"""Turn the Census group tables into the JSON the site's Census pages load.

Inputs  (results/census_overview/):
  obs_groups_<org>.parquet   cell counts, used for the tree skeleton
  groups_<org>.parquet       (tissue, cell type, nucleus/cell, normal) groups
  gene_stats_<org>.parquet   per-group, per-gene n_measured / n_expr / sum_frac
Outputs (site data dir, default ../../../../blogs/a1-receptor-pathway/scgsa-site/data):
  tree_<org>.json            cell-type tree with per-node ADORA stats
  adora_<org>.json           (tissue x cell type) rows for the four receptors
  panel_<org>.json           (tissue) rows for the 60-gene panel
  panel_celltype_<org>.json  (cell type) rows for the 60-gene panel

Cell-type tree: the Cell Ontology is a DAG (a type can have several parents).
Each node's "total" counts every cell whose annotated type is that node or any
descendant, so it is exact. For drawing, each node is placed under ONE parent:
the parent with the most cells. Other parents are kept in `also` and shown in
the tooltip. A node's drawn children therefore need not sum to its total.
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp

R = Path("results/census_overview")
SITE = Path(os.environ.get("SITE_DATA", "../../../../blogs/a1-receptor-pathway/scgsa-site/data")).resolve()
ADORA = {"homo_sapiens": ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"],
         "mus_musculus": ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]}  # stats use canonical upper-case keys for both organisms
COMBOS = [(n, s) for n in ("N", "A") for s in ("B", "C", "U")]   # normal-only / all x both / cell / nucleus

def combo_mask(df, n, s):
    m = np.ones(len(df), bool)
    if n == "N": m &= df.is_normal.to_numpy()
    if s == "C": m &= (df.suspension_type.astype(str) == "cell").to_numpy()
    if s == "U": m &= (df.suspension_type.astype(str) == "nucleus").to_numpy()
    return m

def load_groups(org):
    p = R / f"groups_{org}.parquet"
    if p.exists():
        g = pd.read_parquet(p)
    else:
        d = pd.read_parquet(R / f"obs_groups_{org}.parquet")
        d = d[d.is_primary_data & d.assay.astype(str).str.startswith("10x")].copy()
        d["is_normal"] = d.disease.astype(str).eq("normal")
        g = d.groupby(["tissue_general", "cell_type_ontology_term_id", "cell_type", "suspension_type", "is_normal"],
                      observed=True).n_cells.sum().reset_index()
        g.insert(0, "gid", np.arange(len(g)))
    for c in ["tissue_general", "cell_type_ontology_term_id", "cell_type", "suspension_type"]:
        g[c] = g[c].astype(str)
    return g

def build(org):
    from cellxgene_ontology_guide.ontology_parser import OntologyParser
    op = OntologyParser()
    g = load_groups(org)
    gs_path = R / f"gene_stats_{org}.parquet"
    stats = pd.read_parquet(gs_path) if gs_path.exists() else None
    genes = ADORA[org]
    types = sorted(t for t in g.cell_type_ontology_term_id.unique() if t.startswith("CL:"))
    label = {t: g[g.cell_type_ontology_term_id == t].cell_type.iloc[0] for t in types}
    anc = {}
    for t in types:
        try: a = set(op.get_term_ancestors(t, include_self=True))
        except Exception: a = {t}
        anc[t] = {x for x in a if x.startswith("CL:")}
    nodes = sorted(set().union(*anc.values()))
    nid = {n: i for i, n in enumerate(nodes)}; tid = {t: i for i, t in enumerate(types)}
    A = sp.csr_matrix(([1] * sum(len(a) for a in anc.values()),
                       ([nid[n] for t in types for n in anc[t]], [tid[t] for t in types for n in anc[t]])),
                      shape=(len(nodes), len(types)))
    gt = g[g.cell_type_ontology_term_id.isin(types)]
    own = gt.groupby("cell_type_ontology_term_id").n_cells.sum().reindex(types).fillna(0).to_numpy()
    total = A @ own
    # parents / drawn parent
    par, also = {}, {}
    for n in nodes:
        if n == "CL:0000000": continue
        try: ps = [p for p in op.get_term_parents(n) if p in nid]
        except Exception: ps = []
        if not ps: ps = ["CL:0000000"]
        best = max(ps, key=lambda p: total[nid[p]])
        par[n] = best; also[n] = [op.get_term_label(p) for p in ps if p != best]
    out = {"nodes": {}, "root": "CL:0000000", "genes": genes, "combos": [n + s for n, s in COMBOS]}
    for n in nodes:
        try: nm = op.get_term_label(n)
        except Exception: nm = n
        out["nodes"][n] = {"n": nm, "p": par.get(n), "o": int(own[tid[n]]) if n in tid else 0, "t": int(total[nid[n]]), "a": also.get(n, [])}
    if stats is not None:
        gidx = g.set_index("gid")
        for gene in genes:
            s = stats[stats.gene == gene].set_index("group_id")
            j = g.join(s[["n_measured", "n_expr", "sum_frac"]], on="gid").fillna({"n_measured": 0, "n_expr": 0, "sum_frac": 0})
            j = j[j.cell_type_ontology_term_id.isin(types)]
            for (nn, ss) in COMBOS:
                m = combo_mask(j, nn, ss)
                jj = j[m].groupby("cell_type_ontology_term_id")[["n_cells", "n_measured", "n_expr", "sum_frac"]].sum().reindex(types).fillna(0)
                M = A @ jj.to_numpy()
                for n in nodes:
                    r = M[nid[n]]
                    if r[0] > 0:
                        out["nodes"][n].setdefault("s", {}).setdefault(nn + ss, {})[gene] = [int(r[0]), int(r[1]), int(r[2]), float(f"{r[3]:.5g}")]
    SITE.mkdir(parents=True, exist_ok=True)
    (SITE / f"tree_{org}.json").write_text(json.dumps(out, separators=(",", ":")))
    print(org, "tree nodes", len(nodes), "types", len(types), "stats", stats is not None)

    if stats is None: return
    # --- (tissue x cell type) rows for the receptors
    sub = stats[stats.gene.isin(genes)].pivot_table(index="group_id", columns="gene", values=["n_measured", "n_expr", "sum_frac"], fill_value=0)
    rows = []
    for _, r in g.iterrows():
        if r.n_cells < 100: continue
        gid = r.gid
        vals = []
        for gene in genes:
            if gid in sub.index:
                vals += [int(sub.loc[gid, ("n_measured", gene)]), int(sub.loc[gid, ("n_expr", gene)]), float(f"{sub.loc[gid, ('sum_frac', gene)]:.5g}")]
            else: vals += [0, 0, 0.0]
        rows.append([r.tissue_general, r.cell_type_ontology_term_id, r.cell_type, r.suspension_type[0], int(r.is_normal), int(r.n_cells)] + vals)
    (SITE / f"adora_{org}.json").write_text(json.dumps({"genes": genes, "cols": ["tissue", "ct_id", "ct", "susp", "normal", "n"] + [f"{g_}:{k}" for g_ in genes for k in ("meas", "expr", "sumfrac")], "rows": rows}, separators=(",", ":")))
    print(org, "adora rows", len(rows))

    # --- 60-gene panel by tissue and by cell type
    panel = sorted(stats.gene.unique())
    st = stats.merge(g[["gid", "tissue_general", "cell_type_ontology_term_id", "cell_type", "suspension_type", "is_normal"]], left_on="group_id", right_on="gid")
    def agg(key, name):
        keys = key + ["suspension_type", "is_normal"]
        ncell = g.groupby(keys).n_cells.sum()
        t = st.groupby(keys + ["gene"])[["n_measured", "n_expr", "sum_frac"]].sum()
        recs = []
        for k, nc in ncell.items():
            if nc < 200: continue
            rec = {"k": list(k[:len(key)]), "susp": k[len(key)][0], "normal": int(k[len(key) + 1]), "n": int(nc), "g": {}}
            recs.append((k, rec))
        # fill genes
        tt = t.reset_index()
        d = {}
        nk = len(keys)
        for row in tt[keys + ["gene", "n_measured", "n_expr", "sum_frac"]].itertuples(index=False):
            k = tuple(row[:nk]); gene, a, b, c = row[nk:]
            d.setdefault(k, {})[gene] = [int(a), int(b), float(f"{c:.4g}")]
        for k, rec in recs: rec["g"] = d.get(tuple(k), {})
        (SITE / f"{name}_{org}.json").write_text(json.dumps({"genes": panel, "rows": [r for _, r in recs]}, separators=(",", ":")))
        print(org, name, len(recs))
    agg(["tissue_general"], "panel")
    agg(["cell_type_ontology_term_id", "cell_type"], "panel_celltype")

if __name__ == "__main__":
    for org in (sys.argv[1:] or ["homo_sapiens", "mus_musculus"]):
        build(org)
