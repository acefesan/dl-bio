"""JSON for the 'Gene sets and pathways' page.

kegg_adora.json  (always): the KEGG pathways that contain any adenosine
                  receptor, with size and which receptors; plus how many of
                  all KEGG human pathways each receptor appears in.
pathway_tissue.json (after census_pathway_sample.py has run): per
                  (tissue, nucleus/cell, normal) the weighted mean fraction of
                  each pathway's listed genes detected per cell, plus mean
                  genes detected per cell (the depth proxy).
"""
import json, os
from pathlib import Path
import numpy as np, pandas as pd

R = Path("results/census_overview")
SITE = Path(os.environ.get("SITE_DATA", "../../../../blogs/a1-receptor-pathway/scgsa-site/data")).resolve()
K = json.loads(Path("results/kegg/kegg_human.json").read_text())
P = K["pathways"]; ADORA = ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]
mem = {g: sorted(p for p, d in P.items() if g in d["genes"]) for g in ADORA}
ids = sorted(set().union(*mem.values()))
out = {"fetched": K["fetched"], "n_pathways": len(P), "membership": mem,
       "pathways": {p: {"name": P[p]["name"], "size": len(P[p]["genes"]), "adora": [g for g in ADORA if p in mem[g]]} for p in ids}}
# distribution of pathway sizes for context
sizes = np.array([len(d["genes"]) for d in P.values()])
out["size_quantiles"] = {q: int(np.percentile(sizes, q)) for q in (10, 25, 50, 75, 90)}
SITE.mkdir(parents=True, exist_ok=True)
(SITE / "kegg_adora.json").write_text(json.dumps(out, separators=(",", ":")))
print("kegg_adora.json", len(ids), "pathways")

pg = R / "path_group_homo_sapiens.parquet"
if pg.exists():
    gp = R / "groups_homo_sapiens.parquet"
    g = pd.read_parquet(gp if gp.exists() else R / "groups_meta_homo_sapiens.parquet")
    dep = pd.read_parquet(R / "depth_group_homo_sapiens.parquet").set_index("gid")
    d = pd.read_parquet(pg)
    rows = []
    g = g.set_index("gid")
    d["tissue"] = g.loc[d.gid, "tissue_general"].astype(str).to_numpy()
    d["susp"] = g.loc[d.gid, "suspension_type"].astype(str).str[0].to_numpy()
    d["normal"] = g.loc[d.gid, "is_normal"].astype(int).to_numpy()
    tot = {}
    for (t, s, n), x in d.groupby(["tissue", "susp", "normal"]):
        rec = {"tissue": t, "susp": s, "normal": int(n), "n": 0, "p": {}}
        for p, y in x.groupby("pathway"):
            rec["p"][p] = [round(float(y.w_frac_sum.sum() / y.w_cells.sum()), 5), int(y.n_sampled.sum())]
            rec["n"] = int(y.w_cells.sum())
        sel = g[(g.tissue_general.astype(str) == t) & (g.suspension_type.astype(str).str[0] == s) & (g.is_normal.astype(int) == n)]
        w = dep.loc[sel.index, "n_cells"]
        rec["depth"] = round(float((dep.loc[sel.index, "mean_nnz"] * w).sum() / w.sum()), 1)
        rows.append(rec)
    (SITE / "pathway_tissue.json").write_text(json.dumps({"pathways": {p: out["pathways"][p]["name"] for p in ids}, "rows": rows}, separators=(",", ":")))
    print("pathway_tissue.json", len(rows))
