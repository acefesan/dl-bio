"""Per-(tissue, cell type, nucleus/cell) statistics for a panel of adenosine-
signalling genes, over every primary 10x cell in the local Census.

How it works (so the figures built on it are interpretable):

1. Read the obs columns for all *primary* cells and keep the 10x family
   (assay starts with "10x"). "Primary" avoids counting the same cell twice
   when several datasets re-publish it.
2. A "group" is a unique (tissue_general, cell_type, suspension_type,
   is_normal) combination. Every cell gets an integer group code.
3. X is read only for the panel's gene columns, in Arrow batches. For each
   nonzero entry we add 1 to n_expr[gene, group] and add the value to
   sum_frac[gene, group]. The Census "normalized" layer is count / cell
   total, so cp10k = 1e4 * sum_frac / n_cells.
4. Zeros are not all alike. The Census presence matrix says which datasets
   *included each gene in their feature list at all*. A cell from a dataset
   that did not list a gene cannot have counted it, so its zero means "not
   measured". n_measured[gene, group] counts only cells whose dataset listed
   the gene; use it as the denominator for "% of cells expressing".

Output (results/census_overview/):
  groups_<org>.parquet        group_id + metadata + n_cells
  gene_stats_<org>.parquet    group_id, gene, n_measured, n_expr, sum_frac
"""
import os, sys, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, tiledbsoma as soma

URI = os.environ.get("CENSUS_URI", "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")
OUT = Path("results/census_overview"); OUT.mkdir(parents=True, exist_ok=True)
GENES = json.loads(Path("adenosine_genes.json").read_text())  # {symbol: category}

def run(org):
    print("=====", org, flush=True)
    with soma.open(URI) as c:
        exp = c["census_data"][org]
        ds = c["census_info"]["datasets"].read(column_names=["soma_joinid", "dataset_id"]).concat().to_pandas()
        var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name", "feature_id"]).concat().to_pandas()
        up = var.feature_name.str.upper()
        sel = var[up.isin([g.upper() for g in GENES])].copy()
        sel["symbol"] = sel.feature_name.str.upper().map({g.upper(): g for g in GENES})
        sel = sel.drop_duplicates("symbol").reset_index(drop=True)
        missing = sorted(set(GENES) - set(sel.symbol))
        print("genes found", len(sel), "missing", missing, flush=True)
        gj = sel.soma_joinid.to_numpy(); G = len(sel)
        gi = np.full(var.soma_joinid.max() + 1, -1, np.int32); gi[gj] = np.arange(G)

        obs = exp.obs.read(column_names=["soma_joinid", "dataset_id", "assay", "tissue_general",
                                          "cell_type", "cell_type_ontology_term_id", "disease",
                                          "suspension_type"],
                           value_filter="is_primary_data == True").concat().to_pandas()
        obs = obs[obs.assay.astype(str).str.startswith("10x")].copy()
        obs["is_normal"] = obs.disease.astype(str).eq("normal")
        print("10x primary cells", len(obs), flush=True)
        key = ["tissue_general", "cell_type_ontology_term_id", "cell_type", "suspension_type", "is_normal"]
        obs["gid"] = obs.groupby(key, observed=True, sort=False).ngroup().astype(np.int32)
        groups = obs.drop_duplicates("gid")[["gid"] + key].sort_values("gid").reset_index(drop=True)
        NG = len(groups)
        groups["n_cells"] = np.bincount(obs.gid, minlength=NG)
        print("groups", NG, flush=True)

        code = np.full(int(obs.soma_joinid.max()) + 1, -1, np.int32)
        code[obs.soma_joinid.to_numpy()] = obs.gid.to_numpy()

        # measured denominators via dataset presence
        dsj = ds.set_index("dataset_id").soma_joinid
        obs["dsj"] = obs.dataset_id.astype(str).map(dsj)
        P = exp.ms["RNA"]["feature_dataset_presence_matrix"].read((slice(None), gj.tolist())).tables().concat().to_pandas()
        Pm = sp.csr_matrix((np.ones(len(P), np.int8), (P.soma_dim_0, P.soma_dim_1.map(dict(zip(gj, range(G)))))),
                           shape=(int(ds.soma_joinid.max()) + 1, G))
        M = sp.csr_matrix((np.ones(len(obs), np.int32), (obs.gid, obs.dsj.astype(int))), shape=(NG, Pm.shape[0]))
        n_meas = (M @ Pm).toarray()                      # NG x G

        n_expr = np.zeros(G * NG, np.int64); sum_frac = np.zeros(G * NG, np.float64)
        rows = 0
        for t in exp.ms["RNA"].X["normalized"].read((slice(None), gj.tolist())).tables():
            d0 = t.column("soma_dim_0").to_numpy(); d1 = t.column("soma_dim_1").to_numpy()
            v = t.column("soma_data").to_numpy()
            gc = np.full(len(d0), -1, np.int32); inr = d0 < len(code); gc[inr] = code[d0[inr]]; ok = gc >= 0
            idx = gi[d1[ok]].astype(np.int64) * NG + gc[ok]
            n_expr += np.bincount(idx, minlength=G * NG)
            sum_frac += np.bincount(idx, weights=v[ok], minlength=G * NG)
            rows += len(d0)
            if rows % 50_000_000 < len(d0): print(f"  X rows {rows/1e6:.0f}M", flush=True)
        n_expr = n_expr.reshape(G, NG); sum_frac = sum_frac.reshape(G, NG)

    long = pd.DataFrame({
        "gene": np.repeat(sel.symbol.to_numpy(), NG),
        "group_id": np.tile(np.arange(NG), G),
        "n_measured": n_meas.T.reshape(-1),
        "n_expr": n_expr.reshape(-1), "sum_frac": sum_frac.reshape(-1)})
    long = long[(long.n_measured > 0) | (long.n_expr > 0)]
    groups.to_parquet(OUT / f"groups_{org}.parquet"); long.to_parquet(OUT / f"gene_stats_{org}.parquet")
    print(org, "done", len(groups), len(long), "missing", missing, flush=True)

if __name__ == "__main__":
    for org in (sys.argv[1:] or ["homo_sapiens", "mus_musculus"]):
        run(org)
