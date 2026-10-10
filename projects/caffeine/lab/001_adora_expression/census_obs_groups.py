"""Collapse the whole Census obs table into a compact group-count table.

Why: obs has ~100M+ rows per organism; every overview figure only needs counts
per (assay, tissue, cell type, dataset, ...) combination. Streaming obs in
Arrow batches and counting combinations keeps memory small and gives one file
all the downstream overview scripts can share.

Output: results/census_overview/obs_groups_<organism>.parquet with one row per
unique combination and an `n_cells` column. `n_measured_vars` is the number
of genes the cell's *dataset* measured (the Census records this per cell).
"""
import os, sys
from pathlib import Path
import pandas as pd, pyarrow as pa, tiledbsoma as soma

URI = os.environ.get("CENSUS_URI", "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")
OUT = Path("results/census_overview"); OUT.mkdir(parents=True, exist_ok=True)
COLS = ["dataset_id", "assay", "assay_ontology_term_id", "cell_type",
        "cell_type_ontology_term_id", "tissue", "tissue_ontology_term_id",
        "tissue_general", "disease", "suspension_type", "is_primary_data",
        "tissue_type", "n_measured_vars"]

def run(org):
    acc = None
    with soma.open(URI) as c:
        obs = c["census_data"][org].obs
        n = 0
        for batch in obs.read(column_names=COLS):
            df = batch.to_pandas()
            df["n"] = 1
            g = df.groupby(COLS, observed=True, dropna=False, sort=False)["n"].sum()
            acc = g if acc is None else acc.add(g, fill_value=0)
            n += len(df)
            if True:
                print(org, f"{n/1e6:.0f}M rows, {len(acc)} groups", flush=True)
    out = acc.reset_index().rename(columns={"n": "n_cells"})
    out["n_cells"] = out["n_cells"].astype("int64")
    out.to_parquet(OUT / f"obs_groups_{org}.parquet")
    print(org, "total", out.n_cells.sum(), "groups", len(out))

if __name__ == "__main__":
    for org in (sys.argv[1:] or ["homo_sapiens", "mus_musculus"]):
        run(org)
