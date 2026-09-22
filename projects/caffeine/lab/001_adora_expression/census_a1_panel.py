"""A1 minimal-panel expression across the whole local CELLxGENE Census.

Design note: with ~100M unique human cells, per-cell rank-based scGSA scoring
is not tractable and is also unnecessary here — an 8-gene panel is small
enough to report gene by gene, which is what the method-comparison work
concluded. So this computes exact per-gene statistics over ALL cells rather
than scoring a subsample.

Efficiency: X is read for only the 8 panel gene columns, so only their
nonzeros come back (tens of millions of rows), instead of touching the full
61,497-gene axis.

Outputs per organism, into results/census/:
  <org>_by_cell_type.csv     n_cells, pct_expressing, mean_all, mean_expressing
  <org>_by_tissue.csv        same, grouped by tissue_general
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tiledbsoma as soma

URI = os.environ.get("CENSUS_URI",
                    "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")
OUT = Path(os.environ.get("CENSUS_OUT", "results/census"))
OUT.mkdir(parents=True, exist_ok=True)

PANEL = ["ADORA1", "GNAO1", "GNAI1", "KCNJ3", "KCNJ6", "CACNA1A",
         "CACNA1B", "ADCY5"]

OBS_COLS = ["soma_joinid", "cell_type", "tissue_general", "disease",
            "is_primary_data"]


def run(collection: str, org: str) -> None:
    tag = f"{collection}:{org}"
    print(f"\n========== {tag} ==========", flush=True)
    with soma.open(URI) as census:
        exp = census[collection][org]

        var = exp.ms["RNA"].var.read().concat().to_pandas()
        # Match case-insensitively: mouse orthologs are Adora1/Kcnj3, not
        # ADORA1/KCNJ3, so an exact uppercase match silently finds nothing.
        want = {g.upper(): g for g in PANEL}
        upper = var.feature_name.str.upper()
        hits = var[upper.isin(want)][["soma_joinid", "feature_name"]].copy()
        hits["canonical"] = hits.feature_name.str.upper().map(want)
        found = hits.canonical.tolist()
        missing = sorted(set(PANEL) - set(found))
        print(f"panel genes: {len(found)}/{len(PANEL)}"
              + (f"  missing: {missing}" if missing else ""))
        if hits.empty:
            print("no panel genes; skipping")
            return
        gene_of = dict(zip(hits.soma_joinid, hits.canonical))
        gene_ids = hits.soma_joinid.tolist()

        cols = [c for c in OBS_COLS if c in exp.obs.schema.names]
        obs = exp.obs.read(column_names=cols).concat().to_pandas()
        print(f"obs rows: {len(obs):,}", flush=True)
        if "is_primary_data" in obs:
            obs = obs[obs.is_primary_data]
            print(f"primary-only: {len(obs):,}", flush=True)

        obs = obs.set_index("soma_joinid")

        group_cols = [c for c in ("cell_type", "tissue_general")
                      if c in obs.columns]
        # Categorical codes keep the per-cell label lookup as int16/int32
        # instead of ~100M Python strings.
        lookups, cats, counts = {}, {}, {}
        for gc in group_cols:
            s = obs[gc].astype("category")
            cats[gc] = s.cat.categories
            lookups[gc] = pd.Series(s.cat.codes.values, index=obs.index)
            counts[gc] = s.value_counts().rename("n_cells")
        del obs

        # Only the 8 gene columns are touched, so this returns just their
        # nonzeros rather than a dense 100M x 8 block. Aggregate per chunk so
        # peak memory stays bounded regardless of total nonzero count.
        X = exp.ms["RNA"].X["raw"]
        acc = {gc: [] for gc in group_cols}
        total_nz = kept_nz = 0
        for i, tbl in enumerate(X.read(coords=(slice(None), gene_ids)).tables()):
            cell = tbl.column("soma_dim_0").to_numpy()
            gene = tbl.column("soma_dim_1").to_numpy()
            val = tbl.column("soma_data").to_numpy()
            total_nz += len(cell)
            for gc in group_cols:
                code = lookups[gc].reindex(cell).values  # NaN = non-primary
                ok = ~pd.isna(code)
                if not ok.any():
                    continue
                chunk = pd.DataFrame({
                    "gcode": code[ok].astype(np.int32),
                    "gene": gene[ok],
                    "val": val[ok],
                })
                acc[gc].append(chunk.groupby(["gcode", "gene"], observed=True)
                                    .agg(n_expressing=("val", "size"),
                                         total=("val", "sum")))
                if gc == group_cols[0]:
                    kept_nz += int(ok.sum())
            if i % 20 == 0:
                print(f"  chunk {i}: {total_nz:,} nonzeros seen", flush=True)
        print(f"nonzero entries: {total_nz:,} (primary: {kept_nz:,})",
              flush=True)

        for group_col, label in (("cell_type", "by_cell_type"),
                                 ("tissue_general", "by_tissue")):
            if group_col not in group_cols or not acc[group_col]:
                continue
            agg = (pd.concat(acc[group_col])
                     .groupby(level=["gcode", "gene"], observed=True)
                     .sum().reset_index())
            agg["group"] = cats[group_col][agg.gcode.values]
            agg["gene"] = agg.gene.map(gene_of)
            agg = agg.drop(columns=["gcode"])
            n_cells = counts[group_col]
            agg = agg.merge(n_cells, left_on="group", right_index=True)
            agg["pct_expressing"] = 100 * agg.n_expressing / agg.n_cells
            agg["mean_all"] = agg.total / agg.n_cells
            agg["mean_expressing"] = agg.total / agg.n_expressing
            agg.insert(0, "organism", org)
            agg.insert(0, "collection", collection)

            path = OUT / f"{collection}__{org}__{label}.csv"
            agg.sort_values(["gene", "mean_all"], ascending=[True, False]) \
               .to_csv(path, index=False)
            print(f"wrote {path.name}  ({len(agg):,} rows, "
                  f"{agg.group.nunique():,} groups)", flush=True)

            a1 = agg[agg.gene == "ADORA1"].nlargest(8, "mean_all")
            if not a1.empty:
                print(f"  top ADORA1 {label}:")
                print(a1[["group", "n_cells", "pct_expressing", "mean_all"]]
                      .to_string(index=False))


if __name__ == "__main__":
    targets = [("census_data", "homo_sapiens")]
    if "--all" in sys.argv:
        targets += [("census_data", "macaca_mulatta"),
                    ("census_data", "callithrix_jacchus"),
                    ("census_data", "pan_troglodytes"),
                    ("census_spatial_sequencing", "homo_sapiens"),
                    ("census_spatial_sequencing", "mus_musculus")]
    for coll, org in targets:
        try:
            run(coll, org)
        except Exception as e:
            print(f"FAILED {coll}:{org}: {type(e).__name__}: {e}", flush=True)
