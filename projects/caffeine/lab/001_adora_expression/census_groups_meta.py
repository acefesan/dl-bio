"""Re-derive the group-id -> (tissue, cell type, nucleus/cell, normal) table.

census_gene_groups.py and census_pathway_sample.py both number groups with
pandas ngroup() over the same filtered obs read. This writes that table on its
own so the pathway page can be built before the long gene scan finishes, and
checks it against the group sizes the pathway sample recorded.
"""
import os, numpy as np, pandas as pd, tiledbsoma as soma
URI = "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma"
with soma.open(URI) as c:
    obs = c["census_data"]["homo_sapiens"].obs.read(
        column_names=["soma_joinid", "assay", "tissue_general", "cell_type", "cell_type_ontology_term_id",
                      "disease", "suspension_type", "nnz"], value_filter="is_primary_data == True").concat().to_pandas()
obs = obs[obs.assay.astype(str).str.startswith("10x")].copy()
obs["is_normal"] = obs.disease.astype(str).eq("normal")
key = ["tissue_general", "cell_type_ontology_term_id", "cell_type", "suspension_type", "is_normal"]
obs["gid"] = obs.groupby(key, observed=True, sort=False).ngroup().astype(np.int32)
g = obs.drop_duplicates("gid")[["gid"] + key].sort_values("gid").reset_index(drop=True)
g["n_cells"] = np.bincount(obs.gid, minlength=len(g))
dep = pd.read_parquet("results/census_overview/depth_group_homo_sapiens.parquet").set_index("gid")
assert (dep.loc[g.gid, "n_cells"].to_numpy() == g.n_cells.to_numpy()).all(), "group ids do not match the pathway sample"
print("ok: matches depth_group on", len(g), "groups")
g.to_parquet("results/census_overview/groups_meta_homo_sapiens.parquet")
