"""Overview numbers for the 'What is in the Census' page.

Everything is derived from two sources:
  * results/census_overview/obs_groups_<org>.parquet  (cell counts per
    dataset/assay/tissue/cell-type/... combination; made by census_obs_groups.py)
  * the Census presence matrix (which datasets listed which genes) and
    var table, read directly.

"10x" below means any assay whose name starts with "10x" (3' v1/v2/v3, 5'
v1/v2, multiome, Flex). "primary" means is_primary_data == True, i.e. each
real cell counted once.
"""
import json, os
from pathlib import Path
import numpy as np, pandas as pd, tiledbsoma as soma

URI = os.environ.get("CENSUS_URI", "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")
R = Path("results/census_overview")
out = {"census_build": "2025-11-08"}

def family(a):
    a = str(a)
    if a.startswith("10x"): return "10x"
    if a.startswith("Smart-seq") or a in ("Smart-seq",): return "Smart-seq family"
    if a.startswith("sci-RNA"): return "sci-RNA-seq3"
    return "other"

with soma.open(URI) as c:
    dsinfo = c["census_info"]["datasets"].read(column_names=["soma_joinid", "dataset_id", "dataset_title", "collection_name"]).concat().to_pandas()
    for org in ["homo_sapiens", "mus_musculus"]:
        d = pd.read_parquet(R / f"obs_groups_{org}.parquet")
        d["assay"] = d.assay.astype(str); d["family"] = d.assay.map(family)
        d["tissue_general"] = d.tissue_general.astype(str)
        o = {}
        o["total_cells"] = int(d.n_cells.sum())
        o["primary_cells"] = int(d[d.is_primary_data].n_cells.sum())
        p = d[d.is_primary_data]
        o["primary_by_family"] = p.groupby("family", observed=True).n_cells.sum().astype(int).to_dict()
        o["assays"] = [{"assay": a, "family": family(a), "all": int(g.n_cells.sum()),
                        "primary": int(g[g.is_primary_data].n_cells.sum())}
                       for a, g in d.groupby("assay", observed=True) if g.n_cells.sum() > 0]
        o["assays"].sort(key=lambda r: -r["all"])
        x = p[p.family == "10x"]
        o["tenx_primary"] = int(x.n_cells.sum())
        o["tenx_normal"] = int(x[x.disease.astype(str) == "normal"].n_cells.sum())
        o["tenx_suspension"] = x.groupby(x.suspension_type.astype(str), observed=True).n_cells.sum().astype(int).to_dict()
        o["tenx_cell_types"] = int(x.cell_type_ontology_term_id.nunique())
        o["tenx_datasets"] = int(x.dataset_id.nunique())
        t = x.groupby(["tissue_general", x.suspension_type.astype(str)], observed=True).n_cells.sum().unstack(fill_value=0)
        t["total"] = t.sum(axis=1); t = t.sort_values("total", ascending=False)
        o["tenx_tissue"] = [{"tissue": i, **{k: int(v) for k, v in r.items()}} for i, r in t.iterrows()]
        # dataset-level: genes listed
        x = x.assign(dataset_id=x.dataset_id.astype(str))
        dd = x.groupby("dataset_id", observed=True).agg(n_cells=("n_cells", "sum"), n_measured=("n_measured_vars", "max"),
                                         assay=("assay", "first")).reset_index()
        dd = dd[dd.n_cells > 0]
        dd = dd.merge(dsinfo[["dataset_id", "dataset_title", "collection_name"]], on="dataset_id", how="left")
        o["dataset_gene_counts"] = dd[["n_cells", "n_measured", "assay"]].to_dict("records")
        o["smallest_gene_lists"] = (dd.sort_values("n_measured").head(12)
                                    [["dataset_title", "collection_name", "assay", "n_cells", "n_measured"]].to_dict("records"))
        # gene-level breadth over 10x primary cells
        exp = c["census_data"][org]
        var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name", "feature_type"]).concat().to_pandas()
        w = np.zeros(int(dsinfo.soma_joinid.max()) + 1)
        w[dsinfo.set_index("dataset_id").soma_joinid.reindex(dd.dataset_id).to_numpy()] = dd.n_cells.to_numpy()
        meas = np.zeros(int(var.soma_joinid.max()) + 1)
        for t_ in exp.ms["RNA"]["feature_dataset_presence_matrix"].read().tables():
            d0 = t_.column("soma_dim_0").to_numpy(); d1 = t_.column("soma_dim_1").to_numpy()
            np.add.at(meas, d1, w[d0])
        var["frac_cells_measured"] = meas[var.soma_joinid.to_numpy()] / w.sum()
        bins = [0, 1e-9, .25, .5, .75, .9, .99, 1.0 + 1e-9]
        labels = ["never", "<25%", "25-50%", "50-75%", "75-90%", "90-99%", "99-100%"]
        var["bin"] = pd.cut(var.frac_cells_measured, bins, labels=labels, right=False, include_lowest=True)
        o["gene_breadth"] = {"labels": labels,
                             "all": [int(v) for v in var.bin.value_counts().reindex(labels).fillna(0)],
                             "protein_coding": [int(v) for v in var[var.feature_type == "protein_coding"].bin.value_counts().reindex(labels).fillna(0)]}
        o["n_genes"] = int(len(var)); o["n_protein_coding"] = int((var.feature_type == "protein_coding").sum())
        o["genes_frac_measured"] = var.drop_duplicates("feature_name").set_index("feature_name").frac_cells_measured.reindex(
            ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3", "GAPDH"] if org == "homo_sapiens" else
            ["Adora1", "Adora2a", "Adora2b", "Adora3", "Gapdh"]).round(4).to_dict()
        out[org] = o
        print(org, o["total_cells"], o["primary_cells"], o["tenx_primary"], o["genes_frac_measured"], flush=True)
Path("site_data").mkdir(exist_ok=True)
(R / "overview.json").write_text(json.dumps(out, indent=1, default=lambda v: None if pd.isna(v) else str(v)))
print("wrote overview.json")
