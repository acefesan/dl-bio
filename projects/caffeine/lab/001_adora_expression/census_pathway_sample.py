"""KEGG-pathway detection per cell, on a stratified sample of 10x primary cells.

Scanning ~1,300 pathway genes over 89M cells is heavy and unnecessary for a
"what does a gene set say about a cell" picture, so this samples: for each
tissue a random subset of datasets (about DS_CELLS_PER_TISSUE cells), then up
to N_PER_GROUP random cells per (tissue, cell type, nucleus/cell, normal)
group inside them. Each sampled cell is weighted by (group size within the
chosen datasets) / (sampled in group), so tissue averages describe the chosen
datasets. The datasets are limited because the Census sits on a spinning disk
and scattered reads would cost a full pass over the 538 GB raw layer.

For every sampled cell and every pathway P in PATHWAYS:
    detected_frac(cell, P) = (# genes of P with raw count > 0)
                             / (# genes of P that the cell's dataset listed)
Dividing by listed genes (not all pathway genes) keeps "gene not in the
dataset" from being scored as "gene not detected".

Also saves mean genes-detected-per-cell (obs.nnz) per group, the depth proxy
the earlier A1 analysis flagged as a confound.

Outputs (results/census_overview/):
  path_group_<org>.parquet    gid, pathway, n_sampled, w_cells, w_frac_sum
  depth_group_<org>.parquet   gid, n_cells, mean_nnz
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, tiledbsoma as soma

URI = os.environ.get("CENSUS_URI", "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")
OUT = Path("results/census_overview")
N_PER_GROUP = int(os.environ.get("N_PER_GROUP", 600))
PATHWAYS = ["hsa04015", "hsa04020", "hsa04022", "hsa04024", "hsa04071", "hsa04080", "hsa04082",
            "hsa04270", "hsa04923", "hsa04924", "hsa05012", "hsa05032", "hsa05034"]
kegg = json.loads(Path("results/kegg/kegg_human.json").read_text())["pathways"]

org = "homo_sapiens"
rng = np.random.default_rng(0)
with soma.open(URI) as c:
    exp = c["census_data"][org]
    ds = c["census_info"]["datasets"].read(column_names=["soma_joinid", "dataset_id"]).concat().to_pandas()
    var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name"]).concat().to_pandas()
    genes = sorted(set().union(*[set(kegg[p]["genes"]) for p in PATHWAYS]))
    v = var[var.feature_name.isin(genes)].drop_duplicates("feature_name")
    gj = v.soma_joinid.to_numpy(); gnames = v.feature_name.tolist(); G = len(gj)
    print("pathway genes found in Census var:", G, "of", len(genes), flush=True)
    gidx = {n: i for i, n in enumerate(gnames)}
    gi = np.full(int(var.soma_joinid.max()) + 1, -1, np.int32); gi[gj] = np.arange(G)
    # membership gene -> pathways (CSR)
    rows, cols = [], []
    for pi, p in enumerate(PATHWAYS):
        for g in kegg[p]["genes"]:
            if g in gidx: rows.append(gidx[g]); cols.append(pi)
    Mem = sp.csr_matrix((np.ones(len(rows), np.int32), (rows, cols)), shape=(G, len(PATHWAYS)))
    NP = len(PATHWAYS)

    obs = exp.obs.read(column_names=["soma_joinid", "dataset_id", "assay", "tissue_general",
                                      "cell_type_ontology_term_id", "cell_type", "disease",
                                      "suspension_type", "nnz"],
                       value_filter="is_primary_data == True").concat().to_pandas()
    obs = obs[obs.assay.astype(str).str.startswith("10x")].copy()
    obs["is_normal"] = obs.disease.astype(str).eq("normal")
    key = ["tissue_general", "cell_type_ontology_term_id", "cell_type", "suspension_type", "is_normal"]
    obs["gid"] = obs.groupby(key, observed=True, sort=False).ngroup().astype(np.int32)
    NGR = int(obs.gid.max()) + 1
    chk = pd.read_parquet(OUT / f"groups_{org}.parquet") if (OUT / f"groups_{org}.parquet").exists() else None
    if chk is not None:
        meta = obs.drop_duplicates("gid").sort_values("gid")
        assert (meta.cell_type_ontology_term_id.astype(str).to_numpy() == chk.cell_type_ontology_term_id.astype(str).to_numpy()).all(), "group ids differ from gene_stats groups"
    dep = obs.groupby("gid").agg(n_cells=("nnz", "size"), mean_nnz=("nnz", "mean")).reset_index()
    dep.to_parquet(OUT / f"depth_group_{org}.parquet")

    # Dataset-limited stratified sample. The Census matrix sits on a spinning
    # disk and cells of one dataset are stored contiguously, so reading cells
    # scattered over all 524 datasets means reading the whole 538 GB layer.
    # Instead, for each tissue pick datasets (random order) until DS_CELLS_PER_TISSUE
    # cells are covered, then sample up to N_PER_GROUP cells per group inside
    # them. Estimates therefore describe those datasets, not the whole population.
    DS_CELLS = int(os.environ.get("DS_CELLS_PER_TISSUE", 150_000))
    dsz = obs.groupby(["tissue_general", "dataset_id"], observed=True).size().reset_index(name="n")
    dsz = dsz[dsz.n > 0].sample(frac=1.0, random_state=0)
    dsz["cum"] = dsz.groupby("tissue_general", observed=True).n.cumsum() - dsz.n
    chosen = dsz[dsz.cum < DS_CELLS][["tissue_general", "dataset_id"]]
    obs = obs.merge(chosen.assign(_c=1), on=["tissue_general", "dataset_id"], how="inner")
    print("datasets chosen:", chosen.dataset_id.nunique(), "cells in them:", len(obs), flush=True)
    obs = obs.sample(frac=1.0, random_state=0)
    obs["rank"] = obs.groupby("gid").cumcount()
    samp = obs[obs["rank"] < N_PER_GROUP].sort_values("soma_joinid").reset_index(drop=True)
    gsz = obs.groupby("gid").size(); gsm = samp.groupby("gid").size()
    samp["w"] = (gsz.reindex(samp.gid).to_numpy() / gsm.reindex(samp.gid).to_numpy())
    S = len(samp); print("sampled cells", S, "groups", NGR, flush=True)
    loc = np.full(int(obs.soma_joinid.max()) + 1, -1, np.int64); loc[samp.soma_joinid.to_numpy()] = np.arange(S)

    # listed-gene denominators per sampled cell and pathway
    dsj = ds.set_index("dataset_id").soma_joinid
    samp["dsj"] = samp.dataset_id.astype(str).map(dsj).astype(int)
    P = exp.ms["RNA"]["feature_dataset_presence_matrix"].read((slice(None), gj.tolist())).tables().concat().to_pandas()
    gmap = np.full(int(var.soma_joinid.max()) + 1, -1, np.int64); gmap[gj] = np.arange(G)
    Pm = sp.csr_matrix((np.ones(len(P), np.int32), (P.soma_dim_0, gmap[P.soma_dim_1.to_numpy()])),
                       shape=(int(ds.soma_joinid.max()) + 1, G))
    listed = (Pm @ Mem).toarray()[samp.dsj.to_numpy()]            # S x NP

    det = np.zeros(S * NP, np.int64)
    ptr, ind = Mem.indptr, Mem.indices
    n = 0
    for t in exp.ms["RNA"].X["raw"].read((samp.soma_joinid.to_numpy(), gj.tolist())).tables():
        ci = loc[t.column("soma_dim_0").to_numpy()]; g = gi[t.column("soma_dim_1").to_numpy()]
        k = ptr[g + 1] - ptr[g]
        rep = np.repeat(np.arange(len(g)), k)
        off = np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)
        p = ind[ptr[g][rep] + off]
        det += np.bincount(ci[rep].astype(np.int64) * NP + p, minlength=S * NP)
        n += len(g)
        if n % 100_000_000 < len(g): print(f"  {n/1e6:.0f}M entries", flush=True)
    det = det.reshape(S, NP)

frac = np.divide(det, listed, out=np.zeros_like(det, dtype=float), where=listed > 0)
res = []
for pi, p in enumerate(PATHWAYS):
    d = pd.DataFrame({"gid": samp.gid.to_numpy(), "w": samp.w.to_numpy(), "f": frac[:, pi] * samp.w.to_numpy(),
                      "ok": (listed[:, pi] > 0)})
    g = d.groupby("gid").agg(n_sampled=("w", "size"), w_cells=("w", "sum"), w_frac_sum=("f", "sum")).reset_index()
    g["pathway"] = p; res.append(g)
pd.concat(res).to_parquet(OUT / f"path_group_{org}.parquet")
print("done")
