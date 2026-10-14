"""Sample Census cells per cell type and export full raw counts for scGSA scoring.

For the pathways page's section 4 gene sets (A1 in neurons, A2A in striatum,
vasodilation), score whole transcriptomes with the seven methods, which needs
every gene of each sampled cell (the rank-based methods rank within a cell).
So: for each chosen cell type, take up to N_PER_TYPE cells (10x, primary,
normal) from at most MAX_DS datasets that measured at least 90% of the gene
sets' genes, and read their full raw rows. Datasets are limited because the
Census sits on a spinning disk and scattered reads cost a full pass.

Output (SCGSA_WORK/data, default /mnt/bulk/dl_bio/scgsa_census/data):
  census_counts.h5   10x-style HDF5, genes x cells, raw counts (Read10X_h5)
  census_cells.csv   barcode, group (cell type), dataset_id, tissue_general,
                     assay, suspension_type, nnz (genes detected)
Run: uv run --no-project --python 3.12 --with tiledbsoma --with pandas
     --with pyarrow --with scipy --with h5py python census_scgsa_sample.py
"""
import json, os, time
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp, tiledbsoma as soma

URI = os.environ.get("CENSUS_URI", "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")
OUT = Path(os.environ.get("SCGSA_WORK", "/mnt/bulk/dl_bio/scgsa_census")) / "data"
OUT.mkdir(parents=True, exist_ok=True)
HERE = Path(__file__).resolve().parent
N_PER_TYPE, MAX_DS, MIN_IN_DS, SEED = 200, 5, 40, 0

CELL_TYPES = [
    # vascular: the vasodilation question
    "endothelial cell of artery", "endothelial cell of arteriole", "capillary endothelial cell",
    "vein endothelial cell", "endothelial cell of venule", "endothelial cell of lymphatic vessel",
    "cerebral cortex endothelial cell", "brain vascular cell", "pericyte",
    "vascular associated smooth muscle cell", "blood vessel smooth muscle cell",
    "microcirculation associated smooth muscle cell", "smooth muscle cell of the brain vasculature",
    # neurons: the A1 / A2A diagrams
    "corticothalamic-projecting glutamatergic cortical neuron", "near-projecting glutamatergic cortical neuron",
    "L6b glutamatergic cortical neuron", "L5 extratelencephalic projecting glutamatergic cortical neuron",
    "L2/3 intratelencephalic projecting glutamatergic neuron", "pvalb GABAergic cortical interneuron",
    "sst GABAergic cortical interneuron", "direct pathway medium spiny neuron",
    "indirect pathway medium spiny neuron", "midbrain dopaminergic neuron",
    # glia
    "astrocyte", "microglial cell", "oligodendrocyte", "oligodendrocyte precursor cell",
    # contrast
    "enteric smooth muscle cell", "fibroblast", "regular ventricular cardiac myocyte",
    "regular atrial cardiac myocyte", "hepatocyte", "kidney loop of Henle thick ascending limb epithelial cell",
    "naive B cell", "CD14-positive monocyte", "naive thymus-derived CD4-positive, alpha-beta T cell",
    "CD16-positive, CD56-dim natural killer cell, human", "mast cell",
]

sets = json.loads((HERE / "results/census_scgsa/gene_sets.json").read_text())
set_genes = sorted(set().union(*map(set, sets.values())))
rng = np.random.default_rng(SEED)
t0 = time.time()
with soma.open(URI) as c:
    exp = c["census_data"]["homo_sapiens"]
    var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name"]).concat().to_pandas()
    ds_info = c["census_info"]["datasets"].read(column_names=["soma_joinid", "dataset_id"]).concat().to_pandas()
    # which datasets measured which genes
    sg = var[var.feature_name.isin(set_genes)]
    pres = exp.ms["RNA"]["feature_dataset_presence_matrix"].read((slice(None), sg.soma_joinid.tolist())).tables().concat().to_pandas()
    frac = pres.groupby("soma_dim_0").size() / len(sg)
    ok_ds = set(ds_info.set_index("soma_joinid").loc[frac[frac >= 0.9].index, "dataset_id"])
    print(f"{len(sg)} set genes in Census; {len(ok_ds)} datasets measured >=90% of them", flush=True)

    flt = " or ".join(f"cell_type == '{t}'" for t in CELL_TYPES)
    obs = exp.obs.read(column_names=["soma_joinid", "cell_type", "dataset_id", "tissue_general", "assay",
                                      "suspension_type", "disease", "is_primary_data", "nnz"],
                       value_filter=f"is_primary_data == True and disease == 'normal' and ({flt})").concat().to_pandas()
    obs = obs[obs.assay.str.contains("10x") & obs.dataset_id.isin(ok_ds)]
    print(f"{len(obs):,} candidate cells in {obs.cell_type.nunique()} types ({time.time()-t0:.0f}s)", flush=True)

    picks = []
    for t in CELL_TYPES:
        o = obs[obs.cell_type == t]
        per = o.groupby("dataset_id").size()
        per = per[per >= MIN_IN_DS]
        if per.empty:
            print(f"  skip {t}: no dataset with >= {MIN_IN_DS} cells", flush=True); continue
        dsel = rng.choice(per.index.to_numpy(), size=min(MAX_DS, len(per)), replace=False)
        o = o[o.dataset_id.isin(dsel)]
        k = min(N_PER_TYPE, len(o))
        picks.append(o.sample(n=k, random_state=int(rng.integers(1 << 31))))
        print(f"  {t}: {k} cells from {len(dsel)} datasets", flush=True)
    samp = pd.concat(picks).sort_values("soma_joinid").reset_index(drop=True)
    S = len(samp)
    print(f"reading full rows for {S:,} cells ...", flush=True)

    loc = pd.Series(np.arange(S), index=samp.soma_joinid.to_numpy())
    G = int(var.soma_joinid.max()) + 1
    rows, cols, vals = [], [], []
    for t in exp.ms["RNA"].X["raw"].read((samp.soma_joinid.to_numpy(), slice(None))).tables():
        rows.append(loc.loc[t.column("soma_dim_0").to_numpy()].to_numpy())
        cols.append(t.column("soma_dim_1").to_numpy()); vals.append(t.column("soma_data").to_numpy())
    mat = sp.csr_matrix((np.concatenate(vals).astype(np.float32), (np.concatenate(rows), np.concatenate(cols))), shape=(S, G))
    print(f"read nnz={mat.nnz:,} in {time.time()-t0:.0f}s", flush=True)

symbols = var.set_index("soma_joinid").reindex(range(G)).feature_name.fillna("NA").to_numpy()
# make names unique for Seurat
seen = {}
for i, s in enumerate(symbols):
    if s in seen:
        seen[s] += 1; symbols[i] = f"{s}.{seen[s]}"
    else:
        seen[s] = 0
barcodes = np.array([f"cj{j}" for j in samp.soma_joinid])
with h5py.File(OUT / "census_counts.h5", "w") as o:
    m = o.create_group("matrix")
    m.create_dataset("data", data=mat.data, compression="gzip", compression_opts=1)
    m.create_dataset("indices", data=mat.indices.astype(np.int64), compression="gzip", compression_opts=1)
    m.create_dataset("indptr", data=mat.indptr.astype(np.int64))
    m.create_dataset("shape", data=np.array([G, S], dtype=np.int32))
    m.create_dataset("barcodes", data=barcodes.astype("S"))
    f = m.create_group("features")
    f.create_dataset("id", data=symbols.astype("S")); f.create_dataset("name", data=symbols.astype("S"))
    f.create_dataset("feature_type", data=np.array(["Gene Expression"] * G, dtype="S"))
    f.create_dataset("genome", data=np.array(["GRCh38"] * G, dtype="S"))
samp.assign(barcode=barcodes, group=samp.cell_type)[["barcode", "group", "dataset_id", "tissue_general", "assay",
                                                       "suspension_type", "nnz"]].to_csv(OUT / "census_cells.csv", index=False)
print(f"wrote {OUT}/census_counts.h5 ({G:,} genes x {S:,} cells) in {time.time()-t0:.0f}s", flush=True)
