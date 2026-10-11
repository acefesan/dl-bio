"""ADORA1/2A/2B/3 in the spatial collection of the local CELLxGENE Census.

The Census (2025-11-08) has a separate experiment group,
``census_spatial_sequencing/{homo_sapiens,mus_musculus}``, holding Visium
(spots, ~55 um, several cells each) and Slide-seqV2 (beads, ~10 um) datasets.
Each dataset is one "scene" under ``spatial/<dataset_id>`` with the spot
centres in ``obsl/loc`` (pixel coordinates) and, for Visium, an H&E image in
``img/<name>/hires``. The spatial X only has a ``raw`` (count) layer.

Two steps, both narrow (4 genes only):

  scan    For every spatial dataset, read X for the 4 receptor columns and
          summarise per dataset: spots, spots that could have counted the gene
          (dataset listed it, per the presence matrix), spots with count > 0,
          summed counts and summed CP10K. Visium spots outside the tissue
          (in_tissue == 0) are dropped; Slide-seq has no in_tissue flag so
          every bead is kept.
          -> results/census_spatial/datasets.parquet, overview.json

  export  For a hand-picked list of sections (SECTIONS below), write every
          spot's coordinate plus raw count and CP10K for the 4 genes, and a
          small downscaled H&E JPEG for Visium, into the site's data/ folder.
          -> <site>/data/spatial/*.json, *.jpg, and spatial_index.json

CP10K = 10,000 * count / total counts of that spot (obs.raw_sum), which is
what the Census "normalized" layer is for the non-spatial collection.

Run from this directory (tiledbsoma needs Python 3.12):
  uv run --no-project --python 3.12 --with tiledbsoma --with pandas --with pyarrow \
      --with pillow python census_spatial_adora.py scan
  uv run ... python census_spatial_adora.py export
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tiledbsoma as soma

URI = os.environ.get("CENSUS_URI", "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma")
OUT = Path("results/census_spatial")
SITE = Path(__file__).resolve().parents[4] / "blogs/a1-receptor-pathway/scgsa-site/data"
GENES = ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]
ORGS = ["homo_sapiens", "mus_musculus"]
OBS_COLS = ["soma_joinid", "dataset_id", "assay", "tissue", "tissue_general", "disease",
            "development_stage", "donor_id", "in_tissue", "raw_sum", "nnz", "cell_type"]


def genes_of(exp):
    var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name", "feature_id"]).concat().to_pandas()
    var["up"] = var.feature_name.str.upper()
    sel = var[var.up.isin(GENES)].set_index("up").loc[GENES]
    return sel  # index = upper-case symbol, ordered as GENES


def presence(c, org, exp, var_ids):
    """dataset_id -> bool[4]: did the dataset list each gene at all."""
    ds = c["census_info"]["datasets"].read(column_names=["soma_joinid", "dataset_id"]).concat().to_pandas()
    pm = exp.ms["RNA"]["feature_dataset_presence_matrix"].read((slice(None), var_ids.tolist())).tables().concat().to_pandas()
    col = {v: i for i, v in enumerate(var_ids)}
    j2d = dict(zip(ds.soma_joinid, ds.dataset_id))
    out = {}
    for r in pm.itertuples():
        out.setdefault(j2d[r.soma_dim_0], np.zeros(len(GENES), bool))[col[r.soma_dim_1]] = bool(r.soma_data)
    return out


def read_counts(exp, joinids, var_ids):
    """Dense (len(joinids), 4) raw counts for the 4 genes; only nonzeros are stored, so zeros fill the rest."""
    pos = pd.Series(np.arange(len(joinids)), index=joinids)
    vpos = {v: i for i, v in enumerate(var_ids)}
    M = np.zeros((len(joinids), len(var_ids)), np.float32)
    for t in exp.ms["RNA"].X["raw"].read((joinids, var_ids.tolist())).tables():
        d = t.to_pandas()
        M[pos.loc[d.soma_dim_0].to_numpy(), d.soma_dim_1.map(vpos).to_numpy()] = d.soma_data.to_numpy()
    return M


def scan():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    with soma.open(URI) as c:
        dsinfo = c["census_info"]["datasets"].read(column_names=["dataset_id", "collection_name", "dataset_title", "citation"]).concat().to_pandas().set_index("dataset_id")
        for org in ORGS:
            exp = c["census_spatial_sequencing"][org]
            sel = genes_of(exp); vids = sel.soma_joinid.to_numpy()
            pres = presence(c, org, exp, vids)
            obs = exp.obs.read(column_names=OBS_COLS).concat().to_pandas()
            for k in ["dataset_id", "assay", "tissue", "tissue_general", "disease", "development_stage", "cell_type"]:
                obs[k] = obs[k].astype(str)
            n_all = obs.groupby("dataset_id").size()
            keep = obs.in_tissue.isna() | (obs.in_tissue == 1)
            obs = obs[keep].reset_index(drop=True)
            print(org, "spots kept", len(obs), "of", int(n_all.sum()), flush=True)
            M = read_counts(exp, obs.soma_joinid.to_numpy(), vids)
            cp = 1e4 * M / np.maximum(obs.raw_sum.to_numpy(), 1)[:, None]
            obs["_i"] = np.arange(len(obs))
            for dsid, g in obs.groupby("dataset_id"):
                ii = g._i.to_numpy(); p = pres.get(dsid, np.zeros(4, bool))
                r = dict(org=org, dataset_id=dsid, assay=g.assay.iat[0], tissue=g.tissue.iat[0],
                         tissue_general=g.tissue_general.iat[0], disease=g.disease.iat[0],
                         stage=g.development_stage.iat[0], donor=str(g.donor_id.iat[0]),
                         n_spots_all=int(n_all[dsid]), n_spots=len(ii),
                         med_umi=float(np.median(g.raw_sum)), med_genes=float(np.median(g.nnz)),
                         collection=dsinfo.collection_name.get(dsid, ""), title=dsinfo.dataset_title.get(dsid, ""))
                for j, gene in enumerate(GENES):
                    r[f"{gene}_listed"] = bool(p[j])
                    r[f"{gene}_meas"] = len(ii) if p[j] else 0
                    r[f"{gene}_expr"] = int((M[ii, j] > 0).sum()) if p[j] else 0
                    r[f"{gene}_counts"] = float(M[ii, j].sum()) if p[j] else 0.0
                    r[f"{gene}_cp10k_sum"] = float(cp[ii, j].sum()) if p[j] else 0.0
                rows.append(r)
    df = pd.DataFrame(rows)
    df.to_parquet(OUT / "datasets.parquet")
    print(df.groupby(["org", "assay"]).agg(n=("n_spots", "size"), spots=("n_spots", "sum")))
    for gene in GENES:
        print(gene, df.groupby(["org", "tissue_general"]).apply(
            lambda d: d[f"{gene}_expr"].sum() / max(d[f"{gene}_meas"].sum(), 1)).sort_values().tail(8).round(4).to_dict())


# ---------------------------------------------------------------- export
# Picked from results/census_spatial/datasets.parquet; see the site page for why.
SECTIONS: list[dict] = json.loads(Path(__file__).with_name("census_spatial_sections.json").read_text()) \
    if Path(__file__).with_name("census_spatial_sections.json").exists() else []


def hires_jpeg(scene, out_path, max_px=900):
    """Save a downscaled copy of the Visium H&E 'hires' image; return the
    affine that maps spot (x, y) pixel coordinates into the saved image."""
    from PIL import Image
    imgs = scene["img"]; name = list(imgs.keys())[0]; im = imgs[name]
    # The image's schema has data_axis_permutation [2, 1, 0]: the stored array
    # is (channel, x, y), so reversing the axes gives the usual (y, x, channel).
    # Checked by overlaying in-tissue spots on the H&E: with this order they sit
    # on the tissue; read naively as (C, Y, X) they land on blank glass.
    arr = im["hires"].read().to_numpy()
    perm = json.loads(im.metadata["soma_multiscale_image_schema"])["data_axis_permutation"]
    assert perm == [2, 1, 0], perm
    arr = np.transpose(arr, (2, 1, 0)).astype(np.uint8)
    # loc (full-resolution pixels) -> hires pixels
    tr = scene.get_transform_to_multiscale_image(name, subcollection="img", level="hires")
    sx, sy = None, None
    try:
        A = np.asarray(tr.augmented_matrix, float)
        sx, sy = A[0, 0], A[1, 1]
    except Exception:
        sx = sy = float(getattr(tr, "scale_factors", [1, 1])[0])
    img = Image.fromarray(arr)
    f = min(1.0, max_px / max(img.size))
    if f < 1:
        img = img.resize((round(img.size[0] * f), round(img.size[1] * f)), Image.LANCZOS)
    img.save(out_path, quality=78, optimize=True)
    return dict(w=img.size[0], h=img.size[1], sx=float(sx) * f, sy=float(sy) * f, src=name, transform=str(tr))


def export():
    if not SECTIONS:
        sys.exit("write census_spatial_sections.json first")
    od = SITE / "spatial"; od.mkdir(parents=True, exist_ok=True)
    meta = pd.read_parquet(OUT / "datasets.parquet").set_index("dataset_id")
    index = []
    with soma.open(URI) as c:
        for s in SECTIONS:
            org, dsid = s["org"], s["dataset_id"]
            exp = c["census_spatial_sequencing"][org]
            sel = genes_of(exp); vids = sel.soma_joinid.to_numpy()
            pres = presence(c, org, exp, vids)[dsid]
            scene = exp.spatial[dsid]
            loc = scene["obsl"]["loc"].read().concat().to_pandas()
            obs = exp.obs.read(value_filter=f"dataset_id == '{dsid}'",
                               column_names=["soma_joinid", "in_tissue", "raw_sum", "nnz", "cell_type"]).concat().to_pandas()
            obs = obs[obs.in_tissue.isna() | (obs.in_tissue == 1)]
            d = obs.merge(loc, on="soma_joinid", how="inner").sort_values("soma_joinid").reset_index(drop=True)
            M = read_counts(exp, d.soma_joinid.to_numpy(), vids)
            cp = 1e4 * M / np.maximum(d.raw_sum.to_numpy(), 1)[:, None]
            m = meta.loc[dsid]
            rec = dict(key=s["key"], label=s["label"], why=s["why"], org=org, dataset_id=dsid,
                       assay=m.assay, tissue=m.tissue, tissue_general=m.tissue_general, disease=m.disease,
                       stage=m.stage, donor=m.donor, collection=m.collection, title=m.title,
                       n_spots=len(d), n_spots_all=int(m.n_spots_all), med_umi=float(np.median(d.raw_sum)),
                       med_genes=float(np.median(d.nnz)), genes=GENES, listed=[bool(x) for x in pres],
                       radius=float(scene["obsl"]["loc"].metadata.get("soma_geometry", 0)),
                       x=[round(float(v), 1) for v in d.x], y=[round(float(v), 1) for v in d.y],
                       umi=[int(v) for v in d.raw_sum],
                       counts=[[int(v) for v in M[:, j]] for j in range(4)],
                       cp10k=[[round(float(v), 2) for v in cp[:, j]] for j in range(4)])
            if s.get("cell_type"):
                cats = d.cell_type.astype(str)
                u = sorted(cats.unique()); rec["cell_types"] = u
                rec["ct"] = [u.index(v) for v in cats]
            if "img" in scene and s.get("image", True):
                rec["image"] = hires_jpeg(scene, od / f"{s['key']}.jpg")
                rec["image"]["file"] = f"data/spatial/{s['key']}.jpg"
            (od / f"{s['key']}.json").write_text(json.dumps(rec, separators=(",", ":")))
            print(s["key"], len(d), "spots", [int((M[:, j] > 0).sum()) for j in range(4)], flush=True)
            index.append({k: rec[k] for k in ["key", "label", "why", "org", "assay", "tissue", "disease", "stage",
                                              "n_spots", "collection", "title", "listed", "dataset_id"]})
    (SITE / "spatial_index.json").write_text(json.dumps(index, indent=1))


def overview():
    """Small summary JSON for the page: one row per dataset (no spot data)."""
    df = pd.read_parquet(OUT / "datasets.parquet")
    cols = ["org", "dataset_id", "assay", "tissue", "tissue_general", "disease", "n_spots", "n_spots_all",
            "med_umi", "med_genes", "collection", "title"]
    rows = []
    for r in df.itertuples():
        d = {k: getattr(r, k) for k in cols}
        d["meas"] = [int(getattr(r, f"{g}_meas")) for g in GENES]
        d["expr"] = [int(getattr(r, f"{g}_expr")) for g in GENES]
        d["cp"] = [round(float(getattr(r, f"{g}_cp10k_sum")), 3) for g in GENES]
        rows.append(d)
    (SITE / "spatial_overview.json").write_text(json.dumps({"genes": GENES, "datasets": rows}, separators=(",", ":")))
    print("wrote", len(rows), "datasets")


if __name__ == "__main__":
    {"scan": scan, "export": export, "overview": overview}[sys.argv[1]]()
