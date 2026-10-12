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
          (superseded by export_all; kept for reference)

  export_all
          Every one of the 540 spatial datasets, one compact binary file per
          section (format in _write_bin's docstring), a ~1000 px H&E JPEG per
          Visium section cropped to the tissue, and one index JSON.
          Resumable: a dataset whose file and sidecar already exist is
          skipped. Progress goes to stdout and results/census_spatial/export_all.log.
          -> <site>/data/spatial_all/{index.json, s/<dataset_id>.bin, img/<dataset_id>.jpg}

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


# ---------------------------------------------------------------- export_all
ALL = SITE / "spatial_all"
SIDE = OUT / "all_meta"          # per-dataset metadata sidecars (resume + index)
MAGIC = b"ADSP"
FORMAT_VERSION = 1


def _log(msg, fh=[None]):
    import time
    if fh[0] is None:
        OUT.mkdir(parents=True, exist_ok=True)
        fh[0] = open(OUT / "export_all.log", "a")
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True); fh[0].write(line + "\n"); fh[0].flush()


def _hires_crop(scene, bbox, out_path, max_px=1000, margin=0.04):
    """Visium H&E 'hires' image, cropped to the in-tissue spots' bounding box
    (plus a small margin) and downscaled to <= max_px on the long side.
    Returns the rectangle of *spot coordinates* the JPEG covers, so the page
    can draw it with one drawImage call."""
    from PIL import Image
    imgs = scene["img"]; name = list(imgs.keys())[0]; im = imgs[name]
    arr = im["hires"].read().to_numpy()
    perm = json.loads(im.metadata["soma_multiscale_image_schema"])["data_axis_permutation"]
    assert perm == [2, 1, 0], perm          # stored (channel, x, y): see hires_jpeg()
    arr = np.transpose(arr, (2, 1, 0)).astype(np.uint8)   # -> (y, x, channel)
    tr = scene.get_transform_to_multiscale_image(name, subcollection="img", level="hires")
    A = np.asarray(tr.augmented_matrix, float)
    # Spot -> hires pixel is a pure isotropic scale in every Visium scene checked;
    # assert it so a scene with a translation or anisotropy fails loudly instead
    # of being drawn misaligned.
    assert abs(A[0, 1]) < 1e-9 and abs(A[1, 0]) < 1e-9 and abs(A[0, 2]) < 1e-9 and abs(A[1, 2]) < 1e-9, A
    assert abs(A[0, 0] - A[1, 1]) <= 1e-6 * abs(A[0, 0]), A
    s = A[0, 0]
    H, W = arr.shape[:2]
    x0, y0, x1, y1 = bbox
    m = margin * max(x1 - x0, y1 - y0)
    cx0 = max(0, int(np.floor((x0 - m) * s))); cy0 = max(0, int(np.floor((y0 - m) * s)))
    cx1 = min(W, int(np.ceil((x1 + m) * s))); cy1 = min(H, int(np.ceil((y1 + m) * s)))
    img = Image.fromarray(arr[cy0:cy1, cx0:cx1])
    f = min(1.0, max_px / max(img.size))
    if f < 1:
        img = img.resize((max(1, round(img.size[0] * f)), max(1, round(img.size[1] * f))), Image.LANCZOS)
    img.save(out_path, quality=75, optimize=True, progressive=True)
    return dict(w=img.size[0], h=img.size[1], x0=cx0 / s, y0=cy0 / s, x1=cx1 / s, y1=cy1 / s,
                hires_w=W, hires_h=H, hires_scale=s, src=name)


def _write_bin(path, header, xq, yq, umi, sparse):
    """Binary section file, little-endian, every block 4-byte aligned:

      0   'ADSP'                     magic
      4   uint32 version (1)
      8   uint32 header length L (bytes of UTF-8 JSON, before padding)
      12  header JSON, zero-padded to a multiple of 4
      ..  x    uint16[n]  quantised: x = x0 + q * (x1 - x0) / 65535  (pad to 4)
      ..  y    uint16[n]  same with y0, y1                          (pad to 4)
      ..  umi  uint16[n] or uint32[n] (header.umi_bytes): spot total counts (raw_sum)
      ..  per gene g, for the nnz[g] spots with a count > 0:
            idx   uint32[nnz]   spot index (0..n-1, ascending)
            count uint16[nnz]   raw count                            (pad to 4)
    CP10K is not stored: the page computes 1e4 * count / max(umi, 1), the same
    expression the exporter used everywhere else. A gene the dataset did not
    list has nnz 0 and header.listed[g] false ("not measured").
    """
    def pad(b):
        return b + b"\0" * (-len(b) % 4)
    hj = json.dumps(header, separators=(",", ":")).encode()
    parts = [MAGIC, np.array([FORMAT_VERSION, len(hj)], "<u4").tobytes(), pad(hj),
             pad(xq.astype("<u2").tobytes()), pad(yq.astype("<u2").tobytes()),
             pad(umi.astype("<u2" if header["umi_bytes"] == 2 else "<u4").tobytes())]
    for idx, cnt in sparse:
        parts += [idx.astype("<u4").tobytes(), pad(cnt.astype("<u2").tobytes())]
    tmp = path.with_suffix(".tmp")
    tmp.write_bytes(b"".join(parts)); tmp.replace(path)


def export_all():
    import time
    t_start = time.time()
    (ALL / "s").mkdir(parents=True, exist_ok=True); (ALL / "img").mkdir(parents=True, exist_ok=True)
    SIDE.mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[2:])        # optional: dataset ids to (re)do, for testing
    failed = []
    with soma.open(URI) as c:
        dsinfo = c["census_info"]["datasets"].read(column_names=["dataset_id", "collection_name", "dataset_title", "citation"]).concat().to_pandas().set_index("dataset_id")
        for org in ORGS:
            exp = c["census_spatial_sequencing"][org]
            sel = genes_of(exp); vids = sel.soma_joinid.to_numpy()
            pres = presence(c, org, exp, vids)
            obs = exp.obs.read(column_names=OBS_COLS).concat().to_pandas()
            for k in ["dataset_id", "assay", "tissue", "tissue_general", "disease", "development_stage", "donor_id"]:
                obs[k] = obs[k].astype(str)
            groups = obs.groupby("dataset_id", observed=True)
            dsids = sorted(groups.groups)
            _log(f"{org}: {len(dsids)} datasets, {len(obs):,} obs rows")
            for k, dsid in enumerate(dsids, 1):
                fbin, fside = ALL / "s" / f"{dsid}.bin", SIDE / f"{dsid}.json"
                if (not only and fbin.exists() and fside.exists()) or (only and dsid not in only):
                    continue
                t0 = time.time()
                try:
                    g = groups.get_group(dsid)
                    n_all = len(g)
                    g = g[g.in_tissue.isna() | (g.in_tissue == 1)]
                    scene = exp.spatial[dsid]
                    loc = scene["obsl"]["loc"].read().concat().to_pandas()
                    d = g.merge(loc, on="soma_joinid", how="inner").sort_values("soma_joinid").reset_index(drop=True)
                    assert len(d) == len(g), f"{len(g) - len(d)} in-tissue spots without a location"
                    M = read_counts(exp, d.soma_joinid.to_numpy(), vids).astype(np.int64)
                    assert (M == np.round(M)).all()
                    umi = d.raw_sum.to_numpy()
                    assert (umi == np.round(umi)).all() and umi.min() >= 0
                    umi = umi.astype(np.int64)
                    cp = 1e4 * M / np.maximum(umi, 1)[:, None]
                    p = pres.get(dsid, np.zeros(4, bool))
                    # A gene the dataset did not list must have no counts at all.
                    for j in range(4):
                        if not p[j]:
                            assert (M[:, j] == 0).all(), (dsid, GENES[j])
                    assert M.max() < 65536
                    x, y = d.x.to_numpy(float), d.y.to_numpy(float)
                    x0, x1, y0, y1 = x.min(), x.max(), y.min(), y.max()
                    q = lambda v, a, b: np.round((v - a) / max(b - a, 1e-9) * 65535).astype(np.uint16)
                    xq, yq = q(x, x0, x1), q(y, y0, y1)
                    radius = float(scene["obsl"]["loc"].metadata.get("soma_geometry", 0))
                    hdr = dict(dataset_id=dsid, org=org, n=len(d), genes=GENES, listed=[bool(v) for v in p],
                               nnz=[int((M[:, j] > 0).sum()) for j in range(4)],
                               x0=float(x0), x1=float(x1), y0=float(y0), y1=float(y1), radius=radius,
                               units=json.loads(scene.metadata.get("soma_coordinate_space", "[]")),
                               umi_bytes=2 if umi.max() < 65536 else 4)
                    img_meta = None
                    if "img" in scene:
                        img_meta = _hires_crop(scene, (x0, y0, x1, y1), ALL / "img" / f"{dsid}.jpg")
                        img_meta["file"] = f"img/{dsid}.jpg"
                        hdr["image"] = img_meta
                    sparse = [(np.nonzero(M[:, j] > 0)[0], M[M[:, j] > 0, j]) for j in range(4)]
                    _write_bin(fbin, hdr, xq, yq, umi, sparse)
                    info = dsinfo.loc[dsid] if dsid in dsinfo.index else None
                    nz = [int((M[:, j] > 0).sum()) for j in range(4)]
                    meta = dict(
                        id=dsid, org=org, assay="Visium" if g.assay.iat[0].startswith("Visium") else g.assay.iat[0],
                        tissue_general=g.tissue_general.iat[0], tissue=g.tissue.iat[0], disease=g.disease.iat[0],
                        stage=g.development_stage.iat[0], donor=str(g.donor_id.iat[0]),
                        title=str(info.dataset_title) if info is not None else "",
                        collection=str(info.collection_name) if info is not None else "",
                        n=len(d), n_all=n_all, med_umi=float(np.median(umi)), med_genes=float(np.median(d.nnz)),
                        listed=[bool(v) for v in p],
                        expr=[nz[j] if p[j] else None for j in range(4)],
                        pct=[round(100 * nz[j] / max(len(d), 1), 3) if p[j] else None for j in range(4)],
                        cp_mean=[round(float(cp[:, j].mean()), 4) if p[j] else None for j in range(4)],
                        radius=radius, img=bool(img_meta), bytes=fbin.stat().st_size,
                        img_bytes=(ALL / "img" / f"{dsid}.jpg").stat().st_size if img_meta else 0)
                    fside.write_text(json.dumps(meta))
                    _log(f"{org} {k}/{len(dsids)} {dsid} {meta['assay']} {meta['tissue']!r} n={len(d)} "
                         f"pct={meta['pct']} bin={meta['bytes']//1024}KB img={meta['img_bytes']//1024}KB {time.time()-t0:.1f}s")
                except Exception as e:  # keep going; report at the end
                    failed.append((dsid, repr(e)))
                    _log(f"FAILED {org} {dsid}: {e!r}")
    _build_index()
    _log(f"export_all done in {(time.time()-t_start)/60:.1f} min; {len(failed)} failed: {failed}")


def _build_index():
    rows = [json.loads(p.read_text()) for p in sorted(SIDE.glob("*.json"))]
    rows.sort(key=lambda r: (r["org"], r["tissue_general"], r["tissue"], r["id"]))
    feat = []
    for s in SECTIONS:
        feat.append(dict(id=s["dataset_id"], label=s["label"], why=s["why"]))
    cols = ["id", "org", "assay", "tissue_general", "tissue", "disease", "stage", "donor", "title", "collection",
            "n", "n_all", "med_umi", "med_genes", "listed", "expr", "pct", "cp_mean", "radius", "img", "bytes", "img_bytes"]
    out = dict(census="2025-11-08", genes=GENES, format_version=FORMAT_VERSION, featured=feat,
               cols=cols, rows=[[r[k] for k in cols] for r in rows])
    (ALL / "index.json").write_text(json.dumps(out, separators=(",", ":")))
    tot = sum(r["bytes"] + r["img_bytes"] for r in rows)
    _log(f"index: {len(rows)} datasets, {tot/1e6:.1f} MB of section files + images")


if __name__ == "__main__":
    {"scan": scan, "export": export, "overview": overview, "export_all": export_all,
     "index": _build_index}[sys.argv[1]]()
