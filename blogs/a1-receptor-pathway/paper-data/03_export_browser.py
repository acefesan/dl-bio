"""Pack raw counts for the browser: the naive-B cells (02_export_naive_b.R),
or with argument "ssd" the Splatter pool (05_splatter_pool.R).

Output, under scgsa-site/data/naive_b/ (or data/ssd/):
  counts.bin  per cell: varint(n nonzero), then per nonzero gene
              varint(gene index delta) + varint(count). All cells, all genes.
  meta.json   gene names, per-cell donor, cell count, byte size.
The browser LogNormalizes these itself, exactly as Seurat does
(log1p(count / cell total * 1e4)), so nothing is lost but disk space.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.io import mmread

which = sys.argv[1] if len(sys.argv) > 1 else "naive_b"
src = Path("/mnt/bulk/dl_bio/gse164381")
out = Path(__file__).resolve().parents[1] / "scgsa-site" / "data" / which
out.mkdir(parents=True, exist_ok=True)

X = mmread(src / f"{which}_counts.mtx").tocsc()  # genes x cells
genes = (src / f"{which}_genes.txt").read_text().split()
cells = (src / f"{which}_cells.txt").read_text().split()
assert X.shape == (len(genes), len(cells))

def varint(v, buf):
    while v >= 0x80:
        buf.append((v & 0x7F) | 0x80)
        v >>= 7
    buf.append(v)

buf = bytearray()
for c in range(X.shape[1]):
    lo, hi = X.indptr[c], X.indptr[c + 1]
    idx = X.indices[lo:hi]
    val = X.data[lo:hi].astype(np.int64)
    order = np.argsort(idx)
    idx, val = idx[order], val[order]
    varint(len(idx), buf)
    prev = -1
    for g, v in zip(idx.tolist(), val.tolist()):
        varint(g - prev - 1, buf)
        varint(int(v), buf)
        prev = g

(out / "counts.bin").write_bytes(bytes(buf))
donor_codes = {"GSM5008746": "healthy donor 1", "GSM5008747": "healthy donor 2",
               "GSM5008748": "healthy donor 3", "GSM5008749": "healthy donor 4",
               "SSD": "Splatter-simulated"}
donors = sorted({c.split("_")[0] for c in cells})
meta = {
    "source": ("Splatter pool fitted to the naive-B cells (paper-data/05)" if which == "ssd"
               else "GSE164381, 4 healthy donors, naive-B clusters (see paper-data/01-02)"),
    "nCells": len(cells), "nGenes": len(genes), "nnz": int(X.nnz),
    "bytes": len(buf), "genes": genes,
    "donors": [donor_codes.get(d, d) for d in donors],
    "cellDonor": [donors.index(c.split("_")[0]) for c in cells],
    "barcodes": [c.split("_", 1)[1] for c in cells],
}
(out / "meta.json").write_text(json.dumps(meta, separators=(",", ":")))
print(f"{len(cells)} cells, {len(genes)} genes, nnz {X.nnz}, counts.bin {len(buf)/1e6:.1f} MB, "
      f"meta {(out / 'meta.json').stat().st_size/1e6:.2f} MB")
