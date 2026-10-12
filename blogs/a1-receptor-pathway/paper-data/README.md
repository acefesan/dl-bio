# Paper data for the scGSA playground

Rebuilds the inputs of Wang & Thakar 2024 (NARGAB, lqae124) from GSE164381, so
the reading site's playground can run the paper's scenarios on the same kind
of cells the paper used. Large files live in `/mnt/bulk/dl_bio/gse164381/`;
the packed browser data is in `../scgsa-site/data/{naive_b,ssd}/`.

| Step | What it does |
|---|---|
| `01_naive_b.R` | 4 healthy-donor 10x matrices → QC (<15% mito, <4000 genes) → Seurat clustering → naive-vs-memory marker table |
| `02_export_naive_b.R` | keeps naive-B clusters 0, 1, 3, 4, 7 (IGHD/TCL1A high; CD27/IGHG1/CD3E/JCHAIN low): 9,528 cells; saves `naive_b.rds` + counts |
| `03_export_browser.py [naive_b\|ssd]` | packs raw counts (varint sparse) for the browser |
| `04_authors_rwsd_reference.R` | runs the authors' RWSD script (S1, S3) on `naive_b.rds`, saves its draws and matrices; `../scgsa-site/validate/rwsd_check.js` replays them in the browser simulator and must match exactly |
| `05_splatter_pool.R` | the authors' SSD: `splatEstimate` on the naive-B counts, then a 4,000-cell Splatter pool |
| `06_provenance.R` | records each build step's numbers and the cluster marker table → `scgsa-site/data/provenance.json` (the playground's "Build the dataset" node) |
| `07_authors_results.py` | the authors' published recovery rates from their Zenodo archive → `scgsa-site/data/authors_results.json` (shown beside ours in the recovery and sweep nodes) |
| `08_umap.R` | UMAP of the QC-passed healthy-donor cells (5,000-cell sample: cluster, donor, marker genes) → `data/provenance_umap.json` |
| `09_qc_stats.R` | what each QC rule removed, plus genes-per-cell and %-mitochondrial histograms → `data/provenance_qc.json` |

Runs with `~/.local/share/mamba/envs/scgsa/bin/Rscript` (R 4.4, Seurat 5.5.1,
splatter 1.30.0) and `uv run --python 3.12 --with scipy --with numpy`.

## Source of truth: the authors' code

[Thakar-Lab/scPS](https://github.com/Thakar-Lab/scPS), "Comparison analysis/".
Where it and the paper's text disagree, the playground follows the code:

- **Input:** `Naive Normal B cell` — naive B cells from the healthy donors only.
  The authors' annotation is not published; ours is a marker-based
  reconstruction (see 01/02).
- **Dense genes (S1/S3):** `rows_zeros < Threshold` with Threshold 0.75 (RWSD)
  / 0.9 (SSD), i.e. detected in **more than 25% / 10%** of cells. The text says
  "at least 75% / 90%". On our cells the code's rule gives ~485 genes with a
  median detection of 43% (paper: "at least 500", 44.5%); the text's rule
  would give ~105.
- **"+20%":** `dat[Signal, treatment] * 1.2` on log-normalised values.
- **S3/S4:** 250 genes from a second, independent cell draw, zero in control,
  appended as new genes and added to the non-signal pool; their share in a
  gene set is left to chance. The script as written still multiplies the dense
  genes by 1.2; run that way scPS fires on every S3 set, contradicting the
  paper, so the playground defaults to no boost (`boostNull`), matching the
  authors' "without signal" result folders.
- **Test:** two-sided Wilcoxon (`row_wilcoxon_twosample`, continuity
  correction), BH per replicate, FDR < 0.05.
- **UCell:** irGSEA's default `maxRank = 1500`.

## Known differences on our reconstruction

ssGSEA almost never fires on these cells, including in S3, where the paper
says every method but scPS did.

SSD (Splatter fitted to our naive-B cells) comes out sparser than the
paper's: the median S1 signal gene is detected in ~15% of cells, against the
paper's 32%. Its results match the paper's SSD in direction (UCell and
ssGSEA weak; the S3 trap fires for every method except scPS and UCell) but
its null false-positive rates are lower (JASMINE ~58%, AUCell ~8%, SCSE ~11%
at 100% noise). Splatter gene names (Gene1, ...) are synthetic, not real genes.

## Is the authors' processed dataset published anywhere?

No. Checked on 2026-10-11:

- **The paper's data statement** lists only GEO GSE164381 (and GSE198339 for
  the HIV example), plus code at GitHub and Zenodo.
- **The supplement** (PMC11420841, one PDF) holds figures only. Its captions
  do settle one question: "For SSD and RWSD, no signal was assigned to the
  signal genes in scenarios 3 … and 4". So S3/S4 used no ×1.2 boost, which is
  the playground's default.
- **The Zenodo archive** (doi:10.5281/zenodo.13620619, CC-BY-4.0, 123 MB) is
  a snapshot of the GitHub repo (RuoqiaoWrqqq/scPS v1.0.0): code, the HIV
  example's score tables, and the simulation **results**
  (`Results/{RWSD,SSD}/Sample size N/Scenario X/Two side dat matrix.csv`).
  There are no cell labels and no `Naive Normal B cell.RData`.
- **Reference 25** (Hoehn et al. 2021, J Immunol) has no processed object on
  Zenodo either.

So the naive-B subset stays a reconstruction. Against the published results
(07), our RWSD runs match closely for every method except ssGSEA.
