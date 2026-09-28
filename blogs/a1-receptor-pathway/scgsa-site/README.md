# scGSA from first principles (reading companion)

A static site for reading Wang & Thakar 2024 (`../lqae124.pdf`) and working
through each of its seven single-cell gene set scoring methods from first
principles. Pedagogical and disposable — delete the folder when done.

**Live (tailnet only):** <https://acefsan-ubuntu.taila37484.ts.net:3044>

## Layout

| Path | What |
|---|---|
| `index.html` | overview and reading path |
| `paper.html` | PDF reader (pdf.js from cdnjs) for `lqae124.pdf` |
| `debt.html` | the post's three cognitive-debt questions, with folded answers |
| `basics.html`, `scenarios.html`, `compare.html` | foundations, the four scenarios, side-by-side tables |
| `methods/*.html` | one page per method |
| `playground.html` | in-browser simulator of the four scenarios |
| `ours.html` | what the paper implies for the A1 analysis |
| `static/scgsa.js` | simulator + JS versions of all seven methods (runs in node too) |
| `validate/` | checks of `scgsa.js` against the real R packages |

No build step. KaTeX and pdf.js load from cdnjs, so pages need internet
access for maths and the PDF viewer; everything else is local.

## Serving

Same pattern as `adenosine-viz`: an `nginx:alpine` container in
`~/services/docker-compose.yml` (service `scgsa-site`) mounts this folder
read-only on `127.0.0.1:3044`, plus `../lqae124.pdf` as `/lqae124.pdf`, and
`tailscale serve --https=3044` fronts it. The mount points at the main
checkout (`~/src/dl_bio/blogs/...`), so edits there show up on reload.

To take it down: `docker compose -f ~/services/docker-compose.yml rm -sf scgsa-site`
and `tailscale serve --https=3044 off`.

## Fidelity of the JS methods

`validate/export.js` simulates a matrix and scores one gene set in JS;
`validate/validate.R` scores the same inputs with Seurat 5.5.1, AUCell 1.28.0,
UCell 2.10.1, GSVA 2.0.0, the authors' JASMINE source and an scPS pipeline
following `scPS_2025.R`. Results, all four scenarios:

| Method | Agreement with R |
|---|---|
| UCell, ssGSEA, JASMINE, SCSE | exact (r = 1.000000) |
| AddModuleScore | r ≈ 0.91–0.99 — random control genes differ between RNGs |
| AUCell | r ≈ 0.97–0.996 — random tie-breaking differs between RNGs |
| scPS | r ≥ 0.993 once each PC's sign is aligned; unaligned it ranged −0.39 to 0.97 |

Run: `node validate/export.js <dir> <scenario> [log|counts]` then
`Rscript validate/validate.R <dir> <jasmine_functions.R>` in the `scgsa`
conda env built by `projects/caffeine/lab/001_adora_expression/scgsa_setup_env.sh`
(`jasmine_functions.R` = the authors' JASMINE file with its top-level example
call removed). `validate/sweep.js` and `validate/regime.js` print the toy's
recovery rates next to the paper's Table 1.

## Things found while building this

- The published scPS code differs from the paper's Equation 1: it takes a
  square root, subtracts one global minimum rather than one per PC, and with
  Seurat's default `approx = TRUE` normalises "variance explained" over the
  top 10 PCs only (`approx = FALSE` normalises over all of them).
- `scPS_2024.R` assigns every cell the same PC component when one PC explains
  more than 50% of the variance; `scPS_2025.R` fixes it.
- scPS depends on each PC's arbitrary sign.
- HBCA superclusters differ 4.2x in genes detected per cell — the scenario 3
  confound — and the extremes line up with the extremes of our A1 rankings.
