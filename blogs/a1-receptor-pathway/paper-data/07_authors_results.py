"""Collect the authors' published recovery rates for the playground.

Source: the paper's Zenodo archive (doi:10.5281/zenodo.13620619, CC-BY-4.0),
Results/{RWSD,SSD}/Sample size N/Scenario X[ variant]/Two side dat matrix.csv.
Each file: rows = method x replicate ("AUCell", "AUCell1", ...), columns =
"GS-<size>-<noise>-", values = % of the 100 gene sets identified (two-sided
Wilcoxon, BH FDR < 0.05). We keep mean and SD over replicates.
Output: scgsa-site/data/authors_results.json
  {data: {cells: {scenario: {size: {noise: {method: [mean, sd, n]}}}}}}
"""
import csv, json, re, statistics, sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else
            "/mnt/bulk/dl_bio/scps_zenodo/RuoqiaoWrqqq-scPS-bfc34ac/Results")
out = Path(__file__).resolve().parents[1] / "scgsa-site" / "data" / "authors_results.json"
methods = ["AddModuleScore", "AUCell", "UCell", "ssGSEA", "JASMINE", "SCSE", "scPS"]
res = {}
for f in sorted(root.glob("*/Sample size */Scenario */Two side dat matrix.csv")):
    data = f.parts[-4].lower()
    cells = f.parts[-3].split()[-1]
    scen = f.parts[-2].removeprefix("Scenario ").strip()
    rows = list(csv.reader(f.open()))
    head = rows[0][1:]
    by = {}
    for r in rows[1:]:
        m = re.sub(r"\d+$", "", r[0])
        by.setdefault(m, []).append([float(x) if x not in ("", "NA") else None for x in r[1:]])
    node = res.setdefault(data, {}).setdefault(cells, {}).setdefault(scen, {})
    for j, col in enumerate(head):
        mm = re.match(r"GS-(\d+)-([\d.]+)-", col)
        if not mm:
            continue
        size, noise = mm.group(1), str(float(mm.group(2)))
        for m in methods:
            vals = [row[j] for row in by.get(m, []) if row[j] is not None]
            if vals:
                sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
                node.setdefault(size, {}).setdefault(noise, {})[m] = [round(statistics.mean(vals), 1), round(sd, 1), len(vals)]
out.write_text(json.dumps({"source": "doi:10.5281/zenodo.13620619 (CC-BY-4.0), Results/*/Two side dat matrix.csv",
                           "unit": "% of 100 gene sets identified; mean, SD, replicates", "data": res}, separators=(",", ":")))
for d, v in res.items():
    print(d, {c: sorted(s) for c, s in v.items()})
print(f"{out} {out.stat().st_size/1e3:.0f} kB")
