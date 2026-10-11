// Replays the random draws of paper-data/04_authors_rwsd_reference.R (the
// authors' RWSD script run in R) through SCGSA.simulateRWSD and checks the
// browser simulator builds the same matrix, gene for gene and cell for cell.
// Usage: node validate/rwsd_check.js /mnt/bulk/dl_bio/gse164381
const fs = require("fs"), path = require("path");
const S = require("../static/scgsa.js");
const src = process.argv[2] || "/mnt/bulk/dl_bio/gse164381";
const site = path.join(__dirname, "..");
const meta = JSON.parse(fs.readFileSync(path.join(site, "data/naive_b/meta.json")));
const pool = S.decodePool(fs.readFileSync(path.join(site, "data/naive_b/counts.bin")).buffer, meta);
const names = fs.readFileSync(path.join(src, "naive_b_cells.txt"), "utf8").split(/\s+/).filter(Boolean);
const idx = new Map(names.map((c, i) => [c, i]));
const gidx = new Map(meta.genes.map((g, i) => [g, i]));

function readMM(f) {
  const lines = fs.readFileSync(f, "utf8").split("\n").filter(l => l && l[0] !== "%");
  const [nr, nc] = lines[0].split(/\s+/).map(Number);
  const m = new Map();
  for (const l of lines.slice(1)) { const [r, c, v] = l.split(/\s+/); m.set((+r - 1) * nc + (+c - 1), +v); }
  return { nr, nc, get: (r, c) => m.get(r * nc + c) || 0 };
}

function compare(label, sim, ref, rGenes, rCells, jsCellOf) {
  const jsGene = new Map(sim.geneNames.map((g, j) => [g.replace(" (added)", "add"), j]));
  let maxErr = 0, checked = 0, missing = 0;
  rGenes.forEach((g, r) => {
    const j = jsGene.get(g);
    if (j === undefined) { missing++; return; }
    rCells.forEach((c, k) => {
      const v = sim.data[jsCellOf(k) * sim.G + j], w = ref.get(r, k);
      maxErr = Math.max(maxErr, Math.abs(v - w) / Math.max(1, Math.abs(w)));
      checked++;
    });
  });
  console.log(`${label}: genes R ${rGenes.length} / JS ${sim.G}, missing ${missing}, ` +
    `${checked.toLocaleString()} values compared, max relative error ${maxErr.toExponential(2)}`);
  return missing === 0 && rGenes.length === sim.G && maxErr < 1e-5;
}

const ref = path.join(src, "authors_reference");
let ok = true;
{ // Scenario 1
  const j = JSON.parse(fs.readFileSync(path.join(ref, "s1.json")));
  const n = j.cells.length / 2;
  const ctl = j.cells.filter((_, k) => j.cond[k] === "Condition_1");
  const trt = j.cells.filter((_, k) => j.cond[k] === "Condition_2");
  const order = ctl.concat(trt).map(c => idx.get(c));
  const sim = S.simulateRWSD(pool, { scenario: 1, cellsPerGroup: n, cells: order });
  const sigR = new Set(j.signal), sigJ = new Set(sim.signal.map(x => sim.geneNames[x]));
  const same = sigR.size === sigJ.size && [...sigR].every(g => sigJ.has(g));
  console.log(`S1: signal genes R ${sigR.size} / JS ${sigJ.size}, identical: ${same}`);
  const pos = new Map(ctl.concat(trt).map((c, i) => [c, i]));
  ok = compare("S1", sim, readMM(path.join(ref, "s1_dat.mtx")), j.genes, j.cells, k => pos.get(j.cells[k])) && same && ok;
}
{ // Scenario 3, boost as the script is written (Signal_Level = 1.2)
  const j = JSON.parse(fs.readFileSync(path.join(ref, "s3.json")));
  const n = j.cells.length / 2;
  const sim = S.simulateRWSD(pool, { scenario: 3, cellsPerGroup: n, boostNull: true,
    cells: j.cells.map(c => idx.get(c)), cells2: j.cells2.map(c => idx.get(c)),
    added: j.added.map(g => gidx.get(g)) });
  ok = compare("S3", sim, readMM(path.join(ref, "s3_dat.mtx")), j.genes, j.cells, k => k) && ok;
}
console.log(ok ? "PASS: browser RWSD matches the authors' R script" : "FAIL");
process.exit(ok ? 0 : 1);
