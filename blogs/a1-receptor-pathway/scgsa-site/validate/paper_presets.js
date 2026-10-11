// Recovery rates on the paper's real cells (naive-B RWSD) for the playground's
// "Things to try" settings, averaged over a few simulation seeds.
// Usage: node validate/paper_presets.js [naive_b|ssd]   (RWSD threshold 0.75, SSD 0.9)
const fs = require("fs"), path = require("path");
const S = require("../static/scgsa.js");
const site = path.join(__dirname, "..");
const which = process.argv[2] || "naive_b";
const threshold = which === "ssd" ? 0.9 : 0.75;
const meta = JSON.parse(fs.readFileSync(path.join(site, `data/${which}/meta.json`)));
const pool = S.decodePool(fs.readFileSync(path.join(site, `data/${which}/counts.bin`)).buffer, meta);
const cases = [
  ["S1 signal, size 100, 200 cells", { scenario: 1, cellsPerGroup: 200 }, { size: 100, noise: 0 }],
  ["S1 null (noise 100%), size 100, 200 cells", { scenario: 1, cellsPerGroup: 200 }, { size: 100, noise: 1 }],
  ["S2 signal, size 100, 20 cells", { scenario: 2, cellsPerGroup: 20 }, { size: 100, noise: 0 }],
  ["S2 signal, size 100, 200 cells", { scenario: 2, cellsPerGroup: 200 }, { size: 100, noise: 0 }],
  ["S3 null (noise 0), size 100, 200 cells", { scenario: 3, cellsPerGroup: 200 }, { size: 100, noise: 0 }],
  ["S3 null (noise 0), size 25, 200 cells", { scenario: 3, cellsPerGroup: 200 }, { size: 25, noise: 0 }],
  ["S4 null (noise 0), size 100, 200 cells", { scenario: 4, cellsPerGroup: 200 }, { size: 100, noise: 0 }],
];
const seeds = [7, 8, 9];
for (const [label, sopt, spec] of cases) {
  const acc = {};
  for (const seed of seeds) {
    const sim = S.simulateRWSD(pool, Object.assign({ seed, threshold }, sopt));
    const prep = S.prepare(sim);
    const r = S.recovery(sim, prep, spec, 30, 100 + seed);
    for (const [k, v] of Object.entries(r)) acc[k] = (acc[k] || 0) + v / seeds.length;
  }
  console.log(label.padEnd(44), Object.entries(acc).map(([k, v]) => `${k} ${Math.round(100 * v)}%`).join("  "));
}
