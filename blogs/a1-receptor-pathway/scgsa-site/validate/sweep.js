// Does the toy simulator reproduce the paper's qualitative findings?
// Prints recovery rates (share of gene sets with BH-adjusted Wilcoxon p < 0.05).
const S = require("../static/scgsa.js");

const nSets = Number(process.argv[2] || 30);
const M = Object.keys(S.METHODS);
const fmt = r => M.map(m => String(Math.round(100 * r[m])).padStart(5)).join(" ");
console.log("".padEnd(44) + M.map(m => m.slice(0, 5).padStart(5)).join(" "));

function row(label, sim, prep, spec) {
  const t0 = Date.now();
  const r = S.recovery(sim, prep, spec, nSets, 123);
  console.log(label.padEnd(44) + fmt(r) + `   (${((Date.now() - t0) / 1000).toFixed(1)}s)`);
}

for (const scenario of [1, 2]) {
  const sim = S.simulate({ scenario, cellsPerGroup: 200, seed: 7 });
  const prep = S.prepare(sim);
  const ds = S.dataSummary(sim);
  console.log(`\n-- scenario ${scenario}: median % cells expressing signal genes = ${ds.medianPctCellsExpressingSignal.toFixed(1)}`);
  for (const noise of [0, 0.5, 0.8, 1]) row(`  S${scenario} size100 noise ${noise * 100}%`, sim, prep, { size: 100, noise });
  row(`  S${scenario} size10  noise 0%`, sim, prep, { size: 10, noise: 0 });
}
for (const scenario of [3, 4]) {
  const sim = S.simulate({ scenario, cellsPerGroup: 200, seed: 7 });
  const prep = S.prepare(sim);
  const ds = S.dataSummary(sim);
  console.log(`\n-- scenario ${scenario}: genes detected per cell control/treatment = ${ds.genesPerCell.map(x => x.toFixed(0)).join(" / ")}`);
  for (const size of [10, 50, 100, 200]) row(`  S${scenario} size${size} WITHOUT cond-specific`, sim, prep, { size, withCS: false, csFrac: 0 });
  row(`  S${scenario} size100 WITH 2% cond-specific`, sim, prep, { size: 100, withCS: true, csFrac: 0.02 });
  row(`  S${scenario} size100 WITH 10% cond-specific`, sim, prep, { size: 100, withCS: true, csFrac: 0.1 });
}
for (const n of [20, 50, 200]) {
  const sim = S.simulate({ scenario: 2, cellsPerGroup: n, seed: 7 });
  const prep = S.prepare(sim);
  row(`\n  S2 cells/group ${n} size100 noise 0%`, sim, prep, { size: 100, noise: 0 });
}
