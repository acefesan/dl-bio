// Compare simulator regimes against the paper's Table 1 (size 100, 200 cells,
// sensitivity = noise 0%, false-positive rate = noise 100%).
const S = require("../static/scgsa.js");
const M = Object.keys(S.METHODS);
const nSets = Number(process.argv[2] || 20);

const regimes = [
  { label: "dense toy  G3000 ln(-1.2,1.6)", genes: 3000, meanLog: -1.2, sdLog: 1.6 },
  { label: "sparse     G6000 ln(-3.2,2.2)", genes: 6000, meanLog: -3.2, sdLog: 2.2 },
  { label: "sparser    G6000 ln(-3.6,2.3)", genes: 6000, meanLog: -3.6, sdLog: 2.3 },
];

console.log("".padEnd(52) + M.map(m => m.slice(0, 5).padStart(6)).join(""));
for (const R of regimes) {
  for (const sc of [1, 2]) {
    const t0 = Date.now();
    const sim = S.simulate({ scenario: sc, seed: 7, genes: R.genes, meanLog: R.meanLog, sdLog: R.sdLog });
    const prep = S.prepare(sim);
    const ds = S.dataSummary(sim);
    const tag = `${R.label} S${sc}`;
    for (const noise of [0, 1]) {
      const r = S.recovery(sim, prep, { size: 100, noise }, nSets, 123);
      console.log(`${tag} ${noise ? "FPR " : "sens"}`.padEnd(52) + M.map(m => String(Math.round(100 * r[m])).padStart(6)).join(""));
    }
    console.log(`   nSignal=${ds.nSignal} dense=${ds.nDense} medianDetect(all)=${ds.medianPctCellsExpressingAll.toFixed(1)}% `
      + `medianDetect(signal)=${ds.medianPctCellsExpressingSignal.toFixed(1)}% sparsity=${(100 * ds.sparsity).toFixed(1)}% `
      + `(${((Date.now() - t0) / 1000).toFixed(1)}s)`);
  }
}
console.log("\npaper Table 1 (RWSD) S1 sens: AMS 100 AUC 100 UC 100 ssG 17 JAS 100 SCSE 100 scPS 99.5");
console.log("paper Table 1 (RWSD) S1 FPR : AMS 4.9 AUC 100 UC 11.8 ssG 0.7 JAS 99.7 SCSE 40 scPS 4.1");
console.log("paper Table 1 (RWSD) S2 sens: AMS 98 AUC 99 UC 96 ssG 1.5 JAS 100 SCSE 99 scPS 99");
console.log("paper Table 1 (RWSD) S2 FPR : AMS 2.9 AUC 17 UC 1.8 ssG 0 JAS 25 SCSE 2.7 scPS 3.5");
