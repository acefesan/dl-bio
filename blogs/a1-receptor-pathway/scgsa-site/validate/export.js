// Export one simulated dataset + gene set, and the JS engine's scores, so
// validate.R can score the same inputs with the real R packages.
// usage: node export.js <outdir> <scenario>
const fs = require("fs");
const path = require("path");
const S = require("../static/scgsa.js");

const out = process.argv[2] || "/tmp/scgsa-validate";
const scenario = Number(process.argv[3] || 1);
fs.mkdirSync(out, { recursive: true });

const inject = process.argv[4] || "log";
const sim = S.simulate({ genes: 800, cellsPerGroup: 60, scenario, seed: 3, inject });
const prep = S.prepare(sim, { ucellMaxRankFrac: 60 / 800 });
const rng = S.makeRng(5);
const spec = scenario <= 2 ? { size: 30, noise: 0.2 } : { size: 30, withCS: true, csFrac: 0.1 };
const set = S.makeSet(sim, spec, rng);

const { G, N, counts, data } = sim;
const genes = Array.from({ length: G }, (_, g) => `g${g}`);
const cells = Array.from({ length: N }, (_, c) => `c${c}`);

let csv = "gene," + cells.join(",") + "\n";
for (let g = 0; g < G; g++) {
  const row = [];
  for (let c = 0; c < N; c++) row.push(counts[c * G + g]);
  csv += genes[g] + "," + row.join(",") + "\n";
}
fs.writeFileSync(path.join(out, "counts.csv"), csv);

let dcsv = "gene," + cells.join(",") + "\n";
for (let g = 0; g < G; g++) {
  const row = [];
  for (let c = 0; c < N; c++) row.push(data[c * G + g]);
  dcsv += genes[g] + "," + row.join(",") + "\n";
}
fs.writeFileSync(path.join(out, "data_js.csv"), dcsv);

fs.writeFileSync(path.join(out, "set.txt"), set.map(g => genes[g]).join("\n") + "\n");
fs.writeFileSync(path.join(out, "group.txt"), Array.from(sim.group).join("\n") + "\n");

const res = S.scoreAll(sim, prep, set, { AddModuleScore: { ctrl: 20 } });
let scsv = "cell," + Object.keys(res).join(",") + "\n";
for (let c = 0; c < N; c++) scsv += cells[c] + "," + Object.keys(res).map(k => res[k].scores[c]).join(",") + "\n";
fs.writeFileSync(path.join(out, "scores_js.csv"), scsv);

fs.writeFileSync(path.join(out, "params.json"), JSON.stringify({
  aucMaxRank: prep.aucMaxRank, ucellMaxRank: prep.ucellMaxRank, ctrl: 20, inject,
  scps: res.scPS.meta, setSize: set.length,
}));
console.log("exported", out, "G", G, "N", N, "set", set.length,
  "aucMaxRank", prep.aucMaxRank, "ucellMaxRank", prep.ucellMaxRank);
