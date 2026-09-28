(function () {
  const S = window.SCGSA;
  const $ = id => document.getElementById(id);
  const METHODS = Object.keys(S.METHODS);
  const DESC = {
    AddModuleScore: "vs control genes",
    AUCell: "top-5% ranks",
    UCell: "Mann–Whitney",
    ssGSEA: "random walk",
    JASMINE: "rank + detection",
    SCSE: "share of total",
    scPS: "PCA in the set",
  };
  const REGIMES = {
    dense: { genes: 3000, meanLog: -1.2, sdLog: 1.6 },
    sparse: { genes: 6000, meanLog: -3.2, sdLog: 2.2 },
  };
  let sim = null, prep = null;

  const tick = () => new Promise(r => setTimeout(r, 30));
  const pfmt = p => (p < 1e-4 ? p.toExponential(1) : p.toFixed(4));

  function syncBoxes() {
    const sc = Number($("scenario").value);
    $("noiseBox").hidden = sc >= 3;
    $("csBox").hidden = sc <= 2;
  }

  function enable(on) {
    ["one", "rec", "sw"].forEach(id => { $(id).disabled = !on; });
  }

  async function simulate() {
    const regime = $("regime").value;
    let cells = Number($("cells").value);
    let clamped = false;
    if (regime === "sparse" && cells > 200) {
      cells = 200;
      clamped = true;
      $("cells").value = "200";
    }
    $("simulate").disabled = true;
    enable(false);
    $("simstatus").textContent = "Simulating counts…";
    await tick();
    const t0 = performance.now();
    sim = S.simulate(Object.assign({
      scenario: Number($("scenario").value),
      cellsPerGroup: cells,
      seed: Number($("seed").value) || 7,
      inject: $("inject").value,
    }, REGIMES[regime]));
    $("simstatus").textContent = "Ranking every cell's genes…";
    await tick();
    prep = S.prepare(sim);
    const ds = S.dataSummary(sim);
    const secs = ((performance.now() - t0) / 1000).toFixed(1);
    $("simstatus").textContent = `Done in ${secs}s.` + (clamped ? " The sparse regime is capped at 200 cells per group to keep memory reasonable on phones." : "");
    const sc = sim.opt.scenario;
    const items = [
      [sim.G.toLocaleString(), "genes"],
      [sim.N.toLocaleString(), "cells"],
      [`${(100 * ds.sparsity).toFixed(0)}%`, "of entries are zero"],
      [`${ds.medianPctCellsExpressingAll.toFixed(1)}%`, "median gene detection"],
    ];
    if (sc <= 2) {
      items.push([String(ds.nSignal), "signal genes (+20%)"]);
      items.push([`${ds.medianPctCellsExpressingSignal.toFixed(1)}%`, "median signal-gene detection"]);
    } else {
      items.push([String(sim.nCS), "treatment-only genes"]);
      items.push([`${ds.genesPerCell[0].toFixed(0)} → ${ds.genesPerCell[1].toFixed(0)}`, "genes detected per cell, control → treatment"]);
    }
    $("summary").innerHTML = items.map(([b, s]) => `<div><b>${b}</b><span>${s}</span></div>`).join("");
    $("simulate").disabled = false;
    enable(true);
    ["onePanel", "recPanel", "sweep"].forEach(id => { $(id).hidden = true; });
  }

  function spec(size) {
    const sc = sim.opt.scenario;
    if (sc <= 2) return { size, noise: Number($("noise").value) };
    const cs = Number($("cs").value);
    return { size, withCS: cs > 0, csFrac: cs };
  }

  function isNull() {
    const sc = sim.opt.scenario;
    return sc <= 2 ? Number($("noise").value) === 1 : Number($("cs").value) === 0;
  }

  function strip(scores, group) {
    let lo = Infinity, hi = -Infinity;
    for (const v of scores) { if (v < lo) lo = v; if (v > hi) hi = v; }
    const span = hi - lo || 1;
    const W = 400, pad = 6;
    const x = v => pad + (W - 2 * pad) * (v - lo) / span;
    const ctl = [], trt = [];
    scores.forEach((v, i) => (group[i] ? trt : ctl).push(v));
    const med = a => { const s = a.slice().sort((p, q) => p - q); return s[Math.floor(s.length / 2)]; };
    const pick = (a, n) => a.length <= n ? a : a.filter((_, i) => i % Math.ceil(a.length / n) === 0);
    const jit = i => ((i * 2654435761) % 1000) / 1000;
    let h = `<svg viewBox="0 0 ${W} 46" preserveAspectRatio="none">`;
    pick(ctl, 160).forEach((v, i) => { h += `<circle class="dot-c" cx="${x(v).toFixed(1)}" cy="${(6 + 12 * jit(i)).toFixed(1)}" r="2.2"/>`; });
    pick(trt, 160).forEach((v, i) => { h += `<circle class="dot-t" cx="${x(v).toFixed(1)}" cy="${(27 + 12 * jit(i + 7)).toFixed(1)}" r="2.2"/>`; });
    const mc = x(med(ctl)), mt = x(med(trt));
    h += `<line class="med" x1="${mc}" x2="${mc}" y1="3" y2="21"/><line class="med" x1="${mt}" x2="${mt}" y1="24" y2="43"/></svg>`;
    return h;
  }

  async function scoreOne() {
    enable(false);
    $("setstatus").textContent = "Scoring…";
    await tick();
    const rng = S.makeRng(Math.floor(Math.random() * 1e9));
    const set = S.makeSet(sim, spec(Number($("size").value)), rng);
    const res = S.scoreAll(sim, prep, set);
    $("strips").innerHTML = METHODS.map(m => {
      const r = res[m];
      const hit = r.test.p < 0.05;
      const arrow = r.test.z > 0 ? "↑" : r.test.z < 0 ? "↓" : "";
      const extra = m === "scPS" && r.meta ? `<small>${r.meta.maxComb} of ${r.meta.npcs} PCs</small>` : `<small>${DESC[m]}</small>`;
      return `<div class="mrow"><div class="name">${m}${extra}</div>${strip(r.scores, sim.group)}` +
        `<div class="p">${hit ? `<span class="pill ${isNull() ? "bad" : "good"}">p ${pfmt(r.test.p)} ${arrow}</span>` : `p ${pfmt(r.test.p)}`}</div></div>`;
    }).join("");
    const sc = sim.opt.scenario;
    let note = `A fresh random gene set of ${set.length} genes. `;
    note += isNull()
      ? "This set contains no gene that differs between the groups, so every highlighted p-value is a false positive."
      : sc <= 2 ? "This set contains signal genes, so highlighted p-values are correct detections." : "This set contains treatment-only genes, so a detection is correct.";
    note += " Arrows give the direction of the treatment shift; a downward shift is still a detection, because the test is two-sided.";
    $("onenote").textContent = note;
    $("onePanel").hidden = false;
    $("setstatus").textContent = "";
    enable(true);
  }

  async function recovery() {
    enable(false);
    $("setstatus").textContent = "Scoring 30 gene sets × 7 methods…";
    await tick();
    const size = Number($("size").value);
    const r = S.recovery(sim, prep, spec(size), 30, Math.floor(Math.random() * 1e9));
    const nul = isNull();
    $("recTitle").textContent = nul ? "False-positive rate (these sets contain no real difference)" : "Sensitivity (these sets contain a real difference)";
    $("recTable").innerHTML = `<tr><th>Method</th><th class="num">sets identified</th><th></th></tr>` + METHODS.map(m => {
      const v = Math.round(100 * r[m]);
      const cls = nul ? (v <= 10 ? "good" : v <= 40 ? "mid" : "bad") : (v >= 80 ? "good" : v >= 40 ? "mid" : "bad");
      return `<tr><td>${m}</td><td class="num"><span class="pill ${cls}">${v}%</span></td>` +
        `<td style="width:50%"><span style="display:inline-block;height:9px;border-radius:3px;background:var(--${cls === "good" ? "good" : cls === "mid" ? "mid" : "bad"});width:${v}%"></span></td></tr>`;
    }).join("");
    $("recnote").textContent = `30 random gene sets of size ${size}; each scored by every method and tested with a two-sided Wilcoxon test; p-values adjusted with Benjamini–Hochberg per method; identified = adjusted p < 0.05. Rerun for a different random draw of sets.`;
    $("recPanel").hidden = false;
    $("setstatus").textContent = "";
    enable(true);
  }

  async function sweep() {
    enable(false);
    const sizes = [10, 25, 50, 100];
    const out = {};
    METHODS.forEach(m => { out[m] = []; });
    for (const size of sizes) {
      $("setstatus").textContent = `Sweeping: size ${size}…`;
      await tick();
      const r = S.recovery(sim, prep, spec(size), 20, 1000 + size);
      METHODS.forEach(m => out[m].push(r[m]));
    }
    const nul = isNull();
    $("swTitle").textContent = nul ? "False-positive rate vs. gene set size" : "Sensitivity vs. gene set size";
    const W = 520, H = 230, L = 40, B = 34, T = 12, R = 120;
    const X = i => L + i * (W - L - R) / (sizes.length - 1);
    const Y = v => T + (1 - v) * (H - T - B);
    const colors = ["#3b7bbf", "#d0762b", "#3d9a8b", "#8a5cc2", "#b23b2a", "#b07d12", "#2e7d4f"];
    let h = `<svg viewBox="0 0 ${W} ${H}">`;
    h += `<line class="ax" x1="${L}" y1="${Y(0)}" x2="${W - R}" y2="${Y(0)}"/><line class="ax" x1="${L}" y1="${Y(0)}" x2="${L}" y2="${Y(1)}"/>`;
    [0, 0.5, 1].forEach(v => { h += `<text x="${L - 30}" y="${Y(v) + 4}">${v * 100}%</text>`; });
    sizes.forEach((s, i) => { h += `<text x="${X(i) - 6}" y="${H - 14}">${s}</text>`; });
    h += `<text x="${(L + W - R) / 2 - 40}" y="${H - 1}">gene set size</text>`;
    METHODS.forEach((m, k) => {
      const pts = out[m].map((v, i) => `${X(i).toFixed(1)},${Y(v).toFixed(1)}`);
      h += `<polyline fill="none" stroke="${colors[k]}" stroke-width="2" points="${pts.join(" ")}"/>`;
      out[m].forEach((v, i) => { h += `<circle cx="${X(i)}" cy="${Y(v)}" r="3" fill="${colors[k]}"/>`; });
      h += `<text x="${W - R + 8}" y="${T + 12 + k * 17}" style="fill:${colors[k]};font-weight:600">${m}</text>`;
    });
    h += "</svg>";
    $("swChart").innerHTML = h;
    $("swnote").textContent = nul
      ? "Every point is a rate of detecting differences that are not there. Rising lines are the 'bigger is shadier' effect: averaging more genes shrinks noise but not a shared systematic shift."
      : "Every point is a rate of detecting a real difference. The same averaging that amplifies artifacts also amplifies real signal.";
    $("sweep").hidden = false;
    $("setstatus").textContent = "";
    enable(true);
  }

  const PRESETS = {
    "s1-null": { scenario: 1, regime: "dense", cells: 200, size: 100, noise: "1" },
    "s2-signal": { scenario: 2, regime: "sparse", cells: 20, size: 100, noise: "0" },
    "s2-small": { scenario: 2, regime: "dense", cells: 20, size: 25, noise: "0" },
    "s3-null": { scenario: 3, regime: "dense", cells: 200, size: 100, cs: "0" },
    "s4-null": { scenario: 4, regime: "dense", cells: 200, size: 100, cs: "0" },
  };

  async function applyHash() {
    const p = PRESETS[location.hash.slice(1)];
    if (!p) return;
    $("scenario").value = String(p.scenario);
    $("regime").value = p.regime;
    $("cells").value = String(p.cells);
    $("inject").value = "log";
    $("size").value = String(p.size);
    if (p.noise !== undefined) $("noise").value = p.noise;
    if (p.cs !== undefined) $("cs").value = p.cs;
    syncBoxes();
    await simulate();
    await recovery();
    $("recPanel").scrollIntoView({ behavior: "smooth", block: "center" });
  }

  function init() {
    $("scenario").addEventListener("change", syncBoxes);
    $("simulate").addEventListener("click", simulate);
    $("one").addEventListener("click", scoreOne);
    $("rec").addEventListener("click", recovery);
    $("sw").addEventListener("click", sweep);
    window.addEventListener("hashchange", applyHash);
    syncBoxes();
    if (PRESETS[location.hash.slice(1)]) applyHash();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
