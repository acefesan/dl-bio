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
  let sim = null, prep = null, last = null, selected = null;
  const pools = {};
  const PAPER = { paper: { dir: "naive_b", threshold: 0.75, what: "the real cells (~12 MB, once)" },
                  ssd: { dir: "ssd", threshold: 0.9, what: "the Splatter-simulated cells (once)" } };

  const tick = () => new Promise(r => setTimeout(r, 30));
  const pfmt = p => (p < 1e-4 ? p.toExponential(1) : p.toFixed(4));

  function syncBoxes() {
    const sc = Number($("scenario").value);
    const paper = $("regime").value in PAPER;
    $("noiseBox").hidden = !paper && sc >= 3;
    $("csBox").hidden = paper || sc <= 2;
    $("boostBox").hidden = !paper || sc <= 2;
    $("injectBox").hidden = paper;
  }

  // The real cells: all 9,528 naive B cells from GSE164381's 4 healthy donors,
  // raw counts (paper-data/03_export_browser.py). Fetched once, ~12 MB.
  // SSD: one 4,000-cell Splatter pool fitted to those cells (paper-data/05).
  async function loadPool(kind) {
    if (pools[kind]) return pools[kind];
    $("simstatus").textContent = `Downloading ${PAPER[kind].what}…`;
    await tick();
    const d = PAPER[kind].dir;
    const [meta, buf] = await Promise.all([
      fetch(`data/${d}/meta.json`).then(r => r.json()),
      fetch(`data/${d}/counts.bin`).then(r => r.arrayBuffer()),
    ]);
    pools[kind] = S.decodePool(buf, meta);
    return pools[kind];
  }

  function enable(on) {
    ["one", "rec", "sw"].forEach(id => { $(id).disabled = !on; });
  }

  async function simulate() {
    const regime = $("regime").value;
    let cells = Number($("cells").value);
    let clamped = false;
    closeCell();
    if (regime === "sparse" && cells > 200) {
      cells = 200;
      clamped = true;
      $("cells").value = "200";
    }
    $("simulate").disabled = true;
    enable(false);
    const t0 = performance.now();
    if (regime in PAPER) {
      let pool;
      try { pool = await loadPool(regime); } catch (e) {
        $("simstatus").textContent = "Could not load the real cells: " + e.message;
        $("simulate").disabled = false;
        return;
      }
      $("simstatus").textContent = "Drawing cells and building the scenario (authors' procedure)…";
      await tick();
      sim = S.simulateRWSD(pool, {
        scenario: Number($("scenario").value),
        cellsPerGroup: cells,
        seed: Number($("seed").value) || 7,
        boostNull: $("boostNull").value === "1",
        threshold: PAPER[regime].threshold,
      });
      sim.kind = regime;
    } else {
      $("simstatus").textContent = "Simulating counts…";
      await tick();
      sim = S.simulate(Object.assign({
        scenario: Number($("scenario").value),
        cellsPerGroup: cells,
        seed: Number($("seed").value) || 7,
        inject: $("inject").value,
      }, REGIMES[regime]));
    }
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
    if (sim.paper && sc <= 2) {
      const ssd = sim.kind === "ssd";
      items.push([String(ds.nSignal), sc === 1 ? `signal genes (detected in >${ssd ? 10 : 25}% of cells), ×1.2` : "signal genes (550 random), ×1.2"]);
      items.push([`${ds.medianPctCellsExpressingSignal.toFixed(1)}%`, "median signal-gene detection" + (sc === 1 ? ` (paper: ${ssd ? 32 : 44.5}%)` : "")]);
    } else if (sim.paper) {
      items.push([String(ds.nSignal), `${sc === 3 ? "dense" : "random"} "signal" genes, ${sim.boost === 1 ? "not boosted" : "×1.2"}`]);
      items.push([String(sim.nCS), "added treatment-only genes"]);
      items.push([`${ds.genesPerCell[0].toFixed(0)} → ${ds.genesPerCell[1].toFixed(0)}`, "genes detected per cell, control → treatment"]);
    } else if (sc <= 2) {
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
    if (sim.paper || sc <= 2) return { size, noise: Number($("noise").value) };
    const cs = Number($("cs").value);
    return { size, withCS: cs > 0, csFrac: cs };
  }

  // Paper-mode S3/S4 sets have no fixed condition-specific share, so "null"
  // is decided per set (see setIsNull); this is the setting-level answer.
  function isNull() {
    const sc = sim.opt.scenario;
    if (sim.paper) return sc <= 2 ? Number($("noise").value) === 1 : sim.boost === 1 && Number($("noise").value) === 0;
    return sc <= 2 ? Number($("noise").value) === 1 : Number($("cs").value) === 0;
  }

  function setIsNull(set) {
    if (!sim.paper || sim.opt.scenario <= 2) return isNull();
    const cs = new Set(sim.csGenes);
    const hasCS = set.some(g => cs.has(g));
    const hasBoost = sim.boost !== 1 && set.some(g => sim.isSignal[g]);
    return !hasCS && !hasBoost;
  }

  function strip(scores, group) {
    let lo = Infinity, hi = -Infinity;
    for (const v of scores) { if (v < lo) lo = v; if (v > hi) hi = v; }
    const span = hi - lo || 1;
    const W = 400, pad = 6;
    const x = v => pad + (W - 2 * pad) * (v - lo) / span;
    const ctl = [], trt = [];
    scores.forEach((v, i) => (group[i] ? trt : ctl).push([v, i]));
    const med = a => { const s = a.map(e => e[0]).sort((p, q) => p - q); return s[Math.floor(s.length / 2)]; };
    const pick = (a, n) => a.length <= n ? a : a.filter((_, i) => i % Math.ceil(a.length / n) === 0);
    const jit = i => ((i * 2654435761) % 1000) / 1000;
    // Every drawn dot carries its cell index; click one to open that cell.
    let h = `<svg viewBox="0 0 ${W} 46" preserveAspectRatio="none">`;
    pick(ctl, 160).forEach(([v, c]) => { h += `<circle class="dot-c cell" data-c="${c}" cx="${x(v).toFixed(1)}" cy="${(6 + 12 * jit(c)).toFixed(1)}" r="2.2"/>`; });
    pick(trt, 160).forEach(([v, c]) => { h += `<circle class="dot-t cell" data-c="${c}" cx="${x(v).toFixed(1)}" cy="${(27 + 12 * jit(c + 7)).toFixed(1)}" r="2.2"/>`; });
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
    last = { set, res };
    closeCell();
    $("strips").innerHTML = METHODS.map(m => {
      const r = res[m];
      const hit = r.test.p < 0.05;
      const arrow = r.test.z > 0 ? "↑" : r.test.z < 0 ? "↓" : "";
      const extra = m === "scPS" && r.meta ? `<small>${r.meta.maxComb} of ${r.meta.npcs} PCs</small>` : `<small>${DESC[m]}</small>`;
      return `<div class="mrow"><div class="name">${m}${extra}</div>${strip(r.scores, sim.group)}` +
        `<div class="p">${hit ? `<span class="pill ${setIsNull(set) ? "bad" : "good"}">p ${pfmt(r.test.p)} ${arrow}</span>` : `p ${pfmt(r.test.p)}`}</div></div>`;
    }).join("");
    const sc = sim.opt.scenario;
    let note = `A fresh random gene set of ${set.length} genes. `;
    if (sim.paper && sim.opt.scenario >= 3) {
      const nCS = set.filter(g => sim.csGenes.includes(g)).length;
      note += `It happens to contain ${nCS} added treatment-only gene${nCS === 1 ? "" : "s"} (${(100 * nCS / set.length).toFixed(0)}%): in the authors' design that share is left to chance. `;
    }
    note += setIsNull(set)
      ? "This set contains no gene that differs between the groups, so every highlighted p-value is a false positive."
      : sc <= 2 ? "This set contains signal genes, so highlighted p-values are correct detections." : "This set contains treatment-only genes, so a detection is correct.";
    note += " Arrows give the direction of the treatment shift; a downward shift is still a detection, because the test is two-sided. Click any dot to open that cell.";
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
    const mixed = sim.paper && sim.opt.scenario >= 3 && !nul && sim.boost === 1;
    $("recTitle").textContent = mixed ? "Detection rate (sets hold added treatment-only genes only by chance)"
      : nul ? "False-positive rate (these sets contain no real difference)" : "Sensitivity (these sets contain a real difference)";
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

  // ------------------------------------------------------------ one cell
  // Side panel for a clicked dot: where this cell sits under every method,
  // and the gene set's genes in this cell, boosted ("signal") genes on top,
  // as in the paper's Figure 1.
  const esc = t => String(t).replace(/[&<>"]/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));

  function closeCell() {
    selected = null;
    const d = $("cellPanel");
    if (d) d.hidden = true;
    document.querySelectorAll("circle.cell.sel").forEach(e => e.classList.remove("sel"));
  }

  function openCell(c) {
    if (!last) return;
    selected = c;
    document.querySelectorAll("circle.cell.sel").forEach(e => e.classList.remove("sel"));
    document.querySelectorAll(`circle.cell[data-c="${c}"]`).forEach(e => e.classList.add("sel"));
    const { set, res } = last;
    const G = sim.G, N = sim.N, trt = sim.group[c] === 1;
    const row = g => sim.data[c * G + g];
    const mean = (g, grp) => { let s = 0, k = 0; for (let i = 0; i < N; i++) if (sim.group[i] === grp) { s += sim.data[i * G + g]; k++; } return s / k; };
    const csSet = new Set(sim.csGenes);
    const role = g => csSet.has(g) ? "added" : sim.isSignal[g] ? "signal" : "non-signal";
    const boosted = g => trt && sim.isSignal[g] && (sim.paper ? sim.boost !== 1 : true);
    const order = { signal: 0, added: 1, "non-signal": 2 };
    const genes = set.slice().sort((a, b) => (order[role(a)] - order[role(b)]) || (row(b) - row(a)));
    let detected = 0;
    for (let g = 0; g < G; g++) if (row(g) > 0) detected++;
    const inSet = genes.filter(g => row(g) > 0).length;
    let vmax = 0;
    const rows = genes.map(g => { const r = [row(g), mean(g, 0), mean(g, 1)]; vmax = Math.max(vmax, ...r); return r; });
    const cellStyle = v => `background:color-mix(in srgb, var(--s1) ${Math.round(100 * Math.min(1, v / (vmax || 1)))}%, transparent)`;

    let h = `<button class="cp-close" type="button" aria-label="Close">×</button>`;
    h += `<div class="kicker">${trt ? "Treatment" : "Control"} cell #${c + 1}</div>`;
    if (sim.paper) {
      const ci = sim.cellInfo[c];
      h += `<h3>${esc(ci.donor)} · <code>${esc(ci.barcode)}</code></h3>`;
      if (ci.addedFrom !== null && sim.nCS) h += `<p class="cp-note">Its ${sim.nCS} added genes were copied from another real cell (pool #${ci.addedFrom + 1}), as the authors' script does.</p>`;
    } else h += `<h3>Simulated cell</h3>`;
    h += `<div class="summary cp-sum"><div><b>${detected.toLocaleString()}</b><span>genes detected in this cell (of ${G.toLocaleString()})</span></div>` +
      `<div><b>${inSet} / ${set.length}</b><span>gene-set genes detected here</span></div></div>`;

    h += `<h4>Its score under each method</h4><table class="cp-scores"><tr><th>Method</th><th class="num">score</th><th class="num">rank among ${N}</th></tr>`;
    for (const m of METHODS) {
      const s = res[m].scores, v = s[c];
      let above = 0;
      for (let i = 0; i < N; i++) if (s[i] > v) above++;
      h += `<tr><td>${m}</td><td class="num">${Number(v).toPrecision(3)}</td><td class="num">${above + 1}</td></tr>`;
    }
    h += `</table><p class="cp-note">Rank 1 = highest score of all ${N} cells. The cell is highlighted in every strip.</p>`;

    h += `<h4>The gene set in this cell</h4>`;
    h += `<p class="cp-note">Log-normalised expression, log1p(count ÷ cell total × 10⁴). Darker = higher (scale shared across the three columns). ` +
      `<span class="role signal">signal</span> genes were multiplied by ${sim.paper ? sim.boost : sim.opt.foldChange} in treatment cells` +
      (sim.paper && sim.opt.scenario >= 3 && sim.boost === 1 ? " — not in this run, so here they are just dense genes" : "") +
      `; <span class="role added">added</span> genes exist only in treatment cells. Signal genes first, as in the paper's Figure 1.</p>`;
    h += `<div class="table-wrap"><table class="cp-genes"><tr><th>Gene</th><th></th><th class="num">this cell</th><th class="num">control mean</th><th class="num">treatment mean</th></tr>`;
    genes.forEach((g, k) => {
      const [v, mc, mt] = rows[k], r = role(g);
      const pre = boosted(g) && v > 0 ? `<small title="value before the ×${sim.paper ? sim.boost : sim.opt.foldChange}">${(v / (sim.paper ? sim.boost : sim.opt.foldChange)).toFixed(2)}→</small>` : "";
      h += `<tr><td>${esc(sim.geneNames[g])}</td><td><span class="role ${r === "non-signal" ? "non" : r}">${r}</span></td>` +
        `<td class="num hm" style="${cellStyle(v)}">${pre}${v.toFixed(2)}</td>` +
        `<td class="num hm" style="${cellStyle(mc)}">${mc.toFixed(2)}</td><td class="num hm" style="${cellStyle(mt)}">${mt.toFixed(2)}</td></tr>`;
    });
    h += `</table></div>`;
    const d = $("cellPanel");
    d.innerHTML = h;
    d.hidden = false;
    d.querySelector(".cp-close").addEventListener("click", closeCell);
  }

  const PRESETS = {
    "s1-null": { scenario: 1, regime: "paper", cells: 200, size: 100, noise: "1" },
    "s2-signal": { scenario: 2, regime: "paper", cells: 20, size: 100, noise: "0" },
    "s2-small": { scenario: 2, regime: "dense", cells: 20, size: 25, noise: "0" },
    "s3-null": { scenario: 3, regime: "paper", cells: 200, size: 100, noise: "0" },
    "s4-null": { scenario: 4, regime: "paper", cells: 200, size: 100, noise: "0" },
  };

  async function applyHash() {
    const p = PRESETS[location.hash.slice(1)];
    if (!p) return;
    $("scenario").value = String(p.scenario);
    $("regime").value = p.regime;
    $("cells").value = String(p.cells);
    $("inject").value = "log";
    $("boostNull").value = "0";
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
    $("regime").addEventListener("change", syncBoxes);
    $("strips").addEventListener("click", e => {
      const t = e.target.closest("circle.cell");
      if (t) openCell(Number(t.dataset.c));
    });
    document.addEventListener("keydown", e => { if (e.key === "Escape") closeCell(); });
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
