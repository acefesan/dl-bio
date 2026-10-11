// The playground as a pipeline graph: one node per processing step, from the
// source data (root) to the test results. Each node holds its own options,
// shows what it produced, and links to the code that does it — ours
// (static/scgsa.js) and the authors' R (Thakar-Lab/scPS, pinned commit).
// Changing an option reruns only that node and the ones below it.
(function () {
  const S = window.SCGSA;
  const $ = id => document.getElementById(id);
  const METHODS = Object.keys(S.METHODS);
  const esc = t => String(t).replace(/[&<>"]/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));
  const tick = () => new Promise(r => setTimeout(r, 20));
  const pfmt = p => (p < 1e-4 ? p.toExponential(1) : p.toFixed(4));
  const pct = x => `${(100 * x).toFixed(1)}%`;

  // ------------------------------------------------------------ code references
  const OURS = "https://github.com/acefesan/dl-bio/blob/main/blogs/a1-receptor-pathway/";
  const AUTH_SHA = "e690b450d7b7be204df38c42febc641cf1b25d88";
  const AUTH_BLOB = `https://github.com/Thakar-Lab/scPS/blob/${AUTH_SHA}/`;
  const AUTH_RAW = `https://raw.githubusercontent.com/Thakar-Lab/scPS/${AUTH_SHA}/`;
  const RWSD = "Comparison analysis/Rscript RWSD and GS simulation.rmd";
  const SSD = "Comparison analysis/Rscript SSD and GS simulation.rmd";
  const WILCOX = "Comparison analysis/Rscript Wilcox test.rmd";
  const IRGSEA = "Comparison analysis/Rscript irGSEA.R";
  // ours: a string to find in static/scgsa.js; the snippet runs to the end of
  // that block (or `lines` lines). authors: [file, firstLine, lastLine].
  const CODE = {
    data:     { ours: { find: "function decodePool", label: "decodePool()" },
                pipe: { href: "paper-data/README.md", label: "paper-data/ (R + Python that built the file)" },
                auth: [[RWSD, 29, 29, "load(\"Naive Normal B cell.RData\")"]] },
    toy:      { ours: { find: "function simulate(opt)", label: "simulate()" } },
    draw:     { ours: { find: "const cells = o.cells", lines: 3, label: "simulateRWSD(): draw" },
                auth: [[RWSD, 48, 48, "sim.mc: draw 2n cells"], [RWSD, 76, 77, "split into conditions"]] },
    norm:     { ours: { find: "function poolRow", label: "poolRow()" },
                pipe: { href: "paper-data/02_export_naive_b.R", label: "NormalizeData (paper-data/02)" } },
    filter:   { ours: { find: "const nz = new Uint32Array(PG)", lines: 6, label: "simulateRWSD(): keep detected" },
                auth: [[RWSD, 85, 89, "rows_zeros, Genes"]] },
    signal:   { ours: { find: "let signalPool;", lines: 4, label: "simulateRWSD(): signal" },
                auth: [[RWSD, 92, 95, "S1/S3: rows_zeros < Threshold"], [RWSD, 224, 224, "S2/S4: 550 random"], [SSD, 26, 26, "SSD Threshold = 0.9"]] },
    boost:    { ours: { find: "if (group[i] === 1 && boost !== 1", lines: 2, label: "simulateRWSD(): ×boost" },
                auth: [[RWSD, 44, 46, "Threshold, Signal_Level"], [RWSD, 101, 101, "dat[Signal, cond 2] * Signal_Level"]] },
    added:    { ours: { find: "let added = [], cells2 = null;", lines: 11, label: "simulateRWSD(): added genes" },
                auth: [[RWSD, 350, 359, "S3: Add_gene … rbind"], [RWSD, 464, 464, "S4: Add_gene (any detected)"]] },
    geneset:  { ours: { find: "function makeSet", label: "makeSet()" },
                auth: [[RWSD, 130, 160, "gene sets: size × noise"]] },
    prepare:  { ours: { find: "function prepare(sim, opt)", label: "prepare()" },
                auth: [[IRGSEA, 24, 32, "irGSEA.score(...)"]] },
    AddModuleScore: { ours: { find: "function addModuleScore", label: "addModuleScore()" }, auth: [["Comparison analysis/Rscript AddmoduleScore.R", 28, 28, "Seurat AddModuleScore"]] },
    AUCell:   { ours: { find: "function aucell", label: "aucell()" }, auth: [[IRGSEA, 24, 32, "irGSEA: AUCell"]] },
    UCell:    { ours: { find: "function ucell", label: "ucell()" }, auth: [[IRGSEA, 24, 32, "irGSEA: UCell (maxRank 1500)"]] },
    ssGSEA:   { ours: { find: "function ssgsea", label: "ssgsea()" }, auth: [[IRGSEA, 24, 32, "irGSEA: ssgsea"]] },
    JASMINE:  { ours: { find: "function jasmine", label: "jasmine()" }, auth: [[IRGSEA, 24, 32, "irGSEA: JASMINE (oddsratio)"]] },
    SCSE:     { ours: { find: "function scse", label: "scse()" }, auth: [[IRGSEA, 24, 32, "irGSEA: scSE"]] },
    scPS:     { ours: { find: "function scps", label: "scps()" }, auth: [["Comparison analysis/Rscript scPS.R", 34, 70, "RunPCA(npcs = 10) …"]] },
    test:     { ours: { find: "function wilcox", label: "wilcox()" }, auth: [[WILCOX, 124, 124, "row_wilcoxon_twosample(two.sided)"]] },
    recovery: { ours: { find: "function recovery", label: "recovery() + bh()" }, auth: [[WILCOX, 125, 126, "p.adjust(fdr), FDR < 0.05"]] },
  };

  let srcText = null;
  async function ourSnippet(ref) {
    if (!srcText) srcText = await fetch("static/scgsa.js").then(r => r.text());
    const lines = srcText.split("\n");
    const i = lines.findIndex(l => l.includes(ref.find));
    if (i < 0) return { start: 0, end: 0, text: "(not found)" };
    let end = i;
    if (ref.lines) end = i + ref.lines - 1;
    else {
      let depth = 0, seen = false;
      for (let j = i; j < lines.length; j++) {
        for (const ch of lines[j]) { if (ch === "{") { depth++; seen = true; } else if (ch === "}") depth--; }
        if (seen && depth <= 0) { end = j; break; }
      }
    }
    return { start: i + 1, end: end + 1, text: lines.slice(i, end + 1).join("\n") };
  }
  const authCache = {};
  async function authSnippet([file, a, b]) {
    if (!authCache[file]) authCache[file] = fetch(AUTH_RAW + file.split("/").map(encodeURIComponent).join("/")).then(r => r.ok ? r.text() : Promise.reject(r.status));
    const t = await authCache[file];
    return t.split("\n").slice(a - 1, b).join("\n");
  }

  async function showCode(key) {
    const c = CODE[key];
    const panel = $("sidePanel");
    let h = `<button class="cp-close" type="button" aria-label="Close">×</button><div class="kicker">Code</div><h3>${esc(TITLES[key] || key)}</h3>`;
    if (c.ours) {
      const s = await ourSnippet(c.ours);
      h += `<h4>Ours · <a href="${OURS}scgsa-site/static/scgsa.js#L${s.start}-L${s.end}" target="_blank" rel="noopener">static/scgsa.js L${s.start}–${s.end}</a></h4>` +
        `<pre class="code"><code>${esc(s.text)}</code></pre>`;
    }
    if (c.pipe) h += `<p class="cp-note">Built offline by <a href="${OURS}${c.pipe.href}" target="_blank" rel="noopener">${esc(c.pipe.label)}</a>.</p>`;
    for (const ref of c.auth || []) {
      const [file, a, b, what] = ref;
      h += `<h4>Authors · <a href="${AUTH_BLOB}${file.split("/").map(encodeURIComponent).join("/")}#L${a}-L${b}" target="_blank" rel="noopener">${esc(file.split("/").pop())} L${a}${b > a ? "–" + b : ""}</a></h4><p class="cp-note">${esc(what)}</p>`;
      try { h += `<pre class="code"><code>${esc(await authSnippet(ref))}</code></pre>`; }
      catch (e) { h += `<p class="cp-note">(Couldn't fetch the snippet from GitHub; follow the link.)</p>`; }
    }
    panel.innerHTML = h;
    panel.hidden = false;
    panel.querySelector(".cp-close").addEventListener("click", closePanel);
  }

  function closePanel() {
    $("sidePanel").hidden = true;
    document.querySelectorAll("circle.cell.sel").forEach(e => e.classList.remove("sel"));
  }

  // ------------------------------------------------------------ state
  const SOURCES = {
    rwsd: { label: "RWSD · the paper's real cells", dir: "naive_b", threshold: 0.75, paper: true },
    ssd: { label: "SSD · Splatter fitted to those cells", dir: "ssd", threshold: 0.9, paper: true },
    dense: { label: "Toy · dense, 3,000 made-up genes", toy: { genes: 3000, meanLog: -1.2, sdLog: 1.6 } },
    sparse: { label: "Toy · sparse, 6,000 made-up genes", toy: { genes: 6000, meanLog: -3.2, sdLog: 2.2 } },
  };
  const cfg = { source: "rwsd", cells: 200, seed: 7, scenario: 1, fold: 1.2, boostNull: 0, nAdded: 250,
    inject: "log", size: 100, noise: 0, cs: 0.02, setSeed: 1, nSets: 30 };
  const SIM_KEYS = new Set(["source", "cells", "seed", "scenario", "fold", "boostNull", "nAdded", "inject"]);
  const pools = {};
  let sim = null, prep = null, set = null, res = null, poolNow = null;

  const TITLES = {
    data: "Source data", toy: "Simulate counts (toy model)", draw: "Draw cells and split into two groups",
    norm: "Log-normalise each cell", filter: "Keep genes detected in the drawn cells", signal: "Choose the signal genes",
    boost: "Boost signal genes in treatment cells", added: "Add treatment-only genes (S3/S4)",
    geneset: "Build one gene set", prepare: "Rank genes within every cell", score: "Score every cell — 7 methods",
    test: "Test: control vs. treatment", recovery: "Repeat over many gene sets", sweep: "Sweep gene set size",
  };

  const opt = (key, values, fmt) => `<label>${esc(fmt.label)}<select data-k="${key}">` +
    values.map(([v, t]) => `<option value="${v}"${String(cfg[key]) === String(v) ? " selected" : ""}>${esc(t)}</option>`).join("") + `</select></label>`;
  const codeBtns = key => {
    const c = CODE[key];
    if (!c) return "";
    let h = `<div class="pn-code">`;
    if (c.ours) h += `<button type="button" data-code="${key}">ours: <code>${esc(c.ours.label)}</code></button>`;
    if (c.auth) h += `<button type="button" data-code="${key}">authors' R: ${c.auth.map(a => `L${a[1]}`).join(", ")}</button>`;
    return h + `</div>`;
  };

  // ------------------------------------------------------------ graph layout
  function nodesFor() {
    const paper = !!SOURCES[cfg.source].paper;
    const list = ["data"];
    if (paper) list.push("draw", "norm", "filter", "signal", "boost", "added");
    else list.push("toy");
    list.push("geneset", "prepare", "score", "test", "recovery", "sweep");
    return list;
  }

  function controls(key) {
    const paper = !!SOURCES[cfg.source].paper;
    switch (key) {
      case "data": return opt("source", Object.entries(SOURCES).map(([k, s]) => [k, s.label]), { label: "Data" });
      case "draw": return opt("cells", [[20, "20"], [50, "50"], [200, "200"], [500, "500"]], { label: "Cells per group" }) +
        `<label>Seed<input type="number" data-k="seed" value="${cfg.seed}" min="1"></label>`;
      case "toy": return opt("scenario", [[1, "S1 · ×1.2 on dense genes"], [2, "S2 · ×1.2 on random genes"], [3, "S3 · extra dense genes in treatment"], [4, "S4 · extra sparse genes in treatment"]], { label: "Scenario" }) +
        opt("cells", [[20, "20"], [50, "50"], [200, "200"], [500, "500"]], { label: "Cells per group" }) +
        opt("inject", [["log", "after normalisation (as the paper)"], ["counts", "before normalisation (leaks)"]], { label: "Apply change" }) +
        `<label>Seed<input type="number" data-k="seed" value="${cfg.seed}" min="1"></label>`;
      case "signal": return opt("scenario", [[1, "S1 · dense genes (signal)"], [2, "S2 · 550 random genes (signal)"], [3, "S3 · dense genes, no signal"], [4, "S4 · 550 random genes, no signal"]], { label: "Scenario" });
      case "boost": return opt("fold", [[1.2, "×1.2 (paper)"], [1.1, "×1.1"], [1.5, "×1.5"], [2, "×2"]], { label: "Multiply log values by" }) +
        (cfg.scenario >= 3 ? opt("boostNull", [[0, "not boosted (authors' null runs)"], [1, "boosted (script as literally written)"]], { label: "S3/S4 dense genes" }) : "");
      case "added": return cfg.scenario >= 3 ? opt("nAdded", [[250, "250 (paper)"], [50, "50"], [100, "100"], [500, "500"], [1000, "1000"]], { label: "Genes added" }) : "";
      case "geneset": return opt("size", [[10, "10"], [25, "25"], [50, "50"], [100, "100"], [200, "200"]], { label: "Gene set size" }) +
        ((paper || cfg.scenario <= 2)
          ? opt("noise", [[0, "0% (all signal)"], [0.2, "20%"], [0.5, "50%"], [0.8, "80%"], [1, "100% (no signal)"]], { label: "Noise: non-signal share" })
          : opt("cs", [[0, "0% (null)"], [0.02, "2%"], [0.05, "5%"], [0.1, "10%"], [0.2, "20%"]], { label: "Treatment-only share" })) +
        `<button type="button" class="btn ghost" data-act="reroll">New random set</button>`;
      case "recovery": return opt("nSets", [[30, "30"], [100, "100"]], { label: "Gene sets" }) + `<button type="button" class="btn" data-act="recovery">Run</button>`;
      case "sweep": return `<button type="button" class="btn" data-act="sweep">Run sweep (sizes 10–100)</button>`;
      default: return "";
    }
  }

  const BLURB = {
    data: "The root. Everything below is computed from this, in your browser.",
    draw: "2n random cells; first half control, second half treatment. Same population, so any difference later is put there by the steps below.",
    norm: "Per cell: log1p(count ÷ cell total × 10,000) — Seurat's LogNormalize.",
    filter: "Drop genes with zero counts in every drawn cell.",
    signal: "Which genes the scenario changes (S1/S2) or treats as the gene-set pool (S3/S4).",
    boost: "The signal: their log-normalised values × the factor, in treatment cells only. Nothing else moves.",
    added: "Copy genes from a second, independent draw of cells into treatment cells only (zero in control). They join the non-signal pool.",
    geneset: "size × (1 − noise) genes from the signal pool, the rest from everything else.",
    prepare: "Per-cell work that doesn't depend on the gene set: rankings for AUCell, UCell, ssGSEA, JASMINE; bins for AddModuleScore.",
    score: "Each method turns the gene set into one number per cell. Click any dot to open that cell.",
    test: "Two-sided Wilcoxon rank-sum on the scores of control vs. treatment cells, per method.",
    recovery: "Many random gene sets with these settings; Benjamini–Hochberg per method; identified = adjusted p < 0.05.",
    sweep: "Recovery at sizes 10, 25, 50, 100 (20 sets each).",
    toy: "Negative-binomial counts for made-up genes; the scenario's change is applied inside. For mechanisms, not magnitudes.",
  };

  function render() {
    const g = $("graph");
    g.innerHTML = nodesFor().map((key, i) => {
      const fan = key === "score";
      return `<section class="pnode${fan ? " fan" : ""}" id="n-${key}" data-key="${key}">
        <div class="pn-head"><span class="pn-num">${i + 1}</span><h3>${esc(TITLES[key])}</h3><span class="pn-state" id="st-${key}"></span></div>
        <p class="pn-blurb">${BLURB[key] || ""}</p>
        ${controls(key) ? `<div class="pn-ctrls">${controls(key)}</div>` : ""}
        <div class="pn-out" id="out-${key}"></div>
        ${fan ? `<div class="pn-methods" id="methods"></div>` : codeBtns(key)}
      </section>`;
    }).join("");
  }

  // ------------------------------------------------------------ running
  function setState(keys, text) {
    for (const k of keys) { const e = $("st-" + k); if (e) e.textContent = text; }
  }
  const downstream = from => { const l = nodesFor(); return l.slice(Math.max(0, l.indexOf(from))); };

  async function loadPool(kind) {
    if (pools[kind]) return pools[kind];
    setState(["data"], "downloading…");
    await tick();
    const d = SOURCES[kind].dir;
    const [meta, buf] = await Promise.all([
      fetch(`data/${d}/meta.json`).then(r => r.json()),
      fetch(`data/${d}/counts.bin`).then(r => r.arrayBuffer()),
    ]);
    pools[kind] = S.decodePool(buf, meta);
    pools[kind].bytes = buf.byteLength;
    return pools[kind];
  }

  let running = false, again = null;
  async function run(level) {
    if (running) { again = again === "sim" || level === "sim" ? "sim" : level; return; }
    running = true;
    closePanel();
    try {
      const src = SOURCES[cfg.source];
      if (level === "sim" || !sim) {
        setState(downstream("data").slice(1), "…");
        if (src.paper) {
          poolNow = await loadPool(cfg.source);
          setState(["data"], "");
          await tick();
          sim = S.simulateRWSD(poolNow, { scenario: cfg.scenario, cellsPerGroup: cfg.cells, seed: cfg.seed,
            boostNull: !!cfg.boostNull, threshold: src.threshold, foldChange: cfg.fold, nAdded: cfg.nAdded });
          sim.kind = cfg.source;
        } else {
          poolNow = null;
          await tick();
          sim = S.simulate(Object.assign({ scenario: cfg.scenario, cellsPerGroup: Math.min(cfg.cells, cfg.source === "sparse" ? 200 : 500),
            seed: cfg.seed, inject: cfg.inject, foldChange: cfg.fold }, src.toy));
        }
        describeSim();
        setState(["prepare"], "ranking…");
        await tick();
        const t0 = performance.now();
        prep = S.prepare(sim);
        $("out-prepare").innerHTML = `<div class="kv"><b>${((performance.now() - t0) / 1000).toFixed(1)} s</b><span>to rank ${sim.G.toLocaleString()} genes in each of ${sim.N.toLocaleString()} cells (×3 ways)</span></div>`;
        setState(["prepare"], "");
      }
      newSet();
      await tick();
      scoreNow();
      $("out-recovery").innerHTML = "";
      $("out-sweep").innerHTML = "";
      setState(nodesFor(), "");
    } finally {
      running = false;
      if (again) { const l = again; again = null; run(l); }
    }
  }

  function kv(items) { return `<div class="kvs">${items.map(([b, s]) => `<div class="kv"><b>${b}</b><span>${s}</span></div>`).join("")}</div>`; }

  function describeSim() {
    const src = SOURCES[cfg.source];
    const ds = S.dataSummary(sim);
    if (src.paper) {
      const P = poolNow;
      $("out-data").innerHTML = kv([
        [P.nCells.toLocaleString(), "cells in the pool"], [P.nGenes.toLocaleString(), "genes"],
        [(P.bytes / 1e6).toFixed(1) + " MB", "packed raw counts"],
      ]) + `<p class="pn-note">${cfg.source === "rwsd"
        ? `GSE164381: naive B cells (our marker-based reconstruction) from the 4 healthy donors — the input the authors' script names "Naive Normal B cell". Raw 10x counts, every gene.`
        : `4,000 cells from Splatter (splatEstimate on the naive-B counts, then splatSimulate), as the authors build SSD. Genes are synthetic (Gene1, …).`}</p>`;
      const donors = {};
      sim.cellInfo.forEach(c => { donors[c.donor] = (donors[c.donor] || 0) + 1; });
      $("out-draw").innerHTML = kv([[`${sim.n} + ${sim.n}`, "control + treatment cells"],
        [Object.keys(donors).length, "donors represented"]]) +
        `<p class="pn-note">${Object.entries(donors).map(([d, k]) => `${esc(d)}: ${k}`).join(" · ")}</p>`;
      const libs = sim.cellInfo.map(c => { let s = 0; for (let k = P.ptr[c.pool]; k < P.ptr[c.pool + 1]; k++) s += P.ct[k]; return s; }).sort((a, b) => a - b);
      $("out-norm").innerHTML = kv([[libs[Math.floor(libs.length / 2)].toLocaleString(), "median counts per cell"],
        [`${libs[0].toLocaleString()}–${libs[libs.length - 1].toLocaleString()}`, "range"]]);
      $("out-filter").innerHTML = kv([[sim.G0.toLocaleString(), `genes kept of ${P.nGenes.toLocaleString()}`],
        [`${(100 * ds.sparsity).toFixed(0)}%`, "of the matrix is zero"], [`${ds.medianPctCellsExpressingAll.toFixed(1)}%`, "median gene detection"]]);
      const top = sim.signal.slice().sort((a, b) => sim.zeroFrac[a] - sim.zeroFrac[b]).slice(0, 8).map(g => sim.geneNames[g]);
      const rule = cfg.scenario === 1 || cfg.scenario === 3 ? `detected in >${Math.round(100 * (1 - src.threshold))}% of drawn cells` : "550 random detected genes";
      $("out-signal").innerHTML = kv([[sim.signal.length, rule], [`${ds.medianPctCellsExpressingSignal.toFixed(1)}%`, "their median detection"]]) +
        `<p class="pn-note">Most detected: ${top.map(esc).join(", ")}${cfg.scenario <= 2 ? "" : " — in S3/S4 these are the gene-set pool, not a change"}</p>`;
      $("out-boost").innerHTML = sim.boost === 1
        ? kv([["none", "no gene is changed in this scenario"]])
        : kv([[`×${sim.boost}`, `on ${sim.signal.length} genes × ${sim.n} treatment cells`]]);
      $("out-added").innerHTML = cfg.scenario >= 3
        ? kv([[sim.nCS, "genes added (zero in control)"], [`${ds.genesPerCell[0].toFixed(0)} → ${ds.genesPerCell[1].toFixed(0)}`, "genes detected per cell, control → treatment"]])
        : `<p class="pn-note">Skipped in S1/S2.</p>`;
      document.getElementById("n-added").classList.toggle("off", cfg.scenario < 3);
    } else {
      $("out-data").innerHTML = `<p class="pn-note">No file: the next node generates made-up counts.</p>`;
      $("out-toy").innerHTML = kv([[sim.G.toLocaleString(), "genes"], [sim.N.toLocaleString(), "cells"],
        [`${(100 * ds.sparsity).toFixed(0)}%`, "zeros"], [cfg.scenario <= 2 ? ds.nSignal : sim.nCS, cfg.scenario <= 2 ? "signal genes" : "treatment-only genes"]]);
    }
  }

  function newSet() {
    const spec = (sim.paper || sim.opt.scenario <= 2) ? { size: cfg.size, noise: cfg.noise } : { size: cfg.size, withCS: cfg.cs > 0, csFrac: cfg.cs };
    set = S.makeSet(sim, spec, S.makeRng(cfg.setSeed * 7919 + cfg.size));
    const cs = new Set(sim.csGenes);
    const nSig = set.filter(g => sim.isSignal[g]).length, nCS = set.filter(g => cs.has(g)).length;
    $("out-geneset").innerHTML = kv([[set.length, "genes"], [nSig, sim.boost === 1 && sim.paper ? "from the dense/random pool (unchanged)" : "signal (changed)"],
      [nCS, "treatment-only"], [set.length - nSig - nCS, "other, unchanged"]]) +
      `<p class="pn-note">${set.slice(0, 12).map(g => esc(sim.geneNames[g])).join(", ")}${set.length > 12 ? ", …" : ""}</p>`;
  }

  function setIsNull() {
    const cs = new Set(sim.csGenes);
    const hasCS = set.some(g => cs.has(g));
    const boosted = sim.paper ? sim.boost !== 1 : sim.opt.scenario <= 2;
    return !hasCS && !(boosted && set.some(g => sim.isSignal[g]));
  }

  function strip(scores, group) {
    let lo = Infinity, hi = -Infinity;
    for (const v of scores) { if (v < lo) lo = v; if (v > hi) hi = v; }
    const span = hi - lo || 1, W = 400, pad = 6;
    const x = v => pad + (W - 2 * pad) * (v - lo) / span;
    const ctl = [], trt = [];
    scores.forEach((v, i) => (group[i] ? trt : ctl).push([v, i]));
    const med = a => { const s = a.map(e => e[0]).sort((p, q) => p - q); return s[Math.floor(s.length / 2)]; };
    const pick = (a, n) => a.length <= n ? a : a.filter((_, i) => i % Math.ceil(a.length / n) === 0);
    const jit = i => ((i * 2654435761) % 1000) / 1000;
    let h = `<svg viewBox="0 0 ${W} 46" preserveAspectRatio="none">`;
    pick(ctl, 160).forEach(([v, c]) => { h += `<circle class="dot-c cell" data-c="${c}" cx="${x(v).toFixed(1)}" cy="${(6 + 12 * jit(c)).toFixed(1)}" r="2.2"/>`; });
    pick(trt, 160).forEach(([v, c]) => { h += `<circle class="dot-t cell" data-c="${c}" cx="${x(v).toFixed(1)}" cy="${(27 + 12 * jit(c + 7)).toFixed(1)}" r="2.2"/>`; });
    const mc = x(med(ctl)), mt = x(med(trt));
    return h + `<line class="med" x1="${mc}" x2="${mc}" y1="3" y2="21"/><line class="med" x1="${mt}" x2="${mt}" y1="24" y2="43"/></svg>`;
  }

  function scoreNow() {
    res = S.scoreAll(sim, prep, set);
    const nul = setIsNull();
    $("methods").innerHTML = `<div class="legend2"><span><i style="background:var(--muted)"></i>control</span><span><i style="background:var(--s1)"></i>treatment</span><span>tick = median</span></div>` +
      METHODS.map(m => {
        const r = res[m], hit = r.test.p < 0.05, arrow = r.test.z > 0 ? "↑" : r.test.z < 0 ? "↓" : "";
        return `<div class="mnode"><div class="mrow"><div class="name">${m}<button type="button" class="codechip" data-code="${m}">code</button></div>${strip(r.scores, sim.group)}` +
          `<div class="p">${hit ? `<span class="pill ${nul ? "bad" : "good"}">p ${pfmt(r.test.p)} ${arrow}</span>` : `p ${pfmt(r.test.p)}`}</div></div></div>`;
      }).join("");
    const fired = METHODS.filter(m => res[m].test.p < 0.05);
    $("out-test").innerHTML = kv([[`${fired.length} / 7`, "methods call a difference (p < 0.05)"],
      [nul ? "no" : "yes", "does this set contain a real difference?"]]) +
      `<p class="pn-note">${nul
        ? (fired.length ? `Every call here is a false positive: ${fired.join(", ")}.` : "No false positives on this set.")
        : (fired.length < 7 ? `Missed by: ${METHODS.filter(m => !fired.includes(m)).join(", ")}.` : "Every method detects it.")} Arrows give the direction; the test is two-sided.</p>`;
  }

  async function recoveryRun() {
    const out = $("out-recovery");
    out.innerHTML = `<p class="pn-note">Scoring ${cfg.nSets} gene sets × 7 methods…</p>`;
    await tick();
    const spec = (sim.paper || sim.opt.scenario <= 2) ? { size: cfg.size, noise: cfg.noise } : { size: cfg.size, withCS: cfg.cs > 0, csFrac: cfg.cs };
    const r = S.recovery(sim, prep, spec, cfg.nSets, 1000 + cfg.setSeed);
    const nul = sim.paper ? (cfg.scenario <= 2 ? cfg.noise === 1 : sim.boost === 1 && cfg.noise === 0) : (cfg.scenario <= 2 ? cfg.noise === 1 : cfg.cs === 0);
    out.innerHTML = `<p class="pn-note">${nul ? "These sets contain no real difference: every identification is a false positive." : "Rate of identifying the sets."}</p>` +
      `<table class="rec">` + METHODS.map(m => {
        const v = Math.round(100 * r[m]);
        const cls = nul ? (v <= 10 ? "good" : v <= 40 ? "mid" : "bad") : (v >= 80 ? "good" : v >= 40 ? "mid" : "bad");
        return `<tr><td>${m}</td><td class="num"><span class="pill ${cls}">${v}%</span></td><td class="bar"><span style="width:${v}%;background:var(--${cls})"></span></td></tr>`;
      }).join("") + `</table>`;
  }

  async function sweepRun() {
    const out = $("out-sweep");
    const sizes = [10, 25, 50, 100], acc = {};
    METHODS.forEach(m => { acc[m] = []; });
    const spec = size => (sim.paper || sim.opt.scenario <= 2) ? { size, noise: cfg.noise } : { size, withCS: cfg.cs > 0, csFrac: cfg.cs };
    for (const size of sizes) {
      out.innerHTML = `<p class="pn-note">Size ${size}…</p>`;
      await tick();
      const r = S.recovery(sim, prep, spec(size), 20, 1000 + size);
      METHODS.forEach(m => acc[m].push(r[m]));
    }
    const W = 520, H = 210, L = 40, B = 30, T = 10, R = 120;
    const X = i => L + i * (W - L - R) / (sizes.length - 1), Y = v => T + (1 - v) * (H - T - B);
    const colors = ["#3b7bbf", "#d0762b", "#3d9a8b", "#8a5cc2", "#b23b2a", "#b07d12", "#2e7d4f"];
    let h = `<svg class="sweep" viewBox="0 0 ${W} ${H}"><line class="ax" x1="${L}" y1="${Y(0)}" x2="${W - R}" y2="${Y(0)}"/><line class="ax" x1="${L}" y1="${Y(0)}" x2="${L}" y2="${Y(1)}"/>`;
    [0, 0.5, 1].forEach(v => { h += `<text x="${L - 32}" y="${Y(v) + 4}">${v * 100}%</text>`; });
    sizes.forEach((s, i) => { h += `<text x="${X(i) - 6}" y="${H - 12}">${s}</text>`; });
    METHODS.forEach((m, k) => {
      h += `<polyline fill="none" stroke="${colors[k]}" stroke-width="2" points="${acc[m].map((v, i) => `${X(i).toFixed(1)},${Y(v).toFixed(1)}`).join(" ")}"/>`;
      h += `<text x="${W - R + 8}" y="${T + 12 + k * 17}" style="fill:${colors[k]};font-weight:600">${m}</text>`;
    });
    out.innerHTML = h + `</svg><p class="pn-note">x: gene set size · y: share of 20 sets identified.</p>`;
  }

  // ------------------------------------------------------------ one cell
  function openCell(c) {
    if (!res) return;
    document.querySelectorAll("circle.cell.sel").forEach(e => e.classList.remove("sel"));
    document.querySelectorAll(`circle.cell[data-c="${c}"]`).forEach(e => e.classList.add("sel"));
    const G = sim.G, N = sim.N, trt = sim.group[c] === 1;
    const row = g => sim.data[c * G + g];
    const mean = (g, grp) => { let s = 0, k = 0; for (let i = 0; i < N; i++) if (sim.group[i] === grp) { s += sim.data[i * G + g]; k++; } return s / k; };
    const cs = new Set(sim.csGenes);
    const fold = sim.paper ? sim.boost : sim.opt.foldChange;
    const boostedNow = sim.paper ? sim.boost !== 1 : sim.opt.scenario <= 2;
    const role = g => cs.has(g) ? "added" : sim.isSignal[g] ? "signal" : "non-signal";
    const order = { signal: 0, added: 1, "non-signal": 2 };
    const genes = set.slice().sort((a, b) => (order[role(a)] - order[role(b)]) || (row(b) - row(a)));
    let detected = 0;
    for (let g = 0; g < G; g++) if (row(g) > 0) detected++;
    let vmax = 0;
    const rows = genes.map(g => { const r = [row(g), mean(g, 0), mean(g, 1)]; vmax = Math.max(vmax, ...r); return r; });
    const hm = v => `background:color-mix(in srgb, var(--s1) ${Math.round(100 * Math.min(1, v / (vmax || 1)))}%, transparent)`;
    let h = `<button class="cp-close" type="button" aria-label="Close">×</button><div class="kicker">${trt ? "Treatment" : "Control"} cell #${c + 1}</div>`;
    if (sim.paper) {
      const ci = sim.cellInfo[c];
      h += `<h3>${esc(ci.donor)} · <code>${esc(ci.barcode)}</code></h3>`;
      if (ci.addedFrom !== null && sim.nCS) h += `<p class="cp-note">Its ${sim.nCS} added genes were copied from another cell (pool #${ci.addedFrom + 1}).</p>`;
    } else h += `<h3>Simulated cell</h3>`;
    h += kv([[detected.toLocaleString(), `genes detected (of ${G.toLocaleString()})`], [`${genes.filter(g => row(g) > 0).length} / ${set.length}`, "gene-set genes detected here"]]);
    h += `<h4>Score under each method</h4><table class="cp-scores"><tr><th>Method</th><th class="num">score</th><th class="num">rank of ${N}</th></tr>`;
    for (const m of METHODS) {
      const s = res[m].scores, v = s[c];
      let above = 0;
      for (let i = 0; i < N; i++) if (s[i] > v) above++;
      h += `<tr><td>${m}</td><td class="num">${Number(v).toPrecision(3)}</td><td class="num">${above + 1}</td></tr>`;
    }
    h += `</table><h4>The gene set in this cell</h4><p class="cp-note">Log-normalised values; darker = higher. <span class="role signal">signal</span> genes first` +
      (boostedNow ? ` (×${fold} in treatment cells; before → after shown)` : " (not changed in this run)") + `; <span class="role added">added</span> = treatment-only.</p>`;
    h += `<div class="table-wrap"><table class="cp-genes"><tr><th>Gene</th><th></th><th class="num">this cell</th><th class="num">control mean</th><th class="num">treatment mean</th></tr>`;
    genes.forEach((g, k) => {
      const [v, mc, mt] = rows[k], r = role(g);
      const pre = trt && boostedNow && sim.isSignal[g] && v > 0 ? `<small>${(v / fold).toFixed(2)}→</small>` : "";
      h += `<tr><td>${esc(sim.geneNames[g])}</td><td><span class="role ${r === "non-signal" ? "non" : r}">${r}</span></td>` +
        `<td class="num hm" style="${hm(v)}">${pre}${v.toFixed(2)}</td><td class="num hm" style="${hm(mc)}">${mc.toFixed(2)}</td><td class="num hm" style="${hm(mt)}">${mt.toFixed(2)}</td></tr>`;
    });
    h += `</table></div>`;
    const d = $("sidePanel");
    d.innerHTML = h;
    d.hidden = false;
    d.querySelector(".cp-close").addEventListener("click", closePanel);
  }

  // ------------------------------------------------------------ events
  function onChange(e) {
    const k = e.target.dataset.k;
    if (!k) return;
    const v = e.target.value;
    cfg[k] = k === "source" || k === "inject" ? v : Number(v);
    if (k === "source" || k === "scenario") render();
    if (k === "boostNull" || k === "fold") render();
    run(SIM_KEYS.has(k) ? "sim" : "set");
  }

  function init() {
    applyHash();
    render();
    const g = $("graph");
    g.addEventListener("change", onChange);
    g.addEventListener("click", e => {
      const code = e.target.closest("[data-code]");
      if (code) { showCode(code.dataset.code); return; }
      const act = e.target.closest("[data-act]");
      if (act) {
        if (act.dataset.act === "reroll") { cfg.setSeed++; newSet(); scoreNow(); }
        if (act.dataset.act === "recovery") recoveryRun();
        if (act.dataset.act === "sweep") sweepRun();
        return;
      }
      const dot = e.target.closest("circle.cell");
      if (dot) openCell(Number(dot.dataset.c));
    });
    document.addEventListener("keydown", e => { if (e.key === "Escape") closePanel(); });
    window.addEventListener("hashchange", () => { applyHash(); render(); run("sim").then(() => { if (PRESETS[location.hash.slice(1)]) recoveryRun(); }); });
    run("sim").then(() => { if (PRESETS[location.hash.slice(1)]) { recoveryRun(); $("n-recovery").scrollIntoView({ behavior: "smooth", block: "center" }); } });
  }

  const PRESETS = {
    "s1-null": { source: "rwsd", scenario: 1, cells: 200, size: 100, noise: 1 },
    "s2-signal": { source: "rwsd", scenario: 2, cells: 20, size: 100, noise: 0 },
    "s2-small": { source: "dense", scenario: 2, cells: 20, size: 25, noise: 0 },
    "s3-null": { source: "rwsd", scenario: 3, cells: 200, size: 100, noise: 0 },
    "s4-null": { source: "rwsd", scenario: 4, cells: 200, size: 100, noise: 0 },
  };
  function applyHash() {
    const p = PRESETS[location.hash.slice(1)];
    if (p) Object.assign(cfg, { fold: 1.2, boostNull: 0, nAdded: 250, inject: "log", seed: 7 }, p);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
