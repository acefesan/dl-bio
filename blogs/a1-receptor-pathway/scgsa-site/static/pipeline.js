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
  const OURS_RAW = "https://raw.githubusercontent.com/acefesan/dl-bio/main/blogs/a1-receptor-pathway/";
  const CODE = {
    build:    { files: [["paper-data/01_naive_b.R", "healthy donors → QC → clusters → marker table"],
                        ["paper-data/02_export_naive_b.R", "keep the naive-B clusters"],
                        ["paper-data/03_export_browser.py", "pack raw counts for the browser"],
                        ["paper-data/05_splatter_pool.R", "SSD: splatEstimate + splatSimulate"]],
                auth: [[RWSD, 29, 29, "their input: load(\"Naive Normal B cell.RData\") — not published"], [SSD, 19, 21, "SSD: splatEstimate on that object's counts"]] },
    data:     { ours: { find: "function decodePool", label: "decodePool()" } },
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
    // Fetch the exact (versioned) file this page runs, revalidating, so the
    // snippet can't come from a stale cached copy of an older scgsa.js.
    if (!srcText) {
      const tag = document.querySelector('script[src*="static/scgsa.js"]');
      srcText = await fetch(tag ? tag.src : "static/scgsa.js", { cache: "no-cache" }).then(r => r.text());
    }
    const lines = srcText.split("\n");
    const i = lines.findIndex(l => l.includes(ref.find));
    if (i < 0) return null;
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
    for (const [path, what] of c.files || []) {
      h += `<h4>Ours · <a href="${OURS}${path}" target="_blank" rel="noopener">${esc(path)}</a></h4><p class="cp-note">${esc(what)} — run offline, not in the browser.</p>`;
      try { h += `<pre class="code"><code>${esc(await fetch(OURS_RAW + path).then(r => r.ok ? r.text() : Promise.reject(r.status)))}</code></pre>`; }
      catch (e) { h += `<p class="cp-note">(Couldn't fetch it from GitHub; follow the link.)</p>`; }
    }
    if (c.ours) {
      const s = await ourSnippet(c.ours);
      h += s
        ? `<h4>Ours · <a href="${OURS}scgsa-site/static/scgsa.js#L${s.start}-L${s.end}" target="_blank" rel="noopener">static/scgsa.js L${s.start}–${s.end}</a></h4><pre class="code"><code>${esc(s.text)}</code></pre>`
        : `<h4>Ours · <code>${esc(c.ours.label)}</code></h4><p class="cp-note">Couldn't locate this in static/scgsa.js — please reload the page (an old copy may be cached). <a href="${OURS}scgsa-site/static/scgsa.js" target="_blank" rel="noopener">Open the file on GitHub</a>.</p>`;
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
    build: "Build the dataset", data: "Load the dataset", toy: "Simulate counts (toy model)", draw: "Draw cells and split into two groups",
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
    if (c.files) h += `<button type="button" data-code="${key}">ours: <code>paper-data/</code> scripts</button>`;
    if (c.auth) h += `<button type="button" data-code="${key}">authors' R: ${c.auth.map(a => `L${a[1]}`).join(", ")}</button>`;
    return h + `</div>`;
  };

  // ------------------------------------------------------------ graph layout
  function nodesFor() {
    const paper = !!SOURCES[cfg.source].paper;
    const list = ["build", "data"];
    if (paper) list.push("draw", "norm", "filter", "signal", "boost", "added");
    else list.push("toy");
    list.push("geneset", "prepare", "score", "test", "recovery", "sweep");
    return list;
  }

  function controls(key) {
    const paper = !!SOURCES[cfg.source].paper;
    switch (key) {
      case "build": return opt("source", Object.entries(SOURCES).map(([k, s]) => [k, s.label]), { label: "Dataset" });
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
      case "geneset": return opt("size", [10, 20, 30, 40, 50, 100, 200, 300, 400, 500].map(v => [v, String(v)]), { label: "Gene set size" }) +
        ((paper || cfg.scenario <= 2)
          ? opt("noise", [[0, "0% (all signal)"], [0.2, "20%"], [0.5, "50%"], [0.8, "80%"], [1, "100% (no signal)"]], { label: "Noise: non-signal share" })
          : opt("cs", [[0, "0% (null)"], [0.02, "2%"], [0.05, "5%"], [0.1, "10%"], [0.2, "20%"]], { label: "Treatment-only share" })) +
        `<button type="button" class="btn ghost" data-act="reroll">New random set</button>`;
      case "recovery": return opt("nSets", [[30, "30 (quick)"], [100, "100 (as the paper)"]], { label: "Gene sets" }) + `<button type="button" class="btn" data-act="recovery">Run</button>`;
      case "sweep": return `<button type="button" class="btn" data-act="sweep">Run sweep (sizes 10–100)</button>`;
      default: return "";
    }
  }

  const BLURB = {
    build: "Where the cells come from. These steps ran once, offline, in R and Python; everything from step 2 on runs live in your browser. Click any step to see it explained.",
    data: "Download the saved counts and unpack them.",
    draw: "Pick some cells at random and split them into two equal groups, \"control\" and \"treatment\". They're the same kind of cell, so any difference found later was put there on purpose by the next steps.",
    norm: "Put every cell on the same scale: some cells yield more molecules than others, so each gene's count is divided by the cell's total, then log-transformed.",
    filter: "Drop genes that none of the picked cells expressed — they can't differ between the groups.",
    signal: "Choose which genes the experiment changes (S1/S2), or which genes gene sets are drawn from when nothing is changed (S3/S4).",
    boost: "The planted difference: in treatment cells only, multiply those genes' log-normalised values by a factor. Nothing else changes.",
    added: "S3/S4 only: give treatment cells extra genes that control cells don't have at all, without changing any shared gene.",
    geneset: "Pick a list of genes to test: some share from the changed genes, the rest (\"noise\") from everything else.",
    prepare: "Sort each cell's genes from highest to lowest, once; most methods work from these ranks.",
    score: "Each of the 7 methods boils the gene set down to one number per cell. Click any dot to open that cell.",
    test: "Ask, per method: do treatment cells score differently from control cells? (Wilcoxon test.)",
    recovery: "Repeat with many random gene sets and count how often each method says \"different\". Where the paper ran the same setting, its published number is shown beside ours.",
    sweep: "The same, at several gene set sizes, to see whether bigger sets make mistakes more or less likely.",
    toy: "Made-up counts from a simple statistical model, with the scenario's change built in. Good for seeing mechanisms; not the paper's data.",
  };


  function render() {
    const g = $("graph");
    g.innerHTML = nodesFor().map((key, i) => {
      const fan = key === "score";
      return `<section class="pnode${fan ? " fan" : ""}" id="n-${key}" data-key="${key}">
        <div class="pn-head"><span class="pn-num">${i + 1}</span><h3>${esc(TITLES[key])}</h3><span class="pn-state" id="st-${key}"></span><button type="button" class="explainbtn" data-explain="${key}">explain ›</button></div>
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

  // ------------------------------------------------------------ dataset provenance + published results
  let PROV = null, PUB = null, QC = null, UMAP = null;
  const getJSON = u => fetch(u).then(r => r.ok ? r.json() : null).catch(() => null);
  async function loadSide() {
    if (!PROV) PROV = await getJSON("data/provenance.json");
    if (!QC) QC = await getJSON("data/provenance_qc.json");
    if (!PUB) PUB = await getJSON("data/authors_results.json");
  }

  function describeBuild() {
    const out = $("out-build");
    if (!out) return;
    if (!SOURCES[cfg.source].paper) {
      out.innerHTML = `<p class="pn-note">No dataset: the toy simulator makes up negative-binomial counts for ${SOURCES[cfg.source].toy.genes.toLocaleString()} genes, each with a random mean (log-normal, μ = ${SOURCES[cfg.source].toy.meanLog}, σ = ${SOURCES[cfg.source].toy.sdLog}) and dispersion 0.2.</p>`;
      return;
    }
    if (!PROV) { out.innerHTML = `<p class="pn-note">(provenance.json not found)</p>`; return; }
    const P = PROV, steps = [
      ["b-geo", `${P.geo.samples.length} samples · ${P.geo.cells.toLocaleString()} cells`, `Downloaded from GEO (${P.geo.accession}): one file per person, a table of how many RNA molecules of each gene were counted in each cell. B cells were purified from each blood sample before sequencing, so most cells are B cells. 4 healthy people and 7 with mild COVID-19 (Hoehn et al. 2021).`],
      ["b-healthy", `4 healthy donors · ${P.healthy_cells.toLocaleString()} cells`, `Keep only the 4 healthy people. The authors' input is named "Naive <em>Normal</em> B cell" — "Normal" is GEO's word for healthy.`],
      ["b-qc", `Quality control · ${P.qc.cells_after.toLocaleString()} cells`, `Throw out what can't be trusted: genes seen in fewer than 3 cells; "cells" with fewer than 200 genes (empty droplets); cells with 15% or more of their RNA from mitochondria (dying cells); cells with 4,000+ genes (probably two cells stuck together).`],
      ["b-cluster", `Group similar cells · ${P.clustering.n_clusters} clusters`, `Without telling the computer what any cell is: put cells on a common scale, keep the 2,000 genes that differ most between cells, compress them into 20 summary axes (PCA), link each cell to its most similar cells, and find the tightly linked groups (Louvain).`],
      ["b-naive", `naive B · ${P.naive_b.cells.toLocaleString()} cells × ${P.naive_b.genes.toLocaleString()} genes`, `Name each group by its marker genes and keep the groups that look like naive B cells: high ${P.clustering.naive_markers.slice(0, 3).join(", ")}…, low memory-B markers (${P.clustering.memory_markers.slice(0, 2).join(", ")}…), and no sign of T, NK or plasma cells.`],
    ];
    if (cfg.source === "ssd") steps.push(["b-splatter", `Splatter · ${P.splatter.pool_cells.toLocaleString()} synthetic cells`, `Fit a statistical model to those naive B cells (how common each gene is, how many molecules a cell yields, how noisy genes are), then generate ${P.splatter.pool_cells.toLocaleString()} made-up cells from it. Gene names are lost: Gene1, Gene2, …`]);
    steps.push(["b-pack", cfg.source === "ssd" ? "packed · data/ssd/" : "packed · data/naive_b/", "Save the kept cells' counts (every gene, nothing averaged) in a compact file the browser can read."]);
    let h = `<p class="pn-note">Click any step for the full story.</p><ol class="chain">${steps.map(([k, b, t]) => `<li data-explain="${k}" tabindex="0" role="button"><b>${b} <span class="more">explain ›</span></b><span>${t}</span></li>`).join("")}</ol>`;
    const M = P.clustering.markers_shown;
    h += `<details><summary>Cluster marker table (mean normalised expression)</summary><div class="table-wrap"><table class="clus"><tr><th>cluster</th><th class="num">cells</th><th class="num">naive − memory</th>${M.map(m => `<th class="num">${m}</th>`).join("")}</tr>` +
      P.clustering.clusters.map(c => `<tr class="${c.kept ? "kept" : ""}"><td>${c.cluster}${c.kept ? " ✓" : ""}</td><td class="num">${c.cells.toLocaleString()}</td><td class="num">${c.naive_minus_memory}</td>${M.map(m => `<td class="num">${c[m]}</td>`).join("")}</tr>`).join("") +
      `</table></div></details>`;
    h += `<p class="pn-note"><b>Why rebuilt:</b> the paper says the cells were "normalized, clustered and annotated as described in the original paper", and its data statement points only to GEO ${P.geo.accession}. We checked the paper, its supplement and its Zenodo archive (doi:10.5281/zenodo.13620619): code and result tables, no cell labels or processed object. So the naive-B subset is our marker-based reconstruction. It reproduces the paper's own numbers (≈485 dense genes, 43% median detection vs. their "at least 500" and 44.5%).</p>`;
    out.innerHTML = h;
  }

  // The authors' published recovery for the current setting, if they ran it.
  function published(size, noise) {
    if (!PUB || !SOURCES[cfg.source].paper) return null;
    if (cfg.fold !== 1.2 || cfg.nAdded !== 250 || (cfg.scenario >= 3 && cfg.boostNull)) return null;
    const d = PUB.data[cfg.source];
    const node = d && d[String(cfg.cells)] && d[String(cfg.cells)][String(cfg.scenario)];
    return (node && node[String(size)] && node[String(size)][Number(noise).toFixed(1)]) || null;
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
      ]) + (cfg.source === "ssd" ? `<p class="pn-note">Splatter genes are synthetic (Gene1, …), not real genes.</p>` : "");
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
    const pub = published(cfg.size, cfg.noise);
    const grade = v => nul ? (v <= 10 ? "good" : v <= 40 ? "mid" : "bad") : (v >= 80 ? "good" : v >= 40 ? "mid" : "bad");
    out.innerHTML = `<p class="pn-note">${nul ? "These sets contain no real difference: every identification is a false positive." : "Rate of identifying the sets."}</p>` +
      `<table class="rec"><tr><th>Method</th><th class="num">ours (${cfg.nSets} sets)</th><th></th>${pub ? `<th class="num">authors' published</th>` : ""}</tr>` + METHODS.map(m => {
        const v = Math.round(100 * r[m]), cls = grade(v);
        const a = pub && pub[m];
        return `<tr><td>${m}</td><td class="num"><span class="pill ${cls}">${v}%</span></td><td class="bar"><span style="width:${v}%;background:var(--${cls})"></span></td>` +
          (pub ? `<td class="num">${a ? `<span class="pill ${grade(a[0])}">${a[0]}%</span> <small>± ${a[1]}</small>` : "—"}</td>` : "") + `</tr>`;
      }).join("") + `</table>` +
      `<p class="pn-note">${pub ? `Authors' published: mean ± SD over their ${Object.values(pub)[0][2]} replicates of 100 gene sets each, same data type, cells per group, scenario, size and noise (Zenodo, "Two side dat matrix.csv"). Ours: one replicate of ${cfg.nSets} sets.`
        : "No published result for this exact setting (the paper ran sizes 10–500 at 20/50/200/500 cells for S1/S2, and 200 cells for S3/S4, with ×1.2 and 250 added genes)."}</p>`;
  }

  async function sweepRun() {
    const out = $("out-sweep");
    const sizes = [10, 20, 50, 100], acc = {};
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
    let anyPub = false;
    METHODS.forEach((m, k) => {
      const pubPts = sizes.map((sz, i) => { const p = published(sz, cfg.noise); return p && p[m] ? `${X(i).toFixed(1)},${Y(p[m][0] / 100).toFixed(1)}` : null; }).filter(Boolean);
      if (pubPts.length > 1) { anyPub = true; h += `<polyline fill="none" stroke="${colors[k]}" stroke-width="1.5" stroke-dasharray="4 3" opacity="0.8" points="${pubPts.join(" ")}"/>`; }
      h += `<polyline fill="none" stroke="${colors[k]}" stroke-width="2" points="${acc[m].map((v, i) => `${X(i).toFixed(1)},${Y(v).toFixed(1)}`).join(" ")}"/>`;
      h += `<text x="${W - R + 8}" y="${T + 12 + k * 17}" style="fill:${colors[k]};font-weight:600">${m}</text>`;
    });
    out.innerHTML = h + `</svg><p class="pn-note">x: gene set size · y: share of sets identified. Solid: ours (20 sets per size).${anyPub ? " Dashed: the authors' published means for the same setting." : ""}</p>`;
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

  // ------------------------------------------------------------ explanations
  // Every build step and pipeline node opens a panel: what it does from first
  // principles, the real numbers, and (for the in-browser nodes) a small
  // "sample in → sample out": the same few cells × genes before and after.
  const n0 = v => Number(v).toLocaleString();
  const f2 = v => (v === 0 ? "0" : Number(v).toFixed(2));

  function histSVG(h, lines, xlab) {
    const W = 460, H = 150, L = 34, B = 26, T = 8, R = 8;
    const n = h.counts.length, max = Math.max(...h.counts) || 1;
    const x0 = h.breaks[0], x1 = h.breaks[h.breaks.length - 1];
    const X = v => L + (v - x0) / (x1 - x0) * (W - L - R), bw = (W - L - R) / n;
    let s = `<svg class="hist" viewBox="0 0 ${W} ${H}">`;
    h.counts.forEach((c, i) => {
      const hh = (H - T - B) * c / max;
      s += `<rect x="${(L + i * bw).toFixed(1)}" y="${(H - B - hh).toFixed(1)}" width="${Math.max(0.5, bw - 0.6).toFixed(1)}" height="${hh.toFixed(1)}" class="bar"/>`;
    });
    for (const [v, label] of lines) {
      s += `<line x1="${X(v)}" x2="${X(v)}" y1="${T}" y2="${H - B}" class="cut"/><text x="${X(v) + 3}" y="${T + 10}" class="cutl">${label}</text>`;
    }
    s += `<line x1="${L}" x2="${W - R}" y1="${H - B}" y2="${H - B}" class="ax"/>`;
    for (let k = 0; k <= 4; k++) { const v = x0 + k * (x1 - x0) / 4; s += `<text x="${X(v)}" y="${H - 10}" text-anchor="middle">${n0(Math.round(v))}</text>`; }
    return s + `<text x="${(L + W) / 2}" y="${H}" text-anchor="middle">${xlab}</text></svg>`;
  }

  // small matrix table: rows = genes, cols = cells
  function grid(title, cells, genes, val, opts) {
    opts = opts || {};
    let h = `<div class="sio"><div class="sio-t">${title}</div><div class="table-wrap"><table class="sio-tab"><tr><th></th>` +
      cells.map(c => `<th class="num">${c.label}</th>`).join("") + `</tr>`;
    for (const g of genes) {
      h += `<tr><td>${esc(g.label)}${g.tag ? ` <span class="role ${g.tagcls || "non"}">${g.tag}</span>` : ""}</td>` +
        cells.map(c => { const v = val(c, g); const hi = opts.mark && opts.mark(c, g); return `<td class="num${hi ? " chg" : ""}">${v === null ? "—" : (opts.fmt || f2)(v)}</td>`; }).join("") + `</tr>`;
    }
    return h + `</table></div></div>`;
  }
  const sio = (a, b, note) => `<h4>Sample in → sample out</h4><div class="sio-row">${a}<div class="sio-arrow">→</div>${b}</div>${note ? `<p class="cp-note">${note}</p>` : ""}`;

  // a fixed small sample: 2 control + 2 treatment cells, and genes of each role
  function sample() {
    const cells = [0, 1, sim.n, sim.n + 1].map(i => ({ i, label: `${sim.group[i] ? "T" : "C"}${i < sim.n ? i + 1 : i - sim.n + 1}` }));
    const G = sim.G, has = g => cells.some(c => sim.data[c.i * G + g] > 0);
    const pick = (list, k) => list.filter(has).slice(0, k);
    const sig = pick(sim.signal, 2);
    const non = pick(sim.nonSignal.filter(g => !sim.csGenes.includes(g)), 2);
    const cs = pick(sim.csGenes, 1);
    const genes = [...sig.map(g => ({ g, label: sim.geneNames[g], tag: sim.paper && sim.boost === 1 ? "pool" : "signal", tagcls: "signal" })),
      ...cs.map(g => ({ g, label: sim.geneNames[g], tag: "added", tagcls: "added" })),
      ...non.map(g => ({ g, label: sim.geneNames[g] }))];
    return { cells, genes };
  }
  const rawCount = (c, gName) => {
    if (!poolNow) return null;
    const pc = sim.cellInfo[c.i].pool, gi = poolNow.genes.indexOf(gName);
    for (let k = poolNow.ptr[pc]; k < poolNow.ptr[pc + 1]; k++) if (poolNow.gi[k] === gi) return poolNow.ct[k];
    return 0;
  };
  const libSize = c => { const pc = sim.cellInfo[c.i].pool; let s = 0; for (let k = poolNow.ptr[pc]; k < poolNow.ptr[pc + 1]; k++) s += poolNow.ct[k]; return s; };
  const isBoosted = (c, g) => sim.group[c.i] === 1 && sim.isSignal[g.g] && (sim.paper ? sim.boost !== 1 : sim.opt.scenario <= 2);

  async function umapPanel(colorBy) {
    if (!UMAP) UMAP = await getJSON("data/provenance_umap.json");
    const box = $("umapBox");
    if (!UMAP || !box) return;
    const cv = box.querySelector("canvas"), ctx = cv.getContext("2d");
    const css = getComputedStyle(document.documentElement);
    const W = cv.width = box.clientWidth * 2, H = cv.height = Math.round(W * 0.8);
    cv.style.height = `${H / 2}px`;
    const xs = UMAP.x, ys = UMAP.y;
    const mnx = Math.min(...xs), mxx = Math.max(...xs), mny = Math.min(...ys), mxy = Math.max(...ys);
    const X = v => 12 + (v - mnx) / (mxx - mnx) * (W - 24), Y = v => H - 12 - (v - mny) / (mxy - mny) * (H - 24);
    ctx.clearRect(0, 0, W, H);
    const cat = (i, n) => `hsl(${Math.round(i * 360 / n)}, 60%, 52%)`;
    const kept = new Set(UMAP.kept.map(Number));
    let color, legend = "";
    if (colorBy === "cluster") {
      color = i => cat(UMAP.cluster[i], 14);
      legend = Array.from({ length: 14 }, (_, k) => `<span><i style="background:${cat(k, 14)}"></i>${k}${kept.has(k) ? "✓" : ""}</span>`).join("");
    } else if (colorBy === "kept") {
      const on = css.getPropertyValue("--s3").trim(), off = css.getPropertyValue("--line").trim();
      color = i => kept.has(UMAP.cluster[i]) ? on : off;
      legend = `<span><i style="background:${on}"></i>kept as naive B</span><span><i style="background:${off}"></i>dropped</span>`;
    } else if (colorBy === "donor") {
      color = i => cat(UMAP.donor[i], 4);
      legend = UMAP.donors.map((d, k) => `<span><i style="background:${cat(k, 4)}"></i>${d}</span>`).join("");
    } else {
      const v = UMAP.expr[colorBy], mx = Math.max(...v) || 1, hot = css.getPropertyValue("--s1").trim(), cold = css.getPropertyValue("--line").trim();
      color = i => v[i] > 0 ? `color-mix(in srgb, ${hot} ${Math.round(25 + 75 * v[i] / mx)}%, ${cold})` : cold;
      legend = `<span><i style="background:${cold}"></i>0</span><span><i style="background:${hot}"></i>high ${colorBy}</span>`;
    }
    const order = UMAP.x.map((_, i) => i);
    if (colorBy !== "cluster" && colorBy !== "donor") order.sort((a, b) => (colorBy === "kept" ? Number(kept.has(UMAP.cluster[a])) - Number(kept.has(UMAP.cluster[b])) : UMAP.expr[colorBy][a] - UMAP.expr[colorBy][b]));
    for (const i of order) { ctx.fillStyle = color(i); ctx.beginPath(); ctx.arc(X(xs[i]), Y(ys[i]), 2.6, 0, 7); ctx.fill(); }
    // cluster labels at centroids
    if (colorBy === "cluster" || colorBy === "kept") {
      ctx.font = "bold 22px system-ui"; ctx.fillStyle = css.getPropertyValue("--ink").trim(); ctx.textAlign = "center";
      for (let k = 0; k < 14; k++) {
        let sx = 0, sy = 0, m = 0;
        UMAP.cluster.forEach((c, i) => { if (c === k) { sx += X(xs[i]); sy += Y(ys[i]); m++; } });
        if (m) ctx.fillText(String(k), sx / m, sy / m);
      }
    }
    box.querySelector(".legend2").innerHTML = legend;
  }

  const EXPLAIN = {
    build: () => {
      if (!SOURCES[cfg.source].paper) return EXPLAIN.toy();
      const steps = [["b-geo", "Raw data from GEO"], ["b-healthy", "Healthy donors only"], ["b-qc", "Quality control"], ["b-cluster", "Clustering"], ["b-naive", "Which clusters are naive B cells"]]
        .concat(cfg.source === "ssd" ? [["b-splatter", "Splatter (SSD)"]] : []).concat([["b-pack", "Packing for the browser"]]);
      return `<p>The paper's input is "naive B cells from the healthy donors of GSE164381, annotated as in the original study". That object was never published, so it is rebuilt here from the raw GEO files, offline, once. Each step:</p>
        <ol class="steps">${steps.map(([k, t]) => `<li><button type="button" class="linkbtn" data-explain="${k}">${t} ›</button></li>`).join("")}</ol>
        <p class="cp-note">Check that it worked: on these cells the authors' own S1 rule picks ~485 genes with median detection 43%, against the paper's "at least 500" and 44.5%; and our recovery rates sit beside their published ones in the recovery node.</p>`;
    },
    "b-geo": () => {
      const P = PROV;
      return `<p>A 10x experiment puts each cell in a droplet with a barcoded bead, reads ~90 letters of each captured RNA, and counts how many distinct molecules (UMIs) of each gene each barcode caught. Cell Ranger's <b>filtered</b> matrix keeps only barcodes it judged to be real cells (not empty droplets). Each file is cells × <b>${n0(QC ? QC.genes_total : 36601)} genes</b>, values = UMI counts, mostly zero.</p>
        <p>GSE164381 is from Hoehn et al. 2021 (J Immunol), a study of B-cell receptor repertoires in mild vs severe COVID-19. Per GEO's protocol, B cells were <b>purified</b> from each thawed blood sample with a magnetic "pan-B negative selection" kit (antibodies pull out every cell that is <em>not</em> a B cell), and ~7,000 purified B cells per sample went into the 10x run (5′ gene expression plus a B-cell-receptor library).</p>
        <p>That is why the data is B-cell-heavy: B cells are only ~5–15% of normal blood white cells, but here, after QC, ~84% of the healthy donors' cells are B cells (the rest is T, NK and monocyte carry-over — negative selection removes most non-B cells, not all), and ~79% of those B cells are naive, a little above the usual 60–70% for adults, as expected for mostly young donors (23–31 years for three of the four). Net: 9,528 naive B cells = 55% of the 17,362 raw cells, 66% of those passing QC.</p>
        <table class="cp-scores"><tr><th>GEO sample</th><th>group</th><th class="num">cells</th></tr>${P.geo.samples.map(s => `<tr><td>${s.gsm}</td><td>${s.healthy ? "healthy" : "mild COVID-19"}</td><td class="num">${n0(s.cells)}</td></tr>`).join("")}</table>`;
    },
    "b-healthy": () => `<p>The authors' simulation script loads an object named <code>Naive Normal B cell</code>: "Normal" is GEO's label for the healthy donors ("disease state: Normal healthy"). So the paper's input excludes all 7 COVID-19 samples — none of its simulated "treatment" effects should be confused with disease.</p>
        <table class="cp-scores"><tr><th>donor</th><th>GEO</th><th>age, sex</th><th class="num">cells</th></tr>
        <tr><td>healthy donor 1</td><td>GSM5008746</td><td>53, M</td><td class="num">${n0(PROV.geo.samples[0].cells)}</td></tr><tr><td>healthy donor 2</td><td>GSM5008747</td><td>31, F</td><td class="num">${n0(PROV.geo.samples[1].cells)}</td></tr>
        <tr><td>healthy donor 3</td><td>GSM5008748</td><td>26, M</td><td class="num">${n0(PROV.geo.samples[2].cells)}</td></tr><tr><td>healthy donor 4</td><td>GSM5008749</td><td>23, M</td><td class="num">${n0(PROV.geo.samples[3].cells)}</td></tr></table>`,
    "b-qc": () => {
      if (!QC) return "<p>(provenance_qc.json missing)</p>";
      return `<p>Quality control removes barcodes that are not one healthy cell, and genes too rare to carry information. Four rules (Seurat <code>CreateSeuratObject(min.cells = 3, min.features = 200)</code>, then <code>subset(percent.mt &lt; 15 &amp; nFeature_RNA &lt; 4000)</code>):</p>
        <table class="cp-scores"><tr><th>rule</th><th>why</th><th class="num">removed</th></tr>
        <tr><td>gene detected in ≥ 3 cells</td><td>a gene seen in 0–2 cells can't vary in any useful way; ${n0(QC.genes_in_0_cells)} of them are never seen at all</td><td class="num">${n0(QC.genes_in_lt3_cells)} genes</td></tr>
        <tr><td>cell has ≥ 200 genes</td><td>very few genes = an empty or broken droplet that slipped past the filter</td><td class="num">${n0(QC.cells_lt200_genes)} cells</td></tr>
        <tr><td>mitochondrial reads &lt; 15%</td><td>a dying cell's membrane leaks its cytoplasmic mRNA but keeps the mitochondria, so its reads skew to MT- genes</td><td class="num">${n0(QC.cells_mt_ge15)} cells</td></tr>
        <tr><td>&lt; 4,000 genes</td><td>far more genes than a typical cell = probably two cells in one droplet (a doublet)</td><td class="num">${n0(QC.cells_genes_ge4000)} cells</td></tr></table>
        <p class="cp-note">${n0(QC.cells_total)} cells × ${n0(QC.genes_total)} genes in → ${n0(PROV.qc.cells_after)} cells × ${n0(QC.genes_total - QC.genes_in_lt3_cells)} genes out. Typical cell after QC: ${n0(QC.median_genes)} genes, ${n0(QC.median_umi)} UMIs, ${QC.median_pmt}% mitochondrial. (Seurat counts genes per cell after dropping the rare genes, which moves 3 cells across the lines; hence ${n0(PROV.qc.cells_after)} rather than ${n0(QC.cells_after)}.)</p>
        <h4>Genes detected per cell</h4>${histSVG(QC.hist_genes, [[200, "200"], [4000, "4,000"]], "genes detected per cell")}
        <h4>% of a cell's reads from mitochondrial genes</h4>${histSVG(QC.hist_pmt, [[15, "15%"]], "% mitochondrial")}
        <p class="cp-note">Thresholds are ours (the authors' are not published); 15% sits in the usual 10–20% range and cuts the long right tail.</p>`;
    },
    "b-cluster": () => `<p>Clustering groups cells by similar expression without being told what they are. Step by step:</p>
        <ol class="steps"><li><b>LogNormalize</b> each cell (log1p(count ÷ total × 10⁴)) so deep and shallow cells are comparable.</li>
        <li><b>2,000 variable genes</b>: genes that vary between cells more than their average level predicts. Housekeeping genes vary little; identity genes (IGHD, CD3E, …) vary a lot. Only these carry "which kind of cell" information.</li>
        <li><b>Scale</b> each gene to mean 0, SD 1, so a highly expressed gene doesn't dominate.</li>
        <li><b>PCA</b>: compress 2,000 correlated genes into 20 axes of shared variation; the first ~15 hold the structure, the rest mostly noise.</li>
        <li><b>Neighbour graph</b>: link each cell to the cells closest to it in PC 1–15 space, weighting links by how many neighbours they share.</li>
        <li><b>Louvain</b>: find groups of cells more densely linked to each other than chance would give. Resolution 0.5 sets the granularity; here it gives 14 clusters.</li></ol>
        <p>The map below is a <b>UMAP</b> of the same PCs — a 2-D picture for the eye only; clusters were found on the graph, not on this picture. ${UMAP_NOTE}</p>${umapHTML("cluster")}`,
    "b-naive": () => {
      const P = PROV, M = P.clustering.markers_shown, kept = new Set(P.clustering.clusters.filter(c => c.kept).map(c => c.cluster));
      return `<p>Clusters have numbers, not names; we name them by their <b>marker genes</b>:</p>
        <table class="cp-scores"><tr><th>gene</th><th>marks</th></tr>
        <tr><td>MS4A1 (CD20)</td><td>any B cell</td></tr><tr><td>IGHD, IGHM</td><td>an unswitched antibody receptor (IgD/IgM) — naive B cells carry both</td></tr>
        <tr><td>TCL1A, FCER2 (CD23), IL4R, BACH2</td><td>the naive-B programme</td></tr><tr><td>CD27, AIM2, TNFRSF13B (TACI), CD80</td><td>memory B cells</td></tr>
        <tr><td>IGHG1, IGHA1</td><td>class-switched antibody (memory/plasma)</td></tr><tr><td>CD3E · NKG7 · LYZ · JCHAIN</td><td>T cells · NK/cytotoxic · monocytes · plasma cells</td></tr></table>
        <p><b>The rule.</b> For each cluster, average each marker's normalised expression; score = mean(naive markers) − mean(memory markers). Keep clusters with a high score, MS4A1 present, and no T/NK/monocyte/plasma signal. The kept five score 2.29–3.69; the next, cluster 13 (1.26), is plasma-like (JCHAIN 5.7). The cut falls in a clear gap — but it was chosen by reading this table, not by a pre-set threshold.</p>
        ${umapHTML("kept")}
        <h4>Marker table</h4><div class="table-wrap"><table class="clus"><tr><th>cl.</th><th class="num">cells</th><th class="num">score</th>${M.map(m => `<th class="num">${m}</th>`).join("")}</tr>` +
        P.clustering.clusters.map(c => `<tr class="${c.kept ? "kept" : ""}"><td>${c.cluster}${c.kept ? " ✓" : ""}</td><td class="num">${n0(c.cells)}</td><td class="num">${c.naive_minus_memory}</td>${M.map(m => `<td class="num">${c[m]}</td>`).join("")}</tr>`).join("") + `</table></div>
        <h4>Why five naive clusters? Donors.</h4><p>Each kept cluster is mostly one donor — the clustering separated <em>people</em> (age, sex, sample handling) before cell states, because nothing corrected for donor. It doesn't change which cells are naive B, but it means "5 clusters" is not 5 kinds of naive B cell.</p>
        ${UMAP && UMAP.cluster_donor_pct ? `<div class="table-wrap"><table class="clus"><tr><th>cl.</th>${UMAP.donors.map(d => `<th class="num">${d.replace("healthy ", "")}</th>`).join("")}</tr>${UMAP.cluster_donor_pct.filter(r => kept.has(r.cluster)).map(r => `<tr class="kept"><td>${r.cluster}</td>${r.pct.map(v => `<td class="num">${v}%</td>`).join("")}</tr>`).join("")}</table></div>` : ""}`;
    },
    "b-splatter": () => `<p><b>Splatter</b> fits a statistical model of a scRNA-seq dataset and draws new cells from it. <code>splatEstimate</code> learns, from the naive-B counts: how gene means are spread (Gamma(${PROV.splatter.mean_shape}, ${PROV.splatter.mean_rate})), how library sizes vary (log-normal(${PROV.splatter.lib_loc}, ${PROV.splatter.lib_scale})), how noisy each gene is beyond Poisson (BCV ${PROV.splatter.bcv_common}), and how often a gene is an outlier (${PROV.splatter.out_prob}). <code>splatSimulate</code> then draws ${n0(PROV.splatter.pool_cells)} synthetic cells. They look statistically like the real ones, but their genes are anonymous — Gene1 is not any real gene — and real biology (co-expressed programmes, cell subtypes) is gone.</p>`,
    "b-pack": () => `<p>The browser can't read an R object, so the kept cells' raw counts are packed into one binary file: per cell, the number of nonzero genes, then for each nonzero gene its index (as the gap from the previous one) and its count, each as a variable-length integer (1 byte when small). ${n0(PROV.naive_b.cells)} cells × ${n0(PROV.naive_b.genes)} genes become ~12 MB; nothing is lost. Node 2 decodes it.</p>`,

    data: () => {
      const s0 = sample(), c = s0.cells.slice(0, 2);
      return `<p>Fetch the packed counts and decode them back into, per cell, a list of (gene, count) pairs. Nothing is computed yet — this is the dataset exactly as built in node 1.</p>` +
        (poolNow ? sio(`<div class="sio"><div class="sio-t">in: bytes of counts.bin</div><p class="cp-note">a varint byte stream, ${(poolNow.bytes / 1e6).toFixed(1)} MB</p></div>`,
          grid("out: raw counts (UMIs)", c.map(x => ({ ...x, label: sim.cellInfo[x.i].barcode.slice(0, 8) })), s0.genes, (x, g) => rawCount(x, g.label.replace(" (added)", "")), { fmt: v => String(v) }),
          "Rows: genes; columns: two real cells (barcodes). Most entries in a full row would be 0.") : "");
    },
    draw: () => `<p>The paper's "control" and "treatment" are not two biological groups: they are one population, split at random. ${2 * sim.n} cells are drawn without replacement from the ${poolNow ? n0(poolNow.nCells) : ""} pool; the first ${sim.n} are labelled control, the rest treatment. Because the split is random, the two groups differ only by chance — every real difference later is one the next nodes put there. That is what gives the benchmark its ground truth.</p>
        <p>Cells, not people, are the unit here. For a real experiment, the person is the replicate; for testing whether a method can find a planted difference, cells are fine.</p>` +
      sio(`<div class="sio"><div class="sio-t">in</div><p class="cp-note">${poolNow ? n0(poolNow.nCells) : "?"} cells, unlabelled</p></div>`,
        `<div class="sio"><div class="sio-t">out</div><div class="table-wrap"><table class="sio-tab"><tr><th>cell</th><th>group</th><th>donor</th><th>barcode</th></tr>${sample().cells.map(c => `<tr><td>${c.label}</td><td>${sim.group[c.i] ? "treatment" : "control"}</td><td>${esc(sim.cellInfo[c.i].donor)}</td><td><code>${esc(sim.cellInfo[c.i].barcode)}</code></td></tr>`).join("")}</table></div></div>`),
    norm: () => {
      const s0 = sample(), c = s0.cells[0], T = libSize(c), g = s0.genes.find(x => rawCount(c, x.label.replace(" (added)", "")) > 0) || s0.genes[0], k = rawCount(c, g.label.replace(" (added)", ""));
      return `<p>Cells are captured with very different efficiency: one may yield 600 molecules, another 5,000. Raw counts would mostly measure that. Seurat's LogNormalize divides each gene's count by the cell's total, scales to 10,000, and takes log(1 + x): the log tames the few very high genes, the +1 keeps zeros at zero.</p>
        <p class="formula">x<sub>gc</sub> = log(1 + count<sub>gc</sub> ÷ total<sub>c</sub> × 10,000)</p>
        <p class="cp-note">Worked example — cell ${c.label}, gene ${esc(g.label)}: count ${k}, cell total ${n0(T)} → log(1 + ${k} ÷ ${n0(T)} × 10,000) = log(1 + ${(k / T * 1e4).toFixed(2)}) = <b>${Math.log1p(k / T * 1e4).toFixed(3)}</b>.</p>` +
        sio(grid("in: raw counts", s0.cells, s0.genes, (x, gg) => rawCount(x, gg.label.replace(" (added)", "")), { fmt: v => String(v) }),
          grid("out: log-normalised", s0.cells, s0.genes, (x, gg) => { const r = rawCount(x, gg.label.replace(" (added)", "")); return Math.log1p(r / libSize(x) * 1e4); }),
          "Before any boost. C = control, T = treatment cells.");
    },
    filter: () => {
      const P = poolNow, drawn = sim.cellInfo.map(c => c.pool), N = drawn.length, nz = new Uint32Array(P.nGenes);
      for (const pc of drawn) for (let k = P.ptr[pc]; k < P.ptr[pc + 1]; k++) nz[P.gi[k]]++;
      const bins = [["0 cells (removed)", g => nz[g] === 0], ["1 cell", g => nz[g] === 1], ["2 – 1%", g => nz[g] > 1 && nz[g] / N <= 0.01],
        ["1 – 10%", g => nz[g] / N > 0.01 && nz[g] / N <= 0.1], ["10 – 25%", g => nz[g] / N > 0.1 && nz[g] / N <= 0.25], ["> 25% (S1 signal)", g => nz[g] / N > 0.25]];
      const cnt = bins.map(([, f]) => { let k = 0; for (let g = 0; g < P.nGenes; g++) if (f(g)) k++; return k; });
      const gone = []; for (let g = 0; g < P.nGenes && gone.length < 8; g++) if (nz[g] === 0) gone.push(P.genes[g]);
      const max = Math.max(...cnt);
      return `<p>The pool has ${n0(P.nGenes)} genes, but a random ${N} cells only ever see part of them. A gene that is zero in <b>every</b> drawn cell is identical in both groups and can't be scored or tested, so the authors drop it (<code>rows_zeros &lt; 1</code>). This is re-done for every draw: with more cells, more genes survive.</p>
        <p>It matters more than it looks: the kept genes are the "universe" the rank-based methods rank within (AUCell's top-5% cut-off is 5% of <em>this</em> number) and the pool non-signal genes are drawn from.</p>
        <h4>How many of the drawn ${N} cells detect each gene</h4><table class="rec">${bins.map(([l], i) => `<tr><td>${l}</td><td class="num">${n0(cnt[i])}</td><td class="bar"><span style="width:${(100 * cnt[i] / max).toFixed(0)}%;background:var(--${i === 0 ? "bad" : i === 5 ? "s1" : "muted"})"></span></td></tr>`).join("")}</table>` +
        sio(`<div class="sio"><div class="sio-t">in</div><p class="cp-note">${n0(P.nGenes)} genes × ${N} cells</p></div>`,
          `<div class="sio"><div class="sio-t">out</div><p class="cp-note">${n0(sim.G0)} genes kept; ${n0(P.nGenes - sim.G0)} removed, e.g. ${gone.map(esc).join(", ")}</p></div>`,
          "Most genes in a single-cell experiment are rare: the median kept gene is seen in under 2% of cells.");
    },
    signal: () => `<p><b>S1/S3:</b> every kept gene whose zero fraction is below ${SOURCES[cfg.source].threshold ?? 0.75} — detected in more than ${Math.round(100 * (1 - (SOURCES[cfg.source].threshold ?? 0.75)))}% of the drawn cells — becomes a "signal" gene. <b>S2/S4:</b> 550 kept genes at random. In S1/S2 these are the genes that get changed; in S3/S4 nothing is changed and they are just the pool gene sets are drawn from.</p>
        <p>In these B cells the S1 genes are mostly ribosomal (RPL/RPS), mitochondrial (MT-) and other housekeeping genes, plus a few B-cell genes (CD74, CD52): what every cell expresses. So S1 asks "can a method see a small shift in very commonly detected genes", not "can it see a pathway". S2's random genes are mostly rarely detected, which is why S2 is harder.</p>
        <p class="cp-note">${sim.signal.length} genes this run; their median detection ${S.dataSummary(sim).medianPctCellsExpressingSignal.toFixed(1)}%. First 30: ${sim.signal.slice(0, 30).map(g => esc(sim.geneNames[g])).join(", ")}…</p>`,
    boost: () => {
      const s0 = sample(), b = sim.paper ? sim.boost : sim.opt.foldChange;
      const ex = [0.5, 1, 2, 3].map(x => [x, Math.expm1(x), Math.expm1(x * 1.2)]);
      return `<p>The "signal" is a multiplication of the <b>log-normalised</b> value of each signal gene, in treatment cells only: x → ${b} × x. Control cells and all other genes are untouched.</p>
        <p>×1.2 on the log scale is <em>not</em> +20% expression. Undo the log: x = log(1 + y), so 1.2x means y' = (1 + y)<sup>1.2</sup> − 1. The change grows with the level, and zeros stay zero (a gene not detected in a cell is never "boosted"):</p>
        <table class="cp-scores"><tr><th class="num">log value x</th><th class="num">expression y</th><th class="num">after ×1.2</th><th class="num">change in y</th></tr>${ex.map(([x, y, y2]) => `<tr><td class="num">${x}</td><td class="num">${y.toFixed(2)}</td><td class="num">${y2.toFixed(2)}</td><td class="num">+${(100 * (y2 / y - 1)).toFixed(0)}%</td></tr>`).join("")}</table>` +
        (sim.paper ? sio(grid("in: log-normalised", s0.cells, s0.genes, (c, g) => isBoosted(c, g) ? sim.data[c.i * sim.G + g.g] / b : sim.data[c.i * sim.G + g.g]),
          grid("out", s0.cells, s0.genes, (c, g) => sim.data[c.i * sim.G + g.g], { mark: isBoosted }),
          sim.boost === 1 ? "No boost in this scenario: in = out." : "Highlighted: changed values (signal genes × treatment cells).") : "");
    },
    added: () => {
      if (cfg.scenario < 3) return `<p>Only in S3/S4. Switch the scenario in node "Choose the signal genes".</p>`;
      const s0 = sample(), cs = sim.csGenes.slice(0, 3).map(g => ({ g, label: sim.geneNames[g], tag: "added", tagcls: "added" }));
      return `<p>S3/S4 model two groups that differ in how many genes they express (e.g. cancer vs normal cells), without changing any shared gene. The authors draw a second, independent set of ${2 * sim.n} cells, pick ${sim.nCS} genes from it (S3: genes detected in >25% of those cells; S4: any detected gene), and append them as <b>new</b> genes: in treatment cell <i>i</i>, the value is that gene's value in the second draw's cell <i>i</i>; in every control cell it is 0.</p>
        <p>Result: treatment cells detect more genes (${S.dataSummary(sim).genesPerCell.map(v => v.toFixed(0)).join(" → ")} per cell). A gene set with none of the added genes is truly unchanged — yet rank- and total-based scores shift, because the added genes change each cell's ranks and totals. That's the trap S3 exposes.</p>` +
        sio(`<div class="sio"><div class="sio-t">in</div><p class="cp-note">${n0(sim.G0)} genes</p></div>`,
          grid(`out: ${n0(sim.G)} genes (added rows shown)`, s0.cells, cs, (c, g) => sim.data[c.i * sim.G + g.g], { mark: (c, g) => sim.group[c.i] === 1 }));
    },
    geneset: () => {
      const csSet = new Set(sim.csGenes);
      return `<p>A "gene set" here is just a list of genes to score. The authors build ${cfg.size}-gene sets as ${cfg.size} × (1 − noise) genes drawn from the signal pool plus ${cfg.size} × noise from all other kept genes (in S3/S4 that "other" pool includes the added genes, so a set may pick some up by chance). Noise 0 = the set is pure signal; noise 1 = it contains no signal at all, so any detection is a false positive.</p>
        <p class="cp-note">This set (${set.length} genes): ${set.map(g => `${esc(sim.geneNames[g])}${sim.isSignal[g] ? "*" : csSet.has(g) ? "†" : ""}`).join(", ")}<br>* signal pool · † added (treatment-only)</p>`;
    },
    prepare: () => {
      const s0 = sample(), c = s0.cells[0], G = sim.G;
      const vals = Array.from(sim.data.subarray(c.i * G, (c.i + 1) * G));
      const order = vals.map((v, g) => [v, g]).sort((a, b) => b[0] - a[0]).slice(0, 6);
      return `<p>Rank-based methods don't use the values themselves but each gene's <b>rank within its cell</b> (1 = the cell's highest). Ranks are independent of the gene set, so they are computed once for all ${n0(sim.N)} cells × ${n0(G)} genes: AUCell (ties broken at random), UCell (ties averaged), ssGSEA (its own rank transform), JASMINE (ranks among the cell's detected genes only). AddModuleScore instead bins all genes by average expression, to pick matched control genes later.</p>
        <p>Zeros are the catch: a cell with 600 detected genes out of ${n0(G)} has ${n0(G - 600)} genes tied at the bottom. How a method treats that tie is much of why the methods disagree.</p>` +
        sio(grid(`in: cell ${c.label}, top values`, [c], order.map(([, g]) => ({ g, label: sim.geneNames[g] })), (x, g) => sim.data[x.i * G + g.g]),
          grid("out: rank in the cell (UCell)", [c], order.map(([, g]) => ({ g, label: sim.geneNames[g] })), (x, g) => prep.rU[x.i * G + g.g], { fmt: v => String(Math.round(v * 10) / 10) }));
    },
    score: () => `<p>Each method turns "this cell's values for the set's genes" into one number. In one line each (full first-principles pages under <em>Methods</em> in the sidebar):</p>
        <table class="cp-scores"><tr><td><a href="methods/addmodulescore.html">AddModuleScore</a></td><td>mean of the set's genes minus mean of matched random control genes</td></tr>
        <tr><td><a href="methods/aucell.html">AUCell</a></td><td>how many of the set's genes sit in the cell's top 5% of ranks (area under a recovery curve)</td></tr>
        <tr><td><a href="methods/ucell.html">UCell</a></td><td>a Mann–Whitney U on the set's ranks, capped at rank 1,500</td></tr>
        <tr><td><a href="methods/ssgsea.html">ssGSEA</a></td><td>a running-sum walk down the cell's ranked genes</td></tr>
        <tr><td><a href="methods/jasmine.html">JASMINE</a></td><td>rank among detected genes, plus an odds ratio of detection in vs out of the set</td></tr>
        <tr><td><a href="methods/scse.html">SCSE</a></td><td>the set's share of the cell's total expression</td></tr>
        <tr><td><a href="methods/scps.html">scPS</a></td><td>PCA on the set's genes alone, weighted by the set's mean expression</td></tr></table>` +
      (res ? sio(grid("in: this set's genes", sample().cells, sample().genes.filter(g => set.includes(g.g)).slice(0, 4), (c, g) => sim.data[c.i * sim.G + g.g]),
        `<div class="sio"><div class="sio-t">out: one score per cell</div><div class="table-wrap"><table class="sio-tab"><tr><th></th>${sample().cells.map(c => `<th class="num">${c.label}</th>`).join("")}</tr>${METHODS.map(m => `<tr><td>${m}</td>${sample().cells.map(c => `<td class="num">${Number(res[m].scores[c.i]).toPrecision(3)}</td>`).join("")}</tr>`).join("")}</table></div></div>`,
        "Only the set genes these four cells detect are shown on the left; each score uses all the set's genes (and, for some methods, the rest of the cell).") : ""),
    test: () => `<p>For each method: are treatment cells' scores systematically higher or lower than control cells'? The Wilcoxon rank-sum test pools all ${n0(sim.N)} scores, ranks them, and asks whether the treatment ranks sum to more (or less) than chance would give. Two-sided, with a continuity correction, as the authors' <code>row_wilcoxon_twosample</code>. It sees only the final scores, not why they differ — so an artefact that shifts scores is "detected" just like real signal. <a href="wilcoxon.html">Play with it here</a>.</p>` +
      (res ? `<table class="cp-scores"><tr><th>method</th><th class="num">p</th><th class="num">P(treatment &gt; control)</th></tr>${METHODS.map(m => `<tr><td>${m}</td><td class="num">${res[m].test.p < 1e-4 ? res[m].test.p.toExponential(1) : res[m].test.p.toFixed(4)}</td><td class="num">${(100 * res[m].test.auc).toFixed(0)}%</td></tr>`).join("")}</table>` : ""),
    recovery: () => `<p>One gene set says little. The paper scores 100 random sets per setting, tests each, and — because testing 100 sets would give ~5 false hits at p &lt; 0.05 by chance — adjusts the p-values with <b>Benjamini–Hochberg</b>, which controls the expected share of false discoveries. A set counts as "identified" if its adjusted p &lt; 0.05. The recovery rate is the share identified: on sets with real signal that's sensitivity; on sets with none, it's the false-positive rate. They repeat this for 10 replicate draws; we show one.</p>`,
    sweep: () => `<p>Same as the recovery node, at several gene set sizes. The paper's key observation: with a systematic artefact (S3), bigger sets make false positives <em>more</em> likely, because averaging more genes shrinks random noise but not a shift shared by every gene.</p>`,
    toy: () => `<p>Made-up data: each gene gets a random mean, and each cell's count is drawn from a negative binomial around it (a Poisson whose rate itself varies — the usual model for single-cell counts). The scenario's change is applied inside. Useful for seeing mechanisms on a small, fast dataset; not the paper's data.</p>`,
  };
  const UMAP_NOTE = "Colour by cluster, by kept/dropped, by donor, or by any marker gene.";
  function umapHTML(def) {
    const genes = ["MS4A1", "IGHD", "TCL1A", "FCER2", "CD27", "IGHG1", "CD3E", "NKG7", "LYZ", "JCHAIN"];
    return `<div id="umapBox" class="umap" data-def="${def}"><label class="cp-note">Colour by <select id="umapColor">${[["cluster", "cluster"], ["kept", "kept as naive B"], ["donor", "donor"], ...genes.map(g => [g, g])].map(([v, t]) => `<option value="${v}"${v === def ? " selected" : ""}>${t}</option>`).join("")}</select></label><canvas></canvas><div class="legend2"></div><p class="cp-note">5,000 of the ${PROV ? n0(PROV.qc.cells_after) : ""} QC-passed cells, random sample.</p></div>`;
  }

  async function explain(key) {
    if (key === "b-cluster" || key === "b-naive") { if (!UMAP) UMAP = await getJSON("data/provenance_umap.json"); }
    const f = EXPLAIN[key];
    if (!f) return;
    const needsSim = !key.startsWith("b-") && key !== "toy" && key !== "recovery" && key !== "sweep";
    if (needsSim && !sim) return;
    let body;
    try { body = f(); } catch (e) { body = `<p class="cp-note">(Couldn't build this explanation: ${esc(e.message)})</p>`; }
    const title = key.startsWith("b-") ? { "b-geo": "Raw data from GEO", "b-healthy": "Healthy donors only", "b-qc": "Quality control", "b-cluster": "Clustering", "b-naive": "Which clusters are naive B cells", "b-splatter": "Splatter (SSD)", "b-pack": "Packing for the browser" }[key] : TITLES[key];
    const panel = $("sidePanel");
    panel.innerHTML = `<button class="cp-close" type="button" aria-label="Close">×</button><div class="kicker">${key.startsWith("b-") ? "Build the dataset" : "Node"}</div><h3>${esc(title)}</h3>${body}` +
      (CODE[key] ? `<p class="cp-note"><button type="button" class="codechip" data-codepanel="${key}">show the code</button></p>` : "");
    panel.hidden = false;
    panel.scrollTop = 0;
    panel.querySelector(".cp-close").addEventListener("click", closePanel);
    const cp = panel.querySelector("[data-codepanel]");
    if (cp) cp.addEventListener("click", () => showCode(cp.dataset.codepanel));
    panel.querySelectorAll("[data-explain]").forEach(b => b.addEventListener("click", () => explain(b.dataset.explain)));
    const sel = panel.querySelector("#umapColor");
    if (sel) { sel.addEventListener("change", () => umapPanel(sel.value)); umapPanel(sel.value); }
  }


  // ------------------------------------------------------------ events
  function onChange(e) {
    const k = e.target.dataset.k;
    if (!k) return;
    const v = e.target.value;
    cfg[k] = k === "source" || k === "inject" ? v : Number(v);
    if (k === "source" || k === "scenario") { render(); describeBuild(); }
    if (k === "boostNull" || k === "fold") { render(); describeBuild(); }
    run(SIM_KEYS.has(k) ? "sim" : "set");
  }

  function init() {
    applyHash();
    render();
    loadSide().then(describeBuild);
    const g = $("graph");
    g.addEventListener("change", onChange);
    g.addEventListener("click", e => {
      const code = e.target.closest("[data-code]");
      if (code) { showCode(code.dataset.code); return; }
      const ex = e.target.closest("[data-explain]");
      if (ex) { explain(ex.dataset.explain); return; }
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
    document.addEventListener("keydown", e => {
      if (e.key === "Escape") closePanel();
      if (e.key === "Enter" && e.target.matches && e.target.matches("li[data-explain]")) explain(e.target.dataset.explain);
    });
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
