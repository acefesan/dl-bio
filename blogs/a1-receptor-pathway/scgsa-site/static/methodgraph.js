// "Step by step on real cells": each method page gets the method's own
// computation drawn as a chain of nodes, run on the paper's real naive-B cells
// (the playground's RWSD data). One cell is followed through every step with
// its real intermediate values; the last step checks the result equals the
// score the playground computes (static/scgsa.js) for that cell.
// Usage: <div id="mg" data-method="AUCell"></div> plus scgsa.js and this file.
(function () {
  const S = window.SCGSA;
  const root = document.documentElement.dataset.root || "";
  const esc = t => String(t).replace(/[&<>"]/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));
  const n0 = v => Number(v).toLocaleString();
  const f = (v, d = 3) => (v === 0 ? "0" : Math.abs(v) >= 1000 || Math.abs(v) < 0.001 ? Number(v).toPrecision(3) : Number(v).toFixed(d));
  const OURS = "https://github.com/acefesan/dl-bio/blob/main/blogs/a1-receptor-pathway/scgsa-site/static/scgsa.js";

  const SETTINGS = {
    s1sig: { label: "S1 · set made of boosted genes (a real difference)", scenario: 1, noise: 0 },
    s1null: { label: "S1 · set with no boosted gene (no real difference)", scenario: 1, noise: 1 },
    s3null: { label: "S3 · no change, but treatment cells gained 250 genes", scenario: 3, noise: 0 },
  };
  const cfg = { setting: "s1sig", size: 25, cell: 200, setSeed: 1 };
  let pool = null, sim = null, prep = null, set = null, res = null, simKey = "";

  // ------------------------------------------------------------ code snippets
  let src = null;
  async function snippet(find, lines) {
    if (!src) {
      const tag = document.querySelector('script[src*="static/scgsa.js"]');
      src = await fetch(tag ? tag.src : root + "static/scgsa.js", { cache: "no-cache" }).then(r => r.text());
    }
    const L = src.split("\n"), i = L.findIndex(l => l.includes(find));
    if (i < 0) return null;
    let end = i + (lines || 1) - 1;
    if (!lines) {
      let d = 0, seen = false;
      for (let j = i; j < L.length; j++) { for (const ch of L[j]) { if (ch === "{") { d++; seen = true; } else if (ch === "}") d--; } if (seen && d <= 0) { end = j; break; } }
    }
    return { start: i + 1, end: end + 1, text: L.slice(i, end + 1).join("\n") };
  }

  // same algorithm as scgsa.js sampleWithout (not exported), so control genes match exactly
  function sampleWithout(rng, poolArr, k) {
    const a = Array.from(poolArr);
    k = Math.min(k, a.length);
    for (let i = 0; i < k; i++) { const j = i + rng.int(a.length - i); const t = a[i]; a[i] = a[j]; a[j] = t; }
    return a.slice(0, k);
  }

  // ------------------------------------------------------------ small renderers
  const kv = items => `<div class="kvs">${items.map(([b, s]) => `<div class="kv"><b>${b}</b><span>${s}</span></div>`).join("")}</div>`;
  const table = (head, rows, num) => `<div class="table-wrap"><table class="mg-tab"><tr>${head.map((h, i) => `<th class="${num && num[i] ? "num" : ""}">${h}</th>`).join("")}</tr>` +
    rows.map(r => `<tr>${r.map((c, i) => `<td class="${num && num[i] ? "num" : ""}">${c}</td>`).join("")}</tr>`).join("") + `</table></div>`;
  const io = (a, b) => `<div class="sio-row"><div class="sio"><div class="sio-t">in</div>${a}</div><div class="sio-arrow">→</div><div class="sio"><div class="sio-t">out</div>${b}</div></div>`;
  const role = g => sim.csGenes.includes(g) ? `<span class="role added">added</span>` : sim.isSignal[g] ? `<span class="role signal" title="${sim.boost === 1 ? "drawn from the dense-gene pool; not changed in this scenario" : "a signal gene: ×1.2 in treatment cells only"}">${sim.boost === 1 ? "pool" : "signal"}</span>` : "";
  const gname = g => `${esc(sim.geneNames[g])} ${role(g)}`;
  const cellLabel = c => c < sim.n ? `control cell C${c + 1}` : `treatment cell T${c - sim.n + 1}`;
  const val = (c, g) => sim.data[c * sim.G + g];
  function verify(m, mine) {
    const ref = res[m].scores[cfg.cell];
    const ok = Math.abs(mine - ref) <= 1e-9 * Math.max(1, Math.abs(ref));
    return `<p class="mg-check ${ok ? "ok" : "bad"}">${ok ? "✓" : "✗"} Same as the score the playground computes for this cell (${f(ref, 4)}).</p>`;
  }
  function curve(points, W, H, xs, ys, opts) {
    const X = x => 30 + (x - xs[0]) / (xs[1] - xs[0]) * (W - 40), Y = y => H - 18 - (y - ys[0]) / (ys[1] - ys[0] || 1) * (H - 28);
    let s = `<svg class="mg-plot" viewBox="0 0 ${W} ${H}"><line class="ax" x1="30" x2="${W - 10}" y1="${Y(ys[0])}" y2="${Y(ys[0])}"/><line class="ax" x1="30" x2="30" y1="10" y2="${H - 18}"/>`;
    if (opts.zero) s += `<line class="zero" x1="30" x2="${W - 10}" y1="${Y(0)}" y2="${Y(0)}"/>`;
    for (const [pts, cls] of points) s += `<polyline class="${cls}" points="${pts.map(([x, y]) => `${X(x).toFixed(1)},${Y(y).toFixed(1)}`).join(" ")}"/>`;
    for (const [x, label] of opts.vlines || []) s += `<line class="cut" x1="${X(x)}" x2="${X(x)}" y1="10" y2="${H - 18}"/><text x="${X(x) + 3}" y="20" class="cutl">${label}</text>`;
    s += `<text x="${W - 10}" y="${H - 4}" text-anchor="end">${opts.xlab}</text><text x="32" y="10">${opts.ylab}</text>`;
    return s + `</svg>`;
  }

  // ------------------------------------------------------------ the seven methods, step by step
  const STEPS = {
    AddModuleScore: () => {
      const c = cfg.cell, k = set.length, G = sim.G;
      const rng = S.makeRng(1), ctrl = new Set();
      for (const g of set) for (const x of sampleWithout(rng, prep.bins[prep.bin[g]], Math.min(100, prep.bins[prep.bin[g]].length))) ctrl.add(x);
      const ctl = Array.from(ctrl);
      let a = 0, b = 0; for (const g of set) a += val(c, g); for (const g of ctl) b += val(c, g);
      const sm = a / k, cm = b / ctl.length, score = sm - cm;
      const show = set.slice(0, 6);
      return [
        { t: "Take the gene set's values in this cell", p: "Start from the cell's log-normalised value for each gene in the set.",
          b: io(`<p class="cp-note">${k} genes in the set</p>`, table(["gene", "value in this cell"], show.map(g => [gname(g), f(val(c, g), 2)]), [0, 1]) + `<p class="cp-note">first 6 of ${k} shown</p>`), code: ["let a = 0, b = 0;", 4] },
        { t: "Sort every gene into 24 bins by its average level", p: `Average each of the ${n0(G)} genes over all ${sim.N} cells, sort them, and cut the list into 24 equal bins (~${Math.round(G / 24)} genes each). Bin 1 = the rarest genes, bin 24 = the most abundant. The bins let the method compare like with like.`,
          b: table(["gene", "average over all cells", "bin (1–24)"], show.map(g => [gname(g), f(prep.avg[g], 3), prep.bin[g] + 1]), [0, 1, 1]), code: ["const avg = new Float64Array(G);", 13] },
        { t: "Pick control genes from the same bins", p: "For each set gene, draw 100 random genes from its bin. Together they form a 'control' set with the same overall expression levels as the gene set but no shared biology. It answers: how high would a random set of similarly abundant genes score in this cell?",
          b: kv([[n0(ctl.length), `distinct control genes (${k} genes × up to 100, overlaps merged)`]]) + `<p class="cp-note">e.g. ${ctl.slice(0, 10).map(g => esc(sim.geneNames[g])).join(", ")}, …</p>`, code: ["const ctrlSet = new Set();", 6] },
        { t: "Average both sets in this cell", p: "Mean value of the gene-set genes in this cell, and mean value of the control genes in this cell.",
          b: kv([[f(sm, 4), `gene-set mean (${k} genes)`], [f(cm, 4), `control mean (${n0(ctl.length)} genes)`]]), code: ["out[c] = a / set.length - b / ctl.length;", 1] },
        { t: "Score = gene-set mean − control mean", p: "Above 0: the set is higher in this cell than comparable random genes. Below 0: lower.",
          b: kv([[f(score, 4), `${f(sm, 4)} − ${f(cm, 4)}`]]) + verify("AddModuleScore", score), code: ["function addModuleScore", 0] },
      ];
    },
    AUCell: () => {
      const c = cfg.cell, G = sim.G, k = set.length, thr = Math.round(prep.aucMaxRank);
      const ranks = set.map(g => [g, prep.rAUC[c * G + g]]).sort((x, y) => x[1] - y[1]);
      const xs = ranks.map(r => r[1]).filter(r => r < thr);
      let s = 0; for (let i = 0; i < xs.length; i++) s += ((i + 1 < xs.length ? xs[i + 1] : thr) - xs[i]) * (i + 1);
      let maxA = 0; { const m = []; for (let i = 1; i <= k; i++) if (i < thr) m.push(i); for (let i = 0; i < m.length; i++) maxA += ((i + 1 < m.length ? m[i + 1] : thr) - m[i]) * (i + 1); }
      const score = maxA > 0 ? s / maxA : 0;
      const det = ranks.filter(([g]) => val(c, g) > 0).length;
      const steps = [[0, 0]]; xs.forEach((x, i) => { steps.push([x, i]); steps.push([x, i + 1]); }); steps.push([thr, xs.length]);
      const ideal = [[0, 0]]; for (let i = 1; i <= Math.min(k, thr - 1); i++) { ideal.push([i, i - 1]); ideal.push([i, i]); } ideal.push([thr, Math.min(k, thr - 1)]);
      return [
        { t: "Rank all genes inside this cell", p: `Sort the cell's ${n0(G)} genes from highest to lowest value; rank 1 = the cell's most expressed gene. Genes with value 0 are all tied at the bottom, so AUCell orders them at random.`,
          b: io(`<p class="cp-note">${k} gene-set genes; ${det} of them detected in this cell</p>`, table(["gene", "value", "rank in this cell"], ranks.slice(0, 8).map(([g, r]) => [gname(g), f(val(c, g), 2), n0(r)]), [0, 1, 1]) + `<p class="cp-note">the set's 8 best-ranked genes</p>`), code: ['const ra = rankDesc(v, "random", rng);', 1] },
        { t: "Look only at the top 5% of the ranking", p: `The cut-off is 5% of ${n0(G)} genes = rank ${n0(thr)}. Anything ranked lower is ignored, however close it is.`,
          b: kv([[n0(thr), "cut-off rank"], [`${xs.length} / ${k}`, "set genes inside the top 5%"]]), code: ["const aucMaxRank = Math.ceil(0.05 * G);", 1] },
        { t: "Draw the recovery curve", p: "Walk down the ranking from the top; each time you meet a set gene, the curve steps up by one. A set concentrated at the top climbs early and encloses a big area.",
          b: curve([[ideal, "ideal"], [steps, "rec"]], 460, 170, [0, thr], [0, Math.max(1, Math.min(k, thr - 1))], { xlab: "rank in this cell (top 5%)", ylab: "set genes found" }) + `<p class="cp-note">Solid: this cell. Dashed: the best possible (all set genes at the very top).</p>`, code: ["for (const g of set) { const r = prep.rAUC", 7] },
        { t: "Score = area under the curve ÷ best possible area", p: "0 = no set gene in the top 5%; 1 = the set's genes are the cell's very top genes.",
          b: kv([[f(score, 4), `${n0(s)} ÷ ${n0(maxA)}`]]) + verify("AUCell", score), code: ["function aucell", 0] },
      ];
    },
    UCell: () => {
      const c = cfg.cell, G = sim.G, k = set.length, maxR = prep.ucellMaxRank;
      const ranks = set.map(g => [g, prep.rU[c * G + g]]).sort((x, y) => x[1] - y[1]);
      const capped = ranks.map(([g, r]) => [g, Math.min(r, maxR)]);
      const sum = capped.reduce((a, [, r]) => a + r, 0), minSum = k * (k + 1) / 2;
      const anySig = ranks.some(([, r]) => r < maxR);
      const score = anySig ? 1 - (sum - minSum) / (k * maxR - minSum) : 0;
      const zeroRank = ranks.find(([g]) => val(c, g) === 0);
      return [
        { t: "Rank all genes inside this cell (ties averaged)", p: `Rank 1 = the cell's top gene. Tied genes share the average of their ranks, so every undetected gene gets the same, very low rank${zeroRank ? ` (here ${n0(zeroRank[1])})` : ""}.`,
          b: table(["gene", "value", "rank"], ranks.slice(0, 8).map(([g, r]) => [gname(g), f(val(c, g), 2), n0(r)]), [0, 1, 1]), code: ['const ru = rankDesc(v, "average");', 1] },
        { t: `Cap every rank at ${n0(maxR)}`, p: `A gene ranked below ${n0(maxR)} is treated as rank ${n0(maxR)}: "not near the top" is all that counts. (The authors' tool, irGSEA, uses 1,500.)`,
          b: kv([[`${capped.filter(([, r]) => r === maxR).length} / ${k}`, `set genes capped at ${n0(maxR)}`], [n0(Math.round(sum)), "sum of the capped ranks"]]), code: ["if (r >= maxRank) r = maxRank; else anySig = true;", 1] },
        { t: "Turn the rank sum into a 0–1 score", p: `This is a Mann–Whitney U statistic: subtract the smallest possible sum (${n0(minSum)}, all set genes at the very top), divide by the full range, and flip so 1 = best.`,
          b: kv([[f(score, 4), `1 − (${n0(Math.round(sum))} − ${n0(minSum)}) ÷ (${k} × ${n0(maxR)} − ${n0(minSum)})`]]) + verify("UCell", score), code: ["out[c] = anySig ?", 1] },
      ];
    },
    ssGSEA: () => {
      const c = cfg.cell, G = sim.G, n = G, k = set.length;
      const inSet = new Set(set);
      let wsum = 0, wpos = 0, psum = 0;
      for (const g of set) { const pos = prep.ssPos[c * G + g], w = prep.ssRa[c * G + g]; wsum += w; wpos += w * (n - pos + 1); psum += n - pos + 1; }
      const inStep = wsum > 0 ? wpos / wsum : 0, outStep = (n * (n + 1) / 2 - psum) / (n - k), score = inStep - outStep;
      // running walk, for the picture: by position (1 = top)
      const byPos = new Array(n); for (let g = 0; g < G; g++) byPos[prep.ssPos[c * G + g] - 1] = g;
      let cin = 0, cout = 0; const pts = []; const step = Math.max(1, Math.floor(n / 600));
      for (let p = 0; p < n; p++) { const g = byPos[p]; if (inSet.has(g)) cin += prep.ssRa[c * G + g] / wsum; else cout += 1 / (n - k); if (p % step === 0 || p === n - 1) pts.push([p + 1, cin - cout]); }
      const mn = Math.min(...pts.map(q => q[1])), mx = Math.max(...pts.map(q => q[1]));
      const ranked = set.map(g => [g, prep.ssPos[c * G + g]]).sort((a, b) => a[1] - b[1]);
      return [
        { t: "Rank the cell's genes and weight each one", p: "Every gene gets a position (1 = the cell's highest) and a weight that grows with its rank (rank to the power ¼), so highly expressed set genes count more.",
          b: table(["gene", "value", "position", "weight"], ranked.slice(0, 8).map(([g, p]) => [gname(g), f(val(c, g), 2), n0(p), f(prep.ssRa[c * G + g], 2)]), [0, 1, 1, 1]), code: ["const Rasc = rankAscAvg(v);", 7] },
        { t: "Walk down the ranking", p: "Start at 0 at the top gene. At a set gene, step up by its share of the set's total weight; at any other gene, step down a little. If the set sits near the top, the walk climbs early and stays high.",
          b: curve([[pts, "walk"]], 460, 170, [1, n], [Math.min(0, mn), Math.max(0, mx)], { xlab: `position in this cell (1 … ${n0(n)})`, ylab: "running sum", zero: true }), code: ["let wsum = 0, wpos = 0, psum = 0;", 8] },
        { t: "Score = the total area under the walk", p: "Unlike classic GSEA (which takes the walk's peak), ssGSEA adds up the walk over every position. Positive = the set sits above the cell's typical gene.",
          b: kv([[f(score, 4), `${f(inStep, 1)} − ${f(outStep, 1)} (set side − rest side)`]]) + verify("ssGSEA", score), code: ["function ssgsea", 0] },
      ];
    },
    JASMINE: () => {
      const c = cfg.cell, G = sim.G, k = set.length, N = sim.N;
      const RMall = new Float64Array(N), ORall = new Float64Array(N);
      let mine = null;
      for (let cc = 0; cc < N; cc++) {
        let rs = 0, m = 0; for (const g of set) { const r = prep.rNZ[cc * G + g]; if (r > 0) { rs += r; m++; } }
        const RM = m > 0 ? (rs / m) / prep.nnz[cc] : 0;
        let SE = m, NE = prep.nnz[cc] - SE, SN = k - SE; if (SN === 0) SN = 1; if (NE === 0) NE = 1;
        let NN = G - (NE + SE); NN -= SN; const OR = (SE * NN) / (SN * NE);
        RMall[cc] = RM; ORall[cc] = OR;
        if (cc === c) mine = { rs, m, RM, SE, NE, SN, NN, OR };
      }
      const mm = a => { let lo = Infinity, hi = -Infinity; for (const v of a) { if (v < lo) lo = v; if (v > hi) hi = v; } return [lo, hi]; };
      const [rl, rh] = mm(RMall), [ol, oh] = mm(ORall);
      const rS = rh > rl ? (mine.RM - rl) / (rh - rl) : 0, oS = oh > ol ? (mine.OR - ol) / (oh - ol) : 0, score = (rS + oS) / 2;
      const det = set.filter(g => val(c, g) > 0);
      return [
        { t: "Keep only the genes this cell detected", p: `JASMINE ignores zeros entirely: this cell detected ${n0(prep.nnz[c])} genes; ${det.length} of the ${k} set genes are among them.`,
          b: kv([[n0(prep.nnz[c]), "genes detected in this cell"], [`${det.length} / ${k}`, "set genes detected"]]), code: ["const nzIdx = [];", 9] },
        { t: "Rank the detected genes; average the set's ranks", p: "Rank the detected genes from lowest (1) to highest. Average the set genes' ranks and divide by the number of detected genes: near 1 = the set's genes are among this cell's highest.",
          b: table(["gene", "value", "rank among detected"], det.slice(0, 6).map(g => [gname(g), f(val(c, g), 2), n0(prep.rNZ[c * G + g])]), [0, 1, 1]) + kv([[f(mine.RM, 4), `mean rank ${f(mine.rs / Math.max(1, mine.m), 1)} ÷ ${n0(prep.nnz[c])} detected`]]), code: ["RM[c] = m > 0", 1] },
        { t: "Is the set enriched among detected genes? (odds ratio)", p: "A 2×2 table: set vs not-set genes × detected vs not detected in this cell. The odds ratio says how much more likely a set gene is to be detected than other genes.",
          b: table(["", "detected", "not detected"], [["in the set", n0(mine.SE), n0(mine.SN)], ["not in the set", n0(mine.NE), n0(mine.NN)]], [0, 1, 1]) + kv([[f(mine.OR, 3), "odds ratio = (SE × NN) ÷ (SN × NE)"]]), code: ["let SE = m;", 8] },
        { t: "Rescale both across all cells and average", p: `Each part is min–max scaled over all ${N} cells (0 = lowest cell, 1 = highest) and the two are averaged. So a cell's JASMINE score depends on the other cells in the dataset.`,
          b: kv([[f(rS, 3), `rank part, scaled (${f(mine.RM, 3)} in ${f(rl, 3)}–${f(rh, 3)})`], [f(oS, 3), `odds part, scaled (${f(mine.OR, 2)} in ${f(ol, 2)}–${f(oh, 2)})`], [f(score, 4), "score = their average"]]) + verify("JASMINE", score), code: ["const a = minmax(RM), b = minmax(OR);", 2] },
      ];
    },
    SCSE: () => {
      const c = cfg.cell; let s = 0; for (const g of set) s += val(c, g);
      const tot = prep.colSum[c], score = tot > 0 ? 100 * s / tot : 0;
      return [
        { t: "Add up the set's values in this cell", p: "Sum the log-normalised values of the set's genes in this cell.",
          b: table(["gene", "value"], set.slice(0, 6).map(g => [gname(g), f(val(c, g), 2)]), [0, 1]) + kv([[f(s, 3), `sum over all ${set.length} set genes`]]), code: ["for (const g of set) s += data[c * G + g];", 1] },
        { t: "Add up everything in this cell", p: `The same sum over all ${n0(sim.G)} genes: the cell's total.`,
          b: kv([[f(tot, 1), "sum over all genes in this cell"]]), code: ["colSum[c] = s;", 1] },
        { t: "Score = the set's share of the total, in %", p: "What percentage of this cell's (log) expression goes to the set. Extra genes elsewhere in the cell raise the total and lower the share — the S3 weakness.",
          b: kv([[f(score, 4), `100 × ${f(s, 2)} ÷ ${f(tot, 1)}`]]) + verify("SCSE", score), code: ["function scse", 0] },
      ];
    },
    scPS: () => {
      const c = cfg.cell, G = sim.G, N = sim.N, k = set.length;
      const A = new Float64Array(N * k), mean = [], sd = [];
      for (let j = 0; j < k; j++) {
        const g = set[j]; let s = 0; for (let i = 0; i < N; i++) s += val(i, g); const m = s / N;
        let vv = 0; for (let i = 0; i < N; i++) { const d = val(i, g) - m; vv += d * d; } const sdv = Math.sqrt(vv / (N - 1));
        mean.push(m); sd.push(sdv);
        for (let i = 0; i < N; i++) { let z = sdv > 0 ? (val(i, g) - m) / sdv : 0; if (z > 10) z = 10; A[i * k + j] = z; }
      }
      const C = new Float64Array(k * k);
      for (let i = 0; i < N; i++) for (let a = 0; a < k; a++) { const x = A[i * k + a]; if (x === 0) continue; for (let b = a; b < k; b++) C[a * k + b] += x * A[i * k + b]; }
      for (let a = 0; a < k; a++) for (let b = 0; b < a; b++) C[a * k + b] = C[b * k + a];
      const { values, vectors } = S.jacobiEigen(C, k);
      const order = Array.from({ length: k }, (_, i) => i).sort((a, b) => values[b] - values[a]);
      const npcs = Math.max(1, Math.min(10, k - 1, N - 1));
      const emb = [], sdev = [];
      for (let t = 0; t < npcs; t++) {
        const e = order[t], d = Math.sqrt(Math.max(values[e], 0));
        let big = 0; for (let i = 0; i < k; i++) { const x = vectors[i * k + e]; if (Math.abs(x) > Math.abs(big)) big = x; }
        const sign = big < 0 ? -1 : 1, u = new Float64Array(N);
        for (let i = 0; i < N; i++) { let s = 0; for (let a = 0; a < k; a++) s += A[i * k + a] * vectors[a * k + e]; u[i] = d > 0 ? sign * s / d : 0; }
        emb.push(u); sdev.push(d / Math.sqrt(Math.max(1, N - 1)));
      }
      const tot = sdev.reduce((s, x) => s + x * x, 0), varExp = sdev.map(x => (tot > 0 ? x * x / tot : 0));
      let cum = 0, maxComb = npcs; for (let t = 0; t < npcs; t++) { cum += varExp[t]; if (cum > 0.5) { maxComb = t + 1; break; } }
      let gmin = Infinity; for (const u of emb) for (const x of u) if (x < gmin) gmin = x;
      let agg = 0; for (let t = 0; t < maxComb; t++) agg += (emb[t][c] - gmin) * varExp[t]; agg = Math.sqrt(agg);
      let ms = 0; for (const g of set) ms += val(c, g); const mexp = ms / k, score = agg * mexp;
      return [
        { t: "Standardise each set gene across all cells", p: `Take the ${N} cells × ${k} set-genes matrix and rescale each gene to average 0 and spread 1 (capped at 10), so every gene counts equally. Below: this cell's values before and after.`,
          b: io(table(["gene", "value"], set.slice(0, 5).map(g => [gname(g), f(val(c, g), 2)]), [0, 1]),
            table(["gene", "standardised"], set.slice(0, 5).map((g, j) => [esc(sim.geneNames[g]), f(A[c * k + j], 2)]), [0, 1])), code: ["const A = new Float64Array(N * k);", 15] },
        { t: "PCA on the gene set alone", p: `Find the main patterns of co-variation among the set's genes across cells (principal components). Keep the first ${npcs}; then keep as many as needed to explain more than half of the variation: here ${maxComb}.`,
          b: table(["PC", "share of variation", "kept"], varExp.map((v, t) => [`PC${t + 1}`, `${(100 * v).toFixed(1)}%`, t < maxComb ? "✓" : ""]), [0, 1, 0]), code: ["const { values, vectors } = jacobiEigen(C, k);", 3] },
        { t: "Combine this cell's PC scores", p: `Take this cell's position on each kept PC, shift so the lowest value anywhere is 0, weight by each PC's share of variation, add, and take the square root.`,
          b: table(["PC", "this cell", "shifted", "weight"], emb.slice(0, maxComb).map((u, t) => [`PC${t + 1}`, f(u[c], 3), f(u[c] - gmin, 3), f(varExp[t], 3)]), [0, 1, 1, 1]) + kv([[f(agg, 4), "√(Σ shifted × weight)"]]), code: ["for (let t = 0; t < maxComb; t++) s += (emb[t][c] - gmin) * varExp[t];", 1] },
        { t: "Score = that × the set's mean expression in this cell", p: "Multiplying by the cell's average set expression ties the score to how much the set is actually expressed, not just to the pattern.",
          b: kv([[f(score, 4), `${f(agg, 3)} × ${f(mexp, 3)}`]]) + verify("scPS", score), code: ["out[c] = agg[c] * (m / k);", 1] },
      ];
    },
  };

  // ------------------------------------------------------------ last node: all cells
  function strip(scores) {
    const N = sim.N; let lo = Infinity, hi = -Infinity; for (const v of scores) { if (v < lo) lo = v; if (v > hi) hi = v; }
    const W = 400, x = v => 6 + 388 * (v - lo) / (hi - lo || 1), jit = i => ((i * 2654435761) % 1000) / 1000;
    let h = `<svg class="mg-strip" viewBox="0 0 ${W} 46" preserveAspectRatio="none">`;
    for (let i = 0; i < N; i++) h += `<circle class="${sim.group[i] ? "dot-t" : "dot-c"} cell${i === cfg.cell ? " sel" : ""}" data-c="${i}" cx="${x(scores[i]).toFixed(1)}" cy="${((sim.group[i] ? 27 : 6) + 12 * jit(i)).toFixed(1)}" r="${i === cfg.cell ? 4 : 2.2}"/>`;
    return h + `</svg>`;
  }

  // ------------------------------------------------------------ render
  async function build(el) {
    const m = el.dataset.method;
    el.innerHTML = `<p class="cp-note">Loading the real cells (~12 MB, once)…</p>`;
    if (!pool) {
      const [meta, buf] = await Promise.all([fetch(root + "data/naive_b/meta.json").then(r => r.json()), fetch(root + "data/naive_b/counts.bin").then(r => r.arrayBuffer())]);
      pool = S.decodePool(buf, meta);
    }
    const st = SETTINGS[cfg.setting];
    const key = `${st.scenario}`;
    if (key !== simKey) { sim = S.simulateRWSD(pool, { scenario: st.scenario, cellsPerGroup: 200, seed: 7 }); prep = S.prepare(sim); simKey = key; }
    set = S.makeSet(sim, { size: cfg.size, noise: st.noise }, S.makeRng(cfg.setSeed * 7919 + cfg.size));
    res = S.scoreAll(sim, prep, set);
    const steps = STEPS[m]();
    const r = res[m];
    const cellOpts = [[0, "C1 (control)"], [1, "C2 (control)"], [sim.n, "T1 (treatment)"], [sim.n + 1, "T2 (treatment)"]];
    if (!cellOpts.some(([v]) => v === cfg.cell)) cellOpts.push([cfg.cell, `${cellLabel(cfg.cell)} (clicked)`]);
    const ci = sim.cellInfo[cfg.cell];
    let h = `<div class="mg-ctrls">
      <label>Setting<select data-k="setting">${Object.entries(SETTINGS).map(([k, s]) => `<option value="${k}"${k === cfg.setting ? " selected" : ""}>${s.label}</option>`).join("")}</select></label>
      <label>Gene set size<select data-k="size">${[8, 25, 100].map(v => `<option${v === cfg.size ? " selected" : ""}>${v}</option>`).join("")}</select></label>
      <label>Follow this cell<select data-k="cell">${cellOpts.map(([v, t]) => `<option value="${v}"${v === cfg.cell ? " selected" : ""}>${t}</option>`).join("")}</select></label>
      <button type="button" class="btn ghost" data-act="reroll">New random gene set</button></div>
      <p class="cp-note">Data: the paper's real cells (naive B, GSE164381), 200 control + 200 treatment, as in the <a href="${root}playground.html">playground</a>. Following <b>${cellLabel(cfg.cell)}</b> — ${esc(ci.donor)}, <code>${esc(ci.barcode)}</code>.</p>
      <div class="mg-graph">`;
    steps.forEach((s, i) => {
      h += `<section class="pnode"><div class="pn-head"><span class="pn-num">${i + 1}</span><h3>${esc(s.t)}</h3></div><p class="pn-blurb">${s.p}</p><div class="pn-out">${s.b}</div>` +
        (s.code ? `<div class="pn-code"><button type="button" class="codebtn" data-find="${esc(s.code[0])}" data-lines="${s.code[1]}">code: static/scgsa.js</button></div><div class="mg-code" hidden></div>` : "") + `</section>`;
    });
    h += `<section class="pnode"><div class="pn-head"><span class="pn-num">${steps.length + 1}</span><h3>Do the same for every cell, then compare the groups</h3></div>
      <p class="pn-blurb">Every one of the ${sim.N} cells gets a score this way. The Wilcoxon test then asks whether treatment cells score differently from control cells. Click any dot to follow that cell instead.</p>
      <div class="pn-out"><div class="legend2"><span><i style="background:var(--muted)"></i>control</span><span><i style="background:var(--s1)"></i>treatment</span><span>big dot = the cell followed above</span></div>${strip(r.scores)}
      ${kv([[r.test.p < 1e-4 ? r.test.p.toExponential(1) : r.test.p.toFixed(4), "Wilcoxon p, control vs treatment"], [r.test.p < 0.05 ? (st.noise === 1 || st.scenario === 3 ? "false positive" : "detected") : (st.noise === 1 || st.scenario === 3 ? "correctly quiet" : "missed"), "verdict for this set"]])}</div></section></div>`;
    el.innerHTML = h;
  }

  function init() {
    const el = document.getElementById("mg");
    if (!el || !S) return;
    const rerun = () => build(el).catch(e => { el.innerHTML = `<p class="cp-note">Couldn't load the real cells: ${esc(e.message)}</p>`; });
    el.addEventListener("change", e => {
      const k = e.target.dataset.k; if (!k) return;
      cfg[k] = k === "setting" ? e.target.value : Number(e.target.value);
      rerun();
    });
    el.addEventListener("click", async e => {
      const b = e.target.closest(".codebtn");
      if (b) {
        const box = b.closest(".pnode").querySelector(".mg-code");
        if (!box.hidden) { box.hidden = true; return; }
        const s = await snippet(b.dataset.find.replace(/&quot;/g, '"').replace(/\\n/g, "\n").split("\n")[0], Number(b.dataset.lines));
        box.innerHTML = s ? `<p class="cp-note"><a href="${OURS}#L${s.start}-L${s.end}" target="_blank" rel="noopener">static/scgsa.js L${s.start}–${s.end}</a></p><pre class="code"><code>${esc(s.text)}</code></pre>` : `<p class="cp-note">Couldn't locate this code; reload the page.</p>`;
        box.hidden = false;
        return;
      }
      if (e.target.closest('[data-act="reroll"]')) { cfg.setSeed++; rerun(); return; }
      const dot = e.target.closest("circle.cell");
      if (dot) { cfg.cell = Number(dot.dataset.c); rerun(); el.scrollIntoView({ behavior: "smooth", block: "start" }); }
    });
    rerun();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
