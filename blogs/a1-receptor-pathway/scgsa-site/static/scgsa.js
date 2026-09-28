// Toy re-implementation of the seven scGSA methods benchmarked in
// Wang & Thakar 2024 (NARGAB lqae124), plus a small simulator of the
// paper's four scenarios. Formulas follow the installed package source
// (Seurat 5.5.1, AUCell 1.28.0, UCell 2.10.1, GSVA 2.0.0, irGSEA wrappers,
// scPS_2024/2025.R, JASMINE_V1). This is for intuition, not for results:
// it is validated against the real R packages by validate/ in this folder.
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.SCGSA = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  // ------------------------------------------------------------ random
  function mulberry32(a) {
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function makeRng(seed) {
    const u = mulberry32(seed >>> 0);
    let spare = null;
    const rng = {
      uniform() { let x; do { x = u(); } while (x === 0); return x; },
      normal() {
        if (spare !== null) { const s = spare; spare = null; return s; }
        let a, b, s;
        do { a = 2 * u() - 1; b = 2 * u() - 1; s = a * a + b * b; } while (s >= 1 || s === 0);
        const m = Math.sqrt(-2 * Math.log(s) / s);
        spare = b * m;
        return a * m;
      },
      gamma(shape, scale) {
        if (shape < 1) {
          return rng.gamma(shape + 1, 1) * Math.pow(rng.uniform(), 1 / shape) * scale;
        }
        const d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
        for (;;) {
          let x, v;
          do { x = rng.normal(); v = 1 + c * x; } while (v <= 0);
          v = v * v * v;
          const U = rng.uniform();
          if (U < 1 - 0.0331 * x * x * x * x) return d * v * scale;
          if (Math.log(U) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v * scale;
        }
      },
      poisson(lam) {
        if (lam <= 0) return 0;
        if (lam < 30) {
          const L = Math.exp(-lam);
          let k = 0, p = 1;
          do { k++; p *= rng.uniform(); } while (p > L);
          return k - 1;
        }
        return Math.max(0, Math.round(lam + Math.sqrt(lam) * rng.normal()));
      },
      int(n) { return Math.floor(u() * n); },
    };
    return rng;
  }

  function sampleWithout(rng, pool, k) {
    const a = Array.from(pool);
    k = Math.min(k, a.length);
    for (let i = 0; i < k; i++) {
      const j = i + rng.int(a.length - i);
      const t = a[i]; a[i] = a[j]; a[j] = t;
    }
    return a.slice(0, k);
  }

  // ------------------------------------------------------------ simulate
  // Scenario 1: +20% on densely expressed genes (>= 90% of cells).
  // Scenario 2: +20% on the same number of randomly chosen genes.
  // Scenarios 3/4: no signal; ~2% extra genes exist only in treatment cells,
  // densely (3) or sparsely (4) expressed.
  //
  // inject = "log" (default) follows the paper's real-world simulation: the
  // +20% and the extra genes are applied to the already log-normalised
  // matrix, so no other gene's value moves. inject = "counts" applies them
  // before normalisation instead, which leaks into every other gene through
  // the library size.
  function simulate(opt) {
    const o = Object.assign({
      genes: 3000, cellsPerGroup: 200, scenario: 1, seed: 7,
      denseSignalShare: 1, foldChange: 1.2, csFraction: 0.02,
      meanLog: -1.2, sdLog: 1.6,
      dispersion: 0.2, inject: "log",
    }, opt || {});
    const rng = makeRng(o.seed);
    const G0 = o.genes;
    const nCS = o.scenario >= 3 ? Math.round(o.csFraction * G0) : 0;
    const G = G0 + nCS;
    const n = o.cellsPerGroup;
    const N = 2 * n;
    const phi = o.dispersion;

    const mu = new Float64Array(G);
    for (let g = 0; g < G0; g++) mu[g] = Math.exp(o.meanLog + o.sdLog * rng.normal());
    for (let g = G0; g < G; g++) {
      mu[g] = o.scenario === 3
        ? 2.5 * Math.exp(rng.uniform() * Math.log(8))
        : Math.exp(o.meanLog + o.sdLog * rng.normal());
    }

    // "dense" = expected detection >= 90% (NB zero probability <= 0.1)
    const dense = [];
    for (let g = 0; g < G0; g++) {
      const p0 = Math.pow(1 + mu[g] * phi, -1 / phi);
      if (p0 <= 0.1) dense.push(g);
    }
    const nSignal = Math.max(1, Math.round(o.denseSignalShare * dense.length));
    const isSignal = new Uint8Array(G);
    let signal = [];
    if (o.scenario === 1) signal = sampleWithout(rng, dense, nSignal);
    else if (o.scenario === 2) signal = sampleWithout(rng, Array.from({ length: G0 }, (_, i) => i), nSignal);
    for (const g of signal) isSignal[g] = 1;

    const group = new Uint8Array(N);
    const size = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      group[c] = c < n ? 0 : 1;
      size[c] = Math.exp(0.3 * rng.normal());
    }

    const atCounts = o.inject === "counts";
    const counts = new Float32Array(N * G);
    for (let c = 0; c < N; c++) {
      for (let g = 0; g < G; g++) {
        if (g >= G0 && group[c] === 0) continue;
        let m = mu[g] * size[c];
        if (atCounts && group[c] === 1 && isSignal[g]) m *= o.foldChange;
        counts[c * G + g] = rng.poisson(rng.gamma(1 / phi, m * phi));
      }
    }

    // Seurat LogNormalize: log1p(count / library size * 1e4). In "log" mode
    // the library size is taken over the shared genes only, so the extra
    // treatment-only genes do not rescale anything else.
    const data = new Float32Array(N * G);
    for (let c = 0; c < N; c++) {
      let lib = 0;
      const top = atCounts ? G : G0;
      for (let g = 0; g < top; g++) lib += counts[c * G + g];
      const f = lib > 0 ? 1e4 / lib : 0;
      for (let g = 0; g < G; g++) {
        let v = Math.log1p(counts[c * G + g] * f);
        if (!atCounts && group[c] === 1 && isSignal[g]) v *= o.foldChange;
        data[c * G + g] = v;
      }
    }

    const nonSignal = [];
    for (let g = 0; g < G0; g++) if (!isSignal[g]) nonSignal.push(g);
    const csGenes = [];
    for (let g = G0; g < G; g++) csGenes.push(g);

    return { opt: o, G, G0, N, n, nCS, group, counts, data, mu, signal, nonSignal, csGenes, isSignal, dense };
  }

  // ------------------------------------------------------------ ranking helpers
  // rank descending (1 = highest). ties: "average" or "random"
  function rankDesc(vals, ties, rng) {
    const n = vals.length;
    const idx = Array.from({ length: n }, (_, i) => i);
    let key = null;
    if (ties === "random") {
      key = new Float64Array(n);
      for (let i = 0; i < n; i++) key[i] = rng.uniform();
      idx.sort((a, b) => (vals[b] - vals[a]) || (key[a] - key[b]));
      const r = new Float32Array(n);
      for (let p = 0; p < n; p++) r[idx[p]] = p + 1;
      return r;
    }
    idx.sort((a, b) => vals[b] - vals[a]);
    const r = new Float32Array(n);
    let p = 0;
    while (p < n) {
      let q = p;
      while (q + 1 < n && vals[idx[q + 1]] === vals[idx[p]]) q++;
      const avg = (p + q) / 2 + 1;
      for (let t = p; t <= q; t++) r[idx[t]] = avg;
      p = q + 1;
    }
    return r;
  }

  // ascending average rank (1 = lowest), as colRanks(ties = "average")
  function rankAscAvg(vals) {
    const n = vals.length;
    const idx = Array.from({ length: n }, (_, i) => i);
    idx.sort((a, b) => vals[a] - vals[b]);
    const r = new Float64Array(n);
    let p = 0;
    while (p < n) {
      let q = p;
      while (q + 1 < n && vals[idx[q + 1]] === vals[idx[p]]) q++;
      const avg = (p + q) / 2 + 1;
      for (let t = p; t <= q; t++) r[idx[t]] = avg;
      p = q + 1;
    }
    return r;
  }

  // ------------------------------------------------------------ prepare
  // Per-cell work that does not depend on the gene set.
  function prepare(sim, opt) {
    const o = Object.assign({ tieSeed: 11, ucellMaxRankFrac: 1500 / 20000 }, opt || {});
    const { G, N, data } = sim;
    const rng = makeRng(o.tieSeed);
    const rAUC = new Float32Array(N * G);
    const rU = new Float32Array(N * G);
    const ssPos = new Float32Array(N * G);
    const ssRa = new Float32Array(N * G);
    const rNZ = new Float32Array(N * G);
    const nnz = new Int32Array(N);
    const colSum = new Float64Array(N);

    for (let c = 0; c < N; c++) {
      const v = data.subarray(c * G, (c + 1) * G);

      const ra = rankDesc(v, "random", rng);
      const ru = rankDesc(v, "average");
      rAUC.set(ra, c * G);
      rU.set(ru, c * G);

      // ssGSEA: ascending average ranks, then mode(R) <- "integer" (truncates)
      const Rasc = rankAscAvg(v);
      const Rint = new Int32Array(G);
      for (let g = 0; g < G; g++) Rint[g] = Math.trunc(Rasc[g]);
      const order = Array.from({ length: G }, (_, i) => i)
        .sort((a, b) => (Rint[b] - Rint[a]) || (a - b));
      for (let p = 0; p < G; p++) ssPos[c * G + order[p]] = p + 1;
      for (let g = 0; g < G; g++) ssRa[c * G + g] = Math.pow(Math.abs(Rint[g]), 0.25);

      // JASMINE: rank among non-zero genes only (ascending, average ties)
      const nzIdx = [];
      let s = 0;
      for (let g = 0; g < G; g++) { if (v[g] !== 0) nzIdx.push(g); s += v[g]; }
      colSum[c] = s;
      nnz[c] = nzIdx.length;
      const nzVals = nzIdx.map(g => v[g]);
      const nzRank = rankAscAvg(nzVals);
      for (let i = 0; i < nzIdx.length; i++) rNZ[c * G + nzIdx[i]] = nzRank[i];
    }

    // AddModuleScore: 24 equal-frequency bins on each gene's mean expression
    const avg = new Float64Array(G);
    for (let c = 0; c < N; c++) for (let g = 0; g < G; g++) avg[g] += data[c * G + g];
    for (let g = 0; g < G; g++) avg[g] /= N;
    const nbin = 24;
    const binRng = makeRng(1);
    const jitter = new Float64Array(G);
    for (let g = 0; g < G; g++) jitter[g] = binRng.normal() / 1e30;
    const ord = Array.from({ length: G }, (_, i) => i).sort((a, b) => (avg[a] + jitter[a]) - (avg[b] + jitter[b]));
    const bin = new Int32Array(G);
    for (let p = 0; p < G; p++) bin[ord[p]] = Math.min(nbin - 1, Math.floor(p * nbin / G));
    const bins = Array.from({ length: nbin }, () => []);
    for (let g = 0; g < G; g++) bins[bin[g]].push(g);

    const aucMaxRank = Math.ceil(0.05 * G);
    const ucellMaxRank = Math.max(10, Math.round(o.ucellMaxRankFrac * G));

    return { rAUC, rU, ssPos, ssRa, rNZ, nnz, colSum, avg, bin, bins, aucMaxRank, ucellMaxRank };
  }

  // ------------------------------------------------------------ methods
  function addModuleScore(sim, prep, set, opt) {
    const ctrl = (opt && opt.ctrl) || 100;
    const rng = makeRng(1);
    const { G, N, data } = sim;
    const ctrlSet = new Set();
    for (const g of set) {
      const pool = prep.bins[prep.bin[g]];
      for (const x of sampleWithout(rng, pool, Math.min(ctrl, pool.length))) ctrlSet.add(x);
    }
    const ctl = Array.from(ctrlSet);
    const out = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      let a = 0, b = 0;
      for (const g of set) a += data[c * G + g];
      for (const g of ctl) b += data[c * G + g];
      out[c] = a / set.length - b / ctl.length;
    }
    return out;
  }

  function aucell(sim, prep, set) {
    const { G, N } = sim;
    const thr = Math.round(prep.aucMaxRank);
    const k = set.length;
    let maxAUC = 0;
    {
      const xs = [];
      for (let i = 1; i <= k; i++) if (i < thr) xs.push(i);
      for (let i = 0; i < xs.length; i++) {
        const next = i + 1 < xs.length ? xs[i + 1] : thr;
        maxAUC += (next - xs[i]) * (i + 1);
      }
    }
    const out = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      const xs = [];
      for (const g of set) { const r = prep.rAUC[c * G + g]; if (r < thr) xs.push(r); }
      xs.sort((a, b) => a - b);
      let s = 0;
      for (let i = 0; i < xs.length; i++) {
        const next = i + 1 < xs.length ? xs[i + 1] : thr;
        s += (next - xs[i]) * (i + 1);
      }
      out[c] = maxAUC > 0 ? s / maxAUC : 0;
    }
    return out;
  }

  function ucell(sim, prep, set) {
    const { G, N } = sim;
    const maxRank = prep.ucellMaxRank;
    const k = set.length;
    const minSum = k * (k + 1) / 2;
    const out = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      let sum = 0, anySig = false;
      for (const g of set) {
        let r = prep.rU[c * G + g];
        if (r >= maxRank) r = maxRank; else anySig = true;
        sum += r;
      }
      out[c] = anySig ? 1 - (sum - minSum) / (k * maxRank - minSum) : 0;
    }
    return out;
  }

  function ssgsea(sim, prep, set) {
    const { G, N } = sim;
    const n = G, k = set.length;
    const out = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      let wsum = 0, wpos = 0, psum = 0;
      for (const g of set) {
        const pos = prep.ssPos[c * G + g];
        const w = prep.ssRa[c * G + g];
        wsum += w;
        wpos += w * (n - pos + 1);
        psum += (n - pos + 1);
      }
      const inStep = wsum > 0 ? wpos / wsum : 0;
      const outStep = (n * (n + 1) / 2 - psum) / (n - k);
      out[c] = inStep - outStep;
    }
    return out;
  }

  function minmax(x) {
    let lo = Infinity, hi = -Infinity;
    for (const v of x) { if (v < lo) lo = v; if (v > hi) hi = v; }
    const d = hi - lo;
    return x.map(v => (d > 0 ? (v - lo) / d : 0));
  }

  function jasmine(sim, prep, set) {
    const { G, N } = sim;
    const k = set.length;
    const RM = new Float64Array(N);
    const OR = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      let rs = 0, m = 0;
      for (const g of set) {
        const r = prep.rNZ[c * G + g];
        if (r > 0) { rs += r; m++; }
      }
      RM[c] = m > 0 ? (rs / m) / prep.nnz[c] : 0;

      let SE = m;
      let NE = prep.nnz[c] - SE;
      let SN = k - SE;
      if (SN === 0) SN = 1;
      if (NE === 0) NE = 1;
      let NN = G - (NE + SE);
      NN = NN - SN;
      OR[c] = (SE * NN) / (SN * NE);
    }
    const a = minmax(RM), b = minmax(OR);
    return a.map((v, i) => (v + b[i]) / 2);
  }

  function scse(sim, prep, set) {
    const { G, N, data } = sim;
    const out = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      let s = 0;
      for (const g of set) s += data[c * G + g];
      out[c] = prep.colSum[c] > 0 ? 100 * s / prep.colSum[c] : 0;
    }
    return out;
  }

  // cyclic Jacobi eigen-decomposition of a symmetric n x n matrix
  function jacobiEigen(S, n) {
    const a = Float64Array.from(S);
    const v = new Float64Array(n * n);
    for (let i = 0; i < n; i++) v[i * n + i] = 1;
    for (let sweep = 0; sweep < 80; sweep++) {
      let off = 0, diag = 0;
      for (let p = 0; p < n; p++) {
        diag += a[p * n + p] * a[p * n + p];
        for (let q = p + 1; q < n; q++) off += a[p * n + q] * a[p * n + q];
      }
      if (off <= 1e-22 * Math.max(diag, 1e-300)) break;
      for (let p = 0; p < n - 1; p++) {
        for (let q = p + 1; q < n; q++) {
          const apq = a[p * n + q];
          if (Math.abs(apq) < 1e-300) continue;
          const theta = (a[q * n + q] - a[p * n + p]) / (2 * apq);
          const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          const cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
          for (let r = 0; r < n; r++) {
            const arp = a[r * n + p], arq = a[r * n + q];
            a[r * n + p] = cs * arp - sn * arq;
            a[r * n + q] = sn * arp + cs * arq;
          }
          for (let r = 0; r < n; r++) {
            const apr = a[p * n + r], aqr = a[q * n + r];
            a[p * n + r] = cs * apr - sn * aqr;
            a[q * n + r] = sn * apr + cs * aqr;
          }
          for (let r = 0; r < n; r++) {
            const vrp = v[r * n + p], vrq = v[r * n + q];
            v[r * n + p] = cs * vrp - sn * vrq;
            v[r * n + q] = sn * vrp + cs * vrq;
          }
        }
      }
    }
    const values = new Float64Array(n);
    for (let i = 0; i < n; i++) values[i] = a[i * n + i];
    return { values, vectors: v };
  }

  // scPS: PCA on the scaled gene-set matrix; weight shifted PC scores by
  // variance explained up to 50% cumulative; sqrt; times mean expression.
  function scps(sim, prep, set, opt) {
    const o = Object.assign({ bug2024: false, flipPC1: false }, opt || {});
    const { G, N, data } = sim;
    const k = set.length;
    const A = new Float64Array(N * k);
    for (let j = 0; j < k; j++) {
      const g = set[j];
      let s = 0;
      for (let c = 0; c < N; c++) s += data[c * G + g];
      const m = s / N;
      let vv = 0;
      for (let c = 0; c < N; c++) { const d = data[c * G + g] - m; vv += d * d; }
      const sd = Math.sqrt(vv / (N - 1));
      for (let c = 0; c < N; c++) {
        let z = sd > 0 ? (data[c * G + g] - m) / sd : 0;
        if (z > 10) z = 10;
        A[c * k + j] = z;
      }
    }
    const C = new Float64Array(k * k);
    for (let c = 0; c < N; c++) {
      for (let i = 0; i < k; i++) {
        const ai = A[c * k + i];
        if (ai === 0) continue;
        for (let j = i; j < k; j++) C[i * k + j] += ai * A[c * k + j];
      }
    }
    for (let i = 0; i < k; i++) for (let j = 0; j < i; j++) C[i * k + j] = C[j * k + i];

    const { values, vectors } = jacobiEigen(C, k);
    const order = Array.from({ length: k }, (_, i) => i).sort((a, b) => values[b] - values[a]);
    const npcs = Math.max(1, Math.min(10, k - 1, N - 1));

    const emb = [];
    const sdev = [];
    for (let t = 0; t < npcs; t++) {
      const e = order[t];
      const lam = Math.max(values[e], 0);
      const d = Math.sqrt(lam);
      let big = 0, sign = 1;
      for (let i = 0; i < k; i++) {
        const x = vectors[i * k + e];
        if (Math.abs(x) > Math.abs(big)) big = x;
      }
      sign = big < 0 ? -1 : 1;
      if (t === 0 && o.flipPC1) sign = -sign;
      const u = new Float64Array(N);
      for (let c = 0; c < N; c++) {
        let s = 0;
        for (let i = 0; i < k; i++) s += A[c * k + i] * vectors[i * k + e];
        u[c] = d > 0 ? sign * s / d : 0;
      }
      emb.push(u);
      sdev.push(d / Math.sqrt(Math.max(1, N - 1)));
    }
    const tot = sdev.reduce((s, x) => s + x * x, 0);
    const varExp = sdev.map(x => (tot > 0 ? x * x / tot : 0));
    let cum = 0, maxComb = npcs;
    for (let t = 0; t < npcs; t++) { cum += varExp[t]; if (cum > 0.5) { maxComb = t + 1; break; } }

    let gmin = Infinity;
    for (const u of emb) for (const x of u) if (x < gmin) gmin = x;

    const agg = new Float64Array(N);
    if (o.bug2024 && maxComb === 1) {
      let s = 0;
      for (let c = 0; c < N; c++) s += (emb[0][c] - gmin) * varExp[0];
      agg.fill(Math.sqrt(s));
    } else {
      for (let c = 0; c < N; c++) {
        let s = 0;
        for (let t = 0; t < maxComb; t++) s += (emb[t][c] - gmin) * varExp[t];
        agg[c] = Math.sqrt(s);
      }
    }
    const out = new Float64Array(N);
    for (let c = 0; c < N; c++) {
      let m = 0;
      for (const g of set) m += data[c * G + g];
      out[c] = agg[c] * (m / k);
    }
    out.meta = { npcs, maxComb, varExp };
    return out;
  }

  const METHODS = {
    AddModuleScore: addModuleScore,
    AUCell: aucell,
    UCell: ucell,
    ssGSEA: ssgsea,
    JASMINE: jasmine,
    SCSE: scse,
    scPS: scps,
  };

  // ------------------------------------------------------------ statistics
  function erfc(x) {
    const z = Math.abs(x);
    const t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
      t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
      t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }

  // two-sided Wilcoxon rank-sum, normal approximation with tie and
  // continuity correction (what wilcox.test does at these sample sizes)
  function wilcox(x, group) {
    const N = x.length;
    const idx = Array.from({ length: N }, (_, i) => i).sort((a, b) => x[a] - x[b]);
    const r = new Float64Array(N);
    let tieTerm = 0;
    let p = 0;
    while (p < N) {
      let q = p;
      while (q + 1 < N && x[idx[q + 1]] === x[idx[p]]) q++;
      const avg = (p + q) / 2 + 1;
      for (let t = p; t <= q; t++) r[idx[t]] = avg;
      const tl = q - p + 1;
      if (tl > 1) tieTerm += tl * tl * tl - tl;
      p = q + 1;
    }
    let n1 = 0, R1 = 0;
    for (let i = 0; i < N; i++) if (group[i] === 1) { n1++; R1 += r[i]; }
    const n2 = N - n1;
    const W = R1 - n1 * (n1 + 1) / 2;
    const mean = n1 * n2 / 2;
    const varW = n1 * n2 / 12 * ((N + 1) - tieTerm / (N * (N - 1)));
    if (varW <= 0) return { p: 1, z: 0, auc: 0.5 };
    const dev = W - mean;
    const corr = dev > 0 ? 0.5 : dev < 0 ? -0.5 : 0;
    const z = (dev - corr) / Math.sqrt(varW);
    const pv = Math.min(1, erfc(Math.abs(z) / Math.SQRT2));
    return { p: pv, z, auc: W / (n1 * n2) };
  }

  function bh(ps) {
    const m = ps.length;
    const idx = ps.map((p, i) => [p, i]).sort((a, b) => b[0] - a[0]);
    const out = new Array(m);
    let run = 1;
    idx.forEach(([p, i], j) => {
      const rank = m - j;
      run = Math.min(run, p * m / rank);
      out[i] = Math.min(1, run);
    });
    return out;
  }

  // ------------------------------------------------------------ gene sets
  function makeSet(sim, spec, rng) {
    const size = spec.size;
    if (sim.opt.scenario <= 2) {
      const nSig = Math.min(sim.signal.length, Math.round(size * (1 - spec.noise)));
      return sampleWithout(rng, sim.signal, nSig)
        .concat(sampleWithout(rng, sim.nonSignal, size - nSig));
    }
    const nCS = spec.withCS ? Math.max(1, Math.round(size * spec.csFrac)) : 0;
    return sampleWithout(rng, sim.csGenes, nCS)
      .concat(sampleWithout(rng, sim.nonSignal, size - nCS));
  }

  function scoreAll(sim, prep, set, opt) {
    const res = {};
    for (const [name, fn] of Object.entries(METHODS)) {
      const t0 = Date.now();
      const s = fn(sim, prep, set, opt && opt[name]);
      res[name] = { scores: s, test: wilcox(s, sim.group), ms: Date.now() - t0, meta: s.meta };
    }
    return res;
  }

  function recovery(sim, prep, spec, nSets, seed, opt) {
    const rng = makeRng(seed || 99);
    const ps = {};
    for (const k of Object.keys(METHODS)) ps[k] = [];
    for (let i = 0; i < nSets; i++) {
      const set = makeSet(sim, spec, rng);
      const r = scoreAll(sim, prep, set, opt);
      for (const k of Object.keys(METHODS)) ps[k].push(r[k].test.p);
    }
    const out = {};
    for (const k of Object.keys(METHODS)) {
      const adj = bh(ps[k]);
      out[k] = adj.filter(p => p < 0.05).length / nSets;
    }
    return out;
  }

  function dataSummary(sim) {
    const { G, G0, N, data, group } = sim;
    const detected = [0, 0], cells = [0, 0];
    for (let c = 0; c < N; c++) {
      let d = 0;
      for (let g = 0; g < G; g++) if (data[c * G + g] > 0) d++;
      detected[group[c]] += d;
      cells[group[c]]++;
    }
    const pctExpr = gs => {
      const v = gs.map(g => {
        let e = 0;
        for (let c = 0; c < N; c++) if (data[c * G + g] > 0) e++;
        return 100 * e / N;
      }).sort((a, b) => a - b);
      return v.length ? v[Math.floor(v.length / 2)] : NaN;
    };
    let zeros = 0;
    for (let i = 0; i < data.length; i++) if (data[i] === 0) zeros++;
    return {
      genesPerCell: [detected[0] / cells[0], detected[1] / cells[1]],
      medianPctCellsExpressingSignal: sim.signal.length ? pctExpr(sim.signal) : NaN,
      medianPctCellsExpressingAll: pctExpr(Array.from({ length: G0 }, (_, i) => i)),
      nSignal: sim.signal.length,
      nDense: sim.dense.length,
      sparsity: zeros / data.length,
    };
  }

  return { makeRng, simulate, prepare, METHODS, scoreAll, recovery, makeSet, wilcox, bh, dataSummary, jacobiEigen };
});
