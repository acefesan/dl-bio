// The Wilcoxon rank-sum test as a pipeline of steps, like the playground.
// It follows the interactive chart above it: static/wilcoxon.js sends a
// "wilcoxon:update" event with the pooled, ranked cells and the statistics
// (including the U value of every possible relabelling), and each node here
// shows that step's sample in → out with the live numbers.
(function () {
  const box = document.getElementById("wpipe");
  if (!box) return;
  const fmt = v => (Number.isInteger(v) ? String(v) : v.toFixed(1));
  const pfmt = p => (p < 0.001 ? "< 0.001" : p.toFixed(3));
  const tag = g => `<span class="role ${g === "Treatment" ? "trt" : "ctl"}">${g === "Treatment" ? "T" : "C"}</span>`;
  const table = (head, rows, num) => `<div class="table-wrap"><table class="mg-tab"><tr>${head.map((h, i) => `<th class="${num && num[i] ? "num" : ""}">${h}</th>`).join("")}</tr>` +
    rows.map(r => `<tr>${r.map((c, i) => `<td class="${num && num[i] ? "num" : ""}">${c}</td>`).join("")}</tr>`).join("") + `</table></div>`;
  const io = (a, b) => `<div class="sio-row"><div class="sio"><div class="sio-t">in</div>${a}</div><div class="sio-arrow">→</div><div class="sio"><div class="sio-t">out</div>${b}</div></div>`;
  const kv = items => `<div class="kvs">${items.map(([b, s]) => `<div class="kv"><b>${b}</b><span>${s}</span></div>`).join("")}</div>`;

  function nullHist(st) {
    const counts = new Map();
    for (const u of st.nullUs) counts.set(u, (counts.get(u) || 0) + 1);
    const xs = [...counts.keys()].sort((a, b) => a - b), max = Math.max(...counts.values());
    const W = 460, H = 170, L = 30, B = 30, T = 12, R = 10;
    const lo = 0, hi = st.n1 * st.n2, bw = Math.max(4, (W - L - R) / (hi - lo + 1) - 2);
    const X = u => L + (u - lo) / (hi - lo || 1) * (W - L - R - bw) + bw / 2;
    let s = `<svg class="mg-plot" viewBox="0 0 ${W} ${H}">`;
    for (const u of xs) {
      const h = (H - T - B) * counts.get(u) / max;
      const ext = Math.abs(u - st.center) + 1e-9 >= st.deviation;
      s += `<rect class="${ext ? "bar ext" : "bar"}" x="${(X(u) - bw / 2).toFixed(1)}" y="${(H - B - h).toFixed(1)}" width="${bw.toFixed(1)}" height="${h.toFixed(1)}"><title>U = ${fmt(u)}: ${counts.get(u)} of ${st.total} relabellings</title></rect>`;
    }
    s += `<line class="ax" x1="${L}" x2="${W - R}" y1="${H - B}" y2="${H - B}"/>`;
    s += `<line class="obs" x1="${X(st.u)}" x2="${X(st.u)}" y1="${T - 4}" y2="${H - B}"/><text class="obsl" x="${X(st.u) > W / 2 ? X(st.u) - 4 : X(st.u) + 4}" y="${T + 6}" text-anchor="${X(st.u) > W / 2 ? "end" : "start"}">observed U = ${fmt(st.u)}</text>`;
    s += `<line class="mid" x1="${X(st.center)}" x2="${X(st.center)}" y1="${T}" y2="${H - B}"/>`;
    for (let k = 0; k <= 4; k++) { const u = lo + (hi - lo) * k / 4; s += `<text x="${X(u)}" y="${H - 14}" text-anchor="middle">${fmt(Math.round(u * 2) / 2)}</text>`; }
    return s + `<text x="${(L + W) / 2}" y="${H - 2}" text-anchor="middle">U under every possible relabelling</text></svg>`;
  }

  function render({ pooled, stats: st }) {
    const ctl = pooled.filter(d => d.group === "Control"), trt = pooled.filter(d => d.group === "Treatment");
    const N = pooled.length, hasTies = new Set(pooled.map(d => d.value)).size < N;
    const minSum = st.n1 * (st.n1 + 1) / 2;
    const steps = [
      { t: "Two groups of numbers", p: "The input: one score per cell, in two labelled groups. Nothing else about the cells enters the test.",
        b: kv([[ctl.map(d => fmt(d.value)).join(", "), `control (n2 = ${st.n2})`], [trt.map(d => fmt(d.value)).join(", "), `treatment (n1 = ${st.n1})`]]) },
      { t: "Pool them", p: "Forget the labels for a moment and put all the scores in one list. The labels travel with the scores, but the next step looks only at the values.",
        b: io(kv([[`${st.n2} + ${st.n1}`, "scores in two lists"]]), `<p class="cp-note">${pooled.slice().sort((a, b) => a.id.localeCompare(b.id)).map(d => `${fmt(d.value)} ${tag(d.group)}`).join(" · ")}</p>`) },
      { t: "Sort and rank", p: `Sort the pooled scores from lowest to highest and number them 1 to ${N}. Equal scores share the average of their ranks.${hasTies ? " This data has ties." : ""} Only the order survives; the distances between scores are gone.`,
        b: io(`<p class="cp-note">${N} scores</p>`, table(["score", "group", "rank"], pooled.map(d => [fmt(d.value), tag(d.group), fmt(d.rank)]), [1, 0, 1])) },
      { t: "Add up the treatment ranks", p: "R1 = the sum of the ranks that landed on treatment cells. If treatment scores tend to be high, they hold the high ranks and R1 is large.",
        b: kv([[fmt(st.rankSum), `R1 = ${st.treatmentRanks.map(fmt).join(" + ")}`]]) },
      { t: "Shift R1 to U", p: `Subtract the smallest R1 possible, n1(n1 + 1)/2 = ${fmt(minSum)} (all treatment cells at the bottom). The result, U, runs from 0 to n1 × n2 = ${st.n1 * st.n2}, and counts the control–treatment pairs where the treatment cell scores higher.`,
        b: kv([[fmt(st.u), `U = ${fmt(st.rankSum)} − ${fmt(minSum)}`], [`0 – ${st.n1 * st.n2}`, "possible range"], [fmt(st.center), "midpoint if labels mean nothing (n1 × n2 / 2)"]]) },
      { t: "Build the null distribution", p: `If the labels meant nothing, any ${st.n1} of the ${N} cells could have been the treatment group. List all ${st.total} ways, compute U for each: that is the null distribution. Red bars are as far from the midpoint as the observed U, or farther.`,
        b: nullHist(st) },
      { t: "Two-sided p-value", p: "The share of relabellings at least as extreme as the one observed, in either direction. Small p: these labels are unusual for random labels.",
        b: kv([[pfmt(st.p), `p = ${st.extreme} ÷ ${st.total} relabellings`], [st.p < 0.05 ? "yes" : "no", "below 0.05?"]]) + (hasTies ? `<p class="cp-note">With ties, this exact count still holds (it uses the tied ranks); textbook tables assume no ties.</p>` : "") },
      { t: "How big is the shift?", p: "U ÷ (n1 × n2) is the chance that a randomly picked treatment cell outscores a randomly picked control cell (ties count half). 0.5 = no tendency.",
        b: kv([[(st.u / (st.n1 * st.n2)).toFixed(2), "P(treatment > control)"]]) },
    ];
    box.innerHTML = steps.map((s, i) => `<section class="pnode"><div class="pn-head"><span class="pn-num">${i + 1}</span><h3>${s.t}</h3></div><p class="pn-blurb">${s.p}</p><div class="pn-out">${s.b}</div></section>`).join("");
  }

  document.addEventListener("wilcoxon:update", e => render(e.detail));
})();
