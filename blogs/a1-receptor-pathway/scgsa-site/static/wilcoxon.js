(function () {
  const NS = "http://www.w3.org/2000/svg";
  const scoreSvg = document.getElementById("w-score-chart");
  const rankSvg = document.getElementById("w-rank-chart");
  if (!scoreSvg || !rankSvg) return;

  let serial = 0;
  let control = makeGroup("C", [2.0, 4.8, 7.0]);
  let treatment = makeGroup("T", [3.0, 5.6, 8.0]);
  let drag = null;
  // Each move redraws the chart, replacing the dragged element, so the drag is
  // tracked on the window rather than on the dot (pointer capture would be lost).
  window.addEventListener("pointermove", event => {
    if (drag) moveDot(event.clientX, scoreSvg.getBoundingClientRect());
  });
  window.addEventListener("pointerup", () => { drag = null; });
  window.addEventListener("pointercancel", () => { drag = null; });

  const presets = {
    overlap: [[2.0, 4.8, 7.0, 8.5, 9.2, 9.7], [3.0, 5.6, 8.0, 8.8, 9.4, 9.9]],
    separated: [[1.5, 2.8, 4.1, 4.8, 5.4, 5.8], [6.2, 7.1, 8.0, 8.7, 9.3, 9.8]],
    reverse: [[6.2, 7.1, 8.0, 8.7, 9.3, 9.8], [1.5, 2.8, 4.1, 4.8, 5.4, 5.8]],
    ties: [[2.0, 4.0, 6.0, 6.0, 8.0, 9.0], [2.0, 4.0, 6.0, 7.0, 8.0, 9.0]],
  };

  function makeGroup(prefix, values) {
    return values.map(value => ({ id: `${prefix}${++serial}`, value }));
  }

  function resizeGroup(group, count, prefix) {
    while (group.length < count) {
      const previous = group.length ? group[group.length - 1].value : 5;
      group.push({ id: `${prefix}${++serial}`, value: clamp(previous + 0.7) });
    }
    group.length = count;
  }

  function clamp(value) { return Math.max(0, Math.min(10, Math.round(value * 10) / 10)); }
  const PLOT_L = 116; // left edge of the score axis; room for the row labels
  function sx(value, width) { return PLOT_L + value / 10 * (width - PLOT_L - 30); }
  function el(name, attrs, text) {
    const node = document.createElementNS(NS, name);
    for (const [key, value] of Object.entries(attrs || {})) node.setAttribute(key, value);
    if (text !== undefined) node.textContent = text;
    return node;
  }

  // One shape per group everywhere: control = square, treatment = circle.
  // Fill/stroke are also set as SVG attributes so the groups stay distinct
  // even if the stylesheet is missing or stale (CSS overrides them).
  const LOOK = { Control: { fill: "#e3ecf6", stroke: "#3b7bbf" }, Treatment: { fill: "#f8e6d6", stroke: "#d0762b" } };
  function shape(group, half) {
    const look = LOOK[group];
    const common = { fill: look.fill, stroke: look.stroke, "stroke-width": 2 };
    return group === "Treatment"
      ? el("circle", Object.assign({ r: half }, common))
      : el("rect", Object.assign({ x: -half, y: -half, width: 2 * half, height: 2 * half, rx: Math.round(half / 4) }, common));
  }

  function ranked() {
    const pooled = [
      ...control.map(d => ({ ...d, group: "Control" })),
      ...treatment.map(d => ({ ...d, group: "Treatment" })),
    ].sort((a, b) => a.value - b.value || a.id.localeCompare(b.id));
    for (let start = 0; start < pooled.length;) {
      let end = start + 1;
      while (end < pooled.length && pooled[end].value === pooled[start].value) end++;
      const averageRank = ((start + 1) + end) / 2;
      for (let i = start; i < end; i++) pooled[i].rank = averageRank;
      start = end;
    }
    return pooled;
  }

  function combinations(n, k) {
    const out = [];
    function walk(start, chosen) {
      if (chosen.length === k) { out.push(chosen.slice()); return; }
      for (let i = start; i <= n - (k - chosen.length); i++) {
        chosen.push(i); walk(i + 1, chosen); chosen.pop();
      }
    }
    walk(0, []);
    return out;
  }

  function statistics(pooled) {
    const n1 = treatment.length, n2 = control.length;
    const treatmentRanks = pooled.filter(d => d.group === "Treatment").map(d => d.rank);
    const rankSum = treatmentRanks.reduce((a, b) => a + b, 0);
    const u = rankSum - n1 * (n1 + 1) / 2;
    const center = n1 * n2 / 2;
    const deviation = Math.abs(u - center);
    const assignments = combinations(pooled.length, n1);
    let extreme = 0;
    for (const indices of assignments) {
      const r = indices.reduce((sum, i) => sum + pooled[i].rank, 0);
      const candidate = r - n1 * (n1 + 1) / 2;
      if (Math.abs(candidate - center) + 1e-9 >= deviation) extreme++;
    }
    return { n1, n2, treatmentRanks, rankSum, u, center, p: extreme / assignments.length };
  }

  function drawAxes(svg, width, y, label) {
    svg.append(el("line", { x1: PLOT_L, y1: y, x2: width - 30, y2: y, class: "wil-axis" }));
    for (let v = 0; v <= 10; v += 2) {
      const x = sx(v, width);
      svg.append(el("line", { x1: x, y1: y, x2: x, y2: y + 6, class: "wil-tick" }));
      svg.append(el("text", { x, y: y + 23, class: "wil-tick-label", "text-anchor": "middle" }, String(v)));
    }
    svg.append(el("text", { x: (PLOT_L + width - 30) / 2, y: y + 43, class: "wil-axis-label", "text-anchor": "middle" }, label));
  }

  function drawScoreChart(pooled) {
    const width = Math.max(560, scoreSvg.parentElement.clientWidth || 800);
    const height = 250;
    scoreSvg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    scoreSvg.replaceChildren();
    const rows = [{ name: "Control", data: control, y: 65 }, { name: "Treatment", data: treatment, y: 145 }];
    for (const row of rows) {
      const label = el("text", { x: PLOT_L - 12, y: row.y + 5, class: "wil-row-label", "text-anchor": "end", "font-size": 13 });
      label.append(el("tspan", { fill: LOOK[row.name].stroke, class: `wil-legend-${row.name.toLowerCase()}` }, row.name === "Treatment" ? "● " : "■ "));
      label.append(document.createTextNode(row.name));
      scoreSvg.append(label);
      scoreSvg.append(el("line", { x1: PLOT_L, y1: row.y, x2: width - 30, y2: row.y, class: "wil-row-line" }));
      row.data.forEach((d, index) => {
        const x = sx(d.value, width);
        const g = el("g", {
          class: `wil-dot wil-${row.name.toLowerCase()}`,
          transform: `translate(${x} ${row.y})`,
          role: "slider", tabindex: "0", "aria-label": `${row.name} cell ${index + 1} score`,
          "aria-valuemin": "0", "aria-valuemax": "10", "aria-valuenow": d.value,
          "data-id": d.id,
        });
        g.append(shape(row.name, row.name === "Treatment" ? 22 : 20));
        g.append(el("text", { y: 5, "text-anchor": "middle" }, d.value.toFixed(1)));
        g.addEventListener("pointerdown", event => {
          event.preventDefault();
          drag = { id: d.id, width };
          moveDot(event.clientX, scoreSvg.getBoundingClientRect());
        });
        g.addEventListener("keydown", event => {
          if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
          event.preventDefault();
          d.value = clamp(d.value + (event.key === "ArrowRight" ? 1 : -1) * (event.shiftKey ? 1 : 0.1));
          render(d.id);
        });
        scoreSvg.append(g);
      });
    }
    drawAxes(scoreSvg, width, 190, "Gene-set score");
  }

  function moveDot(clientX, rect) {
    if (!drag) return;
    const local = (clientX - rect.left) / rect.width * drag.width;
    const value = clamp((local - PLOT_L) / (drag.width - PLOT_L - 30) * 10);
    const point = [...control, ...treatment].find(d => d.id === drag.id);
    if (point && point.value !== value) { point.value = value; render(point.id); }
  }

  function drawRankChart(pooled) {
    const width = Math.max(560, rankSvg.parentElement.clientWidth || 800);
    const height = 190;
    rankSvg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    rankSvg.replaceChildren();
    const left = 70, right = width - 35, y = 80;
    rankSvg.append(el("line", { x1: left, y1: y, x2: right, y2: y, class: "wil-axis" }));
    const n = pooled.length;
    pooled.forEach((d, index) => {
      const x = n === 1 ? (left + right) / 2 : left + index / (n - 1) * (right - left);
      const isTreatment = d.group === "Treatment";
      const mark = el("g", { class: `wil-rank-mark wil-${d.group.toLowerCase()}`, transform: `translate(${x} ${y})` });
      mark.append(shape(d.group, isTreatment ? 17 : 15));
      mark.append(el("text", { y: 5, "text-anchor": "middle" }, d.value.toFixed(1)));
      rankSvg.append(mark);
      rankSvg.append(el("text", { x, y: 120, class: "wil-rank-label", "text-anchor": "middle" }, `rank ${format(d.rank)}`));
    });
    rankSvg.append(el("text", { x: left, y: 28, class: "wil-legend-control" }, "■ Control"));
    rankSvg.append(el("text", { x: left + 105, y: 28, class: "wil-legend-treatment" }, "● Treatment"));
    rankSvg.append(el("text", { x: (left + right) / 2, y: 158, class: "wil-axis-label", "text-anchor": "middle" }, "Lowest score → highest score"));
  }

  function format(value) { return Number.isInteger(value) ? String(value) : value.toFixed(1); }
  function render(focusId) {
    const pooled = ranked();
    const stats = statistics(pooled);
    drawScoreChart(pooled);
    drawRankChart(pooled);
    document.getElementById("w-rank-sum").textContent = format(stats.rankSum);
    document.getElementById("w-rank-list").textContent = `Treatment ranks: ${stats.treatmentRanks.map(format).join(" + ")}`;
    document.getElementById("w-u").textContent = format(stats.u);
    document.getElementById("w-u-range").textContent = `Range 0–${stats.n1 * stats.n2}; null midpoint ${format(stats.center)}`;
    document.getElementById("w-p").textContent = stats.p < 0.001 ? "< 0.001" : stats.p.toFixed(3);
    document.getElementById("w-verdict").textContent = stats.p < 0.05 ? "Unusual under random labels" : "Compatible with random labels";
    const probability = stats.u / (stats.n1 * stats.n2);
    let direction = "about equally likely to outrank";
    if (probability > 0.5) direction = "more likely to outrank";
    if (probability < 0.5) direction = "less likely to outrank";
    document.getElementById("w-interpretation").textContent = `A randomly chosen treatment cell is ${direction} a randomly chosen control cell (${Math.round(probability * 100)}% of pairwise comparisons, with ties counting as half).`;
    if (focusId) requestAnimationFrame(() => scoreSvg.querySelector(`[data-id="${focusId}"]`)?.focus());
  }

  document.getElementById("w-control-count").addEventListener("input", event => {
    const n = Number(event.target.value); resizeGroup(control, n, "C");
    document.getElementById("w-n-control").value = n; render();
  });
  document.getElementById("w-treatment-count").addEventListener("input", event => {
    const n = Number(event.target.value); resizeGroup(treatment, n, "T");
    document.getElementById("w-n-treatment").value = n; render();
  });
  document.querySelectorAll("[data-preset]").forEach(button => button.addEventListener("click", () => {
    const [c, t] = presets[button.dataset.preset];
    control.forEach((d, i) => { d.value = c[i]; });
    treatment.forEach((d, i) => { d.value = t[i]; });
    render();
  }));
  new ResizeObserver(() => render()).observe(scoreSvg.parentElement);
  render();
})();
