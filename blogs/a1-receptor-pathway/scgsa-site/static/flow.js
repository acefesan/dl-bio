// Small renderer for signalling ("activity flow") diagrams: boxes for
// molecules/processes, arrows for "activates" (→) and bar-ended lines for
// "inhibits" (⊣). A diagram is a spec { w, h, nodes, edges }; each node has
// id, label, x, y (centre), kind, and optional genes (shown on hover/tap and in
// the gene-set table). Usage: SCFLOW.draw(svgElement, spec).
(function () {
  const NS = "http://www.w3.org/2000/svg";
  const esc = t => String(t).replace(/[&<>"]/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));

  function box(n) {
    const lines = n.label.split("\n");
    const w = n.w || Math.max(...lines.map(l => l.length)) * 6.6 + 20;
    const h = 14 * lines.length + 12;
    return { x0: n.x - w / 2, y0: n.y - h / 2, w, h, lines };
  }
  // point where the segment from a's centre towards (tx, ty) leaves a's box
  function edgePoint(n, b, tx, ty) {
    const dx = tx - n.x, dy = ty - n.y;
    if (!dx && !dy) return [n.x, n.y];
    const sx = dx ? (b.w / 2) / Math.abs(dx) : Infinity, sy = dy ? (b.h / 2) / Math.abs(dy) : Infinity;
    const s = Math.min(sx, sy);
    return [n.x + dx * s, n.y + dy * s];
  }

  function draw(svg, spec) {
    const B = {};
    spec.nodes.forEach(n => { B[n.id] = box(n); });
    const N = Object.fromEntries(spec.nodes.map(n => [n.id, n]));
    let h = `<defs><marker id="${spec.id}-arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" class="fl-head"/></marker></defs>`;
    for (const lane of spec.lanes || []) {
      h += `<rect class="fl-lane" x="${lane.x}" y="${lane.y}" width="${lane.w}" height="${lane.h}" rx="12"/><text class="fl-lanel" x="${lane.x + 10}" y="${lane.y + 18}">${esc(lane.label)}</text>`;
    }
    for (const e of spec.edges) {
      const a = N[e.from], b = N[e.to];
      const [x1, y1] = edgePoint(a, B[e.from], b.x, b.y);
      let [x2, y2] = edgePoint(b, B[e.to], a.x, a.y);
      const len = Math.hypot(x2 - x1, y2 - y1) || 1, ux = (x2 - x1) / len, uy = (y2 - y1) / len;
      x2 -= ux * 3; y2 -= uy * 3;
      if (e.type === "inh") {
        h += `<line class="fl-edge inh" x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}"/>` +
          `<line class="fl-edge inh" x1="${x2 - uy * 8}" y1="${y2 + ux * 8}" x2="${x2 + uy * 8}" y2="${y2 - ux * 8}"/>`;
      } else {
        h += `<line class="fl-edge ${e.type || "act"}" x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" marker-end="url(#${spec.id}-arr)"/>`;
      }
      if (e.label) h += `<text class="fl-elabel" x="${(x1 + x2) / 2 + (e.dx || 6)}" y="${(y1 + y2) / 2 + (e.dy || 0)}">${esc(e.label)}</text>`;
    }
    for (const n of spec.nodes) {
      const b = B[n.id];
      const tip = n.genes ? `${n.label.replace(/\n/g, " ")} — genes: ${n.genes.join(", ")}` : n.label.replace(/\n/g, " ");
      h += `<g class="fl-node ${n.kind}" tabindex="0"><title>${esc(tip)}</title><rect x="${b.x0}" y="${b.y0}" width="${b.w}" height="${b.h}" rx="${n.kind === "effect" ? 4 : 9}"/>` +
        b.lines.map((l, i) => `<text x="${n.x}" y="${b.y0 + 17 + i * 14}" text-anchor="middle">${esc(l)}</text>`).join("") +
        (n.genes ? `<text class="fl-genes" x="${n.x}" y="${b.y0 + b.h + 11}" text-anchor="middle">${esc(n.genes.length > 3 ? n.genes.slice(0, 3).join(" ") + " …" : n.genes.join(" "))}</text>` : "") + `</g>`;
    }
    svg.setAttribute("viewBox", `0 0 ${spec.w} ${spec.h}`);
    svg.innerHTML = h;
  }

  // gene-set table for a diagram: node → genes, marking genes in our Census panel
  function geneTable(spec, panel) {
    const rows = spec.nodes.filter(n => n.genes).map(n => `<tr><td>${esc(n.label.replace(/\n/g, " "))}</td><td>${n.genes.map(g => panel && panel.has(g) ? `<b>${g}</b>` : g).join(", ")}</td></tr>`);
    return `<div class="table-wrap"><table class="fl-tab"><tr><th>step</th><th>genes</th></tr>${rows.join("")}</table></div>`;
  }
  const genesOf = spec => [...new Set(spec.nodes.flatMap(n => n.genes || []))];

  window.SCFLOW = { draw, geneTable, genesOf };
})();
