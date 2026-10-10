/* Tiny chart helpers (plain DOM, no dependencies). Colours come from CSS
   variables so light/dark follows the page. Any element with data-tip="..."
   gets a hover tooltip. */
window.V = (function () {
  const fmtN = n => n >= 1e9 ? (n / 1e9).toFixed(2) + "B" : n >= 1e6 ? (n / 1e6).toFixed(n >= 1e7 ? 1 : 2) + "M" : n >= 1e3 ? (n / 1e3).toFixed(n >= 1e4 ? 0 : 1) + "k" : String(Math.round(n));
  const fmtP = (x, d) => x == null || isNaN(x) ? "–" : (100 * x).toFixed(d == null ? (x >= 0.1 ? 0 : x >= 0.01 ? 1 : 2) : d) + "%";
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  function h(tag, attrs, kids) {
    const e = document.createElement(tag);
    for (const k in (attrs || {})) {
      if (k === "class") e.className = attrs[k]; else if (k === "html") e.innerHTML = attrs[k];
      else if (k === "style") e.style.cssText = attrs[k]; else e.setAttribute(k, attrs[k]);
    }
    (kids || []).forEach(c => e.append(c));
    return e;
  }

  // tooltip
  const tip = h("div", { class: "vtip" }); document.addEventListener("DOMContentLoaded", () => document.body.append(tip));
  document.addEventListener("mousemove", e => {
    const t = e.target.closest && e.target.closest("[data-tip]");
    if (!t) { tip.classList.remove("on"); return; }
    tip.textContent = t.getAttribute("data-tip"); tip.classList.add("on");
    const w = tip.offsetWidth, hh = tip.offsetHeight;
    tip.style.left = Math.min(e.clientX + 14, innerWidth - w - 8) + "px";
    tip.style.top = (e.clientY + 18 + hh > innerHeight ? e.clientY - hh - 10 : e.clientY + 18) + "px";
  });

  // segmented control: opts = [[value,label],...]; calls cb(value)
  function seg(label, opts, init, cb) {
    const g = h("span", { class: "grp" }); let cur = init;
    const bs = opts.map(([v, l]) => { const b = h("button", { type: "button", "aria-pressed": v === cur }); b.textContent = l;
      b.onclick = () => { cur = v; bs.forEach((x, i) => x.setAttribute("aria-pressed", opts[i][0] === cur)); cb(cur); }; g.append(b); return b; });
    const wrap = h("span", {}, [h("span", { class: "lab", html: esc(label) }), g]);
    return wrap;
  }

  // horizontal bars. rows: [{label, segs:[{v,color,tip}], text, tip}]
  function hbars(el, rows, opts) {
    opts = opts || {}; const max = opts.max || Math.max(...rows.map(r => r.segs.reduce((a, s) => a + s.v, 0)), 1e-9);
    const g = h("div", { class: "hb", style: opts.labelW ? `grid-template-columns:${opts.labelW} minmax(0,1fr) 9.6rem` : "" });
    rows.forEach(r => {
      const tr = h("div", { class: "row" });
      const t = h("div", { class: "t" });
      let left = 0;
      r.segs.forEach(s => { const w = 100 * s.v / max;
        t.append(h("div", { class: "b", style: `left:${left}%;width:${w}%;background:${s.color}`, "data-tip": s.tip || r.tip || "" })); left += w; });
      tr.append(h("div", { class: "l", title: r.label, "data-tip": r.tip || "", html: esc(r.label) }), t, h("div", { class: "v", html: esc(r.text) }));
      g.append(tr);
    });
    el.replaceChildren(g);
  }

  // sequential colour: t in [0,1] -> css colour (one hue, light -> dark via mixing)
  const seqColor = t => `color-mix(in oklab, var(--accent) ${Math.round(6 + 94 * Math.max(0, Math.min(1, t)))}%, var(--panel))`;
  const inkOn = t => t > 0.62 ? "var(--bg)" : "var(--ink)";

  // heatmap. cells[r][c] = {t (0-1) or null, text, tip}
  function heat(el, rowLabels, colLabels, cells, opts) {
    opts = opts || {}; const colW = opts.colW || 56;
    const g = h("div", { class: "heat", style: `grid-template-columns: ${opts.rowW || 190}px repeat(${colLabels.length}, minmax(${colW}px, 1fr))` });
    g.append(h("div"));
    colLabels.forEach(c => g.append(h("div", { class: "cl" + (opts.rot ? " rot" : ""), html: esc(c), title: c })));
    rowLabels.forEach((rl, r) => {
      g.append(h("div", { class: "rl", title: rl, html: esc(rl) }));
      colLabels.forEach((_, c) => { const x = cells[r][c];
        if (!x || x.t == null) g.append(h("div", { class: "hc na", "data-tip": x && x.tip || "not measured / no cells", html: x && x.text || "" }));
        else g.append(h("div", { class: "hc", style: `background:${seqColor(x.t)};color:${inkOn(x.t)}`, "data-tip": x.tip || "", html: esc(x.text || "") }));
      });
    });
    el.replaceChildren(h("div", { class: "heat-wrap" }, [g]));
  }

  // vertical-column histogram. bars: [{label, v, tip}]
  function hist(el, bars) {
    const max = Math.max(...bars.map(b => b.v), 1e-9);
    const cols = h("div", { class: "hist" }), xs = h("div", { class: "hist-x" });
    bars.forEach(b => { const c = h("div", { class: "col", "data-tip": b.tip || "" }); c.append(h("div", { class: "b", style: `height:${100 * b.v / max}%` })); cols.append(c);
      xs.append(h("span", { html: esc(b.label) })); });
    el.replaceChildren(cols, xs);
  }

  function load(url) { return fetch(url).then(r => { if (!r.ok) throw new Error(url + " " + r.status); return r.json(); }); }
  function legend(items) { return h("div", { class: "legend" }, items.map(([c, l]) => h("span", { html: `<i style="background:${c}"></i>${esc(l)}` }))); }
  function ready(fn) { if (document.readyState !== "loading") fn(); else document.addEventListener("DOMContentLoaded", fn); }
  return { fmtN, fmtP, esc, h, seg, hbars, heat, hist, load, legend, seqColor, ready };
})();
