// Signalling diagrams for pathways.html section 4 (drawn by static/flow.js).
// kinds: ligand, drug, receptor, gprotein, enzyme, channel, messenger, process, effect.
// Each node's genes double as the diagram's gene set (table under the figure).
(function () {
  const A1 = {
    id: "a1", w: 720, h: 480,
    nodes: [
      { id: "ado", label: "adenosine\n(builds up while awake)", x: 180, y: 32, kind: "ligand" },
      { id: "caf", label: "caffeine", x: 470, y: 32, kind: "drug" },
      { id: "r", label: "A1 receptor", x: 180, y: 112, kind: "receptor", genes: ["ADORA1"] },
      { id: "g", label: "Gi/o protein", x: 180, y: 192, kind: "gprotein", genes: ["GNAI1", "GNAI2", "GNAI3", "GNAO1", "GNB1", "GNG2"] },
      { id: "ac", label: "adenylyl cyclase", x: 100, y: 285, kind: "enzyme", genes: ["ADCY1", "ADCY5", "ADCY8"] },
      { id: "girk", label: "GIRK K⁺ channel", x: 355, y: 285, kind: "channel", genes: ["KCNJ3", "KCNJ6"] },
      { id: "ca", label: "N / P/Q Ca²⁺ channels", x: 590, y: 285, kind: "channel", genes: ["CACNA1A", "CACNA1B"] },
      { id: "camp", label: "less cAMP", x: 100, y: 370, kind: "messenger" },
      { id: "hyp", label: "K⁺ leaves:\ncell hyperpolarised", x: 355, y: 370, kind: "process" },
      { id: "rel", label: "less Ca²⁺ enters:\nless transmitter released", x: 590, y: 370, kind: "process" },
      { id: "out", label: "wake-promoting neurons fire less → sleep pressure", x: 355, y: 450, kind: "effect" },
    ],
    edges: [
      { from: "ado", to: "r" }, { from: "caf", to: "r", type: "inh", label: "blocks" },
      { from: "r", to: "g" }, { from: "g", to: "ac", type: "inh" },
      { from: "g", to: "girk", label: "Gβγ opens" }, { from: "g", to: "ca", type: "inh", label: "Gβγ" },
      { from: "ac", to: "camp" }, { from: "girk", to: "hyp" }, { from: "ca", to: "rel" },
      { from: "camp", to: "out" }, { from: "hyp", to: "out" }, { from: "rel", to: "out" },
    ],
  };
  const A2A = {
    id: "a2a", w: 720, h: 480,
    nodes: [
      { id: "ado", label: "adenosine", x: 160, y: 32, kind: "ligand" },
      { id: "caf", label: "caffeine", x: 380, y: 32, kind: "drug" },
      { id: "dop", label: "dopamine", x: 600, y: 32, kind: "ligand" },
      { id: "r", label: "A2A receptor", x: 160, y: 112, kind: "receptor", genes: ["ADORA2A"] },
      { id: "d2", label: "D2 dopamine receptor", x: 560, y: 112, kind: "receptor", genes: ["DRD2"] },
      { id: "g", label: "Golf protein", x: 160, y: 192, kind: "gprotein", genes: ["GNAL", "GNB1", "GNG7"] },
      { id: "ac", label: "adenylyl cyclase 5", x: 160, y: 270, kind: "enzyme", genes: ["ADCY5"] },
      { id: "camp", label: "more cAMP", x: 160, y: 345, kind: "messenger" },
      { id: "pka", label: "PKA", x: 370, y: 345, kind: "enzyme", genes: ["PRKACA", "PRKACB", "PRKAR1B", "PRKAR2B"] },
      { id: "darpp", label: "DARPP-32", x: 580, y: 345, kind: "enzyme", genes: ["PPP1R1B"] },
      { id: "out", label: "indirect-pathway neuron more active → its output damps arousal → sleepier", x: 370, y: 445, kind: "effect" },
    ],
    edges: [
      { from: "ado", to: "r" }, { from: "caf", to: "r", type: "inh", label: "blocks" }, { from: "dop", to: "d2" },
      { from: "r", to: "d2", type: "inh", label: "opposes", dy: -6 }, { from: "d2", to: "ac", type: "inh", label: "via Gi" },
      { from: "r", to: "g" }, { from: "g", to: "ac" }, { from: "ac", to: "camp" }, { from: "camp", to: "pka" },
      { from: "pka", to: "darpp", label: "phosphorylates", dx: -40, dy: -10 }, { from: "pka", to: "out" }, { from: "darpp", to: "out" },
    ],
  };
  const VASO = {
    id: "vaso", w: 720, h: 610,
    lanes: [
      { x: 6, y: 72, w: 300, h: 270, label: "endothelial cell (vessel lining)" },
      { x: 316, y: 72, w: 398, h: 470, label: "smooth muscle cell (vessel wall)" },
    ],
    nodes: [
      { id: "ado", label: "adenosine (more when tissue\nis busy or short of O₂)", x: 300, y: 32, kind: "ligand" },
      { id: "caf", label: "caffeine", x: 600, y: 32, kind: "drug" },
      { id: "re", label: "A2A / A2B receptor", x: 156, y: 125, kind: "receptor", genes: ["ADORA2A", "ADORA2B"] },
      { id: "enos", label: "eNOS", x: 156, y: 205, kind: "enzyme", genes: ["NOS3"] },
      { id: "no", label: "nitric oxide (NO)", x: 156, y: 300, kind: "messenger" },
      { id: "rs", label: "A2A / A2B receptor", x: 430, y: 125, kind: "receptor", genes: ["ADORA2A", "ADORA2B"] },
      { id: "gs", label: "Gs protein", x: 430, y: 190, kind: "gprotein", genes: ["GNAS"] },
      { id: "ac", label: "adenylyl cyclase", x: 430, y: 255, kind: "enzyme", genes: ["ADCY3", "ADCY5", "ADCY6"] },
      { id: "camp", label: "more cAMP", x: 430, y: 340, kind: "messenger" },
      { id: "pka", label: "PKA", x: 430, y: 400, kind: "enzyme", genes: ["PRKACA"] },
      { id: "sgc", label: "soluble guanylyl cyclase", x: 615, y: 300, kind: "enzyme", genes: ["GUCY1A1", "GUCY1B1"] },
      { id: "cgmp", label: "more cGMP", x: 615, y: 360, kind: "messenger" },
      { id: "pkg", label: "PKG", x: 615, y: 418, kind: "enzyme", genes: ["PRKG1"] },
      { id: "k", label: "K⁺ channels open\n(KATP, BK), so Ca²⁺\nchannels close", x: 410, y: 485, kind: "channel", genes: ["KCNJ8", "ABCC9", "KCNMA1", "CACNA1C"] },
      { id: "mlc", label: "less myosin\nlight-chain activity\n(MLCK off, MLCP on)", x: 620, y: 485, kind: "process", genes: ["MYLK", "PPP1R12A"] },
      { id: "out", label: "smooth muscle relaxes → vessel widens (vasodilation) → more blood flow", x: 380, y: 580, kind: "effect" },
    ],
    edges: [
      { from: "ado", to: "re" }, { from: "ado", to: "rs" },
      { from: "caf", to: "rs", type: "inh" }, { from: "caf", to: "re", type: "inh", label: "blocks both", dx: 40, dy: -16 },
      { from: "re", to: "enos" }, { from: "enos", to: "no" }, { from: "no", to: "sgc", label: "diffuses into the muscle", dx: -150, dy: -7 },
      { from: "rs", to: "gs" }, { from: "gs", to: "ac" }, { from: "ac", to: "camp" }, { from: "camp", to: "pka" },
      { from: "sgc", to: "cgmp" }, { from: "cgmp", to: "pkg" },
      { from: "pka", to: "k" }, { from: "pka", to: "mlc" }, { from: "pkg", to: "mlc" }, { from: "pkg", to: "k" },
      { from: "k", to: "out" }, { from: "mlc", to: "out" },
    ],
  };

  async function init() {
    if (!window.SCFLOW) return;
    let panel = new Set();
    try { panel = new Set(Object.keys(await fetch("data/adenosine_genes.json").then(r => r.json()))); } catch (e) { /* bolding is optional */ }
    for (const spec of [A1, A2A, VASO]) {
      const svg = document.getElementById("fl-" + spec.id);
      if (!svg) continue;
      SCFLOW.draw(svg, spec);
      const t = document.getElementById("gt-" + spec.id);
      if (t) {
        const genes = SCFLOW.genesOf(spec);
        const notIn = genes.filter(g => !panel.has(g));
        t.innerHTML = `<details><summary>Gene set for this diagram (${genes.length} genes)</summary>${SCFLOW.geneTable(spec, panel)}` +
          `<p class="cp-note">As one list: <code>${genes.join(" ")}</code>${notIn.length ? `<br>Not yet in our Census panel (would need a new scan): ${notIn.join(", ")}.` : ""}</p></details>`;
      }
    }
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
