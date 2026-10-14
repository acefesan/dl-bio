(function () {
  const root = document.documentElement.dataset.root || "";
  const here = document.documentElement.dataset.page || "";

  const sections = [
    { title: "Start", items: [
      ["index", "index.html", "Overview"],
      ["paper", "paper.html", "Read the paper"],
      ["debt", "debt.html", "Cognitive debt"],
      ["data", "data.html", "Which data is this?"],
    ]},
    { title: "Foundations", items: [
      ["basics", "basics.html", "Scoring from zero"],
      ["scenarios", "scenarios.html", "The four scenarios"],
      ["compare", "compare.html", "Side by side"],
    ]},
    { title: "Methods", items: [
      ["addmodulescore", "methods/addmodulescore.html", "AddModuleScore"],
      ["scse", "methods/scse.html", "SCSE"],
      ["aucell", "methods/aucell.html", "AUCell"],
      ["ucell", "methods/ucell.html", "UCell"],
      ["ssgsea", "methods/ssgsea.html", "ssGSEA"],
      ["jasmine", "methods/jasmine.html", "JASMINE"],
      ["scps", "methods/scps.html", "scPS"],
    ]},
    { title: "Census (SOMA)", items: [
      ["census", "census.html", "What is in the Census"],
      ["census-tree", "census-tree.html", "Cell type tree"],
      ["adora-census", "adora-census.html", "ADORA across the Census"],
      ["adora-spatial", "adora-spatial.html", "ADORA in tissue sections"],
      ["pathways", "pathways.html", "Gene sets and pathways"],
    ]},
    { title: "Hands on", items: [
      ["wilcoxon", "wilcoxon.html", "Wilcoxon playground"],
      ["playground", "playground.html", "Playground"],
      ["ours", "ours.html", "Our A1 data"],
    ]},
  ];

  const side = document.createElement("nav");
  side.className = "sidebar";
  let html = `<a class="brand" href="${root}index.html">scGSA from first principles</a>
    <div class="sub">Wang &amp; Thakar 2024, taken apart</div>`;
  for (const s of sections) {
    html += `<h4>${s.title}</h4>`;
    for (const [key, href, label] of s.items) {
      html += `<a class="nav${key === here ? " active" : ""}" href="${root}${href}">${label}</a>`;
    }
  }
  html += `<button class="theme-toggle" type="button">Toggle theme</button>`;
  side.innerHTML = html;

  const mob = document.createElement("nav");
  mob.className = "mobile-nav";
  mob.innerHTML = sections.flatMap(s => s.items).map(([key, href, label]) =>
    `<a class="${key === here ? "active" : ""}" href="${root}${href}">${label}</a>`).join("");

  // Data badge: <html data-source="geo|bulk|census|spatial ..."> names the
  // page's data source(s), linked to the matching card on data.html.
  const SRC = {
    geo: ["Paper data · GEO GSE164381", "the paper's own cells, rebuilt from GEO; not the Census"],
    bulk: ["Single-dataset files · HBCA / Tabula", "three CELLxGENE study files with their own labels (e.g. CA4)"],
    census: ["Census (SOMA) · 2025-11-08", "all studies in one harmonised store; standard labels; no batch correction"],
    spatial: ["Census spatial (SOMA) · 2025-11-08", "the Census's separate spatial collection, 540 sections"],
  };
  const srcs = (document.documentElement.dataset.source || "").split(/\s+/).filter(k => SRC[k]);
  if (srcs.length) {
    const kick = document.querySelector("article .kicker");
    const div = document.createElement("div");
    div.className = "data-badges";
    div.innerHTML = srcs.map(k => `<a class="data-badge ${k === "spatial" ? "census" : k}" href="${root}data.html#${k === "spatial" ? "census" : k}" title="${SRC[k][1]}">Data: ${SRC[k][0]}</a>`).join("");
    if (kick) kick.after(div);
  }

  const layout = document.querySelector(".layout");
  if (layout) {
    layout.prepend(side);
    document.body.prepend(mob);
    const act = mob.querySelector("a.active");
    if (act) act.scrollIntoView({ inline: "center", block: "nearest" });
  }

  // Theme: "auto" (default) follows the system via prefers-color-scheme in CSS;
  // "light"/"dark" are explicit overrides set from the button. Old key
  // "scgsa-theme" is ignored so a stale override can't pin the site.
  const KEY = "scgsa-theme-pref";
  const MODES = ["auto", "light", "dark"];
  const btn = side.querySelector(".theme-toggle");
  let mode = "auto";
  try {
    const saved = localStorage.getItem(KEY);
    if (MODES.includes(saved)) mode = saved;
    localStorage.removeItem("scgsa-theme");
  } catch (e) { /* storage unavailable */ }

  function applyTheme() {
    if (mode === "auto") delete document.documentElement.dataset.theme;
    else document.documentElement.dataset.theme = mode;
    const sys = window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
    btn.textContent = mode === "auto" ? `Theme: auto (${sys})` : `Theme: ${mode}`;
  }
  applyTheme();
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", applyTheme);

  btn.addEventListener("click", () => {
    mode = MODES[(MODES.indexOf(mode) + 1) % MODES.length];
    try { localStorage.setItem(KEY, mode); } catch (e) { /* ignore */ }
    applyTheme();
  });

  function renderMath() {
    if (window.renderMathInElement) {
      window.renderMathInElement(document.body, {
        delimiters: [
          { left: "$$", right: "$$", display: true },
          { left: "\\(", right: "\\)", display: false },
        ],
        throwOnError: false,
      });
    }
  }
  if (document.readyState === "complete") renderMath();
  else window.addEventListener("load", renderMath);
})();
