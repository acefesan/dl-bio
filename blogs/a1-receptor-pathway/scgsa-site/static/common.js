(function () {
  const root = document.documentElement.dataset.root || "";
  const here = document.documentElement.dataset.page || "";

  const sections = [
    { title: "Start", items: [
      ["index", "index.html", "Overview"],
      ["paper", "paper.html", "Read the paper"],
      ["debt", "debt.html", "Cognitive debt"],
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
    { title: "Hands on", items: [
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

  const layout = document.querySelector(".layout");
  if (layout) {
    layout.prepend(side);
    document.body.prepend(mob);
    const act = mob.querySelector("a.active");
    if (act) act.scrollIntoView({ inline: "center", block: "nearest" });
  }

  const KEY = "scgsa-theme";
  try {
    const saved = localStorage.getItem(KEY);
    if (saved) document.documentElement.dataset.theme = saved;
  } catch (e) { /* storage unavailable */ }

  side.querySelector(".theme-toggle").addEventListener("click", () => {
    const dark = window.matchMedia("(prefers-color-scheme: dark)").matches;
    const cur = document.documentElement.dataset.theme || (dark ? "dark" : "light");
    const next = cur === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    try { localStorage.setItem(KEY, next); } catch (e) { /* ignore */ }
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
