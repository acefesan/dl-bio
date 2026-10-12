// Plain-language glossary. Any page that loads this gets a dotted underline on
// the first use of each technical term inside the containers listed in
// data-glossary on <html> (a CSS selector; default "main"); hover, focus or
// tap shows a one-line definition. New content (e.g. panels opened later) is
// picked up automatically.
(function () {
  const G = [
    [["GEO"], "Gene Expression Omnibus: the US public archive where researchers deposit sequencing data."],
    [["10x"], "A commercial droplet method that measures the RNA of thousands of single cells at once."],
    [["UMIs", "UMI"], "Unique molecular identifier: a random tag stuck on each captured RNA molecule, so copies of one molecule are counted once."],
    [["barcodes", "barcode"], "A DNA tag shared by every molecule from one droplet; it tells you which cell a molecule came from."],
    [["droplets", "droplet"], "A tiny oil-wrapped bubble that (ideally) holds one cell and one barcoded bead."],
    [["negative selection"], "Purifying one cell type by removing all the others: antibodies against non-B-cell markers pull those cells out with magnets, leaving untouched B cells."],
    [["PBMC"], "Peripheral blood mononuclear cells: the white blood cells with one round nucleus (T, B and NK cells, monocytes)."],
    [["quality control", "QC"], "Removing cells and genes that are broken, empty or uninformative before any analysis."],
    [["mitochondrial", "MT-"], "Genes carried by the mitochondria, the cell's power plants. When a large share of a cell's RNA is mitochondrial, the cell was probably dying."],
    [["doublets", "doublet"], "Two cells caught in one droplet; they look like a single cell with far too many genes."],
    [["LogNormalize", "log-normalised", "log-normalise", "normalisation", "normalise"], "Rescale each cell's counts to the same total (10,000) and take log(1 + x), so cells captured with more or fewer molecules can be compared."],
    [["variable genes"], "Genes whose levels differ a lot from cell to cell; they carry the information about what kind of cell each one is."],
    [["ScaleData"], "Seurat step that shifts and stretches each gene to average 0 and spread 1, so a very highly expressed gene doesn't drown out the others."],
    [["principal components", "PCA"], "Principal component analysis: compresses thousands of genes into a few summary axes, each a weighted mix of genes, that capture the biggest patterns of variation."],
    [["PCs"], "Principal components: the summary axes PCA produces."],
    [["nearest neighbours", "neighbour graph", "neighbours"], "A network linking each cell to the cells most similar to it."],
    [["Louvain"], "An algorithm that finds groups in a network whose members are more connected to each other than to everyone else."],
    [["resolution"], "Louvain's granularity knob: higher gives more, smaller clusters."],
    [["clusters", "cluster"], "A group of cells the algorithm found to be similar. It gets a number, not a name; naming it is a separate step."],
    [["UMAP"], "A method that squashes the data onto a 2-D picture so you can look at it. Good for seeing groups; distances on it are only roughly meaningful."],
    [["marker genes", "markers"], "Genes expressed mainly in one cell type; used to put names on clusters."],
    [["naive B"], "B cells that haven't yet met the thing their antibody recognises; they carry IgD and IgM antibody receptors."],
    [["memory B", "memory"], "B cells that have met their target before; often switched to IgG or IgA antibodies and positive for CD27."],
    [["plasma"], "Antibody-factory cells descended from B cells (marker: JCHAIN)."],
    [["Seurat"], "The standard R package for single-cell analysis."],
    [["Splatter"], "An R package that fits a statistical model to real single-cell data and simulates new, made-up cells from it."],
    [["RWSD"], "Real-world simulated data: real cells, randomly split in two, with a change added to one half."],
    [["SSD"], "Splatter-simulated data: made-up cells from a model fitted to the real ones."],
    [["zero fraction"], "The share of cells in which a gene has a count of 0."],
    [["detected", "detection"], "A gene is detected in a cell if at least one of its molecules was counted (count > 0)."],
    [["signal genes", "signal gene"], "The genes the simulation deliberately changes in treatment cells: the 'right answer' a method should find."],
    [["gene sets", "gene set"], "A list of genes scored together, such as the genes of one pathway."],
    [["ranks", "rank"], "A gene's position when one cell's genes are sorted from highest to lowest expression (1 = highest)."],
    [["Wilcoxon"], "A test of whether one group's values tend to sit higher than another's, using only their order, not their size."],
    [["p-values", "p-value"], "The chance of seeing a difference at least this large if the two groups were really the same. Small = surprising."],
    [["Benjamini–Hochberg", "adjusted p", "FDR"], "A correction for testing many gene sets at once: it keeps the expected share of false discoveries under the cut-off (here 5%)."],
    [["false positives", "false positive"], "Calling a difference that isn't there."],
    [["sensitivity"], "The share of real differences a method finds."],
    [["recovery rate", "recovery"], "The share of gene sets a method calls different: sensitivity when they truly differ, false-positive rate when they don't."],
    [["negative binomial", "negative-binomial"], "The usual statistical model for counts that are noisier than pure chance (Poisson) would give."],
    [["dispersion"], "How much extra spread counts have beyond pure chance."],
    [["varint"], "A compact way to store whole numbers: 1 byte for small ones, more only when needed."],
    [["scGSA"], "Single-cell gene-set analysis: turning a gene set into one score per cell."],
    [["pseudobulk"], "Adding up all cells of a group into one sample before analysis."],
  ];
  const lookup = new Map();
  const alts = [];
  G.forEach(([names, def], i) => names.forEach(n => { lookup.set(n.toLowerCase(), i); alts.push(n); }));
  alts.sort((a, b) => b.length - a.length);
  const escRe = s => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
  const RE = new RegExp(`(^|[^A-Za-z0-9-])(${alts.map(escRe).join("|")})(?![A-Za-z0-9])`, "i");
  const SKIP = "code, pre, select, option, button, svg, canvas, a, h1, h2, h3, h4, th, .term, .glossary-off, label, .kicker";

  function gloss(root) {
    if (!root || root.closest && root.closest(".glossary-off")) return;
    const seen = new Set(Array.from(root.querySelectorAll(".term")).map(t => t.dataset.term));
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
      acceptNode: n => (n.nodeValue.trim() && !(n.parentElement && n.parentElement.closest(SKIP))) ? NodeFilter.FILTER_ACCEPT : NodeFilter.FILTER_REJECT,
    });
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    for (let node of nodes) {
      let guard = 0;
      while (node && guard++ < 20) {
        const m = RE.exec(node.nodeValue);
        if (!m) break;
        const word = m[2], i = lookup.get(word.toLowerCase());
        const start = m.index + m[1].length;
        if (seen.has(String(i))) {
          // already explained in this container: continue after this match
          const rest = node.splitText(start + word.length);
          node = rest;
          continue;
        }
        seen.add(String(i));
        const after = node.splitText(start);
        const rest = after.splitText(word.length);
        const span = document.createElement("span");
        span.className = "term";
        span.tabIndex = 0;
        span.dataset.term = String(i);
        span.textContent = after.nodeValue;
        after.replaceWith(span);
        node = rest;
      }
    }
  }

  let tip;
  function show(el) {
    if (!tip) { tip = document.createElement("div"); tip.className = "gloss-tip"; tip.setAttribute("role", "tooltip"); document.body.appendChild(tip); }
    const def = G[Number(el.dataset.term)][1];
    tip.textContent = "";
    const t = document.createElement("b");
    t.textContent = el.textContent;
    tip.append(t, " — " + def);
    tip.hidden = false;
    const r = el.getBoundingClientRect(), w = Math.min(320, window.innerWidth - 16);
    tip.style.width = w + "px";
    let left = Math.max(8, Math.min(window.innerWidth - w - 8, r.left + r.width / 2 - w / 2));
    tip.style.left = left + window.scrollX + "px";
    const below = r.bottom + 8 + tip.offsetHeight < window.innerHeight;
    tip.style.top = (below ? r.bottom + 6 : r.top - tip.offsetHeight - 6) + window.scrollY + "px";
  }
  const hide = () => { if (tip) tip.hidden = true; };

  function init() {
    const sel = document.documentElement.dataset.glossary || "main";
    const roots = () => document.querySelectorAll(sel);
    let pending = false;
    const run = () => {
      pending = false;
      obs.disconnect();
      roots().forEach(r => (r.matches(".per-section") ? r.querySelectorAll(":scope > *") : [r]).forEach(gloss));
      roots().forEach(r => obs.observe(r, { childList: true, subtree: true }));
    };
    const obs = new MutationObserver(() => { if (!pending) { pending = true; requestAnimationFrame(run); } });
    run();
    document.addEventListener("mouseover", e => { const t = e.target.closest && e.target.closest(".term"); if (t) show(t); });
    document.addEventListener("mouseout", e => { if (e.target.closest && e.target.closest(".term")) hide(); });
    document.addEventListener("focusin", e => { if (e.target.classList && e.target.classList.contains("term")) show(e.target); });
    document.addEventListener("focusout", hide);
    document.addEventListener("click", e => { const t = e.target.closest && e.target.closest(".term"); if (t) show(t); else hide(); });
    window.addEventListener("scroll", hide, { passive: true });
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
