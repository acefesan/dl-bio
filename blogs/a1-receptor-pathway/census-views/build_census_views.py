"""Build compact Census visual fragments from the local SOMA snapshot."""

from __future__ import annotations

import html
import json
import urllib.request
from collections import Counter
from pathlib import Path

import pyarrow.compute as pc
import tiledbsoma as soma


HERE = Path(__file__).resolve().parent
URI = "/mnt/bulk/dl_bio/cellxgene_census/2025-11-08/soma"
OBO_PATH = HERE / "cl.obo"
OBO_URL = "https://purl.obolibrary.org/obo/cl.obo"
ASSAYS = [
    "10x 3' v3", "10x 3' v2", "10x 5' v2", "10x 5' v1",
    "10x multiome", "10x 5' transcription profiling",
    "10x 3' transcription profiling", "10x gene expression flex",
    "10x 3' v1",
]


def census_counts():
    tissues = Counter()
    cell_types = Counter()
    assay_filter = " or ".join(f"assay == {label!r}" for label in ASSAYS)
    with soma.open(URI) as census:
        obs = census["census_data"]["homo_sapiens"].obs
        reader = obs.read(
            column_names=["tissue_general", "cell_type", "cell_type_ontology_term_id"],
            value_filter=f"is_primary_data == True and ({assay_filter})",
        )
        for table in reader:
            for row in pc.value_counts(table.column("tissue_general")).to_pylist():
                tissues[str(row["values"])] += row["counts"]
            grouped = table.group_by(
                ["cell_type_ontology_term_id", "cell_type"]
            ).aggregate([("cell_type", "count")])
            for row in grouped.to_pylist():
                key = (row["cell_type_ontology_term_id"], row["cell_type"])
                cell_types[key] += row["cell_type_count"]
    return tissues, cell_types


def read_obo():
    if not OBO_PATH.exists():
        urllib.request.urlretrieve(OBO_URL, OBO_PATH)
    terms = {}
    term = None
    for line in OBO_PATH.read_text().splitlines():
        if line == "[Term]":
            if term and "id" in term:
                terms[term["id"]] = term
            term = {"parents": []}
        elif line.startswith("["):
            if term and "id" in term:
                terms[term["id"]] = term
            term = None
        elif term is not None:
            if line.startswith("id: "):
                term["id"] = line[4:]
            elif line.startswith("name: "):
                term["name"] = line[6:]
            elif line.startswith("is_a: "):
                term["parents"].append(line[6:].split(" ! ")[0])
    if term and "id" in term:
        terms[term["id"]] = term
    return terms


def make_tissue_view(tissues):
    rows = sorted(tissues.items(), key=lambda kv: (-kv[1], kv[0]))
    total = sum(tissues.values())
    biggest = rows[0][1]
    entries = []
    for name, count in rows:
        count_short = f"{count / 1_000_000:.2f}m" if count >= 1_000_000 else f"{count / 1000:.0f}k" if count >= 1000 else str(count)
        width = max(1.5, 100 * count / biggest)
        entries.append(
            f'<div class="tissue" data-tooltip="{html.escape(name)}: {count:,} cells">'
            f'<div class="label"><span>{html.escape(name)}</span><span class="num">{count_short}</span></div>'
            f'<div class="track"><span style="width:{width:.2f}%"></span></div></div>'
        )
    output = f'''<div id="census-tissues-v1">
<h2>10x human cells by tissue</h2>
<div class="meta">{total:,} primary cells · {len(rows)} broad tissue labels · Census 2025-11-08</div>
<div class="tissue-grid">{''.join(entries)}</div>
<style>
#census-tissues-v1 {{ color:var(--foreground); }}
#census-tissues-v1 .meta {{ color:var(--muted-foreground); margin-bottom:1rem; }}
#census-tissues-v1 .tissue-grid {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:.5rem 1rem; }}
#census-tissues-v1 .tissue {{ min-width:0; }}
#census-tissues-v1 .label {{ display:flex; justify-content:space-between; align-items:baseline; gap:.3rem; font-size:12px; }}
#census-tissues-v1 .label span:first-child {{ min-width:0; overflow-wrap:anywhere; }}
#census-tissues-v1 .num {{ flex:none; color:var(--muted-foreground); font-variant-numeric:tabular-nums; }}
#census-tissues-v1 .track {{ height:4px; background:var(--muted); margin-top:2px; }}
#census-tissues-v1 .track span {{ display:block; height:100%; background:var(--viz-series-1); }}
@media(max-width:600px) {{ #census-tissues-v1 .tissue-grid {{ grid-template-columns:repeat(2,minmax(0,1fr)); }} }}
@media(max-width:420px) {{ #census-tissues-v1 .tissue-grid {{ grid-template-columns:1fr; }} }}
</style>
</div>'''
    (HERE / "census-10x-tissues.html").write_text(output)


def make_cell_tree(cell_types, terms):
    root = "CL:0000000"
    observed = {}
    labels = {}
    for (term_id, label), count in cell_types.items():
        observed[term_id] = observed.get(term_id, 0) + count
        labels[term_id] = label
    cache = {}

    def path_to_root(term_id, visited=frozenset()):
        if term_id == root:
            return [root]
        if term_id in cache:
            return cache[term_id]
        if term_id in visited or term_id not in terms:
            return [root, term_id]
        options = []
        for parent in terms[term_id]["parents"]:
            path = path_to_root(parent, visited | {term_id})
            options.append(path + [term_id])
        # One documented is-a path per term; prefer the more detailed path.
        best = sorted(options, key=lambda p: (-len(p), tuple(p)))[0] if options else [root, term_id]
        cache[term_id] = best
        return best

    parent_of = {}
    nodes = {root}
    for term_id in observed:
        path = path_to_root(term_id)
        for parent, child in zip(path, path[1:]):
            if child not in parent_of:
                parent_of[child] = parent
            nodes.add(parent)
            nodes.add(child)
    children = {term_id: [] for term_id in nodes}
    for child, parent in parent_of.items():
        children[parent].append(child)
    subtree_count = {}

    def n_observed(term_id):
        total = int(term_id in observed)
        for child in children.get(term_id, []):
            total += n_observed(child)
        subtree_count[term_id] = total
        return total

    n_observed(root)
    data = []
    for term_id in nodes:
        data.append({
            "id": term_id,
            "name": labels.get(term_id, terms.get(term_id, {}).get("name", term_id)),
            "parent": parent_of.get(term_id),
            "direct": observed.get(term_id, 0),
            "n": subtree_count.get(term_id, 0),
        })
    data.sort(key=lambda node: node["id"])
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    output = f'''<div id="census-cell-tree-v1">
<h2>10x cell-type hierarchy</h2>
<div class="meta">{len(observed)} observed labels · {sum(observed.values()):,} primary human cells · is-a classification, not a developmental lineage</div>
<label class="form-label" for="census-cell-search">Find a cell type</label>
<input class="form-control" id="census-cell-search" type="search" placeholder="e.g., hippocampal, astrocyte, T cell">
<div id="census-tree-status" class="meta" aria-live="polite"></div>
<div id="census-tree-nodes" role="tree" aria-label="Cell Ontology hierarchy"></div>
<style>
#census-cell-tree-v1 {{ color:var(--foreground); }}
#census-cell-tree-v1 .meta {{ color:var(--muted-foreground); margin-bottom:.65rem; }}
#census-cell-tree-v1 #census-tree-nodes {{ margin-top:.5rem; }}
#census-cell-tree-v1 details {{ margin-left:1rem; border-left:1px solid var(--border); padding-left:.7rem; }}
#census-cell-tree-v1 details:first-child {{ margin-left:0; border:0; padding-left:0; }}
#census-cell-tree-v1 summary, #census-cell-tree-v1 .leaf {{ padding:.22rem 0; }}
#census-cell-tree-v1 summary {{ cursor:pointer; }}
#census-cell-tree-v1 .leaf {{ margin-left:1rem; }}
#census-cell-tree-v1 .name {{ font-weight:500; }}
#census-cell-tree-v1 .count {{ color:var(--muted-foreground); font-variant-numeric:tabular-nums; margin-left:.45rem; }}
#census-cell-tree-v1 mark {{ color:var(--foreground); background:var(--accent); }}
</style>
<script>
(() => {{
  const root = document.getElementById('census-cell-tree-v1');
  const nodes = {payload};
  const byId = new Map(nodes.map(d => [d.id, d]));
  const children = new Map(nodes.map(d => [d.id, []]));
  for (const d of nodes) if (d.parent && children.has(d.parent)) children.get(d.parent).push(d.id);
  for (const list of children.values()) list.sort((a,b) => byId.get(a).name.localeCompare(byId.get(b).name));
  const input = root.querySelector('#census-cell-search');
  const status = root.querySelector('#census-tree-status');
  const target = root.querySelector('#census-tree-nodes');
  const saved = window.openai?.widgetState?.privateContent;
  if (saved && typeof saved.search === 'string') input.value = saved.search;
  function render() {{
    const query = input.value.trim().toLocaleLowerCase();
    const visible = new Set(); let hits = 0;
    if (query) {{
      for (const d of nodes) if (d.direct && d.name.toLocaleLowerCase().includes(query)) {{
        hits++; let p = d;
        while (p) {{ visible.add(p.id); p = p.parent ? byId.get(p.parent) : null; }}
      }}
    }}
    function add(id, depth) {{
      if (query && !visible.has(id)) return null;
      const d = byId.get(id), childIds = children.get(id).filter(x => !query || visible.has(x));
      const node = document.createElement(childIds.length ? 'details' : 'div');
      if (!childIds.length) node.className = 'leaf';
      if (childIds.length && (query || depth === 0)) node.open = true;
      const label = document.createElement(childIds.length ? 'summary' : 'span');
      const name = document.createElement('span'); name.className = 'name'; name.textContent = d.name;
      const count = document.createElement('span'); count.className = 'count';
      count.textContent = d.direct ? `${{d.direct.toLocaleString()}} cells` : `${{d.n}} types`;
      label.append(name, count); node.append(label);
      for (const child of childIds) {{ const el = add(child, depth + 1); if (el) node.append(el); }}
      return node;
    }}
    target.replaceChildren(add('CL:0000000', 0));
    status.textContent = query ? `${{hits}} matching observed cell types` : 'Expand branches or search to explore every observed type.';
  }}
  input.addEventListener('input', () => {{ render(); window.openai?.setWidgetState?.({{privateContent:{{search:input.value}}, modelContent:null}}).catch(() => {{}}); }});
  render();
}})();
</script>
</div>'''
    (HERE / "census-10x-cell-tree.html").write_text(output)
    print(f"tree: {len(observed)} observed terms, {len(nodes)} displayed nodes")


if __name__ == "__main__":
    tissues, cell_types = census_counts()
    make_tissue_view(tissues)
    make_cell_tree(cell_types, read_obo())
    print(f"tissues: {len(tissues)} labels, {sum(tissues.values())} cells")
