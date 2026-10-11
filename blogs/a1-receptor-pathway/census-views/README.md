# Standalone Census 10x views

These are the two standalone views created during the 10x-only exploration of
the local CELLxGENE Census 2025-11-08 SOMA snapshot:

- `census-10x-tissues-preview.html` — open directly in a browser.
- `census-10x-cell-tree-preview.html` — open directly in a browser.
- The versions without `-preview` are embeddable fragments.
- `build_census_views.py` regenerates all four files from the local SOMA store.

The tree is an ontology classification view, not a developmental-lineage tree.
The Cell Ontology is a graph; this visualization chooses one `is_a` parent per
observed term to make it displayable as a tree.

The builder downloads the official Cell Ontology `cl.obo` file when it is not
already present. That generated dependency is intentionally not checked in.

These snapshots are separate from the newer integrated Census pages in
`../scgsa-site/`, which include human and mouse views and are the active work.
