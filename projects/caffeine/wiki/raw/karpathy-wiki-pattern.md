# Karpathy-Style LLM Wiki Pattern

## Source

External web search, 2026-05-29. Useful references found:

- AgentWiki describes a shared knowledge base inspired by Karpathy's LLM Wiki concept, where raw sources are decomposed into atomic pages and cross-referenced.
- Charles Chen's wiki overview describes the pattern as a persistent, compounding Markdown knowledge base maintained by an LLM instead of rediscovering knowledge through stateless RAG on every query.
- Multiple public implementations describe the same practical conventions: one page per concept/entity/source, raw source summaries, synthesis pages, wikilinks, provenance, and a schema or maintenance file that tells future agents how to update the wiki.

## Compiled Pattern

A Karpathy-style wiki is not just a folder of notes. It is a small knowledge system with:

- Raw pages that summarize source documents and preserve provenance.
- Synthesis pages that compile concepts across sources.
- Dense cross-links so related pages are discoverable.
- A concept map or index that helps a reader enter the graph.
- Maintenance rules so future updates refine existing pages before creating duplicate concepts.

## How This Project Uses It

This caffeine wiki uses a lightweight version:

- `raw/` owns source-level summaries for local files and the style reference.
- `concepts/` owns reusable biological, data, and method concepts.
- `labs/` owns active experiment context.
- `README.md`, `concept-map.md`, and `learning-path.md` are the entry points.

Related pages: [maintenance](../maintenance.md), [concept map](../concept-map.md), [proposal source](proposal-source.md)
