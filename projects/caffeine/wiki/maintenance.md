# Maintenance

## Wiki Rules

- Prefer updating an existing page over creating a near-duplicate.
- One page should own one durable concept, source, lab, or decision.
- Every synthesis page should link to at least two related pages.
- Every factual claim that comes from a local file should point back to a raw source page or source path.
- Use wikilinks for concept navigation, even though plain Markdown viewers will not resolve them automatically.
- Keep source summaries in `raw/`; keep interpretation in `concepts/` or `labs/`.

## Page Types

- `raw/*`: source summaries and provenance.
- `concepts/*`: reusable biology, data, and computational concepts.
- `labs/*`: experiment-specific context and outputs.
- root pages: entrypoints, maps, and maintenance rules.

## Ingest Checklist

When adding a new paper, dataset, or lab:

1. Create or update a `raw/` source page.
2. Update any concept pages touched by the source.
3. Add backlinks from the source page to the concepts.
4. Update [concept map](concept-map.md) if the new page changes the mental model.
5. Update [learning path](learning-path.md) only if the ramp-up sequence changes.

## Known Limitations

This initial wiki is compiled from the local project files plus a brief web search for the Karpathy wiki pattern. It does not independently verify every biological citation in `PROPOSAL.md`.
