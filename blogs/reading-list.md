# Reading list: blog post series

Things to read for the dl_bio blog posts (`blogs/`). One entry per item:
date added, link, what it is, and why it is here. Move an item into a post's
own folder when it becomes part of that post.

- **2026-10-20** — [Neural Networks Are Already Clean](https://already-clean.goodfire.com/)
  (Nick Cammarata, Goodfire; informal working-hypothesis essay, no date on the
  page). Argues that a network's own neurons are already mostly interpretable
  ("monosemantic"), and that most apparent polysemanticity comes from weak
  tools, poor architectures and "recursive generalization", not superposition.
  Studies vision models (InceptionV1, ConvNeXt, DINOv2) and language models
  (Qwen3-32B: about 90% of neuron poles with one main job in early audits)
  directly, with feature visualization, dataset examples, steering and
  attribution, and uses sparse autoencoders mainly as a comparison. Nothing on
  biology models.

  **Why it is here — an idea for the A1 / caffeine series:** probe the
  biology through a *model* instead of only through the data. Single-cell
  foundation models (for example Geneformer, scGPT, UCE) are trained on large
  cell corpora such as the CELLxGENE Census; if their neurons are as clean as
  this essay claims for vision and language models, some units could encode
  cell types, vascular zonation, or adenosine / cAMP programs, and could be
  read and steered. Open questions before we try it:
  - which single-cell models are open and small enough to run here, and
    whether interpretability work on them already exists (search first; do
    not assume);
  - a model trained on the Census inherits its problems (study batch effects,
    submitter labels, cells and nuclei mixed, uneven gene coverage), so an
    interpretable unit can encode a study or a chemistry, not biology;
  - a low-copy gene such as ADORA2A is mostly zeros in the training data, so
    a model may barely represent it;
  - the repo's chapter 2 already works with protein language models (ESM2),
    so the tools for reading neurons can be shared.
