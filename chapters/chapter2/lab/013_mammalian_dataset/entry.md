# 013 — Balanced Mammalian Dataset Design

**Date:** 2026-03-22
**Model:** N/A (dataset design)
**Status:** in-progress

## Motivation
The CAFA3 dataset has severe biases for orthology analysis (entries 011-012):
- **Species bias:** 8/22 taxa are vertebrates covering ~60% of proteins
- **HOG fragmentation:** 36,240 root HOGs but 60% are singletons — orphan genes
  because the dataset lacks close relatives for divergent species
- **HOG imbalance:** only HOG 801468 (2,236 proteins) is large enough for sub-HOG
  analysis; at level 5+, classes shrink to 5-18 proteins

A purpose-built mammalian dataset from OMA solves these problems:
- Mammalia-only scope → all HOGs are mammalian, no cross-kingdom noise
- 78 species in OMA → pick 30 spread across all major orders
- Uniform sampling → same N proteins per species, same N HOGs per species
- 56,136 HOGs defined at Mammalia level → rich orthology structure

## Dataset Design

### Species selection (30 species, 14 orders)

| Order             | N | Species |
|-------------------|---|---------|
| Monotremata       | 2 | Platypus, Echidna |
| Metatheria        | 2 | Opossum, Tasmanian devil |
| Afrotheria        | 2 | Elephant, Tenrec |
| Xenarthra         | 1 | Armadillo |
| Eulipotyphla      | 1 | Hedgehog |
| Chiroptera        | 2 | Horseshoe bat, Little brown bat |
| Carnivora         | 3 | Cat, Dog, Giant panda |
| Perissodactyla    | 1 | Horse |
| Artiodactyla      | 3 | Cow, Pig, Dolphin |
| Pholidota         | 1 | Pangolin |
| Primates          | 5 | Human, Macaque, Marmoset, Gibbon, Mouse lemur |
| Rodentia          | 4 | Mouse, Rat, Guinea pig, Naked mole-rat |
| Lagomorpha        | 1 | Rabbit |
| Scandentia        | 1 | Tree shrew |

Smallest species: Hedgehog (14,488 proteins). This sets the upper bound for
uniform sampling at ~14K proteins per species.

### Sampling strategy

1. **Uniform proteins per species:** Sample N proteins per species (N ≤ 14,488)
   - Only sample proteins that belong to a HOG at Mammalia level
   - Ensures every protein has a meaningful orthology label
2. **HOG-balanced sampling:** Within each species, sample to cover as many
   distinct HOGs as possible (1 protein per HOG first, then fill)
3. **Target:** ~12,000 proteins × 30 species = ~360,000 total

### Data to collect per protein

| Field | Source |
|-------|--------|
| omaid | OMA API `/api/genome/{code}/proteins/` |
| canonicalid | same |
| sequence | OMA API `/api/protein/{omaid}/` |
| oma_hog_id | same (full HOG path with sub-levels) |
| roothog_id | same |
| species_code | from genome |
| taxon_id | from genome |

### Data storage

```
dl_bio/assets/proteins/mammalia/
├── species.json             # 30 species metadata
├── sequences/
│   ├── HUMAN.fasta
│   ├── MOUSE.fasta
│   └── ...
├── hog_memberships.feather  # omaid → oma_hog_id, roothog_id
├── mammalia_annotations.feather  # full merged dataset
└── embeddings/
    ├── esm2_150m.feather
    ├── esm2_650m.feather
    └── esm2_3b.feather
```

## Download plan

### Step 1: Species metadata + protein lists (fast)
```bash
# For each species: GET /api/genome/{code}/proteins/?per_page=500
# Returns omaid, oma_hog_id, sequence_length — no sequence (fast)
# Paginate through all proteins, save to JSON/feather
# ~30 species × ~50 pages each = ~1500 API calls
```

### Step 2: HOG-balanced sampling
```python
# For each species:
#   1. Group proteins by roothog_id
#   2. Sample 1 protein per HOG (random)
#   3. If N_hogs < target, sample additional from largest HOGs
#   4. Result: uniform N proteins per species, maximum HOG diversity
```

### Step 3: Fetch sequences for sampled proteins
```bash
# POST /api/protein/bulk_retrieve/ (up to 1000 IDs per call)
# ~360,000 proteins / 1000 = 360 API calls
# Or: download oma-seqs.fa.gz (5.3GB) and extract by omaid
```

### Step 4: Compute ESM2 embeddings
```bash
# Use existing 01_compute_embeddings.py with --model flag
# RTX 5090 (32GB): 150M and 650M fit easily, 3B needs batch_size=4
```

## OMA API details

- **Base:** `https://omabrowser.org/api/`
- **Pagination:** `per_page` (max 500), `x-total-count` header
- **HOGs at Mammalia level:** 56,136
- **Bulk FASTA:** `https://omabrowser.org/All/oma-seqs.fa.gz` (5.3GB)
- **Rate limits:** undocumented, be polite (1 req/sec)

## Next steps
1. Write `fetch_mammalia.py` script implementing the download plan
2. Run protein list collection (Step 1)
3. Analyze HOG distribution across species before sampling
4. Execute sampling, sequence fetch, and embedding computation
