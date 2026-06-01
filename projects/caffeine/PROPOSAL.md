# Mapping caffeine's epigenomic fingerprint across every human cell type

**Caffeine rewires gene regulation far beyond the brain, and public datasets now make it possible to map those effects computationally across the entire body.** This report provides a complete project blueprint — from foundational knowledge through research design to a working analysis pipeline — for building a purely computational model of how caffeine alters the epigenomic landscape of every major human tissue and cell type. The project leverages more than a dozen public databases encompassing millions of single cells, hundreds of reference epigenomes, and dozens of caffeine-related GWAS loci. A critical finding that shapes the entire project: **very few epigenomic experiments have directly tested caffeine treatment in human cells**, meaning the most productive strategy is to integrate baseline epigenomic atlases with caffeine-responsive gene lists, genetic association data, and signaling pathway knowledge to *predict* caffeine's regulatory effects across untested cell types.

---

## Phase 1: Self-study reference — the biology you need to know

This section serves as a standalone reference document. Every technical concept is introduced and explained so that a reader with basic epigenetics and endocrine knowledge can follow the entire project.

### How caffeine works at the molecular level

Caffeine (1,3,7-trimethylxanthine) is a **non-selective adenosine receptor antagonist** — it blocks all four adenosine receptor subtypes (ADORA1, ADORA2A, ADORA2B, ADORA3) at concentrations achieved by normal coffee consumption (**1–50 µM** in plasma). Adenosine is a nucleoside that accumulates during wakefulness and signals tiredness; caffeine blocks this signal.

Each adenosine receptor couples to a different G-protein and triggers distinct intracellular cascades:

- **ADORA1** (A1 receptor): Couples to inhibitory Gi/o proteins. When adenosine binds, it *inhibits* adenylyl cyclase (the enzyme that makes cAMP, a key signaling molecule). Caffeine blocks this inhibition, effectively *increasing* cAMP. ADORA1 is the most abundant adenosine receptor in the brain (cortex, hippocampus, cerebellum), and is also highly expressed in **adipose tissue, heart (atria), kidney, and testis**.

- **ADORA2A** (A2A receptor): Couples to stimulatory Gs/Golf proteins. Normally *increases* cAMP when adenosine binds. Caffeine blocks this, which paradoxically *decreases* cAMP through A2A — but the net effect across the cell depends on which receptors dominate. ADORA2A is extremely concentrated on **striatal medium spiny neurons** (co-expressed with dopamine D2 receptors), and is also expressed on **T cells, macrophages, NK cells, endothelial cells, and platelets**. Under pathological conditions, A2A is robustly upregulated in astrocytes and microglia.

- **ADORA2B** (A2B receptor): A low-affinity receptor activated mainly under pathological conditions (when adenosine levels surge). Couples to Gs. Expressed on **intestinal epithelium, endothelial cells, cardiac fibroblasts, astrocytes, dendritic cells**, and is selectively induced in **liver** during ischemia.

- **ADORA3** (A3 receptor): Couples to inhibitory Gi/o. Expressed on **mast cells, neutrophils, macrophages, eosinophils, synovial tissue, lung**, and at lower levels in heart, brain, and testis.

Beyond adenosine receptor antagonism, caffeine has several secondary molecular targets that become relevant at higher concentrations:

**Phosphodiesterase (PDE) inhibition** (IC50 ~500 µM–1 mM). PDEs are enzymes that break down cAMP and cGMP, terminating their signaling. Caffeine inhibits PDEs nonselectively, prolonging cAMP/cGMP elevation. This is supraphysiological at normal coffee intake but pharmacologically relevant at high doses.

**Ryanodine receptor (RyR) activation** (threshold ~250 µM, robust at 5–20 mM). RyRs are calcium release channels on the endoplasmic/sarcoplasmic reticulum. Caffeine directly enhances calcium-induced calcium release (CICR), primarily in cardiac and skeletal muscle. RYR1 operates in skeletal muscle, RYR2 in cardiac muscle, RYR3 in brain.

**ATM/ATR kinase inhibition** (IC50 ~0.2 mM / ~1.1 mM). These are master regulators of the DNA damage checkpoint. Caffeine's inhibition of these kinases was discovered decades ago and is still used in cancer radiosensitization research.

**mTOR inhibition** (IC50 ~0.4 mM). mTOR (mechanistic target of rapamycin) is the master regulator of cell growth and autophagy (cellular self-cleaning). Caffeine inhibits mTOR signaling through the PI3K/Akt pathway, inducing autophagy — a potentially unifying mechanism for many of caffeine's pleiotropic effects in bone, skin, and neurons.

### The seven signaling cascades caffeine perturbs

Understanding these cascades is essential for interpreting epigenomic data, because each cascade terminates at specific **transcription factors** (TFs) that bind DNA and change gene expression.

**1. cAMP/PKA/CREB pathway.** This is the primary physiological cascade. Adenosine receptor blockade alters adenylyl cyclase activity → changes cAMP levels → activates or inhibits protein kinase A (PKA) → PKA phosphorylates CREB (cAMP response element-binding protein) at serine-133 → phospho-CREB recruits the coactivator CBP/p300 → activates transcription at CRE (cAMP response element) sequences throughout the genome. CREB target genes include **c-fos, BDNF, PER1, PER2** (circadian clock), tyrosine hydroxylase, and Bcl-2. There are an estimated **4,000+ CREB-controlled genes** genome-wide.

**2. PDE inhibition → sustained cAMP/cGMP.** By blocking cAMP and cGMP degradation, caffeine amplifies and prolongs signaling through PKA and protein kinase G (PKG). In cardiac tissue, PDE4D3 is part of the RyR2 macromolecular complex; its inhibition causes PKA-hyperphosphorylation of RyR2, creating "leaky" calcium channels.

**3. Ryanodine receptor → calcium release → CaMKII → HDAC export.** Caffeine-induced calcium release activates calmodulin-dependent protein kinase II (CaMKII), which phosphorylates HDAC5 (histone deacetylase 5), causing it to leave the nucleus. This directly connects caffeine to **epigenetic regulation**: HDAC5 removal from gene promoters allows histone hyperacetylation, opening chromatin and activating transcription. This was demonstrated at the GLUT4 promoter in muscle cells.

**4. MAPK/ERK pathway.** Caffeine has complex, concentration-dependent effects on the RAS→RAF→MEK→ERK cascade. At high concentrations it activates ERK1/2 in neuroblastoma cells; at moderate concentrations it inhibits ERK in osteosarcoma cells. ERK connects to CREB phosphorylation via MSK1 and can activate Nrf2.

**5. PI3K/Akt/mTOR pathway → autophagy.** Caffeine inhibits PI3K/Akt/mTOR/p70S6K signaling dose-dependently, inducing autophagy. This was demonstrated in neuroblastoma (SH-SY5Y) and osteosarcoma (HOS) cells at ≥5 mM caffeine — supraphysiological but with potential implications for chronic exposure.

**6. NF-κB pathway → inflammation.** Caffeine inhibits IKK/NF-κB signaling, reducing phosphorylation of NF-κB p65 at serine-536 and blocking its nuclear translocation. Downstream, this suppresses TNF-α, IL-6, IL-1β, COX-2, and the **NLRP3 inflammasome**. Anti-inflammatory effects are observed from ~100 µM to 5 mM depending on the model system.

**7. Nrf2/ARE pathway → antioxidant defense.** Caffeine activates Nrf2 (nuclear factor erythroid 2-related factor 2), which binds antioxidant response elements (AREs) to upregulate SOD, catalase, NQO1, and HMOX1. However, in PC12 neuronal cells, caffeine *inhibits* Nrf2 — demonstrating **cell-type-specific effects** that are a central theme of this project.

### Essential epigenomics vocabulary

**Chromatin accessibility** refers to how "open" or "closed" a region of DNA is. Open chromatin (accessible) means transcription factors can physically reach the DNA and activate nearby genes. Two assays measure this: **ATAC-seq** (Assay for Transposase-Accessible Chromatin with sequencing) uses a transposase enzyme that preferentially inserts into open regions, and **DNase-seq** uses the enzyme DNase I, which cuts open chromatin.

**ChIP-seq** (Chromatin Immunoprecipitation followed by sequencing) identifies where specific proteins — transcription factors or histone modifications — bind to DNA. Antibodies against the target protein pull down DNA fragments that were bound to it.

**Histone modifications** are chemical tags on histone proteins (which DNA wraps around). Key marks include: **H3K4me3** (trimethylation of lysine 4 on histone H3, marks active promoters), **H3K27ac** (acetylation at lysine 27, marks active enhancers and promoters), **H3K4me1** (marks poised/active enhancers), **H3K27me3** (marks repressed regions), and **H3K9me3** (marks constitutive heterochromatin).

**DNA methylation** is the addition of a methyl group to cytosine bases, predominantly at CpG dinucleotides. Methylation at gene promoters generally silences genes. **WGBS** (Whole-Genome Bisulfite Sequencing) maps methylation genome-wide. **EWAS** (Epigenome-Wide Association Studies) use methylation arrays (Illumina 450K or EPIC) to find CpG sites associated with traits — analogous to GWAS but for epigenetic marks.

**eQTL** (expression quantitative trait locus) is a genetic variant (SNP) that affects gene expression levels. GTEx cataloged eQTLs across 54 tissues, revealing tissue-specific genetic regulation.

**ChromHMM** is a computational method that combines multiple histone mark ChIP-seq datasets to segment the genome into chromatin states (e.g., active promoter, strong enhancer, repressed, quiescent). Roadmap Epigenomics provides 15-state ChromHMM models for 127 reference epigenomes.

**GRN** (Gene Regulatory Network) inference uses computational methods to reconstruct which transcription factors regulate which target genes, typically from single-cell RNA-seq or multi-omics data.

---

## Phase 2: Detailed research proposal — what data exists and what questions to ask

### The public data landscape for caffeine epigenomics

The following databases collectively provide the raw material for this project. A critical finding is that **direct caffeine-treatment epigenomic experiments are extremely rare** — the strategy must therefore integrate reference epigenomes with caffeine-responsive gene lists and genetic data.

**ENCODE** (encodeproject.org) provides over **25,000 experiments** across dozens of tissues and cell lines. For caffeine research, the most valuable data includes ATAC-seq on adrenal gland (ENCSR113MBR, ENCSR548KIL), liver (ENCSR373TDL), and multiple brain regions; ChIP-seq for >140 transcription factors and all major histone marks on key cell lines including **HepG2** (liver cancer, critical for CYP1A2 regulation), **GM12878** (B-lymphoblastoid), and **K562** (myeloid leukemia); and Hi-C data for 3D genome organization. ENCODE contains **no caffeine-treatment experiments** — it provides baseline epigenomic maps only.

**Roadmap Epigenomics** provides **127 consolidated reference epigenomes** (labeled E001–E129) with ChromHMM chromatin state maps. The caffeine-relevant epigenomes span 10 brain regions (E067–E074, E081–E082), liver (E066), four heart samples (E083, E095, E104, E105), fetal adrenal gland (E080), over a dozen immune cell types (E029–E062), and adipose tissue (E023, E063). Each epigenome includes core histone marks (H3K4me1, H3K4me3, H3K36me3, H3K27me3, H3K9me3), and many include DNA methylation and DNase-seq. All data is downloadable via FTP and viewable in the WashU Epigenome Browser.

**GEO** (Gene Expression Omnibus) is where the rare caffeine-treatment datasets live. The most important finding from GEO is the **Findley et al. (2019)** dataset containing paired RNA-seq and ATAC-seq from caffeine-treated human umbilical vein endothelial cells (HUVECs) across 17 donors — this is the single most valuable caffeine epigenomic dataset currently available. Additional GEO datasets include caffeine-treatment gene expression studies, adenosine receptor perturbation experiments, and theophylline treatment studies in airway cells. The major EWAS datasets (Karabegović et al. 2021, ~15,789 subjects; Schellhas et al. 2024, ~3,725 cord blood samples) used Illumina 450K/EPIC arrays deposited in GEO.

**Human Cell Atlas** (data.humancellatlas.org) provides single-cell resolution across **42 organs and >10.8 million cells** with unified cell-type labels. For this project, the most valuable resource is the **Zhang et al. (2021, Cell)** single-cell ATAC-seq atlas covering 30 adult tissue types with 615,998 nuclei and ~1.2 million candidate cis-regulatory elements (cCREs) across 222 cell types. The **Li et al. (2023, Science)** brain-specific snATAC-seq atlas covers 1.1 million cells across 42 brain regions with 107 cell types — essential for mapping ADORA1/ADORA2A regulatory elements in specific neuronal subtypes.

**GTEx** (gtexportal.org, version 8) provides bulk RNA-seq across **54 tissue types** from 838 donors (15,201 samples) plus comprehensive eQTL data. GTEx reveals that **CYP1A2 expression is virtually liver-exclusive** (~99%), ADORA1 is highest in brain and adipose, ADORA2A is enriched in caudate/putamen/nucleus accumbens and immune cells, ADORA2B is expressed in colon/blood vessels/lung, and ADORA3 is highest in whole blood and spleen. GTEx eQTL data for CYP1A2 in liver and ADORA2A in brain tissues are directly relevant to understanding caffeine sensitivity variation.

**CistromeDB** (cistrome.org/db) aggregates ~56,000 uniformly processed ChIP-seq, DNase-seq, and ATAC-seq profiles. Its toolkit allows querying "what factors regulate your gene?" for any gene — inputting ADORA2A or CYP1A2 returns all TFs with ChIP-seq evidence of binding near those genes. Key TF datasets include CREB1, AP-1 (JUN/FOS), NF-κB (RELA), Nrf2 (NFE2L2), AHR (regulates CYP1A2 directly), and HNF4A (liver-specific TF).

**JASPAR 2026** (jaspar.genereg.net, 10th release) and **HOCOMOCO v12** (hocomoco12.autosome.org) provide transcription factor binding motif databases. JASPAR contains curated position frequency matrices (PFMs) for all caffeine-relevant TFs: CREB1 (MA0018), AP-1/JUN/FOS (MA0099, MA0476), NF-κB/RELA (MA0107), Nrf2/NFE2L2 (MA0150), AHR::ARNT (MA0006), HNF4A (MA0114), and MEF2A/C (MA0052). Both databases are accessible programmatically through R (TFBSTools) and Python (pyJASPAR) packages.

**FANTOM5** (fantom.gsc.riken.jp/5) provides CAGE (Cap Analysis Gene Expression) profiles across >1,900 samples and >180 primary cell types, mapping ~201,000 promoters and ~65,000 enhancers. This provides promoter-level expression of ADORA genes and cell-type-specific enhancer usage that complements ENCODE/Roadmap.

**4D Nucleome** (4dnucleome.org) provides Hi-C, HiChIP, SPRITE, and Micro-C data for 3D chromatin organization. Critical for understanding which enhancers physically contact ADORA and CYP1A2 promoters — essential for interpreting non-coding GWAS variants.

**Additional resources** include IHEC (7,000+ reference datasets integrating 7 consortia), ChIP-Atlas (433,000+ experiments with data-mining tools), LINCS L1000/CMap (>1.3 million gene expression profiles for ~19,811 compounds including caffeine), PharmGKB (caffeine pharmacogenomics pathway), and DrugBank (DB00201, listing all known caffeine targets).

### What is already known about caffeine's epigenetic effects

**DNA methylation.** The landmark EWAS by Karabegović et al. (2021, PMID: 33990564) in 15,789 participants identified **11 CpG sites** reaching epigenome-wide significance for coffee consumption. These map to genes including **AHRR** (xenobiotic metabolism), **F2RL3**, **GFI1**, **HDAC4**, and **PHGDH** (linked to fatty liver disease). High coffee intake correlated with lower methylation at AHRR and F2RL3 — though these overlap with smoking-associated CpGs, creating a confounding concern. In cord blood, the Schellhas et al. (2024) meta-analysis found only **one CpG (cg19370043, near PRRX1)** associated with maternal caffeine intake, suggesting limited intrauterine epigenetic effects.

**Histone modifications.** The most detailed evidence comes from Mukwevho et al. (2008, PMID: 18198354), who showed that 10 mM caffeine induced **hyperacetylation of histone H3** at the MEF2 site on the GLUT4 promoter in C2C12 myotubes, through a CaMKII-dependent calcium signaling mechanism. Caffeine reduced nuclear HDAC5 content and increased MEF2A binding ~2.2-fold. In fetal rat adrenal tissue, prenatal caffeine enhanced DNMT1, DNMT3a, HDAC1, and HDAC2 expression while decreasing H3K9 and H3K14 acetylation at the steroidogenic factor-1 (SF-1) promoter. **No genome-wide ChIP-seq studies** comparing caffeine-treated vs. untreated cells have been published — only gene-specific ChIP.

**Chromatin accessibility.** The most important dataset is **Findley et al. (2019, PMID: 31492806)**, who performed paired RNA-seq and ATAC-seq in caffeine-treated HUVECs from 17 donors. They found significant changes in both gene expression and chromatin accessibility at 6 hours, with genes near differentially accessible regions 3.4–6.5× more likely to be differentially expressed. Motif enrichment revealed **INSM1, PLAG1, and ZNF423** as caffeine response factors, and SNPs in caffeine response TF footprints were enriched in GTEx eQTLs from artery tissues and colocalized with coronary artery disease risk. The follow-up study by **Boye et al. (2024, eLife)** used BiT-STARR-Seq (a massively parallel reporter assay) to test >43,000 genetic variants, finding **29 variants** that modulate gene expression specifically in response to caffeine, with **NFAT family transcription factors** (NFATC1, NFATC2, NFATC4) most enriched among differentially active targets.

**3D chromatin and non-coding RNA.** No Hi-C or chromatin conformation capture studies with caffeine treatment have been published — this is a major gap. For non-coding RNAs, caffeine downregulated ~8% of miRNAs in HeLa cells (notably miR-183-5p and miR-33a-5p) and reduced serum miR-9-3p, miR-15b-5p, miR-16-5p, and miR-222-3p in ethanol-exposed rat models.

### GWAS reveals caffeine biology extends far beyond the liver

Over a dozen GWAS have identified approximately **56 genomic loci** associated with caffeine/coffee intake. The established loci include AHR (rs4410790), CYP1A1/CYP1A2 (rs2472297), ABCG2 (rs1481012), POR (rs17685), and GCKR (rs1260326), which relate to caffeine pharmacokinetics. But the **novel pharmacodynamic loci** tell a more interesting story:

- **BDNF** (11p14) — brain-derived neurotrophic factor, linking caffeine to neuroplasticity
- **SLC6A4** (17q11.2) — serotonin transporter, connecting caffeine to serotonergic signaling
- **RORA** (15q22.2) — circadian rhythm regulation via retinoid-related orphan receptor
- **HCN2** (16q12.1) — cardiac/neuronal pacemaker channel
- **CACNA2D2** (12p13.31) — calcium channel subunit involved in neuronal signaling
- **MC4R** (18q21.32) — melanocortin 4 receptor, appetite/energy homeostasis
- **SORCS2** (10q21.1) — brain sorting receptor associated with ADHD and bipolar disorder
- **BRWD1** (21q22) — bromodomain WD repeat protein involved in chromatin remodeling

Tissue enrichment analyses of these GWAS loci show significant concentration in **central nervous system** tissues, indicating that caffeine's genetic architecture involves far more than hepatic metabolism. A comprehensive systematic review identified **59 Mendelian Randomization studies** on coffee/caffeine, revealing a key insight: MR analyses show divergent effects of **coffee intake** (often beneficial) versus **plasma caffeine levels** (potentially harmful), suggesting non-caffeine coffee components drive many observed health benefits.

### Cell-type-specific caffeine responsiveness can be computationally predicted

The central analytical challenge of this project is that most cell types have never been directly exposed to caffeine in an epigenomic experiment. The solution is a **computational imputation strategy** with four layers:

**Layer 1: Expression-based receptor mapping.** Using single-cell RNA-seq atlases (Human Cell Atlas, Tabula Sapiens, Human Protein Atlas), map which cell types express each adenosine receptor. This immediately identifies which cells have the molecular machinery to respond to caffeine. Key findings from existing data show that ADORA2A expression on CD8+ T cells means caffeine can enhance anti-tumor immunity, ADORA1 expression on adipocytes means caffeine promotes lipolysis, and ADORA2B expression on intestinal epithelium means caffeine affects gut chloride secretion.

**Layer 2: Chromatin accessibility at receptor loci.** Using single-cell ATAC-seq atlases, identify cell types with open chromatin at ADORA gene regulatory elements. A cell type with accessible chromatin at ADORA loci — even if current mRNA levels are low — is "primed" for receptor expression and could become responsive under stress or developmental transitions. Published methods validate this approach: pancreatic cancer organoids showed that ATAC-seq peaks predict drug sensitivity, and the **scE2G** tool (2024) classifies enhancer-gene regulation from scATAC-seq data trained on >10,000 CRISPR-validated element-gene pairs.

**Layer 3: Transcription factor activity inference.** Using chromVAR on scATAC-seq data, calculate per-cell TF motif deviation scores. Cells with high CREB motif accessibility are likely cAMP-responsive and therefore caffeine-responsive. TF footprinting tools like **TOBIAS** can identify occupied TF binding sites at nucleotide resolution from ATAC-seq data, revealing which TFs are physically bound at ADORA regulatory regions in each cell type.

**Layer 4: Gene regulatory network perturbation.** Using **CellOracle**, build cell-type-specific GRNs from paired scRNA-seq and scATAC-seq, then simulate what happens when adenosine receptor signaling is blocked. This predicts the full downstream transcriptomic response to caffeine in each cell type, even without direct experimental data.

### Unexpected cell types that respond to caffeine

Beyond the canonical CNS/cardiovascular/hepatic/endocrine systems, several surprising cell types show caffeine responsiveness:

**Oligodendrocyte precursor cells (OPCs)** express all four adenosine receptor subtypes. Adenosine promotes OPC differentiation and myelination — caffeine would antagonize this, with potential implications for white matter development.

**Osteoblasts and osteoclasts** show a biphasic dose response: low-moderate caffeine (3.125–12.5 µg/mL) inhibits osteoclastogenesis via MAPK/NF-κB and promotes bone formation, while high caffeine (50 µg/mL) reverses this, promoting bone resorption.

**Spermatocytes** express adenosine receptors (primarily A1R) within the seminiferous tubule epithelium. Caffeine enhances sperm motility in vitro through PDE inhibition but may increase DNA aneuploidy.

**Vascular smooth muscle cells** undergo osteogenic differentiation driven by ADORA2A → cAMP/CREB1/RUNX2 signaling. Caffeine's A2A antagonism could potentially attenuate vascular calcification — a finding from 2025 with direct clinical implications.

**Gut microbiome interactions** represent an emerging axis: high caffeine consumption associates with increased alpha diversity, higher Faecalibacterium and Roseburia (butyrate producers), and decreased Erysipelatoclostridium (linked to metabolic disease).

**Skin cells** show that caffeine protects against oxidative stress-induced senescence at just ~10 µM through A2A → SIRT3/AMPK → autophagy, and selectively induces apoptosis in melanoma cells while sparing normal melanocytes.

### Twenty specific, answerable computational research questions

These questions are ordered from beginner-accessible to publication-worthy and are all achievable with Python/R notebooks using public data.

**Beginner level (foundational mapping)**

1. *Which human cell types express adenosine receptors at highest levels?* Map ADORA1/2A/2B/3 across all cell types in Human Cell Atlas scRNA-seq data using Scanpy. Data: HCA, Tabula Sapiens.

2. *What are the tissue-specific eQTL effects of caffeine GWAS SNPs?* Perform colocalization analysis between published caffeine GWAS summary statistics and GTEx eQTLs using the coloc R package. Data: GWAS Catalog, GTEx v8.

3. *What gene sets are enriched among caffeine GWAS loci?* Run MAGMA or PASCAL gene-set enrichment on caffeine GWAS summary statistics against MSigDB, KEGG, and GO databases. Data: GWAS summary statistics, MSigDB.

4. *What is caffeine's transcriptomic signature in LINCS L1000, and which other drugs share it?* Query caffeine's gene expression profile in CMap/clue.io and identify drugs with correlated or anti-correlated signatures. Data: LINCS L1000 (GSE92742).

5. *How does CYP1A2 expression vary across liver cell subtypes?* Analyze single-cell RNA-seq data from human liver atlases to determine whether CYP1A2 is exclusively in hepatocytes or present in stellate cells, Kupffer cells, or cholangiocytes. Data: Human Liver Cell Atlas.

**Intermediate level (integrative analyses)**

6. *Do caffeine GWAS variants overlap enhancers active in unexpected cell types?* Test enrichment of caffeine GWAS SNPs in cell-type-specific ChromHMM enhancer states from Roadmap Epigenomics using GREGOR or GARFIELD. Data: GWAS Catalog, Roadmap Epigenomics 127 epigenomes.

7. *What WGCNA co-expression modules contain adenosine receptor genes across tissues?* Run weighted gene co-expression network analysis on GTEx data for brain, liver, heart, and immune tissues; identify hub genes within ADORA-containing modules. Data: GTEx v8 bulk RNA-seq.

8. *What transcription factor regulons are active in ADORA-high cell types?* Run pySCENIC on Human Cell Atlas scRNA-seq to identify TF regulons enriched specifically in cell types with high adenosine receptor expression. Data: HCA scRNA-seq.

9. *Can chromatin accessibility at ADORA loci predict receptor expression across cell types?* Correlate scATAC-seq peak scores at ADORA1/2A/2B/3 promoters and enhancers with matched scRNA-seq expression using multiome data. Data: 10x Multiome datasets from HCA.

10. *How does the NFAT transcription factor family respond to caffeine across vascular cell types?* Reanalyze the Findley (2019) and Boye (2024) HUVEC datasets, focusing on NFAT motif activity and connecting to vascular biology. Data: GEO (Findley/Boye datasets).

11. *Do caffeine-associated CpG sites from EWAS studies overlap with tissue-specific regulatory elements?* Map the 11 CpGs from Karabegović et al. onto ChromHMM states across all 127 Roadmap epigenomes to determine tissue-specific regulatory context. Data: EWAS results, Roadmap Epigenomics.

12. *Is there genetic pleiotropy between caffeine consumption and bone mineral density?* Perform cross-trait LD score regression and colocalization between caffeine GWAS and GEFOS bone mineral density GWAS. Data: GWAS summary statistics.

**Advanced level (publication-worthy)**

13. *Can CellOracle predict cell-type-specific transcriptomic responses to adenosine receptor blockade?* Build cell-type-specific GRNs from HCA multiome data, then simulate ADORA2A perturbation across all cell types to predict caffeine's downstream effects. Data: HCA scRNA-seq + scATAC-seq.

14. *What are the PheWAS signatures of caffeine-metabolizing enzyme variants?* Run phenome-wide association studies for rs2472297 (CYP1A2), rs4410790 (AHR), and novel GWAS loci across thousands of phenotypes in UK Biobank. Data: UK Biobank, FinnGen.

15. *Can multi-omics integration reveal how caffeine GWAS variants affect chromatin in immune cells?* Map caffeine GWAS SNPs onto immune cell-specific enhancers using DICE eQTLs, Blueprint Epigenome data, and scATAC-seq of PBMCs. Data: DICE, Blueprint, ENCODE.

16. *Can network pharmacology predict caffeine's effects on the tumor immune microenvironment?* Map caffeine targets (DrugBank) onto tumor microenvironment cell-type networks using TCGA single-cell data and CellChat. Data: TCGA, TISIDB, DrugBank.

17. *How does genetic variation in caffeine metabolism interact with gene expression programs in the developing brain?* Map caffeine GWAS variant effects onto developmental cell types using fetal brain scRNA-seq atlases and SCENIC regulon analysis. Data: Fetal brain atlases, BrainSpan.

18. *What is the 3D chromatin context of caffeine GWAS variants at ADORA loci?* Use 4D Nucleome Hi-C data and promoter-capture Hi-C to identify enhancer-promoter contacts at ADORA genes and determine which GWAS variants disrupt these contacts. Data: 4D Nucleome, ENCODE Hi-C.

19. *Can a multi-tissue caffeine response model identify causal tissue-mediating pathways for kidney function?* Combine UK Biobank GWAS, CKDGen data, and GTEx to perform tissue-partitioned Mendelian Randomization with mediation analysis for caffeine→kidney outcomes. Data: UK Biobank, CKDGen, GTEx.

20. *Can a graph neural network predict cell-type-specific caffeine responses from molecular features?* Train a model combining drug structure (molecular fingerprints), target profiles (DrugBank), and cell-type GRN features to predict response signatures in untested cell types. Data: LINCS L1000, scRNA-seq atlases, DrugBank.

---

## Phase 3: The computational analysis pipeline

### Pipeline architecture in seven phases

The pipeline follows a modular design where each phase produces defined outputs that feed into subsequent phases. The entire workflow can be orchestrated with **Snakemake** or **Nextflow** for reproducibility.

**Phase 1: Data collection and preprocessing.** Query the ENCODE REST API (`https://www.encodeproject.org/search/?type=Experiment&format=json`) for ATAC-seq and ChIP-seq across target tissues. Download GEO datasets with GEOparse, GTEx expression matrices via the portal API, and HCA data from the Data Portal. Use SRA Toolkit's `prefetch` and `fasterq-dump` for raw sequencing data. Organize into a standardized directory tree: `raw/{scRNA,scATAC,bulkRNA,ChIP,ATAC,WGBS,methylarray}/`. Validate file integrity with MD5 checksums. Expected storage: **5–10 TB** for a comprehensive analysis across 20+ cell types.

**Phase 2: Expression mapping.** For single-cell data, process with **Scanpy** (Python): filter cells (minimum 200 genes, <20% mitochondrial reads), remove doublets with Scrublet, normalize with scran, select highly variable genes, reduce dimensionality with PCA, correct batch effects with scVI (from scvi-tools), cluster with Leiden algorithm, and annotate cell types with **CellTypist**. Build a cell-type × gene expression matrix for all adenosine pathway genes (ADORA1/2A/2B/3, CYP1A2, PDE4B/4D, CREB1, RYR1/2/3, all downstream signaling genes). For bulk data, use STAR alignment → featureCounts → DESeq2. QC metrics: genes/cell, UMIs/cell, mitochondrial percentage, doublet rate, silhouette score for clustering quality.

**Phase 3: Chromatin accessibility analysis.** Process scATAC-seq with **ArchR** or **SnapATAC2**: import fragment files, QC filter (TSS enrichment >4, >1,000 fragments per cell), remove doublets, perform latent semantic indexing (LSI) for dimensionality reduction, cluster, and transfer cell-type labels from scRNA-seq. Call peaks per cell type with **MACS3** (`--nomodel --shift -100 --extsize 200`). Compute gene activity scores at ADORA loci. Run **Cicero** for co-accessibility analysis to link distal regulatory elements to ADORA promoters. For bulk ATAC-seq/DNase-seq, use MACS3 peak calling followed by DiffBind for differential accessibility across cell types. QC metrics: TSS enrichment (>6 ideal), fraction of reads in peaks (FRiP >0.3), fragment size distribution showing nucleosome periodicity.

**Phase 4: Epigenetic modification analysis.** For histone ChIP-seq, call peaks with MACS3 (narrow peaks for H3K4me3/H3K27ac, broad peaks for H3K27me3/H3K36me3). Run **ChromHMM** with 15-state model to segment genomes into chromatin states per cell type. Use DiffBind to compare histone marks at ADORA loci across cell types. Identify bivalent domains (H3K4me3 + H3K27me3) at adenosine receptor genes — these indicate poised genes that could be rapidly activated. For DNA methylation, process Illumina arrays with **minfi** (IDAT → SWAN normalization) and analyze with **RnBeads**; for bisulfite sequencing use **methylKit** to identify differentially methylated regions (≥3 CpGs, ≥25% methylation difference). Call super-enhancers near caffeine-responsive loci with the ROSE algorithm on H3K27ac signal.

**Phase 5: Transcription factor and motif analysis.** Run **HOMER** `findMotifsGenome.pl` on cell-type-specific ATAC-seq peaks near ADORA genes for both known and de novo motif enrichment. Use **FIMO** (from MEME Suite) to scan differentially accessible regions for specific TF motifs — particularly CREB (CRE: TGACGTCA), AP-1 (TGA[CG]TCA), NF-κB, Nrf2/ARE, and AHR::ARNT (XRE). Perform **TOBIAS** TF footprinting on bulk ATAC-seq: ATACorrect for Tn5 bias correction → ScoreBigwig for footprint scoring → BINDetect for genome-wide TF occupancy. For single-cell data, use **chromVAR** to compute per-cell TF motif deviation scores, identifying which cells have high CREB or AP-1 activity. Correlate TF expression (from Phase 2) with motif accessibility to identify activating versus repressing TFs at caffeine-responsive loci.

**Phase 6: Gene regulatory network inference.** Run **pySCENIC**: GRNBoost2 for TF-target co-expression → cisTarget for motif-based regulon pruning → AUCell for per-cell regulon activity scoring. If multiome data is available, use **SCENIC+** to infer enhancer-driven regulatory networks (eRegulons) that connect ATAC-seq peaks to target genes through TF binding. Build cell-type-specific GRNs with **CellOracle**: integrate scATAC-seq co-accessibility (Cicero) with motif scanning to define a base GRN, then use scRNA-seq to infer context-specific edge weights. The key analysis: simulate **in silico perturbation** of adenosine receptor signaling in each cell type to predict caffeine's downstream transcriptomic effects. For atlas-scale analysis, consider **LINGER**, which uses neural networks with bulk data regularization for improved accuracy.

**Phase 7: Integration and visualization.** Create unified cell-type × feature matrices linking expression, accessibility, methylation, histone marks, TF activity, and GRN topology. Run **GSEApy** prerank on differentially expressed genes per cell type against MSigDB Hallmark, KEGG, and Reactome collections. Build a caffeine-response regulatory network with NetworkX/igraph, overlay STRING protein-protein interactions, and run Louvain community detection to identify pathway modules. Generate genome browser views with **pyGenomeTracks** showing multi-layer ATAC-seq, ChIP-seq, methylation, ChromHMM states, TF footprints, and GRN links at ADORA loci across cell types. Create summary heatmaps with **ComplexHeatmap** showing cross-cell-type epigenomic profiles. Export results as UCSC Genome Browser track hubs for interactive exploration.

### Complete toolkit reference

The pipeline requires tools across ten functional categories:

| Category | Tools | Language |
|---|---|---|
| Data access | GEOparse, SRA Toolkit 3.1, ffq, encode-client | Python/CLI |
| Single-cell RNA-seq | Scanpy 1.12, Seurat 5.4, scvi-tools 1.4, CellTypist | Python/R |
| Single-cell ATAC-seq | ArchR 1.0, Signac 1.14, SnapATAC2 2.7, chromVAR 1.26 | R/Python |
| Bulk epigenomics | MACS3 3.0, DiffBind 3.14, DESeq2 1.44, ChromHMM 1.25 | Python/R/Java |
| DNA methylation | methylKit 1.30, RnBeads 2.22, minfi 1.50 | R |
| Motif analysis | MEME Suite 5.5, HOMER 4.11, TOBIAS 0.16 | CLI/Python |
| GRN inference | pySCENIC 0.12, CellOracle 0.18, SCENIC+ 1.0, LINGER | Python |
| Pathway analysis | GSEApy 1.1, clusterProfiler 4.12, g:Profiler | Python/R |
| Genomic operations | pybedtools 0.12, PyRanges 1.1, GenomicRanges 1.56 | Python/R |
| Visualization | pyGenomeTracks 3.9, ComplexHeatmap 2.20, IGV 2.18 | Python/R/Java |

### Hardware and runtime requirements

**Minimum configuration** for small-scale analysis (<50,000 cells): 16-core CPU, 64 GB RAM (128 GB recommended), NVIDIA GPU with ≥8 GB VRAM (for scvi-tools and SCENIC+), 2 TB NVMe SSD, Linux (Ubuntu 22.04+).

**Recommended cloud setup** for full pipeline (>200,000 cells, 20+ cell types): AWS r6i.8xlarge (32 vCPU, 256 GB RAM) for main analysis, g5.2xlarge (A10G GPU, 24 GB VRAM) for deep learning steps, with Snakemake or Nextflow orchestration using spot instances for parallel steps. Total storage: **5–10 TB**. Total runtime: **40–80 hours sequential, 20–40 hours parallelized** (approximately 1–2 days on recommended hardware). Key bottlenecks are TOBIAS footprinting (requires >50M reads per sample), pySCENIC GRNBoost2 (CPU-bound), and WGBS processing (I/O-bound).

Pin all tool versions in conda environment YAML files and use Snakemake's `--use-conda` with per-rule environments. Containerize with Docker/Singularity for full reproducibility.

---

## Where the field stands and where the gaps are

The existing literature reveals a strikingly uneven landscape. **The Findley/Boye endothelial cell studies** (2019, 2024) represent the only high-quality epigenomic datasets with direct caffeine treatment in human cells. The EWAS literature provides DNA methylation associations but exclusively from blood, creating a tissue-specificity gap. The GWAS field has identified ~56 loci but most functional annotation remains unexplored. Single-cell atlases now cover most human tissues but contain zero caffeine perturbation experiments.

**Ten critical gaps** define the opportunity space for this project:

1. No genome-wide ChIP-seq for histone modifications after caffeine treatment in any cell type
2. No Hi-C or 3D genome data with caffeine treatment
3. No single-cell multi-omics (scRNA-seq + scATAC-seq) under caffeine exposure
4. No systematic ATAC-seq time-course for caffeine response
5. No lncRNA profiling specific to caffeine
6. No integrative multi-omics analysis combining methylation + histones + accessibility + transcription for caffeine
7. Most signaling studies used supraphysiological caffeine concentrations (5–25 mM vs. 1–50 µM physiological)
8. EWAS data is blood-only; no liver, brain, or other tissue methylation data for caffeine
9. No computational gene regulatory network inference specific to caffeine has been published
10. The NFAT/calcium signaling axis identified by Boye et al. has not been explored in any non-endothelial cell type

## Conclusion

This project is feasible today because the data landscape has reached a tipping point: single-cell atlases covering millions of cells across dozens of tissues can now be combined with caffeine-specific GWAS, EWAS, and the Findley/Boye ATAC-seq datasets to build the first comprehensive computational model of caffeine's epigenomic impact across the human body. The most impactful contributions would come from three areas. First, **systematically mapping adenosine receptor expression and chromatin accessibility** across all cell types in existing single-cell atlases — work that has simply not been done despite the data being available. Second, **using CellOracle to simulate caffeine perturbation** in every major cell type, predicting downstream transcriptomic changes without new experiments. Third, **integrating GWAS variants with cell-type-specific enhancer maps** to discover which tissues and cell types mediate caffeine's genetic effects on complex traits. The mTOR/autophagy axis and NFAT/calcium signaling pathway emerging from recent work may prove to be unifying mechanisms that explain caffeine's surprisingly broad effects across cell types as different as neurons, osteoclasts, spermatocytes, and melanocytes. The computational tools now exist to test these hypotheses — the field simply needs someone to connect the datasets.
