# Novelty Search Results

## Summary

Searched across 6 rounds (direct competitors, DMS tools, MPRA tools, in silico/AI tools, composable/declarative frameworks, and reviews). No tool was found that replicates PoolParty's combination of composable DAG-based design, multi-assay support, and design card metadata. Three tools and two reviews are worth further consideration.

---

## Papers to Investigate Further

### 1. Oligopool Calculator — Hossain et al. (2024)

**"Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements"**
Hossain A, Cetnar DP, LaFleur TL, McLellan JR, Salis HM.
ACS Synth Biol. 2024;13(12):4218-4232. DOI: 10.1021/acssynbio.4c00661

**What it does:** End-to-end Python suite for designing oligonucleotide pools (barcodes, primer binding sites) and analyzing barcoded sequencing reads. Can design 4 million barcodes in 1.2 hours, primer sites for 1 million oligos in 15 minutes.

**Relevance:** This is the most relevant new find. It designs components of oligo pools (barcodes, primers) computationally. However, from the abstract, it appears focused on the *physical construction* aspects of oligo pools (barcode uniqueness, primer compatibility, sequencing analysis) rather than on *variant content design* (what mutations/motifs to include). PoolParty addresses the latter; Oligopool Calculator addresses the former. They are likely complementary.

**Recommendation:** **Investigate.** If it does not design variant content (mutations, motifs, combinatorial logic), it belongs in a different category than PoolParty and may not need citation. But if it has any library content design features, it could warrant citation as a complementary tool. **I would need to read the paper to confirm.** Please download the PDF if you'd like me to do a full analysis.

### 2. Mutation Maker — Hiraga et al. (2021)

**"Mutation Maker, An Open Source Oligo Design Platform for Protein Engineering"**
Hiraga K et al. ACS Synth Biol. 2021;10(2):357-370. DOI: 10.1021/acssynbio.0c00542

**What it does:** Open-source software for designing mutagenic oligos for large-scale protein engineering. Supports multi-site random mutagenesis, directed mutagenesis, and de novo gene synthesis workflows. Uses optimization and constraint-satisfaction algorithms for primer design.

**Relevance:** Designs mutagenic oligos for protein engineering, which overlaps with PoolParty's DMS use case. However, Mutation Maker focuses on the *physical realization* (designing primers and oligos for wet-lab mutagenesis protocols) rather than on *specifying library content* as a structured object. It is already in the bib file (commented out as `article{Hiraga2021yg`).

**Recommendation:** **Borderline.** Could be cited alongside VaLiAnT as another assay-specific tool for DMS/protein engineering that focuses on oligo/primer design rather than library-level specification. Not essential, but would strengthen the "existing tools are format-specific" argument. **I would need to read the paper to confirm.** Please download if you'd like a full analysis.

### 3. SimDNA — Kundaje Lab (no formal publication)

**GitHub: https://github.com/kundajelab/simdna**

**What it does:** Python library for generating simulated regulatory DNA sequences with embedded motifs. Uses a composable architecture: BackgroundGenerator creates base sequences, Embedders place motifs at specified positions, and SequenceSetGenerator assembles the pipeline.

**Relevance:** The composable architecture (generators + embedders) has conceptual overlap with PoolParty's composable Operations. However, SimDNA is designed for generating *training/test data for ML models*, not for designing *experimental variant libraries*. It has no concept of design cards, lazy evaluation, DMS mutagenesis, barcodes, or library-level metadata. It also lacks a formal publication.

**Recommendation:** **Note but probably do not cite.** No formal publication, different purpose (ML training data vs. experimental libraries), and limited adoption. Worth being aware of if a reviewer asks about composable sequence generation frameworks.

---

## Reviews to Consider Citing

### 4. McEwen et al. (2025)

**"Multiplexed assays of variant effect for clinical variant interpretation"**
Nat Rev Genet. 2025. DOI: 10.1038/s41576-025-00870-x

**What it does:** Comprehensive recent review of MAVEs for clinical variant interpretation.

**Recommendation:** **Consider citing** in the Introduction alongside Starita2017 and Kinney2019 as a more recent MAVE reference. The manuscript already cites these two, but a 2025 review would show the field is still active and growing. Optional but could strengthen timeliness.

### 5. Chen et al. (2023)

**"Deep mutational scanning: A versatile tool in systematically mapping genotypes to phenotypes"**
Frontiers in Genetics. 2023. DOI: 10.3389/fgene.2023.1087267

**What it does:** Review of DMS methodology, applications, and computational analysis.

**Recommendation:** **Low priority.** The manuscript already has adequate DMS references (Fowler2014, Olson2014). Adding another review is not necessary unless you want to cite a more recent one.

---

## Confirmed Not Relevant

| Tool/Paper | Why not relevant |
|------------|-----------------|
| **Ledidi** (Schreiber et al. 2025b) | Optimizes individual sequences for desired ML model outputs; not a library design tool |
| **Oligo Designer Toolsuite** | Designs probes for spatial transcriptomics (MERFISH, etc.); different domain entirely |
| **OMEGA** (2025 bioRxiv) | Gene assembly from oligo pools; about physical construction, not variant design |
| **fastISM / Yuzu / TISM** | In silico saturation mutagenesis for model interpretation; not library design |
| **DYNAMCC_D / DC-Analyzer** | Degenerate codon optimization tools; same category as CodonGenie (already cited) |
| **Golden Mutagenesis** | Primer design for Golden Gate-based mutagenesis; wet-lab protocol tool |
| **DMS_PrimerDesignTool** | Simple GitHub primer design script; not a published tool |
| **Easy DNA** | Simple sequence manipulation utilities; no library concept |
| **Tiled-Region Exchange Mutagenesis** | Wet-lab method, not software |
| **inMOTIFin** (Ferenc 2025) | Regulatory sequence simulator for ML training; not experimental library design |

---

## Bottom Line

The search confirms that **no existing tool replicates PoolParty's combination of**:
- Composable, DAG-based library design
- Support across MPRA, DMS, and in silico experiments
- Design cards as structured metadata
- Lazy evaluation with on-demand sequence generation

The closest prior art remains **VaLiAnT** (already cited). The **Oligopool Calculator** is the most interesting new find but likely addresses a complementary problem (physical oligo pool construction) rather than variant content specification. **Mutation Maker** could optionally be added as another assay-specific tool example.
