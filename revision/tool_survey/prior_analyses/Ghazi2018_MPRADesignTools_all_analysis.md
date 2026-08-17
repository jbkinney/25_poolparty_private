# Ghazi et al. (2018) — Design tools for MPRA experiments

## BibTeX Entry

```bibtex
@article{Ghazi2018aa,
  year     = {2018},
  title    = {Design tools for {MPRA} experiments},
  author   = {Ghazi, Andrew R and Chen, Edward S and Henke, David M and Madan, Namrata and Edelstein, Leonard C and Shaw, Chad A},
  journal  = {Bioinformatics},
  issn     = {1367-4803},
  doi      = {10.1093/bioinformatics/bty150},
  pmid     = {30052913},
  pmcid    = {{PMC}6454564},
  pages    = {2682--2683},
  number   = {15},
  volume   = {34},
  mynotes  = {Web-based tool and R package for interactive MPRA experimental design, focused on statistical power analysis and automated sequence generation from VCF files for SNP-focused MPRAs. Should be cited in the Introduction as an example of an assay-specific tool that automates library design for a particular experimental format (SNP-focused MPRAs) but is not easily adapted to novel designs such as DMS libraries, motif grammar studies, or in silico experiments.}
}
```

## 1. Summary

MPRA Design Tools is a web application (built in R/Shiny) and companion R package for interactive design of MPRA experiments. Its two main contributions are: (1) **statistical power analysis** — users can adjust parameters (barcodes per allele, activity variance) to estimate the power to detect transcriptional shifts of a given effect size, calibrated using real MPRA data from Tewhey et al. (2016) and Ulirsch et al. (2016); (2) **automated sequence generation** — users upload a VCF file of variants of interest and the tool generates barcoded MPRA construct sequences with appropriate genomic context from the hg38 reference genome, including restriction sites and PCR primers following the Melnikov et al. (2012) oligo layout. The tool explicitly cites and differentiates itself from MPRAnator, noting improved ease-of-use, interactivity, and genomic context acquisition from the reference genome rather than requiring user input.

## 2. Relevance to PoolParty Manuscript

MPRA Design Tools and PoolParty both generate sequences for MPRA experiments, but they address different aspects of the design problem:

- **MPRA Design Tools** focuses narrowly on **SNP-centric MPRAs**: take a VCF of variants, pull genomic context, assign barcodes, and output construct sequences. Its distinguishing feature is the power analysis module.
- **PoolParty** provides a general-purpose framework for designing arbitrary variant libraries (MPRA, DMS, in silico, etc.) with composable Operations, design cards, and lazy evaluation.

Key differences:

| Feature | MPRA Design Tools | PoolParty |
|---------|------------------|-----------|
| Scope | SNP-focused MPRAs | General-purpose library design |
| Input | VCF files + reference genome | Arbitrary sequences + composable Operations |
| Power analysis | Yes (core feature) | No (out of scope) |
| Motif grammar designs | No | Yes (insert_multiscan, etc.) |
| DMS support | No | Yes (mutagenize_orf, etc.) |
| Design cards | No | Core feature |
| Architecture | R/Shiny web app | Python DAG-based library |
| Extensibility | Limited | User-defined Operations |

## 3. Citation Recommendations

**Should be cited** in the Introduction (paragraph 3) alongside MPRAnator as an assay-specific MPRA design tool. Its narrow focus on SNP-centric MPRAs exemplifies the manuscript's claim that existing assay-specific tools "are not easily adapted to novel designs." The pair of MPRAnator + MPRA Design Tools makes the case stronger than either alone, as they show two independent tools both limited to specific MPRA sub-types.

No additional citation locations needed.

## 4. Novelty Assessment

**No threat to novelty claims.** MPRA Design Tools is even more narrowly scoped than MPRAnator — it handles only SNP-focused MPRAs with a fixed oligo layout. It has no combinatorial design capabilities, no motif grammar support, no DMS functionality, and no composable architecture. The power analysis module is orthogonal to PoolParty's contributions (it addresses experimental design statistics, not sequence design). No special precautions needed.

## 5. References to Investigate

- **Tewhey et al. (2016)** — "Direct identification of hundreds of expression-modulating variants using a multiplexed reporter assay." Cell 165:1519-1529. Landmark SNP-MPRA paper. Already well-known; not needed for the PoolParty manuscript unless adding more MPRA background references. **No action needed.**
- **Ulirsch et al. (2016)** — "Systematic functional dissection of common genetic variation affecting red blood cell traits." Cell 165:1530-1545. Another major SNP-MPRA study. Same assessment. **No action needed.**

## 6. Citing Papers to Investigate

Searched via Semantic Scholar (DOI: 10.1093/bioinformatics/bty150). 8 citing papers total. All are MPRA applications or methodology papers.

- **Gordon et al. (2020)** — "lentiMPRA and MPRAflow." 118 citations. MPRA protocol/analysis workflow. Not a library design tool. **No action needed.**
- **Zheng & VanDusen (2023)** — MPRA methods review. 15 citations. Could serve as a general MPRA reference, but manuscript already has adequate MPRA citations. **No action needed.**

No citing paper poses a novelty concern or fills a gap in the manuscript.

---
**Status:** Keep citation as-is. Correctly categorized as an assay-specific MPRA design tool. Narrower than MPRAnator (SNP-focused only). No novelty concerns. No changes to manuscript needed.
