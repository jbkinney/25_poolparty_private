# Georgakopoulos-Soares et al. (2016) — MPRAnator: a web-based tool for the design of massively parallel reporter assay experiments

## BibTeX Entry

```bibtex
@article{GeorgakopoulosSoares2017gb,
  year     = {2016},
  title    = {{MPRAnator}: a web-based tool for the design of massively parallel reporter assay experiments},
  author   = {Georgakopoulos-Soares, Ilias and Jain, Naman and Gray, Jesse M and Hemberg, Martin},
  journal  = {Bioinformatics},
  doi      = {10.1093/bioinformatics/btw584},
  pmid     = {27605100},
  pmcid    = {{PMC}5198521},
  pages    = {137--138},
  number   = {1},
  volume   = {33},
  mynotes  = {Web-based tool for designing MPRA experiments, with modules for motif placement, SNP variant design, PWM-based sequence generation, and negative control design (scrambling/mutating). Should be cited in the Introduction as an example of an assay-specific tool that automates library design for a particular experimental format (MPRAs) but is not easily adapted to novel designs such as DMS libraries or in silico experiments.}
}
```

Note: Advance Access published September 2016; journal issue is Bioinformatics 33(1), 2017. The year=2016 matches the online publication date.

## 1. Summary

MPRAnator is a web-based tool set for designing MPRA experiments. It provides four modules: (1) **MPRA Motif design** — systematically generates synthetic sequences with single or combinatorial motifs placed at pre-selected positions to study TF occupancy rules; (2) **MPRA SNP design** — examines the regulatory effects of single or combinations of SNPs at regulatory sequences; (3) **PWM Seq-Gen** — generates sequences from PWMs (probabilistic realizations or exhaustive k-mers above a threshold); (4) **Transmutation** — designs negative controls via scrambling, reversing, complementing, or random mutagenesis. Users can incorporate sub-components (barcodes, adapters, restriction sites) and arrange them via drag-and-drop. Implemented in Python, Perl, and Javascript with a REST API.

## 2. Relevance to PoolParty Manuscript

MPRAnator is the most directly relevant prior tool to PoolParty's MPRA use case. Both tools address the problem of designing sequence libraries for MPRA experiments, and both handle motif placement, combinatorial designs, and sub-components like barcodes.

Key differences:

| Feature | MPRAnator | PoolParty |
|---------|-----------|-----------|
| Scope | MPRA-specific | General-purpose library design (MPRA, DMS, in silico, etc.) |
| Interface | Web-based GUI with drag-and-drop | Python library with composable syntax |
| Architecture | Four separate tools | Unified DAG of composable Operations |
| Assay types | MPRA only | MPRA, DMS, in silico experiments, and more |
| Codon-level operations | None | mutagenize_orf and other codon-aware Operations |
| Design cards | None | Core feature |
| Lazy evaluation | No | Yes (on-demand generation) |
| Sequence metadata | None | Names, styling, design cards |
| Extensibility | Limited (web tool) | User-defined Operations via subclassing |

MPRAnator handles some of the same operations as PoolParty (motif placement, SNP design, barcode addition) but is constrained to MPRA experiments and a web GUI. It cannot handle DMS libraries, in silico experiments, or arbitrary composable designs.

## 3. Citation Recommendations

**Should be cited** in the Introduction (paragraph 3) as the primary example of an assay-specific tool that automates MPRA library design but is not easily adapted to novel designs. This is the strongest reference for that specific claim because it is explicitly MPRA-focused and shares the most functional overlap with PoolParty's MPRA capabilities.

No additional citation locations needed.

## 4. Novelty Assessment

**Low threat to novelty, but worth being precise about.** MPRAnator is the closest prior art for PoolParty's MPRA functionality. However, PoolParty's novelty claims are safe because:

1. **Broader scope**: PoolParty handles DMS, in silico experiments, and arbitrary designs — not just MPRAs.
2. **Composable architecture**: PoolParty's DAG-based design with composable Operations is fundamentally different from MPRAnator's four separate tools.
3. **Design cards and metadata**: MPRAnator generates sequences but not structured metadata about how each sequence was constructed.
4. **Lazy evaluation**: MPRAnator generates all sequences upfront; PoolParty generates on demand.
5. **Python API**: PoolParty's programmatic interface enables integration into computational workflows (e.g., surrogate modeling); MPRAnator is a web tool.

The current manuscript text ("assay-specific tools automate library design for particular experimental formats but are not easily adapted to novel designs") correctly distinguishes PoolParty from MPRAnator. No additional precautions needed.

## 5. References to Investigate

- **Ghazi et al. (2018)** — "Design tools for MPRA experiments." Already in our to_analyze list. Another MPRA design tool; will be analyzed separately.
- **Gordon et al. (2020)** — "lentiMPRA and MPRAflow for high-throughput functional characterization of gene regulatory elements." 118 citations. This is an MPRA protocol/workflow paper, not a library design tool. **Low priority** — could be cited as background on MPRA methodology if needed, but not directly relevant to PoolParty's software contribution.

## 6. Citing Papers to Investigate

Searched via Semantic Scholar (DOI: 10.1093/bioinformatics/btw584). 9 citing papers total.

- **Ghazi et al. (2018)** — "Design tools for MPRA experiments." 8 citations. Already in our to_analyze list; will be analyzed separately.
- **MPRAbase (Zhao et al. 2023)** — Database of 129 MPRA experiments. Not a design tool. **No action needed.**
- **Zheng & VanDusen (2023)** — MPRA methodology review. Could be worth citing as a general MPRA reference, but the manuscript already cites Inoue2015 and Kinney2019 for this purpose. **No action needed.**
- **Gordon et al. (2020)** — lentiMPRA/MPRAflow. MPRA protocol, not design software. **No action needed.**

No citing paper poses a novelty concern or fills a gap in the manuscript.

---
**Status:** Keep citation as-is. Correctly categorized as an assay-specific MPRA design tool. Closest prior art for MPRA functionality, but PoolParty's broader scope, composable architecture, and metadata features clearly distinguish it. No changes to manuscript needed.
