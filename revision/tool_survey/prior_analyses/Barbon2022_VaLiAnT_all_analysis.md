# Barbon et al. (2021/2022) — VaLiAnT: an oligonucleotide library design and annotation tool for saturation genome editing and other deep mutational scanning experiments

## BibTeX Entry

```bibtex
@article{Barbon2022th,
  year     = {2021},
  title    = {Variant Library Annotation Tool ({VaLiAnT}): an oligonucleotide library design and annotation tool for saturation genome editing and other deep mutational scanning experiments},
  author   = {Barbon, Luca and Offord, Victoria and Radford, Elizabeth J and Butler, Adam P and Gerety, Sebastian S and Adams, David J and Tan, Hong Kee and Waters, Andrew J},
  journal  = {Bioinformatics},
  issn     = {1367-4803},
  doi      = {10.1093/bioinformatics/btab776},
  pmid     = {34791067},
  pmcid    = {{PMC}8796380},
  pages    = {892--899},
  number   = {4},
  volume   = {38},
  mynotes  = {Python command-line tool for designing oligonucleotide variant libraries for saturation genome editing (SGE) and DMS experiments. Generates exhaustive single-position variants (SNVs, deletions, alanine scan, etc.) with codon-aware annotation. Should be cited in the Introduction as the most relevant example of an assay-specific tool that automates library design for a particular experimental format (SGE/DMS saturation mutagenesis) but is not easily adapted to novel designs involving combinatorial variation, motif grammar studies, or in silico experiments. This is the closest prior art to PoolParty and the distinction should be drawn carefully.}
}
```

Note: Advance Access published November 2021; journal issue is Bioinformatics 38(4), 2022. The year=2021 matches the online publication date.

## 1. Summary

VaLiAnT (Variant Library Annotation Tool) is a Python command-line tool for designing oligonucleotide variant libraries for saturation genome editing (SGE) and deep mutational scanning (DMS) experiments. Starting from genomic coordinates, it retrieves reference sequences and applies "mutator functions" to generate systematic variant libraries. Mutator functions include: `snv` (all single-nucleotide substitutions), `1del`/`2del0`/`2del1` (single and dinucleotide deletions), `snvre` (SNVs with synonymous codon replacement for controls), `ala` (alanine scan), `stop` (stop codon substitutions), and `inframe` (in-frame deletions). CDS-aware functions classify variants as synonymous/missense/nonsense. Users define "targetons" (genomic regions) with sub-regions that can receive different mutator functions independently. VaLiAnT also supports cDNA DMS and prime editing saturation libraries, handles PAM/protospacer protection edits for CRISPR compatibility, and generates rich metadata (mutation type, amino acid changes, variant classification). Demonstrated on BRCA1 exons. 193 citations as of early 2026 — well-adopted in the SGE community.

## 2. Relevance to PoolParty Manuscript

**VaLiAnT is the closest prior art to PoolParty.** Both are Python command-line tools that design variant libraries with metadata annotation. Key comparisons:

**Shared territory:**
- Both generate variant libraries for DMS/mutagenesis experiments
- Both are Python-based command-line tools
- Both handle codon-aware mutagenesis (synonymous/missense/nonsense classification)
- Both generate metadata for each variant (VaLiAnT: CSV/VCF output; PoolParty: design cards)

**Critical differences:**

| Feature | VaLiAnT | PoolParty |
|---------|---------|-----------|
| Design paradigm | Genome-centric (coordinates + mutator functions) | Sequence-centric (composable DAG of Operations) |
| Mutation types | Exhaustive single-position saturation only | Sequential (exhaustive), random (sampled), and fixed modes |
| Combinatorial design | No (single-position variants only) | Yes (stack, join, pairwise/higher-order mutations) |
| MPRA support | No | Yes (motif insertion, insert_multiscan, barcodes) |
| In silico experiments | Not designed for this | Core use case (design cards as covariates) |
| Motif operations | None | PWM sampling, IUPAC motifs, motif insertion/scanning |
| Barcode generation | No | Yes (get_barcodes) |
| Deletion/insertion scans | Limited (1del, 2del) | Comprehensive (deletion_scan, insertion operations) |
| Lazy evaluation | No (generates all variants upfront) | Yes (on-demand generation) |
| Composable Operations | No (fixed mutator functions) | Yes (DAG of chainable Operations) |
| Sequence styling | No | Yes (colored/formatted output) |
| Extensibility | Limited (fixed mutator set) | User-defined Operations via subclassing |
| SGE-specific features | PAM protection, protospacer edits, prime editing | Not applicable |

VaLiAnT excels at its specific niche: exhaustive saturation mutagenesis with genomic context and CRISPR-compatibility features for SGE. PoolParty is broader and more composable but does not handle SGE-specific features like PAM protection.

## 3. Citation Recommendations

**Should be cited** in the Introduction (paragraph 3) as the most relevant example of an assay-specific tool for DMS/SGE library design. It is the strongest reference supporting the claim that existing assay-specific tools are tailored to particular experimental formats. Among the three assay-specific tools cited (MPRAnator, MPRA Design Tools, VaLiAnT), VaLiAnT has the most functional overlap with PoolParty and thus makes the distinction most meaningful.

No additional citation locations needed, though a brief mention in the Discussion acknowledging VaLiAnT's SGE-specific strengths (PAM protection, prime editing support) while noting that PoolParty provides a more general framework could strengthen the paper. This is optional.

## 4. Novelty Assessment

**Moderate relevance — requires careful positioning but not a threat.** VaLiAnT is the one tool that most overlaps with PoolParty's DMS capabilities. However, PoolParty's novelty claims are secure because:

1. **Composable architecture**: VaLiAnT applies fixed mutator functions to genomic regions; PoolParty chains Operations into a DAG. This is a fundamental architectural difference.
2. **Combinatorial design**: VaLiAnT generates only exhaustive single-position variants. PoolParty supports pairwise mutations, higher-order random mutants, and arbitrary combinations via stack/join.
3. **Broader scope**: PoolParty handles MPRA grammar libraries and in silico experiments — domains VaLiAnT does not address.
4. **Design cards as covariates**: PoolParty's design cards enable downstream analyses (surrogate modeling). VaLiAnT's metadata is variant annotation, not experimental covariates.
5. **Lazy evaluation**: PoolParty's on-demand generation allows exploration before committing to full library generation.

The current manuscript text ("assay-specific tools automate library design for particular experimental formats but are not easily adapted to novel designs") is accurate but could be slightly more specific. Consider adding a clause noting that these tools focus on exhaustive single-position saturation mutagenesis, which would sharpen the contrast with PoolParty's combinatorial and multi-mode capabilities. This is optional but recommended.

## 5. References to Investigate

- **Findlay et al. (2018)** — "Accurate classification of BRCA1 variants with saturation genome editing." Nature 562:217-222. Landmark SGE paper that motivated VaLiAnT. Not a software tool. Already well-known. **No action needed** for PoolParty manuscript.
- **Erwood et al. (2020)** — "Saturation variant interpretation using CRISPR prime editing." Referenced for prime editing saturation. **No action needed.**

## 6. Citing Papers to Investigate

Searched via Semantic Scholar (PMID:34791067). 193 total citations; API returned 7 representative papers. All are SGE applications or DMS methodology papers (BRCA1, RAD51C, BAP1, DDX3X saturation editing). None are library design software tools.

- **Waters et al. (2024)** — "Saturation genome editing of BAP1." 37 citations. SGE application using VaLiAnT. Not a competing tool. **No action needed.**
- **Olvera-León et al. (2024)** — "High-resolution functional mapping of RAD51C by SGE." 29 citations. Same. **No action needed.**
- **Cooper et al. (2024)** — "Analyzing the functional effects of DNA variants with gene editing." 8 citations. Review of gene editing approaches for variant analysis. Could cite PoolParty in the future, but not needed in our manuscript. **No action needed.**

No citing paper poses a novelty concern or needs to be added to the PoolParty manuscript.

---
**Status:** Keep citation as-is. Closest prior art — correctly categorized as assay-specific (SGE/DMS saturation). PoolParty's composable architecture, combinatorial design, broader scope, and design cards clearly distinguish it. Consider optional sharpening of the distinction in the Introduction. No major changes to manuscript needed.
