# Hossain et al. (2024) — Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements

## BibTeX Entry

```bibtex
@article{Hossain2024oc,
  year     = {2024},
  title    = {Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements},
  author   = {Hossain, Ayaan and Cetnar, Daniel P and LaFleur, Travis L and McLellan, James R and Salis, Howard M},
  journal  = {ACS Synth Biol},
  doi      = {10.1021/acssynbio.4c00661},
  pmid     = {39641628},
  pages    = {4218--4232},
  number   = {12},
  volume   = {13},
  mynotes  = {Python suite for designing the physical infrastructure of oligo pools (barcodes, primer binding sites, spacers, overlap regions) and analyzing barcoded MPRA sequencing reads. Assumes variant sequences are already designed as input; does not address variant content design (which mutations, motifs, or combinatorial logic to use). Complementary to PoolParty: PoolParty designs the variant content, Oligopool Calculator designs the surrounding infrastructure. Could optionally be cited in the Introduction or Discussion as a complementary tool focused on oligo pool construction rather than variant library specification.}
}
```

## 1. Summary

The Oligopool Calculator is an end-to-end Python suite for designing and analyzing oligonucleotide pools. In **Design Mode**, it takes pre-designed variant sequences as input and designs optimized barcodes (with orthogonally symmetric selection for maximal Hamming distance), universal primer binding sites (using adaptive decision trees and the Non-Repetitive Parts Calculator), spacers, padding, and overlap regions for multi-fragment assembly. In **Analysis Mode**, it indexes barcoded sequencing reads from MPRAs, packs reads efficiently, and counts variant frequencies to produce count matrices. Demonstrated on three projects: 93,180 promoter variants, 62,120 mRNA stability elements, and 6,232 ribozymes. Performance: 4M barcodes in 1.2 hours, primers for 1M oligos in 15 minutes, ~500M reads/hour analysis. Published December 2024; 2 citations so far.

## 2. Relevance to PoolParty Manuscript

The Oligopool Calculator and PoolParty address **different stages** of the same overall workflow:

- **PoolParty** designs the *variant content* of a library: which sequences to include, what mutations/motifs/combinations to generate, with structured metadata (design cards).
- **Oligopool Calculator** designs the *physical infrastructure* around pre-existing variants: barcodes, primers, spacers, and overlap regions needed for synthesis, PCR, cloning, and sequencing.

Critically, the Oligopool Calculator's Design Mode takes "Designed Variants" as *input* (see their Figure 1a). It does not decide what variants to make — that is the user's (or PoolParty's) job. The tools are complementary, not competing.

| Feature | Oligopool Calculator | PoolParty |
|---------|---------------------|-----------|
| Purpose | Physical oligo pool construction | Variant content design |
| Input | Pre-designed variant sequences | Template sequences + Operations |
| Output | Oligos with barcodes, primers, spacers | Library of variant sequences + metadata |
| Barcode design | Core feature (advanced algorithms) | Supported (get_barcodes) |
| Primer design | Core feature | Not in scope |
| Mutagenesis | Not in scope | Core feature (20+ Operations) |
| Combinatorial logic | Not in scope | Core feature (DAG, stack, join) |
| Design cards | No | Core feature |
| MPRA read analysis | Yes (Analysis Mode) | Not in scope |

## 3. Citation Recommendations

**Optional citation.** Could be cited in:
- The **Introduction** when noting that existing tools address parts of the process but not the unified library design problem — Oligopool Calculator handles downstream construction but not upstream variant specification.
- The **Discussion** as a complementary tool that could be paired with PoolParty in practice.

Not essential, but citing it would demonstrate awareness of the broader oligo pool design ecosystem and strengthen the claim that PoolParty fills a specific gap (variant content design) not addressed by other tools.

## 4. Novelty Assessment

**No threat to novelty claims.** The Oligopool Calculator operates at a different stage of the pipeline. It does not design variant content, has no composable Operations, no DAG architecture, no design cards, no lazy evaluation, and no mutagenesis capabilities. No special precautions needed.

## 5. References to Investigate

None of the references in the Oligopool Calculator paper appear to be missed library design tools. The paper cites many MPRA application papers and primer/barcode design tools (Primer3, FreeBarcodes, BARCOSEL, SiteOut) that are in a different category.

## 6. Citing Papers to Investigate

Only 2 citations so far (very recent paper). Neither is relevant to library design software.

---
**Status:** Not currently cited. Optional citation as a complementary tool. No novelty concerns.
