# Zulkower & Rosser (2020) — DNA Chisel, a versatile sequence optimizer

## BibTeX Entry

```bibtex
@article{Zulkower2020jk,
  year     = {2020},
  title    = {{DNA} Chisel, a versatile sequence optimizer},
  author   = {Zulkower, Valentin and Rosser, Susan},
  journal  = {Bioinformatics},
  issn     = {1367-4803},
  doi      = {10.1093/bioinformatics/btaa558},
  pmid     = {32647895},
  pages    = {4508--4509},
  number   = {16},
  volume   = {36},
  mynotes  = {Python framework for optimizing individual DNA sequences against composable constraints and objectives (codon usage, GC content, restriction site removal, homology avoidance, etc.). Operates on single sequences, not libraries. Should be cited in the Introduction as an example of a sequence optimization tool that addresses synthesis constraints on individual sequences but not at the level of variant libraries, motivating the need for PoolParty.}
}
```

## 1. Summary

DNA Chisel is a Python library, command-line tool, and web application for optimizing DNA sequences against multiple composable constraints and objectives. Users define specifications (hard constraints like "avoid BsaI sites" or "enforce translation in this region," and soft objectives like "match E. coli codon usage" or "enforce 40-60% GC content") that are applied to regions of a starting sequence. The solver uses a two-step algorithm: first resolve all hard constraints, then maximize objectives. Over 15 built-in specification classes handle codon optimization, GC content enforcement, restriction site avoidance, homology removal, and more. Specifications can be defined via Python scripts or Genbank annotations. DNA Chisel is extensible via user-defined specification classes. It produces optimization reports for traceability. Built on Biopython. 57 citations.

## 2. Relevance to PoolParty Manuscript

DNA Chisel and PoolParty address fundamentally different problems:

- **DNA Chisel** optimizes a single DNA sequence to satisfy manufacturing and biological constraints (codon usage, GC content, restriction sites, etc.). It modifies an existing sequence to make it "better" according to defined criteria.
- **PoolParty** designs libraries of many variant sequences for functional assays or in silico experiments. It generates collections of sequences that systematically explore sequence space.

DNA Chisel has no concept of variant libraries, mutagenesis, combinatorial design, or design cards. Its composable specification system is architecturally interesting (specifications can be combined via Python scripts, similar in spirit to PoolParty's composable Operations), but the purpose and output are entirely different: one optimized sequence vs. thousands of variant sequences.

The tools could be complementary in practice: a user might use PoolParty to design a variant library, then use DNA Chisel to post-process individual sequences for synthesis constraints (e.g., removing problematic restriction sites while preserving the intended mutations).

## 3. Citation Recommendations

**Should be cited** in the Introduction (paragraph 3) as an example of a sequence optimization tool that "address[es] synthesis constraints on individual sequences, but not at the level of libraries." This is the current placement, and it's accurate.

No additional citation locations needed.

## 4. Novelty Assessment

**No threat to novelty claims.** DNA Chisel operates on individual sequences for optimization, not on libraries for variant generation. The composable specification system is a shared design philosophy (both tools allow combining modular operations), but applied to different domains. No special precautions needed.

## 5. References to Investigate

- **D-tailor (Guimaraes et al. 2014)** — "Automated analysis and design of DNA sequences." Bioinformatics 30:1087-1094. A predecessor to DNA Chisel that enables programmatic definition and combination of sequence design specifications via Python scripts. Could be relevant as another sequence design tool, but focused on single-sequence optimization, not libraries. **Low priority.**

## 6. Citing Papers to Investigate

Searched via Semantic Scholar (PMID:32647895). 57 total citations; top 20 reviewed.

- **Tycko et al. (2020)** — "High-throughput discovery and characterization of human transcriptional effectors." 135 citations. Used DNA Chisel for codon optimization of DMS library oligos. DNA Chisel served as a post-processing step, not as the library design tool. Illustrates the complementary relationship. **No action needed.**
- **Lund et al. (2023)** — "Highly Parallelized Construction of DNA from Low-Cost Oligonucleotide Mixtures." 14 citations. DNA assembly from oligo pools. About constructing DNA, not designing variant libraries. **No action needed.**
- **James et al. (2024)** — "The design and engineering of synthetic genomes." 19 citations. Genome-scale DNA design review. Tangential. **No action needed.**

No citing paper poses a novelty concern or fills a gap in the PoolParty manuscript.

---
**Status:** Keep citation as-is. Correctly categorized as a sequence optimization tool operating on individual sequences. No novelty concerns. No changes to manuscript needed.
