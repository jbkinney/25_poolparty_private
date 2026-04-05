# Swainston et al. (2017) — CodonGenie: optimised ambiguous codon design tools

## BibTeX Entry

```bibtex
@article{Swainston2017rb,
  year     = {2017},
  title    = {{CodonGenie}: optimised ambiguous codon design tools},
  author   = {Swainston, Neil and Currin, Andrew and Green, Lucy and Breitling, Rainer and Day, Philip J and Kell, Douglas B},
  journal  = {PeerJ Comput Sci},
  doi      = {10.7717/peerj-cs.120},
  pages    = {e120},
  volume   = {3},
  mynotes  = {Web application for designing optimal ambiguous codons (degenerate codons) for protein mutagenesis, ranking them by encoding efficiency and organism-specific codon usage. Operates at the single-codon level, not at the level of full sequence libraries. Should be cited in the Introduction as an example of a sequence optimization tool that addresses codon-level design for individual positions but not the design of complete variant libraries, motivating the need for PoolParty.}
}
```

## 1. Summary

CodonGenie is a web application and RESTful web service for designing optimal ambiguous (degenerate) codons for protein mutagenesis and directed evolution. Given a user-specified set of target amino acids and a host organism, CodonGenie enumerates all ambiguous codons that encode the desired amino acids, then ranks them by a scoring function that balances encoding efficiency (minimizing off-target amino acids and stop codons) with organism-specific codon usage preferences. For example, to encode the five non-polar amino acids F, I, L, M, V, CodonGenie identifies that DTK (6 DNA variants) is more efficient than the naive NTN (16 variants) for E. coli, while DTS is preferred for S. coelicolor. An Analyse module lets users check existing ambiguous codons. The tool also provides a web service API for programmatic integration into larger pipelines (e.g., iterating over a multiple sequence alignment to build a degenerate synthetic gene). Written in Python/Flask, ~700 lines of code. 11 citations.

## 2. Relevance to PoolParty Manuscript

CodonGenie operates at a fundamentally different granularity than PoolParty:

- **CodonGenie** designs a single optimal ambiguous codon for a single position, given a set of desired amino acids. It is a codon-level optimization tool.
- **PoolParty** designs complete libraries of variant sequences across entire genes or regulatory elements, with composable Operations, design cards, and lazy evaluation.

CodonGenie has no concept of full sequence libraries, DAGs, composable operations, design cards, or multi-region designs. It solves one narrow sub-problem (which degenerate codon to use at one position) that could arise as a component within a larger library design workflow. PoolParty's `mutagenize_orf` Operation handles codon-level mutagenesis internally as part of its library generation, but uses a different approach (specifying explicit codons for each amino acid substitution rather than degenerate codons).

## 3. Citation Recommendations

**Should be cited** in the Introduction (paragraph 3) alongside DNA Chisel as a sequence optimization tool that addresses constraints at the individual sequence/codon level but not at the level of libraries. CodonGenie is a good complement to DNA Chisel in this category because it specifically addresses codon design for mutagenesis — a topic directly adjacent to PoolParty's domain — yet operates at a completely different scale (single codon vs. full library).

No additional citation locations needed.

## 4. Novelty Assessment

**No threat to novelty claims.** CodonGenie solves a narrow optimization problem (best degenerate codon for a position) that is orthogonal to PoolParty's contributions. PoolParty does not use degenerate codons — it generates explicit sequences for each variant. No special precautions needed.

## 5. References to Investigate

- **DYNAMCC (Halweg-Edwards et al. 2016)** — "A Web Interface for Codon Compression." ACS Synth Biol 5:1021-1023. Designs sets of ambiguous codons to encode amino acid sets with minimal redundancy. Same domain as CodonGenie, even narrower. **No action needed** — not relevant to PoolParty's library-level design.
- **Hiraga et al. (2021)** — "Mutation Maker, An Open Source Oligo Design Platform for Protein Engineering." ACS Synth Biol 10:357-370. Already in the bib file (commented out as `article{Hiraga2021yg`). Open-source platform for designing mutagenic oligos for protein engineering. **Medium priority** — could be worth investigating as it designs oligos for multi-site mutagenesis, which is closer to library design than CodonGenie. However, it appears focused on primer/oligo design for site-directed mutagenesis, not on designing complete variant libraries for high-throughput assays.
- **GeneGenie (Swainston et al. 2014)** — "Optimized oligomer design for directed evolution." NAR 42:W395-W400. Same group, predecessor tool. Designs oligomers for directed evolution. **Low priority.**

## 6. Citing Papers to Investigate

Searched via Semantic Scholar (DOI: 10.7717/peerj-cs.120). 11 total citations.

- **GGAssembler (Hoch et al. 2024)** — "Precise and economical design and synthesis of combinatorial mutation libraries." 4 citations. Graph-theoretical method for designing DNA fragments for combinatorial variant library assembly. **Medium priority** — this is the most relevant citing paper, as it addresses combinatorial library design. Worth being aware of, but appears focused on the DNA fragment assembly strategy (how to physically construct libraries) rather than on specifying library content (what sequences to include), which is PoolParty's domain.

- **Hiraga et al. (2020)** — "Mutation Maker." 7 citations. Same as noted above.

No citing paper poses a novelty concern or needs to be added to the PoolParty manuscript.

---
**Status:** Keep citation as-is. Correctly categorized as a codon-level optimization tool. No novelty concerns. No changes to manuscript needed.
