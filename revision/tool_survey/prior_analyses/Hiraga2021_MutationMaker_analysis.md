# Hiraga et al. (2021) — Mutation Maker, An Open Source Oligo Design Platform for Protein Engineering

## BibTeX Entry

```bibtex
@article{Hiraga2021yg,
  year     = {2021},
  title    = {Mutation Maker, An Open Source Oligo Design Platform for Protein Engineering},
  author   = {Hiraga, Kaori and Mejzlik, Petr and Marcisin, Matej and Vostrosablin, Nikita and Gromek, Anna and Arnold, Jakub and Wiewiora, Sebastian and Svarba, Rastislav and Prihoda, David and Clarova, Kamila and Klempir, Ondrej and Navratil, Josef and Tupa, Ondrej and Vazquez-Otero, Alejandro and Walas, Marcin W and Holy, Lukas and Spale, Martin and Kotowski, Jakub and Dzamba, David and Temesi, Gergely and Russell, Jay H and Marshall, Nicholas M and Murphy, Grant S and Bitton, Danny A},
  journal  = {ACS Synth Biol},
  doi      = {10.1021/acssynbio.0c00542},
  pmid     = {33433999},
  pages    = {357--370},
  number   = {2},
  volume   = {10},
  mynotes  = {Open-source software for designing mutagenic primers and oligos for protein engineering, supporting site-scanning saturation mutagenesis (SSSM), multi site-directed mutagenesis (MSDM), and de novo gene synthesis (PAS). Focuses on primer design given user-specified mutations, not on deciding which variants to include in a library. Could optionally be cited in the Introduction as another assay-specific tool for DMS/protein engineering that focuses on the physical realization (primer/oligo design) rather than on library-level variant specification.}
}
```

## 1. Summary

Mutation Maker is an open-source Python-based platform for designing mutagenic oligos at industrial scale. It supports three protein engineering workflows: (1) **SSSM** (site-scanning saturation mutagenesis) — designs overlapping mutagenic primer pairs for parallel PCR reactions to introduce random amino acid substitutions at specified positions; (2) **MSDM** (multi site-directed mutagenesis) — designs unidirectional mutagenic primers carrying degenerate codons to introduce specific combinations of mutations at multiple sites simultaneously (using the QCLM kit); (3) **PAS** (PCR-based accurate synthesis) — designs overlapping mutagenic oligos for de novo gene synthesis with user-specified mutations. Uses novel brute-force and fast-approximation algorithms for primer optimization, incorporating constraints on melting temperature, GC content, hairpins, primer-dimers, and coverage. Experimentally validated with ~80% success rates across multiple genes and vectors.

## 2. Relevance to PoolParty Manuscript

Mutation Maker and PoolParty address different aspects of variant library creation:

- **Mutation Maker** designs the *physical primers/oligos* needed to introduce mutations via wet-lab protocols (PCR, QCLM, PAS). It takes user-specified mutation sites as input and outputs optimized primers.
- **PoolParty** designs the *variant content* of a library: which sequences to include, with what combinations of mutations/motifs, generating complete variant sequences with structured metadata.

Key distinction: Mutation Maker does not decide which variants to make — it designs primers to make the variants the user has already specified. PoolParty decides which variants to make and generates their complete sequences.

| Feature | Mutation Maker | PoolParty |
|---------|---------------|-----------|
| Purpose | Design mutagenic primers/oligos | Design variant sequence libraries |
| Input | Gene + list of mutation sites | Template sequences + Operations |
| Output | Primer sequences + PCR protocols | Complete variant sequences + metadata |
| Decides which variants? | No (user specifies) | Yes (Operations generate variants) |
| Wet-lab protocol aware | Yes (SSSM, MSDM, PAS) | No (protocol-agnostic) |
| Combinatorial logic | Limited (MSDM combinations) | Extensive (DAG, stack, join, etc.) |
| MPRA support | No | Yes |
| Design cards | No | Core feature |

## 3. Citation Recommendations

**Optional citation.** Could be cited in the Introduction alongside VaLiAnT as another assay-specific tool for protein engineering that focuses on primer/oligo design. It strengthens the argument that existing tools address specific aspects of the library creation process (primer design, saturation mutagenesis) but not the unified library specification problem that PoolParty solves.

However, since Mutation Maker is about primer design rather than library content design, it is less directly relevant than VaLiAnT (which does generate variant sequences). **Not essential.**

## 4. Novelty Assessment

**No threat to novelty claims.** Mutation Maker operates at a different level (primer design for wet-lab protocols) than PoolParty (library content specification). No overlap in architectural design (no DAG, no composable Operations, no design cards). No special precautions needed.

## 5. References to Investigate

None relevant to PoolParty's library design domain.

## 6. Citing Papers to Investigate

Semantic Scholar lookup failed (rate limit), but web search indicates the paper has ~7 citations, none of which appear to be library design software tools.

---
**Status:** Not currently cited. Optional citation to broaden the "existing tools" landscape. No novelty concerns.
