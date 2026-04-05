# Pereira et al. (2015) — Pydna: a simulation and documentation tool for DNA assembly strategies using python

## BibTeX Entry

```bibtex
@article{Pereira2015wj,
  year     = {2015},
  title    = {Pydna: a simulation and documentation tool for {DNA} assembly strategies using python},
  author   = {Pereira, Filipa and Azevedo, Fl{\'a}vio and Carvalho, {\^A}ngela and Ribeiro, Gabriela F and Budde, Mark W and Johansson, Bj{\"o}rn},
  journal  = {BMC Bioinformatics},
  doi      = {10.1186/s12859-015-0544-x},
  pmid     = {25933606},
  pmcid    = {{PMC}4472420},
  pages    = {142},
  number   = {1},
  volume   = {16},
  mynotes  = {Python library for simulating molecular biology DNA operations (restriction digestion, ligation, PCR, Gibson assembly, homologous recombination). Operates on individual DNA constructs, not variant libraries. Should be cited in the Introduction as an example of a general-purpose Python sequence toolkit that can manipulate individual sequences but has no notion of a variant library as a structured object, motivating the need for PoolParty.}
}
```

## 1. Summary

Pydna is a free, open-source Python library for simulating basic molecular biology DNA unit operations: restriction digestion, ligation, PCR, primer design, Gibson assembly, and homologous recombination. A cloning strategy expressed as a pydna script serves as a complete, unambiguous, human-readable description that can be automatically executed to yield the sequences of final and intermediate constructs. Pydna is built on Biopython and NetworkX, and its central data type (`Dseqrecord`) represents double-stranded DNA molecules. The paper demonstrates pydna through three examples: (1) construction of an expression vector by restriction digestion and ligation, (2) construction of a plasmid by homologous recombination, and (3) assembly of a two-gene metabolic pathway combining cut-and-paste cloning with homologous recombination. The authors compare pydna to GUI-based tools like j5 and RavenCAD, arguing that Python's flexibility and extensibility are advantages over purpose-built languages or graphical interfaces. The paper had ~25 citing papers as of early 2026.

## 2. Relevance to PoolParty Manuscript

Pydna and PoolParty share the philosophy of using Python code as a declarative description of a DNA construction process, but they address fundamentally different problems:

- **Pydna** simulates individual molecular biology operations (cut, ligate, PCR, assemble) to produce and verify the sequence of a single DNA construct. It operates on one molecule at a time.
- **PoolParty** designs *libraries* of variant sequences (thousands to millions) using a DAG of composable Operations that generate sequences on demand with structured metadata (design cards, names, styling).

Pydna has no concept of a variant library, mutagenesis, combinatorial sequence variation, design cards, lazy evaluation, or DAG-based library specification. It is a cloning simulation tool, not a library design tool.

The manuscript correctly categorizes pydna among "general-purpose sequence toolkits" that "can manipulate individual sequences but have no notion of a variant library as a structured object" (Introduction, paragraph 3).

## 3. Citation Recommendations

**Current citation**: Introduction, paragraph 3 — cited alongside Biopython, tangermeme, and SeqPro as a general-purpose sequence toolkit.

**Assessment**: The citation is accurate and appropriately placed. Pydna is a reasonable example of a Python-based DNA tool that operates at the individual-sequence level. No changes needed.

**Optional addition**: Pydna could also be mentioned in the Discussion as an example of the "Python script as documentation" paradigm applied to a different domain (cloning simulation vs. library design). This would acknowledge shared design philosophy while clarifying the distinct scope. However, this is not necessary and might distract from the main points.

## 4. Novelty Assessment

**No threat to novelty claims.** Pydna operates in a completely different problem domain (cloning simulation of individual constructs) than PoolParty (variant library design). Specific distinctions:

| Feature | Pydna | PoolParty |
|---------|-------|-----------|
| Purpose | Simulate cloning of individual constructs | Design libraries of variant sequences |
| Output | One DNA sequence per script | Thousands to millions of sequences |
| Core abstraction | Dseqrecord (one dsDNA molecule) | Pool + Operation DAG |
| Combinatorics | Not supported | Central feature |
| Mutagenesis | Not supported | 20+ built-in Operations |
| Design cards | Not available | Core feature |
| Lazy evaluation | No (immediate execution) | Yes (on-demand generation) |

No special precautions are needed. The current manuscript text adequately distinguishes these tools.

## 5. References to Investigate

From pydna's reference list and comparison section:

- **Hillson et al. (2011)** — "j5 DNA assembly design automation software," ACS Synth Biol 1:14-21. DNA assembly design tool with cost optimization and part reuse. Different domain (assembly automation, not variant libraries), but worth knowing about as a related tool in the DNA design space. **Low priority** — already a different problem domain.
- **Appleton et al. (2014)** — "Interactive assembly algorithms for molecular cloning," Nat Methods 11:657-62 (RavenCAD). Interactive DNA assembly design. Same assessment as j5. **Low priority.**

Neither j5 nor RavenCAD design variant libraries; they automate cloning strategy design. They are not relevant enough to cite in the PoolParty manuscript.

## 6. Citing Papers to Investigate

Searched via Semantic Scholar (DOI: 10.1186/s12859-015-0544-x). ~25 citing papers total. Most relevant:

- **Mori & Yachie (2022)** — "QUEEN: A framework to efficiently describe and share reproducible DNA materials and construction protocols." 13 citations. Extends the idea of scripted DNA documentation with a formal framework for reproducibility. Still focused on individual constructs, not variant libraries. **Low priority** for PoolParty manuscript, but worth being aware of as a recent entry in the "DNA construction as code" space.

- **Zulkower (2021)** — "Computer-Aided Design and Pre-validation of Large Batches of DNA Assemblies." 4 citations. Could be tangentially relevant if it addresses batch/library-level design, but this appears to be about validating many individual assembly designs (e.g., for DNA foundries), not about variant library design for functional assays. **Low priority.**

- **Pereira et al. (2016)** — "Yeast Pathway Kit." 17 citations. Extension of pydna for metabolic pathway assembly. Not relevant to variant library design. **No action needed.**

No citing paper poses a novelty concern or fills a gap in the PoolParty manuscript's literature review.

---
**Status:** Keep citation as-is. Accurately categorized. No novelty concerns. No changes to manuscript needed.
