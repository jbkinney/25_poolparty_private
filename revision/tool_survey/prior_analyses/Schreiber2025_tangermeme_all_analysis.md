# Schreiber (2025) — tangermeme: A toolkit for understanding cis-regulatory logic using deep learning models

## BibTeX Entry

```bibtex
@article{Schreiber2025nd,
  year     = {2025},
  title    = {tangermeme: A toolkit for understanding cis-regulatory logic using deep learning models},
  author   = {Schreiber, Jacob},
  journal  = {{bioRxiv}},
  doi      = {10.1101/2025.08.08.669296},
  pages    = {2025.08.08.669296},
  mynotes  = {Python toolkit for interpreting genomic deep learning models via sequence manipulations (marginalization, ablation, saturation mutagenesis, variant effect prediction) paired with model operations (predictions, DeepLIFT/SHAP attributions). Includes novel seqlet-calling algorithms. Should be cited in the Introduction as an example of a general-purpose sequence toolkit that provides sequence manipulation operations but has no notion of a variant library as a structured object. Could alternatively or additionally be cited alongside the genomic AI references as a tool for in silico probing of models, complementary to the library design role that PoolParty fills.}
}
```

## 1. Summary

tangermeme is a Python toolkit for interpreting trained genomic deep learning models. It implements "everything-but-the-model": optimized sequence manipulations (marginalization, ablation, variant effect prediction, saturation mutagenesis), model operations (predictions, DeepLIFT/SHAP attributions), and novel algorithms for identifying "seqlets" (contiguous spans of high-attribution nucleotides). Seqlets are automatically called, annotated against motif databases, and counted to derive statistics on motif usage, co-occurrence, and spacing. tangermeme's key design principle is the separation of sequence manipulation from model operations, enabling richer analyses (e.g., computing attributions on marginalized sequences rather than just predictions). It is significantly faster than alternatives at core operations like one-hot encoding. The paper demonstrates tangermeme by comparing the learned logics of multiple models (BPNet, Beluga, ChromBPNet, ProCapNet) at the PLD6 promoter.

## 2. Relevance to PoolParty Manuscript

tangermeme and PoolParty overlap in the space of computational sequence manipulation but serve fundamentally different purposes:

- **tangermeme** manipulates sequences to *interpret* trained deep learning models (in silico perturbation experiments). Its sequence operations are model-focused: substitute a motif, ablate a region, compute variant effects — all in the service of understanding what a model has learned.
- **PoolParty** designs *libraries* of variant sequences for experimental or computational use, with structured metadata (design cards, names, styling) and a composable DAG architecture.

tangermeme has no concept of a variant library as a structured object, no DAG-based design workflow, no design cards, no lazy evaluation, and no library-level combinatorics. It operates on individual sequences or batches of sequences for model interpretation.

However, there is a meaningful connection: tangermeme's in silico perturbation experiments are exactly the kind of analyses that PoolParty-designed libraries could feed into. The SpliceAI surrogate modeling example in the PoolParty manuscript is conceptually similar to what tangermeme facilitates. The tools are complementary rather than competing.

The manuscript currently categorizes tangermeme among "general-purpose sequence toolkits" alongside Biopython, pydna, and SeqPro. This is defensible — tangermeme does provide sequence manipulation operations — though tangermeme is more specifically a genomic AI interpretation toolkit. The manuscript's point holds either way: tangermeme can manipulate individual sequences but lacks library-level abstractions.

## 3. Citation Recommendations

**Should be cited** in the Introduction (paragraph 3) as an example of a toolkit that provides sequence manipulation operations but has no notion of a variant library as a structured object. The current placement is appropriate.

**Optional additional citation**: tangermeme could also be mentioned in the Discussion or in the SpliceAI section as an example of a complementary tool for in silico model interpretation, illustrating how PoolParty-designed libraries could be used with interpretation toolkits like tangermeme. This would strengthen the narrative that PoolParty fills a gap in the genomic AI workflow (library design) that existing tools like tangermeme do not address.

## 4. Novelty Assessment

**No threat to novelty claims.** tangermeme operates in a different problem domain (model interpretation) and at a different level of abstraction (individual sequence perturbations, not library design). Specific distinctions:

| Feature | tangermeme | PoolParty |
|---------|-----------|-----------|
| Purpose | Interpret genomic DL models | Design variant sequence libraries |
| Core operations | Marginalization, ablation, variant effect, saturation mutagenesis | 20+ Operations for library construction |
| Output | Model predictions/attributions on perturbed sequences | Libraries of thousands–millions of sequences |
| Library abstraction | None | Central (Pools, DAGs) |
| Design cards | None | Core feature |
| Lazy evaluation | No | Yes |
| Sequence metadata | None | Names, styling, design cards |

No special precautions needed. The tools are complementary.

## 5. References to Investigate

- **Schreiber et al. (2025b)** — "Programmatic design and editing of cis-regulatory elements." bioRxiv 2025.04.22.650035. Same first author. This paper is about *designing* regulatory sequences programmatically, which could be more directly relevant to PoolParty's domain. **Medium priority** — worth checking whether it addresses library-level design or operates on individual sequences.

- **inMOTIFin: Ferenc et al. (2025)** — "A lightweight end-to-end simulation software for regulatory sequences." A simulator for generating synthetic regulatory sequences with controlled motif placement. Cited by tangermeme. Could be relevant as another tool in the sequence generation space. **Low priority** — appears to be a training data simulator, not a library design tool.

## 6. Citing Papers to Investigate

Searched via Semantic Scholar (DOI: 10.1101/2025.08.08.669296). 8 citing papers total (recent preprint, so limited coverage).

- **inMOTIFin (Ferenc 2025)** — "A lightweight end-to-end simulation software for regulatory sequences." Generates synthetic regulatory DNA with precise motif control. Potentially relevant as another sequence generation tool, but aimed at training data simulation rather than experimental library design. **Low priority.**

- **GaugeFixer (Martí-Gómez 2025)** — Resolves parameter non-identifiability in sequence-function models. Not relevant to library design.

No citing paper poses a novelty concern or needs to be added to the PoolParty manuscript.

---
**Status:** Keep citation as-is. Correctly categorized as a sequence toolkit lacking library-level abstractions. No novelty concerns. Consider optional additional citation in Discussion to position PoolParty as complementary to model interpretation tools.
