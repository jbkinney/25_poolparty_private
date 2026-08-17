# Ledidi — evidence memo (tool survey for PoolParty referee response)

**Tool:** Ledidi (`ledidi`), Jacob Schreiber et al.
**One-line:** Gradient-based optimizer that turns any frozen, differentiable PyTorch sequence-to-function model into a *sequence editor*: it designs a compact set of substitutions to a single template sequence so the oracle's prediction matches a desired output.
**Not a library-design tool** in the PoolParty sense; it is a single-sequence, model-driven design method.

## Sources consulted

| Kind | Reference | Date accessed |
|---|---|---|
| PDF (preprint) | `papers/Schreiber2025_ledidi.pdf` — Schreiber, Lorbeer, Heinzl, Lu, Stark, Noble. "Programmatic design and editing of cis-regulatory elements." bioRxiv 2025.04.22.650035, posted 2025-04-23 (40 pp., extracted with PyMuPDF) | 2026-08-10 |
| Prior analysis | `prior_analyses/Schreiber2025_ledidi_analysis.md` | 2026-08-10 |
| Repo | https://github.com/jmschrei/ledidi (README.md, LICENSE, pyproject.toml, `ledidi/*.py`, `docs/*`) via raw.githubusercontent + GitHub API | 2026-08-10 |
| Docs | https://ledidi.readthedocs.io/en/latest/ (live); `docs/index.rst`, `getting_started.rst`, `input_output.rst`, `parameters.rst`, `faq.rst`, `whats_new.rst`, 9 tutorial notebooks | 2026-08-10 |
| PyPI | https://pypi.org/project/ledidi/ (JSON API) | 2026-08-10 |
| Model zoo | https://zenodo.org/records/14604495 (BPNet oracles used by paper + tutorials) | referenced, not downloaded |

Prior-analysis check: the prior notes are broadly **correct** (single-sequence ML-guided editing, not library design; affinity catalogs are repeated independent optimizations; no design cards/metadata). Two refinements: (1) the shipped `ledidi()` API *does* return many sequences per call — a batch (`batch_size`, default 16), `n_samples` (e.g. 10,000), `n_repeats`, and an affinity-catalog axis — so "output: one or a few sequences" understates it; the sequences are still a bare tensor with no identity or metadata. (2) Ledidi is Apache 2.0 on GitHub, actively developed as of June 2026 (the prior note said nothing about maintenance).

## What actually ships (repo/API)

Package layout: `ledidi/{__init__,ledidi,losses,plot,pruning,wrappers}.py`. Public surface:

- `ledidi(model, X, y_bar, n_repeats=1, n_samples=None, return_designer=False, return_history=False, device='cuda', random_state=None, **kwargs)` — the one-call entry point.
- `Ledidi(model, shape, target, input_loss=L1Loss(sum), output_loss=MSELoss(), tau=1, l=0.1, batch_size=16, max_iter=1000, early_stopping_iter=100, lr=1.0, input_mask=None, initial_weights=None, eps=1e-4, ...)` — the `torch.nn.Module` designer.
- `ledidi.wrappers.DesignWrapper([model_a, model_b])` — concatenates several oracles' predictions.
- `ledidi.losses.MinGap(in_mask)` — on-target vs off-target gap loss (cell-type-specific design).
- `ledidi.pruning.greedy_pruning(model, X, X_hat, threshold=1, target=None)`.
- `ledidi.plot.{plot_loss, plot_history, plot_edits}` (plot_edits built on `tangermeme.plot.plot_logo`).

I/O contract (docs/input_output.rst, verbatim): *"Sequences — both the template `X` you pass in and the designs `X_hat` you get back — are **one-hot encoded** tensors of shape `(batch, n_channels, length)` and dtype `torch.float32`."* and *"The template `X` passed to `ledidi` has a batch dimension of **1**."* There is **no** sequence-record, FASTA, or table object anywhere in the API — designs come back as a raw tensor.

Multi-sequence output does exist, but only as tensor axes:
> *"For an **affinity catalog**, pass a *list* of such tensors. The design is repeated once per element and the results are stacked, adding a leading catalog dimension to `X_hat` … `X_hat.shape == (4, batch_size, n_channels, length)`"* (docs/input_output.rst)
> Docstring of `ledidi()`: *"y: torch.Tensor, shape=(*ny, *n_repeats, n_samples, n_channels, length)"*

## Capability-by-capability evidence

### Block A — library specification

- **library_as_object = partial.** A single call returns a *batch* of designed sequences and can add catalog/repeat axes (`n_samples=10000` in one call; README: *"X_hat = ledidi(model, X, y_bar, n_samples=10000)  # Runs almost as fast as the default!"*). But the returned object is a bare `torch.Tensor`; there is no library abstraction, no per-sequence identity, no writer to FASTA/CSV, and heterogeneous design tasks must be scripted by the user. Every design starts from exactly one template (`_validate_input(X, "X", shape=(1, -1, -1), ohe=True)`).
- **dag_chaining = no.** Composition in Ledidi is over *oracles*, not over design steps: `DesignWrapper` *"concatenates their predictions"* (README). There is no pipeline/graph object, no operation chaining, no nesting of design specs. Roadmap item 3 ("Interoperability and end-to-end pipelines … so that one can go from raw data all the way to finished designs using command-line tools alone") is explicitly future work.
- **lazy_evaluation = partial.** `ledidi()` returns fully materialized tensors; nothing is deferred at the specification level (there is no specification). However, the fitted designer is a callable sampler that emits new designs on demand: *"once the Ledidi weight matrix is learned, edited sequences can be sampled extremely quickly using the forward function"* (README); Tutorial 8 is entirely about `return_designer=True` and calling `designer(X)` repeatedly. That is on-demand *generation*, not lazy evaluation of a library graph.
- **mixed_mutagenesis_one_pool = no.** Ledidi performs substitutions only, all driven by one objective; there are no mutagenesis "types" to mix, and no way to state "exhaustive singles + sampled random + WT replicates" in one spec. Insertions/deletions are a stated future item: roadmap 5, *"**Insertions and deletions.** Extending Ledidi beyond substitutions so that designs can change the length of a sequence, not just its content."*
- **combinatorial_motif_place = no.** The closest features are hard masks (`input_mask`), soft priors (`initial_weights`), and in-painting (zeroed columns that Ledidi fills in). Tutorial 3 (In-Painting): *"sections of the sequence are removed entirely and Ledidi is tasked with 'filling in' the region."* The motif content is *optimized*, not enumerated; there is no API for placing a set of motifs across positions/orientations/permutations.
- **barcode_generation = no.** No occurrence of "barcode" anywhere in the repo (README, docs, source). No edit-distance or GC-constrained tag generation, no concatenation of barcodes onto designs.
- **per_sequence_provenance = no.** Output is an unlabeled tensor. To find out what changed, the README tells you to diff it yourself: *"To see exactly which positions changed and how, compare the original sequence to a designed one"* followed by a manual `torch.where` loop. `return_history=True` gives per-*iteration* optimization statistics (input/output/total loss, where edits were proposed), not per-sequence provenance records. The README roadmap lists this as a cross-cutting *future* thread: *"**reproducibility and provenance** — capturing the full recipe behind every design so it can be audited and regenerated."*
- **automatic_naming = no.** Sequences are tensor slices; nothing in the API assigns names or IDs.
- **design_visualization = yes (designs, not a design graph).** `ledidi.plot.plot_edits` renders *"attribution values as a series of stacked logo tracks, with characters colored if they are edited with respect to the initial sequence. The original sequence is prepended as the first track."* `plot_history` plots *"the position of each proposed edit over the course of the optimization"*, `plot_loss` the loss curves. There is no graph/DAG view because there is no design graph; attributions must be computed by the user (tangermeme).

### Block B — assay coverage

- **assay_dms = no.** No codon-aware or amino-acid-level mutagenesis, no exhaustive variant enumeration, no coding-sequence handling. The README notes only that the *method* is modality-agnostic: *"one can also apply Ledidi out-of-the-box to RNA or protein models (or really, to any model with a sequence of categorical inputs)"* — i.e. protein oracles are usable, but nothing DMS-specific ships and no protein example is documented.
- **assay_mpra = partial.** Ledidi designs regulatory sequences *against MPRA/STARR-seq oracles* — Malinois (MPRA, K562) and DeepSTARR appear among the 14 design settings (paper Fig 2A, Supplementary Table 1; paper §3.1: *"such as TF binding, chromatin accessibility, transcription, or massively parallel reporter assay (MPRA) activity"*), and §3.4 validates designed 8-mer affinity catalogs against Reiter et al. STARR-seq data in *Drosophila* enhancers. But nothing about MPRA *library construction* is present: no barcodes, no adapters/primer sites, no oligo-length or synthesis handling, no replicate/control structure, no output file for ordering.
- **assay_insilico = partial.** Every Ledidi design is an in-silico interrogation of a genomic AI model, and affinity catalogs are explicitly framed as a probe of the model's learned code (paper §3.3, Fig 4E: *"an affinity catalog for GATA2 binding reveals a sophisticated usage of a learned cis-regulatory code"*; abstract: catalogs *"provide a resource for understanding the full range of potential sequence edits and their relative strengths"*). Tutorial 7 round-trips designs through independent held-out models (Enformer, an ATAC BPNet) and FIMO. But the sequence sets are produced by repeated *optimization*, not by a declarative specification of systematic in-silico perturbations; ISM/marginalization tooling lives in the sibling package tangermeme, not in ledidi.

### Block C — genomics integration

- **genome_coordinates = no.** The API accepts only one-hot tensors; coordinate/FASTA handling is explicitly delegated: *"The sibling library tangermeme provides `one_hot_encode` and `characters` utilities that handle this — as well as reading FASTA and genomic loci — for you"* (docs/input_output.rst). Paper examples (e.g. the SMYD3 promoter) extract sequence upstream of Ledidi.
- **transcript_models = no.** No GTF/GFF anywhere in repo or docs. Tutorial 2 discusses *"not editing … transcribed regions"* and gene bodies, but the user supplies a boolean position mask by hand; the paper likewise says *"we masked out the first 1,057bp of sequence"* manually.
- **exon_intron_split_codons = no.** No exon/intron or codon machinery of any kind in the source.
- **hgvs_input = no.** No HGVS parsing; input is a tensor.
- **vcf_vep_output = no.** No VCF/BED/GFF writers; the only output is `X_hat` (torch tensor) plus optional history dicts and matplotlib axes.
- **consequence_annotation = no.** No molecular-consequence logic; "consequence" in the paper means predicted regulatory activity from an oracle, not variant effect classification.

### Block D — physical construction

- **primer_design = no.** No primers, oligos, or wet-lab protocol output anywhere in repo/docs/paper.
- **codon_optimization = no.** Named only as future work: roadmap 4, *"Designing for several properties at once while respecting hard constraints such as GC content, restriction sites, codon usage, and the absence of off-target effects."*
- **synthesis_constraints = no.** Same roadmap sentence establishes GC content and restriction sites are *not* yet handled. `initial_weights` priors and `input_mask` constrain *where/what* edits may occur, but no synthesizability check is applied to the resulting sequences.

### Block E — engineering

- **interface:** Python API only (PyTorch). `pyproject.toml` declares no `[project.scripts]`/console entry points; there is no CLI, web service, or GUI. A CLI is roadmap item 3 (future), agentic interfaces roadmap item 8.
- **license:** Apache 2.0 — `LICENSE` file is the Apache License 2.0 text; README badge and closing section say Apache 2.0; GitHub API reports `spdx_id: Apache-2.0`. **Discrepancy:** PyPI metadata for v2.1.0 reports "MIT License" (stale classifier); the repository text is authoritative.
- **maintained:** Actively maintained. Last commit on `master` **2026-06-23** ("UPDATE plot_edits to render via tangermeme.plot_logo (#19)"); repo `pushed_at` 2026-06-23, `updated_at` 2026-07-30; not archived; 110 stars; 0 open issues; CI test workflow present with a 100% coverage badge. Latest PyPI release **2.1.0 (2025-04-24)**; only tagged GitHub release is v2.0.0 (2025-01-06); `docs/whats_new.rst` has an unreleased "Version 2.2.0" section, so master is ahead of PyPI.

## Documented examples / vignettes (for possible PoolParty reproduction)

Nine tutorial notebooks in `docs/tutorials/`, all rendered on readthedocs. Tutorials 0 and 8 run on CPU with no downloads; the rest use BPNet/Enformer/Malinois oracles from https://zenodo.org/records/14604495.

0. **A First Design** — toy parameter-free AP-1 (`TGACTCA`) convolution oracle on a random 50 bp sequence; Ledidi inserts the motif in ~2 edits. Fully reproducible with no data.
1. **Design of Protein Binding Sites** — BPNet GATA2 (K562) oracle; design edits creating GATA2 binding; inspect proposed edits.
2. **Constraints and Priors** — hard `input_mask` (e.g. protect gene body / TATAA box) and soft `initial_weights` priors.
3. **In-Painting** — zero out a span and let Ledidi fill it in ("free" edits); used to insert/move/replace motifs.
4. **Multiple Models** — `DesignWrapper` over several oracles/tasks: cell-type-specific, cross-modal, and "design-property" objectives.
5. **Affinity Catalogs** — pass a *list* of `y_bar` targets to sweep design strength; GATA2 BPNet example.
6. **Custom Loss Functions** — custom output losses, incl. using BPNet base-resolution profiles to control *where* a site is placed.
7. **Validating Your Designs** — attributions, FIMO motif counting, and round-tripping designs through held-out models (Enformer, ATAC BPNet).
8. **The Ledidi Object** — fit the weight matrix once, then sample many designs cheaply (`return_designer`, `n_samples`).

Paper-level case studies (Methods §2.4–2.6, Results §3): 14 design settings across 10 pre-trained models (BPNet GATA2/E2F6, Beluga CTCF/MAX, BPNet-ATAC, ChromBPNet, Enformer DNase/CAGE, Sei, ProCapNet, Puffin, Borzoi, Malinois K562, DeepSTARR); percentile-based design from the 50th to 99.5th percentile on chr12/hg38; affinity catalogs at the SMYD3 promoter (GATA2 binding, accessibility, transcription initiation); 8-mer affinity catalogs at seven sites in the *Drosophila* ced6 and ZnT63C enhancers validated against STARR-seq (r = 0.893–0.925); 8-mer catalogs in yeast 80 bp promoters validated against DREAM/GPRA data via LegNet.

## Notable capabilities not covered by the row list

1. **Model-in-the-loop objective-driven design** (the core of the tool) — no row captures "generate sequences by optimizing against a predictive oracle".
2. **Multi-objective design across several models** (`DesignWrapper`, `MinGap` loss for on-target vs off-target, e.g. cell-type-specific enhancers).
3. **Edit-count minimization + greedy pruning** — an explicit input loss penalizing number of edits, plus post-hoc `greedy_pruning` to revert edits that change the prediction by less than a threshold. Relevant to genome-editing feasibility, absent from all library-design tools.
4. **Hard masks and soft priors over edit space** (`input_mask`, `initial_weights` with `-inf` to forbid specific characters at specific positions) and **in-painting** (free edits in a span).
5. **Affinity catalog** — a swept series of design targets producing a graded set of sequences; the closest thing here to a "library", but produced by independent optimizations.
6. **Reproducible sampling without mutating global RNG** (`random_state` with a private `torch.Generator`).
7. **GPU-native, modality-agnostic** (DNA/RNA/protein/any categorical alphabet, any differentiable PyTorch oracle).

## Stated limitations (author's own words)

- Substitutions only; **no insertions or deletions** yet (roadmap 5).
- **No biological/synthesis constraints yet** — GC content, restriction sites, codon usage, off-target effects are roadmap 4.
- **No CLI / end-to-end pipeline yet** (roadmap 3); **no provenance capture yet** (cross-cutting roadmap thread).
- Designs in a returned batch are **correlated**: *"Every sequence in the returned batch is sampled from the same learned weight matrix … If a motif could have been inserted at several locations, Ledidi commits to one"* (FAQ). Independent designs require `n_repeats` or new `random_state`.
- **Oracle exploitation risk**: *"gradient-based design can produce sequences that score well on the oracle yet are not biologically meaningful"* (Tutorial 7); the paper's affinity catalogs are *"most accurate within the dynamic range of the training data"* (§3.4).
- Requires a differentiable pre-trained model and, practically, a GPU (`device='cuda'` default).
- Self-described as **"research software under active development"** (README roadmap).

## Availability today (2026-08-10)

Installable and runnable: `pip install ledidi` (v2.1.0, 2025-04-24; requires Python ≥3.10, torch ≥2.0, tangermeme ≥1.3.0, numpy, matplotlib). GitHub repo alive and not archived, last commit 2026-06-23, CI configured, 0 open issues, 110 stars. Docs site https://ledidi.readthedocs.io/en/latest/ is live. Oracle models used by the tutorials are on Zenodo (record 14604495). No web service or CLI exists — it is a library. Paper remains a bioRxiv preprint (posted 2025-04-23), CC-BY 4.0; the preprint has no formal "Code availability" section, but the README/docs point to the repo and the preprint DOI.

## Confidence notes

Highest confidence: everything read directly from the shipped source, docs, `pyproject.toml`, LICENSE, and GitHub/PyPI metadata (Blocks C, D, E, and the "no" values in Block A). Lower confidence / judgement calls: `library_as_object` and `lazy_evaluation` ("partial" rather than "no" — batch/`n_samples`/catalog axes and the resamplable fitted designer are real, but there is no library object or deferred spec); `assay_mpra` and `assay_insilico` ("partial" — MPRA/STARR oracles are central to the paper and designs are used to probe models, but no MPRA construct assembly and no declarative in-silico perturbation framework); `design_visualization` ("yes" on the strength of `plot_edits` highlighting edited bases, though there is no library- or graph-level view). Not verified by execution — no installation was performed, per instructions.
