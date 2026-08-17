# Ledidi — v2 capability record (PoolParty referee response)

**Tool:** Ledidi (`ledidi`) · **slug:** `ledidi` · **Citation key:** `Schreiber2025ledidi` · **Tier:** 3
**Row set:** v2 (`ROWS_v2.md`) · **Scored:** 2026-08-10
**Basis:** `final/ledidi.md` (self-contained, post-adversarial-review), re-read in full. No new
extraction was performed; no web re-check was needed — every v2 row was answerable from the memo.
Nothing was installed or executed.

**One-line:** Gradient-based optimizer that turns any frozen, differentiable PyTorch
sequence-to-function model into a *sequence editor*: it learns a per-position/per-character weight
matrix (Gumbel-softmax) so that a small set of substitutions to **one** template sequence drives the
oracle's prediction to a user-specified target. It is a single-template, model-driven design method
that happens to emit batches — not a library-design tool in the PoolParty sense.

**Scoring stance.** Ledidi is a purpose-built tool with a genuine, well-documented capability
surface. Where the v2 wording is broader than v1's, Ledidi is credited for what it actually does
(sampling from a fitted designer, multi-substitution designs, catalog stacking, model-in-the-loop).
Where a row asks for an enumerative or declarative *library* operation, Ledidi is "no" and the
evidence says exactly why, so the row cannot be misread as "Ledidi cannot do X at all".

---

## Block A — library specification

| Key | v2 value | v1 value | Note |
|---|---|---|---|
| `library_first_class_object` | **partial** | (`library_as_object` = partial) | split |
| `composable_operations` | **no** | (`dag_chaining` = no) | rename, value carried |
| `lazy_generation` | **partial** | (`lazy_evaluation` = partial) | rename, value carried and *firmer* |
| `library_algebra` | **partial** | (`library_as_object` = partial) | split |
| `exhaustive_single_scans` | **no** | (`mixed_mutagenesis_one_pool` = no) | split |
| `sampled_random_mutagenesis` | **partial** | (`mixed_mutagenesis_one_pool` = no) | split — **upgraded** |
| `higher_order_combinatorial` | **partial** | (`mixed_mutagenesis_one_pool` = no) | split — **upgraded** |
| `heterogeneous_components_one_library` | **no** | (`mixed_mutagenesis_one_pool` = no) | split |
| `combinatorial_motif_place` | **no** | no | unchanged |
| `barcode_generation` | **no** | no | unchanged |
| `per_sequence_provenance` | **no** | no | unchanged |
| `automatic_naming` | **no** | no | unchanged |
| `design_visualization` | **partial** | partial | unchanged |

### `library_first_class_object` = **partial**

Ledidi is **not** a file-writing tool — everything is in-memory and can be held, inspected,
transformed and passed onward. `ledidi()` returns a `torch.Tensor` of designs (rank-5 when catalog /
repeat / sample axes are used: `(*ny, *n_repeats, n_samples, n_channels, length)`), which the user
can feed to another oracle, to `ledidi.pruning.greedy_pruning(model, X, X_hat, ...)`, or to
`ledidi.plot.plot_edits`. `return_designer=True` additionally returns a fitted `Ledidi`
`torch.nn.Module` — a genuine first-class object that can be stored and re-sampled (Tutorial 8).
What is missing is the *library* half of the row: the returned object is a bare unlabeled tensor
with no per-element identity, no metadata, no naming, and no library abstraction or writer
(no FASTA/CSV/table anywhere in the API); the template is hard-restricted to a **single** sequence
(`ledidi.py:204`, `_validate_input(X,"X",shape=(1,-1,-1),ohe=True,allow_N=True)`).
`partial` is the generous-but-defensible reading and is a *weaker* partial than tangermeme's.
**Source:** `ledidi/ledidi.py:204, 212–252`, `ledidi/pruning.py`, `ledidi/plot.py`;
`docs/input_output.rst`; README quick-start; Tutorial 8.

### `composable_operations` = **no**

The only composition Ledidi ships is over *oracles*, not over design steps: `DesignWrapper`
*"concatenates their predictions"* (a `torch.cat` of model outputs). There is no pipeline object,
graph, step composition, or nesting of design specifications in any of the 6 modules; chaining one
design into the next is a user-written loop. Pre-empting the author objection: paper §2.1 says
masks, in-painting and priors *"can be seamlessly combined with each other"* — that is **constraint
composition inside one design**, not composition of design steps, and it is credited under
"notable capabilities" rather than here. README roadmap item 3 puts pipelines in the future:
*"Interoperability and end-to-end pipelines … so that one can go from raw data all the way to
finished designs using command-line tools alone."*
**Source:** `ledidi/wrappers.py`; README roadmap item 3; PDF §2.1.

### `lazy_generation` = **partial**

The v2 wording ("produced **on demand** rather than fully materialized") fits Ledidi better than
v1's "lazy evaluation" did, so `partial` is firmer here than it was in v1. `ledidi()` by default
returns fully materialized tensors and there is no deferred specification object. But
`return_designer=True` yields a fitted `Ledidi` module that is a real on-demand sampler: README —
*"once the Ledidi weight matrix is learned, edited sequences can be sampled extremely quickly using
the forward function"*; Tutorial 8 — *"Fitting is the expensive part … drawing a fresh batch of
designs afterward is a single forward pass."* Not `yes` because each materialization yields
*different* stochastic sequences (resampling from the fitted weight matrix), so there is no fixed
sequence set being lazily realized, and the default path materializes everything.
**Source:** `ledidi/ledidi.py:467` (`forward`); README; Tutorial 8.

### `library_algebra` = **partial**

Repetition and stacking of design *runs* are first-class parameters inside the tool, not user
scripting: `n_repeats` (independent re-designs), `n_samples` (draw N designs from the fitted weight
matrix; realized as `n_iter = n_samples // batch_size + 1` then `torch.cat`, `ledidi.py:212–252`),
and the **affinity catalog** — `docs/input_output.rst`: *"For an affinity catalog, pass a *list* of
such tensors. The design is repeated once per element and the results are stacked, adding a leading
catalog dimension to `X_hat`."* That is repeat + sample + stack as tool operations. What does not
exist is algebra over *libraries*: no operation combines two independently specified design sets,
no concat/union/subsample of libraries as objects, and the stacked axes carry no identity, so any
real combination (different templates, different objectives, WT spike-ins) is a user-written
`torch.cat` outside the API. Hence `partial`, not `yes`.
**Source:** `docs/input_output.rst` (catalog stacking); `ledidi/ledidi.py:212–252`; `ledidi()`
signature (`n_repeats`, `n_samples`).

### `exhaustive_single_scans` = **no**

Ledidi never enumerates. It optimizes a weight matrix by gradient descent and samples from it; no
saturation-mutagenesis or all-positions × all-characters operation exists in the 6 modules (that
lives in the sibling library tangermeme). Enumeration appears only as the *baseline Ledidi is
compared against*: paper Fig 2G benchmarks Ledidi as *"almost an order of magnitude faster than
greedy substitution and several orders of magnitude faster than a comprehensive motif implantation
method"*, where those baselines implant consensus sequences *"at each possible position"*. Also,
deletions are impossible at all — indels are README roadmap item 5 (*"Insertions and deletions.
Extending Ledidi beyond substitutions…"*).
**Source:** module listing (no enumeration API); PDF §3 / Fig 2G; README roadmap item 5.

### `sampled_random_mutagenesis` = **partial**

Stochastic sampling of variants is native and industrial-scale: designs are drawn from a learned
Gumbel-softmax weight matrix, `batch_size` designs per iteration and `n_samples` designs per call —
README: *"X_hat = ledidi(model, X, y_bar, n_samples=10000)  # Runs almost as fast as the default!"* —
with `random_state` backed by a private `torch.Generator` for reproducible draws. FAQ, verbatim:
*"Every sequence in the returned batch is sampled from the same learned weight matrix."* So this is
genuinely sampling, not enumeration. It is `partial` rather than `yes` because the sampling
distribution is **learned from an objective**, not a user-specified random or per-base mutation
rate: there is no "mutate at rate p" / "k random substitutions" mode, the number of edits is only
steered indirectly through the input-loss weight `l` (mean ~16.54 edits in the paper's chr12
benchmark), and draws within one batch are correlated by construction (FAQ: *"If a motif could have
been inserted at several locations, Ledidi commits to one"*).
**Source:** README quick-start (`n_samples=10000`); `ledidi/ledidi.py:212–252, 467`; `docs/faq.rst`;
`docs/parameters.rst` (`l`, `tau`, `random_state`); PDF §3 (mean edit count).

### `higher_order_combinatorial` = **partial**

Every Ledidi design is a **multi-substitution** variant — combinations of mutations co-occurring in
one sequence are the tool's native output, not a special mode (paper: percentile-based designs on
chr12/hg38 at ~16.54 edits mean; Tutorial 2 knocks out a chosen GATA site with several coordinated
edits; forcing a whole motif in is a block of simultaneous substitutions). Across `n_samples` the
user gets many *distinct* multi-mutant combinations, and `greedy_pruning` then searches subsets of
an edit set by reverting edits whose effect falls below a threshold. It is `partial`, not `yes`,
because nothing **enumerates or specifies** combinations: there is no pairwise/k-wise mutator, no
"all pairs of positions" operation, no control over combination order or identity — which
substitutions co-occur is chosen by the optimizer, and the resulting combinations are unlabeled.
**Source:** PDF §3 (edit counts); Tutorials 2–3; `ledidi/pruning.py` (`greedy_pruning`);
module listing (no combinatorial mutator).

### `heterogeneous_components_one_library` = **no**

One `ledidi()` call performs substitutions only, over one fixed-length one-hot template, under one
edit-space specification. There is no way to express structurally different component types —
"exhaustive singles + sampled higher-order + WT replicates + scrambles" — in a single design
specification, because there are no component or mutagenesis *types* to mix. Output length always
equals template length, so no component can differ structurally. The nearest thing, the affinity
catalog, varies only the numeric *target* across an otherwise identical design (`y_bar` list), and
`DesignWrapper` varies the *objective*, not the component structure.
**Source:** `ledidi/ledidi.py:204` (single template, fixed shape); `docs/input_output.rst`
(output shape == input shape); README roadmap item 5 (indels future).

### `combinatorial_motif_place` = **no** — *no enumeration; placement itself IS supported*

Ledidi **can** force a specified motif into a user-chosen position: Tutorial 2, *"Adding in
Motifs"* — *"we can also use careful usage of a mask to force edits to happen -- even editing in an
entire motif. … If we set all the characters to `-inf` except for some other character, we could
force that edit to be included in the design."* It can also **relocate** motifs: Tutorial 3,
*"Moving Motifs (CTCF + Transcription)"* — in-paint the old site, in-paint a landing pad, hold the
oracle output fixed. Placements **can** be swept, but only in a user-written loop: Tutorial 3,
*"…moving the 20bp region of the BPNet motif downstream between 1 and 800 bp. Every 10th basepair…"*.
What does not exist is any operation enumerating motif sets × positions × orientations × spacings ×
permutations as a library specification. The cell is `no` strictly on absence of combinatorial
enumeration.
**Source:** Tutorial 2 ("Adding in Motifs"); Tutorial 3 ("Moving Motifs", positional sweep);
`ledidi/ledidi.py` (`input_mask`, `initial_weights`).

### `barcode_generation` = **no**

Zero occurrences of `barcode`, `oligo`, or `adapter` across all 6 modules, all 7 `.rst` pages, the
README and the tutorial notebooks (grepped independently by extractor and reviewer). No
edit-distance or GC-constrained tag generation, no attachment of tags to designs — and output length
always equals template length, so a tag could not be appended in the first place.
**Source:** repo-wide grep; `docs/input_output.rst` (output shape == input shape).

### `per_sequence_provenance` = **no**

The output is an unlabeled tensor and **no per-sequence record travels with it**. README tells users
to diff by hand: *"To see exactly which positions changed and how, compare the original sequence to
a designed one"*, followed by a manual `torch.where` loop. Stated precisely so the row cannot be
attacked: `return_history=True` does record `history['edits'].append(torch.where(X_hat != X_))`
(`ledidi.py:603`) — *(batch index, channel, position)* triples **per optimization iteration**, so a
batch-member index does exist. But the returned `best_sequence` is cloned from whichever iteration
last improved the loss and **that iteration index is never recorded**, so nothing maps a returned
design back to its history rows. README roadmap lists *"reproducibility and provenance — capturing
the full recipe behind every design so it can be audited and regenerated"* as a **future**
cross-cutting thread.
**Source:** `ledidi/ledidi.py:603–606`; README (diff recipe + roadmap provenance thread).

### `automatic_naming` = **no**

Designs are slices of a `torch.float32` tensor; nothing in the API assigns names, IDs or labels. The
docs decode designs to strings with a hand-written local `decode()` helper.
**Source:** `docs/getting_started.rst`; `docs/input_output.rst`; all 6 modules.

### `design_visualization` = **partial**

Ledidi ships `plot_loss` (input/output loss curves), `plot_history` (*"the position of each proposed
edit over the course of the optimization"*) and `plot_edits` (*"attribution values as a series of
stacked logo tracks, with characters colored if they are edited with respect to the initial
sequence"*). `partial` not `yes` for three reasons: (1) **cross-tool consistency** — `plot_edits` is
a thin wrapper over `tangermeme.plot.plot_logo` (`plot.py:21` import, `plot.py:193` call) and this
survey scores tangermeme `partial` for a strictly richer plotting module; a tool cannot score above
the library it delegates rendering to. (2) **Scope** — it renders one design run at a time and
requires user-supplied attribution tensors; there is no library-level view, no design-set summary,
no design-graph view. (3) **Release reality** — in the pip-installable v2.1.0 (commit `7adfcd6453`)
`plot.py` imports `pandas` and `logomaker`, neither declared in that build's `install_requires`, and
`plot_loss` does not exist there at all.
**Source:** `ledidi/plot.py:17–21, 24, 85, 118, 193` (master); `ledidi/plot.py:4–8, 11, 44` and
`setup.py` at `7adfcd6453`; tangermeme row of this survey.

---

## Block B — assay coverage

### `assay_dms` = **no**

No codon-aware or amino-acid-level mutagenesis, no reading-frame logic, no exhaustive variant
enumeration, no coding-sequence handling, no protein example documented. The only protein-adjacent
statement is a modality claim about the *method*: README — *"one can also apply Ledidi
out-of-the-box to RNA or protein models (or really, to any model with a sequence of categorical
inputs)"*. Consistent with the survey's own bar: tangermeme scores `no` here despite shipping
`saturation_mutagenesis`.
**Source:** repo-wide grep (codon / amino-acid / reading-frame / exon); README.

### `assay_mpra` = **partial**

Ledidi designs regulatory sequences **against MPRA/STARR-seq oracles**: paper §3.1 lists *"TF
binding, chromatin accessibility, transcription, or massively parallel reporter assay (MPRA)
activity"* among design targets; Malinois (MPRA, K562) and DeepSTARR are among the 10 pre-trained
oracles across 14 design settings (Fig 2A, Suppl. Table 1); §3.4 validates designed 8-mer affinity
catalogs against Reiter et al. STARR-seq data in *Drosophila* enhancers (r = 0.893–0.925, p.25). But
**no MPRA construct assembly ships**: no barcode, adapter, primer, oligo or synthesis handling; no
oligo-length, replicate or control structure; no orderable output file.
**Source:** PDF §3.1, §3.4 (p.25), Fig 2A / Suppl. Table 1; repo-wide grep.

### `assay_insilico` = **partial**

Every Ledidi design is an in-silico interrogation of a genomic AI model, and affinity catalogs are
explicitly framed as probes of the learned code (paper §3.3 / Fig 4E: *"an affinity catalog for
GATA2 binding reveals a sophisticated usage of a learned cis-regulatory code"*; abstract: catalogs
*"provide a resource for understanding the full range of potential sequence edits and their relative
strengths"*). Paper §2.6 additionally uses in-silico marginalization (implanting a motif into 1,000
dinucleotide-shuffled backgrounds) as an evaluation. `partial` because Ledidi ships **no declarative
perturbation operations** — the bar for `yes` is tangermeme's `marginalize`,
`saturation_mutagenesis`, `variant_effect.*`, which live in tangermeme, not ledidi. Ledidi reaches
in-silico interrogation only through repeated optimization; Tutorial 7 is *design validation*, not a
perturbation framework.
**Source:** PDF §2.6, §3.3, Fig 4E, abstract; Tutorial 7; module listing (no perturbation API).

---

## Block C — genomics integration

### `genome_coordinates` = **no**

The API accepts only one-hot `(1, n_channels, length)` tensors; coordinate and FASTA handling is
explicitly delegated: `docs/input_output.rst` — *"The sibling library tangermeme provides
`one_hot_encode` and `characters` utilities that handle this — as well as reading FASTA and genomic
loci — for you."* Named so the row is not misread: tutorials **do** work from hg38 loci (e.g. the
SMYD3 promoter), but via `pyfaidx` + `tangermeme.utils.one_hot_encode` **outside** ledidi's API.
**Source:** `docs/input_output.rst`; `ledidi.py:204`; Tutorials 1–3.

### `transcript_models` = **no**

Zero GTF/GFF hits anywhere in repo, docs or notebooks. Gene-body protection is achieved only by a
**user-supplied boolean position mask**: Tutorial 2, *"Preventing Edits (Gene Bodies)"* (*"one may
wish to prevent edits at the TATAA box or within the coding region"*, a hand-built boolean mask);
paper §2.5 *"we masked out the first 1,057bp of sequence"*.
**Source:** repo-wide grep; Tutorial 2; PDF §2.5.

### `exon_intron_split_codons` = **no**

No exon, intron, reading-frame or codon concept anywhere in the 6 modules or the docs; the sequence
is an unannotated fixed-length categorical array.
**Source:** repo-wide grep over modules + `.rst` + notebooks.

### `vcf_vep_output` = **no**

Reading the full return path of `ledidi()`: outputs are `X_bar` (tensor), optional designer objects,
optional history dicts, plus matplotlib axes. No VCF, BED, GFF, FASTA or CSV writer exists in any
module.
**Source:** `ledidi/ledidi.py` (return path); `ledidi/plot.py`; module listing.

### `consequence_annotation` = **no**

No molecular-consequence classification anywhere. Every readout in the repo and paper is either an
oracle prediction or an edit count/position; "consequence" in the paper means predicted regulatory
activity, not variant-effect classification.
**Source:** repo-wide grep; PDF results sections.

*(v1's `hgvs_input = no` is dropped from v2 per `ROWS_v2.md`.)*

---

## Block D — adjacent / complementary

### `primer_design` = **no**

No primer, oligo or wet-lab protocol output in repo, docs, notebooks or paper. The paper frames the
wet-lab route as CRISPR base/prime editing but designs **no reagents**; its "experimental
validation" is a retrospective lookup against published STARR-seq and DREAM/GPRA data.
**Source:** repo-wide grep; PDF §3.4.

### `codon_optimization` = **no**

The word "codon" appears exactly once in the entire repo, in README roadmap item 4 as **future**
work: *"Designing for several properties at once while respecting hard constraints such as GC
content, restriction sites, codon usage, and the absence of off-target effects."*
**Source:** README roadmap item 4.

### `synthesis_constraints` = **no** — *but a real edit-composition constraint system exists*

There is **no** GC-content target, homopolymer/repeat check, restriction-site avoidance, Tm
calculation or post-hoc synthesizability check anywhere; all are named in roadmap item 4 as future
work. What does exist, and what an author would cite, is a constraint system over *edits*:
Tutorial 2, *"Preventing Edits (No Added Cs)"* — `initial_weights[:, 1, :] = float("-inf")` forbids
ever introducing a C, motivated by editing feasibility (*"you might want to prevent certain
nucleotides from being edited in, potentially due to issues with prime or base editing"*), down to
blocking a single specific substitution. That is an **edit-composition constraint, not a synthesis
constraint** — it constrains which edits the optimizer may propose, never whether the resulting
molecule is orderable.
**Source:** Tutorial 2; README roadmap item 4; `ledidi.py:430`.

### `degenerate_iupac_codons` = **no**

No IUPAC/degenerate-code design, expansion or compression anywhere: there is no codon concept at all
(the word "codon" occurs once, in the roadmap), no ambiguity alphabet, and the alphabet is whatever
one-hot channel set the oracle takes. The one adjacent feature is input tolerance, not design:
`_validate_input(X, "X", shape=(1,-1,-1), ohe=True, allow_N=True)` (`ledidi.py:204`) permits `N`
(an all-zero one-hot column) in the *template*; designs themselves are always concrete one-hot
characters sampled from the weight matrix.
**Source:** `ledidi/ledidi.py:204`; repo-wide grep (codon / IUPAC / degenerate); README roadmap
item 4.

### `negative_control_generation` = **no**

No scramble, shuffle, or reverse/complement control generator is exposed anywhere in the 6 modules
or documented as a feature. Two near-misses, both external to ledidi's API: the paper's **evaluation**
implants motifs into 1,000 **dinucleotide-shuffled** backgrounds (§2.6 marginalization — shuffling
performed with tangermeme, not ledidi), and Tutorial 7 counts FIMO hits for **decoy motifs** as a
validation comparison. Neither generates control *sequences* as part of a design.
**Source:** PDF §2.6; Tutorial 7; module listing / repo-wide grep.

### `ml_model_in_loop` = **yes**

This is the entire tool. Ledidi optimizes a Gumbel-softmax weight matrix by gradient descent through
a **frozen, differentiable, pre-trained PyTorch oracle** so that the oracle's prediction hits a
user-specified numeric target `y_bar`: `ledidi(model, X, y_bar, ...)`, `Ledidi(model, shape,
target=..., input_loss=..., output_loss=...)`. Multi-model objectives are first-class —
`wrappers.DesignWrapper([model_a, model_b])` concatenates several oracles' predictions and
`losses.MinGap(in_mask)` expresses on-target vs off-target gaps (e.g. cell-type-specific enhancers).
14 design settings across 10 pre-trained models are demonstrated (BPNet GATA2/E2F6, Beluga
CTCF/MAX, BPNet-ATAC, ChromBPNet, Enformer DNase/CAGE, Sei, ProCapNet, Puffin, Borzoi, Malinois
K562, DeepSTARR). Design is symmetric: knock-in and knock-down/knock-out (SMYD3 catalogs span −3.0
to +4.0). Requirement, stated as a limitation: a differentiable pre-trained model, and practically a
GPU (`device='cuda'` default).
**Source:** `ledidi/ledidi.py` (entry point + `Ledidi` module); `ledidi/wrappers.py`;
`ledidi/losses.py` (`MinGap`); PDF §2–§3, Fig 2A / Suppl. Table 1; Tutorials 1–6.

### `readout_analysis` = **no**

No analysis of sequencing readout from a built library: no FASTQ/count/amplicon handling, no
demultiplexing (there are no barcodes), no writer or reader of assay data anywhere in the 6 modules.
Ledidi's validation loop is entirely in-silico — attributions, FIMO motif-hit deltas against decoys,
and round-tripping designs through independently trained oracles (Tutorial 7) — and the paper's
"experimental validation" is a retrospective comparison of *predictions* to **already-published**
STARR-seq and DREAM/GPRA measurements (§3.4), not analysis of sequencing from a Ledidi-built
library. No Ledidi library has been synthesized.
**Source:** Tutorial 7; PDF §3.4; module listing / repo-wide grep.

---

## Block E — engineering and availability

| Key | Enum | Descriptor for the table |
|---|---|---|
| `interface` | **yes** | **Python API (library) only — no CLI, GUI or web service** |
| `license` | **yes** | **Apache-2.0** (shipped v2.1.0 `setup.py` says `license='MIT License'` — repo `LICENSE` is authoritative) |
| `installable_today` | **yes** | **Yes — `pip install ledidi` → v2.1.0 (2025-04-24); caveats below** |
| `last_activity` | **yes** | **Last commit 2026-06-23; PyPI 2.1.0 2025-04-24; last tag v2.0.0 2025-01-06 (observed 2026-08-10)** |
| `documented_examples` | **yes** | **9 tutorial notebooks (Tutorial 0–8) on readthedocs + 7 `.rst` doc pages; 2 tutorials CPU-runnable with no downloads** |

### `interface`

`pyproject.toml` (master, 51 lines) declares **no** `[project.scripts]`; the file tree contains no
CLI module, no GUI, no web service and — unlike sibling tangermeme — no bundled agent `SKILL.md`. A
CLI is roadmap item 3 and agentic interfaces roadmap item 8, both future.
**Source:** `pyproject.toml` (master); GitHub tree API; README roadmap items 3 and 8.

### `license`

Apache-2.0 confirmed three ways: the `LICENSE` file is the Apache License 2.0 text at **both**
master and the v2.1.0 commit `7adfcd6453`; the README badge and closing section say Apache 2.0; the
GitHub API reports `spdx_id: Apache-2.0`. The "MIT" string on PyPI is not a stale classifier — the
v2.1.0 `classifiers` list has no license entry; it comes from `setup.py`'s `license='MIT License'`
in the distributed artifact. The shipped package therefore carries **contradictory license
metadata**; the repository `LICENSE` governs.
**Source:** `LICENSE` at master and `7adfcd6453`; `setup.py` at `7adfcd6453`; README badge; GitHub
API `spdx_id`.

### `installable_today`

**Yes** — on PyPI, `pip install ledidi` installs **v2.1.0, uploaded 2025-04-24**; that build is a
`setup.py` build with `install_requires=["torch >= 1.9.0", "matplotlib"]` (no `requires_python`, no
numpy, no tangermeme). Caveats that must travel with this cell: (i) `ledidi/plot.py` in v2.1.0
imports `pandas` and `logomaker`, **neither declared**, so `from ledidi.plot import plot_edits`
fails on a clean install, and `plot_loss` does not exist in that release; (ii) several capabilities
cited in this record exist only on **unreleased master** — input validation via
`tangermeme.utils._validate_input`, the tangermeme-based `plot_edits`, and
`greedy_pruning(threshold=0)` (`docs/whats_new.rst` "Version 2.2.0"; `ledidi/__init__.py` on master
still reads `__version__ = '2.1.0'`). Scoring the capability table from master is right — it credits
the authors with everything they have built — but the referee response must not claim the
installable wheel contains all of it. No installation was performed (per instructions), so this is
metadata-verified, not execution-verified.
**Source:** PyPI JSON API; `setup.py` and `ledidi/plot.py` at `7adfcd6453`; `docs/whats_new.rst`;
`ledidi/__init__.py`.

### `last_activity`

GitHub API on 2026-08-10: last commit `2026-06-23T21:43:47Z` (*"UPDATE plot_edits to render via
tangermeme.plot_logo (#19)"*), `pushed_at` 2026-06-23, `updated_at` 2026-07-30, not archived, 110
stars, 0 open issues, CI test workflow present. PyPI latest 2.1.0 (2025-04-24); only tag/release is
v2.0.0 (2025-01-06). Docs site live. Actively maintained; release cadence lags development. (The
"100% coverage" badge is a static shields.io badge, not a coverage service.)
**Source:** GitHub API (commits/tags/releases/repo); PyPI JSON API; `docs/whats_new.rst`.

### `documented_examples`

**Nine** tutorial notebooks in `docs/tutorials/` (Tutorial 0–8), all rendered on readthedocs, plus 7
`.rst` doc pages (`index`, `getting_started`, `input_output`, `parameters`, `faq`, `whats_new`, +
tutorial index). Tutorials 0 and 8 run on CPU with no downloads; the rest use BPNet/Enformer/Malinois
oracles from Zenodo record 14604495. Topics: 0 getting started (toy AP-1 oracle), 1 protein binding
sites, 2 constraints and priors, 3 in-painting / motif moving, 4 multiple models, 5 affinity
catalogs, 6 custom losses, 7 validating designs, 8 the Ledidi object. Beyond the docs, the paper
reports 14 design settings across 10 pre-trained models as case studies.
**Source:** `docs/tutorials/Tutorial_0..8`; readthedocs; PDF §2.4–2.6, §3.

---

## Changes from the v1 record

| Row | v1 | v2 | Why |
|---|---|---|---|
| `library_as_object` → `library_first_class_object` | partial | **partial** | Split. Ledidi is genuinely object-returning (tensor + fitted `Ledidi` module) and never writes files, but the object is unlabeled and single-template — partial survives the split. |
| `library_as_object` → `library_algebra` | partial | **partial** | Split. Repeat/sample/stack (`n_repeats`, `n_samples`, catalog list → stacked leading axis) are real in-tool operations; combining two independently specified design sets is not. |
| `dag_chaining` → `composable_operations` | no | **no** | Rename; value carried. Composition is over oracles and constraints, not design steps. |
| `lazy_evaluation` → `lazy_generation` | partial | **partial** | Rename; value carried and *strengthened*: v1 flagged `no` as equally defensible under "lazy evaluation", but under "produced on demand" the fitted-designer sampler squarely supports `partial`. |
| `mixed_mutagenesis_one_pool` → `exhaustive_single_scans` | no | **no** | Ledidi optimizes, never enumerates; deletions impossible. |
| `mixed_mutagenesis_one_pool` → `sampled_random_mutagenesis` | no | **partial** (**upgraded**) | The v1 single row was answering "can you mix types", which hid that Ledidi *is* a sampler (`n_samples=10000` from a learned weight matrix). Not `yes`: no user-specified random/rate mode. |
| `mixed_mutagenesis_one_pool` → `higher_order_combinatorial` | no | **partial** (**upgraded**) | Every design is a multi-substitution variant (~16.5 edits mean) and `greedy_pruning` searches subsets. Not `yes`: nothing enumerates or specifies combinations. |
| `mixed_mutagenesis_one_pool` → `heterogeneous_components_one_library` | no | **no** | No component types to mix; fixed length, substitutions only. |
| `hgvs_input` | no | *dropped* | Per `ROWS_v2.md`. |
| `maintained` | yes (2026-06-23) | *replaced by* `installable_today` + `last_activity` | Same underlying evidence, split into the two v2 rows. |
| `degenerate_iupac_codons` | — | **no** | New row; answered from the memo (no codon concept; `allow_N` is input tolerance only). |
| `negative_control_generation` | — | **no** | New row; dinucleotide shuffles are a paper-level evaluation done with tangermeme, not a ledidi feature. |
| `ml_model_in_loop` | — | **yes** | New row; the memo's "notable capabilities" item 1 — this is what Ledidi *is*. Previously uncaptured by any row. |
| `readout_analysis` | — | **no** | New row; validation is in-silico + retrospective lookup of published assay data. |
| `documented_examples` | — | **yes / 9 tutorials** | New row; from the memo's tutorial inventory. |

All other Block B/C/D values carried across unchanged.

## Confidence

**Highest confidence** (read directly from shipped source, docs, `pyproject.toml`, `setup.py`,
`LICENSE`, GitHub/PyPI metadata, verified twice by extractor and reviewer): all of Blocks C and E,
`primer_design`, `codon_optimization`, `synthesis_constraints`, `degenerate_iupac_codons`,
`readout_analysis`, `ml_model_in_loop`, `barcode_generation`, `automatic_naming`,
`per_sequence_provenance`, `exhaustive_single_scans`, `heterogeneous_components_one_library`,
`assay_dms`.

**Judgement calls flagged for human review:** `library_first_class_object`, `library_algebra`,
`lazy_generation`, `sampled_random_mutagenesis`, `higher_order_combinatorial`,
`composable_operations`, `combinatorial_motif_place`, `design_visualization`, `assay_mpra`,
`assay_insilico`, `negative_control_generation`, `installable_today`. See the per-row evidence for
the exact boundary each turns on.

**Two sentences most likely to draw an author objection, and how this record pre-empts them:**
(a) any claim that Ledidi cannot place motifs — it **can**, at a user-chosen position, and can
relocate them; only combinatorial *enumeration* is absent;
(b) any claim that Ledidi has no constraints — it has a real edit-composition constraint system;
what it lacks is *synthesis* constraints.

**Not verified by execution** — no installation was performed, per instructions.
