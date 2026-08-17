# Ledidi — FINAL capability record (PoolParty referee response)

**Tool:** Ledidi (`ledidi`)
**Citation key:** `Schreiber2025ledidi`
**Tier:** 3
**Status of this document:** final, post-review. Merges the extraction memo
(`extractions/ledidi.md ⟨deleted at cleanup — in commit 35d65d8⟩`) with the adversarial review (`reviews/ledidi.md ⟨deleted at cleanup — in commit 35d65d8⟩`) and my own
re-verification pass on 2026-08-10. This file is self-contained — it is the source for the
comparison table.

**One-line:** Gradient-based optimizer that turns any frozen, differentiable PyTorch
sequence-to-function model into a *sequence editor*: it learns a per-position/per-character weight
matrix (Gumbel-softmax) so that a small set of substitutions to **one** template sequence drives the
oracle's prediction to a user-specified target. It is **not** a library-design tool in the PoolParty
sense — it is a single-template, model-driven design method that happens to emit batches.

**Changes from the extraction memo (all reviewer-driven, all verified by me):**
1. `design_visualization` **yes → partial** (the single change of value in the whole record).
2. `combinatorial_motif_place`, `synthesis_constraints`, `per_sequence_provenance` — values unchanged
   (`no`), evidence rewritten so the row cannot be misread as "Ledidi cannot place motifs" / "Ledidi
   has no constraints" / "Ledidi records nothing".
3. `availability_status` — corrected: the extraction quoted **master's** `pyproject.toml` as if it
   described the installable wheel. It does not.
4. `license` — the "MIT" string is not a stale classifier; it is `setup.py`'s `license=` field in the
   shipped v2.1.0 artifact.
5. Four missed capabilities added (knock-down/knock-out symmetry, forced motif insertion +
   relocation, speed-vs-enumeration claim, documented design-validation recipe, plus minor
   engineering surface).

---

## Sources consulted

| Kind | Reference | Accessed |
|---|---|---|
| PDF (preprint) | `../../../lit_review/analyzed/Schreiber2025_ledidi.pdf` — Schreiber, Lorbeer, Heinzl, Lu, Stark, Noble, "Programmatic design and editing of cis-regulatory elements", bioRxiv 2025.04.22.650035, posted 2025-04-23, 40 pp. (PyMuPDF text) | 2026-08-10 |
| Prior analysis | `prior_analyses/Schreiber2025_ledidi_analysis.md ⟨deleted at cleanup — in commit 35d65d8⟩` | 2026-08-10 |
| Repo (master) | https://github.com/jmschrei/ledidi — `README.md`, `LICENSE`, `pyproject.toml`, `ledidi/{__init__,ledidi,losses,plot,pruning,wrappers}.py`, `docs/*` | 2026-08-10 |
| Repo (released tree) | commit `7adfcd6453` = the v2.1.0 source behind the PyPI wheel — `setup.py`, `ledidi/plot.py`, `LICENSE` | 2026-08-10 |
| Docs | repo `docs/` at master — `index/getting_started/input_output/parameters/faq/whats_new.rst` + 9 tutorial notebooks (`docs/tutorials/Tutorial_0..8`); the deployed site https://ledidi.readthedocs.io/en/latest/ is a stale v2.0.0 build that serves only `index`/`whats_new`, the 5 API pages and Tutorials 1–6 | 2026-08-10 |
| PyPI | https://pypi.org/project/ledidi/ (JSON API) | 2026-08-10 |
| GitHub API | commits, tags, releases, repo metadata, tree listing | 2026-08-10 |
| Model zoo | https://zenodo.org/records/14604495 (BPNet oracles used by paper + tutorials) | referenced, not downloaded |

Nothing was installed or executed. Prior-analysis check: the prior notes are broadly **correct**
(single-sequence ML-guided editing, not library design; affinity catalogs are repeated independent
optimizations; no design cards/metadata). Two refinements: (1) one `ledidi()` call *does* return many
sequences (`batch_size`, `n_samples`, `n_repeats`, catalog axis), so "output: one or a few sequences"
understates it — but they are a bare tensor with no identity; (2) Ledidi is Apache-2.0 and actively
developed as of June 2026, which the prior note did not cover.

## What actually ships (master)

Six modules: `ledidi/{__init__,ledidi,losses,plot,pruning,wrappers}.py`. Public surface:

- `ledidi(model, X, y_bar, n_repeats=1, n_samples=None, return_designer=False, return_history=False, device='cuda', random_state=None, **kwargs)` — one-call entry point.
- `Ledidi(model, shape, target=None, input_loss=L1Loss(sum), output_loss=MSELoss(), tau=1, l=0.1, batch_size=16, max_iter=1000, early_stopping_iter=100, lr=1.0, input_mask=None, initial_weights=None, eps=1e-4, ...)` — the `torch.nn.Module` designer (`.fit_transform`, `.forward`).
- `ledidi.wrappers.DesignWrapper([model_a, model_b])` — concatenates several oracles' predictions.
- `ledidi.losses.MinGap(in_mask)` — on-target vs off-target gap loss.
- `ledidi.pruning.greedy_pruning(model, X, X_hat, threshold=1, target=None)`.
- `ledidi.plot.{plot_loss, plot_history, plot_edits}` (`plot_edits` delegates to `tangermeme.plot.plot_logo`, `plot.py:21,193`).

I/O contract (`docs/input_output.rst`, verbatim): *"Sequences — both the template `X` you pass in and
the designs `X_hat` you get back — are **one-hot encoded** tensors of shape `(batch, n_channels,
length)` and dtype `torch.float32`"* and *"The template `X` passed to `ledidi` has a batch dimension
of **1**."* Enforced at `ledidi.py:204`: `_validate_input(X, "X", shape=(1, -1, -1), ohe=True,
allow_N=True)`. There is **no** sequence-record, FASTA, table, or writer anywhere in the API.

Multi-sequence output exists only as tensor axes: *"For an **affinity catalog**, pass a *list* of such
tensors. The design is repeated once per element and the results are stacked, adding a leading catalog
dimension to `X_hat`"*; `ledidi()` docstring: *"y: torch.Tensor, shape=(*ny, *n_repeats, n_samples,
n_channels, length)"*; `n_samples` is realized as `n_iter = n_samples // batch_size + 1` then
`torch.cat` (`ledidi.py:212–252`).

---

## Capability-by-capability record

### Block A — library specification

| Key | Value |
|---|---|
| `library_as_object` | **partial** |

**Evidence.** One call can return a rank-5 tensor `(n_catalog, n_repeats, n_samples, n_channels,
length)` — README: *"X_hat = ledidi(model, X, y_bar, n_samples=10000)  # Runs almost as fast as the
default!"* — so "returns one sequence" would be wrong. Against that: the template is hard-restricted
to a **single** sequence (`ledidi.py:204`, `_validate_input(X,"X",shape=(1,-1,-1),ohe=True,allow_N=True)`);
the return is an unlabeled `torch.Tensor` with no library abstraction, no per-element identity, and no
writer to FASTA/CSV; heterogeneous design tasks must be scripted by the user. Reviewer notes `no`
would also be defensible; `partial` is kept as the generous-but-defensible reading, and it is a
*weaker* partial than tangermeme's (tangermeme's ops are batch-first over many templates).
**Source:** `ledidi/ledidi.py:204, 212–252`; `docs/input_output.rst`; README quick-start.

| Key | Value |
|---|---|
| `dag_chaining` | **no** |

**Evidence.** The only composition in Ledidi is over *oracles*, not over design steps: `DesignWrapper`
*"concatenates their predictions"* (README; `wrappers.py` is a `torch.cat` of model outputs). No
pipeline object, graph, step composition, or nesting of design specs exists in any of the 6 modules.
README roadmap item 3 puts this in the future: *"Interoperability and end-to-end pipelines … so that
one can go from raw data all the way to finished designs using command-line tools alone."*
Pre-empting an author objection: paper §2.1 says masks, in-painting and priors *"can be seamlessly
combined with each other"* — that is **constraint composition inside one design**, not chaining of
design steps. **Source:** `ledidi/wrappers.py`; README roadmap item 3; PDF §2.1.

| Key | Value |
|---|---|
| `lazy_evaluation` | **partial** |

**Evidence.** `ledidi()` returns fully materialized tensors and there is no specification object to
defer. However `return_designer=True` yields a fitted `Ledidi` module that is a genuine on-demand
sampler: README *"once the Ledidi weight matrix is learned, edited sequences can be sampled extremely
quickly using the forward function"*; Tutorial 8 *"Fitting is the expensive part … drawing a fresh
batch of designs afterward is a single forward pass."* **Qualifier that must travel with this cell:**
successive materializations yield *different* stochastic sequences (resampling from the fitted weight
matrix), so this is on-demand generation, not deferred evaluation of a fixed specification. Reviewer
notes `no` is equally defensible; `partial` retained as the safer direction.
**Source:** `ledidi/ledidi.py:467` (`forward`); README; Tutorial 8.

| Key | Value |
|---|---|
| `mixed_mutagenesis_one_pool` | **no** |

**Evidence.** Ledidi performs **substitutions only**, over a fixed-length one-hot tensor, driven by a
single objective. There are no mutagenesis "types" to mix and no way to express "exhaustive singles +
sampled random + WT replicates" in one specification. Indels are roadmap item 5 verbatim:
*"**Insertions and deletions.** Extending Ledidi beyond substitutions so that designs can change the
length of a sequence, not just its content."* **Source:** `ledidi/ledidi.py`; README roadmap item 5.

| Key | Value |
|---|---|
| `combinatorial_motif_place` | **no** (no combinatorial enumeration — but placement *is* supported) |

**Evidence — reworded per review; the old wording was exposed.** Ledidi **can** force a specified
motif into a user-chosen position: Tutorial 2, *"Adding in Motifs"* — *"we can also use careful usage
of a mask to force edits to happen -- even editing in an entire motif. … If we set all the characters
to `-inf` except for some other character, we could force that edit ot [sic] be included in the design."*
It can also **relocate** motifs: Tutorial 3, *"Moving Motifs (CTCF + Transcription)"* — in-paint the
old site, in-paint a landing pad, and hold the oracle output fixed. And placements **can** be swept,
but only in a user-written loop: Tutorial 3, *"let's take a comprehensive look at what the model
predictions are if we try moving the 20bp region of the BPNet motif downstream between 1 and 800 bp.
Every 10th basepair …"*. What does **not** exist is any operation that enumerates motif sets ×
positions × orientations × spacings × permutations as a library specification. The cell is `no`
strictly on the absence of combinatorial enumeration.
**Source:** Tutorial 2 ("Adding in Motifs"), Tutorial 3 ("Moving Motifs", positional sweep);
`ledidi/ledidi.py` (`input_mask`, `initial_weights`).

| Key | Value |
|---|---|
| `barcode_generation` | **no** |

**Evidence.** Zero occurrences of `barcode`, `oligo`, or `adapter` across all 6 modules, all 11 `.rst`
pages, the README, and the tutorial notebooks (grepped independently by extractor and reviewer). No
edit-distance or GC-constrained tag generation, no concatenation of tags onto designs — and output
length is always the template length, so tags could not be appended in the first place.
**Source:** repo-wide grep; `docs/input_output.rst` (output shape == input shape).

| Key | Value |
|---|---|
| `per_sequence_provenance` | **no** |

**Evidence — phrasing tightened per review.** The output is an unlabeled tensor; **no per-sequence
record travels with it**. README tells users to diff by hand: *"To see exactly which positions changed
and how, compare the original sequence to a designed one"* followed by a manual `torch.where` loop.
Nuance to state precisely so the row cannot be attacked: `return_history=True` does record
`history['edits'].append(torch.where(X_hat != X_))` (`ledidi.py:603`), i.e. *(batch index, channel,
position)* triples **per optimization iteration** — so a batch-member index does exist. But the
returned `best_sequence` is a clone from whichever iteration last improved the loss and **that
iteration index is never recorded**, so nothing maps a returned design back to its history rows.
README roadmap lists *"**reproducibility and provenance** — capturing the full recipe behind every
design so it can be audited and regenerated"* as a cross-cutting **future** thread.
**Source:** `ledidi/ledidi.py:603–606`; README (diff recipe + roadmap provenance thread).

| Key | Value |
|---|---|
| `automatic_naming` | **no** |

**Evidence.** Designs are slices of a `torch.float32` tensor; nothing in the API assigns names, IDs,
or labels. The docs decode designs to strings with a hand-written local `decode()` helper.
**Source:** `docs/getting_started.rst`; `docs/input_output.rst`; all 6 modules.

| Key | Value |
|---|---|
| `design_visualization` | **partial** — *corrected from `yes`* |

**Evidence.** Ledidi ships `plot_loss` (input/output loss curves), `plot_history` (*"the position of
each proposed edit over the course of the optimization"*) and `plot_edits` (*"attribution values as a
series of stacked logo tracks, with characters colored if they are edited with respect to the initial
sequence"*). Three reasons this is `partial`, not `yes`:
1. **Shipped surface.** The module calls these *"thin matplotlib helpers"*: `plot_edits` renders
   user-supplied attribution tensors and explicitly leaves attribution computation and broader
   sequence-level plotting to the caller (`plot.py:12–14, 127–129`). It is a plotting helper, not a
   full design-inspection system.
2. **Scope.** It renders edits against **one reference sequence at a time** and requires the user to
   supply attribution tensors — the module docstring says attribution computation and broader
   sequence-level plotting are outside these helpers' scope. There is no library-level view, no
   design-set summary, and no design-graph view (there is no design graph).
3. **Release reality.** The `yes` rested on unreleased master. In the pip-installable v2.1.0 (commit
   `7adfcd6453`) `ledidi/plot.py` imports `pandas` and `logomaker`, **neither declared** in that
   build's `install_requires` (`torch >= 1.9.0`, `matplotlib`), and `plot_loss` **does not exist**
   there at all (verified: v2.1.0 `plot.py` defines only `plot_history` and `plot_edits`).
**Source:** `ledidi/plot.py:12–21, 24, 85, 118, 127–129, 193` (master); `ledidi/plot.py:4–8, 11, 44`
and `setup.py` at commit `7adfcd6453`.

### Block B — assay coverage

| Key | Value |
|---|---|
| `assay_dms` | **no** |

**Evidence.** No codon-aware or amino-acid-level mutagenesis, no reading-frame logic, no exhaustive
variant enumeration, no coding-sequence handling, and no protein example documented. The only
protein-adjacent statement is a modality claim about the *method*: README, *"one can also apply Ledidi
out-of-the-box to RNA or protein models (or really, to any model with a sequence of categorical
inputs such as small molecules)"*, plus one *"amino acids"* mention in a prose docstring. Consistent
with this survey's own bar: tangermeme scores `no` here despite shipping `saturation_mutagenesis`.
**Source:** repo-wide grep (codon/amino-acid/reading-frame/exon); README.

| Key | Value |
|---|---|
| `assay_mpra` | **partial** |

**Evidence.** Ledidi designs regulatory sequences **against MPRA/STARR-seq oracles**: paper §3.1 lists
*"TF binding, chromatin accessibility, transcription, or massively parallel reporter assay (MPRA)
activity"* among design targets; Malinois (MPRA, K562) and DeepSTARR are among the 10 pre-trained
oracles across 14 design settings (Fig 2A, Suppl. Table 1); §3.4 validates designed 8-mer affinity
catalogs against Reiter et al. STARR-seq data in *Drosophila* enhancers (r = 0.893–0.925, p.25). But
**no MPRA construct assembly ships**: grep found no barcode, adapter, primer, oligo, or synthesis
handling; no oligo-length or replicate/control structure; no orderable output file. Consistent with
the tangermeme row.
**Source:** PDF §3.1, §3.4 (p.25), Fig 2A / Suppl. Table 1; repo-wide grep.

| Key | Value |
|---|---|
| `assay_insilico` | **partial** |

**Evidence.** Every Ledidi design is an in-silico interrogation of a genomic AI model, and affinity
catalogs are explicitly framed as probes of the learned code (paper Introduction, p. 4: *"an affinity
catalog for GATA2 binding reveals a sophisticated usage of a learned cis-regulatory code"*; same
paragraph: catalogs *"provide a resource for understanding the full range of potential sequence edits
and their relative strengths"*). Paper §2.6 also uses in-silico marginalization (implanting a motif
into 1,000 dinucleotide-shuffled backgrounds) as an evaluation. But Ledidi ships **no declarative perturbation
operations** — the survey's bar for `yes` is tangermeme's `marginalize`, `saturation_mutagenesis`,
`variant_effect.*`, all of which live in tangermeme, not ledidi. Ledidi reaches in-silico
interrogation only through repeated optimization; Tutorial 7 (attributions → FIMO hit deltas →
Enformer / ATAC-BPNet round-trip) is *design validation*, not a perturbation framework.
**Source:** PDF §2.6, §3.3, Fig 4E, Introduction p. 4; Tutorial 7; module listing (no perturbation API).

### Block C — genomics integration

| Key | Value |
|---|---|
| `genome_coordinates` | **no** |

**Evidence.** The API accepts only one-hot `(1, n_channels, length)` tensors; coordinate and FASTA
handling is explicitly delegated: `docs/input_output.rst`, *"The sibling library tangermeme provides
`one_hot_encode` and `characters` utilities that handle this — as well as reading FASTA and genomic
loci — for you."* Worth naming so the row is not misread: the tutorials **do** work from hg38 loci
(e.g. the SMYD3 promoter), but via `pyfaidx` + `tangermeme.utils.one_hot_encode` **outside** ledidi's
API — which is precisely the point.
**Source:** `docs/input_output.rst`; `ledidi.py:204`; Tutorials 1–3.

| Key | Value |
|---|---|
| `transcript_models` | **no** |

**Evidence.** Zero GTF/GFF hits anywhere in repo, docs, or notebooks. Gene-body protection is achieved
only by a **user-supplied boolean position mask**: Tutorial 2 *"Preventing Edits (Gene Bodies)"*
(*"one may wish to prevent edits at the TATAA box or within the coding region"*, implemented as a
hand-built boolean mask); paper §2.5 *"we masked out the first 1,057bp of sequence"*.
**Source:** repo-wide grep; Tutorial 2; PDF §2.5.

| Key | Value |
|---|---|
| `exon_intron_split_codons` | **no** |

**Evidence.** No exon, intron, reading-frame, or codon concept anywhere in the 6 modules or the docs;
the sequence is an unannotated fixed-length categorical array.
**Source:** repo-wide grep over modules + `.rst` + notebooks.

| Key | Value |
|---|---|
| `hgvs_input` | **no** |

**Evidence.** No HGVS or variant-nomenclature handling. The only accepted input is a one-hot
`torch.float32` tensor validated by `tangermeme.utils._validate_input`; the FAQ's accepted-input list
contains no variant syntax.
**Source:** `ledidi/ledidi.py:26, 204`; `docs/faq.rst`; `docs/input_output.rst`.

| Key | Value |
|---|---|
| `vcf_vep_output` | **no** |

**Evidence.** Reading the full return path of `ledidi()`: outputs are `X_bar` (tensor), optional
designer objects, optional history dicts, plus matplotlib axes from `ledidi.plot`. No VCF, BED, GFF,
FASTA, or CSV writer exists in any module.
**Source:** `ledidi/ledidi.py` (return path); `ledidi/plot.py`; module listing.

| Key | Value |
|---|---|
| `consequence_annotation` | **no** |

**Evidence.** No molecular-consequence classification anywhere. The readouts in the repo and the
paper — oracle predictions, edit counts/positions, attributions, motif-hit counts, correlations with
experimental data, runtimes — never include a variant-consequence class; "consequence" in the paper
means predicted regulatory activity, not variant-effect classification.
**Source:** repo-wide grep; PDF results sections.

### Block D — physical construction

| Key | Value |
|---|---|
| `primer_design` | **no** |

**Evidence.** Grep finds no primer, oligo, or wet-lab protocol output in repo, docs, notebooks, or
paper. The paper frames the wet-lab route as CRISPR base/prime editing but designs **no reagents**;
its "experimental validation" is a retrospective lookup against published STARR-seq and DREAM/GPRA
data, with no wet-lab work of its own.
**Source:** repo-wide grep; PDF §3.4 (Reiter STARR-seq, DREAM/GPRA via LegNet).

| Key | Value |
|---|---|
| `codon_optimization` | **no** |

**Evidence.** The word "codon" appears exactly once in the entire repo, in README roadmap item 4 as
**future** work: *"Designing for several properties at once while respecting hard constraints such as
GC content, restriction sites, codon usage, and the absence of off-target effects."*
**Source:** README roadmap item 4.

| Key | Value |
|---|---|
| `synthesis_constraints` | **no** (but a real *edit-composition* constraint system exists — name it) |

**Evidence — reworded per review.** There is **no** GC-content target, homopolymer/repeat check,
restriction-site avoidance, Tm calculation, or post-hoc synthesizability check anywhere; all are named
in roadmap item 4 as future work. What **does** exist, and what an author would cite, is a constraint
system over *edits*: Tutorial 2, *"Preventing Edits (No Added Cs)"* — `initial_weights[:, 1, :] =
float("-inf")` forbids ever introducing a C, motivated explicitly by editing feasibility (*"you might
want to prevent certain nucleotides from being edited in, potentially due to issues with prime or base
editing"*), down to blocking a single specific substitution (*"the G in each of the GATAA motif might
be blocked from becoming a T"*). That is an **edit-composition constraint, not a synthesis
constraint** — it constrains which edits the optimizer may propose, never whether the resulting
molecule is orderable.
**Source:** Tutorial 2 ("Preventing Edits (No Added Cs)"); README roadmap item 4; `ledidi.py:430`.

### Block E — engineering

| Key | Value (enum) | Descriptor for the table |
|---|---|---|
| `interface` | **yes** | **Python API (library) only — no CLI, GUI, or web service** |

**Evidence.** `pyproject.toml` (master, 51 lines) declares **no** `[project.scripts]`; the file tree
contains no CLI module, no GUI, no web service, and — unlike sibling tangermeme — no bundled agent
`SKILL.md`. A CLI is roadmap item 3 and agentic interfaces roadmap item 8, both future. Note on
encoding: the extraction recorded the bare enum `yes`, which is not an answer to "interface"; the
descriptor above is what belongs in the comparison table (other tools in this survey carry
descriptors such as "CLI only").
**Source:** `pyproject.toml` (master); GitHub tree API; README roadmap items 3 and 8.

| Key | Value (enum) | Descriptor for the table |
|---|---|---|
| `license` | **yes** | **Apache-2.0** (shipped v2.1.0 `setup.py` says `license='MIT License'` — repo `LICENSE` is authoritative) |

**Evidence.** Apache-2.0 confirmed three ways: the `LICENSE` file is the Apache License 2.0 text at
**both** `master` and the v2.1.0 commit `7adfcd6453`; the README badge and closing section say Apache
2.0; the GitHub API reports `spdx_id: Apache-2.0`. **Correction to the extraction:** the "MIT" string
on PyPI is *not* a stale classifier — the v2.1.0 `classifiers` list contains no license entry at all.
It comes from `setup.py`'s `license='MIT License'` field in the distributed artifact (verified
directly: `setup.py` at `7adfcd6453` reads `license='MIT License'`). So the shipped package carries
**contradictory license metadata**; the repository `LICENSE` governs.
**Source:** `LICENSE` at master and `7adfcd6453`; `setup.py` at `7adfcd6453`; README badge; GitHub API
`spdx_id`.

| Key | Value (enum) | Descriptor for the table |
|---|---|---|
| `maintained` | **yes** | **yes — last commit 2026-06-23; PyPI 2.1.0 (2025-04-24), master well ahead** |

**Evidence.** GitHub API on 2026-08-10: last commit `2026-06-23T21:43:47Z` (*"UPDATE plot_edits to
render via tangermeme.plot_logo (#19)"*), `pushed_at` 2026-06-23, `updated_at` 2026-07-30, not
archived, 110 stars, 0 open issues, CI test workflow present. PyPI latest **2.1.0 (2025-04-24)**; the
only tag/release is v2.0.0 (2025-01-06); `docs/whats_new.rst` carries an unreleased "Version 2.2.0"
section while `ledidi/__init__.py` on master still reads `__version__ = '2.1.0'`. Caveat worth
carrying: the "100% coverage" badge is a static shields.io badge, not a coverage service.
**Source:** GitHub API (commits/tags/releases/repo); PyPI JSON API; `docs/whats_new.rst`;
`ledidi/__init__.py`.

---

## Documented examples / vignettes (candidates for PoolParty reproduction)

Nine tutorial notebooks in `docs/tutorials/` at master; the deployed readthedocs build is stale
(v2.0.0) and renders only Tutorials 1–6. Tutorials 0 and 8 run on CPU with no downloads; the
rest use pre-trained oracles — the BPNet models from Zenodo record 14604495, Enformer from
EleutherAI via `enformer-pytorch`, and Malinois/ChromBPNet/Beluga checkpoints from elsewhere.

0. **Tutorial 0 — Getting Started / A First Design.** Toy parameter-free AP-1 (`TGACTCA`) convolution
   oracle on a random 50 bp sequence; Ledidi inserts the motif in ~2 edits. Fully reproducible with no
   data or GPU.
1. **Tutorial 1 — Design of Protein Binding Sites.** BPNet GATA2 (K562) oracle; design edits that
   create GATA2 binding; inspect the proposed edits.
2. **Tutorial 2 — Constraints and Priors.** Hard `input_mask` (protect gene body / TATAA box), soft
   `initial_weights` priors, blocking a whole character ("No Added Cs") or one specific substitution,
   and **forcing an entire motif in** at a chosen position ("Adding in Motifs"). Also shows using
   masks to choose **which** of three GATA sites gets knocked out.
3. **Tutorial 3 — In-Painting.** Zero out a span so Ledidi fills it with "free" edits; add a motif to
   background, **replace** an SP-1 site while holding accessibility constant, and **move** a CTCF site
   (in-paint old site + in-paint landing pad + hold oracle output fixed), plus a user-written sweep of
   the motif's placement over 1–800 bp.
4. **Tutorial 4 — Multiple Models.** `DesignWrapper` over several oracles/tasks: cell-type-specific,
   cross-modal, and "design-property" objectives; also user-written wrappers for oracles with
   different input widths, with the warning that *"Ledidi will try to hide the edits from the models
   with the smaller receptive field"* unless the regions not seen by every model are masked.
5. **Tutorial 5 — Affinity Catalogs.** Pass a *list* of `y_bar` targets to sweep design strength;
   GATA2 BPNet example.
6. **Tutorial 6 — Custom Loss Functions.** Custom output losses, including using BPNet base-resolution
   profiles to control *where* a site is placed.
7. **Tutorial 7 — Validating Your Designs.** Attributions, FIMO motif-hit counting against decoys, and
   round-tripping designs through held-out models (Enformer, ATAC BPNet).
8. **Tutorial 8 — The Ledidi Object.** Fit the weight matrix once, then sample many designs cheaply
   (`return_designer=True`, `n_samples`). CPU-runnable.

**Paper-level case studies** (Methods §2.4–2.6, Results §3): 14 design settings across 10 pre-trained
models (BPNet GATA2/E2F6, Beluga CTCF/MAX, BPNet-ATAC, ChromBPNet, Enformer DNase/CAGE, Sei, ProCapNet,
Puffin, Borzoi, Malinois K562, DeepSTARR); percentile-based design from the 50th to the 99.27th
percentile on chr12/hg38 at ~16.54 edits mean; affinity catalogs at the SMYD3 promoter (GATA2 binding
+0.2…+10.0 in steps of 0.2 with the first 1,057 bp masked; accessibility and transcription catalogs
spanning −3.0 to +4.0); 8-mer affinity catalogs at seven sites in the *Drosophila* ced6 and ZnT63C
enhancers validated against STARR-seq (r = 0.893–0.925); 8-mer catalogs in yeast 80 bp promoters
validated against DREAM/GPRA data via LegNet; runtime benchmark vs greedy substitution and two motif
implantation baselines (Fig 2G).

## Notable capabilities not covered by the row list

1. **Model-in-the-loop, objective-driven design** — the core of the tool. No row captures "generate
   sequences by optimizing against a differentiable predictive oracle to hit a numeric target."
2. **Knock-down / knock-out design, fully symmetric with knock-in** *(added per review)*. Ledidi is
   not a "make it stronger" tool: Tutorial 2 diminishes GATA2 binding and uses masks to choose *which*
   of three binding sites is destroyed (*"we have used the mask to control which binding sites get
   removed by the editing process"*); Tutorial 3 removes an SP-1 site while preserving accessibility;
   the paper's SMYD3 accessibility/transcription catalogs span −3.0 to +4.0. Any table saying Ledidi
   "designs sequences with higher predicted activity" understates it.
3. **Forced motif insertion at a chosen position, and motif relocation** *(added per review)*. Setting
   all characters but one to `-inf` in `initial_weights` forces a whole motif in (Tutorial 2, "Adding
   in Motifs"); in-painting the old site plus a landing pad, with the oracle output held fixed, moves a
   site (Tutorial 3, "Moving Motifs"). This is the nearest analogue Ledidi has to PoolParty's motif
   operations — and it is *placement without enumeration*.
4. **Multi-objective design across several models** — `DesignWrapper`, plus `MinGap` loss for
   on-target vs off-target (e.g. cell-type-specific enhancers). `MinGap` ignores `y_bar` entirely and
   maximizes the gap between the *weakest* on-target and *strongest* off-target output, so it controls
   relative separation rather than absolute activity.
5. **Edit-count minimization + greedy pruning** — an explicit input loss penalizing the number of
   edits (itself replaceable: `input_loss` takes any differentiable callable, so which changes count
   as costly is user-definable), plus post-hoc `greedy_pruning` to revert edits that change the
   prediction by less than a threshold (one designed sequence per call, scored by the summed absolute
   change in the oracle output relative to the fully edited design). Directly aimed at genome-editing
   feasibility; absent from all the library-design tools in this survey.
6. **Hard masks, soft priors, and edit-type blocking over edit space** — `input_mask`,
   `initial_weights` with `-inf`/`+inf`/graded logits, including blocking a *substitution type* for
   base/prime-editor feasibility. Only `-inf` is hard; finite priors *"nudge, but do not force, the
   design"* (`docs/input_output.rst`) and can be optimized away.
7. **In-painting** — zeroed spans treated as free edits (input loss set to zero there), the conceptual
   opposite of masking. It is preferential, not exclusive: positions outside the span stay editable
   unless masked, and any `N`/all-zero column is a candidate for editing and is not preserved.
8. **Affinity catalog** — a swept series of design targets producing a graded set of sequences; the
   closest thing here to a "library", but produced by independent optimizations with no shared
   identity.
9. **Speed as an explicit competitive claim against enumeration** *(added per review)*. Paper §3, Fig
   2G: Ledidi is *"almost an order of magnitude faster than greedy substitution and several orders of
   magnitude faster than a comprehensive motif implantation method"*, where the baselines implant
   consensus sequences *"at each possible position"* from JASPAR2024 (or *"a set of 17 motifs thought
   to be relevant"*). If the referee response contrasts enumerative library design with optimization,
   the authors will cite this figure — it is an argument about *cost*, not about *coverage or
   auditability*.
10. **Documented design-validation example, not a shipped Ledidi capability** — Tutorial 7 documents
    attributions, FIMO motif-hit deltas against decoy motifs, and round-tripping through independently
    trained oracles (Enformer, ATAC BPNet). Those checks are notebook code built from external
    packages and models; Ledidi exposes no validation API, and README roadmap item 2 says first-class
    evaluation remains future work.
11. **Reproducible sampling without mutating global RNG** — `random_state` backed by a private
    `torch.Generator` (master only; the released v2.1.0 signatures do not accept `random_state`
    despite a stale docstring, and sample from the global RNG).
12. **GPU-native and modality-agnostic** — DNA/RNA/protein or any categorical alphabet, any
    differentiable PyTorch oracle.
13. **Minor engineering surface** *(added per review)* — `target=` slices one task out of a multi-task
    oracle; `verbose` gives per-iteration edit/loss logging; the test suite is CPU-runnable with GPU
    tests gated behind a pytest marker (`addopts = "-m 'not gpu'"`, `pyproject.toml`).

## Stated limitations (author's own words, plus verified release caveats)

- **Substitutions only; no insertions or deletions** yet (roadmap 5) — designs cannot change sequence
  length.
- **No global sequence-property or synthesis constraints yet** — GC content, restriction sites, codon
  usage, and off-target effects are all roadmap 4. (Constraints over *edits* do exist; see
  `synthesis_constraints`.)
- **No CLI and no end-to-end pipeline** (roadmap 3); **no provenance capture** (cross-cutting roadmap
  thread); agentic interfaces are roadmap 8.
- **Designs within one returned batch are correlated.** FAQ: *"Every sequence in the returned batch is
  sampled from the same learned weight matrix … If a motif could have been inserted at several
  locations, Ledidi commits to one and all sampled sequences reflect that choice."* Independent
  designs require `n_repeats` or a new `random_state`.
- **Oracle-exploitation risk.** Tutorial 7: *"gradient-based design can produce sequences that score
  well on the oracle yet are not biologically meaningful"*; paper §3.4 notes affinity catalogs are
  *"most accurate within the dynamic range of the training data"*.
- **Requires a differentiable PyTorch oracle** — not necessarily pre-trained; the README quickstart
  uses a parameter-free toy convolution — and practically a GPU (`device='cuda'` default).
- Self-described as **"research software under active development"** (README roadmap).
- **Release lag (my addition, not the authors').** Several capabilities cited above exist only on
  unreleased master: input validation, `plot_loss`, the tangermeme-based `plot_edits`, and
  `random_state`'s private `torch.Generator`. See availability below.

## Availability today (2026-08-10) — *corrected*

**Installable and runnable: yes**, `pip install ledidi` → **v2.1.0, uploaded 2025-04-24**.

**Correction to the extraction memo:** it described the wheel as *"requires Python ≥3.10, torch ≥2.0,
tangermeme ≥1.3.0, numpy, matplotlib"*. Those are **master's `pyproject.toml`** requirements, not the
released package. The published 2.1.0 is a `setup.py` build (verified at commit `7adfcd6453`):

```
license='MIT License',
install_requires=["torch >= 1.9.0", "matplotlib"]
```

— i.e. **no `requires_python`, no numpy, no tangermeme**. Consequences, all verified:

- `ledidi/plot.py` in v2.1.0 imports `pandas` and `logomaker`, neither declared, so
  `from ledidi.plot import plot_edits` fails on a clean install; `plot_loss` does not exist in that
  release at all, and that release's `plot_edits` hard-codes the DNA alphabet `A, C, G, T`.
- Input validation via `tangermeme.utils._validate_input`, the tangermeme-based `plot_edits`, and
  `random_state`'s private `torch.Generator` are **unreleased 2.2.0 work on master**
  (`docs/whats_new.rst`); `ledidi/__init__.py` on master still reads `__version__ = '2.1.0'`.
- Most evidence quotes in this record (`_validate_input`, FAQ, `input_output.rst`, `parameters.rst`,
  the roadmap, Tutorials 0–8) come from **master**. That is the right basis for a capability table —
  it credits the tool with everything its authors have built — but the referee response must not also
  claim the installable wheel contains them.

**Repo/service health:** GitHub repo alive, not archived, last commit 2026-06-23, CI configured, 0
open issues, 110 stars, Apache-2.0. Docs site https://ledidi.readthedocs.io/en/latest/ is live but
stale (a v2.0.0 build; the getting-started, input/output, parameters and FAQ pages and Tutorials 0, 7
and 8 return 404 there). The BPNet oracles used by the tutorials are on Zenodo (record 14604495);
their Enformer, Malinois, ChromBPNet and Beluga models come from elsewhere. **No web service, no CLI,
no GUI** — it is a Python library. The paper remains a bioRxiv preprint (2025.04.22.650035, posted 2025-04-23,
CC-BY 4.0); it has no formal "Code availability" section, but README/docs point to the repo.

## Unresolved disagreements and confidence

**No unresolved extractor/reviewer disagreements remain.** The reviewer proposed exactly one value
change (`design_visualization` yes → partial), which I independently verified and applied; every other
verdict was "supported". No cell is set to `unknown`.

**Two cells where extractor and reviewer agreed on the value but explicitly flagged that the opposite
call is also defensible** — recorded here so the referee response does not over-lean on them:

- `library_as_object = partial`. `no` is defensible (one template only; unlabeled tensor; no writer).
  `partial` was kept as the *generous* reading; it is a strictly weaker partial than tangermeme's.
- `lazy_evaluation = partial`. `no` is equally defensible (each materialization gives *different*
  stochastic sequences, so this is resampling, not deferred evaluation of a fixed spec). `partial` was
  kept because it costs nothing and removes an attack surface.

**Highest confidence:** everything read directly from the shipped source, docs, `pyproject.toml`,
`setup.py`, LICENSE, and GitHub/PyPI metadata — all of Blocks C, D, E, and every `no` in Blocks A and
B.

**Lower confidence / judgement calls:** `library_as_object`, `lazy_evaluation` (see above);
`assay_mpra` and `assay_insilico` (`partial` — MPRA/STARR oracles are central to the paper and designs
are used to probe models, but there is no MPRA construct assembly and no declarative perturbation
framework); `design_visualization` (`partial` — the plotting exists and is real, but it uses one
reference sequence at a time, needs user-computed attributions, and is partly unreleased).

**Two sentences most likely to draw an author objection, and how this record pre-empts them:**
(a) any claim that Ledidi cannot place motifs — it **can**, at a user-chosen position, and can
relocate them; only combinatorial *enumeration* is absent (see `combinatorial_motif_place`);
(b) any claim that Ledidi has no constraints — it has a real edit-composition constraint system; what
it lacks is *synthesis* constraints (see `synthesis_constraints`).

**Not verified by execution** — no installation was performed, per instructions.
