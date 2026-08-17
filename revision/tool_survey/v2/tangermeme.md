# tangermeme — v2 capability record (revised row set)

**Tool:** tangermeme · **Slug:** `tangermeme` · **Citation key:** `Schreiber2025nd` · **Tier:** 1
**Author:** Jacob Schreiber (IMP Vienna / UMass Chan)
**Version at survey:** v1.4.0 (PyPI, 2026-06-25); repo `pushed_at` 2026-07-19
**Basis:** the final v1 record (`final/tangermeme.md`, extraction → adversarial review → final),
re-scored against the v2 row set. Targeted re-checks of the live `main` tree were run this pass
for the seven new rows only (see "Re-checks performed"). Nothing was installed.

## Scope rule applied

tangermeme is a general-purpose "everything-but-the-model" toolkit, not a library-design tool.
Only functionality that participates in **specifying a DNA sequence library** is scored:
`ersatz.*`, `marginalize*`, `ablate*`, `space`, `saturation_mutagenesis`, `variant_effect.*`,
`product.*`, `design/*`, `match`, `io`. Explicitly **excluded and not allowed to inflate any
cell**: DeepLIFT/SHAP, PISA, seqlet calling, seqlet annotation / tomtom-lite, logo plotting as
interpretation, k-mer counting, encoding/tokenization speed, general file I/O. Being general
purpose is not a criticism; it is a scoping fact.

Data model, which drives most of Block A: every sequence collection is a bare
`torch.Tensor` of shape `(-1, len(alphabet), length)`. A whole-package `grep "^class "` returns
**six classes** — five `NamedTuple`s of model outputs (`SpaceResult`,
`SaturationMutagenesisRawResult`, `PerturbationResult`, `PerturbationAnnotationsResult`,
`AttributionReferencesResult`) and one `UserWarning` subclass. **No class represents a sequence,
a set of sequences, a library, or a design.**

## Re-checks performed this pass (new rows only)

Full `main` tarball downloaded to a scratch dir and grepped (55 `.py` incl. tests, 17 `.md`,
`docs/` tree, `pyproject.toml`):

- `iupac` → **zero hits**. `degener` → 2 hits, both in `tests/test_saturation_mutagenesis.py`
  meaning "degenerate input" (a zero-width / single-position span), not degenerate bases.
  `codon` → **zero hits**. `utils.pwm_consensus` returns the **argmax** consensus, one-hot, no
  ambiguity code. The only ambiguity concept in the package is `N` / `[0,0,0,0]` handling in
  encoding (`whats_new.rst` L513) and the `N` padding in `greedy_marginalize`'s construct.
- `fastq`, `demultiplex`, `count table`, `barcode`, `primer` → **zero hits**. (The one `readout`
  hit is the README prose "a model that predicts a readout of interest".)
- `docs/` enumerated: 14 tutorial notebooks, 3 vignettes, 1 how-to notebook, 2 paper-figure
  notebooks, 18 API reference pages; `tangermeme/_skills/data/references/` holds 14 reference
  `.md` docs beside `SKILL.md`.
- `pyproject.toml` re-read: `license = "MIT"`, `requires-python = ">=3.10"`, 12 runtime deps
  (numpy, scipy, pandas, torch>=2.0, pyfaidx, pybigtools, tqdm, numba, joblib, scikit-learn,
  matplotlib, memelite); extras `dev` / `docs` / `interactive` (mpld3) only; `[project.scripts]`
  has exactly one entry, `tangermeme-install-skills`.
- `ersatz.randomize` signature re-read: `probs` is shape `(1, len(alphabet))` **or**
  `(len(X), len(alphabet))` — per-example, **not per-position**; there is no per-site mutation
  rate. `design/screen.py` re-read: generates random sequences via `utils.random_one_hot` (or a
  user-supplied `func=`) and keeps `n_best`.

---

## BLOCK A — library specification

### `library_first_class_object` → **partial**
*(split from v1 `library_as_object` = partial; value carried)*

tangermeme is emphatically **not** a file-writing tool — it is an in-memory API. The sequence
sets it produces (`ersatz.substitute` / `insert` / `delete` / `multisubstitute` / `randomize` /
`shuffle` / `dinucleotide_shuffle`, `beam_substitution`, `screen`) are real `torch.Tensor`
objects the user holds, slices, inspects, transforms and passes onward, and operations are
batch-first so no user loop is needed to apply one operation to many sequences
(`ersatz.py` `substitute` docstring: *"If a motif with a batch size equal to that of X is
provided, there will be 1-1 correspondence between the motifs and the sequence"*). Only one
writer exists in the whole package (`io.one_hot_to_fasta`), and it is optional.

Why not "yes": the object is a bare tensor, not a library. The package defines no class holding
sequences (six classes total, five NamedTuples of model outputs plus a Warning), so the "library"
carries no identity, no names, no membership, no metadata, and cannot be inspected *as a design*
— only as an array. Nothing distinguishes member from member except row index.

*Source:* full-tree `grep "^class "`; `ersatz.py` L20–342 signatures and `substitute` docstring;
`io.py` L540–590 (sole writer); `design/beam_substitution.py` L305–306; `design/screen.py` L142.

### `composable_operations` → **partial**
*(rename of v1 `dag_chaining` = partial; value carried)*

A genuine, deliberate composition mechanism exists: the `func=` plug-point, contract
`func(model, X, args=None, **kwargs) -> torch.Tensor | list[torch.Tensor]`
(`_skills/data/references/func-pattern.md` L11), accepted by `ablate`, `marginalize`, `space`,
`variant_effect.*`, `product.apply_pairwise`, `product.apply_product`. Nesting is documented and
demonstrated (Tutorial B7):

```python
apply_pairwise(marginalize, model, X, args=(cell_states,),
               additional_func_kwargs={'func': deep_lift_shap}, motif="TGACGTCA")
```

Paper p.3 (verbatim): *"An important consequence of tangermeme's design is that operations can be
stacked."* Fig 1A calls it "Stacked Operations". There is no fixed pipeline anywhere; ersatz
operations are pure tensor→tensor functions and chain freely in plain Python.

Why not "yes": what the **tool-level** mechanism composes is *sequence-manipulation → model
operation* — every leaf is a model call returning predictions or attributions, never a set of
sequences. Composition of *design steps into a library specification* exists only as raw
primitives the user chains by hand, which caps this at "partial" under the survey's
raw-primitives rule. **A referee could reasonably argue "yes" under the new wording** (the v2
phrasing dropped the graph-object requirement that carried part of the v1 "partial"); flagged.

*Source:* `_skills/data/references/func-pattern.md`, `product.md`; `product.py`; paper p.3;
Tutorial B7; `ersatz.py`.

### `lazy_generation` → **partial**
*(rename of v1 `lazy_evaluation` = partial; value carried)*

`product.md` L50, verbatim: *"Batches are built iteratively, so the full product is never
materialized in memory"* — confirmed in `product.py`
(`for x in tqdm(itertools.product(X, *args))` accumulating to `batch_size`). `io.extract_loci`
documents `max_jitter` expansion to *"reduce the memory footprint"*.

Counter-evidence: `ersatz.*` returns fully materialized tensors;
`saturation_mutagenesis.py` L212–235 eagerly builds `X_ = X[i].repeat(n_edits, 1, 1)` per example;
`design/greedy_substitution.py` L198 materializes `X_ = X.float().repeat(input_idxs.shape[0], 1, 1)`
every motif, every round. **Streaming exists on the Cartesian-product axis only, never as a
library abstraction.**

*Source:* `product.md` L50; `product.py`; `saturation_mutagenesis.py` L212–235;
`design/greedy_substitution.py` L188–205; `io.py`.

### `library_algebra` → **no**
*(split from v1 `library_as_object`)*

No operation in tangermeme combines, samples, repeats or otherwise treats whole sequence sets as
operands. `grep "^def .*(concat|stack|combine|sample|pool)"` over the package returns nothing
relevant. `apply_product` / `apply_pairwise` cross `X` with **extra model inputs** (the worked
example is DragoNNFruit cell-state × read-depth), not with other libraries. The n-axis on
`randomize` / `shuffle` / `dinucleotide_shuffle` replicates *within one generation call*; it is
not a repeat operator over a library.

Building a library with several design intents requires the user to `torch.cat` tensors
themselves (a `torch` call, not a tangermeme operation) and to track membership outside the tool.
By the row's explicit rule — "if the user must concatenate outputs with an external script, this
is no" — this is **no**. Noted honestly: because the data model is a plain tensor, that external
concatenation is one line rather than a script, so a lenient reader might say "partial"; flagged.

*Source:* full-tree function grep, this pass; `product.py` L30–60 docstring;
`ersatz.py` L343–360, L425–460, L608–672; final memo §3 `library_as_object`.

### `exhaustive_single_scans` → **partial**
*(new sub-row from v1 `mixed_mutagenesis_one_pool` = no; raised)*

The exhaustive scan genuinely exists and is a headline feature.
`saturation_mutagenesis` enumerates **every single substitution at every position**
(docstring L85: *"…of the sequences with an edit distance of one on them"*;
L212–235 builds `X_ = X[i].repeat(n_edits, 1, 1)` covering the full single-mutant set).
Separately, `design/_substitute.py` — whose entire docstring is *"This function takes a motif and
inserts it at all possibilities"* — tiles a motif at **every allowed position** each round of
`greedy_substitution` / `beam_substitution`, both strands when `reverse_complement=True`.

Why not "yes", precisely:
1. **Substitutions only.** There is no exhaustive deletion (or insertion) scan;
   `variant_effect.deletion_effect` / `insertion_effect` apply a *user-supplied* table of events,
   they do not enumerate.
2. **The enumerated sequences are never returned.** `saturation_mutagenesis` builds the mutant
   set internally, predicts, and returns aggregated ISM attributions or
   `SaturationMutagenesisRawResult(y0, y_hat)` when `raw_outputs=True` — per-mutant *values*,
   never the mutant *sequences*. In `design/`, the full {motifs} × {all allowed positions} ×
   {fwd, rc} enumeration is materialized as a real tensor and then discarded except for the
   argmin. **The enumeration exists; the library does not come out.**

So: exhaustive single-substitution mutagenesis as an *in silico assay* — yes; as a *library
specification you can obtain* — no. "partial" is the honest middle.

*Source:* `saturation_mutagenesis.py` L19–35, L85, L205–250; `design/_substitute.py`;
`design/greedy_substitution.py` L188–205; `variant_effect.py` L20–73.

### `sampled_random_mutagenesis` → **partial**
*(new sub-row from v1 `mixed_mutagenesis_one_pool` = no)*

Random variant generation exists and — unusually for this package — **returns sequences**:

- `ersatz.randomize(X, start, end, probs=..., n=1, random_state=None)` replaces a window with
  randomly drawn sequence, *"`n` times for each sequence in X"*, returning
  `shape=(-1, n, len(alphabet), length)` (`ersatz.py` L343–397).
- `utils.random_one_hot(shape, probs=..., random_state=...)` generates random sequences outright.
- `design.screen(model, shape, ..., n_best=k, func=random_one_hot)` samples random sequences in
  batches and keeps the best `k` — the `func=` hook lets a user substitute any generator,
  including a rate-based mutagenizer.

Why not "yes": `randomize`'s `probs` is per-example (`(1, |alphabet|)` or `(len(X), |alphabet|)`),
**not per-position**, so there is **no per-site mutation rate and no notion of sampling variants
of a WT at rate r** — the window is replaced wholesale, which is random-background generation
more than mutagenesis. There is no "sample k variants from the mutational space" operator, no
sampling with or without replacement over an enumerated space, and no record of which sample is
which. Rate-based sampling is reachable only as a user-written `func=`, which is "partial" at
most by the raw-primitives rule.

*Source:* `ersatz.py` L343–397 (signature + docstring re-read this pass); `utils.py` L579–585;
`design/screen.py` L19–47.

### `higher_order_combinatorial` → **partial**
*(new sub-row from v1 `mixed_mutagenesis_one_pool` = no)*

Multiple simultaneous edits in the same sequence are first-class, and one combinatorial axis is
genuinely enumerated:

- `ersatz.multisubstitute(X, motifs, spacing, start=None)` places a **list** of motifs with
  per-gap spacings into one sequence — a higher-order composite edit.
- `space(model, X, motifs, spacing)` takes `spacing` of shape `(n_spacings, n_motifs-1)` —
  *"Each row in this tensor is a different combination of spacings between motifs and each column
  is the spacing between an adjacent pair of motifs"* — and sweeps the whole grid, returning
  `SpaceResult(y_before, y_afters)` with `y_afters` shape `(batch, n_spacings, ...)`. This is the
  one truly enumerated combinatorial axis in the package (Tutorial B3: AP-1 pair cooperativity vs
  distance).
- `variant_effect.substitution_effect` applies **multiple substitution rows to the same example**
  simultaneously; the docstring notes multi-row adjacent positions encode *"entire motifs or just
  multiple characters"*.

Why not "yes": nothing enumerates **combinations of mutations**. `saturation_mutagenesis` is
edit-distance-one only; there is no double-mutant / k-wise combination generator, no
epistasis-style pairwise mutation library, and no combination of a mutation set with itself.
`apply_pairwise` / `apply_product` cross `X` with extra **model inputs**, not with mutations, so
they do not supply the missing axis. Multi-edit sequences must be enumerated by the user.

*Source:* `ersatz.py` L198–260; `space.py` L22, L34–44 + `_validate_input(..., shape=(-1, len(motifs)-1))`;
`variant_effect.py` L41–72; `product.py` L30–60; `saturation_mutagenesis.py` L85; Tutorial B3.

### `heterogeneous_components_one_library` → **no**
*(new sub-row from v1 `mixed_mutagenesis_one_pool` = no; value carried)*

The nearest analogue is real and must be acknowledged: `variant_effect.substitution_effect` takes
a heterogeneous per-example COO edit table `(example_idx, position, alphabet_idx)`, rows grouped
by example — *"one can encode longer variants (e.g., entire motifs or just multiple characters) by
passing in multiple rows with adjacent positions"* — so **one call can mix a SNV on example 0 with
a 10-bp motif swap on example 1**. Likewise `marginalize_annotations` / `ablate_annotations` sweep
a *set* of annotation-defined regions `(example_idx, start, end)` in a single call.

Still "no": there is **no specification object** in which structurally different component types
could coexist. One perturbation *type* per call (substitutions here; indels via
`deletion_effect` / `insertion_effect`; shuffles via `ablate` / `ersatz.shuffle` — never together);
no WT-replicate class; no co-residency of an exhaustive scan with sampled variants with
motif-placed constructs; nothing records which scheme produced which row (indexing is the user's
own row order); and the return is `PerturbationResult` / `PerturbationAnnotationsResult` —
predictions, not a pool of sequences.

*Source:* `variant_effect.py` L20–73; `marginalize.py` L133–205; `ablate.py` L171–215;
`results.py`; full-tree class grep.

### `combinatorial_motif_place` → **partial**
*(unchanged from v1)*

Real combinatorics on two axes: (i) **multiple motifs + spacings** — `multisubstitute` places a
motif list with per-gap spacings, `space` sweeps a 2-D `(n_spacings, n_motifs-1)` grid;
(ii) **motif × position × orientation search** — `design.greedy_substitution` /
`beam_substitution` try every motif at every allowed position each round, `reverse_complement=True`
adds each motif's reverse complement, `input_mask` restricts editable positions.

Why not "yes": `substitute` / `multisubstitute` / `space` take a **scalar** `start`
(`ersatz.py` L97–103, L198–205; `space.py` L34–44), so **user-facing position is not a swept
axis**; and the exhaustive {motifs} × {every allowed position} × {fwd, rc} enumeration that does
exist is a **hidden internal step of an optimizer** (`design/_substitute.py`,
`greedy_substitution.py` L188–205: `_fast_tile_substitute(...)` then `pos = loss_curr.argmin()`)
— only the argmin survives and there is no way to obtain the enumerated set. Permutations of
motif order are not enumerated at all.

*Source:* `space.py`; `ersatz.py` L97–103, L198–205; `design/_substitute.py`;
`design/greedy_substitution.py` L188–205; `_skills/data/references/design.md`.

### `barcode_generation` → **no**
*(unchanged from v1)*

Recursive grep over all 28 package `.py` modules (20 under `tangermeme/`, 6 under
`tangermeme/design/`, 2 under `tangermeme/_skills/`) plus the 15 bundled skill `.md` docs:
**zero hits** for `barcode`, `primer`, `oligo`; re-confirmed this pass over the full tarball
including tests. No minimum-edit-distance code construction, no attachment, no demultiplexing.
The sole `edit distance` hit is `saturation_mutagenesis.py` L85 (single mutants). GC machinery
exists (`utils.gc_content`, `match.py`) but serves background matching.

*Source:* recursive grep (v1 final pass and this pass); `saturation_mutagenesis.py` L85;
`utils.py`; `match.py`.

### `per_sequence_provenance` → **no**
*(unchanged from v1)*

No designed or perturbed sequence carries any record of the operation that made it. Perturbations
return tensors or NamedTuples of model outputs. The only DataFrames the package emits are
`seqlet.recursive_seqlets` and `annotate.annotate_seqlets` — features **discovered in attribution
scores** (out of scope here), not construction records. `AttributionReferencesResult.references`
does hold sequences, but they are DeepLIFT background shuffles, not provenance.

*Source:* `results.py`; `seqlet.py`; `annotate.py`; full-tree class grep.

### `automatic_naming` → **no**
*(unchanged from v1)*

`io.one_hot_to_fasta` is the package's only writer (whole-tree writer grep: every hit is in
`io.py`). Body, `io.py` L577–582:

```python
if headers is None:
    outfile.write("> {}\n".format(i))
else:
    outfile.write("> {}\n".format(headers[i]))
```

The fallback name is the bare integer batch index. No name is ever derived from a design
operation, motif identity, position, or variant.

*Source:* `io.py` L540–590; full-tree writer grep.

### `design_visualization` → **partial**
*(unchanged from v1)*

`plot.plot_logo` (`plot.py` L142ff) accepts an `annotations` DataFrame
(*"motif_name … start … end … strand (optional) … score"*, L222–233) and draws labelled,
non-overlapping boxes under the glyphs (`check_box_overlap` / `place_new_box`, L35–115), with
`annot_cmap`, `score_key`, `n_tracks`, `show_score`. So **a designed construct can be rendered
with its placed motifs named and positioned**. `plot.interactive_logo` (v1.3.0) adds mpld3 hover
tooltips, but needs the optional `[interactive]` extra and is absent from a default install.

Why not "yes": visualization is per-sequence. There is no view of a design or library
*specification*, no pipeline/graph view, no library-level summary, and nothing in `design/` emits
the annotation table — the user must build it. (Attribution-logo plotting as model interpretation
is out of scope and is not counted here.)

*Source:* `plot.py` L35–115, L142–233; `pyproject.toml` optional-dependencies;
`docs/whats_new.rst` v1.3.0; Tutorial C2.

---

## BLOCK B — assay coverage

### `assay_dms` → **no** — *definitional; footnote required*
*(unchanged from v1)*

`saturation_mutagenesis` is explicitly the **in silico** analogue — README: *"This is another form
of attribution method that is conceptually similar to deep mutational scanning but using a
predictive model instead of running an experiment."* The "no" holds only under the row definition
*"designs a DMS library for synthesis"*, and rests on three verified facts: (1) in silico only, no
synthesis output of any kind; (2) nucleotide-level only — zero occurrences of `codon`, `amino`,
`protein`, `orf`, `reading frame` anywhere in the package or docs (re-confirmed this pass); (3)
the mutant sequences are never returned — even `raw_outputs=True` yields
`SaturationMutagenesisRawResult(y0, y_hat)`, per-mutant *values*.

The footnote must state this, because the package ships a module literally named
`saturation_mutagenesis` and `design.md` L48 documents an ISM-style mode
(*"pass `['A','C','G','T']` as `motifs` → greedy single-nucleotide (ISM-style) design"*).

*Source:* README; `saturation_mutagenesis.py` L19–35, L205–250;
`_skills/data/references/design.md` L48.

### `assay_mpra` → **partial** — *generous by design; state the row definition in the caption*
*(unchanged from v1)*

`tangermeme.design` genuinely designs cis-regulatory sequences against a model oracle
(`screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize`; Tutorial B6 designs a
strong-AP-1 sequence against Beluga). `beam_substitution(n_best=k)` and `screen(n_best=k)` return
`shape=(n_best, len(alphabet), length)` — a small ranked **set** of designed sequences. And the
package bulk-generates the standard MPRA scramble controls as materialized sequences
(`randomize` / `shuffle` / `dinucleotide_shuffle`, each with an `n` axis).

Why still "partial": nothing MPRA-specific — no barcodes, adapters, primer sites, oligo-length or
synthesis constraints; no identity or naming for members; no pool object mixing designed members
with their controls; no replicate or collision accounting. This is **regulatory sequence design
against a model oracle, not MPRA library design.**

*Source:* `design/__init__.py`; `design/beam_substitution.py` L173, L305–306;
`design/screen.py` L19–47, L142; `ersatz.py` L343–360, L425–460, L608–672; Tutorial B6; Fig 1A.

### `assay_insilico` → **yes**
*(unchanged from v1)*

The package's stated raison d'être, fully implemented: `marginalize`, `ablate`, `space`,
`saturation_mutagenesis`, `variant_effect.{substitution,deletion,insertion}_effect`,
`marginalize_annotations`, `ablate_annotations`, `apply_pairwise` / `apply_product` — with
built-in batching, device handling and multi-input/multi-output support. Paper p.3: *"These
operations have built-in batching, support for alternative data types and devices, and work
out-of-the-box on multi-input/output models."* Fig 1A "Screening" column.

*Source:* paper abstract and p.3; Fig 1A; module tree.

---

## BLOCK C — genomics integration

### `genome_coordinates` → **yes** *(unchanged)*
`io.extract_loci(loci, sequences, signals=..., chroms=..., in_window=2114, out_window=1000,
max_jitter=0, exclusion_lists=..., min_counts=..., max_counts=..., summits=False, ...)`
(`io.py` L246–265) — *"Either the path to a bed file or a pandas DataFrame object containing three
columns: the chromosome, the start, and the end"*; FASTA via pyfaidx, bigwig via pybigtools; with
exclusion lists, chromosome subsetting, summit-centring, jitter and count thresholds.
`utils.example_to_fasta_coords` maps window-relative spans **back** to genome coordinates by
cross-referencing the originating BED, with zero/one-indexing handling (`utils.py` L233–253) — a
genuine coordinate round-trip. Also `match.extract_matching_loci` (BED + FASTA → GC-matched
negatives).
*Source:* `io.py` L246–265; `utils.py` L233–253; `match.py`; Tutorial C1.

### `transcript_models` → **no** *(unchanged)*
Recursive grep over all modules + skill docs: zero hits for `gtf`, `gff`, `transcript`, `exon`,
`intron` (re-confirmed this pass). `io.py`'s complete public surface is `extract_loci`,
`one_hot_to_fasta`, `read_meme`, `read_vcf`. Genomic input is BED + FASTA + bigwig only.
*Source:* recursive grep; `io.py` public functions.

### `exon_intron_split_codons` → **no** *(unchanged)*
Zero hits for `exon`, `intron`, `codon`, `reading frame`, `orf`, `amino`, `protein`. Sequences are
alphabet-agnostic one-hot tensors (`alphabet: list[str] = ['A','C','G','T']` as a plain default);
there is no coding-biology layer at all.
*Source:* recursive grep, this pass.

### `vcf_vep_output` → **no** *(unchanged)*
`io.read_vcf` is **input-only** (*"returns a pandas DataFrame with the comments filtered out …
only the first 9 columns"*; drops genotype columns, no BCF). The exhaustive whole-tree writer
search finds file-writing code in `io.py` alone, and the only writer is `one_hot_to_fasta`
(FASTA). No VCF writer, no VEP-consumable emitter.
*Source:* `io.py` `read_vcf`; full-tree writer grep.

### `consequence_annotation` → **no** *(unchanged)*
`tangermeme.annotate` means **motif matching of seqlets against a motif database via tomtom-lite**
plus co-occurrence/spacing statistics (`annotate_seqlets`, `count_annotations`,
`pairwise_annotations`, `pairwise_annotations_spacing`); README: *"an annotation is any genomic
span"*. No molecular-consequence logic (stop-gained / synonymous / frameshift / in-frame) exists,
consistent with the zero hits for `codon`, `amino`, `protein`, `orf`.
*Source:* `annotate.py`; README; paper pp.7–8; recursive grep.

---

## BLOCK D — adjacent / complementary

### `primer_design` → **no** *(unchanged)*
Zero hits for `primer`, `oligo`, `melting`, `Tm` across all modules and skill docs (re-confirmed
this pass over the full tarball); the single `adapter` hit is `pisa.py` L174, an `nn.Module` output
adapter. No wet-lab layer in README, docs, tutorials, `whats_new.rst` through v1.4.0, or the paper.
*Source:* recursive grep; `pisa.py` L174; `docs/whats_new.rst`.

### `codon_optimization` → **no** *(unchanged)*
Zero hits for `codon`, `amino`, `protein`, or any translation concept; no codon table, no
translation layer, no expression host. The alphabet-agnosticism claim is the **README's** (*"all
functions are extensible to any alphabet"*) — no protein support is demonstrated and the word
`protein` appears nowhere in the source.
*Source:* recursive grep; README; signatures across `ersatz.py`, `utils.py`, `io.py`.

### `synthesis_constraints` → **no** *(unchanged)*
Zero hits for `homopolymer`, `restriction`, `synthesis`, `secondary structure`, or any repeat /
oligo-length constraint check. The GC machinery serves **background selection**, not
synthesizability: `utils.gc_content`; `match.extract_matching_loci`, whose stated rationale is the
GC "mirage" effect (*"Uniformly random sequences are much higher in GC content than real genomic
DNA …"*, `motif-effects.md`) and whose other filters (`max_n_perc`, optional `bigwig` +
`signal_beta` / `signal_threshold`, `match.py` L285–395) reject *negatives carrying too much
signal*. Still negative-set selection, not synthesis QC.
*Source:* recursive grep; `utils.gc_content`; `match.py` L285–395; `motif-effects.md`.

### `degenerate_iupac_codons` → **no** — *NEW ROW*

Checked directly this pass over the full tree (55 `.py` incl. tests, 17 `.md`, `docs/`,
`pyproject.toml`): **zero hits for `iupac`**; `codon` zero hits; the two `degener` hits are in
`tests/test_saturation_mutagenesis.py` and mean a *degenerate (zero-width / single-position) input
span*, not degenerate bases. `utils.pwm_consensus` explicitly resolves a PWM by **argmax** —
*"the consensus, which is the most likely individual sequence … by taking the argmax at each
position"* — so even the PWM→sequence path produces a single concrete sequence, never an ambiguity
code. There is no expansion of a degenerate string into members and no compression of members into
one.

Two near-misses, neither of which is degenerate-codon design: (1) encoding-level handling of `N` /
`[0,0,0,0]` as an ambiguous genomic position (`whats_new.rst` L513, `utils.characters` /
`_validate_input`); (2) `design.greedy_marginalize`, which returns a variable-width one-hot
*construct* carrying `N` where motif contributions cancel — an output artefact of the optimizer,
not a specifiable degenerate design.

*Source:* own recursive grep this pass; `utils.py` L789–800 (`pwm_consensus` docstring);
`docs/whats_new.rst` L513; final memo §5.

### `negative_control_generation` → **yes** — *NEW ROW*

This is a first-class, documented, sequence-returning feature set — one of the very few places the
package hands back sequences rather than model outputs:

- `ersatz.dinucleotide_shuffle(X, start=0, end=-1, n=20, ...)` → documented return
  `shape=(-1, n, len(alphabet), length)`. Dinucleotide shuffles are *the* standard MPRA/CRE
  scramble control.
- `ersatz.shuffle(X, start=0, end=-1, n=1, random_state=...)` → same n-axis, mononucleotide.
- `ersatz.randomize(X, start, end, probs=..., n=1, ...)` → n random-substitution controls per
  sequence, composition controllable via `probs`.
- `utils.reverse_complement` (`utils.py` L507) → reverse-complement controls.
- `utils.random_one_hot` → random background sequences.
- `match.extract_matching_loci` → **GC-matched genomic negatives** sampled genome-wide, with
  `max_n_perc` and optional bigwig `signal_beta` / `signal_threshold` rejection of
  signal-carrying candidates (`match.py` L285–395). The package documents *why* this matters (the
  GC "mirage" effect, `motif-effects.md`).
- `ablate(model, X, start, end)` = shuffle-out of a window as a built-in negative operation;
  Tutorial A1 and Tutorial B2 demonstrate the whole set.

Caveat recorded for honesty (does not reduce the value): the controls are anonymous tensors — no
label marks a sequence as a control, and they cannot be declared as members of a pool alongside
non-control members. The row asks whether controls are a first-class feature; here they are.

*Source:* `ersatz.py` L343–397, L425–470, L608–672; `utils.py` L507, L579–585; `ablate.py`;
`match.py` L285–395; Tutorials A1, B2; final memo §3 `assay_mpra` and §5.

### `ml_model_in_loop` → **yes** — *NEW ROW*

Design driven by a predictive model's output is the entire content of `tangermeme.design`
(`__all__` = `screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize`):
each takes `model`, a `loss`, an optional target `y`, and iterates edits to minimise the loss.
`greedy_substitution` tiles each motif at every allowed position (both strands), predicts, and
keeps the argmin, for `max_iter` rounds; `beam_substitution` maintains a beam and returns
`n_best` designs; `screen` samples random sequences and keeps the `n_best` by model loss;
`greedy_marginalize` returns a construct of the motifs whose marginal effect survives. Fig 1A's
"Design" column lists Greedy Substitution, Motif Implantation, Construct Marginalization, and
Ledidi; the docs point at `ledidi` for the gradient-based variant
(`design.md`: *"`tangermeme.design` is discrete and greedy. When you want gradient-based, minimal
edits … reach for the `ledidi` library."*). Tutorial B6 is a full worked example against Beluga.

Note for the comparison: the model-in-loop output is a small ranked set of winners, not a library
— this row credits the mechanism, not library-scale output.

*Source:* `design/__init__.py`; `design/greedy_substitution.py` L63, L188–205;
`design/beam_substitution.py` L173, L305–306; `design/screen.py` L19–47, L142;
`_skills/data/references/design.md`; Tutorial B6; paper Fig 1A.

### `readout_analysis` → **no** — *NEW ROW*

No analysis of sequencing readout from a built library exists. Checked this pass: **zero hits** for
`fastq`, `demultiplex`, `count table`, `barcode` anywhere in the tree (the one `readout` hit is
README prose about a model predicting "a readout of interest"). The package's readers are
`extract_loci` (BED/FASTA/bigwig), `read_meme`, `read_vcf` — reference and model inputs, not assay
output. Everything downstream (`seqlet`, `annotate`, `kmers`, statistics) analyses **model
attributions**, not experimental counts, and is out of scope for this row regardless.

*Source:* own recursive grep this pass; `io.py` public functions; `seqlet.py`; `annotate.py`.

---

## BLOCK E — engineering and availability

*(row value is "yes" for all five; the informative answer is the evidence, and the rendered table
cell should print the evidence string, not "yes")*

### `interface` → **yes** — render as: **"Python API only (no CLI, no GUI, no web service)"**
`pyproject.toml` `[project.scripts]` contains **exactly one** entry —
`tangermeme-install-skills = "tangermeme._skills.install:main"` — which only installs the bundled
Claude Code Agent Skill. There is no `[project.entry-points]` group, no web service, no GUI.
Analysis CLIs were removed upstream (README: *"These FIMO and Tomtom command-line tools have been
moved to memesuite-lite … Please use those!"*). Notable extra surface: a bundled **Claude Code
Agent Skill** (`SKILL.md` + 14 reference docs) so an LLM agent can drive the API.
*Source:* `pyproject.toml` (re-read this pass); README; `tangermeme/_skills/` tree.

### `license` → **yes** — render as: **"MIT"**
Verified three ways: `pyproject.toml` → `license = "MIT"`, `license-files = ["LICENSE"]`;
PyPI JSON → `license_expression: MIT`; GitHub API → `license.spdx_id: "MIT"`. Paper p.9:
*"tangermeme is free and open source software available under the MIT license at
https://github.com/jmschrei/tangermeme."*
*Source:* `pyproject.toml`; PyPI JSON; `gh api repos/jmschrei/tangermeme`; paper p.9.

### `installable_today` → **yes** — render as: **"pip install tangermeme (PyPI v1.4.0, 2026-06-25); Python >=3.10, torch>=2.0"**
On PyPI with wheels, hatchling build, `Development Status :: 5 - Production/Stable`, classifiers
for Python 3.10–3.13. Twelve pure-pip runtime dependencies (numpy, scipy, pandas, torch>=2.0,
pyfaidx, pybigtools, tqdm, numba, joblib, scikit-learn, matplotlib, memelite) — no compiled
external tool, no MEME-suite install, no Docker requirement. Optional extras: `dev`, `docs`,
`interactive` (mpld3; note `plot.interactive_logo` is therefore **absent from a default install**).
No conda-forge package observed. GitHub repo is public and `archived: false`.
*Source:* `pyproject.toml` (re-read this pass); PyPI JSON; GitHub API.

### `last_activity` → **yes** — render as: **"v1.4.0 released 2026-06-25; repo pushed 2026-07-19 (8 releases in 12 months)"**
PyPI latest release **1.4.0, uploaded 2026-06-25T18:05:55**; GitHub `pushed_at`
**2026-07-19T18:34:08Z**, `updated_at` 2026-08-07T19:42:00Z, `archived: false`, 308 stars / 32
forks / 8 open issues. Release ladder: 1.0.0 (2025-08-27), 1.0.2 (2026-01-19), 1.0.3 (2026-02-15),
1.0.4 (2026-04-23), 1.1.0 and 1.2.0 (2026-05-27), 1.3.0 (2026-06-23), 1.4.0 (2026-06-25) — **8
releases in 12 months**. CI unit-test badge; Read the Docs live and matching the repo tree.
*Source:* PyPI JSON and GitHub API (fetched live in the v1 final pass); tree re-downloaded this
pass and consistent.

### `documented_examples` → **yes** — render as: **"20 executable notebooks (14 tutorials + 3 vignettes + 1 how-to + 2 paper figures), 18 API pages, 14 agent-skill reference docs"**
Enumerated directly from `docs/` this pass: `docs/tutorials/` = 14 notebooks (A1 Ersatz Sequence
Manipulation, A2 Predictions, A3 DeepLIFT/SHAP, A4 Seqlets, A5 Annotations, B1 Marginalization,
B2 Ablation, B3 Spacing, B4 Saturation Mutagenesis, B5 Variant Effect, B6 Design, B7 Cartesian
Product, C1 IO and Data Loading, C2 Plotting); `docs/vignettes/` = 3; `docs/howto/` = 1
(*How To — Reduce Friction and Save Time with Tangermeme*); `docs/paper/` = 2 figure notebooks;
`docs/api/` = 18 module reference pages; `tangermeme/_skills/data/references/` = 14 reference `.md`
docs beside `SKILL.md`. Design-relevant subset: A1, B1, B2, B3, B4, B5, B6, B7, C1.
*Source:* `docs/` tree enumerated this pass; `tangermeme/_skills/data/references/` listing.

---

## Changes vs. the v1 record

| v2 key | v1 key | v1 → v2 | Why |
|---|---|---|---|
| `library_first_class_object` | `library_as_object` | partial → **partial** | Split. Batch-first in-memory tensors, held/inspected/transformed/passed, never file-only — that credit stays here. Still no library class. |
| `library_algebra` | `library_as_object` | partial → **no** | Split. Nothing in tangermeme combines/samples/repeats libraries; assembly is the user's own `torch.cat`, i.e. outside the tool. |
| `composable_operations` | `dag_chaining` | partial → **partial** | Rename; value carried. New wording drops the graph-object requirement, so "yes" is arguable — held at partial because design-step composition is user-assembled primitives. |
| `lazy_generation` | `lazy_evaluation` | partial → **partial** | Rename; value and evidence carried unchanged. |
| `exhaustive_single_scans` | `mixed_mutagenesis_one_pool` | no → **partial** | Split. `saturation_mutagenesis` really is an every-position, every-substitution scan; substitutions only, and the sequences are never returned. |
| `sampled_random_mutagenesis` | `mixed_mutagenesis_one_pool` | no → **partial** | Split. `randomize(n=)` / `random_one_hot` / `screen` sample random sequences and return them; no per-site rate, no sampling over an enumerated variant space. |
| `higher_order_combinatorial` | `mixed_mutagenesis_one_pool` | no → **partial** | Split. `multisubstitute` + `space`'s `(n_spacings, n_motifs-1)` grid are genuine multi-edit combinatorics; combinations of *mutations* are never enumerated. |
| `heterogeneous_components_one_library` | `mixed_mutagenesis_one_pool` | no → **no** | Split; value carried. The COO edit table is heterogeneous per example, but there is no specification object and one perturbation type per call. |
| — | `hgvs_input` | no → **dropped** | Row removed from v2. |
| `degenerate_iupac_codons` | — | **no** (new) | Zero `iupac` hits; `pwm_consensus` is argmax; `N` handling is encoding-level only. |
| `negative_control_generation` | — | **yes** (new) | `dinucleotide_shuffle` / `shuffle` / `randomize` (n-axis, sequences returned), `reverse_complement`, `ablate`, GC-matched negatives via `match`. |
| `ml_model_in_loop` | — | **yes** (new) | The whole `tangermeme.design` subpackage: `screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize`. |
| `readout_analysis` | — | **no** (new) | Zero `fastq` / `demultiplex` / count-table machinery. |
| `installable_today` | — | **yes** (new) | pip, PyPI v1.4.0, pure-pip deps, Python 3.10–3.13. |
| `last_activity` | — | **yes** (new) | v1.4.0 2026-06-25; repo pushed 2026-07-19; 8 releases/12 months. |
| `documented_examples` | — | **yes** (new) | 20 executable notebooks + 18 API pages + 14 agent-skill docs. |

All Block B, C and D-1..D-3 values carry over from the v1 final record unchanged, with evidence
re-verified where cited above.

## Bottom line (unchanged from v1, re-stated for v2)

The v2 split makes tangermeme's real shape visible rather than diminishing it: it scores
**partial across four of the mutagenesis sub-rows and yes on both new mechanism rows**
(`negative_control_generation`, `ml_model_in_loop`) — more credit than the single v1 row gave it.
The boundary is equally visible: `library_algebra` = no, `heterogeneous_components_one_library` = no,
`per_sequence_provenance` = no, `automatic_naming` = no. **tangermeme's composition target is a
model call, not a library.** Every pipeline terminates in predictions or attributions; the only
sequences it hands back are DeepLIFT backgrounds, generated controls, and design winners — all
anonymous tensors. The honest framing remains **complementary, not competing**: PoolParty specifies
and materializes the library; tangermeme is one of its most natural consumers.
