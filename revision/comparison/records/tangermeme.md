# tangermeme — FINAL capability record for the PoolParty tool comparison

**Tool:** tangermeme
**Slug:** `tangermeme`
**Citation key:** `Schreiber2025nd`
**Tier:** 1
**Author:** Jacob Schreiber (IMP Vienna / UMass Chan)
**Self-description (pyproject / PyPI):** "Biological sequence analysis for the modern age."
**Paper self-description:** "a highly optimized toolkit for 'everything-but-the-model' when it comes to genomic deep learning"
**Version at survey:** v1.4.0 (PyPI, 2026-06-25); repo `pushed_at` 2026-07-19
**Status of this document:** final. Extraction (2026-08-10) → adversarial review (2026-08-10) → this record.
All 24 capability **values** survive adversarial review unchanged; what changed is **evidence and
citations**. Every quote below was re-verified from primary source for this final pass
(see §1 for method). **No unresolved extractor/reviewer disagreements** — see §7.

---

## 1. Sources consulted and verification method for this final pass

| Kind | Reference |
|---|---|
| PDF (local) | `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/../../../lit_review/analyzed/Schreiber2025_tangermeme_all.pdf` — bioRxiv 2025.08.08.669296, 13 pp., extracted with PyMuPDF |
| Repo (full tree, this pass) | `gh api repos/jmschrei/tangermeme/git/trees/main?recursive=1` — all paths enumerated |
| Repo source (this pass) | **All 28 `.py` files** and **all 15 bundled skill `.md` docs** downloaded from `raw.githubusercontent.com/jmschrei/tangermeme/main/` and grepped locally |
| Packaging | `pyproject.toml` (raw, this pass) |
| PyPI | `https://pypi.org/pypi/tangermeme/json` (fetched live, this pass) |
| GitHub metadata | `gh api repos/jmschrei/tangermeme` (fetched live, this pass) |
| Docs site | https://tangermeme.readthedocs.io/en/latest/ ; `docs/` tree enumerated via the trees API |

**CITATION CORRECTION — important, applies to nine "no" cells.** The extraction cited
*"grep over `tangermeme/*.py` (20 files)"*. That shell glob silently excludes
`tangermeme/design/` (6 modules) and `tangermeme/_skills/` (2 modules) — the design
subpackage being exactly where a referee would expect a synthesis/barcode helper to live.
The reviewer flagged this and re-ran recursively, citing "25 modules"; **that count is also
wrong**. The verified figure, from the trees API this pass, is:

> **28 `.py` files = 20 directly under `tangermeme/` + 6 under `tangermeme/design/`
> (`__init__`, `_substitute`, `beam_substitution`, `greedy_marginalize`,
> `greedy_substitution`, `screen`) + 2 under `tangermeme/_skills/` (`__init__`, `install`),
> plus 15 bundled skill `.md` docs (`SKILL.md` + 14 reference files).**

I re-ran the negative greps myself over all 43 of those files. **Result identical to both
prior passes: zero hits** for `barcode`, `primer`, `oligo`, `codon`, `amino`, `protein`,
`orf`, `reading frame`, `hgvs`, `gtf`, `gff`, `transcript`, `exon`, `intron`, `synthesis`,
`homopolymer`, `restriction`, `melting`, `secondary structure`. The only `adapter` hit is
`pisa.py:174` ("adapter that transposes or flattens the output" — an `nn.Module` wrapper).
The only `edit distance` hit is `saturation_mutagenesis.py:85`. **No verdict changes; use
the "28 modules including `design/` and `_skills/`" citation everywhere.**

**Structural fact worth having in hand (verified this pass, stronger than anything in either
prior document):** `grep -rn "^class " --include=*.py` over the whole package returns
**exactly six classes, five of which are `NamedTuple`s** (four of model outputs; the fifth,
`AttributionReferencesResult`, pairs attributions with the reference *sequences* used):

```
space.py:22                    class SpaceResult(NamedTuple)
saturation_mutagenesis.py:19   class SaturationMutagenesisRawResult(NamedTuple)
results.py:16                  class PerturbationResult(NamedTuple)
results.py:28                  class PerturbationAnnotationsResult(NamedTuple)
results.py:41                  class AttributionReferencesResult(NamedTuple)
utils.py:18                    class TangermemeWarning(UserWarning)
```

There is **no class of any kind representing a sequence, a set of sequences, a library, or a
design**. That single grep is the cleanest possible support for the Block-A cells.

---

## 2. What tangermeme is

Paper p.3: *"we introduce tangermeme, a software package which implements 'everything-but-the-model'
when it comes to genomic deep learning. … tangermeme focuses on using trained models to drive
genomics research (Fig 1A), including support for sequence manipulations, model operations,
and efficient combinations of the two (e.g. variant effect prediction)."*

Fig 1A capability grid, tangermeme column: **Design** — Greedy Substitution, Motif Implantation,
Construct Marginalization, Ledidi; **Screening** — Sequence Manipulation, Stacked Operations,
Cartesian Products, Saturation Mutagenesis, Marginalizations, Ablations, Spacings, Variant
Effect Prediction; **Seqlets** — Calling, Annotation, Ablation+Marginalization, Annotation
Counting, Spacing; **Core** — DeepLIFT/SHAP, PISA, Plotting. The same figure marks Model
Training, Data Preprocessing and Model Zoo as **not** provided (unlike Selene, EUGENe, gReLU,
CREsted).

**Data model:** every sequence collection is a bare `torch.Tensor` — shape
`(-1, len(alphabet), length)`, or `(-1, n, len(alphabet), length)` for the replicate-generating
`ersatz.randomize` / `shuffle` / `dinucleotide_shuffle`. See the class grep above — there is no
library, pool, or record type.

---

## 3. Capability-by-capability record

Format: **key → value.** Evidence, then source. Corrections relative to the original
extraction are marked **[CORRECTED]** or **[EVIDENCE STRENGTHENED]**.

### BLOCK A — library specification

---

**`library_as_object` → partial**

Operations are batch-first, so a user does **not** write a loop to apply *one* operation to
many sequences: `substitute(X, motif)` accepts a whole batch, and *"If a motif with a batch
size equal to that of X is provided, there will be 1-1 correspondence between the motifs and
the sequence"* (`ersatz.py` docstring). That batch-first ergonomic is the entire basis for
"partial".

But the batch is a raw tensor. **[CORRECTED]** The extraction said *"`results.py` … all four
of its NamedTuples hold model outputs, not sequences"*. Two errors: (a) `results.py` holds
**three** NamedTuples (`PerturbationResult`, `PerturbationAnnotationsResult`,
`AttributionReferencesResult`); `SpaceResult` is defined in `space.py:22` and
`SaturationMutagenesisRawResult` in `saturation_mutagenesis.py:19` — five in total across the
package, which is likely where the miscount came from (the v1.2.0 changelog lists all five
together). (b) "all hold model outputs, not sequences" is not exact:
`AttributionReferencesResult.references` **is** a tensor of sequences (the DeepLIFT background
shuffles).

The accurate and much stronger statement: **the entire package defines six classes; five are
NamedTuples and one is a Warning subclass; none represents a library, and none carries
identity, name, or metadata for any sequence.** Building a library with several different
design intents requires the user to `torch.cat` tensors themselves and track membership
outside tangermeme.

*Source:* full-tree `grep "^class "` (this pass); `results.py` L16/L28/L41; `space.py` L22;
`saturation_mutagenesis.py` L19; `ersatz.py` `substitute` docstring.
*Note:* "no" would also be defensible. "partial" is the deliberately generous reading so the
row cannot be attacked as understating a competitor.

---

**`dag_chaining` → partial**

There is a genuine, deliberate composition mechanism — the `func=` plug-point.
`_skills/data/references/func-pattern.md` L11: the contract is

```
func(model, X, args=None, **kwargs) -> torch.Tensor | list[torch.Tensor]
```

accepted by `ablate`, `marginalize`, `space`, `variant_effect.*`, `product.apply_pairwise`
and `product.apply_product`. Nesting is documented and demonstrated (`product.md`, Tutorial B7):

```python
apply_pairwise(marginalize, model, X, args=(cell_states,),
               additional_func_kwargs={'func': deep_lift_shap}, motif="TGACGTCA")
```

Paper p.3 (verbatim, verified in the PDF): *"An important consequence of tangermeme's design
is that operations can be stacked."* Fig 1A calls it "Stacked Operations".

**Why partial:** what is composed is *sequence-manipulation → model-operation*. It is plain
Python function nesting — no graph object, no node/edge inspection, no reuse of a node across
branches — and **every leaf is a model call returning predictions or attributions, never a
set of sequences.**

*Source:* `_skills/data/references/func-pattern.md`; `_skills/data/references/product.md`;
`product.py`; paper p.3; Tutorial B7.

---

**`lazy_evaluation` → partial**

`product.md` L50, verbatim: *"Batches are built iteratively, so the full product is never
materialized in memory."* Confirmed in `product.py`: `for x in tqdm(itertools.product(X, *args))`
accumulating to `batch_size` before each call. `io.extract_loci` documents `max_jitter`
expansion to *"reduce the memory footprint"*.

Counter-evidence (all verified): `ersatz.*` returns fully materialized tensors;
`saturation_mutagenesis.py` L212–235 eagerly builds `X_ = X[i].repeat(n_edits, 1, 1)` per
example; `design/greedy_substitution.py` L198 materializes
`X_ = X.float().repeat(input_idxs.shape[0], 1, 1)` every motif, every round.
**Streaming exists only on the Cartesian-product axis, never as a library abstraction.**

*Source:* `_skills/data/references/product.md` L50; `product.py`;
`saturation_mutagenesis.py` L212–235; `design/greedy_substitution.py` L188–205; `io.py`.

---

**`mixed_mutagenesis_one_pool` → no**

No specification object holds several mutagenesis schemes with provenance. Each entry point
implements one perturbation type per call and returns a homogeneous result of *model outputs*.

**[EVIDENCE STRENGTHENED — two things the original evidence omitted that a referee would exploit:]**

1. **`marginalize_annotations` (`marginalize.py` L133) and `ablate_annotations`
   (`ablate.py` L171)** sweep a *set* of annotation-defined regions in one call.
   `marginalize_annotations(model, X, X0, annotations, **kwargs)` takes
   `annotations: torch.Tensor, shape=(n_annotations, 3)` = `(example_idx, start, end)`, and
   its body is literally `for idx, start, end in annotations: seq = X[idx, :, start:end]...;
   marginalize(model, X0, seq, **kwargs)` — i.e. an **n_annotations × n_backgrounds sweep in
   one call**, returning `PerturbationAnnotationsResult(y_befores, y_afters)` stacked per
   annotation. This is the closest thing tangermeme has to a per-element indexed sweep.
2. **`variant_effect.substitution_effect` accepts a heterogeneous per-example COO edit table.**
   `substitutions` is `shape=(-1, 3)` = `(example_idx, position, alphabet_idx)`; rows are
   grouped by example and applied together. Docstring L41–44, verbatim: *"The substitutions
   provided must be individual variants, i.e., each row in the tensor corresponds to a single
   substitution in a single example, but one can encode longer variants (e.g., entire motifs
   or just multiple characters) by passing in multiple rows with adjacent positions."*
   So **one call can mix a SNV on example 0 with a 10-bp motif swap on example 1.**

**Still "no", and here is precisely why:** substitutions only (indels and shuffles require
separate calls — `deletion_effect`, `insertion_effect`, `ablate`); one perturbation *type* per
call; no WT-replicate or sampled-random scheme co-resident; nothing records which scheme
produced which row (indexing is by the user's own row order); and the return is
`PerturbationResult(y_before, y_after)` / `PerturbationAnnotationsResult` — predictions, not a
pool of sequences.

*Source:* `marginalize.py` L133–205; `ablate.py` L171–215; `variant_effect.py` L20–73;
`results.py`.

---

**`combinatorial_motif_place` → partial**

Real combinatorics exist along two axes:

- **Multiple motifs + spacings.** `multisubstitute(X, motifs, spacing, start=None)` places a
  motif list with per-gap spacings; `space(model, X, motifs, spacing)` takes `spacing` of
  shape `(n_spacings, n_motifs-1)` — *"Each row in this tensor is a different combination of
  spacings between motifs and each column is the spacing between an adjacent pair of motifs"*
  — returning `SpaceResult(y_before, y_afters)` with `y_afters` shape `(batch, n_spacings, ...)`.
  This is a user-specified, output-exposed enumerated spacing axis.
- **Motif × position × orientation search.** `design.greedy_substitution` /
  `design.beam_substitution`: each round tries substituting every motif at candidate positions
  and keeps the best edit; `reverse_complement=True` also considers each motif's reverse
  complement; `input_mask` restricts the positions at which a motif may *start*, not the
  positions it may overwrite — a motif substituted at an allowed start runs over masked-out
  positions downstream of it (`design/greedy_substitution.py` L92–97).

**Why "partial" and not "yes":** `substitute` / `multisubstitute` / `space` take a **scalar**
`start` (`ersatz.py` L97–103, L198–205; `space.py` L34–44 — `start: int | None = None` in all
three), so **position is not a swept axis**. And `apply_product` / `apply_pairwise` cross `X`
with **extra model inputs** (`args` — the worked example is DragoNNFruit cell-state × read-depth),
not with motifs, positions or orientations.

**[CORRECTED — the original wording "No API enumerates {motifs} × {positions} × {orientations}
into an output sequence set" is true of outputs but reads as if the enumeration does not exist;
it does, and the author would say so.]** The correct, un-attackable phrasing:

> **The exhaustive {motifs} × {every allowed position} × {fwd, rc} enumeration IS built and
> materialized as a real sequence tensor — it is simply a hidden internal step of an optimizer.
> Only the winner of each round survives (the argmin in `greedy_substitution`, the `beam_size`
> best in `beam_substitution`), and there is no way to obtain the enumerated set.**

Verified in `design/_substitute.py`, whose entire docstring is *"This function takes a motif
and inserts it at all possibilities"*, driven from `design/greedy_substitution.py` L186–214
(abridged below, with the explanatory comments added here and the reverse-complement
augmentation shown at its real location, L161–162):

```python
for idx, motif in enumerate(tqdm(motifs, ...)):          # motifs already += [rc(m) for m in motifs]
    input_idxs = torch.where(input_mask_ == True)[0].numpy()
    X_ = X.float().repeat(input_idxs.shape[0], 1, 1).numpy(force=True)
    _fast_tile_substitute(X_, motif_ohe, input_idxs)     # motif at EVERY allowed position
    y_hat = predict(model, X_, ...)
    pos = loss_curr.argmin()                              # only the winner escapes
```

*Source:* `space.py` docstring + `_validate_input(..., shape=(-1, len(motifs)-1))`;
`ersatz.py` L97–103, L198–205; `product.py` L30–60 docstring;
`design/_substitute.py`; `design/greedy_substitution.py` L188–205;
`_skills/data/references/design.md`.

---

**`barcode_generation` → no**

**[CITATION CORRECTED]** Recursive grep over **all 28 `.py` modules (including `design/` and
`_skills/`) plus the 15 bundled skill `.md` docs**: **zero hits** for `barcode`, `primer`,
`oligo`. The sole `edit distance` hit is `saturation_mutagenesis.py` L85 —
*"…of the sequences with an edit distance of one on them"*, i.e. single mutants, not barcode
separation. No demultiplexing, no minimum-distance code construction, no barcode attachment.

GC machinery exists but for an unrelated purpose: `utils.gc_content` ("Compute the GC content
of one-hot encoded sequences") and `match.py` (genome-wide GC-content calculation and sampling
of GC-matched negatives).

*Source:* own recursive grep, this pass (§1); `saturation_mutagenesis.py` L85; `utils.py`;
`match.py`.

---

**`per_sequence_provenance` → no**

Every perturbation returns tensors or NamedTuples of model outputs. The DataFrames the package
produces are seqlet tables (`seqlet.recursive_seqlets` — example_idx, start, end, attribution,
p-value — and `seqlet.tfmodisco_seqlets`), matched-locus coordinates
(`match.extract_matching_loci`), a coordinate remapping (`utils.example_to_fasta_coords`) and a
parsed VCF (`io.read_vcf`) — features **discovered in attribution scores**, genomic intervals,
or passed-through input, never how a sequence was constructed. (`annotate.annotate_seqlets`
returns two tensors of motif indexes and p-values, not a DataFrame; its `-> tuple[pandas.DataFrame, list[str]]`
annotation at `annotate.py` L25 is stale.) No designed or perturbed sequence carries any record
of the operation that made it.

**[CORRECTED, minor]** As under `library_as_object`: `AttributionReferencesResult.references`
does hold *sequences* — but they are DeepLIFT background shuffles, not provenance.

*Source:* `results.py`; `seqlet.py`; `annotate.py`; full-tree class grep.

---

**`automatic_naming` → no**

**[EVIDENCE STRENGTHENED — verified in source, not just the docstring.]** `io.one_hot_to_fasta`
is the package's **only sequence writer** (confirmed this pass by
`grep -rn "\.write(\|to_csv\|open(...'w')" --include=*.py` over the whole tree — every hit is
in `io.py`). Its body, `io.py` L577–582:

```python
if headers is None:
    outfile.write("> {}\n".format(i))
else:
    outfile.write("> {}\n".format(headers[i]))
```

The fallback name is the bare integer batch index. No name is ever derived from a design
operation, motif identity, position, or variant.

*Source:* `io.py` L540–590; full-tree writer grep, this pass.

---

**`design_visualization` → partial**

**[EVIDENCE CORRECTED — the original cited only the v1.3.0 `interactive_logo` and missed the
stronger static case.]**

`plot.plot_logo` (`plot.py` L142ff) itself takes an `annotations` DataFrame and draws labelled,
non-overlapping boxes under the glyphs. Verbatim parameter doc (`plot.py` L222–233):

> *"annotations: pandas.DataFrame, optional — A set of annotations with the following columns
> in any order except for `motif_name`, which can be called anything but must come first:
> motif_name … start: the start of the hit relative to the window provided … end … strand
> (optional) … score. These will probably come from the output of the hit caller."*

with `annot_cmap`, `score_key`, `n_tracks`, `show_score` controls and helper functions
`check_box_overlap` / `place_new_box` (`plot.py` L35–115) for non-overlapping placement. So a
**designed construct CAN be rendered with its placed motifs named and positioned** — the user
just has to supply the table themselves.

`plot.interactive_logo` (v1.3.0) adds mpld3 hover tooltips over the same annotation model —
but note it requires the **optional `[interactive]` extra** (`pyproject.toml`:
`interactive = ["mpld3>=0.5"]`). The function itself ships and imports in a default install
(`plot.py` L578); it imports `mpld3` lazily and **raises `ImportError` when called** unless
the extra is present.

**Why still "partial":** there is no view of a design or library *specification*, no graph /
pipeline view, and no library-level summary plot. Visualization is per-sequence, and the
annotation table must come from outside the design step (nothing in `design/` emits one).

*Source:* `plot.py` L35–115, L142–233, L578ff; `pyproject.toml` `[project.optional-dependencies]`;
`docs/whats_new.rst` v1.3.0; Tutorial C2.

---

### BLOCK B — assay coverage

---

**`assay_dms` → no** — *definitional risk; footnote required (see §7).*

`saturation_mutagenesis` is explicitly the **in silico** analogue. README: *"This is another
form of attribution method that is conceptually similar to deep mutational scanning but using
a predictive model instead of running an experiment."*

**[FOOTNOTE THE DEFINITION EXPLICITLY — the author is a likely referee and his package
literally contains a module named `saturation_mutagenesis` plus a documented ISM-style design
mode:** `design.md` L48 — *"pass `['A', 'C', 'G', 'T']` as `motifs` → greedy single-nucleotide
(ISM-style) design"`.**]** The "no" holds only under the row definition *"designs a DMS
library for synthesis"*, and rests on three verified facts:

1. **In silico only** — no synthesis output of any kind.
2. **No coding-biology layer** — zero occurrences of `codon`, `amino`, `protein`, `orf`,
   `reading frame` anywhere in the 28 modules or 15 skill docs (own grep, this pass). The
   ISM primitive itself is channel-generic — it enumerates whatever alphabet axis `X` carries
   (`saturation_mutagenesis.py` L222–232) — so the grep establishes the absence of translation
   machinery, not a hard nucleotide restriction.
3. **The mutant sequences are never returned.** `saturation_mutagenesis.py` L212–245 builds
   `X_ = X[i].repeat(n_edits, 1, 1)` internally, predicts, and returns either aggregated
   ISM attributions or `SaturationMutagenesisRawResult(y0, y_hat)` when `raw_outputs=True` —
   i.e. even the raw escape hatch gives you per-mutant *values*, never the mutant sequences.

*Source:* README; `saturation_mutagenesis.py` L19–35, L205–250;
`_skills/data/references/design.md` L48; own recursive grep.

---

**`assay_mpra` → partial** — *definitional note; kept generous (see §7).*

`tangermeme.design` genuinely designs cis-regulatory sequences against a model oracle:
`screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize`
(`design/__init__.py` `__all__`). Tutorial B6 designs a sequence with strong AP-1 binding
using Beluga. Fig 1A's "Design" column is real.

**[EVIDENCE CORRECTED — two overreaches in the original that a referee would catch:]**

1. **"no notion of controls/scrambles as library members" overreaches.** tangermeme *does*
   bulk-generate scramble and random controls as materialized sequences:
   - `ersatz.randomize(X, start, end, probs=..., n=1, random_state=None)` — docstring L353–356:
     *"It will do this `n` times for each sequence in X and so return a tensor with one more
     dimension than `X`."*
   - `ersatz.shuffle(X, start=0, end=-1, n=1, ...)` — same n-axis.
   - `ersatz.dinucleotide_shuffle(X, start=0, end=-1, n=20, ...)` — documented return
     `shape=(-1, n, len(alphabet), length)`. These functions return sequences rather than
     model outputs.
   The accurate statement is: **the control sequences are generated, but they carry no
   identity, no name, and no membership in a pool, and cannot be co-resident with
   motif-placed members in one specification.**
2. **"it designs one sequence at a time" is true only of `greedy_substitution` /
   `greedy_marginalize`** (`greedy_substitution` docstring: `X: torch.tensor,
   shape=(1, len(alphabet), length)`). `beam_substitution(..., n_best=k)` returns
   `X: torch.Tensor, shape=(n_best, len(alphabet), length)`
   (`design/beam_substitution.py` L173, L305–306: `n_best = min(n_best, len(beam));
   return torch.cat([seq for _, seq in beam[:n_best]], dim=0)`), and `screen(..., n_best=k)`
   does likewise (`design/screen.py` L142, L183). **The design layer can return a small ranked
   SET of designed sequences.**

**Why still "partial", not "yes":** nothing MPRA-specific — no barcodes, adapters, primer
sites, oligo-length or synthesis constraints; no identity or naming for members; no pool
object mixing designed members with their controls; no replicate/collision accounting. This is
**regulatory sequence design against a model oracle, not MPRA library design.**

*Source:* `design/__init__.py`; `design/beam_substitution.py` L25–43, L129, L173, L305–306;
`design/screen.py` L27, L96, L142, L183; `design/greedy_substitution.py` L63;
`ersatz.py` L343–360, L425–460, L608–672; Tutorial B6; paper Fig 1A.

---

**`assay_insilico` → yes**

This is the package's stated raison d'être and is fully implemented: `marginalize`, `ablate`,
`space`, `saturation_mutagenesis`, `variant_effect.{substitution,deletion,insertion}_effect`,
`marginalize_annotations`, `ablate_annotations`, `apply_pairwise` / `apply_product`,
`deep_lift_shap`, `pisa` — all present in the tree, all with built-in batching and device
handling. Multi-input/multi-output support is per-operation rather than universal:
`deep_lift_shap` requires a single `(batch, n_targets)` tensor (`deep_lift_shap.py` L255–259)
and the default ISM aggregation raises unless `raw_outputs=True`
(`saturation_mutagenesis.py` L198–201). Paper p.3: *"These operations have built-in
batching, support for alternative data types and devices, and work out-of-the-box on
multi-input/output models."* Fig 1A "Screening" column.

The atomic `ersatz.insert` / `delete` operations return, respectively, longer / shorter edited
sequence tensors at one shared coordinate across the batch; they do not search indel positions.

*Source:* paper abstract and p.3 (verified verbatim in the PDF); Fig 1A; module tree
(`ablate.py`, `marginalize.py`, `space.py`, `saturation_mutagenesis.py`, `variant_effect.py`,
`product.py`, `deep_lift_shap.py`, `pisa.py`); `ersatz.py` L20–94, L300–340.

---

### BLOCK C — genomics integration

---

**`genome_coordinates` → yes**

`io.extract_loci(loci, sequences, signals=None, in_signals=None, chroms=None, in_window=2114,
out_window=1000, max_jitter=0, exclusion_lists=None, min_counts=None, max_counts=None,
summits=False, return_mask=False, ...)` (signature at `io.py` L246–265; the quoted `loci`
parameter doc at L310–315) — *"Either the path to a bed file
or a pandas DataFrame object containing three columns: the chromosome, the start, and the end"*;
`sequences` is a FASTA path read via **pyfaidx**; `signals` are bigwigs read via **pybigtools**.
Supports exclusion-list filtering, chromosome subsetting, summit-centring, jitter, and count
thresholds.

**[ADDED — deserves its own line, per review:]** `utils.example_to_fasta_coords(example_df,
loci_df, window=None, one_indexed=False)` maps window-relative spans (seqlets, hits) **back**
to genome coordinates by cross-referencing the originating BED, handling zero/one-indexing:
*"Frequently, though, one wants to know these coordinates on the original genome or other
entity encoded in a FASTA file. This function cross-references the example indexes and
BED-formatted locus file that the examples were extracted from…"* (`utils.py` L233–253).
This is a genuine coordinate round-trip and the closest tangermeme analogue to emitting
genome-anchored records.

Also `match.extract_matching_loci` (BED + FASTA → GC-matched negative DataFrame).

*Source:* `io.py` L246–265; `utils.py` L233–253; `match.py`; Tutorial C1.

---

**`transcript_models` → no**

**[CITATION CORRECTED]** Recursive grep over **all 28 modules + 15 skill docs**: zero hits for
`gtf`, `gff`, `transcript`, `exon`, `intron`. `io.py`'s complete public surface is
`extract_loci`, `one_hot_to_fasta`, `read_meme`, `read_vcf`. Genomic input is BED intervals +
FASTA + bigwig only; there is no gene model, no transcript structure, no annotation database.

*Source:* own recursive grep, this pass; `io.py` public functions.

---

**`exon_intron_split_codons` → no**

**[CITATION CORRECTED]** Same recursive grep over all 28 modules and 15 skill docs: zero hits
for `exon`, `intron`, `codon`, `reading frame`, `orf`, `amino`, `protein`. Sequences are
alphabet-agnostic one-hot tensors (`alphabet: list[str] = ['A','C','G','T']` as a plain
default argument); there is no coding-biology layer at all.

*Source:* own recursive grep, this pass.

---

**`hgvs_input` → no**

Zero hits for `hgvs` anywhere in the package or skill docs. The actual variant interface is a
COO integer tensor — `variant_effect.py` L63–72: *"A set of variants that should be substituted
into each sequence. This tensor is formatted like a COO-sparse matrix where each row is a
single variant, the first column is the index in `X`, the second index is the position in that
example, and the third index is the index into the alphabet…"* — and
`_skills/data/references/variant-effect.md` warns explicitly that **position is relative to
the extracted window `X`, not a genomic coordinate**. README example:
`substitutions = torch.tensor([[0, 1058, 0]])`. `io.read_vcf` loads a VCF into a DataFrame,
but the user must convert to indices themselves. **No nomenclature parser of any kind.**

*Source:* own recursive grep; `variant_effect.py` L63–72;
`_skills/data/references/variant-effect.md`; `io.py` `read_vcf`; README.

---

**`vcf_vep_output` → no**

`read_vcf` is **input-only**: *"returns a pandas DataFrame with the comments filtered out …
only the first 9 columns"* (drops genotype columns; no BCF support). The writer search over the
whole tree (this pass) finds sequence/data-serializing code in **`io.py` alone**, and the only
writer is `one_hot_to_fasta` (FASTA). **There is no VCF writer and no VEP-consumable
emitter.**

*Source:* `io.py` `read_vcf` docstring; full-tree writer grep, this pass.

---

**`consequence_annotation` → no**

`tangermeme.annotate` exists, but "annotation" here means **motif matching of seqlets against
a motif database via tomtom-lite**, plus co-occurrence and spacing statistics. Its four
functions are `annotate_seqlets`, `count_annotations`, `pairwise_annotations`,
`pairwise_annotations_spacing`; the README's general definition is *"an annotation is any
genomic span"*. Paper Methods pp.7–8 describe exactly that. There is **no molecular-consequence
logic** (stop-gained / synonymous / frameshift / in-frame) anywhere — consistent with the zero
hits for `codon`, `amino`, `protein`, `orf`.

*Source:* `annotate.py`; README; paper pp.7–8; own recursive grep.

---

### BLOCK D — physical construction

---

**`primer_design` → no**

**[CITATION CORRECTED]** Zero hits for `primer`, `oligo`, `melting`, `Tm` across all 28
modules and 15 skill docs; the single `adapter` hit is `pisa.py` L174, an `nn.Module` output
adapter. No wet-lab layer appears in the README, docs site, tutorials, release history
(`whats_new.rst` through v1.4.0), or the paper.

*Source:* own recursive grep, this pass; `pisa.py` L174; `docs/whats_new.rst`.

---

**`codon_optimization` → no**

Zero hits for `codon`, `amino`, `protein`, or any translation concept. The `alphabet` argument
is generic (default `['A','C','G','T']`) with no codon table, translation layer, or
expression-host concept.

**[CAVEAT ADDED, per review]** The extraction's separate claim that alphabet-agnosticism means
"the same operations apply to protein" must be phrased as **the README's claim, not as
demonstrated protein support** — the README says *"although the library was built with
operations on DNA sequences in mind, all functions are extensible to any alphabet"*, and the
word `protein` appears **nowhere** in the source.

*Source:* own recursive grep, this pass; README; function signatures across `ersatz.py`,
`utils.py`, `io.py`.

---

**`synthesis_constraints` → no**

Zero hits for `homopolymer`, `restriction`, `synthesis`, `secondary structure`, or any
repeat/oligo-length constraint check. The GC machinery that exists serves **background
selection**, not synthesizability: `utils.gc_content`, and `match.extract_matching_loci`
whose stated rationale (`_skills/data/references/motif-effects.md`) is the GC "mirage" effect
— *"Uniformly random sequences are much higher in GC content than real genomic DNA …
substituting a motif into the wrong background can produce a 'mirage' effect"*.

**[ADDED, per review]** `extract_matching_loci`'s other filters confirm the same reading:
`max_n_perc` (applied to **both** the input loci and the candidates) and an optional
`bigwig` + `signal_beta` / `signal_threshold` to reject candidate negatives carrying too much
signal (`match.py` L285–352, L377–395). Still background/negative selection, still not
synthesis QC.

*Source:* own recursive grep, this pass; `utils.gc_content`; `match.py` L285–395;
`_skills/data/references/motif-effects.md`.

---

### BLOCK E — engineering

---

**`interface` → yes** — **RENDER THIS CELL AS: "Python API; one installation-only CLI (no analysis CLI, GUI, or web service)"**

**[ENCODING NOTE, per review — mandatory.]** A bare "yes" is uninformative next to
MutationMaker's web GUI + REST API and VaLiAnT's CLI-only row; the table cell must print the
discriminating fact.

Verified: `pyproject.toml` `[project.scripts]` contains **exactly one** entry —
`tangermeme-install-skills = "tangermeme._skills.install:main"` — which only installs the
bundled Claude Code Agent Skill. There is **no `[project.entry-points]` group**, no web
service, no GUI. Analysis CLIs were removed upstream: README — *"These FIMO and Tomtom
command-line tools have been moved to memesuite-lite … Please use those!"*
Notable extra surface: a bundled **Claude Code Agent Skill** (`_skills/data/SKILL.md` + 14
reference docs) so an LLM agent can drive the API. Optional extras are `[dev]`, `[docs]`,
`[interactive]` (mpld3) only.

*Source:* `pyproject.toml` (fetched raw, this pass); README; `tangermeme/_skills/` tree.

---

**`license` → yes** — **RENDER THIS CELL AS: "MIT"**

**[ENCODING NOTE, per review — mandatory.]** Permissive-vs-copyleft is the comparative point
(cf. VaLiAnT's AGPL-3.0), so the cell must print the licence name.

Verified three independent ways this pass: `pyproject.toml` → `license = "MIT"` with
`license-files = ["LICENSE"]`; PyPI JSON → `license_expression: MIT`; GitHub API →
`license.spdx_id: "MIT"`. Plus paper Code Availability (p.9, verbatim): *"tangermeme is free
and open source software available under the MIT license at
https://github.com/jmschrei/tangermeme."*

*Source:* `pyproject.toml`; PyPI JSON; `gh api repos/jmschrei/tangermeme`; paper p.9.

---

**`maintained` → yes (actively)**

Live figures, re-fetched this pass:

| Fact | Value | Source |
|---|---|---|
| Latest release | **1.4.0, uploaded 2026-06-25T18:05:55** | PyPI JSON |
| Repo `pushed_at` | **2026-07-19T18:34:08Z** | GitHub API |
| Repo `updated_at` | 2026-08-07T19:42:00Z | GitHub API |
| `archived` | **false** | GitHub API |
| Stars / forks / open issues+PRs | **308 / 32 / 8** (`open_issues_count`; enumerating gives 5 issues + 3 PRs) | GitHub API |
| Development status classifier | "Development Status :: 5 - Production/Stable" | `pyproject.toml` |
| Python / torch floor | `>=3.10` / `torch>=2.0` | `pyproject.toml` |

Release ladder over the last 12 months (PyPI `upload_time`, exact): 1.0.0 (2025-08-27),
1.0.2 (2026-01-19), 1.0.3 (2026-02-15), 1.0.4 (2026-04-23), 1.1.0 and 1.2.0 (both 2026-05-27),
1.3.0 (2026-06-23), 1.4.0 (2026-06-25) — **8 releases in 12 months**. CI unit-test badge; Read
the Docs live and matching the repo tree.

*Source:* PyPI JSON and GitHub API, both fetched live this pass; `pyproject.toml`;
https://tangermeme.readthedocs.io/en/latest/.

---

## 4. tangermeme's own documented examples (candidates for PoolParty reproduction)

Verified against the `docs/` tree this pass — the list below matches the repository exactly.

| Notebook | Content |
|---|---|
| Tutorial A1 — Ersatz Sequence Manipulation | `substitute` / `multisubstitute` / `insert` / `delete` / `randomize` / `shuffle` / `dinucleotide_shuffle` / `reverse_complement`, incl. inserting JASPAR AP-1 (MA1144.1) into a random background |
| Tutorial A2 — Predictions | batched `predict` |
| Tutorial A3 — DeepLIFT/SHAP | attributions |
| Tutorial A4 — Seqlets | recursive seqlet caller + TF-MoDISco caller |
| Tutorial A5 — Annotations | seqlet annotation, counting, pairwise, spacing |
| Tutorial B1 — Marginalization | motif effect on predictions/attributions in backgrounds |
| Tutorial B2 — Ablation | shuffle-out of a window |
| **Tutorial B3 — Spacing** | **motif-pair spacing sweeps (AP-1 cooperativity fading with distance) — a small combinatorial motif-placement library** |
| Tutorial B4 — Saturation Mutagenesis | in silico ISM |
| Tutorial B5 — Variant Effect | substitution / deletion / insertion effects |
| **Tutorial B6 — Design** | **`greedy_substitution` and `beam_substitution` against Beluga to design a sequence with high AP-1 binding, from one random sequence + a JASPAR motif subset** |
| **Tutorial B7 — Cartesian Product** | **`apply_pairwise` / `apply_product`, incl. nesting `marginalize(func=deep_lift_shap)` inside a product** |
| Tutorial C1 — IO and Data Loading | `extract_loci`, BED/FASTA/bigwig |
| Tutorial C2 — Plotting | logos (incl. `plot_logo(annotations=...)`), attributions |

Vignettes (`docs/vignettes/`): **Inspecting What Cis-Regulatory Features a Model Has Learned**
(README: *"If you only read one vignette, read THIS ONE"*); Attribution Trickiness and
DeepLiftShap Implementations; Wrappers are Productivity Hacks.
How-to: *How To — Reduce Friction and Save Time with Tangermeme*.
Paper figure notebooks: `docs/paper/Fig1_Timings_Examples.ipynb`,
`docs/paper/Fig2_Seqlet_Calling_and_Downstream_Analyses.ipynb` (PLD6 promoter;
BPNet / Beluga / ChromBPNet / ProCapNet).

**Reproduction targets most relevant to PoolParty:** Tutorial B3 (spacing sweep = a small
combinatorial motif-placement library) and Tutorial A1 (motif implantation into
shuffled/GC-matched backgrounds). Tutorial B6 (design) is **search**, not enumeration —
PoolParty would not replace it.

---

## 5. Capabilities not covered by the row list

Original list, corrected and extended with everything the review found missing (new or
corrected items marked **[NEW]** / **[CORRECTED]**).

- DeepLIFT/SHAP with built-in batching and correctness guarantees on convergence deltas
  (paper Fig 1B: *"these implementations can silently fail, causing the convergence deltas to
  become quite large instead of near zero … whereas tangermeme will never do so"*).
- PISA (`pisa.py`) as a saliency-style attribution alternative.
- **Novel recursive, variable-length seqlet caller** — *"This seqlet caller is the first to
  call variable-length seqlets directly"* — plus a TF-MoDISco-style caller.
- Seqlet annotation against a motif database via tomtom-lite (`memelite`), then
  `count_annotations`, `pairwise_annotations`, `pairwise_annotations_spacing` (motif
  co-occurrence and spacing statistics).
- **[NEW] `marginalize_annotations` (`marginalize.py` L133) / `ablate_annotations`
  (`ablate.py` L171)** — perturb each of a *set* of annotation-defined regions individually,
  returning `PerturbationAnnotationsResult(y_befores, y_afters)` stacked per annotation.
  `marginalize_annotations` extracts `X[idx, :, start:end]` as the motif and marginalizes each
  into a background set `X0` (an n_motifs × n_backgrounds sweep in one call). This is the
  closest thing tangermeme has to a per-element indexed sweep.
- **[NEW] `plot_logo(annotations=...)`** — the *static* logo accepts a
  `motif_name/start/end/strand/score` DataFrame and draws labelled non-overlapping boxes under
  the glyphs (`plot.py` L142–233, with `check_box_overlap` / `place_new_box` at L35–115).
  `interactive_logo` adds mpld3 tooltips but needs the optional `[interactive]` extra.
- **[NEW] Replicate-generating control primitives** — `ersatz.randomize(X, start, end, probs, n=)`,
  `ersatz.shuffle(X, start, end, n=)`, `ersatz.dinucleotide_shuffle(X, start, end, n=20)`
  each return an extra sequence axis (`shape=(-1, n, len(alphabet), length)`), i.e. bulk
  materialized scramble/random controls — one of the few places the package returns sequences
  rather than model outputs.
- **[NEW] `variant_effect.substitution_effect`'s heterogeneous per-example COO edit table** —
  rows grouped by `example_idx` are applied together, and multi-row adjacent positions encode
  *"entire motifs or just multiple characters"*; the nearest thing in the package to a
  heterogeneous edit specification.
- **[NEW] Exhaustive internal enumeration in `design/`** — `design/_substitute.py`'s numba
  kernel (*"This function takes a motif and inserts it at all possibilities"*) plus
  `greedy_substitution` L188–205 materialize {motifs} × {all allowed positions} × {fwd, rc}
  each round; only the argmin escapes.
- **[NEW] `beam_substitution(n_best=k)` / `screen(n_best=k)`** return
  `shape=(n_best, len(alphabet), length)` — a small ranked *set* of designed sequences. Note
  the two are different modes: `beam_substitution` edits a template, whereas `screen` generates
  independent random candidate batches de novo (default generator `utils.random_one_hot`,
  overridable via `func=` + `additional_func_kwargs`) and keeps the `n_best` lowest-loss
  candidates in a heap without de-duplicating them (`design/screen.py` L19–44, L155–185).
- **[NEW] `SaturationMutagenesisRawResult(y0, y_hat)`** via `raw_outputs=True` — the only route
  to per-mutant *values* rather than aggregated attributions, and notable precisely because it
  still does not hand back the mutant *sequences*.
- **[CORRECTED] k-mer utilities (`kmers.py`)** — `kmers(X, k, scores=None)` returns a **dense
  `torch.Tensor`** count matrix (`kmers.py` L68–70; the `-> scipy.sparse.csr_matrix` annotation
  at L22 is stale) **optionally weighted by the summed per-position attribution across each
  k-mer's span** (L40–44, L61–66); `gapped_kmers(X, scores=None, ...)` does return a
  `scipy.sparse.csr_matrix` (L250–253), but weights each gapped k-mer by the **average**
  attribution over its characters (`new_gkmer_attr / k`, L136), not the sum its docstring
  claims. The attribution-weighted mode is the non-obvious capability the original list omitted.
- **[NEW] `utils.example_to_fasta_coords`** — window-relative spans → genome coordinates by
  cross-referencing the originating BED, with zero/one-indexing handling (`utils.py` L233ff).
- **[CORRECTED] GC-matched genome-wide negative/background sampling (`match.py`)** — beyond GC
  matching, `extract_matching_loci` optionally takes a `bigwig` + `signal_beta` to reject
  candidate negatives with too much signal, and applies `max_n_perc` to both inputs and
  candidates (`match.py` L285–395).
- Extremely fast one-hot encoding (paper Fig 1C: 1.61 s for hg38 chr1, ~3× the next fastest —
  benchmarked on tangermeme **v0.5.1**, paper Methods p.6, not the v1.4.0 surveyed here).
- `utils.entropy` / `information_content` / `pwm_consensus`.
- Alphabet-agnostic throughout — **README's claim**: *"although the library was built with
  operations on DNA sequences in mind, all functions are extensible to any alphabet."*
  (Note: no protein support is demonstrated anywhere; `protein` appears nowhere in the source;
  and support is not literally uniform — `utils.reverse_complement`'s complement map and
  `utils.gc_content` are DNA-semantic by default, and `greedy_substitution` /
  `beam_substitution` default to `reverse_complement=True`.)
- Model-oracle sequence design: `screen`, `greedy_substitution`, `beam_substitution`,
  `greedy_marginalize` (the last returns a variable-width one-hot *construct* whose internal
  `N`s mark positions no chosen motif covers — it initialises an all-`N` background, writes each
  chosen motif into its slice in order, and strips terminal `N`s, `greedy_marginalize.py`
  L242–246; the bundled `design.md` L96–97 instead says `N` marks "where overlapping motifs
  cancel", which the implementation does not do — later motifs simply overwrite earlier ones).
- Bundled Claude Code Agent Skill (`SKILL.md` + 14 reference docs) for LLM-driven use, with
  the `tangermeme-install-skills` console script.
- Reads bigwig signal tracks, MEME motif files, and VCF (input only).

**Suggested new comparison rows this tool would motivate (if the referees want them):**
"model-oracle sequence optimization / design search"; "attribution & interpretation of a
trained model"; "background/negative-set construction (GC-, signal- or dinucleotide-matched)".

---

## 6. Stated limitations

- Scope is deliberately "everything-but-the-model": Fig 1A marks **Model Training, Data
  Preprocessing and Model Zoo as NOT provided** (unlike Selene, EUGENe, gReLU, CREsted). A
  trained PyTorch model is a prerequisite for most of the package.
- Design is discrete and combinatorial; the bundled docs call it greedy, though
  `design/__init__.py` also exports the non-greedy random `screen`.
  `_skills/data/references/design.md`: *"`tangermeme.design` is
  discrete and greedy. When you want gradient-based, minimal edits to a specific template
  sequence … reach for the `ledidi` library."*
- `greedy_substitution` requires `X` of batch size 1 (`shape=(1, len(alphabet), length)`) and a
  positive `max_iter` (`max_iter=-1` documented as *"a silent no-op"*).
- The recursive seqlet caller assumes constant sequence length: *"In principle, l does not need
  to be constant across examples but in our current implementation it is assumed to be fixed."*
- FIMO/TOMTOM CLIs no longer ship with tangermeme (moved to `memesuite-lite`).
- `read_vcf` drops genotype columns past column 9 and does not support BCF.
- `plot.interactive_logo` requires the optional `[interactive]` extra (mpld3); the function
  ships in a default `pip install tangermeme` but **raises `ImportError` when called** without
  the extra.
- Attribution-based workflows are documented as expensive (marginalize + `deep_lift_shap`
  ≈ 200 forward/backward passes for 5 examples).
- `dinucleotide_shuffle`'s `random_state` is applied as `random_state + i` per batch position,
  so reproducibility depends on batch composition (documented in the docstring).

---

## 7. Definitional risks and disagreements — explicit flags

**No unresolved extractor/reviewer disagreements exist.** The adversarial review returned
"supported" on all 24 verdicts; every correction it raised was to evidence, citation, or
phrasing, and all have been applied above. No cell is set to `unknown`.

Three cells are **generous by design** and should carry table footnotes so no referee can
claim either understatement or overstatement:

1. **`library_as_object = partial` (generous).** A strict reading gives "no" — the package
   defines no class representing a library, and no class carrying identity or metadata for a
   sequence (six classes total, five NamedTuples plus a Warning; the one that does hold
   sequences, `AttributionReferencesResult`, holds DeepLIFT backgrounds). "partial" credits
   only the batch-first ergonomics. Erring generous
   is the right direction for a referee response.
2. **`assay_mpra = partial` (generous).** Under a strict *"designs an MPRA library"* reading
   this is "no". Keep "partial" — it credits the real model-oracle design layer — but state
   the row definition explicitly in the table caption.
3. **`assay_dms = no` (definitional).** Defensible **only** under the row definition *"designs
   a DMS library for synthesis"*. The package contains a module named `saturation_mutagenesis`
   and a documented ISM-style design mode, so the footnote must state the three saving facts:
   in silico only; no coding-biology layer (zero codon/AA/ORF concepts anywhere, though the
   ISM primitive itself is alphabet-generic); and the mutant sequences are never returned.

Two cells must **not** be rendered as bare yes/no:
`interface` → **"Python API; one installation-only CLI (no analysis CLI, GUI, or web service)"**;
`license` → **"MIT"**.

---

## 8. Bottom line for the PoolParty comparison

tangermeme is the Tier-1 tool with the most genuine overlap in *sequence-manipulation
primitives*: `ersatz.substitute` / `multisubstitute` really do implant motifs;
`space` really does sweep a 2-D spacing grid; `randomize` / `shuffle` /
`dinucleotide_shuffle` really do bulk-generate scramble controls; `design.greedy_substitution`
really does enumerate motifs × positions × strands internally; and `func=` / `apply_product`
really do compose operations.

The overlap ends at the abstraction boundary. **tangermeme's composition target is a model
call, not a library.** Every model-operation pipeline terminates in predictions or attributions;
and the sequences it does hand back — the products of the `ersatz` edits (`insert`,
`substitute`, `multisubstitute`, `delete`), DeepLIFT backgrounds, generated controls, and
design winners — are all anonymous tensors with no identity, no provenance, no names, and no
way to hold heterogeneous design intents in one specification (the package defines no such
class). It has
no wet-lab layer (no barcodes, primers, codons, synthesis checks) and no
transcript/variant-nomenclature layer (no GTF, HGVS, VCF output, consequence calls).

The honest framing for the referee response is **complementary, not competing**: PoolParty
specifies and materializes the library; tangermeme is one of its most natural consumers. This
framing survived adversarial review unfalsified.
