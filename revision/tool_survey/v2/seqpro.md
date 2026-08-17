# SeqPro — v2 capability record

**Slug:** `seqpro`
**Full name:** SeqPro — "Genomic sequence preprocessing toolkit"
**Citation key:** `Klie2023kg`
**Tier:** 1
**Authors:** David Laub, Adam Klie — ML4GLand org (UCSD)
**Version assessed:** 0.22.0 (`main`, commit 2026-07-27); metadata re-verified 2026-08-10
**Row set:** v2 (`ROWS_v2.md`)
**Basis:** `final/seqpro.md` (adversarially reviewed; all 24 v1 values returned "supported"), plus two
targeted re-checks of the live repo for the new rows (`rag/__init__.py` `__all__`;
`alphabets/_alphabets.py` IUPAC grep), both performed 2026-08-10.

---

## Scope statement (read before the table)

SeqPro is a general-purpose ML **preprocessing** library for DNA/RNA/protein sequences, and its authors
say so: *"SeqPro is a Python package for processing DNA/RNA sequences… fully functional on its own but
is heavily utilized by other packages including SeqData, MotifData, SeqExplainer, and EUGENe."* It is
not in the library-design problem space, and this record scores **only** functionality that would
participate in specifying a DNA library. Encoding/tokenization, file-format I/O, coverage statistics,
count normalization, and QC plotting are named where a row asks about them and otherwise excluded.
A column of "no"s here is a scope finding, not a quality judgment: SeqPro is one of the most actively
maintained and most rigorously engineered packages in this survey (Rust/PyO3 + maturin, Numba, CodSpeed
CI benchmarking, pandera/narwhals schema validation, `py.typed`, hypothesis property tests, 55 PyPI
releases, 10 releases in the six weeks before assessment). Two of its authors are plausible referees.

---

## Block A — library specification

| Key | Value |
|---|---|
| `library_first_class_object` | **yes** |
| `composable_operations` | **partial** |
| `lazy_generation` | **no** |
| `library_algebra` | **partial** |
| `exhaustive_single_scans` | **no** |
| `sampled_random_mutagenesis` | **no** |
| `higher_order_combinatorial` | **no** |
| `heterogeneous_components_one_library` | **no** |
| `combinatorial_motif_place` | **no** |
| `barcode_generation` | **no** |
| `per_sequence_provenance` | **no** |
| `automatic_naming` | **no** |
| `design_visualization` | **no** |

**`library_first_class_object` — yes.** *(v1 `library_as_object` = partial; raised on the split.)*
SeqPro never writes a library file — everything it produces is an in-memory object the user holds,
inspects, transforms and passes onward. Fixed-length sets are `(N, L)` NumPy `S1`/`uint8` arrays.
Variable-length sets are `sp.rag.Ragged` (`rag/_core.py`:
`class Ragged(NDArrayOperatorsMixin, Generic[RDTYPE_co])`), a Rust-native Arrow-style flat buffer +
offsets container with exactly one ragged axis, NumPy ufunc dispatch, `from_lengths` / `from_fields` /
`.fields`, `to_padded` / `to_packed` / `to_ak()` / `to_chars()` / `to_strings()`, and `np.memmap`
safety. `sp.rag.zip({"seq": …, "score": …})` builds a **record-layout** Ragged carrying multiple
aligned per-sequence fields. The design philosophy is explicitly batch-first (`SKILL.md`: *"No Python
loops over sequences in library code."*).
**Caveat that the split does not absorb, and the reason this is flagged:** the object is a materialized
*sequence batch*, not a *design specification* — it carries no design steps, no build history and no
names. If the table legend reads this row as "the designed library is an object", SeqPro should be
downgraded to partial; if it reads as the ROWS_v2 wording ("an object the user can hold, inspect,
transform and pass onward — as opposed to a tool that only writes files"), SeqPro clears it outright.
*Source:* `python/seqpro/rag/_core.py`, `rag/__init__.py`, `docs/ragged.md`, `skills/seqpro/SKILL.md`.

**`composable_operations` — partial.** *(rename of v1 `dag_chaining` = partial; value carried
unchanged.)* `transforms/augmentation.py` defines `Sequential(*transforms)` whose `__call__` body is
verbatim `for t in self.transforms: x = t(x)` / `return x`, and `Random(p, *transforms, seed=…)`,
which applies a nested chain with probability `p` — so composition and arbitrary nesting genuinely
exist and are not a fixed pipeline. Declared surface:
`transforms/__all__ = ["TMM", "Jitter", "KShuffle", "Random", "ReverseComplement", "Sequential"]`
(`Tokenize` is defined at `augmentation.py:114` but is **not** in `__all__`). This is a
torchvision-style **linear** augmentation chain over already-materialized arrays: no graph, no
branching or merging, no named nodes, and the operands are augmentations rather than design steps.
*Source:* `python/seqpro/transforms/__init__.py`, `python/seqpro/transforms/augmentation.py`.

**`lazy_generation` — no.** *(rename of v1 `lazy_evaluation` = no; value carried unchanged.)*
No sequence-producing function is lazy. `random_seqs`' body is literally
`return seed.choice(alphabet.array, size=shape)`; every encoder and modifier returns a materialized
NumPy array or `Ragged`. Real laziness does exist in the package but never for sequence production:
`sp.gtf.scan(path)` returns a polars `LazyFrame`, and the experimental `seqpro.xr` module runs `ohe`,
`bin_coverage` and `translate` through `xr.apply_ufunc(..., dask="parallelized")` — genuinely lazy,
chunk-parallel, out-of-core evaluation over dask-backed arrays. That is lazy *compute over sequences
that already exist*, not on-demand generation of designed sequences; there is no design step whose
construction could be deferred. (BED reading is eager — `bed.read` dispatches to `_read_bed` /
`_read_narrowpeak` / `_read_broadpeak`, all calling `pl.read_csv`.)
*Source:* `python/seqpro/_modifiers.py:281`, `python/seqpro/xr/__init__.py:42–48, 93–99, 147–152`,
`python/seqpro/gtf.py`, `python/seqpro/bed.py:236, 287, 356`.

**`library_algebra` — partial.** *(new row from the `library_as_object` split.)*
Combining sequence sets does **not** require an external script: `concatenate` is a public member of
`seqpro.rag.__all__` (re-verified 2026-08-10: `__all__ = [OFFSET_TYPE, DTYPE_co, RDTYPE_co, Ragged,
concatenate, hash, is_rag_dtype, lengths_to_offsets, reverse_complement, to_packed, to_padded, zip]`),
`Ragged` implements `__getitem__` (with BREAKING-CHANGE entries in 0.21.x/0.22.0, i.e. a real indexing
API) so an index array subsamples a set, `to_padded`/`pad_seqs` bring ragged sets to a common length so
plain `np.concatenate`/`np.stack` apply, and `rag.zip` merges aligned per-sequence fields.
**Why not "yes":** these are *array* operations, not *library* operations. There is no `sample` or
`repeat`/replicate function (confirmed absent from `rag.__all__`), no weighted or stratified sampling,
and — because SeqPro has no provenance or naming — combining two sets loses all record of which set
each sequence came from. Composition is possible; library-level algebra with semantics is not.
*Source:* `python/seqpro/rag/__init__.py` (`__all__`, re-verified), `rag/_core.py`, `CHANGELOG.md`.

**`exhaustive_single_scans` — no.** *(from the v1 `mixed_mutagenesis_one_pool` split.)*
No mutagenesis of any kind exists in SeqPro. The complete `__all__` (28 names) contains no mutagenesis
symbol; the untruncated recursive git tree (190 blobs, 27 Python files) contains no mutagenesis module;
a grep over all 608 CHANGELOG lines (55 releases) for `mutagen` returns 0 hits. Nothing enumerates
substitutions, deletions or insertions at positions of a template.
*Source:* `python/seqpro/__init__.py`, `python/seqpro/_modifiers.py`, `CHANGELOG.md`, recursive git tree.

**`sampled_random_mutagenesis` — no.** *(from the split.)* There is no mutation-rate or variant-sampling
facility. `sp.random_seqs(shape, alphabet, seed)` draws i.i.d. random sequences **de novo** — it takes
no template and introduces no variants into a parent sequence. `k_shuffle` is a k-let-preserving
permutation of an existing sequence (composition-preserving, introduces no substitutions) and `jitter`
is a random offset crop. None of the three is parameterized by mutation type or mutation rate.
*Source:* `python/seqpro/_modifiers.py:39, 281`, `python/seqpro/__init__.py`.

**`higher_order_combinatorial` — no.** *(from the split.)* Follows from the absence of any mutagenesis
machinery: there are no single mutations, hence no pairwise or higher-order combinations of them in one
sequence, and no combination-order parameter anywhere in the API.
*Source:* `python/seqpro/__init__.py`, recursive git tree, `CHANGELOG.md` grep.

**`heterogeneous_components_one_library` — no.** *(from the split.)* SeqPro has no design components, so
none can be mixed. The nearest capability is real but different in kind: `Ragged` holds sequences of
**different lengths** in one container and `rag.zip` carries several aligned per-sequence fields, so a
structurally heterogeneous *batch* is representable. What is absent is any notion of component *type*
(e.g. exhaustive singles + sampled higher-order + WT replicates) coexisting in one specification, and
any label distinguishing them.
*Source:* `python/seqpro/rag/_core.py`, `rag/__init__.py`, `python/seqpro/__init__.py`.

**`combinatorial_motif_place` — no.** No motif-insertion or implantation function exists. A source-wide
grep for `motif` across every module returns 0 hits; CHANGELOG returns 0 hits; `docs/api/` contains only
`index, alphabets, bed, gtf, ragged, types`. Motif handling in the ML4GLand suite lives in the separate
packages `MotifData` (org repo, last pushed 2023-07-21) and `SeqExplainer`.
*Source:* source-wide grep, `docs/api/` listing, https://github.com/ML4GLand org listing.

**`barcode_generation` — no.** `sp.random_seqs` draws i.i.d. nucleotides with **no constraints** (no
distance, GC, homopolymer or uniqueness filter); `sp.gc_content` only measures and cannot filter or
enforce; nothing attaches a barcode to a sequence. An `edit_distance(seq1, seq2, dual=False)` exists
only in `python/seqpro/experimental/_experimental.py`, which is (a) headed *"Unmaintained experimental
module, not part of the public package surface."*, (b) a **Hamming** distance requiring equal-length
inputs (`assert len(seq1) == len(seq2)`), not an edit distance, and (c) not wired into any pool
construction. `python/seqpro/experimental/` has **no `__init__.py` at all** (raw URL 404; absent from
the recursive tree), so it is not an importable package.
*Source:* `python/seqpro/_modifiers.py:281`, `_analyzers.py:32`, `experimental/_experimental.py`,
recursive git tree, raw-URL 404.

**`per_sequence_provenance` — no.** Nothing records how a sequence was built, and there is no design
step whose provenance could exist. The record-layout `Ragged` (`rag.zip` / `Ragged.from_fields`) is a
generic mechanism letting a *user* carry arbitrary parallel per-sequence fields, and
`rag.hash("sha256"|"md5"|"rapidhash")` can fingerprint sequences — but SeqPro itself writes no
provenance and no operation in the package emits build metadata.
*Source:* `python/seqpro/rag/__init__.py`, `python/seqpro/rag/_ops.py:410`.

**`automatic_naming` — no.** No naming or ID-generation function anywhere; outputs are anonymous array
rows. Two robustness notes so the cell cannot be nitpicked: `bed.py`'s `BEDSchema` *reads* a `name`
column from a BED file, and a record-layout `Ragged` exposes `.fields` (field names, not per-sequence
names) — neither generates per-sequence identifiers.
*Source:* `python/seqpro/__init__.py`, `docs/api/index.md`, `python/seqpro/bed.py`, `rag/_core.py:303`.

**`design_visualization` — no.** No plotting in the public API and no plotting dependency.
`experimental/_visualizers.py` defines `plot_gc_content` and `plot_nucleotide_content`, but its own
header reads *"Unmaintained experimental module, not part of the public package surface."* and
*"NOTE: this module is currently broken at import time — `gc_content_seqs` and
`nucleotide_content_seqs` were renamed…"*; its directory has no `__init__.py`, so it cannot be imported
at all. Even if working, these are QC histograms of a sequence set, not a view of a design.
*Source:* `python/seqpro/experimental/_visualizers.py` (verbatim header), `_analyzers.py`, recursive tree.

---

## Block B — assay coverage

| Key | Value |
|---|---|
| `assay_dms` | **no** |
| `assay_mpra` | **no** |
| `assay_insilico` | **partial** |

**`assay_dms` — no.** No coding-variant or DMS machinery: no variant table, no codon-substitution
enumeration, no mutant-library construction. *Pre-rebuttal (a referee may say "SeqPro handles stop
codons"):* `AA.translate(..., truncate_stop=True)` truncates each output at the first stop codon
inclusive (Ragged input only, after drop-compaction), backed by a dedicated Rust kernel
`_translate_stop_ends` and an OHE path locating the `*` one-hot column. That is genuine
coding-sequence awareness — and still not DMS design: truncating a translated sequence enumerates
nothing, mutates nothing, and produces no variant library.
*Source:* `python/seqpro/alphabets/_alphabets.py:421–720` (docstring 485–487, kernel at 648, 714),
`src/translate.rs`.

**`assay_mpra` — no.** No regulatory-library design of any kind: no promoter/enhancer scaffolding, no
barcoding, no adapter or cloning handles, no motif perturbation, no oligo assembly. Confirmed against
`SKILL.md`'s "When to use" list, `__all__`, and the full source tree.
*Source:* `python/seqpro/__init__.py`, `skills/seqpro/SKILL.md`, recursive git tree.

**`assay_insilico` — partial.** SeqPro supplies the standard *primitives* used when probing
sequence-to-function models: one-hot encode/decode, integer tokenization, padding, `k_shuffle`
(k-let-preserving shuffle — the canonical null/background for in-silico experiments, dispatched to the
Rust kernel `src/kshuffle.rs` and CodSpeed-benchmarked), `random_seqs` (random backgrounds), `jitter`
and `ReverseComplement` augmentations, and `Ragged` batching. Explicitly framed as ML infrastructure:
*"heavily utilized by … SeqExplainer and EUGENe."* It does **not** design in-silico libraries, run
models, do ISM, or score anything — that layer is SeqExplainer/EUGENe, separate packages.
Primitives-only, therefore partial and never yes.
*Source:* `python/seqpro/_modifiers.py:39`, `src/kshuffle.rs`, `README.md`.

---

## Block C — genomics integration

| Key | Value |
|---|---|
| `genome_coordinates` | **partial** |
| `transcript_models` | **partial** |
| `exon_intron_split_codons` | **no** |
| `vcf_vep_output` | **no** |
| `consequence_annotation` | **no** |

**`genome_coordinates` — partial (legend-dependent — see flags).** `python/seqpro/bed.py` provides
`read` (dispatching to `_read_bed`, `_read_narrowpeak`, `_read_broadpeak`), `sort`, `with_len`, and
`to_pyr`/`from_pyr` (pyranges interop); `python/seqpro/_coords.py` provides `CoordSchema`,
`detect_schema`, `set_schema` with **5 built-in schemas** (bed/pb/pr/gtf/gff) covering 0- vs 1-based and
open vs closed conventions; `polars-config-meta` tags frames with `coordinate_system_zero_based`;
`with_len` is **narrowPeak-summit aware**, recentering on `2 * (chromStart + peak)` (falling back to
the interval midpoint when `peak` is null/absent). So SeqPro ingests and manipulates genomic intervals
competently — **but it does not fetch sequence from a reference genome**: a grep across the whole tree
finds no FASTA, 2bit or bigWig reader, and there is no coordinate→sequence resolution (that lives in the
sibling package `SeqData`).
*Source:* `python/seqpro/bed.py:38–164, 219–378`, `python/seqpro/_coords.py`, source-wide grep.

**`transcript_models` — partial.** `sp.gtf.scan(path)` parses GFF/GTF into a polars `LazyFrame` with
columns `seqname, source, feature, start, end, score, strand, frame, attribute`, and `sp.gtf.attr(...)`
returns a polars expression regex-extracting one attribute as Utf8. That is the **entire** GTF module:
`gtf.py` is 1,813 bytes with exactly two public functions and `docs/api/gtf.md` lists only those two.
There is **no transcript object, no exon assembly, no CDS/UTR resolution, no gene-model semantics** — a
user gets a flat annotation table and must build transcript logic themselves. "Partial" is generous, in
the safe direction.
*Source:* `python/seqpro/gtf.py`, `docs/api/gtf.md`.

**`exon_intron_split_codons` — no.** Nothing joins exons or tracks reading frame across exon
boundaries. `pl.col("frame").fill_null(0)` in `gtf.py` is the **only** use of `frame` anywhere in the
package — parsed, defaulted, never consumed. `translate()` operates on a contiguous nucleotide buffer
in stride-3 codons, and `truncate_stop` likewise works within a single buffer, never crossing a
boundary.
*Source:* `python/seqpro/gtf.py`, `alphabets/_alphabets.py`, `src/translate.rs`, grep for `frame`.

**`vcf_vep_output` — no.** No VCF reader/writer, no VEP-compatible export, no variant representation of
any kind. Package I/O is BED / narrowPeak / broadPeak / GTF-GFF in, NumPy / polars / pyranges / awkward
out. Confirmed via source grep and the PyPI `requires_dist` list (`numba, numpy, polars, pyranges,
pandera, pandas, pyarrow, natsort, narwhals, setuptools, awkward, polars-config-meta[polars], attrs`) —
none is a variant-format library.
*Source:* `pyproject.toml`, PyPI JSON, source-wide grep.

**`consequence_annotation` — no.** There is no variant representation anywhere, hence nothing to compare
against a reference allele and nothing to classify. `AA.translate` can emit a stop codon `*` and
`truncate_stop=True` can truncate at the first one, but finding a stop codon in a translated sequence is
not classifying a variant consequence: no reference allele, no variant object, no consequence class.
*Source:* `python/seqpro/alphabets/_alphabets.py`, grep for `variant|consequence`.

*(v1 `hgvs_input` is dropped from the row set. For the record it was "no": 0 hits for `hgvs` in the
recursive tree, every fetched module, the docs, and all 55 CHANGELOG releases.)*

---

## Block D — adjacent / complementary

| Key | Value |
|---|---|
| `primer_design` | **no** |
| `codon_optimization` | **no** |
| `synthesis_constraints` | **no** |
| `degenerate_iupac_codons` | **no** |
| `negative_control_generation` | **yes** |
| `ml_model_in_loop` | **no** |
| `readout_analysis` | **partial** |

**`primer_design` — no.** A source-wide grep for
`primer|oligo|melting|Tm|restriction|overhang|assembl` returns 0 substantive hits; CHANGELOG returns 0.
No primer, oligo, melting-temperature, overhang or assembly functionality.
*Source:* source-wide grep, `CHANGELOG.md`.

**`codon_optimization` — no.** Codon tables exist strictly for **forward** translation:
`AminoAlphabet(codons, amino_acids)` builds `self.codon_to_aa = dict(zip(codons, amino_acids))` and a
64-entry LUT via `_build_translate_lut` (gated by `_can_build_lut`). There is **no reverse map, no
codon-usage frequency table, and no optimizer**. All 6 `codon` hits in the CHANGELOG are
translate-kernel performance work ("O(1) LUT codon→AA lookup — 179× speedup on bench").
*Source:* `python/seqpro/alphabets/_alphabets.py:246, 297–350`, `CHANGELOG.md`.

**`synthesis_constraints` — no.** All analyzers are read-only; there is no homopolymer, repeat,
restriction-site or GC-window check and no enforcement/repair loop. `gc_content`, `nucleotide_content`
and `length` measure only; `_cleaners.py` provides `remove_N_seqs`, `remove_only_N_seqs`,
`remove_whitespace`, `sanitize` — input hygiene, not synthesis feasibility. Three reinforcing facts:
`count_kmers_seq` is **not** public (absent from `__all__`; `sp.count_kmers_seq` raises
`AttributeError`), its batched sibling `_analyzers._count_kmers` is `raise NotImplementedError`
(`# TODO: non-trivial to parallelize/SIMD this`), and `remove_whitespace`'s entire body is `pass`.
The whole `_cleaners` module is private.
*Source:* `python/seqpro/_analyzers.py:138, 165–187`, `_cleaners.py`, `python/seqpro/__init__.py`.

**`degenerate_iupac_codons` — no.** *(new row; re-checked against the live source 2026-08-10.)*
The only module-level nucleotide alphabet is `DNA = NucleotideAlphabet("ACGT", "TGCA")` — canonical
bases only. There is no IUPAC/ambiguity alphabet constant, no expansion of a degenerate sequence into
concrete sequences, and no compression of a sequence set into degenerate codons. The single occurrence
of the string "IUPAC" in the package is a docstring in `_can_build_lut()` noting that *"Non-standard
alphabets (different codon size or extended/IUPAC characters) use the generic `gufunc_translate` scan
instead"* — i.e. a user who **defines their own** extended alphabet will still be able to translate,
which is IUPAC *tolerance* in one code path, not degenerate codon design. `_check_nuc_bytes()`
validates input against allowed characters but never interprets or expands ambiguity codes; `N` is
handled only by the private cleaners `remove_N_seqs` / `remove_only_N_seqs`, i.e. by deletion.
*Source:* `python/seqpro/alphabets/_alphabets.py` (`DNA` constant, `_can_build_lut` docstring,
`_check_nuc_bytes`), `python/seqpro/_cleaners.py` — re-verified 2026-08-10.

**`negative_control_generation` — yes.** *(new row; the strongest cell in Blocks A–D.)*
All three control classes named by the row exist as complete, public, top-level operations, not as
primitives the user must assemble:
- **shuffle** — `sp.k_shuffle(seqs, k)`, k-let-preserving shuffle implemented in Rust
  (`src/kshuffle.rs`, with `src/kshuffle_ref.rs` and `benches/kshuffle.rs`, CodSpeed-benchmarked). This
  is the canonical dinucleotide-shuffle null for regulatory genomics, and the memo's `assay_insilico`
  evidence identifies it as exactly that.
- **reverse / complement** — `sp.reverse_complement` across every representation (`str`, `S1` bytes,
  one-hot), plus `sp.rag.reverse_complement` for Ragged batches (`rag/_ops.py:22`, in `rag/__all__`),
  plus `transforms.ReverseComplement(*types=("dna","track"))` which also reverses coverage tracks.
- **random background** — `sp.random_seqs(shape, alphabet, seed)`, seeded and reproducible.
All three appear in the README/docs Quick Start and are exposed as reusable transform objects
(`KShuffle`, `ReverseComplement`) that compose via `Sequential`/`Random`.
**Honest limit on the "yes":** SeqPro generates control sequences but does not *label* them as controls,
attach them to a designed library, or size a control set relative to a design — it has no library to
attach them to. The capability is the generation, and the generation is first-class and complete.
*Source:* `python/seqpro/_modifiers.py:39, 281`, `python/seqpro/__init__.py` (`__all__`),
`python/seqpro/rag/_ops.py:22`, `python/seqpro/transforms/__init__.py`, `src/kshuffle.rs`, `README.md`.

**`ml_model_in_loop` — no.** *(new row.)* SeqPro runs no model and consumes no model output. There is
no torch/jax/tensorflow/onnx dependency (runtime deps are `numba, numpy, polars, pyranges, pandera,
pandas, pyarrow, natsort, narwhals, setuptools, awkward, polars-config-meta[polars], attrs`), no
scoring function, no ISM, no gradient or optimization loop, and nothing that selects or edits sequences
according to a prediction. The memo is explicit: SeqPro *"does not design in-silico libraries, run
models, do ISM, or score anything — that layer is SeqExplainer/EUGENe"*, which are separate packages.
SeqPro is the *input* side of a model loop (encoding, batching, augmentation), not a participant in one.
*Source:* `pyproject.toml` / PyPI `requires_dist`, `python/seqpro/__init__.py`, `README.md`,
`skills/seqpro/SKILL.md`.

**`readout_analysis` — partial.** *(new row.)* SeqPro reaches genuinely into the assay-count side, but
generically. `python/seqpro/transforms/tmm.py` (205 lines) implements a full **TMM** (edgeR-style
Trimmed Mean of M-values) estimator — a `TMM` class with `fit(counts)` / `transform(counts,
library_size=1e6)`, quantile-based reference-sample selection (`expression_cutoff, log_ratio_trim=0.3,
expression_trim=0.05, apply_weighting, quantile=0.75`), log-ratio and absolute-expression trimming,
inverse-variance weighting, and an `@nb.njit(parallel=True, nogil=True, cache=True)` `_tmm_helper`
kernel — and it **is** exported (`transforms/__all__` includes `"TMM"`). Alongside it:
`bin_coverage(bin_width, normalize)` (public) and `normalize_coverage(method="CPM"|"CPKM",
total_counts)` (private, not in `__all__`).
**Why not "yes":** none of this connects a readout to a designed library. There is no read/FASTQ
ingestion, no demultiplexing, no barcode or variant counting, no assignment of reads to designed
members, and no effect-size or enrichment estimation. TMM normalizes an arbitrary count matrix that the
user must have produced elsewhere.
**Why not "no":** count normalization is a real, exported, well-implemented piece of readout-side
functionality, and omitting it would read to an author-referee as not having read the code.
*Source:* `python/seqpro/transforms/tmm.py`, `python/seqpro/transforms/__init__.py`,
`python/seqpro/_modifiers.py` (`bin_coverage`, `normalize_coverage`), `python/seqpro/__init__.py`.

---

## Block E — engineering and availability

| Key | Value | Real answer |
|---|---|---|
| `interface` | **yes** | Python API only (`import seqpro as sp`) — no CLI, GUI or web service |
| `license` | **yes** | MIT |
| `installable_today` | **yes** | `pip install seqpro` → 0.22.0; prebuilt wheels on PyPI |
| `last_activity` | **yes** | 2026-07-27 (commit + release 0.22.0 + PyPI upload, same day) |
| `documented_examples` | **yes** | 5 shipped example/reference documents; **0** library-design vignettes |

**`interface` — yes (Python API only).** `import seqpro as sp`. `pyproject.toml` declares **no**
`[project.scripts]` and no console entry points — there is no CLI, no GUI, no web service, so there is
none that could be dead. `requires-python = ">=3.10"`; `py.typed` shipped; build backend `maturin`
(Rust extension via PyO3). Docs site (MkDocs/zensical) at https://ml4gland.github.io/SeqPro/, also the
repo's declared GitHub homepage.
*Source:* `pyproject.toml`, GitHub repo metadata.

**`license` — yes (MIT).** GitHub API `license.spdx_id = "MIT"`; `LICENSE` at repo root;
`pyproject.toml` `license = { file = "LICENSE" }`.
*Source:* `api.github.com/repos/ML4GLand/SeqPro`, `pyproject.toml`, `LICENSE`.

**`installable_today` — yes.** `pip install seqpro` installs **0.22.0**, uploaded to PyPI
**2026-07-27T19:07:05**, `requires-python >=3.10`. Wheels are built via maturin for
manylinux/musllinux/windows/macos, so no Rust toolchain is needed (a *source* install would require
one). Runtime deps are all mainstream and pinned loosely: `numba>=0.58.1, numpy>=1.26, polars>=1.21,<2,
pyranges>=0.1.3,<0.2, pandera>=0.31.1, pandas, pyarrow, natsort, narwhals>=2.20, setuptools>=70,
awkward>=2.5, polars-config-meta[polars]>=0.3.2, attrs`. Repo public and **not archived**.
*Source:* PyPI JSON, `pyproject.toml`, GitHub repo metadata (re-verified 2026-08-10).

**`last_activity` — yes: 2026-07-27, the most recent in this survey.** Last commit
`2026-07-27T19:02:07Z` ("bump: version 0.21.2 → 0.22.0"); `pushed_at 2026-07-27T19:07:56Z`; release
0.22.0 tagged `2026-07-27T19:02:25Z`. 55 PyPI releases since `2023-03-26T17:33:00`; **10 GitHub releases
(0.16.0 → 0.22.0) in the ~6 weeks before assessment**. 14 stars, 1 fork, 0 open issues.
*Source:* GitHub repo/commits/releases API, PyPI JSON (re-confirmed 2026-08-10).

**`documented_examples` — yes: 5 shipped, but none is a library-design use case.**
1. **README / docs "Quick Start"** (`README.md`, `docs/index.md`) — `random_seqs`, `pad_seqs`,
   `DNA.ohe`/`decode_ohe`, `tokenize`/`decode_tokens`, `reverse_complement`, `k_shuffle(k=2)`,
   `gc_content`, `nucleotide_content`, `jitter`, `bin_coverage`, `AA.translate`, `AA.ohe`.
2. **"Ragged Arrays" guide** (`docs/ragged.md`) — a worked example over per-position coverage for three
   intervals of different widths: `Ragged.from_lengths`, ufunc arithmetic, `to_padded`, `to_packed`,
   record layouts via `rag.zip`.
3. **API reference** — `docs/api/{index,alphabets,bed,gtf,ragged,types}.md`, mkdocstrings-rendered
   signatures; `docs/api/index.md` lists 14 top-level members.
4. **`skills/seqpro/SKILL.md`** — author-maintained agent-facing cheat sheet with a task→call table plus
   Ragged construction, record-layout and hashing examples; effectively the richest usage document.
5. **Test suite as usage examples** — `SKILL.md` points at `tests/`; ~30 test modules dominated by
   `test_ragged*`, `test_ohe`, `test_tokenize`, `test_translate`, `tests/bed/*`, including
   `hypothesis[numpy]` property-based tests.
Additionally 49 internal design/spec documents under `docs/superpowers/` plus a repo-root `CLAUDE.md`
(intent documents, not user tutorials). **There is no SeqPro tutorial notebook**: `ML4GLand/tutorials`
(last push 2025-02-06) contains notebooks only under `eugene/*` and `seqexplainer/*`, plus
`motifdata/TODO.txt`; `ML4GLand/use_cases` last pushed 2024-01-03.
*Reproducibility note for PoolParty:* none of these is a library-design use case, so there is nothing
to reproduce in PoolParty beyond generating random sequences and computing GC.
*Source:* `README.md`, `docs/index.md`, `docs/ragged.md`, `docs/api/*`, `skills/seqpro/SKILL.md`,
`tests/`, https://github.com/ML4GLand org listing.

---

## v1 → v2 changes

| v1 row | v1 value | v2 row(s) | v2 value(s) | Note |
|---|---|---|---|---|
| `library_as_object` | partial | `library_first_class_object` / `library_algebra` | **yes** / **partial** | Split. The "hold an object" half is unambiguous for SeqPro (`Ragged` + record layout, nothing ever written to file); the "combine libraries" half is real but array-level only (`rag.concatenate`, no sample/repeat, no provenance on merge). |
| `dag_chaining` | partial | `composable_operations` | partial | Rename; value carried unchanged. `Sequential`/`Random` compose and nest arbitrarily but over augmentations, not design steps. |
| `lazy_evaluation` | no | `lazy_generation` | no | Rename; value carried unchanged. Laziness exists (`gtf.scan`, `seqpro.xr` dask) but never for sequence generation. |
| `mixed_mutagenesis_one_pool` | no | `exhaustive_single_scans` / `sampled_random_mutagenesis` / `higher_order_combinatorial` / `heterogeneous_components_one_library` | no / no / no / no | Split resolves to four independent "no"s; SeqPro has no mutagenesis machinery at all, so the split does not separate anything here. |
| `hgvs_input` | no | — | — | Dropped from the row set. |
| `maintained` | yes | `installable_today` / `last_activity` / `documented_examples` | yes / yes / yes | Split into checkable facts: pip-installable 0.22.0; 2026-07-27; 5 examples, 0 design vignettes. |
| — | — | `degenerate_iupac_codons` | no | New; re-checked against live source. |
| — | — | `negative_control_generation` | **yes** | New; `k_shuffle` + `reverse_complement` (incl. Ragged) + `random_seqs`, all public and complete. |
| — | — | `ml_model_in_loop` | no | New; no model dependency or scoring anywhere. |
| — | — | `readout_analysis` | partial | New; exported TMM normalizer + `bin_coverage`/`normalize_coverage`, but no read ingestion, demultiplexing or member-level counting. |

All other v1 values (Blocks B, C, D1–D3, `interface`, `license`) are carried across unchanged.

---

## Flagged cells for human review

1. **`library_first_class_object = yes` — the one value raised in v2, and legend-dependent.** Under the
   ROWS_v2 wording ("an object the user can hold… as opposed to a tool that only writes files") this is
   a clean yes. Under a reading of "library" as *designed* library (with steps, history, names), SeqPro
   should be partial or no. Decide the legend before the table is built; if the legend stays as written,
   note that SeqPro will read as "yes" here while several file-writing design tools read as "no", which
   is correct but may look counterintuitive in the printed table.
2. **`negative_control_generation = yes`.** Defensible and fair — `k_shuffle` is the canonical
   dinucleotide-shuffle null and is a complete, exported, Rust-accelerated operation. The counterargument
   a strict reader could make: SeqPro frames these as *augmentations*, never as controls, and cannot
   attach or label them within a library. If the row is scored as "controls as a first-class feature
   **of a library design**", downgrade to partial. Check this against how tangermeme and ledidi are
   scored so the two are consistent.
3. **`readout_analysis = partial`.** TMM is real, exported and substantial, but it is generic count
   normalization with no link to a designed library. Could reasonably be "no" if the row means
   member-resolved readout analysis (Oligopool Calculator's sense). Consistency with the Oligopool
   Calculator column matters here.
4. **`library_algebra = partial`.** Rests on `rag.concatenate` being public (verified) plus NumPy-level
   stack/index. If the row intends library-level algebra with semantics (labelled provenance on merge,
   sampling/replication as operations), this is "no". There is no `sample` or `repeat`.
5. **`genome_coordinates = partial` — carried over from v1 as the one real framing risk.** Under
   "resolves genomic coordinates to sequence", partial is correct (no FASTA/2bit reader,
   grep-confirmed). Under "accepts genomic intervals", SeqPro is a clean **yes** and the cell would read
   as an understatement an author-referee could attack. **Action: the row legend must state
   coordinate→sequence resolution explicitly.**
6. **`heterogeneous_components_one_library = no`.** SeqPro *can* hold structurally heterogeneous
   (variable-length, multi-field) batches in one `Ragged`. The "no" rests on there being no design
   *component types* to be heterogeneous. Worth a legend sentence so the distinction is visible.
7. **`composable_operations = partial` and `transcript_models = partial` and `assay_insilico = partial`
   — all deliberately generous, in the safe direction** (crediting SeqPro rather than under-crediting).
   A referee could argue "no" for any of them; none can argue "yes".
8. **`degenerate_iupac_codons = no`.** Low risk but flagged: `NucleotideAlphabet` is user-definable and
   `_can_build_lut`'s docstring explicitly contemplates extended/IUPAC alphabets in the translate path.
   That is tolerance, not design/expansion/compression — but an author-referee might cite the docstring.
9. **`documented_examples`** is reported as a count of *shipped example documents* (5). If the survey
   counts *executable* vignettes/notebooks, SeqPro's number is **0** and the other columns must be
   counted the same way.
10. **Method caveat, carried from v1:** source was read from `main` (0.22.0); the package was **not**
    installed or executed. All claims are source-, docs- and metadata-derived. The two v2-specific
    re-checks (`rag/__init__.py` `__all__`; IUPAC in `alphabets/_alphabets.py`) were fetched from
    `raw.githubusercontent.com/ML4GLand/SeqPro/main` on 2026-08-10.
