# SeqPro — FINAL capability record

**Slug:** `seqpro`
**Full name:** SeqPro — "Genomic sequence preprocessing toolkit"
**Citation key:** `Klie2023kg`
**Tier:** 1
**Authors:** David Laub (dlaub@ucsd.edu), Adam Klie (aklie@ucsd.edu) — ML4GLand org (UCSD)
**Paper:** none — software only. SeqPro has no standalone publication; `Klie2023kg` is a legacy
suite-level citation key and is not used as SeqPro evidence here.
**Version assessed:** 0.22.0 (`main`, commit of 2026-07-27)
**Record status:** final. All 24 capability values retain their prior assessment; no value changed.
Five evidence-text defects and five missed capabilities are corrected/incorporated below, and every
correction was independently re-verified for this memo.

---

## Sources consulted

| Kind | URL / ref | Retrieved |
|---|---|---|
| repo | https://github.com/ML4GLand/SeqPro | 2026-08-10 |
| repo API | `api.github.com/repos/ML4GLand/SeqPro` — metadata, license, recursive git tree (`truncated: false`, 190 entries = 162 blobs + 28 trees) | metadata/license 2026-08-14; pinned tree 2026-08-10 |
| repo (raw) | `python/seqpro/__init__.py` (full `__all__`, 28 names) | 2026-08-10 (re-verified) |
| repo (raw) | `python/seqpro/_analyzers.py`, `_modifiers.py`, `_cleaners.py`, `_coords.py`, `_encoders.py` | 2026-08-10 |
| repo (raw) | `python/seqpro/bed.py` (378 lines), `python/seqpro/gtf.py` (1,813 bytes) | 2026-08-10 (re-verified) |
| repo (raw) | `python/seqpro/alphabets/_alphabets.py` | 2026-08-10 (re-verified) |
| repo (raw) | `python/seqpro/transforms/__init__.py`, `transforms/augmentation.py`, `transforms/tmm.py` | 2026-08-10 (re-verified) |
| repo (raw) | `python/seqpro/rag/__init__.py`, `rag/_core.py`, `rag/_ops.py` | 2026-08-10 (re-verified) |
| repo (raw) | `python/seqpro/xr/__init__.py` | 2026-08-10 (re-verified) |
| repo (raw) | `python/seqpro/experimental/_experimental.py`, `_visualizers.py`; **`experimental/__init__.py` → HTTP 404** | 2026-08-10 (re-verified) |
| repo (raw) | `pyproject.toml`, `CHANGELOG.md` (608 lines / 51 release headings), `README.md`, `CLAUDE.md` | 2026-08-10 (re-verified) |
| repo (raw) | `skills/seqpro/SKILL.md` (author-written API cheat-sheet) | 2026-08-10 |
| docs | https://ml4gland.github.io/SeqPro/ ; `docs/index.md`, `docs/ragged.md`, `docs/api/*.md` | 2026-08-10 (re-verified) |
| pypi | https://pypi.org/pypi/seqpro/json | 2026-08-14 (re-verified) |

No installation or clone; evidence comes from the official repository/docs and package metadata.

---

## What SeqPro actually is

From `README.md`, verbatim:

> "SeqPro is a Python package for processing DNA/RNA sequences, with support for protein sequences.
> SeqPro is fully functional on its own but is heavily utilized by other packages including SeqData,
> MotifData, SeqExplainer, and EUGENe."

> "All functions in SeqPro take as input a string, a list of strings, a NumPy array of strings, a NumPy
> array of single character bytes (`S1`), or a NumPy array of one-hot encoded sequences (`uint8`).
> Variable-length sequence collections are supported via the `seqpro.rag` submodule (`Ragged` arrays)."

From the author-written `skills/seqpro/SKILL.md`, "When to use":

> "Encoding/decoding DNA, RNA, or protein sequences (OHE, integer tokens, padding). Sequence
> augmentation: reverse complement, k-mer shuffle, jitter, random draws. Sequence stats: GC content,
> nucleotide composition, length. Variable-length batches (e.g. peaks, transcripts of different sizes)
> → `sp.Ragged`. Genomic interval I/O: `sp.bed`, `sp.gtf`."

**Author-documentation caveat:** the quoted API generalizations are broader than the 0.22.0 code:
`Ragged` is reached as `sp.rag.Ragged`, not root `sp.Ragged`, and accepted input representations vary
by function (`pad_seqs` and `k_shuffle`, for example, have no `Ragged` overload). Here, **root export**
means a name in `seqpro.__all__`; declared submodule surfaces are identified separately.

**SeqPro is an ML-preprocessing library, and its authors say so.** The complete root export list
(`python/seqpro/__init__.py`, `__all__`, 28 names, re-derived verbatim for this memo) is:

```
AA, DNA, RNA, AminoAlphabet, NestedStr, NucleotideAlphabet, PathLike, SeqType, StrSeqType,
alphabets, bed, gtf, rag, transforms,
bin_coverage, cast_seqs, decode_ohe, decode_tokens, gc_content, jitter, k_shuffle,
length, nucleotide_content, ohe, pad_seqs, random_seqs, reverse_complement, tokenize
```

There is no mutagenesis, motif, barcode, primer, variant, codon-optimization, provenance, plotting,
or sequence-naming symbol in that list, and — apart from the two unmaintained plotting helpers in
`experimental/_visualizers.py` (see `design_visualization` below) — none anywhere in the 27 Python
source files of the recursive tree. Keyword greps over all 608 CHANGELOG lines (51 release headings) for
`mutagen|motif|barcode|primer|variant|vcf|hgvs|exon|intron|fasta|2bit|genome|provenance|plot|visual|design`
return **0 substantive hits**; the only `codon` hits (6×) are forward-translation kernel/LUT work
("O(1) LUT codon→AA lookup — 179× speedup on bench", plus LUT-indexing and codon-table-test entries).

**Framing for the referee response (important).** SeqPro is an excellently engineered package that is
not in the library-design problem space. Lead with the authors' own scope statement above, credit the
engineering (Rust/PyO3 + maturin, Numba, CodSpeed CI benchmarking, pandera/narwhals schema validation,
`py.typed`, `hypothesis[numpy]` property tests, 10 releases in six weeks, 49 design/spec documents
under `docs/superpowers/` plus a repo `CLAUDE.md`), and do not let a column of "no"s read as a
deficiency claim.

---

## Capability assessment

### Block A — library specification

**`library_as_object` — partial.**
A *batch* of sequences is genuinely a first-class object, and the design philosophy is
vectorized-over-the-batch rather than user loops (`SKILL.md`: *"No Python loops over sequences in
library code."*). Fixed-length batches are `(N, L)` NumPy `S1`/`uint8` arrays. Variable-length
collections are `sp.rag.Ragged` (`rag/_core.py`: `class Ragged(NDArrayOperatorsMixin, Generic[RDTYPE_co])`,
Rust-native, Arrow-style flat buffer + offsets, one or two ragged axes — R≥3 rejected — NumPy ufunc dispatch), with
`Ragged.from_lengths`, `Ragged.from_fields`, `.fields`, `.to_ak()`, `.to_chars()`, `.to_strings()`.
`sp.rag.zip({"seq": seq_rag, "score": score_rag})` (an alias for `Ragged.from_fields`) builds a
**record-layout** Ragged carrying multiple aligned per-sequence fields. These are real, composable,
first-class batch containers.
**But** they are preprocessing data containers, not *library specifications*: they hold materialized
sequences with no notion of design steps, no per-sequence build history, and no names. Nothing in the
object records or reproduces how the sequences were designed.
*Framing note:* "partial" is a generous reading, generous in the safe direction. Scoring "no" would be
attackable by an author-referee pointing at `Ragged` + record layout.
*Source:* `python/seqpro/rag/_core.py`, `rag/__init__.py`, `skills/seqpro/SKILL.md`, `docs/ragged.md`.

**`dag_chaining` — partial.**
`python/seqpro/transforms/augmentation.py` defines `Sequential(*transforms)` whose `__call__` is
verbatim `for t in self.transforms: x = t(x)` / `return x`, and `Random(p, *transforms, seed=...)`
which applies a nested chain with probability `p`. The declared transform surface is
`transforms/__all__ = ["TMM", "Jitter", "KShuffle", "Random", "ReverseComplement", "Sequential"]`.
This is a **linear** augmentation pipeline over arrays: composition and nesting
exist, but there is no graph, no branching or merging, no named nodes, no design semantics, and the
operands are augmentations applied to already-materialized sequences, not design operations.
**Correction vs. the original extraction:** the transform list is *not* "ReverseComplement, KShuffle,
Jitter, Tokenize". `Tokenize` is defined in `augmentation.py` (line 114) but is **omitted from
`transforms/__all__`**; `TMM` is in `__all__` and was missing from the original list. Re-verified from
`transforms/__init__.py`.
*Source:* `python/seqpro/transforms/__init__.py`, `python/seqpro/transforms/augmentation.py`.

**`lazy_evaluation` — no.** *(evidence rewritten — the original text contained two errors)*
No sequence-producing function is lazy: `random_seqs`' body is literally
`return seed.choice(alphabet.array, size=shape)`, and every encoder/modifier returns a materialized
NumPy array or `Ragged`. There is no design step whose construction could be deferred.
Deferred execution does exist in the package, but only for I/O and array transforms, never for
library construction: `sp.gtf.scan(path)` returns a polars `LazyFrame`, and the experimental
`seqpro.xr` module runs `ohe`, `bin_coverage`, and `translate` through
`xr.apply_ufunc(..., dask="parallelized")` — genuinely lazy, chunk-parallel evaluation over
dask-backed sequence arrays.
**Two corrections vs. the original extraction, both material:** (1) the claim that "the only laziness
in the package is annotation-table I/O" is **false** because of `seqpro.xr` (re-verified: `xr/__init__.py`
lines 42–48, 93–99, 147–152 all pass `dask="parallelized"`); (2) BED reading is **not** lazy —
`bed.read` dispatches to `_read_bed`/`_read_narrowpeak`/`_read_broadpeak`, all of which call the eager
`pl.read_csv` (lines 236, 287, 356). The "no" verdict is unaffected: laziness of *chunked array
compute* is not lazy generation of a designed library.
*Source:* `python/seqpro/_modifiers.py:281` (`random_seqs`), `python/seqpro/xr/__init__.py`,
`python/seqpro/gtf.py`, `python/seqpro/bed.py`.

**`mixed_mutagenesis_one_pool` — no.** No mutagenesis of any kind exists. Re-derived `__all__`, walked
the untruncated recursive tree (190 entries / 162 blobs, 27 Python files), grepped every Python module and all 608
CHANGELOG lines: zero mutagenesis machinery. The nearest things are `k_shuffle` (k-let-preserving
shuffle) and `jitter` (random offset crop) — both augmentations, neither parameterized by mutation
type, and neither able to produce a single annotated pool mixing mutation classes.
*Source:* `python/seqpro/__init__.py`, `python/seqpro/_modifiers.py`, `CHANGELOG.md`, recursive git tree.

**`combinatorial_motif_place` — no.** No motif-insertion/implantation function exists. A source-wide
grep for `motif` across every fetched module returns exactly one hit — the unfinished comment
`# Performing motif anal` in `experimental/_experimental.py`, no code; CHANGELOG returns 0 hits; `docs/api/`
contains only `index, alphabets, bed, gtf, ragged, types`.
*Source:* source-wide grep, `docs/api/` listing.

**`barcode_generation` — no.** `sp.random_seqs(shape, alphabet, seed)` draws i.i.d. random nucleotides
with **no constraints** (no distance, GC, homopolymer, or uniqueness filter); `sp.gc_content` only
measures GC and cannot filter or enforce it; nothing attaches a barcode to a sequence. An
`edit_distance(seq1, seq2, dual=False)` exists only in `python/seqpro/experimental/_experimental.py`,
which is (a) headed `# Unmaintained experimental module, not part of the public package surface.`,
(b) a **Hamming** distance requiring equal-length inputs (`assert len(seq1) == len(seq2), "Both
sequences must be of same length."`), not an edit distance, and (c) not wired into any pool
construction routine.
**Strengthened vs. the original extraction:** `python/seqpro/experimental/` has **no `__init__.py` at
all** — the raw URL returns HTTP 404 and the path is absent from the recursive tree (which contains
only `experimental/_experimental.py` and `experimental/_visualizers.py`). `seqpro.experimental` is
therefore an implicit namespace package — undocumented, not root-exported and outside the declared
public surface — though the 0.22.0 wheel does ship both modules and
`seqpro.experimental._experimental` still imports.
*Source:* `python/seqpro/_modifiers.py:281`, `python/seqpro/_analyzers.py:32`,
`python/seqpro/experimental/_experimental.py`, recursive git tree, raw-URL 404.

**`per_sequence_provenance` — no.** Nothing records how a sequence was built, and there is no design
step whose provenance could exist. The record-layout `Ragged` (`sp.rag.zip` / `Ragged.from_fields`)
is a generic mechanism that lets a *user* carry arbitrary parallel per-sequence fields, and
`rag.hash("sha256"|"md5"|"rapidhash")` can fingerprint sequences — but SeqPro itself writes no
provenance and no operation in the package emits build metadata.
*Source:* `python/seqpro/rag/__init__.py`, `python/seqpro/rag/_ops.py:410`.

**`automatic_naming` — no.** No naming or ID-generation function anywhere; outputs are anonymous array
rows. Two robustness notes so the cell cannot be nitpicked: `bed.py`'s `BEDSchema` *reads* a `name`
column from a BED file, and a record-layout `Ragged` exposes `.fields` (field names, not per-sequence
names) — neither generates per-sequence identifiers.
**Correction vs. the original extraction:** `docs/api/index.md` lists **14** top-level members, not 15
(counted from the rendered file: `bin_coverage, cast_seqs, decode_ohe, decode_tokens, gc_content,
jitter, k_shuffle, length, nucleotide_content, ohe, pad_seqs, random_seqs, reverse_complement,
tokenize`).
*Source:* `python/seqpro/__init__.py`, `docs/api/index.md`, `python/seqpro/bed.py`, `rag/_core.py:303`.

**`design_visualization` — no.** No plotting on the root or declared submodule surfaces and no
plotting dependency.
`python/seqpro/experimental/_visualizers.py` defines `plot_gc_content` and `plot_nucleotide_content`,
but the file's own header reads:
> "Unmaintained experimental module, not part of the public package surface."
> "NOTE: this module is currently broken at import time — `gc_content_seqs` and `nucleotide_content_seqs`
> were renamed to `gc_content` / `nucleotide_content` with different signatures."

Confirmed: it imports two functions that no longer exist in `_analyzers.py`, so importing the module
raises; its containing directory also has no `__init__.py` (404), leaving it an undeclared namespace
package. Even if working, these are QC plots of a sequence set — a GC histogram and a positional
nucleotide-frequency line plot — not a view of a design.
*Source:* `python/seqpro/experimental/_visualizers.py` (verbatim header), `_analyzers.py`, recursive tree.

### Block B — assay coverage

**`assay_dms` — no.** No coding-variant or DMS machinery: no variant table, no codon-substitution
enumeration, no mutant-library construction. `sp.AA.translate()` translates DNA→AA (the built-in `AA`
table is the DNA codon table, so RNA codons containing `U` fall to the `unknown` policy) with a Rust
64-entry codon LUT, `unknown="X"|"drop"` policies, and `validate=True` fast-fail.
**Pre-rebuttal (added; a referee may say "SeqPro handles stop codons"):** `AA.translate(...,
truncate_stop=True)` truncates each output at the first stop codon **inclusive** (Ragged input only,
applied after drop-compaction), backed by a dedicated Rust kernel `_translate_stop_ends` and an OHE
path that locates the `*` one-hot column. That is more coding-sequence awareness than a bare "the `AA`
alphabet includes `*`" — and it is still not DMS design: truncating a translated sequence at a stop
codon enumerates nothing, mutates nothing, and produces no variant library.
*Source:* `python/seqpro/alphabets/_alphabets.py:421–720` (`translate` overloads, `truncate_stop`
docstring at 485–487, `_translate_stop_ends` at 648 and 714), `src/translate.rs`.

**`assay_mpra` — no.** No regulatory-library design of any kind: no promoter/enhancer scaffolding, no
barcoding, no adapter or cloning handles, no motif perturbation, no oligo assembly. Confirmed against
`SKILL.md`'s "When to use" list, `__all__`, and the full source tree.
*Source:* `python/seqpro/__init__.py`, `skills/seqpro/SKILL.md`, recursive git tree.

**`assay_insilico` — partial.** SeqPro supplies preprocessing primitives used when probing
sequence-to-function models: one-hot encoding/decoding, integer tokenization, padding, `k_shuffle`
(k-let-preserving shuffling, dispatched to the Rust kernel `src/kshuffle.rs` and CodSpeed-benchmarked),
`random_seqs` (random backgrounds), `jitter`
(which applies the same random crop across several aligned arrays, keeping a sequence and its per-base
tracks in register) and `ReverseComplement` augmentations, and `Ragged` batching. Explicitly framed as ML infrastructure:
*"heavily utilized by ... SeqExplainer and EUGENe."*
It does **not** design in-silico libraries, run models, do ISM, or score anything: no model-scoring,
ISM, or in-silico library-construction code exists in SeqPro.
*Judgment call*, flagged in confidence notes: primitives-only.
*Source:* `python/seqpro/_modifiers.py:39` (`k_shuffle`), `src/kshuffle.rs`, `README.md`.

### Block C — genomics integration

**`genome_coordinates` — partial.** *(see the flagged legend dependency below)*
`python/seqpro/bed.py` provides `read` (dispatching to `_read_bed`, `_read_narrowpeak`,
`_read_broadpeak`), `sort`, `with_len`, and `to_pyr`/`from_pyr` (pyranges interop);
`python/seqpro/_coords.py` provides `CoordSchema`, `detect_schema`, `set_schema` with **5 built-in
schemas** (bed / pb / pr / gtf / gff), each carrying a 0- vs 1-based flag (there is no open/closed
field, and `set_schema` renames the coordinate columns and tags the frame rather than converting values).
So SeqPro ingests and manipulates genomic intervals competently — **but it does not fetch sequence
from a reference genome**: a grep across the whole tree finds no FASTA, 2bit, or bigWig reader, and
there is no coordinate→sequence resolution.
**Added (missed by the original extraction):** `with_len(bed, length)` is **narrowPeak-summit aware** —
when a `peak` column is present it recenters on `2 * (chromStart + peak)` (falling back to
`chromStart + chromEnd`, the interval midpoint, when `peak` is null or absent), then sets
`chromStart = (double_center - length) // 2`, `chromEnd = (double_center + length) // 2`. This is a
genuine correctness nicety, in the same class as the coordinate-schema handling.
*Source:* `python/seqpro/bed.py:38–164, 219–378`, `python/seqpro/_coords.py`, source-wide grep.

**`transcript_models` — partial.** `sp.gtf.scan(path)` parses GFF/GTF into a polars `LazyFrame` with
columns `seqname, source, feature, start, end, score, strand, frame, attribute`, and
`sp.gtf.attr("gene_id")` returns a polars expression regex-extracting one attribute as Utf8. That is
the **entire** GTF module: `gtf.py` is 1,813 bytes with exactly two public functions, and
`docs/api/gtf.md` lists only those two. There is **no transcript object, no exon assembly, no CDS/UTR
resolution, and no gene-model semantics** — a user gets a flat annotation table and must build any
transcript logic themselves.
*Framing note:* "partial" is generous (a flat annotation table is not a gene model), generous in the
safe direction.
*Source:* `python/seqpro/gtf.py`, `docs/api/gtf.md`.

**`exon_intron_split_codons` — no.** Nothing joins exons or tracks reading frame across exon
boundaries. `pl.col("frame").fill_null(0)` in `gtf.py` is the **only** use of `frame` anywhere in the
package — it is parsed, defaulted, and never consumed. `translate()` operates on a contiguous
nucleotide buffer in stride-3 codons; `truncate_stop` likewise works within a single buffer and never
crosses an exon boundary.
*Source:* `python/seqpro/gtf.py`, `python/seqpro/alphabets/_alphabets.py`, `src/translate.rs`,
source-wide grep for `frame`.

**`hgvs_input` — no.** Zero hits for `hgvs` in the full recursive tree, in every fetched module, in the
docs, and across all 51 release sections of `CHANGELOG.md`. Independently reproduced.
*Source:* source-wide grep, `CHANGELOG.md`, recursive git tree.

**`vcf_vep_output` — no.** No VCF reader/writer, no VEP-compatible export, no variant representation of
any kind. Package I/O is BED / narrowPeak / broadPeak / GTF-GFF in, NumPy / polars / pyranges /
awkward out. Confirmed via source grep and the PyPI `requires_dist` list — the runtime dependencies are
`numba, numpy, polars, pyranges, pandera, pandas, pyarrow, natsort, narwhals, setuptools, awkward,
polars-config-meta[polars], attrs`, none of which is a variant-format library.
**Correction vs. the original extraction:** the dependency list omitted `polars-config-meta[polars]>=0.3.2`
(used to tag BED frames with `coordinate_system_zero_based` metadata). The claim is unaffected.
*Source:* `pyproject.toml`, PyPI JSON, source-wide grep.

**`consequence_annotation` — no.** There is no variant representation anywhere in the package, hence
nothing to compare against a reference allele and nothing to classify. `AA.translate` can emit a stop
codon `*` and `truncate_stop=True` can truncate at the first one, but finding a stop codon in a
translated sequence is not classifying a variant consequence: there is no reference allele, no variant
object, and no consequence class.
*Source:* `python/seqpro/alphabets/_alphabets.py`, source-wide grep for `variant|consequence`.

### Block D — physical construction

**`primer_design` — no.** A source-wide grep for `primer|oligo|melting|Tm|restriction|overhang|assembl`
returns 0 substantive hits; the CHANGELOG returns 0. No primer, oligo, melting-temperature, overhang,
or assembly functionality. Independently reproduced.
*Source:* source-wide grep, `CHANGELOG.md`.

**`codon_optimization` — no.** Codon tables exist strictly for **forward** translation:
`AminoAlphabet(codons, amino_acids)` builds `self.codon_to_aa = dict(zip(codons, amino_acids))` and a
64-entry LUT via `_build_translate_lut` (gated by `_can_build_lut`). There is **no reverse map, no
codon-usage frequency table, and no optimizer**. All 6 `codon` hits in the CHANGELOG are
forward-translation kernel/LUT work (kernels, LUT indexing, a codon-table test, a benchmark) — re-grepped.
*Source:* `python/seqpro/alphabets/_alphabets.py:246, 297–350`, `CHANGELOG.md`.

**`synthesis_constraints` — no.** All analyzers are read-only and there is no homopolymer, repeat,
restriction-site, or GC-window check and no enforcement/repair loop. `gc_content`, `nucleotide_content`,
and `length` measure; `_cleaners.py` provides `remove_N_seqs`, `remove_only_N_seqs`, `remove_whitespace`,
`sanitize` — input hygiene, not synthesis feasibility.
**Two corrections vs. the original extraction:** (1) `count_kmers_seq` is **not** public — it is absent
from `__all__` and not imported by `__init__.py`, so `sp.count_kmers_seq` raises `AttributeError`
(reachable only as `seqpro._analyzers.count_kmers_seq`), and its batched sibling `_analyzers._count_kmers`
is `raise NotImplementedError` with a `# TODO: non-trivial to parallelize/SIMD this` comment above it;
(2) `remove_whitespace`'s entire body is `pass` — it is a no-op stub. Note also that the whole
`_cleaners` module is private (none of its functions appear in `__all__`). All three facts reinforce
the "no".
*Source:* `python/seqpro/_analyzers.py:138, 165–187`, `python/seqpro/_cleaners.py`, `python/seqpro/__init__.py`.

### Block E — engineering

**`interface` — yes (Python API only).** `import seqpro as sp`. `pyproject.toml` declares **no**
`[project.scripts]` and no console entry points — there is no CLI, no GUI, no web service.
`requires-python = ">=3.10"`; `py.typed` shipped; build backend `maturin` (Rust extension via PyO3).
Docs site (MkDocs/zensical) at https://ml4gland.github.io/SeqPro/, which is also the repo's declared
GitHub `homepage`.
*Source:* `pyproject.toml`, GitHub repo metadata.

**`license` — yes (MIT).** GitHub API `license.spdx_id = "MIT"`; `LICENSE` at repo root;
`pyproject.toml` `license = { file = "LICENSE" }`.
*Source:* `api.github.com/repos/ML4GLand/SeqPro`, `pyproject.toml`, `LICENSE`.

**`maintained` — yes (ongoing upstream commit/package-release activity).** Last commit
`2026-07-27T19:02:07Z` ("bump: version 0.21.2 →
0.22.0"); repo `pushed_at 2026-07-27T19:07:56Z`; release 0.22.0 tagged `2026-07-27T19:02:25Z` and
uploaded to PyPI `2026-07-27T19:07:05`. 55 PyPI releases since `2023-03-26T17:33:00`; 10 GitHub
releases 0.16.0 → 0.22.0 in ~6 weeks. Not archived, 14 stars, 1 fork, 0 open issues. All dates
independently re-confirmed for this memo (2026-08-14).
*Source:* GitHub repo/commits/releases API, PyPI JSON.

---

## Documented examples / vignettes

SeqPro ships **no library-design vignettes** — the examples are API demos, not use cases.

1. **README / docs "Quick Start" snippet** (`README.md`; `docs/index.md`) — the repository's main example:
   `sp.random_seqs(shape=(2,9), alphabet=sp.DNA, seed=1234)`, then `pad_seqs`, `DNA.ohe`/`decode_ohe`,
   `tokenize`/`decode_tokens`, `reverse_complement`, `k_shuffle(k=2)`, `gc_content`,
   `nucleotide_content`, `jitter(max_jitter=2)`, `bin_coverage(bin_width=128)`, `AA.translate(cds)`,
   `AA.ohe`. As written, two of those calls do not run: `k_shuffle` omits the required `alphabet`
   argument and `jitter` omits the required `jitter_axes`, so both raise `TypeError`.
2. **"Ragged Arrays" guide** (`docs/ragged.md`; https://ml4gland.github.io/SeqPro/ragged/) — a worked
   example building `Ragged.from_lengths` over "per-position coverage for three intervals of different
   widths", ufunc arithmetic, and record layouts (built in the guide with `ak.zip`, which returns a
   `Ragged`).
3. **API Reference pages** (`docs/api/{index,alphabets,bed,gtf,ragged,types}.md`) — mkdocstrings-rendered
   signatures, docstrings and source (`show_source: true`), with no narrative examples.
   `docs/api/index.md` lists **14** top-level members (corrected from "15").
4. **`skills/seqpro/SKILL.md`** — an author-maintained agent-facing cheat sheet with a task→call quick
   reference table plus Ragged construction, record-layout, and hashing examples. Effectively the
   richest usage document in the repo.
5. **Test suite as usage examples** — `SKILL.md` points to `tests/` ("Tests as usage examples");
   ~30 test modules, dominated by `test_ragged*`, `test_ohe`, `test_tokenize`, `test_translate`,
   `tests/bed/*`; `hypothesis[numpy]` property-based tests among them.
6. **`docs/superpowers/` design & spec documents** — 49 markdown files (plus a repo-root `CLAUDE.md`)
   documenting internal design decisions. Not user tutorials, but they are the repo's densest
   statement of intent and scope.
7. **No repository notebook** — the complete SeqPro tree contains no `.ipynb` file.

**Reproducibility note for PoolParty:** none of these is a library-design use case. There is nothing
here to "reproduce" in PoolParty beyond trivially generating random sequences and computing GC — which
is the honest headline for the comparison table.

---

## Notable capabilities not in the row list

- **`seqpro.transforms.TMM` — TMM count normalization.**
  *(Added; the single biggest omission from the original extraction.)* `python/seqpro/transforms/tmm.py`
  (205 lines) implements a full TMM estimator: a `TMM` class with `fit(counts)` /
  `transform(counts, library_size=1e6)`, quantile-based reference-sample selection
  (`expression_cutoff, log_ratio_trim=0.3, expression_trim=0.05, apply_weighting, quantile=0.75`),
  log-ratio and absolute-expression trimming, inverse-variance weighting, and an
  `@nb.njit(parallel=True, nogil=True, cache=True)` `_tmm_helper` kernel. It **is** exported
  (`transforms/__all__` includes `"TMM"`). It creates no Block A–D "yes" (TMM normalizes count
  matrices; it designs nothing), but it shows the package reaching into assay-count normalization.
- **One-hot encoding / decoding / integer tokenization / padding**, vectorized over batches
  (`ohe`, `decode_ohe`, `tokenize`, `decode_tokens`, `pad_seqs`), Numba-accelerated; alphabet objects
  also carry `.tokenize` / `.decode_tokens` methods (`NucleotideAlphabet`, overloaded).
- **k-let-preserving shuffling** (`k_shuffle`) implemented in Rust (`src/kshuffle.rs`, with
  `src/kshuffle_ref.rs` and `benches/kshuffle.rs`).
- **`Ragged` variable-length array container** — Rust-native (`crates/seqpro-core/src/ragged.rs`,
  `src/ragged.rs`), Arrow-style flat data + offsets, one or two ragged axes (R≥3 rejected),
  `NDArrayOperatorsMixin` ufunc dispatch, `from_lengths`/`from_fields`/`.fields`,
  `to_padded`/`to_packed`/`concatenate`, `np.memmap`-safe, awkward-array interop via `to_ak()`,
  and `to_chars()` / `to_strings()` for opaque-string ↔ `S1`-char conversion.
- **Multi-field record layout** via `sp.rag.zip` (alias of `Ragged.from_fields`) — parallel
  per-sequence fields in one container.
- **Parallel Rust string hashing** of ragged batches: `rag.hash("sha256"|"md5"|"rapidhash")`
  (`src/hashing.rs`).
- **Reverse complement across every representation**: on `str`, `S1` bytes, and one-hot arrays
  (`sp.reverse_complement`), on coverage tracks (`transforms.ReverseComplement(*types=("dna","track"))`),
  **and Ragged-native** via `sp.rag.reverse_complement` (`rag/_ops.py:22`, in `rag/__all__`) —
  *the Ragged variant was missed by the original extraction.*
- **Coverage utilities**: `bin_coverage(bin_width, normalize)` (public) and
  `normalize_coverage(method="CPM"|"CPKM", total_counts)` — the latter is **private**
  (`seqpro._modifiers.normalize_coverage`; not in `__all__`, so `sp.normalize_coverage` raises
  `AttributeError`).
- **k-mer counting** `count_kmers_seq(seq, k)` — also **private** (`seqpro._analyzers`), and its
  batched sibling `_count_kmers` is `raise NotImplementedError`.
- **Custom alphabets**: `NucleotideAlphabet(alphabet, complement)`, `AminoAlphabet(codons, amino_acids)`
  — user-definable, including non-standard codon tables.
- **Translation with explicit policies**: `AA.translate(..., validate=False, unknown="X"|"drop",
  truncate_stop=False)`; `"drop"` returns a `Ragged` because lengths become variable; `truncate_stop=True`
  truncates at the first stop codon inclusive (Ragged only, Rust kernel `_translate_stop_ends`).
- **Coordinate-schema handling** (`_coords.py`): `CoordSchema`, `detect_schema`, `set_schema` with 5
  built-in schemas (bed/pb/pr/gtf/gff) carrying a 0- vs 1-based flag; `polars-config-meta` tags frames
  with `coordinate_system_zero_based`.
- **Peak-format readers**: BED, narrowPeak, broadPeak (pandera/narwhals-validated schemas); pyranges
  and polars/narwhals interop; **summit-aware `with_len`** (recenters on the narrowPeak `peak` column).
- **Transform objects** for ML dataloaders: `Sequential`, `Random`,
  `ReverseComplement`, `KShuffle`, `Jitter`, `TMM` (declared surface); `Tokenize` exists in the module
  but is *not* in `transforms/__all__`.
- **Experimental XArray/dask integration** (`seqpro.xr`): `ohe` and `bin_coverage`
  (in `xr.__all__`) plus `translate` (defined, not exported) all run through
  `xr.apply_ufunc(..., dask="parallelized")` — lazy, chunk-parallel, out-of-core.
- **Engineering rigor**: Rust/PyO3 + maturin (10 `.rs` files incl. a `seqpro-core` crate), Numba,
  CodSpeed continuous benchmarking, pandera/narwhals schema validation, `py.typed`, pyrefly strict
  type checking wired into pre-commit (`pixi run typecheck`) alongside a configured basedpyright,
  commitizen conventional-commit changelog, dev-deps
  `biopython>=1.85` and `hypothesis[numpy]>=6.131.2` (property-based testing), 49 design/spec docs
  under `docs/superpowers/`.

---

## Stated limitations / scope disclaimers

- README's final section heading is a bare **"# More to come!"** (followed only by a
  contributions/code-of-conduct paragraph) — the authors treat the surface as incomplete.
- `python/seqpro/experimental/*` is headed **"Unmaintained experimental module, not part of the public
  package surface."**, and the directory has **no `__init__.py`** (raw URL 404), so it is an implicit
  namespace package outside the declared public surface — the modules still ship in the wheel and
  `_experimental` imports. `_visualizers.py` adds **"NOTE: this module is currently broken at import
  time"** — it imports `gc_content_seqs` / `nucleotide_content_seqs`, which no longer exist.
- `seqpro.xr` is described as an **"experimental integration with XArray"**.
- README frames SeqPro as a dependency of a suite: **"fully functional on its own but is heavily
  utilized by other packages including SeqData, MotifData, SeqExplainer, and EUGENe"**; it does not
  assign responsibilities to those packages or state a formal scope boundary.
- `SKILL.md`: **"Public API: see `python/seqpro/__init__.py` for the full export list. Re-read it before
  assuming a symbol exists."**
- Several functions are stubs or private-only: `remove_whitespace` body is `pass`;
  `_analyzers._count_kmers` is `raise NotImplementedError` (`# TODO: non-trivial to parallelize/SIMD`);
  `count_kmers_seq`, `normalize_coverage`, and the whole `_cleaners` module are not in `__all__`.
- 0.22.0 carries a **BREAKING CHANGE** entry for `Ragged.__getitem__` / record-row-peeling semantics
  (the only such entry in the changelog) — the API is still moving (`major_version_zero = true` in
  commitizen config).

---

## Availability status

**Fully available; `maintained` here means ongoing upstream commit/package-release activity.**
`pip install seqpro` → **0.22.0**, uploaded to PyPI **2026-07-27T19:07:05**, `requires-python >=3.10`.
Runtime deps: `numba>=0.58.1, numpy>=1.26, polars>=1.21,<2, pyranges>=0.1.3,<0.2, pandera>=0.31.1,
pandas, pyarrow, natsort, narwhals>=2.20, setuptools>=70, awkward>=2.5, polars-config-meta[polars]>=0.3.2,
attrs`. Repo public, **not archived**, MIT, **0 open issues**, 14 stars, 1 fork, last commit
**2026-07-27T19:02:07Z**, `pushed_at 2026-07-27T19:07:56Z`. 55 PyPI releases since 2023-03-26; 10 GitHub
releases (0.16.0 → 0.22.0) in the six weeks before assessment. Docs site live at
https://ml4gland.github.io/SeqPro/. No CLI, GUI, or web service exists, so there is none to be alive or
dead. Wheels are built via maturin for manylinux/musllinux/windows/macos — a *source* install requires a
Rust toolchain, but published wheels make `pip install` work without one. (All metadata re-verified
2026-08-14.)

---

## Unresolved disagreements and flagged judgment calls

No capability value is recorded as unknown; all 24 retain their prior assessments, and no value
changed. The following are flagged for the table build:

1. **`genome_coordinates = partial` is legend-dependent (the one real framing risk).** Under the
   intended reading — *resolves genomic coordinates to sequence* — "partial" is correct: SeqPro parses,
   validates, sorts, resizes, and converts intervals but has no FASTA/2bit reader (grep-confirmed) and
   cannot produce sequence from coordinates. Under a looser reading — *accepts genomic intervals* —
   SeqPro is a clean **yes**, and the cell would read as an understatement an author-referee could
   attack. **Action: the row legend must state coordinate→sequence resolution explicitly.**
2. **Four "partial"s are deliberately generous, all in the safe direction** (crediting SeqPro rather
   than under-crediting it): `library_as_object` (`Ragged` + record layout are genuinely first-class
   composable containers, but hold no design steps, history, or names), `dag_chaining` (linear
   `Sequential`/`Random` composition, no graph), `transcript_models` (a flat annotation table is not a
   gene model), and `assay_insilico` (primitives only — no model scoring, ISM, or library
   construction). A referee could argue "no" for any of them; none can argue "yes".
3. **`lazy_evaluation = no` is a scoped claim.** SeqPro does contain real laziness (`gtf.scan` →
   `LazyFrame`; `seqpro.xr` via `xr.apply_ufunc(dask="parallelized")`). The "no" refers strictly to
   lazy *construction of a designed library*, of which there is none — there is no design step to defer
   and `random_seqs` is eager. The original extraction's blanket sentence ("the only laziness in the
   package is annotation-table I/O") was false and has been removed.
4. **Everything in Blocks B (dms/mpra), C (exon/hgvs/vcf/consequence), and D is a high-confidence "no"**,
   supported by four repository checks: the complete root `__all__` export list, the untruncated
   recursive file tree, the rendered docs navigation, and keyword greps over all 51 release sections
   of `CHANGELOG.md` plus every Python module (0 substantive hits for
   mutagen/motif/barcode/primer/oligo/restriction/vcf/hgvs/variant/consequence/library design).
5. Source was read from `main` (0.22.0). The package was **not installed**; all claims are source-,
   docs-, and metadata-derived.
