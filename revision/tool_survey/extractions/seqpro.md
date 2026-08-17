# SeqPro — evidence memo

**Slug:** `seqpro`
**Full name:** SeqPro (Sequence processing toolkit) — "Genomic sequence preprocessing toolkit"
**Authors:** David Laub (dlaub@ucsd.edu), Adam Klie (aklie@ucsd.edu) — ML4GLand org (UCSD)
**Paper:** none. Software only. (Confirmed by web search; SeqPro is a support library of the ML4GLand / EUGENe suite and has no standalone publication.)

## Sources consulted

| Kind | URL / ref | Retrieved |
|---|---|---|
| repo | https://github.com/ML4GLand/SeqPro | 2026-08-10 |
| repo (raw) | https://raw.githubusercontent.com/ML4GLand/SeqPro/main/README.md | 2026-08-10 |
| repo (raw) | `python/seqpro/__init__.py` (full public API surface) | 2026-08-10 |
| repo (raw) | `python/seqpro/_modifiers.py`, `_analyzers.py`, `_cleaners.py`, `_coords.py` | 2026-08-10 |
| repo (raw) | `python/seqpro/gtf.py`, `python/seqpro/bed.py` | 2026-08-10 |
| repo (raw) | `python/seqpro/transforms/augmentation.py` | 2026-08-10 |
| repo (raw) | `python/seqpro/experimental/_experimental.py`, `_visualizers.py` | 2026-08-10 |
| repo (raw) | `python/seqpro/alphabets/_alphabets.py` (signatures) | 2026-08-10 |
| repo (raw) | `skills/seqpro/SKILL.md` (author-written API cheat-sheet) | 2026-08-10 |
| repo (raw) | `CHANGELOG.md`, `pyproject.toml`, `docs/index.md`, `docs/api/*.md`, `docs/ragged.md` | 2026-08-10 |
| repo API | `api.github.com/repos/ML4GLand/SeqPro` (metadata, commits, releases, full file tree) | 2026-08-10 |
| docs | https://ml4gland.github.io/SeqPro/ | 2026-08-10 |
| pypi | https://pypi.org/pypi/seqpro/json | 2026-08-10 |
| web | https://github.com/ML4GLand (org repo listing incl. `tutorials`, `use_cases`) | 2026-08-10 |

## What SeqPro actually is

From `docs/index.md` and README (verbatim):

> "SeqPro is a Python package for processing DNA/RNA sequences, with support for protein sequences. SeqPro is fully functional on its own but is heavily utilized by other packages including SeqData, MotifData, SeqExplainer, and EUGENe."

> "All functions in SeqPro take as input a string, a list of strings, a NumPy array of strings, a NumPy array of single character bytes (`S1`), or a NumPy array of one-hot encoded sequences (`uint8`). Variable-length sequence collections are supported via the `seqpro.rag` submodule (`Ragged` arrays)."

> "Computational bottlenecks or code that is impossible to vectorize with NumPy alone are accelerated with Numba ... and Rust via PyO3/maturin (e.g. k-mer shuffling)."

**It is an ML-preprocessing library, not a library-design tool.** The complete public API (`python/seqpro/__init__.py`, `__all__`) is:

```
AA, DNA, RNA, AminoAlphabet, NucleotideAlphabet, NestedStr, SeqType, StrSeqType, PathLike,
alphabets, bed, gtf, rag, transforms,
bin_coverage, cast_seqs, decode_ohe, decode_tokens, gc_content, jitter, k_shuffle,
length, nucleotide_content, ohe, pad_seqs, random_seqs, reverse_complement, tokenize
```

There is no mutagenesis, motif-insertion, barcode, primer, codon-optimization, variant, or
sequence-naming symbol anywhere in that list. Grepping `CHANGELOG.md` (55 releases of history)
for `mutagen|motif|barcode|primer|codon optim|library design|variant|vcf|hgvs` returns **0 hits**;
the only `codon` hits are translation-kernel optimizations ("O(1) LUT codon→AA lookup — 179× speedup on bench").

The author-written `skills/seqpro/SKILL.md` states the intended scope directly:

> "**When to use** — Encoding/decoding DNA, RNA, or protein sequences (OHE, integer tokens, padding). Sequence augmentation: reverse complement, k-mer shuffle, jitter, random draws. Sequence stats: GC content, nucleotide composition, length. Variable-length batches (e.g. peaks, transcripts of different sizes) → `sp.Ragged`. Genomic interval I/O: `sp.bed`, `sp.gtf`."

## Capability assessment

### Block A — library specification

**library_as_object — partial.**
A *batch* of sequences is genuinely a first-class object, and the whole design philosophy is
vectorized-over-the-batch rather than user loops. `SKILL.md`: *"No Python loops over sequences in
library code."* / *"Don't write Python `for` loops over residues or positions."* Fixed-length batches
are `(N, L)` NumPy `S1`/`uint8` arrays; variable-length collections are `sp.rag.Ragged`
(Arrow-style flat buffer + offsets, exactly one ragged axis, NumPy ufunc dispatch via
`NDArrayOperatorsMixin`). `sp.rag.zip({"seq": seq_rag, "score": score_rag})` builds a *record-layout*
Ragged carrying multiple aligned per-sequence fields.
**But** this is a data container for preprocessing, not a *library specification*: it holds
materialized sequences, has no notion of design steps, no per-sequence build history, no names.
Nothing in the object records or reproduces how the sequences were designed.

**dag_chaining — partial.**
`python/seqpro/transforms/augmentation.py` defines `Sequential(*transforms)` (applies callables
left-to-right) and `Random(p, *transforms, seed=...)` (applies a nested chain with probability `p`),
plus transform objects `ReverseComplement`, `KShuffle`, `Jitter`, `Tokenize`. This is a
torchvision-style **linear** augmentation pipeline over arrays — composition and nesting exist, but
there is no graph, no branching/merging, no named nodes, no design semantics, and the operands are
augmentations applied to already-materialized sequences, not design operations.

**lazy_evaluation — no (for sequences).**
Every sequence-producing function returns a materialized NumPy array; e.g. `random_seqs` is
`return seed.choice(alphabet.array, size=shape)` — eager. The only laziness in the package is
annotation-table I/O: `sp.gtf.scan(path) -> pl.LazyFrame` (polars `scan_csv`), and BED reading is
polars-backed. That is lazy *table* scanning, not lazy sequence generation.

**mixed_mutagenesis_one_pool — no.** No mutagenesis of any kind in the public API. The nearest
things are `k_shuffle` (k-let-preserving shuffle) and `jitter` (random offset crop), both
augmentations, not mutagenesis; neither is parameterized by mutation type and they cannot be mixed
into a single annotated pool.

**combinatorial_motif_place — no.** No motif-insertion/implantation function exists. Motif handling
in the ML4GLand suite lives in the *separate* packages `MotifData` (last push 2023-07-21) and
`SeqExplainer`, not in SeqPro. `docs/api/` contains only `alphabets`, `bed`, `gtf`, `ragged`, `types`.

**barcode_generation — no.** `sp.random_seqs(shape, alphabet, seed)` draws i.i.d. random nucleotides
with no constraints; `sp.gc_content` measures GC but does not filter or enforce it; nothing attaches
a barcode to a sequence. An `edit_distance(seq1, seq2, dual=False)` function exists only in
`python/seqpro/experimental/_experimental.py`, which is (a) not imported in `__init__.py`, (b) headed
`"# Unmaintained experimental module, not part of the public package surface."`, and (c) is a
Hamming distance requiring equal-length inputs (`assert len(seq1) == len(seq2)`), not an edit distance,
and is not wired into any pool-construction routine.

**per_sequence_provenance — no.** Nothing records how a sequence was built. The record-layout
`Ragged` (`sp.rag.zip`) is a generic mechanism that lets a *user* carry arbitrary parallel per-sequence
fields, and `rag.hash("sha256"/"md5"/"rapidhash")` can fingerprint sequences — but SeqPro itself
writes no provenance and no operation emits build metadata.

**automatic_naming — no.** No naming or ID-generation function; outputs are anonymous array rows.

**design_visualization — no.** No plotting in the public API. `python/seqpro/experimental/_visualizers.py`
has `plot_gc_content` and `plot_nucleotide_content`, but the file's own header says:
> "Unmaintained experimental module, not part of the public package surface."
> "NOTE: this module is currently broken at import time — `gc_content_seqs` and `nucleotide_content_seqs`
> were renamed to `gc_content` / `nucleotide_content` with different signatures."
Even if working, these are QC histograms of a sequence set, not a view of a design.

### Block B — assay coverage

**assay_dms — no.** No coding-variant or DMS machinery. `sp.AA.translate()` translates DNA/RNA→AA
(with a Rust codon LUT, `unknown="X"`/`"drop"` policies, `validate=True` fast-fail) and `sp.AA` includes
the stop codon `*`, but there is no codon-substitution enumeration, no CDS handling, no variant table.

**assay_mpra — no.** No regulatory-library design: no promoter/enhancer scaffolding, no barcoding,
no adapter/cloning handles, no motif perturbation.

**assay_insilico — partial.** SeqPro supplies the standard *primitives* used when probing
sequence-to-function models: one-hot encoding/decoding, tokenization, `k_shuffle` (k-let-preserving
dinucleotide shuffle — the canonical null/background for in-silico experiments, Rust-accelerated),
`random_seqs` (random backgrounds), `jitter` and `ReverseComplement` augmentations, and `Ragged`
batching. Explicitly framed as infrastructure for ML: *"heavily utilized by ... SeqExplainer and EUGENe."*
It does **not** design in-silico libraries, run models, do ISM, or score anything — that layer is
SeqExplainer/EUGENe, separate packages.

### Block C — genomics integration

**genome_coordinates — partial.** `python/seqpro/bed.py` provides `read` (dispatching to `_read_bed`,
`_read_narrowpeak`, `_read_broadpeak`), `sort`, `with_len(bed, length)` (resize intervals to a fixed
length), `to_pyr`/`from_pyr` (pyranges interop). `python/seqpro/_coords.py` provides `CoordSchema`,
`detect_schema`, `set_schema` for 0/1-based, open/closed interval conventions. So SeqPro ingests and
manipulates genomic intervals — but it **does not fetch sequence from a reference genome**; there is
no FASTA/2bit reader and no coordinate→sequence resolution (that lives in the sibling `SeqData`).

**transcript_models — partial.** `sp.gtf.scan(path)` parses GFF/GTF into a polars LazyFrame with
columns `seqname, source, feature, start, end, score, strand, frame, attribute`, and
`sp.gtf.attr("gene_id")` returns a polars expression regex-extracting one attribute as Utf8.
That is the entire GTF module (2 public functions, ~75 lines). There is **no transcript object, no
exon assembly, no CDS/UTR resolution, no gene model semantics** — a user gets a flat annotation table
and must build any transcript logic themselves.

**exon_intron_split_codons — no.** Nothing joins exons or tracks reading frame across exon
boundaries. `translate()` operates on a contiguous nucleotide buffer in stride-3 codons; the `frame`
GTF column is merely parsed (`pl.col("frame").fill_null(0)`) and never used.

**hgvs_input — no.** No HGVS parser or reference to HGVS anywhere in source, docs, or changelog.

**vcf_vep_output — no.** No VCF reader/writer, no VEP-compatible export. Package I/O is BED /
narrowPeak / broadPeak / GTF-GFF in, NumPy/polars out.

**consequence_annotation — no.** No variant-consequence logic. `AA.translate` can produce a stop
codon `*` for a given sequence, but nothing compares a variant to a reference or classifies a
consequence class.

### Block D — physical construction

**primer_design — no.** No primer, oligo, Tm, overhang, or assembly functionality.

**codon_optimization — no.** Codon tables exist strictly for translation (`AminoAlphabet(codons, amino_acids)`,
Rust codon LUT). There is no reverse-translation, no usage-frequency table, no optimizer.

**synthesis_constraints — no.** `gc_content`, `nucleotide_content`, `count_kmers_seq` are read-only
analyzers; `_cleaners.py` has `remove_N_seqs`, `remove_only_N_seqs`, `remove_whitespace`, `sanitize` —
input hygiene, not synthesis feasibility (no homopolymer, repeat, restriction-site, or
GC-window checks, no enforcement/repair loop).

### Block E — engineering

**interface — Python API only.** `import seqpro as sp`. `pyproject.toml` declares **no**
`[project.scripts]` — there is no CLI, no GUI, no web service. Docs site (MkDocs/zensical) at
https://ml4gland.github.io/SeqPro/. Requires Python ≥3.10.

**license — MIT.** GitHub API `license.spdx_id = "MIT"`; `pyproject.toml` `license = { file = "LICENSE" }`.

**maintained — yes, very actively.** Last commit `2026-07-27T19:02:07Z` ("bump: version 0.21.2 → 0.22.0");
release 0.22.0 tagged 2026-07-27, uploaded to PyPI 2026-07-27T19:07:05. 55 PyPI releases since
2023-03-26; 10 GitHub releases in the last ~2 months (0.16.0 on 2026-06-14 → 0.22.0 on 2026-07-27).
0 open issues. 14 stars, 1 fork.

## Documented examples / vignettes

SeqPro ships **no library-design vignettes** — the examples are API demos, not use cases.

1. **README / docs "Quick Start" snippet** (README.md; docs/index.md) — the canonical example: generate
   `sp.random_seqs(shape=(2,9), alphabet=sp.DNA, seed=1234)`, then `pad_seqs`, `DNA.ohe`/`decode_ohe`,
   `tokenize`/`decode_tokens`, `reverse_complement`, `k_shuffle(k=2)`, `gc_content`, `nucleotide_content`,
   `jitter(max_jitter=2)`, `bin_coverage(bin_width=128)`, `AA.translate(cds)`, `AA.ohe`.
2. **"Ragged Arrays" guide** (docs/ragged.md; https://ml4gland.github.io/SeqPro/ragged/) — worked example
   building `Ragged.from_lengths` over "per-position coverage for three intervals of different widths",
   ufunc arithmetic, `to_padded`, `to_packed`, record layouts.
3. **API Reference pages** (docs/api/{index,alphabets,bed,gtf,ragged,types}.md) — mkdocstrings-rendered
   signatures only; `docs/api/index.md` lists exactly 15 top-level members.
4. **`skills/seqpro/SKILL.md`** — an author-maintained agent-facing cheat sheet with a task→call quick
   reference table and Ragged construction/record-layout/hashing examples. Effectively the richest
   usage document in the repo.
5. **Test suite as usage examples** — `SKILL.md` points to `tests/` ("Tests as usage examples");
   ~30 test modules, dominated by `test_ragged*`, `test_ohe`, `test_tokenize`, `test_translate`,
   `tests/bed/*`.
6. **Sibling-package notebooks** — the org's `ML4GLand/tutorials` repo (last push 2025-02-06) contains
   notebooks only for `eugene/*` and `seqexplainer/*` (attribution_analysis, filter_interpretation,
   sequence_evolution); **there is no SeqPro tutorial notebook**. `ML4GLand/use_cases` last pushed 2024-01-03.

**Reproducibility note for PoolParty:** none of these are library-design use cases. There is nothing
here to "reproduce" in PoolParty other than trivially generating random sequences and computing GC —
which is the honest headline for the comparison table.

## Notable capabilities not in the row list

- **One-hot encoding / decoding / integer tokenization / padding**, vectorized over batches
  (`ohe`, `decode_ohe`, `tokenize`, `decode_tokens`, `pad_seqs`), Numba-accelerated.
- **k-let-preserving shuffling** (`k_shuffle`) implemented in Rust (`src/kshuffle.rs`) — the standard
  dinucleotide-shuffle null for regulatory-genomics ML.
- **`Ragged` variable-length array container** — Rust-native, Arrow-style flat data + offsets, exactly one
  ragged axis, NumPy ufunc dispatch, `to_padded`/`to_packed`/`to_numpy`/`concatenate`, multi-field
  **record layout** via `sp.rag.zip`, `np.memmap`-safe, awkward-array interop via `to_ak()`.
- **Parallel Rust string hashing** of ragged sequence batches: `rag.hash("sha256"|"md5"|"rapidhash")`.
- **Reverse complement that works on strings, `S1` bytes, and one-hot arrays alike**, and on coverage
  tracks (`transforms.ReverseComplement(*types=("dna","track"))`).
- **Coverage utilities**: `bin_coverage(bin_width, normalize)`, `normalize_coverage(method="CPM"|"CPKM", total_counts)`.
- **k-mer counting** (`count_kmers_seq`) and `length`.
- **Custom alphabets**: `NucleotideAlphabet(alphabet, complement)`, `AminoAlphabet(codons, amino_acids)`
  — user-definable, incl. non-standard codon tables.
- **Translation with explicit unknown-codon policy**: `AA.translate(..., validate=False, unknown="X"|"drop")`;
  `"drop"` returns a `Ragged` because lengths become variable.
- **Coordinate-schema handling** (`_coords.py`): `detect_schema`/`set_schema` for 0- vs 1-based and
  open vs closed intervals — a real correctness feature most tools leave implicit.
- **Peak-format readers**: BED, narrowPeak, broadPeak; pyranges and polars/narwhals interop.
- **torchvision-style transform objects** for ML dataloaders (`Sequential`, `Random`, `ReverseComplement`,
  `KShuffle`, `Jitter`, `Tokenize`).
- **Experimental XArray integration** (`seqpro.xr`) for SeqData interop.
- **Engineering rigor**: Rust/PyO3 + maturin, Numba, CodSpeed continuous benchmarking, pandera/narwhals
  schema validation, `py.typed`, pre-commit, conventional-commit changelog, dedicated benchmark harnesses.

## Stated limitations / scope disclaimers

- README ends with a bare **"# More to come!"** — the authors treat the surface as incomplete.
- `python/seqpro/experimental/*` is headed **"Unmaintained experimental module, not part of the public
  package surface."**, and `_visualizers.py` adds **"this module is currently broken at import time ...
  Left as-is per scope; needs porting if revived."**
- `seqpro.xr` is described as **"experimental integration with XArray"**.
- README frames SeqPro as a dependency of a suite: **"fully functional on its own but is heavily utilized
  by other packages including SeqData, MotifData, SeqExplainer, and EUGENe"** — i.e. sequence *acquisition*
  (SeqData), motifs (MotifData), and model probing (SeqExplainer) are explicitly out of scope.
- `SKILL.md`: **"Public API: see `python/seqpro/__init__.py` for the full export list. Re-read it before
  assuming a symbol exists."**
- 0.22.0 and 0.21.x carry **BREAKING CHANGE** entries for `Ragged.__getitem__` semantics — the API is
  still moving.

## Availability today

**Fully available and actively developed.** `pip install seqpro` → 0.22.0 (PyPI upload 2026-07-27,
requires Python ≥3.10; deps: numpy, numba, polars, pyranges, pandera, narwhals, awkward, pyarrow,
pandas, natsort, attrs). Repo public, not archived, MIT, 0 open issues, last commit 2026-07-27.
Docs site live at https://ml4gland.github.io/SeqPro/. No web service or CLI to be alive or dead.
Wheels are built via maturin (Rust extension) — a source install requires a Rust toolchain, but PyPI
wheels are published.

## Confidence notes

- Least certain: **`library_as_object`** and **`dag_chaining`** (both "partial"). The partials are
  generous readings — `Ragged`/record-layout and `transforms.Sequential` are real, first-class,
  composable objects — but they operate on already-materialized sequences for ML preprocessing, not on a
  design specification. A referee could argue "no" for both on the grounds that no *library* abstraction
  exists; a referee from ML4GLand could argue `Ragged` deserves credit. The evidence quotes support
  either read, so "partial" with the caveat is the defensible cell.
- **`assay_insilico` = partial** is the other judgment call: SeqPro supplies the canonical primitives
  (`k_shuffle`, `random_seqs`, `ohe`) but performs no design or model probing itself.
- **`lazy_evaluation` = no** despite `gtf.scan` returning a polars LazyFrame — the row asks about
  sequences, and sequence generation is unambiguously eager.
- Everything in Blocks B(dms/mpra), C(exon/hgvs/vcf/consequence), and D is a **high-confidence "no"**:
  verified against the complete `__all__` export list, the full repo file tree, the rendered docs
  navigation, and a keyword grep over all 55 releases of CHANGELOG.md (0 hits for
  mutagen/motif/barcode/primer/vcf/hgvs/library design).
- Source was read from `main` (0.22.0). I did not install or execute the package, per instructions.
