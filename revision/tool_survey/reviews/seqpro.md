# Adversarial review — SeqPro extraction

**Reviewed:** `extractions/seqpro.md` + the structured extraction for slug `seqpro`
**Review date:** 2026-08-10
**Reviewer method:** independent re-fetch of `ML4GLand/SeqPro@main` (GitHub API metadata, recursive git tree,
raw source for every Python module, `CHANGELOG.md`, `pyproject.toml`, `docs/`, `skills/seqpro/SKILL.md`),
PyPI JSON, and the `ML4GLand` org + `ML4GLand/tutorials` trees. No install, no clone, no execution.

## Headline

**The extraction survives falsification.** Every "no" I tried to break, I could not break. I independently
re-derived `__all__`, the recursive file tree, the CHANGELOG keyword greps, and the source of every module
the extraction cites, and the verbatim quotes it uses are accurate. All 24 capability values are defensible;
none is an understatement of a real design capability.

What I *did* find are **five small factual defects in the supporting evidence text** and **five genuinely
missed capabilities**. None flips a cell, but two of the defects are the kind of thing an ML4GLand referee
(Laub/Klie are the authors, both at UCSD, both active on the repo *this month*) would notice and use to
impeach the whole row. Fix them before the table ships.

---

## Independently verified (my own retrieval, not the extractor's quotes)

- **`__all__` is exactly as claimed.** Fetched `python/seqpro/__init__.py`. The 28 exported names match the
  extraction character-for-character. No mutagenesis, motif, barcode, primer, variant, naming, or plotting
  symbol.
- **Recursive file tree (`truncated: false`)** — 160-odd paths. Python source is 25 files. There is no
  `fasta`, `vcf`, `variant`, `motif`, `primer`, `design`, or `library` module. Confirmed.
- **`python/seqpro/experimental/` has NO `__init__.py`** (HTTP 404 on the raw URL; absent from the tree).
  This is *stronger* than the extraction's claim that the module "is not imported by `__init__.py`" — it is
  not even an importable package. `seqpro.experimental` raises `ModuleNotFoundError`. Use this; it makes
  the `barcode_generation` and `design_visualization` "no"s unassailable.
- **`_experimental.py` header** verbatim: `# Unmaintained experimental module, not part of the public
  package surface.` and `edit_distance` really does `assert len(seq1) == len(seq2), "Both sequences must be
  of same length."` — it is Hamming, as claimed.
- **`_visualizers.py` header** verbatim, including `NOTE: this module is currently broken at import time`.
  Confirmed: it imports `gc_content_seqs` / `nucleotide_content_seqs`, which do not exist in `_analyzers.py`.
- **`random_seqs`** body is literally `return seed.choice(alphabet.array, size=shape)`. Eager. Confirmed.
- **`gtf.py`** is two public functions, `scan` and `attr`, 1,813 bytes total. `pl.col("frame").fill_null(0)`
  is present and `frame` is never consumed anywhere in the package. Confirmed.
- **`transforms/augmentation.py`** `Sequential.__call__` is `for t in self.transforms: x = t(x); return x`.
  Confirmed verbatim.
- **CHANGELOG grep** over all 608 lines / 55 releases for
  `mutagen|motif|barcode|primer|variant|vcf|hgvs|exon|intron|fasta|2bit|genome|provenance|plot|visual|design`
  → the only hits are `codon` (6×, all translate-kernel perf work), one `bed.sort() temp column name
  collision`, one `rename bed functions`, one `restore designed contract`. **0 substantive hits.** Confirmed.
- **Source-wide grep** across every fetched module for
  `fasta|2bit|bigwig|motif|mutat|barcode|primer|restriction|homopolymer|provenance|hgvs|vcf|consequence|codon.?optim|reverse.?transl|oligo|assembl`
  → 2 hits, both the English word "mutate" in `rag/_ops.py` docstrings about in-place buffer mutation.
- **`pyproject.toml`** has no `[project.scripts]`. No CLI. Confirmed.
- **Repo/PyPI metadata**: MIT (`spdx_id: MIT`), not archived, 14 stars, 1 fork, 0 open issues, last commit
  `2026-07-27T19:02:07Z`, `pushed_at 2026-07-27T19:07:56Z`, 55 PyPI releases, first `0.0.1` at
  `2023-03-26T17:33:00`, latest `0.22.0` at `2026-07-27T19:07:05`. 10 GitHub releases 0.16.0→0.22.0 in the
  window claimed. All confirmed.
- **`ML4GLand/tutorials`** tree re-fetched: notebooks under `eugene/*` and `seqexplainer/*` only;
  `motifdata/TODO.txt`; **no `seqpro/` directory**. The negative finding is correct.

---

## Defects in the evidence text (fix before referees see it)

### D1. "BED reading is polars-backed [lazy]" — WRONG (in `lazy_evaluation` evidence)

The extraction writes: *"The only laziness in the package is annotation-table I/O: `sp.gtf.scan(path) ->
pl.LazyFrame` (polars scan_csv over the GTF) and polars-backed BED reading. That is lazy table scanning."*

`bed.read` → `_read_bed` / `_read_narrowpeak` / `_read_broadpeak` all call **`pl.read_csv`**, which is
**eager**. Only `gtf.scan` is lazy. Delete "and polars-backed BED reading" from that sentence.

### D2. "The only laziness in the package is annotation-table I/O" — FALSE, and referee-baitable

`python/seqpro/xr/__init__.py` implements `ohe`, `bin_coverage`, and `translate` via
`xr.apply_ufunc(..., dask="parallelized")`. On a dask-backed xarray input these are **lazily evaluated and
chunk-parallel** — real out-of-core laziness over sequence data, not table I/O. An ML4GLand referee will know
this module exists (it is their SeqData bridge) and the universal quantifier ("the only laziness") is a free
hit for them.

**The `lazy_evaluation = no` verdict still stands** — the row asks about lazy *generation of a designed
library*, and SeqPro generates nothing lazily (`random_seqs` is eager; there is no design step to defer).
But the evidence sentence must be rewritten to something like: *"No sequence-producing function is lazy —
`random_seqs` materializes immediately. Deferred execution exists only for I/O and array transforms
(`gtf.scan` → `pl.LazyFrame`; the experimental `seqpro.xr` ops run via `xr.apply_ufunc(dask='parallelized')`),
never for library construction."* That version is unattackable.

### D3. `count_kmers_seq` and `normalize_coverage` are NOT public API

Both are listed in `additional_capabilities` and `count_kmers_seq` is also cited in the
`synthesis_constraints` evidence as an "analyzer". Neither is in `__all__`, neither is imported by
`__init__.py`. `sp.count_kmers_seq` and `sp.normalize_coverage` raise `AttributeError`; they are reachable
only as `seqpro._analyzers.count_kmers_seq` / `seqpro._modifiers.normalize_coverage`. Worse, the *batched*
k-mer counter `_analyzers._count_kmers` is `raise NotImplementedError` on its first executable line, with a
`# TODO: non-trivial to parallelize` comment above it. Qualify both, or drop them.

### D4. `remove_whitespace` is a no-op stub

Cited in `synthesis_constraints` and in `additional_capabilities` as one of the "input hygiene cleaners".
Its full body is `pass`. It does nothing. (Doesn't change the "no" — if anything it reinforces it — but a
referee reading the cleaner list will spot it.)

### D5. Member-count arithmetic

`docs/api/index.md` lists **14** members, not 15 (I counted the rendered file: `bin_coverage, cast_seqs,
decode_ohe, decode_tokens, gc_content, jitter, k_shuffle, length, nucleotide_content, ohe, pad_seqs,
random_seqs, reverse_complement, tokenize`). The extraction says "exactly 15" twice — in the
`automatic_naming` evidence and in `documented_examples` — while its own enumeration lists 14. Trivial, but
it is a checkable number in an evidence field.

Minor, sub-defect: the `vcf_vep_output` dependency list omits `polars-config-meta[polars]>=0.3.2` (used to
attach `coordinate_system_zero_based` metadata to BED frames). Doesn't affect the claim — that package is
not a variant-format library either.

---

## Capabilities the extraction missed entirely

### M1. `seqpro.transforms.TMM` — Trimmed Mean of M-values normalization (**biggest miss**)

`python/seqpro/transforms/tmm.py` implements a full edgeR-style TMM estimator: a `TMM` class with
`fit(counts)` / `transform(counts, library_size=1e6)`, reference-sample selection by quantile, log-ratio and
absolute-expression trimming, inverse-variance weighting, and a `@nb.njit(parallel=True)` `_tmm_helper`
kernel. It is **exported**: `transforms/__all__ = ["TMM", "Jitter", "KShuffle", "Random",
"ReverseComplement", "Sequential"]`.

The extraction never mentions it, and its `dag_chaining` evidence enumerates the transform objects as
"ReverseComplement, KShuffle, Jitter, Tokenize" — which is the *wrong* list (see M4). This is an entire
assay-data-normalization capability absent from the memo. It does not create any Block A–D "yes" (TMM
normalizes count matrices, it does not design sequences), but it belongs in `additional_capabilities` and it
is exactly the kind of omission an author-referee reads as "they didn't actually look at our package".

### M2. `AA.translate(..., truncate_stop=True)` — stop-codon truncation

Truncates each output at the first stop codon, inclusive; Ragged input only; applied after drop-compaction.
There is also `_translate_stop_ends` (a dedicated Rust kernel) and an OHE path that locates the `*` column.
This is more coding-sequence awareness than the extraction's `assay_dms` / `consequence_annotation` /
`exon_intron_split_codons` evidence implies (those say only that `AA` "includes the stop codon `*`").

**Verdicts do not change** — truncating at a stop is not variant consequence classification, not exon
assembly, and not mutant-library construction. But if a referee says "SeqPro handles stop codons", the memo
currently has no answer prepared. Add: *"`translate(truncate_stop=True)` truncates at the first stop codon;
there is still no reference allele, no variant representation, and no consequence class."*

### M3. `sp.rag.reverse_complement` — Ragged-native reverse complement

`rag/_ops.py` exports `reverse_complement` alongside `concatenate`, `hash`, `to_packed`, `to_padded`. The
extraction's `additional_capabilities` covers RC on `str`/`S1`/OHE and on coverage tracks, but not on
variable-length Ragged batches.

### M4. `Tokenize` is **not** exported from `seqpro.transforms`

`dag_chaining` evidence and `additional_capabilities` both list `Tokenize` among the "torchvision-style
transform objects". It is defined in `augmentation.py` but omitted from `transforms/__all__` (which instead
contains `TMM`, unmentioned). `sp.transforms.Tokenize` works via attribute access but is not part of the
declared surface. Swap `Tokenize` → `TMM` in that list.

### M5. `bed.with_len` is peak-summit aware

The extraction says `with_len(bed, length)` "resize[s] intervals to a fixed length". It actually centers on
the **narrowPeak summit** when a `peak` column is present (`2 * (chromStart + peak)`), falling back to the
interval midpoint otherwise. A genuine correctness nicety in the same class as the coordinate-schema
handling the extraction (rightly) praises.

Also worth one line: `Ragged.to_chars()` / `to_strings()` (opaque-string ↔ S1-char conversion), the
`NucleotideAlphabet.tokenize` / `.decode_tokens` methods, and dev-deps `biopython` + `hypothesis[numpy]`
(property-based testing) as further engineering-rigor evidence.

---

## Cell-by-cell

| key | value | verdict | note |
|---|---|---|---|
| library_as_object | partial | supported | Generous but the safe direction. `Ragged` + `rag.zip` record layout are genuinely first-class and composable, so "no" would be attackable by an author-referee. Errs toward *crediting* SeqPro, which is the correct error to make here. |
| dag_chaining | partial | supported | `Sequential` / `Random` verified verbatim; both exported. Fix the operand list (M4). Linear composition, no graph — "partial" is right. |
| lazy_evaluation | no | supported | Value correct; **evidence must be rewritten** (D1, D2). `seqpro.xr` dask laziness is real and the "only laziness" sentence is false. |
| mixed_mutagenesis_one_pool | no | supported | Re-verified via `__all__`, full tree, source grep, CHANGELOG grep. Nothing mutagenic exists. |
| combinatorial_motif_place | no | supported | Zero `motif` hits anywhere in source or changelog. Motifs live in the separate, 2023-dormant MotifData. |
| barcode_generation | no | supported | Strengthen with the `experimental/` has-no-`__init__.py` finding — `seqpro.experimental` is not importable at all. |
| per_sequence_provenance | no | supported | `rag.zip` / `rag.hash` are user-driven mechanisms; nothing in the package writes provenance. Correctly characterized. |
| automatic_naming | no | supported | Note for robustness: `bed.py`'s `BEDSchema` *reads* a `name` column and record-Ragged has *field* names — neither generates per-sequence IDs. Fix the "15 members" count (D5). |
| design_visualization | no | supported | Verbatim headers confirmed; module is not importable (no `__init__.py`) and imports two functions that no longer exist. |
| assay_dms | no | supported | Add the `truncate_stop` pre-rebuttal (M2). |
| assay_mpra | no | supported | Nothing regulatory-design anywhere. |
| assay_insilico | partial | supported | Judgment call, correctly flagged. `k_shuffle` (Rust, `src/kshuffle.rs`, CodSpeed-benchmarked) really is the canonical null-background primitive; SeqPro supplies primitives only. |
| genome_coordinates | partial | supported | Verified `read`/`sort`/`with_len`/`to_pyr`/`from_pyr` + `CoordSchema`/`detect_schema`/`set_schema` with 5 built-in schemas (bed/pb/pr/gtf/gff). No FASTA/2bit anywhere — grep-confirmed. **Caveat:** if the row means "accepts genomic intervals" rather than "resolves coordinates to sequence", SeqPro is a clean "yes". Make sure the row legend says coordinate→sequence resolution, or this cell reads as an understatement. |
| transcript_models | partial | supported | `gtf.py` really is 2 functions / ~75 lines. "Partial" is generous (a flat annotation table is not a transcript model); the generosity is in the safe direction. |
| exon_intron_split_codons | no | supported | Confirmed `frame` is parsed, defaulted, never consumed. `truncate_stop` (M2) is single-buffer and does not cross exon boundaries. |
| hgvs_input | no | supported | 0 hits, four independent ways. |
| vcf_vep_output | no | supported | Confirmed; dep list is missing `polars-config-meta` (also not a variant library). |
| consequence_annotation | no | supported | Confirmed. No variant representation exists at all. |
| primer_design | no | supported | 0 hits for primer/oligo/Tm/restriction/assembly. |
| codon_optimization | no | supported | `AminoAlphabet` holds `codon_to_aa` and a 64-entry LUT; there is no reverse map, no usage table, no optimizer. All `codon` CHANGELOG hits are translate perf work — re-verified. |
| synthesis_constraints | no | supported | Confirmed. Trim the cited "analyzers" per D3/D4 (`count_kmers_seq` is private; `remove_whitespace` is `pass`). |
| interface | yes | supported | No `[project.scripts]` in `pyproject.toml`; docs homepage `https://ml4gland.github.io/SeqPro/` confirmed as the repo's declared homepage. Python-library-only is correct. |
| license | yes | supported | GitHub `license.spdx_id = "MIT"`, `LICENSE` at root. |
| maintained | yes | supported | Among the most actively maintained tools in this survey. All dates independently re-confirmed. |

---

## Framing advice for the referee response

The extraction's own `confidence_notes` land on the right posture and should be kept: *SeqPro is an
excellently engineered ML-preprocessing library that is not in the library-design problem space, and its own
README and SKILL.md say so.* Two additions:

1. **Lead with the authors' own scope statement, not with our "no"s.** `SKILL.md`'s "When to use" list
   (encoding/decoding, augmentation, stats, ragged batches, BED/GTF I/O) and the README's "fully functional
   on its own but is heavily utilized by SeqData, MotifData, SeqExplainer, and EUGENe" do the work for us. A
   row of "no"s reads as an attack; a scope quote reads as a fair reading.
2. **Praise the engineering explicitly.** Rust/PyO3 + maturin, Numba, CodSpeed CI benchmarking, pandera
   schema validation, `py.typed`, hypothesis property tests, ~50 design/spec docs under
   `docs/superpowers/`, and 10 releases in six weeks. The comparison must not read as though SeqPro is
   deficient — it is excellent at a different job. Note also `transforms.TMM` (M1), which shows the package
   reaches into assay-count normalization; omitting it invites "you didn't read our code".
