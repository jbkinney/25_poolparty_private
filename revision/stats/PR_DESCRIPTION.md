# Add `pool.stats()`; fix two silent defects in the export path

Answers Reviewer 2 item 1a of the BMC Bioinformatics revision, which asked for
"additional statistical readouts describing pool generation … how many unique
sequences vs duplicates … how unique are the sequences … how frequent are
homopolymers?"

Base `4f125f2`. 16 files, +2040 / −37.

---

## Read this first

The branch carries **three separable subjects at different levels of risk**, and
a reviewer should be able to take them separately.

| Subject | Commits | Risk |
|---|---|---|
| **`stats` readout** — a new read-only function | `437e45a`, most of `f6345d6` | **None.** Pure addition; cannot change any existing output |
| **`num_cycles` export fix** — changes shipped behaviour | `2ddedb1`, `cc82dc0`, part of `f6345d6` | Alters what `to_df`/`to_file` return on filtered pools |
| **FASTA null warning** — new warning on shipped code | `ae5e3d4` | Emits a `UserWarning` where none was emitted before |

If you want only the feature, drop `2ddedb1`, `cc82dc0` and `ae5e3d4`; the
`stats` code touches none of the same behaviour.

**The commits must be read in order.** `2ddedb1` argues for a rule that
`cc82dc0` overturns — see "How the export fix went wrong" below. Reading only
the first commit gives you the wrong rule.

---

## 1. The new readout

`pp.stats(pool)` and `pool.stats()` generate a library and describe it. Report
only: the pool is left exactly as it was, including the internal cursor that
decides which sequence `generate_library` returns next.

```python
>>> print(pool.stats())
pool.stats()  -  5 of 5 sequences in the design

Composition
  design size (num_states)             5
  generated                            5
  filtered out                         2
  unique sequences                     3
  duplicate sequences                  0   (0.0%)
  most-repeated sequence               1 copy

Length
  min / max                        8 / 8

GC content
  min / mean / max          0.000 / 0.167 / 0.500

Homopolymer runs
  longest run                          4
  sequences with a run > 6          0.0%

Repetitiveness (DUST)
  mean / max                 0.33 / 0.33

Pairwise distance (Hamming)
  exact, all 3 pairs
  min / mean / max          4 / 4.7 / 6
```

The return value is a `dict` of 25 keys, so numbers can be read out
individually, saved as JSON, or stacked with
`pd.DataFrame([a.stats(), b.stats()])`. All 25 are documented in `docs/pool.rst`
and cross-checked against the code.

### Why it is worth having: the paper's own libraries

| Library | What the readout found |
|---|---|
| **GB1 deep mutational scan** | 547,230 generated, **546,362 unique, 868 duplicates**. The wild-type ORF appears **130 times** where the design asks for 100 — the extra 30 from the `mutation_rate=0.1` random arm drawing zero mutations |
| **SpliceAI surrogate** | 200,000 generated, 98 duplicates, and **every sequence carries a homopolymer run longer than 6**, inherited from the 201-bp background at position 16. That design filters homopolymers out of its 9-mer *source*; filtering the parts does not constrain the assembled whole |
| **MPRA regulatory grammar** | 24,000 generated, 24,000 unique, 0 duplicates. Clean |

None of this was visible before. The GB1 and MPRA figures were reproduced
independently by a throwaway prototype before the feature was written, and the
shipped code matches them digit for digit.

### How much of the design gets measured

A pool records how to build a library rather than the library itself, and
duplicates cannot be counted in a recipe, so `stats` must generate first.

| Design | `pool.stats()` |
|---|---|
| Fixed size, ≤ 1,000,000 sequences | measures all of it |
| Fixed size, > 1,000,000 | raises, naming the size and how to proceed |
| No fixed size | raises: there is no "all of it" to measure |

An explicit `num_seqs=` is always honoured and never capped. The 1,000,000 limit
reuses `Operation.max_num_sequential_states` rather than introducing a second
threshold.

**Expect a question about the third row.** Many operations default to
`mode='random'` — `mutagenize`, `shuffle_seq`, `from_iupac`, `get_kmers`,
`from_seqs`, the scans — and without `num_states` such an operation draws a fresh
sequence per row, so the design genuinely has no total. `pool.num_states` reports
1 for a design that produces hundreds of distinct sequences. Refusing to guess
was a deliberate choice; `docs/pool.rst` names the affected operations so the
first failure is self-explanatory.

### Cost

Everything except the pairwise distances is exact and linear. Measured against
generation time on the paper's libraries: **7% (MPRA), 19% (SpliceAI), 5%
(GB1)**. Comparing every pair is quadratic, so above `max_hamming_seqs` (default
2,000) a seeded subsample is compared and `hamming_exact` reports `False`.

Checked against an exact all-pairs run on the 24,000-sequence MPRA library: a
2,000-sequence subsample estimates the **mean** to within 0.1%, but gives a
minimum of 1, 2 or 3 against a true 1. The docs, the docstring and the printed
report all state that a sampled minimum is an upper bound and the maximum a lower
bound. The question people actually have — *are any two sequences identical?* —
is answered exactly by `num_duplicate_seqs`, in linear time.

### Deliberately not included

- **No deduplication and no filtering.** Five constraint filters already act on
  GC, homopolymers, repetitiveness and restriction sites, so `stats` measures and
  `filter_*` fixes; a filtering option here would mean two implementations of the
  same predicates. Deduplication cannot be an operation at all — `op.compute`
  sees only its parents and an RNG, so no step can know what has already been
  emitted.
- **No melting temperature, DNA folding, or reference-genome matching.** New
  dependencies or external databases, and primer-only in every comparable tool.
- **No edit distance.** O(L²) per pair, roughly 200× Hamming's cost at 200 bp.
- **Protein pools unsupported.** Five of the statistics are DNA-specific;
  `ProteinPool` has no export methods either. `pp.stats(protein_pool)` raises
  `TypeError`.

The boundary was set by surveying what eight comparable tools report. The
included set covers every statistic with two or more tools behind it — and
notably, **no surveyed tool reports a duplicate count or a library-wide distance
distribution at all**.

---

## 2. `num_cycles` no longer re-emits sequences from a filtered pool

A filter replaces a rejected sequence with `NullSeq` rather than removing it, so
one cycle through a filtered pool yields fewer rows than `num_states`. All four
export loops treated `num_cycles` as a *row* target and made up the shortfall by
asking `generate_library` for more, which restarted the traversal and returned
survivors a second time. On a five-state pool whose filter rejects two:

```
to_df(num_cycles=1)   before: 5 rows holding 3 distinct sequences
                       after:  3 rows
```

Nothing indicated that two were copies. `generate_library` had already warned
that only three valid sequences existed; the loop discarded the warning and asked
again.

`num_cycles` now counts states traversed rather than rows collected, **for a
design whose states determine its sequences**. Cycling itself is unchanged:
`num_cycles=3` on a filtered pool still returns each survivor three times, and an
unfiltered pool still returns `num_cycles × num_states`. `num_seqs` is untouched
— asking for a count legitimately means "keep drawing".

### How the export fix went wrong, and why the history keeps it

`2ddedb1` applied that rule unconditionally. It is only valid when the state
determines the sequence: a random operation built without `num_states` is seeded
from the row counter, so for those designs the top-up had been producing **fresh**
sequences, and removing it lost rows that were never duplicates.

```
from_seqs(S5).filter_gc(...)          num_cycles=10   base 10 rows -> 6
5-state source -> mutagenize(random)  num_cycles=1    base  5 rows -> 3
                                                      (5 distinct: nothing
                                                      was being re-emitted)
seeded random pool, num_cycles=4      chunk 1/2/3/5/1000
                                                      base 4/4/4/4/4 -> 2/2/2/3/3
```

`cc82dc0` corrects the predicate and shares it with `stats` as
`_draws_fresh_sequences`, so the export path and the readout can no longer
disagree about what a design's size means. Base parity is restored on all three
cases and the original defect stays fixed at every chunk size.

The history is kept unsquashed on purpose, so the wrong turn is on the record.
That is why the commits have to be read in order.

---

## 3. FASTA now says when it drops a filtered-out sequence

A FASTA record needs a sequence, so a rejected one has no representation and was
dropped even with `discard_null_seqs=False`. CSV, TSV and JSONL write a blank and
keep the row, so the same library exported two ways gave files of different
lengths and `to_file` returned different counts, with nothing to explain it.

Each format was individually right, so the change makes the loss visible rather
than altering it: one `UserWarning` per call naming how many records were
omitted, plus a sentence on `to_file`'s `discard_null_seqs`. This follows the
`write_style` warning a few lines above it in the same function. It fires only
when something was actually lost — unlike `write_style`, the loss is
data-dependent, so warning on the flag alone would fire on every unfiltered FASTA
export.

Pre-existing behaviour, unrelated to item 1a, kept as its own commit so it can be
dropped independently.

---

## 4. Documentation corrections

Item 1a exposed three statements that were wrong, all now fixed:

- `docs/pool.rst` and `docs/operations/library_size.rst` described `num_states`
  as the number of *distinct sequences* a pool produces. It is a state count, and
  neither an upper nor a lower bound on the sequence count: an open-ended design
  reports 1 while producing hundreds.
- `library_size.rst` listed `filter` as reducing `num_states`. It does not.
- `calc_dust` cited Hancock & Armstrong (1994), which describes a different
  algorithm (SIMPLE). The implementation is the standard DUST score; it now cites
  Morgulis et al. (2006).

---

## How this was checked

- **Full suite: 3,076 passed, 14 xfailed** (base: 2,974). Sphinx builds with
  exit 0 and no warning or error located in any file this branch touches. `ruff
  check` and `ruff format --check` clean on every file added or edited; the
  repository's pre-existing findings are unchanged.
- **Three independent review passes** — adversarial correctness, API design, and
  test quality — each verifying claims by running code rather than reading it.
  Every finding was reproduced before being acted on; one was rejected as
  pre-existing rather than introduced here.
- **Mutation-checked.** Nine plausible mutations that the first version of the
  suite accepted are now caught, among them pinning the export predicate true
  (8 failures), dropping the reproducibility pinning (3), discarding the `seed`
  argument (3), ignoring `num_cycles` (2), and inverting `longest_homopolymer`
  (5). The blocked pairwise-comparison path is exercised explicitly, since no
  default call reaches it.
- **Reproducibility.** `stats` pins both the start state and the seed, so a
  report depends on the design alone. Verified: the same design reports
  identically whether or not sequences were generated from the pool beforehand,
  and the pool's cursor and seed are exactly as they were afterwards.
- **The printed example in `docs/pool.rst` is verified character-for-character**
  against what the code prints.

## Reviewer notes

- `docs/api.rst` gains one `autofunction` entry. A separate piece of work on
  Reviewer 2 item 2c also edits that file to document the five existing
  `filter_*` methods; the two additions are in different sections and should
  merge cleanly, but they are worth landing in a known order.
- `docs/operations/library_size.rst` is edited here, which also overlaps that
  work.
- The `!` in `cc82dc0`'s type is relative to `2ddedb1`, not to base. Nothing on
  this branch is breaking with respect to `origin/main`.
