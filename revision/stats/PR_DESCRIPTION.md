# Add `pool.stats()`; fix two silent defects in the export path

Answers Reviewer 2 item 1a, which asked how many sequences in a pool are unique
versus duplicates, how far apart they are, and how often homopolymers appear.

Base `4f125f2`. 16 files, +2075 / −40, in 7 commits.

## What you are approving

Four separable subjects, at different levels of risk:

| Subject | Commits | Risk |
|---|---|---|
| **`stats` readout** — a new read-only function | `437e45a`, most of `f6345d6`, most of `424a1ae` | **None.** Pure addition; cannot change existing output |
| **`num_cycles` export fix** | `2ddedb1`, `cc82dc0`, part of `f6345d6` | Changes what `to_df`/`to_file` return on filtered pools |
| **FASTA null warning** | `ae5e3d4` | Emits a `UserWarning` where none was emitted before |
| **Documentation** | `d5106ba`, most of `424a1ae` | Docs only, including three corrections to statements that were already wrong before this branch |

For the feature alone, drop `2ddedb1`, `cc82dc0` and `ae5e3d4`; nothing in the
`stats` code depends on them.

**Read the commits in order, and read the last one.** Two commits correct an
earlier commit on this branch: `cc82dc0` overturns a rule `2ddedb1` argues for
(§2), and `424a1ae` fixes a defect `f6345d6` introduced (§5). Reading either
earlier commit alone gives you a wrong rule and code that is not what shipped.
The history is unsquashed on purpose.

## Where to start

`src/poolparty/stats.py` and `src/poolparty/utils/stats_utils.py` are the new
logic (~510 lines). `pool_mixins/export_mixin.py` is the behaviour change (~90
lines). The remaining ~1,400 lines are tests and docs. If you read only one
thing, read `_generate_rows` in `stats.py` — it decides how much of a design
gets measured, and both defects fixed within the branch were in that area.

## 1. The new readout

`pp.stats(pool)` and `pool.stats()` generate a library and describe it. Report
only: the pool is left exactly as it was, including the cursor that decides which
sequence `generate_library` returns next.

```
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

The return value is a plain `dict`, so numbers can be read out individually,
saved as JSON, or stacked with `pd.DataFrame([a.stats(), b.stats()])`:

`num_states`, `open_ended`, `num_generated_seqs`, `frac_design_covered`,
`num_filtered_out_seqs`, `num_valid_seqs`, `num_unique_seqs`,
`num_duplicate_seqs`, `frac_duplicate_seqs`, `max_seq_copies`, `length_min`,
`length_max`, `gc_min`, `gc_mean`, `gc_max`, `longest_homopolymer`,
`frac_seqs_with_long_homopolymer`, `dust_mean`, `dust_max`,
`frac_seqs_with_restriction_site`, `hamming_exact`, `hamming_seqs_compared`,
`hamming_min`, `hamming_mean`, `hamming_max`.

### What it finds in the paper's own libraries

| Library | Result |
|---|---|
| GB1 deep mutational scan | 547,230 generated, 546,362 unique, **868 duplicates**. The wild-type ORF appears **130 times** where the design asks for 100 — the extra 30 from the random arm drawing zero mutations |
| SpliceAI surrogate | 98 duplicates in 200,000, and **every sequence** carries a homopolymer run longer than 6, inherited from the fixed background. That design filters homopolymers out of its 9-mer source |
| MPRA regulatory grammar | 24,000 generated, 24,000 unique, 0 duplicates |

None of this was visible before.

### How much of the design gets measured

A pool records how to build a library rather than the library itself, so `stats`
must generate before it can measure anything.

| Design | `pool.stats()` |
|---|---|
| Fixed size, ≤ 1,000,000 | measures all of it |
| Fixed size, > 1,000,000 | raises, naming the size and how to proceed |
| No fixed size | raises: there is no "all of it" to measure |

An explicit `num_seqs=` is always honoured and never capped. The limit reuses
`Operation.max_num_sequential_states` rather than adding a second threshold.

**Expect a question about the third row.** Many operations default to
`mode='random'` — `mutagenize`, `shuffle_seq`, `from_iupac`, `get_kmers`,
`from_seqs`, the scans — and without `num_states` they draw a fresh sequence per
row, so `num_states` reports 1 for a design producing hundreds. Refusing to guess
was deliberate; `docs/pool.rst` names the affected operations.

Everything except the pairwise distances is exact, and costs 5–19% on top of
generation. Distances are subsampled above `max_hamming_seqs` (default 2,000),
with `hamming_exact` reporting `False`; a sampled minimum is documented as an
upper bound, since measurement showed it varies (1, 2 or 3 against a true 1)
while the mean is accurate to 0.1%.

### Deliberately excluded

No deduplication or filtering — the five existing `filter_*` methods already act
on GC, homopolymers, repetitiveness and restriction sites, so `stats` measures and
`filter_*` fixes. Deduplication could not be an operation anyway: `op.compute`
sees only its parents and an RNG, so no step can know what was already emitted.
No melting temperature or folding (new dependencies, primer-only in every
comparable tool), no reference-genome matching (external database), no edit
distance (~200× Hamming's cost). No protein pools — five statistics are
DNA-specific.

## 2. `num_cycles` no longer re-emits sequences from a filtered pool

A filter replaces a rejected sequence with `NullSeq` rather than removing it, so
one cycle yields fewer rows than `num_states`. The export loops treated
`num_cycles` as a row target and made up the shortfall, restarting the traversal
and returning survivors twice. On a five-state pool rejecting two,
`to_df(num_cycles=1)` returned **5 rows holding 3 distinct sequences**, silently.

`num_cycles` now counts states traversed, for designs whose states determine
their sequences. Cycling is unchanged (`num_cycles=3` still returns each survivor
three times); `num_seqs` is untouched.

`2ddedb1` applied that rule unconditionally, which lost rows on random-mode
designs where the top-up had been producing *fresh* sequences. `cc82dc0` corrects
the predicate and shares it with `stats`, so the export path and the readout
cannot disagree about what a design's size means.

## 3. FASTA now says when it drops a filtered-out sequence

A FASTA record needs a sequence, so a rejected one was dropped even with
`discard_null_seqs=False`, while CSV/TSV/JSONL kept the row — same library, files
of different lengths, no explanation. One `UserWarning` per call naming the count,
following the `write_style` warning in the same function. Pre-existing behaviour;
own commit so it can be dropped.

## 4. Documentation corrections

Four statements in the existing docs were wrong before this branch:

- `num_states` was described as the number of *distinct sequences* a pool
  produces. It is a state count, and neither an upper nor a lower bound.
- `library_size.rst` listed `filter` as reducing `num_states`. It does not.
- The same table listed `sample` as reducing it. Sampling 999 states from 256
  raises the count, so that row now reads "sets a new size".
- `calc_dust` cited Hancock & Armstrong (1994), which describes a different
  algorithm. Now cites Morgulis et al. (2006).

Two statements added earlier on this branch were also wrong and are corrected in
`424a1ae`: `num_states` called a "floor" for an open-ended design, in the same
file that says two paragraphs earlier it is neither bound — and a floor is a
lower bound; and gap characters said to "count as characters", which holds for
length, homopolymer runs and DUST but not for GC, since `calc_gc` ignores
anything outside ACGT.

## Verification

- **3,077 tests pass, 14 xfailed** (base 2,974).
- Sphinx exits 0 with no warning located in any file this branch touches.
- `ruff check` and `format --check` clean on every file added or edited.
- Four independent review passes — correctness, API design, tests, and a
  documentation pass whose voice assessment was written and sealed before the
  diff was visible. Each verified claims by running code rather than reading it.
  Two findings were rejected as pre-existing rather than introduced here, and
  three were declined as adding words readers do not need.
- Nine mutations the first version of the suite accepted are now caught,
  including pinning the export predicate true (8 failures) and discarding the
  `seed` argument (3).
- The printed example above is verified character-for-character against output,
  and all 25 documented keys are cross-checked against what the code returns.

## 5. Two defects were found and fixed within the branch

Both are already fixed in the diff you are reviewing; they are listed because
the commits are unsquashed and you will meet them in the history.

`2ddedb1` scoped the `num_cycles` change too broadly. A random operation built
without `num_states` is seeded from the row counter, so for those designs the
top-up it removed had been producing fresh sequences rather than repeats, and
rows were lost — silently, because the change also suppressed the warning that
would have said so. `cc82dc0` corrects the predicate and shares it with `stats`,
so the export path and the readout cannot disagree about what a design's size
means.

`f6345d6` added the restore that makes `stats()` leave the pool untouched, but
read the saved cursor with a `None` default and put it back only when it was not
`None`. A pool that had never generated therefore kept a cursor it did not
start with, and its first real `generate_library` returned the wrong window.
`424a1ae` distinguishes absent from zero with a sentinel. The test missed it
because it generated before calling `stats`, so it only exercised the
restore-a-value path; the added test fails without the fix.

## Reviewer notes

- `docs/api.rst` and `docs/operations/library_size.rst` are also touched by the
  separate item-2c work. Different sections; worth landing in a known order.
- The `!` in `cc82dc0` is relative to `2ddedb1`, not to base. Nothing here is
  breaking with respect to `origin/main`.
