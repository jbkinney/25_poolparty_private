# Constraint-filter API and documentation plan

Implementation plan for making PoolParty's ready-made sequence filters complete,
auditable, and discoverable.

**Status:** package implementation, pool-stats reconciliation, tests, and
documentation verification complete; PR #21 merged into `main` with merge
commit `d4618fa` on 2026-08-22. Manuscript follow-up remains pending.

**Package worktree:** the temporary
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-constraint-filters` linked
worktree was removed after merge verification. The primary package worktree is
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter`, on `main`.

**Package branch:** `docs/constraint-filters`, based on local `origin/main` at
`4f125f2` when the worktree was created. Its local and remote references are
retained at `eee404f` after the linked worktree's removal.

**Environment constraint:** use system Python for code, examples, and tests. Use
the existing ``poolparty_dev`` Conda environment only for the fully configured
Sphinx build. Do not create or use a virtual environment, and do not install
anything.

---

## 1. Problem

PoolParty already provides five convenience methods that wrap the generic
`filter` operation:

- `filter_gc`
- `filter_homopolymer`
- `filter_complexity`
- `filter_dust`
- `filter_restriction_sites`

They are tested but difficult to discover. They do not appear in the operations
guide or API reference, and unlike generic `filter`, they do not expose its
optional `cards=` argument. The current complexity documentation also promises a
0--1 score even though accepted IUPAC input can produce a value above 1, and the
DUST description overstates how closely the implementation matches NCBI
DustMasker.

The existing `filter` guide already owns predicate behavior, `NullSeq`, design
cards, and final output semantics. The work will improve that page rather than
create a competing constraint-filter page.

---

## 2. Settled decisions

| # | Decision | Status |
|---|---|---|
| 1 | Keep one generic `FilterOp`; the named methods remain thin convenience wrappers around it. | approved |
| 2 | Add an inclusive `filter_length(min_length=None, max_length=None, ...)` convenience method. | approved |
| 3 | Length filtering is sequence-generic, so place it on the existing generic Pool operation path rather than in the DNA-only `FilterMixin`. Do not create a new mixin or framework for it. | approved |
| 4 | Add optional `cards=` passthrough to `filter_length` and all five existing named wrappers. | approved |
| 5 | Cards expose only the existing Boolean `passed` key. They do not record calculated GC, length, DUST, or complexity values. Use `score` when measured values are needed. | approved |
| 6 | Preserve `cards=None` as the default, so existing output and behavior remain unchanged. | approved |
| 7 | Keep both `filter_complexity` and `filter_dust`; they overlap conceptually but measure different forms of repetition and can disagree. | approved |
| 8 | Fix and precisely define PoolParty's short-k linguistic-complexity calculation. Its result must honor the documented 0--1 range, and `k_range` must be validated. | approved |
| 9 | Retain the current DUST calculation, but describe it as a DUST-style whole-sequence triplet-repetition score rather than the complete NCBI windowed masking algorithm. | approved |
| 10 | Fold the named-filter guidance into `docs/operations/filter.rst`; do not create `constraint_filters.rst`. | approved |
| 11 | Treat the six named methods as convenience methods that construct a `filter` operation, not as six additional operations. | approved |
| 12 | Keep one `filter` entry in the State Operations sidebar and toctree. Use anchors and API links for individual convenience methods. | approved |
| 13 | Keep the documentation lean: shared behavior and parameters are explained once, with grouped examples rather than six repetitive operation pages. | approved |
| 14 | Link briefly to `pool.stats()` for inspection; do not duplicate the statistics guide inside `filter.rst`. | approved |
| 15 | Do not add a post-hoc maximum-Hamming selector, duplicate remover, or other library-global operation under `filter`. | approved |
| 16 | Do not implement full symmetric DUST, redesign `FilterOp` to return arbitrary metrics, add a filter registry, or perform the unrelated broad API-reference cleanup in this task. | approved |
| 17 | Keep the complexity fix minimal: count accepted symbols literally, use at least four symbols in the possible-vocabulary denominator, retain the arithmetic mean across requested k values, and add one validation guard for an empty or nonpositive `k_range`. | approved |
| 18 | Do not sort or deduplicate `k_range`; repeated values may intentionally weight a k. Return `1.0` when no requested k fits the sequence, preserving existing short-sequence behavior. | approved |

---

## 3. Required code changes

### 3.1 `filter_length`

Add a convenience method with this conceptual interface:

```python
pool.filter_length(
    min_length=None,
    max_length=None,
    name=None,
    prefix=None,
    cards=None,
)
```

Behavior:

- At least one of `min_length` or `max_length` is required.
- Bounds are inclusive.
- Bounds are nonnegative integers.
- `min_length` cannot exceed `max_length`.
- Exact-length filtering uses equal bounds.
- The method delegates to generic `filter`; it does not introduce a new
  Operation class.
- It is available on generic Pools, including DNA and protein pools.

Examples of boundary behavior:

```python
pool.filter_length(min_length=8)                 # keep length >= 8
pool.filter_length(max_length=12)                # keep length <= 12
pool.filter_length(min_length=8, max_length=12)  # keep 8 <= length <= 12
pool.filter_length(min_length=8, max_length=8)   # keep exactly length 8
```

### 3.2 Cards on all named wrappers

Add `cards: CardsType = None` to:

- `filter_length`
- `filter_gc`
- `filter_homopolymer`
- `filter_complexity`
- `filter_dust`
- `filter_restriction_sites`

Every wrapper passes `cards` unchanged to `self.filter(...)`. The available key
is `passed`, inherited from `FilterOp`.

Do not silently change default operation names. Users who want readable columns
for chained filters can supply explicit names or dict-style card aliases:

```python
pool.filter_gc(
    min_gc=0.4,
    max_gc=0.6,
    name="gc_check",
    cards={"passed": "passed_gc"},
)
```

### 3.3 Short-k complexity

The implementation must satisfy these public guarantees:

- The returned score is between 0 and 1.
- Higher values mean a richer short-k vocabulary and less repetition.
- `k_range` is nonempty and contains only positive k values (its public type is
  `tuple[int, ...]`).
- The exact aggregation rule is stated in the docstring and guide.
- Accepted IUPAC and gap characters are counted literally. The effective
  alphabet size is `max(4, number of distinct symbols in the sequence)`, so
  these symbols cannot accidentally increase the score above 1.
- Empty and very short sequence behavior is explicit and tested.

For each requested k that fits the sequence, divide the observed number of
distinct k-mers by `min(number of k-mer positions, alphabet_size ** k)`, then
take the arithmetic mean of those ratios. Ignore k values longer than the
sequence. If none fit, return `1.0`, preserving current short-sequence behavior.
Do not expand ambiguity codes, strip gaps, sort k values, or deduplicate repeated
k values. This keeps the implementation literal and avoids hidden biological
assumptions or unnecessary normalization machinery.

The documentation must call this PoolParty's **short-k linguistic-complexity
score**, not imply that its default `(1, 2, 3)` calculation or a threshold such
as `0.5` is universal.

### 3.4 DUST description

Keep the existing whole-sequence calculation and threshold direction:

- low score: less repeated triplet content;
- high score: more repeated triplet content;
- `filter_dust` keeps scores at or below `max_score`.

Correct the docstring and guide to say that PoolParty computes a DUST-style
whole-sequence triplet-repetition score. It does not run NCBI DustMasker's
64-base windowing, interval selection, or masking procedure. PoolParty's
`max_score` values must not be equated with DustMasker's command-line `level`.

The `feature/pool-stats` branch already replaces the incorrect SIMPLE citation
with Morgulis et al. (2006). Reconcile that overlapping edit rather than
duplicating or reverting it.

---

## 4. Documentation design

The canonical guide remains `docs/operations/filter.rst`.

### 4.1 Lean page structure

```text
filter
  1. Filtering model
     - generic predicate or convenience method
     - NullSeq and discard_null_seqs=True
  2. Generic filtering
     - common parameters
     - one custom predicate example not covered by a named method
  3. Ready-made checks
     - one method-selection table
     - compact grouped examples
     - one short DUST-versus-complexity explanation
  4. Chaining, cards, and output
     - one chained example
     - passed cards
     - links to stats, score, get_barcodes, generate_library, and API
```

Target a concise guide, approximately 220--280 lines. Do not repeat full common
parameter tables or rendered output for every wrapper.

### 4.2 Selection table

The page will include one table of the following form:

| Need | Method | Pass condition |
|---|---|---|
| Control sequence length | `filter_length` | Length is within inclusive bounds |
| Control GC content | `filter_gc` | GC fraction is within inclusive bounds |
| Avoid long single-character runs | `filter_homopolymer` | No run exceeds the maximum |
| Avoid concentrated triplet repetition | `filter_dust` | Score is at or below the maximum |
| Require a richer short-k vocabulary | `filter_complexity` | Score is at or above the minimum |
| Avoid recognition sites | `filter_restriction_sites` | No selected site is found |

Document `name`, `prefix`, and `cards` once as shared wrapper arguments.

### 4.3 Examples

Use three compact example groups:

1. Length, GC, and homopolymer filters chained together.
2. DUST and short-k complexity on the same small pool, showing that they can
   disagree.
3. Restriction-site filtering plus one pass/fail card example.

All displayed output must come from executed examples. Name pools explicitly so
headers match the output; do not reproduce the repository's known idealized
header mismatch.

### 4.4 Navigation and cross-references

- Keep one `filter` row and toctree entry in `state_operations.rst`.
- Expand that row to mention custom predicates and ready-made checks.
- Add stable anchors for direct links to each wrapper subsection.
- Add only the API entries required for `filter`, `filter_length`, and the five
  existing named methods.
- Update the design-card reference to state that generic and ready-made filters
  expose `passed`.
- Add a short `seealso` link to the incoming `pool.stats()` guide rather than
  restating its report.
- Update the README operations summary and the Unreleased changelog concisely.

---

## 5. Relationship to `pool.stats()`

The implemented `feature/pool-stats` branch reports:

- generated, filtered-out, valid, unique, and duplicate sequence counts;
- length range;
- GC minimum, mean, and maximum;
- longest homopolymer and fraction over a selected limit;
- DUST mean and maximum;
- selected restriction-site frequency;
- pairwise Hamming minimum, mean, and maximum for equal-length sequences.

It deliberately does not report linguistic complexity. Keep that decision: DUST
is the single general repetitiveness summary in `stats`.

The documentation relationship is:

```text
pool.stats()                             inspect a generated library
filter_*()                              enforce per-sequence rules
generate_library(discard_null_seqs=True) emit only passing sequences
```

`filter_length` is the only new filter justified directly by a statistic that
currently lacks a named per-sequence check. Do not add ordinary filters for
duplicates or Hamming distance because those depend on relationships among
multiple sequences.

The known overlaps with the now-merged `feature/pool-stats` work were:

- `utils/seq_properties.py` (DUST citation);
- `docs/api.rst`;
- `CHANGELOG.md`;
- the final `filter.rst` link from the statistics section.

They were reconciled while rebasing onto merge commit `aba0751`: the full
Morgulis et al. citation and all changelog entries were retained, the
`longest_homopolymer` API entry was added, and the filter guide now links to the
existing statistics guide. No statistics implementation was copied.

---

## 6. Tests and verification

### Code tests

- Existing filter behavior remains unchanged when `cards` is omitted.
- `filter_length` covers minimum-only, maximum-only, bounded, exact-length, and
  inclusive boundary cases.
- `filter_length` rejects missing, negative, and reversed bounds using existing
  package validation conventions.
- Every named wrapper accepts `cards=["passed"]` and dict-style aliases.
- Passing rows record `True`; rejected rows record `False` when null rows are
  retained.
- Invalid card keys continue to be rejected by the underlying card machinery.
- Complexity stays within 0--1 for canonical DNA, IUPAC, gaps, short strings,
  and representative repetitive strings.
- `k_range` validation covers empty, zero, and negative values.
- Existing DUST numeric behavior remains unchanged.

Run the existing focused filter and sequence-property tests, then the full suite
because public signatures and a shared helper change.

### Documentation verification

- Execute every displayed example using system Python.
- Build the fully configured Sphinx HTML with the existing ``poolparty_dev``
  Conda environment; do not create or use a venv.
- Confirm the `filter` page is reachable from State Operations.
- Confirm direct links to all six convenience methods resolve.
- Confirm no new dangling references or warnings are introduced.
- Inspect the rendered tables and raw HTML output manually.

Base Python includes Sphinx 7.3.7 but not all configured documentation extras.
The existing ``poolparty_dev`` Conda environment contains Sphinx 8.1.3,
``sphinx-rtd-theme``, ``sphinx-autodoc-typehints``, and ``sphinx-copybutton``.
Use that existing environment for the configured build; no installation or new
virtual environment is needed.

---

## 7. Anticipated package files

Code and tests:

- `poolparty/src/poolparty/pool_mixins/common_ops_mixin.py`
- `poolparty/src/poolparty/pool_mixins/filter_mixin.py`
- `poolparty/src/poolparty/utils/seq_properties.py`
- focused filter and sequence-property test files

Documentation:

- `poolparty/docs/operations/filter.rst`
- `poolparty/docs/operations/state_operations.rst`
- `poolparty/docs/metadata/design_cards.rst`
- `poolparty/docs/api.rst`
- `poolparty/README.md`
- `poolparty/CHANGELOG.md`

No `constraint_filters.rst` file will be created.

---

## 8. Manuscript and response follow-up

After code and documentation are verified:

1. Add the approved wording proposal to `revision/comparison/MAIN_TEX_CHANGES.md`.
2. Resolve that file's open decision #7 before editing the authoritative
   submitted manuscript.
3. Integrate the constraint-checking distinction into the existing Reviewer 3
   limitations rewrite rather than create a competing edit.
4. Finalize the point-by-point response, distinguishing constraint checking from
   sequence redesign and minimum-distance barcode construction from post-hoc
   maximum-diversity selection.
5. Mark `revision/comparison/FINDINGS.md` A1 resolved.
6. Change the PoolParty synthesis-constraint matrix cell from partial to
   supported and update its totals from `14/3/3` to `15/2/3` only after the
   documentation builds successfully.

Only `26.05.18_bmc_submission/latex/main.tex` is authoritative. Do not edit older
manuscript copies.

---

## 9. Execution order

| Step | Work | Status |
|---|---|---|
| 1 | Record decisions and documentation architecture in this plan | complete |
| 2 | Settle and record the IUPAC/gap policy for short-k complexity | complete |
| 3 | Implement complexity validation/correction and DUST wording | complete |
| 4 | Implement generic `filter_length` | complete |
| 5 | Add `cards` passthrough and tests for all six wrappers | complete |
| 6 | Rewrite `filter.rst` using the lean structure | complete |
| 7 | Make targeted navigation, API, metadata, README, and changelog edits | complete; README reviewed and unchanged |
| 8 | Run focused and full tests with system Python | complete after docs |
| 9 | Build and inspect documentation | complete using existing `poolparty_dev` Conda environment |
| 10 | Rebase onto the merged pool-stats main and reconcile the documented overlaps | complete |
| 11 | Update manuscript-change record, response, finding, and matrix | blocked on manuscript ownership |
| 12 | Perform a final cross-repository consistency review | pending |

Code-phase verification on 2026-08-22 used system Python 3.12.2 directly from
this worktree's `src` directory (no virtual environment):

- focused filter and sequence-property suite: 88 passed;
- full package suite: 2,997 passed and 14 expected failures;
- `git diff --check`: clean.

Documentation-phase verification on 2026-08-22 also used system Python with no
virtual environment:

- every displayed `filter.rst` example executed with asserted output;
- focused suite: 88 passed;
- full package suite: 2,997 passed and 14 expected failures;
- fallback structural Sphinx build: succeeded with no warnings from the changed
  filter, API, score, design-card, or navigation pages;
- all six helper links resolved to rendered API anchors, and the rendered tables
  and raw HTML content were inspected;
- the fully configured Read the Docs build succeeded under the existing
  `poolparty_dev` Conda environment with no installation; its 48 warnings are
  pre-existing project or offline-intersphinx issues, and none originate from
  the changed filter, API, score, design-card, or navigation pages.

Post-pool-stats integration verification on 2026-08-22:

- rebased the two existing commits onto `origin/main` at `aba0751` and added one
  focused documentation-integration commit;
- every displayed `filter.rst` example executed with asserted output under
  system Python 3.12.2;
- focused filter, sequence-property, statistics, export, and scan suite: 273
  passed;
- full package suite: 3,233 passed and 14 expected failures;
- fully configured Sphinx HTML build: exit 0 under the existing
  `poolparty_dev` Conda environment, with 49 pre-existing/offline warnings and
  no warning from the filter or sequence-property changes;
- all six ready-made filter links, the `pool-stats` cross-reference, the State
  Operations entry, and the `longest_homopolymer` API target resolve in the
  rendered HTML.

Post-merge verification on 2026-08-22:

- PR #21 merged using a true merge commit, `d4618fa`;
- the merge commit's first parent is the prior `main` tip `aba0751`, and its
  second parent is the exact constraint-filter tip `eee404f`;
- all three constraint-filter commits are present in `main` without rewriting;
- both local and remote `docs/constraint-filters` branch references remain at
  `eee404f`, and the local worktree is clean and synchronized with its upstream.
- after PR #19 merged, local and remote `main` advanced to `b3192f1`; the
  primary package worktree was switched back to that `main`, and the completed
  pool-stats and constraint-filter linked worktrees were removed after confirming
  that they contained no tracked or untracked changes (only generated caches and
  documentation build output).

---

## 10. Completion criteria

This work is complete only when:

- all six convenience methods are documented and API-linked;
- `filter_length` and cards behave as specified;
- complexity honors its stated range and defined input policy;
- DUST wording matches the implementation without changing its results;
- `filter.rst` is the single canonical filtering guide;
- focused tests, the full package suite, and the documentation build pass;
- examples and rendered output are verified;
- manuscript and response wording match the shipped behavior;
- the comparison finding and matrix are updated after, not before, the docs land.
