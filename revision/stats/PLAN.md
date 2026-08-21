# Reviewer 2, item 1a — library statistics

Implementation plan for `pp.stats()` / `pool.stats()`.

**Status: approved. Ready to implement.** Nothing in the package, its docs, or the
manuscript has been modified. All numbers below were measured in the read-only
venv at `poolparty-statecounter/.venv` against `poolparty` 0.1.0.

---

## 1. The comment

> "Additional statistical readouts describing pool generation would be useful.
> For example, how many unique sequences vs duplicates were produced in each pool
> (i.e. out of the universe of all permutations that meet the specifications)?
> How unique are the sequences (min/max/average hamming distance)? How frequent
> are homopolymers? Etc."

PoolParty has no such readout. This plan adds one: `pp.stats(pool)` /
`pool.stats()`, which generates a library and reports on it. Report only — it
never changes a design or a library.

---

## 2. Decisions

| # | Decision | Settled |
|---|---|---|
| 1 | **Name:** `pp.stats(pool)` and `pool.stats()`, per the package rule that every pool-taking function is also a method | yes |
| 2 | **Reported:** the composition funnel, length, GC, homopolymer runs, DUST, pairwise Hamming, restriction sites (opt-in) | yes |
| 3 | **Not reported:** melting temperature, DNA folding, reference-genome matching, edit distance, linguistic complexity, positional base composition | yes |
| 4 | **Returns** a `dict` subclass that prints as a formatted report | yes |
| 5 | **Sampling** only for pairwise Hamming, default 2,000 sequences. Everything else exact | yes |
| 6 | **Report only.** No filtering, no deduplication, and no `filter=`/`fix=`/`drop_duplicates=` parameter, ever | yes |
| 7 | **How many sequences to examine:** closed designs auto up to 1,000,000, open-ended designs require a count. See §5 | yes |
| 8 | **Wording:** `num_states` where meaningful, "unbounded" where not. Never "declared sequences" | yes |
| 9 | **Protein pools:** not supported. `ProteinPool` raises `TypeError` | yes |
| 10 | **Per-filter breakdown:** dropped. `num_filtered_out_seqs` (the total) stays | yes |
| 11 | **Restriction argument names:** keep `enzymes` / `sites`, matching `filter_restriction_sites`, `has_restriction_site` and `get_sites_for_enzymes` | yes |
| 12 | **Homopolymer argument name:** `max_homopolymer_run` | yes |
| 13 | **Printed report** hides sections with no value; `dict(s)` always shows everything | yes |
| 14 | **`to_df(unique=True)`:** not in scope. Revisit after the paper | yes |
| 15 | **`to_df` repeat-rows defect** (§3f): **fix it.** Scoped to the `num_cycles` path — see §13 | yes |
| 16 | **The auto limit is not a parameter.** Hard-coded at 1,000,000; `num_cycles=1` is the override | yes |
| 17 | **Every count key that counts sequences ends in `_seqs`.** See §4 | yes |
| 18 | **Three keys cut** as pure arithmetic nobody quotes: `length_variable`, `hamming_num_pairs`, `hamming_sd` | yes |

### Rationale for the decisions that were debated

**11 — why not `restriction_enzymes` / `restriction_sites`.** `sites` really is
vague inside `stats`, where no function name supplies the context. But the name is
used identically in three shipped functions, and the story we are telling is
*measure with `stats`, fix with `filter_*`* — users will copy arguments across.
Renaming everywhere would also make the filter itself stutter:
`filter_restriction_sites(restriction_sites=[...])`. The vagueness is fixed by the
first line of a docstring, read once; the stutter is in the code, read every time.

**10 — why the per-filter breakdown was dropped.** It duplicates something the
design graph already provides. Every step is a pool you can inspect, so calling
`.stats()` before and after a filter gives the same answer with no new machinery.
Measured, using the plain named filters:

| Pool | rows | filtered out |
|---|---:|---:|
| base | 1,710 | 0 |
| after `filter_gc` | 1,710 | 180 |
| after `filter_homopolymer` | 1,710 | 224 |

`filter_gc` rejected 180; `filter_homopolymer` rejected 224 − 180 = 44. Identical
to what the design-card route produced. This also removes the need to add a
`cards` passthrough to the five named filters, so `stats` touches no existing
behaviour at all.

**16 — why the auto limit is not a parameter.** A `stats()` call on a design
larger than the limit raises, and the escape hatch already exists in the API:
`stats(num_cycles=1)` means "yes, I mean all of it" and is never capped. An
override parameter would add a tenth argument for a case the existing API already
covers.

**18 — why those three keys and no others.** Several keys are derivable from
others — the funnel is redundant by two fields, and `open_ended` is
`num_states is None`. The rule applied: **keep a derived field when it is a number
people quote; cut it when it is arithmetic nobody asks for.** The funnel stays
whole because that chain *is* the answer to the reviewer's question and nobody
should have to reconstruct it by subtraction. `open_ended` and `hamming_exact`
stay because they are named flags carrying meaning, not arithmetic — relying on a
`None` sentinel is what trips people up. Cut: `hamming_num_pairs` (exactly
`C(n,2)`, and the printed line computes it on the fly), `length_variable`
(`length_min != length_max` is right there; still used internally to decide
whether to skip the distance section), `hamming_sd` (never printed, not what the
reviewer asked for, and the sampled mean is measurably near-exact so the spread is
not acting as a caveat).

---

## 3. What exists today (verified)

| Claim | Status | Evidence |
|---|---|---|
| No summary/QC function anywhere | confirmed | `grep -rn "def describe\|def summary\|def stats\|def qc\|def summarize\|def report" src/` → 0 hits |
| Nothing named `stats` in the package | confirmed | `grep -rn "\bstats\b" src/` → 0 hits |
| `_hamming_distance` is private, barcode-only | confirmed | `base_ops/get_barcodes.py:17`. A vectorised numpy Hamming also exists inline at `:176` — reuse that, not the scalar helper |
| No library-level deduplication | confirmed | `dedup` appears only at `utils/scan_utils.py:196-207`, keyed on `(position, insert_idx)` tuples |
| `num_states` counts state slots | confirmed | `pool.py:134` returns `self.state.num_values` |
| Filtering replaces, does not remove | confirmed | `base_ops/filter_seq.py:65` returns `NullSeq()` |
| **`filter` is the only origin of a null** | confirmed | The other four `NullSeq()` sites — `operation.py:313`, `materialize.py:93`, `translate.py:195`, `reverse_translate.py:184` — all fire only when a parent is already null. They propagate; they do not originate. This is why the key is `num_filtered_out_seqs` |
| Property helpers exist | confirmed | `utils/seq_properties.py`: `calc_gc`, `calc_complexity`, `calc_dust`, `has_homopolymer`, `has_restriction_site`, `get_sites_for_enzymes` |
| `ExportMixin`/`FilterMixin` are on `DnaPool` only | confirmed | `dna_pool.py:7` |
| The five named filters do **not** accept `cards` | confirmed | `filter_gc(min_gc, max_gc, name, prefix)` etc.; each calls `self.filter(predicate, name=, prefix=)` and drops cards |
| `filter_gc`/`filter_homopolymer`/`filter_restriction_sites` are undocumented | confirmed | 0 hits across `docs/`, `README.md`, `CHANGELOG.md` |
| 112 enzymes, 7 presets | confirmed | `data/restriction_enzymes.py` |
| `from_seq`/`from_seqs`/`get_kmers` return `DnaPool`; only `translate()` yields `ProteinPool` | confirmed | live check |

### Findings not in the original briefing

**(a) `filter` does not reduce `num_states`.** `docs/operations/library_size.rst`
lists `filter` under "State — reduces". Measured: a 30-state pool through
`filter_gc(min_gc=0.5)` still reports 30.

**(b) Two design kinds, and only one has a size.** *(Load-bearing — see §5.)*
Same design, written two ways:

```python
pp.from_seq("ACGTACGTAC").mutagenize(num_mutations=2, mode="sequential")   # num_states = 405
pp.from_seq("ACGTACGTAC").mutagenize(num_mutations=2, mode="random")       # num_states = 1
```

The second draws afresh for every row, so `num_states` has nothing to count.
Drawing from it:

| Requested | Unique | Duplicates | `num_states` |
|---:|---:|---:|---:|
| 100 | 88 | 12 (12.0%) | 1 |
| 500 | 286 | 214 (42.8%) | 1 |
| 1,000 | 372 | 628 (62.8%) | 1 |
| 5,000 | 405 | 4,595 (91.9%) | 1 |

So "one pass through the design" is meaningless for such a pool — it yields one
row — and duplicate fraction is a property of the draw size, not of the design.

The two kinds are distinguishable from the existing flag
`Operation.action_uniquely_determined_by_state` (`operation.py:161-163`,
docstring *"True if same state value always produces the same output"*), `False`
exactly for random operations built without `num_states`. Verified:

| Design | `determined_by_state` | `num_states` | 3 cycles gives |
|---|---|---:|---|
| random, no `num_states` | False | 1 | 3 rows, 3 unique |
| random, `num_states=50` | True | 50 | 150 rows, **45 unique** |
| 4-state source → random, no `num_states` | False | 4 | 12 rows, **12 unique** |
| `from_motif(pwm)` | False | 1 | open-ended |
| `from_motif(pwm, num_states=2000)` | True | 2000 | closed |

`from_motif` is random-mode-only, so every PWM-based design — including the
paper's SpliceAI example — is open-ended unless `num_states` is passed.

**Terminology used below:** a design is **closed** when every operation is
determined by its state, **open-ended** when any random operation was built
without `num_states`.

**(c) Three documentation sites are wrong.**

- `docs/pool.rst:63-65` — *"Number of distinct sequences this pool produces"*
- `docs/pool.rst:82-85` — *"the total number of distinct states (and therefore distinct sequences) the pool can produce"*
- `docs/operations/library_size.rst:4-5` — *"the number of distinct sequences it can produce"*

For a closed design these overstate by the duplicate count. For an open-ended
design the stated value is **1** while the design produces **405** distinct
sequences. `library_size.rst` contradicts itself — its own "Practical tips"
section ends *"Without `num_states`, each sequence gets a fresh draw and
`pool.num_states` is unchanged."*

**(d) Four duplicate mechanisms, not three.** `sample()` defaults to
`with_replacement=True` (`state_ops/sample.py:23`).

**(e) A second `stats()` call would differ without pinning the start point.**
`generate_library` persists `pool._current_state`. Measured on an open-ended
design, two identical calls gave 157 then 163 unique sequences. `stats` must pass
`init_state=0`, as `print_library` does. Not a user-facing argument.

**(f) `to_df` returns repeat rows on filtered designs.** Traced on a 5-slot
design with a filter that rejects 2:

```
generate_library(num_cycles=1)          to_df(num_cycles=1)
  slot 1: AAAATTTT                        row 1: AAAATTTT
  slot 2: (filtered out)                  row 2: AATTAATT
  slot 3: AATTAATT                        row 3: ACGTACGT
  slot 4: (filtered out)                  row 4: AAAATTTT   ← repeat of row 1
  slot 5: ACGTACGT                        row 5: AATTAATT   ← repeat of row 2
  → 5 rows, 3 real sequences              → 5 rows, 3 unique
```

Mechanism: `to_df` converts `num_cycles=1` into a target row count
(`target_count = num_cycles * state.num_values`), then loops
`while generated < target_count` re-calling `generate_library`
(`export_mixin.py`). `generate_library` finds 3 survivors, stops on state-space
exhaustion (`generate_library.py:173-183`), **raises a warning** —
`"Reached max_iterations (5) with only 3 valid sequences (requested 5)"` — and
persists `pool._current_state` (`:186`). `to_df` ignores the warning and asks
again; `generate_library` restarts through the design and re-emits survivors.

The refill rule is correct for `num_seqs` on a randomly-sampling design ("keep
drawing until you have N"). It is wrong for `num_cycles`, which is about passes
through a design, not a row count. Fix in §13.

**Consequence for `stats` regardless:** it calls `generate_library` directly and
never the `to_df` chunk loop, so it cannot inherit this.

**(g) `calc_dust` cites the wrong paper.** The docstring cites Hancock &
Armstrong (1994), which describes SIMPLE. The implementation is the standard DUST
score and should cite Morgulis et al. (2006).

**(h) A 1,000,000 cap already exists.** `Operation.max_num_sequential_states =
1_000_000` (`operation.py:33`) is why a single sequential operation cannot declare
more than a million states. Reused as the `stats` auto limit in §5, so no new
threshold is invented. The cap is per operation, not per pool: chaining three
single-base mutagenize steps on a 200-mer yields `num_states = 216,000,000`,
built instantly.

---

## 4. API

### Three layers

Mirrors how `generate_library` is layered, with one addition: the counting code
takes **sequences**, not a pool.

| Layer | File | Signature | Why |
|---|---|---|---|
| Counting | `utils/stats_utils.py` (new) | `_stats_from_seqs(seqs, ...) -> dict` | No `Pool`, no `Party`. Most tests become one-liners. Reusable on an imported CSV, and on another tool's output for item 2b |
| Generation | `stats.py` (new) | `stats(source, ...) -> LibraryStats` | Decides how many to generate (§5), generates, delegates |
| Method | `pool_mixins/stats_mixin.py` (new) | `DnaPool.stats(...)` | Discoverable where `to_df`/`to_file` are |

Only `num_states` and the open-ended flag need the pool. Everything else needs
sequences.

### Signature

```python
@beartype
def stats(
    source: Union[DnaPool, Sequence[str]],
    num_seqs: Optional[Integral] = None,
    num_cycles: Optional[Integral] = None,
    seed: Optional[Integral] = None,
    max_hamming_seqs: Optional[Integral] = 2000,
    max_homopolymer_run: Optional[Integral] = 6,
    enzymes: Optional[list[str]] = None,
    sites: Optional[list[str]] = None,
    show_progress: bool = True,
) -> LibraryStats:
```

`pool.stats(...)` takes the same arguments minus `source`. Nine arguments.

### Table 1 — arguments

| # | Argument | Type | Default | What it controls |
|---|---|---|---|---|
| 1 | `source` | `DnaPool \| Sequence[str]` | required | What to measure. A pool is generated first; a list of strings is measured directly |
| 2 | `num_seqs` | `int \| None` | `None` | Generate exactly this many sequences, starting from state 0. Never capped |
| 3 | `num_cycles` | `int \| None` | `None` | Generate `num_cycles × num_states` sequences. Never capped |
| 4 | `seed` | `int \| None` | `None` → 0 | Which sequences random steps produce, **and** which subset the Hamming section draws |
| 5 | `max_hamming_seqs` | `int \| None` | `2000` | At most this many sequences enter the pairwise-distance comparison. `None` skips the section |
| 6 | `max_homopolymer_run` | `int \| None` | `6` | Single-base runs *longer than* this are counted in `frac_seqs_with_long_homopolymer`. `None` omits that key |
| 7 | `enzymes` | `list[str] \| None` | `None` | Enzyme names or preset names (112 enzymes, 7 presets) |
| 8 | `sites` | `list[str] \| None` | `None` | Explicit recognition sequences; IUPAC codes allowed |
| 9 | `show_progress` | `bool` | `True` | Progress bar during generation, matching `to_df` |

`num_seqs` and `num_cycles` are mutually exclusive; passing both is an error, as
in `to_df`. There is deliberately no parameter for the auto limit (decision 16).

**Naming record.** Renamed during review: `hamming_sample_size` →
`max_hamming_seqs`, because as an argument it is a *cap* while the matching key is
the *actual number used* — on a 500-sequence library `stats(max_hamming_seqs=2000)`
compares 500, so one name for two values would repeat the `num_seqs` /
`num_valid_seqs` collision. "Sample size" is also wrong when the calculation is
exact. Removed during review: `auto_seq_limit` (decision 16). Considered and
kept: `source` (over `pool`, which would be a lie when a list is passed);
`enzymes`/`sites` (decision 11); `num_seqs`/`num_cycles`/`seed`/`show_progress`
(fixed by precedent in `generate_library` and `to_df`).

**`max_homopolymer_run` semantics.** The value is the longest run tolerated;
anything longer is counted. Measured:

| Sequence | Longest run | `=3` | `=4` | `=6` |
|---|---:|---|---|---|
| `ACGTACGTAC` | 1 | no | no | no |
| `ACGTAAACGT` | 3 | no | no | no |
| `ACGTAAAACG` | 4 | **yes** | no | no |
| `ACGTAAAAAC` | 5 | **yes** | **yes** | no |
| `ACGTAAAAAAAC` | 7 | **yes** | **yes** | **yes** |

`longest_homopolymer` — the longest run anywhere in the library — is reported
regardless of this setting.

**`enzymes` / `sites` semantics.** Combined through the existing
`get_sites_for_enzymes`; both `None` omits the section. Measured expansions:

```
enzymes=["EcoRI","BamHI"]                 -> ['GAATTC', 'GGATCC']
enzymes=["golden_gate"]                   -> ['GGTCTC','CGTCTC','GAAGAC','GCTCTTC','CACCTGC']   (8 names, 5 unique sites)
enzymes=["BsaI"], sites=["GGCGCGCC"]      -> ['GGTCTC', 'GGCGCGCC']
sites=["RGATCY"]                          -> ['RGATCY']   (R = A/G, Y = C/T)
```

Presets: `golden_gate`, `common`, `mcs`, `gibson`, `frequent_cutters`,
`rare_cutters`, `blunt`.

Reverse complements are always checked, matching
`has_restriction_site(check_rc=True)`. Measured with BsaI (`GGTCTC` /
reverse complement `GAGACC`):

| Sequence | Found? |
|---|---|
| `ACGTGGTCTCACGT` (forward) | yes |
| `ACGTGAGACCACGT` (reverse strand) | yes |
| `ACGTACGTACGTAC` | no |

`check_rc` is not exposed; `True` is right in essentially all cases.

### Fixed internally, deliberately not arguments

| Setting | Value | Why |
|---|---|---|
| `init_state` | `0` | Correctness — finding 3e. A report must be reproducible |
| `discard_null_seqs` | `False` | Filtered-out sequences must be visible to be counted |
| Auto limit | `1_000_000` | Decision 16. Reuses `Operation.max_num_sequential_states` |
| Generation chunk size | 1000 | Matches `to_df`. **Verified**: chunked generation via `init_state` reproduces a single call byte-for-byte, on both closed and open-ended designs |
| Hamming block size | 2000 | Bounds the comparison matrix's memory |
| Hamming warning threshold | 20,000 | ≈22 s, the point where a user should be told |

### Return value

`class LibraryStats(dict)` — a dict subclass with a formatted `__repr__` (~40
lines). `s["num_duplicate_seqs"]` works, `dict(s)` round-trips, it is JSON-
serialisable, and `pd.DataFrame([a.stats(), b.stats()])` gives the supplementary
table. Rejected: a bespoke non-dict class (a new public type for no gain); a bare
`dict` (no formatted output); a `pd.Series` (mixed dtypes print poorly, cannot
group sections).

### Table 2 — return keys

25 keys. Example values are from the paper's GB1 library except where noted.

| # | Key | Type | What it is | Example |
|---|---|---|---|---|
| 1 | `num_states` | `int \| None` | Size of the design; `None` when open-ended. Not a sequence count, so no `_seqs` | `547230` |
| 2 | `open_ended` | `bool` | Whether the design has a fixed size | `False` |
| 3 | `num_generated_seqs` | `int` | Sequences generated, including filtered-out ones | `547230` |
| 4 | `frac_design_covered` | `float \| None` | Fraction of the design examined; `None` when open-ended | `1.0` |
| 5 | `num_filtered_out_seqs` | `int` | Sequences a filter rejected | `0` |
| 6 | `num_valid_seqs` | `int` | Sequences that survived | `547230` |
| 7 | `num_unique_seqs` | `int` | How many different sequences there are | `546362` |
| 8 | `num_duplicate_seqs` | `int` | Excess copies — a sequence appearing 3× contributes 2 | `868` |
| 9 | `frac_duplicate_seqs` | `float` | Duplicates over `num_valid_seqs` | `0.001586` |
| 10 | `max_seq_copies` | `int` | Most copies any single sequence has | `130` |
| 11 | `length_min` | `int` | Shortest sequence | `168` |
| 12 | `length_max` | `int` | Longest sequence | `168` |
| 13 | `gc_min` | `float` | Lowest GC fraction | `0.488` |
| 14 | `gc_mean` | `float` | Mean GC fraction | `0.532` |
| 15 | `gc_max` | `float` | Highest GC fraction | `0.595` |
| 16 | `longest_homopolymer` | `int` | Longest single-base run in the library | `8` |
| 17 | `frac_seqs_with_long_homopolymer` | `float \| None` | Fraction of sequences with a run longer than `max_homopolymer_run`; `None` when that argument is `None` | `0.0000603` |
| 18 | `dust_mean` | `float` | Mean DUST repetitiveness score | `1.88` |
| 19 | `dust_max` | `float` | Highest DUST score | `2.37` |
| 20 | `frac_seqs_with_restriction_site` | `float` | Fraction of sequences containing a site; present only when `enzymes`/`sites` is given | — |
| 21 | `hamming_exact` | `bool` | Whether the distance calculation was exact or sampled | `False` |
| 22 | `hamming_seqs_compared` | `int` | How many sequences were actually compared | `2000` |
| 23 | `hamming_min` | `int` | Smallest pairwise distance | `1` |
| 24 | `hamming_mean` | `float` | Mean pairwise distance | `8.90` |
| 25 | `hamming_max` | `int` | Largest pairwise distance | `44` |

**Naming record.** Every count of sequences ends in `_seqs`, so `num_valid` alone
can never be read as "valid what?". Renamed during review: `num_rows` →
`num_generated_seqs`, `num_null` → `num_filtered_out_seqs`, `num_valid` →
`num_valid_seqs`, `num_distinct` → `num_unique_seqs`, `num_duplicates` →
`num_duplicate_seqs`, `frac_duplicates` → `frac_duplicate_seqs`, `max_copies` →
`max_seq_copies`, `frac_long_homopolymer` → `frac_seqs_with_long_homopolymer`
(the old name read as "fraction *of a sequence* that is homopolymer"),
`frac_with_restriction_site` → `frac_seqs_with_restriction_site`,
`hamming_sample_size` → `hamming_seqs_compared`.

`num_filtered_out_seqs` was chosen over `num_null_seqs`: `NullSeq` is a mechanism,
and a public report should name the user's concept. It is exactly as accurate,
since `filter` is verifiably the only origin of a null (§3), and "null" in most
codebases suggests missing data — i.e. a bug — rather than a deliberate rejection.
The printed line reads "filtered out" to match. *Caveat for the future:* if an
operation is ever changed to null sequences for a non-filter reason, this name
goes stale, and at that point the two causes should be reported separately anyway.

`longest_homopolymer` deliberately breaks the `gc_*`/`length_*`/`dust_*`/
`hamming_*` prefix-grouping convention. The grouped alternative `homopolymer_max`
is confusable with the `max_homopolymer_run` argument, and for two keys
readability beats sort order.

No `dust_min`: low DUST means "not repetitive", so the minimum carries no
information. The asymmetry with `gc_min`/`gc_mean`/`gc_max` is deliberate and
belongs in the docstring.

Cut during review: `length_variable`, `hamming_num_pairs`, `hamming_sd` —
decision 18.

### One new helper

`utils/seq_properties.py` gains `longest_homopolymer(seq: str) -> int`, matching
the key it produces. Today only the boolean `has_homopolymer(seq, max_length)`
exists; the reviewer asks *how frequent*, which needs the length. ~20 lines with
the numpydoc block. Tag stripping reuses `export_mixin._strip_tags`
(`export_mixin.py:17`).

### Where the method lives

`StatsMixin` mixed into `DnaPool`, beside `FilterMixin` and `ExportMixin`. Five
statistics are DNA-only because their matching filters are DNA-only, and
`ProteinPool` has no `to_df` or `to_file` either.

---

## 5. How many sequences to examine

A design holds no sequences, only a recipe. Duplicates cannot be counted in a
recipe, so `stats` must generate first, and something must decide how many.

| Design kind | Bare `pool.stats()` | `pool.stats(num_seqs=N)` |
|---|---|---|
| **Closed**, `num_states` ≤ 1,000,000 | report on all of it | report on N |
| **Closed**, `num_states` > 1,000,000 | `ValueError`, naming the size | report on N |
| **Open-ended** | `ValueError` — there is no "all of them" | report on N |

An explicit count is always honoured and never capped. The limit applies only
when `stats` would otherwise be guessing.

Error messages:

```
ValueError: This design declares 216,000,000 sequences, above the 1,000,000
            limit for an automatic report. Either examine part of it —
            stats(num_seqs=100_000) — or ask for all of it explicitly with
            stats(num_cycles=1).

ValueError: This design samples randomly without a fixed size, so it has no
            total number of sequences. Say how many to examine, e.g.
            stats(num_seqs=10_000). (To give the design a fixed size instead,
            pass num_states=... to the random step.)
```

### Why not the alternatives

| Option | Rejected because |
|---|---|
| Always require a count (matches `to_df`) | `pool.stats()` — the most obvious thing to type — would fail even on a 30-sequence tutorial pool |
| Default to everything, warn when large (matches `generate_library`) | A printed warning does not stop a 2.6-hour job or a 150 GB allocation |
| Always sample, e.g. 10,000 | Silently breaks the headline number. A 10,000-of-547,230 sample of GB1 sees ~2% of sequences, so a duplicated pair appears only if both copies land in it — about 0.03% of the time. It would report **zero duplicates** where the true answer is 868 |
| Derive statistics from the design without generating | Not computable. Random draws, `recombine`, `from_motif` and `filter` make uniqueness a property of realised sequences |

### Cost, measured

Generation from `benchmarks/SUPPLEMENTARY_TEXT.md`; statistics measured here.

| Library | n | generate | statistics | overhead |
|---|---:|---:|---:|---:|
| MPRA regulatory grammar | 24,000 | 8.9 s | 0.6 s | 7% |
| SpliceAI surrogate | 200,000 | 12.4 s | 2.4 s | 19% |
| GB1 deep mutational scan | 547,230 | 121.2 s | 6.1 s | 5% |

Everything except pairwise Hamming looks at each sequence once and is free
relative to generation.

### Pairwise Hamming — the one sampled number

| n | pairs | time |
|---:|---:|---:|
| 1,000 | 5.0 × 10⁵ | 0.15 s |
| 10,000 | 5.0 × 10⁷ | 5.5 s |
| 24,000 (real MPRA library) | 2.9 × 10⁸ | 23.8 s |
| 100,000 | 5.0 × 10⁹ | ~9 min (extrapolated) |
| 547,230 | 1.5 × 10¹¹ | ~4.6 h (extrapolated) |

Default `max_hamming_seqs=2000` costs ~0.3 s. If the library is no larger than
that, the calculation is exact and `hamming_exact` is `True`.

**Sampling accuracy, measured** against an exact all-pairs run on the
24,000-sequence MPRA library:

| | min | max | mean |
|---|---:|---:|---:|
| **exact**, 2.88 × 10⁸ pairs | **1** | **95** | **60.194** |
| m=2000, three seeds | 1, 2, 3 | 90, 93, 91 | 60.160 / 60.132 / 60.239 |
| m=500, three seeds | 4, 3, 5 | 88, 87, 88 | 60.078 / 60.112 / 60.246 |

The mean is estimated essentially perfectly. **The minimum is biased upward and
seed-dependent** — a sample of m from n sees (m/n)² of all pairs, which is 10⁻⁴
at m=2000, n=200,000. So the docstring, docs and printed report must state that
when `hamming_exact` is `False`, `hamming_min` is an upper bound and `hamming_max`
a lower bound, while `hamming_mean` is unbiased. Note also that
`num_duplicate_seqs > 0` settles "is the true minimum zero?" exactly, by hashing,
in O(n).

Deliberately **not** doing: a confidence interval on the mean. The C(m,2) pairs
share sequences, so a naive interval would be wrong by a large factor, and a
correct one (U-statistic variance or a jackknife) is scope creep.

**Edit distance is out.** `_edit_distance` (`get_barcodes.py:22`) is O(L²) per
pair — ~200× Hamming's cost at 200 bp.

**Hamming needs equal lengths.** Most PoolParty libraries have them (deletion
scans pad with gap characters by default). When `length_min != length_max` the
section is omitted with a note rather than computed on a subset.

---

## 6. Printed output

House style follows `print_library` and `print_dag` — plain ASCII, no
dependencies. Sections with no value are omitted. **Every number below is
measured.**

### Closed design — the paper's GB1 library

```
pool.stats()  —  547,230 of 547,230 sequences in the design

Composition
  design size (num_states)         547,230
  generated                        547,230
  filtered out                           0
  unique sequences                 546,362
  duplicate sequences                  868   (0.16%)
  most-repeated sequence           130 copies

Length
  min / max                        168 / 168

GC content
  min / mean / max                 0.488 / 0.532 / 0.595

Homopolymer runs
  longest run                            8
  sequences with run > 6            0.006%

Repetitiveness (DUST)
  mean / max                       1.88 / 2.37

Pairwise distance (Hamming)
  sampled 2,000 of 547,230 sequences
  min / mean / max                 1 / 8.9 / 44
  min is an upper bound on the true minimum; mean is unbiased
```

### Open-ended design, 500 draws

```
pool.stats(num_seqs=500)  —  500 sequences drawn

Composition
  design size                      unbounded  (random sampling, no fixed size)
  generated                              500
  filtered out                             0
  unique sequences                       278
  duplicate sequences                    222   (44.4%)
  most-repeated sequence             4 copies

  Note: this design draws randomly without a fixed size, so the duplicate count
  reflects how many you drew, not a property of the design. Drawing more raises it.

Length
  min / max                          10 / 10

GC content
  min / mean / max                   0.300 / 0.507 / 0.700

Homopolymer runs
  longest run                              3
  sequences with run > 6              0.000%

Repetitiveness (DUST)
  mean / max                         0.08 / 0.50

Pairwise distance (Hamming)
  exact, all 124,750 pairs
  min / mean / max                   0 / 3.5 / 4
```

`hamming_min = 0` because duplicates exist — two identical sequences differ at
zero positions. Consistent by construction whenever `num_duplicate_seqs > 0` and
the calculation is exact. The pair count on the "exact" line is computed on the
fly from `hamming_seqs_compared`; it is not a stored key (decision 18).

### Design with two filters

```
pool.stats()  —  1,710 of 1,710 sequences in the design

Composition
  design size (num_states)             1,710
  generated                            1,710
  filtered out                           224
  unique sequences                     1,486
  duplicate sequences                      0   (0.0%)
  most-repeated sequence               1 copy

Length
  min / max                          20 / 20

GC content
  min / mean / max                   0.450 / 0.511 / 0.600

Homopolymer runs
  longest run                              2
  sequences with run > 6              0.000%

Repetitiveness (DUST)
  mean / max                         0.86 / 1.33

Pairwise distance (Hamming)
  exact, all 1,103,355 pairs
  min / mean / max                   1 / 3.7 / 4
```

`num_states` stays 1,710 after filtering — filters replace rather than remove,
which is the documentation error at `library_size.rst`. GC reads 0.450–0.600
because the GC filter enforced exactly that.

---

## 7. Errors

| Condition | Raises |
|---|---|
| Both `num_seqs` and `num_cycles` | `ValueError: Specify only one of num_seqs or num_cycles` |
| Neither, closed design over 1,000,000 | `ValueError` naming the size, `num_seqs=`, and `num_cycles=1` |
| Neither, open-ended design | `ValueError` explaining there is no fixed size |
| `num_seqs` or `num_cycles` ≤ 0 | `ValueError` |
| `max_hamming_seqs` ≤ 1 | `ValueError` |
| `max_homopolymer_run` < 1 | `ValueError` (matches `has_homopolymer`) |
| `ProteinPool` or any non-`DnaPool` pool | `TypeError: stats() supports DnaPool; got ProteinPool` |
| `num_seqs` or `num_cycles` passed with a list of sequences | `ValueError` — refuses rather than silently ignoring |
| `max_hamming_seqs` > 20,000 | `UserWarning` with an estimated runtime, not an error |
| Zero sequences (empty, or all filtered out) | No error. `num_valid_seqs = 0`; length, GC, homopolymer, DUST and Hamming sections omitted |

---

## 8. What the readout finds (already measured)

Run on all three published libraries. This is the response letter's evidence that
1a was worth doing.

**GB1 deep mutational scan — 547,230 generated, 546,362 unique, 868 duplicates
(0.159%).** By design arm:

| arm | rows | unique | duplicates within |
|---|---:|---:|---:|
| `double` | 536,085 | 536,085 | 0 |
| `random` | 10,000 | 9,957 | 43 |
| `single` | 1,045 | 1,045 | 0 |
| `wt` | 100 | 1 | 99 |

726 further sequences appear in more than one arm. **The wild-type ORF appears
130 times — 100 by design from `repeat(times=100)`, plus 30 more because the
`mutation_rate=0.1` random arm drew zero mutations.** Verified byte-identical to
the input ORF, distributed `{wt: 100, random: 30}`. This is the finding
`max_seq_copies` surfaces.

**SpliceAI surrogate — 200,000 generated, 199,902 unique, 98 duplicates
(0.049%). And `frac_seqs_with_long_homopolymer` (runs > 6) = 1.000.** Every
sequence carries a long run, because the 201-bp background contains `TTTTTTTT` at
position 16 and the scan positions (51–89, 107–167) never overwrite it. This
library's own design *filters* homopolymers out of the 9-mer source — filtering at
the source does not constrain the assembled product.

**MPRA regulatory grammar — 24,000 generated, 24,000 unique, 0 duplicates.**
Clean. `frac_seqs_with_long_homopolymer` 0.0205, longest run 10, mean GC 0.543,
exact mean pairwise Hamming 60.19 (min 1, max 95).

Three libraries, three different stories, none visible today.

---

## 9. Scope: what "Etc." covers

The boundary is set by a survey of what comparable tools report, using the
verified records in `revision/comparison/records/` (13 tools) and the pre-scored
`scored/synthesis_constraint_checking.md` row.

| Statistic | Tools exposing it | In? |
|---|---|---|
| Length | Oligopool Calculator `lenstat`; VaLiAnT min/max + discard counts; MPRAnator `lengthsSame`; MPRA Design Tools `max_construct_size`; SeqPro `length`; DNA Chisel `SequenceLengthBounds` | **yes** (6/8) |
| GC content | Oligopool Calculator (primers); MPRAnator (barcodes); DNA Chisel `EnforceGCContent`; SeqPro; tangermeme | **yes** (5/8) |
| Pairwise Hamming | Oligopool Calculator `minimum_hamming_distance`; MPRAnator barcode distance; DNA Chisel min-Hamming tag sets; MPRA Design Tools freebarcodes; SeqPro (experimental) | **yes** (5/8) — every one uses it as a *barcode design constraint*; none reports it over a library |
| Restriction / excluded motifs | MPRAnator `DUPLICATE_RESTRICTION_SITES`; MPRA Design Tools `countDigSites`; DNA Chisel `AvoidPattern`; Oligopool Calculator `excluded_motifs` | **yes**, opt-in (4/8) |
| Homopolymer runs | Oligopool Calculator; DNA Chisel `AvoidPattern('6xA')`; MPRA Design Tools (audited runs ≥4) | **yes** (3/8) |
| Repeats / low complexity | Oligopool Calculator `maximum_repeat_length`; DNA Chisel `RepeatedKmerPattern` | **yes**, DUST only (2/8) |
| Duplicate count | none. VaLiAnT collapses silently to `_unique.csv`; MPRAnator has a keep-or-drop toggle | **yes** — nobody reports it |
| Per-rejection reporting | MPRA Design Tools `failed` tibble; VaLiAnT `_meta_excluded.csv`; Oligopool Calculator `verify` | **total only** — decision 10 |
| Tm / secondary structure | Oligopool Calculator, Mutation Maker, DNA Chisel, pydna | **no** — primer-only in every case; new dependencies |
| Reference-genome matching | Oligopool Calculator `background_directory`; DNA Chisel `AvoidBlastMatches` | **no** — external database; out of scope |
| Positional base composition | SeqPro only | **no** (1/8) |
| Edit distance | — | **no** — O(L²) per pair |
| Linguistic complexity | — | **no** — DUST measures the same property under a name reviewers recognise. `filter_complexity` still exists and is still documented; the report and the filters need not be symmetric |
| Free-space / synthesis budget | Oligopool Calculator `lenstat` | **no** — needs an `oligo_length_limit` concept PoolParty lacks |
| Cross-pool overlap | — | **no** — item 2b, addressed in `comparison/PLAN.md` §D3 |
| Stats for every node in a design | — | **no** — would generate every intermediate pool |

Two points for the response letter: the reported set covers every statistic with
two or more tools behind it; and **no surveyed tool reports a duplicate count or a
library-wide distance distribution**, so both are new contributions rather than
catch-up.

Corroboration worth noting: Oligopool Calculator's `verify` checks motif
occurrence as **excess over a library-wide baseline** — a motif present in every
oligo is baseline, not a conflict. That is exactly the SpliceAI case
(`frac_seqs_with_long_homopolymer = 1.000` from the shared background), and it
supports reporting a fraction rather than a pass/fail.

### Report only — the reasoning

1. **The write path already exists.** Five constraint filters ship today
   (`filter_mixin.py:20-278`). A `filter=` parameter on `stats` would mean two
   implementations of the same five predicates.
2. **A query that mutates is misnamed.** `s = pool.stats()` silently returning a
   different library than the user asked about is a durable source of bug reports.
3. **Exactly one cleanup the filters cannot express: deduplication.** `filter`
   takes `Callable[[str], bool]` — stateless, per sequence. A predicate cannot
   know it has seen a sequence before.
4. **And deduplication cannot be a PoolParty operation at all.**
   `_compute_one(pool, sorted_ops, global_state, ...)` sets
   `pool.state.value = global_state % num_values` and computes the row from that
   state plus the master seed; `op.compute(parents, op_rng)` receives only the
   parent sequences and an RNG, with no access to anything previously emitted.
   `docs/pool.rst:81-84` states the invariant: *"each sequence is identified by a
   state — an integer that, together with a random seed, uniquely determines the
   sequence content."* A stateful dedup step would make output depend on iteration
   order, break `init_state` random access, and reset across chunked generation.

Deduplication could therefore only live on the export path, where VaLiAnT puts
it. **Decision 14: not in this revision.** Item 2c is answered by documenting the
five filters nobody knew existed, which is a solid answer on its own; adding a
feature to match a reviewer's choice of word is the wrong reason to add a feature.

**Guard against:** a `filter=`, `drop_duplicates=` or `fix=` parameter on `stats`.
It is one line and looks helpful. The two-step form reads fine:

```python
s = pool.stats()
if s["frac_seqs_with_long_homopolymer"] > 0.05:
    pool = pool.filter_homopolymer(max_length=6)
```

That is also the response-letter sentence: **PoolParty measures with `stats` and
fixes with the `filter_*` methods**, and every number in the report has a
documented operation that addresses it.

---

## 10. Tests

New `tests/test_stats.py`. Conventions from `tests/test_seq_properties.py`: one
`Test<Thing>` class per concern, a docstring on every test, the autouse
`pp.init()` fixture from `conftest.py`. ~40 tests. Most target
`_stats_from_seqs` directly and need no `Party`.

| Class | Covers |
|---|---|
| `TestStatsFromSeqs` | counting on plain lists — `["ACGT","ACGT"]` → `num_duplicate_seqs == 1`; empty list; single sequence; all-identical |
| `TestComposition` | the funnel across all four duplicate mechanisms: random-mode collision, `repeat`, `sample(with_replacement=True)`, `filter`. Also that `num_generated_seqs == num_filtered_out_seqs + num_valid_seqs` and `num_valid_seqs == num_unique_seqs + num_duplicate_seqs` |
| `TestDesignKind` | closed → `open_ended is False`, `num_states` populated, `frac_design_covered` computed. Open-ended → `open_ended is True`, `num_states is None`, `frac_design_covered is None`, bare call raises. The three verified cases in finding 3b become three tests |
| `TestAutoLimit` | closed and under 1,000,000 → bare call works; closed and over → `ValueError` naming the size; explicit `num_seqs` above the limit → honoured, not capped |
| `TestExactHamming` | hand-computable pool (`from_seqs(["AAAA","AAAC","ACGT"])`) → known min/max/mean; `hamming_exact is True`; `hamming_min == 0` whenever `num_duplicate_seqs > 0` |
| `TestSampledHamming` | `max_hamming_seqs=2` on a 10-sequence pool → `hamming_exact is False`, `hamming_seqs_compared == 2`; same seed reproduces; value > 20,000 warns |
| `TestDnaStats` | `gc_mean` agrees with `calc_gc`; `frac_seqs_with_long_homopolymer` with `has_homopolymer`; `longest_homopolymer` with the new helper; `dust_mean` with `calc_dust`; `frac_seqs_with_restriction_site` with `has_restriction_site` — cross-checked against the existing helpers rather than restated constants |
| `TestRestrictionArgs` | `enzymes=["golden_gate"]` expands to 5 sites; `sites=["RGATCY"]` matches all four expansions; reverse-complement hits are found; neither argument → key absent |
| `TestReproducibility` | two identical calls return identical results (guards finding 3e); chunked generation equals a single call |
| `TestEdgeCases` | `length_min != length_max` → Hamming section absent, no crash; all-filtered-out pool → `num_valid_seqs == 0`, no `ZeroDivisionError`; `max_hamming_seqs=None` → section absent; `max_homopolymer_run=None` → `frac_seqs_with_long_homopolymer` absent |
| `TestFormatting` | `LibraryStats` is a `dict` (`s["num_duplicate_seqs"]`, `dict(s)` round-trip, JSON-serialisable); `repr` contains the section headings; empty sections absent from `repr`; open-ended `repr` contains the caveat note; the cut keys are absent |
| `TestRegression` | the funnel comes from `generate_library`, not `to_df` — a filtered pool must report the real duplicate count, not the refill artefact of finding 3f |
| `TestNeverMutates` | `pool.num_states`, `pool.parents` and a re-generated library are byte-identical before and after `stats()` |

Add one line to `tests/test_types.py::TestAllExports` for `stats`. The suite is
2,929 tests and passes today; that must stay true.

---

## 11. Documentation

| File | Change | Owner |
|---|---|---|
| `docs/pool.rst` | New section **"Summarising a library — `stats(...)`"**, after `to_file(...)` (ends line 415) and before `print_dag(...)` (line 418). The page already has one section per Pool readout | this work |
| `docs/pool.rst:63-65, 82-85` | Correct `num_states` per finding 3c: state slots, not distinct sequences; note the open-ended case; cross-link `stats` | this work |
| `docs/operations/library_size.rst:4-5` | Same correction; reconcile the headline with the page's own correct "Practical tips" note | this work |
| `docs/operations/library_size.rst` per-category table | `filter` row: "reduces" → "unchanged (rejected sequences become `NullSeq`)" — finding 3a | this work |
| `src/poolparty/utils/seq_properties.py` | Fix the `calc_dust` citation — finding 3g | this work |
| `docs/api.rst` | New "Library Statistics" section with `.. autofunction:: poolparty.stats`, after "Library Generation" (line ~189) | **shared — see below** |
| `CHANGELOG.md` | `[Unreleased]/Added` for `stats`; `[Unreleased]/Fixed` for §13 | this work |

No new page under `docs/operations/`: `stats` returns a dict, not a Pool, so it is
not an operation. `docs/pool.rst` is where the other non-operation readouts live.

### Overlap with the item-2c agent

They are documenting the five existing `filter_*` methods; their work touches
`docs/api.rst` and `docs/operations/`.

- **`docs/api.rst` — real collision.** They add the five filter methods; this work
  adds one `autofunction`. **Proposal: they own the file; this work hands them a
  three-line block to paste.**
- **`docs/operations/filter.rst` — theirs.** It should gain a "see also `stats`"
  pointer, in their voice.
- **`docs/operations/index.rst` toctree — theirs.** This work adds nothing.
- **`docs/operations/library_size.rst` — this work.** Confirm before either starts.
- **`utils/seq_properties.py` — this work**, adding `longest_homopolymer` and
  fixing the DUST citation.
- **No source overlap otherwise.** Their change is docs-only. Decision 10 removed
  the need to touch `filter_mixin.py` at all.

Framing point to agree before either drafts: **1a and 2c are the same
response-letter paragraph** — measure with `stats`, fix with `filter_*`. Whoever
writes it should write both halves.

---

## 12. Cost

| Item | Lines |
|---|---:|
| `utils/stats_utils.py` — the counting layer | ~200 |
| `stats.py` — generation, size decision, `LibraryStats` | ~180 (≈40 is the formatted `__repr__`) |
| `utils/seq_properties.py` — `longest_homopolymer` + citation fix | ~25 |
| `pool_mixins/stats_mixin.py` + `pool_mixins/__init__.py` + `dna_pool.py` + `__init__.py` | ~40 |
| `tests/test_stats.py` | ~280 |
| `docs/pool.rst` new section | ~100 |
| Doc corrections + CHANGELOG + the `api.rst` block | ~35 |
| `to_df`/`to_file` repeat-rows fix + tests (§13) | ~40 |
| **Total** | **~900** |

ruff, line-length 100, py310. About a day of work plus review.

Manuscript side: one paragraph (Results or Methods), one supplementary table
produced by running `stats` over the three designs already defined in
`benchmarks/paper_benchmark.py`, and one response-letter paragraph. The §8 numbers
are in hand, so the table is a script run.

**If the schedule slips:** publish the §8 table plus two paragraphs, computed by a
script in `revision/stats/`. That analysis is already done. Everything above turns
a one-off analysis into a feature users get — a better answer, but not the only
one.

---

## 13. Decision 15 — fixing the `to_df` repeat-rows defect

Independent of `stats`, which never touches that code path. Ships alongside it as
a separate, self-contained change.

**The defect** is traced in finding 3f: `to_df(num_cycles=1)` on a 5-slot design
whose filter rejects 2 returns 5 rows containing only 3 unique sequences, two of
them silent repeats. `generate_library` already detects the condition and warns —
`"Reached max_iterations (5) with only 3 valid sequences (requested 5)"` — and
`to_df` discards the warning and asks again.

**The fix.** In `ExportMixin.to_df` / `to_file`, when `generate_library` returns
fewer rows than the chunk requested, the design is exhausted: stop and return what
exists, instead of looping for more. `to_df(num_cycles=1)` on the traced design
then returns 3 rows, 3 unique — matching what `generate_library(num_cycles=1)`
already reports. Roughly five lines in one loop.

**Scope it to the `num_cycles` path.** The refill rule is **correct** for
`num_seqs` on a randomly-sampling design: "give me 1,000 sequences" reasonably
means "keep drawing until you have 1,000", and an open-ended design can always
supply more. It is wrong only where `to_df` converts a *pass count*
(`target_count = num_cycles * state.num_values`) into a *row count* and then
applies the same rule. Fixing both paths would break legitimate use; fixing
neither leaves the defect. Fix the `num_cycles` path only.

**Tests.** Add to `tests/test_export.py` (not `test_stats.py` — different
subject): the 5-slot traced design returns 3 rows and 3 unique sequences under
`num_cycles=1`; an unfiltered design still returns exactly
`num_cycles × num_states` rows; a randomly-sampling design with `num_seqs=N` still
returns N rows, refill intact.

**Documentation.** `CHANGELOG.md` gains a `[Unreleased]/Fixed` entry. If
`docs/pool.rst`'s `to_df` section states or implies the returned row count, update
it to say that a design containing a filter yields fewer rows than
`num_cycles × num_states`.

---

## 14. Also worth tracking elsewhere

`comparison/PLAN.md` lists R2 1a under "Out of scope for this document" — correct
for that file (the tool-comparison cluster), but nothing else tracked it. Point
its "Remaining work" table at this file.
