# Findings about the PoolParty package itself

Things the tool survey turned up about **PoolParty**, as distinct from decisions
about the comparison table. Each is verified against the v2 repository
(`jbkinney/poolparty-statetracker`, local checkout at commit `1bb0179`, 2026-04-07)
or by running code in the read-only venv.

Grouped by what fixing them requires. Nothing here has been changed in the package —
this is a to-do list, not a changelog.

---

## A. Documentation gaps — fixable by writing docs, no code change

### A1. Five constraint filters are undocumented *(caused a table correction)*

`pool_mixins/filter_mixin.py` defines five named, public filter operations, all live
on `DnaPool`:

| Method | Line |
|---|---|
| `filter_gc(min_gc=0.0, max_gc=1.0, …)` | `filter_mixin.py:20` |
| `filter_homopolymer(…)` | `:65` |
| `filter_complexity(…)` | `:106` |
| `filter_dust(…)` | `:151` |
| `filter_restriction_sites(enzymes=None, …)` | `:194` |

**None appears anywhere in the documentation.** A grep for all five across
`docs/**/*.rst`, `README.md` and `CHANGELOG.md` returns nothing.
`docs/operations/filter.rst` documents only the generic `pp.filter(pool, predicate)`,
and its worked example is a hand-written lambda:

```python
high_gc = pp.filter(seqs, lambda s: s.count("G") + s.count("C") >= 3)
```

**Consequence for the comparison table.** `synthesis_constraint_checking` was
initially scored `yes` for PoolParty on the strength of these methods, then
**corrected to `partial`** — the row requires constraint types modelled as named
**and documented** checks. Documenting them would legitimately restore `yes`.

This is the cheapest high-value fix in this file: five methods already written and
working, invisible to users, and costing a capability rating.

### A2. The API reference is incomplete despite its own heading

`docs/api.rst` omits public interfaces including `get_barcodes`, `annotate_orf`,
`translate`, `reverse_translate`, `stylize_orf`, `set_genetic_code`,
`set_progress_mode`, `fixed_operation`, `DnaPool`, `ProteinPool`, the DnaPool
convenience filters (A1), `to_df`, `to_file`, `clear_tags`, and the
sequence-property and restriction-enzyme helpers.

### A3. No tutorial for the in-silico application

`docs/index.rst` advertises "in silico analysis of genomic models" on the landing
page. `docs/tutorials/index.rst` ships two tutorials, neither in-silico.
`figure4a.drawio.svg` and `figure4b_g.drawio.svg` are committed and referenced by
nothing.

---

## B. Capability gaps — require code

### B1. Codon-aware indels — present, undocumented *(RECLASSIFIED 2026-08-17)*

**This entry previously said codon-aware indels were absent. That was wrong, and the
correction matters: it removed a blocker from the VaLiAnT reproduction experiment.**

What remains true: `deletion_scan` and `insertion_scan` take no `frame` parameter and
contain **zero** references to codons, ORFs, frames or triplets. `grep -rniE
"inframe|in_frame|in-frame"` over `src/poolparty/` returns nothing.

What was wrong: the conclusion that the user must compute codon-aligned positions
themselves. `positions` is typed `Sequence[Integral] | slice | None`, and the shared
`validate_positions` (`utils/seq_utils.py`) resolves a slice through
`positions.indices(...)`, **honouring the step**. So a codon scan is declarative and
tool-resolved, not hand-assembled.

Verified by execution against the installed package:

| Call | Result |
|---|---|
| `deletion_scan(deletion_length=3, positions=slice(0,None,3), deletion_marker=None, mode='sequential')` on a 6-codon ORF | **6 states, each removing one complete codon**, frame preserved |
| same with a 2-nt leader and `slice(2,None,3)` | correct — leader intact, codon boundaries respected |
| `deletion_length=6`, step 3 | works — 5 states, all 12 nt (two-codon deletions) |
| `insertion_scan(ins_pool, positions=slice(0,None,3), mode='sequential')` | **7 states**, insertions at codon boundaries |

So VaLiAnT's `inframe` mutator is reproducible today with no code change.

**The real gap is documentation, and it is inconsistent.**
`docs/operations/deletion_scan.rst` types `positions` as `list[int] | None` and
describes it as an *"Explicit list of window start positions"* — a slice is never
mentioned or demonstrated. Yet `mutagenize_orf`'s `codon_positions=slice(1, 56)` **is**
documented, in `docs/tutorials/dms_gb1.rst`. The same capability is documented for one
operation and hidden for another. Same failure mode as A1.

**Optional code, much smaller than a port:** a `frame=` argument on `deletion_scan` /
`insertion_scan` that derives the offset (and, where a region is an `OrfRegion`, reads
the frame from it) would turn a two-part idiom into one named parameter. That is sugar,
not a missing capability.

**v1 comparison, for the record.** `poolparty/deletion_scan_orf_pool.py` defined
`DeletionScanORFPool`: *"Scan deletions across an ORF at the codon level … Operations
work at the codon level to maintain reading frame integrity"*, with `deletion_size` and
`start`/`end`/`step_size` in **codon units**, plus per-position weights.
`InsertionScanORFPool` likewise. v1 expressed this in codon units with a named class;
v2 expresses it in nucleotide units with a stepped slice. The capability survived; the
vocabulary did not.

Corrected v1 → v2 mapping:

| v1 class | v2 status |
|---|---|
| `KMutationORFPool` | **absorbed** into `mutagenize_orf(num_mutations=k)` — verified: k=1/2/3 → 171 / 12,996 / 576,156 states |
| `RandomMutationORFPool` | **absorbed** into `mutagenize_orf(mutation_rate=…, mode='random')` — verified |
| `DeletionScanORFPool` | **expressible** via `deletion_scan(deletion_length=3n, positions=slice(offset,None,3))` — no named codon API |
| `InsertionScanORFPool` | **expressible** via `insertion_scan(positions=slice(offset,None,3))` — no named codon API |

Nothing was lost. Two classes were generalised and two were re-expressed as idioms.

**Table status: unchanged.** Row 8 (`Insertions and deletions`) does not distinguish
nucleotide-level from codon-aware indels, so PoolParty's `yes` stands either way.

**Do not** describe this gap in the Discussion as a structural consequence of the DAG
architecture. It never was, and now it is not even a gap — only an undocumented idiom.
That framing remains available for B2 and for the non-uniform enumeration limit in C1.

### B2. Genomic coordinates are not composed for designed variants

`from_fasta` in batch mode emits a coordinate-bearing name (`chr1:10-30(+)`), but
nothing composes it into a per-variant position. Verified behaviourally:

- Four insertion variants at different offsets all carried the **identical** name; output columns after a design operation are `['name', 'seq']` only.
- `from_fasta`'s valid card keys are `['seq', 'seq_index', 'seq_name', 'state']`; requesting `chrom`/`start`/`stop`/`strand` is rejected. The origin exists only as a packed string needing regex parsing.
- Offsets are emitted but operation-namespaced (`op[2]:insertion_scan(region_scan).start`) and **identical across strands**: `chr1:10-30(+)` and `chr1:10-30(-)` both give `0,1,2`, yet genomic position is `start+offset` on plus and `stop-1-offset` on minus. **The tool gives no warning; composing them wrongly fails silently.**
- With a single coordinate tuple rather than a list, the name is `None` entirely.

**Two tiers of fix.** *Cheap:* give `from_fasta` real `chrom`/`start`/`stop`/`strand`
card keys instead of only the packed `seq_name` — removes the regex step and makes
strand explicit as data. This would **not** change the table score, which needs a
per-variant coordinate. *Expensive:* propagate a coordinate mapping through the
operation set — and for `recombine` (two parents), `stack` (mixed ancestry), `join`
(synthetic adapters) and `shuffle` the answer is genuinely undefined. Unlike B1, this
one **is** a structural consequence of composability and can be framed as such.

---

### B3. `deletion_scan` cannot be used twice with different deletion lengths

`deletion_scan` hard-codes the name of its internal marker region —
`marker_name = "_del"` (`scan_ops/deletion_scan.py:81`), passed as `tag_name` to
`region_scan`, with **no parameter to override it**. Region lengths are Party-global,
so a second `deletion_scan` with a different `deletion_length` collides:

```
ValueError: Region '_del' already registered with seq_length=3,
cannot re-register with seq_length=6. Region lengths must be consistent
within a Party.
```

Reproduced directly: `deletion_length=3` then `deletion_length=6` in one process,
on independent pools. Each works alone; together the second raises.

**Why it matters for the comparison.** VaLiAnT routinely combines deletion spans in a
single library — `1del` with `2del0` or `2del1` appears in **three of its five shipped
examples** (`brca1_nuc`, `cdna`, `brca1_prime_editing`). None of those is expressible
in one PoolParty Party today. It also blocks any mixed 1-codon / 2-codon in-frame
deletion series.

`brca1_pep` is unaffected — it uses a single deletion length — which is part of why it
is the right primary reproduction target.

**Fix:** expose the marker name, or derive it from `deletion_length` / the operation's
prefix, so distinct lengths get distinct regions.

### B4. Prefixed insertion scans cannot be copied or deep-copied

Verified on the Phase 2 ORF-scan branch `feature/orf-indel-scans` (2026-08-21).
Both generic `insertion_scan` and the new `insertion_scan_orf` add a
`PassthroughOp` when any of `prefix`, `prefix_position`, or `prefix_insert` is
provided. The resulting pool generates and composes correctly, but `copy()` and
`deepcopy()` fail:

```python
import poolparty as pp

with pp.Party():
    pool = pp.insertion_scan(
        "AAACCC", "TAG", mode="sequential", prefix="variant"
    )
    pool.copy()
```

```text
TypeError: PassthroughOp.__init__() got an unexpected keyword argument 'name'
```

The same failure occurs with `insertion_scan_orf` and with either of the other
two prefix arguments. Unprefixed insertion scans copy and deep-copy successfully.
Normal generation, molecular output, design cards, state counts, naming, and
downstream composition are unaffected.

**Root cause.** The generic operation-copy machinery reconstructs an operation
with `name=`, but `PassthroughOp.__init__` does not accept that argument. Adding
the argument alone is not a complete fix: the insertion naming callback is a
closure over the original position and insertion state objects, so a deep copy
could still calculate names from the original DAG rather than its copied states.

**Fix direction.** Replace the closure-based naming wrapper with one reusable,
state-rebindable naming operation used by both generic and ORF insertion scans.
Its position and insertion state dependencies must be explicit so `copy()` and
`deepcopy()` can rebind them to the copied DAG. Fix this once in generic
infrastructure; do not add an ORF-only workaround.

**Status:** unresolved and intentionally out of scope for the ORF indel-scan
feature. It is release-blocking only if copying a prefixed insertion pool is a
required supported workflow.

## C. Behaviour worth documenting as limitations

Verified live; all are correct behaviour that will surprise users.

1. **Synonymous variants cannot be exhaustively enumerated, and the reason is architectural.** `mutagenize_orf` refuses `mode='sequential'` for any non-uniform mutation type (`mutagenize_orf.py:205-207`); `synonymous`, `missense_only_random` and `nonsynonymous_random` all raise. The mechanism is explicit in `codon_table.py:8-13`: `UNIFORM_MUTATION_TYPES` is a dict of **fixed alternatives per codon** — `any_codon: 63`, `nonsynonymous_first: 20`, `missense_only_first: 19`, `nonsense: 3` — consumed as `self.uniform_num_alts`. Synonymous substitution has no fixed count (Met and Trp have none; Leu, Ser and Arg have five), so the per-operation state count is not a constant and the product arithmetic behind lazy state decomposition does not apply. A DMS user building a synonymous-control arm must fall back to random sampling, which does not guarantee coverage. **Unlike B1, this one may honestly be described as a consequence of the architecture** — it blocks VaLiAnT's `snvre` mutator, which appears in four of its five shipped examples.
2. **Codon-level "exhaustive" means exhaustive over amino acids, not codons.** `missense_only_first` enumerates 19 substitutions per codon using one representative codon each, not all synonymous DNA variants.
3. **Deletion scans emit gap characters by default** and do not shorten the sequence; `deletion_marker=None` gives true shortening. Verified: `deletion_length=3` on a 12-mer returns `---ATTTTGGGG` by default and `ATTTTGGGG` with the marker disabled.
4. **`from_motif` is random-mode only** (`from_motif.py:67`), so a PWM-derived library has no exhaustive traversal.
5. **`materialize` severs the DAG** — verified: `parents` goes from 1 to 0, losing construction history and state-space arithmetic.
6. **`num_states` counts state slots, not distinct sequences.** Random draws can duplicate, `repeat` duplicates by design, `filter` leaves nulls. Observed directly: two minus-strand shuffle rows produced identical sequences from different offsets.
7. **Dinucleotide shuffling fixes the first and last characters** (Euler-path constraint); a perfect tandem repeat has exactly one valid shuffle and returns unchanged.
8. **There is no "substitute to a specific amino acid" mutation type.** `VALID_MUTATION_TYPES` (`codon_table.py:16-24`) offers only *classes* of substitution — `any_codon`, `nonsynonymous_first/random`, `missense_only_first/random`, `synonymous`, `nonsense`. Targeting one named residue at every codon is not expressible, so VaLiAnT's `ala` mutator (alanine scanning, a standard DMS design) has no PoolParty equivalent. `nonsense` is the nearest analogue to VaLiAnT's `stop`, but yields 3 alternatives per codon rather than one representative, so counts will not match.

---

## D. Distribution and packaging

1. **`__version__` reports `0.1.0`** (`src/poolparty/__init__.py:3`) while PyPI ships `0.1.1`.
2. **All executable notebooks are gitignored** (`.gitignore:70`) and therefore not distributed. The DMS and MPRA notebooks are regenerable from the shipped `.rst` tutorials; the SpliceAI notebook — the only in-silico worked example — does not ship at all. This bears on the manuscript claim at `main.tex:148` that "example notebooks are available on the ReadTheDocs and GitHub project websites".
3. **The assessed checkout is behind canonical `main`** — `1bb0179` (2026-04-07) vs `9d6a205` (2026-04-13). Both intervening commits touch `.github/workflows/test.yml` only, so no capability evidence moves, but the record's CI description derives from the stale file.

---

## E. Manuscript-facing, listed here for completeness

Tracked in `PLAN.md`; repeated so this file is a single package to-do list.

1. **Figure 2B omits `.named()` and `cards=`**, so it cannot produce panels C and D beside it. Run as printed it yields `pool[0]`…`pool[5]` and `['name','seq']`.
2. **Figure 3B does not run**: `insertion_multiscan` missing required `num_insertions`; `replace_region(region=…)` should be `region_name=`.
3. **GB1 codon 2 = `CAG` (Gln) is CORRECT** — Olson et al.'s construct carries T2Q, verified against 530,737 double mutants (Q at position 2 appears 511,397 times). Do not "fix" it.
