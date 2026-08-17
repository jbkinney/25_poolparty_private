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

### B1. No codon-aware indels *(a v1 → v2 regression)*

`deletion_scan` and `insertion_scan` take no `frame` parameter and contain zero
references to codons or ORFs. `grep -rn "inframe|in_frame|in-frame"` over the
package returns nothing. An in-frame deletion requires the user to pass
`deletion_length=3` with codon-aligned `positions` computed themselves.

The codon machinery exists — `mutagenize_orf` takes an explicit `frame` (±1/2/3) —
it is simply not wired to the indel operations.

**v1 had this.** `poolparty/deletion_scan_orf_pool.py` defined
`DeletionScanORFPool`: *"Scan deletions across an ORF at the codon level …
Operations work at the codon level to maintain reading frame integrity"*, with
`deletion_size` and `start`/`end`/`step_size` all in **codon units**, plus a
position-based interface with per-position weights. `InsertionScanORFPool` likewise.
That is finer control than VaLiAnT's `inframe` mutator.

Corrected v1 → v2 mapping (an earlier note overstated the loss):

| v1 class | v2 status |
|---|---|
| `KMutationORFPool` | **absorbed** into `mutagenize_orf(num_mutations=k)` — verified: k=1/2/3 → 171 / 12,996 / 576,156 states |
| `RandomMutationORFPool` | **absorbed** into `mutagenize_orf(mutation_rate=…, mode='random')` — verified |
| `DeletionScanORFPool` | **lost** |
| `InsertionScanORFPool` | **lost** |

Two of four were consolidated into a more general operation — a refactor, arguably an
improvement. Only the two indel classes were lost.

**Table status: deferred, non-blocking.** Row 8 (`Insertions and deletions`) does not
distinguish nucleotide-level from codon-aware indels, so PoolParty's `yes` is correct
under the wording used. Options: restore in v2 (row stays `yes` honestly), sharpen
row 8 (PoolParty likely → `partial`), or leave undifferentiated.

**Do not** leave it undifferentiated *and* describe the gap in the Discussion as a
structural consequence of the DAG architecture. v1 proves it is not, and the v1 code
is public history. That framing is available for B2; it is not available here.

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

## C. Behaviour worth documenting as limitations

Verified live; all are correct behaviour that will surprise users.

1. **Synonymous variants cannot be exhaustively enumerated.** `mutagenize_orf` refuses `mode='sequential'` for any non-uniform mutation type (`mutagenize_orf.py:205-207`); `synonymous`, `missense_only_random` and `nonsynonymous_random` all raise. A DMS user building a synonymous-control arm must fall back to random sampling, which does not guarantee coverage.
2. **Codon-level "exhaustive" means exhaustive over amino acids, not codons.** `missense_only_first` enumerates 19 substitutions per codon using one representative codon each, not all synonymous DNA variants.
3. **Deletion scans emit gap characters by default** and do not shorten the sequence; `deletion_marker=None` gives true shortening. Verified: `deletion_length=3` on a 12-mer returns `---ATTTTGGGG` by default and `ATTTTGGGG` with the marker disabled.
4. **`from_motif` is random-mode only** (`from_motif.py:67`), so a PWM-derived library has no exhaustive traversal.
5. **`materialize` severs the DAG** — verified: `parents` goes from 1 to 0, losing construction history and state-space arithmetic.
6. **`num_states` counts state slots, not distinct sequences.** Random draws can duplicate, `repeat` duplicates by design, `filter` leaves nulls. Observed directly: two minus-strand shuffle rows produced identical sequences from different offsets.
7. **Dinucleotide shuffling fixes the first and last characters** (Euler-path constraint); a perfect tandem repeat has exactly one valid shuffle and returns unchanged.

---

## D. Distribution and packaging

1. **`__version__` reports `0.1.0`** (`src/poolparty/__init__.py:3`) while PyPI ships `0.1.1`.
2. **All executable notebooks are gitignored** (`.gitignore:70`) and therefore not distributed. The DMS and MPRA notebooks are regenerable from the shipped `.rst` tutorials; the SpliceAI notebook — the only in-silico worked example — does not ship at all. This bears on the manuscript claim at `main.tex:148` that "example notebooks are available on the ReadTheDocs and GitHub project websites".
3. **The assessed checkout is behind canonical `main`** — `1bb0179` (2026-04-07) vs `9d6a205` (2026-04-13). Both intervening commits touch `.github/workflows/test.yml` only, so no capability evidence moves, but the record's CI description derives from the stale file.

---

## E. Manuscript-facing, listed here for completeness

Tracked in `comparison_plan.md`; repeated so this file is a single package to-do list.

1. **Figure 2B omits `.named()` and `cards=`**, so it cannot produce panels C and D beside it. Run as printed it yields `pool[0]`…`pool[5]` and `['name','seq']`.
2. **Figure 3B does not run**: `insertion_multiscan` missing required `num_insertions`; `replace_region(region=…)` should be `region_name=`.
3. **GB1 codon 2 = `CAG` (Gln) is CORRECT** — Olson et al.'s construct carries T2Q, verified against 530,737 double mutants (Q at position 2 appears 511,397 times). Do not "fix" it.
