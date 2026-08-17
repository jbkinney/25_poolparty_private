# Capability matrix — 20 rows x 8 tools

**Authoritative scored matrix.** Every cell traces to an audit in `scored/`; row
definitions are in `ROW_DEFINITIONS.md`. Assembled 2026-08-17; grid/totals
reconciled 2026-08-17 (see *Provenance of corrections*).

`●` supported · `◐` partially supported · `○` not supported

```
                                            PoolP   VaLi  MPRAn  MPRAd  Oligo  MutMk  DNAch   tang
LIBRARY SPECIFICATION
  Library object                                ●      ○      ○      ◐      ◐      ◐      ○      ◐
  Composable operations                         ●      ○      ◐      ○      ●      ○      ◐      ◐
  Mixed variant types in one library            ●      ●      ●      ●      ○      ◐      ○      ○
VARIANT GENERATION
  Saturation mutagenesis                        ●      ●      ○      ○      ◐      ○      ○      ●
  Randomly sampled variants                     ●      ○      ●      ◐      ○      ○      ●      ○
  Pairwise and higher-order variants            ●      ○      ●      ○      ◐      ●      ●      ◐
  Model-guided variants                         ◐      ○      ○      ○      ○      ○      ●      ●
VARIANT TYPES
  Codon / amino-acid substitutions              ●      ●      ○      ○      ○      ◐      ●      ○
  Insertions and deletions                      ●      ●      ◐      ●      ○      ○      ○      ●
  Combinatorial multi-motif placement           ●      ○      ●      ○      ○      ○      ○      ●
  Recombination                                 ●      ○      ○      ○      ○      ○      ○      ○
  Shuffling                                     ●      ○      ◐      ○      ○      ○      ○      ●
PHYSICAL CONSTRUCTION
  Synthesis-constraint checking                 ◐      ◐      ◐      ◐      ●      ○      ●      ○
  Constraint-based optimisation                 ○      ○      ◐      ◐      ●      ●      ●      ○
  Primer design                                 ○      ○      ○      ○      ●      ●      ◐      ○
METADATA AND INSPECTION
  Per-sequence design cards                     ●      ●      ◐      ●      ●      ●      ●      ○
  Automatic naming                              ●      ●      ●      ◐      ◐      ●      ○      ◐
  Sequence styling                              ●      ○      ◐      ○      ○      ○      ○      ○
GENOMIC INTEGRATION
  Genome coordinates                            ◐      ●      ●      ◐      ○      ○      ○      ◐
  Transcript / annotation aware                 ○      ●      ○      ○      ○      ○      ○      ○
```

## Totals

| | PoolP | VaLi | MPRAn | MPRAd | Oligo | MutMk | DNAch | tang |
|---|---|---|---|---|---|---|---|---|
| ● | 14 | 8 | 6 | 3 | 5 | 5 | 7 | 5 |
| ◐ | 3 | 1 | 7 | 6 | 4 | 3 | 2 | 5 |
| ○ | 3 | 11 | 7 | 11 | 11 | 12 | 11 | 10 |

## Balance

**PoolParty sole ● (3 rows):** Library object, Recombination, Sequence styling

**PoolParty not best (6 rows):**

| Row | PoolParty | Ahead of it |
|---|---|---|
| Model-guided variants | `partial` | DNAch, tang |
| Synthesis-constraint checking | `partial` | Oligo, DNAch |
| Constraint-based optimisation | `no` | Oligo, MutMk, DNAch |
| Primer design | `no` | Oligo, MutMk |
| Genome coordinates | `partial` | VaLi, MPRAn |
| Transcript / annotation aware | `no` | VaLi |

Concessions outnumber unique wins. Both are concentrated in blocks whose headings
explain them — Physical construction and Genomic integration — rather than
scattered across the table.

## Provenance of corrections

Nine cells differ from their scoring agent's first output. Each is documented in
the relevant `scored/*.md` under a CORRECTION banner:

| Row / tool | Change | Basis |
|---|---|---|
| `saturation_mutagenesis` / MPRAnator, MPRA-DT, Mutation Maker | partial -> no | Cited the struck "exhaustiveness is the user's responsibility" branch by name |
| `mixed_variant_types` / tangermeme | partial -> no | Cited the struck "separate runs the user merges" branch by name |
| `saturation_mutagenesis` / DNA Chisel | partial -> no | Batch-7 rescore under corrected wording |
| `mixed_variant_types` / DNA Chisel, Oligopool | partial -> no | Batch-7 rescore under corrected wording |
| `genome_coordinates` / PoolParty | yes -> partial | Behavioural test: four insertion variants share one name; no coordinate cards; offsets identical across strands where arithmetic inverts |
| `synthesis_constraint_checking` / PoolParty | yes -> partial | The five named checks (`filter_gc`, `filter_homopolymer`, `filter_complexity`, `filter_dust`, `filter_restriction_sites`, `filter_mixin.py:20-194`) exist in source but appear **nowhere** in `docs/`, README or CHANGELOG. Row 19's `yes` requires named **and documented** checks; `filter.rst` documents only the generic `pp.filter(pool, predicate)`. |
| `transcript_annotation_aware` / VaLiAnT | partial -> yes | Restriction cited (one transcript per gene) matches neither `partial` branch |

Seven of the nine moved against PoolParty.

**Transcription fix, 2026-08-17.** The `synthesis_constraint_checking` / PoolParty
correction above had been applied to the totals and to *Balance*, but not to the
grid, which still showed `●`. The grid summed to 15/2/3 for PoolParty against a
stated 14/3/3. The grid has been corrected to `◐`; all eight columns now sum to
their stated totals. No other cell was affected, and the error inflated PoolParty.

## Contested-cell review — CLOSED 2026-08-17

Four cells were held for a single review pass once all rows were scored, asking of
each: does this cell earn `partial`, or is `partial` absorbing something that should
be `no`? **One of four moved.**

| Cell | Verdict |
|---|---|
| `library_object` / tangermeme ◐ | **No change.** The objection is sound — `beam_substitution(n_best=10)` returns what a search converged on, not a declared design space. But adopting that test moves only this cell (the other `partial`s all have user-determined membership), it moves *in our favour*, and `◐` already conveys "returns a set, but not a library object". `○` would overstate. |
| `synthesis_constraint_checking` / PoolParty ● | **Corrected to ◐.** The five named filters are real but undocumented — see `FINDINGS.md` A1. Row 19 requires named **and documented** checks. |
| `randomly_sampled_variants` / DNA Chisel ● (medium) | **No change.** `pick_random_mutations(n_mutations, sequence)` is documented public API taking a count parameter. Borderline against the row's exclusion of "stochastic search during optimisation", but exposed independently of the solver. |
| `composable_operations` / Oligopool ● (medium) | **No change.** Docs state "Chainable: Output of one module feeds into the next"; `primer()` and `motif()` share the DataFrame contract in either order. Inferred from documentation rather than a worked example, hence medium. |

A scan of all 160 cells found 12 scored below `high` confidence. Two were already
superseded by batch-7 rescores; of the remaining ten, none changes PoolParty's
standing — they are competitor cells where the evidence is documentary rather than
executed, concentrated in the two tools that cannot be installed and run.

## Scope

Scoring covers only the 8 tools carried into Table 2. The five tools in Table 1 only
— CodonGenie, ledidi, Biopython, pydna, SeqPro — are not scored here.
