# MATRIX_FINAL — 20 rows x 8 tools

**Authoritative scored matrix.** Supersedes `MATRIX_verified.md` and `verified/`.
Every cell traces to an audit in `scored/`; row definitions are in
`../tables/ROW_DEFINITIONS.md`. Assembled 2026-08-17.

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
  Synthesis-constraint checking                 ●      ◐      ◐      ◐      ●      ○      ●      ○
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

## Open for review

Two `●` cells carried into this matrix are contested; see `../tables/PARTIAL_REVIEW.md`.

1. **`library_object` / tangermeme = ◐** — REVIEWED 2026-08-17, no change. — returns `n_best` designed sequences via
   its own API, but as a bare tensor. Contested on the grounds that a beam-search
   result set is not a *library* in the sense the row intends. Would move to `○`.
2. **`synthesis_constraint_checking` / PoolParty** — REVIEWED 2026-08-17, **corrected to ◐** (see above). Superseded text: — earned on the wording
   ("several constraint types modelled as named documented checks": `calc_gc`,
   `calc_complexity`/`calc_dust`, `has_homopolymer`, `has_restriction_site`), but we
   had expected `partial` given the absent length caps, Tm, secondary structure and
   background screening.

## Scope

Scoring covers only the 8 tools carried into Table 2. The five tools in Table 1 only
— CodonGenie, ledidi, Biopython, pydna, SeqPro — are not scored here.
