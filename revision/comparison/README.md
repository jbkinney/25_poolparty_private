# `comparison/` — tool comparison (Tables 1 and 2)

Answers the referee cluster asking how PoolParty compares to existing tools:
**R1 #3/#5** (feature table, code conciseness), **R2 2b** (output comparison),
**R3** (scope honesty vs VaLiAnT), and the **editor's** ask for a well-defined
tool-selection process.

Sibling of `../benchmarks/`, which answers the runtime/memory asks. No shared
artifacts.

## Files

| Path | What it is |
|---|---|
| `PLAN.md` | Referee asks → deliverables, current status, remaining work, manuscript-facing findings not yet acted on |
| `MATRIX.md` | **Authoritative cell values.** 20 rows x 8 tools; totals, balance analysis, correction provenance, closed contested-cell review |
| `ROW_DEFINITIONS.md` | **The scoring instrument.** Operational test for each of the 20 rows. The only copy — and the only file given to scoring agents |
| `TERMINOLOGY.md` | Controlled vocabulary, chosen by measured usage across the 11 tool papers. Binding on all table text |
| `CONVENTIONS.md` | Survey of four real bioinformatics comparison tables; logs every convention we deviate from, and why |
| `FINDINGS.md` | Findings about the PoolParty **package** that the survey turned up. A to-do list, not a changelog — nothing has been changed in the package |
| `table1.md` | Table 1 — tool overview, 8 rows x 6 columns |
| `table2.md` | Table 2 — display layout, substitution protocol, change log |
| `scored/` | **Table 2's evidence.** One audit per row, one agent per row across all 8 tools. `_recheck/` holds batch-7 per-cell rescores, already folded in |
| `records/` | **Table 1's evidence.** Per-tool records for all 13 tools. The only evidence for the five tools in Table 1 that were never scored per row — CodonGenie, ledidi, Biopython, pydna, SeqPro |
| `examples/` | Expressiveness experiment: rebuilding other tools' published examples in PoolParty (R2 2b, R1 #5) |

Where `MATRIX.md` and `table2.md` disagree on a value, **`MATRIX.md` wins**.

## Method

**One agent per row, not per tool.** Scoring a whole row in one pass applies a
single threshold across all eight tools. The earlier per-tool pass produced rater
variance that made cells incomparable, and is why `records/` is evidence for
Table 1 but not the source of Table 2's values.

**Scorers are never shown expected outcomes.** `table2.md` protocol rule 6 forbids
it, and `ROW_DEFINITIONS.md` is sanitised accordingly. An earlier version leaked a
"rows flagged as substitution risks" list into twelve agents, which echoed the flags
back; that was then misreported as independent nominations. Do not repeat it.

**Corrections are logged with their direction.** Nine cells moved after first
scoring; **seven moved against PoolParty**. That ratio is the argument that the
matrix was not tuned, and it only holds if it keeps being recorded.

## Reading pre-cleanup paths in `scored/` and `records/`

Those files record **what each agent actually read**, at the paths that existed when
they read it. Those statements are the audit trail — rewriting them to match the
current layout would falsify the record — so they were left alone. Translate as:

| Path as written | Now |
|---|---|
| `revision/tables/ROW_DEFINITIONS.md` | `ROW_DEFINITIONS.md` |
| `revision/tool_survey/final/` | `records/` |
| `revision/tool_survey/scored/` | `scored/` |
| `revision/tool_survey/papers/` | `../../lit_review/analyzed/` |
| `extractions/`, `reviews/`, `prior_analyses/`, `verified/`, `citation_audit/`, `factcheck/`, `fixes/`, `MATRIX_verified.md`, `ROWS_v2.md`, `ROWS_v3.md` | deleted — commit `35d65d8`, branch `bmc-revision` |

PDF citations **were** rewritten, since the files still exist and a live pointer is
more useful than a historical one. All 19 resolve.

## Standing evidence policies

Adopted while ruling on 62 escalations from the record-repair pass; still binding.

- **Sourcing — primary only.** The tool's repository, documentation site, package
  index page, or paper. Third-party aggregators (bio.tools) not admitted.
  Cross-citation between sibling records not admitted: each record stands alone.
- **Version drift — the header version governs.** A record describes the version
  named in its header. Where the current release differs materially, add one
  parenthetical; do not renumber.
- **Omissions — materiality.** Add an omitted capability only if it bears on a
  Table 2 row or a Table 1 prose cell.
- **Balance — out of scope.** These records justify a table; they are not
  deliverables and no referee reads them.

## Known open items

- **`Recombination` is a one-tool row** — PoolParty ● against seven ○. **Resolved
  2026-08-17: kept, on currency.** It carries a drafting obligation — the Results
  must state that no compared tool provides it, and why it is in the table.
  Defence in `table2.md` *Open risks*.
- **Codon-aware indels** — a v1 → v2 regression, deferred by decision. Does not
  change any row as currently worded. See `FINDINGS.md` B1.
- **Table 1 tool set** does not match Table 2's eight columns; tangermeme's
  placement is undecided.
