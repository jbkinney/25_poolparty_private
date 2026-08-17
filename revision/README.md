# `revision/` — what is what

Working directory for the BMC Bioinformatics referee response (report of 26.07.29).
Reference manuscript: `../26.05.18_bmc_submission/latex/main.tex` — **not**
`../manuscript/main.tex`, which is an older preprint draft and disagrees with the
submission on the MPRA library size and on operation names.

Three statuses:
**LIVE** — current, cite it, edit it.
**SUPERSEDED** — kept for provenance only. Do not cite, do not edit, do not hand
to an agent.
**AUDIT TRAIL** — read only when tracing how a decision was reached.

---

## Top level

| Path | Status | What it is |
|---|---|---|
| `comparison_plan.md` | LIVE | Referee asks → deliverables, current status, manuscript-facing findings not yet acted on |
| `README.md` | LIVE | this file |

## `tables/`

| Path | Status | What it is |
|---|---|---|
| `table2_candidate_17row.md` | **LIVE — the measurement instrument** | Locked 17-row set, operational definition of every row, substitution protocol, candidate rows. **This is the file scoring agents are given.** Contains no expected outcomes, by protocol rule 6 — keep it that way. |
| `table1_scope.md` | LIVE | Table 1 draft v3 — 8 rows x 6 columns (Tool, Purpose, Key features, Output, Availability, Reference). Known issue: its tool list no longer matches Table 2's 8 columns. |
| `table2_capabilities.md` | LIVE | The *current shipping* Table 2 — 8 rows x 8 columns. Remains authoritative until the 17-row set is scored. |
| `TERMINOLOGY.md` | LIVE | Controlled vocabulary, chosen by measured usage across the 11 tool papers and the manuscript. Binding on all table text. |
| `TABLE_CONVENTIONS.md` | LIVE | Survey of four real comparison tables; records every deviation we make and why |
| `table1_scope.tex` | SUPERSEDED | Two design generations behind `table1_scope.md`. Regenerate, do not edit. |

## `tool_survey/`

| Path | Status | What it is |
|---|---|---|
| `papers/` (11) | LIVE | Source PDFs. SeqPro has none — no paper exists. PoolParty's paper is the BMC submission. |
| `final/` (13) | **LIVE — the evidence base** | Per-tool capability records. Extracted → adversarially reviewed → corrected → policies applied. Best evidence we have. **Leads only when given to agents: never cite as evidence.** |
| `ESCALATION_DECISIONS.md` | LIVE | The four standing policies (balance / omissions / sourcing / version drift) and the rulings on all 62 escalations |
| `MATRIX_FINAL.md` | **LIVE — the authoritative matrix** | 20 rows x 8 tools, all 160 cells traced to `scored/`. Includes the correction log and the two contested cells. |
| `MATRIX_verified.md` | SUPERSEDED | 8 x 13 rating grid. Predates the 338 corrections; being replaced by `scored/`. Its five pre-registered rulings have been ported into the row definitions. |
| `verified/` (8) | SUPERSEDED | Per-row audits behind `MATRIX_verified.md`. Same caveat. Three of these eight rows were scored under wordings since replaced. |
| `scored/` (20 + 7 rechecks) | LIVE | Per-row audits, one file per capability, plus `_recheck/` for individually rescored cells. Assembled into `MATRIX_FINAL.md`. |
| `ROWS_v3.md` | SUPERSEDED | Defines only 8 of the 17 rows. Replaced by `table2_candidate_17row.md`. |
| `ROWS_v2.md` | SUPERSEDED | Already carries its own supersession banner |
| `v2/` (12) | SUPERSEDED | Tools re-scored against the v2 row set, which no longer exists |
| `extractions/` (13) | SUPERSEDED | Stage 1 of the original survey; absorbed into `final/` |
| `reviews/` (13) | SUPERSEDED | Stage 2 adversarial reviews; absorbed into `final/` |
| `prior_analyses/` (11) | SUPERSEDED | The pre-survey lit-review analyses. Ancestors of `final/`; none of the corrections propagated back. |
| `citation_audit/` (13) | AUDIT TRAIL | Did every quote, line reference, URL and count resolve? |
| `factcheck/` (13) | AUDIT TRAIL | What was factually wrong, what was missing, was the record balanced |
| `fixes/` (13) | AUDIT TRAIL | Change logs. Pass 1 = repair; Pass 2 = policy application. Includes what was *rejected* and why. |

---

## What to hand an agent

**Scoring / verification** — `table2_candidate_17row.md` (definitions), `papers/`,
the tool's live repo and docs, and `final/` **as leads only**.
Never give `verified/`, `MATRIX_verified.md`, `v2/`, or any prior rating: the point
is independent replication, and the diff against the old matrix is done afterwards
by hand.

**Record repair** — `final/<tool>.md`, plus that tool's `citation_audit/`,
`factcheck/` and `fixes/` entries, plus `ESCALATION_DECISIONS.md`.

**Never hand any agent** a document containing values we expect. See protocol rule
6 in `table2_candidate_17row.md` for why — an earlier version leaked our
suspicions into twelve agents at once and the echo was misread as independent
corroboration.

## Known open items

1. Table 1's tool list does not match Table 2's 8 columns (Mutation Maker and DNA Chisel are individual there, grouped here), and `sequence_styling` should move into its Key features cell.
2. Two contested cells pending review — see `tables/PARTIAL_REVIEW.md`.
3. `table1_scope.tex` needs regenerating; there is no `.tex` for Table 2.
4. Manuscript patches and the response letter are not drafted.
5. Deferred: codon-aware indel gap (restore in v2 / sharpen row 8 / leave) — `tables/PARTIAL_REVIEW.md` section 3.
