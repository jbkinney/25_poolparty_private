# `revision/` — BMC Bioinformatics referee response

Working directory for the referee report of 26.07.29.

**Reference manuscript:** `../26.05.18_bmc_submission/latex/main.tex` — **not**
`../manuscript/main.tex`, which is an older preprint draft and disagrees with the
submission on the MPRA library size and on operation names.

**Nothing outside `revision/` is modified without explicit instruction.**

## Two independent tasks

| Directory | Referee asks | Owner |
|---|---|---|
| `benchmarks/` | R1 #1, #4 — runtime and peak memory across library sizes, with environment spec | separate session |
| `comparison/` | R1 #3, #5 · R2 2b · R3 · editor — feature comparison, output comparison, scope honesty, tool-selection process | this work |

They are siblings and do not share artifacts. See each directory's `README.md`.

## Conventions

- Tool PDFs are **not** duplicated here. All 11 are tracked at
  `../lit_review/analyzed/`.
- Superseded material is **deleted, not archived**. It is recoverable from
  commit `35d65d8` on branch `bmc-revision`, which snapshots the working directory
  as it stood before cleanup. Version control does this job; a parallel `_archive/`
  tree would only duplicate it.
