# Referee response: the tool-comparison cluster

Working plan for the 26.07.29 BMC Bioinformatics referee report.
Reference manuscript: `26.05.18_bmc_submission/latex/main.tex`.
**Not** `manuscript/main.tex` — that is an older preprint draft which disagrees
with the submission on the MPRA library size (10,000/30,000 vs 8,000/24,000) and
misnames `insertion_multiscan`. Both are already correct in the submission.

**Positioning (settled).** PoolParty's claim is **generality and flexibility**.
Code conciseness and speed are supporting detail, not the argument. Design cards
and provenance are secondary in this cluster. Reviewer 1 did explicitly ask for a
conciseness demonstration (LL3.2), so report it as a minor column — do not build
the case on it.

## The five asks

| # | Source | Ask | Status |
|---|---|---|---|
| A1 | R1 #1, #4 | Runtime + peak memory across library sizes, with environment spec | **Done** (see D4) |
| A2 | R1 #3, #5 | Feature table vs existing tools | Scoring complete (20 rows x 8 tools, `MATRIX.md`); `.tex` not yet drafted |
| A3 | R1 #3 | Code conciseness on a common use case | Folded into D3 as a minor column |
| A4 | R2 2b | Overlap of pool elements; investigate unique elements | **Not feasible as asked** — see D3 |
| A5 | R3 | Scope honesty vs VaLiAnT; qualify "replace" (L256) and "in contrast" (L274) | Evidence complete; text not drafted |
| A6 | Editor | Well-defined tool-selection process; benchmarking + features + outputs | Framework settled; see "Selection" |

A2 and A5 are the same table. A3 and A4 are the same experiment.

## Selection framework (answers A6)

No PRISMA-style screen. The Background already contains a three-tier taxonomy;
formalise it and state the inclusion rule per tier. That explains *why* n is what
it is, which is what the editor objected to ("arbitrarily picking two or three").

- **Tier 1 — general-purpose sequence toolkits.** Biopython, pydna, SeqPro, tangermeme.
- **Tier 2 — assay-specific library designers.** VaLiAnT, MPRAnator, MPRA Design Tools.
- **Tier 3 — adjacent problems, complementary.** CodonGenie, DNA Chisel, ledidi,
  Mutation Maker, and **Oligopool Calculator (`Hossain2024oc`, currently in the
  .bib but cited nowhere)**.

**Decided: no new tools go into the manuscript.** The screening record (~20
candidates with exclusion reasons, from `lit_review/novelty_search_results.md`)
goes in the **point-by-point response letter**, not the paper. Referees and editor
see the rigour; the citation list does not grow.

## D1 — evidence base: COMPLETE

Two evidence bases, each backing one table:

- **`scored/`** — Table 2. One audit per row, **one agent per row across all 8
  tools** (not per tool), which removes the rater variance that per-tool scoring
  produced. 20 rows x 8 tools = 160 cells, every one traceable.
  `scored/_recheck/` holds the batch-7 per-cell rescores, already folded in.
- **`records/`** — Table 1. Per-tool capability records for all 13 tools
  (12 competitors + PoolParty), built by a 3-stage extract -> adversarial review ->
  fix pass. This is the **only** evidence for the five tools that appear in Table 1
  but were never scored per row: CodonGenie, ledidi, Biopython, pydna, SeqPro.

Tool PDFs are **not** duplicated here — all 11 are tracked at
`../../lit_review/analyzed/` (SeqPro has no paper; software only).

**Central result:** PoolParty scores ● on 14 of 20 rows and is sole-● on three —
`Library object`, `Recombination`, `Sequence styling`. It is *not* best on six,
concentrated in Physical construction and Genomic integration. Nine cells were
corrected during scoring; **seven moved against PoolParty**. See `MATRIX.md`.

The superseded intermediate artifacts (`extractions/`, `reviews/`, `v2/`,
`verified/`, the audit-trail directories, and the earlier row sets and matrices)
were deleted from the working tree after commit `35d65d8` on branch `bmc-revision`,
where they remain recoverable.

## D2 — screening table: CUT from the manuscript

Moved to the response letter. See "Selection" above.

## D3 — expressiveness, replacing output concordance

**A4 as literally asked is not feasible.** Availability, verified 2026-08-10:

| Tool | Status |
|---|---|
| VaLiAnT | Installable (pip / Docker `quay.io/wtsicgp/valiant` / Singularity, Python >=3.11) but **dormant** — last push 2024-04-22, zero GitHub releases, not on PyPI or bioconda, AGPL-3.0 |
| MPRAnator | **Not installable** — web service only, no package. Pages serve (HTTP 200) but the design endpoint returns **HTTP 500** on its own documented use case, checked 2026-08-22. Repo has 1 commit (2015-12-27), 0 releases, **no LICENSE file** |
| MPRA Design Tools | **Not installable** — Shiny app 404, CRAN 404, and Bioconductor serves *"Removed Packages"*. GitHub reachable, last push 2017-09-26. Checked 2026-08-22 |
| Mutation Maker | The URL printed in its paper (`github.com/Merck/Mutation_Maker`) is **HTTP 404**; a live author-maintained fork exists at `github.com/ra100/Mutation_Maker` |

So a pool-overlap diff against MPRAnator / MPRA Design Tools cannot be run. Say so
in one line — it is a legitimate finding about the state of tooling, not a gap in
our work.

**Replacement experiment (agreed):** take each Tier-2 tool's **own documented
example libraries** and express them as PoolParty DAGs. This demonstrates
generality directly and is well supplied — VaLiAnT ships 7 examples, MPRAnator 6,
MPRA Design Tools 5, Oligopool Calculator 13, tangermeme 20 (inventoried per tool
in `records/*.md`, under each record's "documented examples / vignettes" heading —
several are explicitly labelled *candidates for PoolParty reproduction*). Build
these in `examples/`, one directory per attempt.

Report per attempt: can express (Y/N) · invocations needed · lines of code (minor)
· any capability that blocked it. The mirror direction — VaLiAnT **cannot** express
the GB1 library, because it has no sampled and no higher-order mutagenesis — is
the strongest single result and needs no installation to state, since it follows
from source-verified facts about `MutatorType`.

## D4 — benchmarking: DONE (another session)

`poolparty-statecounter/poolparty/benchmarks/` now contains
`paper_benchmark.py`, `paper_benchmark.json`, `make_figure.py`,
`figS_scaling.{pdf,png}`, `tableS1.{pdf,png}`, `SUPPLEMENTARY_TEXT.md`.
Coverage: **dms n=1 -> 547,230** (7 points), spliceai -> 200,000, mpra -> 24,000,
runtime and peak memory. The gap flagged earlier (DMS stopping at n=1,000) is closed.

Remaining check: confirm the environment spec in `SUPPLEMENTARY_TEXT.md` matches
the machine that produced `paper_benchmark.json`, since R1 asked for it explicitly.

## Also done by another session

`poolparty-statecounter/poolparty/examples/` now contains
`spliceai_surrogate.ipynb` (38 KB), `spliceai_design_cards_published.csv.gz`,
`build_spliceai_nb.py` and the flank files — **the Example 3 / SpliceAI migration,
previously the top open item.** Its README gained *Performance* and *Example 3
migration notes* sections (v1 -> v2 API mapping, verified equivalence, four API
traps, MaxEntScan preloading).

**Coordination risk:** multiple sessions are active on this revision. Confirm
ownership before drafting the Results subsection, or two sessions will write it.

## Remaining work

| ID | Deliverable | Blocked by |
|---|---|---|
| T1 | Main-text feature table — general, not technical; Tier 2 in full, Tiers 1/3 collapsed; include a block of capabilities PoolParty lacks | nothing — all 20 rows scored; Table 1 still needs restructuring to Table 2's tool set |
| T2 | Expressiveness results (D3) | nothing |
| T3 | Background rewrite — tier framing; decide on Oligopool Calculator citation | T1 |
| T4 | Discussion rewrite — qualify L256 "replace" and L274 "in contrast"; add the genomics-scope limitation R3 asked for | T1 |
| T5 | Point-by-point response letter, incl. the screening record | all |

## Manuscript-facing findings, not yet acted on

Recorded here only; **nothing outside `revision/` is to be modified without
explicit instruction.**

1. **R3's HGVS premise is wrong.** `hgvs_input` is "no" for all 13 tools including
   VaLiAnT. Concede genome coordinates, transcript models, split codons, VCF/VEP
   output and consequence annotation — all genuinely VaLiAnT-only — but do not
   concede HGVS.
2. **L148 claim.** "Documentation, tutorials, and example notebooks are available
   on the ReadTheDocs and GitHub project websites" — `poolparty/examples/` is
   gitignored and not distributed, and there is still no SpliceAI docs tutorial
   even though `docs/_static/images/figure4a.drawio.svg` and `figure4b_g.drawio.svg`
   are committed and referenced by nothing.
3. **Figure code panels.** Fig 2B omits `.named()` and `cards=`, so it cannot
   produce panels C and D beside it. Fig 3B does not run: `insertion_multiscan`
   missing required `num_insertions`, and `replace_region(region=...)` should be
   `region_name=`. Consider generating panel code from the executed notebooks
   rather than retyping into drawio.
4. **PoolParty limitations for R3's paragraph** (verified live): synonymous
   variants cannot be exhaustively enumerated (`mutagenize_orf` refuses
   `mode='sequential'` for non-uniform mutation types); codon-level "exhaustive"
   means exhaustive over amino acids, not codons; deletion scans emit gap
   characters and do not shorten sequences.
5. **GB1 codon 2 = `CAG` (Gln) is CORRECT** — Olson et al.'s construct carries
   T2Q. Verified against the Olson dataset itself (530,737 double mutants; Q at
   position 2 appears 511,397 times). Do not "fix" it. Full note in
   `poolparty/examples/README.md`.

## Out of scope for this document

R1 #2/#7 (design-card incremental predictive value), R1 #6 (four suggested
references — check relevance honestly; the editor wrote that including them is not
required), R2 1a (pool statistics), R2 1b (LLM test-harness disclosure), R2 2c
(new operations), R2 2d (name collision with the PoolParty popgen software,
doi:10.1111/1755-0998.12784), R3's line-level language edits.
