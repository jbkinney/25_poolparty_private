# Changes required in `main.tex`

**Nothing here has been applied.** `main.tex` is unmodified. This file is the
complete list of edits the comparison-table work implies, so the manuscript can be
changed in one deliberate pass rather than piecemeal.

Target file: `../../26.05.18_bmc_submission/latex/main.tex`
(**not** `../../manuscript/main.tex`, an older preprint draft that disagrees with
the submission on the MPRA library size and on operation names.)

Line numbers are as of commit `da98f61` and will drift once editing starts. Each
item quotes enough text to relocate.

---

## A. Preamble — 4 new packages

Insert after `\usepackage{graphicx}` (line 8):

```latex
\usepackage{booktabs}                  % \toprule \midrule \botrule
\usepackage{array}                     % >{\raggedright\arraybackslash}p{} columns
\usepackage{rotating}                  % sidewaystable, Table 1 is landscape
\usepackage[nointegrals]{wasysym}      % \CIRCLE \LEFTcircle \Circle harvey balls
```

`nointegrals` prevents wasysym from clashing with `amssymb`, already loaded at
line 9. Springer's own template (`sn-article-template/sn-article.tex`) loads
booktabs and multirow, so none of this is unusual for the class.

**Verified:** all four compile against `sn-jnl.cls` in both `referee` and
production modes with zero errors, no oversized floats, and no overfull boxes
above 3pt. Harness: `preview_tables.tex`, which mirrors this preamble exactly.

## B0. Table 1 groups DNA Chisel — decided 2026-08-18

Table 1 has **nine** rows, not ten: seven library-design tools individually, then
two grouped rows. **DNA Chisel is grouped**, under *Sequence optimization tools
(CodonGenie, DNA Chisel, ledidi)*, but **keeps its own column in Table 2**.

Why the group is named that: `main.tex` line 235 already says *"PoolParty is not a
sequence optimization tool, e.g., for optimizing codon usage [CodonGenie],
satisfying synthesis constraints on individual sequences [DNA Chisel], or
designing sequences guided by machine learning models [ledidi]"* -- the same three
tools, in that order. Table and Discussion now describe them identically.

Why the earlier label failed: *Single-sequence design tools* excluded DNA Chisel,
whose own abstract calls it *"a sequence optimization framework"* and whose Purpose
cell read *"Optimize one sequence at a time"*. Any accurate category name for
CodonGenie and ledidi also describes DNA Chisel, so the name could not be fixed
while it sat outside the group.

**Why it stays in Table 2.** DNA Chisel is ahead of PoolParty on **four of the six
rows where PoolParty is not best, the only tool in all four** -- model-guided
variants, synthesis-constraint checking, constraint-based optimization, primer
design. Removing it would drop `Model-guided variants` to tangermeme alone, a third
near-uniform row, and would show the competitor that beat us most disappearing
between submission and revision. It also scores 7 supported, third highest in the
matrix.

**Consequence to expect:** a reader tracing DNA Chisel's Table 2 column finds it
inside a Table 1 group rather than on its own row. The group label names it, so it
is findable. An earlier caption draft called this out explicitly; that was cut --
the label already shows it, and naming one competitor in a caption gives it odd
prominence.

**Open:** the Discussion positions DNA Chisel and Mutation Maker as complementary
(line 235), while Table 2 scores them as competitors across 20 library-design
capabilities where they show 11 and 12 unsupported. Both framings are defensible.
Reconcile deliberately in the D3 rewrite rather than leaving it for a referee.

## B. Insert the two tables

`\input{table1.tex}` and `\input{table2.tex}`, or paste their contents.

| Table | Label | Where | Why |
|---|---|---|---|
| Table 1 — tool overview | `tab:tools` | Background, near line 88 | It is the landscape of existing tools; the paragraph at 88 is where they are introduced |
| Table 2 — capability matrix | `tab:capabilities` | Results, or Background right after Table 1 | Depends on whether the comparison is framed as background or as a result. **Undecided.** |

Table 1 is a `sidewaystable` (landscape). This was measured, not chosen for
convenience: in portrait the float exceeds the page by 99pt in referee mode and
83pt in production, driven entirely by the *Key features* column. Portrait is
possible only by dropping a column — see `table1_options.pdf` for the three
rendered options and `table1.md` for the reasoning.

## C. One new reference

`Hossain2024oc` (Oligopool Calculator) is cited by both tables and by nothing in
the current manuscript. Accepted by decision, 2026-08-17: the tool is a scored row
and column, and showing a tool without citing it would be worse. The reference
list goes from 61 keys to 62.

`Ghazi2018aa` (MPRA Design Tools) is **already cited**, at line 88 as
`\citep[see also][]{Ghazi2018aa}`. No action.

This does not reopen the earlier decision that the citation list should not grow.
That applies to the ~20 screened candidates, which go in the response letter, not
the paper. It never applied to tools inside the tables.

## D. Discussion — deliverable T4

R3 asked for two claims to be qualified. `main.tex` contains exactly two candidate
sentences, so the mapping from the referee's PDF line numbers is unambiguous.

### D1. Line 227 — "replaces" (R3's L256)

> "We have presented PoolParty, a Python package that **replaces** the ad hoc
> scripting commonly used to design oligonucleotide libraries..."

Qualify. PoolParty does not replace ad hoc scripting for the genomics-anchored
workflows VaLiAnT serves.

### D2. Line 233 — "In contrast to" (R3's L274)

> "**In contrast to** assay-specific tools such as VaLiAnT and MPRAnator,
> PoolParty provides a single framework... Because the DAG encodes the full
> construction logic, structured design-card metadata are generated automatically
> as a byproduct."

Two separate problems in one paragraph.

**The contrast itself** needs qualifying per R3. The framework claim is defensible;
the sentence currently reads as a blanket contrast.

**The design-card clause creates a false implicature.** Placing design cards
immediately after "in contrast to VaLiAnT" invites the reading that cards are part
of the contrast. They are not — VaLiAnT scores ● on that row, and so do six of
eight tools. VaLiAnT's `META_CSV_FIELDS` is 32 typed columns including PAM
protection, sgRNA IDs and background variants; for an SGE library it records more
per-oligo fields than PoolParty does.

There is no explicit false claim anywhere in the manuscript — an earlier note in
this project overstated this. The distinction that survives scrutiny is
**composability, not completeness**: PoolParty's card schema is generated from the
user's DAG, VaLiAnT's is fixed because its pipeline is fixed. Line 233 already says
this ("Because the DAG encodes the full construction logic"); the surrounding text
is what misleads.

**Status: user deferred, 2026-08-17** — reviewers did not raise it directly, so it
is not a must-fix on its own. But this paragraph is being rewritten anyway for R3,
so the clause's wording is decided as part of T4. Do not add any claim of greater
completeness; a referee who has run VaLiAnT can refute it from the README.

### D3. Add the genomics-scope limitation R3 asked for

Line 235 already lists limitations (not a sequence optimiser, not codon
optimization, not synthesis constraints, not primer design). Extend it with the
genomics gap, which Table 2 makes visible: `genome_coordinates` ◐ and
`transcript_annotation_aware` ○ against VaLiAnT's ● on both.

**Concede:** genome coordinates, transcript models, split codons, VCF/VEP output,
consequence annotation. All genuinely VaLiAnT-only.

**Also concede here, rather than as table rows** (decided 2026-08-17):

- **GC- and length-matched background sequences drawn from a reference.** PoolParty
  has no such operation -- verified by search over the package: its control
  operations are `shuffle_seq`, `shuffle_scan`, `shuffle_states`, `rc` and
  `sample`. tangermeme has `match.extract_matching_loci` (BED + FASTA -> GC-matched
  negatives, Tutorial C1). Matched negatives are standard MPRA practice, so the
  omission is worth naming.

  *Not added as a row* because it would score roughly tangermeme ● against seven ○,
  giving the table a **third** near-uniform row and doubling its exposure to the
  "this row does not discriminate" objection we are already defending on
  `Recombination`. The symmetry it would have bought is already available:
  `Transcript / annotation aware` is VaLiAnT-alone, the exact mirror of
  `Recombination`. It would also split control generation across two rows, with
  PoolParty winning one and losing the other, which reads as arbitrary.

- **Codon-aware insertions and deletions.** Row 8 does not distinguish
  nucleotide-level from codon-aware indels, so PoolParty's ● is correct as worded;
  this is a scope concession, not a mis-scored cell. **Do not frame it as a
  structural consequence of the DAG architecture** -- v1 shipped
  `DeletionScanORFPool` and `InsertionScanORFPool` and the code is public history.
  See `FINDINGS.md` B1. If the two classes are ported back to v2 before
  resubmission, delete this concession; the table needs no change either way.

**Do not concede HGVS.** R3's premise is wrong — `hgvs_input` is `no` for all 13
tools surveyed, VaLiAnT included.

The genomic-coordinate gap **may** be framed as a structural consequence of
composability: for `recombine` (two parents), `stack` (mixed ancestry), `join`
(synthetic adapters) and `shuffle`, a per-variant coordinate is genuinely
undefined. See `FINDINGS.md` B2. That framing is **not** available for the
codon-aware indel gap — v1 shipped `DeletionScanORFPool` and `InsertionScanORFPool`
and the code is public history. See `FINDINGS.md` B1.

## E. Results — the recombination clause (required)

`Recombination` is a one-tool row: PoolParty ● against seven ○. Kept by decision,
2026-08-17, on currency grounds. That decision carries a drafting obligation.

**State plainly in the Results that recombination is the one capability no
compared tool provides, and why it is in the table.** A referee who notices an
uncontested row first reads it as padding; naming it first costs a clause and
removes the objection.

Two supporting facts for the response letter, and optionally the text:

1. **The table is symmetric on one-tool rows.** `Recombination` is PoolParty-alone;
   `Transcript / annotation aware` is **VaLiAnT-alone** (`○●○○○○○○`). Same
   structure, opposite direction. These are the only two near-uniform rows in the
   matrix. A referee arguing we invented a row only we can do has to explain the
   mirror image two blocks below.
2. **Corrections ran against us.** Nine cells moved after their scoring agent's
   first output; **seven moved against PoolParty**. That ratio is the evidence the
   instrument was not tuned.

## F. Background — narrow one sentence

Line 88:

> "General-purpose sequence toolkits [Biopython, pydna, SeqPro, tangermeme] can
> manipulate individual sequences but **have no notion of a variant library as a
> structured object**."

The citation group includes `Schreiber2025nd` (tangermeme), and the audit found
tangermeme satisfies saturation mutagenesis and composable design operations, and
scores ◐ on `library_object`. The blanket claim does not hold for it.

Either narrow the sentence, or move tangermeme out of that citation group. Table 1
already gives tangermeme its own row rather than grouping it — **that grouping is
an open decision**, see `table1.md`.

## G. Supplementary material — none from this task

**Decided 2026-08-17: no Additional file.** The 20 operational row definitions go
in the **point-by-point response letter**, not the paper.

Reasons, in order of weight:

1. The **editor** asked for the well-defined process. The response letter is
   addressed to the editor, so that is where the answer belongs.
2. It matches a decision already taken: the ~20-candidate screening record goes in
   the response letter rather than the paper, on the same logic.
3. **The caption already carries the part that matters** -- the binding global rule
   that a capability counts only where the tool supplies an operation, parameter,
   or mode for it. The per-row tests are adjudication detail.
4. Nothing in our own convention survey suggests published comparison tables ship
   operational definitions. The DMS review's columns (*Scoring approach*, *Scoring
   statistic*) are not defined anywhere either.

**Consequence for the other revision task:** the "Additional file N" vs "Table S1"
question is no longer shared. Whatever `benchmarks/` needs, it decides alone.

The full definitions live in `ROW_DEFINITIONS.md` and reach the referees via the
response letter. Keep that file current -- it is now the only copy.

## H. Captions -- DONE

Both captions were rewritten 2026-08-17 and are already in `table1.tex` and
`table2.tex`. No further action.

- Table 2's promise of "Additional file 1" is gone, per G.
- Neither caption explains the bold on PoolParty. It was added, then cut: a
  reader of this paper does not need telling which tool is ours, and Table 1's
  Reference column already reads *This work*. Figure 3's caption explains bold
  because barcodes against TFBSs in a colored sequence display are genuinely
  ambiguous; this is not that.
- Spelling switched to US throughout both tables to match `main.tex`, which uses
  *modeling* (x10), *optimiz-* (x2), *behavior* (x3), *color* (x7), *center* (x6),
  and *analyzed*. Six British forms were corrected, including the row label
  `Constraint-based optimisation` -> `Constraint-based optimization`. The internal
  slug in `scored/` was left alone; it is an identifier, not display text.
- Serial commas added, matching `main.tex` (*"structured, transparent, and
  reproducible"*; *"source, transformation, composition, and state"*).

## I. Findings unrelated to the tables

Surfaced by the tool survey; listed so the manuscript pass is complete. Detail in
`FINDINGS.md`.

| # | Line | Issue |
|---|---|---|
| I1 | 148 | "example notebooks are available on the ReadTheDocs and GitHub project websites" — `poolparty/examples/` is gitignored and not distributed. The SpliceAI notebook, the only in-silico worked example, does not ship at all. |
| I2 | Fig 2B | Omits `.named()` and `cards=`, so the printed code cannot produce panels C and D beside it. Run as printed it yields `pool[0]`…`pool[5]` and columns `['name','seq']`. |
| I3 | Fig 3B | Does not run: `insertion_multiscan` is missing the required `num_insertions`, and `replace_region(region=...)` should be `region_name=`. |
| I4 | 235 | Limitations paragraph could absorb three verified behaviours: synonymous variants cannot be exhaustively enumerated; codon-level "exhaustive" means exhaustive over amino acids, not codons; deletion scans emit gap characters rather than shortening. |
| I5 | — | **GB1 codon 2 = `CAG` (Gln) is CORRECT.** Olson et al.'s construct carries T2Q; verified against the dataset (Q at position 2 in 511,397 of 530,737 double mutants). **Do not "fix" it.** |

---

## Decisions still open

Every one of these changes the text, so none should be drafted until answered.

| # | Decision | Blocks |
|---|---|---|
| 1 | Supplementary naming — "Additional file N" vs "S1" | G, H, and the benchmarks session |
| 2 | Where Table 2 goes — Background or Results | B |
| 3 | Alphabetical order, or PoolParty first | both tables |
| 4 | tangermeme — own row, or back in the general-purpose group | F, Table 1 |
| 5 | Design-card wording inside the D2 rewrite | D2 |
| 6 | Codon-aware indels — restore in v2, sharpen row 8, or leave | D3 framing |
| 7 | Who owns `main.tex` — this session or the benchmarks session | **all of it** |

**Item 7 first.** Two sessions editing one file will collide.
