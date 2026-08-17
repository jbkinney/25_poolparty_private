# Table 2 — LOCKED row set, 20 rows x 8 tools

> **STATUS: row set locked 2026-08-14; all 20 rows scored 2026-08-17.** Rows may
> still be *substituted* — see the protocol below — but the set is otherwise
> closed, and the count stays at 20. Substitution is now expensive: every row has
> a scored column that a swap discards.

This is the shipping row set. The earlier 8-row draft (`table2_capabilities.md`)
was superseded and deleted; it is recoverable from commit `35d65d8` on branch
`bmc-revision`. Cell values are authoritative in `MATRIX.md`; row definitions in
`ROW_DEFINITIONS.md`. This file owns the **display layout**, the **substitution
protocol**, and the **change log**.

## Substitution protocol

A better candidate row may displace an existing one. To keep that from becoming
uncontrolled churn, and from wasting scored cells:

1. **One in, one out.** The set stays at 20. A candidate must name the row it
   replaces.
   *Amended 2026-08-15: the set was 17; `shuffling`,
   `synthesis_constraint_checking` and `primer_design` were each added without a removal by explicit
   decision, so the ceiling is now 20.*
2. **A candidate must beat the incumbent on one of three grounds, with evidence:**
   - *discrimination* — it produces more distinct values across the 8 tools;
   - *honesty* — it is a row where PoolParty is not best, replacing one where it is; or
   - *currency* — the incumbent measures a capability that no longer matters to
     library design as practised now, and the candidate does.
   Intuition alone is not sufficient; log the evidence that motivated it.

   **Currency is scoped to the paper's three declared applications** — DMS, MPRA,
   and in-silico probing of genomic models. A capability that is marginal for two
   of them but current for the third still passes. This is discriminating, not
   merely permissive: it retains `Recombination` (marginal for wet-lab MPRA, but
   standard practice for in-silico attribution work) while rejecting
   *degenerate-codon compression* (a synthesis-cost optimisation — dead in silico,
   where nothing is synthesised, and largely obsoleted elsewhere by falling pool
   prices). Compression also sits in direct tension with the paper's own thesis:
   a compressed library cannot carry per-variant design cards, because you do not
   know which molecule is which until you sequence it.
3. **Substitutions are cheap before verification and expensive after.** Scoring
   the outstanding cells costs roughly 88 cell-assessments. A row swapped after
   that pass discards its scored column and needs a fresh one.
   **Propose substitutions before the verification pass launches.**
4. **The verification pass may itself propose substitutions.** Agents scoring a
   row will be told they may return "this row does not discriminate — replace it
   with X", with evidence. That is how three of the previous eight rows were found
   defective, so it is the likeliest source of good candidates.
5. **Log every substitution here**, with the date, the row removed, the row added,
   and the ground under rule 2. A row set that changes without a record is how the
   pydna drift happened.
6. **This document is handed to scoring agents, so it must not contain expected
   outcomes.** No note may say which rows we suspect are uniform, which values we
   expect, or which way a cell is likely to fall. Learned the hard way: an earlier
   version carried a "rows already flagged as substitution risks" section; the
   twelve Pass-2 agents read it and echoed the flags back, which was then
   misreported as six independent nominations. Predictions about data we are about
   to collect belong nowhere near the people collecting it.

### Substitution log

| Date | Change | Ground |
|---|---|---|
| 2026-08-15 | Row 11 renamed `constraint_repaired_variants` -> `constraint_based_optimisation`. Label only; the operational test was unchanged, so scores stand. | Terminology: "optimisation" appears 110x in 11/11 tool papers, "repair" 7x in 3 and never in the sequence-design sense. |
| 2026-08-17 | Block restructure, no rescoring: `How variants are generated` -> **Variant generation**; `What variants can be generated` -> **Variant types**. Moved `mixed_variant_types` to Library specification (it describes what one specification can express, not how variants are chosen) and `model_guided_variants` to Variant generation (it is a way of arriving at variants, not a kind of change). | Consistency — all six block names are now noun phrases, and each block answers one question. |
| 2026-08-17 | Merged batch-7 rechecks: `saturation_mutagenesis`/DNA Chisel `partial`->`no`; `mixed_variant_types`/DNA Chisel `partial`->`no`; `mixed_variant_types`/Oligopool `partial`->`no`. Three cells unchanged. | Rescored under corrected definitions after the struck "user's responsibility" and "separate runs" clauses. |
| 2026-08-15 | **Removed from display** `Library algebra` and `On-demand generation`. They were in an early 19-row proposal, dropped when the set was trimmed to 17, and were never defined in `ROW_DEFINITIONS.md` nor scored in any batch. They were reinstated by mistake when the display table was rewritten by hand to add the Physical construction block; this entry records the correction, not a new decision. | Bookkeeping. Both are rows PoolParty held alone in the original 8-row matrix, so reinstating them would have added two sole-● rows unnoticed. |
| 2026-08-15 | **Added** row 20 `primer_design`; created a **Physical construction** block holding rows 19, 11 and 20; moved `model_guided_variants` out of it into variant generation. | Coherence and honesty. The block matches the Discussion's own scoping sentence; `model_guided_variants` concerns biological function, not manufacturability. Label measured: `primer` 270 uses in 5/11 papers vs `oligo design` 14 in 1. |
| 2026-08-15 | **Added** row 19 `synthesis_constraint_checking` to the variant block, adjacent to row 11. No removal. | Honesty and discrimination. Three tools (Oligopool, Mutation Maker, DNA Chisel) at `yes` against PoolParty `partial`; a genuine three-level spread across all eight — the best discrimination of any candidate identified. |
| 2026-08-15 | **Added** row 18 `shuffling` to the variant block. `sequence_styling` retained rather than displaced, so the set grows 17 -> 18. | Discrimination and currency. Shuffling is a variant type with a real threshold (dinucleotide-preserving), current for both MPRA controls and in-silico nulls, and PoolParty ties rather than wins. |
| 2026-08-17 | **Retained** row 10 `recombination` despite failing rule 2's discrimination ground (PoolParty ● against seven ○). No substitution; the set stays at 20. Decision by the user. | Currency, under rule 2's scoping to the paper's three declared applications: marginal for wet-lab MPRA, standard for in-silico attribution. Logged because rule 5 covers decisions *not* to change the set as well as changes — an uncontested row that nobody recorded a reason for is how a padding criticism lands unanswered. Defence and required Results text in *Open risks* below. |

### Candidate rows — documented, NOT integrated

Held here so they are not lost and not adopted by drift. Neither has displaced a
locked row; neither should be scored in the verification pass unless explicitly
promoted first.

Both already have unverified values in the 28-row v2 matrix, which is why they were
noticed — they are not new inventions, they are existing rows that may deserve
promotion into the shipping set.

---

#### Candidate A — Control / background sequence generation  **[PARTLY PROMOTED]**

> The shuffling sub-case was promoted 2026-08-15 as row 18 `shuffling`, narrowed to
> permutation of a sequence's own residues with a dinucleotide-preserving threshold.
> What remains a candidate is the rest of this row: **GC- and length-matched
> background sequences drawn from a reference** (e.g. tangermeme's `match.py`), which
> is not shuffling and is not scored anywhere in the set.

**Proposed definition.** Generation of sequences intended as negatives or baselines
rather than as designed variants: scrambles; shuffles that preserve composition
(mononucleotide or, more stringently, dinucleotide frequency); reverse and
reverse-complement controls; or GC- and length-matched background sequences drawn
from a reference. `partial` = one control type only, or shuffling that does not
preserve composition.

*Dinucleotide-preserving shuffling is the discriminating sub-case* — it is the
standard null for motif and model-attribution work, because dinucleotide frequency
alone drives much of a model's output.

**Existing v2 values (unverified), across the 8 shipping columns:**
`PoolParty ● · VaLiAnT ○ · MPRAnator ● · MPRA Design Tools ○ · Oligopool ○ ·
Mutation Maker ○ · DNA Chisel ○ · tangermeme ●` — 3 ● / 5 ○.

| Criterion | Assessment |
|---|---|
| Discrimination | adequate — 3 ● against 5 ○ |
| Honesty | **PoolParty ties, does not win** |
| Currency (scoped) | **strong** — matched negatives are mandatory in MPRA practice, and shuffled nulls are standard in in-silico model probing. Passes on two of the paper's three applications. |

---

#### Candidate B — Synthesis-constraint checking  **[PROMOTED 2026-08-15 → row 19]**

> Promoted in full. The definition below was tightened on promotion: `yes` now
> requires several constraint types modelled as named documented checks, and a
> boundary against row 11 was added (row 19 detects, row 11 fixes).

**Proposed definition.** The tool checks designed sequences against physical
synthesis or cloning constraints — oligo length caps, homopolymer runs, restriction
site occurrences, extreme GC — and reports or filters violations. `partial` = a
generic predicate mechanism through which a user can express such a check, without
the constraint types being modelled by the tool.

**Existing v2 values (unverified), across the 8 shipping columns:**
`PoolParty ◐ · VaLiAnT ◐ · MPRAnator ◐ · MPRA Design Tools ◐ · Oligopool ● ·
Mutation Maker ● · DNA Chisel ● · tangermeme ○` — 3 ● / 4 ◐ / 1 ○.

| Criterion | Assessment |
|---|---|
| Discrimination | good — a genuine three-level spread |
| Honesty | **PoolParty is ◐ against three tools at ●** |
| Currency (scoped) | **strong for wet-lab scopes** — oligo length caps and cloning-site constraints still bind hard regardless of synthesis price. Weak for in-silico, where nothing is synthesised. |

---

#### Why neither is integrated

The locked set is at 18 and the protocol is one-in-one-out. Promoting either
requires naming the row it replaces, and no incumbent has yet been shown to fail
on discrimination, honesty or currency. Score the locked set first, then decide
from the resulting distributions what — if anything — either of these displaces.


**What changed and why.** The shipping draft has PoolParty scoring ● on all eight
rows — the single most attackable thing in it. This expansion adds rows where
PoolParty is *not* best, so the table reads as a comparison rather than a feature
list. It also adds the variant-type axis, which is where PoolParty's breadth
actually lives.

**Columns: 8.** General-purpose toolkits (Biopython, pydna, SeqPro) are dropped —
they do not design libraries, Table 1 places them, and the same rule that excluded
ledidi excludes them. Verified consequence: dropping the column removes no
concession. PoolParty gains no row where it becomes uniquely ●.

20 x 8 = 160 cells, against the DMS review's 120.

---

## The table

```
                                          Pool  VaLi  MPRA  MPRA   Oligo  Muta   DNA    tanger
                                          Party AnT   nator Design pool   tion   Chisel meme
                                                            Tools  Calc   Maker
───────────────────────────────────────────────────────────────────────────────────────────────
LIBRARY SPECIFICATION
  Library object                            ●     ○     ○     ◐      ◐      ◐      ○      ◐
  Composable operations                     ●     ○     ◐     ○      ●      ○      ◐      ◐
  Mixed variant types in one library        ●     ●     ●     ●      ○      ◐      ○      ○

VARIANT GENERATION
  Saturation mutagenesis                    ●     ●     ○     ○      ◐      ○      ○      ●
  Randomly sampled variants                 ●     ○     ●     ◐      ○      ○      ●      ○
  Pairwise and higher-order variants        ●     ○     ●     ○      ◐      ●      ●      ◐
  Model-guided variants                     ◐     ○     ○     ○      ○      ○      ●      ●

VARIANT TYPES
  Codon / amino-acid substitutions          ●     ●     ○     ○      ○      ◐      ●      ○
  Insertions and deletions                  ●     ●     ◐     ●      ○      ○      ○      ●
  Combinatorial multi-motif placement       ●     ○     ●     ○      ○      ○      ○      ●
  Recombination                             ●     ○     ○     ○      ○      ○      ○      ○
  Shuffling                                 ●     ○     ◐     ○      ○      ○      ○      ●

PHYSICAL CONSTRUCTION
  Synthesis-constraint checking             ◐     ◐     ◐     ◐      ●      ○      ●      ○
  Constraint-based optimisation             ○     ○     ◐     ◐      ●      ●      ●      ○
  Primer design                             ○     ○     ○     ○      ●      ●      ◐      ○

METADATA AND INSPECTION
  Per-sequence design cards                 ●     ●     ◐     ●      ●      ●      ●      ○
  Automatic naming                          ●     ●     ●     ◐      ◐      ●      ○      ◐
  Sequence styling                          ●     ○     ◐     ○      ○      ○      ○      ○

GENOMIC INTEGRATION
  Genome coordinates                        ◐     ●     ●     ◐      ○      ○      ○      ◐
  Transcript / annotation aware             ○     ●     ○     ○      ○      ○      ○      ○
```

**Key:** ● supported · ◐ partially supported · ○ not supported
`†` order of design operations is fixed in source

**All 20 rows are scored.** This grid is the display layout; `MATRIX.md` is
authoritative for cell values and carries the totals, the balance analysis and the
correction provenance. If the two ever disagree, `MATRIX.md` wins.

---

# Row definitions

**Moved.** The operational tests live in `ROW_DEFINITIONS.md`, which is the single
authoritative instrument and the only version given to scoring agents. This file
previously carried a second copy covering 17 of the 20 rows, annotated with
per-row scoring status (*verified*, *NOT YET SCORED*, *rescore required*). Those
labels described the mid-scoring state and were all obsolete once scoring
completed; the copy has been removed rather than re-synchronised, so there is one
place a definition can be changed.

## Evidence status — all 20 rows scored

All 160 cells are scored and traceable to a per-row audit in `scored/`. Scoring
used **one agent per row across all eight tools**, rather than one agent per tool,
which is what removed the rater variance that made the earlier per-tool values
untrustworthy.

Nine cells were corrected after their agent's first output; **seven of the nine
moved against PoolParty**. Provenance for each is in `MATRIX.md`.

## Open risks

**`Recombination` is a one-tool row — RESOLVED 2026-08-17: kept, on currency.**
PoolParty scores ● and all seven others ○, so the row has zero discrimination and
is the likeliest target of a padding criticism. It stays. The defence, which the
Results text must carry rather than leave implicit:

1. **The table is visibly not stacked.** PoolParty is *not* best on six of twenty
   rows, two of them to VaLiAnT. Sole-● on three. A padded table does not concede
   six rows.
2. **The exclusion is principled, not gerrymandered.** The attack is that row 10's
   "concatenating user-supplied fragments does not count" was worded to exclude
   competitors. It was not — it is the global rule applied uniformly, and the same
   rule cost PoolParty cells elsewhere. **Seven of the nine scoring corrections
   moved against PoolParty**; that ratio is the evidence the instrument was not
   tuned, and it belongs in the response letter.
3. **Currency comes from the paper's third declared application.** Recombination is
   marginal for wet-lab MPRA but standard for in-silico attribution — swapping
   segments between parents to ask which region drives a model's prediction.

**Required in the Results:** state plainly that recombination is the one capability
no compared tool provides, and why it is included. A referee who notices it first
reads it as padding; naming it first costs a clause and removes the objection.

**Two rows still carry medium-confidence competitor cells** —
`randomly_sampled_variants`/DNA Chisel and `composable_operations`/Oligopool were
reviewed and stand, but rest on documentation rather than executed evidence. Both
are in the two tools that could not be installed and run. See the contested-cell
review in `MATRIX.md`.
