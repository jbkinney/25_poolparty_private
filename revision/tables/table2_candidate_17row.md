# Table 2 — LOCKED row set, 20 rows x 8 tools

> **STATUS: row set locked 2026-08-14.** This is the baseline against which
> Policy B filtering, reconciliation, and the verification pass are all scoped.
> Rows may still be *substituted* — see the protocol below — but the set is
> otherwise closed, and the count stays at 20.

Parallel draft. Does **not** replace `table2_capabilities.md` (8 rows x 8 columns),
which remains the shipping version until this one is verified.

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
  Synthesis-constraint checking             ?     ?     ?     ?      ?      ?      ?      ?
  Constraint-based optimisation             ○     ○     ◐     ◐      ●      ●      ●      ○
  Primer design                             ?     ?     ?     ?      ?      ?      ?      ?

METADATA AND INSPECTION
  Per-sequence design cards                 ●     ●     ◐     ●      ●      ●      ●      ○
  Automatic naming                          ●     ●     ●     ◐      ◐      ●      ○      ◐
  Sequence styling                          ●     ○     ◐     ○      ○      ○      ○      ○

GENOMIC INTEGRATION
  Genome coordinates                        ◐     ●     ●     ◐      ○      ○      ○      ◐
  Transcript / annotation aware             ○     ●     ○     ○      ○      ○      ○      ○
```

**Key:** ● supported · ◐ partially supported · ○ not supported · ? **not yet scored**
`†` order of design operations is fixed in source

---

# Row definitions

**Binding on any agent scoring these rows.** These are the operational tests, not
descriptions. Where a row was verified under an earlier wording, that wording is
reproduced verbatim so scoring stays comparable.

### 1. Library object — *verified*
A documented API call returns a value representing the whole multi-sequence library
as an instance of a type the tool defines for that purpose (carrying member
identity, design provenance, or membership operations), which the user can bind to
a variable and inspect — as opposed to the tool only writing files, or handing back
a general-purpose container (ndarray, DataFrame, ragged array, torch tensor,
list/dict of sequences) whose library meaning lives only in the caller's head.
`partial` = a durable, documented, in-memory value representing the **whole
library** exists and is **returned by the tool's own design API**, but its type is
general-purpose (DataFrame, list, tensor). A collection the **caller assembles** —
appending results to a list across repeated calls — is `no`, not `partial`.
An object representing a single sequence or a single optimisation problem is also
`no`, however purpose-built its type.

*Tightened after the first scoring batch: the previous wording ("an in-memory
library value exists but its type is general-purpose") had no floor — almost any
tool returns something in memory — and five of eight cells moved, four of them in
the same direction. The two added requirements codify the distinction the batch-1
scorer already applied when arguing DNA Chisel down to `no`.*

### 2. Composable operations — *verified*
Two or more distinct **library-design** operations share a carrier type and compose
in either order, with fan-out (one intermediate feeding two downstream operations)
expressible. General sequence utilities (parse, revcomp, translate, pad, one-hot,
write) and cloning-simulation primitives are **not** design operations.
`partial` = design operations exist and chain, but only in a tool-fixed order, only
for a small documented subset, or only by re-feeding output artefacts.
`no†` = the tool fixes the order of its design operations in source.
**Shared-node memoisation is deliberately excluded from the threshold** as
self-serving to PoolParty; the row measures expressibility only.

### 3. Saturation mutagenesis — *verified*
From one specification, the tool enumerates every single substitution at every
eligible position of a region. `partial` = exhaustive over positions but only for
one narrow event type, or exhaustiveness is the user's responsibility, or the
capability exists only in an example script.

### 4. Randomly sampled variants — *verified*
A documented rate, count or RNG parameter produces a *sample* of the variant space.
**Random sequence generation is not sampled mutagenesis** — generating random
barcodes, shuffled controls or stochastic optimisation does not count.

### 5. Pairwise and higher-order variants — *verified*
Two or more **independent** mutation events co-occurring in a single output
sequence, with the combinations enumerated or sampled **by the tool**.

Explicitly excluded: hand-authored multi-variant input (a user-written VCF);
uniform context edits applied to every oligo (e.g. PAM-protection edits); one
multi-base single event (a 3 bp deletion is one event); combinations the user
produces by running the tool repeatedly and merging.

**Boundary rules — pre-registered, and binding.** This row sits next to two others
that describe adjacent capabilities. Without these rules the same capability gets
credited twice, once here and once in its own row.

- **vs. row 9 (Combinatorial multi-motif placement): inserting multiple motifs is
  *placement*, not mutagenesis.** A tool that positions two or more motifs across a
  background sequence scores in row 9 and **not** here, however many combinations it
  produces. Score here only if the tool combines *mutation events* — substitutions,
  insertions, deletions applied to a parent sequence.
- **vs. degenerate templates: IUPAC expansion does not count.** Expanding
  `ATG NNK GGC` into its concrete sequences yields multi-position variation, but the
  user authored the combinatorial specification and the tool merely decompressed it.
  The tool must generate the combinations from a design specification, not enumerate
  a template the user already wrote.
- **vs. row 12 (Model-guided variants): objective-driven optimisers score `partial`
  at most.** A tool that emits a sequence carrying several edits because an
  optimiser converged there has not enumerated or sampled a *declarable*
  combinatorial space. Score such tools `partial` here and full credit in row 12.

*Provenance: these three were pre-registered before the earlier scoring pass and are
carried forward verbatim. They matter more in this row set than they did in the
previous one, because rows 9 and 12 now ship — previously neither was in the table,
so there was no adjacent row to double-count into.*

### 6. Mixed variant types in one library — *verified*
Two or more structurally different component types declared in **one**
specification and emitted as **one** pooled output. `partial` = achievable only by
separate runs the user merges, or limited to two component types of the same kind.

### 7. Codon / amino-acid substitutions — **NOT YET SCORED**
Substitution specified at the codon or amino-acid level, with the tool handling the
codon↔residue mapping and offering a choice of replacement policy (e.g. all
residues, most-frequent codon, degenerate codon). Nucleotide substitution that
merely happens to fall inside a coding region does **not** count — the tool must be
codon-aware. `partial` = codon awareness exists but the substitution set is fixed
or the user must supply the residue list.

### 8. Insertions and deletions — **NOT YET SCORED**
Generation of variants that insert or remove bases as a **designed variant type**,
via a supported operation. Supplying pre-made indel sequences as input does not
count. Gap-marked deletions count provided the tool either shortens the sequence or
explicitly marks the deletion in the output. `partial` = one of insertion or
deletion only, or fixed length only.

### 9. Combinatorial multi-motif placement — *v2, unverified*
Two or more distinct motifs placed across multiple positions and/or orientations,
with the combinations enumerated or sampled by the tool. `partial` = single-motif
placement, or positions enumerated but orientations fixed.

### 10. Recombination — **NOT YET SCORED**
*Gloss for the table: segments joined from two or more parent sequences.*

Construction of variants by joining segments from two or more parent sequences at
breakpoints the tool chooses or enumerates. Concatenating user-supplied fragments
does not count; the tool must generate the breakpoint combinations.

*Named "Recombination", not "chimeras" — the latter is protein-engineering
register, and this row must also cover the in-silico sense (swapping segments
between parents to ask which region drives a model's prediction).*

### 11. Constraint-based sequence optimisation — **NOT YET SCORED**
The tool **modifies** a sequence so it satisfies declared constraints, rather than
only rejecting sequences that violate them. **Rejection-only filtering does not
count** — this row is about the tool altering the sequence to comply. `partial` =
optimisation against a narrow constraint class only.
*Note:* the nearest existing key, `synthesis_constraints`, measures checking rather
than repair. Its values must **not** be lifted into this row.

### 12. Model-guided variants — *v2, unverified*
Design driven by the output of a predictive model, where the model's prediction
determines the edit. `partial` = the tool accepts an arbitrary scoring callable so
a model can be attached, but performs no optimisation against it. `yes` = an
optimisation loop against the model's output.

### 13. Per-sequence design cards — *v2, unverified*
Structured, machine-readable metadata attached to each output sequence recording
**how it was constructed**, beyond naming the mutation it carries. `partial` = some
provenance recorded, but not per sequence, not structured, or limited to the
variant identity.

### 14. Automatic naming — *v2, unverified*
Informative identifiers generated **by the tool** for each output sequence,
composed from the construction history. **Carrying a user-supplied identifier
through to the output does not count.** `partial` = names generated but not
composed from construction history, or only for some operations.

### 15. Sequence styling — *v2 values, DEFINITION NARROWED — rescore required*
Visual marking of **what changed inside the sequence itself** — highlighting,
casing, colour — so a user can audit a variant by eye. Plots, graph views and
report documents are **not** this row.
*Was `design_visualization`, which included graph and plot output. Existing v2
values were scored against the broader wording and are not valid here.*

### 16. Genome coordinates — *v2, unverified*
Accepts reference-genome coordinates as input, and/or emits coordinates that locate
each designed sequence in a reference. `partial` = coordinates accepted for
extraction only, with no coordinate output for designed variants.

### 17. Transcript / annotation aware — *v2 values, DEFINITION BROADENED — rescore required*
Represents transcript models, exon/intron structure, or annotation files (GTF/GFF),
such that design respects them.
*Was `transcript_models`, which asked only about GTF/GFF transcript models.
Broadening it to include exon/intron structure absorbs the former
`exon_intron_split_codons` row; existing values need re-checking against the wider
wording.*

---

## Evidence status — 6 of 17 rows are publishable today

| Status | Rows | Count |
|---|---|---|
| Cross-tool verified | 1, 2, 3, 4, 5, 6 | 6 |
| v2 per-tool values, unverified | 9, 12, 13, 14, 16 | 5 |
| v2 values but definition changed — rescore | 15, 17 | 2 |
| Never scored | 7, 8, 10, 11 | 4 |

**Work to ship this table:** 4 rows scored from scratch (32 cells), 2 rows rescored
under new definitions (16 cells), 5 rows verified cross-tool (40 cells). The pass
that verified the existing six changed 25 of 104 cells and found three rows
defective as worded, so none of the unverified values should be trusted.

## Open risks

**`Recombination and chimeras` may be a one-tool row.** The lexical probe suggests
PoolParty alone. If scoring confirms that, it is a row only we score — the padding
criticism, inside the block designed to avoid it. Decide in advance whether it
survives that outcome.

**Four `?` rows are guesses about what will discriminate.** They were chosen from a
word-frequency probe of the records, not from assessment. Any of them could turn
out uniform and have to be dropped.
