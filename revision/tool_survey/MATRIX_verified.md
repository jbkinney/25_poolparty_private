# Verified capability matrix — 8 argument-critical rows, 13 tools

Row definitions: `ROWS_v3.md`. Per-row audits with quoted evidence: `verified/`.
Every cell below was scored by a single max-effort agent that owned that row
across all 13 tools, from primary sources only (memos were leads, not evidence).

`Y` = yes · `~` = partial · `.` = no · `†` order fixed in source · `‡` no design ops to compose

```
                                       PoolP VaLiA MPRAn MPRAD Oligo MutMk Codon DNACh ledid tange Biopy pydna SeqPr
library_first_class_object   (reworded)   Y     .     .     .     ~     .     .     Y     .     .     .     Y     .
user-composable design ops   (reworded)   Y     .†    ~     .†    Y     .†    .†    ~     ~     Y     .‡    .‡    .‡
unmaterialized_lib_addressing (replaced)  Y     .     .     .     .     .     .     ~     .     .     .     .     .
library_algebra                           Y     .     .     .     ~     .     .     .     .     .     .     .     .
exhaustive_single_scans                   Y     Y     ~     .     ~     ~     .     .     .     Y     .     .     .
sampled_random_mutagenesis                Y     .     Y     ~     .     ~     ~     ~     ~     ~     .     .     .
higher_order_combinatorial                Y     .     Y     .     .     Y     .     ~     ~     ~     .     .     .
heterogeneous_components_one_library      Y     Y     ~     ~     .     ~     .     .     .     .     .     .     .
```

## The central result

**PoolParty is the only tool scoring `yes` on all four mutagenesis rows.**

| Row | Tools scoring yes |
|---|---|
| `exhaustive_single_scans` | VaLiAnT, tangermeme, **PoolParty** |
| `sampled_random_mutagenesis` | MPRAnator, **PoolParty** |
| `higher_order_combinatorial` | MPRAnator, Mutation Maker, **PoolParty** |
| `heterogeneous_components_one_library` | VaLiAnT, **PoolParty** |

No other tool reaches even 3 of 4. This is a **generality** claim, not a novelty
claim: every competitor covers a different subset; only PoolParty covers the union.

It also explains the GB1 result mechanically. The GB1 library needs exhaustive
singles **+** exhaustive pairs **+** sampled higher-order **+** WT replicates
simultaneously. VaLiAnT has exhaustive + heterogeneous but not sampled or
higher-order (source-verified: `MutatorType` is a closed enum of `DEL, SNV,
SNV_RE, IN_FRAME, ALA, STOP, AA`, no stochastic member, no plugin mechanism).
MPRAnator has sampled + higher-order but is MPRA-only (`assay_dms` = no).

The result survived a pass that changed 25 cells and demoted PoolParty on one
row — which is what makes it worth publishing.

## Rows PoolParty holds alone

`unmaterialized_library_addressing` (DnaChisel partial) and `library_algebra`
(Oligopool Calculator partial). Everything else it scores `yes` on, at least one
other tool also does.

That modesty is an asset: a column of unique wins would read as marketing and
invite exactly the overstatement criticism the editor raised.

## Corrections that change the competitive picture

| Tool | Change | Consequence |
|---|---|---|
| **tangermeme** | `yes` on `user-composable design ops` and on `exhaustive_single_scans` | A more serious competitor than the Background's "general-purpose toolkit" framing implies. `saturation_mutagenesis` does the positional loop itself over the whole sequence. |
| **DNA Chisel** | `no` -> `yes` on `library_first_class_object` | The earlier pass missed `MutationSpace` entirely because DnaChisel never calls it a library. Also the strictest competitor on laziness. |
| **Oligopool Calculator** | `yes` -> `partial` on library object; `no` -> `partial` on `library_algebra` | `oligopool.join` is a documented Design Mode op over two variant tables. Still uncited in the manuscript. |
| **Mutation Maker** | `yes` -> `partial` on `exhaustive_single_scans` | All three workflows require a user-supplied `mutations` list, so exhaustiveness is the user's job. |
| **MPRA Design Tools** | `no` -> `partial` on two rows | One `processVCF()` call dispatches SNV/insertion/deletion to three distinct construct generators. |

## Rulings recorded for the response letter

Referees may contest these; each was pre-registered before scoring and applied
identically to PoolParty and every competitor.

1. **Combinatorial motif placement does not count as higher-order mutagenesis** —
   a separate `combinatorial_motif_place` row scores it. Verified to change no
   cell, since PoolParty and MPRAnator both qualify on strictly mutational engines.
2. **IUPAC-degenerate template expansion does not count as higher-order** — the
   user hand-authors the combinatorial spec and the tool decompresses it. Forced
   by consistency: the identical product-over-IUPAC-positions function exists in
   both Oligopool (`core_degenerate.py:858`) and pydna (`primer_screen.py:178`).
3. **Objective-driven optimisers that emit multi-edit sequences without a
   declarable combinatorial space score `partial`, not `yes`** (DnaChisel, Ledidi,
   tangermeme).
4. **A capability existing only in an example script caps at `partial`**, never `yes`.
5. **Shared-node memoisation was excluded** from the composability threshold as
   self-serving to PoolParty.

## Lowest-confidence cell

`MPRAnator` on `user-composable design ops` = `partial` (medium confidence). It
rests on user-orderable sequence-element slots (`part1.html:81` ->
`part1.py:241`) and an undocumented but mechanically valid FASTA-out/FASTA-in
re-feed (`oligo.py:74` vs `mycustom.py:32-33`). Strictly above VaLiAnT on both
counts, so anchor ordering holds, but a referee who requires the re-feed to be
*documented* would move it to `no†`. Authors' choice.

## Scope limit

Only these 8 rows are verified. The other 20 rows (Blocks B–E) still carry
per-tool-agent values with known rater variance. Audit any row the same way
before putting it in the main-text table.
