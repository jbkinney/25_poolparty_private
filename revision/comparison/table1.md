# Table 1 — tool overview (main text)

**Rendered by `table1.tex`.** That file is the deliverable; this one records the
design decisions behind it. Keep them in step.

**Shape follows surveyed convention** (see `CONVENTIONS.md`): tools are rows in
4/4 surveyed tables; column headers are full words with no invented acronyms; the
dominant descriptive pattern is *Tool · Purpose · Key features · Availability*.
Two columns are added beyond that pattern — *Output*, because the editor
explicitly asked for "similarities and differences in program outputs", and
*Reference*, which is near-universal.

Terminology follows `TERMINOLOGY.md`.

## Structure — revised 2026-08-17

**Rows now match Table 2's columns, one-to-one and in the same order.** DNA Chisel
and Mutation Maker were previously folded into a grouped *Adjacent design tools*
row while appearing as individual columns in Table 2. A reader could not map one
table onto the other. Both are now individual rows.

The residue of that group is CodonGenie and ledidi, renamed **Single-sequence
design tools**, which describes what they have in common rather than their
distance from us.

Ten rows: the eight scored tools, then two grouped rows.

**Order is alphabetical in both tables**, with PoolParty bold rather than first —
the surveyed DMS review orders chronologically so its own subject is not
privileged, and the same reasoning applies here.

**Key features are held to 2-3 clauses per tool.** PoolParty's cell previously ran
to five clauses against two for everyone else. Whatever its accuracy, a cell three
times longer than its neighbours reads as self-promotion.

## Layout

`sidewaystable` (landscape). Measured, not assumed: in portrait the float exceeds
the page by 99pt in referee mode and 83pt in production mode, in both cases
because the prose columns wrap heavily at 13.1cm. Landscape gives 19.4cm and the
float then fits with no overfull boxes in either mode.

Requires `booktabs`, `array`, `rotating`. Springer's own template loads booktabs
and multirow, so this is unremarkable, but all three are additions to `main.tex`.

A `\singlespacing` override was tried to defeat referee-mode double spacing and
**measured to have no effect inside the float** — it made the table 8pt taller.
Removed rather than left in as decoration.

## The table

**Not duplicated here.** The rendered table is `table1.tex`; keeping a second copy
in markdown is how the two fell out of step before (the markdown still showed
`Adjacent design tools` after Table 2 had been scoring DNA Chisel and Mutation
Maker as separate columns for days). Read the cell contents there.

Row order, matching Table 2: DNA Chisel · MPRA Design Tools · MPRAnator ·
Mutation Maker · Oligopool Calculator · **PoolParty** · tangermeme · VaLiAnT ·
General-purpose toolkits *(Biopython, pydna, SeqPro)* · Single-sequence design
tools *(CodonGenie, ledidi)*.

---

## Design decisions

**Row order is alphabetical within the comparable set, with PoolParty in place
rather than first.** The DMS review orders chronologically so its own subject is
not privileged. Bolding identifies PoolParty without promoting it. Trivially
changed if you prefer it first.

**Output earns its column.** It is the most direct response to the editor's ask
about program outputs, and reading down it makes the complementarity argument
structurally: VaLiAnT emits VCF with consequence annotation, Oligopool Calculator
emits a synthesis-ready library from variants a user already has, DNA Chisel and
the single-sequence tools emit one sequence. Nobody has to be told these are
different problems.

**No Pros/Cons column**, though two of the seven CRISPR-review tables use one.
Writing "Cons" about tools whose authors may referee this paper is precisely the
overstatement risk the editor flagged. Purpose and Output let the reader conclude.

**No genomic-context column.** It dissolves into VaLiAnT's Purpose
("from genomic coordinates", "transcript-aware") and Output ("VCF carrying variant
consequence annotation"), and into PoolParty's Key features, which names its input
types without claiming genomic integration. R3's concession lands without a
bespoke column a referee would notice we invented.

**No version, license or documentation-link columns.** Absent from all four
surveyed tables. Version is handled by the caption footnote; license belongs in
the supplement.

## On the two grouped rows

Convention lists every tool individually (4/4 surveyed). We deviate for the two
final rows only, because the table spans tool *categories*, which none of the
surveyed tables does — the CRISPR systematic review, which faces the same problem,
answers it with seven separate tables, which is worse here.

The deviation is now much smaller than it was: every tool scored in Table 2 has
its own row, so grouping applies only to five tools that Table 2 does not score.

Three mitigations: members are named in the Tool cell, all members are cited in
Reference, and the footnote states the grouping rationale.

**Known instability:** tangermeme was moved out of the general-purpose group after
the audit found it satisfies saturation mutagenesis and composable design
operations. That grouping is a judgement call, not a fact, and the caption should
own it. It also means the Background sentence claiming general-purpose toolkits
"have no notion of a variant library as a structured object" needs narrowing.

## Open decision, unchanged

**tangermeme's row.** Given its own row here. Alternative: return it to the group
and narrow the Background sentence instead.
