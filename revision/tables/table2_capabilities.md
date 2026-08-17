# Table 2 draft v5 — library-design capabilities (main text)

Matrix style, **pivoted**: capabilities as rows, tools as columns. Rows carry a
short gloss so the operational definitions sit in the table rather than in a
caption or a separate file.

**Columns are exactly Table 1's tool column**, same eight entries, same order,
same spelling — including the two grouped entries. A reader can map between the
tables by position.

Every cell comes from the verified row audits (`tool_survey/verified/`), where one
max-effort agent owned each capability across all 13 tools and scored from primary
sources only. These eight rows are the only ones audited to that standard.

---

## The table

```
                                         MPRA     MPRA   Pool           Oligopool           General-  Adjacent
                                         Design   nator  Party  VaLiAnT Calculator tanger-  purpose   design
                                         Tools                                     meme     toolkits  tools
──────────────────────────────────────────────────────────────────────────────────────────────────────────────
Library specification
  Library object                            ○        ○      ●       ○       ◐         ○        ● a       ● b
    held, inspected, passed on
  Composable operations                     ○        ◐      ●       ○       ●         ●        ○         ◐ c
    chained and branched by the user
  Library algebra                           ○        ○      ●       ○       ◐         ○        ○         ○
    libraries combined, sampled, repeated
  On-demand generation                      ○        ○      ●       ○       ○         ○        ○         ◐ b
    any slice, without materialising

Variant generation
  Saturation mutagenesis                    ○        ◐      ●       ●       ◐         ●        ○         ◐ d
    every substitution, every position
  Randomly sampled variants                 ◐        ●      ●       ○       ○         ◐        ○         ◐
    drawn by rate or count
  Pairwise and higher-order variants        ○        ●      ●       ○       ○         ◐        ○         ● d
    two or more mutations per sequence
  Mixed variant types in one library        ◐        ◐      ●       ●       ○         ○        ○         ◐ d
    different components, one specification
```

**Key:** ● supported · ◐ partially supported · ○ not supported.
*Tools examined August 2026.*

**Grouped columns** show the highest level any member reaches; the footnote names
which member. Members not named are ○.

- **a** — pydna. Biopython and SeqPro ○.
- **b** — DNA Chisel. CodonGenie, ledidi and Mutation Maker ○.
- **c** — DNA Chisel and ledidi. CodonGenie and Mutation Maker ○.
- **d** — Mutation Maker. On *pairwise and higher-order variants*, DNA Chisel and ledidi are ◐ and CodonGenie ○.

Cell distribution: 18 ● · 15 ◐ · 31 ○ (64 cells).

Group members in full: **General-purpose toolkits** = Biopython, pydna, SeqPro.
**Adjacent design tools** = CodonGenie, DNA Chisel, ledidi, Mutation Maker.

---

## Consistency with Table 1: now exact

All eight column names are byte-identical to Table 1's tool column, in the same
order. Nothing appears in one table and not the other.

The previous draft listed 13 individual tools, which put CodonGenie, DNA Chisel,
ledidi, Mutation Maker, Biopython, pydna and SeqPro in Table 2 with no
corresponding Table 1 row, and abbreviated several names incorrectly in the header
(`DNA Chisl`, `Mut Makr`, `Bio-pyth`, `Codon Genie`, `Seq-Pro`). Both problems are
gone: eight columns, full names, no truncation.

**Cost of grouping, stated plainly.** Per-tool resolution inside the two grouped
columns is lost from the face of the table and survives only in the footnotes. The
supplementary matrix carries all 13 tools individually, so nothing is
unrecoverable, but a reader of the main text cannot see that DNA Chisel reaches ●
on `Library object` without reading footnote b.

## What the grouped columns still preserve

The evidence that matters most survives. On `Library object`, three columns reach ●
— PoolParty, General-purpose toolkits (pydna) and Adjacent design tools (DNA
Chisel). On `Composable operations`, three reach ● — PoolParty, Oligopool
Calculator and tangermeme. PoolParty is visibly not alone on either row, which is
the strongest available evidence the table is not rigged.

## Rendering: check this compiles before the table depends on it

**Harvey balls are not core LaTeX.** The usual route is:

```latex
\usepackage{wasysym}     % \CIRCLE  \LEFTcircle  \Circle
```

Three known risks, in order of likelihood:

1. **Springer `sn-jnl` package conflicts.** The class is fussy about added
   packages. `wasysym` redefines several symbols (notably `\Box`, `\diamond`,
   `\iint`) and can clash with `amssymb`, which `main.tex` already loads. If both
   are needed: `\usepackage[nointegrals]{wasysym}`.
2. **BMC production re-typesetting.** BMC may re-set the article from source; a
   symbol that renders in our PDF is not guaranteed to survive their pipeline.
   Keep the legend spelled out in words so the table stays readable if symbols
   degrade.
3. **Half-fill orientation.** `\LEFTcircle` fills the left half. It does not matter
   semantically, but be consistent and do not mix with any other half-symbol in
   the paper.

**Fallbacks, in preference order:**

| Fallback | Gives | Cost |
|---|---|---|
| `\usepackage{pifont}` -> `\ding{108}` `\ding{109}` | solid and open circles only | no half state; would force resolving 15 partials |
| TikZ-drawn balls | full control, no package clash | verbose; ~6 lines of preamble macros |
| `\usepackage{harveyballs}` | purpose-built | less widely installed; may not be on the journal's TeX system |
| ✓ / ~ / ✗ via `amssymb` | no new package | loses the ordinal reading; ✗ reads as a verdict against other people's software |

At eight columns the headers fit horizontally with two- or three-line wrapping, so
`\rotatebox` is not required — one advantage of matching Table 1 rather than
listing all 13 tools.

**Action before submission:** compile `main.tex` with the chosen package and check
the symbols at final print size, not just on screen.

## Why circles rather than ✓ / ~ / ✗

- **The data is ordinal**, so the encoding should be. ● ◐ ○ are one object at three fill levels and read as a scale without the legend; ✓ / ~ / ✗ are unrelated glyphs, and `~` does not sit intuitively between a check and a cross.
- **○ is neutral where ✗ is a verdict.** Thirty-one cells here are "not supported". Empty circles read as absence; crosses read as an attack on other people's software, in a table whose referees may be its authors.

**Limitation, stated honestly:** all four tables in `TABLE_CONVENTIONS.md` use
binary encodings, so they show checkmarks suit *binary* data, not that they suit
ours. Harvey balls are the general convention for graded comparison, but no
in-domain graded example has been located.

## The "not applicable" level is not used

The legend previously reserved `—` for capabilities that do not apply to a class of
tool. There are **no such cells**: the audit returned a real value for every tool
on every row. Biopython genuinely does not provide saturation mutagenesis, which is
`○`. Three symbols ship; `—` is retired.

## Why pivoted

Full-word capability names as column headers would need three-line stacked headers
or the acronyms we ruled out. Pivoting is what makes the full-word decision
affordable, and it creates room for the glosses.

**Logged deviation:** all four surveyed tables put tools in rows, including the
codon-optimization feature matrix — which affords criteria-as-columns only because
its criteria are short acronyms (CAI, ICU, CC). Recorded in
`TABLE_CONVENTIONS.md`.

## The result this table carries

**PoolParty is the only tool satisfying all four variant-generation capabilities.
No other tool satisfies three.** Each competitor covers a different subset:

- **VaLiAnT** — saturation mutagenesis and mixed types; no random sampling, no higher-order variants. Source-verified: `MutatorType` is a closed enum of seven members, no stochastic element, no combination mutator.
- **MPRAnator** — random sampling and higher-order variants, but MPRA only.
- **Mutation Maker** (in Adjacent design tools) — higher-order variants, but the user supplies the mutation list.

A generality claim, not a novelty claim. It also explains the GB1 result
mechanically: that library needs all four at once.

## PoolParty's column is all-●

That is what the verified data says, and it survived a pass that changed 25 cells
and demoted PoolParty on one row before that row was replaced. But it is the visual
most likely to draw the overstatement criticism. Mitigations, all already true:

1. Table 1 shows plainly where PoolParty does not compete — no genomic coordinates, no VCF, no consequence annotation, no primer design.
2. The caption should state these axes were chosen because the paper's three examples exercise them, not chosen post hoc.
3. The supplement carries the full 28-row matrix over all 13 tools individually, where PoolParty is not uniformly best.

If it still reads too strong, the defensible trim is dropping `Library algebra` and
`On-demand generation` — PoolParty is sole ● on both, they are the least intuitive
of the eight, and they survive in the supplement.

## Flagged deviation: three symbol levels

Both surveyed matrix tables use two (✓ / blank). Dropping ◐ is not viable:
**15 of 64 cells (23%) are partial.** Promoting them to ● overstates competitors
and dilutes our own claim; demoting to ○ understates tools whose authors may
referee, which is the more damaging direction.
