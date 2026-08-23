# `code_comparison/` — expressing other tools' libraries in PoolParty

Answers **R2 2b** (output overlap with existing tools) and supports **R1 #3/#5**
and **R3**. Sibling of `../comparison/`, which builds Tables 1 and 2, and of
`../benchmarks/`, which answers the runtime and memory asks.

One directory per tool. Each holds a `comparison.md` (the finding and how it was
verified), a `figure.md` (the design record, and the source of the caption text),
and a `repro.py` that executes the comparison.

| Directory | Tool | Status |
|---|---|---|
| `valiant/` | VaLiAnT | Verified — 357/357 sequences including adaptors. Figure specified, not drawn |
| `tangermeme/` | tangermeme | Verified — 200,000 rows, identical set and order against the published SpliceAI GT-arm design cards. Figure specified, not drawn |
| `mpranator/` | MPRAnator | Not comparable — the design endpoint returns HTTP 500 on its own documented use case. No figure and no numerical claim; see that directory's `comparison.md` |

## Shared drawing spec

Both figures are drawn to match Supplementary Figure S2
(`../mechanism_figure/fig_s2.drawio`), the only figure in this revision that
already contains a code block. Reuse its vocabulary rather than inventing a
second one.

**Upright, full width, wide aspect — not a rotated page.** Every figure in the
submission is a plain `\includegraphics[width=\textwidth]` (`main.tex` lines 100,
170, 189, 209); the manuscript contains no sideways float. A side-by-side code
comparison does not need one, and an upright figure spares the reader turning the
page.

**Width budget.** `../comparison/figure_s2.tex` records that BMC scales
full-width figures to 170 mm = 482 pt. Figure S2 fits 138 characters of Consolas
across that width — its content bounding box is 710 units at 5.16 units per
character, exported 1:1 to a 511.9 x 523.0 pt PDF.

| | Characters | pt/char | Effective type |
|---|---|---|---|
| Figure S2, as drawn | 138 | 3.49 | ~6.4 pt |
| Full width, our figures | 130 | 3.71 | ~6.7 pt |
| One of two side-by-side columns | 64 | 3.71 | ~6.7 pt |

Both figures therefore sit at slightly larger type than the figure already
accepted into the revision. Current layouts are 130 characters (`tangermeme/`)
and 83 (`valiant/`), so neither needs re-flowing.

**Code block cell style**, copied from `fig_s2.drawio`:

```
text;html=1;strokeColor=none;fillColor=none;align=left;verticalAlign=top;
whiteSpace=wrap;rounded=0;fontSize=9;fontFamily=Consolas;spacing=5;
spacingLeft=8;spacingTop=6;
```

Consolas at `font-size: 8.3px` inline, one `<div>` per line, `<br>` for blank
lines.

**Syntax palette**, verified against the tokens each colour actually wraps in
`fig_s2.drawio`:

| Colour | Role |
|---|---|
| `rgb(0, 102, 204)`, bold | keywords — `import`, `as` |
| `rgb(106, 115, 125)` | comments |
| `rgb(163, 21, 21)` | string literals |

`fig_s2.drawio` is plain uncompressed XML, so it can be edited by script.

Figure S2's own code block ends `df = Final.generate_library(num_cycles=1)`, which
is why both comparison figures use `generate_library` rather than `to_df`.

## Open

- **170 mm vs 552 pt.** `figure_s2.tex` gives 170 mm as the production width for a
  full-width figure. An earlier version of `valiant/figure.md` assumed a landscape
  page at 552 pt. The character budget is nearly identical either way, so no layout
  depends on the answer, but the number should be confirmed against BMC's author
  instructions before the figures are drawn.
- **Neither figure is drawn yet.** Both `figure.md` files are specifications.
