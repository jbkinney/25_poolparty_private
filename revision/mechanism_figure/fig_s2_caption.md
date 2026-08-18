# Supplementary Figure S2 — caption

Paste-ready for `main.tex`. 78 words.

```latex
\caption{\textbf{A worked example of sequence generation.}
\textbf{(A)} Python code defining the library.
\textbf{(B)} State counts.
\textbf{(C)} Generating the sequence for $s_{Final} = 9$. Left: state assignment (blue).
Right: sequence construction (magenta). Inactive branches are shown in grey.
\textbf{(D)} Style and name for the same sequence. The \texttt{style} and \texttt{prefix}
arguments are supplied to Operations; the resulting \texttt{.style} and name are held by
Pools. Sequence characters are colored by their \texttt{.style} entries; where entries
overlap, the last one applies.}
\label{fig:figureS2}
```

## Plain reading

**A worked example of sequence generation.**
**(A)** Python code defining the library.
**(B)** State counts.
**(C)** Generating the sequence for s_Final = 9. Left: state assignment (blue). Right:
sequence construction (magenta). Inactive branches are shown in grey.
**(D)** Style and name for the same sequence. The `style` and `prefix` arguments are supplied
to Operations; the resulting `.style` and name are held by Pools. Sequence characters are
colored by their `.style` entries; where entries overlap, the last one applies.

## Why the caption is this short

The figure is self-annotating. Everything below is already printed **in** the figure, so the
caption does not repeat it:

| Already in the figure | Where |
|---|---|
| "Chaining an Operation multiplies the state count; stack adds it." | panel B footnote, verbatim |
| "S = a number of states; s = one particular state." | panel B footnote, verbatim |
| "records active branch, not counted" | panel B, stack row |
| `generate_library` sets s_Final = 9 | panel C, step 1 |
| `stack` / `mutagenize` decompose ... with the mod/div and if/else rules | panel C, steps 2-3 |
| Pool A emits TGC; mutagenize applies the mutation; stack returns the active branch | panel C, steps 4-6 |
| `s_C = None` and `inactive` | panel C, Pool C nodes |
| branch structure of the DAG | panel A, code comments |
| `style=`, `.style`, `prefix=`, `contributes`, `name` row labels | panel D gutter |

Only four things were **not** in the figure, and those are what the caption carries: the
blue / magenta / grey conventions, the `.style` render precedence rule, and panel D's
argument-on-Operations / result-on-Pools row structure.

The title covers all four panels. Justified by the main text: naming happens "when sequences
are generated" (Background), and styles "combine in a natural manner as sequences propagate
through the DAG" (Sequence metadata). The two-step definition in Sequence generation scopes
how the sequence string is built, not what generation as a whole produces.

## Fits the BMC height budget

BMC scales full-width figures to 170 mm and caps figure + legend at 225 mm.

| | |
|---|---|
| PDF as exported | 180.6 x 184.5 mm, cropped to content, all fonts embedded, fully vector |
| scale to 170 mm wide | x 0.9413 |
| printed figure height | 173.7 mm |
| legend, 78 words | ~5.3 lines ~ 18 mm |
| **total** | **~192 mm** (cap 225) |
| headroom | ~33 mm |

The headroom is what makes the font-size fix below feasible.

## Still open on the figure itself

1. **Body text prints at 5.4-5.6 pt** at 170 mm - about 1,570 of ~1,660 characters. BMC
   requires text legible at final size. Fix: raise 8 px fonts to 10-11 px in the drawio and
   spend some of the 33 mm headroom.
2. **221 vector paths at 0.01-0.02 pt**, against BMC's 0.25 pt floor. Likely table and box
   borders exported at near-zero width; they may drop out in print.
3. **Panel A cannot produce panel D.** `generate_library(num_cycles=1)` returns
   `['name', 'seq']` only - verified by execution. Add `_include_inline_styles=True`, or the
   caption needs a clause saying styles require it (~12 words).

If (3) is not fixed, append to panel D: "Styles are returned by \texttt{generate\_library}
when \texttt{\_include\_inline\_styles=True}."
