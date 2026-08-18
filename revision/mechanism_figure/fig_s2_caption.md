# Supplementary Figure S2 — caption

Paste-ready for `main.tex`. 115 words.

```latex
\caption{\textbf{A worked example of sequence generation.}
\textbf{(A)} Python code defining the library.
\textbf{(B)} Number of states of each Pool and Operation.
\textbf{(C)} Generating the sequence for $s_{Final} = 9$. Left: state assignment (blue).
Right: sequence construction (magenta). The lists beside Pool A and \texttt{mutagenize} give
the sequences and mutations indexed by their states. Inactive branches are shown in grey.
\textbf{(D)} Style and name for the same sequence. Each generated sequence is a \texttt{Seq}
object that bundles the sequence string with its \texttt{.style}. The \texttt{style} and
\texttt{prefix} arguments are supplied to Operations; the resulting \texttt{.style} and name
are held by Pools. Sequence characters are colored by their \texttt{.style} entries; where
entries overlap, the last one applies.}
\label{fig:figureS2}
```

## What the caption carries, and why only this

The figure is self-annotating, so the caption states only what appears nowhere else.

Already **in the figure** — not repeated in the caption:

| | Where |
|---|---|
| "Chaining an Operation multiplies the state count; stack adds it." | panel B footnote, verbatim |
| "S = a number of states; s = one particular state." | panel B footnote, verbatim |
| "records active branch, not counted" | panel B, stack row |
| `generate_library` sets s_Final = 9 | panel C, step 1 |
| `stack` / `mutagenize` decompose, with the if/else and mod/div rules | panel C, steps 2-3 |
| Pool A emits TGC; mutagenize applies the mutation; stack returns the active branch | panel C, steps 4-6 |
| `s_C = None`, `inactive` | panel C, Pool C nodes |
| branch structure of the DAG | panel A, code comments |
| `style=`, `.style`, `prefix=`, `contributes`, `name` row labels | panel D gutter |

Already **in the main text** — also not repeated:

| | Where |
|---|---|
| state decomposes into an internal state plus component states for upstream Pools | Sequence generation, verbatim |
| Cartesian product / disjoint union / mixed-radix decomposition | StateTracker |
| each Operation constructs a sequence from its internal state and upstream sequences | Sequence generation |
| sequences are automatically assigned an informative name | Background |
| the `style` parameter passed to `from_seq` | Example 2 |

So the caption carries exactly six things, each of which is otherwise undefined:

1. blue = state assignment, magenta = sequence construction, grey = inactive
2. what the unlabelled lists beside Pool A and `mutagenize` are
3. that a generated sequence is a `Seq` object bundling the string with `.style`
4. `style=` / `prefix=` are arguments on Operations
5. `.style` / name are results held by Pools
6. the render precedence rule (last overlapping entry wins)

Item 3 closes a real gap: `.style` appears in panel D as a labelled row, but the main text
never names an object or an attribute — Sequence metadata only says styles "combine in a
natural manner as sequences propagate through the DAG".

## Title scope

"Sequence generation" covers all four panels, justified by the main text: naming happens
"when sequences are generated" (Background) and styles "combine ... as sequences propagate
through the DAG" (Sequence metadata). The two-step definition in Sequence generation scopes
how the sequence string is built, not what generation as a whole produces.

## Fits the BMC height budget

BMC scales full-width figures to 170 mm and caps figure + legend at 225 mm.

| | |
|---|---|
| PDF as exported | 180.6 x 184.5 mm, cropped to content, fonts embedded, fully vector |
| scale to 170 mm wide | x 0.9413 |
| printed figure height | 173.7 mm |
| legend, 115 words | ~7.6 lines ~ 27 mm |
| **total** | **~200 mm** (cap 225) |
| headroom | ~25 mm |

## Suggested main-text change

Cite this figure from **StateTracker**, after the mixed-radix sentence. That paragraph states
the rules abstractly ("Cartesian product", "disjoint union", "mixed-radix decomposition") and
currently cites no figure; panel C is the worked instance of it. A reader meeting the abstract
term there is exactly the reader who wants the example.

This is better than naming mixed-radix in the caption: the reader looking at `mod` / `div` in
the figure does not need the term, while the reader hitting the term in the text does need the
example. Captions elsewhere in this paper never cross-reference sections.

## Still open on the figure itself

1. **Body text prints at 5.4-5.6 pt** at 170 mm - about 1,570 of ~1,660 characters. BMC
   requires text legible at final size. Fix: raise 8 px fonts to 10-11 px and spend some of
   the 25 mm headroom.
2. **221 vector paths at 0.01-0.02 pt**, against BMC's 0.25 pt floor. Likely table and box
   borders exported at near-zero width; they may drop out in print.
3. **Panel A cannot produce panel D.** `generate_library(num_cycles=1)` returns
   `['name', 'seq']` only - verified by execution. Add `_include_inline_styles=True`, or the
   caption needs a clause saying styles require it (~12 words).
4. **Grey does double duty.** The lookup lists beside Pool A and `mutagenize` are set in grey
   `#666666`, and grey also marks the inactive branch - two meanings, one colour. Setting the
   lists to black would reserve grey for "inactive" and let the panel C caption clause drop
   its workaround phrasing.
5. **Figure and text use different words for the same rules**: the figure says "multiplies" /
   "adds", StateTracker says "Cartesian product" / "disjoint union". Both correct; plain
   language in the figure is reasonable. Noted in case a reviewer asks.
