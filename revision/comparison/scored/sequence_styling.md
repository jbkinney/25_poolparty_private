# Audit: `sequence_styling`

Audit date: 2026-08-15

## Measurement instrument and operational test

I read the preamble and binding global rule in
`25_poolparty_private/revision/tables/ROW_DEFINITIONS.md:1-34` and the
`sequence_styling` definition at lines 169-172 before inspecting any tool. The
binding rule requires a tool-provided operation, parameter, or mode; styling
assembled by the caller does not count. The row is deliberately limited to
visual marking inside sequence text itself. Plots, graph/feature views, and
report documents do not count, even if they color positions or bases.

**Operational test declared before tool inspection:** For each tool, inspect its
documented API and implementation/output paths for a tool-provided operation,
parameter, or mode that renders the changed bases within sequence text itself
using highlighting, case, or colour; full documented support is `yes`, a
restricted or example-only implementation is `partial`, absence after targeted
source/docs searches is `no`, and plots, graphs, reports, mutation labels, or
user-written styling do not qualify.

The same test and exclusions were applied to all eight tools.

## Result

| Tool | Value | Confidence |
|---|---|---|
| poolparty | **yes** | high |
| valiant | **no** | high |
| mpranator | **partial** | high |
| mpradesign | **no** | high |
| oligopoolcalc | **no** | high |
| mutationmaker | **no** | high |
| dnachisel | **no** | high |
| tangermeme | **no** | high |

## Source snapshots

Repository evidence was examined at these commits:

- PoolParty local repository: `1bb0179e1c3720b1fffd471802b3040f9336de28`
- VaLiAnT `develop`: `8796cc112dafd4919fec59913f58cd2be87c45eb`
- VaLiAnT official wiki: `5f4d9023a55b2c2ee7ee65ad785865e99cc86576`
- MPRAnator: `9969790d62410138d4281b7955da6d085f07b1bc`
- mpradesigntools: `afd386ef12051bb0a37ad63a6f92acd555246584`
- designMPRA: `0cf56ee602fc86dde705906d071a39cbdf6e99a8`
- Oligopool Calculator: `b88fa394ca67ed4c48ec41127b5707694ee7cf0a`
- Mutation Maker: `396c1c0ede7529f3dedf65e821e8c1f20c9a7043`
- DnaChisel: `68c09304341c3656f3dfe63eda37757d6a7b3917`
- tangermeme: `2006b310cd72a28c56c3ea4ba67f738fff74bb89`

The seven external repositories were shallow-cloned read-only into `/tmp` for
full-tree searches. The per-tool survey records were consulted only to locate
official source URLs. Their ratings and conclusions were not used as evidence.

## Per-tool audit

### poolparty — `yes` (high confidence)

**Qualifying evidence.** The dedicated documentation states:

> “PoolParty can apply **inline styles** to sequences so that mutations,
> regions, and other features are visually highlighted when printed. Styles
> are tracked per-character through the entire DAG and rendered as ANSI escape
> codes in terminal output.”

Source: `poolparty-statecounter/poolparty/docs/metadata/styling.rst:4-7`.

The same page says, “Many operations accept a `style` parameter that
automatically highlights the modified positions” and demonstrates
`wt.mutagenize(num_mutations=1, style="red")` followed by
`muts.print_library()` (`docs/metadata/styling.rst:130-146`). The documented
`mutagenize` parameter is explicitly “Named display style applied to mutated
bases” (`docs/operations/mutagenize.rst:52-55`).

This is implemented as a direct API parameter:

> `style : Optional[str], default=None`
>
> `    Style to apply to mutated positions (e.g., 'red', 'blue bold').`

Source: `src/poolparty/base_ops/mutagenize.py:18-24,51-52`. The mutation
operation adds the style at the mutated raw positions and returns it with the
sequence (`mutagenize.py:582-596`); `print_library` applies those inline styles
to the sequence before printing (`src/poolparty/pool.py:361-376`). This is
marked-up sequence text, not a plot or report.

**Behavioral confirmation.** As required for PoolParty, I ran its local,
read-only Python 3.12 environment with bytecode disabled:

```text
PYTHONDONTWRITEBYTECODE=1 poolparty-statecounter/.venv/bin/python -c \
  'import sys; sys.path.insert(0,"poolparty-statecounter/poolparty/src"); \
   import poolparty as pp; pp.init(); \
   p=pp.from_seq("AAAA").mutagenize(num_mutations=1,mode="sequential",style="red"); \
   p.print_library(num_seqs=2)'
```

The sequence lines were:

```text
\x1b[91mC\x1b[0mAAA
\x1b[91mG\x1b[0mAAA
```

Thus the changed first base, and only that base, was wrapped in red ANSI codes.
A separate generation call returned `_inline_styles = SeqStyle(1 style,
length=4)` with the plain sequence, confirming that the operation carries the
per-character styling metadata into its renderer.

**Search and disconfirmation attempt.** I searched the complete local
`src/poolparty` and `docs` trees for
`style|highlight|color|colour|case|uppercase|lowercase|ansi|html|changed|mutation`,
then inspected `docs/metadata/styling.rst`,
`docs/operations/mutagenize.rst`, `base_ops/mutagenize.py`, `pool.py`,
`fixed_ops/stylize.py`, `utils/style_utils.py`, and `text_viz/`. I also searched
all operation signatures for `style=` and ran the behavior above. Evidence that
would have lowered the result would be styling available only in documentation
figures/reports, only through a caller-built diff, or only in an example; I
searched those alternatives and instead found a documented operation parameter,
implementation, and direct runtime output. No evidence supports a value below
`yes`.

### valiant — `no` (high confidence)

**Evidence.** VaLiAnT documents its visualizable construction facts as separate
plain CSV fields: `mut_position`, `ref`, `new`, and `mseq` are respectively the
mutation position, reference bases, mutated bases, and full oligo sequence
(`README.md:737-772`). Its sequence-only output is described as:

> “Comma-separated values (CSV) file containing only the label and the sequence
> of the oligonucleotides generated for any given targeton...”

Source: `README.md:805-807`. The implementation writes those values without
markup:

> `fh.write(names[0])`
>
> `fh.write(',')`
>
> `fh.write(seq)`

Source: `src/valiant/meta_table.py:676-686`. Metadata lets a caller reconstruct
a visual diff, but the global rule says that caller reconstruction is `no`.

**Search and disconfirmation attempt.** I searched all 426 tracked repository
files, all 11 pages of the separately cloned official wiki, then the full
`src/valiant` tree and `README.md` specifically, with the
families `style|styling|stylize`, `highlight`, `color|colour`, `ANSI`,
`font/background-color`, `uppercase|lowercase|swapcase`, `sequence view/display`,
and `mutation view/display`; I separately searched the output/write paths for
`SeqIO|write(|to_csv|FASTA|oligo`. The only source/docs candidates were a README
warning that dependency installers may “highlight” missing packages and calls to
`.upper()`/`.lower()` used to normalize variants, logging, SQL, and fetched
reference sequence (`common_cli.py:37`, `queries.py:135`, `sge_utils.py:32`,
`custom_variant.py:88,97`). Checked-in GenBank *input examples* contain `/Color`
feature qualifiers, but the design output path neither creates them nor styles
sequence text. Wiki color/highlight hits occur in the captions for the
“pegRNA design schematic,” tiling schematic, and cDNA design schematic
(`Saturation-prime-editing-example.md:71-73`, `Advanced-usage.md:37`,
`cDNA-example.md:17`), not in generated output. The extracted paper's only
color hit likewise describes colored rectangles in Figure 3. All are excluded
figures/schematics.

Evidence for `yes` would have been a documented output mode or parameter that
renders `mseq` with the `new` bases highlighted/cased/colored. Example-only or
single-mode support would have produced `partial`. I searched the README,
source, examples, and paper for both and found neither; the complete plain CSV
writer is positive evidence for absence, so `unknown` is not warranted.

### mpranator — `partial` (high confidence)

**Qualifying evidence.** MPRAnator has tool-provided inline sequence styling in
two web workflows. Its Transmutation handler calls the mutation routine, receives
the changed positions, and invokes the tool's own highlighter:

> `mutatedString, positionMutated = part3.mutateString(seq, numToMutate)`
>
> `outputSequenceHTMLL.append(myf.highlightString(mutatedString, positionMutated))`

Source: `iliasApp/views.py:199-207`. `highlightString` replaces each changed
character with `<span style='color:red;'>BASE</span>` by default
(`myfunctions.py:79-85`), and the result template inserts that marked-up sequence
with `{{ sequence|safe }}` (`templates/iliasApp/part3Results.html:18-20`). This is
changed-base coloring inside sequence text.

Motif Design also lowercases the background, uppercases each inserted motif,
wraps that motif in a red `<span>`, and returns both `sequence` and
`sequenceHTML` (`oligo.py:39-58,70-76`). The view passes that HTML sequence to
the result template (`iliasApp/views.py:102-123`; `templates/iliasApp/results.html:23-25`).

**Why restricted.** Styling is not a general output parameter and is not present
across the tool's workflows. Downloads deliberately use the plain sequence list:
the Transmutation handler builds `forDownload` from `finalOutputSequenceL`, not
`outputSequenceHTMLL` (`views.py:218-224`), and Motif download responses are
`text/plain` (`views.py:109-115`). More importantly, SNP code creates an
`htmlSequences` list in the combinations branch (`parseSNPs.py:35,56,95-130`)
but returns only `ExtraNamesWithArrow, ExtraSequences`
(`parseSNPs.py:199-206`). The SNP view then constructs items with only
`"header"` and plain `"sequence"` (`views.py:360-373`), so the result assembler's
fallback displays plain sequence (`part1.py:253-259`). Non-combination SNPs do
not build styled strings at all.

This is a genuine tool-provided web rendering path, not caller composition, but
it is restricted to Motif and Transmutation result pages and absent from SNP and
download outputs. Under the global definition of `partial` as a restricted
tool-provided capability, the score is `partial`.

**Search and disconfirmation attempt.** I searched all 41 tracked files for
`style|highlight|color|colour|case|ANSI|<span|<font|sequence view/display`, then
followed every `htmlSequences`, `SequencedHTML`, `sequenceHTML`,
`highlightString`, `color:red`, and `make_sequence_copies` reference through
`oligo.py`, `part1.py`, `parseSNPs.py`, `myfunctions.py`, `iliasApp/views.py`,
and both result templates. The bootstrap CSS bundle was excluded after checking
that it only styles the site generally. I also searched the extracted paper;
it does not document sequence styling.

Evidence for `yes` would have been a general documented styling mode, or the
same direct sequence marking in SNP and downloadable results; I explicitly
traced those paths and found plain output. Evidence for `no` would have been the
HTML-generating functions being dead/example code with no user-facing call;
the view-to-template calls above disconfirm that. Therefore `partial`, rather
than `yes`, `no`, or `unknown`, is the reproducible result.

### mpradesign — `no` (high confidence)

This score covers both the `mpradesigntools` package and its `designMPRA` Shiny
wrapper.

**Evidence.** The documented package result is:

> “A list of two data_frames. The first, named 'result', is a data_frame
> containing the labeled MPRA sequences.”

Source: `mpradesigntools/man/processVCF.Rd:78-86`. In the implementation,
`processSnp` builds the output `sequence` as an ordinary `paste0(...)` string
(`R/processVCFfast.R:392-418`), selects it as a plain table column, and writes the
table with `write_tsv` (`processVCFfast.R:1248-1258`). The Shiny app offers the
result only through `Download Output` and writes `vcfOut()$result` as TSV
(`designMPRA/server.R:136-165`). Its rendered tables are the failed-SNP list and
input preview, not a styled output sequence.

**Search and disconfirmation attempt.** Across all 69 package files and 21 Shiny
files I searched `R/`, `man/`, both READMEs, `ui.R`, `server.R`, scripts, and the
source R Markdown for `style|highlight|color|colour|ANSI|<span|font`,
`toupper|tolower|case`, plus `sequence|oligo|output|write|export|VCF|BED`. The
package produced no styling hit: its only `case` matches were commented function
names and prose (“three ... cases”). Shiny hits were layout
`style='text-align: center'`, a progress-notification style, and ggplot colors in
the statistical supplement. The extracted paper likewise uses colors only in
power-analysis plots. None marks changed bases in sequence text.

Evidence for `yes` would have been a package/Shiny option or returned field
containing marked-up sequence text; an example-only formatter would have yielded
at most `partial`. I searched the complete package API/manual, Shiny render and
download paths, examples/supplement, and paper and found neither. The plain TSV
construction and export settle the result as `no`.

### oligopoolcalc — `no` (high confidence)

**Evidence.** The public API documents finalization as returning a DataFrame
whose output contains only `ID`, `CompleteOligo`, and `OligoLength`, with optional
CSV output (`docs/api.md:590-616`). The implementation is equally explicit:

> “Output contains only `CompleteOligo` and `OligoLength` (sequence annotations
> are not preserved).”

Source: `oligopool/final.py:30-38`. It concatenates the sequence columns into a
plain string and writes the DataFrame as CSV (`final.py:100-116`).

**Search and disconfirmation attempt.** I searched all 76 tracked files and the
complete `oligopool/`, `docs/`, README, and examples trees for
`style|styling|stylize`, `highlight`, `color|colour`, exact ANSI escape/color
libraries, `uppercase|lowercase|swapcase`, `HTML|CSS|span|font`, and
`sequence/mutation view/display`; output schemas for every documented operation
were also inspected. The direct candidates were: README navigation CSS, example
and parser `.upper()` calls that normalize DNA, comments saying “return style”
meaning DataFrame-index convention, and an `_ANSI_RE` helper in `cli.py:91-102`
that strips terminal control codes from CLI help text. None records changed
positions or styles sequence text. The paper keyword scan found figures and
ordinary reporting language only.

Evidence for `yes` would have been an operation or parameter attaching changed-
position markup to a returned sequence or printing it with inline highlighting;
restricted/example-only support would have produced `partial`. I searched the
API docs, implementation, CLI, examples, and paper for both. There is no such
operation, and the finalizer explicitly discards annotations, so the score is
`no`, not `unknown`.

### mutationmaker — `no` (high confidence)

**Tempting but excluded evidence.** Mutation Maker does color mutation positions,
but only in outputs the row excludes. Its interactive Feature Viewer adds
mutation locations as colored rectangular features (`SSMFeatureViewer.tsx:83-103`)
and instantiates an interactive viewer with axis, zoom, and toolbar
(`shared/components/FeatureViewerComponent/index.tsx:33-40,65-80`). That is a
graph/feature view, not marked-up sequence text.

The PAS spreadsheet export also constructs Excel rich text in which codons at
`oligo.reds` and `oligo.blues` positions receive red and blue fonts
(`shared/components/SaveFile/index.tsx:226-305`), then writes an `.xlsx` workbook
and labels the button “Download as XLSX” (`SaveFile/index.tsx:326-357`). This is
a report document and is explicitly excluded by the row.

The actual on-page result tables render sequence values as plain strings. For
example, the SSM “Sequence” columns return
`record.forward_primer.sequence`/`record.reverse_primer.sequence`
(`SSMResultTable/index.tsx:88-126`), and QCLM's “Primer” column returns
`record.sequence` (`QCLMResultTable.tsx:45-65`).

**Search and disconfirmation attempt.** I searched all 300 tracked files, the
backend, frontend source, READMEs, tests, and extracted paper for
`style|highlight|color|colour|case|ANSI|HTML|CSS|span|font|richText`, then traced
all result-table, `FeatureViewer`, `withSelectedAndHighlighted`, `SaveFile`,
`reds`, and `blues` paths. The remaining color hits were validation icons,
out-of-range number styling, input normalization, and the Feature Viewer.
The paper's “interactive visualization” discussion points to the customized
Feature Viewer, not a sequence-text renderer.

Evidence for `yes` would have been the on-page `record.sequence`/primer sequence
renderers wrapping changed bases in spans or case/color markup. Example-only
support would have yielded `partial`. I searched each SSM, QCLM, and PAS result
table and all shared display components; no such renderer exists. Since the only
positive visuals are an excluded graph and excluded spreadsheet report, the
row's instruction requires `no`.

### dnachisel — `no` (high confidence)

**Tempting but excluded evidence.** DNA Chisel can emit edit annotations, but it
does not style changed bases inside sequence text. Its README distinguishes the
plain final string from an annotated record:

> `final_sequence = problem.sequence  # string`
>
> `final_record = problem.to_record(with_sequence_edits=True)`

Source: `README.rst:79-82`. The `with_sequence_edits` option is documented as
adding “annotations representing each nucleotide change”
(`RecordRepresentationMixin.py:128-130`). The implementation converts difference
segments to red Biopython/GenBank **features** labeled `before=>after`
(`RecordRepresentationMixin.py:186-206`) and either returns that record or writes
a GenBank file (`:189-192`). The CLI likewise says this option annotates edits in
the GenBank file (`dnachisel/cli.py:7-14,63-72`). These are feature annotations,
not inline sequence-text styling.

DNA Chisel's other visual path is explicitly report/plot based: the report
implementation says it writes “a PDF of plots of the sequence ... and an
annotated genbank” (`dnachisel/reports/optimization_reports.py:82-92`) and uses
Matplotlib/GraphicRecord plotting (`:117-140,170-184`). The paper likewise calls
it a “PDF report ... indicating modified nucleotides.” Reports and plots are
excluded from this row.

**Search and disconfirmation attempt.** I searched all 281 tracked files and the
full `dnachisel/`, `docs/`, README, examples, and paper for
`style|highlight|color|colour|upper/lower/swapcase|ANSI|HTML|CSS|span|font`,
`with_sequence_edits`, and `sequence_edits`. I inspected
`biotools/sequences_differences.py`, `RecordRepresentationMixin.py`, `cli.py`,
the report implementation, GenBank API docs, and README. Sequence-difference
utilities return a Boolean array, a count, or coordinate segments
(`sequences_differences.py:6-42`), never marked text. Colors in checked-in
GenBank examples are feature colors.

Evidence for `yes` would have been a renderer/parameter producing highlighted,
cased, or colored bases in sequence text itself; restricted or example-only
inline rendering would have produced `partial`. I searched both the direct
string path and every edit/report/record path. Only excluded annotations,
figures, and reports exist, so `no` is determinate.

### tangermeme — `no` (high confidence)

**Tempting but excluded evidence.** tangermeme has sophisticated per-position
coloring, but its API names and types make clear that it is a plot:

> `def plot_logo(... color=... ) -> matplotlib.axes.Axes:`
>
> “Make a logo plot and optionally annotate it.”

Source: `tangermeme/plot.py:142-167`. The `color` parameter may color each
position independently (`plot.py:186-206`), and the implementation draws glyph
paths with Matplotlib face colors and adds them to an Axes
(`plot.py:330-380`). The current release notes also describe
`interactive_logo` as an interactive counterpart to `plot_logo`, with colored
boxes behind logo glyphs (`docs/whats_new.rst:71-72`). These are precisely the
plots excluded by the row, not rendered variant sequence text.

**Search and disconfirmation attempt.** I searched all 141 tracked files and the
full package, docs, README, skills references, tests, and extracted paper for
`style|highlight|color|colour|ANSI|HTML|CSS|span|font|upper/lower/swapcase`,
`sequence/mutation view/display`, and all `read|write|plot` definitions. I then
repeated the styling search excluding `plot.py`, plot tests, and generated
notebooks to look specifically for a non-plot renderer. Remaining candidates
were conceptual statements that attribution values “highlight” important input
characters, `.upper()` calls for genome input normalization, documentation
headings, and plot documentation. The public I/O module has readers but no
styled-sequence writer; `one_hot_encode` returns a tensor
(`tangermeme/utils.py:427-465`).

Evidence for `yes` would have been a design/result operation returning or
printing variant sequence text with changed positions marked. Example-only or
restricted non-plot rendering would have produced `partial`. I searched all
non-plot source and documentation paths and found none. Because the only
color-capable result is a Matplotlib/mpld3 figure, the mandated value is `no`.

## Row-level finding

The row discriminates on one consistent threshold: one tool has full direct
sequence-text styling (`poolparty`), one has a restricted web-page implementation
(`mpranator`), and six have none under the stated exclusions. The exclusions are
substantive rather than cosmetic: `mutationmaker`, `dnachisel`, and `tangermeme`
all visually encode positions, but only in spreadsheets/reports, annotated
records/plots, or sequence-logo/feature views. Counting those would collapse the
row into a general “has visualization” measure and contradict its definition.

## Complete search log

The following records every search performed. Line-range reads used to inspect
hits are represented by the source citations above; they did not discover
additional files.

1. Read the complete definitions preamble and row with
   `sed -n '1,240p' .../ROW_DEFINITIONS.md`.
2. Enumerated the survey tree and possible local repository directories with
   `find ... -type f` / `find ... -type d`. Searched only URL/source headings in
   the eight `final/<slug>.md` lead files using
   `rg -n 'https?://|github\.com|Source|Repository|Documentation|PyPI'`; this was
   source discovery only, and no lead conclusion was used.
3. PoolParty first-pass search over all source/docs:
   `rg -n -i 'style|highlight|color|colour|case|uppercase|lowercase|ansi|html|changed|mutation'`.
   Follow-up definition search:
   `rg -n 'def (stylize|mutagenize|print_library)|style='`. Inspected the
   styling/mutagenize docs and relevant implementations and ran five behavior
   invocations. The first omitted `pp.init()` and did not finish within the
   10-second yield; the second incorrectly treated `p[0]` as a sequence and
   raised `AttributeError`; the third exposed a pandas formatting incompatibility
   when `DataFrame.to_string()` tried to index `SeqStyle`; and the final two
   inspected the sequence/style cells directly and established the ANSI output
   reported above.
4. Attempted shallow clones in the sandbox; DNS was unavailable. Repeated the
   same official `git clone --depth 1` operations with approved network access
   for VaLiAnT (`develop`), MPRAnator, mpradesigntools, designMPRA, oligopool,
   Mutation_Maker, DnaChisel, and tangermeme. Used `git rev-parse`, `git log -1`,
   and `git ls-files | wc -l` to record commits and search scope.
5. Extracted all eight supplied PDFs with PyMuPDF exactly as directed, then
   searched each text for
   `style|highlight|color|colour|uppercase|lowercase|case|sequence view/display|displayed sequence|visualize/visualise|figure|report`.
   For MPRAnator I also searched `highlight|red|color|case|display|result page`.
6. Cross-repository first pass over all seven cloned repositories used
   `style/styling/stylize`, `highlight`, `color|colour`, `uppercase|lowercase`,
   `ANSI`, `font/background-color`, `text-decoration`,
   `sequence/mutation view/display`, and `diff`, excluding `.git`, locks, large
   data formats, SVGs, and notebooks. A companion `rg -l` pass listed every
   candidate file. Data/example formats were then searched or inspected
   separately where relevant (notably VaLiAnT GenBank `/Color` fields and
   notebook/report artifacts).
7. VaLiAnT targeted searches: README `output|FASTA|metadata|sequence`; source
   `SeqIO|write(|to_csv|FASTA|fasta|oligo`; and source/docs
   `style|highlight|color|colour|upper(|lower(|swapcase|ANSI|HTML|CSS`.
   A final exact candidate pass used
   `style|highlight|color|colour|ANSI|<span|font|.upper(|.lower(|swapcase`.
   The official wiki was cloned separately, all 11 Markdown pages were listed,
   and the complete wiki was searched for
   `style|highlight|color|colour|uppercase|lowercase|ANSI|font|background-color|sequence display|mutation display|output|FASTA|metadata`.
8. MPRAnator targeted searches: templates enumerated with `find`; all Python/app
   files searched for `FASTA|result|oligo|sequence|download|render|HttpResponse`;
   non-static source searched for
   `style|highlight|color|colour|upper(|lower(|swapcase|ANSI|<span|<font`;
   follow-up reference searches used
   `htmlSequences|SequencedHTML|color:red|make_sequence_copies(` and
   `def highlightString|def createMPRAResultOutput|backgroundHTML|sequenceHTML`.
9. mpradesigntools targeted searches: `R/`, `man/`, and README for
   `style|highlight|color|colour|upper/lower|toupper|tolower|ANSI|HTML|CSS|span|font`,
   then manual/docs for `sequence|oligo|output|write|export|VCF|BED`, and
   `processVCFfast.R` for
   `sequence=|sequences|outPath|write.table|write_tsv|return(`. designMPRA was
   searched separately for the same styling terms plus Shiny
   `render*`/`downloadHandler`, excluding only the generated Supplement HTML
   after checking its candidates came from embedded dependencies/images.
10. Oligopool Calculator targeted searches: complete package/docs/README/examples
    for styling, highlighting, colors, case, ANSI, HTML/CSS/span/font, and output
    terms; then a narrower `highlight|color|colour|ANSI|<span|font|swapcase|toupper|tolower`
    pass. I attempted a literal escape-code/library regex that `rg` rejected as
    an unsupported backreference because the shell-transmitted `\033` was
    malformed; I did not treat that failed search as evidence. The successful
    searches nevertheless found and inspected the only ANSI candidate,
    `_ANSI_RE`/`_strip_ansi` in `cli.py:91-102`, as well as every public output
    schema and `final.py`.
11. Mutation Maker targeted searches: all result-table and Feature Viewer files
    for `sequence|oligo|primer|highlight|style|color|FeatureViewer`; all
    frontend/backend/docs for
    `style=.color|color:|<span|dangerouslySetInnerHTML|innerHTML|mark>|backgroundColor|textDecoration|toUpperCase|toLowerCase`;
    and follow-ups for `richText`, `reds`, `blues`, `SaveFile`, and each SSM,
    QCLM, PAS result component.
12. DNA Chisel targeted searches: source/docs/README/examples for
    `style|highlight|color|colour|upper(|lower(|swapcase|ANSI|HTML|CSS|span|font|with_sequence_edits|sequence_edits`;
    followed by `sequence_edits|feature|color|label` in the record mixin and CLI.
    Inspected the complete sequence-difference utility and relevant report paths.
13. tangermeme targeted searches: full package/docs/README for styling/case terms;
    `plot_logo` and all color arguments/implementation; all public
    `read|write|plot` definitions; then the same styling search excluding
    `plot.py`, plot tests, and notebooks to test for any non-plot renderer. The
    latter found only attribution prose, normalization, and plot documentation.
14. Finally, I reran compact exact-candidate passes for VaLiAnT,
    mpradesigntools/designMPRA, Oligopool Calculator, and non-plot tangermeme to
    ensure broad-regex noise had not hidden a qualifying match. Candidate
    dispositions are reported in their per-tool sections above.
