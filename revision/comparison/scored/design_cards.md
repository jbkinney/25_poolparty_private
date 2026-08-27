# Audit: per-sequence design cards (`design_cards`)

Date of audit: 2026-08-15.

## Measurement instrument and operational test

I read the preamble and the complete **Global rule — binding on every row below**
in `25_poolparty_private/revision/tables/ROW_DEFINITIONS.md`, then row 13. I used
the following operational test, stated before inspecting any tool:

> For each tool, check whether its documented design operation emits every
> output sequence with tool-attached, structured, machine-readable metadata
> that records construction provenance beyond merely the variant/mutation
> identity; score `yes` if all elements hold, `partial` if provenance exists but
> is non-per-sequence, unstructured, or identity-only, and `no` only after
> checking the design outputs, schemas/models, serializers/exporters, examples,
> docs, and paper for such a tool-provided field or object.

I treated a native DataFrame/CSV field, JSON field, or sequence-feature
qualifier as structured. A descriptive string does not become a structured
record merely because its producer uses delimiters; a documented delimited
header is therefore evidence of provenance but is `partial` when the header is
the only carrier. The test is format-neutral: an in-memory value and a file both
qualify. Optional metadata qualifies because `cards`, export, or annotation is a
tool-provided parameter/mode, not user reconstruction. Conversely, values kept
only in local variables and discarded before return do not qualify.

The seven remote repositories were shallow-cloned solely for read-only search
and pinned at the commits below. PoolParty was inspected in the supplied local,
read-only repository and was also run in the supplied venv with
`PYTHONDONTWRITEBYTECODE=1`. The files under `tool_survey/final/` were used only
as navigation leads and are not evidence here.

Common search log: I first read `ROW_DEFINITIONS.md` with `sed`, listed the
available survey/paper files and possible local source directories with `find`,
and searched the eight `tool_survey/final/<slug>.md` lead records only for
repository URLs and likely primary-source filenames (terms including
`Repository|GitHub|docs|Source|URL|file|metadata|output|schema|model|class|CSV|FASTA`).
No proposition or score below comes from those lead records. I then cloned the
seven cited repositories, recorded every `git rev-parse HEAD`, searched each
complete tree as documented in its section, and extracted all eight author PDFs
with PyMuPDF. The first PDF pass counted the same common terms in every paper:
`metadata`, `provenance`, `header`, `output`, `CSV`, `FASTA`, `data frame`,
`table`, `annotation`, `label`, `barcode`, `primer`, and `oligo`; the targeted
follow-up paper searches are listed under each tool. File listings (`find`) and
line-numbered reads (`nl`/`sed`) were used to make sure a grep hit was understood
in its output/schema context rather than scored from the hit alone.

## Result summary

| tool | value | confidence |
|---|---|---|
| PoolParty | **yes** | high |
| VaLiAnT | **yes** | high |
| MPRAnator | **partial** | high |
| MPRA Design Tools | **yes** | high |
| Oligopool Calculator | **yes** | medium |
| Mutation Maker | **yes** | high |
| DNA Chisel | **yes** | high |
| tangermeme | **no** | high |

## PoolParty — **yes** (confidence: high)

### Evidence

Primary source: supplied local repository
`poolparty-statecounter/poolparty/`.

The documentation states the capability directly:

> “PoolParty can automatically pair each generated sequence with a **design
> card**, a DataFrame row that records how the sequence was constructed.”

It immediately specifies construction fields beyond a mutation label:

> “Columns report the changes applied by each operation: mutation positions,
> substituted characters, scores, orientations, and more.”

Source: `poolparty-statecounter/poolparty/docs/metadata/design_cards.rst:4-9`.
The same page says the mode is tool-provided and opt-in (`cards`), not a user
join: “Design cards are opt-in: unless you pass the `cards` parameter, the
output contains only `name` and `seq`” (`:11-12`). It documents list and dict
forms, automatic operation-prefixed column names, and custom column names
(`:41-69`).

The documented schema includes both per-operation state and intermediate
sequence (`:73-95`) and construction-specific fields, including
`mutagenize.positions/wt_chars/mut_chars`, codon and amino-acid changes,
`stack.active_parent`, `repeat.repeat_index`, orientation, recombination
breakpoints and source-pool assignments, shuffle permutation, and source
sequence index (`:112-168`). In particular, recombination provenance is
described as:

> “`breakpoints`, `pool_assignments` — Breakpoint positions and which source
> pool contributed each segment.”

Source: `docs/metadata/design_cards.rst:147-150`.

The return contract says:

> “DataFrame with columns: name, seq, plus any requested design card columns.”

Source: `src/poolparty/generate_library.py:48-55`. At row generation time the
implementation topologically visits every operation, calls
`output_seq, raw_card = op.compute(...)`, filters the requested card, and writes
each key/value into the same row as the final `seq` (`generate_library.py:253-258,
270-319,321-330`). The operation API itself returns `tuple[Seq, dict]`, “Output
Seq (with string and style) and design card dict”
(`src/poolparty/operation.py:286-301`), with native dict filtering and column
mapping at `operation.py:338-386`.

Behavioral check, run with the supplied read-only Python 3.12 venv:

```text
name      seq  construction_state mut_pos     wt    mut
None ATCGATAA                   0  (6, 7) (C, G) (A, A)
None TTCGAACG                   1  (0, 5) (A, T) (T, A)
```

The code was a documented `from_seq(...).mutagenize(...,
cards={...}).generate_library(seed=7)` call. This establishes actual per-row,
structured construction fields rather than relying only on prose.

### Searches and disconfirmation attempt

Searched the full `src/poolparty` tree and `docs/metadata` with case-insensitive
terms `design.?card|provenance|metadata|construction|history|lineage|origin|created|operation|variant`.
Read `docs/metadata/design_cards.rst`, `src/poolparty/generate_library.py`, and
`src/poolparty/operation.py`, and searched the author PDF `main.pdf` for
`metadata`, `provenance`, `header`, `output`, `CSV`, `FASTA`, `data frame`,
`table`, `annotation`, `label`, `barcode`, `primer`, and `oligo` (three
`metadata` hits; the repository was treated as authoritative as instructed).

A different value would have resulted if cards were only names, an example-only
user join, an unstructured report, or global rather than per-row metadata. I
specifically searched the implementation for where raw cards are returned,
filtered, and attached to the final row, and ran the API. Those disconfirming
conditions are absent. `yes` follows even though cards are opt-in: `cards` is a
documented tool parameter, exactly what the global rule requires.

## VaLiAnT — **yes** (confidence: high)

### Evidence

Primary repository: <https://github.com/cancerit/VaLiAnT>, branch `develop`,
commit `8796cc112dafd4919fec59913f58cd2be87c45eb`.

The README defines an explicit per-oligo structured output:

> “Comma-separated values (CSV) file containing name, label, and all metadata
> of the oligonucleotides generated for any given targeton.”

Source: `README.md:737-740`. The schema is not restricted to variant identity.
Each row includes reference interval/strand/reverse-complement status and source
sequences (`:749-761`), variant fields and the `mutator` that generated it
(`:762-770`), the full sequence with and without adaptors (`:771-773`), applied
PAM-protection annotations and sgRNA IDs (`:774-775`), and background variants
and the background-altered reference sequence (`:776-780`). The literal CSV
header and example row appear at `README.md:782-787`.

The code fixes those fields in a native schema. `META_CSV_FIELDS` contains 32
named columns including `revc`, `ref_seq`, `pam_seq`, `mutator`, `mseq`,
`mseq_no_adapt`, `pam_mut_annot`, `pam_mut_sgrna_id`, `background_variants`, and
`background_seq` (`src/valiant/meta_table.py:49-82`). `write_meta_record` accepts
these as distinct typed arguments (`meta_table.py:175-208`) and writes one
record. The intermediate `MetaRow` is a slotted dataclass with separate variant,
mutator, oligo, mutation-type, exon/codon, PAM-edit, and sgRNA fields
(`src/valiant/meta_row.py:54-76`).

The authors’ paper corroborates the intended machine-readable output. Its
Supplementary Table 3 calls `*_meta.csv` a “Metadata and sequence annotation
file” for bioinformatic analyses and distinguishes it from `*_unique.csv`, the
synthesis-only name/sequence file (Barbon et al. 2022 PDF, supplementary p. 11).
The paper also says a VCF alias and optional identifier are assigned “to preserve
variant provenance” (main article p. 5). These are additional confirmation; the
repository schema is the decisive evidence.

### Searches and disconfirmation attempt

Searched the complete repository with
`metadata|meta.?table|header|provenance|construction|history|design|oligo.?id|name|sequence.?id|variant|mutation|mutator|to_csv|output`.
Read the README sections “Output files,” “Oligonucleotide metadata file,”
“Variant file,” and “Unique oligonucleotides file”; read
`src/valiant/meta_row.py`, `meta_table.py`, `cdna_proc.py`, and `sge_proc.py`;
and searched all repository docs/examples/tests. The paper was searched for the
common audit term set above and specifically `provenance`, `metadata file`,
`CSV file`, and `oligonucleotide metadata`.

A different value would have resulted if the only metadata were mutation/HGVS
identity, if metadata were targeton-level only, or if the synthesis sequence
could not be linked to the rich row. I searched for exactly those failure modes:
the metadata table contains `mseq` itself and `oligo_name`, plus construction
fields such as adaptors, PAM protection, reverse complementation and background
edits. This is per output oligo and beyond mutation identity, so `yes`.

## MPRAnator — **partial** (confidence: high)

### Evidence

Primary repository: <https://github.com/hemberg-lab/MPRAnator>, commit
`9969790d62410138d4281b7955da6d085f07b1bc`.

There is real per-sequence provenance, so `no` would understate the tool. The
shipped documentation says:

> “The description line (header) has information about the options chosen by
> the user during submission. A header is composed of one or more DESCRIPTORs
> and each DESCRIPTOR is composed of a LABEL and INFO. The descriptors are
> delimited by a |, i.e. a ‘pipe’.”

Source: `iliasApp/templates/iliasApp/docs.html:45-51`; the MPRA-motif page repeats
the contract and says it applies “for each sequence” (`:84-111`). The worked
header records two motif placements and two restriction elements, not merely a
mutation name (`:125-136`). The Transmutation header records mutation count and
whether scrambling occurred (`:150-160`).

The implementation appends `BARCODE`, `RESTRICTION`, and `ADAPTER` descriptors
to each sequence header and flags duplicate restriction sites
(`part1.py:210-292`). The Transmutation view emits, per sequence,
`Mutated_nucleotides`, `Scrambled`, `Reversed`, and `Complemented`
(`iliasApp/views.py:178-221`). SNP/indel generation also records length
normalisation, e.g. “removed ... nucleotides from the left edge” or “added ...
nucleotides from the right edge” (`parseSNPs.py:395-428`). These are genuine
construction facts beyond identity.

The restriction is format/structure. The only carrier is a FASTA description
string assembled by concatenation. No native per-sequence dict/table/JSON field
is returned. Syntax is also not uniform: the documented general form is
`<LABEL> - <INFO> |` (`views.py:232-252`), while the SNP header is assembled as
`... | SNP ... | Nucleotide change ...` and edge-fix phrases
(`parseSNPs.py:395-428`). Thus the header is semi-structured and human-readable,
but downstream code must parse operation-specific prose. Under the stated test,
this is the row definition’s “some provenance recorded, but ... not structured”
case: `partial`.

### Searches and disconfirmation attempt

Searched every repository file with
`header|metadata|descriptor|label|info|provenance|construct|history|fasta|output|oligo|barcode|restriction|mutated|scrambled|reversed|complemented|prob`.
Read `part1.py`, `oligo.py`, `parseSNPs.py`, `iliasApp/views.py`, the entire
shipped `docs.html`, result templates, downloadable clients, models, and tests.
The two-page author paper was searched with the common term set; it has no
`metadata`, `provenance`, `header`, `CSV`, `FASTA`, `data frame`, `annotation`,
or `label` hit, so it supplies no richer output contract than the repository.

Evidence that would have produced `yes` was a native per-sequence record or
sidecar (CSV/TSV/JSON/database field) carrying these descriptors, or a single
uniform typed header grammar. I explicitly searched for `dict`, models,
serializers, CSV/TSV/JSON, and all output templates and found none. Evidence that
would have produced `no` was headers carrying only input names or mutation
identity; the construction-option and length-repair descriptors disconfirm that.

## MPRA Design Tools — **yes** (confidence: high)

### Evidence

Primary repositories:

- <https://github.com/andrewGhazi/mpradesigntools>, commit
  `afd386ef12051bb0a37ad63a6f92acd555246584` (decisive package source);
- <https://github.com/andrewGhazi/designMPRA>, commit
  `0cf56ee602fc86dde705906d071a39cbdf6e99a8` (companion web app).

The documented design API return for `processSnp` is:

> “a data_frame of labeled sequences with appropriate information on the
> changes made”

Source: `mpradesigntools/man/processSnp.Rd:44-50`, generated from
`R/processVCFfast.R:187-205`. `processVCF` returns a named list whose `result`
DataFrame contains labeled MPRA sequences and whose `failed` DataFrame records
failures and reasons (`man/processVCF.Rd:78-87`; source
`R/processVCFfast.R:1068-1071`).

The implementation attaches native fields to each successful sequence:

```r
select(ID, type, allele, snpIndex, totIndex, barcode, sequence,
       site_fix_info, notes)
```

Source: `R/processVCFfast.R:1243-1254` and, for mixed success/failure,
`:1281-1292`. This is more than variant identity. `notes` records, for example,
that input alleles were flipped due to the VCF `RV` tag (`:236-287`) and exactly
how much context was shortened on each side to meet `max_construct_size`
(`:317-324`). When an aberrant restriction site is repaired, distinct fields
`aberrant_pattern`, `fixed_pattern`, `fixed_pattern_index`, and
`constrseq_fixed` are generated (`:146-172`) and nested as per-row
`site_fix_info` (`:452-481`). The sidecar can also be unnested and written as TSV
(`:1299-1329`). These are explicit construction-repair records, not names.

The paper describes the public output as “a tab-separated file containing the
MPRA sequences for their experiment” and the R package as providing more
customizable sequence generation (Ghazi et al. 2018, main article p. 2). The
repository’s current return schema supplies the stronger provenance evidence.

### Searches and disconfirmation attempt

Searched both complete repositories with
`metadata|provenance|construct|history|header|name|id|variant|mutation|barcode|adapter|sequence|output|data.frame|dataframe|table|csv|fasta|bed`.
Read both READMEs, all R package functions and `.Rd` pages, the Shiny server,
scripts, and Supplement Rmd. In the author PDF I searched the common term set
and specifically `data frame`, `labeled`, `output file`, `barcode sequence`,
`oligo sequence`, and `annotat`.

A different value would have resulted if the returned table stopped at
`ID/type/allele/sequence`, if repair information were only printed globally, or
if rich fields existed only in a repo example. I checked the exported package
function and generated manual: `site_fix_info` and `notes` survive into the
official `result` DataFrame and TSV. The fields are per sequence and structured,
so `yes`.

## Oligopool Calculator — **yes** (confidence: medium)

### Evidence

Primary repository: <https://github.com/ayaanhossain/oligopool>, commit
`b88fa394ca67ed4c48ec41127b5707694ee7cf0a` (tag visible in the shallow clone:
`v2026.02.22.1`).

The current README describes the official data contract:

> “**DataFrame-centric:** modules operate on CSV/DataFrames and return updated
> tables plus `stats`.”

Source: `README.md:37-48`. Its API synopsis says “Modules return (DataFrame,
stats). Chain them iteratively” (`README.md:151-157`). The design tables are not
just sequence/name pairs: they retain one row per `ID` and one named column per
constructed molecular element. The documented workflow inserts a designed
forward primer as `FwdPrimer`, a barcode as `BC1`, and a Tm-matched reverse
primer as `RevPrimer`, all into the same per-variant DataFrame
(`docs/docs.md:1304-1358`), then explicitly saves it as
`library_design.csv` (`:1367-1374`).

Column placement itself records construction order. The barcode documentation
shows:

```text
Input:  Left, Core, Right
barcode(..., barcode_column='BC', left_context_column='Left',
        right_context_column='Core')
Output: Left, BC, Core, Right
```

Source: `docs/docs.md:317-345`. It explains that when both context columns are
given, “newly designed columns are inserted between those context columns in the
DataFrame” and that this resolves final linear column ordering
(`docs/docs.md:255-275`). Thus the per-ID row records the actual primer, motif,
barcode, spacer, core sequence, and their construction order—well beyond a
variant label.

The API docstring calls `input_data` a CSV/DataFrame “with annotated oligo pool
variants,” requires a unique `ID`, takes the tool-created `barcode_column`, and
returns a pandas DataFrame (`oligopool/barcode.py:35-68`). Equivalent contracts
exist in `primer.py`, `motif.py`, and `spacer.py`.

Important caveat checked rather than hidden: `final()` concatenates the row into
`CompleteOligo` and `OligoLength`, then drops the component annotation columns
(`docs/docs.md:739-765`; `oligopool/final.py:15-38`). This does not reduce the
score under the operational test because the tool’s documented design output is
the annotated, per-ID DataFrame and the docs explicitly instruct saving it; the
synthesis-only projection is an additional output mode, not the only carrier.
No user reconstruction is needed to obtain the annotated table.

### Searches and disconfirmation attempt

Searched the full repository (package, docs, examples) with
`metadata|provenance|construct|history|header|name|id|variant|mutation|barcode|adapter|sequence|output|dataframe|data.frame|table|csv|fasta|dictionary|dict|column|context`.
Read README, `docs/api.md`, the relevant User Guide sections, all design-module
docstrings and returns, `final.py`, `merge.py`, `join.py`, `verify.py`, and the
design assembly example. The author PDF was searched with the common term set
and specifically `design pars`, `input table`, `annotated`, `output`,
`oligo pool variants`, and `column`.

A different value would have resulted if the component columns were global
statistics, if operations returned only their newly generated vector detached
from `ID`, if the annotated table existed only in the example parser, or if
`final` were the sole design output. I searched the documented API and source
for all four core element operations: each returns an updated per-ID DataFrame.
The final projection creates some interpretive risk about “attached,” hence
medium rather than high confidence, but the documented annotated output meets
every element of the test and supports `yes`.

## Mutation Maker — **yes** (confidence: high)

### Evidence

Primary repository: <https://github.com/ra100/Mutation_Maker>, commit
`396c1c0ede7529f3dedf65e821e8c1f20c9a7043`. (The paper’s Merck URL is no
longer available; this is the supplied successor source.)

Mutation Maker defines native JSON schemas containing each output sequence and
its design properties. For SSSM, `PrimerOutput` includes direction, sequence,
start, three-prime temperature, GC content, and whether parameters are in range;
it is attached as forward and reverse primers to a mutation result that also
records non-optimality and overlap design (`backend/mutation_maker/ssm_types.py:
221-245`).

For QCLM/MSDM, each `PrimerOutput` carries:

> `sequence`, `start`, `length`, `temperature`, `gc_content`,
> `degenerate_codons`, and `overlap_with_following`

Source: `backend/mutation_maker/qclm_types.py:110-117`. These structured primer
objects are grouped with the mutations they cover (`:119-126`) inside a
`QCLMOutput` JSON object that also records input data, full sequence, GOI offset,
and results (`:129-133`). Length, thermodynamics, start, degeneracy and overlap
are construction facts beyond mutation identity.

For PAS gene synthesis, every `PASOligoOutput` has sequence, mixture ratio,
mutation annotations and positional highlight lists (`pas_types.py:168-183`). It
is nested in a fragment result carrying fragment/start/end/length and overlap
sequence, Tm, GC, and length (`pas_types.py:186-196`), and the complete output
also retains the typed input (`:213-216`). The generator actually creates a
separate `PASOligo(sequence=dna, ratio=concentration)` for each mutation
combination and normalizes ratios (`generate_oligos.py:122-142`).

The authors’ paper independently describes “interactive table ... and Excel
export,” says results plus input parameters can be exported, and states:

> “Designed oligos can be exported as annotated XLSX file.”

Source: Hiraga et al. 2021 PDF, p. 9. The schema above shows what that annotation
contains.

### Searches and disconfirmation attempt

Searched the entire repository with
`metadata|provenance|construct|history|header|name|id|variant|mutation|output|result|sequence|primer|oligo|json|dict|model|schema`.
Read all `ssm_types.py`, `qclm_types.py`, `pas_types.py`, `pas_output.py`,
solution/generation files, API serializers and XLSX export paths, frontend result
tables, tests, and README. Searched the author PDF with the common term set and
specifically `output file`, `output`, `download`, `export`, `mutation`,
`oligo sequence`, and `primer sequence`.

A different value would have resulted if the JSON output paired sequences only
with mutation labels, or if thermodynamic/fragment/overlap information were
global rather than attached to the primer/oligo result. The three workflow
schemas and export code were searched explicitly; each embeds rich typed fields
at the relevant sequence object or its containing fragment result. Therefore
`yes`.

## DNA Chisel — **yes** (confidence: high)

### Evidence

Primary repository:
<https://github.com/Edinburgh-Genome-Foundry/DnaChisel>, commit
`68c09304341c3656f3dfe63eda37757d6a7b3917`.

The documented API returns the optimized sequence as an annotated Biopython
record:

```python
final_sequence = problem.sequence
final_record = problem.to_record(with_sequence_edits=True)
```

Source: `README.rst:79-82`. `to_record` is specifically documented to
“Return/write record representing the final sequence and problem”; its options
add annotations for constraints, objectives, and every nucleotide edit
(`dnachisel/DnaOptimizationProblem/mixins/RecordRepresentationMixin.py:76-103,
118-141`). The default is `with_constraints=True, with_objectives=True`, while
edits are an explicit tool parameter (`:76-88`).

The implementation creates one `SeqRecord` for the final sequence and attaches
constraint and objective features, original features, and optionally edit
features (`RecordRepresentationMixin.py:153-192`). These are structured
Biopython `SeqFeature` objects. Each specification feature has a genomic/sequence
location and a `qualifiers` dictionary containing role, label, and color
(`dnachisel/Specification/FeatureRepresentationMixin.py:13-54`). Each edit
feature records its location, `before=>after` label, and `is_edit="true"`
(`RecordRepresentationMixin.py:194-207`). Thus the record holds both what
changed and why/how the sequence was constructed (the applied constraints and
objectives), not just edit identity.

The report mode also writes structured before/after constraint and objective
CSV tables and `final_sequence_with_edits.gb`
(`dnachisel/reports/optimization_reports.py:366-380,404-420`). The author paper
calls the reports “comprehensive optimization reports for traceability and
troubleshooting” and says the optimized GenBank sequence contains annotations
indicating modified nucleotides (Zulkower & Rosser 2020 PDF, p. 1). Figure 1C is
described as summarizing specifications, constraint resolutions, and objective
score improvements (p. 2).

### Searches and disconfirmation attempt

Searched README, all docs, examples, and the full package with
`metadata|provenance|construct|history|sequence_edits|with_sequence_edits|annotate_differences|SeqFeature|qualifiers|constraints|objectives|to_record|write_record|final_sequence|report`.
Read `RecordRepresentationMixin.py`, `FeatureRepresentationMixin.py`, GenBank
helpers, report writers, CLI, GenBank API docs and examples. Searched the author
paper with the common term set and specifically `annotation`, `annotated`,
`GenBank`, and `report`.

A different value would have resulted if edits were the only annotations
(identity-only), if constraint/objective reports were detached global prose, or
if the annotated record were an example-only helper. I checked the public
`to_record` defaults and implementation: typed specification features are
attached directly to the final sequence record, and the README documents the
call. Although DNA Chisel normally designs one sequence per problem, “per
sequence” still holds for every sequence it emits. Score: `yes`.

## tangermeme — **no** (confidence: high)

### Evidence

Primary repository: <https://github.com/jmschrei/tangermeme>, commit
`2006b310cd72a28c56c3ea4ba67f738fff74bb89`.

The decisive documented contract is:

> “all data are PyTorch tensors”

Source: `README.md:80-82`. The design example binds only `X_hat` from
`greedy_substitution` (`README.md:268-280`). The public function return is a
bare tensor:

> “`X: torch.Tensor ... The edited sequence.`”

Source: `tangermeme/design/greedy_substitution.py:146-150`; implementation
returns only `X` (`:233-249`). `beam_substitution` likewise documents and
returns only a tensor of designed sequences (`design/beam_substitution.py:
173-176,279-314`), and `screen` returns only its tensor of best sequences
(`design/screen.py:140-143,183-185`).

The strongest attempted disconfirmation is in `greedy_marginalize`: it does
internally accumulate `chosen_motifs` and `chosen_pos`, exactly the raw material
for construction provenance (`design/greedy_marginalize.py:214-229`), but it
then builds and returns only a one-hot-encoded sequence (`:242-246`). Internal
bookkeeping discarded before return is excluded by the global rule.

Complete package class search (`rg '^class ' --glob '*.py' tangermeme`) found
exactly six classes: five output `NamedTuple`s and one warning. The result types
hold predictions/attributions/references, not designed sequences with identity
or provenance: `PerturbationResult`, `PerturbationAnnotationsResult`, and
`AttributionReferencesResult` contain only tensor fields
(`tangermeme/results.py:16-52`); `SpaceResult` and
`SaturationMutagenesisRawResult` similarly hold model outputs. No sequence,
library, record, metadata, design-result, or edit-history class exists.

### Searches and disconfirmation attempt

Searched every file in `tangermeme/`, all docs/tutorials/vignettes, README, and
tests with
`metadata|provenance|construct|history|header|name|id|variant|mutation|output|result|sequence|annotation|namedtuple|dataclass|dict|dataframe|record|label`.
Read every file under `tangermeme/design/`, `results.py`, `space.py`,
`saturation_mutagenesis.py`, `ersatz.py`, design docs/tutorial, and ran the
complete class grep above. Searched the author PDF with the common term set:
zero hits for `metadata`, `provenance`, `header`, `CSV`, `FASTA`, or `data
frame`; its annotation discussion concerns genomic spans used for perturbation,
not metadata attached to designed sequences.

Evidence that would have produced `yes` or `partial` was any public design return
type/field, DataFrame, record, FASTA header, serializer, or documented sidecar
carrying chosen motif/position, objective score, seed, iteration, or source ID
with each returned sequence. I searched all those terms, all four design return
paths, all package classes, docs, examples, tests and the paper. None exists.
Because the internal motif/position lists are discarded and users would have to
reconstruct changes by comparing tensors and adding their own bookkeeping, the
global rule requires `no`, not `partial`.

## Row-level discrimination finding

The row is applicable on one consistent scale and does discriminate, but the
distribution is top-heavy: six `yes`, one `partial`, one `no`. It separates
native per-sequence construction records (DataFrame/CSV, JSON, or sequence
features) from MPRAnator’s semi-structured header and tangermeme’s sequence-only
tensor returns. Its weakest edge is Oligopool Calculator: the annotated design
DataFrame clearly qualifies, but the separate synthesis projection deliberately
drops annotations. That caveat does not require a replacement rule and did not
change the score under the operation-based global rule.
