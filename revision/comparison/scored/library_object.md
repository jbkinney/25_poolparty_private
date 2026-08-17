# `library_object` audit

Date of audit: 2026-08-15.  The eight scores were derived independently from
primary sources at the revisions listed below.  The files under
`tool_survey/final/` were used only as leads for locating repositories and likely
source files; no statement or score in those records is evidence here.  I did not
open or use `library_object.SUPERSEDED.md`.

## Binding instrument and operational test

I read both the preamble (including **Global rule — binding on every row below**)
and `### 1. Library object` in `revision/tables/ROW_DEFINITIONS.md` before looking
at any tool.

**Operational test (stated before tool inspection):** for each tool, identify a
documented design API call and check whether one call returns the entire
multi-sequence result as (yes) a tool-defined library type carrying
identity/provenance/membership behavior, (partial) a durable documented in-memory
general-purpose container, or (no) only files, single-sequence/problem objects, or
no whole-library return; then search its API, examples, source return paths, and
paper for the missing stronger return form.

This applies the global rule literally: an internal purpose-built object that is
serialized before it crosses the documented interface is not credited as the
returned *type* (a documented whole-result JSON dictionary can still be
`partial`), and a caller-written loop plus list accumulation is `no`, not
`partial`.

## Result summary

| Tool | Value | Confidence | Short basis |
|---|---|---:|---|
| PoolParty | **yes** | high | documented design operations return an inspectable `DnaPool`, which represents the whole lazy library and its provenance DAG |
| VaLiAnT | **no** | high | documented CLI produces CSV/VCF files; its export method returns only counts |
| MPRAnator | **no** | high | documented web/REST result is FASTA text; internal lists are not a documented API value |
| MPRA Design Tools | **partial** | high | documented `processVCF()` returns a list whose `result` member is a whole-library data frame |
| Oligopool Calculator | **partial** | high | documented design functions return the complete output `DataFrame` plus a stats dict |
| Mutation Maker | **partial** | medium | the asynchronous design API's result call returns the complete result as a JSON/dictionary, while its purpose-built internal output type does not cross the interface |
| DNA Chisel | **no** | high | its documented purpose-built object contains one sequence/problem; collection examples build a Python list in a caller loop |
| tangermeme | **partial** | high | documented design calls return all `n_best` designed sequences as a raw `torch.Tensor` |

## Per-tool audits

### PoolParty — **yes** (high confidence)

Primary revision: local read-only repository
`poolparty-statecounter`, commit
`1bb0179e1c3720b1fffd471802b3040f9336de28`.

Evidence:

* The documentation defines the carrier, rather than leaving its meaning to the
  caller: “**A Pool represents a designed collection of DNA sequences. Pools are
  lazy: they record the rules for generating sequences**” and “**every operation
  returns a new Pool**” (`poolparty/docs/quickstart.rst:32-40`).
* The dedicated Pool page makes the whole-library boundary explicit: “**The final
  Pool ... is called the root Pool**”; its rooted DAG describes the library and
  internal bookkeeping (`poolparty/docs/pool.rst:16-19`).
* The object carries identity and membership/design inspection: `name`,
  `num_states`, `parents`, `seq_length`, `regions`, and `has_region` are defined on
  `Pool` (`src/poolparty/pool.py:94-155`); `operation` records the operation that
  created it (`docs/pool.rst:75-81`).
* A documented multi-sequence source has an explicit purpose-built return:
  `from_seqs` says “**A Pool object yielding the provided sequences**” and returns
  `DnaPool(operation=op)` (`src/poolparty/base_ops/from_seqs.py:60-63,95-96`).

Behavioral check, run read-only as required:

```text
PYTHONDONTWRITEBYTECODE=1 poolparty-statecounter/.venv/bin/python -c "..."
poolparty.dna_pool.DnaPool
audit_library 12 4 1 MutagenizeOp
pandas.core.frame.DataFrame 12
```

The test constructed `from_seq('ACGT').mutagenize(num_mutations=1,
mode='sequential').named('audit_library')`; the bound result was a `DnaPool` with
12 members, one parent and a `MutagenizeOp`.  Materializing it separately returned
a 12-row DataFrame.  Thus the qualifying value is the design call's `DnaPool`, not
the weaker export DataFrame.

Disconfirmation attempt: I searched all package source and documentation with
`rg -n "def from_seqs|return .*Pool|class .*Pool|to_df|generate_library"` and read
`quickstart.rst`, `pool.rst`, `pool.py`, `from_seqs.py`, and the generation/export
references.  I specifically checked whether `Pool` was merely a single-sequence
wrapper or whether only `generate_library()` returned the collection.  The cited
root-Pool documentation, exact `num_states`, parent link, and operation provenance
disconfirm both alternatives.

What would change the score: only a general-purpose returned collection would
have produced `partial`; only files or a single-sequence/problem object would have
produced `no`.  I searched both materialization return paths and the object model;
neither weaker condition holds.

### VaLiAnT — **no** (high confidence)

Primary revision: `cancerit/VaLiAnT` `develop`, commit
`8796cc112dafd4919fec59913f58cd2be87c45eb`.

Evidence:

* The official README describes the result exclusively as files: “**Output
  files:**” followed by QC CSV, oligonucleotide metadata CSV, VCF, unique
  oligonucleotide CSV and configuration JSON (`README.md:96-102`).
* The official wiki's **Output files** page likewise enumerates metadata CSV,
  excluded metadata CSV, VCF, and unique-oligonucleotide CSV (wiki
  `Output-files.md:1-14`), and says the unique output is “**A file with
  oligonucleotide name and generated sequence**” (`:32-36`).
* The current implementation confirms that boundary. `MetaTable.to_csv(...)` is
  declared to return `OligoGenerationInfo` (`src/valiant/meta_table.py:340`), opens
  and writes the metadata and variant files (`:441-449`), writes unique sequences
  to CSV (`:676-686`), and returns `info` (`:688`). `OligoGenerationInfo` contains
  only `too_short`, `in_range`, and `too_long` counters
  (`src/valiant/oligo_generation_info.py:26-47`), not sequences.
* The developers' paper says “**VaLiAnT is a command line tool written in Python**”
  and describes the command as generating output files; Figure 3's text says the
  outputs include full-library metadata, VCF, and a unique CSV for synthesis
  (Barbon et al. 2022, abstract and Fig. 3 caption, local PDF
  `Barbon2022_VaLiAnT_all.pdf`).

Disconfirmation attempt: I enumerated the complete `src/valiant/` tree and
stream-grepped every Python file for class declarations, with special attention to
`library|pool|collection|table|output|oligo|experiment|meta`. The plausible hits
were `MetaTable`, `MutatorCollection`, `OligoGenerationInfo`, and `OligoSeq`.
`MutatorCollection` is a collection of mutator algorithms, not output members
(`src/valiant/mutator.py:76-118`); `OligoSeq` is one oligo
(`src/valiant/oligo_seq.py:34-45`); `MetaTable` is the file writer above. I also
checked README **Usage**, README **Output files**, wiki **Home**, wiki **Output
files**, `sge_proc.py`, `cdna_proc.py`, `meta_table.py`, and the paper's
Implementation/Output sections. `run_sge` and `run_cdna` return `None`; per-targeton
processing returns the count object.

What would change the score: a documented Python design call returning all oligos
as a DataFrame/list would produce `partial`; returning a purpose-built library
instance with member access would produce `yes`. I explicitly searched the CLI
entry points, processing returns, export class, wiki, README, and all class names
for either form; neither exists.

### MPRAnator — **no** (high confidence)

Primary revision: `hemberg-lab/MPRAnator` `master`, commit
`9969790d62410138d4281b7955da6d085f07b1bc`.

Evidence:

* The shipped documentation says: “**The result page (plain text view) displays
  the synthesized oligonucleotides in FASTA format**”
  (`iliasApp/templates/iliasApp/docs.html:45-50`; repeated for motif design at
  `:84-93`).
* The web return path constructs `finalOutput` by caller-side list accumulation
  over backgrounds and combinations (`iliasApp/views.py:92-100`), converts it to
  `mpraOutput`, and either writes that text into an `HttpResponse` or places it in
  a download context (`:106-123`). The SNP path similarly writes `mpraOutput` as a
  plain-text response (`:365-387`).
* There is an internal `oligo(...)` helper which returns `allResults`, a Python
  list of dicts (`oligo.py:35-40,74-78`). It is not presented as a documented API;
  the documented interface converts/merges these lists and returns FASTA text.
* Although the developers' paper states that a REST API provides programmatic
  access, the current repository's only named `MpraSnpApiView` returns the literal
  string `testing` (`iliasApp/views.py:419-424`); the actual documented result is
  the FASTA text above (Georgakopoulos-Soares et al. 2017, abstract, local PDF
  `Georgakopoulos-Soares2017_MPRAnator.pdf`).

Disconfirmation attempt: I enumerated the repository and streamed a recursive
`^class ` grep over all top-level and `iliasApp` Python files. The output contained
Django form/model classes and parsing classes (`FastaFile`, `Sequence`,
`SnpFile`), but no result/library class. I read `README.md`, the full shipped
`docs.html` result descriptions, `views.py`, `viewsCore.py`, `oligo.py`,
`iliasApp/tasks.py`, both downloadable-script paths in the tree, and searched the
paper for `download|output|FASTA|result|web|standalone|script`. The internal list
was the strongest in-memory near miss, but it fails the row's documented-API
condition, and combining it across inputs is exactly the caller accumulation
excluded by the global rule.

What would change the score: a documented callable or REST response returning the
whole design as a list/DataFrame/JSON collection would produce `partial`; a
documented purpose-built result/library type would produce `yes`. I searched the
web view, purported REST view, downloadable scripts, docs, paper and all defined
classes for these alternatives; only plain-text FASTA crosses the documented
interface.

### MPRA Design Tools — **partial** (high confidence)

Primary revisions: `andrewGhazi/mpradesigntools` `master`, commit
`afd386ef12051bb0a37ad63a6f92acd555246584`, and companion Shiny repository
`andrewGhazi/designMPRA`, commit `0cf56ee602fc86dde705906d071a39cbdf6e99a8`.

Evidence:

* The package documentation defines `processVCF` as the main whole-design call: it
  turns a VCF into a labeled MPRA sequence set (`man/processVCF.Rd:84-87`).
* Its documented value is exactly the row's `partial` case: “**A list of two
  data_frames. The first, named 'result', is a data_frame containing the labeled
  MPRA sequences. The second, named 'failed' ...**” (`man/processVCF.Rd:78-82`).
* The implementation creates `res = list(result = successes, failed = ...)`
  (`R/processVCFfast.R:1254,1270,1292`) and returns the complete value with
  `return(res)` (`:1341`). Writing a TSV is optional (`outPath = NULL` in the
  signature, `:1099-1115`; optional writes at `:1256-1259,1294-1329`). Thus a user
  can bind and inspect the whole result without a file.
* `NAMESPACE` exports `processVCF` (`NAMESPACE:3`), so this is a documented package
  API, not an example-only helper.

Disconfirmation attempt: I enumerated all R and man-page files and read
`README.md`, `man/processVCF.Rd`, `man/processSnp.Rd`, `NAMESPACE`, and the full
`processVCFfast.R` return construction. I recursively searched both repositories'
R sources for `setClass|setRefClass|R6Class|structure\(|class\s*\(`. No S3/S4/R6
library/result class is constructed; the only `class(...)` uses in the package are
type checks. The returned list and tibbles contain useful labels, but their types
are general-purpose and their library meaning is not embodied in a tool-defined
class.

What would change the score: a documented `MPRA_library`/equivalent class carrying
membership or provenance would produce `yes`; only an output TSV with no returned
whole DataFrame/list would produce `no`. I searched package exports,
documentation, both codebases, and the actual final expression for both forms.

### Oligopool Calculator — **partial** (high confidence)

Primary revision: `ayaanhossain/oligopool` `master`, commit
`b88fa394ca67ed4c48ec41127b5707694ee7cf0a`.

Evidence:

* The official API reference documents the whole barcode-design call as
  `df, stats = op.barcode(...)` (`docs/api.md:62-92`) and states “**Returns:
  `(DataFrame, stats_dict)`**” (`:117`). Its input is the annotated whole oligo-pool
  DataFrame and its output adds the designed barcode column to all variants
  (`:96-101`).
* The same documented general-container return is present for primer, motif,
  spacer, merge, reverse-complement, join, final, split, pad, compress and expand
  (the `**Returns**` declarations at `docs/api.md:200,279,353,456,507,560,613,
  675,743,810,868`). This is a deliberate DataFrame pipeline, not an accidental
  internal variable.
* The barcode implementation documents “**A pandas DataFrame of generated
  barcodes**” and a stats dict (`oligopool/barcode.py:62-64`), annotates the return
  as `Tuple[pd.DataFrame, dict]` (`:17-34`), and returns `(outdf_return, stats)` on
  success (`:1074-1083`). `output_file` is “optional in library usage” (`:47-50`),
  so the in-memory whole result does not depend on file output.

Disconfirmation attempt: I read README **Getting Started**, the API reference's
complete table of contents and every `**Returns**` declaration, plus
`oligopool/barcode.py`'s validation, success, failure, write and return paths. I
enumerated all `oligopool/` and `oligopool/base/` Python files and streamed a
recursive `^class ` grep. The defined classes were `Scry`, `SafeCounter`,
`SafeQueue`, `FailureSampler`, `vectorDB`, and CLI parser/formatter classes; none is
a sequence-library/result class. (`vectorDB` is a k-mer screening database, not the
designed oligo collection.) An attempted GitHub code-search query for `class` was
unauthorized/null, so I replaced it with the complete streamed archive grep rather
than treating the failed query as negative evidence.

What would change the score: a documented tool-defined library result class would
produce `yes`; stats-only/file-only design would produce `no`. I checked all
documented module return schemas, the representative design implementation, and
all class declarations for these alternatives. The documented DataFrame is exactly
the binding definition of `partial`.

### Mutation Maker — **partial** (medium confidence)

Primary revision: `ra100/Mutation_Maker` `master`, commit
`396c1c0ede7529f3dedf65e821e8c1f20c9a7043` (the paper's printed Merck repository
URL is dead; this is the supplied surviving repository).

Evidence:

* The documented product is an application: the README calls it an
  “**Application for mutagenic primer design**” (`README.md:1-4`), documents the
  webserver and API service addresses (`:33-52`), and describes the API as the
  server which “accepts tasks from the frontend and schedules them”
  (`:207-218`). It does not document a Python design call returning a library
  instance.
* The REST submit endpoint `/qclm` calls `start_celery_task`
  (`api/api.py:48-53`), which returns a general dictionary containing a
  `result_url` (`:89-105`). Feeding this tool-produced result URL into the shipped
  result endpoint is tool-provided composition under the global rule. That result
  endpoint returns the Celery task's complete serialized result (`:123-129`).
* Internally, `QCLMOutput(JsonObject)` holds `input_data`,
  `full_sequence`, `goi_offset`, and a list of mutation/primer results
  (`backend/mutation_maker/qclm_types.py:119-133`). The workflow entry point
  returns `output.to_json()` (`backend/mutation_maker/qclm.py:92-97`), and the
  project's own tests index that returned whole value as a dictionary, including
  `qclm_result["results"]` (`backend/tests/unit_tests/test_qclm_overlap.py:52,
  91,139,515`). This establishes the general-purpose whole-result container
  required for `partial`, but the `QCLMOutput` instance itself does not cross the
  interface, preventing `yes`.
* The frontend binds the result call's `response.data` in memory
  (`frontend/src/services/api.ts:353-354`) and casts the complete response to
  `QCLMResponseData` (`:394-407`), whose `results` member is a list of mutation
  records (`frontend/src/shared/lib/Api.ts:227-238`). It then converts that JSON
  into ordinary TypeScript records and arrays
  (`frontend/src/shared/lib/ResultData.ts:35-64,84-106`). The developers' paper
  describes the API as communication between UI and worker and Redis as temporary
  storage for the task/results, while describing visualization, printing, export
  and result-sharing—not a documented user Python object (Hiraga et al. 2021,
  Discussion and Container-Based Architecture, local PDF
  `Hiraga2021_MutationMaker.pdf`).

Disconfirmation attempt: I enumerated the complete repository tree and read
`README.md`, `backend/README.md`, `api/api.py`, `backend/tasks.py`, `qclm.py`,
`qclm_types.py`, `ssm.py`, `ssm_types.py`, `pas_types.py`, and the frontend API
result types, `frontend/src/services/api.ts`, and the QCLM/PAS/SSM tests. I searched
those sources and the paper for `API|output|result|return|JsonObject|Output|results|
primers|FASTA|Excel|download|response.data|to_json`. I attacked `partial` in both
directions: the result is stronger than a mere file because the documented/shipped
API and frontend bind the complete dictionary in memory, but weaker than `yes`
because the purpose-built `QCLMOutput`, `SSMOutput`, and `PASOutput` schemas are
serialized before the client receives them.

What would change the score: a documented design boundary returning the actual
`QCLMOutput`/`SSMOutput`/`PASOutput` instance would produce `yes`; only task-status
URLs or XLSX/FASTA output, with no complete result response, would produce `no`. I
searched the main README, backend README, route implementations/docstrings, solver
entry points, frontend retrieval/contract, tests and paper for both alternatives.
The general dictionary return is affirmative, so `no` is not tenable; no
purpose-built instance crosses the boundary, so `yes` is not tenable. Confidence
is medium because the REST contract is documented mostly by shipped route code and
frontend types rather than a standalone API manual.

### DNA Chisel — **no** (high confidence)

Primary revision: `Edinburgh-Genome-Foundry/DnaChisel` `master`, commit
`68c09304341c3656f3dfe63eda37757d6a7b3917`.

Evidence:

* The dedicated class is purpose-built but explicitly singular. Its source module
  says the whole problem consists of “**sequence, constraints, objectives**”
  (`dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py:1-5`); the class
  accepts `sequence`, “**A string of ATGC characters**” (`:47-52`), and stores one
  `self.sequence` (`:115-139`).
* The official usage gets one final sequence: “**GET THE FINAL SEQUENCE (AS STRING
  OR ANNOTATED BIOPYTHON RECORDS)**”, followed by `problem.sequence` and
  `problem.to_record(...)` (`README.rst:79-82`). The row definition explicitly
  says a single-sequence or single-optimization object is `no`, however
  purpose-built.
* The official genome-wide example proves that a collection is caller
  reconstruction: `optimized_records = []`, a `for feature` loop creates one
  `DnaOptimizationProblem`, and the caller executes
  `optimized_records.append(problem.to_record(...))`
  (`examples/common_scenarios/genome-wide-optimization.py:9-29`).
* The primer-collection example is even more explicit: “**DNA Chisel is not
  originally meant for creating collections of sequences**” and creates them “one
  after the other” (`examples/common_scenarios/primers_collection.py:1-8`); the
  caller initializes `existing_primers = []`, loops, and appends each returned
  string (`:39-46`). Under the global rule, that cannot earn `partial`.

Disconfirmation attempt: I enumerated the full repository tree and read
`README.rst`, `docs/ref/core_classes.rst`, the complete `DnaOptimizationProblem`
constructor/object model, `genome-wide-optimization.py`, and
`primers_collection.py`. I searched official docs/source/examples for
`collection|library|batch|DnaOptimizationProblem|return` and checked the reports
API; `optimize_with_report(target='@memory')` can return report archive bytes, but
still for one problem/sequence (`DnaOptimizationProblem.py:202-219`). No API call
returns the caller-accumulated collection.

What would change the score: a documented batch call returning a list/DataFrame of
all optimized sequences would produce `partial`; a documented purpose-built
collection/library class would produce `yes`. I searched the core class docs,
source tree, collection examples and report return for both. The official example's
manual list append affirmatively establishes their absence from the tool-provided
interface.

### tangermeme — **partial** (high confidence)

Primary revision: `jmschrei/tangermeme` `main`, commit
`2006b310cd72a28c56c3ea4ba67f738fff74bb89`.

Evidence:

* The official API page documents `screen`, `greedy_substitution`,
  `beam_substitution`, and `greedy_marginalize` as the design API
  (`docs/api/design.rst:1-14`).
* `beam_substitution` is a whole multi-sequence design return: `n_best` is “**The
  number of sequences to return at the end**”
  (`tangermeme/design/beam_substitution.py:125-135`), and its documented return is
  “**X: torch.Tensor, shape=(n_best, len(alphabet), length) — The designed
  sequences, ranked from lowest loss to highest loss.**” (`:173-176`). The
  implementation returns one concatenated tensor containing all retained designs
  (`:313-314`).
* Independently, `screen` describes screening a “large pool” and keeping `n_best`
  hits (`tangermeme/design/screen.py:37-48,96-104`); it documents and returns a
  tensor of shape `(n_best, len(alphabet), length)` (`:140-143,183-185`).
* The tensor has no member identities, provenance, or library methods. It is the
  exact general-purpose container allowed by `partial`.

Disconfirmation attempt: I read `docs/api/design.rst` and the complete return and
input sections/implementations of `greedy_substitution.py`,
`beam_substitution.py`, and `screen.py`. The single-input restriction of greedy
design (`greedy_substitution.py:63-67,155-159`) does not erase the stronger
multi-output `beam_substitution` and `screen` calls. I enumerated the full tree and
stream-grepped all 28 package Python files for `^class `. Exactly six classes were
found: five analysis/result `NamedTuple`s (`PerturbationResult`,
`PerturbationAnnotationsResult`, `AttributionReferencesResult`,
`SaturationMutagenesisRawResult`, `SpaceResult`) and `TangermemeWarning`
(`results.py:16,28,41`; `saturation_mutagenesis.py:19`; `space.py:22`;
`utils.py:18`). None represents designed sequences or a library.

What would change the score: a documented tool-defined designed-library class
carrying identity/provenance/member operations would produce `yes`; a design call
returning only one sequence or only predictions/files would produce `no`. I
searched every class plus every documented design return. Multiple designed
sequences are genuinely returned, but only as a tensor, fixing the value at
`partial`.

## Search ledger

This ledger records all searches used for the row, including unsuccessful ones,
so a referee can rerun or widen them.

1. Read `revision/tables/ROW_DEFINITIONS.md:1-43` (preamble/global rule and row 1)
   before inspecting tools, then stated the operational test above.
2. Listed `revision/tool_survey/final/*.md` and read the eight named lead records
   only to locate primary repositories/files. Listed local PDFs and searched for
   local repository clones; only PoolParty was local.
3. PoolParty: `rg` over `src/poolparty` and `docs` for
   `def from_seqs|return .*Pool|class .*Pool|to_df|generate_library`; read the files
   named in its audit and ran the recorded venv behavior check with bytecode writes
   disabled.
4. For every remote repository, queried GitHub's repository endpoint for the
   default branch, enumerated the recursive Git tree, and pinned the branch commit
   SHA. Relevant tree filters were `README|docs|api|output|\.py$|\.R$|man/|vignette`.
5. VaLiAnT: read/grepped README, `sge_proc.py`, `cdna_proc.py`, `meta_table.py`,
   `oligo_generation_info.py`, `oligo_seq.py`, `mutator.py`; read official wiki
   **Home** and **Output files**; streamed all `src/valiant/*.py` class declarations
   and filtered for `library|pool|collection|table|output|oligo|experiment|meta`.
6. MPRAnator: read/grepped README, shipped `docs.html`, `views.py`, `viewsCore.py`,
   `tasks.py`, `oligo.py`; enumerated downloadable scripts; streamed all top-level
   and `iliasApp` class declarations. The first archive glob incorrectly included
   an extra `MPRAnator/` path and returned no files; it was corrected to `*/*.py`
   and `*/iliasApp/*.py`, producing the class list summarized above.
7. MPRA Design Tools: read README, `processVCF.Rd`, `processSnp.Rd`, `NAMESPACE`,
   and `processVCFfast.R`; searched the implementation for
   `processVCF|results = list|return(res)|outPath|write_tsv`; recursively searched
   both R repositories for `setClass|setRefClass|R6Class|structure\(|class\s*\(`.
8. Oligopool Calculator: read README, full API contents/return declarations, and
   `barcode.py`; enumerated source; archive-grepped all Python class declarations.
   GitHub code search for `repo:ayaanhossain/oligopool "class " language:Python`
   returned an unauthorized/null result and was not used; the complete archive
   stream grep replaced it.
9. Mutation Maker: read/grepped main and backend READMEs, `api.py`, `tasks.py`,
   QCLM/SSM/PAS solver and type files, frontend `Api.ts`, `ResultData.ts`,
   `services/api.ts`, and QCLM/PAS/SSM tests; tree filters included backend,
   frontend, docs, examples and all Python files; searched for
   `API|output|result|return|JsonObject|Output|results|primers|FASTA|Excel|download|
   response.data|to_json`.
10. DNA Chisel: read/grepped README, core class docs/source, genome-wide and primer
    collection examples; recursive tree filters covered `dnachisel/`, `docs/`, and
    `examples/`; search terms were
    `collection|library|batch|DnaOptimizationProblem|return|sequence`.
11. tangermeme: read/grepped design API and the three design implementations;
    enumerated `docs/` and all package files; streamed `^class ` across all Python
    modules with file/line labels.
12. Extracted all three especially relevant local PDFs (VaLiAnT, MPRAnator,
    Mutation Maker) with PyMuPDF as instructed and searched respectively for
    `output|CSV|command.line|Python|library`,
    `download|output|FASTA|result|web|standalone|script`, and
    `download|output|FASTA|Excel|result|API|web interface`. Repository/docs evidence
    was preferred where it directly settled the return type. The remaining papers
    were not needed to override a documented implementation return.

## Row-level finding

The row is applicable on one consistent scale and does discriminate this set:
**1 yes / 4 partial / 3 no**. Its useful separation is not “Python versus web” or
“object-oriented versus functional”; it is whether the documented design boundary
returns (i) a purpose-built whole-library carrier, (ii) only a general batch
container, or (iii) no whole-library in-memory result. The strict single-sequence
exclusion and the global ban on caller accumulation are essential to keeping DNA
Chisel and MPRAnator on the same threshold as the other tools.
