# `library_object` audit

> **SUPERSEDED — do not use these values.** Scored against the earlier wording of
> row 1's `partial` clause, which had no floor and produced a threshold shift:
> five of eight cells moved versus the prior scoring, four of them the same way.
> The clause has since been tightened (whole library + returned by the tool's own
> API; caller-assembled collections are `no`). This row is being **rescored in
> batch 3**; use that result instead.
>
> Retained for provenance, and because its DNA Chisel reasoning is what the
> tightened wording codifies.

Assessed 2026-08-15. This audit applies only the definition in
`revision/tables/ROW_DEFINITIONS.md` §1. The records under
`revision/tool_survey/final/` were treated as leads only; none is evidence here and
none of their ratings was used as a baseline.

## Operational test

For each tool, locate a documented API call that creates or returns the complete
multi-sequence design result, then verify from its documented return contract and
implementation whether the bindable in-memory value is (a) an instance of a
tool-defined library-specific type carrying member identity, design provenance, or
membership operations (`yes`), (b) a durable general-purpose container
(`partial`), or (c) only files/no such return (`no`); search the API, documentation,
and source for a qualifying alternative before assigning `no`.

This is one threshold for all tools. In particular, a custom object for **one**
sequence is not a whole-library object; a DataFrame/tibble, tensor, list/dict, HTTP
response, or decoded JSON object is `partial` even when it has excellent columns or
metadata; and an internal custom result model does not earn `yes` if the documented
user-facing call serializes it to a general-purpose value.

## Result

| tool | value | confidence |
|---|---|---|
| PoolParty | **yes** | high |
| VaLiAnT | **no** | high |
| MPRAnator | **partial** | medium |
| MPRA Design Tools | **partial** | high |
| Oligopool Calculator | **partial** | high |
| Mutation Maker | **partial** | medium |
| DNA Chisel | **no** | high |
| tangermeme | **partial** | high |

Repository snapshots used were PoolParty `1bb0179` (the supplied read-only local
checkout), VaLiAnT `8796cc1` (`develop`), MPRAnator `9969790`, mpradesigntools
`afd386e`, designMPRA `0cf56ee`, oligopool `b88fa39`, Mutation Maker `396c1c0`,
DNA Chisel `68c0930`, and tangermeme `2006b31`. External repositories were shallow
cloned read-only into `/tmp`; nothing was installed.

## Per-tool audit

### PoolParty — **yes** (high confidence)

**Qualifying return.** The documentation defines the carrier for precisely this
purpose:

> “A **Pool** represents a designed collection of DNA sequences. A Pool can
> represent the final library you wish to generate, or an intermediate set of
> sequences used to construct it.”

Source: `poolparty-statecounter/poolparty/docs/pool.rst:4-6`.

The documented multi-member constructor is bindable and returns that type:

```python
def from_seqs(...) -> DnaPool:
    """Create a Pool containing the specified sequences."""
    ...
    result_pool = DnaPool(operation=op)
    return result_pool
```

Source: `poolparty-statecounter/poolparty/src/poolparty/base_ops/from_seqs.py:16-30,95-96`.
Its return documentation says:

> “A Pool object yielding the provided sequences using the specified selection
> mode.”

Source: `poolparty-statecounter/poolparty/src/poolparty/base_ops/from_seqs.py:60-63`.

The object carries member identity and inspectable design state rather than merely
wrapping a list. `from_seqs` accepts `seq_names` and design-card keys `seq_name` and
`seq_index` (`from_seqs.py:21-22,42-45,57-58,99-103`). `Pool` exposes `name`,
`num_states`, `parents`, `seq_length`, and `regions`
(`src/poolparty/pool.py:94-155`), and its `operation` records the operation that
created it (`docs/pool.rst:75-81`). The docs explicitly demonstrate binding and
inspection:

```python
library = (
    pp.from_seq("ATCGATCG")
    .mutagenize(num_mutations=1, mode="sequential")
    .named("mutants")
    .print_library(num_seqs=3)
    .repeat(times=2)
)
print(library.num_states)  # 48
```

Source: `docs/quickstart.rst:127-162`; the same page states “At any point you can
check `pool.num_states`, `pool.seq_length`, and `pool.regions`” and use
`pool.print_dag()` (`:176-178`).

**Behavioral check.** Run read-only with the supplied interpreter and
`PYTHONDONTWRITEBYTECODE=1`:

```python
pp.init()
library = pp.from_seqs(
    ['AAA', 'CCC'], seq_names=['a', 'c'], mode='sequential',
    cards=['seq_name', 'seq_index']).named('audit_library')
```

Observed type and inspection output:

```text
poolparty.dna_pool.DnaPool
DnaPool(id=0, name='audit_library', op='op[0]:from_seqs', num_states=2)
audit_library 2 3 [] from_seqs
```

Materializing it showed two named rows and the requested per-member cards. This
settles that the value bound to `library` is the tool-defined carrier, not the later
DataFrame export.

**Disconfirmation attempt and counterfactual.** A `partial` would have resulted if
the documented construction call returned only `pd.DataFrame`, list, or dict; a
`no` would have resulted if it only wrote a file. I searched all documented Pool
creation and export paths for those alternatives. The export call does return a
DataFrame, but it is downstream of the already-held `Pool`; docs say “every
operation returns a new Pool” (`docs/quickstart.rst:91-92,132-144`). No contrary
constructor was found.

**Searches performed.** `rg -n` over `docs/{index,quickstart,pool}.rst` and
`src/poolparty/` for `class (Pool|DnaPool)`, `from_seq`, `from_seqs`, `return .*Pool`,
`num_states`, `generate_library`, and `to_df`; read `pool.py:1-220`,
`from_seqs.py:1-105`, `quickstart.rst:40-190`, and `pool.rst:1-125`; ran the
behavioral check above. The supplied repository remained read-only.

### VaLiAnT — **no** (high confidence)

**Observed public contract.** VaLiAnT documents a command/file interface. Its
README enumerates CSV, VCF, and JSON **output files** (`README.md:96-102`) and makes
`OUTPUT` an output-directory argument (`README.md:211-216`). The authors' paper is
equally explicit:

> “VaLiAnT is run from the command line to generate output files, including
> ‘meta’ full library metadata, ‘VCF’ of all variants generated and a ‘unique’ csv
> file for easy ordering of sequence synthesis.”

Source: `papers/Barbon2022_VaLiAnT_all.pdf`, Figure 1 caption; see also §2.1,
“VaLiAnT is run from the command line,” and §2.3, “standalone executable Python
package exposing a command line interface.”

The source confirms that the apparent in-memory `MetaTable` is an internal writer,
not a returned library. `run_sge` is annotated `-> None`
(`src/valiant/sge_proc.py:527`), builds an in-memory database, processes targetons,
and finalizes (`:534-568`). Per targeton it does:

```python
# Write metadata files
return mt.to_csv(conn, targeton_name)
```

Source: `src/valiant/sge_proc.py:351-368`. `MetaTable.to_csv` returns only
`OligoGenerationInfo` statistics, opens the metadata/VCF paths for writing, writes
the unique CSV, and returns `info` (`src/valiant/meta_table.py:340-355,440-449,
676-688`). `run_cdna` is also `-> None` (`src/valiant/cdna_proc.py:137`). Thus the
complete set of oligos crosses the user boundary only as files.

**Disconfirmation attempt and counterfactual.** `partial` would require a documented
call returning the complete set as a list/dict/DataFrame; `yes` would require a
returned VaLiAnT library class. I searched the complete `src/valiant` tree and
README for `class .*Library`, `class .*Pool`, `class .*Collection`, `def .*library`,
`return .*library`, `DataFrame`, `pandas`, `to_csv`, and `output files`. The only
name-like hit was `MutatorCollection`, which is a collection of mutator actions,
not output oligos. The only whole-result path found was `MetaTable.to_csv`, whose
return is statistics. Therefore this is an observed absence, not a failure to look.

**Searches/files/pages checked.** Read `README.md:1-115,200-235`,
`src/valiant/{sge_cli,sge_proc,cdna_proc,meta_table,db}.py`, `pyproject.toml`, and the
paper abstract, §§2.1, 2.3, 2.5.3-2.5.4, Figure 1 caption, and Supplementary Table
2. Ran `rg -ni` over README, `src`, examples, and configuration for `library`,
`return`, `to_csv`, `write`, `output`, `class .*library`, `dataframe`, `pandas`,
`sqlite`, and `api`, followed by the type/class searches listed above.

### MPRAnator — **partial** (medium confidence)

**Qualifying general-purpose value.** The paper documents programmatic access:

> “The REST API allows programmatic access to MPRAnator using simple URLs.”

Source: `papers/Georgakopoulos-Soares2017_MPRAnator.pdf`, abstract.

The repository ships the client script served for download by the application
(`iliasApp/views.py:393-412`). It binds the complete response before writing it:

```python
result = requests.post(url = url, data = data)

# writing the results
openfile = open(fileToWriteResults,"w")
openfile.write(result.content)
```

Source: `downloadables/MpraMotifs_script.py:68-74`; the SNP client is identical in
this respect (`downloadables/MpraSnps_script.py:58-64`). On the server, download
mode returns a plain-text `HttpResponse` containing `mpraOutput`
(`iliasApp/views.py:106-115` for motif designs and `:372-381` for SNP designs).
Consequently a user can retain the complete FASTA-like library in the ordinary
`requests.Response.content` bytes value. The carrier is general-purpose and has no
MPRAnator membership/provenance methods, so the row's explicit `partial` rule
applies. Because the positive API evidence is a shipped/downloadable example
script rather than a separately documented Python package API, it could not score
higher than `partial` in any event.

**Disconfirmation attempt and counterfactual.** `yes` would require the API to
return a tool-defined result/library type. I searched every Python class and return
path. `mycustom.FastaFile` is tool-defined and has membership operations, names,
and lengths (`mycustom.py:16-120`), but its constructor documentation identifies it
as a FASTA file/input parser and the views instantiate it only for input
(`iliasApp/views.py:70-83,172-178,352-357`); no design call returns it. Internally,
`oligo.oligo` returns `allResults`, a list of dictionaries with `header`, `sequence`,
and `sequenceHTML` (`oligo.py:35-78`), and final output is flattened to text.
Conversely, `no` would require there to be no bindable whole result; the shipped
REST client disproves that by binding `result.content`. I therefore searched for
both counterfactuals and assigned `partial`.

**Searches/files/pages checked.** Read the complete one-line `README.md`,
`downloadables/{MpraMotifs_script,MpraSnps_script}.py`, `iliasApp/{urls,views}.py`,
`viewsCore.py`, `oligo.py`, `part1.py`, and `mycustom.py`, plus the paper abstract.
Ran `rg -ni` across all Python files for `library`, `output`, `return`, `DataFrame`,
`pandas`, `class`, `write`, `download`, `api`, and `usage`, then targeted
`MpraSnpApiView`, `HttpResponse`, `createMPRAResultOutput`, and `oligo(`.

### MPRA Design Tools — **partial** (high confidence)

**Qualifying general-purpose value.** `processVCF` is the documented package API:

> “Currently the main function of MPRA Design Tools package is to design a set of
> barcoded sequences for MPRA experiments ... This is done with the `processVCF`
> function.”

Source: `mpradesigntools/README.md:41-47`.

Its generated R manual gives the return contract verbatim:

> “A list of two data_frames. The first, named 'result', is a data_frame
> containing the labeled MPRA sequences. The second, named 'failed', is a
> data_frame listing the SNPs that are not able to have MPRA sequences generated
> and the reason why.”

Source: `mpradesigntools/man/processVCF.Rd:78-82`; same Roxygen text at
`R/processVCFfast.R:1068-1071`.

The output path is optional (`man/processVCF.Rd:76`), and implementation constructs
`res = list(result = successes, failed = ...)` and `return(res)`
(`R/processVCFfast.R:1241-1270,1281-1297,1333-1341`). Thus the complete labeled
design persists in memory even without a file, but its carrier is an ordinary R
list of tibbles/data frames. That is exactly `partial` under the instrument.

**Disconfirmation attempt and counterfactual.** `yes` would require a package-
defined MPRA library class rather than a list/tibble; `no` would require file-only
output. I searched both the package and companion Shiny repository for S3/S4 class
definitions and constructors (`setClass`, `setOldClass`, `structure`, class
assignment, `S3method`) and for all documented returns. None defines a library
class. `NAMESPACE:3-4` exports only `processVCF` and `spread_and_fix_indels`; the
companion `designMPRA/scripts/processVCFfast.R:430-432` likewise returns a general
list while optionally writing. The optional `outPath` and explicit `return(res)`
rule out `no`.

**Searches/files/pages checked.** Read `mpradesigntools/README.md`, `DESCRIPTION`,
`NAMESPACE`, `man/processVCF.Rd`, and `R/processVCFfast.R:1000-1100,1230-1350`.
Searched all `R/`, `man/`, and metadata for `library`, `return`, `data.frame`,
`list(`, `output`, `write`, `class`, `setClass`, `@return`, and `usage`. In the
companion `designMPRA` repo, searched `server.R`, `ui.R`, and `scripts/*.R` for the
same terms plus `reactive`, `renderTable`, and `download`; inspected its
`scripts/processVCF{,fast}.R` return paths.

### Oligopool Calculator — **partial** (high confidence)

**Qualifying general-purpose value.** The package-level manual makes the public
contract unambiguous:

```python
>>> import oligopool as op
>>> df, stats = op.barcode(input_data='variants.csv', ...)
```

and:

> “Modules return (DataFrame, stats). Chain them iteratively; use patch_mode=True
> to extend pools without overwriting existing designs.”

Source: `oligopool/__init__.py:130-136`.

For a concrete design call, `barcode` is annotated
`-> Tuple[pd.DataFrame, dict]` (`oligopool/barcode.py:17-34`) and documents:

> “A pandas DataFrame of generated barcodes; saves to `output_file` if specified.”

Source: `oligopool/barcode.py:62-64`. `output_file` is optional for library use
(`:47-49`), the input must carry a unique `ID` (`:67-68`), and the implementation
returns `(outdf_return, stats)` (`:1074-1083`). The DataFrame therefore is durable,
bindable, member-identified, and represents the whole current oligo-pool design,
but remains a general pandas type. This satisfies `partial`, not `yes`.

**Disconfirmation attempt and counterfactual.** `yes` would require the design
functions to return an Oligopool-defined library carrier. I searched all public API
modules and docs for `class .*Pool`, `class .*Library`, dataclasses, NamedTuples,
Protocols, and return annotations. The only public tool-defined classes found were
`vectorDB` and `Scry` (`oligopool/base/vectordb.py`, `base/scry.py`); they are a
k-mer database and classifier, not the carrier returned by design functions.
Every documented design/assembly return in `docs/api.md` is a DataFrame plus stats
(for example the repeated contracts at `docs/api.md:117,200,279,353,456,507,560,
613,676,743`). `no` would require only files, which the optional output path and
actual DataFrame return disprove.

**Searches/files/pages checked.** Read `README.md`, `docs/{README,docs,api}.md`,
`oligopool/__init__.py`, and the public design modules `barcode`, `primer`, `motif`,
`spacer`, `merge`, `join`, `revcomp`, `final`, `split`, `pad`, `compress`, and
`expand`. Ran `rg -ni` over repository docs/source for `library`, `return`,
`DataFrame`, `pandas`, `dict`, `list`, `output`, `write`, `class .*pool`, and `API`,
then targeted all class/NamedTuple/dataclass definitions and all
`Tuple[pd.DataFrame` annotations.

### Mutation Maker — **partial** (medium confidence)

**Qualifying general-purpose value.** Mutation Maker's documented architecture has
an API between the UI and worker and a Redis in-memory store for “the design task
and its respective results”; its results page provides tables and exports
(`papers/Hiraga2021_MutationMaker.pdf`, Methods, “Software Architecture” and “User
Interface”). The repository's API returns a result URL for each task
(`api/server_fastapi.py:43-50,53-65`) and the result endpoint returns the completed
Celery result (`:111-117`). The frontend binds it directly:

```typescript
export const fetchJobResultData = (jobId: string): Promise<any> =>
  Axios.get(`${API_PREFIX}/result/${jobId}`).then((response) => response.data)
```

Source: `frontend/src/services/api.ts:353-354`.

The ordinary returned object does represent the whole design and carries rich
identity/provenance. For example, the public response shape has `input_data`,
`results: QCLMMutationData[]`, `full_sequence`, and `goi_offset`
(`frontend/src/shared/lib/Api.ts:227-238`); each result carries mutation identities
and primers (`:227-230`). PAS similarly exposes `results` whose fragments contain
their oligos (`:287-328`). These are TypeScript structural types over decoded JSON,
not a durable Mutation Maker runtime library class, hence `partial`.

**Why the internal custom types do not make this `yes`.** Internally there really
are purpose-defined schemas:

```python
class QCLMOutput(JsonObject):
    input_data = ObjectProperty(QCLMInput, required=True)
    full_sequence = StringProperty(required=True)
    goi_offset = IntegerProperty(required=True)
    results = ListProperty(QCLMMutationOutput, required=True)
```

Source: `backend/mutation_maker/qclm_types.py:129-133`; corresponding `SSMOutput`
at `ssm_types.py:248-260` and `PASOutput` at `pas_types.py:213-216`. The internal
solver returns `QCLMOutput` (`qclm.py:285,597-675`), but the worker-facing callable
immediately returns `output.to_json()` (`qclm.py:92-97`; SSM does the same at
`ssm.py:55-70`). The tests consume this value by general dictionary operations such
as `result["results"]` and `.keys()`
(`backend/tests/unit_tests/test_pas_output.py:45-67`). The tool-defined class does
not survive the documented API boundary.

**Disconfirmation attempt and counterfactual.** `yes` would result if the
documented user-facing call returned `QCLMOutput`, `SSMOutput`, or `PASOutput`
itself so callers could inspect that purpose-defined instance. I searched all
output schemas, solver return annotations, Celery tasks, HTTP routes, frontend API
types, and result-fetch code; the serialization step above prevents that. `no`
would result if only XLSX/PDF files were exposed, but `/result/{task_id}` and the
bindable decoded response disprove it. Confidence is medium rather than high
because the repository contains both an internal Python object boundary and a
public JSON boundary; the operational test explicitly selects the latter.

**Searches/files/pages checked.** Read `README.md:180-223`, `api/{api,
server_fastapi}.py`, `backend/tasks.py`, `backend/mutation_maker/{qclm,qclm_types,
ssm,ssm_types,pas,pas_types,pas_output}.py`, `frontend/src/services/api.ts`,
`frontend/src/shared/lib/Api.ts`, and output tests. Searched the repository for
`library`, `return`, `output`, `result`, `class .*Output`, `class .*Result`,
`NamedTuple`, `dataclass`, `BaseModel`, `API`, `json`, `to_json`, and output model
constructors. Also checked the paper's architecture, UI/output, and export
descriptions.

### DNA Chisel — **no** (high confidence)

**Observed object is single-sequence, not whole-library.** DNA Chisel has an
excellent purpose-defined object, but it represents one optimization problem for
one sequence. Its class documentation says:

> “Problem specifications: sequence, constraints, optimization objectives.”

and its `sequence` parameter is “A string of ATGC characters”
(`dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py:20-30,47-52`). The
main documented example binds `problem = DnaOptimizationProblem(...)`, then obtains
only:

```python
final_sequence = problem.sequence  # string
final_record = problem.to_record(with_sequence_edits=True)
```

Source: `README.rst:57-82`. No call in this path returns a multi-sequence design.

The strongest contrary example is explicit that collection handling is outside
the abstraction:

> “DNA Chisel is not originally meant for creating collections of sequences ...
> it is still possible to create collections of inter-compatible sequences.”

Source: `examples/common_scenarios/primers_collection.py:1-10`. It realizes that
possibility with a caller-authored `existing_primers = []` loop and `.append`
(`:39-46`). The genome-wide example likewise creates `optimized_records = []`,
constructs one `DnaOptimizationProblem` per feature, appends records, and writes
the list (`examples/common_scenarios/genome-wide-optimization.py:9-32`). A
caller-created list does not satisfy the row.

**Disconfirmation attempt and counterfactual.** `partial` would require a
documented batch call returning the complete collection in a list/DataFrame or
other general container; `yes` would require a tool-defined collection/library
object. I searched README, full docs, all `dnachisel/`, and examples for
`class .*Library`, `class .*Pool`, `class .*Collection`, `def .*library`,
`return .*library`, `list[DnaOptimizationProblem]`, `collection of sequences`,
`multiple sequences`, and collection/batch return paths. No such call or type was
found. The two collection examples above are the only strong hits and both require
the caller to perform the loop and own the list. A related-project note says an
external project, `dnachisel-dtailor-mode`, brings generation of large collections
to DNA Chisel (`README.rst:212-218`), corroborating rather than filling the gap.

**Searches/files/pages checked.** Read `README.rst`, `docs/index.rst`,
`docs/ref/{core_classes,biotools,reports}.rst`, the full
`DnaOptimizationProblem` class and record-representation paths, and the collection,
genome-wide, and manuscript examples. Ran the class/library/batch searches above
and a second search for `DnaOptimizationProblem(`, `problem.sequence`, `to_record`,
`SeqIO.write`, and user loops. Checked the paper's Usage and Implementation
sections, which define “an optimization problem” around “a starting linear or
circular sequence.”

### tangermeme — **partial** (high confidence)

**Qualifying general-purpose value.** `beam_substitution` is a documented design
call that can return several ranked designs via `n_best`. Its signature returns
`torch.Tensor` (`tangermeme/design/beam_substitution.py:25-43`), and its return
contract is:

> “X: torch.Tensor, shape=(n_best, len(alphabet), length) — The designed
> sequences, ranked from lowest loss to highest loss.”

Source: `tangermeme/design/beam_substitution.py:131-135,173-176`. Implementation
returns `torch.cat([...], dim=0)` (`:313-314`). Thus one call yields the complete
multi-sequence result in memory, but only as a general PyTorch tensor with no
member IDs or construction provenance. The independent `screen` design API has
the same contract: `-> torch.Tensor` and returns the `n_best` screened examples
(`tangermeme/design/screen.py:20-44,96-105,140-143,183-185`). This meets the
instrument's exact `partial` condition.

**Disconfirmation attempt and counterfactual.** `yes` would require a tangermeme
library/design-result type containing the designed sequences plus identity,
provenance, or membership operations. I searched all class, NamedTuple, dataclass,
and return definitions. Tangermeme does define `PerturbationResult`,
`PerturbationAnnotationsResult`, `AttributionReferencesResult`, `SpaceResult`, and
`SaturationMutagenesisRawResult`, but their fields hold predictions,
attributions/references, or perturbation outputs; for example
`PerturbationResult` contains only `y_before` and `y_after`
(`tangermeme/results.py:16-25`). None is returned by the design functions as a
sequence-library carrier. `no` would require no multi-sequence in-memory result;
`beam_substitution(n_best>1)` and `screen(n_best>1)` directly disprove that.

**Searches/files/pages checked.** Read `README.md`, `docs/api/{design,ersatz,
saturation_mutagenesis}.rst`, the design package (`screen`, `beam_substitution`,
`greedy_substitution`, `greedy_marginalize`), `ersatz.py`, `results.py`, and
`saturation_mutagenesis.py`. Searched all repository source/docs for `library`,
`return`, `torch.Tensor`, `tensor`, `array`, `batch`, `variants`, `mutagen`,
`design`, `class .*Library`, `class .*Dataset`, `substitute`, `randomize`,
`insert`, and `delete`, then enumerated every class/NamedTuple/dataclass definition.

## Row-level finding

The row **does discriminate** on the definition as written: one `yes`, five
`partial`, and two `no`. Its strongest separation is architectural. PoolParty alone
returns a dedicated whole-library carrier; VaLiAnT exposes only files, and DNA
Chisel's purpose-defined object stops at one sequence. The broad `partial` bucket
then deliberately groups five quite different general-purpose carriers: an HTTP
response body (MPRAnator), R list/tibbles (MPRA Design Tools), pandas DataFrames
(Oligopool Calculator), decoded JSON/dictionaries (Mutation Maker), and PyTorch
tensors (tangermeme). Therefore the row has limited discrimination *within* tools
that return complete results in memory, but it remains consistently applicable;
no replacement row was used or proposed.
