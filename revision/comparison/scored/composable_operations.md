# Audit: composable operations

Assessed 2026-08-15. Primary sources only. The records under
`revision/tool_survey/final/` were used only as navigation leads; no statement or
rating from those records is evidence in this audit.

## Measurement instrument and operational test

The binding definition is `revision/tables/ROW_DEFINITIONS.md`, row 2:

> Two or more distinct **library-design** operations share a carrier type and
> compose in either order, with fan-out (one intermediate feeding two downstream
> operations) expressible.

The definition excludes general sequence utilities, cloning simulation, and
internal caching; it assigns `partial` to fixed-order, small-subset, or
output-artifact-only composition, and `no` when source fixes the order or there
are no design operations.

**Concrete operational test (stated before any tool was examined):** for each
tool, determine whether its documented public design interface exposes at least
two distinct library-design operations that accept and return the same carrier
type, can be applied as A→B and B→A, and allow one saved intermediate to feed
both operations; score `partial` only for fixed-order/subset/artifact-refeeding
composition and `no` when source fixes the design order or exposes no design
operations.

## Results

| Tool | Score | Confidence | Short basis |
|---|---|---:|---|
| PoolParty | **yes** | high | `Pool -> Pool` design operations work in both orders; a saved `Pool` fans out in one DAG |
| VaLiAnT | **no** | high | source fixes preprocessing, PPE, mutator, assembly, and output order |
| MPRAnator | **partial** | high | authors document PWM-output → Motifs input, but only as one-way FASTA re-feeding among separate web tools |
| MPRADesignTools | **no** | high | one exported design pipeline; its component order is hard-coded inside `processVCF`/`processSnp` |
| Oligopool Calculator | **yes** | medium | documented in-memory DataFrame flow supports reorderable `primer` and `motif`; DAG example proves fan-out/rejoin |
| Mutation Maker | **no** | high | SSM, QCLM, and PAS are separate endpoints/tasks with incompatible job types and no hand-off |
| DNA Chisel | **partial** | high | two in-place design phases share a problem object but are constrained to resolve-then-optimize |
| tangermeme | **partial** | high | two design algorithms share tensors and can reverse only for singleton outputs; batch/library output is rejected downstream |

## Per-tool audit

### PoolParty — **yes** (confidence: high)

**Positive evidence.** Both tested operations are design operations, not excluded
utilities. `base_ops/mutagenize.py:18-32,67-70` declares
`mutagenize(pool: Union[Pool, str], ...) -> Pool` and documents:

> A Pool that generates mutated sequences.

`scan_ops/deletion_scan.py:10-23,55-59` likewise declares
`deletion_scan(pool: Union[Pool, str], ...) -> Pool` and documents a pool yielding
the designed deletions. Their method forms preserve the carrier:
`pool_mixins/common_ops_mixin.py:25-37` and
`pool_mixins/scan_ops_mixin.py:99-111` both return `Self`.

I ran the operational test read-only with the repository venv and
`PYTHONDONTWRITEBYTECODE=1`:

```python
base = pp.from_seq("ACGTAC")
ab = base.mutagenize(num_mutations=1, mode="sequential", num_states=2).deletion_scan(
    deletion_length=1, deletion_marker=None, mode="sequential", num_states=2)
ba = base.deletion_scan(
    deletion_length=1, deletion_marker=None, mode="sequential", num_states=2
).mutagenize(num_mutations=1, mode="sequential", num_states=2)

shared = base.mutagenize(num_mutations=1, mode="sequential", num_states=2)
fan = pp.stack([
    shared.deletion_scan(deletion_length=1, deletion_marker=None,
                         mode="sequential", num_states=2),
    shared.insertion_scan("TT", positions=[1], mode="sequential")
])
```

Observed: `ab` and `ba` were both `DnaPool`, each generated a `(4, 2)` DataFrame;
`fan` was a `DnaPool` with 6 states and generated `(6, 2)`. The shared node is
expressible as one object, not as duplicated output: `_topo_sort_operations`
deduplicates operations through `visited` (`generate_library.py:220-237`), while
each row's downstream operations retrieve their parents from `seq_cache`
(`:260,270-292`). Caching is not counted as capability evidence; it only confirms
that the syntax denotes a shared node.

**Disconfirmation and different-value condition.** I searched
`src/poolparty/{base_ops,scan_ops,pool_mixins,region_ops,multiscan_ops,state_ops}`
and `docs/operations/` for `mutagenize`, `deletion_scan`, `insertion_scan`,
`Callable[[Pool], Pool]`, `_topo_sort_operations`, `visited`, and `seq_cache`, and
ran A→B, B→A, and fan-out. `partial` would have resulted if one ordering failed,
only serialized output could be re-fed, or only an excluded utility supplied the
second operation; `no` would have resulted if operation order were hard-coded.
None occurred.

### VaLiAnT — **no** (confidence: high)

**Evidence of a source-fixed design order.** `src/valiant/mutator_type.py:22-32`
defines the closed seven-value `MutatorType` enum. The only declared mutator
dependency is one-way (`:58-61`):

```python
# Mutator types based on the output of other mutators
DEPENDENT_MUTATOR_TYPES = {MutatorType.SNV_RE: {MutatorType.SNV}}
```

The public job does not expose those as carrier-preserving operations.
`sge_proc.py:306-337` applies PAM-protection edits and then invokes
`Targeton.process`; `targeton.py:300-312` fixes the latter order as:

```python
self._process_custom_variants(...)
self._process_pattern_variants(...)
```

Only after that does `sge_proc.py:341-368` build the `MetaTable` and write CSV.
At the outer level `sge_proc.py:409-503` fixes loading, annotation/background
lifting, targeton processing, and output traversal. Packaging exposes one console
script (`pyproject.toml:23-24`), not operation entry points.

**Search and disconfirmation.** I enumerated the full `develop` tree, read
`mutator_type.py`, `targeton.py`, `sge_proc.py`, and `pyproject.toml`, and
stream-grepped all `.py/.md/.rst/.toml` files for
`compos|chain|pipeline|workflow|plugin|entry.?point|fan.?out|DAG|dependenc`.
The only capability-relevant dependency match was the internal one-way
`SNV_RE <- SNV`; other `chain` hits were `itertools.chain`. I also checked for a
plugin/entry-point interface and found none. `partial` would require at least a
documented artifact hand-off between design operations; `yes` would require two
public operations with both orders plus fan-out. I searched for both and found
neither. This is `no` specifically because source fixes the design order, not
because VaLiAnT lacks design functionality.

### MPRAnator — **partial** (confidence: high)

**Positive but bounded evidence.** The authors' official documentation states:

> The output of PWM Seq-Gen can be inputted into MPRAs Motifs tool, therefore
> allowing for the design of MPRA experiments using PWMs.

Source: <https://www.genomegeek.com/MPRA/documentation/>, section `PWM-SeqGen`.
That page also says PWM Seq-Gen results are generated motif sequences in FASTA,
and the MPRA result is displayed/downloaded as FASTA. Thus two distinct design
operations genuinely chain, but by downloading/copying an output artifact and
submitting it to another web form, in one documented direction.

The paper independently identifies the distinct operations: PWM Seq-Gen performs
probabilistic or thresholded PWM realization, while Transmutation designs
negative controls (local PDF `Georgakopoulos-Soares2017_MPRAnator.pdf`, extracted
with PyMuPDF; text lines matching `PWM Seq-Gen` and `Transmutation`). The repository
confirms separation: `iliasApp/urls.py:7-31` maps independent Motif, SNP,
Transmutation, documentation, and API routes, and `views.py` has independent
result handlers that render/download complete response artifacts.

**Search and disconfirmation.** I searched the official documentation for
`output of PWM`, `inputted`, `Transmutation`, and `negative controls`; extracted
the full paper with PyMuPDF and searched the same concepts; enumerated the
50-entry repository tree; and read `iliasApp/urls.py`, `views.py`, and
`templates/iliasApp/docs.html`. No pipeline object, common in-memory carrier,
reverse Motifs→PWM order, or fan-out/rejoin API appears. `yes` would require those
features. `no` would require the documented PWM→Motifs hand-off to be absent;
I explicitly searched for it and found the quoted author documentation. The
artifact-only rule therefore makes this `partial`.

### MPRADesignTools / designMPRA — **no** (confidence: high)

**Evidence of one hard-coded design pipeline.** The package namespace has only
two exports (`NAMESPACE:3-4`):

```r
export(processVCF)
export(spread_and_fix_indels)
```

The second is excluded input preprocessing, not a library-design operation. Its
own documentation says it reads/spreads/fixes VCF rows and writes `*_fixed.vcf`
"ready to be fed into processVCF()" (`R/processVCFfast.R:1395-1408`); this is parse/
format normalization, not design.

The sole exported design operation is the monolith `processVCF`.
`R/processVCFfast.R:1099-1115` declares it; `:1121-1202` fixes VCF loading and
barcode-set filtering, `:1212-1237` allocates barcodes and calls the internal
`processSnp`, and `:1241-1297` materializes/writes the result. Inside
`processSnp` (`:210-993`) the SNV/insertion/deletion branches, context creation,
primer/enzyme/barcode concatenation, site repair, and return order are source
code, not independently callable documented operations. For example `:940-961`
hard-codes the concatenation order with `paste0(fwprimer, ..., enzyme1, enzyme2,
barcodes, ..., revprimer)`.

**Search and disconfirmation.** I enumerated both official repositories. The
package tree contains `NAMESPACE`, three `R/` files, README, data, and man pages;
the companion `designMPRA` tree contains Shiny `ui.R/server.R`, scripts, and the
supplement. I stream-grepped all `.R/.Rd/.Rmd/.md` files in both repositories for
`compos|chain|pipeline|workflow|plugin|fan.?out|DAG|dependenc`. The only relevant
companion hit was the paper's static "MPRA Design Tools Workflow" figure; there
was no operation interface. I also read both exports and the complete
`processVCF`/`processSnp` path. `partial` would require two genuine design
operations even in fixed order or via artifact re-feeding; `yes` would require
both orders and fan-out. I searched for each. Because preprocessing is excluded
and the design components are fixed inside one source pipeline, the row's
explicit `no` rule applies.

### Oligopool Calculator — **yes** (confidence: medium)

**Positive evidence.** The official user guide explicitly defines a common
in-memory carrier (`docs/docs.md:136-158`):

> Most modules return both a DataFrame and a stats dictionary

and

> Chainable: Output of one module feeds into the next

Two distinct qualifying design operations have the same documented contract.
`docs/api.md:143-200` documents `primer(input_data: str | pd.DataFrame, ...)` and
returns `(DataFrame, stats_dict)`; `:230-279` documents
`motif(input_data: str | pd.DataFrame, ...)` and the same return, with the purpose
"Insert sequence motifs (per-variant or constant anchors) with constraint
satisfaction." Both can use an original sequence column as context and add a
differently named column. Therefore a saved `df` permits both
`motif(primer(df))` and `primer(motif(df))`, and can independently feed both
calls. This is an inference directly from the documented input/return contracts;
the package was not installed or executed.

Fan-out is also a shipped, documented design feature. The runnable
`examples/cli-yaml-pipeline/mpra_design_parallel.yaml:5-25` says:

> two independent barcode branches are designed from the same backbone, then
> recombined with `join` into a single design table for downstream steps.

It encodes `primer -> {barcode_a, barcode_b} -> join -> spacer -> final`.
`docs/docs.md:1635-1677` defines serial and explicit dependency-DAG formats, and
`docs/api.md:526-570` documents `join` for reconciling parallel design branches.
The example's branches use CSV artifacts, but that is not the only carrier: the
documented Python API accepts/returns the DataFrame itself.

**Restriction search and different-value condition.** I read `docs/docs.md`
sections "The DataFrame Flow", `join`, and "Parallel Pipeline Execution";
`docs/api.md` design modules `barcode`, `primer`, `motif`, `spacer`, and `join`;
the example README and both serial/parallel design YAML files. I searched for
`Design order`, `Chainable`, `after`, `DAG`, `join`, and all design-operation
signatures. The docs do impose order for some pairs: barcode is "after
primers/motifs, before spacers" (`api.md:119-124`) and spacer is last
(`:353-358`). Those pairs would only be `partial`; they are not the tested pair.
No ordering prohibition exists between `primer` and `motif`, and both accept the
same original context. `partial` would have resulted if every pair were fixed in
order or fan-out were artifact-only; `no` if there were only a fixed monolith.
The medium rather than high confidence records the lack of a local behavioral
run, not uncertainty about the documented contracts.

### Mutation Maker — **no** (confidence: high)

**Evidence of isolated workflows.** The legacy API exposes separate endpoints
for `/ssm`, `/qclm`, and `/pas`, each dispatching a different task
(`api/api.py:40-61`). The newer FastAPI port preserves exactly those isolated
endpoints (`api/server_fastapi.py:53-65`). The task implementations instantiate
different input types and immediately invoke independent solvers
(`backend/tasks.py:42-60,68-73`):

```python
input = SSMInput(data); return ssm_solve(...)
input = QCLMInput(data); return qclm_solve(input)
input = PASInput(data); output = pas_solve(input)
```

The common response is a Celery task URL/result artifact, not a design carrier
accepted by another workflow (`api/api.py:89-105,123-129`). No endpoint accepts
another workflow's output.

**Search and disconfirmation.** I enumerated the current `master` tree and read
`AGENTS.md`, root/backend READMEs, both API implementations, `backend/tasks.py`,
and the SSM/QCLM/PAS solver and type file locations. I stream-grepped all
`.py/.md/.rst/.toml/.ts/.tsx` files for
`compos|chain|pipeline|workflow|plugin|entry.?point|fan.?out|DAG|dependenc`.
The `compose` matches were Docker Compose or React's `recompose`; workflow
matches describe the three separate screens/solvers. I specifically checked the
new `server_fastapi.py` so the conclusion does not rely on the legacy Hug API.
`partial` would require a documented SSM→QCLM/PAS (or other design-operation)
hand-off; `yes` would require a shared carrier, both orders, and fan-out. Neither
exists. This is `no` because the public design workflows are isolated, not
because the tool lacks design operations.

### DNA Chisel — **partial** (confidence: high)

**Positive but fixed-order evidence.** The public `DnaOptimizationProblem`
carrier supports two distinct in-place design phases. The official README gives
the canonical sequence (`README.rst:69-82`):

```python
problem.resolve_constraints()
problem.optimize()
final_sequence = problem.sequence
```

These phases both modify the same problem/sequence, so they chain. They do not
compose in either order for a nontrivial problem.
`ObjectivesMaximizerMixin.py:26-32,62-68` explicitly rejects optimization before
constraints pass:

> Optimization can only be done when all constraints are verified.

The combined public operation hard-codes the same order:
`DnaOptimizationProblem.py:202-253` documents `optimize_with_report` as "Resolve
constraints, optimize objectives" and calls `resolve_constraints()` before
`optimize()` at `:232-252`. The methods mutate one object; no documented fan-out
or branch/rejoin carrier exists.

**Search and disconfirmation.** I read the README, both solver mixins,
`DnaOptimizationProblem.py`, `SpecificationSet.py`, core-class/CLI/example doc
paths, and enumerated the full `master` tree. I stream-grepped every
`.py/.md/.rst` file for the alternatives `fan.?out`, `DAG`, `pipeline`, `chain`,
`compose`, `resolve_constraints()`, and `optimize()`. Numerous examples repeat
resolve-then-optimize; none documents the
reverse or fan-out/rejoin. `SpecificationSet.py:1-13` composes constraint
specifications, but specifications are not sequence-producing design operations,
so that cannot produce `yes`. `yes` would require a second order plus fan-out;
`no` would require the two modifying phases not to chain at all. The fixed-order
clause makes the observed capability `partial`.

### tangermeme — **partial** (confidence: high)

**Positive evidence.** `greedy_substitution` and `beam_substitution` are distinct
model-guided design algorithms and share a `torch.Tensor` carrier.
`tangermeme/design/greedy_substitution.py:24-48` has
`X: torch.Tensor -> torch.Tensor` and describes itself as a design function;
`:146-149,249` documents and returns the edited tensor.
`design/beam_substitution.py:25-67` also has
`X: torch.Tensor -> torch.Tensor` and describes beam search through edit space;
`:173-176,313-314` returns the best designed tensor(s). With
`beam_substitution(..., n_best=1)`, each result has shape `(1, alphabet, length)`,
so `greedy(beam(X))` and `beam(greedy(X))` are expressible, and a saved singleton
tensor can feed both functions.

**Why only partial.** Both functions explicitly reject a carrier batch larger
than one. `greedy_substitution.py:155-159` says it "designs a single sequence at
a time and requires X.shape[0] == 1"; `beam_substitution.py:182-186` says the
same. Beam can return a genuine multi-sequence batch when `n_best > 1`
(`:131-135,173-176`), but that output cannot feed either design operation. Thus
the reverse ordering works only for the documented singleton subset of the
carrier, exactly the row's small-subset `partial` case. The restriction is not
hidden: the official changelog states both now raise `ValueError` for batch size
other than one.

**Search and disconfirmation.** I read all five files in `tangermeme/design/`,
the design API page/changelog and bundled design reference, and stream-grepped
those sources/docs for `greedy_substitution`, `beam_substitution`, `compose`,
`chain`, `pipeline`, `fan-out`, `stack`, `return Tensor`, and `single sequence`.
I also checked the package
was not installed locally, so no behavioral run was claimed. `yes` would require
the batch output of one documented design operation to remain acceptable to the
other (or another unrestricted pair); `no` would require no same-carrier design
pair. I explicitly searched all public design algorithms: `screen` takes a shape,
not a prior `X`; `greedy_marginalize` documents a rank-changing return; only the
singleton greedy/beam pair passes. Hence `partial`.

## Row-level finding

The row **does discriminate** on one consistent scale: 2 `yes`, 3 `partial`, and
3 `no`. It also exposes an important interpretation boundary. A `yes` need not be
a graph object: PoolParty's `Pool` DAG and Oligopool Calculator's ordinary
DataFrame flow both pass because the instrument measures expressibility, not
carrier richness or caching. Conversely, DNA Chisel's declarative specification
composition does not pass as full composable operations because the actual
sequence-modifying phases are fixed-order, and tangermeme's tensor compatibility
is only singleton-wide. Those distinctions follow the supplied definition; no
replacement row was introduced.

## Search log (audit trail)

All repository access was read-only; nothing was installed. Remote sources were
read from official GitHub repositories through raw files, GitHub tree/metadata
endpoints, or streamed `codeload` archives (archives were not written to disk).

- Read the complete binding definition with
  `sed -n '1,240p' revision/tables/ROW_DEFINITIONS.md` before looking at tools.
- Inventoried the workspace with `find` to locate the supplied papers, local
  PoolParty repository, and lead files; no other tool repository was present
  locally.
- Inspected `final/{poolparty,valiant,mpranator,mpradesign,oligopoolcalc,
  mutationmaker,dnachisel,tangermeme}.md` only to locate primary files/pages;
  no claims or ratings were reused.
- Initial discovery web queries were `site:github.com/ayaanhossain/oligopool
  "Parallel Pipeline Execution"`, `site:github.com/Edinburgh-Genome-Foundry/
  DnaChisel "resolve_constraints" "optimize"`, `site:github.com/jmschrei/
  tangermeme "greedy_substitution" "beam_substitution"`, and
  `site:github.com/cancerit/VaLiAnT "DEPENDENT_MUTATOR_TYPES"`; direct official
  repository files, not search-result summaries, supplied the evidence.
- PoolParty: `rg` over the source/doc paths and symbols listed in its section;
  line-numbered reads of the three operations, mixins, and generator; one
  read-only venv execution of A→B, B→A, and fan-out.
- VaLiAnT: GitHub metadata/default branch and recursive tree; raw reads of
  `mutator_type.py`, `targeton.py`, `sge_proc.py`, `pyproject.toml`; full text-file
  archive grep with the exact regex reported above.
- MPRAnator: full PDF extraction with
  `python3 -c "import fitz; ..."` and grep for `PWM Seq-Gen|Transmutation|negative
  controls|inputted`; GitHub recursive tree; raw route/view/docs-template reads;
  web searches for the exact official sentence `The output of PWM Seq-Gen can be
  inputted` and the official GenomeGeek documentation page.
- MPRADesignTools: metadata/default branch and recursive trees for both
  `andrewGhazi/mpradesigntools` and `andrewGhazi/designMPRA`; raw reads of
  `NAMESPACE` and line ranges covering both exports plus the complete
  `processSnp/processVCF` route; full text archive greps for both repositories.
- Oligopool Calculator: metadata/default branch and recursive tree; raw reads of
  `docs/docs.md`, `docs/api.md`, both example README/YAML paths; searches for
  `DataFrame Flow|Chainable|Design order|Parallel Pipeline Execution|after|join`;
  local import probe returned `None` and no install was attempted.
- Mutation Maker: metadata/default branch and recursive tree; raw reads of both
  APIs and `backend/tasks.py`; full source/docs/frontend text archive grep with
  the exact regex reported above.
- DNA Chisel: metadata/default branch and recursive tree; raw reads of README,
  both solver mixins, core problem, and `SpecificationSet`; full source/docs
  archive grep with the exact regex reported above.
- tangermeme: metadata/default branch; raw reads of all public design algorithm
  files and design docs/changelog; full design-source/docs archive grep with the
  exact regex reported above; local import probe returned `None` and no install
  was attempted.
