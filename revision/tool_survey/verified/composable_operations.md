# Verified audit — row `composable_operations`

**Auditor pass:** single-rater, all 13 tools, one threshold.
**Date:** 2026-08-10
**Row question:** Do design steps compose/chain/nest arbitrarily in a user-chosen order, or is the pipeline fixed by the tool?

**Row verdict: `reword`** (cells are still filled and internally consistent — see §4 for why the row as
worded is unsafe to print, and §5 for the concrete replacement wording).

---

## 1. Operational test (fixed BEFORE any tool was opened)

Written to `verified/OPTEST_composable_operations.txt` before the first source file was read.

> For each tool, enumerate its **library-design operations** — operations whose output is candidate
> library/variant sequence content (mutagenesis, recoding, insertion/replacement of a user-supplied
> payload, barcoding/indexing, combinatorial or degenerate expansion, constraint-directed sequence
> optimisation, attribution-guided sequence editing). General sequence utilities (parse / revcomp /
> translate / pad / one-hot / write) and cloning-simulation primitives do **not** count.
>
> Then, from primary sources only (code `file:line` or exact quoted official docs):
>
> **(a) COMPOSITION** — is there a carrier type `T` such that at least **two distinct** design
> operations accept `T` and return `T` (or return something another design op accepts), so that a user
> can write `op_b(op_a(x))` **and** `op_a(op_b(x))` — i.e. the order is the *user's* choice, not the
> tool's?
>
> **(b) BRANCHING** — can one intermediate result be bound once and fed to **two different**
> downstream design operations (fan-out), as shown in source, tests, or official examples?
>
> - **yes** = (a) demonstrated with ≥2 user-orderable design ops on a shared carrier **and** (b) expressible.
> - **partial** = design ops exist and some user-expressible chaining is possible, but **either** the order
>   is set by a fixed driver/pipeline with configurable stages, **or** only a small documented subset
>   chains, **or** composition works only by manually re-feeding output files between invocations.
> - **no** = the topology is fixed in the tool's source (the user only toggles hard-coded steps);
>   inserting or reordering a design step requires editing the source; **or** there is only one design
>   op, so there is nothing to compose.
> - **unknown** = source genuinely inaccessible.

### 1.1 Three interpretation rulings, fixed once and applied to all 13

These were needed to keep the anchors (PoolParty=yes, VaLiAnT=no) satisfied while staying re-runnable.

**R1 — "some chaining is possible" means *expressible by the user*.** A tool whose internal stages chain
(e.g. VaLiAnT applies PAM/background variants to the reference, then runs mutators on that background)
does **not** earn `partial` for that, because the chain lives inside a single hard-coded driver and the
user cannot insert, remove, or reorder a link. This is what keeps VaLiAnT at the `no` anchor.

**R2 — a design op is named for a design act, not inherited from a container API.** `tangermeme.ersatz.insert`
("Insert a motif into a set of sequences at a defined position") counts. `Bio.Seq.replace` ("Return a copy
with all occurrences of subsequence old replaced by new" — a mirror of `str.replace`) does not. This is the
line that stops the row from collapsing into "is this a pure-functional array library?".

**R3 — "manual re-feeding" requires the tool's own design *output* to be accepted by the tool's own
design *input*, in a way that deepens the library.** MPRAnator emits FASTA and consumes FASTA, so round 2
inserts a second motif set into round-1 oligos → counts. VaLiAnT emits `*_meta.csv` / `*_unique.csv` /
a variants VCF and consumes an indexed reference FASTA + targeton TSV; re-feeding its output VCF as
`custom variants` merely re-declares the same variants against the same reference and adds no design
layer → does not count. (Checked: `README.md:96–101` output list; `src/valiant/meta_table.py:679–686`
writes `oligo_name,mseq` CSV.)

**Deliberately NOT part of the threshold:** whether a shared node is *computed once* (memoised). PoolParty
memoises; requiring that would be a self-serving criterion, since it is an efficiency property, not a
composition property. PoolParty scores `yes` on expressibility alone.

---

## 2. Cells, with primary-source evidence

### 2.1 `oligopoolcalc` = **yes** (prior: yes — confirmed, and it is a strong yes)

Design-mode ops all share carrier `T = pandas.DataFrame` of an annotated oligo pool.

`oligopool/barcode.py:17-34`:
```python
def barcode(
    input_data:str|pd.DataFrame,
    ...
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
```
`oligopool/join.py:13-22`:
```python
def join(
    input_data:str|pd.DataFrame,
    other_data:str|pd.DataFrame,
    ...
    Useful for recombining parallel branch outputs into a single design table.
```
`oligopool/__init__.py` package docstring: *"Modules return (DataFrame, stats). Chain them iteratively; use
patch_mode=True to extend pools without overwriting existing designs."*

`docs/docs.md:157`: *"- **Chainable**: Output of one module feeds into the next"*
`docs/docs.md:1711`: *"If two design steps are truly independent … you can run them as parallel branches and
then recombine their outputs with `join`."*
`docs/docs.md:2475` (section **"Composability Cheat Sheet"**): *"Compound element | chain
`motif`/`barcode`/`primer`/`spacer` → `merge`/`revcomp` → repeat"*

Fan-out, verbatim from `examples/cli-yaml-pipeline/mpra_design_parallel.yaml:12-22`:
```yaml
    - name: primer
      command: primer
    - name: barcode_a
      command: barcode
      after: [primer]
    - name: barcode_b
      command: barcode
      after: [primer]
    - name: rejoin
      command: join
      after: [barcode_a, barcode_b]
```
One node (`primer`) feeding two downstream design branches, then a fan-in. (a) and (b) both satisfied.

**Disconfirmation:** a `no`/`partial` would have required either that the design modules only accept file
paths and not in-memory frames (they accept `str|pd.DataFrame`), or a documented mandatory order. I grepped
`docs/docs.md` for `any order|order of|arbitrary order|order matters|in order` — **zero hits**, i.e. no
mandated order. The only order constraints found are advisory (`docs/docs.md:2395,2497`: treat `compress`
output as a *terminal* branch), which is a documented dead-end for one module, not a fixed pipeline.

### 2.2 `pydna` = **no** (prior: **yes** → CHANGED)

pydna's composable carrier (`Dseqrecord`/`Dseq`) is real and excellent, but nothing composed on it is a
library-design operation. Full public surface, `src/pydna/all.py:45-59`:
```python
from pydna.amplify import Anneal
from pydna.amplify import pcr
from pydna.assembly import Assembly
from pydna.genbank import genbank
from pydna.genbank import Genbank
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.readers import read
from pydna.readers import read_primer
from pydna.parsers import parse
from pydna.parsers import parse_primers
from pydna.design import primer_design
from pydna.design import assembly_fragments
from pydna.utils import eq
from pydna.genbankfixer import gbtext_clean
```
Every documented chain is cloning/assembly simulation, e.g. `README.md:161-205`:
```python
a, payload, c = pcr_product.cut (BamHI)
...
rec_vector= (linear_vector_bgl + payload).looped().synced(vector)
```
This is exactly the case the row's CRITICAL DISTINCTION excludes.

**Disconfirmation (searched, absent):** `grep -rniE "mutagen|mutation|variant|library|degenerate|saturat|randomiz" src/` over all 39
modules returns only (i) `snapgene_history_parser.py:345` `elif node.operation == "primerDirectedMutagenesis":`
— *parsing* a SnapGene history log, not generating a library; (ii) `recombinase.py:23,72` and
`sequence_regex.py:17` — IUPAC degenerate bases used for *searching*; (iii) `primer_screen.py:144
expand_iupac_to_dna` — a primer-screening helper. `grep -rniE "\blibrar(y|ies)\b"` over `.py/.rst/.md/.ipynb`
returns only "full library documentation" boilerplate plus `docs/example_gallery.rst:24,33`, which link to
**third-party** notebooks (`cazyme_primer_design`, `teemi`) for library design — i.e. pydna itself supplies
no library-design op. Module list checked individually: `design.py` (primer design + assembly linkers),
`codon.py`, `primer_screen.py`, `crispr.py`, `gateway.py`, `recombinase.py`, `fusionpcr.py`, `cre_lox.py`.
None generates a variant library.

⚠️ **This is the cell that most needs the row reworded.** `no` is defined in the rubric as "the tool fixes
the pipeline order", which is flatly false for pydna — pydna fixes no pipeline at all. See §4.

### 2.3 `tangermeme` = **yes** (prior: **partial** → CHANGED)

Carrier `T = torch.Tensor` one-hot, shape `(N, 4, L)`. Four ersatz ops are pure `T → T`:
`tangermeme/ersatz.py:20-25` `def insert(X, motif, start=None, alphabet=...) -> torch.Tensor` with
`:94 return torch.cat([X[:, :, :start], motif, X[:, :, start:]], dim=-1)`;
`:97 def substitute(...) -> torch.Tensor` with `:186 X = torch.clone(X)` … `:195 return X`;
`:198 def multisubstitute(...)` → `:297 return X`; `:300 def delete(X, start, end) -> torch.Tensor` →
`:340 return torch.cat([X[:, :, :start], X[:, :, end:]], dim=-1)`.
Purity is documented: `ersatz.py:31-32` *"It will then return a copy of the data with the insertion,
leaving the original data unperturbed."*
`README.md:86`: *"tangermeme implements **atomic sequence operations** to help you ask 'what if?' questions
of your data."*
Design subpackage also returns the carrier: `tangermeme/design/greedy_substitution.py:24-40`
`… ) -> torch.Tensor` with `:241 return X`; `design/screen.py:185 return X`;
`design/beam_substitution.py:306 return torch.cat([seq for _, seq in beam[:n_best]], dim=0)`.

**(a) documented, two DISTINCT design ops, output→input** — `docs/tutorials/Tutorial_B6_Design.ipynb`
cell 31 then cell 33:
```python
from tangermeme.ersatz import substitute
X = random_one_hot((1, 4, 2000), random_state=0).type(torch.float32)
X = substitute(X, "GTGACTCATC")
```
```python
X_hat3 = greedy_substitution(model, X, y, motifs, max_iter=1, verbose=True)
```
**(b) documented fan-out** — same notebook, cell 43: one `X` feeding two different downstream design ops:
```python
X_hat_greedy = greedy_substitution(model, X, y, motifs, output_mask=idxs, max_iter=5)
X_hat_beam = beam_substitution(model, X, y, motifs, output_mask=idxs, beam_size=5, max_iter=5, verbose=True)
```
The library's own higher-order functions compose ersatz ops (`marginalize.py:15,126`;
`space.py:16,151`; `ablate.py:14`), confirming this is the intended pattern, not an accident.

**Disconfirmation:** a `partial` would have required a fixed driver, a mandated order, or type
incompatibility. There is no driver — `tangermeme/` is 21 flat modules with no pipeline object
(`grep -rn "class .*Pipeline\|def run_pipeline"` → nothing). One genuine wrinkle found and weighed:
`randomize` (`:422 return torch.stack(X_rands).permute(1, 0, 2, 3)`) and `shuffle` (`:498`) return a
**4-D** tensor, so `substitute(shuffle(X))` needs a reshape; that limits two of seven ersatz ops, not the
four 3-D ops that carry the `yes`. I also searched for nested calls anywhere in the repo
(`grep -rnE "(insert|substitute|delete|randomize|shuffle)\((insert|substitute|delete|randomize|shuffle)\("`
over `.py/.ipynb/.rst`) — **zero hits**; the composition is by named intermediate (cell 31→33), which
satisfies (a) equally.

### 2.4 `poolparty` = **yes** (prior: yes — confirmed; same standard as competitors)

51 public methods across the 10 `pool_mixins` modules, of which **47 are typed `-> Self`**, i.e. `Pool → Pool`
(`grep -rn "^    def [a-z]" src/poolparty/pool_mixins/*.py | grep -v "def _"` → 51;
`grep -rn ") -> Self:" src/poolparty/pool_mixins/*.py` → 47). The four non-`Self` ones are the terminals
`score`, `materialize`, `to_df`, `to_file`. Full `-> Self` list: `mutagenize, shuffle_seq, recombine, filter,
insert_from_iupac, insert_from_motif, insert_kmers, flip, rc, annotate_orf, stylize_orf, mutagenize_orf,
translate, slice_seq, add_prefix, swapcase, upper, lower, clear_gaps, clear_annotation, stylize,
reverse_translate, filter_gc, filter_homopolymer, filter_complexity, filter_dust, filter_restriction_sites,
annotate_region, apply_at_region, extract_region, insert_tags, remove_tags, replace_region, clear_tags,
mutagenize_scan, deletion_scan, insertion_scan, replacement_scan, shuffle_scan, subseq_scan,
deletion_multiscan, insertion_multiscan, replacement_multiscan, repeat, sample, shuffle_states,
slice_states` (signatures at `pool_mixins/common_ops_mixin.py:25-37`, `scan_ops_mixin.py:24-36`,
`region_ops_mixin.py:256`).

Documented: `docs/pool.rst:12-14` — *"Pools are also *immutable*: every operation returns a new Pool, leaving
the original unchanged. You can branch a pipeline at any point and apply different operations to each branch
without interference."* `docs/pool.rst:501-503`:
```python
base = pp.from_iupac("NNNN", mode="sequential")
branch_a = base.mutagenize(num_mutations=1).named("branch_a")
branch_b = base.copy().mutagenize(num_mutations=2).named("branch_b")
```
`README.md:53-62` shows fan-out then fan-in (`template` → `mutagenize` / `deletion_scan(...).repeat(...)`
→ `pp.stack([mut_pool, del_pool])`); `state_ops/stack.py:19` is the fan-in combinator.

**Behavioural verification** (read-only venv, `PYTHONDONTWRITEBYTECODE=1`):
```
t = pp.from_seq("ACGTACGTACGT")
t.rc().mutagenize(num_mutations=1, mode="sequential")   -> DnaPool, 36 states
t.mutagenize(num_mutations=1, mode="sequential").rc()   -> DnaPool, 36 states   # both orders legal
shared = t.mutagenize(...); br1 = shared.rc(); br2 = shared.repeat(times=2)
shared: 36  br1: 36  br2: 72  pp.stack([br1, br2]): 108                        # fan-out + fan-in
```
Both orders are legal (a) and one node feeds two branches (b).

**Disconfirmation:** a `partial` would have required a driver that fixed the order. `generate_library.py:105,220`
`sorted_ops = _topo_sort_operations(pool)` sorts whatever DAG the *user* built; it does not impose an order.
I specifically tried the harsher criteria and rejected them as self-serving (see §1.1, "Deliberately NOT part
of the threshold"). PoolParty would still score `yes` if memoisation were removed.

### 2.5 `valiant` = **no** (prior: no — confirmed; anchor holds)

`src/valiant/mutator_type.py:22-32` — closed 7-member enum:
```python
class MutatorType(str, Enum):
    DEL = 'del'; SNV = 'snv'; SNV_RE = 'snvre'; IN_FRAME = 'inframe'
    ALA = 'ala'; STOP = 'stop'; AA = 'aa'
```
`src/valiant/mutator.py:43-49` is a closed dispatch dict. All mutators are applied to the **same** input
region, in parallel, in one hard-coded call: `src/valiant/cdna_proc.py:94-96`
```python
    # Generate pattern variants
    vars, annot_vars = targeton_cfg.mutator_collection.get_variants(
        codon_table, r_seq)
```
No mutator's output is any other mutator's input. Step order lives in `proc_targeton`
(`cdna_proc.py:62`, `sge_proc.py:267`).

**Disconfirmation, three probes:**
(i) plugin extension — `pyproject.toml:23-24` has only `[project.scripts] valiant = "valiant.__main__:main"`;
no plugin entry-point group. (ii) re-feeding — outputs are CSV/VCF/JSON (`README.md:96-101`;
`meta_table.py:679-686` writes `oligo_name,mseq`), inputs are an indexed reference FASTA + targeton TSV +
VCF; the designed oligos are never an accepted input, and re-feeding the output VCF as `custom variants`
re-declares the same variants against the same reference without adding a design layer (ruling **R3**).
(iii) documented multi-pass — `grep -n -iE "two-pass|re-run|rerun|pipeline|chain|compose" README.md` →
**zero hits**. VaLiAnT *does* chain internally (background/PAM variants → mutators on the background
sequence); ruling **R1** excludes that, which is what preserves this anchor.

### 2.6 `ledidi` = **partial** (prior: **no** → CHANGED)

Two design ops on carrier `T = torch.Tensor (N,4,L)`, and the chain is in the official tutorial.
`ledidi/ledidi.py:78` `def ledidi(model, X, y_bar, ...)` → `:262 return ledidi_output[0] if len(...)==1 else ...`
where `ledidi_output = [X_bar]` (`:254`).
`ledidi/pruning.py:25` `def greedy_pruning(model, X, X_hat, threshold=1, target=None, verbose=False)` →
`:129 return X_hat`. Docstring `pruning.py:33`: *"the procedure will stop and return the remaining edits."*

`docs/tutorials/Tutorial_7_-_Validating_Your_Designs.ipynb`, cells 9 → 10:
```python
X_bar = ledidi(model, X, y_bar, n_samples=50, verbose=True)
```
```python
X_bar_p = torch.cat([greedy_pruning(model, X, xb[None], threshold=0.25, verbose=False) for xb in tqdm(X_bar)])
```
Output of design op 1 is the input of design op 2 — genuine user-expressible chaining, which is more than
`no` allows.

**Why not `yes`:** the order is semantically forced. `greedy_pruning(model, X, X_hat)` needs both the
original and the edited sequence, so `ledidi(greedy_pruning(...))` is not a meaningful second ordering —
criterion (a)'s "`op_a(op_b(x))` also legal" fails. Two ops with one legal order = "chaining limited to a
small documented subset" = `partial`.

**Disconfirmation:** I looked for a third design op or a documented re-entry of `ledidi` into itself.
`grep -rn "^def \|^class " ledidi/*.py` gives only `ledidi`, `Ledidi`, `greedy_pruning`, `DesignWrapper`
(a multi-model `torch.nn.Module` wrapper, not a sequence op), `MinGap`/losses, and plotting. I read
Tutorials 3 and 7 for repeated `ledidi` calls: `Tutorial_3` calls `ledidi` many times but always on
freshly-constructed in-painted inputs (`X_bar1 = ledidi(model, Xip, y_bar)`, `X_bar2 = ledidi(chrombpnet,
Xip3, y_bar)`), never on a previous `X_bar`. So the ledidi→ledidi chain is type-legal but undocumented,
which does not lift the cell to `yes`.

### 2.7 `dnachisel` = **partial** (prior: partial — confirmed)

DnaChisel has **one** design operation — solve a `DnaOptimizationProblem` — parameterised by a user-composed
set of declarative specifications. Specifications are *aggregated*, not chained: `README.rst:59-66`
```python
    problem = DnaOptimizationProblem(
        sequence=...,
        constraints=[...],
        objectives=[CodonOptimize(species='e_coli', location=(500, 1400))])
```
`builtin_specifications/AllowPrimer.py:18` `class AllowPrimer(SpecificationSet)` nests five other
specifications — but a `SpecificationSet` is still a scoring/repair declaration bundle, not an op whose
output is another op's input. Resolution is a fixed two-stage pipeline (`resolve_constraints()` then
`optimize()`; `CircularDnaOptimizationProblem.py:137,149`) and within a stage the solver treats all
constraints simultaneously (`mixins/ConstraintsSolverMixin.py:61-75`).

What earns `partial` is the officially documented iterative re-feed, `examples/common_scenarios/primers_collection.py:23-43`:
```python
def create_new_primer(existing_primers):
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(length=20),
        constraints=[AvoidHeterodimerization(existing_primers, tmax=3), ...])
    problem.resolve_constraints(); problem.optimize()
    return problem.sequence
...
for i in range(20):
    new_primer = create_new_primer(existing_primers)
    existing_primers.append(new_primer)
```
plus the record round-trip (`mixins/RecordRepresentationMixin.py:17 from_record`, `:76 to_record`), which
makes problem→problem re-feeding a first-class path.

**Disconfirmation:** I searched for a second, distinct design op and for a documented multi-problem chain.
`grep` for files instantiating `DnaOptimizationProblem(` **twice** across `examples/` and `docs/` → **zero
files**. `grep -rn -iE "successive|sequentially|multi-step|two steps|second problem|new problem|then feed|pipeline"
docs/*.rst README.rst` → one hit, `README.rst:213`, and it points at an external webinar video, not an API.
So chaining exists but is single-op, re-feed-only → `partial`, not `yes`.

### 2.8 `mpranator` = **partial** (prior: partial — confirmed, but for narrower reasons than the prior)

The design pipeline is hard-coded in the Django views: `iliasApp/views.py:90-107` runs
`oligo.oligo(...)` (motif insertion) → `part1.getBarCodes(...)` → `part1.createMPRAResultOutput(...)`,
with the SNP view (`views.py:341-373`) hard-coding the analogous `parseSNPs.make_sequence_copies(...)` →
`getBarCodes` → `createMPRAResultOutput`. There is no API; the two tools cannot feed each other.

Two things keep it above `no`:
1. **Documented user-chosen element order.** `templates/iliasApp/part1.html:81` *"Now Pick your ordering of
   the different sequences"*; the drag order is posted as `ordering` (`part1.html:118`) and consumed at
   `part1.py:241` `for orderItem in ordering:` which concatenates adapter1 / restriction1 / Background /
   restriction2 / adapter2 / Barcode in the user's order. This is ordering of *elements within the oligo*,
   not of operations — but the user, not the tool, sets it.
2. **Re-feeding works (ruling R3).** Output records are FASTA (`oligo.py:74`
   `"header": "> Background-%s|%s" % (...)`), and the design input demands FASTA
   (`mycustom.py:32-33` `if ">" not in fasta: raise ValueError("It is not in Fasta Format")`), so a second
   round can insert a second motif set into round-1 oligos.

**Disconfirmation:** I checked for any programmatic composition surface — `iliasApp/views.py`,
`iliasApp/urls.py`, `iliasApp/tasks.py`, `part1.py`, `part3.py`, `oligo.py`, `parseSNPs.py`, `mycustom.py`,
`myfunctions.py` (2482 lines total). There is no library entry point, no ordering of *operations*, and no
`setup.py`/`pyproject.toml`. Had `ordering` permuted design *operations* rather than sequence elements, this
would have been `yes`; it does not.

### 2.9 `mpradesign` (mpradesigntools + designMPRA) = **no** (prior: no — confirmed)

`NAMESPACE` exports exactly two functions:
```
export(processVCF)
export(spread_and_fix_indels)
```
`processVCF` (`R/processVCFfast.R:1099`) is a monolithic driver: it reads the VCF file
(`:1122 vcf_lines = readLines(con = vcf)`), filters and shuffles a barcode pool (`:1184-1215`), then calls
the unexported `processSnp` per row (`:1221 print('Processing SNPs...')`). The Shiny front end is a thin
wrapper: `designMPRA/server.R:17 source('scripts/processVCFfast.R')`, `:114 res = processVCF(vcf = inVCF(), ...)`,
`:162 filename = function() {'out.tsv'}`.

The one documented chain is not a chain of design ops. `R/processVCFfast.R:1405-1406`:
> *"The output is written to the same directory as the input named `"*_fixed.vcf"`. This is ready to be fed
> into processVCF()"*

`spread_and_fix_indels` normalises a **VCF** (splits multi-allelic rows, reformats indels) — its output is
input metadata, not candidate library sequence content, so it is not a design op under the test. One design
op ⇒ nothing to compose.

**Disconfirmation:** re-feeding fails (ruling R3): design output is a TSV of oligos (`out.tsv`), design input
is a VCF file path. Internal helpers that *could* have been composable (`processSnp:210`, `randomly_fix:124`,
`generateInsConstruct:36`, `generateDelConstruct:63`, `change_pattern:90`) are **not exported** — using them
in a different order requires `:::` or editing the package, which is the definition of `no`. `man/` documents
`processSnp.Rd` and `randomly_fix.Rd` but `NAMESPACE` does not export them.

### 2.10 `mutationmaker` = **no** (prior: no — confirmed)

Three independent monolithic solvers, one per workflow, each a single shot from JSON input to a
primer/oligo table. `backend/tasks.py:42-72`:
```python
@celery.task(name='tasks.ssm')
def ssm(ssm_input):  ...  return ssm_solve(input, main_generator, secondary_generator)
@celery.task(name='tasks.qclm')
def qclm(qclm_input): ... return qclm_solve(input)
@celery.task(name='tasks.pas')
def pas(pas_input):  ...  output = pas_solve(input); return output
```
Exposed as three separate endpoints, `api/server_fastapi.py:53,58,63` (`/v1/ssm`, `/v1/qclm`, `/v1/pas`).

**Disconfirmation:** for `partial` I would have needed either a workflow whose output is another workflow's
input, or a user-orderable stage list. Outputs are result tables retrieved by task id and exported as xlsx
(`server_fastapi.py:111-148`); inputs are gene sequence + mutation-site JSON — no re-feed path (R3).
No plugin surface: `backend/mutation_maker/` is 31 modules with the stage order fixed inside
`ssm_solve`/`qclm_solve`/`pas_solve`; there is no PyPI package and no entry-point group.

### 2.11 `codongenie` = **no** (prior: no — confirmed)

One design operation. `main.py:61-66`:
```python
@app.route('/codons')
def get_codons():
    if 'aminoAcids' in request.args:
        codons = _CODON_SELECTOR.optimise_codons(request.args['aminoAcids'],
                                                 request.args['organism'])
```
`codon_genie/codon_selector.py:72 def optimise_codons(self, amino_acids, organism_id)` is the only op that
produces sequence content; `:83 def analyse_codon(self, ambig_codon, tax_id)` is an *analysis* op (it scores
a codon you supply). Client mirror: `codon_genie/client.py:28 def get_codons(self, amino_acids, taxonomy_id)`.

**Disconfirmation:** I noted the one superficially chain-like pair — `optimise_codons` emits ambiguous codons
and `analyse_codon` consumes an ambiguous codon — and rejected it: `analyse_codon` returns analysis, not
library content, so nothing composes into a deeper design. `CodonOptimiser.get_codon_optim_seq` /
`.mutate` exist in `codon_genie/codon_utils.py:62,119` but are not reachable from any route or from
`CodonGenieClient`; using them requires importing internals. Single design op ⇒ `no`.

### 2.12 `biopython` = **no** (prior: no — confirmed)

Biopython supplies no library-design operation, so there is nothing in scope to compose. `Bio/Seq.py`
public methods are container-API mirrors — `count, find, index, split, strip, removeprefix, upper, lower,
translate, complement, reverse_complement, transcribe, back_transcribe, join, replace`, and on `MutableSeq`
`append, insert, pop, remove, reverse, extend`. The docstring makes the provenance explicit,
`Bio/Seq.py:1951-1952`:
```python
    def replace(self, old, new, inplace=False):
        """Return a copy with all occurrences of subsequence old replaced by new.
```
Ruling **R2** excludes these.

**Disconfirmation (searched, absent):**
`grep -rniE "def .*mutagen|def .*mutate|def .*randomize|def .*shuffle|def .*barcode|def .*degenerate|def .*library" --include=*.py Bio/`
→ only `Bio/Nexus/Trees.py:569 def randomize` (phylogenetic trees), `Bio/Phylo/BaseTree.py:758 def randomized`
(trees), and `Bio/motifs/matrix.py:165` / `Bio/motifs/__init__.py:423 def degenerate_consensus` (a consensus
*string*, not a library).
`grep -rniE "def (mutagenize|make_library|generate_variants|design_)" --include=*.py Bio/` → **zero hits**.
`Bio/SeqUtils/` contains only `CheckSum, IsoelectricPoint, MeltingTemp, ProtParam, ProtParamData, lcc` — all
analysis. Tutorial chapters (`Doc/Tutorial/chapter_*.rst`, 27 files) have no library-design chapter.

⚠️ Same false-implicature problem as pydna: Biopython fixes no pipeline. See §4.

### 2.13 `seqpro` = **no** (prior: **partial** → CHANGED)

SeqPro has a genuine, explicit, user-ordered composition combinator —
`python/seqpro/transforms/augmentation.py:12-19`:
```python
class Sequential:
    def __init__(self, *transforms: Callable[..., Any]):
        self.transforms = transforms
    def __call__(self, x: Any) -> Any:
        for t in self.transforms:
            x = t(x)
        return x
```
— but the composable units are not library-design operations. Full public surface,
`python/seqpro/__init__.py` `__all__`: `AA, DNA, RNA, AminoAlphabet, NestedStr, NucleotideAlphabet,
PathLike, SeqType, StrSeqType, alphabets, bed, bin_coverage, cast_seqs, decode_ohe, decode_tokens,
gc_content, gtf, jitter, k_shuffle, length, nucleotide_content, ohe, pad_seqs, rag, random_seqs,
reverse_complement, tokenize, transforms`. The transform classes are
`Jitter, KShuffle, Random, ReverseComplement, Sequential, TMM` (`transforms/__init__.py:4`) — ML
data-augmentation and normalisation. `CLAUDE.md:57` describes `transforms/` as *"Composable transform objects
(`KShuffle`, `ReverseComplement`, `Jitter`, `Sequential`, `TMM`)"*. `README.md` frames the package as
*"a Python package for processing DNA/RNA sequences"*. Under the row's CRITICAL DISTINCTION and ruling **R2**,
composing `pad → ohe → revcomp → tokenize → jitter` is composing general utilities.

**Disconfirmation (searched, absent):** `grep -rniE "library design|variant librar|mutagenesis|barcode" docs/ README.md`
→ **zero hits**. No user-payload insertion/substitution op exists (nothing analogous to
`tangermeme.ersatz.insert(X, motif)`), which is the exact op that lifts tangermeme to `yes`. The nearest
candidates — `k_shuffle`, `jitter`, `random_seqs` — take no design payload and are documented as processing /
augmentation, so even taken together they do not give two composable *design* ops. Had SeqPro shipped one
motif-insertion or mutagenesis op, `Sequential` would have made this an immediate `yes`.

⚠️ Same false-implicature problem as pydna/biopython. See §4.

---

## 3. Result table

| tool | prior | verified | one-line basis |
|---|---|---|---|
| poolparty | yes | **yes** | 47 `Pool -> Pool` methods; both orders verified in venv; `stack` fan-in; `docs/pool.rst:12-14` |
| valiant | no | **no** | closed 7-member enum, all mutators on one shared input, order in `proc_targeton`, no re-feed path |
| mpranator | partial | **partial** | fixed 3-step view pipeline, but user-chosen element `ordering` + FASTA-out/FASTA-in re-feed |
| mpradesign | no | **no** | 2 exported fns, only one a design op; internals unexported; TSV-out / VCF-in |
| oligopoolcalc | yes | **yes** | 9 DataFrame→DataFrame design modules; docs "Chainable"; runnable fan-out+`join` YAML |
| mutationmaker | no | **no** | three independent monolithic solvers, three endpoints, no re-feed |
| codongenie | no | **no** | one design op (`optimise_codons`); second op is analysis-only |
| dnachisel | partial | **partial** | one solve op; specs aggregate not chain; documented in-memory re-feed loop |
| ledidi | no | **partial** ← changed | `ledidi -> greedy_pruning` chain in Tutorial 7; only one legal order |
| tangermeme | partial | **yes** ← changed | pure `(N,4,L)` ersatz ops; `substitute -> greedy_substitution` and greedy/beam fan-out both in Tutorial B6 |
| biopython | no | **no** | no library-design op exists; `Seq` methods are `str`-API mirrors |
| pydna | yes | **no** ← changed | zero library-design ops in `pydna.all`; all composition is cloning simulation |
| seqpro | partial | **no** ← changed | real `Sequential` combinator, but composes processing/augmentation utilities only |

Four changes. Distribution: 3 yes, 3 partial, 7 no.

---

## 4. Why the row needs rewording

The row's `no` is *defined* as "the tool fixes the pipeline order; extending it requires editing the source."
That definition is true of valiant, mpradesign, mutationmaker, and codongenie. It is **false** of pydna,
biopython, and seqpro, which fix no pipeline whatsoever and are among the most compositional packages in
the table — they simply have no library-design operations to compose. Three of the seven `no` cells
therefore carry a claim the primary sources contradict.

This is not a scoring failure that better rating fixes; it is a definitional collision. The row silently
conflates two independent axes:

- **(A)** does the tool expose library-design operations as first-class units at all?
- **(B)** if so, is their composition order under user control?

`no` currently means "(B) is false", but is being used for tools where (A) is false and (B) is undefined.
A referee who authored pydna, Biopython, or SeqPro will read the printed cell as an incorrect statement
about their package, and they will be right.

Secondary hazard: because `yes` reduces to "≥2 pure `T -> T` design functions", any purely functional
numpy/torch perturbation library scores `yes` automatically. That is defensible (tangermeme genuinely does
compose), but it means the row measures *function purity plus presence of design ops*, not the DAG
machinery the anchor text emphasises. If the authors intend the stronger claim, the row must say so — and
must then say it in a way that does not smuggle in PoolParty's memoisation as the discriminator.

## 5. Concrete proposed replacement

**Minimal fix (recommended).** Keep one row, keep the 13 values above, and change the label + legend:

> **Row label:** `user-composable design operations`
> **Question:** Are library-design operations first-class values that the user can chain and branch in an
> order of their choosing?
> **yes** — ≥2 distinct design operations share a carrier type and are composable in either order, with
> fan-out expressible.
> **partial** — design operations exist and chain, but only in a tool-fixed order, only for a small
> documented subset, or only by re-feeding output artefacts.
> **no†** — the tool fixes the order of its design operations in source.
> **no‡** — the tool exposes no library-design operations, so there is nothing to compose (it may still be
> highly composable over general sequence utilities).
>
> Assignment: `no†` = valiant, mpradesign, mutationmaker, codongenie. `no‡` = pydna, biopython, seqpro.

The dagger/double-dagger split costs one legend line and removes the only claim in the row that the
sources contradict. If a two-symbol legend is unacceptable, the alternative is to **split** into
`design_ops_first_class` (yes/no) and `composition_order_user_chosen` (yes/partial/no, scored only where
the first is yes) — cleaner, but it adds a row and partly duplicates the table's scope row.
