# Verified row audit — `library_algebra`

**Auditor pass:** one rater, one threshold, all 13 tools.
**Date:** 2026-08-10
**Row question:** Can whole libraries be combined, sampled, or repeated AS OPERATIONS INSIDE the tool?
**Verdict on the row:** `keep` (scoreable as defined — see §6 for a caveat the caption should carry).

---

## 1. Operational test (fixed BEFORE any tool was opened)

Written to `verified/_library_algebra_test.txt` before inspection:

> For each tool, search its own source + official docs for a **named, documented operation of the tool's
> own design API** whose **operand(s) and result are the tool's own library object** (the type the tool
> produces as its design output: pool / library / variant set / oligo table) and which belongs to one of
> three families:
> **(a) COMBINE** — concatenate / stack / merge / union two or more libraries → one library;
> **(b) SAMPLE** — sample / subsample / downsample a library → smaller library;
> **(c) REPEAT** — repeat / replicate / multiply a library → larger library.
>
> `yes` = all three families. `partial` = one or two families, or a family present only under a
> restricted type/condition. `no` = none; the user must combine or subsample with an **external script**.
>
> **Exclusion, applied uniformly:** generic container operations do not count — Python `list +`,
> `np.concatenate`, `pandas.concat`, `Bio.Align.MultipleSeqAlignment.extend`, `torch.cat`, `dict.update`.
> The operand must be a **library object of the tool** and the operation must be **part of the design API**.

Two gates therefore decide every cell: **(1)** is the operand a library object of the tool's design API?
**(2)** is the operation in one of the three families? A `no` may be earned at either gate; where a tool
fails both I say so, because that makes the cell robust to a referee who disputes one gate.

A third, purely mechanical criterion settles the hardest cases: **does the operation change the member
set (the number of designed sequences)?** Combine adds members, sample removes members, repeat duplicates
members. The anchors force this reading — the anchor NO (VaLiAnT) is missing exactly the *member-count*
concatenation ("the final concatenated library is 62 754"), and the anchor YES (PoolParty `stack`) is
exactly a member-count sum. Operations that leave the member set untouched (column-wise joins,
within-member concatenation) are recorded as **near-misses** rather than silently ignored.

---

## 2. Scores

| tool | prior | **verified** | one-line reason |
|---|---|---|---|
| poolparty | yes | **yes** | `stack` + `sample` + `repeat` as `Operation` subclasses; `+`, `*`, `[]` on `Pool` |
| oligopoolcalc | no | **partial** ⬆ | `join` merges two oligo-pool tables (design API), but column-wise only — never adds rows; no sample, no repeat |
| valiant | no | **no** | 2 CLI commands; per-targeton CSVs; user concatenates |
| mpranator | no | **no** | 4 form endpoints, one FASTA blob each; no library object |
| mpradesign | no | **no** | 2 exported R functions; combining runs is `bind_rows` by the user |
| mutationmaker | no | **no** | 4 celery tasks, independent jobs; no operation takes two designs |
| codongenie | no | **no** | 4 REST routes; one codon at a time; no library object |
| dnachisel | no | **no** | 55-name `__all__`, no library object; own example uses `list.append` + `SeqIO.write` |
| ledidi | partial | **no** ⬇ | no library object (output is a tensor); `n_repeats`/`n_samples` are generation-time params of one call |
| tangermeme | no | **no** | no such function in the entire public API |
| biopython | no | **no** | no library object anywhere in `Bio/`; 0 hits for any sampling function |
| pydna | no | **no** | 15-name headline API; all classes single-sequence |
| seqpro | partial | **no** ⬇ | `rag.concatenate` is a generic-container op that *cannot* change the sequence count; no sample, no repeat |

**Changes from the prior per-tool pass: 3** (seqpro `partial`→`no`, ledidi `partial`→`no`,
oligopoolcalc `no`→`partial`).

---

## 3. The three flagged cells, resolved

### 3.1 seqpro — `partial` → **no**

Fails **both** gates, which is why this is high confidence.

*Gate 2 (decisive).* The only candidate in the entire package is `sp.rag.concatenate`
(`python/seqpro/rag/_ops.py:341`). Its docstring:

> `"""Concatenate Ragged arrays along the ragged axis.`
> `Given N ``Ragged`` arrays that share the same number of groups and the same leading fixed dimensions,`
> `concatenate their per-group element sequences so that each output group contains the elements of`
> `` rags[0] `` followed by `` rags[1] ``, …, `` rags[-1] ``.`

and the body enforces it (`_ops.py:377-380`):

```python
ax = axis % len(ref.shape)
if ax != ref.rag_dim:
    raise ValueError(
        f"concatenate only supports the ragged axis ({ref.rag_dim}), got {axis}"
    )
```

So the operation **cannot** be applied to the group axis: it appends elements *inside* each sequence and
requires the two inputs to have the same number of groups. It can never increase the number of sequences.
It is not a library concatenation in any reading — not even a generic-container one. There is no `sample`
and no `repeat`: the complete function list of `python/seqpro/_modifiers.py` is `reverse_complement`,
`k_shuffle`, `bin_coverage`, `jitter`, `random_seqs`, `normalize_coverage`.

*Gate 1.* `Ragged` is a generic array container, not a library:
`rag/_core.py`: `class Ragged(NDArrayOperatorsMixin, Generic[RDTYPE_co])` with the one-line docstring
`"""A non-branching ragged array with a single ragged axis (Spec A)."""`, and `concatenate` "Supports any
fixed-itemsize numeric dtype (e.g. ``int32``, ``float32``)" — it carries coverage tracks and scores as
readily as sequences. SeqPro also has no design API for it to be part of.

*Documentation status.* `concatenate` is **not** in the API reference. `docs/api/index.md` documents 15
members (`bin_coverage, cast_seqs, decode_ohe, decode_tokens, gc_content, jitter, k_shuffle, length,
nucleotide_content, ohe, pad_seqs, random_seqs, reverse_complement, tokenize`) and `docs/api/ragged.md`
documents only `Ragged`, `is_rag_dtype`, `lengths_to_offsets`. Conversely the docs use NumPy for the user
side: `docs/ragged.md:115` — `ohe_data = np.concatenate([sp.DNA.ohe(sp.cast_seqs(s)) for s in cds_seqs])`.

*Checked and absent:* `grep -rniE "^\s*(def|fn|pub fn) [a-z_]*(concat|stack|sampl|repeat|replicat|tile|
merge|union|combine)[a-z_]*\("` over `python/` and `crates/` returns exactly two hits, `_ops.py:341
concatenate` and `crates/seqpro-core/src/ragged.rs:134 ragged_concat` (its Rust kernel). Also read
`python/seqpro/__init__.py` `__all__` (28 names), `rag/__init__.py` `__all__` (12 names),
`transforms/__init__.py` `__all__` (`TMM, Jitter, KShuffle, Random, ReverseComplement, Sequential`).

### 3.2 ledidi — `partial` → **no**

Fails **both** gates.

*Gate 1.* There is no library object. `ledidi()` returns a `torch.Tensor`
(`ledidi/ledidi.py:186-190`: `y: torch.Tensor, shape=(*ny, *n_repeats, n_samples, n_channels, length)`).

*Gate 2.* The prior `partial` rested on `n_repeats`, `n_samples` and the affinity-catalog axis. All three
are **parameters of a single generation call**, not operations applied to a library:

- `n_repeats`: *"The number of times to run the Ledidi procedure."* Implementation
  (`ledidi.py:222-227`) constructs a **fresh** `Ledidi(...)` designer per repeat with a distinct seed and
  runs `fit_transform` again. This produces *k independent designs*; it does not replicate a library.
- `n_samples`: *"The number of samples to draw from Ledidi after the optimization process."*
  Implementation (`ledidi.py:234-236`): `X_bar_ = torch.cat([designer(X) for _ in range(n_iter)], dim=0)
  [:n_samples]` — forward passes drawing new sequences from a learned Gumbel-softmax distribution. This is
  *generation*, not sampling of an existing library, and the operand (a `torch.nn.Module`) is not a library
  and is not returned.
- the catalog axis: `X_bar[i] = torch.stack(X_bar[i])` / `X_bar = torch.stack(X_bar)`
  (`ledidi.py:240,242`) — internal `torch.stack`, an excluded generic container op.

*No combine.* The one thing in Ledidi that "concatenates" is `DesignWrapper`
(`ledidi/wrappers.py:67`): `return torch.cat([model(X).clone() for model in self.models], dim=-1)` —
it concatenates **model predictions**, and its module docstring says so: *"takes a list of models, runs the
input sequence through each, and concatenates their predictions into a single output tensor."* Nothing
combines two design outputs.

*Checked and absent:* the complete documented API is 8 names, read from `docs/api/*.rst`
(`automodule` `:members:` lists): `ledidi`, `Ledidi`, `MinGap`, `plot_loss`, `plot_history`, `plot_edits`,
`greedy_pruning`, `DesignWrapper`. `grep -rn "torch.cat\|torch.stack\|\.repeat(\|repeat_interleave\|
torch.tile" ledidi/*.py` returns 6 hits, all internal (listed above, plus two in `plot.py`).

### 3.3 biopython — **no** (prior `no` upheld, and now consistent with the two demotions)

The prior pass's worry — that `biopython = no` was inconsistent with seqpro/ledidi `partial` — is resolved
by demoting those two rather than promoting Biopython. Biopython's own position:

- **Gate 1 fails absolutely.** `grep -rn "^class [A-Za-z]*\(Library\|Pool\)" --include=*.py Bio/` → **0
  hits**. There is no design layer and therefore no library object.
- **Gate 2, sampling.** `grep -rn "^\s*def [a-z_]*sampl" --include=*.py Bio/` → **0 hits** across all
  2 440 files. Same for `^\s*def [a-z_]*\(repeat\|replicat\)` → 0 hits.
- **Nearest primitives, both excluded and both named in the row definition.**
  `Bio/SeqRecord.py:915 def __add__` — *"Add another sequence or string to this sequence"*: concatenates a
  record with a record/`Seq`/`str`, i.e. **sequence-level**, member set of size one throughout.
  `Bio/Align/__init__.py:438 def extend(self, records)` — *"Add more SeqRecord objects to the alignment as
  rows. They must all have the same length as the original alignment."* This *does* increase the row count,
  which makes it a stronger member-set operation than anything in seqpro; but the operand is an
  **alignment**, not a library, it requires equal-length aligned records, and the row definition names it
  explicitly as not qualifying. `Bio/Align/__init__.py:562 def __add__` — *"Combine two alignments with the
  same number of rows by adding them"* (by column) — the same-cardinality horizontal merge shape again.

---

## 4. The anchors, verified rather than assumed

### PoolParty — **yes** (anchor holds; not softened)

Three first-class `Operation` subclasses, one per family, all exported from the top level
(`src/poolparty/__init__.py:124-126, 225-233`) and each with its own documentation page:

| family | factory | class | file:line | doc page |
|---|---|---|---|---|
| combine | `stack(pools, ...) -> T` | `StackOp` | `state_ops/stack.py:19`, `:51` | `docs/operations/stack.rst` |
| sample | `sample(pool, num_seqs=…, seq_states=…, with_replacement=…) -> T` | `SampleOp` | `state_ops/sample.py:18`, `:74` | `docs/operations/sample.rst` |
| repeat | `repeat(pool, times, ...) -> T` | `RepeatOp` | `state_ops/repeat.py:16`, `:55` | `docs/operations/repeat.rst` |

Docstrings, verbatim: `stack` — *"Create a pool by stacking multiple input pools state-wise (disjoint
union)."* returning *"A Pool whose states are the disjoint union of all input pools' states"*;
`sample` — *"Sample states from a pool."* returning *"A Pool containing the sampled states from the input
Pool"*; `repeat` — *"Repeat a pool's states a specified number of times."* returning *"A new Pool with
``times`` as many states as the input pool."*

Operator sugar on the library object itself (`src/poolparty/pool.py:169-189`):

```python
def __add__(self, other: Pool_type) -> Self:      # "Stack two pools (union of states via sum_counters)."
    from .state_ops.stack import stack
    return stack([self, other])
def __mul__(self, n: int) -> Self:                # "Repeat this pool n times (repeat states)."
    from .state_ops.repeat import repeat
    return repeat(self, n)
def __getitem__(self, key: Union[int, slice]) -> Self:   # "Slice this pool's states (not sequences)."
    from .state_ops.state_slice import state_slice
    return state_slice(self, key)
```

Closure and type-genericity checked, not assumed: every factory is `@beartype`d `(T, ...) -> T` with
`T = TypeVar("T", bound=Pool)` and returns `pool_class = type(pool)`, so the operations are **not**
restricted to `DnaPool` — the rubric's "restricted type" downgrade does not apply.

Behaviour confirmed in the read-only venv (`PYTHONDONTWRITEBYTECODE=1
/mnt/c/.../poolparty-statecounter/.venv/bin/python`):

```
a 3   b 2
stack   -> 5      a+b     -> 5
repeat3 -> 9      a*3     -> 9
sample2 -> 2      a[0:2]  -> 2
generate_library(a+b) -> DataFrame with 5 rows: AAAA CCCC GGGG TTTT ACGT
```

and the docs state the arithmetic independently (`docs/operations/library_size.rst`): *"``stack`` (or the
``+`` operator) places its input pools side by side … the resulting ``num_states`` is the sum of the
inputs' ``num_states``."*

*Disconfirmation sought:* I looked for the failure modes that would demote this — ops implemented only on
one pool subclass (no: `TypeVar` + `type(pool)`), ops that exist but are undocumented (no: three `.rst`
pages plus `library_size.rst`), or ops that are helpers rather than graph nodes (no: `StackOp`, `SampleOp`,
`RepeatOp` are `Operation` subclasses with `factory_name` and `build_pool_counter`, so they are nodes in
the pool DAG). None applied. `yes` stands.

### VaLiAnT — **no** (anchor holds)

- Two subcommands only: `src/valiant/cli.py:110-111` — `main.add_command(sge)` / `main.add_command(cdna)`.
  Their full option lists (`sge_cli.py:26-90`, `cdna_cli.py`) contain nothing that takes a second library.
- Output is **per targeton**: `src/valiant/meta_table.py:340-352` —
  `def get_path(suffix, ext): return self.cfg.get_output_file_path(f"{targeton_name}_{suffix}.{ext}")`,
  then `unique_fp = get_csv_path('unique')`, `meta_fp = get_csv_path('meta')`.
  README §"Unique oligonucleotides file" (line 807): *"file containing only the label and the sequence of
  the oligonucleotides generated for **any given targeton**, where the sequences are unique."*
- `grep -rniE "^\s*(def|class) [a-z_]*(concat|stack|sampl|repeat|replicat|tile|merge|union|combine)"`
  over `src/valiant` → **0 hits**.
- Near-miss recorded: the `_unique.csv` de-duplication (`meta_table.py:359,591,676-681`) is a
  generation-time dedup within one run's own output, not an operation on library objects.

---

## 5. The remaining nine cells

**oligopoolcalc — `no` → partial.** This is the row's one genuine judgment call; both readings are given
so the authors can flip it with the reasoning already written.

*Why partial.* `oligopool.join` is a documented **Design Mode** function whose two operands are both
oligo-pool variant tables — the tool's own library representation — and whose result is one such table
(`oligopool/join.py:13-46`):

> *"Join two oligo tables on `ID` and insert new columns from `other_data` into the `input_data` column
> order. **Useful for recombining parallel branch outputs into a single design table.**"*

`docs/docs.md:702-704`: *"**When to use it**: Recombining parallel design branches (for example, from a
YAML CLI DAG) back into one design table, or reconciling two independently designed DataFrames."*
It is in `__api__` (`oligopool/__init__.py:6-32`) and has its own `docs/api.md` entry. It therefore clears
the row's BEWARE gate — the operands *are* library objects and the operation *is* part of the design API —
and "merge" is one of the three verbs the rubric names. One family, restricted → `partial`.

*Why the cell must say what it does not do.* `join` is **column-wise**, and the docs are explicit
(`docs/docs.md:726-728`):

> *"`join` is an inter-table operation (contrast with `merge`/`revcomp`, which operate within one
> DataFrame)."*
> *"`input_data` and `other_data` must contain the same `ID` set; **`join` never creates or drops rows**
> (ID mismatches are an error)."*

So it cannot union two variant sets. Union of members is done by the user with pandas — in the official
docs (`docs/docs.md:219-224`, Patch Mode):

```python
# Later: append new variants; mark missing barcodes with '-'
df = pd.concat([df, pd.DataFrame({'ID': ['V3'], 'Variant': ['TTAA' * 10], 'BC': ['-']})],
               ignore_index=True)
```

There is no sample and no repeat: the 22-name `__api__` and the 20 documented functions in
`docs/api.md` (`barcode, primer, motif, spacer, background, merge, revcomp, join, final, split, pad,
compress, expand, index, pack, acount, xcount, lenstat, verify, inspect`) contain neither, and
`grep -n -i "subsample\|def sample" docs/api.md` → 0. Other near-misses recorded: `merge` (*"collapse
columns into single element"*) and `final` (*"concatenate into synthesis-ready oligos"*) are
within-variant; `compress`/`expand` change the row count but are IUPAC de/re-generation, not
combine/sample/repeat; the `pipeline` YAML DAG (`docs/docs.md:1631+`) only schedules the same 20 commands.
*If the authors prefer the strict member-set reading, this cell becomes `no`* — the deciding fact either
way is the quoted "never creates or drops rows".

**mpranator — no.** `iliasApp/urls.py` exposes four design endpoints (`MPRAResults/`,
`MPRAResults/SNPs/`, `TransmutationResults/`, plus `API/MPRA/SNP`); each POST returns one FASTA blob and
none accepts a second design. No library object exists to operate on. Live documentation
(genomegeek.com/MPRA/documentation/, fetched 2026-08-10) describes five sections and no combine, merge,
subsample or library-repeat feature; the only replication is the generation-time form field *"Number of
barcodes per sequence (replicates)."* — recorded as a near-miss of the same class as ledidi's `n_repeats`.
`grep -rniE "^\s*(def|class) [a-z_]*(concat|stack|sampl|repeat|merge|union|combine)"` over the repo → the
only hit is `iliasApp/models.py:5 class SampleCounting(models.Model)`, a Django visit counter.

**mpradesign — no.** `NAMESPACE` exports exactly two functions: `export(processVCF)` and
`export(spread_and_fix_indels)`. `R/processVCFfast.R:1099 processVCF = function(vcf, nper, …, outPath)`
runs one VCF to one output; roxygen (`:195`) *"@param nper The number of barcoded sequences to be
generated per allele per SNP"* is a generation-time replication parameter (near-miss). `bind_rows` appears
only as `importFrom(dplyr,bind_rows)` used internally on per-SNP tibbles — the generic container op the
row excludes; combining two runs is the user's `bind_rows`. No sample of a library (the sampling in
`randomly_fix` draws **barcodes** from a barcode pool: *"@param snp a tibble containing the VCF
information for one SNP as well as a barcode pool to sample from"*, `:193`) and no repeat operation.

**mutationmaker — no.** The whole compute surface is four Celery tasks (`backend/tasks.py:42,55,62,68`):
`ssm`, `qclm`, `species_table`, `pas`. Jobs are independent and there is no task, function or endpoint
that takes two designs. Near-miss recorded and inspected: normalised **physical** mixing ratios,
`backend/mutation_maker/generate_oligos.py:136-139` — `total_ratio = sum(oligo.ratio for oligo in
oligos)` … `oligo.ratio /= total_ratio` — a wet-lab mixture specification, not an operation on library
objects. `pas_output.py:90 combine_oligos_list` is an internal output-formatting helper that pairs oligos
with mutation indices and returns a plain `List`. The `union`/`combine`/`merge` hits in
`degenerate_codon.py:43,395`, `qclm.py:697`, `site_split.py:55`, `pas_solution.py:202` all operate on
codons/amino-acid subsets/primer solutions, not libraries.

**codongenie — no.** Four routes, read from `main.py:40-62`: `/`, `/organisms/`, `/organisms/<term>`,
`/codons`. `get_codons()` calls `_CODON_SELECTOR.optimise_codons(request.args['aminoAcids'],
request.args['organism'])` or `analyse_codon(...)` — one codon at a time. No library object, so nothing to
combine, sample or repeat; `grep -rniE "^\s*(def|class) [a-z_]*(concat|stack|sampl|repeat|merge|union|
combine)"` over the repo → 0 hits.

**dnachisel — no.** Complete public API read from `dnachisel/__init__.py:75-129` (`__all__`, 55 names):
one problem class (`DnaOptimizationProblem`, `CircularDnaOptimizationProblem`), specifications, patterns
and biotools — no library object and no concat/sample/repeat. The shipped example proves the combining is
external (`examples/common_scenarios/genome-wide-optimization.py:10,29,32`):

```python
optimized_records = []
...
        optimized_records.append(problem.to_record(record_id=protein_id))
SeqIO.write(optimized_records, "genome_wide_domestication.fa", format="fasta")
```

Near-misses recorded: `dnachisel/reports/constraints_reports/constraints_breaches_dataframe.py:30
constraints_breaches_dataframe(constraints, sequences)` consumes a user-held set for QC but neither builds
nor transforms a library; `Location.merge_overlapping_locations` and `MutationChoice.merge_with` merge
locations/choices inside one sequence.

**tangermeme — no.** `tangermeme/__init__.py` contains only `__version__`; the documented surface is the
18 pages in `docs/api/`. A grep of every `^def ` in `tangermeme/*.py` plus `tangermeme/design/*.py` for
`sampl|concat|stack|repeat|merge|combin|pool` returns exactly one hit, `deep_lift_shap.py:172 _maxpool`
(a private backward hook). Near-misses inspected: `product.py:30 apply_pairwise` / `:197 apply_product`
apply *a function* over the cartesian product of `X` and `args` and return **predictions**, not a library
(*"Apply a function on the cartesian product between X and args"*); `utils.py:649 chunk` tiles one long
sequence into fixed-length blocks; `match.py` selects GC-matched background loci from a **genome**, not
from a designed library; `design/_substitute.py:8 _fast_tile_substitute` is a private kernel. Combining
tensors is the user's `torch.cat`.

**pydna — no.** Headline API read from `src/pydna/all.py:26-41` (`__all__`, 15 names): `Anneal, pcr,
Assembly, genbank, Genbank, Dseqrecord, Dseq, read, read_primer, parse, parse_primers, primer_design,
assembly_fragments, eq, gbtext_clean`. `grep -rn "^\s*def [a-z_]*\(sampl\|repeat\|concat\|stack\|subset\)"
--include=*.py src/pydna/` → one hit, `genbankfixer.py:288 concat_dict(dlist)` (GenBank text repair).
Every sequence class is single-sequence (`Dseqrecord`, `Amplicon`, `Contig`, `Primer` — full class list
grepped); `Dseqrecord.__add__` is ligation of two records, and `Assembly` recombines **fragments** into
constructs, so both are construct assembly, not member-set algebra. `assembly_fragments(list) -> list`
adds primer tails and passes a plain Python list through.

---

## 6. Row-level judgement, stated for the authors

`row_verdict = keep`. Every cell fell out of the fixed test without a coin-flip, including the two
demotions and the one promotion, so the row is scoreable on one scale.

Two caveats the authors should carry into the caption or the response letter, because a referee will
raise them:

1. **The row is partly entailed by `library_first_class_object`.** For nine tools (valiant, mpranator,
   mpradesign, mutationmaker, codongenie, dnachisel, ledidi, biopython, pydna) the `no` follows from there
   being no library object at all — the cell is a corollary, not independent evidence. The row carries
   genuinely independent information only for the tools that *do* hold a library-like object: poolparty
   (yes), oligopoolcalc (partial), seqpro (no), and marginally tangermeme (no). Saying this once in the
   text is cheaper than being told it by a referee.
2. **The distribution is 1 yes / 1 partial / 11 no.** That is what the sources say, but it means the row
   reads as a PoolParty self-portrait unless the *mechanism* is named. The mechanism worth stating is
   sharp and checkable: the near-misses across the field are all **same-cardinality** merges — oligopool
   `join` ("never creates or drops rows"), seqpro `rag.concatenate` (ragged axis only), Biopython
   `MultipleSeqAlignment.__add__` ("same number of rows"), plus **generation-time** replication knobs
   (ledidi `n_repeats`, MPRAnator "barcodes per sequence", mpradesigntools `nper`). What is absent
   field-wide is an operation that *changes the member set of a held library*. If a shorter row is wanted,
   `library_member_set_algebra` with that phrasing would say the same thing less contestably.

---

## 7. Files, greps and doc pages consulted

Source snapshots downloaded 2026-08-10 to a scratch dir (read-only analysis; nothing installed):
SeqPro `main`, ledidi `master`, tangermeme `master`, VaLiAnT `develop`, oligopool `master`, DnaChisel
`master`, pydna `master`, CodonGenie `master`, ra100/Mutation_Maker `master`, MPRAnator `master`,
mpradesigntools `master`, designMPRA `master`, biopython `master`.

Read directly (beyond the file:line citations above): `poolparty/src/poolparty/state_ops/__init__.py`,
`pool.py:150-215`, `__init__.py:120-360`, `docs/operations/` (60 pages listed; `stack.rst`, `sample.rst`,
`repeat.rst`, `library_size.rst` read), `seqpro/python/seqpro/{__init__,_modifiers}.py`,
`rag/{__init__,_core,_ops}.py`, `docs/api/{index,ragged}.md`, `docs/ragged.md`,
`ledidi/ledidi/{ledidi,wrappers,plot,pruning,losses}.py`, `ledidi/docs/api/*.rst`, `docs/index.rst`,
`valiant/src/valiant/{cli,sge_cli,meta_table}.py`, `valiant/README.md`, `oligopool/oligopool/{__init__,
join,merge,descriptions}.py`, `oligopool/docs/{docs,api}.md`, `dnachisel/dnachisel/__init__.py`,
`dnachisel/examples/`, `pydna/src/pydna/all.py`, `tangermeme/tangermeme/{__init__,product,utils,match,
ersatz}.py`, `tangermeme/tangermeme/design/__init__.py`, `tangermeme/docs/api/`,
`biopython/Bio/{SeqRecord,Align/__init__,motifs/__init__}.py`, `mutationmaker/backend/tasks.py`,
`backend/mutation_maker/{pas_output,generate_oligos}.py`, `codongenie/main.py`,
`mpranator/iliasApp/{urls,forms,views}.py`, `mpranator/part1.py`,
`mpradesigntools/{NAMESPACE,R/processVCFfast.R}`.

Named absence greps (each returning zero relevant hits, i.e. every `no` is "I searched and it is absent"):

```
# per repo, function/class declarations in any of the three families
grep -rniE "^\s*(def|class|fn|pub fn) [a-z_]*(concat|stack|sampl|repeat|replicat|tile|merge|union|combine)[a-z_]*\(" <repo>
# biopython, recursive over Bio/ (2440 files)
grep -rn "^\s*def [a-z_]*sampl"                       --include=*.py Bio/   -> 0
grep -rn "^class [A-Za-z]*\(Library\|Pool\)"          --include=*.py Bio/   -> 0
grep -rn "^\s*def [a-z_]*\(repeat\|replicat\)"        --include=*.py Bio/   -> 0
# oligopool documented API
grep -n -i "subsample\|def sample" docs/api.md                              -> 0
```

Live pages fetched: genomegeek.com/MPRA/documentation/ (2026-08-10).
Behavioural check: PoolParty only, in the existing read-only venv, `PYTHONDONTWRITEBYTECODE=1`.
