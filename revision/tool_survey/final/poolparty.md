# PoolParty — FINAL capability record

**Slug:** `poolparty`
**Full name:** PoolParty — "streamlined design of DNA sequence libraries in Python"
**Tier:** **subject of the paper** (this is the authors' own tool; every value below is written to
survive a hostile referee, not to flatter)
**Authors:** Zhihan Liu, Aidan Cordero, Justin B. Kinney (Cold Spring Harbor Laboratory)
**Version assessed:** PyPI `poolparty` 0.1.1 (wheel + sdist uploaded 2026-04-06T21:10 UTC);
working tree at commit `1bb0179` "fix RTD citation", 2026-04-07.
**Canonical repository:** `https://github.com/jbkinney/poolparty-statetracker`
(the local working copy is named `poolparty-statecounter/` — see §Unresolved / housekeeping).
**Assessed:** 2026-08-10.

**Record status:** FINAL. Built from `extractions/poolparty.md` + `reviews/poolparty.md`, with every
disputed fact re-derived by me from source, docs, `.gitignore`, `pyproject.toml`, the CI workflow,
the live PyPI JSON API, and **24 read-only executions** against
`poolparty-statecounter/.venv/bin/python` (run with `PYTHONDONTWRITEBYTECODE=1`; nothing written
outside `revision/tool_survey/`; the repository was not modified).

**Adversarial review outcome:** 26 of the 30 values it reviewed `supported`, 0 `overstated`, 0
`wrong` (the review's own tally, `reviews/poolparty.md:394-397`: `supported` 26 · `understated` 1 ·
`unsupported as encoded` 3).
The non-`supported` verdicts were all *encoding* or *understatement* defects, not value errors:
`interface` / `last_activity` (values were the uninformative bare string `yes`), and
`documented_examples` (bare `yes` **and** two counts wrong low). The reviewer additionally scored
`combinatorial_motif_place` and `negative_control_generation` `supported` while judging their
evidence *too modest*. All are corrected below.
**No capability value changed.** No genuine extractor/reviewer disagreement over a capability value
exists; the disagreements were over line numbers and file counts, and I adjudicated each one
myself — see §Adjudicated disagreements.

---

## Sources consulted

| Kind | Path / ref | Note |
|---|---|---|
| source | `poolparty/src/poolparty/` — **102** `.py` files | re-counted |
| source | `pool.py`, `operation.py`, `generate_library.py` (336 lines) | read in full |
| source | `base_ops/{mutagenize,get_barcodes,flip,shuffle_seq,from_iupac,filter_seq}.py` | read |
| source | `orf_ops/{mutagenize_orf,reverse_translate,translate}.py`, `codon_table.py` | read |
| source | `state_ops/{stack,sample,repeat,state_slice,state_shuffle,sync}.py`, `fixed_ops/join.py` | read |
| source | `region_ops/{apply_at_region,replace_region,region_scan,region_multiscan}.py`, `multiscan_ops/` | read |
| source | `fixed_ops/{from_fasta,score}.py`, `pool_mixins/{export_mixin,filter_mixin}.py` | read |
| readme | `poolparty/README.md` | |
| docs | `docs/index.rst`, `quickstart.rst`, `pool.rst`, `regions.rst`, `metadata/{design_cards,naming,styling}.rst` | |
| docs | `docs/operations/` — **62** `.rst` files (**50** per-operation + 12 concept/index); 74 `.rst` total under `docs/` | re-counted |
| docs | `docs/tutorials/{index,dms_gb1,mpra_regulatory_grammar}.rst` | numbers re-derived |
| tests | `poolparty/tests/` — **77** `.py` files, **2,906** `test_` functions | re-counted |
| changelog | `poolparty/CHANGELOG.md` (`[Unreleased]`, 0.1.1, 0.1.0) | |
| examples | `poolparty/examples/{README.md, dms_gb1.ipynb, mpra_regulatory_grammar.ipynb, spliceai_surrogate.ipynb}` | **gitignored (`.gitignore:70`), NOT distributed** — `git ls-files` returns nothing |
| other | `poolparty/pyproject.toml`, `poolparty/LICENSE`; **monorepo root**: `LICENSE`, `CITATION.cff`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `.pre-commit-config.yaml`, `.readthedocs.yaml`, `.github/workflows/test.yml` | path corrected |
| other | companion package `statetracker/` (state algebra engine, MIT, on PyPI) | |
| pypi | `https://pypi.org/pypi/poolparty/json` — pulled live 2026-08-10 | |
| live | 24 read-only executions against the repo `.venv` (Python 3.12) | this memo |
| prior | `extractions/poolparty.md`, `reviews/poolparty.md` | |

---

## What PoolParty is

A Python package in which a DNA sequence library is specified as a **directed acyclic graph** of
`Pool` nodes and `Operation` edges. `pool.num_states` is an exact state count and `pool.seq_length`
is exact whenever it is statically determined (`None` for variable-length pools); both are available
at DAG-construction cost; sequences are produced one row at a time only when
`generate_library`/`to_df`/`to_file` is called. There are 39 operation factories with `Pool`/`DnaPool`
method forms and ~50 documented operations in total; the state algebra is delegated to the companion
`statetracker` package.

---

# Capability assessment

Format: **`key` — value.** Evidence. *Source.*

## BLOCK A — library specification

### `library_first_class_object` — **yes**

`Pool` is a returned, holdable, inspectable, passable object. Verified at file:line in `pool.py`:
identity and naming (`named()` 195, `copy()` 200, `deepcopy()` 225 — the latter recursively copies
the whole upstream DAG); inspection (`num_states` 134, `parents` 139, `seq_length` 144, `regions`
149, `has_region` 153); operator overloads (`__add__`→stack 169, `__mul__`→repeat 175,
`__getitem__`→state_slice 185); presentation (`print_library` 282, `print_dag` 385).

Pools travel as **arguments to other operations**, not merely as chain receivers:
`insertion_multiscan(..., insertion_pools=[m1, m2])`, `replace_region(..., content_pool=barcode_pool)`,
`replacement_scan(cryptic, ...)`, `stack([...])`, `sync([...])` — all exercised live for this memo.

*Caveat, stated but do not over-concede it.* `pool.py:48` raises
`RuntimeError("Pools must be created inside a Party context…")`. **A default `Party` is created at
import time** (`__init__.py:269 _init_default_party()`, verified live:
`pp.get_active_party() is not None` → `True`), and the failure is **not** reachable by either of the
routes previously claimed here: `pp.init()` creates, enters and activates a fresh default `Party`
(`party.py:62-89`), and exiting a `with pp.Party()` block restores the previously active party
(`party.py:207-226`) — both re-verified live, `pp.from_seq` succeeds after each. The real cost is architectural, not
ergonomic: there is an implicit module-level registry, which complicates concurrent or
library-embedded use.

*Source:* `src/poolparty/pool.py`, `src/poolparty/__init__.py:269`, `src/poolparty/party.py:62-89,207-226`,
live execution.

### `composable_operations` — **yes**

Most operations are `Pool(s) → Pool` and are exposed both as a module function and as a `Pool` method
(`__init__.py:331 _POOL_FACTORY_MAP` (31 entries, one of which is `generate_library`, not an
operation) / `:373 _DNAPOOL_FACTORY_MAP` (9 entries) — **39 operation factories with method forms**,
verified live; the loops at `:390/395` copy docstrings onto methods the mixins already supply, they
do not bind factories). Source and multi-input operations — `from_seq`, `from_seqs`, `from_fasta`,
`get_barcodes`, `stack`, `join`, `sync`, `region_scan`, `region_multiscan` — are module functions
only, and `sync` returns `None` rather than a `Pool`. The result is an
arbitrary DAG, topologically sorted once per `generate_library` call
(`generate_library.py:220-237 _topo_sort_operations`; `visit()` dedupes via a `visited` set, so a
shared node appears exactly once in `sorted_ops`).

Three properties, in increasing order of what a referee will actually test:

1. **Nesting.** `region_ops/apply_at_region.py` takes a `transform_fn` and applies an arbitrary
   *sub-pipeline* — not a single op — to the contents of a named region. A two-op sub-pipeline
   (`deletion_scan(mutagenize(p, …), …)`) yields 168 states with composed names `m_00.d_0`,
   `m_00.d_1`, … *Cite the docstring (`Callable[[Pool], Pool]`), not the annotation, which is a bare
   `Callable`.*
2. **Branching.** A three-branch `stack` over `mutagenize`, `deletion_scan`+`replace_region`, and
   `get_barcodes` renders correctly under `print_dag`.
3. **TRUE SHARED DAG NODE — this is the row's whole point and it is now in the record.**
   One `mutagenize` pool feeding *both* an `rc` branch and a second `mutagenize` branch, stacked.
   Reproduced by me verbatim:

   ```
   base = pp.from_seq('ACGT')
   m    = base.mutagenize(num_mutations=1, mode='sequential', prefix='m')   # 12 states
   lib  = pp.stack([pp.rc(m), m.mutagenize(num_mutations=1, mode='sequential', prefix='m2')])
   # -> 12 + 144 = 156 states
   ```
   `print_dag` shows `pool[1]` under both branches; `_topo_sort_operations`' `visited` set puts
   `op[1]` in the list once, and `seq_cache` (`generate_library.py:260`, populated at `:292`,
   consumed at `:273`) computes it **once per row**. A tree cannot do this. Nesting and branching
   alone would not have distinguished PoolParty from a linear pipeline with a fan-out.

*Undocumented composability restrictions (belong in limitations, not here):* `mode='sequential'` is
refused with `mutation_rate` (`base_ops/mutagenize.py:157-159`) and with non-uniform `mutation_type`
(`orf_ops/mutagenize_orf.py:205-207`); `flip` refuses `mode='fixed'` (`base_ops/flip.py:122-125`);
`shuffle_seq` supports random mode only (docstring, `shuffle_seq.py:44`).

*Source:* `generate_library.py:220-237, 260, 273, 292`; `region_ops/apply_at_region.py`;
`__init__.py:331,373,390,395`; live execution.

### `lazy_generation` — **yes** (with a precisely bounded caveat — evidence amended)

**Mechanism, read line by line.** `generate_library.py:115` `while len(rows) < num_seqs:` →
`_compute_one` → `:265` `pool.state.value = global_state % pool.state.num_values`, which propagates
through the state DAG; then a single pass over `sorted_ops` calling `op.compute()`. Nothing
enumerates the library. `ExportMixin.to_file` (`pool_mixins/export_mixin.py:222`) streams in
`chunk_size` blocks straight to csv/tsv/fasta/jsonl (+`.gz`), docstring: *"Generates sequences in
chunks and writes them incrementally to avoid loading the entire library into memory."*

**Measured (authors' own table, `examples/README.md`):** DAG construction 0.03 / 0.05 / 0.24 s
against 15.08 / 10.12 / 170.32 s to generate the same libraries — a **202-710x** gap (SpliceAI 503x,
MPRA 202x, DMS 710x). *The README's own "500-3000x" summary does not follow from its own table — it
cross-pairs the fastest construction with the slowest generation. Quote 202-710x.*

**The evidence sentence "Nothing materialises the library" was an overreach and has been removed.**
What is eager, stated exhaustively so no referee finds it first:

- **`get_barcodes` pre-generates and stores every barcode at DAG-construction time.** The authors'
  own class docstring says so (`base_ops/get_barcodes.py:317-318`: *"All barcodes are pre-generated
  at construction time and stored"*). I measured `get_barcodes(num_barcodes=24000, length=8)` at
  **0.043 s** with **24,000 strings** held in `op._barcode_strings` before a single library row
  exists.
- **Sequential `mutagenize` / `mutagenize_orf` build the complete list of position/substitution index
  tuples at construction** (`base_ops/mutagenize.py:332-381 _build_caches`; the count formula
  `comb(num_positions, num_mutations) * (alpha_size-1)**num_mutations` at 351-353). For the GB1
  double-mutant pool that is 536,085 tuples.
- **Source pools hold their sources:** `from_seqs`, `from_fasta` — but **not** `get_kmers`, which
  stores only `length`, alphabet size and the total count (`get_kmers.py:146`) and derives each
  k-mer from the state index on demand (`_value_to_kmer:198`, `_compute_core:212`); verified live,
  `get_kmers(length=8, mode='sequential')` holds 65,536 states and no k-mer collection.
- **`generate_library()` itself accumulates every row in a Python list** (`:108 rows = []`) before
  building the DataFrame; **`to_df` chunks but then `pd.concat`s** (`export_mixin.py:220`), so it
  also materialises. **Only `to_file` truly streams.** State the distinction explicitly: *the pool
  is lazy; a materialising call is not.*

**Positive flip side, previously unclaimed:** eagerness makes constraint infeasibility **fail fast,
at DAG-construction time, before any generation cost**. Live:
`get_barcodes(24000, length=8, gc_range=(0.4,0.6), min_edit_distance=3, max_homopolymer=2)` raises
`ValueError: Could only generate N barcodes satisfying constraints within 100000 attempts
(requested 24000)` immediately. *The call is unseeded, so N is not reproducible evidence — 419 when
first measured, 424 on re-execution. Quote the fail-fast behaviour, not the number.*

*Value stays `yes`:* the *library members* are never materialised, which is the row's question and
the paper's central mechanism. This is nonetheless the most attackable `yes` in Block A; the
fallback position if pressed is *"yes for library members; source pools and sequential index caches
are eager, and this is documented."*

*Source:* `generate_library.py:108,115,265`; `pool_mixins/export_mixin.py:220,222`;
`base_ops/get_barcodes.py:317-318`; `base_ops/mutagenize.py:332-381`; `examples/README.md`
(performance table); live execution.

### `library_algebra` — **yes**

All primitives but `sync` are first-class `Operation` subclasses, not user recipes:

| Primitive | File | Note |
|---|---|---|
| `stack` | `state_ops/stack.py` | disjoint union; delegates to `st.stack(parent_states)` at `:89`; `seq_length=None` when parents disagree (`:66-71`); card `active_parent` (`:55`) records which branch produced each row |
| `join` | `fixed_ops/join.py` | end-to-end concatenation (`spacer_str.join(seqs)`, `:89`), optional spacer (`:45`), length accounting (`:78-82`) |
| `sample` | `state_ops/sample.py` | `num_seqs`/`seq_states`/`seed`/`with_replacement` |
| `repeat` | `state_ops/repeat.py` | card `repeat_index` |
| `state_slice` | `state_ops/state_slice.py` | `pool[a:b]` |
| `state_shuffle` | `state_ops/state_shuffle.py` | |
| **`sync`** | `state_ops/sync.py` | **not an `Operation` subclass** — a module-level function that rewires the pools' states in place and returns `None` (`:26-68`). **See below — the under-evidenced differentiator** |

Operator sugar `pool_a + pool_b`, `pool * 3`, `pool[a:b]` verified live (`pool.py:169,175,185`).

**`sync()` deserves its own sentence and previously did not get one.** It couples state spaces so
coupled pools advance *together* rather than crossing, producing **1:1 pairing instead of a
Cartesian product**. It is what makes "give each of 24,000 CRE variants exactly one barcode" a
single operation rather than a manual zip, and it is now the **default** for `replace_region`
(`region_ops/replace_region.py:19-20 sync: bool = True`, a `[Unreleased]` **BREAKING** change from
`False`). Tools that lack it force the user to align two libraries by hand outside the tool.

*Source:* `src/poolparty/state_ops/`, `fixed_ops/join.py`, `region_ops/replace_region.py:19-20`,
`CHANGELOG.md [Unreleased]`, live execution.

### `exhaustive_single_scans` — **yes**

`mode="sequential"` is the exhaustive mode across the whole op set: `mutagenize`, `mutagenize_orf`,
`deletion_scan`, `insertion_scan`, `replacement_scan`, `mutagenize_scan`, `subseq_scan`,
`shuffle_scan`, `region_scan`, `get_kmers`, `from_iupac`.

Reproduced live: `pp.from_seq('ACGTACGT').mutagenize(num_mutations=1, mode='sequential')` →
`num_states=24` (8 × 3), names `mut_00`…`mut_23`, first rows `CCGTACGT`, `GCGTACGT`, `TCGTACGT`.
`deletion_scan(deletion_length=3)` on a 12-mer → 10 states. `mutagenize_orf` on a 5-codon ORF:
`missense_only_first` → 95 = 5 × 19; `nonsynonymous_first` → 100 = 5 × 20; `any_codon` → 315 = 5 × 63;
`nonsense` → 15 = 5 × 3.

**Scope precision (does not change the value, but write it this way).** At codon level "exhaustive"
means exhaustive **over amino acids, one representative codon each** — say *"all amino-acid
substitutions"*, not *"all codon substitutions"*. And the synonymous class **cannot be enumerated at
all**: `orf_ops/mutagenize_orf.py:205-207` refuses `mode='sequential'` for any non-uniform
`mutation_type`. Verified live — `synonymous`, `missense_only_random` and `nonsynonymous_random` all
raise *"mode='sequential' requires a uniform mutation type"*. This is a real limitation and is
recorded as one.

*Source:* `base_ops/mutagenize.py`, `orf_ops/mutagenize_orf.py:205-207`, `scan_ops/`, live execution.

### `sampled_random_mutagenesis` — **yes**

`mutagenize(mutation_rate=…, mode="random", num_states=N)`. The rate is a genuine **per-position
mutation probability**: `num_mut = rng.binomial(num_mutable, self.mutation_rate)`
(`base_ops/mutagenize.py:458`). Mutually exclusive with `num_mutations`
(`mutagenize.py:144-147`); `mode='sequential'` refused with `mutation_rate` (`:157-159`).
*(Citation drift corrected — the extraction memo cited "26-27, 149-155"; line 26 is a parameter
default, not a validation.)*

Codon-level equivalent: `mutagenize_orf(mutation_rate=0.1, mode="random", num_states=10000)`
(`docs/tutorials/dms_gb1.rst:126-141`).

**Reproducible by construction.** Every random operation is seeded deterministically per operation
per row: `np.random.SeedSequence([pool._master_seed, op.id, state_val])`
(`generate_library.py:283`, confirmed verbatim). The authors report byte-identical notebook output
across Python 3.10 and 3.12.

*Source:* `base_ops/mutagenize.py:144-147, 157-159, 458`; `generate_library.py:283`;
`docs/tutorials/dms_gb1.rst`.

### `higher_order_combinatorial` — **yes**

`base_ops/mutagenize.py:332-381 _build_caches` enumerates
`comb(num_positions, num_mutations) * (alpha_size-1)**num_mutations` (`:351-353`) via
`itertools.combinations`, with a genuine **non-uniform** branch (`:365-378`) for per-position IUPAC
masks where positions have unequal alternative counts.

Arithmetic re-derived, not taken on trust: 55 × 19 = **1,045** singles;
C(55,2) × 19² = 1,485 × 361 = **536,085** doubles. Both appear in `docs/tutorials/dms_gb1.rst`
(lines 91, 116) with the composition table totalling **547,230** at lines 263-275.
Orders beyond exhaustive range are covered by `mutation_rate` sampling.

*Source:* `base_ops/mutagenize.py:332-381`, `orf_ops/mutagenize_orf.py`, `docs/tutorials/dms_gb1.rst`.

### `heterogeneous_components_one_library` — **yes** *(best-evidenced row; likely held alone in the table)*

Structurally **and dimensionally** different components coexist in one library specification.
Reproduced end to end for this memo:

```
p1 = pp.from_seq('ACGTACGT').mutagenize(num_mutations=1, mode='sequential', prefix='mut')  #  8 bp, 24 states
p2 = pp.from_seq('T'*12).deletion_scan(deletion_length=3, mode='sequential', prefix='del')  # 12 bp, 10 states
p3 = pp.get_barcodes(num_barcodes=4, length=6, seed=2)                                      #  6 bp,  4 states
lib = pp.stack([p1, p2, p3])
# lib.num_states == 38 ; lib.seq_length is None
# lib.generate_library(num_seqs=38) -> 38 rows, realised lengths {6, 8, 12}
```

`state_ops/stack.py:66-71` explicitly sets `seq_length=None` when member lengths disagree — this is
a deliberate design decision, not an accident. `stack`'s `active_parent` design card
(`stack.py:55`) records **which component produced each row**, so heterogeneity is recoverable
downstream. `print_dag` renders the three-branch DAG.

At scale: the GB1 tutorial stacks four semantically different sub-libraries (exhaustive singles,
exhaustive doubles, rate-sampled higher-order, WT replicates) into one 547,230-member pool; the MPRA
tutorial combines a motif-placement branch with a barcode branch inside one construct.

*Source:* `state_ops/stack.py:55, 66-71, 89`; `docs/tutorials/{dms_gb1,mpra_regulatory_grammar}.rst`;
live execution.

### `combinatorial_motif_place` — **yes** *(evidence upgraded from asserted to demonstrated)*

`multiscan_ops/insertion_multiscan.py` (and `replacement_multiscan`) is a purpose-built combinatorial
motif placer: `insertion_pools: Sequence[Pool]` (several distinct motifs at once), `num_insertions`,
flat or per-insertion `positions`, `min_spacing` / `max_spacing`,
`insertion_mode: Literal["ordered","unordered"]`, and `replace=True` to hold total length constant.

**Demonstrated, not asserted.** Two distinct motif pools, four candidate positions, `min_spacing=6`:

| `insertion_mode` | states | interpretation |
|---|---|---|
| `'unordered'` | **12** | 6 spacing-legal position pairs × 2 pool→position permutations |
| `'ordered'` | **6** | 6 spacing-legal position pairs, fixed pool order |

Design cards returned `combination_index`, `starts` (`[0,10]`, `[0,20]`, `[0,30]`, `[10,20]`,
`[10,30]`, `[20,30]`) and `names`, which flip from `[_ins_0, _ins_1]` to `[_ins_1, _ins_0]` exactly
at the permutation boundary. So **spacing constraints, permutation enumeration and per-site
provenance are all directly demonstrable.** Orientation comes from `flip(mode='sequential')` → 2
states, verified `GGGCCCAAA` → `TTTGGGCCC`.

*Attribution nit, worth getting right:* those card keys are declared on `RegionMultiscanOp`
(`region_ops/region_multiscan.py:129`), which `insertion_multiscan` composes — the emitted column is
literally `op[3]:insertion_multiscan(region_multiscan).starts`.

At scale (`docs/tutorials/mpra_regulatory_grammar.rst:110-122`): 3 TFBSs × 1,000 sampled position
configurations × 2³ orientation combinations = **8,000 CRE states**. *Say "states", not "unique
CREs" as the tutorial does: the 1,000 configurations are drawn with replacement, so a few coincide —
re-running the documented pipeline gives 7,968-7,992 distinct sequences across seeds.*

*Source:* `multiscan_ops/insertion_multiscan.py`, `region_ops/region_multiscan.py:129`,
`base_ops/flip.py`, `docs/tutorials/mpra_regulatory_grammar.rst`, live execution.

### `barcode_generation` — **yes**

`base_ops/get_barcodes.py` (579 lines, genuinely vectorised: `_gc_filter_batch:40`,
`_homopolymer_filter_batch:47` with a cumsum sliding window, batched candidate draw at `:74`).
Constraints implemented and **enforced during generation, not post hoc**: `min_edit_distance`
(Levenshtein), `min_hamming_distance` (vectorised numpy against the accepted set), `max_homopolymer`,
`gc_range`, `avoid_sequences` + `avoid_min_distance`, variable `lengths` with `length_proportions`,
`seed`, and `max_attempts` with an explicit `_raise_insufficient` failure mode.

**Cite the operation's capability, not the tutorial's setting.** Exercised at a harder setting than
the shipped tutorial: `get_barcodes(num_barcodes=8, length=8, gc_range=(0.4,0.6),
min_edit_distance=3, max_homopolymer=2, seed=7)` → measured **minimum pairwise Levenshtein distance
exactly 3**, every GC exactly 0.50, no run > 2. **Presentational warning:** the shipped MPRA
tutorial uses `min_edit_distance=1` (`mpra_regulatory_grammar.rst:139`), which for distinct barcodes
is a no-op (distance ≥ 1 is just uniqueness). A referee will open the tutorial and find the distance
constraint doing nothing — do not lean on that line.

**Attachment is a supported operation, not a manual step.** `replace_region(region_name="bc",
content_pool=barcode_pool)` with the new default `sync=True` pairs each of 24,000 CRE variants with
exactly one barcode (`mpra_regulatory_grammar.rst:135-150`). Tests: `tests/test_get_barcodes.py`.

*Source:* `base_ops/get_barcodes.py:40,47,74,317-322`; `region_ops/replace_region.py:19-20`;
`docs/tutorials/mpra_regulatory_grammar.rst:135-150`; live execution.

### `per_sequence_provenance` — **yes**

**Design cards.** Each `Operation` subclass declares `design_card_keys`. I enumerated all
declarations myself: **29 sites, 17 of them non-empty.** The complete non-empty set:

| Op | Keys | File:line |
|---|---|---|
| `MutagenizeOp` | `positions`, `wt_chars`, `mut_chars` | `base_ops/mutagenize.py:113` |
| `MutagenizeOrfOp` | `codon_positions`, `wt_codons`, `mut_codons`, `wt_aas`, `mut_aas` | `orf_ops/mutagenize_orf.py:162` |
| `RegionMultiscanOp` | `combination_index`, `starts`, `ends`, `names`, `region_seqs` | `region_ops/region_multiscan.py:129` |
| `RegionScanOp` | `position_index`, `start`, `end`, `name`, `region_seq` | `region_ops/region_scan.py:98-104` |
| `FlipOp` | **`flip`** (not `strand` — extraction memo corrected) | `base_ops/flip.py:103` |
| `StackOp` | `active_parent` | `state_ops/stack.py:55` |
| `GetBarcodesOp` | `barcode_index`, `barcode` | `base_ops/get_barcodes.py:322` |
| `GetKmersOp` | `kmer_index`, `kmer` | `base_ops/get_kmers.py:107` |
| `RecombineOp` | `breakpoints`, `pool_assignments` | `base_ops/recombine.py:151` |
| `SeqShuffleOp` | `permutation` | `base_ops/shuffle_seq.py:89` |
| `FromSeqsOp` | `seq_name`, `seq_index` | `base_ops/from_seqs.py:103` (inherited by `from_fasta`) |
| `FromIupacOp` | `iupac_state` | `base_ops/from_iupac.py:89` |
| `FromMotifOp` | `prob_state` | `base_ops/from_motif.py:95` |
| `MaterializeOp` | `seq_index`, `seq_name` | `base_ops/materialize.py:25` |
| `RepeatOp` | `repeat_index` | `state_ops/repeat.py:59` |
| `FilterOp` | `passed` | `base_ops/filter_seq.py:21` |
| `ScoreOp` | the user's `card_key` | `fixed_ops/score.py:95` |

`generate_library.py:298-319` filters the raw card per the `cards` spec and writes one DataFrame
column per requested key, with either `op.name`-prefixed or user-supplied column names; a card from
an **inactive** DAG branch is written as `None` (`:316-319`). Two universal keys (`"seq"`, `"state"`)
work on every operation (`docs/metadata/design_cards.rst`).

**Caveat (kept — it is honest and cheap):** cards are **opt-in**. With no `cards=` argument the
DataFrame is exactly `['name', 'seq']` — confirmed live. Provenance is structured and
machine-readable, but the user must ask each operation for it.

*Source:* the 29 `design_card_keys` declarations above; `generate_library.py:298-319`;
`docs/metadata/design_cards.rst`; live execution.

### `automatic_naming` — **partial** *(honest downgrade; keep it — it pre-empts the obvious attack)*

Names **are** auto-composed: each operation contributes a segment
(`operation.py:413-449 compute_name_contributions`), segments are collected in topological order and
dot-joined (`generate_library.py:328-330`), with zero-padding sized to the state count. Verified
live: `m_00.d_0`, `m_00.d_1`, … for a chained `mutagenize`→`deletion_scan`.

Two reasons this is not `yes`:

1. **Naming is off unless the user opts in per operation.** `operation.py:436-437`:
   `if self.prefix is None: return []`. `docs/metadata/naming.rst:27-29` states it in the authors'
   own words: *"If no operation in the pipeline sets a prefix, the `name` column is `None`."*
   Verified live — `pp.from_seq('ACGT').mutagenize(num_mutations=1, mode='sequential')` gives
   `name=None` on every row, columns `['name','seq']`.
2. **Names are state indices, not variant descriptions.** `mut_00` identifies a variant but does not
   say *what* it is; there is no `p.Gln1Phe` or `chr:pos:ref>alt` informative name. Semantic content
   lives in design cards instead. Exceptions: `from_fasta` **in batch (list-of-coordinates) mode**
   names sequences `{chrom}:{start}-{stop}({strand})` (`fixed_ops/from_fasta.py:145-153`; a single
   coordinate tuple delegates to `from_seq` at `:109-125` and gets no such name), and
   `from_seqs(seq_names=[...])` accepts explicit names.

*Source:* `operation.py:413-449`, `generate_library.py:328-330`, `docs/metadata/naming.rst:27-29`,
live execution.

### `design_visualization` — **yes** *(text/ANSI only — keep the qualifier in the table cell)*

Both halves of the row's own wording ("graph view, highlighted sequences") exist, and both were run:

- **Graph view.** `pool.print_dag()` (`pool.py:385` → `text_viz/pool_op_tree.py`, `text_viz/graph.py`)
  renders an ASCII tree of the pool/operation DAG with per-node modes and state counts. The
  shared-node DAG above is real `print_dag` output.
- **Highlighted sequences.** `pool.print_library()` (`pool.py:282`) plus inline styling
  (`utils/style_utils.py`, `fixed_ops/stylize.py`, `stylize_orf`, per-op `style=`) renders mutations,
  deletions, region tags and inserted motifs in colour — ANSI in a terminal (`\x1b[91mC\x1b[0m` on a
  mutated base, observed live), with the HTML equivalent embedded in the docs (the MPRA output shows
  three TFBSs in blue/purple/orange plus a bold barcode). `docs/metadata/styling.rst`;
  `tests/test_inline_styles.py`.

**No graphical output.** `grep -riE 'matplotlib|pyplot|seaborn|plotly|logomaker'` over
`poolparty/src/**/*.py` and `poolparty/docs/**/*.rst` returns **0** hits (matches elsewhere in the
tree are producer metadata inside SVG assets — the untracked `docs/_build/html/` render **and** the
committed `docs/_static/images/figure4b_g.drawio.svg`, which carries a literal `Matplotlib v3.10.7`
string). No sequence logos, no rendered DAG image, no interactive viewer.

*Calibration:* `pydna` scores `yes` partly on ASCII `figure()` output; `ledidi` scores `yes` on plots
with no graph view. PoolParty has the graph view **and** the highlighted sequences — literally both
halves of the row. `yes` is consistent with the survey's own bar.

*Source:* `pool.py:282,385`; `text_viz/`; `utils/style_utils.py`; `docs/metadata/styling.rst`;
grep; live execution.

---

## BLOCK B — assay coverage

### `assay_dms` — **yes**

`orf_ops/mutagenize_orf.py` with `mutation_type` ∈ {`any_codon`, `nonsynonymous_first`,
`nonsynonymous_random`, `missense_only_first`, `missense_only_random`, `synonymous`, `nonsense`}
(`codon_table.py:16-24 VALID_MUTATION_TYPES`), `codon_positions` (list **or** `slice`), `frame` ∈
{±1,±2,±3}. `ProteinPool`, `translate` (style-preserving) and `reverse_translate` close the loop
between DNA and protein.

Shipped tutorial `docs/tutorials/dms_gb1.rst` — protein GB1, all numbers re-derived from the docs:
1,045 exhaustive singles + 536,085 exhaustive doubles + 10,000 rate-sampled higher-order + 100 WT
replicates = **547,230 library members**.

**Two precisions that must survive into the paper:**

1. **547,230 is a MEMBER count, not a distinct-sequence count.** With `mutation_rate=0.1` over 55
   codons, Binomial(55, 0.1) gives P(0 mut)=0.0030, P(1)=0.0186, P(2)=0.0558 — roughly **774 of the
   10,000 random draws are expected to duplicate the WT, single- or double-mutant arms.** Write
   *"547,230 library members"*; never *"547,230 distinct sequences"*.
2. **A synonymous-control arm cannot be exhaustively enumerated** (`mutagenize_orf.py:205-207`,
   verified live). This is directly relevant to DMS practice and is in the limitations.

*Source:* `orf_ops/mutagenize_orf.py`, `codon_table.py:16-24`, `docs/tutorials/dms_gb1.rst:91,116,263-275`.

### `assay_mpra` — **yes**

Shipped tutorial `docs/tutorials/mpra_regulatory_grammar.rst`, verified line by line: Melnikov-layout
construct (5′ adaptor / `<cre>` / KpnI+XbaI junction / `<bc>` / sequencing adapter) at `:33-35`;
`num_states=1000` sampled position configurations (`:110-115`) × 2³ TFBS orientation combinations =
8,000 CRE states (`:120`; not 8,000 *unique* CREs — see `combinatorial_motif_place`);
`repeat(times=3)` → 24,000 (`:113,122`); `get_barcodes` +
`replace_region` with default `sync=True` (`:135-150`); final `mpra_pool seq_length=153,
num_states=24000` (`:165`). Independently corroborated by the authors' own reproduction record in
`examples/README.md`. Backed by `insertion_multiscan`, `get_barcodes`, `replace_region`, `flip`,
`repeat`; design cards become the design factors for downstream regression.

*Caveat:* the tutorial's `min_edit_distance=1` is effectively a no-op — see `barcode_generation`.

*Source:* `docs/tutorials/mpra_regulatory_grammar.rst`, `examples/README.md`.

### `assay_insilico` — **partial** *(honest and correct)*

The mechanism exists: `fixed_ops/score.py` attaches any `Callable[[str], Any]` output as a
design-card covariate without altering the sequence, and `base_ops/filter_seq.py` gates on any
predicate. `docs/index.rst:9-10` advertises *"in silico analysis of genomic models"* on the landing
page. `docs/metadata/design_cards.rst` pitches cards as covariates for downstream models.

Why not `yes`:

- **No shipped example.** `docs/tutorials/index.rst` lists exactly two tutorials (`dms_gb1`,
  `mpra_regulatory_grammar`). `docs/_static/images/figure4a.drawio.svg` and `figure4b_g.drawio.svg`
  are committed and referenced by **zero** `.rst` files (grep confirmed) — a tutorial was staged and
  never written. The authors say so themselves (`examples/README.md:27-33`).
- **The only in-silico worked example is not distributed.** `examples/spliceai_surrogate.ipynb` is
  gitignored (`.gitignore:70`; `git ls-files` returns nothing). It covers **library construction
  only** — SpliceAI scoring and the GAM analysis "remain the standalone scripts in
  `poolparty/experiment/`", which is untracked and exists on one machine.
- **No batched inference path.** `ScoreOp._compute_core` (`fixed_ops/score.py:107-119`) calls `fn`
  once per sequence.

*Size figure, disambiguated:* the SpliceAI example builds **400,000 sequences across both arms**
(`examples/README.md`, "Stage C, both arms (400,000 seqs)", 29.8 s); the performance table's
"Library size 200,000" row refers to **one arm**. Quote whichever you mean, explicitly.

*Source:* `fixed_ops/score.py:107-119`, `base_ops/filter_seq.py`, `docs/index.rst:9-10`,
`docs/tutorials/index.rst`, `.gitignore:70`, `examples/README.md:27-33`, grep.

---

## BLOCK C — genomics integration

*This block is what buys the table its credibility. All five re-derived by me independently. My own
greps over `poolparty/src/poolparty/**/*.py`: `exon|intron|splice|transcript|gtf|gff|ensembl|refseq`
→ **0 hits**; `vcf|vep|hgvs` → **0 hits**.*

### `genome_coordinates` — **partial** *(do not let it be pushed to `yes`; do not concede it to `no`)*

`fixed_ops/from_fasta.py` (166 lines, read in full) is real genome-coordinate support at the **input**
end: `from_fasta(fasta_path, coordinates=(chrom, start, stop, strand))` or a list of such tuples,
indexed with `pyfaidx` (a declared runtime dependency). 0-based half-open; `strand='-'`
reverse-complements (`_extract_sequence:16-25`); `start > stop` handles circular-genome wrap-around
(`:27-29`); batch mode names sequences `{chrom}:{start}-{stop}({strand})` (`:145-153`); the extracted
sequence drops straight into a named region of a background pool. Tested
(`tests/test_from_fasta.py`), documented (`docs/operations/from_fasta.rst`).

What is absent, and why this is not `yes`:

- No BED / GTF / GFF / interval-file input — coordinates must be typed as Python tuples.
- Coordinates are **not** carried as structured metadata; they appear only inside the *name string*.
  `from_fasta` inherits `FromSeqsOp`'s card keys `seq_name`, `seq_index` (`from_seqs.py:103`) and
  declares none of its own.
- No coordinate tracking through downstream operations and **no way to map a designed variant back
  to a genomic position** — after an insertion or deletion the offsets are gone.
- No assembly/build awareness, no chromosome-name normalisation, no liftover.

*Cross-tool calibration (this is what makes `partial` defensible):* `tangermeme` = `yes` (BED + FASTA
via pyfaidx); `seqpro` = `partial` (interval I/O **without** sequence fetch); `biopython` = `partial`;
`dnachisel` = `no`. **PoolParty is seqpro's mirror image — fetch without interval files —** so
`partial` is exactly where it belongs.

*Source:* `fixed_ops/from_fasta.py`, `base_ops/from_seqs.py:103`, `pyproject.toml`,
`tests/test_from_fasta.py`.

### `transcript_models` — **no**

Nothing in the package reads or represents a transcript annotation.
`grep -riE 'exon|intron|splice|transcript|gtf|gff|ensembl|refseq'` over
`poolparty/src/poolparty/ --include=*.py` returns **0 hits** (run by me). The nearest structure,
`OrfRegion` (`region.py:69-95`), is `(name, seq_length, frame)` — a single contiguous span plus a
frame integer, not a transcript with exons, UTRs or a CDS model. No GTF/GFF parser exists.

*Source:* grep; `src/poolparty/region.py:69-95`.

### `exon_intron_split_codons` — **no**

`OrfRegion` is contiguous by construction, and `mutagenize_orf` resolves a **single** frame for a
**single** region (`orf_ops/mutagenize_orf.py:21-55 _resolve_frame`). There is no exon list, no
junction concept, and therefore no handling of a codon split across an exon boundary. Reading frames
±1/±2/±3 handle offset and orientation **within one contiguous ORF only**. Confirmed by the same
zero-hit grep.

*Source:* `orf_ops/mutagenize_orf.py:21-55`, `src/poolparty/region.py:69-95`, grep.

### `vcf_vep_output` — **no**

Export formats are exhaustively `Literal["csv", "fasta", "tsv", "jsonl"]` — declared on `to_file` at
`pool_mixins/export_mixin.py:225` and re-validated at `:351-353`, with gzip handled by extension.
`grep -riE 'vcf|vep|hgvs'` over `src/poolparty/` returns **0** substantive hits. No VCF writer, no
VEP/ANN field, no HGVS.

*Source:* `pool_mixins/export_mixin.py:225,351-353`; grep.

### `consequence_annotation` — **no** *(name the near-miss, then deny the capability — this paragraph is what a referee will read most carefully)*

The near-miss is real and must be stated precisely. `mutagenize_orf`'s `mutation_type`
(`synonymous`, `missense_only_*`, `nonsynonymous_*`, `nonsense`) **specifies which class of change to
generate**, and its design cards record `wt_codons` / `mut_codons` / `wt_aas` / `mut_aas`. That is
**bookkeeping of changes PoolParty itself introduced** into a user-supplied contiguous ORF.

It is **not** consequence annotation: there is no classification of arbitrary input variants, no
transcript context, no consequence vocabulary (`splice_donor`, `frameshift`, `stop_retained`, …), no
impact ranking, and nothing that would consume a VCF. Confirmed by grep — no classifier and no
consequence terms exist anywhere in the source.

*Source:* `orf_ops/mutagenize_orf.py:162`, `codon_table.py:16-24`, grep.

---

## BLOCK D — adjacent / complementary

### `primer_design` — **no**

`grep -riE 'primer|melting|anneal|hairpin'` over `poolparty/src/poolparty/ --include=*.py` returns
**exactly one hit, and it is a code comment**:
`data/restriction_enzymes.py:202  # Gibson assembly - common enzymes to avoid in Gibson primers`.
A lexical near-miss inside the restriction-enzyme *avoidance* presets, not a design capability. No
amplicon, overhang, Tm, or assembly-junction design.

*Source:* grep; `data/restriction_enzymes.py:202`.

### `codon_optimization` — **partial**

`orf_ops/reverse_translate.py` offers `codon_selection="first"` against a frequency-ordered table.
`codon_table.py:33-56 STANDARD_GENETIC_CODE` is the Kazusa *Homo sapiens* usage table with codons
sorted high→low frequency, and the source comment at `:26-32` documents that the ordering is load-
bearing: *"This ordering is important for mutation types like 'missense_only_first' and
'nonsynonymous_first' which select the first codon in each list (the most frequent when using the
built-in standard genetic code)."* So `"first"` = most-frequent-codon optimisation. Verified live:
`pp.reverse_translate('MKV', codon_selection='first')` → `ATGAAGGTG` (M→ATG, K→AAG, V→GTG — the head
of each list). A custom table can be supplied **per call** (`genetic_code: Union[str, dict]`), so
another host's preference order can be installed — but note that `pp.set_genetic_code` does *not*
reach `reverse_translate`, which builds its own `CodonTable` from its own argument (default
`"standard"`, `reverse_translate.py:143`); verified live.
`codon_selection="random"` samples synonymous codons uniformly.

Why not `yes`: the **entire strategy set is {take the first codon, pick uniformly at random}**. No
CAI or CAI-target optimisation, no GC-content targeting, no codon harmonisation, no rare-codon or
ramp handling, no avoidance of hairpins / repeats / restriction sites *during* codon choice, no
windowed or constrained optimisation. Constraint checking is a separate, after-the-fact `filter`.

*Minor note:* the live check works because `reverse_translate` string-promotes to a `ProteinPool`
(`reverse_translate.py:88-89`); passing a `DnaPool` raises `TypeError`.

*Source:* `orf_ops/reverse_translate.py:88-89`, `codon_table.py:26-56`, live execution.

### `synthesis_constraints` — **partial** *(justification rewritten — lead with missing constraint TYPES, not with "reject-only")*

Implemented and documented. `pool_mixins/filter_mixin.py` provides five filters:
`filter_gc` (`:20`), `filter_homopolymer` (`:65`), `filter_complexity` (linguistic complexity, `:106`),
`filter_dust` (NCBI DUST, `:151`), `filter_restriction_sites(enzymes=…, sites=…, check_rc=True)`
(`:194`) — backed by `data/restriction_enzymes.py`, a curated, IUPAC-aware, reverse-complement-aware
enzyme database with seven cloning presets (`golden_gate`, `common`, `mcs`, `gibson`,
`frequent_cutters`, `rare_cutters`, `blunt`). `get_barcodes` enforces GC range and homopolymer limits
at generation time. The underlying predicates are also **exported as standalone public API** —
`calc_gc`, `calc_complexity`, `calc_dust`, `has_homopolymer`, `has_restriction_site` are all on the
`poolparty` namespace (verified live) and usable outside the DAG. Live: a GC filter over
`from_iupac('NNN')` kept 32 of 64 and emitted 32 `NullSeq` rows with the acceptance-rate warning.

**Why not `yes` — the missing constraint types:** no *built-in* oligo-length or vendor-capacity limit
(a length check is expressible as a generic `filter` predicate, and `docs/operations/filter.rst:97-122`
demonstrates one), no Tm, no
secondary-structure / hairpin check, no pool-level repeat-content or cross-hybridisation analysis, no
background k-mer screen, and no constructive `split`/`pad` for oligos exceeding synthesis length.

*Explicitly NOT the reason:* "enforcement is reject-only". **Rejection is enforcement**, and the row
asks for "constraint checking/enforcement", not repair. Do not concede that framing — but do record
the consequence separately (a filter-heavy design burns states and can come up short, with an
explicit warning from `generate_library`), which belongs in limitations.

*Calibration:* `oligopoolcalc` earns `yes` on exactly the missing types plus constructive split/pad;
`mpradesign` and `mpranator` earn `partial` on **narrower** check sets than PoolParty's.
**PoolParty sits at the top of `partial`.**

*Source:* `pool_mixins/filter_mixin.py:20,65,106,151,194`; `data/restriction_enzymes.py`;
`utils/seq_properties.py`; live execution.

### `degenerate_iupac_codons` — **yes**

`base_ops/from_iupac.py` accepts any IUPAC string and either enumerates it exhaustively
(`mode="sequential"`) or samples it (`mode="random"`). Reproduced live:
`pp.from_iupac('NNK', mode='sequential')` → **`num_states=32`** (4 × 4 × 2), first rows
`AAG, AAT, ACG, ACT` — the standard NNK degenerate codon handled exactly. Additionally,
`mutagenize(allowed_chars=…)` takes a per-position IUPAC mask restricting which bases each position
may mutate to (`base_ops/mutagenize.py:161-185`), with correct **non-uniform** combinatorics in
`_build_caches` (`:365-378`); restriction-site matching also expands IUPAC codes.

*Caveat, worth keeping:* IUPAC is **expanded** into explicit sequences. PoolParty does not emit a
degenerate string (a literal `NNK` oligo order) as output — which matters because the row's peers
(CodonGenie, MPRAnator) do.

*Source:* `base_ops/from_iupac.py`, `base_ops/mutagenize.py:161-185, 365-378`, live execution.

### `negative_control_generation` — **yes** *(evidence rebuilt; the old demo was misleading)*

`base_ops/shuffle_seq.py` implements `shuffle_type="dinuc"` — an **Euler-path shuffle preserving
dinucleotide frequencies** (the standard scrambled-control construction; first and last characters
are fixed, a mathematical constraint of the algorithm, `shuffle_seq.py:40-43`) — as well as
`"mono"`. Region-scoped, card-recording (`permutation`, `:89`). `shuffle_scan` applies it
positionally. Tested: `tests/test_dinuc_shuffle.py`.

**Live, on a non-cyclic 24-mer `ACGGTTACCGATTGCAACGTTAGC`:** 5 requested dinuc shuffles → 5 **distinct**
sequences, every one with a dinucleotide-frequency vector identical to the source (Counter-checked)
and mononucleotide composition preserved; 3 distinct `mono` shuffles.

> **⚠ Do not reuse the extraction memo's demo.** It used `ACGT`×4, whose dinucleotide graph is a
> single Euler cycle — that sequence has **exactly one** valid dinuc shuffle and the operation
> returns the input unchanged. I reproduced this: all three "shuffles" came back as
> `ACGTACGTACGTACGT`. If that demo reached a referee it would read as a bug. It is not one; it is
> the correct output for a degenerate input.

**The row definition explicitly names "reverse/complement controls", and PoolParty has them as
first-class operations** — `rc` (fixed reverse complement; `GGGCCCAAA` → `TTTGGGCCC`) and `flip`
(`mode='sequential'` → 2 states, forward + reverse complement). Both verified live. WT/replicate
controls come from `repeat` (GB1 tutorial `wt_pool`, 100 copies); matched-variant control arms from
`replacement_scan`/`mutagenize` (the SpliceAI notebook builds a `GT`→`GA` disrupted twin for every
variant).

*The "no operation is NAMED control" hedge has been dropped* — it conceded more than the facts
require. Controls are assembled from first-class named operations (`shuffle_seq`, `rc`, `flip`,
`repeat`, `stack`), which is what the row asks for.

> **Framing note: this row is NOT a PoolParty exclusive.** `extractions/mpranator.md:169` records
> scramble/reverse/complement as first-class in MPRAnator, and `seqpro`'s Rust `k_shuffle` is the
> same construction. **Do not present it as a differentiator in the response letter.**

*Source:* `base_ops/shuffle_seq.py:40-43,89`; `base_ops/flip.py`; `rc`; `state_ops/repeat.py`;
`tests/test_dinuc_shuffle.py`; live execution.

### `ml_model_in_loop` — **partial**

`score(pool, fn, card_key=…)` and `filter(pool, predicate)` accept arbitrary Python callables, so a
model can annotate or gate a library from inside the DAG. Demonstrated in
`examples/spliceai_surrogate.ipynb`: `pp.from_motif(donor_pwm).filter(no_long_homopolymer)
.score(score_5ss, card_key="maxent_score", …)` with MaxEntScan preloaded.

Why not `yes` — `ScoreOp` is a **literal passthrough**. `fixed_ops/score.py:107-119`:

```python
def _compute_core(self, parents, rng=None):
    parent_seq = parents[0]
    ...
    result = self._fn(clean_seq)
    return parent_seq, {self._card_key: result}
```

The sequence is returned unchanged; the model's output influences only what is *recorded*, never
what is *generated*. `mode="fixed"`, `num_states=1`. `fn` is invoked **one sequence at a time** with
no batching and no tensor/GPU path (the notebook's own headline finding is that a naive per-call
model costs 3.8 ms/call). Only `filter` can act on a model, and only by rejection. There is no
gradient-based or iterative design loop, no active learning, no model-guided search — and no shipped
example in the distributed docs.

*Calibration:* `ledidi` and `tangermeme` earn `yes` on gradient/greedy optimisation. `partial` is the
correct distance from them.

*Source:* `fixed_ops/score.py:95,107-119`; `base_ops/filter_seq.py`; `examples/README.md`.

### `readout_analysis` — **no**

Design and export only. `grep -riE 'enrichment|count table|readout|fastq|bam'` over
`poolparty/src/poolparty/ --include=*.py` returns **7 hits, every one of them the string `BamHI`**
(`data/restriction_enzymes.py:28,169,195,205`; `pool_mixins/filter_mixin.py:217,249`;
`utils/seq_properties.py:319`) — lexical false positives from `bam`, nothing readout-related
(re-run by me; the previously recorded "0 hits" was wrong). Export is one-way
(`to_df`/`to_file` → csv/tsv/fasta/jsonl). `docs/metadata/design_cards.rst` explicitly frames design
cards as covariates for analysis performed **elsewhere**; in the authors' own SpliceAI example the
GAM analysis is an external, untracked script (`examples/README.md`).

*This row marks complementarity, not deficiency* — Oligopool Calculator's `yes` here is the reason
the row exists.

*Source:* grep; `pool_mixins/export_mixin.py`; `docs/metadata/design_cards.rst`; `examples/README.md`.

---

## BLOCK E — engineering and availability

> **Encoding fix (reviewer Finding 2, applied).** All four of these were the bare string `yes` in the
> extracted record, which answers nothing for a column asking "which interface?", "which license?",
> "what date?", "how many?". Below, each row gives the **exact table cell** to typeset — the survey
> convention already used in `final/seqpro.md` and `final/valiant.md`. (The structured record's
> `value` field is constrained to a yes/partial/no/unknown enum, so the qualified string is carried
> at the head of each evidence field. **Typeset the qualified string, not the enum.**)

### `interface` — table cell: **`yes (Python API only)`**

Python library / API only. `pyproject.toml` declares **no `[project.scripts]` and no entry points**
(grep confirmed: zero matches for `scripts`, `entry-points`, `gui`). No GUI, no web service, no R
package. Most operations are available both as `pp.op(pool, …)` and as `pool.op(…)`
(`__init__.py:331/373` factory maps; **39** operation factories with method forms — source and
multi-input operations such as `from_seq`, `from_fasta`, `get_barcodes`, `stack`, `join` and `sync`
are module functions only),
so pipelines read left-to-right. Runtime type checking via `beartype` on public factories (visible in
every traceback triggered during these checks). `requires-python = ">=3.10"`.

*Source:* `poolparty/pyproject.toml`, `src/poolparty/__init__.py:331,373,390,395`, live execution.

### `license` — table cell: **`yes (MIT)`**

Confirmed three ways: `poolparty/LICENSE` — "MIT License / Copyright (c) 2025 Justin B. Kinney"
(identical file at monorepo root); `poolparty/pyproject.toml:11 license = "MIT"`; PyPI
`info.license_expression = "MIT"` (**cite `license_expression` — the legacy `info.license` field is
`null`**). The companion `statetracker` package is also MIT.

*Source:* `LICENSE`, `poolparty/LICENSE`, `poolparty/pyproject.toml:11`, live PyPI JSON.

### `installable_today` — table cell: **`yes (pip, PyPI)`**

`pip install poolparty`. Live PyPI JSON (pulled 2026-08-10): version **0.1.1**, wheel
`2026-04-06T21:10:37` + sdist `21:10:39`; 0.1.0 at `21:03` the same day; `requires_python >=3.10`.
Runtime dependencies exactly `numpy>=1.20, pandas>=1.3, beartype>=0.22.9, statetracker>=0.1.0,
pyfaidx>=0.8.1, typing_extensions>=4.0`. **PoolParty itself ships no compiled extension, but `numpy`
and `pandas` do** — 63 `.so` files under `numpy`/`pandas` in this environment, e.g.
`numpy/_core/_multiarray_umath.cpython-312-x86_64-linux-gnu.so`. *(The earlier "all pure-Python, no
compiled extensions" claim was wrong and is falsifiable in one `find`.)*
CI (`.github/workflows/test.yml`, re-read): ubuntu-latest × {3.10, 3.11, 3.12}, plus macOS-latest and
windows-latest on 3.11, running the `statetracker` **and** `poolparty` suites (**2,906** `test_`
functions across **77** test files, both counts reproduced by me). Docs on Read the Docs
(`https://poolparty.readthedocs.io`). Import and execution verified in this environment.

*Source:* live PyPI JSON; `.github/workflows/test.yml`; `poolparty/pyproject.toml`; own counts.

### `last_activity` — table cell: **`2026-04-07 (commit 1bb0179); PyPI 0.1.1, 2026-04-06`**

Last repo commit `1bb0179` "fix RTD citation", **2026-04-07** (`git log -1`). PyPI 0.1.0 and 0.1.1
both uploaded **2026-04-06**. `CHANGELOG.md [Unreleased]` carries a **BREAKING** change not yet on
PyPI (`replace_region` now defaults to `sync=True`, `keep_tags=True`; previously both `False`), plus
the two tutorial pages and vectorised `get_barcodes` internals. Classifier:
`Development Status :: 3 - Alpha`. Local (gitignored) example notebooks were re-executed 2026-08-07.

*Packaging defect, stated because it is user-visible:* `src/poolparty/__init__.py:3` reads
`__version__ = "0.1.0"` while `pyproject.toml` and PyPI say `0.1.1`. In this environment
`importlib.metadata.version("poolparty")` **also** returns `0.1.0`, so the skew is present in the
editable install as well as in the module constant. A user who `pip install poolparty==0.1.1` and
prints `pp.__version__` sees the wrong number. Trivial to fix; fix it before resubmission.

*Source:* `git log`; live PyPI JSON; `CHANGELOG.md`; `src/poolparty/__init__.py:3`; live execution.

### `documented_examples` — table cell: **`3 worked examples (quickstart + 2 tutorials) + 50 operation reference pages`**

Counts re-derived by me (the extraction memo was wrong low on two of them):

| Item | Count | Note |
|---|---|---|
| Full worked examples shipped | **3** | quickstart, DMS GB1 tutorial, MPRA tutorial |
| `docs/operations/` `.rst` files | **62** | *not 57* |
| — of which per-operation reference pages | **50** | each with parameter table + runnable example + printed output |
| — of which concept/index pages | **12** | `modes`, `library_size`, `scanning`, `mutagenesis`, `region_operations`, `state_operations`, `source_operations`, `composition_operations`, `orf_operations`, `sequence_utilities`, `case_ops`, `index` |
| All `.rst` under `docs/` (excl. `_build`) | 74 | |
| Test files / test functions | **77** / **2,906** | *not 78 files* — the 78th directory entry is `__pycache__` |
| Executable notebooks | 3 | **gitignored, NOT distributed** |

**Shipped examples:**

| Name | What it is | Where |
|---|---|---|
| Quickstart | Step-by-step walkthrough of Pool, operations, regions, cards | `docs/quickstart.rst` |
| Deep Mutational Scanning: Protein GB1 | Full DMS design — 1,045 singles + 536,085 doubles + 10,000 rate-sampled + 100 WT reps = 547,230 members | `docs/tutorials/dms_gb1.rst` |
| MPRA Library for Regulatory Grammar | 3 TFBSs × 1,000 position configs × 8 orientations → 8,000 CREs × 3 barcoded replicates = 24,000 seqs, `seq_length=153` | `docs/tutorials/mpra_regulatory_grammar.rst` |
| Operation reference pages | 50 pages, parameter tables + runnable examples + printed output | `docs/operations/` |
| Concept pages | Modes, Library Size, Regions, Design Cards, Naming, Styling, Pool | `docs/operations/*.rst`, `docs/regions.rst`, `docs/metadata/*.rst`, `docs/pool.rst` |
| API reference | Autodoc of *selected* public functions/classes — **not complete**: omits `get_barcodes`, `translate`, `reverse_translate`, `annotate_orf`, `annotate_region`, `stylize_orf`, `set_genetic_code`, `set_progress_mode`, `fixed_operation`, `DnaPool`/`ProteinPool`, `to_df`/`to_file`, `clear_tags` and the sequence-property/restriction helpers (checked name by name) | `docs/api.rst` |
| README quick example | Branch + stack + design cards, with output | `README.md` |

**NOT shipped** (`.gitignore:70 /poolparty/examples/`; also `:69 /poolparty/notebooks/`,
`:68 /poolparty/benchmarks/`; `git ls-files` returns nothing for any of them):
`dms_gb1.ipynb`, `mpra_regulatory_grammar.ipynb`, `spliceai_surrogate.ipynb`.
The first two are generated verbatim from the shipped `.rst` tutorials, so their content *is*
distributed; the third — the only in-silico worked example — is not, and no docs tutorial replaces
it.

*Source:* own file counts; `.gitignore:68-70`; `git ls-files`; `docs/tutorials/index.rst`.

---

## Additional capabilities not covered by the row set

1. **Region tagging.** XML-style inline tags (`<cre>…</cre>`, `<bc/>`) travel with the sequence
   through operations, so any **region-aware** operation can target a named region **by name rather
   than by index** (`region.py`, `region_ops/`, `docs/regions.rst`). The state-level operations —
   `stack`, `sample`, `repeat`, `state_slice`, `state_shuffle`, `sync`, `join` — take no `region`
   argument (signatures checked live). Includes region-length validation, zero-length
   insertion points, and `insert_tags` / `remove_tags` / `extract_region` / `replace_region` /
   `annotate_region` / `annotate_orf` / `clear_annotation`.
2. **Addressable, random-access state space.** Backed by the companion `statetracker` package (MIT,
   on PyPI): every library index maps deterministically to exactly one sequence, so a library can be
   sliced, shuffled or sampled **without generating it**. (Not *split* — `split` exists only on the
   companion `statetracker` `State` API; PoolParty exposes neither `pp.split` nor `Pool.split`.)
3. **Library size known before generation.** See `library_first_class_object` for the exact semantics
   and caveats. The consequence not recorded there: the count is available at DAG-construction cost
   (0.03-0.24 s for libraries of 24k-547k), which is what makes synthesis-budget checking possible
   before committing to a design.
4. **`sync()` — 1:1 coupling instead of Cartesian product.** See `library_algebra`. Now the default
   for `replace_region`. A genuine differentiator against tools that force manual zipping.
5. **Shared DAG nodes with per-row memoisation.** See `composable_operations`.
6. **Fail-fast constraint infeasibility.** `get_barcodes` raises at **DAG-construction time** if the
   requested constraints cannot be met, before any generation cost is paid.
7. **Deterministic reproducibility.** Per-operation, per-row seeding via
   `np.random.SeedSequence([master_seed, op.id, state_val])`; the authors report byte-identical
   notebook output across Python 3.10 and 3.12.
8. **Constraint predicates as standalone public API.** See `synthesis_constraints` for the predicate
   list. Additionally `get_preset_enzymes`, `get_enzyme_site`, `ENZYME_PRESETS`, `ENZYME_SITES`.
9. **Recombination / chimeragenesis.** `recombine(sources=[…], num_breakpoints=k, positions=…)` with
   `breakpoints` / `pool_assignments` design cards. Requires two or more sources of equal length.
10. **k-mer tiling and sub-sequence scanning.** `get_kmers`, `subseq_scan`.
11. **Protein-level pools.** `ProteinPool`, `translate` (style-preserving), `reverse_translate`,
    `set_genetic_code`. There is no public protein-sequence source and no direct amino-acid
    substitution operation: protein variants are designed as codon-mutated DNA and then translated.
12. **Streaming export** to csv/tsv/fasta/jsonl with gzip, chunking, tqdm progress, and FASTA
    description templates/callables — on **`DnaPool` only** (`to_df`/`to_file` come from
    `ExportMixin`, which `ProteinPool` does not inherit; verified live).
13. **Performance toggles.** `toggle_styles`, `toggle_cards`, `set_progress_mode` — with the authors'
    own measurements of what each is worth (styles −26% on style-heavy DAGs; cards −2 to −6%).
14. **Engineering.** 2,906 tests across 77 files, `beartype` runtime typing, pre-commit config,
    multi-OS CI, `CITATION.cff` with ORCIDs, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, Read the Docs.
15. **User-defined operations.** New Operation types can be added by subclassing `Operation` and
    implementing a few required methods; custom Operations inherit state tracking, design-card
    generation and visualisation without reimplementing core infrastructure. The exported
    `fixed_operation` helper wraps a plain callable as a fixed-mode Operation. This is the
    extensibility the manuscript itself claims (`latex/main.tex:146-148,239`) and it is supporting
    evidence for `composable_operations`.

---

## Stated limitations

*(This section feeds the limitations paragraph Reviewer 3 asked for. It is deliberately longer, more
specific and more damaging than anything a referee is likely to produce unaided. Every item is
checkable at file:line or was measured. Do not soften it — its length is the point.)*

**Scope: genomics integration**

1. **No genome or transcript integration beyond FASTA extraction.** No GTF/GFF/BED input, no
   transcript or exon models, no exon-spanning codons, no VCF/VEP output, no consequence annotation
   of arbitrary input variants. `from_fasta` extracts by `(chrom, start, stop, strand)` tuples typed
   in Python, and **the coordinates survive only inside the sequence name string** — designed
   variants cannot be mapped back to genomic positions after an insertion or deletion.
   *(`fixed_ops/from_fasta.py`; grep over `src/poolparty/` for
   `exon|intron|splice|transcript|gtf|gff|ensembl|refseq` and `vcf|vep|hgvs` returns 0 hits.)*

**Scope: mutagenesis**

2. **Synonymous variants cannot be exhaustively enumerated.** `orf_ops/mutagenize_orf.py:205-207`
   refuses `mode='sequential'` for any non-uniform `mutation_type`; verified live, `synonymous`,
   `missense_only_random` and `nonsynonymous_random` all raise *"mode='sequential' requires a uniform
   mutation type"*. **A DMS user building a synonymous-control arm hits this immediately** — the
   workaround is random sampling, which does not guarantee coverage.
3. **Codon-level "exhaustive" means exhaustive over amino acids, not over codons.**
   `missense_only_first` enumerates 19 substitutions per codon using **one representative
   (most-frequent) codon** per amino acid, not all synonymous DNA variants.
4. **Deletion scans emit gap characters by default.** `deletion_scan` on a
   12-mer with `deletion_length=3` returns 12-character strings such as `TTTTTTTT---T`. Passing
   `deletion_marker=None` produces genuinely shorter sequences instead
   (`scan_ops/deletion_scan.py:100-107`, with a worked "True deletion" example at
   `docs/operations/deletion_scan.rst:111-134`); `clear_gaps` is a third route. Anyone expecting
   shorter oligos *by default* for fixed-length synthesis will be surprised.
5. **Undocumented mode restrictions.** `mode='sequential'` is refused with `mutation_rate`
   (`base_ops/mutagenize.py:157-159`); `flip` refuses `mode='fixed'` (`base_ops/flip.py:122-125`);
   `shuffle_seq` supports random mode only; and **a single sequential operation is hard-capped at
   1,000,000 states** (`Operation.max_num_sequential_states`, `operation.py:28-55`) — verified live,
   `get_kmers(length=11, mode='sequential')` raises `ValueError: Number of states (4194304) exceeds
   max_num_sequential_states (1000000). Use mode='random' instead.` rather than enumerating 4^11,
   and the ceiling appears nowhere in the `.rst` docs. None is fatal; all are ten-minute discoveries.
6. **Dinucleotide shuffling fixes the first and last characters** (a mathematical constraint of the
   Euler-path algorithm), and a sequence whose dinucleotide graph is a single cycle (e.g. a perfect
   tandem repeat) has exactly one valid shuffle — the operation correctly returns the input
   unchanged, which can look like a failure.

**Metadata and naming**

7. **Naming is opt-in and non-descriptive.** With no `prefix` anywhere in the pipeline the `name`
   column is `None` (`operation.py:436-437`; `docs/metadata/naming.rst:27-29`, the authors' own
   words). Where present, names are **state indices** (`mut_00.del_3`), not variant descriptions —
   there is no `p.Gln1Phe` or `chr:pos:ref>alt` form.
8. **Design cards are opt-in too.** Without an explicit `cards=` argument the output DataFrame is
   exactly `['name','seq']`. Provenance is structured and machine-readable, but the user must ask
   each operation for it.

**Design constraints and adjacent tooling**

9. **Constraints reject; they never repair.** A violating sequence becomes `NullSeq` and is dropped
   (`discard_null_seqs=True`). This is **neither silent nor an automatic shrink**: with `num_seqs=N`
   `generate_library` advances to further states until N valid rows are collected, and warns
   explicitly (`min_acceptance_rate`, `max_iterations`, state-space exhaustion) when it comes up
   short (`generate_library.py:115-184`; verified live — a 64-state pool filtered to 8 survivors
   returns 8 rows for `num_seqs=8` and warns for `num_seqs=20`). With `num_cycles` the output does
   shrink to the valid subset. Nothing is repaired, re-drawn *under* the constraint, or optimised.
10. **Missing synthesis-constraint types.** No *built-in* oligo-length or vendor-capacity limit — a
    length check is expressible as a generic `filter` predicate and the docs demonstrate one
    (`docs/operations/filter.rst:97-122`) — no Tm, no
    secondary-structure/hairpin check, no pool-level repeat-content or cross-hybridisation analysis,
    no background k-mer screen, no constructive split/pad for over-length oligos.
11. **Codon optimisation is rudimentary.** Strategy set is exactly {most-frequent codon, uniform
    random synonymous}. No CAI or CAI target, no GC targeting, no harmonisation, no rare-codon/ramp
    handling, no constraint-aware codon selection.
12. **No primer, adapter, or assembly design.** The only lexical hit for `primer` in the entire
    source is a comment in the restriction-enzyme avoidance presets.
13. **No readout or analysis side.** PoolParty designs and exports; it does not process sequencing
    counts, compute enrichment, or fit models. Design cards are explicitly framed as covariates for
    analysis performed elsewhere — in the authors' own SpliceAI example, by an external untracked
    script.
14. **Models are observers, not drivers.** `ScoreOp` is a literal passthrough
    (`fixed_ops/score.py:107-119`): the model's output is recorded, never fed back into generation.
    Calls are **one sequence at a time** with no batching or tensor/GPU path; there is no
    gradient-based or iterative optimisation loop and no active learning.
15. **Visualisation is text-only.** ASCII DAG tree and ANSI/HTML sequence styling. No plots, no
    sequence logos, no rendered graphics, no interactive viewer — zero plotting dependencies by
    design.

**Performance and implementation**

16. **`stack` evaluates dead branches — a measured 3.03x self-inflicted penalty.** `_compute_one`
    walks the **entire** topologically sorted operation list for **every** row, not just the active
    branch; only the *design-card write* is gated on `op.state.is_active`
    (`generate_library.py:316-319`). On the GB1 library this costs **170.32 s stacked against 56.20 s
    for the sum of its parts (+114.1 s, 3.03x)**; cProfile at 50,000 sequences shows the stacked pool
    calling the inner `mutagenize_orf` 150,000 times (3/seq) against 50,000 (1/seq) bare. The
    penalty reproduces at similar magnitude on Python 3.10 and 3.12, so it is a property of `stack`,
    not of the interpreter. **The authors' own pre-publication audit lists branch pruning (item H3,
    est. ~70% speedup) as unimplemented**, alongside style suppression (H2, ~34%).
17. **Measured throughput is 2,400-13,300 seq/s** (MPRA 2,372; DMS 3,213; SpliceAI 13,263) with peak
    memory +19 to +332 MB. Generation is 1.24-1.62x slower on Python 3.10 than 3.12, and the penalty
    is **operation-dependent, not a flat tax** (`mutagenize_orf` degrades far more than
    `replacement_scan`) — always quote the interpreter alongside a generation number.
18. **Laziness has exceptions, and they are real.** `get_barcodes` pre-generates and stores **every**
    barcode at DAG-construction time (the authors' own docstring, `get_barcodes.py:317-318`;
    measured 24,000 strings in 0.043 s). Sequential `mutagenize`/`mutagenize_orf` materialise the
    complete list of position/substitution index tuples at construction (536,085 for GB1 doubles).
    `from_seqs`/`from_fasta` hold their sources (`get_kmers` does **not** — it derives each k-mer
    from the state index on demand). `generate_library()` accumulates all
    rows in a list, and `to_df` chunks but then `pd.concat`s — **only `to_file` truly streams.**
    *The pool is lazy; a materialising call is not.*
19. **Global `Party` context.** Pools cannot exist outside a `Party`, and `pp.init()` mutates
    module-level state (`pool.py:48`, `__init__.py:269`). A default Party is created on import so the
    error is rarely seen, but the implicit global registry complicates concurrent and
    library-embedded use.

**Maturity and distribution**

20. **Version 0.1.1, `Development Status :: 3 - Alpha`, first released 2026-04-06**, with a
    **BREAKING** change already queued in the unreleased changelog (`replace_region` sync/keep_tags
    defaults). `src/poolparty/__init__.py:3` still reports `__version__ = "0.1.0"`.
21. **The in-silico application advertised on the docs landing page has no shipped tutorial.**
    `docs/index.rst:9-10` promises "in silico analysis of genomic models"; `docs/tutorials/index.rst`
    ships two tutorials, neither of them in-silico; `figure4a.drawio.svg` and `figure4b_g.drawio.svg`
    are committed and referenced by nothing.
22. **All three executable notebooks are gitignored and therefore not distributed**
    (`.gitignore:70`). The DMS and MPRA notebooks are regenerated verbatim from the shipped `.rst`
    tutorials so their content survives; the SpliceAI notebook — the only in-silico worked example —
    does not ship at all.
23. **Docs tutorials are not covered by CI.** `docs/conf.py` loads no `nbsphinx`/`myst`, so tutorial
    outputs are pasted by hand and `test.yml` runs only pytest. The undistributed notebooks are
    currently the only check that the tutorials still execute.
24. **`from_motif` is random-mode only.** A PWM source pool cannot be enumerated in sequential mode,
    so a motif-derived library has no exhaustive traversal and no exact enumeration of the sequences
    a PWM can produce above a threshold.
25. **`materialize` severs the DAG.** Materialising a pool replaces the upstream graph with a fixed
    sequence collection, so downstream operations lose access to the construction history and to the
    state-space arithmetic that makes size known before generation.

---

## Availability

- **Status:** installable today. `pip install poolparty` → 0.1.1 (wheel + sdist, PyPI, 2026-04-06).
- **Python:** `>=3.10`; CI on 3.10/3.11/3.12 (ubuntu) + 3.11 (macOS, Windows).
- **Dependencies:** `numpy>=1.20`, `pandas>=1.3`, `beartype>=0.22.9`, `statetracker>=0.1.0`,
  `pyfaidx>=0.8.1`, `typing_extensions>=4.0`. PoolParty itself has no compiled extension; `numpy`
  and `pandas` do.
- **License:** MIT (both `poolparty` and `statetracker`).
- **Docs:** https://poolparty.readthedocs.io
- **Repository:** https://github.com/jbkinney/poolparty-statetracker
- **Maturity:** `Development Status :: 3 - Alpha`; 2,906 tests across 77 files.

---

## Adjudicated disagreements (extractor vs. reviewer)

**No capability VALUE was disputed.** All 33 values agree between the extraction memo, the
adversarial review, and this final record. The disagreements were factual details, and I resolved
each one against the source myself. **No row is set to `unknown`** — none of the disputes was
genuinely unresolvable.

| Item | Extractor | Reviewer | **My finding** | Winner |
|---|---|---|---|---|
| `_topo_sort_operations` span | 220-237 | 219-236 | **220-237** (`def` on 220, `return result` on 237) | **extractor** |
| `seq_cache` line | — | 245 | **260** (populated 292, consumed 273) | neither; corrected |
| card `is_active` gate | 298-319 (block) | 314 | block **298-319**; gate at **316-319** | **extractor** |
| `mutagenize` mutual exclusivity | 26-27, 149-155 | 144-147 / 157-159 | **144-147** and **157-159** | **reviewer** |
| `docs/operations/` `.rst` files | 57 | 62 | **62** (50 operation + 12 concept/index) | **reviewer** |
| `tests/` `.py` files | 78 | 77 | **77** `.py` (78 = directory entries incl. `__pycache__`); 2,906 `test_` fns exact | **reviewer** |
| `FlipOp` design card key | `strand` | `flip` | **`flip`** (`flip.py:103`) | **reviewer** |
| Multiscan card keys owner | `InsertionMultiscanOp` | `RegionMultiscanOp` | **`RegionMultiscanOp`** (`region_multiscan.py:129`); column reads `op[N]:insertion_multiscan(region_multiscan).starts` | **reviewer** |
| `design_card_keys` sites | not counted | "20 sites" | **29 declarations, 17 non-empty** (full table above) | neither; corrected |
| `codon_table` mutation-type set | 16-25 | 16-25 | **16-24** (`VALID_MUTATION_TYPES`); Kazusa table 33-56, its comment 26-32 | corrected |
| `stack.py` `seq_length=None` | 68-73 | 66-71 | **66-71** | **reviewer** |
| SpliceAI library size | "400,000 seqs" | "defensible; disambiguate" | **400,000 across both arms**; the perf table's "200,000" is **one arm** | both, disambiguated |
| Party caveat severity | "pools cannot exist outside a Party" | "less damaging — default Party on import" | **reviewer is right**: `pp.get_active_party() is not None` → `True` at import (`__init__.py:269`) | **reviewer** |
| `matplotlib` grep | "0 hits" | "0 hits" | **0 hits in `src/**/*.py` and `docs/**/*.rst`**; matches elsewhere are SVG producer metadata — in `docs/_build/html/` **and** in the committed `docs/_static/images/figure4b_g.drawio.svg` (`Matplotlib v3.10.7`) — scope the grep or the claim is falsifiable | corrected |
| `synthesis_constraints` justification | "reject-only" | "missing constraint types" | **reviewer** — rejection *is* enforcement; the row does not ask for repair | **reviewer** |
| `negative_control_generation` demo | `ACGT`×4 | "bad demo — single Euler cycle" | **reviewer** — reproduced: all three "shuffles" returned the input unchanged | **reviewer** |
| `lazy_generation` evidence | "Nothing materialises the library" | overreach | **reviewer** — `get_barcodes` docstring contradicts it; caveat rewritten and extended to `to_df`/`generate_library` | **reviewer** |

## Housekeeping items to fix before resubmission

*(Not capability findings. Listed because each is visible to anyone who clones the repository.)*

1. **`__version__` skew** — `src/poolparty/__init__.py:3` says `0.1.0`; PyPI/`pyproject.toml` say
   `0.1.1`. `importlib.metadata.version()` also returns `0.1.0` in the editable install.
2. **Repo-name skew** — the local working copy is `poolparty-statecounter/` and contains **both** a
   `statetracker/` and a `statecounter/` directory, while PyPI `project_urls` and `CITATION.cff`
   give the canonical repository as `github.com/jbkinney/poolparty-statetracker`. **Cite the
   published URL in the manuscript, never the local directory name.**
3. **`src/poolparty/marker_ops/` contains nothing but `__pycache__`** — a dead package directory.
4. **Path corrections for the sources table** — `CITATION.cff`, `CONTRIBUTING.md`,
   `CODE_OF_CONDUCT.md`, `.pre-commit-config.yaml` and `.readthedocs.yaml` live at the **monorepo
   root**, not inside `poolparty/`. (`LICENSE` and `.readthedocs.yaml` exist in both places.)
5. **`CITATION.cff` still says `version: "0.1.0"`, `date-released: "2026-04-03"`** — update to 0.1.1
   / 2026-04-06, and replace the `https://doi.org/XXXX` placeholder in `CHANGELOG.md`.
6. **MPRA tutorial's `min_edit_distance=1` is a no-op.** Either raise it to a value that constrains
   (≥ 2) or add a sentence saying it is set to 1 deliberately for uniqueness only.
7. **Docs tutorials are not in CI.** Consider adding an nbsphinx/myst execution check so the shipped
   tutorials cannot silently rot.
