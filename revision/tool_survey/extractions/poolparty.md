# PoolParty — capability extraction memo

**Slug:** `poolparty`
**Assessed:** 2026-08-10
**Version assessed:** PyPI `poolparty` 0.1.1 (uploaded 2026-04-06); repo working tree at commit `1bb0179` (2026-04-07).
**Note on bias:** this is the authors' own tool. Every "yes" below is anchored to implementing
code or an executed check, not to a docs sentence. Where a capability exists only as a manual
recipe or only as an unshipped notebook, it is scored `partial`. Block C is scored `no`
throughout except one honest `partial`.

## Sources consulted

| Kind | Path / ref |
|---|---|
| source | `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/src/poolparty/` (102 `.py` files) |
| readme | `poolparty/README.md` |
| docs | `poolparty/docs/index.rst`, `quickstart.rst`, `pool.rst`, `regions.rst`, `metadata/{design_cards,naming,styling}.rst`, `operations/` (57 `.rst` pages), `tutorials/{dms_gb1,mpra_regulatory_grammar}.rst` |
| changelog | `poolparty/CHANGELOG.md` |
| examples | `poolparty/examples/{README.md,dms_gb1.ipynb,mpra_regulatory_grammar.ipynb,spliceai_surrogate.ipynb}` (**gitignored — not distributed**, see `.gitignore:70 /poolparty/examples/`) |
| tests | `poolparty/tests/` — 78 files, 2,906 `test_` functions |
| other | `poolparty/pyproject.toml`, `LICENSE`, `.github/workflows/test.yml`, `CITATION.cff` |
| other | companion package `statetracker/` (state algebra engine) |
| pypi | `https://pypi.org/pypi/poolparty/json` — 0.1.0 and 0.1.1, both 2026-04-06 |
| other | live checks run against `poolparty-statecounter/.venv/bin/python` (read-only, no writes outside `tool_survey/`) |

---

## BLOCK A — library specification

### `library_first_class_object` — **yes**

`Pool` is a returned, holdable object. `pool.py:33-247` gives it identity (`_id`, `name`,
`named()`), inspectable state (`num_states`, `seq_length`, `parents`, `regions`,
`has_region()`), value semantics (`copy()`, `deepcopy()` — the latter recursively copies the
whole upstream DAG), and operator overloads (`__add__` = stack, `__mul__` = repeat,
`__getitem__` = state slice).

Pools are passed onward as *arguments to other operations*, not just chained:
`insertion_multiscan(..., insertion_pools=[hnf4a, ppara, xbp1])`
(`multiscan_ops/insertion_multiscan.py:15`), `replace_region(..., content_pool=barcode_pool)`,
`replacement_scan(cryptic, ...)`.

Caveat (recorded, not disqualifying): pools cannot exist outside a `Party` context —
`pool.py:45-49` raises `RuntimeError("Pools must be created inside a Party context")`. A
module-level default Party is created on import (`__init__.py: _init_default_party()`), and
`pp.init()` resets it. So there is an implicit global registry behind the object.

### `composable_operations` — **yes**

Every operation is `Pool(s) -> Pool` and is exposed both as a module function and as a Pool
method (`__init__.py:_POOL_FACTORY_MAP` / `_DNAPOOL_FACTORY_MAP` copy factory docstrings onto
methods). The result is an arbitrary DAG, topologically sorted at generation time
(`generate_library.py:220-237 _topo_sort_operations`).

Nesting is genuine, not just left-to-right chaining:
`region_ops/apply_at_region.py:10-34` takes `transform_fn: Callable[[Pool], Pool]` and applies
an arbitrary sub-pipeline to the content of a named region, e.g.
`pp.apply_at_region(bg, 'region', lambda p: pp.mutagenize(p, num_mutations=1))`.

Verified live: a three-branch DAG printed by `print_dag` shows `stack` over `mutagenize`,
`deletion_scan(region_scan)+replace_region`, and `get_barcodes` subtrees.

### `lazy_generation` — **yes** (with one honest caveat)

`generate_library.py:115-186` computes **one row at a time**: it assigns
`pool.state.value = global_state % pool.state.num_values` (line 265), which propagates through
the state DAG, then walks `sorted_ops` calling `op.compute(...)` per row. Nothing materialises
the library. `ExportMixin.to_file` (`pool_mixins/export_mixin.py:222-246`) streams in
`chunk_size` blocks straight to csv/tsv/fasta/jsonl (+`.gz`) "to avoid loading the entire
library into memory".

Measured (authors' own numbers, `examples/README.md:104-114`): DAG construction 0.03–0.24 s vs
10–170 s to generate the same libraries — a 500–3000x gap.

**Caveat:** sequential-mode `mutagenize` / `mutagenize_orf` build a *complete Python list of
index combinations* at construction time (`base_ops/mutagenize.py:_build_caches`, lines
331-380: `for positions in combinations(range(num_positions), self.num_mutations)`). For the
GB1 double-mutant pool that is 536,085 tuples held in memory before a single sequence exists.
Sequences are lazy; the sequential index enumeration for these two ops is not.

### `library_algebra` — **yes**

All are first-class operations in `state_ops/`:
- `stack` (`state_ops/stack.py`) — disjoint union of state spaces, delegating to
  `statetracker.stack`; `seq_length=None` when members differ.
- `join` (`fixed_ops/join.py`) — end-to-end concatenation of pools with optional spacer.
- `sample(pool, num_seqs=..., seq_states=..., seed=..., with_replacement=True)`.
- `repeat(pool, times=...)`; `state_slice` (`pool[0:100]`); `state_shuffle`; `sync(pools)`
  (couples state spaces so pools advance together instead of crossing).
- Operator sugar: `pool_a + pool_b`, `pool * 3`, `pool[a:b]` (`pool.py:169-189`).

### `exhaustive_single_scans` — **yes**

`mode="sequential"` is the exhaustive mode across the op set:
`mutagenize(num_mutations=1, mode="sequential")` (all 3L substitutions — README example yields
9 for a 3 bp region), `mutagenize_orf` (all 19 missense per codon — GB1 tutorial: 1,045 =
55x19), `deletion_scan`, `insertion_scan`, `replacement_scan`, `mutagenize_scan`,
`subseq_scan`, `shuffle_scan`, `region_scan`, `get_kmers`.

Verified live: `pp.from_seq('ACGTACGT').mutagenize(num_mutations=1, mode='sequential')` →
`num_states=24` (8x3), names `mut_00`…`mut_23`.

### `sampled_random_mutagenesis` — **yes**

`mutagenize(mutation_rate=..., mode="random", num_states=N)` —
`base_ops/mutagenize.py:26-27, 149-155`; rate is a per-position mutation probability
(binomial), explicitly mutually exclusive with `num_mutations`.
Codon-level equivalent: `mutagenize_orf(mutation_rate=0.1, mode="random", num_states=10000)`
(GB1 tutorial, "Random higher-order mutants"). Every random op is seeded deterministically per
operation via `np.random.SeedSequence([master_seed, op.id, state_val])`
(`generate_library.py:283`), so random libraries are reproducible.

### `higher_order_combinatorial` — **yes**

`base_ops/mutagenize.py:_build_caches` enumerates
`comb(num_positions, num_mutations) * (alpha_size-1)**num_mutations` via
`itertools.combinations` — exhaustive pairwise and k-wise. At codon level,
`orf_ops/mutagenize_orf.py` with `num_mutations=2` gives C(55,2)x19² = **536,085** double
mutants (GB1 tutorial, verified by the authors' notebook run). Higher orders beyond exhaustive
range are covered by `mutation_rate` sampling.

### `heterogeneous_components_one_library` — **yes**

Verified live (not just documented). Stacking three structurally *and dimensionally* different
pools into one library:

```
a = pp.from_seq('ACGTACGT').mutagenize(num_mutations=1, mode='sequential', prefix='mut')  #  8 bp
b = pp.from_seq('TTTTTTTTTTTT').deletion_scan(deletion_length=3, mode='sequential', ...)  # 12 bp
c = pp.get_barcodes(num_barcodes=4, length=6, ...)                                        #  6 bp
lib = pp.stack([a, b, c])   # -> seq_length=None, num_states=38, generates all 38 rows
```

`state_ops/stack.py:68-73` explicitly sets `seq_length=None` when members disagree. The GB1
tutorial stacks four semantically different sub-libraries (singles, doubles, rate-sampled
higher-order, WT replicates) into one 547,230-member pool; the MPRA tutorial combines a
motif-placement branch with a barcode branch in one construct.

### `combinatorial_motif_place` — **yes**

`multiscan_ops/insertion_multiscan.py` (and `replacement_multiscan`) is a purpose-built
combinatorial motif placer:
- `insertion_pools: Sequence[Pool]` — several distinct motifs at once;
- `positions` — flat or per-insertion position sets; `min_spacing` / `max_spacing`;
- `insertion_mode: Literal["ordered","unordered"]` — "'unordered' uses all permutations";
- `replace=True` keeps total length constant;
- orientation via `flip(mode="sequential")`, which enumerates forward + reverse complement.

MPRA tutorial (`docs/tutorials/mpra_regulatory_grammar.rst`): 3 TFBSs x 1,000 sampled position
configurations x 2³ orientation combinations = 8,000 unique CREs, with design cards recording
`positions`, spatial order (`tfbs`), and per-site strand.

### `barcode_generation` — **yes**

`base_ops/get_barcodes.py` (579 lines, vectorised). Constraints implemented and enforced at
generation, not post-hoc: `min_edit_distance` (Levenshtein, `_edit_distance`),
`min_hamming_distance` (vectorised numpy against the accepted set), `max_homopolymer`
(`_homopolymer_filter_batch`, sliding-window cumsum), `gc_range` (`_gc_filter_batch`),
`avoid_sequences` + `avoid_min_distance`, variable `lengths` with `length_proportions`, `seed`,
`max_attempts` with an explicit `_raise_insufficient` failure mode.

Attachment is a supported operation, not a manual step: MPRA tutorial does
`cre_pool.replace_region(region_name="bc", content_pool=barcode_pool)` with default `sync=True`,
pairing each of 24,000 CRE variants with exactly one barcode. Tests: `tests/test_get_barcodes.py`.

### `per_sequence_provenance` — **yes**

Design cards. Each `Operation` subclass declares `design_card_keys` (e.g. `MutagenizeOp`:
`["positions","wt_chars","mut_chars"]`; `MutagenizeOrfOp`:
`["codon_positions","wt_codons","mut_codons","wt_aas","mut_aas"]`; `InsertionMultiscanOp`:
`["combination_index","starts","ends","names","region_seqs"]`; `FlipOp`: strand; `FilterOp`:
`passed`; `ScoreOp`: the user's `card_key`). `generate_library.py:298-319` filters the raw card
per the `cards` spec and writes one DataFrame column per requested key, with either
`op.name`-prefixed or user-supplied column names. Two universal keys (`"seq"`, `"state"`) work
on every operation (`docs/metadata/design_cards.rst`).

**Caveat:** cards are **opt-in**. Without a `cards=` argument the output is `name` + `seq` only
(`generate_library` docstring, lines 53-55). Provenance is structured and machine-readable, but
the user must ask each operation for it.

### `automatic_naming` — **partial**

Names are auto-composed: each operation contributes a segment
(`operation.py:413-449 compute_name_contributions`), segments are collected in topological
order and dot-joined (`generate_library.py:295, 329`), with zero-padding sized to the state
count. Verified live: `mut_00.del_0`, `mut_00.del_1`, … for a chained mutagenize→deletion_scan.

Two reasons this is not `yes`:
1. **Naming is off unless the user opts in per operation.** `operation.py:436-437`:
   `if self.prefix is None: return []`. `docs/metadata/naming.rst:26-29` states plainly: "If no
   operation in the pipeline sets a prefix, the `name` column is `None`." Verified live —
   `pp.from_seq('ACGT').mutagenize(num_mutations=1, mode='sequential')` produces `name=None` for
   every row.
2. **Names are state indices, not descriptions.** `mut_00` identifies a variant but does not say
   *what* it is; there is no `p.Gln1Phe` / `chr:pos:ref>alt` style informative name. The
   semantic content lives in design cards instead. Exceptions: `from_fasta` names encode
   `{chrom}:{start}-{stop}({strand})`, and `from_seqs(seq_names=[...])` accepts explicit names.

### `design_visualization` — **yes** (text/ANSI, not graphical)

Two implemented, documented mechanisms:
- `pool.print_dag()` (`pool.py:385-390` → `text_viz/pool_op_tree.py`, `text_viz/graph.py`) —
  ASCII tree of the pool/operation DAG with modes and state counts. Verified live.
- Inline styling (`utils/style_utils.py`, `fixed_ops/stylize.py`, `stylize_orf`, per-op `style=`
  parameter) — `print_library` renders mutations, deletions, region tags and inserted motifs in
  colour (ANSI in a terminal; the docs embed the HTML equivalent, e.g. the MPRA output showing
  three TFBSs in blue/purple/orange plus a bold barcode).
  `docs/metadata/styling.rst`; `tests/test_inline_styles.py`.

Recorded honestly: there is **no** graphical output — no matplotlib/plotting dependency (grep
for `plot`/`matplotlib`/`logo` across `src/` and `docs/*.rst` returns 0 hits), no sequence
logos, no rendered DAG image, no interactive viewer.

---

## BLOCK B — assay coverage

### `assay_dms` — **yes**
`orf_ops/mutagenize_orf.py` with `mutation_type` ∈ {`any_codon`, `nonsynonymous_first`,
`nonsynonymous_random`, `missense_only_first`, `missense_only_random`, `synonymous`, `nonsense`}
(`codon_table.py:16-25`), `codon_positions` (list or `slice`), `frame` ∈ {±1,±2,±3}.
`ProteinPool`, `translate`, `reverse_translate` complete the loop.
Shipped tutorial: `docs/tutorials/dms_gb1.rst` — GB1, 1,045 singles + 536,085 doubles + 10,000
rate-sampled higher-order + 100 WT replicates = 547,230.

### `assay_mpra` — **yes**
Shipped tutorial `docs/tutorials/mpra_regulatory_grammar.rst`: Melnikov-layout construct
(5' adaptor / `<cre>` / KpnI+XbaI junction / `<bc>` / sequencing adapter), TFBS placement,
constrained barcodes, replicate structure, 24,000 sequences, design cards as design factors.
Backed by `insertion_multiscan`, `get_barcodes`, `replace_region`, `flip`, `repeat`.

### `assay_insilico` — **partial**
The mechanism exists: `fixed_ops/score.py` attaches any `Callable[[str], Any]` output as a
design-card covariate without altering the sequence, and `base_ops/filter_seq.py` gates on any
predicate. `docs/index.rst:10` advertises "in silico analysis of genomic models" and
`docs/metadata/design_cards.rst:35-38` pitches cards as "covariates in regression models".

Why not `yes`:
- **No shipped example.** `docs/tutorials/` contains only DMS and MPRA.
  `examples/README.md:27-33` (authors' own note): "**No SpliceAI docs tutorial.**
  `docs/_static/images/figure4a.drawio.svg` and `figure4b_g.drawio.svg` are committed but
  referenced by nothing — a tutorial was staged and never written."
- The only in-silico worked example, `examples/spliceai_surrogate.ipynb`, is **gitignored**
  (`.gitignore:70`) and therefore not distributed; it also covers library construction only —
  "SpliceAI scoring and the GAM analysis remain the standalone scripts in
  `poolparty/experiment/`" (untracked, on one machine).
- `score` calls `fn` once per sequence (`ScoreOp._compute_core`), so there is no batched
  inference path for a neural model.

---

## BLOCK C — genomics integration

### `genome_coordinates` — **partial** (do not read this as more than it is)
`fixed_ops/from_fasta.py` is real genome-coordinate support at the *input* end:
`from_fasta(fasta_path, coordinates=(chrom, start, stop, strand))` or a list of such tuples,
indexed with `pyfaidx` (a declared runtime dependency, `pyproject.toml`). 0-based half-open;
`strand='-'` reverse-complements (`_extract_sequence`); `start > stop` handles circular-genome
wrap-around; batch mode names sequences `{chrom}:{start}-{stop}({strand})`; extracted sequence
can be dropped into a named region of a background pool. Tested (`tests/test_from_fasta.py`),
documented with worked examples (`docs/operations/from_fasta.rst`).

What is **absent**, and why this is not `yes`:
- No BED / GTF / GFF / interval-file input — coordinates must be typed as Python tuples.
- Coordinates are **not** carried as structured metadata. They appear only inside the *name
  string*; `from_fasta`'s design-card keys are `'seq_name'`, `'seq_index'` only.
- No coordinate tracking through downstream operations, and **no way to map a designed variant
  back to a genomic position** — after an insertion or deletion the offsets are gone.
- No assembly/build awareness, no chromosome-name normalisation, no liftover.

### `transcript_models` — **no**
Nothing in the package reads or represents a transcript annotation. A `grep -riE
"exon|intron|splice|transcript|gtf|gff|ensembl|refseq"` over `src/poolparty/` returns **zero**
hits (the only `chrom` hits are `from_fasta`'s local variables). The nearest structure,
`OrfRegion` (`region.py:69-95`), is a single contiguous span plus a `frame` integer — not a
transcript with exons, UTRs, or a CDS model.

### `exon_intron_split_codons` — **no**
`OrfRegion` is contiguous by construction, and `mutagenize_orf` resolves a single `frame` for a
single region (`orf_ops/mutagenize_orf.py:21-55 _resolve_frame`). There is no exon list, no
junction concept, and therefore no handling of a codon split across an exon boundary. Reading
frames ±1/±2/±3 handle offset and orientation within one contiguous ORF only.

### `vcf_vep_output` — **no**
Export formats are exhaustively `{csv, tsv, fasta, jsonl}` (+ `.gz`) —
`pool_mixins/export_mixin.py:_detect_file_type` mapping, and the `file_type` Literal on
`to_file`. No VCF writer, no VEP/ANN field, no HGVS.

### `consequence_annotation` — **no**
Near-miss worth stating precisely, because a referee will look here. `mutagenize_orf`'s
`mutation_type` (`synonymous`, `missense_only_*`, `nonsense`, `nonsynonymous_*`) *specifies
which class of change to generate*, and its design cards record `wt_codons`/`mut_codons`/
`wt_aas`/`mut_aas`. That is bookkeeping of changes PoolParty itself introduced into a
user-supplied contiguous ORF. It is **not** consequence annotation: there is no classification
of arbitrary input variants, no transcript context, no consequence vocabulary (splice_donor,
frameshift, stop_retained, …), no impact ranking, and nothing that would consume a VCF.

---

## BLOCK D — adjacent, complementary

### `primer_design` — **no**
Zero hits for `primer`, `melting`, `Tm`, `anneal` across `src/poolparty/`. No amplicon,
overhang, or assembly-junction design. (The only near-miss is the `"gibson"` /
`"golden_gate"` entries in the restriction-enzyme *avoidance* presets.)

### `codon_optimization` — **partial**
`orf_ops/reverse_translate.py` offers `codon_selection="first"` against a frequency-ordered
table — `codon_table.py:26-56` is the Kazusa *Homo sapiens* usage table with codons sorted
high→low frequency, so "first" = most-frequent-codon optimisation. Verified live:
`pp.reverse_translate('MKV', codon_selection='first')` → `ATGAAGGTG`. A custom table can be
supplied (`genetic_code: Union[str, dict]`, also `pp.set_genetic_code`), so a different host's
preference order can be installed. `codon_selection="random"` samples synonymous codons.

Why not `yes`: the only strategies are "take the first codon" and "pick uniformly at random".
There is no CAI or CAI-target optimisation, no GC-content targeting, no codon harmonisation, no
rare-codon or ramp handling, no avoidance of hairpins / repeats / restriction sites *during*
codon choice, no windowed or constrained optimisation. Constraint checking is a separate,
after-the-fact `filter`.

### `synthesis_constraints` — **partial**
Implemented and documented: `pool_mixins/filter_mixin.py` provides `filter_gc`,
`filter_homopolymer`, `filter_complexity` (linguistic complexity), `filter_dust` (NCBI DUST),
and `filter_restriction_sites(enzymes=..., sites=..., check_rc=True)` backed by
`data/restriction_enzymes.py` (a curated enzyme database with presets `golden_gate`, `common`,
`mcs`, `gibson`, `frequent_cutters`, `rare_cutters`, `blunt`; IUPAC-aware site matching).
`get_barcodes` enforces GC range and homopolymer limits at generation time.

Why not `yes`: these operations only **reject** (a failing sequence becomes `NullSeq` and is
dropped with `discard_null_seqs=True`) — nothing repairs or re-designs a violating sequence, so
a filter-heavy design can silently shrink (hence `min_acceptance_rate` warnings in
`generate_library`). There are also no oligo-length or vendor-capacity limits, no Tm or
secondary-structure/hairpin checks, no pool-level repeat-content or cross-hybridisation
analysis.

### `degenerate_iupac_codons` — **yes**
`base_ops/from_iupac.py` accepts any IUPAC string and either enumerates it exhaustively
(`mode="sequential"`) or samples it (`mode="random"`). Verified live:
`pp.from_iupac('NNK', mode='sequential')` → `num_states=32` (4x4x2), first rows `AAG, AAT, ACG,
ACT`. Additionally, `mutagenize(allowed_chars=...)` takes a per-position IUPAC mask that
restricts which bases each position may mutate to (`base_ops/mutagenize.py:161-185`), and
restriction-site matching expands IUPAC codes.

Caveat recorded: IUPAC is *expanded* into explicit sequences. PoolParty does not emit a
degenerate string (e.g. a literal `NNK` oligo order) as output.

### `negative_control_generation` — **yes**
`base_ops/shuffle_seq.py` implements `shuffle_type="dinuc"` — an Euler-path shuffle preserving
dinucleotide frequencies (the standard scrambled-control construction) as well as `"mono"`;
region-scoped and card-recording (`permutation`). Tested: `tests/test_dinuc_shuffle.py`.
Verified live on a 16-mer, 3 distinct dinuc-preserving shuffles. `shuffle_scan` applies it
positionally. WT/replicate controls come from `repeat` (GB1 tutorial `wt_pool`, 100 copies) and
matched-variant control arms from `replacement_scan`/`mutagenize` (SpliceAI notebook builds a
`GT`→`GA` disrupted twin for every variant).
Caveat: there is no operation *named* "control"; controls are assembled from
`shuffle_seq`/`repeat`/`stack`.

### `ml_model_in_loop` — **partial**
`score(pool, fn, card_key=...)` and `filter(pool, predicate)` accept arbitrary Python
callables, so a model can annotate or gate a library from inside the DAG. Demonstrated in
`examples/spliceai_surrogate.ipynb`: `pp.from_motif(donor_pwm).filter(no_long_homopolymer)
.score(score_5ss, card_key="maxent_score", ...)` with MaxEntScan preloaded.

Why not `yes`: `ScoreOp` is a **passthrough** — the model's output never influences what is
generated, only what is recorded; `fn` is invoked one sequence at a time with no batching or
tensor/GPU path (the notebook's own headline finding is that a naive per-call model costs
3.8 ms/call); there is no gradient-based or iterative design loop, no active learning, no
model-guided search, and no shipped example in the distributed docs.

### `readout_analysis` — **no**
Design and export only. Zero hits for `enrichment`, `count table`, `readout`, `fastq`, `bam`
across `src/` and the `.rst` docs. Design cards are explicitly framed as covariates for
analysis performed *elsewhere* (`docs/metadata/design_cards.rst`); in the authors' own SpliceAI
example the GAM analysis is an external, untracked script.

---

## BLOCK E — engineering and availability

- **interface** — Python library / API only. No CLI (`pyproject.toml` declares no
  `[project.scripts]` / entry points), no GUI, no web service. Every operation is available both
  as `pp.op(pool, ...)` and as `pool.op(...)` (`__init__.py` factory maps), so pipelines read
  left-to-right. Runtime type checking via `beartype`.
- **license** — MIT. `poolparty/LICENSE` ("MIT License / Copyright (c) 2025 Justin B. Kinney"),
  `pyproject.toml: license = "MIT"`, confirmed in PyPI metadata.
- **installable_today** — Yes. `pip install poolparty`; PyPI has 0.1.1 as sdist + wheel
  (`requires_python >=3.10`). Runtime deps: numpy>=1.20, pandas>=1.3, beartype>=0.22.9,
  statetracker>=0.1.0, pyfaidx>=0.8.1, typing_extensions>=4.0. CI (`.github/workflows/test.yml`)
  runs 2,906 tests on Python 3.10/3.11/3.12 on ubuntu, plus macOS and Windows on 3.11. Docs on
  Read the Docs. Verified importable and functional in this environment.
- **last_activity** — PyPI 0.1.0 and 0.1.1 both uploaded **2026-04-06**; last repo commit
  `1bb0179` "fix RTD citation" **2026-04-07**. `CHANGELOG.md` carries an `[Unreleased]` section
  with a BREAKING change (`replace_region` defaults) not yet on PyPI. Local (uncommitted,
  gitignored) example notebooks were re-executed 2026-08-07. Classifier:
  `Development Status :: 3 - Alpha`.
  Minor packaging nit noticed: `src/poolparty/__init__.py:3` still reads `__version__ = "0.1.0"`
  while `pyproject.toml` and PyPI say 0.1.1.
- **documented_examples** — see the list below.

### Documented examples shipped

| Name | What it is | Where |
|---|---|---|
| Quickstart | Step-by-step walkthrough of Pool, operations, regions, cards | `docs/quickstart.rst` |
| Deep Mutational Scanning: Protein GB1 | Full DMS design — 1,045 singles, 536,085 doubles, 10,000 rate-sampled, 100 WT reps; 547,230 total | `docs/tutorials/dms_gb1.rst` |
| MPRA Library for Regulatory Grammar | 3 TFBSs x positions x orientations in a 100 bp CRE, barcoded, 24,000 seqs | `docs/tutorials/mpra_regulatory_grammar.rst` |
| Operation reference pages | 57 `.rst` pages; ~50 per-operation pages each with parameter table + runnable examples + printed output | `docs/operations/` |
| Concept pages | Modes, Library Size, Sequence Regions, Design Cards, Sequence Names, Styling, Pool | `docs/operations/{modes,library_size}.rst`, `docs/regions.rst`, `docs/metadata/*.rst`, `docs/pool.rst` |
| API reference | Autodoc of all public functions/classes | `docs/api.rst` |
| README quick example | Branch + stack + design cards, with output | `README.md` |
| `dms_gb1.ipynb` | Executable notebook, generated verbatim from the DMS tutorial | `examples/` — **gitignored, not distributed** |
| `mpra_regulatory_grammar.ipynb` | Executable notebook, generated verbatim from the MPRA tutorial | `examples/` — **gitignored, not distributed** |
| `spliceai_surrogate.ipynb` | In-silico surrogate-modelling design, 400,000 seqs, MaxEntScan via `score` | `examples/` — **gitignored, not distributed; no docs tutorial exists** |

---

## Additional capabilities not covered by the key list

- **Region tagging.** XML-style inline tags (`<cre>...</cre>`, `<bc/>`) travel with the sequence
  through operations, so any op can target a named region by name rather than by index
  (`region.py`, `region_ops/`, `docs/regions.rst`). Region-length validation, zero-length
  insertion points, `insert_tags`/`remove_tags`/`extract_region`/`replace_region`.
- **Addressable, random-access state space.** Backed by the companion `statetracker` package
  (also MIT, on PyPI): every library index maps deterministically to one sequence, so a library
  can be sliced, shuffled, sampled or split without generating it.
- **Library size known before generation.** `pool.num_states` and `pool.seq_length` are exact
  and available at DAG-construction cost (0.03–0.24 s for libraries of 24k–547k), which is what
  makes synthesis-budget checks possible before committing.
- **Deterministic reproducibility.** Per-operation seeding via
  `np.random.SeedSequence([master_seed, op.id, state_val])`; the authors report byte-identical
  notebook output across Python 3.10 and 3.12.
- **Recombination / chimeragenesis.** `recombine(sources=[...], num_breakpoints=k, positions=...)`.
- **k-mer tiling.** `get_kmers`, `subseq_scan`.
- **Protein-level pools.** `ProteinPool`, `translate` (style-preserving), `reverse_translate`.
- **Streaming export** to csv/tsv/fasta/jsonl with gzip, chunking, tqdm progress, FASTA
  description templates/callables.
- **Performance toggles.** `toggle_styles`, `toggle_cards`, `set_progress_mode`.
- **Engineering.** 2,906 tests, beartype runtime typing, pre-commit config, multi-OS CI,
  `CITATION.cff`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`.

---

## Stated limitations (write these into the paper; an all-green column will not survive review)

1. **No genome/transcript integration beyond FASTA extraction.** No GTF/GFF/BED input, no
   transcript or exon models, no exon-spanning codons, no VCF/VEP output, no consequence
   annotation. `from_fasta` extracts by `(chrom, start, stop, strand)` but the coordinates are
   preserved only in the sequence *name* — designed variants cannot be mapped back to genomic
   positions.
2. **Naming is opt-in and non-descriptive.** With no `prefix` anywhere in the pipeline the
   `name` column is `None`; where present, names are state indices (`mut_00.del_3`), not
   variant descriptions. Provenance lives in design cards, which are themselves opt-in.
3. **Constraints reject; they never repair.** All synthesis-constraint handling is
   filter-based — a violating sequence becomes null and is discarded. No sequence is fixed,
   re-drawn under constraint, or optimised. No oligo-length/vendor limits, Tm, secondary
   structure, or pool-level cross-hybridisation checks.
4. **Codon optimisation is rudimentary** — most-frequent-codon or uniform-random synonymous
   choice only; no CAI/GC targets, harmonisation, or constraint-aware codon selection.
5. **No primer, adapter, or assembly design.**
6. **No readout/analysis side.** PoolParty designs and exports; it does not process sequencing
   counts, enrichment, or model fitting.
7. **Models are observers, not drivers.** `score`/`filter` accept arbitrary callables, but the
   model output cannot steer generation; calls are per-sequence with no batching, and there is
   no optimisation loop.
8. **Visualisation is text-only** — ASCII DAG tree and ANSI/HTML sequence styling; no plots,
   logos, or rendered graphics.
9. **Performance.** Measured throughput 2,400–13,300 seq/s. `stack` is the dominant cost
   because `_compute_one` evaluates *every* branch of the DAG for *every* row, not just the
   active branch — a measured **3.0x** penalty on the GB1 library (170 s stacked vs 56 s as the
   sum of its parts). The authors' own audit lists branch pruning as an unimplemented ~70%
   speedup.
10. **Sequential mutagenesis is not fully lazy.** `mutagenize`/`mutagenize_orf` in sequential
    mode materialise the complete list of position/substitution index tuples at pool
    construction (536,085 entries for the GB1 double-mutant pool).
11. **Global Party context.** Pools cannot exist outside a `Party`; `pp.init()` mutates
    module-level state, which complicates use in concurrent or library-embedded settings.
12. **Maturity.** Version 0.1.1, `Development Status :: 3 - Alpha`, first release 2026-04-06,
    with a BREAKING change already queued in the unreleased changelog. The in-silico
    application advertised on the docs landing page has no shipped tutorial, and the three
    executable notebooks are gitignored and therefore not distributed to users.
