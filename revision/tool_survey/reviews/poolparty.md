# Adversarial review — PoolParty capability record

**Reviewed:** 2026-08-10
**Record under review:** `revision/tool_survey/extractions/poolparty.md` + the structured record
**Reviewer stance:** hostile. PoolParty is the authors' own tool; the dominant failure mode is
self-flattery, and the second is false modesty that concedes ground the paper does not have to concede.
**Method:** independent re-derivation from `poolparty/src/poolparty/`, `docs/`, `examples/`,
`.gitignore`, `pyproject.toml`, `.github/workflows/test.yml`, the live PyPI JSON API, and
**16 read-only executions** against `poolparty-statecounter/.venv/bin/python` run with
`PYTHONDONTWRITEBYTECODE=1` from `/tmp`. No file was written outside
`revision/tool_survey/`. The repo was not modified.

---

## Headline

**The record survives the attack.** 26 of 30 capability values are supported as written; none is
`wrong`; none is a manual recipe dressed up as a feature. Every load-bearing Block A "yes" was
reproduced by independent execution, not merely re-read. Block C is clean: my own greps returned
**zero** hits for `exon|intron|splice|transcript|gtf|gff|ensembl|refseq` and zero for
`vcf|vep|hgvs` across `src/poolparty/`, exactly as claimed. The self-flattery risk did not
materialise — if anything the record is harsh on itself in two places.

**Four findings require action before typesetting:**

1. **`lazy_generation`'s evidence overreaches.** "Nothing materialises the library" is false as an
   absolute. `get_barcodes` pre-generates and stores **every** barcode at DAG-construction time —
   the class docstring says so outright (`base_ops/get_barcodes.py:317-318`: *"All barcodes are
   pre-generated at construction time and stored"*). I measured `pp.get_barcodes(num_barcodes=24000,
   length=8, ...)` at **0.07 s and 24,000 stored strings before a single library row exists**. The
   value `yes` is still correct — the *library* is not materialised — but the caveat currently names
   only `mutagenize`/`mutagenize_orf` and must name `get_barcodes` (and, for completeness,
   `from_seqs`/`from_fasta`/`get_kmers`, which hold their source sequences). A referee who greps
   `get_barcodes.py` and finds the authors' own docstring contradicting their own "nothing
   materialises" sentence will use it.
2. **Four Block E values are uninformative as encoded.** `interface`, `license`, `last_activity` and
   `documented_examples` are all recorded as the bare string `"yes"`. `last_activity` asks for a date
   and `documented_examples` asks for a count; `"yes"` answers neither. This is the *same* defect the
   VaLiAnT review already flagged (`final/valiant.md` §8, "column-semantics problem"), and the
   established survey convention is the qualified form used in `final/seqpro.md` and
   `final/valiant.md`: `yes (Python API only)`, `yes (MIT)`, `yes (AGPL-3.0-or-later)`. Re-encode.
3. **A real, checkable limitation is missing: synonymous variants cannot be exhaustively enumerated.**
   `orf_ops/mutagenize_orf.py:205-207` refuses `mode="sequential"` for any non-uniform
   `mutation_type`. Verified live — `synonymous`, `missense_only_random` and `nonsynonymous_random`
   all raise *"mode='sequential' requires a uniform mutation type"*, while `missense_only_first`,
   `any_codon` and `nonsense` enumerate fine. A DMS referee building a synonymous-control arm hits
   this immediately. It belongs in `stated_limitations`; it does **not** downgrade
   `exhaustive_single_scans`.
4. **The record's strongest available evidence for `composable_operations` is not in it.** The record
   demonstrates *nesting* (`apply_at_region`) and *branching* (a three-branch `stack`), but never
   demonstrates a **shared DAG node** — the property that separates a DAG from a tree, and the one a
   referee will actually test. I verified it: a single `mutagenize` pool feeding both an `rc` branch
   and a second `mutagenize` branch, stacked, gives `12 + 144 = 156` states with `pool[8]` appearing
   once under each branch in `print_dag` and computed once per row (`generate_library.py:245`
   `seq_cache`). Put this in the record. It is free and it is the row's whole point.

---

## 1. Block A — the rows the paper's argument rests on

### `library_first_class_object` = yes — **supported**

Re-read `pool.py` in full. Everything claimed is at the cited location: `num_states` (135),
`parents` (141), `seq_length` (147), `regions`/`has_region` (152-156), `named()` (197),
`copy()` (201), `deepcopy()` (224, recursively copies the upstream DAG via `Operation.deepcopy`),
`__add__`→stack (168), `__mul__`/`__rmul__`→repeat (175/181), `__getitem__`→state_slice (186).
Pools genuinely travel as arguments, not just as chain receivers — confirmed in
`insertion_multiscan(insertion_pools=[...])` (`multiscan_ops/insertion_multiscan.py:15`),
`replace_region(content_pool=...)` (`region_ops/replace_region.py:16`), and my own live runs.

The `Party` caveat is stated honestly and is real (`pool.py:45-49`). One refinement: the record says
"pools cannot exist outside a Party context"; in practice a default Party is created on import
(`_init_default_party`), so the failure mode is only hit after an explicit `pp.init()` or inside a
`with pp.Party()` block that has exited. That makes the caveat *less* damaging than written — do not
over-concede it.

### `composable_operations` = yes — **supported (evidence should be strengthened)**

`_topo_sort_operations` is at `generate_library.py:219-236` (record says 220-237; off by one, cosmetic).
Verified live, beyond what the record claims:

- **Two-operation nested sub-pipeline** through `apply_at_region`:
  `apply_at_region(bg, 'r', lambda p: deletion_scan(mutagenize(p, ...), ...))` → 168 states,
  names `m_00.d_0`, `m_00.d_1`, … So `transform_fn` takes an arbitrary sub-DAG, not one op.
- **Shared node / true DAG** (see Finding 4 above) — 156 states, `pool[8]` reused.
- `apply_at_region`'s *annotation* is bare `Callable`; `Callable[[Pool], Pool]` is the docstring.
  Quote the docstring, not a type annotation, if this is cited in the paper.

Composability restrictions that exist and are not stated anywhere in the record:
`mode="sequential"` is refused with `mutation_rate` (`base_ops/mutagenize.py:157-159`), with
non-uniform `mutation_type` (`orf_ops/mutagenize_orf.py:205-207`, Finding 3), and `flip` refuses
`mode="fixed"` (`base_ops/flip.py:122-125`). None is fatal; all are the kind of thing a referee
finds in ten minutes.

### `lazy_generation` = yes — **supported, evidence must be amended** (Finding 1)

The core mechanism is exactly as described and I read it line by line.
`generate_library.py:115` `while len(rows) < num_seqs:` → `_compute_one` → line 265
`pool.state.value = global_state % pool.state.num_values`, which propagates through the state DAG,
then a single pass over `sorted_ops` calling `op.compute()`. Nothing enumerates the library.
`ExportMixin.to_file` (`pool_mixins/export_mixin.py:222`) streams in `chunk_size` blocks —
docstring: *"Generates sequences in chunks and writes them incrementally to avoid loading the entire
library into memory."* The authors' performance table (`examples/README.md`) is real and I read it:
DAG construction 0.03/0.05/0.24 s against 15.08/10.12/170.32 s generation.

Two required corrections to the *evidence*, not the value:

- `get_barcodes` eager pre-generation (Finding 1). Note this also has a **positive** side the record
  misses: constraint infeasibility fails fast, at DAG-construction time, before any generation cost.
  Live: `get_barcodes(num_barcodes=24000, length=8, gc_range=(0.4,0.6), min_edit_distance=3,
  max_homopolymer=2)` raises `ValueError: Could only generate 11538 barcodes satisfying constraints
  within 100000 attempts` immediately. That is good engineering and is worth stating.
- `generate_library()` itself accumulates every row in a Python list before building the DataFrame.
  Only `to_file`/`to_df` chunk. The distinction "the *pool* is lazy; a materialising *call* is not"
  should be made explicitly rather than left for a referee to make.

The `stack` branch-evaluation cost (limitation 9) is real and independently confirmed in the code:
`_compute_one` iterates the whole topo-sorted op list per row with no branch pruning; only the
*design-card* write is gated on `op.state.is_active` (`generate_library.py:314`). The authors'
measured 3.03x penalty and their own audit item H3 (~70% speedup) are in `examples/README.md` as
quoted. Reporting this is the single most credibility-buying thing in the record. Keep it.

### `library_algebra` = yes — **supported**

All six primitives confirmed as first-class ops with their own `Operation` subclasses:
`state_ops/{stack,sample,repeat,state_slice,state_shuffle,sync}.py`, `fixed_ops/join.py` (spacer
handling at `join.py:13-33`). `stack` delegates to `statetracker.stack` (`state_ops/stack.py:86-90`)
for a genuine disjoint union of state spaces. Operator sugar verified live.

### `exhaustive_single_scans` = yes — **supported**

Reproduced exactly: `pp.from_seq('ACGTACGT').mutagenize(num_mutations=1, mode='sequential')` →
`num_states=24`, names `mut_00`…`mut_23`, first rows `CCGTACGT`, `GCGTACGT`, `TCGTACGT`.
`deletion_scan(deletion_length=3)` on a 12-mer → 10 states. Codon level: `mutagenize_orf` on a 6-codon
ORF with `missense_only_first` → 114 = 6×19. All the named scan ops exist in `scan_ops/`.

Scope note (does not change the value): "exhaustive" at codon level means exhaustive over *amino
acids*, one representative codon each (`missense_only_first`), not over all synonymous DNA variants —
and the synonymous class cannot be enumerated at all (Finding 3). Say "all amino-acid substitutions",
not "all codon substitutions".

### `sampled_random_mutagenesis` = yes — **supported**

`rng.binomial(num_mutable, self.mutation_rate)` at `base_ops/mutagenize.py:458` — per-position
probability, as claimed. Mutual exclusivity with `num_mutations` is at **144-147**, and the
sequential refusal at **157-159**; the record cites "26-27,138-155", which is drift of ~10 lines.
Fix the citations — a referee who opens line 26 finds a parameter default, not a validation.
Per-op seeding at `generate_library.py:283` (`np.random.SeedSequence([master_seed, op.id, state_val])`)
confirmed verbatim.

### `higher_order_combinatorial` = yes — **supported**

`_build_caches` at `base_ops/mutagenize.py:332-381`, `comb(num_positions, num_mutations) *
(alpha_size-1)**num_mutations` at 351-353, with the genuine non-uniform IUPAC branch at 365-378.
Arithmetic checks: 55×19 = 1,045; C(55,2)×19² = 1,485×361 = **536,085**. Both numbers appear in
`docs/tutorials/dms_gb1.rst` (lines 91, 113, 116) and in the composition table at 263-275, total
**547,230** — all reproduced from the docs, not taken on trust.

### `heterogeneous_components_one_library` = yes — **supported, and independently reproduced**

I ran the record's own construction from scratch and got its exact numbers: three pools of
`seq_length` 8 / 12 / 6 with `num_states` 24 / 10 / 4 → `pp.stack([...])` gives
`num_states=38, seq_length=None`, and `generate_library()` emits 38 rows whose realised lengths are
`{6, 8, 12}`. `state_ops/stack.py:66-71` sets `seq_length=None` when parents disagree. `print_dag`
renders the three-branch DAG. This row is the one PoolParty is most likely to hold alone in the
table, and it is the best-evidenced row in the record. No change needed.

### `combinatorial_motif_place` = yes — **supported; evidence is weaker than reality**

Verified live and it does more than the record shows. Two distinct motif pools, four candidate
positions, `min_spacing=6`:
`insertion_mode='unordered'` → **12** states (6 position pairs × 2 pool→position permutations);
`insertion_mode='ordered'` → **6**. Design cards returned `combination_index`, `starts` (`[0, 10]`,
`[0, 20]`, …) and `names` (`[_ins_0, _ins_1]` flipping to `[_ins_1, _ins_0]` at the permutation
boundary). `flip(mode='sequential')` → 2 states, forward + reverse complement, confirmed on
`GGGCCCAAA` → `TTTGGGCCC`. Spacing constraints, permutation enumeration and per-site provenance are
therefore all demonstrable, not just documented. Attribution nit: those card keys are declared on
`RegionMultiscanOp` (`region_ops/region_multiscan.py:129`), which `insertion_multiscan` composes —
the emitted column is literally `op[3]:insertion_multiscan(region_multiscan).starts`.

### `barcode_generation` = yes — **supported**

Independently exercised at a *harder* setting than the shipped tutorial:
`get_barcodes(num_barcodes=8, length=8, gc_range=(0.4,0.6), min_edit_distance=3, max_homopolymer=2,
seed=7)` → 8 barcodes, measured minimum pairwise Levenshtein distance **3**, every GC exactly 0.50,
no run > 2. Constraints are enforced during generation, not post hoc. File is 579 lines as claimed,
genuinely vectorised (`_gc_filter_batch:40`, `_homopolymer_filter_batch:47` with a cumsum sliding
window, batch candidate draw at 74).

**Presentational warning.** The shipped MPRA tutorial uses `min_edit_distance=1`
(`docs/tutorials/mpra_regulatory_grammar.rst:139`), which for distinct barcodes is a no-op — edit
distance ≥ 1 is just uniqueness. The record's `assay_mpra` evidence says "constraint-checked 8 bp
barcodes ... (gc_range, min_edit_distance)". True but thin. Cite the **operation's** capability
(verified above) rather than the tutorial's setting, or a referee will open the tutorial and find the
distance constraint doing nothing.

### `per_sequence_provenance` = yes — **supported**

I enumerated every `design_card_keys` declaration in the package (20 sites) and the record's list is
accurate except: **`FlipOp`'s key is `"flip"`, not `"strand"`** (`base_ops/flip.py:103`), and the
multiscan keys belong to `RegionMultiscanOp` as noted. `RegionScanOp` additionally declares
`["position_index","start","end","name","region_seq"]` (`region_ops/region_scan.py:98-104`), which
the record omits and which is good evidence. Card filtering and column naming at
`generate_library.py:298-319` as cited. The opt-in caveat is correct and confirmed live: with no
`cards=` the DataFrame has exactly `['name','seq']`.

### `automatic_naming` = partial — **supported (honest downgrade, correctly reasoned)**

Both stated reasons hold. Verified live: `pp.from_seq('ACGT').mutagenize(num_mutations=1,
mode='sequential')` → `name` is `None` on every row, columns `['name','seq']`. With prefixes, names
compose in topological order (`m_00.d_0`, `m_00.d_1`, …). `docs/metadata/naming.rst:27-29` states
the `None` behaviour in the authors' own words. Names are state indices, not variant descriptions.
`partial` is the right call and pre-empts the obvious attack.

### `design_visualization` = yes — **supported**

Both halves of the row's own wording exist and I ran both: `print_dag()` (`pool.py:385`) rendered the
ASCII DAG shown above; `print_library()` (`pool.py:282`) emitted ANSI-coloured mutated bases
(`\x1b[91mC\x1b[0m` on the mutated position). Grep for `matplotlib|pyplot|seaborn|plotly|logomaker`
across `src/` and `docs/` returns **0** — the record's "text only" disclaimer is accurate.

Calibration against peers: `pydna` scores `yes` in part on *ASCII* `figure()` output
(`extractions/pydna.md:113-118`) and `ledidi` on plots alone with no graph view. PoolParty has the
graph view *and* highlighted sequences, i.e. literally both halves of the row. `yes` is consistent
with the survey's own bar. Keep the honest qualifier in the cell.

---

## 2. Block C — the rows that would discredit the whole table

All five re-derived independently. **No quiet claim of genome awareness anywhere.**

| Key | Record | My check | Verdict |
|---|---|---|---|
| `genome_coordinates` | partial | `fixed_ops/from_fasta.py` (166 lines) takes `(chrom, start, stop, strand)` tuples, pyfaidx-indexed, `-` strand reverse-complements, `start > stop` wraps; names carry `{chrom}:{start}-{stop}({strand})`; card keys are `seq_name`/`seq_index` only. No BED/GTF/interval file, no build awareness, no back-mapping. | **supported** |
| `transcript_models` | no | `grep -riE 'exon\|intron\|splice\|transcript\|gtf\|gff\|ensembl\|refseq' src/poolparty/ --include=*.py` → **0 hits** | **supported** |
| `exon_intron_split_codons` | no | `OrfRegion` is `(name, seq_length, frame)`; `_resolve_frame` (`mutagenize_orf.py:21`) resolves one frame for one contiguous region | **supported** |
| `vcf_vep_output` | no | export `Literal["csv","fasta","tsv","jsonl"]` (`export_mixin.py:225`, re-validated at 351); `grep -riE 'vcf\|vep\|hgvs'` → 0 | **supported** |
| `consequence_annotation` | no | `mutation_type` specifies what to *generate*; no classifier, no vocabulary, no VCF consumer | **supported** |

**Cross-tool calibration of `genome_coordinates = partial` holds.** `tangermeme` scores `yes` on
BED + FASTA via pyfaidx; `seqpro` scores `partial` for interval I/O *without* sequence fetch;
`biopython` `partial`; `dnachisel` `no`. PoolParty is seqpro's mirror image — fetch without interval
files — so `partial` is exactly where it belongs. Do not let anyone push it to `yes`, and do not
concede it to `no`.

The `consequence_annotation` entry is the best-written cell in the record: it names the near-miss
(`synonymous`/`missense`/`nonsense` generation), explains precisely why it is not the capability, and
denies it anyway. That is the paragraph that will buy the table its credibility.

---

## 3. Block D

- **`primer_design` = no — supported.** `grep -riE 'primer|melting|anneal|hairpin'` over `src/` returns
  one hit: a code *comment* at `data/restriction_enzymes.py:202` ("enzymes to avoid in Gibson
  primers"). Correctly identified as a lexical near-miss.
- **`codon_optimization` = partial — supported.** `codon_table.py:33-56` is the Kazusa *H. sapiens*
  table, frequency-ordered high→low, with the ordering's significance documented in the source comment
  at 28-32. Live: `pp.reverse_translate('MKV', codon_selection='first')` → `ATGAAGGTG` (M→ATG, K→AAG,
  V→GTG = the first, most-frequent codon in each list). Note the record's own live check works only
  because `reverse_translate` string-promotes to a `ProteinPool` (`reverse_translate.py:88-89`);
  passing a `DnaPool` raises `TypeError`. Strategy set really is {first, uniform random} — no CAI, no
  GC target, no harmonisation, no constraint-aware selection.
- **`synthesis_constraints` = partial — supported, but for a different reason than given.** All five
  filters exist (`filter_mixin.py:20,65,106,151,194`) with a curated IUPAC- and RC-aware enzyme
  database and seven cloning presets. Live: a GC filter over `from_iupac('NNN')` kept 32 of 64 and
  emitted 32 `NullSeq` rows with the acceptance-rate warning. The record justifies `partial` by
  "enforcement is reject-only" — but *rejection is enforcement*, and the row asks for "constraint
  checking/enforcement", not repair. The defensible justification is the one the record also gives:
  **missing constraint types** — no oligo-length/vendor cap, no Tm, no secondary structure, no
  pool-level repeat or cross-hybridisation analysis, no background k-mer screen. Calibration:
  `oligopoolcalc` earns `yes` on exactly those missing types plus constructive `split`/`pad`;
  `mpradesign` and `mpranator` earn `partial` on *narrower* check sets than PoolParty's. PoolParty
  sits at the top of `partial`. Lead with the missing types, not with "reject-only".
- **`degenerate_iupac_codons` = yes — supported.** `pp.from_iupac('NNK', mode='sequential')` →
  **32** states, first rows `AAG, AAT, ACG, ACT`. NNK handled exactly. The "expands, does not emit a
  degenerate oligo string" caveat is correct and worth keeping, since the row's peers (CodonGenie,
  MPRAnator) emit degenerate strings.
- **`negative_control_generation` = yes — supported, and the evidence is unnecessarily modest.**
  Live on `ACGGTTACCGATTGCAACGTTAGC`: 5 distinct `dinuc` shuffles, all with dinucleotide
  frequency vectors identical to the source (checked with a Counter), plus 3 distinct `mono`
  shuffles. The Euler-path claim is the authors' own (`shuffle_seq.py:41-43`). Two notes: (a) the
  record's own live check used a perfect tandem repeat (`ACGT`×4) whose dinucleotide graph is a single
  cycle — that sequence has exactly **one** valid dinuc shuffle, so it returns the input unchanged;
  it is a bad demo sequence and if it reached a referee it would look like a bug. Use a non-cyclic
  sequence. (b) The row definition names "reverse/complement controls" explicitly and PoolParty has
  `rc` and `flip` as first-class ops — add them; the record's "no operation is NAMED control" caveat
  gives away more than it needs to. Finally, note that this row is **not** a PoolParty exclusive:
  `extractions/mpranator.md:169` records scramble/reverse/complement as first-class in MPRAnator and
  `seqpro`'s Rust `k_shuffle` is the same construction. Do not frame it as a differentiator.
- **`ml_model_in_loop` = partial — supported.** `ScoreOp._compute_core` (`fixed_ops/score.py:113-119`)
  returns `parent_seq` unchanged with `{card_key: fn(clean_seq)}` — a literal passthrough, one
  sequence at a time, `mode="fixed"`, `num_states=1`. Model output cannot steer generation; only
  `filter` can reject. `partial` is right and the reasoning is exact.
- **`readout_analysis` = no — supported.** `grep -riE 'enrichment|readout|fastq|bam'` over `src/` → 0.

---

## 4. Block E — encoding defects (Finding 2)

Facts all check out; the *values* do not carry them.

- **`interface` — value `"yes"` is unsupported as encoded.** Facts confirmed: `pyproject.toml`
  declares no `[project.scripts]` and no entry points; `_POOL_FACTORY_MAP`/`_DNAPOOL_FACTORY_MAP`
  exist at `__init__.py:331`/`373` and are bound onto `Pool` at 390/395; `beartype` wraps public
  factories (visible in every traceback I triggered); `requires-python = ">=3.10"`.
  **Re-encode as `yes (Python API only)`** — the form used by `final/seqpro.md` and required by
  `final/valiant.md` §8.
- **`license` — supported, re-encode as `yes (MIT)`.** `poolparty/LICENSE` "MIT License / Copyright
  (c) 2025 Justin B. Kinney"; `pyproject.toml:11 license = "MIT"`; PyPI `license_expression: "MIT"`
  (the legacy `info.license` field is `null`, so cite `license_expression`); `statetracker` also MIT.
- **`installable_today` = yes — supported.** PyPI JSON pulled live: `poolparty` 0.1.1, wheel + sdist
  both `2026-04-06T21:10`, 0.1.0 at `21:03` the same day, `requires_python >=3.10`, runtime deps
  exactly `numpy>=1.20, pandas>=1.3, beartype>=0.22.9, statetracker>=0.1.0, pyfaidx>=0.8.1,
  typing_extensions>=4.0`. CI matrix confirmed: ubuntu × {3.10, 3.11, 3.12} + macOS 3.11 + Windows
  3.11, running both `statetracker` and `poolparty` test suites. Import and execution verified.
- **`last_activity` — value `"yes"` is unsupported as encoded.** The row asks for a date.
  **Re-encode as `2026-04-07 (last commit 1bb0179); PyPI 0.1.1 2026-04-06`.** The `__version__ =
  "0.1.0"` / distribution-0.1.1 mismatch is real and I reproduced it: the installed distribution
  reports `0.1.0` from `importlib.metadata` too, so the skew is in the editable install as well as
  the module constant. Keep it as the minor-hygiene note it is.
- **`documented_examples` — value `"yes"` is unsupported as encoded.** The row asks for a **count**.
  Also two counts in the record are wrong low: `docs/operations/` holds **62** `.rst` files, not 57
  (of which ~50 are per-operation pages and the rest are concept/grouping pages); `tests/` holds
  **77** `.py` files (the 78th entry is `__pycache__`) — the 2,906 `test_` function count is exact and
  I reproduced it. **Re-encode as e.g. `3 (quickstart + 2 full tutorials); plus ~50 per-operation
  reference pages and 7 concept pages`.** The gitignore finding is correct and important: `.gitignore`
  line 70 is `/poolparty/examples/` (also `/poolparty/notebooks/` at 69 and `/poolparty/benchmarks/`
  at 68), `git ls-files` returns nothing for any of them, so none of the three notebooks ships.

---

## 5. Things the record missed

1. **`get_barcodes` is eager** (Finding 1) — belongs in the `lazy_generation` caveat, *and* its
   fail-fast constraint-infeasibility error is a positive worth claiming.
2. **Shared DAG nodes** (Finding 4) — the best evidence for `composable_operations`, unclaimed.
3. **Synonymous variants cannot be exhaustively enumerated** (Finding 3) — belongs in
   `stated_limitations`.
4. **`min_spacing`/`max_spacing` layout constraints and permutation enumeration are demonstrable**,
   not merely documented (12 vs 6 states, measured) — strengthen `combinatorial_motif_place`.
5. **Constraint predicates are exported as standalone public functions** — `calc_gc`,
   `calc_complexity`, `calc_dust`, `has_homopolymer`, `has_restriction_site` are in `__all__` and
   usable outside the DAG. Small, but it is API surface the record does not credit.
6. **Deletion semantics.** `deletion_scan` emits gap characters (`TTTTTTTT---T`), it does not shorten
   the sequence; `clear_gaps` is a separate op. Anyone reading "deletion scan" and expecting shorter
   oligos will be surprised. Worth one clause in the record, since it interacts with fixed-length
   oligo synthesis.
7. **`sync()` is under-evidenced.** It is what makes 1:1 barcode pairing work instead of a Cartesian
   product, and it is now the *default* for `replace_region` (`replace_region.py:19-20`, the
   `[Unreleased]` BREAKING change). It is listed under `library_algebra` but never explained; it is a
   genuine differentiator against tools that force the user to zip libraries by hand.
8. **`assay_dms`'s 547,230 is a member count, not a distinct-sequence count.** With
   `mutation_rate=0.1` over 55 codons, Binomial(55, 0.1) puts ~0.3% of the 10,000 random draws at 0
   mutations, ~1.9% at 1 and ~5.6% at 2 — roughly **774 expected draws that duplicate the WT, single-
   or double-mutant arms**. Say "547,230 library members"; do not say "547,230 distinct sequences".
9. **Path corrections.** `CITATION.cff`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` and
   `.pre-commit-config.yaml` live at the **monorepo root**, not in `poolparty/` as the record's
   sources table says. The claims hold; the paths do not.
10. **Repo-name skew.** The working copy is `poolparty-statecounter/` and contains *both* a
    `statetracker/` and a `statecounter/` directory, while PyPI `project_urls` gives the canonical
    repo as `github.com/jbkinney/poolparty-statetracker`. Make sure the manuscript and `CITATION.cff`
    cite the published URL, not the local directory name.
11. **`src/poolparty/marker_ops/` contains nothing but `__pycache__`** — a dead package directory.
    Cosmetic, but a referee cloning the repo may notice.

---

## 6. Where the record is too modest

- **`negative_control_generation`** — see §3. `rc` and `flip` are first-class and are named in the row
  definition; the "no operation named control" hedge concedes more than the facts require.
- **`library_first_class_object`** — the Party caveat reads as harsher than it is, since a default
  Party exists from import.
- **`combinatorial_motif_place`** — spacing constraints and per-site provenance are demonstrable and
  the record only asserts them.
- **`documented_examples`** — 62 operation pages counted as 57.

## 7. Where a referee will still push, and the answer

| Attack | Answer already in the record? |
|---|---|
| "Your DAG is lazy but `mutagenize` isn't" | Yes, stated. Extend to `get_barcodes`. |
| "`stack` costs 3x because you evaluate dead branches" | Yes, with the authors' own measurement and audit item. Strongest possible framing. |
| "Naming is not actually automatic" | Yes — scored `partial` pre-emptively. |
| "You advertise in-silico use with no shipped tutorial" | Yes — scored `partial`, quoting the authors' own `examples/README.md`. |
| "Your notebooks aren't distributed" | Yes, with the `.gitignore` line number. |
| "You can't enumerate synonymous variants" | **No — add it.** |
| "Your DAG is really a tree" | **No — add the shared-node demonstration.** |

---

## Verdict tally

`supported` 26 · `understated` 1 (`documented_examples` counts) · `unsupported as encoded` 3
(`interface`, `last_activity`, `documented_examples` value fields) · `overstated` 0 · `wrong` 0.

The record is publishable after the four findings are addressed. Its most valuable property is that
the limitations section is longer, more specific and more damaging than anything a referee is likely
to produce unaided — including a measured 3x self-inflicted performance penalty and a pointer to the
authors' own unimplemented fix. Keep all of it.
