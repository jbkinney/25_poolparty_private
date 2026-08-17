# PoolParty — repair change log

**Record repaired:** `revision/tool_survey/final/poolparty.md` (edited in place)
**Audits processed:** `citation_audit/poolparty.md` (7 NOT-FOUND + 16 other), `factcheck/poolparty.md`
(10 in §A, 14 in §B, 1 in §C)
**Repaired:** 2026-08-14
**Method:** every finding re-verified against primary source before any edit — repository source read
directly, ~30 read-only executions against `poolparty-statecounter/.venv/bin/python` (Python 3.12,
`PYTHONDONTWRITEBYTECODE=1`, all runs from `/tmp` or the session scratchpad), live `pypi.org` JSON and
`api.github.com` fetches, and `git`/`grep`/`find` over the read-only checkout. Nothing was written
outside `revision/tool_survey/`; the PoolParty repository was not modified.

**Outcome: 23 findings applied · 6 rejected · 15 escalated (as 4 questions).** No capability value
was changed. Several findings appear in both audits (citation NOT-FOUND 1/2/4/6/7 = factcheck
A3/A1/A4/A6/A9; citation other 11 = factcheck A2) and are recorded once, under the applied item that
carries the edit. The 11 §B coverage findings and the §C balance finding are escalated as two grouped
questions rather than twelve separate ones, because they are one decision each.

---

## APPLIED

### 1. Party-error reachability (citation NOT-FOUND 1 = factcheck A3)

* **Was:** "the failure is reachable only after an explicit `pp.init()` or outside an exited
  `with pp.Party()` block."
* **Now:** the failure is **not** reachable by either route, with `party.py:62-89` / `:207-226` cited.
* **Verified:** read `party.py` — `init()` creates a new default Party, calls
  `_state_manager.__enter__()` and sets `_active_party` (`:62-89`); `Party.__exit__` restores
  `_previous_party` (`:218-226`). Live: `pp.init(); pp.from_seq('ACGT')` → OK;
  `with pp.Party(): pass` then `pp.from_seq('ACGTA')` → OK.
* Source line also extended with `party.py:62-89,207-226`.

### 2. "40 method-bound factories" / "every operation is a `Pool` method" (citation NOT-FOUND 2 = factcheck A1)

Corrected in **three** places (§What PoolParty is, `composable_operations`, `interface`).

* **Was:** "Every operation is `Pool(s) → Pool`, exposed both as a module function and as a `Pool`
  method … bound onto `Pool` at 390/395 — **40 distinct method-bound factories**".
* **Now:** "Most operations…", **39** operation factories with method forms, the loops at `:390/395`
  copy docstrings rather than bind, and the module-function-only operations are named.
* **Verified live:** `len(_POOL_FACTORY_MAP)` = 31 (one entry is `generate_library`, not an
  operation), `len(_DNAPOOL_FACTORY_MAP)` = 9 → 39 operation factories. Read
  `__init__.py:390-399`: both loops are `if hasattr(...): …__doc__ = …`. `hasattr` checks confirm
  `stack`, `join`, `sync`, `region_scan`, `region_multiscan`, `from_seq`, `from_seqs`, `from_fasta`,
  `get_barcodes` are **not** methods of `Pool` or `DnaPool`.

### 3. `sync` is not an `Operation` subclass (citation NOT-FOUND 3)

* **Was:** "All primitives are first-class `Operation` subclasses"; table note for `sync` said only
  "see below".
* **Now:** "All primitives **but `sync`**…"; the `sync` table cell records that it is a module-level
  function that rewires pool states in place and returns `None` (`sync.py:26-68`).
* **Verified:** `grep 'class .*Op' src/poolparty/state_ops/*.py` returns `RepeatOp`, `SampleOp`,
  `StackOp`, `StateShuffleOp`, `StateSliceOp` — no `SyncOp`. `sync.py` is a single `@beartype def
  sync(pools) -> None` ending `pool.state = w` at `:68`.
* The same fix also removes the "Every operation is `Pool(s) → Pool`" absolute (see item 2).

### 4. `get_kmers` does not hold its source (citation NOT-FOUND 4 = factcheck A4)

Corrected in **two** places (`lazy_generation` bullet, limitation 18).

* **Was:** "Source pools hold their sources: `from_seqs`, `from_fasta`, `get_kmers`."
* **Now:** `from_seqs`/`from_fasta` only; `get_kmers` stores length/alphabet/total (`:146`) and
  derives each k-mer from the state index (`_value_to_kmer:198`, `_compute_core:212`).
* **Verified live:** `pp.get_kmers(length=8, mode='sequential')` → `num_states=65536`; dumping
  `vars(op)` shows no collection attribute of length > 4.

### 5. "Any operation can target a named region" (citation NOT-FOUND 5)

* **Was:** "so any operation can target a named region by name rather than by index".
* **Now:** "any **region-aware** operation…", with the state-level operations that take no `region`
  argument named.
* **Verified live** with `inspect.signature`: `region` absent from `stack`, `sample`, `repeat`,
  `state_slice`, `state_shuffle`, `sync`, `join`, `filter`, `get_barcodes`, `materialize`; present on
  `mutagenize`, `rc`, `flip`, `score`, `recombine`, `shuffle_seq`, `translate`.

### 6. "Silently shrinks" / "nothing re-drawn" (citation NOT-FOUND 6 = factcheck A6)

Corrected in **two** places (`synthesis_constraints` closing note, limitation 9).

* **Was:** "a filter-heavy design **silently shrinks** — hence the `min_acceptance_rate` warning …
  Nothing is fixed, re-drawn under constraint, or optimised."
* **Now:** neither silent nor an automatic shrink — with `num_seqs=N` generation advances to further
  states until N valid rows are collected and warns explicitly when short; with `num_cycles` the
  output does shrink. The "nothing repaired / re-drawn *under* the constraint" point is kept.
* **Verified:** read `generate_library.py:100-190` (the `while len(rows) < num_seqs` loop appends
  only non-null rows and keeps advancing `state`); three warning paths at `:150-184`. Live on a
  64-state `from_iupac('NNN')` pool filtered to 8 survivors: `num_seqs=8, discard_null_seqs=True`
  → **8 rows, no warning**; `num_seqs=20` → 8 rows **with** `Reached max_iterations (64) …`;
  `num_cycles=1` → 8 rows with the same warning.

### 7. "All dependencies pure-Python, no compiled extensions" (citation NOT-FOUND 7 = factcheck A9)

Corrected in **two** places (`installable_today`, §Availability).

* **Now:** PoolParty itself ships no compiled extension, but `numpy` and `pandas` do.
* **Verified:** `find .venv/.../{numpy,pandas} -name '*.so'` → **63** files, e.g.
  `numpy/_core/_multiarray_umath.cpython-312-x86_64-linux-gnu.so`,
  `pandas/_libs/algos.cpython-312-x86_64-linux-gnu.so`.

### 8. Adversarial-review tally (citation other 1)

* **Was:** "29 of 33 values `supported` … The four non-`supported` verdicts … plus two rows the
  reviewer judged *too modest* (`combinatorial_motif_place`, `negative_control_generation`)."
* **Now:** the review's own tally, quoted: 26 `supported` · 1 `understated` · 3 `unsupported as
  encoded` · 0 `overstated` · 0 `wrong`, over the 30 values it reviewed; the two "too modest" rows
  are stated to have been scored **`supported`**.
* **Verified:** `reviews/poolparty.md:394-397` gives exactly those counts; `:14` says "26 of 30
  capability values are supported"; `:168` and `:280` label `combinatorial_motif_place` and
  `negative_control_generation` **supported**.

### 9. "500-3000x" construction/generation gap (citation other 2)

* **Now:** **202-710x** (SpliceAI 503x, MPRA 202x, DMS 710x), with a note that the README's own
  "500-3000x" cross-pairs the fastest construction with the slowest generation.
* **Verified:** `examples/README.md:105-113` table — construction 0.05 / 0.03 / 0.24 s, generation
  10.12 / 15.08 / 170.32 s. Ratios recomputed: 202.4, 502.7, 709.7. The record's own paired figures
  do not yield 3000x.

### 10. `join` evidence location (citation other 3)

* **Was:** "end-to-end concatenation, optional spacer (`:13-33`)".
* **Now:** concatenation `spacer_str.join(seqs)` at `:89`, spacer argument at `:45`, length
  accounting at `:78-82`.
* **Verified:** read `fixed_ops/join.py` in full — `:13-39` is `_make_join_style_combiner`
  (style handling only).

### 11. "8,000 unique CREs" (citation other 4)

Corrected in **two** places (`combinatorial_motif_place`, `assay_mpra`).

* **Now:** "8,000 CRE **states**", with a note that the tutorial's word *unique* is not warranted
  because the 1,000 position configurations are drawn with replacement.
* **Verified:** reproduced the documented pre-repeat MPRA pipeline verbatim from
  `docs/tutorials/mpra_regulatory_grammar.rst:41-122`. `num_states` = 8000 exactly, but distinct
  `seq` strings = 7992 / 7992 / 7968 / 7976 at seeds 0 / 7 / 42 / 123. The exact deficit is
  stochastic, so no specific number was written into the record beyond the observed range.
* *(The audit's companion claim that "collided arrangements receive six not three barcodes" was not
  applied: the record never makes the three-barcodes-per-arrangement claim — only the shipped
  tutorial does.)*

### 12. `readout_analysis` grep "0 hits" (citation other 6)

* **Now:** **7 hits, all `BamHI`**, cited at `data/restriction_enzymes.py:28,169,195,205`;
  `pool_mixins/filter_mixin.py:217,249`; `utils/seq_properties.py:319` — lexical false positives from
  `bam`, nothing readout-related.
* **Verified:** re-ran the record's own regex; 7 hits, listed above. The `no` value is untouched and
  unaffected.

### 13. `matplotlib` "only matches anywhere" (citation other 7)

Corrected in **two** places (`design_visualization`, the adjudication table row) — the same false
statement appears twice and it would be incoherent to fix one and not the other.

* **Now:** matches elsewhere are SVG producer metadata in `docs/_build/html/` **and** in the
  **committed** `docs/_static/images/figure4b_g.drawio.svg`.
* **Verified:** `grep -oiE 'matplotlib…' docs/_static/images/figure4b_g.drawio.svg` → `Matplotlib
  v3.10.7`; `git ls-files docs/_static/images/` lists that file; `git ls-files docs/_build` is empty.
  The narrowed claim (0 hits in `src/**/*.py` and `docs/**/*.rst`) re-verified: 0 and 0.

### 14. Unseeded "419 barcodes" (citation other 8, first half)

* **Now:** the quoted `ValueError` carries `N`, with a note that the call is unseeded so N is not
  reproducible evidence (419 first measured, **424** on re-execution here).
* **Verified live:** three consecutive runs of the exact displayed call all returned 424 within one
  session; the audit reports 414 and the prior review a third value.
* *The 0.043 s timing in the same paragraph was **not** changed — see REJECTED 5.*

### 15. `from_fasta` coordinate names are batch-mode only (citation other 9)

* **Now (in `automatic_naming`):** "`from_fasta` **in batch (list-of-coordinates) mode**…", with the
  single-coordinate path noted as delegating to `from_seq` at `:109-125` with no such name.
* **Verified:** read `fixed_ops/from_fasta.py:107-164` — `is_single` branch returns `from_seq(...)`
  with no `seq_names`; the batch branch builds `{chrom}:{start}-{stop}({strand})` at `:145-153`.
* The `genome_coordinates` section already qualified this correctly and was left alone.

### 16. `pp.set_genetic_code` does not reach `reverse_translate` (citation other 10)

* **Now (in `codon_optimization`):** a custom table can be supplied **per call**; `pp.set_genetic_code`
  does *not* reach `reverse_translate`, which builds its own `CodonTable` from its own argument
  (default `"standard"`, `reverse_translate.py:143`).
* **Verified live:** after `pp.set_genetic_code({'K': ['AAA','AAG'], …})`, default
  `pp.reverse_translate('K', codon_selection='first')` still yields `AAG`; passing
  `genetic_code={'K':['AAA','AAG']}` explicitly yields `AAA`.

### 17. `num_states` / `seq_length` "are exact" (citation other 11 = factcheck A2)

Corrected in **two** places (§What PoolParty is, Additional capability 3).

* **Now:** `num_states` is an exact **state** count — not necessarily a distinct-sequence count —
  and `seq_length` is exact **when statically determined** (`None` for variable-length pools).
* **Verified:** `pool.py:143-146` returns `None` for variable-length pools; the record's own
  heterogeneous-stack example already shows `lib.seq_length is None`; `repeat`, random draws and
  filtering all break the state-count-equals-distinct-sequence identity (`state_ops/repeat.py`,
  `generate_library.py:29-46`).

### 18. "sliced, shuffled, sampled or **split** without generating it" (citation other 12)

* **Now:** "sliced, shuffled or sampled", with a parenthetical that `split` exists only on the
  companion `statetracker` `State` API.
* **Verified live:** `hasattr(pp,'split')`, `hasattr(Pool,'split')`, `hasattr(DnaPool,'split')` all
  `False`; `hasattr(statetracker,'split')` `True` (`statetracker/ops/split_op.py:8`).

### 19. Deletion scans and `deletion_marker=None` (factcheck A5)

* **Was (limitation 4):** "Deletion scans emit gap characters, **they do not shorten the sequence**
  … `clear_gaps` is a separate operation the user must remember."
* **Now:** "emit gap characters **by default**"; `deletion_marker=None` produces genuinely shorter
  sequences (`scan_ops/deletion_scan.py:100-107`), with the documented "True deletion" example at
  `docs/operations/deletion_scan.rst:111-134`.
* **Verified:** read both files; the docs page shows `deletion_scan(deletion_length=2,
  deletion_marker=None)` on an 8-mer giving `seq_length=6, num_states=7`.

### 20. "API reference — Autodoc of all public functions/classes" (factcheck A8)

* **Now:** "Autodoc of *selected* public functions/classes — **not complete**", listing the omissions.
* **Verified:** `docs/api.rst` uses explicit `autofunction`/`autoclass` directives (no blanket
  `automodule`). Name-by-name grep returns **0** occurrences for `get_barcodes`, `annotate_orf`,
  `annotate_region`, `translate`, `reverse_translate`, `stylize_orf`, `set_genetic_code`,
  `set_progress_mode`, `fixed_operation`, `DnaPool`, `ProteinPool`, `to_df`, `to_file`, `clear_tags`,
  `calc_gc`, `calc_complexity`, `calc_dust`, `has_homopolymer`, `has_restriction_site`,
  `ENZYME_SITES`, `get_preset_enzymes`.

### 21. §B additions applied as brief clauses

Three §B findings were verified and are concrete restrictions rather than coverage judgments, so they
were added as a clause in the relevant existing entry (no new sections):

* **B1 — the 1,000,000-state sequential ceiling.** Added to limitation 5 ("Undocumented mode
  restrictions"). Verified: `Operation.max_num_sequential_states = 1_000_000` (`operation.py:33`),
  raising in `validate_num_states` (`:37-50`). Live: `pp.get_kmers(length=11, mode='sequential')` and
  `pp.from_iupac('N'*11, mode='sequential')` both raise `Number of states (4194304) exceeds
  max_num_sequential_states (1000000). Use mode='random' instead.` Grep confirms the ceiling appears
  in no `.rst` file, so limitation 5 is the correct home.
* **B11 — protein-pool boundaries.** Added to Additional capability 11: no public protein-sequence
  source and no direct amino-acid substitution operation. Verified live:
  `set(dir(ProteinPool)) - set(dir(Pool))` == `{'reverse_translate'}`.
* **B12 — export is `DnaPool`-only.** Added to Additional capability 12. Verified:
  `DnaPool(Pool, DnaMixin, FilterMixin, ExportMixin)` vs `ProteinPool(Pool, ProteinMixin)`
  (`dna_pool.py:7`, `protein_pool.py:13`); live `hasattr(ProteinPool,'to_df')` → `False`.

---

## REJECTED

### 1. citation other 14 — "168-state nested sub-pipeline is uncited"

**Why wrong:** the number reproduces exactly under the natural reading of the record's own text.
**Evidence:** `pp.apply_at_region(pp.from_seq('AAAA<r>ACGTACGT</r>TTTT'), 'r', lambda p:
p.mutagenize(num_mutations=1, mode='sequential', prefix='m').deletion_scan(deletion_length=2,
mode='sequential', prefix='d'))` → `num_states = 168` (24 × 7), first rows named `m_00.d_0`,
`m_00.d_1`, `m_00.d_2` — precisely what the record shows. The audit's complaint is about a missing
transcript, not a wrong fact. Record left untouched.

### 2. citation other 15 — "GC filter kept 32 of 64 is uncited"

**Why wrong:** reproduces exactly. **Evidence:** `pp.from_iupac('NNN', mode='sequential')` → 64
states; `.filter_gc(min_gc=0.5)` → 32 kept / 32 `NullSeq`; `.filter_gc(max_gc=0.5)` → likewise 32/32.
The partition is forced by the combinatorics (8 + 24 + 24 + 8 by G/C count), not by a lucky
threshold. Record left untouched.

### 3. citation other 13 — "24 read-only executions and repository non-modification are uncited"

**Why not applied:** this asserts no factual error in the record; it is a methodological-transparency
complaint about a first-person provenance claim I cannot falsify (and the audit itself concedes the
individual runtime outputs are reproducible). There is nothing to correct surgically. Record left
untouched. *(Noted for the authors: if a referee-facing artefact is wanted, keep a command transcript
next time.)*

### 4. citation other 16 — "dead-link DOI placeholder"

**Why not applied:** confirmed dead (`https://doi.org/XXXX` → HTTP 404), but the record is already
correct about it: it names it a placeholder and housekeeping item 5 prescribes replacing it. Nothing
in the record is wrong. Record left untouched. *(The other two URLs check out: `poolparty.readthedocs.io`
→ 200, `github.com/jbkinney/poolparty-statetracker` → 200.)*

### 5. citation other 8, second half — "0.043 s barcode timing is not reproducible"

**Why wrong here:** the timing reproduces to within noise in this environment. **Evidence:**
`get_barcodes(num_barcodes=24000, length=8)` measured **0.039 s** with exactly 24,000 strings in
`op._barcode_strings`. The audit's 0.091 s is a different machine; a one-off wall-clock number is not
an error. Only the unseeded barcode *count* was corrected (APPLIED 14).

### 6. factcheck A10 — "'first released 2026-04-06' contradicts the project record"

**Why wrong:** the record's statement is accurate about PyPI, which is what "released" means for a
distributed package. **Evidence:** live `pypi.org/pypi/poolparty/json` — 0.1.0 wheel
`2026-04-06T21:03:21Z` + sdist `21:03:22Z`; 0.1.1 wheel `21:10:37Z` + sdist `21:10:39Z`. There are no
earlier artefacts. `CHANGELOG.md:34` does date `[0.1.0]` to 2026-04-03 and `CITATION.cff` says the
same, but the record already flags that skew (housekeeping item 5) rather than repeating it. Record
left untouched. *(Flagging for the authors: housekeeping item 5 names `CITATION.cff` only —
`CHANGELOG.md:34` carries the same 2026-04-03 date and would need the same decision.)*

### 7. factcheck C — see ESCALATED 4 (not rejected; referred)

---

## ESCALATED

### 1. `last_activity` — the local checkout is stale; the canonical repository has moved

* **Finding:** citation other 5. **Verified true.**
* **Evidence:** live `api.github.com/repos/jbkinney/poolparty-statetracker/commits?sha=main` — tip is
  `9d6a205` "Fix workflow file to properly test different Python versions", **2026-04-13T20:55Z**,
  preceded by `b3306a7` "Add more Python versions to testing workflow", 2026-04-13T20:51Z, then
  `1bb0179` 2026-04-07. The repo `pushed_at` is 2026-06-19, not archived. Locally,
  `git cat-file -t 9d6a205` fails and local `origin/main` still points at `1bb0179` — the checkout
  has not been fetched since 2026-04-08.
* **Why escalated, not edited:** the `last_activity` row's **table cell is its value**
  (`2026-04-07 (commit 1bb0179); PyPI 0.1.1, 2026-04-06`). Correcting it changes what gets typeset
  into the survey table, which the instructions reserve for the ratings pass. The record was left
  untouched at that point.
* **Question for the authors:** should `last_activity` report the canonical `main` tip
  (2026-04-13, `9d6a205`) or the assessed working tree (2026-04-07, `1bb0179`)?
* **Options:** (a) re-cell as `2026-04-13 (commit 9d6a205, canonical main); PyPI 0.1.1, 2026-04-06`
  and keep the `1bb0179` reference in the "Version assessed" header where it correctly describes the
  tree that was read; (b) keep 2026-04-07 and add "as assessed; canonical `main` has since advanced
  to `9d6a205`, 2026-04-13 (CI-workflow commits only)"; (c) re-fetch the checkout and re-run the
  affected checks before deciding. Note the two new commits are CI-workflow-only, so no capability
  evidence is affected — but the record's CI description (`test.yml`: ubuntu × 3.10/3.11/3.12 +
  macOS/Windows on 3.11) is derived from the stale tree and both new commits touch that file.

### 2. `synthesis_constraints` — "no oligo-length or vendor-capacity limit"

* **Finding:** factcheck A7. **Partly verified.**
* **Evidence:** there is no dedicated oligo-length/vendor-capacity constraint object — but the
  generic `filter` accepts any predicate and `docs/operations/filter.rst:97-122` ships a worked
  example that is *exactly* a length constraint (`pp.filter(seqs, lambda s: len(s) == 8)`, with the
  printed output). So a categorical "no oligo-length limit" is falsifiable from the project's own
  documentation.
* **Why escalated, not edited:** the phrase is load-bearing. It is the **first item** in the
  `synthesis_constraints` "Why not `yes` — the missing constraint types" list, and the same list is
  repeated as limitation 10; that list *is* the justification for `partial`. Weakening it changes
  what the value rests on, which the instructions reserve for the ratings pass. Record left untouched
  in both places.
* **Question:** does a user-supplied `len(s)` predicate count as oligo-length constraint checking for
  this row, or does the row require a vendor-capacity/synthesis-length concept the tool models itself?
* **Options:** (a) keep the categorical denial and accept the referee risk; (b) narrow to "no
  *built-in* oligo-length or vendor-capacity limit — a length check is expressible as a generic
  `filter` predicate, and the docs demonstrate one — and no constructive split/pad"; (c) drop the
  item from the missing-types list and lean on the remaining five (Tm, secondary structure,
  repeat-content/cross-hybridisation, background k-mer screen, split/pad), which are unaffected.
  Note (b) and (c) both leave `partial` intact on the other five types.

### 3. §B coverage additions — 11 of the 14 "incomplete" findings

* **Findings:** factcheck B2, B3, B4, B5, B6, B7, B8, B9, B10, B13, B14. **All verified as facts**;
  what is escalated is whether the omissions are material.
* **What they ask for:** add `from_motif` PWM sampling (random-only); `from_seqs`/`get_kmers` as
  source-library interfaces; `recombine`'s constraints; `mutagenize_scan` / `shuffle_scan` /
  `subseq_scan` as distinct geometries; `deletion_multiscan` / `replacement_multiscan`;
  `region_scan` / `region_multiscan` as design operations; `clear_tags` vs `clear_annotation` vs
  `clear_gaps`; `materialize`; user extension via `Operation` subclassing and `fixed_operation`;
  `slice_seq` / case ops; `stylize_orf`.
* **Why escalated, not edited:** each is a judgment about whether an omission is material, which the
  instructions reserve for the authors — and they interact directly with factcheck §C (below), which
  says the record is *already* capability-heavy. Applying all eleven would add roughly a dozen new
  capability claims and push the balance further in the direction §C objects to. They also fall
  outside the row set the record is explicitly scoped to. Record left untouched.
* **Question:** which, if any, of these eleven should enter the record, and in which section?
* **Options:** (a) none — the record is scoped to the row set plus a deliberately short "Additional
  capabilities" list; (b) add only **B10** (extensibility via `Operation` subclassing and the exported
  `fixed_operation`), which is the one the paper itself claims (`latex/main.tex:146-148,239`) and
  which materially strengthens `composable_operations`; (c) add B10 plus the three that name
  *restrictions* rather than features (B2 `from_motif` is random-only; B9 `materialize` severs the
  DAG; B4 `recombine` requires ≥2 equal-length sources); (d) add all eleven and accept the §C balance
  cost. *(B1, B11 and B12 were applied — see APPLIED 21 — because each states a restriction with an
  obvious existing home.)*

### 4. §C — balance between capabilities and limitations

* **Finding:** factcheck §C: "capability-heavy overall, not limitation-heavy"; capabilities get long
  subsections, demonstrations and arithmetic plus repetition in "Additional capabilities", while
  limitations occupy a shorter closing inventory.
* **Why escalated:** explicitly a judgment about emphasis and allocation of space, reserved for the
  authors. No edit made.
* **Question:** is the current allocation acceptable for a subject-of-the-paper record, given the
  record's own stated design (limitations "deliberately longer, more specific and more damaging than
  anything a referee is likely to produce unaided — do not soften it, its length is the point")?
* **Options:** (a) accept as is — the record's stated intent is that limitations be *sharper*, not
  *longer*, than the capability sections; (b) trim the duplication between the capability blocks and
  "Additional capabilities", which is where most of the extra capability space actually goes;
  (c) expand limitations to match. Note this decision constrains escalation 3: choosing (d) there
  makes (b)/(c) here more pressing.

---

## Local roughness left in place (deliberate, per the surgical-editing rule)

1. **Limitation 5 still ends "None is fatal; all are ten-minute discoveries."** The newly added
   1,000,000-state sequential ceiling is arguably more than a ten-minute discovery for anyone wanting
   exhaustive 11-mers. The closing sentence was left alone rather than rewritten.
2. **§Adjudicated disagreements still opens "All 33 values agree between the extraction memo, the
   adversarial review, and this final record."** The review covered 30 values, not 33 (verified,
   `reviews/poolparty.md:14,394-397`). The audits flagged the tally in the header (fixed, APPLIED 8)
   but not this sentence, so it was left untouched.
3. **`composable_operations` now opens "Most operations are `Pool(s) → Pool`" while the row heading
   and value remain `yes`.** The value is untouched and still correct — the qualification is about
   the module-function/method surface, not about composability.
4. **`lazy_generation` keeps "0.043 s" for the 24,000-barcode construction** (re-measured 0.039 s
   here) while the barcode *count* in the adjacent infeasibility example was changed to `N`. The two
   numbers came from the same paragraph but only one is unreproducible.
5. **`design_visualization`'s corrected matplotlib parenthetical now sits next to the unchanged
   sentence "No sequence logos, no rendered DAG image, no interactive viewer,"** which remains true —
   the committed SVG is a hand-drawn figure asset that happens to embed a Matplotlib producer tag,
   not plotting code.

---

## Files touched

| Path | Action |
|---|---|
| `25_poolparty_private/revision/tool_survey/final/poolparty.md` | edited in place — 24 surgical edits across 20 applied findings |
| `25_poolparty_private/revision/tool_survey/fixes/poolparty.md` | created (this change log) |
| `poolparty-statecounter/` | **not modified** — read-only throughout; all execution from `/tmp` with `PYTHONDONTWRITEBYTECODE=1` |

**Repository non-modification, checked rather than asserted:** `git rev-parse --short HEAD` is still
`1bb0179`. `git status --porcelain` shows ` M poolparty/README.md` and two untracked
`.docs_buildhtml/` directories, all of which **predate this session** (mtimes 2026-04-09, 2026-04-04,
2026-04-06). `find . -newermt '2026-08-14 00:00'` outside `.git`, `.venv`, `__pycache__` and
`.pytest_cache` returns **nothing**; the two cache paths that did change are timestamped 03:29 and
03:36, i.e. the earlier fact-check audit's session, not this one (this session began ~12:00).
