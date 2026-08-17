# SeqPro — repair change log

Record repaired: `revision/tool_survey/final/seqpro.md` (edited in place).
Audits processed: `citation_audit/seqpro.md` (19 findings), `factcheck/seqpro.md` (A×8, B×10, C×1).
Repair date: 2026-08-14.

**No capability value was changed.** Where a verified correction touched the evidence under a value,
the value's stated rationale was re-checked and found to survive; anything doubtful is escalated below.

## How verification was done

Primary source only. The repository tarball for the exact audited commit
`63a843985d96dd3f5a7bc8cc20e8bd03f1dabdd9` (0.22.0, 2026-07-27) was downloaded and read directly
(no clone into the project, no install, no execution of SeqPro). In addition:

- GitHub recursive-tree API for that commit (`?recursive=1`), commits API (July 2026 window and
  per-author), contributors API, releases API, repo metadata.
- PyPI JSON (`pypi.org/pypi/seqpro/json`) and the published wheel
  `seqpro-0.22.0-cp39-abi3-manylinux_2_17_x86_64...whl`, opened as a zip and listed (not installed).
- Rendered docs page `https://ml4gland.github.io/SeqPro/ragged/`.
- Python import semantics tested locally on a synthetic package and then directly against the
  extracted wheel tree via `importlib.machinery.PathFinder`.
- Signature behaviour of the README quick-start calls reproduced with stub functions carrying the
  real signatures.

---

## APPLIED (verified, corrected)

| # | Finding (source) | Verified how | Correction |
|---|---|---|---|
| 1 | Citation 1 / Factcheck A2 — `seqpro.experimental` *is* importable | Wheel ships `seqpro/experimental/_experimental.py` + `_visualizers.py` and no `__init__.py`; `PathFinder.find_spec('experimental', [<wheel>/seqpro])` returns a namespace spec; loading `_experimental` succeeds and `edit_distance("AC","AT") == 1`; `_experimental.py` has zero imports | `barcode_generation` entry: "not an importable package" → "not a declared subpackage but an implicit namespace directory — undocumented, unexported and unsupported — though the 0.22.0 wheel does ship both modules and `seqpro.experimental._experimental` still imports". Parallel fixes in `design_visualization` and the limitations bullet |
| 2 | Citation 2 — "190 blobs" | GitHub recursive tree at the pinned commit: `truncated: false`, 190 entries = **162 blobs + 28 trees**; 27 `.py` files under `python/seqpro/` confirmed | Sources table and `mixed_mutagenesis_one_pool`: "190 blobs" → "190 entries = 162 blobs + 28 trees" / "190 entries / 162 blobs" |
| 3 | Citation 3 — "608 lines / 55 releases" of `CHANGELOG.md` | `wc -l` = 608; `grep -c '^## '` = **51**; PyPI has **55** release keys | Changelog counts → "51 release headings" / "51 release sections" (3 places). PyPI "55 releases" statements left untouched — those are correct |
| 4 | Citation 4 / Factcheck A1 / B1 — "exactly one ragged axis" | `rag/_core.py`: `from_offsets` rejects only `shape.count(None) >= 3`; `from_lengths(data, (outer, inner))` builds a two-`None` shape; `from_fields` docstring admits R=2 fields; CHANGELOG "validate R=2 nested ragged layouts; cap at R<=2"; dedicated R=2 test modules | Both occurrences: "exactly one ragged axis" → "one or two ragged axes (R≥3 rejected)" |
| 5 | Citation 5 / Factcheck A4 — `sp.AA.translate()` "DNA/RNA→AA" | `alphabets/__init__.py`: built-in `AA` is built from `canonical_codons_to_aas`, a pure-`T` DNA table; `translate` docstring defines the unknown policy over bytes outside `{A,C,G,T,a,c,g,t}` | `assay_dms`: "translates DNA/RNA→AA" → "translates DNA→AA (the built-in `AA` table is the DNA codon table, so RNA codons containing `U` fall to the `unknown` policy)" |
| 6 | Citation 6 / Factcheck A7 — "0.22.0 and 0.21.x carry BREAKING CHANGE entries" | Only one `### BREAKING CHANGE` heading in `CHANGELOG.md`, under 0.22.0; 0.21.0/0.21.1/0.21.2 contain only Fix/Feat/Perf | "0.22.0 and 0.21.x carry **BREAKING CHANGE** entries" → "0.22.0 carries a **BREAKING CHANGE** entry … (the only such entry in the changelog)" |
| 7 | Citation 7 — "Two of the paper's authors (Laub, Klie) were committing … July 2026" | Commits API for 2026-07-01→07-31: 17 commits, all `d-laub` / `David Laub` plus `github-actions[bot]`. `?author=adamklie` newest commit **2023-10-31** | → "One of the two named authors, David Laub (UCSD), was committing to this repo in July 2026 (Klie's last commits date to 2023-10-31); both are plausible referees" |
| 8 | Citation 8 — `docs/ragged.md` example inventory | `to_padded` / `to_packed` appear **0 times** in `docs/ragged.md` and 0 times on the rendered page; the record-layout example at lines 133–146 uses **`ak.zip`** ("ak.zip returns a Ragged automatically"); `rag.zip` appears 0 times | Vignette 2: dropped `to_padded`/`to_packed`; "record layouts via `rag.zip`" → "record layouts (built in the guide with `ak.zip`, which returns a `Ragged`)" |
| 9 | Citation 9 — "grep for `motif` … returns 0 hits" | `grep -rni motif python/` returns exactly one hit: `experimental/_experimental.py:134` `# Performing motif anal`. CHANGELOG: 0 hits (record correct there) | `combinatorial_motif_place`: "returns 0 hits" → "returns exactly one hit — the unfinished comment `# Performing motif anal` in `experimental/_experimental.py`, no code" |
| 10 | Citation 10 — all 6 changelog `codon` hits are "performance work" | 6 hits confirmed; two are not performance entries — line 280 "single source of truth for codon LUT index" and line 513 "test full codon table"; the 179× entry is real at line 285 | Both occurrences: "translate-kernel performance work" → "forward-translation kernel/LUT work (…, plus LUT-indexing and codon-table-test entries)" / "(kernels, LUT indexing, a codon-table test, a benchmark)" |
| 11 | Citation 11 / Factcheck A5 — schemas cover "open vs closed" | `_coords.py:10–27`: `CoordSchema` fields are `chrom, start, end, zero_based, strand` — no closure field. `set_schema` (100–143) renames columns and sets `coordinate_system_zero_based` metadata; no arithmetic on start/end | `genome_coordinates`: → "each carrying a 0- vs 1-based flag (there is no open/closed field, and `set_schema` renames the coordinate columns and tags the frame rather than converting values)". Same trim in the notable-capabilities bullet |
| 12 | Citation 12 — "From `README.md` / `docs/index.md`, verbatim" | Both quoted passages are verbatim `README.md:10` and `README.md:12` (links rendered away). `docs/index.md:6,8` is reworded ("It makes almost zero compromises on speed…", "All functions accept strings, …") | → "From `README.md`, verbatim:" |
| 13 | Citation 13 — API pages are "signatures only" | `zensical.toml:67` and `docs/api/index.md:21` set `show_source = true`; `docstring_style = "numpy"`, so parameter/return docs render. No `Examples` section exists in any `python/seqpro` docstring | → "mkdocstrings-rendered signatures, docstrings and source (`show_source: true`), with no narrative examples" |
| 14 | Citation 14 — "README ends with a bare '# More to come!'" | `README.md:81` is the heading; the file continues to line 83 with a contributions / code-of-conduct paragraph | → "README's final section heading is a bare **\"# More to come!\"** (followed only by a contributions/code-of-conduct paragraph)" |
| 15 | Citation 15 — both visualizers are "histograms" | `_visualizers.py:16` `ax.hist`; `_visualizers.py:26` `ax.plot(nuc_contents.T)` — a positional line plot | `design_visualization`: "QC histograms" → "QC plots … — a GC histogram and a positional nucleotide-frequency line plot" |
| 16 | Citation 16 — README makes sibling scope "explicit" | `README.md:10` says only "heavily utilized by other packages including …"; no exclusion statement anywhere in the README | Limitations bullet: "are explicitly out of scope" → "sit with those packages, a division the README implies rather than states" |
| 17 | Citation 19 — "basedpyright + pyrefly strict type checking" | `pyproject.toml:66–67` `[tool.pyrefly] preset = "strict"`; `pixi.toml:25` `typecheck = "pyrefly check python"`; `.pre-commit-config.yaml:24–29` runs `pixi run typecheck`. `[tool.basedpyright]` sets only `include`/`enableTypeIgnoreComments`/`reportMissingTypeArgument`; no task or workflow invokes it | → "pyrefly strict type checking wired into pre-commit (`pixi run typecheck`) alongside a configured basedpyright" |
| 18 | Factcheck A6 — plotting symbols do exist in the tree | `experimental/_visualizers.py` imports matplotlib and defines `plot_gc_content`, `plot_nucleotide_content` | Scope paragraph: "nor anywhere in the 27 Python source files" → "and — apart from the two unmaintained plotting helpers in `experimental/_visualizers.py` (see `design_visualization` below) — none anywhere in the 27 Python source files" |
| 19 | Factcheck B3 — the canonical quick start contains two non-running calls | `_modifiers.py:39–46` — `alphabet` is a required positional of `k_shuffle`; `_modifiers.py:144–149` — `jitter_axes` is a required keyword-only arg. Reproduced with stub signatures: both README calls raise `TypeError: missing 1 required … argument` | Vignette 1, one appended sentence: "As written, two of those calls do not run: `k_shuffle` omits the required `alphabet` argument and `jitter` omits the required `jitter_axes`, so both raise `TypeError`." |
| 20 | Factcheck B4 — synchronized multi-array jitter omitted | `_modifiers.py:144–207`: `jitter(*arrays, …)`, docstring "using the same jitter across arrays"; output length `L - 2*max_jitter` | `assay_insilico`, one clause: "`jitter` (which applies the same random crop across several aligned arrays, keeping a sequence and its per-base tracks in register)" |

## REJECTED (finding not confirmed, or not an error in the record)

| # | Finding (source) | Why rejected | Evidence |
|---|---|---|---|
| R1 | Citation 17 — "canonical"/"standard" dinucleotide-shuffle null is uncited | Status is "uncited", not "wrong"; the statement is accurate field knowledge and the repo itself treats it as the standard comparator. Removing or hedging a true statement would degrade the record; citation policy for background claims is an author decision | `pixi.toml:85–90` adds `tangermeme` and `dinuc-shuf` as benchmark dependencies — SeqPro benchmarks `k_shuffle` against the field's reference dinucleotide-shuffle implementations |
| R2 | Citation 18 — "edgeR-style" TMM is uncited | The characterization is correct, not merely unsupported. `transforms/tmm.py` implements TMM with edgeR `calcNormFactors` defaults: `log_ratio_trim=0.3`, `expression_trim=0.05`, `apply_weighting=True`, `expression_cutoff=-1e10`, 0.75-quantile reference-sample selection, log-ratio + absolute-expression trimming, inverse-variance weighting | `transforms/tmm.py:8–63`, `:162` (`@nb.njit(parallel=True, nogil=True, cache=True)`) |
| R3 | Factcheck A3 — "`sp.Ragged` is not a root API name" | Verified true of the *package* (`Ragged` is absent from root `__all__`), but it is not an error in the record: the string `sp.Ragged` appears only inside the verbatim `SKILL.md` quotation; the record's own prose consistently writes `sp.rag.Ragged`. Annotating a quotation is a framing decision — escalated as E1 | `python/seqpro/__init__.py:13–41` (28 names, no `Ragged`); `SKILL.md:15,52`; record lines 60 (quote) vs. 97, 401 (record's own prose) |
| R4 | Factcheck A8 — the "All functions" quote is unqualified | Same class as R3. The record presents it as an attributed verbatim README quotation and makes no such claim in its own voice. The underlying README statement *is* imprecise (verified: `pad_seqs` and `k_shuffle` have no `Ragged` overload), but adding a caveat inside a quoted block is a framing decision — escalated as E1 | `README.md:12`; `_encoders.py:24` (`pad_seqs`, dense only), `_modifiers.py:39` (`k_shuffle`, `SeqType` only) |
| R5 | Factcheck B9 — omission of `seq2kmer` / `kmer2seq` | Verified they exist (`experimental/_experimental.py:137–170`), but they are unsupported private helpers in a module the record already characterises, and they bear on no capability row. Adding them would be inventory padding, not repair | `experimental/_experimental.py` — both are pure string helpers, unexported, undocumented |

## ESCALATED (human decision required; record left untouched)

**E1 — Should verbatim author-doc quotations carry inline corrections?**
Two audit findings (Factcheck A3 and A8) are not errors in the record but errors *inside quotations the
record reproduces faithfully*: SKILL.md's `sp.Ragged` (the root package exports no `Ragged`; the correct
path is `sp.rag.Ragged`) and the README's "All functions in SeqPro take as input …" (acceptance is
function-specific; `pad_seqs` and `k_shuffle` have no `Ragged` overload). Both verified against
`python/seqpro/__init__.py`, `_encoders.py` and `_modifiers.py`.
*Options:* (a) leave the quotations bare, as now; (b) add a one-line bracketed editor's note after each
quote; (c) add a single footnote under "What SeqPro actually is" noting that both author documents
overstate the uniformity of the API. Not decided here — annotating quoted source is an editorial call,
and one of the quotes is the record's headline scope statement.

**E2 — How much API detail should the record carry? (Factcheck B2, B5, B6, B7, B8, B10)**
Six further "incomplete" findings were checked and are true statements about the code: `k_shuffle`'s full
call contract and seed/batch-size reproducibility semantics; padding/truncation modes and all-zero OHE of
out-of-alphabet characters; `rag.reverse_complement`'s per-row mask and `copy=False` in-place mode;
advanced `Ragged` selection/conversion surface (`empty`, `from_offsets`, fancy indexing, `to_numpy`,
method- vs module-level `to_padded`); custom-alphabet construction constraints (complement must be the
reverse of the alphabet, equal-width codons, single-character amino symbols); and peripheral limits
(xarray is an **undeclared** optional dependency — absent from `pyproject.toml` and `pixi.toml`, with
`xr/__init__.py:10–13` raising `ImportError`; `xr.translate` is defined at line 107 but absent from
`xr.__all__`; an xarray `jitter` is fully commented out from `xr/__init__.py:163`; and `bed.with_len`
has no clamp, so it can emit negative starts).
None of them changes a capability value, and adding all six would roughly double the record's API prose.
*Options:* (a) add none — the record is a capability memo, not an API reference; (b) add only the two with
referee-facing bite (undeclared xarray dependency; `with_len` produces negative starts); (c) add all six as
one-clause notes. Whether these omissions are "material" is exactly the judgment reserved for the authors.

**E3 — Factcheck section C (balance).**
The fact-check argues the record is "limitation-heavy rather than balanced": that limitations are documented
in substantially greater depth than capabilities, and that the (now corrected) one-ragged-axis error had been
suppressing a major current strength. This is an emphasis judgment about a memo whose stated purpose is to
defend a column of "no"s to an author-referee, so it is not something a repair pass should re-decide.
*Options:* (a) accept the current proportions; (b) promote the "Notable capabilities" list above the
capability assessment; (c) expand Block A/B positive evidence. Note that fix #4 has already restored the
R=2 Ragged strength, which is the one concrete instance the balance finding cites.

## Local roughness left in place (per the surgical-editing rule)

- `barcode_generation` still labels its paragraph "**Strengthened vs. the original extraction:**", but after
  fix #1 that paragraph now partly *retracts* rather than strengthens (the directory is still not a declared
  subpackage, yet the module does import). The label was left untouched; a one-word change to "Refined" would
  resolve it if the authors want it.
- `mixed_mutagenesis_one_pool` still describes `jitter` as a bare "random offset crop" while `assay_insilico`
  now records its synchronized multi-array behaviour (fix #20). Both statements are true; the terser one was
  left alone.
- The flagged-judgment list (item 3) still says the original extraction's blanket laziness sentence "has been
  removed" and item 4 still summarises Blocks B–D as verified four ways. Both re-checked and still accurate
  after the changelog-count corrections.

## Re-verified and left unchanged (record was right)

`__all__` = 28 names; 27 Python files; `gtf.py` 1,813 bytes; `bed.py` 378 lines; `tmm.py` 205 lines; 49
`docs/superpowers/` markdown files; 30 `test*.py` modules; 10 `.rs` files; `docs/api/` contains exactly
`index, alphabets, bed, gtf, ragged, types` and `index.md` lists 14 members; `transforms/__all__` =
`["TMM","Jitter","KShuffle","Random","ReverseComplement","Sequential"]` with `Tokenize` defined at
`augmentation.py:114` but unexported; `xr` `dask="parallelized"` at lines 48/98/152 with
`__all__ = ["bin_coverage","ohe"]`; `bed.read` dispatches to eager `pl.read_csv` at 236/287/356;
`with_len` summit recentering exactly as described; `remove_whitespace` body is `pass`;
`_analyzers._count_kmers` is `raise NotImplementedError` under `# TODO: non-trivial to parallelize/SIMD`;
`count_kmers_seq` and `normalize_coverage` absent from `__all__`; `random_seqs` returns
`seed.choice(alphabet.array, size=shape)`; 64-entry codon LUT gated by `_can_build_lut`;
`_translate_stop_ends` at `_alphabets.py:648`; 0 hits for `hgvs` anywhere; MIT licence, not archived,
14 stars, 1 fork, 0 open issues, last commit `2026-07-27T19:02:07Z`, `pushed_at 2026-07-27T19:07:56Z`,
10 GitHub releases 0.16.0→0.22.0 (2026-06-14 → 2026-07-27), 55 PyPI releases since `2023-03-26T17:33:00`,
`requires-python >=3.10`, dependency list matches `pyproject.toml` including
`polars-config-meta[polars]>=0.3.2`.

## Pass 2 — policy application

**Date:** 2026-08-14

**Baseline and counts.** Every factual item was rechecked against the official SeqPro repository at
`63a843985d96dd3f5a7bc8cc20e8bd03f1dabdd9`, its documentation, or PyPI. Counts use the 11 atomic
policy items below: **4 applied · 7 declined-by-policy · 0 rejected-unverifiable · 0 still
escalated.** No capability value changed, and no Table 1 or Table 2 file was edited.

| Open item | Outcome | Edit and primary-source verification |
|---|---|---|
| **E1 / Factcheck A3+A8 — author-documentation errors in quotations** | **applied** | Kept both quotations verbatim and added one short caveat in the existing scope entry: `Ragged` is reached as `sp.rag.Ragged`, not root `sp.Ragged`, and accepted input representations are function-specific. Defined **root export** as a name in `seqpro.__all__`, with submodule surfaces treated separately. Verified in `python/seqpro/__init__.py`, `rag/__init__.py`, `_encoders.py:24–123`, and `_modifiers.py:39–108`; `pad_seqs` and `k_shuffle` have no `Ragged` overload. |
| **E2 / B2 — full `k_shuffle` call and reproducibility contract** | **declined-by-policy** | Verified the mandatory alphabet, explicit axes, fixed-seed per-row derivation, batch-size dependence, and `k >= length` behavior in `_modifiers.py:39–108` and `src/kshuffle.rs:135–203`. These details do not alter Table 1's grouped Purpose, transformation Key-features, transformed-sequence Output, or PyPI Availability cells. |
| **E2 / B5 — padding/truncation and unknown-symbol encoding** | **declined-by-policy** | Verified left/right/both behavior, odd-change placement, zero OHE padding, dense `pad_seqs` restriction, and all-zero unknown OHE behavior in `_encoders.py:24–228`. These are lower-level transformation semantics and do not correct a Table 1 cell. |
| **E2 / B6 — Ragged reverse-complement modes** | **declined-by-policy** | Verified the per-row Boolean mask, `copy=False` mutation, S1/final-axis restriction, and 256-entry LUT contract in `rag/_ops.py:22–100`. The grouped Table 1 Key-features cell already accurately states transformation support; these modes do not change it. |
| **E2 / B7 — advanced Ragged selection/conversion surface** | **declined-by-policy** | Verified `empty`, `from_offsets`, indexed/masked/fancy selection, `to_numpy`, and method/module packing and padding behavior in `rag/_core.py` and `rag/_ops.py:113–338`. The record already credits one- and two-axis `Ragged`, packing/padding, fields, and interop; the omitted method inventory changes no Table 1 cell. |
| **E2 / B8 — custom-alphabet construction restrictions** | **declined-by-policy** | Verified complement-equals-reverse validation, equal-width codons, single-character amino symbols, the generic nonstandard-codon path, and the sparse-table `unknown="drop"` caveat in `alphabets/_alphabets.py:25–77, 297–361, 527–568`. These restrictions do not change Table 1's sequence-manipulation Purpose or broad Key-features description. |
| **E2 / B10 — XArray and interval edge limits** | **declined-by-policy** | Verified that XArray is absent from `pyproject.toml` and `pixi.toml`, `seqpro.xr` raises a targeted `ImportError` when absent, `translate` is outside `xr.__all__`, XArray jitter is commented out, and `bed.with_len` does not clamp coordinates (`xr/__init__.py:1–228`, `bed.py:63–99`). The default PyPI package remains available as Table 1 states, and the grouped Output cell does not promise bounded genomic intervals; no clause qualifies for shipping. |
| **E3 / fact-check section C — balance and emphasis** | **declined-by-policy** | Policy A expressly declines rebalancing or proportional-emphasis edits. No record text changed for this item. |
| **Citation audit 17 — “standard/canonical” shuffle-null assertions** | **applied** | Removed `standard` and `canonical null/background`; retained only the repository-verifiable facts that `k_shuffle` preserves k-lets, uses `src/kshuffle.rs`, and is benchmarked. No SeqPro repository/docs/package source or SeqPro paper establishes the field-wide practice claim; there is no SeqPro paper. |
| **Citation audit 18 — “edgeR-style” TMM equivalence** | **applied** | Changed the existing heading to **TMM count normalization** and retained the directly verified estimator mechanics. `transforms/tmm.py` contains no edgeR mention or citation, so the external equivalence label was dropped rather than inferred. |
| **Policy C — provenance, anecdotes, and ambiguous labels** | **applied** | Removed the organization/sibling-repository and prior-analysis rows from Sources; removed sibling-package responsibility claims, cross-tool superlatives, and the unsupported torchvision analogy; replaced the sibling-notebook entry in place with the official SeqPro-tree fact that it contains no `.ipynb`; dropped the plausible-referee anecdote; clarified that `Klie2023kg` is not SeqPro evidence; defined root exports and `maintained`; described `experimental` consistently as an importable implicit namespace package outside the declared public surface; and removed the stale no-execution claim. Primary checks were the official tree, README/docs, package source/wheel, GitHub metadata, and PyPI. |

### Version drift (Policy D)

No edit was required. Read-only `git ls-remote` on 2026-08-14 returned
`63a843985d96dd3f5a7bc8cc20e8bd03f1dabdd9` for both `refs/heads/main` and `refs/tags/0.22.0`;
PyPI still reports **0.22.0**, 55 release keys, Python `>=3.10`, and the 2026-07-27 upload. The GitHub
API still reports the same `pushed_at`, MIT license, archive state, stars/forks/issues, and the same
0.16.0→0.22.0 ten-release sequence. The header therefore governs without a drift parenthetical.

### Value-basis check

The quotation and provenance corrections do not weaken a capability value. `assay_insilico = partial`
still rests on SeqPro's own preprocessing operations and absence of model scoring/design;
`genome_coordinates = partial` still rests on official BED/schema support plus the absence of a
reference-sequence reader; and `combinatorial_motif_place = no` still rests on the complete root and
submodule tree. TMM and shuffle wording lost only unsupported external equivalence labels. No value
escalation is required.

### Rejections, escalations, and row-substitution report

- **Rejected-unverifiable findings:** none. The two unverifiable *underlying assertions* were the
  field-wide shuffle-null label and edgeR equivalence; Policy C required removing them, and the audit
  finding that SeqPro's admitted primary sources do not support them was verified.
- **Remaining escalations:** none. E1–E3 are resolved, and no value basis failed.
- **Uniform-row evidence:** none can be established from this SeqPro-only pass. The locked reference
  already flags `Codon / amino-acid substitutions` as potentially near-uniform and `Recombination and
  chimeras` as potentially one-positive, but their eight shipping columns remain unscored here.
- **Candidate considered, report only:** replace `Codon / amino-acid substitutions` with
  **Variable-length/nested sequence-batch representation**. SeqPro's official source supplies unusually
  concrete evidence: R=1/R=2 `Ragged`, record fields, packing/padding, ufunc dispatch, and awkward/NumPy
  interop (`rag/_core.py`, `rag/_ops.py`, `docs/ragged.md`). It is **not nominated for substitution**:
  SeqPro has no Table 2 column, the locked table compares eight library-design tools, and this pass has
  no primary-source cross-tool values proving discrimination. No locked row was changed.

### Neighbouring tensions left intact

The historical Pass-1 `## ESCALATED` section still says E1–E3 were untouched; this Pass-2 section
supersedes those outcomes without rewriting the audit trail. `barcode_generation` still says
**Strengthened** although the namespace-package correction is partly a retraction. Table 1's grouped
`no library abstraction` wording remains consistent with `Ragged` being a general-purpose batch
container rather than a library specification, but its broad `Transformed sequences` output cell does
not enumerate interval tables, arrays, or TMM-normalized matrices; that compression is logged, not
edited.
