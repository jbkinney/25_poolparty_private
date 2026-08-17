# v2 capability record — Oligopool Calculator

**Slug:** `oligopoolcalc`
**Full name:** Oligopool Calculator (`oligopool` Python package, v2026.2.22.1)
**Citation key:** `Hossain2024oc`
**Tier:** 3
**Row set:** v2 (`ROWS_v2.md`) — 33 rows. Supersedes the 24-row v1 scoring in `final/oligopoolcalc.md`.
**Basis:** `final/oligopoolcalc.md` read in full (the adversarially-reviewed, thrice-verified record), plus targeted re-greps of the cached raw sources for the seven NEW rows only (`README.md` 16,084 B, `docs/docs.md` 98,900 B, `docs/api.md` 76,438 B, `docs/agent-skills.md` 23,019 B, `examples/library-compressor/mutant_generator.py` 8,431 B). No re-extraction was performed.

**Scoring stance:** this is a purpose-built tool in the oligopool-design space and is scored on its own terms. Referees may include Hossain or Salis. Where the v2 wording genuinely supports a stronger value than v1's coarser row, the stronger value is given (`library_first_class_object`, `exhaustive_single_scans`, `higher_order_combinatorial`, `degenerate_iupac_codons`, `readout_analysis`). Where it does not, it is not.

**Block B caption rule (carried from v1, must appear in any published table):** `yes` = demonstrated end-to-end on published libraries of that assay type, with shipped tooling for it; `partial` = documented recipe only, with the assay-specific variant content supplied by the user.

---

## Block A — library specification

### `library_first_class_object` = **yes**  *(v1 `library_as_object` was `partial`; the split resolves it upward)*
*Evidence:* Every one of the 22 modules has the signature `(input_data, ..., output_file=None) -> (DataFrame, stats)` — a pandas DataFrame goes in, a DataFrame comes out, and `output_file` is optional, so nothing need touch disk. docs.md "Core Concepts → The DataFrame Flow", verbatim: *"**Input**: CSV path or pandas DataFrame with an `ID` column … **Chainable**: Output of one module feeds into the next."* The user holds the whole library in memory, inspects it, transforms it, and passes it onward; every module operates on the whole table at once, so the user never loops over sequences. This is emphatically **not** a tool that only writes files. The honest qualifier: it is a *generic* DataFrame, not a dedicated library class carrying design semantics (no library class exists anywhere in `oligopool/**/*.py`), the "library" is the accumulated set of named columns, and `final` reduces it to `ID`/`CompleteOligo`/`OligoLength` only (api.md L613), after which docs.md L764 warns verbatim: *"Save your annotated design DataFrame/CSV before `final` if you plan to `index`, `pack`, and count later."*
*Source:* docs.md "The DataFrame Flow", L764; api.md L613 and all 22 module signatures.
*v1→v2:* v1 scored `partial` on the conflated row. The conflation was precisely between "can you hold a library" (yes) and "can you combine libraries" (no) — which is what the split exists to separate. Holding/inspecting/transforming/passing-onward is unambiguous here.

### `composable_operations` = **yes**  *(rename of `dag_chaining`; value carried)*
*Evidence:* Two shipped mechanisms. (a) CLI YAML pipelines with an explicit DAG — docs.md §"Parallel Pipeline Execution" (L1630 ff.), verbatim: *"**Serial format**: a list of command names … **Parallel/DAG format**: a list of step objects with `name`, `command`, and optional `after`, which defines dependencies explicitly. **Mixed format (supported)**: you can mix strings and step objects in one list; if any dict-style step appears, DAG parsing is used…"*; execution model: *"Steps with no `after` dependencies form level 1 (eligible to run concurrently). Each subsequent level waits for dependencies … `--dry-run` shows execution levels and parallelism."* (b) A runnable **design-mode** DAG, `examples/cli-yaml-pipeline/mpra_design_parallel.yaml`: `primer` → (`barcode_a` ∥ `barcode_b`, both `after: [primer]`) → `rejoin` (`command: join`) → `spacer` → `final`. The pipeline is not fixed by the tool: modules may be ordered and re-invoked freely.
*Caveat to keep verbatim:* the DAG is an **execution/scheduling graph over design steps applied to one flat table**, not a **content-composition graph over library sub-populations** (`join` requires *"`ID` sets must match exactly across both inputs"* and *"never creates or drops rows; mismatched IDs are an error"*, api.md L563/L567 — re-verified). There is also no *nesting* of operations. The authors concede the design-side limit: docs.md L1634 *"design workflows are often dependency-dense and mostly sequential"*; design branch-join is called **"Rare"**. Two documented ordering constraints exist: `compress` is a *"terminal branch in design mode"* (docs.md L2395) and `pad` must be run per fragment (api.md L749).
*Source:* docs.md L1630–L1700, L1634, L2395; api.md L563, L567, L749; `examples/cli-yaml-pipeline/mpra_design_parallel.yaml`.

### `lazy_generation` = **no**  *(rename of `lazy_evaluation`; value carried)*
*Evidence:* All 22 modules return materialized DataFrames/CSVs; the evaluation model is eager throughout. Documentation-wide grep for `lazy|generator|yield|iterator|on-demand` returns no feature hits. The one enumerating module, `expand`, is explicitly eager with a safety cap — api.md `expand` Notes, verbatim: *"Expansion can be exponential; use `expansion_limit` for safety with highly degenerate sequences"* (re-verified at source this pass).
*Completeness note (state before a referee does):* the shipped example `examples/library-compressor/mutant_generator.py` **is** generator-based — re-verified this pass: `generate_single_mutants`, `generate_codon_variants`, `generate_multi_position_variants` are all typed `Generator[...]` with `yield`. But it is example code feeding a DataFrame *into* the package, not the package's evaluation model.
*Source:* api.md (all module Returns lines, `expand` Notes); documentation-wide grep; `mutant_generator.py` L43–L172.

### `library_algebra` = **no**
*Evidence:* **The load-bearing structural fact, verified and unrebuttable: across all 22 modules there is no row-union / sub-library concatenation operator, and no library-level sample or repeat operator.** `merge` "Concatenate[s] contiguous **columns** into a single column" (api.md L423–L427, re-verified verbatim this pass: *"Purpose: Concatenate contiguous columns into a single column"*). `join` requires *"`ID` sets must match exactly across both inputs (order can differ)"* and *"`join` never creates or drops rows; mismatched IDs are an error"* (api.md L556–L567, re-verified verbatim). **An Oligopool Calculator library is therefore always exactly the rows the user supplied.** Combining two sub-libraries requires an external `pd.concat`. No sampling operator and no replicate/repeat operator exists on any module.
*Honest counter-note:* two modules do change row counts — `compress` collapses variants into fewer degenerate rows and `expand` enumerates degenerate rows into more concrete rows. These are compression/decompression of one library, not algebra over two libraries, and `compress` output is documented as a terminal branch (docs.md L2395).
*Source:* api.md L423–L427 (`merge`), L556–L567 (`join`), `compress`/`expand` sections; docs.md L729, L2395.

### `exhaustive_single_scans` = **partial**
*Evidence:* **Not in the shipped API; present in shipped author-written example code.** None of the 22 `op.*` modules generates variants — variant content is a required user input (docs.md Quick Start *"# Start with your variants"*; docs.md L1304 *"Start with your variants in a CSV (must have an `ID` column)"*; agent-skills.md L371 lists *"generate all substitutions → `compress` → order `synthesis_df`"* where "generate all substitutions" is the user's step; the tutorial notebook: *"designing exact ribozymes is beyond the scope here"*). **However**, the repo ships `examples/library-compressor/mutant_generator.py` (8,431 B, header `Author: Ayaan Hossain`) whose `generate_single_mutants(sequence, include_wildtype=True)` enumerates **all 3L single-nucleotide substitutions** at every position, auto-minting HGVS-style IDs (`variant_id = f'{ref_base}{pos+1}{alt_base}'` → `A5T`) and optionally yielding the WT first. Re-verified at source this pass (L43–L79).
*Scope of the `partial`:* substitutions only — **no deletion scan, no insertion scan, no amino-acid-level scan** (there is no translation table or reading frame anywhere in the package). And it is example code, not `op.*` API.
*Source:* `mutant_generator.py` L43–L79; docs.md Quick Start, L1304; agent-skills.md L371; api.md (22 module signatures).
*Judgement:* `no` is defensible on a strict shipped-API reading and is what v1's coarser row implied. `partial` is given because the row asks about the capability, the code is in the repo, and the fairness requirement forbids understating. **Flagged — must be scored consistently with how example/vignette code was treated for the other 11 tools.**

### `sampled_random_mutagenesis` = **no**
*Evidence:* Nothing anywhere samples variants at a rate or at random. No mutation-rate, sampling-fraction, or library-size-sampling parameter exists on any of the 22 module signatures. The `random_seed` parameter present on `barcode`/`primer`/`motif`/`spacer`/`split`/`pad`/`compress` seeds the **stochastic design of scaffold elements** (and `compress`'s Monte-Carlo rollouts, api.md `compress`: *"`rollout_simulations` … Monte Carlo simulations per decision"*), not variant sampling. The shipped example `mutant_generator.py` is exhaustive-only: its three generators enumerate 3L, 64, and 4^N variants respectively; its only random function is `generate_random_sequence(length, seed)`, which makes a starting sequence, not a variant sample. Re-verified at source this pass.
*Source:* api.md all 22 signatures, `compress`; docs.md §"Reproducibility" L279–L289; `mutant_generator.py` L27, L43–L172.

### `higher_order_combinatorial` = **partial**
*Evidence:* Same structure as `exhaustive_single_scans` — not in the shipped API, present in shipped author-written example code. `mutant_generator.py::generate_multi_position_variants(sequence, positions, include_wildtype=True)` enumerates **all `4**len(positions)` combinations** across a user-chosen set of positions in the same sequence, IDs `'_'.join(f'p{p+1}{b}' ...)` → `p1A_p5C`. `generate_codon_variants` similarly enumerates all 64 nucleotide triplets at one codon position (i.e. up to 3 simultaneous changes within a triplet). Re-verified at source this pass (L82–L172).
*The limits:* exhaustive-only (no sampled higher-order); combinations are over **positions the user names**, with no ceiling on order but also no order parameter, no distance/linkage model, and no way to mix orders in one call; and it is example code, not `op.*`. Within the package proper, higher-order combinatorial *content* has no representation at all — the only cross-product anywhere in `op.*` is `expand`, which enumerates concrete sequences from an IUPAC string and is documented as *"Primarily used as a verification tool to confirm `compress` output"* (api.md `expand` Notes, re-verified verbatim).
*Source:* `mutant_generator.py` L82–L172; api.md `expand`; api.md 22 module signatures.
*Judgement:* flagged, same consistency caveat as `exhaustive_single_scans`.

### `heterogeneous_components_one_library` = **no**
*Evidence:* Read as the row intends it (structurally different **variant classes** coexisting as rows in one library specification — e.g. exhaustive singles + sampled higher-order + WT replicates), this is `no`, twice over. (i) The package never derives variants, so there is no specification language for variant classes at all. (ii) Even in the example code, `generate_variant_dataframe(sequence, strategy=...)` takes **`strategy` as a single string — `'single'` | `'codon'` | `'multi'`, dispatched by an `if/elif/elif/else` that raises on anything else** (re-verified at source, `mutant_generator.py` L200–L231) — exactly one strategy per call, with no mixed specification expressible. (iii) And two differently-generated sets **cannot be unioned by any shipped operation**, because there is no row-union operator (`merge` concatenates columns; `join` "never creates or drops rows"). A user can of course hand-build a heterogeneous CSV externally and hand it to `oligopool`, but that is composition outside the tool.
*Alternate reading, disclosed:* if the row were read as "structurally different **element types** coexisting in one architecture", Oligopool would score `yes` — a single design table carries variant, primer, barcode, motif, spacer, split and pad columns simultaneously, and `examples/design-assembly-parser` exposes exactly this as a declarative spec with `element_type ∈ {variant, primer, barcode, motif, spacer}`. The v2 row text ("e.g. exhaustive singles + sampled higher-order + WT replicates") is unambiguously about row-wise variant heterogeneity, so `no` stands — but the alternate reading is noted so a referee cannot spring it.
*Source:* `mutant_generator.py` L200–L231; api.md L423–L427, L556–L567; docs.md L1304; `examples/design-assembly-parser`.

### `combinatorial_motif_place` = **no**  *(carried from v1)*
*Evidence:* `motif` inserts **one** motif into **one** named column per call — re-verified verbatim at api.md L230–L234: *"Purpose: Insert sequence motifs (per-variant or constant anchors) with constraint satisfaction"*, with a single required `motif_sequence_constraint` and a single `motif_column`; `motif_type` ∈ {`0`/`'variable'`, `1`/`'constant'`}. **No position, orientation, or motif-set enumeration parameter exists** on the signature. Multiple motifs require multiple chained `motif` calls at fixed column positions; there is no cross-product over motif sets. Paper §"Restriction Sites, Motifs, and Spacers": the pathfinder *"identifies an oligonucleotide-specific sequence solution that satisfies design rule #1"* — one solution per oligo, not a combinatorial expansion. `revcomp` is global, not per-variant orientation enumeration — re-verified verbatim at api.md L480: *"Purpose: Reverse complement a range of columns and reverse their order"*, parameterized only by `left_context_column`/`right_context_column`, i.e. a whole-table column range.
*Source:* api.md L230–L262 (`motif`), L476–L504 (`revcomp`), `expand` Notes; paper §"Restriction Sites, Motifs, and Spacers"; docs.md L970–L977.

### `barcode_generation` = **yes**  *(carried from v1 — flagship)*
*Evidence:* Full verified signature (api.md L76–L120): `op.barcode(input_data, oligo_length_limit, barcode_length, minimum_hamming_distance, maximum_repeat_length, barcode_column, output_file=None, barcode_type=0, left_context_column=None, right_context_column=None, patch_mode=False, cross_barcode_columns=None, minimum_cross_distance=None, excluded_motifs=None, background_directory=None, random_seed=None, verbose=True)`. Constraints: minimum pairwise **Hamming distance**; `maximum_repeat_length` (*"Maximum shared repeat length between barcode and context/oligos"*); `excluded_motifs` (list, CSV/FASTA, DataFrame with `Exmotif` column, comma-string, or named multi-source dict); context/edge-effect awareness (*"At least one of `left_context_column` or `right_context_column` must be specified"*); `background_directory` k-mer screening, *"junction-aware when context columns are provided"*; cross-set orthogonality (`cross_barcode_columns` + `minimum_cross_distance`); `patch_mode` incremental fill; `random_seed`. **Two distinct algorithms** via `barcode_type`: *"`0`/`'terminus'`=fast terminus-optimized … `1`/`'spectrum'`=thorough spectrum-optimized (targets k-mer saturation)."* The barcode is written as a column at its correct position in the oligo — genuinely attached, not merely enumerated. Paper: *"over four million highly unique and compact barcodes in 1.2 h."*
***One honest gap — keep it:*** **GC content is not constrainable for barcodes** (no GC or sequence-constraint argument on `op.barcode`; the paper's four barcode criteria never mention GC). Phrase it as "not constrainable for barcodes", never "absent from the tool" — GC **is** computed and reported for primers (`primer_guanine_cytosine_content`, 6 occurrences in `oligopool/primer.py`).
*Source:* api.md L76–L125; `oligopool/primer.py`; paper ¶"Barcodes"; docs.md L471.

### `per_sequence_provenance` = **partial**  *(carried from v1)*
*Evidence:* Three real per-sequence mechanisms, one real gap. (i) Each designed element lives in its own named column, so a row structurally records which barcode/primer/motif/spacer composes that oligo. (ii) `verify` emits genuine **per-row, per-COLUMN attribution** (api.md L1369–L1405): `LengthConflictDetails.column_lengths` = *"Per-column contribution (gap-stripped lengths) for attribution"*; `ExmotifConflictDetails` carries `motif`, `occurrences`, `library_baseline`, `excess_occurrences`, `positions` (*"0-based start offsets in `CompleteOligo`"*) and `columns` (*"Column attribution for each position"*); `BackgroundConflictDetails` carries the same `matched_kmers`/`positions`/`columns` triple. (iii) `compress` returns `mapping_df` (`ID`, `Sequence`, `DegenerateID`) and `expand(mapping_file=...)` *"restores original variant IDs in output"* — re-verified verbatim at api.md `expand` Optional Parameters this pass — real per-sequence provenance across the compression step.
*What does not exist:* any record of **how a variant was derived** — no mutation history, no design history — because the package never derives variants. The per-run `stats` dict is aggregate (docs.md L176: *"`stats` is **aggregate**: programmatic pass/fail decisions and summary reporting"*).
*v2 wording note:* the v2 row says "beyond naming the mutation". Oligopool does not name the mutation either (IDs are user-supplied), but it does record per-column composition and conflict attribution, which is a different and genuine kind of provenance. `partial` remains the right cell.
*Source:* api.md L1369–L1405, L811, L862, `expand`; docs.md L176.

### `automatic_naming` = **no**  *(carried from v1)*
*Evidence:* **Scoped to library-member (variant) naming, where it is airtight.** The `ID` column is a **required user input** on every module (docs.md L1304: *"Start with your variants in a CSV (must have an `ID` column)"*), and no module renames, enumerates, or invents variant IDs. Output *file* names get auto-suffixes (`.oligopool.<module>.csv`).
*Two things this cell deliberately does NOT claim (both are easy to disprove):* (1) Output **column** names are not all user-specified — five modules auto-create columns: `split` → `Split1..SplitN` (api.md L676), `pad` → `5primeSpacer`/`ForwardPrimer`/`<split_column>`/`ReversePrimer`/`3primeSpacer` (L743), `final` → `CompleteOligo`/`OligoLength` (L613), `expand` → `ExpandedSeq`/`OligoLength` (L868), `compress` → `DegenerateID`/`DegenerateSeq`/`Degeneracy` (L811–L812). (2) `compress` **does** mint identifiers for rows the user never named — `synthesis_df` rows are identified solely by auto-minted `DegenerateID`s. Those name **degenerate groups**, not designed library members, so the value holds. Separately, the example `mutant_generator.py` auto-mints informative HGVS-style variant IDs (`A5T`, `ATG_AAA`, `p1A_p5C`) — again example code, not `op.*`.
*Source:* docs.md L1304; api.md L613, L676, L743, L811–L812, L868; `mutant_generator.py`.

### `design_visualization` = **no**  *(carried from v1)*
*Evidence:* Airtight, verified by positive-absence at source level. All **40** `oligopool/**/*.py` files grepped for `seaborn|matplotlib|pyplot`: **zero hits** (`seaborn>=0.13.2` is a stale declared PyPI dependency that no shipped module imports). Documentation-wide grep for `visuali|plot|matplotlib|graph` returns no feature hits. All inspection is **textual**: `lenstat` prints a per-element length report, `verify` emits a conflict table, `inspect` prints binary-artifact metadata, CLI `--dry-run` prints DAG execution levels as text. No plotting, no graph view, no sequence highlighting. The tutorial notebook imports matplotlib/seaborn, but those are the **author's ad-hoc plots**, not tool functionality; architecture figures in the docs are hand-drawn SVGs.
*Source:* source-tree grep over 40 `.py` files; `oligopool/core_lenstat.py`; README/docs/api/agent-skills grep.

---

## Block B — assay coverage

*(Caption rule: `yes` = demonstrated end-to-end on published libraries with shipped tooling; `partial` = documented recipe only, variant content user-supplied.)*

### `assay_dms` = **partial**  *(carried from v1)*
*Evidence:* Infrastructure yes, content no — **no protein-level or reading-frame model in any module**. What ships is a documented DMS **recipe that takes the variants as input**: docs.md §"Saturation Mutagenesis Compression" (L2333 ff.), verbatim: *"**The ask**: 'I have 1000 single-amino-acid substitution variants. Can I synthesize fewer oligos?' … **The decomposition**: Generate all substitution variants, then `compress` into IUPAC-degenerate oligos."* — with a worked `op.compress` → `op.expand` losslessness check and *"# 2. Order synthesis_df for synthesis (6-20x fewer oligos than individual variants)"* (docs.md L2383; agent-skills.md L371 repeats the figure). agent-skills.md lists both *"Saturation mutagenesis (barcoded)"* and *"Saturation mutagenesis (degenerate)"* recipes.
*The limit:* no translation table, codon-usage model, reading frame, or amino-acid concept in any module. It builds the pool infrastructure around a DMS library **enumerated elsewhere** and compresses it 6–20×; it does not enumerate the substitutions. *Pre-empt the rebuttal:* `mutant_generator.py::generate_codon_variants` enumerates all 64 nucleotide triplets at a position — but it is example code with no protein model (it never translates).
*Source:* docs.md L2325–L2400 (esp. L2333, L2383); agent-skills.md L371; `mutant_generator.py`.

### `assay_mpra` = **yes**  *(carried from v1 — strongest cell in the record)*
*Evidence:* Paper title/keywords include *"massively parallel reporter assay"*; **three real libraries designed, built, and characterized**: 93,180 promoters (Fig. 2a/b), 62,120 5′UTR mRNA-stability elements (Fig. 2c/d), 6,232 hammerhead ribozymes (Fig. 2e/f). docs.md ships a **"Promoter MPRA Library"** application template verbatim (L2403–L2420): `background(host.fasta) → primer(fwd) → motif(constant, anchor) → motif(constant, anchor) → barcode → primer(rev, paired) → spacer(auto) → lenstat → verify → final`, followed by *"This same single-barcode MPRA pattern also applies to other regulatory/stability-element libraries (e.g., UTR variants, degrons, structured elements)."* Analysis Mode (`index`/`pack`/`acount`/`xcount`) then counts the MPRA reads, closing the loop.
***Caveat that must be stated:*** it designs the **synthesis-ready MPRA oligo** (primers, indexing anchors, barcodes, spacers, splits/pads) and analyzes the resulting sequencing data — **the regulatory variants themselves are user input.**
*Source:* paper title/abstract, ¶192/¶216/¶230, Fig. 2; docs.md L2403–L2420; api.md Analysis Mode.

### `assay_insilico` = **no**  *(carried from v1)*
*Evidence:* No sequence-to-function model, no in-silico perturbation, no genomic-AI integration, and no model-loading parameter on any of the 22 module signatures. The paper mentions ML only **downstream of the wet lab**: count matrices *"can be combined with biophysics and/or machine learning to develop predictive sequence-to-function models"*. **The README's "AI" section is about LLM agent tooling, not genomic AI**: *"Oligopool Calculator is optimized for AI-assisted workflows. Either share the `docs/agent-skills.md` file with your agent…"*. Re-verified this pass: a grep for `machine learning|\bML\b|predictive model|neural|deep learning|pytorch|tensorflow` across all four documentation files returns **exactly one hit** — docs.md L895: *"ML-generated libraries and saturation mutagenesis libraries often compress well"* — which places the model upstream and external.
*Source:* README.md AI section; docs.md L895; paper Discussion; api.md (all 22 signatures).

---

## Block C — genomics integration (all `no`)

Verified by independent positive-absence greps over the four documentation files: **zero hits** for `GTF|GFF|HGVS|VCF|chrom|coordinate|genome coord|exon|intron|reading frame`. The only `codon`/`amino` hits in the docs are two incidental ones (api.md L1686 IUPAC-table row *"M | A, C | Amino"*; docs.md L2333 *"single-amino-acid"*). The entire I/O contract is a CSV/DataFrame of `ID` + DNA-string columns — docs.md L1195: *"It assumes all non-ID columns are DNA strings."*

- **`genome_coordinates` = no.** No coordinate is ever parsed, tracked, or emitted. The only genome-shaped input is `background()`, whose `input_data` is *"list of DNA strings, CSV/DataFrame with a `Sequence` column, or FASTA path"* and whose only output is a k-mer DB directory. Paper: background is *"e.g., the genome sequence of E. coli"* — purely off-target k-mer screening. *Source:* api.md `background`; documentation-wide grep.
- **`transcript_models` = no.** Zero `GTF|GFF|transcript-model` hits across all four documentation files; no annotation-file parameter on any of the 22 signatures; inputs are CSV/DataFrame/FASTA sequence only. *Source:* api.md; documentation-wide grep.
- **`exon_intron_split_codons` = no.** Zero `exon|intron|reading frame` hits **across the four documentation files** (grep scoped explicitly, because the repo does ship a 64-entry `CODONS` table in `examples/library-compressor/mutant_generator.py`). No module has a reading-frame, exon, intron, or codon-boundary parameter. docs.md L1195: *"It assumes all non-ID columns are DNA strings."* *Source:* documentation-wide grep; api.md; docs.md L1195.
- **`vcf_vep_output` = no.** Zero VCF hits. Outputs are CSV design tables (`final` → `ID`/`CompleteOligo`/`OligoLength`), CSV count matrices from `acount`/`xcount`, and binary artifacts (`.oligopool.index`/`.pack`/`.background`) readable only via `inspect`. No VCF or VEP writer exists. *Source:* api.md `final`, `acount`, `xcount`, `inspect`.
- **`consequence_annotation` = no.** No molecular-consequence vocabulary anywhere. The only per-row annotations are **physical-design conflicts** from `verify`: `HasIntegrityConflict`, `HasLengthConflict`, `HasExmotifConflict`, `HasBackgroundConflict`, `HasAnyConflicts` (api.md L1369–L1377). Structurally impossible without a protein model, which does not exist. *Source:* api.md L1369–L1377; documentation-wide grep.

*(v1's `hgvs_input = no` is dropped per ROWS_v2. For the record, in case it is needed as prose: zero HGVS hits anywhere; `mutant_generator.py` **mints** HGVS-*style* IDs (`A5T`) but no HGVS **parser** exists in the package or the examples.)*

---

## Block D — adjacent / complementary

### `primer_design` = **yes**  *(carried from v1, with the scope caveat)*
*Evidence:* `op.primer(input_data, oligo_length_limit, primer_sequence_constraint, primer_type, minimum_melting_temperature, maximum_melting_temperature, maximum_repeat_length, primer_column, output_file=None, paired_primer_column=None, oligo_set=None, left/right_context_column=None, patch_mode=False, excluded_motifs=None, background_directory=None, random_seed=None, verbose=True)`. Verified behaviours (docs.md L468–L474, verbatim): *"If `paired_primer_column` is provided, the paired primer type is inferred and Tm matching is applied within 1 °C."* / *"When `oligo_set` is provided, primers are designed per set and screened for cross-set compatibility."* / *"`primer_sequence_constraint` is flexible by design; mix fixed and degenerate bases to enforce patterns (for example, GC clamp-like starts such as `'SS' + 'N'*18`)."* / *"Structural screening is first-class: candidates with strong hairpin/homodimer/heterodimer behavior … are rejected during optimization."* GC content is computed and reported (`primer_guanine_cytosine_content`, 6× in `primer.py`). Paper: 7 explicit primer design rules; benchmark *"universal primer binding sites for one million 200-mer oligos in 15 min."* `split` additionally designs Tm- and Hamming-constrained overlap regions; `pad` designs pad primers with Type IIS sites from a **34-enzyme built-in catalogue** with 3′ cut-offset models.
***Caveat — essential, keep it:*** these are **amplification / assembly** primers, **not mutagenic** primers. There is no site-directed-mutagenesis primer mode (contrast Mutation Maker). Under the v2 row text's first clause ("**Mutagenic** primer") the answer is `no`; under its second ("oligo design for wet-lab protocols") it is clearly `yes`. `yes` is recorded because the row text is disjunctive and the capability is genuinely, heavily present.
*Source:* api.md `primer`, `split`, `pad`, §"Type IIS Enzymes" (L1724–L1730); docs.md L468–L474; `oligopool/primer.py`; paper ¶"Primers".

### `codon_optimization` = **no**  *(carried from v1)*
*Evidence:* Zero `codon` hits **in the four documentation files** (grep scoped explicitly; the repo ships a 64-entry `CODONS` table in `mutant_generator.py`, but that is codon **enumeration** for a compression demo, not codon optimization). No translation table, no codon-usage table, no CDS handling, no reading-frame parameter on any of the 22 signatures. Nothing in the package computes or optimizes synonymous recoding.
*Source:* documentation-wide grep; api.md (all 22 signatures); `mutant_generator.py`.

### `synthesis_constraints` = **yes**  *(carried from v1 — one of its strongest cells)*
*Evidence:* Enforced per sequence, at design time and again at QC time.
*Length:* `oligo_length_limit` is a **required** parameter on every design module; `lenstat` reports length statistics and free space; over-length constructs are resolved **constructively** by `split` (Tm- and Hamming-constrained overlaps for Gibson/homology assembly) → `pad` (Type IIS pads for scarless Golden Gate), not merely flagged.
*Forbidden sequences:* homopolymers and restriction sites are banned via `excluded_motifs`, which accepts a list, CSV/FASTA path, DataFrame with an `Exmotif` column, comma-string, or multiple named sources. Documented usage, verbatim: docs.md L360 `op.barcode(..., excluded_motifs={'cutsites': 'cutsites.csv', 'homo': ['AAAA', 'TTTT']})`; docs.md L2223 `'homopolymers': ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']`. Motifs must be strict ATGC (no IUPAC).
*Repetitiveness (a different mechanism):* `maximum_repeat_length` caps the maximum **shared** repeat length between the designed element and the input oligos / flanking context, k = `maximum_repeat_length + 1` (api.md L397, L409). docs.md L468: *"`maximum_repeat_length` controls non-repetitiveness against `input_data` only; screening against a background requires `background_directory`."*
*Off-target:* `background_directory` k-mer screening against host genome/plasmid/phiX DBs, junction-aware when context columns are given, layerable across multiple DBs, re-checked in `verify`.
*Thermodynamics:* Tm windows and hairpin/homodimer/heterodimer/cross-dimer rejection for primers and overlaps.
*Dedicated pre-synthesis QC:* `verify` emits per-row conflict flags plus four JSON detail columns with per-position, per-column attribution. **Motif emergence at junctions is checked** — `ExmotifConflictDetails.library_baseline` and `excess_occurrences = occurrences - library_baseline` flag motifs that *emerge from concatenation* rather than counting absolute occurrences. Most competing tools do not do this. Cost is addressed separately via `compress`.
*Source:* api.md L100, L184, L264, L338, L397, L409, L1369–L1405, `split`/`pad`/`lenstat`; docs.md L360, L468, L2223; paper ¶"Barcodes", ¶"Primers".

### `degenerate_iupac_codons` = **yes**  *(NEW row)*
*Evidence:* Degenerate/IUPAC handling is a **whole shipped mode** (Degenerate Mode), in both directions, plus IUPAC constraint syntax across the design modules. Re-verified at source this pass.
(i) **Compression:** `op.compress(input_data, mapping_file=None, synthesis_file=None, rollout_simulations=100, rollout_horizon=4, random_seed=None, verbose=True)` collapses a concrete A/T/G/C variant library into IUPAC-degenerate synthesis oligos by Monte-Carlo rollout, returning `mapping_df` (`ID`, `Sequence`, `DegenerateID`) and `synthesis_df` (`DegenerateID`, `DegenerateSeq`, `Degeneracy`, `OligoLength`) plus `compression_ratio`/`min|max|mean_degeneracy` stats. api.md `compress` Notes, verbatim: *"Core guarantee (lossless): `degeneracy(prefix) <= count(compatible variants)` - no invented sequences"*; docs.md/agent-skills.md document **6–20× fewer oligos**.
(ii) **Expansion:** `op.expand(input_data, sequence_column, mapping_file=None, output_file=None, expansion_limit=None)` — api.md L838, verbatim: *"**Purpose**: Expand IUPAC-degenerate sequences into all concrete A/T/G/C sequences"*, with the full code table documented in-line: *"IUPAC codes: N=any, R=A/G, Y=C/T, S=C/G, W=A/T, K=G/T, M=A/C, B=C/G/T, D=A/G/T, H=A/C/T, V=A/C/G"*, an `expansion_limit` safety cap, and CPU-parallel expansion.
(iii) **Design-time IUPAC constraints:** `motif_sequence_constraint` (paper: motifs are *"defined using the IUPAC degenerate nucleotide code"*) and `primer_sequence_constraint` (docs.md L471: *"mix fixed and degenerate bases … `'SS' + 'N'*18`"*); api.md L1686 ships a full IUPAC code table. `IUPAC|degenerate` occurs 29× in docs.md, 28× in api.md, 8× in README.md, 5× in agent-skills.md.
*The one qualifier:* this is **nucleotide-level** IUPAC degeneracy. There is no codon/reading-frame model, so "degenerate **codon** design" in the classical NNK/NNS protein-library sense is only reachable by the user supplying triplets (as `mutant_generator.py::generate_codon_variants` does by enumerating all 64 nucleotide triplets, never translating). The row's "expansion, or compression" clauses are satisfied outright.
*Source:* api.md L774–L832 (`compress`), L834–L878 (`expand`), L230–L262 (`motif`), L1686; docs.md L471, L895, L2333–L2400; paper §"Restriction Sites, Motifs, and Spacers".

### `negative_control_generation` = **no**  *(NEW row)*
*Evidence:* Verified by positive-absence this pass: a grep for `scrambl|shuffl|negative control|control sequence|randomi[sz]|dinucleotide` across all four documentation files (README.md, docs.md, api.md, agent-skills.md) returns **zero hits**. No module generates scrambles, dinucleotide-preserving shuffles, or reverse/complement negative controls, and no module has a control-generation parameter.
*Pre-empt the two obvious rebuttals:* (1) `revcomp` exists, but api.md L480 is explicit that its **purpose** is *"Reverse complement a range of columns and reverse their order"* — a whole-table architecture operation (parameterized only by a column range) used to flip oligo orientation, not a per-variant control generator. (2) `spacer` generates filler sequence to pad oligos to a target length under constraints — that is length management, not a designed negative control. Control sequences, like all variant content, are the user's input.
*Source:* documentation-wide grep (zero hits, re-run 2026-08-10); api.md L476–L504 (`revcomp`), L306–L374 (`spacer`); docs.md L1304.

### `ml_model_in_loop` = **no**  *(NEW row)*
*Evidence:* No design step is driven by a predictive model's output. No model, checkpoint, scorer, or predictor parameter appears on any of the 22 module signatures; there is no torch/tensorflow/sklearn dependency in the shipped package. Re-verified this pass: `machine learning|\bML\b|predictive model|neural|deep learning|pytorch|tensorflow` across all four documentation files returns **exactly one hit** — docs.md L895, verbatim: *"ML-generated libraries and saturation mutagenesis libraries often compress well"* — i.e. the model runs **upstream and outside** the tool, and Oligopool merely compresses what it produced. The paper's only ML mention is **downstream**: count matrices *"can be combined with biophysics and/or machine learning to develop predictive sequence-to-function models"*. The README's "AI" section refers to **LLM agent tooling** (`docs/agent-skills.md`) for driving the CLI, not to a model in the design loop — this distinction must be kept, because the word "AI" in the README is otherwise easy to misread.
*Source:* docs.md L895; README.md AI section; paper Discussion; api.md (all 22 signatures); documentation-wide grep.

### `readout_analysis` = **yes**  *(NEW row — a headline strength; complementarity, not competition)*
*Evidence:* Half the tool is a sequencing-readout analyzer, and **no other tool in this survey closes the design→measurement loop.** Analysis Mode ships four modules: `index` (builds a barcode index with an explicit anchor model — `barcode_prefix_column`/`barcode_suffix_column` with `barcode_prefix_gap`/`barcode_suffix_gap` = *"Bases between prefix anchor and barcode in read"*, plus a parallel associate set; api.md L906–L956), `pack` (efficient chunked read packing), `acount` (association counting), `xcount` (combinatorial multi-barcode counting). Supporting machinery: the **Scry** barcode classifier; built-in **phiX spike-in filtering** (`oligopool/base/phiX.py` k-mer spectrum, `phix_match` failure category *"PhiX contamination detected; first k-mer (k=30) match"*, `phix_reads` counter); **index-hopping detection**; **failed-read diagnostics** — docs.md L2613, verbatim: *"Use `failed_reads_file` to sample discarded reads by failure category (anchor missing, barcode absent, barcode ambiguous, callback rejected, etc.)"*; a user-extensible **`callback(r1, r2, ID, count, coreid)` hook** in `acount`/`xcount` (api.md L1154–L1162); throughput ~500 M reads/hour on 8 cores. The paper's three published libraries were counted with this half of the tool. Analysis workflows are the documented sweet spot for the DAG runner (docs.md L1634).
*Source:* api.md L889–L1288 (`index`/`pack`/`acount`/`xcount`), L906–L956, L1113, L1154–L1162, L1235, L1666; docs.md L1634, L2483, L2613; `oligopool/base/phiX.py`, `core_count.py`; paper Analysis Mode / Fig. 2.

---

## Block E — engineering and availability

*(Per instruction, `value = "yes"` for all five; the real answer is in the evidence.)*

### `interface`
**Python API + CLI (two console scripts) + YAML pipelines (serial / DAG / mixed) + Jupyter + Docker. No web server, no GUI.**
README.md L20, verbatim: *"`Oligopool Calculator` is a Swiss-army knife for oligo pool libraries: a unified toolkit for high-throughput design, assembly, compression, and analysis of massively parallel assays, designed to integrate seamlessly with Python, the CLI, Jupyter, containers, and AI-assisted workflows."* Installs **two console entry points**, `oligopool` and `op` (README L166). CLI supports YAML config for single commands and multi-step pipelines, `--dry-run`, `--stats-json` / `--stats-file path.json`, `--quiet` (docs.md L1477–L1479), YAML merge-key parameter sharing, and documented CLI > config > defaults precedence.
***Warning to preserve:*** **`oligopool.com` is NOT this tool** — an unrelated generic browser-side oligo calculator citing SantaLucia 1998 / Breslauer 1986, with no mention of Salis, Hossain, or the ACS paper. A web-search summariser misidentified it. **No hosted web service exists and none is claimed.**
*Source:* README.md L20, L166; docs.md L1477–L1479, L1630–L1700; live fetch of oligopool.com.

### `license`
**GPL-3.0-only.** Confirmed three independent ways: PyPI JSON `license_expression: "GPL-3.0-only"`; GitHub repo metadata `license.spdx_id: "GPL-3.0"`; paper Data Availability (p.4220): *"available at https://github.com/hsalis/SalisLabCode and https://github.com/ayaanhossain/oligopool under a GPLv3 open-source license."* `LICENSE` file (35,155 B) in repo root.
*Source:* PyPI JSON API; GitHub REST API; paper p.4220.

### `installable_today`
**Yes — pip, git, or Docker.** PyPI `oligopool` is alive: latest **2026.2.22.1** (wheel + sdist), `requires_python: <4,>=3.10`. Re-verified at source this pass, README L69/L75/L80: `pip install --upgrade oligopool`; `pip install git+https://github.com/ayaanhossain/oligopool.git`; *"If you are on `Windows` or simply prefer to, `Oligopool Calculator` can also be used via `Docker`"* (`Dockerfile` in repo root + `docs/docker-notes.md`). GitHub repo public and unarchived. Platform: **Linux, macOS, WSL — not native Windows**; Docker is the documented Windows path.
***Not verified by execution*** — the survey ran under a no-install constraint. Nothing was installed or run; installability on a current toolchain is inferred from PyPI metadata and the README, not tested.
*Source:* PyPI JSON API (re-fetched 2026-08-10); README.md L69, L75, L80, L106; repo root `Dockerfile`.

### `last_activity`
**2026-02-22 — about 6 months before this survey (2026-08-10), with no activity since.** Last repo push **2026-02-22T02:46:27Z**; latest PyPI release **2026.2.22.1** uploaded **2026-02-22T02:48:39Z**; `open_issues_count = 0`.
***Do NOT characterize this as "rapid iteration."*** There have only ever been **two PyPI releases in total** — `2024.12.2` (2024-12-02) and `2026.2.22.1` (2026-02-22) — **14 months apart**, followed by ~6 months of silence. What happened within four days in Feb 2026 were two *in-repo version bumps* (v2026.02.18 → v2026.02.22), only one of which was published. State the dates and stop. The tool is current and complete, not abandoned; it is simply quiet.
*Source:* PyPI JSON API (full release list); GitHub REST API `pushed_at` / commit list / `open_issues_count`.

### `documented_examples`
**Rich — 4 runnable example projects + 1 tutorial notebook + 5 in-docs application templates + 3 published libraries in the paper.** Re-verified at source this pass, README L98: the `examples` directory includes *"a design parser, a library compressor, an analysis pipeline, and a complete CLI YAML pipeline"*; README L100 points to `examples/OligopoolCalculatorInAction.ipynb` for *"the full end-to-end walkthrough"* (284 cells, re-executed 2026-02-18, building a 6,232-variant ribozyme oligopool through Design → Assembly → Degenerate → Analysis Mode). docs.md §"Application Templates" (L2403 ff.) adds five worked templates: Promoter MPRA Library, CRISPR Guide Library, Ribozyme Library (pre/post-cleavage dual barcode), Saturation Mutagenesis Compression, Selection-Based Discovery Without Barcodes. Documentation totals ~214 kB across `docs.md` / `api.md` / `agent-skills.md`, plus `docker-notes.md` and `examples/README.md`, all consistent with v2026.2.22.1.
*Highest-value reproduction target for PoolParty:* `examples/design-assembly-parser/` — a declarative spec-based parser over 4,351 real promoter sequences, architecture `[Primer1]-[Cut1]-[Promoter]-[Barcode]-[Primer2]-[Cut2]-[Primer3]-[Filler]`, `element_type ∈ {variant, primer, barcode, motif, spacer}`, with `get_primer_order()` topologically resolving the paired-primer dependency graph. Closest analogue in the surveyed ecosystem to PoolParty's declarative library spec — **but a worked example, not a shipped API.**
*Source:* README.md L98–L106; `examples/README.md`; docs.md L2403 ff.

---

## Summary table

| Key | v2 value | v1 value | Changed? |
|---|---|---|---|
| `library_first_class_object` | **yes** | `library_as_object` = partial | **↑ yes** (split resolves the conflation) |
| `composable_operations` | yes | `dag_chaining` = yes | rename, carried |
| `lazy_generation` | no | `lazy_evaluation` = no | rename, carried |
| `library_algebra` | **no** | `library_as_object` = partial | **↓ no** (no row-union operator) |
| `exhaustive_single_scans` | **partial** | `mixed_mutagenesis_one_pool` = no | **↑ partial** (shipped example enumerates all 3L singles) |
| `sampled_random_mutagenesis` | no | `mixed_mutagenesis_one_pool` = no | carried |
| `higher_order_combinatorial` | **partial** | `mixed_mutagenesis_one_pool` = no | **↑ partial** (shipped example enumerates 4^N multi-position) |
| `heterogeneous_components_one_library` | no | `mixed_mutagenesis_one_pool` = no | carried (single-strategy dispatch + no row-union) |
| `combinatorial_motif_place` | no | no | carried |
| `barcode_generation` | yes | yes | carried |
| `per_sequence_provenance` | partial | partial | carried |
| `automatic_naming` | no | no | carried |
| `design_visualization` | no | no | carried |
| `assay_dms` | partial | partial | carried |
| `assay_mpra` | yes | yes | carried |
| `assay_insilico` | no | no | carried |
| `genome_coordinates` | no | no | carried |
| `transcript_models` | no | no | carried |
| `exon_intron_split_codons` | no | no | carried |
| `vcf_vep_output` | no | no | carried |
| `consequence_annotation` | no | no | carried |
| `primer_design` | yes | yes | carried (amplification/assembly, not mutagenic) |
| `codon_optimization` | no | no | carried |
| `synthesis_constraints` | yes | yes | carried |
| `degenerate_iupac_codons` | **yes** | — | NEW |
| `negative_control_generation` | **no** | — | NEW |
| `ml_model_in_loop` | **no** | — | NEW |
| `readout_analysis` | **yes** | — | NEW |
| `interface` | yes (Python+CLI+YAML+Jupyter+Docker; no web/GUI) | yes | carried |
| `license` | yes (GPL-3.0-only) | yes | carried |
| `installable_today` | yes (pip/git/Docker; not executed) | `maintained` = yes | NEW (derived) |
| `last_activity` | yes (2026-02-22; ~6 mo quiet) | `maintained` = yes | NEW (derived) |
| `documented_examples` | yes (4 projects + notebook + 5 templates) | — | NEW |

**`hgvs_input` is dropped per ROWS_v2 and is not reported.**

---

## Cells flagged for human review

1. **`library_first_class_object` = yes** — upgraded from v1 `partial`. Defensible as `partial` if the row is read to require a *dedicated* library class with design semantics; it is a generic pandas DataFrame, and `final` discards the annotation columns.
2. **`exhaustive_single_scans` = partial** — the capability lives in `examples/library-compressor/mutant_generator.py`, not in `op.*`. `no` is defensible on a strict shipped-API reading. **Must be scored consistently with how example/vignette code was treated for the other 11 tools.** Also: substitutions only, no deletion/insertion scan.
3. **`higher_order_combinatorial` = partial** — identical caveat to (2); same file, same consistency question.
4. **`heterogeneous_components_one_library` = no** — hinges on reading the row as row-wise *variant-class* heterogeneity (per the ROWS_v2 gloss). Under an element-type reading it would be `yes` (`design-assembly-parser` mixes variant/primer/barcode/motif/spacer in one spec).
5. **`library_algebra` = no** — airtight for combine/sample/repeat, but `compress`/`expand` do change row counts, and a referee may point at them.
6. **`composable_operations` = yes** — the DAG is a scheduling graph over one flat table, not content composition; no nesting; `compress` is a terminal branch and `pad` is per-fragment. `partial` is arguable.
7. **`per_sequence_provenance` = partial** — judgement call, deliberately generous. Per-column attribution and `mapping_df` are real; derivation history does not exist.
8. **`assay_dms` = partial vs `assay_mpra` = yes** — the asymmetry is deliberate and *only* holds under the stated caption rule. If a single uniform rule is preferred, both become `partial (infrastructure only)`. Do not leave the asymmetry unexplained — it is the first thing a referee will poke.
9. **`primer_design` = yes** — only under the row's second clause ("oligo design for wet-lab protocols"). Under "**mutagenic** primer" it is `no`. The row text is disjunctive, hence `yes`; if the published row text is narrowed, this cell flips.
10. **`degenerate_iupac_codons` = yes** — unambiguous for IUPAC degeneracy, design/expansion/compression. But there is no codon or reading-frame model, so "degenerate **codon**" in the NNK/NNS protein-library sense is not natively expressible.
11. **`automatic_naming` = no** — true for variant IDs (the intended scope), but `compress` auto-mints `DegenerateID`s and five modules auto-create column names. If the row is read as "any automatic naming", this flips.
12. **`negative_control_generation` = no** — zero-hit grep is strong, but `revcomp` exists and could be waved at; the evidence pre-empts that.
13. **`barcode_generation` = yes** — with the verified gap that **GC is not constrainable for barcodes**. Footnote it if the published row text names GC.
14. **`installable_today` = yes** — inferred from PyPI/README metadata only. **Nothing was installed or run** (no-install constraint), so this is unverified by execution.
15. **`last_activity` = 2026-02-22** — the GitHub REST API rate-limited the third verification pass; the figure rests on two prior independent fetches that agree exactly, corroborated by the PyPI upload timestamp two minutes after the final push.
16. **`documented_examples`** — the count depends on what is counted (4 projects? +1 notebook? +5 docs templates? +3 paper libraries?). Normalize the counting rule across all 12 tools before publishing.
