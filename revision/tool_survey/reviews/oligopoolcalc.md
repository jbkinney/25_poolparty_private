# Adversarial review — Oligopool Calculator (`oligopoolcalc`)

**Reviewer stance:** try to falsify the extraction. Referees may include Hossain / Salis.
**Review date:** 2026-08-10. **No installs; document + web research only.**

## What I did independently

Downloaded fresh copies of `README.md`, `docs/docs.md` (98,900 B), `docs/api.md` (76,438 B),
`docs/agent-skills.md` from `raw.githubusercontent.com/ayaanhossain/oligopool/master`; pulled the
full repo file tree via the GitHub trees API; pulled `pypi.org/pypi/oligopool/json`, the GitHub repo
metadata and commit list; extracted the paper PDF with PyMuPDF (15 pp., 92,921 chars); read
`examples/design-assembly-parser/{README.md,run_design_assembly_parser.py,design_assembly_parser.py}`,
`examples/library-compressor/mutant_generator.py`, `examples/cli-yaml-pipeline/mpra_design_parallel.yaml`;
grepped all 40 `oligopool/**/*.py` source files for plotting imports; independently checked
`oligopool.com`.

**Headline:** the extraction's *values* are almost all correct and the framing
("infrastructure yes, variant content no") survives scrutiny. But **four evidence strings contain
statements that a hostile author-referee can falsify**, including **one fabricated quotation**. Those
must be fixed before this reaches a referee response. And one repo file the extractor never opened —
`examples/library-compressor/mutant_generator.py` — is the single biggest attack surface in the whole
extraction.

---

## The one thing the extractor missed that matters most

`examples/library-compressor/mutant_generator.py` (8,431 B, header `Author: Ayaan Hossain`) ships in
the repo and implements:

- `generate_single_mutants()` — all 3L single-nucleotide mutants, **auto-minting HGVS-style variant IDs**
  (`f'{ref_base}{pos+1}{alt_base}'` → `A5T`)
- `generate_codon_variants()` — **all 64 codons at a chosen codon position** (a `CODONS` table built from
  `itertools.product(DNA_BASES, repeat=3)`), IDs like `ATG_AAA`
- `generate_multi_position_variants()` — **all 4^N combinatorial variants at an arbitrary position list**,
  IDs like `p1A_p5C`
- all three are **Python generators (`yield`)**, and `generate_variant_dataframe(strategy=...)` wraps them
  into a `DataFrame` with `ID`/`Sequence`

It is an *example helper*, not part of the `oligopool` package API (`op.*` has no such function), and
`generate_variant_dataframe` takes **one** strategy per call, so the capability values `no` still stand.
But four evidence strings are worded as blanket absences that this file directly contradicts:

| Row | Wording that is falsifiable | Fix |
|---|---|---|
| `mixed_mutagenesis_one_pool` | "The tool does not generate mutagenized variants **at all**" | "The `oligopool` package has no mutagenesis module; the repo ships a *worked-example helper* (`examples/library-compressor/mutant_generator.py`) that enumerates single, codon-saturation and multi-position variants, but it is example code, exposes one strategy per call, and cannot express a mixed-mutagenesis spec." |
| `automatic_naming` | "no module invents or enumerates sequence names"; "Output COLUMN names are user-specified" | The second clause is simply wrong (see below). The example helper also auto-mints variant IDs. |
| `exon_intron_split_codons` / `codon_optimization` | "Grep for 'codon' … returns zero hits" | True **for the four documentation files only** — say so explicitly. The repo has a `CODONS` table and `generate_codon_variants`. Codon *optimization* (synonymous recoding / usage tables) genuinely does not exist anywhere. |
| `assay_dms` | "zero codon/amino-acid awareness anywhere in the documentation" | `docs.md` L2333 literally says "single-amino-acid substitution variants". Rephrase to "no protein-level or reading-frame model in any module". |

---

## Verified-correct core (no change needed)

- **Zero hits** for `GTF|GFF|HGVS|VCF|chrom|coordinate|exon|intron|reading frame` across README + docs.md
  + api.md + agent-skills.md. `codon`/`amino` hits: exactly two, both incidental
  (`api.md` L1686 IUPAC table row "M | A, C | Amino"; `docs.md` L2333 "single-amino-acid"). Blocks C and
  `codon_optimization` are safe.
- **No GC-content parameter for barcodes.** I read the full `op.barcode` signature: `oligo_length_limit,
  barcode_length, minimum_hamming_distance, maximum_repeat_length, barcode_column, barcode_type,
  left/right_context_column, patch_mode, cross_barcode_columns, minimum_cross_distance, excluded_motifs,
  background_directory, random_seed, verbose`. No sequence-constraint and no GC argument. The paper's four
  barcode criteria (PDF ¶"Barcodes") are: ≥1 barcode per variant; minimum pairwise Hamming distance;
  minimal length; no homopolymer tracts + researcher-defined excluded sequences. **GC is never mentioned.**
  The extractor's "one honest gap" is correct and should be kept. *Refinement:* GC **is** computed and
  reported for primers — `primer_guanine_cytosine_content` appears 6× in `oligopool/primer.py` stats — so
  say "not constrainable for barcodes", not "absent".
- **No visualization.** I grepped **all 40 `oligopool/**/*.py` files** for `seaborn|matplotlib|pyplot`:
  **zero hits.** (`seaborn>=0.13.2` is a declared PyPI dependency but no shipped module imports it — a
  stale requirement.) `lenstat` prints a per-element textual length report
  (`core_lenstat.py`: `' Element {}: Occupies {:,} Base Pair(s)'`), which is textual, not graphical.
  `design_visualization = no` is airtight.
- **`motif` is not combinatorial.** One motif → one named column per call; `motif_type` ∈
  {`0`/`variable` = unique per row, `1`/`constant` = one anchor for all rows}. No position, orientation or
  motif-set enumeration. `revcomp` flips a contiguous **column range** for the whole table. Verified.
- **`join` does not compose sub-libraries.** api.md `join` notes: "`ID` sets must match exactly across both
  inputs"; "`join` never creates or drops rows; mismatched IDs are an error". Combined with the module
  inventory (no row-union operator anywhere), this is the **sharpest and safest** contrast with PoolParty
  and the extraction only implies it. Promote it to an explicit sentence.
- **License / maintenance facts.** PyPI `license_expression: GPL-3.0-only`, `requires_python <4,>=3.10`,
  version `2026.2.22.1` uploaded `2026-02-22T02:48:39`. GitHub `license.spdx_id: GPL-3.0`,
  `pushed_at 2026-02-22T02:46:27Z`, `open_issues_count 0`, 6 stars, created 2020-08-06. Paper: "under a
  GPLv3 open-source license" (PDF p.4220). All confirmed.
- **`oligopool.com` is genuinely unrelated.** Fetched it: "Free, professional tools for oligonucleotide
  design and analysis… peer-reviewed methods, including SantaLucia 1998 and Breslauer 1986", 19 generic
  browser-side calculators, **no mention of Salis, Hossain, or the ACS paper**. A web-search summariser
  misidentified it as a hosted version of the tool — the extractor's caution here was correct and worth
  keeping in the memo so nobody re-makes that mistake.
- **Paper benchmarks** all reproduce verbatim: "over four million highly unique and compact barcodes in
  1.2 h", "universal primer binding sites for one million 200-mer oligos in 15 min", "about 500 million
  deep sequencing reads per hour", "one million sequence variants takes about 35 min".
- **`examples/design-assembly-parser`** verified: `promoters.txt` = **4,351 lines** exactly (README says
  "~4,500"); `element_type` ∈ {variant, primer, barcode, motif, spacer}; architecture
  `[Primer1]-[Cut1]-[Promoter]-[Barcode]-[Primer2]-[Cut2]-[Primer3]-[Filler]`; `Cut1` constraint
  `'NNNGGATCCNNN'`; `Barcode` length 16, Hmin 3; `oligo_length_limit` 250. Also — the extractor missed it —
  the parser accepts `split_spec` and `padding_spec` in the same declarative call, and
  `get_primer_order()` *topologically resolves the paired-primer dependency graph*. The "closest analogue
  to PoolParty's declarative spec" characterisation is if anything conservative.

---

## Problems found, row by row

### 1. `synthesis_constraints = yes` — **contains a fabricated quotation**

The evidence says:

> `maximum_repeat_length` bans homopolymers and shared repeats … (docs.md: set it **'to 4 to reject AAAA'**)

**That string does not exist.** `grep -niE "reject|to 4 "` over docs.md returns only four hits, none of
them about homopolymers (they concern primer structural screening, alias collisions, a callback code
comment, and failed-read categories). The same grep over api.md / README.md / agent-skills.md is empty.

The mechanism claim is also wrong. Homopolymers are banned via **`excluded_motifs`**, not
`maximum_repeat_length`:
- docs.md L360: `op.barcode(..., excluded_motifs={'cutsites': 'cutsites.csv', 'homo': ['AAAA','TTTT']})`
- docs.md L2223: `'homopolymers': ['AAAAA','CCCCC','GGGGG','TTTTT']`
- docs.md L468: "`maximum_repeat_length` controls non-repetitiveness against `input_data` only"
- paper: barcodes must "not contain homopolymer tracts … while also excluding a set of researcher-defined
  sequences"

`maximum_repeat_length` limits *shared repeat length against the input oligos / context / background DB*
(k-mer size = `maximum_repeat_length + 1`), which caps homopolymer runs only incidentally and at ~k+1, not
at 4. **The value `yes` is correct and richly supported by the rest of the evidence** (`oligo_length_limit`
on every module, `excluded_motifs` from list/CSV/FASTA/DataFrame/dict, `background_directory`, Tm +
hairpin/homodimer/heterodimer/cross-dimer screening, `verify`'s five `Has*Conflict` columns, `lenstat`,
`split`→`pad`). **Delete the fabricated quote and fix the mechanism sentence.**

### 2. `automatic_naming = no` — evidence contains a false blanket statement

"Output COLUMN names are user-specified (e.g. `barcode_column='BC1'`)" is wrong for five modules:
- `split` → auto-creates `Split1`, `Split2`, … (count varies per oligo; even-numbered splits are
  auto-reverse-complemented)
- `pad` → auto-creates `5primeSpacer`, `ForwardPrimer`, `<split_column>`, `ReversePrimer`, `3primeSpacer`
- `final` → `CompleteOligo`, `OligoLength`
- `expand` → `ExpandedSeq`, `OligoLength`
- `compress` → `DegenerateID`, `DegenerateSeq`, `Degeneracy` — and `synthesis_df` is a table of **new rows
  the user never named**, identified solely by auto-minted `DegenerateID`s

The value `no` still holds for *library-member* naming (the input `ID` column is required and
user-supplied, and no module renames or enumerates variant IDs), but the cell should say that precisely
rather than making a claim about column names that is easy to disprove.

### 3. `mixed_mutagenesis_one_pool = no` — see the `mutant_generator.py` section above

Value stands; the "at all" absolute must go.

### 4. `maintained = yes` / `availability_status` — factual error in the spin

`availability_status` says "two releases within four days in Feb 2026, i.e. rapid active iteration".
**There have only ever been two PyPI releases in total**: `2024.12.2` (2024-12-02) and `2026.2.22.1`
(2026-02-22) — 14 months apart. What happened within four days were two *in-repo version bumps*
(v2026.02.18 → v2026.02.22), only one of which was published. And the repo has had **zero commits in the
~6 months since 2026-02-22**. `maintained = yes` is still fair (last touched 6 months before the survey,
0 open issues, docs current), but "rapid active iteration" is a claim a referee could invert
("actually it's been quiet since February"). State the dates and stop.

### 5. Two paraphrases presented as verbatim quotes (low severity, fix anyway)

- `interface` evidence: README.md quoted as *"seamlessly integrating with Python, CLI, Jupyter, containers,
  and AI workflows"*. Actual text: *"designed to integrate seamlessly with Python, the CLI, Jupyter,
  containers, and AI-assisted workflows."*
- `library_as_object` / `per_sequence_provenance`: *"Keep a saved annotated design DataFrame/CSV before
  `final` if you plan to `index`, `pack`, and count later"*. Actual docs.md L762: *"Save your annotated
  design DataFrame/CSV before `final` if you plan to `index`, `pack`, and count later."* (There is a
  separate, similar note at L821 about `split`.)

Substantively harmless, but if a referee checks one quote and finds it isn't verbatim, they will check
them all.

### 6. `per_sequence_provenance = partial` — value right, evidence undersells

`verify` gives genuinely **per-row, per-element** attribution that the evidence doesn't mention:
- `LengthConflictDetails.column_lengths` — per-column bp contribution for that oligo
- `ExmotifConflictDetails` — a list of `{motif, occurrences, library_baseline, excess_occurrences,
  positions (0-based offsets), columns (which element contributed each occurrence)}`
- `compress` `mapping_df` (`ID`, `Sequence`, `DegenerateID`) plus `expand(mapping_file=…)` restoring
  original `ID`s is real per-sequence provenance across the compression step

`partial` is still the right call (no derivation/mutation history, and `stats` is aggregate — docs.md L176
"`stats` is **aggregate**: programmatic pass/fail decisions and summary reporting" — verbatim, confirmed).
But if the referee is Hossain, quote the `verify` detail schema rather than let him supply it.

### 7. `dag_chaining = yes` — value right, evidence cites the weaker example

The extractor cited only `analysis_multi.yaml`. There is a **runnable Design-Mode DAG**:
`examples/cli-yaml-pipeline/mpra_design_parallel.yaml` — `primer` → (`barcode_a` ∥ `barcode_b`) →
`join` → `spacer` → `final`, with a comment acknowledging it is "intentionally 'mostly sequential' with
one parallel fan-out". Cite this one; it pre-empts "you only found our analysis DAG".

The extractor's caveat is exactly right and airtight: the DAG is an **execution/scheduling** graph over
design *steps* on one flat table (verified via `join`'s ID-set constraint), not a **content-composition**
graph over sub-populations. Keep that sentence verbatim. The prior-analysis correction ("no DAG
architecture" is outdated) is also correct and important.

### 8. `assay_mpra = yes` vs `assay_dms = partial` — the asymmetry is real

Both are "infrastructure yes, content no". The extractor already flags this in `confidence_notes`. My view:
keep `assay_mpra = yes` (three real published MPRA libraries; a dedicated shipped template; Analysis Mode
closes the loop) and keep `assay_dms = partial`, but make the *rule* explicit in the caption — e.g.
"yes = demonstrated end-to-end on published libraries of this type; partial = documented recipe only" —
otherwise the inconsistency is the first thing a referee will poke.

---

## Capabilities the extraction missed entirely

1. **`examples/library-compressor/mutant_generator.py`** — see top section. Highest priority.
2. **Reproducibility as a first-class feature.** docs.md §"Reproducibility": "All stochastic design modules
   support `random_seed` … Same seed + same inputs = same outputs. Great for debugging and publications."
   `random_seed` is echoed back in the stats dict for `barcode`/`primer`/`motif`/`spacer`/`split`/`pad`/
   `compress`. PoolParty presumably claims reproducibility too — do not let this look like a differentiator.
3. **Machine-readable output contract.** CLI `--stats-json` (stats dict to stdout) and `--stats-file
   path.json`, plus `--quiet`; api.md ships a formal **Stats Dict Schema** (`status`, `basis`, `step`,
   `step_name`, `vars`, `warns`, `module`, `input_rows`, `output_rows`) with a documented `basis`
   vocabulary (`solved`/`complete`/`unsolved`/`infeasible`/`verified`/`conflicts`/`corrupted`).
4. **User-extensible counting via `callback`** — Python-only per-read accept/reject hook in
   `acount`/`xcount` ("Returns: bool - True to accept, False to reject"). An extension point.
5. **Failed-read diagnostics** — `failed_reads_file` with per-category attribution (anchor missing,
   barcode absent, `barcode_ambiguous`, callback rejected) plus `DiagnosticDetails` JSON columns.
6. **Built-in phiX spike-in filtering** — `oligopool/base/phiX.py` ships the phiX k-mer spectrum;
   `core_count.py` imports it and exposes `phix_match` / `phix_reads` counters.
7. **Built-in enzyme knowledge base** — 34 Type IIS enzymes with recognition motif *and* 3′ cut-offset
   models (`BsaI: GGTCTC + N*5`), chosen so pad excision is scarless; plus a common-restriction-site
   appendix. Not just a string list.
8. **`index` anchor model** — `barcode_prefix_column` / `barcode_suffix_column` with
   `barcode_prefix_gap` / `barcode_suffix_gap`, and a parallel set for *associates*. This is what lets
   barcodes be found at variable offsets in reads; it is the design↔analysis contract.
9. **Two distinct barcode algorithms** — `barcode_type` `0`/`terminus` (fast, distinctive 5′/3′ ends) vs
   `1`/`spectrum` (thorough, k-mer saturation).
10. **GC content is reported (not constrained) for primers** — `primer_guanine_cytosine_content` in stats.
11. **The absence of any row-union / sub-library concatenation operator.** Across all 22 modules there is
    no way to vertically combine two libraries with different `ID`s into one pool. `merge` concatenates
    *columns*; `join` requires identical `ID` sets and "never creates or drops rows". **This is the single
    cleanest, most defensible statement of what PoolParty does that Oligopool Calculator does not** — much
    safer than any "no DAG" or "no composability" framing. Recommend making it the load-bearing sentence.
12. **Runnable Design-Mode DAG example** (`mpra_design_parallel.yaml`) and the parser's topological
    primer-dependency resolution (`get_primer_order`) — see §7.

---

## Bottom line

- **0 values are wrong.** The extraction's capability grid is honest and, if anything, slightly generous
  to Oligopool Calculator in places (`assay_dms = partial`, `library_as_object = partial`) — which is the
  safe direction.
- **4 evidence strings need rewriting before use** (`synthesis_constraints` — fabricated quote;
  `automatic_naming` — false claim about column names; `mixed_mutagenesis_one_pool` — falsifiable absolute;
  `maintained`/`availability_status` — "two releases in four days" is factually wrong).
- **2 quotes are paraphrases masquerading as verbatim** and should be corrected.
- **The biggest unexamined artefact is `examples/library-compressor/mutant_generator.py`.** Read it before
  writing the referee response.
- **Recommended load-bearing contrast sentence:** *"Oligopool Calculator has no operation that combines
  two variant sets into one pool — `merge` concatenates columns and `join` requires identical `ID` sets and
  'never creates or drops rows' — so a library is always exactly the rows the user supplied."* This is
  verifiable from the shipped API reference and cannot be rebutted by pointing at an example script.
