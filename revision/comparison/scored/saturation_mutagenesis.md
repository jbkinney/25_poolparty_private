# Saturation mutagenesis capability audit

> ## CORRECTION applied after scoring — three cells demoted
>
> `ROW_DEFINITIONS.md` gained a **global rule**: a capability counts only where the
> tool provides an operation for it, never where the user reconstructs it. The
> `partial` branch this row previously offered — *"exhaustiveness is the user's
> responsibility"* — was struck as a consequence.
>
> Three cells were scored on that branch **by name**, so their demotion follows
> deterministically from their own stated reasoning. No rescoring was performed and
> no other cell was touched.
>
> | Tool | Was | Now | The rationale it rested on |
> |---|---|---|---|
> | MPRAnator | partial | **no** | *"meets the definition's 'exhaustiveness is the user's responsibility' branch"* |
> | MPRA Design Tools | partial | **no** | *"the consistent score is `partial`, specifically 'exhaustiveness is the user's responsibility'"* |
> | Mutation Maker | partial | **no** | *"exhausting the eligible region remains the user's responsibility"* |
>
> **Still unresolved — needs targeted rescoring, not demotion:** Oligopool
> Calculator and DNA Chisel are `partial` on other grounds (an `examples/`-only
> helper, and a Cartesian space requiring user filtering). Neither cites the struck
> branch, so neither can be demoted by inspection.
>
> The per-tool sections below retain their original headings and reasoning as
> written at scoring time; the Results table has been updated.

Audit date: 2026-08-15. Measurement instrument: `revision/tables/ROW_DEFINITIONS.md`, row 3, read before any tool-specific source. The per-tool records under `tool_survey/final/` were used only to locate canonical repositories/documentation; no statement, rating, or conclusion in those records is evidence here.

## Operational test

Given one documented design specification containing a parent sequence and an eligible region, does the tool itself emit or evaluate every one-event substitution—every non-reference allowed alternative at every eligible position—in one run (`yes`), with `partial` reserved for a position-exhaustive but narrowly fixed substitution, a scan whose exhaustive sites/alleles or Hamming-distance-one filtering remain the user's responsibility, or support confined to an example rather than the documented/exported API?

This test intentionally does **not** credit random mutagenesis, a tool merely consuming pre-made variants as full support, a full Cartesian multi-position expansion as an exact single-mutant scan, or the enumeration of primer candidates rather than substitution variants.

## Results

| Tool | Score | Confidence | Short reason |
|---|---|---|---|
| PoolParty | **yes** | high | Documented `mutagenize(..., num_mutations=1, mode="sequential")` covers all positions and all non-WT bases; reproduced as 12 variants for a 4-mer. |
| VaLiAnT | **yes** | high | Documented `snv` action replaces each nucleotide by all three alternatives throughout the specified target region. |
| MPRAnator | **no** [corrected] | high | It applies every ALT in a user-authored SNP table, but it does not derive all positions/alternatives from a region; its other mutator is random. |
| MPRA Design Tools | **no** [corrected] | high | It expands ALT alleles already present in a user-authored VCF and constructs ref/alt MPRA sequences; it does not generate the exhaustive VCF. |
| Oligopool Calculator | **partial** | high | An official example helper generates exactly all single mutants, but that helper is under `examples/` and is absent from the documented/exported package API. |
| Mutation Maker | **no** [corrected] | high | SSSM saturates user-listed amino-acid sites with a degenerate codon; choosing every eligible position is explicitly the user's responsibility and concrete variants are represented degenerately. |
| DNA Chisel | **partial** | medium | Documented `MutationSpace.all_variants()` exhausts a Cartesian mutation space, but an exact single-substitution library requires the user to supply/filter the space to Hamming distance one; the normal optimizer returns only the best sequence. |
| tangermeme | **yes** | high | Its documented `saturation_mutagenesis` builds every true single-base edit in the requested `[start,end)` region and exposes per-edit predictions. |

## Per-tool audit records

### PoolParty — **yes** (high confidence)

**Positive evidence.** The documented API states:

> “With `mode='sequential'`, every single-base variant within the region is enumerated.”

Source: canonical local repository at commit `1bb0179e1c3720b1fffd471802b3040f9336de28`, `poolparty/docs/operations/mutagenize.rst:147-158`.

It states even more explicitly:

> “`mode='sequential'` with `num_mutations=1` enumerates every single-point variant in deterministic order, covering all positions and non-wild-type bases.”

Source: `poolparty/docs/operations/mutagenize.rst:201-210`. The example's reported size is 12 for `ACGT` (`:213-220`), exactly 4 positions × 3 alternatives.

The implementation builds the product of position combinations and substitution indices. For the uniform DNA case it computes

```python
num_combinations = comb(num_positions, self.num_mutations) * (
    alpha_minus_1**self.num_mutations
)
```

and iterates every `positions` combination and every `mut_pattern` (`poolparty/src/poolparty/base_ops/mutagenize.py:332-379`). At generation time the chosen cache entry is mapped back to positions and alternative characters (`:563-584`).

**Behavioral check required for this repository.** With `PYTHONDONTWRITEBYTECODE=1` and the read-only Python 3.12 environment:

```text
pp.from_seq("ACGT").mutagenize(num_mutations=1, mode="sequential")
num_states = 12
CCGT GCGT TCGT AAGT AGGT ATGT ACAT ACCT ACTT ACGA ACGC ACGG
```

All sequences are unique and are precisely the 12 Hamming-distance-one substitutions.

**Disconfirmation attempt / alternative values.** A lower value would have required evidence that sequential mode sampled, omitted some non-WT bases/positions by default, required a caller loop, or existed only in a tutorial. I searched the operation documentation, the DMS tutorial, implementation, and tests for `mutagenize`, `sequential`, `substitution`, `allowed_chars`, `num_mutations`, and region restrictions. The implementation and live run confirm the full default alphabet and all positions. A user may deliberately narrow eligibility with `region` or `allowed_chars`, but that is the specification of eligible positions/alternatives, not an unexplained omission. No disconfirming restriction was found.

### VaLiAnT — **yes** (high confidence)

**Positive evidence.** The canonical README defines the documented `snv` action:

> “Each nucleotide is replaced with all the alternatives.”

It then shows target `AA` producing `CA, GA, TA, AC, AG, AT`.

Source: `https://github.com/cancerit/VaLiAnT`, branch `develop`, commit `8796cc112dafd4919fec59913f58cd2be87c45eb`, `README.md:294-307`.

The target specification is a single row containing target-region coordinates and an action vector; the README calls it a file describing “the types of mutation to be applied to the three target regions” and documents `action_vector` as mutation labels grouped by region (`README.md:630-656`). For cDNA, one row similarly contains the target region and action labels (`:659-679`).

The source independently establishes completeness. `SnvMutator.get_variants` obtains all one-nucleotide windows and flattens `Variant.get_snvs(...)` over them (`src/valiant/mutators/snv.py:37-49`). `Variant.get_snvs` is documented “Create the three SNV's for a given position” and constructs one variant for every member of `NT_SNVS[nt]` (`src/valiant/variant.py:107-113`); `NT_SNVS` is `NTS - {nt}` over `ACGT` (`src/valiant/constants.py:19-24`).

The authors' paper is direct corroboration:

> “`snv` results in all possible nucleotide substitutions at each position in the targeton”

Source: Barbon et al. 2022, supplied `Barbon2022_VaLiAnT_all.pdf`, p. 3, Fig. 2 caption.

**Disconfirmation attempt / alternative values.** `partial` would have applied if `snv` meant one fixed alternative, only CDS positions, or a caller-authored VCF; `no` would have applied if VaLiAnT merely tiled oligos around supplied variants. I searched the complete README, examples, mutator classes, loaders, target-region code, and paper for `saturat`, `mutator`, `snv`, `all alternatives`, `all possible`, `each nucleotide/position`, `substitut`, and custom VCF handling. The non-CDS `snv` mutator is built-in, action-vector driven, and covers all three alternatives per base. Custom VCF support is a separate mode and does not limit `snv`. No contrary restriction was found.

### MPRAnator — **partial** (high confidence)

**Positive but limited evidence.** The SNP module supports all user-listed ALT alleles in one submission. Its official downloadable client makes the caller provide the variants directly:

```python
snpS = "1  42991  rs43434  T  G\n 1  43030  rs121212  A  G"
...
data = {"sequenceS": sequenceS, "SnpS": snpS, ...}
```

Source: `https://github.com/hemberg-lab/MPRAnator`, commit `9969790d62410138d4281b7955da6d085f07b1bc`, `downloadables/MpraSnps_script.py:12-15,41-55`.

The server parses that user field as `SnpFile` and passes it to `make_sequence_copies` (`iliasApp/views.py:319-363`). In the substitution path, position, REF, and ALT are read straight from each input record; comma-separated ALT values are split and each is emitted (`parseSNPs.py:137-158,171-183`). The documentation describes the page only as studying supplied SNPs and optionally combining them:

> “This page allows the user to synthesize oligonucleotides for MPRA experiments to study the effects of SNPs. The user can select to include or exclude combinations of SNPs...”

Source: `iliasApp/templates/iliasApp/docs.html:40-43`.

This can realize an exhaustive single-substitution scan if the caller authors every position and all three ALT alleles, so it meets the definition's “exhaustiveness is the user's responsibility” branch, but not the one-specification automatic scan required for `yes`.

**Disconfirmation attempt / alternative values.** `yes` would require a documented region-level switch that derives positions and all non-reference alleles; `no` would require absence of even a user-supplied exhaustive route. I searched every Python function, the Django forms/views/templates, downloadable scripts, README, and supplied paper for `saturat`, `mutagen`, `mutation`, `SNP`, `substitut`, `variant`, `position`, `all possible`, and `every/each position`. No region-wide exhaustive SNP generator exists. The only sequence-derived mutator is `part3.mutateString`, which selects positions and alternatives with `random.choice` (`part3.py:13-48`), and the docs call it “Mutate random positions” (`docs.html:138-159`). That rules out `yes`, while the explicit multi-ALT SNP application rules out `no` under the instrument's user-responsibility clause.

### MPRA Design Tools — **partial** (high confidence)

This record covers the authors' `mpradesigntools` package and its companion `designMPRA` repository.

**Positive but limited evidence.** The package's documented primary operation is:

> “`processVCF` takes a VCF of SNPs ... and turns them into a set of labeled MPRA sequences barcoded with inert n-mers”

Source: `https://github.com/andrewGhazi/mpradesigntools`, commit `afd386ef12051bb0a37ad63a6f92acd555246584`, `man/processVCF.Rd:78-87`.

The README makes clear that position and alleles come from the user's VCF:

> “Only the CHROM, POS, REF, and ALT columns are used.”

and

> “Multiple alternate alleles should be separated in the ALT field by a comma and no spaces”

Source: `README.md:50-59`.

The code spreads comma-delimited ALT values into rows (`R/processVCFfast.R:5-16`), identifies an SNV solely from the row's `REF` and `ALT` (`:210-234`), and constructs the alternate sequence with `replaceLetterAt(..., snp$ALT)` (`:327-386`). Thus a caller-authored exhaustive VCF is processed, but the tool does not derive all eligible positions or all alternatives from a parent region.

**Disconfirmation attempt / alternative values.** `yes` would require a documented/exported function accepting a parent/region and generating all position × alternative substitutions. I searched all package R sources, generated `.Rd` documentation, README, `NAMESPACE`, the complete companion `designMPRA` scripts/server/UI, and the supplied paper for `saturat`, `mutagen`, `mutation`, `substitut`, `variant`, `snv`, `snp`, `allele`, `position`, `every position`, and `all possible`. `NAMESPACE:3-4` exports only `processVCF` and `spread_and_fix_indels`; the function inventory contains VCF processing, barcoding, and assay-power code, not a saturation generator. `no` would require inability to realize an exhaustive list at all, contradicted by comma-separated ALT expansion and arbitrary VCF rows. Therefore the consistent score is `partial`, specifically “exhaustiveness is the user's responsibility.”

### Oligopool Calculator — **partial** (high confidence)

**Positive but limited evidence.** The canonical repository includes an author-documented example helper whose docstring says:

> “Generate all single-nucleotide mutants of a DNA sequence. For a sequence of length L, this generates up to 3*L mutants.”

The code loops `for pos in range(len(sequence))`, then `for alt_base in DNA_BASES`, excluding the reference base and yielding each mutant (`examples/library-compressor/mutant_generator.py:43-79`).

Source: `https://github.com/ayaanhossain/oligopool`, commit `b88fa394ca67ed4c48ec41127b5707694ee7cf0a`.

A no-install behavioral check of this helper on `ACGT` returned the WT plus exactly 12 unique single mutants, matching 4 × 3. The example README explicitly labels `mutant_generator.py` a “Helper utility for generating test variants” (`examples/library-compressor/README.md:5-9`) and presents saturation as a compressor use case, not as an exported design function (`:69-85`).

This exact capability therefore exists, but only under `examples/`; the measurement instrument caps example-only support at `partial`.

**Disconfirmation attempt / alternative values.** `yes` would require the same generator in the documented/exported `oligopool` API. I searched the full source, docs, README, notebooks, and supplied paper for `saturat`, `mutagen`, `mutation`, `single mutant`, `substitut`, `variant`, and `all/every position`, then checked the package export table and API table of contents. `oligopool/__init__.py:6-31,42-69` exports design, assembly, degenerate, and analysis functions but no mutant generator. The documented `compress` API requires a pre-built DataFrame of concrete sequences (`docs/api.md:770-819`), while `expand` expands a caller-supplied IUPAC sequence into its full Cartesian space (`:834-874`), not an exact one-edit scan over a multi-position region. `no` would require absence even from examples; the exact helper and live 12-mutant result disprove that. The 2024 paper contains no saturation-mutagenesis API; its `saturates` usage concerns barcode k-mer coverage, not variant scanning.

### Mutation Maker — **partial** (high confidence)

**Positive but limited evidence.** The authors' paper describes the SSSM workflow and its input boundary precisely:

> “The minimal input for both SSSM algorithms is the plasmid sequence containing the parental gene of interest (GOI), the forward and reverse flanking primers and a list of mutation sites with their respective degenerate codon (default, ‘NNK’).”

Source: Hiraga et al. 2021, supplied `Hiraga2021_MutationMaker.pdf`, p. 3.

The same paper says target sites may be selected rationally “or alternatively through an exhaustive scan of the entire sequence,” but Fig. 1's caption identifies the actual user input as “a list of amino acid positions to be changed” (p. 2). Thus saturation at each listed site is supported, while exhausting the eligible region remains the user's responsibility. It is also a narrow codon/amino-acid event represented by a degenerate codon rather than a concrete enumeration.

Repository source matches the paper: `SSMInput` requires `mutations = ListProperty(str, required=True)` and parses every supplied string; its degenerate codon defaults to `NNS` (`backend/mutation_maker/ssm_types.py:211-218`). The UI explicitly requires users to enter mutations in `[Amino Acid Codon][Location][X]` form (`frontend/src/scenes/SSM/components/SSMForm.tsx:217-235`). The `/ssm` endpoint is a documented task type (`api/api.py:40-45`).

**Disconfirmation attempt / alternative values.** `yes` would require a region/full-gene option that automatically derives all mutation sites and concretely emits each single alternative; `no` would require no exhaustive-per-site mechanism. I searched the README, all backend and API Python, frontend SSM form/result code, feature tests, and paper for `site saturation`, `site-scanning`, `SSM/SSSM`, `mutation input/list`, `all/every position`, `amino acid`, `degenerate codon`, and `substitution`. No automatic “all positions” input switch was found: both paper and code require a mutation-site list. Conversely, NNK/NNS saturation and the SSSM endpoint are explicit. This is the definition's `partial` case (user-managed regional exhaustiveness and a narrow codon event), not `yes` or `no`.

### DNA Chisel — **partial** (medium confidence)

**Positive but limited evidence.** `MutationSpace.all_variants(sequence)` is documented in the API through the `MutationSpace` automodule (`docs/ref/core_classes.rst:45-49`) and says:

> “Iterate through all sequence variants in this mutation space.”

Its implementation constructs variant choices for each mutable segment and iterates `itertools.product(*variants_slots)` (`dnachisel/MutationSpace/MutationSpace.py:132-164`). A default problem initially gives every nucleotide position all four DNA choices (`:166-180`). Therefore every Hamming-distance-one substitution is present in the Cartesian enumeration.

However, the method also emits the WT and all higher-order combinations. To obtain **exactly** the single-substitution scan measured here, the caller must filter yielded strings by Hamming distance one or hand-construct a mutation space containing only those region variants. That is not a documented saturation operation; it is user-managed extraction from a more general enumerator. The normal documented workflow instead solves an optimization problem and returns one final sequence (`README.rst:44-82`). Internally, exhaustive optimization loops over `all_variants` only to retain `current_best_sequence`, then assigns that one value back to `self.sequence` (`DnaOptimizationProblem/mixins/ObjectivesMaximizerMixin.py:26-60`).

**Why `partial`, not `yes` or `no`.** Under the instrument, the documented exhaustive generator is stronger than absence: a caller can obtain the exact scan without inventing substitution chemistry, but is responsible for restricting the exhaustive Cartesian output to one edit. `yes` would require a documented call/parameter producing only all single substitutions for a declared region. `no` would require no exhaustive sequence-variant enumerator. I searched the entire source, docs, README, examples, and paper for `saturat`, `mutagen`, `single mutation/substitution`, `all_variants`, `MutationSpace`, `pick/apply_random_mutations`, `every/each position`, and `EnforceChanges`. There are zero `saturat` or `mutagen` hits in source/docs and no scan tutorial/example; the only exhaustive primitive is the general Cartesian `all_variants`, while `pick_random_mutations(n_mutations=1)` is sampled (`MutationSpace.py:106-130`). The paper likewise says exhaustive search through “all possible sequence variants” is an optimizer search mode, not a saturation-library output (supplied `Zulkower2020_dnachisel.pdf`, p. 2). A no-install behavioral import was attempted but stopped at the missing runtime dependency `Bio`; per instruction nothing was installed. Source control flow is unambiguous, but the boundary between a general Cartesian enumerator plus caller filtering and a saturation API is interpretive, hence medium rather than high confidence.

### tangermeme — **yes** (high confidence)

**Positive evidence.** The canonical README defines the operation as comparing the original with sequences that:

> “comprehensively each contain one mutation with respect to that original sequence.”

Source: `https://github.com/jmschrei/tangermeme`, commit `2006b310cd72a28c56c3ea4ba67f738fff74bb89`, `README.md:224-234`.

The documented API docstring says it returns predictions on the original and:

> “each of the sequences with an edit distance of one on them.”

Source: `tangermeme/saturation_mutagenesis.py:81-102`, included in the API reference by `docs/api/saturation_mutagenesis.rst:1-5`. `start` and `end` declare the eligible region (`:121-131`).

Implementation confirms exact enumeration. For each example it identifies every non-identity character/position pair (`edit_chars, edit_positions = torch.where(~identity)`), repeats the parent once per edit, zeroes that position and sets the alternative character (`:212-232`), predicts every edit, and scatters results into the full `[character, position]` grid (`:240-269`). For clean one-hot DNA there are three true edits per eligible position; identity slots are filled from the reference prediction rather than redundantly evaluated. `raw_outputs=True` returns the per-edit predictions in `SaturationMutagenesisRawResult` (`:19-30,180-187,269`).

**Disconfirmation attempt / alternative values.** `partial` would apply if the operation sampled edits, handled only one fixed alternative, depended on a user loop, or existed only in a notebook; `no` would apply if “saturation mutagenesis” were merely a plotting label. I searched the README, source, API docs, tutorial inventory, skill reference, tests, and supplied paper for `saturation_mutagenesis`, `every possible`, `each position`, `all substitution`, `edit distance`, `n_edits`, and start/end handling. Source enumerates all true alternatives over the requested interval in one call, and API docs expose it. The paper also lists and uses windowed saturation mutagenesis (supplied `Schreiber2025_tangermeme_all.pdf`, pp. 2, 7). No contrary sampling or example-only restriction was found.

## Search and verification ledger

This section records every evidence search performed, including negative searches. Commands are normalized for readability; paths were rooted at `/mnt/c/Users/zhliu/Desktop/KinneyLab/` (case-equivalent mounted path used by the shell).

1. **Instrument check before tool review:** `sed -n '1,240p' 25_poolparty_private/revision/tables/ROW_DEFINITIONS.md`. The file existed; row 3 was used verbatim.
2. **Lead-only lookup:** listed `tool_survey/final/<slug>.md` and used repository/documentation URL and source-filename strings only to locate primary material. These records were not cited and their ratings/conclusions were disregarded.
3. **Primary-source acquisition:** shallow-cloned the seven external canonical repositories plus the separate `designMPRA` companion repo into `/tmp/saturation-audit.H0l57X/`; recorded each `git rev-parse HEAD`. PoolParty was read only from `poolparty-statecounter/` and its supplied PDF.
4. **Uniform first-pass repository grep:** for every clone, recursively searched `saturat|mutagen|substitut|single[- ]nucleotide|all variants|every position|each position|alanine scan|deep mutation|DMS`, excluding lock/minified/map files, and inventoried matching files.
5. **PoolParty searches:** recursively searched package source, docs, tutorials, tests, and examples for `sequential|exhaust|mutagen|substitut|every position|each position|all possible`; then narrowed to `mutagenize.py`, `docs/operations/mutagenize.rst`, and `docs/tutorials/dms_gb1.rst` with `def mutagenize|mode.*sequential|substitution|mutation_type|num_mutations|_build_caches`. Inspected source around lines 320-390 and 550-590 and docs around lines 147-220. Ran the 4-mer behavior check in the mandated read-only venv with bytecode disabled.
6. **VaLiAnT searches:** recursively searched README, examples, and `src/` for `saturation|snv|single nucleotide|mutator|all amino|all.*substitut|all (three|possible)|every/each base/position`. Inspected mutation-type docs, SGE/cDNA target formats, `mutators/snv.py`, `variant.py`, `constants.py`, `mutators/__init__.py`, and `mutator.py`.
7. **MPRAnator searches:** inventoried the complete repository tree and every top-level/Django Python function. Recursively searched all Python, HTML, and Markdown for `saturat|mutagen|mutation|snp|substitut|transmut|random|position|variant`. Inspected `docs.html`, `MpraSnps_script.py`, `views.py`, `parseSNPs.py`, `part3.py`, forms, templates, and README.
8. **MPRA Design Tools searches:** inventoried both repository trees; recursively searched `.R`, `.Rmd`, `.Rd`, `.md`, and `.html` for `saturat|mutagen|mutation|substitut|variant|snv|snp|allele|position|sequence`. Inspected README, `processVCF.Rd`, `R/processVCFfast.R`, companion `scripts/processVCFfast.R`, every R function definition, and `NAMESPACE` exports.
9. **Oligopool Calculator searches:** recursively searched source, docs, README, examples, and notebooks for `saturat|mutagen|mutation|single.?mut|substitut|variant|degenerate|compress|all.*position`; counted hits as a cross-check. Inspected the example README and both example scripts, documented API sections for `compress`/`expand`, and the package export/lazy-import table. Ran `generate_single_mutants('ACGT')` without importing/installing the core package; observed WT + 12 exact mutants.
10. **Mutation Maker searches:** recursively searched README, API/backend Python, frontend TypeScript/TSX, feature files, tests, and docs for `site.?saturation|saturation mutagenesis|single.?site|mutation.*input|input.*mutation|mutations|amino acid|SSM|SSSM|scan|every position|all positions|all.*substitut|degenerate_codon`. Inspected `/ssm`, `SSMInput`, `ssm_solve`, the SSM UI mutation field, and paper pp. 2-3.
11. **DNA Chisel searches:** inventoried docs; recursively searched source, docs, README, examples, and tests for `all_variants|mutation.?space|mutagen|saturat|single.?mutation|single.?substitut|all.*substitut|every position|each position|pick_random_mutations|apply_random_mutations|EnforceChanges`, plus term-count cross-checks (zero `saturat` and zero `mutagen` hits). Inspected `MutationSpace.py`, core-class docs, README workflow, `DnaOptimizationProblem`, optimizer mixins, `DnaNotationPattern`, built-in-specification index, and examples. Attempted a no-install runtime check; import failed only because `Bio` was unavailable, and no dependency was installed.
12. **tangermeme searches:** recursively searched README, source, docs, tests, tutorial inventory, and bundled reference docs for `saturation mutagenesis|saturation_mutagenesis|every possible|each position|all.*substitut|edit distance|n_edits`. Inspected the full API docstring and implementation, API toctree page, README example, tests, and window/start/end behavior.
13. **Paper searches:** using `python3 -c "import fitz; ..."` as required, extracted all eight supplied PDFs and searched case-insensitively for `saturat`, `mutagen`, `substitut`, `single nucleotide`, `every position`, `all possible`, `deep mutational`; for the three MPRA/oligopool papers also searched `SNP`, `variant`, and `input`. Mutation Maker pp. 2-3 were then extracted in full to resolve whether the tool or user selects exhaustive sites.
14. **Repository-state checks:** recorded commit hashes and checked the PoolParty worktree without modifying it. Existing PoolParty changes/build artifacts were left untouched. No packages were installed and no source repository was edited.

## Row-level finding

The row is applicable on one scale and does discriminate, but only into two occupied levels in this set: three tools have a documented automatic exact scan (`yes`), while five expose a materially real but restricted route (`partial`). There are no `no` or `unknown` cells. The partial group is heterogeneous—user-authored variant/site lists (MPRAnator, MPRA Design Tools, Mutation Maker), example-only exact generation (Oligopool Calculator), and caller-filtered Cartesian enumeration (DNA Chisel)—so the row distinguishes automatic saturation well but compresses several different limitations into the definition's single `partial` category.
