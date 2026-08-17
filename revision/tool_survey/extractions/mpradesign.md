# MPRA Design Tools (Ghazi et al. 2018) — evidence memo

**Slug:** `mpradesign`
**Full name:** MPRA Design Tools (Shiny app `designMPRA` + R package `mpradesigntools`)
**Citation:** Ghazi AR, Chen ES, Henke DM, Madan N, Edelstein LC, Shaw CA. "Design tools for MPRA experiments." *Bioinformatics* 34(15):2682–2683 (2018). doi:10.1093/bioinformatics/bty150

---

## 1. Sources consulted

| Kind | Reference |
|---|---|
| PDF (paper + full supplement, 28 pp) | `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/papers/Ghazi2018_MPRADesignTools_all.pdf` |
| Prior analysis | `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/prior_analyses/Ghazi2018_MPRADesignTools_all_analysis.md` |
| Repo (R package) | https://github.com/andrewGhazi/mpradesigntools |
| Repo (Shiny app) | https://github.com/andrewGhazi/designMPRA |
| Source read verbatim | `R/processVCFfast.R` (1459 lines), `DESCRIPTION`, `NAMESPACE`, `README.md`, `man/processVCF.Rd`, `man/processSnp.Rd`, `man/twelvemers.Rd`, `designMPRA/ui.R`, `designMPRA/server.R` |
| Web service | https://andrewghazi.shinyapps.io/designmpra/ — **HTTP 404, "Not Found"** (checked 2026-08-10) |
| Docs mirror | https://rdrr.io/github/andrewGhazi/mpradesigntools/ |
| CRAN | Not present. rdrr.io indexes it under `/github/`; confirmed "GitHub only (not on CRAN)". |

Notable: the prior analysis says "no DMS support, no motif grammar, no combinatorial design, SNP-focused only, R/Shiny." All of that is **confirmed** by reading the source. The prior analysis understates two things: (a) the barcode machinery is more substantial than a one-liner (28 barcode sets incl. indel-error-correcting `freebarcodes` sets), and (b) the tool does real construct-level restriction-site checking and repair. It also does not mention that the web app is now dead.

---

## 2. What the tool actually is

Two artifacts:

1. **`designMPRA`** — an R/Shiny web app with three tabs: "MPRA" (static explanatory text + two PNG diagrams), "Power" (interactive `pwr.t.test` power curve), "Create sequences" (upload VCF → download TSV). Last commit **2017-09-26**. The hosted app is **dead (404)**.
2. **`mpradesigntools`** — the companion R package. Version 0.2.0. Last commit **2020-12-17**. Two exported functions: `processVCF()` and `spread_and_fix_indels()` (plus 28 barcode datasets).

The headline scientific contribution is the **power analysis**, not the sequence design. Abstract: *"an online web-tool and R package that allows for interactive MPRA experimental design encompassing both power analysis and design of constructs."* The app itself says of its own sequence generator: *"we would like to emphasize that this sequence generation functionality of the web tool is mainly for simple demonstrative purposes"* (`ui.R`).

### The one design operation

`processVCF()` is a single monolithic function. Signature (`man/processVCF.Rd`):

```r
processVCF(vcf, nper, upstreamContextRange, downstreamContextRange,
           fwprimer, revprimer,
           enzyme1 = "GGTACC", enzyme2 = "TCTAGA", enzyme3 = "GGCCNNNNNGGCC",
           filterPatterns = "AATAAA", alter_aberrant = FALSE, extra_elements = FALSE,
           max_construct_size = NULL, barcode_set = "twelvemers",
           ensure_all_4_nuc = TRUE, flip_RV = TRUE, outPath = NULL)
```

Semantics: *"takes a VCF of SNPs (preferably from dbSNP) and turns them into a set of labeled MPRA sequences barcoded with inert n-mers."* For each variant it pulls hg38 context from `BSgenome.Hsapiens.UCSC.hg38`, builds ref and alt versions, and emits `nper` barcoded copies of each.

Fixed construct layout, hard-coded in `processVCFfast.R`:

```r
sequence = paste0(fwprimer, constrseq, enzyme1, enzyme2, barcodes, revprimer)
# with extra_elements = TRUE:
sequence = paste0(fwprimer, 'TG', constrseq, enzyme1, enzyme2, barcodes, 'GGC', revprimer)
```

Supplement S1.3: *"enzyme1 and enzyme2 default to KpnI and XbaI as in Melnikov et al., Nature Biotechnology 2012. The forward and reverse primers default to ACTGGCCAG and CTCGGCGGCC respectively."*

### Output schema

`processVCF()` returns `list(result = <tibble>, failed = <tibble>)`. The result columns are exactly (from the `select()` at the end of `processVCF`):

```r
select(ID, type, allele, snpIndex, totIndex, barcode, sequence, site_fix_info, notes)
```

`failed` carries a per-SNP `reason` string, e.g. `'Failed - Context contained a digestion site'`, `'Failed - the variant is located within an aberrant digestion site'`, `'Failed - SNP sequence could not be generated without an aberrant digestion site even after barcode resampling'`.

---

## 3. Capability-by-capability

### Block A — library specification

**`library_as_object` = partial.** One `processVCF()` call materializes the whole library, so the user does not loop over variants themselves — but the product is a plain tibble (`list(result = <tibble>, failed = <tibble>)`), not a library object with a type or methods. There is no S3/S4 class, no `NAMESPACE` export of any class, and no operation that consumes a library and returns a library. The only compositional pairing in the package is the preprocessing helper `spread_and_fix_indels()` → writes `*_fixed.vcf` → feed to `processVCF()`.

**`dag_chaining` = no.** Two exported functions total (`NAMESPACE`: `export(processVCF)`, `export(spread_and_fix_indels)`). The design is a single fixed pipeline hard-coded inside `processSnp()`; there is no way to compose, nest, or reorder design steps. The only user-visible "chaining" is the file handoff described above, which is a VCF-normalization prestep, not composition of design operations.

**`lazy_evaluation` = no.** Everything is materialized eagerly. The entire barcode set is loaded into memory and filtered up front (`mers = get(barcode_set)`, then `mers[purrr::map_lgl(mers, matches_all_nucleotides)]`), the whole VCF is read to a tibble, and all constructs are built and `rbind`-reduced. The README quantifies the cost: *"5 + .01 * Number of barcodes per allele * Number of SNPs in VCF * 2."* `processVCF` hard-errors if the design exceeds the materialized barcode pool: `'Your design requires more barcodes than is possible with the selected barcode_set. Try a bigger set.'`

**`mixed_mutagenesis_one_pool` = no.** The only sequence variation produced is ref vs. alt for each VCF row: `type = rep(c('ref', 'alt'), each = nper)`. There is no mutagenesis scheme at all — no exhaustive singles, no pairs, no random sampling, no scanning. Multi-allelic and indel variants are supported (`spreadAllelesAcrossRows`, `fix_indels`), but that is variant parsing, not mutagenesis design. The ref allele always accompanies the alt, which is the closest thing to a WT control, but the user cannot specify a different mix.

**`combinatorial_motif_place` = no.** The word "motif" does not appear in the package. There is no concept of an element, a position, an orientation, or a placement. The only positional control is symmetric flanking-context width (`upstreamContextRange`, `downstreamContextRange`) around the variant, and the only orientation control is a per-variant strand flag (`MPRAREV` / `RV` tags in the VCF INFO field).

**`barcode_generation` = yes.** This is genuinely a strength. Supplement S5 documents the original set: *"We generated the set of all possible DNA 12-mers then screened these according to the design parameters in Melnikov et al. ... each nucleotide occurs at least once; no runs of single nucleotides greater than length 3; do not contain restriction sites for KpnI/XbaI; do not start with TCT (this creates a restriction site with XbaI with the preceding sequence); they do not match any human miR seed sequences."* Yielding 1,140,292 inert twelvemers (`man/twelvemers.Rd`: "1140292 inert twelvemers"). Version 0.2.0 adds 27 further sets from `freebarcodes` (Hawkins et al., PNAS 2018) indexed as `barcodes<length>-<n_correctable_errors>`, i.e. **edit-distance-constrained** barcodes (10–17 bp, 1–2 indel errors correctable). At runtime `processVCF` re-filters for all-four-nucleotides (`ensure_all_4_nuc`) and removes barcodes matching `filterPatterns` (default `'AATAAA'`, the polyA signal) plus all three enzyme sites and their reverse complements/reverses. Barcodes are then shuffled, partitioned into per-SNP pools, and sampled without replacement (`sample(snp$bcPools %>% unlist, 2*nper, replace = FALSE)`) and pasted into the construct. Custom barcode vectors are accepted. **Caveat:** there is no GC-content constraint anywhere.

**`per_sequence_provenance` = partial.** Each output row carries structured, machine-readable build metadata beyond the mutation name: `type` (ref/alt), `allele`, `barcode`, `snpIndex`, `totIndex`, plus two genuine provenance fields — `site_fix_info` (a nested data frame recording exactly which bases were changed to destroy an aberrant restriction site, unnested into the TSV) and `notes`. Example note strings from the source: `'Shortened context by ', amount_to_remove, ' bp on each side to account for input maximum construct size. '` and `'The alleles for this SNP were flipped from the input VCF because of the presence of the RV tag in the INFO field.'` But the schema is a fixed flat set of columns for one design operation; there is no general record of "how it was built" because there is only one way anything is ever built.

**`automatic_naming` = partial.** The docs say `processVCF` produces *"labeled MPRA sequences"* and `processSnp` returns *"a data_frame of labeled sequences with appropriate information on the changes made."* But the labelling is a set of identifier **columns** (`ID` from the VCF, `type`, `snpIndex`, `totIndex`), not a composed name string. No name-generating code (no `paste0` producing an identifier) exists in `processVCFfast.R`. If you export FASTA you must build names yourself.

**`design_visualization` = no.** Searched `designMPRA/server.R`: the only outputs are `powerPlot` (a `pwr.t.test` power curve — a statistics plot, not a view of the library), `timeText`, `warnText`, `conditionalDownload`, `failed` (a table of failed SNPs) and `inputHead` (a table of the first few VCF rows). The construct-layout picture in the "Create sequences" tab is a **static PNG** (`img(src = 'MPRAseqDiagram.png', width = 955, height = 20)`), identical regardless of the design. The R package has no plotting function at all (no `graphics`/`ggplot2` import in `NAMESPACE`); the only visual feedback is `print()`ing a construct-size distribution table.

### Block B — assay coverage

**`assay_dms` = no.** No codon, amino-acid, ORF or reading-frame concept anywhere in the package. `NAMESPACE` imports only `Biostrings` nucleotide utilities (`replaceLetterAt`, `complement`, `reverseComplement`, `subseq`, `countPattern`); nothing translational.

**`assay_mpra` = yes.** This is the tool's entire purpose. Title: "Design tools for MPRA experiments." `DESCRIPTION`: `Description: Design MPRA sequences from input VCF's`. Scope is narrow within MPRA: SNP/indel allelic-contrast MPRAs with the Melnikov et al. 2012 oligo layout.

**`assay_insilico` = no.** No model, no scoring, no ML dependency; nothing in the paper, README, or source mentions predictive models.

### Block C — genomics integration

**`genome_coordinates` = yes.** VCF `CHROM`/`POS` are used to slice the reference directly: `snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], ...)`, where `genome` is `BSgenome.Hsapiens.UCSC.hg38` (a hard `Depends`/`import` in `DESCRIPTION`/`NAMESPACE`). The paper explicitly frames this as its advantage over MPRAnator: *"Our tool acquires genomic context from the hg38 reference genome, rather than requiring input by the user."* **Constraint:** hg38 human only — the README's "Planned Features" list still has an unchecked *"mm10 genomic context"*.

**`transcript_models` = no.** No GTF/GFF anywhere; no genome-annotation dependency in `DESCRIPTION`. Strandedness is not derived from annotation — the user must hand-annotate it: *"If the user wishes to generate SNPs for genes that normally are read from the reverse strand, add a string containing 'MPRAREV' to the INFO field of the VCF"* and *"the MPRAREV tag will need to be added by the user (where appropriate) because the VCF's do not always specify which strand the relevant gene is on."*

**`exon_intron_split_codons` = no.** Follows from the two above: no transcript structure, no codons.

**`hgvs_input` = no.** The only accepted input is a VCF path. README: *"Only the CHROM, POS, REF, and ALT columns are used. The INFO column is used only for detecting reverse strand constructs."* No HGVS parser in the source.

**`vcf_vep_output` = no.** VCF is the *input* format. Output is a TSV of sequences (`write_tsv(successes, path = outPath)`, optionally an `.RData`). Paper: *"users provide a VCF containing variants to receive a tab-separated file containing the MPRA sequences for their experiment."* Nothing VEP-consumable is emitted.

**`consequence_annotation` = no.** No molecular-consequence vocabulary. The `notes` and `failed$reason` fields describe *construction* events (context shortened, alleles flipped, restriction site hit), not variant consequence.

### Block D — physical construction

**`primer_design` = no.** `fwprimer` and `revprimer` are required user-supplied strings ("a string containing the forward PCR primer to be used") that are simply pasted onto the construct ends; defaults are the fixed Melnikov-derived `ACTGGCCAG` / `CTCGGCGGCC`. There is no Tm calculation, no primer3, no mutagenic-primer generation.

**`codon_optimization` = no.** No codon table, no translation, no optimization objective.

**`synthesis_constraints` = partial.** Real but narrow, and cloning-oriented rather than synthesis-oriented. What it does check on the assembled construct: `countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)` counts aberrant KpnI/XbaI/SfiI sites in each full sequence, including sites created *at the junctions between elements* (`'Failed - aberrant site generated at interface between sequence elements'`). Ambiguous IUPAC N is supported in the enzyme patterns ("The three enzyme arguments may contain ambiguous nucleotides"). SfiI is checked but never inserted, deliberately: *"The sequence for enzyme3 does not show up in the output sequences, however it is necessary to check for it's presence in the output sequences as it is used when preparing the plasmid library."* On failure it retries with resampled barcodes, and with `alter_aberrant = TRUE` it will *"randomly alter aberrant digestion sites across barcodes"* (`randomly_fix`), recording the edits in `site_fix_info`. `max_construct_size` trims flanking context symmetrically until the oligo fits. What it does **not** check: GC content, homopolymer runs, repeats, secondary structure, or synthesis complexity of the construct as a whole (homopolymer/miR/polyA screening is applied to the **barcode** only). `alter_aberrant`, `extra_elements` and `max_construct_size` are all flagged "under development" in `man/processVCF.Rd`.

### Block E — engineering

**`interface`:** R package API (single main function `processVCF()`) + an R/Shiny web GUI. **No CLI.** The GUI additionally caps output at 20,000 barcoded sequences ("To avoid overloading our Shiny server ... at most 20000") and locks barcodes to 12 bp; the app itself directs serious users to the package.

**`license`:** `GPL-2`, declared only in `DESCRIPTION`. There is no `LICENSE`/`COPYING` file in the repository and the GitHub API reports `license: null` for both `mpradesigntools` and `designMPRA`.

**`maintained`:** `mpradesigntools` last commit **2020-12-17** ("replace data_frame with tibble"); `designMPRA` last commit **2017-09-26** ("final chad changes"). Both repos are public, unarchived, 0 open issues, 3 and 1 stars respectively.

---

## 4. The tool's own documented examples / vignettes

There is **no `vignettes/` directory, no `tests/`, no `inst/extdata`, and no example VCF shipped with the package**. What exists:

1. **README usage example** (the only complete invocation anywhere):
   ```r
   processVCF(vcf = '/path/to/the.vcf', nper = 14,
              upstreamContextRange = 55, downstreamContextRange = 55,
              outPath = '/path/to/the/output.tsv',
              fwprimer = 'ACTGGCCGCTTCACTG', revprimer = 'AGATCGGAAGAGCGTCG',
              alter_aberrant = TRUE, extra_elements = FALSE,
              max_construct_size = 170, barcode_set = 'barcodes14-1',
              ensure_all_4_nuc = TRUE)
   ```
   A 170 bp oligo, 14 barcodes per allele, 55 bp of hg38 context each side, error-correcting 14-mer barcodes. **This is the single reproduction target.**

2. **Shiny "Create sequences" defaults, stated to mirror a real published library:** *"The default values are comparable to those used in Ulirsch et al. 2016"* — 14 barcodes per SNP, 75 bp genomic context each side, KpnI/XbaI/SfiI, primers `ACTGGCCAG` / `CTCGGCGGCC`.

3. **Melnikov et al. 2012 construct layout** as the fixed default architecture (Supplement S1.3, plus the `MPRAseqDiagram.png` schematic).

4. **`generate12mers.R`** (Supplement S5, full code reproduced in the paper's supplement): generate all 12-mers → require all four nucleotides → drop runs ≥4 → drop KpnI/XbaI sites → drop 12-mers starting `TCT` → drop human miRBase seed matches → 1,140,292 inert barcodes. A reproducible barcode-design use case in its own right.

5. **Power-analysis vignette (not library design):** Supplement S2/S3 re-analyse Tewhey et al. 2016 and Ulirsch et al. 2016 barcode counts to estimate activity SD (mean ≈ 0.926) and justify the t-test; S4 gives the `pwr.t.test` call driving the app.

---

## 5. Notable capabilities not in the row list

- **Statistical power analysis** — the paper's headline feature and entirely orthogonal to sequence design. Interactive `pwr.t.test` power curves over transcriptional shift, parameterized by barcodes/allele, transfection replicates, number of variants (Bonferroni), alpha, and assumed activity SD, and *calibrated on real MPRA data* (Tewhey 2016, Ulirsch 2016). Validated against Monte Carlo simulation (max abs. difference 0.00341).
- **Error-correcting barcode sets.** 27 `freebarcodes` sets indexed by length (3–17 bp) and number of correctable indel errors (1–2), selectable via `barcode_set`.
- **Automatic repair of aberrant restriction sites** with per-base provenance (`randomly_fix` / `site_fix_info`), plus retry-by-barcode-resampling.
- **Structured failure reporting.** A second `failed` tibble with an explicit machine-readable `reason` per rejected variant, rather than silently dropping.
- **VCF normalization utility.** `spread_and_fix_indels()` splits multi-allelic rows and rewrites indels to the dash convention, emitting `*_fixed.vcf`.
- **Automatic strand handling.** `flip_RV = TRUE` reverse-complements alleles carrying dbSNP's `RV` tag; `MPRAREV` flips the construct orientation and swaps the up/downstream context semantics.
- **Runtime/feasibility estimation.** The app pre-computes and displays expected job time and warns when the design exceeds the available barcode pool.
- **Named downstream analysis companion**, `malacoda` (Bayesian MPRA analysis, PLOS Comput Biol 2020) — a design→analysis pairing PoolParty does not attempt.

---

## 6. Availability today (checked 2026-08-10)

| Component | Status |
|---|---|
| Web app https://andrewghazi.shinyapps.io/designmpra/ | **DEAD** — HTTP 404, page body renders shinyapps.io's `404 / Not Found`. The URL is the primary availability link in the *Bioinformatics* paper. |
| `github.com/andrewGhazi/mpradesigntools` | Alive, public, unarchived. Last commit 2020-12-17. v0.2.0. 3 stars, 1 fork, 0 open issues. |
| `github.com/andrewGhazi/designMPRA` | Alive, public, unarchived. Last commit 2017-09-26. 1 star. |
| CRAN / Bioconductor | **Not published on either.** GitHub-only, installed via `devtools::install_github()`. |
| Installability risk | Non-trivial. README still instructs `source("https://bioconductor.org/biocLite.R"); biocLite(...)` — `biocLite` was removed from Bioconductor in 2018 (superseded by `BiocManager`), so the documented install path is broken. The code also uses `tidyr::unnest_legacy()`, `purrrlyr`, and `readr::write_tsv(path=)` (deprecated in favour of `file=`), and requires the ~700 MB `BSgenome.Hsapiens.UCSC.hg38`. Likely still runnable with effort on a pinned toolchain, but not verified — no install was attempted this round per instructions. |

**Bottom line:** the R package is fetchable and the source is intact, but the interactive web tool that the paper advertises is gone, the package was never distributed through CRAN/Bioconductor, its documented installation instructions are stale, and it has not been touched in ~5.5 years.

---

## 7. Confidence

Highest confidence on Blocks B/C/D and on `dag_chaining`, `combinatorial_motif_place`, `assay_dms`, `barcode_generation` — read directly from the full 1459-line source, `NAMESPACE`, and man pages, which are short enough to be exhaustive.

Least confident:
- **`library_as_object` (partial)** — judgement call. One call yields the whole library so the user writes no loop, but the product is a bare tibble with no class or methods. A referee could argue "yes" (it does make the whole library at once) or "no" (there is no library abstraction).
- **`per_sequence_provenance` (partial)** and **`automatic_naming` (partial)** — the metadata is real and structured (`site_fix_info`, `notes`), but fixed-schema and single-purpose; naming is by identifier columns rather than a composed name string.
- **`synthesis_constraints` (partial)** — restriction-site and max-length enforcement on the assembled construct is genuine, but there is no GC/homopolymer/repeat/structure checking of the oligo, and three of the relevant arguments are marked "under development".
- **Installability** — asserted from stale README instructions and deprecated dependency calls, not from an actual install (installation was out of scope).
