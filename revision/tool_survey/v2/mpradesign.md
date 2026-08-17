# v2 capability record — MPRA Design Tools

**Slug:** `mpradesign`
**Full name:** MPRA Design Tools (Shiny app `designMPRA` + companion R package `mpradesigntools` v0.2.0)
**Citation key:** `Ghazi2018aa` — Ghazi AR *et al.*, *Bioinformatics* 34(15):2682–2683 (2018). doi:10.1093/bioinformatics/bty150
**Tier:** 2
**Basis:** re-scored against ROWS_v2 from `final/mpradesign.md` (FINAL, adversarially reviewed, re-verified 2026-08-10). No new extraction was performed; no cell required a web re-check — every v2 row, including all seven new ones, is answerable from the final memo's §3, §4, §5 and §7.
**Date:** 2026-08-10

---

## Block A — library specification

| Key | Value |
|---|---|
| `library_first_class_object` | **partial** |
| `composable_operations` | **no** |
| `lazy_generation` | **no** |
| `library_algebra` | **no** |
| `exhaustive_single_scans` | **no** |
| `sampled_random_mutagenesis` | **no** |
| `higher_order_combinatorial` | **no** |
| `heterogeneous_components_one_library` | **no** |
| `combinatorial_motif_place` | **no** |
| `barcode_generation` | **yes** |
| `per_sequence_provenance` | **partial** |
| `automatic_naming` | **partial** |
| `design_visualization` | **no** |

**`library_first_class_object` = partial.**
`processVCF()` returns an in-memory `list(result = <tibble>, failed = <tibble>)`; `outPath` is optional, so the tool is genuinely *not* a file-writer-only — the user can hold and inspect the whole library in R and manipulate it with any dplyr verb. That earns real credit. It stops short of `yes` because there is no library abstraction: no S3/S4 class, no `methods` import, a hard-coded terminal `select(ID, type, allele, snpIndex, totIndex, barcode, sequence, site_fix_info, notes)`, and — decisively for "pass onward" — **no operation anywhere consumes a library and returns a library**. `NAMESPACE` exports exactly two functions (`processVCF`, `spread_and_fix_indels`) and zero datasets; the 28 barcode sets arrive via `LazyData: true`. `spread_and_fix_indels()` operates on a VCF *file*, upstream of the library, not on the library.
*Source:* final memo §3 `library_as_object`, §2; `NAMESPACE`; `DESCRIPTION`; `R/processVCFfast.R` L1252, L1290.

**`composable_operations` = no.** *(rename of `dag_chaining`; value and evidence carried unchanged)*
Two exported functions total; the design is one fixed pipeline hard-coded inside `processSnp()` (~890 lines, L210–L1098) with no composition, nesting or reordering. Strongest independent demonstration: the only fork of the package (`mmandelia/mpradesigntools`, 2 commits, 2023-03-24) had to **edit `R/processVCFfast.R` in place** to swap the KpnI/XbaI/SfiI trio for AgeI/SbfI/I-SceI and to drop the barcode. Changing construct architecture requires forking the source.
*Source:* final memo §3 `dag_chaining`; `NAMESPACE`; GitHub API `/repos/andrewGhazi/mpradesigntools/forks`.

**`lazy_generation` = no.** *(rename of `lazy_evaluation`; value and evidence carried unchanged)*
Everything is materialized eagerly. `mers = get(barcode_set)` loads an entire set (`data/twelvemers.rda`, 3.92 MB / 1,140,292 entries) and filters it up front; the whole VCF is read to a tibble; constructs are `Reduce('rbind')`-ed. `processVCF` hard-errors *before* generating anything if `nrow(vcf) * 2 * nper > length(mers)`. No streaming, generators or chunking; README quantifies the eager cost as ≈ `5 + .01 * barcodes-per-allele * SNPs * 2` seconds.
*Source:* final memo §3 `lazy_evaluation`; `R/processVCFfast.R` L1180–L1220; `README.md`.

**`library_algebra` = no.** *(new half of the `library_as_object` split)*
No stack/concat/sample/repeat operation exists over libraries. The output is a bare tibble pair; merging the products of two `processVCF()` runs is an external `rbind`/`bind_rows` the user writes themselves, and doing so would break the tool's own guarantee of library-wide barcode disjointness (pools are partitioned *within* a single call: `split(shuffled_mers, ceiling(1:length(shuffled_mers) %% nrow(vcf)))`). Per the v2 rule — external script ⇒ `no`.
*Source:* final memo §3 `library_as_object`, §5 "Library-wide barcode uniqueness"; `NAMESPACE`; `R/processVCFfast.R` L1180–L1220.

**`exhaustive_single_scans` = no.** *(from the `mixed_mutagenesis_one_pool` split)*
There are **no mutagenesis schemes at all** — not "only one kind". The sole per-variant variation is `type = rep(c('ref','alt'), each = nper)`: the alleles come from the user's VCF, they are not enumerated by the tool. No saturation scan, no positional walk, no deletion scan. Positive-absence verified over the full 1459-line source and the 28-page PDF.
*Source:* final memo §3 `mixed_mutagenesis_one_pool`, §6; `R/processVCFfast.R` L5–L21, L36–L71, L210–L1098.

**`sampled_random_mutagenesis` = no.** *(from the split)*
No mutation rate, no sampled variant mode, no random mutagenesis parameter. Randomness in the tool serves two unrelated purposes: barcode assignment (`sample(..., replace = FALSE)` within disjoint per-variant pools) and constraint repair (`randomly_fix()` samples `nrow(res_df)/2` different single-base alterations of an offending restriction site, matched across ref and alt). Neither samples variants into the design; the second is site repair whose edits are logged as damage in `site_fix_info`.
*Source:* final memo §3 `mixed_mutagenesis_one_pool`, `synthesis_constraints`, §5; `R/processVCFfast.R` L502–L540, L723–L760, L926–L960, L1180–L1220.

**`higher_order_combinatorial` = no.** *(from the split)*
Nothing pairs or combines mutations within one sequence. Each output construct carries exactly one allele of one VCF record; a multi-allelic site is expanded to *separate* rows by `spreadAllelesAcrossRows()`, never combined. `grep -E 'combinator|permut'` returns 0 hits across source, README, app and full PDF.
*Source:* final memo §3 `mixed_mutagenesis_one_pool`, `combinatorial_motif_place`; `R/processVCFfast.R` L5–L21.

**`heterogeneous_components_one_library` = no.** *(from the split; see confidence flag)*
Giving the tool its due: one pool **can** contain SNVs, insertions, deletions and multi-allelic expansions simultaneously, each routed through its own ~200-line construct-generation branch (`generateInsConstruct`, `generateDelConstruct`), and every alt construct is automatically accompanied by a matched ref arm — the closest thing to a WT control. But these are mixed variant **classes** flowing through one scheme, not structurally different component **types**: every construct has the identical hard-coded layout `paste0(fwprimer, constrseq, enzyme1, enzyme2, barcodes, revprimer)`, and the user cannot specify a mix (no exhaustive-singles-plus-sampled-higher-order-plus-controls composition exists to specify). The mix is a consequence of what is in the VCF, not a design the user writes.
*Source:* final memo §3 `mixed_mutagenesis_one_pool`, §2 "The one design operation"; `R/processVCFfast.R` L5–L21, L36–L71, L210–L1098.

**`combinatorial_motif_place` = no.**
Positive-absence confirmed: `grep -E 'motif|spacer|permut|combinator'` returns **0 hits** across `README.md`, `ui.R`, `server.R`, `R/processVCFfast.R`, `man/processVCF.Rd` and the full 28-page paper + supplement. No concept of an element, position, orientation or placement. The only positional control is flanking-context width (`upstreamContextRange`/`downstreamContextRange`); the only orientation control is a per-variant strand flag (`MPRAREV`/`RV`).
*Source:* final memo §3 `combinatorial_motif_place` (full-source and full-PDF grep).

**`barcode_generation` = yes.** *(genuine strength — do not understate)*
(a) Original set: all 12-mers screened per Melnikov *et al.* — all four nucleotides present, no homopolymer run > 3, no KpnI/XbaI sites, no leading `TCT`, no human miR seed matches → **1,140,292 inert twelvemers** (Supplement S5, reproducible code `generate12mers.R`). (b) v0.2.0 adds **27** error-correcting `freebarcodes` sets (Hawkins *et al.* PNAS 2018), `barcodes3-1` … `barcodes17-2`, spanning **3–17 bp** with 1–2 correctable indel errors. (c) Runtime re-filtering: `ensure_all_4_nuc`, `filterPatterns` (default `AATAAA` polyA), and removal of all three enzyme sites plus reverse complements and reverses. (d) **Library-wide uniqueness by construction**: disjoint per-variant barcode pools + `sample(replace = FALSE)` make the barcode→variant map injective with no post-hoc dedup. (e) Documented quality audit of the borrowed sets in `R/import_freebarcodes.R`. Custom barcode vectors accepted. *Real caveat:* no GC-content constraint anywhere.
*Source:* final memo §3 `barcode_generation`, §5; Supplement S5 (merged PDF pp. 25–26); `man/twelvemers.Rd`; `R/processVCFfast.R` L394, L1180–L1220; `R/import_freebarcodes.R`; repo `data/` tree.

**`per_sequence_provenance` = partial.**
Each row carries structured build metadata beyond the mutation name: `type` (ref/alt), `allele`, `barcode`, `snpIndex`, `totIndex`, `site_fix_info` (a nested data frame recording exactly which bases were changed to destroy an aberrant restriction site) and `notes`. Three limits keep it below `yes`: (i) fixed flat schema for one design operation — there is no general "how it was built" record because there is only one way anything is built; (ii) the terminal `select()` **drops `CHROM`/`POS`**, so the delivered library carries no genomic coordinates and traceability depends on the VCF `ID` column being populated; (iii) `notes` is assigned once per variant (`res$notes = notes`, L987) and broadcast to all `2*nper` rows — variant-level, not per-sequence. Only `site_fix_info` genuinely varies row to row.
*Source:* final memo §3 `per_sequence_provenance`; `R/processVCFfast.R` L236, L287, L323, L986–L989, L1252, L1290.

**`automatic_naming` = partial.**
Docs claim `processVCF` yields *"labeled MPRA sequences"*, and identifier columns (`ID`, `type`, `snpIndex`, `totIndex`) are emitted automatically. But labelling is a set of **columns**, not a composed name string — no name-composition code exists anywhere in `processVCFfast.R` (no `paste0` producing an identifier) — and `ID` is carried verbatim from the input VCF, so a VCF without rsIDs yields rows identified only by integer indices. FASTA export requires the user to build names themselves. `partial` is the generous reading; a stricter matrix would say `no`.
*Source:* final memo §3 `automatic_naming`; `man/processVCF.Rd`, `man/processSnp.Rd`.

**`design_visualization` = no.**
`designMPRA/server.R` has eight `output$` assignments (six `render*`); the only plot is `powerPlot`, a `pwr.t.test` power curve — a statistics plot, not a view of the library. `testHead` is a leftover `head(mtcars)` stub. The construct-layout picture is a **static PNG** identical regardless of design, and `www/*` is gitignored so it is absent from the distributed source. The R package has no plotting function at all (no `graphics`/`ggplot2` in `NAMESPACE`); the only visual feedback is a `print()`ed construct-size distribution table.
*Source:* final memo §3 `design_visualization`; `designMPRA/server.R`, `ui.R` L144, `.gitignore`; `NAMESPACE`.

---

## Block B — assay coverage

| Key | Value |
|---|---|
| `assay_dms` | **no** |
| `assay_mpra` | **yes** |
| `assay_insilico` | **no** |

**`assay_dms` = no.** No codon, amino-acid, ORF or reading-frame concept. `Biostrings` imports are strictly nucleotide-level (`DNAString`, `complement`, `replaceLetterAt`, `reverseComplement`, `subseq`, …) — nothing translational. `grep 'codon|amino|ORF|frame'` over full source and 28-page PDF: 0 hits. *Source:* final memo §3; `NAMESPACE`.

**`assay_mpra` = yes.** The tool's entire purpose — title "Design tools for MPRA experiments"; `DESCRIPTION`: *"Design MPRA sequences from input VCF's"*. Scope is narrow *within* MPRA: SNP/indel allelic-contrast libraries in the Melnikov *et al.* 2012 layout; the README's Planned Features leaves *"bed file to Sharpr-MPRA library oligo design"* (tiling MPRA) **unchecked**. *Source:* final memo §3; paper title/abstract; `DESCRIPTION`; `README.md`.

**`assay_insilico` = no.** `DESCRIPTION` Imports are dplyr/purrr/readr/tidyr/tibble/stringr/purrrlyr/magrittr/Biostrings/BSgenome.Hsapiens.UCSC.hg38 only — no model, scorer or ML dependency. "Model" in the paper refers solely to the `pwr.t.test` power model. *Source:* final memo §3; `DESCRIPTION`; `NAMESPACE`.

---

## Block C — genomics integration

| Key | Value |
|---|---|
| `genome_coordinates` | **yes** |
| `transcript_models` | **no** |
| `exon_intron_split_codons` | **no** |
| `vcf_vep_output` | **no** |
| `consequence_annotation` | **no** |

**`genome_coordinates` = yes.** VCF `CHROM`/`POS` slice the reference directly: `subseq(genome[[paste0('chr', snp$CHROM)]], ...)` with `genome = BSgenome.Hsapiens.UCSC.hg38` assigned inside `processSnp` (L227, L766) and imported wholesale in `NAMESPACE`. The paper frames this as its advantage over MPRAnator: *"Our tool acquires genomic context from the hg38 reference genome, rather than requiring input by the user."* *Constraint:* hg38 human only; "mm10 genomic context" is an unchecked Planned Feature; and note the delivered `result` table itself drops `CHROM`/`POS`. *Source:* final memo §3; `R/processVCFfast.R` L227, L337, L555, L766, L777.

**`transcript_models` = no.** No GTF/GFF parsing, no annotation package in Imports. Strandedness is **user-supplied**, not derived: *"the MPRAREV tag will need to be added by the user (where appropriate) because the VCF's do not always specify which strand the relevant gene is on"*; in code `reverseGene = grepl('MPRAREV', INFO)`. *Source:* final memo §3; `README.md`; `R/processVCFfast.R` L1217.

**`exon_intron_split_codons` = no.** Context extraction is a flat genomic window `subseq(genome[[chr]], POS - upstreamContextRange, POS + downstreamContextRange)` — no splice-junction and no reading-frame awareness. Follows from `transcript_models = no` and `assay_dms = no`. *Source:* final memo §3; `R/processVCFfast.R` L337, L555, L777.

**`vcf_vep_output` = no.** VCF is the **input** format. Output is `write_tsv()` (+ optional `.RData`); `readr::write_tsv` is the only writer in `NAMESPACE`, and there is no FASTA or GenBank writer. The only VCF ever *written* is `spread_and_fix_indels()`'s normalized `*_fixed.vcf` copy of the **input** — a prestep, not a result. *Source:* final memo §3; `R/processVCFfast.R` L1418–L1450.

**`consequence_annotation` = no.** No molecular-consequence vocabulary. The only annotation-like strings are **construction events**: `notes` (e.g. *"Shortened context by … bp on each side"*, *"alleles … flipped … because of the presence of the RV tag"*) and five `failed$reason` values, all describing restriction-site handling. *Source:* final memo §3; `R/processVCFfast.R` L136, L177, L287, L323, L506, L727, L930.

---

## Block D — adjacent / complementary

| Key | Value |
|---|---|
| `primer_design` | **no** |
| `codon_optimization` | **no** |
| `synthesis_constraints` | **partial** |
| `degenerate_iupac_codons` | **no** |
| `negative_control_generation` | **no** |
| `ml_model_in_loop` | **no** |
| `readout_analysis` | **no** |

**`primer_design` = no.** `fwprimer`/`revprimer` are **required user-supplied strings** pasted onto the construct ends (defaults `ACTGGCCAG` / `CTCGGCGGCC` per Melnikov *et al.*). No Tm calculation, no primer3, no mutagenic-primer generation. *Source:* final memo §3; Supplement S1.3 (merged PDF p.5 = supplement p.3); `man/processVCF.Rd`.

**`codon_optimization` = no.** No codon table, no translation, no optimization objective, and nothing in the dependency graph capable of providing one. *Source:* final memo §3; `DESCRIPTION`; `NAMESPACE`; full-source grep.

**`synthesis_constraints` = partial.** Real and more sophisticated than a one-liner suggests, but cloning-oriented rather than synthesis-oriented. *Does:* `countDigSites()` counts aberrant KpnI/XbaI/SfiI sites in each **fully assembled** construct, including sites created **at element junctions**; enzyme patterns are IUPAC-aware (`fixed = FALSE`); SfiI is checked but never inserted (needed for plasmid-library prep); failure triggers barcode resampling under a **bounded 40-attempt budget** per variant; `alter_aberrant = TRUE` invokes `randomly_fix()` (different single-base repair per barcode replicate, **matched** between ref and alt so the allelic contrast is not confounded), `nonrandomly_replace_N()` protects N positions of ambiguous sites, `check_cross_center()` hard-fails variants inside the aberrant site, and per-base edits are recorded in `site_fix_info`; `max_construct_size` trims context symmetrically until the oligo fits. *Does not:* GC content, homopolymer runs, repeats, secondary structure, or synthesis-complexity scoring **of the construct** (homopolymer/miR/polyA screening applies to the barcode only). `alter_aberrant`, `extra_elements`, `max_construct_size` are all tagged "under development". *Source:* final memo §3, §5; `R/processVCFfast.R` L22–L35, L72–L184, L502–L540, L723–L760, L926–L960, L1451–L1459; `man/processVCF.Rd`.

**`degenerate_iupac_codons` = no.** *(new row)* IUPAC ambiguity appears in exactly one place and it is not codon design: the three **enzyme recognition patterns** may contain ambiguous nucleotides (`enzyme3 = "GGCCNNNNNGGCC"`, SfiI, matched with `fixed = FALSE`), and `nonrandomly_replace_N()` exists only to protect the N positions of such a site during repair. There is no degenerate-codon specification, no expansion of a degenerate string into sequences, and no compression of sequences into a degenerate string — and no codon concept at all (`assay_dms = no`; `grep 'codon|amino|ORF|frame'` = 0 hits over full source and full PDF). Output sequences are fully specified A/C/G/T constructs.
*Source:* final memo §3 `synthesis_constraints`, `assay_dms`, §2 (function signature); `R/processVCFfast.R` L124–L184; `man/processVCF.Rd`.

**`negative_control_generation` = no.** *(new row)* No scramble, shuffle, dinucleotide-preserving permutation, or reverse/complement control generator — `grep 'motif|spacer|permut|combinator'` returns 0 hits and nothing in the source, README, app, paper or supplement offers a control-sequence feature. What exists is adjacent but different: every alt construct is automatically paired with a matched **ref (reference-allele) arm**, `type = rep(c('ref','alt'), each = nper)`, which functions as the experiment's internal comparator; and `flip_RV`/`MPRAREV` reverse-complement constructs for genes on the minus strand — an orientation correction, not a control. Neither is a user-requestable negative-control class, and the user cannot add controls to the pool through the tool.
*Source:* final memo §3 `mixed_mutagenesis_one_pool`, `combinatorial_motif_place`, §5 "Automatic strand handling"; `R/processVCFfast.R` L5–L21, L36–L71, L1217.

**`ml_model_in_loop` = no.** *(new row)* Nothing in the dependency graph can score a sequence: `DESCRIPTION` Imports are dplyr/purrr/readr/tidyr/tibble/stringr/purrrlyr/magrittr/Biostrings/BSgenome.Hsapiens.UCSC.hg38, plus `pwr` in the app. No predictive model, no scoring function, no selection or optimization step over generated sequences — construct generation is deterministic given the VCF, barcode set and enzyme parameters. The only "model" in the paper is `pwr.t.test`, a statistical power model applied *before* design, not a sequence predictor in the loop. *Source:* final memo §3 `assay_insilico`, §5; `DESCRIPTION`; `NAMESPACE`.

**`readout_analysis` = no.** *(new row)* The tool ends at the TSV of oligos; it never ingests sequencing counts. Two things could be mistaken for a yes and are not: (i) the **power analysis** — the paper's headline feature — is prospective, taking user-entered barcodes/allele, replicates, variant count, alpha and assumed activity SD, and it was *calibrated* on published Tewhey 2016 / Ulirsch 2016 barcode counts in the supplement (combined mean activity SD ≈ .95), but no count-data ingestion is exposed to the user; (ii) `malacoda` (Bayesian MPRA analysis, PLOS Comput Biol 2020) is a **separate package by the same author**, cited as a companion, not a feature of `mpradesigntools`/`designMPRA`. Scored on this tool's own code, `no`.
*Source:* final memo §3 `assay_insilico`, §4 item 5, §5 (power analysis; `malacoda`); `designMPRA/server.R` (`powerPlot`); Supplement S2–S4.

---

## Block E — engineering and availability

*(Per the v2 instrument, value is `yes` = row addressed; the real answer is in the evidence.)*

| Key | Value | Short answer |
|---|---|---|
| `interface` | **yes** | R package API (one main function) + R/Shiny web GUI; **no CLI** |
| `license` | **yes** | GPL-2, declared in `DESCRIPTION` only |
| `installable_today` | **yes** | `devtools::install_github()` only — never on CRAN/Bioconductor/r-universe; documented install path broken since 2018 |
| `last_activity` | **yes** | package 2020-12-17; app 2017-09-26; **0 tags, 0 releases** |
| `documented_examples` | **yes** | 0 vignettes, 0 tests, 0 example data; 1 README invocation + 1 supplement script + power-analysis supplement |

**`interface`.** R package API (single main function `processVCF()`) **plus** an R/Shiny web GUI. **No CLI** — verified absence of `inst/`, `exec/` and any `Rscript` entry point in the recursive git tree. GUI limits: ≤ 20,000 barcoded sequences ("to avoid overloading our Shiny server"), barcodes locked to 12 bp, and the app directs serious users to the package. Two verifiable caveats: (1) the **hosted** instance https://andrewghazi.shinyapps.io/designmpra/ returns **HTTP 404** (re-verified 2026-08-10; case variant and account root also 404) — the source is still runnable locally via `shiny::runApp()`; (2) a local run is degraded and stale — `www/*` is gitignored so all three referenced images are absent from the distributed source, and `server.R` sources an **18,609-byte 2017 snapshot** of the engine (vs the package's 62,351 B) whose `processVCF` has no `barcode_set`, `alter_aberrant`, `max_construct_size`, `flip_RV`, `filterPatterns` or `ensure_all_4_nuc`, and only a symmetric `seqwidth`. The GUI cannot reach any v0.2.0 feature. *(The earlier "GUI is non-functional" claim was overstated and is withdrawn.)*
*Source:* final memo §3 `interface`, §2, §7; `designMPRA/ui.R` L108, L126, L144; `server.R` L17, L114; `.gitignore`; recursive git trees; live HTTP checks.

**`license`.** `DESCRIPTION` declares `License: GPL-2`. Caveat: declared **only** there — the GitHub REST API returns `license: null` for **both** repos, and neither recursive file tree contains a `LICENSE`/`COPYING` blob. *Source:* final memo §3 `license`; `DESCRIPTION`; GitHub API.

**`installable_today`.** Not on CRAN (**HTTP 404, never published**), not on Bioconductor (URL 200 but serves the "Removed Packages" fallback — never indexed), 0 results on r-universe; no conda recipe, no Docker image, no pip equivalent (R package). The only path is `devtools::install_github('andrewGhazi/mpradesigntools')`. Risk is **non-trivial and unverified** (no install attempted, per survey instructions): the README's documented install still says `source("https://bioconductor.org/biocLite.R"); biocLite(...)` and `biocLite` was removed in 2018; the code uses `tibble::data_frame` (deprecated), `tidyr::unnest_legacy()`, `purrrlyr`, and `readr::write_tsv(path=)` (deprecated); and it requires the ~700 MB `BSgenome.Hsapiens.UCSC.hg38`. Likely still runnable on a pinned 2020-era toolchain. The web app the paper advertises as the primary availability link is **dead**.
*Source:* final memo §6, §7.

**`last_activity`.** `mpradesigntools` last commit **2020-12-17T19:27:54Z** ("replace data_frame with tibble") — 5 yr 8 mo dormant as of 2026-08-10; **0 tags, 0 releases**, 1 branch, 3 stars, 1 fork, 0 open issues. `designMPRA` last commit **2017-09-26** — 8 yr 10 mo dormant, 1 star. Both public and unarchived. The only fork (`mmandelia/mpradesigntools`, 2 commits, 2023-03-24) is a private hack, not a successor; `goldenac/MPRADesignGenerator` is confirmed independent (`fork: false`). No maintained descendant. *(This row supersedes v1 `maintained` = no.)*
*Source:* final memo §3 `maintained`, §7; GitHub API `/repos`, `/tags`, `/releases`, `/forks`.

**`documented_examples`.** Verified against the full recursive git tree: **no `vignettes/`, no `tests/`, no `inst/extdata`, no example VCF**. What exists: (1) **one** complete README invocation — the single reproduction target (170 bp oligo, `nper = 14`, 55 bp context each side, `barcode_set = 'barcodes14-1'`); (2) Shiny "Create sequences" defaults stated to mirror Ulirsch *et al.* 2016 (14 barcodes/SNP, 75 bp context, KpnI/XbaI/SfiI); (3) `generate12mers.R`, a reproducible barcode-design script, in Supplement S5 (merged PDF pp. 25–26) and in `designMPRA/scripts/`; (4) the power-analysis worked examples in Supplement S2–S4 (re-analysis of Tewhey 2016 / Ulirsch 2016) — documentation of the power feature, not of library design. Count as **1 shipped runnable design example + 1 shipped auxiliary script; 0 vignettes, 0 tests**.
*Source:* final memo §4, §7; `README.md`; recursive git tree; Supplement S2–S5.

---

## Changes vs. the v1 record

| v2 key | v1 key | v1 → v2 | Why |
|---|---|---|---|
| `composable_operations` | `dag_chaining` | no → no | pure rename; new wording ("compose/chain/nest vs fixed pipeline") if anything strengthens the `no` — one hard-coded pipeline. |
| `lazy_generation` | `lazy_evaluation` | no → no | pure rename; eager throughout. |
| `library_first_class_object` | `library_as_object` (half) | partial → **partial** | returns an in-memory tibble pair, not files only — but no class, no methods, nothing to pass it onward to. |
| `library_algebra` | `library_as_object` (half) | partial → **no** | no stack/concat/sample over libraries; combining runs is an external `bind_rows` and would break barcode disjointness. |
| `exhaustive_single_scans` | `mixed_mutagenesis_one_pool` | no → no | no mutagenesis scheme of any kind; alleles come from the VCF. |
| `sampled_random_mutagenesis` | `mixed_mutagenesis_one_pool` | no → no | randomness is barcode assignment and site repair, not variant sampling. |
| `higher_order_combinatorial` | `mixed_mutagenesis_one_pool` | no → no | one allele per construct; multi-allelic sites are split, never combined. |
| `heterogeneous_components_one_library` | `mixed_mutagenesis_one_pool` | no → no | mixed variant *classes* (SNV/ins/del/multi-allelic + paired ref arms) coexist, but one fixed construct layout and no user-specified mix. |
| `degenerate_iupac_codons` | — | new → **no** | IUPAC appears only in enzyme recognition patterns; no codon concept. |
| `negative_control_generation` | — | new → **no** | no scramble/shuffle/RC generator; the automatic ref arm is a comparator, not a control feature. |
| `ml_model_in_loop` | — | new → **no** | no scorer in the dependency graph; the only "model" is `pwr.t.test`, pre-design. |
| `readout_analysis` | — | new → **no** | tool ends at the oligo TSV; `malacoda` is a separate package. |
| `installable_today` | — | new → yes (github-only, risky) | never on CRAN/Bioconductor/r-universe; documented install broken since 2018. |
| `last_activity` | (`maintained` = no) | new → yes (2020-12-17 / 2017-09-26) | same underlying evidence, restated as a date; 0 tags, 0 releases. |
| `documented_examples` | — | new → yes (1 + 1; 0 vignettes/tests) | verified against the recursive git tree. |
| `hgvs_input` | dropped | — | v1 value was `no` (`readLines` + `read_tsv` on a VCF path; 0 HGVS hits). Not output. |
| `maintained` | not in v2 | — | superseded by `last_activity` + `installable_today`. |

Unchanged and carried verbatim: `barcode_generation` (yes), `per_sequence_provenance` (partial), `automatic_naming` (partial), `design_visualization` (no), `combinatorial_motif_place` (no), all of Block B, all of Block C, `primer_design`/`codon_optimization` (no), `synthesis_constraints` (partial), `interface` (yes), `license` (yes).

## Confidence flags (for human review)

- `library_first_class_object` = partial — the judgement call. A referee could argue `yes` (a tibble *is* a first-class R object you can hold, inspect and transform with dplyr; `outPath` is optional) or `no` (there is no library abstraction, no class, nothing to pass it to). `partial` is the defensible middle and the generous-to-competitor reading.
- `heterogeneous_components_one_library` = no — the strongest candidate for `partial` in this record. SNVs, insertions, deletions and multi-allelic expansions genuinely do coexist in one pool via distinct code branches, and every alt is paired with a ref arm. Scored `no` because the construct layout is single and fixed and the mix is dictated by the input VCF, not composed by the user. Worth a human decision, since the same reasoning must be applied identically to the other 11 tools.
- `per_sequence_provenance` = partial — real, structured metadata, but fixed-schema, coordinate-free (`CHROM`/`POS` dropped), and `notes` is variant-level. Left generous.
- `automatic_naming` = partial — identifier columns, not a composed name; `ID` verbatim from input. A stricter matrix would say `no`.
- `synthesis_constraints` = partial — well calibrated but sits near the `yes`/`partial` line: the restriction-site machinery is sophisticated (junction checking, IUPAC awareness, matched per-replicate repair, 40-attempt retry) while GC/homopolymer/repeat/structure checks on the construct are absent and three relevant arguments are "under development".
- `negative_control_generation` = no — depends on whether the automatic ref-allele arm counts as a first-class control. I read it as the experiment's comparator arm rather than a control-generation feature; a referee (possibly this tool's author) could push back.
- `readout_analysis` = no — depends on the scoring boundary. If the survey credits a tool for a *named companion* analysis package by the same author, `malacoda` would make this `partial`. Scored strictly on this tool's own code. Needs a consistent rule across all 12 tools.
- `installable_today` = yes-with-caveats — **asserted, not tested**. No install was attempted (out of scope this round); the broken-`biocLite`/deprecated-call/700 MB-BSgenome reasoning is inference from source, not a failed build.
- `documented_examples` = the count depends on convention: 1 if only shipped runnable package examples count, 2 including the supplement's `generate12mers.R`, 4 if Shiny defaults and the power-analysis supplement count. Stated all four so the survey can normalize.
- `sampled_random_mutagenesis` = no — noted only because `randomly_fix()` does sample random single-base alterations; it is constraint repair, not design, so the `no` is firm, but the word "random" in the source could mislead a re-reader.

**Highest confidence (read directly from exhaustive source/API/PDF checks):** `composable_operations`, `lazy_generation`, `library_algebra`, `exhaustive_single_scans`, `higher_order_combinatorial`, `combinatorial_motif_place`, `barcode_generation`, `design_visualization`, all of Block B, all of Block C, `primer_design`, `codon_optimization`, `degenerate_iupac_codons`, `ml_model_in_loop`, `interface`, `license`, `last_activity`.
