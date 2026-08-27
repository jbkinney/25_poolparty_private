# FINAL capability record — MPRA Design Tools

**Slug:** `mpradesign`
**Full name:** MPRA Design Tools (Shiny app `designMPRA` + companion R package `mpradesigntools` v0.2.0)
**Citation key:** `Ghazi2018aa`
**Citation:** Ghazi AR, Chen ES, Henke DM, Madan N, Edelstein LC, Shaw CA. "Design tools for MPRA experiments." *Bioinformatics* 34(15):2682–2683 (2018). doi:10.1093/bioinformatics/bty150
**Tier:** 2
**Status of this document:** FINAL. Supersedes `extractions/mpradesign.md ⟨deleted at cleanup — in commit 35d65d8⟩`. Incorporates every correction from `reviews/mpradesign.md ⟨deleted at cleanup — in commit 35d65d8⟩` plus an independent third re-verification pass (2026-08-10) against GitHub raw, the GitHub REST API, and the 28-page PDF.

---

## 0. Changes relative to the original extraction

| Cell | Was | Now | Why |
|---|---|---|---|
| `maintained` | `yes` | **`no`** | The original value contradicted its own evidence string. Re-verified: last commit 2020-12-17, **0 tags, 0 releases**, 1 branch, 0 open issues, 3 stars, not on CRAN (404), not on Bioconductor, not on r-universe. 5 yr 8 mo dormant. |
| `interface` | `yes` (sub-claim "the GUI is currently non-functional") | **`yes`** (sub-claim replaced) | Value correct. The broad sub-claim was narrowed: the *hosted* instance is 404, while a checkout can launch the UI/power tab but cannot run sequence generation without an absent, gitignored barcode object. The record also notes the missing `www/*` images and the 2017 engine snapshot. |
| `barcode_generation` evidence | "10–17 bp" | **"3–17 bp"** | Shipped sets are `barcodes3-1` … `barcodes17-2`; verified in the repo file tree and the README table. Only `barcodes17-2` exists at 17 bp. |
| `library_as_object` evidence | "NAMESPACE exports two functions and 28 datasets" | **"two exported functions; 28 lazy-loaded datasets"** | Re-read NAMESPACE: exactly `export(processVCF)` + `export(spread_and_fix_indels)`; **zero** dataset exports. Datasets reach the user via `LazyData: true`. |
| `design_visualization` evidence | "exactly six outputs" | **"eight `output$` assignments, six of them `render*`"** | Verified: `powerPlot`, `testHead`, `timeText`, `warnText`, `conditionalDownload`, `failed`, `inputHead`, `downloadSequences`. Verdict unaffected. |
| Nine capabilities the reviewer flagged as missed | — | added | See §5. |

**Direction of error check:** every remaining `partial` (`library_as_object`, `per_sequence_provenance`, `automatic_naming`, `synthesis_constraints`) is calibrated *generously toward the competitor*. That is the correct direction for a referee response.

---

## 1. Sources consulted

| Kind | Reference |
|---|---|
| PDF (paper + full supplement, 28 pp) | `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/../../../lit_review/analyzed/Ghazi2018_MPRADesignTools_all.pdf` |
| Prior analysis | `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/prior_analyses/Ghazi2018_MPRADesignTools_all_analysis.md` |
| Repo (R package) | https://github.com/andrewGhazi/mpradesigntools (re-fetched raw: `DESCRIPTION`, `NAMESPACE`, `R/processVCFfast.R` 1459 lines / 62,351 B, `R/import_freebarcodes.R`, `R/data.R`, `README.md`, all 32 man pages, full recursive git tree) |
| Repo (Shiny app) | https://github.com/andrewGhazi/designMPRA (re-fetched raw: `ui.R`, `server.R`, `.gitignore`, full recursive git tree) |
| GitHub REST API | `/repos/andrewGhazi/{mpradesigntools,designMPRA}` + `/tags`, `/releases`, `/forks` |
| Web service | https://andrewghazi.shinyapps.io/designmpra/ — **HTTP 404** (re-checked 2026-08-14) |
| CRAN | https://cran.r-project.org/package=mpradesigntools — **HTTP 404** |
| Bioconductor | https://bioconductor.org/packages/mpradesigntools/ — HTTP 200 but redirects to the generic **"Removed Packages" page**, which does not name the package; never indexed established separately by its absence from the release package indexes (checked Bioc 3.6–3.14) and from the removed-packages list |
| r-universe | 0 search results |
| Docs mirror | https://rdrr.io/github/andrewGhazi/mpradesigntools/ |

---

## 2. What the tool actually is

Two artifacts, and — a fact worth stating plainly — **two different programs**, not one engine with two front-ends:

1. **`designMPRA`** — an R/Shiny web app (two files at repo root, `ui.R` + `server.R`). Three tabs: "MPRA" (static text + PNG diagrams), "Power" (interactive `pwr.t.test` curve), "Create sequences" (upload VCF → download TSV). Last commit **2017-09-26**. `server.R` does `source('scripts/processVCFfast.R')` — an **18,609-byte 2017 snapshot**, not the package's 62,351-byte `R/processVCFfast.R` — and calls `processVCF(vcf, nper, seqwidth=, fwprimer, revprimer, enzyme1..3, updateProgress)`. That signature has **no** `barcode_set`, `alter_aberrant`, `max_construct_size`, `flip_RV`, `filterPatterns`, or `ensure_all_4_nuc`, and only a **symmetric** `seqwidth` instead of separate up/downstream ranges. The GUI therefore cannot reach any v0.2.0 feature.
2. **`mpradesigntools`** — the companion R package, v0.2.0, last commit **2020-12-17**. Two exported functions: `processVCF()` and `spread_and_fix_indels()`. 28 barcode datasets reach the user through `LazyData: true`, not through `export()`.

The headline scientific contribution is the **power analysis**, not the sequence design. Abstract: *"an online web-tool and R package that allows for interactive MPRA experimental design encompassing both power analysis and design of constructs."* The app says of its own generator: *"We would like to emphasize that this sequence generation functionality of the web tool is mainly for simple demonstrative purposes"* (`ui.R` L108).

### The one design operation

`processVCF()` is a single monolithic function (`man/processVCF.Rd`):

```r
processVCF(vcf, nper, upstreamContextRange, downstreamContextRange,
           fwprimer, revprimer,
           enzyme1 = "GGTACC", enzyme2 = "TCTAGA", enzyme3 = "GGCCNNNNNGGCC",
           filterPatterns = "AATAAA", alter_aberrant = FALSE, extra_elements = FALSE,
           max_construct_size = NULL, barcode_set = "twelvemers",
           ensure_all_4_nuc = TRUE, flip_RV = TRUE, outPath = NULL)
```

Fixed construct layout, hard-coded:

```r
sequence = paste0(fwprimer, constrseq, enzyme1, enzyme2, barcodes, revprimer)
# extra_elements = TRUE:
sequence = paste0(fwprimer, 'TG', constrseq, enzyme1, enzyme2, barcodes, 'GGC', revprimer)
```

Output: `list(result = <tibble>, failed = <tibble>)`, with `result` columns fixed by a terminal
`select(ID, type, allele, snpIndex, totIndex, barcode, sequence, site_fix_info, notes)`
(`R/processVCFfast.R` L1252 and L1290) — **`CHROM` and `POS` are dropped**.

---

## 3. Capability-by-capability (FINAL values)

### Block A — library specification

**`library_as_object` = partial.**
*Evidence:* One `processVCF()` call materializes the whole library, so the user writes no per-variant loop — but the product is a plain tibble pair, `list(result = <tibble>, failed = <tibble>)`, with a hard-coded terminal `select()`. There is no S3/S4 class and no `methods` import. `NAMESPACE` (re-read verbatim) contains exactly **two** export lines — `export(processVCF)`, `export(spread_and_fix_indels)` — and **zero** dataset exports; the 28 barcode datasets are lazy-loaded via `LazyData: true`, not exported. No operation consumes a library and returns a library. The only compositional pairing is the file handoff `spread_and_fix_indels()` → `*_fixed.vcf` → `processVCF()`.
*Source:* `NAMESPACE`; `DESCRIPTION`; `R/processVCFfast.R` L1252, L1290.
*Note:* a stricter matrix could justify `no`; `partial` is deliberately the generous reading.

**`dag_chaining` = no.**
*Evidence:* Two exported functions total; the design is a single fixed pipeline hard-coded inside `processSnp()` (a 784-line function, L210–L993) with no composition, nesting, or reordering. The strongest independent demonstration: the **only fork of the package** (`mmandelia/mpradesigntools`, 2 commits, pushed 2023-03-24) had to **edit `R/processVCFfast.R` in place** to swap the KpnI/XbaI/SfiI enzyme trio for AgeI/SbfI/I-SceI and to drop the barcode from the construct. Changing the construct architecture requires forking the source.
*Source:* `NAMESPACE`; `R/processVCFfast.R`; GitHub API `/repos/andrewGhazi/mpradesigntools/forks`.

**`lazy_evaluation` = no.**
*Evidence:* Everything is materialized eagerly. `mers = get(barcode_set)` loads a whole set into memory (`data/twelvemers.rda` is 3.92 MB / 1,140,292 entries) and filters it up front; the whole VCF is read to a tibble; all constructs are `Reduce('rbind')`-ed. `processVCF` hard-errors *before* generating anything if `nrow(vcf) * 2 * nper > length(mers)`: `'Your design requires more barcodes than is possible with the selected barcode_set. Try a bigger set.'` (That capacity test runs against the unfiltered pool, *before* `ensure_all_4_nuc` and `filterPatterns` shrink it, so a design that clears the advertised set size can still run short of usable barcodes.) No streaming, no generators, no chunking. README quantifies the eager cost: runtime ≈ `5 + .01 * barcodes-per-allele * SNPs * 2` seconds.
*Source:* `R/processVCFfast.R` L1121–L1145 (VCF read), L1159 (`get(barcode_set)`), L1167–L1173 (capacity error), L1175–L1202 (filters), L1243–L1290 (result assembly / `Reduce('rbind')`); `README.md`; `data/` tree.

**`mixed_mutagenesis_one_pool` = no.**
*Evidence:* **There are no mutagenesis schemes at all** — not "only one kind of sequence". The only per-variant variation is `type = rep(c('ref','alt'), each = nper)`. No exhaustive singles, no double mutants, no random sampling, no scanning, no saturation. *Pre-empting the obvious rebuttal:* one pool **can** contain SNVs, insertions, deletions and multi-allelic expansions simultaneously, each with its own construct-generation branch (three ~200-line blocks inside `processSnp`, plus `generateInsConstruct`/`generateDelConstruct`), and `spreadAllelesAcrossRows()` duplicates the ref arm per ALT. Even those classes are narrowly defined: only canonical single-base SNVs (`REF` and `ALT` both in A/C/G/T) and dash-coded insertions/deletions are recognised, and anything else — MNVs, complex alleles — is rejected with `'Failed - Not identifiable as SNV, insertion, or deletion.'` That is mixed variant **classes**, not mixed mutagenesis **schemes**, so the value stands at `no`.
*Source:* `R/processVCFfast.R` L5–L21 (`spreadAllelesAcrossRows`), L36–L71, L210–L993, L232–L234, L967–L974.

**`combinatorial_motif_place` = no.**
*Evidence:* Positive-absence confirmed: `grep -E 'motif|spacer|permut|combinator'` returns **0 hits** across `README.md`, `ui.R`, `server.R`, `R/processVCFfast.R`, `man/processVCF.Rd`, and the full 28-page paper + supplement. There is no concept of an element, a position, an orientation, or a placement. The only positional control is flanking-context width (`upstreamContextRange`, `downstreamContextRange`) around the variant; the only orientation control is a per-variant strand flag (`MPRAREV`/`RV` in the VCF INFO field).
*Source:* full-source and full-PDF grep.

**`barcode_generation` = yes.**
*Evidence:* Genuinely a strength; the tool does more than a one-line "yes" suggests.
(a) **Original set.** Supplement S5: *"We generated the set of all possible DNA 12-mers then screened these according to the design parameters in Melnikov et al. … each nucleotide occurs at least once; no runs of single nucleotides greater than length 3; do not contain restriction sites for KpnI/XbaI; do not start with TCT …; they do not match any human miR seed sequences."* → 1,140,292 inert twelvemers (`man/twelvemers.Rd`). The KpnI/XbaI clause is prose only: the reproduced `generate12mers.R` has no `GGTACC`/`TCTAGA` removal step, and **3,767 of the shipped 1,140,292 twelvemers contain `GGTACC`** (none contain `TCTAGA`). The enzyme filtering `processVCF()` does at runtime is a separate, later operation.
(b) **Error-correcting sets.** v0.2.0 adds **27** `freebarcodes` sets (Hawkins et al., PNAS 2018) named `barcodes<length>-<n_correctable_errors>`, spanning **3–17 bp** with 1–2 correctable indel errors (`barcodes3-1` … `barcodes17-2`; only `barcodes17-2` exists at 17 bp). Verified against the repo `data/` tree, `man/`, and the README table. **(The original memo's "10–17 bp" was wrong.)**
(c) **Runtime re-filtering.** All-four-nucleotides (`ensure_all_4_nuc`), removal of `filterPatterns` (default `'AATAAA'`, the polyA signal), and removal of all three enzyme sites plus their reverse complements and reverses — though the removal step is guarded by `if (nrow(barcodeFilter) > 1)` (L1200), so a pool in which exactly one barcode matches is left unfiltered.
(d) **Library-wide uniqueness by construction.** `vcf %<>% mutate(bcPools = split(shuffled_mers, ceiling(1:length(shuffled_mers) %% nrow(vcf))))` gives every variant a **disjoint** barcode pool, then `sample(snp$bcPools %>% unlist, 2*nper, replace = FALSE)` draws within it. No barcode is reused anywhere in the library and the barcode→variant map is injective **without any post-hoc dedup step** — for the shipped sets, which are duplicate-free (verified). `split()` partitions positions rather than values, and the only `unique()` check is per-variant (L978), so duplicates inside a user-supplied custom vector would not be caught.
(e) **Documented quality audit of the borrowed sets.** `R/import_freebarcodes.R` contains commented-out code — not re-runnable as shipped: it reads hard-coded absolute paths (`/mnt/bigData2/...`, `~/plateletMPRA/...`) to files absent from the repo — measuring, per freebarcodes set, the fraction containing all four nucleotides, homopolymer runs ≥4, and human miR seeds, concluding: *"So none of the new sets have nucleotide runs, and the larger ones have all 4 nucleotides ~85% - 98% of the time"* — the basis of the man page's "varying degree" caveat.
Custom barcode vectors are accepted.
*Caveat (real):* `mpradesigntools` applies no GC-content filter of its own at any stage, and the `twelvemers` set it ships spans GC 2/12–10/12. The 27 borrowed freebarcodes sets are GC-balanced upstream, though — every barcode in every set sits within one or two bases of half GC — so for those the constraint is inherited rather than absent.
*Source:* Supplement S5 (merged PDF pp. 27–28 = supplement pp. 25–26); `man/twelvemers.Rd`; `R/processVCFfast.R` L1180–L1220, L394; `R/import_freebarcodes.R`; `README.md` barcode table; repo `data/` tree.

**`per_sequence_provenance` = partial.**
*Evidence:* Each output row carries structured build metadata beyond the mutation name: `type` (ref/alt), `allele`, `barcode`, `snpIndex`, `totIndex`, plus `site_fix_info` (a nested data frame recording exactly which bases were changed to destroy an aberrant restriction site, unnested into the TSV) and `notes`. But three limits keep it below `yes`: (i) the schema is a fixed flat column set for one design operation — there is no general "how it was built" record because there is only one way anything is ever built; (ii) the terminal `select()` **drops `CHROM` and `POS`**, so the delivered library carries **no genomic coordinates** and traceability back to the locus depends entirely on the VCF `ID` column being populated (`CHROM`/`POS` survive only in the `failed` table); (iii) `notes` is assigned once per variant (`res$notes = notes`, L987) and recycled across all `2*nper` rows, making it **variant-level, not per-sequence** — only `site_fix_info` genuinely varies across rows.
*Source:* `R/processVCFfast.R` L236, L287, L323, L986–L989, L1252, L1290.

**`automatic_naming` = partial.**
*Evidence:* The docs claim `processVCF` produces *"labeled MPRA sequences"* and `processSnp` returns *"a tibble of labeled sequences with appropriate information on the changes made"*. But labelling is a set of identifier **columns** (`ID`, `type`, `snpIndex`, `totIndex`), not a composed name string, and **no name-composition code exists anywhere** in `processVCFfast.R` (no `paste0` producing an identifier). `ID` is carried **verbatim** from the input VCF's ID column, so a VCF without rsIDs yields a library whose rows are identified only by integer indices. Exporting FASTA requires the user to build names themselves.
*Source:* `man/processVCF.Rd`, `man/processSnp.Rd`; `R/processVCFfast.R`.
*Note:* `partial` is the generous reading of "labeled"; a stricter matrix would say `no`.

**`design_visualization` = no.**
*Evidence:* `designMPRA/server.R` has **eight** `output$` assignments — `powerPlot`, `testHead`, `timeText`, `warnText`, `conditionalDownload`, `failed`, `inputHead`, `downloadSequences` — of which six are `render*`. The only plot is `powerPlot`, a `pwr.t.test` power curve: a statistics plot, not a view of the library. `testHead` is a leftover `head(mtcars)` debug stub and `downloadSequences` is the TSV `downloadHandler`. The construct-layout picture is a **static PNG** (`img(src = 'MPRAseqDiagram.png', width = 955, height = 20)`, `ui.R` L144), identical regardless of design — and that PNG is **not in the repo**, since `www/*` is gitignored. The R package has **no plotting function at all** (no `graphics`/`ggplot2` in `NAMESPACE`); the only visual feedback is a `print()`ed construct-size distribution table.
*Source:* `designMPRA/server.R` L24, 63, 67, 128, 136, 148, 152, 161; `designMPRA/ui.R` L144; `designMPRA/.gitignore`; `NAMESPACE`.

### Block B — assay coverage

**`assay_dms` = no.**
*Evidence:* No codon, amino-acid, ORF or reading-frame concept anywhere. `NAMESPACE` `Biostrings` imports are strictly nucleotide-level: `DNAString`, `DNAStringSet`, `complement`, `countPattern`, `replaceLetterAt`, `reverse`, `reverseComplement`, `subseq`, `toString` — nothing translational. `grep -E 'codon|amino acid|\bORF\b|reading[- ]?frame'` over the full source and the 28-page PDF: 0 hits. (The looser pattern `frame` does hit, but only as `data_frame` / `data.frame` / "data frame" / "web application framework".)
*Source:* `NAMESPACE`; full-source and full-PDF grep.

**`assay_mpra` = yes.**
*Evidence:* The tool's entire purpose. Title: "Design tools for MPRA experiments." `DESCRIPTION`: `Description: Design MPRA sequences from input VCF's`. Scope is narrow *within* MPRA: SNP/indel allelic-contrast MPRAs with the Melnikov et al. 2012 oligo layout. Notably, the README's "Planned Features" list leaves *"bed file to Sharpr-MPRA library oligo design"* **unchecked** — tiling/Sharpr-style MPRA design was contemplated and never implemented. (Also unchecked: mm10 genomic context, parallelization.)
*Source:* paper title/abstract; `DESCRIPTION`; `README.md` Planned Features.

**`assay_insilico` = no.**
*Evidence:* `DESCRIPTION` Imports are dplyr / purrr / readr / tidyr / tibble / stringr / purrrlyr / magrittr / Biostrings / BSgenome.Hsapiens.UCSC.hg38 only. No model, scorer, or ML dependency in the graph. The word "model" in the paper refers solely to the `pwr.t.test` power model.
*Source:* `DESCRIPTION`; `NAMESPACE`; PDF.

### Block C — genomics integration

**`genome_coordinates` = yes.**
*Evidence:* VCF `CHROM`/`POS` slice the reference directly: `snpseq = subseq(genome[[paste0('chr', as.character(snp$CHROM))]], ...)` where `genome = BSgenome.Hsapiens.UCSC.hg38` is assigned **inside `processSnp` (L227, L766) — hard-coded, not a parameter** — and `import(BSgenome.Hsapiens.UCSC.hg38)` is a wholesale NAMESPACE import. The paper frames this as its advantage over MPRAnator: *"Our tool acquires genomic context from the hg38 reference genome, rather than requiring input by the user."*
*Constraint:* **hg38 human only**; the README's Planned Features still lists *"mm10 genomic context"* unchecked.
*Source:* `R/processVCFfast.R` L227, L337, L555, L766, L777; `NAMESPACE`; `README.md`; PDF.

**`transcript_models` = no.**
*Evidence:* No GTF/GFF parsing, no annotation package in `DESCRIPTION` Imports. Strandedness is **user-supplied**, not derived: *"If the user wishes to generate SNPs for genes that normally are read from the reverse strand, add a string containing 'MPRAREV' to the INFO field of the VCF"* and *"the MPRAREV tag will need to be added by the user (where appropriate) because the VCF's do not always specify which strand the relevant gene is on."* In code: `reverseGene = grepl('MPRAREV', INFO)`.
*Source:* `README.md`; `R/processVCFfast.R` L1217.

**`exon_intron_split_codons` = no.**
*Evidence:* Context extraction is a flat genomic window — `subseq(genome[[chr]], POS - upstreamContextRange, POS + downstreamContextRange)` — with no splice-junction and no reading-frame awareness. Follows from `transcript_models = no` and `assay_dms = no`.
*Source:* `R/processVCFfast.R` L337, L555, L777.

**`hgvs_input` = no.**
*Evidence:* `processVCF`'s first act is `readLines(con = vcf)` then `readr::read_tsv` — a **file path** only. README: *"Only the CHROM, POS, REF, and ALT columns are used. The INFO column is used only for detecting reverse strand constructs."* Zero `HGVS` hits in code or paper. (The Shiny app's older `processVCF` takes an in-memory data frame instead, but it is still VCF-shaped.)
*Source:* `R/processVCFfast.R` L1099–L1140; `README.md`; `designMPRA/server.R` L114.

**`vcf_vep_output` = no.**
*Evidence:* VCF is the **input** format. Output is `write_tsv(successes, path = outPath)`, optionally a `save()`d `.RData`. Paper: *"users provide a VCF containing variants to receive a tab-separated file containing the MPRA sequences for their experiment."* The only VCF ever *written* is `spread_and_fix_indels()`'s normalized `*_fixed.vcf` copy of the **input** — a prestep, not a result. No FASTA or GenBank writer either.
*Source:* `R/processVCFfast.R` L1256–L1259 and L1294–L1330 (result writer / `.RData`), L1418–L1450 (`spread_and_fix_indels`, the `*_fixed.vcf` writer); `NAMESPACE` (`readr::write_tsv` is the only writer); PDF.

**`consequence_annotation` = no.**
*Evidence:* No molecular-consequence vocabulary. The only annotation-like strings are **construction events**: the `notes` strings (`'Shortened context by … bp on each side …'`, `'The alleles for this SNP were flipped from the input VCF because of the presence of the RV tag …'`) and the **nine** distinct `failed$reason` values (`'Failed - Context contained a digestion site'`, `'Failed - the variant is located within an aberrant digestion site'`, `'Failed - More than one aberrant digestion site in context'`, `'Failed - aberrant site generated at interface between sequence elements'`, `'Failed - SNP sequence could not be generated without an aberrant digestion site even after barcode resampling'`, `'Failed - SNP sequence could not be generated without multiple aberrant digestion sites'`, `'Failed - SNP sequence could not be generated without an aberrant digestion site in all constructs'`, `'Failed - SNP sequence could not be generated without an aberrant digestion site'`, `'Failed - Not identifiable as SNV, insertion, or deletion.'`). The nine failure strings are verbatim; the two `notes` strings above are elided — the first is a `paste0` template with an interpolated width.
*Source:* `R/processVCFfast.R` L287, L323 (`notes`); L137, L179, L361, L371, L448, L489, L513, L733, L973 (one location per distinct failure string).

### Block D — physical construction

**`primer_design` = no.**
*Evidence:* `fwprimer` and `revprimer` are **required user-supplied strings** ("a string containing the forward PCR primer to be used") pasted onto the construct ends. Supplement S1.3: *"enzyme1 and enzyme2 default to KpnI and XbaI as in Melnikov et al., Nature Biotechnology 2012. The forward and reverse primers default to ACTGGCCAG and CTCGGCGGCC respectively."* No Tm calculation, no primer3, no mutagenic-primer generation.
*Source:* Supplement S1.3 (PDF **p.5** of the merged 28-page file = supplement's own p.3); `man/processVCF.Rd`.

**`codon_optimization` = no.**
*Evidence:* No codon table, no translation, no optimization objective, and nothing in the dependency graph capable of providing one.
*Source:* `DESCRIPTION`; `NAMESPACE`; full-source grep.

**`synthesis_constraints` = partial.**
*Evidence:* Real, and more sophisticated than a one-liner suggests — but cloning-oriented rather than synthesis-oriented.
*What it does:* `countDigSites(DNAString(.x), enzyme1, enzyme2, enzyme3)` counts aberrant KpnI/XbaI/SfiI sites in each **fully assembled** construct, including sites created **at the junctions between elements** (`'Failed - aberrant site generated at interface between sequence elements'`). Enzyme patterns are IUPAC-aware (`fixed = FALSE`; *"The three enzyme arguments may contain ambiguous nucleotides"*). SfiI is checked but never inserted, deliberately: *"The sequence for enzyme3 does not show up in the output sequences, however it is necessary to check for it's presence in the output sequences as it is used when preparing the plasmid library."* On failure it retries with resampled barcodes under a **bounded budget of 40 attempts** per variant (`resample_attempts > 40`, three copies of the loop at L502, L723, L926, one per variant class). With `alter_aberrant = TRUE`, `randomly_fix()` samples `nrow(res_df)/2` single-base alterations of the offending site and hands one to each barcode replicate — the source comment says *"This assures that the changes are unique, if possible"*, and `sample_n(., nrow(res_df)/2, replace = ((nrow(res_df)/2) > nrow(.)))` deliberately reuses alterations once `nper` exceeds the number of possible one-base changes — while giving **matched** repairs to the ref and alt constructs (`fixed_pattern = rep(altered_patterns$altered_pattern, times = 2) # give the same altered patterns to both alleles`) so the allelic contrast is not confounded by the repair; `nonrandomly_replace_N()` protects the N positions of ambiguous (SfiI) sites; `check_cross_center()` hard-fails variants sitting inside the aberrant site, though its test is `start < center & (start + width) > center`, so a site beginning exactly at the variant position is not caught. Edits are recorded per base in `site_fix_info`. `max_construct_size` trims flanking context symmetrically until the oligo fits.
*What it does not do:* no GC content, no homopolymer runs, no repeats, no secondary structure, no synthesis-complexity scoring **of the construct** (homopolymer/miR/polyA screening applies to the **barcode** only). `alter_aberrant`, `extra_elements` and `max_construct_size` are all tagged **"under development"** in `man/processVCF.Rd`.
*Source:* `R/processVCFfast.R` L22–L35, L72–L123, L124–L184, L323, L502–L540, L723–L760, L926–L960, L1451–L1459; `man/processVCF.Rd`.

### Block E — engineering

**`interface` = yes.**
*Evidence:* R package API (single main function `processVCF()`) **plus** an R/Shiny web GUI. **No CLI** — verified absence of `inst/`, `exec/`, and any `Rscript` entry point in the recursive git tree. GUI limits: output capped just below 20,000 barcoded sequences — the UI prose says *"To avoid overloading our Shiny server … can be at most 20000"*, but `server.R` L97 validates `input$nBCperSNP*nrow(inVCF())*2 < 20000`, so exactly 20,000 is rejected — barcodes locked to 12 bp, and the app itself directs serious users to the package.
*Two verifiable caveats (the broad GUI-level claim is narrowed to the specific failures below):*
1. The **hosted** instance at https://andrewghazi.shinyapps.io/designmpra/ returns HTTP 404 (re-verified 2026-08-14; case variant and the shinyapps.io account root also 404). The app source can still be **launched** locally via `shiny::runApp()`, but its sequence-generation path cannot run from the distributed checkout: `scripts/processVCFfast.R` L383 does `load('outputs/inertTwelveMersChar.RData')`, `outputs/*` is gitignored, and that barcode object is absent from the repository tree — `scripts/shinyAppsDeploy.R` lists it as a separately supplied deployment asset.
2. But a local run is degraded and stale: **`www/*` is in `.gitignore`**, and the repo tree confirms there is no `www/` directory — so `MPRAdiagram.png`, `MPRAseqDiagram.png` and `formulas.png` (every image the UI references, including the construct-layout diagram) are **absent from the distributed source**; a local run renders broken images. And `server.R` sources a **2017 snapshot** of the engine (`scripts/processVCFfast.R`, 18,609 B vs the package's 62,351 B) whose `processVCF` signature has no `barcode_set`, `alter_aberrant`, `max_construct_size`, `flip_RV`, `filterPatterns` or `ensure_all_4_nuc`, and only a symmetric `seqwidth`. The GUI cannot reach any v0.2.0 feature.
*Source:* `designMPRA/ui.R` L108, L126, L144; `designMPRA/server.R` L17, L97, L114; `designMPRA/scripts/processVCFfast.R` L383; `designMPRA/scripts/shinyAppsDeploy.R`; `designMPRA/.gitignore`; GitHub recursive trees for both repos; live HTTP checks.

**`license` = yes.**
*Evidence:* `DESCRIPTION` declares `License: GPL-2`. Caveat: declared **only** there — the GitHub REST API returns `license: null` for **both** `mpradesigntools` and `designMPRA`, and neither recursive file tree contains a `LICENSE` or `COPYING` blob.
*Source:* `DESCRIPTION`; GitHub API `/repos/...`; recursive git trees.

**`maintained` = no.** **(CHANGED from `yes`.)**
*Evidence:* Here, **maintained** means ongoing code/release maintenance, not occasional issue responses. `mpradesigntools` last commit **2020-12-17T19:27:27Z** ("replace data_frame with tibble"; the GitHub API's `pushed_at` of 19:27:54Z is a push, not a commit, timestamp) — **5 years 8 months** dormant as of 2026-08-14. **Zero tags, zero releases**, one branch, 3 stars, 1 fork, 0 open issues — but five issues have been filed (#1–#5, 2018 through 2024), all now closed; the author replied to every one, closed four himself, and answered as recently as **2024-10-27**. In that reply he wrote: *"I should probably archive this repository soon as I don't have adequate bandwidth to maintain it at this point unfortunately."* `designMPRA` last commit **2017-09-26** — 8 years 10 months dormant, 1 star. Never published to CRAN (both the package page and the CRAN Archive path 404), never on Bioconductor (absent from the release package indexes checked 3.6–3.14 and from the removed-packages list), 0 results on r-universe. The only fork of the package, `mmandelia/mpradesigntools` (2 commits, 2023-03-24), is a private hack, not a successor; `goldenac/MPRADesignGenerator` is confirmed independent (`fork: false`, no shared history). No maintained descendant exists. The advertised web service is dead. Both repos are public and unarchived, so **"source still retrievable" is true — but that is a different claim from "maintained."**
*Source:* GitHub API `/repos` for both repos and `/tags`, `/releases`, `/forks`, `/issues?state=all`, `/issues/5/comments` for `mpradesigntools`; local clone at `afd386e` for the HEAD commit timestamp; CRAN/Bioconductor/r-universe checks.

---

## 4. The tool's own documented examples / vignettes

There is **no `vignettes/`, no `tests/`, no `inst/extdata`, and no example VCF** shipped with the package (verified against the full recursive git tree). What exists:

1. **README usage example** — the only complete invocation anywhere, and **the single reproduction target**:
   ```r
   processVCF(vcf = '/path/to/the.vcf', nper = 14,
              upstreamContextRange = 55, downstreamContextRange = 55,
              outPath = '/path/to/the/output.tsv',
              fwprimer = 'ACTGGCCGCTTCACTG', revprimer = 'AGATCGGAAGAGCGTCG',
              alter_aberrant = TRUE, extra_elements = FALSE,
              max_construct_size = 170, barcode_set = 'barcodes14-1',
              ensure_all_4_nuc = TRUE)
   ```
   A 170 bp oligo, 14 barcodes per allele, 55 bp of hg38 context each side, error-correcting 14-mer barcodes (`barcodes14-1`, **157,197** barcodes; the README table says 157,196 — it was generated with `wc -l` over source files lacking a trailing newline, so all 27 of its counts are one short).
2. **Shiny "Create sequences" defaults, stated to mirror a real published library:** *"The default values are comparable to those used in Ulirsch et al. 2016"* — 14 barcodes per **allele** (the control is labelled "Number barcodes per SNP", but the value is passed as `nper` and the engine emits `nper` ref plus `nper` alt, i.e. 28 per ordinary biallelic SNP), 75 bp genomic context each side, KpnI/XbaI/SfiI, primers `ACTGGCCAG` / `CTCGGCGGCC`.
3. **The paper's Sequence Element Layout**, with Melnikov et al. 2012 enzyme/primer defaults (Supplement S1.3, merged PDF p.5, plus the `MPRAseqDiagram.png` schematic — which is itself missing from the repo). The published layout inserts `TG` before the genomic context and `GGC` before the reverse primer; the 2017 web worker always emits those, but in the package `extra_elements` defaults to `FALSE`, so they are **not** part of the package's default construct — only `extra_elements = TRUE` reproduces the published layout.
4. **`generate12mers.R`** (Supplement S5, full code reproduced in the paper's supplement, merged PDF pp. 27–28 = supplement pp. 25–26; also present as `designMPRA/scripts/generate12mers.R`): all 12-mers → require all four nucleotides → drop runs ≥4 → drop 12-mers starting `TCT` → drop human miRBase seed matches → 1,140,292 barcodes. There is **no** KpnI/XbaI removal step in the script, despite the S5 prose (see `barcode_generation`). Not re-runnable as shipped: it reads `~/plateletMPRA/mirBaseSpecies.txt` and `~/plateletMPRA/mature.fa` and writes under `~/designMPRA/outputs/`, none of which are in the repo.
5. **Power-analysis vignette (not library design):** Supplement S2/S3 re-analyse Tewhey et al. 2016 and Ulirsch et al. 2016 barcode counts to estimate activity SD and justify the t-test; S4 gives the `pwr.t.test` call driving the app. **Correct number to quote: the combined-study mean activity SD is ≈ .95** (paper p.2: *"on average around .95"*; supplement p.14: `## [1] 0.9523771`). The figure **.926 is Tewhey-only** (supplement p.9: *"the mean of every allele in every transfection in this study is .926"*); Ulirsch-only is **1.1221643** (supplement p.14). Note the authors themselves mislabel .926 in the S3.3.1 simulation code (merged PDF p.23, `sd = .926) # using the average SD from Ulirsch and Tewhey together`) — quote **.95** for the combined figure.

---

## 5. Notable capabilities not in the row list

- **Statistical power analysis** — the paper's headline feature, entirely orthogonal to sequence design. Interactive `pwr.t.test` power curves over transcriptional shift, parameterized by barcodes/allele, transfection replicates, number of variants (Bonferroni), alpha, and assumed activity SD, **calibrated on real MPRA data** (Tewhey 2016, Ulirsch 2016; combined mean allele activity SD ≈ .95, central 95% interval 0.453–1.587). Validated against Monte Carlo simulation: max abs. theoretical-vs-simulated power difference **0.0400** for the most non-normal allele (S3.3.1, merged PDF p.24) and **0.00341** for the normal control allele (S3.3.2, merged PDF p.26).
- **Error-correcting barcode sets.** 27 `freebarcodes` sets indexed by length (**3–17 bp**) and number of correctable indel errors (1–2), selectable via `barcode_set`; plus the original 1,140,292 inert twelvemers.
- **Library-wide barcode uniqueness by construction.** Disjoint per-variant barcode pools (`split(shuffled_mers, ceiling(1:length(shuffled_mers) %% nrow(vcf)))`) + `sample(replace = FALSE)` within each pool make the barcode→variant map injective across the whole library with **no post-hoc dedup step** — for the shipped sets, which are duplicate-free; `split()` partitions positions rather than values, so duplicates in a user-supplied custom vector are not caught. Distinct from "generating a barcode set."
- **Per-replicate randomised constraint repair.** `randomly_fix()` samples `nrow(res_df)/2` single-base alterations of the offending restriction site and gives one to each barcode replicate — unique *"if possible"*, falling back to sampling with replacement once `nper` exceeds the number of possible one-base changes — while giving **matched** repairs to ref and alt so the allelic contrast is not confounded; `nonrandomly_replace_N()` protects N positions of ambiguous (SfiI) sites; `check_cross_center()` hard-fails variants inside the aberrant site (except one beginning exactly at the variant position, which its `start < center` test misses); per-base edits recorded in `site_fix_info`. Materially more sophisticated than "repairs restriction sites."
- **Bounded retry budget.** Up to **40** barcode-resampling attempts per variant (three copies of the loop, one per variant class) before emitting the resampling-exhausted failure.
- **Structured failure reporting.** A second `failed` tibble with an explicit machine-readable `reason` per rejected variant (nine distinct strings), rather than silently dropping. `failed` retains `CHROM`/`POS`, which `result` does not.
- **Documented barcode-quality audit of the borrowed sets.** `R/import_freebarcodes.R` measures, per freebarcodes set, the fraction with all four nucleotides, homopolymer runs ≥4, and human miR seeds — *"none of the new sets have nucleotide runs, and the larger ones have all 4 nucleotides ~85% - 98% of the time"* — the basis of the man page's "varying degree" caveat. Cite this if the referee response characterises their barcode quality; it deserves respect, not dismissal.
- **VCF normalization utility.** `spread_and_fix_indels()` splits multi-allelic rows and rewrites indels to the dash convention, emitting `*_fixed.vcf`.
- **Automatic strand handling.** `flip_RV = TRUE` **complements** alleles carrying dbSNP's `RV` tag — the man page says "reverse complement", but the implementation (L291–L313) calls `Biostrings::complement()`, never `reverseComplement()`, so multibase indel alleles are complemented without reversing base order; `MPRAREV` flips construct orientation and swaps the up/downstream context semantics.
- **Runtime/feasibility estimation.** The app pre-computes and displays an expected job time (`nSnp * nBCperSNP * 10/1000 + 5` seconds — which omits the ×2 for ref/alt that both the generator and the cap check apply) and warns when the request exceeds the 20,000-sequence Shiny-server threshold. It does **not** warn about barcode-pool exhaustion: the worker's own 1,140,292 check (`scripts/processVCFfast.R` L379) sits behind the far lower GUI cap and is unreachable.
- **Named downstream analysis companion**, `malacoda` (QC and statistical analysis, per the package README) — a design→analysis pairing PoolParty does not attempt.

### Behaviour worth knowing but undocumented by the authors

- **Multi-allelic expansion duplicates the reference arm.** `spreadAllelesAcrossRows()` emits one row per ALT, and each row independently generates `nper` ref + `nper` alt constructs — so a tri-allelic site contributes **two full reference sets**. This affects library-size accounting and is documented nowhere.
- **App/package divergence** (see §2 and `interface`) — the two published artefacts are two different programs.

---

## 6. Stated limitations (the authors' own, plus verified structural ones)

**Authors' own:**
- *"We would like to emphasize that this sequence generation functionality of the web tool is mainly for simple demonstrative purposes"* (`ui.R`) — the GUI is explicitly not for production libraries.
- GUI hard caps: ≤ 20,000 barcoded sequences ("to avoid overloading our Shiny server"), barcodes locked to 12 bp.
- `alter_aberrant`, `extra_elements`, `max_construct_size` all tagged **"under development"** in `man/processVCF.Rd`.
- *"freebarcodes sets only meet the traditional MPRA barcode requirements to varying degree"* (`man/processVCF.Rd`).
- README "Planned Features", still **unchecked**: mm10 genomic context, parallelization, bed-file → Sharpr-MPRA library oligo design.

**Structural, verified:**
- hg38 human only — `BSgenome.Hsapiens.UCSC.hg38` is a hard NAMESPACE `import()` and is assigned inside `processSnp`, not exposed as a parameter.
- One fixed construct architecture; changing it requires forking the source (as the only fork in fact did).
- No mutagenesis schemes of any kind; ref-vs-alt only.
- No motif/element/spacer/combinatorial vocabulary anywhere.
- No GC, Tm, secondary-structure or repeat checking of the assembled construct, and no GC or Tm check anywhere. Barcode repeat screening does exist, but only upstream of the tool: homopolymer runs ≥4 were removed when `twelvemers` was built and are absent from all 27 freebarcodes sets. Nothing re-checks at design time.
- No FASTA/GenBank writer; TSV (+ optional `.RData`) only.
- Delivered `result` table carries no genomic coordinates.
- No validation that the input VCF's `REF` allele matches the hg38 base at `CHROM`/`POS`; the fetched reference sequence is labelled `ref` regardless.
- Output is stochastic and not reproducible on re-run: barcode shuffling, per-variant pool assignment, barcode sampling and `randomly_fix()` all sample, and `processVCF()` exposes no seed argument (there is no `set.seed` anywhere in the source).
- No CLI; no plugin/extension mechanism; no S4 generics; no user hooks.
- Documented install path is stale: README still instructs `source("https://bioconductor.org/biocLite.R"); biocLite(...)`; that URL now returns HTTP 404, while [current Bioconductor instructions](https://bioconductor.org/install/) use `BiocManager::install()`. Active code calls `tidyr::unnest_legacy()` ([legacy interface](https://tidyr.tidyverse.org/reference/nest_legacy.html)), `purrrlyr::by_row()`, and `readr::write_tsv(path=)` (`path` is [deprecated in favour of `file`](https://readr.tidyverse.org/news/index.html)); `BSgenome.Hsapiens.UCSC.hg38` is an imported dependency.

---

## 7. Availability today (re-checked 2026-08-14)

| Component | Status |
|---|---|
| Web app https://andrewghazi.shinyapps.io/designmpra/ | **DEAD — HTTP 404.** This is the primary availability link in the *Bioinformatics* paper. Case variant and the shinyapps.io account root also 404. |
| `github.com/andrewGhazi/mpradesigntools` | Alive, public, unarchived. Last commit **2020-12-17**. v0.2.0. **0 tags, 0 releases**, 1 branch, 3 stars, 1 fork, 0 open issues. |
| `github.com/andrewGhazi/designMPRA` | Alive, public, unarchived. Last commit **2017-09-26**. 0 tags, 0 releases, 1 star, 1 fork. Ships **without** its `www/` image assets (gitignored). |
| CRAN | **HTTP 404 — never published.** The package page and the CRAN Archive path `/src/contrib/Archive/mpradesigntools/` both 404. |
| Bioconductor | **Never indexed** — the package URL returns 200 only because it redirects to the generic "Removed Packages" page, which does not name it. Established instead by absence from the Bioc release package indexes (checked 3.6–3.14) and from the removed-packages list. |
| r-universe | 0 results. |
| Only documented install command | `devtools::install_github('andrewGhazi/mpradesigntools')` (README). Distribution is GitHub-only; a checkout can equally be installed with `remotes`/`pak` or `R CMD INSTALL`. |
| Installability risk | **Non-trivial, not verified** (no install attempted, per instructions). The documented `biocLite.R` URL returns 404; active code uses the legacy `unnest_legacy()` interface and deprecated `write_tsv(path=)` argument; `BSgenome.Hsapiens.UCSC.hg38` is an imported dependency. Primary documentation is cited in §6. |

**Bottom line:** the R package source is fetchable and intact, but the interactive web tool the paper advertises is gone, the package was never distributed through CRAN or Bioconductor, its documented Bioconductor bootstrap URL now returns 404, it has zero releases, and it has not been touched in 5 years 8 months.

---

## 8. Unresolved disagreements and confidence

**No capability cell is set to `unknown`.** Extractor and reviewer agreed on every capability value except `maintained` (reviewer right — flipped to `no`) and one sub-claim of `interface` (the broad claim was narrowed to the verified hosted and local sequence-generation failures). Both changes are in the direction of internal consistency, and the `maintained` change is *generous-to-competitor* in no direction that creates an attack surface.

**Two reviewer corrections that I checked and did NOT adopt** (flagged here so nobody re-litigates them):

1. **Supplement S1.3 page number.** The reviewer said S1.3 is on PDF p.4, not p.5. Re-extracted with PyMuPDF: PDF page 4 contains only "S1 - MPRA diagrams / S1.1 MPRA experimental diagram"; the S1.3 heading and the primer/enzyme default sentence are on **PDF page 5** (which is the supplement's own page 3). **The original extraction's "p.5" was correct**, given the merged 28-page file. Cite as "Supplement S1.3 (merged PDF p.5 = supplement p.3)" to be unambiguous.
2. **Activity-SD attribution, partially.** The reviewer is right that **.95 is the correct combined figure** and that **.926 is Tewhey-only** (supplement p.9). But the reviewer's claim that ".952 is the Ulirsch figure" is wrong: **Ulirsch-only is 1.1221643** (supplement p.14) and **0.9523771 is the combined mean of both studies** (supplement p.14). Adopted: quote **≈ .95 combined**. Note the authors' own S3.3.1 simulation-code comment (merged PDF p.23) mislabels `.926` as "the average SD from Ulirsch and Tewhey together" — an internal inconsistency in their supplement, not in ours.

**Least-confident cells** (all deliberately generous toward the competitor):
- **`library_as_object` = partial** — judgement call. One call yields the whole library so the user writes no loop, but the product is a bare tibble with no class or methods. A referee could argue `yes` ("it does make the whole library at once") or `no` ("there is no library abstraction"). `partial` is the defensible middle.
- **`per_sequence_provenance` = partial** — metadata is real and structured, but fixed-schema, coordinate-free, and `notes` is variant-level. Arguably could be tightened; left generous.
- **`automatic_naming` = partial** — labelling is identifier columns, not a composed name; `ID` is verbatim from the input VCF. A stricter matrix would say `no`.
- **`synthesis_constraints` = partial** — well calibrated. The restriction-site machinery is genuinely sophisticated (matched per-replicate repair, 40-attempt retry, junction checking, IUPAC awareness); the gap is that nothing checks GC, homopolymers, repeats or structure, and three relevant arguments are "under development."
- **Installability** — asserted from stale README instructions and deprecated dependency calls, **not** from an actual install (out of scope this round).

**Highest confidence:** Blocks B/C/D, `dag_chaining`, `combinatorial_motif_place`, `assay_dms`, `barcode_generation`, `maintained`, `license` — read directly from the full 1459-line source, `NAMESPACE`, all man pages, the recursive git trees, the GitHub REST API, and live HTTP checks, all of which are small enough to be exhaustive. Every `no` in this record survived a positive-absence check by two independent readers.
