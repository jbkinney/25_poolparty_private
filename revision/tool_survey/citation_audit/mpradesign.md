# Citation-integrity audit — `mpradesign`

Record audited: `25_poolparty_private/revision/tool_survey/final/mpradesign.md`

Repository revisions checked:

- `andrewGhazi/mpradesigntools` at `afd386ef12051bb0a37ad63a6f92acd555246584`
- `andrewGhazi/designMPRA` at `0cf56ee602fc86dde705906d071a39cbdf6e99a8`
- The cited merged 28-page PDF, extracted with PyMuPDF
- Live URLs and GitHub/r-universe APIs on 2026-08-14

## NOT FOUND IN ANY SOURCE

No wholly invented quotation was found. The most serious implementation-evidence failure is item 1: the record attributes KpnI/XbaI removal to `generate12mers.R`, but that operation is absent from both copies of the script and the shipped data contradict the KpnI claim.

## Every item not verified

1. **misquoted** — Record L107 repeats the supplement's restriction-site claim, and L222 says the reproduced `generate12mers.R` workflow includes “drop KpnI/XbaI sites.” Neither `designMPRA/scripts/generate12mers.R` nor the code reproduced on merged-PDF pp. 27–28 contains a step that removes either motif. The script removes missing nucleotides, homopolymer runs, a leading `TCT`, and human miRNA seeds. Independent parsing of `data/twelvemers.rda` found 1,140,292 entries, including **3,767 containing `GGTACC` (KpnI)**. It contains zero `TCTAGA`, but no XbaI-removal operation appears in the script. The supplement prose lists absence of KpnI/XbaI sites as a desired parameter; the reproduced implementation does not implement it. Runtime `processVCF()` later filters enzyme motifs, which is a separate operation.

2. **wrong-location** — Record L114 and L222 cite Supplement S5 as “PDF pp. 25–26.” In the merged 28-page file, S5 starts on **PDF p. 27** and continues on **PDF p. 28**. Pages 25–26 contain the normal-allele Monte Carlo analysis and the start of S4. “25–26” are the supplement's printed internal page numbers, not merged-PDF page numbers; elsewhere the record explicitly distinguishes these numbering systems.

3. **number-wrong** — Record L219 says `barcodes14-1` contains **157,196** barcodes. The packaged object and `man/barcodes14-1.Rd` both contain **157,197**. This is systematic: every one of the 27 freebarcodes counts in the package README is one too small. `R/import_freebarcodes.R` generated that table with `wc -l`, while every upstream text file lacks a final newline; record count (`awk ... NR`), each packaged RData object, and every generated man page are all README count + 1. For example, `barcodes3-1` has 2, not 1, and `barcodes17-2` has 23,025, not 23,024.

4. **minor-discrepancy** — Record L113 says there is “no GC-content constraint anywhere.” The borrowed freebarcodes sets are documented by their cited upstream repository as being generated with a **balanced-GC filter of 40–60%**. `mpradesigntools` itself does not perform a runtime GC filter, and the original `twelvemers` script has no GC filter, but “anywhere” is false for 27 shipped sets.

5. **minor-discrepancy** — Record L109 says runtime filtering removes every `filterPatterns`/enzyme match. At `R/processVCFfast.R:1200`, removal is guarded by `if (nrow(barcodeFilter) > 1)`. If exactly one barcode matches, it is not removed. The general filtering description is accurate for zero or at least two matches but not for the one-match case.

6. **minor-discrepancy** — Record L110 and L231 claim library-wide barcode injectivity “by construction” for the tool without qualification. The guarantee holds for the unique shipped sets, but custom character vectors are accepted unchanged (`R/processVCFfast.R:1158–1164`). `split()` separates positions, not values; duplicate strings in a custom vector can enter different variant pools. The only uniqueness check is inside each individual `processSnp()` result (`R/processVCFfast.R:977–980`), not across the final library. There is no global `unique()`/dedup check.

7. **number-wrong** — Record L90 calls `processSnp()` a “~890-line function, L210–L1098,” and L99 repeats that range. The function starts at L210 and closes at L993: **784 lines inclusive**, not ~890. L994–L1098 are blank lines and the roxygen documentation for the next function, `processVCF()`.

8. **wrong-location** — The source range for `lazy_evaluation` at record L95 (`R/processVCFfast.R` L1180–L1220) does not contain several cited facts in L94. `mers = get(barcode_set)` is L1158–L1159; the quoted insufficient-barcode error is L1167–L1172; VCF reading/materialization is L1121–L1145; and result `Reduce('rbind')` calls are L1243–L1290. L1180–L1220 covers the tail of barcode filtering and pool construction only.

9. **minor-discrepancy** — Record L117 and L298 say `notes` is assigned once per variant and broadcast uniformly to all `2*nper` rows. If both the RV note (L287) and construct-shortening note (L323) apply, `notes` is a length-2 vector. `res$notes = notes` at L987 recycles those two messages alternately across rows; it does not broadcast one combined variant-level value to every row.

10. **number-wrong** — Record L132 states that `grep` for `codon|amino|ORF|frame` over the full source and PDF gives 0 hits. The literal stated pattern has many hits: `data_frame`, `as.data.frame`, the PDF phrase “This data frame is then plotted,” and the paper phrase “web application framework for R.” A narrower semantic pattern such as `codon|amino|\bORF\b|reading[- ]?frame` has zero relevant hits, so the biological conclusion is supported, but the reported grep count is not.

11. **misquoted** — Record L155 presents `subseq(genome[[chr]], POS - upstreamContextRange, POS + downstreamContextRange)` as code and cites L337/L555/L777. That text does not occur at those locations. The implementation first assigns `rangestart` and `rangeend`, optionally swaps their context-width semantics, and then calls `subseq(genome[[paste0('chr', as.character(snp$CHROM))]], start = rangestart, end = rangeend)`. The record's expression is accurate pseudocode but is not a verbatim code quotation.

12. **wrong-location** — Record L164 cites `R/processVCFfast.R` L1418–L1450 for the result writer described in L163 as `write_tsv(successes, path = outPath)`. That exact result writer is at L1256–L1259 and again in L1294–L1330. L1418–L1450 is `spread_and_fix_indels()` and writes the normalized input VCF, not the MPRA result table. The PDF and NAMESPACE independently support TSV output.

13. **misquoted** — Record L167 says its two abbreviated `notes` strings and five failure strings were “all read verbatim from source.” The note strings containing the Unicode ellipsis (`…`) do not occur verbatim; the source contains full sentences at L287 and L323. Elliptical paraphrase is reasonable, but it is not verbatim as claimed.

14. **number-wrong** — Record L167 says there are five `failed$reason` values. There are **nine distinct `Failed - ...` strings** in `R/processVCFfast.R`: context site; more than one site in context; variant inside site; interface-created site; multiple aberrant sites; aberrant site in all constructs; explicit exhaustion after resampling; generic aberrant-site failure after insertion/deletion resampling; and unidentifiable variant class.

15. **minor-discrepancy** — Record L182 and L232 say `randomly_fix()` gives a different single-base repair to every barcode replicate. Source L150–L158 says “unique, if possible” and explicitly sets `replace = TRUE` when `nrow(res_df)/2` exceeds the number of possible one-base alterations. Matched ref/alt repairs are verified; unconditional uniqueness across replicates is not.

16. **number-wrong** — Record L189 and L252 describe the GUI limit as at most/≤ **20,000** sequences. `server.R:97` validates `input$nBCperSNP*nrow(inVCF())*2 < 20000`, so a design of exactly 20,000 is rejected. The UI prose and warning use ≤/> 20,000, so the upstream app itself is internally inconsistent. For ordinary integer inputs the generated total is even, making 19,998 the largest attainable total below the strict threshold.

17. **minor-discrepancy** — Record L191 says the distributed app source is “still runnable locally via `shiny::runApp()`,” with only missing images described as degradation. The UI/power portion may launch if dependencies are installed, but the “Create sequences” path calls `load('outputs/inertTwelveMersChar.RData')` at `designMPRA/scripts/processVCFfast.R:383`. `outputs/*` is gitignored and that RData file is absent from the recursive repository tree. Sequence generation is therefore not runnable from the distributed source alone, independently of the missing images.

18. **minor-discrepancy** — Record L200 labels `2020-12-17T19:27:54Z` as the package's “last commit” time. `19:27:54Z` is GitHub REST's repository `pushed_at` value. The actual HEAD author and committer timestamp for `afd386e` is **2020-12-17T19:27:27Z**. The date and commit subject are correct; the exact timestamp is a push timestamp, not a commit timestamp.

19. **minor-discrepancy** — Record L200 adds “0 open issues (nobody files any).” Zero open issues is verified, but “nobody files any” is false: the repository has **five issues**, all closed (#1–#5), filed from 2018 through 2024.

20. **number-wrong** — Record L220 calls the Shiny default “14 barcodes per SNP.” The UI labels it that way, but passes the value as `nper`, and the engine emits `nper` ref plus `nper` alt constructs. The cap calculation also multiplies it by 2. The behavior is **14 per allele / 28 per biallelic SNP**, not 14 total per SNP.

21. **wrong-location** — Record L223 and L294 say the authors' `.926` “average SD from Ulirsch and Tewhey together” comment is in S4. It is on merged PDF p. 23 in **S3.3.1, Non-normal allele** (and in `Rmd/Supplement.Rmd:627`). S4 begins on merged PDF p. 26 and continues on p. 27.

22. **number-wrong** — Record L229 says the Monte Carlo validation has “max abs. difference 0.00341.” `0.00341` is the result for the deliberately normally distributed control allele (S3.3.2, merged PDF p. 26). The primary highly non-normal-allele simulation (S3.3.1) reports **0.0400** on merged PDF p. 24. Across the two reported validation simulations, the maximum is 0.0400, not 0.00341.

23. **number-wrong** — Record L234 repeats that structured failure reporting has “five distinct strings.” The source contains **nine** distinct failure strings, as enumerated in item 14.

24. **misquoted** — Record L237 says `flip_RV = TRUE` “reverse-complements alleles” carrying `RV`. The documentation says reverse complement, but the implementation at `R/processVCFfast.R:291–313` calls `Biostrings::complement()`, not `reverseComplement()`, for SNVs, insertions, and deletions. Those operations are equivalent only for one-base alleles; multibase RV indels are complemented without reversing base order.

25. **minor-discrepancy** — Record L238 says the app “warns when the design exceeds the available barcode pool.” `warnText` at `server.R:128–134` warns only when the requested total exceeds the 20,000 server-workload threshold. The old engine separately stops above 1,140,292 barcodes at `scripts/processVCFfast.R:379–381`; that is not the displayed warning and is unreachable after the GUI's much lower validation cap.

26. **minor-discrepancy** — Record L262 says there is no repeat checking “on either barcode or construct.” Barcode sources do perform repeat checks: `generate12mers.R:29–33` removes homopolymer runs of length ≥4, and the cited freebarcodes generator documents a no-homopolymer-triples filter. No general repeat/secondary-structure check is applied to assembled constructs, but the blanket barcode statement is false.

27. **wrong-location** — Record L38/L200/L278 infers “never indexed” on Bioconductor from `https://bioconductor.org/packages/mpradesigntools/`. That URL redirects to the generic “Removed Packages” page. The page establishes that no current package page exists; it does **not** say that `mpradesigntools` was never indexed, and a generic removed-package fallback cannot by itself distinguish never-indexed from formerly indexed-and-removed. The historical claim needs a release/archive index citation.

28. **minor-discrepancy** — Record L280 calls `devtools::install_github('andrewGhazi/mpradesigntools')` the “Only install path.” It is the only installation command documented in the package README, but it is not the only technical path: the live source can be cloned/downloaded and installed with standard R tooling, or installed from GitHub with other clients such as `remotes`/`pak`. “Only documented install command” would match the evidence.

29. **uncited** — Record L274 says “No alternative hosted instance found.” No search query, search-result URL, registry, or dated candidate list is cited, so this negative web-wide claim is not independently reproducible from the record. The three specifically cited Shiny URLs do all return 404.

30. **uncited** — Record L281 says the package is “Likely still runnable on a pinned 2020-era toolchain” while also stating that no install was attempted. No environment, dependency lockfile, installation log, or execution result supports that claim.

31. **minor-discrepancy** — The app quotation at record L51 and L251 begins *“we would like to emphasize...”*; `ui.R:108` begins **“We would like to emphasize...”**. The text is otherwise verbatim. The unmarked capitalization change is minor but fails a strict verbatim check.

32. **minor-discrepancy** — The first README quotation at record L151 changes the source's double quotes around **"MPRAREV"** to single quotes **'MPRAREV'**. The words are otherwise verbatim; this is a punctuation-only discrepancy.

33. **uncited** — Record L239 identifies `malacoda` as a Bayesian MPRA-analysis package published in *PLOS Computational Biology* in 2020. The claim is true—the live `malacoda` README gives that description and cites DOI `10.1371/journal.pcbi.1007504`—but neither that repository nor the paper/DOI is cited in this record. The cited `mpradesigntools` README supports only the weaker statement that `malacoda` provides downstream QC and statistical analysis.

34. **uncited** — Record L266 gives several external, version-sensitive facts without citations: removal of `biocLite` in 2018 and replacement by `BiocManager`; deprecation of `tibble::data_frame`; deprecation of `readr::write_tsv(path=)` in favor of `file=`; and the ~700 MB size of `BSgenome.Hsapiens.UCSC.hg38`. The repository verifies that these APIs/dependencies are used, but it cannot verify their external deprecation dates/status or download size. The size is currently 731,245,539 bytes and the Bioconductor 3.8 materials support the 2018 BiocManager transition, but those primary sources are absent from the record.

35. **wrong-location** — Record L37/L200/L277 uses the current CRAN package URL's HTTP 404 to support “never published.” That URL verifies current absence only; a 404 package page does not itself state historical publication status. The canonical CRAN archive path also currently returns 404, which supports the inference, but the archive check is not cited in the record.

## URL check result

No unexpected dead link was found. The URLs the record labels dead still return 404 (Shiny app, CRAN package page, and `biocLite.R`); the two GitHub repositories, rdrr mirror, freebarcodes repository, malacoda repository, and NCBI documentation URL are live. The DOI resolves to the Oxford Academic article (the automated request receives HTTP 403 at the publisher, not a dead-link response). The r-universe API returned an empty `results` array for `mpradesigntools`.
