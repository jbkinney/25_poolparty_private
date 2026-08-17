# Repair change log — `mpradesign`

Record repaired in place: `final/mpradesign.md`
Audits adjudicated: `citation_audit/mpradesign.md` (35 items), `factcheck/mpradesign.md` (A1–A13, B1–B7, C)
Date: 2026-08-14

## How verification was done

Every finding was checked against primary source before any edit. No finding was
applied on the audit's assertion alone.

- **Repos**: both cloned fresh to scratch and confirmed at the revisions the audits cite —
  `andrewGhazi/mpradesigntools` @ `afd386ef12051bb0a37ad63a6f92acd555246584`
  (HEAD author/committer date `2020-12-17T14:27:27-05:00` = `19:27:27Z`),
  `andrewGhazi/designMPRA` @ `0cf56ee602fc86dde705906d071a39cbdf6e99a8`.
  All line numbers below were read out of those working trees.
- **PDF**: the supplied 28-page merged file re-extracted page-by-page with PyMuPDF, so every
  "merged PDF p. N" claim was checked against an explicit page boundary. Supplement page 5
  was additionally rendered to an image, because S1.3's sequence layout is a figure and does
  not survive text extraction.
- **Packaged data**: R is not installed here, so a minimal R-serialization (RDX2/XDR) reader was
  written in Python to open `data/*.RData` and `data/twelvemers.rda` directly. Every barcode
  count, GC range, homopolymer count and duplicate check below is a direct read of the shipped
  binary objects, not a repetition of the audit or of the README.
- **Live checks**: GitHub REST (`/repos`, `/tags`, `/releases`, `/forks`, `/issues?state=all`,
  `/issues/N/comments`), the hosted Shiny URL, CRAN package + archive paths, Bioconductor
  release `PACKAGES` indexes 3.6–3.14 and the removed-packages page, r-universe search,
  `biocLite.R`, the BSgenome tarball `Content-Length`, and the malacoda README/DOI.

Nothing outside `revision/tool_survey/` was written. Nothing was installed.

---

## APPLIED

Citation-audit item numbers are `C-n`; fact-check items are `F-An` / `F-Bn`.

| # | Finding | Correction made | How verified |
|---|---|---|---|
| C-1 / F-A1 | Record §3(a) and §4(4) credit `generate12mers.R` with a "drop KpnI/XbaI sites" step | §4(4): removed that step from the workflow chain and stated plainly that no such step exists despite the S5 prose. §3(a): appended one sentence recording that the KpnI/XbaI clause is prose only and that the shipped set does not honour it. The S5 **quotation itself was left untouched** — it is verbatim. | Read `designMPRA/scripts/generate12mers.R` in full and the code reproduced on merged PDF pp. 27–28: filters are missing-nucleotide, runs of 4, leading `TCT`, human miR seeds. No `GGTACC`/`TCTAGA` step in either copy. Parsed `data/twelvemers.rda`: 1,140,292 entries, **3,767 contain `GGTACC`**, 0 contain `TCTAGA`, 0 start with `TCT`, 0 runs ≥4, 0 missing a nucleotide. |
| C-2 | S5 cited as "PDF pp. 25–26" in §3 source line and §4(4) | Both changed to "merged PDF pp. 27–28 = supplement pp. 25–26" | Page-boundary extraction: `S5 - Barcode design` begins inside the merged-p.27 block and runs into p.28. Merged pp. 25–26 hold S3.3.2 and the start of S4. 25/26 are the supplement's own printed folios (visible at page foot). |
| C-3 / F-A9(1) | `barcodes14-1` given as 157,196 | Corrected to **157,197**, with the cause of the off-by-one noted | Parsed `data/barcodes14-1.RData` → length 157197; `man/barcodes14-1.Rd` → "length 157197". Parsed all 27 sets: **every one is exactly README + 1** (e.g. `barcodes3-1` = 2 not 1, `barcodes17-2` = 23,025 not 23,024, `barcodes16-1` = 1,636,418 not 1,636,417). `R/import_freebarcodes.R` L31–37 builds the README table with `wc -l`. |
| C-4 | "no GC-content constraint anywhere" | Caveat narrowed: no GC filter in `mpradesigntools` itself, `twelvemers` spans GC 2/12–10/12, but the 27 freebarcodes sets are GC-balanced upstream | Computed GC per barcode across all 28 shipped objects. `twelvemers`: GC 0.167–0.833, 377,563 outside 40–60%. Every freebarcodes set is confined to a narrow band around 50% (e.g. `barcodes14-1` 0.429–0.571, `barcodes10-1` 0.400–0.600, `barcodes16-1` 0.375–0.625). Verified empirically, **not** from the upstream doc the audit cited. |
| C-5 | Runtime pattern filtering described as unconditional | Added the guard clause | `R/processVCFfast.R` L1200: `if (nrow(barcodeFilter) > 1) { mers %<>% .[-barcodeFilter$removeIndex] }`. Exactly one match ⇒ no removal. |
| C-6 | Barcode injectivity asserted "by construction" without qualification (§3(d) and §5) | Both places qualified: holds for the shipped sets, which are duplicate-free; `split()` partitions positions not values; only per-variant `unique()` check | L1215–1216 `split(shuffled_mers, ceiling(...%% nrow(vcf)))`; L1158–1164 custom character vectors accepted unchanged; L978 `if (length(unique(res$barcodes)) != nrow(res))` is inside `processSnp`, so per-variant only. I additionally parsed all 28 shipped objects and confirmed **none contains a duplicate**, which is why the qualifier is "for the shipped sets" rather than a retraction. |
| C-7 | `processSnp()` called "a ~890-line function, L210–L1098" (twice) | Corrected to "784-line function, L210–L993" in both places | `processSnp = function(` at L210; its closing `}` at L993 (`return(res)` L992). L994–996 blank, L997 begins the roxygen block for `processVCF`, whose `function(` is L1099. 993−210+1 = 784. |
| C-8 | `lazy_evaluation` source range L1180–L1220 does not contain the cited facts | Source line list replaced with the ranges that actually hold each fact | VCF read L1121–L1145; `mers = get(barcode_set)` L1159; capacity `stop()` L1167–L1173; filters L1175–L1202; result `Reduce('rbind')` L1246/L1284 within L1243–L1290. |
| C-9 | `notes` said to be "broadcast to all `2*nper` rows" | One word: "broadcast" → "recycled" | `notes = c(notes, ...)` at both L287 and L323, so `notes` can be length 2; `res$notes = notes` (L987) then recycles. "Recycled" is correct for both the length-1 and length-2 cases. The record's substantive conclusion (variant-level, not per-sequence) is unaffected and was left alone. |
| C-10 | `grep 'codon\|amino\|ORF\|frame'` claimed to give 0 hits | Pattern replaced with the semantic one that genuinely returns 0, with a parenthetical explaining the loose-pattern hits | Literal pattern: 4 hits in `R/processVCFfast.R`, 1 in `NAMESPACE`, 2 in `server.R`, 9 in the PDF — all the substring `frame` (`data_frame`, `data.frame`, "data frame", "web application framework"). Narrow pattern `codon\|amino acid\|\bORF\b\|reading[- ]?frame`: **0** across source and PDF. Conclusion (`assay_dms = no`) stands. I also re-ran the neighbouring `motif\|spacer\|permut\|combinator` grep and confirmed it really is 0 hits, so the two claims are now consistent. |
| C-12 | `vcf_vep_output` source cited L1418–L1450 for the result writer | Source split: L1256–L1259 and L1294–L1330 for the result writer / `.RData`; L1418–L1450 retained and labelled as `spread_and_fix_indels`' `*_fixed.vcf` writer | `write_tsv(successes, path = outPath)` at L1258 and L1329; `save(res, ...)` L1296; `spread_and_fix_indels` starts L1418 and writes at L1447. The record's prose made both claims; only the citation was wrong. |
| C-13 | "all read verbatim from source" | Corrected: the failure strings are verbatim, the two `notes` strings are elided (the first is a `paste0` template) | Compared each quoted string to source. Five quoted failure strings verbatim at L137/L179/L361/L371/L513. `notes` strings at L287/L323 are longer sentences, and the shortening one interpolates `amount_to_remove`. |
| C-14 / C-23 / F-A6 | "five distinct" `failed$reason` strings (§3 and §5) | Corrected to **nine** in both places; the four missing strings added to the §3 list; source line list rebuilt to one line per distinct string | `grep -o "Failed - [^']*" R/processVCFfast.R \| sort -u` → exactly 9 distinct values across 18 occurrences. First occurrences L137, L179, L361, L371, L448, L489, L513, L733, L973. |
| C-15 / F-A5 | `randomly_fix()` said to give a **different** repair to every replicate; `check_cross_center()` said to hard-fail variants inside the site | Both qualified in §3 and §5: uniqueness is "if possible" with an explicit fallback to sampling with replacement; the cross-center test misses a site beginning exactly at the variant position | L150 comment "This assures that the changes are unique, if possible"; L156–158 `sample_n(., nrow(res_df)/2, replace = ((nrow(res_df)/2) > nrow(.)))`. L1451–1459: `any(start < center_point & (start + width) > center_point)` — `start == center_point` is not caught. The **matched** ref/alt claim was separately re-verified and left intact (L167 `rep(..., times = 2)` with the authors' own comment). |
| C-16 / F-A4 | GUI "capped at 20,000" | §3 `interface` changed to "capped just below 20,000", naming the actual predicate | `server.R` L97: `need(input$nBCperSNP*nrow(inVCF())*2 < 20000, ...)`. Strict `<`, so exactly 20,000 is rejected. **§6 left unchanged** — that bullet sits under "Authors' own" and accurately reports the UI's own prose ("can be at most 20000", `ui.R` L126). The app is internally inconsistent; the record now shows both sides. |
| C-17 / F-A3 | "The app source is still runnable locally via `shiny::runApp()`" | Narrowed: it can be *launched*, but the sequence-generation path cannot run from the distributed checkout | `designMPRA/scripts/processVCFfast.R` L383 `load('outputs/inertTwelveMersChar.RData')`; `.gitignore` L5 `outputs/*`; `git ls-files` confirms no `outputs/` and no `www/` tracked; `scripts/shinyAppsDeploy.R` L2–8 lists that RData and the three PNGs as separately supplied deployment assets. **See escalation 2.** |
| C-18 | `2020-12-17T19:27:54Z` labelled "last commit" | Corrected to `19:27:27Z`, with `19:27:54Z` named as the API `pushed_at` | Local clone HEAD `afd386e`: author and committer date `2020-12-17T14:27:27-05:00`. GitHub API `pushed_at` = `2020-12-17T19:27:54Z`. Date and commit subject were already right. |
| C-19 | "0 open issues (nobody files any)" | Replaced with the verified issue history | `/issues?state=all`: five issues #1–#5, created 2018-08-30, 2019-11-29, 2023-07-04, 2023-07-05, 2024-10-23, all closed. Comment threads show the author replied on every one, closed four himself, last reply 2024-10-27. **See escalation 1** — this is the one correction that plausibly bears on a capability value. |
| C-20 / F-A9(2) | Shiny default described as "14 barcodes per SNP" | Corrected to 14 per **allele**, with the UI's own label noted and 28-per-biallelic-SNP spelled out | `ui.R` L35–37 `numericInput('nBCperSNP', label = 'Number barcodes per SNP', value = 14)`; passed as `nper` (`server.R` L115); engine emits `nper` ref + `nper` alt; `ui.R` L126 itself defines the total as "number of SNPs times the number of barcodes per SNP times 2 for reference/alternate alleles". The record's §4(1) already said "per allele" for the README example, so this also removes an internal inconsistency. |
| C-21 | `.926` mislabel located "in S4" (§4 and §8) | Both corrected to the S3.3.1 simulation code, merged PDF p.23, with the comment quoted | Page-boundary extraction: `sd = .926) # using the average SD from Ulirsch and Tewhey together` sits in the merged-p.23 block, inside S3.3.1. S4 begins on merged p.26. |
| C-22 / F-A11 | Monte Carlo "max abs. difference 0.00341" | Both figures now given and attributed | Merged p.24 (S3.3.1, most non-normal allele): "maximum absolute difference between theoretical and simulated power of **0.0400**". Merged p.26 (S3.3.2, normal allele): "maximum absolute difference of **0.00341**". |
| C-24 / F-A7 | `flip_RV` said to "reverse-complement" | Corrected to "complements", with the man-page/implementation divergence noted | L286–L315: all three branches (SNV, insertion ALT, deletion REF) call `Biostrings::complement()`; `reverseComplement()` never appears in the RV block. Equivalent for one-base alleles, not for multibase indels. |
| C-25 / F-A13 | "warns when the design exceeds the available barcode pool" | Corrected: it warns on the 20,000 Shiny-server threshold; the pool check is unreachable. Also added that the displayed time estimate omits the ref/alt ×2 | `server.R` L128–133 `warnText` fires only on `> 20000`; L74 `timeGuess = round(nSnp*input$nBCperSNP*10/1000 + 5, 3)` with prose "This yields a total of `nSnp*input$nBCperSNP` sequences" — no ×2, while L97 validation and the generator both apply it. `scripts/processVCFfast.R` L379 caps at 1,140,292, far above the GUI cap. |
| C-26 / F-A8 | "No GC-content, Tm, secondary-structure or repeat checking on either barcode or construct" | Rewritten: no construct-level checking and no GC/Tm anywhere, but barcode homopolymer screening exists upstream | `generate12mers.R` L29–33 removes runs ≥4; measured 0 barcodes with a run ≥4 in `twelvemers` and in all 27 freebarcodes sets; `R/import_freebarcodes.R` L93–104 + L148 audit comment. This also removes a self-contradiction: §3 `synthesis_constraints` already said homopolymer screening "applies to the **barcode** only". |
| C-27 | Bioconductor "never indexed" inferred from a generic fallback page (§1 and §7) | Both entries now state that the URL redirects to the generic Removed Packages page, which does not name the package, and cite the checks that do establish it | `https://bioconductor.org/packages/mpradesigntools/` → 200 via redirect to `/about/removed-packages/`, `<title>Bioconductor - Removed Packages</title>`; that page contains 0 occurrences of "mpradesigntools". `Package: mpradesigntools` absent from the `PACKAGES` index of Bioc 3.6, 3.7, 3.8, 3.9, 3.10, 3.12, 3.14 and from the current release `VIEWS`. Conclusion unchanged; only its basis was wrong. |
| C-28 / F-A12 | Table row "Only install path" | Relabelled "Only documented install command", with GitHub-only distribution stated | README L27–37 documents only `devtools::install_github`. That is a claim about the README, not about technical possibility; `remotes`/`pak`/`R CMD INSTALL` on a checkout all work, and rdrr publishes the `remotes` form. |
| C-31 | App quotation began lowercase "we" (§2 and §6) | Capital "We" in both | `ui.R` L108 begins "We would like to emphasize". One-character fidelity fix, no nesting conflict. |
| C-35 | CRAN "never published" rested on the package-page 404 alone | Archive-path check added to the §7 row and to the `maintained` evidence | `https://cran.r-project.org/package=mpradesigntools` → 404; `https://cran.r-project.org/src/contrib/Archive/mpradesigntools/` → 404. |
| F-A2 | `generate12mers.R` called "a reproducible barcode-design use case"; `import_freebarcodes.R` called "reproducible (commented-out) code" | Both corrected to say they are not re-runnable as shipped, naming the absent hard-coded inputs | `generate12mers.R` L42/L48 read `~/plateletMPRA/mirBaseSpecies.txt` and `~/plateletMPRA/mature.fa`, L55/L62 write under `~/designMPRA/outputs/`. `import_freebarcodes.R` L6 reads `/mnt/bigData2/resources/freebarcodes/...`, L108 reads `~/plateletMPRA/data/mature.fa`, and the audit pipeline L130–146 is commented out. None of those inputs are in either repo. |
| F-A10 | "Melnikov et al. 2012 construct layout as the fixed default architecture" | §4(3) rewritten to name it as the paper's Sequence Element Layout and to state that `extra_elements` defaults to `FALSE`, so `TG`/`GGC` are not in the package's default construct | Rendered merged PDF p.5: the S1.3 figure reads `Forward PCR primer | TG | Genomic context with appropriate allele | enzyme1 | enzyme2 | Barcode | GGC | Reverse PCR primer`. Package signature L1110 `extra_elements = FALSE`; assembly L401–L419 adds `TG`/`GGC` only in the `TRUE` branch. Web worker `scripts/processVCFfast.R` L125–133 always emits them. S1.3's Melnikov attribution is specifically about the enzyme/primer defaults, which is why the heading was re-worded rather than deleted. |
| F-B1 | Variant-class limits never stated | One clause added to the existing `mixed_mutagenesis_one_pool` evidence sentence | L232 `isSNV` requires both `REF` and `ALT` in A/C/G/T; L233–234 `isINS`/`isDEL` require a literal `-`; L967–974 everything else returns `'Failed - Not identifiable as SNV, insertion, or deletion.'`. `fix_indels` L1391 `stop('cannot repair two multibase alleles')`. |
| F-B2 | No REF-vs-hg38 validation never stated | One bullet added to §6 "Structural, verified" | `snp$REF` is never compared against the fetched `snpseq`; grep over all 21 `snp$REF` uses shows classification, sizing, RV complementing and failure-row construction only. The record already covered the `ID`/`INFO`/`chr`-prefix half of this finding (§3 `automatic_naming`, `transcript_models`, `genome_coordinates`), so only the missing REF check was added. |
| F-B4 | Capacity check runs before the pool-shrinking filters | One clause added to the existing `lazy_evaluation` evidence | L1167 capacity `stop()` precedes L1175 `ensure_all_4_nuc` and L1184–1202 `filterPatterns`/enzyme filtering. The second half of this finding (freebarcodes inertness caveats) was **not** added — already present at §3(e), §6 and §8. |
| F-B5 | Stochastic output with no seed argument | One bullet added to §6 "Structural, verified" | No `set.seed` anywhere in `R/processVCFfast.R`; no seed parameter in the L1099–L1115 signature; sampling at L73, L77, L156, L394, L519, L531, L609, L629, L739, L751, L830, L941, L953 and the L1213 `runif` shuffle. |

---

## REJECTED

| # | Finding | Why not applied | Evidence |
|---|---|---|---|
| C-11 | "misquoted" — the `subseq(genome[[chr]], POS - upstream, POS + downstream)` expression is not a verbatim code quotation | The audit concedes the expression is accurate pseudocode, and the substantive claim (flat genomic window, no splice or frame awareness) verifies. The expression is self-evidently schematic — `chr` is not a variable anywhere in the source — so it does not read as a quotation. Rewriting it would create fresh unaudited text for zero factual gain. | Actual code at L337/L555/L777 is `subseq(genome[[paste0('chr', as.character(snp$CHROM))]], start = rangestart, end = rangeend)` with `rangestart`/`rangeend` set at L329–L335 and swapped when `snp$reverseGene`. The window is flat; nothing splice- or frame-aware exists. |
| C-32 | "minor-discrepancy" — README's `"MPRAREV"` rendered with single quotes in the record | Purely typographic, and an artifact of nesting: the quotation already sits inside a double-quoted italic span, so inner single quotes avoid a nested-quote collision. The words are verbatim. Changing it would degrade the markup for no factual gain. | `README.md` L58 and `ui.R` L125 both use `"MPRAREV"`; the surrounding sentence in the record is otherwise word-for-word identical to source. |
| F-B3 | "incomplete" — `alter_aberrant`'s one-site, cross-variant and element-junction limits could be missed | The omission does not exist. Every fact the finding lists is already in the record, just distributed rather than collected in one place. | Record §3 `consequence_annotation` quotes `'Failed - More than one aberrant digestion site in context'` and `'Failed - aberrant site generated at interface between sequence elements'`; §3 `synthesis_constraints` covers `check_cross_center()` (now with the corrected caveat from C-15) and tags `alter_aberrant` "under development"; §6 repeats the under-development tag. |

---

## ESCALATED — human decision required, record left as-is at these points

**1. `maintained = no` may now rest on a weaker footing than the record states.**
The correction at C-19 was factual and has been applied, but it introduces information the record
did not previously contain. Verified: `mpradesigntools` has five issues (#1–#5) filed between
2018-08-30 and 2024-10-23; all are closed; the author replied on every thread, closed four of
them himself, and his most recent reply is **2024-10-27** — about 22 months before this record's
as-of date, not 5 years 8 months. In #1 he shipped a fix the same day and pointed the reporter at
a newly added `alter_aberrant`; in #2 he debugged a user's VCF over five weeks; in #3–#5 he
diagnosed failures from stack traces.
The code is genuinely frozen at 2020-12-17 and everything else in the evidence string (0 releases,
dead hosted app, no CRAN/Bioconductor, no descendant) is verified and unchanged.
**Question:** does "maintained" mean commits, or maintainer responsiveness? Under the first the
value stands; under the second it is contestable. Options: (a) keep `no` and let the corrected
evidence stand as a documented tension; (b) keep `no` but add an explicit "code frozen, author
still responsive through 2024" gloss; (c) revisit the value. **I did not change the value.**

**2. `interface`: the correction at C-17 partly rehabilitates a sub-claim a previous pass withdrew.**
§0 row 2 and §8 both record that "the GUI is currently non-functional" was judged "overstated and
rebuttable — only the *hosted* instance is 404", and that the sub-claim was withdrawn. The verified
fact is narrower than the withdrawn claim but points back toward it: `outputs/inertTwelveMersChar.RData`
is gitignored and absent, so *Create sequences* cannot run from the distributed checkout at all,
independent of the missing images. I corrected the specific false sentence in §3 `interface` and
**left §0 and §8 untouched**, per the instruction to leave nearby unfixed sentences alone.
**This is a live inconsistency in the record right now:** §3 says the sequence generator cannot run
from source; §0 and §8 say the non-functionality claim was withdrawn as overstated.
**Question:** how should §0's change-history row and §8's "reviewer right — sub-claim withdrawn"
line be reconciled? Note `interface = yes` does not appear to be at risk either way — it rests on
the package API plus a GUI, and the package API is intact. **I did not change the value.**

**3. `license = yes` and the two-artifact scope (F-B7).**
Verified: `mpradesigntools/DESCRIPTION` L8 declares `License: GPL-2`; `designMPRA` has **no**
license file (`git ls-files` shows none) and GitHub returns `license: null` for both repos. The
record already states the API result and the absence of a `LICENSE` blob in either tree, but never
draws the conclusion the audit asks for: GPL-2 covers the package only, and the Shiny app —
half of what this record defines as "the tool" — carries no declared license.
**Question:** stating that plainly is a one-sentence addition, but it changes what `license = yes`
rests on for a subject defined as both artifacts. Add the sentence and keep `yes`, or send the
value to the reconciliation pass? **I did not change the value or add the sentence.**

**4. Three further example-only barcode scripts (F-B6).**
Verified present: `designMPRA/scripts/dnaBarcodeTest.R`, `dnaBarcodeTimingTest.R`,
`dnaBarcodeTwelverMers.R` — third-party `DNABarcodes` experiments creating a 10-mer set,
benchmarking Conway-heuristic sets of length 4–12, and generating an Ashlock-heuristic 12-mer set.
All use hard-coded local paths and none is reachable from either public API.
**Question:** whether omitting them is material is a framing call, and adding them cuts both ways —
they could read as broader barcode-design capability (the audit's own stated worry) even though
they are developer scratch work adjacent to `barcode_generation = yes`. **Not added.**

**5. True-but-uncited external facts (C-29, C-33, C-34).**
I verified all three sets of claims independently and found them **correct**:
`malacoda` is Bayesian MPRA analysis published in PLOS Computational Biology, 2020-07-21,
DOI `10.1371/journal.pcbi.1007504` (resolves); `biocLite.R` now returns 404;
`BSgenome.Hsapiens.UCSC.hg38` 1.4.5 is 731,245,539 bytes, so "~700 MB" is fair; and my own search
turned up no alternative hosted designMPRA instance. What the audit objects to is that none of this
is sourced in the record — the cited `mpradesigntools` README supports only the weaker "downstream
QC and statistical analysis" for malacoda, and no query is recorded for the negative web claim.
**Question:** whether this record's standard requires citations for external, version-sensitive
facts is a documentation-policy call. The facts are sound; only their provenance is thin.
**No citations added.**

**6. "Likely still runnable on a pinned 2020-era toolchain" (C-30).**
Confirmed unsupported: no install was attempted (correctly, per instructions), and there is no
lockfile, log or execution result behind it. The record does already label the whole row
"Non-trivial, not verified", so the speculation is disclosed rather than hidden.
**Question:** keep the disclosed prediction, soften it, or drop it? Deciding would be re-deciding an
author's hedge. **Left as written.**

**7. Fact-check section C (balance).**
Explicitly reserved: judgments about emphasis and proportion are for the authors. For what it is
worth, the four strength-overstatements section C names (reproducible barcode generation,
universally distinct repairs, local Shiny usability, and the contrary blanket denial of barcode
homopolymer checking) have all now been corrected above, which is the factual half of its
complaint. The proportion question is untouched.

---

## Known local roughness left in place, by design

- **§0 change table and §8** still describe the `interface` sub-claim as withdrawn (escalation 2),
  and §0 still lists the `barcode_generation` "3–17 bp" fix without the new KpnI finding. §0 is a
  history of a *previous* pass; rewriting it would be re-narrating work I did not do.
- **§6 "Authors' own"** still says "≤ 20,000" while §3 now says "just below 20,000". Both are
  correct: §6 reports the authors' prose, §3 reports the implemented predicate. The app itself is
  inconsistent.
- **§3(a)** still calls the shipped set "1,140,292 inert twelvemers" (the wording of
  `man/twelvemers.Rd`'s own title) immediately before the new sentence noting 3,767 of them carry a
  KpnI site.
- **§8 "Highest confidence"** still lists `barcode_generation` and `maintained` among the cells read
  "exhaustively" by two independent readers. Several corrections above land in exactly those two
  cells. Not rewritten — that paragraph is a prior pass's self-assessment.

## Files touched

| File | Action |
|---|---|
| `revision/tool_survey/final/mpradesign.md` | Edited in place — 34 findings applied across 35 surgical edits (several findings recur in both §3 and §5) |
| `revision/tool_survey/fixes/mpradesign.md` | Created (this change log) |

## Pass 2 — policy application

**Outcome count (seven open items): 4 applied · 2 declined-by-policy ·
0 rejected-unverifiable · 1 escalated.** No capability value was changed.

### 1. `maintained = no` value basis — applied

- **Edit:** defined `maintained` as ongoing code/release maintenance rather than
  occasional issue responses, added the maintainer's decisive 2024-10-27 statement,
  and added `/issues/5/comments` to the source line.
- **Verification:** the current GitHub API reports five closed issues and the author
  commented in all five. In issue #5 he wrote, *"I should probably archive this
  repository soon as I don't have adequate bandwidth to maintain it at this point
  unfortunately."* The API and current clone independently give HEAD
  `afd386ef12051bb0a37ad63a6f92acd555246584`, commit time
  `2020-12-17T19:27:27Z`; the API gives zero tags and zero releases.
- **Decision:** the corrected evidence still supports `no`. The recent reply is
  retained as counterevidence rather than hidden, but the author's own statement
  resolves the supposed responsiveness/maintenance tension.

### 2. `interface` change-history inconsistency — applied

- **Edit:** changed only the stale wording in §0, the caveat lead-in in §3, and the
  §8 summary. They now distinguish a launchable UI/power tab from a sequence-generation
  path that cannot run from the distributed checkout. `interface = yes` is unchanged.
- **Verification:** `designMPRA/scripts/processVCFfast.R` loads
  `outputs/inertTwelveMersChar.RData`; `.gitignore` excludes `outputs/*`; `git ls-files`
  contains no such object; and `scripts/shinyAppsDeploy.R` lists it as a separately
  supplied deployment file. The published endpoint again returned HTTP 404 on
  2026-08-14. This bears on Table 1's Availability cell and is now stated consistently.

### 3. `license = yes` across two artifacts — escalated

- **Edit:** none; the value and prose were left unchanged.
- **Verification:** `mpradesigntools/DESCRIPTION` declares `License: GPL-2`.
  Complete tracked-file lists for both repositories contain no `LICENSE` or `COPYING`,
  `designMPRA` contains no other license declaration, and GitHub's repository API
  returns `license: null` for both.
- **Reason:** the record defines its subject as the package plus the Shiny app. The
  package declaration does not support an unqualified `yes` for both artifacts.
  Per the value-basis ruling, this remains an escalation and the value was not changed.

### 4. Three example-only `DNABarcodes` scripts (F-B6) — declined-by-policy

- **Edit:** none.
- **Verification:** all three scripts were read in full at the current `designMPRA`
  head. `dnaBarcodeTest.R` creates a 10-mer set and reads a hard-coded laboratory CSV;
  `dnaBarcodeTimingTest.R` benchmarks Conway sets of length 4–12 and writes under
  `outputs/`; `dnaBarcodeTwelverMers.R` creates an Ashlock 12-mer set and writes under
  `outputs/`. None is called by `ui.R`, `server.R`, the historical worker, or either
  exported package function.
- **Policy-B test:** these scratch experiments do not change any of the 17 locked-row
  scores. In particular, they generate barcode candidates rather than modify a designed
  variant, so they do not satisfy `Constraint-repaired variants`. They also do not change
  Table 1's Key-features wording, which already says barcode **assignment** with
  error-correcting sets.

### 5. External sourcing and provenance (C-29, C-33, C-34) — applied

- **C-29:** dropped *"No alternative hosted instance found"*. The specific published
  endpoint's 404 is reproducible; the web-wide negative was not.
- **C-33:** replaced the uncited `malacoda` publication/algorithm parenthetical with
  the narrower primary-source statement supported by `mpradesigntools/README.md`:
  QC and statistical analysis.
- **C-34:** replaced the uncited date/size assertions with plain verified statements:
  the README's `biocLite.R` URL returns 404; current Bioconductor instructions use
  `BiocManager::install()`; active source calls `unnest_legacy()`, `purrrlyr::by_row()`
  and `write_tsv(path=)`; and the hg38 BSgenome is imported. Added links to the official
  Bioconductor, tidyr and readr documentation. The audit premise that active code uses
  `tibble::data_frame` did **not** verify: HEAD uses `tibble()`; only an unused NAMESPACE
  import and commented developer-script lines remain, so that claim was removed.

### 6. Pinned-2020-toolchain hedge (C-30) — applied

- **Edit:** dropped *"Likely still runnable on a pinned 2020-era toolchain."* The
  surrounding `Non-trivial, not verified` assessment remains.
- **Verification:** neither tracked tree contains a lockfile, environment manifest,
  installation log or test log, and no installation was performed in either pass.
  The surviving risk statements are direct repo/docs facts with primary sources.

### 7. Balance and emphasis (fact-check §C) — declined-by-policy

- **Edit:** none. Policy A declines all rebalancing and proportional-emphasis requests
  for these working records.
- **Verification:** no factual residue from §C was used as a reason to rebalance; its
  four concrete factual points had already been independently adjudicated in Pass 1.

### Version-drift check

No edit. The header governs v0.2.0, and current `DESCRIPTION` at the current remote HEAD
still says `Version: 0.2.0`; GitHub has no tags or releases. There is therefore no later
release requiring a Policy-D parenthetical. The record now also defines commit time versus
GitHub `pushed_at` explicitly and consistently.

### Row-substitution candidate — report only

**Provisional candidate:** replace **`Recombination and chimeras`** with **`MPRA
statistical-power / assay-sizing analysis`** if cross-tool scoring confirms the locked
file's flagged risk that recombination is a PoolParty-only row.

- **Ground:** honesty. This would be a row where PoolParty is not best, replacing one
  where it is expected to be best.
- **Evidence:** `designMPRA/server.R` computes a Bonferroni-adjusted `pwr.t.test` curve
  from barcode replicates, transfection replicates, number of variants, alpha and activity
  SD; `ui.R` exposes those inputs; the paper abstract and §2 make power analysis a headline
  capability. A repo-wide search of the current PoolParty package, README and tests found
  no power-analysis or assay-sizing API. PoolParty's documented `recombine` operation,
  by contrast, constructs breakpoint-defined chimeras. No substitution was made; the
  other six tool cells would still require primary-source scoring.
