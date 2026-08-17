# Adversarial review — `mpradesign` (MPRA Design Tools, Ghazi et al. 2018)

Reviewed 2026-08-10. All source re-fetched independently from GitHub raw + the GitHub REST API;
paper/supplement re-extracted from the PDF with PyMuPDF. No install attempted (per instructions).

## What I independently verified

| Check | Result |
|---|---|
| `DESCRIPTION` | v0.2.0, `License: GPL-2`, Imports = dplyr/purrr/readr/tidyr/tibble/stringr/purrrlyr/magrittr/Biostrings/BSgenome.Hsapiens.UCSC.hg38. **Confirmed.** |
| `NAMESPACE` | exactly `export(processVCF)` + `export(spread_and_fix_indels)`; only nucleotide-level Biostrings imports. **Confirmed.** |
| `R/processVCFfast.R` | 1459 lines, 16 top-level functions, 2 exported. Read the barcode-filter block, `processVCF` body, `processSnp`, `randomly_fix`. **Confirmed.** |
| `man/processVCF.Rd` | full signature and the "under development" tags on `alter_aberrant`/`extra_elements`/`max_construct_size`; "freebarcodes sets only meet the traditional MPRA barcode requirements to varying degree". **Confirmed.** |
| README | biocLite install path, 5+0.01·n runtime formula, barcode table, Planned Features (mm10, parallelization, Sharpr-MPRA all unchecked). **Confirmed.** |
| `designMPRA/ui.R`, `server.R` | 20 000-sequence cap, 12 bp barcode lock, "mainly for simple demonstrative purposes", accuracy disclaimer, static `MPRAseqDiagram.png`. **Confirmed.** |
| GitHub API | `mpradesigntools` pushed_at 2020-12-17T19:27:54Z, 3★, 1 fork, 0 issues, unarchived, **0 tags/releases, 1 branch**; `designMPRA` pushed_at 2017-09-26, 1★. `license: null` on both. **Confirmed.** |
| `https://andrewghazi.shinyapps.io/designmpra/` | HTTP **404**, 3407-byte shinyapps error page. Case variant `designMPRA/` also 404. shinyapps.io account root 404. **Confirmed.** |
| CRAN | `cran.r-project.org/package=mpradesigntools` → **404**. r-universe search → 0 results. **Confirmed.** |
| Bioconductor | `bioconductor.org/packages/mpradesigntools/` returns HTTP 200 but the page is Bioconductor's **"Removed Packages"** fallback, i.e. not indexed. **Confirmed not on Bioconductor.** |
| Paper Supplement S1.3 (PDF p.4), S5 (pp.25-26) | primer/enzyme defaults and the 12-mer screening criteria quoted verbatim. **Confirmed.** |
| grep for `motif|spacer|permut|combinator` across README, ui.R, server.R, processVCFfast.R, processVCF.Rd, and the whole paper+supplement | **0 hits everywhere.** |

Net: the extraction is unusually accurate. Every "no" I could attack held up. Two cells and a handful of
evidence strings need fixing.

---

## Cell-by-cell

### Cells that survive attack unchanged

`dag_chaining`, `lazy_evaluation`, `combinatorial_motif_place`, `assay_dms`, `assay_insilico`,
`genome_coordinates`, `transcript_models`, `exon_intron_split_codons`, `hgvs_input`, `vcf_vep_output`,
`consequence_annotation`, `primer_design`, `codon_optimization`, `assay_mpra`, `license`,
`barcode_generation`, `design_visualization`, `synthesis_constraints`, `library_as_object`,
`per_sequence_provenance`, `automatic_naming`, `mixed_mutagenesis_one_pool`.

Independent corroboration worth adding for `dag_chaining = no`: the **only** fork of the package
(`mmandelia/mpradesigntools`, 2 commits, 2023-03-24) had to **edit `R/processVCFfast.R` in place** to swap
KpnI/XbaI/SfiI for AgeI/SbfI/I-SceI and to drop the barcode from the construct. That is the strongest possible
demonstration that the construct architecture is not composable or parameterisable — the only way to change it
is to fork the source. Recommend citing this in the referee response.

### 1. `maintained = "yes"` — **WRONG**

This is the one cell I would change outright. The value contradicts its own evidence string
("Effectively unmaintained for ~5.5 years as of 2026-08-10; never published to CRAN or Bioconductor").
Independently: last commit 2020-12-17, **zero tags, zero releases, one branch**, 0 open issues (because nobody
files any), 3 stars, not on CRAN, not on Bioconductor, not on r-universe, and the only fork is a private
2023 hack. Nothing has happened in 5 years 8 months.

**Suggested value: `no`** (or `partial`/`dormant` if the matrix has such a value — repos are public and
unarchived, so "source still retrievable" is true; "maintained" is not).

Risk of leaving it: a matrix row reading "maintained: yes" for a 5.7-year-dormant package is the kind of
internal inconsistency a referee will use to argue the whole survey was done carelessly. And it is
inconsistent in the *generous* direction, so fixing it does not create an attack surface.

### 2. `interface = "yes"` — value right, one sub-claim **OVERSTATED**

The row asserts "**The GUI is currently non-functional** (see availability)." That is too strong and is
exactly the sentence Ghazi would rebut.

What is actually true:
- The *hosted* instance at shinyapps.io is 404 (verified).
- The app itself is **still runnable locally**: `designMPRA` is a standard two-file Shiny app
  (`ui.R` + `server.R` at repo root) that a user can clone and launch with `shiny::runApp()`.
- Two real caveats that are *better* attacks than "non-functional", because they are verifiable:
  1. **`www/*` is in `.gitignore`.** The repo contains no `www/` directory, so `MPRAdiagram.png`,
     `MPRAseqDiagram.png` and `formulas.png` — every image the UI references, including the construct-layout
     diagram — are **absent from the distributed source**. A local run renders broken images.
  2. **The app is a divergent 2017 code path.** `server.R` does `source('scripts/processVCFfast.R')`
     (18.6 kB) — *not* the package's `R/processVCFfast.R` (62.4 kB) — and calls
     `processVCF(vcf, nper, seqwidth=, fwprimer, revprimer, enzyme1..3, updateProgress)`. That signature has
     **no** `barcode_set`, `alter_aberrant`, `max_construct_size`, `flip_RV`, `filterPatterns`,
     `ensure_all_4_nuc`, and only a **symmetric** `seqwidth` instead of separate up/downstream ranges. So the
     GUI cannot reach any 0.2.0 feature, and the two artefacts are not two front-ends on one engine.

**Suggested rewording:** keep `interface = yes`; replace "currently non-functional" with "the hosted instance
is dead (404); the app source is still runnable locally but ships without its `www/` image assets and wraps a
2017 snapshot of the engine that predates every v0.2.0 option."

### 3. Evidence-string errors (values unaffected, but fix before quoting)

- **`barcode_generation`**: "named `barcodes<length>-<n_correctable_errors>` (**10-17 bp**, 1-2 indel errors)".
  Wrong — the shipped sets run **3-17 bp** (`barcodes3-1` … `barcodes17-2`; verified in `available_sets`,
  `data/`, `man/`, and the README table). The extraction's own `additional_capabilities` bullet says 3-17,
  so this is an internal inconsistency. Also note only `barcodes17-2` exists at 17 bp (no `17-1`).
- **`library_as_object`**: "NAMESPACE exports only two functions **and 28 datasets**". NAMESPACE exports
  **no** datasets — there are zero `export()`/`exportData` lines for them; they are available via
  `LazyData: true`. Say "two exported functions; 28 lazy-loaded barcode datasets".
- **`design_visualization`**: "server.R has **exactly six** outputs". There are **eight** `output$`
  assignments (`powerPlot`, `testHead`, `timeText`, `warnText`, `conditionalDownload`, `failed`, `inputHead`,
  `downloadSequences`). The two omitted are a leftover `head(mtcars)` debug stub and the TSV download handler —
  neither is a visualization, so **the `no` stands**, but do not quote a falsifiable count.
- **Supplement page cites**: S1.3 is on PDF **p.4**, not p.5.
- **Activity SD**: the extraction says "mean 0.926". The *paper* says "**on average around .95**"; .926 is the
  Tewhey-only figure (supplement line 440), .952 is the Ulirsch figure. Quote ".95 across both studies" or
  attribute .926 to Tewhey.

### 4. Findings that make two "partial" cells *safer*, not weaker

`per_sequence_provenance = partial` is if anything **generous**, which is the right direction. Two things the
extractor missed that support it:
- The final `select(ID, type, allele, snpIndex, totIndex, barcode, sequence, site_fix_info, notes)`
  **drops `CHROM` and `POS`** (they exist in the intermediate tibble and are carried in the `failed` table, but
  not in `result`). The delivered library therefore contains **no genomic coordinates** — provenance back to
  the locus depends entirely on the VCF `ID` column being populated.
- **`notes` is variant-level, not sequence-level**: `res$notes = notes` assigns one string to all `2*nper`
  rows of a SNP. Only `site_fix_info` genuinely varies per sequence.

`automatic_naming = partial` similarly errs generous: `ID` is verbatim carried from the input VCF's ID column,
so a VCF without rsIDs yields a library whose rows are identified only by integer indices.

---

## Capabilities the extractor missed

1. **Library-wide barcode uniqueness by construction.** `split(shuffled_mers, ceiling(1:length(shuffled_mers)
   %% nrow(vcf)))` gives every variant a **disjoint** barcode pool, then `sample(..., replace = FALSE)` draws
   within it — so no barcode is ever reused anywhere in the library and the barcode→variant map is injective by
   construction, not by post-hoc dedup. This is a distinct capability from "generating a barcode set" and is
   the kind of thing a referee would say the matrix ignored.
2. **Per-replicate randomised constraint repair.** `randomly_fix()` does not apply one fix to a variant — it
   samples `nrow(res_df)/2` **different** single-base alterations of the offending site and gives a different
   one to each barcode replicate (`rep(altered_patterns$altered_pattern, times = 2)`), while giving **matched**
   fixes to the ref and alt constructs so the allelic contrast is not confounded by the repair. It also refuses
   to touch N positions of an ambiguous site (`nonrandomly_replace_N`) and hard-fails if the variant itself
   falls inside the aberrant site (`check_cross_center`). This is materially more sophisticated than "repairs
   restriction sites".
3. **Bounded retry budget.** Up to **40** barcode-resampling attempts per variant before emitting
   `'Failed - SNP sequence could not be generated without an aberrant digestion site even after barcode
   resampling'` (three copies of this loop, one per variant class).
4. **Multi-allelic expansion duplicates the reference arm.** `spreadAllelesAcrossRows()` emits one row per ALT,
   and each row independently generates `nper` ref + `nper` alt constructs — so a tri-allelic site contributes
   two full ref sets. Affects library-size accounting; not documented anywhere by the authors.
5. **Documented barcode-quality audit of the borrowed sets.** `R/import_freebarcodes.R` contains (commented-out
   but reproducible) code computing, per freebarcodes set, the fraction with all four nucleotides, with
   homopolymer runs ≥4, and with human miR seeds, concluding "none of the new sets have nucleotide runs, and
   the larger ones have all 4 nucleotides ~85%-98% of the time". This is the basis of the man page's
   "varying degree" caveat and is worth citing if the response characterises their barcode quality.
6. **App/package divergence** (see §2 above) — arguably its own row-level fact: the two published artefacts are
   two different programs, not one engine with two front-ends.
7. **Resolution of the extractor's open question #6.** The only fork of `mpradesigntools` is
   `mmandelia/mpradesigntools` (2 commits, 2023-03-24), a private modification swapping the enzyme trio and
   *removing* the barcode from the construct. It is not a successor, not published, and not referenced by the
   authors. `goldenac/MPRADesignGenerator` is independently confirmed **not a fork** (`fork: false`, no shared
   history, no citation link). No successor or maintained descendant exists.
8. **Never-on-Bioconductor is now positively verified** (the package URL resolves to Bioconductor's "Removed
   Packages" fallback page), which the extraction had only inferred from rdrr.io.

---

## Things I tried to falsify and could not

- No hidden CLI: no `inst/`, no `exec/`, no `Rscript` entry point in the file tree.
- No plugin/extension mechanism: no S4 generics, no `methods` import, no `...` pass-through to user hooks.
- No FASTA/GenBank writer: only `readr::write_tsv` and `save()`.
- No motif/element/spacer/permutation vocabulary anywhere (0 grep hits across code, docs, and the 28-page PDF).
- No codon/amino-acid/reading-frame concept; Biostrings imports are strictly nucleotide-level.
- No GC-content, Tm, secondary-structure, or repeat check on either the barcode or the construct.
- No mm10/other genome: `BSgenome.Hsapiens.UCSC.hg38` is a hard `import()` in NAMESPACE and is hard-coded as
  `genome = BSgenome.Hsapiens.UCSC.hg38` inside `processSnp` — not a parameter.
- No alternative hosted instance of the Shiny app on shinyapps.io, posit.cloud, or elsewhere that I could find.

## Bottom line

The extraction is defensible as written except for `maintained`, which should flip to `no`, and the
"GUI is currently non-functional" phrasing in `interface`, which should be narrowed to the hosted instance.
Fix the 10-17 bp → 3-17 bp barcode range, the "28 datasets in NAMESPACE" claim, and the "exactly six outputs"
count before any of those strings are quoted in the referee response. The four flagged "partial" cells all err
in the competitor's favour, which is the correct direction for this document.
