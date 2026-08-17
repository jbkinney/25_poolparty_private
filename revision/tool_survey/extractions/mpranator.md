# MPRAnator — evidence memo

**Tool:** MPRAnator (MPRA Motifs, MPRA SNPs, Transmutation, PWM Seq-Gen)
**Authors:** Ilias Georgakopoulos-Soares, Naman Jain, Jesse M. Gray, Martin Hemberg (Wellcome Trust Sanger Institute)
**Paper:** *Bioinformatics* 33(1):137–138, 2017 (Advance Access 6 Sep 2016). doi:10.1093/bioinformatics/btw584, PMID 27605100, PMCID PMC5198521
**Survey date:** 2026-08-10

---

## 1. Sources consulted

| Kind | Reference |
|---|---|
| PDF | `revision/tool_survey/papers/Georgakopoulos-Soares2017_MPRAnator.pdf` (2 pp., Application Note) |
| PDF (supp.) | `btw584_supp.docx`, retrieved via Europe PMC `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5198521/supplementaryFiles` (OUP/silverchair direct link is 403; PMC `/bin/` link returned HTML, not the docx) |
| Prior analysis | `revision/tool_survey/prior_analyses/Georgakopoulos-Soares2017_MPRAnator_analysis.md` |
| Repo | https://github.com/hemberg-lab/MPRAnator (+ GitHub REST API for commits/contents/license) |
| Repo source read | `part1.py`, `oligo.py`, `parseSNPs.py`, `part3.py`, `mycustom.py`, `myfunctions.py`, `iliasApp/forms.py`, `iliasApp/views.py`, `iliasApp/urls.py`, `iliasApp/templates/iliasApp/docs.html`, `downloadables/MpraMotifs_script.py`, `downloadables/MpraSnps_script.py` |
| Web service | https://www.genomegeek.com/ and its pages `/MPRA/Motifs/`, `/MPRA/SNPs/`, `/MPRA/Transmutation/`, `/MPRA/PWM/`, `/MPRA/documentation/` |
| Web (mirror) | https://www.sanger.ac.uk/tool/mpranator/ (the paper's URL `www.sanger.ac.uk/science/tools/mpranator` 302-redirects here) |
| PyPI | `https://pypi.org/pypi/mpranator/json` → **404**; no package |
| GitHub search | `q=MPRAnator` → exactly 1 repository (`hemberg-lab/MPRAnator`); no newer fork/successor |

---

## 2. Availability status (checked 2026-08-10)

### Repository — alive but frozen and skeletal
- **Exactly one commit:** `2015-12-27T14:57:53Z`, author Naman Jain, message `"first commit"`. (GitHub API `/commits`.)
- GitHub API metadata: `created_at` 2015-12-27, `pushed_at` 2015-12-27, `updated_at` 2016-01-20, size 79 KB, 1 star, 1 fork, 0 open issues, `archived: false`.
- **`README.md` is 12 bytes and contains exactly one line: `# MpraNator`.** No install instructions, no usage docs, no dependency list, no example commands.
- **No `LICENSE` file; GitHub API reports `license: null`** — despite the paper's claim of MIT (see §3).
- Code is **Python 2** (`print "…"` statements, `xrange`) on **Django 1.x** (`django.conf.urls.patterns`, `staticfiles` template tag). It would not run on a modern Python/Django without porting.
- The committed snapshot is a **development snapshot**: live `pdb.set_trace()` calls remain in `parseSNPs.py` (e.g. lines 406, 443), and `oligo.py` carries a hard-coded interpreter path `#!/Users/Naman/miniconda/bin/python`.
- **The repo is older than the published tool.** It contains only 3 of the 4 published modules — there is no PWM Seq-Gen code, and `Part3Form` in the repo lacks the reverse/complement options that the live site and the published Supplementary describe. So the code actually deployed at genomegeek.com is **not** what is in the public repo.

### Web service — alive, but the flagship module is broken
- `https://www.genomegeek.com/` → HTTP 200. Homepage text: *"We have 4 tools: MPRA (Motifs) construction … MPRA (SNPs) construction … Transmutation … PWM-SeqGen …"*. Contact addresses (Naman Jain, Ilias Georgakopoulos-Soares) still listed. Note: the URL paths changed since the repo (now `/MPRA/Motifs/`, `/MPRA/SNPs/`, `/MPRA/Transmutation/`, `/MPRA/PWM/`, `/MPRA/documentation/`); the repo's paths 404.
- **Transmutation — WORKS.** POST to `/MPRA/TransmutationResults/` with the site's own sample data returned HTTP 200 and real output:
  `>sequence1| Mutated_nucleotides - 3 | Scrambled - Yes | Reversed - No | Complemented - No` / `CGGTAAGGAGAAGATAGTTGCAGTGTCC`
- **MPRA SNPs — WORKS.** POST to `/MPRA/SNPs/Results/` with the site's sample data returned HTTP 200 and real output, including a WT/reference entry:
  ```
  >1 42990 43041| SNP rs43434 | Position 42991 | Nucleotide change T:G   AGAGATTAGTCAGT...
  >1 42990 43041| REFERENCE                                              ACAGATTAGTCAGT...
  >1 42990 43041| SNP rs121212 | Position 43030 | Nucleotide change A:G  ACAGATTAGTCAGT...
  ```
- **PWM Seq-Gen — reachable and validating** (returned the form with `"This field is required"` for the threshold field on my under-specified POST, i.e. the view is live; I did not push a fully valid submission through).
- **MPRA Motifs — BROKEN (HTTP 500).** Every *valid* submission to `/MPRA/Motifs/Results/` returned `Server Error (500)`. Tested: (a) the site's own "Load Sample Data" values verbatim; (b) with and without barcode parameters; (c) four different `ordering` strings; (d) a single-motif minimal case. Form validation itself still works — an intentionally invalid submission (empty motifs) returned HTTP 200 with the inline error `"Enter motifs!"` — so the 500 originates in the design/generation code, not in request parsing. **This is the module the paper leads with.** (Single observation window from one host; a referee could conceivably see different behaviour, so phrase any manuscript claim as "observed to return HTTP 500 on 2026-08-10".)
- `https://www.sanger.ac.uk/tool/mpranator/` → HTTP 200, still describes the tool and points to `http://genomegeek.com/`. It is a landing page only; the tool is not hosted at Sanger.

### Installability
There is **no installable artifact**: no PyPI/conda/CRAN package, no Dockerfile, no `requirements.txt`/`setup.py` in the repo, no install section in the README. Programmatic use is HTTP POST against genomegeek.com (see §4, `interface`). The tool is therefore usable today only as a hosted web service, and only for 3 of its 4 modules.

---

## 3. License discrepancy (worth flagging carefully)

Paper, Availability section (verbatim): *"The source code is available on www.github.com/hemberg-lab/MPRAnator/ under the MIT license."*

Repository reality: no `LICENSE`/`COPYING` file anywhere in the tree (`GET /repos/hemberg-lab/MPRAnator/contents/` and the recursive tree listing both show none), and the GitHub API reports `"license": null`. Record the row as "MIT per the paper; no license file in the repository."

---

## 4. Capability-by-capability evidence

### Block A — library specification

**`library_as_object` — partial.**
One form submission (or one HTTP POST) yields an entire multi-sequence library as a FASTA stream; the user never loops over sequences. But there is no library object the user can hold, inspect, subset, or pass on. The only user-facing representation is a text blob. Internally `mycustom.py` defines `FastaFile`/`Sequence` classes, but these are input parsers inside a Django app, not an exposed API. Documented output (docs page): *"The result page displays the synthesized oligonucleotides in the FASTA format. The user is able to download ( plain text file ) the generated oligonucleotides."*

**`dag_chaining` — partial.**
Manual, documented hand-off only. Docs page, PWM-SeqGen section: *"The output of PWM Seq-Gen can be inputted into MPRAs Motifs tool, therefore allowing for the design of MPRA experiments using PWMs."* Supplementary §2.3: the Transmutation tool serves *"to generate negative controls for the MPRA SNP and Motif design tools."* But there are four separate web forms with four separate result pages; chaining means copy-pasting FASTA between them. No pipeline object, no graph, no nesting, no composition primitive anywhere in the source (`iliasApp/urls.py` exposes independent views). The one intra-tool "composition" is the drag-and-drop *ordering* of sub-components (adapter / restriction site / background / barcode) — that is linear concatenation, not a DAG.

**`lazy_evaluation` — no.**
Everything is materialized eagerly. `oligo.py::oligo()` builds `allResults = []` by enumerating `itertools.product(range(distanceFromLeftEdge, len(BackgroundS)), repeat=len(motifsL))` and appending every passing sequence; `part1.py::createMPRAResultOutput()` concatenates the whole library into a single `response` string that is then rendered/downloaded. Consequently the tool must impose hard size caps: `forms.py` `maximumNumberOfMotifsTimesSequences = 100` with the error *"Only a maximum of %s motifs^4*sequences allowed"*, `sequenceLimit` of 1000 (motifs) / 50000 (SNPs), and `sequenceLengthLimit = 20000`.

**`mixed_mutagenesis_one_pool` — partial.**
*Within* a module, heterogeneous designs do coexist. Supplementary §2.1: *"If multiple motifs are inputted, sequences will be generated for each motif separately and for the motif combinations. For instance, if three motifs are inputted, output sequences for single, pairs and triplets of motifs will be generated."* (`part1.py::generateCombinations` + `generatePermutations`.) The SNP module emits exhaustive single-SNP variants, all SNP combinations, **and** a `REFERENCE` (WT) sequence in one output — confirmed in the live run above — plus per-sequence barcode replicates (`numOfBarCodesPerSequence`, max 5). *Across* modules it is not possible: motif placement, SNP substitution and random mutagenesis live in three separate tools with three separate submissions and three separate output files; there is no single specification that mixes them. Sampled/random mutagenesis (Transmutation) can never coexist in the same specification as motif or SNP design.

**`combinatorial_motif_place` — yes.**
Supplementary Fig. 1 caption: *"Up to 4 motifs are substituted into each of the background sequences using all possible permutations. Distance from the edges and minimum and maximum spacing between motifs restrict the positioning of motifs in the sequences. Interval of substitution … determines the frequency of motif substitutions."* Code: `part1.py::generateCombinations` (all subsets), `generatePermutations` (all orderings within a subset), `equaliseSize(..., 4)` (*"Adds empty strings instead of motifs, if less than 4 motifs are used"*), `oligo.py::oligo()` enumerating position tuples with `itertools.product`, filtered by `testIfWithin` (min/max spacing) and `testIfFarFromRightEdge`. Orientation is limited: there is a single global `reverseComplement` checkbox (*"Reverse complement sequence before motif substitution?"*) applied to the scaffold, not per-motif orientation enumeration.

**`barcode_generation` — yes.**
Form fields (`forms.py`, and live on both `/MPRA/Motifs/` and `/MPRA/SNPs/`): `barCodeLength` (min 10), `minimumGCContent`, `maximumGCContent`, `barCodeDistance`, `numOfBarCodesPerSequence` (1–5). Docs tooltips: *"Barcode edit distance: The Levenshtein distance between each barcode. The default is 2."*; *"Number of barcodes per sequence: This specifies the number of barcodes inserted per sequence (replicates)."* Implementation `part1.py::barcode_generator(number, mingc, maxgc, barCodeLength, diffs=2)` — *"generates barcodes that differ from existing barcodes by a minimum number of bases, each sequence has a unique barcode"*. Barcodes are attached by `createMPRAResultOutput()` at the user-chosen position in the ordering. **Caveat for precision:** the implementation counts position-wise identities (`if short_term[pos] == bc[pos]: count += 1`), i.e. a Hamming-style criterion, even though the docs and Supplementary call it Levenshtein. The row is still "yes"; the caveat is an accuracy note, not a capability denial.

**`per_sequence_provenance` — partial.**
Every output sequence carries a structured, machine-parseable header. Docs: *"A header is composed of one or more DESCRIPTORs and each DESCRIPTOR is composed of a LABEL and INFO. The descriptors are delimited by a |, i.e. a 'pipe'."* Documented LABEL vocabulary — Motifs: `<MOTIF>` (motif + substitution position), `BARCODE` (*"The variant of the same sequence"*), `RESTRICTION`, `ADAPTER`, `DUPLICATE_RESTRICTION_SITES`. SNPs: `<SNP>`, `<NUCLEOTIDE>` (REF/ALT), `BARCODE`, `RESTRICTION`, `ADAPTER`, `DUPLICATE_RESTRICTION_SITES`, plus indel-normalization notes (`"| added N nucleotides from the left"`, `"| removed N nucleotides from the right edge"`, from `parseSNPs.py::trimSequenceAndMakeHeader`). Transmutation: `Mutated_nucleotides - 3 | Scrambled - No | Reversed - No | Complemented - Yes`. So it goes **beyond naming the mutation** — it records which sub-components were added, in what role, and which construction-time warnings fired. But it is a delimited string inside the FASTA header, not a separate structured record/table, and there is no machine-readable side-car (no TSV/JSON manifest). Hence "partial".

**`automatic_naming` — yes.**
Headers are generated entirely by the tool, not supplied by the user. `oligo.py`: `headerString = "|".join(["-".join([i[0], str(i[1])]) for i in header])`; documented example *"> ATGTG - 53|AAAAA-61|RESTRICTION - 1|RESTRICTION - 2"* with the explanation *"ATGTG - 53 is the motif starting at position 54 in the background sequence."* Confirmed live for SNPs and Transmutation (§2).

**`design_visualization` — partial.**
Highlighted sequences only, no design/graph view. Docs, Motifs results: *"Image below is an excerpt of the Results page ( showing 8 nucleotides ). The nucleotides in red are the substituted motifs."* Supplementary §2.1: *"Substituted motifs are colour-marked for visualization purposes."*; §2.3: *"Mutated nucleotides are colour-marked for visualization purposes."* Implementation: `oligo.py` builds a parallel `sequenceHTML` with `<span style='color:red'>` around motif spans; `myfunctions.py::highlightString`. There is nothing resembling a graph/DAG view or a library-level summary plot.

### Block B — assay coverage

**`assay_dms` — no.** Searched paper, Supplementary, docs page and the whole repo tree: no codon, no amino acid, no ORF, no reading-frame, no translation-table concept anywhere. Mutagenesis is strictly nucleotide-level (motif substitution, SNP substitution, random point mutation, scramble/reverse/complement). Title and every module name are MPRA-scoped.

**`assay_mpra` — yes.** Title: *"MPRAnator: a web-based tool for the design of massively parallel reporter assay experiments."* Abstract: *"We introduce MPRAnator, a set of tools that facilitate rapid design of MPRA experiments."* The Supplementary also notes portability to related protocols: *"the output sequences are not restricted to MPRA experiments, but they can easily be linked to similar protocols such as BunDLE-seq (Levo et al., 2015)."*

**`assay_insilico` — no.** No model, predictor, scoring function, or ML integration of any kind in paper, Supplementary, docs, or repo. PWM Seq-Gen is a sequence *generator* from a PWM, not an in-silico probe of a genomic AI model.

### Block C — genomics integration

**`genome_coordinates` — partial.**
The SNP module *requires* coordinates, but only as text in the FASTA header. `mycustom.py::FastaFileWithHeaderRange.getHeaderRange` parses `chr(\w+):(\d+)-(\d+)` and raises *"Header does not contain chromosome and positions in the proper format. eg. 'chr1:start-end'"*; the live sample header is `>sequence1chr1:42990-43041,testest`. `parseSNPs.py::snpInChromosomeRegion` then checks each VCF record against that interval. There is also a `BedFile` parser in `mycustom.py` and a homepage utility: *"Click here for a Python script that converts original coordinates to equal sized, centered coordinates with range X nucleotides on each side."* **However, MPRAnator never reads a reference genome** — the user must paste/upload the actual sequence; there is no genome build selector, no FASTA index, no coordinate→sequence retrieval. Chromosome validation is hard-coded human-only (`[str(i) for i in range(1,23)] + ["X","Y"]`, with the error *"Chromosome numbers are between 1 and 22"* / *"Only 'X' and 'Y' are accepted for chromosome letters"*).

**`transcript_models` — no.** No GTF/GFF parsing, no gene/transcript/exon concept anywhere in the repo, docs, paper, or Supplementary. Input formats are FASTA, VCF-like 5-column text, BED intervals, and PWM matrices.

**`exon_intron_split_codons` — no.** Follows from the above: no exon/intron model and no codon model at all.

**`hgvs_input` — no.** The only variant input is a whitespace-delimited VCF-like table. `mycustom.py::SnpFile` requires ≥5 columns and validates columns 4/5 as `ATGCU.` only, raising *"less than 5 columns, probably not a VCF file"*. No `c.`/`g.`/`p.` parsing.

**`vcf_vep_output` — no.** VCF is an *input* format only (form label: *"Enter your SNPs (VCF Format)"*). Every documented output is FASTA plain text: *"The result page (plain text view) displays the synthesized oligonucleotides in FASTA format."* No VCF, GFF, or annotation-tool-consumable emission; nothing about VEP anywhere.

**`consequence_annotation` — no.** Headers record the nucleotide change (`Nucleotide change T:G`) and, for indels, the length-normalization performed — but never a molecular consequence. There is no gene model, so a consequence could not be computed. Confirmed by the exhaustive documented LABEL vocabulary in §4/`per_sequence_provenance`, which contains no consequence term.

### Block D — physical construction

**`primer_design` — no.** No primer, Tm/melting-temperature, annealing, PCR, or oligo-assembly-primer functionality in paper, Supplementary, docs, or repo. The only wet-lab-facing elements are user-supplied fixed adapter and restriction-site strings that are concatenated verbatim.

**`codon_optimization` — no.** No codon usage table, no organism selector, no synonymous-recoding logic. (Same absence of any codon concept as `assay_dms`.)

**`synthesis_constraints` — partial.**
Two real, documented checks/enforcements, both narrow:
1. *Restriction-site collision detection in the final oligo.* Docs LABEL: *"DUPLICATE_RESTRICTION_SITES : The restriction site which has multiple copies present."* Supplementary §2.1: *"If any of the restriction sites are identified in any generated oligonucleotide sequence they are reported in the header of the corresponding sequence."* Implementation in `part1.py::createMPRAResultOutput`: `if len(myf.findMatch(sequenceS=outputSequence.lower(), motifS=restriction1.lower())) > 1: outputHeader += "|DUPLICATE_RESTRICTION_SITES - RESTRICTION1"`. Note this **reports**, it does not repair or filter.
2. *Uniform oligo length, explicitly motivated by synthesis.* Supplementary §2.2: *"Since several methods of oligonucleotide synthesis work best when all oligos have similar lengths, the instances with an insertion must be trimmed while the sequences with a deletion must be expanded. MPRAnator solves this problem by adding adenines to one end of the sequences that are too short, and by trimming one end of the sequences that are too long."* (`parseSNPs.py::trimSequenceAndMakeHeader`.) The Motifs module separately enforces `lengthsSame()` on inputs (*"Sizes of the sequences should be the same"*).

What is **absent**: no GC-content check on the full oligo (GC bounds apply to barcodes only), no homopolymer/repeat limit, no secondary-structure or synthesis-difficulty scoring, no vendor-specific screening, no length cap tied to a synthesis platform.

### Block E — engineering

**`interface`** — Web GUI (Django app, jQuery-UI drag-and-drop for sub-component ordering) plus an HTTP POST "REST API". Paper: *"The REST API allows programmatic access to MPRAnator using simple URLs."* In practice the API is two downloadable Python 2 scripts that `requests.post()` the same form parameters to `https://genomegeek.com/MPRAResults/` and `https://genomegeek.com/MPRAResults/SNPs/` and write the FASTA response to a file (`downloadables/MpraMotifs_script.py`, `MpraSnps_script.py`; the live site serves them from `/MPRA/ScriptDownload/MpraMotifs/` and `/MPRA/ScriptDownload/MpraSnps/`, and the live endpoints are `/MPRA/Motifs/Results/` and `/MPRA/SNPs/Results/`). There is **no library API, no CLI, and no installable package**; the repo has no entry point other than Django's `manage.py`. Note the repo's `MpraSnpApiView` is a stub returning the literal string `"testing"`, and `restView` just redirects home. Only 2 of 4 modules have documented programmatic access.

**`license`** — Paper claims MIT; the repository contains no license file and GitHub reports `license: null`. See §3.

**`maintained`** — Last (and only) commit `2015-12-27`; GitHub `updated_at` 2016-01-20. No releases, no tags, no issues activity. ~10.6 years stale as of 2026-08-10. The deployed web app is newer than the repo but its source is not public.

---

## 5. The tool's own documented example libraries / use cases

Reproducing these in PoolParty is feasible for #1 and #3–5; #2 (the full worked motif use case) is the only one with a stated expected library size.

1. **MPRA Motifs "Load Sample Data"** (live at `/MPRA/Motifs/`, and mirrored in `downloadables/MpraMotifs_script.py`):
   scaffold `>sequence1 AGAGATTAGTCAGTACGGCTAGCTAGCTACGTCTATATTATAGCGATACGGG` (52 nt); motifs `>motif1 GGCTT`, `>motif2 AAACGG`; `minSpacing=2`, `maxSpacing=4`, `leftDistance=3`, `rightDistance=3`, `frequencyOfInsertion=5`; barcodes length 12, GC 20–80 %, edit distance 2, 4 barcodes/sequence. **This is the input that currently triggers the HTTP 500.**

2. **"MPRA Motif Use Case Example"** — the fully worked design on the live documentation page, and the strongest reproduction target because it states an expected count:
   > *"Inspired by the results presented by Nguyen et al, we would like to investigate the effects of AP1 (TGACTCA), ELK1 (ACCGGAAGT) and RFX (CGTTGCTAGGCAACG) on gene expression."*
   Backgrounds: two non-regulatory mouse (mm9) tiles, `>bg1_chr6:77195320` and `>bg2_chr9:37271330`, 90 nt each (sequences given verbatim on the page). Total tile length 117 bp; two flanking restriction sequences `CACGTG` and `CAATTG` *"rearranged so that they flank the background"*; barcodes 15 bp, GC 40–60 %, edit distance 3, 6 barcodes per sequence; minimum 15 bp to each edge; min/max spacing 6 and 24; interval of substitution 12.
   > *"This results in a total of 5856 sequences."*

3. **MPRA SNPs "Load Sample Data"** (live at `/MPRA/SNPs/`; mirrored in `downloadables/MpraSnps_script.py`):
   `>sequence1chr1:42990-43041,testest` + 52 nt scaffold; SNPs `1 42991 rs43434 T G` and `1 43030 rs121212 A G`; barcode length 12, GC 20–80 %, edit distance 2, 4 barcodes/sequence. **Verified working today** — output shown in §2.

4. **Transmutation "Load Sample Data"**: `>sequence1 AGAGAAGGACATATTCGATCGGATCGAT`, `>sequence2 AGAGACCCATAGAACACA`, `numToMutate = 3`. **Verified working today.**

5. **PWM Seq-Gen "Load Sample Data"**: two JASPAR matrices pasted as text — `>MA0056.1 MZF1_1-4` and `>MA0089.1 NFE2L1::MafG` (4 rows × 6 columns each, counts) — with `simNum = 3`. Documented output header example: `> MA0056.1 MZF1_1-4 | Simulation number - 4`.

6. **Sub-component ordering example** (Supplementary Fig. 3): *"Modular design for final output. Sub-components can be placed in any order. Replicates of each sequence can be generated, each with distinct barcode sequence."* Ordering vocabulary in the API scripts: `Barcode, Background, restriction site 1, restriction site 2, adapter site 1, adapter site 2`.

---

## 6. Notable capabilities NOT covered by the current row list

These are candidates for new matrix rows, or at least for a footnote so a referee cannot say we ignored what the tool does well.

- **PWM → sequence generation.** Two modes: probabilistic realizations of a PWM, or exhaustive enumeration of all k-mers exceeding a user probability threshold (Supplementary §2.4). Nothing in Blocks A–E captures "motif sampling/enumeration from a PWM".
- **Negative-control / scramble generation as a first-class feature.** Scramble, reverse, complement, and N random point mutations (Transmutation), with the choices recorded in the header. A row like `negative_control_generation` would be a fair addition — this is genuinely a capability PoolParty should be compared on.
- **Modular oligo assembly with user-orderable sub-components** (adapters ×2, restriction sites ×2, barcode, background) via drag-and-drop. This is "flanking/scaffold assembly", distinct from `combinatorial_motif_place`.
- **Restriction-site collision reporting in the assembled oligo** (`DUPLICATE_RESTRICTION_SITES`) — a cloning-compatibility check rather than a synthesis check; partially credited under `synthesis_constraints`.
- **Indel handling with length normalization** — insertions/deletions up to 10 nt from the VCF input, with automatic A-padding or edge-trimming to keep all oligos equal length, and the operation recorded in the header.
- **Barcode replicates / sequence multiplicity** — N distinct-barcode copies of every designed sequence in one pass (up to 5 on the web form).
- **Reference/WT sequences emitted alongside variants** in the SNP module (`| REFERENCE` header), without the user asking for it.
- **BED coordinate-centering helper script** — converts arbitrary intervals into equal-sized windows centered on the original, X nt each side (offered on the homepage as a download).
- **SNP module doubles as a plain assembly tool**: Supplementary §2.2 — *"the MPRAnator SNP design tool can accept a set of FASTA sequences without a list of associated SNPs and perform their incorporation into the user-designed constructs. This is particularly useful if the regulatory role of existing sequences is being investigated."*
- **Hard library-size caps** (a *limitation* worth a row if we ever add one on scalability): `motifs^4 × sequences ≤ 100`, ≤1000 input sequences for Motifs / ≤50000 for SNPs, ≤20000 nt per sequence, ≤5 barcodes/sequence, indels ≤10 nt, ≤4 motifs, ≤2 restriction sites, ≤2 adapters.

---

## 7. Corrections to the prior analysis

The prior analysis (`prior_analyses/Georgakopoulos-Soares2017_MPRAnator_analysis.md`) is directionally right but has three claims that would not survive referee scrutiny:

1. **"Sequence metadata | None"** — wrong. MPRAnator emits structured pipe-delimited DESCRIPTOR headers recording motif identity + position, barcode replicate index, restriction/adapter presence, restriction-site collisions, SNP id/position/change, indel length-normalization, and scramble/reverse/complement flags. The honest distinction is *string-encoded header provenance* vs. PoolParty's *structured design cards / programmatic metadata*, not "none vs. some".
2. **"Design cards | None" and the implication that MPRAnator has no visual inspection** — the results page colour-marks substituted motifs and mutated nucleotides in red. Again a difference of kind (highlighted output sequences vs. a design/graph view), not a total absence.
3. **Implementation attribution** — the prior note says "Implemented in Python, Perl and Javascript". That is the paper's claim; the public repo contains no Perl at all (Python 2 + Django templates + jQuery only). Do not repeat the Perl claim as if verified.

Two further points the prior analysis missed entirely and that strengthen the manuscript's positioning if stated carefully: the repository has been untouched since a single commit in **December 2015** and ships no license file despite the paper's MIT claim, and the **MPRA Motifs module of the live service returns HTTP 500 today**.
