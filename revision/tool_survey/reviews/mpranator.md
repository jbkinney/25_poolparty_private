# Adversarial review — MPRAnator extraction

**Reviewed:** `extractions/mpranator.md` + the structured extraction JSON
**Review date:** 2026-08-10
**Posture:** assume the referee is Georgakopoulos-Soares, Jain or Hemberg. Try to falsify every cell.
**Method:** independent re-derivation from the local PDF (PyMuPDF), the Supplementary `.docx` re-downloaded
from Europe PMC and unzipped/parsed myself, the GitHub REST API, the raw source of 16 repo files fetched
from `raw.githubusercontent.com`, the live `genomegeek.com` HTML/JS, `sanger.ac.uk`, PyPI, Bioconda, and
**~25 fresh live HTTP POSTs** to all four modules. Nothing was installed or cloned.

---

## 0. Bottom line

The Block-A/B/C/D science content of this extraction is **very good**. I re-derived every quotation from
the paper, the Supplementary and the live docs page and they are all verbatim-accurate; I re-read the
source files and every code citation checks out. Twenty of the twenty-four cells I would leave alone.

The problems are concentrated in **Block E**, where all three cells (`interface`, `license`, `maintained`)
carry the value `"yes"`. That is not the convention used by the sibling extractions
(`valiant.md`, `dnachisel.md`, `ledidi.md`, `seqpro.md`, `tangermeme.md` all use descriptive strings there),
and for `maintained` the literal cell value **contradicts its own evidence string** and is flatly false.
A matrix row reading `maintained: yes` for a repo with one commit dated 2015-12-27 is the single most
attackable thing in this file — and it is attackable in the direction that *helps* MPRAnator, so it is a
credibility problem rather than a fairness problem, but it must be fixed.

Second, the **availability narrative is now provably worse than the extraction says**, and the extraction
should be corrected *upward in severity* (with care). The extractor wrote "MPRA SNPs — WORKS". It does not
work in its advertised configuration. See §2.

Third, there are **three evidence-level factual errors** a tool-author referee could use to chip at us. See §3.

---

## 1. What I independently confirmed (no change needed)

| Claim | Independent check |
|---|---|
| Repo metadata | GitHub API: `created_at` 2015-12-27T14:57:02Z, `pushed_at` **2015-12-27T15:29:04Z**, `updated_at` 2016-01-20, size 79 KB, 1 star, 1 fork, 0 issues, `archived:false`, `license: null`, description "A webtool to create DNA sequence variants", homepage `http://genomegeek.com` — all as claimed |
| Single commit | `/commits` returns exactly 1: `2015-12-27T14:57:53Z`, Naman Jain, "first commit" — confirmed |
| README 12 bytes | confirmed (`# MpraNator`) |
| No LICENSE anywhere | confirmed via recursive `git/trees/master?recursive=1` (44 entries, no LICENSE/COPYING) |
| Paper claims MIT | verbatim from PDF: *"The source code is available on www.github.com/hemberg-lab/MPRAnator/ under the MIT license."* — confirmed |
| Paper claims Perl | verbatim: *"implemented in Python, Perl and Javascript"*; I grepped the whole tree — **no Perl**. Extractor's caveat is right |
| No PyPI / conda | `pypi.org/pypi/mpranator/json` → 404; `bioconda.github.io/recipes/mpranator/` → 404 |
| GitHub search → 1 repo | confirmed (`total_count: 1`); no successor |
| All Supplementary quotes (§2.1–2.4, Figs 1–3) | I re-downloaded `btw584_supp.docx` from Europe PMC and parsed `word/document.xml` directly. **Every quoted sentence is verbatim-correct**, including the BunDLE-seq sentence, the adenine-padding sentence, the "single, pairs and triplets" sentence, and the Fig.1 "Up to 4 motifs" caption |
| All live-docs quotes | re-scraped `genomegeek.com/MPRA/documentation/` and rendered to text — all verbatim, including the DESCRIPTOR/LABEL vocabulary, the `> ATGTG - 53|AAAAA-61|…` example, the Levenshtein tooltip, and the 5856-sequence use case |
| `oligo.py` eager `allResults`, `itertools.product`, `testIfWithin`, `testIfFarFromRightEdge`, `sequenceHTML` red spans | confirmed line-by-line |
| `part1.py` `generateCombinations` / `generatePermutations` / `equaliseSize(...,4)` / `barcode_generator` Hamming-style `if short_term[pos]==bc[pos]` / `DUPLICATE_RESTRICTION_SITES` append | confirmed line-by-line |
| `mycustom.py` `chr(\w+):(\d+)-(\d+)`, `SnpFile` ≥5 cols + `ATGCU.` letter set, human-only chromosome list `range(1,23)+["X","Y"]`, `BedFile` | confirmed |
| `views.py` `MpraSnpApiView` returns literal `"testing"`, `restView` redirects home, all downloads `text/plain` | confirmed |
| No codon / amino / ORF / primer / Tm / PCR / GTF / GFF / exon / intron / splice / HGVS / homopolymer / secondary-structure token anywhere | I grepped the full tree for all of these — **zero hits**. The Block-B/C/D "no"s are genuinely absences, not unfound |
| Sanger page | `sanger.ac.uk/science/tools/mpranator` → 302 → `sanger.ac.uk/tool/mpranator/`, HTTP 200, landing page only, points to `http://genomegeek.com/` and to the 2015 bioRxiv preprint; its "example codes" links are dead. It describes only **three** modules (no PWM Seq-Gen) |
| Live service alive | `genomegeek.com` HTTP 200, still advertises 4 tools, contacts still listed |
| Sample data for all four modules | re-extracted from the `loadData_id` JS on each live page — matches the extraction exactly |
| Indel length normalization | **verified live end-to-end**, which the extractor only cited from the Supplementary. Insertion `T→GGG` returned `… \| removed 2 nucleotides from right edge`; deletion `TTT→G` returned `… \| added 2 Adenines bases in the right edge` |
| Red highlighting | **verified live** in raw HTML of a Transmutation run: `<span style='color:red;'>T</span>GAGAAGGACATATTC<span style='color:red;'>C</span>ATCGG…`. `design_visualization = partial` is correct and now has a live receipt |
| `lengthsSame()` | verified live: unequal input lengths → `"Sizes of the sequences should be the same"` |

---

## 2. The availability picture is materially worse than the extraction states

This is the most important correction and it must be worded carefully, because a referee may re-test.

### 2.1 MPRA Motifs — 500 confirmed independently

I reproduced the HTTP 500 from a **different host and session** than the extractor, across 8 more payloads:
the site's own Load-Sample-Data values verbatim; barcodes off; single motif; 1, 2, 3 and 4 motifs against a
52 nt scaffold; and the exact payload shipped in the live `MpraMotifs` API script. **All 500.** An empty-motifs
POST returns 200 with inline validation errors, and an unequal-length POST returns 200 with
`"Sizes of the sequences should be the same"` — so request parsing and form validation are healthy and the
fault is downstream in generation. This confirms the extraction.

### 2.2 MPRA SNPs — the extraction's "WORKS" is **overstated**

The SNP module works only in a stripped-down configuration. Verified today:

| SNP configuration | Result |
|---|---|
| Sample data, **no barcodes**, `numOfBarCodesPerSequence=0` | **HTTP 200**, real FASTA (singles + `\| REFERENCE`) |
| + restriction site, no barcodes | **HTTP 200**, `\|RESTRICTION - 1` header, site prepended |
| + insertion / deletion, no barcodes | **HTTP 200**, correct length-normalization headers |
| Site's own **Load Sample Data** (barcodes: len 12, GC 20–80, dist 2, 4/seq) | **HTTP 500** |
| Any barcode config at all (`numOfBarCodesPerSequence≥1` + `barCodeLength`) | **HTTP 500** |
| `makeSnpCombinations=on` (barcodes off) | **HTTP 500** |
| The **exact payload in the live `MpraSnps` API script** (ships with barcodes on) | **HTTP 500** |

So the live failure is not confined to one module. As of 2026-08-10:

- **Barcode generation is non-functional on the live service in both modules.** That is the tool's signature
  MPRA feature and the one PoolParty most directly overlaps.
- **SNP-combination design is non-functional.**
- **Both shipped "REST API" example scripts 500 out of the box**, because both ship a configuration the
  server can no longer execute. The paper's `"The REST API allows programmatic access to MPRAnator using
  simple URLs"` therefore has no working demonstrator today.
- The **BED coordinate-centering helper script download** (`POST /MPRA/ScriptDownload/Bed/`) also returns
  **HTTP 500**, while `MpraMotifs` and `MpraSnps` download fine. The extraction lists the BED helper as an
  available capability without noting it cannot currently be obtained.

### 2.3 PWM Seq-Gen — the extraction **understates** it

The extractor wrote *"live and validating … I did not push a fully valid submission through."* I did. Both
documented modes work:

- Simulation mode (`simNum=3`): `>MA0056.1 MZF1_1-4 | Simulation number - 1` / `TGGGGA`, etc. HTTP 200.
- Exhaustive-threshold mode (`simOrAllMotifs=on`, `threshold=0.1`): returns every k-mer above threshold
  **with its probability in the header** — `> MA0056.1 MZF1_1-4 | Prob : 0.10965375000000001` / `AGGGGA`, etc.

So the correct net statement is: **2 of 4 modules fully functional (Transmutation, PWM Seq-Gen);
1 partially functional (SNPs, only without barcodes and without combinations); 1 non-functional (Motifs).**
Not "3 of 4".

### 2.4 Suggested replacement wording for `availability_status`

> Usable today only as a hosted web service, and only in part. Repo `github.com/hemberg-lab/MPRAnator` is
> reachable but frozen at a single commit (2015-12-27), has no LICENSE despite the paper's MIT claim, is
> Python 2 / Django 1.x with `pdb.set_trace()` left in, and predates the published tool (no PWM Seq-Gen).
> No PyPI, conda, Docker, `setup.py` or `requirements.txt`. On 2026-08-10 the live service at genomegeek.com
> returned working output for Transmutation and PWM Seq-Gen; the SNP module worked only with barcodes and
> SNP-combinations disabled; every barcode-enabled request and every MPRA Motifs request returned HTTP 500,
> including the site's own "Load Sample Data" values and the two example API scripts the site itself ships.

---

## 3. Evidence-level errors a referee could use

### 3.1 `barcode_generation` — the "1–5" range is wrong for the live service

The evidence string says, attributing to the live pages:

> "Form fields on both /MPRA/Motifs/ and /MPRA/SNPs/: barCodeLength (min 10), … numOfBarCodesPerSequence (1-5)."

The live HTML says otherwise:

- `/MPRA/Motifs/`: `<input type="number" name="numOfBarCodesPerSequence" min="0" max="100">`, and
  `<input type="number" name="barCodeLength" min="10" max="25">`
- `/MPRA/SNPs/`: `<input type="number" name="numOfBarCodesPerSequence" min="0" max="100">`, `barCodeLength min="10"` (no max)

`1–5` is the **2015 repo** value (`forms.py: max_value=5`), not the deployed value. This error propagates
into `stated_limitations` ("≤5 barcodes per sequence") and into the two `documented_examples`. It matters
because the docs' own worked use case asks for **6 barcodes per sequence** — i.e. the extraction's stated
cap is contradicted by MPRAnator's own flagship example, which an author-referee would spot instantly.
Fix: say "up to 100 per sequence on the live forms (the 2015 repo capped it at 5)".

Related: the extraction's whole `stated_limitations` block is repo-derived and should be labelled as such.
I probed the live Motifs validator with 1–5 motifs; the `motifs^4 × sequences ≤ 100` rule does not fire
where the repo code says it should (4 motifs × 1 sequence = 256 still reaches the generator and 500s).

### 3.2 `assay_insilico` — one sentence in the evidence is factually false

> "PWM Seq-Gen generates sequences FROM a position weight matrix; **it does not score** or probe any genomic model."

It does score. In exhaustive mode it computes and reports a PWM probability per generated k-mer:
`> MA0056.1 MZF1_1-4 | Prob : 0.10965375000000001` (verified live), and Supplementary §2.4 says the header
*"contains information regarding the probability of each k-mer occurring."* The **cell value `no` is still
correct** — a PWM is not a sequence-to-function model and there is no model-in-the-loop design step — but
delete the "it does not score" clause. As written it is a one-line falsification for a referee.

### 3.3 `interface` — the API scripts are Python 3 on the live site, not Python 2

The evidence says *"two downloadable Python 2 scripts"*. True of the repo copies (`print "…"`, old URL
`https://genomegeek.com/MPRAResults/`). The scripts **actually served today** are Python 3 (`print(...)`)
and point at the current endpoints `http://www.genomegeek.com/MPRA/Motifs/Results/` and
`/MPRA/SNPs/Results/`. Say "Python 2 in the repo; the live site serves Python 3 versions", or just drop the
version. (Both still 500 as shipped — see §2.2.)

### 3.4 Two smaller ones

- `mixed_mutagenesis_one_pool` cites the form's opt-in `"Make Snp Combinations?"`. Supplementary §2.2 states
  the **opposite** default: *"Although the default is to substitute all combinations of SNPs in the input
  sequences, there is also the option to substitute only a single SNP at a time."* Paper-vs-deployment
  disagreement; worth one clause so a referee cannot claim we misread their Supplementary.
- `genome_coordinates` says chromosome validation is hard-coded `[1..22]+X,Y`. Exactly right for
  `SnpFile`; `FastaFileWithHeaderRange.getHeaderRange` actually accepts up to 23 (`int(chrNumber) > 23`).
  Trivial, but it is a direct code quotation, so make it precise.

---

## 4. Block E: all three cells are mis-valued

The sibling extractions record these as descriptive strings (`valiant.md`: "interface = CLI only",
"license = AGPL-3.0-or-later", "maintained = last commit 2024-04-22"; `dnachisel.md`, `ledidi.md`,
`seqpro.md`, `tangermeme.md` likewise). MPRAnator's three cells all say `"yes"`.

- **`maintained: "yes"` is wrong**, not merely mis-formatted. The evidence string in the same object says
  *"~10.6 years stale"*, *"Last (and only) commit 2015-12-27"*, *"No releases and no tags."* A rendered
  matrix cell reading "yes" is false and self-contradicting. Suggested: `no — single commit 2015-12-27;
  unmaintained for ~10.6 y; deployed service newer than the repo but closed`.
- **`license: "yes"` is wrong.** There is no license. Suggested: `MIT per the paper; no LICENSE file in the
  repository (GitHub API: license = null)`.
- **`interface: "yes"` is wrong and is the one that would most flatter MPRAnator.** In a comparison whose
  whole point is "web forms vs. a composable Python API", a cell reading "yes" reads as parity. Suggested:
  `Web GUI (Django + jQuery-UI drag-and-drop) + HTTP-POST "REST API" for 2 of 4 modules; no library API,
  no CLI, no installable package`.

Note the codongenie review flagged the same `maintained` convention problem. This should be settled once,
cross-tool.

---

## 5. Capabilities the extraction missed

The extractor's `additional_capabilities` list is unusually thorough (PWM generation, negative controls,
modular ordering, restriction collision, indel normalization, barcode replicates, WT emission, BED helper,
SNP-module-as-assembler, size caps). Three genuine gaps remain:

1. **IUPAC / degenerate-base awareness throughout.** `myfunctions.regex_mapper` expands R/Y/S/W/K/M/B/D/H/V/N
   into regex character classes, so restriction-site and motif matching are degenerate-code aware;
   `myfunctions.mapperDict` defines legal substitutions for every IUPAC code, so `part3.mutateString`
   mutates degenerate positions to a chemically consistent alternative and **explicitly skips `N` positions**
   (`while sequenceS[posChosen] == "N"`). The live docs say so: *"The input will accept all nucleic acid
   IUPAC letters."* Nothing in the extraction records that MPRAnator handles degenerate bases at all.
2. **Duplicate removal in PWM Seq-Gen.** Live form field `dupMotifs`; docs: *"Remove duplicates: Allows users
   to either keep only unique Sequence Motifs or include Duplicates as well"*; Supplementary §2.4 same.
   A library-level dedup control — relevant if the matrix ever gains a dedup/uniqueness row.
3. **Per-k-mer PWM probability reported in the output header** (`| Prob : 0.1096…`). This is a *score*
   attached to a designed sequence, i.e. a (weak) form of per-sequence quantitative annotation. Worth one
   line in `per_sequence_provenance`, and it is the thing that makes the `assay_insilico` evidence sentence
   in §3.2 falsifiable.

Two lower-value observations, useful only as footnotes:

4. `mycustom.FastaFile.areDuplicatesPresent()` and `lengthsSame()` are input-QC checks; the latter is
   user-visible (`"Sizes of the sequences should be the same"`, verified live).
5. **Implementation bug worth knowing about, not worth putting in the manuscript:**
   `myfunctions.revcompl()` is `complement(backgroundS[::1])` — `[::1]` is the identity slice, so it returns
   the **complement, not the reverse complement**. It is used by `findMatch(getRevCompMatch=True)`
   (restriction-site collision detection) and by the Motifs `reverseComplement` option. Do not put this in a
   referee response — attacking a competitor's correctness invites a fight we do not need — but it means the
   `DUPLICATE_RESTRICTION_SITES` check does not do what its name implies.

---

## 6. Judgement calls I would leave alone

The extractor flagged three cells as least-sure. I would keep all three as `partial`:

- **`library_as_object = partial`.** A strict "is there a library object?" reading gives `no`. But one POST
  does yield a whole library with no user looping, and the extractor's choice avoids an inflated negative
  against a possible author-referee. Keep `partial`, keep the "no object, only a text blob" clause — that
  clause is what carries the PoolParty contrast.
- **`dag_chaining = partial`.** The hand-off is authors-documented (PWM → Motifs; Transmutation → SNP/Motif
  controls) even though it is copy-paste. Keep `partial`; make sure other web-only tools in the survey are
  scored the same way.
- **`per_sequence_provenance = partial` / `design_visualization = partial`.** Both are correct and both
  correctly overturn the prior analysis's "Sequence metadata: None" / "Design cards: None". I verified the
  red-span highlighting live in raw HTML. Do **not** revert to the prior analysis's wording — those two
  lines are the easiest things in the whole survey for an author-referee to rebut.

---

## 7. Priority fix list

| # | Item | Change |
|---|---|---|
| 1 | `maintained` cell | `"yes"` → `"no — single commit 2015-12-27; ~10.6 y stale"`. Internal contradiction, must go |
| 2 | `interface` cell | `"yes"` → descriptive string; "yes" reads as parity with PoolParty's API |
| 3 | `license` cell | `"yes"` → `"MIT per paper; no LICENSE file in repo"` |
| 4 | `availability_status` | "SNPs WORKS" is false as advertised; barcodes and SNP-combinations 500 in **both** modules, as do both shipped API scripts and the BED helper download. PWM Seq-Gen is fully working (upgrade). Net = 2 working / 1 partial / 1 broken |
| 5 | `barcode_generation` evidence + `stated_limitations` | "1–5 barcodes/sequence" is the 2015 repo value; live forms allow 0–100, and MPRAnator's own worked example uses 6. Contradicts itself as written |
| 6 | `assay_insilico` evidence | delete "it does not score" — PWM mode emits `\| Prob : …` per k-mer. Keep the `no` |
| 7 | `interface` evidence | live API scripts are Python 3 with current URLs, not Python 2 |
| 8 | `additional_capabilities` | add IUPAC/degenerate-base handling, PWM duplicate removal, per-k-mer probability in header |
| 9 | `mixed_mutagenesis_one_pool` evidence | note Supp §2.2 states combinations are the *default*, contradicting the opt-in checkbox |
| 10 | HTTP-500 phrasing | keep the extractor's discipline ("observed to return HTTP 500 on 2026-08-10"). I reproduced it from a second host/session across ~12 payloads, so it is a server-side defect, but the phrasing should stay date-stamped |
