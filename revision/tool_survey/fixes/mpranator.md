# MPRAnator — repair change log

**Record repaired:** `final/mpranator.md` (edited in place)
**Audits processed:** `citation_audit/mpranator.md` (20 findings + link check), `factcheck/mpranator.md` (A1–A9, B1–B8, C)
**Repair date:** 2026-08-14
**Verification base:** local PDF via PyMuPDF; `btw584_supp.docx` re-downloaded from Europe PMC (`PMC5198521/supplementaryFiles`) and parsed from `word/document.xml`; the repository cloned at commit `9969790d62410138d4281b7955da6d085f07b1bc` and read/grepped locally; GitHub REST (`git/trees/master?recursive=1`, `/forks`, `/search/repositories`, repo metadata); bio.tools API; PyPI; Bioconda; ~25 fresh live HTTP requests to `genomegeek.com` on 2026-08-14, including running both distributed API scripts verbatim under Python 3 + `requests` 2.32.2.

**No capability value was changed.** Two findings that would have required changing a value are escalated below.

---

## Applied (verified, corrected)

### CIT-1 / "44-entry recursive tree" → 50
**Verified:** `git ls-tree -r -t --name-only HEAD | wc -l` = **50** (41 blobs + 9 trees); `git/trees/master?recursive=1` returns 50 entries, `truncated: false`, same 41/9 split.
**Correction:** `44` → `50` at all four occurrences (§0 row 2, `assay_dms` source column — now "50 entries: 41 files + 9 directories" —, `interface`, `license`). The absence of `setup.py`/`requirements.txt`/`Dockerfile`/`LICENSE`/`COPYING` re-confirmed and left unchanged.

### CIT-2 / `additional_capabilities` "+5" → "+4"
**Verified:** `extractions/mpranator.md` §6 has 10 bullets; the final has 14 numbered items; exactly 4 carry *(ADDED on review)*. The extraction's first bullet already reads "**PWM → sequence generation.** Two modes: probabilistic realizations of a PWM, or exhaustive enumeration…", so dual-mode operation was not an addition.
**Correction:** §0 row 10: `+5` → `+4`, and "PWM Seq-Gen dual-mode operation," struck from the list of additions.

### CIT-3 / FC-A1 / "the four live forms submit a JS-generated `ordering` field"
**Verified:** fetched all four live form pages. `ordering` occurs twice in `/MPRA/Motifs/` and `/MPRA/SNPs/` (the `#sortable` drag list plus `$("#mpraForm").append("<input type='hidden' … name='ordering'>")`) and **zero times** in `/MPRA/Transmutation/` and `/MPRA/PWM/`.
**Correction:** §2 "the four live forms" → "the live Motifs and SNPs forms … (Transmutation and PWM Seq-Gen have no such field and no such JavaScript)". The missing-`ordering` → HTTP 500 claim was independently reproduced (live SNP POST omitting `ordering` → 500; mechanism visible in `views.py`, `request.POST.get('ordering').strip()` on `None`) and left as written.

### CIT-4 / "no fork or successor"
**Verified:** `/repos/hemberg-lab/MPRAnator/forks` returns **one** fork, `xtmgah/MPRAnator`, `pushed_at` 2015-12-27T15:29:04Z — identical to the canonical repo, i.e. no newer commit. `search/repositories?q=MPRAnator` returns `total_count: 1` (default repo search excludes forks).
**Correction:** §1 source row and the `maintained` cell now read "no successor" and name the single fork with its identical `pushed_at`. `maintained = no` is untouched and unaffected.

### CIT-5 / FC-A3 / "`Part3Form` lacks the reverse/complement options"
**Verified:** `iliasApp/forms.py` `Part3Form` declares `reverse = forms.BooleanField(label="Reverse motif?" …)` and `complement = forms.BooleanField(label="Complement motif?" …)`; `views.py::part3RresultsView` reads `request.POST.get('reverse')`/`('complement')` and applies `seq[::-1]` / `myf.complement(seq)`. The claim is contradicted by its own cited source.
**Correction:** §6 — clause deleted; the bullet now reads "a **development snapshot that predates the published tool** — no PWM Seq-Gen code. The deployed source is not public." (The PWM absence is real: no PWM route in `urls.py`, no PWM code in the tree.)

### CIT-6 / FC-A2 / "Headers are generated entirely by the tool; the user never supplies them"
**Verified:** `oligo.py` line 74 emits `"> Background-%s|%s" % (backgroundSequenceHeaderS, headerString)` — `backgroundSequenceHeaderS` is the user's FASTA header. `views.py` line 210 builds Transmutation headers as `">"+seq.name+"| Mutated_nucleotides …"`. The record's own live examples show it (`>sequence1| Mutated…`, `>MA0056.1 MZF1_1-4 | Simulation number - 1`).
**Correction:** `automatic_naming` evidence now says the *descriptor* part is tool-generated and the user never writes an output name, but the user's input identifier is carried through as the header prefix. **Value `yes` untouched** — see escalation ESC-3.

### CIT-8 / FC-A5 / "both shipped scripts 500 as distributed"
**Verified by running the distributed scripts verbatim** (downloaded 2026-08-14 from `POST /MPRA/ScriptDownload/{MpraMotifs,MpraSnps}/`, both HTTP 200):
- `requests.post('http://www.genomegeek.com/MPRA/Motifs/Results/', …)` → history `[301 → https://…/Results/, 302 → /MPRA/Motifs/]`, final **200 text/html** (the input form). Identical for SNPs.
- the next line, `openfile.write(result.content)`, then raises `TypeError: write() argument must be str, not bytes` under Python 3 — reproduced for both scripts.
- the same payloads POSTed directly to the **HTTPS** result endpoints → **500** for both (Motifs with empty barcode fields; SNPs with its shipped `numOfBarCodesPerSequence = "4"`).
**Correction:** the mechanism was rewritten in the `interface` cell, the §7 "Shipped REST API scripts" row and the §9 REST-API row; §0 row 9 and `barcode_generation` now say "the payloads of both shipped API scripts" / "the barcode payload shipped in the `MpraSnps` API script". The record's conclusion — no working demonstrator for the paper's REST-API claim — survives and is unchanged.

### CIT-9 / Motifs sample "mirrored in `MpraMotifs_script.py`"
**Verified:** live `#loadData_id` handler sets `leftEdge=3`, `rightEdge=3`, `frequencyOfInsertion=5`, barcode 12/20/80/2/4; the served `MpraMotifs_script.py` ships `leftDistance=2`, `rightDistance=2`, `frequencyOfInsertion=1` and empty barcode fields. Scaffold and both motifs do match. (The SNP preset *is* mirrored accurately — left unchanged.)
**Correction:** §8 item 1 now says "only partly mirrored … which ships `leftDistance=2`, `rightDistance=2`, `frequencyOfInsertion=1` and empty barcode fields".

### CIT-10 / `generatePermutations` / `equaliseSize` cited as the implementation path
**Verified:** `grep -rn` over the clone — `equaliseSize` has **no call site at all**; `generatePermutations` is called only at `oligo.py:95`, inside `if __name__ == "__main__"`. `views.py:79` calls `part1.generateCombinations(motifsL)` and passes each combination straight to `oligo.oligo()`.
**Correction:** `combinatorial_motif_place` evidence now credits the position enumeration and explicitly labels the two functions as off the Django path; the `mixed_mutagenesis_one_pool` source column drops `/generatePermutations`. **Value `yes` untouched** — the Supp. Fig. 1 "all possible permutations" quotation and the `itertools.product` position enumeration (which is what actually realises distinct spatial orders) both stand.

### CIT-11 / FC-A7 / "indels ≤10 nt"
**Verified:** Supp. §2.2 does say "Deletions or insertions (up to 10nt)". `parseSNPs.py` generates only when `abs(difference) < 10`; `MpraSnpsForm.clean` raises "SNP sizes cannot be greater than 10!" from `SnpFile.getMaximumTargetSize()`, which measures **ALT allele length**, not indel size. Live 2026-08-14: `T → TGGGGGGGGG` (9-nt difference) emitted "removed 9 nucleotides from left edge"; an 11-character ALT was rejected at validation.
**Correction:** §6 bullet now states the Supplementary claim and the implemented `< 10` boundary plus the separate ALT-length check, with the live result.
**Local roughness accepted:** §2 module list and §4 item 7 still say "indels ≤10 nt" / "up to 10 nt". Left alone deliberately (they paraphrase the Supplementary). A §9 "claim vs. deployment" row would be the natural home if the authors want the mismatch surfaced there too.

### CIT-13 / "reviewer marked 21 of 24 verdicts supported"
**Verified:** `reviews/mpranator.md` line 17: "**Twenty of the twenty-four cells I would leave alone.**" No "21", "21/24" or "twenty-one" occurs anywhere under `reviews/`, `extractions/`, `v2/` or `verified/`.
**Correction:** §10 now reads "The reviewer would leave 20 of the 24 cells alone".

### CIT-14 / "zero-hit grep across the full recursive tree"
**Verified — but the audit's stated mechanism is wrong.** Over the clone, `codon|amino|protein|translat|orf|reading.frame` matches **one file**: `iliasApp/static/iliasApp/bootstrap.min.css`, 14 times, all via **`translat`** matching CSS `translate(…)`. The `orf` alternative the audit blames matches **nothing** (`transform` does not contain `orf`). Every biological term (`codon`, `amino`, `protein`, `orf`, `reading.frame`) is zero-hit everywhere, and the whole pattern is zero-hit across all `.py`/template files.
**Correction:** `assay_dms` now says "across every source file in the tree (the pattern's only match anywhere is `translate` inside the vendored `bootstrap.min.css`)"; `codon_optimization` now says "the same source-file zero-hit grep". The correction is based on my own grep, not the audit's reasoning.
**Also checked and left alone:** the analogous grep claims for `gtf|gff|exon|intron|splice|transcript` (`transcript_models`) and `primer|melting|anneal|PCR|\bTm\b` (`primer_design`) are genuinely zero-hit over the full tree — every term returns 0 matches.

### CIT-15 / "`mycustom.py` contains exactly three parsers"
**Verified:** `mycustom.py` defines `GenericFile`, `FastaFile`, `Sequence`, `FastaFileWithHeaderRange`, `SequenceWithRange`, `BedFile`, `SnpFile`. `FastaFileWithHeaderRange` is a distinct parser (it overrides `getSequences` and adds `getHeaderRange`) and is the class the SNP module actually uses.
**Correction:** `transcript_models` evidence now reads "parses exactly three input-format families: FASTA (`FastaFile`, plus the `FastaFileWithHeaderRange` subclass), BED (`BedFile`) and VCF-like SNP tables (`SnpFile`)". Value `no` untouched.

### CIT-16 / FC-A6 / IUPAC claims
Three sub-claims checked; two corrected.
- **"motif and restriction-site matching are degenerate-code aware" — wrong.** `regex_mapper` is called only by `myf.findMatch`, and `findMatch` is reached only from `part1.py:280/285` (the duplicate-restriction-site check) and from `Invalid_data`/`oligosynthesizer`, neither of which has a call site. Motif placement in `oligo.py` is plain slice substitution with no matching. **Corrected** to "restriction-site collision matching is degenerate-code aware", with the reason.
- **"`myfunctions.complement` carries the full IUPAC complement map" — wrong.** `complementD` maps `"S": "W"` and `"W": "S"`; correct IUPAC complements are S→S and W→W. **Corrected** to "a complement map for all 15 DNA ambiguity codes, though it maps `S`↔`W` rather than each to itself".
- **`U` is not accepted.** `mapperDict` and `complementD` have no `U`; live Transmutation POST of `>s1 / AUGCAUGCAUGC` returned "Invalid letters!" on 2026-08-14, against the docs' "The input will accept all nucleic acid IUPAC letters". **Added** as one closing sentence to the same item.
- **Not changed:** the word "throughout" in the §4 item 3 heading and in §11 item 5. That is an emphasis judgement, not a checkable fact; the concrete corrections now sit immediately beneath it. Noted as tension — §11 item 5 still reads "IUPAC-aware throughout (regex expansion for matching …)".
- **Not applied:** the audit's point that `'R': '[A|G]'` also matches a literal `|`. True as read, but it has no effect on any claim the record makes, and editing it would add unaudited text for no gain.

### CIT-17 / FC-A9 / "LABEL vocabulary … is exhaustive"
**Verified:** the live documentation page has two LABEL lists. Motif section: `<MOTIF>`, `BARCODE`, `RESTRICTION`, `ADAPTER`, `DUPLICATE_RESTRICTION_SITES`. SNP section (same page): `<SNP>  :  Name of SNP`, `<NUCLEOTIDE>  : REF / ALT`, `BARCODE`, `RESTRICTION`, `ADAPTER`, `DUPLICATE_RESTRICTION_SITES`.
**Correction:** `consequence_annotation` evidence now names both lists instead of calling the Motif list exhaustive. Value `no` untouched — neither list contains a consequence term.

### CIT-19 / FC-A8 / "no importable module"
**Verified:** the tree contains ordinary Python modules (`mycustom.py`, `part1.py`, `part3.py`, `oligo.py`, `parseSNPs.py`, `myfunctions.py`) and classes `FastaFile`, `Sequence`, `SnpFile`, `BedFile`. They are undocumented Django internals with no packaging, and Python-2-only (`xrange`, `print` statements), but the absolute wording is false.
**Correction:** §2 → "no importable library API (the repo's modules are undocumented Django internals)"; §10 item 4 → "no importable library API". Values untouched.

---

## Applied — section B (missing facts added as a clause or one sentence)

| # | Fact added | Where | How verified |
|---|---|---|---|
| B1 | The 2015 result view truncates silently to the first 20 background records and first 800 nt of each | §6, repo-caps paragraph | `views.py`: `numOfSequencesToUse = 20`; `[i.sequence[:800] for i in sequenceO[:numOfSequencesToUse]]` |
| B2 | One random realisation per input; mononucleotide shuffle; fixed order scramble → reverse → complement → mutation; no seed, no replicate count | §4 item 4 | `part3.scramble_motifs` = `rd.shuffle`; `views.py:184–205` applies the four options in that order; no `seed(` anywhere in the tree |
| B3 | Live PWM bounds: threshold 0–1 (default 0.1), 1–1000 simulations (default 1) | §4 item 1 | live `/MPRA/PWM/` HTML: `threshold value="0.1" min="0" max="1"`, `simNum value="1" min="1" max="1000"` |
| B4 | Neither API client exposes its module's full form (Motifs omits `reverseComplement`, SNPs omits `makeSnpCombinations`) | `interface` cell | data dictionaries of both scripts downloaded 2026-08-14 — neither key present |
| B5 | The barcode filter compares the *maximum matching-position count* against the requested distance, never rechecks after regeneration, and skips barcoding unless length + both GC bounds + replicate count are all supplied | `barcode_generation` caveat (i) | `part1.barcode_generator` lines 52–88 (`countlist` accumulated across candidates, `if max(countlist) < diffs`, retry loop tests only membership/GC); `part1.getBarCodes` `if barCodeLength and minimumGCContent and maximumGCContent and numOfBarCodesPerSequence` |
| B6 | SNP parser wants ≥5 whitespace columns with a bare chromosome token and no `#` header; REF is never checked against the background; a background with no in-range variant is dropped entirely | §6 structural bullets | `mycustom.SnpFile.getSnps`; `parseSNPs.py:145–190`; live 2026-08-14 — two backgrounds submitted, variant matching only the first, output contained only the first |
| B7 | The form permits a substitution interval of 0, which `oligo.py` uses as a modulus divisor | §6, repo-caps paragraph | live `frequencyOfInsertion` has `min="0"`; `oligo.py:70` `sorted(positionsL)[0] % frequencyOfInsertion == 0` |
| B8 (part) | "no published version number or release date" | §7 Repository row | GitHub: 0 tags, 0 releases, no version declaration. The bio.tools half of this finding is escalated (ESC-4) |

---

## Rejected

### FC-A4 — "SNP module doubles as a plain assembly tool" is not a current capability
**The audit is wrong about the live tool.** It asserts "the current live form also rejects an empty SNP list".
**Evidence (live, 2026-08-14):** a POST to `/MPRA/SNPs/Results/` with `SnpS=''`, `restriction1=CACGTG`, `adapter1=TTTTAAAA`, `usingDownload=True` returned **HTTP 200, `text/plain`**:
```
>1 42990 43041|RESTRICTION - 1|ADAPTER - 1
ACAGATTAGTCAGTACGGCTAGCTAGCTACGTCTATATTATAGCGATACGGGCACGTGTTTTAAAA
```
That is exactly the Supp. §2.2 behaviour the record describes — FASTA sequences with no associated SNPs, incorporated into the user-designed construct. The audit's repo citation (`forms.py` "Enter SNPs!") is real but applies only to the 2015 snapshot, which the record already flags as predating the deployed tool. Record left untouched.

### CIT-12 — "5856 sequences" not reproducible
The record quotes "This results in a total of 5856 sequences" as the tool's own documented claim (verified verbatim on the live documentation page) and nowhere asserts that it recomputed or confirmed the number. The audit's own alternative (6648) is derived from repository code that the audit simultaneously acknowledges is not the deployed source, and the deployed Motifs module returns 500, so no independent count is obtainable. Nothing in the record is wrong; record left untouched.

### CIT-14 (mechanism) — see the applied entry
The audit's claim that the `orf` alternative matched `transform`/`text-transform` is false (`orf` is not a substring of `transform`; 0 matches tree-wide). The finding's *conclusion* was applied on the strength of my own grep, not the audit's reasoning.

---

## Escalated (record left untouched at these points)

### ESC-1 — CIT-7: "Only **2 of 4** modules have any programmatic access at all"
**What I verified:** all four live modules accept scripted POSTs. Motifs and SNPs are `@csrf_exempt` and take a bare `requests.post`. Transmutation and PWM Seq-Gen require the normal CSRF cookie + token (403 without) but then work: `/MPRA/TransmutationResults/` → 200 with results, `/MPRA/PWM/Results/` → 200 in both modes. Only two modules have a *shipped API script*.
**Why escalated:** "2 of 4 modules" sits inside the `interface` cell's **descriptive value for the rendered matrix**, so correcting it is changing a value, not evidence.
**Question:** should the `interface` descriptive value read "HTTP-POST form API (2 of 4 modules)" — i.e. counting offered/documented programmatic access — or be restated as "shipped API scripts for 2 of 4 modules; all four form endpoints are POST-scriptable, two of them only with a CSRF session"?
**Options:** (a) keep "2 of 4" and change the *evidence* sentence to "only 2 of 4 modules ship an API client"; (b) restate the descriptive value; (c) leave both as-is and add a footnote.

### ESC-2 — CIT-18: FASTA-as-format vs plain-text-as-transport
**What I verified:** a normal browser-style POST to SNPs/Transmutation/PWM returns an **HTML** results page (e.g. SNP 200, `text/html`, 3021 bytes); adding `usingDownload=True` returns **`text/plain`** FASTA (325 bytes) — both reproduced 2026-08-14. `views.py` uses `content_type='text/plain'` only on the download branch.
**Why escalated:** the wording spans four separate places in three sections (§2 opening sentence, §2 "structural fact", `library_as_object`, `vcf_vep_output`), so a fix is more than a single paragraph of new text; and I cannot verify which request mode produced the earlier pass's observations, so "every live output obtained this pass was plain-text FASTA" may simply be true of that pass. No capability value depends on it (the output is FASTA either way; `vcf_vep_output = no` is unaffected).
**Question:** do you want the FASTA/HTML-vs-plain-text distinction spelled out at those four points?
**Options:** (a) leave as-is; (b) add "(plain text in download mode; embedded in an HTML results page otherwise)" once in §2 and let the other three inherit it; (c) qualify all four.

### ESC-3 — value check arising from CIT-6 (`automatic_naming = yes`)
**What changed:** the evidence for `automatic_naming` no longer says headers are generated entirely by the tool; it now records that the user's input identifier is carried through as the header prefix.
**Why escalated:** I judged `yes` still correct (the user never writes an output name and all descriptors are machine-generated) and therefore did **not** touch the value — but the cell's evidence string is now materially different from the one the value was set against.
**Question:** does `automatic_naming` stay `yes` under the corrected evidence?
**Options:** (a) `yes` unchanged; (b) `partial` on the grounds that identifiers are user-seeded.

### ESC-4 — FC-B8 (bio.tools registry material)
**What I verified:** `https://bio.tools/api/tool/mpranator?format=json` → `version: []`, `download: []`, `license: "MIT"`, `toolType: ["Command-line tool", "Web application"]`, `operatingSystem: ["Linux","Windows","Mac"]`, `maturity: "Mature"`, `lastUpdate: 2024-11-24T21:03:06Z`. All of the audit's claims check out, and the "Command-line tool" + three-platform listing with no download is a genuine registry-vs-reality mismatch.
**Why escalated:** bio.tools is **not** among the record's §1 "Sources consulted". Adding these facts means either citing a source the record did not consult, or adding a source row asserting a consultation that this pass — not the original — performed.
**Question:** admit bio.tools as a source for this record?
**Options:** (a) no — drop the finding; (b) yes — add a §1 `registry` row and a §9 row for "bio.tools lists a command-line tool and three desktop platforms with no downloadable package"; (c) yes, but confine it to §7.

### ESC-5 — CIT-20: unverifiable provenance and request counts
**What I verified:** nothing survives to check. There is no request log, no response artefact and no host/session identifier anywhere under `revision/tool_survey/`. I did independently reproduce the *substance* of the observations (Motifs 500, barcode-enabled SNPs 500, combinations 500, missing-`ordering` 500, BED download 500, barcode-free SNPs 200, Transmutation 200, both PWM modes 200), but "third independent pass", "~30 fresh live HTTP requests", "~10 live POSTs", "three independent hosts/sessions" and "~20 payloads" cannot be recomputed.
**Why escalated:** these appear at 11 points including the record-status line; they are the record's own methodological attestations, which only the original author can confirm, weaken or withdraw. Removing them is a voice-level rewrite, not a factual correction.
**Question:** keep, soften, or drop the quantified provenance claims — and should this pass's reproductions be logged as an artefact so they are re-derivable?

---

## Findings deliberately not acted on

- **FC-C (balance).** The audit itself concludes "the imbalance is not material". Judgement reserved for the authors; no edit.
- **CIT link check.** Re-verified independently: PyPI 404, Bioconda 404, `sanger.ac.uk/science/tools/mpranator` → 302 → `/tool/mpranator/`, GenomeGeek pages 200, `ScriptDownload/{MpraMotifs,MpraSnps}` 200, `ScriptDownload/Bed` 500. Nothing to fix.

## Observation found during verification, not in either audit, not edited

`genome_coordinates` attributes two error strings to `mycustom.SnpFile`: "Chromosome numbers are between 1 and 22" and "Only 'X' and 'Y' are accepted for chromosome letters". Both strings exist verbatim in `mycustom.py`, but at lines 183 and 186, inside **`FastaFileWithHeaderRange.getHeaderRange`** — the very method the same sentence goes on to distinguish. `SnpFile`'s own message is "Chromosome position should be between 1-22 or X or Y" (line 270). The chromosome *list* attribution (`[str(i) for i in range(1,23)] + ["X","Y"]` in `SnpFile`, line 250) and the "`getHeaderRange` separately accepts up to 23" note are both correct. Reported rather than edited, since it is outside the scope of both audits; the cell value `partial` does not depend on which class raises which message.

---

## Supplementary quotations re-verified verbatim (no change needed)

Re-parsed `word/document.xml` from a freshly downloaded `btw584_supp.docx`. All of the record's Supplementary quotations match the source exactly: §2.1 motif-combination sentence, "select the barcode size (set to zero to exclude barcodes)…", "the user can select the multiplicity of each sequence…", the BunDLE-seq sentence, the restriction-site reporting sentence, both "colour-marked for visualization purposes" sentences; §2.2 "Although the default is to substitute all combinations of SNPs…", the oligo-length-normalisation passage, the no-SNP-list passage; §2.3 "therefore serving to generate negative controls…"; §2.4 "the header contains information regarding the probability of each k-mer occurring" and the duplicate-removal sentence; Fig. 1 and Fig. 3 captions. Live-output receipts also reproduced exactly, including `> MA0056.1 MZF1_1-4 | Prob : 0.10965375000000001`, `>MA0056.1 MZF1_1-4 | Simulation number - 1`, and the SNP `| REFERENCE` entry.

---

## Pass 2 — policy application

### ESC-1 — CIT-7 programmatic-access count — **applied**

**Edit:** In §0, `interface`, and §9, changed the claim from programmatic access for only 2 of 4 modules to the narrower verified claim that **2 of 4 modules ship API clients**. The `interface` evidence now also states that all four live form endpoints are POST-scriptable, with Transmutation and PWM requiring a CSRF session; its source column now cites the live HTML for all four forms. Value `partial` was not changed.

**Verification:** The official repository's complete `downloadables/` listing contains only `MpraMotifs_script.py` and `MpraSnps_script.py`. All four official live forms declare POST result endpoints. Scripted valid submissions reached all four endpoints: SNPs, Transmutation, and PWM returned results, while Motifs reached its downstream HTTP 500; bare Transmutation/PWM POSTs returned 403, and the same submissions with the official CSRF cookie/token flow returned 200 HTML. The corrected evidence still supports `interface = partial`: two shipped clients and a GUI exist, but no library API, CLI, or installable package does.

### ESC-2 — CIT-18 FASTA format versus HTTP transport — **applied**

**Edit:** Surgically qualified the five transport statements in §2, `library_as_object`, `lazy_evaluation`, and `vcf_vep_output`: FASTA is embedded in an HTML results page by default and returned as plain text on download. This corrects evidence only; no capability value changed. This was a factual correction, not a Policy-B omission.

**Verification:** Valid official sample submissions to SNPs, Transmutation, and PWM each returned HTTP 200 `text/html` with FASTA in the result page. SNP `usingDownload=True` returned `text/plain`; the result pages' official `/MPRA/TableDownload/` action returned the submitted FASTA as a `text/plain` attachment. Repository `views.py` likewise renders the results template outside its download branch.

### ESC-3 — `automatic_naming = yes` value-basis check — **applied**

**Edit:** No further prose edit was needed; the corrected evidence already names both halves of the behavior. Value `yes` remains unchanged.

**Verification:** Locked row 14 requires tool-generated informative identifiers composed from construction history and says that merely carrying through a user identifier does not count. Official `oligo.py` carries the user's scaffold header as a prefix but generates motif/position descriptors; `part1.py` adds barcode/restriction descriptors; `parseSNPs.py` generates SNP/position/nucleotide-change descriptors; and `views.py` generates Transmutation operation descriptors. The generated construction-history suffixes independently satisfy the locked test, so the user-seeded prefix does not undermine `yes`.

### ESC-4 — FC-B8 bio.tools registry material — **declined-by-policy**

**Edit:** None. The final record already contains no bio.tools citation or registry-derived claim.

**Verification/policy:** Policy C excludes third-party aggregators and permits only the repository, official documentation/package page, or paper. A search of `final/mpranator.md` confirmed that it has no reliance on bio.tools. The admissible repository evidence already supports the separate statement that there are no tags, releases, or published version number; the registry-only remainder was dropped.

### ESC-5 — CIT-20 provenance and request-count anecdotes — **applied**

**Edit:** Dropped the uncheckable independent-pass/host/session/request/payload counts and first-round narrative from the status line, source table, capability sources, §7, §10, and §11. Replaced process phrasing such as "this pass" with plain evidence labels, retained the verified substantive observations, and refreshed the §11 live-failure statement to the 2026-08-14 recheck.

**Verification:** A survey-wide filename/content search found no retained MPRAnator request log, response artefact, host/session identifier, or matching receipt from which the historical counts could be reconstructed. Current primary-source rechecks returned: Motifs sample with and without barcodes, HTTP 500; SNP combinations and barcode-enabled SNPs, HTTP 500; barcode-free SNPs, HTTP 200; Transmutation and PWM samples, HTTP 200; BED helper POST, HTTP 500. These checks support the plain status statements but cannot verify the deleted historical counts or independence claims.

**Pass-2 counts:** 4 applied · 1 declined-by-policy · 0 rejected-unverifiable · 0 escalated.

**Row-substitution candidates:** none observed.
