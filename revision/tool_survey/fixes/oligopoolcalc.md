# Repair change log — `oligopoolcalc`

**Record repaired:** `final/oligopoolcalc.md` (edited in place)
**Audits processed:** `citation_audit/oligopoolcalc.md` (28 findings), `factcheck/oligopoolcalc.md` (A×10, B×9, C×1)
**Date:** 2026-08-14

## Verification basis

Every finding was re-checked against primary source before any edit. Sources used:

| Source | How obtained |
|---|---|
| Paper | `papers/Hossain2024_OligopoolCalculator.pdf`, text extracted with PyMuPDF (15 pp., 92,921 chars) |
| Repository | `codeload.github.com/ayaanhossain/oligopool/tar.gz/b88fa394ca67ed4c48ec41127b5707694ee7cf0a` — the exact commit the fact-check pinned; unpacked and read directly |
| PyPI | `https://pypi.org/pypi/oligopool/json` (live fetch) |
| GitHub REST API | repo metadata + commit list for `ayaanhossain/oligopool` and repo metadata + top-level contents for `hsalis/SalisLabCode` (live fetch) |
| oligopool.com | live fetch, HTML stripped to text |

**No capability value was changed.** All edits touch evidence and prose only.

---

## 1. APPLIED

### A. Identity and citation

**1. Wrong paper title.** *(citation audit "Other" #1; fact-check A1)*
- **Correction:** record L6 title replaced with *"Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements"*.
- **Verified how:** PDF page 1 carries that exact title. The repository README's own citation block and BibTeX entry (`README.md` L288–L306) repeat it verbatim with the same authors, journal, year, volume 13, number 12, pages 4218–4232, and DOI `10.1021/acssynbio.4c00661`. The displaced title belongs to the Non-Repetitive Parts Calculator paper — consistent with `NRP_Calculator` being a sibling directory in `hsalis/SalisLabCode`.
- Everything else on L6 (authors, journal, year, volume, pages, DOI, PMID 39641628, PMCID PMC11669329) checked out and was left alone.

**2. "Paper title/keywords include 'massively parallel reporter assay'."** *(citation audit #2)*
- **Correction:** `assay_mpra` evidence now says the phrase is a paper **keyword**, and quotes the actual title fragment.
- **Verified how:** PDF page 1: the title reads *"…Rapid Analysis of Massively Parallel Barcoded Measurements"*; the phrase appears only in the `KEYWORDS:` line.

**3. `hsalis/SalisLabCode` described as a "mirror."** *(citation audit #23)*
- **Correction:** §1 Repo row and §8 row now describe it as a multi-project Salis-lab monorepo containing an `Oligopool_Calculator` subdirectory, explicitly not a mirror, with currency unverified.
- **Verified how:** GitHub API — `description: "Models, design algorithms, and other software related to Salis lab publications"`, `fork: false`; top-level contents are `ELSA_Calculator`, `ModelTestSystem`, `NRP_Calculator`, `Oligopool_Calculator`, `Promoter_Calculator`, `Synthesis_Success_Calculator`, `T7_Promoter_Calculator`, `mRNA_Stability_Calculator`.

### B. Paper attributions

**4. Data Availability "lists only the two GitHub repos."** *(citation audit #3)*
- **Correction:** `interface` evidence now says the Data Availability Statement lists three locations, naming PyPI.
- **Verified how:** PDF Data Availability Statement: *"A Python implementation of Oligopool Calculator can be found at (1) github.com/hsalis/SalisLabCode, (2) github.com/ayaanhossain/oligopool, and (3) pypi.org/project/oligopool."*

**5. GPL quotation attributed to Data Availability.** *(citation audit #4)*
- **Correction:** `license` evidence now attributes the sentence to the end of the Introduction (printed pp. 4219–4220) and says explicitly that it is *not* the Data Availability Statement; the `*Source:*` line changed from "paper p.4220" to "paper pp. 4219–4220".
- **Verified how:** the sentence begins on PDF p. 2 (*"The Oligopool Calculator is available at https://github.com/hsalis/SalisLabCode and"*) and completes at the top of PDF p. 3 (*"https://github.com/ayaanhossain/oligopool under a GPLv3 open-source license."*), immediately before the `RESULTS AND DISCUSSION` heading. The Data Availability Statement sits near the end of the paper and contains no license wording.

**6. 6.8-second benchmark misidentified.** *(citation audit #5; fact-check A6)*
- **Correction:** §5 item 1 now quotes the paper and states the figure times the **spacer** design step.
- **Verified how:** paper, Motifs/Spacers section: *"The design process is very fast; for example, designing the spacer regions for 93,180 oligonucleotides required only 6.8 s (Table S1)."*

### C. Module-surface claims

**7. "Every module signature is `(input_data, ..., output_file) -> (DataFrame, stats)`."** *(citation audit #7; fact-check A3)*
- **Correction:** `library_as_object` evidence now scopes the signature to design/transform modules and names the five stats-only modules.
- **Verified how:** `docs/docs.md` L140–L156 distinguishes three return shapes: design/transform → `(out_df, stats)`; stats-only → `stats` (`background`, `lenstat`, `inspect`, `index`, `pack`); counting → `(counts_df, stats)`. `api.md` L383–L390 shows `op.background(input_data, maximum_repeat_length, output_directory)` returning `stats`; `api.md` L1471–L1478 shows `op.inspect(target, kind, verbose)` returning `stats`. `oligopool/__init__.py` `__api__` lists 22 public names — 20 functions plus the `vectorDB` and `Scry` classes — so the record's "22 modules" count is itself correct and was left alone.

**8. "All 22 modules return materialized DataFrames/CSVs; the package evaluation model is eager throughout."** *(citation audit #8; fact-check A3)*
- **Correction:** `lazy_evaluation` evidence rewritten to state that every module materializes its result in one of three documented return shapes and that the **public API** exposes no deferred evaluation; the grep parenthetical now names the three real `yield` hits; a closing clause notes internal generators and `__getattr__` lazy imports are not user-facing.
- **Verified how:** `grep -niE '\byield' docs/{docs,api,agent-skills}.md README.md` returns exactly three hits — `api.md:955`, `api.md:1257`, `agent-skills.md:247` — all ordinary-language "yield different barcodes" prose, i.e. non-feature hits as the record claimed. `grep -rln '^\s*yield ' oligopool --include=*.py` returns `base/vectordb.py`, `base/core_pack.py`, `base/scry.py`, `base/core_barcode.py`, `base/core_count.py`, `base/utils.py`, `cli.py`. `oligopool/__init__.py` L73 defines `__getattr__` using `importlib.import_module`.
- Value `no` is unaffected: no user-facing lazy or deferred API exists.

**9. "An Oligopool Calculator library is therefore always exactly the rows the user supplied."** *(citation audit #9; fact-check A2)*
- **Correction:** §2 sentence replaced with "No shipped operation can therefore union two independently designed tables into a single pool of rows", plus a clause naming `compress`/`expand` as the only row-count changes, both within one library.
- **Verified how:** `api.md` L812 — `compress` returns `synthesis_df` (`DegenerateID`, `DegenerateSeq`, `Degeneracy`, `OligoLength`), one row per degenerate oligo; `api.md` L818 — *"If sequences are too diverse to compress, returns 1:1 mapping"*, so it is normally fewer rows. `api.md` L873 — *"Expansion can be exponential; use `expansion_limit` for safety"*, so `expand` multiplies rows.
- **The `merge`/`join` quotations and the no-row-union fact itself were re-verified and are intact** (`api.md` L427, L563, L567). `mixed_mutagenesis_one_pool = no` rests on that fact, which survives, so no value moves.

**10. "None of the 22 `op.*` modules generates variants."** *(citation audit #10)*
- **Correction:** narrowed to "no `op.*` module derives variants from a mutation specification".
- **Verified how:** `api.md` L838 — `expand`'s purpose is *"Expand IUPAC-degenerate sequences into all concrete A/T/G/C sequences"*, so it does emit sequences the user did not enumerate row-by-row; but it infers no mutagenesis design. `motif(motif_type='variable')` likewise generates per-row DNA. The narrower wording is what the sources support and is what the value rests on.

**11. "The `ID` column is a required user input on every module."** *(citation audit #11; fact-check A4)*
- **Correction:** `automatic_naming` evidence now says "every table-based design module" and names the counterexamples in a parenthetical.
- **Verified how:** `api.md` L396 — `background` takes *"list of DNA strings, CSV/DataFrame with `Sequence` column, or FASTA path"*; `api.md` L1486 — `inspect` takes `target`, an artifact path; `pack` takes FastQ paths (`api.md` L1008–L1011); `vectorDB`/`Scry` are constructors, not DataFrame modules (`api.md` L1526, L1574).

**12. "`oligo_length_limit` is a required parameter on every design module."** *(citation audit #14; fact-check A5, first half)*
- **Correction:** `synthesis_constraints` *Length:* clause now names the modules that require it and the four Design-inventory modules that do not.
- **Verified how:** `grep -n 'oligo_length_limit' docs/api.md` plus the package signatures. Present in `barcode` (`api.md` L73), `primer` (L154), `motif` (L241), `spacer` (L317), `join` (L539), `pad` (L713), `lenstat` (L1301), `verify` (L1345). Absent from `background`, `merge`, `revcomp`, `final` — confirmed against `oligopool/{background,merge,revcomp,final}.py` signatures.

### D. Mechanism claims

**13. `pad` described as "Type IIS pads for scarless Golden Gate."** *(fact-check A5, second half)*
- **Correction:** both occurrences (`synthesis_constraints` *Length:* and §6 item 4) now say the pads are excised scarlessly before assembly, with the enzyme removing pads and the `split` overlaps doing the joining.
- **Verified how:** `docs/docs.md` L878, verbatim: *"The Type IIS enzyme is for **pad removal**, not fragment-to-fragment ligation - the 15-30 bp overlaps from `split` drive the actual assembly."* The post-synthesis workflow (docs.md L873–L876) is PCR → Type IIS digest → optional mung bean nuclease → *"Assemble the cleaned fragments via their split-designed overlaps (Gibson, overlap-extension PCR, etc.)"*. **"Golden Gate" appears zero times** in `docs/` or `README.md`.

**14. "It caps homopolymer runs only incidentally, at ~k+1."** *(citation audit #15)*
- **Correction:** clause replaced — the constraint bites only by coincidence of corpus membership; a homopolymer run is not intrinsically capped.
- **Verified how:** `oligopool/base/core_barcode.py` L549–L604, `is_nonrepetitive()`: it streams the canonical `(maxreplen+1)`-mer spectrum of `left_context + barcode + right_context` and returns `False` only `if kmer in corpus`, where `corpus = oligorepeats[index]` is the indexed input/context repeat set. Nothing tests run length directly.

**15. Causal over-reading of `verify` motif emergence.** *(citation audit #16)*
- **Correction:** three sites. (a) `per_sequence_provenance`: "the element column that produced it" → "the element column containing its start offset", with a note that a junction-spanning motif is attributed to the column it starts in. (b) `synthesis_constraints`: heading changed to "Motif emergence is checked as excess over a library-wide baseline", the causal "emerge from concatenation" clause replaced with the excess-over-minimum test plus an explicit note that the excess may sit wholly inside one column. (c) §6 item 20: same start-offset wording.
- **Verified how:** `oligopool/base/core_verify.py` L199 — `baselines[motif] = min(counts) if counts else 0`, where `counts = [oligo.count(motif) for oligo in oligo_sequences]`; L210 — flag when `count > baselines[motif]`; L220 — `columns.append(ut.get_boundary_column(pos, boundaries))`, and `utils.get_boundary_column` is documented as *"Map a flat character position to its enclosing column name"*. `api.md` L1439 states the semantics directly: *"**Emergence**: A motif has emergence when its count exceeds the library-wide minimum (baseline)"*.
- Design-time junction awareness (context columns, junction-aware background screening) is a separate documented mechanism and was left untouched. Values unaffected.

**16. `get_primer_order()` "topologically resolves the paired-primer dependency graph."** *(citation audit #20; fact-check A9)*
- **Correction:** both occurrences (§5 item 5, §6 item 8) now describe a one-pass `deque` insertion heuristic and say explicitly it is not a general topological sort.
- **Verified how:** `examples/design-assembly-parser/design_assembly_parser.py` L222–L234 — a single `for` loop over `elements_spec` that appends each primer and inserts an unseen `paired_primer_column` at `order.index(element_name)`. No graph construction, no indegree bookkeeping, no cycle detection.

**17. Analysis Mode credited with index-hopping detection and mismatch-distribution reporting.** *(fact-check A7)*
- **Correction:** §6 item 1 now attributes both to the paper by quotation and states that neither term appears in the current documentation, describing what `xcount`/`acount` actually emit.
- **Verified how:** `grep -rniE 'index.?hop|hopping|mismatch distribution' docs/ README.md` → **zero hits**. `api.md` L1226 — `xcount` returns *"one `<indexname>.ID` column per index, plus `CombinatorialCount`. Missing barcodes shown as `'-'`"*. `api.md` L1105 — `acount` returns `<indexname>.ID`, `BarcodeCount`, `AssociationCount`. The paper does make both claims (index hopping: *"identified during mapping and removed during counting"*; *"reporting of mismatch distributions"*).
- §5 item 3 already attributed index hopping to the paper's ribozyme example, which is where the paper makes the claim — left unchanged and correct.

### E. Counts and file facts

**18. "All 40 `oligopool/**/*.py` files" and `oligopool/core_lenstat.py`.** *(citation audit #12; fact-check A10)*
- **Correction:** 40 → **41** (with the 23/18 split) in the `design_visualization` evidence and its `*Source:*` line; the path citation corrected to `oligopool/base/core_lenstat.py` (2 occurrences).
- **Verified how:** at commit `b88fa394…`, `find oligopool -type f -name '*.py' | wc -l` = 41 (23 at top level, 18 nested); `find . -name 'core_lenstat*'` resolves to `oligopool/base/core_lenstat.py`. The zero-hit plotting-import result the record relies on was re-run and **is correct**: `grep -rniE 'seaborn|matplotlib|pyplot' oligopool --include=*.py` returns nothing.

**19. `primer_guanine_cytosine_content` "6 occurrences".** *(citation audit #13)*
- **Correction:** four sites reworded from "6 occurrences"/"6×" to "6 lines"; the precise site adds "7 occurrences in all because L1245 carries it twice".
- **Verified how:** `grep -c` = 6 lines; `grep -o … | wc -l` = 7 tokens; lines are L542, L1070, L1176, L1245, L1321, L1347 — exactly the six the record listed, with L1245 containing the token twice.

**20. "13 missed capabilities."** *(citation audit #6)*
- **Correction:** 13 → **11**.
- **Verified how:** `extractions/oligopoolcalc.md` §7 contains exactly 11 numbered items; `final` §6 contains 22; items 1–11 of §6 correspond one-to-one to the extraction's 1–11, so 11 were added.

### F. Example-file facts

**21. Parser "over 4,351 real promoter sequences".** *(citation audit #19)*
- **Correction:** §5 item 5 now states the file holds 4,351 lines but `run_design_assembly_parser.py` drops 848 first, so `design_parser` runs on **3,503**.
- **Verified how:** `wc -l promoters.txt` = 4,351 (all non-empty). `run_design_assembly_parser.py` L21–L25 filters `excluded_motifs = ['GGATCC','TCTAGA','GGTCTC','GAGACC','CCCCC','AAAAA','TTTTT','GGGGG']` before `dp.design_parser(pool_size=len(promoter_list), …)`. Recomputing that filter in Python: 3,503 kept, 848 removed. The README's "~4,500" wording was also confirmed (`examples/design-assembly-parser/README.md` L9).

**22. Notebook "re-executed 2026-02-18".** *(citation audit #17)*
- **Correction:** §5 item 4 now says the markdown header reads *"Updated February 18, 2026"* and gives the actual last commit touching the file.
- **Verified how:** notebook cell 0 markdown contains `**Updated** February 18, 2026`. GitHub commits API filtered to `path=examples/OligopoolCalculatorInAction.ipynb` returns `b88fa394ca67`, `2026-02-22T02:46:06Z`, subject `chore: rerun notebook` — which is also repo HEAD. "Updated" is not evidence of an execution date. (Cell count 284 was re-verified and is correct.)

**23. "3–17-bp spacers."** *(citation audit #18)*
- **Correction:** §5 item 4 now says "spacers padding to 170 bp", with a parenthetical giving both the notebook's prose figure and the executed `lenstat` output.
- **Verified how:** notebook prose: *"Variable length spacers ranging from 3 to 17bp to pad the oligos to 170bp."* Executed output in the same notebook: `Element Spacer: Occupies 0 to 14 Base Pair(s)` on a `DataFrame w/ 6,232 Record(s)` at `Oligo Limit: At most 170 Base Pair(s)`. The three extra degenerate bases live in the separate `AatIIPadded` motif (`primer_sequence_constraint='GACGTC'+'NNN'`).

**24. "HGVS-style" variant IDs.** *(citation audit #22)*
- **Correction:** three sites reworded to "ref-position-alt"; §3 and the `hgvs_input` note add that this is not HGVS syntax and that no source calls it that.
- **Verified how:** `mutant_generator.py` L78 — `variant_id = f'{ref_base}{pos + 1}{alt_base}'`. `grep -rniE '\bhgvs\b'` over `docs/`, `README.md`, `oligopool/`, and the example `.py` files returns **zero hits**. HGVS DNA-substitution syntax requires a coordinate-type prefix and `>` (`c.5A>T`). The record's `hgvs_input = no` value is unaffected — it rests on the zero-hit grep and the absence of any parser, both re-confirmed.

### G. Quantitative claim

**25. "6–20×" presented as fact.** *(citation audit #21)*
- **Correction:** only the one site that asserted it as an outcome (`assay_dms` *The limit:* sentence) was hedged to "the documented 6–20×", with a parenthetical noting no shipped executed example demonstrates it.
- **Verified how:** the figure is genuinely documented — `docs.md` L2383 (`# 2. Order synthesis_df for synthesis (6-20x fewer oligos than individual variants)`) and `agent-skills.md` L371. Executed evidence in the repo: the notebook prints `Compression Ratio : 1.2x` (12 → 10) and `Compression ratio: 1.4x` (500 → 360); `run_degenerate_demo.py` prints a table of theoretical outcomes (64×, 4,096×, and no compression for single mutants). The three other sites already said "documented" or quoted the docs, and were left alone.

### H. Hosted-service claim

**26. "No hosted web service exists."** *(citation audit #24)*
- **Correction:** folded into the L195 fix and the §8 table row — both now say "no hosted web service **is claimed** anywhere in the paper, README, or docs".
- **Verified how:** the universal negative is unverifiable; the claim-scoped version is verifiable and was verified — the Data Availability Statement lists only two GitHub repos and PyPI, and no web-service claim appears in `README.md` or `docs/`.

### I. Section-B additions (incomplete)

**27. `verify` skips degenerate columns.** *(fact-check B6)*
- **Correction:** one bullet added to §7, next to the existing `verify` bullet.
- **Verified how:** `api.md` L1438, verbatim: *"IUPAC/degenerate columns are skipped silently during DNA column detection"*.
- **Why applied rather than escalated:** the record presents `verify` as the "Dedicated pre-synthesis QC" underpinning `synthesis_constraints = yes` and separately promotes Degenerate Mode; without this, a reader would reasonably assume compressed output receives the same QC. It does not.

**28. `compress` applicability limits.** *(fact-check B2)*
- **Correction:** one bullet added to §7, next to the existing `compress` terminal-branch bullet.
- **Verified how:** `api.md` L798 (*"Must contain only A/T/G/C (no degenerate codes)"*), L817 (*"Sequences of different lengths are compressed independently"*), L818 (*"If sequences are too diverse to compress, returns 1:1 mapping"*), L805–L806 (`rollout_simulations`, `rollout_horizon`).
- **Why applied rather than escalated:** these restrictions directly qualify the record's own "6–20×" and "a synthesis-cost capability no other surveyed tool has" claims.

---

## 2. REJECTED

**R1. Citation audit, "NOT FOUND IN ANY SOURCE" #1 — the `oligopool.com` quotation.**
- **Finding:** the quoted wording *"Free, professional tools for oligonucleotide design and analysis… peer-reviewed methods, including SantaLucia 1998 and Breslauer 1986"* does not occur on the live page.
- **Why the finding is wrong:** both quoted fragments occur **verbatim** on the live page, and the record joins them with an ellipsis, which correctly signals the elision.
- **Evidence:** live fetch of `https://oligopool.com` (HTTP 200, 82,348 B), HTML stripped to text. Footer: *"OligoPool.com Free, professional tools for oligonucleotide design and analysis. Built for the molecular biology research community."* Body, under "Why Use These Tools → Published Method Basis": *"Tm and related calculations reference peer-reviewed methods, including SantaLucia 1998 and Breslauer 1986, with assumptions documented in the method pages."* The audit found the second fragment but appears to have missed the first, which sits in the page footer.
- The record's surrounding claims were re-verified and are correct: the page advertises **"19 Tools"**, and `Salis`, `Hossain`, `acssynbio`, and `oligopool calculator` each return **zero** hits in both the visible text and the raw HTML.
- **No edit made.**

**R2. Citation audit #25 — "runnable" / "installable and runnable today" uncited.**
- **Finding:** artifact availability does not verify the stronger runtime claim.
- **Why the finding is inaccurate as stated:** the record already carries exactly the qualification the audit asks for, twice, including in the same table. §8 "Install paths" row ends *"**Not attempted — no-install constraint.**"*, and §9 states *"**Nothing was installed or run.** Runtime behaviour, actual benchmark reproduction, and installability on a current toolchain are all unverified."*
- The underlying artifact facts were independently re-verified: PyPI serves `2026.2.22.1` (wheel + sdist, `requires_python: <4,>=3.10`, `license_expression: GPL-3.0-only`), the GitHub repo is live and unarchived, and a `Dockerfile` plus `docs/docker-notes.md` are present at the pinned commit.
- **No edit made.**

**R3. Fact-check section C — balance.**
- The audit's own conclusion is that the record is *"broadly balanced rather than systematically understating the tool"* and that its problems are *"discrete factual overclaims, a few paper/current-release conflations, and missing lower-level modes and restrictions"* — all of which are handled above. No corrective action is called for, and balance judgments are reserved for the authors in any case.
- **No edit made.**

---

## 3. ESCALATED — left untouched, authors to decide

**E1. Fact-check A8 — three paper-vs-current-release divergences are not flagged.**
- All three are real and were verified:
  - Paper primer rule #2: *"Pairs of primers must have matching DNA melting temperatures (e.g., with a difference smaller than 0.3 °C)"*; current contract, `docs.md` L469: *"Tm matching is applied within 1 °C."*
  - Paper: *"pairs of primers and primer binding sites with target length (14 < L < 23)"* (i.e. 15–22 nt); `oligopool/primer.py` L167 validates `minlenval=15` with no maximum — no 22/23 cap exists in the file.
  - Paper: the callback can *"tally the function's results in dictionary format… for example, identifying transcription start sites, RNA cleavage locations, or mutated amino acids"*; current contract, `api.md` L1163: *"bool - True to accept, False to reject"*.
- **Why escalated, not applied:** the record states **none** of these incorrectly — it never quotes 0.3 °C or a primer-length bound, and §6 item 14 already describes the callback correctly as Boolean-only. The finding is therefore an omission, and whether it is material is an authors' call. It would also cut against the framing of the record's own scoping note at §1 ("the shipped software has moved substantially beyond the 2024 paper"), because two of the three are cases where the release is **narrower** than the paper — reframing that note is more than a clause.
- **Question:** should §1's scoping note (or §7) gain a sentence recording that in these three respects the current release differs from, and in two cases is narrower than, the paper?

**E2. Fact-check B1, B3, B4, B5, B7, B8, B9 — lower-level modes and restrictions not covered.**
- Facts verified at source where checkable: spacer length modes (`api.md` L344: `None`=auto-fill, `int`=fixed, `list`=per-row, DataFrame/CSV with `ID` + `SpacerLength`); `pack` modes and filters (`api.md` L1008–L1023: single/paired FastQ incl. `.gz`, read orientation, minimum length and minimum average Phred, `pack_type` concatenate vs merge, dedup, reusable `.oligopool.pack`); `acount` vs `xcount` output semantics (`api.md` L1105, L1226); `index` anchor requirements (`api.md` L953, L956); `vectorDB`/`Scry` direct method surfaces (`api.md` L1526–L1545, L1574–L1600); `background` bounds (`api.md` L397: `maximum_repeat_length` 6–20).
- **Why escalated, not applied:** none of these contradicts or materially qualifies an existing claim in the record — they are depth-of-coverage judgments about a capability-matrix record rather than a manual. Several also have no unambiguous existing entry, so placement would itself be an editorial decision. B7 (`inspect`) and B9 (background input formats) are already partly covered at `design_visualization`, §8, and the `genome_coordinates` evidence.
- **Question:** which, if any, of these seven should be added, and to which entries? Options: (a) add none — the record is a capability record, not an API reference; (b) add only the two that touch survey-relevant behaviour (`acount` vs `xcount` semantics; spacer length modes); (c) add all seven as brief clauses in §6 and §7.

**E3. Citation audit #26 — "A web-search summariser misidentified it."**
- **Why escalated:** this is a claim about the prior session's own process, not about the tool. Nothing in any primary source can confirm or refute it, and the record cites no transcript or URL. Deleting a true-but-uncited note and keeping a false one are both wrong, and I cannot tell which it is.
- The substantive warning it supports — that `oligopool.com` is unrelated and must not be cited — is independently verified and correct (see R1).
- **Question:** keep the misidentification anecdote, drop it and keep only the verified "unrelated site" warning, or keep it with an explicit "process note, not source-backed" marker?

**E4. Citation audit #27 — third-pass provenance and rate-limit narrative.**
- Affects the `**Status:**` line (§0 header), §1's "Sources consulted" GitHub row, §9's "Corrections adopted from the reviewer, re-verified here", and §9's "Third-pass limitation to note".
- **Why escalated:** these are assertions about how the record was produced, not about the tool; no logs are cited and none can be recovered.
- **Relevant new fact for the decision:** the GitHub REST API was fetched successfully during this repair pass and **the disputed figures are confirmed exactly** — `pushed_at 2026-02-22T02:46:27Z`, `open_issues_count 0`, `license.spdx_id GPL-3.0`, `archived: false`. The commit list also confirms the three commits §9 names, to the second: `b88fa394ca67` 02:46:06 `chore: rerun notebook`; `7c056c87d315` 02:43:15 `chore: bump version to v2026.02.22.1`; `d6bdb8ec6a3a` 02:40:15 `chore: update stale v2026.02.18 references to v2026.02.22`. PyPI likewise confirms exactly two releases ever: `2024.12.2` (2024-12-02T18:43:35Z) and `2026.2.22.1` (2026-02-22T02:48:39Z).
- **Question:** the §9 "Third-pass limitation to note" paragraph is now stale — the data it hedges has since been fetched directly. Replace it with a plain statement of the verified figures, delete it, or leave the provenance narrative alone as a record of what the earlier passes did?

---

## 4. Tensions left in place (deliberately not smoothed)

Per the surgical-editing rule, these nearby sentences were left alone even though a corrected fact now sits next to them:

1. **§2 still frames the row-union fact as "unrebuttable" and "cannot be rebutted by pointing at an example script."** Only the over-broad inference was replaced; the no-row-union fact itself was re-verified and stands, so the framing is still defensible — but a reader may notice the rhetoric is louder than the (now narrower) sentence.
2. **`assay_mpra`'s `*Source:*` line still reads "paper title/abstract"** although the MPRA phrase is a keyword. The evidence body above it is now precise; the source pointer was not touched.
3. **`library_as_object`'s `*Source:*` line still reads "all 22 module signatures in api.md."** The evidence body now states the three return shapes; the pointer is a locator, not a claim.
4. **§6 item 22 and §9 still recommend the row-union sentence as "the recommended load-bearing sentence for the referee response."** They describe it correctly as "the absence of any row-union / sub-library concatenation operator", which survives the correction unchanged.
5. **§0 still says "No capability VALUE changed" and "All 24 values survived."** Still true — this pass changed no values.
6. **§0's change table rows describe what the *previous* pass changed.** Only the "13 missed capabilities" count was corrected; the rest were left as a historical record of that pass, not restated for this one.

---

## 5. Summary

| Outcome | Count |
|---|---|
| Applied | 28 corrections across 26 record locations |
| Rejected | 3 (1 citation-audit finding demonstrably wrong; 1 already qualified in the record; 1 balance judgment needing no action) |
| Escalated | 4 questions (A8; B1/B3/B4/B5/B7/B8/B9 as a group; citation audit #26; citation audit #27) |
| Capability values changed | **0** |

## Pass 2 — policy application

**Date:** 2026-08-14

**Verification basis:** the paper was read from `papers/Hossain2024_OligopoolCalculator.pdf` with the required PyMuPDF command; repository documentation and source were checked at commit `b88fa394ca67ed4c48ec41127b5707694ee7cf0a`; live GitHub repository metadata and the commit list were fetched on 2026-08-14; the live PyPI project page/JSON reports release `2026.2.22.1`. The live repository and PyPI version still match the record header, so Policy D requires no version-drift parenthetical.

### Policy A — balance and emphasis

**Outcome: declined-by-policy.** No balance or proportional-emphasis edit was made. Fact-check §C itself concludes that the record is broadly balanced; the authors' ruling independently declines rebalancing.

### E1 — fact-check A8, paper/current-release divergences

**Outcome: declined-by-policy.** No edit.

- **Verified:** the paper gives the paired-primer example as a Tm difference below 0.3 °C, gives `14 < L < 23`, and describes callbacks tallying arbitrary properties. Current `docs/docs.md` documents paired-primer matching within 1 °C; current `oligopool/primer.py` validates `primer_sequence_constraint` with `minlenval=15` and supplies no maximum-length validator; current `docs/api.md` defines the callback return as Boolean (`True` accept, `False` reject).
- **Shipping-row test:** none of the three differences changes any of the 17 locked row scores. They also do not change Table 1's Purpose, Key features, Output, or Availability cells: current primer design remains constraint-aware and the callback remains an analysis-time filter. The record already states the current 1 °C and Boolean contracts correctly.
- **Value-basis check:** current primary evidence continues to support `primer_design = yes`; no value escalation is warranted. Header release `2026.2.22.1` is still current, so this is paper/current behavior, not header/current-release drift.

### E2 — lower-level modes and restrictions

The grouped escalation is resolved item by item under the shipping-row filter.

**B1, spacer length modes — outcome: declined-by-policy.** Verified in `docs/docs.md` / `docs/api.md`: `None` auto-fills to the oligo limit, an integer fixes length, a list supplies per-row lengths, and CSV/DataFrame input uses `ID` + `SpacerLength`; `'-'` marks no spacer. These modes do not change a locked-row score or the Table 1 description "constraint-aware ... spacers."

**B3, read packing modes and filters — outcome: declined-by-policy.** Verified in `docs/docs.md` / `docs/api.md`: single- or paired-end FastQ, read orientation, length and mean-Phred filters, concatenate/merge modes, deduplication, and reusable `.oligopool.pack` artifacts. These preprocessing details do not change a locked-row score or a Table 1 cell.

**B4, `acount` versus `xcount` semantics — outcome: applied.** In final §6 item 1, added the brief output clause that `xcount` retains reads when at least one barcode maps and fills missing positions with `'-'`, while `acount` separately emits `BarcodeCount` and stricter `AssociationCount`. Verified in `docs/docs.md` (`acount`/`xcount` output sketches and notes) and `docs/api.md` return contracts. This qualifies through Table 1's **Output** cell because it changes "variant counts" into a materially more precise account of the two shipped count products.

**B5, index restrictions — outcome: declined-by-policy.** Verified in `docs/api.md`, `oligopool/index.py`, and `oligopool/base/core_index.py`: barcode indices require an anchor side, constant anchors, unique fixed-length barcodes, and associate data/anchors for association indexing. These prerequisites do not change a locked-row score or Table 1 cell.

**B7, artifact inspection — outcome: declined-by-policy.** Verified in `docs/docs.md` / `docs/api.md`: `inspect` is read-only and stats-only, auto-detects background/index/pack artifacts, reports metadata, and returns `Valid`, `Corrupted`, or `Invalid`. This utility does not change a locked-row score or Table 1 cell.

**B8, direct advanced-object interfaces — outcome: declined-by-policy.** Verified in `docs/docs.md` / `docs/api.md`: `vectorDB` exposes add/check/iterate/remove/clear/close/drop operations with U→T normalization, while `Scry` exposes train/prime/predict and fast/sensitive modes. These programmable surfaces do not change a locked-row score or Table 1 cell.

**B9, background bounds and formats — outcome: declined-by-policy.** Verified in `docs/api.md`: `background` accepts DNA iterables, FASTA, or CSV/DataFrame with `Sequence`, and `maximum_repeat_length` is 6–20. These details do not alter the `genome_coordinates = no` basis, any locked-row score, or Table 1's off-target-screening description.

**Value-basis check for E2:** the corrected B4 evidence supports the existing account of Analysis Mode and requires no capability-value change. The six declined omissions likewise do not undermine a stated value.

### E3 — citation audit #26, web-search-summariser anecdote

**Outcome: rejected-unverifiable.** The anecdote has no admissible primary source, transcript, or retained search result. Per Policy C it was removed, without an "unverified" marker, from both occurrences. The non-primary `oligopool.com` description/source row was also removed. The record now retains only the primary-source-supported statement: the paper's Data Availability Statement names the two GitHub repositories and PyPI, and neither the paper, README, nor docs claims a hosted service.

### E4 — citation audit #27, provenance and stale rate-limit hedge

**Outcome: applied.** Unsupported reviewer/independent-pass provenance was removed from the status, source table, and §9; the sibling prior-analysis row and attribution were removed; process wording was recast as primary-source evidence; and the stale rate-limit paragraph was replaced with the freshly fetched facts.

- **Verified:** GitHub REST on 2026-08-14 reports `pushed_at 2026-02-22T02:46:27Z`, `open_issues_count 0`, `archived false`, and `license.spdx_id GPL-3.0`. The latest commit is `b88fa394ca67` at 2026-02-22T02:46:06Z (`chore: rerun notebook`), followed by `7c056c87d315` and `d6bdb8ec6a3a` with the messages/times previously stated. PyPI still lists only `2024.12.2` and `2026.2.22.1`; the latter was uploaded 2026-02-22T02:48:39Z.
- **Definition applied consistently:** "last repository activity" now means the latest code push (`pushed_at`), not issue, comment, or repository-metadata activity; neighboring "last touched" wording was narrowed to "last repository push."
- **Value-basis check:** the refreshed evidence still supports the retained `maintained = yes` account (available, public, unarchived release/source/docs), with the existing six-month inactivity caveat. No value escalation is warranted.

### Resulting neighboring tension

The new B4 output clause sits beside the existing comparative sentence, *"No other tool in this survey closes the design→measurement loop."* This pass verified Oligopool Calculator's output contracts, not that cross-tool superlative; it was left untouched under the surgical-editing rule.

### Row-substitution candidate — report only, no action

**Candidate:** replace the locked **`Codon / amino-acid substitutions`** row (already flagged in the locked-row file as a likely near-uniform row) with **`NGS read preprocessing and barcode/association/combinatorial counting`**. Oligopool Calculator primary evidence is unusually concrete: `pack` filters/merges/deduplicates FastQ reads; `acount` emits barcode and association counts; `xcount` emits single- or multi-index combination counts with partial tuples. Cross-tool discrimination was not re-scored from primary sources in this pass, so this remains a candidate only and no substitution was made.

### Pass 2 counts

Counts below treat the seven E2 findings separately; at the original grouped-escalation level, E2 has a mixed outcome.

| Outcome | Count |
|---|---:|
| Applied | 2 (`E2/B4`, `E4`) |
| Declined-by-policy | 7 (`E1`; `E2/B1`, `B3`, `B5`, `B7`, `B8`, `B9`) |
| Rejected-unverifiable | 1 (`E3`) |
| Still escalated | **0** |
| Additional Policy-A decline (not an open escalation) | 1 |
| Capability values changed | **0** |
