# Mutation Maker — repair pass change log

**Record repaired:** `final/mutationmaker.md` (edited in place)
**Audits processed:** `citation_audit/mutationmaker.md` (20 findings), `factcheck/mutationmaker.md` (6 in A, 15 in B, 1 in C)
**Date:** 2026-08-14

## Verification baseline

Every finding below was checked against primary source before any edit. Nothing was
applied on the audits' assertion alone.

- **Repository:** cloned `https://github.com/ra100/Mutation_Maker` to a scratch directory.
  HEAD = `396c1c0ede7529f3dedf65e821e8c1f20c9a7043`, branch `master`, 85 commits — the same
  commit both audits used as their baseline. All file/line claims were re-read from that
  working tree, not from the audits' quotations.
- **Paper:** `papers/Hiraga2021_MutationMaker.pdf` re-extracted with PyMuPDF (14 pp.).
  *(Note: a stale `paper.txt` in the session scratchpad held a different paper — VaLiAnT.
  The extraction was redone to a fresh file before any paper claim was checked.)*
- **Live network:** `gh api` for repository metadata; `curl` for the bio.tools API;
  WebFetch/WebSearch for the Merck publication page and the hosted-instance/docs negatives.

**No capability value was changed.** All 24 values in §2 and §3 are byte-identical to the
pre-repair record; only evidence and prose were touched.

---

## APPLIED (18)

### From the citation-integrity audit

**C1 — `generate_random_SSM_input` mis-cited to the scaling benchmark** (§4 item 6)
Applied. `grep -rn generate_random_SSM_input` over the tree returns exactly three hits:
the definition at `backend/tests/test_support.py:309`, an import at
`backend/tests/unit_tests/test_ssm.py:34`, and a call at `test_ssm.py:498`. There is no
reference in `backend/ssm_algorithm_comparison.py`, which defines and calls its own
`generate_inputs()` at lines 91 and 112. Re-pointed the citation to the test module and
noted that the benchmark uses its own generator.
*Audit sub-claim NOT carried over:* the audit added that "the enclosing
`SsmSolveTestConfigRange` class has `num_of_tests = 0`, so even that test case is
disabled." That is wrong. `num_of_tests = 0` on the base class, but the subclasses
`SsmSolveTestSpecificConfigsRandom2` (line 667) and `...Random3` (line 680) set
`num_of_tests = 1`. The helper *is* exercised. The corrected record says only that the
tests use it, and does not repeat the audit's "disabled" claim.

**C2 — "no scoring model anywhere" over-broad** (`assay_insilico`, §2 and §3)
Applied. `backend/mutation_maker/primer_scoring.py` and `translation_scoring.py` both
exist; `PrimerScoring` is imported by `qclm.py:52` and `qclm_solution.py:36`,
`TranslationScoring` by `generate_oligos.py:23`, `pas_output.py:48`,
`pas_degeneracy_recursion.py:21` and `reverse_translation.py:20`. Narrowed the claim to
"no predictive sequence-function model and no model-in-the-loop code", naming the two
scoring modules as thermodynamic/codon-usage penalty functions.
*Value basis unchanged:* `assay_insilico = no` rests on the absence of sequence-function
prediction, which the correction preserves — the record already argued Fig. 2 is algorithm
benchmarking. Deterministic thermodynamic penalty functions are not predictive models.

**C3 — "Gene X / Gene Y are proprietary Merck genes"** (§4 heading)
Applied. Searched the full extracted paper text for `proprietar|confidential|not
disclos|undisclos|withheld`. The single hit is p. 8 line 838, "several proprietary and
freely available codon [optimization tools]" — about third-party software, not the genes.
"Gene X" occurs 20 times and "Gene Y" 4 times, always as a bare anonymised label. Replaced
the unsupported explanation with what the paper actually does.

**C4 — no XLSX template artifact** (§4 item 11)
Applied. `find` for `*.xlsx`/`*.xls` over the whole tree (excluding `node_modules`) returns
nothing. `FileUploadMutations.tsx:26` sets `FILE_EXTENSIONS = '.xlsx'` and lines 54-58
validate `['site','mt','mt%']`, so the *format* is real and enforced in code, but no
template workbook ships. Changed "template" to "input format" and stated the absence.

**C5 — `PASOligoOutput.mutations` typing/containment** (`per_sequence_provenance`)
Applied. `pas_types.py:170` declares `mutations = ListProperty(required=False)` — untyped.
`pas_output.py:107` populates it with `mutations.append(index)`, an integer index. The
typed `ListProperty(PASMutationFormatted, required=False)` is at `pas_types.py:195`, on the
enclosing `PASResult`. Rewrote the clause to describe indices resolving into the enclosing
list.
*Value basis unchanged:* the per-sequence record still exists for every emitted oligo, by
reference rather than by embedding, so both the "kept at partial" grounds and the "lenient
reading = yes" flag survive intact.

**C6 — `reds`/`blues` are not codon indices** (`per_sequence_provenance`)
Applied, in the same sentence as C5. `pas_output.py:110-114` computes
`mut_rel_position = 3 * (mutation.position - 1) + goi_offset - pas_fragment.get_start()`
and appends that. `PASOligoOutput.make_reverse_complement` reindexes with
`abs(len(sequence) - ind - 3)`, which only makes sense for nucleotide offsets. Corrected to
"nucleotide start offsets within the oligo", with the formula, and kept the record's
(correct) reds = mutated / blues = wild-type mapping.

**C7 — bio.tools called a metadata-free placeholder** (§0 sources)
Applied. `curl https://bio.tools/api/tool/mutation_maker/?format=json` → HTTP 200 with
`description`, `toolType`, five EDAM `topic` entries, a `function` block with three
operations, `language: ["Python","JavaScript"]`, `license: "GPL-3.0"`, `publication`,
`credit`, and dates. Replaced "placeholder, no usable metadata" with an accurate
description that keeps the record's skepticism where it is warranted: the `homepage` is the
dead `Merck/Mutation_Maker` URL, `toolType` is `["Command-line tool"]` (wrong — the record
correctly states elsewhere that there is no CLI), and the linked publication DOI is
`10.1101/2020.06.26.171819`, the 2020 bioRxiv preprint, not the ACS paper.

**C8 / factcheck A1 — the active server is FastAPI, not hug/falcon** (`interface`, §7)
Applied. `api/server.py:19` is `from server_fastapi import app`; `api/Dockerfile:20` is
`CMD gunicorn server:app --bind=0.0.0.0:8000`; `Makefile:42` is `cd api && gunicorn
server:app`. `git log --all` shows the commit `2d586c0 ✨ Migrate API from hug to FastAPI`.
Corrected the framework attribution in the `interface` entry and in §7 item 3.
*Scope control:* rather than sweep-edit the six other `api/api.py` citations in the record,
I added one sentence in the `interface` entry stating that `api/api.py` is retained legacy
and that its route inventory is identical to the FastAPI one — which I verified by
comparing both files. This re-frames the remaining citations without generating six more
passages of unaudited text. See "Residual tensions" below.
I did re-point two citations where the file identity is load-bearing for a
*does-not-exist* claim, and verified each against `server_fastapi.py` directly:
- `interface` item 1 (`/v1/export_pas`): `server_fastapi.py` defines only
  `/v1/export_qclm/{task_id}.xlsx` (line 120), `/v1/export_ssm/...` (134),
  `/v1/export_species_table/...` (148). `"export_pas"` appears only at line 65 as the
  `export_name` argument that `get_urls` (line 43) interpolates into a URL. Confirmed.
- §6 "Dead server-side export endpoints": those three routes call
  `CELERY.send_task("tasks.export_qclm" / "tasks.export_ssm" /
  "tasks.export_species_table")`, and `backend/tasks.py` registers only `tasks.ssm` (42),
  `tasks.qclm` (55), `tasks.species_table` (62), `tasks.pas` (68). Confirmed dead.

**C9 / factcheck A4 — configurable Tm models are Biopython, not primer3-py** (`primer_design`)
Applied. `temperature_calculator.py:21` imports `Bio.SeqUtils.MeltingTemp`; line 23 imports
`calcHairpin, calcHomodimer, calcTm, calcHeterodimer` from `primer3`. The
`TemperatureConfig.create_calculator` dispatch (lines 252-284) builds Wallace, GC and NN
calculators from `MeltingTemp.Tm_Wallace` / `Tm_GC` / `Tm_NN`, and `get_nn_table` maps to
`MeltingTemp.DNA_NN1..NN4`. Only `create_neb_calculator` reaches primer3's `calcTm`.
Corrected the attribution and added a one-clause note that the paper (p. 366: "Melting
temperature was computed using Primer3 Python library") describes the 2021 implementation.

**C10 — "every source file carries" the GPL header** (`license`)
Applied. Of the 72 `.ts`/`.tsx` files under `frontend/src`, 66 contain "GNU General Public
License". The six without include `SSMForm.tsx`, `QCLMForm.tsx` and
`PASForm/index.tsx`; `backend/mutation_maker/__init__.py` is 0 bytes. Changed "every" to
"nearly every" and stated the counterexamples.
*Value basis unchanged:* `license = yes` rests on the repository LICENSE and the GitHub
API, not on header universality.

**C11 — FeatureViewer glob matches no files** (`design_visualization` source line)
Applied. `ls frontend/src/scenes/*/components/*FeatureViewer.tsx` → no such file. The three
components are one level deeper, under `*/components/*Result/components/`. Corrected the
glob to `frontend/src/scenes/{SSM,QCLM,PAS}/components/*Result/components/*FeatureViewer.tsx`,
which does resolve.

**C12 — plate/well-aware naming is SSM-only** (`automatic_naming`, §3)
Applied. In `SaveFile/index.tsx`, `PLATE_SIZE` (35) and `getPosition` (48) appear only in
`createSSMTable` (106, 135-142). QCLM names are `${oligoPrefix}-${siteIndex+1}${letter}`
(183) and PAS names are `${oligoPrefix}-Fr${index+1}-${oIndex+1}` (301) — no plate index,
no well. The record's own quoted formulas already showed this; only the summarising
sentence overstated it. Narrowed the lead sentence and the §3 basis cell.
*Value basis unchanged:* all three workflows are still automatically named, which is what
`automatic_naming = yes` asserts.

**C13 — Table 1 clone range** (§4 item 1)
Applied. Paper Table 1 `clones tested` column, read directly from the PDF text: 48, 24, 24,
24, 16, 16, 24. Range is 16–48. Changed "24–48" to "16–48". (The record's other Table 1
figures — 7 libraries, 31–96 target sites, pAVE011-Kan / pET30a(+), W3110 / BL21(DE3),
79.9 % mean — all check out against the same table and were left alone.)

**C14 — PAS fixture `ind=3` site count** (§4 item 8)
Applied. `sample_pas_mutations(ind=3)` at `test_support.py:766-808` defines
`PASMutationSite` entries at positions 1, 2, 10, 12, 30, 36, 39, 63, 66, 73, 84, 87, 90 —
13 sites. Changed "(11 sites)" to "(13 sites)".

**C15 — `generate_qclm_input` preset count** (§4 item 7)
Applied. `generate_qclm_input` (`test_support.py:359`) has reachable branches
`if ind == 0` through `elif ind == 11` (lines 368-520), i.e. 12 cases, plus a second,
unreachable `elif ind == 11` at line 521. Changed "`ind=0..9`  — 10 preset MSDM cases" to
"`ind=0..11` — 12 preset MSDM cases". The three examples the record cites (`ind=2`,
`ind=3`, `ind=4`) are unaffected and were left as written.

**C16 — `generate_pas_input` range** (§4 item 8)
Applied. `sample_pas_mutations` defines `elif ind == 6` at line 835 and
`sample_pas_sequences` defines `elif ind == 6` at line 909. `sample_pas_config` has no
`ind == 6` branch but ends in an unconditional `return PASConfig(...)` fallthrough (line
~984), so `generate_pas_input(6)` is fully constructible. Changed "`ind=1..5`" to
"`ind=1..6`". The record's note that `ind=5` uses a protein sequence is confirmed by
`generate_pas_input` line 704 (`if ind == 5: is_dna_seq = False`) and was left alone.

**C17 — NEB provenance of the enzyme list** (§5 item 11)
Applied. `frontend/src/shared/motifs.json` is a bare JSON array of 162 strings with no
attribution; `backend/mutation_maker/motifs.py` hardcodes the name→pattern dict with no
source comment. Repository-wide `NEB` hits are only the `NEB_like` Tm calculator and its
link to `tmcalculator.neb.com`, plus unrelated `package-lock.json` noise. Replaced "NEB
enzyme names" with "162 restriction-enzyme names" and noted the absent provenance.

**C19 — Merck publication page not identifiable** (§0 sources)
Applied. Located the page at
`https://www.merck.com/publication/mutation-maker-an-open-source-oligo-design-platform-for-protein-engineering/`
and fetched it: it carries no repository, download or web-application link, so the record's
parenthetical is correct. Added the URL so the claim is now checkable. No claim changed.

### From the fact-check audit

**A2 — PDF export and share links presented as current** (`design_visualization`, §3, §5 item 13)
Applied. The paper does say it (p. 366, "parameters can be exported to PDF or as a
temporary sharable link"), so the record's source line was right — but the feature is gone.
Searched `frontend/src` and `frontend/package.json` for `jspdf|html2pdf|pdfmake|\.pdf|
topdf|exportpdf|generatepdf` → zero hits; searched `frontend/src`, `backend` and `api` for
`shareab|sharable|share_link|shareLink|copyLink|permalink` → zero hits. What exists is
print-oriented CSS (`App.css:85,97`) and `Print-Only` classes in the three result scenes.
Attributed the claim to the paper and stated its absence from the current code; removed the
two derived mentions in §3 and §5 item 13.
*Value basis unchanged:* `design_visualization = yes` rests on the neXtProt Feature Viewer
fork, which is present (`frontend/feature-viewer/` plus the three scene components).

**A3 — GC clamp listed as an active constraint** (`primer_design`)
Applied. `gc_clamp` is declared at `ssm_types.py:158` and `qclm_types.py:48`, but
`grep` shows no solver reads `config.gc_clamp`; the only assignment in a solver path is
`ssm.py:603`, `primer3_config.gc_clamp(0)`, which *disables* Primer3's GC-clamp
requirement, matching `PRIMER_GC_CLAMP=0` in `lambda/test.conf` and
`lambda/deploy-package/config_string`. `Primer.get_gc_clamp()` (`primer.py:119`) is called
only from `backend/tests/` and `backend/features/`. The paper's Discussion says Mutation
Maker "by design does not optimize for the presence of GC clamps" — which the record
already quotes verbatim in §6, so the list item also contradicted the record itself.
Removed "GC clamp" from the constraints list.
*Value basis unchanged:* `primer_design = yes` rests on ~13 solver modules and a long
constraint set; one inert field does not bear on it.

**B2 — SSSM fast/Primer3 mode interaction** (`primer_design`)
Applied as a clause, because it corrects an existing imprecise sentence rather than adding
a new topic. `backend/tasks.py:47-50` selects `primer3` or `NullPrimerGenerator()` from
`input.config.use_primer3`, independent of the fast-approximation flag; the UI
(`ParametersFormSection/index.tsx:66-77`) makes them mutually exclusive —
`usePrimerGrowingAlgorithmToggle` sets `usePrimer3` false, `handlePrimer3Change` sets
`usePrimerGrowingAlgorithm` false, and the Primer3 switch is `disabled={usePrimerGrowingAlgorithm}`.
Changed "brute force with Primer3" to "a brute-force search, with Primer3 as an optional
candidate generator via `SSMConfig.use_primer3` … the UI makes the two mutually exclusive".

---

## REJECTED (4)

**factcheck A6 — "the license identifier is too narrow"**
Rejected: the record is not wrong as written. It attributes "GPL-3.0" explicitly to the
GitHub API, and `gh api repos/ra100/Mutation_Maker` really does return
`"license":{"spdx_id":"GPL-3.0"}` — I ran it. The record then quotes the or-later grant
verbatim in the very same sentence, so it nowhere asserts a GPL-3.0-only licence. The
audit's preferred `GPL-3.0-or-later` is a more precise SPDX expression, and it is attested
in the repository (`frontend/package.json:5` `"license": "GPL-3.0-or-later"`; `README.md:326`
"licensed under [GPLv3+]"), but choosing between an accurately-cited registry identifier
and the stricter expression of the grant is a labelling preference, not a correction of a
falsehood. Applying it would have rewritten a true, sourced sentence. Recorded here so the
authors can make the swap deliberately if they prefer the precise expression.

**C18 — "No public hosted instance found" is uncited**
Rejected: this is a complaint about the record not showing a search log, not a claim that
the record is wrong — and the claim is right. I searched for a hosted Mutation Maker
instance and found only the paper, the bioRxiv preprint, bio.tools, the Merck publication
page and the two GitHub repositories; no public deployment. The record's wording is already
hedged ("found"). Adding a search transcript would be fresh text with no factual gain.

**C20 — "No docs site" / "No tutorial or vignette site" are uncited**
Rejected: same reasoning, and the claims verify. `gh api` confirms `has_pages=false`; my own
search surfaced no documentation or tutorial site anywhere. The fact-check's Phase 1 notes
that the FastAPI app serves autogenerated `//docs` when deployed, but that is per-deployment
Swagger, not a documentation site, and its workflow bodies are typed `data: Any` — it does
not contradict the record.

**C1 sub-claim — "even that test case is disabled"**
Rejected and deliberately kept out of the record: `num_of_tests = 0` holds for the base
`SsmSolveTestConfigRange`, but `SsmSolveTestSpecificConfigsRandom2` and `...Random3` set
`num_of_tests = 1`, so `generate_random_SSM_input` is genuinely exercised. Applying the
audit's wording would have introduced a new error. (The rest of C1 was applied — see above.)

---

## ESCALATED (5)

These need an author decision. The record is untouched at every point below.

1. **QCLM custom-organism wiring looks broken in current code.** Verified: the frontend
   sends `codon_usage` and `taxonomy_id` (`frontend/src/services/api.ts:167-168`), while
   `QCLMConfig` declares only `organism = StringProperty(default="e-coli")`
   (`qclm_types.py:79`) and no `codon_usage`/`taxonomy_id` property. Does this belong in
   §6, and does it qualify `codon_optimization` (whose evidence is PAS-centric) or the
   "35,799 species" claim in §5 item 7? I did not edit, because the answer decides what a
   `yes` rests on.

2. **Three further current-code wiring defects in the PAS path.** All three are the kind of
   evidence that could qualify `synthesis_constraints = yes`, whose stated basis is
   "Verified directly in `PASConfig`": (a) the generated-oligo motif check tests the
   original `dna` string rather than each mutated oligo, so a mutation that *creates* a
   forbidden motif may not be rejected; (b) the frontend sends `compute_hairpin_homodimer`
   but `PASConfig` declares no such property; (c) backend `PASConfig` defaults
   (40/56/90 nt) differ from the web-form defaults (40/40/60), so "default settings" means
   different things by interface. Whether these are material enough to state — and whether
   they touch the value — is an authors' call.

3. **The remaining section-B omissions.** Eleven findings that add facts without
   contradicting anything the record says: circular-plasmid wrap-around and the
   exactly-once uniqueness requirement (B1); QCLM site grouping and multi-site primers
   (B4); QCLM non-overlap vs. staged-overlap modes (B5); QCLM degeneracy one-minute timeout
   with silent fallback to random non-degenerate codons (B7); PAS alternating-strand tiling
   with mutations excluded from overlaps (B8); PAS running with zero mutations (B9); PAS
   degeneracy merging only equal-probability alternatives (B10); concrete solver bounds —
   ten-minute reverse-translation timeout, 600-non-improvement stop, five seconds per Tm
   attempt (B13); and the FastAPI `/docs` schemas being typed `Any` (B15). Each is a
   materiality judgment about an omission, which the brief reserves for the authors. I
   verified B10's condition and B2 (applied above) but did not verify the rest in depth,
   since none of them corrects an existing statement.

4. **The SSM export plate ceiling and first-codon display (B3).** Verified:
   `createSSMTable` opens "Second" sheets only at `index === PLATE_SIZE` (96) with no third
   transition, so beyond 192 targets `index % PLATE_SIZE` repeats well coordinates on the
   second sheet; and `ssm.py:780` uses `degenerate_codon.split(",")[0]`, with the code
   comment "The actual degenerate outputs for the export are created on the frontend. This
   is simply for the resulting table, which can only show one primer at a time." This is a
   real defect in exactly the component `automatic_naming` cites, and it also bears on §7
   item 2 (job-wide codon *set*). I edited that entry only for C12 and left the defect
   out, because adding it is a materiality call on a `yes` row.

5. **Fact-check section C (balance).** The audit's overall judgement — "breadth is fair,
   but present-day capability/readiness is modestly overstated" — is explicitly a balance
   and emphasis call. Several of its specific complaints (paper-only PDF/share features, the
   inert GC-clamp field) are now fixed above; whether the residue warrants rebalancing the
   record is for the authors.

---

## Residual tensions (local roughness left deliberately unfixed)

- **`api/api.py` citations elsewhere in the record.** §1, `dag_chaining`, `lazy_evaluation`
  and `vcf_vep_output` still cite `api/api.py`. I verified that every claim they make holds
  equally of `api/server_fastapi.py` (same four workflow routes, same lifecycle routes, same
  three export routes, no composition primitive), so none is false — but they now name the
  superseded file. The `interface` entry carries one sentence saying so. Sweep-editing six
  more citations would have generated more unaudited text than the imprecision costs; a
  future pass may want to do it uniformly.
- **§4 line "Items 5–8 and 11 are concrete and self-contained."** Item 11 is now described
  as a code-defined format with no shipped template (C4). "Concrete and self-contained" is
  still defensible — the column contract is enforced in `FileUploadMutations.tsx` — but the
  sentence was written when item 11 was believed to be a shipped artifact. Left alone.
- **`interface` heading still reads "(two factual errors from the extraction corrected)".**
  That count refers to the earlier review pass, not this one, and this pass corrected a
  third thing in that entry (the server framework). Left alone: it is a historical note
  about a prior pass, and renumbering it is an authors' bookkeeping decision.
- **§8 "Unresolved disagreements and confidence"** still says the reviewer independently
  re-ran the Block C code searches and lists `primer_design`, `interface`,
  `automatic_naming` and `design_visualization` under "High confidence". Three of those
  four entries had evidence corrected in this pass (though no value moved). Left alone; it
  is a confidence statement about values, not about the corrected evidence.

---

## Pass 2 — policy application

**Baseline.** Every finding below was rechecked against Hiraga et al. 2021 and official
`ra100/Mutation_Maker` `master` at
`396c1c0ede7529f3dedf65e821e8c1f20c9a7043`. Read-only `git ls-remote` verified that SHA
was still the branch head on 2026-08-14. Counts use the 19 atomic policy findings below:
**6 applied, 13 declined-by-policy, 0 rejected-unverifiable, 0 escalated**. No capability
value changed.

1. **Open item 1 — QCLM custom-organism wiring: applied.** `api.ts:145-175` sends
   `codon_usage` and `taxonomy_id`; `QCLMConfig` instead declares `organism` with an
   E. coli default; `qclm.py:196` calls `Degeneracy(config)` without an organism. The
   `codon_optimization = yes` value nevertheless verifies independently in PAS:
   `pas.py:160-166` passes `PASConfig.organism` into `Translator`, mutates the GOI to the
   reverse-translated DNA, and returns it in `PASOutput.input_data`. **Edit:** the
   `codon_optimization` entry now says explicitly that its value rests on PAS and records
   the QCLM wiring limit.

2. **Open item 2 — three PAS-path defects: applied (three findings).** All verify in the
   current source: `generate_oligos.py:323-325` checks motifs against input `dna`, not the
   generated `PASOligo.sequence`; `api.ts:297` sends `compute_hairpin_homodimer`, absent
   from `PASConfig`; and `PASConfig`'s 40/56/90 nt defaults differ from the PAS form's
   40/40/60. The value remains supported by active fragment length, overlap Tm/length, GC,
   and structural scoring. **Edit:** one caveat sentence and exact sources were added to
   `synthesis_constraints`; the value was not changed.

3. **Open item 3 — remaining section-B omissions: mixed.**

   - **B1 circular wrap-around and exactly-once checks — declined-by-policy.** Verified at
     `ssm_types.py:40-84` (including both wrap branches and `occurrences != 1`) and the
     form's 60-nt primer cap. These restrictions do not change the locked `Saturation
     mutagenesis` value (`partial`) under its operational test and do not alter a Table 1
     cell.
   - **B4 QCLM grouping/multi-site primers — applied.** `QCLMInput.parse_mutations`
     groups same-position substitutions; `qclm.py:348-355,398-413,621-638` groups nearby
     sites and emits all mutations carried by a primer. The paper's Figure 1B likewise says
     primers cover one or more sites. This is direct support for the locked `Pairwise and
     higher-order variants = yes` value. **Edit:** a brief clause was added to the existing
     QCLM entry.
   - **B5 non-overlap versus staged-overlap modes — declined-by-policy.** Verified in
     `qclm_types.py:52-53`, `qclm.py:358-369,640-665`, and paper p. 361: the non-overlap
     option supports simultaneous annealing; default overlapping neighbours require
     successive MSDM rounds. It is an experimental-layout option, not a composable design
     operation, so no locked score or Table 1 cell changes.
   - **B7 one-minute degeneracy timeout and fallback — declined-by-policy.** Verified at
     `qclm.py:315-342`. The fallback is random non-degenerate codons, but codon-aware
     substitution and a choice of degenerate/non-degenerate policy remain supported, so
     the `Codon / amino-acid substitutions` score would not change.
   - **B8 alternating PAS tiling with mutations outside overlaps — declined-by-policy.**
     Verified by paper Figure 1C and PAS Results plus `pas_output.py:323-328` (odd-indexed
     fragments are reverse-complemented) and fragment-boundary constraints. The existing
     Table 1 output “an optimised sequence or a primer set” remains accurate.
   - **B9 zero-mutation PAS — declined-by-policy.** The backend path verifies:
     `PASInput.mutations` accepts an explicit empty list, `pas.py:197-203` handles no
     mutation protofragments, and PAS can still tile the gene; protein input is
     reverse-translated first. The current GUI, however, requires at least one mutation
     (`InputMutations/index.tsx:56-59,183-188`). This source/API-only mode changes neither
     a locked score nor the already accurate grouped Table 1 purpose/output cells.
   - **B10 equal-probability condition for PAS degeneracy — declined-by-policy.** Verified
     at `generate_oligos.py:75-119,289-307`. It limits when alternatives can be compressed
     together but does not remove codon-aware substitution or the replacement-policy
     choice, so the relevant locked score does not change.
   - **B13 reverse-translation ten-minute timeout — declined-by-policy.** Verified at
     `reverse_translation.py:71-74`; it changes no locked score or Table 1 cell.
   - **B13 600-non-improvement stop — declined-by-policy.** Verified at
     `pas.py:160-166` and `reverse_translation.py:37,87-88`; same policy result.
   - **B13 five seconds per fragment-Tm attempt — declined-by-policy.** Verified at
     `pas.py:32-33,254-274`; same policy result.
   - **B15 FastAPI schemas typed `Any` — declined-by-policy.** Verified at
     `server_fastapi.py:36-65`. Autogenerated `/docs` does not expose workflow schemas,
     but the REST service still exists, so neither the locked rows nor Table 1 Availability
     changes.

4. **Open item 4 — SSM first-codon display and plate ceiling: declined-by-policy.** Both
   verify: `ssm.py:775-803` substitutes only the first comma-separated codon in the result
   display, while `SaveFile/index.tsx:92-149` expands every codon but creates a second
   sheet only at index 96. Above 192 targets, well coordinates repeat on those second
   sheets. Generated identifiers nevertheless remain informative and unique because their
   plate number is `floor(index/96)+1`; the XLSX defect therefore does not change
   `Automatic naming` (`yes`), and the first-codon display does not change `Sequence
   styling` (`yes`). No edit.

5. **Open item 5 / finding A — balance and emphasis: declined-by-policy.** The authors'
   ruling expressly declines rebalancing. No edit.

6. **Finding C — provenance: applied.** Removed the sibling prior-analysis row, third-party
   mirror, bio.tools aggregator, and Merck web-page row. Replaced broad “no docs site” and
   “no public hosted instance found” claims with repository-verifiable facts, removed the
   referee anecdote, and defined both “current” and the SSM/SSSM and QCLM/MSDM aliases.
   The record now relies on the paper, official repository/docs, and PyPI checks only.

7. **Finding D — version drift: declined-by-policy (no edit required).** The header's
   repository snapshot still equals the official branch head; there is no newer tag or
   GitHub release and thus no materially different current release requiring the single
   allowed parenthetical. The scoping note was not renumbered or restructured.

**Escalations:** none. **Rejections:** none; every audit fact was independently verifiable,
including B5's “successive rounds” wording in the paper and B9's GUI/backend distinction.

**Row-substitution candidate (reported, not applied).** Replace the already flagged,
likely-uniform **`Codon / amino-acid substitutions`** row with **`Degenerate-codon
compression`**. Mutation Maker has direct primary evidence: the paper describes heuristic
degenerate-codon generation and benchmarks it against CodonGenie (Fig. 2E), while current
source implements SSSM set cover in `frontend/src/scenes/SSM/components/amino.ts` and
MSDM/PAS set cover in `qclm.py`, `generate_oligos.py`, and
`pas_degeneracy_recursion.py`. This is a better Mutation Maker discriminator than generic
codon-aware substitution; cross-tool values still need verification before any table
substitution.

**Neighbouring tension left intact:** the earlier `## ESCALATED (5)` section is the
historical Pass 1 record and still says those points were untouched. This Pass 2 section
supersedes its outcomes; the earlier section was not rewritten so the audit trail remains
legible.
