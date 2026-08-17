# Citation-integrity audit — Mutation Maker

Record audited: `final/mutationmaker.md`

Repository baseline: `ra100/Mutation_Maker`, default branch `master`, commit `396c1c0ede7529f3dedf65e821e8c1f20c9a7043` (2026-02-14). Paper baseline: `papers/Hiraga2021_MutationMaker.pdf` (14 pages). Only non-verified items are reported below.

## NOT FOUND IN ANY SOURCE

### 1. `generate_random_SSM_input` is not used by the scaling benchmark

- **Record claim (line 223):** “`generate_random_SSM_input(...)` — randomised site generator used for the scaling benchmarks (`backend/ssm_algorithm_comparison.py`).”
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** Repository-wide search finds the function definition at `backend/tests/test_support.py:309`, but no reference to it in `backend/ssm_algorithm_comparison.py`. That benchmark defines and calls its own `generate_inputs()` at lines 91 and 112. The only call to `generate_random_SSM_input` is in `backend/tests/unit_tests/test_ssm.py`; the enclosing `SsmSolveTestConfigRange` class has `num_of_tests = 0`, so even that test case is disabled. Full-history search likewise finds no version in which the cited benchmark used this helper.

### 2. “No scoring model ... anywhere” is contradicted by the cited repository

- **Record claim (lines 112 and 195):** “No predictive model, scoring model, or model-in-the-loop code anywhere in `backend/mutation_maker/`”; “no predictive/scoring model anywhere.”
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** The repository contains `backend/mutation_maker/primer_scoring.py` and `backend/mutation_maker/translation_scoring.py`, as well as numerous score and penalty functions. `qclm.py` and `qclm_solution.py` import and use `PrimerScoring`; PAS modules import and use `TranslationScoring`. The narrower claim that there is no *sequence-function predictive* model may be supportable, but the cited, unqualified assertion that there is no “scoring model” anywhere is not.

### 3. The paper does not identify Gene X or Gene Y as “proprietary Merck genes”

- **Record claim (line 215):** “sequences NOT disclosed — ‘Gene X’ / ‘Gene Y’ are proprietary Merck genes.”
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** The paper uses the anonymized labels “Gene X” and “Gene Y” and does not publish their sequences, but searches of the full PDF find no statement that the genes are proprietary, confidential, Merck-owned, or withheld for that reason. The repository also supplies no such attribution. “Proprietary Merck genes” is an unsupported explanation presented as fact.

### 4. No XLSX “template” artifact exists

- **Record claim (lines 228 and 230):** “**XLSX mutation-table template** for PAS ...”; item 11 is said to be “concrete and self-contained.”
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** The complete repository contains no `.xlsx` or `.xls` file and no downloadable mutation-table template. `FileUploadMutations.tsx:26,54-58` verifies only that the UI accepts `.xlsx` and requires headers `site`, `mt`, and `mt%`. Thus the accepted format exists, but the claimed template artifact does not.

## Other non-verified evidence

### 5. `PASOligoOutput.mutations` is not `List[PASMutationFormatted]`

- **Record claim (lines 84, 190, and 288):** `PASOligoOutput.mutations` is a `List[PASMutationFormatted]` embedded in every emitted oligo.
- **Status:** **wrong-location**
- **Finding:** `backend/mutation_maker/pas_types.py:168-173` declares `PASOligoOutput.mutations = ListProperty(required=False)`, with no element type. `pas_output.py:100-129` populates it with integer indices (`mutations.append(index)`), not `PASMutationFormatted` objects. The typed `ListProperty(PASMutationFormatted, ...)` is instead on the enclosing `PASResult` at `pas_types.py:195`. The integer indices can be resolved against that enclosing list, but the cited per-oligo schema is at the wrong object level.

### 6. `reds` and `blues` are not codon-index lists

- **Record claim (line 84):** `reds`/`blues` are “the mutated vs. wild-type codon index lists.”
- **Status:** **misquoted**
- **Finding:** `pas_output.py:110-119` computes `mut_rel_position = 3 * (mutation.position - 1) + goi_offset - pas_fragment.get_start()` and appends that value to `reds` or `blues`. These are nucleotide start offsets relative to the fragment/oligo, used for display colouring—not codon indices.

### 7. The bio.tools record is not a metadata-free placeholder

- **Record claim (line 25):** `https://bio.tools/mutation_maker` is a “placeholder, no usable metadata.”
- **Status:** **misquoted**
- **Finding:** The URL is live. Its bio.tools API record contains a description, homepage, tool type, topics, functions/operations, Python and JavaScript language entries, GPL-3.0 license, the paper DOI/publication, credits, and dates. Some metadata is stale or questionable (including the old Merck homepage and a CLI label), but “placeholder, no usable metadata” does not describe the page.

### 8. The active REST server is FastAPI, not the cited Hug/Falcon module

- **Record claims (lines 33, 57, 162-168, and 278):** the REST API is defined by `api/api.py`, runs on Hug/Falcon, while `server_fastapi.py` merely “also exists alongside” it.
- **Status:** **wrong-location**
- **Finding:** On the audited default branch, `api/server.py:19` imports `app` from `server_fastapi`; the Makefile and `api/Dockerfile` launch `gunicorn server:app`. `api/server_fastapi.py:25,36` constructs the active FastAPI app and defines the routes. Git history includes the migration commit titled “Migrate API from hug to FastAPI.” `api/api.py` is a retained legacy implementation, so it is the wrong source for describing the current server. The endpoint inventory and missing export-task observations remain substantially applicable to the FastAPI implementation.

### 9. Configurable Tm models are implemented with Biopython, not primer3-py

- **Record claim (line 147):** “Tm via `primer3-py` with configurable model (`TemperatureConfig`: SantaLucia 1997/2004, Breslauer 1986, Sugimoto 1996; ...).”
- **Status:** **minor-discrepancy**
- **Finding:** `temperature_calculator.py:21-29,238-299` implements the selectable Wallace, GC, nearest-neighbour tables, and salt corrections with `Bio.SeqUtils.MeltingTemp`. Primer3 functions are used for hairpin/homodimer/heterodimer calculations and for the separate fixed `NEB_like` Tm calculation (`calcTm` with fixed SantaLucia/Owczarzy settings). The thermodynamic options exist, but the library attribution in the citation is wrong for the configurable models.

### 10. Not every source file carries the quoted GPL header

- **Record claim (lines 171-172):** “every source file carries” the quoted Merck/GPL header, citing `frontend/src/**` among other paths.
- **Status:** **misquoted**
- **Finding:** The repository license and quoted header are real, but the universal per-file claim is false. Counterexamples with no copyright/GPL header include `frontend/src/scenes/QCLM/components/QCLMForm.tsx`, `frontend/src/scenes/PAS/components/PASForm/index.tsx`, `frontend/src/scenes/SSM/components/SSMForm.tsx`, and the empty `backend/mutation_maker/__init__.py`. This finding concerns only the evidence claim, not the repository's GPL-3.0 license.

### 11. The cited FeatureViewer glob resolves to no files

- **Record citation (line 99):** `frontend/src/scenes/*/components/*FeatureViewer.tsx`.
- **Status:** **wrong-location**
- **Finding:** No file matches that path. The three cited components exist one level deeper, at `frontend/src/scenes/SSM/components/SSMResult/components/SSMFeatureViewer.tsx`, `frontend/src/scenes/QCLM/components/QCLMResult/components/QCLMFeatureViewer.tsx`, and `frontend/src/scenes/PAS/components/PASResult/components/PASFeatureViewer.tsx`. The component names are real; the location given is not.

### 12. Plate/well-aware naming applies to SSM, not to all three workflows

- **Record claims (lines 89-95 and 191):** generated names “are plate/well aware”; `SaveFile/index.tsx` provides “plate-aware naming for SSM/QCLM/PAS, 96-well positions.”
- **Status:** **minor-discrepancy**
- **Finding:** The cited formulas show that only SSM names contain a plate index and `A01`-style well. QCLM names use site index plus an optional letter suffix; PAS names use fragment and oligo indices. All three are automatically named, but only SSM naming is plate/well-aware.

### 13. Table 1's clone-count range is 16–48, not 24–48

- **Record claim (line 216):** “24–48 clones Sanger-sequenced.”
- **Status:** **number-wrong**
- **Finding:** The seven `clones tested` entries in paper Table 1 are 48, 24, 24, 24, 16, 16, and 24. The recomputed range is **16–48**.

### 14. PAS fixture `ind=3` has 13 sites, not 11

- **Record claim (line 225):** `generate_pas_input(ind=3)` contains “... (11 sites).”
- **Status:** **number-wrong**
- **Finding:** `backend/tests/test_support.py:766-808` defines mutation positions 1, 2, 10, 12, 30, 36, 39, 63, 66, 73, 84, 87, and 90: **13 distinct sites**.

### 15. `generate_qclm_input` exposes 12 reachable preset indices, not 10

- **Record claim (line 224):** “`generate_qclm_input(ind=0..9)` — 10 preset MSDM cases.”
- **Status:** **number-wrong**
- **Finding:** `backend/tests/test_support.py:359-520` contains reachable cases for every integer index from 0 through 11, i.e. **12 cases**. There is a second, duplicate `elif ind == 11` at lines 521-545 that is unreachable, but the first index-11 case is valid and is exercised by tests.

### 16. `generate_pas_input` has a sixth defined case

- **Record claim (line 225):** “`generate_pas_input(ind=1..5)`.”
- **Status:** **number-wrong**
- **Finding:** `sample_pas_mutations` defines `ind == 6` at `backend/tests/test_support.py:835-876`, and `generate_pas_input(6)` also receives a defined sequence and configuration. The supported fixture range is 1–6, not 1–5.

### 17. The NEB provenance of the restriction-enzyme list is uncited

- **Record claim (line 246):** `frontend/src/shared/motifs.json` contains “NEB enzyme names.”
- **Status:** **uncited**
- **Finding:** `motifs.json` does contain 162 named restriction enzymes, but neither that file nor the motif code identifies NEB as the source/vendor of the list. Repository-wide `NEB` searches find only unrelated Tm-calculator documentation and organism names. No NEB page or other provenance citation is supplied.

### 18. “No public hosted instance found” has no supporting search record

- **Record claim (line 27):** “No public hosted instance found.”
- **Status:** **uncited**
- **Finding:** No search query, deployment-directory check, web result, or URL is cited for this external negative claim. The self-deployment instructions are real but do not establish the absence of every public hosted instance.

### 19. The “Merck publication page” cannot be checked from the record

- **Record claim (line 25):** “Merck publication page (no repo/webapp link).”
- **Status:** **uncited**
- **Finding:** The record supplies no URL or page title for the alleged Merck publication page. Its existence and the claimed absence of links therefore cannot be mechanically checked from the citation.

### 20. Broad claims that no documentation/tutorial site exists are not fully cited

- **Record claims (lines 23 and 213):** “No docs site”; “No tutorial or vignette site.”
- **Status:** **uncited**
- **Finding:** `has_pages=false` verifies that this GitHub repository has no GitHub Pages site, and the named repository Markdown files can be inspected. Those facts do not by themselves establish that no documentation/tutorial site exists anywhere. No external-search evidence is provided for the broader negative claims.
