# Mutation Maker — evidence memo

**Tool:** Mutation Maker (Merck / MSD)
**Paper:** Hiraga K, Mejzlik P, ... Bitton DA. "Mutation Maker, An Open Source Oligo Design Platform for Protein Engineering." *ACS Synth Biol* 2021;10(2):357–370. doi:10.1021/acssynbio.0c00542 (CC-BY; PMID 33433999)
**Local PDF:** `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/papers/Hiraga2021_MutationMaker.pdf`
**Compiled:** 2026-08-10

---

## 0. Sources consulted

| Kind | Reference |
|---|---|
| pdf | `papers/Hiraga2021_MutationMaker.pdf` (14 pp., extracted with PyMuPDF) |
| prior_analysis | `prior_analyses/Hiraga2021_MutationMaker_analysis.md` |
| repo (dead) | https://github.com/Merck/Mutation_Maker — **HTTP 404** as of 2026-08-10 (both WebFetch and `gh api repos/Merck/Mutation_Maker` return "Not Found"). This is the URL printed in the paper. |
| repo (live) | https://github.com/ra100/Mutation_Maker — surviving copy / network root; retains the original Merck commit history ("Merge pull request #9 from Merck/feature/add-readme-for-backend-code…", commits by Merck authors *Martin Špale* and *prihoda* = paper authors Martin Spale, David Prihoda). GPL-3.0, 19 stars, 13 forks, last commit **2026-02-14**, CI green 2026-02-14. |
| repo (mirror) | https://github.com/ljm176/MutationMaker — description literally "Forked from Merck", last push 2023-01-16 |
| docs | Repo `README.md`, `backend/README.md`, `AGENTS.md`, `e2e-tests/README.md` (no separate docs site, no GitHub Pages; `has_pages=false`) |
| pypi | No PyPI package: `mutation-maker`, `mutation_maker`, `mutationmaker` all return 404 from pypi.org |
| web | bio.tools entry https://bio.tools/mutation_maker (dev registry served a placeholder, no usable metadata); Merck publication page https://www.merck.com/publication/mutation-maker-an-open-source-oligo-design-platform-for-protein-engineering/ (no repo/webapp link) |

**Availability headline:** the repository cited *in the paper* is gone. The code survives only in a community-maintained copy (`ra100/Mutation_Maker`), which is in fact actively maintained (React 18 / antd 4 migration, Python 3.11, Biopython-compat shims, Playwright e2e suite, all in Feb 2026). No public hosted instance was found — you must run Docker Compose yourself.

---

## 1. What the tool is

Three fixed workflows, each a separate solver and a separate REST endpoint (`api/api.py`):

- **SSM / SSSM** (site-scanning saturation mutagenesis) — `POST /v1/ssm`. Input: plasmid sequence + GOI + flanking primers + list of positions + one degenerate codon. Output: overlapping forward/reverse mutagenic **primer pairs** for parallel PCR + SOEing.
- **QCLM / MSDM** (multi site-directed mutagenesis, QuikChange Lightning Multi kit) — `POST /v1/qclm`. Input: GOI + list of `E52W`-style substitutions. Output: unidirectional mutagenic **primers** (degenerate where useful) + their mix ratios.
- **PAS** (PCR-based accurate synthesis, de novo gene synthesis) — `POST /v1/pas`. Input: GOI as DNA *or protein* + per-site mutation lists with per-site **frequencies** + avoided motifs + host organism. Output: set of overlapping **oligos** with mix ratios.

Paper, abstract: *"we developed Mutation Maker, an open source mutagenic oligo design software for large-scale protein engineering experiments."*

Critically: **the outputs are oligos/primers, never the enumerated library-member sequences.** There is no object anywhere in the codebase representing "the set of variant sequences in this library".

---

## 2. Capability-by-capability evidence

### BLOCK A — library specification

**library_as_object — PARTIAL.**
One job covers many sites at once, so the user does not loop; but the "library" is only ever a *design job*, not a collection of sequences. `QCLMInput` (`backend/mutation_maker/qclm_types.py`):
```python
class QCLMInput(JsonObject):
    sequences = ObjectProperty(QCLMSequences, required=True)
    config = ObjectProperty(QCLMConfig, required=True)
    mutations = ListProperty(str, required=True)
```
and the output is `QCLMOutput.results: List[QCLMMutationOutput]`, where each record is `{mutations: [str], result_found: bool, primers: [PrimerOutput]}`. No variant-sequence list, no library container type.

**dag_chaining — NO.**
Three monolithic solvers (`ssm_solve`, `qclm_solve`, `pas_solve` in `backend/tasks.py`), reached by three unrelated endpoints (`@hug.post('/ssm')`, `/qclm`, `/pas`). Nothing composes, nests or chains design steps; you cannot feed one workflow's output into another programmatically. `AGENTS.md`: *"Mutation Maker is a mutagenic primer design application with three main workflows: SSM …, QCLM …, PAS …"*

**lazy_evaluation — NO.**
Eager throughout. Paper: *"the brute-force algorithm precomputes all possible forward and reverse primers around each mutation site"*; the Celery task returns a fully materialized JSON result stored in Redis (`api/api.py` `/v1/result/{task_id}`).

**mixed_mutagenesis_one_pool — PARTIAL.**
Within *one* QCLM or PAS job you can mix heterogeneous per-site specifications, including deliberate wild-type/parental retention and per-site frequencies:
- QCLM shipped fixture (`backend/tests/test_support.py`, `generate_qclm_input(ind=3)`): `["P269P","P269D","A271A","A271R","V272V","V272N","Q376Q","Q376R", …]` — the `X→X` entries are parental codons, so WT and mutant coexist at each site.
- PAS: `PASMutationSite(position=10, mutations=[PASMutation(mutation="I,V,F,S,R,T", frequency=0.5)])` alongside `position=12, mutation="Y", frequency=0.1` — different substitution sets *and different frequencies* per site in one design. Paper: *"the algorithm computes the oligo ratios in the reaction mixture that are necessary to achieve the user-deﬁned mutation frequency in the ﬁnal library."*

But: (a) SSM applies a **single** degenerate codon to the whole job — `SSMInput.degenerate_codon = StringProperty(required=True, default="NNS")`, one field, not per-site; (b) the three workflows cannot be combined in one specification; (c) there is no notion of sampled/random variants or WT replicate spike-ins as first-class library components.

**combinatorial_motif_place — NO.**
"Motif" in Mutation Maker means *sequence to avoid*, not element to place. `PASConfig.avoided_motifs = ListProperty(str)`; `backend/mutation_maker/motifs.py`; `frontend/src/shared/motifs.json` is a list of restriction enzyme names (`AarI, AatII, Acc65I, AccI, …`). Paper: *"exclusion of hairpins, homodimers and certain motifs"*. There is no insertion, no orientation/permutation/position enumeration of elements.

**barcode_generation — NO.**
GitHub code search for `barcode` in the repo returns 1 hit, in `frontend/feature-viewer/css/bootstrap.css` (irrelevant CSS). No barcode/UMI/edit-distance machinery anywhere.

**per_sequence_provenance — PARTIAL.**
Each designed oligo carries structured metadata, but it is physical-property metadata plus mutation identity, not a build history. E.g. `PrimerOutput` = `{sequence, start, length, temperature, gc_content, degenerate_codons, overlap_with_following}`; `SSMMutationOutput` = `{mutation, non_optimality, parameters_in_range, result_found, forward_primer, reverse_primer, overlap}`; `PASOligoOutput` = `{sequence, mix_ratio, mutations (indices), reds, blues}` where `blues` marks wild-type codon positions and `reds` mutated ones (`backend/mutation_maker/pas_output.py`). Excel export columns (`frontend/src/shared/components/SaveFile/index.tsx`): QCLM `['Name','Primer Sequence','Scale','Purification','Mutation Syntax','Overlap']`, PAS `['Name','Fragment Sequence','Mutation Syntax','Mix ratio','Target and MT%']`. No provenance for library members (they are never emitted) and no record of which operations produced a given sequence.

**automatic_naming — YES.**
Names are generated automatically, and are plate/well aware. `SaveFile/index.tsx`:
```ts
`${oligoPrefix}-${Math.floor(index / PLATE_SIZE) + 1}F-${getPosition(index % PLATE_SIZE)}`   // SSM forward
`${oligoPrefix}-${siteIndex + 1}${site.length > 1 ? String.fromCharCode(CHAR_LOWER_A + index) : ''}` // QCLM
```
with `PLATE_SIZE = 96` and `getPosition` producing `A01`-style wells. A separate "Mutation Syntax" column carries `m.mutations.map(prop('identifier')).join(',')` plus the degenerate codons.

**design_visualization — YES.**
Paper, Methods: *"A JavaScript (JS) User Interface (UI) where users can submit oligo design tasks and visualize the results using a customized version of the neXtProt viewer"*, and *"The following features were added: GC content graph, primer directionality, mutation alignment, and integration with interactive results table and with React user interface."* Results export to PDF or a temporary shareable link. Repo confirms: `frontend/feature-viewer` is a fork of `calipho-sib/feature-viewer` with an `arrow` object type "for primer visualisation". It visualizes the *primer layout on the gene*, not a design graph.

### BLOCK B — assay coverage

**assay_dms — YES (with caveat).**
This is the tool's core purpose: coding-sequence mutagenesis libraries for protein engineering. Paper: *"SSSM can generate a diverse gene library … a set of amino acids in a given protein sequence are subjected to random mutagenesis (at a single site out of deﬁned positions per clone)"*, and SSSM *"enable[s] parallel and deep mutational scans"*. Caveat for the table: it designs the **reagents** (primers/oligos) for such a library, and requires the user to hand it the position list; it never emits the variant sequences. Positions are amino-acid indices in the supplied gene (`AminoMutation`, `mutation_string[0]` = old AA, `[-1]` = new AA or `X`).

**assay_mpra — NO.**
Every entry point is codon/amino-acid indexed on a coding gene (`parse_codon_mutation` computes `zero_based_base_position = (one_based_codon_position - 1) * 3`). There is no promoter/enhancer/regulatory-element concept, no TF motif insertion, no scrambles, no MPRA vocabulary anywhere in paper or repo.

**assay_insilico — NO.**
No model-in-the-loop functionality; the paper's stated purpose is wet-lab oligo design and its validation is PCR/Sanger. Nothing about scoring sequences with predictive models.

### BLOCK C — genomics integration

All six are **NO**. Verified both in the paper and by GitHub code search over `ra100/Mutation_Maker`:

- **genome_coordinates — NO.** Inputs are raw sequence strings: `SSMSequences{forward_primer, reverse_primer, gene_of_interest, five_end_flanking_sequence, three_end_flanking_sequence, plasmid}`, `Plasmid{gene_start_in_plasmid, gene_end_in_plasmid, plasmid_sequence}`. Positions are 1-based *codon* indices in the GOI. Code search `chromosome` → 0 hits, `genome` → 0 hits.
- **transcript_models — NO.** Code search `GTF` → 0 hits, `gff` → 0 hits. No annotation ingestion of any kind.
- **exon_intron_split_codons — NO.** Code search `exon` → 0 hits, `intron` → 0 hits. `AminoMutation` hard-codes `self.length = 3` and contiguous codons: `zero_based_base_position = (one_based_codon_position - 1) * 3`.
- **hgvs_input — NO.** Code search `hgvs` → 0 hits. Mutation syntax is a bespoke one-letter format `"E52W"` (`parse_codon_mutation`: `mutation_string[0]`, `int(mutation_string[1:-1])`, `mutation_string[-1]`), which is *not* HGVS (`p.Glu52Trp` / `p.E52W` prefixes are not parsed and would raise `ValueError`).
- **vcf_vep_output — NO.** Code search `vcf` → 0 hits. Outputs are JSON (REST), XLSX (client-side via `exceljs`/`xlsx`), and PDF (browser print).
- **consequence_annotation — NO.** Nothing labels stop-gained / synonymous / in-frame. The nearest things are (i) stop codons being *excluded by default* from SSSM degeneracy solutions (paper: *"only stop codons are excluded by default"*) and (ii) a boolean `wild_type` flag on PAS mutations used only to color the viewer (`reds`/`blues` in `pas_output.py`). Neither is molecular-consequence annotation.

### BLOCK D — physical construction

**primer_design — YES.** This is the entire tool. Two SSSM algorithms (brute force with Primer3, and a "fast-approximation" primer-growth algorithm, the default), an MSDM partition/extension algorithm, and a PAS constraint-satisfaction + backtracking gene-splitter. Constraints: primer length, 5′/3′ end length, 3′-end Tm, overlap size and Tm, GC content, GC clamp, hairpin/homodimer/heterodimer Tm, Tm consistency across the whole set. Tm via `primer3-py` with configurable model (`TemperatureConfig`: SantaLucia_1997/2004, Breselauer_1986, Sugimoto_1996; Owczarzy/SantaLucia salt corrections; Na/K/Tris/Mg/dNTP/primer concentrations). Weighted-least-squares scoring with user-adjustable weights and explicit reporting of violated constraints. Experimentally validated: SSSM 79.9 % mean success across 7 libraries (Table 1), MSDM 72 % with 2.5 mutations/clone (Table 2), PAS full-length 1.2 kb gene, 32/32 clones correct.

**codon_optimization — YES.** `backend/mutation_maker/reverse_translation.py`, `codon_usage_table.py`, `usage_table.py`, `translation_scoring.py`. Paper: *"Mutation Maker's reverse translation algorithm generates a DNA sequence from a set of randomly selected codons of the host organism (excluding rare codons) … the algorithm computes the codon frequencies product (CFP) which is proportional to codon adaptation index (CAI)"*, benchmarked against GeneGenie's CAI optimizer (Fig. 2D). Ships codon-usage data for **35,799 species** and **33 taxa-specific genetic code tables** (`backend/mutation_maker/codon_usage_data/*.spsum`, `species.table`). `PASConfig.organism = "e-coli"`, `codon_usage_frequency_threshold = 0.1` (rare-codon cutoff).

**synthesis_constraints — YES.** Each designed oligo/fragment is checked against GC-content bounds (`min_gc_content=40`, `max_gc_content=60`), length bounds (`min/opt/max_oligo_size = 40/56/90`), overlap length and Tm, hairpin and homodimer formation temperature (penalized when within `safe_temp_difference` of reaction Tm), and a user-supplied excluded-motif list (restriction sites). Paper: *"a constraint-satisfaction algorithm with a backtracking component that splits a given gene sequence into an even number of overlapping DNA fragments that conform to length, Tm, mutation location, GC content, hairpins, homodimers and motif exclusion constraints"*. Note these are constraints on the **oligos**, not on final assembled library members.

### BLOCK E — engineering

**interface — YES.** Self-hosted **web GUI** (React + TypeScript + Ant Design, `localhost:3000`) is the primary interface; a **REST API** (`hug`/`falcon` + Gunicorn, `localhost:8000`) exposes `POST /v1/ssm`, `/v1/qclm`, `/v1/pas`, `GET /v1/get_species`, plus `/v1/check|result|cancel|forget/{task_id}` and `/v1/export_{ssm,qclm,pas}/{task_id}.xlsx`; the solvers are an importable Python package `mutation_maker` (`from mutation_maker.qclm import qclm_solve`). **No CLI and no PyPI/conda distribution** — deployment is Docker Compose or manual (Redis + Celery + Gunicorn + nginx + optional AWS Lambda for Primer3 parallelism). Version string in `backend/tasks.py`: `print("Mutation Maker version: 1.0.0")`.

**license — YES.** GPL-3.0 (GitHub API `license.spdx_id = "GPL-3.0"`; every source file carries the MSD GPLv3 header: *"Mutation Maker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version."*). The article itself is CC-BY.

**maintained — YES (but not at the published URL).** `Merck/Mutation_Maker` → 404. `ra100/Mutation_Maker`: last commit **2026-02-14** ("⚡ Reduce backend test count for faster CI execution"), CI + Primer3 workflows both green on 2026-02-14, 0 open issues, not archived. Recent work is genuine modernization (React 18 / antd 4 migration plans in `docs/plans/2026-02-*`, Python 3.11 in the Ubuntu deployment instructions, Biopython `GC`→`gc_fraction` compat shim, Playwright e2e suite). **No tagged releases and no GitHub releases** (`/tags` and `/releases` both empty).

---

## 3. The tool's own documented examples / vignettes

There is no tutorial or vignette site. What exists is (a) the paper's experimental validations and (b) shipped test fixtures.

### From the paper (sequences NOT disclosed — "Gene X" / "Gene Y" are proprietary Merck genes)
1. **SSSM libraries, Table 1** — 7 libraries on Gene X (~700 bp) and Gene Y (~1.2 kbp), 31–96 target sites each, vectors pAVE011-Kan / pET30a(+), hosts *E. coli* W3110 / BL21(DE3); 24–48 clones Sanger-sequenced; mean success 79.9 %.
2. **MSDM/QCLM libraries, Table 2** — 13 designs on Gene X, 13–18 target mutations each, primer:template ratios Low/Medium/High, pET30a(+) / BL21(DE3); mean success 72 %, 2.5 mutations per clone.
3. **PAS de novo synthesis, Figure 3** — ~1.2 kb gene split into 32 overlapping fragments (default settings), of which 6 fragments carried **16 mutations at 8 sites**, 50 oligos total; assembled in two PCR rounds via 3-way (A1–A3) or 4-way (B1–B4) segment splits; 32/32 clones full-length.
4. **In-silico benchmarks** — 50 simulated genes × 10 mutations for brute-force vs fast-approximation (Fig. 2A/B); 10 genes with 2–10 sites and ≤5 mutations/site for MSDM coverage weights (Fig. 2C); 5 random 250-aa proteins reverse-translated, MM-CFP vs GeneGenie-CAI (Fig. 2D); random amino-acid groups, MM degeneracy vs CodonGenie (Fig. 2E).

### From the repository (`backend/tests/test_support.py`) — these ARE runnable
5. `generate_SSM_input(mut_ind)` — 10 preset SSM mutation sets over a bundled plasmid (`sample_ssm_sequence()`), e.g. `["L50X","E60X","W70X","R80X","E100X","L120X","G160X","G200X","I220X","L240X"]`; `X` = saturate with the job's degenerate codon.
6. `generate_random_SSM_input(...)` — randomized site generator used for the scaling benchmarks (`backend/ssm_algorithm_comparison.py`).
7. `generate_qclm_input(ind=0..9)` — 10 preset MSDM cases, e.g. `ind=2`: `E52W/L/F/A, V73L/G/F/I/S, A165P/E/R/K, R169I/K/L/Y` (4 sites × 4–5 AAs); `ind=3`: 21 entries including parental codons (`P269P, A271A, V272V, Q376Q, M379M, K383K, T408T, H411H, L412L`); `ind=4`: a 13-site cysteine scan (`G310C, H348C, I6C, L9C, S11C, G15C, A16C, A20C, E21C, V26C, G27C, L31C, G32C`) with degeneracy off.
8. `generate_pas_input(ind=1..5)` — PAS cases with per-site frequencies, e.g. `ind=3`: `pos 1 → "R,L" @0.25`, `pos 2 → "R" @0.9`, `pos 10 → "I,V,F,S,R,T" @0.5`, `pos 12 → "Y" @0.1`, … (11 sites). `ind=5` uses a **protein** input sequence (triggers reverse translation).
9. `e2e-tests/tests/app.spec.ts` (32 Playwright UI tests) and `api.spec.ts` (9 API tests) exercising all three workflow forms and endpoints.

**Reproducibility note for PoolParty:** items 1–3 cannot be reproduced literally because the gene sequences are withheld. Items 5–8 are concrete and self-contained — the natural PoolParty comparison target is #7 (`ind=2`: 4 sites × 4–5 defined substitutions + WT) and #8 (per-site substitution sets with per-site frequencies), which PoolParty can express as a library specification even though PoolParty does not design the primers.

---

## 4. Notable capabilities NOT covered by the current row list

1. **Degenerate-codon set-cover solvers.** Two distinct heuristics: (a) SSSM — given amino-acid *inclusion* and *exclusion* lists, find a minimal degenerate codon/alphabet; (b) MSDM/PAS — degeneracy strictly restricted to the specified substitutions *plus the parental amino acid*, ranked by codon-frequency product. Benchmarked against **CodonGenie** (Fig. 2E). This is a whole capability axis (degenerate-codon compression) that neither PoolParty nor most other surveyed tools have.
2. **Oligo mixing-ratio calculation to hit a target variant frequency.** PAS computes the molar ratio of each oligo in the reaction so the *final library* has the user-specified per-mutation frequency (`PASOligo.ratio`, exported as "Mix ratio"). This is library-composition control expressed physically rather than as an enumerated pool.
3. **Species-specific codon usage at scale** — 35,799 species (CUTG/EBI `.spsum` tables) × 33 NCBI genetic code tables, with a rare-codon frequency threshold.
4. **Configurable thermodynamics** — nearest-neighbour model, salt correction, and full buffer composition (Na⁺/K⁺/Tris/Mg²⁺/dNTP/primer conc.) are user-settable; hairpin/homo-/heterodimer Tm penalties relative to the reaction temperature.
5. **96-well plate layout assignment** in the Excel export (well positions `A01`…, plate index, per-plate sheets, "Scale"/"Purification" columns ready for a synthesis vendor order form).
6. **Constraint-violation reporting with user-tunable weights** — soft vs hard constraints; the UI reports *which* constraints a returned primer violates and its penalty score (`non_optimality`, `parameters_in_range`).
7. **Restriction-site avoidance from a named-enzyme list** (`frontend/src/shared/motifs.json`, NEB enzyme names) — a synthesis-constraint feature keyed by enzyme name rather than raw sequence.
8. **Asynchronous job infrastructure** — Celery + Redis queue, cancel/forget/poll endpoints, optional AWS Lambda fan-out for Primer3, shareable temporary result links, PDF export.

---

## 5. Stated limitations (verbatim from the Discussion)

- *"the SSSM brute-force approach is only exhaustive within fairly large Tm windows that slide along the user-speciﬁed Tm range. Thus, in essence, not all possible primers across the speciﬁed Tm range are interrogated."*
- *"Mutation Maker by design does not optimize for the presence of GC clamps throughout, since this constraint often eliminated better scoring solutions."*
- *"The scoring functions … have not been optimized using exhaustive laboratory testing and therefore are subject to improvement over time."*
- *"Mutation Maker's gene split algorithm can be also further optimized, since only simpliﬁed scores were used to evaluate fragments during backtracking as a trade-off for reduced computational time."*
- *"The gene synthesis algorithm may not guarantee an optimal solution for complex cases when it hits the hardcoded timeout for each Tm iteration."*
- *"the algorithms for reverse translation and for degeneracy are solving these problems using heuristics and randomization and therefore may not output the optimal solution."*
- *"the gene synthesis algorithm only checks for secondary structures at the fragment level and does not scan the entire reverse-translated sequence for resulting mRNA secondary structures that may affect translation efficiency."*
- MSDM coverage is not guaranteed: *"The MSDM algorithm may report an incomplete primer set that covers only some of the speciﬁed mutation sites"* (coverage enters the score with weight `mutation_coverage_weight = 160`).

Undocumented issue observed in the live code: `api/api.py` schedules Celery tasks `tasks.export_ssm` / `tasks.export_qclm` / `tasks.export_pas` for its `/v1/export_*.xlsx` endpoints, but `backend/tasks.py` defines no `export*` task (0 occurrences of "export"). Excel export in practice is done client-side in the React app (`frontend/src/shared/components/SaveFile/index.tsx`), so the server-side export endpoints appear vestigial/broken.

---

## 6. Corrections to the prior analysis

The prior note (`prior_analyses/Hiraga2021_MutationMaker_analysis.md`) is broadly correct — "designs primers for the variants the user has already specified" is accurate. Three refinements:

1. It says the repo is `https://github.com/Merck/Mutation_Maker`; **that URL is now 404**. Any citation should point at the surviving copy or be worded carefully.
2. "Combinatorial logic: Limited (MSDM combinations)" understates PAS: PAS does combine mutations across sites within a fragment (`cartesian_dictionary` / `itertools.product` in `generate_oligos.py`) and computes oligo mix ratios to hit per-site target frequencies. The combinatorics are real, but they live at the oligo level and are never enumerated as library-member sequences.
3. "Design cards: No" and "MPRA support: No" are confirmed; add that the tool has genuinely strong per-oligo metadata and automatic plate-aware naming, so a "no metadata / no naming" framing would be wrong and would be easy for the authors to rebut.

---

## 7. Confidence

- **High confidence:** all BLOCK C rows (no genome/transcript/HGVS/VCF/consequence support) — checked by code search and by the input type definitions; `primer_design`, `codon_optimization`, `license`, `maintained`, `interface`, `barcode_generation`, `combinatorial_motif_place`, `dag_chaining`.
- **Medium confidence / judgement calls:**
  - `assay_dms` = yes: defensible (SSSM *is* deep-mutational-scan library construction) but it designs reagents, not variant sequences. Anyone wanting a stricter reading could call it "partial".
  - `library_as_object` = partial: depends on whether "library" means the design job (yes) or the sequence collection (no).
  - `mixed_mutagenesis_one_pool` = partial: PAS/QCLM genuinely mix WT + heterogeneous per-site substitution sets with frequencies in one job; they cannot mix *mutagenesis types* (saturation + defined + sampled) in one specification.
  - `per_sequence_provenance` = partial: rich per-oligo metadata, but no compositional build history and no per-library-member records.
  - `lazy_evaluation` = no: inferred from eager precomputation and materialized JSON results; there is no generator/iterator API, but I did not exhaustively audit every solver for internal streaming.
