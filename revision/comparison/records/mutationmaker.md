# Mutation Maker — FINAL capability record

**Tool:** Mutation Maker
**Slug:** `mutation_maker`
**Citation key:** `Hiraga2021yg`
**Tier:** 3
**Paper:** Hiraga K, Mejzlik P, Marcisin M, Vostrosablin N, Gromek A, Arnold J, Wiewiora S, Svarba R, Prihoda D, Clarova K, Klempir O, Navratil J, Tupa O, Vazquez-Otero A, Walas MW, Holy L, Spale M, Kotowski J, Dzamba D, Temesi G, Russell JH, Marshall NM, Murphy GS, Bitton DA. "Mutation Maker, An Open Source Oligo Design Platform for Protein Engineering." *ACS Synth Biol* 2021;10(2):357–370. doi:10.1021/acssynbio.0c00542 (article CC-BY; PMID 33433999)
**Local PDF:** `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/../../../lit_review/analyzed/Hiraga2021_MutationMaker.pdf`
**Status:** final — supersedes `extractions/mutationmaker.md ⟨deleted at cleanup — in commit 35d65d8⟩`; incorporates all corrections from `reviews/mutationmaker.md ⟨deleted at cleanup — in commit 35d65d8⟩`.
**Compiled:** 2026-08-10

---

## 0. Sources consulted

| Kind | Reference |
|---|---|
| pdf | `../../../lit_review/analyzed/Hiraga2021_MutationMaker.pdf` (14 pp., extracted with PyMuPDF) |
| repo (dead) | https://github.com/Merck/Mutation_Maker — **HTTP 404** as of 2026-08-10 (`gh api repos/Merck/Mutation_Maker` → `Not Found`). This is the URL printed in the paper. |
| repo (live, author-maintained) | https://github.com/ra100/Mutation_Maker — network root, retains original Merck history. Default-branch commit `396c1c0ede7529f3dedf65e821e8c1f20c9a7043` names Rastislav Švarba as committer; **Rastislav Svarba is the 8th listed author of Hiraga et al. 2021** and MSD is the copyright holder. GPL-3.0, last commit **2026-02-14**, CI + Primer3 workflows both green 2026-02-14. |
| docs | Repo `README.md`, `backend/README.md`, `AGENTS.md`, `e2e-tests/README.md`; the repository has no GitHub Pages (`has_pages=false`). Journal Supporting Information is only *"Figures S1−S4 illustrate Mutation Maker algorithms (PDF)"* — no user guide or tutorial. |
| pypi | No PyPI package (`mutation-maker`, `mutation_maker`, `mutationmaker` JSON endpoints all 404 on 2026-08-14). Stronger: **no `setup.py`, `pyproject.toml` or `setup.cfg` anywhere in the 367-path repo tree** (verified via `gh api .../git/trees/master?recursive=1`), so it is not pip-installable even from a git URL. |

**Source scope and labels.** Repository claims refer to the official `master` branch at `396c1c0ede7529f3dedf65e821e8c1f20c9a7043`, verified as its head on 2026-08-14; “current” means that snapshot. `SSM`/`SSSM` and `QCLM`/`MSDM` are the code/paper labels for the same two workflows, respectively.

**Availability headline.** The repository URL printed in the paper (`github.com/Merck/Mutation_Maker`) now 404s. The code is maintained by **co-author Rastislav Svarba** at `github.com/ra100/Mutation_Maker`. **Do not describe this as a "community fork" or "community-maintained copy"** — it is an author-maintained repository. Maintenance is continuous rather than a single burst: commits by month 2020-07 (9), 2020-11 (3), 2021-01 (4), 2021-02 (10), 2022-08/09/10 (20), 2025-07 (7), 2026-02 (32). Not archived, 0 open issues, but **0 tags and 0 releases**. The repository documents self-deployment with Docker Compose (Redis + Celery + Gunicorn + nginx, optional AWS Lambda Primer3 fan-out).

---

## 1. What the tool is

Three fixed workflows, each a separate solver and a separate REST endpoint (`api/api.py`):

- **SSM / SSSM** (site-scanning saturation mutagenesis) — `POST /v1/ssm`. Input: plasmid + GOI + flanking primers + list of positions + a degenerate-codon field. Output: overlapping forward/reverse mutagenic **primer pairs** for parallel PCR + SOEing.
- **QCLM / MSDM** (multi site-directed mutagenesis, QuikChange Lightning Multi) — `POST /v1/qclm`. Input: GOI + list of `E52W`-style substitutions. Output: unidirectional mutagenic **primers** (degenerate where useful) + mix ratios; nearby sites can be grouped so one returned primer carries substitutions at multiple sites.
- **PAS** (PCR-based accurate synthesis / de novo gene synthesis) — `POST /v1/pas`. Input: GOI as DNA *or protein*, per-site mutation lists with per-site **frequencies**, avoided motifs, host organism. Output: overlapping **fragments**, each carrying an enumerated set of **variant oligo sequences** with normalised mixing ratios, plus the reverse-translated full-length gene echoed back in `PASOutput.input_data`.

Paper, abstract: *"we developed Mutation Maker, an open source mutagenic oligo design software for large-scale protein engineering experiments."*

**Correct framing of the core gap vs. PoolParty** (this replaces the earlier, falsifiable claim that the tool "never emits enumerated library members"): Mutation Maker **enumerates variant sequences only at the oligo/fragment level, and only in PAS**. It never emits full-length assembled library members, and SSM and QCLM emit primers only; the distinction is PAS fragment-level enumeration versus full-length member enumeration.

---

## 2. Capability-by-capability record

### BLOCK A — library specification

**`library_as_object` — PARTIAL.**
One design job covers many sites at once, so the user does not loop. There is a job-level input object but no full-library container. `QCLMInput` (`backend/mutation_maker/qclm_types.py`) is `{sequences, config, mutations: List[str]}` and `QCLMOutput.results` is `List[QCLMMutationOutput]` = `{mutations, result_found, primers}` — no variant-sequence collection anywhere in SSM or QCLM.
*Correction from review:* enumeration **does** exist in PAS. `backend/mutation_maker/generate_oligos.py::generate_oligos_from_combinations()` iterates `cartesian_dictionary(*mutations_on_site_with_prob)` (= `itertools.product` over the per-site codon dictionaries, folding per-site frequencies by multiplication), writes each codon combination into the fragment DNA via `Codons.replace_codon`, and appends one `PASOligo(sequence=dna, ratio=concentration)` per combination; ratios are then normalised (`oligo.ratio /= total_ratio`, asserted to sum to 1). Those *are* enumerated variant sequences with per-variant frequency — scoped to a fragment.
Library specification can also be **declarative and tabular**: PAS mutations upload as XLSX with columns `Site` / `MT` / `MT%` (`frontend/src/scenes/PAS/components/PASForm/components/InputMutations/FileUploadMutations.tsx`, `const header = ['site','mt','mt%']`, `FILE_EXTENSIONS = '.xlsx'`), and sequences upload as `.fa/.fasta/.gb`.
**Partial, not yes, because:** enumeration exists at oligo/fragment level in PAS only, never as full-length assembled library members, and not at all in SSM/QCLM; there is no first-class object representing "the library" that can be inspected, filtered, or composed.
*Source:* `backend/mutation_maker/generate_oligos.py`, `qclm_types.py`, `pas_types.py`, `FileUploadMutations.tsx`.

**`dag_chaining` — NO.**
`api/api.py` defines only `@hug.post('/ssm')`, `/qclm`, `/pas` (+ `/species_table`), each dispatching to a separate Celery task (`tasks.ssm`, `tasks.qclm`, `tasks.pas`, `tasks.species_table` in `backend/tasks.py`). No composition, nesting, or chaining primitive; there is no way to feed one workflow's output into another. `AGENTS.md`: *"three main workflows: SSM …, QCLM …, PAS …"*.
*Source:* `api/api.py`, `backend/tasks.py`, `AGENTS.md`.

**`lazy_evaluation` — NO.**
Fully eager. Solvers return materialised JSON (`PASOutput(...).to_json()`, `SSMOutput`, `QCLMOutput`) stored in the Celery/Redis result backend and retrieved by `GET /v1/result/{task_id}`. No generator/iterator/streaming API at any user-facing level; internal `itertools.product` is eagerly listified (`self.concrete_mutations = list(itertools.product(...))` in `mutation.py`). Paper: *"the brute-force algorithm precomputes all possible forward and reverse primers around each mutation site."*
*Source:* `api/api.py`, `backend/mutation_maker/{pas,ssm,qclm,mutation}.py`, paper p. 359.

**`mixed_mutagenesis_one_pool` — PARTIAL.**
Within one QCLM or PAS job you can mix heterogeneous per-site specifications, and **wild-type retention is automatic, not a user trick**:
- QCLM: `MutationSite.__init__` (`backend/mutation_maker/mutation.py:147`) does `self.new_aminos = frozenset(old_aminos.union({m.new_amino for m in mutations}))` — the parental amino acid is unioned into **every** site's target set unconditionally.
- PAS: in `OligoGenerator.generate_solution`, `if sum_of_probabilities != 1:` the wild-type codon at that site is added with `wild_type_prob = 1 - sum_of_probabilities`, so WT is a first-class member of the enumerated pool at any site whose specified frequencies do not sum to 1.
- PAS accepts different substitution sets **and different frequencies** per site in one design (fixture: `position=10, "I,V,F,S,R,T" @0.5` alongside `position=12, "Y" @0.1`). Paper: *"the algorithm computes the oligo ratios in the reaction mixture that are necessary to achieve the user-defined mutation frequency in the final library."*
- Mutations may be given as **raw DNA codons** rather than amino acids (`PASInput.is_mutations_as_codons = BooleanProperty(required=True)`; `PASMutationSite` docstring "Each mutation is a amino acid IUPAC code, or a DNA codon"; parse path `(mutation.position, str(mut)[:3], float(str(mut)[3:]))`), giving nucleotide-level control of the exact codon placed.

**Downgraded to partial because:** (a) SSM applies one **job-wide** `degenerate_codon` field (`ssm_types.py:215`, `StringProperty(required=True, default="NNS")`) — the UI's "Amino acids" tab writes back into that single field; (b) the three workflows cannot be combined in one specification; (c) there is no sampled/random-variant or WT-replicate spike-in concept.
*Source:* `backend/mutation_maker/mutation.py`, `generate_oligos.py`, `pas_types.py`, `ssm_types.py`, `frontend/src/scenes/SSM/components/SSMForm.tsx`, paper p. 362.

**`combinatorial_motif_place` — NO.**
"Motif" in Mutation Maker means *sequence to avoid*: `PASConfig.avoided_motifs = ListProperty(str)`, `backend/mutation_maker/motifs.py`, `frontend/src/shared/motifs.json` is a restriction-enzyme name list (`AarI, AatII, Acc65I, …`). Paper: *"exclusion of hairpins, homodimers and certain motifs."* No element placement, orientation, permutation, or position enumeration.
*Stronger point added on review:* **Mutation Maker supports no insertions or deletions at all.** `AminoMutation.__init__` hardcodes `self.length = 3` (`mutation.py:67`; also `ConcreteTripletMutation`, `MutationSite`), and every mutation is a codon-for-codon substitution. GitHub code search: `insertion` → 1 irrelevant hit, `deletion` → 0.
*Source:* `backend/mutation_maker/{pas_types,motifs,mutation}.py`, `frontend/src/shared/motifs.json`, paper p. 360.

**`barcode_generation` — NO.**
GitHub code search `barcode` over `ra100/Mutation_Maker` → 1 hit, `frontend/feature-viewer/css/bootstrap.css` (irrelevant CSS). No UMI, edit-distance, or tag-set machinery in code or paper.
*Source:* GitHub code search; paper (full text).

**`per_sequence_provenance` — PARTIAL (evidence corrected; a lenient reading would support "yes").**
*The earlier claim "no records exist for library-member sequences (they are never emitted)" was false and has been removed.* Every sequence PAS emits carries a per-sequence record: `PASOligoOutput{sequence, mix_ratio, mutations, reds, blues}` where `mutations` is an untyped `ListProperty` of integer indices into the enclosing `PASResult.mutations`, which is the `ListProperty(PASMutationFormatted)` = `{position, mutated_amino, wild_type_amino, wild_type_codon, mutated_codon, frequency, wild_type}`, and `reds`/`blues` are the mutated vs. wild-type nucleotide start offsets within the oligo (`3*(position-1) + goi_offset - fragment_start`, used for display colouring); the enclosing `PASResult` adds `{fragment, start, end, length, overlap, overlap_Tm, overlap_GC, overlap_length}`. Primers likewise carry `PrimerOutput{sequence, start, length, temperature, gc_content, degenerate_codons, overlap_with_following}` and `SSMMutationOutput{mutation, non_optimality, parameters_in_range, result_found, forward_primer, reverse_primer, overlap}`. Excel exports carry it through (`frontend/src/shared/components/SaveFile/index.tsx` headers: SSM `['Well Position','Name','Sequence','Notes']`, QCLM `['Name','Primer Sequence','Scale','Purification','Mutation Syntax','Overlap']`, PAS `['Name','Fragment Sequence','Mutation Syntax','Mix ratio','Target and MT%']`).
**Kept at partial** on two narrow grounds only: it is a *state* record (what mutation is present, at what frequency), not a *compositional build history* (no record of which design operations produced the sequence), and no records exist for full-length assembled members because those are never emitted. **Flagged:** a lenient reading of this row would score Mutation Maker "yes"; if the row is defined as "each emitted sequence carries a machine-readable record of what was applied to it", Mutation Maker meets it.
*Source:* `backend/mutation_maker/pas_types.py:158-196`, `pas_output.py`, `ssm_types.py`, `frontend/src/shared/components/SaveFile/index.tsx`.

**`automatic_naming` — YES.**
Names are generated automatically for all three workflows; only the SSM names are plate/well aware. `frontend/src/shared/components/SaveFile/index.tsx`: `PLATE_SIZE = 96`, `getPosition` yields `A01`-style wells, and
- SSM: `` `${oligoPrefix}-${Math.floor(index/PLATE_SIZE)+1}F-${getPosition(index%PLATE_SIZE)}` `` (and `…R-…` for reverse)
- QCLM: `` `${oligoPrefix}-${siteIndex+1}${site.length>1 ? String.fromCharCode(CHAR_LOWER_A+index) : ''}` ``
- PAS: `` `${oligoPrefix}-Fr${index+1}-${oIndex+1}` ``

with prefix from `input.config.oligo_prefix || 'prefix'`. A separate "Mutation Syntax" column carries `m.mutations.map(prop('identifier')).join(',')` plus degenerate codons.
*Source:* `frontend/src/shared/components/SaveFile/index.tsx:35,48,55-57,136,142,183,301,314`.

**`design_visualization` — YES.**
Paper, Methods: *"A JavaScript (JS) User Interface (UI) where users can submit oligo design tasks and visualize the results using a customized version of the neXtProt viewer"*; *"The following features were added: GC content graph, primer directionality, mutation alignment, and integration with interactive results table and with React user interface."* The paper also states that results *"can be exported to PDF or as a temporary sharable link"*, but neither exists in the current repository — there is no PDF-generation and no share-link code, only print-oriented CSS. Repo: `frontend/feature-viewer` (fork of `calipho-sib/feature-viewer`, adds an `arrow` object type for primers) plus `SSMFeatureViewer.tsx`, `QCLMFeatureViewer.tsx`, `PASFeatureViewer.tsx`. **Caveat kept:** it visualises the *primer/fragment layout on the gene*, not a design graph or library composition plot.
*Source:* paper Methods p. 366; `frontend/feature-viewer/`, `frontend/src/scenes/{SSM,QCLM,PAS}/components/*Result/components/*FeatureViewer.tsx`.

### BLOCK B — assay coverage

**`assay_dms` — YES (with caveat).**
Core purpose. Paper: SSSM *"can generate a diverse gene library … a set of amino acids in a given protein sequence are subjected to random mutagenesis (at a single site out of defined positions per clone)"*, and the primers *"enable parallel and deep mutational scans."* Validation is 7 real DMS-style libraries with 31–96 target sites each (Table 1), 79.9 % mean success. **Caveat kept:** it designs the *reagents* for such a library from a user-supplied position list; it does not enumerate the resulting full-length variants.
*Source:* paper pp. 358-359, Table 1; `backend/mutation_maker/ssm.py`.

**`assay_mpra` — NO.**
Every mutation entry point is codon-indexed: `parse_codon_mutation` computes `zero_based_base_position = (one_based_codon_position - 1) * 3`; `AminoMutation` validates against IUPAC one-letter amino acids; PAS filters mutations onto fragments via `goi_offset + mut[0]*3`. No promoter/enhancer/regulatory vocabulary, no TF-motif insertion, no scrambles or shuffles anywhere in paper or repo. *Hedge for honesty:* PAS with `is_dna_sequence=True` will synthesise arbitrary DNA, so non-coding synthesis is technically possible — but mutation positions remain codon-indexed and there is no regulatory-element concept, so the "no" stands.
*Source:* `backend/mutation_maker/mutation.py`, `generate_oligos.py`, `pas_types.py`.

**`assay_insilico` — NO.**
No predictive sequence-function model and no model-in-the-loop code anywhere in `backend/mutation_maker/` (`primer_scoring.py` and `translation_scoring.py` do exist, but they are thermodynamic and codon-usage penalty functions, not predictive models). The paper's in-silico work (Fig. 2) is *algorithm* benchmarking — runtime, penalty scores, MM-CFP vs. GeneGenie CAI, degeneracy vs. CodonGenie — not sequence-function prediction. Experimental validation is PCR + Sanger.
*Source:* paper Fig. 2 and Results; repo file listing.

### BLOCK C — genomics integration

All six **NO**, verified in the paper and by GitHub code search over `ra100/Mutation_Maker` (`hgvs` 0, `vcf` 0, `gtf` 0, `gff` 0, `exon` 0, `intron` 0, `chromosome` 0, `genome` 0).

**`genome_coordinates` — NO.**
Inputs are sequences plus intra-gene indices: `SSMSequences{forward_primer, reverse_primer, gene_of_interest, five_end_flanking_sequence, three_end_flanking_sequence, plasmid}`, `Plasmid{gene_start_in_plasmid, gene_end_in_plasmid, plasmid_sequence}`; positions are 1-based **codon** indices in the GOI. Code search `chromosome` → 0, `genome` → 0.
*Pre-empting the obvious rebuttal ("it reads GenBank"):* `frontend/src/shared/components/FileUploadInput/index.tsx` does accept `FILE_EXTENSIONS = '.fa, .fasta, .gb'` — but the "GenBank parser" is `text.split('ORIGIN')[1].split('').filter(letter => /[ATCGatcgFLIMVSYHQNKDEWRPflimvsyhqnkdewrp]/.test(letter)).join('').toUpperCase()`. It discards LOCUS, FEATURES, every coordinate and every annotation, keeping only raw sequence letters. No genomic coordinate ever enters the system.
*Source:* `backend/mutation_maker/ssm_types.py:35-38,96-102`; `frontend/src/shared/components/FileUploadInput/index.tsx:25,35-37`.

**`transcript_models` — NO.**
Code search `GTF` → 0, `gff` → 0. The only external data files are codon-usage `.spsum` tables, `species.table` / `speciesInfo.table`, and NCBI genetic-code tables. Same GenBank caveat as above: `.gb` files are accepted by the GUI but all annotation is stripped before use, so no transcript model is ever ingested.
*Source:* GitHub code search; `backend/mutation_maker/codon_usage_data/`; `FileUploadInput/index.tsx`.

**`exon_intron_split_codons` — NO.**
Code search `exon` → 0, `intron` → 0. `AminoMutation.__init__` sets `self.length = 3` unconditionally, `parse_codon_mutation` computes `zero_based_base_position = (one_based_codon_position - 1) * 3`, and `mutations_on_fragments` in `generate_oligos.py` assumes `goi_offset + mut[0]*3` contiguity. Codons are contiguous by construction with no escape hatch.
*Source:* `backend/mutation_maker/mutation.py:28,67`, `generate_oligos.py:60-65`.

**`hgvs_input` — NO.**
Code search `hgvs` → 0. Mutation syntax is a bespoke one-letter format `"E52W"`: `parse_codon_mutation` reads `mutation_string[0]`, `int(mutation_string[1:-1])`, `mutation_string[-1]` and raises `ValueError('Position must be positive number')` otherwise; the SSM form regex is `/^\s*([A-Z][0-9]*[A-Z])(\s+[A-Z][0-9]*[A-Z]\s*)*$/i`. `p.Glu52Trp` and `c.155A>G` both fail. The format resembles HGVS protein shorthand *minus the prefix*, so the distinction is deliberate, not an oversight.
*Source:* `backend/mutation_maker/mutation.py:28-42`; `frontend/src/scenes/SSM/components/SSMForm.tsx`.

**`vcf_vep_output` — NO.**
Code search `vcf` → 0. Outputs are JSON over REST, client-side XLSX (`exceljs`/`xlsx` in `SaveFile/index.tsx`), and browser-print PDF.
*Source:* `api/api.py`, `frontend/src/shared/components/SaveFile/index.tsx`.

**`consequence_annotation` — NO.**
Nothing labels stop-gained, synonymous, missense, or frameshift. The two near-misses are not consequence calls: (i) stop codons are *excluded by default* from degeneracy solutions (paper: *"only stop codons are excluded by default"*; `defaultAvoid = [...aminoToCodonMap.STOP]` in `frontend/src/scenes/SSM/components/amino.ts`), and (ii) `PASMutationFormatted.wild_type = BooleanProperty(required=True)` only drives the `reds`/`blues` colouring in `pas_output.py`.
*Source:* paper p. 360; `backend/mutation_maker/pas_types.py:165`, `pas_output.py`; `frontend/.../amino.ts`.

### BLOCK D — physical construction

**`primer_design` — YES.**
The entire point of the tool, and if anything under-sold. Two SSSM algorithms (a brute-force search, with Primer3 as an optional candidate generator via `SSMConfig.use_primer3`, and a default "fast-approximation" primer-growth algorithm; the UI makes the two mutually exclusive), an MSDM partition/extension algorithm, and a PAS constraint-satisfaction + backtracking gene splitter. Repo modules: `ssm.py`, `ssm_fast_approximation.py`, `qclm.py`, `qclm_solution.py`, `pas.py`, `pas_back_track.py`, `pas_solution.py`, `primer_scoring.py`, `primer_filtering.py`, `primer3_interoperability.py`, `temperature_calculator.py`, `site_split.py`, plus an AWS Lambda Primer3 fan-out (`lambda/`) with a full thermodynamic parameter set. Constraints: primer length, 5′/3′ end length, 3′-end Tm, overlap size/Tm, GC content, hairpin/homodimer/heterodimer Tm, Tm consistency across the set; weighted-least-squares scoring with user-adjustable weights and explicit reporting of violated constraints (`non_optimality`, `parameters_in_range`). Tm via Biopython `Bio.SeqUtils.MeltingTemp` with configurable model (`TemperatureConfig`: SantaLucia 1997/2004, Breslauer 1986, Sugimoto 1996; Owczarzy/SantaLucia salt corrections; Na/K/Tris/Mg/dNTP/primer concentrations); `primer3-py` supplies the hairpin/homodimer/heterodimer calculations and the separate fixed `NEB_like` Tm. (The paper, p. 366, attributes Tm computation to `primer3-py`; the current code does not.) Wet-lab validated: SSSM 79.9 % mean success across 7 libraries (Table 1), MSDM 72 % with 2.5 mutations/clone (Table 2), PAS 1.2 kb gene 32/32 clones full-length (Fig. 3).
*Source:* repo `backend/mutation_maker/*`, `lambda/`; paper Tables 1-2, Fig. 3.

**`codon_optimization` — YES.**
`reverse_translation.py`, `codon_usage_table.py`, `usage_table.py`, `translation_scoring.py`, `degenerate_codon.py`. Paper: *"Mutation Maker's reverse translation algorithm generates a DNA sequence from a set of randomly selected codons of the host organism (excluding rare codons) … the algorithm computes the codon frequencies product (CFP) which is proportional to codon adaptation index (CAI)"*, benchmarked against GeneGenie's CAI optimizer (Fig. 2D). Ships codon usage for **35,799 species** and **33 taxa-specific genetic code tables** (paper p. 366; `codon_usage_data/` has 9 `.spsum` division tables plus `species.table`, 35,793 lines). `PASConfig.organism = "e-coli"`, `codon_usage_frequency_threshold = 0.1` (rare-codon cutoff).
*Added on review:* the optimised gene is actually **delivered to the user**, not just used internally — `pas.py::find_solution` calls `sequences.translate_goi_sequences(translator)`, which does `self.gene_of_interest = translator(self.gene_of_interest)` **in place** on `workflow_input.sequences`, and `pas_solve` then returns `PASOutput(input_data=workflow_input, results=results).to_json()`.
This **yes rests on PAS**. In the current QCLM web path, the frontend sends `codon_usage`/`taxonomy_id` while `QCLMConfig` reads `organism` (default `"e-coli"`), and QCLM invokes its degeneracy helper without an organism, so the selected custom host is not honoured there.
*Source:* `backend/mutation_maker/pas.py:83-99,160-166`, `pas_types.py:79-86`, `reverse_translation.py`; `frontend/src/services/api.ts:145-175`, `qclm_types.py:48-79`, `qclm.py:196,267-275`; paper pp. 361, 366, Fig. 2D.

**`synthesis_constraints` — YES.**
Verified directly in `PASConfig`: `min/max/opt_oligo_size = 40/90/56`, `min/max/opt_overlap_tm = 50/65/56`, `min/max/opt_overlap_length = 15/25/21`, `min/max_gc_content = 40/60`, `avoided_motifs`, `codon_usage_frequency_threshold = 0.1`, `safe_temp_difference = 10`, `hairpin_homodimer_weight = 2`, plus a full `TemperatureConfig`. Paper: *"a constraint-satisfaction algorithm with a backtracking component that splits a given gene sequence into an even number of overlapping DNA fragments that conform to length, Tm, mutation location, GC content, hairpins, homodimers and motif exclusion constraints."* **Caveat kept:** these constrain the **oligos/fragments**, not final assembled library members.
Current-code caveats do not remove the active length/Tm/GC/structure constraints that sustain the value: the mutated-oligo motif check searches the original fragment rather than each candidate, the frontend's `compute_hairpin_homodimer` field is not a `PASConfig` switch, and 40/56/90 nt are backend defaults (the web form uses 40/40/60).
*Source:* `backend/mutation_maker/pas_types.py:29-69`, `generate_oligos.py:311-328`; `frontend/src/services/api.ts:275-303`, `frontend/src/scenes/PAS/components/PASForm/index.tsx:47-55,202-212`; paper p. 361.

### BLOCK E — engineering

**`interface` — YES (two factual errors from the extraction corrected).**
Self-hosted **web GUI** (React + TypeScript + Ant Design, `localhost:3000`) is the primary interface; a **REST API** (FastAPI + Gunicorn, `localhost:8000`) exposes `POST /v1/ssm`, `/v1/qclm`, `/v1/pas`, `POST /v1/species_table`, `GET /v1/get_species`, `GET /v1/{check,cancel,forget,result}/{task_id}`, and the XLSX export endpoints.
1. **There is no `/v1/export_pas` endpoint.** `api/server_fastapi.py` defines only `GET /v1/export_qclm/{task_id}.xlsx`, `/v1/export_ssm/{task_id}.xlsx`, `/v1/export_species_table/{task_id}.xlsx` (the string `"export_pas"` appears only as an `export_name` argument used to build a URL in `get_urls`). All three export endpoints schedule Celery tasks (`tasks.export_ssm`, `tasks.export_qclm`, `tasks.export_species_table`) that **`backend/tasks.py` never defines** — they are dead. Excel export in practice is client-side in `SaveFile/index.tsx`.
2. **"Importable Python package" overstates it.** There is no `setup.py`, `pyproject.toml` or `setup.cfg` anywhere in the 367-path tree, so `mutation_maker` is importable only as a source directory on `PYTHONPATH` and cannot be pip-installed even from a git URL.

**No CLI** (`argparse` code search hits are only `LICENSES_THIRD_PARTY` and `frontend/package-lock.json`). Deployment is Docker Compose or manual (Redis + Celery + Gunicorn + nginx, optional AWS Lambda). The API was migrated post-publication from the paper's hug/falcon server to FastAPI (commit *"Migrate API from hug to FastAPI"*): `api/server.py` imports `app` from `api/server_fastapi.py`, and both `api/Dockerfile` and the `Makefile` run `gunicorn server:app`. `api/api.py` is retained as a legacy file but is no longer the entry point, so citations to it elsewhere in this record describe the superseded implementation (its route inventory is identical to the FastAPI one). Version string in `backend/tasks.py`: `print("Mutation Maker version: 1.0.0")`.
*Footnote relevant to programmatic use:* the SSSM degenerate-codon set-cover solver runs **client-side in TypeScript** (`frontend/src/scenes/SSM/components/amino.ts`: `modifyGrammar`, `joinGrammars`, `generateCombinationsFromDegenerateCodon`, `MAX_NUMBER_OF_COMBINATIONS = 100`, `MAX_DURATION = 10 min`), not in the Python backend — REST API users do not get it.
*Source:* `api/server_fastapi.py:36-159`, `api/server.py:19`, `api/Dockerfile:20`, `backend/tasks.py`, `gh api .../git/trees/master?recursive=1`, `frontend/.../amino.ts`.

**`license` — YES.**
GPL-3.0. `gh api repos/ra100/Mutation_Maker` → `license.spdx_id = "GPL-3.0"`; nearly every source file carries *"Copyright (c) 2020 Merck Sharp & Dohme Corp. … Mutation Maker is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License … either version 3 of the License, or (at your option) any later version."* Coverage is not universal: 66 of the 72 `.ts`/`.tsx` files under `frontend/src` carry the header, but the SSM/QCLM/PAS form components do not, and `backend/mutation_maker/__init__.py` is empty. The article itself is CC-BY.
*Source:* GitHub API; per-file headers in `api/api.py`, `backend/tasks.py`, `backend/mutation_maker/*.py`, `frontend/src/**`.

**`maintained` — YES (reframed: author-maintained, not a community fork).**
The paper's URL `github.com/Merck/Mutation_Maker` 404s. The live repository `github.com/ra100/Mutation_Maker` is maintained by **co-author Rastislav Svarba**: the current head's commit metadata names Rastislav Švarba as committer, and Hiraga et al. list him as an author affiliated with MSD Czech Republic. Commits by month: 2020-07 (9), 2020-11 (3), 2021-01 (4), 2021-02 (10), 2022-08/09/10 (20), 2025-07 (7), 2026-02 (32) — continuous, not a single burst. Last commit 2026-02-14; CI and Primer3 workflows both `success` 2026-02-14; not archived; 0 open issues. Recent work is genuine modernisation (React 18 / antd 4 migration plans in `docs/plans/2026-02-*`, Python 3.11, Biopython `GC`→`gc_fraction` compat shim, Playwright e2e suite). **No tagged releases and no GitHub releases** (`/tags` and `/releases` both empty).
*Source:* official repository metadata, commit history, actions runs, tags and releases; author list of Hiraga et al. 2021.

---

## 3. Capability summary table

| Key | Value | One-line basis |
|---|---|---|
| library_as_object | partial | PAS enumerates fragment-level variant oligos with normalised ratios; XLSX/FASTA/GenBank spec upload; but no full-library object and no enumeration in SSM/QCLM |
| dag_chaining | no | three independent endpoints/Celery tasks; no composition primitive |
| lazy_evaluation | no | eager solve, materialised JSON in Redis; no generator/streaming API |
| mixed_mutagenesis_one_pool | partial | per-site heterogeneous specs + frequencies + automatic WT union (QCLM) / WT remainder (PAS) / codon-level specs; but SSM has one job-wide degenerate codon, workflows cannot be combined, no sampled/random component |
| combinatorial_motif_place | no | motifs are only things to avoid; no element placement; and no indels at all (`self.length = 3`) |
| barcode_generation | no | `barcode` code search → 1 hit, in bootstrap.css |
| per_sequence_provenance | partial | full per-emitted-sequence records (`PASOligoOutput`/`PASMutationFormatted`/`PrimerOutput`); not a compositional build history, no full-length members. **Lenient reading = yes** |
| automatic_naming | yes | `SaveFile/index.tsx` automatic naming for SSM/QCLM/PAS; 96-well plate positions in SSM only |
| design_visualization | yes | neXtProt Feature Viewer fork with GC graph, primer directionality, mutation alignment |
| assay_dms | yes | SSSM is deep-mutational-scan library construction; 7 validated libraries, 31–96 sites |
| assay_mpra | no | all positions codon-indexed; no regulatory-element vocabulary |
| assay_insilico | no | no predictive sequence-function model; Fig. 2 is algorithm benchmarking |
| genome_coordinates | no | GOI-relative codon indices only; `.gb` upload strips all annotation |
| transcript_models | no | `GTF`/`gff` → 0 hits; only codon-usage and genetic-code tables |
| exon_intron_split_codons | no | `exon`/`intron` → 0; codons contiguous by construction |
| hgvs_input | no | bespoke `E52W` parser; HGVS prefixes raise ValueError |
| vcf_vep_output | no | JSON / XLSX / print-PDF only |
| consequence_annotation | no | stop-codon exclusion and a display-only `wild_type` flag are not consequence calls |
| primer_design | yes | four solvers, Primer3 + configurable thermodynamics, wet-lab validated |
| codon_optimization | yes | reverse translation with CFP scoring, 35,799 species; optimised gene returned in output |
| synthesis_constraints | yes | full `PASConfig` length/Tm/GC/hairpin/homodimer/motif constraint set |
| interface | yes | web GUI + REST API; no CLI, no pip-installable package |
| license | yes | GPL-3.0 (MSD), article CC-BY |
| maintained | yes | author-maintained (current head committed by paper co-author R. Svarba), commits through 2026-02-14, CI green |

---

## 4. The tool's own documented examples / vignettes

The repository documentation and Supporting Information (only "Figures S1−S4 illustrate Mutation Maker algorithms") contain no tutorial or vignette. What exists is (a) the paper's experimental validations and (b) shipped test fixtures and executable specs.

### From the paper (sequences NOT disclosed — the paper uses the anonymised labels "Gene X" / "Gene Y" and gives neither the sequences nor a stated reason for withholding them)
1. **SSSM libraries, Table 1** — 7 libraries on Gene X (~700 bp) and Gene Y (~1.2 kbp), 31–96 target sites each, vectors pAVE011-Kan / pET30a(+), hosts *E. coli* W3110 / BL21(DE3); 16–48 clones Sanger-sequenced; mean success 79.9 %.
2. **MSDM/QCLM libraries, Table 2** — 13 designs on Gene X, 13–18 target mutations each, primer:template ratios Low/Medium/High, pET30a(+) / BL21(DE3); mean success 72 %, 2.5 mutations per clone.
3. **PAS de novo synthesis, Figure 3** — ~1.2 kb gene split into 32 overlapping fragments (default settings), 6 fragments carrying 16 mutations at 8 sites, 50 oligos total; assembled in two PCR rounds via 3-way (A1–A3) or 4-way (B1–B4) segment splits; 32/32 clones full-length.
4. **In-silico benchmarks** — 50 simulated genes × 10 mutations, brute force vs. fast approximation (Fig. 2A/B); 10 genes with 2–10 sites and ≤5 mutations/site for MSDM coverage weights (Fig. 2C); 5 random 250-aa proteins reverse-translated, MM-CFP vs. GeneGenie CAI (Fig. 2D); random amino-acid groups, MM degeneracy vs. CodonGenie (Fig. 2E).

### From the repository — these ARE runnable
5. `backend/tests/test_support.py::generate_SSM_input(mut_ind)` — 10 preset SSM mutation sets over a bundled plasmid (`sample_ssm_sequence()`), e.g. `["L50X","E60X","W70X","R80X","E100X","L120X","G160X","G200X","I220X","L240X"]`; `X` = saturate with the job's degenerate codon.
6. `generate_random_SSM_input(...)` (`backend/tests/test_support.py:309`) — randomised site generator, used by the SSM solver tests in `backend/tests/unit_tests/test_ssm.py`; the scaling benchmark `backend/ssm_algorithm_comparison.py` does not call it and defines its own `generate_inputs()`.
7. `generate_qclm_input(ind=0..11)` — 12 preset MSDM cases, e.g. `ind=2`: `E52W/L/F/A, V73L/G/F/I/S, A165P/E/R/K, R169I/K/L/Y` (4 sites × 4–5 AAs); `ind=3`: 21 entries including explicit parental codons (`P269P, A271A, V272V, Q376Q, M379M, K383K, T408T, H411H, L412L`) — note WT is included automatically anyway by `MutationSite.__init__`; `ind=4`: a 13-site cysteine scan (`G310C, H348C, I6C, L9C, S11C, G15C, A16C, A20C, E21C, V26C, G27C, L31C, G32C`) with degeneracy off.
8. `generate_pas_input(ind=1..6)` — PAS cases with per-site frequencies, e.g. `ind=3`: `pos 1 → "R,L" @0.25`, `pos 2 → "R" @0.9`, `pos 10 → "I,V,F,S,R,T" @0.5`, `pos 12 → "Y" @0.1`, … (13 sites). `ind=5` uses a **protein** input sequence (triggers reverse translation).
9. **BDD / Gherkin executable specs** (missed by the original extraction): `backend/features/*.feature` — `degenerate_codon.feature`, `mutation.feature`, `primer.feature`, `temperature_calculator.feature`, `possible_mutagenic_primer_position.feature`, `extract_sequence_from_plasmid.feature`, `primer3_interoperability.feature` — with `backend/features/steps/*.py` (`steps_common.py`, `steps_degeneracy.py`, `steps_extract_sequence.py`, `steps_mutagenic_primer_position.py`, `steps_mutation.py`, `steps_primer.py`, `steps_primer3interoperability.py`, `steps_temperature_calculator.py`).
10. `e2e-tests/tests/app.spec.ts` (32 Playwright UI tests) and `api.spec.ts` (9 API tests) exercising all three workflow forms and endpoints.
11. **XLSX mutation-table input format** for PAS (`Site` / `MT` / `MT%` columns) accepted and validated by `FileUploadMutations.tsx` — a tabular input format defined in code; no template workbook ships in the repository (it contains no `.xlsx`/`.xls` file at all).

**Reproducibility note for PoolParty:** items 1–3 cannot be reproduced literally (gene sequences withheld). Items 5–8 and 11 are concrete and self-contained. The natural PoolParty comparison targets are #7 `ind=2` (4 sites × 4–5 defined substitutions + automatic WT) and #8 (per-site substitution sets with per-site frequencies), which PoolParty can express as a library specification — with the honest note that PoolParty enumerates the members while Mutation Maker designs the oligos.

---

## 5. Notable capabilities NOT covered by the current row list

1. **Degenerate-codon set-cover solvers.** Two distinct heuristics: (a) SSSM — given amino-acid inclusion and exclusion lists, find a minimal degenerate codon/alphabet (runs **client-side in TypeScript**, `frontend/src/scenes/SSM/components/amino.ts`, so REST-API users do not get it); (b) MSDM/PAS — degeneracy restricted to the specified substitutions plus the parental amino acid, ranked by codon-frequency product (`generate_oligos.py::find_candidates_for_set_cover` / `Codons.solve_set_cover` / `pas_degeneracy_recursion.py`). Benchmarked against **CodonGenie** (Fig. 2E). A whole capability axis (degenerate-codon compression) that PoolParty and most surveyed tools lack.
2. **PAS combinatorial variant enumeration with normalised mixing probabilities.** `cartesian_dictionary` (itertools.product over per-site codon dicts, folding frequencies by multiplication) → `generate_oligos_from_combinations` → one `PASOligo(sequence, ratio)` per codon combination, ratios normalised to sum to 1. The closest thing in Mutation Maker to an enumerated pool.
3. **Oligo mixing-ratio calculation to hit a target variant frequency.** Library-composition control expressed *physically* (molar ratios in the reaction) rather than as an enumerated pool; exported as the "Mix ratio" column.
4. **Automatic wild-type inclusion.** QCLM: `new_aminos = frozenset(old_aminos.union(...))` unconditionally. PAS: if the user's per-site frequencies sum to <1, the wild-type codon is added with the remaining probability.
5. **Mutations specifiable as raw DNA codons, not only amino acids** (`PASInput.is_mutations_as_codons`) — direct nucleotide-level control of the exact codon placed.
6. **Reverse-translated full-length gene returned in the output payload** — a designed full-length DNA sequence *is* emitted (one codon-optimised parent, not a library).
7. **Species-specific codon usage at scale** — 35,799 species (CUTG/EBI `.spsum`) × 33 NCBI genetic-code tables, with a rare-codon frequency threshold.
8. **Configurable thermodynamics** — nearest-neighbour model, salt correction, and full buffer composition (Na⁺/K⁺/Tris/Mg²⁺/dNTP/primer conc.) user-settable; hairpin/homo-/heterodimer Tm penalties relative to reaction temperature.
9. **96-well plate layout assignment** in the Excel export (well positions `A01`…, plate index, per-plate sheets, "Scale"/"Purification" columns ready for a synthesis vendor order form).
10. **Constraint-violation reporting with user-tunable weights** — soft vs. hard constraints; the UI reports *which* constraints a returned primer violates and its penalty score (`non_optimality`, `parameters_in_range`).
11. **Restriction-site avoidance from a named-enzyme list** (`frontend/src/shared/motifs.json`, 162 restriction-enzyme names; the repository gives no vendor or catalogue provenance for the list) — synthesis constraint keyed by enzyme name rather than raw sequence.
12. **File-based / tabular specification input** — PAS mutations from `.xlsx` (`Site`/`MT`/`MT%`), sequences from `.fa`/`.fasta`/`.gb`.
13. **Asynchronous job infrastructure** — Celery + Redis, cancel/forget/poll endpoints, optional AWS Lambda fan-out for Primer3.

---

## 6. Stated limitations

Verbatim from the Discussion:
- *"the SSSM brute-force approach is only exhaustive within fairly large Tm windows that slide along the user-specified Tm range. Thus, in essence, not all possible primers across the specified Tm range are interrogated."*
- *"Mutation Maker by design does not optimize for the presence of GC clamps throughout, since this constraint often eliminated better scoring solutions."*
- *"The scoring functions … have not been optimized using exhaustive laboratory testing and therefore are subject to improvement over time."*
- *"Mutation Maker's gene split algorithm can be also further optimized, since only simplified scores were used to evaluate fragments during backtracking as a trade-off for reduced computational time."*
- *"The gene synthesis algorithm may not guarantee an optimal solution for complex cases when it hits the hardcoded timeout for each Tm iteration."*
- *"the algorithms for reverse translation and for degeneracy are solving these problems using heuristics and randomization and therefore may not output the optimal solution."*
- *"the gene synthesis algorithm only checks for secondary structures at the fragment level and does not scan the entire reverse-translated sequence for resulting mRNA secondary structures that may affect translation efficiency."*
- MSDM coverage is not guaranteed: *"The MSDM algorithm may report an incomplete primer set that covers only some of the specified mutation sites"* (coverage enters the score with `mutation_coverage_weight = 160`).

Limitations observed in the code but not stated by the authors:
- **No insertions or deletions of any kind.** `AminoMutation.__init__` hardcodes `self.length = 3`; every mutation is a codon-for-codon substitution. Code search `insertion` → 1 irrelevant hit, `deletion` → 0.
- **Enumeration of variant sequences is fragment-scoped and PAS-only.** Full-length assembled library members are never emitted; SSM and QCLM emit primers only.
- **Dead server-side export endpoints.** `api/server_fastapi.py` (like the legacy `api/api.py`) schedules `tasks.export_ssm` / `tasks.export_qclm` / `tasks.export_species_table`, but `backend/tasks.py` defines only `tasks.ssm`, `tasks.qclm`, `tasks.pas`, `tasks.species_table`. Excel export works only through the React client.
- **Not installable as a package** — no `setup.py`/`pyproject.toml`/`setup.cfg` in the tree.
- **SSSM degeneracy solving is browser-only** — a REST-API-only consumer must supply the degenerate codon(s) themselves.
- **No tagged releases** — no version pinning is possible beyond a commit SHA.

---

## 7. Paper/code discrepancies (footnote material; state them, do not rely on them)

1. **Default degenerate codon.** Paper p. 359: *"a list of mutation sites with their respective degenerate codon (default, "NNK")"*; and Methods: forward `"NNK"` / reverse `"MNN"` primers "with default settings". Code: `ssm_types.py:215` `degenerate_codon = StringProperty(required=True, default="NNS")`.
2. **Per-site vs. job-wide degenerate codon.** The paper's phrase *"a list of mutation sites with their **respective** degenerate codon"* could be read as per-site. The implementation is a **single job-wide field**: `SSMInput.degenerate_codon`, and `frontend/src/scenes/SSM/components/SSMForm.tsx` has one `degenerateCodon` input plus an "Amino acids" tab that calls `getDegenerateCodons(include, avoid)` and writes the result back into that same one field. The field does accept a **comma-separated set** of codons (`degenerate_codon.split(",")` in `ssm.py:780`; `joinGrammars` in `amino.ts` emits e.g. `"NDT,VHG"`), so it is a job-wide codon *set*, not strictly one codon — but it is not per-site. **We follow the shipped implementation's one-field contract.**
3. **Server description is dated.** The hug/falcon server described in the paper has been replaced by FastAPI (`api/server_fastapi.py`, run as `gunicorn server:app`); `api/api.py` is retained in the tree as a legacy file but is not the entry point.
4. **`/v1/export_pas` does not exist**, despite the symmetric naming implied by the other two workflows.

---

## 8. Unresolved disagreements and confidence

**No capability VALUE is disputed.** The adversarial review attempted to falsify every yes/no/partial and could not change a single one. Two verdicts were "understated" (`per_sequence_provenance`, `maintained`) — the values are unchanged and the evidence has been rewritten to give the tool proper credit. No cell is set to `unknown`.

**One flagged judgement call (not a value change, but declare it if challenged):**
- `per_sequence_provenance = partial`. Under a lenient reading of the row ("each emitted sequence carries a machine-readable record of what was applied to it"), Mutation Maker earns a **yes** — `PASOligoOutput` + `PASMutationFormatted` is a complete per-sequence record. We keep `partial` only because it is a state record rather than a compositional build history, and because full-length assembled members are never emitted. Anyone reading the row the lenient way should upgrade this cell.

**Other judgement calls, stated openly:**
- `assay_dms = yes` — defensible and should not be softened; SSSM *is* deep-mutational-scan library construction, validated on 7 real libraries. The caveat (it designs reagents from a user-supplied position list) belongs in the footnote, not in the value.
- `library_as_object = partial` — hinges on whether "library" means the design job and its fragment-level enumerated oligo pool (yes) or a full-length variant-sequence collection (no).
- `mixed_mutagenesis_one_pool = partial` — PAS/QCLM genuinely mix WT + heterogeneous per-site substitution sets with frequencies in one job; they cannot mix *mutagenesis types* (saturation + defined + sampled) in one specification, and SSM's degenerate codon is job-wide.
- `lazy_evaluation = no` — inferred from eager precomputation plus fully materialised Redis-backed JSON results; no generator/iterator API is exposed at any level, though I did not audit every internal solver loop for streaming.

**High confidence:** all Block C rows (code search `hgvs` 0, `vcf` 0, `gtf` 0, `gff` 0, `exon` 0, `intron` 0, `chromosome` 0, `genome` 0, `barcode` 1-in-CSS, independently re-run by the reviewer); `primer_design`, `codon_optimization`, `synthesis_constraints`, `license`, `maintained`, `interface`, `automatic_naming`, `design_visualization`, `dag_chaining`, `combinatorial_motif_place`, `barcode_generation`.

**Three text-level errors from the original extraction, now fixed — do not reintroduce them:**
1. Calling `ra100/Mutation_Maker` a "community copy"/"community-maintained copy". It is maintained by co-author Rastislav Svarba.
2. The absolute "the tool emits primers/oligos, never the enumerated library members". False for PAS.
3. "Inputs are raw sequence strings only" as the Block C justification. The GUI accepts `.fa/.fasta/.gb`; the point is that the GenBank reader discards all annotation.
