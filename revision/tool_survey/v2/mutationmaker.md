# Mutation Maker — v2 capability record

**Tool:** Mutation Maker
**Slug:** `mutationmaker` (v1 slug `mutation_muker`/`mutation_maker`)
**Citation key:** `Hiraga2021yg`
**Paper:** Hiraga et al., "Mutation Maker, An Open Source Oligo Design Platform for Protein Engineering," *ACS Synth Biol* 2021;10(2):357–370. doi:10.1021/acssynbio.0c00542
**Primary source for this re-score:** `final/mutationmaker.md` (adversarially reviewed FINAL record), read in full.
**Live repo:** https://github.com/ra100/Mutation_Maker (author/MSD-maintained; the paper's `github.com/Merck/Mutation_Maker` URL 404s)
**Row set:** v2 (`ROWS_v2.md`) — 33 keys. `hgvs_input` dropped.
**Compiled:** 2026-08-10

**Scope note honored:** Mutation Maker is a purpose-built mutagenesis-library design tool. It is scored on its own terms, and every cell where the memo's evidence supports a stronger value under the v2 wording has been upgraded (see §Changed).

---

## Block A — library specification

### `library_first_class_object` — **partial**
There *is* a typed, machine-retrievable result object, but not a library container the user can transform.
- Typed schematics classes on both sides: `SSMInput/SSMOutput`, `QCLMInput/QCLMOutput`, `PASInput/PASOutput` (`backend/mutation_maker/{ssm,qclm,pas}_types.py`). `PASOutput = {input_data, results: List[PASResult]}`, `PASResult = {fragment, start, end, length, overlap, overlap_Tm, overlap_GC, overlap_length, oligos: List[PASOligoOutput]}`. So the tool does **not** merely write files: results are fetched programmatically via `GET /v1/result/{task_id}` as structured JSON, and the input job object is echoed back in `PASOutput.input_data`.
- **But:** no object representing "the library" that can be inspected, filtered, sliced, or passed onward as an operand; the only collections are per-fragment oligo lists and per-mutation primer lists. There is no importable Python surface either — no `setup.py`/`pyproject.toml`/`setup.cfg` anywhere in the 367-path tree, so `mutation_maker` is usable only as a source directory on `PYTHONPATH`; the primary interface is a self-hosted React GUI. `SSMOutput`/`QCLMOutput` contain no variant-sequence collection at all.
*Source:* `final/mutationmaker.md` §2 (`library_as_object`), §Block E `interface`; `backend/mutation_maker/{qclm,pas,ssm}_types.py`, `pas_output.py`, `api/api.py`.
*Judgement:* the v1 `library_as_object = partial` splits as partial here + `no` on `library_algebra`. Flagged — a strict reading ("hold, transform, pass onward") would make this `no`.

### `composable_operations` — **no**  *(rename of `dag_chaining`; value carried unchanged)*
`api/api.py` defines only `@hug.post('/ssm')`, `/qclm`, `/pas` (+ `/species_table`), each dispatching to a separate Celery task (`tasks.ssm`, `tasks.qclm`, `tasks.pas` in `backend/tasks.py`). No composition, nesting, or chaining primitive, and no way to feed one workflow's output into another. `AGENTS.md`: *"three main workflows: SSM …, QCLM …, PAS …"*. The pipeline is fixed by the tool.
*Source:* `final/mutationmaker.md` §2 `dag_chaining`; `api/api.py`, `backend/tasks.py`, `AGENTS.md`.

### `lazy_generation` — **no**  *(rename of `lazy_evaluation`; value carried unchanged)*
Fully eager. Solvers return materialised JSON (`PASOutput(...).to_json()`, `SSMOutput`, `QCLMOutput`) stored in the Celery/Redis result backend. No generator/iterator/streaming API at any user-facing level; internal products are eagerly listified (`self.concrete_mutations = list(itertools.product(...))`, `mutation.py`). Paper: *"the brute-force algorithm precomputes all possible forward and reverse primers around each mutation site."*
*Source:* `final/mutationmaker.md` §2 `lazy_evaluation`; `backend/mutation_maker/{pas,ssm,qclm,mutation}.py`; paper p. 359.

### `library_algebra` — **no**
No stack/concat/sample/repeat operation over libraries exists anywhere in the tool. The three workflows cannot be combined in one specification; each job is independent, results are keyed by `task_id`, and there is no primitive that takes two designs and returns one. Combining two designs (e.g. an SSSM primer set with a QCLM set, or two PAS jobs) means merging the exported XLSX/JSON by hand or with a user script outside the tool. What Mutation Maker *does* provide instead is **physical** composition control — normalised oligo mixing ratios (`oligo.ratio /= total_ratio`, asserted to sum to 1) so that a single reaction achieves a user-defined variant frequency — which is a wet-lab mixture, not a library-algebra operation on library objects.
*Source:* `final/mutationmaker.md` §2 (`dag_chaining`, `mixed_mutagenesis_one_pool`), §5 items 2–3; `api/api.py`, `backend/mutation_maker/generate_oligos.py`.

### `exhaustive_single_scans` — **yes**
SSM/SSSM is site-scanning **saturation** mutagenesis: at each targeted position the degenerate codon (`NNS` default in code, `NNK` per the paper) covers all 20 amino acids, one site mutated per clone. Paper: SSSM *"can generate a diverse gene library … a set of amino acids in a given protein sequence are subjected to random mutagenesis (at a single site out of defined positions per clone)"*, and the designed primers *"enable parallel and deep mutational scans."* Wet-lab validated on 7 libraries with **31–96 target sites each**, 79.9 % mean success (Table 1). Shipped fixture `generate_SSM_input` uses `["L50X","E60X",…]` where `X` = saturate at that site. QCLM/PAS additionally accept explicit per-site substitution sets (e.g. a 13-site cysteine scan, fixture `ind=4`).
**Caveats (kept out of the value, per the same reasoning that keeps `assay_dms = yes`):** (i) the position list is user-supplied — there is no "scan every position in the ORF" switch; (ii) **substitutions only — Mutation Maker supports no insertions or deletions at all** (`AminoMutation.__init__` hardcodes `self.length = 3`; code search `deletion` → 0), so the "every deletion" half of the row is not met; (iii) the exhaustive single-site set is realised as a degenerate-codon primer rather than as enumerated variant sequences.
*Source:* `final/mutationmaker.md` §2 (`assay_dms`, `mixed_mutagenesis_one_pool`, `combinatorial_motif_place`), §4 items 1, 5, 7; paper pp. 358–359, Table 1; `backend/mutation_maker/{ssm,mutation}.py`.
*Flagged:* a reviewer who reads the row as requiring both substitution *and* deletion scans, or requiring automatic whole-sequence position enumeration, would score `partial`.

### `sampled_random_mutagenesis` — **partial**
Rate-based library composition is first-class, but the tool never draws a random *sample* of variants.
- **For it:** PAS takes a **per-site mutation frequency** (`MT%` column in the uploaded XLSX; fixture `ind=3`: `pos 10 → "I,V,F,S,R,T" @0.5`, `pos 12 → "Y" @0.1`) and computes oligo mixing ratios to hit it. Paper: *"the algorithm computes the oligo ratios in the reaction mixture that are necessary to achieve the user-defined mutation frequency in the final library."* Per-combination probabilities are folded multiplicatively (`cartesian_dictionary(*mutations_on_site_with_prob)`) and normalised. Wild type is auto-added with the leftover probability when a site's frequencies sum to <1 (`wild_type_prob = 1 - sum_of_probabilities`). The physical library that results *is* a random draw at the specified rates, and SSSM's degenerate codon is described by the authors as *"random mutagenesis … at a single site per clone."*
- **Against:** computationally the enumeration is exhaustive-with-weights, not sampled — no `random.sample`/`choice` over a variant space, no draw-N-variants control, no per-base error-rate / error-prone-PCR mutagenesis model, and no sampled component that can be requested as a library element. The v1 record's stated basis for downgrading `mixed_mutagenesis_one_pool` was precisely *"there is no sampled/random-variant or WT-replicate spike-in concept."*
*Source:* `final/mutationmaker.md` §2 `mixed_mutagenesis_one_pool`, §5 items 2–4, §4 item 8; paper p. 362; `backend/mutation_maker/generate_oligos.py`, `pas_types.py`.
*Flagged:* genuinely borderline. `no` is defensible if the row requires an explicit sampling operation; `yes` is defensible if frequency-targeted library composition counts.

### `higher_order_combinatorial` — **yes**
Multiple mutations in the same molecule is a core, wet-lab-validated capability, in two independent workflows.
- **QCLM/MSDM** exists to install several substitutions in one clone: Table 2 reports 13 designs with 13–18 target mutations each, 72 % mean success and **2.5 mutations per clone**; the solver partitions/extends primers over multiple sites with a `mutation_coverage_weight = 160` term.
- **PAS** explicitly enumerates cross-site combinations: `generate_oligos.py::generate_oligos_from_combinations()` iterates `cartesian_dictionary(*mutations_on_site_with_prob)` (= `itertools.product` over per-site codon dictionaries), writes each combination into the fragment DNA via `Codons.replace_codon`, and appends one `PASOligo(sequence, ratio)` per combination. Figure 3: 6 fragments carrying **16 mutations at 8 sites**, 50 oligos total, 32/32 clones full-length.
**Caveat:** the enumerated combinations are fragment-scoped — combinations spanning different fragments are realised physically at assembly rather than enumerated as full-length sequences.
*Source:* `final/mutationmaker.md` §1, §2 (`library_as_object` correction), §5 item 2, §4 items 2–3; `backend/mutation_maker/generate_oligos.py`, `qclm.py`; paper Table 2, Fig. 3.

### `heterogeneous_components_one_library` — **partial**
One job can carry genuinely heterogeneous **per-site specifications**, but not structurally different component *types*.
- **For it:** PAS/QCLM accept a different substitution set *and a different frequency* at every site in one design (fixture: `position=10, "I,V,F,S,R,T" @0.5` alongside `position=12, "Y" @0.1`); mutations may be given as amino acids **or as raw DNA codons** in the same job (`PASInput.is_mutations_as_codons`; `PASMutationSite` docstring "Each mutation is a amino acid IUPAC code, or a DNA codon"), so codon-level and residue-level specs coexist; wild type is unioned into every QCLM site unconditionally (`new_aminos = frozenset(old_aminos.union({m.new_amino for m in mutations}))`, `mutation.py:147`) and added with the residual probability in PAS. A saturated site and a single-defined-substitution site can therefore sit side by side in one specification.
- **Against:** every component is the same structural kind — a codon-for-codon substitution at a codon-indexed site. No indels, no motif/element insertion, no barcode component, no WT-replicate spike-ins as distinct members, and no sampled component. The three workflows (SSSM / MSDM / PAS) cannot be combined in one specification, and SSM applies a **single job-wide** `degenerate_codon` field (`ssm_types.py:215`, `default="NNS"`; the UI's "Amino acids" tab writes back into that same field, though it accepts a comma-separated codon *set*).
*Source:* `final/mutationmaker.md` §2 `mixed_mutagenesis_one_pool`, §5 items 4–5, §7 item 2; `backend/mutation_maker/{mutation,pas_types,ssm_types}.py`, `frontend/src/scenes/SSM/components/SSMForm.tsx`.
*Flagged:* hinges on whether "structurally different component types" means different specification shapes per site (→ `yes`) or different kinds of library component (→ `no`).

### `combinatorial_motif_place` — **no**
"Motif" in Mutation Maker means *sequence to avoid*: `PASConfig.avoided_motifs = ListProperty(str)`, `backend/mutation_maker/motifs.py`, and `frontend/src/shared/motifs.json` is a restriction-enzyme name list (`AarI, AatII, Acc65I, …`). Paper: *"exclusion of hairpins, homodimers and certain motifs."* No element placement, orientation, permutation, or position enumeration. Reinforcing: **no insertions or deletions of any kind** — `AminoMutation.__init__` hardcodes `self.length = 3` (also `ConcreteTripletMutation`, `MutationSite`); code search `insertion` → 1 irrelevant hit, `deletion` → 0. Inserting a motif is structurally impossible.
*Source:* `final/mutationmaker.md` §2 `combinatorial_motif_place`; `backend/mutation_maker/{pas_types,motifs,mutation}.py`, `frontend/src/shared/motifs.json`; paper p. 360.

### `barcode_generation` — **no**
GitHub code search `barcode` over `ra100/Mutation_Maker` → 1 hit, in `frontend/feature-viewer/css/bootstrap.css` (irrelevant CSS). No UMI, edit-distance, or tag-set machinery in code or paper; no barcode attachment step in any of the three workflows.
*Source:* `final/mutationmaker.md` §2 `barcode_generation`; GitHub code search; paper (full text).

### `per_sequence_provenance` — **yes**  *(upgraded from v1 `partial`)*
Every sequence Mutation Maker emits carries a structured, machine-readable record of what was applied to it, well beyond naming the mutation — which is exactly the v2 row's bar ("structured per-sequence metadata recording HOW it was built, **beyond naming the mutation**"). The v1 FINAL record explicitly invites this upgrade: *"a lenient reading of this row would score Mutation Maker 'yes'; if the row is defined as 'each emitted sequence carries a machine-readable record of what was applied to it', Mutation Maker meets it."*
- `PASOligoOutput = {sequence, mix_ratio, mutations, reds, blues}`, where `mutations: List[PASMutationFormatted] = {position, mutated_amino, wild_type_amino, wild_type_codon, mutated_codon, frequency, wild_type}` and `reds`/`blues` are the mutated vs. wild-type codon index lists. So each emitted oligo records, per site, the exact codon replaced, the exact codon installed, the residue change, and the intended frequency.
- The enclosing `PASResult` adds `{fragment, start, end, length, overlap, overlap_Tm, overlap_GC, overlap_length}` — where in the parent the sequence came from.
- Primers likewise: `PrimerOutput{sequence, start, length, temperature, gc_content, degenerate_codons, overlap_with_following}`; `SSMMutationOutput{mutation, non_optimality, parameters_in_range, result_found, forward_primer, reverse_primer, overlap}` — including *which constraints were violated* in producing it.
- It survives export: XLSX headers carry "Mutation Syntax", "Mix ratio", "Target and MT%", plate/well position (`frontend/src/shared/components/SaveFile/index.tsx`).
**Residual caveats (do not change the value under the v2 wording):** it is a *state* record (what mutation is present, at what frequency) rather than a *compositional build history* (no record of which design operations were composed), and no records exist for full-length assembled members because those are never emitted.
*Source:* `final/mutationmaker.md` §2 `per_sequence_provenance` + §8 flagged judgement call; `backend/mutation_maker/pas_types.py:158-196`, `pas_output.py`, `ssm_types.py`; `frontend/src/shared/components/SaveFile/index.tsx`.
*Flagged:* value changed from v1. If the row is scored strictly as "compositional build history", revert to `partial`.

### `automatic_naming` — **yes**
Names are generated automatically and are plate/well aware. `frontend/src/shared/components/SaveFile/index.tsx`: `PLATE_SIZE = 96`, `getPosition` yields `A01`-style wells, prefix from `input.config.oligo_prefix || 'prefix'`:
- SSM `` `${oligoPrefix}-${Math.floor(index/PLATE_SIZE)+1}F-${getPosition(index%PLATE_SIZE)}` `` (and `…R-…`)
- QCLM `` `${oligoPrefix}-${siteIndex+1}${site.length>1 ? String.fromCharCode(CHAR_LOWER_A+index) : ''}` ``
- PAS `` `${oligoPrefix}-Fr${index+1}-${oIndex+1}` ``
A separate "Mutation Syntax" column carries `m.mutations.map(prop('identifier')).join(',')` plus degenerate codons, so names plus adjacent columns are informative and vendor-order-ready.
*Source:* `final/mutationmaker.md` §2 `automatic_naming`; `frontend/src/shared/components/SaveFile/index.tsx:35,48,55-57,136,142,183,301,314`.

### `design_visualization` — **yes**
Paper, Methods: *"A JavaScript (JS) User Interface (UI) where users can submit oligo design tasks and visualize the results using a customized version of the neXtProt viewer"*; *"The following features were added: GC content graph, primer directionality, mutation alignment, and integration with interactive results table and with React user interface."* Repo: `frontend/feature-viewer` (fork of `calipho-sib/feature-viewer`, adds an `arrow` object type for primers) plus `SSMFeatureViewer.tsx`, `QCLMFeatureViewer.tsx`, `PASFeatureViewer.tsx`. Results export to PDF or a temporary shareable link.
**Caveat kept:** it visualises the primer/fragment layout on the gene, not a design graph or library-composition plot.
*Source:* `final/mutationmaker.md` §2 `design_visualization`; paper Methods p. 366; `frontend/feature-viewer/`, `frontend/src/scenes/*/components/*FeatureViewer.tsx`.

---

## Block B — assay coverage

### `assay_dms` — **yes**
Core purpose, and the v1 record marks softening this cell as the easiest thing for an author-referee to rebut. SSSM designs deep-mutational-scan libraries; the paper says the primers *"enable parallel and deep mutational scans."* Validated on 7 real DMS-style libraries with 31–96 target sites each, 79.9 % mean success (Table 1), plus 13 MSDM designs (Table 2). Caveat belongs in the footnote, not the value: it designs the *reagents* from a user-supplied position list and does not enumerate the resulting full-length variants.
*Source:* `final/mutationmaker.md` §2 `assay_dms`, §8; paper pp. 358–359, Table 1; `backend/mutation_maker/ssm.py`.

### `assay_mpra` — **no**
Every mutation entry point is codon-indexed: `parse_codon_mutation` computes `zero_based_base_position = (one_based_codon_position - 1) * 3`; `AminoMutation` validates against IUPAC one-letter amino acids; PAS filters mutations onto fragments via `goi_offset + mut[0]*3`. No promoter/enhancer/regulatory vocabulary, no TF-motif insertion, no scrambles or shuffles anywhere in paper or repo (re-verified: code search `shuffle` → 0, `scramble` → 0). *Honest hedge:* PAS with `is_dna_sequence=True` will synthesise arbitrary DNA, so non-coding synthesis is technically possible — but mutation positions remain codon-indexed and there is no regulatory-element concept.
*Source:* `final/mutationmaker.md` §2 `assay_mpra`; `backend/mutation_maker/{mutation,pas_types}.py`, `generate_oligos.py`; GitHub code search 2026-08-10.

### `assay_insilico` — **no**
No predictive model, scoring model, or model-in-the-loop code anywhere in `backend/mutation_maker/`. The paper's in-silico work (Fig. 2) is *algorithm* benchmarking — runtime, penalty scores, MM-CFP vs. GeneGenie CAI, degeneracy vs. CodonGenie — not sequence-function prediction. Experimental validation is PCR + Sanger.
*Source:* `final/mutationmaker.md` §2 `assay_insilico`; paper Fig. 2 and Results; repo file listing.

---

## Block C — genomics integration (all no)

All verified in the paper and by GitHub code search over `ra100/Mutation_Maker`: `hgvs` 0, `vcf` 0, `gtf` 0, `gff` 0, `exon` 0, `intron` 0, `chromosome` 0, `genome` 0 — independently re-run by the adversarial reviewer.

### `genome_coordinates` — **no**
Inputs are sequences plus intra-gene indices: `SSMSequences{forward_primer, reverse_primer, gene_of_interest, five_end_flanking_sequence, three_end_flanking_sequence, plasmid}`, `Plasmid{gene_start_in_plasmid, gene_end_in_plasmid, plasmid_sequence}`; positions are 1-based **codon** indices in the GOI. Pre-empting the "but it reads GenBank" rebuttal: `frontend/src/shared/components/FileUploadInput/index.tsx` accepts `.fa, .fasta, .gb`, but the "GenBank parser" is `text.split('ORIGIN')[1].split('').filter(letter => /[ATCGatcg…]/.test(letter)).join('').toUpperCase()` — it discards LOCUS, FEATURES, every coordinate and every annotation. No genomic coordinate ever enters the system.
*Source:* `final/mutationmaker.md` §2 `genome_coordinates`; `backend/mutation_maker/ssm_types.py:35-38,96-102`; `frontend/src/shared/components/FileUploadInput/index.tsx:25,35-37`.

### `transcript_models` — **no**
`GTF` → 0, `gff` → 0. The only external data files are codon-usage `.spsum` tables, `species.table`/`speciesInfo.table`, and NCBI genetic-code tables. `.gb` upload strips all annotation, so no transcript model is ever ingested.
*Source:* `final/mutationmaker.md` §2 `transcript_models`; GitHub code search; `backend/mutation_maker/codon_usage_data/`.

### `exon_intron_split_codons` — **no**
`exon` → 0, `intron` → 0. `AminoMutation.__init__` sets `self.length = 3` unconditionally, `parse_codon_mutation` computes `zero_based_base_position = (one_based_codon_position - 1) * 3`, and `mutations_on_fragments` assumes `goi_offset + mut[0]*3` contiguity. Codons are contiguous by construction with no escape hatch.
*Source:* `final/mutationmaker.md` §2; `backend/mutation_maker/mutation.py:28,67`, `generate_oligos.py:60-65`.

### `vcf_vep_output` — **no**
`vcf` → 0. Outputs are JSON over REST, client-side XLSX (`exceljs`/`xlsx` in `SaveFile/index.tsx`), and browser-print PDF.
*Source:* `final/mutationmaker.md` §2; `api/api.py`, `frontend/src/shared/components/SaveFile/index.tsx`.

### `consequence_annotation` — **no**
Nothing labels stop-gained, synonymous, missense, or frameshift. The two near-misses are not consequence calls: stop codons are *excluded by default* from degeneracy solutions (paper: *"only stop codons are excluded by default"*; `defaultAvoid = [...aminoToCodonMap.STOP]` in `frontend/src/scenes/SSM/components/amino.ts`), and `PASMutationFormatted.wild_type` only drives `reds`/`blues` colouring in `pas_output.py`.
*Source:* `final/mutationmaker.md` §2; paper p. 360; `backend/mutation_maker/pas_types.py:165`, `pas_output.py`.

*(`hgvs_input` dropped per ROWS_v2. For the record, it was `no`: bespoke `E52W` parser, `p.Glu52Trp` and `c.155A>G` both raise `ValueError`.)*

---

## Block D — adjacent / complementary

### `primer_design` — **yes**
The entire point of the tool, and if anything under-sold. Two SSSM algorithms (brute force with Primer3, and a default "fast-approximation" primer-growth algorithm), an MSDM partition/extension algorithm, and a PAS constraint-satisfaction + backtracking gene splitter: `ssm.py`, `ssm_fast_approximation.py`, `qclm.py`, `qclm_solution.py`, `pas.py`, `pas_back_track.py`, `pas_solution.py`, `primer_scoring.py`, `primer_filtering.py`, `primer3_interoperability.py`, `temperature_calculator.py`, `site_split.py`, plus an AWS Lambda Primer3 fan-out (`lambda/`). Constraints: primer length, 5′/3′ end length, 3′-end Tm, overlap size/Tm, GC content, GC clamp, hairpin/homodimer/heterodimer Tm, Tm consistency; weighted-least-squares scoring with user-adjustable weights and explicit reporting of violated constraints (`non_optimality`, `parameters_in_range`). Tm via `primer3-py` with configurable model (SantaLucia 1997/2004, Breslauer 1986, Sugimoto 1996; Owczarzy/SantaLucia salt corrections; Na/K/Tris/Mg/dNTP/primer concentrations). Wet-lab validated: SSSM 79.9 % across 7 libraries, MSDM 72 % at 2.5 mutations/clone, PAS 1.2 kb gene 32/32 clones full-length.
*Source:* `final/mutationmaker.md` §2 `primer_design`; repo `backend/mutation_maker/*`, `lambda/`; paper Tables 1–2, Fig. 3.

### `codon_optimization` — **yes**
`reverse_translation.py`, `codon_usage_table.py`, `usage_table.py`, `translation_scoring.py`, `degenerate_codon.py`. Paper: *"Mutation Maker's reverse translation algorithm generates a DNA sequence from a set of randomly selected codons of the host organism (excluding rare codons) … the algorithm computes the codon frequencies product (CFP) which is proportional to codon adaptation index (CAI)"*, benchmarked against GeneGenie's CAI optimizer (Fig. 2D). Ships codon usage for **35,799 species** and **33 taxa-specific genetic code tables**; `PASConfig.organism = "e-coli"`, `codon_usage_frequency_threshold = 0.1` (rare-codon cutoff). The optimised gene is **delivered to the user**, not just used internally: `pas.py::find_solution` calls `sequences.translate_goi_sequences(translator)`, mutating `workflow_input.sequences.gene_of_interest` in place, and `pas_solve` returns `PASOutput(input_data=workflow_input, results=results)`.
*Source:* `final/mutationmaker.md` §2 `codon_optimization`; `backend/mutation_maker/pas.py:83-99,160-166`, `pas_types.py:79-86`, `reverse_translation.py`; paper pp. 361, 366, Fig. 2D.

### `synthesis_constraints` — **yes**
Verified directly in `PASConfig`: `min/max/opt_oligo_size = 40/90/56`, `min/max/opt_overlap_tm = 50/65/56`, `min/max/opt_overlap_length = 15/25/21`, `min/max_gc_content = 40/60`, `avoided_motifs`, `codon_usage_frequency_threshold = 0.1`, `safe_temp_difference = 10`, `hairpin_homodimer_weight = 2`, plus a full `TemperatureConfig`. Paper: *"a constraint-satisfaction algorithm with a backtracking component that splits a given gene sequence into an even number of overlapping DNA fragments that conform to length, Tm, mutation location, GC content, hairpins, homodimers and motif exclusion constraints."* Restriction-site avoidance is keyed by **enzyme name** (`frontend/src/shared/motifs.json`, NEB names). **Caveat kept:** these constrain the oligos/fragments, not final assembled library members.
*Source:* `final/mutationmaker.md` §2 `synthesis_constraints`, §5 item 11; `backend/mutation_maker/pas_types.py:29-69`; paper p. 361.

### `degenerate_iupac_codons` — **yes** *(new row; strongly supported — a headline capability)*
Two distinct degenerate-codon set-cover solvers, both directions of the problem (design **and** expansion), plus published benchmarking:
- **SSSM direction (design/compression):** given amino-acid inclusion and exclusion lists, find a minimal degenerate codon or codon set. `frontend/src/scenes/SSM/components/amino.ts`: `getDegenerateCodons(include, avoid)`, `modifyGrammar`, `joinGrammars` (emits e.g. `"NDT,VHG"`), `generateCombinationsFromDegenerateCodon` (expansion), `MAX_NUMBER_OF_COMBINATIONS = 100`, `MAX_DURATION = 10 min`, `defaultAvoid = [...aminoToCodonMap.STOP]`.
- **MSDM/PAS direction:** degeneracy restricted to the specified substitutions plus the parental amino acid, ranked by codon-frequency product — `generate_oligos.py::find_candidates_for_set_cover`, `Codons.solve_set_cover`, `pas_degeneracy_recursion.py`, `degenerate_codon.py`; `backend/features/degenerate_codon.feature` is an executable Gherkin spec for it.
- Job input is itself a degenerate IUPAC codon field: `SSMInput.degenerate_codon = StringProperty(required=True, default="NNS")` (paper says `NNK`), accepting a comma-separated codon set (`degenerate_codon.split(",")`, `ssm.py:780`); outputs report `PrimerOutput.degenerate_codons`.
- Benchmarked against **CodonGenie** (Fig. 2E). The v1 record calls this *"a whole capability axis (degenerate-codon compression) that PoolParty and most surveyed tools lack."*
**Caveat:** the SSSM set-cover solver runs client-side in TypeScript, so REST-API-only consumers must supply the degenerate codon(s) themselves.
*Source:* `final/mutationmaker.md` §5 item 1, §2 `interface` footnote, §7 items 1–2; `frontend/src/scenes/SSM/components/amino.ts`, `backend/mutation_maker/{degenerate_codon,generate_oligos,pas_degeneracy_recursion}.py`; paper Fig. 2E.

### `negative_control_generation` — **no** *(new row)*
No scramble, shuffle, reverse, or reverse-complement control generation. Re-verified by targeted GitHub code search over `ra100/Mutation_Maker` on 2026-08-10: `shuffle` → 0 hits, `scramble` → 0 hits; the v1 record independently states *"no scrambles or shuffles anywhere in paper or repo"*, and since no insertions/deletions/arbitrary-sequence substitutions are possible (`self.length = 3`, codon-for-codon only), a scrambled control cannot be expressed as a mutation either.
**Closest analogue, stated for fairness:** **automatic wild-type inclusion** is first-class — QCLM unions the parental amino acid into every site unconditionally (`new_aminos = frozenset(old_aminos.union(...))`, `mutation.py:147`) and PAS adds the WT codon with `wild_type_prob = 1 - sum_of_probabilities`. The parental sequence is the natural negative control in protein engineering, so a reviewer could argue `partial` here. We score `no` because the row asks for *generated* randomization controls, and WT is retained rather than generated.
*Source:* `final/mutationmaker.md` §2 (`assay_mpra`, `mixed_mutagenesis_one_pool`), §5 item 4, §6; GitHub code search 2026-08-10.
*Flagged.*

### `ml_model_in_loop` — **no** *(new row)*
No predictive model anywhere. Nothing in `backend/mutation_maker/` scores sequences with a learned model; the paper's stated purpose is wet-lab oligo design and its validation is PCR/Sanger. The tool's "scoring" is a hand-specified weighted-least-squares penalty over thermodynamic/geometric primer constraints (`primer_scoring.py`, user-tunable weights) plus codon-frequency products — deterministic objective functions, not model predictions. Fig. 2 is algorithm benchmarking, not sequence-function prediction.
*Source:* `final/mutationmaker.md` §2 `assay_insilico`, §5 items 8, 10; `extractions/mutationmaker.md` (*"No model-in-the-loop functionality … Nothing about scoring sequences with predictive models"*).

### `readout_analysis` — **no** *(new row)*
No analysis of sequencing readout from the built library. Targeted GitHub code search over `ra100/Mutation_Maker` on 2026-08-10: `fastq` → 0 hits, `demultiplex` → 0 hits; the repo's 367-path tree contains no read-processing, counting, variant-calling, or enrichment-scoring module (module inventory in `final/mutationmaker.md` §2 is entirely design-side: solvers, primer scoring, thermodynamics, codon usage, degeneracy). Validation sequencing in the paper (24–48 Sanger clones per library, Tables 1–2) was done by the authors with external tooling; no analysis code ships. Outputs are design artefacts only: JSON, XLSX oligo order sheets, PDF.
*Source:* GitHub code search 2026-08-10; `final/mutationmaker.md` §2, §4 items 1–3, §Block E.

---

## Block E — engineering and availability (value `yes`; substance in evidence)

### `interface` — **yes**
Self-hosted **web GUI** (React + TypeScript + Ant Design, `localhost:3000`) is the primary interface; a **REST API** (`hug`/`falcon` + Gunicorn, `localhost:8000`) exposes `POST /v1/ssm`, `/v1/qclm`, `/v1/pas`, `POST /v1/species_table`, `GET /v1/get_species`, `GET /v1/{check,cancel,forget,result}/{task_id}`. `api/server_fastapi.py` now also exists (post-publication addition). **No CLI** (`argparse` hits only in `LICENSES_THIRD_PARTY` and `frontend/package-lock.json`). **No importable/installable Python package** — no `setup.py`/`pyproject.toml`/`setup.cfg` in the 367-path tree, so `mutation_maker` works only as a source dir on `PYTHONPATH`. Two dead ends to note: the three `GET /v1/export_*.xlsx` endpoints schedule Celery tasks that `backend/tasks.py` never defines (Excel export is in practice client-side), and there is no `/v1/export_pas` endpoint at all. Version string: `print("Mutation Maker version: 1.0.0")`.
**Summary value: web GUI + REST API; no CLI, no library API.**
*Source:* `final/mutationmaker.md` §2 `interface`; `api/api.py:40-152`, `backend/tasks.py`, `gh api .../git/trees/master?recursive=1`.

### `license` — **yes**
**GPL-3.0.** Re-verified 2026-08-10: `gh api repos/ra100/Mutation_Maker` → `license.spdx_id = "GPL-3.0"`, `archived: false`. Every source file carries *"Copyright (c) 2020 Merck Sharp & Dohme Corp. … free software … under the terms of the GNU General Public License … version 3 … or (at your option) any later version."* The article itself is CC-BY. Note for the survey: GPL-3.0 is copyleft, unlike the permissive licenses of most tools in this table.
*Source:* GitHub API 2026-08-10; per-file headers in `api/api.py`, `backend/tasks.py`, `backend/mutation_maker/*.py`, `frontend/src/**`.

### `installable_today` — **yes**
**Deployable today via Docker Compose from source; not installable as a package.** Re-verified against the live tree on 2026-08-10: `docker-compose.yml` at root, `Makefile`, and per-service Dockerfiles (`api/Dockerfile`, `backend/Dockerfile`, `frontend/Dockerfile`, `webserver/Dockerfile`, `lambda/Dockerfile`, `lambda/primer3-build-env/Dockerfile`) with `api/requirements.txt`, `backend/requirements.txt`, `lambda/requirements.txt`. The full stack is Redis + Celery + Gunicorn + nginx, with an optional AWS Lambda Primer3 fan-out. CI and the Primer3 workflow were both green on 2026-02-14, and a 41-test Playwright/API e2e suite passes, so the Docker path is live rather than nominal.
**Qualifications:** (i) **no PyPI or conda package** (`mutation-maker`, `mutation_maker`, `mutationmaker` all 404) and no `setup.py`/`pyproject.toml`/`setup.cfg` anywhere, so not pip-installable even from a git URL; (ii) **no tagged releases and no GitHub releases** (`/tags`, `/releases` both empty) — version pinning is only by commit SHA; (iii) **no public hosted instance**; (iv) the repository URL printed in the paper, `github.com/Merck/Mutation_Maker`, **404s** — you must know to use `github.com/ra100/Mutation_Maker`, maintained by co-author Rastislav Švarba (MSD IT).
*Source:* `gh api repos/ra100/Mutation_Maker/git/trees/master?recursive=1` and `repos/ra100/Mutation_Maker` (2026-08-10); `final/mutationmaker.md` §0, §2 `maintained`, §6.
*Flagged:* value depends on whether "installable today" admits self-hosted Docker Compose from an unreleased default branch.

### `last_activity` — **yes**
**Last commit 2026-02-14** on `github.com/ra100/Mutation_Maker` (repo `pushed_at = 2026-02-14T20:01:33Z`, re-verified 2026-08-10); CI and Primer3 workflows both `success` the same day; not archived; 0 open issues; 19 stars, 13 forks; **0 tags, 0 releases**. Maintenance is continuous rather than a single burst — commits by month: 2020-07 (9), 2020-11 (3), 2021-01 (4), 2021-02 (10), 2022-08/09/10 (20), 2025-07 (7), 2026-02 (32) — and recent work is genuine modernisation (React 18 / antd 4 migration plans in `docs/plans/2026-02-*`, Python 3.11, Biopython `GC`→`gc_fraction` compat shim, Playwright e2e suite).
**Must be stated correctly:** this is **author- and copyright-holder-maintained**, not a community fork. `gh api users/ra100` → `{"login":"ra100","name":"Rastislav Švarba","company":"MSD IT"}`, the 8th listed author of Hiraga et al. 2021; co-authors `spalemartin` (Martin Spale) and `prihoda` (David Prihoda) also contribute, and the fork network contains `prihoda/` and `swiewiora/`. The paper's `Merck/Mutation_Maker` URL 404s.
*Source:* GitHub API 2026-08-10 (`repos/ra100/Mutation_Maker`, `users/ra100`, commits, actions runs, tags, releases); `final/mutationmaker.md` §0, §2 `maintained`.

### `documented_examples` — **yes**
**No tutorial, vignette, or docs site; ~25 runnable preset input fixtures plus 7 executable Gherkin specs and 41 e2e tests.** Breakdown:
- **Narrative docs: none.** No docs site, no GitHub Pages (`has_pages = false`); documentation is `README.md`, `backend/README.md`, `AGENTS.md`, `e2e-tests/README.md`. Journal Supporting Information is only *"Figures S1−S4 illustrate Mutation Maker algorithms (PDF)"* — no user guide or tutorial.
- **Runnable shipped examples:** `backend/tests/test_support.py::generate_SSM_input(mut_ind)` — **10** preset SSM mutation sets over a bundled plasmid (e.g. `["L50X","E60X","W70X","R80X","E100X","L120X","G160X","G200X","I220X","L240X"]`); `generate_qclm_input(ind=0..9)` — **10** preset MSDM cases (e.g. `ind=2`: `E52W/L/F/A, V73L/G/F/I/S, A165P/E/R/K, R169I/K/L/Y`; `ind=4`: a 13-site cysteine scan); `generate_pas_input(ind=1..5)` — **5** PAS cases with per-site frequencies (`ind=5` uses a protein input, triggering reverse translation); `generate_random_SSM_input(...)` for the scaling benchmarks (`backend/ssm_algorithm_comparison.py`).
- **Executable specs:** **7** Gherkin `.feature` files re-verified in the live tree (`degenerate_codon`, `mutation`, `primer`, `temperature_calculator`, `possible_mutagenic_primer_position`, `extract_sequence_from_plasmid`, `primer3_interoperability`) with `backend/features/steps/*.py`.
- **e2e:** `e2e-tests/tests/app.spec.ts` (32 Playwright UI tests) + `api.spec.ts` (9 API tests) exercising all three workflow forms and endpoints.
- **Input template:** an XLSX PAS mutation-table format (`Site`/`MT`/`MT%`) accepted by `FileUploadMutations.tsx`; sequences via `.fa`/`.fasta`/`.gb`.
- **Paper worked examples (not reproducible literally — "Gene X"/"Gene Y" sequences withheld as proprietary Merck genes):** 7 SSSM libraries (Table 1), 13 MSDM designs (Table 2), 1 PAS de novo synthesis (Fig. 3), 4 in-silico benchmarks (Fig. 2).
*Source:* `final/mutationmaker.md` §0, §4; `gh api repos/ra100/Mutation_Maker/git/trees/master?recursive=1` (2026-08-10).

---

## Summary table (v2)

| Block | Key | Value |
|---|---|---|
| A | library_first_class_object | partial |
| A | composable_operations | no |
| A | lazy_generation | no |
| A | library_algebra | no |
| A | exhaustive_single_scans | yes |
| A | sampled_random_mutagenesis | partial |
| A | higher_order_combinatorial | yes |
| A | heterogeneous_components_one_library | partial |
| A | combinatorial_motif_place | no |
| A | barcode_generation | no |
| A | per_sequence_provenance | **yes** (upgraded) |
| A | automatic_naming | yes |
| A | design_visualization | yes |
| B | assay_dms | yes |
| B | assay_mpra | no |
| B | assay_insilico | no |
| C | genome_coordinates | no |
| C | transcript_models | no |
| C | exon_intron_split_codons | no |
| C | vcf_vep_output | no |
| C | consequence_annotation | no |
| D | primer_design | yes |
| D | codon_optimization | yes |
| D | synthesis_constraints | yes |
| D | degenerate_iupac_codons | **yes** (new) |
| D | negative_control_generation | no (new) |
| D | ml_model_in_loop | no (new) |
| D | readout_analysis | no (new) |
| E | interface | web GUI + REST API; no CLI, no library API |
| E | license | GPL-3.0 (MSD); article CC-BY |
| E | installable_today | Docker Compose from source; no pip/conda, no releases |
| E | last_activity | 2026-02-14 (author-maintained `ra100/Mutation_Maker`); paper URL 404s |
| E | documented_examples | 0 tutorials; ~25 preset fixtures + 7 Gherkin specs + 41 e2e tests |

---

## How the v1 rows map into v2

| v1 row | v1 value | v2 outcome |
|---|---|---|
| `library_as_object` | partial | SPLIT → `library_first_class_object` = partial; `library_algebra` = **no** (combining designs requires an external script) |
| `dag_chaining` | no | RENAME → `composable_operations` = no (unchanged) |
| `lazy_evaluation` | no | RENAME → `lazy_generation` = no (unchanged) |
| `mixed_mutagenesis_one_pool` | partial | SPLIT → `exhaustive_single_scans` = **yes**; `sampled_random_mutagenesis` = partial; `higher_order_combinatorial` = **yes**; `heterogeneous_components_one_library` = partial |
| `per_sequence_provenance` | partial | **yes** — the v2 wording ("structured metadata recording how it was built, beyond naming the mutation") is the lenient reading the v1 record explicitly said would earn a yes |
| `hgvs_input` | no | DROPPED |
| `maintained` | yes | superseded by `installable_today` + `last_activity` (+ `documented_examples`) |
| all others | — | unchanged |

**Net effect of the v2 restructure on Mutation Maker: it scores higher.** The split of `mixed_mutagenesis_one_pool` turns one hedged `partial` into two clear `yes` values (exhaustive saturation scans; higher-order combinatorial), the new `degenerate_iupac_codons` row captures a headline capability the v1 rows missed entirely, and `per_sequence_provenance` upgrades to `yes`. The one place the restructure is stricter is `library_algebra`, which is an honest `no`.

## Cells flagged for human review

`library_first_class_object` (partial vs. no), `exhaustive_single_scans` (yes vs. partial — no deletions, user-supplied positions), `sampled_random_mutagenesis` (partial; genuinely between no and yes), `heterogeneous_components_one_library` (partial; depends on reading of "structurally different component types"), `per_sequence_provenance` (upgraded to yes; revert to partial if the row means compositional build history), `negative_control_generation` (no; automatic WT inclusion could argue partial), `installable_today` (Docker-from-source counts?), `documented_examples` (how to count test fixtures as "examples"), `assay_dms` (kept yes deliberately — do not soften; easiest cell for an author-referee to rebut).
