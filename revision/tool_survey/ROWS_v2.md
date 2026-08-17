# Capability matrix rows, v2

> **SUPERSEDED by `ROWS_v3.md`.** Kept for provenance. Three of the eight
> argument-critical rows defined here were found defective by the per-row
> verification pass (`verified/`) and are reworded in v3:
> `library_first_class_object` (re-merged two other rows' criteria),
> `composable_operations` ("no" asserted something false for pydna/Biopython/SeqPro),
> and `lazy_generation` (no tool reached "yes"; replaced by
> `unmaterialized_library_addressing`). Verified cell values are in
> `MATRIX_verified.md`.

Restructured from the v1 24-row draft using the completed 12-tool survey.
Changes are justified by measured discriminating power, not by preference.

## What changed and why

| Change | Row(s) | Evidence |
|---|---|---|
| **DROP** | `hgvs_input` | "no" for all 12 tools including VaLiAnT. Zero discriminating power. Move to limitations prose. Also: do NOT concede HGVS to VaLiAnT in the response letter — it does not have it. |
| **SPLIT** | `mixed_mutagenesis_one_pool` -> A5/A6/A7/A8 | VaLiAnT scored "yes", which hides the gap the paper turns on. Its mutators are *exhaustive single-event scans only* — no sampled/random mode, no pairwise or higher-order mutator (code-verified in `src/valiant/mutator_type.py`, a closed 7-member Enum). As one row this reads as parity; split, it shows the real boundary. |
| **SPLIT** | `library_as_object` -> A1/A4 | 8 of 12 scored "partial" — the row was too coarse to discriminate. VaLiAnT's own paper says the final concatenated library is assembled by the user outside the tool, so "can hold a library object" and "can combine libraries" are distinct claims. |
| **ADD** | D4-D7 | Recurring capabilities across `additional_capabilities` that the v1 rows missed entirely. |
| **ADD** | E3-E5 | 2 of 3 named competitors are not installable. Availability is a legitimate, checkable comparison axis and the referees asked about program features and outputs. |

## Row set

### Block A — library specification (differentiating block)

| Key | Question |
|---|---|
| `library_first_class_object` | Is the library an object the user can hold, inspect, transform and pass onward — as opposed to a tool that only writes files? |
| `composable_operations` | Do design steps compose/chain/nest arbitrarily, or is the pipeline fixed by the tool? |
| `lazy_generation` | Are sequences produced on demand rather than fully materialized? |
| `library_algebra` | Can whole libraries be combined/sampled/repeated as operations (stack, concat, sample), inside the tool? |
| `exhaustive_single_scans` | Exhaustive single-event mutagenesis (every substitution/deletion at every position) |
| `sampled_random_mutagenesis` | Random or rate-based sampling of variants |
| `higher_order_combinatorial` | Pairwise and higher-order combinations of mutations |
| `heterogeneous_components_one_library` | Can structurally DIFFERENT component types (e.g. exhaustive singles + sampled higher-order + WT replicates) coexist in one library specification? |
| `combinatorial_motif_place` | Combinatorial placement of motifs (multiple motifs, positions, orientations, permutations) |
| `barcode_generation` | Constraint-aware barcode generation and attachment |
| `per_sequence_provenance` | Structured per-sequence metadata recording HOW it was built, beyond naming the mutation |
| `automatic_naming` | Informative sequence names generated automatically |
| `design_visualization` | Visual inspection of design or library (graph view, highlighted sequences) |

### Block B — assay coverage

`assay_dms` · `assay_mpra` · `assay_insilico`

### Block C — genomics integration (PoolParty expected "no" throughout)

`genome_coordinates` · `transcript_models` · `exon_intron_split_codons` ·
`vcf_vep_output` · `consequence_annotation`

### Block D — adjacent / complementary

| Key | Question |
|---|---|
| `primer_design` | Mutagenic primer / oligo design for wet-lab protocols |
| `codon_optimization` | Codon usage optimization |
| `synthesis_constraints` | Constraint checking/enforcement on individual sequences |
| `degenerate_iupac_codons` | **NEW** — degenerate/IUPAC codon design or expansion (present in MPRAnator, CodonGenie, Mutation Maker, Oligopool Calculator) |
| `negative_control_generation` | **NEW** — scrambles, shuffles, reverse/complement controls as a first-class feature |
| `ml_model_in_loop` | **NEW** — design driven by a predictive model's output (ledidi, tangermeme) |
| `readout_analysis` | **NEW** — analysis of sequencing readout from the built library (Oligopool Calculator). Marks complementarity, not competition. |

### Block E — engineering and availability

| Key | Question |
|---|---|
| `interface` | Python API / R package / CLI / web / GUI |
| `license` | License |
| `installable_today` | **NEW** — installable via pip/conda/Docker today, or not at all |
| `last_activity` | **NEW** — last release or commit observed |
| `documented_examples` | **NEW** — count of shipped examples/vignettes |

## Note on main-text vs. supplement

This full set is the *survey* instrument. The main-text table should be a
general subset — Blocks A (collapsed), B, and C, plus one availability row.
Blocks D and E belong in supplement or in the response letter. Decide after
PoolParty's column is filled.
