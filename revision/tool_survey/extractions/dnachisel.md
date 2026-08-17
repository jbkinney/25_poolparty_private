# DNA Chisel — evidence memo

**Tool:** DNA Chisel (`dnachisel`)
**Authors / group:** Valentin Zulkower & Susan Rosser, Edinburgh Genome Foundry (EGF), University of Edinburgh
**Paper:** Zulkower V, Rosser S. "DNA Chisel, a versatile sequence optimizer." *Bioinformatics* 36(16):4508–4509 (2020). doi:10.1093/bioinformatics/btaa558, PMID 32647895
**Survey date:** 2026-08-10

## Sources consulted

| Kind | Reference |
|---|---|
| PDF | `/mnt/c/.../tool_survey/papers/Zulkower2020_dnachisel.pdf` (2 pp., full text extracted with PyMuPDF) |
| Prior analysis | `/mnt/c/.../tool_survey/prior_analyses/Zulkower2020_dnachisel_analysis.md` |
| Repo | https://github.com/Edinburgh-Genome-Foundry/DnaChisel (README.rst; source tree; GitHub API for commits/releases/metadata) |
| Docs | https://edinburgh-genome-foundry.github.io/DnaChisel/ (index, `ref/builtin_specifications.html`, `ref/core_classes.html`, `ref/biotools.html`, `ref/reports.html`, `ref/constraints_reports.html`, `genbank/genbank_api.html`, `examples.html`) |
| Repo source files | `dnachisel/builtin_specifications/*` (directory listing), `EnforceTranslation.py`, `AllowPrimer.py`, `EnforceChoice.py`, `DnaOptimizationProblem/mixins/RecordRepresentationMixin.py`, `Specification/SpecificationSet.py`, `examples/common_scenarios/*` |
| PyPI | https://pypi.org/project/dnachisel/ (JSON API) |
| Web app | https://cuba.genomefoundry.org/sculpt_a_sequence (HTTP 200 as of 2026-08-10) |

## What the tool actually is

From the abstract: *"DNA Chisel is an easy-to-use, easy-to-extend sequence optimization framework allowing to freely define and combine optimization specifications via Python scripts or Genbank annotations."*

From the Usage section: *"An optimization problem is defined in DNA Chisel by a list of global or local specifications against which a starting linear or circular sequence will be optimized. A specification can be either a hard constraint, which must be satisfied in the final sequence, or an optimization objective, whose score must be maximized."*

Solver (Implementation section): *"(i) resolution of all hard constraints, ignoring optimization objectives, and (ii) objectives maximization with respect to the constraints… A region's sequence is optimized via either a stochastic search (using random mutations of the sequence) or an exhaustive search through all possible sequence variants, depending on the number of variants."*

The unit of work is therefore **one sequence in, one optimized sequence out**. The word "variants" in the paper refers to candidate solutions explored *inside* the solver, not to library members emitted to the user. There is no library, pool, or variant-collection object anywhere in the API (`dnachisel/` top level contains only `DnaOptimizationProblem/`, `MutationSpace/`, `SequencePattern/`, `Specification/`, `biotools/`, `builtin_specifications/`, `reports/`, `utils/`, `Location.py`, `cli.py`).

The prior analysis's core claim ("operates on single sequences, not libraries") is **confirmed**. Two prior-analysis details need correction/nuance:
- The prior note says DNA Chisel has "no concept of… combinatorial design". More precisely, it *does* have an internal combinatorial mutation space (`MutationSpace` / `MutationChoice`, and an `EnforceChoice` spec taking a list of allowed alternatives), but the solver collapses it to a single chosen solution; it never enumerates the space as output.
- The prior note omits that DNA Chisel ships primer-oriented functionality (`AllowPrimer` spec, `AvoidHeterodimerization`, `EnforceMeltingTemperature`, and a documented "collection of compatible primers" example). This is de-novo primer *sequence* generation, not mutagenic-primer design for library construction — but a referee could reasonably object if we wrote a flat "no" on primer design.

## Built-in specifications (docs: ref/builtin_specifications.html + source listing)

Pattern: `AvoidPattern`, `AvoidHairpins`, `AvoidStopCodons`, `AvoidRareCodons`, `AvoidChanges`, `AvoidMatches` (Bowtie), `AvoidBlastMatches`, `UniquifyAllKmers`
Sequence: `EnforceSequence` (degenerate allowed), `EnforceTranslation`, `EnforceChanges`, `EnforceChoice`, `EnforcePatternOccurence`, `SequenceLengthBounds`, `EnforceRegionsCompatibility`
Codon: `CodonOptimize` (wrapper), `MaximizeCAI`, `MatchTargetCodonUsage`, `HarmonizeRCA`
Physicochemical: `EnforceGCContent`, `EnforceTerminalGCContent`, `EnforceMeltingTemperature`, `AllowPrimer`, `AvoidHeterodimerization`

Notably **absent** from the specification set and the whole source tree: anything named barcode, variant, mutagenesis, library, HGVS, VCF, GFF/GTF.

## Documented examples / vignettes (docs `examples.html`, repo `examples/common_scenarios/`)

1. **Optimization with report** (`optimization_with_report.py`) — load a Genbank record, optimize, emit a multi-file report.
2. **Optimization of a circular sequence** (`circular_sequence.py`) — cross-origin specs on a circular molecule.
3. **Sequence standardization / domestication** (`gene_domestication.py`, `sequence_without_restriction_sites.py`) — remove restriction sites while codon-optimizing and holding GC.
4. **Programmatic optimization of a full genome** (`genome-wide-optimization.py`) — downloads the E. coli genome via `genome_collector`, loops over CDS features in Python, builds one `DnaOptimizationProblem` per gene, writes a multi-record FASTA. *This is the closest thing to a multi-sequence workflow and it is an explicit user-written `for` loop.*
5. **Creating a collection of compatible primers** (`primers_collection.py`) — iteratively builds 20 mutually non-heterodimerizing 20-mers with ~60% GC; each new primer is constrained against all previously accepted ones. Docs note DNA Chisel is "not originally designed for sequence collection creation" but can generate "inter-compatible sequences".
6. **Writing a custom Specification: 9-mer score** (`gen9_9mer_score_minimization.py`) — user-defined objective class.
7. **Sequence without E. coli TF binding sites** (`sequence_without_tf_binding_sites.py`) — downloads TFBS data, avoids all motifs.
8. **Constraints reports example** (`constraints_breaches_report/`) — `constraints_breaches_dataframe(constraints, sequences)` over a *list/dict of many sequences*, exported to Excel + `breaches_records_to_pdf`.
9. **Plasmid optimization** (`plasmid_optimization.py`), **E. coli genes optimization** (`ecoli_genes_optimization.py`).
10. **Manuscript examples** (`examples/manuscript_examples/`): `B_multiobjectives/`, `C_sequence_without_sites.py`, `D_genome-wide-optimization.py` — correspond to Supplementary Sections S1B–S1D.

**Reproducibility target for PoolParty:** #5 (compatible primer/oligo collection under mutual-heterodimer constraints) and #8 (batch synthesis-constraint QC over a set of sequences) are the two examples that touch multi-sequence territory and are the most interesting to attempt in PoolParty. #4 (per-gene domestication across a genome) is a natural "PoolParty expresses this as one library object instead of a for-loop" demonstration.

## Capability-by-capability evidence

### Block A — library specification

- **library_as_object — no.** No class in the API represents a set of designed sequences. `DnaOptimizationProblem(sequence=…)` takes a single sequence; docs core-classes page lists only `DnaOptimizationProblem`, `CircularDnaOptimizationProblem`, `Specification`, `SpecEvaluation`, `Location`, `MutationSpace`/`MutationChoice`, pattern classes. Multi-sequence work is an explicit Python loop, e.g. `genome-wide-optimization.py`: `for feature in genome_record.features: … problem = DnaOptimizationProblem(...); optimized_records.append(problem.to_record(...))`. The only multi-sequence API is a *reporting* utility, `constraints_breaches_dataframe(constraints, sequences)`, which validates sequences the user already has.
- **dag_chaining — partial.** Specifications compose freely and can nest one level via `SpecificationSet` ("Generic class for writing Specs which are actually made of more specs… the initialization actually creates a dictionary of standard Specifications in the DNAOptimizationProblem"; `AllowPrimer` is such a set). But this composes *constraints on one sequence*, not design steps producing sequences. There is no pipeline/graph object; chaining problems (output of one feeding the next) is done by hand in user Python, as in `primers_collection.py`.
- **lazy_evaluation — no.** The problem holds a concrete `sequence` string in memory and mutates it in place (`_replace_sequence`, `sequence_edits_as_array`). Solver enumerates candidate variants internally but nothing is emitted lazily to the user; there are no generators/iterators over designed sequences.
- **mixed_mutagenesis_one_pool — no.** There is no mutagenesis abstraction at all (no exhaustive-singles/pairs/random-sampling/WT-replicate concepts). `EnforceChanges` forces the sequence to differ from the original in a region, and `EnforceChoice` restricts a region to a list of alternatives, but the solver returns one sequence satisfying them, not a pool of different mutagenesis regimes.
- **combinatorial_motif_place — no.** `EnforcePatternOccurence` sets a required *count* of a pattern and the solver may `insert(...)` it ("an `insert(CGTCTC)` constraint will attempt to place the pattern 'CGTCTC' at different locations of the sequence"), but the alternative placements are search states, not enumerated library members. No API for multiple motifs × positions × orientations × permutations as a product.
- **barcode_generation — no.** No barcode/UMI functionality: the string "barcode" appears nowhere in the built-in specification list, the biotools reference, the docs TOC, or the `builtin_specifications/` file listing. Closest analogue is the primer-collection example, which builds mutually compatible short sequences under Tm/GC/heterodimer constraints but has no edit-distance/decoding notion and no mechanism to attach the results to library members.
- **per_sequence_provenance — partial.** For the single output sequence, `problem.to_record(..., with_constraints=True, with_objectives=True, with_sequence_edits=True)` writes GenBank features for every specification applied and for each nucleotide edit ("original=>modified"), and `optimize_with_report()` emits a PDF/plots/GenBank bundle described in the paper as "comprehensive optimization reports for traceability and troubleshooting". This is provenance of *how one sequence was edited*, not structured metadata about a library member's construction lineage.
- **automatic_naming — no.** `to_record(..., record_id=None)`: the id is whatever the user passes (in the genome example the user supplies `protein_id`), otherwise inherited from the input record. No generated names encoding what was designed.
- **design_visualization — partial.** Reports render the sequence with annotated constraints, objectives and edits as PDF plots (paper Fig. 1C; `write_optimization_report` produces "a PDF summary, plots, and genbanks"; `max_features_in_plots` option), plus `breaches_records_to_pdf` for constraint breaches across many records, plus a troubleshooting figure for unresolvable regions. But there is no view of a design graph/pipeline and no library-level view — because there is no library.

### Block B — assay coverage

- **assay_dms — no.** No DMS/variant-library design capability; no codon-substitution enumeration, no per-position saturation. (Third parties, e.g. Tycko et al. 2020, use DNA Chisel to codon-optimize oligos *after* designing a DMS library elsewhere — that is post-processing, not library design.)
- **assay_mpra — no.** No regulatory/MPRA library design. Motif-related functionality is *avoidance* (`sequence_without_tf_binding_sites.py`, `MotifPssmPattern` used to forbid matches), never combinatorial motif installation across a library.
- **assay_insilico — no.** Nothing about model-probing libraries; the tool predates and does not mention genomic AI models. Objectives are biochemical/manufacturing scores only.

### Block C — genomics integration

- **genome_coordinates — no.** `Location(start, end, strand)` is defined relative to the problem's own sequence: docs state "Location(5, 10) represents sequence[5, 6, 7, 8, 9]". A genome can be *loaded* as a Biopython record and features extracted by the user (`feature.location.extract(genome_record)` in the genome-wide example), but there is no reference-assembly coordinate system, no chrom:pos input, no coordinate liftover back to the genome.
- **transcript_models — no.** Input is FASTA/GenBank/Snapgene (`load_record`). The GenBank parser reads *specification annotations* (`@`/`~`-prefixed labels), not gene models. No GFF/GTF support anywhere in docs or source; no transcript/isoform concept.
- **exon_intron_split_codons — no.** `EnforceTranslation.set_location` validates that "Location length in Codon Specifications should be a 3x" and iterates codons contiguously via `codon_index_to_location(i)`. Location is a single contiguous (start, end, strand); no CompoundLocation/exon-array handling and no logic for codons spanning an intron.
- **hgvs_input — no.** No HGVS parsing anywhere in the API, docs, or annotation shorthand vocabulary (`@no()`, `@insert()`, `@keep`, `@change`, `@cds`, `@sequence`, `@choice`, `@gc`, `~use_best_codon`, …).
- **vcf_vep_output — no.** Outputs are GenBank/FASTA records, PDF reports, and Excel/dataframe constraint tables (`write_record` supports "genbank or fasta"). No VCF writer.
- **consequence_annotation — no.** `annotate_differences` / `sequences_differences` mark where sequences differ, and `EnforceTranslation` guarantees changes are synonymous, but nothing classifies molecular consequence (stop-gained, missense, frameshift, in-frame indel). Consequence is *constrained away*, never annotated.

### Block D — physical construction

- **primer_design — partial.** `AllowPrimer` ("Enforce various specifications for enabling primers at the location… useful for making sure that you will be able to conduct a PCR or sanger sequencing with a primer annealing at a particular location") with `tmin/tmax`, `max_homology_length`, `avoid_heterodim_with`, `max_heterodim_tm`, `avoided_repeats`; plus `AvoidHeterodimerization` and `EnforceMeltingTemperature`; plus the documented 20-primer compatible-collection example. **However**: these design/validate primer *annealing sites and sequences* de novo; there is no mutagenic-primer design (no NNK/degenerate mutagenesis primers, no assembly-oligo tiling, no protocol-aware oligo splitting for a variant library).
- **codon_optimization — yes.** Core feature. `CodonOptimize` wrapper plus `MaximizeCAI`, `MatchTargetCodonUsage`, `HarmonizeRCA`, `AvoidRareCodons`, `EnforceTranslation`; codon tables from the Codon Usage Database (Nakamura 2000). Paper: "Gene NUM1 of Saccharomyces cerevisiae will be codon optimized for Escherichia coli".
- **synthesis_constraints — yes.** This is the tool's raison d'être. README: "meeting commercial DNA provider requirements"; specs `AvoidPattern` (restriction sites, homopolymers, repeated k-mers, regex), `AvoidHairpins` ("per IDT guidelines"), `EnforceGCContent` (global/windowed), `UniquifyAllKmers`, `AvoidMatches`/`AvoidBlastMatches`, `SequenceLengthBounds`, `EnforceTerminalGCContent`, plus batch QC via `constraints_breaches_dataframe`. Applies to individual sequences (checked/enforced one sequence at a time).

### Block E — engineering

- **interface — API + CLI + web app.** Paper: "It can be used as Python library, a command-line interface, or a web application". CLI: `dnachisel annotated_record.gb optimized_record.gb` (`dnachisel/cli.py`; docs "Command Line Interface" section). Web app: https://cuba.genomefoundry.org/sculpt_a_sequence (part of EGF CUBA), with GenBank drag-and-drop and an annotation editor (Fig. 1B). No desktop GUI.
- **license — MIT.** "MIT license, Copyright 2017 Edinburgh Genome Foundry, University of Edinburgh" (README.rst); PyPI metadata `License: MIT`.
- **maintained — yes, actively.** Latest PyPI release **v3.2.16, 2025-05-10**; latest commit on `master` **2025-05-10** ("Include pkg data, exclude unwanted folders", author veghp); GitHub `pushed_at` 2025-11-05 (branch activity), repo not archived, 278 stars, 22 open issues, created 2017-09-07. Release cadence through 2024–2025 (v3.2.12 Nov 2024, v3.2.13 Jan 2025, v3.2.14/15 Apr 2025, v3.2.16 May 2025). Maintenance has passed from Zulkower to veghp (Peter Vegh, EGF).

## Availability today (2026-08-10)

Installable and runnable. `pip install dnachisel` / `pip install 'dnachisel[reports]'`; latest release 3.2.16 (2025-05-10) on PyPI. Repo is public, MIT, not archived, last master commit 2025-05-10, repo pushed_at 2025-11-05. Docs site https://edinburgh-genome-foundry.github.io/DnaChisel/ returns HTTP 200. Web app https://cuba.genomefoundry.org/sculpt_a_sequence returns HTTP 200 (a Vue SPA, so a plain fetch shows only the "EGF CUBA" shell — server is up; interactive optimization was not exercised because this round is document research only). Optional external dependencies: Bowtie (`AvoidMatches`), BLAST (`AvoidBlastMatches`), Matplotlib/pdf-reports/sequenticon (reports extra).

## Stated limitations / scope statements

- The problem is defined on "a starting linear or circular sequence" (singular) — the framework has no library-level scope.
- Docs for the primer-collection example concede DNA Chisel is not designed for collection creation, only that constraint-based iteration can produce "inter-compatible sequences".
- Solver is heuristic (stochastic or exhaustive local search); constraints may be unresolvable, in which case the tool emits "a troubleshooting figure representing the problematic region" rather than a solution.
- `EnforceTranslation` requires a contiguous, in-frame, length-divisible-by-3 location (no spliced CDS).
- GenBank annotation interface does not support all specs — docs note `AvoidHeterodimerization` and `EnforceRegionsCompatibility` cannot be declared as annotations.
- Report plotting degrades with many features ("plots with thousands of features may take ages to plot").

## Capabilities worth considering as new matrix rows

1. **Constraint/objective solving on an existing sequence** (repair/optimize rather than generate) — DNA Chisel's defining capability, orthogonal to every current row.
2. **Batch synthesis-constraint QC over a set of sequences** — `constraints_breaches_dataframe(constraints, sequences)` accepting `[("name","seq")…]` / dict / records, exporting Excel + PDF. This is a genuine multi-sequence feature and the fairest thing to credit DNA Chisel with in a library-oriented table.
3. **Extensibility by user-defined design primitives** — new `Specification` subclasses with custom `evaluate`/`localized` methods, optionally registered with the GenBank parser. Directly comparable to PoolParty's user-defined Operations.
4. **Annotation-driven (no-code) design specification** — problems declared as GenBank feature labels (`@keep`, `@no(BsaI_site)`, `~use_best_codon`) editable in SnapGene/Benchling.
5. **Circular-sequence / cross-origin support** — `CircularDnaOptimizationProblem`.
6. **Homology / repeat-uniqueness control across a construct** — `UniquifyAllKmers`, `AvoidMatches` (Bowtie), `AvoidBlastMatches`.
7. **Oligo/primer cross-reactivity control** — `AvoidHeterodimerization`, `max_heterodim_tm`; relevant if PoolParty claims anything about pool-wide oligo compatibility.

## Confidence notes

- Highest confidence: absence of library/mutagenesis/barcode/HGVS/VCF/GFF machinery (checked docs TOC, full built-in-spec list, and the actual `builtin_specifications/` and `dnachisel/` file listings — nothing exists under those names).
- `dag_chaining = partial` is the judgement call most open to challenge. Specifications compose and nest (`SpecificationSet`), which the authors legitimately call "freely define and combine"; what is absent is a graph of *design steps generating sequences*. Evidence is worded to make that distinction explicit.
- `primer_design = partial` is the second judgement call. `AllowPrimer` and the primer-collection example are real, so a flat "no" would be attackable; but no mutagenic/assembly oligo design exists.
- `per_sequence_provenance = partial` and `design_visualization = partial` both hinge on reading "per sequence" as "the one output sequence". If the matrix intends library-member provenance/visualization, both are effectively "no".
- Not verified by execution (document research only): actual behaviour of the CUBA web app, and whether any undocumented helper handles compound locations. The `EnforceTranslation` source read makes spliced-CDS support very unlikely.
- Supplementary Sections S1A–S1E of the paper were not obtainable (supplementary PDF not in the local papers directory); their content is inferred from the main-text descriptions and the matching `examples/manuscript_examples/` scripts in the repo.
