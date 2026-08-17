# Oligopool Calculator — evidence memo

**Slug:** `oligopoolcalc`
**Tool:** Oligopool Calculator (`oligopool` Python package)
**Authors / lab:** Ayaan Hossain, Daniel P. Cetnar, Travis L. LaFleur, James R. McLellan, Howard M. Salis (Salis lab, Penn State)
**Paper:** Hossain et al., *ACS Synth. Biol.* 2024, 13(12), 4218–4232. DOI 10.1021/acssynbio.4c00661 (PMID 39641628, PMCID PMC11669329)
**Survey date:** 2026-08-10

---

## 1. Sources consulted

| Kind | Reference |
|---|---|
| PDF | `revision/tool_survey/papers/Hossain2024_OligopoolCalculator.pdf` (15 pp., text extracted with PyMuPDF) |
| Prior analysis | `revision/tool_survey/prior_analyses/Hossain2024_OligopoolCalculator_analysis.md` |
| Repo | https://github.com/ayaanhossain/oligopool (canonical; also mirrored at https://github.com/hsalis/SalisLabCode) |
| Docs | https://github.com/ayaanhossain/oligopool/blob/master/docs/docs.md (user guide, ~99 KB) |
| Docs | https://github.com/ayaanhossain/oligopool/blob/master/docs/api.md (API reference, ~76 KB) |
| Docs | https://github.com/ayaanhossain/oligopool/blob/master/docs/agent-skills.md (LLM-agent guide) |
| Repo | https://github.com/ayaanhossain/oligopool/blob/master/README.md |
| Repo | https://github.com/ayaanhossain/oligopool/tree/master/examples (+ `examples/README.md`, `examples/design-assembly-parser/README.md`, `examples/OligopoolCalculatorInAction.ipynb`) |
| PyPI | https://pypi.org/project/oligopool/ (JSON API) |
| Web | GitHub commits API for `ayaanhossain/oligopool` |

**Important:** the shipped software has moved substantially beyond the 2024 paper. The paper describes Design Mode (barcode/primer/spacer/motif/split/pad) + Analysis Mode. The current release (v2026.02.22.1) adds a **CLI with YAML pipelines including explicit DAG execution**, **Degenerate Mode (`compress`/`expand`)**, **QC Mode (`lenstat`/`verify`/`inspect`)**, **Patch Mode**, `join`/`merge`/`revcomp`, Docker, and an LLM-agent guide. Any claim in the referee response should be checked against the repo, not only the paper.

---

## 2. Verdict on the prior analysis

The prior analysis is **broadly correct** on the core framing ("designs the physical infrastructure of oligo pools; assumes variant sequences are already designed as input"). Two corrections:

1. It says "no DAG architecture". **That is now wrong as stated.** The CLI YAML pipeline explicitly supports DAG step graphs with `after:` dependencies and level-based parallel execution, and the Python API has a `join` module for recombining parallel design branches. The DAG is an *execution/scheduling* graph over design steps, not a *content-composition* graph over library sub-populations — that distinction is the one to make in the referee response, not "no DAG".
2. It underspecifies the current feature set (missing Degenerate Mode, Patch Mode, QC Mode, CLI/YAML, declarative design parser example).

The "complementary, not competing" conclusion holds and is well supported by the docs.

---

## 3. Module inventory (v2026.02.22.1)

| Mode | Modules |
|---|---|
| Design | `barcode`, `primer`, `motif`, `spacer`, `background`, `merge`, `revcomp`, `join`, `final` |
| Assembly | `split`, `pad` |
| Degenerate | `compress`, `expand` |
| Analysis | `index`, `pack`, `acount`, `xcount` |
| QC | `lenstat`, `verify`, `inspect` |
| Advanced | `vectorDB` (k-mer background store), `Scry` (barcode classifier) |

---

## 4. Capability assessment with evidence

### Block A — library specification

**`library_as_object` = partial.**
The library is a plain **pandas DataFrame / CSV**, not a dedicated library class. docs.md "Core Concepts → The DataFrame Flow":
> "**Input**: CSV path or pandas DataFrame with an `ID` column … **Chainable**: Output of one module feeds into the next"

Every module operates on the *whole table at once*, so the user does not loop over sequences — that part is genuinely first-class. But there is no library object carrying design semantics; the "library" is the accumulated set of columns, and `final` "Concatenate[s] all columns into synthesis-ready oligos" (api.md), after which non-sequence columns are dropped (docs.md advises: *"Keep a saved annotated design DataFrame/CSV before `final` if you plan to `index`, `pack`, and count later."*). Table-as-library, not library-as-object.

**`dag_chaining` = yes.**
Two mechanisms.
(a) CLI YAML pipelines, docs.md §"Parallel Pipeline Execution":
> "**Parallel/DAG format**: a list of step objects with `name`, `command`, and optional `after`, which defines dependencies explicitly."
> "Steps with no `after` dependencies form level 1 (eligible to run concurrently). Each subsequent level waits for dependencies, and independent steps are grouped automatically for concurrent execution."
(b) Python-level branch recombination via `join`, docs.md:
> "Most Design Mode workflows are sequential, because element constraints tend to couple steps. If two design steps are truly independent … you can run them as parallel branches and then recombine their outputs with `join`."
Caveat for the response: docs.md admits *"design workflows are often dependency-dense and mostly sequential"* and calls design branch-join "**Rare** parallel design branches". The DAG schedules *design steps applied to one flat table*; it does not compose heterogeneous sub-libraries into one pool.

**`lazy_evaluation` = no.**
Everything is materialized in DataFrames/CSVs. The one enumerating module, `expand`, is explicitly eager with a safety cap: docs.md — *"Expansion can be exponential (N positions = 4^N for all-N sequence) … Use `expansion_limit` as a safety cap for highly degenerate sequences."* No generator/iterator/on-demand API anywhere in api.md.

**`mixed_mutagenesis_one_pool` = no.**
The tool does not generate mutagenized variants at all, so mixing mutagenesis types in one spec is out of scope. docs.md Quick Start: *"# Start with your variants"*. docs.md §Workflows: *"Start with your variants in a CSV (must have an `ID` column)"*. agent-skills.md recipe table: *"Cost-efficient saturation mutagenesis | generate all substitutions → `compress` → order `synthesis_df`"* — the "generate all substitutions" step is the user's, external to the tool. The notebook is explicit that variant content is out of scope: *"let us first generate a random library of 59-73 mers (**designing exact ribozymes is beyond the scope here**)"*.

**`combinatorial_motif_place` = no.**
`motif` inserts **one** motif into **one** named column position per call — either per-variant (`motif_type=0`, a distinct sequence per row satisfying an IUPAC constraint) or a single constant anchor (`motif_type=1`). api.md: *"Insert sequence motifs (per-variant or constant anchors) with constraint satisfaction"*. Multiple motifs are placed by chaining multiple `motif` calls at fixed column positions; there is no enumeration over positions, orientations, or motif combinations, and no cross-product of motif sets. Paper §"Restriction Sites, Motifs, and Spacers": motifs are *"defined using the IUPAC degenerate nucleotide code"* and the pathfinder *"identifies an oligonucleotide-specific sequence solution that satisfies design rule #1"* — one solution per oligo, not a combinatorial expansion. (`revcomp` flips a **range of columns** globally for the whole table, not per-variant orientation enumeration.) The nearest thing to combinatorial content generation is `expand`, which enumerates all concrete sequences from an IUPAC-degenerate string — but it is documented as a *verification* step for `compress`, not a design primitive.

**`barcode_generation` = yes.**
`op.barcode(input_data, oligo_length_limit, barcode_length, minimum_hamming_distance, maximum_repeat_length, barcode_column, barcode_type, left/right_context_column, patch_mode, cross_barcode_columns, minimum_cross_distance, excluded_motifs, background_directory, random_seed)`. Constraints: minimum pairwise **Hamming distance**, max shared repeat length (kills homopolymers), excluded motifs (restriction sites; list/CSV/FASTA/DataFrame/dict of sources), context/edge-effect awareness, background k-mer screening, cross-set orthogonality to previously designed barcode columns. The barcode is inserted as a column at the correct position in the oligo, i.e. genuinely attached. Paper: *"4 million highly unique and compact barcodes in 1.2 h"*; barcode design is its flagship feature ("orthogonally symmetric barcode design").
**One honest gap:** there is **no GC-content parameter** for barcodes — the row's "(GC, edit distance)" is only half satisfied. GC is controlled only indirectly (excluded motifs, repeat length) and directly only for primers via the IUPAC sequence constraint (api.md: *"`'SS' + 'N'*18` for GC clamp"*) and Tm windows.

**`per_sequence_provenance` = partial.**
Each designed element lives in its own named column, so a row structurally records *which* barcode/primer/motif/spacer composes that oligo, and `verify` emits JSON-serialized diagnostic columns (`IntegrityConflictDetails`, `LengthConflictDetails`, `ExmotifConflictDetails`, `BackgroundConflictDetails`). But there is no record of *how a variant was derived* (no mutation/design history), because the tool never derives variants. The per-run `stats` dict (`status`, `basis`, `step`, `step_name`, `vars`, `warns`, `module`, `input_rows`, `output_rows`) is **aggregate, not per-sequence** — docs.md: *"`stats` is **aggregate**: programmatic pass/fail decisions and summary reporting."*

**`automatic_naming` = no.**
docs.md: input *"must have an `ID` column"*; the ID is user-supplied and required, and modules never invent or enumerate sequence names. Output *column* names are user-specified (`barcode_column='BC1'`); output *file* names get auto-suffixes (`.oligopool.<module>.csv`). `compress` mints `DegenerateID` labels for degenerate groups, which is the only auto-generated identifier found, and it names groups, not designed variants.

**`design_visualization` = no.**
No plotting, no graph view, no sequence-highlighting anywhere in README/docs/api. Grep for `visuali|plot|matplotlib|graph` across README.md, docs.md, api.md, agent-skills.md returns no feature hits. The tutorial notebook imports `matplotlib`/`seaborn`, but those are the *author's* ad-hoc plots of length histograms and count distributions, not tool functionality. Inspection is textual: `lenstat` (length ruler), `verify` (conflict table), `inspect` (binary-artifact metadata), CLI `--dry-run` (prints DAG execution levels as text). Architecture figures in the docs are hand-drawn SVGs.

### Block B — assay coverage

**`assay_dms` = partial.**
No DMS content generation, no codon or amino-acid awareness anywhere in the codebase docs. But DMS pools are an explicitly documented *input*: docs.md §"Saturation Mutagenesis Compression" — *"**The ask**: 'I have 1000 single-amino-acid substitution variants. Can I synthesize fewer oligos?' … **The decomposition**: Generate all substitution variants, then `compress` into IUPAC-degenerate oligos."* agent-skills.md lists both *"Saturation mutagenesis (barcoded) | generate all substitutions → standard design pipeline for individual tracking"* and *"Saturation mutagenesis (degenerate)"*. So: it will build the pool infrastructure around a DMS library you enumerated elsewhere, and compress it 6–20×; it will not enumerate the substitutions.

**`assay_mpra` = yes.**
This is the tool's home turf. Paper title/keywords: *"massively parallel reporter assay"*; three real MPRA libraries built and characterized (93,180 promoters; 62,120 5′UTR mRNA-stability elements; 6,232 ribozymes). docs.md ships a *"Promoter MPRA Library"* application template (`background → primer(fwd) → motif(anchor) → motif(anchor) → barcode → primer(rev,paired) → spacer → lenstat → verify → final`) and notes *"This same single-barcode MPRA pattern also applies to other regulatory/stability-element libraries (e.g., UTR variants, degrons, structured elements)."* Plus Analysis Mode counts the MPRA reads. **Caveat to state in the response:** it designs the synthesis-ready MPRA oligo (primers, anchors, barcodes, spacers, splits/pads) and analyzes the data, but the regulatory variants themselves are user input.

**`assay_insilico` = no.**
No sequence-to-function model, no in-silico perturbation, no genomic AI integration. The paper mentions ML only downstream: count matrices *"can be combined with biophysics and/or machine learning to develop predictive sequence-to-function models"* (i.e. the user builds a model from the wet-lab data afterward). The README's "AI workflows" refers to **LLM agent skills for driving the tool** (`docs/agent-skills.md`), not to probing genomic AI models. Degenerate Mode mentions *"ML-generated libraries … often compress well"* — again, models are upstream and external.

### Block C — genomics integration (all no)

Grep across README.md, docs.md, api.md, agent-skills.md for `GTF|GFF|HGVS|VCF|genome coord|chrom|codon` returns **zero** hits. The entire I/O contract is a CSV/DataFrame of `ID` + DNA-string columns.

- **`genome_coordinates` = no.** No coordinate input anywhere. The only genome-shaped input is `background(input_data=<FASTA/sequence>)`, which builds a k-mer database purely for off-target screening — paper: *"e.g., the genome sequence of E. coli, which can appear as"* background. Coordinates are never parsed or emitted.
- **`transcript_models` = no.** No GTF/GFF anywhere.
- **`exon_intron_split_codons` = no.** No reading frame, exon, intron, or codon concept in any module signature.
- **`hgvs_input` = no.** No HGVS parsing; variants arrive as raw DNA strings.
- **`vcf_vep_output` = no.** Outputs are CSV (design tables, count matrices) and binary `.oligopool.index` / `.oligopool.pack` artifacts. No VCF writer.
- **`consequence_annotation` = no.** No molecular-consequence vocabulary anywhere; the tool has no protein-level model of its sequences.

### Block D — physical construction

**`primer_design` = yes** (with a scope caveat).
`op.primer(...)` designs PCR/universal primer binding sites with `primer_sequence_constraint` (IUPAC, e.g. GC clamp), `minimum/maximum_melting_temperature`, `maximum_repeat_length`, `paired_primer_column` for Tm matching, `oligo_set` for per-sub-pool primers, `excluded_motifs`, and `background_directory`. Paper: primers are designed for *"maximum specificity and efficiency during PCR, which includes precise matching of primer melting temperatures, eliminating all off-target binding, eliminating inhibitory DNA structures, and greatly reducing primer dimer formation"*; docs list screening for *"hairpin, homodimer, heterodimer, cross-primer dimer"*. Benchmark: *"the design of universal primer binding sites for one million 200-mer oligos in 15 min"*. `split` designs overlap regions (Tm- and Hamming-constrained) and `pad` adds Type IIS pads for scarless Golden Gate assembly.
**Caveat:** these are **amplification/assembly** primers, not **mutagenic** primers. There is no site-directed-mutagenesis primer mode (contrast Mutation Maker). If the row is read strictly as "mutagenic primers", the answer is no; as "designs oligos/primers for wet-lab protocols", clearly yes.

**`codon_optimization` = no.**
Zero hits for "codon" in the entire documentation set. No translation table, no CDS handling, no reading frame.

**`synthesis_constraints` = yes.**
Enforced per sequence: `oligo_length_limit` on every module; `maximum_repeat_length` (bans homopolymers and shared repeats, including junction/edge effects from flanking context); `excluded_motifs` (restriction sites and their reverse complements, palindromes, polymeric runs); `background_directory` k-mer screening against host genome/plasmid; Tm/structure for primers and overlaps. `verify` is a dedicated pre-synthesis QC pass emitting per-row `HasIntegrityConflict`, `HasLengthConflict`, `HasExmotifConflict`, `HasBackgroundConflict`, `HasAnyConflicts` plus JSON detail columns; `lenstat` reports *"length statistics and free space remaining (non-destructive)"*. Over-length constructs are handled constructively via `split` → `pad` rather than merely flagged; cost is handled via `compress`. Notably, **motif emergence at junctions** is checked (*"Motif emergence = count exceeds library-wide minimum"*), which most competing tools do not do.

### Block E — engineering

**`interface` = Python API + CLI + YAML pipelines + Docker + Jupyter; no web server, no GUI.**
README: *"seamlessly integrating with Python, CLI, Jupyter, containers, and AI workflows"*; installs two console entry points, `oligopool` and `op`. CLI supports YAML config for single commands and multi-step pipelines (serial, DAG, or mixed) with `--dry-run` and CLI flag overrides. Docker image documented in `docs/docker-notes.md`. Searched for a hosted web server: none exists (the paper's Data Availability lists only the two GitHub repos and PyPI; `oligopool.com` is an unrelated commercial site, not this tool).

**`license` = GPL-3.0-only.**
PyPI metadata: `GPL-3.0-only`. Paper: *"available at https://github.com/hsalis/SalisLabCode and https://github.com/ayaanhossain/oligopool under a GPLv3 open-source license."*

**`maintained` = last commit 2026-02-22; last PyPI release 2026-02-22 (v2026.2.22.1).**
GitHub commits API: most recent commits 2026-02-22 02:46 UTC (`chore: rerun notebook`), 02:43 (`chore: bump version to v2026.02.22.1`), 02:40 (docs update v2026.02.18 → v2026.02.22). Actively maintained ~6 months before this survey, and under rapid iteration (two version bumps within four days in Feb 2026).

---

## 5. Availability today

**Installable and runnable today: yes.**
- PyPI package `oligopool`, latest **2026.2.22.1** released **2026-02-22**; `requires-python >=3.10,<4`. First release 2024-12-02 (v2024.12.2), contemporaneous with the paper.
- `pip install --upgrade oligopool`, or `pip install git+https://github.com/ayaanhossain/oligopool.git`, or Docker.
- Repo https://github.com/ayaanhossain/oligopool alive, last commit 2026-02-22; mirror at https://github.com/hsalis/SalisLabCode.
- Docs, API reference, agent guide, four runnable example projects and a re-executed tutorial notebook all live in-repo.
- Platform note from README: Linux, macOS, and **Windows Subsystem for Linux** (not native Windows).
- No hosted web service exists — and none is claimed.

---

## 6. The tool's own documented examples (candidates for PoolParty reproduction)

**From the paper (real, published libraries):**
1. **93,180 designed promoter variants** (Fig. 2A/B) — architecture: universal primer binding site, BamHI + HindIII flanks, promoter variant, unique barcode (distinguishable up to 6 mismatches), EcoRI + SpeI sites for a 2-step cloning procedure that shifts the barcode into the 3′UTR; read out as mRFP1 transcription rate. Used to train the Promoter Calculator / Multi-Sigma Promoter Calculator.
2. **62,120 5′UTR mRNA-stability elements** (Fig. 2C/D) — same architecture family with longer primer binding sites and shorter barcodes, different cloning enzymes; sfGFP readout of mRNA decay.
3. **6,232 hammerhead ribozymes** (Fig. 2E/F) — double barcode (BC1 uncleaved-proportional, BC2 total-mRNA-proportional) + triple primer binding sites; ribozymes themselves designed by the Non-Repetitive Parts Calculator to ≤16-nt max shared repeat. Self-cleavage efficiency by combinatorial counting; index hopping detected and removed.

**From the repo (runnable):**
4. **`examples/OligopoolCalculatorInAction.ipynb`** — full walkthrough (284 cells, updated 2026-02-18) building a 6,232-variant ribozyme oligopool: variants simulated as random 59–73-mers screened against EcoRI/AatII/BsaI; 5 kb random plasmid stored as background; three 20-bp primers; two 11-bp barcodes (Hmin 3, cross-set separated); 3–17-bp spacers padding to 170 bp; then Assembly Mode (`split`/`pad`), Degenerate Mode (`compress`/`expand`, incl. "Realistic Example: 500 ML-style variants"), and Analysis Mode (`index`/`pack`/`acount`/`xcount` with a synchronized callback).
5. **`examples/design-assembly-parser/`** — declarative spec-based parser over **4,351 real promoter sequences** (`promoters.txt`), architecture `[Primer1]-[Cut1]-[Promoter]-[Barcode]-[Primer2]-[Cut2]-[Primer3]-[Filler]`, element types `variant | primer | barcode | motif | spacer`. **This is the closest thing in the ecosystem to PoolParty's declarative library spec and is the highest-value reproduction target.**
6. **`examples/library-compressor/`** — Degenerate Mode: compress a variant library into IUPAC oligos.
7. **`examples/analysis-pipeline/`** — NGS indexing, packing, combinatorial counting.
8. **`examples/cli-yaml-pipeline/`** — `analysis_single.yaml`, `analysis_multi.yaml` (multi-sample DAG: index once → per-sample pack → parallel `xcount`/`acount`).

**Application templates in docs.md §"Application Templates":**
9. **Promoter MPRA Library** — `background(host.fasta) → primer(fwd) → motif(const anchor) → motif(const anchor) → barcode → primer(rev,paired) → spacer → lenstat → verify → final`.
10. **CRISPR Guide Library** — `background → primer(fwd) → barcode → spacer → lenstat → verify → final`; *"Guides are the variant column. Add `excluded_motifs` for cloning-system recognition sites (e.g., BsmBI: `CGTCTC`/`GAGACG`, BbsI: `GAAGAC`/`GTCTTC`)."*
11. **Ribozyme Library (Pre/Post Cleavage Readout)** — dual barcode with `cross_barcode_columns=['BC_molecule']`.
12. **Saturation Mutagenesis Compression** — 1,000 single-amino-acid substitution variants → `compress` → order `synthesis_df` (6–20× fewer oligos).

---

## 7. Notable capabilities NOT covered by the current row list

These are real strengths of the tool that the matrix would otherwise miss; several may deserve new rows (or an explicit "out of scope for this comparison" note in the response).

1. **NGS analysis / demultiplexing half (Analysis Mode).** `index`, `pack`, `acount` (association counting), `xcount` (combinatorial multi-barcode counting), the **Scry** barcode classifier, efficient read packing; ~500M reads/hour on 8 cores; index-hopping detection and mismatch-distribution reporting. No other tool in this survey closes the design→measurement loop.
2. **Degenerate-oligo compression (`compress`/`expand`).** Lossless IUPAC compression of a concrete variant library into far fewer synthesis oligos (6–20× documented), plus a `mapping_df` that maps sequenced survivors back to original variant IDs. This is a *synthesis-cost* capability no other surveyed tool has.
3. **Background off-target screening against a host genome/plasmid** (`background` / `vectorDB` k-mer DB), applied junction-aware to barcodes, primers, motifs, spacers, and re-checked in `verify`. Multiple background DBs can be layered.
4. **Length-overflow resolution by construction**: `split` (Tm- and Hamming-constrained overlaps for Gibson/homology assembly) and `pad` (Type IIS pads for scarless Golden Gate).
5. **Patch Mode** — incrementally extend an already-ordered library (append rows, fill only missing elements) without redesigning or perturbing existing barcodes/primers. Directly relevant to real multi-round library projects.
6. **Cross-set barcode orthogonality** (`cross_barcode_columns`, `minimum_cross_distance`) and **tandem/multi-barcode architectures** (CombiGEM-style), designed *and* counted combinatorially.
7. **Per-sub-pool primer design** (`oligo_set`) with cross-set dimer screening and per-set Tm matching — multiplexed sub-libraries amplifiable independently from one synthesis pool.
8. **Declarative architecture spec** (`examples/design-assembly-parser`): an element-list + element-spec dictionary rather than imperative calls. The nearest analogue to PoolParty's declarative style, though it is a worked example, not a shipped API.
9. **CLI + YAML config with DAG scheduling, `--dry-run`, YAML merge-key parameter sharing, and config precedence.**
10. **Scale.** 4M barcodes in 1.2 h; universal primer binding sites for 1M 200-mers in 15 min; ~35 min for a 1M-variant oligopool — on an 8-core desktop. If the referee cares about scale, this is the benchmark to beat or concede.
11. **LLM-agent-ready documentation** (`docs/agent-skills.md`) — a design decision worth noting given current trends.

---

## 8. Stated limitations / scope boundaries (author's own words)

- Variant content is the user's responsibility: *"Start with your variants in a CSV (must have an `ID` column)"*; notebook: *"designing exact ribozymes is beyond the scope here."*
- Excluded motifs are **global**: docs.md — *"Excluded motifs are applied globally to all variants (no per-variant exclusion)."*
- Design pipelines are mostly sequential: *"design workflows are often dependency-dense and mostly sequential"*; parallel design branches are described as *"Rare"*.
- `compress` is a terminal branch: *"treat `compress` output (`synthesis_df`) as a terminal branch in design mode. You usually don't chain it into `barcode`, `primer`, or `spacer`."*
- `expand` blows up exponentially and needs an `expansion_limit`.
- Motif matching in `verify` is literal substring matching — *"IUPAC bases are not expanded as wildcards."*
- `final` concatenates all columns and drops metadata; users are told to save an annotated DataFrame beforehand.
- No native Windows support (WSL required).

---

## 9. Confidence notes

- **High confidence (multiple independent sources, direct quotes):** absence of genome/GTF/HGVS/VCF/consequence/codon features; DataFrame-based data model; barcode and primer capabilities; synthesis-constraint checking; license; maintenance dates; DAG/YAML pipelines; variants-as-input.
- **Judgement calls, flag if challenged:**
  - `library_as_object = partial` — one could argue "yes" (a DataFrame *is* the whole library and no user looping is required) or "no" (there is no library abstraction, just a table). "Partial" is the defensible middle; the evidence string states both readings.
  - `assay_mpra = yes` vs `assay_dms = partial` — both are "infrastructure yes, content no". MPRA gets "yes" because three real MPRA libraries were built end-to-end and a dedicated MPRA template ships; DMS gets "partial" because only a compression recipe exists and there is zero protein/codon awareness. If the response wants a single consistent rule, both could be reported as "partial (infrastructure only)".
  - `primer_design = yes` depends on reading the row as "designs primers/oligos for wet-lab protocols". Under the strict reading "**mutagenic** primers", it is no. The evidence string makes the distinction explicit.
  - `barcode_generation = yes` despite **no GC-content constraint** — the row text names GC explicitly and that specific sub-feature is genuinely absent.
- **Not verified by execution** (per the no-install constraint): all behavior claims come from the paper, README, docs.md, api.md, agent-skills.md, and the example READMEs/notebook. Nothing was installed or run.
