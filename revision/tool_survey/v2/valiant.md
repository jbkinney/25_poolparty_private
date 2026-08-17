# VaLiAnT — v2 capability record

**Slug:** `valiant`
**Full name:** VaLiAnT — Variant Library Annotation Tool
**Citation key:** `Barbon2022th`
**Version assessed:** **4.0.0** (tag and last commit 2024-04-22). Re-scored **2026-08-10**.
**Row set:** v2 (`ROWS_v2.md`) — 33 keys. `hgvs_input` dropped; `library_as_object` split into
`library_first_class_object` + `library_algebra`; `mixed_mutagenesis_one_pool` split into four;
`dag_chaining` → `composable_operations`, `lazy_evaluation` → `lazy_generation` (renames).

**Basis:** `final/valiant.md` (FINAL, adversarially reviewed, re-verified against the live repo,
v4.0.0 source, GitHub API, quay.io, PyPI and anaconda.org on 2026-08-10). No new extraction was
performed for this pass; every value below traces to evidence already in that record. Four v2 rows
(`degenerate_iupac_codons`, `negative_control_generation`, `ml_model_in_loop`, `readout_analysis`)
were not asked as rows in v1 but are answered by the FINAL record's exhaustive-absence evidence base
(875-line README, CHANGELOG, `pyproject.toml`, 8 wiki pages, all shipped examples, and the full
461-path recursive git tree of `develop`).

> **Framing note carried forward.** VaLiAnT is the closest tool in this survey to PoolParty's problem
> space and its authors are plausible referees. Nothing here is falsely denied. Where the v2 wording
> lets VaLiAnT score *higher* than v1 it does (`heterogeneous_components_one_library` = yes,
> `exhaustive_single_scans` = yes). Where the split exposes a genuine boundary it says so with
> code-level evidence. The table must be labelled **VaLiAnT 4.0.0**, not the published 2022 version.

---

## Summary table

| Block | Key | v2 value |
|---|---|---|
| A | `library_first_class_object` | no |
| A | `composable_operations` | no |
| A | `lazy_generation` | no |
| A | `library_algebra` | no |
| A | `exhaustive_single_scans` | **yes** |
| A | `sampled_random_mutagenesis` | no |
| A | `higher_order_combinatorial` | no |
| A | `heterogeneous_components_one_library` | **yes** |
| A | `combinatorial_motif_place` | no |
| A | `barcode_generation` | no |
| A | `per_sequence_provenance` | **yes** |
| A | `automatic_naming` | **yes** |
| A | `design_visualization` | no |
| B | `assay_dms` | **yes** |
| B | `assay_mpra` | partial |
| B | `assay_insilico` | no |
| C | `genome_coordinates` | **yes** |
| C | `transcript_models` | **yes** (one transcript per execution) |
| C | `exon_intron_split_codons` | **yes** |
| C | `vcf_vep_output` | **yes** (SGE mode only) |
| C | `consequence_annotation` | **yes** |
| D | `primer_design` | no |
| D | `codon_optimization` | partial |
| D | `synthesis_constraints` | partial (length only) |
| D | `degenerate_iupac_codons` | no |
| D | `negative_control_generation` | no |
| D | `ml_model_in_loop` | no |
| D | `readout_analysis` | no |
| E | `interface` | yes — CLI only (Docker/Singularity; no Python API) |
| E | `license` | yes — AGPL-3.0-or-later |
| E | `installable_today` | yes — source + `quay.io/wtsicgp/valiant:4.0.0`; not on PyPI/bioconda |
| E | `last_activity` | yes — 2024-04-22 (commit + tag 4.0.0); ~27.6 months dormant |
| E | `documented_examples` | yes — 5 runnable examples (md5-validated) + 8 wiki pages |

---

## Block A — library specification

### `library_first_class_object` — **no**
*(v1 `library_as_object` = partial; the "partial" was carried entirely by the declarative-spec half,
which now lives in `heterogeneous_components_one_library`, and by the library-combining half, which
now lives in `library_algebra`.)*

There is no library object a user can hold, inspect, transform or pass onward. VaLiAnT is a CLI that
writes files: `MetaTable.to_csv(conn, targeton_name)` writes **one CSV per targeton**, and Supp.
Table 2 marks `_meta.csv`, `_meta_excluded.csv`, `_unique.csv` and both VCFs "targeton-specific",
with only `ref_sequences.csv` "execution-specific". The package is importable but exposes **no
documented Python API** — `pyproject.toml` has exactly one console script
(`valiant = "valiant.__main__:main"`) and **no `[project.entry-points]`**; there is no `docs/`
directory in the 461-path recursive tree, no API reference in the 875-line README, and no API page
among the 8 wiki pages. Internally the library exists only as rows in an **in-memory SQLite database**
that is drained to CSV and discarded (`db.py: with sqlite3.connect(':memory:') as conn`).
*Credit where due:* the *specification* is genuinely declarative and inspectable — one BED-like TSV,
one row per targeton — and one invocation processes many targetons across many contigs
(`run_sge` → `proc_contig` → `proc_targeton`), with no user loops. That strength is scored under
`heterogeneous_components_one_library`, not here; this row asks about the library as a manipulable
object, and the answer is no.
*Source:* `src/valiant/meta_table.py`; `src/valiant/db.py:131-132`; `src/valiant/sge_proc.py`;
`pyproject.toml`; paper Supp. Table 2; recursive git tree of `develop`.

### `composable_operations` — **no**  *(rename of `dag_chaining`; value and evidence carried unchanged — code-verified)*
No mechanism exists to chain, nest or compose design steps. `src/valiant/mutator_type.py` is a
**closed `Enum` of exactly seven mutator types**, whose only composition in the entire system is the
hard-coded `DEPENDENT_MUTATOR_TYPES = {MutatorType.SNV_RE: {MutatorType.SNV}}` (`snvre` internally
consumes `snv`). The pipeline order is fixed in `proc_contig` / `proc_targeton`: background variants →
PAM/protospacer protection edits → mutators on r1/r2/r3 → custom VCF variants → `--adaptor-5`/`-3`
appended → min/max-length filter → order-invariant dedup. Mutators act independently on disjoint
sub-regions of a single reference: "Regions 1-3 (r1-3) can be changed by mutator functions — **each
independently of each other**" (§2.4.1). No plugin hook; the authors state extension means editing the
source: "We envision that, as opensource software, VaLiAnT will be further modified by the community.
For example, by the addition of new mutator actions" (Discussion).
*Source:* `src/valiant/mutator_type.py`; `src/valiant/sge_proc.py`; `pyproject.toml`; paper §2.4.1,
§2.5.3, Discussion.

### `lazy_generation` — **no**  *(rename of `lazy_evaluation`; value and evidence carried unchanged — code-verified)*
`run_sge` opens an in-memory SQLite connection, inserts exons / PAM / custom / background variants,
calls `targeton.process(...)` to **materialise every mutation into the database**, then
`MetaTable.to_csv(...)` drains it to file. No generator, no streaming, no partial-materialisation
option, and no API through which a subset could be requested. The only short-circuit in the whole flow
is `--sequences-only`, which writes the reference-QC file and `sys.exit(0)`. Full library counts are
reported per run (paper Table 2, Supp. Table 3).
*Source:* `src/valiant/sge_proc.py` (~lines 535–568); `src/valiant/db.py:131-132`;
`src/valiant/data/ddl.sql`; `src/valiant/meta_table.py`; paper §2.5.3, Table 2.

### `library_algebra` — **no**  *(new row from the `library_as_object` split)*
Whole libraries cannot be combined, sampled or repeated as operations inside the tool. Outputs are
per-targeton CSV/VCF files; the paper's own phrase for the cDNA design is "The total number of variants
to cover all of the CDS and making up the **final concatenated library** is 62 754" (§3.2) —
**concatenation is the user's job, done outside the tool**, which is exactly the condition the row
definition names for "no". There is no stack/concat/merge operation, no sampling operation (see
`sampled_random_mutagenesis`), and no replicate/repeat operation: `--include-no-op-oligo` adds a
**single** WT sequence (CHANGELOG 4.0.0), not a repeat count. The only pooling that happens inside the
tool is *within* one targeton, where all mutators plus custom variants land in one `_meta.csv` and one
order-invariantly deduplicated `_unique.csv` — that is generation into one pool, not an algebra over
libraries.
*Source:* paper §3.2, Supp. Table 2; `src/valiant/meta_table.py`; CHANGELOG 4.0.0, 2.0.0; README
§*Oligonucleotide metadata file*.

### `exhaustive_single_scans` — **yes**  *(new row from the `mixed_mutagenesis_one_pool` split — VaLiAnT's core competence)*
This is precisely what every VaLiAnT mutator does. `snv` performs saturation single-nucleotide
substitution across the target region; `del` is **fully parametric** —
`<SPAN>del[<OFFSET>]`, "Non-overlapping stretches of nucleotides of a given length are deleted starting
from a given offset", with span and offset settable to "any valid value" (CHANGELOG 4.0.0;
`IntPatternBuilder(offset, span)` with `span > 0` the only constraint), giving arbitrary tiled deletion
scans (`1del`, `2del0`, `2del1`, `3del1`, `6del0`, `10del2` …); `aa` replaces each codon with the
top-ranking codon of every other amino acid — "**19 mutated sequences for each codon** mapping to an
amino acid ... and 20 for each stop codon" (README); `ala`, `stop`, `inframe` and `snvre` are likewise
exhaustive over every eligible codon. Measured outputs confirm exhaustiveness: BRCA1 exon 2 = 1000
total / 583 unique sequences from ten mutator/region combinations (paper Table 2); unique counts for
exons 2–5 of 583/740/825/1185 (nucleotide library) and 340/500/580/918 (amino-acid library, Supp.
Table 3); 62 754 variants across the 40-targeton cDNA design.
*Source:* README §*Parametric deletion*, §*All amino acid codon substitution*, §*SNVRE*;
`src/valiant/mutator_type.py`, `int_pattern_builder.py`; CHANGELOG 4.0.0; paper Tables 1–2, Supp.
Table 3, §3.2.

### `sampled_random_mutagenesis` — **no**  *(new row from the split)*
Every available mutation type is an exhaustive, deterministic, single-event scan. There is **no
sampled or random mutagenesis mode**: no mutation-rate parameter, no sampling depth, no RNG or seed
anywhere in the CLI option tables, the JSON configuration schema (`minOligoLength`, `maxOligoLength`,
`reverseComplementOnMinusStrand`, `forceBackground*`, file paths), the 8 wiki pages, or the 461-path
recursive tree. The closed 7-member `MutatorType` Enum contains no stochastic member. Library size is
determined entirely by the target-region length and the chosen mutators — reproducibility is in fact an
explicit design goal (order-invariant deduplication, md5-validated expected outputs).
*Source:* `src/valiant/mutator_type.py`; README CLI tables; shipped
`chr17_31226400_31226678_plus_sgRNA_1146241047_valiant_config.json`; recursive git tree; CHANGELOG
2.0.0; `examples/README.md`.

### `higher_order_combinatorial` — **no**  *(new row from the split)*
There is **no pairwise or higher-order mutator**. Each generated oligo carries exactly one designed
mutation event drawn from one mutator; multi-base and combinatorial variants can only enter via a
**hand-authored custom-variant VCF**, and even there README restricts to "simple variants ...
substitutions, insertions, deletions, indels" classified from `POS`/`REF`/`ALT` alone. No enumeration
of mutation × mutation combinations exists in `src/valiant/mutators/`
(`codon.py, deletion.py, snv.py, snv_re.py, utils.py`) or anywhere in the recursive tree.
*Honest adjacent fact, stated pre-emptively:* an oligo can carry more than one sequence change relative
to the reference — PAM/protospacer protection edits are applied to **every** oligo of a tagged
targeton, and background variants (v4.0.0) are laid into the reference frame before mutation. But those
are **fixed, uniform context edits**, not combinations being enumerated: the protection edit set is the
same across the whole targeton, and mutations that *overlap* background variants are **discarded**, not
combined ("Mutations that overlap background variants are discarded; warnings are always raised",
README). A parametric `3del` is one multi-base event, not a combination of events. The row asks for
pairwise and higher-order combinations *of mutations* in the same sequence; VaLiAnT does not generate
them.
*Source:* `src/valiant/mutator_type.py`, `src/valiant/mutators/` listing; README §*Custom variants*,
§*Background variants*, §*PAM protection*; paper Table 2; recursive git tree.

### `heterogeneous_components_one_library` — **yes**  *(new row from the split — this is where VaLiAnT's headline strength lands; score it generously)*
Structurally different component types genuinely coexist in one specification and are emitted into one
pool. A single targeton row carries a three-part action vector, so different mutator types run on
different sub-regions and several mutators run per region. Verbatim from the shipped
`brca1_nuc_targeton_input.txt` (exon 2 row):

```
"(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"
```

and from the prime-editing example (`input/input_a.txt`), six mutator types in r2 and three in r1:

```
"(1del,snv,2del0),(1del,snvre,inframe,stop,ala,aa),()"
```

The component types so mixed are structurally distinct: single-nucleotide substitutions, parametric
multi-base deletions at chosen span/offset, codon-level replacements (`ala`, `stop`, `aa`), in-frame
codon deletions, synonymous-recoding substitutions (`snvre`) — **plus** custom variants ingested from
real catalogues (ClinVar 2020-11-07 and gnomAD v3.0 via a manifest CSV), which merge into the same pool
across the whole targeton **including constant regions** (flagged per oligo by `vcf_var_in_const`),
**plus** an optional single WT no-op oligo (`--include-no-op-oligo`, CHANGELOG 4.0.0). Paper Table 2
shows the pooled composition for BRCA1 exon 2: 1000 total sequences, 583 unique, drawn from ten
mutator/region combinations plus custom variants. One file, one invocation, one pooled output — the
`brca1_nuc` input has four targeton rows and the cDNA input has 40, all in single runs.
*Boundary to record honestly, since the v2 row's illustrative example names it:* the heterogeneity is
over *exhaustive scan types and catalogue variants*. It cannot include sampled components or
higher-order components, because VaLiAnT has neither (rows above). Also, pooling is per targeton —
outputs are one CSV per targeton, concatenated by the user (see `library_algebra`).
*Source:* `examples/sge/brca1/parameter_input_files/brca1_nuc_targeton_input.txt`;
`examples/sge/brca1_prime_editing/input/input_a.txt`, `input_b.txt`;
`examples/cdna/input/cdna_targeton.tsv` (40 data rows); paper Tables 1–2, §2.4.1; CHANGELOG 4.0.0;
`src/valiant/meta_table.py` (`vcf_var_in_const`).

### `combinatorial_motif_place` — **no**  *(unchanged from FINAL, including its corrected evidence)*
VaLiAnT has no motif or element concept and no insertion-scan machinery. Fixed sequence enters an oligo
two ways: (1) `--adaptor-5` / `--adaptor-3`, appended once to **every** oligo in the run — "targetons
need to be processed individually as this function **appends to all sequences in the library**"
(§2.4); (2) as an **insertion in a user-supplied custom-variant VCF** — README §*Custom variants*:
"Only simple variants such as the following are supported: substitutions, **insertions**, deletions,
indels", so any sequence including a TF-binding motif can be inserted at any position and VaLiAnT will
build, name, annotate and emit an oligo for it. The value is nonetheless **no** because the
*combinatorial* half is absent: no motif × position enumeration, no orientation or permutation control,
no spacing or copy-number sweep, no element-library concept — the user must author every insertion by
hand. Confirmed against the full 461-path recursive tree: no motif, element or insertion-scan module
anywhere in `src/valiant/` (93 files).
*Source:* README §*Custom variants*, §*Adaptor sequences*; paper §2.4; recursive git tree of `develop`.

### `barcode_generation` — **no**
No barcode functionality of any kind. Absent from the 32-field `META_CSV_FIELDS` schema, from both CLI
option tables, from the JSON configuration schema, from all 8 wiki pages, and from the full recursive
tree. There is **no GC, edit-distance, Tm or uniqueness machinery** anywhere in the source. VaLiAnT's
answer to library identity is the variant itself plus the amplicon/adaptor structure, not a barcode.
*Source:* `src/valiant/meta_table.py`; README metadata + CLI tables; recursive git tree.

### `per_sequence_provenance` — **yes**  *(strong; provenance is bidirectional)*
Every oligo gets a row in `_meta.csv` with **32** structured columns (README, verbatim header):

```
oligo_name, species, assembly, gene_id, transcript_id, src_type, ref_chr, ref_strand, ref_start,
ref_end, revc, ref_seq, pam_seq, vcf_alias, vcf_var_id, mut_position, ref, new, ref_aa, alt_aa,
mut_type, mutator, oligo_length, mseq, mseq_no_adapt, pam_mut_annot, pam_mut_sgrna_id, mave_nt,
mave_nt_ref, vcf_var_in_const, background_variants, background_seq
```

That is provenance well beyond naming the mutation: `mutator` (which mutator built it),
`vcf_alias`/`vcf_var_id` (which custom VCF and which record — "To preserve variant provenance, an alias
is assigned to each VCF file", §2.5.2), `pam_mut_annot`/`pam_mut_sgrna_id` (which PAM edits, which
sgRNA), `background_variants`/`background_seq`, `ref_seq`/`pam_seq` (the reference and PAM-protected
reference the oligo was built from), `mave_nt`/`mave_nt_ref` (MAVE-HGVS), and `vcf_var_in_const`.
**Bidirectional:** the output VCFs carry `INFO` tags linking each variant record *back* to its oligo —
`SGE_OLIGO`, `SGE_SRC`, `SGE_VCF_ALIAS`, `SGE_VCF_VAR_ID`. Since v3.0.0 each execution also writes a
`config.json` recording the full parameter set.
*Source:* `src/valiant/meta_table.py` (`META_CSV_FIELDS`); README §*Oligonucleotide metadata file*,
§*Variant call files*; paper §2.5.2; CHANGELOG 3.0.0; `examples/sge/brca1_prime_editing/output_a_exp/`.

### `automatic_naming` — **yes**
`oligo_name` is generated by the tool; README example
`ENST00000357654.9.ENSG00000012048.23_chr17:43104102_1del_rc`
(transcript.gene_chr:pos_mutator_rc), implemented as `get_sge_oligo_name`, `get_cdna_oligo_name`,
`get_sge_oligo_no_op_name` in `src/valiant/meta_table.py`, with a documented `NO_TRANSCRIPT`
placeholder (wiki *Advanced usage*). Output **filenames** are auto-generated too: SGE mode from the
design — "The file names report chromosome, coordinates, strand and the sgRNA IDs associated with the
targeton" (§2.5.3), e.g. `chr17_43115722_43115768_minus_sgRNA_ex2_a_meta.csv`; cDNA mode hash-derived,
`<seq_id>_<md5>_meta.csv`. Duplicate sequences resolve deterministically: "When multiple
oligonucleotides have the same sequence, the first name in **lexicographic order** is chosen".
*Source:* `src/valiant/meta_table.py`; README; wiki *Advanced usage*; paper §2.5.3, Supp. Table 3;
shipped `output_a_exp/`, `examples/cdna/output_exp_v3/`.

### `design_visualization` — **no**
CSV / VCF / JSON output only. The complete pinned dependency set is
`charset-normalizer==3.3.2, click==8.1.7, pydantic==2.7.0, pysam==0.22.0` — **no plotting library** —
and there is no plotting or rendering module in the 461-path tree. The paper's own design figure was
made in third-party software (Fig. 3 legend: "sequence information modified from **Geneious Prime**
(version 2019.04) visualization"), and the Discussion places visualization downstream and outside the
tool. *Pre-emptive clause:* the repo ships static PNGs under `examples/images/wiki/` and a hand-built
GenBank plasmid model (`examples/cdna/construct_generation/BRCA1_cdna_pCW57.1_model.gb`) — these are
hand-authored documentation assets, not tool output.
*Source:* `pyproject.toml`; recursive git tree; paper Fig. 3 legend, Discussion.

---

## Block B — assay coverage

### `assay_dms` — **yes**
The tool's entire purpose, stated in its title: "an oligonucleotide library design and annotation tool
for **saturation genome editing** and other **deep mutational scanning** experiments." Both subcommands
are DMS modes (`sge` for endogenous-locus SGE / prime editing, `cdna` for expression-cassette DMS), and
the repo ships four worked BRCA1 DMS/SGE libraries with md5-validated expected outputs.
*Source:* paper title, abstract, §2.1–§3.4; `examples/`.

### `assay_mpra` — **partial**  *(attacked from both directions in review; survived both)*
No MPRA-specific support and no MPRA claim — MPRA appears in the paper only as background. **But**
non-coding targetons genuinely work: wiki *Advanced usage* — "UTR regions and non-coding sequence are
treated as if intronic sequence. Targetons that do not overlap any GTF/GFF2 feature won't be associated
with a gene or transcript..."; wiki *Input files* — "If CDS-specific mutations are not required, then a
GTF/GFF2 file need not be supplied"; and `sge_proc.proc_contig` has an explicit `if config.gff_fp: ...
else:` branch for the un-annotated case. So `snv` / parametric-deletion saturation of a regulatory
element produces usable MPRA-style oligos — but there are no barcodes, no reporter-cassette assembly,
no motif insertion or scanning, no promoter/enhancer grammar, and no MPRA example in the docs.
*Source:* paper §1; wiki *Advanced usage*, *Input files*; `src/valiant/sge_proc.py:448`.

### `assay_insilico` — **no**
Nothing in the paper, README, wiki or source relates to designing libraries to probe sequence-to-
function models: no model interface, no scorer, no covariate/design-matrix export. The single ML
mention in the Discussion runs the other way — DMS data later *training* models.
*Source:* paper Discussion; README; all 8 wiki pages; recursive git tree.

---

## Block C — genomics integration

### `genome_coordinates` — **yes**
`valiant sge` is coordinate-driven end to end: "Targeton ranges are defined by 1-based indexing of
genomic coordinates using 'ref_start' and 'ref_end' fields, with reference chromosome and strand also
defined, as 'ref_chr' and 'ref_strand'" (§2.1). Confirmed in the shipped targeton TSV
(`chr17  -  43115634  43115878  43115726  43115779 ...`). Sequences are fetched from a local FASTA+FAI
via pysam. Species and assembly are required positional CLI arguments; "The tool is species and genome-
build agnostic" (§2.3).
*Source:* paper §2.1, §2.3; `brca1_nuc_targeton_input.txt`; `src/valiant/sge_proc.py`;
`examples/sge/brca1/run_brca1_nuc.sh`.

### `transcript_models` — **yes (one transcript per execution)**
"To collect gene and transcript information and to apply CDS-specific mutator functions, appropriate
transcript annotation must be provided via a GTF/GFF2 file; only CDS, UTR and stop features are taken
into consideration" (§2.4.2). Wiki *Input files* is stricter and matches the code: "**One specific
transcript must be supplied for an execution.** GTF/GFF files should not contain multiple
transcripts... CDS and UTR features are required. The last CDS feature gets extended by one codon to
include STOP in computation." Code: `proc_contig` builds a single `Transcript` per contig/strand and
carries `# TODO (future multi-transcript)`; `run_sge` warns "Annotation provided for targetons on
multiple contigs: verify only one transcript is used!" Real exon/CDS/UTR-aware handling, but one
transcript per run — not a full transcript-model engine.
*Source:* wiki *Input files*; paper §2.4.2; `src/valiant/sge_proc.py:448-470`, `:549`;
`src/valiant/transcript.py`, `loaders/gtf.py`.

### `exon_intron_split_codons` — **yes**
Explicitly handled, with a published formula (eqns 1–2) and Supp. Fig. S2: "Retrieval of additional
positions from the reference sequence is necessary to obtain the context of partial codons... **The
extra bases required at either end of the target may come from the same or an adjacent exon**." SNVs
reach partial liminal codons and are annotated using adjacent exonic bases, while codon-replacing or
codon-deleting mutators apply only to complete codons (§2.5.1). v4.0.0 goes beyond the paper —
CHANGELOG 4.0.0 *Added*: "Support for out-of-frame CDS targeton regions whose 5' and 3' extensions span
**both adjacent and distal bases**." Documented restriction (§2.4.2): no single region r1/r2/r3 may
span both coding and non-coding sequence, though a targeton as a whole may.
*Source:* paper §2.5.1, §2.4.2, eqns 1–2, Supp. Fig. S2; CHANGELOG 4.0.0; `src/valiant/cds_seq.py`,
`exon.py`, `annot_variant.py`.

### `vcf_vep_output` — **yes (SGE mode only)**
Since v3.0.0 there are **two VCFs per targeton**: README §*Variant call files* — "The `REF` field
reports the reference sequence including (`*_pam.vcf`) or excluding (`*_ref.vcf`) PAM protection
edits", confirmed in shipped expected outputs
(`chr17_43115722_43115768_minus_sgRNA_ex2_a_pam.vcf` and `..._ref.vcf`). Each record carries
`SGE_OLIGO`, `SGE_SRC`, `SGE_VCF_ALIAS`, `SGE_VCF_VAR_ID` `INFO` tags, and **all** generated variants
are emitted, not just custom ones. *Caveats kept:* **cDNA mode emits no VCF** (relative coordinates),
the metadata CSV deliberately breaks VCF convention (§2.5.4), and **VEP itself is never named by the
authors** — the claim is standard-VCF compatibility, not a VEP integration.
*Source:* README §*Variant call files*; CHANGELOG 3.0.0;
`examples/sge/brca1_prime_editing/output_a_exp/`; paper §2.5.4, Supp. Table 2.

### `consequence_annotation` — **yes**
Codon-level consequence is computed and stored in `ref_aa`, `alt_aa`, `mut_type`: "SNVs carry extra
information when introduced into CDS targets ... namely whether the SNV results in a **synonymous,
missense or nonsense** change" (§2.5.1); strand-aware, precomputed as a codon-indexed table at start-up
(1152 codon/SNV combinations). v4.0.0 widens this beyond SNVs — CHANGELOG: "report the type of mutation
for the `aa`, `ala`, and `stop` mutators (`mut_type` field)" and "in cDNA mode, mutations preserving
stop codons are now annotated as nonsense". `pam_mut_annot` classifies each PAM edit as
`syn|mis|non|ncd`; background variants are classified and gated in code
(`is_variant_nonsynonymous`, `is_variant_frame_shifting`, with `--force-bg-ns` / `--force-bg-indels`).
*Limits (kept):* no splice-site consequence prediction, no frameshift/NMD calling, no external
annotator invoked — annotation is codon arithmetic.
*Source:* paper §2.5.1; CHANGELOG 4.0.0; README metadata table (field 26), §*Background variants*;
`src/valiant/sge_proc.py`, `targeton.py`, `annot_variant.py`, `codon_table.py`.

---

## Block D — adjacent / complementary

### `primer_design` — **no**
VaLiAnT appends **user-supplied** constant sequences (`--adaptor-5`, `--adaptor-3`) verbatim and leaves
flanking constant regions for user-chosen primers; it never designs one. In the prime-editing example
the P5/P7, Golden Gate Type IIS ends, sgRNA spacer, pegRNA scaffold and PBS are hand-written into the
adaptor strings (`run_prime_a.sh`: a 120-nt `--adaptor-5` literal). cDNA example: "The 21-base flanking
regions at targeton boundaries remain unmutated to facilitate primer binding" — the user picks them. No
Tm, specificity or dimer calculation in the tree; no primer3-class dependency among the four pinned
deps.
*Source:* README CLI table; `examples/sge/brca1_prime_editing/run_prime_a.sh`; wiki *cDNA example*;
`pyproject.toml`; recursive git tree.

### `codon_optimization` — **partial**  *(deliberately generous toward VaLiAnT — the safe direction)*
Codon-usage-aware but does not optimise a sequence. Ships a human codon-frequency table
(`src/valiant/data/default_codon_table.csv`; triplet / amino acid / frequency / rank) and accepts
`--codon-table` for other species. Ranks drive **per-codon substitution**: `ala`/`stop` use the
top-ranking triplet; `aa` replaces each codon with the top-ranking codon of every other amino acid;
`snvre` swaps a missense codon for the next most frequent triplet giving the same amino-acid change,
"allow[ing] for insights into the effect of codon sequence on missense changes". That is per-codon
frequency-ranked substitution — **not** whole-CDS optimisation, harmonisation, CAI maximisation, or
avoidance of structure/repeats/restriction sites.
*Source:* README §*All amino acid codon substitution*, §*SNVRE*, §*Alanine/stop codon substitution*;
`src/valiant/codon_table*.py`; `src/valiant/data/default_codon_table.csv`; CHANGELOG 4.0.0.

### `synthesis_constraints` — **partial (length only)**
Only oligo length is enforced: `--max-length` (default 300 bp) — "Oligonucleotides exceeding a given
length ... will not be included in the unique oligonucleotide files and their metadata will be stored
in separate files marked as 'excluded'"; `--min-length` added in v3.0.0; both present in the JSON
schema as `minOligoLength` / `maxOligoLength` (shipped NXRL config: `1` / `300`). The Discussion frames
it as a synthesis-platform constraint, and `_unique.csv` is described as ready for "Synthesis
submission" (Supp. Table 2). One adjacent *input-validity* rule, mentioned pre-emptively: "Ambiguous
nucleotides are not allowed in the reference sequence. Soft-masking is ignored." Still **no GC-content,
homopolymer, repeat, secondary-structure or restriction-site checking** anywhere in the tree.
*Source:* README CLI tables and §*Reference sequence*; paper Discussion, Supp. Table 2;
`examples/sge/nxrl/chr17_31226400_31226678_plus_sgRNA_1146241047_valiant_config.json`.

### `degenerate_iupac_codons` — **no**  *(NEW row; resolved from the FINAL record's absence evidence)*
VaLiAnT emits only concrete, fully specified sequences. The strongest positive evidence that it is
*not* a degenerate-codon designer is the `aa` mutator's own documented behaviour: it enumerates "**19
mutated sequences for each codon** mapping to an amino acid ... and 20 for each stop codon" (README) —
i.e. one explicit oligo per amino-acid substitution, exactly where a degenerate design (NNK/NNS) would
emit a single mixed-base oligo. Ambiguity is affirmatively rejected on input: "**Ambiguous nucleotides
are not allowed in the reference sequence.** Soft-masking is ignored." There is no IUPAC alphabet,
expansion or compression facility in the 7-member `MutatorType` Enum, in `src/valiant/mutators/`
(`codon.py, deletion.py, snv.py, snv_re.py, utils.py`), in the CLI option tables, in the JSON schema,
in the 32-field metadata header, in the 8 wiki pages, or in the 461-path recursive tree. Custom
variants are restricted to concrete "substitutions, insertions, deletions, indels" classified from
`POS`/`REF`/`ALT`.
*Source:* README §*All amino acid codon substitution*, §*Reference sequence*, §*Custom variants*;
`src/valiant/mutator_type.py`, `src/valiant/mutators/`; recursive git tree of `develop`.

### `negative_control_generation` — **no**  *(NEW row)*
No scramble, shuffle, dinucleotide-preserving permutation, or reverse/complement control generator
exists. There is no control-generation module in the 461-path tree, no such option in either CLI option
table or the JSON configuration schema, and no `src_type` / control category in the 32-field metadata
schema beyond the mutation sources actually designed.
*Two adjacent features, named honestly so a referee cannot claim they were missed:*
(1) `--include-no-op-oligo` (CHANGELOG 4.0.0) adds a **single** unmutated WT reference oligo — a
reference control, not a scramble/shuffle, and not replicated; (2) `--revcomp-minus-strand` emits
reverse-complemented oligos for minus-strand targets — that is strand-correct construction for
ssODN/SGE orientation, not a reverse-complement *control* arm. Neither is the capability this row asks
about.
*Source:* CHANGELOG 4.0.0; README CLI tables, §*Strand*; `src/valiant/meta_table.py`; recursive git
tree; shipped NXRL `*_valiant_config.json`.

### `ml_model_in_loop` — **no**  *(NEW row; same evidence as `assay_insilico`)*
No predictive model participates in design at any point. No model interface, scorer, weights file,
inference dependency (the complete pinned set is `charset-normalizer, click, pydantic, pysam`), or
covariate/design-matrix export anywhere in the paper, README, 8 wiki pages, or the recursive tree.
Sequence selection is deterministic enumeration from the action vector. The one ML mention in the
Discussion runs the other way — DMS data later *training* models ("as predicative models and machine
learning increase in utility ... training datasets rich in biological context will be important").
*Source:* paper Discussion; `pyproject.toml`; README; all 8 wiki pages; recursive git tree.

### `readout_analysis` — **no**  *(NEW row)*
VaLiAnT is a design-and-annotation tool only: it stops at `_unique.csv` for "Synthesis submission"
(Supp. Table 2). There is no FASTQ/BAM/count ingestion, no variant-count or enrichment/score
computation, and no analysis module in the recursive tree; `pysam` is present solely to read the
reference FASTA/FAI and VCFs on the *input* side. The authors place readout analysis explicitly
downstream and outside the tool: "VaLiAnT may be combined with downstream analysis software tools ...
and to produce annotated visualizations of variant effect" (Discussion). The MAVE-HGVS output strings
(`mave_nt`, `mave_nt_ref`) exist precisely to hand off to MaveDB-style downstream analysis — an
interoperability hook, not analysis.
*Source:* paper Discussion, Supp. Table 2; README §*Oligonucleotide metadata file*; `pyproject.toml`;
recursive git tree.

---

## Block E — engineering and availability
*(cell value is "yes" per the v2 convention; the real answer is in the evidence)*

### `interface` — yes → **CLI only (Docker/Singularity; no Python API)**  — flagged, see below
"The tool is implemented as a standalone executable Python package exposing a command line interface"
(§2.3). `pyproject.toml` `[project.scripts]` has **exactly one** entry
(`valiant = "valiant.__main__:main"`) and **no `[project.entry-points]`**. Two subcommands (`sge`,
`cdna`) plus a whole-run JSON configuration entry point (`valiant -c config.json`, demonstrated by the
NXRL example). Distributed as source + Dockerfile + prebuilt `quay.io/wtsicgp/valiant`, with
**Singularity a first-class documented path** (its own README section with a full
`singularity exec --cleanenv -B ...` invocation). **No documented Python API** (no `docs/` directory,
no API section in the 875-line README, no API page among the 8 wiki pages), **no web service, no GUI**.
*Column-semantics carry-over from FINAL §8:* if the column means "ships a documented, usable
interface" → yes (matching the `seqpro` record's `yes (Python API only)` convention); if it means
"exposes a programmatic / library API" → the honest value is **no**. **Never print a bare `yes`**;
recommended cell text: `CLI only (Docker/Singularity; no Python API)`.
*Source:* paper §2.3; `pyproject.toml`; README §*Docker image*, §*Singularity image*;
`examples/sge/nxrl/run.sh`; recursive git tree; PyPI/anaconda API checks.

### `license` — yes → **AGPL-3.0-or-later**
Confirmed three ways: `gh api /repos/cancerit/VaLiAnT` → `"license": "AGPL-3.0"`; README LICENCE block
and per-file headers ("GNU Affero General Public License ... version 3 of the License, or (at your
option) any later version"); paper Availability — "VaLiAnT is licensed under AGPLv3."
*Presentation warning:* put the SPDX identifier in the cell. AGPL is strong copyleft with a network-use
clause, materially unlike the MIT licences of `seqpro` and `tangermeme` in the same survey.
*Source:* GitHub API; README; source-file headers; paper Availability.

### `installable_today` — yes → **source install + container; not on PyPI or bioconda**
- **Containers pull and run today:** `quay.io/wtsicgp/valiant:4.0.0` and `:latest`, both pushed
  2024-04-22 (older tags 3.0.1 = 2023-07-19, 3.0.0 = 2022-10-12, 2.0.0 = 2021-07-12,
  1.0.0 = 2020-12-21). Singularity documented as a first-class path.
- **From source:** `python3.11 -m venv .env && source .env/bin/activate && pip install .` Requires
  Python ≥ 3.11; Linux/macOS native, Windows via WSL/Docker/Singularity. Dependencies fully pinned
  (`charset-normalizer==3.3.2, click==8.1.7, pydantic==2.7.0, pysam==0.22.0`) — so it should still
  build, though pinned `pysam==0.22.0` may need a source build on very new interpreters.
- **Not on PyPI** — `pypi.org/project/valiant` is Duncan Dickinson's unrelated dependency-audit tool
  (v0.2.3, 2021). **Not on bioconda** — `api.anaconda.org/package/bioconda/valiant` returns
  `{"error": "\"valiant\" could not be found"}`. **No web service.**
*Source:* quay.io tag API, PyPI JSON API, anaconda.org API (all re-queried 2026-08-10); README install
sections; `pyproject.toml`.

### `last_activity` — yes → **2024-04-22 (commit + tag 4.0.0); ~27.6 months dormant at 2026-08-10**
Last commit on `develop` 2024-04-22T13:51:54Z ("Merge tag '4.0.0' into develop — Background variant
support", Luca Barbon); `pushed_at = 2024-04-22T13:52:12Z`. Tags 1.0.0 / 2.0.0 / 3.0.0 / 3.0.1 / 4.0.0
but **zero GitHub Release objects**. Docker `latest` and `4.0.0` both last pushed 2024-04-22. Repo is
public and **not archived**; 6 stars, 3 forks; both open issues are unattended Snyk dependency-bot PRs
(#10 2024-03-03, #11 2024-04-03) and there is no human issue traffic at all.
*Carry-over from FINAL:* under the v1 `maintained` row this record scored **partial**. The survey
should publish its criterion (suggested: commit within 24 months = yes; older but installable and not
archived = partial; archived/uninstallable = no) so the column is defensible across all tools.
*Source:* GitHub repo/commits/releases/issues API; quay.io tag API — re-queried 2026-08-10.

### `documented_examples` — yes → **5 runnable examples (md5-validated) + 8 wiki pages + 875-line README**
All under `examples/`, each with `run_*.sh`, `check_*.sh` and a committed expected-output directory;
validation is by `md5sum` against `output_exp/` — a genuinely strong reproducibility claim. Inputs must
be unpacked first (`sge/unpack_reference.sh`, `brca1_prime_editing/unpack_inputs.sh`).
1. **`brca1_nuc`** — BRCA1 exons 2–5 SGE nucleotide library, 4 targeton rows in one file / one
   invocation; exons 2–4 action vector `"(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"`,
   exon 5 `"(1del,2del0,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"` (its r1 is 20 bp, even).
   Flags actually passed: `--pam`, `--vcf`, `--revcomp-minus-strand`, `--adaptor-5`, `--adaptor-3`,
   `--gff` — **not** `--max-length 300` (300 bp is the default). Unique counts 583/740/825/1185.
2. **`brca1_pep`** — same four targetons, action vector `"(),(aa,inframe),()"`; unique counts
   340/500/580/918.
3. **`brca1_prime_editing`** — two pegRNAs over BRCA1 exon 2 (the RTT *is* the targeton);
   `"(1del,snv,2del0),(1del,snvre,inframe,stop,ala,aa),()"` and its mirror; `--max-length 250` **is**
   passed; 120-nt `--adaptor-5` and 38-nt `--adaptor-3` literals; 454 / 483 unique.
4. **`examples/cdna`** — 5592 bp BRCA1 CDS in pCW57.1; `cdna_targeton.tsv` has **exactly 40 data
   rows**, all with action vector `"snv,1del,snvre,ala,stop,inframe,aa"`; 62 754 variants in the "final
   concatenated library"; ships `output_exp/` and `output_exp_v3/` plus a hand-authored
   `construct_generation/` (evidence the 40-targeton tiling plan was written by hand).
5. **`examples/sge/nxrl`** — background-variant example, repo only; the **only worked demonstration of
   the JSON-config entry point** (`valiant -c ..._valiant_config.json`), setting
   `forceBackgroundNonSynonymous`/`forceBackgroundFrameShifting` true, `minOligoLength: 1`,
   `maxOligoLength: 300`, and **no** `maskBackgroundFilePath`.
Plus wiki *Docker-usage-example* (**stale** — still pins `VERSION=2.0.0`, four releases behind).
*Source:* `examples/README.md` and all five example directories (scripts, inputs, expected outputs);
8 wiki pages; README.

---

## Changes vs the v1 record

| v1 key | v1 value | v2 key(s) | v2 value(s) | Why |
|---|---|---|---|---|
| `library_as_object` | partial | `library_first_class_object` | **no** | The v1 "partial" was carried by the declarative-spec strength, which the split reassigns. Under the narrow wording — an object the user can hold/inspect/transform/pass onward — VaLiAnT writes per-targeton files and exposes no Python API. |
| `library_as_object` | partial | `library_algebra` | **no** | The paper itself says the "final concatenated library" is assembled by the user outside the tool; no stack/concat/sample/repeat operation exists in-tool. |
| `mixed_mutagenesis_one_pool` | yes | `exhaustive_single_scans` | **yes** | Carries the strength: `snv`, parametric `del`, `aa` (19/codon), `ala`, `stop`, `inframe`, `snvre` are all exhaustive. |
| `mixed_mutagenesis_one_pool` | yes | `sampled_random_mutagenesis` | **no** | v1 memo already recorded verbatim: "There is no sampled or random mutagenesis mode." No rate, depth, RNG or seed anywhere. |
| `mixed_mutagenesis_one_pool` | yes | `higher_order_combinatorial` | **no** | v1 memo already recorded: "no pairwise/higher-order mutator"; multi-base/combinatorial variants only via hand-authored VCF. |
| `mixed_mutagenesis_one_pool` | yes | `heterogeneous_components_one_library` | **yes** | Verified verbatim in shipped inputs: different mutator types across r1/r2/r3, several per region, plus ClinVar/gnomAD custom variants and an optional WT no-op, pooled into one output. |
| `dag_chaining` | no | `composable_operations` | **no** | Rename; value and evidence unchanged. |
| `lazy_evaluation` | no | `lazy_generation` | **no** | Rename; value and evidence unchanged. |
| `hgvs_input` | no | — | *(dropped)* | Zero discriminating power across the 12 tools. Note for the response letter: do **not** concede HGVS input to VaLiAnT — `mave_hgvs.py` is write-only. |
| `maintained` | partial | `last_activity` | yes → 2024-04-22 | Reframed from a judgement label to an observed date; the underlying facts are identical. |
| — | — | `installable_today` | yes (source + container) | NEW; answered from FINAL §9. |
| — | — | `documented_examples` | yes (5 + 8 wiki pages) | NEW; answered from FINAL §5. |
| — | — | `degenerate_iupac_codons` | **no** | NEW; answered from FINAL absence evidence. |
| — | — | `negative_control_generation` | **no** | NEW; `--include-no-op-oligo` (single WT) and `--revcomp-minus-strand` are adjacent but are not this capability. |
| — | — | `ml_model_in_loop` | **no** | NEW; same evidence base as `assay_insilico`. |
| — | — | `readout_analysis` | **no** | NEW; design/annotation only, analysis explicitly placed downstream. |

All other keys (`combinatorial_motif_place`, `barcode_generation`, `per_sequence_provenance`,
`automatic_naming`, `design_visualization`, Block B, Block C, `primer_design`, `codon_optimization`,
`synthesis_constraints`, `interface`, `license`) are unchanged from the FINAL record.

---

## Confidence flags (for human review)

1. **`library_first_class_object` = no** — the label moved from `partial` to `no` because the split
   reassigned the declarative-spec credit. Facts are certain (per-targeton CSV output; no Python API;
   in-memory SQLite discarded at exit); the label is the judgement call. A referee-author could argue
   the targeton TSV *is* the library. Recommend the cell text read `no (CLI; per-targeton files)`.
2. **`heterogeneous_components_one_library` = yes** — a generous but, I believe, correct reading:
   different *scan types* plus catalogue variants plus an optional WT no-op genuinely coexist in one
   specification and one pooled output. The v2 row's illustrative example ("exhaustive singles +
   sampled higher-order + WT replicates") lists two things VaLiAnT cannot do. If the survey intends the
   row to require *those specific* component types, the value would be `partial`. Decide once, apply to
   all 12 tools.
3. **`higher_order_combinatorial` = no** — high confidence on generation, but note the honest adjacent
   fact: oligos do carry PAM-protection edits and (v4) background variants alongside the designed
   mutation, so a finished oligo can differ from reference at more than one position. Those are fixed
   uniform context edits, not enumerated combinations, and overlapping mutations are *discarded*. Worth
   a footnote in the response letter to pre-empt the objection.
4. **`degenerate_iupac_codons` = no** — NEW row, resolved by absence rather than by a statement in the
   memo. The absence evidence is strong (full recursive tree, README, wiki, mutator Enum) and the `aa`
   mutator's 19-sequences-per-codon behaviour is affirmative counter-evidence, but no source explicitly
   says "no IUPAC support". Low residual risk.
5. **`negative_control_generation` = no** — NEW row, and the only one where a "partial" is arguable:
   `--include-no-op-oligo` is a first-class flag that adds a WT control sequence. It is a *reference*
   control, not a scramble/shuffle/revcomp control, so `no` is right under this row's wording — but if
   the survey later broadens the row to "control sequences as a first-class feature", VaLiAnT becomes
   `partial`. Flagging so the definition is applied consistently.
6. **`readout_analysis` = no** — NEW row, resolved by absence. Very high confidence (the Discussion
   explicitly places analysis downstream), but flagged as not directly asked in the original
   extraction.
7. **`ml_model_in_loop` = no** — NEW row, resolved by absence; effectively the same evidence as
   `assay_insilico = no`, which survived adversarial review. Low risk, flagged for completeness.
8. **`interface` = yes** — unchanged column-semantics problem from FINAL §8. Never print a bare `yes`;
   define the column in the caption.
9. **`assay_mpra` = partial** — carried unchanged; survived attack from both directions in review, but
   it remains the second genuine judgement call in the record.
10. **`vcf_vep_output` = yes** — VEP is never named by the authors, and cDNA mode emits no VCF. The row
    key overstates what the tool claims; consider renaming the row to `vcf_output` or captioning it.
11. **`last_activity`** — the value depends on the assessment date (2026-08-10). Re-check before
    submission if the response letter is filed much later.
12. **`codon_optimization` = partial and `synthesis_constraints` = partial** — both are deliberately
    generous partials erring toward VaLiAnT. Underlying facts are certain; the labels are chosen for
    safety, not precision.

*Version-drift warning (carried forward):* background variants, parametric deletions, two output VCFs,
MAVE-HGVS, `--include-no-op-oligo`, `mut_type` for `aa`/`ala`/`stop` and the JSON config all **postdate
the 2022 paper**. The table must be explicit that it describes **VaLiAnT 4.0.0**.
