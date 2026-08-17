# VaLiAnT — FINAL capability record

**Slug:** `valiant`
**Full name:** VaLiAnT — Variant Library Annotation Tool
**Citation key:** `Barbon2022th`
**Tier:** 2
**Authors / institution:** Luca Barbon, Victoria Offord, Elizabeth J. Radford, Adam P. Butler,
Sebastian S. Gerety, David J. Adams, Hong Kee Tan, Andrew J. Waters — Wellcome Sanger Institute
**Paper:** Barbon L. *et al.*, "VaLiAnT: an oligonucleotide library design and annotation tool for
saturation genome editing and other deep mutational scanning experiments," *Bioinformatics*
38(4):892–899 (2022), doi:10.1093/bioinformatics/btab776 (Advance Access 16 Nov 2021)
**Repo:** https://github.com/cancerit/VaLiAnT (default branch `develop`)
**Manual:** https://github.com/cancerit/VaLiAnT/wiki (11 content pages)
**Container:** `quay.io/wtsicgp/valiant` (Docker + Singularity)
**License:** AGPL-3.0-or-later
**Version assessed:** **4.0.0** (tag and last commit both 2024-04-22). Assessment date **2026-08-10**.

**Record status: FINAL.** The adversarial review returned **"supported" for 21 of 24 capability
values**, "understated" for 1 (`combinatorial_motif_place` — value correct, evidence false),
"unsupported as encoded" for 1 (`interface` — column-semantics problem), and "overstated" for 1
(`maintained`). **One value changed** in this pass (`maintained`: yes → partial). Every correction
below was independently re-verified against the live repo, the v4.0.0 source, the GitHub API, quay.io,
PyPI and anaconda.org on 2026-08-10 — not merely copied from the review.

> **Critical framing note for the referee response.** VaLiAnT is the closest tool in this survey to
> PoolParty's problem space, it is well engineered, and its authors are plausible referees. **No
> capability is falsely denied in this record**: every one of the eight "no" cells was attacked
> against the current v4.0.0 source (not just the 2022 paper) and held. Where VaLiAnT is strong
> (mixed mutagenesis in one pool, per-oligo provenance, automatic naming, CDS/split-codon geometry,
> VCF output) this record says so in full. The comparison table must be labelled as describing
> **VaLiAnT 4.0.0**, not the published 2022 version — several capabilities grew substantially after
> the paper (background variants, two output VCFs, parametric deletions, `mut_type` for `aa`/`ala`/
> `stop`), and the paper's own §2.3 dependency list (pandas, PyRanges) is already obsolete.

---

## 1. Sources consulted

| Kind | Ref | Retrieved |
|---|---|---|
| pdf | `../../../lit_review/analyzed/Barbon2022_VaLiAnT_all.pdf` (11 pp incl. Supp. Tables S1–S3, Supp. Figs S1–S2), text via PyMuPDF; all quoted passages re-extracted and re-matched for this memo | 2026-08-10 |
| prior_analysis | `revision/tool_survey/extractions/valiant.md ⟨deleted at cleanup — in commit 35d65d8⟩` (original extraction memo) | 2026-08-10 |
| prior_analysis | `revision/tool_survey/reviews/valiant.md ⟨deleted at cleanup — in commit 35d65d8⟩` (adversarial review) | 2026-08-10 |
| repo (raw) | `README.md` (875 lines, read in full) | 2026-08-10 (re-verified) |
| repo (raw) | `CHANGELOG.md` (all five releases), `pyproject.toml` | 2026-08-10 (re-verified) |
| repo (raw) | `src/valiant/mutator_type.py`, `src/valiant/int_pattern_builder.py`, `src/valiant/sge_proc.py` (568 lines), `src/valiant/db.py` | 2026-08-10 (re-verified) |
| repo (API) | `gh api /repos/cancerit/VaLiAnT` — `pushed_at`, `archived`, stars, forks, issues, license | 2026-08-10 (re-verified) |
| repo (API) | `gh api /repos/cancerit/VaLiAnT/commits?sha=develop`, `/releases`, `/issues` | 2026-08-10 (re-verified) |
| repo (API) | full recursive git tree of `develop` (`/git/trees/develop?recursive=1`, 461 paths) | 2026-08-10 (re-verified) |
| repo (raw) | `examples/README.md`; `examples/sge/brca1/run_brca1_nuc.sh`, `run_brca1_pep.sh` and both `parameter_input_files/*_targeton_input.txt`; `examples/sge/brca1_prime_editing/run_prime_a.sh`, `input/input_a.txt`, `input/input_b.txt`; `examples/sge/nxrl/run.sh` and its `*_valiant_config.json`; `examples/cdna/run.sh`, `input/cdna_targeton.tsv`, `input/cdna_annot.tsv` | 2026-08-10 (re-verified) |
| docs | Wiki raw markdown, all 11 content pages: `Home`, `Input-files`, `Output-files`, `Advanced-usage`, `Configuration`, `cDNA-example`, `cDNA-DMS-file-formats`, `Saturation-prime-editing`, `Saturation-prime-editing-example`, `Software-and-dependencies`, `Docker-usage-example` | 2026-08-10 (last 3 pages added 2026-08-14) |
| web | `quay.io/api/v1/repository/wtsicgp/valiant/tag/` — all six tags + timestamps | 2026-08-10 (re-verified) |
| pypi | `pypi.org/pypi/valiant/json` → **not this tool** (Duncan Dickinson's "Audit tool to help investigate Python dependencies", v0.2.3) | 2026-08-10 (re-verified) |
| web | `api.anaconda.org/package/bioconda/valiant` → `{"error": "\"valiant\" could not be found"}` | 2026-08-10 (re-verified) |

No install, no clone, no execution — document/web research only.

---

## 2. What VaLiAnT actually is

A command-line Python package (Python ≥ 3.11) that turns a declarative table of genomic *targetons*
into a synthesis-ready oligonucleotide pool with a per-oligo metadata row. Two subcommands:

- `valiant sge` — reads a local reference FASTA; **absolute genomic coordinates**; up to five regions
  per targeton (`c1-r1-r2-r3-c2`), where r1/r3 are derived from r2 by an extension vector.
- `valiant cdna` — reads a user-supplied multi-FASTA; **relative coordinates**; a single target
  region (r2) per row; several cDNAs (`seq_id`) each with several targeton rows in one run.

> "VaLiAnT is run from the command line." (paper §2.1)
> "The tool is implemented as a standalone executable Python package exposing a command line
> interface." (§2.3)

**Mutator functions** (README, `develop`, v4.0.0) — a closed set of seven types
(`src/valiant/mutator_type.py`, condensed — comments added, type annotations and the
`ANNOTABLE_MUTATOR_TYPES` set omitted):

```
class MutatorType(str, Enum):
    DEL = 'del'      # parametric span:      <SPAN>del[<OFFSET>]
    SNV = 'snv'
    SNV_RE = 'snvre' # codon span, CDS only
    IN_FRAME = 'inframe'
    ALA = 'ala'
    STOP = 'stop'
    AA = 'aa'
PARAMETRIC_MUTATOR_TYPES = {MutatorType.DEL}
CDS_MUTATOR_TYPES = {SNV_RE, ALA, STOP, AA, IN_FRAME}
DEPENDENT_MUTATOR_TYPES = {MutatorType.SNV_RE: {MutatorType.SNV}}
```

plus **custom variants imported from VCF** (labelled `custom`). `src/valiant/mutators/` contains only
`__init__.py, codon.py, deletion.py, snv.py, snv_re.py, utils.py`.

**Fixed processing pipeline** (hard-coded across `src/valiant/sge_proc.py`, `targeton.py` and
`meta_table.py`, not user-composable): background variants → PAM/protospacer protection edits →
custom VCF variants (incl. constant regions) → mutator functions on r1/r2/r3 →
`--adaptor-5`/`--adaptor-3` appended → min/max-length filter →
order-invariant dedup to `_unique.csv`.

**Engine:** an **in-memory SQLite database** (`db.py: with sqlite3.connect(':memory:') as conn`,
`data/ddl.sql`, `queries.py`, `sql_gen.py`). Note that the paper's §2.3 says tabular data is handled
by pandas and GTF by PyRanges — **both dependencies are gone in v4.0.0**, whose complete pinned
dependency set is `charset-normalizer==3.3.2, click==8.1.7, pydantic==2.7.0, pysam==0.22.0`.

---

## 3. Changes made in this final pass

| Key | Extraction | Review verdict | FINAL | What changed |
|---|---|---|---|---|
| `maintained` | yes | overstated | **partial** | **Value changed.** ~27.6 months dormant (last commit 2024-04-22), zero GitHub Release objects, both open issues are unattended Snyk bot PRs. Still not archived; source and container images remain published. |
| `combinatorial_motif_place` | no | understated (evidence false) | **no** (unchanged) | **Evidence rewritten.** The sentence "the only insertable fixed sequences are `--adaptor-5`/`--adaptor-3`" is **false** — custom-variant VCFs support insertions of arbitrary sequence at arbitrary positions. Value stays "no" because the *combinatorial* half is genuinely absent. |
| `interface` | yes | unsupported as encoded | **yes (CLI only)** — flagged | Value kept for cross-table consistency with the `seqpro` record's `yes (Python API only)` convention, but the cell text must carry the qualifier. See §8 *Unresolved / flagged*. |
| `library_as_object` | partial | supported | partial | Evidence strengthened: one invocation genuinely handles many targetons across many contigs. |
| `lazy_evaluation` | no | supported | no | Hedge dropped — now code-verified, not inferred. |
| `vcf_vep_output` | yes | supported (stale description) | yes | Refreshed from v4.0.0: **two** VCFs per targeton with oligo back-links, not one. |
| `exon_intron_split_codons` | yes | supported | yes | Added CHANGELOG 4.0.0 distal-base extension support. |
| `consequence_annotation` | yes | supported | yes | Added v4 `mut_type` for `aa`/`ala`/`stop`, cDNA stop-preserving → nonsense, `pam_mut_annot`, background-variant classification. |
| `per_sequence_provenance` | yes | supported | yes | Added bidirectional provenance (VCF `INFO` back-links) and `vcf_var_in_const`. |
| `transcript_models` | yes | supported | yes | Re-anchored on the wiki's stricter wording + the `# TODO (future multi-transcript)` code comment. |
| `license` | yes | supported (presentation warning) | yes (AGPL-3.0-or-later) | SPDX identifier moved into the cell. |
| all other 13 keys | — | supported | unchanged | Evidence enriched where the review supplied verified detail. |
| `documented_examples` | — | 3 checkable errors | corrected | NXRL flags, `--max-length 300`, exon 5 action vector — all three re-verified against the shipped files and fixed. |

---

## 4. Capability assessment

### Block A — library specification

**`library_as_object` — partial.**
The declarative half is real and stronger than the extraction stated. The whole library is declared in
one BED-like targeton TSV — one row per targeton, each carrying an `action_vector` naming the mutators
for r1/r2/r3 — and a **single invocation processes many targetons across many contigs**:
`run_sge` → `for contig, targeton_configs in exp_config.targeton_configs.items(): proc_contig(...)`
→ `for targeton_config in tcs: proc_targeton(...)` (`src/valiant/sge_proc.py`). The shipped
`brca1_nuc_targeton_input.txt` is one file with four targeton rows executed by one command; the cDNA
example's `cdna_targeton.tsv` has **40 data rows** in one run. The user writes no loops.
**But** there is no library object a user can hold, transform, inspect or pass onward.
`MetaTable.to_csv(conn, targeton_name)` writes **per-targeton files, not one library-wide file**;
Supp. Table 2 marks
`_meta.csv`, `_meta_excluded.csv`, `_unique.csv` and the VCFs "targeton-specific" and only
`ref_sequences.csv` "execution-specific". The paper's own phrase for the cDNA design is "The total
number of variants to cover all of the CDS and making up the **final concatenated library** is 62 754"
(§3.2) — concatenation is the user's job, done outside the tool.
*Source:* paper §2.4.1, §3.2; `src/valiant/sge_proc.py`; `src/valiant/meta_table.py`; wiki
`Input-files`, `Output-files`; `examples/sge/brca1/parameter_input_files/brca1_nuc_targeton_input.txt`;
`examples/cdna/input/cdna_targeton.tsv`.

**`dag_chaining` — no.** *(code-verified)*
No mechanism exists to chain, nest, or compose design steps. `src/valiant/mutator_type.py` is a
**closed `Enum` of exactly seven mutator types** whose only composition in the entire system is the
hard-coded `DEPENDENT_MUTATOR_TYPES = {MutatorType.SNV_RE: {MutatorType.SNV}}` (i.e. `snvre`
internally consumes `snv`). `pyproject.toml` declares **no `[project.entry-points]`** and no plugin
hook — only one console script. The pipeline order (background → PAM protection → custom VCF →
mutators → adaptors → length filter → dedup) is fixed in
`proc_contig`/`proc_targeton`/`Targeton.process`/`MetaTable.to_csv`, not user-specified.
Mutators are applied independently to disjoint sub-regions of a single reference: "Regions 1-3 (r1-3)
can be changed by mutator functions — **each independently of each other** — by detailing corresponding
mutator lists in the BED-like input file" (§2.4.1). The authors are explicit that extension means
editing the source: "We envision that, as opensource software, VaLiAnT will be further modified by the
community. For example, by the addition of new mutator actions" (Discussion).
*Source:* `src/valiant/mutator_type.py`; `src/valiant/sge_proc.py`; `pyproject.toml`; paper §2.4.1,
§2.5.3, Discussion.

**`lazy_evaluation` — no.** *(now code-verified — the extraction's "not verified in code" hedge is dropped)*
`run_sge` opens an **in-memory SQLite** connection (`db.py: with sqlite3.connect(':memory:') as conn`;
`init_db(conn)`; `src/valiant/data/ddl.sql`), inserts exons / PAM / custom / background variants,
calls `targeton.process(...)` to **materialise every mutation into the database**, and then
`MetaTable.to_csv(conn, targeton_name)` drains it to file. The only short-circuit in the whole flow is
`--sequences-only`, which writes the reference-QC file and `sys.exit(0)`. There is no partial
materialisation option and no public Python API through which a subset could be requested. (The only
row-at-a-time work is draining the finished database to file: `MetaTable.to_csv` walks the result set
with `fetchone`.)
*Source:* `src/valiant/sge_proc.py` (lines ~535–568), `src/valiant/db.py:131-132`,
`src/valiant/data/ddl.sql`, `src/valiant/meta_table.py`; paper §2.5.3, Table 2.

**`mixed_mutagenesis_one_pool` — yes.** *(VaLiAnT's headline strength — state it generously)*
A single targeton row carries a three-part action vector, so **different mutator types run on
different sub-regions and several mutators run per region**, all pooled into one output. Verified
against the shipped input file, not just paper Table 1 —
`brca1_nuc_targeton_input.txt`, exon 2 row, verbatim:

```
"(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"
```

The prime-editing example mixes six mutator types in r2 and three in r1 in one row:
`"(1del,snv,2del0),(1del,snvre,inframe,stop,ala,aa),()"` (`input/input_a.txt`). Custom ClinVar/gnomAD
variants are merged into the same pool across the whole targeton, including constant regions. Paper
Table 2 shows the pooled composition for BRCA1 exon 2: 1000 total sequences, 583 unique, drawn from
twelve generated mutator/region combinations (3 in r1, 6 in r2, 3 in r3) plus custom variants.
*Caveat to record honestly (keep verbatim — it is the load-bearing PoolParty contrast):* every
available mutation type is an **exhaustive, single-event** scan; `snv` emits all three non-reference
bases at every eligible position, one event per oligo. There is **no sampled or random
mutagenesis mode**, and **no pairwise/higher-order mutator** to mix in; combinatorial
variants can only enter via a hand-authored VCF. A no-op (no pattern or custom mutation) oligo can be
added via `--include-no-op-oligo` (CHANGELOG 4.0.0) as a **single** sequence, not replicates — and
only when a background or PAM edit applies to the targeton, in which case what it carries is that
edited reference rather than true WT (`meta_table.py`; CLI help: "Include no-op oligonucleotide if any
background or PAM protection variant is applied to the targeton").
*Source:* `examples/sge/brca1/parameter_input_files/brca1_nuc_targeton_input.txt`;
`examples/sge/brca1_prime_editing/input/input_a.txt`, `input_b.txt`; paper Tables 1–2; CHANGELOG 4.0.0.

**`combinatorial_motif_place` — no.** *(evidence rewritten — the extraction contained a false sentence)*
VaLiAnT has **no motif or element concept and no insertion-scan machinery**. Fixed sequence can enter
an oligo two ways, and the extraction's claim that adaptors are the only route was **wrong**:
1. `--adaptor-5` / `--adaptor-3`, appended once to **every** oligo in the run — "targetons need to be
   processed individually as this function **appends to all sequences in the library**" (§2.4);
2. as an **insertion in a user-supplied custom-variant VCF**. README §*Custom variants*: "Only simple
   variants such as the following are supported:" — then a list: substitutions, **insertions**,
   deletions, indels.
   Any sequence — including a TF-binding motif — can therefore be inserted at any position inside a
   targeton, and VaLiAnT will build, name, annotate and emit one oligo per such insertion.

The value remains **no** because the *combinatorial* half is genuinely absent: there is no motif ×
position enumeration, no orientation or permutation control, no spacing/copy-number sweep, no
element library concept — **the user must author every insertion by hand in a VCF**. Confirmed against
the full recursive tree of `develop` (461 paths): no motif, element, or insertion-scan module anywhere
in `src/valiant/`.
*Source:* README §*Custom variants* and §*Adaptor sequences*; paper §2.4; recursive git tree of
`develop`; `src/valiant/` module listing (87 files).

**`barcode_generation` — no.**
No barcode functionality of any kind. Absent from the 32-field `META_CSV_FIELDS` schema
(`src/valiant/meta_table.py`, re-derived verbatim below), from both CLI option tables, from the JSON
configuration schema, from all 11 wiki pages, and from the full recursive tree. There is **no GC,
edit-distance or Tm machinery** anywhere in the source, and no sequence-uniqueness *design* step —
the only uniqueness handling is post-hoc deduplication of identical oligos into `_unique.csv`.
VaLiAnT's answer to library identity is the variant itself plus the amplicon/adaptor structure, not a
barcode.
*Source:* `src/valiant/meta_table.py`; README metadata + CLI tables; recursive git tree.

**`per_sequence_provenance` — yes.** *(stronger than the extraction stated — provenance is bidirectional)*
Every oligo gets a row with **32** structured columns — in `_meta.csv`, or in `_meta_excluded.csv` if
the length filter rejected it (README, verbatim header):

```
oligo_name, species, assembly, gene_id, transcript_id, src_type, ref_chr, ref_strand, ref_start,
ref_end, revc, ref_seq, pam_seq, vcf_alias, vcf_var_id, mut_position, ref, new, ref_aa, alt_aa,
mut_type, mutator, oligo_length, mseq, mseq_no_adapt, pam_mut_annot, pam_mut_sgrna_id, mave_nt,
mave_nt_ref, vcf_var_in_const, background_variants, background_seq
```

That is provenance well beyond a mutation name: `mutator` (which mutator built it), `vcf_alias` /
`vcf_var_id` (which custom VCF and which record — "To preserve variant provenance, an alias is
assigned to each VCF file", §2.5.2), `pam_mut_annot` / `pam_mut_sgrna_id` (which PAM edits, which
sgRNA), `background_variants` / `background_seq`, `ref_seq` / `pam_seq` (the reference and
PAM-protected reference the oligo was built from), `mave_nt` / `mave_nt_ref` (MAVE-HGVS), and
`vcf_var_in_const` (whether a custom variant landed in a **non-mutable constant region**).
**Provenance is bidirectional:** the output VCFs carry `INFO` tags linking each variant record *back*
to its oligo — `SGE_OLIGO` (= `oligo_name`) and `SGE_SRC` (= `mutator`) on every record, plus
`SGE_VCF_ALIAS` and `SGE_VCF_VAR_ID` on custom-variant records only (README: "The variants can be
linked to the corresponding oligonucleotides via the `SGE_OLIGO` tag, and, for custom variants, to the
original VCF files via the `SGE_VCF_ALIAS` and `SGE_VCF_VAR_ID` tags"). Since v3.0.0 each
library-generating execution also writes a `config.json` recording the full parameter set (present in
the shipped `output_a_exp/config.json`); `--sequences-only` runs exit before that step.
*Source:* `src/valiant/meta_table.py` (`META_CSV_FIELDS`); README §*Oligonucleotide metadata file*,
§*Variant call files*; paper §2.5.2; CHANGELOG 3.0.0; `examples/sge/brca1_prime_editing/output_a_exp/`.

**`automatic_naming` — yes.**
`oligo_name` is generated automatically by the tool; README example:
`ENST00000357654.9.ENSG00000012048.23_chr17:43104102_1del_rc`
(transcript.gene_chr:pos_mutator_rc). Implemented in `src/valiant/meta_table.py` as
`get_sge_oligo_name`, `get_cdna_oligo_name`, `get_sge_oligo_no_op_name`, with the documented
`NO_TRANSCRIPT` placeholder — wiki *Advanced usage*: "When the transcript and gene ID cannot be
determined the oligonucleotide name contains a placeholder string: `NO_TRANSCRIPT`."
Output **filenames** are auto-generated too: in SGE mode from the design — "The file names report
chromosome, coordinates, strand and the sgRNA IDs associated with the targeton" (§2.5.3), e.g. the
shipped `chr17_43115722_43115768_minus_sgRNA_ex2_a_meta.csv`; in cDNA mode they are hash-derived,
`<seq_id>_<md5>_meta.csv`, e.g. `BRCA1_NP_009225.1_sense_0168090b9a169f735c3f685fca7f8063_meta.csv`
— still automatic, just not coordinate-derived.
*Source:* `src/valiant/meta_table.py`; README; wiki *Advanced usage*; paper §2.5.3, Supp. Table 3;
shipped `output_a_exp/` and `examples/cdna/output_exp_v3/` filenames.

**`design_visualization` — no.**
VaLiAnT emits CSV / VCF / JSON only. The complete pinned dependency set is
`charset-normalizer, click, pydantic, pysam` — **no plotting library** — and there is no plotting or
rendering module anywhere in the 461-path recursive tree. The paper's own design figure was made in
third-party software: Fig. 3 legend — "sequence information modified from **Geneious Prime**
(version 2019.04) visualization." The Discussion places visualization downstream and outside the tool:
"VaLiAnT may be combined with downstream analysis software tools ... and to produce annotated
visualizations of variant effect."
*Pre-emptive clause (so a referee cannot say "there are figures in the repo"):* the repo does ship
static PNGs under `examples/images/wiki/` and a hand-built GenBank plasmid model
(`examples/cdna/construct_generation/BRCA1_cdna_pCW57.1_model.gb`) — these are **documentation assets
authored by hand, not tool output**.
*Source:* `pyproject.toml`; recursive git tree; paper Fig. 3 legend, Discussion.

### Block B — assay coverage

**`assay_dms` — yes.** The tool's entire purpose, stated in the title: "an oligonucleotide library
design and annotation tool for **saturation genome editing** and other **deep mutational scanning**
experiments." Both modes are DMS modes (`sge` for endogenous-locus SGE / prime editing, `cdna` for
expression-cassette DMS), and the repo ships four worked BRCA1 DMS/SGE libraries with expected outputs.
*Source:* paper title, abstract, §2.1–§3.4; `examples/`.

**`assay_mpra` — partial.** *(the review attacked this in both directions and could not break it — keep "partial")*
There is **no MPRA-specific support and no MPRA claim**; MPRA appears in the paper only as background
("MAVEs can assess coding and non-coding loci through approaches such as deep mutational scanning
(DMS) and massively parallel reporter assays (MPRAs), respectively"). **However**, non-coding targetons
genuinely work: wiki *Advanced usage* — "UTR regions and non-coding sequence are treated as if
intronic sequence. Targetons that do not overlap any GTF/GFF2 feature won't be associated with a gene
or transcript and their target regions will be considered intronic for the purposes of mutation
generation"; wiki *Input files* — "If CDS-specific mutations are not required, then a GTF/GFF2 file
need not be supplied"; and `sge_proc.proc_contig` has an explicit `if config.gff_fp: ... else:` branch
for the un-annotated case. So `snv` / parametric-deletion saturation of a regulatory element is
achievable and produces usable MPRA-style oligos — **but** there are no barcodes, no reporter-cassette
assembly, no motif insertion/scanning, no promoter/enhancer element grammar, and no MPRA example
anywhere in the docs.
*Source:* paper §1; wiki *Advanced usage*, *Input files*; `src/valiant/sge_proc.py:448`.

**`assay_insilico` — no.** Nothing in the paper, README, wiki, or source relates to designing
libraries to probe sequence-to-function models. There is no model interface, no scorer, no
covariate/design-matrix export for surrogate modelling anywhere in the tree. The single ML mention in
the Discussion runs the other way — DMS data later *training* models ("as predicative models and
machine learning increase in utility ... training datasets rich in biological context will be
important") — not designing libraries to probe a model.
*Source:* paper Discussion; recursive git tree; README; all 11 wiki pages.

### Block C — genomics integration

**`genome_coordinates` — yes.** `valiant sge` is coordinate-driven end to end. "Targeton ranges are
defined by 1-based indexing of genomic coordinates using 'ref_start' and 'ref_end' fields, with
reference chromosome and strand also defined, as 'ref_chr' and 'ref_strand'" (§2.1). Confirmed in the
real shipped targeton TSV (`chr17  -  43115634  43115878  43115726  43115779 ...`). Sequences are
retrieved from a local FASTA + FAI via pysam (`fetch_sequence`, `open_fasta` in `sge_proc.py`).
Species and assembly are **required positional CLI arguments**; "The tool is species and genome-build
agnostic" (§2.3).
*Source:* paper §2.1, §2.3; `examples/sge/brca1/parameter_input_files/brca1_nuc_targeton_input.txt`;
`src/valiant/sge_proc.py`; `examples/sge/brca1/run_brca1_nuc.sh`.

**`transcript_models` — yes (one transcript per execution).**
"To collect gene and transcript information and to apply CDS-specific mutator functions, appropriate
transcript annotation must be provided via a GTF/GFF2 file; only CDS, UTR and stop features are taken
into consideration" (§2.4.2 — paper-era: v4's loader recognises only `CDS` and `UTR` rows and ignores
all other feature types; README and `loaders/gtf.py`). The operative constraint is best anchored on
the **wiki**, which is stricter than the paper's "one transcript per gene": wiki *Input files* —
"**One specific transcript must be supplied for an execution.** GTF/GFF files should not contain
multiple transcripts... CDS and UTR features are required. The last CDS feature gets extended by one codon to include STOP in
computation." The documentation's execution-level rule is **not globally enforced by the code**:
`GtfLoader(contig, strand)` rejects multiple gene/transcript pairs only within that group, while
`run_sge` merely warns when annotation spans multiple contigs, so different transcripts on different
contigs can proceed. `proc_contig` still builds one `Transcript` per contig/strand and carries the
comment `# TODO (future multi-transcript): clear existing exons from the database`. So: real
exon/CDS/UTR-aware transcript handling; the documented rule is one transcript per run, while the code
enforces one per contig/strand — not a full transcript-model engine.
*Source:* wiki *Input files*; paper §2.4.2; `src/valiant/sge_proc.py:448-470`, `:549`;
`src/valiant/transcript.py`, `loaders/gtf.py`.

**`exon_intron_split_codons` — yes.** *(refreshed from v4.0.0 — the extraction cited only the paper)*
Explicitly handled, with a published formula (eqns 1–2) and Supp. Fig. S2: "Retrieval of additional
positions from the reference sequence is necessary to obtain the context of partial codons... **The
extra bases required at either end of the target may come from the same or an adjacent exon**." And:
"SNVs can be introduced into partial, liminal codons in addition to complete codons; the annotation of
the mutations in partial codons is informed by the exonic bases adjacent to the target", while
"CDS-specific mutator functions that result in codon replacement or deletion are only applied to
complete codons within CDS targets (partial, liminal codons are ignored)" (§2.5.1).
**v4.0.0 goes further than the published description** — CHANGELOG 4.0.0, *Added*: "Support for
out-of-frame CDS targeton regions whose 5' and 3' extensions span **both adjacent and distal bases**."
Documented restriction (§2.4.2): "No discrete target region (r1, r2, r3) within a targeton can span
both coding and non-coding sequences; although the complete targeton can be divided into coding and
non-coding regions."
*Source:* paper §2.5.1, §2.4.2, eqns 1–2, Supp. Fig. S2; CHANGELOG 4.0.0;
`src/valiant/cds_seq.py`, `exon.py`, `annot_variant.py`.

**`hgvs_input` — no.** *(source-level verified)*
Every input is a table or a standard genomics file: targeton TSV, reference FASTA + FAI, GTF/GFF2,
PAM-protection VCF, custom-variant VCF manifest CSV, background VCF + optional BED mask,
codon-frequency CSV, JSON config, cDNA multi-FASTA + annotation TSV. Variants enter as VCF rows
(POS/REF/ALT) — "The classification of the variants is based exclusively on the `POS`, `REF`, and
`ALT` fields to be agnostic with respect to the VCF source" (README). Every format loader in
`src/valiant/loaders/` is bed/csv/fasta/gtf/vcf/tsv (the rest of the directory is config, error and
utility modules); **there is no HGVS parser anywhere**.
`src/valiant/mave_hgvs.py` is **write-only** (`get_mave_nt`, consumed by `MetaTable`) — HGVS appears
only on the **output** side (`mave_nt`, `mave_nt_ref`).
*Source:* `src/valiant/loaders/` listing (16 files); `src/valiant/mave_hgvs.py`; README §*Custom
variants*, §*Input files*; wiki *Input files*, *cDNA DMS file formats*.

**`vcf_vep_output` — yes (SGE mode only).** *(refreshed from v4.0.0 — the extraction's single "`_.vcf`" came from paper-era Supp. Table 2)*
Since **v3.0.0** there are **two VCFs per targeton**, not one (CHANGELOG 3.0.0: "Output files:
generate additional VCF file where `REF` does not include PAM variants"). README §*Variant call
files*: "VCF files containing a subset of the metadata in VCF format. The metadata are stored in the
`INFO` field. The `REF` field reports the reference sequence including (`*_pam.vcf`) or excluding
(`*_ref.vcf`) PAM protection edits." Confirmed in shipped expected outputs, e.g.
`chr17_43115722_43115768_minus_sgRNA_ex2_a_pam.vcf` **and** `..._ref.vcf`. Every record carries
`SGE_OLIGO` and `SGE_SRC` `INFO` tags; `SGE_VCF_ALIAS` and `SGE_VCF_VAR_ID` are added for custom
variants only. **All in-range** generated variants are emitted, not just custom ones (oligos rejected
by the length filter appear in neither VCF).
*Caveats to keep:* **cDNA mode emits no VCF** (relative coordinates), and the metadata CSV
deliberately breaks VCF convention — "the output metadata format does not follow VCF convention in
reporting positions, reference and alternative sequences to favour streamlining of downstream
processing" (§2.5.4). **VEP itself is never named by the authors**; the claim is standard-VCF
compatibility, not a VEP integration.
*Source:* README §*Variant call files*; CHANGELOG 3.0.0;
`examples/sge/brca1_prime_editing/output_a_exp/` file listing; paper §2.5.4, Supp. Table 2.

**`consequence_annotation` — yes.** *(strengthened with two v4 facts)*
Codon-level consequence is computed and stored in the metadata fields `ref_aa`, `alt_aa`, `mut_type`.
"SNVs carry extra information when introduced into CDS targets ... namely whether the SNV results in a
**synonymous, missense or nonsense** change" (§2.5.1); the classification is strand-aware ("For the
same codon on different strands, the nucleotide changes are the same but the amino acid changes, and
therefore their classification ... differ") and is precomputed as a codon-indexed table at start-up
(1152 codon/SNV combinations). In-frame deletions are labelled via the `mutator` field.
**v4.0.0 widens this beyond SNVs** — CHANGELOG 4.0.0, *Changed*: "Metadata table: report the type of
mutation for the `aa`, `ala`, and `stop` mutators (`mut_type` field)" and "in cDNA mode, mutations
preserving stop codons are now annotated as nonsense (*vs.* synonymous)". Additionally
`pam_mut_annot` classifies each PAM protection edit as `syn|mis|non|ncd` (README field 26), and
background variants are classified and gated in code — `is_variant_nonsynonymous`,
`is_variant_frame_shifting` (`sge_proc.py`, imported from `targeton.py`), with README: "By default,
errors are raised when background variants are not synonymous or shift the reading frame; the
`force-bg-ns` and `force-bg-indels` flags may be passed to allow them."
*Limits (correct, keep):* no splice-site consequence prediction, no frameshift/NMD calling, no VEP or
other annotator invoked — annotation is codon arithmetic only.
*Source:* paper §2.5.1; CHANGELOG 4.0.0; README metadata table (field 26), §*Background variants*;
`src/valiant/sge_proc.py`, `targeton.py`, `annot_variant.py`, `codon_table.py`.

### Block D — physical construction

**`primer_design` — no.**
VaLiAnT appends **user-supplied** constant sequences (`--adaptor-5`: "DNA sequence to be added at the
5' end of the oligonucleotide"; `--adaptor-3`: the same "at the 3' end") verbatim, and leaves flanking
constant regions for user-chosen primers; it never designs a primer. In the prime-editing example the
P5/P7, Golden Gate Type IIS ends, sgRNA spacer, pegRNA scaffold and PBS are all hand-written into the
adaptor strings by the user (`run_prime_a.sh`: a 121-nt `--adaptor-5` literal). In the cDNA example
the 21 bp at each targeton boundary lie outside r2 and are therefore not subjected to mutators, which
the wiki says "allows for primer binding to amplify and clone specific targetons" — the user picks
them.
No Tm, specificity, or dimer calculation exists anywhere in the tree, and there is no primer3-class
dependency among the four pinned deps.
*Source:* README CLI table; `examples/sge/brca1_prime_editing/run_prime_a.sh`; wiki *cDNA example*;
`pyproject.toml`; recursive git tree.

**`codon_optimization` — partial.** *(deliberately generous toward VaLiAnT — the safe direction)*
VaLiAnT is codon-usage-aware but does not optimise a sequence. It ships a human codon-frequency table
(`src/valiant/data/default_codon_table.csv`; columns triplet / amino acid / frequency / rank
`RANKT, RANK2, ..., RANKU`) and accepts `--codon-table` for other species
(`codon_table.py`, `codon_table_builder.py`, `codon_table_loader.py`, `codon_table_row.py`). Ranks
drive **per-codon substitution**: `ala`/`stop` use the top-ranking triplet for Ala/stop and skip a
codon already equal to it; `inframe` deletes each complete reading-frame triplet; `aa`
"Replace each codon with the top-ranking codon of all amino acids. Given the default codon table, this
results in **19 mutated sequences for each codon mapping to an amino acid** (the reference amino acid
being excluded) **and 20 for each stop codon**" (README); `snvre` swaps a missense codon for the
**top-ranking** synonymous triplet of the SNV's amino-acid product — the second-ranked one if that
triplet is already top-ranking — and expands a synonymous SNV to every other synonymous triplet;
nonsense SNVs follow the analogous stop-codon rule. Reference and triggering triplets and duplicate
products are excluded, so no additional synonymous triplet means no added design (README §*SNVRE*;
`mutators/snv_re.py`). This is what the paper says "allows for insights into the effect of codon
sequence on missense changes" (§2.1).
That is per-codon frequency-ranked substitution — **not** whole-CDS codon optimisation, codon
harmonisation, CAI maximisation, or avoidance of structure/repeats/restriction sites.
*Source:* README §*All amino acid codon substitution*, §*SNVRE*, §*Alanine/stop codon substitution*;
paper §2.1; `src/valiant/codon_table*.py`; `src/valiant/mutators/codon.py`, `deletion.py`, `snv_re.py`;
`src/valiant/data/default_codon_table.csv`; CHANGELOG 4.0.0 (codon-table loading fix).

**`synthesis_constraints` — partial.** *(length only)*
Only oligo length is enforced. `--max-length` (default 300 bp): "Oligonucleotides exceeding a given
length (`max-length` option) will not be included in the unique oligonucleotide files and their
metadata will be stored in separate files marked as 'excluded'." `--min-length` was added in v3.0.0.
The measured length **includes the adaptors**, which are appended after any reverse complement
(`meta_table.py`).
Confirmed in the JSON configuration schema as `minOligoLength` / `maxOligoLength` (shipped NXRL config:
`"minOligoLength": 1, "maxOligoLength": 300`). The Discussion frames it as a synthesis-platform
constraint: "At the time of writing, chemical synthesis at the scale of ~300 bp is possible on a large
scale ... We have provided an optional filter to remove sequences above a user-defined length, which
is configurable to accommodate for future increases in synthesis capability." `_unique.csv` is
described as ready for "Synthesis submission" (Supp. Table 2).
*One adjacent input-level check, mentioned pre-emptively:* README — "**Ambiguous nucleotides are not
allowed in the reference sequence. Soft-masking is ignored.**" That is an input-validity rule, not a
synthesis constraint.
Still **no GC-content, homopolymer, repeat, secondary-structure, or restriction-site checking**
anywhere in the source tree.
*Source:* README CLI tables and §*Reference sequence*; paper Discussion, Supp. Table 2;
`examples/sge/nxrl/chr17_31226400_31226678_plus_sgRNA_1146241047_valiant_config.json`.

### Block E — engineering

**`interface` — yes (CLI only; no Python API, no web service, no GUI).** — **see §8, flagged**
"The tool is implemented as a standalone executable Python package exposing a command line interface"
(§2.3). `pyproject.toml` `[project.scripts]` has **exactly one** entry
(`valiant = "valiant.__main__:main"`) and **no `[project.entry-points]`**. Two subcommands (`sge`,
`cdna`) plus a whole-run JSON configuration entry point (`valiant -c config.json`), which the NXRL
example demonstrates. Distributed as source + Dockerfile + prebuilt `quay.io/wtsicgp/valiant`, with
**Singularity as a first-class documented path** (its own README section with a full
`singularity exec --cleanenv -B ...` invocation).
**No documented Python API** — the package is importable but there is no `docs/` directory in the
tree, no API reference in the 875-line README, and no API page among the 11 wiki pages. **No web
service, no GUI. Not on PyPI, not on bioconda.**
*Column-semantics note:* the value `yes` is used here in the same sense as the `seqpro` record's
`yes (Python API only)` — i.e. "ships a documented, usable interface" — and the cell **must carry the
"(CLI only)" qualifier**. If the survey instead defines this column as "exposes a programmatic /
library API", the honest value for VaLiAnT is **no**. See §8.
*Source:* paper §2.3; `pyproject.toml`; README §*Docker image*, §*Singularity image*;
`examples/sge/nxrl/run.sh`; recursive git tree; pypi/anaconda API checks.

**`license` — yes (AGPL-3.0-or-later).**
Confirmed three ways: `gh api /repos/cancerit/VaLiAnT` → `"license": "AGPL-3.0"`; README LICENCE block
and per-file headers ("GNU Affero General Public License ... version 3 of the License, or (at your
option) any later version"); paper Availability — "VaLiAnT is licensed under AGPLv3."
*Presentation warning for the table:* put the **SPDX identifier in the cell**, not a bare `yes`. AGPL
is strong copyleft with a network-use clause, materially unlike the MIT licences of `seqpro` and
`tangermeme` in the same survey; a single `yes` column would erase a real difference.
*Source:* GitHub API; README; source-file headers; paper Availability.

**`maintained` — partial.** — **VALUE CHANGED (was `yes`)**
Dormant but alive. Re-verified 2026-08-10:
- Last commit on `develop`: **2024-04-22T13:51:54Z**, "Merge tag '4.0.0' into develop" with body
  "Background variant support" (Luca Barbon). `pushed_at = 2024-04-22T13:52:12Z`. **~27.6 months** with zero
  commits at the assessment date.
- Tags 1.0.0 / 2.0.0 / 3.0.0 / 3.0.1 / 4.0.0, but **zero GitHub Release objects**
  (`gh api .../releases` → length 0).
- Docker images `latest` and `4.0.0` both last pushed **2024-04-22**; no newer tag on quay.io.
- **Both open issues are automated Snyk dependency-bot PRs** — #10 (2024-03-03) and #11 (2024-04-03),
  both "[Snyk] Fix for 5 vulnerabilities", user `superjw` — left unattended. External issue traffic is
  minimal but not zero: #8 (2023-06-27, `delagee`, `author_association: NONE`) is a substantive human
  bug report on MAVE-HGVS strings, answered and closed by a maintainer on 2023-07-18 — the day 3.0.1
  shipped with a MAVE-HGVS fix. Every other issue was authored by the maintainers.
- Not archived; 6 stars, 3 forks; documentation and examples with expected outputs present.

`yes` claimed active maintenance the evidence does not show, so this record downgrades to `partial`
("dormant since 2024-04-22; not archived; source and container images remain published").
**Note for the survey:** this correction *flatters nothing* — it moves VaLiAnT down, so it carries no
referee-attack risk, but the survey should **state its maintenance criterion** (suggested: "commit
within 24 months of the assessment date = yes; older but with published artifacts and not archived =
partial; archived / no published artifacts = no") so the column is defensible across all eleven tools.
*Source:* GitHub repo/commits/releases/issues API; quay.io tag API — all re-queried 2026-08-10.

---

## 5. VaLiAnT's own documented example libraries

All under https://github.com/cancerit/VaLiAnT/tree/develop/examples, each with `run*.sh` and
`check*.sh` and a committed expected-output directory. **Validation is by `md5sum` against that
directory** ("The validation scripts depend on the `md5sum` utility (with `md5` as a fallback)",
`examples/README.md`) — a genuinely strong reproducibility claim. The shared SGE reference and the
prime-editing inputs must be **unpacked first**
(`sge/unpack_reference.sh`; `brca1_prime_editing/unpack_inputs.sh`).

1. **BRCA1 nucleotide SGE library (`brca1_nuc`)** — exons 2–5 of BRCA1 (ENST00000357654.9 /
   NM_007294.4, GRCh38), **4 targeton rows in one file, one invocation**, 245–251 bp with 25/25 bp (or
   20/41 bp for exon 5) intronic extensions. **Corrected action vectors, quoted verbatim from
   `brca1_nuc_targeton_input.txt`:** exons 2–4 use
   `"(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"`, but **exon 5 uses
   `"(1del,2del0,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"`** — its r1 is 20 bp (even), so
   the offset-1 two-base deletion pattern is replaced by the offset-0 one. Custom variants from
   ClinVar (2020-11-07) and gnomAD v3.0 via a manifest CSV; PAM/protospacer protection edits at one
   sgRNA per targeton; P5/P7 adaptors; `--revcomp-minus-strand`.
   **Corrected flags:** `run_brca1_nuc.sh` passes **only** `--pam`, `--vcf`, `--revcomp-minus-strand`,
   `--adaptor-5`, `--adaptor-3`, `--gff`. It does **not** pass `--max-length 300` — 300 bp is simply
   the default. (The paper describes the option as used; the shipped script does not set it.)
   Exon 2 alone: 1000 total sequences, 583 unique (paper Table 2). Unique counts exons 2–5:
   583 / 740 / 825 / 1185 (Supp. Table 3). (The committed v4.0.0 outputs instead contain 1000
   in-range plus 1 excluded exon-2 metadata row and `_unique.csv` counts of 583 / 740 / 825 / 1217
   for nucleotide and 340 / 500 / 580 / 919 for peptide.)
2. **BRCA1 peptide/amino-acid SGE library (`brca1_pep`)** — same four targetons; action vector
   verbatim `"(),(aa,inframe),()"` for all four rows (r1/r3 empty). Same flags as `brca1_nuc` minus
   `--vcf`. Unique counts exons 2–5: 340 / 500 / 580 / 918 (Supp. Table 3).
3. **BRCA1 saturation prime editing library (`brca1_prime_editing`)** — two pegRNAs over BRCA1 exon 2;
   the pegRNA RTT *is* the targeton. `input_a.txt`: chr17:43115722–43115768, minus strand, ext `"4,0"`,
   action vector `"(1del,snv,2del0),(1del,snvre,inframe,stop,ala,aa),()"`; `input_b.txt`:
   chr17:43115735–43115783, ext `"0,4"`, `"(),(1del,snvre,inframe,stop,ala,aa),(1del,snv,2del0)"`.
   `--max-length 250` **is** passed here; `--adaptor-5` is a 121-nt literal (P5 + Golden-Gate Type IIS
   end + sgRNA spacer + pegRNA scaffold) and `--adaptor-3` a 38-nt literal (PBS + Type IIS end + P7),
   both hand-written by the user. ClinVar + gnomAD custom variants. 454 / 483 unique sequences.
   Docs: wiki *Saturation-prime-editing-example*.
4. **BRCA1 cDNA DMS library (`examples/cdna`)** — the 5592 bp BRCA1 CDS (NP_009225.1, 22 exons)
   inserted *in silico* into pCW57.1 (Addgene #41393) downstream of a TRE promoter (SalI/NheI digest —
   both the paper and the cDNA wiki page spell the second enzyme `NehI`).
   `cdna_targeton.tsv` has **exactly 40 data rows**, targetons 123–237 bp (the paper says 132–237;
   the shortest shipped row, 901–1023, is 123 bp), r2 ranges 81–195 bp, tiled
   with overlap; **every row uses the same action vector `"snv,1del,snvre,ala,stop,inframe,aa"`**;
   `cdna_annot.tsv` supplies `cds_start=331, cds_end=5922`. Adaptors P5/P7; no `--max-length` passed.
   858–2092 unique oligos per targeton; **62,754 variants** in the "final concatenated library"
   (§3.2). The directory also ships **`output_exp_v3/`** (v3-era expected outputs, filenames
   `<seq_id>_<md5>_meta.csv`) alongside `output_exp/`, and **`construct_generation/`** containing
   `BRCA1_cdna_pCW57.1_model.gb` (hand-built GenBank plasmid model), `BRCA1_NP_009225_1_CDS.fa` and
   `cdna_targeton_plan.txt` — direct evidence that **the 40-targeton tiling plan was authored by hand**,
   reinforcing that tiling is a documented manual recipe, not a function.
5. **NXRL background-variant example (`examples/sge/nxrl`, repo only, not in the paper)** — targeton
   `chr17_31226400_31226678_plus_sgRNA_1146241047`, transcript ENST00000358273.9. **Corrected flags,
   read from the shipped
   `chr17_31226400_31226678_plus_sgRNA_1146241047_valiant_config.json`:** it sets
   `"backgroundVCFFilePath": "input/background_GATK_NXRL_variants.vcf"`,
   `"forceBackgroundNonSynonymous": true`, `"forceBackgroundFrameShifting": true`,
   `"reverseComplementOnMinusStrand": false`, `minOligoLength: 1`, `maxOligoLength: 300`, plus PAM VCF
   and custom-variant manifest — and **`maskBackgroundFilePath` is absent entirely (no BED mask is
   used)**, contrary to the original extraction. It is invoked as
   `valiant -c chr17_31226400_31226678_plus_sgRNA_1146241047_valiant_config.json`, making it the
   repo's **only worked demonstration of the JSON-config entry point**.
6. **Docker usage example** — wiki *Docker-usage-example*, mounting inputs into
   `quay.io/wtsicgp/valiant`. **Stale:** the page still pins `VERSION=2.0.0`, three tagged releases
   behind (`3.0.0`, `3.0.1`, `4.0.0`).

---

## 6. Notable capabilities not covered by the row list

1. **The deletion mutator is FULLY PARAMETRIC in v4.0.0** — the single most under-described capability
   in the original extraction, and exactly what an author-referee would notice. It is **not** a fixed
   `1del`/`2del0`/`2del1` triple. README §*Parametric deletion*: "Non-overlapping stretches of
   nucleotides of a given length are deleted starting from a given offset. No partial deletions are
   performed at the end of the target regions. Format: `<SPAN>del[<OFFSET>]` (the offset is assumed to
   be zero if not set)." CHANGELOG 4.0.0: "Mutator: deletion pattern span and offset can be set to
   **any valid value**." Confirmed in `src/valiant/int_pattern_builder.py`
   (`IntPatternBuilder(offset, span)`, `span > 0` the only constraint) and `mutator_type.py`
   (`PARAMETRIC_MUTATOR_TYPES = {MutatorType.DEL}`). So `3del1`, `6del0`, `10del2` … all work —
   **arbitrary tiled deletion scans**. (For backwards compatibility the metadata reports `1del0` as
   `1del`.) Changes no cell value.
2. **PAM / protospacer protection edits** — user-defined SNVs keyed to an
   sgRNA ID via a VCF `INFO` tag, applied to **every** oligo of a tagged targeton so the HDR template
   is refractory to Cas9 re-cutting. The paper describes them as synonymous or non-coding, but that is
   design intent, not an input restriction: v4's parser (`pam_variant.py`) requires only a one-base
   REF/ALT and an `SGRNA` tag, which is why `pam_mut_annot` can report `mis` and `non` as well.
   Multiple edits per sgRNA and multiple sgRNAs per targeton.
   Recorded per oligo as `pam_mut_annot` (`syn|mis|non|ncd`) and `pam_mut_sgrna_id`. Supp. Table S1
   gives the authors' selection heuristics — used by hand, **not automated by the tool**.
3. **Saturation prime editing (pegRNA RTT) design** — treats the RTT as the targeton, with
   `--revcomp-minus-strand` for strand-specific RTTs and adaptors carrying spacer/scaffold/PBS.
4. **Background variants (v4.0.0)** with **computed, not user-specified, scoping.** README: background
   variants "are applied in the minimal range of positions that spans at least the entire CDS, further
   extended to the boundaries of any targeton overlapping the CDS, and finally to the start and end
   position of the first and last background variant intersecting the resulting range" — with an exon
   lift-over over the shifted coordinate frame (`GenomicPositionOffsets`, `transcript.lift_exons(gpo)`
   in `sge_proc.py`). Non-trivial coordinate machinery.
5. **Automatic conflict detection between generated mutations and background variants.** README:
   "Mutations that overlap background variants are discarded; warnings are always raised to identify
   them, with their positions expressed in absolute genomic coordinates for custom variants, and
   targeton-relative coordinates for pattern variants." In v4.0.0 this discard check is limited in
   code to coordinate-shifting background variants (indels); coincident background SNVs are not
   masked by it (`GenomicPositionOffsets`). Plus hard errors on non-synonymous or frame-shifting
   background variants unless `--force-bg-ns` / `--force-bg-indels` are passed.
6. **Two output VCFs per targeton with oligo back-links** (v3.0.0+): `*_pam.vcf` (REF **includes** PAM
   edits) and `*_ref.vcf` (REF **excludes** them), every record carrying `SGE_OLIGO` and `SGE_SRC`
   `INFO` tags, and custom-variant records additionally `SGE_VCF_ALIAS` and `SGE_VCF_VAR_ID` —
   provenance is bidirectional.
7. **Order-invariant deduplication as an explicit reproducibility guarantee.** README: "When multiple
   oligonucleotides have the same sequence, the first name in **lexicographic order** is chosen";
   CHANGELOG 2.0.0: "Make unique oligonucleotide name selection order-invariant."
8. **In-memory SQLite is the generation engine** (`db.py`, `data/ddl.sql`, `queries.py`, `sql_gen.py`)
   — relevant context for any scalability comparison, and what makes `lazy_evaluation = no` directly
   verifiable rather than inferred. Note this replaced the paper-era pandas/PyRanges stack.
9. **`vcf_var_in_const` metadata field** — flags custom variants that landed in a constant
   (non-mutable) region, confirming the paper's "custom variants are still incorporated into constant
   regions."
10. **Ingesting real variant catalogues** — ClinVar / gnomAD VCFs turned into oligos alongside the
    systematic scans, with per-file aliases and original IDs preserved; the Discussion offers this as a
    standalone use ("conversion of VCF annotation files to oligonucleotide sequences").
11. **MAVE-HGVS output strings** (`mave_nt`, `mave_nt_ref`) for MaveDB-style interoperability —
    linear-genomic, targeton-relative, reference-free; the README documents the known position-zero
    edge case for liminal insertions.
12. **Strand handling** — `--revcomp-minus-strand` emits reverse-complemented oligos for minus-strand
    targets / ssODN-orientation-specific SGE; codon consequence classification is strand-aware.
13. **`aa` mutator quantified** — 19 substitutions per amino-acid codon, 20 per stop codon (README),
    useful for library-size comparisons against PoolParty.
14. **cDNA mode multiplexes** — the multi-FASTA can hold several cDNAs (`seq_id`), each with several
    targeton rows, in one run.
15. **Singularity is a first-class documented deployment path** — its own README section with a full
    `singularity exec --cleanenv -B ...` invocation, not merely a Windows workaround.
16. **Deduplication + synthesis-ready export** (`_unique.csv` for "Synthesis submission") and
    `_meta_excluded.csv` for length-filtered sequences.
17. **Reference-retrieval QC** — `ref_sequences.csv` plus `--sequences-only`, which writes the QC file
    and exits before generating anything.
18. **JSON run configuration** — `valiant -c config.json` as an entry point, and a `config.json`
    written per library-generating execution recording the full parameter set (v3.0.0+).
19. **Species/assembly agnostic** with pluggable codon-frequency tables (`--codon-table`).
20. **Targeton tiling recipe** for exons longer than one oligo (wiki *Advanced usage*), including
    in-frame overlapping tiles — a documented **manual pattern**, not an automated tiling function
    (corroborated by the hand-authored `cdna_targeton_plan.txt`).

---

## 7. Stated limitations (from the paper, README and wiki — verbatim or near-verbatim)

- **One transcript per execution (documented, not globally enforced).** Paper: "One transcript per
  gene is allowed to avoid ambiguities in matching target regions." Wiki (stricter): "One specific
  transcript must be supplied for an execution. GTF/GFF files should not contain multiple
  transcripts." The code instead enforces one gene/transcript pair per contig/strand and only warns
  when annotation spans multiple contigs; it carries `# TODO (future multi-transcript)`.
- "No discrete target region (r1, r2, r3) within a targeton can span both coding and non-coding
  sequences; although the complete targeton can be divided into coding and non-coding regions." A
  region overlapping **more than one CDS exon** is likewise rejected
  (`targeton.py: InvalidTargetonRegion("overlaps multiple exons")`).
- CDS-specific mutators are rejected on non-CDS targetons: `CRITICAL:root:Invalid mutator 'snvre' for
  targeton!` (wiki *Advanced usage*).
- Partial / liminal codons are ignored by codon-replacement mutators; SNVs, parametric deletions and
  custom variants can still reach them.
- **cDNA mode:** "All mutator functions except 'custom vcf' are currently supported" (§2.2), and no
  VCF output. `valiant cdna` also exposes no `--pam`, `--bg`/`--bg-mask`, `--revcomp-minus-strand` or
  `--sequences-only` (`cdna_cli.py` takes only `--annot` plus the common options).
- **Only "simple" custom variants:** README §*Custom variants* — "Only simple variants such as the
  following are supported:" (substitutions, insertions, deletions, indels), classified from
  `POS`/`REF`/`ALT` alone. v4 also reads **only the first ALT** of a multiallelic record
  (`custom_variant.py`, `r.alts[0]`, with a `# TODO: handle multiple ALT's?`) and skips monomorphic
  records, so a catalogue VCF can lose alleles silently apart from a log line.
- **Reference input rule:** "Ambiguous nucleotides are not allowed in the reference sequence.
  Soft-masking is ignored."
- Oligos above `--max-length` (default 300 bp) are dropped from the synthesis file into
  `_meta_excluded.csv`; the authors tie this to synthesis-platform limits.
- **Adaptors are appended to *all* sequences in a run**, so per-targeton adaptor differences require
  separate runs: "targetons need to be processed individually as this function appends to all
  sequences in the library."
- The README states that mutations overlapping background variants are **discarded** and warnings are
  always raised; v4 code applies this discard test only to coordinate-shifting background variants
  (indels), not coincident background SNVs.
- **Extension requires editing the source** — no plugin system: "We envision that, as opensource
  software, VaLiAnT will be further modified by the community. For example, by the addition of new
  mutator actions."
- Upstream design choices are the user's: "There is also scope for VaLiAnT to be combined with other
  software, such as upstream heuristic functions to select appropriate input information." Targeton
  ranges, sgRNA choice and protection-edit selection are all done by hand (Supp. Table S1 documents the
  authors' heuristics, which the tool does not implement).
- Not installable natively on Windows (WSL / Docker / Singularity required).

---

## 8. Unresolved / flagged

**`interface` — flagged, value kept at `yes` with a mandatory qualifier.**
This is **not a factual dispute** — extractor and reviewer agree on every fact (one console script,
`valiant = valiant.__main__:main`; no `[project.entry-points]`; no `docs/` directory; no API page in
the wiki; not on PyPI or bioconda; Docker + Singularity images). It is a **column-semantics** problem:

- If the column means *"ships a documented, usable interface"* → **yes**, and the cell must read
  **`yes (CLI only)`**. This is the reading used by the `seqpro` final record, whose cell reads
  `yes (Python API only)`; consistency across the survey requires the same reading here.
- If the column means *"exposes a programmatic / library API"* — the natural reading in a table whose
  thesis is that PoolParty is a composable Python library — the honest value for VaLiAnT is **no**.

**Required action for the table builder:** state the column definition in the table caption, and never
print a bare `yes` for VaLiAnT. A bare `yes` sitting beside evidence that says "no Python API" is an
internal contradiction a referee will quote back. Recommended cell text:
`CLI only (Docker/Singularity; no Python API)`.

**Nothing else is unresolved.** No capability value in this record is set to `unknown`: the reviewer
attacked every "no" against the current v4.0.0 source and none of them broke, and the two "partial"
judgement calls (`library_as_object`, `assay_mpra`) survived attack from both directions.

---

## 9. Availability today (re-checked 2026-08-10)

- **Repo alive, public, not archived:** https://github.com/cancerit/VaLiAnT — default branch
  `develop`, 6 stars, 3 forks, 2 open issues (both Snyk bot PRs).
- **Last commit 2024-04-22** (`pushed_at 2024-04-22T13:52:12Z`); latest tag **4.0.0**, same date.
  **Zero GitHub Release objects** — tags only. ~27.6 months dormant.
- **Container images remain published:** `quay.io/wtsicgp/valiant:4.0.0` and `:latest`, both pushed
  2024-04-22 (older tags: 3.0.1 = 2023-07-19, 3.0.0 = 2022-10-12, 2.0.0 = 2021-07-12,
  1.0.0 = 2020-12-21). Singularity documented as a first-class path.
- **Documented source-install path:** `python3.11 -m venv .env && source .env/bin/activate && pip install .`
  Requires Python ≥ 3.11; Linux/macOS native, Windows via WSL / Docker / Singularity. Dependencies are
  fully pinned: `charset-normalizer==3.3.2, click==8.1.7, pydantic==2.7.0, pysam==0.22.0`.
- **Not on PyPI** — `pypi.org/project/valiant` is Duncan Dickinson's unrelated dependency-audit tool
  (v0.2.3, 2021). **Not on bioconda** — `api.anaconda.org/package/bioconda/valiant` returns
  `{"error": "\"valiant\" could not be found"}`. **No web service.**
- **Documentation:** 875-line README plus an 11-page GitHub wiki; examples with md5-validated
  expected outputs. One stale page: *Docker-usage-example* still pins `VERSION=2.0.0`.
- **Verdict: source and container artifacts remain published, but the project is dormant** — no
  commits in ~2.3 years, almost no external issue traffic (one substantive report, #8 in 2023,
  answered), unattended security-bot PRs.

## 10. Confidence notes

- **Highest confidence:** seven of the eight "no" cells (`dag_chaining`, `combinatorial_motif_place`,
  `barcode_generation`, `design_visualization`, `assay_insilico`, `hgvs_input`, `primer_design`; the
  eighth, `lazy_evaluation`, is covered in the next bullet).
  Each was checked against the paper end to end, the 875-line README, the CHANGELOG, `pyproject.toml`,
  all 11 wiki pages, the shipped example scripts and input files, and the **full 461-path recursive git
  tree** — plus direct reads of `mutator_type.py`, `int_pattern_builder.py`, `meta_table.py`, `db.py`
  and all 568 lines of `sge_proc.py`. These are verified absences, not failures to find.
- **`lazy_evaluation` (no)** — upgraded from *inferred* to *code-verified* in this pass
  (`sqlite3.connect(':memory:')` → full materialisation → `MetaTable.to_csv`).
- **`library_as_object` (partial)** — the one genuine judgement call. The factual anchors are solid on
  both sides (one declarative file drives many targetons across many contigs in one run; outputs are
  per-targeton files the user concatenates, in the paper's own words). The *label* is arguable; the
  facts are not.
- **`assay_mpra` (partial)** — attacked from both directions and survived both. Non-coding saturation
  genuinely works (wiki + the `else:` branch in `proc_contig`); there is genuinely no MPRA machinery,
  claim, or example. "Partial" is the only label that survives.
- **`interface` (yes, CLI only)** — see §8. Semantics, not facts.
- **`maintained` (partial)** — the only value changed in this pass; the survey should publish its
  maintenance criterion so the column is consistent across all eleven tools.
- **`vcf_vep_output` (yes)** — VEP is **never named** by the authors; the claim is standard-VCF output
  (two files per targeton in SGE mode), not a VEP integration. cDNA mode emits no VCF at all.
- **`codon_optimization` / `synthesis_constraints` (partial)** — both are deliberately generous
  "partials" that err toward VaLiAnT (the safe direction). The underlying facts (rank-driven per-codon
  substitution; min/max length filter only) are certain.
- **Version drift is the main residual risk in any comparison table.** Everything about background
  variants, parametric deletions, two output VCFs, MAVE-HGVS, `--include-no-op-oligo`, `mut_type` for
  `aa`/`ala`/`stop`, and the JSON config comes from the repo/CHANGELOG only — it **postdates the 2022
  paper**. The paper's §2.3 already misdescribes v4 (pandas and PyRanges are no longer dependencies).
  The table must be explicit that it describes **VaLiAnT 4.0.0**.
