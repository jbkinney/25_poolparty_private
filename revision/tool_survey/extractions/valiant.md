# VaLiAnT — evidence memo

**Full name:** Variant Library Annotation Tool
**Authors / institution:** Barbon, Offord, Radford, Butler, Gerety, Adams, Tan, Waters — Wellcome Sanger Institute
**Paper:** Barbon L. et al., *Bioinformatics* 38(4):892–899 (2022), doi:10.1093/bioinformatics/btab776 (Advance Access 16 Nov 2021)
**Repo:** https://github.com/cancerit/VaLiAnT (default branch `develop`)
**Wiki / manual:** https://github.com/cancerit/VaLiAnT/wiki
**Container:** quay.io/wtsicgp/valiant
**License:** AGPL-3.0-or-later

---

## 1. Sources consulted

| Kind | Ref |
|---|---|
| pdf | `papers/Barbon2022_VaLiAnT_all.pdf` (11 pp incl. supplementary), text extracted with PyMuPDF |
| prior_analysis | `prior_analyses/Barbon2022_VaLiAnT_all_analysis.md` |
| repo | https://github.com/cancerit/VaLiAnT (README on `develop`) |
| repo | https://raw.githubusercontent.com/cancerit/VaLiAnT/develop/CHANGELOG.md |
| repo | https://raw.githubusercontent.com/cancerit/VaLiAnT/develop/pyproject.toml |
| repo | GitHub API: `/repos/cancerit/VaLiAnT`, `/tags`, `/releases`, `/commits?sha=develop`, `/contents/src/valiant`, `/contents/src/valiant/mutators`, `/contents/examples/sge/nxrl` |
| docs | Wiki raw markdown: `Input-files.md`, `Output-files.md`, `Advanced-usage.md`, `cDNA-example.md`, `cDNA-DMS-file-formats.md`, `Saturation-prime-editing-example.md` |
| docs | https://raw.githubusercontent.com/cancerit/VaLiAnT/develop/examples/README.md |
| pypi | https://pypi.org/pypi/valiant/json — **NOT this tool** (an unrelated "audit tool to help investigate Python dependencies" by Duncan Dickinson, v0.2.3, 2021-04-04). VaLiAnT is **not on PyPI**. |
| web | https://quay.io/api/v1/repository/wtsicgp/valiant/tag/ |

### Corrections to the prior analysis
The prior notes are broadly right about positioning but contain errors of fact:
- Prior note says mutators include `'ala'` etc. and lists `snv/1del/2del0/2del1/snvre/ala/stop/inframe` — correct — but **omits `aa`** (all-amino-acid codon substitution), which is the protein-level saturation mutator and is a real, documented feature.
- Prior note says "MPRA support: No". Strictly, VaLiAnT *can* saturate non-coding/intergenic targetons with basic mutators (no GTF needed), which is an MPRA-adjacent design; it just has no MPRA-specific machinery. Recorded here as **partial**, not a flat no.
- Prior note does not mention **background variants** (`--bg`, added v4.0.0, 2024), **MAVE-HGVS output**, the **JSON config file**, or `--include-no-op-oligo`. These postdate the paper and only exist in the repo.
- Prior note's "193 citations" and Semantic Scholar claims were not re-verified here.

---

## 2. What VaLiAnT is

Command-line Python tool. Two subcommands:
- `valiant sge` — genomic reference sequences, **absolute (genome) coordinates**, up to five regions per targeton (`c1-r1-r2-r3-c2`).
- `valiant cdna` — user-supplied multi-FASTA, **relative coordinates**, single target region (r2).

> "VaLiAnT is run from the command line. Inputs, outputs and broad processes are summarized in Figure 1" (paper §2.1)

> "SGE and cDNA DMS functions are mapped to separate subcommands: 'sge' and 'cdna'. Library type is listed in the source field ('src_type') in metadata output files, this can aid coordinate interpretation, which is absolute (i.e. genomic location) for SGE/saturation prime editing and relative (i.e. reference sequence position) for cDNA." (paper §2.5)

### Mutator functions (README, `develop`)
1. `1del` / `2del0` / `2del1` — parametric deletions (frame-agnostic)
2. `snv` — all single-nucleotide substitutions at every position
3. `inframe` — codon (triplet) deletions, CDS only
4. `ala` — alanine scan, CDS only
5. `stop` — stop-codon substitution scan, CDS only
6. `aa` — all-amino-acid codon substitution, CDS only
7. `snvre` — `snv` plus synonymous-codon expansion and alternative-codon redundant oligos, CDS only
8. custom variants imported from VCF

Source tree confirms the mutator set is closed and small: `src/valiant/mutators/` contains only `codon.py, deletion.py, snv.py, snv_re.py, utils.py`.

### Fixed processing pipeline (not user-composable)
background variants → PAM/protospacer protection edits → mutator functions on r1/r2/r3 → custom VCF variants (incl. constant regions) → adaptor-5/adaptor-3 appended → max-length filter → dedup to `_unique.csv`.

> "The final oligonucleotide sequences are then assembled from: invariant region sequences (optional), adaptor sequences (optional) and target region mutated sequences." (paper §2.5.3)

---

## 3. Capability-by-capability evidence

### BLOCK A — library specification

**library_as_object = partial.** The whole library is declared *declaratively* in one BED-like targeton TSV — one row per targeton, with an `action_vector` naming the mutators for each of r1/r2/r3 — so the user does not write loops. But there is no library object a user can hold, transform, or pass around; outputs are **per-targeton files** that the user concatenates. Wiki *Output files* / Supp. Table 2: `_meta.csv`, `_meta_excluded.csv`, `_unique.csv`, `_.vcf` are all "Targeton-specific"; only `ref_sequences.csv` is "Execution-specific". The cDNA example designs "40 targetons ... The total number of variants to cover all of the CDS and making up the **final concatenated library** is 62 754" (paper §3.2) — i.e. concatenation is the user's job. Source: paper §2.4.1, §3.2; wiki Input-files.md, Output-files.md.

**dag_chaining = no.** There is no mechanism to chain, nest, or compose design steps. Mutators are applied independently to disjoint sub-regions of a single reference; the order of operations (background → PAM → mutator → custom VCF → adaptors) is hard-coded in the tool, not user-specified. Quote: "Regions 1-3 (r1-3) can be changed by mutator functions—**each independently of each other**—by detailing corresponding mutator lists in the BED-like input file." (paper §2.4.1). The only "composition" is that `snvre` internally calls `snv`. Source: paper §2.4.1, §2.5.3; README mutator table; `src/valiant/mutators/` listing.

**lazy_evaluation = no.** Batch CLI: one invocation materialises every oligo, writes the complete metadata CSV, unique CSV, and VCF, and reports full library counts (Table 2 of the paper enumerates totals per mutator per region). There is no generator/streaming/on-demand interface, and no documented Python API through which to request a subset. Internally the tool builds full tables (`src/valiant/db.py`, `queries.py`, `sql_gen.py`, `meta_table.py`). Source: paper §2.5.3, Table 2, Supp. Table 3; repo source listing.

**mixed_mutagenesis_one_pool = yes.** This is one of VaLiAnT's headline features. A single targeton row carries a three-part action vector so different mutator types run on different sub-regions, and multiple mutators can be listed per region. BRCA1 exon 2 (paper Table 1): "mutators for r1: `2del1, snv, 1del`; mutators for r2: `snvre, inframe, ala, stop, 1del`; mutators for r3: `2del0, snv, 1del`" — plus ClinVar and gnomAD custom variants incorporated across the entire targeton, all in one output pool (Table 2 shows the pooled composition: 12+12+104+17+312+140+17+17+42+327 = 1000 sequences, 583 unique).
*Caveat to record honestly:* every available mutation type is an **exhaustive, single-event** scan. There is no sampled/random mutagenesis mode and no pairwise/higher-order mutator to mix in; multi-base variants can only enter via user-supplied VCF. WT ("no-op") oligo inclusion exists only as a single sequence via `--include-no-op-oligo` (CHANGELOG v4.0.0), not as replicates.

**combinatorial_motif_place = no.** VaLiAnT has no concept of a motif or a placeable element. The only insertable fixed sequences are `--adaptor-5` / `--adaptor-3`, which are appended once to *every* oligo in the run: "generic adapters for universal amplification of all oligonucleotides in a synthesized pool can be specified on the command line and are appended after mutator function actions ... **this function appends to all sequences in the library**" (paper §2.4/§2.1). No positional, orientational, or permutational placement. Searched README, all 11 wiki pages, and the full `src/valiant/` module listing — no motif/element/insertion-scan module exists.

**barcode_generation = no.** No barcode functionality anywhere: not in the README feature list, not in the CLI option table, not in the 32-column metadata schema, not in the source module listing. Design constraints such as GC content and edit distance are never mentioned. (VaLiAnT's answer to library identity is the variant itself plus adaptor/amplicon structure, not barcodes.)

**per_sequence_provenance = yes.** Every oligo gets a row in `_meta.csv` with 32 structured columns (README "Oligonucleotide metadata file"), including provenance beyond the mutation name: `mutator` ("Label of the mutator type that generated the oligonucleotide"), `vcf_alias` and `vcf_var_id` (which custom VCF a variant came from — "To preserve variant provenance, an alias is assigned to each VCF file", paper §2.5.2), `pam_mut_annot` and `pam_mut_sgrna_id` (which PAM-protection edits and which sgRNA were applied), `background_variants` and `background_seq`, `ref_seq` / `pam_seq` (the reference and PAM-protected reference the oligo was built from), `src_type`, `species`, `assembly`, `gene_id`, `transcript_id`, `mave_nt` / `mave_nt_ref` (MAVE-HGVS strings), `mseq` / `mseq_no_adapt`. v3.0.0 additionally writes a JSON config file recording the full parameter set for reproducibility.

**automatic_naming = yes.** `oligo_name` is generated automatically; README example: `ENST00000357654.9.ENSG00000012048.23_chr17:43104102_1del_rc` (transcript.gene_chr:pos_mutator_rc). Output *filenames* are also auto-generated from the design: "The file names report chromosome, coordinates, strand and the sgRNA IDs associated with the targeton" (paper §2.5.3), e.g. `chr17_43115634_43115878_minus_sgRNA_ex2` (Supp. Table 3). Wiki *Advanced usage*: "When the transcript and gene ID cannot be determined the oligonucleotide name contains a placeholder string: `NO_TRANSCRIPT`."

**design_visualization = no.** VaLiAnT emits only CSV/VCF/JSON. Dependencies are `charset-normalizer, click, pydantic, pysam` (pyproject.toml) — no plotting library; no plotting/rendering module in `src/valiant/`. The design figure in the paper was made in third-party software: Fig. 3 legend — "sequence information modified from **Geneious Prime (version 2019.04) visualization**". The Discussion places visualization downstream and outside the tool: "VaLiAnT may be combined with downstream analysis software tools ... and to produce annotated visualizations of variant effect."

### BLOCK B — assay coverage

**assay_dms = yes.** This is the tool's entire purpose. Title: "an oligonucleotide library design and annotation tool for saturation genome editing and other deep mutational scanning experiments." Both SGE (endogenous locus, HDR templates) and cDNA DMS (expression-cassette) modes are demonstrated on BRCA1.

**assay_mpra = partial.** No MPRA-specific support and no MPRA claim; MPRA appears in the paper only as background ("MAVEs can assess coding and non-coding loci through approaches such as deep mutational scanning (DMS) and massively parallel reporter assays (MPRAs), respectively"). However, non-coding targetons *are* supported: wiki *Advanced usage* — "UTR regions and non-coding sequence are treated as if intronic sequence. Targetons that do not overlap any GTF/GFF2 feature won't be associated with a gene or transcript and their target regions will be considered intronic for the purposes of mutation generation"; and "If CDS-specific mutations are not required, then a GTF/GFF2 file need not be supplied" (wiki Input-files). So `snv`/`1del`/`2del*` saturation of a regulatory element is achievable and would produce usable MPRA-style oligos — but there are no barcodes, no reporter-cassette assembly, no motif insertion/scanning, no promoter/enhancer element grammar, and no MPRA example in the docs.

**assay_insilico = no.** Nothing in the paper, README, wiki, or source relates to designing libraries to probe sequence-to-function models. The one ML mention in the Discussion is about DMS data later *training* models ("as predicative models and machine learning increase in utility ... training datasets rich in biological context will be important"), not about in-silico perturbation experiments. No model interface, no scoring, no covariate export for surrogate modelling.

### BLOCK C — genomics integration

**genome_coordinates = yes.** `valiant sge` is coordinate-driven. "Targeton ranges are defined by 1-based indexing of genomic coordinates using 'ref_start' and 'ref_end' fields, with reference chromosome and strand also defined, as 'ref_chr' and 'ref_strand'" (paper §2.1); sequences retrieved from a local FASTA + FAI via pysam. Species and assembly are required CLI arguments; "The tool is species and genome-build agnostic" (§2.3). Example: BRCA1 exon 2, chr17:43115634–43115878, minus strand, GRCh38.

**transcript_models = yes.** "To collect gene and transcript information and to apply CDS-specific mutator functions, appropriate transcript annotation must be provided via a GTF/GFF2 file; only CDS, UTR and stop features are taken into consideration. **One transcript per gene is allowed** to avoid ambiguities in matching target regions." (paper §2.4.2). Wiki Input-files: "GTF/GFF files should not contain multiple transcripts... CDS and UTR features are required. The last CDS feature gets extended by one codon to include STOP in computation." Constraint worth quoting for the referee response: it is one transcript per run, not a full transcript-model engine.

**exon_intron_split_codons = yes.** Explicitly handled, with a published formula (paper eqns 1–2) and Supplementary Fig. S2: "Retrieval of additional positions from the reference sequence is necessary to obtain the context of partial codons... **The extra bases required at either end of the target may come from the same or an adjacent exon**." And: "SNVs can be introduced into partial, liminal codons in addition to complete codons; the annotation of the mutations in partial codons is informed by the exonic bases adjacent to the target" while "CDS-specific mutator functions that result in codon replacement or deletion are only applied to complete codons within CDS targets (partial, liminal codons are ignored)" (§2.5.1). Documented restriction: "No discrete target region (r1, r2, r3) within a targeton can span both coding and non-coding sequences; although the complete targeton can be divided into coding and non-coding regions" (§2.4.2).

**hgvs_input = no.** Every input is a table or a standard genomics file: targeton TSV, reference FASTA + FAI, GTF/GFF2, PAM-protection VCF, custom-variant VCF manifest CSV, codon-frequency CSV, JSON config, cDNA multi-FASTA + annotation TSV. Variants enter as VCF rows (POS/REF/ALT), not HGVS. HGVS appears only on the **output** side (`mave_nt`, `mave_nt_ref`). Checked README input table, wiki Input-files.md and cDNA-DMS-file-formats.md.

**vcf_vep_output = yes (SGE mode).** "`_.vcf` — Variant Call Format of generated and custom variants — Purpose: Bioinformatic analyses" (Supp. Table 2); README lists "Variant calls (VCF, SGE only)". Every generated variant, not just the custom ones, is emitted. Caveat: **cDNA mode emits no VCF** (relative coordinates), and the metadata CSV deliberately breaks VCF conventions — "the output metadata format does not follow VCF convention in reporting positions, reference and alternative sequences to favour streamlining of downstream processing" (§2.5.4). VEP itself is never named; the claim is VCF-compatibility, not a VEP integration.

**consequence_annotation = yes.** Codon-level consequence is computed and stored: metadata fields `ref_aa`, `alt_aa`, `mut_type`. "SNVs carry extra information when introduced into CDS targets ... namely whether the SNV results in a **synonymous, missense or nonsense** change" (§2.5.1); the classification is strand-aware ("For the same codon on different strands, the nucleotide changes are the same but the amino acid changes, and therefore their classification ... differ"), and is precomputed as a codon-indexed table at start-up. In-frame deletions are labelled via the `mutator` field (`inframe`). Limits: no splice-site consequence prediction, no frameshift/NMD calling, no VEP invocation — mutation-type annotation is codon-arithmetic only.

### BLOCK D — physical construction

**primer_design = no.** VaLiAnT appends *user-supplied* constant sequences (`--adaptor-5`, `--adaptor-3`: "DNA sequence to be added at the 5'/3' end of the oligonucleotide") and leaves flanking constant regions for user-chosen primers; it never designs a primer. In the prime-editing example the P5/P7, Golden Gate overhangs, sgRNA spacer, scaffold and PBS are all written by the user into the adaptor strings. cDNA example: "The 21-base flanking regions at targeton boundaries remain unmutated to facilitate primer binding" — the user picks them. No Tm/specificity calculation; no primer3 or equivalent in the dependency list (`click, pydantic, pysam, charset-normalizer`).

**codon_optimization = partial.** VaLiAnT is codon-usage-aware but does not optimise a sequence. It ships a human codon-frequency table (`src/valiant/data/default_codon_table.csv`, columns: triplet, amino acid, frequency, rank "RANKT/RANK2/.../RANKU") and accepts `--codon-table` for other species. It uses ranks to pick codons: `ala`/`stop` use "the most frequent triplet code" for Ala/stop; `aa` "exchanges each wild-type codon for the most frequent triplet code of each other amino acid"; `snvre` swaps a missense codon for "the next most frequent triplet code for the same missense change ... this allows for insights into the effect of codon sequence on missense changes." That is per-codon frequency-ranked substitution — not whole-CDS codon optimisation, harmonisation, CAI maximisation, or avoidance of structure/repeats.

**synthesis_constraints = partial.** Only length is enforced. `--max-length` (default 300 bp): "Oligonucleotides exceeding a given length (`max-length` option) will not be included in the unique oligonucleotide files and their metadata will be stored in separate files marked as 'excluded'." `--min-length` was added in v3.0.0. Discussion frames this as a synthesis-platform constraint: "At the time of writing, chemical synthesis at the scale of ~300 bp is possible on a large scale ... We have provided an optional filter to remove sequences above a user-defined length, which is configurable to accommodate for future increases in synthesis capability." No GC-content, homopolymer, repeat, secondary-structure, or restriction-site checking. The `_unique.csv` output is described as ready for "Synthesis submission" (Supp. Table 2).

### BLOCK E — engineering

**interface = CLI only.** "The tool is implemented as a standalone executable Python package exposing a command line interface" (§2.3). `pyproject.toml` declares one console script: `valiant = valiant.__main__:main`. Two subcommands (`sge`, `cdna`); `valiant -c config.json` for a JSON config. Distributed as source + Dockerfile + prebuilt image `quay.io/wtsicgp/valiant`. **No documented Python API** (the package is importable but no API reference exists in the README or any of the 11 wiki pages), **no web service, no GUI**. **Not on PyPI** — `pypi.org/project/valiant` is an unrelated dependency-audit package.

**license = AGPL-3.0-or-later.** Paper: "VaLiAnT is licensed under AGPLv3." GitHub API license field: "GNU Affero General Public License v3.0"; README: "GNU Affero General Public License v3.0 or later". Note for the referee response: AGPL is a copyleft license, unlike a permissive MIT/BSD.

**maintained = last commit 2024-04-22.** GitHub API `pushed_at` = 2024-04-22; latest commits on `develop` are "Merge tag '4.0.0' into develop" (2024-04-22, Luca Barbon). Tags: 1.0.0, 2.0.0, 3.0.0, 3.0.1, **4.0.0** (no GitHub "Releases" objects — tags only). Docker: `quay.io/wtsicgp/valiant:4.0.0` and `:latest` both pushed 2024-04-22 (older tags 3.0.1 = 2023-07-19, 3.0.0 = 2022-10-12, 2.0.0 = 2021-07-12, 1.0.0 = 2020-12-21). Repo not archived; 2 open issues; 6 stars. So: **~2.3 years since the last commit as of Aug 2026**, but fully installable and runnable today.

---

## 4. VaLiAnT's own documented example libraries

All under https://github.com/cancerit/VaLiAnT/tree/develop/examples with `run_*.sh` and `check_*.sh` scripts and expected outputs (`output_exp/`). Reference FASTA must be unpacked first (`cd sge && ./unpack_reference.sh`).

1. **BRCA1 nucleotide SGE library (`brca1_nuc`)** — exons 2–5 of BRCA1 (transcript NM_007294.4 / ENST00000357654.9, GRCh38), 4 targetons of 245–251 bp with 25 bp (or 20/41 bp) intronic extensions. Mutators: r1 `2del1, snv, 1del`; r2 `snvre, inframe, ala, stop, 1del`; r3 `2del0, snv, 1del`. Custom variants from ClinVar (2020-11-07) and gnomAD v3.0 merged in. PAM/protospacer protection edits at one sgRNA per targeton. P5/P7 adaptors, `--max-length 300`, `--revcomp-minus-strand`. Exon 2 alone: 1000 total sequences, 583 unique. Scripts: `run_brca1_nuc.sh`, `check_brca1_nuc.sh`.
2. **BRCA1 peptide/amino-acid SGE library (`brca1_pep`)** — same targetons, mutators `aa` + `inframe` only. Exon 2: 340 sequences; exons 2–5 total 340/500/580/920 (Supp. Table 3). Scripts: `run_brca1_pep.sh`, `check_brca1_pep.sh`.
3. **BRCA1 saturation prime editing library** — two pegRNAs over BRCA1 exon 2; variant RTT regions of 47 bp (pegRNA_a) and 49 bp (pegRNA_b) covering the 54 bp CDS + 4 bp flanking intron. 454 and 483 unique sequences. Adaptor-5 = P5 + Golden-Gate Type IIS end + sgRNA spacer + pegRNA scaffold; adaptor-3 = PBS + Type IIS end + P7. `--max-length 250`. Docs: wiki `Saturation-prime-editing-example`; scripts `run_prime_a.sh`, `run_prime_b.sh`.
4. **BRCA1 cDNA DMS library** — the 5592 bp BRCA1 CDS (protein NP_009225.1, 22 exons) inserted in silico into pCW57.1 (Addgene #41393) downstream of a TRE promoter (SalI/NheI digest). 40 targetons of 132–237 bp, r2 ranges 81–195 bp, tiled with overlap to cover all CDS. Mutators per targeton: `snv, 1del, snvre, ala, stop, inframe, aa`. 858–2092 unique oligos per targeton; **62,754 variants** in the concatenated library. Docs: wiki `cDNA-example`; scripts `examples/cdna/run.sh`, `check.sh`.
5. **NXRL example (repo only, not in the paper)** — targeton `chr17_31226400_31226678_plus_sgRNA_1146241047`; "This example applies background variants" (v4.0.0 feature: `--bg`, `--bg-mask`, `--force-bg-ns`, `--force-bg-indels`). Scripts: `examples/sge/nxrl/run.sh`, `check.sh`; includes a `*_valiant_config.json`.
6. **Docker usage example** — wiki `Docker-usage-example`, mounting inputs into `quay.io/wtsicgp/valiant`.

---

## 5. Notable capabilities not covered by the row list

- **PAM / protospacer protection edits** — user-defined synonymous or non-coding SNVs, keyed to an sgRNA ID via a VCF INFO tag, applied to *every* oligo of a tagged targeton so the HDR template is refractory to Cas9 re-cutting. Multiple edits per sgRNA and multiple sgRNAs per targeton supported. Supplementary Table S1 gives the sgRNA/edit selection heuristics (used by the authors, not automated by the tool).
- **Saturation prime editing (pegRNA RTT) design** — treats the RTT as the targeton, with `--revcomp-minus-strand` for strand-specific RTTs and adaptors for spacer/scaffold/PBS.
- **Background variants (v4.0.0)** — apply a haplotype/background VCF beneath the mutagenesis, with a BED mask to exclude regions and explicit guards `--force-bg-ns` / `--force-bg-indels` against non-synonymous or frame-shifting backgrounds; recorded per oligo in `background_variants` / `background_seq`.
- **Ingesting real variant catalogues** — ClinVar/gnomAD VCFs turned into oligos alongside the systematic scans, with per-file aliases and ID tags preserved; also usable standalone as "conversion of VCF annotation files to oligonucleotide sequences" (Discussion).
- **MAVE-HGVS output strings** (`mave_nt`, `mave_nt_ref`) for MaveDB-style interoperability.
- **Strand handling** — `--revcomp-minus-strand` produces reverse-complemented oligos for minus-strand targets / ssODN-orientation-specific SGE.
- **Deduplication and synthesis-ready export** — `_unique.csv` (unique sequence + name) for "Synthesis submission", plus `_meta_excluded.csv` for over-length sequences.
- **Reference-retrieval QC file** (`ref_sequences.csv`) and `--sequences-only` mode to check retrieved sequences before generating anything.
- **JSON run configuration** (`-c config.json`, and a config JSON written per execution) for reproducibility.
- **Species/assembly agnostic** with pluggable codon-frequency tables.
- **Targeton tiling recipe** for exons longer than one oligo (wiki *Advanced usage*), including in-frame overlapping tiles — but this is a documented manual pattern, not an automated tiling function.

## 6. Stated limitations

From the paper and wiki, verbatim or near-verbatim:
- One transcript per execution: "One transcript per gene is allowed to avoid ambiguities in matching target regions"; wiki: "GTF/GFF files should not contain multiple transcripts."
- "No discrete target region (r1, r2, r3) within a targeton can span both coding and non-coding sequences."
- CDS-specific mutators are rejected on non-CDS targetons: "CRITICAL:root:Invalid mutator 'snvre' for targeton!" (wiki *Advanced usage*).
- Partial/liminal codons are ignored by codon-replacement mutators (only SNVs reach them).
- cDNA mode: "All mutator functions except 'custom vcf' are currently supported" (paper §2.2) and no VCF output.
- Only "simple" custom variants: "Simple variants are currently supported, including substitutions, insertions, deletions and indels."
- Oligos above `--max-length` (default 300 bp) are dropped from the synthesis file; the authors tie this to synthesis-platform limits.
- Adaptors are appended to *all* sequences in a run, so per-targeton adaptor differences require separate runs: "targetons need to be processed individually as this function appends to all sequences in the library."
- Extension is expected to come from the community rather than from a plugin system: "We envision that, as opensource software, VaLiAnT will be further modified by the community. For example, by the addition of new mutator actions" — i.e. new mutation types require editing the source.
- Upstream design choices (targeton ranges, sgRNAs, protection edits) are the user's: "There is also scope for VaLiAnT to be combined with other software, such as upstream heuristic functions to select appropriate input information."

## 7. Availability today (checked Aug 2026)

- Repo **alive and public**, not archived: https://github.com/cancerit/VaLiAnT (default branch `develop`, 6 stars, 2 open issues).
- Last commit **2024-04-22**; latest tag **4.0.0** (same date). No GitHub Release objects, tags only.
- Docker image **quay.io/wtsicgp/valiant:4.0.0** and `:latest`, both pushed **2024-04-22** — pull-and-run today.
- Install from source: `python3.11 -m venv .env && source .env/bin/activate && pip install .` (requires Python ≥3.11; Linux/macOS native, Windows via WSL/Docker/Singularity). Deps: `charset-normalizer==3.3.2, click==8.1.7, pydantic==2.7.0, pysam==0.22.0` (all pinned).
- **Not on PyPI**, no conda package found, no web service.
- Documentation (11-page GitHub wiki) and runnable examples with expected outputs are present and current.
- Verdict: **installable and runnable today, but dormant** — no commits in ~2.3 years. Pinned deps mean it should still build; the Python ≥3.11 floor is fine on current interpreters.

## 8. Confidence notes

- Lowest confidence: `library_as_object` (partial) — this is a judgement call about whether a declarative TSV of targetons counts as a first-class library object. The factual anchors (declarative spec in one file; per-targeton outputs requiring concatenation) are solid; the label is arguable.
- `assay_mpra` (partial) — VaLiAnT never claims MPRA and has no MPRA-specific features, but non-coding saturation genuinely works. A referee could reasonably argue either "no" (no MPRA support) or "partial" (non-coding scans are possible). "Partial" is the conservative, defensible choice.
- `lazy_evaluation` (no) — asserted from the absence of any streaming/on-demand interface in a batch CLI plus the internal table-building modules; I did not run the tool or read the execution code line by line.
- `vcf_vep_output` (yes) — VEP is never named by the authors; the claim rests on a standard VCF of all generated variants being emitted in SGE mode. cDNA mode has no VCF.
- `codon_optimization` / `synthesis_constraints` (partial) — both are deliberately weak "partials"; the underlying facts (codon-frequency-ranked substitution; max/min length filter only) are certain.
- `dag_chaining`, `combinatorial_motif_place`, `barcode_generation`, `design_visualization`, `assay_insilico`, `primer_design`, `hgvs_input` — all "no" after checking the paper end to end, the full README (feature list, mutator table, CLI option table, metadata schema), all 11 wiki pages, the dependency list, and the complete `src/valiant/` and `src/valiant/mutators/` module listings. These are genuine absences, not failures to find.
- Everything about v4.0.0 (background variants, MAVE-HGVS, `--include-no-op-oligo`, JSON config) comes from the repo/CHANGELOG only — it postdates the 2022 paper. Any comparison table should be explicit that it describes VaLiAnT 4.0.0, not the published version.
