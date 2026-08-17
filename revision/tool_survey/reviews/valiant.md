# Adversarial review — VaLiAnT extraction

Reviewer pass, Aug 2026. Goal: falsify the extraction. Sources re-checked independently
(not merely re-read from the extractor's memo):

- `papers/Barbon2022_VaLiAnT_all.pdf` — full text pages 0–10 extracted with PyMuPDF (main text + Supp. Tables S1–S3, Supp. Figs S1–S2).
- `https://raw.githubusercontent.com/cancerit/VaLiAnT/develop/README.md` (875 lines, read in full).
- `CHANGELOG.md`, `pyproject.toml`, `examples/README.md`, `examples/sge/brca1/run_brca1_nuc.sh`,
  `examples/sge/brca1/parameter_input_files/brca1_nuc_targeton_input.txt`,
  `examples/sge/nxrl/run.sh` and its `*_valiant_config.json`.
- Source read directly, not just listed: `src/valiant/mutator_type.py`, `src/valiant/int_pattern_builder.py`,
  `src/valiant/meta_table.py` (head), `src/valiant/sge_proc.py` (full, 568 lines).
- Full recursive git tree of `develop` (`gh api .../git/trees/develop?recursive=1`).
- GitHub API via `gh`: repo metadata, commits on `develop`, tags, releases, all issues (open + closed).
- Wiki raw markdown: `Home`, `Input-files`, `Output-files`, `Advanced-usage`, `cDNA-example`,
  `cDNA-DMS-file-formats`, `Saturation-prime-editing-example`, `Docker-usage-example`.
- `quay.io` tag API, `api.anaconda.org/package/bioconda/valiant`, `pypi.org/pypi/valiant/json`.

**Headline:** the extraction is unusually good. Every "no" I tried to break held up against the
current v4.0.0 source, not just the 2022 paper. There is **no capability that is falsely denied**.
The problems are (a) two ambiguous yes/no cells whose column semantics are undefined
(`interface`, `maintained`), (b) one factually wrong sentence inside an otherwise-correct "no"
(`combinatorial_motif_place`), and (c) several concrete factual errors in
`documented_examples` that a VaLiAnT author would spot immediately.

---

## Independently re-verified facts (all confirmed)

| Claim | Status |
|---|---|
| Last commit `develop` 2024-04-22, "Merge tag '4.0.0' into develop", Luca Barbon | confirmed |
| `pushed_at` 2024-04-22; archived=false; default branch `develop`; 6 stars; 2 open issues | confirmed (forks=3, not stated) |
| Tags 1.0.0 / 2.0.0 / 3.0.0 / 3.0.1 / 4.0.0; **zero** GitHub Release objects | confirmed (`gh api .../releases` → length 0) |
| quay.io tags `latest` and `4.0.0` both pushed 2024-04-22 | confirmed |
| Not on PyPI (`pypi.org/project/valiant` = Duncan Dickinson's dependency-audit tool v0.2.3) | confirmed |
| Not on bioconda (`{"error": "\"valiant\" could not be found"}`) | confirmed |
| Deps exactly `charset-normalizer==3.3.2, click==8.1.7, pydantic==2.7.0, pysam==0.22.0`; `requires-python >=3.11`; one console script `valiant = valiant.__main__:main` | confirmed |
| AGPL-3.0 (`gh api` → `"license": "AGPL-3.0"`; README LICENSE block; paper Availability) | confirmed |
| 32-column metadata schema, exactly as listed | confirmed against `META_CSV_FIELDS` in `src/valiant/meta_table.py` |
| Supp. Table 2: only `ref_sequences.csv` is execution-specific | confirmed |
| Supp. Table 3 counts (583/740/825/1185 nt; 340/500/580/918 aa unique) | confirmed |
| `src/valiant/mutators/` = codon.py, deletion.py, snv.py, snv_re.py, utils.py | confirmed |
| No plotting / barcode / primer / motif module anywhere in `src/valiant/` | confirmed against full tree |

---

## Capability-by-capability

### 1. `library_as_object = partial` — **supported**
Verified stronger than the extraction states, and this matters. A single `valiant sge`
invocation processes **many targetons across many contigs**: `run_sge` → `for contig,
targeton_configs in exp_config.targeton_configs.items(): proc_contig(...)` →
`for targeton_config in tcs: proc_targeton(...)`. The shipped `brca1_nuc_targeton_input.txt`
is one file with four targeton rows, run by one command. So the declarative-spec half of
"partial" is real, not generous. The other half is equally real: `MetaTable.to_csv(conn,
targeton_name)` writes one CSV per targeton, and the paper's own phrase for the cDNA design is
"the **final concatenated library**" (§3.2).

**Action:** add "one invocation handles multiple targetons across multiple contigs" to the
evidence. Without it, a referee can call the cell dismissive.

### 2. `dag_chaining = no` — **supported (code-verified)**
`src/valiant/mutator_type.py` is a closed `Enum` of exactly seven types
(`DEL, SNV, SNV_RE, IN_FRAME, ALA, STOP, AA`) with
`DEPENDENT_MUTATOR_TYPES: dict = {MutatorType.SNV_RE: {MutatorType.SNV}}` — literally the only
composition in the system, and it is hard-coded, exactly as the extraction says. No registry,
no entry points, no plugin hook in `pyproject.toml`. Pipeline order is fixed in `proc_contig`
(background → PPE → custom → per-targeton mutators → adaptors).

### 3. `lazy_evaluation = no` — **supported, and the caveat can now be dropped**
The extraction hedges: "inferred from documentation and module structure, not from running the
code." I read the code. `run_sge` opens an **in-memory SQLite** connection (`get_db_conn`,
`init_db`, `src/valiant/data/ddl.sql`), inserts exons/PAM/custom/background variants, calls
`targeton.process(...)` to materialise every mutation into the DB, then `MetaTable.to_csv`
drains it to file. `sequences_only` short-circuits with `sys.exit(0)` after writing the QC file.
There is no generator, no partial-materialisation option, no public API. The "no" is now
directly verified, not inferred — **strengthen the confidence note accordingly**.

### 4. `mixed_mutagenesis_one_pool = yes` — **supported**
Verified against the real input file (not just paper Table 1):
`"(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"`. Multiple mutator types per
region, different types per region, one pooled output. The caveat in the evidence (exhaustive
single-event scans only; no sampled/random mode; no pairwise mutator; `--include-no-op-oligo`
is one WT oligo, not replicates) is accurate and should be kept verbatim.

### 5. `combinatorial_motif_place = no` — **value correct, but the evidence contains a false sentence**
This is the one place a VaLiAnT author could land a hit. The evidence says:

> "The only insertable fixed sequences are `--adaptor-5` / `--adaptor-3`"

That is **not true**. README §*Custom variants*: "Only simple variants such as the following are
supported: substitutions, **insertions**, deletions, indels." A user can supply a VCF with an
insertion of any sequence at any position inside a targeton — including a TF-binding motif —
and VaLiAnT will build one oligo per such insertion, name it, annotate it, and emit it. That is
sequence insertion at user-chosen positions.

The value `no` still stands, because the *combinatorial* half is absent: VaLiAnT does not
enumerate motif × position, does not permute or orient elements, has no motif/element concept,
and requires the user to have already written every insertion out longhand in a VCF. But the
evidence must say that, not the false exclusivity claim.

**Suggested replacement wording:** "VaLiAnT has no motif or element concept and no
insertion-scan machinery. Fixed sequence can enter an oligo two ways: `--adaptor-5`/`--adaptor-3`,
appended once to every oligo in the run; or as an insertion in a user-supplied custom-variant
VCF, which VaLiAnT applies verbatim at the given position. Neither enumerates placements —
there is no motif × position combinatorics, no orientation or permutation control, and the
user must author each insertion by hand. Checked the full `src/valiant/` tree: no motif,
element, or insertion-scan module."

### 6. `barcode_generation = no` — **supported**
Nothing in the 32-field schema, the CLI tables, the config schema, all 8 wiki pages, or the
source tree. No GC/edit-distance/Tm machinery of any kind. Held under attack.

### 7. `per_sequence_provenance = yes` — **supported, and understated in detail**
All 32 fields confirmed in `META_CSV_FIELDS`. Two things the extraction missed that make the
"yes" stronger:
- The output VCF carries `INFO` tags that link **back** to the oligo: `SGE_OLIGO` (oligo_name),
  `SGE_SRC` (mutator), `SGE_VCF_ALIAS`, `SGE_VCF_VAR_ID`. Provenance is bidirectional.
- `vcf_var_in_const` records whether a custom variant landed in a non-mutable constant region.

### 8. `automatic_naming = yes` — **supported**
`get_sge_oligo_name` / `get_cdna_oligo_name` / `get_sge_oligo_no_op_name` in `meta_table.py`;
`NO_TRANSCRIPT` placeholder confirmed in code and wiki. Note cDNA output filenames are
`<seq_id>_<md5>_meta.csv` (hash-derived) rather than coordinate-derived — still automatic.

### 9. `design_visualization = no` — **supported**
No plotting dependency, no rendering module, and the paper's own Fig. 3 legend credits Geneious
Prime. The repo does ship static PNGs under `examples/images/wiki/` and a hand-built GenBank
model (`examples/cdna/construct_generation/BRCA1_cdna_pCW57.1_model.gb`) — documentation
assets, not tool output. Worth one clause in the evidence so a referee cannot say "there are
figures in the repo".

### 10. `assay_dms = yes` — **supported.** Trivially so.

### 11. `assay_mpra = partial` — **supported**
I tried to break this both ways and could not. Wiki *Advanced usage* genuinely permits non-coding
targetons ("Targetons that do not overlap any GTF/GFF2 feature ... will be considered intronic
for the purposes of mutation generation"), and `sge_proc.proc_contig` has an explicit
`else:` branch for `config.gff_fp` being unset. Equally, there is no barcode, no reporter
cassette, no MPRA example. "partial" is the only label that survives both attacks. Keep it.

### 12. `assay_insilico = no` — **supported.** No model, scorer, or covariate export.

### 13. `genome_coordinates = yes` — **supported.** Confirmed in the real targeton TSV and `fetch_sequence`/pysam FASTA path.

### 14. `transcript_models = yes` — **supported**
One nuance for the record: the paper says "One transcript per **gene** is allowed" and the README
says the GTF "should only contain features for one transcript per gene", but the wiki is
stricter — "One specific transcript must be supplied for an execution. *GTF/GFF files should not
contain multiple transcripts.*" The code confirms the wiki: `sge_proc.proc_contig` builds a
single `Transcript` per contig/strand and carries the comment
`# TODO (future multi-transcript): clear existing exons from the database`. The extraction's
framing ("one transcript per run") is the operative one and is correct; keep the wiki quote as
the anchor rather than the paper quote, which is looser.

### 15. `exon_intron_split_codons = yes` — **supported; one v4 improvement missed**
CHANGELOG 4.0.0 adds "Support for out-of-frame CDS targeton regions whose 5' and 3' extensions
span both adjacent and distal bases" — i.e. v4 handles harder split-codon geometry than the
published eqns (1)–(2) description. The extraction cites only the paper here. Add it; it makes
the "yes" unimpeachable.

### 16. `hgvs_input = no` — **supported**
`src/valiant/mave_hgvs.py` is write-only (`get_mave_nt`, consumed by `MetaTable`). Every loader
in `src/valiant/loaders/` is bed/csv/fasta/gtf/vcf/tsv. No HGVS parser exists.

### 17. `vcf_vep_output = yes` — **supported; the file description is stale**
Since v3.0.0 there are **two** VCFs per targeton, not one: `*_pam.vcf` (REF includes PAM
protection edits) and `*_ref.vcf` (REF excludes them), each carrying the `SGE_*` INFO tags. The
extraction's `_.vcf` wording comes from Supp. Table 2 (paper era). Update it — the current
behaviour is richer than the cell describes, and this is exactly the "describes the 2022 paper,
not v4.0.0" trap the extraction itself warns about. The "VEP is never named by the authors"
caveat is correct and should stay.

### 18. `consequence_annotation = yes` — **supported**
Plus two v4 facts worth adding: CHANGELOG 4.0.0 "report the type of mutation for the `aa`,
`ala`, and `stop` mutators (`mut_type` field)" (so consequence is no longer SNV-only), and
"in cDNA mode, mutations preserving stop codons are now annotated as nonsense". Also
`pam_mut_annot` classifies each PAM edit as `syn|mis|non|ncd`, and background variants are
classified synonymous / non-synonymous / frame-shifting with a hard error unless forced
(`is_variant_nonsynonymous`, `is_variant_frame_shifting` in `sge_proc.py`). The stated limits
(no splice prediction, no NMD, no VEP) remain correct.

### 19. `primer_design = no` — **supported.** No Tm/specificity/dimer code; no primer3-class dependency.

### 20. `codon_optimization = partial` — **supported**
Rank-driven codon choice only (`codon_table.py`, `codon_table_builder.py`, ranks
`RANKT/RANK2/.../RANKU`). Note the README quantifies `aa`: "19 mutated sequences for each codon
mapping to an amino acid ... and 20 for each stop codon" — worth citing. "partial" is generous
toward VaLiAnT, which is the safe direction.

### 21. `synthesis_constraints = partial` — **supported**
`min-length`/`max-length` only, confirmed in the CLI table and the config schema
(`minOligoLength`, `maxOligoLength`). One extra input-level check exists that is not a synthesis
constraint but is adjacent: "Ambiguous nucleotides are not allowed in the reference sequence.
Soft-masking is ignored." Mention it or a referee will. Still no GC / homopolymer / repeat /
secondary-structure / restriction-site checking anywhere in the tree.

### 22. `interface = yes` — **UNSUPPORTED AS ENCODED (column semantics undefined)**
The evidence text is entirely accurate — "CLI ONLY (plus Docker) ... NO documented Python API,
NO web service, NO GUI" — but the value `yes` is not derivable from it without knowing what the
column asks. If the column means "has a usable interface", `yes`. If it means "has a
programmatic / library API" — which is the natural reading in a table whose point is that
PoolParty is a Python library — the honest answer for VaLiAnT is **no**. Publishing `yes` next
to evidence that says "CLI only, no Python API" is an internal contradiction a referee will
quote back.

**Action:** define the column, and make the cell carry the qualifier, e.g.
`CLI only (Docker/Singularity images; no Python API, no web service, no GUI)`. Do not ship a
bare `yes`. Confirmed independently: `pyproject.toml` `[project.scripts]` has exactly one entry;
no `[project.entry-points]`; no `docs/` directory in the tree; no API page among the wiki pages.

### 23. `license = yes` — **supported, with a presentation warning**
AGPL-3.0-or-later confirmed three ways. `yes` is fine if the column means "open-source
licensed". But AGPL copyleft is materially unlike the MIT/BSD licences of several other tools
in this survey (seqpro, tangermeme are MIT). If the table has a single `license` column, put the
SPDX identifier in the cell, not `yes` — otherwise the table erases a real difference and looks
careless.

### 24. `maintained = yes` — **OVERSTATED (but in VaLiAnT's favour, so low attack risk)**
Verified: zero commits since 2024-04-22, i.e. ~28 months dormant at the Aug 2026 assessment
date. Zero GitHub Release objects. Docker images last pushed 2024-04-22. And a detail the
extraction did not surface: the **two open issues are both automated Snyk dependency-bot PRs**
("[Snyk] Fix for 5 vulnerabilities", #10 and #11) — left unattended. There is no human issue
traffic at all; the newest closed issue is #9 "Release 3.0.1" (2023). Every closed issue is
authored by the maintainers themselves.

`yes` therefore claims active maintenance the evidence does not show. Suggested `partial`
("dormant since 2024-04-22; not archived; installable and containerised"). Note this error
flatters VaLiAnT, so no referee will attack it — but the survey should state its criterion
(e.g. "commit within 24 months") so the column is defensible across all eleven tools.

---

## Errors in `documented_examples` (fix these — they are checkable in 30 seconds)

1. **NXRL example — wrong flags.** The extraction says it uses "`--bg` with `--bg-mask` BED,
   `--force-bg-ns`, `--force-bg-indels`". The actual
   `chr17_31226400_31226678_plus_sgRNA_1146241047_valiant_config.json` has
   `"maskBackgroundFilePath"` **absent entirely** (no BED mask is used), and sets
   `"forceBackgroundNonSynonymous": true`, `"forceBackgroundFrameShifting": true`. Also, the
   example is invoked as `valiant -c <config>.json` — so it doubles as the repo's only worked
   demonstration of the JSON-config entry point, which is worth saying.
2. **brca1_nuc — `--max-length 300` is not passed.** `run_brca1_nuc.sh` passes only `--pam`,
   `--vcf`, `--revcomp-minus-strand`, `--adaptor-5`, `--adaptor-3`, `--gff`. 300 bp is simply
   the default. The paper says the option was "used"; the shipped script does not set it.
3. **brca1_nuc — exon 5 mutators differ.** The extraction gives one action vector for all four
   targetons. In the real input file, exons 2–4 use r1 = `(1del,2del1,snv)` but **exon 5 uses
   r1 = `(1del,2del0,snv)`** (its r1 is 20 bp, even). Also the r2 vector in the file is
   `(1del,snvre,inframe,stop,ala)` — no `2del*`, matching the paper, but the extraction's
   ordering should be quoted from the file, not paraphrased.
4. **cDNA example** also ships an `output_exp_v3/` directory of v3-era expected outputs
   alongside `output_exp/`, and a `construct_generation/` directory with the GenBank plasmid
   model and `cdna_targeton_plan.txt` — evidence the 40-targeton tiling plan was authored by
   hand, which reinforces "tiling is a manual recipe, not a function".
5. **Docker usage wiki page is stale** — it pins `VERSION=2.0.0`, four versions behind. Minor,
   but if the survey cites documentation quality, this is the counterexample.
6. The examples require an unpack step first (`sge/unpack_reference.sh`,
   `brca1_prime_editing/unpack_inputs.sh`); validation is by `md5sum` against `output_exp/`.
   The extraction says "expected outputs with check scripts", which is right, but the md5
   comparison is a stronger reproducibility claim worth making.

---

## Capabilities the extraction missed entirely

1. **The deletion mutator is fully parametric in v4.0.0.** README: "Non-overlapping stretches of
   nucleotides of a given length are deleted starting from a given offset. Format:
   `<SPAN>del[<OFFSET>]`", and CHANGELOG 4.0.0: "Mutator: deletion pattern span and offset can
   be set to **any valid value**". Confirmed in `int_pattern_builder.py` (`IntPatternBuilder(offset,
   span)`, span > 0 the only constraint) and `mutator_type.py`
   (`PARAMETRIC_MUTATOR_TYPES = {MutatorType.DEL}`). So `3del1`, `6del0`, `10del2` … all work —
   arbitrary tiled deletion scans, not just the `1del/2del0/2del1` triple the extraction implies
   throughout. This is the single most under-described capability in the file, and it is the
   kind of thing an author-referee would notice. It does not change any cell value, but it
   belongs in `additional_capabilities`.
2. **Two output VCFs with back-links** (`*_pam.vcf` / `*_ref.vcf`; `SGE_OLIGO`, `SGE_SRC`,
   `SGE_VCF_ALIAS`, `SGE_VCF_VAR_ID`) — see §17 and §7.
3. **Order-invariant deduplication.** "When multiple oligonucleotides have the same sequence, the
   first name in lexicographic order is chosen" (README), and CHANGELOG 2.0.0 "Make unique
   oligonucleotide name selection order-invariant" — a deliberate reproducibility guarantee.
4. **Automatic conflict detection between mutations and background variants.** README:
   "Mutations that overlap background variants are discarded; warnings are always raised to
   identify them, with their positions expressed in absolute genomic coordinates for custom
   variants, and targeton-relative coordinates for pattern variants."
5. **Background-variant scoping is computed, not user-specified.** README describes an automatic
   range expansion: background variants are applied "in the minimal range of positions that
   spans at least the entire CDS, further extended to the boundaries of any targeton overlapping
   the CDS, and finally to the start and end position of the first and last background variant
   intersecting the resulting range" — with an exon lift-over (`transcript.lift_exons(gpo)`,
   `GenomicPositionOffsets`). Non-trivial coordinate machinery worth one sentence.
6. **In-memory SQLite as the generation engine** (`db.py`, `data/ddl.sql`, `queries.py`,
   `sql_gen.py`). Relevant context for any scalability comparison, and it is what makes the
   `lazy_evaluation = no` finding verifiable.
7. **`vcf_var_in_const`** — flags custom variants that landed in a constant (non-mutable) region;
   confirms the paper's "custom variants are still incorporated into constant regions".
8. **Singularity is a first-class documented deployment path**, with its own README section and
   `singularity exec --cleanenv -B ...` invocation — not merely a Windows workaround as the
   availability text implies.
9. **cDNA mode multiplexes** — the multi-FASTA can hold several cDNAs (`seq_id`), each with
   several targeton rows, in one run.
10. **`aa` mutator is quantified**: 19 substitutions per amino-acid codon, 20 per stop codon
    (README) — useful for library-size comparisons against PoolParty.
11. **Repo has 3 forks** (not stated); no bioconda package (checked and confirmed absent, which
    the extraction asserted without showing the check — it is correct).

---

## Bottom line

No capability is falsely denied. Do not soften any of the seven "no" cells. Required fixes,
in priority order:

1. Rewrite the `combinatorial_motif_place` evidence — the "only insertable sequences are the
   adaptors" sentence is false and is the most attackable line in the file.
2. Define and qualify the `interface` cell; a bare `yes` beside "CLI only, no Python API" is
   self-contradicting.
3. Fix the three `documented_examples` errors (NXRL flags, `--max-length`, exon 5 action vector).
4. Reconsider `maintained = yes` → `partial`, and state the survey's criterion.
5. Refresh `vcf_vep_output` and `exon_intron_split_codons` from v4.0.0 rather than the paper.
6. Add the parametric-deletion generality to `additional_capabilities`.
7. Drop the "not verified in code" hedge on `lazy_evaluation` — it is now verified.
