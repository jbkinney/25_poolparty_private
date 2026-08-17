# Repair log — `final/valiant.md`

**Date:** 2026-08-14
**Record repaired:** `revision/tool_survey/final/valiant.md` (edited in place)
**Audits processed:** `citation_audit/valiant.md` (34 findings) + `factcheck/valiant.md` (A: 11, B: 11, C: 1)

**Totals:** 42 findings applied · 3 rejected · 12 escalated.
The 42 applied findings became **37 distinct corrections** — 5 findings (fact-check A1, A3, A7, A9,
A11) restate citation-audit findings 5, 10, 20, 25 and 14 and were fixed once each.
7 of the 12 escalations are audit-section-B items grouped into a single question (E5).

## Verification method

Every finding was checked against primary source before any edit. No finding was applied on the
audit's assertion alone.

- **Repo:** `git clone --branch develop https://github.com/cancerit/VaLiAnT` into a scratch dir.
  Head = `8796cc112dafd4919fec59913f58cd2be87c45eb`, 2024-04-22 — the same commit both audits cite.
  Files read directly, counts recomputed with `git ls-tree`.
- **Wiki:** `git clone https://github.com/cancerit/VaLiAnT.wiki.git` — all 11 pages read.
- **Paper:** `papers/Barbon2022_VaLiAnT_all.pdf` extracted with PyMuPDF (11 pp), quotes matched by
  string search, section boundaries located by heading line numbers.
- **Live APIs:** `gh api` for repo metadata / releases / issues; `curl` for quay.io tags, PyPI JSON,
  anaconda.org.
- **Shipped examples:** row counts, adaptor-literal lengths, VCF `INFO` fields and expected-output
  files all recomputed from the cloned tree.
- No install, no execution of VaLiAnT itself (the record's own "no install, no execution" constraint
  is preserved — see escalation E3).

Nothing was written outside `revision/tool_survey/`.

---

## APPLIED (42 findings → 37 corrections)

### From the citation audit

| # | Finding | Correction made | How verified |
|---|---|---|---|
| C1 | `mutator_type.py` code block labelled "verbatim" but is edited | Label changed to "condensed — comments added, type annotations and the `ANNOTABLE_MUTATOR_TYPES` set omitted". Block content left as-is (substantively accurate). | Read `src/valiant/mutator_type.py` in full: source has `# Parametric span` etc. on separate lines, carries `set[MutatorType]` annotations, qualified members, and a fourth set the block omits. |
| C2 | Wiki has 11 content pages, not 8 | "(8 content pages)" → "(11 content pages)"; the four "all 8 wiki pages" / "8-page wiki" claims → 11; sources table row now lists all 11 page names with an honest retrieval note. | Cloned the wiki: exactly 11 `.md` pages. **I then read the three previously unconsulted pages** (`Configuration`, `Saturation-prime-editing`, `Software-and-dependencies`) and grepped them for barcode / MPRA / plot / visuali / HGVS / primer3 / Tm / "python api" — **zero hits**, so every absence claim the record rests on those pages survives at 11. |
| C3 | Record encodes 8 `no` values, not 7 | "every one of the seven 'no' cells" → "eight"; §10 "the seven 'no' cells (…)" → "seven of the eight 'no' cells (…; the eighth, `lazy_evaluation`, is covered in the next bullet)". | Counted the 24 capability entries in §4: `no` = `dag_chaining`, `lazy_evaluation`, `combinatorial_motif_place`, `barcode_generation`, `design_visualization`, `assay_insilico`, `hgvs_input`, `primer_design` = 8. |
| C4 | `mutators/` has 6 files, not 5 | Added `__init__.py` to the list. | `ls src/valiant/mutators/`. |
| C5 / A1 | Pipeline order reversed: custom variants run **before** pattern mutators | §2 and the `dag_chaining` entry both re-ordered to "background → PAM → custom VCF → mutators → adaptors → length filter → dedup"; attribution widened from `sge_proc.py` alone to `sge_proc.py` / `targeton.py` / `meta_table.py`. | `targeton.py`, `Targeton.process()`: `self._process_custom_variants(...)` then `self._process_pattern_variants(...)`. Adaptor concat (`get_full_oligo`), length filter (`info.eval_in_range`) and dedup (`unique_oligos`) all live in `meta_table.py`. |
| C6 | `MetaTable.to_csv` does not write "one CSV per targeton" | → "writes **per-targeton files, not one library-wide file**". | `meta_table.py:340–355` names `_unique.csv`, `_meta.csv`, `_meta_excluded.csv`, `_ref.vcf`, `_pam.vcf`; `:666–686` unlinks the empty ones. The record's very next sentence already listed three CSVs — an internal contradiction now removed. |
| C7 | "no generator, no streaming" is false | Both clauses deleted; replaced by a parenthetical that the only row-at-a-time work is `MetaTable.to_csv` draining the finished DB with `fetchone`. | `loaders/csv.py:53` `yield r` (real generator); `meta_table.py:488–489` `it = cur.execute(...)` / `while r := it.fetchone()`. The `lazy_evaluation = no` basis (full materialisation into SQLite first, no partial-materialisation option, no API) is untouched. |
| C8 | "Full library counts are reported per run (paper Table 2, Supp. Table 3)" is unsupported | Sentence deleted. | Paper Table 2 / Supp. Table 3 are static example summaries. In v4, `common_cli.finalise()` calls only `log_excluded_oligo_counts`; grep of every `logging.info/warning` in `src/valiant/` found no total-library count anywhere. |
| C9 | "ten mutator/region combinations" wrong | → "twelve generated mutator/region combinations (3 in r1, 6 in r2, 3 in r3)". | Counted non-zero generated cells in paper Table 2: r1 = 2del1/1del/snv; r2 = 1del/inframe/snv/snvre/ala/stop; r3 = 2del0/1del/snv. 3+6+3 = 12. |
| C10 / A3 | `--include-no-op-oligo` is not a general WT switch | Sentence replaced: still "a **single** sequence, not replicates", but now states it fires only when a background or PAM edit applies and carries that edited reference rather than true WT, with the CLI help quoted. | `meta_table.py:422–425`: `include_no_op_oligo = self.opt.include_no_op_oligo and (len(self.ppe_mut_types) > 0 or len(self.bg_variants) > 0)`; the written sequence is `self.alt_seq.s`. `sge_cli.py` help string quoted verbatim. |
| C11 | 93 files → 87 | Corrected. | `git ls-tree -r --name-only HEAD src/valiant \| wc -l` = 87. (461 recursive entries incl. dirs is correct and left alone.) |
| C12 | "no … uniqueness machinery" contradicted by `_unique.csv` dedup | → "no GC, edit-distance or Tm machinery … and no sequence-uniqueness *design* step — the only uniqueness handling is post-hoc deduplication of identical oligos into `_unique.csv`". | `grep -rniE "melting\|\btm\b\|hamming\|levenshtein\|edit_distance\|barcode\|gc_content\|gc_frac" src/valiant/` → **0 hits**, so the GC/Tm/edit-distance half stands. `meta_table.py` `unique_oligos` defaultdict + lexicographic name selection is the dedup the record itself describes in §6.7. |
| C13 | Not every oligo gets a `_meta.csv` row | → "in `_meta.csv`, or in `_meta_excluded.csv` if the length filter rejected it". | Code branches `fh = meta_fh` / `fh = meta_excl_fh` on `info.eval_in_range`. Shipped exon-2 `_meta.csv` = 1000 rows, `_meta_excluded.csv` = 1 row, and that oligo name is absent from `_meta.csv` (grep = 0). |
| C14 / A11 | `SGE_VCF_ALIAS` / `SGE_VCF_VAR_ID` are not on every VCF record | Corrected in all three places (`per_sequence_provenance`, `vcf_vep_output`, §6.6): `SGE_OLIGO` + `SGE_SRC` on every record, alias/ID for custom variants only. README quote extended to its full, accurate form. | `vcf_writer.py:150–160`: `if vcf_alias and vcf_var_id:` / `if vcf_var_id:`. README §*Variant file* marks both "(Optional) … only for `custom` variants". Committed exon-2 `_pam.vcf`: 651 records with `SGE_SRC;SGE_OLIGO` only, 349 with alias/ID. |
| C15 | "All generated variants are emitted" false under the length filter | → "**All in-range** generated variants … (oligos rejected by the length filter appear in neither VCF)". | VCF writes sit inside the `info.eval_in_range(...)` branch (`meta_table.py:585–625`). The exon-2 excluded oligo is absent from `_pam.vcf` (grep = 0). |
| C16 | `loaders/` has 16 files, and not all are format loaders | "(17 files)" → "(16 files)"; "Every loader …" → "Every **format** loader … (the rest of the directory is config, error and utility modules)". | `git ls-tree` = 16; directory contains `cdna_experiment_config.py`, `mutator_config.py`, `errors.py`, `utils.py`, `experiment.py`, etc. No HGVS parser found, so the underlying `hgvs_input = no` evidence is unaffected. |
| C17 | Fused `5'/3'` adaptor "quotation" | Split into the two real option descriptions. | README:223–224 are two separate rows: "…at the 5' end of the oligonucleotide." / "…at the 3' end…". |
| C18 | `--adaptor-5` literal is 121 nt, not 120 | Corrected in both places (`primer_design`, §5 item 3). | Measured the literal in `run_prime_a.sh` and `run_prime_b.sh` with a regex + `len()`: 121 in both. `--adaptor-3` = 38, as the record says. |
| C19 | "21-base flanking regions…" presented as a quotation but is a paraphrase | De-quoted; replaced with the fact plus the wiki's actual quoted fragment "allows for primer binding to amplify and clone specific targetons". | `wiki/cDNA-example.md:13` — the page says `targeton_start` is 21 bp upstream of `r2_start` etc.; the record's sentence appears nowhere. |
| C20 / A7 | SNVRE rule misstated and the "insights" quote mis-located | → "the **top-ranking** synonymous triplet of the SNV's amino-acid product — the second-ranked one if that triplet is already top-ranking (README §*SNVRE*) — which the paper says '…' (§2.1)". `paper §2.1` added to the *Source:* line. | README:383–395 states the top-ranked / second-ranked rule. The quoted fragment is at paper line 159, inside §2.1 (lines 119–321), not in any README section the paragraph cited. |
| C21 | cDNA targetons are 123–237 bp, not 132–237 | → "targetons 123–237 bp (the paper says 132–237; the shortest shipped row, 901–1023, is 123 bp)". | Recomputed inclusive lengths for all 40 rows of `examples/cdna/input/cdna_targeton.tsv`: sorted min = 123, second = 132, max = 237. The paper (§3.2) and cDNA wiki both say 132. |
| C23 | `VERSION=2.0.0` is 3 releases behind, not 4 | → "three tagged releases behind (`3.0.0`, `3.0.1`, `4.0.0`)". | `git tag -l` = 1.0.0, 2.0.0, 3.0.0, 3.0.1, 4.0.0. |
| C25 / A9 | "only SNVs reach" partial/liminal codons | → "SNVs, parametric deletions and custom variants can still reach them". | Paper §2.5.1 scopes the restriction to CDS-specific codon-replacement/deletion mutators. `mutators/deletion.py` `DeletionMutator.get_variants` has no codon gating; custom variants are applied across the whole targeton (paper §2.1, README §*Custom variants*). |
| C26 | "silently discarded (with warnings)" self-contradictory and contrary to the README | → "**discarded**; the README states that warnings are always raised to identify them". | README:120 verbatim: "warnings are **always** raised to identify them". |
| C27 | Commit message not verbatim (em dash inserted) | → subject in quotes, "with body 'Background variant support'". | `git log -1 --format='%B'` on the develop head. |
| C28 | "no human issue traffic at all" is false | Replaced in both places (§4 `maintained`, §9 verdict) with the accurate statement: external traffic is minimal but not zero — issue #8, 2023-06-27, `delagee`, `author_association: NONE`, answered and closed by a maintainer 2023-07-18, the day 3.0.1 shipped a MAVE-HGVS fix. | `gh api /repos/cancerit/VaLiAnT/issues?state=all` — #8 is not a PR and its author is external; its single comment is from `lbarbon` (COLLABORATOR) at 2023-07-18T14:36:58Z; CHANGELOG dates 3.0.1 to 2023-07-18 with a `mave_nt` fix. |
| C29 | `run_*.sh` / `check_*.sh` / `output_exp/` glob and naming wrong | Globs → `run*.sh` / `check*.sh`; "against `output_exp/`" → "against that directory"; "Inputs must be unpacked first" → "The shared SGE reference and the prime-editing inputs must be unpacked first". | `find examples -name "*.sh"`: cdna and nxrl use `run.sh`/`check.sh`. Expected dirs are `brca1_nuc_output_exp`, `brca1_pep_output_exp`, `output_a_exp`, `output_b_exp`, `cdna/output_exp`, `nxrl/output_exp`. Only two unpack scripts exist and neither is in `examples/cdna`. |
| C33 | "Only simple variants…" quote is a compression of prose + a bullet list | Fixed in both places: quotation now ends at the colon, the four types follow unquoted. | README:467–472 — introductory sentence then a four-item Markdown list. (Paper §2.5.2 uses different wording again.) |
| C34 | `SalI/NheI` does not match the cited sources | → "SalI/NheI digest — both the paper and the cDNA wiki page spell the second enzyme `NehI`". | Paper line 611 and `wiki/cDNA-example.md:1` both read "SalI and NehI". The record's silent correction is now documented rather than hidden. |

### From the fact-check, section A

| # | Finding | Correction made | How verified |
|---|---|---|---|
| A2 | "multi-base and combinatorial variants can only enter via a hand-authored VCF" — the multi-base half is false | Struck "multi-base and". Sentence now reads "combinatorial variants can only enter via a hand-authored VCF". | Parametric deletions of arbitrary span, `inframe`, `ala`, `stop`, `aa` and `snvre` all emit multi-nucleotide events (README:272–463; `mutators/`). The record lists `2del1`, `inframe`, `ala`, `stop` in the same paragraph — an internal contradiction now removed. The PoolParty contrast (single-event, no sampling, no pairwise mutator) is preserved intact. |
| A4 | Not every execution writes `config.json` | Both places → "each **library-generating** execution", plus "`--sequences-only` runs exit before that step". | `sge_proc.run_sge`: `if sequences_only: sys.exit(0)` occurs before `finalise(config, stats)`, and `common_cli.finalise` is what writes `OUTPUT_CONFIG_FILE_NAME`. |
| A5 | v4 does not consume GTF `stop` feature rows | Short clause added to the existing paper quote: "(§2.4.2 — paper-era: v4's loader recognises only `CDS` and `UTR` rows and ignores all other feature types; README and `loaders/gtf.py`)". The paper quote itself is accurate and was left verbatim. | `loaders/gtf.py`: `GtfFeatureType` has only `CDS` and `UTR`; the `match` statement drops everything else; stop context comes from `cds_features_to_exons` extending the terminal CDS by 3 bases. README:106: "Any features of type other than `CDS` and `UTR` are ignored." Wiki *Input-files*:91 says the same. |
| A10 | PAM edits are not restricted to synonymous/non-coding | → "user-defined SNVs … The paper describes them as synonymous or non-coding, but that is design intent, not an input restriction: v4's parser (`pam_variant.py`) requires only a one-base REF/ALT and an `SGRNA` tag, which is why `pam_mut_annot` can report `mis` and `non` as well." | `pam_variant.py:47`: `assert self.ref_len == 1 and self.alt_len == 1` plus a required `SGRNA` INFO tag — no consequence check. README field 26 allows `syn\|mis\|non\|ncd`. The record already listed all four codes two lines later. |

### From the fact-check, section B (added as brief clauses in existing entries)

| # | Missing fact added | Where | How verified |
|---|---|---|---|
| B4 | v4 reads only the **first ALT** of a multiallelic record and skips monomorphic records | §7, existing "Only simple custom variants" bullet | `custom_variant.py:97` `DnaStr(r.alts[0].upper())` with `# TODO: handle multiple ALT's?`; `_vcf_filter` drops `MONOMORPHIC` with a `logging.info`. Material because §6.10 sells ClinVar/gnomAD ingestion. |
| B7 | cDNA mode also lacks `--pam`, `--bg`/`--bg-mask`, `--revcomp-minus-strand`, `--sequences-only` | §7, existing cDNA bullet | `cdna_cli.py:26–28`: `@common_params` plus only `--annot`. |
| B8 | A mutable region overlapping **more than one CDS exon** is rejected | §7, existing coding/non-coding bullet | `targeton.py:96–103` `get_targeton_region_exon_id` → `InvalidTargetonRegion("overlaps multiple exons")`. |
| B10 | Length filtering measures the oligo **including adaptors**, appended after any reverse complement | `synthesis_constraints` entry | `meta_table.py:456–457` / `498–499`: `reverse_complement(oligo_no_adapt)` then `get_full_oligo(...)`, then `oligo_length = len(oligo)` fed to `info.eval_in_range`. |

---

## REJECTED (3)

| # | Finding | Why rejected | Evidence |
|---|---|---|---|
| C24 | "`span > 0` the only constraint" omits the `^(\d+)del(\d*)$` parser constraint | The record's claim is explicitly scoped to the named file, and within that file it is exactly right. The conclusion it supports ("arbitrary tiled deletion scans"; `3del1`, `6del0`, `10del2` all work) is unaffected — the regex admits any non-negative decimal span/offset, which is what "any valid value" means. Applying this would replace a true statement with a longer one. | `int_pattern_builder.py:__post_init__` contains only `if self.span <= 0: raise ValueError`. `loaders/mutator_config.py` regex admits `\d+` / `\d*`. |
| C30 | "authored by hand" / "hand-built" are uncited inferences | The claims are correct and I verified them directly rather than by inference. Leaving the record's wording is right; weakening it would make it less accurate. | `cdna_targeton_plan.txt` carries Excel thousands separators (`"1,402"`, `"1,593"`) and a hand-labelled `targeton` column with a duplicated `9.1` where `9.10` belongs — spreadsheet artefacts. `BRCA1_cdna_pCW57.1_model.gb` header reads `DEFINITION Ligation of 2 sequences` / `UNA 17-MAY-2021` with `/note="Geneious type: Editing History Insertion"` features. VaLiAnT emits neither TSV-plan nor GenBank output. (Strictly, the `.gb` was built interactively **in Geneious**, not literally by hand — a nuance not worth an edit.) |
| C31 | Cross-tool `seqpro` / `tangermeme` claims uncited | Not uncited-and-wrong — I checked them against the sibling records and they are accurate. No correction is warranted; adding cross-references is a formatting preference for the survey editor. | `final/seqpro.md:323` = "**`interface` — yes (Python API only).**"; `final/seqpro.md:330` = "**`license` — yes (MIT).**"; `final/tangermeme.md:603` = "**`license` → yes** — **RENDER THIS CELL AS: 'MIT'**". |

---

## ESCALATED (12) — record left untouched at these points

### E1 — `transcript_models`: "the code confirms the wiki" / "one transcript per run" (fact-check A6)

**Verified:** `loaders/gtf.py:131–155` raises `"Multiple transcript in annotation: not supported!"` only
for the `(gene_id, transcript_id)` pairs found **within the selected contig+strand group** —
`GtfLoader(contig, strand)` is constructed fresh per group inside `proc_contig`. `run_sge` only
*warns* ("Annotation provided for targetons on multiple contigs: verify only one transcript is used!")
when annotation accompanies multiple contigs. So different transcripts on different contigs can
proceed in a single invocation; the enforced unit is one transcript per contig/strand, not per
execution.

**Why escalated, not fixed:** the value's own qualifier is "**yes (one transcript per execution)**",
and correcting the evidence changes what that qualifier rests on. The wiki and paper do state the
one-transcript-per-execution rule as documentation, so the qualifier is not simply wrong — it is
documented-but-not-enforced. Per the no-value-changes rule this is yours to decide.

**Question:** should the qualifier read "one transcript per execution (documented; the code enforces
one per contig/strand and only warns otherwise)", or should the value text change?
**Options:** (a) keep the value, narrow the evidence sentence to the contig/strand scope;
(b) keep both as-is and accept that the code claim overstates; (c) change the qualifier.

### E2 — Background-overlap discard is indel-gated (fact-check A8)

**Verified:** `genomic_position_offsets.py:_compute_ref_del_mask` and `_compute_ref_offsets` populate
`_shift_mask` / `del_mask` / `pos_offset` **only** where `variant.alt_ref_delta != 0`. Both discard
paths (`ref_var_overlaps_var` for custom variants, `alt_var_overlaps_var` for pattern variants) read
those masks, and `alt_var_overlaps_var` returns `False` immediately for any variant with
`ref_len <= 1`. A generated mutation coincident with a background **SNV** is therefore not discarded.

**Why escalated:** the record's §6.5 sentence is a *verbatim README quote*, and the README states the
rule without the indel qualifier. Refining it means the record would contradict the tool's own
documentation on a code-derived nuance — in a document whose likely referees are the tool's authors.
That is an editorial call, not a repair.
**Options:** (a) leave the README quote and add a code-derived footnote; (b) leave entirely as-is;
(c) state the code behaviour in the record's own voice and flag the docs as imprecise.

### E3 — Installability / "pulls and runs today" claims (citation audit C32)

**Verified:** the record states "No install, no clone, no execution — document/web research only"
(§1), yet asserts "**Containers pull and run today**", "it should still build", "**installable and
runnable today**", "still containerised and installable", and speculates that `pysam==0.22.0` "may
need a source build on very new interpreters". The GitHub and quay.io APIs establish that the tags and
images **exist** with the stated timestamps (I re-confirmed all six quay tags and their dates); they
establish nothing about executability.

**Why escalated:** five separate locations across §3, §4 (`maintained`) and §9, and the claims are
part of what distinguishes `maintained = partial` from `no`. Fixing them is a hedging pass over more
than one paragraph and touches a capability value's basis.
**Options:** (a) soften to existence claims ("images published, both tags dated 2024-04-22; not
verified by execution"); (b) actually pull and run the container once and cite that; (c) keep the
claims and add an explicit caveat that they are inferred from artefact availability.

### E4 — Stale BRCA1 unique counts vs. v4.0.0 outputs (citation audit C22)

**Verified both sides.** Paper Supp. Table 3: nucleotide exons 2–5 = 583 / 740 / 825 / **1185**;
peptide = 340 / 500 / 580 / **918**. Recount of the committed v4 `_unique.csv` files:
nucleotide = 583 / 740 / 825 / **1217**; peptide = 340 / 500 / 580 / **919**. Exon 2 currently has
1000 in-range metadata rows **plus** 1 excluded row (1001 generated), whereas paper Table 2 labels
1000 as total and 1 as excluded. The cDNA figures do *not* drift: 40 files, 858–2092 per targeton,
sum exactly 62 754.

**Why escalated:** the record's numbers are correctly attributed to the paper, so nothing is
mis-cited. The question is what a record declared to describe **VaLiAnT 4.0.0** should print — an
authorial decision, not a factual one.
**Options:** (a) keep the paper numbers and add the v4 recount in parentheses; (b) replace with the v4
recount and cite the shipped outputs; (c) keep as-is since the attribution is explicit.

### E5 — Fact-check section B, remaining 7 items (B1, B2, B3, B5, B6, B9, B11)

Grouped because they share one question. Each asserts a real, verifiable omission; each is also a
materiality judgment, and together they would require substantial new prose (the record has no
mutator-semantics section to attach B1–B3 to).

- **B1** core `snv` semantics (all three non-reference bases at every eligible position, one per oligo)
- **B2** `inframe` / `ala` / `stop` semantics and the no-op skip when a codon already equals the
  selected replacement (`mutators/codon.py:45–88` — verified; it changes exact counts)
- **B3** full SNVRE behaviour (synonymous expansion, nonsense handling, exclusion of
  reference/triggering/duplicate codons, no extra variant where no additional synonymous codon exists)
- **B5** `--force-bg-indels` requires `--force-bg-ns`; non-coding background indels warn; PAM/background
  overlap in one coding codon is fatal
- **B6** duplicate-position PAM edits and multiple edits in one coding codon are rejected; a linked
  edit outside the targeton warns
- **B9** output directory must pre-exist; `.fai` required; contig names must match; `SPECIES`/`ASSEMBLY`
  are metadata strings only
- **B11** duplicate mutator labels and whitespace are ignored; a version-mismatched replayed config
  warns rather than guaranteeing compatibility

**Question:** which of these belong in a capability-comparison record, versus being tool-manual detail?
Adding all seven would lengthen §7 by roughly half and shift the record's altitude.
**Options:** (a) add B2 and B3 only (they change library-size arithmetic, which is the comparison axis
against PoolParty); (b) add all seven to §7; (c) add none — the record's job is capability presence,
not parameter semantics.

### E6 — Fact-check section C (balance)

The audit judges the record "mostly fair, with a capability rather than limitation emphasis", the main
imbalance being that absence claims get exhaustive treatment while ordinary `snv` / `inframe` / `ala` /
`stop` / `snvre` semantics are name-only. This is explicitly an emphasis judgment reserved for the
authors; no edit made. It is the same underlying question as E5.

---

## Tensions left in place deliberately

Per the "leave the nearby sentence alone" rule, these were **not** smoothed over:

1. **`lazy_evaluation` *Source:* line still cites "paper §2.5.3, Table 2"** although the Table 2
   sentence it supported was deleted (C8). Left as-is rather than editing an unaudited line.
2. **§10's `lazy_evaluation` bullet** now partly overlaps the reworded "seven of the eight" bullet
   above it (C3). Both are accurate; the redundancy is cosmetic.
3. **§5's preamble says "All under …/examples"** but item 6 is a wiki page, not an example directory,
   and has no expected-output dir. The item labels itself as a wiki page, so the reader is not misled;
   fixing it would mean rewriting the preamble.
4. **`mixed_mutagenesis_one_pool`'s caveat block is marked "keep verbatim"** in the record, yet three
   corrections landed inside it (A2, C10/A3, B1). All were factual omissions or errors, so the
   instruction was overridden — flagging it here so the authors can re-approve the new wording.
5. **§2's mutator code block is still edited code**, now honestly labelled "condensed" rather than
   replaced with true verbatim source (C1). If you would rather it be genuinely verbatim, that is a
   one-block swap.

## Errors observed but NOT in either audit — not applied, reported only

Found while verifying; outside the repair mandate, so the record is untouched at these points.

1. **Wrong section number.** The quote "targetons need to be processed individually as this function
   appends to all sequences in the library" is cited as **§2.4** in `combinatorial_motif_place`. It is
   at paper line 304, inside **§2.1** (lines 119–321); §2.4 begins at line 358.
2. **README section names that do not exist.** The record cites README "§*Variant call files*"
   (twice), "§*Adaptor sequences*", "§*Reference sequence*", "§*Background variants*" and
   "§*Input files*". The actual headings are `### Variant file`, and for the others the material lives
   under `## Command line interface` / `## File formats` / `## Usage` with no such heading. The
   underlying quotes are all genuine and were verified; only the section labels are wrong.
3. **Wiki staleness beyond the Docker page.** `Software-and-dependencies` states "Python 3.7+ is
   required"; `pyproject.toml` requires `>=3.11`. The record flags only *Docker-usage-example* as
   stale.
4. **README INFO-tag descriptions are swapped upstream** (`SGE_VCF_ALIAS` described as "VCF variant
   identifier", `SGE_VCF_VAR_ID` as "VCF variant source file alias"). A VaLiAnT documentation bug, not
   a record error; noted in case a referee raises it.

## Re-verified and found correct (no change needed)

Spot-checked because adjacent findings touched them: the 32-column metadata header (verbatim match);
the `oligo_name` example; the `aa` 19/20 quote; the lexicographic-dedup quote; the "Ambiguous
nucleotides…" quote; the background-variant scoping quote; the `force-bg-ns`/`force-bg-indels` quote;
all four CHANGELOG 4.0.0 quotes and the CHANGELOG 3.0.0 two-VCF entry; all BRCA1 and prime-editing
action vectors, coordinates and extension vectors; `run_brca1_nuc.sh` passing no `--max-length`;
`run_prime_a.sh` passing `--max-length 250`; the NXRL config's absent `maskBackgroundFilePath`;
`pyproject.toml` (one console script, no entry-points, four pinned deps, `requires-python >=3.11`);
README = 875 lines; `sge_proc.py` = 568 lines; recursive tree = 461 entries; GitHub metadata (not
archived, 6 stars, 3 forks, 2 open issues, `pushed_at` 2024-04-22T13:52:12Z, **0** Release objects);
all six quay.io tags and dates; PyPI = Duncan Dickinson's unrelated v0.2.3; anaconda 404 JSON.

## Pass 2 — policy application

**Date:** 2026-08-14

**Outcome counts:** 7 applied · 5 declined-by-policy · 0 rejected-unverifiable · 0 escalated.

**Version/source check:** the `develop` head remains
`8796cc112dafd4919fec59913f58cd2be87c45eb`; annotated tag `4.0.0` peels to
`822153257712d7732ec6fb948f984aecb5ccde4a`. Both commits have tree
`b40d87b53388f64e1c2f3d5dc8bb17c358232c94`, so the primary code and committed outputs checked on
`develop` are exactly the assessed v4.0.0 tree. Paper checks used PyMuPDF text extracted directly from
`papers/Barbon2022_VaLiAnT_all.pdf`.

| Open item | Outcome | Edit and verification |
|---|---|---|
| E1 / A6 — transcript scope and value basis | **applied** | Kept the capability value `yes` and the documented one-transcript-per-execution qualifier, but removed the false claim that code globally confirms it. The entry and §7 now say that `GtfLoader(contig, strand)` rejects multiple gene/transcript pairs within that group, while `run_sge` loops over contigs and only warns for annotated multi-contig input. Verified in the wiki *Input-files* page, `loaders/gtf.py:124–155`, and `sge_proc.py:409–482, 548–563`. The corrected evidence still supports locked row 17: the tool consumes GTF/GFF transcript/exon structure and makes design respect it. |
| E2 / A8 — background-overlap discard | **applied** | Kept the accurate README quotation and added the v4 implementation qualifier in §6.5 and §7: the discard masks are populated only for `alt_ref_delta != 0`, so the test covers coordinate-shifting background variants (indels), not coincident background SNVs. Verified in `genomic_position_offsets.py:108–146, 246–277` and both call sites in `targeton.py:204–227, 256–296`. |
| E3 / C32 — untested install/runtime claims and `maintained` basis | **applied** | Dropped “pull and run today,” “should still build,” the `pysam` source-build speculation, “runnable examples,” and “installable and runnable today.” Replaced them only with verified statements that source and container images remain published and that the README documents a source-install path. GitHub's primary API reports `archived=false`, `pushed_at=2024-04-22T13:52:12Z`, six stars, three forks, two open items, and zero Release objects; the Quay API lists `latest` and `4.0.0`, both dated 2024-04-22. `maintained = partial` remains supported by the dated inactivity, non-archived repository and published artifacts; no value changed. “Last commit” and `pushed_at` remain explicitly separate definitions. |
| E4 / C22 — paper/v4 BRCA1 count drift | **applied** | Added one parenthetical to the existing BRCA1 entry, as Policy D requires. Paper Supp. Table 3 gives nucleotide 583/740/825/1185 and peptide 340/500/580/918; direct row counts of the committed v4.0.0 `_unique.csv` files are 583/740/825/1217 and 340/500/580/919. The committed exon-2 metadata files contain 1000 in-range rows plus 1 excluded row. The paper figures remain attributed rather than replaced. |
| E5 / B1 — core `snv` semantics | **applied** | Added a short clause to `mixed_mutagenesis_one_pool`: all three non-reference bases are emitted at every eligible position, one event per oligo. Verified independently in README lines 294–309, `mutators/snv.py:44–49`, and paper Fig. 2. This directly bears on locked row 3's every-substitution-at-every-position test. |
| E5 / B2 — `inframe`/`ala`/`stop` semantics | **applied** | Added brief clauses to the existing codon entry: `inframe` deletes complete reading-frame triplets, while `ala`/`stop` select the top-ranked codon and skip a codon already equal to it. Verified in README lines 311–350, `mutators/deletion.py:51–59`, and `mutators/codon.py:41–88`. This bears directly on locked row 7's codon-level policy test. |
| E5 / B3 — full SNVRE behavior | **applied** | Expanded the existing SNVRE clause, without creating a section: synonymous expansion, analogous nonsense handling, exclusion of reference/triggering/duplicate products, and no design when no additional synonymous codon exists. Verified in README lines 383–463 and `mutators/snv_re.py:35–80`. This bears directly on locked row 7's replacement-policy test. |
| E5 / B5 — background-mask/validation edges | **declined-by-policy** | Verified: the BED mask filters positions (`sge_proc.py:131–168`); `--force-bg-indels` requires `--force-bg-ns` (`sge_config.py:63–69`); non-coding background indels warn and coding PAM/background overlap is fatal (`sge_proc.py:109–128, 313–320`). These rejection/failure edges do not change any locked-row score: row 11 explicitly excludes rejection-only checking, and row 6's pooled mixed-variant result is unchanged. They also do not alter a Table 1 cell. |
| E5 / B6 — PAM-edit conflict edges | **declined-by-policy** | Verified in `sge_proc.py:214–260`: duplicate positions and multiple edits in a coding codon fail; an edit outside the targeton warns. These rejection rules do not change a locked operational score or Table 1 Purpose, Key features, Output or Availability. |
| E5 / B9 — operational input contract | **declined-by-policy** | Verified in README lines 78–106, 211–217 and 630–638 plus `common_cli.py:32–46`: local FASTA/FAI and a pre-existing output directory are required, names must match, and `SPECIES`/`ASSEMBLY` populate metadata. None changes row 16's coordinate-input result or Table 1's “source, Docker” availability and coordinate-driven key-feature cells. |
| E5 / B11 — parser/replay restrictions | **declined-by-policy** | Verified in README lines 519–536 and 630–636: version mismatch warns; duplicate mutator labels and spacing are ignored. These parser details do not change any locked row score or Table 1 cell. |
| E6 / section C — balance/emphasis | **declined-by-policy** | Policy A declines all balance and proportional-emphasis edits to working records. No record text changed for this item. |

**Remaining escalations:** none. Both value-basis checks (E1 and E3) remain supported after their
evidence was corrected.

**Rejected-unverifiable findings:** none. Every factual open item was verified in an admitted primary
source; the five non-applied items were declined solely by the shipping-row/Table 1 policy or Policy A.

**Row-substitution candidates:** none observed. No verified VaLiAnT behavior here shows a locked row
to be uniform across tools, and these findings provide no cross-tool evidence sufficient to nominate a
better-discriminating replacement.
