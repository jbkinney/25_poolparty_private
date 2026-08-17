# Citation-integrity audit — VaLiAnT

Record audited: `final/valiant.md`

Reference source checked: VaLiAnT tag `4.0.0` (`822153257712d7732ec6fb948f984aecb5ccde4a`) and current `develop` head (`8796cc112dafd4919fec59913f58cd2be87c45eb`), the live GitHub and Quay APIs, the live wiki, PyPI and Anaconda APIs, and all 11 pages of `papers/Barbon2022_VaLiAnT_all.pdf` extracted with PyMuPDF. This audit concerns citation accuracy only.

## NOT FOUND IN ANY SOURCE

None. Several passages presented as quotations are not verbatim, but each has an identifiable underlying source; they are classified `misquoted` below rather than `NOT FOUND IN ANY SOURCE`.

## Findings (every item not verified)

1. **Status: misquoted. Record lines 74–89.** The code block is labelled “`src/valiant/mutator_type.py`, verbatim,” but it is not verbatim. The source puts `# Parametric span`, `# Single-nucleotide span`, and `# Codon span` on separate preceding lines; it does not contain the invented inline comments `# parametric span: <SPAN>del[<OFFSET>]` or `# codon span, CDS only`; the sets have type annotations; and `CDS_MUTATOR_TYPES` contains qualified members such as `MutatorType.SNV_RE`, not bare `SNV_RE`. The seven enum values and set membership are substantively accurate, but the displayed quotation is edited code.

2. **Status: number-wrong. Record lines 13, 218, 305, 459, 704, and 714.** The GitHub wiki has **11**, not 8, Markdown content pages: `Advanced-usage`, `Configuration`, `Docker-usage-example`, `Home`, `Input-files`, `Output-files`, `Saturation-prime-editing-example`, `Saturation-prime-editing`, `Software-and-dependencies`, `cDNA-DMS-file-formats`, and `cDNA-example`. The source table lists only eight consulted pages and omits `Configuration`, `Saturation-prime-editing`, and `Software-and-dependencies`; therefore the repeated claims that “all 8 wiki pages” were checked are also false.

3. **Status: number-wrong. Record lines 27 and 711–718.** The final record encodes **eight** `no` values, not seven: the seven listed at lines 711–712 plus `lazy_evaluation = no`. Treating `lazy_evaluation` in a separate bullet does not make the claim “every one of the seven ‘no’ cells” numerically correct.

4. **Status: minor-discrepancy. Record lines 91–92.** `src/valiant/mutators/` contains six files, not only the five listed: it also contains `__init__.py`. The five named implementation modules do exist.

5. **Status: wrong-location. Record lines 94–97 and 153–154.** The claimed fixed order “mutators → custom VCF” is reversed in the cited v4 source. `Targeton.process()` calls `_process_custom_variants(...)` first and `_process_pattern_variants(...)` second (`src/valiant/targeton.py:300–312`). The relevant order is not wholly hard-coded in `sge_proc.py`; adaptor addition, length filtering, and deduplication are in `meta_table.py`. Background application precedes PAM protection as claimed.

6. **Status: number-wrong. Record line 139.** `MetaTable.to_csv(conn, targeton_name)` does not write “one CSV per targeton.” It names and writes `_meta.csv` and `_unique.csv` and conditionally retains `_meta_excluded.csv` (`src/valiant/meta_table.py:340–355, 440–447, 666–686`). Thus an ordinary targeton has two CSV outputs and a length-filtered targeton has three, in addition to the two SGE VCFs.

7. **Status: minor-discrepancy. Record lines 168–170.** The blanket claims “There is no generator, no streaming” are false. The source contains generator/context-manager functions in `db.py`, `loaders/csv.py`, `loaders/fasta.py`, `loaders/vcf.py`, `meta_table.py`, and `vcf_writer.py`; `loaders/csv.load_csv()` yields rows, and `MetaTable.to_csv()` fetches one SQLite row at a time and writes it immediately. The narrower claim that VaLiAnT exposes no lazy oligo-library API remains supported, but that is not what these two blanket statements say.

8. **Status: wrong-location. Record line 170.** Paper Table 2 and Supplementary Table 3 contain static, targeton/example-library counts; they do not show that “Full library counts are reported per run” by the software. In v4, `finalise()` logs only counts excluded by the length filter (`src/valiant/common_cli.py:94–109`), then writes `config.json`.

9. **Status: number-wrong. Record lines 187–188.** Paper Table 2 has **12** nonzero generated-mutator/region combinations: 3 in r1, 6 in r2, and 3 in r3. Alternatively it has 8 distinct generated mutator rows plus 2 custom-source rows. Neither reading supports “ten mutator/region combinations plus custom variants.” The totals 1000 and 583 do appear in the paper table.

10. **Status: minor-discrepancy. Record lines 192–193.** `--include-no-op-oligo` is not an unconditional WT-inclusion switch. The CLI says it includes a no-op oligo **if a background or PAM-protection variant is applied**, and `meta_table.py:422–425` enforces that condition. The emitted “no-op” sequence is the background/PAM-altered reference, so calling it WT without qualification is inaccurate. It is a single row when emitted.

11. **Status: number-wrong. Record line 213.** There are **87 files** (Git blobs) under `src/valiant/` at v4.0.0, not 93. The full recursive GitHub tree has 461 entries when directories are included, which is verified. The apparent “93” comes from counting directory/tree entries and/or the `src/valiant` directory itself as files.

12. **Status: minor-discrepancy. Record lines 218–220.** The source does have sequence-uniqueness machinery: `_unique.csv`, a `unique_oligos` map, sequence deduplication, and lexicographic representative-name selection. The evidence supports “no barcode generation or barcode uniqueness design,” but not the unqualified phrase “no ... uniqueness machinery anywhere in the source.”

13. **Status: misquoted. Record line 224.** Not every generated oligo gets a row in `_meta.csv`. Length-filtered oligos are written only to `_meta_excluded.csv`; the committed BRCA1 exon-2 example contains one such row. The 32-column header itself is verified verbatim and is used for both metadata files.

14. **Status: misquoted. Record lines 239–242, 368–370, and 591–593.** `SGE_VCF_ALIAS` and `SGE_VCF_VAR_ID` are not carried by every VCF record. `SGE_SRC` and `SGE_OLIGO` are universal; alias/ID are conditional and appear for custom variants when the relevant values exist (`src/valiant/vcf_writer.py:134–160`). The committed VCFs visibly omit alias/ID on ordinary `snv`, `del`, `aa`, etc. records.

15. **Status: misquoted. Record lines 368–370.** “All generated variants are emitted” is false under a length filter. VCF records are written only inside the `info.eval_in_range(...)` branch (`src/valiant/meta_table.py:585–625`). For example, the committed exon-2 `_meta_excluded.csv` oligo is absent from both output VCFs. The accurate statement is that all **in-range** generated and custom variants are emitted in SGE mode.

16. **Status: number-wrong. Record line 359.** `src/valiant/loaders/` contains **16 files**, not 17. Also, the literal assertion that every module there “is bed/csv/fasta/gtf/vcf/tsv” is false: the directory includes experiment/configuration, error, and utility modules such as `cdna_experiment_config.py`, `mutator_config.py`, `errors.py`, and `utils.py`. No HGVS parser was found, so the underlying no-HGVS-input evidence is otherwise supported.

17. **Status: misquoted. Record lines 402–403.** The quoted README text “DNA sequence to be added at the 5'/3' end of the oligonucleotide” does not occur. The README has two separate option descriptions, one saying “at the 5' end” and one saying “at the 3' end.” Combining them as `5'/3'` inside quotation marks is a paraphrase presented as a quotation.

18. **Status: number-wrong. Record lines 406 and 529–530.** The literal passed to `--adaptor-5` in both `run_prime_a.sh` and `run_prime_b.sh` is **121 nt**, not 120 nt. The `--adaptor-3` literals are 38 nt, as claimed.

19. **Status: misquoted. Record lines 406–407.** The sentence “The 21-base flanking regions at targeton boundaries remain unmutated to facilitate primer binding” is not verbatim in the cited cDNA wiki page. The page separately states that `targeton_start`/`targeton_end` are 21 bp outside r2, that those bases are not subjected to mutators, and that this allows primer binding. The record is an accurate paraphrase, but quotation marks incorrectly present it as source text.

20. **Status: wrong-location. Record lines 421–423 and source list 426–428.** The quotation fragment “allow[ing] for insights into the effect of codon sequence on missense changes” comes from paper §2.1, not from any of the cited README sections, code files, default table, or CHANGELOG entry listed for this paragraph. In addition, the current README describes SNVRE as choosing the top-ranking synonymous codon of the SNV product (or second-ranked if it is already top), which is more specific than the paragraph’s blanket “next most frequent triplet.”

21. **Status: number-wrong. Record lines 535–536.** Recomputing inclusive lengths from all 40 rows of `examples/cdna/input/cdna_targeton.tsv` gives targeton lengths **123–237 bp**, not 132–237 bp. Row 5 is `901..1023`, which is 123 bp. The paper and wiki themselves repeat the 132-bp typo. The r2 range 81–195 bp, 40-row count, common action vector, and 62,754 summed unique sequences are verified.

22. **Status: minor-discrepancy. Record lines 520–524.** The stated counts are genuinely present in the 2022 paper, but they are stale relative to the record’s declared v4.0.0 scope and the cited current shipped examples. Recounting current v4 `_unique.csv` files gives nucleotide exons 2–5 = **583 / 740 / 825 / 1217**, not 583 / 740 / 825 / 1185, and peptide exons 2–5 = **340 / 500 / 580 / 919**, not 340 / 500 / 580 / 918. Current exon-2 outputs contain 1000 in-range metadata rows plus 1 excluded row (1001 generated), whereas paper Table 2 labels 1000 as total and 1 as excluded. These are version-drift discrepancies, not invented paper citations.

23. **Status: number-wrong. Record line 556.** Wiki `VERSION=2.0.0` is **three later software releases/tags behind** 4.0.0 (`3.0.0`, `3.0.1`, `4.0.0`), not four releases behind. It is the second of five tagged releases.

24. **Status: minor-discrepancy. Record lines 568–570.** `span > 0` is not the only CLI-level constraint on parametric deletion parameters. `MutatorConfig` requires the form `^(\d+)del(\d*)$`, so span and explicit offset are nonnegative decimal integers, and `IntPatternBuilder` additionally rejects span 0. The examples `3del1`, `6del0`, and `10del2` parse, but the “only constraint” assertion omits the parser constraints.

25. **Status: minor-discrepancy. Record line 640.** “Only SNVs reach” partial/liminal codons is too absolute. It is true for generated codon-replacement mutators, but custom VCF variants are applied across the whole targeton and can affect partial codons. The cited paper sentence concerns generated CDS-specific mutators and SNVs, not all possible variant sources.

26. **Status: misquoted. Record line 652.** “Silently discarded (with warnings)” contradicts the cited README, which says warnings are **always raised**. The variants are discarded, but not silently.

27. **Status: misquoted. Record lines 479–480.** The quoted commit message is not verbatim. The Git commit message is `Merge tag '4.0.0' into develop\n\nBackground variant support`; the record substitutes an em dash between subject and body. The timestamp, author, and message content are otherwise verified.

28. **Status: number-wrong. Record lines 485–488 and 706–707.** “There is no human issue traffic at all; every closed issue was authored by the maintainers” is contradicted by GitHub issue **#8**, opened by external user `delagee` (`author_association: NONE`) on 2023-06-27, with a detailed human bug report and a maintainer response. The two currently open PRs are indeed automatically created Snyk PRs #10 and #11, and the other live repository counts are correct.

29. **Status: minor-discrepancy. Record lines 502–506 and 555–556.** Not every example uses filenames matching `run_*.sh` and `check_*.sh`: cDNA and NXRL use `run.sh` and `check.sh`. Nor is every expected directory literally named `output_exp/` (the BRCA1 examples use `brca1_nuc_output_exp`, `brca1_pep_output_exp`, etc.). The validation scripts do use `md5sum` with `md5` fallback against their corresponding committed expected-output directories. The blanket “Inputs must be unpacked first” does not apply to the cDNA example; it applies to the shared SGE reference and the prime-editing archives identified by the two scripts. The sixth listed item, Docker usage, is a wiki page rather than an example under the cited `examples/` tree and has no committed expected-output directory of its own.

30. **Status: uncited. Record lines 270–273, 540–543, and 624–626.** File existence proves that PNG, GenBank, and `cdna_targeton_plan.txt` assets are committed and that VaLiAnT does not emit those file types. It does not by itself prove the stronger authorship claims “authored by hand,” “hand-built,” or “the 40-targeton tiling plan was authored by hand.” Those are plausible inferences, but no cited source states them.

31. **Status: uncited. Record lines 112, 461–464, 472–474, and 671–675.** Claims about the `seqpro` record’s exact interface cell and the MIT licences of `seqpro` and `tangermeme` are cross-tool facts with no supporting citation in this VaLiAnT record. They were outside the VaLiAnT repository/documentation/paper source set and are not verified here.

32. **Status: uncited. Record lines 110, 488–491, and 694–706.** The cited GitHub and Quay APIs establish source/tag/image existence and timestamps, not that the package or containers “pull and run today,” “should still build,” or are “installable and runnable today.” The record explicitly says “No install, no clone, no execution” at line 56. The speculation that `pysam==0.22.0` may require a source build on new interpreters is also uncited. These runtime/installability claims require an actual installation/container execution or an external build result.

33. **Status: misquoted. Record lines 202–203 and 643–644.** The quoted sentence “Only simple variants such as the following are supported: substitutions, insertions, deletions, indels” does not appear verbatim. The README has the introductory sentence followed by a four-item Markdown list; the paper instead says “Simple variants are currently supported, including substitutions, insertions, deletions and indels.” The supported variant types are accurate, but the record compresses source prose and list items into a new sentence inside quotation marks.

34. **Status: minor-discrepancy. Record line 534.** The record says the in-silico digest used `SalI/NheI`, while both the cited paper text and the cDNA wiki page spell the second enzyme `NehI`. `NheI` is likely a sensible correction of the authors’ typo, but it does not match the stated source and the correction is not documented in the record.

## Link check

No dead links found among the explicit record links. The repository, wiki, `develop/examples`, and Quay repository resolve; the DOI resolves to the publisher (the publisher returns 403 to automated retrieval, but the DOI itself is live); the PyPI JSON endpoint is live and is the unrelated Duncan Dickinson package as claimed; and the Anaconda endpoint returns the claimed 404 JSON error.

## Verified counts and anchors relevant to the findings

- Paper PDF: 11 pages; Supplementary Tables S1–S3 and Figures S1–S2 are present.
- `README.md`: 875 lines; `src/valiant/sge_proc.py`: 568 lines.
- Recursive Git tree: 461 entries including directories; `src/valiant/`: 87 files; `src/valiant/loaders/`: 16 files.
- Metadata schema: 32 columns, matching the displayed header.
- Quay: six tags with the dates stated in the record.
- GitHub as checked 2026-08-14: default branch `develop`, not archived, 6 stars, 3 forks, 2 open PRs, last commit 2024-04-22T13:51:54Z, `pushed_at` 2024-04-22T13:52:12Z, and zero GitHub Release objects.
