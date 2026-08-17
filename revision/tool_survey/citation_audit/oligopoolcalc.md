# Citation-integrity audit — `oligopoolcalc`

Record audited: `final/oligopoolcalc.md`. Only evidence items that are not `verified` are reported below. Capability ratings and evaluative comparisons were not assessed.

## NOT FOUND IN ANY SOURCE

1. **NOT FOUND IN ANY SOURCE — record L195 and L294, quotation attributed to `https://oligopool.com`.** The quoted wording, *"Free, professional tools for oligonucleotide design and analysis… peer-reviewed methods, including SantaLucia 1998 and Breslauer 1986"*, does not occur verbatim on the live page or in an exact-phrase web search, and no retained 2026-08-10 page snapshot is cited. The current page instead says *"Oligo Pool Design & QC Before Vendor Upload"* and separately says that Tm calculations reference SantaLucia 1998 and Breslauer 1986. The live page does support the separate claims that it has 19 tools and is unrelated to the Hossain/Salis package, but it does not support this purported quotation.

## Other non-verified evidence

1. **misquoted — record L6, paper title.** The cited PDF/DOI `10.1021/acssynbio.4c00661` is titled *"Automated Design of Oligopools and Rapid Analysis of Massively Parallel Barcoded Measurements."* The record substitutes *"Automated Design of Thousands of Nonrepetitive Parts for Engineering Stable Genetic Systems,"* which is a different 2020 *Nature Biotechnology* paper about the Non-Repetitive Parts Calculator (`10.1038/s41587-020-0584-2`). The authors, journal, year, volume, pages, DOI, PMID, and PMCID that follow in L6 belong to the Oligopool Calculator paper.

2. **minor-discrepancy — record L151, “Paper title/keywords include ‘massively parallel reporter assay’.”** That exact phrase is a paper keyword, but it is not in the paper title; the title says *"Massively Parallel Barcoded Measurements."*

3. **misquoted — record L195, paper Data Availability.** The record says Data Availability “lists only the two GitHub repos.” The PDF's Data Availability Statement lists three locations: `(1) github.com/hsalis/SalisLabCode`, `(2) github.com/ayaanhossain/oligopool`, and `(3) pypi.org/project/oligopool`.

4. **wrong-location — record L199, GPL quotation.** The GPL sentence is real, but it is not in the paper's Data Availability Statement. It starts at the end of the Introduction on printed p.4219 and finishes on p.4220. The actual Data Availability Statement is near the end of the paper and does not include the license wording.

5. **misquoted — record L213, 6.8-second benchmark.** The paper says, *"designing the spacer regions for 93,180 oligonucleotides required only 6.8 s (Table S1)."* The record calls this a “Verify/design step,” which changes the operation being timed.

6. **number-wrong — record L29, “13 missed capabilities added; see §6.”** The original extraction's notable-capability list has 11 items, and final §6 has 22. Final items 12–22 are 11 newly enumerated items, not 13. The final section itself likewise contains no group of 13 additions.

7. **misquoted — record L94–L108 and repeated references to “all 22 module signatures.”** The cited API does not say that every public item has signature `(input_data, ..., output_file) -> (DataFrame, stats)` or returns a materialized DataFrame/CSV. `docs/docs.md` L140–L156 explicitly distinguishes transform functions, five stats-only functions (`background`, `lenstat`, `inspect`, `index`, `pack`), and two counting functions. The 22 public names are 20 functions plus two advanced classes, `vectorDB` and `Scry`; the classes do not have module-style DataFrame signatures. Several functions also lack `input_data` or `output_file` under those names.

8. **minor-discrepancy — record L105–L108, lazy-evaluation evidence.** The public transform/count results are eager, but “the package evaluation model is eager throughout” overstates what the cited checks establish. Package source contains generators/streaming in `base/core_count.py`, `base/core_pack.py`, `base/vectordb.py`, `base/scry.py`, `base/utils.py`, and `base/core_barcode.py`; `oligopool/__init__.py` also lazily imports public attributes. None of that creates a user-facing lazy DataFrame API, but it contradicts the implementation-wide wording. In addition, the stated documentation grep has three literal `yield` hits (`agent-skills.md:247`, `api.md:955`, `api.md:1257`), not “only chunked/downstream prose”; they are ordinary-language, non-feature hits.

9. **misquoted — record L66, “An Oligopool Calculator library is therefore always exactly the rows the user supplied.”** The `join` and `merge` quotations are accurate, and no row-union operator was found, but the inference is false across the cited inventory. `compress` creates a smaller table of new synthesis rows, while `expand` can create exponentially more rows. The record itself acknowledges both behaviors at L106, L136, and L270.

10. **minor-discrepancy — record L111, “none of the 22 `op.*` modules generates variants.”** The narrower claim that the package has no mutation-specification module is supported. The blanket wording is not: `op.expand` enumerates every concrete sequence represented by IUPAC-degenerate input and can multiply rows. It does not infer a mutagenesis design, which is the supported narrower distinction.

11. **wrong-location — record L133, “The `ID` column is a required user input on every module (docs.md L1304).”** `docs.md` L1304 is one workflow instruction, *"Start with your variants in a CSV (must have an `ID` column),"* not a package-wide signature rule. Counterexamples include `background` (sequence list/FASTA/`Sequence` table), `inspect` (artifact path), `pack` (FASTQ/index artifacts), and the `vectorDB`/`Scry` classes.

12. **number-wrong — record L140/L141/L317, “all 40 `oligopool/**/*.py` files.”** At release tag `v2026.02.22.1` / commit `b88fa394ca67ed4c48ec41127b5707694ee7cf0a`, `find oligopool -type f -name '*.py'` returns 41 files (23 top-level and 18 nested). The zero-match plotting-import result is correct; the file count is not.

13. **number-wrong — record L26/L123/L173/L259, `primer_guanine_cytosine_content` “6 occurrences.”** The literal token occurs 7 times in `oligopool/primer.py`. It appears on six source lines because line 1245 contains it twice; the other lines are 542, 1070, 1176, 1321, and 1347.

14. **misquoted — record L183, “`oligo_length_limit` is a required parameter on every design module.”** Under the record's own Design inventory, `background`, `merge`, `revcomp`, and `final` do not accept `oligo_length_limit`; `barcode`, `primer`, `motif`, `spacer`, and `join` do. The cited signatures therefore do not support “every.”

15. **minor-discrepancy — record L185, “It caps homopolymer runs only incidentally, at ~k+1.”** The cited source establishes only that `maximum_repeat_length=k` rejects a candidate `(k+1)`-mer when that k-mer is present in the input/context/background corpus. It does not intrinsically cap a homopolymer run in a candidate. Without a matching corpus k-mer or an explicit homopolymer in `excluded_motifs`, such a run is not prohibited by this mechanism.

16. **minor-discrepancy — record L128/L188/L258, causal interpretation of `verify`.** Source computes each motif's baseline as the minimum absolute count across completed oligos and flags a row when its count exceeds that minimum. An excess can occur wholly inside one variant column and need not have emerged from concatenation; conversely, a junction motif present equally in every row is baseline and is not flagged. `columns` is assigned from each match's start offset, so saying it identifies the element that “produced” a junction-spanning motif is also stronger than the implementation supports. The quoted schema fields themselves are real.

17. **number-wrong — record L219, notebook “re-executed 2026-02-18.”** The notebook markdown says *"Updated February 18, 2026,"* but the current file's latest Git commit is `b88fa394...`, dated 2026-02-22T02:46:06Z with subject `chore: rerun notebook`. “Updated” is not evidence of the claimed execution date, and repository history gives the later rerun date.

18. **number-wrong — record L219, “3–17-bp spacers.”** Early notebook prose states that target, but the executed pipeline puts three degenerate bases into the separate 9-bp `AatIIPadded` motif and then creates a `Spacer` column of 0–14 bp. The executed `lenstat` output says `Element Spacer: Occupies 0 to 14 Base Pair(s)`. A 3–17-bp range describes motif-flank-plus-spacer padding, not the actual spacer column.

19. **number-wrong — record L220/L246, parser over “4,351 real promoter sequences.”** `promoters.txt` does contain 4,351 lines, but `run_design_assembly_parser.py` filters out any sequence containing eight excluded motifs before calling `design_parser`. Recomputing the code leaves 3,503 sequences and removes 848. Thus 4,351 is the raw-file count, not the parser invocation's pool size.

20. **minor-discrepancy — record L220/L246, `get_primer_order()` “topologically resolves the paired-primer dependency graph.”** The cited function is a one-pass `deque` routine: for each primer, it inserts an unseen paired primer immediately before it. It has no graph construction, indegree calculation, cycle detection, or general topological-sort algorithm. It produces a usable order for the shipped reciprocal-pair example, but the source does not support the broader algorithm description.

21. **minor-discrepancy — record L146/L147/L230/L240, quantitative “6–20×” compression claim.** The `6-20x` text genuinely occurs in `docs.md` and `agent-skills.md`, but the cited saturation-mutagenesis template uses placeholder input and supplies no counts from which the range can be recomputed. Shipped executed examples do not demonstrate that range: the notebook shows 12→10 (1.2×) and 500→360 (1.4×), while `run_degenerate_demo.py` states theoretical/example outcomes of 64×, 4,096×, and 1×. The record can verify that the authors document the range, but not present 6–20× as a rechecked empirical result.

22. **uncited — record L77/L136/L166, “HGVS-style” IDs.** `mutant_generator.py` verifies the format `f'{ref_base}{pos + 1}{alt_base}'` (for example `A5T`), but neither source nor documentation calls it HGVS. The official HGVS DNA-substitution syntax includes a coordinate-type prefix and `>` (for example `c.5A>T`; see `https://hgvs-nomenclature.org/21.1.0/recommendations/DNA/substitution/`), so the “HGVS-style” characterization has no cited support in the record.

23. **minor-discrepancy — record L41/L288, `hsalis/SalisLabCode` described as a “mirror.”** The URL is live and the paper lists it. The repository is a multi-project Salis Lab monorepo containing an `Oligopool_Calculator` subdirectory, not a repository mirror of `ayaanhossain/oligopool`. No evidence is supplied that its contents are identical to, or current with, release `v2026.02.22.1`.

24. **uncited — record L193/L294, “no hosted web service exists.”** Repository and paper sources verify that no hosted service is claimed and that `oligopool.com` is unrelated. They do not establish the stronger universal existence claim; no exhaustive service search or author statement is cited.

25. **uncited — record L203/L217/L292/L296, “runnable” / “installable and runnable today.”** PyPI artifacts, install commands, Docker files, examples, and version metadata exist, but no installation or example execution is cited. L289 says installation was not attempted, and L319 explicitly says runtime behavior and current-toolchain installability are unverified. Availability of artifacts does not verify the stronger runtime claim.

26. **uncited — record L195/L294, “A web-search summariser misidentified it.”** No search result, transcript, URL, or retained artifact is cited for this process claim. The unrelatedness of the live site is independently verifiable; the alleged summarizer error is not.

27. **uncited — record L8/L14/L47/L302/L321, independent-pass provenance and rate-limit narrative.** The current files and live metadata support many underlying facts, but the record supplies no logs or retained responses proving a third independent fetch, a rate-limit event, or two prior independent API fetches that agreed exactly. These are unsupported process claims rather than source-backed evidence.
