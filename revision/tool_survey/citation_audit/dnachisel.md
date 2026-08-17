# DNA Chisel citation-integrity audit

Record audited: `revision/tool_survey/final/dnachisel.md`

Repository evidence was checked at `Edinburgh-Genome-Foundry/DnaChisel` commit `68c09304341c3656f3dfe63eda37757d6a7b3917` (the cited 2025-05-10 `master` commit). Only non-verified items are reported below.

## NOT FOUND IN ANY SOURCE

1. **Status: NOT FOUND IN ANY SOURCE — record line 58.** The claim that anything named `variant` or `library` is absent “from the entire source tree” is contradicted by that tree. `variant`/`variants` occurs throughout ten files under `dnachisel/`, including `MutationChoice.variants`, `MutationSpace.all_variants`, and both solver mixins; `library` occurs in three files under `dnachisel/`. This also conflicts with the record’s own later citations to `all_variants`.

2. **Status: NOT FOUND IN ANY SOURCE — record line 86.** “The only multi-sequence API is” `constraints_breaches_dataframe` is false. At minimum, `AvoidMatches(..., sequences=None)` accepts and indexes a set of sequences (`dnachisel/builtin_specifications/AvoidMatches.py:13-18,65-87`), `AvoidBlastMatches(..., sequences=None)` searches supplied subject sequences (`AvoidBlastMatches.py:45-88`), `AvoidHeterodimerization(other_primers_sequences)` consumes a list of other sequences, and `MotifPssmPattern.from_sequences(sequences)` consumes a sequence collection.

3. **Status: NOT FOUND IN ANY SOURCE — record lines 102 and 129.** The claims that `MotifPssmPattern` exists only for `AvoidPattern`-style constraints and that motif functionality is avoidance-only are contradicted by the cited sources. `docs/ref/core_classes.rst:54-55` says both `AvoidPattern` and `EnforcePatternOccurence` accept a Pattern. `EnforcePatternOccurence.__init__` accepts a `SequencePattern`, and `evaluate` calls `self.pattern.find_matches` (`dnachisel/builtin_specifications/EnforcePatternOccurence.py:54-67,88-94`). The official Supplementary Information also explicitly lists a PSSM as a pattern usable by “specifications such as AvoidPattern and EnforcePattern” (SI1, Note 1).

4. **Status: NOT FOUND IN ANY SOURCE — record line 133.** “Every objective is a biochemical or manufacturing score” is not supported and is contradicted by the API. `AvoidChanges` and `EnforceChanges` provide generic edit-count objectives, and the paper/source explicitly permit arbitrary user-defined `Specification` subclasses with custom evaluation functions. The cited finite list cannot support the absolute “every objective” claim.

5. **Status: NOT FOUND IN ANY SOURCE — record line 183.** PyPI JSON does not contain an OSI MIT classifier. At audit time, `info.classifiers` is an empty list. `info.license_expression = "MIT"` is real, as are the repository `LICENSE` and README sentence, but the claimed fourth piece of metadata is absent.

6. **Status: NOT FOUND IN ANY SOURCE — record line 267.** “Supplementary Sections S1A-S1E ... remain unobtainable” is false. The University of Edinburgh hosts a live 21-page author PDF containing the paper and the complete 16-page SI1; its contents page lists A through E, and the following pages contain those sections: <https://www.pure.ed.ac.uk/ws/portalfiles/portal/154025952/DnaChisel_with_SI.pdf>. Absence from the local `papers/` directory does not make the supplement unobtainable.

## Other non-verified evidence

1. **Status: misquoted — record line 18.** The meta-claim that “all quotes below [were] verified verbatim” is false. The non-verbatim quotations/signatures are itemized below.

2. **Status: minor-discrepancy — record line 44.** The purported exhaustive top-level listing after “`dnachisel/` contains only” omits `README.md`, `__init__.py`, and `version.py`, all present at the cited commit.

3. **Status: misquoted — record lines 46, 106, 224, and 260.** The re-verified `random_compatible_dna_sequence` signature is not reproduced verbatim. The source signature ends `max_random_iters=5000, logger="bar", **kwargs`; line 46 omits `**kwargs`, while lines 106 and 224 omit both `logger` and `**kwargs`. Line 260 nevertheless labels the signature “re-verified verbatim.”

4. **Status: number-wrong — record lines 73 and 86.** `dnachisel-dtailor-mode` was not “released once on PyPI.” PyPI has four releases—0.0.1, 0.0.2, 0.1.0, and 0.1.1—all uploaded on 2020-01-14. “All four releases ... on 2020-01-14” elsewhere in the record is verified; “released once” is the contradictory count.

5. **Status: wrong-location — record lines 86 and 217.** The record attributes “one `DnaOptimizationProblem` per CDS, results appended to a list and written as a multi-record FASTA” to both `ecoli_genes_optimization.py` and `genome-wide-optimization.py`. Only the latter does that (`genome-wide-optimization.py:9-32`). `ecoli_genes_optimization.py:23-40` accumulates per-gene constraints, constructs one problem over the whole GenBank record, and writes one GenBank record; it is not a per-CDS problem/output loop.

6. **Status: misquoted — record line 98.** The purported `changes.md` quotation adds a word. The cited text is exactly “AvoidChanges and EnforceChanges can now tunable,” not “can now **be** tunable” (`changes.md:19`).

7. **Status: misquoted — record line 102.** The second quoted documentation phrase silently corrects the source. `docs/genbank/genbank_api.rst:149-150` says new occurrences “will be **be** preferentially placed towards the center”; the record presents “will be preferentially placed ...” as a verbatim quotation. This is a harmless source typo correction, but it is not verbatim.

8. **Status: minor-discrepancy — record line 106.** The evidence cited for “Generate N short sequences de novo” only generates one sequence per call: `random_compatible_dna_sequence` returns the single value `problem.sequence` (`dnachisel/utils/utils.py:40-46`). None of the cited barcode sources implements an N-sequence generator or returns the claimed bare list; a user-written loop is required.

9. **Status: minor-discrepancy — record line 111.** `to_record(..., with_constraints=True)` does not write a feature for “every specification applied.” Constraints are filtered to objects with a truthy singular `location` (`RecordRepresentationMixin.py:158-168`); for example, `EnforceRegionsCompatibility` has `locations` and no `location`, so it is omitted.

10. **Status: number-wrong — record line 111.** It does not create “one feature per nucleotide edit.” `sequence_edits_as_features` calls `sequences_differences_segments`, which merges adjacent differing nucleotides into one `(start, end)` segment and makes one feature per contiguous segment (`RecordRepresentationMixin.py:194-206`; `biotools/sequences_differences.py:28-42`).

11. **Status: wrong-location — record lines 111, 155, 235, and 260.** The report filenames conflate two different report functions. A successful `write_optimization_report` writes `Report.pdf`, the two before/after CSVs, `final_sequence_with_edits.gb`, and `final_sequence.gb` (`optimization_reports.py:286-423`); it does **not** write `logs.txt` or `plots.pdf`. Those two files belong to `write_no_solution_report`, whose other principal artifact is `local_constraints_breaches.gb` (`optimization_reports.py:82-190`). No single report bundle contains the combined artifact list stated on line 155, and `logs.txt` is not part of the successful provenance bundle claimed on lines 111 and 235.

12. **Status: minor-discrepancy — record lines 111 and 235.** The source-copy/hash evidence is overstated. Both are emitted only when `file_path is not None`, and `file_hash` is only the first eight hexadecimal characters of the MD5 digest (`optimization_reports.py:349-357`), not an unconditional source copy with a full MD5 hash.

13. **Status: wrong-location — record line 119.** Current `write_optimization_report` does not render the sequence with annotated constraints and objectives “as PDF plots” through `dna_features_viewer`. Its PDF contains summary tables, sequenticons, and an optional GeneBlocks edit-diff figure. The whole-sequence and local annotated plots belong to `write_no_solution_report`. The cited paper describes the historical/high-level report, but the cited current implementation does not perform the behavior attributed to this function.

14. **Status: minor-discrepancy — record lines 119 and 246.** `max_features_in_plots` exists in the `write_optimization_report` signature and its quoted warning is present in the docstring, but the parameter is never read anywhere in the function after the docstring. It is not an operative “option” in the cited implementation.

15. **Status: wrong-location — record lines 70, 147, 244, and 260.** `Location.from_biopython_location` would indeed reduce any object directly passed to it to outer `start/end/strand`, but that is not the normal `EnforceTranslation`/record path claimed. `EnforceTranslation.__init__` calls `Location.from_data`; that accepts a Biopython `FeatureLocation` (`SimpleLocation`) but not a sibling `CompoundLocation`, so a directly supplied compound location falls through to `None` rather than being collapsed (`Location.py:163-186`; Biopython `Bio.SeqFeature` API). In addition, `DnaOptimizationProblem.from_record` ignores CDS features and parses only `misc_feature` specification labels (`RecordRepresentationMixin.py:63-74`). The cited helper supports a hypothetical explicit call, not the stated end-to-end claim that a spliced CDS is silently treated as its intron-containing outer span.

16. **Status: minor-discrepancy — record line 151.** The list called the “complete GenBank annotation vocabulary” is incomplete. The cited documentation also shows objective forms such as `~no`, `~keep`, `~change`, `~gc`, `~tm`, and `~CodonOptimize`; the parser accepts registered long class names as well as shorthands, and `AllowPrimer` supplies the supported `primer` shorthand. The list is evidence of no HGVS syntax, but it is not complete.

17. **Status: misquoted — record line 165.** The shown runtime exception text is truncated. The actual concatenated message is `Using avoid_heterodimerization requires primer3 installed (pip install primer3-py)` (`AvoidHeterodimerization.py:48-53`), not the complete `ImportError("Using avoid_heterodimerization requires primer3 installed")` expression printed in the record.

18. **Status: wrong-location — record lines 173 and 226.** `breaches_records_to_pdf` does not export the breaches DataFrame to Excel. The caller exports Excel via `dataframe.to_excel(...)`, then separately calls `records_from_breaches_dataframe(...)` and passes those records to `breaches_records_to_pdf(...)` (`constraints_reports.py:4-10`; example script lines 27-30). Line 173’s phrase “exporting to Excel with `breaches_records_to_pdf`” attributes the Excel operation to the wrong function; line 226 omits the required record-conversion step.

19. **Status: number-wrong — record lines 187 and 254.** “22 open issues” is the GitHub repository API’s `open_issues_count`, which includes pull requests. Recomputing the typed records gives 19 open issues plus 3 open pull requests = 22 aggregate open issue/PR objects. Reporting the aggregate as 22 issues is a category/count error.

20. **Status: uncited — record lines 7 and 187.** “Current maintainer Peter Vegh” and “Maintenance has passed from Zulkower to Peter Vegh” are not established by any cited source. `veghp` authored the latest `master` commit, and Peter Vegh is associated with EGF, but neither the repository nor PyPI metadata designates him maintainer or documents a handoff. This is an inference presented as a sourced fact.

21. **Status: dead-link — record lines 129 and 204.** The cited `sequence_without_tf_binding_sites.py` says it downloads data from `http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt`. The literal URL used by the script no longer returns the dataset; current RegulonDB download infrastructure is at different routes and may require access through its current portal. Consequently the present-tense “downloads” claim and the example’s data-fetch step are no longer live as cited.

22. **Status: wrong-location — record line 225.** `MutationSpace.plot(ax, ...)` exists in source, but it is not on the rendered public core-classes API page claimed by “all on the public API reference page.” The page lists the other named documented methods through `string_representation`; `plot` has no docstring and is absent from the rendered member list.

23. **Status: number-wrong — record lines 94 and 225.** `space_size` is not an exact cardinality in general. Its implementation returns `exp(min(100, sum(log(choice_counts))))`, capping any true product above `exp(100)`, and returns 0 when there are no multichoices even though the empty product has one sequence (`MutationSpace.py:96-104`). Calling it cardinality/product without this cap is numerically wrong for those cases.

24. **Status: minor-discrepancy — record line 235.** The GeneBlocks diff figure is conditional, not a guaranteed provenance artifact. It is generated only when `GENEBLOCKS_AVAILABLE and plot_figure`; `geneblocks` is in the `tests` extra, not the `reports` extra (`pyproject.toml:26-42`; `optimization_reports.py:359-364`). When unavailable, the PDF template prints an installation note instead (`optimization_report.pug:43-49`).

25. **Status: misquoted — record line 242.** The primer-example source does not say DNA Chisel is “not originally designed for sequence collection creation.” Its actual wording is: “DNA Chisel is not originally meant for creating collections of sequences ... it is still possible to create collections of inter-compatible sequences” (`examples/common_scenarios/primers_collection.py:1-10`). The record presents a paraphrase as a quotation.

26. **Status: uncited — record lines 254 and 266.** “Installable and runnable today” is stronger than the evidence collected. PyPI availability supports “installable,” but the record explicitly says no installation was performed and no runtime behavior was exercised. “Runnable today” therefore has no execution evidence in this audit trail.
