# DNA Chisel — repair change log

**Record repaired:** `revision/tool_survey/final/dnachisel.md` (edited in place)
**Audits processed:** `citation_audit/dnachisel.md` (6 "not found" + 26 "other"), `factcheck/dnachisel.md` (7 section A + 12 section B + 1 section C)
**Repair date:** 2026-08-14
**No capability value was changed.** All 24 rows hold their original yes / partial / no.

## How verification was done

Primary sources used for every finding below — the audits' assertions were never taken on trust:

| Source | How obtained |
|---|---|
| Repository at the exact audited commit `68c09304341c3656f3dfe63eda37757d6a7b3917` | tarball downloaded from GitHub and read directly (`dnachisel/`, `docs/`, `examples/`, `tests/`, `pyproject.toml`, `changes.md`, `README.rst`) |
| Paper | `papers/Zulkower2020_dnachisel.pdf`, full text via PyMuPDF |
| Supplementary Information (SI1) | local 18-page `poolparty-statecounter/literature/Zulkower2020_DNAChisel_all.pdf` **and** the 21-page Edinburgh author PDF at `https://www.pure.ed.ac.uk/ws/portalfiles/portal/154025952/DnaChisel_with_SI.pdf` (HTTP 200) |
| PyPI | live `https://pypi.org/pypi/dnachisel/json` and `.../dnachisel-dtailor-mode/json` |
| GitHub API | repo metadata, `commits/master`, issue-vs-PR search counts |
| Rendered docs | live `https://edinburgh-genome-foundry.github.io/DnaChisel/ref/core_classes.html` |
| Biopython class hierarchy | `Bio/SeqFeature.py` at `master` and at tag `biopython-179` |
| Live URLs | RegulonDB dataset URL, CUBA web app, docs site |

Nothing was installed. No file outside `revision/tool_survey/` was written.

---

## APPLIED

### From the citation audit — "NOT FOUND IN ANY SOURCE"

**CA-N1 — line 58, `variant` / `library` absent from the source tree.** *Applied.*
Verified: `grep -rli variant dnachisel/` hits 10 files, including the public `MutationChoice.variants`, `MutationSpace.all_variants` and `DnaNotationPattern.all_variants`. The audit is right about `variant`. It over-reaches on `library`: the three hits are prose ("the snapgene_reader library", "the genome_collector library", "used throughout the library"), not anything *named* library. Correction reflects both facts. Corroborated by fact-check A1.

**CA-N2 — line 86, "the only multi-sequence API".** *Applied.*
Verified in source: `AvoidMatches(match_length, bowtie_index=None, sequences=None, …)` builds a temporary Bowtie index from a supplied set; `AvoidBlastMatches(blast_db=None, sequences=None, …)`; `AvoidHeterodimerization(other_primers_sequences, …)`; `MotifPssmPattern.from_sequences`. Reworded to distinguish specifications that *consume* a sequence set as reference from the one API that takes a collection as its object of work.

**CA-N3 — lines 102 and 129, `MotifPssmPattern` avoidance-only.** *Applied.*
Verified three ways: `EnforcePatternOccurence.__init__` accepts any `SequencePattern` and `evaluate` calls `self.pattern.find_matches`; `MotifPssmPattern` implements `find_matches_in_string`, so it works there; `docs/ref/core_classes.rst:54-55` states "Specifications ``AvoidPattern``, ``EnforcePatternOccurence`` accept a Pattern as argument"; and SI1 Note 1 (read in the Edinburgh PDF) lists a PSSM among patterns for "specifications such as AvoidPattern and EnforcePattern". Corroborated by fact-check A2. The construction methods `from_sequences` / `from_file` (fact-check B9) were folded into the same clause.

**CA-N4 / FC-A7 — line 133, "every objective is a biochemical or manufacturing score".** *Applied.*
Verified: `EnforceChanges` and `AvoidChanges` score nucleotide edit counts, not biochemistry (`EnforceChanges.py`, `AvoidChanges.py`); the paper states DNA Chisel "allows to define any new sequence specification that the Python language can express". Narrowed to shipped objectives and added the extensibility caveat. **See escalation E4** — `assay_insilico = no` now rests solely on the absence of any shipped model-probing workflow.

**CA-N5 — line 183, "the OSI MIT classifier".** *Applied.*
Verified on live PyPI JSON: `info.classifiers` is `[]`. `info.license_expression = "MIT"` is real, as are `LICENSE` and the README sentence, so "confirmed three ways" still stands.

**CA-N6 — line 267, Supplementary Sections "unobtainable".** *Applied.*
Verified twice over: a local 18-page `Zulkower2020_DNAChisel_all.pdf` contains SI1 pages 1–16 with a contents page listing A–E, and the Edinburgh URL returns HTTP 200, `application/pdf`, 21 pages. Rewritten to say the sections were absent from `papers/` and that the record's SI-derived statements have not been re-checked against the actual supplement. **See escalation E3.**

### From the citation audit — "Other non-verified evidence"

**CA-2 — line 44, `dnachisel/` "contains only".** *Applied.* Directory listing at the cited commit also has `__init__.py`, `version.py`, `README.md`. Added.

**CA-3 — lines 46 and 106, `random_compatible_dna_sequence` signature.** *Applied.*
Verified in `dnachisel/utils/utils.py:6-14`: the signature ends `max_random_iters=5000, logger="bar", **kwargs`. Both places that present the signature with defaults now reproduce it in full. Line 224 was left alone — it lists bare parameter names with no defaults and does not read as a quotation.

**CA-4 — lines 73 and 86, dtailor-mode "released once".** *Applied.*
Verified on PyPI: releases 0.0.1, 0.0.2, 0.1.0, 0.1.1, all uploaded 2020-01-14. The record contradicted itself (line 73 already said "all four releases"). Both the evidence text and its self-quotation in the corrections table were fixed.

**CA-5 — line 86, per-CDS loop attributed to two scripts.** *Applied.*
Verified by reading both scripts. `genome-wide-optimization.py:11-32` does build one problem per CDS, append to `optimized_records`, and `SeqIO.write(...)` a multi-record FASTA. `ecoli_genes_optimization.py:24-40` loops over `gene` features only to accumulate constraints, builds one problem over the whole record, and writes one GenBank. Attribution corrected.

**CA-6 — line 98, `changes.md` quotation.** *Applied.* Source reads "AvoidChanges and EnforceChanges can now tunable"; the record had added "be". Now quoted exactly with `[sic]`.

**CA-7 — line 102, "will be preferentially placed".** *Applied.* `docs/genbank/genbank_api.rst:149-150` reads "will be be preferentially placed towards the center of the selected region". The quotation boundary was moved so what is inside the quote marks is now an exact substring of the source, without importing the typo.

**CA-8 — line 106, "Generate N short sequences".** *Applied.* `utils.py:40-46` returns the single value `problem.sequence`. Retitled to one sequence per call, with the user loop stated.

**CA-9 — line 111, "features for every specification applied".** *Applied.* `RecordRepresentationMixin.py` filters constraints with `if cst.__dict__.get("location", False)`; objectives are unfiltered. `EnforceRegionsCompatibility.__init__` sets `self.locations` and never `self.location`, so it is skipped. Corrected precisely.

**CA-10 / FC-A4 — line 111, "one feature per nucleotide edit".** *Applied.* `sequence_edits_as_features` calls `sequences_differences_segments`, which merges consecutive differing positions into `(start, end)` runs and emits one feature per run. Corrected to "per contiguous run of edited nucleotides".

**CA-11 / FC-A3 — lines 111, 155, 235, report-bundle conflation.** *Applied in three places.*
Verified by reading both functions. `write_optimization_report` (`optimization_reports.py:286-423`) writes `constraints_before_and_after.csv`, `objectives_before_and_after.csv`, `Report.pdf`, `final_sequence_with_edits.gb`, `final_sequence.gb`, and the source copy — **no** `logs.txt`, **no** `plots.pdf`. `write_no_solution_report` (`:82-190`) writes `plots.pdf`, `local_constraints_breaches.gb` and `logs.txt`. Line 111 now names the success bundle and says `logs.txt` is not in it; line 155 splits the outputs into success and failure bundles; line 235 drops `logs.txt`.

**CA-12 — lines 111 and 235, source hash/copy overstated.** *Applied.* Both functions emit the copy only `if file_path is not None`, and `file_hash = hashlib.md5(file_content).hexdigest()[:8]`. Both places now say "first 8 hex digits" and "only when the problem was built from a file".

**CA-13 — line 119, what `write_optimization_report` plots.** *Applied.*
Verified against `optimization_report.pug`: the PDF holds a success/fail table with before/after sequenticons, an optional GeneBlocks diff image (with an install note in its place otherwise), and the two before/after dataframes. No `SpecAnnotationsTranslator` whole-sequence plot is produced by this function; that path is in `write_no_solution_report`. Note the diff figure *does* pass `translator_class=SpecAnnotationsTranslator` to GeneBlocks, and the replacement text says so rather than denying it outright.

**CA-14 — line 119, `max_features_in_plots`.** *Applied.* `grep` finds the name only at `optimization_reports.py:294` (signature) and `:326` (docstring) — never in the body. Now stated as accepted but never read. Line 246's §8 bullet was left alone (see "left with noted tension").

**CA-16 — line 151, "the complete GenBank annotation vocabulary".** *Applied.*
Verified by extracting every `@`/`~` token from `docs/genbank/genbank_api.rst`: the file also documents `~no`, `~change`, `~CodonOptimize` and the long class-name forms `@AvoidPattern`, `@AvoidChanges`, `@EnforceTranslation`, `@EnforceSequence`, `@EnforceChoice`, `@EnforcePatternOccurence`. "Complete" dropped and the omitted forms added. The list's purpose — showing there is no HGVS syntax — is unaffected.

**CA-17 — line 165, `ImportError` text.** *Applied.* `AvoidHeterodimerization.py:48-53` concatenates to "Using avoid_heterodimerization requires primer3 installed (pip install primer3-py)". Completed.

**CA-18 — line 173, Excel attributed to `breaches_records_to_pdf`.** *Applied.*
Verified in the shipped example (`constraints_breaches_report.py:27-30`): `dataframe.to_excel(...)`, then `records_from_breaches_dataframe(...)`, then `breaches_records_to_pdf(...)`. Corrected. Line 226 was left alone — it lists the two artifacts with a "+" and does not attribute the Excel export to the PDF function.

**CA-19 — lines 187 and 254, "22 open issues".** *Applied.* GitHub API: `open_issues_count` 22; issue search gives 19 open issues and 3 open PRs. Both places now state the split.

**CA-20 — lines 7 and 187, "current maintainer Peter Vegh".** *Applied.*
Verified: `pyproject.toml` has `authors = [{ name = "Zulko" }]`; PyPI `info.author` is "Zulko" and `info.maintainer` is `null`; no CODEOWNERS. Against that, the GitHub user `veghp` is "Peter Vegh", company `@Edinburgh-Genome-Foundry`, and all 30 most recent `master` commits (back to 2024-12-17) are his. So the identification is solid but the *designation* is an inference. Both places now say so.

**CA-21 — line 129, RegulonDB download.** *Applied (two words).*
Verified: the script's literal URL is `http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt`; port 80 refuses connection, and the https equivalent returns a 1653-byte RegulonDB Browser SPA shell, not the dataset. "downloads" → "attempts to download". Line 204 (a table row naming the script) left alone.

**CA-22 — line 225, `MutationSpace.plot` on the public API page.** *Applied.*
Verified against the live rendered `ref/core_classes.html`: the `MutationSpace` member anchors are exactly `all_variants`, `apply_random_mutations`, `choices_span`, `constrain_sequence`, `from_optimization_problem`, `localized`, `pick_random_mutations`, `space_size`, `string_representation`. `plot` has no docstring in source and does not appear. Moved out of the "all on the public API reference page" scope. Lines 66 and 94, which list exactly the nine real members, were already correct and were not touched.

**CA-23 — lines 94 and 225, `space_size` as exact cardinality.** *Applied.*
Verified `MutationSpace.py:96-104`: returns `0` when `len(self.multichoices) == 0`, otherwise `np.exp(min(100, np.log(choices).sum()))` — a saturating float, not an exact product. Both places corrected.

**CA-24 — line 235, GeneBlocks figure conditional.** *Applied.* `if GENEBLOCKS_AVAILABLE and plot_figure` in source; `pyproject.toml` puts `geneblocks` under `tests`, not `reports`; the pug template prints an install note otherwise. Now stated as conditional at both line 119 and line 235.

**CA-25 — line 242, primer-example quotation.** *Applied.* `primers_collection.py:1-10` reads "DNA Chisel is not originally meant for creating collections of sequences (frameworks such as D-tailor were written with this purpose in mind), but it is still possible to create collections of inter-compatible sequences." The paraphrase-as-quotation was replaced with the real sentence.

**CA-26 — lines 254 and 266, "Installable and runnable today".** *Applied.* The record itself states no installation or execution was performed. "and runnable" removed, with a pointer to §10.

### From the fact-check — section B (incomplete)

Each added as a brief clause inside an existing entry; no new sections.

**B2 — `EnforceSequence` / `EnforceChoice` work-in-progress warning.** *Applied.* Both module headers in v3.2.16 literally read `"""Implement EnforceSequence (DO NOT USE YET: Work in progress, stabilizing)"""`. Added as one §8 bullet.

**B3 — `SequenceLengthBounds` is verification-only.** *Applied.* Docstring: "Quite an uncommon specification as it can't really be solved or optimized. But practical as part of a list of constraints to verify." Added as a parenthetical where the class is listed at line 173.

**B5 — `AvoidChanges` edit budgets.** *Applied.* `AvoidChanges(max_edits=0, max_edits_percent=None, …)` verified in source. Added to the primitives list at line 98, which already enumerates the mutagenesis primitives.

**B6 + B11 + FC-A5 — default parser registry gaps.** *Applied as one clause.*
Verified against `builtin_specifications/__init__.py:32-58`: `DEFAULT_SPECIFICATIONS_DICT` omits `AvoidMatches` (whose `shorthand_name` is `no_match`, not `avoid_matches`), `AvoidStopCodons`, `EnforceTerminalGCContent` and `SequenceLengthBounds`. `EnforceTerminalGCContent` is additionally absent from `docs/ref/builtin_specifications.rst`. Folded into the existing §8 GenBank bullet.

**B7 — `UniquifyAllKmers` hard-constraint warning.** *Applied.* Class docstring: "For sequences with subsequences appearing more than 2 times, the specification may not work as a problem constraint, but will work as a problem optimization objective." Added where the class is listed at line 173.

**B8 — `snapgene_reader` undeclared.** *Applied.* `genbank_operations.py:22-28` wraps the import in try/except and raises `ImportError` otherwise; the package appears nowhere in `pyproject.toml`. Added to the §8 dependency bullet.

**B9 — PSSM construction.** *Applied*, folded into the CA-N3 / FC-A2 correction at line 102 (`from_sequences` / `from_file`).

**B10 — `AvoidBlastMatches` local-only and deprecated in favour of `AvoidMatches`.** *Applied.* Source docstring opens "WARNING: try using AvoidMatches instead, it is much better!!" and "Only local BLAST is supported/tested as for now". Added to the §8 dependency bullet.

---

## REJECTED

**CA-Other-1 — line 18, "all quotes below verified verbatim" is false.** *Rejected.*
The sentence sits inside the **PDF row** of the sources table, so "all quotes below" naturally scopes to quotations from the paper. I checked every paper quotation in the record against the PyMuPDF extraction — the abstract sentence, the two Usage sentences, the Implementation two-step and stochastic/exhaustive sentences, the `insert(CGTCTC)` sentence, "comprehensive optimization reports for traceability and troubleshooting", the Fig. 1A NUM1 sentence, and "It can be used as Python library, a command-line interface, or a web application". **All are verbatim.** The audit reached its conclusion by reading the claim as covering repository and docs quotations too, which is a stretch. The genuine repository/docs misquotations it found (CA-6, CA-7, CA-25) are all fixed above on their own merits.

**FC-A5 — the record presents a stale `@avoid_matches` as current default vocabulary.** *Rejected as stated; underlying fact applied elsewhere.*
The record's sentence attributes the vocabulary to `docs/genbank/genbank_api.rst`, and the docs do document `avoid_matches` (line 283, with the caution "Only works if supported by the server"). So the record is not wrong about its cited source. The audit's *code* claim is correct and verified, and I recorded it in the §8 GenBank bullet (B6/B11 above) rather than in the `hgvs_input` entry, where it has no bearing.

**FC-B12 — a shipped manuscript example "collects the resulting sequences" from repeated boost profiles.** *Rejected.*
I read `examples/manuscript_examples/B_multiobjectives/competing_objectives_example.py` in full. It builds 11 boost profiles, solves each, and stores `results[name] = get_scores_from_problem(problem)` — **objective scores only**. The `problem` variable is overwritten each iteration and the optimized sequences are never retained; the script ends by writing `optimization_data.xlsx` and a bar chart SVG. The example's own README confirms it "explores different variations of a problem" and plots scores. The audit's premise — that this is a hand-built way to obtain a set of different designed sequences — is false, so no addition was made.

---

## ESCALATED

**E1 — CA-15: the `CompoundLocation` mechanism (record lines 70, 147, 244, and referenced at 260).**
The record says `Location.from_biopython_location` reads only `.start/.end/.strand`, "so a Biopython `CompoundLocation` **silently collapses to its outer span**", making a spliced CDS one contiguous intron-including region.

What I could verify: the description of `from_biopython_location` itself is exactly right. But `CompoundLocation` is **not** a subclass of `FeatureLocation` in any Biopython version (in ≥1.81 `SimpleLocation`/`FeatureLocation` and `CompoundLocation` are siblings under the `Location` ABC; in ≤1.80 `CompoundLocation` is a standalone class) — checked against `Bio/SeqFeature.py` at `master` and tag `biopython-179`. Every specification standardizes its location through `Location.from_data`, which tests `isinstance(location_data, FeatureLocation)` and, finding no match, falls off the end of the function and returns `None`. `list_from_biopython_feature` calls `from_biopython_location` only to format an error message; it hands the raw `feature.location` to the spec, which again goes through `from_data`. So on the normal paths the location is silently *dropped*, not collapsed.

Why I did not fix it: the audit is right that the record's mechanism is not the normal path, but I cannot verify the *replacement* mechanism without executing DnaChisel (Biopython is not installed here and I may not install). A dropped location makes `initialize_location_from_problem` true, which routes to `_copy_with_full_span_if_no_location` and then back through `set_location`'s "should be a 3x" check — whether that raises or silently spans the whole record depends on runtime state I cannot observe. The claim appears in three places and underpins the `exon_intron_split_codons` evidence, so guessing would mean writing three passages of unaudited text about behaviour I have not seen.

Note two verified facts that bear on the answer: the shipped `ecoli_genes_optimization.py:25` explicitly guards with `len(feature.location.parts) == 1`, i.e. it skips compound locations; and `from_record` parses only `misc_feature` labels, never CDS features.

**Question for the authors:** should the mechanism be restated as "the location is silently dropped by `Location.from_data`" (verifiable from source), or should someone install DnaChisel and run the spliced-CDS case to determine the actual end-to-end behaviour before rewording? Either way the value `no` looks safe — there is no exon-array handling on any reading — but the evidence sentence should not stay as written.

**E2 — FC-B1: reverse translation and degenerate-pattern expansion.**
Verified present and exported: `reverse_translate(protein_sequence, randomize_codons=False, table="Standard")`, `random_protein_sequence(length, seed=None)`, and `DnaNotationPattern.all_variants()`, which returns *all* concrete ATGC sequences of an IUPAC-degenerate pattern via `itertools.product`.

Not applied because the third one is not a neutral addition: a documented public method that expands a degenerate pattern into its full concrete variant set touches what `lazy_evaluation = partial` and `mixed_mutagenesis_one_pool = no` rest on, and the record's current `lazy_evaluation` argument is built entirely around `MutationSpace.all_variants`. **Question:** should `DnaNotationPattern.all_variants()` be added to the `lazy_evaluation` evidence, and does its existence change either rating?

**E3 — the Supplementary Information is now in hand and has not been mined.**
Following CA-N6 I read SI1. One item is directly relevant and no audit raised it: the SI's built-in-specification table gives, as the *use case* for `AvoidMatches`, "Creation of 'barcode' sequences which are mutually orthogonal, or orthogonal to given construsts or genomes (as is Casini et al, 2014)". The record's `barcode_generation` entry says the string "barcode" appears nowhere — but scopes that claim to the built-in specification list, the biotools reference, the docs TOC and the `builtin_specifications/` listing, none of which is the SI, so the sentence is not falsified and I left it alone.

**Question:** `barcode_generation` is already flagged in the record as the one definitional fork. Does the authors' own SI advertising `AvoidMatches` for orthogonal barcode-set creation change that fork's resolution? This is a rating decision, not an evidence fix. The rest of SI1 (sections B–E, the worked examples, and the algorithm description in S1E) has also not been checked against the record's inferences.

**E4 — `assay_insilico` evidence base after the CA-N4/FC-A7 fix.**
The universal claim "every objective is a biochemical or manufacturing score" was false and is corrected, but that claim was one of the row's two supports. The row now rests entirely on "nothing in the paper, docs, README, source tree or examples concerns probing predictive or generative genomic AI models", plus the new sentence noting that a user-written `Specification.evaluate` is unconstrained. **Question:** is the surviving support sufficient for `no`, given the framework's explicit extensibility?

**E5 — FC-B4: passive objectives and priority scheduling.**
Verified: `Specification.as_passive_objective()` exists with the documented semantics (contributes to the global score during other objectives' passes without getting its own pass), and specifications carry a `priority` attribute. Not applied because no existing entry in the record is a natural home — the rows are about library-scale design, not solver scheduling — and whether this omission is material is an emphasis judgment. **Question:** does this belong anywhere, e.g. as a candidate row in §7?

**E6 — FC section C: balance.**
The fact-check calls the record "broadly fair but somewhat limitation-heavy", with the single-sequence boundary repeated across many rows. That is exactly the emphasis judgment reserved for the authors, and acting on it would mean rewriting far more than a paragraph. No change made. Worth noting that this repair pass has, on net, moved the record slightly *toward* the tool — positive PSSM enforcement, additional multi-sequence-consuming specifications, `variant` naming, `AvoidChanges` edit budgets — while also adding limitations (undeclared `snapgene_reader` and `geneblocks`, the WIP headers, the parser-registry gaps, the inoperative `max_features_in_plots`).

---

## Left in place, with noted tension

These are places where a corrected fact makes a nearby unfixed sentence read slightly oddly. Per the surgical-editing rule they were not rewritten.

1. **Line 260** ("Highest confidence — re-verified verbatim from raw source this round") still lists "the report output filenames" among items re-verified. Several of those filenames were wrong and are now fixed, so that meta-claim about the *previous* round is inaccurate as to that item. The CSV columns it also lists **were** verified correct against `constraints_before_after_dataframe` / `objectives_before_after_dataframe`, as was the `random_compatible_dna_sequence` docstring, `sequence_to_biopython_record(sequence, id="<unknown id>", …)`, `load_record` → `SeqIO.read`, the full text of `example_EnforceRegionsCompatibility.py`, the `requires_dist` placement of `primer3-py`, the README sentences and the `docs/notes.rst` determinism contract.
2. **Line 217** still lists "#4/#11 per-gene domestication across a genome" as the demonstration where PoolParty would replace a user `for` loop. After the CA-5 fix, only #11 (`genome-wide-optimization.py`) is a per-CDS problem/output loop; #4 (`ecoli_genes_optimization.py`) produces a single record. The pairing is now looser than the corrected line 86.
3. **Line 226** (§7 item 4) still reads "exporting Excel + `breaches_records_to_pdf`", which compresses away the `records_from_breaches_dataframe` step now spelled out at line 173.
4. **Line 246** still carries the §8 bullet "Report plotting degrades with many features ('plots with thousands of features may take ages to plot')". The quotation is verbatim from the `write_optimization_report` docstring, but line 119 now records that the `max_features_in_plots` knob it describes is never read.
5. **Line 169** ("Provenance modernized — the extraction cited 2020 sources for 2025 software") now sits next to the corrected chronology showing that the `EnforceTranslation` additions predate the paper's submission rather than following publication. The broader point — the paper cites Nakamura 2000 while current DnaChisel depends on `python_codon_tables`, verified in `requires_dist` — is unaffected.
6. **Line 261** still refers to "the complete GenBank annotation vocabulary" in the confidence summary, while line 151 no longer claims completeness.

---

## Pass 2 — policy application

**Date:** 2026-08-14

**Outcome count (six open items): 3 applied · 2 declined-by-Policy-B · 1 declined-by-Policy-A · 0 rejected-unverifiable.** No capability value was changed. E4 also produces one locked-row value escalation; two additional locked-row reconciliation escalations exposed by the corrected evidence are reported below.

**Verification basis:** the main paper and the University of Edinburgh author PDF containing SI1 were read with the required PyMuPDF command. Official source was checked at tag `v3.2.16` / commit `68c09304341c3656f3dfe63eda37757d6a7b3917` from a source archive, without installation. Current official docs, PyPI, GitHub repository/commit/issue metadata, the CUBA route and the EGF software-suite page were fetched on 2026-08-14. Biopython's official API documentation was used only to verify that `FeatureLocation` is the simple-location alias while `CompoundLocation` is a sibling location class.

| Open item | Outcome | Edit and verification |
|---|---|---|
| **E1 / CA-15 — `CompoundLocation` mechanism** | **applied** | Corrected the mechanism in the corrections table, `exon_intron_split_codons`, §8 and §10. In v3.2.16, `EnforceTranslation.__init__` calls `Location.from_data`; that function handles tuple/list, Biopython `FeatureLocation` and DNA Chisel `Location`, but has no `CompoundLocation` branch and falls through to `None`. `EnforceTranslation` consequently marks the location for initialization, and `Specification._copy_with_full_span_if_no_location` supplies `Location(0, len(problem.sequence))`. `from_record` independently parses only labelled `misc_feature` entries, not CDS features. Official Biopython docs confirm `FeatureLocation` is the `SimpleLocation` alias and `CompoundLocation` is its sibling under the location base. The original outer-span-collapse claim was false on the normal path. The corrected evidence still supports `exon_intron_split_codons = no`; value unchanged. |
| **E2 / FC-B1 — reverse translation and degenerate expansion** | **applied where score-relevant; `random_protein_sequence` declined-by-policy within the item** | Added brief clauses to existing entries only. `reverse_translate(protein_sequence, randomize_codons=False, table="Standard")` uses the first valid codon or a random synonymous codon and emits one DNA sequence; this bears on locked row 7, although it is not a substitution-library generator. `DnaNotationPattern.all_variants()` uses `itertools.product` and returns an eager list of every concrete ATGC sequence represented by one IUPAC pattern; this bears on locked rows 1 and 5. `random_protein_sequence` merely emits one toy M…`*` protein and changes no locked-row score or Table 1 cell, so Policy B excludes it. No capability value changed. |
| **E3 — SI1 `AvoidMatches` barcode use case** | **declined-by-policy** | No capability prose was added. SI1 p. 2 was verified directly: its use case for `AvoidMatches` is creation of mutually orthogonal "barcode" sequences or sequences orthogonal to constructs/genomes. That supports the already-recorded `barcode_generation = partial` fork but does not change any of the 17 locked rows, none of which is a barcode row, and does not alter Table 1 Purpose, Key features, Output or Availability. The SI availability hedge in §10 was refreshed because the supplement was fetched, but the omitted use case was not inserted. |
| **E4 — `assay_insilico` value basis / custom objectives** | **evidence applied; locked-row score escalated** | Replaced the weak extensibility caveat with the verified behavior: `Specification.__init__(evaluate=...)` accepts an evaluation callable, user subclasses may implement `evaluate`, and `ObjectivesMaximizerMixin` optimizes the resulting total score by exhaustive or random mutation. Thus a predictive model can be wrapped and optimized, although no shipped specification, example or documentation does so. The legacy `assay_insilico = no` remains grounded in the absence of a shipped model-probing assay workflow. Locked row 12 has a different binding test and is escalated below; no value was changed. |
| **E5 / FC-B4 — passive objectives and priority scheduling** | **declined-by-policy** | Verified in `Specification.py`: `priority = 0` is the class default and `as_passive_objective()` returns a copy with `optimize_passively=True`; `ObjectivesMaximizerMixin.optimize()` omits passive objectives from their own passes while their scores remain in the global objective sum. These are solver-scheduling controls on one sequence, not two library-design operations sharing a carrier, so they do not change locked row 2. They also do not alter a Table 1 cell. No edit. |
| **E6 / fact-check §C — balance and emphasis** | **declined-by-policy** | Policy A declines balance and proportional-emphasis edits. No factual finding was inferred from the balance judgment and no record text changed for it. |

### Policy C — primary-source boundary and current facts

**Outcome: applied.** Removed the prior-analysis source row, the related-project package/repository discussion, the DMS-use and Benchling adoption anecdotes, and all corresponding confidence/availability residue. The official-repository primer quotation was shortened so it no longer relies on the sibling-project parenthetical. Capability conclusions now use the paper, official v3.2.16 repository/docs, DNA Chisel's PyPI page and first-party project metadata.

Stale live-status language was also replaced rather than hedged. On 2026-08-14 PyPI still reported v3.2.16; GitHub still reported 278 stars, 53 forks, 19 open issues plus 3 open PRs, `archived=false`, `pushed_at=2025-11-05T10:28:23Z`, and `master` at `68c0930` from 2025-05-10. The advertised CUBA URL returned only the 1,035-byte Vue shell, while EGF's current official suite page says the website is "not public". The record now distinguishes a documented web interface from verified current public usability. `maintained` is consistently defined as maintenance through May 2025, not active ongoing development in August 2026; `interface = yes` remains supported independently by the Python API and CLI.

### Policy D — version drift

**Outcome: no edit.** The assessed source and package release are v3.2.16, and live PyPI on 2026-08-14 still lists v3.2.16 as current. There is no materially different current release requiring a parenthetical.

### Value and locked-row escalations — report only

No value was changed in either the record or the locked table.

1. **Locked row 12, `Model-guided variants`: current DNA Chisel cell `○`; corrected evidence supports `●` under the binding definition.** The definition gives `●` when an optimization loop acts on a model output and gives `◐` for an arbitrary scoring callable without optimization. DNA Chisel accepts an arbitrary `evaluate` function returning a score and its normal exhaustive/random objective loop optimizes that score. Requiring the user to wrap a predictor in `Specification.evaluate` is extension work, but the row's text does not exclude extension work. This is the E4 value escalation.
2. **Locked row 1, `Library object`: current DNA Chisel cell `●`; corrected evidence supports at most `◐`.** No DNA Chisel-defined library type exists. `DnaNotationPattern.all_variants()` returns a general Python list and `MutationSpace.all_variants()` returns a generator of bare strings. Under the row's explicit rule, a general-purpose container can earn only `partial`; the full `yes` cell conflicts with the final record's `library_as_object = no` evidence. Escalated without changing the cell.
3. **Locked row 5, `Pairwise and higher-order variants`: current DNA Chisel cell `◐`; the binding test may require `●`.** Both `MutationSpace.all_variants()` and `DnaNotationPattern.all_variants()` enumerate Cartesian products across independently variable positions, so outputs can contain two or more co-occurring events and the tool enumerates the combinations. The current `partial` seems to import a managed-library/provenance requirement that row 5 does not state. This requires the table authors' ruling; the cell was not changed.

**Pending row-7 scoring note:** `reverse_translate` is relevant to `Codon / amino-acid substitutions` but emits one back-translation of a user-supplied protein, not a substitution set. It supports at most `partial` under the locked definition; the current cell is `?`, so this is not a value change.

### Resulting neighbouring tension

The `assay_insilico = no` heading now sits beside an explicit statement that a custom predictive-model score can be optimized. That is intentional only under the legacy row's shipped-workflow reading; locked row 12 expressly scores extensibility and therefore conflicts, as escalated above. The qualified `interface` heading similarly retains `yes` because API + CLI are verified while no longer presenting the reachable CUBA shell as proof of a usable public web optimizer.

### Row-substitution candidate — report only, no action

**Provisional candidate:** replace the already-flagged near-uniform **`Codon / amino-acid substitutions`** row with **`Annotation-driven design specification`**.

- **Proposed test:** executable constraints and/or objectives can be encoded in sequence-feature annotations and parsed directly into a design problem without user-written orchestration. `partial` = annotation export only, or declarative input covering only a narrow fixed subset.
- **Ground:** discrimination and honesty. The incumbent is already flagged in the locked file as likely near-uniform. The candidate is a first-class DNA Chisel strength on which PoolParty is not best.
- **Evidence:** the paper describes `@` constraints and `~` objectives in GenBank; official docs expose the annotation grammar and CLI; `RecordRepresentationMixin.from_record` parses labelled `misc_feature` entries into constraints/objectives; and the stock registry covers multiple pattern, sequence, codon and physicochemical specifications. A tracked-source search of the current PoolParty repository found no `GenBank`, `SeqIO`, `misc_feature` or `from_record` parser path. Other tool cells still require primary-source scoring, so no substitution was made.

### Pass-2 totals

| Outcome | Count |
|---|---:|
| Applied open items | 3 (`E1`, score-relevant parts of `E2`, `E4` evidence) |
| Declined by Policy B | 2 (`E3`, `E5`; plus the `random_protein_sequence` subfinding inside E2) |
| Declined by Policy A | 1 (`E6`) |
| Rejected as unverifiable | 0 |
| Open-item score escalations | 1 (`E4` → locked row 12) |
| Additional locked-row reconciliation escalations | 2 (rows 1 and 5) |
| Capability values changed | **0** |
