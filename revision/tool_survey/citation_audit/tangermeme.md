# Citation-integrity audit — Tangermeme

Record audited: `25_poolparty_private/revision/tool_survey/final/tangermeme.md`

Repository baseline: `jmschrei/tangermeme` main commit `8d732b8c08764057b7ae5faa644d48664f36b44b` (2026-07-19). The complete repository tree, live GitHub/PyPI/Read the Docs endpoints, and all 13 pages of `papers/Schreiber2025_tangermeme_all.pdf` extracted with PyMuPDF were checked. This audit concerns citation integrity only.

## NOT FOUND IN ANY SOURCE

### 1. “Nucleotide-level only” is not established by the cited grep and is contradicted by the implementation

- **Record claim (lines 363–365 and 788–792):** `saturation_mutagenesis` is “Nucleotide-level only,” with zero occurrences of `codon`, `amino`, `protein`, `orf`, or `reading frame` offered as the evidence.
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** No source states that `saturation_mutagenesis` is restricted to nucleotides. Its implementation uses the input alphabet dimension dynamically (`X.shape[1]`) and obtains candidate symbols with `torch.where(~identity)`; it does not hard-code four nucleotide channels (`saturation_mutagenesis.py:203–232`). The README separately states that all functions are extensible to any alphabet. The grep verifies the absence of codon/translation-specific machinery, but it cannot support the stronger “nucleotide-level only” assertion.

### 2. “No CLI” is contradicted by the cited packaging source

- **Record claim (lines 584 and 794–795):** render the interface as “Python API only (no CLI, no GUI, no web service).”
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** `pyproject.toml:71–72` defines the console script `tangermeme-install-skills`, and `tangermeme/_skills/install.py:108–149` implements its `argparse` command-line interface. The record itself acknowledges this entry point at lines 590–592 and 739–740. A narrower statement such as “no analysis CLI” would match the evidence; the literal “no CLI” does not.

### 3. The claimed exhaustive list of sequences returned by the package is false

- **Record claim (lines 808–812):** “the only sequences it ever hands back are DeepLIFT backgrounds, generated controls, and design winners.”
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** The ordinary public `ersatz` operations `insert`, `substitute`, `multisubstitute`, `delete`, and `reverse_complement` all directly return sequence tensors, in addition to `randomize`, `shuffle`, and `dinucleotide_shuffle`. For example, `ersatz.insert` explicitly documents and returns an edited sequence tensor at `ersatz.py:20–94`, and `substitute` does so at `ersatz.py:97–195`. These outputs are omitted from the purported exhaustive list even though the record cites them elsewhere.

No wholly fabricated external-source quotation was found. The three items above are placed here because the claimed evidence is absent and the repository directly supplies counterexamples.

## Other non-verified evidence

### 4. Three abbreviated local source paths do not resolve

- **Record citations (lines 23–25):** `/mnt/c/.../prior_analyses/Schreiber2025_tangermeme_all_analysis.md`, `/mnt/c/.../tool_survey/extractions/tangermeme.md`, and `/mnt/c/.../tool_survey/reviews/tangermeme.md`.
- **Status:** **wrong-location**
- **Finding:** None of those literal paths exists. The corresponding files do exist under the full workspace path `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/`, but the ellipsized references are not resolvable file citations.

### 5. The raw GitHub directory URL is dead

- **Record citation (line 27):** `raw.githubusercontent.com/jmschrei/tangermeme/main/`.
- **Status:** **dead-link**
- **Finding:** The HTTPS form of this directory-prefix URL returns HTTP 404. Raw URLs for individual files below that prefix are live, but the cited prefix itself is not a live page or directory listing.

### 6. “Six ‘no’ cells” is not reproducible

- **Record claim (line 33):** the recursive-grep citation correction “applies to six ‘no’ cells.”
- **Status:** **number-wrong**
- **Finding:** The original extraction explicitly invokes the narrow grep in seven `no` cells: `barcode_generation`, `assay_dms`, `transcript_models`, `exon_intron_split_codons`, `hgvs_input`, `primer_design`, and `codon_optimization`. In the final record, the recursive absence search is used for at least nine `no` sections after also supporting `consequence_annotation` and `synthesis_constraints`. The count six cannot be reproduced from either document.

### 7. The five NamedTuples are not all simply “of model outputs”

- **Record claims (lines 53–55 and 781–783):** five of the six classes are “NamedTuples of model outputs.”
- **Status:** **minor-discrepancy**
- **Finding:** The class count and class locations are correct, but `AttributionReferencesResult` contains `references`, explicitly documented as the reference/background **sequence tensor** used for attribution (`results.py:41–52`). The record correctly notes this exception at lines 115–117, making the repeated unqualified characterization internally inconsistent.

### 8. Not every sequence collection has the stated three-dimensional shape

- **Record claim (lines 86–88):** “every sequence collection is a bare `torch.Tensor` of shape `(-1, len(alphabet), length)`.”
- **Status:** **minor-discrepancy**
- **Finding:** The bare-tensor portion is supported, but `randomize`, `shuffle`, and `dinucleotide_shuffle` return four-dimensional tensors with an added replicate axis, documented as `(-1, n, len(alphabet), length)`. The record itself gives that shape at lines 387–394 and 702–706.

### 9. Streaming is not exclusive to the Cartesian-product implementation

- **Record claim (lines 164–176):** “Streaming exists only on the Cartesian-product axis.”
- **Status:** **minor-discrepancy**
- **Finding:** `match.py` contains multiple genuine streaming generators, including `_sequence_generator_stream` (“only one sequence is kept in memory at a time”), `_count_generator_stream`, coordinate generators, and `numpy.fromiter` consumers (`match.py:22–145, 194–245`). The narrower claim that the package has no streaming *library abstraction* is supportable; the literal exclusivity claim is not.

### 10. Spacing is not the package's only enumerated axis

- **Record claim (line 224):** the spacing grid “is the only genuinely enumerated axis in the package.”
- **Status:** **minor-discrepancy**
- **Finding:** `randomize`, `shuffle`, and `dinucleotide_shuffle` expose an enumerated `n` replicate axis, while `saturation_mutagenesis` internally enumerates alphabet symbols by positions and returns outputs on that grid. The cited `space` source verifies the spacing axis, not its claimed uniqueness.

### 11. The greedy-substitution code block is edited pseudocode, not source code

- **Record evidence (lines 244–254):** a fenced Python block attributed to `design/greedy_substitution.py:188–205`.
- **Status:** **misquoted**
- **Finding:** The block is not verbatim. The source has `tqdm(motifs, disable=not verbose)`, not `tqdm(motifs, ...)`; reverse-complement augmentation is at lines 161–162, outside the cited range; `predict` arguments are replaced with `...`; and the three explanatory inline comments (“motifs already…,” “motif at EVERY…,” and “only the winner escapes”) do not exist in the file. The underlying enumeration behavior is real, but the displayed code is an annotated paraphrase presented as code from that location.

### 12. The DataFrame inventory is wrong, and `annotate_seqlets` does not return a DataFrame

- **Record claim (lines 282–291):** “The only DataFrames the package produces are `seqlet.recursive_seqlets` ... and `annotate.annotate_seqlets`.”
- **Status:** **wrong-location**
- **Finding:** `annotate_seqlets` actually returns two `torch.Tensor` objects (`annotate.py:68–89`); its return annotation at line 25 is stale and contradicts both the return documentation and implementation. Conversely, several omitted public functions do return DataFrames: `seqlet.tfmodisco_seqlets` (`seqlet.py:211–353`), `utils.example_to_fasta_coords` (`utils.py:233–335`), `match.extract_matching_loci` (`match.py:372–621`), and `io.read_vcf` (`io.py:622–656`), in addition to `recursive_seqlets`.

### 13. `one_hot_to_fasta` is not literally the package's only writer

- **Record claims (lines 297–300 and 504–507):** `one_hot_to_fasta` is the package's “only writer,” and all file-writing code is in `io.py`.
- **Status:** **minor-discrepancy**
- **Finding:** The narrow grep misses the installer: `tangermeme/_skills/install.py:95–104` removes an existing destination, creates directories, and copies the bundled skill tree with `shutil.copytree`. The evidence does support the narrower conclusion that `one_hot_to_fasta` is the only sequence/data serializer and that there is no VCF writer, but not the literal package-wide writer inventory.

### 14. `interactive_logo` is present in a default installation; its optional dependency is not

- **Record claims (lines 334–336 and 761–764):** `plot.interactive_logo` “is absent from a default install.”
- **Status:** **minor-discrepancy**
- **Finding:** The function ships in `tangermeme/plot.py:578` and remains importable after a default installation. What is absent is the optional `mpld3` dependency; the function imports it lazily and raises an error when invoked without the `[interactive]` extra. The dependency claim is verified, but the function-absence wording is not.

### 15. The “standard MPRA scramble control” assertion has no citation

- **Record claim (lines 391–394):** “Dinucleotide shuffles are *the* standard MPRA scramble control.”
- **Status:** **uncited**
- **Finding:** The cited `ersatz.py` source establishes how Tangermeme generates dinucleotide shuffles and their return shape. It does not establish this external claim about standard MPRA experimental practice, and no paper or documentation citation is supplied for it.

### 16. The `extract_loci` quotation is about 45 lines beyond the cited range

- **Record evidence (lines 440–459):** the `loci` quotation is cited to `io.py:246–265`.
- **Status:** **wrong-location**
- **Finding:** Lines 246–265 contain the function signature. The quoted text beginning “Either the path to a bed file...” occurs at `io.py:310–315`, not within a few lines of the cited range. The quotation itself is verbatim at the later location.

### 17. The reported eight “open issues” include three pull requests

- **Record claim (line 628):** “Stars / forks / open issues — 308 / 32 / 8.”
- **Status:** **number-wrong**
- **Finding:** GitHub's `open_issues_count` field is 8 but includes pull requests. Enumerating the open endpoints gives **5 open issues and 3 open pull requests**. The stars and forks values, and the raw API field value 8, were verified; labeling all eight as issues is inaccurate.

### 18. The k-mer return type and gapped-k-mer weighting are described incorrectly

- **Record claim (lines 720–724):** `kmers()` returns a SciPy sparse matrix, and `gapped_kmers()` “does the same” summed per-position weighting.
- **Status:** **misquoted**
- **Finding:** Despite a stale SciPy return annotation, `kmers.py:68–70` constructs and returns a dense `torch.Tensor`. `gapped_kmers` does return a SciPy CSR matrix, but its implementation accumulates `new_gkmer_attr / k` (`kmers.py:127–136`), i.e. the average across selected characters, not the claimed sum. Its docstring says “sum,” so the source documentation and implementation are themselves inconsistent; the cited lines 21 and 40–64 do not verify the actual return object or the gapped implementation.

### 19. `N` does not specifically mark positions “where motifs cancel”

- **Record claim (lines 736–738):** `greedy_marginalize` returns a construct “with `N` where motifs cancel.”
- **Status:** **minor-discrepancy**
- **Finding:** The bundled skill documentation does make that claim (`design.md:95–97`), but the implementation initializes the entire background as `N`, writes each selected motif into its chosen slice in order, strips terminal `N`s, and returns the result (`greedy_marginalize.py:242–246`). Thus internal `N`s ordinarily denote uncovered gaps; overlapping motif characters are overwritten by later assignments, not converted to `N` by a cancellation rule.

## Link check

The PyPI JSON endpoint, GitHub repository/API, Read the Docs root, and GitHub URL quoted from the paper were live when checked on 2026-08-14. The only dead URL-like citation is the raw GitHub directory prefix in item 5.

## Verified counts and anchors relevant to the findings

- The local paper is 13 pages, and the record's quoted paper passages on pp. 1, 3, and 9 were found.
- The repository contains 28 Python files: 20 directly under `tangermeme/`, six under `tangermeme/design/`, and two under `tangermeme/_skills/`.
- The bundled skill contains 15 Markdown files (`SKILL.md` plus 14 references).
- The class inventory is exactly six classes: five NamedTuples and `TangermemeWarning`.
- The stated recursive negative grep over those 43 files is reproducible, including the sole `adapter` hit in `pisa.py` and sole `edit distance` hit in `saturation_mutagenesis.py`.
- The docs tree contains 14 tutorial notebooks, three vignettes, one how-to notebook, and two paper notebooks, matching the inventory in §4.
- PyPI version 1.4.0, its 2026-06-25 upload time, the eight listed releases in the preceding 12 months, MIT metadata, Python `>=3.10`, repository dates, archive state, stars, and forks were verified.
