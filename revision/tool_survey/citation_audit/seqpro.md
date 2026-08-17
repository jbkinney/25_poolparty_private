# Citation-integrity audit — SeqPro

Record audited: `25_poolparty_private/revision/tool_survey/final/seqpro.md`

Repository baseline: `ML4GLand/SeqPro` main commit `63a843985d96dd3f5a7bc8cc20e8bd03f1dabdd9` (2026-07-27, version 0.22.0). The complete repository tree, the published 0.22.0 wheel, live GitHub/PyPI/documentation endpoints, GitHub commit and release history, and the local `Klie2023kg` bibliography entry were checked. There is no SeqPro paper PDF in the supplied `papers/` directory, consistent with the record's software-only citation. This audit concerns citation integrity only.

## NOT FOUND IN ANY SOURCE

### 1. Absence of `experimental/__init__.py` does not make the namespace unimportable

- **Record claims (lines 163–168, 194–197, and 434–437):** because `python/seqpro/experimental/` has no `__init__.py`, “`seqpro.experimental` is therefore not an importable package” and “cannot be imported at all.”
- **Status:** **NOT FOUND IN ANY SOURCE**
- **Finding:** No repository, packaging, or Python source supports that inference. Python supports implicit namespace packages without an `__init__.py`. The published 0.22.0 wheel contains both `seqpro/experimental/_experimental.py` and `seqpro/experimental/_visualizers.py`; importing `seqpro.experimental` yields a namespace module, and importing `seqpro.experimental._experimental` succeeds (`edit_distance("AC", "AT") == 1`). `_visualizers.py` separately fails at import time because it imports the removed names `gc_content_seqs` and `nucleotide_content_seqs`, exactly as its header says. The raw `experimental/__init__.py` 404 is real, but it cannot support the stronger non-importability claim.

## Other non-verified evidence

### 2. The recursive Git tree has 190 entries, not 190 blobs

- **Record claims (lines 22 and 142–147):** the untruncated recursive tree contains “190 blobs.”
- **Status:** **number-wrong**
- **Finding:** Recomputing from the cited GitHub recursive-tree API gives `truncated: false` and **190 total tree entries: 162 blobs and 28 trees**. The separate count of 27 Python files under `python/seqpro/` is correct.

### 3. `CHANGELOG.md` contains 51 releases, not 55

- **Record claims (lines 31, 72–77, 266–268, and 492–496):** `CHANGELOG.md` is “608 lines / 55 releases,” and greps are described as spanning “all 55 releases” of the changelog.
- **Status:** **number-wrong**
- **Finding:** The file is exactly 608 lines, but it contains **51** version headings matching `^## `. The separate PyPI count is **55 release keys**, so the record has conflated the PyPI release count with the number represented in the changelog. The whole-file greps still cover all 608 changelog lines.

### 4. `Ragged` supports two nested ragged axes

- **Record claims (lines 97–99 and 387–391):** the Rust-native `Ragged` container has “exactly one ragged axis.”
- **Status:** **number-wrong**
- **Finding:** Current 0.22.0 code supports **R=1 and R=2**. `rag/_core.py:139–140` rejects only `shape.count(None) >= 3`; `Ragged.from_lengths` accepts a tuple and constructs a shape containing two `None` axes at `rag/_core.py:190–209`; and `from_fields` explicitly accepts R=2 fields at `rag/_core.py:212–245`. `CHANGELOG.md:95–113`, the Ragged roadmap, and numerous R=2 tests independently document the two-level implementation and the `R <= 2` cap. The bundled `SKILL.md` sentence saying “exactly one” is stale and contradicted by the implementation and newer project documentation.

### 5. The built-in `sp.AA` translator does not translate RNA codons

- **Record claim (lines 201–203):** “`sp.AA.translate()` translates DNA/RNA→AA.”
- **Status:** **misquoted**
- **Finding:** The built-in `AA` object is created from a DNA codon table containing `T`, not `U` (`alphabets/__init__.py:4–69`). Its own translation documentation says the accepted canonical nucleotide bytes are `{A,C,G,T,a,c,g,t}` and treats all others as unknown (`alphabets/_alphabets.py:496–501`). Thus an RNA codon such as `AUG` is rejected with `validate=True` or mapped through the unknown policy otherwise; it is not translated as methionine by `sp.AA`. A user can construct a custom `AminoAlphabet` with RNA codons, but that does not make the cited built-in `sp.AA` claim true. The README/docs quick-start comment “DNA/RNA → amino acids” is itself inconsistent with the implementation.

### 6. Only version 0.22.0 carries the cited breaking-change entry

- **Record claim (lines 448–449):** “0.22.0 and 0.21.x carry **BREAKING CHANGE** entries for `Ragged.__getitem__` semantics.”
- **Status:** **number-wrong**
- **Finding:** In `CHANGELOG.md`, the only `### BREAKING CHANGE` heading is under **0.22.0** at lines 1–9. The 0.21.0, 0.21.1, and 0.21.2 sections contain no breaking-change entry. The corresponding GitHub release bodies likewise do not add such an entry for any 0.21.x release.

### 7. Only David Laub, not both named authors, committed in July 2026

- **Record claim (lines 82–85):** “Two of the paper's authors (Laub, Klie, UCSD) were committing to this repo in July 2026.”
- **Status:** **number-wrong**
- **Finding:** The GitHub commits API for 2026-07-01 through 2026-07-31 shows human commits by **David Laub / d-laub** and automated version bumps by `github-actions[bot]`; it shows no Adam Klie commit. The API's latest commits attributed to Adam Klie are dated **2023-10-31**. The asserted contributor count is therefore one, not two, and the record supplies no citation for the claim.

### 8. The Ragged guide inventory includes examples that are not in that guide

- **Record claim (lines 344–346):** `docs/ragged.md` contains a worked coverage example, ufunc arithmetic, `to_padded`, `to_packed`, and record layouts via `rag.zip`.
- **Status:** **misquoted**
- **Finding:** The coverage construction and ufunc arithmetic examples are present. Neither `to_padded` nor `to_packed` occurs anywhere in `docs/ragged.md` or on the cited rendered page. The record-layout example at `docs/ragged.md:133–146` uses **`ak.zip`**, not `rag.zip`; the file even says “`ak.zip` returns a Ragged automatically.” The methods and `seqpro.rag.zip` do exist elsewhere in the package, but not as the cited guide examples.

### 9. The source-wide `motif` grep has one hit, not zero

- **Record claim (lines 149–153):** a source-wide grep for `motif` “across every fetched module returns 0 hits.”
- **Status:** **number-wrong**
- **Finding:** A case-insensitive grep over `python/seqpro/` returns one hit, `python/seqpro/experimental/_experimental.py:134`: `# Performing motif anal`. It is only an unfinished comment and supplies no motif functionality, so the substantive capability evidence is unaffected. The later, narrower wording “0 substantive hits” is supportable; the literal zero-hit claim is not.

### 10. The six changelog `codon` hits are not all performance work

- **Record claims (lines 74–77 and 293–298):** all six `codon` occurrences in `CHANGELOG.md` are “translate-kernel performance work.”
- **Status:** **minor-discrepancy**
- **Finding:** The count of six is correct, but the occurrences include `test full codon table` at line 513 and `single source of truth for codon LUT index` at line 280, in addition to kernel-routing, compaction, and benchmark/speedup entries. The cited `179×` performance entry exists verbatim at line 285; it does not characterize all six hits.

### 11. The coordinate schema does not encode open-versus-closed intervals

- **Record claims (lines 232–245 and 411–414):** five schemas cover “0- vs 1-based and open vs closed interval conventions,” described as a correctness feature.
- **Status:** **minor-discrepancy**
- **Finding:** The five schema names are real, and `CoordSchema` stores a `zero_based: bool`. It has no field for interval closure. `set_schema` renames columns and attaches `coordinate_system_zero_based` metadata for Polars frames; it neither models an open/closed flag nor adjusts start/end values when moving between a 0-based schema and a 1-based schema (`_coords.py:10–27, 100–143`). The standard BED/GTF format names imply conventional closure to a knowledgeable reader, but the cited implementation does not explicitly handle that dimension as claimed.

### 12. The opening “verbatim” quotations are from the README, not both named sources

- **Record evidence (lines 45–53):** two passages are introduced as “From `README.md` / `docs/index.md`, verbatim.”
- **Status:** **wrong-location**
- **Finding:** Both quoted passages exist in `README.md:10–12` (with Markdown links rendered away), so the quotation text is real. They do not occur verbatim in `docs/index.md`: that page instead says “It makes almost zero compromises on speed...” and begins the input statement “All functions accept...,” with additional XArray wording. The combined attribution is overbroad; only the README is the verbatim location.

### 13. The rendered API pages are not “signatures only”

- **Record claim (lines 347–348):** the six API reference pages are “mkdocstrings-rendered signatures only.”
- **Status:** **minor-discrepancy**
- **Finding:** The live rendered pages include docstrings, parameter and return descriptions, and source-code disclosures in addition to signatures. The documentation configuration explicitly sets `show_source = true` (`zensical.toml:67`). The count of six files and the count of 14 top-level members in `docs/api/index.md` are correct.

### 14. The README does not end at `# More to come!`

- **Record claim (line 433):** “README ends with a bare `# More to come!`.”
- **Status:** **minor-discrepancy**
- **Finding:** The heading exists at `README.md:81`, but the file continues at line 83 with a full contributions/community paragraph and code-of-conduct link. It is the last heading, not the end of the README and not a bare terminal line.

### 15. The two experimental visualizers are not both histograms

- **Record claim (lines 194–196):** the two functions in `_visualizers.py` are “QC histograms.”
- **Status:** **minor-discrepancy**
- **Finding:** `plot_gc_content` calls `ax.hist` (`_visualizers.py:16`), whereas `plot_nucleotide_content` calls `ax.plot` (`_visualizers.py:26`) and is a positional line plot. The count and names of the two functions, their broken imports, and their non-public status are otherwise supported.

### 16. The README does not explicitly assign the sibling-package boundaries claimed

- **Record claim (lines 439–442):** the README suite sentence means reference acquisition, motif work, and model probing “are explicitly out of scope” for SeqPro.
- **Status:** **misquoted**
- **Finding:** The quoted README text says only that SeqPro is “heavily utilized” by SeqData, MotifData, SeqExplainer, and EUGENe. It does not state an explicit scope exclusion or assign those three responsibilities to those packages. The sibling repositories and SeqPro's own absence of those functions can support the division as an inference, but the cited README quotation does not make it explicit.

### 17. The “canonical” or “standard” shuffle-null assertions are uncited

- **Record claims (lines 218–221 and 384–386):** k-let/dinucleotide shuffling is “the canonical null/background for in-silico experiments” and “the standard dinucleotide-shuffle null for regulatory-genomics ML.”
- **Status:** **uncited**
- **Finding:** The repository verifies that `k_shuffle` preserves k-lets, dispatches to the Rust implementation, and is benchmarked. No SeqPro source, documentation, paper, or other cited authority establishes the external field-practice claims “canonical” or “standard.”

### 18. “edgeR-style” TMM normalization is uncited

- **Record claim (line 371):** `seqpro.transforms.TMM` performs “edgeR-style” count normalization.
- **Status:** **uncited**
- **Finding:** `transforms/tmm.py` verifies the 205-line implementation, `TMM` class, parameters, trimming, weighting, and Numba kernel. Neither that file nor any SeqPro documentation mentions edgeR or cites the edgeR/TMM literature. The comparison to edgeR is an external equivalence claim with no supporting citation.

### 19. Only pyrefly is configured and invoked in strict mode

- **Record claim (lines 423–427):** the engineering inventory includes “basedpyright + pyrefly strict type checking.”
- **Status:** **minor-discrepancy**
- **Finding:** `pyproject.toml:66–67` explicitly gives pyrefly `preset = "strict"`, `pixi.toml:25` defines the `typecheck` task as `pyrefly check python`, and the pre-commit hook invokes that task. basedpyright is installed and has a configuration section, but it has no `typeCheckingMode = "strict"` setting and no repository task or workflow invokes it. The evidence supports strict pyrefly checking plus an available/configured basedpyright dependency, not strict checking by both tools.

## Link check

The GitHub repository and organization pages, GitHub API endpoints, PyPI JSON endpoint, documentation root, rendered Ragged page, and sibling-repository URLs were live when checked on 2026-08-14. The raw `experimental/__init__.py` URL returns HTTP 404 exactly as the record reports; that intentional negative lookup is not an unanticipated dead citation.

## Verified counts and anchors relevant to the findings

- `python/seqpro/__init__.py` exports exactly 28 names, and there are exactly 27 tracked Python files under `python/seqpro/`.
- `bed.py` is 378 lines, `gtf.py` is 1,813 bytes, and `transforms/tmm.py` is 205 lines.
- `docs/superpowers/` contains 49 Markdown files; the repository has 30 `test*.py` modules and 10 Rust `.rs` files.
- `docs/api/` contains exactly `index`, `alphabets`, `bed`, `gtf`, `ragged`, and `types`; `docs/api/index.md` lists 14 top-level members.
- PyPI contains 55 release keys, reports version 0.22.0 and Python `>=3.10`, and records the cited first/current upload dates. The dependency list in the record matches the package metadata.
- The current GitHub metadata, release count/date range, license, archive state, stars, forks, open-issue count, commit/tag dates, and sibling-repository push dates were reproduced.
- Apart from the source-attribution qualification in item 12, the displayed README and `SKILL.md` quotations were found. The `Sequential`, experimental-header, Hamming assertion, lazy-XArray, BED eager-read, translation-policy, stop-truncation, hashing, GTF, cleaner-stub, and private-symbol source anchors all matched within the stated lines or acceptable line drift.
