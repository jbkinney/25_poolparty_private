# Literature Review Workflow

## Goal

Round out the literature search for the PoolParty manuscript (`manuscript/poolparty_main.tex`) by systematically analyzing each paper in `lit_review/to_analyze/`, producing structured analysis files, and updating the bibliography.

## Papers to Analyze

Eight PDFs in `lit_review/to_analyze/`:

1. `Schreiber2025_tangermeme_all.pdf` — tangermeme toolkit (general-purpose sequence toolkit)
2. `Cock2009_biopython.pdf` — Biopython (general-purpose sequence toolkit)
3. `Ghazi2018_MPRADesignTools_all.pdf` — MPRA design tools (assay-specific tool)
4. `Barbon2022_VaLiAnT_all.pdf` — VaLiAnT (assay-specific tool)
5. `Zulkower2020_dnachisel.pdf` — DNA Chisel (sequence optimization tool)
6. `Georgakopoulos-Soares2017_MPRAnator.pdf` — MPRAnator (assay-specific tool)
7. `Pereira2015_pydna.pdf` — pydna (general-purpose sequence toolkit)
8. `Swainston2017_codongenie.pdf` — CodonGenie (sequence optimization tool)

All are already cited in the manuscript. The analysis will verify citation accuracy, assess novelty implications, and identify additional papers to consider.

## Citation Contexts in the Manuscript

The manuscript cites related work in several distinct contexts:

1. **General-purpose sequence toolkits** (Introduction, paragraph 3): Cock2009 (Biopython), Schreiber2025 (tangermeme), Pereira2015 (pydna), Klie2023 (SeqPro) — cited as tools that can manipulate individual sequences but lack a notion of a variant library as a structured object.

2. **Assay-specific tools** (Introduction, paragraph 3): GeorgakopoulosSoares2017 (MPRAnator), Ghazi2018 (MPRADesignTools), Barbon2022 (VaLiAnT) — cited as tools that automate library design for particular experimental formats but are not easily adapted to novel designs.

3. **Sequence optimization tools** (Introduction, paragraph 3): Zulkower2020 (DNA Chisel), Swainston2017 (CodonGenie) — cited as tools that address synthesis constraints on individual sequences but not at the level of libraries.

## PoolParty's Novelty Claims

- First unified, declarative, composable framework for DNA sequence library design
- DAG-based design with lazy evaluation (no sequences generated until requested)
- Design cards providing structured metadata for each sequence
- Automatic informative sequence naming and styling
- Extensible architecture (custom Operations via subclassing)
- Over 20 built-in Operation types covering nucleotide- and codon-level mutagenesis, motif insertion, scans, barcodes, etc.

## Analysis File Format

Each analysis file (`[pdf_file_name]_analysis.md`) follows this structure:

0. **Validated BibTeX entry** with `mynotes=` field summarizing relevance and citation location
1. **Summary** — brief description of what the paper does
2. **Relevance to PoolParty manuscript** — specific connections
3. **Citation recommendations** — where to cite, if at all
4. **Novelty assessment** — does this threaten any claims? Precautions needed?
5. **References to investigate** — papers cited by this reference worth looking into
6. **Citing papers to investigate** — papers citing this work (via Semantic Scholar)

## Per-Paper Workflow

For each PDF:

1. Read the PDF to understand its content
2. Compare with the existing BibTeX entry in `poolparty_refs.bib` — validate fields (title, authors, year, journal, DOI, pages, volume, number)
3. Write the analysis following the format above
4. Search Semantic Scholar (by DOI or title) for papers citing this work; summarize notable citers
5. Save analysis to `lit_review/analyzed/`
6. Move PDF from `to_analyze/` to `analyzed/`
7. Update the BibTeX entry in `poolparty_refs.bib` with the validated entry including `mynotes=`

## Citing-Papers Search Strategy

**Primary tool: Semantic Scholar API** (via web fetch of `https://api.semanticscholar.org/graph/v1/paper/DOI:{doi}?fields=citations.title,citations.year,citations.authors,citations.citationCount,citations.abstract`)

- Semantic Scholar provides structured data (title, abstract, year, citation count) and is programmatically accessible
- Google Scholar is not reliably accessible programmatically (CAPTCHA blocks)
- For very recent preprints where Semantic Scholar coverage may be thin, supplement with general web searches

**Known limitation:** Semantic Scholar can lag behind Google Scholar by weeks/months for newly published papers. This mainly affects `Schreiber2025_tangermeme_all.pdf`. Any gaps will be noted in the analysis.
