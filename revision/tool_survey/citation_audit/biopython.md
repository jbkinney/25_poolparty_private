# Biopython citation-integrity audit

Audited record: `final/biopython.md`  
Tool/version baseline: official Biopython `biopython-188` tag (`d7e4b8b`) plus live master where the record explicitly relies on master  
Audit date: 2026-08-14

This audit checks only whether factual evidence, quotations, locations, URLs, and counts are accurate. It does not assess any capability rating or conclusion.

## NOT FOUND IN ANY SOURCE

1. **NOT FOUND IN ANY SOURCE** — `per_sequence_provenance`, record line 109: “`annotations` ... `features` ... plus `letter_annotations` ... all round-trip to GenBank.” This is contradicted by the GenBank writer. `Bio/SeqIO/InsdcIO.py` serializes recognized GenBank annotation keys and `record.features`; it never reads `record.letter_annotations`, and arbitrary keys in `record.annotations` are not written. Thus arbitrary provenance annotations and per-letter provenance do not round-trip to GenBank.

2. **NOT FOUND IN ANY SOURCE** — `transcript_models`, record line 174: “The same exon-block model is read and written for ... `bigmaf`.” `Bio/Align/bigmaf.py:32-56` declares bigMaf as BED3 plus one `mafBlock` field (`chrom`, `chromStart`, `chromEnd`, `mafBlock`). It has no BED12 `blockCount`/`blockSizes`/`blockStarts` transcript model and no `thickStart`/`thickEnd` CDS bounds. BigMaf read/write support is real, but the claimed transcript/exon-block evidence is not.

3. **NOT FOUND IN ANY SOURCE** — `exon_intron_split_codons`, record line 182: “BED12/PSL exon blocks with `thickStart`/`thickEnd`.” Plain PSL has block coordinates, but neither `Bio/Align/psl.py` nor the PSL format has `thickStart` or `thickEnd`; those attributes occur in BED/bigBed and bigPsl. The cited PSL source does not support the stated evidence.

4. **NOT FOUND IN ANY SOURCE** — §5 item 1, record line 256: “Modern `Bio.Align` genomic file I/O — read *and write* for the full `formats` tuple.” Four members of that tuple are reader-only: `emboss`, `hhr`, `msf`, and `tabular` each define `AlignmentIterator` but no `AlignmentWriter`. `Bio.Align.write` explicitly raises `ValueError("File writing has not yet been implemented ...")` when a module lacks `AlignmentWriter` (`Bio/Align/__init__.py:4870-4875`). Sixteen of the 20 tuple members are writable, not all 20.

5. **NOT FOUND IN ANY SOURCE** — §4 assessment, record line 248: the Tutorial recipe “Producing randomized genomes” is described as “dinucleotide-preserving shuffling.” The recipe converts the sequence to a list of individual nucleotides and calls `random.shuffle(nuc_list)` (`Doc/Tutorial/chapter_cookbook.rst:119-123, 164-175, 188-190`). It preserves mononucleotide composition only; no dinucleotide-preserving algorithm appears in the recipe.

6. **NOT FOUND IN ANY SOURCE** — §6, record line 274: “Biopython no longer drives external command-line tools.” Removal of the generic `Bio.Application` wrapper framework is verified, but the broader claim is false in 1.88. Current source still launches external executables, including DSSP/mkdssp (`Bio/PDB/DSSP.py:168-195`), NACCESS (`Bio/PDB/NACCESS.py:55-68`), PSEA (`Bio/PDB/PSEA.py:28-38`), MSMS (`Bio/PDB/ResidueDepth.py:510-535`), and PAML programs (`Bio/Phylo/PAML/_paml.py:108-121`).

## Misquoted evidence

7. **misquoted** — `per_sequence_provenance`, record lines 109-110: the purported verbatim `SeqRecord.__getitem__` quotation omits its leading qualification. The actual text is: “With the exception of any molecule type, the annotations dictionary and the dbxrefs list are not used for the new SeqRecord” (`Bio/SeqRecord.py:433-437`). The record quotes only “the annotations dictionary and the dbxrefs list are not used ...” and then says slicing silently drops annotations. Slicing preserves `annotations["molecule_type"]` when present (`Bio/SeqRecord.py:572-574`), so the omission is material.

8. **misquoted** — `genome_coordinates`, record lines 143-150: the fenced block introduced by “The Tutorial does exactly this” is a rewritten excerpt, not the Tutorial code verbatim. The Tutorial uses a separate `filename = f"{name}.2bit"`, parses `filename`, asserts matching lengths, assigns the genome sequence, and then resets each record ID (`Doc/Tutorial/chapter_align.rst:1592-1599`). The record folds the filename into `SeqIO.parse(...)` and omits both the assertion and ID assignment. The central operation is supported, but the displayed code is a paraphrase presented as exact source code.

## Wrong numbers and dates

9. **number-wrong** — record line 61: “Since 1.78 it also supports undefined and partially-defined sequences.” Release 1.79 introduced `Seq(None, length=...)` and `_UndefinedSequenceData` (`NEWS.rst`, 1.79 section, lines 554 and 600-632). The `.defined` property was then added in 1.80 (`NEWS.rst:438, 482-483`). The 1.78 start version is wrong.

10. **number-wrong** — record lines 69 and 126: “19-item Tutorial cookbook.” The 1.88 Tutorial cookbook has 17 recipe subsections plus the BioSQL section, for 18 enumerated items. The record's own list at lines 238-241 also contains 13 sequence-file recipes + 4 plotting recipes + BioSQL = 18.

11. **number-wrong** — `maintained`, record line 227: “PyPI first file upload `2026-08-06T12:30:38Z`.” The earliest 1.88 file is `biopython-1.88.tar.gz` at `2026-08-06T12:13:36.003082Z`. `12:30:38.131147Z` is the first cp310 macOS wheel, not the first release file.

## Minor discrepancies

12. **minor-discrepancy** — `Bio.motifs`, record line 100: the listed methods/properties are described as the “complete documented surface.” The named members do exist, but the public `Motif` surface also includes slicing via `__getitem__`, configurable `mask`, `pseudocounts`, and `background`, and public attributes such as `alignment`, `counts`, `length`, `alphabet`, and `name` (`Bio/motifs/__init__.py:187-345`). The no-sampler/no-placer evidence is unaffected; only the claim of completeness is inaccurate.

13. **minor-discrepancy** — package inventory, record line 67: the list is labeled “Package list from `pyproject.toml`” but omits four literal entries from `[tool.setuptools].packages`: `Bio.Align.substitution_matrices.data`, `Bio.Entrez.DTDs`, `Bio.Entrez.XSDs`, and `Bio.SearchIO._model` (`pyproject.toml:45-96`). These are data/internal subpackages, but an exhaustive package-list claim should include them or state that such packages were collapsed.

14. **minor-discrepancy** — SeqIO format correction, record lines 22 and 176: the record says the old list omitted `twobit` and auto-registered AlignIO formats, but it also omitted six direct `_FormatToIterator` keys: `abi-trim`, `fasta-blast`, `fasta-pearson`, `embl-cds`, `genbank-cds`, and `sff-trim` (`Bio/SeqIO/__init__.py:414-452`). The correction remains incomplete.

15. **minor-discrepancy** — §5 item 7, record line 262: “~30 `Bio.SeqIO` formats read/write” is not reproducible from the cited registries without an unstated alias-grouping rule. In 1.88 there are 38 direct iterator keys and 18 direct writer keys; AlignIO auto-registration adds 11 iterator keys and 8 writer keys, yielding 49 readable format identifiers and 26 writable identifiers. “About 30 writable identifiers” is reasonable, but “~30 formats read/write” is ambiguous and understates the registered read surface.

16. **minor-discrepancy** — `design_visualization`, record line 120: “`pip install biopython` therefore draws nothing” is broader than its evidence. A clean environment without ReportLab cannot import `Bio.Graphics`, which is correctly documented. But `Bio.motifs.Motif.weblogo()` is a separate visualization path that POSTs to the remote WebLogo service and has no ReportLab dependency (`Bio/motifs/__init__.py:467-569`); the record itself later identifies that path. The statement is accurate only if restricted to `Bio.Graphics` local rendering.

17. **minor-discrepancy** — load-bearing absence claim, record line 69: “There is no ... function ... anywhere in Biopython concerned with designing ... a mutagenesis scheme” is too absolute. `Bio.SeqUtils.MeltingTemp.Tm_GC(valueset=2)` documents the QuikChange formula as “Recommended ... for the design of primers for QuikChange mutagenesis” (`Bio/SeqUtils/MeltingTemp.py:719-727`). This is only a Tm calculation, not a mutagenesis or primer-design operation, so the narrower functional claim elsewhere in the record remains supported.

18. **minor-discrepancy** — license quotation, record line 223: the words and omissions match `LICENSE.rst`, but it is not character-for-character verbatim. The source uses double quotation marks around “Biopython License Agreement” and “BSD 3-Clause License”; the record changes those internal marks to single quotation marks. The meaning is unchanged.

## Uncited claims, wrong locations, and links

No additional purely uncited factual claim was identified after treating the global source table and the per-capability `Source:` lines as citations. No cited file was missing, no material file:line reference was at the wrong location, and all eight URLs in §1 were live (the two `http://biopython.org/wiki/...` links redirect to live HTTPS pages). The `Bio/SeqUtils/Mapper.py` 404 is also real and is used as evidence of absence, not as a live-source citation.

## Time-scoped metadata note (verified; not a finding)

The record's “31 files” count for PyPI 1.88 was correct on its 2026-08-10 finalization date: 30 wheels plus the sdist had been uploaded by then. PyPI now lists 41 files because ten cp315/cp315t wheels were added on 2026-08-12. The count should not be treated as an extraction error; only the claimed first-upload timestamp is wrong.
