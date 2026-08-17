# Biopython — FINAL corrected capability record

**Tool name:** Biopython
**Tool slug:** `biopython`
**Citation key:** `Cock2009df` — Cock PJA *et al.* (2009) "Biopython: freely available Python tools for computational molecular biology and bioinformatics." *Bioinformatics* 25(11):1422–1423. doi:10.1093/bioinformatics/btp163
**Tier:** 1
**Category:** general-purpose bioinformatics library (sequence/alignment I/O + manipulation + genomic-coordinate arithmetic). **Not** a library-design tool.
**Record finalized:** 2026-08-10 (extraction + adversarial review + independent re-verification of both contested cells)
**Assessed version:** 1.88 (released 2026-08-06), tag `biopython-188`; “current” below means that tag unless an explicit retrieval date is given

---

## 0. What changed from the original extraction

The adversarial review upheld 22 of 24 cells verbatim and flagged two as **understated**, both in the same direction: the extraction assessed Biopython as though `Bio.Align` were still the 2010-era `AlignIO` wrapper. I independently re-verified both, plus every artifact the reviewer added, against tag `biopython-188` over `raw.githubusercontent.com`.

| Cell | Was | Now | Why |
|---|---|---|---|
| `genome_coordinates` | partial | **yes** | The extraction's evidence asserted three absences (no chromosome-name resolution against a reference, no coordinate-system conversion/liftover, no `hg38 chr7:1000-2000` entry point). All three are false as of 1.88 and I confirmed the refutations in source. Corrected, not split. |
| `transcript_models` | partial | **partial** (value unchanged; evidence rewritten) | The GTF/GFF-absent fact is correct, but the reasoning generalized illegitimately from the `Bio.SeqIO` format table. BED12/bigBed/bigPsl exon-block transcript models with `thickStart`/`thickEnd` CDS bounds — and PSL exon blocks without CDS bounds — are read **and written** natively via `Bio.Align`. |

Also corrected: the quoted `Bio.SeqIO` format list (omitted `twobit`, six other direct `_FormatToIterator` keys, and the auto-registered `AlignIO` formats); the `additional_capabilities` inventory (missed the entire modern `Bio.Align` genomic-file stack); `design_visualization` (ReportLab caveat that *weakens* the cell); `availability_status` (ReportLab; second CVE); `documented_examples` (second, user-contributed wiki cookbook enumerated).

No cell was set to `unknown` — the one place extractor and reviewer differed on value (`genome_coordinates`) is a factual question with a verified answer, not a judgment call. Borderline judgment calls are flagged in §6 rather than hidden in an `unknown`.

---

## 1. Sources consulted

| Kind | Reference |
|---|---|
| pdf | `papers/Cock2009_biopython.pdf` — Cock et al. (2009) *Bioinformatics* 25(11):1422–1423 |
| repo | https://github.com/biopython/biopython (tag `biopython-188`, re-verified 2026-08-14) |
| repo | `Bio/Align/__init__.py` (`map`, `mapall`, `formats` tuple, `CodonAligner`) — re-read |
| repo | `Bio/Align/bigbed.py` (`search`), `Bio/Align/bed.py`, `Bio/Align/chain.py`, `Bio/Align/analysis.py` — re-read |
| repo | `Bio/SeqIO/TwoBitIO.py`, `Bio/SeqIO/__init__.py` (`_FormatToIterator`, AlignIO auto-registration) — re-read |
| repo | `Bio/SeqRecord.py` (`__add__`, `__getitem__`), `Bio/Seq.py` (`MutableSeq`, `defined`, `defined_ranges`) — re-read |
| repo | `Bio/Graphics/__init__.py` (ReportLab guard), `Bio/SeqUtils/__init__.py`, `Bio/Emboss/Primer3.py`, `Bio/Restriction/__init__.py` |
| repo | `pyproject.toml` (package list, deps, license, absence of `[project.scripts]`), `NEWS.rst`, `DEPRECATED.rst`, `LICENSE.rst` |
| docs | https://biopython.org/docs/latest/Tutorial/index.html (Tutorial & Cookbook, 1.88) |
| docs | https://biopython.org/docs/latest/Tutorial/chapter_align.html (liftover + twobit worked example) |
| docs | https://biopython.org/docs/latest/api/Bio.html (API index) |
| docs | http://biopython.org/wiki/GFF_Parsing |
| docs | http://biopython.org/wiki/Category:Cookbook (second, user-contributed cookbook; 16 entries) |
| pypi | https://pypi.org/pypi/biopython/json |
| web | https://biopython.org/ (homepage) |

**Nothing was installed or executed.** All source read over HTTPS; metadata from the PyPI JSON API.

---

## 2. What Biopython is

The 2009 applications note is 2 pages and describes ~1.49. Its own CONCLUSIONS section calls it *"a large open-source application programming interface (API) used in both bioinformatics software development and in everyday scripts"* and warns *"The features described herein are only a subset; potential users should refer to the tutorial and API documentation for further information."* Everything below is therefore anchored on **1.88 source and docs**, not the paper.

Core abstractions:

- `Bio.Seq.Seq` — immutable sequence with `translate()`, `reverse_complement()`, `transcribe()`. Since 1.79 it also supports **undefined and partially-defined** sequences: `Seq(None, length=248956422)`, `.defined` (added in 1.80), `.defined_ranges` — a chromosome can be represented by name + length with no bases in memory.
- `Bio.Seq.MutableSeq` — *"An editable sequence object. … the MutableSeq lets you edit the sequence in place"* (`Bio/Seq.py:2173`). Per-position assignment; the raw primitive a user would loop over to make variants.
- `Bio.SeqRecord.SeqRecord` — sequence + `id`/`name`/`description` + `annotations` dict + `features` list + `letter_annotations`. Supports slicing (shifts in-range features) and **concatenation with annotation bookkeeping** (`Bio/SeqRecord.py:915 __add__`).
- `Bio.SeqFeature` — `SimpleLocation` / `CompoundLocation` (join locations) + `.extract()` (trans-splicing since 1.78).
- `Bio.SeqIO` / `Bio.AlignIO` / **`Bio.Align`** — the last is the actively developed one and now carries a full UCSC genomic-format stack (see §5).

Package list from `pyproject.toml` (`[tool.setuptools] packages`; data and internal subpackages — `Bio.Align.substitution_matrices.data`, `Bio.Entrez.DTDs`, `Bio.Entrez.XSDs`, `Bio.SearchIO._model` — collapsed): Bio, Affy, Align (+substitution_matrices), AlignIO, Alphabet (an import-failing compatibility stub; removed in 1.78), Blast, CAPS, Cluster, codonalign, Compass, Data, Emboss, Entrez, ExPASy, GenBank, Geo, Graphics(+GenomeDiagram), HMM, KEGG(+Compound/Enzyme/Gene/Map/KGML), Medline, motifs(+jaspar), Nexus, NMR, Pathway(+Rep), PDB(+mmtf), phenotype, PopGen(+GenePop only; FDist and SimCoal were removed in 1.70), Restriction, SCOP, SearchIO(+BlastIO/HHsuiteIO/HmmerIO/**InfernalIO**/ExonerateIO/InterproscanIO), SeqIO, SeqUtils, Sequencing, SVDSuperimposer, SwissProt, TogoWS, Phylo(+PAML), UniGene, UniProt, BioSQL.

**There is no module, class, function, tutorial chapter, or cookbook recipe in Biopython that generates a variant library, a mutagenesis scheme, an oligo pool, or a combinatorial construct.** Verified against the 1.88 API index, the full 27-chapter Tutorial TOC, the 18-item Tutorial cookbook, and the separate 16-entry wiki cookbook. This is the load-bearing claim for the comparison table and it survived both review passes.

---

## 3. Capability-by-capability record

### Block A — library specification

**`library_as_object` = no**
The only multi-sequence containers are (a) a plain Python `list`/generator of `SeqRecord`s, (b) `Bio.Align.Alignment` / `Bio.Align.Alignments` / `MultipleSeqAlignment` — *alignments*: only the legacy `MultipleSeqAlignment` stores equal-length gapped rows, a modern `Alignment` stores the original (possibly unequal-length) sequences plus a coordinate array, and `Alignments` is a container of those; the 2009 paper describes `Bio.SeqIO` as interpreting multiple-sequence-alignment formats as "collections of equal length (gapped) sequences" — and (c) the in-memory dict eagerly built by `SeqIO.to_dict()` together with the dict-like `index()` / `index_db()` views over a file. None carries design semantics: no wild-type reference, no variant set, no record of the operation that produced a member. Building a library means writing your own loop.
*Source:* `pyproject.toml` package list; `Bio/Align/__init__.py` class list; `Bio.SeqIO` API page; Cock 2009 p.1423. (Review note: the extraction missed `Bio.Align.Alignments`; it is a multi-alignment container, still not a design object.)

**`dag_chaining` = no**
No pipeline/graph/workflow/composition abstraction anywhere in `Bio/`. Chaining is ordinary Python function composition written by the user. `Bio.Application` — the only layer that ever resembled a pipeline (chaining external command-line tools) — was *"Declared obsolete in release 1.79, deprecated in release 1.82, and removed in release 1.86."*
*Source:* `DEPRECATED.rst:228–232` (verified verbatim, twice); API index.

**`lazy_evaluation` = partial**
Lazy *reading* is real, first-class, and now extends to sub-record granularity:
- `SeqIO.parse()` returns an iterator; `SeqIO.index()` gives a dict-like object where *"When you access a particular record via the dictionary methods, the code will jump to the appropriate part of the file and then parse that section into a SeqRecord"*; `index_db()` persists that index in SQLite. Cookbook recipe: "Indexing a FASTQ file".
- **Region-level laziness**: `Bio/SeqIO/TwoBitIO.py` module docstring — *"Using the information in the index, the `__getitem__` method calculates the file position at which the requested region starts, and only reads the requested sequence region. Note that the full sequence of a record is loaded only if specifically requested, making the parser memory-efficient. The TwoBitIterator object implements the `__getitem__`, keys, and `__len__` methods that allow it to be used as a dictionary."* (verified verbatim). `Bio/SeqIO/__init__.py:637` explicitly documents "lazy-loading file formats such as twobit".
- `Seq(None, length)` / `.defined_ranges` let a record carry chromosome-scale coordinates with no sequence loaded at all.
- `Bio.bgzf` gives random access into BGZF-compressed files.

**Partial, not yes**, because this is lazy *consumption of sequences that already exist on disk*. There is no lazy *generation* of designed sequences, because there is no design step to defer.
*Source:* `Bio/SeqIO/TwoBitIO.py` lines 7–22 (read directly); `Bio/SeqIO/__init__.py:450, 637`; Tutorial `chapter_seqio.html` ("Sequence files as Dictionaries"), `chapter_cookbook.html`.

**`mixed_mutagenesis_one_pool` = no**
Biopython has no mutagenesis operation of any kind — not exhaustive single, double, random, or WT-replicate; a fortiori no mixing of mutagenesis modes in one pool. Nothing in the API index, the Tutorial, or either cookbook. The closest primitives are `MutableSeq`'s in-place edits — item and slice assignment/deletion, `insert`, `append`, `extend`, `pop`, `remove`, plus the shared `replace`, `join`, `+` and `*`: *"An editable sequence object. … the MutableSeq lets you edit the sequence in place"*.
*Source:* `Bio/Seq.py:2173` (verbatim); full Tutorial + both cookbook TOCs.

**`combinatorial_motif_place` = no**
`Bio.motifs` is analysis-only. Its documented surface includes `pwm`, `pssm`, `consensus`, `anticonsensus`, `degenerate_consensus`, `relative_entropy`, `reverse_complement()`, `weblogo()`, `format()`, plus parsers (MEME, TRANSFAC, JASPAR, …) and PSSM searching. There is **no** function to sample a sequence from a PWM and **no** function to insert or place a motif into a host sequence, hence no combinatorial placement machinery.
*Wording correction from review:* the original memo said insertion "would be hand-written string surgery on a MutableSeq". That undersells the available primitive — `SeqRecord` slicing plus concatenation **does** shift feature coordinates correctly (`Bio/SeqRecord.py:915 __add__`, documented with the pBAD30 plasmid re-linearization doctest: `new = right + left`, `len(new.features) == len(left.features) + len(right.features)`, annotations merged when unambiguous). Insertion is `rec[:i] + insert + rec[i:]` with features maintained — but only those falling entirely within a slice survive, so a feature spanning position `i` is lost. Still no motif sampler, no placement combinatorics, no library.
*Source:* `Bio/motifs/__init__.py` method list; `Bio/SeqRecord.py:915–975` (read directly).

**`barcode_generation` = no**
No barcode/UMI/index generator. `Bio/SeqUtils/` contains exactly `CheckSum`, `IsoelectricPoint`, `MeltingTemp`, `ProtParam`, `ProtParamData`, `lcc`, `__init__` — verified; the `Bio/SeqUtils/Mapper.py` referenced by a stale wiki cookbook entry is absent from the 1.88 repository tree. Constraint primitives a user could build a generator from exist (`gc_fraction()`, `MeltingTemp`, `Bio.Align.PairwiseAligner` for distance, `lcc` for complexity), but there is no generator, no edit-distance-constrained set construction, and no attach-to-library step.
*Source:* `Bio/SeqUtils/` directory listing at tag `biopython-188`; official wiki cookbook.

**`per_sequence_provenance` = partial**
The container exists and is idiomatic: `SeqRecord.annotations` is *"A dictionary of additional information about the sequence"*, plus a structured `features` list of `SeqFeature`s with locations and qualifiers, plus `letter_annotations` for per-position data; `features` and the recognised INSDC `annotations` keys round-trip to GenBank, but arbitrary `annotations` keys and `letter_annotations` are not written by the GenBank writer (`Bio/SeqIO/InsdcIO.py` never reads `letter_annotations`). **But** Biopython never constructs a designed sequence, so nothing is ever populated automatically — provenance is 100% user-written. And the container leaks: `SeqRecord.__getitem__` docstring states *"With the exception of any molecule type, the annotations dictionary and the dbxrefs list are not used for the new SeqRecord"*, so slicing silently drops annotations apart from `molecule_type` (features in range and letter_annotations are kept). Concatenation is better behaved — `__add__` merges annotations when unambiguous and *"for any ambiguities (e.g. different names) it defaults to omitting that annotation."*
*Source:* `Bio/SeqRecord.py` `__getitem__` and `__add__` docstrings (both verified verbatim); `Bio/SeqIO/InsdcIO.py` (GenBank writer).

**`automatic_naming` = no**
`SeqRecord.id` / `.name` / `.description` are passive slots: user-supplied or read from a file. No variant-naming scheme of any kind — no HGVS-style, position-based, or operation-derived name generation anywhere in the package.
*Source:* `Bio/SeqRecord.py`; API index; Cock 2009 p.1422 (SeqRecord description).

**`design_visualization` = partial**
`Bio.Graphics.GenomeDiagram` is a genuine, well-documented visualization layer — *"designed for drawing whole genomes, in particular prokaryotic genomes, either as linear diagrams (optionally broken up into fragments to fit better) or as circular wheel diagrams"* — with stacked tracks, feature sets, sigils (BOX/ARROW/OCTO/JAGGY/BIGARROW), cross-links between tracks, and PDF/EPS/SVG/PNG output; plus `BasicChromosome` and `motifs.weblogo()`. Cookbook adds "Nucleotide dot plots", "Histogram of sequence lengths", "Plot of sequence GC%". So *annotated sequences* can be drawn well.
Two things cut against the cell:
1. There is no **design** view — no library/pool rendering, no operation-graph view, no highlighting of designed variants against a reference, because those objects do not exist.
2. **It does not work out of the box.** `Bio/Graphics/__init__.py` docstring is *"Bio.Graphics offers several graphical outputs, all using ReportLab"* and the module raises `MissingPythonDependencyError("Please install ReportLab if you want to use Bio.Graphics…")` on import without it — while PyPI `requires_dist` is `["numpy"]` only. `pip install biopython` therefore renders no `Bio.Graphics` output.
*Source:* `Bio/Graphics/__init__.py:9–24` (read directly); Tutorial `chapter_graphics.html`; PyPI JSON `requires_dist`.

### Block B — assay coverage

**`assay_dms` = no**
No deep-mutational-scanning or saturation-mutagenesis functionality, and no such term in the 27-chapter Tutorial TOC, the 18-item Tutorial cookbook, or the separate 16-entry wiki cookbook (all three enumerated).
*Source:* Tutorial index; `chapter_cookbook.html`; `biopython.org/wiki/Category:Cookbook`.

**`assay_mpra` = no**
No regulatory/MPRA library design. `Bio.motifs` is TFBS *analysis* (PWM/PSSM scanning and parsing), not reporter-library construction; no tiling, no element scrambling, no promoter/enhancer library builder in either cookbook.
*Source:* `Bio/motifs/__init__.py` surface; both cookbooks.

**`assay_insilico` = no**
No model-probing or in-silico perturbation design; no sequence-to-function model interface of any kind. Biopython's only ML/statistical content is `Bio.Cluster` and `Bio.PopGen` — classical methods unrelated to genomic deep-learning models; the `Bio.HMM` implementation modules (`DynamicProgramming`, `Trainer`, `MarkovModel`, `Utilities`), together with `Bio.kNN`, `Bio.LogisticRegression`, `Bio.NaiveBayes` and `Bio.MaxEntropy`, were removed in 1.86, leaving `Bio.HMM` an empty package.
*Source:* API index; Cock 2009 p.1422 feature list (re-extracted from PDF, matches verbatim).

### Block C — genomics integration

**`genome_coordinates` = yes** *(corrected from `partial`; this is the one value change)*

Biopython handles genomic coordinates against real assemblies at three documented levels:

1. **Chromosome-name resolution against a reference, lazily.** `twobit` (UCSC `.2bit`, the format UCSC ships hg38 in) is a registered `Bio.SeqIO` format — `_FormatToIterator["twobit"] = TwoBitIO.TwoBitIterator` at `Bio/SeqIO/__init__.py:450` — and the parser is dict-like and region-lazy (docstring quoted under `lazy_evaluation`). The Tutorial does exactly this for six assemblies:
   ```python
   names = ("panTro6", "hg38", "rheMac10", "calJac4", "mm39", "rn7")
   for i, name in enumerate(names):
       filename = f"{name}.2bit"
       genome = SeqIO.parse(filename, "twobit")
       chromosome = genome_alignment.sequences[i].id
       assert len(genome_alignment.sequences[i]) == len(genome[chromosome])
       genome_alignment.sequences[i] = genome[chromosome]
       genome_alignment.sequences[i].id = f"{name}.{chromosome}"
   ```
   So `SeqIO.parse("hg38.2bit","twobit")["chr7"][1000:2000]` is a one-liner that reads only that byte range. The original extraction's claim that this "requires the user to load the FASTA and slice it themselves" is false.
2. **Assembly-to-assembly liftover via UCSC chain files.** `Bio/Align/__init__.py:3396` (`Alignment.map`): *"The map method can also be used to lift over an alignment between different genome assemblies. In this case, self is a DNA alignment between two genome assemblies, and the argument is an alignment of a transcript against one of the genome assemblies"* — with a doctest reading `panTro5ToPanTro6.over.chain` and printing the mapped coordinate arrays. `Bio/Align/chain.py` is a complete read/write implementation of the UCSC chain format, and `"chain"` is in the `Bio.Align.formats` tuple (`Bio/Align/__init__.py:4818`). The Tutorial works a 6-species liftover with `mapall()` over `hg19ToHg38.chain`, `mm10ToMm39.chain`, `panTro5ToPanTro6.chain`, `rheMac8ToRheMac10.chain`, `calJac3ToCalJac4.chain`, `rn6ToRn7.chain`.
3. **Indexed region query.** `Bio/Align/bigbed.py:1042`:
   ```python
   def search(self, chromosome=None, start=None, end=None):
       """Iterate over alignments overlapping the specified chromosome region..
       This method searches the index to find alignments to the specified
       chromosome that fully or partially overlap the chromosome region
       between start and end."""
   ```
   (verified verbatim, including the stray double period). This is an R-tree index lookup, not a linear scan, and the same method is exposed on bigBed, bigPsl and bigMaf — Tutorial: `alignments.search("chr3", 48000000, 49000000)`.

Assembly identity is also carried in the data model: MAF/bigMaf records are assembly-qualified (`hg19.chr1`, `mm10.chr3`), and `Seq(None, length=…)` lets a chromosome be a pure name+length coordinate frame (`SeqRecord(Seq(None, length=1575), id="chr1")` in the Tutorial, to supply a PSL target size).

**Honest riders** (these are why the cell is not simply "best in class"): there is no built-in assembly registry or reference-download service — the user supplies the `.2bit` and `.chain` files; and, as everywhere in this record, there is no *design* object for coordinates to be attached to.
*Source:* `Bio/SeqIO/TwoBitIO.py`; `Bio/SeqIO/__init__.py:450`; `Bio/Align/__init__.py:3339 map`, `:3396`, `:3584 mapall`, `:4815–4832 formats`; `Bio/Align/bigbed.py:1042`; `Bio/Align/chain.py`; Tutorial `chapter_align.html`. All re-read directly at finalization.

**`transcript_models` = partial** *(value unchanged; evidence replaced)*

Two halves, and the original memo only had the first:

*Absent:* GTF/GFF is **not** in core. `pyproject.toml` lists no GFF package; the project's own wiki is explicit — *"GFF parsing is not yet integrated into Biopython. This documentation is work towards making it ready for inclusion"*. GTF/GFF support therefore remains external to core.

*Present:* Exon-block transcript models **are** read and written natively via `Bio.Align`. `Bio/Align/bed.py` module docstring: *"The Browser Extensible Data (BED) format, stores a series of pairwise alignments in a single file. Typically they are used for transcript to genome alignments."* The reader parses `blockCount`/`blockSizes`/`blockStarts` into `Alignment.coordinates` (i.e. exons) and the writer emits `thickStart`/`thickEnd` — the UCSC CDS boundaries ("start codon"/"stop codon"). The same exon-block model is read and written for `bigbed` and `bigpsl` — and, without `thickStart`/`thickEnd`, for `psl` — and parsed by `Bio.SearchIO.ExonerateIO`. GenBank/EMBL `mRNA`/`CDS` features with `join(...)` locations are parsed into `CompoundLocation` — a transcript model in all but file format.

*Correction to the original memo's format table:* it listed the `Bio.SeqIO` formats as "abi, ace, cif-atom, cif-seqres, embl, fasta, fasta-2line, fastq(+variants), gck, genbank/gb, gfa1, gfa2, ig, imgt, nib, pdb-seqres, pdb-atom, phd, pir, seqxml, sff, snapgene, swiss, tab, qual, uniprot-xml, xdna". That list is incomplete: it omits **`twobit`** (`Bio/SeqIO/__init__.py:450`), the six other direct `_FormatToIterator` keys `abi-trim`, `fasta-blast`, `fasta-pearson`, `embl-cds`, `genbank-cds` and `sff-trim`, and the `AlignIO` formats auto-registered into `SeqIO` at `Bio/SeqIO/__init__.py:509–519`.

*Why `partial` and not `yes`:* no GTF/GFF, and no transcript-model *object* with named exons/UTRs/CDS that a design tool could target ("place this variant in exon 4 of ENST00000…"). What exists is the alignment-coordinate representation, which is sufficient arithmetic but not a gene-model API.
*Source:* `Bio/Align/bed.py:7–12` and `:82–218` (read directly); `Bio/Align/__init__.py:4815–4832`; `Bio/SeqIO/__init__.py:414–519`; `pyproject.toml` package list; `biopython.org/wiki/GFF_Parsing`.

**`exon_intron_split_codons` = partial**
The mechanism is correct and documented: `CompoundLocation` represents joined exons, and *"The `SeqFeature` object has an `extract` method to take care of all this (and since Biopython 1.78 can handle trans-splicing by supplying a dictionary of referenced sequences)"* — exons are concatenated in transcript order with strand handled, so a codon split across an exon boundary translates correctly via `feature.extract(genome).translate()`. Cookbook: "Translating a FASTA file of CDS entries". A second route to the same arithmetic exists through `Bio.Align`: BED12 exon blocks with `thickStart`/`thickEnd`, plus `Bio.Align.CodonAligner` (confirmed present in `Bio/Align/__init__.py`) for codon-aware alignment.
**But** because Biopython performs no mutagenesis, it never has to *place* an edit into a codon that straddles an intron — the hard part of this capability in a design tool is simply never exercised.
*Source:* Tutorial `chapter_seq_annot.html`; `Bio/SeqFeature.py`; `Bio/Align/__init__.py` (`CodonAligner`); `Bio/Align/bed.py`.

**`hgvs_input` = no**
No HGVS module, class, parser, or Tutorial section. The full `Bio/` package listing in `pyproject.toml` contains nothing HGVS-related.
*Source:* `pyproject.toml` package list; API index; Tutorial TOC.

**`vcf_vep_output` = no**
No `Bio.VCF`/variant module; VCF is not a `Bio.SeqIO` or `Bio.AlignIO` or `Bio.Align` format. No VEP interoperability of any kind.
*Wording correction from review:* the original phrase "writes sequence and alignment formats only" is now inaccurate — Biopython also writes genomic-interval formats (BED, PSL, SAM, bigBed, chain, MAF). The precise statement is: **Biopython writes no variant-call formats.**
*Source:* `Bio/SeqIO/__init__.py` `_FormatToWriter`; `Bio/Align/__init__.py:4815–4832`; API index.

**`consequence_annotation` = no**
No molecular-consequence classifier. `Seq.translate()` + `Bio.Data.CodonTable` let a user compare two translations by hand, but there is no stop-gained / missense / synonymous / frameshift / in-frame-indel annotation function. The nearest thing in the package is `Bio/Align/analysis.py calculate_dn_ds(alignment, method="NG86"|"LWL85"|"YN00"|"ML")`, which partitions synonymous vs non-synonymous *sites* across a codon alignment — a population-genetics statistic, not a per-variant consequence call. (`Bio.CAPS` is restriction-site polymorphism detection in an alignment, a different thing again.)
*Source:* `Bio/Align/analysis.py` (read directly); `Bio/CAPS/`; API index.

### Block D — physical construction

**`primer_design` = no**
`Bio.Emboss.Primer3` is explicitly a **parser**: module docstring *"Code to parse output from the EMBOSS eprimer3 program."* `Bio.Emboss.PrimerSearch` likewise writes input for / parses output from EMBOSS primersearch. `Bio/Emboss/` contains exactly `Primer3.py`, `PrimerSearch.py`, `__init__.py`. The wrappers that used to *launch* eprimer3 (`Bio.Emboss.Applications`, built on `Bio.Application`) were **removed in 1.86**. `Bio.SeqUtils.MeltingTemp` computes Tm — including `Tm_GC(valueset=2)`, documented as the QuikChange formula recommended for mutagenesis-primer design — and the Cookbook has "Trimming off primer sequences" (trimming, not design). The formula accepts an already supplied sequence and returns only a temperature; it does not select a primer. Net: Biopython reads someone else's primer designs and scores a Tm; it designs no primers and no mutagenic oligos.
*Source:* `Bio/Emboss/Primer3.py` docstring (verbatim); `Bio/Emboss/` listing; `Bio/SeqUtils/MeltingTemp.py:719–727`; `DEPRECATED.rst`.

**`codon_optimization` = yes (basic — see caveat, this is the second-most-debatable cell)**
`Bio.SeqUtils.CodonAdaptationIndex.optimize()` exists and is documented: *"Return a new DNA sequence with preferred codons only. Uses the codon adaptiveness table defined by the CodonAdaptationIndex object to generate DNA sequences with only preferred codons. May be useful when designing DNA sequences for transgenic protein expression or codon-optimized proteins like fluorophores."* It accepts DNA/RNA/protein input and a `strict` flag for ties; `calculate()` scores CAI.
**Caveats that must travel with the cell:** the CAI table is learned from a *user-supplied* set of coding sequences (`CodonAdaptationIndex.__init__(sequences, table=standard_dna_table)`) — no host presets ship, and the old `CodonUsageIndices.py` / `SharpEcoliIndex` were removed. The implementation is literally `"".join(pref_codons[aa] for aa in aa_seq)`: replace every codon with the single most-preferred one. No codon harmonization, no usage-frequency sampling, no GC-window / repeat / restriction-site / secondary-structure constraints.
*Source:* `Bio/SeqUtils/__init__.py:578–725` (read in full by both extractor and reviewer; quote verbatim).

**`synthesis_constraints` = partial**
`Bio.Restriction` is a first-class restriction-site analyser over the ~1000-enzyme REBASE set — `Analysis(AllEnzymes, seq)`, `.print_that()`, `.blunt()`, enzyme-set filtering, doctest on the pBluescript MCS including *"Enzymes which do not cut the sequence"* — which provides unwanted-site checking. Add `gc_fraction()`, `MeltingTemp`, `molecular_weight()`, `lcc`, `nt_search()`.
**But** there is no synthesis-constraint feature as such: no homopolymer or repeat rules, no GC-window enforcement, no oligo-length or vendor limits, no "reject or repair this sequence" API. The user assembles checks from primitives, and nothing enforces anything on a library (there being no library object to enforce it on).
*Source:* `Bio/Restriction/__init__.py` doctest (verified verbatim); `Bio/SeqUtils/__init__.py`.

### Block E — engineering

**`interface` = yes** — *Python library API only: no CLI, no GUI, no web service.*
The bare value `yes` means "has a usable interface"; the modality matters and must be carried into the table cell. `pyproject.toml` declares **no `[project.scripts]`** and no console entry points (grep-verified), so `pip install biopython` installs zero command-line tools. The repo's `Scripts/` directory (GenBank, PDB, Performance, Restriction, SeqGui, `query_pubmed.py`, `scop_pdb.py`, `update_ncbi_codon_table.py`, xbbtools) holds contributed demos — including two Tk GUI demos, `SeqGui` and `xbbtools` — that are **not installed by the package**. No hosted web service. The paper's CONCLUSIONS section is consistent: *"a large open-source application programming interface (API) used in both bioinformatics software development and in everyday scripts."*
*Caveat to the "offline library" reading:* `Bio.motifs.weblogo()` POSTs to the remote WebLogo web service, and `Bio.Entrez`/`ExPASy`/`KEGG`/`TogoWS`/`UniProt` are network clients. Biopython consumes web services; it does not expose one.
*Source:* `pyproject.toml`; `Scripts/`; Cock 2009 p.1423 CONCLUSIONS (re-extracted from PDF, verbatim match).

**`license` = yes (open, permissive)**
`pyproject.toml`: `license = "LicenseRef-Biopython-License-Agreement"`, `license-files = ["LICENSE.rst"]`. `LICENSE.rst`: *"Biopython is currently released under the 'Biopython License Agreement' … Some files are explicitly dual licensed under your choice of the 'Biopython License Agreement' or the 'BSD 3-Clause License' … with the intention of later offering all of Biopython under this dual licensing approach."* The Biopython License Agreement is an MIT/BSD-style permissive grant. GitHub's licence detector reports NOASSERTION only because the text is custom.
*Source:* `pyproject.toml:11–12`; `LICENSE.rst`.

**`maintained` = yes (very actively)**
Latest release **1.88, 2026-08-06** — PyPI first file upload `2026-08-06T12:13:36.003082Z` (the sdist), 41 files, `requires-python >=3.10`, `requires_dist ["numpy"]` (re-verified against the PyPI JSON on 2026-08-14). Wheels cover cp310–cp315 including free-threaded builds, manylinux x86_64 + aarch64, macOS, and Windows. `NEWS.rst` confirms *"6 August 2026: Biopython 1.88"* and an open *"(In progress, not yet released): Biopython 1.89"* section. Release cadence: 1.83 (2024-01-10), 1.84 (2024-06-28), 1.85 (2025-01-15), 1.86 (2025-10-28), 1.87 (2026-03-30), 1.88 (2026-08-06). Ships `py.typed`; type annotations expanding in 1.88.
**Security:** two consecutive security releases — 1.87 *"Addressed security issue CVE-2025-68463 in `Bio.Entrez.Parser`"*; 1.88 removed use of Python's built-in `eval` reachable via the `Bio.Nexus` NEXUS parser (arbitrary code execution), following an email report.
*Minor caveat:* on 2026-08-14, biopython.org's homepage advertised 1.87, i.e. the homepage lagged PyPI and the 1.88 Tutorial by one release.
*Source:* PyPI JSON and official homepage (re-queried 2026-08-14); `NEWS.rst` lines 13–52 (read directly).

---

## 4. Biopython's own documented examples

Biopython ships a **Tutorial & Cookbook** (27 chapters: Introduction; Quick Start; Sequence objects; Sequence annotation objects; Sequence Input/Output; Sequence alignments; Pairwise sequence alignment; Multiple Sequence Alignment objects; pairwise2; BLAST (new); BLAST (old); BLAST and other search tools; Entrez; Swiss-Prot and ExPASy; Bio.PDB; Bio.PopGen; Bio.Phylo; Bio.motifs; Cluster analysis; Graphics incl. GenomeDiagram; KEGG; Bio.phenotype; Cookbook; testing framework; contributing; Python appendix; Bibliography).

**Tutorial cookbook recipes** (the closest thing to "example libraries"):
*Working with sequence files* — Filtering a sequence file; Producing randomized genomes; Translating a FASTA file of CDS entries; Making the sequences in a FASTA file upper case; Sorting a sequence file; Simple quality filtering for FASTQ files; Trimming off primer sequences; Trimming off adaptor sequences; Converting FASTQ files; Converting FASTA and QUAL files into FASTQ files; Indexing a FASTQ file; Converting SFF files; Identifying open reading frames.
*Sequence parsing plus simple plots* — Histogram of sequence lengths; Plot of sequence GC%; Nucleotide dot plots; Plotting the quality scores of sequencing read data.
*BioSQL* — storing sequences in a relational database.

**Second, user-contributed cookbook** at `biopython.org/wiki/Category:Cookbook` (16 entries, enumerated at finalization so the "I checked the full cookbook" claim is airtight): Plotting ABI traces; Represent an alignment from contig archived in ACE files; Retrieve and annotate Entrez Gene IDs; Concatenating multiple alignments NEXUS files; Converting sequence files with Bio.SeqIO; Mapping genetic coordinates with the Bio.SeqUtils.Mapper module; Methods for Degenerated Codons; From gene sequence to predicted protein with the GFF module; Workflow to extract intergenic regions; Bio.Phylo Cookbook; Reading data from UNIX pipes; Reading large PDB files; Remove PDB disordered atoms; Retrieve nonmatching blast queries; Sequence Cleaner; Split large file.
Two of these are stale or external: `Bio.SeqUtils.Mapper` does not exist in the 1.88 repository tree, and "From gene sequence to predicted protein with the GFF module" relies on software outside Biopython core.

**Third, doctests as worked examples**: the liftover walkthrough in `chapter_align.html` (6-species MAF + UCSC chain files + `.2bit` genomes) and the pBAD30 plasmid re-linearization doctest in `Bio/SeqRecord.py` are the two most design-adjacent things in the corpus.

**Assessment for the referee response:** none of these is a *library design* example. Not one recipe produces a designed set of variant sequences. The closest adjacents are "Producing randomized genomes" (`random.shuffle` over the list of individual nucleotides of an existing genome, preserving mononucleotide composition only — an analysis control, not a designed library), "Translating a FASTA file of CDS entries", "Trimming off primer sequences", "Methods for Degenerated Codons" (*n*-fold degeneracy *analysis*), and the pBAD30 concatenation doctest (construct assembly of two fragments, not a library). **There is nothing to reproduce head-to-head.** The honest framing for the referee is: Biopython supplies primitives PoolParty builds on; it is infrastructure, not a competing design tool.

---

## 5. Notable capabilities outside the row list

*(Corrected and expanded per the review — the original inventory read like 2010-era Biopython.)*

1. **Modern `Bio.Align` genomic file I/O — read for the full `formats` tuple, write for 16 of its 20 members** (`Bio/Align/__init__.py:4812–4833`): `a2m, bed, bigbed, bigmaf, bigpsl, chain, clustal, emboss, exonerate, fasta, hhr, maf, mauve, msf, nexus, phylip, psl, sam, stockholm, tabular` — `emboss`, `hhr`, `msf` and `tabular` define no `AlignmentWriter`, and `Bio.Align.write` raises `ValueError` for them. This is the newest and most actively developed part of the package and the original extraction omitted it entirely. **Candidate for a "genomic-interval / alignment format interoperability" row where Biopython is the ceiling.**
2. **Genome-assembly liftover** via UCSC chain files: `Alignment.map()` / `Alignment.mapall()`, documented with hg19→hg38, mm10→mm39, panTro5→panTro6, rheMac8→rheMac10, calJac3→calJac4, rn6→rn7.
3. **Indexed genomic-region query**: `alignments.search(chromosome, start, end)` on bigBed/bigPsl/bigMaf (`Bio/Align/bigbed.py:1042`) — R-tree index lookup.
4. **`twobit` lazy genome reader** in `Bio.SeqIO` — chromosome-keyed dict access with byte-range reads of only the requested region.
5. **Undefined / partially-defined sequences**: `Seq(None, length)`, `.defined`, `.defined_ranges` — chromosome-scale coordinates with no sequence in memory.
6. **`SeqRecord` arithmetic with annotation bookkeeping**: `record_a + record_b` merges annotations and shifts feature locations (pBAD30 doctest). The real construct-assembly primitive.
7. **Sequence file I/O breadth** — 49 readable `Bio.SeqIO` format identifiers and 26 writable (GenBank, EMBL, FASTQ variants and xDNA read/write; SnapGene, GCK, GFA1/2 and **twobit** read-only; plus auto-registered AlignIO formats).
8. **Restriction-enzyme analysis** over the full REBASE set (`Bio.Restriction.Analysis`), blunt/sticky/isoschizomer filtering.
9. **Alignment engine**: `Bio.Align.PairwiseAligner` (C-accelerated global/local/affine, substitution matrices), `Bio.Align.CodonAligner`, `Bio.Align.analysis.calculate_dn_ds` (NG86/LWL85/YN00/ML), plus `Bio.SearchIO`/`Bio.Blast` parsers (Infernal parsing added in 1.86).
10. **Remote database clients**: `Bio.Entrez` (incl. coordinate-limited `efetch(seq_start=, seq_stop=)`), `Bio.ExPASy`, `Bio.KEGG`, `Bio.UniProt`, `Bio.TogoWS`; and `Bio.motifs.weblogo()` POSTs to the remote WebLogo service.
11. **3D structure**: `Bio.PDB` (PDB/mmCIF parsing, superposition, DSSP/NACCESS interfaces).
12. **Phylogenetics** (`Bio.Phylo`, `Bio.Nexus`), **population genetics** (`Bio.PopGen.GenePop` parsing/manipulation only; the paper's FDist and SimCoal modules were removed in 1.70), **clustering** (`Bio.Cluster`), **protein physicochemistry** (`Bio.SeqUtils.ProtParam`, `IsoelectricPoint`), **checksums** (`Bio.SeqUtils.CheckSum`), **BGZF random access** (`Bio.bgzf`), **BioSQL ORM**.
13. **Typed**: ships `py.typed`; annotations expanding in 1.88.

---

## 6. Stated limitations, unresolved disagreements, and confidence

**Stated limitations (the tool's own).** The 2009 paper concedes its scope: *"The features described herein are only a subset; potential users should refer to the tutorial and API documentation for further information."* The project's wiki concedes the GFF gap: *"GFF parsing is not yet integrated into Biopython. This documentation is work towards making it ready for inclusion."* `Bio.Application` (external-tool wrappers, incl. everything that launched EMBOSS eprimer3) was declared obsolete in 1.79, deprecated in 1.82, and **removed in 1.86** — Biopython no longer ships a generic framework for driving external command-line tools, though individual modules still launch executables directly via `subprocess` (`Bio.PDB.DSSP`, `NACCESS`, `PSEA`, `ResidueDepth`/MSMS, `Bio.Phylo.PAML`). `Bio.Graphics` hard-requires ReportLab, which the package does not depend on. Two security fixes in consecutive releases (CVE-2025-68463 in `Bio.Entrez.Parser`, 1.87; `eval()`-based RCE in the `Bio.Nexus` parser, 1.88) indicate a parser attack surface that is being actively hardened.

**The overarching limitation for this matrix, stated plainly:** Biopython has no design layer. There is no library object, no mutagenesis operation, no barcode generator, no motif placer, no HGVS or VCF interface, no consequence classifier, and no primer designer. Every `no` in Block A/B and most of Block C/D trace to that single absence rather than to twenty independent gaps — the table should say so once rather than repeating the reasoning per cell.

**Unresolved disagreements — flagged, none hidden:**

1. **`genome_coordinates`: resolved against the extractor.** Extractor said `partial` on the grounds that chromosome-name resolution, coordinate conversion, and liftover do not exist. Reviewer said those three claims are false; I re-read `Bio/Align/__init__.py:3396`, `Bio/Align/bigbed.py:1042`, `Bio/SeqIO/TwoBitIO.py`, `Bio/SeqIO/__init__.py:450`, and `Bio/Align/__init__.py:4815–4832` myself and the reviewer is right on every point. Value corrected to **`yes`**. This is not a judgment call and is not marked `unknown`. **Do not ship the original evidence sentence** because the cited 1.88 code directly contradicts it.
2. **`transcript_models`: value agreed (`partial`), reasoning replaced.** Both parties reach `partial`; the extractor's route there (GTF/GFF absent ⇒ no transcript models) is invalid, since BED12/bigBed/PSL/bigPsl exon+CDS models are read and written natively. Value kept, evidence rewritten. A referee could argue `yes` on the exon-block evidence or `no` on the strict GTF/GFF reading; `partial` is the defensible middle **only with the new evidence attached**.
3. **`codon_optimization` = `yes` vs `partial (naive)`: genuine judgment call, kept at `yes`.** The reviewer verified the `optimize()` docstring verbatim and marked the cell supported, so `yes` is quotable. But the implementation is one line of most-preferred-codon substitution over a user-supplied CAI table. **Consistency risk:** if any other tool in the matrix receives `partial` for a richer optimizer (harmonization, frequency sampling, constraint-aware), Biopython at `yes` will look wrong. Either ship the docstring quote in the cell or downgrade this and re-baseline the column. Flagging rather than silently splitting.
4. **`synthesis_constraints` = `partial`** hinges on counting `Bio.Restriction` site-scanning as a synthesis check. A stricter reading gives `no`. Both extractor and reviewer accept `partial`; note it as the third-most-debatable cell.
5. **`interface` = `yes`** is uninformative as a bare token for a key whose interesting content is "library only". The table cell must carry the modality: *Python library API only — no CLI, no GUI, no hosted web service.*

**Method caveat:** nothing was installed or executed at any stage. All statements derive from published documentation, the PyPI JSON API, the 2009 PDF (text extracted with PyMuPDF), and tag `biopython-188` source read over HTTPS. Quoted strings marked "verbatim" were matched character-for-character against the source at least twice by independent passes.

---

## 7. Bottom line

Biopython is infrastructure, not a competitor. Of 24 capability keys it scores `yes` on four (`genome_coordinates`, `codon_optimization`, plus the two engineering descriptors `license` and `maintained`; `interface` is a fifth `yes` describing modality rather than power), `partial` on six, and `no` on the rest. It is maximally available and maintained — release 1.88 on 2026-08-06, wheels for Python 3.10–3.15 on every major platform, one dependency. That makes it the safe, uncontroversial anchor of the "general-purpose toolkit" row.

The one thing this record must get right for the referee response is the corrected coordinate story: **Biopython does chromosome-name resolution, region-lazy genome access, indexed region query, and UCSC chain-file liftover, and it does them well — it just has nothing to design.** That is a stronger and far more defensible framing than the original memo's, and it concedes exactly what a Biopython author would insist on conceding.
