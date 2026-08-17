# Biopython — v2 capability record (revised row set)

**Tool slug:** `biopython`
**Tool name:** Biopython
**Citation key:** `Cock2009df` — Cock PJA *et al.* (2009) *Bioinformatics* 25(11):1422–1423. doi:10.1093/bioinformatics/btp163
**Tier:** 1
**Category:** general-purpose bioinformatics library (sequence/alignment I/O, sequence manipulation, genomic-coordinate arithmetic). **Not** a library-design tool.
**Assessed version:** 1.88 (released 2026-08-06), master-branch source
**v2 scored:** 2026-08-10, from `final/biopython.md` (extraction + adversarial review + independent re-verification) plus five targeted re-checks listed in §4.

**Scope rule applied.** Biopython is a general-purpose toolkit; that is a description, not a criticism. Only functionality that could participate in *specifying a DNA sequence library* is scored. Out of scope and deliberately excluded from every cell: `Bio.PDB`, `Bio.Phylo`, `Bio.PopGen`, `Bio.Cluster`, `Bio.KEGG`, `Bio.Entrez`/`ExPASy`/`TogoWS`/`UniProt` clients, `Bio.SearchIO`/`Bio.Blast` parsers, `Bio.SeqUtils.ProtParam`/`IsoelectricPoint`/`CheckSum`, dN/dS statistics, phylogenetic and clustering machinery, and pure file-format breadth. Where a capability exists only as raw primitives the user must assemble, the cell is `partial` at most, never `yes`.

**The one structural fact behind most cells.** Biopython has no design layer: no library object, no mutagenesis operation, no barcode generator, no motif placer, no variant-consequence classifier, no primer designer. Verified against the 1.88 API index, the full 27-chapter Tutorial TOC, the Tutorial cookbook chapter, and the separate 16-entry wiki cookbook. Most `no` values in Blocks A, B and D trace to that single absence rather than to twenty independent gaps.

---

## 1. Value table (v2)

### Block A — library specification

| Key | Value |
|---|---|
| `library_first_class_object` | **partial** |
| `composable_operations` | **no** |
| `lazy_generation` | **partial** |
| `library_algebra` | **no** |
| `exhaustive_single_scans` | **no** |
| `sampled_random_mutagenesis` | **no** |
| `higher_order_combinatorial` | **no** |
| `heterogeneous_components_one_library` | **no** |
| `combinatorial_motif_place` | **no** |
| `barcode_generation` | **no** |
| `per_sequence_provenance` | **partial** |
| `automatic_naming` | **no** |
| `design_visualization` | **partial** |

### Block B — assay coverage

| Key | Value |
|---|---|
| `assay_dms` | **no** |
| `assay_mpra` | **no** |
| `assay_insilico` | **no** |

### Block C — genomics integration

| Key | Value |
|---|---|
| `genome_coordinates` | **yes** |
| `transcript_models` | **partial** |
| `exon_intron_split_codons` | **partial** |
| `vcf_vep_output` | **no** |
| `consequence_annotation` | **no** |

### Block D — adjacent / complementary

| Key | Value |
|---|---|
| `primer_design` | **no** |
| `codon_optimization` | **yes** (basic) |
| `synthesis_constraints` | **partial** |
| `degenerate_iupac_codons` | **partial** |
| `negative_control_generation` | **no** |
| `ml_model_in_loop` | **no** |
| `readout_analysis` | **partial** |

### Block E — engineering and availability

| Key | Value | Real answer |
|---|---|---|
| `interface` | yes | Python library API only — no CLI, no installed GUI, no hosted web service |
| `license` | yes | Biopython License Agreement (MIT/BSD-style permissive); parts dual-licensed BSD-3-Clause |
| `installable_today` | yes | `pip install biopython` → 1.88, wheels cp310–cp314 (incl. free-threaded cp314t) on manylinux x86_64/aarch64, macOS arm64, win32/win_amd64; one dependency (`numpy`) |
| `last_activity` | yes | Release **1.88 on 2026-08-06**, four days before this survey; `NEWS.rst` has an open "1.89 in progress" section |
| `documented_examples` | yes | 27-chapter Tutorial & Cookbook + ~17-recipe cookbook chapter + separate 16-entry wiki cookbook + doctests — **none of them a library-design example** |

---

## 2. Evidence, cell by cell

### Block A

**`library_first_class_object` = partial**
*Changed from v1 `library_as_object` = no.* The split separates two claims that v1 fused, and Biopython satisfies one of them. Biopython is an API, not a file-writing program: `Seq`, `MutableSeq`, `SeqRecord`, `SeqFeature`, `Bio.Align.Alignment`, `Bio.Align.Alignments` and `MultipleSeqAlignment` are all genuine in-memory objects the user holds, inspects, transforms (slice, concatenate, translate, reverse-complement, sort) and passes onward to writers. `SeqRecord.__add__` even shifts feature coordinates and merges annotations (`Bio/SeqRecord.py:915`, pBAD30 re-linearization doctest: `new = right + left`, `len(new.features) == len(left.features) + len(right.features)`).
What does **not** exist is a *library* object. The only multi-record containers are (a) a plain Python `list`/generator of `SeqRecord`s, (b) alignment objects — equal-length gapped rows; `Bio.SeqIO` "interprets multiple sequence alignment file formats as collections of equal length (gapped) sequences" — and (c) the dict-like `SeqIO.to_dict()` / `index()` / `index_db()` views over a file. None carries a wild-type reference, a variant set, or a record of the operation that produced a member. Building a library means writing your own loop.
*Source:* `Bio/SeqRecord.py:915–975` (read directly); `Bio/Align/__init__.py` class list; `pyproject.toml` package list; `Bio.SeqIO` API page.

**`composable_operations` = no**
*Rename of v1 `dag_chaining` = no; value carried.* There is no pipeline, graph, workflow, or operation-composition abstraction anywhere in `Bio/`. The only layer that ever resembled one — `Bio.Application`, which chained external command-line tools — was *"Declared obsolete in release 1.79, deprecated in release 1.82, and removed in release 1.86."* The new wording ("design steps compose/chain/nest") does not change the answer, because there are no design steps to compose. Biopython's primitives do compose freely, but as ordinary Python function composition written by the user, with no tool-level semantics carried through the chain.
*Source:* `DEPRECATED.rst:228–232` (verbatim, verified twice); 1.88 API index.

**`lazy_generation` = partial**
*Rename of v1 `lazy_evaluation` = partial; value carried.* Lazy, non-materializing access is real, first-class, and extends to sub-record granularity:
- `SeqIO.parse()` returns an iterator; `SeqIO.index()` yields a dict-like object where *"When you access a particular record via the dictionary methods, the code will jump to the appropriate part of the file and then parse that section into a SeqRecord"*; `index_db()` persists the index in SQLite.
- Region-level laziness: `Bio/SeqIO/TwoBitIO.py` — *"Using the information in the index, the `__getitem__` method calculates the file position at which the requested region starts, and only reads the requested sequence region. Note that the full sequence of a record is loaded only if specifically requested, making the parser memory-efficient."* (verbatim). `Bio/SeqIO/__init__.py:637` names "lazy-loading file formats such as twobit".
- `Seq(None, length)` / `.defined_ranges` carry chromosome-scale coordinates with no bases in memory; `Bio.bgzf` gives random access into BGZF files.

`partial` and not `yes` because this is lazy *consumption of sequences that already exist on disk*. Nothing is *generated* on demand, because there is no generative design step to defer.
*Source:* `Bio/SeqIO/TwoBitIO.py:7–22`; `Bio/SeqIO/__init__.py:450, 637`; Tutorial `chapter_seqio.html`, `chapter_cookbook.html`.

**`library_algebra` = no**
*New half of the v1 `library_as_object` split.* No operation combines, samples, or repeats whole libraries. The nearest in-tool set operations are alignment-level and were re-checked in source for this pass: `MultipleSeqAlignment.append()` (*"Add one more SeqRecord object to the alignment as a new row"*) and `.extend()` (*"Add more SeqRecord objects to the alignment as rows. They must all have the same length as the original alignment"*) stack rows; `aln1 + aln2` concatenates **by column** (*"Combine two alignments with the same number of rows by adding them… Using the extend method adds by row. Using the addition operator adds by column"*); `.sort()` reorders rows. All require equal-length gapped rows, i.e. an alignment, not a library. **There is no random-sampling method and no replication/repeat method on any alignment or record container** (verified over `Bio/Align/__init__.py` for this pass). Combining two designed sets is plain-Python `list1 + list2` — host-language glue outside the tool's semantics — which the row definition scores as `no`.
*Source:* `Bio/Align/__init__.py` (`MultipleSeqAlignment.append`/`extend`/`__add__`/`sort` docstrings, re-read 2026-08-10); wiki cookbook "Concatenating multiple alignments NEXUS files" (file-level alignment concatenation).

**`exhaustive_single_scans` = no**
*One of four rows replacing v1 `mixed_mutagenesis_one_pool` = no.* Biopython has no mutagenesis operation of any kind, so a fortiori no exhaustive single-substitution or single-deletion scan. Nothing in the API index, the Tutorial, or either cookbook. The closest primitive is per-position assignment on a `MutableSeq` — *"An editable sequence object. … the MutableSeq lets you edit the sequence in place"* — i.e. the raw thing a user would loop over themselves.
*Source:* `Bio/Seq.py:2173` (verbatim); full Tutorial TOC + both cookbook TOCs.

**`sampled_random_mutagenesis` = no**
No random or rate-based variant sampler. The only randomization in the corpus is the Tutorial cookbook recipe "Producing randomized genomes", and on re-reading it for this pass it is explicitly **not** a Biopython feature: *"I would use the built-in Python `random` module for this, in particular the function `random.shuffle`"* — the recipe converts the sequence to a list, calls `random.shuffle`, and rebuilds a `Seq`. That is a composition-preserving shuffled control, hand-written in user code, not rate-based mutagenesis. No mutation-rate parameter, no per-position sampling, no RNG-seeded variant generator anywhere in `Bio/`.
*Source:* Tutorial `chapter_cookbook.html` "Producing randomized genomes" (re-fetched 2026-08-10, quote verbatim); API index.

**`higher_order_combinatorial` = no**
No pairwise or higher-order mutation combination, because there is no mutagenesis operation to combine. No cartesian-product, k-of-n, or multi-edit-per-sequence machinery. `SeqRecord` slicing + concatenation lets a user apply two independent edits to one record by hand (`rec[:i] + a + rec[i:j] + b + rec[j:]`), but the combinatorics, the enumeration and the bookkeeping are entirely user code.
*Source:* API index; `Bio/SeqRecord.py:915–975`; both cookbook TOCs.

**`heterogeneous_components_one_library` = no**
Vacuously absent: there is no library specification, hence no notion of a component within one, hence no mixing of structurally different component types. The nearest object, `MultipleSeqAlignment`, requires the opposite — all rows the same length.
*Source:* `Bio/Align/__init__.py` (`append`/`extend` length constraint); `pyproject.toml` package list.

**`combinatorial_motif_place` = no**
`Bio.motifs` is analysis-only. Its complete documented surface is `pwm`, `pssm`, `consensus`, `anticonsensus`, `degenerate_consensus`, `relative_entropy`, `reverse_complement()`, `weblogo()`, `format()`, plus parsers (MEME, TRANSFAC, JASPAR, …) and PSSM searching. There is **no** function to sample a sequence from a PWM and **no** function to insert or place a motif into a host sequence, hence no placement combinatorics over positions, orientations or permutations. The available primitive is honest to state: insertion *is* `rec[:i] + insert + rec[i:]` and feature coordinates are maintained correctly — but that is one insertion, written by the user.
*Source:* `Bio/motifs/__init__.py` method list; `Bio/SeqRecord.py:915`.

**`barcode_generation` = no**
No barcode/UMI/index generator. `Bio/SeqUtils/` contains exactly `CheckSum`, `IsoelectricPoint`, `MeltingTemp`, `ProtParam`, `ProtParamData`, `lcc`, `__init__` (verified; the wiki-referenced `Bio/SeqUtils/Mapper.py` returns HTTP 404 — never merged). Constraint primitives a user could build a generator from exist (`gc_fraction()`, `MeltingTemp`, `PairwiseAligner` for distance, `lcc` for complexity), but there is no generator, no edit-distance-constrained set construction, and no attach-to-library step.
*Source:* `Bio/SeqUtils/` directory listing; direct HTTP check of `Bio/SeqUtils/Mapper.py` (404).

**`per_sequence_provenance` = partial**
The container exists and is idiomatic: `SeqRecord.annotations` is *"A dictionary of additional information about the sequence"*, plus a structured `features` list of `SeqFeature`s with locations and qualifiers, plus `letter_annotations` for per-position data; all round-trip to GenBank. **But** Biopython never constructs a designed sequence, so nothing is ever populated automatically — provenance is 100% user-written. The container also leaks: `SeqRecord.__getitem__` states *"the annotations dictionary and the dbxrefs list are not used for the new SeqRecord"*, so slicing silently drops annotations (in-range features and `letter_annotations` survive). Concatenation is better behaved: `__add__` merges annotations and *"for any ambiguities (e.g. different names) it defaults to omitting that annotation."*
*Source:* `Bio/SeqRecord.py` `__getitem__` and `__add__` docstrings (both verbatim).

**`automatic_naming` = no**
`SeqRecord.id` / `.name` / `.description` are passive slots — user-supplied or read from a file. No variant-naming scheme of any kind: no HGVS-style, position-based, or operation-derived name generation anywhere in the package.
*Source:* `Bio/SeqRecord.py`; API index; Cock 2009 p.1422.

**`design_visualization` = partial**
`Bio.Graphics.GenomeDiagram` is a genuine, well-documented visualization layer — *"designed for drawing whole genomes… either as linear diagrams (optionally broken up into fragments to fit better) or as circular wheel diagrams"* — with stacked tracks, feature sets, sigils (BOX/ARROW/OCTO/JAGGY/BIGARROW), cross-links between tracks, and PDF/EPS/SVG/PNG output; plus `BasicChromosome` and `motifs.weblogo()`. Cookbook adds nucleotide dot plots, length histograms, GC% plots. So *annotated sequences* draw well.
Two things cut against the cell: (1) there is no **design** view — no library/pool rendering, no operation graph, no highlighting of designed variants against a reference, because those objects do not exist; (2) it does not work out of the box — `Bio/Graphics/__init__.py` raises `MissingPythonDependencyError("Please install ReportLab if you want to use Bio.Graphics…")`, while PyPI `requires_dist` is `["numpy"]` only, so `pip install biopython` draws nothing.
*Source:* `Bio/Graphics/__init__.py:9–24`; Tutorial `chapter_graphics.html`; PyPI JSON `requires_dist`.

### Block B

**`assay_dms` = no** — No deep-mutational-scanning or saturation-mutagenesis functionality; the term appears nowhere in the 27-chapter Tutorial TOC, the cookbook chapter, or the 16-entry wiki cookbook (all three enumerated).
*Source:* Tutorial index; `chapter_cookbook.html`; `biopython.org/wiki/Category:Cookbook`.

**`assay_mpra` = no** — No regulatory/MPRA library design. `Bio.motifs` is TFBS *analysis* (PWM/PSSM scanning and parsing), not reporter-library construction; no tiling, no element scrambling, no promoter/enhancer library builder in either cookbook.
*Source:* `Bio/motifs/__init__.py` surface; both cookbooks.

**`assay_insilico` = no** — No model-probing or in-silico perturbation design, and no sequence-to-function model interface of any kind. Biopython's only ML/statistical content is the legacy `Bio.HMM`, `Bio.Cluster`, `Bio.PopGen` — classical methods unrelated to genomic deep-learning models.
*Source:* API index; Cock 2009 p.1422 feature list.

### Block C

**`genome_coordinates` = yes**
Biopython handles coordinates against real assemblies at three levels, all documented:
1. **Chromosome-name resolution against a reference, lazily.** `twobit` (UCSC `.2bit`, the format UCSC ships hg38 in) is a registered `Bio.SeqIO` format (`_FormatToIterator["twobit"]`, `Bio/SeqIO/__init__.py:450`) and the parser is dict-like and region-lazy, so `SeqIO.parse("hg38.2bit","twobit")["chr7"][1000:2000]` reads only that byte range. The Tutorial does this for six assemblies (`panTro6, hg38, rheMac10, calJac4, mm39, rn7`).
2. **Assembly-to-assembly liftover via UCSC chain files.** `Bio/Align/__init__.py:3396` (`Alignment.map`): *"The map method can also be used to lift over an alignment between different genome assemblies…"*, with a doctest over `panTro5ToPanTro6.over.chain`. `Bio/Align/chain.py` is a complete read/write chain implementation and `"chain"` is in the `Bio.Align.formats` tuple. The Tutorial works a 6-species liftover via `mapall()`.
3. **Indexed region query.** `Bio/Align/bigbed.py:1042` `search(chromosome, start, end)` — *"This method searches the index to find alignments to the specified chromosome that fully or partially overlap the chromosome region between start and end"* — an R-tree lookup, exposed on bigBed, bigPsl and bigMaf.
Assembly identity is carried in the data model (MAF/bigMaf records are `hg19.chr1`, `mm10.chr3`), and `Seq(None, length=…)` makes a chromosome a pure name+length coordinate frame.
*Honest riders:* no built-in assembly registry or reference-download service (the user supplies `.2bit` / `.chain`), and no design object for coordinates to attach to.
*Source:* `Bio/SeqIO/TwoBitIO.py`; `Bio/SeqIO/__init__.py:450`; `Bio/Align/__init__.py:3339, 3396, 3584, 4815–4832`; `Bio/Align/bigbed.py:1042`; `Bio/Align/chain.py`; Tutorial `chapter_align.html`.

**`transcript_models` = partial**
*Absent:* GTF/GFF is not in core. The project wiki: *"GFF parsing is not yet integrated into Biopython. This documentation is work towards making it ready for inclusion"*, directing users to `BCBio.GFF` or gffutils.
*Present:* exon-block transcript models are read **and written** natively via `Bio.Align`. `Bio/Align/bed.py`: *"The Browser Extensible Data (BED) format, stores a series of pairwise alignments in a single file. Typically they are used for transcript to genome alignments."* The reader parses `blockCount`/`blockSizes`/`blockStarts` into `Alignment.coordinates` (exons); the writer emits `thickStart`/`thickEnd` (UCSC CDS bounds). Same model for `psl`, `bigpsl`, `bigbed`, `bigmaf`, plus `Bio.SearchIO.ExonerateIO`. GenBank/EMBL `mRNA`/`CDS` features with `join(...)` parse to `CompoundLocation`.
*Why `partial`:* no GTF/GFF, and no transcript-model *object* with named exons/UTRs/CDS that a design step could target ("place this variant in exon 4 of ENST…"). What exists is sufficient coordinate arithmetic, not a gene-model API.
*Source:* `Bio/Align/bed.py:7–12, 82–218`; `Bio/Align/__init__.py:4815–4832`; `Bio/SeqIO/__init__.py:414–519`; `biopython.org/wiki/GFF_Parsing`.

**`exon_intron_split_codons` = partial**
The mechanism is correct and documented: `CompoundLocation` represents joined exons and *"The `SeqFeature` object has an `extract` method to take care of all this (and since Biopython 1.78 can handle trans-splicing by supplying a dictionary of referenced sequences)"* — exons concatenate in transcript order with strand handled, so a codon straddling an exon boundary translates correctly via `feature.extract(genome).translate()`. Cookbook: "Translating a FASTA file of CDS entries". A second route exists through `Bio.Align`: BED12/PSL exon blocks with `thickStart`/`thickEnd`, plus `Bio.Align.CodonAligner`.
**But** because Biopython performs no mutagenesis, it never has to *place* an edit into a codon that straddles an intron — the hard part of this capability in a design tool is never exercised.
*Source:* Tutorial `chapter_seq_annot.html`; `Bio/SeqFeature.py`; `Bio/Align/__init__.py` (`CodonAligner`); `Bio/Align/bed.py`.

**`vcf_vep_output` = no**
No `Bio.VCF` or variant module; VCF is not a `Bio.SeqIO`, `Bio.AlignIO` or `Bio.Align` format; no VEP interoperability. Stated precisely (the v1 phrasing was too broad — Biopython *does* write genomic-interval formats BED/PSL/SAM/bigBed/chain/MAF): **Biopython writes no variant-call formats.**
*Source:* `Bio/SeqIO/__init__.py` `_FormatToWriter`; `Bio/Align/__init__.py:4815–4832`.

**`consequence_annotation` = no**
No molecular-consequence classifier. `Seq.translate()` + `Bio.Data.CodonTable` let a user compare two translations by hand, but there is no stop-gained / missense / synonymous / frameshift / in-frame-indel annotation function. The nearest thing is `Bio/Align/analysis.py calculate_dn_ds(alignment, method="NG86"|"LWL85"|"YN00"|"ML")`, which partitions synonymous vs non-synonymous *sites* across a codon alignment — a population-genetics statistic, not a per-variant consequence call. (`Bio.CAPS` is restriction-site polymorphism detection in an alignment, different again.)
*Source:* `Bio/Align/analysis.py`; `Bio/CAPS/`; API index.

### Block D

**`primer_design` = no**
`Bio.Emboss.Primer3` is explicitly a **parser**: *"Code to parse output from the EMBOSS eprimer3 program."* `Bio.Emboss.PrimerSearch` likewise writes input for / parses output from EMBOSS primersearch. `Bio/Emboss/` contains exactly `Primer3.py`, `PrimerSearch.py`, `__init__.py`. The wrappers that used to *launch* eprimer3 (`Bio.Emboss.Applications`, built on `Bio.Application`) were removed in 1.86. `Bio.SeqUtils.MeltingTemp` computes Tm; the cookbook has "Trimming off primer sequences" (trimming, not design). Net: Biopython reads someone else's primer designs and scores a Tm; it designs no primers and no mutagenic oligos.
*Source:* `Bio/Emboss/Primer3.py` docstring (verbatim); `Bio/Emboss/` listing; `DEPRECATED.rst`.

**`codon_optimization` = yes (basic)**
`Bio.SeqUtils.CodonAdaptationIndex.optimize()` exists and is documented: *"Return a new DNA sequence with preferred codons only. Uses the codon adaptiveness table defined by the CodonAdaptationIndex object to generate DNA sequences with only preferred codons. May be useful when designing DNA sequences for transgenic protein expression or codon-optimized proteins like fluorophores."* It accepts DNA/RNA/protein input and a `strict` flag for ties; `calculate()` scores CAI.
**Caveats that must travel with the cell:** the CAI table is learned from a *user-supplied* set of coding sequences (`__init__(sequences, table=standard_dna_table)`) — no host presets ship, and `CodonUsageIndices.py` / `SharpEcoliIndex` were removed. The implementation is literally `"".join(pref_codons[aa] for aa in aa_seq)`: replace every codon with the single most-preferred one. No harmonization, no usage-frequency sampling, no GC-window / repeat / restriction-site / secondary-structure constraints.
*Source:* `Bio/SeqUtils/__init__.py:578–725` (read in full by extractor and reviewer; quote verbatim).

**`synthesis_constraints` = partial**
`Bio.Restriction` is a first-class restriction-site analyser over the ~1000-enzyme REBASE set — `Analysis(AllEnzymes, seq)`, `.print_that()`, `.blunt()`, enzyme-set filtering, doctest on the pBluescript MCS including *"Enzymes which do not cut the sequence"* — covering the most common synthesis/cloning check (unwanted sites). Add `gc_fraction()`, `MeltingTemp`, `molecular_weight()`, `lcc`, `nt_search()`.
**But** there is no synthesis-constraint feature as such: no homopolymer or repeat rules, no GC-window enforcement, no oligo-length or vendor limits, no reject-or-repair API. The user assembles checks from primitives, and nothing enforces anything on a library.
*Source:* `Bio/Restriction/__init__.py` doctest (verbatim); `Bio/SeqUtils/__init__.py`.

**`degenerate_iupac_codons` = partial** *(new row; re-checked in source for this pass)*
Biopython ships real IUPAC/degenerate-codon machinery, but built for *interpreting* degenerate input rather than *designing* a degenerate library:
- `Bio.Data.IUPACData` defines `ambiguous_dna_values` (`M:AC, R:AG, W:AT, S:CG, Y:CT, K:GT, V:ACG, H:ACT, D:AGT, B:CGT, X:GATC, N:GATC`), `ambiguous_rna_values`, `ambiguous_dna_complement`, and the ambiguous/unambiguous letter alphabets — i.e. the expansion table.
- `Bio.Data.CodonTable` provides `AmbiguousCodonTable`, `AmbiguousForwardTable`, `ambiguous_dna_by_name`/`by_id` (+ RNA and generic variants) and two public helpers: `list_possible_proteins` (degenerate codon → the amino acids it can encode) and `list_ambiguous_codons`, docstring *"Extend a codon list to include all possible ambiguous codons. e.g. `['TAG', 'TAA'] -> ['TAG', 'TAA', 'TAR']`"* — that is genuine codon-set → degenerate-codon compression.
- `Seq.translate()` handles ambiguous codons through these tables; `Bio.SeqUtils.nt_search()` and `Bio.Restriction` expand ambiguous patterns for searching; `Bio.motifs.degenerate_consensus` compresses an aligned motif to an IUPAC string.
**Not `yes`:** there is no degenerate-codon *design* function — no NNK/NNS/NNB or reduced-codon scheme, no "pick the minimal IUPAC codon covering this amino-acid set" solver, no expansion of a degenerate template into the member sequences of a library, and no attachment of any of this to a library specification. Everything above is a table or a translation helper the user must assemble into a design.
*Source:* `Bio/Data/IUPACData.py` and `Bio/Data/CodonTable.py` (both re-read 2026-08-10, quotes verbatim); `Bio/motifs/__init__.py`; wiki cookbook "Methods for Degenerated Codons" (*n*-fold degeneracy *analysis*).

**`negative_control_generation` = no** *(new row)*
There is no control-generation feature and, on re-checking, **no shuffle or scramble function anywhere in `Bio/`**. The Tutorial cookbook's "Producing randomized genomes" recipe explicitly delegates to the standard library — *"I would use the built-in Python `random` module for this, in particular the function `random.shuffle`"* — so even the one randomized-control example is user code, not a Biopython capability; there is no dinucleotide-preserving or *k*-mer-preserving shuffle. `Seq.reverse_complement()` and `Seq.complement()` do exist as first-class methods, but they are generic strand operations used throughout the package, not offered as control generation, and nothing attaches controls to a designed set (there being no designed set).
*Source:* Tutorial `chapter_cookbook.html` "Producing randomized genomes" (re-fetched 2026-08-10, quote verbatim); `Bio/Seq.py`; API index.

**`ml_model_in_loop` = no** *(new row)*
No predictive model participates in anything. There is no sequence-to-function model interface, no scoring/optimization loop, no gradient or search over sequences against a model objective, and no dependency on any ML framework (PyPI `requires_dist` is `["numpy"]`). Biopython's only model-flavoured content is the legacy `Bio.HMM`, `Bio.Cluster` and `Bio.PopGen` — classical statistics applied to existing data, not design.
*Source:* API index; PyPI JSON `requires_dist`; `pyproject.toml` package list.

**`readout_analysis` = partial** *(new row)*
Biopython has substantial, documented sequencing-read handling, but nothing that closes the loop against a built library:
- `Bio.SeqIO` parses FASTQ in all three quality encodings with per-base `letter_annotations["phred_quality"]`; cookbook recipes cover "Simple quality filtering for FASTQ files", "Trimming off primer sequences", "Trimming off adaptor sequences", "Indexing a FASTQ file", "Converting FASTQ files", "Plotting the quality scores of sequencing read data"; `Bio.bgzf` and `index_db()` handle large read files; `Bio.Align` reads SAM; `PairwiseAligner` can match reads to references.
**Not `yes`:** there is no demultiplexing by barcode, no counting of reads against a designed variant set, no count matrix, no enrichment/effect-size estimation, and no concept of "the library that was built" for reads to be assigned to. These are raw primitives the user assembles, which caps the cell at `partial`.
*Source:* Tutorial `chapter_cookbook.html` recipe list (re-fetched 2026-08-10); `Bio/SeqIO/QualityIO.py`; `Bio/Align/sam.py`; `Bio/bgzf.py`.

### Block E

**`interface` = yes** — *Python library API only: no CLI, no GUI, no web service.*
`pyproject.toml` declares **no `[project.scripts]`** and no console entry points (grep-verified), so `pip install biopython` installs zero command-line tools. The repo's `Scripts/` directory (including two Tk GUI demos, `SeqGui` and `xbbtools`) holds contributed demos that are **not installed by the package**. No hosted web service. The paper's CONCLUSIONS agrees: *"a large open-source application programming interface (API) used in both bioinformatics software development and in everyday scripts."*
*Caveat to an "offline library" reading:* `Bio.motifs.weblogo()` POSTs to the remote WebLogo service, and `Bio.Entrez`/`ExPASy`/`KEGG`/`TogoWS`/`UniProt` are network clients. Biopython consumes web services; it does not expose one.
*Source:* `pyproject.toml`; `Scripts/`; Cock 2009 p.1423.

**`license` = yes** — Biopython License Agreement, an MIT/BSD-style permissive grant.
`pyproject.toml`: `license = "LicenseRef-Biopython-License-Agreement"`, `license-files = ["LICENSE.rst"]`. `LICENSE.rst`: *"Biopython is currently released under the 'Biopython License Agreement' … Some files are explicitly dual licensed under your choice of the 'Biopython License Agreement' or the 'BSD 3-Clause License' … with the intention of later offering all of Biopython under this dual licensing approach."* GitHub's detector reports NOASSERTION only because the text is custom.
*Source:* `pyproject.toml:11–12`; `LICENSE.rst`.

**`installable_today` = yes** *(new row)* — `pip install biopython`, cleanly, today.
PyPI 1.88: first file upload `2026-08-06T12:30:38Z`, 31 files, `requires-python >=3.10`, `requires_dist ["numpy"]`. Wheels cover cp310–cp314 **including free-threaded cp314t**, manylinux x86_64 + aarch64, macosx_11_0_arm64, win32 + win_amd64 — i.e. no compiler needed on any mainstream platform, and one runtime dependency. Ships `py.typed`.
*The one install caveat:* `Bio.Graphics` hard-requires ReportLab, which is not a declared dependency, so plotting needs a second install.
*Source:* PyPI JSON (re-queried 2026-08-10); `pyproject.toml`.

**`last_activity` = yes** *(new row; supersedes v1 `maintained`)* — release **1.88 on 2026-08-06**, four days before this survey.
`NEWS.rst`: *"6 August 2026: Biopython 1.88"* plus an open *"(In progress, not yet released): Biopython 1.89"* section. Cadence: 1.83 (2024-01-10), 1.84 (2024-06-28), 1.85 (2025-01-15), 1.86 (2025-10-28), 1.87 (2026-03-30), 1.88 (2026-08-06). Last commit to master observed at extraction time: `2026-08-06T14:20:26Z` ("Post Biopython 1.88 release version bump"). Security posture is active: 1.87 fixed CVE-2025-68463 in `Bio.Entrez.Parser`; 1.88 removed an `eval()`-reachable RCE in the `Bio.Nexus` parser.
*Minor caveat:* biopython.org's homepage still advertises 1.87, i.e. the website lags PyPI/GitHub by one release.
*Source:* PyPI JSON (re-queried 2026-08-10); `NEWS.rst:13–52`.

**`documented_examples` = yes** *(new row)* — extensive, and none of it library design.
Three corpora: (1) the **Tutorial & Cookbook**, 27 chapters; (2) its **cookbook chapter**, ~17 recipes (13 sequence-file recipes: filtering, randomized genomes, translating CDS, upper-casing, sorting, FASTQ quality filtering, primer trimming, adaptor trimming, FASTQ/QUAL/SFF conversion, FASTQ indexing, ORF finding; plus 4 plotting recipes: length histogram, GC% plot, dot plots, quality-score plots; plus a BioSQL section); (3) a **separate user-contributed wiki cookbook**, 16 entries (of which "Mapping genetic coordinates with `Bio.SeqUtils.Mapper`" is stale — module 404s — and "From gene sequence to predicted protein with the GFF module" uses the external `BCBio.GFF`). Plus doctests as worked examples, notably the 6-species liftover walkthrough in `chapter_align.html` and the pBAD30 plasmid re-linearization doctest in `Bio/SeqRecord.py`.
**Assessment for the referee response:** not one example produces a designed set of variant sequences. The closest adjacents are "Producing randomized genomes" (a shuffled control), "Methods for Degenerated Codons" (*n*-fold degeneracy analysis), and the pBAD30 two-fragment concatenation. **There is nothing to reproduce head-to-head** — Biopython supplies primitives PoolParty builds on.
*Source:* Tutorial index and `chapter_cookbook.html` (re-fetched 2026-08-10); `biopython.org/wiki/Category:Cookbook`; `Bio/SeqRecord.py`.

---

## 3. What changed from the v1 record

| Key | v1 | v2 | Why |
|---|---|---|---|
| `library_as_object` → `library_first_class_object` | no | **partial** | Split. Biopython genuinely hands the user in-memory `SeqRecord`/`Alignment` objects to inspect, transform and pass onward — it is not a file-writing program. What it lacks is a *library* object with design semantics, which is now the other row. |
| `library_as_object` → `library_algebra` | (no, fused) | **no** | Split. No combine/sample/repeat operation on libraries; alignment `append`/`extend`/`+` stack rows or columns of an alignment, and combining designed sets is plain-Python list concatenation. No sampling or replication method exists anywhere. |
| `dag_chaining` → `composable_operations` | no | **no** | Rename; value and evidence carried. New wording does not change the answer — there are no design steps to compose, and `Bio.Application` was removed in 1.86. |
| `lazy_evaluation` → `lazy_generation` | partial | **partial** | Rename; value carried. Lazy *reading* is first-class and region-level (`SeqIO.parse`/`index`/`index_db`, `twobit` byte-range reads, `Seq(None, length)`, `bgzf`); lazy *generation* is vacuously absent. |
| `mixed_mutagenesis_one_pool` → `exhaustive_single_scans` | no | **no** | Split into four; all four resolve to `no` because Biopython has no mutagenesis operation at all. |
| → `sampled_random_mutagenesis` | (no) | **no** | Re-checked specifically: the one randomization example delegates to Python's `random.shuffle` in user code. |
| → `higher_order_combinatorial` | (no) | **no** | No multi-edit enumeration or combination machinery. |
| → `heterogeneous_components_one_library` | (no) | **no** | Vacuous — no library specification to hold components. |
| `hgvs_input` | no | *dropped* | Per v2 row set. The fact (no HGVS anywhere in `Bio/`) moves to limitations prose. |
| `maintained` | yes | → `last_activity` **yes** | Reframed as a checkable date: 1.88 on 2026-08-06, 1.89 in progress. |
| `degenerate_iupac_codons` | — | **partial** *(new)* | `IUPACData.ambiguous_dna_values`, `AmbiguousCodonTable`/`AmbiguousForwardTable`, `list_possible_proteins`, `list_ambiguous_codons` (codon-set → degenerate codon) are real public machinery — but for interpretation, not degenerate-library design. |
| `negative_control_generation` | — | **no** *(new)* | No shuffle/scramble function in `Bio/`; the cookbook recipe uses Python's `random.shuffle`. `reverse_complement()` exists as a generic strand operation, not as a control feature. |
| `ml_model_in_loop` | — | **no** *(new)* | No model interface; sole dependency is `numpy`. |
| `readout_analysis` | — | **partial** *(new)* | Strong FASTQ/SAM primitives and six read-handling cookbook recipes, but no demultiplexing, variant counting, or enrichment against a designed library. |
| `installable_today` | — | **yes** *(new)* | `pip install biopython`, wheels cp310–cp314 on all mainstream platforms, one dependency. |
| `documented_examples` | — | **yes** *(new)* | Three documented corpora, none of them a library-design example. |

Unchanged: `genome_coordinates` **yes**, `transcript_models` **partial**, `exon_intron_split_codons` **partial**, `vcf_vep_output` **no**, `consequence_annotation` **no**, `combinatorial_motif_place` **no**, `barcode_generation` **no**, `per_sequence_provenance` **partial**, `automatic_naming` **no**, `design_visualization` **partial**, all three `assay_*` **no**, `primer_design` **no**, `codon_optimization` **yes**, `synthesis_constraints` **partial**, `interface` **yes**, `license` **yes**.

---

## 4. Re-checks performed for this pass

Everything else is carried from `final/biopython.md`, which was already extracted, adversarially reviewed, and independently re-verified. For v2 I re-read five things over `raw.githubusercontent.com` / `biopython.org` (nothing installed, nothing executed):

1. `Bio/Data/IUPACData.py` — full top-level name list; `ambiguous_dna_values` definition. → `degenerate_iupac_codons`
2. `Bio/Data/CodonTable.py` — `AmbiguousCodonTable`, `AmbiguousForwardTable`, `ambiguous_dna_by_name`, `list_possible_proteins`, `list_ambiguous_codons` docstring. → `degenerate_iupac_codons`
3. `Bio/Align/__init__.py` — `MultipleSeqAlignment.append`/`extend`/`__add__`/`sort`; confirmed absence of any sampling or replication method. → `library_algebra`
4. Tutorial `chapter_cookbook.html` — full recipe list; the "Producing randomized genomes" recipe's use of Python's `random.shuffle`. → `sampled_random_mutagenesis`, `negative_control_generation`, `readout_analysis`, `documented_examples`
5. PyPI JSON + `pyproject.toml` figures carried forward for `installable_today` / `last_activity`.

One discrepancy against `final/biopython.md` worth recording: the final memo describes "Producing randomized genomes" as *dinucleotide-preserving* shuffling and counts the Tutorial cookbook at 19 items. The re-fetched chapter shows a plain `random.shuffle` (composition-preserving, not dinucleotide-preserving) and ~17 recipe headings plus a BioSQL section. Neither correction changes any cell value; the v2 text above uses the re-verified wording.

---

## 5. Confidence and open judgment calls

Flagged for human review, most consequential first:

1. **`library_first_class_object` = partial** — the deliberate upgrade from v1 `no`. Defensible reading: Biopython hands you real objects, so the row's contrast ("vs a tool that only writes files") puts it above `no`. Stricter reading: the row says *library*, and Biopython has no library object, so `no`. This is the cell most likely to need a matrix-wide consistency decision.
2. **`library_algebra` = no** — borderline. `MultipleSeqAlignment.append`/`extend`/`+` are genuine in-tool set operations; they just operate on alignments, and there is no sample or repeat. A lenient reader could argue `partial`.
3. **`readout_analysis` = partial** — borderline. Biopython's FASTQ/SAM handling is real and documented but has no connection to a designed library; a stricter reading gives `no`.
4. **`degenerate_iupac_codons` = partial** — `list_ambiguous_codons` is closer to degenerate-codon *compression* than I expected, so a lenient reader could push toward `yes`; but nothing designs a degenerate library, so `partial` is the ceiling under the scope rule.
5. **`negative_control_generation` = no** — a lenient reader could give `partial` on `reverse_complement()`/`complement()` alone. I scored `no` because there is no shuffle function at all and no control-generation feature; if the matrix credits generic `reverse_complement` for other tools, this must be re-baselined together.
6. **`lazy_generation` = partial** — carried from v1 by the rename rule, but the new key word is *generation*, and Biopython generates nothing. If the column is read strictly as "lazy production of designed sequences", this is `no` for Biopython (and probably for several other tools) and the whole column should be re-baselined at once.
7. **`composable_operations` = no** — carried. A referee could argue Biopython is maximally composable *as Python* (`rec[:i] + insert + rec[i:]` with feature bookkeeping). The `no` is about tool-level design-step composition, and the cell text should say so.
8. **`codon_optimization` = yes (basic)** — inherited judgment call, still live. The docstring supports `yes`, but the implementation is one line of most-preferred-codon substitution over a user-supplied CAI table with no host presets. **Consistency risk:** if any other tool gets `partial` for a richer optimizer, Biopython at `yes` looks wrong. Ship the docstring quote with the cell or downgrade and re-baseline the column.
9. **`synthesis_constraints` = partial** — inherited; hinges on counting `Bio.Restriction` site scanning as a synthesis check. Stricter reading gives `no`.
10. **`transcript_models` = partial** — inherited; a referee could argue `yes` on BED12/bigBed/PSL exon+CDS models or `no` on the strict GTF/GFF reading. `partial` is defensible **only with the `Bio.Align` evidence attached**.
11. **`per_sequence_provenance` = partial** and **`exon_intron_split_codons` = partial** — inherited; both credit correct machinery that is never exercised for design.
12. **`documented_examples`** — the recipe count (~17 vs the v1 memo's 19) is soft; the load-bearing claim ("none is a library-design example") is verified across all three corpora.
13. **`last_activity`** — the release date and `NEWS.rst` are re-verified; the last-commit SHA/timestamp and the GitHub star/issue counts come from the extraction pass and could not be re-verified (GitHub API rate-limited). Do not quote star/issue numbers.

**Method caveat:** nothing was installed or executed at any stage of v1 or v2. All statements derive from published documentation, the PyPI JSON API, the 2009 PDF (text extracted with PyMuPDF), and master-branch source read over HTTPS.

---

## 6. Bottom line for the referee response

Biopython is infrastructure, not a competitor, and the v2 row set makes that sharper rather than harsher. Under v2 it scores `yes` on five cells (`genome_coordinates`, `codon_optimization`, plus the three availability descriptors `installable_today`, `last_activity`, `documented_examples`; `interface` and `license` are two further `yes` values describing modality and terms rather than design power), `partial` on eight, and `no` on the remaining fifteen — with every Block-A/B `no` tracing to one absence: there is no design layer.

The framing to keep is the one the finalized record fought for: **Biopython does chromosome-name resolution, region-lazy genome access, indexed region query, and UCSC chain-file liftover, and it does them well — it just has nothing to design.** Under v2 that is joined by two more concessions worth making explicitly, because they are real and a Biopython author would insist on them: Biopython hands the user first-class, transformable sequence objects (`library_first_class_object` = partial), and it ships genuine IUPAC/degenerate-codon tables and helpers (`degenerate_iupac_codons` = partial). Neither amounts to library design, and saying so plainly is stronger than pretending they are absent.
