# Biopython — evidence memo for the PoolParty tool-comparison table

**Surveyed:** 2026-08-10
**Tool slug:** `biopython`
**Category:** general-purpose bioinformatics library (sequence I/O + manipulation), *not* a library-design tool.

---

## 1. Sources consulted

| Kind | Reference |
|---|---|
| PDF | `papers/Cock2009_biopython.pdf` — Cock et al. (2009) *Bioinformatics* 25(11):1422–1423, doi:10.1093/bioinformatics/btp163 |
| prior analysis | `prior_analyses/Cock2009_biopython_analysis.md` |
| repo | https://github.com/biopython/biopython (master; GitHub API queried 2026-08-10) |
| repo (raw source read) | `Bio/SeqUtils/__init__.py`, `Bio/Seq.py`, `Bio/Emboss/Primer3.py`, `Bio/Restriction/__init__.py`, `LICENSE.rst`, `pyproject.toml`, `NEWS.rst`, `DEPRECATED.rst` |
| docs | https://biopython.org/docs/latest/Tutorial/index.html (Tutorial & Cookbook, 1.88) |
| docs | https://biopython.org/docs/latest/api/Bio.html (API index), `Bio.SeqIO`, `Bio.SeqUtils`, `Bio.motifs` API pages |
| docs | https://biopython.org/docs/latest/Tutorial/chapter_cookbook.html, `chapter_seq_annot.html`, `chapter_graphics.html` |
| wiki | http://biopython.org/wiki/GFF_Parsing |
| pypi | https://pypi.org/pypi/biopython/json |
| web | https://biopython.org/ (homepage) |

**Nothing was installed.** Source files were read via `raw.githubusercontent.com`; metadata via the GitHub and PyPI JSON APIs.

---

## 2. What Biopython actually is

The 2009 applications note is 2 pages and describes the library as of 1.49-ish. Its own conclusion says: *"The features described herein are only a subset; potential users should refer to the tutorial and API documentation for further information."* So the assessment below is anchored on the **current 1.88 docs and source**, per the task hint.

Core abstractions (unchanged in kind since 2009):

- `Bio.Seq.Seq` — immutable sequence, string-like, with `translate()`, `reverse_complement()`, `transcribe()`.
- `Bio.Seq.MutableSeq` — *"An editable sequence object. … the MutableSeq lets you edit the sequence in place"* (`Bio/Seq.py:2174`). This is per-position assignment (`my_seq[5] = "A"`), i.e. the raw primitive a user would loop over to make variants themselves.
- `Bio.SeqRecord.SeqRecord` — sequence + `id`/`name`/`description` + free-form `annotations` dict + `features` list + `letter_annotations`.
- `Bio.SeqFeature` — `SimpleLocation` / `CompoundLocation` (the latter *"handles 'join' locations in EMBL/GenBank files"*, since 1.62) + `.extract()`.
- `Bio.SeqIO` / `Bio.AlignIO` / `Bio.Align` — I/O and alignment.

Full current subpackage list (from the `Bio` API index): Affy, Align, AlignIO, Blast, CAPS, Cluster, Compass, Data, Emboss, Entrez, ExPASy, GenBank, Geo, Graphics, HMM, KEGG, Medline, NMR, Nexus, PDB, Pathway, Phylo, PopGen, Restriction, SCOP, SVDSuperimposer, SearchIO, SeqIO, SeqUtils, Sequencing, SwissProt, TogoWS, UniGene, UniProt, codonalign, motifs, phenotype; submodules File, Seq, SeqFeature, SeqRecord, bgzf, pairwise2.

**There is no module, class, function, tutorial chapter, or cookbook recipe anywhere in Biopython concerned with designing a variant library, a mutagenesis scheme, an oligo pool, or a combinatorial construct.** The prior analysis's core claim is correct and I confirmed it against the 1.88 API index and the full Tutorial + Cookbook table of contents.

---

## 3. Capability-by-capability evidence

### Block A — library specification

**`library_as_object` = no.**
The only multi-sequence containers are (a) a plain Python `list`/generator of `SeqRecord`s, (b) `Bio.Align.Alignment` / `MultipleSeqAlignment` (an *alignment*, i.e. equal-length gapped rows — Bio.SeqIO "interprets multiple sequence alignment file formats as collections of equal length (gapped) sequences"), and (c) the dict-like `SeqIO.to_dict()` / `SeqIO.index()` / `SeqIO.index_db()` views over a file. None of these carries design semantics — no notion of a wild-type reference, a variant set, or an operation that produced the members. Building a library means writing your own loop over `MutableSeq`.

**`dag_chaining` = no.**
No pipeline/graph/composition abstraction exists in the API index. Chaining in Biopython is ordinary Python function composition written by the user. (`Bio.Application` — the only thing that ever resembled a pipeline layer, for chaining external command-line tools — *"Declared obsolete in release 1.79, deprecated in release 1.82, and removed in release 1.86"* per `DEPRECATED.rst`.)

**`lazy_evaluation` = partial.**
Lazy *parsing* is real and first-class: `SeqIO.parse()` returns an iterator; `SeqIO.index()` gives a dict-like object where *"When you access a particular record via the dictionary methods, the code will jump to the appropriate part of the file and then parse that section into a SeqRecord"*; `index_db()` stores that index in SQLite. Cookbook has a dedicated recipe "Indexing a FASTQ file". **But** this is lazy reading of sequences that already exist on disk — there is no lazy *generation* of designed sequences, because there is no design step.

**`mixed_mutagenesis_one_pool` = no.**
Biopython has no mutagenesis operation at all (exhaustive single, double, random, or WT-replicate). Nothing in the API index or Tutorial/Cookbook TOC. The closest primitive is per-base assignment on `MutableSeq`.

**`combinatorial_motif_place` = no.**
`Bio.motifs` is analysis-only. Its documented surface is `consensus`, `anticonsensus`, `degenerate_consensus`, `pwm`, `pssm`, `relative_entropy`, `reverse_complement()`, `weblogo()`, plus parsers (MEME, TRANSFAC, JASPAR, …) and PSSM searching. There is **no** function to sample a sequence from a PWM or to insert/place a motif into a host sequence, and no combinatorial placement machinery. (Insertion would be hand-written string surgery on a `MutableSeq`.)

**`barcode_generation` = no.**
No barcode/UMI/index generator anywhere. Constraint primitives a user could build one from exist — `Bio.SeqUtils.gc_fraction()`, `Bio.SeqUtils.MeltingTemp`, `Bio.Align.PairwiseAligner` (for distance), `Bio.SeqUtils.lcc` (complexity) — but there is no generator, no edit-distance-constrained set construction, and no attach-to-library step.

**`per_sequence_provenance` = partial.**
`SeqRecord.annotations` is *"A dictionary of additional information about the sequence. The keys are the name of the information, and the information is contained in the value"* — a free-form container, plus a structured `features` list of `SeqFeature`s with locations and qualifiers, and `letter_annotations` for per-position data; these round-trip to GenBank. So a structured metadata slot exists and is idiomatic. **However** Biopython never constructs a designed sequence, so nothing is ever recorded automatically: any provenance is 100% user-written. Note also that slicing a `SeqRecord` deliberately *drops* the `annotations` dict (it keeps id/name/description, in-range features, and letter_annotations), so even simple edits do not propagate metadata.

**`automatic_naming` = no.**
`SeqRecord.id` / `.name` / `.description` are user-supplied or read from a file. No variant-naming scheme (no HGVS-style, position-based, or operation-derived name generation) exists.

**`design_visualization` = partial.**
`Bio.Graphics.GenomeDiagram` is a genuine, well-documented visualisation layer: *"designed for drawing whole genomes, in particular prokaryotic genomes, either as linear diagrams (optionally broken up into fragments to fit better) or as circular wheel diagrams"*, with stacked tracks, feature sets, sigils (BOX/ARROW/OCTO/JAGGY/BIGARROW), cross-links between tracks, and PDF/EPS/SVG/PNG output via ReportLab; plus `BasicChromosome` and `motifs.weblogo()`. Cookbook adds "Nucleotide dot plots", "Histogram of sequence lengths", "Plot of sequence GC%". So *annotated sequences* can be drawn. What cannot be inspected is a **design**: there is no library/pool view, no operation-graph view, and no highlighting of designed variants relative to a reference, because those objects do not exist.

### Block B — assay coverage

**`assay_dms` = no.** No DMS/saturation-mutagenesis functionality; no such term in the Tutorial/Cookbook TOC.
**`assay_mpra` = no.** No regulatory/MPRA library design; `Bio.motifs` is TFBS *analysis*, not reporter-library construction.
**`assay_insilico` = no.** No model-probing, no in-silico perturbation design. Biopython's ML modules are the legacy `Bio.HMM`, `Bio.Cluster`, `Bio.PopGen` — classical statistics, nothing to do with genomic deep-learning models.

### Block C — genomics integration

**`genome_coordinates` = partial.**
Coordinates are a first-class idea at the *record* level: `SimpleLocation`/`CompoundLocation` carry start/end/strand, `SeqFeature.extract()` pulls the subsequence, `SeqRecord` supports slicing, and `Bio.Entrez.efetch` accepts `seq_start`/`seq_stop` to retrieve a genomic region remotely. **But** there is no assembly/build concept, no chromosome-name resolution service, no coordinate-system conversion (0- vs 1-based, lift-over), and no "give me hg38 chr7:1000-2000" entry point against a local reference beyond loading the FASTA yourself and slicing.

**`transcript_models` = partial (GTF/GFF specifically: no).**
GFF/GTF is **not** supported by `Bio.SeqIO` — the format list is abi, ace, cif-atom, cif-seqres, embl, fasta, fasta-2line, fastq(+variants), gck, genbank/gb, gfa1, gfa2, ig, imgt, nib, pdb-seqres, pdb-atom, phd, pir, seqxml, sff, snapgene, swiss, tab, qual, uniprot-xml, xdna. The project's own wiki is explicit: *"GFF parsing is not yet integrated into Biopython. This documentation is work towards making it ready for inclusion"*, and directs users to the separate `BCBio.GFF` (bcbio-gff) package / gffutils. What Biopython *does* understand natively is the GenBank/EMBL feature table: `mRNA`/`CDS` features with `join(...)` locations are parsed into `CompoundLocation`, which is a transcript model in all but file format.

**`exon_intron_split_codons` = partial.**
The mechanism is correct and documented: `CompoundLocation` represents joined exons, and *"The `SeqFeature` object has an `extract` method to take care of all this (and since Biopython 1.78 can handle trans-splicing by supplying a dictionary of referenced sequences)"* — i.e. exons are concatenated in transcript order with strand handled, so a codon split across an exon boundary translates correctly via `feature.extract(genome).translate()`. Cookbook recipe "Translating a FASTA file of CDS entries" covers the CDS case. **But** because Biopython performs no mutagenesis, it never has to *place* an edit into a codon that straddles an intron — the hard part of this capability in a design tool is simply not exercised.

**`hgvs_input` = no.**
No HGVS module, class, or parser in the `Bio` API index; no HGVS section in the Tutorial. HGVS in the Python ecosystem is the separate biocommons `hgvs` package, not Biopython.

**`vcf_vep_output` = no.**
VCF is absent from the `Bio.SeqIO` format table and there is no `Bio.VCF`/variant module. Biopython writes sequence/alignment formats (FASTA, GenBank, EMBL, Clustal, PHYLIP, Stockholm, NEXUS, …), not variant-call formats. No VEP interoperability.

**`consequence_annotation` = no.**
No molecular-consequence classifier. `Seq.translate()` + `Bio.Data.CodonTable` let a user compare two translations by hand, but there is no stop-gained / synonymous / frameshift / in-frame-indel annotation function. (`Bio.CAPS` is about restriction-site polymorphisms in an alignment — a different thing.)

### Block D — physical construction

**`primer_design` = no.**
`Bio.Emboss.Primer3` is explicitly a **parser**, not a designer: its module docstring is *"Code to parse output from the EMBOSS eprimer3 program."* `Bio.Emboss.PrimerSearch` similarly writes input for / parses output from EMBOSS primersearch. The command-line wrappers that used to *launch* eprimer3 (`Bio.Emboss.Applications`, built on `Bio.Application`) were **removed in release 1.86**. `Bio.SeqUtils.MeltingTemp` computes Tm and Cookbook has "Trimming off primer sequences" (trimming, not design). Net: Biopython can read someone else's primer designs and score a Tm; it designs no primers or mutagenic oligos itself.

**`codon_optimization` = yes (basic).**
`Bio.SeqUtils.CodonAdaptationIndex.optimize()` exists and is documented: *"Return a new DNA sequence with preferred codons only. Uses the codon adaptiveness table defined by the CodonAdaptationIndex object to generate DNA sequences with only preferred codons. May be useful when designing DNA sequences for transgenic protein expression or codon-optimized proteins like fluorophores."* It accepts DNA/RNA/protein input and a `strict` flag for ties. The CAI table itself is learned from a user-supplied set of coding sequences (`CodonAdaptationIndex.__init__(sequences, table=standard_dna_table)`), and `calculate()` scores CAI. **Caveat for the table:** this is the naive "replace every codon with the single most-preferred codon" strategy — no codon harmonisation, no usage-frequency sampling, no GC-window / repeat / restriction-site / secondary-structure constraints, no host presets shipped.

**`synthesis_constraints` = partial.**
`Bio.Restriction` is a first-class restriction-site analyser (`Analysis(AllEnzymes, seq)`, `.print_that()`, `.blunt()`, and enzyme-set filtering over ~1000 REBASE enzymes), which covers the single most common synthesis/cloning constraint check (unwanted restriction sites). Add `gc_fraction()`, `MeltingTemp`, `molecular_weight()`, `lcc` (local composition complexity), `nt_search()`. **But** there is no synthesis-constraint feature as such: no homopolymer/repeat rules, no GC-window enforcement, no oligo length/vendor limits, no "reject or repair this sequence" API — the user assembles checks from primitives, and nothing enforces anything on a library.

### Block E — engineering

**`interface`** — Python library API only. `pyproject.toml` declares no `[project.scripts]` / console entry points, so `pip install biopython` installs **no** command-line tool. The repo's `Scripts/` directory (GenBank, PDB, Performance, Restriction, SeqGui, query_pubmed.py, scop_pdb.py, update_ncbi_codon_table.py, xbbtools) holds contributed demo scripts — including two Tk GUI demos, `SeqGui` and `xbbtools` — that are **not installed by the package**. No web service.

**`license`** — Biopython License Agreement (a permissive MIT/BSD-style grant). `pyproject.toml`: `license = "LicenseRef-Biopython-License-Agreement"`. `LICENSE.rst`: *"Biopython is currently released under the 'Biopython License Agreement' … Some files are explicitly dual licensed under your choice of the 'Biopython License Agreement' or the 'BSD 3-Clause License' … with the intention of later offering all of Biopython under this dual licensing approach."* (GitHub's licence detector reports NOASSERTION because of the custom text.)

**`maintained`** — extremely active. Latest release **1.88 on 2026-08-06** (PyPI upload 2026-08-06T12:30:38; `NEWS.rst`: *"6 August 2026: Biopython 1.88"*). Last commit to master 2026-08-06T14:20:26Z ("Post Biopython 1.88 release version bump"); repo `pushed_at` 2026-08-06. 1.89 already "in progress" in NEWS. Recent release cadence: 1.83 (2024-01-10), 1.84 (2024-06-28), 1.85 (2025-01-15), 1.86 (2025-10-28), 1.87 (2026-03-30), 1.88 (2026-08-06). 5,158 GitHub stars, 607 open issues, not archived. Supports Python 3.10–3.14 (+3.15rc), only dependency is numpy. Note: biopython.org's homepage still advertises 1.87 (30 March 2026) — the website lags PyPI/GitHub by one release.

---

## 4. Biopython's own documented examples (candidates for PoolParty reproduction)

Biopython ships a Tutorial & Cookbook. Chapters: Introduction; Quick Start; Sequence objects; Sequence annotation objects; Sequence Input/Output; Sequence alignments; Pairwise sequence alignment; Multiple Sequence Alignment objects; Pairwise alignments using pairwise2; BLAST (new); BLAST (old); BLAST and other sequence search tools; Accessing NCBI's Entrez databases; Swiss-Prot and ExPASy; Going 3D: The PDB module; Bio.PopGen; Phylogenetics with Bio.Phylo; Sequence motif analysis using Bio.motifs; Cluster analysis; Graphics including GenomeDiagram; KEGG; Bio.phenotype; Cookbook; testing framework; contributing; Python appendix; Bibliography.

Cookbook recipes (the closest thing to "example libraries"):

*Working with sequence files* — Filtering a sequence file; Producing randomized genomes; Translating a FASTA file of CDS entries; Making the sequences in a FASTA file upper case; Sorting a sequence file; Simple quality filtering for FASTQ files; Trimming off primer sequences; Trimming off adaptor sequences; Converting FASTQ files; Converting FASTA and QUAL files into FASTQ files; Indexing a FASTQ file; Converting SFF files; Identifying open reading frames.

*Sequence parsing plus simple plots* — Histogram of sequence lengths; Plot of sequence GC%; Nucleotide dot plots; Plotting the quality scores of sequencing read data.

*BioSQL* — storing sequences in a relational database.

**Assessment for the referee response:** none of these is a *library design* example. Not one recipe produces a designed set of variant sequences. The only ones remotely adjacent to PoolParty's domain are "Producing randomized genomes" (dinucleotide-preserving shuffling of an existing genome, an analysis control — not a designed library), "Translating a FASTA file of CDS entries", and "Trimming off primer sequences". There is nothing to reproduce head-to-head; the honest framing is that Biopython supplies primitives PoolParty builds on, not competing designs.

---

## 5. Notable capabilities outside the row list

- **Sequence file I/O breadth** (~28 formats read/write incl. GenBank, EMBL, FASTQ variants, SnapGene, GCK, xDNA, GFA1/2) — the practical reason every design tool depends on it; PoolParty's own export story is judged against this.
- **Restriction-enzyme analysis** over the full REBASE enzyme set (`Bio.Restriction.Analysis`), with blunt/sticky/isoschizomer filtering.
- **Alignment engine**: `Bio.Align.PairwiseAligner` (C-accelerated global/local/affine, substitution matrices) plus `Bio.AlignIO`/`Bio.Align` MSA objects and `Bio.SearchIO`/`Bio.Blast` parsers.
- **Remote database access**: `Bio.Entrez` (E-utilities, incl. coordinate-limited `efetch`), `Bio.ExPASy`, `Bio.KEGG`, `Bio.UniProt`, `Bio.TogoWS`.
- **3D structure**: `Bio.PDB` (PDB/mmCIF parsing, superposition, DSSP/NACCESS interfaces) — irrelevant to library design but a major part of the package's footprint.
- **Phylogenetics** (`Bio.Phylo`, `Bio.Nexus`), **population genetics** (`Bio.PopGen`), **clustering** (`Bio.Cluster`), **protein physicochemistry** (`Bio.SeqUtils.ProtParam`, `IsoelectricPoint`), **checksums** (`Bio.SeqUtils.CheckSum`), **codon alignments / dN-dS** (`Bio.codonalign`), **bgzip random access** (`Bio.bgzf`), **BioSQL ORM**.
- **Typed**: ships `py.typed`; type annotations are being expanded (1.88 news).

None of these argues for a new row in the comparison matrix except possibly a "sequence file format breadth / interoperability" row, where Biopython would be the ceiling.

---

## 6. Bottom line and confidence

Biopython is infrastructure, not a competitor. Of 24 capability keys it scores `yes` on exactly three (codon_optimization, plus the two engineering descriptors), `partial` on seven, and `no` on the rest. It is maximally available and maintained (release 4 days before this survey), which makes it a safe, uncontroversial anchor in the "general-purpose toolkit" row of the matrix.

**Least-confident cells** (state them carefully in the paper):
1. `codon_optimization` = yes. Defensible from the docstring, but a referee could equally argue "this is a one-line most-frequent-codon substitution, not codon optimisation". The evidence quote should travel with the cell.
2. `synthesis_constraints` = partial. Hinges on counting `Bio.Restriction` site-scanning as a synthesis check. A stricter reading gives `no`.
3. `transcript_models` = partial. `no` under a literal GTF/GFF reading (and the wiki backs that: GFF is not in core); `partial` credits GenBank/EMBL `join()` CompoundLocations.
4. `per_sequence_provenance` = partial. Credits the existence of `SeqRecord.annotations`/`features` as a metadata container; nothing is populated automatically, and `SeqRecord` slicing discards `annotations`.
5. `exon_intron_split_codons` = partial. The splicing machinery is genuinely correct; it is just never used for variant placement.
6. `lazy_evaluation` = partial. Lazy *parsing* is real (`SeqIO.parse`/`index`/`index_db`); lazy *generation* is vacuously absent.

I did **not** verify by execution (no installs permitted); all statements come from published docs and the master-branch source read over HTTPS.
