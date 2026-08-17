# Adversarial review — Biopython extraction

**Reviewer pass:** 2026-08-10. Independent re-verification against master-branch source (`raw.githubusercontent.com`), the live 1.88 Tutorial (`biopython.org/docs/latest`), PyPI JSON, and the 2009 PDF. Nothing installed or executed.

**Headline:** the extraction is careful, well-sourced, and correct on all 18 Block-A/B/D/E cells. It is **materially wrong on one cell (`genome_coordinates`) and incomplete on a second (`transcript_models`)**, both in the same direction — it understates modern `Bio.Align`, which is exactly the part of Biopython that a referee named **Michiel de Hoon** wrote. Two of the three most-quoted "no such thing exists" phrases in the `genome_coordinates` evidence are falsifiable in one line of documented code.

---

## 1. The falsifications

### 1.1 `genome_coordinates = partial` — evidence is demonstrably wrong

The extraction's evidence says, verbatim:

> "there is no assembly/genome-build concept, no chromosome-name resolution against a reference, no coordinate-system conversion or liftover: 'give me hg38 chr7:1000-2000' requires the user to load the FASTA and slice it themselves."

Every clause after the colon is false as of 1.88 (and has been since ~1.80–1.83).

**(a) Liftover between genome assemblies is a documented, first-class feature.**
`Bio/Align/__init__.py:3339 Alignment.map()` docstring:

> "The map method can also be used to lift over an alignment between different genome assemblies. In this case, self is a DNA alignment between two genome assemblies, and the argument is an alignment of a transcript against one of the genome assemblies"

The Tutorial chapter *Sequence alignments* (`chapter_align.html`) has a worked example that reads UCSC `.chain` files and lifts a 6-species MAF alignment onto new builds:

> "These are provided by UCSC as chain files, typically used for UCSC's liftOver tool. … `paths = ["Blat/panTro5ToPanTro6.chain", "Blat/hg19ToHg38.chain", "Blat/rheMac8ToRheMac10.chain", "Blat/calJac3ToCalJac4.chain", "Blat/mm10ToMm39.chain", "Blat/rn6ToRn7.chain"]` … `genome_alignment = genome_alignment.mapall(liftover_alignments)`"

`Bio/Align/chain.py` is a full read/write implementation of the UCSC chain format ("Bio.Align support for the 'chain' pairwise alignment format … See https://genome.ucsc.edu/goldenPath/help/chain.html").

**(b) Indexed query by chromosome region exists.**
`Bio/Align/bigbed.py:1042`:

```python
def search(self, chromosome=None, start=None, end=None):
    """Iterate over alignments overlapping the specified chromosome region.
    This method searches the index to find alignments to the specified
    chromosome that fully or partially overlap the chromosome region
    between start and end."""
```

Tutorial: `selected_alignments = alignments.search("chr3", 48000000, 49000000)` — with the same method available on bigPsl and bigMaf ("we can look for regions on chr10 between positions 3018000 and 3019000 … `alignments.search("mm9.chr10", 3018000, 3019000)`"). This is an R-tree index lookup, not a linear scan.

**(c) Chromosome-name resolution against a genome, lazily, is one line.**
`Bio/SeqIO/TwoBitIO.py` (UCSC `.2bit`, i.e. the format UCSC actually distributes hg38 in) is a **registered SeqIO format** (`_FormatToIterator["twobit"]`, `Bio/SeqIO/__init__.py:450`) whose parser is dict-like and region-lazy:

> "Using the information in the index, the `__getitem__` method calculates the file position at which the requested region starts, and only reads the requested sequence region. Note that the full sequence of a record is loaded only if specifically requested, making the parser memory-efficient. The TwoBitIterator object implements the `__getitem__`, keys, and `__len__` methods that allow it to be used as a dictionary."

Tutorial (`chapter_align.html`) uses exactly the hg38 case: `genome = SeqIO.parse("hg38.2bit", "twobit"); genome[chromosome]`. So "give me hg38 chr7:1000-2000" is `SeqIO.parse("hg38.2bit","twobit")["chr7"][1000:2000]`, without loading chr7.

**(d) Assembly identity is carried in the data model.** MAF/bigMaf records are keyed `hg19.chr1`, `mm10.chr3`, etc.; `Seq(None, length=248956422)` (undefined sequences, `.defined`, `.defined_ranges`) lets a chromosome be represented by name+length alone for pure coordinate work — the Tutorial demonstrates `SeqRecord(Seq(None, length=1575), id="chr1")` precisely to supply a target size for PSL output.

**Verdict: understated.** The cell value should be `yes` for coordinate handling (with the honest rider that there is no *design* step to attach coordinates to). At minimum the evidence string must be rewritten — as written it is the single most attackable sentence in the extraction.

### 1.2 `transcript_models = partial` — right value, wrong reasoning

The "GTF/GFF is not in core" fact is correct and I re-verified it (`Bio/` has no GFF module; the wiki still says GFF parsing "is not yet integrated"). But the evidence then reasons *from* the `Bio.SeqIO` format table to a conclusion about transcript annotation, and that inference does not hold:

- `Bio/Align/bed.py` is a **read and write** implementation of BED, described in its own docstring as: *"The Browser Extensible Data (BED) format … Typically they are used for **transcript to genome alignments**."* It parses `blockCount`/`blockSizes`/`blockStarts` (exons) into `Alignment.coordinates`, and exposes `alignment.thickStart` / `alignment.thickEnd` — which the UCSC spec the Tutorial quotes defines as *"Start of where display should be thick (start codon)"* / *"(stop codon)"*, i.e. the CDS boundaries. That is a transcript model with exons and a CDS.
- The same is true of `psl.py` / `bigpsl.py` (blockSizes/qStarts/tStarts), `bigbed.py`, and `Bio/SearchIO/ExonerateIO`.
- Separately, the quoted SeqIO format table is itself incomplete: it omits **`twobit`**, and it omits the AlignIO formats that `Bio/SeqIO/__init__.py:509–519` auto-registers into `SeqIO`.

**Verdict: understated.** Keep `partial` if you like, but the evidence must say "GTF/GFF absent from core (use BCBio.GFF); BED12/bigBed/PSL/bigPsl exon-block transcript models are read *and written* natively via `Bio.Align`."

---

## 2. Cells I tried and failed to break (all confirmed)

| Cell | What I checked independently | Result |
|---|---|---|
| `library_as_object` = no | Full `Bio/` dir listing, `pyproject.toml` package list, `Bio/Align/__init__.py` class list | Confirmed. (Extractor missed the newer `Bio.Align.Alignments` collection class — still not a design object.) |
| `dag_chaining` = no | `DEPRECATED.rst:228–232` quoted verbatim, correct; no pipeline module anywhere | Confirmed |
| `mixed_mutagenesis_one_pool` = no | `Bio/Seq.py:2173 MutableSeq` docstring quoted verbatim, correct; no mutagenesis anywhere | Confirmed |
| `combinatorial_motif_place` = no | `Bio/motifs/__init__.py` full method list: `pwm, pssm, consensus, anticonsensus, degenerate_consensus, relative_entropy, reverse_complement, weblogo, format` + parsers only. No sampler, no placer. | Confirmed |
| `barcode_generation` = no | `Bio/SeqUtils/` = `CheckSum, IsoelectricPoint, MeltingTemp, ProtParam, ProtParamData, lcc, __init__` | Confirmed |
| `per_sequence_provenance` = partial | `Bio/SeqRecord.py __getitem__` docstring: *"the annotations dictionary and the dbxrefs list are not used for the new SeqRecord"* | Confirmed verbatim |
| `automatic_naming` = no | — | Confirmed |
| `design_visualization` = partial | GenomeDiagram real; **but** `Bio/Graphics/__init__.py` raises `MissingPythonDependencyError` without ReportLab, and PyPI `requires_dist == ["numpy"]` — so `pip install biopython` cannot draw anything | Confirmed, with an added caveat that *weakens* the cell |
| `assay_dms` / `assay_mpra` / `assay_insilico` = no | Full Tutorial + Cookbook TOC re-fetched; also checked the *second*, user-contributed wiki cookbook (16 entries) the extractor did not enumerate — nothing relevant | Confirmed |
| `exon_intron_split_codons` = partial | `CompoundLocation` + `extract` + trans-splicing since 1.78 | Confirmed |
| `hgvs_input` = no | No HGVS anywhere in `Bio/` | Confirmed |
| `vcf_vep_output` = no | No VCF module, not a SeqIO/AlignIO format | Confirmed (though "writes sequence and alignment formats only" now also covers BED/PSL/SAM/bigBed genomic-interval output — phrase it as "no variant-call formats") |
| `consequence_annotation` = no | Nearest thing is `Bio/Align/analysis.py calculate_dn_ds` (NG86/LWL85/YN00), which classifies synonymous vs non-synonymous *sites* in a codon alignment — not per-variant consequence calls | Confirmed |
| `primer_design` = no | `Bio/Emboss/` = `Primer3.py, PrimerSearch.py, __init__.py` only; `Primer3.py` docstring is verbatim *"Code to parse output from the EMBOSS eprimer3 program."*; `Bio.Application` removal confirmed in `DEPRECATED.rst` | Confirmed |
| `codon_optimization` = yes | `Bio/SeqUtils/__init__.py:667 optimize()` read in full — quote is verbatim and complete. Also confirms the caveat: `__init__(sequences, ...)` requires user CDS, and `CodonUsageIndices.py` / `SharpEcoliIndex` no longer ship, so "no host tables" is right. Implementation is literally `"".join(pref_codons[aa] for aa in aa_seq)`. | Confirmed. `partial` is defensible; keep the docstring quote attached to the cell. |
| `synthesis_constraints` = partial | `Bio/Restriction/__init__.py` doctest quoted verbatim, correct | Confirmed |
| `interface` = yes | `pyproject.toml` re-read: no `[project.scripts]`, no console entry points; deps `["numpy"]` | Confirmed. Flag: the bare value `yes` doesn't encode *which* interface — a reader could take it as "has a CLI". Suggest the value carry the modality. |
| `license` = yes | `pyproject.toml`: `license = "LicenseRef-Biopython-License-Agreement"`, `license-files = ["LICENSE.rst"]` | Confirmed |
| `maintained` = yes | PyPI: version **1.88**, first file uploaded **2026-08-06T12:30:38Z**, 31 files, `requires-python >=3.10`, `requires_dist ["numpy"]`. Wheel list confirms cp310–cp314 incl. **cp314t**, manylinux x86_64 + aarch64, macosx_11_0_arm64, win32 + win_amd64. `NEWS.rst` confirms the 1.88 date, the `Bio.Nexus` `eval()` RCE fix, and the "In progress: 1.89" section. | Confirmed. (GitHub API was rate-limited at review time, so the 5,158 stars / 607 issues / commit-SHA figures are unverified — low risk, but they are the only unchecked numbers.) |
| Paper quotes | All four PDF quotes re-extracted with PyMuPDF and matched verbatim, incl. *"a large open-source application programming interface (API) used in both bioinformatics software development and in everyday scripts"* (CONCLUSIONS) and *"The features described herein are only a subset"* | Confirmed |

---

## 3. Capabilities missed entirely

1. **Modern `Bio.Align` genomic file I/O** — read *and write* for `bed, bigbed, psl, bigpsl, maf, bigmaf, chain, sam, exonerate, a2m, mauve, hhr, msf, tabular`. The extraction's `additional_capabilities` mentions "Bio.AlignIO (clustal, phylip, stockholm, nexus, maf, ...)" and stops there, which reads like 2010-era Biopython. This is the biggest omission and it is the newest, most actively developed part of the package.
2. **Assembly liftover** via UCSC chain files (`Alignment.map` / `mapall`), documented with hg19→hg38, mm10→mm39, panTro5→panTro6.
3. **Indexed genomic-region query** `alignments.search(chrom, start, end)` over bigBed/bigPsl/bigMaf.
4. **`twobit` lazy genome reader** in `Bio.SeqIO` — chromosome-keyed dict access with byte-range reads. Missing from the extraction's SeqIO format list entirely, and it is the single strongest piece of evidence for both `lazy_evaluation` and `genome_coordinates`.
5. **Undefined / partially-defined sequences** — `Seq(None, length)`, `.defined`, `.defined_ranges`; lets records carry chromosome-scale coordinates with no sequence in memory.
6. **`SeqRecord` concatenation with feature bookkeeping** — `record_a + record_b` merges annotations and shifts feature locations (documented pBAD30 plasmid re-linearization doctest). The extraction twice says insertion "would be hand-written string surgery on MutableSeq"; that undersells the real construct-assembly primitive and is the kind of sentence a Biopython author would push back on. It does not change any `no`.
7. **`Bio.Align.CodonAligner` + `Bio.Align.analysis.calculate_dn_ds`** (NG86, LWL85, YN00, ML) — codon-aware alignment now lives in `Bio.Align`, not only the legacy `Bio.codonalign`.
8. **`Bio.Graphics` requires ReportLab**, which `pip install biopython` does not install. Belongs in `availability_status` and as a caveat on `design_visualization`.
9. **`Bio.motifs.weblogo()` POSTs to the remote WebLogo server** — a network dependency inside an otherwise offline library.
10. **Security posture**: 1.87 fixed **CVE-2025-68463** in `Bio.Entrez.Parser`; 1.88 fixed the `Bio.Nexus` `eval()` RCE. The extraction notes only the second. Two consecutive security releases is worth one clause if the table has a maintenance column.
11. **Second cookbook** — `biopython.org/wiki/Category:Cookbook`, 16 user-contributed entries, is a separate corpus from the Tutorial cookbook. I enumerated it: nothing is library design (closest are "Methods for Degenerated Codons" = *n*-fold degeneracy analysis, and "From gene sequence to predicted protein with the GFF module" = external BCBio.GFF). Worth one sentence so the "I checked the full cookbook" claim is airtight.
12. `Bio.SearchIO` gained **Infernal** parsing in 1.86 — trivial, listed only for completeness.

---

## 4. Wording risks (no verdict change)

- `interface = "yes"`: uninformative for a key whose interesting content is "library only, no CLI, no web service". Encode the modality in the value.
- `codon_optimization = "yes"`: I confirmed the docstring verbatim, so "yes" is quotable. But the implementation is a one-line most-preferred-codon substitution over a user-supplied CAI table. If any other tool in the matrix gets `partial` for a richer optimizer, Biopython at `yes` will look inconsistent. Ship the quote with the cell, or downgrade to `partial (naive)`.
- The repeated construction "*X does not exist because Biopython has no design step*" is correct but appears in five cells. A referee will read it as a template. Vary it, and in the two coordinate cells replace it with what Biopython *does* do.

## 5. Bottom line

Do not ship the `genome_coordinates` evidence as written — it asserts three absences that a five-line doctest in Biopython's own Tutorial refutes, in a module written by one of the paper's authors. Fix that cell and the `transcript_models` reasoning, add `twobit` and the `Bio.Align` UCSC formats to the additional-capabilities list, and the extraction is defensible everywhere else.
