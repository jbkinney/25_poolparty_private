# Audit: `genome_coordinates`

> ## CORRECTION applied after scoring — PoolParty demoted
>
> **`poolparty`: `yes` → `partial`.** Verified behaviourally, not by re-reading docs.
>
> `from_fasta` in batch mode does emit a coordinate-bearing name (`chr1:10-30(+)`),
> but nothing composes it into a per-variant position:
>
> * Four insertion variants at different offsets all carried the **identical** name
>   `chr1:10-30(+)`; output columns after a design op are `['name','seq']` only.
> * `from_fasta`'s valid card keys are `['seq','seq_index','seq_name','state']` —
>   requesting `chrom`/`start`/`stop`/`strand` is rejected. The origin exists only as
>   a packed string requiring regex parsing.
> * Offsets are emitted, but operation-namespaced
>   (`op[2]:insertion_scan(region_scan).start`) and **identical across strands**:
>   `chr1:10-30(+)` and `chr1:10-30(-)` both give `0,1,2`, yet genomic position is
>   `start+offset` on plus and `stop-1-offset` on minus. The tool gives no warning;
>   composing them wrongly fails silently.
> * With a single coordinate tuple rather than a list, the name is `None` entirely.
>
> This is the row's `partial` branch verbatim: *coordinates accepted for extraction
> only, with no coordinate output for designed variants.* Consistency check: VaLiAnT
> and MPRAnator keep `yes` because both emit a per-variant position that moves with
> the variant (`mut_position`; absolute `Position` in the header). MPRA Design Tools
> keeps `partial` because its final output drops the coordinate fields.


Audit date: 2026-08-15

## Measurement instrument and operational test

I read the preamble and scoring vocabulary in
`25_poolparty_private/revision/tables/ROW_DEFINITIONS.md:1-12`, the binding
global rule at lines 14-34, and the `genome_coordinates` definition at lines
174-177 before inspecting any tool. The global rule requires a tool-provided
operation, parameter, or mode and excludes coordinate recovery assembled by the
caller from unrelated primitives and their own bookkeeping. The row assigns
`partial` specifically when reference coordinates are accepted only to extract
source sequence and no coordinate output accompanies designed variants.

**Operational test declared before tool inspection:** For each tool, check
whether a documented tool operation accepts a chromosome/contig plus reference
positions and whether its designed-sequence output contains reference-genome
locations; score `partial` only when coordinates are accepted solely to extract
source sequence with no coordinate output, `yes` when the tool emits a
reference location for each designed sequence (with or without coordinate
input), and `no` only after repository, documentation, and paper searches show
that neither tool-provided interface exists.

For consistent application:

- A coordinate-bearing CSV field, DataFrame field, or FASTA/name field counts as
  emitted coordinates; the row does not require a particular file format.
- A source interval carried by the tool into each mutated construct counts as
  locating that designed construct's reference origin. It need not imply that
  the mutated bases themselves occur in the unmodified reference.
- Retaining the input BED/VCF beside an unlabelled tensor and joining it back by
  row number is caller bookkeeping and does not count under the global rule.
- Transcript/exon/GTF/GFF behavior is not credited here; it belongs to
  `transcript_annotation_aware`.

## Result

| Tool | Value | Confidence |
|---|---|---|
| poolparty | **partial** [corrected] | high |
| valiant | **yes** | high |
| mpranator | **yes** | high |
| mpradesign | **partial** | high |
| oligopoolcalc | **no** | high |
| mutationmaker | **no** | high |
| dnachisel | **no** | high |
| tangermeme | **partial** | medium |

No cell is `unknown`.

## Primary-source snapshots

Repository evidence was examined at these commits:

- PoolParty local repository: `1bb0179e1c3720b1fffd471802b3040f9336de28`
  (the worktree had a pre-existing README modification; the coordinate source,
  documentation, and tests cited here were unmodified)
- VaLiAnT `develop`: `8796cc112dafd4919fec59913f58cd2be87c45eb`
- MPRAnator: `9969790d62410138d4281b7955da6d085f07b1bc`
- mpradesigntools: `afd386ef12051bb0a37ad63a6f92acd555246584`
- designMPRA: `0cf56ee602fc86dde705906d071a39cbdf6e99a8`
- Oligopool Calculator: `b88fa394ca67ed4c48ec41127b5707694ee7cf0a`
- Mutation Maker: `396c1c0ede7529f3dedf65e821e8c1f20c9a7043`
- DNA Chisel: `68c09304341c3656f3dfe63eda37757d6a7b3917`
- tangermeme: `2006b310cd72a28c56c3ea4ba67f738fff74bb89`

Official repositories:

- <https://github.com/cancerit/VaLiAnT>
- <https://github.com/hemberg-lab/MPRAnator>
- <https://github.com/andrewGhazi/mpradesigntools>
- <https://github.com/andrewGhazi/designMPRA>
- <https://github.com/ayaanhossain/oligopool>
- <https://github.com/ra100/Mutation_Maker>
- <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>
- <https://github.com/jmschrei/tangermeme>

The per-tool records under `tool_survey/final/` were used only as path leads.
No quotation, score, or conclusion from them is evidence in this audit.

## Per-tool audit

### poolparty — `yes` (high confidence)

**Evidence.** PoolParty has a documented reference-coordinate source operation:
“Extract one or more genomic regions from a FASTA file and create a pool” and
“Coordinates are 0-based half-open intervals `[start, stop)` following the
convention `(chrom, start, stop, strand)`”
(`poolparty-statecounter/poolparty/docs/operations/from_fasta.rst:4-8`). The
documented parameter is “A single tuple `(chrom, start, stop, strand)` or a list
of such tuples,” and minus-strand intervals are reverse-complemented (lines
33-39). The implementation accepts that tuple at
`src/poolparty/fixed_ops/from_fasta.py:48-69`, indexes the FASTA via `pyfaidx`,
and extracts `fasta[chrom][start:stop]` at lines 16-34.

Crucially, this is not extraction-only with the coordinate link discarded.
Batch mode itself generates output names
`"{chrom}:{start}-{stop}({strand})"` and passes them into the pool
(`from_fasta.py:145-165`). `FromSeqsOp` returns the generated name as
`seq_name` and contributes it to downstream sequence names
(`src/poolparty/base_ops/from_seqs.py:225-237`). The repository test asserts
that generated library output has a `name` column containing
`chr1:0-4(+)` and `chr1:4-8(-)`
(`tests/test_from_fasta.py:208-222`). Thus each derived member can retain a
tool-generated reference-origin location.

**Behavioral confirmation.** Per the task requirement, I ran the complete
`tests/test_from_fasta.py` module in the read-only Python 3.12 environment with
bytecode and pytest caching disabled; all 23 tests passed. I then ran a targeted
composition from one coordinate-extracted pool into documented sequential
`mutagenize`. The emitted designed variants were:

```text
PYTHONDONTWRITEBYTECODE=1 poolparty-statecounter/.venv/bin/python \
  -m pytest -q -p no:cacheprovider \
  poolparty-statecounter/poolparty/tests/test_from_fasta.py

23 passed
```

```text
       name  seq
chr1:0-4(+) CCGT
chr1:0-4(+) GCGT
chr1:0-4(+) TCGT
```

The test command used `PYTHONDONTWRITEBYTECODE=1`, the provided
`poolparty-statecounter/.venv/bin/python`, an ephemeral FASTA, and removed the
temporary FASTA/index afterward. This settles the otherwise ambiguous question
of whether the coordinate name survives an actual design operation.

**Search and disconfirmation attempt.** I searched the entire local
`src/poolparty`, `docs`, and README trees (excluding binary/image files) with:

```text
(chromosome|genome.coordinates|genomic.coordinates|reference.coordinates|\bBED\b|\bVCF\b|\bGTF\b|\bGFF\b|\blocus\b|\bloci\b|coordinates)
```

The only reference-genome interface found was `from_fasta`; other `coordinate`
hits were sequence-relative molecular/tag coordinate systems. I inspected the
full `from_fasta.py`, its operation page, `from_seqs.py` naming path, and all
`from_fasta` tests. I also searched the authors' `main.pdf` with
`genom|coordinate|chromosome|BED|VCF|FASTA|loci|locus|position|reference sequence`;
it does not document the newer `from_fasta` operation, so the current repository
is the controlling primary source.

Evidence that would have produced `partial` was coordinate input used for
extraction but a designed output that lost the coordinate name; evidence for
`no` was absence of a reference-coordinate operation. I explicitly tested the
first and searched for the second. The preserved downstream `name` rules out
both alternatives.

### valiant — `yes` (high confidence)

**Evidence.** The authors' paper states: “Coordinates for a genomic range are
provided by the user to retrieve a corresponding oligonucleotide reference
sequence” (local `papers/Barbon2022_VaLiAnT_all.pdf`, p. 1). It further says the
BED-like targeton ranges use 1-based `ref_start`/`ref_end`, `ref_chr`, and
`ref_strand` (p. 3).

The current README documents those exact input fields and gives the concrete
row `chrX + 41334132 41334320 ...`
(`README.md:638-656`). The implementation logs and fetches the interval with
`fa.fetch(reference=contig, start=r.start - 1, end=r.end)`
(`src/valiant/sge_utils.py:28-34`).

VaLiAnT also emits the reference linkage for each designed oligonucleotide. Its
documented metadata schema includes `assembly`, `ref_chr`, `ref_strand`,
`ref_start`, `ref_end`, and `mut_position` alongside `oligo_name` and `mseq`
(`README.md:747-780`). The checked-in example at lines 782-786 contains
`GRCh38, ..., chr17, -, 43104080, 43104330, ..., 43104102, ...` in the same row
as the designed oligo sequence. The source declares those fields in
`src/valiant/meta_table.py:48-59` and writes `ref_start`/`ref_end` at lines
234-238.

**Search and disconfirmation attempt.** I searched all 426 files in the
official `develop` archive in memory for
`ref_start|ref_end|ref_chr|ref_strand|fetch_sequence|meta_headers|oligo_name|assembly|vcf.*output|output.*vcf`,
then inspected the README input schema, output schema and example, `sge_utils.py`,
`sge_proc.py`, `meta_table.py`, `vcf_writer.py`, and checked-in SGE output CSVs.
The paper received the same general PDF term scan listed in the audit trail.

Evidence for `partial` would have been coordinate-driven extraction followed by
sequence-only output; evidence for `no` would have been sequence-only input.
I specifically searched the final CSV/VCF write paths and found per-output
assembly and reference coordinates, so neither alternative applies.

### mpranator — `yes` (high confidence)

**Evidence.** MPRAnator's SNP module requires coordinate-bearing FASTA headers.
The primary source says, “The headers in the FASTA file MUST contain the
following information chr:start-end” (`parseSNPs.py:5-6`).
`FastaFileWithHeaderRange.getHeaderRange` parses
`chr(\w+):(\d+)-(\d+)` and raises an error unless the header has that form
(`mycustom.py:167-188`). The parsed chromosome/start/end are used to select VCF
records within each sequence interval (`parseSNPs.py:31-50` and 367-380).

The designed output retains reference coordinates rather than discarding them.
In the actual generation path, when `Named[index]` is the parsed coordinate
tuple, MPRAnator builds `stemHeader` by joining that tuple and appending the SNP
ID, absolute `Position`, and nucleotide change
(`parseSNPs.py:62-88`). The web view places each returned header directly beside
its generated sequence in `finalOutput`
(`iliasApp/views.py:360-367`). A second header builder likewise emits
`chromosome_number`, interval `start`/`end`, SNP ID, absolute `Position`, REF and
ALT (`parseSNPs.py:395-405`). These coordinate-bearing FASTA headers locate each
designed SNP sequence in its reference interval.

**Search and disconfirmation attempt.** I searched all 41 files in the official
repository archive for
`getHeaderRange|chr\w|chromosome|snpInChromosomeRegion|bed|coordinates|positions|header.*range|reference genome|pyfaidx|pysam|twobit|genome build|hg19|hg38`,
then inspected `mycustom.py`, all of `parseSNPs.py`, the called view path in
`iliasApp/views.py`, `forms.py`, and the repository documentation template. I
also searched the authors' two-page paper for the general PDF terms; it
describes SNP design but does not specify this interface.

Evidence for `partial` would have been coordinate headers used only for
filtering/extraction and then removed from the generated FASTA headers. Evidence
for `no` would have been purely sequence-relative SNP input. I traced the
headers through the web result construction and found the original interval
plus absolute SNP position in output, establishing `yes`. The absence of an
internal reference-genome fetcher does not lower the score because the row is
input/output-coordinate based, not a genome-download row.

### mpradesign — `partial` (high confidence)

This score covers both the `mpradesigntools` R package and its `designMPRA`
Shiny wrapper.

**Evidence.** The authors' paper says users “upload VCFs of their variants to
obtain MPRA construct sequences based on the hg38 genome” and, “Our tool
acquires genomic context from the hg38 reference genome, rather than requiring
input by the user” (`papers/Ghazi2018_MPRADesignTools_all.pdf`, pp. 1-2).
The package README says only VCF `CHROM`, `POS`, `REF`, and `ALT` are used
(`README.md:52`) and that context is pulled from the reference genome (lines
58-60).

The code confirms coordinate-only extraction: `processSnp` hard-codes
`genome = BSgenome.Hsapiens.UCSC.hg38`, computes `rangestart`/`rangeend` from
`snp$POS`, and calls
`subseq(genome[[paste0('chr', snp$CHROM)]], start=..., end=...)`
(`R/processVCFfast.R:227`, 326-339). The dependency is declared in
`DESCRIPTION:22` and imported in `NAMESPACE:5`.

The final designed-construct output drops the coordinate fields. Although an
intermediate tibble contains `CHROM = snp$CHROM`
(`R/processVCFfast.R:389-397`), the final package output explicitly selects only
`ID, type, allele, snpIndex, totIndex, barcode, sequence, site_fix_info, notes`
and writes that table (`R/processVCFfast.R:1248-1258`). The Shiny repository does
the same, selecting only `ID, type, allele, snpIndex, totIndex, barcode,
sequence` (`designMPRA/scripts/processVCFfast.R:419-430`) and writing the result
unchanged (`designMPRA/server.R:161-165`). Recovering `CHROM`/`POS` afterward
requires joining the caller's VCF by ID, which the global rule excludes.

**Search and disconfirmation attempt.** I searched all 69 files in
`mpradesigntools` for
`BSgenome|hg38|CHROM|POS|snpseq|subseq\(|write|output|Construct|variant|context`
and all 21 files in `designMPRA` for
`CHROM|POS|hg38|BSgenome|processVCF|download|write.table|write_tsv|output`.
I inspected package input parsing, every SNV/insertion/deletion context fetch,
the package's final select/write path, and the Shiny `processVCF`, download
handler, and UI output description. The author PDF received the same full-text
term scan.

Evidence for `yes` would have been `CHROM`/`POS` or equivalent reference
intervals preserved in the exported construct table; I searched both codebases'
output selection and write paths and found they are explicitly removed.
Evidence for `no` would have been no coordinate-driven reference operation; the
VCF/hg38 `subseq` path rules that out. The row's stated extraction-only case
therefore applies exactly.

### oligopoolcalc — `no` (high confidence)

**Evidence.** The public API is sequence-table based. The documentation defines
ordinary `input_data` as “CSV path or DataFrame with `ID` column + DNA sequence
columns” (`docs/api.md:96`; repeated for the design modules), and `final`
returns only `ID`, `CompleteOligo`, and `OligoLength` (`docs/api.md:596-616`).
The only genome-shaped operation is `background()`, whose input is “list of DNA
strings, CSV/DataFrame with `Sequence` column, or FASTA path”
(`docs/api.md:383-408`). It creates a k-mer screening database, not a
coordinate-indexed extractor. The paper likewise describes a genome sequence
only as template DNA “that must be avoided” during non-repetitiveness screening
(`papers/Hossain2024_OligopoolCalculator.pdf`, p. 6).

**Search and disconfirmation attempt.** I searched all 64 source/documentation
text files in the 76-file official archive for each of:

```text
coordinate | chromosome | \bchrom\b | \bbed\b | \bvcf\b | \bgtf\b |
\bgff\b | \blocus\b | \bloci\b | genomic position | reference position
```

There were zero hits for every term except generic `coordinate` (45 hits).
I inspected all 45: they are barcode-search stream coordinates, read-extraction
offsets, fragment/split coordinates, and internal caches in
`core_barcode.py`, `core_count.py`, `core_pack.py`, `core_split.py`, and
`base/utils.py`; none has a chromosome/contig or reference-genome contract. I
also inspected `README.md`, `docs/README.md`, `docs/api.md`, `docs/docs.md`,
`docs/agent-skills.md`, all public operation signatures and output schemas, and
searched the author paper with the general PDF term scan.

Evidence for `partial` would have been a documented BED/chromosome interval
input that extracts sequence; evidence for `yes` would additionally have been a
reference-location output for designed oligos. Those are the exact terms,
signatures, and output schemas searched, and neither exists. A whole-genome
FASTA used as a k-mer rejection background is explicitly not a coordinate
interface, so the supported value is `no`, not `unknown`.

### mutationmaker — `no` (high confidence)

**Evidence.** Mutation Maker's inputs are raw molecular sequences and
gene/plasmid-relative positions. The paper's SSSM input is “plasmid sequence ...
with gene of interest, a list of amino acid positions to be changed ... forward
and reverse flanking primers”
(`papers/Hiraga2021_MutationMaker.pdf`, p. 2). In source, `Plasmid` has
`gene_start_in_plasmid`, `gene_end_in_plasmid`, and `plasmid_sequence`
(`backend/mutation_maker/ssm_types.py:35-38`), while `SSMSequences` requires
primer and gene strings (lines 96-102). Mutation positions are 1-based codon
indices converted to a gene-relative base offset by
`(one_based_codon_position - 1) * 3`
(`backend/mutation_maker/mutation.py:28-41`).

The GUI accepting `.gb` does not preserve a reference coordinate. Its parser
splits at `ORIGIN`, filters sequence letters, and returns only the resulting
string (`frontend/src/shared/components/FileUploadInput/index.tsx:25-41`);
LOCUS, FEATURES, chromosome, and annotation coordinates are discarded.

**Search and disconfirmation attempt.** Across 170 source/documentation text
files in the official archive (excluding the vendored `feature-viewer` and
`node_modules`), exact case-insensitive searches returned zero hits for each of
`genome`, `chromosome`, `genomic coordinate`, `reference coordinate`, whole-word
`BED`, `VCF`, `GTF`, `GFF`, `locus`, and `loci`. I additionally searched the
complete archive for `gene_start_in_plasmid|gene_end_in_plasmid|codon_position|ORIGIN|FILE_EXTENSIONS`,
inspected `ssm_types.py`, `mutation.py`, the upload component and author paper,
and applied the general PDF scan.

Evidence for `partial` would have been a chromosome/contig plus positions used
to extract a target sequence; evidence for `yes` would have been coordinates in
primer/oligo output. I searched the input schemas, parser, output schemas, and
the explicit genomic-format vocabulary above and found neither. The plasmid and
codon offsets are sequence-relative and therefore do not satisfy this row.

### dnachisel — `no` (high confidence)

**Evidence.** DNA Chisel's `Location` is explicitly sequence-relative:
“Location(5, 10) represents sequence[5, 6, 7, 8, 9]”
(`dnachisel/Location.py:15-24`). `Location.from_data` accepts only a tuple/list,
a Biopython `FeatureLocation`, or another DNA Chisel `Location`, and copies only
start/end/strand (`Location.py:145-186`). There is no contig, assembly, or
reference identifier.

The apparently contrary “genome-wide” example demonstrates caller-side
reconstruction rather than a tool coordinate interface. The user's script
loads a genome, loops over its features, calls
`feature.location.extract(genome_record)`, supplies that raw sequence to one
`DnaOptimizationProblem`, appends records to a list, and uses Biopython
`SeqIO.write` (`examples/common_scenarios/genome-wide-optimization.py:1-32`).
DNA Chisel never receives the reference chromosome/location and cannot emit it.
The paper's phrase “genome-scale gene domestication”
(`papers/Zulkower2020_dnachisel.pdf`, p. 2) points to this example and does not
change that interface fact.

**Search and disconfirmation attempt.** Across all 178 `.py`, `.md`, `.rst`,
`.txt`, YAML, and TOML files in the official archive, exact case-insensitive
searches returned zero hits for `genomic coordinate`, `reference coordinate`,
whole-word `BED`, `VCF`, `GTF`, `GFF`, `locus`, `loci`, and `chromosome`. I also
searched all 281 archive files for
`Location\(|from_data|FeatureLocation|extract\(genome_record\)` plus the
coordinate vocabulary, inspected all of `Location.py`, the core-class docs,
the genome-wide example and the author paper.

Evidence for `partial` would have been a DNA Chisel loader or operation taking
reference chromosome coordinates for extraction. Evidence for `yes` would
have been a tool-produced record carrying reference location after optimization.
I searched for both interfaces and found only sequence-local `Location` plus
the user's Biopython extraction/write loop. Under the global rule that is `no`,
not `partial`.

### tangermeme — `partial` (medium confidence)

**Evidence.** `io.extract_loci` is a genuine documented reference-coordinate
extraction operation. Its `loci` parameter is “Either the path to a bed file or
a pandas DataFrame object containing three columns: the chromosome, the start,
and the end” and `sequences` may be a FASTA path keyed by chromosome
(`tangermeme/io.py:308-320`). It extracts fixed windows around those loci
(lines 269-284), returning sequence and optional signal tensors in locus order
(lines 401-418).

Tangermeme also ships `utils.example_to_fasta_coords`, but its documented
subject is relative spans such as “seqlets or other spans,” which it
cross-references with a separately supplied BED-formatted locus table
(`tangermeme/utils.py:233-258`). Its inputs are two caller-provided DataFrames
(lines 267-282), and it returns corrected `chrom/start/end` spans
(lines 299-308). This is useful coordinate output for analysis annotations, but
it is not connected to designed sequences by a shipped interface.

The design operations return bare edited tensors. For example,
`design.greedy_substitution` documents and returns only
`X: torch.Tensor ... The edited sequence`
(`tangermeme/design/greedy_substitution.py:146-150,244-249`); the other design
modules have the same tensor-only return contract. Supplying a newly authored
`example_df` that identifies designed tensors and reconciling it with the
original loci would be the caller's bookkeeping. The global rule therefore
prevents promoting the generic span mapper to designed-variant coordinate
output.

**Search and disconfirmation attempt.** I searched all 141 files in the
official repository archive for
`example_to_fasta_coords|extract_loci|chromosome|genome.coordinates|genomic.coordinates|reference.coordinates|\bBED\b|\bVCF\b|\bGTF\b|\bGFF\b|\blocus\b|\bloci\b`.
I inspected `io.py`, `utils.py`, `match.py`, `read_vcf`, the I/O API page,
Tutorial C1, `docs/whats_new.rst`, all files under `tangermeme/design/`, and the
author paper. I separately searched every design module's definitions,
`Returns` blocks and return statements; `beam_substitution`,
`greedy_marginalize`, `greedy_substitution`, and `screen` all return tensors,
not coordinate-bearing records.

Evidence for `yes` would have been a design operation returning a compatible
example/span table or directly carrying chromosome/start/end for each edited
sequence; I searched every design return contract and found none. Evidence for
`no` would have been absence of a tool-provided coordinate input; `extract_loci`
rules that out. The binding extraction-only `partial` case is therefore the
best-supported score. Confidence is medium rather than high only because the
generic `example_to_fasta_coords` mapper is adjacent to the threshold; the
global rule resolves it against `yes` because no design output satisfies its
input contract.

## Complete search/audit trail

The following searches were performed in addition to the per-tool searches
reported above:

1. Read `ROW_DEFINITIONS.md:1-177` before tool inspection and declared the
   operational test above.
2. Used the eight `tool_survey/final/<slug>.md` files only to identify likely
   primary-source paths. Their conclusions and ratings were disregarded and
   are not cited.
3. Extracted every author PDF with PyMuPDF, as required, and searched the text
   with
   `genom|coordinate|chromosome|BED|VCF|FASTA|loci|locus|position|reference sequence`.
   This covered PoolParty `main.pdf` plus the seven PDFs in
   `tool_survey/papers/`. Targeted follow-ups searched the exact phrases
   `Coordinates for a genomic range`, `Targeton ranges are defined`,
   `acquires genomic context`, `genome sequence of E. coli`,
   `User's input includes plasmid sequence`, and
   `genome-scale gene domestication` to recover page-local quotations.
4. Locator-only web searches were:
   `site:github.com/cancerit/VaLiAnT ref_start ref_end ref_chr`,
   `site:github.com/hemberg-lab/MPRAnator FastaFileWithHeaderRange getHeaderRange`,
   `site:github.com/andrewGhazi/mpradesigntools BSgenome.Hsapiens.UCSC.hg38 processSnp`,
   and `site:github.com/jmschrei/tangermeme example_to_fasta_coords extract_loci`.
   Search indexing was incomplete, so evidence was read from the official
   repositories directly.
5. Direct raw-file opens succeeded for DNA Chisel `Location.py` and tangermeme
   `io.py`/`utils.py`; several other raw URLs returned cache misses. To avoid
   treating search-engine snippets or mirrors as evidence, I fetched official
   GitHub branch archives into memory, performed the listed greps/counts, and
   wrote no repository copies to disk.
6. `git ls-remote ... HEAD` supplied the source snapshot identifiers above.
   PoolParty's identifier came from its local read-only checkout.
7. The failed first Oligopool archive fetch used the nonexistent `main` branch;
   the official GitHub metadata endpoint reported `master`, after which the
   complete 76-file archive was searched successfully.

## Row-level finding

The row discriminates on one consistent scale: **3 yes / 2 partial / 3 no / 0
unknown**. The decisive distinction is not whether a tool mentions genomic
data, nor whether it has sequence-relative positions. It is whether a shipped
interface consumes reference coordinates and, for full credit, preserves or
emits the reference linkage with designed sequences.

This audit also exposes two easy misclassifications. First, a coordinate string
in a generated FASTA/name field is still coordinate output; applying the same
threshold to MPRAnator and PoolParty makes both `yes`. Second, a generic
coordinate conversion utility does not by itself locate designed variants when
the design API returns an incompatible bare tensor; applying the global rule to
tangermeme keeps it `partial`. MPRA Design Tools is the row's clearest literal
`partial`: VCF coordinates drive hg38 extraction, but both exported construct
tables explicitly drop `CHROM` and `POS`.
