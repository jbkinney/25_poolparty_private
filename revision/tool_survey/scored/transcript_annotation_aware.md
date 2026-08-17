# Audit: `transcript_annotation_aware`

Date of audit: 2026-08-15


> ## CORRECTION applied after scoring — VaLiAnT promoted
>
> **`VaLiAnT`: `partial` → `yes`.** The audit's own evidence establishes full
> annotation-aware design: it loads GTF/GFF2, builds
> `Annotation(TranscriptInfo, utr_ranges, exons)`, and the design pipeline consumes
> it — `Transcript.get_cds_seq` retrieves partial-codon bases from preceding and
> following exons (`src/valiant/transcript.py:84-146`). The scorer's own summary:
> *"This is direct annotation-aware design, not coordinate extraction."*
>
> The only restriction cited is one transcript per gene
> (`loaders/gtf.py:148-156`). That matches **neither** `partial` branch in the row
> definition — annotation is not consumed for a single purpose (it drives exon
> detection, reading frame, split codons and CDS-specific mutators), and exon
> structure is not handled without a transcript model. A documented scope limit is
> not a partial capability.
>
> Direction of the correction is worth noting: it moves **against** PoolParty, which
> scores `no` on this row.

## Measurement instrument and operational test

I read the preamble, including **“Global rule — binding on every row below,”**
and row 17 in `25_poolparty_private/revision/tables/ROW_DEFINITIONS.md` before
examining any tool. In particular, I applied the rule that a capability must be
a tool-provided operation, parameter, or mode; a user-created composition of
generic primitives and bookkeeping is `no`, not `partial`.

**Operational test stated before examining the tools:** for each tool, check
whether a documented design operation, parameter, or mode directly consumes or
represents transcript models, exon/intron structure, or GTF/GFF annotation and
uses that information to determine the design; score `yes` for general
transcript/annotation-aware design, `partial` when the tool-provided support is
limited to one annotation purpose, exon structure without transcript models, an
example-only interface, or another material scope/type restriction, and `no`
only after source/docs/paper searches show no such operation (coordinates alone
do not count).

## Result

| Tool | Value | Confidence |
|---|---|---|
| PoolParty | **no** | high |
| VaLiAnT | **yes** [corrected] | high |
| MPRAnator | **no** | high |
| MPRA Design Tools | **no** | high |
| Oligopool Calculator | **no** | high |
| Mutation Maker | **no** | high |
| DNA Chisel | **no** | high |
| tangermeme | **no** | high |

VaLiAnT is the sole positive cell, but it is `partial` because its documented
GTF/GFF2 interface and loader reject more than one transcript per gene. That is
a material restriction on transcript-model support under the global scoring
rule. It nevertheless does substantially more than consume coordinates: it
constructs a `Transcript` with ordered exons and uses adjacent-exon bases,
strand, frame, CDS/UTR identity, and transcript lifting during design.

## Source versions and common search protocol

The per-tool `final/*.md` records were used only as URL/file-location leads.
None is cited as evidence and none of their ratings was used.

Repository snapshots inspected (all official repositories) were:

| Tool | Repository and branch | Commit searched |
|---|---|---|
| PoolParty | `jbkinney/poolparty-statetracker`, `main` | protected local checkout `1bb0179e1c3720b1fffd471802b3040f9336de28`; current remote head `9d6a2058a5d0d4742d0874775b2f3a3de2f0706a` was also cloned read-only to `/tmp` and gave the same row-specific search result |
| VaLiAnT | `cancerit/VaLiAnT`, `develop` | `8796cc112dafd4919fec59913f58cd2be87c45eb` |
| MPRAnator | `hemberg-lab/MPRAnator`, `master` | `9969790d62410138d4281b7955da6d085f07b1bc` |
| MPRA Design Tools | `andrewGhazi/mpradesigntools`, `master`; `andrewGhazi/designMPRA`, `master` | `afd386ef12051bb0a37ad63a6f92acd555246584`; `0cf56ee602fc86dde705906d071a39cbdf6e99a8` |
| Oligopool Calculator | `ayaanhossain/oligopool`, `master` | `b88fa394ca67ed4c48ec41127b5707694ee7cf0a` |
| Mutation Maker | `ra100/Mutation_Maker`, `master` | `396c1c0ede7529f3dedf65e821e8c1f20c9a7043` |
| DNA Chisel | `Edinburgh-Genome-Foundry/DnaChisel`, `master` | `68c09304341c3656f3dfe63eda37757d6a7b3917` |
| tangermeme | `jmschrei/tangermeme`, `main` | current remote head `2006b310cd72a28c56c3ea4ba67f738fff74bb89` |

`git ls-remote` was used to verify the official branch heads. Seven original
snapshots matched. PoolParty and tangermeme had advanced; fresh read-only clones
of those two heads were searched. The protected PoolParty checkout was not
modified.

Every PDF was converted to text with the instructed PyMuPDF method:

```bash
python3 -c "import fitz; d=fitz.open('PATH'); print(chr(10).join(p.get_text() for p in d))"
```

The same structural-term search was run over every tracked repository tree
(`-I` skips binary contents; generated minified JS/CSS, SVG, and lock files were
excluded only from the displayed results):

```bash
git grep -n -I -P '(?i)\b(gtf|gff2?|exons?|introns?|transcripts?|splice|splicing|spliced)\b|gene[ _-]?models?'
```

A second full-tree search was run to catch tools using a different sense of
“annotation”:

```bash
git grep -n -I -P '(?i)\bannotat(?:e|ed|es|ing|ion|ions)\b'
```

The extracted papers were searched with:

```bash
rg -n -i '\b(GTF|GFF2?|exons?|introns?|transcripts?|splice|splicing|spliced)\b|gene[ _-]?models?' PAPER.txt
```

I also searched each documented input/API surface for format and carrier terms
(`FASTA|VCF|BED|GTF|GFF|input|sequence|region|annotation`) and inspected the
matching definitions rather than treating a keyword hit as support. The
per-tool sections below give the complete result of these searches, the files
and documentation pages inspected, and the disconfirmation test.

## Per-tool audit

### PoolParty — **no** (confidence: high)

**Quoted evidence.** The documented region annotation is caller-authored by
numeric extent:

> “Tag a region by position range”

and its `extent` is:

> “`(start, stop)` tuple (0-based, exclusive stop).”

Source: official repository,
`poolparty/docs/operations/annotate_region.rst:4-7,34-40` at protected checkout
`1bb0179e...` (the same text is present at current remote head `9d6a2058...`).

The underlying scan interface is likewise generic:

> “`positions`: Valid insertion positions (0-based). If None, all positions are valid.”
>
> “`region`: Region to constrain the scan to. Can be region name (str) or [start, stop].”

Source: `poolparty/src/poolparty/region_ops/region_scan.py:15-46` at the current
official head. PoolParty does have `annotate_orf`, but its docs require the user
to supply one contiguous `extent` and a `frame` (`annotate_orf.rst:4-8,34-47`);
that is ORF/codon handling, not a transcript model or exon structure.

The authors' paper contains a splicing application, but describes caller-chosen
positions around one already-known boundary:

> “Using the 5’ss of SMN2 exon 7 as the canonical site, we used PoolParty to
> construct a library in which cryptic 5’ss 9-mers of varying strength were
> inserted at varying positions on either side of the canonical 5’ss.”

Source: `25_poolparty_private/26.05.18_bmc_submission/latex/main.pdf`, §3.3,
manuscript lines 234-239 (PyMuPDF extraction). The same section says the GAMs
were separately fit for exonic and intronic positions; it does not expose a
PoolParty transcript/annotation input. Under the global rule, knowing the SMN2
boundary and labeling numeric positions is the user's biological bookkeeping.

**Searches and disconfirmation.** I searched all 356 tracked files at the
current head with both common patterns above. The structural query hit only
`poolparty/tests/test_integration_libraries.py` (`exon`, `intron`, and a splice
screen test); it had zero hits in `src/`, `README.md`, or textual `docs/`.
Annotation hits were inspected in `README.md`, `docs/api.rst`,
`docs/operations/{annotate_region,annotate_orf,clear_annotation,translate}.rst`,
`src/poolparty/{region_ops,orf_ops}/`, and related mixins: they are XML-like
user-defined regions, ORFs, display styles, and design-card annotations. I also
searched the entire PDF; its row terms are confined to the SpliceAI use case and
references, with no GTF/GFF/transcript-model interface. I inspected the region
operation signatures and the paper's complete §2.8/§3.3 discussion.

Evidence that would have changed this cell: a documented PoolParty operation
that parses GTF/GFF or a transcript/exon carrier and automatically constrains
design would yield `yes` (or `partial` if single-purpose/example-only); a tool
operation that at least derives exon structure from supplied annotation would
yield `partial`. I explicitly searched for both and found neither. Generic
regions plus user-known exon coordinates do not satisfy the test.

### VaLiAnT — **partial** (confidence: high)

**Quoted evidence.** The documented SGE input includes:

> “features file (GTF/GFF2, optional)”

and the manual states:

> “The features file (SGE-only `gff` option) is required to detect exonic
> regions in the targeton ... The files should only contain features for one
> transcript per gene ... Any features of type other than `CDS` and `UTR` are
> ignored.”

Source: official `develop` README,
`README.md:76-106`, commit `8796cc112daf...`. The CLI table repeats:

> “Path to GTF/GFF2 file containing CDS and UTR features; one transcript per
> gene only.”

Source: `README.md:227-239`.

The loader genuinely builds a transcript model. It collects CDS/UTR features,
validates `(gene_id, transcript_id)`, converts CDS features to exons, and returns
`Annotation(TranscriptInfo(...), utr_ranges, exons)`
(`src/valiant/loaders/gtf.py:126-169`). The hard restriction is executable code:

> `if len(id_pair) > 1:`
>
> `    raise ValueError("Multiple transcript in annotation: not supported!")`

Source: `src/valiant/loaders/gtf.py:148-156`.

The design pipeline then uses the result rather than merely retaining it:

> `annot = GtfLoader(contig, strand).load_gtf(config.gff_fp)`
>
> `transcript = Transcript(annot.transcript_info, annot.cds)`
>
> `ctx = get_ctx_range(transcript.exons + [x.ref for x in tcs])`

and, with background variants, lifts and stores the exon model
(`src/valiant/sge_proc.py:448-473`). `Transcript.get_cds_seq` retrieves missing
partial-codon bases from preceding/following exons
(`src/valiant/transcript.py:84-146`). This is direct annotation-aware design,
not coordinate extraction.

The paper confirms the same behavior:

> “To collect gene and transcript information and to apply CDS-specific mutator
> functions, appropriate transcript annotation must be provided via a GTF/GFF2
> file; only CDS, UTR and stop features are taken into consideration. One
> transcript per gene is allowed ... Targeton region reading frame is computed
> using user-specified transcript feature annotations.”

Source:
`25_poolparty_private/revision/tool_survey/papers/Barbon2022_VaLiAnT_all.pdf`,
§2.4.2, journal p. 895.

**Searches and disconfirmation.** Both common searches were run across all 426
tracked files. Real structural hits occur in the README, CLI, GTF loader,
`annotation.py`, `exon.py`, `transcript.py`, `sge_proc.py`, `cdna_proc.py`, SQL
schema/queries, tests, and shipped GTF examples. I inspected those implementation
paths plus the paper's §2.4.2, CDS-frame/partial-codon discussion, negative-strand
handling, and BRCA1 examples. The PDF has extensive GTF/GFF, transcript, exon,
intron, and splice hits, not merely references.

Evidence that would have changed the value: removal of the single-transcript
guard, with documented support for multiple alternative transcript models per
gene, would produce `yes`; absence of a transcript object or use of GTF only for
coordinate extraction would produce `no` or at most a narrower `partial`. I
searched the loader, CLI/manual, model class, and design call sites and confirmed
both the substantial support and the explicit restriction. Therefore `partial`
is not a hedge: it records a concrete accepted-model limitation.

### MPRAnator — **no** (confidence: high)

**Quoted evidence.** The input parser's actual purpose-built classes are
`FastaFile`, `FastaFileWithHeaderRange`, `BedFile`, and `SnpFile`. For example:

> “This class represents a FastaFile.”

and `FastaFileWithHeaderRange` parses headers of the form:

> `chr15:12900-12950`

Source: official repository `mycustom.py:16-39,145-188`, commit
`9969790d6241...`. The user-facing SNP form is explicitly:

> `SnpS = SnpTextField(label="Enter your SNPs (VCF Format)")`

Source: `iliasApp/forms.py:501-502`. These are sequences/coordinates/variants;
none represents transcript or exon structure.

**Searches and disconfirmation.** The two common queries were run over all 41
tracked files, including all Python, templates, downloadable scripts, and docs;
both returned zero row-relevant hits. I specifically inspected
`mycustom.py`, `iliasApp/forms.py`, `iliasApp/templates/iliasApp/docs.html`,
`part1.py`, `parseSNPs.py`, `part3.py`, both downloadable scripts, and the full
paper. The docs/input search found FASTA, BED, and VCF interfaces only. The PDF
structural-term query returned zero hits.

Evidence that would have changed the value: a GTF/GFF upload/parameter, a
transcript/exon model in `mycustom.py`, or a documented mode that chooses edits
using exon/intron identity would give `yes` or restricted `partial`. I searched
the complete tree, documentation template, forms, parser classes, and paper for
exactly those interfaces and confirmed their absence.

### MPRA Design Tools — **no** (confidence: high)

**Quoted evidence.** The package documents its input semantics narrowly:

> “Only the CHROM, POS, REF, and ALT columns are used. The INFO column is used
> only for detecting reverse strand constructs.”

It further requires caller bookkeeping:

> “If the user wishes to generate SNPs for genes that normally are read from
> the reverse strand, add a string containing `MPRAREV` to the INFO field...”

and:

> “the MPRAREV tag will need to be added by the user (where appropriate) because
> the VCF's do not always specify which strand the relevant gene is on.”

Source: official `mpradesigntools` `README.md:45-61`, commit
`afd386ef1205...`.

The public R function takes a VCF plus numeric upstream/downstream context; its
documented signature has no annotation carrier (`R/processVCFfast.R:997-1043,
1099-1115`). The paper similarly says:

> “Users can adjust experimental parameters ... as well as upload VCFs for
> automated construct sequence generation.”

Source:
`25_poolparty_private/revision/tool_survey/papers/Ghazi2018_MPRADesignTools_all.pdf`,
abstract and §2 (PDF pp. 1-2).

**Searches and disconfirmation.** I searched all 69 tracked files in
`mpradesigntools` and all 21 in `designMPRA`, including every `.R`, `.Rd`, `.Rmd`,
README, UI/server file, and supplement. There were zero structural hits in the R
package. The Shiny repository's sole hit was “transcripts” in prose explaining
MPRA measurement variance (`ui.R:78`), not a transcript model. Annotation-term
searches returned zero hits in both repos. The complete paper search found no
GTF/GFF/exon/intron/splice/transcript-model interface; its relevant input terms
are VCF and hg38 reference context.

Evidence that would have changed the value: automatic gene-strand/exon handling
from a transcript annotation would yield `yes` or restricted `partial`.
I searched both repositories, all package manuals, Shiny UI/server, supplement,
and paper for it. The explicit need for the caller to add `MPRAREV` confirms that
gene orientation is not annotation-derived.

### Oligopool Calculator — **no** (confidence: high)

**Quoted evidence.** Its documented design carrier is sequence columns, not
genome annotation:

> “`input_data` (`str` / `pd.DataFrame`): Path to a CSV file or DataFrame with
> annotated oligo pool variants.”

Here “annotated” means table columns; the notes define the contents:

> “`input_data` must contain a unique 'ID' column; all other columns must be
> non-empty DNA strings.”

Source: official repository `oligopool/barcode.py:17-45,62-76`, commit
`b88fa394ca67...`. The core docs agree:

> “Input: CSV path or pandas DataFrame with an `ID` column”

Source: `docs/docs.md:136-158`.

**Searches and disconfirmation.** Both common searches were run over all 76
tracked files, including `README.md`, `docs/{README,api,docs,agent-skills}.md`,
all public modules under `oligopool/`, examples, CLI/YAML docs, and the paper.
The structural query returned zero repository hits. Annotation hits were all
phrases such as “annotated oligo pool variants” in DataFrame-oriented APIs;
inspection showed no biological feature parser. The paper's only structural
hits were application background (“alternative transcript polyadenylation”) and
a Salmon citation, not a design input or operation.

Evidence that would have changed the value: a design module accepting GTF/GFF,
an exon/transcript column schema that changes design, or a transcript-model
object would give `yes` or restricted `partial`. I searched the complete source
and docs tree, public signatures, CLI configuration documentation, examples, and
paper and found none. A general DataFrame called “annotated” does not meet the
row definition.

### Mutation Maker — **no** (confidence: high)

**Quoted evidence.** The documented UI collects caller-supplied contiguous
sequences:

> “Enter or Upload Plasmid sequence containing flanking primers and Gene of
> Interest sequence.”
>
> “Enter or Upload Gene of Interest sequence.”

and mutations are position strings:

> “Enter mutations in following format [Amino Acid Codon][Location][X]”

Source: official author-maintained repository
`frontend/src/scenes/SSM/components/SSMForm.tsx:156-218`, commit
`396c1c0ede75...`. Backend parsing converts the caller's one-based codon number
directly to a base offset:

> `zero_based_base_position = (one_based_codon_position - 1) * 3`

Source: `backend/mutation_maker/mutation.py:28-41`.

**Searches and disconfirmation.** Both common searches were run over all 300
tracked files, including backend Python, frontend TS/TSX/JS/JSX, READMEs, and
docs/plans. The only repository structural hit was JavaScript array
`newData.splice(index, 1, ...)`, unrelated to RNA splicing. Annotation hits were
an SSM docstring and third-party licensing text; neither is a genome annotation
interface. I inspected the SSM/QCLM/PAS form inputs, backend mutation types and
API request types, top-level/backend READMEs, and the full paper. Paper “splicing”
hits mean SOEing (“splicing by overlap extension”), not transcript splicing; no
GTF/GFF/exon/intron/transcript model occurs.

Evidence that would have changed the value: accepting a genomic annotation and
mapping targets through its exon/transcript structure would give `yes`, while a
documented exon-only mode would give `partial`. I searched all implementation
languages, UI inputs, backend API/types, documentation, and paper for those and
confirmed that the tool instead requires a pre-built gene sequence and explicit
codon locations.

### DNA Chisel — **no** (confidence: high)

**Quoted evidence.** DNA Chisel does read “annotations,” but its parser defines
them as DNA Chisel commands in GenBank `misc_feature`s:

> “Create a DnaOptimizationProblem by parsing a record's annotations.”

The implementation then does:

> `if feature.type != "misc_feature":`
>
> `    continue`

and only keeps labels recognized as specifications
(`DnaOptimizationProblem/mixins/RecordRepresentationMixin.py:16-74`, official
commit `68c09304341c...`). The label parser requires the label/note to start with
`@` or `~` (`dnachisel/biotools/genbank_operations.py:195-209`). Thus ordinary
biological `exon`/`CDS` annotations are not interpreted as transcript structure.

Its `@cds` documentation is likewise a caller-marked contiguous span:

> “To indicate that a region is a CDS and the protein sequence should be
> conserved ... use @cds ... on a region whose span is a multiple of 3.”

Source: `docs/genbank/genbank_api.rst:168-186`. This is codon-aware constraint
placement, not an exon model.

**Searches and disconfirmation.** Both common queries were run over all 281
tracked files, including `README.rst`, all `docs/**/*.rst`, all package modules,
examples, and the paper. There were no GTF/GFF/exon/intron/transcript-model hits
in textual documentation or package code. Code's lone “splicing” hit is Python
slice terminology in `Location.py:19`. A shipped GenBank example contains an
ordinary `exon` feature (`examples/common_scenarios/constraints_breaches_report/
genbanks/HC_Amp_ccdB.gb:21`), but the quoted parser skips it because it is not a
`misc_feature`; this is direct disconfirmation, not support. Other transcript
hits are free text inside a whole-E. coli GenBank example. The paper search had
zero structural hits. I inspected `docs/genbank/genbank_api.rst`,
`docs/ref/{core_classes,builtin_specifications,biotools,cli}.rst`, the record
parser, feature-label parser, built-in specification set, README, and paper.

Evidence that would have changed the value: parsing ordinary exon/CDS joins or
GTF/GFF into a transcript/exon carrier and using it to scope optimization would
give `yes`; a supported single-exon/exon-only mode would give `partial`. I
searched those exact terms and inspected the only biological exon feature found.
GenBank as a container for user-authored optimization commands is not transcript
annotation awareness under the global rule.

### tangermeme — **no** (confidence: high)

**Quoted evidence.** Its coordinate input is explicitly BED plus FASTA:

> “Either the path to a bed file or a pandas DataFrame object containing three
> columns: the chromosome, the start, and the end...”
>
> “Either the path to a fasta file ... or a dictionary where the keys are ...
> chromosomes...”

Source: current official repository `tangermeme/io.py:246-320`, commit
`2006b310cd72...`. Its VCF helper simply returns the first nine VCF columns
(`io.py:622-655`); it does not interpret transcript consequences.

Tangermeme also has an `annotate` module, but its own definition is motif
annotation:

> “Annotate a set of seqlets according to a motif database using TOMTOM.”

Its “BED-formatted” seqlet coordinates are example-relative, and the third input
is a motif dictionary or MEME file (`tangermeme/annotate.py:18-54`). This is not
GTF/GFF or exon/transcript structure.

**Searches and disconfirmation.** Both common searches were run over all 141
tracked files at the current head, including package Python, `README.md`, all
textual docs/API pages, bundled skill docs, examples/notebooks, tests, and the
paper. The structural query returned zero repository hits. Annotation hits were
inspected in `annotate.py`, `docs/api/annotate.rst`, `docs/index.rst`, annotation
tutorials, result types, plotting, and design modules; all denote motifs,
seqlets, attribution spans, or generic span metadata. The complete paper's only
row-like term was “alternative splicing” in scientific background/references;
its annotations are explicitly motif matches. I also inspected
`io.extract_loci`, `io.read_vcf`, `variant_effect`, and the documented design
operations to test whether coordinate/VCF inputs were joined to transcript data;
they were not.

Evidence that would have changed the value: a documented design operation that
accepts GTF/GFF/transcript objects or automatically masks/changes edits by
exon/intron identity would give `yes` or restricted `partial`. I searched all
source/docs/paper surfaces for that and found none. BED extraction and motif
annotation belong to adjacent coordinate/motif capabilities and are not credited
again here.

## Row-level discrimination finding

The row is applicable on one consistent scale and does discriminate, although
sparsely: seven tools are `no`, while VaLiAnT is a substantive but restricted
`partial`. The important discriminant is not whether a paper discusses genes,
splicing, “annotations,” or genomic loci; it is whether the tool itself carries
and uses transcript/exon structure in design. Applying the global rule prevents
PoolParty's caller-defined regions, MPRA Design Tools' caller-added `MPRAREV`,
DNA Chisel's command annotations, and tangermeme's motif annotations from being
mistaken for transcript awareness. No replacement row was needed.
