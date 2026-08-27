# Audit: recombination

Audited 2026-08-15. This audit was performed independently from the binding
definition and primary sources. Files in `tool_survey/final/` were used only in
an initial locator pass; no statement or score from them is evidence here.

## Measurement instrument and operational test

I read the preamble and **Global rule — binding on every row below** in
`revision/tables/ROW_DEFINITIONS.md`, then row 10 and row 5's binding boundary
rules. The global rule requires a tool-provided operation, parameter, or mode;
caller-written composition and bookkeeping are `no`, not `partial`. Row 10 is:

> “Construction of variants by joining segments from two or more parent
> sequences at breakpoints the tool chooses or enumerates. Concatenating
> user-supplied fragments does not count; the tool must generate the breakpoint
> combinations. This covers both the wet-lab sense (chimeric constructs) and
> the in-silico sense (swapping segments between parents to ask which region
> drives a model's prediction).”

Row 5's adjacent-capability rules were also enforced: multi-motif placement is
placement rather than mutational combination; IUPAC expansion is template
decompression; and an optimiser's multi-edit output is not an enumerated or
sampled declarable combinatorial mutation space. Those rules prevent motif
placement, degenerate templates, generic multi-mutation libraries, and
optimizer-found edits from being repurposed here as “recombination.”

**Operational test stated before inspecting any tool:** For each tool, does a
documented tool-provided operation, parameter, or mode accept at least two
parent sequences and construct chimeric variants by itself choosing or
enumerating segment-joining breakpoints, excluding user-authored concatenation,
pre-specified fragments, generic multi-mutation combinations, motif placement,
IUPAC expansion, and caller-written bookkeeping, with a restricted or
example-only implementation worth at most `partial`?

I applied that sentence unchanged to all eight tools.

## Results

| Tool | Score | Confidence | Decisive basis |
|---|---:|---:|---|
| PoolParty | **yes** | high | Documented `recombine` accepts at least two equal-length parent pools and enumerates or samples breakpoint × source-assignment chimeras. |
| VaLiAnT | **no** | high | Its documented mutators act on one retrieved reference/cDNA sequence; no operation takes multiple parents or generates crossover breakpoints. |
| MPRAnator | **no** | high | It places motifs or SNP combinations into a background and concatenates user-supplied construct parts; neither path is parent-segment recombination. |
| MPRA Design Tools | **no** | high | `processVCF` builds ref/alt constructs one VCF row at a time from one genomic context and fixed supplied parts. |
| Oligopool Calculator | **no** | high | `join` reconciles columns for identical row IDs and `split` fragments one existing oligo for assembly; neither constructs chimeras from parent sequences. |
| Mutation Maker | **no** | high | SSSM, QCLM, and PAS each expose one gene of interest; PAS chooses fragment boundaries within that single sequence. |
| DNA Chisel | **no** | high | A problem optimizes one starting sequence; region compatibility and choice specifications do not join multiple parents or enumerate crossovers. |
| tangermeme | **no** | high | `insert`/`substitute` accept a caller-supplied motif and one placement; multi-motif design is placement, not breakpoint generation between parents. |

## Tool audits

### PoolParty — yes (high confidence)

**Positive evidence.** The documented operation says:

> “Produce chimeric sequences by slicing multiple source pools at breakpoints
> and stitching together alternating segments. Requires at least two source
> pools of equal sequence length. Use `mode='sequential'` to enumerate all
> breakpoint combinations deterministically, or `mode='random'` (default) to
> draw chimeras stochastically.”

(`poolparty-statecounter/poolparty/docs/operations/recombine.rst:4-8`, commit
`1bb0179e1c3720b1fffd471802b3040f9336de28`). The same page documents `sources`
as at least two pools, `num_breakpoints`, optional `positions`, and modes where
random “selects breakpoints randomly” and sequential “enumerates all breakpoint
positions” (`:34-48,65-74`). Its worked example calls
`pp.recombine(sources=[src_a, src_b], mode="sequential")` and reports 18
single-crossover chimeras for two ten-base parents (`:95-119`).

The implementation confirms that the tool, not the caller, creates the
combinations. `RecombineOp` states:

> “In sequential mode, enumerates all breakpoint positions × pool assignment
> combinations. In random mode, randomly selects breakpoints and pool
> assignments.”

(`src/poolparty/base_ops/recombine.py:143-148`). When `positions` is omitted the
operation supplies every internal position (`:205-209`); sequential mode counts
`C(P,K) × N × (N-1)^K` states (`:241-252`), loops over
`combinations(self.positions, self.num_breakpoints)` and every valid source-pool
assignment (`:301-319`), and slices the selected source at each breakpoint
before adding the final segment (`:421-451`). This is exactly the row's
two-parent, tool-generated-breakpoint construction.

**Behavioral check.** Required read-only execution used the repository venv
with `PYTHONDONTWRITEBYTECODE=1`. For parents `AAAA` and `TTTT`, one breakpoint,
and sequential mode, the operation materialized six rows:

```text
ATTT  breakpoint=(0,) sources=(0,1)
TAAA  breakpoint=(0,) sources=(1,0)
AATT  breakpoint=(1,) sources=(0,1)
TTAA  breakpoint=(1,) sources=(1,0)
AAAT  breakpoint=(2,) sources=(0,1)
TTTA  breakpoint=(2,) sources=(1,0)
num_rows=6
```

The complete targeted test module also passed: `44 passed in 0.77s` for
`tests/test_recombine.py`.

**Disconfirmation attempt.** I searched the public docs, `src/`, and `tests/`
for `recombine`, `breakpoint`, `pool_assign`, `chimera`, and the adjacent
`join`/concatenation operation, then read the complete constructor,
enumeration, sampling, and segment-assembly paths. A different value would have
resulted if `recombine` merely concatenated caller-supplied fragments, required
the caller to enumerate each breakpoint construct, or existed only in an
example. Those possibilities were explicitly searched and contradicted by the
documented API, source, live output, and passing operation tests.

### VaLiAnT — no (high confidence)

**Evidence of implemented scope.** VaLiAnT describes itself as a library design
tool for saturation editing (`README.md:1-5`, commit
`8796cc112dafd4919fec59913f58cd2be87c45eb`). Its two modes obtain sequence from
one reference source per target: SGE uses reference-genome coordinates and cDNA
DMS uses user-provided cDNA sequence (`README.md:189-196`). Its complete
documented mutation-type list is parametric deletion, SNV, in-frame deletion,
alanine/stop/all-amino-acid codon substitution, SNVRE, and custom variants
(`README.md:255-270`). The targeton specification declares mutation types to be
applied to regions of that reference sequence (`README.md:630-649`). None is a
multi-parent or crossover operation.

The author paper likewise says that for each targeton **reference sequence**
mutator functions are applied (Barbon et al. 2022, §2.2), and more decisively:

> “Discrete mutation events are applied per oligonucleotide so that each
> generated sequence contains specific PAM/protospacer protection edits and
> the desired variant composed of one or more base changes.”

(Barbon et al. 2022, §2.4.2). It explains that final variant oligos are made by
replacing target regions in a PAM/protospacer-protected reference and assembling
that result with invariant regions/adaptors (§2.5.2). This is mutation and fixed
construct assembly around one parent, not recombination between parents.

One lexical hit in the repository does not change the score: the example
GenBank file `examples/cdna/construct_generation/BRCA1_cdna_pCW57.1_model.gb`
contains the annotation “recombination site for the Gateway(R) LR reaction” at
line 2121. It is input biological annotation, not a VaLiAnT operation.

**Disconfirmation attempt / absence search.** I searched all repository source,
tests, README, changelog, examples, and the full paper for `recombin*`,
`chimer*`, `cross-over`/`crossover`, `breakpoint`, `mosaic`, `parent sequence`,
and segment/swap combinations; after excluding the GenBank annotation, there
were no operation or documentation hits. I separately searched adjacent terms
`join`, `merge`, `concat`, `splice`, `assemble`, `fragment`, `overlap`,
`replace`, and `insert`, and read the CLI inputs, mutation-type, targeton,
custom-variant, and output sections. A different score required a supported
parameter or mode accepting at least two reference/cDNA parents and generating
crossover positions. I searched source, documentation, examples, and paper for
that exact mechanism and found only one-reference mutagenesis and fixed construct
parts.

### MPRAnator — no (high confidence)

**Evidence of implemented scope.** The paper enumerates the supported scope as
four investigations: motif placement, SNP effects, PWM realization, and
Transmutation negative controls (Georgakopoulos-Soares et al. 2017, §2). Its
motif and SNP combinatorics are adjacent capabilities, not recombination: motif
design places single or combinations of motifs at preselected positions; SNP
design applies single or combinations of SNPs to each provided background
sequence. The shipped documentation is equally explicit that the motif user

> “is able to specify the locations to substitute the motifs, insert
> restriction sites, adapter sites and barcodes.”

(`iliasApp/templates/iliasApp/docs.html:79-82`, commit
`9969790d62410138d4281b7955da6d085f07b1bc`). Thus even the positions are
caller-specified motif placement, which row 5's boundary assigns away from
recombination.

The SNP option can enumerate combinations of mutation records
(`parseSNPs.py:20-30,36-62`), but these are multiple events applied within one
background, not segments drawn from multiple parents. Construct assembly also
does not pass: `createMPRAResultOutput` loops through a caller-selected
`ordering` and appends adapter, restriction-site, background, and barcode
strings (`part1.py:210-292`). That is the row's expressly excluded
concatenation of user-supplied fragments.

**Disconfirmation attempt / absence search.** I searched every repository
Python and HTML file plus the paper for `recombin*`, `chimer*`, `crossover`,
`breakpoint`, `mosaic`, `parent sequence`, and segment/swap phrases; there were
no hits. I then searched `join`, `merge`, `concat`, `splice`, `assemble`,
`fragment`, `overlap`, `substitute`, `background`, and the four module names,
and read `oligo.py`, `part1.py`, `parseSNPs.py`, `part3.py`, the shipped docs,
and the paper's module descriptions. A different score required an operation
that accepts at least two background parents and itself selects or enumerates
where to switch from one to the other. I explicitly searched for it and found
only SNP-event combinations, motif placement, and caller-directed construct
concatenation.

### MPRA Design Tools — no (high confidence)

This score covers both primary repositories named for the tool:
`andrewGhazi/mpradesigntools` commit
`afd386ef12051bb0a37ad63a6f92acd555246584` and
`andrewGhazi/designMPRA` commit
`0cf56ee602fc86dde705906d071a39cbdf6e99a8`.

**Evidence of implemented scope.** The package README says its main design
function is `processVCF`, and sizes the output as barcodes per allele × SNPs ×
two reference/alternate alleles (`README.md:39-48`). The Rd page says:

> “`processVCF` takes a VCF of SNPs (preferably from dbSNP) and turns them into
> a set of labeled MPRA sequences barcoded with inert n-mers”

and documents the output multiplicity as “per allele per SNP”
(`man/processVCF.Rd:28-37,78-87`). The only exported design-facing functions are
`processVCF` and `spread_and_fix_indels` (`NAMESPACE:1-4`). The internal helper
is expressly “process an individual SNP”; it obtains that SNP's genomic context
and concatenates the parts (`man/processSnp.Rd:3-50`).

The paper similarly describes constructs containing the reference and alternate
alleles “along with surrounding genomic sequence,” and says users upload a VCF
to obtain construct sequences (Ghazi et al. 2018, pp. 2682-2683). The source's
`paste0`/concatenation constructs combine fixed primer, context, allele,
restriction, and barcode elements. No second sequence parent or crossover
enumerator enters the API.

**Disconfirmation attempt / absence search.** Across both repositories I
searched every R/Rd/README/workflow/Shiny file and the complete paper for
`recombin*`, `chimer*`, `crossover`, `breakpoint`, `mosaic`, `parent sequence`,
and segment/swap phrases; there were no hits. I then searched for adjacent
`join`, `merge`, `concat`, `paste`, `splice`, `assemble`, `fragment`, `overlap`,
`allele`, and `construct` code and read both `processVCF` implementations, the
package's `processSnp` and `processVCF` documentation, exports, and README. A
different score required a supported operation grouping two parent sequences
and choosing/enumerating crossover positions. I searched for precisely that,
including generic product/combination mechanisms, and found one-variant genomic
context construction only.

### Oligopool Calculator — no (high confidence)

**Evidence of implemented scope.** The current README lists design modules for
primers, barcodes, motifs/anchors, and spacers; assembly modules split/pad long
constructs; degenerate mode compresses/expands supplied variants
(`README.md:20-22,34-42`, commit
`b88fa394ca67ed4c48ec41127b5707694ee7cf0a`). Two operations initially look
close but fail the operational test for different reasons.

First, `join` is described with recombination vocabulary, but it joins **table
columns**, not biological parent segments:

> “Join two oligo tables on `ID` and reconcile branch outputs into one design
> table.”

(`docs/api.md:526-531`). Both inputs must have exactly the same ID set; shared
columns must be identical; only new columns are inserted; and, decisively,
“`join` never creates or drops rows” (`docs/api.md:548-567`). The source copies
the backbone DataFrame and inserts each new column from the other table
(`oligopool/join.py:260-309`). This reconciles two processing branches over the
same members; it does not create chimeric member sequences or enumerate
breakpoints.

Second, `merge` merely “Concatenate[s] contiguous columns into a single column”
(`docs/api.md:423-460`), while `split` fragments each long input oligo into
overlapping pieces for later assembly (`README.md:39-40`; paper, Hossain et al.
2024, pp. 4219-4220). The latter can choose split coordinates, but it does so
inside one existing sequence and returns pieces of that one construct. It does
not join segments from two or more parents.

**Disconfirmation attempt / absence search.** I searched the full package,
docs, changelog, examples, and paper for `recombin*`, `chimer*`, `crossover`,
`breakpoint`, `mosaic`, `parent sequence`, and segment/swap phrases. The only
recombination-word hits were `join`'s parallel-table-branch wording. I then read
the complete `join` API/source, `merge`, `split`, `final`, the YAML branch-join
example, assembly-parser example, degenerate-mode example, README workflow, and
paper assembly description, while searching `join`, `merge`, `split`,
`concatenate`, `fragment`, `overlap`, and `assembly`. A different score required
an operation whose two inputs are sequence parents and whose outputs cross from
one parent to another at generated breakpoints. I specifically checked every
similarly named path and found table reconciliation, fixed column
concatenation, and single-construct fragmentation instead.

### Mutation Maker — no (high confidence)

**Evidence of implemented scope.** Mutation Maker supports SSSM, QCLM/MSDM,
and PAS. Each API input schema has exactly one `gene_of_interest`:
`SSMSequences.gene_of_interest` (`backend/mutation_maker/ssm_types.py:96-102`),
`QCLMSequences.gene_of_interest` (`qclm_types.py:82-93`), and
`PASSequences.gene_of_interest` (`pas_types.py:72-102`), commit
`396c1c0ede7529f3dedf65e821e8c1f20c9a7043`. Their job inputs wrap that one
sequence object plus mutations and configuration (`ssm_types.py:211-218`,
`qclm_types.py:96-107`, `pas_types.py:199-210`). There is no list of parent
sequences or crossover parameter.

PAS is the strongest apparent counterexample, because it chooses fragment
boundaries. The author paper states that its algorithm

> “splits a given gene sequence into an even number of overlapping DNA
> fragments”

and that it optimizes fragment positions until “the entire gene sequence is
covered” (Hiraga et al. 2021, pp. 362-363). Those fragments can contain mutation
combinations and are assembled to full-length genes, but all are derived from
one given gene. This is single-parent gene synthesis/assembly, not joining
segments from two parent sequences. QCLM's combinatorial mutation library is
likewise row-5 mutation combination, not crossover.

**Disconfirmation attempt / absence search.** I searched backend source,
frontend source, README, tests/features, and the complete paper for
`recombin*`, `chimer*`, `crossover`, `breakpoint`, `mosaic`, `parent sequence`,
and segment/swap phrases. After excluding organism/codon-usage data entries
whose species names contain “mosaic,” “chimeric,” or “recombinant,” no design
operation hit remained. I searched and read adjacent `SSSM`, `QCLM`, `PAS`,
`assembly`, `fragment`, `overlap`, `primer`, `gene sequence`, and input-schema
paths. A different score required an input accepting multiple gene parents and
a tool-generated choice of crossover coordinates among them. The full schemas,
workflow source, README, and paper were searched for that mechanism and expose
only one gene of interest.

### DNA Chisel — no (high confidence)

**Evidence of implemented scope.** DNA Chisel documents a problem as one DNA
sequence plus constraints/objectives. The README example constructs
`DnaOptimizationProblem(sequence=...)`, solves it, and retrieves the singular
`problem.sequence` (`README.rst:19-31,44-82`, commit
`68c09304341c3656f3dfe63eda37757d6a7b3917`). The paper is equally explicit:

> “An optimization problem is defined in DNA Chisel by a list of global or
> local specifications against which a starting linear or circular sequence
> will be optimized.”

(Zulkower and Rosser 2020, §2). Exhaustive local variant search therefore edits
one starting sequence; it is not enumeration of chimeras between parents.

I also checked `EnforceRegionsCompatibility`, the closest-sounding built-in.
It takes locations and a caller-provided compatibility condition, evaluates all
pairs of locations **within `problem.sequence`**, and reports incompatible
locations (`dnachisel/builtin_specifications/EnforceRegionsCompatibility.py:8-50`).
It neither takes parent sequences nor constructs joined-segment variants.
`EnforceChoice` similarly lets a region match one of caller-supplied
same-length alternatives; that is a local mutation-space choice, not a
breakpoint generator between parents.

**Disconfirmation attempt / absence search.** I searched all source, docs,
examples, tests, README, and the paper for `recombin*`, `chimer*`, `crossover`,
`breakpoint`, `mosaic`, `parent sequence`, and segment/swap phrases. All
`recombination` hits were biological annotations in the large example E. coli
GenBank record, not code or documentation. I separately searched built-in
specifications and core optimization code for `join`, `merge`, `concat`,
`splice`, `assemble`, `fragment`, `overlap`, `region`, `choice`, and `replace`,
and read `EnforceRegionsCompatibility`, `EnforceChoice`, the README/API index,
and paper. A different score required a shipped specification or operation
accepting at least two parent sequences and producing tool-selected crossover
constructs. That evidence was explicitly sought and is absent.

### tangermeme — no (high confidence)

**Evidence of implemented scope.** The documented atomic sequence operations
are insertions/substitutions of caller-supplied motifs. `insert` accepts `X`, a
`motif`, and one optional `start`; it inserts that motif at the defined position
or the center (`tangermeme/ersatz.py:20-76`, commit
`2006b310cd72a28c56c3ea4ba67f738fff74bb89`). `substitute` has the same
caller-supplied motif/start contract (`ersatz.py:97-166`). Although a caller
could manually slice a second sequence and supply the slice as `motif`, that
would be exactly the global rule's prohibited user reconstruction: tangermeme
does not accept two parents or generate breakpoint combinations.

`multisubstitute` also does not rescue the row. It takes a list of motifs and
caller-provided spacing, then repeatedly calls `substitute`
(`ersatz.py:198-242,266-297`). The README describes `space` as measuring
interactions among motifs at supplied spacings (`README.md:210-220`). Under row
5's binding boundary, that is combinatorial multi-motif placement, not
recombination.

The newer model-guided `beam_substitution` tiles every supplied motif at every
allowed position and retains high-scoring complete sequences
(`tangermeme/design/beam_substitution.py:25-67,82-118`). It is objective-driven
motif substitution from one base sequence, not crossover between parent
sequences; row 5 assigns such output to the model-guided/motif-placement
boundaries rather than this row. The author paper describes the toolkit in the
same vocabulary—sequence manipulations, motif substitution/marginalization,
variant effects, and design—and never documents a recombination operation.

**Disconfirmation attempt / absence search.** I searched all package source,
docs, bundled skill documentation, tests, README, and the full paper for
`recombin*`, `chimer*`, `crossover`, `breakpoint`, `mosaic`, `parent sequence`,
and segment/swap phrases; no hits were returned. I then searched and read every
adjacent API: `insert`, `substitute`, `multisubstitute`, `delete`, `space`,
`product`, variant-effect functions, and all design modules (`screen`, greedy
and beam substitution/marginalization). A different score required a documented
operation taking two or more sequence parents and extracting/combining their
segments across generated breakpoint choices. I specifically checked whether
motif tensors, batch broadcasting, Cartesian-product helpers, or model-guided
design supplied that behavior; they do not.

## Search ledger

The following records every evidence-search pass used for this row.

1. **Instrument.** Read `revision/tables/ROW_DEFINITIONS.md:1-138`, including
   the preamble, global rule, complete row 5 boundary rules, and row 10.
2. **Lead-only locator pass.** Searched each of the eight
   `tool_survey/final/<slug>.md` files for
   `recombin|chimera|crossover|breakpoint|segment|swap|join|parent`; these files
   were not cited and their ratings/conclusions were disregarded.
3. **Repository inventory.** Inventoried local PoolParty docs/source/tests and
   all Markdown/RST/Python/R/Rd/HTML/YAML files in fresh read-only audit clones
   of the seven cited public repositories. Recorded the exact commits shown in
   the tool sections/table above; also inspected `designMPRA` alongside the
   package repository as requested.
4. **Uniform direct-capability grep.** Across every repository, ran the same
   case-insensitive family:
   `recombin|chim(ae|e)r|cross[- ]?over|breakpoint|mosaic|parent(al)? (sequence|allele)|segment.{0,30}swap|swap.{0,30}segment`.
   Data-only GenBank, species-name, and codon-usage hits were inspected and
   classified explicitly above rather than silently discarded.
5. **Uniform adjacent-operation grep.** Across every repository, searched
   `join|merge|concat|splice|assemble|assembly|substitute|insert|replace|segment|fragment|overlap`,
   then read the operation/API paths that could plausibly implement the test.
   Per-tool additions were: VaLiAnT mutation/targeton/custom/background paths;
   MPRAnator's four modules, SNP combinations, and construct ordering; MPRA
   Design Tools' VCF/SNP/allele/construct paths; Oligopool's `join`, `merge`,
   `split`, `final`, branch DAG, assembly and degenerate examples; Mutation
   Maker's SSSM/QCLM/PAS schemas and workflows; DNA Chisel's problem,
   region/choice and specification paths; tangermeme's ersatz, space, product,
   variant-effect, and design modules.
6. **Papers.** Extracted the seven non-PoolParty author PDFs with PyMuPDF as
   instructed and searched the full text with the same direct family plus each
   tool's adjacent vocabulary (mutator/reference/targeton;
   motif/SNP/PWM/Transmutation; VCF/allele/context; split/assembly/fragment;
   SSSM/QCLM/PAS; constraint/objective/single sequence; and
   motif/substitution/variant/design). PoolParty's repository was treated as
   decisive primary evidence, and its author PDF was not used as evidence, as
   required.
7. **PoolParty behavior.** Read the complete documented operation and source
   paths, materialized a two-parent sequential recombination example with the
   specified repository venv and `PYTHONDONTWRITEBYTECODE=1`, and ran the full
   targeted `tests/test_recombine.py` module (44/44 passed).

## Row-level finding

The row is applicable on one consistent scale and does discriminate, though
narrowly: PoolParty is the sole positive (1 `yes`, 7 `no`). The narrowness is
substantive rather than a search artifact. Several tools contain tempting
near-neighbours—MPRAnator motif/SNP combinations, Mutation Maker and Oligopool
assembly fragments, Oligopool's table-branch `join`, DNA Chisel/tangermeme edit
spaces—but the binding requirement is specifically two or more **parent
sequences** plus **tool-generated crossover breakpoints**. Applying that same
threshold prevents both over-crediting generic assembly/placement and
under-crediting PoolParty's explicit chimeragenesis operation.
