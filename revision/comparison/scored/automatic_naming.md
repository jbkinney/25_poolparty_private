# Audit: `automatic_naming`

Audit date: 2026-08-15

## Measurement instrument and operational test

I read the preamble and binding global rule in
`25_poolparty_private/revision/tables/ROW_DEFINITIONS.md:1-34`, and the
`automatic_naming` definition at lines 163-167, before inspecting any tool. The
binding points are that the capability must be a tool-provided operation,
parameter, or mode; caller reconstruction does not count; and a user-supplied ID
carried through the output does not count. The row explicitly assigns `partial`
when the tool generates names that do not encode construction history, or does
so for only some operations.

**Operational test declared before tool inspection (with “identifier” made
explicit here):** For each tool, inspect the actual output of every documented
design/export operation and ask whether the tool itself assigns each emitted
sequence an identifier whose components encode that sequence's construction
path; score `yes` when it does across the tool's design operations, `partial`
when the tool emits an identifier/header but it is generic/indexed or naming is
limited to some operations, and `no` when the only identity is caller-supplied,
an incidental container index, absent, or requires caller-written naming logic.

An exported field explicitly presented as `Name`, `ID`/index label, or a FASTA
header is an identifier for this test. A bare DataFrame/Tensor row number is not.
This is why a tool-generated numeric FASTA header can be `partial`, while an
unnamed tensor cannot.

## Result

| Tool | Value | Confidence |
|---|---|---|
| poolparty | **yes** | high |
| valiant | **yes** | high |
| mpranator | **yes** | high |
| mpradesign | **partial** | medium |
| oligopoolcalc | **partial** | high |
| mutationmaker | **yes** | medium |
| dnachisel | **no** | high |
| tangermeme | **partial** | high |

## Source snapshots

Repository evidence was examined at these commits:

- PoolParty local repository: `1bb0179e1c3720b1fffd471802b3040f9336de28`
- VaLiAnT `develop`: `8796cc112dafd4919fec59913f58cd2be87c45eb`
- MPRAnator: `9969790d62410138d4281b7955da6d085f07b1bc`
- mpradesigntools: `afd386ef12051bb0a37ad63a6f92acd555246584`
- designMPRA: `0cf56ee602fc86dde705906d071a39cbdf6e99a8`
- Oligopool Calculator: `b88fa394ca67ed4c48ec41127b5707694ee7cf0a`
- Mutation Maker: `396c1c0ede7529f3dedf65e821e8c1f20c9a7043`
- DnaChisel: `68c09304341c3656f3dfe63eda37757d6a7b3917`
- tangermeme: `2006b310cd72a28c56c3ea4ba67f738fff74bb89`

The per-tool survey records were not used as evidence or as assessments.

## Per-tool audit

### poolparty — `yes` (high confidence)

**Evidence.** The dedicated naming documentation says, “Names let you trace
exactly how a sequence was constructed”
(`poolparty-statecounter/poolparty/docs/metadata/naming.rst:4-6`). It documents
operation-contributed segments joined from source to downstream at lines 20-29,
stable variant indices at lines 57-62, and chained names such as `bg.mut_0` at
lines 119-137. Scan operations add distinct operation, position, and inserted-
variant segments; the example is `scan_00.pos_0.base_0`
(`docs/metadata/naming.rst:190-205`).

The implementation collects each operation's contributions
(`src/poolparty/generate_library.py:294-295`) and assigns
`row["name"] = final_name` after joining them in topological order (lines
328-330). `Operation.compute_name_contributions` generates fixed, state-indexed,
or global-state-indexed segments from the operation prefix and state
(`src/poolparty/operation.py:413-449`). The caller supplies a meaningful prefix,
but not the final per-sequence ID: PoolParty appends state/variant indices and
composes the path. That is not mere ID carry-through.

**Behavioral confirmation.** Per the task requirement, I ran the documented
chain from `poolparty-statecounter/poolparty` in the read-only PoolParty
environment with bytecode disabled:

```text
PYTHONDONTWRITEBYTECODE=1 ../.venv/bin/python -c 'import poolparty as pp; pp.init(); wt=pp.from_seq("ATCGATCG", prefix="bg"); muts=wt.mutagenize(num_mutations=1, num_states=3, prefix="mut"); print(muts.generate_library()[["name","seq"]].to_string(index=False))'
```

It returned `bg.mut_0`, `bg.mut_1`, and `bg.mut_2`, paired with the three emitted
sequences.

**Search and disconfirmation attempt.** I searched the full Python and
documentation trees with
`rg -n -i --glob '*.py' --glob '*.md' --glob '*.rst' --glob '*.ipynb' '(automatic.*nam|compute_name|name_contribution|seq_name|oligo_name|record[_ ]?id|identifier|FASTA|write.*fasta|name=|prefix)' poolparty-statecounter/poolparty`, then inspected the complete naming page,
the generator's output-row construction, the base operation's naming method,
and ran the chain above. Evidence that would lower this to `partial` would be a
documented design operation whose emitted sequences cannot receive the composed
name, or names consisting only of a tool-generated ordinal; evidence for `no`
would be only a carried caller ID. I specifically searched for those cases and
found instead the common contribution collector and compound scan naming.

### valiant — `yes` (high confidence)

**Evidence.** VaLiAnT explicitly constructs oligo names. For SGE designs,
`get_sge_oligo_name` combines transcript/gene, contig, a generated variant
fragment, source, and reverse-complement suffix
(`src/valiant/meta_table.py:91-108`); cDNA naming similarly combines sequence,
transcript/gene, variant fragment, and source at lines 111-115. The variant
fragment encodes position/range and the deletion, substitution (`ref>alt`), or
insertion (`alt`) event (`src/valiant/variant.py:141-153`). The name is generated
from the actual protected/altered variant at `meta_table.py:512-514` and written
as `oligo_name` at lines 629-634.

The README calls `oligo_name` the “Name of the oligonucleotide”
(`README.md:737-750`) and shows history-rich examples such as
`...chr3:10146513_G>A_snv` in the unique-oligo output (`README.md:805-815`). When
multiple construction paths yield the same sequence, the unique output chooses
one generated name deterministically after sorting (`meta_table.py:676-685`).

**Search and disconfirmation attempt.** I ran `find` over the repository and
`rg -n -i --hidden --glob '*.py' --glob '*.md' --glob '*.rst' --glob '*.csv' '(oligo.?name|identifier|sequence.?name|record.?id|get_.*name|name.*oligo|fasta|metadata)' /tmp/automatic-naming-valiant`, then inspected the name constructors, variant
formatter, metadata write path, unique-sequence write path, README output schema,
and checked-in expected-output CSVs. Evidence for `partial` would be generic
generated names or a design mode without construction-derived names; evidence
for `no` would be names copied solely from the input VCF ID. I searched the
custom-VCF fields as well: `vcf_var_id` is separate metadata, while `oligo_name`
is constructed independently. No contrary mode was found.

### mpranator — `yes` (high confidence)

**Evidence.** In Motif Design, the source sorts motif/position pairs, formats
each as `motif-position`, joins them, and creates the FASTA header
`> Background-{background}|{motif-position...}` (`oligo.py:35-76`). The final
assembly adds barcode replicate, restriction-site, adapter, and duplicate-site
information to that header (`part1.py:210-291`). SNP Design constructs headers
containing the input sequence identity plus SNP ID, position, nucleotide change,
edge trimming/padding, or `REFERENCE` (`parseSNPs.py:160-189`). Transmutation
appends mutation count and scramble/reverse/complement flags to each output
header (`iliasApp/views.py:175-211`).

For the fourth module, the authors' supplement states, “The header contains
information regarding the probability ... or the simulation number” for PWM
Seq-Gen (`btw584_supp.docx`, §2.4, obtained from the official Europe PMC
supplementary-files endpoint for PMC5198521). Sections 2.1-2.3 of the same
supplement independently describe the motif/subcomponent, SNP/change, and
transmutation information stored in output headers. Thus every documented
module generates headers from what it did to that sequence; the inherited
scaffold name is only one component.

**Search and disconfirmation attempt.** I enumerated the repository with
`find /tmp/automatic-naming-mpranator -type f` and searched it using
`rg -n -i --glob '*.py' --glob '*.md' '(name|identifier|record.?id|fasta|header|barcode|write|output)' /tmp/automatic-naming-mpranator`, followed by targeted searches for
`Prob`, `PWM`, `Seq.Gen`, and `header`. I inspected `oligo.py`, `part1.py`,
`parseSNPs.py`, `iliasApp/views.py`, and the README. Because the repository does
not contain the PWM Seq-Gen implementation, I downloaded the authors' primary
supplement from
`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5198521/supplementaryFiles`,
extracted `btw584_supp.docx`, converted `word/document.xml` to text, and searched
it with `rg -n -i -C 3 '(header|probability|simulation|oligo name|output file|FASTA)'`.
Evidence for `partial` would be one module emitting only generic or inherited
headers; evidence for `no` would be caller-only headers throughout. The
four-module check found neither.

### mpradesign — `partial` (medium confidence)

This score covers the `mpradesigntools` package and its `designMPRA` web wrapper.

**Evidence.** The documented return is a data frame of “labeled MPRA sequences”
(`man/processVCF.Rd:78-87`). For each SNP, `processSnp` copies the caller's
`ID = snp$ID`, but also generates `snpIndex`, ref/alt `type`, allele, and sampled
barcodes (`R/processVCFfast.R:389-399`). The completed output adds a unique
`totIndex = 1:nrow(.)` and exports the fields `ID`, `type`, `allele`, `snpIndex`,
`totIndex`, `barcode`, and `sequence` (`R/processVCFfast.R:1243-1258`). The web
wrapper writes that table unchanged (`designMPRA/server.R:161-165`).

`ID` does not count because it is carried from the VCF. `totIndex` does count as
a tool-generated per-output label, corroborated by the API's “labeled” wording,
but it is only an ordinal and does not compose the adjacent construction facts
(`type`, `allele`, barcode, repair notes) into an identifier. This is precisely
the row's `partial` case: a generated name/identifier not composed from history.
The barcode is also generated per construct, but is a sequence component and is
not presented as an informative construction name.

**Search and disconfirmation attempt.** I enumerated and searched both primary
repositories with
`rg -n -i --glob '*.R' --glob '*.Rd' --glob '*.md' --glob '*.html' '(name|identifier|fasta|header|record.?id|write|output|barcode|totIndex|snpIndex)' /tmp/automatic-naming-mpradesigntools /tmp/automatic-naming-designmpra`.
I inspected `processVCFfast.R`, the `processVCF` and `processSnp` manuals, README,
tests/examples containing output tables, and the Shiny server's download path.
Evidence for `yes` would be a tool-built identifier combining, for example, SNP,
allele, construct index, and/or barcode for every sequence. Evidence for `no`
would be absence of any output identifier beyond the carried VCF ID or an
incidental DataFrame index. I searched for both: no composed name exists, but
the selected/exported `totIndex` field and “labeled” API contract rule out `no`.

### oligopoolcalc — `partial` (high confidence)

**Evidence.** The normal data contract requires the caller to supply a unique
`ID` primary key and says output CSVs include that explicit ID
(`docs/agent-skills.md:82-90`). Input validation resets the DataFrame onto that
`ID`, verifies uniqueness, and preserves it as a string
(`oligopool/base/validation_parsing.py:870-915`). That carry-through does not
count.

One documented operation does generate sequence IDs: `compress` returns a
synthesis table with `DegenerateID` and `DegenerateSeq`
(`oligopool/compress.py:39-54`). Its implementation numbers outputs as
`deg_id = 'D{}'.format(degenerate_idx)` and assigns that ID to every emitted
degenerate oligo (`oligopool/base/core_degenerate.py:1041-1067`). These are
generic generated IDs (`D1`, `D2`, ...), not construction-history names, and the
behavior is limited to compression. Either restriction independently warrants
`partial`.

**Search and disconfirmation attempt.** I searched the full source and docs with
`rg -n -i --glob '*.py' --glob '*.md' --glob '*.rst' '(required.*ID|ID.*required|ID column|preserv|id_from_index|generate.*ID|assign.*ID|DegenerateID|Fragment|identifier|fasta|header)' /tmp/automatic-naming-oligopool` and inspected the core data contract, validation/parser,
all public operation return schemas, `compress.py`, and the degenerate synthesis
builder. Evidence for `yes` would be construction-derived IDs generated across
the normal design operations. Evidence for `no` would be only carried input IDs.
I confirmed carry-through for ordinary operations but found the explicit `D#`
generation in `compress`, so neither `yes` nor `no` fits.

### mutationmaker — `yes` (medium confidence)

**Evidence.** The product's export tables explicitly contain a `Name` column for
all three workflows (`frontend/src/shared/components/SaveFile/index.tsx:54-57`).
The export code generates names as follows:

- SSM: user prefix plus plate number, forward/reverse construction direction,
  and well position (`SaveFile/index.tsx:92-149`), e.g. the expression at lines
  136 and 142.
- QCLM: the tool first groups outputs by mutation source and position, then names
  each site/alternative with a generated site ordinal and optional letter
  (`SaveFile/index.tsx:152-195`).
- PAS: the tool names each emitted oligo from its generated fragment index and
  within-fragment oligo index (`SaveFile/index.tsx:198-310`), specifically
  `${oligoPrefix}-Fr${index + 1}-${oIndex + 1}` at line 301.

The user can choose the prefix (for example,
`frontend/src/scenes/SSM/components/SSMForm.tsx:348-352`), but does not provide
the final identifiers. Mutation Maker composes workflow-specific construction
coordinates onto it and emits one `Name` per export row. Plate/direction/well,
site/alternative, and fragment/oligo are construction history under the same
threshold used for PoolParty's operation/state indices.

**Search and disconfirmation attempt.** I searched the full TypeScript, Python,
and documentation trees with
`rg -n -i --glob '*.py' --glob '*.ts' --glob '*.tsx' --glob '*.md' '(automatic.*nam|oligo.?name|oligo_prefix|identifier|fasta|header|export|SaveFile|prefix)' /tmp/automatic-naming-mutationmaker`, then inspected the shared export component,
all three form/result flows, API types, and backend config types. Evidence for
`partial` would be a workflow with only a generic ordinal or no generated suffix;
evidence for `no` would be the prefix copied unchanged. I searched each SSM,
QCLM, and PAS export branch and found generated, workflow-specific suffixes in
all three. Confidence is medium rather than high because the identifiers encode
construction organization rather than spelling out the full mutation, although
the separate export columns retain that mutation detail.

### dnachisel — `no` (high confidence)

**Evidence.** The core output method has a caller parameter `record_id=None` and
sets `record.id = record_id` only when supplied
(`dnachisel/DnaOptimizationProblem/mixins/RecordRepresentationMixin.py:76-155`).
Without it, the helper defaults to `id="<unknown id>"` and
`name="<unknown name>"` (`dnachisel/biotools/genbank_operations.py:177-190`).
The documented genome-wide example performs the bookkeeping itself: it extracts
the input feature's `protein_id`, calls `problem.to_record(record_id=protein_id)`,
collects the records, and writes FASTA
(`examples/common_scenarios/genome-wide-optimization.py:9-32`). This is carried
input identity, not a name generated from the optimization history.

The only other apparent hit, `_sequences_to_new_records`, accepts caller-provided
`('name', 'sequence')`, dictionaries, or already named records and copies that
name into `id`/`name`
(`dnachisel/reports/constraints_reports/constraints_reports.py:38-55`). It is
therefore not counterevidence.

**Search and disconfirmation attempt.** This `no` is based on an explicit
absence audit. I searched all `dnachisel/`, `docs/`, and `examples/` Python/RST/
Markdown/notebook files with
`rg -n -i --glob '*.py' --glob '*.rst' --glob '*.md' --glob '*.ipynb' '(automatic.?nam|construction.?history|record_id|record\.id|sequence_to_biopython_record|fasta|identifier)'`, then separately searched
`rg -n -i 'record_id\s*=|record\.id\s*=|sequence_to_biopython_record'` over the
same scopes. I inspected `RecordRepresentationMixin.to_record`, GenBank/FASTA
load/write helpers, optimization and constraint reports, CLI, common-scenario
examples, manuscript examples, and README/docs hits. Evidence for `partial`
would be any tool-generated generic sequence ID; evidence for `yes` would be a
record ID derived from applied constraints/objectives/edits. I searched for both
and found only unknown defaults, caller arguments, copied record IDs, and
caller-named report inputs. Thus `no` means inspected-and-absent, not not-found.

### tangermeme — `partial` (high confidence)

**Evidence.** tangermeme's design routines return tensors, not named sequence
records: `greedy_substitution` returns an edited `torch.Tensor`
(`tangermeme/design/greedy_substitution.py:146-159,249`), `beam_substitution`
returns ranked designed sequences as a tensor
(`tangermeme/design/beam_substitution.py:173-186`), and
`greedy_marginalize` returns a tensor construct
(`tangermeme/design/greedy_marginalize.py:144-148`). The README likewise says
“all data are PyTorch tensors” (`README.md:80-86`). A tensor's position alone is
not an output identifier.

However, the tool provides a documented matching export interface:
`one_hot_to_fasta` accepts a tensor and says that when caller headers are absent,
“the numeric index is used” (`tangermeme/io.py:540-568`). The implementation
writes `> {i}` for each sequence at lines 575-582. Feeding the design routine's
documented tensor output to the tool's own documented tensor-to-FASTA input is
allowed by the global composition rule. The result is therefore a tool-generated
per-sequence name, but a generic ordinal with no edit/motif/position history:
`partial`, not `yes`.

**Search and disconfirmation attempt.** I enumerated the package and searched
source, docs, examples, tests, and README with
`rg -n -i --glob '*.py' --glob '*.md' --glob '*.rst' --glob '*.ipynb' '(name|identifier|record.?id|fasta|header|write|output)' /tmp/automatic-naming-tangermeme`, then searched the design modules and `io.py` specifically for `headers`,
`return`, `greedy`, `beam`, and FASTA conversion. Evidence for `yes` would be a
design/export identifier encoding selected motif, insertion/substitution
position, iteration, or design trajectory. Evidence for `no` would be an
exporter that required caller headers or supplied no identifier. I confirmed
that design outputs have no such history metadata, but also confirmed the
exporter's automatic numeric headers; hence `partial`.

## Row-level finding

The row does discriminate on one consistent scale in this set: four tools build
identifiers from construction information, three expose restricted or generic
tool-generated identity, and one exposes no tool-generated sequence identity.
The `partial` bucket is deliberately heterogeneous because the supplied
instrument defines both generic names and operation-limited generation as
`partial`. The closest boundary cases are mpradesign (an explicit exported
`totIndex`) and tangermeme (an explicit numeric FASTA header). Both receive the
same value because both are tool-emitted per-sequence identifiers lacking
construction history; incidental DataFrame/Tensor positions alone were not
credited.
