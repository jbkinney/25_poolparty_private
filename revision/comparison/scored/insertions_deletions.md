# Insertions and deletions — independent capability audit

Date of audit: 2026-08-15

## Instrument and operational test

I read both the preamble/global rule and row 8 in
`25_poolparty_private/revision/tables/ROW_DEFINITIONS.md` before examining any
tool. The binding parts are:

> “A capability counts only where the tool provides an operation, parameter or
> mode for it.” (lines 16–18)

> “Generation of variants that insert or remove bases as a **designed variant
> type**, via a supported operation. Supplying pre-made indel sequences as input
> does not count.” (lines 127–129)

> “`partial` = one of insertion or deletion only, or fixed length only.”
> (lines 130–131)

The global rule also says that a restricted tool-provided capability is
`partial` (lines 20–24).

**Operational test, fixed before inspecting the tools:** Starting from a parent
sequence and an indel design request, can a documented tool-provided operation
generate output variants containing both inserted and deleted bases (`yes`),
only one event class or a restricted/fixed event length (`partial`), or neither
(`no`), with pre-made mutant sequences and caller-written composition excluded?

For this test, an event-level specification such as VCF `POS/REF/ALT` is not a
“pre-made indel sequence”: the tool still has to retrieve/use the parent,
apply the event, and construct the output sequence. A synthesis-ready mutant
FASTA/CSV row supplied as input *is* a pre-made indel sequence and does not
count. This is the only interpretation that treats an event-level API such as
`insert(parent, bases, position)` and the equivalent event-level VCF input
consistently.

## Results

| Tool | Value | Confidence |
|---|---|---|
| PoolParty | **yes** | high |
| VaLiAnT | **yes** | high |
| MPRAnator | **partial** | high |
| MPRA Design Tools | **yes** | high |
| Oligopool Calculator | **no** | high |
| Mutation Maker | **no** | high |
| DNA Chisel | **no** | high |
| tangermeme | **yes** | high |

Repository snapshots used (the records under `tool_survey/final/` were used
only as location leads and are not evidence):

- PoolParty: `jbkinney/poolparty-statetracker`, local read-only checkout at
  `1bb0179e1c3720b1fffd471802b3040f9336de28`.
- VaLiAnT: `cancerit/VaLiAnT`, `develop`,
  `8796cc112dafd4919fec59913f58cd2be87c45eb`.
- MPRAnator: `hemberg-lab/MPRAnator`, `master`,
  `9969790d62410138d4281b7955da6d085f07b1bc`.
- MPRA Design Tools: `andrewGhazi/mpradesigntools`,
  `afd386ef12051bb0a37ad63a6f92acd555246584`, and
  `andrewGhazi/designMPRA`, `0cf56ee602fc86dde705906d071a39cbdf6e99a8`.
- Oligopool Calculator: `ayaanhossain/oligopool`,
  `b88fa394ca67ed4c48ec41127b5707694ee7cf0a`.
- Mutation Maker: `ra100/Mutation_Maker`,
  `396c1c0ede7529f3dedf65e821e8c1f20c9a7043`.
- DNA Chisel: `Edinburgh-Genome-Foundry/DnaChisel`,
  `68c09304341c3656f3dfe63eda37757d6a7b3917`.
- tangermeme: `jmschrei/tangermeme`,
  `8d732b8c08764057b7ae5faa644d48664f36b44b`.

## Per-tool evidence and disconfirmation

### PoolParty — **yes** (high confidence)

**Insertion evidence.** The operation documentation says:

> “Insert sequences from `insertion_pool` at every position along the
> background sequence … no background bases are removed, so output sequences
> are longer than the input.”

Source: `poolparty-statecounter/poolparty/docs/operations/insertion_scan.rst:4-6`.
The `insertion_pool` may be a `Pool | str`, and `positions=None` means all valid
positions (`insertion_scan.rst:34-47`). The implementation sets a zero-width
marker in insertion mode and passes the insertion pool to `replace_region`
(`src/poolparty/scan_ops/insertion_scan.py:96-130`).

**Deletion evidence.** The documentation says:

> “Slide a deletion window of fixed length across the sequence … and, at each
> position, remove those bases.”

and

> “Pass `deletion_marker=None` to produce shorter sequences instead.”

Source: `poolparty-statecounter/poolparty/docs/operations/deletion_scan.rst:4-8`.
`deletion_length` is a caller-selected positive integer and `positions=None`
means all valid starts (`deletion_scan.rst:32-50`). In the implementation,
`deletion_marker=None` becomes `marker_content = ""`, which is passed to
`replace_region` (`src/poolparty/scan_ops/deletion_scan.py:74-120`). Thus the
operation supports multiple deletion lengths rather than one fixed length.

**Required behavioral check.** I ran the local read-only venv with
`PYTHONDONTWRITEBYTECODE=1` and fresh `pp.init()` calls. From parent `ACGT`:

```text
ins2 seq_length= 6 num_states= 5 seqs= ['GGACGT', 'AGGCGT', 'ACGGGT', 'ACGGGT', 'ACGTGG']
del1 seq_length= 3 num_states= 4 seqs= ['CGT', 'AGT', 'ACT', 'ACG']
del2 seq_length= 2 num_states= 3 seqs= ['GT', 'AT', 'AC']
```

This directly confirms lengthening insertions and genuinely shorter 1- and
2-base deletions. (The duplicated `ACGGGT` is expected because inserting `GG`
on either side of the existing `G` produces the same sequence.)

**Disconfirmation attempt.** I recursively searched `src/`, `docs/`, and
`README.md` for the case-insensitive biological terms `indel`, `insert`,
`insertion`, `delete`, and `deletion`; then read both operation pages and both
implementations above. I also checked whether insertion was really replacement,
whether deletion only emitted an unchanged sequence, and whether either event
length was hard-coded. The docs explicitly distinguish insertion from
replacement, and live execution showed length changes and multiple deletion
lengths. Evidence that would change the value: only one of the two operation
classes, an invariant event length, a gap that neither shortened nor marked the
output, or only pre-constructed mutant inputs. I searched for each and found the
opposite.

### VaLiAnT — **yes** (high confidence)

**Designed deletion evidence.** The documented mutation types include
parametric and in-frame deletions (`README.md:255-270`). The parametric deletion
documentation says:

> “Non-overlapping stretches of nucleotides of a given length are deleted
> starting from a given offset.”

Source: `README.md:272-284`. The accepted form is `<SPAN>del[<OFFSET>]`, so the
length is not fixed. In source, `DeletionMutator.get_variants` obtains every
reference window and returns `Variant.get_del(...)`
(`src/valiant/mutators/deletion.py:33-48`); the mutator parser accepts
`^(\d+)del(\d*)$` and constructs the requested span and offset
(`src/valiant/loaders/mutator_config.py:30-57`).

**Insertion and deletion from event specifications.** The custom-variant
documentation says:

> “Applied to the targeton reference sequence as a whole. Only simple variants
> such as the following are supported: substitutions; insertions; deletions;
> indels.”

Source: `README.md:465-476`. These are event-level VCF records, not pre-made
oligo sequences: classification is based on VCF `POS`, `REF`, and `ALT`
(`README.md:474-476`). `CustomVariant.from_record_with_id` normalizes the VCF
event and classifies a positive ALT–REF length delta as `INSERTION` and a
negative delta as `DELETION` (`src/valiant/custom_variant.py:79-123`). The
actual library construction is also explicit: `OligoSeq.from_ref` calls
`alter_seq`, which calls `seq.alter(variant.ref_range, variant.is_insertion,
variant.alt)` (`src/valiant/oligo_seq.py:30-46`); `Seq.alter` uses
`insert_substr` for insertion and `replace_substr` otherwise
(`src/valiant/seq.py:79-91`).

**Disconfirmation attempt.** I recursively searched the entire `develop`
checkout (source, README, changelog, tests, and examples; binary assets and
literal VCF/FASTA data excluded from the first pass) for `indel(s)`,
`insert(ion)(s)`, and `delet(ion)(s)`. I then read the complete documented
mutation-type section, custom-variant section, `mutator_type.py`,
`mutator.py`, `mutators/deletion.py`, `loaders/mutator_config.py`,
`custom_variant.py`, `variant.py`, `seq.py`, `oligo_seq.py`, and the targeton
custom-variant path. This attacked the alternatives that insertion might be
metadata-only or that VCF events might merely be passed through. The code
constructs altered oligo sequences. Evidence that would change the value:
deletion alone would be `partial`; VCF events passed through without sequence
construction, or only already-mutated input sequences, would not count. Both
were searched for and contradicted by the cited implementation.

### MPRAnator — **partial** (high confidence)

MPRAnator provides both events but restricts them to short indels.

**Documentation evidence.** The official paper supplement, §2.2, states:

> “Deletions or insertions (up to 10nt) in VCF format are also allowed.”

It further says that inserted sequences are trimmed and deleted sequences are
expanded with adenines to keep oligo lengths similar. Source: official
supplement `btw584_supp.docx`, §2.2 (local primary-source cache
`/tmp/mpranator-supp/btw584_supp.docx`, `word/document.xml`); the Fig. 2 legend
likewise says variants are generated, deletions receive edge adenines, and
insertions are trimmed. This explicitly supports both local event types while
preserving overall oligo length.

**Implementation evidence.** The core function documents its inputs:

> “Inputs of users are SNPs in VCF format and FASTA file with sequences of
> equal size.”

Source: `parseSNPs.py:5-6`. It computes `difference = len(REF) - len(ALT)` and,
for `difference < 0`, inserts the longer ALT and trims the opposite edge; for
`difference > 0`, it skips deleted parent bases and pads an edge with adenines
(`parseSNPs.py:75-123` for combination mode and `:151-186` for single-event
mode). The supported-operation path is the web app itself:
`MpraSnpsForm` labels the input “Enter your SNPs (VCF Format)” and rejects sizes
greater than 10 (`iliasApp/forms.py:501-545`), and `mpraSnpResults` passes the
validated event objects to `parseSNPs.make_sequence_copies`
(`iliasApp/views.py:319-363`).

**Why partial, not yes.** Both event classes are implemented, but the official
documentation states the indel-length cap and the implementation gates on
`abs(difference) < 10` (`parseSNPs.py:82-84`, `:163-164`). This is a stated
scope restriction under the global partial rule. (The form says `> 10` is
invalid while the older core uses `< 10`; this boundary discrepancy does not
affect the score because either source imposes a short-indel cap.)

**Disconfirmation attempt.** I searched every repository file for the same
indel/insertion/deletion terms, read `parseSNPs.py`, the MPRA SNP form/view,
the downloadable API script, app documentation template, the supplied main
paper PDF (text extracted with PyMuPDF), and the official supplementary
material. I specifically searched for an unrestricted or parametric indel mode
that would warrant `yes`, and for evidence that indels were only pre-made full
sequences (which would warrant `no`). No unrestricted mode exists; instead,
the tool takes VCF events and itself constructs the outputs under a 10-nt cap.
Evidence that would change the value: removal of the cap in the supported app
and implementation would produce `yes`; absence of either event or absence of
the event-application operation would produce `partial` for the former or `no`
for the latter. All were searched.

### MPRA Design Tools — **yes** (high confidence)

**Documented operation.** The package calls its principal operation
`processVCF`; the manual page is titled “Process VCF into MPRA sequences” and
says it:

> “takes a VCF of SNPs … and turns them into a set of labeled MPRA sequences
> barcoded with inert n-mers.”

Source: `mpradesigntools/man/processVCF.Rd:3-31,78-87`. The README explicitly
documents both types:

> “Insertions and deletions must encode the reference and alternate alleles
> (respectively) as a dash character '-'.”

Source: `mpradesigntools/README.md:50-59`. This is an event representation
constraint, not a requirement to supply complete mutant oligos.

**Insertion implementation.** `processSnp` classifies `REF == '-'` as insertion
and `ALT == '-'` as deletion (`R/processVCFfast.R:232-242`).
`generateInsConstruct` takes the genomic context and insertion allele and
returns context-left + inserted allele + context-right
(`R/processVCFfast.R:36-60`). The insertion branch retrieves hg38 context,
calls that function, builds separate reference and alternate `constrseq`
values, and then assembles complete barcoded `sequence` values
(`R/processVCFfast.R:544-640`).

**Deletion implementation.** `generateDelConstruct` concatenates the parent
context before the deleted range with the context after it
(`R/processVCFfast.R:63-70`). The deletion branch derives `refwidth`, retrieves
the parent genomic sequence, calls `generateDelConstruct`, emits separate
reference/alternate rows, and assembles complete output sequences
(`R/processVCFfast.R:764-853`). Event lengths come from `nchar(ALT)` or
`nchar(REF)`, not a single fixed constant.

The paper corroborates the operation at a higher level: users “upload VCFs for
automated construct sequence generation” (supplied
`Ghazi2018_MPRADesignTools_all.pdf`, extracted text lines 17–23) and receive
MPRA sequences after providing a VCF (`paper`, Results/sequence-generation
section, extracted lines 116–121).

**Disconfirmation attempt.** I recursively searched both official repositories
(`mpradesigntools` and `designMPRA`) and the supplied paper for the biological
indel/insertion/deletion terms, then read the README VCF constraints, generated
`processVCF` manual page, both construct functions, all three event branches,
and the old Shiny implementation. I tested the adverse interpretation that
the dash-formatted VCF might simply be copied to output; it is not—the code
retrieves hg38 context and constructs altered and complete oligo sequences.
Evidence that would change the value: only one branch, a fixed event width, or
acceptance only of pre-built mutant sequences. The searches found explicit
variable-width insertion and deletion branches instead.

### Oligopool Calculator — **no** (high confidence)

**What the supported API actually provides.** The package's complete public API
list is `barcode`, `primer`, `motif`, `spacer`, `background`, `split`, `pad`,
`merge`, `revcomp`, `join`, `lenstat`, `verify`, `inspect`, `final`, `compress`,
`expand`, `index`, `pack`, `acount`, and `xcount`
(`oligopool/__init__.py:6-31`; the mode grouping is at `:42-69` and
`README.md:117-144`). There is no indel-variant generator.

The superficially positive `motif` and `join` hits do not pass the operational
test. `motif` takes an input DataFrame “with annotated oligo pool variants” and
adds one designed motif column per existing row
(`oligopool/motif.py:17-41,56-75`); it does not fan out indel alternatives from
a parent. `join` joins two existing tables with exactly the same IDs and inserts
*columns* into column order (`oligopool/join.py:13-45,88-95`). Those are
architecture assembly, not insertion/deletion variant generation. The docs
also direct users to “Start with your variants” (`docs/docs.md:106`) and, for
degenerate mode, “Design variants computationally” before passing them to the
tool (`docs/docs.md:893-902`).

The only actual biological `Indels` hits are in read processing, where the code
says:

> “We don't care about Indels, Skip this Read.”

Source: `oligopool/base/core_pack.py:549-557` and `:779-787`.

**Disconfirmation attempt.** I recursively searched the complete checkout,
including `README.md`, `docs/docs.md`, `docs/api.md`, `docs/agent-skills.md`,
all top-level modules, all `base/` modules, and examples for `indel(s)`,
`insert(ion)(s)`, `delete/deletion(s)`, `variant`, `mutation`, `allele`, and
`VCF`. I then read the complete public API registration plus the apparent
positives `motif`, `join`, `compress`, `expand`, and the indel read-processing
code. No operation creates an insertion or deletion alternative from a parent;
existing variants are input rows. Evidence that would change the value: a
public operation taking a parent plus an insertion/deletion event (both would
yield `yes`, one/fixed would yield `partial`) and emitting newly altered rows.
I searched the API, docs, source, and examples for exactly that and found none.

### Mutation Maker — **no** (high confidence)

**Positive evidence for the narrower mutation model.** The paper defines SSSM
as a library in which a residue is “randomly substituted with all possible
amino acids” and says the methods use designated-position “codon
substitutions” (supplied `Hiraga2021_MutationMaker.pdf`, extracted text lines
83–107). The program inputs confirm this. `SSMInput.mutations` is a required
list of strings parsed by `parse_codon_mutation`
(`backend/mutation_maker/ssm_types.py:211-218`), and `QCLMInput` documents
strings such as `E32W E32L E49K`
(`backend/mutation_maker/qclm_types.py:96-107`). In PAS output generation,
every selected mutation calls `Codons.replace_codon`; that operation replaces
exactly the three parent bases with a codon
(`backend/mutation_maker/generate_oligos.py:122-169`).

**No indel operation.** A recursive biological-term search found no source or
documentation occurrence of `indel`, `insertion`, or `deletion`. The apparent
hits were unrelated: deleting a UI table row, a test comment about matrix-row
insertion, and deleting a JavaScript viewer instance. The internal degenerate
alphabet's `_`/`GAP` entry has an empty base set
(`backend/mutation_maker/degenerate_codon.py:77,95,117`) and is not used as a
sequence-deletion operation; all output mutation paths use equal-length codons.

**Disconfirmation attempt.** I searched `README.md`, every Python file under
`backend/mutation_maker/`, all frontend `.js/.jsx` files, tests, and the full
paper text for the indel terms and for `gap`, then read the SSM/QCLM/PAS input
types, mutation parser paths, oligo generator, codon replacement, and the
paper's descriptions of all three workflows. Evidence that would change the
value: any supported event type or operation that inserts/removes DNA bases;
one event class or fixed-length only would yield `partial`, both would yield
`yes`. I searched all documented workflows and implementations and found only
same-length codon/amino-acid substitutions.

### DNA Chisel — **no** (high confidence)

**Why the apparent “insert” is not an indel.** The public
`EnforcePatternOccurence` specification uses the shorthand `insert`, but its
implementation creates an `EnforceSequence` over a location of exactly the
pattern length and constrains the existing sequence at that location
(`dnachisel/builtin_specifications/EnforcePatternOccurence.py:124-160`). It
therefore substitutes bases in a fixed-length sequence; it does not insert new
bases. The docs describe it as controlling how many times a pattern appears
(`docs/genbank/genbank_api.rst:128-140`), not as an indel-variant operation.

The most direct primary-source negative evidence is the public length
specification:

> “Quite an uncommon specification as it can't really be solved or optimized.
> But practical as part of a list of constraints to verify.”

Source:
`dnachisel/builtin_specifications/SequenceLengthBounds.py:6-20`. If the mutation
machinery could insert/delete bases, length bounds could in principle be
repaired; the class instead only evaluates the current length (`:29-36`). The
mutation-space examples and implementation likewise replace indexed spans
with variants (`dnachisel/MutationSpace/MutationChoice.py:6-35`;
`MutationSpace.py:65-86,124-130`) and expose no indel event operation.

**Disconfirmation attempt.** I recursively searched the whole repository—
README, `docs/`, `dnachisel/`, examples, and tests, excluding binary assets—for
`indel(s)`, `insert(ion)(s)`, `delete/deletion(s)`, and
`SequenceLengthBounds`. I read the complete `EnforcePatternOccurence` and
`SequenceLengthBounds` classes, the GenBank `@insert` docs, mutation-space
classes, built-in-specification index, and public `__init__` exports. Dataset
annotations mentioning biological insertion sequences/deletions were also
inspected and are merely example GenBank content. Evidence that would change
the value: a supported mutation choice or specification that lengthens or
shortens the working sequence (one event class/fixed length would yield
`partial`, both variable would yield `yes`). No such operation, parameter, or
mode exists in the checked API/docs/source.

### tangermeme — **yes** (high confidence)

**Documented insertion.** The README calls these “atomic sequence operations”
and demonstrates `insert(seq, "GCGC")`, turning `AAAAAA` into
`AAAGCGCAAA` (`README.md:84-100`). The function docstring is unambiguous:

> “the returned sequence will be longer than the original sequence.”

Source: `tangermeme/ersatz.py:20-38`. It accepts a motif of arbitrary
`motif_length`, and returns shape `length + motif_length`
(`ersatz.py:48-76`); implementation is the direct concatenation
`[left, motif, right]` (`:79-94`).

**Documented deletion.** `delete(X, start, end)` says:

> “The sequence returned from this function will be shorter than the original
> sequence.”

Source: `tangermeme/ersatz.py:300-329`. It validates caller-selected `start`
and `end` and returns the concatenated left and right flanks (`:331-340`), so
deletion length is not fixed.

There is also a higher-level, documented variant-effect interface. The README
says tangermeme “can also handle deletions and insertions” and shows
`deletion_effect` and `insertion_effect` (`README.md:258-265`). Their docstrings
say multiple single-character events can encode longer deletions/insertions
(`tangermeme/variant_effect.py:124-149,278-332`). This is independent
corroboration that the atomic operations are supported indel variant types,
not incidental tensor tricks.

**Disconfirmation attempt.** I recursively searched the entire repository—
README, source, API docs, tutorials/skills, changelog, examples, and tests—for
the indel/insertion/deletion terms; then read `ersatz.insert`, `ersatz.delete`,
both variant-effect functions, API export docs, and tests. I checked whether
“insert” actually meant same-length substitution, whether deletion merely
masked bases, and whether lengths were fixed. The README explicitly warns that
substitution is different (`README.md:100`), while the cited implementations
change tensor length and accept arbitrary motif/range lengths. Evidence that
would change the value: only one event type, a fixed event width, or only a
caller-supplied complete mutant sequence. All were searched and contradicted.

## Search log

The following searches were performed. They are recorded here so every
positive and every `no` can be reproduced. `rg` searches were case-insensitive
unless noted; large binary/image/data files were excluded from biological-term
passes and then relevant apparent hits were opened manually.

1. Read `ROW_DEFINITIONS.md` lines 1–240, then line-numbered lines 1–35 and
   125–142.
2. Used the eight `tool_survey/final/<slug>.md` records only to locate canonical
   repos/files by searching for `insert|delet|indel|github|documentation|source`;
   no statement or score from those records was used as evidence.
3. Located local primary-source checkouts and recorded each `git remote`, HEAD,
   and relevant branch.
4. PoolParty: recursive `rg` over `src`, `docs`, and `README*` for the five
   biological terms; opened both scan docs and implementations; searched for
   `generate_library`/`print_library`; ran the three behavioral cases quoted
   above with the mandated venv and bytecode disabled.
5. VaLiAnT: recursive `rg` over the full checkout for the biological terms;
   opened the README mutation/custom-variant sections and the source files
   enumerated in its disconfirmation paragraph; separately searched for an
   insertion mutator (`class .*Insertion.*Mutator`, `MutatorType.INS`,
   `insertion.*mutator`, and reversals), and searched application paths for
   `apply.*variant`, `mutat.*variant`, `variant.*seq`, `alt_ref_delta`, and
   `replace.*variant`.
6. MPRAnator: recursive repo search for the biological terms; searches for
   `SNP|VCF|variant|mutation|indel|insert|delet` in README, app, forms, views,
   templates, and downloadable script; extracted the supplied PDF with
   PyMuPDF and searched for indel and SNP terms; extracted official DOCX
   `word/document.xml` with Python `zipfile`, stripped XML tags, and searched
   for the same terms and `10 nt` variants.
7. MPRA Design Tools: recursive searches of both repositories for the
   biological terms; searched all generated `.Rd` pages for `processVCF`,
   `design`, `variant`, and `allele`; searched `processVCFfast.R` for
   `generateInsConstruct|generateDelConstruct|isDEL|isINS`; extracted the
   supplied paper PDF and searched for `indel|insertion|deletion|VCF|variant`.
8. Oligopool Calculator: recursive biological-term search; then recursive
   `variant|mutat|parent|allele|vcf|reference.*alternate|insert.*delete|delete.*insert`
   search across README, all docs, and all Python; opened the complete public
   API list and apparent positives listed in its disconfirmation paragraph.
9. Mutation Maker: recursive biological-term search excluding binary assets;
   enumerated all backend source files; recursively searched README/backend for
   `site-directed|saturation|multi-site|mutagenesis|mutation|amino|codon|oligo`;
   extracted the paper PDF and searched for all indel and mutagenesis terms;
   separately searched source/frontend for `GAP`, underscore codons, and every
   indel spelling.
10. DNA Chisel: recursive search of the whole checkout for all biological terms
    and `SequenceLengthBounds`; opened both apparent-positive classes, GenBank
    docs, mutation-space source, public exports, examples, and tests.
11. tangermeme: recursive full-repo search for all biological terms; opened the
    README, `ersatz.py`, `variant_effect.py`, API export docs, tutorial/skill
    references, changelog, and corresponding tests.

One PoolParty test attempt initially constructed deletion lengths 1 and 2 in
the same `Party` and raised the documented region-length consistency error;
the reported behavioral check therefore used a fresh `pp.init()` per length.
No repository or environment was modified by these checks.

## Row-level finding

The row **does discriminate** on one consistent scale: PoolParty, VaLiAnT,
MPRA Design Tools, and tangermeme provide both event classes without a fixed or
short-event-only restriction (`yes`); MPRAnator provides both but explicitly
caps them at short indels (`partial`); Oligopool Calculator, Mutation Maker,
and DNA Chisel provide no indel-variant operation (`no`). The crucial
cross-tool boundary is not “does the word insertion appear?” but whether a
tool operation constructs a new sequence by applying an insertion/deletion
event to a parent. This excludes Oligopool's architecture columns and DNA
Chisel's fixed-length pattern “insertion” while crediting event-level VCF
construct generation in VaLiAnT, MPRAnator, and MPRA Design Tools.
