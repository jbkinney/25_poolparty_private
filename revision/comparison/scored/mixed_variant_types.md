# Audit: mixed variant types in one library

> ## CORRECTION applied after scoring — one cell demoted
>
> `ROW_DEFINITIONS.md` gained a **global rule**: a capability counts only where the
> tool provides an operation for it, never where the user reconstructs it. The
> `partial` branch this row previously offered — *"achievable only by separate runs
> the user merges"* — was struck as a consequence, and replaced by a `partial` that
> requires the **tool** to emit one pooled output.
>
> | Tool | Was | Now | The rationale it rested on |
> |---|---|---|---|
> | tangermeme | partial | **no** | *"A user can run each operation and manually normalize/concatenate the tensors, which is precisely the separate-run, user-merge `partial` case."* |
>
> That cell cites the struck branch by name, so the demotion follows from its own
> stated reasoning; no rescoring was performed.
>
> **Deliberately NOT demoted — these need targeted rescoring:**
> - **Mutation Maker** — cites *two* grounds, one surviving: *"both independent routes to `partial` in the definition apply: same-kind components and separate runs for different workflows."* Same-kind components remains a valid `partial`, so this cannot be resolved by inspection.
> - **Oligopool Calculator**, **DNA Chisel** — `partial` on grounds that do not cite the struck branch.
> - **MPRAnator**, **MPRA Design Tools** — scored `yes`; both need checking that the pooling is tool-side rather than user-side.
>
> The per-tool sections below retain their original headings and reasoning as
> written at scoring time; the Scores table has been updated.

Audit date: 2026-08-15

## Measurement instrument and operational test

I used row 6 of `revision/tables/ROW_DEFINITIONS.md` exactly as written:

> "Two or more structurally different component types declared in **one**
> specification and emitted as **one** pooled output. `partial` = achievable only
> by separate runs the user merges, or limited to two component types of the same
> kind."

Before inspecting any tool, I fixed the following one-sentence operational test:

> **For each tool, check whether one documented design specification can declare at
> least two structurally different library component types and one documented
> execution emits them together as a single pooled output; score `partial` only
> when this requires separate user-merged runs or supports merely two components
> of the same kind, and score `no` only after searching the relevant API,
> examples, and documentation for such a declaration-and-pooled-emission path.**

I treated a VCF as a declaration rather than a pre-made sequence collection when
the tool reads REF/ALT and itself constructs the corresponding reference and
alternate oligos. Conversely, a CSV of already-complete sequences is not a
declaration of structural variant types. This distinction is applied identically
to MPRAnator, MPRA Design Tools, and Oligopool Calculator.

## Scores

| tool | value | confidence |
|---|---|---|
| poolparty | **yes** | high |
| valiant | **yes** | high |
| mpranator | **yes** | high |
| mpradesign | **yes** | high |
| oligopoolcalc | **partial** | high |
| mutationmaker | **partial** | high |
| dnachisel | **partial** | medium |
| tangermeme | **no** [corrected] | high |

## Per-tool audit

### poolparty — yes (high confidence)

**Why.** The documented design makes substitution mutants and deletion variants
from the same template, passes both branch objects in one `stack([...])`
specification, and emits the result through one final `Pool`.

**Quoted primary evidence.** The package's landing-page example is explicitly
introduced as:

> "Stack different variant types into a single barcoded library:"

and the code is:

```python
mutations = template.mutagenize(region="cre", num_mutations=1)
deletions = template.deletion_scan(region="cre", deletion_length=5)
combined = pp.stack([mutations, deletions])
library = combined.insert_kmers(region="bc", length=10).named("library")
library.print_library(num_seqs=6, seed=0)
```

Source: [`poolparty/docs/index.rst:61-77`](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/docs/index.rst#L61-L77).
The rendered example output immediately below contains both a highlighted base
substitution and a five-base deletion in that one library
([lines 81-88](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/docs/index.rst#L81-L88)).

The implementation defines `stack` as a disjoint-union constructor:

> "Sequence of Pool objects to stack into a single Pool."

> "A Pool whose states are the disjoint union of all input pools' states."

Source: [`poolparty/src/poolparty/state_ops/stack.py:19-48`](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/src/poolparty/state_ops/stack.py#L19-L48).

**Behavioral check.** Per the task instruction, I ran the source with the provided
read-only Python 3.12 environment and `PYTHONDONTWRITEBYTECODE=1`. From a four-base
tagged template, `mutagenize(..., mode="sequential")` produced 9 members,
`deletion_scan(..., deletion_length=2, mode="sequential")` produced 2, and
`pp.stack([m, d]).generate_library()` returned 11 rows. The rows included nine
same-length substitutions plus `A<tag>--T</tag>` and
`A<tag>C--</tag>` in the same DataFrame.

Reproducible command (run from `/tmp`):

```bash
PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/src \
/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/.venv/bin/python -u -c \
'import poolparty as pp; pp.init(); t=pp.from_seq("A<tag>CGT</tag>"); m=t.mutagenize(region="tag",num_mutations=1,mode="sequential"); d=t.deletion_scan(region="tag",deletion_length=2,mode="sequential"); c=pp.stack([m,d]); print(m.num_states,d.num_states,c.num_states); print(c.generate_library().to_string(index=False))'
```

**Disconfirmation attempt.** I searched all public docs and package source under
`poolparty/docs` and `poolparty/src/poolparty` for
`stack|mix|pool|wild.?type|control|variant|mutation|shuffle|scan|specification`,
then inspected `docs/index.rst`, `docs/quickstart.rst`, and
`state_ops/stack.py`. I specifically looked for a restriction that `stack`
requires homogeneous operation ancestry or equal sequence length; none exists.
The source instead assigns `seq_length=None` when parent lengths differ
(`stack.py:66-71`).

**What would change the value.** `partial` would require evidence that the two
branches must be generated in separate executions and merged outside the API, or
that `stack` accepts only components of one structural kind; `no` would require
the displayed branches not to be emitted together. I searched for those
restrictions and ran the mixed branch case; the evidence falsifies both.

### VaLiAnT — yes (high confidence)

**Why.** One targeton row has one `action_vector` that declares amino-acid codon
substitutions (`aa`) and in-frame deletions (`inframe`); the corresponding one
targeton output CSV contains both member types.

**Quoted primary evidence.** The input format documentation states:

> "Multiple types of mutations can be applied to each target region."

and defines `action_vector` as mutation-type labels grouped by target region
([`README.md:630-650`](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/README.md#L630-L650)).
The documented example itself combines deletion and SNV labels in one action
vector:

```text
(1del), (1del, snv), (1del)
```

Source: [`README.md:652-657`](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/README.md#L652-L657).

The shipped BRCA1 peptide example is even more direct. Every input row declares:

```text
"(),(aa,inframe),()"
```

Source: [`examples/sge/brca1/parameter_input_files/brca1_pep_targeton_input.txt:1-5`](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/examples/sge/brca1/parameter_input_files/brca1_pep_targeton_input.txt#L1-L5).

In the single output file
`chr17_43106355_43106599_minus_sgRNA_ex3_meta.csv`, line 2 is an `..._aa_rc`
member whose `mutator` field is `aa`, while line 321 is an
`..._inframe_rc` member whose `mutator` field is `inframe`. Thus the two
structurally different types declared on input row 3 are pooled in one targeton
metadata/output table. Source:
[`...ex3_meta.csv:2`](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/examples/sge/brca1/brca1_pep_output_exp/chr17_43106355_43106599_minus_sgRNA_ex3_meta.csv#L2)
and
[`...ex3_meta.csv:321`](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/examples/sge/brca1/brca1_pep_output_exp/chr17_43106355_43106599_minus_sgRNA_ex3_meta.csv#L321).

**Disconfirmation attempt.** I searched `README.md`, `examples/`, and `src/` for
`mut_type|mutation type|variant type|wild.?type|control|deletion|insertion|custom|output|input`,
inspected the mutation-type section, both targeton formats, the shipped targeton
inputs, and expected output files, and searched source/output for `action`,
`mutator`, `1del`, `inframe`, and `custom`. I looked specifically for a rule that
an action vector may contain only one label or that each label is written to a
separate output; the quoted docs and shipped paired input/output refute both.

**What would change the value.** `partial` would require `aa` and `inframe` to be
separate invocations or separate user-merged pools; `no` would require only one
label per targeton. I searched for both restrictions and found the opposite in
the format definition and canonical example. A future finding that the expected
CSV is assembled outside VaLiAnT would also lower the score; the repository
presents it as the expected output of the one example run.

### MPRAnator — yes (high confidence)

**Why.** A single MPRA-SNP POST accepts a multi-row VCF-like `SnpS` declaration.
Each row's REF and ALT lengths can independently specify substitution, insertion,
or deletion; the view processes every row into one `finalOutput` and serializes
all members into one text/FASTA-like response.

**Quoted primary evidence.** The input parser treats `SnpS` as multiple lines:

```python
firstSplitL = self.fileS.strip().split("\n")
eachByLineL = [re.split(r"\s+", i.strip()) for i in firstSplitL]
...
self.snps = eachByLineL
```

Source: [`mycustom.py:256-280`](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/mycustom.py#L256-L280).
The mutation builder states:

> "if the SNP contains a small deletion or insertion (smaller than 10nt) we
> either remove part of the sequence or we insert adenines in one edge"

and, for every supplied record, computes `difference = initial_size-SNP_size`
before taking distinct negative, zero, and positive branches
([`parseSNPs.py:29-31,137-186`](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/parseSNPs.py#L29-L31)).

The one request path shows pooled emission:

```python
snpO = mc.SnpFile(snpS, fileName=False)
subbedSequenceHeadersL, subbedSequencesL = parseSNPs.make_sequence_copies(...)
finalOutput = []
for headerS, seqS in zip(...):
    finalOutput.append({"header": headerS, "sequence": seqS})
mpraOutput, mpraOutputHTML = part1.createMPRAResultOutput(finalOutput, ...)
...
response.write(mpraOutput)
return response
```

Source: [`iliasApp/views.py:319-381`](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/iliasApp/views.py#L319-L381).
The serializer appends each item to the same `response` string and returns it
([`part1.py:210-292`](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/part1.py#L210-L292)).

**Disconfirmation attempt.** I inventoried the complete repository and searched
all Python/templates/docs for
`SNP|indel|insert|delet|variant|wild.?type|negative control|transmutation|motif|FASTA|result`.
I inspected `SnpFile`, `make_sequence_copies`, the MPRA-SNP form, request view,
serializer, and downloadable client. I looked for a per-request single-variant
restriction and for separate response objects by type; neither exists. The only
relevant restriction found is maximum ALT length 10 in `forms.py:542-545`, which
limits size, not mixing.

**What would change the value.** `partial` would require one POST per structural
type or an external concatenation step; `no` would require the parser to reject
mixed REF/ALT length relations. I searched validation and request code for both;
the parser iterates the common list and the request writes one combined response.

### MPRA Design Tools — yes (high confidence)

**Why.** One documented `processVCF(vcf=...)` run accepts a VCF whose rows may be
SNVs, insertions, or deletions, dispatches each row by REF/ALT structure, row-binds
the resulting constructs, and optionally writes one TSV.

**Quoted primary evidence.** The README says `processVCF` designs the set of
barcoded sequences and documents input handling for indels:

> "Currently the main function of MPRA Design Tools package is to design a set of
> barcoded sequences for MPRA experiments ... This is done with the `processVCF`
> function."

> "Insertions and deletions must encode the reference and alternate alleles
> (respectively) as a dash character '-'."

Source: [`README.md:39-58`](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/README.md#L39-L58).
The documented API example is one `processVCF` call with one input VCF and one
output TSV ([`README.md:106-120`](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/README.md#L106-L120)).

The row processor explicitly distinguishes three structural types:

```r
isSNV = (snp$REF %in% c('A', 'C', 'G', 'T') && snp$ALT %in% c('A', 'C', 'G', 'T'))
isINS = snp$REF == '-'
isDEL = snp$ALT == '-'
```

and notes that the implementation contains three branches because there are
"subtle differences in each of the types of variants"
([`R/processVCFfast.R:232-242`](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L232-L242)).
`processVCF` applies `processSnp` row-wise to the whole VCF, then performs
`Reduce('rbind', .)` and `write_tsv(successes, path = outPath)`
([`R/processVCFfast.R:1221-1259`](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L1221-L1259)).

**Disconfirmation attempt.** I searched both official repositories
(`andrewGhazi/mpradesigntools` and `andrewGhazi/designMPRA`) across R, Rmd,
Markdown, and Rd files for
`variant|mutation|mutagen|control|shuffle|random|motif|SNP|indel|insert|delet|pool|library|output|construct`.
I inspected the README input contract, `processSnp`'s SNV/insertion/deletion
branches, and `processVCF` aggregation/output. I looked for homogeneous-VCF
validation or a separate output path per type; none was found.

**What would change the value.** `partial` would require users to split a mixed
VCF and concatenate outputs; `no` would require only one of SNV/insertion/deletion
to be supported. I searched the validation, dispatch, and aggregation code for
those conditions; the common row-wise dispatch and one `Reduce('rbind', .)` path
show the opposite.

### Oligopool Calculator — partial (high confidence)

**Why.** The documented API consumes already-authored complete variant sequences;
it does not declare different mutation structures. The official example-only
mutant generator offers several substitution-library strategies, but its API
selects exactly one `strategy` per call. Mixed use therefore requires separate
calls and a caller merge, and the available strategies are all substitution
regimes. That is exactly the row definition's `partial` case, with the additional
rule that example-only support cannot exceed `partial`.

**Quoted primary evidence.** The core `compress` API requires pre-existing
variants:

> "`input_data` ... DataFrame with annotated oligo pool variants."

> "All non-'ID' columns are concatenated (left-to-right) to form the sequence per
> variant."

Source: [`oligopool/compress.py:16-55`](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/compress.py#L16-L55).
The official design-parser example is similarly explicit that variants are
supplied, not structurally declared:

```python
if elements_spec[element_name]['element_type'] == 'variant':
    dataframe[element_name] = elements_spec[element_name]['sequences']
```

Source: [`examples/design-assembly-parser/design_assembly_parser.py:141-170`](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/examples/design-assembly-parser/design_assembly_parser.py#L141-L170).

The separate official example script does generate variants, but only one
substitution strategy per call:

```python
def generate_variant_dataframe(..., strategy: str = 'codon', ...):
    ...
    if strategy == 'single': ...
    elif strategy == 'codon': ...
    elif strategy == 'multi': ...
```

Source: [`examples/library-compressor/mutant_generator.py:200-233`](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/examples/library-compressor/mutant_generator.py#L200-L233).
The file describes these as single-nucleotide, codon, and multi-position variants
([lines 43-176](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/examples/library-compressor/mutant_generator.py#L43-L176)); all are substitution spaces.

**Disconfirmation attempt.** I searched README, all docs, all examples, and all
package modules for
`variant|mutation|mutagen|wild.?type|control|pool|library|input|output|element|insert|delet|mix|hetero`,
then separately for
`mutat|substitut|delet|insert|variant type|control|wild.?type`. I inspected the
YAML pipeline, `variants.csv`, design parser, compressor, generator example,
`join`, and `merge`. `join` inserts columns by matching ID and `merge` combines
columns within rows; neither stacks heterogeneous row populations. No public
function or YAML command declares mutation-event types or unions variant
populations.

**What would change the value.** `yes` would require a documented core API/YAML
field accepting two structural regimes (for example substitution plus deletion)
and producing their row union in one output. I searched every design/degenerate
command and example for such a dispatcher/union and found none. `no` would require
even separate generation of two regimes to be absent; the official example does
provide multiple same-kind strategies that can be called separately and merged,
so `no` would understate the tool under the row's explicit partial rule.

### Mutation Maker — partial (high confidence)

**Why.** Mutation Maker exposes three separate workflow endpoints, not one mixed
workflow. Within SSM/QCLM, one input can list multiple amino-acid substitutions,
but they are components of the same structural kind. Thus both independent routes
to `partial` in the definition apply: same-kind components and separate runs for
different workflows.

**Quoted primary evidence.** The API defines separate endpoints:

> "Return task if type Site Saturation Mutagenesis (SSM)"

> "Return task if type QuikChange Lightning Multi Site-Directed Mutagenesis
> (QCLM)"

> "Return task if type PCR-based Accurate Synthesis task (PAS)"

Source: [`api/api.py:40-61`](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/api/api.py#L40-L61).
The worker mirrors them as distinct task/input types
([`backend/tasks.py:42-73`](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/tasks.py#L42-L73)).

Within one SSM input, `mutations` is a list, but every entry is parsed by the same
codon-mutation parser:

```python
mutations = ListProperty(str, required=True)
...
return [parse_codon_mutation(mutation, goi_offset) for mutation in self.mutations]
```

Source: [`backend/mutation_maker/ssm_types.py:211-218`](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/ssm_types.py#L211-L218).
QCLM likewise parses strings such as `E32W E32L E49K` and combines amino-acid
substitution choices by site
([`qclm_types.py:101-107`](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/qclm_types.py#L101-L107)).

**Disconfirmation attempt.** I searched README, docs, API, backend, and tests for
`SSM|MSDM|PAS|mutation type|variant|insert|delet|substitut|wild.?type|control|pool|library|output|workflow|method`,
then inspected all API routes, Celery tasks, and SSM/QCLM/PAS input/output types.
I looked for a union request that selects multiple workflows and for
nucleotide-length-changing indel mutation declarations. The routes accept one
workflow-specific input class each; no mixed-workflow request or indel design
type was found. Incidental `insert`/`delete` hits were data-structure or UI terms,
not mutation-event support.

**What would change the value.** `yes` would require one request object and one
output spanning structurally different mutation events or workflow types; I
searched API routes, task dispatch, and all three input schemas for that and found
only separate paths. `no` would require the tool not to emit multiple declared
components even within one kind; the SSM mutation list and QCLM multi-amino input
show that it does, so the row's same-kind `partial` clause applies.

### DNA Chisel — partial (medium confidence)

**Why.** A single `MutationSpace` can enumerate a pooled stream formed from the
Cartesian product of several segment choices, but those choices are replacement
subsequences—the same structural event kind. The ordinary documented optimizer
collapses a composed specification to one final sequence, and multiple optimized
design regimes would require caller-run problems and caller aggregation. This is
`partial`, not `yes`.

**Quoted primary evidence.** DNA Chisel describes itself as optimizing one DNA
sequence against a composed list of specifications:

> "DNA Chisel ... is a Python library for optimizing DNA sequences with respect
> to a set of constraints and optimization objectives."

> "The library comes with over 15 classes of sequence specifications which can be
> composed ... or all of this at once!"

Source: [`README.rst:18-31`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/README.rst#L18-L31).
`DnaOptimizationProblem` takes one `sequence`, a constraint list, and an objective
list ([`DnaOptimizationProblem.py:47-71`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py#L47-L71)).

The restricted pooled capability is real:

```python
def all_variants(self, sequence):
    """Iterate through all sequence variants in this mutation space."""
    ...
    for variants in itertools.product(*variants_slots):
        ...
        new_sequence[start:end] = variant
        yield new_sequence.decode()
```

Source: [`MutationSpace/MutationSpace.py:132-164`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/MutationSpace/MutationSpace.py#L132-L164).
The mutation model is a segment plus replacement variants
([`MutationChoice.py:6-35`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/MutationSpace/MutationChoice.py#L6-L35)); the optimization-problem constructor's own example is nucleotide
replacement choices and codon replacement choices
([`DnaOptimizationProblem.py:100-106`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py#L100-L106)).

**Disconfirmation attempt.** I searched README, complete docs, all examples, and
the package tree for
`variant|mutation|insert|delet|indel|substitut|library|pool|batch|collection|all_variants|random_compatible|optimize`,
then made a narrower pass for
`deletion|deletions|insertion|insertions|indel|library|pool|variant type|variants`.
I inspected `DnaOptimizationProblem`, `MutationChoice`, `MutationSpace`, the CLI,
built-in specification reference, and every documented example category. The
`insert` shorthand is `EnforcePatternOccurence`, which replaces bases while
preserving the one problem sequence; no deletion event type, heterogeneous
library union, pool writer, or type-labelled member output was found. The CLI is
explicitly one `<source>` to one `<target>` (`dnachisel/cli.py:3-14`).

**What would change the value.** `yes` would require a documented specification
that labels at least two structurally different member types and emits them as
one collection; I searched all specs, docs, examples, output/report APIs, and the
CLI for that and found none. `no` would require even a same-kind multi-component
collection to be absent; the public `MutationSpace.all_variants()` Cartesian
iterator supplies that restricted case. Confidence is medium rather than high
because `all_variants` is principally a solver primitive rather than a named
library-design workflow, but it is public and directly documented in source, so
ignoring it would understate the tool under the row's explicit same-kind clause.

### tangermeme — partial (high confidence)

**Why.** Tangermeme has documented substitution, insertion, and deletion
operations, but they are separate calls with incompatible output lengths and no
documented heterogeneous library-union API. A user can run each operation and
manually normalize/concatenate the tensors, which is precisely the separate-run,
user-merge `partial` case.

**Quoted primary evidence.** The README presents the variant-effect operations as
separate calls:

```python
y, y_var = substitution_effect(model, X, substitutions)
...
y, y_var = deletion_effect(model, X, deletions)
y, y_var = insertion_effect(model, X, insertions)
```

and states that deletion/insertion require trimming because before/after
sequences differ in length. Source: [`README.md:238-265`](https://github.com/jmschrei/tangermeme/blob/8d732b8c08764057b7ae5faa644d48664f36b44b/README.md#L238-L265).

The sequence-returning atomic APIs are likewise separate. `insert` returns
`length + motif_length`
([`ersatz.py:20-94`](https://github.com/jmschrei/tangermeme/blob/8d732b8c08764057b7ae5faa644d48664f36b44b/tangermeme/ersatz.py#L20-L94)), `substitute` returns the original length
([`ersatz.py:97-166`](https://github.com/jmschrei/tangermeme/blob/8d732b8c08764057b7ae5faa644d48664f36b44b/tangermeme/ersatz.py#L97-L166)), and `delete` returns
`length-(end-start)`
([`ersatz.py:300-340`](https://github.com/jmschrei/tangermeme/blob/8d732b8c08764057b7ae5faa644d48664f36b44b/tangermeme/ersatz.py#L300-L340)). There is no call in these APIs that takes a tagged union of operation kinds.

**Disconfirmation attempt.** I searched README, all RST docs/tutorial indices,
all package modules, and bundled skill docs for
`substitute|insert|delet|variant|mutation|mutagen|product|design|stack|concat|cat\(|wild.?type|control|mixed|hetero`,
then specifically for
`mixed|heterogeneous|combine|combined|concatenate|concat|stack|pool|library`.
I inspected `ersatz`, all three `variant_effect` functions, `product`, `space`,
`screen`, and design modules. `product` combines sequence/model-input axes, not
variant-operation types; `torch.cat` hits are internal batch assembly or require
the caller to supply already-compatible tensors. No heterogeneous
variant-specification dispatcher or pooled sequence output was found.

**What would change the value.** `yes` would require a documented call accepting,
for example, substitution and deletion declarations together and returning one
type-labelled pooled sequence tensor/table; I searched the variant-effect,
ersatz, product, and design APIs for it and found none. `no` would require the
structurally different outputs not to be achievable at all; the separate
documented APIs produce them, and caller-side trimming/concatenation can merge
them, so the definition mandates `partial`.

## Row-level finding

The row **does discriminate** on the supplied tool set: four tools (`poolparty`,
`valiant`, `mpranator`, `mpradesign`) natively declare and emit heterogeneous
member types in one execution, while four receive only restricted/manual credit.
The distinction is reproducible at the API boundary: native union/aggregation is
`yes`; one-kind enumeration or separate API calls followed by caller aggregation
is `partial`.

There are no `no` or `unknown` cells. This is not a failure of the row: the
definition deliberately makes separate-run user merging a `partial`, and all four
restricted tools provide enough variant-generation or variant-transformation
machinery to meet that clause. The row's main limitation is that its `partial`
bucket spans a wide range—from a public same-kind Cartesian iterator (DNA Chisel)
to multiple genuinely structural operations with only caller-side union
(tangermeme).

## Search and provenance log

No conclusion or rating from any prior assessment was used. The files under
`revision/tool_survey/final/` were grepped only once, before tool inspection, for
repository/documentation URLs as permitted leads; they are not evidence anywhere
above. I did not inspect `MATRIX_verified.md`, `ROWS_v2.md`, `ROWS_v3.md`, or any
previous heterogeneous-row audit.

### Instrument and source discovery

1. Read `revision/tables/ROW_DEFINITIONS.md` with `sed -n '1,240p'`; fixed the
   operational test before inspecting any tool.
2. Inventoried `revision/tool_survey` with `find ... -maxdepth 3 -type f` and
   checked top-level directories with `find . -maxdepth 2 -type d`.
3. Lead-only URL search: for the eight `final/<slug>.md` files, searched
   `https?://|github|docs|readme|source|paper|repository`.
4. Located already-cloned primary repositories under `/tmp` with
   `find /tmp -maxdepth 2/3 -type d | rg -i
   '(valiant|mpranator|mpradesign|designmpra|oligopool|mutation.?maker|dnachisel|tangermeme)'`.
5. For every repository, recorded `remote get-url origin`, `rev-parse HEAD`, and
   branch. Audited commits: PoolParty `1bb0179`; VaLiAnT `8796cc1` (`develop`);
   MPRAnator `9969790`; mpradesigntools `afd386e`; designMPRA `0cf56ee`;
   Oligopool `b88fa39`; Mutation Maker `396c1c0`; DNA Chisel `68c0930`;
   tangermeme `8d732b8` (`main`).

### Tool searches (exact search expressions and scopes)

- **PoolParty:** `rg -n -i
  'stack|mix|pool|wild.?type|control|variant|mutation|shuffle|scan|specification'
  poolparty/docs poolparty/src/poolparty`; inspected
  `docs/index.rst:55-95`, `docs/quickstart.rst:256-285`,
  `docs/quickstart.rst:335-442`, `docs/pool.rst`, and
  `src/poolparty/state_ops/stack.py:1-125`. Ran the behavioral command described
  in the PoolParty section from `/tmp`, with the repository read-only.
- **VaLiAnT:** inventoried files with `find /tmp/VaLiAnT -maxdepth 3 -type f`;
  searched README/examples/src for
  `mut_type|mutation type|variant type|wild.?type|control|deletion|insertion|custom|output|input`;
  searched non-VCF examples for `action|mutator|1del|inframe|custom`; inspected
  `README.md:245-340,455-490,626-675,688-735`, targeton inputs, and expected
  `_meta.csv`/`_unique.csv` outputs; located representative output rows with
  `rg -n ',aa,'` and `rg -n ',inframe,'`.
- **MPRAnator:** inventoried files with `find /tmp/MPRAnator -maxdepth 3 -type f`;
  searched Python/templates/docs for
  `SNP|indel|insert|delet|variant|wild.?type|negative control|transmutation|motif|FASTA|result`;
  inspected `mycustom.py:241-305`, `parseSNPs.py:1-205,395-456`,
  `forms.py:501-545`, `views.py:225-390`, `part1.py:210-330`, the SNP form,
  results template, and downloadable SNP client. Also searched `views.py` for
  `def.*Snp|ViewMpraSnp|forDownload|sequenceHTML` and `parseSNPs.py` for
  `difference|deletion|insertion|small`.
- **MPRA Design Tools:** inventoried both official repositories; searched all
  R/Rmd/Markdown/Rd sources for
  `variant|mutation|mutagen|control|shuffle|random|motif|SNP|indel|insert|delet|pool|library|output|construct`;
  inspected README, package manual, `processVCF`, all three `processSnp` branches,
  output aggregation, and the Shiny server.
- **Oligopool Calculator:** inventoried docs/examples/package modules; searched
  README/docs/examples/package for
  `variant|mutation|mutagen|wild.?type|control|pool|library|input|output|element|insert|delet|mix|hetero`
  and then
  `mutat|substitut|delet|insert|variant type|control|wild.?type`; inspected
  `compress.py`, `expand.py`, `join.py`, `merge.py`, the design parser, library
  compressor example, serial YAML pipeline, and `variants.csv`.
- **Mutation Maker:** inventoried docs/backend/API; searched them for
  `SSM|MSDM|PAS|mutation type|variant|insert|delet|substitut|wild.?type|control|pool|library|output|workflow|method`;
  inspected README workflow statements, `api.py`, `server_fastapi.py`,
  `tasks.py`, and `ssm_types.py`, `qclm_types.py`, `pas_types.py`, plus tests for
  the corresponding input/output structures.
- **DNA Chisel:** searched README/docs/examples/package for
  `variant|mutation|insert|delet|indel|substitut|library|pool|batch|collection|all_variants|random_compatible|optimize`
  and separately
  `deletion|deletions|insertion|insertions|indel|library|pool|variant type|variants`;
  inspected README, `DnaOptimizationProblem`, `MutationChoice`, `MutationSpace`,
  `EnforcePatternOccurence`, CLI, built-in-specification docs, and all categories
  in `docs/examples.rst`.
- **tangermeme:** inventoried package/docs; searched README/docs/package for
  `substitute|insert|delet|variant|mutation|mutagen|product|design|stack|concat|cat\(|wild.?type|control|mixed|hetero`
  and then
  `mixed|heterogeneous|combine|combined|concatenate|concat|stack|pool|library`;
  inspected README variant-effect section, `ersatz.py`, all three functions in
  `variant_effect.py`, `product.py`, `space.py`, and design modules.

All quoted evidence is from these primary repositories. Local PDFs were not
needed to settle this row, so no PDF quotation is used.
