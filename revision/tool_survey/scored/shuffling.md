# Row 18 audit — Shuffling

Audit date: 2026-08-17.

## Binding operational test

This audit applies row 18 and the global rule in `ROW_DEFINITIONS.md` literally:

- `yes`: the tool itself exposes a dinucleotide-preserving shuffle of an input sequence;
- `partial`: the tool itself exposes only a mononucleotide shuffle or scramble;
- `no`: no tool-provided shuffling operation exists; and
- `unknown`: the primary sources genuinely do not permit a determination.

The operation must permute residues of an input sequence. Reverse complementing or
orientation flipping does not count. Nor do random sequence generation, random search
inside an optimizer, random assignment of barcodes, or user-written composition of
unrelated primitives. A source-level helper counts only when it backs a tool-provided
operation, parameter, or mode.

## Results

| Tool | Value | Confidence | Operational finding |
|---|---|---|---|
| PoolParty | `yes` | high | Public `shuffle_seq(..., shuffle_type="dinuc")` performs an Euler-path, dinucleotide-preserving shuffle. |
| VaLiAnT | `no` | high | The complete mutator enum/dispatcher has deletion, SNV, and codon mutators, but no sequence shuffle. |
| MPRAnator | `partial` | high | Its Transmutation “Scramble” option applies Python's in-place character shuffle, preserving mononucleotide counts only. |
| MPRA Design Tools | `no` | high | Its public design API has no input-sequence shuffle; an internal shuffle only randomizes barcode assignment. |
| Oligopool Calculator | `no` | high | The complete Python/CLI operation registry has no shuffle; internal shuffles concern search/counting state, not input-sequence variants. |
| Mutation Maker | `no` | medium-high | Its complete web API exposes only SSM, QCLM/MSDM, and PAS design workflows; the surviving source contains no shuffle/scramble operation. |
| DNA Chisel | `no` | high | Neither the complete public API nor built-in specifications expose shuffling; the sole source hit shuffles a list of optimization locations. |
| tangermeme | `yes` | high | Public `dinucleotide_shuffle` returns shuffled input sequences and its tests verify exact dinucleotide-count preservation. |

## Per-tool evidence and disconfirmation

### PoolParty — `yes` (high confidence)

**Quoted primary-source evidence.** The public function signature contains
`shuffle_type: Literal["mono", "dinuc"]`, and its API documentation defines `"dinuc"`
as an “Euler-path shuffle preserving dinucleotide frequencies.” In the implementation,
the `dinuc` branch calls `dinucleotide_shuffle` on the molecular characters and writes
the returned permutation back into the sequence. The dedicated test repeatedly asserts
`dinuc_counts(result) == dinuc_counts(seq)`.

**Location.** PoolParty v0.1.1, commit
[`1bb0179e1c3720b1fffd471802b3040f9336de28`](https://github.com/jbkinney/poolparty-statetracker/commit/1bb0179e1c3720b1fffd471802b3040f9336de28):
[`poolparty/src/poolparty/base_ops/shuffle_seq.py`, lines 14–60 and 177–207](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/src/poolparty/base_ops/shuffle_seq.py#L14-L207),
[`poolparty/src/poolparty/utils/shuffle_utils.py`, lines 8–28](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/src/poolparty/utils/shuffle_utils.py#L8-L28), and
[`poolparty/tests/test_dinuc_shuffle.py`, lines 37–58](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/tests/test_dinuc_shuffle.py#L37-L58).
Repository-local equivalents were inspected at
`poolparty-statecounter/poolparty/src/poolparty/base_ops/shuffle_seq.py`,
`poolparty-statecounter/poolparty/src/poolparty/utils/shuffle_utils.py`, and
`poolparty-statecounter/poolparty/tests/test_dinuc_shuffle.py`.

**Reasoning.** This is a public tool operation that accepts an input sequence/pool,
performs the permutation itself, and explicitly preserves dinucleotide frequencies.
It meets the `yes` threshold directly.

**Disconfirmation attempt.** I searched the package source, public operation docs, and
tests for `shuffle`, `scramble`, `dinucleotide`, and related permutation terms, then
checked whether `dinuc` was merely an internal utility, example-only behavior, or a
mislabelled mononucleotide shuffle. It is wired into the public `shuffle_seq` operation
and backed by an exact composition test. Evidence that `shuffle_type="dinuc"` was not
public, did not act on the input, or failed to preserve dinucleotide counts would have
changed this to `partial` or `no`.

### VaLiAnT — `no` (high confidence)

**Quoted primary-source evidence.** The documented mutation types are “parametric
deletion,” “single-nucleotide variant,” “in-frame deletion,” and codon substitutions.
The source's exhaustive `MutatorType` enum contains only `del`, `snv`, `snvre`,
`inframe`, `ala`, `stop`, and `aa`; `MutatorBuilder.MUTATOR_CLASSES` dispatches exactly
those seven types and otherwise raises `NotImplementedError`.

**Location.** VaLiAnT v4.0.0, commit
[`8796cc112dafd4919fec59913f58cd2be87c45eb`](https://github.com/cancerit/VaLiAnT/commit/8796cc112dafd4919fec59913f58cd2be87c45eb):
[`README.md`, lines 255–270](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/README.md#L255-L270),
[`src/valiant/mutator_type.py`, lines 22–32](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/src/valiant/mutator_type.py#L22-L32), and
[`src/valiant/mutator.py`, lines 40–73](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/src/valiant/mutator.py#L40-L73).

**Reasoning.** VaLiAnT has no operation or mode that permutes the residues of an input
sequence. Its reverse-complement support is excluded by the row boundary and cannot be
reconstructed into shuffling under the global rule.

**Disconfirmation attempt.** I searched the complete repository (README, CLI/config
code, mutator loaders/enum/dispatcher, all `src/valiant`, examples, and tests) for
`shuffle`, `shuffling`, `scramble`, `permutation`, and `dinucleotide`. The only
dinucleotide hits are descriptive INFO fields in a supplied example VCF; no shuffle or
scramble implementation was found. A registered mutator or documented mode that
permuted an input targeton would have changed the value (to `yes` if dinucleotide
preserving, otherwise `partial`).

### MPRAnator — `partial` (high confidence)

**Quoted primary-source evidence.** The Transmutation documentation lists
“Scrambling sequences.” Its form exposes `scramble = forms.BooleanField(...)`. When
selected, the view calls `part3.scramble_motifs(seq)` for every input sequence. That
function converts the sequence to a character list, calls `rd.shuffle(motif)`, and
joins the characters back into a string.

**Location.** MPRAnator commit
[`9969790d62410138d4281b7955da6d085f07b1bc`](https://github.com/hemberg-lab/MPRAnator/commit/9969790d62410138d4281b7955da6d085f07b1bc):
[`iliasApp/templates/iliasApp/docs.html`, lines 138–159](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/iliasApp/templates/iliasApp/docs.html#L138-L159),
[`iliasApp/forms.py`, lines 340–349](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/iliasApp/forms.py#L340-L349),
[`iliasApp/views.py`, lines 161–194](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/iliasApp/views.py#L161-L194), and
[`part3.py`, lines 1–9](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/part3.py#L1-L9).

**Reasoning.** Python's `random.shuffle` permutes the character list, so the tool itself
produces a mononucleotide-composition-preserving scramble of each supplied sequence.
There is no dinucleotide-preservation logic. This is exactly `partial`.

**Disconfirmation attempt.** I searched the whole repository and the first-party tool
documentation for `shuffle`, `scramble`, `dinucleotide`, and permutation algorithms,
then traced the form control through the view to the implementation. The unrelated
`itertools.permutations` code enumerates orders of motifs for the motif-placement module
and is excluded by the row-9 boundary. A documented or implemented dinucleotide mode,
transition-count preservation, Euler-path algorithm, or exact dinucleotide-count test
would have changed the value to `yes`; evidence that the form did not reach the helper
would have changed it to `no`.

### MPRA Design Tools — `no` (high confidence)

**Quoted primary-source evidence.** The first-party README states that the package's
“main function” is `processVCF`. The package namespace exports only `processVCF` and
`spread_and_fix_indels`. The apparent shuffle hit is explicitly under “Create a pool of
barcodes for each snp” and permutes `mers` before splitting those barcodes among VCF
rows.

**Location.** `mpradesigntools` v0.2.0, commit
[`afd386ef12051bb0a37ad63a6f92acd555246584`](https://github.com/andrewGhazi/mpradesigntools/commit/afd386ef12051bb0a37ad63a6f92acd555246584):
[`README.md`, lines 39–59](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/README.md#L39-L59),
[`NAMESPACE`, lines 1–14](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/NAMESPACE#L1-L14), and
[`R/processVCFfast.R`, lines 1211–1219](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L1211-L1219).
The companion Shiny source was also checked at commit
[`0cf56ee602fc86dde705906d071a39cbdf6e99a8`](https://github.com/andrewGhazi/designMPRA/commit/0cf56ee602fc86dde705906d071a39cbdf6e99a8).

**Reasoning.** Randomizing which pre-existing barcode is assigned to a variant does not
generate a shuffled variant of an input sequence. Neither exported operation performs a
sequence shuffle, and the user cannot receive credit by composing Biostrings utilities.

**Disconfirmation attempt.** I searched all R/package documentation, exported
functions, the full package source, and the complete companion Shiny repository for
`shuffle`, `scramble`, `dinucleotide`, `permutation`, and sequence-randomization terms.
I inspected the sole relevant hit rather than treating its name as capability evidence;
it shuffles the barcode vector only. A package/app control that permuted the REF/ALT
sequence context itself would have changed the result.

### Oligopool Calculator — `no` (high confidence)

**Quoted primary-source evidence.** The public registry lists the complete commands:
`background`, `barcode`, `primer`, `motif`, `spacer`, `split`, `pad`, `merge`,
`revcomp`, `join`, QC/degenerate operations, and read-analysis operations. The README
describes Design Mode as design of “barcodes, primers, motifs/anchors, and spacers.” No
shuffle operation appears in either inventory.

**Location.** Oligopool Calculator v2026.02.22.1, commit
[`b88fa394ca67ed4c48ec41127b5707694ee7cf0a`](https://github.com/ayaanhossain/oligopool/commit/b88fa394ca67ed4c48ec41127b5707694ee7cf0a):
[`oligopool/cli.py`, lines 54–87](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/cli.py#L54-L87),
[`oligopool/__init__.py`, lines 40–70](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/__init__.py#L40-L70), and
[`README.md`, lines 34–48](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/README.md#L34-L48).

**Reasoning.** The tool can reverse-complement, assemble, design, compress, and analyze
sequences, but none of those operations permutes an input sequence. `revcomp` is
explicitly excluded by the operational boundary.

**Disconfirmation attempt.** I searched the complete source, API guide, user guide,
examples, and command registry for `shuffle`, `scramble`, `dinucleotide`, `permute`,
and sequence randomization. I followed all ambiguous hits: Fisher–Yates shuffles an
internal degenerate-design action space; `np.random.shuffle(indexvec)` randomizes read
counting order; “Dinucleotide Match” is a low-complexity read classifier. None emits a
shuffled variant of an input sequence. A registered API/CLI operation doing so would
have changed the value.

### Mutation Maker — `no` (medium-high confidence)

**Quoted primary-source evidence.** The complete FastAPI design surface exposes
`POST /v1/ssm`, `POST /v1/qclm`, and `POST /v1/pas`; the frontend's complete workflow
enum likewise contains only `ssm`, `qclm`, and `pas`. The backend tasks dispatch these
to `ssm_solve`, `qclm_solve`, and `pas_solve`.

**Location.** Surviving Mutation Maker source, commit
[`396c1c0ede7529f3dedf65e821e8c1f20c9a7043`](https://github.com/ra100/Mutation_Maker/commit/396c1c0ede7529f3dedf65e821e8c1f20c9a7043), API version 1.0.0:
[`api/server_fastapi.py`, lines 36–65](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/api/server_fastapi.py#L36-L65),
[`backend/tasks.py`, lines 23–73](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/tasks.py#L23-L73), and
[`frontend/src/shared/workflow.ts`, lines 19–23](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/frontend/src/shared/workflow.ts#L19-L23).

**Reasoning.** SSM/QCLM/PAS design mutagenic primers and gene-synthesis oligos. None is
an operation that shuffles an input sequence. Random or degenerate mutagenesis is not a
shuffle merely because its outputs are stochastic or diverse.

**Disconfirmation attempt.** I searched the full backend, API, frontend, tests/features,
and README for `shuffle`, `shuffling`, `scramble`, `dinucleotide`, `permutation`, and
input-sequence randomization, and inspected the complete workflow/API dispatch surface.
There were no relevant hits. A fourth workflow, API parameter, or solver operation that
permuted residues would have changed the value.

**Evidence limitation.** The paper's original `Merck/Mutation_Maker` repository URL is
no longer available. The inspected `ra100/Mutation_Maker` repository is the surviving
source with the original Merck history plus later maintenance. That provenance gap
reduces confidence slightly, but it does not make the result `unknown`: the complete
surviving implementation and its public API are determinate and contain no shuffle.

### DNA Chisel — `no` (high confidence)

**Quoted primary-source evidence.** The complete package `__all__` exports its
optimization problem, built-in specifications, sequence patterns, and biotools. It
includes `random_dna_sequence`, `reverse_complement`, and `reverse_translate`, but no
shuffle or scramble operation. The only `np.random.shuffle` source hit is
`np.random.shuffle(locations)` inside evaluation of codon-usage optimization.

**Location.** DNA Chisel v3.2.16, commit
[`68c09304341c3656f3dfe63eda37757d6a7b3917`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/commit/68c09304341c3656f3dfe63eda37757d6a7b3917):
[`dnachisel/__init__.py`, lines 1–68 and 75–132](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/__init__.py#L1-L132),
[`dnachisel/biotools/__init__.py`, lines 16–46 and 58–93](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/biotools/__init__.py#L16-L93), and
[`MatchTargetCodonUsage.py`, lines 125–135](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/builtin_specifications/codon_optimization/MatchTargetCodonUsage.py#L125-L135).

**Reasoning.** Shuffling the order in which an optimizer considers violation locations
does not shuffle sequence residues. Random sequence generation is also explicitly not a
shuffle of an input. A user might write a custom specification or their own permutation
logic, but the global rule forbids credit for user reconstruction.

**Disconfirmation attempt.** I searched the complete source tree, built-in
specification inventory, biotools API, examples, and first-party docs for `shuffle`,
`scramble`, `dinucleotide`, `permutation`, and sequence randomization. I inspected the
sole actual shuffle call and confirmed its target is a Python list of locations, not a
sequence. A built-in specification or biotool accepting an input sequence and returning
its residue permutation would have changed the value.

### tangermeme — `yes` (high confidence)

**Quoted primary-source evidence.** Public `dinucleotide_shuffle` is documented as:
“Given a one-hot encoded sequence, dinucleotide shuffle it.” It accepts `start`, `end`,
`n`, and `random_state`, constructs an Eulerian path, and returns shuffled sequences.
The composition test counts every adjacent pair in original and output and asserts
`dinucs[key] == value`.

**Location.** tangermeme v1.4.1, commit
[`2006b310cd72a28c56c3ea4ba67f738fff74bb89`](https://github.com/jmschrei/tangermeme/commit/2006b310cd72a28c56c3ea4ba67f738fff74bb89):
[`tangermeme/ersatz.py`, lines 529–605 and 608–647](https://github.com/jmschrei/tangermeme/blob/2006b310cd72a28c56c3ea4ba67f738fff74bb89/tangermeme/ersatz.py#L529-L647),
[`README.md`, lines 114–140](https://github.com/jmschrei/tangermeme/blob/2006b310cd72a28c56c3ea4ba67f738fff74bb89/README.md#L114-L140), and
[`tests/test_ersatz.py`, lines 651–711](https://github.com/jmschrei/tangermeme/blob/2006b310cd72a28c56c3ea4ba67f738fff74bb89/tests/test_ersatz.py#L651-L711).

**Reasoning.** This is a documented public operation that performs the shuffle itself,
accepts input sequences, can generate multiple shuffled variants, and preserves exact
dinucleotide counts. It meets `yes` directly.

**Disconfirmation attempt.** I searched the public API, README/tutorial material,
implementation, and tests for both ordinary and dinucleotide shuffling, then checked
whether the latter was only a private reference generator or example. It is a public
function, is directly documented, and has exact-count and regional-shuffle tests.
Evidence that only the ordinary `shuffle` function were public would have reduced the
value to `partial`; evidence that callers had to supply the permutation would have made
it `no`.

## Source and version notes

- Scores describe the source snapshots above, not capabilities inferred from earlier
  survey matrices, scored files, reviews, or analyses. Those materials were not used as
  evidence.
- Full-tree searches included source, first-party documentation, public API/command
  registries, examples, and tests where present. Generated binaries and `.git` metadata
  were excluded.
- MPRAnator has no tagged/package version in its one-commit repository; the commit hash
  is therefore the precise version identifier.
- MPRA Design Tools comprises the v0.2.0 R package and its companion Shiny repository;
  both were checked because either could have supplied the operation.
- Mutation Maker's original paper-cited repository is unavailable; the surviving source
  provenance limitation is recorded in its section rather than converted into an
  `unknown` score.
