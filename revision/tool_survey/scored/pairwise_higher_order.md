# Audit: pairwise_higher_order

Audited 2026-08-15. This audit was performed independently from the row
definition and primary sources. The files in `tool_survey/final/` were used only
as locators; no statement or score from those files is evidence here.

## Measurement instrument and operational test

The binding definition in `revision/tables/ROW_DEFINITIONS.md` is: “Two or more
**independent** mutation events co-occurring in a single output sequence, with
the combinations enumerated or sampled **by the tool**.” It excludes
hand-authored multi-variant inputs, uniform context edits, a multi-base single
event, repeated-run merging, multi-motif placement, and IUPAC expansion; an
objective-driven optimizer can score `partial` at most.

**Operational test stated before inspecting any tool:** Given one parent
sequence and a design specification declaring at least two independent
substitution, insertion, or deletion events, can the documented tool itself
enumerate or sample outputs in which multiple events co-occur—excluding
user-authored multi-variant inputs, uniform context edits, multi-base single
events, repeated-run merging, multi-motif placement, and IUPAC expansion, while
optimizer-found multi-edit sequences can score only partial?

I applied that sentence unchanged to all eight tools. A VCF containing a list of
candidate single events does not by itself pass. It passes only when the tool,
rather than the VCF author, constructs combinations of those events (the
distinction that matters for MPRAnator versus VaLiAnT and MPRA Design Tools).

## Results

| Tool | Score | Confidence | Decisive basis |
|---|---:|---:|---|
| PoolParty | **yes** | high | Public mutagenesis operation exhaustively enumerates fixed-order combinations and samples variable-order mutants. |
| VaLiAnT | **no** | high | Generated mutators emit a flat union of one-event variants; multi-edit context comes only from excluded background/PAM edits or user VCFs. |
| MPRAnator | **yes** | high | Documented SNP-combination option enumerates every nonempty subset and applies all events in each subset. |
| MPRA Design Tools | **no** | high | Processes the VCF row-by-row and emits reference/alternate replicates for one input variant at a time. |
| Oligopool Calculator | **partial** | high | A bundled example helper enumerates multi-position variants, but it is absent from the documented package API. |
| Mutation Maker | **yes** | high | QCLM enumerates Cartesian products of codons across mutation sites into concrete multi-site primer sequences. |
| DNA Chisel | **yes** | high | Documented `MutationSpace.all_variants()` enumerates the Cartesian product of independent segment choices. |
| tangermeme | **partial** | high | Beam search returns model-selected multi-edit sequences; the binding optimizer rule caps this at partial. |

## Tool audits

### PoolParty — yes (high confidence)

**Positive evidence.** `mutagenize` documents `num_mutations` as a “Fixed
number of mutations to apply” and exposes both random and sequential selection
modes (`poolparty/src/poolparty/base_ops/mutagenize.py:43-56`, commit
`1bb0179e1c3720b1fffd471802b3040f9336de28`). The sequential implementation
computes

> `comb(num_positions, self.num_mutations) * (alpha_minus_1**self.num_mutations)`

and loops over

> `for positions in combinations(range(num_positions), self.num_mutations):`

before storing each position/substitution tuple
(`src/poolparty/base_ops/mutagenize.py:332-381`). This is tool-generated
enumeration, not template decompression.

The public DMS tutorial is explicit:

> “The same operation with `num_mutations=2` enumerates every possible pair of
> single amino acid changes.”

It then calls `mutagenize_orf(num_mutations=2, mode="sequential")` and reports
`C(55, 2) x 19^2 = 536,085` combinations
(`poolparty/docs/tutorials/dms_gb1.rst:96-116`). The same page documents random
higher-order mutants: a per-codon `mutation_rate` makes “each sequence” receive
a variable number of changes, with `num_states=10000` random draws
(`dms_gb1.rst:118-140`).

**Behavioral check.** Required read-only execution used the repository venv and
`PYTHONDONTWRITEBYTECODE=1`:

```text
wt = 'ACGT'
p = pp.from_seq(wt).mutagenize(num_mutations=2, mode='sequential')
states 54
rows 54
hamming_counts [2]
first rows ['CAGT', 'CGGT', 'CTGT', 'GAGT']
```

The count is `C(4,2) * 3^2 = 54`; every output differs from the parent at
exactly two positions.

**Disconfirmation attempt.** I searched the public operation docs, DMS
tutorial, and complete `base_ops/mutagenize.py` for `num_mutations`,
`mutation_rate`, `combinations`, `sequential`, and restrictions that would make
the feature motif placement, IUPAC expansion, user enumeration, or example-only.
The operation is in the documented API and the implementation itself generates
the combinations. A different value would have resulted if the two-change
outputs existed only in the tutorial, if the caller supplied every combined
sequence, or if only an optimizer found a single multi-edit sequence; those
possibilities were searched and contradicted by the operation source and live
run.

### VaLiAnT — no (high confidence)

**Evidence of the implemented scope.** The manual's mutation-type list contains
parametric deletion, SNV, in-frame deletion, alanine/stop/all-amino-acid codon
substitution, and SNVRE (`README.md:255-270`, commit
`8796cc112dafd4919fec59913f58cd2be87c45eb`). Its SNV example is explicitly a
single-event scan:

> “Each nucleotide is replaced with all the alternatives.”

and for `AA` it lists `CA, GA, TA, AC, AG, AT`, never a double mutant
(`README.md:294-307`). `SnvMutator.get_variants` constructs a list of individual
`Variant.get_snvs(...)` values (`src/valiant/mutators/snv.py:44-49`), while
`DeletionMutator.get_variants` returns one `Variant.get_del(...)` per reference
span (`mutators/deletion.py:43-48`). Most decisively, a collection of requested
mutators is flattened as

```python
return [
    PatternVariant.from_variant(m.as_str(), v)
    for m in mutators
    for v in m.get_variants(seq)
]
```

(`src/valiant/mutator.py:133-138`). That is a union of one-event outputs, not a
product over events.

VaLiAnT can apply background and PAM-protection variants to all targeton oligos,
but those are uniform context edits, explicitly excluded by the row. It can also
ingest user variants; the manual says custom variants are “Applied to the
targeton reference sequence as a whole” and supports substitutions, insertions,
deletions, and indels (`README.md:465-476`). The paper says such variants
generate “one oligonucleotide sequence each” (Barbon et al. 2022, §2.5.2). This
is the excluded user-VCF route, not tool-generated event combination.

**Disconfirmation attempt / absence search.** I searched all of `README.md`,
`CHANGELOG.md`, `src/valiant/`, `docs/`, and `examples/`, plus the full paper,
for `pairwise`, `higher[- ]order`, `multiple mutations`, `combinatorial`,
`mutations per`, `num.*mut`, `itertools.combinations`, `itertools.product`,
`mutation_rate`, and `num_mutations`. The combination/product and higher-order
mechanism searches returned no hits. I then read the generator path
`mutator.py`, `mutators/snv.py`, `mutators/deletion.py`, `custom_variant.py`, and
the README's mutation/custom-variant/background/PAM sections. A different score
would have required a parameter or implementation that takes two candidate
events and emits/samples their co-occurrence; I explicitly searched the source,
manual, examples, changelog, and paper for it and found only the flat one-event
union and excluded context/input routes.

### MPRAnator — yes (high confidence)

**Positive evidence.** The repository's user documentation says:

> “The user can select to include or exclude combinations of SNPs when
> designing MPRA experiments.”

(`iliasApp/templates/iliasApp/docs.html:40-43`, commit
`9969790d62410138d4281b7955da6d085f07b1bc`). The form exposes the supported
control `makeSnpCombinations = forms.BooleanField(label="Make Snp
Combinations?")` (`iliasApp/forms.py:501-505`), and the view sends it to
`make_sequence_copies` (`iliasApp/views.py:342,351-362`). The paper independently
describes examining “the functional effects of single or combinations of
single-nucleotide polymorphisms” (Georgakopoulos-Soares et al. 2017, abstract).

The source establishes who performs the combination. `generateCombinations`
loops subset sizes 1 through N and calls `itertools.combinations`
(`parseSNPs.py:20-26`). When the option is true, it groups candidate events by
parent sequence, calls that enumerator, then loops `for snped in comb`. For each
event it iterates over the current `Sequenced` list and appends an edited copy,
so the list expands from parent to singles and co-occurring edits
(`parseSNPs.py:36-68,84-123`). There is an implementation quirk worth recording:
the final `ExtraSequences += Sequenced` is indented outside `for comb`
(`parseSNPs.py:124-132`), so earlier subset iterations are discarded. However,
`generateCombinations` orders the full N-event subset last, and the evolving-list
expansion of that final subset itself contains all subcombinations, including
multi-event sequences. Thus the input rows specify candidate single events, but
the user does not pre-author the multi-event outputs: MPRAnator constructs them.

**Disconfirmation attempt.** I searched every repository Python/template file
and the paper for `SNP`, `combination`, `makeSnpCombinations`, `product`,
`multiple mutation`, and `Transmutation`; I read the form, view, complete
combination branch (including its indentation/control flow), and no-combination
branch (`parseSNPs.py:134-185`). A
different value would have followed if “combinations” meant motif placement, if
the VCF already contained the combined haplotypes, or if the checkbox merely
processed each row independently. The source instead generates nonempty subsets
and applies all `snped` members of a subset, while the documented form makes it
a supported user path rather than an example-only script.

### MPRA Design Tools — no (high confidence)

This score covers both primary repositories named for the tool:
`andrewGhazi/mpradesigntools` commit
`afd386ef12051bb0a37ad63a6f92acd555246584` and
`andrewGhazi/designMPRA` commit
`0cf56ee602fc86dde705906d071a39cbdf6e99a8`.

**Evidence of the implemented scope.** The R package documentation says
`processVCF` “takes a VCF of SNPs ... and turns them into a set of labeled MPRA
sequences barcoded with inert n-mers” and defines `nper` as the number of
barcoded sequences “per allele per SNP” (`man/processVCF.Rd:28-31,78-87`). The
core helper is documented as:

> “process an individual SNP”

and its `snp` parameter is “a tibble containing the VCF information for one
SNP” (`R/processVCFfast.R:187-204`). It classifies that one row as SNV,
insertion, or deletion (`:210-242`); for an SNV it constructs exactly one
`altseq` with `replaceLetterAt` and emits replicated `ref`/`alt` rows
(`:326-399`). At the outer level, `processVCF` explicitly does

```r
processed = vcf %>%
    rowwise %>%
    do(seqs = processSnp(., ...))
```

(`R/processVCFfast.R:1221-1237`). The paper likewise describes “reference and
alternate alleles of variants” and VCF-driven construct generation (Ghazi et
al. 2018, pp. 1-2), not event products.

The optional `alter_aberrant` branch repairs restriction sites in the context
(`R/processVCFfast.R:347-353`), which is a context/constraint repair rather than
enumeration of declared mutation combinations and cannot rescue this row.

**Disconfirmation attempt / absence search.** Across both repositories I
searched `R/`, `man/`, both READMEs, all workflow scripts, `server.R`, `ui.R`,
and the paper for `pairwise`, `higher[- ]order`, `haplotype`, `multiple
mutation`, `combinations`, `combn(`, `expand.grid`, `crossing`, `product(`,
`num_mutations`, and `mutation_rate`. The only `crossing` hits were commented
restriction-site bookkeeping; no event-combination generator was found. I read
`processVCF`, `processSnp`, all three variant-type construction branches, the Rd
page, and the Shiny VCF path. A different value would require a supported option
that groups two VCF rows or creates an event product before construct emission;
I searched for exactly those mechanisms and the source instead enforces
row-by-row one-variant processing.

### Oligopool Calculator — partial (high confidence)

**Qualifying evidence, but only in an example.** The current repository includes
`examples/library-compressor/mutant_generator.py` (commit
`b88fa394ca67ed4c48ec41127b5707694ee7cf0a`). Its module docstring says it
includes “Combinatorial variants at multiple positions” (`:1-10`). The helper
`generate_multi_position_variants(sequence, positions, ...)` says:

> “Generate all combinatorial variants at specified positions.”
>
> “Each position can be A, C, G, or T, giving 4^len(positions) variants.”

It calls `itertools.product(DNA_BASES, repeat=len(positions))`, applies every
base in the tuple to the named positions, and yields each concrete sequence
(`mutant_generator.py:129-176`). A read-only execution on parent `AAAA`,
positions `[1,3]`, and `include_wildtype=False` returned 15 sequences with
observed Hamming distances `[1, 2]`; therefore concrete double mutants are among
the outputs. The example README explicitly labels “Multi-position Variants” and
says it “Generates variants with mutations at multiple positions”
(`examples/library-compressor/README.md:44-47`).

**Why not yes.** This helper is outside `oligopool/`; it is not imported or
advertised by the package. The public API list contains barcode, primer, motif,
spacer, assembly/QC, `compress`, `expand`, and counting operations, but no
mutant generator (`oligopool/__init__.py:6-31,42-70`). The main docs describe the
workflow as “Design variants computationally” and then compress them
(`docs/docs.md:889-908`); `compress` accepts already concrete input variants,
while `expand` expands IUPAC sequences. The latter is explicitly excluded by
the row and is not the reason for partial credit.

**Disconfirmation attempt.** I searched all package modules, `README.md`,
`docs/api.md`, `docs/docs.md`, `docs/agent-skills.md`, examples, changelog, and
the paper for `generate_multi_position_variants`, `mutant_generator`,
`pairwise`, `higher-order`, `multiple mutation`, `combinatorial variant`,
`mutation`, `mutagen`, `degenerate`, and `IUPAC`. The qualifying symbol occurs
only in the example helper; searches of package and main API/docs return no
public mutagenesis operation. A `yes` would require the same generator in the
documented API; a `no` would require absence even from examples. Both
alternatives were explicitly checked. The binding “example script” rule
therefore fixes the score at `partial`.

### Mutation Maker — yes (high confidence)

**Positive evidence.** The QCLM/MSDM workflow is not merely multiple separate
primers. The paper states that Mutation Maker outputs primers that “cover one or
more target sites” and that the resulting primers yield a “diverse
combinatorial library” (Hiraga et al. 2021, Fig. 1 caption). It also states that
the algorithm groups neighboring target sites covered by one primer and
“combines selected codons that encode for the desired amino acids” (Results,
MSDM section).

The current source makes the concrete enumeration explicit. The supported QCLM
configuration exposes `use_degeneracy_codon` (`backend/mutation_maker/qclm_types.py:45-53`,
commit `396c1c0ede7529f3dedf65e821e8c1f20c9a7043`), and the UI labels the switch
“Use Degeneracy Codon” (`frontend/src/scenes/QCLM/components/QCLMForm.tsx:308-317`).
When degeneracy is disabled (or no valid degenerate cover is found), QCLM picks
ordinary codons (`qclm.py:315-341`). For each multi-site sequence it builds
codon choices and calls `create_subsets_for_primers` (`qclm.py:398-413`). That
function's contract is decisive:

> “The number of items in the result is always the number of all possible
> combinations of mutations on all sites at the single primer. This
> combination is found as cartesian product of all possible mutations.”

The implementation is `itertools.product(*degenerate_codons_lists)` and emits
one codon list per combination
(`backend/mutation_maker/degenerate_codon.py:334-379`). `create_new_output` then
renders each primer with `primer.spec.get_sequence(...)` and returns its
sequence plus its list of mutations/codons (`qclm.py:597-668`). A single
concrete primer can therefore carry two or more independently chosen codon
substitutions, and the tool enumerates their product.

**Disconfirmation attempt.** I searched the paper, README, backend workflows,
QCLM/PAS types and solutions, API routes, frontend form, and examples/tests for
`multi-site`, `combinatorial`, `multiple mutation`, `mutation sites`,
`degenerate`, `codon combination`, `itertools.product`, and output sequence
construction. A different score would have resulted if the tool only output an
ambiguous IUPAC primer (excluded), only relied on wet-lab stochastic assembly,
or merely designed one primer per site. The non-degenerate switch/fallback,
Cartesian-product implementation, multi-site grouping, and concrete
`PrimerOutput.sequence` contradict all three alternatives.

### DNA Chisel — yes (high confidence)

**Positive evidence.** `MutationSpace` is a documented public core class: the
API reference uses `automodule:: dnachisel.MutationSpace` with `:members:`
(`docs/ref/core_classes.rst:45-49`, commit
`68c09304341c3656f3dfe63eda37757d6a7b3917`). Its class example constructs two
independent `MutationChoice` segments with their own variant sets
(`dnachisel/MutationSpace/MutationSpace.py:9-31`). `MutationChoice` and
`MutationSpace` are both public exports
(`dnachisel/MutationSpace/__init__.py:1-4`).

The enumerator is direct rather than optimizer-mediated:

> `def all_variants(self, sequence):`
>
> `    """Iterate through all sequence variants in this mutation space."""`

It creates one variant slot for every independent multichoice, loops over
`itertools.product(*variants_slots)`, applies every selected segment to the
same output buffer, and yields the concrete sequence
(`MutationSpace.py:132-164`). A mutation space can also be derived from a
declared optimization problem; `from_optimization_problem` starts with
per-position nucleotide choices and intersects them with the constraints'
declared choices (`MutationSpace.py:166-208`). The paper corroborates that the
solver can use “an exhaustive search through all possible sequence variants”
(Zulkower and Rosser 2020, §3), but full credit here comes from the separately
callable enumerator, not from an optimizer's final sequence.

**Disconfirmation attempt.** I searched the entire package/docs/paper for
`MutationSpace`, `MutationChoice`, `all_variants`, `pick_random_mutations`,
`apply_random_mutations`, `product`, `combinatorial`, `higher-order`, and
`multiple mutations`, then read the full `MutationSpace.py`, public exports,
core API page, and solver call sites. A different score would have followed if
only the objective solver returned one converged multi-edit sequence (which
would be partial), or if `all_variants` only decompressed IUPAC input. Instead,
the documented API enumerates the product of independently declared segment
choices and yields concrete sequences, so neither limiting condition holds.

### tangermeme — partial (high confidence)

**Qualifying optimizer evidence.** Current tangermeme documents
`beam_substitution`, which starts from one complete parent sequence and expands
each beam member every round by substituting each candidate at every allowed
position. Its docstring says larger beams

> “can recover good multi-edit combinations that the greedy method prunes
> away”

(`tangermeme/design/beam_substitution.py:44-67`, commit
`2006b310cd72a28c56c3ea4ba67f738fff74bb89`). With `motifs=None`, the candidates
are the alphabet itself, so each step changes one character rather than placing
a multi-base motif (`beam_substitution.py:97-101,192-195`). The implementation
tiles edits for every current beam sequence and carries the selected edited
sequence into the next round (`:218-299`). `n_best` returns one or more complete
designed sequences (`:125-135,173-176`). The changelog describes the same
feature as recovery of “good multi-edit combinations” (`docs/whats_new.rst:34-46`).

This meets the factual co-occurrence condition, but the edits are selected by
model loss. The row's binding boundary says objective-driven optimizers score
`partial` at most; therefore this is `partial`, not `yes`.

**Disconfirmation attempt.** I searched all 28 Python modules, all docs and
bundled skill references, README, changelog, and paper for `pairwise`,
`higher-order`, `combin`, `multiple mutation`, `n_edits`, `max_iter`,
`mutation`, `variant`, and `substitution`. I specifically read the three
plausible alternatives. (1) `saturation_mutagenesis` promises sequences “with
an edit distance of one” and constructs only true single-base edits
(`saturation_mutagenesis.py:81-89,212-232`). (2) `variant_effect.substitution_effect`
can apply multiple rows at once, but the caller supplies every row; its docs say
the tensor rows are the individual variants (`variant_effect.py:29-74`), so it
is excluded hand-authored input. (3) `space` concerns motif spacing and is
excluded placement. A `yes` would require a declarable, non-objective
enumerator/sampler of mutation-event combinations; I searched every design,
ISM, variant-effect, ersatz, product, and spacing entry point and found none. A
`no` would require absence of tool-generated multi-edit output; beam search
contradicts that. Hence `partial` is forced by the optimizer boundary.

## Search and verification log

The following are the searches actually used for the audit (case-insensitive
unless shown otherwise). Repository-wide searches excluded binary images and
PDFs; papers were converted with PyMuPDF exactly as instructed and searched as
plain text.

- **Definition:** read all of
  `revision/tables/ROW_DEFINITIONS.md` before any tool inspection.
- **Repository acquisition:** shallow read-only clones in `/tmp` of the seven
  public repositories/branches named in the task; PoolParty used the supplied
  read-only local checkout. Recorded the full commit hashes shown above.
- **Paper extraction:** `python3 -c "import fitz; ... p.get_text() ..."` for all
  seven supplied public-tool PDFs and the supplied PoolParty PDF. No PDF quote
  was taken from a secondary analysis.
- **PoolParty:** searched `docs/tutorials/dms_gb1.rst`,
  `docs/operations/mutagenize.rst`, and complete `base_ops/mutagenize.py` for
  `num_mutations|combinations|mutation_rate|higher.order|double`; read source
  lines 320-390 and tutorial lines 70-150; executed the 4-nt double-mutant test.
- **VaLiAnT:** repository/manual/paper search for
  `pairwise|higher[- ]order|multiple mutations|combinatorial|mutations per|num.*mut|custom variant|VCF|background|pam protection|mutator`;
  a second source search for
  `itertools.(combinations|product)|num_mutations|mutation_rate|higher.order|pairwise`;
  read `mutator.py`, SNV/deletion mutators, mutation/custom-variant README
  sections, changelog, examples, and relevant paper sections.
- **MPRAnator:** full-tree search for
  `pairwise|higher[- ]order|combin|multiple.*(mut|snp|variant)|mutation|SNP|motif|product`;
  focused paper search for `combinations of SNP|SNP combinations|combination`;
  read `parseSNPs.py`, the form, view, and user-doc template.
- **MPRA Design Tools:** both repositories searched for
  `pairwise|higher[- ]order|combin|multiple.*(mut|variant|snp)|mutation|SNP|variant|oligo`;
  second negative search for
  `itertools|expand.grid|crossing|combn|combinations|product|num_mutations|mutation_rate|higher.order|pairwise|haplotype`;
  read `processVCFfast.R`, its Rd page, both READMEs, Shiny server/UI, scripts,
  and the paper.
- **Oligopool Calculator:** full tree/docs/paper searched for
  `pairwise|higher[- ]order|combin|multiple.*(mut|variant|edit)|mutation|mutagen|variant|degenerate|IUPAC|VCF`;
  exact-symbol search for `generate_multi_position_variants|mutant_generator`
  over package, main docs, pyproject, and example index; read the example helper,
  example README, package exports, main degenerate API/docs, and executed the
  helper on a two-position parent.
- **Mutation Maker:** full tree/paper search for
  `pairwise|higher[- ]order|combin|multiple.*(mut|variant|edit)|multi[- ]site|number.*mutation|mutation|variant|degenerate|IUPAC|NNK`;
  focused searches for `QCLM|PAS|mutations.*primer|sites.*primer|codon combination|itertools.product|use_degeneracy_codon`;
  read QCLM solver/types/solution, degenerate-codon product code, API route,
  frontend selector, README, tests, and the MSDM/PAS paper sections.
- **DNA Chisel:** full package/docs/paper search for
  `all_variants|apply_random_mutations|pick_random_mutations|MutationSpace|mutation space|combinatorial|higher.order|pairwise|multiple mutations|variants`;
  read the complete `MutationSpace.py`, `MutationChoice.py`, public exports,
  core API page, and solver call sites.
- **tangermeme:** all Python modules/docs/README/paper searched for
  `pairwise|higher[- ]order|combin|multiple.*(mut|variant|edit)|n_edits|num.*edit|rounds|mutation|mutagen|variant|substitution`;
  read complete beam-search control flow and returns, variant-effect mutation
  application, saturation-mutagenesis construction, spacing/product candidates,
  public design docs, bundled design/variant references, and changelog.

## Row-level finding

The row works on one consistent scale and does discriminate: four `yes`, two
`partial`, and two `no`. It distinguishes (a) tool-generated explicit event
products/samples, (b) limited implementations confined to examples or
objective-driven search, and (c) flat one-event processing. The two boundary
rules that materially affect the result are valuable rather than fatal:
Oligopool Calculator's example-only generator is capped at `partial`, and
tangermeme's model-guided multi-edit search is capped at `partial`. Motif
placement, degenerate expansion, VCF pass-through, and uniform repair edits were
consistently denied credit across all tools.
