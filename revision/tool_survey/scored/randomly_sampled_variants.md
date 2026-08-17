# Audit: randomly sampled variants

Date: 2026-08-15

## Measurement instrument and operational test

The binding row definition in `25_poolparty_private/revision/tables/ROW_DEFINITIONS.md` is:

> A documented rate, count or RNG parameter produces a *sample* of the variant space. **Random sequence generation is not sampled mutagenesis** — generating random barcodes, shuffled controls, or stochastic search during optimisation does not count.

**Operational test applied unchanged to all eight tools:** look for a documented library-design API in which a user-supplied mutation rate, requested sample count, mutation-event count, or RNG/seed parameter causes the tool to return a non-exhaustive sample of variants derived from one or more parent sequences; reject de novo random sequence or random element generation, shuffled controls, and candidates visited only inside an optimization search.

The term “count” therefore includes a count of mutation events (for example, exactly three randomly positioned substitutions), not only a requested number of output members. A one-member random draw is still a sample. A low-level sampler is credited if it is included in the tool's documented API and directly returns the sampled variant; merely observing random candidates inside an optimizer is not.

## Results

| Tool | Score | Confidence |
|---|---|---|
| PoolParty | **yes** | high |
| VaLiAnT | **no** | high |
| MPRAnator | **yes** | high |
| mpradesigntools / designMPRA | **partial** | medium |
| Oligopool Calculator | **no** | high |
| Mutation Maker | **no** | high |
| DNA Chisel | **yes** | medium |
| tangermeme | **no** | medium-high |

The row discriminates: three tools pass directly, one has a narrowly restricted implementation, and four lack qualifying sampled mutagenesis. The two most important boundary decisions are documented below: DNA Chisel's directly callable and documented `MutationSpace.apply_random_mutations` is credited separately from its optimizer loop, while tangermeme's random replacement backgrounds are excluded as random sequence generation rather than mutation-position sampling.

## Source acquisition and common search protocol

Primary repositories were inspected at these revisions:

- PoolParty local repository `poolparty-statecounter`, commit `1bb0179`.
- VaLiAnT `develop`, commit `8796cc1`.
- MPRAnator `master`, commit `9969790`.
- mpradesigntools `master`, commit `afd386e`, and designMPRA `master`, commit `0cf56ee`.
- Oligopool Calculator `master`, commit `b88fa39`.
- Mutation Maker `master`, commit `396c1c0`.
- DNA Chisel `master`, commit `68c0930`.
- tangermeme `main`, commit `2006b31`.

All remote repositories were shallow-cloned from the GitHub URLs specified in the task. For each tool I recursively searched source and documentation with combinations of:

```text
rg -n -i '(mutation[_ -]?rate|num[_ -]?(samples|variants|mutations)|n[_ -]?(samples|variants|mutations)|random[_ -]?(mut|variant|sample)|sample[_ -]?(mut|variant)|randomly|seed|rng)'
rg -n -i '(random|sample|seed|rng|rate|frequency|mutagen|mutation|variant|shuffle)'
rg -n '^def |^class '
```

Large generated data, locks, SVGs, and notebooks were excluded from the first source-code pass where appropriate. Notebook source was then searched separately when it was first-party documentation. Each local paper was converted to text, page by page, with the required PyMuPDF procedure:

```text
python3 -c "import fitz; d=fitz.open('PATH'); print(chr(10).join(p.get_text() for p in d))"
```

and searched for `random`, `sample`, `seed`, `rng`, `rate`, `mutation`, `mutagenesis`, and tool-specific near matches. The `final/<slug>.md` records were searched only as navigation leads before repository inspection; no statement, score, or citation from those records is evidence here.

## Per-tool audits

### PoolParty — **yes** (high confidence)

**Qualifying evidence.** The public documentation is explicit:

> “`mode="random"` draws a possibility from the design space each time the library is generated.”

and:

> “A single random mutant is drawn. To generate multiple fixed random designs, use the `num_states` parameter.”

Source: local primary repository, `poolparty/docs/operations/modes.rst:86-110`.

The same page gives the concrete multi-draw API:

```python
wt      = pp.from_seq("ACG")
mutants = pp.mutagenize(wt, num_mutations=1, mode="random", num_states=5)
```

and states:

> “Five randomly chosen single-point mutants are generated. Duplicates are possible because each draw is independent.”

Source: `poolparty/docs/operations/modes.rst:210-234`.

The `mutagenize` API also documents `mutation_rate` as the “Probability of mutation at each position,” `mode` as random or sequential, and `num_states` as defaulting to one in random mode (“pure random sampling”). Source: `poolparty/src/poolparty/base_ops/mutagenize.py:35-70`. `generate_library` documents `num_seqs` and `seed` (`poolparty/src/poolparty/generate_library.py:15-43`), and the implementation creates a per-operation NumPy generator from the master seed, operation ID, and state (`generate_library.py:275-284`).

The authors' paper independently describes the same mechanism:

> “In random mode, the Operation samples sequence designs stochastically. For example, when mutagenize operates in random mode, it randomly selects which positions to mutate and which mutations to introduce.”

Source: `25_poolparty_private/26.05.18_bmc_submission/latex/main.pdf`, p. 5. The repository is treated as authoritative, as required.

**Behavioral verification.** Using the required read-only venv and `PYTHONDONTWRITEBYTECODE=1`, I ran `mutagenize('AAA', num_mutations=1, mode='random', num_states=5)` three times. Seed 7 returned `['AAC', 'GAA', 'AAT', 'ACA', 'AAT']`; a second seed-7 construction returned the same list, and seed 8 differed. Five rows were drawn from the nine-member single-substitution space. This settles that `num_states` controls sampled variants and `seed` controls reproducibility, rather than randomness being incidental to some other operation.

**Disconfirmation attempt / what would change the value.** I searched the complete docs, README, `base_ops/mutagenize.py`, `orf_ops/mutagenize_orf.py`, and `generate_library.py` for `random`, `mutation_rate`, `num_states`, `num_seqs`, `seed`, `sample`, and `RNG`; I also ran the API. Evidence that random mode generated unrelated de novo sequences, that `num_states` merely truncated a deterministic enumeration, or that the RNG existed only inside an optimizer would have changed this to `no` or `partial`. The quoted docs, implementation, and run rule all three out.

### VaLiAnT — **no** (high confidence)

**Evidence of the implemented behavior.** The documented SNV operation is exhaustive:

> “Each nucleotide is replaced with all the alternatives.”

The worked `AA` example lists all six single substitutions. Source: `README.md:294-307` at `cancerit/VaLiAnT@8796cc1`.

The implementation agrees. `SnvMutator.get_variants` flattens `Variant.get_snvs(...)` for every reference position (`src/valiant/mutators/snv.py:44-49`); `DeletionMutator` returns a deletion for every reference segment (`mutators/deletion.py:43-48`); amino-acid substitution nests over every codon reference and every eligible alternative codon (`mutators/codon.py:95-103`). `MutatorCollection._get_variants` then returns every variant from every selected mutator via nested list comprehensions (`src/valiant/mutator.py:133-138`). There is no sampling stage.

The paper likewise describes the operations as saturation: nucleotide-level functions perform saturation and the `aa` function “exchanges each wild-type codon for the most frequent triplet code of each other amino acid.” Source: `Barbon2022_VaLiAnT_all.pdf`, p. 2.

**Disconfirmation attempt / what would change the value.** I searched all 76 Python files under `src/valiant`, the complete 875-line `README.md`, `CHANGELOG.md`, examples, tests, and every page of the paper for `random`, `sample`, `seed`, `rng`, `rate`, `mutation`, `variant`, and mutator parameters. The only source hit resembling `rng` is a local variable meaning genomic/codon *range* (`targeton.py:189`); there is no Python random/NumPy RNG import in the design code, no seed option, no mutation-rate option, and no sample-count option. A documented rate/count/seed on any mutator, or a sampling call between the exhaustive mutators and output, would have produced `yes` (or `partial` if restricted to one mutator). None was found. This `no` means the source and documented interface were positively checked, not that documentation was unavailable.

### MPRAnator — **yes** (high confidence)

**Qualifying evidence.** MPRAnator's first-party web documentation distinguishes mutation from shuffling:

> “Mutate random positions in the input sequence.”

and describes the output as:

> “`Mutated_nucleotides - 3` This is the number of randomly chosen mutated nucleotides.”

Source: `iliasApp/templates/iliasApp/docs.html:138-159` at `hemberg-lab/MPRAnator@9969790`.

The documented form exposes the count as `numToMutate = PositiveIntegerField(label="Number of nucleotides to mutate", initial=0)` (`iliasApp/forms.py:340-377`). The request handler passes that value to `part3.mutateString` for each input parent sequence and returns the mutated sequences in FASTA (`iliasApp/views.py:170-224`). The home-page description is unusually direct:

> “Creates a specified number of random mutations on a sequence”

Source: `iliasApp/views.py:136-143`.

Implementation confirms actual sampled mutagenesis. `mutateString` randomly chooses `noOfPos` distinct, non-`N` positions and then randomly chooses a replacement from `mapperDict` at each position (`part3.py:12-48`). For A/C/G/T, `mapperDict` contains only the three non-parent bases (`myfunctions.py:7-8`), so these are genuine substitutions rather than random draws that may silently reproduce the parent.

The paper also states:

> “the Transmutation tool allows for the design of negative controls by permitting scrambling, reversing, complementing or introducing multiple random mutations in the input sequences or motifs.”

Source: `Georgakopoulos-Soares2017_MPRAnator.pdf`, p. 1. The mutation path qualifies; the separate scramble/reverse/complement paths do not.

**Disconfirmation attempt / what would change the value.** I searched all repository Python, templates, README, downloadable clients, tests, and the paper for `random`, `sample`, `seed`, `rng`, `rate`, `combination`, `SNP`, `Transmutation`, and `mutate`. I separately traced `Part3Form -> part3RresultsView -> mutateString -> mapperDict` to rule out the possibility that “random mutation” meant only scrambling or random sequence generation. If the handler only shuffled controls, if `numToMutate` did not influence edit selection, or if this existed only as an example script, the value would be `no` or at most `partial`. It is a documented, full web-tool module with a direct count parameter, so `yes` is consistent with PoolParty's one-draw random mode.

### mpradesigntools / designMPRA — **partial** (medium confidence)

**Qualifying but restricted evidence.** The exported public entry point is `processVCF` (`NAMESPACE:3-4`). Its manual page documents:

> “`alter_aberrant`: under development - logical indicating whether to randomly alter aberrant digestion sites across barcodes”

Source: `man/processVCF.Rd:55-56` at `andrewGhazi/mpradesigntools@afd386e`.

The same call exposes `nper`, “the number of barcoded sequences to be generated per allele per SNP” (`man/processVCF.Rd:7-31`). In implementation, `nper` creates `2*nper` reference/alternate rows (`R/processVCFfast.R:389-399`). When `alter_aberrant` is enabled and a context digestion site is present, `randomly_fix` builds every possible one-base disruption of the site, samples `nrow(res_df)/2` of those choices (with replacement only when necessary), assigns the same sampled repair set to both alleles, and emits the repaired constructs (`R/processVCFfast.R:112-185`, particularly 150-172; call sites at 426-481). Its own generated manual says:

> “For a SNP with aberrant digestion sites in the context, randomly change bases in the site across barcodes”

Source: `man/randomly_fix.Rd:22-24`.

This is real sampled mutation of parent-derived construct context, but only as a restriction-site repair attached to barcode replication. It is not a general random-mutagenesis specification; it is off by default, accepts no mutation rate or seed, fails when more than one aberrant site is present, and the public option is explicitly “under development.” Those stated restrictions warrant `partial`, not `yes`.

**Disconfirmation attempt / what would change the value.** I searched every R source and `.Rd` page in mpradesigntools, its README/DESCRIPTION/NAMESPACE, the designMPRA Shiny repository and scripts, and all pages of `Ghazi2018_MPRADesignTools_all.pdf` for `random`, `sample`, `seed`, `rng`, `mutate`, `variant`, `alter_aberrant`, and `nper`. I separated random barcode selection and statistical resampling in analysis scripts from the qualifying `randomly_fix` context edits. A documented general API that sampled substitutions/indels from parent sequences, or a seed/rate/sample-count parameter independent of digestion-site repair, would have produced `yes`; removal of the repair code from the exported `processVCF` path would have produced `no`. Neither condition holds.

### Oligopool Calculator — **no** (high confidence)

**Evidence of scope.** The current first-party README defines the design surface as:

> “Design modules generate primers, barcodes, motifs/anchors, and spacers; assembly modules split/pad long constructs; Degenerate Mode compresses similar sequences into IUPAC-degenerate oligos...”

Source: `README.md:20-22` at `ayaanhossain/oligopool@b88fa39`; the exhaustive feature/module list is at `README.md:34-48` and `README.md:125-157`. There is no mutagenesis generator in that list.

The tool does expose `random_seed` on barcode, primer, motif, spacer, split, pad, and compression operations. These are excluded near matches: barcode/spacer/motif/primer generation creates new design elements, while compression's Monte Carlo rollouts choose a compression grouping but do not create sampled parent variants. The compression docs explicitly require already concrete sequences as input and guarantee “no invented sequences” (`docs/api.md:778-818`).

Most decisively, the authors' agent-ready workflow table assigns variant generation to the caller:

> “Cost-efficient saturation mutagenesis | generate all substitutions -> `compress`...”

and repeats “generate all substitutions” for both barcoded and degenerate saturation workflows. Source: `docs/agent-skills.md:371-378`.

**Disconfirmation attempt / what would change the value.** I searched all top-level and `oligopool/base` Python modules, `README.md`, `docs/api.md`, `docs/docs.md`, `docs/agent-skills.md`, examples (including the library-compressor demo), and every paper page for `random`, `sample`, `seed`, `rng`, `mutation`, `mutagenesis`, `variant`, and `shuffle`; I also enumerated every public `def`/`class`. A public function taking parent sequences plus a mutation rate/event count/output count and producing sampled derived variants would have changed the value to `yes` or `partial`. The only random parameters belong to excluded element design, assembly heuristics, compression heuristics, or failed-read sampling; documented mutagenesis workflows require externally generated variants. Hence `no`.

### Mutation Maker — **no** (high confidence)

**Evidence and excluded near matches.** The authors describe SSSM input as a parental gene plus “a list of mutation sites with their respective degenerate codon (default, `NNK`)” and describe the brute-force algorithm as precomputing “all potential primers around all mutation sites.” Source: `Hiraga2021_MutationMaker.pdf`, p. 3. The tool designs the mutagenic primers/degenerate oligos; it does not draw concrete library members from the encoded wet-lab distribution.

For multisite work, source constructs all concrete codon combinations with `itertools.product` (`backend/mutation_maker/mutation.py:210-238`). For PAS, the documented mutation `frequency` is carried into `PASMutation`, oligo metadata, and mixture ratios (`pas_types.py:105-170`; `pas_output.py:205-242`). The paper is explicit that the algorithm “computes the oligo ratios in the reaction mixture that are necessary to achieve the user-defined mutation frequency in the final library” (`Hiraga2021_MutationMaker.pdf`, p. 7). That is a physical mixing prescription, not a software sample of concrete variants.

Two source-level random operations were checked and excluded:

- `Translator.generate_dna` “randomly generates DNA from protein sequence based on codon usage frequency” (`reverse_translation.py:53-65`). This is de novo reverse translation/random sequence generation, expressly outside the row.
- `QCLMSolver.pick_random_codon` chooses one permissible synonymous codon during primer/solution construction (`qclm.py:844-851`). That is an internal heuristic choice, not an output sample of a declared variant space.

**Disconfirmation attempt / what would change the value.** I searched every backend Python module, API server, frontend TypeScript/TSX input model, README, tests/features, and all paper pages for `random`, `sample`, `seed`, `rng`, `rate`, `frequency`, `mutation`, and `mutagenesis`, and traced each RNG call (`np.random.choice`, `random.choice`). A documented API that expanded a degenerate design and returned N sampled concrete parent variants, or a mutation-rate/count parameter that selected concrete edited sequences, would have produced `yes`/`partial`. No such API or seed/sample field exists. Claims of “random gene libraries” refer to wet-lab degeneracy; reverse translation and heuristic codon choice are excluded by the instrument. Therefore `no`.

### DNA Chisel — **yes** (medium confidence)

**Qualifying evidence.** DNA Chisel documents `MutationSpace` as a core class (`docs/ref/core_classes.rst:45-49`) and explains that it represents possible mutations at different locations and that its methods “describe how to extract variants from the mutation space” (`dnachisel/README.md:26`). Its documented methods include:

```python
def pick_random_mutations(self, n_mutations, sequence):
    """Draw N random mutations."""
```

and:

```python
def apply_random_mutations(self, n_mutations, sequence):
    """Return a sequence with n random mutations applied."""
```

Source: `dnachisel/MutationSpace/MutationSpace.py:106-130` at `Edinburgh-Genome-Foundry/DnaChisel@68c0930`. The method samples mutation locations without replacement and calls `MutationChoice.random_variant`, which explicitly excludes the current subsequence before drawing an alternative (`MutationChoice.py:38-46`). Thus the count parameter produces a genuine parent-derived random mutant, not a random de novo sequence.

RNG reproducibility is also documented as user-controlled:

> “The randomness in DNA Chisel is entirely determined by the Numpy random generator. As a consequence, setting the Numpy seed at the very beginning of a script ensures that the result will always be the same.”

Source: `docs/notes.rst:4-20`.

**Boundary with excluded optimization search.** `DnaOptimizationProblem.optimize_by_random_mutations` does use this sampler internally (`DnaOptimizationProblem/mixins/ObjectivesMaximizerMixin.py:62-100`), and candidates visited only by that optimization loop do not count. The score instead rests on the separately documented core method that directly returns the sampled sequence. Users can obtain the problem's constructed `mutation_space` and call the method without invoking optimization. The restriction is that each call returns one sampled sequence rather than a whole N-member library; the row definition asks for “a sample” and accepts a count parameter, so the same one-draw threshold used for PoolParty and MPRAnator yields `yes`, not `partial`.

**Disconfirmation attempt / what would change the value.** I searched the full package, documentation reference, README/PyPI README, examples, changes, and paper for `random`, `sample`, `seed`, `rng`, `mutation_rate`, `variant`, `pick_random_mutations`, and `apply_random_mutations`. I distinguished `random_dna_sequence`, `random_compatible_dna_sequence`, randomized codon back-translation, and solver-internal stochastic search as excluded. I attempted a direct runtime check from the shallow clone; it could not start because Biopython is not installed in the environment, and the task forbids installing dependencies, so no behavioral claim relies on that attempt. If `apply_random_mutations` were private/undocumented or only returned an optimizer's winning sequence, the value would be `no`; if it were documented only in an example, it would be at most `partial`. It is a member of the documented core `MutationSpace` API and directly returns the mutant, so `yes`, with medium rather than high confidence because its primary documented role is solver support rather than library export.

### tangermeme — **no** (medium-high confidence)

**Evidence and excluded near matches.** The closest API is `ersatz.randomize`:

> “Replace a region of the provided loci with randomly drawn sequence.”

It accepts `n` and `random_state` and returns `(-1, n, alphabet, length)` tensors (`tangermeme/ersatz.py:343-397` at `jmschrei/tangermeme@2006b31`). However, its implementation independently generates an entirely new one-hot sequence for every position in the selected interval and substitutes that block (`ersatz.py:414-422`). The first-party tutorial frames it exactly as the excluded case:

> “An alternative approach to deleting positions is to replace those positions with randomly generated characters. This would keep the sequence the same length but remove the motif.”

Source: `docs/tutorials/Tutorial_A1_Ersatz_Sequence_Manipulation.ipynb`, markdown cell 23. Cell 25 says `n` creates many “randomizations” to average over random sequence backgrounds. This is random block/background generation, not selection of mutation events relative to the parent, and is therefore barred by “Random sequence generation is not sampled mutagenesis.”

The explicit `shuffle` and `dinucleotide_shuffle` APIs are also barred shuffled controls (`ersatz.py:425-498`, `608-688`). `design.screen` “randomly generate[s] a batch of examples” and retains the model's best hits (`design/screen.py:25-44`), which combines de novo generation with model-guided screening rather than sampled mutagenesis. In contrast, the actual saturation-mutagenesis function creates all true single-base edits at every position (`saturation_mutagenesis.py:212-240`); it has no sample-count/rate/RNG option.

**Disconfirmation attempt / what would change the value.** I searched all 28 package Python modules including `design/` and `_skills/`, all RST/README documentation, every tutorial notebook source, tests, and the paper for `random`, `randomize`, `sample`, `seed`, `rng`, `mutation`, `variant`, `shuffle`, `screen`, and `saturation`. A function that chose a requested number of mutation positions/events from a parent, sampled alternatives different from the parent, or subsampled the saturation variant set would have yielded `yes`/`partial`. The repository instead contains exhaustive substitution, random block generation/backgrounds, shuffling, and model-based screening. Under the binding exclusion, none qualifies, so `no`.

## Row-level finding

The row works on one consistent scale, but only if “sampled mutagenesis” is kept distinct from generic stochastic sequence manipulation. That distinction is empirically useful here:

- PoolParty, MPRAnator, and DNA Chisel explicitly choose mutation events/alternatives relative to a parent and return the resulting draw.
- mpradesigntools does so only as a narrow, under-development restriction-site repair across barcode replicates.
- VaLiAnT enumerates its supported mutation spaces.
- Oligopool Calculator expects variants to be supplied before its design/compression stages.
- Mutation Maker encodes diversity in degenerate oligos and physical mixture ratios but does not sample concrete members in software.
- tangermeme creates random replacement backgrounds or shuffled controls, while its actual mutagenesis routine is exhaustive.

Thus the exclusion sentence in the measurement instrument is load-bearing, but it can be applied consistently without replacing or redefining the row.
