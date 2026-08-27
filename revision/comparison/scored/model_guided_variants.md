# `model_guided_variants` audit

Audit date: 2026-08-15. This is an independent application of
`revision/tables/ROW_DEFINITIONS.md`; no prior score was used. The files in
`tool_survey/final/` were consulted only as locator leads and are not evidence
for any value below.

## Binding operational test

**I score `yes` only when a tool-provided operation runs an optimization loop in
which predictive-model outputs determine sequence edits, `partial` when it
merely accepts an arbitrary scoring callable without optimizing it (or the
capability exists only in an example), and `no` only after primary-source
searches show neither interface.**

This applies the global rule: a user-written reconstruction from unrelated
primitives is `no`, while implementing a model behind an explicitly documented
scoring/objective interface is not reconstruction. It also applies row 5's
boundary: a multi-edit sequence reached by an objective-driven optimizer takes
full credit here and at most `partial` in `pairwise_higher_order`; it is not
double-counted here as declarable combinatorial enumeration.

## Results

| tool | value | confidence |
|---|---|---|
| poolparty | **partial** | high |
| valiant | **no** | high |
| mpranator | **no** | high |
| mpradesign | **no** | high |
| oligopoolcalc | **no** | high |
| mutationmaker | **no** | high |
| dnachisel | **yes** | high |
| tangermeme | **yes** | high |

### poolparty — **partial** (high confidence)

The documented operation accepts exactly the model-shaped hook required for
`partial`, but it is deliberately a passthrough rather than an optimizer:

> “Evaluate a user-supplied function on each sequence and record the result as a
> design card column. The sequence passes through unchanged.”

Source: local canonical repository commit
`1bb0179e1c3720b1fffd471802b3040f9336de28`,
`poolparty/docs/operations/score.rst:4-8`.

The API and implementation make both halves concrete:

> `fn: Callable[[str], Any]`

> “Scoring function. Receives a clean (tag-free) sequence string and returns a
> scalar value to record.”

> `result = self._fn(clean_seq)`
>
> `return parent_seq, {self._card_key: result}`

Source: same snapshot, `poolparty/src/poolparty/fixed_ops/score.py:12-26,35-37,107-117`.

The authors explicitly disconfirm an optimization loop:

> “PoolParty is not a sequence optimization tool, e.g., for optimizing codon
> usage, satisfying synthesis constraints on individual sequences, or designing
> sequences guided by machine learning models.”

Source: author manuscript `25_poolparty_private/26.05.18_bmc_submission/latex/main.pdf`,
PDF p. 13 (manuscript lines 278-283); matching source text in
`25_poolparty_private/26.05.18_bmc_submission/latex/main.tex:235`.

Behavioral check, using the mandated read-only venv and
`PYTHONDONTWRITEBYTECODE=1`: a mock model returned `0.1` for `AAAA` and `0.9` for
`CCCC`; `pp.score(...).generate_library()` returned the same two sequences with
those scores attached. No candidate was selected and no base changed.

**Disconfirmation attempt.** I searched all 102 package `.py` files and all
operation docs for `model|predict|optim|objective|score|loss|argmax|argmin|anneal|
gradient|callable`, then inspected every hit relevant to sequence scoring. The
only generic per-sequence scoring interface is `ScoreOp`, whose source returns
the parent unchanged. A `yes` would have required a tool-provided loop that
proposes edits, calls `fn`, and retains an edit based on its output; that exact
counterevidence was searched for and is absent. A `no` would have required no
arbitrary scoring hook, contradicted by `score(fn=...)`. Hence `partial`.

### valiant — **no** (high confidence)

VaLiAnT documents a fixed library-design scope and fixed mutation modes:

> “The Variant Library Annotation Tool (VaLiAnT) is an oligonucleotide library
> design and annotation tool for Saturation Genome Editing and other Deep
> Mutational Scanning experiments.”

Source: official repository commit
`8796cc112dafd4919fec59913f58cd2be87c45eb`, `README.md:1-5`. The documented
mutation-type list is parametric deletion, SNV, in-frame deletion, alanine,
stop, all-amino-acid, SNVRE, and custom variants (`README.md:20-28`).

The paper's only predictive-model passage puts models downstream of the
generated experimental data, not inside the design loop:

> “as predicative models and machine learning increase in utility ... training
> datasets rich in biological context will be important.”

Source: `Barbon2022_VaLiAnT_all.pdf`, PDF p. 7, Discussion. The next sentence
connects VaLiAnT to producing such reference-aware variant libraries, not to
loading or optimizing against a model.

**Disconfirmation attempt.** I searched 85 Python source files, the complete
README, all 11 official wiki Markdown pages, and the full paper. The exact
case-insensitive primary search was
`predictive|prediction|machine[ -]?learning|deep[ -]?learning|sequence[- ]to[- ]function|model[-_ ]guided|model[-_ ]in[-_ ]the[-_ ]loop|tensorflow|torch|keras|sklearn|pytorch|model[-_ ]loading|scor(e|ing)[-_ ]?(callable|function)`:
zero hits occurred in the 97 repository/wiki text-source files searched. A
broader `callback|callable|objective|loss|score|model` search produced only
Pydantic data models, GTF's `SCORE` column, Click's logger callback, and internal
typing callables for grouping/filtering/oligo retrieval; none is a public
sequence scorer. I also inspected the documented CLI parameters and source tree.
A `partial` would have required a public arbitrary sequence-scoring callable; a
`yes` would additionally have required a candidate-edit optimization loop.
Both were explicitly searched for and absent, so `no` is an observed absence.

### mpranator — **no** (high confidence)

The paper exhaustively names the four provided investigations:

> “The MPRA Motif design tool can be used to systematically generate synthetic
> sequences with single motifs or combinations of motifs placed at pre-selected
> positions. The MPRA SNP design tool can be used to examine the regulatory
> effects of single or combinations of SNPs ... The PWM Seq-Gen tool performs
> probabilistic realizations of PWMs or generates all the corresponding k-mer
> motifs exceeding a probability threshold. The Transmutation tool allows for
> the design of different types of negative controls.”

Source: `Georgakopoulos-Soares2017_MPRAnator.pdf`, PDF p. 2, Materials and
methods.

The source routes likewise expose MPRA motif/SNP, Transmutation, documentation,
and one SNP API (`MPRAnator` commit
`9969790d62410138d4281b7955da6d085f07b1bc`, `iliasApp/urls.py:7-31`). The
implementation performs user-declared combinatorics and random mutation, e.g.
`generatePermutations`, `generateCombinations`, and `mutateString` select
positions/bases with Python `random` (`part1.py:6-21`; `part3.py:13-48`).

The PWM probability threshold is not credited as model-guided editing: it
generates/decompresses motifs from the user-supplied PWM and neither predicts a
sequence property nor feeds a prediction back to choose an edit on a parent.
It therefore fails the same operational test applied to every tool.

**Disconfirmation attempt.** I searched every one of the repository's 36
Python/HTML/Markdown primary files plus both PDF pages. The predictive/model
interface regex used for VaLiAnT returned zero source/doc hits. The broader
`callback|callable|objective|loss|score|model` grep returned only Django database
`models` and WSGI's framework “callable”; it revealed no sequence score hook.
I inspected every function in `part1.py`, `part3.py`, `oligo.py`,
`myfunctions.py`, `viewsCore.py`, and `iliasApp/views.py`, plus every URL route.
A different value required a documented arbitrary predictor/scorer parameter
(`partial`) or an edit/search loop calling it (`yes`); both forms were searched
for and are absent. `no` is therefore not a failure-to-find hedge.

### mpradesign — **no** (high confidence)

The paper describes two functions, neither a sequence-to-function model loop:

> “an online web-tool and R package that allows for interactive MPRA
> experimental design encompassing both power analysis and design of
> constructs.”

> “Users can adjust experimental parameters to examine the predicted effect on
> assay power as well as upload VCFs for automated construct sequence
> generation.”

Source: `Ghazi2018_MPRADesignTools_all.pdf`, PDF p. 1, Abstract. PDF p. 2 says
that after parameter selection, users upload variants and receive construct
sequences. Thus the prediction is statistical power from assay parameters; it
does not determine a DNA edit.

The package documentation says its main design function is `processVCF` and
describes its fixed inputs (`mpradesigntools` commit
`afd386ef12051bb0a37ad63a6f92acd555246584`, `README.md:39-47,106-119`). The
namespace exports only `processVCF` and `spread_and_fix_indels`
(`NAMESPACE:1-4`).

**Disconfirmation attempt.** I searched all 53 R/Rmd/Rd/Markdown files across
both official repositories (`mpradesigntools` above and `designMPRA` commit
`0cf56ee602fc86dde705906d071a39cbdf6e99a8`) plus the paper and supplement. The
predictive/model-interface regex returned zero hits. A broader
`callback|callable|objective|loss|score|model` grep found only the statistical
power model in `Rmd/Supplement.Rmd:376`, a DESeq2 analysis comment, and the
unrelated word “INFO” in VCF documentation. I inspected `NAMESPACE`, all public
R function definitions, `server.R`, and `ui.R`. A `partial` required a callable
that scores DNA sequences; a `yes` required that score to guide edits. I
searched for both and found neither. Predicted assay power does not satisfy the
concrete edit test, so the score is `no`.

### oligopoolcalc — **no** (high confidence)

The official README separates constraint-based sequence design from downstream
analysis:

> “In `Design Mode`, `Oligopool Calculator` generates optimized `barcode`s,
> `primer`s, `spacer`s, and `motif`s.”

> “Once the library is assembled and cloned ... `Analysis Mode` proceeds by
> first `index`ing ... `pack`ing ... and then producing count matrices.”

Source: official repository commit
`b88fa394ca67ed4c48ec41127b5707694ee7cf0a`, `README.md:34-48`.

The paper places sequence-to-function modeling after measurement:

> “The resulting count matrices are analyzed to measure genetic system activity
> across conditions of interest and can be combined with biophysics and/or
> machine learning to develop predictive sequence-to-function models.”

Source: `Hossain2024_OligopoolCalculator.pdf`, PDF p. 2.

There are two tempting but excluded matches. First, the public `Scry` model is a
“1-nearest-neighbor barcode classifier” that “Powers `acount`/`xcount`
internally” (`docs/api.md:1567-1595`); the paper likewise says it identifies
barcodes in sequencing reads (PDF p. 7). It predicts barcode labels during
analysis, not edits during design. Second, the only public arbitrary callbacks
are `acount`/`xcount` “Custom read filter function[s]”
(`docs/api.md:1068-1099,1200-1219`; `oligopool/acount.py:15-60`). They filter
sequencing reads and cannot score design candidates.

**Disconfirmation attempt.** I inspected every documented API signature in
`docs/api.md`, all top-level API modules, all design engines, the complete
README/docs, and the paper (46 Python/Markdown files in the focused source/doc
scope). The predictive-interface regex yielded seven hits, all `Scry` barcode
classification. The broader hook search found only the two analysis callbacks
above, progress-printer `liner` callables, and internal fixed constraint
objectives (`barcode_objectives`, `motif_objectives`, `primer_objectives`). No
design command accepts `model`, `predictor`, `score_fn`, `objective_fn`, or a
design callback. A `partial` required such a public sequence scorer and a `yes`
required feedback from it into edits; both were explicitly sought. The shipped
Scry classifier and “AI-assisted” agent documentation (`README.md:49-59`) do not
change the result. Score: `no`.

### mutationmaker — **no** (high confidence)

Mutation Maker does optimize, but its score functions are fixed primer,
fragment, codon-usage, and GC calculations rather than outputs of an attachable
predictive model. The paper is explicit:

> “All workflows utilize score functions that are based on deviations of
> primers properties from the optimal values.”

Source: `Hiraga2021_MutationMaker.pdf`, PDF p. 9, Methods: Score Functions. It
later states, “The weights used in the score functions are hardcoded” (PDF p.
10).

The primary source matches that description. `PrimerScoring` is a tool-defined
function object whose return is a root-sum-square of melting-temperature, GC,
length, end-size, hairpin, and dimer deviations (official repository commit
`396c1c0ede7529f3dedf65e821e8c1f20c9a7043`,
`backend/mutation_maker/primer_scoring.py:34-94`). `TranslationScoring` computes
CAI and GC (`translation_scoring.py:32-44,57-91`). Fragment scoring is likewise
hard-coded from length and self-binding terms (`pas_solution.py:385-422`).

**Disconfirmation attempt.** I searched all 32 core backend package `.py` files,
the API/task entry points, frontend source, repo README/docs, and the full paper
(152 source/doc files in the wide text scope). The predictive/model-interface
regex returned only two textual “score function” hits and no predictive-model,
ML-framework, or model-guided hit. The broader hook grep was inspected in full:
it found Pydantic request models and the fixed `PrimerScoring`,
`TranslationScoring`, PAS/QCLM/SSM calculations above, but no public callback or
replaceable scorer. A `partial` required a user-supplied arbitrary sequence
scorer; a `yes` required a loop using its result to choose sequence edits. Both
were explicitly searched for and absent. Fixed physical/codon objective
optimization belongs to other rows, not this one. Score: `no`.

### dnachisel — **yes** (high confidence)

DNA Chisel provides a documented generic scoring interface and a genuine edit
optimization loop. The base class says:

> “New types of specifications are defined by subclassing `Specification` and
> providing a custom `evaluate` and `localized` methods.”

> `evaluate` — “function (sequence) => SpecEvaluation”

Source: official repository commit
`68c09304341c3656f3dfe63eda37757d6a7b3917`,
`dnachisel/Specification/Specification.py:14-31`. Its constructor installs a
provided evaluator (`Specification.py:68-72`). Custom specifications are a
documented feature, not a repo-only example: `docs/examples.rst:47-50` includes
“Writing a custom Specification: 9-mer score,” and `README.rst:27-31` says users
can define their own specifications and an objective's score is maximized.

The optimizer feeds scores back into edits. Exhaustive search iterates mutation-
space variants, evaluates `objective_scores_sum()`, and retains the sequence
when its score improves:

> `for variant ...: self.sequence = variant`
>
> `score = self.objective_scores_sum()`
>
> `if score > current_best_score: ... current_best_sequence = self.sequence`

Source: `dnachisel/DnaOptimizationProblem/mixins/ObjectivesMaximizerMixin.py:43-60`.
The random search analog applies random mutations and retains only score
improvements (`:70-103`); `optimize_objective` chooses exhaustive versus random
search (`:105-187`), and `optimize()` invokes that loop for each objective
(`:189-202`). A predictive model placed behind the documented `evaluate`
interface therefore directly determines which edit survives.

The paper independently states:

> “DNA Chisel allows to define any new sequence specification that the Python
> language can express”

and describes “objectives maximization with respect to the constraints.”

Source: `Zulkower2020_dnachisel.pdf`, PDF pp. 2 and 1-2 respectively.

**Disconfirmation attempt.** I searched all 134 Python/RST/Markdown files in
the package, docs, and examples for model/predict/custom-specification/objective/
score/optimization hooks and inspected the complete base specification and
objective maximizer. A downgrade to `partial` would require the evaluator's
score merely to be recorded, but the cited exhaustive and stochastic loops use
it to retain/reject edited candidates. A `no` would require user reconstruction,
but `Specification(evaluate=...)`, `objectives=[...]`, and `problem.optimize()`
are purpose-built documented interfaces. No shipped genomic model is required
by the row. I attempted a no-install runtime check with a four-base mock
predictive objective; it could not start because Biopython is not installed
(`ModuleNotFoundError: Bio`). Per the no-install rule I did not install it;
source and paper settle the behavior directly. Score: `yes`.

### tangermeme — **yes** (high confidence)

Tangermeme exposes model-guided edit search directly. The documented
`greedy_substitution` operation says:

> “This design function will greedily add motifs to achieve a desired output
> from the model. Each round, the function will iterate through all possible
> motifs, substitute each one ... and keep the one whose loss function is the
> smallest.”

Source: official repository commit
`2006b310cd72a28c56c3ea4ba67f738fff74bb89`,
`tangermeme/design/greedy_substitution.py:24-53`. Its parameters are an arbitrary
PyTorch `model`, target `y`, motifs, and a loss callable (`:24-39,58-88`).

The loop materializes every allowed candidate placement, predicts each one,
computes model-output loss, takes `argmin`, and applies the winning edit:

> `y_hat = predict(model, X_, ...)`
>
> `loss_curr = loss(..., y_hat[...])`
>
> `pos = loss_curr.argmin()`
>
> `X = substitute(X, motifs[best_motif_idx], start=best_pos, ...)`

Source: `greedy_substitution.py:194-247`. The companion beam search explicitly
searches “trajectories through edit-space,” scores every complete candidate by
loss, and can recover multi-edit combinations (`design/beam_substitution.py:44-72`).
The API is documented through `docs/api/design.rst:7-10`, and the executed
design tutorial calls `greedy_substitution(model, X, y, motifs, ...,
max_iter=5)` (`docs/tutorials/Tutorial_B6_Design.ipynb`, source cells containing
that call and the explanation that `y` is the target predictive-model output).

The paper says tangermeme focuses on “using trained models to drive genomics
research” and Fig. 1 lists Design / Greedy Substitution (PDF pp. 2-3), consistent
with the source.

**Disconfirmation attempt.** I searched all 66 Python/RST/Markdown files in the
package and docs and inspected all five `tangermeme/design/` algorithms plus the
design tutorial. A downgrade to `partial` would require model output merely to
be reported; the cited loop instead changes `X` based on prediction loss. A
`no` would require user-authored bookkeeping around prediction and editing; the
single documented operation performs enumeration, prediction, comparison, and
edit application itself. The default `greedy_substitution(max_iter=-1)` is a
documented current no-op, but the supported positive `max_iter` parameter runs
the loop and the tutorial supplies it; this is a parameter caveat, not absence
or example-only capability. Score: `yes`.

## Row-level finding

The row discriminates on one consistent scale: five tools are `no`, PoolParty is
the exact documented `partial` case (attachable scoring with no feedback), and
DNA Chisel/tangermeme are `yes` because scores/predictions select edits inside a
tool-provided optimizer. The important boundary is not whether a tool contains
something called a “model” or “optimization”: Oligopool Calculator's Scry model
is downstream barcode decoding, Mutation Maker's objectives are fixed physical
scores, MPRA Design Tools predicts assay power, and MPRAnator's PWM generator
realizes an input motif distribution. None makes a predictive output choose a
sequence edit under the operational test.

DNA Chisel demonstrates a consequence of the binding definition: a tool need
not ship a trained genomic model to score `yes`; a documented arbitrary
objective evaluator plus the tool's own score-driven edit loop is sufficient.
This is applied symmetrically to all eight tools.

## Primary-source snapshots and complete search log

Repository snapshots used (all cloned read-only into `/tmp`, except the mandated
local PoolParty repository): PoolParty `1bb0179...`; VaLiAnT `8796cc1...`
(`develop`) plus its official wiki; MPRAnator `9969790...`; mpradesigntools
`afd386e...`; designMPRA `0cf56ee...`; oligopool `b88fa39...`; Mutation Maker
`396c1c0...`; DNA Chisel `68c0930...`; tangermeme `2006b31...`. The eight
author papers listed in the task were read with PyMuPDF (`fitz`) page text.

Searches/actions performed, in order:

1. Read `revision/tables/ROW_DEFINITIONS.md:1-170`, including the global rule,
   row 5 boundary rules, and row 12 definition; stated the operational test
   before opening any tool record or source.
2. Listed `tool_survey/` and local repository directories with `find`; searched
   the eight `tool_survey/final/*.md` lead files for
   `model|predict|optim|score|objective|github|repository|source|callable` only
   to locate likely primary files. No lead text is cited as evidence.
3. Extracted every supplied PDF in memory with `fitz`; first searched
   `model|predict|optim|objective|score|machine learning|deep learning|callable`,
   then printed page-numbered context for tool-specific hits (including
   `predicative models`, `PWM`, `combined with biophysics`, `score functions`,
   `optimization objective`, `trained models to drive`, and PoolParty's explicit
   limitation). Printed the first two full pages for MPRAnator and MPRA Design
   Tools to verify their complete stated scopes.
4. Searched official GitHub results on the web for tangermeme design,
   DNA Chisel custom objectives, VaLiAnT model/scoring, and MPRAnator
   model/scoring. Only official repository pages were opened; the cloned
   commit-pinned sources above provide the quoted evidence.
5. Cloned the seven external canonical repositories, the second MPRA Design
   repository, and the VaLiAnT official wiki into `/tmp`; recorded each commit
   with `git rev-parse HEAD`.
6. For every clone, listed the complete primary-source tree with `rg --files` or
   `find`, counted relevant source/doc files, and ran a broad first-pass grep for
   `predictive|prediction|machine learning|deep learning|model guided|
   sequence-to-function|objective|optimize|score|callable|anneal|loss` while
   excluding binary/generated assets.
7. Ran the same focused regex across the five proposed `no` tools:
   `predictive|prediction|machine[ -]?learning|deep[ -]?learning|
   sequence[- ]to[- ]function|model[-_ ]guided|model[-_ ]in[-_ ]the[-_ ]loop|
   tensorflow|torch|keras|sklearn|pytorch|model[-_ ]loading|
   scor(e|ing)[-_ ]?(callable|function)`. Results: VaLiAnT 0/97 text files;
   MPRAnator 0/36; MPRA Design 0/53; Oligopool Calculator 7/48, all
   Scry prediction; Mutation Maker 2/152, both fixed “score function” text.
8. Ran `callback|callable|objective|loss|score|model` separately across each
   proposed `no` source/doc scope and manually classified every relevant hit as
   described in its section. Inspected public signatures/exports and route
   lists (`NAMESPACE`, top-level `def` lists, `urls.py`, CLI/config docs) so that
   a scoring hook with different terminology would not be missed.
9. PoolParty: searched all 102 source modules and docs for model/scoring and all
   optimization-loop synonyms; opened `fixed_ops/score.py`, the score operation
   docs, manuscript source, and ran the mock-model passthrough test with the
   specified read-only Python 3.12 venv and bytecode disabled.
10. DNA Chisel: searched all 134 package/doc/example text files; opened the full
    `Specification` evaluator interface, custom-spec example/docs, and all of
    `ObjectivesMaximizerMixin.py`. Attempted the dependency-free mock runtime
    check; stopped at missing Biopython and installed nothing.
11. Tangermeme: searched all 66 package/doc text files; opened all design module
    files, the API page, design tutorial source cells, README, paper, and the
    complete greedy/beam optimizer implementations.

No value is `unknown`: the primary sources were sufficient to determine every
cell.
