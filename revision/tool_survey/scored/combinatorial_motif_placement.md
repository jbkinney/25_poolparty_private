# Audit: combinatorial multi-motif placement

## Measurement instrument and operational test

I read `revision/tables/ROW_DEFINITIONS.md` in full through the relevant
definitions before inspecting any tool. I applied its global rule: a capability
counts only when the tool supplies an operation, parameter, or mode for it; a
caller-written loop, manual list merge, or other reconstruction is `no`, not
`partial`. I also applied row 5's binding boundary: co-placement of motifs is
credited here (row 9), not as pairwise/higher-order mutagenesis. Conversely,
combinations of substitutions, insertions, or deletions relative to a parent are
mutation events and do not establish this row. IUPAC expansion alone is likewise
not motif-placement design.

**Operational test (fixed before tool inspection):** for each tool, look for a
documented operation, parameter, or mode that takes a motif-placement
specification and itself enumerates or samples concrete sequences placing at
least two distinct motifs across position and/or orientation combinations;
score `partial` only for a tool-provided single-motif placement or position
enumeration with fixed orientation, and `no` when users must construct or merge
the combinations themselves.

The test treats model-selected enumeration as enumeration: the definition says
"enumerated or sampled," not necessarily exhaustively returned. It does not
credit transient internal candidates unless the documented operation can return
designed sequence members. This distinction matters for DnaChisel and
tangermeme below.

## Results

| Tool | Value | Confidence |
|---|---|---|
| poolparty | **yes** | high |
| valiant | **no** | high |
| mpranator | **yes** | high |
| mpradesign | **no** | high |
| oligopoolcalc | **no** | high |
| mutationmaker | **no** | high |
| dnachisel | **no** | high |
| tangermeme | **yes** | high |

## Source inventory and common search protocol

All evidence below is primary: authors' repositories, repository documentation,
and authors' papers. The `tool_survey/final/*.md` records were used only as
navigation leads and are not evidence for any decision.

Repository snapshots searched:

| Tool | Repository snapshot |
|---|---|
| poolparty | local `poolparty-statecounter`, commit `1bb0179e1c37`, origin `https://github.com/jbkinney/poolparty-statetracker.git` |
| VaLiAnT | `cancerit/VaLiAnT`, commit `8796cc112daf` |
| MPRAnator | `hemberg-lab/MPRAnator`, commit `9969790d6241` |
| mpradesigntools | `andrewGhazi/mpradesigntools`, commit `afd386ef1205` |
| designMPRA | `andrewGhazi/designMPRA`, commit `0cf56ee602fc` |
| Oligopool Calculator | `ayaanhossain/oligopool`, commit `b88fa394ca67` |
| Mutation Maker | `ra100/Mutation_Maker`, commit `396c1c0ede75` |
| DnaChisel | `Edinburgh-Genome-Foundry/DnaChisel`, commit `68c09304341c` |
| tangermeme | `jmschrei/tangermeme`, commit `8d732b8c0876` |

The non-poolparty clones were clean according to `git status --short`.
Poolparty had unrelated pre-existing changes to `poolparty/README.md` and
untracked built documentation, but `git diff --quiet` confirmed that every file
cited here (`insertion_multiscan.py`, `region_multiscan.py`, the regulatory
grammar tutorial, and the insertion-multiscan documentation) matched the cited
commit.

For every repository I searched, case-insensitively, for the concepts
`motif`, `binding.?site`, `position`, `spacing`, `orientation`,
`reverse.?complement`, `permut`, `product`, `combin`, `insert`, and the names of
any tool-specific nearby APIs. Searches covered source, README/documentation,
examples, and tests while excluding binary images and vendored/minified files.
I also extracted the supplied PDFs with PyMuPDF and searched each for `motif`,
`binding site`, `position`, `spacing`, `orientation`, `permutation`, and
`combinatorial`. Per-tool search logs below identify the closest hits and the
negative spaces used to support every `no`.

## poolparty — **yes** (high confidence)

### Evidence

Poolparty has a purpose-built public operation rather than a caller-assembled
loop. Its signature accepts the background pool, the number of simultaneous
insertions, one pool per distinct inserted motif, candidate positions,
permutation policy, spacing limits, and enumeration/sampling mode:

> `def insertion_multiscan(... num_insertions, insertion_pools, positions=...,
> insertion_mode=..., min_spacing=..., max_spacing=..., mode=...,
> num_states=...) -> Pool`

Source: `poolparty/src/poolparty/multiscan_ops/insertion_multiscan.py:12-30`.

The public docstring is explicit:

> "Insert or replace sequences at multiple positions simultaneously."

> "If a Sequence of Pools is provided, its length must equal
> `num_insertions`."

> "`'unordered'` uses all permutations."

> "Selection mode: 'random' or 'sequential'."

Source: the same file, lines 31-82, especially 32, 40-47, 58-72, and 79-82.
The rendered-source documentation independently says that a list assigns one
pool per site, `None` positions allow any valid non-overlapping arrangement,
and unordered mode allows any permutation
(`docs/operations/insertion_multiscan.rst:30-43,61-73`).

The authors' regulatory-grammar tutorial is a direct, documented multi-motif
use:

```python
cre_pool = template.insertion_multiscan(
    region="cre",
    insertion_pools=[hnf4a, ppara, xbp1],
    insertion_mode="unordered",
    replace=True,
    min_spacing=0,
    num_insertions=3,
    mode="random",
    num_states=1000,
    names=["hnf4a", "ppara", "xbp1"],
    cards={"starts": "positions", "names": "tfbs"},
).repeat(times=3)
```

Source: `docs/tutorials/mpra_regulatory_grammar.rst:90-113`. The tutorial states:

> "Because `flip` uses sequential mode, it exhaustively enumerates both
> orientations for each TFBS rather than sampling. With three TFBSs, this gives
> 2^3 = 8 orientation combinations per position configuration."

Source: lines 115-120. It later shows position, spatial-order, and strand columns
changing independently and says exactly that at lines 191-219.

### Behavioural check

Per the assignment, I ran code using the supplied read-only Python 3.12
environment with `PYTHONDONTWRITEBYTECODE=1`. The test used a ten-base
background, two distinct non-palindromic motifs (`AAC`, `CCG`), each motif's
documented sequential `flip`, two allowed insertion positions, unordered
assignment, and sequential enumeration. The operation returned eight rows and
eight unique sequences: two motif-to-position permutations times two
orientations for each motif. Representative rows were:

```text
AACAAAAACCGAAAAA [0, 5] [_ins_0, _ins_1] forward forward
CCGAAAAAAACAAAAA [0, 5] [_ins_1, _ins_0] forward forward
GTTAAAAACGGAAAAA [0, 5] [_ins_0, _ins_1] rc      rc
CGGAAAAAGTTAAAAA [0, 5] [_ins_1, _ins_0] rc      rc
rows 8 unique_seqs 8
```

This settles that the documented parameters create concrete library members,
not merely metadata or internal search states.

### Disconfirmation attempt and counterfactual

I searched the complete `poolparty/src`, `poolparty/docs`, and tests for
`insertion_multiscan|replacement_multiscan|motif|orientation|flip|spacing|positions|insertion_mode`,
then inspected the implementation's pool normalization and the underlying
`region_multiscan` enumeration (`region_multiscan.py:25-122,125-219`). Evidence
that the operation accepted only one motif, fixed the motif order and strand,
or returned only one solver-selected sequence would have reduced the score to
`partial` or `no`; the code, tutorial, and behavioural result showed the
opposite. I also checked for the row-5 boundary: although the mechanism is named
"insertion," the documented objects inserted are distinct TFBS motifs, so this
is placement here and is not double-counted as higher-order mutagenesis.

## valiant — **no** (high confidence)

### Evidence and closest functionality

VaLiAnT documents a closed set of mutation generators:

> parametric deletion, single-nucleotide variant, in-frame deletion, alanine
> codon substitution, stop codon substitution, all amino acid codon
> substitution, and SNVRE.

Source: `README.md:255-270`. The code enum contains exactly `DEL`, `SNV`,
`SNV_RE`, `IN_FRAME`, `ALA`, `STOP`, and `AA`
(`src/valiant/mutator_type.py:22-33`). None is motif placement.

The closest positive-looking capability is caller-authored insertion in a VCF:

> "Applied to the targeton reference sequence as a whole. Only simple variants
> such as the following are supported: substitutions, insertions, deletions,
> indels."

Source: `README.md:465-474`. The paper confirms that these are already specified
variants, not generated arrangements:

> "Simple variants are currently supported, including substitutions,
> insertions, deletions and indels."

> "Variants that start and end within the targeton are applied, generating one
> oligonucleotide sequence each."

Source: Barbon et al. 2022, PDF p. 5, section 2.5.2.

Thus a user could hand-author a VCF insertion whose ALT happens to be a TFBS,
but VaLiAnT does not enumerate motif × position, order, spacing, or orientation.
The global rule makes that `no`.

### Search and disconfirmation attempt

I ran `rg -i` for
`motif|binding.?site|spacing|permut|combin|orient` over `src/valiant`,
`README.md`, and any documentation. File-hit counts were respectively
`0, 0, 1, 0, 0, 0`; the one spacing hit is an input-parser statement that
whitespace is ignored, not biological spacing. I separately searched
`insertion|insert|reverse.?complement` and inspected `variant.py`,
`custom_variant.py`, `loaders/vcf.py`, `seq.py`, `sge_cli.py`, the full
`mutators/` directory, README mutation/custom-variant sections, and the supplied
paper (especially pp. 3-5). I also listed all files under `src/valiant` to check
for a differently named placement module.

Evidence that would have changed the value was a documented mutator or CLI
parameter accepting a motif set plus placement positions, strand/order, or
spacing and generating the combinations. I explicitly searched for each such
concept and found neither a single-motif sweep (`partial`) nor a multi-motif
one (`yes`).

## mpranator — **yes** (high confidence)

### Evidence

The authors' paper states the row almost verbatim:

> "The MPRA Motif design tool can be used to systematically generate synthetic
> sequences with single motifs or combinations of motifs placed at pre-selected
> positions."

Source: Georgakopoulos-Soares et al. 2017, PDF p. 2, Materials and methods.

The deployed repository's documented home-page description says MPRA Motifs:

> "Allows the substitution of motifs in a background sequence at all possible
> positions within the constraints of the parameters."

Source: `iliasApp/views.py:136-143`.

This is implemented as one tool-provided workflow. The view parses the motif
FASTA into `motifsL`, calls `generateCombinations(motifsL)`, then calls
`oligo.oligo` for every tool-generated motif subset
(`iliasApp/views.py:70-99`). The combination generator is:

```python
for i in xrange(1, len(theList) + 1):
    for acomb in itertools.combinations(theList, i):
        combinationsStore.append(acomb)
```

Source: `part1.py:15-21`.

Within one such subset, the placement operation enumerates the Cartesian
product of all candidate positions, one position axis per motif:

```python
positionsForMotifsL = [p for p in
    itertools.product(range(distanceFromLeftEdge, len(BackgroundS)),
                      repeat=len(motifsL))]
```

It writes each motif at its assigned position, filters by min/max spacing and
right-edge constraints, and appends a concrete sequence and header containing
motif-position pairs (`oligo.py:35-78`). Because each distinct motif has an
independent position axis, spatial order permutations arise without caller
bookkeeping.

### Disconfirmation attempt and counterfactual

I searched every tracked repository file for
`motif|permut|product|position|orientation|reverse|combin`, then inspected
`part1.py`, `oligo.py`, `viewsCore.py`, `iliasApp/views.py`, `forms.py`, the
templates/documentation, tests, and the paper. An older helper in `part1.py`
contains a problematic reverse-background branch, but it is not needed for this
score: multiple distinct motifs and their position/order combinations are
already enumerated by the live view and `oligo.py` with fixed motif orientation.
The row's full criterion is position **and/or** orientation combinations, so
orientation enumeration is not required when multi-position/multi-order
enumeration is present.

Evidence that positions were supplied one complete arrangement at a time, or
that callers had to loop over motif subsets and concatenate results, would have
made this `no` (or at best `partial` for a genuine one-motif tool mode). Instead,
the documented MPRA Motifs submission performs both subset and positional
enumeration itself.

## mpradesign — **no** (high confidence)

### Evidence and closest functionality

The package describes its scope as:

> "Currently the main function of MPRA Design Tools package is to design a set
> of barcoded sequences for MPRA experiments ... This is done with the
> `processVCF` function."

Source: `mpradesigntools/README.md:39-48`. Its only exported package operations
are `processVCF` and `spread_and_fix_indels`
(`NAMESPACE:1-4`). The API documentation says:

> "`processVCF` takes a VCF of SNPs ... and turns them into a set of labeled
> MPRA sequences barcoded with inert n-mers."

Source: `man/processVCF.Rd:78-87`. Its signature contains VCF, context lengths,
primers, restriction-enzyme patterns, barcode options, and output path, but no
motif set, candidate positions, strand sweep, spacing, order, or placement count
(`man/processVCF.Rd:5-76`).

The paper likewise says users "upload VCFs for automated construct sequence
generation" and describes constructs around reference/alternate alleles
(Ghazi et al. 2018 PDF, pp. 1-2). The Shiny UI says only CHROM, POS, REF, ALT and
the strand-detection INFO field are used (`designMPRA/ui.R:117-128`).

### Search and disconfirmation attempt

I searched both author repositories (`mpradesigntools` and `designMPRA`) for
`motif|binding.?site|spacing|permut|combin|orient|insert`, covering all R source,
man pages, README, Shiny UI/server, scripts, and the complete paper supplement.
In the package, file-hit counts for the first six concepts were
`0,0,0,0,0,1`; the orientation hit is only genomic-context strand handling. In
the app they were `0,0,1,0,3,2`; the spacing/combine hits concern prose or
statistical replicate analysis, not sequence placement. The paper itself had
zero occurrences of `motif`, `binding site`, `spacing`, `permutation`,
`orientation`, and `combinatorial` in extracted text.

I inspected the insertion hits in `R/processVCFfast.R` and
`scripts/processVCF*.R`; all are VCF indel handling. Evidence that would have
changed the score was a documented motif/element argument plus automatic
position/orientation/spacing enumeration. I searched the exported API, both app
and package implementations, and paper for it; no single-motif placement mode
(`partial`) or multi-motif mode (`yes`) exists.

## oligopoolcalc — **no** (high confidence)

### Evidence and closest functionality

Oligopool Calculator has real motif design, so this cell cannot be decided from
the word "motif" alone. Its documented operation is one motif/anchor output
column at a caller-chosen architectural location:

> "Purpose: Insert sequence motifs (per-variant or constant anchors) with
> constraint satisfaction."

Its signature accepts one `motif_sequence_constraint`, one `motif_column`, and
left/right context columns. `motif_type` selects per-row unique motifs or one
constant motif shared by all rows; it has no positions list, orientation mode,
spacing sweep, copy count, or motif-order parameter.

Source: `docs/api.md:230-287`; implementation signature and matching docstring
at `oligopool/motif.py:17-76`.

The implementation designs `targetcount` motifs and writes that one result
vector into `motifcol` adjacent to the specified context columns
(`motif.py:758-825`). The paper is consistent:

> "The Oligopool Calculator's design mode can introduce restriction sites,
> researcher-defined sequence motifs, and spacer regions into the oligopool
> design ... These sequences are defined using the IUPAC degenerate nucleotide
> code ... making them useful as padding elements or anchor sequences used to
> define barcode positions."

Source: Hossain et al. 2024, PDF p. 6, "Restriction Sites, Motifs, and
Spacers." This is constraint-aware component design at a fixed architecture,
not a library of placement configurations.

The repository example loops over the user's `elements_spec` and calls
`op.motif` once per named element (`examples/design-assembly-parser/
design_assembly_parser.py:180-200,574-581`). Under the global rule, an example
script repeatedly calling a single-column primitive does not become a
tool-provided placement enumerator; in any event, that script adds fixed columns
and does not sweep their locations, orders, or orientations.

### Search and disconfirmation attempt

I searched the full package, README, `docs/api.md`, `docs/docs.md`, examples, and
paper for
`motif|binding.?site|spacing|permut|combin|orient|position|insert`. I inspected
every public module in the API table of contents, especially `motif`, `spacer`,
`merge`, `join`, `revcomp`, `compress`, `expand`, `xcount`, and `verify`.
Potentially misleading hits were resolved as follows:

- `revcomp` reverses an already selected column range; it does not enumerate
  strand variants (`oligopool/revcomp.py:20-43`).
- `join` inserts DataFrame columns in an order; it does not place motifs in
  sequence coordinate combinations (`join.py:21-44`).
- `xcount` performs combinatorial barcode **read counting**, not design
  (`xcount.py:28-62`).
- IUPAC constraints design or encode a sequence in the single requested column;
  degenerate expansion is excluded by the binding definition and is not a
  position/orientation sweep.

Evidence that would have produced `partial` was a public one-motif scan across
candidate positions with fixed strand. Evidence for `yes` was a motif-list
operation generating multi-motif coordinate/strand/order combinations. I
searched public signatures, docs, examples, and the paper for those parameters;
neither exists. The fixed-column motif designer therefore remains `no`, not
`partial`.

## mutationmaker — **no** (high confidence)

### Evidence and row-5 boundary

Mutation Maker does generate combinatorial **mutation** libraries, but that is
not motif placement. Its paper describes inputs as a gene plus "a list of target
residues and mutations" and constraints including "motifs to exclude"; its
outputs are oligos carrying mutation/parental codon combinations (Hiraga et al.
2021, PDF p. 2). Later it calls the QCLM result a "combinatorial gene library"
and the PAS workflow a mutagenic library (pp. 5-8). Those are substitution-event
combinations governed by row 5, not two distinct motifs positioned across a
background.

In production code the only motif configuration is explicitly negative:

```python
class PASConfig(JsonObject):
    ...
    avoided_motifs = ListProperty(str)
```

Source: `backend/mutation_maker/pas_types.py:29-50`. The generator compiles
`config.avoided_motifs`, rejects sequences containing any, and raises "Not
possible to avoid specified combination of motifs!"
(`generate_oligos.py:233-234,323-328`). Reverse translation does the same
(`reverse_translation.py:37-78`).

### Search and disconfirmation attempt

I searched `backend/mutation_maker`, `api`, `frontend/src`, README, tests, and
the paper for
`motif|binding.?site|spacing|permut|combin|orient|insert`, excluding biological
species tables and vendored/static data when counting hits. Every production
`motif` hit was the `Motifs` restriction-enzyme/degenerate-pattern parser or an
`avoided_motifs` check; the frontend labels the same field "Avoid Motifs."
The API endpoint list exposes only SSM, QCLM, and PAS jobs
(`api/server_fastapi.py:53-65`). I separately inspected the combination code in
`mutation.py`, `generate_oligos.py`, `pas_degeneracy_recursion.py`, and the SSM
frontend: all combinations are amino-acid/codon mutations, degenerate codons,
primer choices, or oligo reaction mixtures.

Evidence that would have changed the value was a workflow input accepting
positive motifs and a coordinate/strand/spacing design space, with generated
placement members. I explicitly searched the three workflows, their schemas,
frontend forms, and paper; positive motif placement is absent. Motif avoidance
is not even a restricted placement operation, so the value is `no` rather than
`partial`.

## dnachisel — **no** (high confidence)

### Evidence and closest functionality

DnaChisel provides a positive pattern constraint, but it solves one sequence
rather than enumerating a placement library. The public API says:

> "Enforce a number of occurrences of the given pattern in the sequence."

`pattern` is singular, `occurences` is a desired count, `location` is one DNA
segment, and strand may be one or both (`EnforcePatternOccurence.py:14-48`). The
GenBank documentation similarly says:

> "You can control how many times a pattern should appear in a sequence region
> with the `@insert()` specification."

> "new occurrences of the pattern will be preferentially placed towards the
> center of the selected region."

Source: `docs/genbank/genbank_api.rst:128-150`.

The implementation makes the decisive distinction. It tries candidate starts
in center-priority order, builds a constrained single
`DnaOptimizationProblem`, and on the first solvable placement assigns
`problem.sequence = new_problem.sequence` and returns
(`EnforcePatternOccurence.py:124-170`, especially 134-160). The paper describes
the same heuristic:

> "an `insert(CGTCTC)` constraint will attempt to place the pattern 'CGTCTC' at
> different locations of the sequence"

Source: Zulkower and Rosser 2020, PDF p. 2. "Attempt" here is solver search; the
alternative locations are not emitted as distinct designed members.

Multiple `EnforcePatternOccurence` specifications can be manually combined in
one optimization problem, but that is composing single-pattern constraints and
still returns one optimized sequence. Crediting it would violate the global
rule against reconstructing the row from primitives and bookkeeping.

### Search and disconfirmation attempt

I searched all `dnachisel/`, docs, examples, README, and the paper for
`EnforcePatternOccur|motif|spacing|permut|combin|orient|insert`, inspected the
complete built-in specification index, `SequencePattern` and
`MotifPssmPattern`, the GenBank API, and every positive insertion example.
There were no permutation APIs and no library/pool motif-placement operation.
I also checked `MutationSpace.all_variants`: it can enumerate a mutation space,
but feeding separately composed insertion constraints into a generic mutation
space and then interpreting outputs as a placement library would be caller
reconstruction, not a documented motif-placement mode.

Evidence that would have yielded `partial` was a documented operation returning
the alternative placements of one pattern as distinct sequences; evidence for
`yes` was the same for two distinct patterns and position/orientation
combinations. The implementation explicitly returns after selecting one
placement, and the docs/examples never present alternative placements as
members, so neither condition was found.

## tangermeme — **yes** (high confidence)

### Evidence

Two adjacent APIs were considered. `space` alone is restricted: it accepts a
list of motifs and a two-dimensional spacing table, where each row is a spacing
combination, but uses the motifs as supplied (fixed orientations):

> "Each row in this tensor is a different combination of spacings between
> motifs and each column is the spacing between an adjacent pair of motifs."

Source: `tangermeme/space.py:34-81`, especially 64-81. Its underlying
`multisubstitute` is genuinely tool-provided multi-motif placement and returns
edited sequence tensors, but for one supplied spacing arrangement
(`tangermeme/ersatz.py:198-264`). Considered alone, this path would establish at
most `partial` under the row's fixed-orientation example.

Current tangermeme also documents `beam_substitution`, which crosses the full
threshold. It is a sequence-design operation whose input includes a motif list,
an allowed-position mask, a reverse-complement mode, edit count, beam size, and
number of returned members. Its documentation states:

> "Each round, every beam member is expanded by tiling every motif at every
> allowed position ... all resulting complete sequences are scored ... larger
> beams hedge across multiple trajectories and can recover good multi-edit
> combinations."

Source: `tangermeme/design/beam_substitution.py:44-67`.

The strand parameter is explicit:

> "Whether to augment the provided list of motifs with their reverse
> complements."

Source: lines 95-116, especially 107-109. The output is also a multi-member
sequence result, not scores or invisible candidates:

> "`n_best`: The number of sequences to return at the end, ranked from the
> lowest loss to the highest loss."

Source: lines 123-133. The implementation adds reverse complements to the motif
set (lines 184-190), tiles every motif into every allowed position for every
current beam sequence (243-269), reconstructs/deduplicates complete edited
sequences (271-290), iterates across edits, and returns
`torch.cat([...beam[:n_best]])` (305-306). With two or more distinct supplied
motifs and `max_iter >= 2`, the second expansion is over combinations of motif
identity, start position, orientation, and the preceding placement trajectory;
`n_best` returns sampled/selected concrete combinations.

The private kernel confirms what "tiling" means:

> "This function takes a motif and inserts it at all possibilities."

Source: `tangermeme/design/_substitute.py:7-16`.

This is motif placement, not row-5 mutagenesis, even though the implementation
uses substitution to implant motifs. It is also model-guided, but row 9 does not
exclude model-ranked enumeration; the operation explicitly enumerates complete
motif-placement candidates and returns multiple designed sequences.

### Disconfirmation attempt and counterfactual

I searched source, public API docs, tutorials, tests, README, changelog, and the
paper for
`multisubstitute|space|spacing|greedy_substitution|beam_substitution|motif|position|reverse_complement|all possibilities|n_best`.
I inspected `space.py`, `ersatz.py`, the complete `design/` package,
`tests/test_space.py`, `tests/test_greedy_substitution.py`, and
`tests/test_beam_substitution.py`. The tests include an explicit `n_best=3`
beam call, and the design tutorial describes beam members as complete sequences.
The paper establishes motif implantation/design as a package capability but
predates or does not detail every current beam-search parameter; the repository
is the decisive primary source for the current API.

I attempted a small behavioural run with a toy model to force two distinct
motifs into different positions, but the audit environment has no `torch`
installation (`ModuleNotFoundError: No module named 'torch'`), and installation
was prohibited. This does not make the API genuinely undeterminable: the public
docstring and straight-line implementation explicitly expose all four required
axes and the returned tensor.

Evidence that would have reduced the value was (a) only `space`, with fixed
motif orientations; (b) a greedy mode returning only one sequence; or (c)
multi-position candidates used internally but never returnable. I explicitly
searched for and found those restricted paths, but also found the documented
beam mode's reverse complements, multi-edit trajectories, and `n_best` concrete
sequence output. Therefore current tangermeme is `yes`, not `partial`.

## Row-level finding

The row discriminates on one consistent scale: three tools satisfy the test
(poolparty, MPRAnator, tangermeme) and five do not. No case is genuinely
unknown. The row separates purpose-built regulatory-grammar/design search from
variant-centric tools, fixed-architecture component design, and single-sequence
constraint repair. It does not require a replacement.

No tool landed at `partial` in the final current-version assessment. Tangermeme's
spacing API alone would have been a clean fixed-orientation `partial`, but its
current documented beam-search API independently reaches `yes`. Oligopool
Calculator's single fixed-column motif designer and DnaChisel's one chosen
pattern placement do not reach `partial`, because neither emits an enumerated
position/orientation placement set; calling either `partial` would relax the
same threshold that was applied to the other six tools.
