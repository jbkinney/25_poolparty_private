# Row audit: `exhaustive_single_scans`

**Row question.** Does the tool enumerate ALL single mutation events of a class at ALL eligible
positions, from one specification?

**Verdict on the row itself:** `keep`. The row is scoreable on one scale across all 13 tools provided
the operational test is stated medium-agnostically (see "Threshold disambiguations" below).

---

## 1. Operational test (fixed before any tool was opened)

> Does a single documented API call / CLI invocation / config specification, given a target region,
> cause the tool itself to loop over positions and produce the complete set of one-event-per-position
> variants of a defined class — every substitution, every deletion, or every codon replacement — at
> every eligible position in that region, with no requirement that the user enumerate the positions?

### Threshold disambiguations (needed because the supplied rubric is internally ambiguous)

The rubric's YES clause names *substitution* as a qualifying class, while its PARTIAL clause demotes
"only one narrow event type". Taken literally these collide, since the ANCHOR YES (VaLiAnT `snv`) is
substitution-only. Resolution actually applied, uniformly:

- **"a class"** = the complete set of alternatives of one kind at a position (all 3 non-WT bases; the
  deletion; all 19/20/63 alternative codons). A substitution-only tool therefore *can* score YES.
- **"one narrow event type"** (→ PARTIAL) = *narrower than a complete class*: a single fixed
  replacement allele/motif per position, or a degenerate-codon primer standing in for the variants
  rather than the variants themselves.
- **Medium-agnostic.** The row's verb is *enumerate*, not *emit as oligos*. A tool that constructs the
  complete variant set internally and hands back per-variant measurements (rather than sequence
  strings) still enumerates. This matters for exactly one cell (tangermeme) and is flagged there; if
  the authors want output medium scored, that belongs in a separate column, not smuggled into this row.
- **Example-script-only → at most PARTIAL** (rubric rule, applied to oligopool).
- **User must list positions → PARTIAL** (rubric rule, applied to Mutation Maker).

---

## 2. Results

| tool | prior | **verified** | one-line basis |
|---|---|---|---|
| poolparty | yes | **yes** | `mutagenize(num_mutations=1, mode='sequential')`; behaviourally verified 3L complete |
| valiant | yes | **yes** | `snv` = "Each nucleotide is replaced with all the alternatives"; `aa` = 19/codon |
| mpranator | partial | **partial** | exhaustive positional scan, but of ONE user-given motif; SNV path is random-position |
| mpradesign | no | **no** | variant set is the input VCF; no enumerator anywhere |
| oligopoolcalc | partial | **partial** | complete 3L scan exists only in `examples/library-compressor/mutant_generator.py` |
| mutationmaker | yes | **partial → CHANGED** | `mutations` is a required user-supplied list in all three workflows |
| codongenie | no | **no** | no sequence input at all; input is a set of amino acids |
| dnachisel | no | **no** | `all_variants` is a full cartesian product, unexported; `optimize()` returns one sequence |
| ledidi | no | **no** | gradient sampling; pruning scans proposed edits, not all positions |
| tangermeme | partial | **yes → CHANGED** | one call mutates every position to every base over `start:end` |
| biopython | no | **no** | `MutableSeq` per-position assignment only; no enumerator in tree |
| pydna | no | **no** | no mutagenesis module; only a SnapGene parser test fixture matches |
| seqpro | no | **no** | `__all__` contains no mutation function of any kind |

---

## 3. Cell-by-cell evidence

### 3.1 poolparty — **yes** (prior: yes; unchanged)

Read-only source and docs at
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/`, plus behavioural checks in
the existing read-only venv (`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python`).

**Documented specification** — `docs/operations/mutagenize.rst:204-206`:

> ``mode='sequential'`` with ``num_mutations=1`` enumerates every single-point
> variant in deterministic order, covering all positions and non-wild-type bases.

`docs/operations/mutagenize.rst:151-152` (region-restricted form):

> ``region`` confines all mutations to the tagged segment; flanks are returned
> unchanged. With ``mode='sequential'``, every single-base variant within the
> region is enumerated.

`docs/operations/modes.rst:58-60`:

> ``mode="sequential"`` enumerates every possibility of the design space in a deterministic
> order. … ``operation.num_states`` equals the number of possibilities.

**Behavioural verification (not taken on trust).** On `ACGTACGTAC` (L=10):

```
pp.from_seq('ACGTACGTAC').mutagenize(num_mutations=1, mode='sequential')
  -> num_states 30                      # 10 positions x 3 alternative bases
  -> 30 generated sequences; each differs from WT at exactly one position
  -> set{(position, alt_base)} == full expected 3L set   => COMPLETE SATURATION SNV: True
```

Deletion class, same one-liner form:

```
pp.from_seq('ACGTACGTAC').deletion_scan(deletion_length=1, mode='sequential')
  -> num_states 10
  -> ['-CGTACGTAC','A-GTACGTAC','AC-TACGTAC','ACG-ACGTAC','ACGT-CGTAC',
      'ACGTA-GTAC','ACGTAC-TAC','ACGTACG-AC','ACGTACGT-C','ACGTACGTA-']
```

Codon-replacement class, on the 5-codon CDS `ATGAAATTTGGGCCC`:

```
mutagenize_orf(..., mutation_type='any_codon',            mode='sequential') -> num_states 315  (= 5 x 63)
mutagenize_orf(..., mutation_type='nonsynonymous_first',  mode='sequential') -> num_states 100  (= 5 x 20)
mutagenize_orf(..., mutation_type='missense_only_first',  mode='sequential') -> num_states  95  (= 5 x 19)
mutagenize_orf(..., mutation_type='nonsense',             mode='sequential') -> num_states  15  (= 5 x 3)
```

`src/poolparty/codon_table.py:8-13` supplies the per-codon alternative counts
(`"any_codon": 63,  # All 64 - 1 (self)`).

**The flagged caveat, verified and reported plainly.** `src/poolparty/orf_ops/mutagenize_orf.py:205-207`:

```python
if mode == "sequential" and mutation_type not in UNIFORM_MUTATION_TYPES:
    raise ValueError(
        f"mode='sequential' requires a uniform mutation type, got '{mutation_type}'"
    )
```

Reproduced in the venv: `mutation_type='synonymous'` with `mode='sequential'` raises
`ValueError: mode='sequential' requires a uniform mutation type, got 'synonymous'`. I also checked the
obvious escape hatch — `reverse_translate` has **no** `mode` parameter
(`TypeError: reverse_translate() got an unexpected keyword argument 'mode'`), and
`docs/operations/reverse_translate.rst:41` describes its `codon_selection` as deterministic-best or
random sampling, not enumeration. So the **synonymous-only** class genuinely cannot be exhaustively
enumerated by label.

**Why this does not demote the cell.** (a) The row asks for exhaustiveness over *some* class; PoolParty
is exhaustive over substitutions, deletions, and codon replacements. (b) The synonymous variants are
not missing from the design space — they are a subset of the `any_codon` enumeration (63 alternative
codons per codon, i.e. all synonymous DNA variants included); what is missing is a *label* that
isolates them. (c) The "exhaustive over amino acids using one representative codon" concern applies
**identically to the anchor**: VaLiAnT's `aa` mutator is "Replace each codon with the **top-ranking
codon** of all amino acids" (README.md:356) — the same 19-per-codon, one-representative-codon design as
PoolParty's `missense_only_first`. PoolParty additionally offers the strictly larger `any_codon`
enumeration that the anchor does not have. On this row PoolParty is at or above the anchor.

**Disconfirmation.** Evidence that would have forced `partial`/`no`: (i) a default that silently
truncates — checked `docs/operations/modes.rst:143-152` and the quick-reference table at 245-266
(`sequential` + `num_states=None` → "all possibilities"); (ii) a requirement that the user list
positions — checked `scan_ops/mutagenize_scan.py:39-42` ("If None, all valid positions are used") and
`docs/operations/deletion_scan.rst:50` ("``None`` = all valid positions"); (iii) an incomplete
enumeration despite the claim — ruled out by set-equality against the full 3L SNV set above, not by
reading the docstring.

---

### 3.2 valiant — **yes** (prior: yes; unchanged) — ANCHOR

`README.md` (branch `develop`), line 297 (`snv` mutator):

> Each nucleotide is replaced with all the alternatives.

with the worked listing for target `AA` giving all six variants (`CA GA TA AC AG AT`).

Line 356 (`aa` mutator):

> Replace each codon with the top-ranking codon of all amino acids. Given the default codon table,
> this results in 19 mutated sequences for each codon mapping to an amino acid (the reference amino
> acid being excluded) and 20 for each stop codon.

Line 281 (parametric deletion, whole-target scan):

> Non-overlapping stretches of nucleotides of a given length are deleted starting from a given offset.
> No partial deletions are performed at the end of the target regions.

**One specification** = the mutator label in the targeton file; README.md:678 shows a real row whose
mutator column is `snv,1del,snvre` for a transcript span. Positions are never listed by the user.

**Disconfirmation.** Looked for a per-position input requirement (there is one, but only for the
separate `custom` mutator sourced from VCFs, README.md:270 "Variants imported from VCF files … are
labelled as `custom`") and for any statement that `snv`/`aa` cover only part of the target — none.

---

### 3.3 mpranator — **partial** (prior: partial; unchanged, but on primary-source grounds)

Repo `hemberg-lab/MPRAnator` (1 commit). Files read in full: `oligo.py`, `part1.py`, `part3.py`,
`parseSNPs.py`, `myfunctions.py`, `mycustom.py`, `iliasApp/views.py`, `iliasApp/forms.py`,
`iliasApp/templates/iliasApp/docs.html`.

**What IS exhaustive over positions** — `oligo.py:39-40`:

```python
positionsForMotifsL = [p for p in
                       itertools.product(range(distanceFromLeftEdge, len(BackgroundS)), repeat=len(motifsL))]
```

with the emission filter at `oligo.py:70` and the per-position output at `oligo.py:74-76`:

```python
if sorted(positionsL)[0] % frequencyOfInsertion == 0:
    ...
    allResults.append({"header": ..., "sequence": "".join(backgroundL), ...})
```

So with `frequencyOfInsertion = 1` the tool itself walks every eligible position of the background and
emits one sequence per position. **But the event is substitution of a window by ONE user-supplied
motif** — the form field is `motifS = FastaTextField(label="Enter your motifs")` (`forms.py:404`) and
the interval field is `frequencyOfInsertion = PositiveIntegerField(label="Interval of substitution of
motifs")` (`forms.py:428`). That is one narrow event type, not a class of alternatives → PARTIAL.

**Why not YES.** No single-nucleotide/codon saturation path exists. `part3.py:32-36` (the
"Transmutation" mutator) chooses positions **randomly**:

```python
for i in range(noOfPos):
    posChosen = rd.choice(LengthList)
    while sequenceS[posChosen] == "N" or posChosen in positionsChosenL:
        posChosen = rd.choice(LengthList)
```

and the substituted base is also random (`part3.py:44`, `therandomchoice = rd.choice(mapperDict[stringL[i]])`).
`iliasApp/templates/iliasApp/docs.html:142` describes it as "Mutate random positions in the input
sequence."; `views.py:142` labels the page "Creates a specified number of random mutations on a
sequence". The only other variant source is a user VCF: `parseSNPs.py:9` `def read_SNPs(SNP_file_path)`,
form field `SnpS = SnpTextField(label="Enter your SNPs (VCF Format)")` (`forms.py:502`).

**Disconfirmation.** Checked the deployed documentation at `genomegeek.com/MPRA/documentation/` for a
saturation option: it lists five modes (MPRAs-Motifs, MPRAs-SNPs, Transmutation, PWM Seq-Gen,
Documentation) and contains no single-/double-nucleotide scanning, "all positions", or saturation
language. Grepped every `.py` and the forms/templates for `scan|mutat|single nucl|double nucl` — the
only hits are the random Transmutation path and dplyr-style incidental words.

---

### 3.4 mpradesign — **no** (prior: no; unchanged)

`andrewGhazi/mpradesigntools` `R/processVCFfast.R` (1459 lines) read; `README.md` and
`andrewGhazi/designMPRA/ui.R` grepped.

The entry point takes the variant set as input — `R/processVCFfast.R:1099`:

```r
processVCF = function(vcf,
                      nper,
                      upstreamContextRange,
                      downstreamContextRange,
                      fwprimer,
                      revprimer, ...)
```

and the roxygen `@return` (`R/processVCFfast.R:1068-1071`) is "A list of two tibbles. The first, named
'result', is a tibble containing the labeled MPRA sequences. The second, named 'failed', is a tibble
listing the SNPs that are not able to have MPRA sequences generated". `spreadAllelesAcrossRows`
(`:5`) only splits multi-allelic VCF records; `fix_indels` (`:1358`) and `spread_and_fix_indels`
(`:1418`) normalise REF/ALT strings.

**Disconfirmation.** Grepped `R/processVCFfast.R`, `README.md`, and `designMPRA/ui.R` for
`mutagen|saturat|scan|all possible|every position|substitut` — zero feature hits (the sole `substitut`
hits are unrelated prose in the Shiny UI, and every one of the 40+ `mutate(` occurrences is dplyr's
dataframe verb, e.g. `:402 res %<>% mutate(sequence = paste0(fwprimer, …))`). The one function that
does change bases, `randomly_fix` (`:124`), exists to repair aberrant restriction sites and does so
randomly (`altered_pattern = map2_chr(pos_to_change, possible_alleles, …)`), not exhaustively.

---

### 3.5 oligopoolcalc — **partial** (prior: partial; unchanged, and the flagged reason confirmed)

**The exhaustive code exists, and only in `examples/`.**
`examples/library-compressor/mutant_generator.py:43-79`:

```python
def generate_single_mutants(
    sequence: str,
    include_wildtype: bool = True,
) -> Generator[Tuple[str, str, int, str, str], None, None]:
    '''
    Generate all single-nucleotide mutants of a DNA sequence.

    For a sequence of length L, this generates up to 3*L mutants.
    ...
    seq_list = list(sequence)
    for pos in range(len(sequence)):
        ref_base = sequence[pos]
        for alt_base in DNA_BASES:
            if alt_base != ref_base:
                ...
                yield (variant_id, mutant_seq, pos + 1, ref_base, alt_base)
```

That is a genuine complete 3L scan driven by the tool. It is example code, reached via
`generate_variant_dataframe(sequence, strategy='single')` (`:200-233`), not via the `op.*` API.

**The package has no mutagenesis module.** GitHub contents listing of `oligopool/` returns exactly:
`__init__, acount, background, barcode, cli, compress, descriptions, expand, final, index, inspect,
join, lenstat, merge, motif, pack, pad, primer, revcomp, spacer, split, verify, xcount` (+ `base/`).
Fetched `docs/api.md` (76 KB, all 22 module signatures) and searched it: **zero** occurrences of
"mutation", "substitution", "saturation", or "SNV"; the only enumerator is `expand`, which expands an
IUPAC-degenerate string (and therefore produces the combinatorial 4^N set for an all-N region, never a
one-event-per-position set).

Per the rubric's example-script rule → PARTIAL.

**Disconfirmation.** Would have moved to YES: a `op.mutate`/`op.saturate` module in the package, or a
documented API-level saturation recipe. Checked the package file listing, `docs/api.md`, and
`examples/` (which contains only `OligopoolCalculatorInAction.ipynb`, `README.md`, and four
sub-directories, of which `library-compressor/` holds `README.md`, `mutant_generator.py`,
`run_degenerate_demo.py`). GitHub code search is auth-walled, so the check was done against the fetched
file listings and files themselves.

---

### 3.6 mutationmaker — **partial** (prior: **yes** → CHANGED)

Live author-maintained fork `ra100/Mutation_Maker`. Files read: `backend/mutation_maker/ssm.py`,
`ssm_types.py`, `qclm_types.py`, `pas_types.py`, `mutation.py`, `README.md`,
`frontend/src/scenes/SSM/components/SSMForm.tsx`, `frontend/src/shared/components/MutationsInput.tsx`.

**Every workflow requires the user to enumerate the sites.**

`backend/mutation_maker/ssm_types.py:213` and `:217-218`:

```python
class SSMInput(JsonObject):
    sequences = ObjectProperty(SSMSequences, required=True)
    config = ObjectProperty(SSMConfig, required=True)
    mutations = ListProperty(str, required=True)
    degenerate_codon = StringProperty(required=True, default="NNS")

    def parse_mutations(self, goi_offset):
        return [parse_codon_mutation(mutation, goi_offset) for mutation in self.mutations]
```

`backend/mutation_maker/qclm_types.py:99-105` — the docstring names the source explicitly:

```python
    mutations = ListProperty(str, required=True)

    def parse_mutations(self, goi_offset: int) -> List[MutationSite]:
        """
        Parses the user input in the format "E32W E32L E49K" and produces multi-amino
        mutations in the format of "E32WL E49K"
        """
```

`backend/mutation_maker/pas_types.py:117-120`:

```python
    """ List of mutations at a given amino acid position """
    position = IntegerProperty(required=True)
    ...
    mutations = ListProperty(PASMutation, required=True)
```

`backend/mutation_maker/mutation.py:28-35` — a position integer is mandatory, so there is no wildcard:

```python
def parse_codon_mutation(mutation_string, gene_of_interest_offset=0) -> "AminoMutation":
    try:
        one_based_codon_position = int(mutation_string[1:-1])
    except Exception:
        raise ValueError("Position must be positive number")
```

Frontend confirms the same: `SSMForm.tsx:218-231` renders a free-text area tooltipped "Enter mutations
in following format [Amino Acid Codon][Location][X]" validated by
`pattern: { value: /^\s*([A-Z][0-9]*[A-Z])(\s+[A-Z][0-9]*[A-Z]\s*)*$/i, message: 'Invalid mutation format' }`
with `required: 'Mutations are required'` — a whitespace-separated explicit list, no range or
all-positions syntax. `MutationsInput.tsx` merely reflows that list into columns.

**Why partial and not no.** Per listed site, SSM does deliver a full amino-acid saturation, because the
mutagenic primer carries a degenerate codon (`degenerate_codon = StringProperty(required=True,
default="NNS")`, `ssm_types.py:214`). So the *event* is saturating; the *positional* exhaustiveness is
the user's job — exactly the rubric's PARTIAL clause ("exhaustiveness is the user's responsibility
(they must enumerate positions themselves)"). Note also that what is emitted is primer pairs
(`create_primer_output`, `ssm.py:775`), not variant sequences.

**Disconfirmation.** Would have kept YES: any code path deriving the site list from the gene, or a UI
"all positions" control. Searched `ssm.py` (all `def`s and every `mutations` reference), `ssm_types.py`,
`qclm_types.py`, `pas_types.py` for any auto-enumeration; searched `SSMForm.tsx` and
`MutationsInput.tsx` for a select-all/range control; grepped `README.md` for
`SSSM|MSDM|PAS|saturat|mutagenesis|all posi|each posi|every` — the only workflow mentions are naming
(lines 6, 230, 287), with no auto-scan claim. `frontend/src/scenes/SSM/components/amino.ts` computes
degenerate-codon coverage, not site lists.

---

### 3.7 codongenie — **no** (prior: no; unchanged)

Canonical repo redirects `synbiochem/CodonGenie` → `neilswainston/CodonGenie` (verified via the GitHub
API `full_name`). Files read: `codon_genie/codon_selector.py`, `main.py`,
`static/codonGenie/codonGenie.ctrl.js`, `README.md`.

**There is no sequence input, hence no positions.** `codon_genie/codon_selector.py:72-79`:

```python
    def optimise_codons(self, amino_acids, organism_id):
        ...
        req_amino_acids = set(amino_acids.upper())
        codons = [CODONS[amino_acid] for amino_acid in req_amino_acids]
        results = [self.__analyse(combo, organism_id, req_amino_acids)
                   for combo in itertools.product(*codons)]
```

REST surface, `main.py:61-72`:

```python
@app.route('/codons')
def get_codons():
    '''Gets codons.'''
    if 'aminoAcids' in request.args:
        codons = _CODON_SELECTOR.optimise_codons(request.args['aminoAcids'], request.args['organism'])
    else:
        codons = _CODON_SELECTOR.analyse_codon(request.args['codon'], request.args['organism'])
```

Client query object, `static/codonGenie/codonGenie.ctrl.js`:
`self.query = {"mode": "aminoAcids", "aminoAcids": [], "codon": ""};`

The tool maps {amino acids} → degenerate codons, or analyses one ambiguous codon. It never receives a
target region, so "all eligible positions" is undefined for it.

**Disconfirmation.** Would have moved off `no`: a sequence or CDS parameter on any route, or a
positional loop. Enumerated all four routes in `main.py` (`/`, `/organisms/`, `/organisms/<term>`,
`/codons`), read the whole of `codon_selector.py`, and checked the Angular controller for any sequence
field. `seq_utils.py` and `ncbi_tax_utils.py` are taxonomy/codon-usage helpers. This is a genuine
absence, not a failure to find.

---

### 3.8 dnachisel — **no** (prior: no; unchanged)

The near-miss is `MutationSpace.all_variants`, and it is the wrong shape.
`dnachisel/MutationSpace/MutationSpace.py:132` and `:159-164`:

```python
    def all_variants(self, sequence):
        """Iterate through all sequence variants in this mutation space."""
        ...
        for variants in itertools.product(*variants_slots):
            new_sequence[choice_start:choice_end] = encoded_segment
            for (start, end), variant in variants:
                new_sequence[start:end] = variant
            yield new_sequence.decode()
```

`itertools.product` over every mutation slot means each yielded sequence carries a choice at **every**
slot simultaneously — the combinatorial space, not the one-event-per-position set. It is also not part
of the public API: `dnachisel/__init__.py`'s `__all__` (lines 75+) lists `DnaOptimizationProblem`,
`CircularDnaOptimizationProblem`, the ~20 `Avoid*`/`Enforce*` specifications, patterns, and biotools —
`MutationSpace` and `MutationChoice` are absent.

The nearest specification-level feature is a constraint, not an enumeration —
`dnachisel/builtin_specifications/EnforceChanges.py:17-31`:

> Specify that some locations of the sequence should be changed/different. … 3 nucleotide changes, or
> EnforceChanges(minimum_percent=70) to enforce … EnforceChanges(amount=10) to aim at a … number of
> nucleotide changes

and DnaChisel's execution model returns a single optimised sequence per problem.

**Disconfirmation.** Would have moved off `no`: a documented library/variant-set output, or a
single-mutant enumerator. Grepped `README.rst` for `saturat|library|all variants|enumerat|scan` — the
only "library" hits describe DnaChisel itself as "a Python library". Enumerated the full non-test
module list (69 modules); the closest names are `EnforceChoice`, `EnforceChanges`, `AvoidChanges`, all
constraint specifications. Read `all_variants` in full rather than trusting its docstring.

---

### 3.9 ledidi — **no** (prior: no; unchanged)

`docs/index.rst` (official docs, in-repo) states the mechanism:

> Ledidi works by phrasing the design process as a continuous optimization problem and then solving
> this problem using off-the-shelf techniques. … Ledidi circumvents this challenge by learning a
> continuous weight matrix from which edits to an initial sequence are sampled

The one function that iterates over positions is post-hoc and scoped to the edits already made —
`ledidi/pruning.py:25-33`:

```python
def greedy_pruning(model, X, X_hat, threshold=1, target=None, verbose=False):
	"""A method for pruning edits to remove those that are irrelevant.

	This method will greedily go through all of the proposed edits and evaluate
	the effect of removing them, one at a time. As a greedy method, this will
	iteratively scan over all edits and remove the one with the smallest change
	in model output …
```

"all of the **proposed edits**" ≠ all eligible positions × all alternative bases.

**Disconfirmation.** Would have moved off `no`: an exhaustive-scan entry point. Enumerated the whole
package (`ledidi/__init__.py`, `ledidi.py`, `losses.py`, `plot.py`, `pruning.py`, `wrappers.py`) and the
docs tree (`getting_started.rst`, `parameters.rst`, `input_output.rst`, `faq.rst`, four tutorials); the
design loop is sampling-based and `pruning` is the only positional iterator. Note ledidi depends on
tangermeme, but that dependency is "for input validation and … the sequence and model utilities used
throughout the tutorials" (`docs/index.rst`), which does not make ISM a ledidi feature.

---

### 3.10 tangermeme — **yes** (prior: **partial** → CHANGED)

One call, tool-side loop, complete substitution class over a user-chosen window.

`tangermeme/saturation_mutagenesis.py:66-79` (signature) and `:126-133` (window semantics):

```python
def saturation_mutagenesis(
	model: torch.nn.Module,
	X: torch.Tensor,
	args: tuple | None = None,
	start: int = 0,
	end: int = -1,
	...
	end: int, optional
		The end of where to make perturbations to the sequence. … `end=-1` maps to `length` and the
		*last* nucleotide is included. … Default is -1, meaning the entire sequence.
```

The enumeration itself, `tangermeme/saturation_mutagenesis.py` (inside the per-example loop):

```python
		ref = X[i, :, start:end]
		identity = (ref == 1) & (ref.sum(dim=0, keepdim=True) == 1)
		edits = ~identity.reshape(-1)

		edit_chars, edit_positions = torch.where(~identity)
		edit_positions = edit_positions + start
		n_edits = edit_chars.shape[0]

		X_ = X[i].repeat(n_edits, 1, 1)
		rows = torch.arange(n_edits)
		X_[rows, :, edit_positions] = 0
		X_[rows, edit_chars, edit_positions] = 1
```

`torch.where(~identity)` over the `(alphabet, window)` grid is precisely "every non-wild-type character
at every position"; the accompanying comment confirms the only omissions are identity edits, which are
filled back in from `y0` ("it is skipped here and filled from y0 below"). Official in-repo doc,
`tangermeme/_skills/data/references/saturation_mutagenesis.md`:

> In-silico saturation mutagenesis (ISM) mutates every position to every base and measures the change in
> the model's output.

and, on windowing, "Windowed ISM is **bit-identical** to running the full sequence and slicing (edits
are independent)".

**Why this is a promotion, and the caveat that goes with it.** Under the medium-agnostic reading fixed
in §1, this is a complete positional scan of the substitution class from one specification — the same
class and the same completeness as the ANCHOR (VaLiAnT `snv`), so scoring it below the anchor would be
understating a competitor. The honest caveats, which belong in adjacent columns rather than this one:
(a) what is returned is per-edit model output / attribution
(`y_hat` shaped `(batch, len(alphabet), end-start, …)`), not the variant sequences — the full
edit-distance-one tensor `X_` is materialised but local; (b) substitution only — tangermeme has no
deletion scan and no codon/reading-frame model (`ersatz.py`'s `insert`/`substitute`/`delete` act at a
single caller-given position).

**Disconfirmation.** Would have kept `partial`: an alphabet or position subset that the function
silently skips, or a requirement that the caller supply positions. Read the function body in full and
checked that `start`/`end` default to the whole sequence, that `positions` is not a parameter, and that
the only skipped cells are identity edits (reconstructed from `y0`). Also enumerated the package (28
modules) and `docs/api/*` to confirm no *other* module is the intended scanner and that
`saturation_mutagenesis` is a first-class public API (`docs/api/saturation_mutagenesis.rst`,
`docs/tutorials/Tutorial_B4_Saturation_Mutagenesis.ipynb`), not an example script — which is what
separates it from oligopool.

---

### 3.11 biopython — **no** (prior: no; unchanged)

Assessed as the current package, per instruction.

Whole-repository path search for `mutag|mutat|satur|variant|ism` over the recursive git tree returned
**three** paths, none of them library code: `Bio/Entrez/DTDs/NCBI_Organism.dtd`,
`Bio/Entrez/DTDs/NCBI_Organism.mod.dtd`, `Tests/MAF/length_coords_mismatch.maf`.

The editing facility is manual and per-position. `Doc/Tutorial/chapter_seq_objects.rst:876-891`:

```
   >>> from Bio.Seq import MutableSeq
   >>> mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
   ...
   MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')
```

i.e. one assignment per edit, written by the user. `Bio/Seq.py`'s public method list is string
operations plus `translate/complement/reverse_complement/transcribe/back_transcribe/join/replace` and
the `MutableSeq` mutators `append/insert/pop/remove/reverse/extend` — no enumerator.
`Bio/SeqUtils/__init__.py` exposes `gc_fraction, GC123, GC_skew, xGC_skew, nt_search, seq3, seq1,
molecular_weight, six_frame_translations`.

**Disconfirmation.** Would have moved off `no`: any module/function producing a variant set. Searched
the full repo tree by path; grepped `Doc/Tutorial/chapter_seq_objects.rst`, `chapter_cookbook.rst`, and
`chapter_seq_annot.rst` for `mutagen|saturat|all possible substitut|every position` (only `MutableSeq`
hits); enumerated `Bio/SeqUtils/`. Searched, absent.

---

### 3.12 pydna — **no** (prior: no; unchanged)

Full module list of `src/pydna/` (40 modules) contains no mutagenesis module: alphabet, amplicon,
amplify, assembly, assembly2, codon, common_sub_strings, contig, cre_lox, crispr, design, dseq,
dseqrecord, fakeseq, fusionpcr, gateway, gel, genbank, genbankfixer, ladders,
oligonucleotide_hybridization, opencloning_models, parsers, primer, primer_screen, readers,
recombinase, seq, seqrecord, sequence_picker, sequence_regex, snapgene_history_parser, tm, utils, …

Grepped `design.py`, `dseqrecord.py`, `utils.py`, `codon.py`, `seq.py` for
`mutagen|mutate|saturat|all possible|every position|substitut`: two hits total, both unrelated prose
(`design.py:97`, `:126`, "can be substituted for a custom made function", about the Tm callable).
Repo-wide path search for `mutag|mutat|satur|scan` returned exactly one file:
`tests/snapgene_history_files/mutagenesis_primer_1primer.dna` — a SnapGene-history parser test fixture.

`docs/example_gallery.rst` does list "02 Primer directed mutagenesis with pydna", but the link target is
`github.com/hiyama341/teemi/.../02_primer_directed_mutagenesis.ipynb` under the heading "Bioengineering/SynBio
projects where you can apply pydna" — an external third-party notebook, not pydna API, and in any case
primer-directed (single site), not a scan.

**Disconfirmation.** Would have moved off `no`: a `mutate`/`saturate` function or a QuikChange-style
scanner. Enumerated all 40 modules, grepped the five most likely, path-searched the whole tree, and read
the docs gallery and docs README for feature claims. Searched, absent.

---

### 3.13 seqpro — **no** (prior: no; unchanged)

The complete public API, `python/seqpro/__init__.py`:

```python
from ._analyzers import gc_content, length, nucleotide_content
from ._encoders import decode_ohe, decode_tokens, ohe, pad_seqs, tokenize
from ._modifiers import bin_coverage, jitter, k_shuffle, random_seqs, reverse_complement
```

with `__all__ = [AA, DNA, RNA, AminoAlphabet, NestedStr, NucleotideAlphabet, PathLike, SeqType,
StrSeqType, alphabets, bed, bin_coverage, cast_seqs, decode_ohe, decode_tokens, gc_content, gtf,
jitter, k_shuffle, length, nucleotide_content, ohe, pad_seqs, rag, random_seqs, reverse_complement,
tokenize, transforms]`. No mutation function of any kind.

`python/seqpro/_modifiers.py` public functions: `reverse_complement` (:14), `k_shuffle` (:39),
`bin_coverage` (:111), `jitter` (:144), `random_seqs` (:281).
`python/seqpro/transforms/augmentation.py` classes: `Sequential, Random, ReverseComplement, KShuffle,
Jitter, Tokenize`. `python/seqpro/experimental/_experimental.py`: `edit_distance, kmer2seq, seq2kmer`.

**Disconfirmation.** Would have moved off `no`: an ISM/saturation helper (plausible, given SeqPro's
ML4GLand sibling packages). Grepped `_analyzers.py`, `_cleaners.py`, `_coords.py`, `_encoders.py`,
`_utils.py` for `mutat|mutag|satur|substitut` — zero hits; enumerated `_modifiers.py`,
`transforms/augmentation.py`, and `experimental/_experimental.py` by symbol. Searched, absent.

---

## 4. Changes from the prior per-tool pass

1. **mutationmaker: yes → partial.** The prior `yes` cannot be sustained against the source: `mutations`
   is `required=True` on `SSMInput`, `QCLMInput`, and `PASMutationSite`, `parse_codon_mutation` demands
   an integer position, and the frontend validates a literal whitespace-separated list. Nothing in the
   backend or frontend derives the site list from the gene, so positional exhaustiveness is the user's
   responsibility — the rubric's explicit PARTIAL condition. (It stays above `no` because each listed
   site does get a full-saturation degenerate codon.)
2. **tangermeme: partial → yes.** `saturation_mutagenesis` performs the complete positional scan itself
   over `start:end` (default: the entire sequence), constructing every non-wild-type character at every
   position; its own documentation says it "mutates every position to every base". That is the same
   class and completeness as the anchor. Demoting it would understate a competitor on the row as worded.

## 5. Notes for the authors

- **PoolParty is not demoted, and the reason matters for the rebuttal.** The `mode='sequential'`
  refusal for non-uniform `mutation_type` is real and reproducible, but it removes only the
  *synonymous-only label*, whose members are still enumerated by `any_codon` (63 alternatives per
  codon). And the "amino-acid-level, one representative codon" pattern is exactly what VaLiAnT's `aa`
  mutator does ("the **top-ranking codon** of all amino acids", README.md:356). PoolParty is at or
  above the anchor on this row; the honest caveat is a missing label, not a missing capability.
- **This row does not score output medium.** Three cells (mutationmaker → primers, tangermeme →
  attribution tensors, oligopoolcalc → example DataFrame) have completeness properties that differ from
  their delivery properties. If the table needs the distinction, add a separate column ("emits
  synthesis-ready oligo sequences") rather than overloading this one; overloading is a likely source of
  the original rater variance.
- **The `yes` set on this row is {poolparty, valiant, tangermeme}** — three tools, from three different
  communities, meeting the same bar. That is a more credible row than one where only the anchor and the
  paper's own tool clear it.
