# Verified audit — row `higher_order_combinatorial`

**Row question.** Does the tool itself enumerate or sample pairwise-and-higher combinations of
mutations co-occurring in one output sequence?

**Auditor note.** One auditor, one threshold, all 13 tools. Primary sources only (file:line or exact
quoted line). The previous per-tool memos were used only to decide where to look and are not cited as
evidence anywhere below.

---

## 1. Operational test (fixed before any tool was opened)

Recorded verbatim in `verified/OPTEST_higher_order_combinatorial.txt` (written before the first tool
was inspected):

> For each tool, search its own source for a construct that, from ONE user specification in ONE
> invocation, ENUMERATES (`itertools.product` / `itertools.combinations` / nested per-site loops over
> per-site variant sets) or RANDOMLY SAMPLES which >=2 independently-positioned sequence-altering
> events co-occur inside a SINGLE output sequence. The decisive, re-runnable signature is CARDINALITY
> GROWTH: the number of distinct output sequences the tool can emit from that one specification grows
> multiplicatively / binomially in the number of candidate sites (`prod_i |V_i|`, or
> `comb(n,k)*|A|^k`), not linearly (one event per output, `sum_i |V_i|`).

Exclusions applied (restating the row's rubric operationally):

| id | excluded |
|---|---|
| X1 | hand-authored multi-variant input — the user supplies the combination |
| X2 | uniform context edits applied to every oligo (PAM protection, background variants, restriction-site repair) |
| X3 | one multi-base single event (a 3 bp parametric deletion, a whole-codon swap) |
| X4 | combinations obtained by re-running the tool and merging outputs |
| X5 | a generic sequence container with no enumeration logic of its own |

### 1.1 Declared sub-decisions (the row as handed to me does not settle these; I settle them here and apply them identically to PoolParty and to every competitor)

**(a) Combinatorial MOTIF placement does NOT count for this row.** A motif substitution is a
sequence-altering event, but the row set already contains a dedicated row,
`combinatorial_motif_place` ("Combinatorial placement of motifs (multiple motifs, positions,
orientations, permutations)"). Scoring motif placement here would make the two rows redundant and
would let a tool bank the same capability twice. This row is therefore restricted to combinations of
**mutations of a reference sequence** (substitutions, indels, codon changes at >=2 distinct sites).
*Consequence, stated for transparency:* this decision costs PoolParty its `insert_from_motif` /
`insert_kmers` evidence and costs MPRAnator its `oligo.py` / `part1.py` motif-product evidence.
Both still score `yes` on strictly mutational engines, so the decision changes no cell — I verified
this rather than assumed it.

**(b) Degenerate/IUPAC template expansion does NOT count for this row.** Expanding `ACGTRCGTRCGT`
into its cartesian product does put >=2 varying positions in one output sequence, but the user has
hand-authored the entire combinatorial specification inline; the tool performs a pure decompression
and never selects which sites co-vary, and there is no reference sequence or notion of "mutation" in
the operation. This is the nearest thing to exclusion X1, and the row set again scores it separately
(`degenerate_iupac_codons`). **This decision is forced by consistency, not preference:** the identical
cartesian-product-over-IUPAC-positions function exists in Oligopool Calculator
(`core_degenerate.py:858`) *and* in pydna (`primer_screen.py:178`). Any standard that gave Oligopool
"partial" for it would have to give pydna "partial" too, which no referee would accept. Both are
therefore scored on their non-IUPAC evidence.

**(c) An objective-driven optimiser/sampler that emits multi-edit sequences scores `partial`, not
`yes`.** Where the tool selects a set of co-occurring edits by searching against a model or a
constraint set, the higher-order output is real, but the combinatorial space is neither declared by
the user nor enumerated/uniformly sampled by the tool, and the user cannot ask for "the doubles".
This is the rubric's PARTIAL ("combinations possible but only across a restricted axis"). Affects
DnaChisel, Ledidi, tangermeme. **The hinge is recorded per cell** — a referee who reads "or sample"
in the row question purely mechanically would promote all three to `yes`; that reading is stated
explicitly in each `disconfirmation` so the authors can flip them deliberately rather than by
accident. Note also that the row set already carries `ml_model_in_loop` for Ledidi/tangermeme.

---

## 2. Result

| tool | prior | **verified** | one-line basis |
|---|---|---|---|
| poolparty | yes | **yes** | `comb(n,k)` x alphabet enumeration; behaviourally verified, 405/405 outputs at Hamming 2 |
| valiant | no | **no** | closed 7-member single-event mutator enum; `OligoSeq` holds exactly one variant |
| mpranator | yes | **yes** | power set of in-window SNVs applied cumulatively to one sequence |
| mpradesign | no | **no** | oligo count is *linear* in #SNPs x 2 alleles, per its own README |
| oligopoolcalc | partial | **no (CHANGED)** | zero mutagenesis in the package; only IUPAC decompression (decision (b)) |
| mutationmaker | yes | **yes** | cartesian product over per-site codon dicts; every PAS oligo carries a codon at every site |
| codongenie | no | **no** | product is within a single 3 bp codon (X3); emits a codon, not a sequence |
| dnachisel | partial | **partial** | `MutationSpace.all_variants` / `apply_random_mutations(n,·)` are real but are solver primitives |
| ledidi | partial | **partial** | samples multi-edit sets, but oracle-driven; no declarable combinatorial space |
| tangermeme | partial | **partial** | greedy/beam accumulate k substitutions into one sequence; `saturation_mutagenesis` is explicitly single-base |
| biopython | no | **no** | no mutagenesis code of any kind in `Bio/` |
| pydna | no | **no** | all combinatorics is fragment assembly; IUPAC expander excluded by (b) |
| seqpro | no | **no** | complete public API contains no mutation operation |

**One change: `oligopoolcalc` partial -> no.**

---

## 3. Cell-by-cell evidence

### 3.1 poolparty — **yes**

Source (`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/src/poolparty/base_ops/mutagenize.py:351-363`):

```python
            num_combinations = comb(num_positions, self.num_mutations) * (
                alpha_minus_1**self.num_mutations
            )
            cache = []
            for positions in combinations(range(num_positions), self.num_mutations):
                num_mut_patterns = alpha_minus_1**self.num_mutations
```

The tool selects both *which* set of `num_mutations` positions co-mutate (`combinations(...)`) and
the substitution pattern at each. Non-uniform branch, same file:367-369:

```python
            for positions in combinations(range(num_positions), self.num_mutations):
                counts_for_positions = [mutation_counts[p] for p in positions]
                num_mut_patterns = prod(counts_for_positions)
```

**Behavioural verification** (read-only venv, `PYTHONDONTWRITEBYTECODE=1
/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/.venv/bin/python`):

```
pp.mutagenize(pp.from_seq("ACGTACGTAC"), num_mutations=2, mode="sequential")
  num_states: 405        expected comb(10,2)*3**2 = 405
  Hamming-distance-to-WT distribution over all 405 emitted seqs: Counter({2: 405})

pp.mutagenize_orf(pp.from_seq("ATGAAACCCGGGTTTTAA"), num_mutations=2, mode="sequential")
  num_states: 5415       differing-codon distribution: Counter({2: 5415})

pp.deletion_multiscan(<30 nt>, deletion_length=2, num_deletions=2, mode="sequential")
  num_states: 378        '-' marker count per seq: Counter({4: 378})   # 2 deletions x 2 bp
```

Every emitted sequence carries exactly two independent mutation events. Three independent
higher-order engines exist (`base_ops/mutagenize.py`, `orf_ops/mutagenize_orf.py:427-439`,
`multiscan_ops/*`; also `base_ops/recombine.py:310` for breakpoint combinations).
`multiscan_ops/deletion_multiscan.py:31` docstring: `"""Delete segments at multiple positions
simultaneously."""` with parameter `num_deletions : Integral / Number of simultaneous deletions to
make.`

*Disconfirmation.* A demotion would require (i) the enumeration to be motif/IUPAC-only (decisions (a),
(b)) — it is not; `mutagenize` operates on a reference with `MUTATIONS_DICT[wt]` targets
(`mutagenize.py:426`); (ii) the emitted sequences to carry one event each — measured, they carry two;
(iii) the combinatorics to live only in a notebook — it is in `base_ops/`, exported on `DnaPool`
(`pool_mixins`), and reachable as `pp.mutagenize(...)`. Checked all three. No demotion.

### 3.2 valiant — **no**

Whole-source search of `src/valiant/` at `develop` (merge of tag 4.0.0, commit `8796cc1`):

```
$ grep -rn "\bproduct(\|combinations(\|permutations(" --include=*.py src/
   (no output)
$ grep -rn "itertools" --include=*.py src/
   src/valiant/utils.py:23:from itertools import groupby
   src/valiant/mutators/snv.py:20:from itertools import chain
   src/valiant/loaders/gtf.py:23:from itertools import chain
   src/valiant/queries.py:19:from itertools import chain
   src/valiant/uint_range.py:23:from itertools import chain
```

There is no cartesian/binomial construct anywhere in the package. The mutator set is a closed enum
(`src/valiant/mutator_type.py:22-30`):

```python
class MutatorType(str, Enum):
    # Parametric span
    DEL = 'del'
    # Single-nucleotide span
    SNV = 'snv'
    # Codon span
    SNV_RE = 'snvre'
    IN_FRAME = 'inframe'
    ALA = 'ala'
    STOP = 'stop'
    AA = 'aa'
```

`DEL` is a parametric span (X3), the four codon-span mutators are one codon per variant (X3). Each
oligo carries exactly one variant (`src/valiant/oligo_seq.py:35,44-45`):

```python
class OligoSeq(Sized, Generic[VariantT]):
    variant: VariantT
...
        s = alter_seq(seq, var).s
```

and `alter_seq` (line 30) applies a single `Variant`: `return seq.alter(variant.ref_range,
variant.is_insertion, variant.alt)`. Also `grep -rn -i "combinator\|pairwise\|higher.order\|epistas\|
double.mutant\|multi.variant" README.md CHANGELOG.md src/` -> no output.

*Disconfirmation.* A promotion would require a multi-variant oligo builder. The only plural applier
is `seq_converter.apply_variants(ref_seq, alt_length, variants)`, and `grep -rn "apply_variants"
--include=*.py src/` returns exactly two call sites: `sge_proc.py:258` (`ppes_in_range`, PAM-protection
edits) and `sge_proc.py:406` (`bg_vars`, background variants). Both are exclusion X2 — uniform edits
applied to every oligo, not enumerated combinations. Searched and absent.

### 3.3 mpranator — **yes** (basis corrected: SNV power set, not motif placement)

`parseSNPs.py:20-26` and `:31,36,53-54,62,84,99`:

```python
def generateCombinations(theList):  # Used to make all combinations
    combinationsStore = []
    upTo = len(theList) + 1
    for i in xrange(1, upTo):
        for acomb in itertools.combinations(theList, i):
            combinationsStore.append(acomb)
    return combinationsStore
...
def make_sequence_copies(SNPs, NamesL, SequencesL, CombinationsB):
    if CombinationsB == True:
...
                combs=generateCombinations(All_groups[f])
                for comb in combs: #for every SNP combination
...
                    for snped in comb:
...
                                        for index in range(len(Sequenced)):
...
                                                    Sequenced+=[Sequenced[index][:position_to_change] + str(change) + Sequenced[index][position_to_change + 1:]]
```

The input VCF lists *single* SNVs (`i[0..4]` = chrom, pos, id, REF, ALT, read at `:43-46`); the tool
groups the SNVs falling inside each oligo window (`:47-50`) and then expands. Because the inner
`for index in range(len(Sequenced))` re-applies each SNV to every sequence already accumulated, the
list doubles per SNV — a power set, not a set of singles.

**Behavioural verification** — the SNV branch (difference==0 path) transcribed line-for-line into
Python 3 (the original is Python 2: `xrange`, `print` statements) and run with 3 SNVs in one window:

```
total emitted: 26
hamming distribution: {0: 7, 1: 12, 2: 6, 3: 1}
distinct seqs: 8
distinct with >=2 SNVs: ['AAAAAGAATA', 'AACAAAAATA', 'AACAAGAAAA', 'AACAAGAATA']
```

All three pairwise doubles and the triple are produced. Not X1: the user authored singles, the tool
authored the combinations. User-facing, in the shipped Django app:
`iliasApp/forms.py:504` — `makeSnpCombinations = forms.BooleanField(label="Make Snp Combinations?",`;
wired at `iliasApp/views.py:342,362` — `parseSNPs.make_sequence_copies(..., CombinationsB=makeSnpCombinations)`;
documented at `iliasApp/templates/iliasApp/docs.html:42-43` — "The user can select to include or
exclude combinations of SNPs when designing MPRA experiments."

*Disconfirmation.* Three ways this could have failed and did not: (i) *is it only motifs?* — no; the
motif engines (`oligo.py:39-48` `itertools.product(range(...), repeat=len(motifsL))`; `part1.py:93-127`
4-motif spacing sweep) are excluded by decision (a), and the SNV power set stands independently;
(ii) *is it dead code?* — no, `views.py:362` is the live call and `:363` shows the *alternative*
(`naman_make_sequence_copies`) commented out, so the power-set function is the one in service;
(iii) *is it X1?* — no, the VCF holds singles. Caveat recorded but not scored here: the code is
Python 2 and will not execute on a modern interpreter (that belongs to `installable_today`).

### 3.4 mpradesign — **no**

Exported surface is two functions (`NAMESPACE`): `export(processVCF)`, `export(spread_and_fix_indels)`.
Oligo construction, `R/processVCFfast.R:392-393` (repeated at `:607-608`, `:828-829`):

```r
                     type = rep(c('ref', 'alt'), each = nper),
                     mid = ifelse(type == 'ref', snp$REF, snp$ALT),
```

Two alleles x `nper` barcodes per SNP — one designed event per output. The tool's own README states
the cardinality (`README.md:47`):

> 5 + .01 * Number of barcodes per allele * Number of SNPs in VCF * 2 (for ref/alt alleles)

Linear in the number of SNPs: the exact `sum_i |V_i|` signature my test uses to reject.

*Disconfirmation.* Searched both repos for any haplotype/multi-SNP construct:
`grep -rn -i "combn\|expand.grid\|haplotype\|permut\|combinat\|double" R/` -> only two hits, both the
English word "double-check" (`processVCFfast.R:287`, `import_freebarcodes.R:64`). Same grep over the
Shiny app `andrewGhazi/designMPRA` (`ui.R`, `server.R`, `Rmd/`, `scripts/`) -> **no output at all**.
The one multi-position writer, `nonrandomly_replace_N` (`:100-109`, loop over `bases_to_replace`), sits
under the header "Randomly correct aberrant digestion sites" (`:111-114`) — restriction-site repair,
exclusion X2. Searched and absent.

### 3.5 oligopoolcalc — **no** (CHANGED from partial)

There is no mutagenesis anywhere in the package:

```
$ grep -rn -i "mutagen\|mutate\|point mutation\|substitut" --include=*.py oligopool/
   oligopool/base/scry.py:584:            # Nope, Matches Overmutated
   oligopool/base/utils.py:3126:       desc - module stats dict (mutated in place)
```

Both hits are unrelated comments. The documented workflow has the user bring the variants
(`docs/docs.md:312`): "1. Start with your core sequences (variants, promoters, genes, etc.)".
The design modules (`barcode`, `primer`, `motif`, `spacer`, `pad`, `split`, `final`) add elements to a
core sequence; none varies a core at >=2 chosen sites.

The prior "partial" rests on `expand`, which is genuine multi-position cartesian expansion
(`oligopool/base/core_degenerate.py:858`):

```python
    return [''.join(combo) for combo in itertools.product(*choices_per_position)]
```

exposed top-level (`oligopool/__init__.py:22 'expand'`) and documented ("Expand IUPAC-degenerate
sequences into all concrete A/T/G/C sequences", `oligopool/expand.py:23`; "Expansion can be
exponential (e.g., 10 N's = 4^10 sequences)", `:46`). **Decision (b) excludes it**, and the
consistency argument is decisive rather than aesthetic: pydna ships the identical function
(`src/pydna/primer_screen.py:144-178`, `return ["".join(tup) for tup in product(*choices_per_pos)]`,
docstring example `expand_iupac_to_dna("ATNG") -> ['ATGG','ATAG','ATTG','ATCG']`). A threshold that
scores Oligopool "partial" for IUPAC expansion must score pydna "partial" too. It cannot, so neither
does. Oligopool's own framing agrees that `expand` is not a design operation
(`docs/docs.md:945`): "**When to use it**: Verifying that `compress` output covers all your original
variants."

*Disconfirmation.* What would have kept "partial": (i) a mutagenesis module — searched, absent (grep
above, plus `grep -n -i "mutagen\|saturation\|variant\|combinatorial\|double mutant" docs/docs.md`,
whose "variant" hits are all *user-supplied* variants, e.g. `:106` "# Start with your variants");
(ii) `motif`/`barcode` enumerating multi-element combinations per oligo — `docs/docs.md:483` describes
`motif` as "Inserts sequence motifs (per-variant or constant anchors) with constraint satisfaction",
one motif per column, and motif placement is excluded by (a) anyway; (iii) a per-oligo variable-region
combinatorics — a barcode is one contiguous designed element (X3). **This demotion is not a criticism
of Oligopool Calculator**: it is a library-assembly + NGS-readout tool that does not claim to design
mutants, and its real strengths land in `barcode_generation`, `degenerate_iupac_codons` and
`readout_analysis`.

### 3.6 mutationmaker — **yes**

PAS. `backend/mutation_maker/generate_oligos.py:68-72`:

```python
def cartesian_dictionary(*args, fold=mul) -> Dict:
    """ Cartesian product of an arbitrary number of dictionaries. """

    return {ks: reduce(fold, starmap(getitem, zip(args, ks)))
            for ks in product(*args)}
```

The arguments are **one dict per mutation site** — built in a `for site in mutations_sites:` loop and
appended at `:298-301` (`mutations_on_site_with_prob.append(aminos_with_probabilities)`;
`chosen_codons_on_sites.append(aminos_with_codons)`) — and the product is taken across sites
(`:303-305`). The stated algorithm, `generate_solution` docstring (`:238-247`):

> 6. Generate combinations with concetrations for different sites with the help of cartesian
> multiplication 7. For every combination replace aminos on mutation sites with selected previousely
> codons (degenerete or normal ones)

Each combination becomes ONE oligo carrying a codon at EVERY site (`:129-134`):

```python
    for combination in mutations_combinations_with_probabilitites:
        concentration = mutations_combinations_with_probabilitites[combination]
        for position, mutation in enumerate(combination):
            codon = chosen_codons_on_sites[position][mutation]
            dna = Codons.replace_codon(dna, mutations_sites[position], codon, start, goi_offset)
        oligos.append(PASOligo(sequence=dna, ratio=concentration))
```

**Behavioural verification** — `cartesian_dictionary` + `generate_oligos_from_combinations`
transcribed and run with two sites (codon 2: W/L/E; codon 4: A/G):

```
n combos: 6 = 3*2 (multiplicative over sites)
('W','A') ATGTGGCCCGCGTTT codons changed vs WT: 2
('W','G') ATGTGGCCCGGCTTT codons changed vs WT: 2
('L','A') ATGCTGCCCGCGTTT codons changed vs WT: 2
... (all six)
```

QCLM independently. `backend/mutation_maker/mutation.py:210-213,231-234`:

```python
class QCLMMutationSiteSequence:
    """
    Combines multiple mutation sites together, with utility function
    for generating options for AAs mutation combinations and mutation codons.
...
        # All n-tuple combinations of triplets at mutation sites which are close to each other.
        self.concrete_mutations = list(
            itertools.product(*[m.get_all_concrete_triplets(codon_usage_table, frequency_threshold)
                                for m in self.ordered_mutations]))
```

plus `:269` `return list(itertools.product(*[m.new_aminos for m in self.ordered_mutations]))`, and
`:227-228` `self.aminos_count = self.aminos_count * len(mutation.new_aminos)` — the tool computes the
multiplicative cardinality itself.

Reachability: `api/api.py:48,56` — `@hug.post('/qclm', versions=1)` and `@hug.post('/pas', versions=1)`;
`OligoGenerator` instantiated in the shipped output path at `backend/mutation_maker/pas_output.py:271`;
`QCLMMutationSiteSequence` consumed in `qclm.py:78,698,707-708`, `qclm_types.py:137`, `site_split.py:131`.
`README.md:287`: "All Mutation Maker workflows (SSSM, MSDM, PAS) are covered with unit tests of the
backend features."

*Disconfirmation.* Three failure modes checked. (i) *Is the product within one codon rather than
across sites?* No — `mutations_sites = Mutations.list_of_mutation_sites(mutations_list)` is
`list(set([item[0] for item in mutations_list]))` (`:190-193`), i.e. distinct residue positions, and
`replace_codon` indexes by `3 * (position - 1) + goi_offset - start` (`:165-168`). (ii) *Is it X1 —
did the user supply the combination?* No; the user supplies per-site point mutations (`E42W,E42L,E42K`
form, see `get_mutation_strings` docstring at `:271-275`) and the tool crosses them. (iii) *Is it
primers only, not sequences?* No; PAS emits `PASOligo(sequence=dna, ...)`, a fragment carrying all the
substitutions. Value stands.

### 3.7 codongenie — **no**

`codon_genie/codon_selector.py:71-79`:

```python
    def optimise_codons(self, amino_acids, organism_id):
        '''Optimises codon selection.'''
        req_amino_acids = set(amino_acids.upper())

        codons = [CODONS[amino_acid] for amino_acid in req_amino_acids]

        results = [self.__analyse(combo, organism_id, req_amino_acids)
                   for combo in itertools.product(*codons)]
```

Every `itertools.product` in the file (`:79`, `:97`, `:112`, `:159`) crosses alternatives **within one
3-base codon**: `CODONS` maps each amino acid to 3-element nucleotide-choice lists (`:20-36`, e.g.
`'L': [['C','T','ACGT'], ['T','T','AG']]`), and `__analyse` transposes to exactly three positions
(`:88-92`, `transpose[:2] + [_optimise_pos_3(transpose[2])]`). The deliverable is an ambiguous codon
string such as `NNK`, not a sequence: `ambig_codons = [''.join([NUCL_CODES[term] for term in cdn]) ...]`
(`:95-97`). One contiguous 3 bp site = one event (X3), and no reference sequence is involved at all.

*Disconfirmation.* A promotion would need a multi-codon or whole-sequence mode. Enumerated the entire
Python surface (`main.py`, `codon_genie/{seq_utils,codon_selector,codon_utils,ncbi_tax_utils,client}.py`);
`analyse_codon(self, ambig_codon, tax_id)` (`:82`) also takes a single codon. No function accepts a
sequence plus a site list. Searched and absent.

### 3.8 dnachisel — **partial**

The machinery is real and I want the authors to see how real it is.
`dnachisel/MutationSpace/MutationSpace.py:132,157-163`:

```python
    def all_variants(self, sequence):
        """Iterate through all sequence variants in this mutation space."""
...
        variants_slots = [
            [
                (choice_.segment, v.encode())
                for v in sort_variants_by_distance_to_current(choice_)
            ]
            for choice_ in self.multichoices
        ]
        for variants in itertools.product(*variants_slots):
```

`self.multichoices` are independent choice slots at **distinct positions**, all applied to the same
sequence before it is yielded (`:161-163`). There is also true multi-site sampling
(`:105-130`):

```python
    def pick_random_mutations(self, n_mutations, sequence):
        """Draw N random mutations."""
...
                for i in np.random.choice(
                    len(self.multichoices), n_mutations, replace=False
                )
...
    def apply_random_mutations(self, n_mutations, sequence):
        """Return a sequence with n random mutations applied."""
```

and the tool computes the multiplicative cardinality itself (`:95-103`, `space_size`, "Return the
number of possible mutations", `np.log(choices).sum()`). The space is derived from a *specification*
relative to a reference — `MutationSpace.from_optimization_problem(problem, new_constraints=None)`
(`:165-172`) — and `MutationSpace` is in the documented API reference
(`docs/ref/core_classes.rst:45-48`, `.. automodule:: dnachisel.MutationSpace`).

**Why `partial` and not `yes`.** These are solver primitives, not design operations. The only callers
of `all_variants` inside the package are the two search mixins —
`DnaOptimizationProblem/mixins/ConstraintsSolverMixin.py:63` and
`ObjectivesMaximizerMixin.py:45` — which consume the generator to find ONE constraint-satisfying
sequence; `problem.sequence` is what the user gets. `MutationSpace` is not re-exported from
`dnachisel/__init__.py` (`grep -n "MutationSpace\|MutationChoice" dnachisel/__init__.py` -> no output),
so a user must reach into the solver object, and no builtin specification, CLI command
(`dnachisel/cli.py`) or documented example emits a multi-mutant library. Restricted axis => PARTIAL.

*Disconfirmation / the hinge, stated for the authors.* A referee reading the row question's "enumerate
or sample" mechanically has a strong case for `yes` here: `problem.mutation_space.all_variants(problem.sequence)`
is a documented public generator that yields pairwise-and-higher mutation combinations, and
`apply_random_mutations(k, seq)` is functionally PoolParty's random k-mutation mode. My demotion rests
entirely on *purpose and delivered output* (one optimised sequence), not on absence of machinery. I
searched for the thing that would force `yes` — a user-facing library-emitting entry point — via
`grep -rn "all_variants" dnachisel/` (4 hits, all internal or IUPAC-pattern expansion at
`SequencePattern/DnaNotationPattern.py:35`) and `grep -rn "all_variants\|MutationSpace" docs/`
(2 hits, both the API-reference stub). It is not there. If the authors prefer the mechanical reading,
this is the cell to flip, and it should be flipped together with Ledidi and tangermeme.

### 3.9 ledidi — **partial**

Ledidi genuinely produces multi-edit sequences and I credit that fully.
`ledidi/ledidi.py:80-89` (function docstring):

> Ledidi is a method for designing compact sets of edits to categorical sequences, such as DNA, to
> make them exhibit desired characteristics as predicted by an oracle model.

`fit_transform` returns a batch of them (`:532-536`):

> y: torch.Tensor, shape=(batch_size, n_channels, length)
>     A tensor containing a batch of one-hot encoded sequences which
>     may contain one or more edits compared to the sequence that was
>     passed in.

Sampling is real (Gumbel-softmax straight-through, `:29 _gumbel_softmax_hard`, described at `:91-97`),
`n_samples` draws arbitrarily many (`:154-157`), and the edit count is the object of the input loss
("which corresponds to the number of positions that have been edited", `:299-301`), with masking via
`input_mask` (`:341-345`, `:551-553`).

**Why `partial`.** The combinatorial space is not declarable and is never enumerated or uniformly
sampled: there is no k-mutation parameter — only `l`, "The mixing weight parameter between the input
loss and the output loss" (`:316-320`) — and the sampling distribution is fitted to a model target,
so the user cannot request "the pairwise set". Per-decision (c) that is a restricted axis. The
capability itself is scored by the row set's `ml_model_in_loop`.

*Disconfirmation / the hinge.* Promotion evidence I looked for and did **not** find: a parameter fixing
the number of edits, or an enumeration mode. Checked the full `Ledidi.__init__` signature
(`:377-381`: `input_loss, output_loss, tau, l, batch_size, max_iter, early_stopping_iter, report_iter,
lr, input_mask, initial_weights, eps` — no `n_edits`), the `ledidi()` wrapper signature (`:78`), and
`ledidi/pruning.py:25 greedy_pruning(model, X, X_hat, threshold=1, ...)`, which *reduces* an existing
multi-edit set rather than enumerating one. Conversely, the argument for `yes` is straightforward and
recorded: the row says "or sample", and Ledidi samples multi-position edit sets by design and returns
a batch of distinct ones. Flip together with DnaChisel and tangermeme if the authors take the
mechanical reading.

### 3.10 tangermeme — **partial**

Multi-mutation accumulation is real, in `tangermeme/design/`. `greedy_substitution.py:44-48,79-81`:

> Greedily add motifs to achieve a desired goal. ... This process will continue until either the
> maximum number of iterations is reached (at which point, `max_iter` motifs will have been inserted
> into the sequence) or the loss falls below `tol`.
>
> motifs: list of strings or None ... If None, use the provided alphabet as the motifs to
> only change one character at a time.

With `motifs=None` this is exactly iterative point mutagenesis: `max_iter` single-base substitutions
accumulate in ONE output sequence, with the tool choosing each position/base
(`greedy_substitution.py:186`, `while iteration < max_iter:`). `beam_substitution.py:47,61-64` keeps
`beam_size` distinct multi-edit trajectories ("the global top-`beam_size` are kept as the next beam").

**Why `partial`.** (i) It is an objective-driven search — decision (c); the user cannot specify the
combinatorial space, and the deliverable is 1 (or `beam_size`) sequences. (ii) The parts that *are*
enumerative are motif placement, excluded by decision (a): `ersatz.multisubstitute(X, motifs, spacing,
...)` — "Substitute a set of motifs into sequences with provided spacings" (`ersatz.py:206`), positions
fixed by the user's `spacing`, one output per call; `space.py:34` sweeps a user-supplied spacing grid
("Each row in this tensor is a different combination of spacings between motifs", `:63-65`).
(iii) The mutational scan is explicitly single-event: `saturation_mutagenesis.py:214` —
`# Build only the true single-base edits. One substitution per position`, and `:88`
"For each single-character substitution, the change in ...".

*Disconfirmation.* I searched specifically for a pairwise/double-mutant mode. `grep -n "def
\|pairwise\|double" tangermeme/saturation_mutagenesis.py` -> single-base only. The package's
`product.py` (`itertools.product` at `:145`, `:300`) is **not** mutation combinatorics: `apply_pairwise`
"Apply a function on the cartesian product between X and args" (`:40`), where args are model inputs
(the docstring's own example is "cell state and read depth", `:52-58`). Checked
`ablate.py`, `marginalize.py`, `variant_effect.py`, `design/screen.py` (`:37` "Screen randomly
generated sequences and choose the best one" — i.i.d. random sequences, not mutation combinations).
No enumerated multi-mutation mode exists; the multi-edit output is search-derived. Same hinge as
DnaChisel/Ledidi.

### 3.11 biopython — **no**

Current package, sparse checkout of `Bio/` + `Doc/` at HEAD:

```
$ grep -rn -i "def mutate\|all_variants\|saturation" --include=*.py Bio/
   (no output)
$ grep -rln -i "mutagen" Bio/
   Bio/Entrez/DTDs/NCBI_Seqfeat.mod.dtd   Bio/Entrez/DTDs/NCBI_BioSource.mod.dtd
   Bio/Entrez/DTDs/NCBI_Sequence.mod.dtd  Bio/Entrez/DTDs/Docsum_3_4.mod.dtd
   Bio/Entrez/DTDs/Docsum_3_3.mod.dtd     Bio/SeqUtils/MeltingTemp.py
   Bio/SwissProt/__init__.py
```

The only hits are database schemas and the SwissProt `MUTAGEN` feature key — parsers, not designers.
Biopython offers `MutableSeq` (a mutable container the *user* edits, exclusion X5); there is no
operation that enumerates or samples mutation combinations.

*Disconfirmation.* Because decision (b) turns on IUPAC expanders, I checked specifically whether
Biopython ships one — it does not: `grep -rn "ambiguous_dna_values" --include=*.py Bio/` yields only
`Bio/Data/IUPACData.py:82` (the table), `Bio/Data/CodonTable.py` (ambiguous stop-codon tables),
`Bio/SeqUtils/__init__.py:285` (molecular-weight ranges) and `Bio/Nexus/Nexus.py:776`; and
`grep -rn "itertools.product" Bio/Seq.py Bio/SeqUtils/ Bio/motifs/` -> no output. So Biopython is `no`
under the strict *and* the mechanical reading. Searched and absent.

### 3.12 pydna — **no**

All of pydna's combinatorics is fragment assembly and restriction/recombination pairing, not
mutation: `assembly.py:154-155` (`itertools.combinations('ABCD', 2)` over fragments), `:359`, `:429`;
`assembly2.py:275,1282-1287,1510,1918,1982-1983`; `fusionpcr.py:49`; `recombinase.py:302`;
`snapgene_history_parser.py:249`. Line `assembly2.py:1500` states the object: "for example two
overlaps between 1 and 2, and single overlap between 2 and 3 should return 3 assemblies" — assemblies
of fragments, not co-occurring mutations of a reference. Site-directed mutagenesis appears only as
*replay* of a foreign file format (`snapgene_history_parser.py:345`, `elif node.operation ==
"primerDirectedMutagenesis":`) and as an external third-party notebook link
(`docs/example_gallery.rst:31`, a teemi notebook).

pydna does ship a full IUPAC expander, `src/pydna/primer_screen.py:144-178`:

```python
def expand_iupac_to_dna(seq: str) -> list[str]:
    """
    Expand an extended IUPAC DNA string to unambiguous IUPAC nucleotide alphabet.
...
    return ["".join(tup) for tup in product(*choices_per_pos)]
```

This is excluded by decision (b) — and it is the reason (b) exists, since it is functionally identical
to Oligopool's `expand`. Here it is a primer-matching helper (module `primer_screen`, feeding
`make_automaton`/`forward_primers` at `:181,290`; module docstring context at `:204` "The primers can
contain ambiguous bases from the extended IUPAC DNA alphabet").

*Disconfirmation.* Searched for a designed-mutant path: `grep -rn -i "def mutate\|mutagen" --include=*.py
src/` -> only the SnapGene replay hits above; and the full module listing of `src/pydna/` (41 modules,
including `design.py`, `primer.py`, `crispr.py`, `codon.py`) contains no mutagenesis designer.
Searched and absent.

### 3.13 seqpro — **no**

The complete public API (`python/seqpro/__init__.py:3-9`):

```python
from . import alphabets, bed, gtf, rag, transforms
from ._analyzers import gc_content, length, nucleotide_content
from ._encoders import decode_ohe, decode_tokens, ohe, pad_seqs, tokenize
from ._modifiers import bin_coverage, jitter, k_shuffle, random_seqs, reverse_complement
from ._types import PathLike
from ._utils import NestedStr, SeqType, StrSeqType, cast_seqs
from .alphabets import AA, DNA, RNA, AminoAlphabet, NucleotideAlphabet
```

Encoders, analyzers, ragged arrays, shuffles and `random_seqs` (i.i.d. sequences, not mutations of a
reference). No mutation operation of any order. This is exclusion X5 territory: a sequence-array
toolkit.

*Disconfirmation.* `grep -rn -i "mutat\|combinatorial\|itertools.product" --include=*.py --include=*.rs .`
over the whole repo (Rust `src/`, `crates/`, `python/`, `tests/`) returns only in-place-buffer-mutation
comments (`tests/test_ragged_rc.py:137,150`, `python/seqpro/rag/_ops.py:52,245`, `scratch_bench_rc.py:117`).
Searched and absent, in both the Python and Rust halves.

---

## 4. Row verdict

**keep** — the row scored cleanly and consistently across all 13 tools once decisions (a), (b) and (c)
were applied, and it discriminates (4 yes / 3 partial / 6 no).

**Recommended wording tightening**, because a referee should not have to reconstruct (a)–(c):

> `higher_order_combinatorial` — Pairwise and higher-order combinations of **mutations of a reference
> sequence** (>=2 independently-positioned substitution/indel/codon events) enumerated or sampled by
> the tool in one output sequence. Combinatorial motif placement is scored in
> `combinatorial_motif_place`; expansion of a user-authored IUPAC-degenerate template is scored in
> `degenerate_iupac_codons`; multi-edit sequences produced by an objective-driven optimiser without a
> declarable combinatorial space score *partial*.

Adding that sentence converts three judgement calls into rule applications and makes the row
re-runnable by a third party.

## 5. Files and commands used

- PoolParty source (read-only): `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/src/poolparty/base_ops/mutagenize.py`, `orf_ops/mutagenize_orf.py`, `multiscan_ops/deletion_multiscan.py`, `base_ops/recombine.py`
- PoolParty behavioural runs: `PYTHONDONTWRITEBYTECODE=1 /mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/.venv/bin/python` (import only; no writes)
- Competitor sources cloned read-only into the session scratchpad
  (`.../scratchpad/repos/{valiant,MPRAnator,Mutation_Maker,mpradesigntools,designMPRA,oligopool,CodonGenie,DnaChisel,ledidi,tangermeme,biopython,pydna,seqpro}`); nothing installed.
- Nothing outside `revision/tool_survey/verified/` was written; `poolparty-statecounter` was read only.
