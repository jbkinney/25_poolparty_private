# Row audit: `sampled_random_mutagenesis`

Single-auditor pass across all 13 tools, one threshold, primary sources only.
Date: 2026-08-10.

---

## 1. Operational test (fixed BEFORE any tool was opened)

> **Is there a documented, user-settable parameter — a mutation rate, a count of
> mutations/variants, or an RNG seed — such that the tool DRAWS SUBSTITUTED CHARACTERS INTO
> POSITIONS OF A USER-SUPPLIED TARGET SEQUENCE and emits the resulting randomized variants,
> so that two runs with different seeds yield different variant sets (i.e. a proper random
> subset of the variant space rather than the deterministic full enumeration)?**

Scoring, applied identically to all 13:

| value | condition |
|---|---|
| **yes** | Such a parameter exists AND the drawn-edit variants are the deliverable of a documented user-facing entry point (function / CLI flag / web-form field). |
| **partial** | The tool's RNG *does* substitute drawn characters into positions of a user-supplied target and those sequences reach the output, but the draw is not governed as a mutagenesis sample: it redraws a whole designated window WT-independently, or it only redraws **synonymous** codons (a non-mutational component), or it is a move operator inside an optimizer whose deliverable is a single optimized sequence. |
| **no** | No such randomness anywhere: exhaustive enumeration only, OR the RNG only produces sequence **de novo** (random DNA, barcodes, spacers, primers, pads), OR only **permutes** existing characters (k-let shuffles), OR only affects non-sequence choices (tree topology, bootstrap columns, plot jitter). |

### Two explicit discriminators a referee can re-run

1. **De-novo generation vs. target editing.** `pydna.utils.randomDNA(length)` and
   `seqpro.random_seqs(shape, alphabet, seed)` are documented random *sequence generators* with
   seeds. Neither takes a target sequence. If "has an RNG that emits DNA" earned partial credit,
   pydna and SeqPro would score partial and the row would carry no information. Therefore
   de-novo generation ⇒ **no**. Corollary: random **barcode/spacer/primer/pad** generation, even
   with a documented `random_seed`, is de-novo generation of *auxiliary* elements added around
   user-supplied variants ⇒ **no**.
2. **Permutation vs. substitution.** A k-let-preserving shuffle changes position but introduces
   no new character and produces the same composition; it is the canonical negative control, not
   a mutant ⇒ **no**. Drawing a *new* character into a position of the target ⇒ at least
   **partial**.

I note where this test is close to flipping a cell, because two cells (DnaChisel, tangermeme)
are genuinely arguable and the paper should be prepared to defend them.

---

## 2. Result

| tool | prior | **audited** | one-line basis |
|---|---|---|---|
| poolparty | yes | **yes** | `mutagenize(mutation_rate=…, mode='random')`; `rng.binomial(num_mutable, mutation_rate)`; seeded via `generate_library(seed=)`; verified by execution |
| valiant | no | **no** | zero `random`/`seed`/`sample` tokens in 85 `.py` files; 10 CLI options, none stochastic; 17 JSON config `params` keys, none stochastic |
| mpranator | yes | **yes** | `Part3Form.numToMutate` → `part3.mutateString(seq, noOfPos)`; `rd.choice` on positions and on non-WT bases; UI text "Creates a specified number of random mutations on a sequence" |
| mpradesign | no | **partial → CHANGED** | `randomly_fix()` draws (position, non-WT allele) pairs with `dplyr::sample_n()` and writes them into `constrseq` — but as restriction-site repair, not variant sampling |
| oligopoolcalc | no | **no** | `random_seed` documented in `barcode`/`primer`/`motif`/`spacer`/`split`/`pad` only; variants are user-supplied input; "no invented sequences" guarantee |
| mutationmaker | partial | **partial** | RNG only picks **synonymous** codons by usage frequency, plus a 500-random-combination heuristic; amino-acid variant space is user-specified |
| codongenie | partial | **partial** | `CodonOptimiser.mutate(protein_seq, dna_seq, mutation_rate)` is a real per-codon rate, but each draw is a **synonymous** codon for the same amino acid |
| dnachisel | partial | **partial** (closest to yes) | `MutationSpace.apply_random_mutations(n_mutations, sequence)` draws genuine non-WT substitutions, and `randomization_threshold` documents the exhaustive-vs-random switch — but they are hill-climber moves; `optimize()` returns one sequence |
| ledidi | partial | **partial** | `n_samples` + `random_state` draw N WT-relative edited sequences, but from the **fitted** design distribution after optimization; no rate/count of edits |
| tangermeme | partial | **partial** | `ersatz.randomize(X, start, end, probs, n, random_state)` redraws a whole window WT-independently (control/ablation primitive); `saturation_mutagenesis` is exhaustive with no sampling parameter |
| biopython | no | **no** | RNG confined to `Bio/Phylo`, `Bio/Nexus`, `Bio/Graphics`; no sequence mutator or random-sequence function exists |
| pydna | no | **no** | `randomDNA/randomRNA/randomORF/randomprot` are de-novo generators; no mutagenesis function in 41 modules |
| seqpro | no | **no** | `random_seqs` (de novo), `k_shuffle` (permutation), `jitter` (coordinate shift), `transforms.Random` (probability of applying a transform); no substitution into a target |

Spread: **2 yes / 6 partial / 5 no.** One change from the prior pass (mpradesign no → partial).

---

## 3. Per-tool evidence

### 3.1 poolparty — **yes**

`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/src/poolparty/base_ops/mutagenize.py`

Signature (L18–32):

```python
def mutagenize(
    pool: Union[Pool, str],
    region: RegionType = None,
    num_mutations: Optional[Integral] = None,
    mutation_rate: Optional[Real] = None,
    ...
    mode: ModeType = "random",
    num_states: Optional[Integral] = None,
```

Docstring (L45–46, L55–60):

> `mutation_rate : Optional[Real], default=None`
> `    Probability of mutation at each position (mutually exclusive with num_mutations).`
> `mode : ModeType, default='random'`
> `    Selection mode: 'random' or 'sequential'. Sequential only available with num_mutations.`
> `num_states : Optional[int], default=None`
> `    … In random mode, if None defaults to 1 (pure random sampling).`

The draw itself, `_random_mutation`, L455–477:

```python
if self.num_mutations is not None:
    num_mut = min(self.num_mutations, num_mutable)
else:
    num_mut = rng.binomial(num_mutable, self.mutation_rate)
...
chosen_indices = rng.choice(num_mutable, size=num_mut, replace=False)
...
mut_indices = (rng.random(num_mut) * counts).astype(np.intp)
mut_bytes = mutation_options_arr[chosen_indices, mut_indices]
```

`mutation_options` excludes the WT base (L426 `valid_muts = dna_utils.MUTATIONS_DICT[wt]`;
L421–423 `... if b != wt_upper`), so every draw is a genuine mutation. Sequential (exhaustive)
mode is a separate branch: `_build_caches` enumerates `combinations(range(num_positions),
num_mutations)` × `(alpha_size-1)**num_mutations` (L351–363). L551–553 requires an RNG for random
mode: `raise RuntimeError(f"{self.mode} mode requires RNG - use Party.generate(seed=...)")`.

Docs, `poolparty/docs/operations/mutagenize.rst` L128–137:

> `Per-base mutation rate (mutation_rate=0.1)`
> `mutation_rate applies an independent per-position probability; the number`
> `    mutants = wt.mutagenize(mutation_rate=0.1, mode="random", style="red")`

Seed: `poolparty/src/poolparty/generate_library.py` L20/L35 `seed: Optional[Integral] = None` /
"Random seed for reproducibility"; L283 `seed_seq = np.random.SeedSequence([pool._master_seed, op.id, state_val])`.

**Behavioural verification** (read-only venv, `PYTHONDONTWRITEBYTECODE=1`,
`/mnt/c/.../poolparty-statecounter/.venv/bin/python`), WT = 32 nt:

```
seed1: ['ACGTACGTACGAACGTACGGACGTACGTACGT', 'ACGTACGTACGTACGTGTATAGGTACGTAGGT', ...]
seed2: ['ATGTTCATCCGTACGTACGTACGTCCGTACGT', 'ACGCACGCACGTACGTACGTACGTACGTACGT', ...]
same set: False
HD seed1 (rate=0.1, L=32, E[HD]=3.2): [2, 5, 4, 1, 2, 4, 7, 4]
reproducible with same seed: True
num_mutations=3 HD: [3, 3, 3, 3, 3, 3, 3, 3]
mode="sequential" first 3: ['CAATACGT…', 'CACTACGT…', 'CATTACGT…']   (exhaustive)
```

Binomial spread around the expected 3.2 confirms a true per-position rate, not a fixed count.

**Disconfirmation.** A "partial" would have followed if the rate had turned out to govern only
barcode/degenerate filling, if `mode='random'` had been an alias for enumeration order, or if the
mutation options had included the WT base (making it a randomization rather than a mutagenesis).
Checked: `_random_mutation` L455–477, `_get_position_mutations` L384–432 (WT excluded),
`_build_caches` L332–382 (the separate exhaustive path), and the execution above. None holds.
PoolParty was scored with exactly the same test as the competitors, including the requirement that
the sample be the deliverable; `generate_library` returns the variant table, so it is.

### 3.2 valiant — **no** (anchor, re-verified independently)

Clone: `github.com/cancerit/VaLiAnT`, branch `develop` (`git branch --show-current` → `develop`),
85 `.py` files under `src/`.

`src/valiant/mutator_type.py` L22–33:

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

Three independent surfaces checked and all negative:

- `grep -rniE "\brandom|\bseed\b|shuffle|choice\(|sample" src/ --include=*.py` → **one** hit,
  `src/valiant/common_cli.py:53: type=click.Choice(list(logging._nameToLevel.keys()), …)` — the
  log-level CLI type. No RNG anywhere in the package.
- Full CLI option inventory (all `click.option` long flags in `src/valiant/*.py`):
  `--adaptor-3 --adaptor-5 --annot --bg --bg-mask --codon-table --config --gff --pam --vcf`.
  No `--seed`, `--rate`, `--sample`, `--n`.
- JSON config `params` keys (`examples/sge/nxrl/chr17_31226400_31226678_plus_sgRNA_1146241047_valiant_config.json`):
  `GFFFilePath, PAMProtectionVCFFilePath, adaptor3, adaptor5, assembly, backgroundVCFFilePath,
  codonTableFilePath, customVCFManifestFilePath, forceBackgroundFrameShifting,
  forceBackgroundNonSynonymous, maxOligoLength, minOligoLength, oligoInfoFilePath, outputDirPath,
  refFASTAFilePath, reverseComplementOnMinusStrand, species`. None stochastic.
- `grep -rniE "\brate\b|probab|stochast" src/ --include=*.py` → zero hits.

**Disconfirmation.** A rate/seed/count in the CLI, in the config `params`, or a stochastic
`MutatorType` member would have moved this to yes/partial. All three enumerated exhaustively above.

### 3.3 mpranator — **yes**

`part3.py` L13–48 (`mutateString`):

```python
def mutateString(sequenceS, noOfPos, mapperDict=myf.mapperDict):
    """
    Args: This function mutates a string
        …
        noOfPos will choose the number of mutations of helper function
    """
    ...
    for i in range(noOfPos):
        posChosen = rd.choice(LengthList)
        while sequenceS[posChosen] == "N" or posChosen in positionsChosenL:
            posChosen = rd.choice(LengthList)
        positionsChosenL.append(posChosen)
    ...
    for i in positionsChosenL:
        therandomchoice = rd.choice(mapperDict[stringL[i]])
        stringL[i] = therandomchoice
```

`myfunctions.py` L7: `mapperDict = {"A": 'TGC', "G": 'ACT', "C": 'TGA', "T": 'GCA', …}` — the WT
base is excluded from every option list, so each draw is a true substitution.

User-facing parameter, `iliasApp/forms.py` L341:

```python
numToMutate = PositiveIntegerField(label="Number of nucleotides to mutate", initial=0)
```

Tool's own catalogue text, `iliasApp/views.py` L142:

```python
["Transmutation", urlresolvers.reverse("iliasApp:ViewPart3"), "Creates a specified number of random mutations on a sequence"],
```

Wired at `iliasApp/views.py` L204 and emitted as the download
(`context['fileName'] = "Transmutation_results.txt"`, L224).

**Disconfirmation.** This would have been partial if `mapperDict` had included the WT base
(randomization, not mutagenesis), if `numToMutate` had turned out to drive only the barcode
mutator, or if `mutateString` were unreachable from the app. Checked `myfunctions.py:7`,
`forms.py:340–352`, `views.py:150–229`. The barcode RNG is separate and lives in `part1.py`
L48/L75–76/L84–85 (`barcode.append(rd.choice(l1))`), which on its own would have been "no" under
discriminator 1; part3 is what earns the yes.

*Recorded limitation (not disqualifying):* one random mutant per input sequence per run, and no
seed parameter — reproducibility is not offered. The rubric requires "a documented rate, count, or
RNG parameter"; `numToMutate` is a documented count.

### 3.4 mpradesign — **partial (CHANGED from `no`)**

`mpradesigntools/R/processVCFfast.R`. Documented function header L112–115:

```r
#' Randomly correct aberrant digestion sites
#'
#' For a SNP with aberrant digestion sites in the context, randomly change bases
#' in the site across barcodes
```

The draw, L149–159 inside `randomly_fix()`:

```r
altered_patterns = tibble(pos_to_change = str_locate_all(aberrant_pattern, '[^N]')[[1]][,1],
                              possible_alleles = map(pos_to_change, ~dplyr::setdiff(c('A','C','G','T'),
                                                                                    substr(aberrant_pattern, .x, .x)))) %>%
  tidyr::unnest_legacy() %>%
  {dplyr::sample_n(., nrow(res_df) / 2, replace = ((nrow(res_df) / 2) > nrow(.)))} %>%
  mutate(altered_pattern = map2_chr(pos_to_change, possible_alleles, change_pattern, ab_pattern = aberrant_pattern))
```

and the write-back into the emitted construct, L167–169:

```r
constrseq_fixed = map2_chr(constrseq, fixed_pattern, reassign_pattern, aberrant_site_loc = aberrant_site)
```

`dplyr::setdiff(c('A','C','G','T'), substr(...))` excludes the WT base, so these are genuine
single-nucleotide substitutions, drawn at random positions, written into the user's oligo, and
carried to the output table. Also `randomly_change_pattern` L72–81 (`sample(1:nchar(dig_pattern),
size = 1)`; `new_allele = sample(allele_options, size = 1)`).

Why partial and not yes: the count of substitutions is not a user parameter — it is
`nrow(res_df)/2`, i.e. one distinct repair per barcode, and its purpose is to destroy an unwanted
restriction site, not to sample a variant library. The variant space itself is exactly the input
VCF (exhaustive over the user's SNPs). Why partial and not no: the RNG demonstrably substitutes
drawn non-WT bases into the target and those sequences are shipped in the output; that is more than
"random barcode assignment" (which is also present, L394/L519/L609 `sample(snp$bcPools %>% unlist, …)`
and would alone have been "no").

This code is live in both distributions: `designMPRA/server.R:17` is
`source('scripts/processVCFfast.R')`.

**Disconfirmation.** A `no` would have been correct if the only `sample()` calls had been barcode
assignment and the `runif`-based 12-mer shuffling at L1213
(`shuffled_mers = mers[sort(runif(length(mers)), index.return = TRUE)$ix]`). I grepped every
`sample(|runif|rbinom|set.seed|rnorm|random` occurrence in `R/processVCFfast.R` (1459 lines) and in
`designMPRA/server.R` + `ui.R`; `randomly_fix`/`randomly_change_pattern` are the ones that write
drawn bases into the construct. A `yes` would have required a user-settable rate/count of
mutations or an option to subsample the SNP set — no such argument exists in `randomly_fix`'s
signature (`snp, res_df, dig_patterns, dig_site_locations`) or in the Shiny inputs.

### 3.5 oligopoolcalc — **no**

`docs/api.md` documents `random_seed` for exactly six modules, mapped by section heading:

```
### barcode -> random_seed at api.md:89     ### primer -> :171     ### motif  -> :254
### spacer  -> :329                        ### split  -> :654     ### pad    -> :722, :791
```

each described identically (`api.md:114`, and at each other site):

> `- random_seed (int / None, default=None): RNG seed for reproducibility`

Source confirms these are local RNGs for element synthesis, e.g. `oligopool/barcode.py:33`
`random_seed:int|None=None`, `oligopool/motif.py:293` `rng = np.random.default_rng(random_seed)`.

Oligopool does not create variants at all — they are input (`api.md:824` `--input-data variants.csv`;
`api.md:1746` header `ID,Variant`), and `api.md:816` states the guarantee:

> `- Core guarantee (lossless): degeneracy(prefix) <= count(compatible variants) - no invented sequences`

`grep -niE "mutagen|random mutat|mutation rate|\bmutate\b" docs/docs.md` (98.9 KB) → three hits, all
describing *input* libraries to be compressed (`docs.md:895` "saturation mutagenesis libraries often
compress well"; `docs.md:2331` "Saturation Mutagenesis Compression"; `docs.md:2485` "Cost-efficient
mutagenesis | `compress` → order `synthesis_df`"). Module inventory (`oligopool/`): acount,
background, barcode, compress, expand, final, index, inspect, join, lenstat, merge, motif, pack,
pad, primer, revcomp, spacer, split, verify, xcount — no mutagenesis module.

**Disconfirmation.** A `partial` would have required a `random_seed` that redraws bases inside the
user's variant sequence. I checked all six documented `random_seed` modules: every one designs or
completes an element that is *added around* the variant (barcode, primer, spacer, motif insert, pad),
and the lossless guarantee above forbids altering variant bases. A `yes` would have required a
mutagenesis entry point; the 22-module API signature list in `docs/api.md` has none.

### 3.6 mutationmaker — **partial**

Every RNG call site in the backend
(`grep -rniE "\brandom\b|random\.|np\.random|shuffle|\bseed\b|sample\(" backend/mutation_maker/*.py backend/*.py api/*.py`):

- `backend/mutation_maker/generate_oligos.py:156–162` —
  `""" Returns random codon for a given amino acid based on codon's usage frequency. """`
  `return np.random.choice(list(codons_with_frequency.keys()), 1, replace=True, p=renormalized_probabilities)[0]`
- `backend/mutation_maker/reverse_translation.py:62` — identical draw.
- `backend/mutation_maker/pas_degeneracy_recursion.py:138–144` — identical draw; and L40
  `1. 500 random combinations of codons are created`, L74 `… (500 random solutions)`.
- `backend/mutation_maker/qclm.py:833, 851` —
  `is to just do a random selection of a codon for each amino. This way we won't` /
  `return random.choice(candidate_codons)`
- `backend/mutation_maker/degenerate_codon.py:347–350` — commented-out `random.choice`.
- `backend/ssm_algorithm_comparison.py:92, 102` — `np.random.choice(["A","C","T","G"], 7000, True)`
  and `mutation_sites = np.random.choice(...)`. **This is a benchmark harness, not a feature**:
  `def generate_inputs()` feeds `def main(): for i in tqdm(range(50)): …` for timing
  (`stats_dataframe_for_solution(solution, end-start)`), and it is not imported by `api/api.py`,
  `api/server.py`, or `api/server_fastapi.py`.

So the randomness redraws **synonymous** codons in the user's gene of interest and drives a
500-sample heuristic. The amino-acid-level variant space is exactly what the user submitted
(`E{site}X` style mutation strings). That is the rubric's "sampling only of a non-mutational
component" plus "stochastic optimization" ⇒ partial.

**Disconfirmation.** A `yes` would have needed a rate/count/seed in the documented input surface.
`grep -rniE "random|seed|rate" README.md docs/*.md` → the only hits are `README.md:54`
("interactive shell"), `:143` ("dependencies … requirements.txt"), `:194` ("separate terminals"),
`:296` (`PYTHONHASHSEED=0` for tests) — no user-facing stochastic parameter. `docs/` contains only
`plans/`. A `no` would have followed if the codon draw did not reach the emitted oligos; it does
(`generate_oligos.py`).

### 3.7 codongenie — **partial**

`codon_genie/codon_utils.py` L119–125:

```python
def mutate(self, protein_seq, dna_seq, mutation_rate):
    '''Mutate a protein-encoding DNA sequence according to a
    supplied mutation rate.'''
    return ''.join([self.get_random_codon(amino_acid)
                    if random.random() < mutation_rate
                    else dna_seq[3 * i:3 * (i + 1)]
                    for i, amino_acid in enumerate(protein_seq)])
```

and L135–137:

```python
def get_random_codon(self, amino_acid, excl_codons=None):
    '''Returns a random codon for a given amino acid,
    based on codon probability from the codon usage table.'''
```

This is a genuine per-codon Bernoulli rate applied to a user-supplied `dna_seq` — the strongest
"rate" wording in the row outside PoolParty — but the replacement is a codon for **the same amino
acid**, i.e. a synonymous re-encoding. No protein-level variant is created. The rubric's
"sampling only of a non-mutational component" ⇒ partial.

CodonGenie's actual product is exhaustive: `codon_genie/codon_selector.py` uses
`import itertools` over the hard-coded `CODONS` / `NUCL_CODES` IUPAC tables (L16–54) to enumerate
and score degenerate codons covering a requested amino-acid set. No RNG there.

`mutate()` is also **not reachable from the tool's documented interface**: `main.py` exposes only
`@app.route('/')`, `/organisms/`, `/organisms/<term>`, `/codons` (L40–62) and imports only
`get_codon_usage_organisms` from `codon_utils` (L22); `mutate` appears nowhere else in the repo
(`grep -rn "mutate\|mutation_rate" main.py codon_genie/*.py static/codonGenie/*.js`). It is,
however, a public method of a class in a packaged distribution (`setup.py:17 name='CodonGenie'`,
`:25 packages=setuptools.find_packages()`, `pypi.sh` = `twine upload dist/*`), so I did not use
unreachability to demote it — the synonymous-only nature is the reason for partial, and that
holds either way.

**Disconfirmation.** A `yes` would have followed if `get_random_codon` had drawn a codon for a
*different* amino acid, or if `/codons` accepted a rate/count/seed. Checked
`codon_utils.py:127–157` (all three getters are per-amino-acid) and `main.py:40–91`.

### 3.8 dnachisel — **partial** *(highest-risk cell in this row; see §4)*

DnaChisel has a first-class, explicitly named mutation space and both an exhaustive and a random
traversal of it.

`dnachisel/MutationSpace/MutationSpace.py` L106–130:

```python
def pick_random_mutations(self, n_mutations, sequence):
    """Draw N random mutations."""
    n_mutations = min(len(self.multichoices), n_mutations)
    ...
    return [ (choice_.segment, choice_.random_variant(sequence=sequence))
             for choice_ in [ self.multichoices[i]
                              for i in np.random.choice(len(self.multichoices), n_mutations, replace=False) ] ]

def apply_random_mutations(self, n_mutations, sequence):
    """Return a sequence with n random mutations applied."""
```

`dnachisel/MutationSpace/MutationChoice.py` L38–46:

```python
def random_variant(self, sequence):
    """Return one of the variants, randomly."""
    subsequence = sequence[self.start : self.end]
    variants = [v for v in self.variants if v != subsequence]
    variants = sorted(variants)
    return variants[np.random.randint(len(variants))]
```

The WT subsequence is excluded, so these are genuine mutations — structurally the same draw as
PoolParty's `num_mutations` mode. And the documented class attributes
(`dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py` L76–90) name the very dichotomy this
row asks about:

> `randomization_threshold`
> `  The algorithm will use an exhaustive search when the size of the mutation space (=the number of possible variants) is above this threshold, and a (guided) random search when it is above.`
> `max_random_iters`
> `  When using a random search, stop after this many iterations`
> `mutations_per_iteration`
> `  When using a random search, produce this many sequence mutations each iteration.`

(defaults L109–111: `randomization_threshold = 10000`, `max_random_iters = 1000`,
`mutations_per_iteration = 2`). `dnachisel.MutationSpace` is in the published API reference:
`docs/ref/core_classes.rst` L45–48 `MutationSpace` / `.. automodule:: dnachisel.MutationSpace` /
`:members:`.

**Why partial anyway.** The random draws are hill-climbing moves, not output. In
`ObjectivesMaximizerMixin.optimize_by_random_mutations` (L62–101) each proposal is reverted unless
it improves the score:

```python
previous_sequence = self.sequence
self.sequence = self.mutation_space.apply_random_mutations(
    n_mutations=self.mutations_per_iteration, sequence=self.sequence)
if self.all_constraints_pass():
    new_score = self.objective_scores_sum()
    if new_score > score: score = new_score; stagnating_iterations = 0
    else: self.sequence = previous_sequence
else:
    self.sequence = previous_sequence
```

Same structure in `ConstraintsSolverMixin.resolve_constraints_by_random_mutations` (L84–127, whose
docstring is explicit: *"Solve all constraints by successive sets of random mutations… If all
constraints pass, the new sequence becomes the problem's new sequence"*). A DnaChisel run yields
**one** optimized sequence; the sampled variants are never returned as a set, and there is no
parameter that asks for a sampled variant library. That is the rubric's "stochastic optimization"
⇒ partial. Also present and separately insufficient: `random_dna_sequence(length, gc_share, probas,
seed)` and `random_protein_sequence(length, seed)` (`dnachisel/biotools/random_sequences.py`
L10–65) are de-novo generators (discriminator 1), and
`reverse_translate(protein_sequence, randomize_codons=False, …)`
(`dnachisel/biotools/sequences_operations.py` L38–59) is the synonymous draw seen in Mutation
Maker and CodonGenie.

**Disconfirmation.** I looked specifically for a way to make DnaChisel emit its sampled variants —
`grep -rn "pick_random_mutations\|random_variant\|all_variants\|def optimize" dnachisel/ --include=*.py`.
The only consumers of `apply_random_mutations` are the two optimizer loops above; the only consumers
of `all_variants` are `optimize_by_exhaustive_search` (L45–48) and
`ConstraintsSolverMixin` L63–66 — both also single-output. `grep -rniE "random_mutations|randomization_threshold|max_random_iters|apply_random" docs/ examples/`
returns exactly one hit, `examples/builtin_specifications_examples/example_max_random_iters.py:29
problem.max_random_iters = 10000` — a solver-tuning knob, not a library request. If a referee
requires that a documented public method returning `"a sequence with n random mutations applied"`
be sufficient, this cell becomes `yes`; that is the one defensible objection to this row and §4
records it.

### 3.9 ledidi — **partial**

`ledidi/ledidi.py` L78–79:

```python
def ledidi(model, X, y_bar, n_repeats=1, n_samples=None, return_designer=False,
    return_history=False, device='cuda', random_state=None, **kwargs):
```

L153–155 and L172–177:

> `n_samples: int or None, positive, optional`
> `    The number of samples to draw from Ledidi after the optimization process.`
> `    If None, draw one batch as defined by batch_size. Otherwise, draw the number of sequences specified.`
> `random_state: int or None, optional`
> `    A seed for the Gumbel-softmax sampling that makes the design procedure reproducible without mutating the global torch RNG.`

`docs/parameters.rst` L63–64:

> `- n_samples -- after optimization, draw this many designed sequences from the fitted weight matrix. Sampling is extremely fast once the weights are learned, so this is nearly free.`
> `- n_repeats -- run the whole design procedure this many times and stack the results, giving genuinely independent sets of edits in a single call.`

The draw is WT-relative — `Designer.forward` L494–499:

```python
logits = torch.log(X + self.eps) + self.weights
logits = logits.expand(self.batch_size, *(-1 for i in range(X.ndim-1)))
if self.random_state is None:
    return torch.nn.functional.gumbel_softmax(logits, tau=self.tau, hard=True, dim=1)
```

so `n_samples` + `random_state` genuinely return N random edited variants of the target `X`. What
makes it partial rather than yes: the distribution is `self.weights`, **fitted** by gradient descent
against `y_bar` (`docs/parameters.rst` L19 describes the only edit-count control as
*"Weight on the edit (input) loss. **The main knob to tune** -- lower values prioritize hitting the
target output, higher values prioritize making fewer edits."*). There is no mutation rate and no
count of edits — the number of edits is an emergent property of a penalty weight — and the sample
is drawn from the optimizer's posterior, not from the variant space. Rubric: stochastic
optimization ⇒ partial.

**Disconfirmation.** A `yes` would have required a documented rate or edit-count parameter, or
sampling from an un-fitted (uniform / rate-parameterized) distribution. I read the full parameter
table (`docs/parameters.rst`), the `ledidi()` and `Ledidi` docstrings (L78–380), and grepped
`ledidi/*.py` for `random|sample|seed|multinomial|gumbel|categorical`: the only RNG is the
Gumbel-softmax draw and `_gumbel_softmax_hard` (L29–70). A `no` would have been wrong — the drawn
edits are WT-relative and are the returned tensor.

### 3.10 tangermeme — **partial**

The package's mutagenesis function is exhaustive and has no sampling parameter.
`tangermeme/saturation_mutagenesis.py` L66–80:

```python
def saturation_mutagenesis(model, X, args=None, start=0, end=-1, batch_size=32,
    target=None, hypothetical=False, raw_outputs=False, dtype=None, device=None,
    verbose=False, func=None):
    """Performs in-silico saturation mutagenesis on a set of sequences.
    … return the predictions on the original sequences and each
    of the sequences with an edit distance of one on them.
```

No `n`, no rate, no `random_state`. `grep -rniE "mutation_rate|mut_rate|binomial|n_mutations|num_mutations" tangermeme/ --include=*.py --include=*.md --include=*.ipynb`
→ eight hits, all `test_plot.py` `..._per_position_...` test names. There is no random-mutagenesis
API.

What earns partial is `tangermeme/ersatz.py` L343–397:

```python
def randomize(X, start, end, probs=[[0.25, 0.25, 0.25, 0.25]], n=1, random_state=None):
    """Replace a region of the provided loci with randomly drawn sequence.
    … It will do this `n` times for each sequence in X …
    Importantly, this function does not shuffle the sequence in the specified
    region but replaces it with a random substitution.
```

This does draw new characters into a designated window of the user's target, with a count `n` and a
seed. Why it is not `yes`: the *whole* window is redrawn and the draw is WT-independent
(`probs` is a fixed alphabet distribution; nothing references the original base), so neither the
number of changes nor which positions change is a controllable mutational quantity, and a drawn
base may equal the WT. It is documented and used as a control/ablation primitive — `ablate.py`
L27/L85–88/L147–153 wires `func(X, start, end, n, random_state)` (default shuffle) through the same
plug, and `docs/api/ersatz.rst:5` lists it beside `shuffle, dinucleotide_shuffle`.

Everything else stochastic in the package fails discriminator 1 or 2:
`shuffle` (L430–492, `random_state.shuffle(idxs)`) and `dinucleotide_shuffle` (L529–591) are
permutations; `design/screen.py` L19–44 `func=random_one_hot` / *"Screen randomly generated
sequences and choose the best one"* takes a `shape`, not a target — de-novo generation;
`utils.random_one_hot` likewise (`docs/api/utils.rst:5`); `deep_lift_shap`/`pisa` `random_state`
seeds attribution backgrounds; `match.py` samples GC-matched genomic loci.

**Disconfirmation.** I searched specifically for a way to subsample the ISM variant space
(`saturation_mutagenesis` signature above; `grep` for `mutation_rate|n_mutations|binomial`) and for
a per-position mutation probability anywhere in the 28 modules — absent. Had `randomize` accepted a
per-position probability of change, or had `saturation_mutagenesis` accepted `n`/`random_state`, this
cell would be `yes`. Had `randomize` not existed (leaving only shuffles and de-novo generation),
this cell would be `no`.

### 3.11 biopython — **no**

Two exhaustive greps over the current package tree:

- `grep -rniE "^\s*def .*(mutat|randomiz|random_seq|shuffle)" --include=*.py Bio/ BioSQL/` → two
  hits, both phylogenetic: `Bio/Nexus/Trees.py:569 def randomize(` and
  `Bio/Phylo/BaseTree.py:758 def randomized(cls, taxa, branch_length=1.0, branch_stdev=None)`.
- `grep -rn "random\.choice\|np\.random\|numpy\.random\|random\.random\|random\.sample\|random\.shuffle\|random\.randint" --include=*.py Bio/`
  → nine hits in five files only: `Bio/Nexus/Nexus.py:1950` (bootstrap resampling of alignment
  columns), `Bio/Phylo/BaseTree.py:779, 789` and `Bio/Nexus/Trees.py:595, 607` (random tree
  construction), `Bio/Phylo/Consensus.py:557, 582, 592` (bootstrap column indices),
  `Bio/Graphics/ColorSpiral.py:115` (`jitter = random.random() * 2 * self._jitter - self._jitter`,
  a plotting colour jitter).
- `grep -rln "^import random\|^from random" --include=*.py Bio/` → exactly those five files.
- `grep -rniE "^\s*class .*(Mutat|Random)" --include=*.py Bio/` → only `Bio/SeqIO/_index.py`
  `*RandomAccess` proxies (random *file* access for indexed parsers).

Biopython has no sequence mutator and no random-sequence generator at all. `MutableSeq` supports
manual in-place editing, which is deterministic user code, not a sampling parameter.

**Disconfirmation.** A `partial` would have required a random-sequence or random-substitution
helper (the pydna/SeqPro situation); a `yes` a rate/count mutagenesis function. The two greps above
enumerate every RNG import and every RNG call in `Bio/`; nothing touches residue identity.

### 3.12 pydna — **no**

`src/pydna/utils.py` L444–542 is the entire RNG surface:

```python
def randomRNA(length, maxlength=None):
    ...
    return "".join([random.choice("GAUC") for x in range(length)])

def randomDNA(length, maxlength=None):
    ...
    return "".join([random.choice("GATC") for x in range(length)])

def randomORF(length, maxlength=None):
...
def randomprot(length, maxlength=None):
```

All four take a length, not a target sequence — de-novo generation (discriminator 1). This is the
control case that fixes the row's threshold: pydna literally ships "random DNA generation" and must
score `no`, otherwise the row degenerates.

**Disconfirmation.** `grep -rniE "mutagen|mutant|quikchange|site.directed" src/ docs` over the 41
modules in `src/pydna/` returns no mutagenesis function: the hits are a `loxP_mutant` string
constant (`cre_lox.py:72`), a SnapGene history-parser branch for a file field
(`snapgene_history_parser.py:345 elif node.operation == "primerDirectedMutagenesis":`, L356
`pcr_product.name = "mutagenesis_pcr_product"` — parsing someone else's cloning history, not
generating variants), one external notebook link in `docs/example_gallery.rst:31`, and GenBank
annotation text. pydna is an assembly/cloning simulator; it has no variant generator, hence nothing
to sample. A `partial` would have required a random substitution into a supplied `Dseqrecord`.

### 3.13 seqpro — **no**

Public API (`python/seqpro/__init__.py` L6, and `__all__` L13–42) exposes three RNG entry points,
all failing the test:

- `python/seqpro/_modifiers.py` L281–304 —
  `def random_seqs(shape, alphabet, seed=None)` / `"""Generate random nucleotide sequences."""` /
  `return seed.choice(alphabet.array, size=shape)`. Takes a `shape`, not a target: de-novo
  generation (discriminator 1).
- L39–108 — `def k_shuffle(..., seed: int | np.random.Generator | None = None)` /
  `"""Shuffle sequences while preserving k-let frequencies."""`. A permutation of the target's own
  characters; introduces no new character (discriminator 2).
- L144–162 — `def jitter(..., max_jitter, jitter_axes, seed=None)` /
  `"""Randomly jitter data from arrays, using the same jitter across arrays."""`. A coordinate
  shift, not a sequence edit.
- `python/seqpro/transforms/augmentation.py` L22–34 — `class Random: def __init__(self, p, *transforms, seed=None)`;
  `if self.rng.random() < self.p:` applies the wrapped transforms. `p` is the probability of
  applying a transform, not a per-position mutation rate.

**Disconfirmation.** `grep -rniE "\brandom|\bseed\b|shuffle|mutat|jitter"` over `python/seqpro/`
and `src/*.rs`: the Rust side is `kshuffle.rs` / `kshuffle_ref.rs` (`SmallRng::seed_from_u64`,
Wilson loop-erased random walk for the Eulerian-path shuffle) and `hashing.rs` (`rapidhash` seed) —
shuffle machinery and a hash seed. The `mutat` hits are `rag/_ops.py:52, 245` about mutating arrays
in place. No mutagenesis anywhere. Had `random_seqs` accepted a template sequence plus a rate, or
had a `mutate`/`substitute` function existed, this cell would move.

---

## 4. Risk register for the authors

1. **DnaChisel is the cell a referee is most likely to challenge.** Its own documentation contains
   the sentence *"The algorithm will use an exhaustive search when the size of the mutation space
   (=the number of possible variants) is above this threshold, and a (guided) random search when it
   is above"*, plus `mutations_per_iteration` = *"produce this many sequence mutations each
   iteration"*, and `MutationSpace.apply_random_mutations(n_mutations, sequence)` →
   *"Return a sequence with n random mutations applied"* is in the published API reference. The
   defensible distinction — and it must appear in the paper, not just here — is that these are
   search moves inside an optimizer that returns a single sequence; DnaChisel cannot be asked for a
   sampled variant library. Recommend a footnote saying exactly that.
2. **Six of thirteen cells are `partial`, for four different mechanisms.** They are not
   interchangeable and a bare "partial" under-informs:
   - synonymous-codon resampling at a rate or by usage frequency — codongenie, mutationmaker
   - random mutations as optimizer moves / post-optimization sampling — dnachisel, ledidi
   - WT-independent randomization of a designated window (control primitive) — tangermeme
   - random single-base repair of unwanted sites, written into the output oligo — mpradesign

   Recommend either a footnote key or a two-column split (see §5).
3. **MPRAnator's `yes` rests on a 1-commit, unreleased, unlicensed repo** whose documentation site
   (`genomegeek.com/MPRA/documentation/`) I did not rely on. The evidence above is entirely
   in-repo (`part3.py`, `forms.py:341`, `views.py:142/204`) and is solid, but the paper should
   expect an "is this maintained?" question rather than an "is this true?" one.
4. **MPRAnator offers no seed**, so its sampled output is not reproducible; PoolParty's is
   (`generate_library(seed=…)`, verified). If the paper cares about reproducible sampling, that is
   a separate, sharper row than this one.

## 5. Row verdict

**keep.** The row is scoreable on one scale: a single test applied to 13 heterogeneous tools
produced a 3-level spread with two unambiguous `yes` cells and five unambiguous `no` cells, and
every cell traces to quoted primary source. I did not need to invent tool-specific criteria.

Two definitional additions I recommend the paper adopt verbatim, because they are what actually
carried the row and without them the column is not reproducible by a referee:

- *"Random generation of new sequence (barcodes, spacers, primers, random backgrounds) and
  k-let-preserving shuffles are not counted; the parameter must draw substituted characters into
  positions of a user-supplied target."*
- *"Random mutations used only as move operators inside an optimizer that returns a single sequence
  are counted as partial, not yes."*

Optional stronger alternative, if the six partials are judged too coarse for the table: split into
two binary rows — **(a) "sampled random mutagenesis of a target sequence (rate / count / seed)"**:
yes = poolparty, mpranator; no = all others; and **(b) "any stochastic component affecting emitted
sequence"**: yes = poolparty, mpranator, mpradesign, mutationmaker, codongenie, dnachisel, ledidi,
tangermeme, oligopoolcalc, pydna, seqpro; no = valiant, biopython. Two clean binaries carry more
information than one row with six partials, at the cost of a column.

## 6. Provenance

Primary sources read this session. Shallow clones (`git clone --depth 1`) into the session
scratchpad `…/13e50e09-5cee-425a-a895-0f176887e7a6/scratchpad/src/`; nothing installed; nothing
outside this file written.

| tool | source read | notes |
|---|---|---|
| poolparty | `/mnt/c/…/poolparty-statecounter/poolparty/src/poolparty/base_ops/mutagenize.py`, `generate_library.py`, `docs/operations/mutagenize.rst`, `docs/tutorials/dms_gb1.rst` | plus read-only execution in `.venv` (`PYTHONDONTWRITEBYTECODE=1`) |
| valiant | `cancerit/VaLiAnT` @ `develop` — 85 `.py`, all `click.option` flags, `examples/sge/nxrl/*_config.json` | |
| mpranator | `hemberg-lab/MPRAnator` — `part3.py`, `part1.py`, `myfunctions.py`, `iliasApp/forms.py`, `iliasApp/views.py` | |
| mpradesign | `andrewGhazi/mpradesigntools` — `R/processVCFfast.R` (1459 L); `andrewGhazi/designMPRA` — `server.R`, `ui.R` | `server.R:17` sources the same script |
| oligopoolcalc | `ayaanhossain/oligopool` — `docs/api.md`, `docs/docs.md`, `oligopool/*.py` module inventory | |
| mutationmaker | `ra100/Mutation_Maker` — `backend/mutation_maker/*.py`, `backend/ssm_algorithm_comparison.py`, `api/*.py`, `README.md`, `docs/` | Merck URL is 404 |
| codongenie | `synbiochem/CodonGenie` — `codon_genie/codon_utils.py`, `codon_genie/codon_selector.py`, `main.py`, `setup.py`, `pypi.sh` | |
| dnachisel | `Edinburgh-Genome-Foundry/DnaChisel` — `MutationSpace/*.py`, `DnaOptimizationProblem/DnaOptimizationProblem.py`, `mixins/ObjectivesMaximizerMixin.py`, `mixins/ConstraintsSolverMixin.py`, `biotools/random_sequences.py`, `biotools/sequences_operations.py`, `docs/ref/core_classes.rst`, `examples/` | |
| ledidi | `jmschrei/ledidi` — `ledidi/ledidi.py`, `docs/parameters.rst` | |
| tangermeme | `jmschrei/tangermeme` — `ersatz.py`, `saturation_mutagenesis.py`, `ablate.py`, `design/screen.py`, `docs/api/ersatz.rst`, `docs/api/utils.rst`, `docs/whats_new.rst` | |
| biopython | `biopython/biopython` — full `Bio/` + `BioSQL/` RNG grep | |
| pydna | `BjornFJohansson/pydna` — `src/pydna/utils.py`, 41-module grep, `docs/` | |
| seqpro | `ML4GLand/SeqPro` — `python/seqpro/__init__.py`, `_modifiers.py`, `transforms/augmentation.py`, `src/*.rs`, `CHANGELOG.md` | |
