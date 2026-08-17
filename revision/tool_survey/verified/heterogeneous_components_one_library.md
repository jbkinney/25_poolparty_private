# Verified row audit: `heterogeneous_components_one_library`

**Auditor scope:** one row, all 13 tools, one threshold.
**Date:** 2026-08-10.
**Question:** Can structurally DIFFERENT component types coexist in ONE library specification producing ONE pooled output?

---

## 1. Operational test (fixed BEFORE any tool was opened)

Written to `verified/_optest_heterogeneous.txt` before the first source file was read:

> For each tool, locate a SINGLE invocation — one config/input file, one CLI run, or one script that uses
> only the tool's own API for BOTH generation AND combination — in which the tool's own primitives declare
> >= 2 STRUCTURALLY DIFFERENT component classes (differing in mutation mechanism or sampling regime, not
> merely in numeric parameters of one rule; e.g. exhaustive single substitution vs. sampled/combinatorial
> higher-order vs. unmutated wild-type control vs. indel/deletion class vs. a different mutator on a
> different sub-region), AND the tool's own code path emits them as ONE pooled output carrying shared
> per-oligo bookkeeping.
>
> - **YES** = a shipped example, shipped input file, or official docs page demonstrates exactly this.
> - **PARTIAL** = heterogeneity is reachable but only across SEPARATE invocations whose outputs the user
>   concatenates himself, OR the in-spec mixing is confined to variants of a single mutator class.
> - **NO** = the tool offers no tool-side mechanism to combine different component classes into one pooled
>   output, so any mixing would rest entirely on a generic language container (a Python list, an R `rbind`,
>   a shell `cat`). *A generic array container is explicitly NOT sufficient for YES.*
> - **UNKNOWN** = source and docs genuinely inaccessible.
>
> This test is applied identically to PoolParty. If PoolParty's only combiner were a bare Python list, it
> would score NO by the same clause that scores Biopython NO.

### 1.1 One sharpening, made during the mpradesign adjudication, then re-applied to all 13

Adjudicating mpradesign forced me to make the word **"declare"** operational, because mpradesign's mix
arises from the *contents* of a user-supplied VCF rather than from any keyword the user writes. The
sharpened form — the form actually applied to all 13 tools, including both anchors:

> **The specification must contain >= 2 *separately addressable declarations* of component generation**
> — i.e. the user can write down two different generative rules and change one without changing the other —
> **and the rules so declared must differ in kind, not only in their parameters.** A component class that
> the tool adds *automatically and unconditionally* (an obligatory wild-type twin) is not a declaration.
> A single global rule whose output happens to contain members of several shapes (a power set; a variant
> list that happens to include indels) is one declaration, not several.

This is a direct operationalisation of the pre-registered clause "not merely in numeric parameters of one
rule", and it coincides with the task rubric's own PARTIAL clause, *"limited to two component types of the
same kind."* It is stated here explicitly so a referee can re-run it. **Every one of the 13 cells below was
scored with this same sentence**, and the two anchors were re-checked against it after it was sharpened.

---

## 2. Result

| tool | prior | **verified** | one-line reason |
|---|---|---|---|
| poolparty | yes | **yes** | 5 named component pools -> `pp.stack` -> one 169-row table with `active_parent` provenance (run) |
| valiant | yes | **yes** | one shipped targeton row declares 3 sub-region action vectors; shipped `_meta.csv` pools 9 `mutator` classes |
| mpranator | partial | **partial** | one global combinations switch + automatic `REFERENCE`; cross-mutator mixing needs separate endpoints |
| **mpradesign** | **no** | **partial** | **3 distinct construct generators + automatic WT arm pooled into one TSV, but zero declarable rules** |
| oligopoolcalc | no | **no** | no member-generating primitive; `merge`/`join` are horizontal, `join` "never creates or drops rows" |
| mutationmaker | partial | **partial** | per-site declarations differ only in amino set/frequency; WT automatic; 3 workflows = 3 endpoints |
| codongenie | no | **no** | emits candidate ambiguous codons, never a sequence pool |
| dnachisel | no | **no** | optimises a single `self.sequence`; `MutationSpace` is a search space, not an emitted library |
| ledidi | no | **no** | one edit class (fixed-length Gumbel-softmax substitution); no pooling primitive |
| tangermeme | no | **no** | has 7 perturbation primitives but no tool-side pooling; `apply_product` crosses ONE func with model args |
| biopython | no | **no** | `MutableSeq` edits one sequence; no library primitive to pool |
| pydna | no | **no** | assembly-topology enumeration; no mutagenesis/variant-class primitive |
| seqpro | no | **no** | array transforms only; no library object, no component-class declaration |

**Changes: 1** (`mpradesign` no -> partial). Both anchors survive unchanged.

---

## 3. Cell-by-cell evidence

### 3.1 poolparty = **yes**

*Declaration (shipped docs, `docs/tutorials/dms_gb1.rst`):* four component classes, each separately
addressable, then one tool-side combiner:

```
L64    single_pool = orf_pool.mutagenize_orf(num_mutations=1, mutation_type="missense_only_first", ...)
L104   double_pool = orf_pool.mutagenize_orf(num_mutations=2, mutation_type="missense_only_first", ...)
L132   random_pool = orf_pool.mutagenize_orf(mutation_rate=0.1,  mutation_type="missense_only_first", ...)
L152   wt_pool     = orf_pool.repeat(times=100, prefix="wt").named("wt_pool")
L163   dms_pool    = pp.stack([single_pool, double_pool, random_pool, wt_pool])
```

`examples/dms_gb1.ipynb:403` — `DnaPool(id=5, name='pool[5]', op='op[5]:stack', num_states=547230)`;
`examples/dms_gb1.ipynb:555` — `| **Total** |  | **547,230** |`.
Exhaustive singles / exhaustive pairwise / rate-sampled random / unmutated replicates are four rules
differing in kind (enumeration vs. rate sampling vs. no mutation), not in parameters.

*The combiner is a tool primitive, not a Python list.* `src/poolparty/state_ops/stack.py:19`
`def stack(pools: Sequence[T], prefix=None, iter_order=None, cards=None) -> T` with docstring
*"Create a pool by stacking multiple input pools state-wise (disjoint union)."*, and
`stack.py:55` `design_card_keys = ["active_parent"]` — the pooled output records which component each
member came from.

*Behavioural verification (read-only venv, `PYTHONDONTWRITEBYTECODE=1`, five classes incl. two indel
classes):*

```
component sizes: [18, 135, 6, 7, 3] sum: 169
stacked num_states: 169
type: DataFrame rows: 169
columns: ['name', 'seq', 'mut_pos', 'mut', 'op[11]:stack.active_parent']
           name                                 seq mut_pos   mut  op[11]:stack.active_parent
0     single_00     TCCGACT<tag>ACATTG</tag>ATTCGGA    (0,)  (A,)                           0
18   double_000     TCCGACT<tag>AAATTG</tag>ATTCGGA    None  None                           1
153       del_0     TCCGACT<tag>-CATTG</tag>ATTCGGA    None  None                           2
159       ins_0  TCCGACT<tag>AAAGCATTG</tag>ATTCGGA    None  None                           3
168        wt_2     TCCGACT<tag>GCATTG</tag>ATTCGGA    None  None                           4
```

Five rules differing in kind (exhaustive single substitution, exhaustive pairwise substitution, deletion
scan, insertion scan, unmutated replicate), variable member length tolerated, one `DataFrame`, per-member
provenance column. Also in `README.md` lines 51-66 as the package's headline example.

*Disconfirmation:* a NO would have required `pp.stack` to be absent or to be a thin wrapper over a bare
Python list with no bookkeeping — I read `state_ops/stack.py` in full and ran the mixed stack; it is a
first-class `Operation` subclass (`class StackOp(Operation)`, `factory_name = "stack"`) that carries
`active_parent` cards. A PARTIAL would have required the components to be emitted as separate objects the
user concatenates — they are not; `generate_library()` returns one frame. **PoolParty is not softened
here: it is scored yes on the same clause that VaLiAnT is, and the clause that would have demoted it
(bare-container combination) was tested for explicitly.**

### 3.2 valiant = **yes** (anchor holds, and is stronger than stated)

*Declaration — shipped input file, verbatim,
`examples/sge/brca1/parameter_input_files/brca1_nuc_targeton_input.txt` line 2:*

```
chr17	-	43115634	43115878	43115726	43115779	"25,25"	"(1del,2del1,snv),(1del,snvre,inframe,stop,ala),(1del,2del0,snv)"	sgRNA_ex2
```

Three parenthesised action vectors = three sub-regions, each with its own separately addressable list of
mutator classes; `src/valiant/mutator_type.py` `class MutatorType(str, Enum)` gives the closed vocabulary
(`DEL='del'`, `SNV='snv'`, `SNV_RE='snvre'`, `IN_FRAME='inframe'`, `ALA='ala'`, `STOP='stop'`, `AA='aa'`).
The companion `brca1_pep_targeton_input.txt` proves the vectors are independently editable:
`"(),(aa,inframe),()"`.

*One pooled output — shipped expected output*
`examples/sge/brca1/brca1_nuc_output_exp_v3/chr17_43115634_43115878_minus_sgRNA_ex2_meta.csv`,
1000 data rows, 30 columns, `mutator` column values counted directly:

```
mutator counts: Counter({'custom': 369, 'snv': 312, 'snvre': 140, '1del': 104,
                         'inframe': 17, 'ala': 17, 'stop': 17, '2del1': 12, '2del0': 12})
vcf_alias counts: Counter({'': 631, 'clinvar_20201107': 327, 'gnomad_v3': 42})
```

Nine mutator classes in one file, each row tagged with `mutator` — and a ninth class, `custom`, is
user-supplied variants merged in from two external VCFs via `--vcf
reference_input_files/brca1_custom_variants_manifest.csv` (`run_brca1_nuc.sh`), tagged by `vcf_alias`.
Note for the paper: the pooled unit is **per targeton** (four `_meta.csv` files for the four spec rows),
not per run. That does not weaken the cell — one spec row is one specification unit — but it should be
stated accurately if the row is footnoted.

*Disconfirmation:* a PARTIAL would have required the action-vector products to land in separate files per
mutator, or the `--vcf` customs to be emitted separately. I enumerated the output directory and counted
the `mutator`/`vcf_alias` columns of a single file; all nine classes are in one table.

### 3.3 mpranator = **partial** (unchanged, uniform justification)

*What genuinely coexists in one submission* — `parseSNPs.py:30`
`def make_sequence_copies(SNPs, NamesL, SequencesL, CombinationsB):`, header comment `parseSNPs.py:29`:
*"if the SNP contains a small deletion or insertion (smaller than 10nt) we either remove part of the
sequence or we insert adenines in one edge"*. In one call the output lists accumulate substitution members
(`difference==0`), insertion members (`difference<0`, `" | removed "+...+" nucleotides from left edge"`),
deletion members (`difference>0`, `"| added "+str(difference)+" Adenines bases in the right edge"`), and an
unmutated reference member — `parseSNPs.py:189`:

```python
ExtraNames += [" ".join( [str(iS) for iS in NamesL[j]] )+"| REFERENCE "]
```

and in the combinations branch `parseSNPs.py:124` `Named[0] = " ".join(...) + "| REFERENCE "`.
All of it leaves as ONE FASTA body: `part1.py:210 def createMPRAResultOutput(...)`, `part1.py:212
response = ""`, `part1.py:289-292 response += outputHeader / response += outputSequence / return response`,
returned as a single `HttpResponse(content_type="text/plain")` (`iliasApp/views.py:379`).

*Why not yes.* There is exactly **one** declarable generative rule per submission and it is a global
boolean, not a per-region declaration: `generateCombinations` (`parseSNPs.py:18-26`, and identically
`part1.py:15-21`) is `for i in xrange(1, upTo): for acomb in itertools.combinations(theList, i)` — one
power-set rule of which the singles are simply the size-1 subsets. The `REFERENCE` member is appended
unconditionally by the tool, so it is not a declaration either. Mixing across mutator *kinds* — motif
placement vs. SNP substitution vs. random Transmutation — is impossible in one specification: they are
three independent Django views with three independent forms and three independent `HttpResponse` bodies
(`iliasApp/urls.py`: `^MPRA/$`/`^MPRAResults/$`, `^MPRA/SNPs/$`/`^MPRAResults/SNPs/$`,
`^Transmutation/$`/`^TransmutationResults/$`; `views.py:33 resultsView`, `views.py:150
part3RresultsView`, `views.py:319 mpraSnpResults`). That is the rubric's PARTIAL clause on both counts:
in-spec mixing confined to one rule, and true heterogeneity requiring separate runs the user merges.

*Disconfirmation:* a YES would have required either (a) a per-region or per-SNP action-vector field
analogous to VaLiAnT's, or (b) one endpoint that accepts both a motif list and a SNP list. I read
`urls.py`, all three result views, `forms.py`, `part1.py`, `part3.py` and `parseSNPs.py`; no such field and
no such endpoint exists. A NO would have required the single output to be homogeneous — it is not; the
`REFERENCE` and the indel branches are in the same response.

### 3.4 mpradesign = **partial** — CHANGED from `no` (the flagged closest call)

**Full credit first.** One `processVCF()` call on one VCF genuinely emits a pooled table containing three
structurally distinct construct classes plus a matched unmutated arm, with tool-side library bookkeeping.

*Dispatch to three distinct generators* — `R/processVCFfast.R:232-234`:

```r
isSNV = (snp$REF %in% c('A', 'C', 'G', 'T') && snp$ALT %in% c('A', 'C', 'G', 'T'))
isINS = snp$REF == '-'
isDEL = snp$ALT == '-'
```

branching at `L327 if (isSNV) {`, `L544 } else if (isINS) {`, `L764 } else if (isDEL) {`, with dedicated
generators `L36 generateInsConstruct = function(snpseq, mid, reverseGene, upstreamContextRange,
downstreamContextRange)` and `L63 generateDelConstruct = function(snpseq, refwidth,
upstreamContextRange)`. These are not cosmetic: `L252` comments *"# insertions are 1bp longer because the
ref is empty"*, and the INS branch overrides the allele semantics (`L608` *"mid = ifelse(type == 'ref', '',
snp$ALT), #This line is line is unique to the isINS block"*).

*Automatic matched wild-type arm* — `L392-393`:

```r
type = rep(c('ref', 'alt'), each = nper),
mid = ifelse(type == 'ref', snp$REF, snp$ALT),
```

so every variant is accompanied by `nper` separately barcoded unmutated constructs, labelled in a `type`
column that survives to the output (`L1252 select(ID, type, allele, snpIndex, totIndex, barcode, sequence,
site_fix_info, notes)`).

*One pooled output with tool-side bookkeeping* — `L1246-1249` `Reduce('rbind', .) %>% mutate(..., totIndex
= 1:nrow(.))`, one TSV at `L1258 write_tsv(successes, path = outPath)`, and library-wide barcode
disjointness established across the whole pool before generation (`L1213-1215` `shuffled_mers = mers[...]`
then `mutate(bcPools = split(shuffled_mers, ...))`).

*Documented as one input* — `README.md:56` verbatim: *"Insertions and deletions must encode the reference
and alternate alleles (respectively) as a dash character '-'."*, in the same VCF-format list as SNPs; and
`designMPRA/ui.R:123` verbatim: *"Insertions and deletions must encode the reference and mutant alleles
(respectively) as a dash character '-'."*.

**Why partial and not yes.** The specification contains **zero separately addressable declarations of a
generative rule.** The VCF is a list of individual variants, not of rules; `processVCF()`'s 17 arguments
(`L1099-1115`) contain no mutator, scheme, or action field; the ref arm is appended unconditionally, so it
is not something the user declares or can decline; and all three branches instantiate one and the same rule
— *place the listed allele into genomic context, compensating length* — under one hard-coded architecture,
identical for every member. `NAMESPACE` exports exactly two functions (`processVCF`,
`spread_and_fix_indels`), and no operation consumes a library and returns a library, so combining the
products of two `processVCF()` runs would be a user-written `rbind` — and would break the tool's own
library-wide barcode-disjointness guarantee, which is established per call at `L1213-1215`. That is
precisely the rubric's *"limited to two component types of the same kind"* -> **partial**.

**Why not `no`.** The prior `no` rested on the distinction "mixed variant *classes*, not mixed mutagenesis
*schemes*". That distinction is real, but it is what the split rows `exhaustive_single_scans`,
`sampled_random_mutagenesis` and `higher_order_combinatorial` exist to measure. *This* row asks only
whether structurally different components can coexist in one specification and one pooled output. Three
distinct construct generators plus an automatic WT arm, pooled by the tool with shared `totIndex` and
pool-wide barcode disjointness, is more than "one component type per run", which is what NO asserts.
Scoring it `no` while scoring mpranator `partial` for a materially weaker version of the same behaviour
(one global power-set switch plus an automatic `REFERENCE`) was the inconsistency this audit exists to fix.

*Disconfirmation:* a YES would have required a declarable rule field — a mutator column, a scheme argument,
a saturation/scan option, or an operation that composes two libraries. I read all 1459 lines of
`R/processVCFfast.R`, `NAMESPACE`, `DESCRIPTION`, `man/processVCF.Rd`, `README.md`, `designMPRA/ui.R` and
`designMPRA/server.R`, and searched for `motif|spacer|permut|combinator|saturat|scan|mutagen`; none exists.
A NO would have required either a single construct generator or per-class separate outputs; there are three
generators (`L36`, `L63`, plus the in-line SNV branch at `L327`) and one `Reduce('rbind', .)`.
R is not installed and installation is forbidden, so this cell is source- and docs-based, not behavioural —
hence confidence medium rather than high.

### 3.5 oligopoolcalc = **no**

The package contains **no primitive that generates library members**. Complete `### ` heading list from
`docs/api.md`: `barcode`, `primer`, `motif`, `spacer`, `background`, `merge`, `revcomp`, `join`, `final`,
`split`, `pad`, `compress`, `expand`, `index`, `pack`, `acount`, `xcount`, `lenstat`, `verify`, `inspect` —
architecture elements, degenerate compression, and analysis. Members arrive in `input_data` as *"CSV path or
DataFrame with `ID` column + DNA sequence columns"*.

Both combination primitives are **horizontal**, so they cannot pool different component classes:
`merge` — *"**Purpose**: Concatenate contiguous columns into a single column."* (`docs/api.md:427`);
`join` — *"**Purpose**: Join two oligo tables on `ID` and reconcile branch outputs into one design table."*
(`docs/api.md:530`), with the decisive note *"`join` never creates or drops rows; mismatched IDs are an
error"* and *"`ID` sets must match exactly across both inputs"*.

*Fairness note, recorded deliberately.* The repo ships
`examples/library-compressor/mutant_generator.py` (header: *"Author: Ayaan Hossain"*) with
`generate_single_mutants`, `generate_codon_variants`, `generate_multi_position_variants`,
`generate_variant_dataframe`, and `run_degenerate_demo.py:62 include_wildtype=True`. That is real
heterogeneous library generation — but it is a plain-Python demo helper in `examples/`, **not part of the
package**: `grep -rn "mutant_generator" oligopool/ pyproject.toml` returns 0 hits, and it is imported as
`from mutant_generator import (...)` from the example directory, not from `oligopool`. Under the
pre-registered clause that a generic container is not sufficient, this is the generic-container route and
does not lift the cell.

*Disconfirmation:* I searched for a vertical pooling primitive (`grep -rniE "def (stack|pool|concat_pools|
combine_pools|merge_pools|union)\("` over the package -> 0 hits) and read the `merge`, `join` and `final`
doc sections in full. Absent.

### 3.6 mutationmaker = **partial** (unchanged, uniform justification)

*Full credit.* Within one PAS job the spec carries genuinely per-site, separately addressable declarations —
`backend/mutation_maker/pas_types.py:116-120`:

```python
class PASMutationSite(JsonObject):
    """ List of mutations at a given amino acid position """
    position = IntegerProperty(required=True)
    # List of mutation. Each mutation is a amino acid IUPAC code, or a DNA codon.
    mutations = ListProperty(PASMutation, required=True)
```

and wild type is first-class in the emitted pool, not a user trick:
`generate_oligos.py:279-287` `if sum_of_probabilities != 1: wild_type_prob = 1 - sum_of_probabilities ... 
aminos_with_probabilities[wild_type_amino] = wild_type_prob`, surfaced per member as
`pas_types.py:165 wild_type = BooleanProperty(required=True)`. In QCLM the parental residue is unioned into
every site unconditionally — `mutation.py:147`:

```python
self.new_aminos = frozenset(old_aminos.union({m.new_amino for m in mutations}))
```

*Why not yes.* The per-site declarations do not differ **in kind** — every site declares the same rule
("this defined substitution set, at this frequency"), differing only in the amino set and the number. The WT
member is automatic, not declared. Two would-be class switches are **job-wide, not per-site**:
`ssm_types.py:215 degenerate_codon = StringProperty(required=True, default="NNS")` and
`pas_types.py:210 is_mutations_as_codons = BooleanProperty(required=True)`. And the three mutagenesis
workflows cannot be combined in one specification at all — three separate endpoints:
`api/api.py:40 @hug.post('/ssm', versions=1)`, `:48 @hug.post('/qclm', ...)`, `:56 @hug.post('/pas', ...)`,
mirrored at `api/server_fastapi.py:53/58/63` (`@app.post("/v1/ssm")`, `"/v1/qclm"`, `"/v1/pas"`).
Combining SSM saturation with QCLM defined substitutions therefore requires separate runs the user merges
= the rubric's PARTIAL clause.

*Disconfirmation:* a YES would have required a per-site or per-region *scheme* selector, or one endpoint
accepting more than one workflow. I read `pas_types.py`, `qclm_types.py`, `ssm_types.py`, `mutation.py`,
`generate_oligos.py`, `api/api.py` and `api/server_fastapi.py`; the only per-site freedom is the amino/codon
set and its frequency.

### 3.7 codongenie = **no**

`codon_genie/codon_selector.py:72 def optimise_codons(self, amino_acids, organism_id)` — the tool returns
candidate ambiguous codons scored for one amino-acid set at one position. No sequence library is emitted at
all, so there is no pooled output that could be heterogeneous. Complete Flask surface, `main.py`:
`:40 @app.route('/')`, `:46 @app.route('/organisms/')`, `:53 @app.route('/organisms/<term>')`,
`:61 @app.route('/codons')`.

*Disconfirmation:* a PARTIAL would have required some sequence-emitting route or a library builder. I
listed every `def`/route in `main.py`, `codon_genie/codon_selector.py`, `codon_utils.py`, `seq_utils.py`,
`client.py`; there is none, and no pooling primitive (`def stack|pool|...` -> 0 hits).

### 3.8 dnachisel = **no**

`DnaOptimizationProblem` holds and rewrites a **single** sequence:
`DnaOptimizationProblem/DnaOptimizationProblem.py:126 self.sequence = str(sequence.seq).upper()`,
`:129 self.sequence = sequence.upper()`, `:189 self.sequence = new_sequence`. `README.rst:20-31` describes
the package as *"a Python library for optimizing DNA sequences with respect to a set of constraints and
optimization objectives"*, composed of *"over 15 classes of sequence specifications"* — those are
constraints on one sequence, not component classes of a pool. `MutationSpace` is the optimiser's internal
search space, self-described `MutationSpace/MutationSpace.py:6-8` as *"Class for mutation space (set of
sequence segments and their variants)"* — never emitted as a library.

*Disconfirmation:* a PARTIAL would have required an operation returning a *set* of designed sequences plus
a way to combine sets. Searched for `def .*librar|variant library|enumerate.*variants` across
`dnachisel/dnachisel/` and for pooling primitives (0 hits). Absent.

### 3.9 ledidi = **no**

One edit class only. `ledidi/ledidi.py:78 def ledidi(model, X, y_bar, n_repeats=1, n_samples=None, ...)`;
the mechanism is fixed-length categorical resampling, `README.md:19` verbatim: *"at each iteration, Ledidi
draws one-hot encoded sequences one position at a time from a Gumbel-softmax distribution defined by
`log(X + eps) + W` ... When the drawn character is different from the character in `X` ... an "edit" is
being made."* That admits substitutions only — no insertion, deletion, or unmutated-replicate class.
The closest thing to heterogeneity is the affinity catalog (docstring `ledidi.py:115-118`: *"one can design
an affinity catalog by passing in a list of target values in `y_bar` ... an additional dimension is added to
the front of the returned tensor"*) — the same rule under different objectives, returned as a tensor axis,
not a pooled library with component-class bookkeeping.

*Disconfirmation:* a PARTIAL would have required a second edit class or a pooling primitive. Read
`ledidi/ledidi.py`, `wrappers.py`, `losses.py`, `pruning.py` and the full `README.md`; searched
`wild[_ -]?type|replicate` -> 0 hits, pooling primitives -> 0 hits.

### 3.10 tangermeme = **no**

Tangermeme *does* ship structurally different perturbation primitives — `ersatz.py` exports `insert`,
`substitute`, `multisubstitute`, `delete`, `randomize`, `shuffle`, `dinucleotide_shuffle` — and this cell
gives it full credit for them. What is absent is any **tool-side pooling** of their outputs into one
library. Each returns its own tensor (`ersatz.py:94 return torch.cat([X[:, :, :start], motif, X[:, :,
start:]], dim=-1)`; `ersatz.py:340 return torch.cat([X[:, :, :start], X[:, :, end:]], dim=-1)` — those
`torch.cat`s splice *within* one sequence, they do not pool classes). The only composition primitives,
`product.py`'s `apply_pairwise` and `apply_product`, apply **one** `func` across a cartesian product of
sequences and *model arguments*, not across several perturbation classes — docstring `product.py:41`:
*"Apply a function on the cartesian product between X and args."*, with the worked example being cell state
crossed with read depth. `saturation_mutagenesis` returns attributions, not a sequence pool ("Returns
------- attr: torch.Tensor / Processed attribution values"). Combining `insert` output with `delete` output
would be a bare `torch.cat` by the user = the generic-container route, explicitly excluded.

*Disconfirmation:* a PARTIAL would have required a function that returns one object holding members of >= 2
perturbation classes with class labels. I enumerated every public function in all 18 modules of
`tangermeme/tangermeme/`, read `product.py` and `saturation_mutagenesis.py`, and searched for pooling
primitives (0 hits) and for `wild[_ -]?type|replicate` (4 hits, all in agent-skill prose about comparing
models). Absent.

### 3.11 biopython = **no**

The only mutation facility is `Bio/Seq.py:2173 class MutableSeq(_SeqAbstractBaseClass)` — in-place editing
of one sequence. There is no mutagenesis, scan, variant-enumeration, or library primitive anywhere in
`Bio/`, hence nothing to pool; any mixing would be a plain Python list.

*Disconfirmation:* searched `Bio/` for `saturation.mutagenesis|mutagenize|deletion scan|variant librar`
(hits only in `Bio/Entrez/DTDs/*.dtd` — NCBI schema text, not code), for pooling primitives (0 hits), and
checked the current Tutorial & Cookbook tree (`Doc/Tutorial/`) for a library-design chapter — the only
`mutagenesis` hit is `Doc/Tutorial/chapter_motifs.rst`. Assessed as the current package (repo at
`c948960 Post Biopython 1.88 release version bump`), not the 2009 paper.

### 3.12 pydna = **no**

Pydna designs cloning, not variant libraries. `src/pydna/design.py` public surface: `primer_design`,
`assembly_fragments`, `circular_assembly_fragments`, `user_assembly_design`. `assembly.py` /`assembly2.py`
enumerate assembly *topologies* from a fragment graph (`assemble_linear`, `assemble_circular`,
`get_linear_assemblies`, `get_circular_assemblies`) — one generative rule (graph traversal), and the
products are assembly outcomes, not declared component classes. No mutagenesis primitive exists:
`grep -niE "def .*librar|combinatorial|mutagen"` over `src/pydna/*.py` returns only
`snapgene_history_parser.py:345` (`elif node.operation == "primerDirectedMutagenesis":` — parsing someone
else's SnapGene log) and one comment, `assembly2.py:1566`.

*Disconfirmation:* a PARTIAL would have required a variant-class primitive or a pooling operation over
sub-libraries. Listed all 39 modules, greped as above, and searched pooling primitives (0 hits) and
`wild[_ -]?type|replicate` (1 hit). Absent.

### 3.13 seqpro = **no**

Preprocessing only. `python/seqpro/__init__.py:4-6` exports `gc_content, length, nucleotide_content`,
`decode_ohe, decode_tokens, ohe, pad_seqs, tokenize`, and `bin_coverage, jitter, k_shuffle, random_seqs,
reverse_complement`. `k_shuffle` and `random_seqs` are two different sequence-producing operations, so the
raw material for heterogeneity exists — but there is no library object, no component-class declaration, and
no pooling primitive, so combining them is `np.concatenate` by the user = the excluded generic-container
route.

*Disconfirmation:* searched for pooling primitives (0 hits) and read `_modifiers.py`, `_encoders.py` and
`__init__.py`'s export list; there is no design/library concept at any level.

---

## 4. Row verdict: **keep**

The row is scoreable on one consistent scale. Applying one sentence to all 13 tools yields 2 yes / 3
partial / 8 no, both anchors survive verbatim, and the boundaries are drawn by evidence a referee can
re-check in one sitting.

**One recommendation for the paper, not a change to the row.** The wording "structurally DIFFERENT
component types" is what produced the original rater variance, because a VCF that happens to contain SNVs
and indels satisfies it on a literal reading. If the row is footnoted, tighten the criterion to the
sentence actually used here — *"two or more separately declared generative rules, differing in kind, pooled
by the tool into one output"* — and cite the two anchors as the operational definition. That wording makes
the mpradesign and mpranator cells self-evidently partial rather than arguable, and it is the wording
under which PoolParty's `pp.stack` and VaLiAnT's action vector are both unambiguously yes.

**Note on Block-A interaction.** `partial` for mpradesign is not a claim about mutagenesis power. Its
sibling rows `exhaustive_single_scans`, `sampled_random_mutagenesis` and `higher_order_combinatorial` are
all `no` for mpradesign (independently verified here: 17 `processVCF()` arguments contain no scan/rate/
combination field; `grep -E 'motif|spacer|permut|combinator|saturat'` over the source is empty), so the
table as a whole still shows the gap. The paper loses nothing by conceding this cell, and gains a row
whose threshold is stated and reproducible.

---

## 5. Everything checked (for reproducibility)

**Local, read-only:** `poolparty/src/poolparty/state_ops/stack.py`; `poolparty/docs/tutorials/dms_gb1.rst`;
`poolparty/examples/dms_gb1.ipynb`; `poolparty/examples/README.md`; `poolparty/README.md`; three scripted
runs against `poolparty-statecounter/.venv/bin/python` with `PYTHONDONTWRITEBYTECODE=1`
(`pp.stack` over 5 heterogeneous pools; `inspect.signature` on `insertion_scan`/`generate_library`).

**Repositories read (all cloned locally, commit pinned):** valiant `8796cc1` (develop) —
`examples/sge/brca1/parameter_input_files/*.txt`, `run_brca1_nuc.sh`,
`brca1_nuc_output_exp_v3/*_meta.csv` (parsed with csv.DictReader), `src/valiant/mutator_type.py`.
mpradesigntools `afd386e` — all 1459 lines of `R/processVCFfast.R`, `NAMESPACE`, `man/processVCF.Rd`,
`README.md`, `data/`. designMPRA `0cf56ee` — `ui.R`, `server.R`, `scripts/processVCF*.R`.
mpranator `9969790` — `parseSNPs.py`, `part1.py`, `part3.py`, `oligo.py`, `myfunctions.py`,
`iliasApp/urls.py`, `iliasApp/views.py`. mutationmaker `396c1c0` (ra100 fork) —
`backend/mutation_maker/{mutation,generate_oligos,pas,pas_types,pas_output,ssm_types,qclm_types}.py`,
`api/api.py`, `api/server_fastapi.py`. oligopool `b88fa39` — `docs/api.md` (all `###` sections;
`merge`/`join`/`final` in full), `oligopool/__init__.py`, module listing,
`examples/library-compressor/{mutant_generator.py,run_degenerate_demo.py}`.
codongenie `c26439f`; dnachisel `68c0930`; ledidi `adbca70`; tangermeme `8d732b8`; biopython `c948960`;
pydna `4e02f81`; seqpro `63a8439` — as cited per cell.

**Absence greps run identically across the 8 `no` tools:**
`grep -rniE "def (stack|pool|concat_pools|combine_pools|merge_pools|union)\(" <tool> --include=*.py`
-> 0 hits in every one of codongenie, dnachisel, ledidi, tangermeme, biopython, pydna, seqpro, oligopool;
`grep -rniE "wild[_ -]?type|wt_control|replicate" <tool>` -> every non-zero result inspected by hand
(tangermeme 4 = agent-skill prose; oligopool 32 = the examples-only `mutant_generator.py`; biopython 42 =
GenBank/PDB corpus text; pydna 1; seqpro 1).

**Leads used only to locate evidence, never cited as evidence:**
`revision/tool_survey/final/*.md`, `revision/tool_survey/v2/*.md`, `revision/tool_survey/ROWS_v2.md`.

**Not obtained:** no R interpreter is available and installation is forbidden, so mpradesign is scored from
source and docs, not from a run. That is the one cell where a behavioural check could still move the
evidence (though not, in my judgement, the value).
