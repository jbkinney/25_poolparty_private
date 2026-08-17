# Verified audit — row `lazy_generation`

**Row question:** Are library members produced on demand rather than fully materialized?
**Auditor pass:** single-rater, one threshold applied to all 13 tools.
**Date:** 2026-08-10

---

## 1. Operational test (fixed BEFORE any tool was inspected)

Written to `verified/_lazy_generation_test.txt` before opening a single source file:

> Using the tool's normal documented API, request the FIRST k=10 members of a library whose declared
> size is N (N large, e.g. 1e5–1e6). Count, in the primary source:
> (i) how many **member-level objects** the tool constructs before the 10th member is readable — where
> "member-level object" includes finished sequences/oligos/rows AND any per-member enumeration entry
> (index tuple, mutation record, coordinate row) that is one-per-library-member; and
> (ii) whether peak resident memory of that path is bounded independent of N.
>
> - **YES** = both O(k): no per-member structure of size N is built first.
> - **PARTIAL** = a lazy/incremental path exists but is gated on an O(N) per-member enumeration
>   (precomputed index / mutation table), OR the whole library is computed and only the WRITE is
>   incremental, OR laziness covers only a subset of the tool's design ops.
> - **NO** = no documented path yields any member before all N members are materialized.
> - **UNKNOWN** = source and docs genuinely inaccessible.

Two clarifying rules I had to add mid-audit (both stated so a referee can re-run them, both applied to
all 13 tools including PoolParty):

- **R1 — the tool must produce the members.** Generic lazy *reading* of sequences the user already has
  (a FASTA iterator) is not lazy library generation. Streaming credit requires the tool to have derived
  the members it emits.
- **R2 — degenerate N.** A tool that *samples* from an unenumerated space (no finite N) satisfies "O(k)"
  vacuously, because k *is* N. Such tools get PARTIAL, not YES: on-demand sampling is real partial
  laziness, but there is no k-of-N economy to demonstrate.

---

## 2. Symmetry test resolved explicitly: PoolParty is DEMOTED to `partial`

### 2a. PoolParty's lazy path is real and very strong

`generate_library(num_seqs=k)` / `to_file(num_seqs=k)` / `print_library(num_seqs=k)` are documented
k-of-N parameters, and the engine is genuinely per-row:

`src/poolparty/generate_library.py:115` and `:265`

```python
    while len(rows) < num_seqs:
...
    pool.state.value = global_state % pool.state.num_values
```

Most operations derive `num_states` arithmetically, never by enumeration:

`src/poolparty/utils/scan_utils.py:9-65` — `build_scan_cache(...) -> int` returns
`num_states = max_start + 1`; despite the name it builds no cache.

`src/poolparty/base_ops/get_kmers.py:200-206` — kmers are **unranked on the fly**, no table:

```python
    def _value_to_kmer(self, value: int) -> str:
        result = []
        remaining = value
        for _ in range(self.length):
            result.append(dna_utils.BASES[remaining % self.alpha_size])
            remaining //= self.alpha_size
        return "".join(reversed(result))
```

Measured (read-only venv `/mnt/c/.../poolparty-statecounter/.venv/bin/python`, `PYTHONDONTWRITEBYTECODE=1`,
`tracemalloc`, script `scratchpad/pp_lazy2.py`, `pp_final2.py`):

| pool | N | construct | peak at construct | k=10 | k=1000 |
|---|---|---|---|---|---|
| `deletion_scan` (300 nt, sequential) | 300 | 0.001 s | 0.0 MB | O(k) | — |
| `del_scan × ins_scan` (3000 nt) | 9,003,000 | 0.014 s | 0.2 MB | 0.07 s | — |
| `del_scan × ins_scan × ins_scan` | **27,036,009,000** | **0.001 s** | **0.0 MB** | 0.0045 s | 0.2198 s |

`state.num_values == 27036009000` is an **exact Python `int`**, and random access deep into the space is
O(1):

```
init_state=0:           0.0368s  names=['d_0000.i1_0.i2_0', ...]
init_state=27036008990: 0.0017s  names=['d_2999.i1_3000.i2_2993', ...]
```

`to_file` streams (`src/poolparty/pool_mixins/export_mixin.py:442-447`):

```python
            while written < target_count:
                remaining = target_count - written
                this_chunk = min(chunk_size, remaining)
...
                df = self.generate_library(
                    num_seqs=this_chunk,
```

### 2b. But five documented operations build one entry per member at construction

`src/poolparty/base_ops/mutagenize.py:354-379` (`_build_caches`, called from `__init__` at `:234`/`:241`
whenever `mode="sequential"`):

```python
            cache = []
            for positions in combinations(range(num_positions), self.num_mutations):
                num_mut_patterns = alpha_minus_1**self.num_mutations
                for mut_pattern in range(num_mut_patterns):
...
                    cache.append((positions, tuple(reversed(mut_indices))))
...
        self._sequential_cache = cache
```

and `_compute_core` then just indexes it (`:569`):

```python
            rel_positions, mut_indices = self._sequential_cache[state % len(self._sequential_cache)]
```

This is, verbatim, the rubric's PARTIAL clause: *"a generator that iterates over a PRE-COMPUTED index of
all members."*

Same pattern in four more places:
- `src/poolparty/orf_ops/mutagenize_orf.py:423-438` (`_build_caches`, called at `:302`)
- `src/poolparty/base_ops/recombine.py:300-319` (`_build_cache`) — `cache.append((breakpoint_combo, pool_assignments))`
- `src/poolparty/region_ops/region_multiscan.py:272-288` — `self._sequential_cache = enumerate_multiscan_combinations(...)`; `return len(self._sequential_cache)`, i.e. `num_states` is obtained **by enumerating every member**
- `src/poolparty/base_ops/get_barcodes.py:314-317` class docstring, and `:412-414`:

```python
    """Generate constrained DNA barcodes via greedy random algorithm.

    All barcodes are pre-generated at construction time and stored.
```
```python
        self._barcode_strings: list[str] = [
            self._pad(bc, max_length) for bc in raw_barcodes
        ]
```

Measured cost of these caches (script `scratchpad/pp_gb1.py`, `pp_lazy.py`):

| pool | N | construct | peak at construct | `len(_sequential_cache)` |
|---|---|---|---|---|
| GB1 `mutagenize_orf(num_mutations=1, sequential)` — docs tutorial | 1,045 | 0.005 s | 0.1 MB | **1045** |
| GB1 `mutagenize_orf(num_mutations=2, sequential)` — docs tutorial | 536,085 | **1.892 s** | **64.8 MB** | **536085** |
| `mutagenize(num_mutations=2, sequential)` 300 nt | 403,650 | 1.767 s | 51.0 MB | 403650 |
| `get_barcodes(num_barcodes=24000, length=10)` | 24,000 | 0.803 s | 6.7 MB | 24000 strings |

`len(_sequential_cache) == N` exactly. For the paper's own DMS tutorial
(`docs/tutorials/dms_gb1.rst:104-114`, `print(double_pool.num_states)  # 536085`) requesting k=10 members
costs 1.89 s and 64.8 MB **before row 1 exists** — i.e. O(N), not O(k).

`mode="sequential"` is not a fringe path: it is used in the front-page README quick example
(`README.md:57-60`, `template.mutagenize(num_mutations=1, region="tag", ..., mode="sequential", ...)`) and
in both shipped tutorials.

Note also that `docs/pool.rst:6-9` states *"no sequences are generated until you explicitly request
them"*. That is accurate for `mutagenize` (an index, not sequences, is precomputed) but **not** for
`get_barcodes`, which pre-generates 24,000 barcode strings at construction. Worth a one-line correction.

### 2c. The threshold decision, stated once

I chose the rubric-literal reading: **any documented core enumeration operation that builds one entry per
library member before the first member is readable demotes the tool to `partial`.** PoolParty has five such
operations, one of which is the flagship DMS use case, measured at 536,085 entries / 64.8 MB. **PoolParty =
partial.** This is a demotion from the prior `yes`.

The same threshold applied to DnaChisel also yields `partial` (§3), so the row does **not** become an
unjustified win for a competitor — but PoolParty no longer wins it either. See §6 for the consequence and
for the ~10-line change that would legitimately earn `yes`.

---

## 3. dnachisel — `partial` (unchanged, but for corrected reasons)

`MutationSpace.all_variants` is a **true generator**, not an index walk
(`dnachisel/MutationSpace/MutationSpace.py:132-164`):

```python
        variants_slots = [
            [
                (choice_.segment, v.encode())
                for v in sort_variants_by_distance_to_current(choice_)
            ]
            for choice_ in self.multichoices
        ]
        for variants in itertools.product(*variants_slots):
            new_sequence[choice_start:choice_end] = encoded_segment
            for (start, end), variant in variants:
                new_sequence[start:end] = variant
            yield new_sequence.decode()
```

`variants_slots` is O(L·|alphabet|) — per **position**, not per member — and `itertools.product` is lazy.
`islice(space.all_variants(seq), 10)` is therefore strictly O(k). On this specific axis DnaChisel is
**more** lazy than PoolParty's sequential `mutagenize`, which materializes the equivalent list.

Size is arithmetic, not enumerated (`MutationSpace.py:96-104`):

```python
    @property
    def space_size(self):
        """Return the number of possible mutations."""
        if len(self.multichoices) == 0:
            return 0
        choices = [len(choice.variants) for choice in self.multichoices]
        return np.exp(min(100, np.log(choices).sum()))
```

It is documented API: `docs/ref/core_classes.rst:45-48` publishes
`.. automodule:: dnachisel.MutationSpace` with `:members:`. And lazy *consumption with early exit* is in
the normal solver path (`DnaOptimizationProblem/mixins/ConstraintsSolverMixin.py:63-79`): it iterates
`all_variants` and `return`s on the first variant that passes.

**Why not `yes`:** DnaChisel has its own eager per-member materialization, in the same shape as
PoolParty's. `dnachisel/SequencePattern/DnaNotationPattern.py:34-41`:

```python
    def all_variants(self):
        """Return all ATGC sequence variants of a sequence"""
        return [
            "".join(nucleotides)
            for nucleotides in itertools.product(
                *[IUPAC_NOTATION[n] for n in self.sequence]
            )
        ]
```

called with an explicit comment in `dnachisel/builtin_specifications/EnforceChoice.py:40-41`:

```python
        # PRECOMPUTE ALL VARIANTS
        choices = [variant for choice in choices for variant in choice.all_variants()]
```

Also: `space_size` is a **float capped at `exp(100)`**, not an exact count, and there is no documented
k parameter — the user must wrap the generator in `islice` themselves. Lazy core enumeration + eager
per-member materialization in a documented op ⇒ **partial**, exactly as for PoolParty.

Checked and did not find: any `num_seqs`/`max_variants`/`limit` parameter on a public library-producing
call (`grep -rn "all_variants\|space_size" --include=*.py`), and any generator-returning public API other
than `MutationSpace.all_variants`.

---

## 4. tangermeme, ledidi, biopython — the other three flagged cells

### tangermeme — `partial` (unchanged)

Eager side: every direct variant-constructing function materializes the full set.
`tangermeme/saturation_mutagenesis.py:~218-222`:

```python
		X_ = X[i].repeat(n_edits, 1, 1)
		rows = torch.arange(n_edits)
		X_[rows, :, edit_positions] = 0
		X_[rows, edit_chars, edit_positions] = 1
```

`n_edits` = every single-base edit of the sequence, i.e. all N members in one tensor before any prediction.
`ersatz.py:489-498` and `:679-688` end in `torch.stack(X_shufs)`; `ersatz.py:422`
`return torch.stack(X_rands).permute(1, 0, 2, 3)`.

Lazy side (genuine, and it is what earns partial): `tangermeme/product.py:300-317` iterates a Cartesian
design space holding only `batch_size` members at a time —

```python
		for x in tqdm(itertools.product(X, *args), disable=not verbose):
			X_.append(x[0])
...
			if len(X_) == batch_size:
				y_ = _apply(func, model, X_, args=args_, batch_size=batch_size,
...
				X_, args_ = [], [[] for _ in args]
```

and `tangermeme/design/screen.py:155-183` keeps only `n_best` candidates in a heap while drawing fresh
batches, so peak memory is O(batch_size + n_best). But there is no k-of-N parameter: `apply_product`
always traverses the whole product, and outputs accumulate O(N). ⇒ **partial**.

Checked and did not find (`grep -rn "yield \|Dataset\|DataLoader\|itertools" --include=*.py tangermeme/`):
any `torch.utils.data.Dataset`/`DataLoader`, any generator-returning public variant function.

### ledidi — `partial` (unchanged; rule R2)

Documented on-demand sampling, `README.md:146`:

> "once the Ledidi weight matrix is learned, edited sequences can be sampled extremely quickly using the
> forward function."

`ledidi/ledidi.py:467-505` — `forward(X)` expands one sequence to `batch_size` and draws Gumbel-softmax
samples; nothing is precomputed. The `n_samples` path is O(k) by construction
(`ledidi/ledidi.py:236`):

```python
				X_bar_ = torch.cat([designer(X) for _ in range(n_iter)], dim=0)[:n_samples]
```

But there is **no enumerable library**: the design space is sampled, never indexed, so "k of N" is
degenerate (k = N). Every requested batch is materialized in full as one tensor
(`fit_transform` returns `best_sequence`, shape `(batch_size, ...)`), and a fixed `max_iter` gradient
optimization runs first. Under rule R2 ⇒ **partial**. Scoring this `yes` would rank a sampler with no
library above tools that actually enumerate one.

### biopython — `no` (CHANGED from `partial`)

Biopython has **no library-generation operation at all**. Searched:
`grep -rn "def mutate\|def mutagen\|combinatorial\|def all_variants\|def variants\|def enumerate" --include=*.py Bio/`
→ only `Bio/PDB/PICIO.py:909 def enumerate_atoms` and a docstring in `Bio/PDB/cealign.py`.
`grep -rn "itertools.product\|from itertools import product" --include=*.py Bio/` → one hit,
`Bio/PDB/mmtf/mmtfio.py:79` (chain-ID generation). Nothing mutagenic, nothing combinatorial over sequences.

What Biopython *does* have, and what I deliberately refused to count under rule R1:

- `Bio/SeqIO/Interfaces.py:228-237` — a genuinely incremental writer:
  ```python
    def write_records(self, records):
        """Write records to the output file, and return the number of records.

        records - A list or iterator returning SeqRecord objects
        """
        count = 0
        for record in records:
            self.write_record(record)
            count += 1
        return count
  ```
  This streams, but the records are produced by *user code*, not by Biopython.
- `Bio/SeqIO/__init__.py:620` — `parse` "Turn a sequence file into an iterator returning SeqRecords."
  Lazy *reading*, not generation.
- `Bio/Align/__init__.py:4076-4087` — the strongest laziness in the package, and it is not a sequence
  library:
  > "Implements an iterator over pairwise alignments returned by the aligner. … Note that pairwise
  > aligners can return an astronomical number of alignments … we therefore recommend to first check the
  > number of alignments, accessible as len(alignments), which can be calculated quickly even if the
  > number of alignments is very large."

  Structurally this is the PoolParty pattern (lazy iteration + cheap exact size) applied to alignments.

The prior `partial` was, I believe, credit for `SeqIO.parse`/`SeqIO.write` — generic sequence I/O. Under
rule R1 that is the same over-credit this audit exists to remove. **`no`: Biopython produces no library
members, therefore produces none on demand.** (The vocabulary has no "n/a"; `no` is the honest carrier.)

---

## 5. Remaining tools

### valiant — `partial` (CHANGED from `no`)

All variants are materialized in a Python list first:
`src/valiant/mutators/snv.py:44-49` `def get_variants(self, seq: Seq) -> list[Variant]: ... return list(chain.from_iterable(...))`;
`src/valiant/mutator.py:133-141` `def _get_variants(...) -> list[PatternVariant]: return [ ... ]`;
`src/valiant/cdna_proc.py:106` `insert_pattern_variants(conn, list(map(get_oligo, vars)))`;
`src/valiant/queries.py:192-205` `cur.executemany(sql_insert_pattern_variants, [...])`.

**But** the export is genuinely row-at-a-time out of SQLite, not `DataFrame.to_csv` of a materialized
object — `src/valiant/meta_table.py:488-490`:

```python
                it = cur.execute(sql_select_meta)
                while r := it.fetchone():
                    mr = MetaRow(*r)
```

with each record written immediately inside that loop into the already-open `meta_fh`. That is precisely
the rubric's PARTIAL clause *"streaming EXPORT only (whole library computed, written incrementally)"*.
There is no k-of-N knob (the CLI always processes everything), and `unique_oligos = defaultdict(list)`
(`meta_table.py:358`) still accumulates O(N) inside the loop. ⇒ **partial**, not `no`.

Checked: `grep -rn "yield " --include=*.py src/valiant/` → all other hits are `@contextmanager` file/db
handles (`db.py:133,140`, `vcf_writer.py:111,131`, `loaders/fasta.py:43`, `loaders/vcf.py:39`) plus
`loaders/csv.py:53` (input row reader).

### oligopoolcalc — `no`

Every module signature is whole-DataFrame in / whole-DataFrame out, e.g. `oligopool/final.py:14-17`
`def final(input_data: str|pd.DataFrame, output_file: str|None=None, verbose: bool=True) -> Tuple[pd.DataFrame, dict]`,
and `:101-103` `outdf = pd.DataFrame(index=indf.index); outdf['CompleteOligo'] = ut.get_df_concat(df=indf)`.
Barcodes are preallocated for the whole target count —
`oligopool/base/core_barcode.py:1000-1002`:

```python
    store = np.zeros(            # Store Encoded Barcodes
        (targetcount, barcodelen),
        dtype=np.float64)
```

filled at `:1095` `store[count, :] = barcode`. Checked for laziness:
`grep -rn "yield |generator|chunk" --include=*.py oligopool/*.py` → only `cli.py:1058` (a filename-suffix
helper) and `revcomp.py` column slicing. The FASTQ *counting* modules stream reads, but that is readout
analysis, not library generation.

### pydna — `no`

`src/pydna/assembly.py:339-345` materializes all graph paths
(`linearpaths = list(itertools.chain(nx.all_simple_paths(...), ...))`), accumulates every product in
`lps`, and `:387-400` returns `sorted((Contig.from_string(...) for lp in lps.values()), key=len, reverse=True)`
— a sort forces full materialization. `src/pydna/assembly2.py:1820-1825`
`return [assemble(self.fragments, a) for a in assemblies]`.

Notably pydna computes the product count cheaply (`assembly2.py:1554-1564`
`get_possible_assembly_number`) and then **refuses** rather than streams
(`:1473-1477`): `raise ValueError(f"Too many assemblies ({possible_assemblies} pre-validation) to assemble")`.
`max_assemblies: int = 50` is a guard, not a k-of-N parameter — the opposite of laziness.

### seqpro — `no`

`python/seqpro/_modifiers.py:281-304` `def random_seqs(shape, alphabet, seed) -> NDArray[np.bytes_]: ... return seed.choice(alphabet.array, size=shape)`
— the full (N, L) array. `k_shuffle` (`:39-47`) returns `NDArray`; `:108 return shuffled`.
Disconfirmation: SeqPro *does* have a dask-capable path, `python/seqpro/xr/__init__.py:42-48`
`out = xr.apply_ufunc(..., dask="parallelized", ...)`, but `__all__ = ["bin_coverage", "ohe"]` — encoding
and coverage binning only, no design op; and the commented-out `jitter` at `:158-160` records "This
function doesn't seem possible to implement with xarray.apply_ufunc". No lazy design path.

### mutationmaker — `no`

`backend/mutation_maker/degenerate_codon.py:372-377` materializes exactly N and asserts it:

```python
        degenerate_site_sequence_codons: List[List[str]] = []
        all_combinations = itertools.product(*degenerate_codons_lists)
        for comb in all_combinations:
            # itertools product creates tuples-> need to convert it to list
            degenerate_site_sequence_codons.append(list(comb))
        assert len(degenerate_site_sequence_codons) == np.prod([len(site) for site in degenerate_codons])
```

The lazy `itertools.product` is consumed immediately into a list. Backend returns JSON payloads; no
streaming, no k parameter.

### codongenie — `no`

`codon_genie/codon_selector.py:71-76`:

```python
    def optimise_codons(self, amino_acids, organism_id):
        '''Optimises codon selection.'''
        req_amino_acids = set(amino_acids.upper())

        codons = [CODONS[amino_acid] for amino_acid in req_amino_acids]

        results = [self.__analyse(combo, organism_id, req_amino_acids)
                   for combo in itertools.product(*codons)]
```

List comprehension over the full product; all results analysed and returned before anything is readable.
(Source is on GitHub, so this is `no`, not `unknown`.)

### mpranator — `no`

`part1.py:6-12`:

```python
def generatePermutations(theListOfLists):  # Used to make all permutations
    allpermutations = []
    for i in theListOfLists:
        for subset in itertools.permutations(i):
            if subset not in allpermutations:
                allpermutations.append(subset)
    return allpermutations
```

(the `if subset not in allpermutations` membership test is itself O(N) per insert). `oligo.py:74-78`
accumulates every designed oligo in `allResults` and `:78 return allResults`.
`part1.py:52-88 barcode_generator(...)` ends `return barcode_storage`.

### mpradesign — `no`

`R/processVCFfast.R:1242-1249`:

```r
    successes = processed %>%
      filter(!failed) %>%
      .$seqs %>%
      Reduce('rbind', .) %>%
      mutate(.,
             constrseq = constrseq %>% unlist,
             totIndex = 1:nrow(.)) %>%
```

`Reduce('rbind', .)` binds the entire library into one data frame, and `totIndex = 1:nrow(.)` requires the
full row count, before `write_tsv(successes, path = outPath)` at `:1259`. Same pattern at `:1145`,
`:1284`.

---

## 6. Consequence for the paper, and `row_verdict = reword`

Under one consistent threshold the column becomes **5 × partial / 8 × no, with zero `yes`**:

| tool | prior | verified |
|---|---|---|
| poolparty | yes | **partial** |
| valiant | no | **partial** |
| biopython | partial | **no** |
| dnachisel | partial | partial |
| ledidi | partial | partial |
| tangermeme | partial | partial |
| mpranator, mpradesign, oligopoolcalc, mutationmaker, codongenie, pydna, seqpro | no | no |

That is internally consistent and I stand behind every cell — but as a table row it is weak: no tool
attains the top value, and PoolParty ties DnaChisel, tangermeme, ledidi and VaLiAnT. It also required me
to invent rule R2, because the row as written cannot distinguish "lazy enumeration of a declared finite
library" from "sampling an unenumerated space" — a sampler satisfies O(k) vacuously.

**Two options for the authors.**

**(A) Reword the row** so it measures what PoolParty demonstrably and uniquely does. Proposed replacement:

> `unmaterialized_library_addressing` — *Can a library far larger than memory be declared, its exact size
> reported, and an arbitrary k-member slice generated at any offset, with work and peak memory O(k) and no
> per-member structure of size N, through a documented parameter of the tool's primary generation call?*

Scoring under that wording, with evidence already in hand:
- **poolparty = yes.** `state.num_values == 27036009000` (exact `int`); DAG construction 0.001 s / 0.0 MB;
  `generate_library(num_seqs=1000)` = 0.22 s; `init_state=27036008990` returns members
  #27,036,008,990–992 in 0.0017 s. Caveat to state in prose: sequential `mutagenize`/`mutagenize_orf`/
  `recombine`/`region_multiscan` and `get_barcodes` precompute their state index at construction.
- **dnachisel = partial.** `islice(all_variants(seq), k)` is O(k), but `space_size` is a float capped at
  `exp(100)` (not an exact size), there is no documented k parameter, and there is no offset/random access.
- **all others = no.** None reports an exact unenumerated size *and* addresses an arbitrary slice; pydna
  raises rather than slices; the rest materialize.

That row is discriminating, fair, and survives a DnaChisel author as referee.

**(B) Or earn the `yes` in code.** The fix is small and PoolParty already contains the technique:
`get_kmers._value_to_kmer` (`base_ops/get_kmers.py:200-206`) unranks a state integer into a member with no
table. Doing the same for `mutagenize`/`mutagenize_orf` — combinatorial unranking of
`(positions, mut_indices)` from `state`, instead of `cache.append(...)` over `combinations(...)` — removes
the only O(N) structure on the DMS path. The state count is *already* computed arithmetically
(`mutagenize.py:351-353`, `num_combinations = comb(num_positions, self.num_mutations) * (alpha_minus_1**self.num_mutations)`),
so the cache is redundant for sizing. `get_barcodes` cannot be fixed this way (greedy constrained search
is inherently sequential) and should simply be documented as eager.

Either option is honest. What the paper should **not** do is ship `lazy_generation: poolparty = yes,
dnachisel = partial` — `MutationSpace.all_variants` is a pure generator where PoolParty's
`_sequential_cache` is a materialized list, and a referee who opens both files will see it.

---

## 7. Files, greps and commands used

Local, read-only:
- `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/src/poolparty/` —
  `generate_library.py`, `base_ops/mutagenize.py`, `base_ops/get_barcodes.py`, `base_ops/get_kmers.py`,
  `base_ops/recombine.py`, `orf_ops/mutagenize_orf.py`, `region_ops/region_scan.py`,
  `region_ops/region_multiscan.py`, `utils/scan_utils.py`, `pool_mixins/export_mixin.py`,
  `state_ops/state_shuffle.py`, `scan_ops/*`, `multiscan_ops/*`, `base_ops/from_iupac.py`,
  `base_ops/from_seqs.py`, `base_ops/shuffle_seq.py`
- greps: `_build_caches|_sequential_cache|combinations(|product(|_cache = []`;
  `num_states\s*=|_cache|= \[\]|list(|combinations|product(` per-op
- docs: `README.md`, `docs/index.rst`, `docs/pool.rst`, `docs/quickstart.rst`, `docs/tutorials/dms_gb1.rst`
- behavioural runs with `/mnt/c/.../poolparty-statecounter/.venv/bin/python`, `PYTHONDONTWRITEBYTECODE=1`:
  `scratchpad/pp_lazy.py`, `pp_lazy2.py`, `pp_gb1.py`, `pp_final.py`, `pp_final2.py` (tracemalloc +
  perf_counter; no writes outside scratchpad)

Cloned to `/tmp` (`git clone --depth 1`) and read as source:
DnaChisel, tangermeme, ledidi, oligopool, VaLiAnT, pydna, SeqPro, MPRAnator, mpradesigntools,
Mutation_Maker (ra100 fork), CodonGenie, biopython.
