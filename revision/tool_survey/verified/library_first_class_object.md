# Verified audit — row `library_first_class_object`

Auditor pass: single-rater, one threshold applied to all 13 tools.
Date: 2026-08-10.

**Row verdict: REWORD.** Cells are scored below under the rubric exactly as written (so the
authors have a usable row if they keep it), but the rubric's criteria (ii) and (iii) are *other
rows of this same table* — `library_algebra` and `lazy_generation` — folded into this one. That
re-merges a split the authors themselves made deliberately, triple-counts PoolParty's advantage,
and produces two demonstrable false equivalences. The fix is to narrow this row back to its own
one-line definition, not to invent new rows. See
[§4 Why the row must be reworded](#4-why-the-row-must-be-reworded).

---

## 1. Operational test (fixed before any tool was opened)

For each tool, using ONLY its public programmatic interface as documented and as implemented in
its primary source: attempt to bind the result of the library-construction call to a variable,
then check three things against primary sources:

- **(i)** the bound value is an instance of a purpose-built library/pool/design type that carries
  library semantics — NOT a bare ndarray / DataFrame / ragged array / torch tensor /
  list-or-dict of sequences, and NOT a single-sequence record;
- **(ii)** at least one OTHER operation in the same library-design API is documented or
  type-annotated to ACCEPT that same type as an argument (library → library, or library →
  derived library composition);
- **(iii)** that object can report the library's TOTAL MEMBER COUNT from the design specification
  alone — no file written, and without the members having been materialised/enumerated.

`yes` = all three. `partial` = an object exists but one or two of (i)–(iii) fail.
`no` = no bindable library object at all, OR the object is a generic container (rubric's explicit
exclusion). `unknown` = source genuinely inaccessible.

A referee re-runs this by, per tool: (a) grep `^class ` over the package for a type whose
instances denote a set of designed sequences; (b) grep the package for a function/method
signature that takes that type as a parameter; (c) grep for a size/count accessor and read its
body to confirm it does no I/O and no enumeration.

## 2. Results

| tool | prior | **audited** | (i) type | (ii) composes | (iii) count w/o generating |
|---|---|---|---|---|---|
| poolparty | yes | **yes** | ✓ `DnaPool` | ✓ `insertion_pools=`, `content_pool=` | ✓ `num_states` |
| dnachisel | no | **yes** | ✓ `MutationSpace` | ✓ `DnaOptimizationProblem(mutation_space=)` | ✓ `space_size` |
| pydna | no | **partial** | ✓ `Assembly` | ✗ | ✓ `get_possible_assembly_number` |
| oligopoolcalc | yes | **no** | ✗ `pd.DataFrame` | ✓ (strongly) | ✓ in substance (`Degeneracy` col) |
| seqpro | yes | **no** | ✗ `Ragged` array | ✓ (container ops only) | ✗ |
| tangermeme | partial | **no** | ✗ `torch.Tensor` | ✓ (tensor ops) | ✗ |
| ledidi | partial | **no** | ✗ `torch.Tensor` | ✗ | ✗ |
| biopython | partial | **no** | ✗ (no library type) | ✗ | ✗ |
| mpradesign | partial | **no** | ✗ `data_frame` in `list` | ✗ | ✗ |
| codongenie | partial | **no** | ✗ `list` of `dict` | ✗ | ✗ |
| mutationmaker | partial | **no** | ✗ (no Python package) | ✗ | ✗ |
| valiant | no | **no** | ✗ (CLI only) | ✗ | ✗ |
| mpranator | no | **no** | ✗ (`dict`/`list`) | ✗ | ✗ |

11 of 13 cells changed. Only `valiant` and `mpranator` are unchanged.

---

## 3. Per-tool evidence

### poolparty — yes (unchanged)

Behavioural run (read-only venv, `PYTHONDONTWRITEBYTECODE=1`,
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/.venv/bin/python`):

```
A repr: DnaPool(id=1, name='pool[1]', op='op[1]:mutagenize', num_states=48) num_states: 48
B repr: DnaPool(id=2, name='pool[2]', op='op[2]:from_iupac', num_states=64) num_states: 64
C repr: DnaPool(id=6, name='pool[6]', op='op[6]:insertion_multiscan(replace_region)', num_states=26738688) num_states: 26738688
D stack: 112
E has .to_fasta? False  .materialize? True
```

**(i)** `poolparty/pool.py:25` — `class Pool(CommonOpsMixin, ScanOpsMixin, GenericFixedOpsMixin,
StateOpsMixin, RegionOpsMixin):` with docstring *"Pool provides generic operations that work on
any sequence type."*; `poolparty/dna_pool.py:7` — `class DnaPool(Pool, DnaMixin, FilterMixin,
ExportMixin):`. Bound to a variable in the run above.

**(ii)** `poolparty/pool_mixins/scan_ops_mixin.py:519`:
`insertion_pools: Union[Pool_type, Sequence[Pool_type]],` — documented as *"Pool(s) providing
content"* (line 540–543). `poolparty/pool_mixins/region_ops_mixin.py:210`:
`content_pool: Union[Pool_type, str],` — *"content_pool : Pool or str / Pool or sequence string to
insert at the region position."* Confirmed working in line C of the run.

**(iii)** `poolparty/pool.py:134-137`:

```python
    @property
    def num_states(self) -> int:
        """Number of states for this pool."""
        return self.state.num_values
```

Docs, `poolparty/docs/pool.rst:84-86`: *"``pool.num_states`` is the total number of distinct
states (and therefore distinct sequences) the pool can produce."*
`poolparty/docs/operations/library_size.rst:4-5`: *"Every pool has a ``num_states`` property — the
number of distinct sequences it can produce."*

Counting is provably separate from generation: generation lives in a distinct op
(`poolparty/base_ops/materialize.py:15` `class MaterializeOp(Operation):`), and a pool that is
*too large to generate* still reports its size —
`pp.from_iupac("N"*20, mode="sequential")` raised
`ValueError: Number of states (1099511627776) exceeds max_num_sequential_states (1000000).`
i.e. 4^20 was computed symbolically with nothing materialised
(`poolparty/operation.py:43-46`).

**Disconfirmation.** A "no"/"partial" would have required either (a) `num_states` to be backed by
a materialised list — refuted by reading `pool.py:134` (`self.state.num_values`, delegated to
`statetracker`) and by the 1.1e12-state error above; or (b) no `Pool`-typed parameter anywhere —
refuted by `grep -rn "insertion_pools\|content_pool" --include=*.py` across the 102-file package,
which returns typed signatures in `scan_ops_mixin.py`, `region_ops_mixin.py`,
`multiscan_ops/insertion_multiscan.py`, `region_ops/replace_region.py`. I applied no softening:
PoolParty passes on the same reading that demotes seqpro and oligopool.

### dnachisel — no → **yes** (upgrade)

Source: `Edinburgh-Genome-Foundry/DnaChisel@68c0930`.

**(i)** `dnachisel/MutationSpace/MutationSpace.py:9-10`:

```python
class MutationSpace:
    """Class for mutation space (set of sequence segments and their variants).
```

Its own docstring binds it to a variable (`MutationSpace.py:24`):
`>>> space = MutationSpace([`. Exported at
`dnachisel/MutationSpace/__init__.py`: `__all__ = ["MutationSpace", "MutationChoice"]`.
Documented in the official reference — `docs/ref/core_classes.rst:45-48`:

```
MutationSpace
-----------------------

.. automodule:: dnachisel.MutationSpace
   :members:
```

`dnachisel/README.md:26`: *"**MutationSpace** is a class to represent the possible mutations at
different locations of the sequence for a given problem."*

**(ii)** `dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py:115-122` —
`def __init__(self, sequence, constraints=None, objectives=None, logger="bar",
mutation_space=None,)`; documented at lines 69-72: *"mutations_space — A MutationSpace indicating
the possible mutations. In most case the mutation space will be left to None … however some core
DNA Chisel methods will create optimization problems with a provided mutation_space to save
computing time."* Used exactly that way in-tree at
`DnaOptimizationProblem/mixins/ConstraintsSolverMixin.py:277` and
`mixins/ObjectivesMaximizerMixin.py:157` (`mutation_space=mutation_space,`). Library→library
transforms also exist: `MutationSpace.py:88-93` `def localized(self, location): … return
MutationSpace(self.choices_index[start:end], left_padding=start)` and
`MutationSpace.py:165` `def from_optimization_problem(problem, new_constraints=None)`.

**(iii)** `dnachisel/MutationSpace/MutationSpace.py:95-104`:

```python
    @property
    def space_size(self):
        """Return the number of possible mutations."""
        if len(self.multichoices) == 0:
            return 0
        choices = [len(choice.variants) for choice in self.multichoices]
        # np.prod(choices) can create overflows and warnings, so instead we use
        # the mechanism below with log/exp and a min.
        return np.exp(min(100, np.log(choices).sum()))
```

A product of per-segment variant counts: no I/O, no enumeration. Generation is a separate
method — `MutationSpace.py:132-133` `def all_variants(self, sequence): """Iterate through all
sequence variants in this mutation space."""` (an `itertools.product` generator, line 159). The
`DnaOptimizationProblem` docstring itself glosses `space_size` as the library size
(lines 76-78): *"The algorithm will use an exhaustive search when the size of the mutation space
(=the number of possible variants) is above this threshold."*

**Caveats recorded for the authors** (they do not defeat (iii) as written, but a referee should
see them): `space_size` returns a **float capped at e^100** (≈2.7e43), not an exact integer, so it
is not an exact count for very large spaces; and it returns `0`, not `1`, for a fully determined
space. Substantively, `MutationSpace` is the *search space* of a single scaffold, not a shipped
multi-member deliverable — DnaChisel's output is one optimised sequence. Under a "the object must
be the designed library you ship" reading of criterion (i), DnaChisel would be `partial`. The
rubric as written does not say that, so the cell is `yes`.

**Disconfirmation.** I looked for the evidence that would keep the prior "no": (a) that no
multi-sequence type exists — refuted by `grep -rn "^class " --include=*.py dnachisel/` which
surfaced `MutationSpace` and `MutationChoice`; (b) that `mutation_space` is private/internal —
refuted by `docs/ref/core_classes.rst` documenting it under "Core Classes", by its `__all__`, and
by the public `DnaOptimizationProblem` docstring; (c) that `space_size` enumerates — refuted by
reading its body. I also checked `grep -rli "variant librar|combinatorial|library of sequences"
--include=*.py dnachisel/` → zero hits, confirming DnaChisel never *calls* this a library, which
is why the prior pass missed it.

### pydna — no → **partial** (upgrade)

Source: `pydna-group/pydna@4e02f81`.

**(i)** `src/pydna/assembly.py:72-77`:

```python
class Assembly(object):
    """Assembly of a list of linear DNA fragments into linear or circular
    constructs. … Accepts a list of Dseqrecords (source fragments) to
    initiate an Assembly object.
```

Its own docstring binds it (`assembly.py:96`): `>>> x = Assembly((a,b,c), limit=14)` and shows it
denoting *multiple* products: `>>> x.assemble_circular()` → `[Contig(o59), Contig(o59)]`.
`README.md:229` binds it too — `asm = Assembly(fragments, limit=10)` — with the comment
(`README.md:231`) *"# From the assembly object, which can generate all possible products, get a
circular"*. `Assembly` is in the top-level namespace: `src/pydna/all.py` `__all__` includes
`"Assembly"`.

**(iii)** `src/pydna/assembly2.py:1550-1562`:

```python
    def get_possible_assembly_number(self, paths: list[list[int]]) -> int:
        """
        Get the number of possible assemblies from a list of node paths. Basically, for each path
        passed as a list of integers / nodes, we calculate the number of paths possible connecting
        the nodes in that order, given the graph (all the edges connecting them).
        """
```

Called *before* any sequence is built — `assembly2.py:1473-1477`:
`possible_assemblies = self.get_possible_assembly_number(unique_linear_paths)` then
`raise ValueError(f"Too many assemblies ({possible_assemblies} pre-validation) to assemble")`.
Same pattern in `_validate_max_assemblies` (`assembly2.py:1574-1581`). So the product count is
obtainable without generating products and without files. Caveat: this method is on
`assembly2.Assembly` (`assembly2.py:1168`), not on the `pydna.all`-exported
`assembly.Assembly` (`assembly.py:72`), whose only methods are `assemble_linear`,
`assemble_circular`, `__repr__`.

**(ii) FAILS.** No operation in pydna accepts an `Assembly` instance.
`grep -rn ": *Assembly\|:Assembly\|assembly: " --include=*.py src/pydna/` returns only
`assembly: EdgeRepresentationAssembly` parameters (a *tuple of graph edges*, not the object) and
`AssemblyAlgorithmType`. `PCRAssembly`/`SingleFragmentAssembly` subclass `Assembly`; subclassing
is not passing. The things that flow onward are the *members*
(`assemble_circular() -> list[Dseqrecord]`), not the library.

**Disconfirmation.** To keep the prior "no" I needed either no bindable multi-product design type
(refuted by `assembly.py:72` + README:229) or no pre-generation count (refuted by
`assembly2.py:1550`). To promote to "yes" I needed an `Assembly`-typed parameter somewhere; the
grep above found none, and `src/pydna/design.py:assembly_fragments` takes a plain `list`. Hence
`partial`, not `yes`.

### oligopoolcalc — yes → **no** (demotion; the rubric's DataFrame exclusion is decisive)

Source: `ayaanhossain/oligopool@b88fa39` (v`2026.02.22.1`).

**(i) FAILS.** There is no library class. `grep -rn "^class " --include=*.py oligopool/` returns
exactly: `Scry`, `vectorDB`, `OligopoolParser(argparse.ArgumentParser)`,
`OligopoolFormatter(argparse.RawTextHelpFormatter)`, `SafeCounter`, `SafeQueue`,
`FailureSampler` — no library type. Every design op is typed to return a bare DataFrame, e.g.
`oligopool/barcode.py:17-34`:

```python
def barcode(
    input_data:str|pd.DataFrame,
    …
    verbose:bool=True) -> Tuple[pd.DataFrame, dict]:
```

`docs/docs.md:13`: *"Most modules take CSV/DataFrames in and return `(out_df, stats)`."*
The rubric's NO clause names `DataFrame` explicitly, so this gates the cell to `no`.

**(ii) PASSES, strongly** — and this is why the cell is contentious. `docs/docs.md:1915`:
*"Because modules accept and return DataFrames, small `pandas` transforms fit naturally between
steps."* `docs/docs.md:140-152`: *"Most modules return both a DataFrame and a stats dictionary:
`out_df, stats = op.barcode(input_data=df, ...)` … **Chainable**: Output of one module feeds into
the next"*. A fully in-memory, file-free chain is a documented example (`docs/docs.md:199-232`):
`df = pd.DataFrame({...})` → `df, _ = op.barcode(input_data=df, …)` → `pd.concat` → `df, _ =
op.barcode(input_data=df, …)`. `output_file` is `None` by default in library usage
(`barcode.py:48-49`).

**(iii) PASSES in substance.** `oligopool/base/core_degenerate.py:860-864`:

```python
def get_expansion_count(sequence):
    '''
    Compute degeneracy (number of concrete
    sequences) without expanding.
```

(marked *"Internal use only."*, line 864). It is reachable through two public paths with no file
written: `compress()` returns a `synthesis_df` whose documented columns are
*"`DegenerateID`, `DegenerateSeq`, `Degeneracy`, `OligoLength`"* (`compress.py:42-43`) — summing
`Degeneracy` gives the library size without expanding; and `expand()`'s `expansion_limit`
pre-check returns `'estimated_expansion': int(total_estimated)` in its stats
(`expand.py:159-171`), computed by `get_parsed_degenerate_info` (`core_degenerate.py:1106-1108`)
before any expansion. Documented: *"`expansion_limit` … Safety cap for maximum total expanded
sequences; if estimated expansion exceeds this limit, expansion fails"* (`expand.py:32-33`).

**Disconfirmation.** To keep "yes" I searched for a dedicated library type
(`grep -rn "^class "`, `oligopool/__init__.py` `__api__` list — 22 entries, all functions plus
`vectorDB`/`Scry`) and for a `.num_variants`-style accessor on a returned object
(`grep -rn "def get_expansion|def .*degeneracy"` → only module-level helpers). Neither exists:
the library really is a `pd.DataFrame`. To demote all the way to the VaLiAnT anchor I would need
oligopool to be file-only, which it demonstrably is not (see (ii)). **This cell is the strongest
single piece of evidence that the row is mis-specified — see §4.**

### seqpro — yes → **no** (demotion; the prior "yes" rested on the ragged container)

Source: `ML4GLand/SeqPro@63a8439`.

**(i) FAILS on the rubric's named exclusion.**
`python/seqpro/rag/_core.py:85-86`:

```python
class Ragged(NDArrayOperatorsMixin, Generic[RDTYPE_co]):
    """A non-branching ragged array with a single ragged axis (Spec A)."""
```

That is a general-purpose array container, not a library type. seqpro has **no library-design
concept at all**: `grep -rni "librar|mutagen|saturation|variant|oligo|pool" --include=*.py
python/seqpro/` returns only RNA-seq *library-size* normalisation in `transforms/tmm.py` and
`TypeVar` names — zero library-design hits. `python/seqpro/__init__.py` `__all__` is
`gc_content, length, nucleotide_content, decode_ohe, decode_tokens, ohe, pad_seqs, tokenize,
bin_coverage, jitter, k_shuffle, random_seqs, reverse_complement, cast_seqs, …` — all sequence
*processing*. Its one generator returns a plain array:
`_modifiers.py:281-285` `def random_seqs(shape…, alphabet…, seed…) -> NDArray[np.bytes_]:`.
README: *"All functions in SeqPro take as input a string, a list of strings, a NumPy array of
strings, a NumPy array of single character bytes (`S1`), or a NumPy array of one-hot encoded
sequences (`uint8`)."*

**(ii) technically passes, for container ops only** — `rag/_ops.py:22-23`
`def reverse_complement(rag: Ragged[np.bytes_],`; `_ops.py:341`
`def concatenate(rags: Any, axis: int) -> "Ragged[Any]":`; `_ops.py:262-263`
`def to_padded(rag: Ragged[Any],`. None of these is a library-design operation.

**(iii) FAILS.** `rag/_core.py:263-268`:

```python
    def __len__(self) -> int:
        """Return the size of the outermost dimension (shape[0])."""
        s = self._layout.shape[0]
```

and `_core.py:287-289` `def shape(self) -> tuple[int | None, ...]: return self._layout.shape`.
Both read a layout that only exists once the data array does — `Ragged.__init__` requires
`data`, and `from_offsets`/`from_lengths` both take a materialised `data: NDArray[Any]`.

**Disconfirmation.** For "yes" I would need a design-semantics type or a pre-generation count.
I grepped `^class ` across `python/seqpro/` (only `Ragged`, layouts, alphabets, sklearn-style
transforms), grepped for `space_size|num_states|n_variants|__len__` in `rag/`, and read
`docs_api_ragged.md`/`docs_ragged.md` notes present in the working copy. Nothing. The prior memo's
own words (`final/seqpro.md:94`, a lead only) — "A *batch* of sequences is genuinely a first-class
object" — describe the batch container, which the rubric excludes by name.

### tangermeme — partial → **no**

Source: `jmschrei/tangermeme@8d732b8`.

Every multi-sequence value is a bare tensor. `tangermeme/ersatz.py:20-25`
`def insert(…) -> torch.Tensor:`; `:97-103 substitute → torch.Tensor`;
`:198-205 multisubstitute → torch.Tensor`; `:300 delete(X: torch.Tensor, …)`;
`:343-350 randomize → torch.Tensor`; `:425-431 shuffle → torch.Tensor`;
`:608-615 dinucleotide_shuffle → torch.Tensor`. Doc string of the parameter:
`product.py:233-234` *"X: torch.tensor, shape=(-1, len(alphabet), length) — A one-hot encoded set
of sequences to make predictions for."*

`grep -rn "^class " --include=*.py tangermeme/` returns only `SaturationMutagenesisRawResult`,
`TangermemeWarning`, `PerturbationResult`, `PerturbationAnnotationsResult`,
`AttributionReferencesResult`, `SpaceResult` — all `NamedTuple` result bundles holding
*prediction* tensors (`results.py:24-26`: `y_before: torch.Tensor | list[torch.Tensor]` /
`y_after: …`), not sequence libraries. **(i) fails.** Cartesian products are built and consumed
inside a call — `product.py:197-208 apply_product(func, model, X, args, …) -> torch.Tensor |
list[torch.Tensor]`, docstring *"Apply a function on the cartesian product between X and each
args"* — never handed back as an object. **(iii) fails:**
`grep -rn "def n_\|def size\|def __len__\|space_size\|num_variants" --include=*.py tangermeme/`
→ **zero hits**.

**Disconfirmation.** I specifically opened the three most promising modules (`space.py`,
`product.py`, `design/`) looking for a library-space object; `design/` is four plain functions
(`screen`, `beam_substitution`, `greedy_marginalize`, `greedy_substitution`) and `space.py` has
only `SpaceResult`. **(ii)** does pass in the weak sense that tensors flow between ersatz ops,
but with (i) and (iii) both failing the rubric's generic-structure clause gates this to `no`.

### ledidi — partial → **no**

Source: `jmschrei/ledidi@adbca70`.

The designed library is a tensor, and the README says so — `README.md:82`: *"The returned `X_hat`
has shape `(batch_size, n_channels, length)`: a batch of independently sampled designed
sequences (`batch_size` defaults to 16)."* `ledidi/ledidi.py:508` `def fit_transform(self, X,
y_bar):` returns `best_sequence` / `(best_sequence, history)` (`ledidi.py:632-633`). Even the
"catalog" mode stays a tensor — `ledidi.py:116-119`: *"one can design an affinity catalog by
passing in a list of target values in `y_bar` … an additional dimension is added to the front of
the returned tensor of designed sequences."*

The one bindable object is explicitly disclaimed as user-facing — `README.md:118`: *"Note that
there is also a `Ledidi` object which is an internal wrapper that you typically won't need to use
directly."* It is a designer (`class Ledidi(torch.nn.Module)`, `ledidi.py:265`), not a library.
**(i) fails.** **(iii) fails:** `grep -rn "def __len__\|def size\|num_seq" --include=*.py ledidi/`
→ **zero hits**. **(ii) fails:** the only other classes are `MinGap` (`losses.py:23`) and
`DesignWrapper` (`wrappers.py:23`), and `DesignWrapper.__init__` takes `models`, raising
`TypeError` for anything that is not a list/tuple of `torch.nn.Module` (`wrappers.py:48-62`) —
nothing accepts a `Ledidi` or a library.

**Disconfirmation.** For "partial" I needed the `Ledidi` object to qualify as the library object;
the README disclaims it and it has no member-count accessor. `greedy_pruning(model, X, X_hat,
…)` (`pruning.py:25`) takes tensors, not a designer.

### biopython — partial → **no**

Source: `biopython/biopython@c948960` (current package, per instruction).

Biopython has **no library-design layer at all**, so there is no library object to hold.
Disconfirmation greps, all over `Bio/`:

- `grep -rli "combinatorial librar|variant librar|oligo pool|saturation mutagenesis|degenerate librar" --include=*.py Bio/` → **zero files**.
- `grep -rn "def space_size|def num_states|def n_variants|def num_variants|def library_size|def n_members" --include=*.py Bio/` → **zero hits**.
- `grep -rn "^class .*(Library|Pool|Design|Variants|Space)" --include=*.py Bio/` → one hit, `Bio/Graphics/BasicChromosome.py:806 class SpacerSegment(ChromosomeSegment):` (a drawing primitive).
- `grep -rn "itertools.product" --include=*.py Bio/` → one hit, `Bio/PDB/mmtf/mmtfio.py:79`, generating PDB chain-ID strings. No sequence-library enumeration anywhere.

The nearest multi-sequence types are containers, which the rubric excludes by name.
`Bio/Align/__init__.py:59-65`:

```python
class MultipleSeqAlignment:
    """Represents a classical multiple sequence alignment (MSA).

    By this we mean a collection of sequences (usually shown as rows) …
```

and lines 82-87 *"In some respects you can treat these objects as lists of SeqRecord objects …
`>>> len(align)` `7`"* — `len()` over materialised records. **(i) and (iii) fail; (ii) has no
library-design operation to pass into.**

Note for the authors: this `no` is a statement about the *absence of a library-design layer*, not
about Biopython's API quality — `Seq`, `SeqRecord`, `Bio.Restriction.Analysis`, and lazy
`SeqIO.index` are all genuinely first-class; none of them is a designed multi-sequence library.
Under the row as written that distinction is invisible, which is part of the case in §4.

### mpradesign — partial → **no**

Sources: `andrewGhazi/mpradesigntools@afd386e`, `andrewGhazi/designMPRA@0cf56ee`.

`NAMESPACE` exports exactly two functions: `export(processVCF)`, `export(spread_and_fix_indels)`.
The documented return value is a generic structure — `man/processVCF.Rd`, `\value{}`:

> *"A list of two data_frames. The first, named 'result', is a data_frame containing the labeled
> MPRA sequences. The second, named 'failed', is a data_frame listing the SNPs that are not able
> to have MPRA sequences generated and the reason why."*

Confirmed in source: `R/processVCFfast.R:1341` `return(res)` where `res$result` is a tibble
(line 1335 `res$result %>% dplyr::mutate(n_bp = nchar(sequence))`). **(i) fails.**
**(ii) fails:** the only other export is `spread_and_fix_indels = function(vcf_path)`
(`R/processVCFfast.R:1418`) — it takes a *file path*, not the result. **(iii) fails:** no count
accessor; the member count is only knowable from the user's own `nper` arithmetic.
The Shiny front end is file-output-only: `designMPRA/server.R:161-164`
`output$downloadSequences = downloadHandler(… write_tsv(vcfOut()$result, file))`.

**Disconfirmation.** For "partial" I needed an R class carrying library semantics:
`grep -rn "setClass|R6Class|setRefClass|structure(class" R/` → **zero hits**. The package defines
no classes at all.

### codongenie — partial → **no**

Source: `synbiochem/CodonGenie@c26439f` (repo `neilswainston/CodonGenie` after redirect).

`codon_genie/codon_selector.py:59-60` `class CodonSelector(): '''Class to optimise codon
selection.'''` is a *selector*, not a library, and its outputs are plain Python collections:
`:81 return _format_results(results)` where `_format_results` (`:184-187`) is
`return sorted([codon for result in results for codon in result], …)` — a list of dicts. Each dict
carries `'ambiguous_codon_expansion': tuple(codons)` (`:126`), i.e. the expansion is **already
materialised** by `codons = [''.join(c) for c in itertools.product(*ambig_codon_nucls)]`
(`:110`). **(i) and (iii) fail.** **(ii) fails:** the only other public entry points are
`analyse_codon` (`:83`) and the REST layer in `client.py`/`main.py`; nothing accepts the result.
Scope note: CodonGenie designs **one degenerate codon** at a time — there is no multi-sequence
library in the tool at all. `README.md` documents only `bash start_server.sh 80`.

**Disconfirmation.** `find … -name "*.py"` over the whole repo returns 13 files; I read
`codon_selector.py` in full and grepped for any class besides `CodonSelector` — none in the
design path. Not `unknown`: the source is fully available, so the cell is scored on source.

### mutationmaker — partial → **no**

Source: `ra100/Mutation_Maker@396c1c0` (the live author-maintained fork).

There is **no installable Python package and no documented Python API**:
`find … -maxdepth 2 -name setup.py -o -name pyproject.toml -o -name setup.cfg` → **zero hits**.
`README.md:157` installs requirements only (`pip install -r backend/requirements.txt`). The user
interface is HTTP — `README.md:47-48` *"Webserver: http://localhost:3000 / API:
http://localhost:8000"* — served by hug + Celery: `api/api.py:41,49,57`
`def find_ssm_primers(body, response, hug_celery)`, `find_qclm_primers`, `find_pas_primers`, with
downloads via `def export_task(task_id, hug_celery)` (`api/api.py:133,142,151`).

Internal backend classes exist (`backend/mutation_maker/pas_solution.py:271 class PASSolution:`
*"A (partial) PAS workflow solution"*; `qclm_solution.py:153 class QCLMSolution:`;
`generate_oligos.py:196 class OligoGenerator(object):`) but they are not a published API surface,
and the oligo set is emitted as a plain `List` — `pas_output.py:90-96
def combine_oligos_list(…) -> List:`, `:252 … -> [PASResult]`. **(i)/(ii)/(iii) all fail** for a
user of the tool.

**Disconfirmation.** For "partial" I needed a documented importable object. I checked for packaging
files (none), read `README.md` install/usage sections (Docker + REST only), listed `docs/` (only
`docs/plans`), and inspected the API layer. `PASSolution.add_fragment(fragment: PASFragment)`
(`pas_solution.py:306`) does compose internally, but nothing in that path is reachable as a
documented library.

### valiant — no (unchanged; anchor confirmed)

Source: `cancerit/VaLiAnT@8796cc1` (develop, v4.0.0).

CLI-only: `pyproject.toml:23-24` `[project.scripts]` / `valiant = "valiant.__main__:main"`.
`src/valiant/__init__.py` contains the licence header and `__version__ = '4.0.0'` — **nothing
else is exported**. Output is written field-by-field to a file handle:
`src/valiant/meta_table.py:85-88`

```python
def _write_field(fh, s: str | None) -> None:
    if s:
        fh.write(s)
    fh.write(',')
```

with `write_meta_record` (`:175`) and `write_no_op_meta_record` (`:118`). The intermediate store is
a discarded in-memory database: `src/valiant/db.py:132` `with sqlite3.connect(':memory:') as
conn:`. README frames everything as input files → *"Output files: … oligonucleotide metadata file
(CSV), variant file (VCF, SGE-only)"* (`README.md:96-100`).

**Disconfirmation.** The only collection-like class in the package is
`src/valiant/mutator.py:77 class MutatorCollection:` — a collection of *mutator types*, not
sequences. Confirmed by `grep -rn "^class .*Librar|^class .*Pool|^class .*Collection|^class .*Set\b" --include=*.py src/`.

### mpranator — no (unchanged)

Source: `hemberg-lab/MPRAnator@9969790`.

Module-level functions over plain builtins, inside a Django app. Returns are dicts, lists and
HTTP responses: `part1.py:187 return {"storeID": storeID, "Sequences": Sequences}`,
`:196 return Synthetics`, `:292 return response, sequenceHTMLL`, `:319 return barCodes,
numOfBarCodesPerSequence`. `grep -n "^class " myfunctions.py oligo.py part1.py part3.py` → **no
classes**. The distributable script only POSTs to the web service and writes the reply to disk:
`downloadables/MpraMotifs_script.py:72-73`

```python
openfile = open(fileToWriteResults,"w")
openfile.write(result.content)
```

`README.md` is a single line, `# MpraNator`. **(i)/(ii)/(iii) all fail.**

---

## 4. Why the row must be reworded

### 4.1 The decisive finding: criteria (ii) and (iii) are other rows of this same table

The rubric I was handed is three criteria wide. Two of the three are already separate rows in the
authors' own row set. From `revision/tool_survey/ROWS_v2.md:22-25`:

| Key | Question |
|---|---|
| `library_first_class_object` | Is the library an object the user can hold, inspect, transform and pass onward — as opposed to a tool that only writes files? |
| `composable_operations` | Do design steps compose/chain/nest arbitrarily, or is the pipeline fixed by the tool? |
| `lazy_generation` | Are sequences produced on demand rather than fully materialized? |
| `library_algebra` | Can whole libraries be combined/sampled/repeated as operations (stack, concat, sample), inside the tool? |

- Rubric criterion **(iii)** ("size queryable without generating all members") *is* `lazy_generation`.
- Rubric criterion **(ii)** ("the object can be passed into another library-design operation") *is*
  `library_algebra` (and overlaps `composable_operations`).

Worse, the authors did not merely happen to have those rows — they created them by splitting this
exact row, and they wrote down why. `ROWS_v2.md:13`:

> **SPLIT** | `library_as_object` -> A1/A4 | 8 of 12 scored "partial" — the row was too coarse to
> discriminate. VaLiAnT's own paper says the final concatenated library is assembled by the user
> outside the tool, so **"can hold a library object" and "can combine libraries" are distinct
> claims.**

The rubric under audit re-merges precisely the claims that sentence separates. Two consequences,
both fatal to the row as written:

1. **Triple counting.** PoolParty wins the same architectural advantage in
   `library_first_class_object`, `lazy_generation`, and `library_algebra`. A referee who notices
   this reads the matrix as padded, which damages the paper more than a lost cell would.
2. **Double penalty.** oligopool loses `library_first_class_object` for a lazy-generation deficit
   that `lazy_generation` already scores, and biopython/seqpro lose it for an algebra deficit that
   `library_algebra` already scores. Competitors are marked down twice for one shortcoming.

This is not a judgement call about thresholds; it is a structural defect, checkable against
`ROWS_v2.md` in ten seconds.

### 4.2 The two false equivalences the merged rubric produces

Applying one consistent threshold produces two false equivalences that a referee will find
immediately, and they point in *opposite* directions — the signature of a compound row.

**False equivalence A: `oligopoolcalc == valiant`.** The row's own headline question is *"a
first-class object the user can hold, inspect, transform and pass onward — as opposed to a tool
that only writes files."* oligopool answers that question affirmatively in its own documentation
(`docs/docs.md:140-152`, `:199-232`, `:1915`): a fully in-memory chain through 12+ design
operations with `output_file=None`. VaLiAnT has no Python API whatsoever
(`src/valiant/__init__.py` exports only `__version__`). The rubric nonetheless collapses them to
the same value, purely because oligopool's in-memory value has the Python type `pd.DataFrame`.
That is a claim about type nominalism, not about whether the library can be held and passed
onward. Understating a competitor this way is precisely the failure mode the audit brief warns
against.

**False equivalence B: `dnachisel == poolparty`.** DnaChisel satisfies all three criteria
literally — `MutationSpace` (purpose-built, documented in `docs/ref/core_classes.rst`),
`DnaOptimizationProblem(mutation_space=…)` (documented parameter), `space_size` (combinatorial,
no I/O). But DnaChisel is a single-sequence optimiser: its deliverable is one sequence, and
`MutationSpace` is an internal search space that is never shipped. Awarding it the same value as
PoolParty removes the row's power to support the paper's claim and overstates DnaChisel.

Both failures trace to the same cause: criterion (i) is a **type** property, criterion (ii) is a
**composability** property, criterion (iii) is a **lazy-evaluation** property. They are
independent, they each already have a home in the row set, and compressing them into one
yes/partial/no forces arbitrary tie-breaks. Note also that criterion (ii)'s phrase "another
LIBRARY-DESIGN operation" is undefined for tools that have no operation anyone would call library
design (DnaChisel, pydna, biopython) — that is where my own judgement had to do the most work, and
it is a bad sign for reproducibility across raters.

### 4.3 Proposed reword

Do **not** add rows. Restore this row to the authors' own one-line definition and delete criteria
(ii) and (iii), which `library_algebra` and `lazy_generation` already carry:

> **`library_first_class_object`** — Does a documented API call return a value that represents the
> whole multi-sequence library as an instance of a **type the tool defines for that purpose**
> (carrying member identity, design provenance, or membership operations), which the user can bind
> to a variable and inspect — as opposed to the tool only writing files, or handing back a
> general-purpose container (ndarray, DataFrame, ragged array, torch tensor, list/dict of
> sequences) whose library meaning lives only in the caller's head?
>
> `yes` = such a type exists and is documented. `partial` = a durable in-memory library value
> exists and is documented, but its type is general-purpose. `no` = the library exists only as
> files, or only as a transient return value with no library identity.

Scores under the reworded row (same primary sources, same evidence as §3; only the threshold
changes):

| tool | reworded value | basis |
|---|---|---|
| poolparty | **yes** | `DnaPool` / `Pool` (`pool.py:25`, `dna_pool.py:7`) |
| dnachisel | **yes** | `MutationSpace` (`MutationSpace.py:9`), documented in `docs/ref/core_classes.rst:45` |
| pydna | **yes** | `Assembly` (`assembly.py:72`), in `pydna.all.__all__` |
| oligopoolcalc | **partial** | durable, documented, file-free in-memory library — but it is a `pd.DataFrame` |
| seqpro | **no** | `Ragged` is "a non-branching ragged array"; no library concept in the package |
| tangermeme | **no** | `torch.Tensor` only |
| ledidi | **no** | `torch.Tensor`; `Ledidi` object disclaimed as internal in README |
| biopython | **no** | no library-design layer; `MultipleSeqAlignment` is "a collection of sequences" |
| mpradesign | **no** | `list` of two `data_frame`s |
| codongenie | **no** | `list` of `dict`; single codon, no multi-sequence library |
| mutationmaker | **no** | no installable/documented Python API; REST + JSON |
| valiant | **no** | CLI only; `__init__.py` exports `__version__` |
| mpranator | **no** | plain `dict`/`list`; no classes in the design path |

Why this is better for the paper, not worse:

- It is the claim the authors themselves wrote, so it needs no defending.
- It ends the double-counting: PoolParty's distinctiveness then comes from being `yes` on
  `library_first_class_object` **and** `library_algebra` **and** `lazy_generation` — three
  independent facts, which is far more persuasive than one compound cell, and it is a claim
  DnaChisel and pydna cannot match (both fail `library_algebra`; both fail or barely reach
  `lazy_generation`).
- It stops understating oligopool, whose in-memory DataFrame pipeline is real and documented. A
  `partial` here is defensible to its author; a `no` is not.
- It stops overstating DnaChisel: `yes` on this row is honest (the type exists and is documented),
  while its failure on the other two rows correctly records that it is a single-sequence optimiser.

If the authors keep the row exactly as the rubric specifies, the cells in §2 are the honest
scoring, but §4.1 and §4.2 should be resolved first — the `ROWS_v2.md:13` contradiction is
discoverable by any referee who reads the authors' own methods.

---

## 5. Primary-source working copies used

Fetched via `curl` from the canonical repositories (read-only; nothing installed). Commit-pinned
tarball roots under
`/tmp/claude-1000/-mnt-c-Users-zhliu-Desktop-KinneyLab/13e50e09-5cee-425a-a895-0f176887e7a6/scratchpad/src2/`:

| tool | root | ref |
|---|---|---|
| oligopoolcalc | `ayaanhossain-oligopool-b88fa39/` | master @ `b88fa39` |
| dnachisel | `Edinburgh-Genome-Foundry-DnaChisel-68c0930/` | master @ `68c0930` |
| pydna | `pydna-group-pydna-4e02f81/` | master @ `4e02f81` |
| ledidi | `jmschrei-ledidi-adbca70/` | master @ `adbca70` |
| tangermeme | `jmschrei-tangermeme-8d732b8/` | main @ `8d732b8` |
| biopython | `biopython-biopython-c948960/` | master @ `c948960` |
| valiant | `cancerit-VaLiAnT-8796cc1/` | develop @ `8796cc1` |
| codongenie | `neilswainston-CodonGenie-c26439f/` | master @ `c26439f` |
| mutationmaker | `ra100-Mutation_Maker-396c1c0/` | master @ `396c1c0` |
| mpranator | `hemberg-lab-MPRAnator-9969790/` | master @ `9969790` |
| mpradesign | `andrewGhazi-mpradesigntools-afd386e/`, `andrewGhazi-designMPRA-0cf56ee/` | master |
| seqpro | `../src/ML4GLand-SeqPro-63a8439/` | main @ `63a8439` |

poolparty was read in place from
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/` (strictly read-only; the
only writes in this audit are this file and its directory).
