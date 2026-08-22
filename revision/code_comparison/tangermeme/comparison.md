# The SpliceAI GT arm, in PoolParty and in tangermeme

**Both tools produce the identical 200,000 sequences — same set, same order —
verified against the published design cards.**

This is **not** a conciseness comparison. Three statements each side. What differs
is the shape of what comes back, and both shapes are right for their tool's purpose.

## Scope

The **GT arm** of the manuscript's Example 3: 2,000 cryptic 5′ splice sites of
graded MaxEntScan strength, substituted at each of 100 positions across a 201 bp
region carrying the canonical 5′ splice site of SMN2 exon 7. 2,000 × 100 = 200,000.

Not the whole example. The full Example 3 also builds the 9-mer ladder from a PWM
(`from_motif` → `filter` → `score` → quantile sampling), slices the region and two
5 kb flanks from GRCh38, and adds a matched **GA control arm** in which the cryptic
GT is disrupted — 400,000 sequences in total, with a seven-column card table. Only
two of those columns come from operations; the rest are derived in analysis.

`SCAN_POSITIONS = list(range(51, 90)) + list(range(107, 168))` — positions 90–106
are omitted because the canonical site sits at 100.

Inputs are committed: `revision/benchmarks/spliceai_background_201bp.txt` and
`spliceai_cryptic_seqs.csv`, whose 2,000 nine-mers match the published GT set
exactly. Ground truth is `poolparty/examples/spliceai_design_cards_published.csv.gz`
filtered to `condition == "GT"`. So this comparison is re-runnable from the repo
with no genome slice.

## PoolParty

```python
cryptic_pool = pp.from_seqs(cryptic_sites, mode="sequential",
                            cards={"seq": "cryptic_sequence"})
library = pp.from_seq(target_region).replacement_scan(
    cryptic_pool, positions=SCAN_POSITIONS, mode="sequential",
    cards={"start": "cryptic_position"})
df = library.to_df(num_cycles=1)
```

`num_states` is 200,000 before anything is generated. The card columns are named by
the operations that produce them — the pool contributes the 9-mer, the scan
contributes the position.

## tangermeme

```python
target_ohe  = one_hot_encode(target_region).unsqueeze(0).expand(len(cryptic_sites), -1, -1)
cryptic_ohe = torch.stack([one_hot_encode(site) for site in cryptic_sites])
library_ohe = torch.stack([substitute(target_ohe, cryptic_ohe, start=p)
                           for p in SCAN_POSITIONS])       # (100, 2000, 4, 201)
```

`expand` rather than `repeat` because `substitute` clones internally
(`ersatz.py`: `X = torch.clone(X)`), so the caller needs no copy. One batched
substitution per position; the motif tensor pairs one-to-one with the sequence
batch. `substitute` rather than `insert` because `replacement_scan` preserves
length — `insert` would lengthen the sequence.

## What differs

| | PoolParty | tangermeme |
|---|---|---|
| Returns | 200,000-row table: `seq`, `cryptic_sequence`, `cryptic_position` | `(100, 2000, 4, 201)` one-hot tensor |
| Generated | lazily; `num_states` known before generation | eagerly; **161 MB** materialised |
| Row order | site-major — each 9-mer's 100 positions contiguous | position-major, the natural loop order |
| Provenance | columns on the output | the axis layout |

Neither shape is better. tangermeme's is what you want when the next call is a
model, which is what the package is for. PoolParty's is what you want when the next
step is analysis with the design parameters as covariates — which is what Example 3
does, using cryptic-site strength and position directly.

## Verification, and what the script does beyond the panels

`repro.py poolparty` and `repro.py tangermeme` each compare independently against
the published cards, so matching both implies matching each other. Both report
identical set and identical order, 0 rows unique to either side.

Three steps in the script are **comparison scaffolding, not how tangermeme is
used**, and they are marked as such:

- `transpose(0, 1).flatten(0, 1)` — reorders position-major to site-major. Without
  it the set still matches but the order does not. This is a real trap: the natural
  tangermeme loop and PoolParty's enumeration nest the two axes oppositely.
- `[characters(seq) for seq in flat]` — decodes one-hot to strings, one sequence per
  call. `characters` accepts a 3-D input only when its first axis is 1, so there is
  no batch decoder.
- the `cards` DataFrame — a tangermeme user feeding a model would index the tensor
  instead.

The two tools need separate environments (`torch` versus `poolparty`), so the
script takes the side as an argument rather than importing both.

## Corrections made while building this

Recorded because each was wrong in a way that would have produced a misleading
figure.

| Claim | Correction |
|---|---|
| Tutorial B7 "Cartesian Product" is a good comparison target | No. `apply_product` crosses sequences against **model side-inputs** — cell state, read depth — not motifs against positions |
| Tutorial B3 "Spacing" is a good target | No. `space()` requires a model and returns predictions; the sequences it builds are local and discarded |
| tangermeme cannot enumerate | It can, in `space`, `saturation_mutagenesis` and all four `design/` functions. It never **returns** the enumeration — `design/` returns the argmax, the others return predictions |
| PoolParty's lazy generation distinguishes it | Not from tangermeme. Its B7 documentation gives the same reason for not materialising: *"only the model predictions are stored in CPU memory as opposed to the (usually much larger) inputs"* |
| The comparison shows PoolParty is more concise | No. Three statements each |
| ~640 MB materialised | **161 MB.** `one_hot_encode` returns `int8`, not `float32` |
| The first tangermeme snippet drafted | A comparison harness, not idiomatic usage — it re-encoded inside the loop and maintained a side table a user would get from the axis layout |

## Why there is no second panel

A model-guided-design panel was considered, to concede that `beam_substitution`
optimises against a model where PoolParty only accepts a scoring callable
(`model_guided_variants`: tangermeme ●, PoolParty ◐). With the comparison above
already even-handed — three statements each, neither shape privileged — the
concession panel is not needed to balance it. The concession itself stands in
Table 2.

## Source review

Three independent agents inventoried the package (20 modules, including the
`design/` subpackage), read the documentation, and checked specifically for any
function that enumerates and returns sequences. All three concluded the same:

- **No public function enumerates positions, spacings or motif combinations and
  returns the resulting sequences.** Every `ersatz` primitive takes a single scalar
  `start`.
- The motif axis is broadcast (1→N) or paired (N→N), never N×M. Their own test
  asserts the failure: `assert_raises(RuntimeError, substitute, X_batch, motif, 10)`.
- `design/_substitute.py::_fast_tile_substitute` is the closest thing to a
  combinatorial generator — private, single-motif, mutates in place, returns `None`.
- No per-sequence metadata for generated sequences anywhere. `example_idx` is
  provenance for *extracted* windows. `io.one_hot_to_fasta(..., headers=)` is the
  only per-sequence naming mechanism and is absent from the API docs and every
  tutorial.
- Only 2 of 14 tutorials return sequences at all: A1 (`ersatz`) and B6 (`design`).

None of that is a deficiency. The README says the package *"aims to be as low-level
and simple as possible"*, and the docs add that the functions *"can easily be built
on top of if you'd like to customize your analyses."* Library construction is
stated as the user's job.
