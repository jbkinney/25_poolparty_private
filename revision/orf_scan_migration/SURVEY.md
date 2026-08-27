# ORF-aware indel scans: v1 → v2 migration survey

**Question.** PoolParty v1 shipped `DeletionScanORFPool` and `InsertionScanORFPool`,
codon-level frame-preserving indel scans. v2's `deletion_scan` / `insertion_scan` have
no ORF vocabulary at all. Should the v1 capability be migrated?

**Answer.** **No port. Add a thin ergonomics layer, and fix B3 in the same change.**

The capability is not missing. All ten v1 behaviours tested here reproduce
**byte-identically** in v2 today with no code change. But the recipe recorded in
`FINDINGS.md` B1 is **not safe as written**, and v2 currently holds the reading frame as
validated, immutable data and then silently ignores it — for frames +2 and +3 the
recorded idiom deletes six triplets, **not one of which is a codon**, with no warning.
That is a correctness hazard that documentation alone cannot remove.

**Scope.** Survey and design only. Nothing in either PoolParty repo was modified.

---

## 0. Contradictions with the briefing

Everything I was handed was re-derived by execution. Four items in the briefing need
correcting, one claim raised against this survey is half right, and one
new hazard was found that had not been reported at all.

| # | Claim as briefed | Verdict | Where |
|---|---|---|---|
| 1 | Bare `positions=slice(frame-1, None, 3)` expresses in-frame codon deletion | **Correct but unsafe.** Silently over-runs the ORF 3′ end whenever a downstream flank exists. `region=` is required, not optional | [§2.2](#22-the-recorded-idiom-is-unsafe-without-region) |
| 2 | `insertion_multiscan` exposes `names=`; precedent for the B3 fix | **Understated.** `deletion_multiscan` exposes it too — the precedent is exact, not analogous | [§5](#5-b3-the-_del-collision) |
| 3 | ORF-awareness is "a property of exactly three operations … per `OrfRegion`'s own docstring" | **Count right, attribution wrong.** The docstring names two (`stylize_orf`, `mutagenize_orf`); the code has three (`translate` as well) | [§3.1](#31-which-operations-actually-read-the-frame) |
| 4 | `write_tags` defaults False on `generate_library` | **Wrong function.** `generate_library` has no such parameter; `write_tags=False` is on the export mixin (`to_df`/`to_file`). Tags **do** leak through `generate_library` | [§3.2](#32-region-does-constrain-but-does-not-frame) |
| 5 | *(raised against this survey by `kinneylab-bb`)* `_resolve_frame` is three distinct implementations, one carrying a circular-import workaround | **Half right.** My "character-for-character identical" was wrong and is corrected. But it is two AST-identical copies plus one with a **vestigial** local import — `party.py` imports nothing from `orf_ops`, and hoisting **both** sites (`stylize_orf.py:22` and `:142`) leaves the suite at 2915 passed / 14 xfailed, unchanged. Raised, tested, and withdrawn by `kinneylab-bb` on independent verification | [§3.1](#31-which-operations-actually-read-the-frame), [§6.4](#64-effort-and-risk) |
| — | *(not reported)* | **New hazard.** `deletion_scan` output is bit-identical across `frame=+1/+3/−1/−3`. For `frame=+3` the recorded idiom deletes `AAT GCC CGG GTT TAA ACC` where the true codons are `ATG CCC GGG TTT AAA CCC` | [§3.3](#33-the-frame-is-annotated-validated-and-ignored) |

---

## 1. What v1 actually did

Grounded in `poolparty/deletion_scan_orf_pool.py`, `insertion_scan_orf_pool.py`,
`orf_pool.py`, and the two test suites, which pass:

```
$ cd /mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty
$ python3 -m pytest tests/test_deletion_scan_orf_pool.py tests/test_insertion_scan_orf_pool.py -q
153 passed in 0.76s
```

### 1.1 How frame is established

Not by a `frame` parameter. `ORFPool.__init__` takes **nucleotide** indices
`orf_start` / `orf_end`, splits the input into three literal strings, and requires the
middle one to be a multiple of three:

```python
self.upstream_flank = seq[:orf_start]
orf_seq             = seq[orf_start:orf_end]
self.downstream_flank = seq[orf_end:]
if len(orf_seq) % 3 != 0:
    raise ValueError(...)          # orf_pool.py:110-117
```

Consequences, all of them structural:

- **Partial codons are impossible.** v1 raises rather than tolerating a remainder.
- **There is no frame offset.** Frame ≡ `orf_start`. Frames 2 and 3 are expressed by
  moving `orf_start`, not by an offset within a region.
- **There is no strand.** v1 has no negative-frame concept anywhere in these two classes.
- **Flanks are inert.** They are sliced off before any operation and re-attached by
  `_reassemble_with_flanks`, so an indel can never touch a UTR.

### 1.2 Parameters, units, semantics

Both classes offer two mutually exclusive interfaces. **Every positional parameter is in
codon units.**

| Parameter | Units | Default | Semantics |
|---|---|---|---|
| `deletion_size` / `insert_seq` | codons / DNA ×3 | required | `>0`, `≤ num_codons`; insert must be `%3 == 0` |
| `start` | codons | `0` | first window start |
| `end` | codons | `num_codons` | exclusive; last valid start is `end − W` |
| `step_size` | codons | `1` | stride |
| `positions` | codons | `None` | explicit list; no duplicates; window must fit |
| `position_weights` | — | uniform | **random mode only**; raises with `mode='sequential'` |
| `mark_changes` | — | `True` (del) / `False` (ins) | deletion: gap-mark vs remove. insertion: `swapcase()` the insert |
| `deletion_character` | — | `'-'` | one char, repeated ×3 per codon |
| `orf_start` / `orf_end` | nt | `0` / `len(seq)` | ORF extent within the construct |
| `insert_or_overwrite` | — | `'overwrite'` | `'insert'` splices between codons; `'overwrite'` replaces them |

**State counts.**

```python
# deletion, range interface
self.positions = list(range(start, min(end, L) - W + 1, step_size))
# insertion, overwrite
self.positions = list(range(start, min(end, L) - W + 1, step_size))
# insertion, insert  -- note the +1: L+1 insertion points for L codons
self.positions = list(range(start, min(end, L) + 1, step_size))
```

`num_internal_states = len(self.positions)` in every case — a plain count. This matters
for §6: it is already a uniform branching factor.

**Output length.** Deletion with `mark_changes=True` is length-preserving (each deleted
codon becomes `'---'`); with `False` the sequence shortens by `3W`. Insertion in
`overwrite` mode preserves length; in `insert` mode it grows by `len(insert_seq)`.

**Design cards.** `codon_pos` (codon units), `codon_pos_abs` (absolute nt, flank-aware),
`del_codons` / `insert`, and `del_aa` / `insert_aa` — the last two are always `None` in
practice, because both classes pass `codon_table=None` to `ORFPool`:

```python
codon_table=None,  # Scan pools don't need codon lookups
```

so the `if hasattr(self, 'codon_table') and self.codon_table:` branch that would populate
them is dead code. **v1's amino-acid card fields never fired.** They are not a capability
to migrate.

---

## 2. What v2 can do today

All of the following was executed against the working tree at
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty`
(poolparty 0.1.0). Progress-bar lines are stripped from the transcripts.

### 2.1 The mechanism

`PositionsType = Sequence[Integral] | slice | None`, resolved by the shared helper
(`utils/seq_utils.py:38-40`):

```python
if isinstance(positions, slice):
    start, stop, step = positions.indices(num_positions)
    return [min_position + i for i in range(start, stop, step)]
```

The step is honoured. `build_scan_cache` then sets
`num_states = len(indices)` (`utils/scan_utils.py:54-62`). Both confirmed. So a codon
stride is declarative and tool-resolved, exactly as B1 records.

### 2.2 The recorded idiom is unsafe without `region=`

**This is the correction that matters most.** `positions` is resolved against the
**whole pool**, not against the ORF. With a 3′ flank present, the last state deletes into
the UTR:

```python
seq = "GG" + "ATGAAACCCGGGTTTTAA" + "CCCC"      # ORF at 2..20
pp.deletion_scan(seq, deletion_length=3, deletion_marker=None,
                 positions=slice(2, None, 3), mode="sequential")
```

```
num_states: 7          <-- a 6-codon ORF; the 7th state is not a codon deletion
  0: GGAAACCCGGGTTTTAACCCC
  ...
  5: GGATGAAACCCGGGTTTCCCC
  6: GGATGAAACCCGGGTTTTAAC   <-- deleted "CCC" from the 3'UTR
```

It gets worse for multi-codon deletions, where the naive 3′ bound is also wrong:

```
deletion_length=6, positions=slice(2, None, 3)  -> num_states: 6   (correct: 5)
deletion_length=6, positions=slice(2, 20,   3)  -> num_states: 6   (still wrong)
deletion_length=6, positions=slice(2, 15,   3)  -> num_states: 5   (correct)
```

The correct stop is `orf_end − deletion_length + 1`, not `orf_end`. The intuitive choice
is silently wrong, and it is wrong by producing *extra* library members, which is the
failure mode least likely to be noticed.

**`region=` removes the arithmetic entirely.** Positions become region-relative, and the
region bounds the window for any `deletion_length`:

```python
base = pp.from_seq(seq)
p0   = pp.annotate_orf(base, "orf", extent=(2, 20), frame=1)
pp.deletion_scan(p0, deletion_length=3, deletion_marker=None,
                 region="orf", positions=slice(0, None, 3), mode="sequential")
```

```
num_states: 6
  GGAAACCCGGGTTTTAACCCC     GGATGAAACCCTTTTAACCCC
  GGATGCCCGGGTTTTAACCCC     GGATGAAACCCGGGTAACCCC
  GGATGAAAGGGTTTTAACCCC     GGATGAAACCCGGGTTTCCCC
```

and `deletion_length=6` under the same call gives 5, with no bound to compute.

**The canonical recipe should be `annotate_orf` + `region=` + `slice(offset, None, 3)`,
not a bare slice.** With `region=`, the slice start is the frame offset *within the
region*, which is `0` for frame +1 regardless of where the ORF sits in the construct.

### 2.3 v1 ↔ v2: ten cases, byte-identical

Ten v1 configurations were enumerated in one process and the v2 equivalents in another,
then diffed. Construct: `GGGGG` + `ATGCCCGGGTTTAAACCC` (6 codons) + `TTTTT`,
`orf_start=5`, `orf_end=23`.

| case | v1 call | v2 call (all with `region="orf"`, `mode="sequential"`) | n | identical |
|---|---|---|---|---|
| D1 | `deletion_size=1, mark_changes=True` | `deletion_length=3, deletion_marker='-', positions=slice(0,None,3)` | 6 | **YES** |
| D2 | `deletion_size=1, mark_changes=False` | `deletion_length=3, deletion_marker=None, positions=slice(0,None,3)` | 6 | **YES** |
| D3 | `deletion_size=2, mark_changes=False` | `deletion_length=6, deletion_marker=None, positions=slice(0,None,3)` | 5 | **YES** |
| D4 | `deletion_size=1, step_size=2` | `deletion_length=3, positions=slice(0,None,6)` | 3 | **YES** |
| D5 | `deletion_size=1, start=2, end=5` | `deletion_length=3, positions=slice(6,15,3)` | 3 | **YES** |
| D6 | `deletion_size=1, positions=[1,3]` | `deletion_length=3, deletion_marker='-', positions=[3,9]` | 2 | **YES** |
| I1 | `insert_seq="TAG", 'insert'` | `insertion_scan("TAG", positions=slice(0,None,3), replace=False)` | 7 | **YES** |
| I2 | `insert_seq="TAG", 'overwrite'` | `insertion_scan("TAG", positions=slice(0,None,3), replace=True)` | 6 | **YES** |
| I3 | `insert_seq="TAGTAG", 'overwrite'` | `insertion_scan("TAGTAG", positions=slice(0,None,3), replace=True)` | 5 | **YES** |
| I4 | `insert_seq="TAG", mark_changes=True` | `insertion_scan("tag", positions=slice(0,None,3), replace=False)` | 7 | **YES** |

```
case                    v1_n  v2_n  identical?
----------------------------------------------------------
D1_marked_1codon           6     6  YES
D2_unmarked_1codon         6     6  YES
D3_unmarked_2codon         5     5  YES
D4_step2_1codon            3     3  YES
D5_start2_end5             3     3  YES
D6_positions_1_3           2     2  YES
I1_insert_TAG              7     7  YES
I2_overwrite_TAG           6     6  YES
I3_overwrite_2codon        5     5  YES
I4_insert_marked           7     7  YES
----------------------------------------------------------
ALL 10 CASES IDENTICAL: True
```

Scripts: `v1_gen.py`, `v2_gen.py` (scratchpad; reproduced in [Appendix A](#appendix-a-reproduction)).

Note the unit conversions the user has to perform by hand: codon → nt for
`deletion_size` (×3), `step_size` (×3), `positions` (×3), and `start`/`end`
(×3, with the `−W+1` correction). Every one is a place to get it wrong.

### 2.4 Random mode and protein-level confirmation

Random mode honours the slice — 200 draws yield exactly the six in-frame deletions:

```
distinct sequences sampled: 6
    50  GGGGGATGCCCGGGAAACCCTTTTT       32  GGGGGCCCGGGTTTAAACCCTTTTT
    35  GGGGGATGCCCGGGTTTCCCTTTTT       31  GGGGGATGCCCGGGTTTAAATTTTT
    26  GGGGGATGCCCTTTAAACCCTTTTT       26  GGGGGATGGGGTTTAAACCCTTTTT
```

Composing with `translate` gives the decisive check. WT protein is `MPGFKP`:

```
codon idiom (positions=slice(0,None,3))     frameless scan (positions=None)
  PGFKP                                       PGFKP  TGFKP  IGFKP  MGFKP
  MGFKP                                       MRFKP  MPFKP  MPFKP  MPVKP
  MPFKP                                       MPGKP  MPGKP  MPG*P  MPGLP
  MPGKP                                       MPGFP  MPGFT  MPGFN  MPGFK
  MPGFP
  MPGFK
```

Six single-residue deletions on the left; sixteen sequences on the right including
substitutions and a premature stop. The idiom is genuinely frame-preserving.

---

## 3. The delta

### 3.1 Which operations actually read the frame

```
$ grep -rn "\.frame\b" src/poolparty/ --include=*.py
src/poolparty/orf_ops/translate.py:42:        return registered_region.frame
src/poolparty/orf_ops/stylize_orf.py:46:     return registered_region.frame
src/poolparty/orf_ops/mutagenize_orf.py:50:  return registered_region.frame
```

Three operations, each via a private `_resolve_frame(region, frame)` helper defined
separately in all three files. **The logic is identical in all three; the text is not.**
Parsed and AST-normalized with docstrings stripped:

| file | AST hash | difference from `mutagenize_orf` |
|---|---|---|
| `mutagenize_orf.py:21-56` | `24fa85ec76` | — (reference) |
| `translate.py:16-48` | `24fa85ec76` | **AST-identical**; three explanatory comments stripped |
| `stylize_orf.py:15-52` | `30c44df929` | one added statement: a function-local `from ..party import get_active_party` |

```
$ diff mutagenize_orf/_resolve_frame  stylize_orf/_resolve_frame
7a8,9
>     from ..party import get_active_party
>
```

That single line is the whole difference, and it is **redundant, not load-bearing**
(see §6.4). Two identical copies plus one with two vestigial import lines is still the strongest
argument that a fourth and fifth consumer belong in the same family rather than bolted
onto `scan_ops`.

`OrfRegion`'s docstring names only two ("`stylize_orf()` and `mutagenize_orf()`"). The
briefing's count of three was right; its attribution to the docstring was not.

And `scan_ops/` has **zero** ORF vocabulary:

```
$ grep -rniE "codon|frame|orf|triplet|inframe" src/poolparty/scan_ops/*.py
$ echo $?
1
```

### 3.2 `region=` does constrain, but does not frame

```python
p0 = pp.annotate_orf(base, "orf", extent=(2,20), frame=1)
pp.deletion_scan(p0, deletion_length=3, deletion_marker=None, region="orf", mode="sequential")
```

```
num_states: 16          <-- 18 − 3 + 1; a nucleotide scan, correctly bounded to the ORF
```

So `region=` solves the **boundary** problem and leaves the **frame** problem untouched;
the slice solves the frame problem and leaves the boundary problem untouched. Both are
needed, and nothing tells the user that.

On tags: `generate_library` **does** emit them —

```
0  GG<orf>AAACCCGGGTTTTAA</orf>CCCC
1  GG<orf>ATGCCCGGGTTTTAA</orf>CCCC
```

— while `to_df` / `to_file` strip them via `write_tags: bool = False`
(`pool_mixins/export_mixin.py:93`). `generate_library` has no `write_tags` parameter at
all; its signature is `(pool, num_cycles, num_seqs, seed, init_state, seqs_only,
_include_inline_styles, discard_null_seqs, max_iterations, min_acceptance_rate,
attempts_per_rate_assessment)`.

### 3.3 The frame is annotated, validated, and ignored

Region `AATGCCCGGGTTTAAACCC` (19 nt) inside a larger construct, annotated at three
different frames. `mutagenize_orf` is frame-aware, so it is the ground truth for "which
triplets are codons":

```
--- FRAME-AWARE ground truth (mutagenize_orf, num_mutations=1) ---
  frame=+3 offset=1: codons = [ATG, CCC, GGG, TTT, AAA]          (+ CCC, deduped)
  frame=-3 offset=1: codons = [GGT, TTA, AAC, CCG, GGC, ATT]     (revcomp, 3'->5')

--- deletion_scan with the recorded slice(0,None,3) ---
  frame=+3: n=6 deleted = ['AAT', 'GCC', 'CGG', 'GTT', 'TAA', 'ACC']
  frame=-3: n=6 deleted = ['AAT', 'GCC', 'CGG', 'GTT', 'TAA', 'ACC']
```

Two things here.

**One.** For `frame=+3`, **none of the six deleted triplets is a codon.** Every state is
an out-of-frame 3-nt deletion. The correct offset is `(4 − |frame|) % 3 = 1`:

```
  naive   slice(0,None,3): ['AAT', 'GCC', 'CGG', 'GTT', 'TAA', 'ACC']
  correct slice(1,None,3): ['ATG', 'CCC', 'GGG', 'TTT', 'AAA', 'CCC']   <-- matches truth
```

**Two.** The `+3` and `−3` rows are byte-identical. `deletion_scan` cannot see the frame,
including its **sign**. (For negative frames the offset is taken off the 3′ end, so
`slice(0, None, 3)` happens to be right there — but the enumeration *order* is reversed
relative to `mutagenize_orf`, which numbers codon 0 as the rightmost. Two ORF operations
in the same Party will disagree about which one is "codon 5".)

The user has annotated the frame. v2 has validated it, stored it in an immutable
`OrfRegion`, and refuses to let it be changed. Then it throws it away. No documentation
fixes that.

### 3.4 The delta, classified

**(a) v1 capabilities v2 cannot express at all**

| Capability | Notes |
|---|---|
| `position_weights` — non-uniform sampling over explicit positions | The **only** true gap. `grep -rn "weight" src/poolparty/` returns nothing. Not emulable: `validate_positions` rejects duplicates (`ValueError: Positions must not contain duplicates`), so a position cannot be repeated to bias it. **Not ORF-specific** — it is a general scan-sampling feature, and it should not be scoped into an ORF change |

`del_aa` / `insert_aa` are *not* in this column: §1.2 shows they were dead code in v1.

**(b) Expressible, but awkwardly or undocumented**

| Capability | Cost today |
|---|---|
| Codon units | Every quantity ×3 by hand; `start`/`end` also need `−W+1` |
| Frame offset | Hand-computed `(4−|frame|)%3`, with a different rule per strand sign; **silently wrong if omitted** (§3.3) |
| Boundary safety | Requires knowing `region=` is mandatory, not optional (§2.2) |
| The slice idiom itself | `docs/operations/deletion_scan.rst:48` types `positions` as `list[int] \| None` and line 50 calls it an *"Explicit list of window start positions"*. `grep -rn "slice" docs/operations/deletion_scan.rst docs/operations/insertion_scan.rst` → **no matches**. Meanwhile `codon_positions=slice(1, 56)` is documented in `docs/tutorials/dms_gb1.rst:44`. Same capability, documented for one operation and hidden for two |
| Two deletion lengths in one Party | Impossible — B3, §5 |
| Codon-unit design cards | `position_index` coincides with the codon index only because the slice makes it so; `start`/`end` are region-relative nt |

**(c) Already equivalent**

Everything else: both interfaces, both marking modes, insert vs overwrite, multi-codon,
flank preservation, both selection modes, state counts, and — per §2.3 — the actual output
bytes.

---

## 4. Recommendation

**Add ergonomics. Do not port. Fix B3 in the same change. Then document.**

Justification, in order of weight:

1. **A port cannot be justified on capability.** Ten of ten cases are byte-identical
   today (§2.3). Any claim that v2 "lost" codon-aware indels is false and would not
   survive review.
2. **Documentation alone is not sufficient**, which is where I part company with the
   simplest reading of B1. §3.3 is a silent-wrong-answer path with the frame sitting
   right there in the Party, unused. §2.2 is a second silent-wrong-answer path in the
   exact recipe currently recorded in `FINDINGS.md`. Prose that says "remember to pass
   `region=` and to compute `(4−|frame|)%3`, and note the rule inverts on the minus
   strand" is a worse artefact than ten lines of code that do it.
3. **The change is small and purely compositional** — no new `Operation` subclass, no new
   state mechanism, no change to any existing signature's behaviour.

### 4.1 For the manuscript

The honest framing, which the evidence supports:

> Codon-aware indel scanning is expressible in PoolParty by composing `annotate_orf` with
> a strided `positions` slice; the v1 classes' behaviour is reproduced exactly. What v1
> provided and v2 does not is the *vocabulary* — codon units and an explicit reading
> frame on the scan operations themselves.

**Do not** describe this as a structural consequence of the DAG architecture. §6 shows
the state count is a plain position count, decomposes as a product, and fits the lazy
mechanism without strain. B1 is right about that, and it remains the reason the framing
matters: the v1 code is public history.

---

## 5. B3: the `_del` collision

Reproduced exactly:

```python
a = pp.deletion_scan(p0, deletion_length=3, ...)
b = pp.deletion_scan(a,  deletion_length=6, ...)
```
```
ValueError: Region '_del' already registered with seq_length=3, cannot re-register
with seq_length=6. Region lengths must be consistent within a Party.
```

Same length twice composes fine, and the product decomposition holds:

```
num_states: 36        # 6 × 6
```

The precedent is stronger than briefed — **both** multiscan operations already expose the
parameter:

```
insertion_multiscan: [..., 'positions', 'region', 'names', 'replace', ...]
deletion_multiscan : [..., 'positions', 'region', 'names', 'min_spacing', ...]
```

**Fix it in the same change.** Two reasons beyond tidiness. v1's `deletion_size` was a
first-class parameter, so a codon wrapper makes multi-length deletion scans the natural
thing to reach for — shipping the wrapper without the fix ships a wall. And it is already
constraining published work: `code_comparison/valiant_brca1_pep.md` records that
`brca1_pep` was chosen partly because *"`deletion_scan` hard-codes its internal region
name, so two different `deletion_length` values collide in one Party … The other examples
combine `1del` with `2del0`/`2del1`."* Four of VaLiAnT's five examples are currently out
of reach for this reason.

The fix is ~3 lines, mirroring `deletion_multiscan`:

```python
def deletion_scan(..., names: Optional[str] = None, ...):
    marker_name = names if names is not None else "_del"
```

`insertion_scan`'s `_ins` / `_rep` have the same latent problem and should get the same
parameter.

---

## 6. Design

### 6.1 Does it fit the lazy state mechanism?

Yes, and without strain.

- `build_scan_cache` returns `num_states = len(resolved_positions)` — a plain count
  (`utils/scan_utils.py:54-62`). A codon scan just resolves a different list.
- The branching factor is uniform by construction, so the `synonymous` /
  `UNIFORM_MUTATION_TYPES` problem does not arise. `UNIFORM_MUTATION_TYPES` maps
  `any_codon: 63, nonsynonymous_first: 20, missense_only_first: 19, nonsense: 3`;
  `synonymous` is absent because its count varies by residue. **A positional scan has no
  per-position variation at all** — every codon offers exactly one deletion.
- Verified empirically: two composed deletion scans give `num_states: 36 = 6 × 6`.

### 6.2 Evaluating `frame=` on `deletion_scan` / `insertion_scan`

This was the briefing's prior guess. **I recommend against it**, on two grounds.

**It creates a hidden units switch.** `positions` is nucleotide-indexed. If `frame=` leaves
it that way, the user still writes the `×3` stride and the parameter buys almost nothing.
If `frame=` re-interprets `positions` as codon-indexed, then one argument silently changes
another's units — invisible to `beartype`, invisible at the call site, and a migration
hazard for existing code. `mutagenize_orf` deliberately avoids this by giving its
codon-indexed parameter a **different name**: `codon_positions`, never `positions`. And
`deletion_length` would remain nucleotides, so a single call would carry both units.

**It puts codon semantics in the nucleotide layer.** `scan_ops/` currently has zero ORF
vocabulary (§3.1). That split is not an accident — it is the same split as
`orf_ops/` vs the rest, and it is what makes "ORF-awareness is a property of specific
operations" a true and teachable statement. Adding `frame=` to `scan_ops` erases it.

Also considered and rejected: **making `region=` implicitly frame-aware** when the region
is an `OrfRegion`. It silently changes existing behaviour and existing state counts, and
it forecloses the legitimate case of a deliberate nucleotide-level scan inside an ORF.

### 6.3 Recommended: `deletion_scan_orf` / `insertion_scan_orf` in `orf_ops/`

Two new functions, named the way v2 already names ORF operations
(`mutagenize_orf`, `stylize_orf`, `annotate_orf`) — and, conveniently, the way v1 named
these classes.

```python
@beartype
def deletion_scan_orf(
    pool: Union[Pool, str],
    num_codons: Integral = 1,                 # codon units (v1: deletion_size)
    region: RegionType = None,
    codon_positions: Union[Sequence[Integral], slice, None] = None,
    frame: Optional[int] = None,              # None -> read from OrfRegion
    deletion_marker: Optional[str] = "-",
    names: Optional[str] = None,              # B3
    prefix: Optional[str] = None,
    mode: ModeType = "random",
    num_states: Optional[Integral] = None,
    style: Optional[str] = None,
    iter_order: Optional[Real] = None,
    cards: CardsType = None,
    _factory_name: Optional[str] = "deletion_scan_orf",
) -> Pool:
    resolved_frame = _resolve_frame(region, frame)          # shared helper, verbatim
    offset   = (4 - abs(resolved_frame)) % 3
    reverse  = resolved_frame < 0
    span     = _region_length(pool, region)
    n_codons = (span - offset) // 3

    codon_idx = validate_positions(codon_positions,         # codon units, slice honoured
                                   max_position=n_codons - num_codons,
                                   min_position=0)
    nt_start  = 0 if reverse else offset                    # negative frames trim the 3' end
    nt_positions = [nt_start + 3 * c for c in codon_idx]

    return deletion_scan(pool, deletion_length=3 * int(num_codons),
                         deletion_marker=deletion_marker, positions=nt_positions,
                         region=region, names=names, prefix=prefix, mode=mode,
                         num_states=num_states, style=style, iter_order=iter_order,
                         cards=cards, _factory_name=_factory_name)
```

`insertion_scan_orf` is the same shape over `insertion_scan`, with `replace: bool = False`
and `max_position = n_codons - ins_codons` when replacing, `n_codons` when inserting
(the `+1` that gives v1's `L+1` insertion points falls out of `validate_positions`'
inclusive `max_position`).

**Why this shape.**

| v2 principle | How it is respected |
|---|---|
| Lazy state decomposition | `num_states = len(codon_idx)`. Nothing new; `deletion_scan` still owns the state |
| Operations are functions returning Pools, composed into a DAG | Pure delegation. **No new `Operation` subclass**, no `Pool` subclass — the v1 class-per-operation shape is not reintroduced |
| Prefer extending composition over parallel implementation | It composes over `deletion_scan`, which already composes `region_scan` + `replace_region`. Three layers, one mechanism |
| `(Seq, dict)` with a design card | Cards pass through unchanged; add `codon_position` (codon units) and `frame`, alongside the existing nt `start`/`end` |
| Regions are Party-global with fixed `seq_length` | Untouched, and `names=` makes the collision addressable |
| ORF-awareness belongs to specific operations | It joins that set explicitly, in `orf_ops/`, sharing `_resolve_frame` with its three siblings |

**Reverse-strand semantics are the one real decision.** `mutagenize_orf` numbers codon 0
as the **rightmost** complete codon for negative frames (`codon_region_end = mol_end −
frame_offset`, `mutagenize_orf.py:541-548`). The new operations **must** adopt that
convention exactly, or two ORF operations in one Party will disagree about which codon is
which. This is the highest-risk detail in the whole design and needs a dedicated test
(§7).

**Do not scope `position_weights` into this change.** It is the only true gap (§3.4a), but
it is a general scan-sampling feature, not an ORF one. It belongs on `positions` across
all scan ops, or nowhere.

### 6.4 Effort and risk

| Item | Size |
|---|---|
| `orf_ops/deletion_scan_orf.py` | ~110 lines, mostly docstring (cf. `annotate_orf.py` at 172) |
| `orf_ops/insertion_scan_orf.py` | ~130 lines |
| Extract `_resolve_frame` to `orf_ops/_frame.py`, import in all 5 | ~40 lines net; reconciles two identical copies plus one with two vestigial local imports (§6.4). Verified safe |
| `orf_ops/__init__.py`, `poolparty/__init__.py` exports | ~8 lines |
| `dna_mixin` methods | ~60 lines |
| `names=` on `deletion_scan` / `insertion_scan` (B3) | ~6 lines |
| Docs: 2 new `.rst`; fix `positions` type + description on the 2 existing scan pages | ~150 lines |
| Tests | ~40 cases |

**Roughly one focused day**, with the reverse-strand convention as the only part that
warrants slow, careful work.

**Is the `_resolve_frame` extraction safe?** Yes — tested, not assumed. `stylize_orf`'s
function-local `from ..party import get_active_party` looks like a circular-import
workaround, but it is not one:

- `party.py` imports **nothing** from `orf_ops`, at module scope or deferred
  (`grep -n "orf_ops\|mutagenize_orf\|stylize_orf\|translate" party.py` → no matches).
  There is no cycle through `orf_ops` to dodge.
- `mutagenize_orf.py` and `translate.py` — its two siblings in the same package, importing
  the same symbol for the same purpose — **already do it at module scope** and always have.

`stylize_orf.py` defers the import in **two** places, not one:

| line | site |
|---|---|
| `stylize_orf.py:22` | inside `_resolve_frame` |
| `stylize_orf.py:142` | inside `StylizeOrfOp.__init__`, guarding `get_active_party()  # Ensure we're in a Party context` |

Confirmed empirically, and for both sites. The package was copied to a scratch tree, both
deferred imports removed and replaced with a single module-level import in the copy, and
the suite run against the copy:

```
$ grep -n "get_active_party" <scratch>/orf_ops/stylize_orf.py   # fully patched copy
9:from ..party import get_active_party
35:    party = get_active_party()
141:        get_active_party()  # Ensure we're in a Party context

$ PYTHONPATH=<scratch> python3 -c "import poolparty"
IMPORT OK

$ PYTHONPATH=<scratch> python3 -m pytest tests/ -q          # both sites hoisted
2915 passed, 14 xfailed in 26.17s

$ PYTHONPATH=<scratch1> python3 -m pytest tests/ -q         # only _resolve_frame hoisted
2915 passed, 14 xfailed in 28.86s

$ python3 -m pytest tests/ -q                               # unpatched baseline
2915 passed, 14 xfailed in 27.47s
```

The same result in all three. Both deferred imports are vestigial. The extraction is a
reconciliation of two identical copies and one carrying two redundant lines — not a
load-order minefield. *(Only scratch copies were patched; the repo was not.)*

**Risks.**

1. **Reverse-strand codon numbering disagreeing with `mutagenize_orf`.** Highest. Test
   directly against it, not against a hand-derived expectation.
2. **Partial codons at region ends.** v1 raised; v2 regions permit remainders and
   `mutagenize_orf` tolerates them via `frame_offset`. Follow `mutagenize_orf`, not v1 —
   consistency inside v2 beats fidelity to v1.
3. **`deletion_length` validation.** `deletion_scan` currently validates against
   `pool.seq_length`, not the region length (`deletion_scan.py:76-80`). A codon wrapper
   makes region-relative validation the thing users will expect.
4. **Marker-name collisions** beyond `_del` — `_ins` / `_rep` need `names=` too.
5. **Scope creep** toward `position_weights`. Resist (§6.3).

---

## 7. Test strategy

**Transfers directly from v1** — the 153 passing tests are a usable spec. In particular:

- All of `TestMarkedDeletionMode`, `TestUnmarkedDeletionMode`, `TestOverwriteMode`,
  `TestInsertMode`, `TestMultiCodonOperations`, `TestReadingFrameIntegrity`,
  `TestSequenceLength`, `TestFlankingRegions` — after rewriting codon parameters onto the
  new signature.
- `TestStateSpaceCalculation` and the `num_states_*` tests transfer as-is in intent.
- **The ten cases in §2.3 are a ready-made golden set**, and they are already known to be
  byte-identical, so they can be committed as regression fixtures on day one.

**Does not transfer:**

- `TestPositionBasedInterface::test_weighted_positions*` and the `position_weights`
  validation tests (~5) — no v2 equivalent, and out of scope.
- `TestRepr` — v2 has no Pool subclass to repr.
- `test_metadata_levels` — v2 uses `cards`, a different contract.
- `generate_seqs` seed-reproducibility tests — different RNG.
- `TestPoolChaining::test_pool_as_seq_input` — v2 pools chain by construction.

**New tests required, none of which exist in v1** (v1 had no frame concept, so none of
this could have been tested there):

1. **Frame agreement with `mutagenize_orf`** for all six frames `±1, ±2, ±3` on the same
   annotated region: the codons deleted by `deletion_scan_orf` must equal the codons
   `mutagenize_orf` reports in `wt_codons`, **in the same order**. This is the §3.3
   regression and the §6.4 risk-1 test in one.
2. **Frame is read, not assumed** — `annotate_orf(frame=3)` followed by
   `deletion_scan_orf(...)` with no explicit `frame` must produce the frame-3 result, and
   must differ from the frame-1 result.
3. **Region boundary non-overrun** for `num_codons` ∈ {1,2,3} with both 5′ and 3′ flanks
   present — the §2.2 regression. Assert flanks are byte-identical in every state.
4. **State-count product** — two composed ORF scans give `n₁ × n₂`.
5. **`names=` collision** — two `deletion_scan_orf` calls with different `num_codons` and
   distinct `names` succeed; with colliding names, raise.
6. **Protein round-trip** — `translate` of each state has exactly `num_codons` residues
   fewer than WT, with all others unchanged (the §2.4 check, as an assertion).
7. **Partial-codon regions** — a region whose length is not `0 mod 3` behaves per
   `mutagenize_orf`, for both strand signs.
8. **Plain `Region` (not `OrfRegion`) with no explicit `frame`** raises the same error
   `mutagenize_orf` raises.
9. **VaLiAnT `brca1_pep`** — the existing reproduction must still yield 2,339/2,339 when
   rewritten onto `deletion_scan_orf`. It is a minus-strand library, so it exercises
   risk 1 against published ground truth rather than against ourselves.

---

## Appendix A: reproduction

Environment: `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/.venv`, or plain
`python3` — poolparty 0.1.0 resolves from
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty/src`.
v1 is imported by prepending `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty` to
`sys.path`; the two packages share the name `poolparty` and must be run in separate
processes.

```python
# v1 side
import sys; sys.path.insert(0, '/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty')
from poolparty import DeletionScanORFPool
FULL = "GGGGG" + "ATGCCCGGGTTTAAACCC" + "TTTTT"
p = DeletionScanORFPool(FULL, deletion_size=1, orf_start=5, orf_end=23,
                        mark_changes=False, mode='sequential')
seqs_v1 = [(p.set_state(s), p.seq)[1] for s in range(p.num_internal_states)]

# v2 side (separate process)
import poolparty as pp
with pp.Party():
    base = pp.from_seq(FULL)
    p0 = pp.annotate_orf(base, "orf", extent=(5, 23), frame=1)
    q  = pp.deletion_scan(p0, deletion_length=3, deletion_marker=None, region="orf",
                          positions=slice(0, None, 3), mode="sequential")
    seqs_v2 = list(q.to_df(num_cycles=1)['seq'])

assert seqs_v1 == seqs_v2
```

Full generators for all ten cases: `v1_gen.py` / `v2_gen.py`, session scratchpad
`/tmp/claude-1000/-mnt-c-Users-zhliu-Desktop-KinneyLab/626bd452-cf38-4426-8103-bd76df8c926b/scratchpad/`.

---

*Survey by session `kinneylab-4a`, 2026-08-18. Neither PoolParty repo was modified.*
