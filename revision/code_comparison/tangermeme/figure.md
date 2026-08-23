# Supplementary figure — SpliceAI GT arm, PoolParty and tangermeme

Design record. See `comparison.md` for the verification and the corrections made
along the way.

**Orientation:** undecided. The widest line is ~88 characters, which fits full
width at 7 pt (89 chars) in portrait. Since (B) and (C) became code-only they are
~6 lines each, so landscape side-by-side is now viable as well — at the cost of
reflowing the code to roughly 66 characters per column.

**This is not a conciseness figure.** Three statements each side. If the layout
invites a line count, it is the wrong layout.

---

## Layout

```
┌──────────────────────────────────────────────────────────────────────────┐
│ (A)  The library                                                         │
│                                                                          │
│      201 bp region, canonical 5'ss at position 100                       │
│      ┌────────────┬───┬──────────────────┐                               │
│      │  51 - 89   │ ● │    107 - 167     │   100 scan positions          │
│      └────────────┴───┴──────────────────┘                               │
│                     ▲                                                    │
│              canonical 5'ss (skipped)                                    │
│                                                                          │
│      2,000 cryptic 5'ss  x  100 positions  =  200,000 sequences          │
└──────────────────────────────────────────────────────────────────────────┘
┌──────────────────────────────────────────────────────────────────────────┐
│ (B)  PoolParty                                                           │
│                                                                          │
│  cryptic_pool = pp.from_seqs(cryptic_sites, mode="sequential",           │
│                              cards={"seq": "cryptic_sequence"})          │
│  library = pp.from_seq(target_region).replacement_scan(                  │
│      cryptic_pool, positions=SCAN_POSITIONS, mode="sequential",          │
│      cards={"start": "cryptic_position"})                                │
│  df = library.generate_library()                                         │
└──────────────────────────────────────────────────────────────────────────┘
┌──────────────────────────────────────────────────────────────────────────┐
│ (C)  tangermeme                                                          │
│                                                                          │
│  target_ohe  = one_hot_encode(target_region).unsqueeze(0) \              │
│                    .expand(len(cryptic_sites), -1, -1)                   │
│  cryptic_ohe = torch.stack([one_hot_encode(s) for s in cryptic_sites])   │
│  library_ohe = torch.stack([substitute(target_ohe, cryptic_ohe, start=p) │
│                             for p in SCAN_POSITIONS])                    │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## Panel notes

**(A)** The gap in the position bar is the point: positions 90–106 are skipped
because the canonical site sits at 100. `SCAN_POSITIONS` is an explicit
non-contiguous list, not a range.

**(B) and (C)** Code only. Neither panel displays output. An earlier version showed
PoolParty's output table against tangermeme's returned tensor shape; that made the
PoolParty panel read as the heavier one, because a table with headers and rows
outweighs a 4-tuple even when both are faithful. Panel (A) already states what the
library contains, so the panels need only show how each tool builds it.

The return types stay legible from the code — `generate_library()` against
`torch.stack` — and the `cards=` arguments name the columns PoolParty attaches.
Those arguments are load-bearing: without them the figure claims only that both
tools produce 200,000 sequences, which is the conciseness reading this design
rejects. `mode="sequential"` is load-bearing for the same reason — without it the
pools sample rather than enumerate.

`generate_library()` rather than `to_df()`: it is the documented idiom (78 uses in
`docs/`, 8 in `examples/`, against 6 and 0 for `to_df`), and it is shorter, since
`num_cycles=1` is its default. Both return the same columns.

The dict form of `cards=` is deliberate. The list shorthand `cards=["start"]` also
works, but names the column `op[6]:replacement_scan(region_scan).start`. The dict
costs a few characters and buys the reader the meaning of what is recorded.

Do **not** put the `transpose`/`flatten` reorder, the `characters` decode, or the
DataFrame construction in panel C. All three are comparison scaffolding, and
including them would make the panel misrepresent how the package is used. They
live in `repro.py`, marked as such.

Both panels use the shared input names `target_region`, `cryptic_sites`,
`SCAN_POSITIONS` so the correspondence is followable. Panel C departs from
tangermeme's own `X` convention deliberately — a reader following two panels should
not need one tool's variable idiom. The `_ohe` suffix carries useful information for
free: it makes visible that tangermeme works in tensors where PoolParty works in
strings.

---

## Caption

> **Supplementary Figure N. One perturbation library, in PoolParty and in
> tangermeme.**
>
> **(A)** The library. Each of 2,000 cryptic 5′ splice sites of graded strength is
> substituted at each of 100 positions across a 201 bp region carrying the canonical
> 5′ splice site of SMN2 exon 7, giving 200,000 sequences. Positions immediately
> around the canonical site are excluded. This is the GT arm of the surrogate-modelling
> library of Fig. N; the matched control arm is built the same way.
>
> **(B)** In PoolParty, the cryptic site and its position are returned as columns
> beside each sequence, named by the `cards` arguments. No sequence is generated
> until requested.
>
> **(C)** In tangermeme, the substitutions are applied one batch per position and
> returned as a one-hot tensor. The site and position are not returned; each is
> recoverable from a sequence's index along one axis.
>
> The two libraries contain the identical 200,000 sequences, verified against the
> published design cards. The tools differ in what they return rather than in what
> they can express: a table with the parameters attached, suited to analysis that
> uses them as covariates, or a tensor in the shape a predictive model consumes.

---

## Rejected alternatives

**Tutorial B3 (Spacing) as the target.** `space()` requires a trained model and
returns predictions; the perturbed sequences are built in a loop and discarded. It
also needs a Beluga checkpoint that is not in the repository.

**Tutorial B7 (Cartesian Product) as the target.** `apply_product` crosses sequences
against model side-inputs — cell state, read depth — not motifs against positions.

**A conciseness framing.** Three statements each side.

**A second panel conceding model-guided design.** Considered so the figure would not
be one-sided. Unnecessary once the comparison is even-handed; the concession stands
in Table 2 (`model_guided_variants`: tangermeme ●, PoolParty ◐).

**The full Example 3.** Would need the PWM ladder, the genome slice, both arms and
the derived card columns — most of which is not a tool comparison.
