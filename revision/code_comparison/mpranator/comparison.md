# MPRAnator "MPRA Motif Use Case" — exact reproduction in PoolParty

**Result: 976 of 976 arrangement sequences reproduced byte-for-byte, on the first
attempt, with no parameter tuning.** MPRAnator's own code was run to establish the
ground truth and returned the **5,856** total its documentation states.

Answers referee **R2 2b** (overlap of pool elements, and investigation of the unique
elements) and contributes to **R1 #5**.

---

## 1. The example, and why this one

MPRAnator's **"MPRA Motif Use Case Example"** — the fully worked design on its live
documentation page, chosen because it is the only one of MPRAnator's five documented
examples that **states an expected count**.

From the page verbatim: *"Inspired by the results presented by Nguyen et al, we would
like to investigate the effects of AP1 (TGACTCA), ELK1 (ACCGGAAGT) and RFX
(CGTTGCTAGGCAACG) on gene expression."*

This is MPRA motif-arrangement design — the paper's second declared application, and
a capability set almost disjoint from VaLiAnT's (which has no motif or element
concept at all).

---

## 2. Exact source

**Documentation page** — the design, the parameters and both background sequences:

```
https://www.genomegeek.com/MPRA/documentation/        (section "MPRA Motif Use Case Example")
```

**Source code** — the Django application that *is* the web service:

```
https://github.com/hemberg-lab/MPRAnator
  oligo.py          position enumeration and the four filters
  part1.py          generateCombinations, barcode_generator, getBarCodes,
                    createMPRAResultOutput
  myfunctions.py    complement / revcompl / findMatch helpers
```

| Fact | Value |
|---|---|
| Commits | **1** — `9969790d`, 2015-12-27 |
| Last push | 2015-12-27 (**predates the 2017 paper by ~18 months**) |
| Licence | **none** (all rights reserved) |
| PyPI / bioconda | absent |

Because the repository predates the paper, matching their documented total is what
validates the source as the code that produced it — see section 5.

**Licence note.** With no LICENCE file the code is all-rights-reserved. It was read
and executed locally to establish ground truth; none of it is vendored into this
repository. The algorithm is described here and re-expressed in PoolParty.

### Design parameters

| Parameter | Value |
|---|---|
| Motifs | AP1 `TGACTCA` (7 nt) · ELK1 `ACCGGAAGT` (9 nt) · RFX `CGTTGCTAGGCAACG` (15 nt) |
| Backgrounds | two non-regulatory mm9 tiles, **90 nt** each |
| Edge margin | >= 15 nt from each end |
| Spacing between motifs | min **6**, max **24** |
| Interval of substitution | **12** |
| Barcodes | 15 nt, GC 40-60%, edit distance 3, **6 per sequence** |
| Restriction sites | `CACGTG` and `CAATTG`, flanking the background |
| Tile | **117 nt** = 90 background + 6 + 6 + 15 barcode |
| **Stated total** | **5,856 sequences** |

Backgrounds, verbatim from the documentation page:

```
>bg1_chr6:77195320
tgtgtcttaaaaaaacaaacaaacaaacaaaatcccgaaataaaacacaacaaaaaaaaccccaccccataatcttcaggacagtctgtc
>bg2_chr9:37271330
gtatctactctctgcccttacaacctcctcccagaaagaataaaatgtttctcatcctggaagctacagtgtgtcacacagtatactctt
```

---

## 3. Mutation types, and how the total is computed

MPRAnator **substitutes** motifs into the background — it overwrites, so sequences
stay 90 nt. There is no insertion, deletion or point-mutation mutator here.

`part1.generateCombinations` takes all **subsets** of the three motifs: 3 singles +
3 pairs + 1 triplet = **7**. For each subset, `oligo.oligo()` enumerates position
tuples with `itertools.product(range(15, 90), repeat=k)` and applies four filters:

| # | Filter | Effect |
|---|---|---|
| 1 | placed motif must not extend past the end | start <= 90 - n |
| 2 | `testIfWithin` | gap between consecutive motifs in **[6, 24]** |
| 3 | `testIfFarFromRightEdge` | rightmost motif **end <= 74** (= 90-1-15), i.e. start <= 75 - n |
| 4 | `sorted(pos)[0] % 12 == 0` | the **leftmost** motif only sits on a multiple of 12 |

Filter 4 is the subtle one: it gates **only the leftmost** motif, acting as an anchor
at 24, 36, 48, 60. Reimplementing it as "every motif on a multiple of 12" changes
every count.

### Singles: 4 each

| Motif | n | start <= 75-n | multiples of 12 in [15, bound] | Count |
|---|---|---|---|---|
| AP1 | 7 | 68 | 24, 36, 48, 60 | **4** |
| ELK1 | 9 | 66 | 24, 36, 48, 60 | **4** |
| RFX | 15 | 60 | 24, 36, 48, 60 | **4** |

All three tie because the next anchor, 72, exceeds every bound. Motif length is
invisible here.

### Pairs: closed form

Anchor the left motif at `a`; the second motif's start lies in
`[a+n1+6, a+n1+24]` (19 gap choices) and must satisfy `start <= 75-n2`. So

```
count at anchor a  =  max(0, min(19, 70 - n1 - n2 - a))
```

which depends only on `n1+n2`, so it is **identical with the two motifs swapped** —
every anchor therefore contributes twice:

| Pair | sum n | a=24 | a=36 | a=48 | Total |
|---|---|---|---|---|---|
| AP1+ELK1 | 16 | 2x19 | 2x18 | 2x6 | **86** |
| AP1+RFX | 22 | 2x19 | 2x12 | - | **62** |
| ELK1+RFX | 24 | 2x19 | 2x10 | - | **58** |

The ordering 86 > 62 > 58 is purely total motif length: longer motifs hit the
end <= 74 cap sooner, and at a=48 only the shortest pair survives.

### Triplet: 270

Factors exactly as **6 left-to-right orderings x 45 position configurations**,
verified for all six orderings.

### Total

| Component | Count |
|---|---|
| singles 4+4+4 | 12 |
| pairs 86+62+58 | 206 |
| triplet | 270 |
| **per background** | **488** |
| x 2 backgrounds | **976** |
| x 6 barcodes | **5,856** |

---

## 4. Running MPRAnator's own code

The source is Python 2 (`xrange`, `print` statements). It was translated mechanically
with `2to3 -w -n` — no logic edits — and driven directly:

```python
subsets = part1.generateCombinations(["TGACTCA", "ACCGGAAGT", "CGTTGCTAGGCAACG"])
finalOutput = []
for header, bg in BACKGROUNDS.items():
    for s in subsets:
        finalOutput += oligo.oligo(bg, 6, 24, list(s), 15, 15, 12, header)
barCodes, nbc = part1.getBarCodes(15, 40, 60, 6, 3, finalOutput)
response, _   = part1.createMPRAResultOutput(finalOutput, nbc, barCodes,
                                             "", "", "", "", ["Background", "Barcode"])
```

| Stage | Result |
|---|---|
| `generateCombinations` | 7 subsets |
| `oligo()` over 2 backgrounds x 7 subsets | **976 arrangements** |
| `getBarCodes(15, 40, 60, 6, 3)` | 5,856 barcodes, all 15 nt, all distinct |
| `createMPRAResultOutput` | **5,856 sequences**, all 105 nt, all distinct, 976 distinct arrangements |

**5,856 reproduced from their code, matching their documentation exactly.**

### Three defects found in their code while doing so

1. **The documented design cannot complete as specified.** With either restriction
   site set, `createMPRAResultOutput` calls
   `myf.findMatch(sequenceS=outputSequence.lower(), ...)`, which reaches
   `complement()`, whose `complementD` has **uppercase keys only** ->
   `KeyError: 'c'`. The use case specifies `CACGTG` and `CAATTG`, so the flagship
   example crashes in its own pipeline. The run above therefore omits the
   restriction sites, giving 105 nt (90 + 15) rather than the documented 117 nt tile.
   This is a concrete, reproducible cause for the HTTP 500 responses observed on the
   live service.
2. **`revcompl` does not reverse.** `myfunctions.revcompl` returns
   `complement(backgroundS[::1])` — step **1**, a no-op slice, not `[::-1]`. It
   computes a complement, not a reverse complement.
3. **The barcode GC constraint is never enforced.** Barcodes come back at
   **GC 6.7-93.3%** against a documented 40-60%. `barcode_generator`'s repair loop
   tests `gc_cont(...) * 100 > maxgc == True`, which Python chains as
   `(gc > maxgc) and (maxgc == True)`; with `maxgc = 60`, `60 == True` is `False`, so
   the clause is dead. `countlist` is also never reset between barcodes.

Only defect 1 affects the comparison, and only by changing the tile length. The
**976 arrangements are unaffected** — that stage runs clean.

---

## 5. PoolParty reproduction

`repro.py` in this directory. The structural mapping:

| MPRAnator | PoolParty |
|---|---|
| substitute, keeping length | `insertion_multiscan(..., replace=True)` |
| `generateCombinations` subsets | `itertools.combinations` over the motif dict |
| gap in [6, 24] | `min_spacing=6, max_spacing=24` |
| `sorted(pos)[0] % 12 == 0` | `insertion_mode="ordered"` (positions strictly increasing) + first insert's candidates restricted to multiples of 12 |
| motif permutations, from `itertools.product` over ordered tuples | the `k!` orderings, built as separate Pools and `stack`ed |
| overflow filter | automatic — PoolParty bounds each insert's positions by its own length |

```python
for k in (1, 2, 3):
    for subset in itertools.combinations(MOT, k):
        for order in itertools.permutations(subset):
            seqs = [MOT[n] for n in order]
            positions = []
            for i, s in enumerate(seqs):
                cand = range(EDGE, RIGHT + 1 - len(s))
                positions.append([p for p in cand if p % INTERVAL == 0] if i == 0
                                 else list(cand))
            pools.append(pp.insertion_multiscan(
                bg, k, [pp.from_seq(s) for s in seqs],
                positions=positions, replace=True,
                insertion_mode="ordered", mode="sequential",
                min_spacing=MIN, max_spacing=MAX,          # k > 1 only
                names=[f"{'_'.join(order)}_{i}" for i in range(k)]))
lib = pp.stack(pools)
```

**The `names=` argument is required, not optional.** Region lengths are Party-global,
so without unique names the second call collides:
`Region '_rep_0' already registered with seq_length=7, cannot re-register with
seq_length=9`. `insertion_multiscan` exposes `names=` for exactly this;
`deletion_scan` does not, which is `../comparison/FINDINGS.md` **B3** — the fix that
entry asks for already exists elsewhere in the package.

### Results

| | MPRAnator | PoolParty |
|---|---|---|
| bg1 arrangements | 488 | **488** |
| bg2 arrangements | 488 | **488** |
| Combined distinct | 976 | **976** |
| **Exact sequence overlap** | | **976 / 976 (100.0%)** |
| Unique to either tool | | **0** |

First attempt, no tuning. Contrast the VaLiAnT reproduction, which needed one codon
table aligned before reaching 100%: here there is no codon choice to disagree about,
because motifs are supplied literally.

### Independent confirmation of the enumeration

`insertion_mode="unordered"` with unrestricted per-insert ranges gives **1,140** for
AP1+ELK1, and MPRAnator's algorithm with filter 4 removed gives **1,140**. Two
independent implementations agreeing on an unconstrained count is strong evidence the
position/spacing semantics match, separately from the anchored counts.

---

## 6. What this does and does not show

**Shows.** PoolParty expresses MPRAnator's full combinatorial motif-placement design
— all subsets, all orderings, all admissible positions under edge, spacing and
interval constraints — and produces the identical 976 sequences.

**Does not show.** The **barcode stage is not compared.** MPRAnator's barcodes are
randomly generated with no published seed, its GC constraint is not enforced
(defect 3), and the live endpoint returns HTTP 500. PoolParty's `get_barcodes` would
produce a different but equally valid set. The comparison is therefore at the
arrangement level — 976 sequences — with the x6 barcode multiplier taken as
documented rather than verified.

**Capability PoolParty has and MPRAnator does not.** Orientation. MPRAnator has a
single **global** scaffold-level `reverseComplement` checkbox, not per-motif
orientation enumeration — it cannot place a motif in both orientations. Our own MPRA
example uses `flip` to do exactly that. Both tools score full marks on
`combinatorial_multi_motif_placement` in Table 2, but only one enumerates
orientation.

---

## 7. Status of this series

| Tool | Target | Result |
|---|---|---|
| **VaLiAnT** | `brca1_pep`, BRCA1 exons 2-5 SGE | **2,339 / 2,339 exact** — see `valiant/comparison.md` |
| **MPRAnator** | MPRA Motif Use Case | **976 / 976 exact** (arrangement level) |
| **tangermeme** | Tutorial B3 — motif-pair spacing sweeps | not started; pip-installable, so runnable head-to-head |

The mirror direction, for the same three tools: **VaLiAnT expresses one of the four
components of our GB1 library** (1,045 of 547,230 sequences, 0.19%) and neither of
our other two examples. **MPRAnator** cannot express the GB1 library at all — it has
no codon or amino-acid concept.
