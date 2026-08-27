# Plan: ORF-aware insertion / deletion scans for PoolParty v2

**Status:** PR #17 was merged into `origin/main` at `7b990ae` and the Phase 2 worktree
`/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-orf-scans`, branch
`feature/orf-indel-scans`, was rebased onto it. The branch contains five unpublished scan
commits, beginning with B3 at `bb102b7`; deletion is `6f94583` and insertion is `83162e3`.
Both wrappers passed adversarial, real-world interaction, and simplicity/robustness review.
The approved post-rebase helper consolidation (`complete_codon_count` from `_frame.py`) is
currently an uncommitted three-file source diff for review. The focused insertion/deletion
suite passes at **94 passed**, the combined ORF-operation suite at **319 passed**, and the
full suite at **3107 passed, 14 xfailed**.

**Audience:** a reviewer with no prior context on this work. Everything needed to evaluate
the plan is in this document. Section 13 lists the specific questions on which review input
is wanted.

**Provenance:** this plan follows a survey (`SURVEY.md`, same directory, 750 lines) which
established the factual base by execution. Where this document states a fact as *verified*,
it was produced by running code, not by reading it. Section 12 gives reproduction commands.

---

## 1. Orientation

### 1.1 The two codebases

| | path | role |
|---|---|---|
| **v1** | `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty` | Previous generation. Class-per-operation (`Pool` subclasses). Public history. |
| **v2** | `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty` | Current. Functions returning `Pool`s, composed into a DAG. `pyproject.toml` says 0.1.1; `poolparty.__version__` still says 0.1.0 (a separate metadata mismatch). |

Both import as `poolparty`, so they cannot be loaded in one process. v1 requires
`sys.path.insert(0, '/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty')`; v2 is importable
with plain `python3` from anywhere.

The pre-fix v2 baseline was **2915 passed, 14 xfailed**. The Phase 1 worktree currently
passes **2965 passed, 14 xfailed** with `PYTHONPATH=src python3 -m pytest tests/ -q` from
its `poolparty/` package directory. Historical defect measurements in Sections 3–4 refer to
the pre-fix baseline; Phase 2 must be developed against the Phase 1 branch or its merge.

### 1.2 The relevant v2 files

```
src/poolparty/
  scan_ops/deletion_scan.py        the non-ORF deletion scan
  scan_ops/insertion_scan.py       the non-ORF insertion scan (+ replacement_scan)
  region_ops/region_scan.py        both of the above are built on this
  orf_ops/annotate_orf.py          creates an OrfRegion (stores frame)
  orf_ops/mutagenize_orf.py        frame-aware
  orf_ops/stylize_orf.py           frame-aware
  orf_ops/translate.py             frame-aware
  region.py                        Region, OrfRegion (frame lives here)
  party.py                         register_region, upgrade_to_orf_region
  codon_table.py                   UNIFORM_MUTATION_TYPES
  utils/seq_utils.py               validate_positions
  utils/scan_utils.py              build_scan_cache
```

### 1.3 v2 design constraints any proposal must respect

1. **Lazy state decomposition.** `num_states` is a *product* of per-operation state counts;
   a backward pass decomposes a global state index into per-operation indices. Any design
   whose state count is not a computable product breaks the core mechanism.
2. **Uniform branching factors.** `codon_table.py:UNIFORM_MUTATION_TYPES` maps a mutation
   type to a *fixed* number of alternatives per codon (`any_codon` 63, `nonsynonymous_first`
   20, `missense_only_first` 19, `nonsense` 3). `synonymous` is deliberately absent because
   its count varies by residue, and `mutagenize_orf` raises for `mode='sequential'` with a
   non-uniform type. A codon-level indel scan has a uniform branching factor of 1 per
   position, so it fits.
3. **Operations are functions**, not classes. A port that reintroduces `Pool` subclasses is
   against the idiom.
4. **Prefer composition.** `deletion_scan` is itself implemented as
   `region_scan` + `replace_region`. Extend that rather than writing a parallel path.
5. Operations return `(Seq, dict)` — sequence plus design card.
6. Regions are registered **Party-globally** with a fixed `seq_length`.

---

## 2. The question, and the answer

**Question as posed.** v1 shipped `DeletionScanORFPool` and `InsertionScanORFPool` —
codon-level, frame-preserving indel scans. v2 has `deletion_scan` / `insertion_scan`, which
are not ORF-aware. Should the v1 capability be migrated?

**Answer.** Not as a port. v2 can already express every v1 behaviour tested — verified, ten
of ten cases byte-identical. But the generic scans ignore `OrfRegion.frame`, so callers must
manually derive codon starts and can silently target the wrong frame. The recommendation is
two compositional wrapper functions with codon-unit positions, explicit strand semantics,
and ORF-level cards, plus documentation.

The frame inconsistency discovered during design was a blocking prerequisite. Phase 1 has
resolved it and is included in the Phase 2 base (Sections 4–5). The internal-marker
collision in Section 7 is isolated in the first Phase 2 commit, before either wrapper.

---

## 3. Evidence base

### 3.1 v1 behaviour (source + 153 passing tests)

v1 files: `poolparty/deletion_scan_orf_pool.py` (506 lines),
`poolparty/insertion_scan_orf_pool.py` (553), `poolparty/orf_pool.py` (415).
Tests: `tests/test_deletion_scan_orf_pool.py` (1203 lines),
`tests/test_insertion_scan_orf_pool.py` (1495). Both suites pass: **153 passed in 0.76s**.

**Frame establishment.** v1 has no `frame` parameter. Frame is implied by nucleotide
indices `orf_start` / `orf_end`; `upstream_flank = seq[:orf_start]`, and the ORF span must
be divisible by 3 (hard error otherwise). There is no strand concept — v1 is plus-strand
only. **Consequence: v1 provides no reverse-strand behaviour to port.**

**`DeletionScanORFPool`** — all position parameters in *codon* units:

| parameter | default | semantics |
|---|---|---|
| `deletion_size` | required | codons deleted per variant |
| `start`/`end`/`step_size` | `0`/`num_codons`/`1` | range interface; `positions = range(start, min(end,L) - W + 1, step_size)` |
| `positions` | `None` | explicit codon positions; mutually exclusive with the range interface |
| `position_weights` | `None` | sampling weights, `mode='random'` only |
| `mark_changes` | `True` | `True` → replace each codon with `deletion_character*3` (length preserved); `False` → excise (length shortens by `3*W`) |
| `deletion_character` | `'-'` | |
| `orf_start`/`orf_end` | `0`/`len(seq)` | nucleotide bounds |

`num_internal_states = len(positions)`.

**`InsertionScanORFPool`** — same range/position interfaces, plus:

| parameter | default | semantics |
|---|---|---|
| `insert_seq` | required | must be ACGT and `len % 3 == 0` (**enforced**) |
| `insert_or_overwrite` | `'overwrite'` | `'overwrite'` → replace `W` codons, `L-W+1` positions; `'insert'` → splice between codons, **`L+1`** positions, independent of insert size |
| `mark_changes` | `False` | `True` → `swapcase()` the inserted codons |

**Dead code in v1:** the `del_aa` and `insert_aa` design-card fields are always `None`,
because both classes pass `codon_table=None` to `ORFPool.__init__`. They do not need porting.

### 3.2 What v2 can do today (verified by execution)

**Mechanism.** `utils/seq_utils.py:validate_positions` accepts
`Sequence[Integral] | slice | None` and resolves a slice via `positions.indices(n)`,
**honouring the step**. `utils/scan_utils.py:build_scan_cache` then sets
`num_states = len(resolved_positions)`. Positions are indices into the list of valid start
positions, and are **region-relative** when `region=` is supplied.

**The working idiom.** With an annotated region, `region=` plus a step-3 slice gives exact
codon-level behaviour:

```python
pp.deletion_scan(p0, deletion_length=3, deletion_marker=None,
                 region="orf", positions=slice(0, None, 3), mode="sequential")
```

**Ten of ten v1 cases reproduce byte-identically** (`SURVEY.md` §2.3):

| case | v1 call | v2 equivalent (all `region="orf"`, `mode="sequential"`) | n | identical |
|---|---|---|---|---|
| D1 | `deletion_size=1, mark_changes=True` | `deletion_length=3, deletion_marker='-', positions=slice(0,None,3)` | 6 | YES |
| D2 | `deletion_size=1, mark_changes=False` | `deletion_length=3, deletion_marker=None, positions=slice(0,None,3)` | 6 | YES |
| D3 | `deletion_size=2, mark_changes=False` | `deletion_length=6, deletion_marker=None, positions=slice(0,None,3)` | 5 | YES |
| D4 | `deletion_size=1, step_size=2` | `deletion_length=3, positions=slice(0,None,6)` | 3 | YES |
| D5 | `deletion_size=1, start=2, end=5` | `deletion_length=3, positions=slice(6,15,3)` | 3 | YES |
| D6 | `deletion_size=1, positions=[1,3]` | `deletion_length=3, deletion_marker='-', positions=[3,9]` | 2 | YES |
| I1 | `insert_seq="TAG", 'insert'` | `insertion_scan("TAG", positions=slice(0,None,3), replace=False)` | 7 | YES |
| I2 | `insert_seq="TAG", 'overwrite'` | `insertion_scan("TAG", positions=slice(0,None,3), replace=True)` | 6 | YES |
| I3 | `insert_seq="TAGTAG", 'overwrite'` | `insertion_scan("TAGTAG", positions=slice(0,None,3), replace=True)` | 5 | YES |
| I4 | `insert_seq="TAG", mark_changes=True` | `insertion_scan("tag", positions=slice(0,None,3), replace=False)` | 7 | YES |

Note I1: `replace=False` yields `n_codons + 1` positions regardless of insert size,
matching v1's `'insert'` mode exactly. I4: v1's `swapcase` marking is reproduced by passing
a lowercase insert.

**v2 also exceeds v1 here.** The insert may be a multi-state `Pool`, giving a true
position × insert product with structured naming — verified 18 states (6 positions × 3
inserts) with `prefix_position`/`prefix_insert`.

### 3.3 The three real gaps (verified)

**Gap 1 — the frame is annotated, validated, and then ignored.** `deletion_scan` output is
bit-identical across `frame=+1/+3/-1/-3`. Passing `region="orf"` to `deletion_scan` treats
the OrfRegion as a plain region: 16 nucleotide-level states on an 18-nt ORF, not 6.
Only three operations consume `OrfRegion.frame` — `translate`, `mutagenize_orf`,
`stylize_orf` — each via its own private `_resolve_frame`. **This is what makes
documentation-only insufficient.**

**Gap 2 — the bare-slice idiom silently overruns the ORF.** Without `region=`, a step-3
slice runs off the 3' end whenever a downstream flank exists. Measured on a 6-codon ORF
with a 4-nt trailer:

```
deletion_length=3, positions=slice(2,None,3)  ->  7 states (should be 6); state 6 deletes trailer bases
deletion_length=6, positions=slice(2,None,3)  ->  6 states (should be 5)
```

The correct bare-slice stop is `orf_end - deletion_length + 1`, **not** `orf_end` — the
intuitive choice is wrong for any multi-codon deletion. Passing `region=` removes the
arithmetic entirely. Any documentation of this idiom must present `region=` as mandatory.

**Gap 3 — `insertion_scan` accepts frameshifting inserts.** v1 raised on
`len(insert) % 3 != 0`; v2 does not:

```
pp.insertion_scan(p0, "TA", region="orf", positions=slice(0,None,3), replace=False)
  -> accepted, 7 states, first = GGGGGTAATGCCCGGGTTTAAACCCTTTTT   (frameshift, silent)
```

**Not a gap:** XML region tags leak into `generate_library()` output but are stripped by
`to_df()` / `to_file()` (via `write_tags`, default `False`). `generate_library` has no
`write_tags` parameter. Cosmetic, out of scope.

### 3.4 The only true capability gap

`position_weights` — v1's weighted random sampling over positions. v2 has no equivalent.
**Deliberately out of scope** (Section 9).

---

## 4. Historical blocking prerequisite: incompatible definitions of `frame`

This was found while designing the wrappers and had to be resolved first. The measurements
below describe the pre-fix baseline. Phase 1 is now implemented on the dedicated branch;
Phase 2 must use its shared convention rather than reproduce any arithmetic locally.

### 4.1 The defect

A reading frame says where codon boundaries fall. v2's three frame-aware operations do not
agree:

| file | line | expression | bases skipped |
|---|---|---|---|
| `orf_ops/translate.py` | 213 (also 161) | `frame_offset = abs(frame) - 1` | `\|f\|-1` |
| `orf_ops/mutagenize_orf.py` | 223 | `(4 - abs(frame)) % 3` | `(4-\|f\|)%3` |
| `orf_ops/stylize_orf.py` | 173, 260 | `region_frame = abs(frame) - 1`, applied as an index **shift** | `(4-\|f\|)%3` |

`stylize_orf` uses the same numeric value as `translate` but applies it as a grid shift
rather than a skip, which lands it on `mutagenize_orf`'s convention. This is the detail that
makes the bug hard to see by reading.

| frame | translate | mutagenize_orf | stylize_orf |
|---|---|---|---|
| ±1 | 0 | 0 | 0 |
| ±2 | **1** | **2** | **2** |
| ±3 | **2** | **1** | **1** |

**Frames ±2 and ±3 are swapped between the two camps.** Negative frames are affected
identically, mirrored from the 3' end.

### 4.2 Minimal reproduction

```
ATGAA with frame=2
A T G A A
0 1 2 3 4

mutagenize_orf operates on  GAA   (nt 2,3,4)
translate       reads       TGA   (nt 1,2,3)
```

Same sequence, same annotation, disjoint codons.

### 4.3 Observable damage

`mutagenize_orf(mutation_type="nonsense")` installs a stop codon at a boundary `translate`
does not read. Variants gaining a stop absent from the wild type:

| frame | +1 | +2 | +3 | -1 | -2 | -3 |
|---|---|---|---|---|---|---|
| gained a new stop | 15/18 | **0/15** | **0/18** | 18/18 | **0/15** | **0/18** |

At frames ±2/±3, the mutation introduces **no new stop at the mutated residue** when
translated by the package's own translator. (Some wild-type translations, notably the
chosen -3 fixture, already contain a stop elsewhere, so "no stop codons at all" would be
false.) Re-translating at the swapped frame restores 100% of eligible states — the mutation
logic is correct; only the offset disagrees.

Secondary, independent defects in `stylize_orf`:
- **Codon numbering** diverges from `mutagenize_orf` even where boundaries agree:
  `stylize_orf` counts the leading partial stub as codon 0, so `mutagenize_orf`'s codon *i*
  receives `style_codons[(i+1) % n]` at frames ±2/±3.
- `stylize_orf.py:297-301` derives the style *group* from an unshifted index but the
  position-within-codon from a shifted one, so with 6 styles a codon straddles two groups.
- `stylize_orf` → `translate(preserve_codon_styles=True)` yields 6 styled residues at
  frame +1 and **0** at ±2/±3 — codon colours silently vanish off-frame.

### 4.4 Why the test suite did not catch it

All 2915 tests pass with the bug present. Both parametrized frame tests in
`tests/test_audit_orf_ops_bravo.py` **re-derive the implementation's own formula in the test
body** (`frame_offset = frame - 1` for translate; `(4 - frame) % 3` for mutagenize), so they
cannot disagree with the code. Of 15 test functions using two or more frame-aware
operations, every one is at `frame=1`, where all three coincide. And
`tests/test_mutagenize_orf.py:176` carries a comment in `translate`'s convention while the
code computes `mutagenize_orf`'s — the two coincide for that 5 bp input, so the assertion
cannot discriminate.

### 4.5 Chosen resolution

**Keep `translate`. Change `mutagenize_orf` and `stylize_orf`.** Canonical:
`skip = abs(frame) - 1`.

Rationale:
- `region.py:84` — the docstring on the shared `OrfRegion` object all three read — says
  "the absolute value indicates the frame offset (1-indexed)", matching `translate`.
- `+1/+2/+3` conventionally means starting translation at bases 1/2/3. This is what outside
  bioinformatics tooling means and what users will expect.
- One `OrfRegion.frame` must identify one set of codons for translation, mutation and styling.
- The alternative concept (`mutagenize_orf`'s reading: "base 1 sits at codon position N",
  i.e. GFF *phase*) is legitimate but different, and if needed should be exposed separately
  as `phase`, not overloaded onto `frame`. Not needed now.

**Orphan-base policy** (bases belonging to no complete codon, at either end):

| operation | policy |
|---|---|
| `translate` | ignore (already does) |
| `mutagenize_orf` | preserve in DNA, assign no codon position (already does) |
| `stylize_orf` | **leave unstyled** (change) |

### 4.6 Measured blast radius

- **`mutagenize_orf` offset fix: 0 test failures.** Applied to a scratch copy →
  2915 passed, 14 xfailed, identical to baseline. The suite offers no protection at all here.
- **`stylize_orf` rewrite: 4 test failures**, all encoding behavior intentionally changed by
  the frame/orphan contract:
  `TestStylizeOrfBasic::test_style_codons_basic` plus
  `TestStylizeOrfFrame::{test_frame_2, test_frame_3, test_frame_affects_style_codons}`.
  The first is the frame +1 trailing-orphan case; all four expectations must be rewritten.
- **Library size changes, not only sequence content.** On a 19-nt region:
  `frame=2` offset 2→1, codons 5→6, `num_states` 315→378; `frame=3` the reverse.
  This happens **only when the region length ≡ 1 (mod 3)**:

  | region length | codons at skip 1 | at skip 2 | size changes? |
  |---|---|---|---|
  | L ≡ 0 (mod 3) | 18 → 5 | 5 | no |
  | L ≡ 1 (mod 3) | 19 → 6 | 5 | **yes** |
  | L ≡ 2 (mod 3) | 20 → 6 | 6 | no |

  For two-thirds of lengths the count is unchanged and only the sequences move — the easier
  case to miss when diffing old against new output.
- `revision/code_comparison/repro_valiant_brca1_pep.py:67` uses `frame=1`, so the frame
  offset change should not affect it. The full 2,339-state VaLiAnT reproduction was **not
  rerun during this review** and must not be presented as newly verified evidence.
- The old `.rst` pages did not pin down the offset. Phase 1 nevertheless updates the public
  ORF pages so the now-canonical convention is explicit rather than merely implicit.

---

## 5. Phase 1 — frame convention fix (implemented, not merged)

### 5.1 New file: `src/poolparty/orf_ops/_frame.py`

```python
"""Shared reading-frame helpers for ORF operations."""

from ..party import get_active_party
from ..region import VALID_FRAMES, OrfRegion
from ..types import Optional, RegionType


def frame_offset(frame: int) -> int:
    """Number of bases skipped before the first complete codon.

    frame=+N skips N-1 bases at the 5' end of the region.
    frame=-N skips N-1 bases at the 3' end (reading 3'->5').
    """
    return abs(frame) - 1


def resolve_frame(region: RegionType, frame: Optional[int]) -> int:
    ...  # moved verbatim from mutagenize_orf.py:21-56
```

`resolve_frame` consolidates the three existing private copies. `mutagenize_orf` and
`translate` are AST-identical (they differ by three stripped comments only);
`stylize_orf` differs solely by two vestigial function-local
`from ..party import get_active_party` statements at lines 22 and 142.

The original sketch omitted the required module-scope `get_active_party` import. The
implemented `_frame.py` includes it. Tests that imported private `_resolve_frame` symbols
were migrated in the extraction commit to import the shared `resolve_frame`; the audit test
also asserts that all three ORF modules reference the same function object.

**The extraction is verified safe.** `party.py` imports nothing from `orf_ops` at module or
deferred scope, so there is no cycle; `mutagenize_orf.py:12` and `translate.py:7` already
import that symbol at module scope. Hoisting both deferred sites in a scratch copy:
package imports, and the suite stays at 2915 passed / 14 xfailed.

### 5.2 `orf_ops/mutagenize_orf.py`

| location | change |
|---|---|
| 219–223 | `self.frame_offset = (4 - abs(frame)) % 3` → `frame_offset(frame)`; delete the now-wrong 3-line comment |
| ~101 | docstring: "the frame of the boundary base (1-indexed)" → "the first complete codon begins at base \|frame\| of the region" |
| 21–56 | delete `_resolve_frame`, import from `._frame` |

### 5.3 `orf_ops/stylize_orf.py`

Delete `self.region_frame` (line 173) entirely — it is the root of both bugs. Rewrite both
style methods so every index derives from one adjusted sequence:

```python
def _compute_codon_styles(self, molecular_positions):
    if len(molecular_positions) == 0:
        return []

    positions = molecular_positions[::-1] if self.reverse else molecular_positions

    # Drop orphans: frame offset at the start, incomplete codon at the end
    positions = positions[frame_offset(self.frame):]
    positions = positions[: (len(positions) // 3) * 3]

    num_styles = len(self.style_codons)
    style_positions = {s: [] for s in self.style_codons}
    for idx, pos in enumerate(positions):
        codon_index = idx // 3
        style_positions[self.style_codons[codon_index % num_styles]].append(pos)

    return [
        (s, np.array(sorted(style_positions[s]), dtype=np.int64))
        for s in self.style_codons
        if style_positions[s]
    ]
```

```python
def _compute_frame_styles(self, molecular_positions):
    if len(molecular_positions) == 0:
        return []

    positions = molecular_positions[::-1] if self.reverse else molecular_positions
    positions = positions[frame_offset(self.frame):]
    positions = positions[: (len(positions) // 3) * 3]

    num_style_groups = len(self.style_frames) // 3
    style_positions: dict[str, list[int]] = {}
    for idx, pos in enumerate(positions):
        codon_index = idx // 3           # both from the same idx, so a codon
        position_in_codon = idx % 3      # can no longer straddle two groups
        style = self.style_frames[(codon_index % num_style_groups) * 3 + position_in_codon]
        style_positions.setdefault(style, []).append(pos)
    ...
```

This closes four defects at once: the frame convention, the codon-numbering off-by-one, the
`style_frames` group/position mismatch, and the lost codon colours off-frame.
Also delete the two vestigial imports (lines 22, 142) and fix the docstring at line 89,
which currently documents `mutagenize_orf`'s convention while the code implements
`translate`'s.

**Note:** the orphan policy changes `stylize_orf` output at **frame ±1 too** — a trailing
partial codon is styled today and will not be. `translate` and `mutagenize_orf` already
ignore orphans and are genuinely unaffected at ±1.

### 5.4 `orf_ops/translate.py`

Lines 161 and 213: `abs(frame) - 1` → `frame_offset(frame)`. **Value unchanged**; the point
is that all three callers share one definition. Delete its `_resolve_frame` copy.

### 5.5 `region.py:84`

```
- The absolute value indicates the frame offset (1-indexed).
+ The absolute value indicates where the first complete codon begins:
+ frame=1 starts at the region's first base, frame=2 at its second,
+ frame=3 at its third. For negative frames the same offset is applied
+ from the 3' end, reading 3'->5'.
```

### 5.6 Tests

Implemented coverage includes:

- rewriting all **4** affected `stylize_orf` tests, including the frame +1 trailing orphan;
- replacing the misleading 5-bp `mutagenize_orf` comment with a 7-bp count discriminator;
- `tests/test_orf_frame_consistency.py`, using a region length ≡ 1 (mod 3) so codon
  *counts* differ between offsets — the case that discriminates;
- literal anchors for all six frames, cross-operation agreement, complete-product
  translation anchors, orphan handling, and end-to-end nonsense mutagenesis.

```python
REGION = "AATGCCCGGGTTTAAACCC"        # 19 nt, 19 % 3 == 1
FULL   = "GGGGG" + REGION + "TTTTT"
EXTENT = (5, 24)

# Hand-derived. NOT computed from any implementation formula.
EXPECTED_CODON0 = {
     1: (0, 1, 2),    2: (1, 2, 3),    3: (2, 3, 4),
    -1: (16, 17, 18), -2: (15, 16, 17), -3: (14, 15, 16),
}
```

The core assertions are all needed:
1. **Literal anchor** — each operation's codon 0 equals `EXPECTED_CODON0[frame]`. Catches
   all three drifting *together*, which is precisely how the existing tests fail.
2. **Cross-operation** — for all six frames, `translate`, `mutagenize_orf` and
   `stylize_orf` agree on codon 0's nucleotide indices.
3. **End-to-end** — annotate at each frame → `mutagenize_orf(mutation_type="nonsense")` →
   `translate` → assert a stop appears where the wild type has none. This is the original
   observable symptom.

**Known missing invariant:** no current test directly proves
`negative_op(S) == rc(positive_op(rc(S)))`, nor the corresponding card-coordinate mirror.
Section 10 makes those explicit pre-merge requirements rather than treating the existing
anchors as substitutes.

### 5.7 Commit sequence

| commit | contents | behaviour |
|---|---|---|
| `fdea858` | add `_frame.py`; route all three ops and affected tests through it | none |
| `be107e9` | add literal/cross-operation tests as `xfail` | none; documents the bug |
| `2ad79c8` | align `mutagenize_orf` | **changes** |
| `8104be3` | align `stylize_orf`; drop orphans | **changes** |
| `ebeede3` | add orphan and nonsense end-to-end tests | tests only |
| `80619ef` | document the convention and compatibility impact | docs only |
| `8791e36` | strengthen all-six-frame anchors and correct branch documentation | tests/docs |

The failing test precedes the behavior fixes in history, so it was not authored by copying
the patched output. The branch is green at **2965 passed, 14 xfailed**.

### 5.8 Release note — three separate statements

1. Frames ±2 and ±3 are swapped for `mutagenize_orf` and `stylize_orf`. Old `+2` ≡ new `+3`,
   old `+3` ≡ new `+2`, and the same for `-2`/`-3`.
2. Variant **counts** also change, but only when the annotated region's length ≡ 1 (mod 3).
   Otherwise counts hold and only sequences move.
3. `stylize_orf` no longer styles orphan bases, **at every frame including ±1**.

---

## 6. Phase 2 — the two new operations

Both live in `src/poolparty/orf_ops/`, are exported from `poolparty/__init__.py`, and are
added to `pool_mixins/dna_mixin.py` as methods, matching their three siblings.

### 6.1 Why wrappers rather than `frame=` on the existing scan ops

The rejected alternative was adding `frame=` to `deletion_scan` / `insertion_scan`.
It fails on **units**: `positions` in those functions is nucleotide-indexed. A `frame=`
argument would either do nothing to `positions`, or silently re-interpret it as
codon-indexed — a units change invisible to `beartype` and to the reader.
`mutagenize_orf` already solved exactly this by naming its codon-indexed parameter
`codon_positions` rather than overloading `positions`. The wrappers follow that convention
instead of bending it, and they join `orf_ops/` where the other three frame-aware
operations live.

### 6.2 `deletion_scan_orf`

```python
@beartype
def deletion_scan_orf(
    pool: Union[Pool, str],
    deletion_codons: Integral,
    deletion_marker: Optional[str] = "-",
    codon_positions: PositionsType = None,
    region: RegionType = None,
    frame: Optional[int] = None,
    prefix: Optional[str] = None,
    mode: ModeType = "random",
    num_states: Optional[Integral] = None,
    style: Optional[str] = None,
    iter_order: Optional[Real] = None,
    cards: CardsType = None,
) -> Pool:
```

- `deletion_codons` — required number of codons removed per variant (v1
  `deletion_size`), matching generic `deletion_scan`'s required width argument.
- `deletion_marker` — `None` excises (v1 `mark_changes=False`); a character gap-marks,
  length-preserving (v1 `mark_changes=True`). Default `'-'` matches `deletion_scan`.
  A non-`None` marker must be exactly one character so marked ORF deletion remains
  length-preserving.
- `codon_positions` — **codon units**, `slice | Sequence[Integral] | None`; `None` = all.
- `frame` — `None` resolves from the `OrfRegion` via `resolve_frame`.
- Scan controls and prefix/style parameters forward to the composed operations. `cards`
  follows the ORF-level ownership contract in Section 6.6 rather than being passed blindly
  to `region_scan`.

### 6.3 `insertion_scan_orf`

```python
@beartype
def insertion_scan_orf(
    pool: Union[Pool, str],
    insertion_pool: Union[Pool, str],
    codon_positions: PositionsType = None,
    region: RegionType = None,
    frame: Optional[int] = None,
    replace: bool = False,
    style: Optional[str] = None,
    prefix: Optional[str] = None,
    prefix_position: Optional[str] = None,
    prefix_insert: Optional[str] = None,
    mode: ModeType = "random",
    num_states: Optional[Integral] = None,
    iter_order: Optional[Real] = None,
    cards: CardsType = None,
) -> Pool:
```

`replace=False` splices between codons (v1 `'insert'`); `replace=True` overwrites whole
codons (v1 `'overwrite'`). `insertion_pool` is interpreted in **coding orientation** for
both strands. Position counts:

| mode | eligible codon slots |
|---|---|
| `replace=False` | `n_codons + 1` — independent of insert size |
| `replace=True` | `n_codons - item_codons + 1`, where `item_codons = insertion_pool.seq_length // 3` |

A multi-state `insertion_pool` is passed straight through, giving the
position × insert product with `prefix_*` naming (verified working today).

For a negative frame, keep the background sequence in stored plus/reference orientation,
map the requested coding slot to its physical coordinate, and reverse-complement the
**entire selected insertion sequence** before insertion. Do not reverse-complement codons
individually without reversing their order. For example, coding insert `ATGGAA` is written
physically as `TTCCAT`, not `CATTTC`.

This choice intentionally differs from generic `insertion_scan`, whose insert is literal
plus-strand content. It matches the ORF abstraction: supplying `TAG` means "insert a stop
codon" on either strand. Users who want literal reference-strand insertion should use the
generic nucleotide operation.

### 6.4 Shared codon → nucleotide translation

```python
def _codon_slots_to_nt(codon_positions, *, span, frame, item_codons, splice):
    """Resolve codon-unit positions to region-relative nucleotide start positions."""
    offset   = frame_offset(frame)             # abs(frame) - 1 after Phase 1
    n_codons = (span - offset) // 3
    n_slots  = n_codons + 1 if splice else n_codons - item_codons + 1
    if n_slots <= 0:
        raise ValueError(...)
    codon_idx = validate_positions(codon_positions, max_position=n_slots - 1, min_position=0)

    if frame > 0:
        return [offset + 3 * k for k in codon_idx]
    # negative frame: coding block is right-aligned; codon 0 is the RIGHTMOST
    # complete codon, matching mutagenize_orf.
    right = span - offset
    if splice:
        return [right - 3 * k for k in codon_idx]
    return [right - 3 * (k + item_codons) for k in codon_idx]
```

**`num_states = len(codon_idx)`** — a plain count, so the product decomposition
(constraint 1) is untouched.

The arithmetic was independently derived and spot-checked against patched
`mutagenize_orf` for one- and two-codon negative windows. That is encouraging but not a
substitute for the general RC-equivalence tests specified in Section 10.

### 6.5 Negative-strand normalization contract

The whole construct remains in stored plus/reference orientation. Each operation normalizes
only the affected biological content:

| operation | negative-frame behavior |
|---|---|
| `annotate_orf` | records frame/tags/styles; molecular DNA bases are unchanged |
| `translate` | traverses physical codons right-to-left, reverse-complements each, translates |
| `mutagenize_orf` | uses coding-oriented WT/mutant codons; reverse-complements mutants before writing |
| `stylize_orf` | traverses positions in coding order and styles the corresponding physical bases |
| `deletion_scan_orf` | maps coding codons to a physical interval and deletes/marks it; no RC is needed |
| `insertion_scan_orf` | maps the coding slot/window and writes the reverse complement of the entire insert |

The strongest executable contract is normalization equivalence. For `N in {1,2,3}` and
the same coding-oriented arguments:

```text
negative_edit(S, frame=-N)
    == RC(positive_edit(RC(S), frame=+N))

translate(S, frame=-N)
    == translate(RC(S), frame=+N)
```

Comparisons use molecular/tag-stripped output where representation-only XML tags differ.
For styles, the negative result's physical style positions must be the mirror image of the
positive-on-RC result. Section 10 requires this property for current and new ORF operations.

### 6.6 Design cards

Cards separate **coding semantics** from **physical coordinates**:

| key | operations | exact meaning |
|---|---|---|
| `codon_positions` | deletion; overwrite insertion | tuple of affected 0-based codon indices in coding order; codon 0 is rightmost for negative frames |
| `codon_slot` | splice insertion | boundary index in coding order, `0..n_codons`; this is not a codon index |
| `wt_codons` | deletion; overwrite insertion | tuple of affected WT codons in coding orientation, aligned with `codon_positions` |
| `start` | both | inclusive physical coordinate in the operation's input region, stored plus/reference orientation |
| `end` | both | exclusive physical coordinate in that input region; equals `start` for a splice |

Do **not** expose the old proposed `codon_position_nt`: "codon start" is ambiguous on the
negative strand and for a splice. `start`/`end` retain v2's existing region-relative,
half-open scan convention. Do **not** rename generic `region_scan.position_index` to a codon
position: it is only the ordinal within the resolved list of allowed starts.

Do **not** duplicate insertion content as `insert_codons`. The insertion pool owns its
identity and provenance. Cards from an insertion pool already propagate through
`insertion_scan`, including a composed/mutagenized pool. Recommended construction:

```python
inserts = pp.from_seqs(
    ["TAG", "TAA", "TGA"],
    mode="sequential",
    cards={"seq": "coding_insert_seq", "seq_index": "insert_index"},
)
```

For a negative ORF, `coding_insert_seq` stays `TAG` while the stored plus-strand output
contains `CTA`. A literal string insert is constant call metadata and is not repeated in
every row by default; users wanting a self-contained export can wrap the one string in a
fixed one-state insertion pool with the same `seq` card.

The generic scan cards cannot manufacture the ORF fields above: `position_index` is an
ordinal, `region_seq` contains internal tags, and a splice marker is empty. Add an internal
fixed, one-state card-producing operation at the point where the selected position and WT
content are still available. Its state-space factor is one, preserving the product
invariant. It owns only the ORF edit cards; insertion-content cards remain upstream.

Card normalization is part of the RC equivalence contract. Direct-negative and
RC-normalized-positive runs must have identical `codon_positions`/`codon_slot` and
`wt_codons`. For an input-region length `L`, physical cards mirror as:

```text
interval [start, end)  <->  [L - end, L - start)
splice point p         <->  L - p
```

Coordinates always describe the operation's **input** region, before a length-changing
edit. v1's `del_aa` / `insert_aa` are not ported; they are always `None` in v1.

### 6.7 Validation owned by the wrappers

1. `insertion_pool.seq_length` must be known and divisible by 3; otherwise raise before
   construction. Every selected insert state must contain only `ACGTacgt` molecular bases;
   validate at compute time when content is state-dependent. This matches v1, which rejects
   non-ACGT content as well as non-triplet length.
2. For a negative frame, reverse-complement the entire selected insertion sequence once,
   including reversal of codon order; preserve the insertion pool's upstream coding cards.
3. `deletion_codons <= 0` or `> n_codons` → raise.
4. `deletion_marker` must be `None` or exactly one character; otherwise raise rather than
   silently changing the length of a marked codon window.
5. `replace=True` and `item_codons > n_codons` → raise.
6. `region` is a plain `Region` and `frame is None` → raise (inherited from `resolve_frame`).
7. Unknown region span when geometry must be computed → raise.
8. Orphan bases (`span - offset` not divisible by 3) are **allowed** and never deleted or
   overwritten, consistent with `translate` / `mutagenize_orf`. Documented, not an error.
   (v1 hard-errored, but v1 had no frames and no orphan concept.)

---

## 7. Implemented prerequisite commit — internal marker collisions (`SURVEY.md` B3)

Before `bb102b7`, `deletion_scan` hard-coded its internal region name `_del`. Region lengths
are Party-global, so two calls with different `deletion_length` in one Party collided:

```
Region '_del' already registered with seq_length=3, cannot re-register with seq_length=6
```

Replacement scans had the same defect because every replacement used `_rep`, whose
registered length is the insertion-pool length. Splice insertion uses `_ins` with a
zero-width marker, so ordinary splice calls do not have the length mismatch.

**Implemented fix:** encode the registered marker length in the existing internal
`tag_name` passed to `region_scan`:

```text
deletion length N      -> _del_lenN
replacement length N   -> _rep_lenN
splice insertion       -> _ins       (marker length is always zero)
```

This makes each internal name imply exactly one Party-global length. Calls of the same
length safely reuse one registration because each composed scan consumes its temporary
tags before returning. The spelling includes `len` so it does not collide with multiscan's
default positional names (`_del_0`, `_rep_0`, and so on).

No public parameter was added: users should not manage high-level scan implementation
tags. The two source changes are in the generic DNA scan operations, so the ORF wrappers
inherit the fix. Four new regression tests cover chained and branched different-length
calls, marked and true deletion, functional and mixin APIs, molecular output, registry
lengths, and the intentionally changed `name` / `region_seq` cards. Two read-only reviews
(one adversarial) approved the final diff. After rebasing onto PR #17, commit `bb102b7` is
the first commit on `feature/orf-indel-scans`; its original pre-rebase hash was `c461833`.

---

## 8. Phase 2 documentation

| file | change |
|---|---|
| `docs/operations/deletion_scan.rst` | `positions` is described as "Explicit list of window start positions" and never mentions slices. Document the `slice` form, and present `region=` as **required** for ORF work, with the overrun (Section 3.3) as a warning. |
| `docs/operations/insertion_scan.rst` | same; also show that insertion-pool cards (`seq`, `seq_index`, provenance cards) propagate to the final library |
| `docs/operations/deletion_scan_orf.rst` | new; define coding-oriented codon numbering, plus/reference `start`/`end`, negative normalization, and orphans |
| `docs/operations/insertion_scan_orf.rst` | new; define coding-oriented inserts, whole-insert RC on negative frames, splice slots, overwrite windows, and insertion-pool card ownership |

Precedent for documenting the slice form already exists:
`mutagenize_orf`'s `codon_positions=slice(1, 56)` appears in `docs/tutorials/dms_gb1.rst`.

---

## 9. Explicitly out of scope

| item | why |
|---|---|
| `position_weights` | v1's only true capability gap, but it is general scan *sampling*, not ORF semantics. Belongs in its own change, applied to `deletion_scan`/`insertion_scan` uniformly. |
| Redesigning generic scan cards | Existing `position_index`/`start`/`end` remain valid low-level cards. Internal `name`/`region_seq` leakage can be cleaned up separately; Phase 2 must not silently change their public meanings. |
| Porting v1's class hierarchy | Against v2's idiom (constraint 3). `ORFPool`, `DeletionScanORFPool`, `InsertionScanORFPool` are not reintroduced in any form. |
| A `phase` parameter | The discarded frame interpretation may be worth exposing separately one day (GFF import). No demonstrated need now. |
| XML tag leakage in `generate_library` | Cosmetic; `to_df`/`to_file` strip tags. |
| Nested annotations inside an edited ORF | Generic `region_scan` can cross nested XML boundaries. Phase 2 rejects these targets before scanning rather than inheriting malformed XML; a general annotation-overlap policy belongs in `region_scan`. |
| Repeated true deletions on one named region | Party-level named-region lengths are immutable and remain at the pre-edit span. Phase 2 documents this broader metadata limitation. Translation and runtime-geometry random mutagenesis remain supported; sequential geometry-dependent ORF operations on that shortened named region do not. Whole-sequence true-delete chains work. Marked deletion does not stale the length, but its gapped output remains outside the wrapper's input scope. |
| Copying prefixed insertion scans | The shared generic `PassthroughOp` naming path cannot currently be reconstructed by `copy()` / `deepcopy()` and captures original state references. Generic `insertion_scan` has the identical limitation. Unprefixed ORF insertion copies are tested; fix prefixed copies once in the generic naming operation rather than adding an ORF-only workaround. |
| v1 `repr` formats, `iteration_order`, Pool-chaining tests | v2 has different mechanisms. |

---

## 10. Test strategy

**Transfers from v1** (as behaviour specs, rewritten to v2 idiom): `TestMarkedDeletionMode`,
`TestUnmarkedDeletionMode`, `TestSequenceLength`, `TestEdgeCases`,
`TestMultiCodonOperations`, `TestReadingFrameIntegrity`, `TestFlankingRegions`,
`TestOverwriteMode`, `TestInsertMode`, `TestMarkChanges`.

**Does not transfer:** anything touching `position_weights` (out of scope), `TestRepr`,
`mode`/`iteration_order` internals, `TestPoolChaining`.

### 10.1 What Phase 1 already covers

`test_orf_frame_consistency.py` now has hand-derived all-six-frame anchors, whole-protein
translation expectations, cross-operation codon-0 agreement, orphan styling/mutation
coverage, and nonsense-mutagenesis end-to-end checks. These strongly pin the selected
frame convention.

They do **not** directly test reverse-complement normalization equivalence, and they do not
test card-coordinate mirroring. A source search found no existing ORF test asserting either
property. An ad hoc executed check on three non-palindromic sequences × all three offsets
passed for `translate` and for `mutagenize_orf` sequence plus
`codon_positions`/`wt_codons`/`mut_codons` cards, but it is not permanent regression
coverage; `stylize_orf` was not included. Add the following Phase 1 follow-up before merge:

1. `translate(S, frame=-N) == translate(RC(S), frame=+N)` for `N=1,2,3`, using
   non-palindromic regions and every length residue class modulo 3.
2. Direct-negative `mutagenize_orf` equals RC-normalized-positive mutagenesis for the same
   `codon_positions` and mutation state. Assert both molecular output and identical
   `codon_positions`, `wt_codons`, `mut_codons`, `wt_aas`, and `mut_aas` cards.
3. Direct-negative `stylize_orf` equals the mirrored styles from positive styling on the RC
   input, for both `style_codons` and multi-group `style_frames`.

These are additional invariants, not replacements for the literal anchors: two operations
could be RC-equivalent while sharing the same wrong offset.

### 10.2 Phase 2 behavior tests

Implementation status: all items below are covered by the committed 42-case deletion suite
and committed 52-case insertion suite. The insertion tests include all-six-frame
anchors, whole-insert reverse complementation, interval boundaries, multistate insertion
products and provenance, B3 coexistence, wrapper validation, public function/mixin exposure,
downstream translation, first/middle/last RC equivalence, style mirroring, and disabled-style
behavior. After rebasing onto PR #17, the combined insertion/deletion suite passes at
94 passed, the combined ORF-operation suite at 319 passed, and the full suite at
3107 passed / 14 xfailed.

1. The ten equivalence cases of Section 3.2 as regression fixtures — they are already
   byte-identical, so they can land on day one.
2. All six frames and all three offsets, with one- and two-codon deletion/overwrite windows.
3. Orphan-base handling: assert leading/trailing orphans are never deleted or overwritten.
4. Insertion validation: unknown length, non-triplet length, and non-`ACGTacgt` content all
   raise; include a multi-state pool whose later state is invalid so runtime validation is
   exercised.
5. **ORF-boundary non-overrun:** reproduce Section 3.3 and assert both flanks are
   byte-identical in every state.
6. Splice slots `0`, middle, and `n_codons`; overwrite windows at first and last eligible
   positions; one- and multi-codon inserts.
7. Multi-state insertion pools: assert the exact position × insert state count, literal
   sequences, names, and propagation of `coding_insert_seq`, index, and provenance cards.
8. B3 integration: two wrapper calls with different deletion/overwrite lengths in one Party
   succeed by using the length-qualified generic marker names from `bb102b7`.

### 10.3 Strongest negative-frame equivalence tests

For each `N in {1,2,3}`, compare a direct negative operation with the positive operation on
the reverse-complemented input. Use the same coding-oriented codon index/slot and insertion
pool state:

```text
negative_delete(S, -N, k, width)
    == RC(positive_delete(RC(S), +N, k, width))

negative_insert(S, -N, slot, coding_insert)
    == RC(positive_insert(RC(S), +N, slot, coding_insert))
```

Cover marked and true deletion, splice and overwrite, one/two-codon windows, first/middle/
last slots, orphan-containing regions, and multi-state insertion pools. Compare molecular
DNA after stripping representation-only tags. For insertions, this test proves that the
negative implementation reverse-complements the **whole insert**, including codon order.

### 10.4 Design-card equivalence tests

For every paired direct-negative / normalized-positive case:

- `codon_positions` or `codon_slot` must be identical;
- `wt_codons` must be identical and coding-oriented;
- insertion-pool `coding_insert_seq`, state/index, and provenance cards must be identical;
- no wrapper-generated `insert_codons` card is expected;
- for input-region length `L`, deletion/overwrite interval cards must satisfy
  `[s_neg, e_neg) == [L - e_pos, L - s_pos)`;
- splice coordinates must satisfy `p_neg == L - p_pos` and `start == end`;
- coordinates must refer to the pre-edit input region, including for true deletion and
  length-changing splice insertion.

Also test explicit, sliced, reordered, and sparse `codon_positions`. This prevents the
internal `region_scan.position_index` (an ordinal within the allowed-position list) from
being mislabeled as the actual codon index.

---

## 11. Risks

1. **Reverse-strand normalization (highest).** Codon 0 is the rightmost complete codon for
   negative frames; splice boundaries run in the opposite physical direction; a multi-codon
   insert must be reverse-complemented as one sequence. Hand examples are insufficient.
   The RC-equivalence and card-mirror tests in Sections 10.3–10.4 are merge requirements.
2. **Card-producing composition.** A pure wrapper around `region_scan` cannot derive truthful
   ORF cards. The fixed one-state card operation must see the selected position and WT
   content without duplicating state factors or leaking internal XML tags. Prototype this
   with sparse/reordered positions and multi-state insertion pools before freezing the API.
3. **Insertion-pool orientation and provenance.** On negative frames the physical insert is
   RC'd, while upstream cards intentionally remain coding-oriented. Docs and tests must call
   the field `coding_insert_seq`; styles and state identity must survive the transformation.
4. **Phase ordering.** Satisfied: the Phase 2 branch starts from merged Phase 1 and its first
   commit is the B3 prerequisite.
5. **Orphan policy for ORF regions** whose span is not `offset + 3n`. v1 hard-errored; the
   proposal allows it. If a user expects an exact CDS, silence may be the wrong default —
   an optional strict mode was suggested and deferred.
6. **Observable generic scan cards.** B3 intentionally changes the low-level `name` and
   `region_seq` card values from `_del` / `_rep` to deterministic length-qualified names.
   Regression tests pin this behavior; molecular output, state counts, and public
   signatures are unchanged.

---

## 12. Reproduction

```bash
V1=/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty
V2=/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-statecounter/poolparty
FRAME_WT=/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-frame-fix/poolparty

# historical pre-fix baseline
cd $V2 && python3 -m pytest tests/ -q                 # 2915 passed, 14 xfailed
cd $V1 && python3 -m pytest tests/test_deletion_scan_orf_pool.py \
                            tests/test_insertion_scan_orf_pool.py -q   # 153 passed

# implemented Phase 1 branch (verified 2026-08-18)
cd $FRAME_WT && PYTHONPATH=src python3 -m pytest tests/ -q  # 2965 passed, 14 xfailed

# the historical frame split on the pre-fix tree
grep -n "frame_offset = " $V2/src/poolparty/orf_ops/translate.py        # abs(frame) - 1
grep -n "self.frame_offset" $V2/src/poolparty/orf_ops/mutagenize_orf.py # (4 - abs(frame)) % 3
grep -n "self.region_frame" $V2/src/poolparty/orf_ops/stylize_orf.py    # abs(frame) - 1, as a shift

# no cycle through orf_ops (the _resolve_frame extraction is safe)
grep -n "orf_ops\|mutagenize_orf\|stylize_orf\|translate" $V2/src/poolparty/party.py   # no matches
```

Minimal frame contradiction:

```python
import poolparty as pp
with pp.Party():
    p0 = pp.annotate_orf(pp.from_seq("ATGAA"), "orf", frame=2)
    m = pp.mutagenize_orf(p0, region="orf", num_mutations=1,
                          mutation_type="any_codon", mode="sequential",
                          cards=['wt_codons'])
    print(m.to_df(num_cycles=1))                       # mutates GAA (nt 2,3,4)
with pp.Party():
    p0 = pp.annotate_orf(pp.from_seq("ATGAA"), "orf", frame=2)
    print(pp.generate_library(pp.translate(p0, region="orf"), num_cycles=1))  # '*', from TGA
```

To measure the Phase 1 blast radius without touching the repo: copy
`$V2/src/poolparty` to a scratch directory, edit the copy, and run
`PYTHONPATH=<scratch> python3 -m pytest tests/ -q` from `$V2`.

---

## 13. Resolved decisions and remaining questions

### 13.1 Resolved

1. Keep `translate`'s `abs(frame)-1` convention; Phase 1 implements it.
2. Allow and preserve orphan bases; never translate, mutate, style, delete, or overwrite
   them as codons.
3. Accept the frame ±1 `stylize_orf` orphan behavior change in Phase 1.
4. Support negative frames in both new wrappers.
5. Define ORF insert sequences in coding orientation; reverse-complement the entire selected
   insert before physical placement on a negative frame.
6. Define codon/card fields in coding orientation and `start`/`end` in stored plus/reference
   orientation, relative to the operation's input region.
7. Let the insertion pool own insert sequence/index/provenance cards. The wrapper does not
   duplicate `insert_codons`.
8. Make direct-negative versus RC-normalized-positive sequence and card equivalence a merge
   requirement; existing tests do not yet cover this property directly.
9. Keep B3 isolated as the first prerequisite commit on the Phase 2 branch; implemented in
   `bb102b7` (rebased from `c461833`) with no public naming parameter.
10. Keep `position_weights` and generic scan-card redesign out of scope.
11. Match generic `deletion_scan` structurally: require `deletion_codons`, retain
    `deletion_marker='-'`, and require a non-`None` marker to be one character.
12. Place a fixed one-state ORF card operation after `region_scan` marks the selected WT
    interval and before `replace_region` removes it. This operation sees the active physical
    position and WT sequence, owns only the four ORF deletion cards, and contributes a
    state-space factor of one. The deletion implementation verifies the state product with
    multi-state parents and preserves the operation through `copy()` / `deepcopy()`.
13. Reject nested annotations in the target ORF with a fixed one-state validation operation
    before `region_scan`; do not allow the inherited low-level tag-slicing defect to emit
    malformed XML. Keep the general annotation-overlap policy and immutable named-region
    length redesign out of Phase 2.

### 13.2 Insertion implementation decisions now resolved

1. **Insertion card operation reuse:** implemented by adapting the proven deletion
   composition while leaving sequence/index/provenance cards with `insertion_pool`. Tests
   prove the complete position × multi-state-insert state product and all three generic
   naming prefixes.
2. **State-dependent insertion validation:** fixed-length validation occurs at construction;
   molecular content is validated per selected state at compute time because random or very
   large pools cannot be exhaustively enumerated during construction. Tests exercise a later
   invalid state and require a clear error.
3. **Shared scan helpers:** deletion and insertion use one private `_scan.py` for region-span
   resolution, codon geometry, target extraction, and fixed one-state target validation.
   The deletion refactor is behavior-neutral: all 42 focused deletion tests still pass.
4. **Review gate:** satisfied. Adversarial, real-world interaction, and
   simplicity/robustness reviews approved the final implementation. The adversarial review
   found one negative-frame crash when inline styles were disabled; the fix and regression
   tests across frames -1/-2/-3 were re-reviewed and approved. The inherited generic
   prefix-plus-copy limitation is explicitly deferred above.

---

## 14. Effort

| phase | work | estimate |
|---|---|---|
| 1 | frame convention fix, 7 commits, 4 old tests rewritten, extensive new test file | implemented; RC-equivalence follow-up still required |
| B3 prerequisite | length-qualified internal names in generic deletion/replacement scans; 4 regression tests | implemented as `bb102b7`; full suite green |
| 2 | two wrappers, shared geometry, negative insertion orientation, fixed card operation | implemented; three insertion reviews approved |
| 2c | docs: 2 generic pages fixed, 2 ORF pages new | implemented |

The original one-to-one-and-a-half-day estimate was optimistic because the pure-wrapper
card promise was not implementable and the negative-strand contract lacked property tests.
Sequence/card RC equivalence and multi-state insertion-pool behavior are the parts that need
slow, careful work.
