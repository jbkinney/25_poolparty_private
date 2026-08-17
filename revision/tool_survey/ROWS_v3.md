# Capability matrix rows, v3 (post-verification)

Supersedes `ROWS_v2.md`. The eight argument-critical rows below were audited by
one max-effort agent per row across all 13 tools (one rater per row, so one
threshold per row). Audits with quoted primary-source evidence, pre-registered
operational tests and mandatory disconfirmation are in `verified/`.

Three rows were returned `reword` because the **rubric** was defective, not
because the cells were mis-scored. Those rewordings are adopted here.

## Verification outcome

| Row | Verdict | Cells changed |
|---|---|---|
| `library_first_class_object` | reword | 10 |
| `composable_operations` | reword | 4 |
| `lazy_generation` | reword -> replaced | 3 |
| `library_algebra` | keep | 3 |
| `exhaustive_single_scans` | keep | 2 |
| `sampled_random_mutagenesis` | keep | 1 |
| `higher_order_combinatorial` | keep | 1 |
| `heterogeneous_components_one_library` | keep | 1 |

25 changes across 8 rows. **PoolParty was demoted on one row** (`lazy_generation`,
yes -> partial) before that row was replaced; see below.

---

## Adopted rewordings

### 1. `library_first_class_object` — criteria (ii) and (iii) deleted

**Defect.** The v2 rubric required three things, but criterion (ii) *is* the
`library_algebra` row and criterion (iii) *is* the `lazy_generation` row. The
rubric re-merged exactly the claims `ROWS_v2.md:13` had split apart. Effects:
PoolParty won one architectural advantage three times (reads as padding);
Oligopool Calculator lost this row for a laziness deficit already scored
elsewhere; and it forced two false equivalences in opposite directions —
Oligopool == VaLiAnT (Oligopool documents a file-free in-memory chain,
`docs.md:140-152, 199-232, 1915`; VaLiAnT exports only `__version__`), and
DnaChisel == PoolParty (DnaChisel satisfies all three literally but optimises a
single sequence and its `MutationSpace` is a search space that is never shipped).

**Adopted wording.** "Does a documented API call return a value representing the
whole multi-sequence library as an instance of a type the tool defines for that
purpose (carrying member identity, design provenance, or membership operations),
which the user can bind to a variable and inspect — as opposed to the tool only
writing files, or handing back a general-purpose container (ndarray, DataFrame,
ragged array, torch tensor, list/dict of sequences) whose library meaning lives
only in the caller's head?"

- `yes` — such a purpose-built type exists and is documented
- `partial` — a durable, documented, in-memory library value exists but its type is general-purpose
- `no` — the library exists only as files, or as a transient return value with no library identity

**Why this is better for the paper, not worse:** it is the claim the authors
actually wrote, so it needs no defending; it ends the triple counting, so
PoolParty's distinctiveness becomes *three independent facts* — `yes` on this row
**and** `library_algebra` **and** `unmaterialized_library_addressing` — which
neither DnaChisel nor pydna can match (both fail `library_algebra`). Three
independent wins are more persuasive than one compound cell.

### 2. `composable_operations` — "no" split into two footnoted variants

**Defect.** v2 defined `no` as "the tool fixes the pipeline order". That is true
of VaLiAnT, MPRA Design Tools, Mutation Maker and CodonGenie. It is **false** of
pydna, Biopython and SeqPro, which fix no pipeline at all — they simply expose no
library-design operations to compose. Three of seven `no` cells asserted
something their own source contradicts, and referees who authored those packages
would be right to object.

**Adopted wording.** Label: *user-composable design operations*. Question: "Are
library-design operations first-class values the user can chain and branch in an
order of their choosing?"

- `yes` — >=2 distinct design operations share a carrier type and compose in either order, with fan-out expressible
- `partial` — design operations exist and chain, but only in a tool-fixed order, only for a small documented subset, or only by re-feeding output artefacts
- `no†` — the tool fixes the order of its design operations in source
- `no‡` — the tool exposes no library-design operations, so there is nothing to compose (it may still be highly composable over general sequence utilities)

Costs one legend line. The alternative (splitting into `design_ops_first_class` +
`composition_order_user_chosen`) is cleaner but adds a row and partly duplicates
the table's scope row — not adopted.

Note: the auditor deliberately **excluded shared-node memoisation** from the
threshold as self-serving to PoolParty. PoolParty earns its `yes` on
expressibility alone, verified behaviourally.

### 3. `lazy_generation` -> replaced by `unmaterialized_library_addressing`

**Defect.** Under one consistent threshold **no tool reached `yes`** — the column
collapsed to 5 partial / 8 no, and PoolParty tied four competitors. Worse, it
exposed the most attackable cell in the table: DnaChisel's
`MutationSpace.all_variants` is a **pure generator** (`yield` over
`itertools.product`) whereas PoolParty's equivalent is a **materialised list**
(`cache.append` over `combinations`). Publishing `poolparty=yes, dnachisel=partial`
would be indefensible to any referee who opened both files. The row also could
not distinguish lazy enumeration of a declared finite library from sampling an
unenumerated space (Ledidi has no N at all).

**Adopted replacement.** `unmaterialized_library_addressing` — "Can a library far
larger than memory be declared, its EXACT size reported, and an arbitrary
k-member slice generated AT ANY OFFSET, with work and peak memory O(k) and no
per-member structure of size N, through a documented parameter of the tool's
primary generation call?"

- **poolparty = yes.** `state.num_values == 27,036,009,000` (exact Python int); DAG construction 0.001 s / 0.0 MB; `generate_library(num_seqs=1000)` = 0.2198 s; `init_state=27036008990` returns members #27,036,008,990-992 in 0.0017 s.
  *Caveat to state in prose:* sequential `mutagenize` / `mutagenize_orf` / `recombine` / `region_multiscan`, and `get_barcodes` always, precompute their state index at construction.
- **dnachisel = partial.** `islice(all_variants(seq), k)` is O(k), but `space_size` is `np.exp(min(100, np.log(choices).sum()))` — a float **capped at e^100**, not an exact count — and there is no documented k parameter and no offset/random access.
- **all others = no.** None reports an exact unenumerated size and addresses an arbitrary slice. pydna explicitly raises instead: `raise ValueError(f"Too many assemblies ({possible_assemblies} pre-validation) to assemble")`.

This wording is discriminating, honest about PoolParty's eager paths, and
winnable on real architecture rather than on a word.

---

## Rows unchanged from v2

Blocks B (assay coverage), C (genomics integration), D (adjacent/complementary)
and E (engineering/availability) keep their v2 definitions.

**Caveat: those 20 rows have NOT been through the per-row verification pass.**
They still carry per-tool-agent values with the known rater-variance problem. Do
not put an unverified row in the main-text table without auditing it the same way.
