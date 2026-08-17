# Contested cells — REVIEW CLOSED 2026-08-17

> **Outcome: one of four cells moved.** `synthesis_constraint_checking`/PoolParty
> `yes` -> `partial`. The other three stand. Verdicts recorded below; the matrix in
> `../tool_survey/MATRIX_FINAL.md` reflects them.
>
> | Cell | Verdict |
> |---|---|
> | `library_object` / tangermeme `partial` | **no change.** The objection is sound — `beam_substitution(n_best=10)` returns what a search converged on, not a declared design space. But adopting that test moves only this one cell (the other `partial`s all have user-determined membership), it moves in our favour, and `partial` already conveys "returns a set, but not a library object". `no` would overstate. |
> | `synthesis_constraint_checking` / PoolParty `yes` | **CORRECTED to `partial`.** The five named filters are real but undocumented — see `../PACKAGE_FINDINGS.md` A1. Row 19 requires named **and documented** checks. |
> | `randomly_sampled_variants` / DNA Chisel `yes` (medium) | **no change.** `pick_random_mutations(n_mutations, sequence)` is documented public API taking a count parameter. Borderline against the row's exclusion of "stochastic search during optimisation", but exposed independently of the solver. |
> | `composable_operations` / Oligopool `yes` (medium) | **no change.** Docs state "Chainable: Output of one module feeds into the next"; `primer()` and `motif()` share the DataFrame contract in either order. Inferred from documentation rather than a worked example, hence medium. |
>
> A scan of all 160 cells found 12 scored below `high` confidence. Two were already
> superseded by batch-7 rescores; of the remaining ten, none changes PoolParty's
> standing — they are competitor cells where the evidence is documentary rather than
> executed, concentrated in the two tools that cannot be installed and run.

## Original contents follow

Held for a single review pass once all 17 rows are scored. The question for each is
narrow: **does this cell earn `partial`, or is `partial` absorbing something that
should be `no`?**

Do not act on these individually. Reviewing them together is the point — a
threshold applied to one cell and not its neighbours is how the earlier
inconsistencies arose.

**This file must not be given to scoring agents.** It records our doubts, and
protocol rule 6 in `table2_candidate_17row.md` forbids handing expected outcomes to
scorers.

---

## 1. tangermeme — `library_object` = partial

**Raised by:** the user, 2026-08-15.

**Objection.** "Library" may not mean the same thing for tangermeme as for the other
tools. Its `beam_substitution(n_best=10)` returns the ten sequences a beam search
converged on; PoolParty's `Pool` is a declared design space representing every
variant the specification produces. One is *search output*, the other is *the
library you intend to build*. Membership in the first is determined by how the
optimiser ran; in the second, by the user's specification.

**Why it currently scores partial.** Row 1 tests whether the tool returns the whole
collection through its own API, and in what container type. tangermeme does return
the set — `n_best` and the documented return shape `(n_best, alphabet, length)` are
tool-provided, not caller-assembled — but as a bare tensor with no member identity,
provenance or library operations. That is the definition's `partial` case as
written.

**What is missing from row 1.** It never says what makes a collection a *library*
rather than a set of sequences. A candidate test: *is membership determined by the
user's design specification, or by what an algorithm found?*

**Caution recorded at the time.** Adopting that test moves tangermeme to `no`, which
**improves PoolParty's relative standing**. Every earlier definitional fix cut
against us or was neutral; this one does not. Scrutinise accordingly.

**Also affected if adopted.** DNA Chisel's `MutationSpace` is a declared space but a
solver search space rather than a deliverable; Oligopool's DataFrame is the user's
variant table plus tool-added infrastructure. Both would need re-examining under the
same test, so this cannot be applied to tangermeme alone.

---

## 2. Cells needing targeted rescoring, not review

Distinct from the above — these are unresolved because a definition changed under
them, not because the definition is doubted.

| Row | Tools | Why |
|---|---|---|
| `saturation_mutagenesis` | Oligopool Calculator, DNA Chisel | `partial` on grounds unrelated to the struck "user's responsibility" branch |
| `mixed_variant_types` | Mutation Maker | cites two grounds — "same-kind components" survives, "separate runs" was struck |
| `mixed_variant_types` | Oligopool Calculator, DNA Chisel | `partial` on unstated grounds |
| `mixed_variant_types` | MPRAnator, MPRA Design Tools | scored `yes`; confirm the pooling is tool-side |

Scheduled as batch 7.

---

## 3. Open TODO — codon-aware indels (does NOT block the table)

**Status: deferred by decision, 2026-08-17. No row change required now.**

**Finding.** PoolParty v2 has no codon-aware insertion or deletion. `deletion_scan`
and `insertion_scan` take no `frame` parameter and contain zero references to
codons or ORFs; `grep -rn "inframe|in_frame|in-frame"` over the package returns
nothing. An in-frame deletion requires the user to pass `deletion_length=3` with
codon-aligned `positions` they compute themselves — user reconstruction under the
global rule.

By contrast `mutagenize_orf` *does* take an explicit `frame` parameter (±1/2/3), so
the codon machinery exists in the package; it is simply not wired to the indel
operations. VaLiAnT has `inframe` as a first-class mutator, and in-frame deletions
are a standard DMS design.

**This is a regression, not an architectural limit.** PoolParty v1 shipped
`DeletionScanORFPool` and `InsertionScanORFPool`: *"Scan deletions across an ORF at
the codon level … Operations work at the codon level to maintain reading frame
integrity"*, with `deletion_size` and `start`/`end`/`step_size` all in **codon
units**, plus a position-based interface with per-position weights. That is finer
control than VaLiAnT's `inframe`.

Corrected v1 -> v2 mapping (an earlier note overstated the loss as three classes):

| v1 class | v2 status |
|---|---|
| `KMutationORFPool` | **absorbed** into `mutagenize_orf(num_mutations=k)` — verified: k=1/2/3 gives 171 / 12,996 / 576,156 states |
| `RandomMutationORFPool` | **absorbed** into `mutagenize_orf(mutation_rate=…, mode='random')` — verified |
| `DeletionScanORFPool` | **lost** |
| `InsertionScanORFPool` | **lost** |

So two of four were consolidated into a more general operation — a refactor, arguably
an improvement. Only the two indel classes were genuinely lost.

**Effect on the table: none as it stands.** Row 8's definition does not distinguish
nucleotide-level from codon-aware indels, so PoolParty's `yes` is correct under the
wording actually used. Row 7 is unaffected — it is substitutions only.

**The three options, and their triggers:**

1. **Restore in v2** before resubmission. Bounded work — porting two v1 classes onto
   the v2 operation model, not rebuilding a subsystem. Row 8 then stays `yes`
   honestly, and nothing in the table changes.
2. **Sharpen row 8** to distinguish codon-aware indels. Costs a rescore of all eight
   tools on that row; PoolParty likely drops to `partial`, VaLiAnT holds `yes`.
3. **Leave row 8 undifferentiated** and say nothing. Defensible on the wording, but
   a referee who knows VaLiAnT's `inframe` may ask.

**Do NOT** take option 3 *and* describe the gap in the Discussion as a structural
consequence of the DAG architecture. v1 proves it is not, and the v1 code is public
history. That framing is available for the genomic-coordinate gap; it is not
available here.

## 4. Rows flagged during scoring

| Row | Finding | Status |
|---|---|---|
| `recombination` | Scored PoolParty `yes`, all seven others `no`. Fails the discrimination ground of the substitution protocol. | **Decision needed** — drop to 16 rows, substitute, or keep and defend |
