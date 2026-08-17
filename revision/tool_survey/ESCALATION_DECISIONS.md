# Escalation rulings — repair pass `wrjc3a1sw`

The repair pass returned **338 applied · 45 rejected · 62 escalated**. PoolParty's four were
decided separately and applied. This records the rulings on the remaining 58.

**Headline: only three of the 58 could move a cell that ships, and none of them did.**

---

## The three cell-affecting escalations

Table 2 ships eight rows. Only these three escalations touched one.

### 1. pydna `library_as_object` — **record annotated, shipping value unchanged**

The record's justification rested on "no collection type appears among the 820 documented objects",
which the citation audit withdrew (the index does contain `RecombinaseCollection`).

Verified position: under the **original** three-criterion rubric pydna is **partial** — criteria (i)
and (iii) pass via `Assembly` and `assembly2.get_possible_assembly_number`; only (ii) fails. Under
the **adopted reworded** rubric in `ROWS_v3.md`, which deleted criteria (ii) and (iii) because they
duplicated other rows, pydna is **yes**.

So the record's `no` was two generations stale. **Ruling:** annotate the entry as superseded and
point to `verified/library_first_class_object.md`; do not re-litigate. Table 2 already shows ● for
pydna via footnote a. **Applied.**

### 2. DNA Chisel `lazy_evaluation` — **no change; already correct**

The escalation asked whether `DnaNotationPattern.all_variants()` should be added. It is already in
the record at line 58, and the record already carries the `no → partial` correction with the note
"Extraction understated the tool". Under the v3 replacement row
(`unmaterialized_library_addressing`) DNA Chisel remains **partial** because `space_size` is a float
capped at e^100 and there is no offset addressing — a degenerate-pattern expander fixes neither.
**No edit needed.**

### 3. tangermeme `lazy_evaluation` — **no change; already correct**

Ruling on the definitional question: `ROWS_v3` defines the row as "an arbitrary k-member slice of a
library larger than memory". **Genome-scanning streaming does not qualify.** The record already
states exactly this — "Streaming exists only on the Cartesian-product axis, never as a library
abstraction" — and rests the `partial` on `product.py`, not `match.py`. **No edit needed.**

---

## Policies for the remaining 55

### Policy A — balance and emphasis · 9 items · **decline all**

One per tool, always the same question: should the record be rebalanced? These are working
documents that justify a table; they are not deliverables and no referee sees them. Record-level
emphasis is out of scope. The one place balance genuinely mattered — PoolParty's own record — was
decided separately and the duplication there has been trimmed.

*Covers:* valiant §C, mpradesign, mutationmaker, codongenie, dnachisel, tangermeme, biopython,
pydna, seqpro.

### Policy B — materiality of omissions · 10 items, ~60 sub-items · **rule, not case-by-case**

**Add an omission only if it bears on one of Table 2's eight shipping rows or on a Table 1 prose
cell. Otherwise decline.**

tangermeme alone proposed 15 and pydna 14; nearly all are API details with no home in the row set.
The repair agents already applied the subset that had an obvious existing home, which is the right
outcome under this rule.

### Policy C — sourcing and provenance · 10 items · **primary sources only**

Admissible: the tool's repository, its documentation site, its package index page, its paper.

- Third-party aggregators (bio.tools) — **not admitted**.
- Cross-citation between sibling records — **not admitted**; each record must stand alone.
- Unverifiable anecdotes — **dropped**, not marked.
- Live reproductions performed during a pass — need not be retained as artefacts.

### Policy D — version drift · 4 items · **header version governs**

A record describes the version named in its header. Where the current release differs materially,
add one parenthetical; do not renumber. Same principle applied to PoolParty's stale checkout, where
the assessed tree was deliberately kept.

### Remaining ~22 tool-specific items — **leave**

None touches a shipping cell. They are quality improvements to documents whose purpose is to
justify a table that is now traceable to verified evidence. Further polishing does not reach the
paper.

---

## Outstanding

**Reconciliation has not been done.** 338 corrections were applied to the evidence records, and the
ratings in `verified/*.md` and `MATRIX_verified.md` have not been re-checked against the corrected
evidence. The three cell-affecting escalations were checked individually and none moved, which is
reassuring but not the same as a systematic pass.
