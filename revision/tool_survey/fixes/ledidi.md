# Ledidi — repair pass change log

**Record repaired:** `final/ledidi.md` (edited in place, 2026-08-14)
**Audits processed:** `citation_audit/ledidi.md` (12 findings), `factcheck/ledidi.md` (A: 8, B: 8, C: 1)
**Outcome counts:** 29 findings → 24 applied · 2 rejected · 3 escalated · 1 no-action (§C balance).
C7 appears twice: applied at one of its two locations, escalated at the other. Two findings are
duplicated across the audits (C5≡A8, C12≡A6) and were fixed by a single edit each, so the 24 applied
findings correspond to 22 edit sites.

**No capability value was changed.** All 24 values (`library_as_object` … `maintained`) are
byte-identical to the pre-repair record; only evidence and prose were touched.

## Primary sources I re-derived everything from

| Source | How accessed |
|---|---|
| Repo master `adbca708d45340fb7f375e4d0e2438d3cfa7c852` | `git clone`; identical to the commit both audits cite |
| Release tree `7adfcd6453a6851e2145bec328d5063bd19f4c4f` (v2.1.0) | `git worktree` of that commit (`setup.py` says `version='2.1.0'`) |
| Paper | `papers/Schreiber2025_ledidi.pdf`, 40 pp., PyMuPDF text dump (title verified: "Programmatic design and editing of cis-regulatory elements") |
| Live docs | `curl` HTTP status for every page cited by the record |
| Zenodo 14604495 | `https://zenodo.org/api/records/14604495` file inventory |

Nothing was installed; no code was executed beyond `git`, `curl` and PyMuPDF text extraction.

---

## Citation-integrity audit

### C1 — Tutorial 2 quotation silently fixes the source's typo · **APPLIED**
- **Verified:** `docs/tutorials/Tutorial_2_-_Constraints_and_Priors.ipynb`, cell 26 reads "we could
  force that edit **ot** be included in the design."
- **Correction:** quotation restored to the source spelling with `[sic]` (`combinatorial_motif_place`
  evidence). Meaning unchanged.

### C2 — README modality quotation drops a clause · **APPLIED**
- **Verified:** `README.md:27` ends "…any model with a sequence of categorical inputs **such as small
  molecules**)."
- **Correction:** the deleted words restored inside the quotation (`assay_dms` evidence).

### C3 — Affinity-catalog quotation attributed to §3.3 / Fig 4E · **APPLIED**
- **Verified:** the sentence occurs exactly once, at the end of the **Introduction** on PDF p. 4
  (paper's own line 103); §2 Methods begins immediately after it on the same page. The Abstract
  (pp. 1–2) contains neither quotation.
- **Correction:** locator changed to "paper Introduction, p. 4". I also corrected the adjacent
  locator in the same parenthetical — "abstract:" — because I verified the second quotation
  ("provide a resource for understanding the full range of potential sequence edits…") occurs
  exactly once, in the *same* Introduction paragraph; it now reads "same paragraph:". The **Source:**
  line's "abstract" was changed to "Introduction p. 4"; its "§3.3, Fig 4E" was left (§3.3 is titled
  "Affinity catalogs reveal how Ledidi uses cis-regulatory logic", so it does support the claim the
  sentence makes — just not the quotation).
- **Note:** the audit called p. 4 "the abstract"; that part of the audit is wrong, but its core claim
  (not §3.3/Fig 4E) is right.

### C4 — ReadTheDocs source inventory · **APPLIED**
- **Verified by HTTP:** `en/latest/` root, `index.html`, `whats_new.html`, the 5 `api/*.html` pages
  and Tutorials 1–6 return 200; `getting_started.html`, `input_output.html`, `parameters.html`,
  `faq.html` and Tutorials 0, 7, 8 return **404**. The rendered page's version banner reads `v2.0.0`.
- **Correction:** the Sources-table "Docs" row now cites `docs/` at repo master as the location of
  those files, with the stale deployed build named as such.

### C5 / A8 — "Nine tutorial notebooks … all rendered on readthedocs" · **APPLIED**
- **Verified:** nine notebooks exist at master; the deployed build indexes and serves only 1–6 (see C4).
- **Correction:** "all rendered on readthedocs" replaced with the stale-build fact.

### C6 — Zenodo provenance for tutorial oracles · **APPLIED**
- **Verified:** record 14604495 is titled "BPNet models used as Ledidi examples" and contains 24 files,
  all BPNet `.torch`/`.fit.json` (GATA2, E2F3, ATAC, CTCF, MAX, MYC, SOX6, RAD21, OTX1, MECOM, JUNB,
  ELF4). No Enformer, Malinois, ChromBPNet or Beluga file. Tutorial 4/7 load Enformer via
  `from_pretrained('EleutherAI/enformer-official-rough')` (`enformer_pytorch`); Malinois, ChromBPNet
  and Beluga are local checkpoint paths with no stated source. Tutorial 1 says the link provides "the
  GATA2 model … as well as the other **BPNet** models used in the paper".
- **Correction:** the two overreaching sentences (vignette preamble; "Repo/service health") now
  attribute only the BPNet oracles to Zenodo. The Sources-table row (line 43) already said "BPNet
  oracles" and was left untouched — the audit's blanket claim about that row is not correct.

### C7 — `.rst` page count · **APPLIED (one location) / ESCALATED (the other)**
- **Verified:** master has **11** `.rst` files (6 doc pages + 5 API stubs); the v2.1.0 tree has 7
  (+1 template). Re-ran the grep myself: zero hits for `barcode|oligo|adapter` across all 6 modules,
  all 11 `.rst` files, the README and all 9 notebooks — the finding's conclusion is unaffected.
- **Correction:** `barcode_generation` evidence now says "all 11 `.rst` pages".
- **Escalated:** the identical number in the confidence section describes *the reviewer's* search
  universe in an earlier pass — see E1.

### C8 — "Every readout in the repo and the paper is either an oracle prediction or an edit count/position" · **APPLIED**
- **Verified:** the repo and paper also report attribution scores (Tutorial 7, `plot_edits`), FIMO
  motif-hit counts (Tutorial 7 cells 19–24), correlations with experimental data (paper "ranging
  between 0.893 and 0.925") and runtime benchmarks (Fig 2G). The exhaustive form is false.
- **Correction:** replaced with a list of the readouts that do occur plus the narrow negative claim
  ("never include a variant-consequence class"), which is what the `no` rests on. Value untouched.

### C9 — E2F6 vs E2F3 · **REJECTED**
- **Why wrong:** the record's cited source is the paper's Methods/Results, and the main text and Fig 2
  say **E2F6** throughout (paper lines 795, 850, 935, 971, 1198, 1213–1219); the tutorials also load
  `E2F6.torch`. Only Supplementary Table 1 (PDF p. 39) says E2F3, and Zenodo hosts `E2F3.torch`. The
  discrepancy is inside the paper/model zoo, not in the record; "correcting" the record to E2F3 would
  contradict its own cited source.
- **Worth knowing:** `E2F6.torch`, which Tutorials 4 and 6 load, is **not** in Zenodo 14604495.

### C10 — "17 curated motifs" · **APPLIED**
- **Verified:** the paper says "a set of 17 motifs thought to be relevant"; the string "curated" does
  not occur in the paper text.
- **Correction:** parenthetical now quotes the paper's own wording.

### C11 — FAQ quotation truncated after "one" · **APPLIED**
- **Verified:** `docs/faq.rst:57` continues "…Ledidi commits to one **and all sampled sequences
  reflect that choice**."
- **Correction:** quotation completed (the omitted clause is what makes the correlation claim).

### C12 / A6 — `greedy_pruning(threshold=0)` listed as unreleased master work · **APPLIED**
- **Verified in the v2.1.0 source:** `pruning.py` has no validation at all and its only stopping test
  is `if best_score < threshold`, so `threshold=0` is accepted and simply reverts nothing. The master
  changelog entry reads "Fixed … to accept a `threshold` of 0 … **again**", and commit `fbb277c`
  describes the intervening rejection as "a back-incompatible regression from 2.1.0".
- **Correction:** removed from both release-lag lists (Stated limitations; Availability). The slot was
  filled by `random_state`'s private generator, which *is* master-only (see A7).

---

## Fact-check §A (factually wrong)

### A1 — "No codon-aware or amino-acid-level mutagenesis" is too broad · **REJECTED**
- **Why wrong:** the clause is one item in a list of *absent operations* ("no reading-frame logic, no
  exhaustive variant enumeration, no coding-sequence handling…"), and the record's very next sentence
  states the README modality claim verbatim and calls it "a modality claim about the *method*". Item
  12 of the capabilities list also credits Ledidi as "modality-agnostic — DNA/RNA/protein". Nothing in
  the record asserts that Ledidi cannot edit amino-acid sequences.
- **Evidence:** v2.1.0 `ledidi.py:182–208` says sequences may be "comprised of nucleotides or amino
  acids" — a statement about the categorical alphabet, not about a mutagenesis operation, of which
  Ledidi ships none. Rewriting here would also alter what `assay_dms = no` rests on.

### A2 — "It renders one design run at a time" is false of `plot_edits` · **APPLIED**
- **Verified:** `plot_edits(X_orig, X_attrs, …)` takes `X_attrs` as an `(n, 4, length)` tensor *or a
  list that is concatenated*, and draws one track per row; nothing ties those rows to a single run
  (v2.1.0 `plot.py:75–79`, master `plot.py:173–177`). What is actually fixed is the single reference
  `X_orig` of shape `(1, 4, len)` and a common length. (`plot_loss`/`plot_history` do take a single
  history dict.)
- **Correction:** the clause now reads "renders edits against **one reference sequence at a time**".
- **Tension noted (two items):**
  1. This is leg 2 of three supporting `design_visualization = partial`. The scope argument it exists
     to make ("no library-level view, no design-set summary") is untouched and the value is unchanged,
     but the authors may want to re-read that leg.
  2. The confidence section still summarizes the same cell as "the plotting … is one design at a
     time". The audit flagged only the evidence cell, and the instruction is to leave nearby unfixed
     sentences alone, so I did — but that summary now sits slightly out of step with the corrected
     evidence and should be reconciled in the same pass that settles item 1.

### A3 — "Documented design-validation recipe as a capability, not just an example" · **ESCALATED** (E2)

### A4 — "No biological or synthesis constraints yet" · **APPLIED**
- **Verified:** the record's own `synthesis_constraints` row and capability 6 document a real
  constraint system over edits (`input_mask`; `initial_weights = -inf` blocking a character, a
  substitution type, or forcing a motif) — implemented at master `ledidi.py:551–555` and documented in
  `docs/input_output.rst` "Masks and priors". Roadmap item 4 covers only GC content, restriction
  sites, codon usage and off-target effects.
- **Correction:** the bullet header is now "No **global sequence-property** or synthesis constraints
  yet", with a pointer to `synthesis_constraints`. This aligns the summary bullet with the row it
  summarizes; the `no` value is unchanged.

### A5 — "Requires a differentiable pre-trained model" · **APPLIED**
- **Verified:** `README.md:44–47` — "a complete, runnable example that uses a tiny parameter-free
  oracle, so there is nothing to download"; the oracle is a `conv1d` over a fixed AP-1 weight tensor
  and the call passes `device="cpu"`.
- **Correction:** "pre-trained" dropped; the quickstart named as the counter-example.

### A6 — see C12. **APPLIED**

### A7 — `random_state` presented without release qualification · **APPLIED**
- **Verified:** v2.1.0 `def ledidi(model, X, y_bar, n_repeats, n_samples, return_designer,
  return_history, device, **kwargs)` and `Ledidi.__init__(…, return_history, verbose)` accept no
  `random_state` (though the 2.1.0 class docstring documents one — stale); its `forward` calls
  `torch.nn.functional.gumbel_softmax`, i.e. the global RNG. Master adds `random_state=` to both
  signatures and a private `torch.Generator` (`ledidi.py:497–503`).
- **Correction:** capability 11 now carries "(master only; the released v2.1.0 signatures do not accept
  `random_state` despite a stale docstring…)", and `random_state` was added to both release-lag lists.
- **Left alone (local roughness):** the FAQ-derived limitation "Independent designs require
  `n_repeats` or a new `random_state`" quotes master's FAQ and is correct for master.

### A8 — see C5. **APPLIED**

---

## Fact-check §B (incomplete) — each added as one clause/sentence inside an existing entry

| # | Fact added | Where | Verified against |
|---|---|---|---|
| B1 | `input_loss` is user-replaceable, so which changes count as costly is definable | capability 5 | `Ledidi.__init__(input_loss=…)` and its use at master `ledidi.py:583` (`self.input_loss(X_hat[…], X_[…])`) |
| B2 | oracle-adapter wrappers for differing input widths, and the "hide the edits from the models with the smaller receptive field" hazard unless unshared flanks are masked | Tutorial 4 vignette | Tutorial 4, cell 59 (quoted verbatim) |
| B3 | finite priors "nudge, but do not force, the design" and can be optimized away; only `-inf` is hard | capability 6 | `docs/input_output.rst` "Masks and priors" (quoted verbatim) |
| B4 | in-painting is preferential, not exclusive: outside positions stay editable unless masked; `N` columns are not preserved | capability 7 | `inpainting_mask = X[0].sum(dim=0) == 1` gates only the *input loss* (master `ledidi.py:557, 583`); `docs/input_output.rst:18` |
| B5 | `MinGap` ignores `y_bar` and maximizes the weakest-on-target/strongest-off-target gap — relative separation, not absolute activity | capability 4 | `losses.py:67–94` and its docstring; `docs/input_output.rst` ("ignore `y_bar` entirely") |
| B6 | `greedy_pruning` handles one designed sequence per call, scored by summed absolute output change vs the fully edited design | capability 5 | `pruning.py:18` ("Only one sequence is pruned at a time"), `pruning.py:66, 77–78` |
| B7 | v2.1.0 `plot_edits` hard-codes the DNA alphabet `A, C, G, T` | Availability bullet | v2.1.0 `plot.py:83, 100` |
| B8 | paper's language-model composition suggestions | — | **ESCALATED (E3)** |

## Fact-check §C (balance)

No action. The audit's own verdict is that the record is balanced; balance is an authors' call in any
case.

---

## Escalations (record left untouched at these points)

**E1 — the reviewer's `.rst` search universe (confidence section).**
The sentence "The reviewer independently re-derived every `no` from primary sources (6 modules, 7
`.rst` pages, README, 6 of 9 notebooks, …)" reports what a previous pass read. Master has 11 `.rst`
files; the v2.1.0 tree has 7. I cannot verify which tree that reviewer actually walked, and rewriting
the sentence would put an unverifiable claim about someone else's process into the record.
*Options:* (a) change 7 → 11 if the reviewer worked from master; (b) write "the release-era 7 `.rst`
pages" if they worked from the v2.1.0 tree; (c) drop the count. Note that I re-ran the decisive grep
over all 11 master `.rst` files myself and it is clean, so no `no` is at risk either way.

**E2 — "Documented design-validation recipe **as a capability, not just an example**" (item 10).**
The fact-check says this misattributes a tutorial workflow to Ledidi: Tutorial 7's attribution, FIMO
and held-out-oracle steps are user/third-party code (`tangermeme`, `memelite.fimo`, Enformer), Ledidi
exposes no validation API, and the README roadmap promises first-class evaluation as future work — all
of which I confirmed. But the record's phrasing was added deliberately "*(added per review)*", it
already ends with "roadmap item 2 promises to make it first-class", and the `assay_insilico` row
already calls Tutorial 7 "*design validation*, not a perturbation framework". Whether the heading
overclaims is an emphasis judgement between two prior passes.
*Options:* (a) leave as is; (b) retitle to "Documented design-validation recipe (documentation, not a
shipped API)"; (c) drop "as a capability, not just an example".

**E3 — the paper's language-model compositions (fact-check B8).**
Confirmed in the paper's Discussion: adding a differentiable genomic language model as a realism term
in a multi-model objective, and using a language model upstream to generate starting templates that
Ledidi then refines. Neither is shipped or demonstrated. There is no existing entry this belongs to —
it is not a shipped capability, not a limitation, and not a vignette — and the instruction is not to
add sections.
*Options:* (a) omit (status quo); (b) one sentence in "Notable capabilities not covered by the row
list" marked as paper-suggested only; (c) one sentence in "Stated limitations" as the authors' stated
route beyond a hand-chosen natural template.

---

## Edits made, in file order

| Record location | Change | Finding |
|---|---|---|
| Sources table, "Docs" row | re-pointed to `docs/` at master; stale live build named | C4 |
| `combinatorial_motif_place` evidence | `ot [sic]` restored in Tutorial 2 quote | C1 |
| `barcode_generation` evidence | "7 `.rst` pages" → "11 `.rst` pages" | C7 |
| `design_visualization` reason 2 | "one design run at a time" → "one reference sequence at a time" | A2 |
| `assay_dms` evidence | README quote completed ("such as small molecules") | C2 |
| `assay_insilico` evidence + Source | quote locators corrected to Introduction p. 4 | C3 |
| `consequence_annotation` evidence | exhaustive readout claim narrowed | C8 |
| Vignettes preamble | readthedocs staleness; oracle provenance split out | C5, C6, A8 |
| Tutorial 4 vignette | adapters + receptive-field hazard clause | B2 |
| Capability 4 | `MinGap` semantics sentence | B5 |
| Capability 5 | `input_loss` replaceable; pruning scope/metric | B1, B6 |
| Capability 6 | finite priors nudge, not force | B3 |
| Capability 7 | in-painting preferential, not exclusive | B4 |
| Capability 9 | "17 curated motifs" → paper's wording | C10 |
| Capability 11 | `random_state` marked master-only | A7 |
| Limitations: constraints bullet | narrowed to global sequence-property/synthesis | A4 |
| Limitations: correlated-batch bullet | FAQ quote completed | C11 |
| Limitations: oracle bullet | "pre-trained" dropped | A5 |
| Limitations: release-lag bullet | `threshold=0` → `random_state` | C12, A6, A7 |
| Availability: v2.1.0 plotting bullet | DNA-alphabet hard-coding added | B7 |
| Availability: unreleased-work bullet | `threshold=0` → `random_state` | C12, A6, A7 |
| Availability: repo/service health | docs staleness; Zenodo scoped to BPNet | C4, C6 |

Whitespace-only re-wraps were made on the touched paragraphs to preserve the file's ~100-column
convention; no other text was altered.

## Pass 2 — policy application

**Open-item outcome:** 3 prior escalations → **2 applied/resolved · 1 rejected/declined · 0 still
escalated**. Two supporting policy corrections were also applied at neighbouring/value-basis text.
There are four edit sites in `final/ledidi.md`; **no capability value changed**, no entry or section
was added, and Table 1/Table 2 were not edited.

### E1 — unverifiable reviewer search-universe anecdote · **APPLIED / RESOLVED**

- **Outcome:** applied provenance ruling C: drop the anecdote rather than guess whether the earlier
  reviewer searched the 7 release-era or 11 master `.rst` files.
- **Edit:** deleted the confidence-section sentence beginning *"The reviewer independently
  re-derived every `no`…"*. The surrounding primary-source confidence statement remains.
- **Verification:** repository master `adbca708d45340fb7f375e4d0e2438d3cfa7c852` contains 6 package
  modules, 11 `.rst` files and 9 notebooks; the v2.1.0 tree is separately identifiable at
  `7adfcd6453a6851e2145bec328d5063bd19f4c4f`. Those trees verify the software evidence, but neither
  can verify what a previous reviewer personally read. Deletion is therefore the only supported
  correction; no capability evidence or value depends on that anecdote.

### E2 / fact-check A3 — Tutorial 7 labelled as a capability · **APPLIED / RESOLVED**

- **Outcome:** corrected a factual category error; this is documentation, not a shipped Ledidi
  validation capability. It was not promoted into a Table 1 cell.
- **Edit:** item 10 now reads *"Documented design-validation example, not a shipped Ledidi
  capability"* and says the checks are notebook code using external packages/models; Ledidi exposes
  no validation API, while README roadmap item 2 describes first-class evaluation as future work.
- **Verification:** the primary Tutorial 7 notebook imports and calls `tangermeme` attribution/
  prediction utilities, `memelite.fimo`, and separately loaded Enformer/BPNet models. The complete
  public definitions in the six Ledidi modules are `ledidi`, `Ledidi`, `MinGap`, `plot_loss`,
  `plot_history`, `plot_edits`, `greedy_pruning`, and `DesignWrapper`; none is a validation API.
  README roadmap item 2 says the authors *"want first-class evaluation built in"*. These first-party
  sources support the corrected label.

### E3 / fact-check B8 — paper-suggested language-model compositions · **REJECTED / DECLINED**

- **Outcome:** verified but not added. Under ruling B, an omission ships only if it bears on a Table
  1 Purpose, Key features, Output, or Availability cell. These two Discussion suggestions are neither
  implemented nor demonstrated and do not change the grouped adjacent-tools row's description of
  Ledidi as model-guided editing.
- **Edit:** none.
- **Verification:** paper Discussion, lines 681–690 (PDF pp. 29–30), proposes (1) a genomic language
  model as a “genome-like” term in a multi-term objective and (2) a language model generating initial
  sequences for later Ledidi refinement. The paper presents both as possible future compositions;
  the repository exposes no bundled language-model integration.

### Balance/emphasis ruling A · **DECLINED / NO EDIT**

- **Outcome:** no balance rewrite. The fact-check's own §C conclusion is that the record is balanced,
  and the ruling expressly declines emphasis-only edits.
- **Verification:** no factual claim was at issue; all specific factual A findings were already
  disposed of in Pass 1, with E2 above the only one left open.

### Value basis, provenance, and neighbour tension · **APPLIED**

- **Outcome:** `design_visualization = partial` remains unchanged; Ledidi's own repo/release evidence
  supports it without relying on a sibling tool's survey score.
- **Edits:** replaced the tangermeme cross-tool comparison with the Ledidi module's own description of
  its functions as *"thin matplotlib helpers"*, user-supplied attribution input, and out-of-scope
  attribution computation/broader plotting. Removed the sibling-row citation. Reconciled the nearby
  confidence summary from *"one design at a time"* to *"one reference sequence at a time"* and
  removed sibling delegation as a scoring premise.
- **Verification:** master `plot.py:12–14, 118–193` exposes three plotting helpers, accepts one
  `X_orig` reference, requires caller-computed `X_attrs`, and has no library/design-set view. Release
  v2.1.0 defines only `plot_history`/`plot_edits`; its `plot.py` imports undeclared `pandas` and
  `logomaker`, and `setup.py` declares only PyTorch and matplotlib. This directly supports `partial`.
  The neighbour wording now matches the already-corrected capability evidence.

### Version drift ruling D · **VERIFIED / NO EDIT**

- **Outcome:** the one-parenthetical drift allowance was not used because there is no material drift.
- **Verification:** live PyPI JSON on 2026-08-14 still reports **2.1.0** as latest. The live GitHub API
  still reports default-branch head `adbca708d45340fb7f375e4d0e2438d3cfa7c852` (2026-06-23), whose
  `ledidi/__init__.py` remains `2.1.0` and whose docs still label 2.2.0 unreleased. This matches the
  record's release/master split.

### Locked-row / substitution report · **NO ACTION**

- **Confirmed uniform rows:** none can be established from a Ledidi-only pass. The locked reference
  already flags `Codon / amino-acid substitutions` as a possible near-uniform row and
  `Recombination and chimeras` as a possible one-tool row; both remain unscored across the eight
  shipping columns, so neither is changed here.
- **Candidate considered:** `Target-directed minimal editing` could replace `Recombination and
  chimeras`; paper Equation 1 combines distance to an oracle target with an input loss penalizing
  edits, and the repo adds `greedy_pruning`. **Rejected as a substitution candidate in this pass:**
  Ledidi has no Table 2 column, `Model-guided variants` already captures the relevant axis, and there
  is no cross-tool evidence that the narrower row would discriminate among the eight included tools.
  No locked row was edited and no substitution was logged in the table file.
