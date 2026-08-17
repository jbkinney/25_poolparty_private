# Adversarial review — Ledidi extraction

**Reviewer pass:** 2026-08-10. Independent re-verification of every claim against the live repo
(`raw.githubusercontent.com/jmschrei/ledidi@master`, GitHub API), the readthedocs sources, the PyPI
JSON API, the published v2.1.0 tree (commit `7adfcd6453`), and the preprint PDF
(`papers/Schreiber2025_ledidi.pdf`, PyMuPDF text). Nothing was installed or executed.

**Headline:** the extraction is unusually accurate. I re-derived every "no" from primary sources and
found **no understatement of a real capability** — the most dangerous failure mode — with one
partial exception (`combinatorial_motif_place`, where a documented ability is missing from the
evidence text even though the verdict itself survives). I found **one overstatement**
(`design_visualization = yes`), a handful of quote/metadata imprecisions, and one materially wrong
statement in `availability_status` about what `pip install ledidi` actually gives you.

---

## 1. Independently re-verified facts

| Fact | Extraction | My check | OK? |
|---|---|---|---|
| Package layout | `ledidi/{__init__,ledidi,losses,plot,pruning,wrappers}.py` | GitHub tree API: exactly those 6 files + `tests/`, `docs/`, `imgs/`, `slides/` | yes |
| Template restricted to 1 sequence | yes | `ledidi.py:204` `_validate_input(X, "X", shape=(1, -1, -1), ohe=True, allow_N=True)` | yes (quote missing `allow_N=True`) |
| Catalog / repeat / sample axes | yes | `ledidi.py:212–252`; `n_iter = n_samples // batch_size + 1` then `torch.cat` | yes |
| Return type is a bare tensor | yes | `ledidi()` returns `X_bar` (+ optional designers/histories lists) — no record class | yes |
| No `[project.scripts]` | yes | `pyproject.toml` (51 lines) has none; no `SKILL.md` either (unlike sibling tangermeme) | yes |
| LICENSE = Apache-2.0 | yes | `LICENSE` = Apache 2.0 text at both `master` and the v2.1.0 commit; GitHub `spdx_id: Apache-2.0` | yes |
| Last commit 2026-06-23 | yes | `UPDATE plot_edits to render via tangermeme.plot_logo (#19)`, 2026-06-23T21:43:47Z; `pushed_at` 2026-06-23, `updated_at` 2026-07-30, not archived, 110 stars, 0 open issues | yes |
| Only tag/release = v2.0.0 | yes | tags `['v2.0.0']`; releases `[('v2.0.0','2025-01-06')]` | yes |
| PyPI latest 2.1.0, 2025-04-24 | yes | `upload_time` 2025-04-24T13:29:19 | yes |
| Zero hits for barcode/codon/exon/intron/HGVS/VCF/GTF/primer/oligo/restriction | yes | grepped all 6 modules, all 7 `.rst` pages, README, and 6 tutorial notebooks: only hit is the roadmap sentence + one "amino acids" in a docstring | yes |
| Paper: 14 settings, 16.54 edits, 50th→99.27th pct | yes | p.12–13 of PDF | yes |
| Paper: r = 0.893–0.925 at 7 *Drosophila* 8-mer sites | yes | p.25: "ranging between 0.893 and 0.925" | yes |
| Paper: SMYD3, first 1,057 bp masked, GATA2 +0.2…+10.0 step 0.2 | yes | §2.5 verbatim | yes |
| Roadmap quotes (items 3, 4, 5; provenance thread) | yes | README lines 258–270, verbatim | yes |

---

## 2. Findings that change or qualify a value

### 2.1 `design_visualization = yes` → should be **partial** (overstatement)

Three reasons, in decreasing order of force:

1. **Cross-tool inconsistency.** The tangermeme extraction in this same survey scores
   `design_visualization → partial` for a plotting module that is a strict superset of Ledidi's
   (`plot_logo`, `interactive_logo` with annotation boxes/tooltips, `plot_attributions`, `plot_pwm`,
   `plot_categorical_scatter`). Ledidi's `plot_edits` is literally a thin wrapper *on top of*
   `tangermeme.plot.plot_logo` (`plot.py:21,193`). A tool cannot reasonably score higher than the
   library it delegates rendering to. If the referee table shows tangermeme = partial and Ledidi =
   yes, that asymmetry is itself attackable.
2. **It visualizes one design at a time, and only with user-supplied attributions.** `plot_edits`
   takes *attribution* tensors, not sequences: "The calculation of attributions is left to the user"
   (`plot.py` docstring). There is no library-level view, no summary of a design set, no graph view.
3. **The pip-installable release is weaker than master.** In v2.1.0 (commit `7adfcd6453`, the wheel
   PyPI serves today) `ledidi/plot.py` imports `pandas` and `logomaker`, **neither of which is in
   `install_requires` (`torch>=1.9.0`, `matplotlib` only)** — so `from ledidi.plot import plot_edits`
   fails out of the box on a clean install; and `plot_loss` does not exist there at all (added in the
   unreleased 2.2.0). The "yes" therefore rests on unreleased master.

Recommended cell: **partial** — "ships `plot_edits` / `plot_history` / `plot_loss` for a single
design run (edited bases highlighted on user-computed attribution logos, edit positions over
iterations, loss curves); no library-level or design-graph view."

### 2.2 `combinatorial_motif_place = no` — verdict survives, but the evidence is incomplete and exposed

The memo says "there is no API to place a set of motifs across positions, orientations, spacings, or
permutations" — true. But it omits two documented facts a Ledidi author would immediately cite:

- **Ledidi can force a specified motif into a specified position.** Tutorial 2, section *"Adding in
  Motifs"*: "we can also use careful usage of a mask to force edits to happen — even editing in an
  entire motif. … If we set all the characters to `-inf` except for some other character, we could
  force that edit to be included in the design." So motif *placement* is supported; only
  *combinatorial enumeration* is not.
- **Tutorial 3 demonstrates a positional sweep.** "let's take a comprehensive look at what the model
  predictions are if we try moving the 20bp region … downstream between 1 and 800 bp. Every 10th
  basepair … let's try using Ledidi to do the same thing by setting two in-painting regions" — i.e.
  the *user* writes the loop over placements; Ledidi has no operation for it.

Keep "no", but the referee-facing text must read something like: "motifs can be forced in at a
user-chosen position via `-inf` priors, and placements can be swept in a user-written loop, but there
is no operation that enumerates motif sets × positions × orientations × spacings."

### 2.3 `synthesis_constraints = no` — verdict survives, evidence should preempt one objection

Tutorial 2 has a whole section *"Preventing Edits (No Added Cs)"* — set `initial_weights[:,1,:] =
-inf` to forbid ever introducing a C, and the more surgical variant "prevent certain nucleotides from
being edited in, potentially due to issues with prime or base editing", down to blocking a single
specific substitution ("the G in each of the GATAA motif might be blocked from becoming a T"). That
is a real *edit-composition* constraint system. It is still not a synthesis constraint: there is no
GC-content target, no homopolymer/repeat check, no restriction-site avoidance, no Tm, and no
post-hoc synthesizability check on the design (all named as roadmap item 4). Recommend keeping "no"
while naming the `-inf` character/edit-type constraints explicitly, so the row cannot be read as
"Ledidi has no constraints at all".

### 2.4 `per_sequence_provenance = no` — verdict survives; one nuance to be precise about

`history['edits'].append(torch.where(X_hat != X_))` (`ledidi.py:603`) records a *(batch index,
channel, position)* triple per iteration, so the history does carry a batch-member index, not only an
iteration index. It still is not provenance for the returned designs: the returned `best_sequence` is
a clone from whichever iteration last improved the loss, and that iteration index is never recorded,
so you cannot map a returned sequence back to its history rows. Phrase the cell as "no per-sequence
record travels with the output" rather than "history is per-iteration only".

### 2.5 `availability_status` — one materially wrong sentence

The memo says `pip install ledidi` gives v2.1.0 "requires Python >= 3.10 and torch >= 2.0, plus
numpy, matplotlib, tangermeme >= 1.3.0". Those are **master's** `pyproject.toml` requirements. The
published 2.1.0 is a `setup.py` build with `install_requires = ["torch >= 1.9.0", "matplotlib"]`,
`license='MIT License'`, and **no** `requires_python`, **no** numpy, **no** tangermeme. Consequences
worth one clause in the referee response, if availability is quoted at all:

- input validation, `plot_loss`, the tangermeme-based `plot_edits`, and `greedy_pruning(threshold=0)`
  are **not** in any released version — they are unreleased 2.2.0 work on master (`docs/whats_new.rst`);
  `ledidi/__init__.py` on master still reads `__version__ = '2.1.0'`;
- several evidence quotes in the extraction (`_validate_input`, the FAQ, `input_output.rst`,
  `parameters.rst`, Tutorials 0/8, the whole roadmap) come from master, not from the installable
  release. That is the right choice for a capability table — just do not also claim the released
  wheel has them.

### 2.6 `license` — right answer, wrong mechanism

The MIT string is not a "stale classifier": the v2.1.0 `classifiers` list contains **no** license
entry at all; the string comes from `setup.py`'s `license='MIT License'` field in the distributed
artifact. The repo `LICENSE` file was already Apache-2.0 at that same commit, so the package as
shipped carries contradictory license metadata. Conclusion (Apache-2.0, per LICENSE + README badge +
GitHub `spdx_id`) is unchanged; the footnote should be worded accurately in case a Ledidi author
reads it.

### 2.7 Block E value encoding

`interface`, `license`, `maintained` are all recorded with `value: "yes"`. For the other tools in
this survey these rows carry descriptors ("CLI only", "AGPL-3.0-or-later", "last commit 2024-04-22").
"yes" is not an answer to "interface". Suggested cells: `interface = Python API (library) only — no
CLI, GUI, or web service`; `license = Apache-2.0`; `maintained = yes (last commit 2026-06-23)`.

---

## 3. Judgement calls I probed and let stand

- **`library_as_object = partial`** — generous but defensible. One call really can return
  `(n_catalog, n_repeats, n_samples, 4, L)`; against that, the template is hard-restricted to one
  sequence, the return is an unlabeled tensor, and there is no writer. "No" would also be defensible;
  partial is the safer direction for a referee response. Note it is *weaker* than tangermeme's
  partial (tangermeme's ops are batch-first over many templates; Ledidi's is one template).
- **`lazy_evaluation = partial`** — the fitted `Ledidi` module is a genuine on-demand sampler
  (Tutorial 8: "drawing a fresh batch of designs afterward is a single forward pass"; verified in the
  notebook). But successive materializations give *different* stochastic sequences, so it is
  resampling, not deferred evaluation of a fixed spec. "No" is equally defensible; keeping partial
  costs nothing and removes an attack surface.
- **`assay_insilico = partial`** — consistent with the survey's own bar: tangermeme earns "yes"
  because it ships declarative perturbation operations (`marginalize`, `saturation_mutagenesis`,
  `variant_effect.*`); Ledidi ships none of those and reaches in-silico interrogation only through
  repeated optimization. Tutorial 7 (attributions → FIMO → held-out Enformer/ATAC-BPNet round-trip)
  is validation of designs, not a perturbation framework. Partial is right.
- **`assay_mpra = partial`**, **`assay_dms = no`**, **`dag_chaining = no`**, **`genome_coordinates =
  no`** — all re-verified. On `dag_chaining`, note the paper's §2.1 line that masks, in-painting and
  priors "can be seamlessly combined with each other" is constraint composition inside one design,
  not composition of design steps; on `genome_coordinates`, note the tutorials do use `pyfaidx` +
  `tangermeme.utils.one_hot_encode` on hg38 loci — outside Ledidi's API, which is the point.

---

## 4. Capabilities the extraction missed

1. **Knock-down / knock-out design is fully symmetric with knock-in.** The memo frames everything as
   raising activity. Tutorial 2 diminishes GATA2 binding and uses masks to choose *which* of three
   binding sites is destroyed; Tutorial 3 removes an SP-1 site while holding accessibility constant;
   the SMYD3 accessibility/transcription catalogs span −3.0 to +4.0. Any competitor table saying
   Ledidi "designs sequences with higher predicted activity" understates it.
2. **Forced motif insertion via `-inf` priors** (Tutorial 2, "Adding in Motifs") and **motif
   relocation** (Tutorial 3, "Moving Motifs": in-paint the old site + in-paint a landing pad + hold
   the oracle output fixed). This is design-with-placement, and it is the nearest thing Ledidi has to
   PoolParty's motif operations. It belongs in `additional_capabilities`.
3. **Speed as an explicit competitive claim** (paper Fig 2G): Ledidi is "almost an order of magnitude
   faster than greedy substitution and several orders of magnitude faster than a comprehensive motif
   implantation method" — where exhaustive motif implantation (JASPAR2024 at every position) is the
   baseline. If the referee response contrasts enumerative library design with optimization, the
   authors will cite this figure.
4. **Built-in design-validation recipe** — attributions → FIMO motif-hit deltas with decoy motifs →
   round-trip through independently trained models. It is a documented workflow (Tutorial 7) and
   roadmap item 2 promises to make it first-class. Listed as an example in the memo but not as a
   capability.
5. Minor: `target=` for slicing one task out of a multi-task oracle; `verbose` per-iteration edit/loss
   logging; a CPU-runnable test suite with GPU tests behind a pytest marker (`addopts = "-m 'not
   gpu'"`).

---

## 5. Bottom line for the referee response

Ledidi's "no" column is safe: I could not find a plugin, subcommand, notebook, or API that supplies
barcodes, naming, provenance records, coordinates, transcript/codon awareness, HGVS/VCF, primers,
codon optimization, or synthesis constraints. The two sentences most likely to draw an author
objection are (a) any claim that Ledidi cannot place motifs — it can, at a user-chosen position, just
not combinatorially, and (b) `design_visualization = yes`, which is inconsistent with the same
survey's tangermeme row and should be softened to partial.
