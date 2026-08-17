# Adversarial review — tangermeme extraction

**Reviewed:** 2026-08-10. **Reviewer method:** independent re-download of the full `main` tree
(GitHub trees API), all 25 `.py` files under `tangermeme/` **including `design/` and `_skills/`**,
all 14 bundled skill reference docs, `README.md`, `pyproject.toml`, `docs/whats_new.rst`,
the live PyPI JSON, the GitHub repo metadata, and the paper PDF (PyMuPDF, all 13 pp.).
Nothing was installed. I did **not** re-use the extractor's quotes; every quote below was
re-read from source.

**Bottom line:** the extraction is unusually solid. All 25 capability values survive
falsification. No verdict needs to change. But there are **five evidence gaps** that a
referee — plausibly Schreiber himself — could use to make the extraction look careless,
and **one small factual error** in a cited source. Fix those and the row is defensible.

---

## 1. Facts I independently confirmed (no change needed)

| Claim | Confirmed how |
|---|---|
| Only sequence writer in the package is `io.one_hot_to_fasta`; header falls back to numeric index | `grep -n "\.write(\|to_csv" tangermeme/**/*.py` → hits only in `io.py` 580–587; body is `outfile.write("> {}\n".format(i))` |
| `substitute` / `multisubstitute` / `space` take a **scalar** `start` (position is not a swept axis) | `ersatz.py` L97–103, L198–205; `space.py` L34–44 — `start: int \| None = None` in all three |
| `space(spacing)` is 2-D `(n_spacings, n_motifs-1)` and is the only real enumerated axis | `space.py` docstring verbatim + `_validate_input(..., shape=(-1, len(motifs)-1))` |
| `apply_pairwise` / `apply_product` cross `X` with **extra model inputs**, not motifs/positions | `product.py` L30–60 docstring: "the first being sequence, and the second being anything else"; the worked example is DragoNNFruit cell-state × read-depth |
| Product streaming is real | `product.md` L50: "Batches are built iteratively, so the full product is never materialized in memory" |
| Zero hits, whole package incl. `design/` and `_skills/`, for: barcode, primer, codon, hgvs, gtf, gff, transcript, exon, intron, synthesis, homopolymer, restriction, amino, protein, orf, reading frame, oligo | my own recursive grep (see §3.5 — the extractor's grep was narrower than it claims) |
| `read_vcf` input-only; no VCF/VEP writer | `io.py` public readers are exactly `extract_loci`, `one_hot_to_fasta`, `read_meme`, `read_vcf` |
| `annotate.*` = TOMTOM motif matching of seqlets, not molecular consequence | `annotate.py` functions: `annotate_seqlets`, `count_annotations`, `pairwise_annotations`, `pairwise_annotations_spacing` |
| MIT; `[project.scripts]` = only `tangermeme-install-skills`; requires-python >=3.10; torch>=2.0 | `pyproject.toml` verbatim; PyPI `license_expression: MIT` |
| v1.4.0 uploaded **2026-06-25**; repo `pushed_at` **2026-07-19T18:34:08Z**; `archived: false`; 308 stars / 32 forks; "Development Status :: 5 - Production/Stable" | PyPI JSON + GitHub `/repos` JSON, both fetched fresh |
| Release cadence 1.0.0 (2025-08-27) → 1.4.0 (2026-06-25), 8 releases in 12 months | PyPI `releases` upload_time, exact match to the extraction |
| Paper quotes ("everything-but-the-model", "operations can be stacked", Fig 1A grid incl. Model Training / Data Preprocessing / Model Zoo as *not* tangermeme, Code Availability MIT) | PDF pp. 1–3, 9 |

---

## 2. The five things a referee could hit — evidence gaps, not wrong values

These do **not** change any value. They change what the evidence must say so the value
cannot be attacked.

### 2.1 `ersatz.randomize` / `shuffle` / `dinucleotide_shuffle` DO generate replicate control sequences (n-axis)

This is the sharpest under-statement risk in the whole extraction. `assay_mpra`'s evidence
asserts "no notion of controls/scrambles as library members". But:

```python
def randomize(X, start, end, probs=[[0.25,0.25,0.25,0.25]], n=1, random_state=None)
    """... It will do this `n` times for each sequence in X and so return a
    tensor with one more dimension than `X`."""
```
(`ersatz.py` L343–360.) `shuffle(X, start, end, n=...)` and `dinucleotide_shuffle(X,
n_shuffles=...)` behave the same way. Dinucleotide shuffles **are** the standard MPRA
scramble control, tangermeme generates them in bulk, and it returns them as a materialized
`(batch, n, alphabet, length)` tensor of *sequences* — one of the very few places the
package hands you sequences rather than model outputs.

**What is still true:** they carry no identity, no name, no membership in a pool, and are
not composable with motif-placed members into one specification. Say that instead of
"no notion of controls/scrambles".

### 2.2 `marginalize_annotations` and `ablate_annotations` sweep a *set* of regions/motifs

`mixed_mutagenesis_one_pool`'s evidence says "marginalize (one motif)". That is incomplete.
`marginalize.py` L133 and `ablate.py` L171 define:

```python
marginalize_annotations(model, X, X0, annotations, **kwargs)   # annotations: (n_annotations, 3) = (example_idx, start, end)
ablate_annotations(model, X, annotations, **kwargs)
```

`marginalize_annotations` loops the annotation table, extracts `X[idx, :, start:end]` as the
motif, and marginalizes each one into the background set `X0` — i.e. an **n_motifs ×
n_backgrounds sweep in one call**, returning `PerturbationAnnotationsResult(y_befores,
y_afters)` stacked per annotation. This is genuinely the closest thing tangermeme has to a
per-element indexed sweep, and the extraction never mentions either function. Add them and
then say why they still don't make a pool: one perturbation *type* per call, position fixed
at centre, and the results are indexed only by row order in the user's own DataFrame.

### 2.3 `variant_effect.substitution_effect` accepts a heterogeneous per-example edit table

`substitution_effect(model, X, substitutions)` takes a COO tensor `(example_idx, position,
char_idx)`; rows are grouped by `example_idx` and all rows for one example are applied
together. The docstring adds: *"one can encode longer variants (e.g., entire motifs or just
multiple characters) by passing in multiple rows with adjacent positions."* So a **single
call** can carry example 0 = one SNV, example 1 = a 10-bp motif swap, example 2 = three
scattered edits. That is a real heterogeneous edit specification.

Still "no" for the row, because: substitutions only (indels and shuffles need separate
calls), no WT-replicate or sampled-random scheme, no record of which scheme produced which
row, and the output is `PerturbationResult(y_before, y_after)` — predictions, not a pool.
But the evidence must acknowledge the table or a referee will say the extractor never
opened `variant_effect.py`.

### 2.4 `greedy_substitution` internally materializes the full {motifs} × {positions} × {orientations} enumeration

The extraction says "No API enumerates {motifs} × {positions} × {orientations} into an
output sequence set" — true of the *output*, but the enumeration exists and is exhaustive.
`design/_substitute.py` is a numba kernel whose entire docstring is *"This function takes a
motif and inserts it at all possibilities"*, and `design/greedy_substitution.py` L188–202:

```python
if reverse_complement:
    motifs = motifs + [rc(motif) for motif in motifs]
...
for idx, motif in enumerate(motifs):
    input_idxs = torch.where(input_mask_ == True)[0].numpy()
    X_ = X.float().repeat(input_idxs.shape[0], 1, 1).numpy(force=True)
    _fast_tile_substitute(X_, motif_ohe, input_idxs)      # motif at EVERY allowed position
    y_hat = predict(model, X_, ...)
```

Every allowed position, every motif, both strands, materialized as a real sequence tensor —
then thrown away except for the argmin. Rephrase as: *the enumeration is a hidden internal
step of an optimizer; only the winner is returned, and there is no way to obtain the
enumerated set.* That is a stronger and un-attackable statement.

Relatedly, `assay_mpra`'s "it designs one sequence at a time" is only true of
`greedy_substitution`/`greedy_marginalize`. `beam_substitution(..., n_best=k)` returns
`(n_best, len(alphabet), length)` — k designed sequences, ranked — and `screen(..., n_best=k)`
does the same from random generation. Small point, but the author will notice.

### 2.5 The grep glob cited for six "no" verdicts does not cover `tangermeme/design/`

Six cells cite *"grep over tangermeme/\*.py (20 files)"* / *"all 20 modules"*. There are
indeed 20 `.py` files directly under `tangermeme/`, but the package also contains
`tangermeme/design/{__init__,_substitute,beam_substitution,greedy_marginalize,greedy_substitution,screen}.py`
and `tangermeme/_skills/`. A shell `tangermeme/*.py` glob misses all of them — including the
entire design subpackage, which is the most plausible place a synthesis or barcode helper
would live.

**I re-ran the grep recursively over the whole package (25 `.py` files + 15 `.md` skill docs)
and the result is identical: zero hits for barcode, primer, codon, hgvs, gtf, gff, transcript,
exon, intron, synthesis, homopolymer, restriction, amino, protein, orf, reading frame, oligo.**
The verdicts stand. Only the citation needs fixing — say "all 25 modules including
`design/`", or a referee gets a free "you didn't even look at our design package".

---

## 3. One factual error in cited evidence

`library_as_object` evidence states: *"tangermeme/results.py is the only shared-container
module and **all four of its NamedTuples** hold MODEL OUTPUTS, not sequences
(PerturbationResult(y_before, y_after); SpaceResult; PerturbationAnnotationsResult;
AttributionReferencesResult)."*

`results.py` contains **three** NamedTuples, not four:
`PerturbationResult`, `PerturbationAnnotationsResult`, `AttributionReferencesResult`.
`SpaceResult` is defined in `space.py` L22, and a fifth, `SaturationMutagenesisRawResult`,
is in `saturation_mutagenesis.py` L19. (The v1.2.0 changelog lists all five together, which
is probably where the miscount came from.)

Also strictly: `AttributionReferencesResult.references` **is** a tensor of sequences (the
DeepLIFT background shuffles), so "all … hold model outputs, not sequences" is not exactly
right. Neither point changes the value — no container is a library, none carries identity or
metadata — but both are the kind of detail the package author will check first.

---

## 4. Definitional risks (flag in the table footnotes, don't change the value)

- **`assay_dms = no`.** Defensible only under the reading "designs a DMS *library for
  synthesis*". tangermeme ships a module literally named `saturation_mutagenesis`, and
  `design.md` documents a second ISM path: *"pass `['A','C','G','T']` as `motifs` → greedy
  single-nucleotide (ISM-style) design"*. The distinction that saves the "no" is precise and
  should be in the footnote: **in silico only, nucleotide-level only, no codon/AA layer, and
  the mutant sequences are never returned** (`saturation_mutagenesis` builds `X_` internally
  and returns only aggregated attributions or `SaturationMutagenesisRawResult(y0, y_hat)`).
  Verified in `saturation_mutagenesis.py` L212–245.
- **`assay_mpra = partial`.** The extractor's own confidence note is right that a strict
  "library design for an MPRA experiment" reading makes this "no". Keep "partial" — it is the
  generous reading and cannot be attacked as understating a competitor — but state the row
  definition explicitly in the caption.
- **`library_as_object = partial`.** Generous. There is no class of any kind in the package
  (only NamedTuples). A "no" would also be defensible; "partial" is the safe choice and I
  would keep it.
- **`interface = yes` / `license = yes`.** The yes/no encoding loses the discriminating fact.
  In the rendered table these should read **"Python API only (no CLI, no GUI, no web
  service)"** and **"MIT"** — otherwise "interface: yes" is uninformative next to
  MutationMaker's web GUI + REST API and VaLiAnT's CLI. Facts themselves fully verified.

---

## 5. Capabilities the extraction missed entirely

1. **`marginalize_annotations` / `ablate_annotations`** — per-annotation perturbation sweeps
   returning `PerturbationAnnotationsResult`. See §2.2.
2. **`plot_logo(annotations=...)`** — the *static* logo already accepts a DataFrame with
   `motif_name, start, end, strand, score` and draws non-overlapping labelled boxes under the
   glyphs. The extraction cites only `interactive_logo` for `design_visualization`. This
   static path is the more relevant one: it means a designed construct *can* be rendered with
   its placed motifs named and positioned — if the user supplies that table themselves.
   (`plot.py` L142–160 and the `annotations` parameter docs.)
3. **Replicate-generating control primitives** (`randomize(n=)`, `shuffle(n=)`,
   `dinucleotide_shuffle(n_shuffles=)`) returning an extra sequence axis. See §2.1.
4. **`utils.example_to_fasta_coords`** — maps window-relative spans (seqlets, hits) back to
   genome coordinates using the originating BED. Mentioned in passing under
   `genome_coordinates` but it deserves a line: it is a genuine coordinate round-trip and the
   nearest tangermeme analogue to PoolParty emitting genome-anchored records.
5. **`extract_matching_loci` signal-aware negatives** — beyond GC matching it optionally takes
   a `bigwig` + `signal_beta` to reject candidate negatives with too much signal
   (`match.py` L372–395). Strengthens the "background/negative-set construction" suggested
   comparison row.
6. **`SaturationMutagenesisRawResult(y0, y_hat)`** — the `raw_outputs=True` escape hatch that
   returns per-mutant predictions rather than aggregated attributions. Relevant because it is
   the only way to get per-mutant *values* out, and it still does not give you the mutant
   *sequences*.
7. **`kmers` / `gapped_kmers`** (`kmers.py`) returning a sparse count matrix optionally
   weighted by per-position attribution scores — listed in the extraction's additional
   capabilities as "k-mer utilities" without noting the attribution-weighted mode, which is
   the non-obvious part.
8. **Optional-extra structure** (`[interactive]` → mpld3, `[docs]`, `[dev]` incl. captum) —
   trivial, but the `interactive` extra means `plot.interactive_logo` is *not* available in a
   default install. Worth a parenthetical if design_visualization is discussed.

---

## 6. Nothing found that would force a value change

I specifically hunted for, and did **not** find:

- any class, registry, or record type representing a set of sequences;
- any plugin, entry-point group, or optional dependency adding wet-lab, transcript, or
  nomenclature functionality (`[project.entry-points]` does not exist; the sole console
  script installs a Claude Code skill);
- any writer other than `one_hot_to_fasta`;
- any post-paper release adding library/pool/synthesis machinery (whats_new through v1.4.0
  reviewed in full: v1.3.0 = skill + plot + ISM fixes, v1.4.0 = beam search + design
  subpackage split, both purely in-silico);
- any docs-site content absent from the repo (readthedocs `latest` landing page matches
  `docs/index.rst`; tutorial/vignette list matches the tree exactly).

The extraction's "Bottom line" framing — *complementary, not competing; tangermeme's
composition target is a model call, not a library* — survives adversarial review and is the
right thing to say to a referee.
