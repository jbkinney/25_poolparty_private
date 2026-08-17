# tangermeme — evidence memo for PoolParty tool comparison

**Surveyed:** 2026-08-10 (document/web research only; nothing installed)
**Tool:** tangermeme, Jacob Schreiber (IMP Vienna / UMass Chan)
**One-line self-description (pyproject/PyPI):** "Biological sequence analysis for the modern age."
**Paper self-description:** "a highly optimized toolkit for 'everything-but-the-model' when it comes to genomic deep learning"

---

## 1. Sources consulted

| Kind | Reference |
|---|---|
| PDF (local) | `/mnt/c/.../papers/Schreiber2025_tangermeme_all.pdf` — bioRxiv 2025.08.08.669296, v2 posted 2025-08-23 |
| Prior analysis | `/mnt/c/.../prior_analyses/Schreiber2025_tangermeme_all_analysis.md` |
| Repo | https://github.com/jmschrei/tangermeme (README, file tree via GitHub API, commits, repo metadata) |
| Repo source | `tangermeme/ersatz.py`, `product.py`, `space.py`, `io.py`, `utils.py`, `match.py`, `annotate.py`, `plot.py`, `results.py`, `saturation_mutagenesis.py`, `design/{__init__,screen}.py`, `pyproject.toml`, `docs/whats_new.rst`, `docs/index.rst` |
| Bundled agent-skill docs (in-repo, authoritative prose) | `tangermeme/_skills/data/references/{func-pattern,motif-effects,product,design,io-loci,variant-effect}.md` |
| Docs site | https://tangermeme.readthedocs.io/en/latest/ |
| PyPI | https://pypi.org/pypi/tangermeme/json (v1.4.0) |

**Verdict on the prior analysis:** broadly correct in its framing ("no notion of a variant library as a structured object"), but it is **incomplete and in two places wrong/outdated**:
- It says design/library-relevant machinery is absent. In fact tangermeme ships a whole `tangermeme.design` subpackage (`screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize`) and a `tangermeme.product` module with Cartesian-product application. Fig 1A of the paper explicitly lists "Design" and "Screening" columns including "Cartesian Products", "Motif Implantation", "Construct Marginalization".
- It says "Sequence metadata: None" — true for designed sequences, but note `tangermeme.seqlet` / `tangermeme.annotate` do produce per-region DataFrames (example_idx, start, end, p-value, motif annotation). That is *analysis* metadata, not construction provenance.
- It was written against the ~v0.5/1.0 era; the package is now at **v1.4.0 (2026-06-25)** with beam search and a bundled Claude Code Agent Skill.

---

## 2. What tangermeme actually is

From the paper (p.3): *"we introduce tangermeme, a software package which implements 'everything-but-the-model' when it comes to genomic deep learning. … tangermeme focuses on using trained models to drive genomics research (Fig 1A), including support for sequence manipulations, model operations, and efficient combinations of the two (e.g. variant effect prediction)."*

Fig 1A capability matrix (tangermeme column) lists: Greedy Substitution, Motif Implantation, Construct Marginalization, Ledidi (Design); Sequence Manipulation, Stacked Operations, Cartesian Products, Saturation Mutagenesis, Marginalizations, Ablations, Spacings, Variant Effect Prediction (Screening); Seqlet Calling/Annotation/Ablation+Marg./Annotation Counting/Spacing (Seqlets); DeepLIFT/SHAP, PISA, Plotting (Core). It explicitly does **not** claim Model Training, Data Preprocessing, or a Model Zoo.

Data model: **every sequence collection is a bare `torch.Tensor` of shape `(-1, len(alphabet), length)`.** There is no library, pool, or record class. `tangermeme/results.py` — the only "container" module — defines four NamedTuples, and every one of them holds *model outputs*, not sequences:

```python
class PerturbationResult(NamedTuple):
    y_before: torch.Tensor | list[torch.Tensor]
    y_after:  torch.Tensor | list[torch.Tensor]
```

---

## 3. Capability-by-capability evidence

### BLOCK A — library specification

**library_as_object → partial.**
Operations are batch-first, so a user does not loop to apply *one* operation to many sequences: e.g. `substitute(X, motif)` where *"If a motif with a batch size equal to that of X is provided, there will be 1-1 correspondence between the motifs and the sequence"* (`ersatz.py` docstring). But the batch is a raw tensor: no class, no per-sequence identity, no composition metadata, no heterogeneous membership. Building a library with several different design intents requires the user to build and `torch.cat` tensors themselves and to track what is what outside the library. `results.py` confirms the only shared return types are model-output pairs.

**dag_chaining → partial.**
There is a genuine, deliberate composition mechanism — the `func=` plug-point. From `_skills/data/references/func-pattern.md`: *"tangermeme's perturbation functions are built around one contract so any analysis can be swapped for any other. … Any function passed as `func=` must satisfy: `func(model, X, args=None, **kwargs) -> torch.Tensor | list[torch.Tensor]`"*, accepted by `ablate`, `marginalize`, `space`, `variant_effect.*`, `product.apply_pairwise`, `product.apply_product`. Nesting is documented and demonstrated:

```python
apply_pairwise(marginalize, model, X, args=(cell_states,),
               additional_func_kwargs={'func': deep_lift_shap}, motif="TGACGTCA")
```
(`product.md`; also Tutorial B7). The paper calls this "Stacked Operations" (Fig 1A) and argues (p.3) *"An important consequence of tangermeme's design is that operations can be stacked."*

**But**: what is composed is *sequence-manipulation → model-operation*, i.e. an analysis pipeline whose leaf is always a model call. It is Python function nesting, not a reified graph: there is no DAG object, no node/edge inspection, no reuse of a node in two branches, and the composition does not yield a library of sequences (it yields tensors of predictions/attributions).

**lazy_evaluation → partial.**
`product.md`: *"Batches are built iteratively, so the full product is never materialized in memory."* Confirmed in `product.py`: `for x in tqdm(itertools.product(X, *args))` accumulating up to `batch_size` before calling. `io.extract_loci` has `max_jitter` expansion explicitly to *"reduce the memory footprint"*. Counter-evidence: `ersatz.*` returns fully materialized tensors, and `saturation_mutagenesis` eagerly materializes every single-base mutant of one example at a time (`X_ = X[i].repeat(n_edits, 1, 1)`). So streaming exists, but only over the *model-input* product axis, never as a library abstraction.

**mixed_mutagenesis_one_pool → no.**
Each entry point implements exactly one perturbation type and returns a homogeneous result: `marginalize` (one motif), `ablate` (one interval, n shuffles), `space` (one motif list across a spacing grid), `saturation_mutagenesis` (all single substitutions), `variant_effect.substitution_effect` (a supplied substitution table). There is no specification object in which exhaustive singles + exhaustive pairs + sampled randoms + WT replicates coexist, and nothing in the API records which scheme produced which sequence. Nothing in README, docs, tutorials, or the bundled skill references describes such a mixture.

**combinatorial_motif_place → partial.**
Real combinatorics exist along two axes:
- *Multiple motifs + spacings*: `multisubstitute(X, motifs, spacing, start=None)` places a list of motifs with per-gap spacings; `space(model, X, motifs, spacing)` takes `spacing` of shape `(n_spacings, n_motifs-1)` — *"Each row in this tensor is a different combination of spacings between motifs and each column is the spacing between an adjacent pair of motifs"* — and returns `SpaceResult(y_before, y_afters)` with `y_afters` shape `(batch, n_spacings, ...)`.
- *Motif × position × orientation search*: `greedy_substitution` / `beam_substitution` — *"Each round, tries substituting every motif at candidate positions and keeps the single edit that most reduces `loss`"*; *"`reverse_complement=True` also considers each motif's reverse complement"*; `input_mask` restricts editable positions (`design.md`).

**Limits that make this "partial" not "yes":** `substitute`/`multisubstitute`/`space` take a single scalar `start` — position is *not* a swept axis; you would loop yourself. `apply_product` crosses `X` with **extra model inputs** (`args`), not with motifs/positions/orientations — `product.md` is explicit that args are model inputs (cell state, read depth). And the greedy/beam routines are *optimizers* that return one (or `n_best`) designed sequence(s), plus `greedy_substitution` *"`X` must have batch size 1 — it designs one sequence at a time."* No API enumerates {motifs} × {positions} × {orientations} into an output sequence set.

**barcode_generation → no.**
Grepped all 20 modules for `barcode`, `edit distance`, `primer`: zero hits for barcode/primer; the single "edit distance" hit is `saturation_mutagenesis.py` docstring *"sequences with an edit distance of one on them"* (i.e. single mutants, not barcode separation). GC machinery does exist but for a different purpose: `utils.gc_content` ("Compute the GC content of one-hot encoded sequences") and `match.py` ("the calculation of GC-content genome-wide and the sampling of GC-matched negatives"). Nothing generates or attaches barcodes.

**per_sequence_provenance → no.**
Outputs are tensors. `results.py` NamedTuples carry only `y_before`/`y_after`/`attributions`/`references`. The only DataFrames produced are `seqlet.recursive_seqlets` (called seqlets: example_idx, start, end, attribution, p-value) and `annotate.annotate_seqlets` (motif index + p-value) — these describe *discovered* features in attribution scores, not how a designed sequence was constructed.

**automatic_naming → no.**
The only writer is `io.one_hot_to_fasta(X, filename, mode='w', headers=None, alphabet=...)`: *"If headers are provided for each sequence, these are used, otherwise the numeric index is used."* No name is derived from the design operation.

**design_visualization → partial.**
`tangermeme.plot` provides `plot_logo`, `interactive_logo`, `plot_attributions`, `plot_pwm`, `plot_categorical_scatter`. v1.3.0 changelog: *"Adds `plot.interactive_logo` … Annotations are drawn as translucent, pastel boxes (colored by any `annot_cmap`) behind the logo glyphs, with the motif name in the box corner and a hover tooltip listing the length and every column of the annotation (e.g. seqlet p-value, annotation p-value, summed attribution)."* So individual sequences *can* be rendered with highlighted regions. There is no view of a design/library specification, no graph view, no library-level summary plot.

### BLOCK B — assay coverage

**assay_dms → no.** `saturation_mutagenesis` is explicitly the *in silico* analogue: README — *"This is another form of attribution method that is conceptually similar to deep mutational scanning but using a predictive model instead of running an experiment."* But it operates on nucleotides only, returns predictions/attributions rather than a sequence set, and the package contains **no** codon, amino-acid, ORF or reading-frame concept (grep: zero hits for `codon` across all modules). Nothing supports designing a coding/DMS library for synthesis.

**assay_mpra → partial.** `tangermeme.design` genuinely designs cis-regulatory sequences against a model oracle (`screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize`); Tutorial B6 designs a sequence with strong AP-1 binding using Beluga; the paper's Fig 1A "Design" row is real. However there is nothing MPRA-specific: no barcodes, adapters, primer sites, oligo-length/synthesis constraints, no notion of controls/scrambles as library members, and the greedy designers work one sequence at a time (`X` batch size 1). This is *regulatory sequence design*, not *MPRA library design*.

**assay_insilico → yes.** This is the tool's raison d'être. `marginalize`, `ablate`, `space`, `saturation_mutagenesis`, `variant_effect.{substitution,deletion,insertion}_effect`, `apply_pairwise`/`apply_product`, `deep_lift_shap`, `pisa`, all with built-in batching, device handling and multi-input/multi-output support. Paper p.3: *"These operations have built-in batching, support for alternative data types and devices, and work out-of-the-box on multi-input/output models."*

### BLOCK C — genomics integration

**genome_coordinates → yes.** `io.extract_loci(loci, sequences, signals=None, in_signals=None, chroms=None, in_window=2114, out_window=1000, max_jitter=0, ...)` — *"Either the path to a bed file or a pandas DataFrame object containing three columns: the chromosome, the start, and the end"*, with `sequences` a FASTA path (pyfaidx) and `signals` bigwigs (pybigtools). Also exclusion-list filtering, chromosome subsetting, `utils.example_to_fasta_coords`, and `match.py` for genome-wide GC-matched negative sampling.

**transcript_models → no.** `io.py` defines exactly three public readers: `read_meme`, `read_vcf`, and the BED/FASTA/bigwig path inside `extract_loci`. Grep across all modules for `gtf`, `gff`, `transcript`, `exon`, `intron`: zero hits.

**exon_intron_split_codons → no.** Same grep; no exon/intron/codon concepts anywhere.

**hgvs_input → no.** Grep for `hgvs`: zero hits. Variants are specified as integer tensors, e.g. README: `substitutions = torch.tensor([[0, 1058, 0]])` (example index, position, character index). `io.read_vcf` can bring in a VCF as a DataFrame, but the user must convert it to indices themselves.

**vcf_vep_output → no.** `read_vcf` is input-only (*"returns a pandas DataFrame with the comments filtered out … only the first 9 columns"*). The only writer in the package is `one_hot_to_fasta`. No VCF/VEP emitter.

**consequence_annotation → no.** `tangermeme.annotate` exists but "annotation" here means *motif matching of seqlets via TOMTOM* (`annotate_seqlets`, `count_annotations`, `pairwise_annotations`, `pairwise_annotations_spacing`) and, more generally, *"an annotation is any genomic span"* (README). There is no molecular-consequence (stop-gained / synonymous / frameshift) logic.

### BLOCK D — physical construction

**primer_design → no.** Zero hits for `primer` across the source tree; no wet-lab protocol support anywhere in README/docs/paper.

**codon_optimization → no.** Zero hits for `codon`.

**synthesis_constraints → no.** No homopolymer, restriction-site, repeat, secondary-structure, or length-constraint checking. `utils.gc_content` measures GC and `match.py` matches GC between peak and background regions — these serve background selection for marginalization (*"Uniformly random sequences are much higher in GC content than real genomic DNA … substituting a motif into the wrong background can produce a 'mirage' effect"*, `motif-effects.md`), not synthesizability screening.

### BLOCK E — engineering

**interface → Python API (library) only.** No web service, no GUI. The one console script is `tangermeme-install-skills` (`pyproject.toml [project.scripts]`), which only installs the bundled Claude Code Agent Skill. Analysis CLIs were removed: README warning — *"These FIMO and Tomtom command-line tools have been moved to memesuite-lite … Please use those!"* Notable extra surface: a bundled **Claude Code Agent Skill** (`SKILL.md` + 14 reference files) so an LLM agent can drive the API.

**license → MIT.** `pyproject.toml`: `license = "MIT"`; GitHub API `license.spdx_id = "MIT"`; paper Code Availability: *"tangermeme is free and open source software available under the MIT license at https://github.com/jmschrei/tangermeme."*

**maintained → yes, actively.** Last commit on `main`: **2026-07-19** ("ADD code of conduct"); previous substantive commit 2026-06-25 (v1.4.0 release). Latest PyPI release **1.4.0, 2026-06-25**. Release cadence over the last year: 1.0.0 (2025-08-27), 1.0.2 (2026-01-19), 1.0.3 (2026-02-15), 1.0.4 (2026-04-23), 1.1.0/1.2.0 (2026-05-27), 1.3.0 (2026-06-23), 1.4.0 (2026-06-25). Repo not archived; 308 stars, 32 forks; CI unit-test badge; Read the Docs live.

---

## 4. tangermeme's own documented examples (candidates for PoolParty reproduction)

Tutorials (`docs/tutorials/`, all rendered on the docs site):

| Notebook | Content |
|---|---|
| Tutorial A1 — Ersatz Sequence Manipulation | substitute / multisubstitute / insert / delete / randomize / shuffle / dinucleotide_shuffle / reverse_complement, incl. inserting JASPAR AP-1 (MA1144.1) into random background |
| Tutorial A2 — Predictions | batched `predict` |
| Tutorial A3 — DeepLIFT/SHAP | attributions |
| Tutorial A4 — Seqlets | recursive seqlet caller + TF-MoDISco caller |
| Tutorial A5 — Annotations | seqlet annotation, counting, pairwise, spacing |
| Tutorial B1 — Marginalization | motif effect on predictions/attributions in backgrounds |
| Tutorial B2 — Ablation | shuffle-out of a window |
| Tutorial B3 — Spacing | motif-pair spacing sweeps (AP-1 cooperativity fading with distance) |
| Tutorial B4 — Saturation Mutagenesis | in silico ISM |
| Tutorial B5 — Variant Effect | substitution / deletion / insertion effects |
| **Tutorial B6 — Design** | **greedy_substitution and beam_substitution against Beluga to design a sequence with high AP-1 binding, starting from one random sequence + a JASPAR motif subset** |
| **Tutorial B7 — Cartesian Product** | **`apply_pairwise` / `apply_product`, including nesting `marginalize(func=deep_lift_shap)` inside a product** |
| Tutorial C1 — IO and Data Loading | `extract_loci`, BED/FASTA/bigwig |
| Tutorial C2 — Plotting | logos, attributions |

Vignettes (`docs/vignettes/`):
- **Inspecting What Cis-Regulatory Features a Model Has Learned** — README: *"If you only read one vignette, read THIS ONE"*
- Attribution Trickiness and DeepLiftShap Implementations
- Wrappers are Productivity Hacks

How-to: *How To — Reduce Friction and Save Time with Tangermeme*.
Paper figure notebooks: `docs/paper/Fig1_Timings_Examples.ipynb`, `docs/paper/Fig2_Seqlet_Calling_and_Downstream_Analyses.ipynb` (PLD6 promoter, BPNet/Beluga/ChromBPNet/ProCapNet).

**Reproduction targets most relevant to PoolParty:** Tutorial B3 (spacing sweep = a small combinatorial motif-placement library) and Tutorial A1 (motif implantation into shuffled/GC-matched backgrounds) are the two use cases a PoolParty DAG would plausibly express as a library, with the model call left to tangermeme. Tutorial B6 (design) is *search*, not enumeration — PoolParty would not replace it.

---

## 5. Capabilities not covered by the row list

- DeepLIFT/SHAP with built-in batching and correctness guarantees on convergence deltas (paper Fig 1B: *"these implementations can silently fail, causing the convergence deltas to become quite large instead of near zero … whereas tangermeme will never do so"*).
- PISA (`pisa.py`), saliency-style attribution alternatives.
- **Novel recursive, variable-length seqlet caller** — *"This seqlet caller is the first to call variable-length seqlets directly"* — plus a TF-MoDISco-style caller.
- Seqlet annotation against a motif database via tomtom-lite, then `count_annotations`, `pairwise_annotations`, `pairwise_annotations_spacing` (motif co-occurrence and spacing statistics).
- GC-matched genome-wide negative/background sampling (`match.py`).
- Extremely fast one-hot encoding (paper Fig 1C: 1.61 s for hg38 chr1, ~3× the next fastest).
- k-mer utilities (`kmers.py`), `utils.entropy` / `information_content` / `pwm_consensus`.
- Alphabet-agnostic throughout — README: *"although the library was built with operations on DNA sequences in mind, all functions are extensible to any alphabet"*.
- Model-oracle sequence design: `screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize` (the last returns a variable-width one-hot *construct* with `N` where motifs cancel).
- Interactive mpld3 sequence logos with hover tooltips.
- Bundled Claude Code Agent Skill for LLM-driven use.
- Reads bigwig signal tracks, MEME motif files, VCF.

**Suggested new comparison rows this tool would motivate (if the referees want them):** "model-oracle sequence optimization / design search", "attribution & interpretation of a trained model", "background/negative-set construction (GC- or dinucleotide-matched)".

---

## 6. Stated limitations

- Scope is deliberately "everything-but-the-model": the paper's own Fig 1A marks Model Training, Data Preprocessing and Model Zoo as **not** provided by tangermeme (unlike Selene, EUGENe, gReLU, CREsted). A trained PyTorch model is a prerequisite for most of the package.
- Design is discrete and greedy; `design.md` states: *"`tangermeme.design` is discrete and greedy. When you want gradient-based, minimal edits to a specific template sequence … reach for the `ledidi` library."*
- `greedy_substitution` requires `X` batch size 1 and a positive `max_iter` (`max_iter=-1` is documented as *"a silent no-op"*).
- The recursive seqlet caller assumes constant sequence length: *"In principle, l does not need to be constant across examples but in our current implementation it is assumed to be fixed."*
- FIMO/TOMTOM CLIs no longer ship with tangermeme (moved to memesuite-lite).
- `read_vcf` drops genotype columns past column 9 and does not support BCF.
- Attribution-based workflows are documented as expensive (marginalize + deep_lift_shap ≈ 200 forward/backward passes for 5 examples).

---

## 7. Bottom line for the PoolParty comparison

tangermeme is the Tier-1 tool with the most genuine overlap in *sequence manipulation primitives* — `ersatz.substitute`/`multisubstitute` really do implant motifs, `space` really does sweep spacing combinatorics, and `func=`/`apply_product` really do compose operations. The overlap ends at the abstraction boundary: **tangermeme's composition target is a model call, not a library.** Every pipeline terminates in predictions or attributions; sequences are anonymous tensors with no identity, no provenance, no names, and no way to hold heterogeneous design intents in one specification. It has no wet-lab layer at all (no barcodes, primers, codons, synthesis checks) and no transcript/variant-nomenclature layer (no GTF, HGVS, VCF output, consequence calls). The honest framing for the referee response is *complementary, not competing*: PoolParty specifies and materializes the library; tangermeme is one natural consumer of it.
