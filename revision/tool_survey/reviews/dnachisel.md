# Adversarial review — DNA Chisel extraction

**Reviewer round:** document/web research only (no installation), 2026-08-10
**Target:** `/mnt/c/.../tool_survey/extractions/dnachisel.md` + the structured extraction
**Independent sources re-checked:** paper PDF (PyMuPDF), GitHub tree API for `master` (323 paths), raw source of `dnachisel/__init__.py`, `MutationSpace/MutationSpace.py`, `Location.py`, `utils/utils.py`, `biotools/{__init__,genbank_operations,random_sequences}.py`, `builtin_specifications/{AllowPrimer,AvoidHeterodimerization,EnforceChanges,EnforceSequence,EnforceRegionsCompatibility}.py`, `reports/optimization_reports.py`, `reports/constraints_reports/constraints_breaches_dataframe.py`, `cli.py`, `README.rst`, `changes.md`, `docs/{index,examples,notes}.rst`, `docs/ref/{core_classes,builtin_specifications}.rst`, `docs/genbank/genbank_api.rst`, all `examples/` scripts of interest, PyPI JSON API, GitHub repo/commits/branches API, live HTTP checks of the docs site and CUBA.

**Bottom line:** the extraction is careful and mostly right. Twenty-one of twenty-four rows survive. But it contains **two falsifiable absolute statements** that an author-referee (Zulkower/Rosser/Vegh know this codebase) could refute in one sentence each, and one of them forces a value change.

---

## FINDING 1 (most serious) — `lazy_evaluation = "no"` is falsified by a documented public API

The evidence says:

> "there are no generators/iterators over designed sequences because there is no multi-sequence output at all."

This is wrong. `dnachisel/MutationSpace/MutationSpace.py` defines:

```python
def all_variants(self, sequence):
    """Iterate through all sequence variants in this mutation space."""
    ...
    for variants in itertools.product(*variants_slots):
        new_sequence[choice_start:choice_end] = encoded_segment
        for (start, end), variant in variants:
            new_sequence[start:end] = variant
        yield new_sequence.decode()
```

That is a Python **generator** (`yield`) performing an `itertools.product` over a per-segment combinatorial choice space — the textbook definition of lazily enumerating a combinatorial design space. Alongside it:

- `@property space_size` — "Return the number of possible mutations", computed as a product of per-choice cardinalities **without materializing** anything (log/exp guarded against overflow).
- `pick_random_mutations(n_mutations, sequence)` / `apply_random_mutations(n_mutations, sequence)` — random sampling from the space.
- `MutationSpace.from_optimization_problem(problem, new_constraints=None)` — build the space from any problem's constraints (e.g. `EnforceTranslation` → the space of all synonymous variants of a CDS).
- `string_representation(separator="|")` and `plot(ax)` — text and matplotlib views of the whole allowed-variant space.

And this is **not** an undocumented internal: `docs/ref/core_classes.rst` contains

```rst
MutationSpace
-----------------------

.. automodule:: dnachisel.MutationSpace
   :members:
```

so `all_variants`, `space_size` and `apply_random_mutations` are all rendered on the public API reference page — the same page the extraction itself cites for `library_as_object`. The extraction even lists `MutationSpace/MutationChoice` among the core classes while simultaneously asserting no iterators exist.

**Fair counter-argument (keep it in the evidence):** `all_variants` yields bare sequence strings with no identity, metadata, or naming; it is a solver helper for exhaustive local search; it operates on the mutation space of one problem over one sequence; and no example, docs page prose, or the paper ever presents it as a way to produce a library. So this is not PoolParty-style lazy library materialization.

**Verdict: understated. Change `lazy_evaluation` to `partial`** with evidence along the lines of: *"`MutationSpace.all_variants(sequence)` is a documented public generator that lazily yields every variant of the constrained mutation space via `itertools.product`, with `space_size` giving the cardinality without materialization and `apply_random_mutations(n, seq)` sampling it. However it is a solver primitive: it yields unlabelled sequence strings for one problem's single sequence, is used internally for exhaustive local search, and no example or docs prose presents it as a way to emit a designed library."*

A flat "no" here is the single easiest thing for a DNA Chisel author to shoot down.

---

## FINDING 2 — `barcode_generation = "no"`: the stated reason is factually wrong

The evidence says the primer-collection example "has no edit-distance or demultiplexing notion". DNA Chisel ships a **documented, tested** example that is exactly minimum-Hamming-distance design over a set of short tags — `examples/builtin_specifications_examples/example_EnforceRegionsCompatibility.py` (test: `tests/builtin_specifications/test_EnforceRegionsCompatibility.py`):

```python
def compatibility_condition(location1, location2, problem):
    seq1 = location1.extract_sequence(problem.sequence)
    seq2 = location2.extract_sequence(problem.sequence)
    return sequences_differences(seq1, seq2) >= 2

problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(200, seed=123),
    constraints=[
        EnforceRegionsCompatibility(
            locations=[(0, 4), (50, 54), (100, 104), (150, 154)],
            compatibility_condition=compatibility_condition,
            condition_label="2bp difference"),
        EnforceGCContent(mini=0.4, maxi=0.6, window=40)])
```

`EnforceRegionsCompatibility.evaluate` runs `itertools.combinations(self.locations, 2)` and scores `-len(incompatible_pairs)`, i.e. it is a general all-pairs constraint over a set of subsequences with a user-supplied predicate. `sequences_differences` is Hamming distance and is a top-level export of `dnachisel`. `EnforceRegionsCompatibility` is listed in `docs/ref/builtin_specifications.rst`.

Add to that:
- `dnachisel.random_compatible_dna_sequence(sequence_length, constraints, probas=None, seed=None, ...)` — a **top-level exported** helper that generates a random sequence satisfying arbitrary constraints from scratch (`dnachisel/utils/utils.py`).
- The documented `primers_collection.py` recipe generates 20 short sequences under `AvoidHeterodimerization(existing, tmax=3)`, `AvoidPattern("3x3mer")`, `AvoidPattern("4xG")`, `EnforceGCContent(target=0.6)` — homopolymer/repeat filters and GC balance, the usual barcode filters.

So of the four things a barcode row normally means — (a) generate N short sequences, (b) enforce minimum pairwise distance, (c) apply GC/homopolymer/repeat filters, (d) attach them to library members — DNA Chisel has documented mechanisms for (a), (b) and (c), and only (d) is genuinely absent.

**Verdict: understated. Change to `partial`.** Suggested evidence: *"No barcode/UMI/demultiplexing abstraction exists ('barcode' appears nowhere in the source tree or docs), and there is no mechanism to attach generated tags to library members. But the primitives are documented and shipped: `random_compatible_dna_sequence(length, constraints)` generates constraint-satisfying sequences de novo; `EnforceRegionsCompatibility(locations, compatibility_condition)` scores all pairs of subregions under a user predicate, and the shipped example uses `sequences_differences(seq1, seq2) >= 2` — i.e. an explicit minimum-Hamming-distance constraint over four 4-nt tags; the primer-collection example adds the usual GC/homopolymer/repeat filters. Barcode-set design is therefore expressible but must be hand-assembled, and the result is a bare list of strings with no member association."*

If the matrix's `barcode_generation` is defined narrowly as "first-class barcode generator producing a demux-ready set bound to library members", keeping "no" is defensible — but **the evidence sentence must be rewritten regardless**, because "no edit-distance notion" is refutable by one shipped example.

---

## FINDING 3 — `mixed_mutagenesis_one_pool`: value survives, evidence does not

The evidence asserts: "There is no mutagenesis abstraction anywhere in the tool: no exhaustive-singles, exhaustive-pairs, random-sampling, or WT-replicate concepts."

Three of those four clauses are wrong at the primitive level:
- exhaustive enumeration: `MutationSpace.all_variants()`;
- random sampling: `MutationSpace.pick_random_mutations(n_mutations, seq)` / `apply_random_mutations(n_mutations, seq)`, plus `DnaOptimizationProblem.mutations_per_iteration`;
- an explicit "force at least N changes" primitive: `EnforceChanges(minimum=3)` / `EnforceChanges(minimum_percent=70)`, with GenBank shorthands `@change`, `@change(minimum=40%)`, `~change(3)` (documented in `docs/genbank/genbank_api.rst` under "Forcing changes"; `changes.md` v3.0.0 lists "AvoidChanges and EnforceChanges can now [be] tunable").

What is genuinely absent — and what the evidence should say instead — is the **pool**: no way to emit multiple mutagenesis regimes into one labelled output, no WT/control replicates, no per-variant records. The solver collapses the space to one sequence, and `all_variants` yields unlabelled strings from a single regime.

**Verdict: unsupported (keep the value "no", rewrite the evidence).**

---

## FINDING 4 — factual error in `automatic_naming` evidence

The evidence says that when `record_id` is not supplied, "the id is inherited from the input record". It is not. `to_record` builds a fresh record:

```python
record = sequence_to_biopython_record(self.sequence)
if record_id is not None:
    record.id = record_id
```

and `sequence_to_biopython_record(sequence, id="<unknown id>", name="<unknown name>", features=())`. So the default output id is the literal string `"<unknown id>"`. The conclusion (`automatic_naming = no`) is **strengthened**, but the sentence as written is wrong and must be fixed.

Also worth a footnote: `write_optimization_report` embeds `sequenticon(seq, output_format="html_image", size=24)` for the before/after sequences — an auto-generated *visual* identicon, not a name. Naming it in the evidence pre-empts a "we do auto-identify sequences" objection.

---

## FINDING 5 — the README advertises a collection-generating plugin

`README.rst` § "Related projects":

> `dnachisel-dtailor-mode <https://github.com/Lix1993/dnachisel_dtailor_mode>`_ brings features from D-tailor to DNA Chisel, **in particular for the generation of large collection of sequences covering the objectives fitness landscape** (i.e. with sequences which are good at some objectives and bad at others, and vice versa).

I checked it (repo has since moved to `HealthCodon/dnachisel_dtailor_mode`, repo id 233731735). It exposes `Design.DesignSpace`, `Design.FullFactorial`, `Design.Optimization`, `DnaDesignProblem`, `SequenceDesigner`, and an SQLite result store — a real multi-sequence design-space generator built on DnaChisel `Specification` objects.

Mitigating facts, all verified: created 2020-01-14, **last push 2020-01-15**, 3 stars, no PyPI package, no releases, Travis-CI badge (defunct service), third-party (not EGF).

This does **not** change `library_as_object = "no"` — it is not part of DNA Chisel and is abandoned. But the prompt specifically flags "documented plugin" as an attack vector, and this one is documented **in the tool's own README**. The `library_as_object` (and/or `dag_chaining`) evidence should pre-empt it with one clause: *"DNA Chisel's README points to a third-party `dnachisel-dtailor-mode` package for 'generation of large collection of sequences', but it is external, unreleased on PyPI, and unmaintained since January 2020."* Leaving this unaddressed is a gift to a referee.

---

## Smaller corrections and strengthenings (value unchanged)

**`library_as_object` — add the strongest available fact.** `load_record` uses `SeqIO.read` (not `parse`), so DNA Chisel's own loader accepts **exactly one record** and raises on a multi-record FASTA/GenBank. Input as well as output is single-sequence. Also the CLI signature is `dnachisel <source> <target>` — one file in, one out, with the docstring conceding "this CLI will be developed on a per-request basis". Both are stronger than the arguments currently used.

**`per_sequence_provenance` — undersold.** `write_optimization_report` emits, beyond the GenBank/PDF the extraction mentions: `constraints_before_and_after.csv` and `objectives_before_and_after.csv` (per-specification rows with `before`/`after` PASS-FAIL or score, `edits` count and `% edited` for that specification's span), an md5 hash of the source file plus a copy of the source file itself, `logs.txt` from the problem logger, and before/after sequenticons. This is machine-readable, per-specification, edit-attributed provenance — richer than "GenBank features + a PDF". "partial" still holds under a library-member reading, but say what is actually there, or an author will.

**`design_visualization` — undersold.** Add: GeneBlocks `DiffBlocks.from_sequences(before, after).plot(...)` renders the edit diff (`plot_optimization_changes`), plots use `dna_features_viewer` via `SpecAnnotationsTranslator`, and `MutationSpace.plot(ax)` / `string_representation()` visualize the allowed-variant space itself. The claim "there is no ... view of the design" is loose given `MutationSpace.plot`; say "no design-graph/pipeline view and no library-level view".

**`primer_design` — add a real limitation the extraction missed.** `primer3-py` is **not a runtime dependency**: PyPI `requires_dist` lists `primer3-py` only under `extra == "tests"`. `AvoidHeterodimerization.evaluate` raises `ImportError("Using avoid_heterodimerization requires primer3 installed")` on a default `pip install dnachisel`. Also `AllowPrimer` is a `SpecificationSet` expanding to `UniquifyAllKmers(k=max_homology_length)` + `EnforceMeltingTemperature(tmin,tmax)` + `AvoidPattern(RepeatedKmerPattern(k,n))` for `((2,5),(3,4),(4,3))` + optional `AvoidHeterodimerization` — worth stating precisely. And `AllowPrimer`/`AvoidHeterodimerization` cannot be declared as GenBank annotations (docs "Specifications not yet supported as Genbank annotations": `AvoidHeterodimerization`, `EnforceRegionsCompatibility`) — the extraction has this in `stated_limitations` but not in the primer row.

**`exon_intron_split_codons` — a stronger and more damning fact.** `Location.from_biopython_location(location)` reads only `location.start`, `location.end`, `location.strand`. A Biopython `CompoundLocation` therefore **silently collapses to its outer span**, so a spliced CDS would be treated as one contiguous region including its introns rather than raising. That is worse than "not supported" and is the sharpest evidence for this row.

**`codon_optimization` — modernize the provenance.** The paper cites Nakamura 2000 / Codon Usage Database, but the current package depends on EGF's `python_codon_tables`; species can be given by name (`e_coli`, `h_sapiens`, `s_cerevisiae`, …) **or by NCBI TaxID** (`species=1423`), in which case frequencies are downloaded from Kazusa. `EnforceTranslation` supports alternative genetic codes (`genetic_table=`) and non-ATG start-codon policies (`start_codon=keep|ATG`) — both added in v3.0.0 after the paper. Value "yes" is unaffected; the row currently cites 2020 provenance for 2025 software.

**`maintained` — value fine, wording slightly generous.** All figures verified exactly (v3.2.16 on 2025-05-10; last `master` commit 2025-05-10 09:40:51Z by veghp; `pushed_at` 2025-11-05; 278 stars, 53 forks, 22 open issues; not archived; ten branches incl. `dev`, `palindromic`, `issue53`). But as of the survey date there has been **no commit to `master` for ~15 months and no release for ~15 months**. "Actively maintained" plus "release cadence remained steady" reads as more momentum than the record shows; "maintained, last release 2025-05-10, no commits since" is safer and still supports "yes".

**Documented-examples list — minor inaccuracies.**
- `docs/examples.rst` puts `ecoli_genes_optimization.py` (not `genome-wide-optimization.py`) under the heading "Programmatic optimization of a full genome". The extraction sources the genome-wide script to `docs examples.html`; the docs actually include the other script. The repo file exists, so the substance is right, but the citation is off.
- `sequence_without_restriction_sites.py`, `plasmid_optimization.py` and the manuscript examples are **repo-only**, not in `docs/examples.rst`. The extraction's `source` fields imply otherwise for the domestication row.
- "A per-specification example directory illustrating each of the 15+ built-in specification classes individually" overstates `examples/builtin_specifications_examples/`: it holds 12 scripts covering ~9 specs (`AvoidBlastMatches`, `AvoidChanges`, `AvoidHairpins`, `AvoidMatches`, `AvoidPattern`, `CodonOptimize` ×2, `EnforcePattern`, `EnforceRegionsCompatibility`, plus `custom_specification`, `gene_PCR-tagging`, `max_random_iters`) — not one per class.
- Two examples are missing from the list, one of them material: **`example_EnforceRegionsCompatibility.py`** (see Finding 2) and `example_gene_PCR-tagging.py` (which despite its name demonstrates `AvoidBlastMatches` against an E. coli BLAST DB).

**Paper-derived quotes — all verified.** "An optimization problem is defined … against which a starting linear or circular sequence will be optimized", the two-phase solver, "an insert(CGTCTC) constraint will attempt to place the pattern 'CGTCTC' at different locations", "produces comprehensive optimization reports for traceability and troubleshooting", Fig. 1A/1B/1C legends, "Python library, a command-line interface, or a web application" — all present verbatim in the PDF. The Supplementary S1A–S1E caveat in `confidence_notes` is honest and correct (the supplement is not in `papers/`).

**Live checks re-run 2026-08-10:** `https://cuba.genomefoundry.org/sculpt_a_sequence` → 200 (Vue SPA shell, `<title>EGF CUBA</title>`, app bundles served — matches the extraction's caveat exactly); `https://edinburgh-genome-foundry.github.io/DnaChisel/` → 200; PyPI `license_expression: MIT`, latest 3.2.16.

---

## Capabilities the extractor missed entirely

1. **`MutationSpace` public API** — `all_variants()` (lazy generator over the combinatorial space), `space_size`, `pick_random_mutations` / `apply_random_mutations`, `from_optimization_problem`, `string_representation()`, `plot(ax)`; documented via `automodule` in `docs/ref/core_classes.rst`. Bears on `lazy_evaluation`, `mixed_mutagenesis_one_pool`, `combinatorial_motif_place`, `design_visualization`.
2. **`random_compatible_dna_sequence(length, constraints, probas, seed, max_random_iters)`** — top-level export that *generates* a constraint-satisfying sequence from scratch. This contradicts the memo's framing of DNA Chisel as "repair-and-optimize rather than generate"; the README's own first example also starts from `random_dna_sequence(10000)`. The additional_capabilities bullet should be reworded.
3. **`EnforceRegionsCompatibility` with arbitrary user predicates** — all-pairs constraint across a set of subregions; the shipped example is a min-Hamming-distance tag set. Comparable to PoolParty-style cross-member constraints, and the closest thing DNA Chisel has to a set-level design primitive.
4. **README-advertised `dnachisel-dtailor-mode`** — full-factorial / fitness-landscape *collection* generation on top of DnaChisel specs (external, 3 stars, dead since 2020-01-15, no PyPI).
5. **Report provenance artifacts** — per-spec before/after CSVs with edit counts and `% edited`, source-file md5 + embedded copy, `logs.txt`, sequenticons, GeneBlocks diff figure.
6. **Documented reproducibility contract** — `docs/notes.rst` guarantees that setting the NumPy seed before any DnaChisel call makes results deterministic ("If you have reproducibility issues despite setting the seed, we would consider it a bug"), backed by `tests/test_determinism/`. If the PoolParty matrix or rebuttal claims reproducibility as a differentiator, DNA Chisel has an explicit one.
7. **`load_record` = `SeqIO.read`** — single-record input only; the loader itself refuses multi-record files. Strong support for `library_as_object = no`.
8. **`primer3-py` is not a runtime dependency** (only the `tests` extra) — `AvoidHeterodimerization` and the heterodimer half of `AllowPrimer` raise `ImportError` on a default install. Belongs in `stated_limitations` next to Bowtie/BLAST.
9. **Codon-table modernization** — `python_codon_tables`, species-by-name or NCBI TaxID with Kazusa download; alternative genetic codes and start-codon policies (post-paper, v3.0.0).
10. **Adoption evidence** — README: "Benchling uses DNA Chisel as part of its sequence optimization pipeline according to this webinar video". Useful if the rebuttal characterizes tool maturity, and dangerous to omit if a referee cites it.
11. **`EnforceChanges(minimum=…)` / `@change(minimum=40%)`** — an explicit "force at least N (or X%) nucleotide changes" primitive with GenBank shorthand. The nearest thing to a mutagenesis knob and currently mentioned only in passing.

---

## Verdict summary

| key | claimed | verdict | action |
|---|---|---|---|
| library_as_object | no | supported | add `SeqIO.read` single-record input + pre-empt the dtailor-mode README pointer |
| dag_chaining | partial | supported | generous but defensible; "no" would also be defensible |
| lazy_evaluation | no | **understated** | **→ partial** (`MutationSpace.all_variants` is a documented public generator) |
| mixed_mutagenesis_one_pool | no | **unsupported** | keep "no", rewrite evidence (random sampling + exhaustive enumeration + `EnforceChanges(minimum=)` all exist) |
| combinatorial_motif_place | no | supported | soften "search states" wording given `all_variants` |
| barcode_generation | no | **understated** | **→ partial** (`EnforceRegionsCompatibility` + `sequences_differences>=2` shipped example; `random_compatible_dna_sequence`) |
| per_sequence_provenance | partial | supported | add before/after CSVs, md5, logs, sequenticon |
| automatic_naming | no | supported | **fix factual error**: default id is `"<unknown id>"`, not inherited |
| design_visualization | partial | supported | add GeneBlocks diff, dna_features_viewer, `MutationSpace.plot` |
| assay_dms | no | supported | — |
| assay_mpra | no | supported | — |
| assay_insilico | no | supported | — |
| genome_coordinates | no | supported | — |
| transcript_models | no | supported | add `SeqIO.read` detail |
| exon_intron_split_codons | no | supported | add `from_biopython_location` silently collapsing `CompoundLocation` |
| hgvs_input | no | supported | — |
| vcf_vep_output | no | supported | — |
| consequence_annotation | no | supported | — |
| primer_design | partial | supported | add primer3-py-not-installed-by-default; spell out the `AllowPrimer` spec set |
| codon_optimization | yes | supported | update provenance to `python_codon_tables` + TaxID |
| synthesis_constraints | yes | supported | — |
| interface | yes | supported | note the CLI is one-file-in/one-file-out by design |
| license | yes | supported | — |
| maintained | yes | supported | soften "actively"; note ~15 months with no commit or release |
