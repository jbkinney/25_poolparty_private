# tangermeme — repair change log

**Record repaired:** `final/tangermeme.md` (edited in place)
**Audits adjudicated:** `citation_audit/tangermeme.md` (19 findings), `factcheck/tangermeme.md` (A1–A9, B1–B17, C)
**Repair date:** 2026-08-14

## Verification method

Every finding was checked against primary source before any edit. Primary sources used:

- **Repository at the audited commit.** `jmschrei/tangermeme` tarball for
  `8d732b8c08764057b7ae5faa644d48664f36b44b` downloaded and extracted to
  `<scratch>/tmrepo/tangermeme-8d732b8c…/`; all line numbers below were read from that tree.
- **Paper.** `papers/Schreiber2025_tangermeme_all.pdf` (13 pp.) extracted with PyMuPDF.
- **Live endpoints.** `api.github.com/repos/jmschrei/tangermeme`,
  `api.github.com/repos/jmschrei/tangermeme/issues?state=open`,
  `raw.githubusercontent.com/…` (HTTP status checks).
- **Local workspace files** for the extraction memo and adversarial review.

**No capability value was changed.** All 24 `key → value` lines are byte-identical to the
pre-repair record (re-verified by grep after editing). Where a verified correction bears on
what a value rests on, the finding was escalated rather than applied.

---

## APPLIED (19)

### Citation audit

**CA-3 / FC-A9 — §8's "only sequences it ever hands back" list is false.**
*Verified:* `ersatz.insert` returns `torch.cat([X[:,:,:start], motif, X[:,:,start:]], dim=-1)`
(`ersatz.py` L94, documented return `(-1, len(alphabet), length + motif_length)`);
`substitute` returns `X` (L195); `multisubstitute` L297; `delete` L340. All four hand back
sequence tensors and all four are named approvingly two paragraphs earlier in the same section.
*Correction:* the clause now reads "the sequences it does hand back — the products of the
`ersatz` edits (`insert`, `substitute`, `multisubstitute`, `delete`), DeepLIFT backgrounds,
generated controls, and design winners — are all anonymous tensors…". "Every pipeline
terminates in predictions or attributions" was narrowed to "Every model-operation pipeline".

**CA-6 — "applies to six 'no' cells" is not reproducible.**
*Verified by enumeration.* In the **final record**, nine `no` cells rest on the recursive
negative keyword grep: `barcode_generation`, `assay_dms`, `transcript_models`,
`exon_intron_split_codons`, `hgvs_input`, `consequence_annotation`, `primer_design`,
`codon_optimization`, `synthesis_constraints`. (Excluded: `per_sequence_provenance` and
`automatic_naming`, which cite the *class* grep and the *writer* grep respectively, and
`vcf_vep_output`, which cites the writer grep.) In the **extraction memo**, seven `no` cells
invoke the narrow grep (`extractions/tangermeme.md` L77, 90, 100, 102, 104, and the
`primer_design` / `codon_optimization` entries). Six is reproducible from neither.
*Correction:* "six" → "nine" (the count that matches the record's own closing instruction,
"use the '28 modules…' citation everywhere").
*Not fixed, noted:* the string the record attributes to the extraction —
*"grep over `tangermeme/*.py` (20 files)"* — is not literal. The extraction says "Grepped all
20 modules". The record inherited the paraphrase from the adversarial review (`reviews/…` §2.5),
which also presented it as a quotation. Neither audit raised this; fixing it is out of scope.

**CA-7 — the five NamedTuples are not all "of model outputs".**
*Verified:* `results.py` L41–52 — `AttributionReferencesResult(attributions, references)`,
where `references` is documented as "the reference sequences it used". Confirmed against the
full-tree class grep (six classes, unchanged).
*Correction, two sites:* §1 now reads "exactly six classes, five of which are `NamedTuple`s
(four of model outputs; the fifth, `AttributionReferencesResult`, pairs attributions with the
reference *sequences* used)". §7 item 1 now reads "defines no class representing a library, and
no class carrying identity or metadata for a sequence (… the one that does hold sequences,
`AttributionReferencesResult`, holds DeepLIFT backgrounds)".
*Residual tension, deliberately left:* §1's closing line still says "no class of any kind
representing a sequence…". That sentence was not among the audit's cited locations, and the
record's own Block-A phrasing ("none represents a library, and none carries identity, name, or
metadata for any sequence") already states the accurate version.

**CA-8 / FC-A1 — not every sequence collection is 3-D.**
*Verified:* `ersatz.randomize` returns `torch.stack(X_rands).permute(1,0,2,3)` (L422), documented
`shape=(-1, n, len(alphabet), length)` (L394); `shuffle` L498 / L471; `dinucleotide_shuffle`
L688 / L665. The record already gives that 4-D shape twice elsewhere.
*Correction:* §2 Data model now reads "a bare `torch.Tensor` — shape `(-1, len(alphabet),
length)`, or `(-1, n, len(alphabet), length)` for the replicate-generating `ersatz.randomize` /
`shuffle` / `dinucleotide_shuffle`".

**CA-11 — the greedy-substitution code block is annotated paraphrase, not source.**
*Verified against `design/greedy_substitution.py`:* the loop is at L190
(`for idx, motif in enumerate(tqdm(motifs, disable=not verbose)):`, not `tqdm(motifs, ...)`);
`predict` at L202–203 has its arguments elided in the record; `pos = loss_curr.argmin()` is at
**L214**, outside the cited L188–205; reverse-complement augmentation is at **L161–162**; the
block omits L193–195 and L200; and none of the three inline comments exists in the file. The
enumeration behaviour itself is real (`_fast_tile_substitute`, `design/_substitute.py`, whose
whole docstring — *"This function takes a motif and inserts it at all possibilities"* — the
record quotes verbatim).
*Correction:* the introducer now reads "`design/greedy_substitution.py` L186–214 (abridged
below, with the explanatory comments added here and the reverse-complement augmentation shown
at its real location, L161–162)". The block itself was left untouched.
*Residual tension, deliberately left:* the `*Source:*` line four lines below still cites
L188–205.

**CA-12 — DataFrame inventory wrong; `annotate_seqlets` returns tensors.**
*Verified:* `annotate.py` L89 returns `torch.from_numpy(idxs).type(torch.int32),
torch.from_numpy(p_values)`; its docstring Returns block documents two tensors; only the
`-> tuple[pandas.DataFrame, list[str]]` annotation at L25 says otherwise. Actual DataFrame
producers, by `pandas.DataFrame(` grep: `seqlet.py` L353 (`tfmodisco_seqlets`), `seqlet.py`
L576 (`recursive_seqlets`), `match.py` L601 (`extract_matching_loci`); plus
`utils.example_to_fasta_coords` (`-> pandas.DataFrame`, L238, returns `coords_df` L335) and
`io.read_vcf` (`-> pandas.DataFrame`, L622).
*Correction:* the sentence was replaced with the corrected inventory plus a parenthetical
noting the stale annotation. The cell's argument ("never how a sequence was constructed") is
preserved.

**CA-13 — `one_hot_to_fasta` is not the package's only writer.**
*Verified:* the record's own grep pattern is accurate (`.write(` hits are `io.py` L580, 582,
585, 587 only), but `_skills/install.py` L100–103 calls `shutil.rmtree`, `Path.mkdir` and
`shutil.copytree`. Nothing else in the tree writes.
*Correction, two sites:* `automatic_naming` — "the package's **only writer**" → "the package's
**only sequence writer**". `vcf_vep_output` — "The exhaustive writer search … finds file-writing
code in **`io.py` alone**" → "The writer search … finds sequence/data-serializing code in
**`io.py` alone**". Both cells' conclusions (integer-index fallback header; no VCF writer) are
untouched and remain correct.

**CA-14 — `interactive_logo` ships; only its dependency is optional.**
*Verified:* `plot.py` L578 defines `interactive_logo`; `plot.py` has no top-level `mpld3`
import; the function does `from mpld3 import plugins` inside a `try/except ImportError` that
re-raises with "interactive_logo requires mpld3. Install it with `pip install mpld3` or
`pip install tangermeme[interactive]`."
*Correction, two sites:* `design_visualization` — "absent from a default install" → "The
function itself ships and imports in a default install (`plot.py` L578); it imports `mpld3`
lazily and **raises `ImportError` when called** unless the extra is present." §6 limitation
reworded to match.

**CA-16 — the `extract_loci` quotation is ~45 lines beyond the cited range.**
*Verified:* `def extract_loci(` at `io.py` L246, `) -> tuple:` at L265 — so L246–265 is the
signature, correctly cited. The quoted `loci` parameter doc ("Either the path to a bed file or a
pandas DataFrame object containing three columns…") is at L310–315, verbatim.
*Correction:* "(`io.py` L246–265)" → "(signature at `io.py` L246–265; the quoted `loci`
parameter doc at L310–315)".

**CA-17 — the eight "open issues" include three PRs.**
*Verified live (2026-08-14):* `open_issues_count: 8`; enumerating
`/issues?state=open` returns 8 items, of which 3 carry a `pull_request` key (#80, #54, #42) and
5 do not (#49, #48, #27, #25, #3). Stars 308 and forks 32 re-confirmed.
*Correction:* row label → "Stars / forks / open issues+PRs"; value → "**308 / 32 / 8**
(`open_issues_count`; enumerating gives 5 issues + 3 PRs)".

**CA-18 — k-mer return type and gapped weighting described incorrectly.**
*Verified:* `kmers.py` L68–70 builds `X_kmers = torch.zeros((X.shape[0], n**k))` and returns it
— a **dense torch.Tensor**; the `-> scipy.sparse.csr_matrix` annotation at L22 is stale (the
docstring Returns block correctly says `torch.Tensor`). `kmers`' weighting **is** a sum
(`conv1d` with a width-`k` ones kernel, L63–64). `gapped_kmers` **does** return a CSR (L250–253)
but accumulates `new_gkmer_attr / k` (L136) — the **average**, contradicting its own docstring.
*Correction:* the §5 bullet was rewritten with the corrected return types, corrected line
anchors, and the sum-vs-average distinction.

**CA-19 — `N` does not mark "where motifs cancel".**
*Verified:* `greedy_marginalize.py` L242–246 — `Xn = ['N' for _ in range(X.shape[-1])]`, then
`Xn[pos:pos+len(motif)] = list(motif)` for each chosen motif in order, then `''.join(Xn).strip("N")`.
Internal `N`s are therefore uncovered gaps; overlapping motifs are overwritten by the later
assignment. The bundled `design.md` L96–97 does make the "cancel" claim, so the record faithfully
repeated a documentation error.
*Correction:* the §5 bullet now describes the actual mechanism and flags the docs/implementation
divergence.

### Fact-check section A

**FC-A2 — `input_mask` restricts starts, not edited positions.**
*Verified verbatim* in both `design/greedy_substitution.py` L92–97 and
`design/beam_substitution.py` L111–116: *"A mask on input positions that can be the start of
substitution. Any motif can be substituted in starting at each allowed position even if the
contiguous span of the mask is shorter than the motif."* The record's source for the old claim
was the bundled `design.md` L45–46, which is the imprecise one.
*Correction:* "`input_mask` restricts editable positions (`design.md`)" → "`input_mask` restricts
the positions at which a motif may *start*, not the positions it may overwrite — a motif
substituted at an allowed start runs over masked-out positions downstream of it
(`design/greedy_substitution.py` L92–97)".

**FC-A3 — "Only the argmin survives" is greedy-specific.**
*Verified:* `greedy_substitution.py` L214 `pos = loss_curr.argmin()`. `beam_substitution.py`
L246–271 builds the same `{motif} × {allowed start}` candidate tensor per beam member, prunes
with `torch.topk` into a size-`beam_size` heap, and carries `beam_size` sequences forward. The
record's second clause ("no way to obtain the enumerated set") is **true of both** — the
candidate tensor `X_` is discarded in both algorithms, and `beam` holds complete sequences, not
the enumeration.
*Correction:* "Only the argmin survives" → "Only the winner of each round survives (the argmin
in `greedy_substitution`, the `beam_size` best in `beam_substitution`)". The load-bearing second
clause was left intact because it verifies.

**FC-A5 — multi-input/multi-output support is not universal.**
*Verified:* `deep_lift_shap.py` L255–259 — *"NOTE: predictions MUST yield a
`(batch_size, n_targets)` tensor… If your model yields something more complicated you must wrap
the model"*; `saturation_mutagenesis.py` L198–201 raises `ValueError("raw_outputs=True is
required for models that return multiple output tensors…")`; `predict.py` L86–94 does handle
lists. The record had restated a paper-level generalisation (paper p.3, quoted verbatim and
correctly) as a verified property of every listed operation.
*Correction:* "all with built-in batching, device handling and multi-input/multi-output support"
→ "all with built-in batching and device handling. Multi-input/multi-output support is
per-operation rather than universal: [two cited counterexamples]".
*Residual tension, deliberately left:* the paper quote immediately following still makes the
broad claim — but it is explicitly attributed to the paper, so the juxtaposition reads as
"record finding, then paper claim" rather than a contradiction.

**FC-A7 — "alphabet-agnostic throughout" overstates.**
*Verified:* `utils.gc_content` L896–915 takes an `alphabet=` argument but is DNA-semantic ("Used
only to locate G and C"); `utils.reverse_complement` L507–511 defaults to
`{"A":"T","C":"G","G":"C","T":"A"}` and accepts a single 2-D tensor;
`greedy_substitution` / `beam_substitution` both default `reverse_complement=True`. The README
quote itself is verbatim (`README.md` L10).
*Correction:* a clause was added to the bullet's **existing** parenthetical Note. The bullet
header and the README attribution were left alone — the record already framed the claim as the
README's, not as a verified finding.
*Not applied:* the audit's further claim that `greedy_marginalize` fails to forward a custom
alphabet on replay is **true** (`chosen_motifs` holds strings via `best_motif = motifs[idx]`
L217/L224, and the replay at L193 calls `substitute(X_, motif, start=pos)` with no `alphabet=`),
but it is an implementation-bug observation beyond the scope of this record.

**FC-A8 — "Design is discrete and greedy".**
*Verified:* the quoted sentence is verbatim at `_skills/data/references/design.md` L119, so the
record reports the bundled docs accurately — but `design/__init__.py` `__all__` exports
`screen`, `greedy_substitution`, `beam_substitution`, `greedy_marginalize`, and `screen` samples
independent random batches with no edit trajectory at all (`screen.py` L37–48, L155–185). The
record already lists all four elsewhere, so the §6 header was internally inconsistent.
*Correction:* the §6 lead clause now reads "Design is discrete and combinatorial; the bundled
docs call it greedy, though `design/__init__.py` also exports the non-greedy random `screen`."
The verbatim quote is unchanged.

**FC-A4 / CA-1 — "Nucleotide-level only" is not established by the grep.**
*Verified:* `saturation_mutagenesis.py` L222–232 derives candidate edits from
`identity = (ref == 1) & (ref.sum(dim=0, keepdim=True) == 1)` and `torch.where(~identity)` over
`ref = X[i, :, start:end]`, i.e. it enumerates whatever alphabet axis `X` carries — no
four-channel hard-coding anywhere. The grep result itself (zero `codon` / `amino` / `protein` /
`orf` / `reading frame`) is reproducible and unchallenged.
*Correction, two sites:* Block-B `assay_dms` fact 2 relabelled "**No coding-biology layer**" with
the channel-generic implementation noted and the grep's scope stated correctly; §7 item 3
reworded to match.
*Why this was applied rather than escalated:* `assay_dms = no` is stated to rest on three facts;
facts 1 (in silico only, no synthesis output) and 3 (mutant sequences never returned — confirmed
at `saturation_mutagenesis.py` L245–250) are each independently dispositive under the row
definition *"designs a DMS library for synthesis"*, and both are untouched. The value does not
move.

### Fact-check section B (only two of seventeen applied — see ESCALATED)

**FC-B1 — `screen` is random de novo generation, not an edit search.**
*Verified:* `screen.py` L19–36 signature (`func: Callable = random_one_hot`,
`additional_func_kwargs`, `n_best`, `max_iter=-1`); docstring L37–48 — *"randomly generate a
batch of examples and evaluate them using the provided model, keeping only the `n_best` top hits…
each batch is independent from the others"*; implementation L155–185 pushes into a size-`n_best`
heap with no de-duplication.
*Why applied rather than escalated:* the record's §5 bullet paired `screen(n_best=k)` with
`beam_substitution(n_best=k)` as if both were template-edit searches, which is inaccurate by
juxtaposition rather than merely incomplete.
*Correction:* a clause distinguishing the two modes was appended to that existing bullet.

**FC-B16 — the paper benchmarked v0.5.1, not the surveyed v1.4.0.**
*Verified in the PDF, Methods p.6:* "We compared the speed of one-hot encoding the entirety of
hg38 chr1 using several software packages. These included tangermeme (**v0.5.1**), SeqPro
(v0.6.1), CREsted (v1.4.0), gReLU (v1.0.7), and Selene SDK (0.6.0)."
*Why applied rather than escalated:* the record cites the Fig 1C timing as a present-tense
property of the surveyed version, which is a citation-integrity problem rather than a coverage
judgment.
*Correction:* the §5 timing bullet now carries "— benchmarked on tangermeme **v0.5.1**, paper
Methods p.6, not the v1.4.0 surveyed here".

---

## REJECTED (2)

**CA-4 — "three abbreviated local source paths do not resolve" (status: wrong-location).**
*Why wrong:* the three files exist and are at the location the record indicates —
`prior_analyses/Schreiber2025_tangermeme_all_analysis.md`, `extractions/tangermeme.md`,
`reviews/tangermeme.md`, all confirmed present under
`/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/`. The `/mnt/c/.../`
form is an ellipsis of a prefix the **immediately preceding table row spells out in full**
(row 1 gives the complete path to the PDF). "wrong-location" is not the right status for a
correctly-located, deliberately abbreviated path, and no reader is misdirected. Record left
untouched.

**CA-5 — "the raw GitHub directory URL is dead" (status: dead-link).**
*The 404 is real* — `https://raw.githubusercontent.com/jmschrei/tangermeme/main/` returns 404
(checked). *But it identifies no error in the record:* line 27 does not cite the prefix as a
page consulted; it states that individual files were "downloaded from
`raw.githubusercontent.com/jmschrei/tangermeme/main/`", i.e. as the path stem for per-file
fetches. Files below the prefix resolve normally (`.../main/README.md` → 200, checked).
raw.githubusercontent has never served directory listings, so a bare-prefix 404 is expected
behaviour, not a broken citation. Record left untouched.

---

## ESCALATED (5 items + 2 groups)

**E1 — `interface` render string: "no CLI" is false, but the render string is the cell value.**
(CA-2 and FC-A6, same finding.)
*Verified:* `pyproject.toml` L71–72 declares `tangermeme-install-skills =
"tangermeme._skills.install:main"`; `_skills/install.py` L122–131 builds an
`argparse.ArgumentParser(prog="tangermeme-install-skills", …)` with `--dest`, `--force`,
`--print-path`. A CLI exists. The record's own body already says so.
*Question for the authors:* the `interface` row's mandated render string —
**"Python API only (no CLI, no GUI, no web service)"** — *is* the table cell content for this
row (§7 says this cell "must not be rendered as bare yes/no"). Changing it changes the cell,
which the repair pass is not permitted to do. How should it read?
*Options:* (a) "Python API only (no analysis CLI, no GUI, no web service)" — smallest change,
keeps the discriminating contrast with MutationMaker/VaLiAnT; (b) "Python API; one
installation-only CLI (no analysis CLI, GUI, or web service)" — the fact-check's proposal, more
literal but longer in a table cell; (c) leave as-is and put the console script in a footnote.
Same decision applies to the echo in §7. Record left untouched at both sites.

**E2 — `lazy_evaluation`: "Streaming exists only on the Cartesian-product axis" is false.**
(CA-9.)
*Verified:* `match.py` contains real streaming generators — `_sequence_generator_stream`
(L63–67, docstring *"Returns a generator, so that only one sequence is kept in memory at a
time"*), `_count_generator_stream` (L119), `_chrom_coords_generator` (L22),
`_resize_coords_generator` (L31), consumed via `numpy.fromiter(..., count=num_regions)` at
L197, L245, L281. The exclusivity claim does not hold.
*Question for the authors:* does genome-scanning streaming in `match.py` count toward the
`lazy_evaluation` row, or does the row mean lazy evaluation *of a designed library* (in which
case the record's trailing clause "never as a library abstraction" already carries the argument
and only the word "only" needs to go)? This determines whether the correction is a one-word
edit or a re-weighing of the `partial`. Record left untouched.

**E3 — `combinatorial_motif_place`: "the only genuinely enumerated axis in the package".**
(CA-10.)
*Verified, partially:* `saturation_mutagenesis` does expose an enumerated
{character} × {position} grid in its output (`y_hat` reshaped to
`(X.shape[0], X.shape[1], end-start, …)`, L258–261), and the record itself documents a second
enumeration three paragraphs later (design/'s internal {motifs} × {positions} × {fwd, rc}). So
"only" is not sustainable as written. *The audit is however wrong about one of its two
counterexamples:* the `n` axis of `randomize` / `shuffle` / `dinucleotide_shuffle` is a
replicate count of random draws, not an enumeration of anything.
*Question for the authors:* what does "genuinely enumerated axis" mean for this row —
user-specified and output-exposed (which would keep `space` unique among motif-placement
operations), or any exhaustive materialised grid (which would not)? Any wording I choose here
decides a definition inside the cell. Record left untouched.

**E4 — "Dinucleotide shuffles are *the* standard MPRA scramble control" is uncited.**
(CA-15.)
*Verified:* the cited `ersatz.py` source establishes the shuffle mechanics and the return shape,
nothing about MPRA experimental practice. The claim is about the field, and neither the
tangermeme paper nor the repository can settle it, so I could not verify it from primary source
and will not restate it in my own words either.
*Question for the authors:* cite an MPRA methods reference, or hedge to "a standard scramble
control in MPRA work"? This sits inside `assay_mpra`, which §7 flags as "generous by design", so
the strength of the wording is an emphasis decision. Record left untouched.

**E5 — fact-check section B, items 2–15 and 17 (fifteen coverage findings).**
Each asks the record to say more about a real capability; I spot-verified several and found the
underlying facts sound (e.g. `greedy_marginalize` optimises `f(X+construct) − f(X)` over a
background *set*, `greedy_marginalize.py` L157–210; `utils.chunk`/`unchunk` exist at
`utils.py` L649–786 and are absent from the record; `shuffle` applies one positional permutation
per replicate across the whole batch; `io.one_hot_to_fasta` calls `characters` without exposing
`force`/`allow_N`). But each is a judgment about whether an omission is *material* to a
capability comparison — explicitly reserved for the authors — and applying all fifteen would add
roughly fifteen sentences of fresh unaudited text to a record whose §5 already runs to ~25
bullets.
*Question for the authors:* which of B2–B15 and B17 should be added? My own ranking of what most
affects a library-design comparison: **B2** (`greedy_marginalize`'s background-averaged mode),
**B10** (`chunk`/`unchunk`), **B11** (FASTA export fails on ambiguous/all-zero columns — which
interacts with `greedy_marginalize`'s `N`-containing construct and so touches `automatic_naming`),
**B12** (`read_meme` → `pwm_consensus` → `characters` is the supported PWM-to-design-motif path),
**B15** (matched loci are not one-to-one paired with inputs). B17 (installation commands) reads
to me as out of scope for a capability record; that is a scope call, not a fact call.
Record left untouched for all fifteen.

**E6 — fact-check section C (balance).**
The audit judges the record "not fully balanced" — absent pooled-library and wet-lab features
documented at greater depth than the four present design algorithms. That is precisely the
emphasis judgment reserved for the authors. Its one concrete, checkable sub-claim (the §8
absolute about returned sequences) was verified and fixed under CA-3/FC-A9 above. Record
otherwise left untouched.

---

## Residual roughness introduced (deliberate, per the surgical-editing rule)

1. `combinatorial_motif_place`: the code-block introducer now says L186–214 while the
   `*Source:*` line four lines below still says L188–205.
2. `assay_insilico`: the record's own per-operation finding now sits immediately before the
   paper's broader (correctly attributed) quote.
3. §1's closing "no class of any kind representing a sequence…" was outside the audited
   locations for CA-7 and is unchanged, though §1's opening and §7 item 1 now carry the
   `AttributionReferencesResult` exception.
4. `assay_dms` fact 2's heading changed from "Nucleotide-level only" to "No coding-biology
   layer"; the surrounding "three verified facts" framing is unchanged and still holds.

## Snapshot note (not an audit finding, no edit made)

The record's live-figures table is dated to its own pass. Re-fetching on 2026-08-14 returns
`pushed_at: 2026-08-14T14:22:12Z` and `updated_at: 2026-08-14T14:23:29Z`, later than the
2026-07-19 / 2026-08-07 values recorded. Stars (308), forks (32), `archived: false` and
`open_issues_count` (8) are unchanged. The record presents these as a dated snapshot, and
refreshing them was not in either audit, so they were left alone.

---

## Pass 2 — policy application

**Date:** 2026-08-14

**Baseline and counts.** Every open factual item was rechecked against the supplied paper and
official `jmschrei/tangermeme` repository/docs at commit
`8d732b8c08764057b7ae5faa644d48664f36b44b`; release/version checks used GitHub Releases and
PyPI. The PDF was read directly with PyMuPDF. Counts use the 21 atomic open policy items below:
**5 applied · 16 declined-by-policy · 0 rejected-unverifiable · 0 unresolved record
escalations.** One downstream Table-2 value implication is reported separately and was not
applied. FC-B1 and FC-B16 were already closed in Pass 1 and are not recounted.

**Version/source check.** The header's surveyed release, **1.4.0**, remains the current PyPI and
GitHub release (2026-06-25), so Policy D requires no current-release parenthetical. The record
already contains its single material paper-snapshot qualification: the Fig. 1C benchmark used
tangermeme 0.5.1. The header version continues to govern.

| Open item | Outcome | Edit and primary-source verification |
|---|---|---|
| E1 / CA-2 / FC-A6 — `interface` render string | **applied** | Kept `interface = yes` and changed both render instructions to **“Python API; one installation-only CLI (no analysis CLI, GUI, or web service)”**. `pyproject.toml:71–72` declares only `tangermeme-install-skills`; `_skills/install.py:108–149` implements that installer with `--dest`, `--force`, and `--print-path`; the README says the former FIMO/Tomtom analysis CLIs moved out. The corrected evidence supports the existing value. |
| E2 / CA-9 — `lazy_evaluation` scope | **declined-by-policy (row definition)** | Verified both sides: `product.py` iterates and batches the Cartesian product without materializing it, while `match.py:22–145` separately streams genome coordinates, FASTA slices, and bigWig counts. The row asks about lazy evaluation of a designed library, so genome scanning does not enter its test; the existing `partial` already rests on `product.py` and “never as a library abstraction.” No edit or value change. |
| E3 / CA-10 — “only genuinely enumerated axis” | **applied** | Replaced the undefined uniqueness claim with **“a user-specified, output-exposed enumerated spacing axis.”** `space.py:25–181` accepts a 2-D row set of spacing combinations and exposes `n_spacings` in `y_afters`; `saturation_mutagenesis.py:212–269` separately exposes a character × position grid, and greedy/beam design materialize motif × position candidates internally. The existing `combinatorial_motif_place = partial` remains supported. |
| E4 / CA-15 — “the standard MPRA scramble control” | **applied** | Dropped the field-wide anecdote and retained only the verified package fact: the three control functions return sequences rather than model outputs. `ersatz.py:343–688` verifies the mechanics and four-dimensional return tensors; neither the paper nor official repository/docs establishes an MPRA-practice standard. `assay_mpra = partial` still rests on the documented design and control APIs. |
| E5 / FC-B2 — background-averaged `greedy_marginalize` | **declined-by-policy** | Verified in `design/greedy_marginalize.py:22–61,157–246` and Tutorial B6: it optimizes `f(X + construct) - f(X)` over background sequences and searches motif/orientation/spacing choices. This adds detail but changes neither locked row 9 (`partial`) nor row 12 (`yes`), and it does not correct a Table-1 Purpose, Key features, Output, or Availability cell. |
| E5 / FC-B3 — design objectives and restrictions | **declined-by-policy** | Verified in `design/greedy_substitution.py:24–239` and `beam_substitution.py:25–306`: `y=None`, custom loss, masks, single-character mode, fixed-length template substitution, and beam-only multiple winners are real. They do not change locked row 12's `yes`, row 9's `partial`, or a Table-1 cell. |
| E5 / FC-B4 — true insertion and deletion | **applied (evidence clause)** | Added one sentence to the existing `assay_insilico` entry: atomic `ersatz.insert` / `delete` return longer / shorter edited tensors at one shared batch coordinate and do not search positions. Verified in `ersatz.py:20–94,300–340`. This directly bears on locked row 8 under its operational definition; the score implication is reported below and was not applied. |
| E5 / FC-B5 — paired/masked substitution semantics | **declined-by-policy** | Verified in `ersatz.py:97–297`: motifs broadcast or pair by batch, all-zero/default-`N` columns preserve the original, and `multisubstitute` places a supplied arrangement. These semantics do not change row 9's `partial` because the arrangement is supplied rather than combinatorially enumerated, nor any Table-1 cell. |
| E5 / FC-B6 — control independence/diversity | **declined-by-policy** | Verified in `ersatz.py:425–688`: each `shuffle` replicate uses one positional permutation across the batch; dinucleotide shuffles preserve endpoints and can warn/fail for insufficient diversity. Locked row 4 expressly excludes shuffled controls, so no score or Table-1 cell changes. |
| E5 / FC-B7 — batched `predict` details | **declined-by-policy** | Verified in `predict.py:17–165`: batching, source-dtype retention, auxiliary inputs, post-processing, tensor/list outputs, and model-state restoration are implemented; auxiliary inputs are cast to model dtype. This does not change model-guided row 12 (`yes`), and Table 1 already says the output includes model predictions. |
| E5 / FC-B8 — indel-effect normalization | **declined-by-policy** | Verified in `variant_effect.py:124–404`: deletion inputs include padding, insertion outputs are flank-trimmed, and the functions return before/after model outputs. Row 8 is already determined by the direct sequence-returning `ersatz` operations in B4; these scoring restrictions do not alter it or Table 1. |
| E5 / FC-B9 — pairwise versus Cartesian product | **declined-by-policy** | Verified in `product.py:30–349`: `apply_pairwise` index-zips auxiliary contexts while `apply_product` independently crosses them. Those axes are model inputs, not independent mutation events under locked row 5; composable row 2 is already `yes`. No Table-1 cell changes. |
| E5 / FC-B10 — `chunk` / `unchunk` | **declined-by-policy** | Verified in `utils.py:649–786`: variable-length inputs are tiled into fixed windows, incomplete terminal windows are discarded, and `unchunk` currently needs original lengths. No locked row or Table-1 cell measures window tiling. |
| E5 / FC-B11 — FASTA ambiguity restriction | **declined-by-policy** | Verified in `io.py:540–587` and `utils.py:338–410`: `one_hot_to_fasta` calls `characters` without `force` or `allow_N`, so ties/all-zero columns can fail. This export edge does not change a locked score or Table 1's accurate output description (“Perturbed sequences and model predictions”). |
| E5 / FC-B12 — MEME/PWM-to-motif path | **declined-by-policy** | Verified in `io.py:590–619` and `utils.py:789–822`: `read_meme` loads PWMs and `pwm_consensus` makes an argmax discrete motif, choosing the earlier channel on ties and preserving all-zero columns. It does not change row 9's `partial` or a Table-1 cell. |
| E5 / FC-B13 — diagnostics and seeding | **declined-by-policy** | Verified in `utils.py:874–1016`: `gc_content`, entropy/information content, and `set_seed` behave as described, but none repairs a sequence to meet a constraint and the finding therefore does not satisfy locked row 11. No Table-1 cell changes. |
| E5 / FC-B14 — seqlet restrictions | **declined-by-policy** | Verified in `seqlet.py:211–577`: both public callers take projected 2-D attribution tracks; `recursive_seqlets` calls positive seqlets and prevents core overlap, while `tfmodisco_seqlets` is only the calling stage. No locked library-design row or Table-1 cell measures seqlet discovery. |
| E5 / FC-B15 — matched-background output semantics | **declined-by-policy** | Verified in `match.py:372–621`: the function returns sorted BED coordinates, explicitly not one-to-one input matches; high-N loci may be dropped and scarce GC bins spill to neighbours. Locked row 16 remains `yes` because it requires coordinate input/output, not pairing, and Table 1 is unchanged. |
| E5 / FC-B17 — installation commands | **declined-by-policy** | Verified in the README installation section and `pyproject.toml:5–35`: `pip install tangermeme` / `uv add tangermeme`, Python `>=3.10`, Torch `>=2.0`. Table 1 already says “Python package (PyPI),” so Availability does not change. |
| E6 / fact-check section C — balance | **declined-by-policy** | Policy A declines record-level balance and proportional-emphasis edits. No record text changed. |
| Policy C — sibling-source provenance | **applied** | Removed the Prior analysis, Extraction memo, and Adversarial review rows from §1's source table. The remaining evidence sources are only the paper, official repository/docs, GitHub repository metadata, and PyPI. Editorial-history wording and correction labels remain as history, not as factual support. |

**Value-basis checks and downstream escalation.** No value in `final/tangermeme.md` changed:
primary evidence still supports `interface = yes`, `lazy_evaluation = partial`,
`combinatorial_motif_place = partial`, and `assay_mpra = partial` after their prose corrections.
However, locked Table-2 row 8 is presently `?` for tangermeme. Under its operational definition,
the verified public `ersatz.insert` **and** `ersatz.delete` operations support **`Insertions and
deletions = ●`** (both true lengthening and shortening are returned). Per the no-value-change rule,
this **`? → ●` value implication is escalated and reported, not edited**.

**Rejected-unverifiable findings:** none. The only unverifiable underlying assertion was the
“standard MPRA scramble control” anecdote; Policy C required dropping the assertion, which was
done. Every audit finding used for an edit was independently verified in an admitted primary
source.

**Remaining record escalations:** none. **Downstream value escalations:** the single locked-row-8
score implication above.

**Row-substitution candidates:** none established. The locked table already includes tangermeme's
best discriminating library-design capability, `Model-guided variants` (tangermeme `●`, PoolParty
`◐`). The open `Recombination and chimeras` and `Codon / amino-acid substitutions` uniformity risks
cannot be resolved from one tool's evidence, and the suggested attribution/interpretation axis lacks
the required verified eight-tool comparison; no substitution is nominated.

**Neighbouring tension logged.** §1 still describes the extraction/review workflow and entries
retain `[CORRECTED]` labels, but those sibling documents are no longer listed or used as evidence.
Table 2 remains `?` for row 8 until the separate scoring authority acts. The existing Table-1 phrase
“insertion and removal of motifs” also remains broader/less explicit than the verified distinction
between true `ersatz.insert`/`delete` and the substitution-based oracle-guided designers.
