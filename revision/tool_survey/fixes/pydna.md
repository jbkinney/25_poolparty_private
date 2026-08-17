# pydna — repair change log

Record repaired: `revision/tool_survey/final/pydna.md` (edited in place)
Audits processed: `citation_audit/pydna.md` (19 findings), `factcheck/pydna.md` (A:7, B:19, C:1 = 27 findings)
Date: 2026-08-14

## Primary sources used for verification

Every finding below was checked against primary source, not against the audit's assertion.

| Source | How obtained |
|---|---|
| `pydna-group/pydna` @ `4e02f81dfd39629c68d2a93dfbb434bfe7d56120` (master, 2026-08-01) | full tarball via `codeload.github.com`, read locally |
| `pydna-group/pydna` @ tag `v5.5.16` | full tarball, read locally |
| PyPI `pydna` 5.5.16 wheel | downloaded, inspected inside the zip (this is the *released* artefact, not a git tag) |
| Sphinx `searchindex.js` | `https://pydna-group.github.io/pydna/searchindex.js` (HTTP 200, 257,700 bytes) |
| GitHub REST API | `/repos/pydna-group/pydna` and `/search/issues` (issue vs PR split) |
| Paper | `papers/Pereira2015_pydna.pdf` extracted with PyMuPDF |
| Live URL checks | `code.google.com/p/pydna/`, its archive JSON, `pydna-shell.appspot.com`, old GitHub URL |

**Capability values: none changed.** All edits are to evidence and prose only.

---

## APPLIED (21)

### Citation audit

**CA-1 — "no collection type appears" among the 820 documented objects.** APPLIED.
Verified: `class RecombinaseCollection` at `src/pydna/recombinase.py:391`; `class OpenDNACollectionsSource` at
`src/pydna/opencloning_models.py:732`. Both appear in the hosted `searchindex.js` (12 and 8 hits
respectively), so both are in the cited 820-object documented inventory. The absolute claim is false as
written.
Correction: narrowed the claim to *sequence* library/pool/collection classes and named the three types that
do exist (`RecombinaseCollection` = a set of recombinase enzymes; `OpenDNACollectionsSource` = a repository-ID
provenance record; `CloningStrategy` = one strategy's sequences/primers/steps), stating that none is a
designed set of sequences. Summary-table row and the "Stated limitations" bullet narrowed the same way.
`library_as_object = no` untouched — the value rests on the absence of a *variant-library* abstraction, which
survives the correction intact.

**CA-2 — "`CloningStrategy` is export-only; its entire surface is …".** APPLIED.
Verified in `opencloning_models.py`: `to_dseqrecords()` (L1322), `normalize()` (L1388), `validate()` (L1395),
`get_ids_of_sequences_that_are_inputs()` / `…_are_not_inputs()`. `to_dseqrecords` is in `searchindex.js`, so
it is documented API. `to_dseqrecords()` walks the dependency graph and rebuilds `Dseqrecord`s **from the
sequences embedded in the strategy itself** (`seq_by_id` is built from `self.sequences`), so the record's
operative reasoning — no re-execution on *new* inputs — is correct and was preserved.
Correction: replaced "export-only / its entire surface is" with a round-trip description naming the import
methods, then kept the original "no re-execution … on new inputs" sentence unchanged. Summary-table row
updated to match. `dag_chaining = partial` untouched.

**CA-3 — `utils.cai()` quotation invented.** APPLIED.
Verified: `src/pydna/utils.py:225-226` docstring is literally `"""docstring."""`. The string "Calculates codon
adaptation index" returns 0 hits across the whole repo tree and 0 hits in `searchindex.js`.
Correction: replaced the fabricated quotation with an unquoted description of what the function does.

**CA-4 — `utils.rarecodons()` quotation invented.** APPLIED.
Verified: `src/pydna/utils.py:235-236` docstring is `"""docstring."""`; "Identifies rare codons" = 0 hits repo-wide
and 0 in `searchindex.js`.
Correction: same treatment; added that all three `utils` codon functions carry only the placeholder docstring.

**CA-other-1 — "using Python" should be "using python".** APPLIED.
Verified against the extracted PDF title line: "Pydna: a simulation and documentation tool for DNA assembly
strategies using python". One character changed on record line 4.

**CA-other-2 — `homologous_recombination_inversion` does not exist.** APPLIED (partially; audit half-wrong).
Verified by listing every top-level `def` in `assembly2.py`: the real names are
`homologous_recombination_excision_or_inversion` (L2876) and `cre_lox_excision_or_inversion` (L2987);
`homologous_recombination_excision` (L2915) and `cre_lox_excision` (L3040) exist but emit
`DeprecationWarning`. Also verified `crispr_integration` (L3145), `recombinase_integration` (L3088),
`recombinase_excision_or_inversion` (L3050) exist and were missing from the inventory.
Correction: fixed the invented `/inversion` names, added the three missing back-ends, and noted the plain
`*_excision` names are deprecated wrappers.
**The audit is wrong on one half:** it claims the record wrote `cre_lox_.../inversion`. The record wrote
`cre_lox_integration/excision`, which was literally accurate. Only the `homologous_recombination` name was
invented.

**CA-other-3 / FC-A5 — unqualified "all" paths / "all valid orders and orientations".** APPLIED.
Verified in `assembly2.py`: `Assembly.__init__` defaults `use_fragment_order: bool = True`, whose docstring
says "only assemblies that start with the first fragment and end with the last are considered";
`get_linear_assemblies(only_adjacent_edges=False, max_assemblies: int = 50)` at L1451 raises `ValueError`
above 50; cycle enumeration wrapped in `limit_iterator(nx.cycles.simple_cycles(self.G), 10000)`.
Also verified the counterweight: the class docstring states that for **circular** assemblies
"use_fragment_order is ignored", so circular enumeration genuinely is exhaustive over orders and
orientations — and the notebook quotation the record leans on is a *circular* example, so it is unaffected.
Correction: two clauses. First occurrence now states the circular/linear asymmetry and the two caps; the
"all valid orders and orientations" sentence now carries "(circularly, or linearly with
`use_fragment_order=False`)". `combinatorial_motif_place = partial` untouched — the "why not `yes`" reasons
(homology-driven, no placement API, simulation not specification) are unaffected.

**CA-other-4 — PCR feature-label inheritance.** APPLIED.
Verified `src/pydna/amplify.py:364-368`: the selector is
`f.location.start == 0 and f.location.end == len(self.template)` — a feature spanning the whole **template**,
not one spanning the amplicon.
Correction: replaced "when a feature fully spans the amplicon" with "when a feature spans the entire
template", with the predicate quoted.

**CA-other-5 — `gel.gel()` does not use matplotlib.** APPLIED.
Verified `src/pydna/gel.py`: imports are `scipy.interpolate.CubicSpline` (L18), `numpy` (L37),
`PIL.Image`/`PIL.ImageDraw` (L38-39). No matplotlib import or call anywhere in the file. (Matplotlib is in the
`gel` extra in `pyproject.toml`, but it is `Contig.figure_mpl()` that uses it.)
Correction: "via PIL/matplotlib" → "via NumPy and PIL, with SciPy for the migration interpolation", plus one
clause pointing matplotlib at `figure_mpl()`.

**CA-other-6 — "single repo-wide code-search hit" for `hgvs`.** APPLIED.
Verified with a case-insensitive count over the extracted commit tree: 31 occurrences in two tracked files —
`tests/pUC_LAC4_correct_rotation.gb` (1) and `docs/notebooks/gibson_eg.gb` (30). Inspected the contexts: all
are `HGVS` inside protein `/translation=` qualifiers, i.e. incidental amino-acid substrings. The record's
separate claim of 0 `hgvs` hits in `searchindex.js` was independently re-confirmed (count = 0).
Correction: replaced the count and named both files and the actual context. `hgvs_input = no` untouched.

**CA-other-7 — `isorf()` paraphrase presented as a quotation.** APPLIED.
Verified `src/pydna/seqrecord.py:145`: "Detect if sequence is an open reading frame (orf) in the 5'-3'.
direction."
Correction: replaced `"is this an ORF?"` with the real docstring text.

**CA-other-8 — `primer_design()` quotation not verbatim.** APPLIED.
Verified `src/pydna/design.py:83-84`: "This function designs a forward primer and a reverse primer for PCR
amplification of a given template sequence."
Correction: quotation replaced with the verbatim text.

**CA-other-9 / CA-other-10 — `utils.express` misquoted; "surfaced as" wrong for `express`.** APPLIED.
Verified `src/pydna/utils.py:248-252`: docstring is `"""docstring.\n\n    **NOT IMPLEMENTED YET**\n    """`
and the body ends `raise NotImplementedError`. The string "Not yet implemented" is absent from the repo and
from `searchindex.js`. Separately verified `Seq.express()` at `src/pydna/seq.py:172-202` is a fully
implemented method that builds and returns a `PrettyTable`; it does **not** call `utils.express()`.
Correction: three occurrences of the misquotation fixed (Block D entry, summary table, Stated limitations),
and the "surfaced as" sentence now says `Seq.express()` is separately implemented rather than a wrapper.
`codon_optimization = no` untouched — no back-translation or CDS-rewriting function exists either way.

**CA-other-11 / FC-A7 — SPDX headers.** APPLIED.
Verified on `master`: 40 files match `src/pydna/*.py`; all 40 carry an `SPDX-License-Identifier` line; 39 are
`BSD-3-Clause`; `src/pydna/threading_timer_decorator_exit.py:3` is `SPDX-License-Identifier: MIT` with the
full MIT text, and it is live code (`src/pydna/assembly.py:61` imports `exit_after` from it).
Also verified on tag `v5.5.16`: 40 files, **0** carry any SPDX line.
Correction: one sentence replaced with the accurate breakdown, including that the released 5.5.16 sources
carry no SPDX headers. Summary-table basis updated. `license = yes (BSD 3-Clause)` untouched — the package
licence rests on `LICENSE.txt`, `pyproject.toml` and the PyPI classifier, all re-verified, and the single
vendored MIT file is permissive and compatible.

**CA-other-12 — "67 open issues".** APPLIED.
Verified: `api.github.com/repos/pydna-group/pydna` → `open_issues_count = 67`; GitHub issue search →
61 open issues; PR search → 6 open PRs. 61 + 6 = 67, confirming 67 is the aggregate.
Correction: both occurrences (Block E `maintained` entry and the Availability section) now give the split and
name the aggregate field.

**CA-other-13 — lactose scripts are not in `docs/cookbook/`.** APPLIED.
Verified: `docs/cookbook/` on `master` contains YEp24PGK_XK and pGUP1 assets (`YEp24PGK_XK.gb`, `pGUP1.gb`,
`GUP1_locus.gb`, `cookbook.ipynb`, …) but no lactose scripts; a repo-wide search finds no set of nineteen
scripts and no `pYPK0_PDC1tp_KlLAC4_…` file. Verified in the paper text: "The lactose pathway was described by
nineteen short pydna scripts (Additional file 2)".
Correction: scope note appended to the lactose bullet only; the heading and the other two bullets, which the
audit correctly says are accurate, were left alone.

**CA-other-14 — Google Code link classified as dead.** APPLIED.
Verified: `curl -L https://code.google.com/p/pydna/` → final URL `https://code.google.com/archive/p/pydna`,
HTTP 200. The archive shell is JS-driven, so I also fetched the underlying archive record at
`storage.googleapis.com/google-code-archive/v2/code.google.com/pydna/project.json` (HTTP 200) and confirmed it
holds the real pydna project description. Also re-confirmed `pydna-shell.appspot.com` → 404 and the old
GitHub URL → 301 to `pydna-group/pydna`.
Correction: Google Code moved out of the dead-links list and described as obsolete-but-archived.

### Fact-check section A

**FC-A2 — provenance wrongly dated to the v6 development line.** APPLIED.
Verified against the *released artefact*, not just the git tag: downloaded the PyPI
`pydna-5.5.16-py3-none-any.whl` and read inside it — `pydna/dseqrecord.py` contains `def history`,
`def validate_history`, `def normalize_history`; `pydna/opencloning_models.py` contains `class CloningStrategy`,
`from_dseqrecords`, `to_dseqrecords`. So the whole provenance/DAG machinery ships in the current stable
release.
Correction: "since the v6 development line" replaced with "already in the current stable release, PyPI
5.5.16, not only on the v6 development line".

**FC-A3 — "every derived `Dseqrecord` carries a typed `source` object".** APPLIED (same sentence as FC-A2).
Verified: `opencloning_models.py` module docstring says "Not all methods generate sources so far" — a line the
record already quotes verbatim in its own `per_sequence_provenance` entry, so the record contradicted itself.
Verified concretely: `Dseqrecord.reverse_complement` (L984) assigns
`ReverseComplementSource(...)`, but `Dseqrecord.__add__` (L810-830) returns a new `Dseqrecord` with no source
assignment, and `dseqrecord.py:257` explicitly sets `obj.source = None`.
Correction: "every derived `Dseqrecord`" → "a `Dseqrecord` produced by a modelled operation", with pydna's own
caveat and the `__add__` counter-example named. `per_sequence_provenance = yes` untouched.

### Fact-check section B

**FC-B7 — `primer_screen` operational restrictions omitted.** APPLIED (partially).
Verified `src/pydna/primer_screen.py:62-68`: a module-level
`warnings.warn("The primer_screen module is experimental and not yet extensively tested. api may change in
future versions.", category=FutureWarning)` fires on import. Verified `pyproject.toml`:
`pyahocorasick = { version = ">=2.2.0", optional = true }` and `primer_screen = ["pyahocorasick"]`.
Correction: one sentence added to the existing `primer_design` entry, because the record leans on
`primer_screen` as its headline reason for widening that cell and the experimental status materially
qualifies it.
Not applied from this finding: the savable/reusable 3′-suffix automaton, the unique-binding-site filter, and
the IUPAC-expansion blowup — see escalation E-3.

**FC-B8 — assembly-engine bounds omitted.** APPLIED (partially, via CA-other-3).
The record now names `use_fragment_order`, `max_assemblies=50` and the 10,000-path cap. Not applied:
`use_all_fragments`, `only_adjacent_edges`, `circular_only`, one-use-per-fragment — see escalation E-3.

**FC-B19 — codon-diagnostic organism scope.** APPLIED.
Verified `src/pydna/codon.py`: `weights = {"sce": _sce_weights}` (no `eco` key);
`rare_codons = {"sce": [...], "eco": [...]}` — **both populated**; `start`, `stop` and `n_end` each have an
`"eco"` key mapping to `{}`.
Correction: the parenthetical "(yeast `sce` is the shipped organism)" replaced with the exact per-table
breakdown.

---

## REJECTED (5)

**CA-other-15a — recorded `pushed_at` is wrong.** REJECTED (unverifiable in the record's favour or against).
The audit reports the live API now shows `pushed_at = 2026-08-11T14:09:17Z`, which I confirmed. But
`pushed_at` is mutable and the record explicitly stamps its snapshot as 2026-08-11. A read taken on 2026-08-11
before 14:09 UTC would legitimately have returned 2026-08-10T08:46:14Z. I cannot establish that the record was
wrong at retrieval time, so under the standing rule I left the dated snapshot untouched. (The second half of
this finding is escalated as E-1.)

**FC-A7b — "master's `threading_timer_decorator_exit.py` lacks an SPDX line".** REJECTED — the audit is wrong.
`src/pydna/threading_timer_decorator_exit.py:3` on `master` reads `# SPDX-License-Identifier: MIT`. It has an
SPDX line; it is simply not BSD-3-Clause. The correct version of this point was applied via CA-other-11.

**FC-B1 — random-generator positive API not recorded.** REJECTED — no defect shown.
The record states these functions "take only `length` (and optional `maxlength`) — unconstrained". Verified
`src/pydna/utils.py:444-544`: `randomRNA(length, maxlength=None)`, `randomDNA(...)`, `randomORF(...)`,
`randomprot(...)`. The record's claim is exactly correct. What the audit wants added (triangular length
selection, ATG-start/random-stop ORF construction) is extra detail, not a correction, and the record's use of
these functions is solely to show they are not constrained barcode generators — which stands.

**FC-B4 — bulk/mixed-format ingestion not explained.** REJECTED — already in the record.
The `lazy_evaluation` entry already states that `parsers.parse()` is annotated `-> list[Dseqrecord |
SeqRecord]` and documented as "a greedy function, use carefully", and the "Round-trip with lab file formats"
item already lists GenBank/EMBL/FASTA/SnapGene parsing. The audit's own point is about which *section* it
appears in, which is placement, not fact.

**CA-other-2b — record used the name `cre_lox_..._inversion`.** REJECTED — the audit misread the record.
The record's inventory reads `cre_lox_integration/excision`; both of those functions exist
(`cre_lox_integration` L2927, `cre_lox_excision` L3040). Only the `homologous_recombination` entry contained
an invented name, and that half was applied.

---

## ESCALATED (5)

**E-1 — "Last activity: 10 days before this survey".** (citation audit other-15, second half)
Both dates in the record are correct: master head commit 2026-08-01, repo push 2026-08-10, survey 2026-08-11.
2026-08-01 → 2026-08-11 is exactly 10 days, so the figure is right *if* "activity" means the last commit; it
is 1 day if a push counts. Deciding which one the survey means is an editorial definition, and it sits next to
a `maintained = yes` value it cannot change in either direction. Left untouched.
**Question:** does "Last activity" in this survey mean last commit to `master`, or last repository push?
**Options:** (a) keep "10 days" and add "(last commit)"; (b) change to "1 day (last push)"; (c) give both.

**E-2 — "The true and only limit: there is no mutagenic-primer designer" (`primer_design`).** (fact-check A6)
The audit's underlying facts check out: `primer_design` rejects two supplied primers (`design.py:74-247`);
`assembly_fragments` requires at least every second non-link fragment to be an `Amplicon`; `primer_screen`
selects from an existing stock, needs the optional `pyahocorasick` extra, and self-declares experimental. But
whether those count as *material* limits, such that "only" is wrong, is a judgment about emphasis, and it
attaches directly to a `yes` cell. I applied only the one verifiable factual omission (the experimental
warning, FC-B7) and left the "only" sentence alone.
**Note the resulting local roughness:** the entry now says `primer_screen` is experimental two paragraphs
above the claim that the mutagenic-primer gap is the "only" limit. I did not smooth this over.
**Question:** should "The true and only limit" be softened to "The main limit", or does the sentence mean
"only limit *on primer design as such*", in which case it stands?
**Options:** (a) soften to "The main limit"; (b) keep "only" and scope it explicitly to design capability;
(c) keep as is.

**E-3 — the remaining 14 fact-check section-B omissions.** (FC-B2, B3, B5, B6, B9, B10, B11, B12, B13, B14,
B15, B16, B17, B18)
I verified the substance of the ones I spot-checked — `ProteinSeq` (`seq.py:265`) with `molecular_weight()`,
`pI()`, `instability_index()`, and `ProteinSeqRecord` (`seqrecord.py:701`) exist and are absent from the
record (B2); `utils.eq` (L545) and `utils.deduplicate` (L990) exist and are absent (B3);
`golden_gate_assembly` (L2508) literally `return restriction_ligation_assembly(...)` and its docstring says
"This is the same as restriction ligation assembly, but with a different name", `in_fusion_assembly` says
"This is the same as Gibson assembly, but with a different name", and README line 25 reads "Golden gate
assembly (in progress)" (B9). So these are true omissions, not audit errors.
I did **not** add them. Each is an API detail with no home in any of the 24 capability keys; adding fourteen
clauses would create a large block of fresh unaudited text and would itself shift the record's balance, which
is the authors' call. Nothing here changes a value.
**Question:** which of these fourteen, if any, should be folded into the record?
**Options:** (a) none — the record is a 24-key capability memo, not an API reference; (b) add only the ones
that qualify a claim the record already makes (B9 Golden Gate / In-Fusion are aliases; B13 CRISPR has no guide
design or off-target search; B16 the SnapGene unsupported-operation list); (c) add all fourteen.

**E-4 — fact-check section C (balance).** Explicitly a balance and emphasis judgment: whether the record
documents library-design absences more deeply than native collection workflows and operational limits.
Reserved for the authors; no edit made.
**Question:** should the record be rebalanced toward native collection workflows and operational caveats?
**Options:** (a) leave as is — the record's job is the 24 keys; (b) commission a balance pass after the value
reconciliation; (c) act on E-3 option (b) and treat that as sufficient.

**E-5 — `library_as_object` evidence is now weaker than it was.** (consequence of CA-1 / FC-A1)
Not a finding either audit raised as such, but a direct consequence of an applied edit, and the brief says to
escalate when a correction changes what a value rests on. The cell's headline argument used to be the clean
absolute "no collection type appears among the 820 documented objects". That absolute is false and is gone.
The value still holds on the substantive ground (no *variant-library* abstraction; assemblies return plain
`list[Dseqrecord]`; building many variants means writing your own loop), and I left `no` in place — but a
referee can now point at `RecombinaseCollection`, `OpenDNACollectionsSource` and `CloningStrategy` and say the
survey's own text concedes collection types exist.
**Question:** is `library_as_object = no` still the right value now that the absolute inventory claim is
withdrawn, and does the cell need a stronger substantive lead sentence?
**Options:** (a) keep `no`, no further change; (b) keep `no` but lead with the variant-library argument rather
than the object-index count; (c) reopen the value in the ratings-reconciliation pass.

---

## Tensions left unsmoothed (deliberate)

Per the brief, nearby sentences that now read slightly oddly were left alone rather than rewritten:

1. `primer_design` — the new experimental-status clause sits above the untouched "The true and only limit"
   sentence. See E-2.
2. `dag_chaining` — the "Stated limitations" bullet still says "`CloningStrategy` serialises", which
   understates the round-trip now described in Block A. It is not false (it does serialise, and the operative
   claim "cannot be re-executed on new inputs" is verified correct), so it was not touched.
3. `combinatorial_motif_place` — the entry now states the linear-enumeration default constraint and then, two
   sentences later, retains "**all valid orders and orientations**" with a parenthetical qualifier. The
   parenthetical is accurate; the emphasis is the original author's.
4. Record header — the "Record status" paragraph still asserts the memo was fully re-verified from source.
   That framing is now partly contradicted by this pass. Left untouched; it is an authorial statement.

---

## Artifacts

| Path | Action |
|---|---|
| `revision/tool_survey/final/pydna.md` | modified in place — 21 surgical edits, no capability value changed |
| `revision/tool_survey/fixes/pydna.md` | created — this change log |
| `revision/tool_survey/fixes/` | created (directory) |

## Pass 2 — policy application

**Date:** 2026-08-14

**Verification basis and counts.** The paper was read directly from
`papers/Pereira2015_pydna.pdf` with PyMuPDF. Repository claims were checked against the official
`pydna-group/pydna` stable tag `v5.5.16` and current `master` head
`4e02f81dfd39629c68d2a93dfbb434bfe7d56120`; release facts were refreshed from PyPI and repository facts
from GitHub's official API. The 18 open findings are counted atomically (E1, E2, fourteen E3 omissions, E4
and E5): **7 applied · 11 declined-by-policy · 0 whole findings rejected · 0 still escalated**. One false
subclaim inside applied finding B17 was independently rejected and is counted separately below. Policy C was
also applied; Policy D required no edit. No capability value changed.

**Version/source check.** PyPI still reports **5.5.16** (uploaded 2026-06-09) as the current stable release,
and `master` still points to the assessed 2026-08-01 commit. GitHub now reports
`pushed_at=2026-08-11T14:09:17Z`. The header versions therefore still govern and no Policy-D current-release
parenthetical is warranted.

| Open item | Outcome | Edit and primary-source verification |
|---|---|---|
| **E1 / citation-audit other-15 — ambiguous “Last activity”** | **applied** | Replaced the ambiguous label with two definitions: “Last `master` commit” (2026-08-01, ten days before the 2026-08-11 survey) and “Last repository push” (`2026-08-11T14:09:17Z`, the survey date). The maintenance entry and availability facts now use the freshly fetched API value. GitHub's commit API confirms the head SHA/date; the repository API confirms `pushed_at`. `maintained = yes` remains supported. |
| **E2 / FC-A6 — “true and only” primer-design limit** | **declined-by-Policy-A** | The underlying facts all verify in `design.py`, `primer_screen.py` and `pyproject.toml`: two supplied primers are rejected by `primer_design`; `assembly_fragments` enforces the Amplicon adjacency rule; primer screening selects from stock, is optional/experimental and filters pair routines to unique binding sites. The authors' ruling classifies the requested “only”→“main” rewrite as balance/emphasis, so the sentence was not edited. |
| **E3 / FC-B2 — protein objects and property analysis** | **applied (Table-1 Purpose/Output)** | Added one clause to existing additional-capability item 16: translation returns `ProteinSeq`/`ProteinSeqRecord`, and `ProteinSeq` reports molecular weight, pI and instability index. Verified in stable `seq.py:265–340` and `seqrecord.py:191–194,704+`. These numeric analysis outputs materially qualify the grouped Table-1 description rather than a Table-2 score. |
| **E3 / FC-B3 — equality and deduplication** | **applied (Table-1 Key features)** | Added a brief clause to existing round-trip/collection item 5: `utils.eq()` handles case, reverse complement and optional circular rotation, while `utils.deduplicate()` removes repeats in order. Verified in stable `utils.py:545–614,990–1012`. This bears on general sequence manipulation/collection handling without creating a library abstraction. |
| **E3 / FC-B5 — PCR product/annealing semantics** | **declined-by-Policy-B** | Verified in `amplify.py`: `_annealing_positions` requires a perfect 3′ suffix of at least `limit`; `Anneal.products` crosses forward/reverse sites; `pcr()` raises unless exactly one product exists. These details do not change Table 1's broad “simulate cloning” purpose or transformed-sequence output, and pydna has no Table-2 column. |
| **E3 / FC-B6 — ordered assembly-primer rules** | **declined-by-Policy-B** | Verified in `design.py:249–790`: `overlap`, `maxlink`, circular handling, absorption of short linkers into primer tails, and the at-least-every-second-Amplicon rule are real. They refine one cloning workflow but do not alter a Table-1 cell. |
| **E3 / FC-B9 — Golden Gate/In-Fusion wrappers** | **declined-by-Policy-B** | Verified: stable `golden_gate_assembly` directly calls `restriction_ligation_assembly`; `in_fusion_assembly` uses the common Gibson-style product path with different provenance; README still labels Golden Gate “in progress.” These implementation qualifications do not alter Table 1 Purpose, Key features, Output or Availability. |
| **E3 / FC-B10 — ligation mode switches** | **declined-by-Policy-B** | Verified in `assembly2.py:2365–2601`: partial-digest filtering, `allow_blunt`, `allow_partial_overlap` and `circular_only` behave as reported. They change which reaction products are returned but not the high-level Table-1 cells. |
| **E3 / FC-B11 — Gateway modes/filtering** | **declined-by-Policy-B** | Verified in `assembly2.py:2615–2727` and `gateway.py`: BP/LR, conservative/greedy, circular-only and multisite filtering exist; the greedy consensus may produce false positives and multisite filtering removes one-site intermediates. No Table-1 cell changes. |
| **E3 / FC-B12 — homologous/site-specific recombination modes** | **declined-by-Policy-B** | Verified in `assembly2.py:2728–3140` and `recombinase.py`: multiple inserts, integration/excision/inversion/general assembly, and the uppercase-flank + contiguous lowercase shared-core format are implemented. This is back-end detail inside the already accurate cloning-simulation purpose. |
| **E3 / FC-B13 — CRISPR scope/limits** | **declined-by-Policy-B** | Verified in `crispr.py` and `assembly2.crispr_integration`: `cas9` is a 20-nt protospacer plus `.GG`; `protospacer()` extracts from a supplied construct; integration requires caller-supplied guides; `search` warns by contract that it ignores `Dseq.circular` unless `linear=False`. The module/source exposes no guide designer or genome-wide off-target scorer. These limits do not alter a Table-1 cell. |
| **E3 / FC-B14 — two-oligo hybridization enumeration** | **declined-by-Policy-B** | Verified in `oligonucleotide_hybridization.py`: exactly two `Primer`s and a minimum annealing length yield every allowed overhang placement, and a mismatch-capable alignment raises. This returns transformed sequences already covered by Table 1. |
| **E3 / FC-B15 — additional physical-end transformations** | **declined-by-Policy-B** | Verified in `dseq.py`: `user`, `melt`, `melt_ss_dna`, `shed_ss_dna`, and nucleotide-selective `T4` behavior exist. These are lower-level transformations within the already accurate grouped Key-features cell. |
| **E3 / FC-B16 — SnapGene-history support boundaries** | **applied (Table-1 Key features/file-format support)** | Added one limitation clause to existing provenance-import item 4. `snapgene_history_parser.py` explicitly rejects GC/TA/TOPO cloning and destroyed restriction fragments, lacks blunting-dependent removal, and stops on listed manual insertion/replacement/reverse-translation/codon/site edits. These boundaries materially qualify the grouped row's file-format support. |
| **E3 / FC-B17 — Sanger-trace mapping** | **applied with one subclaim rejected (Table-1 Key features/Output)** | Added the verified interface to item 5: `map_trace_files()` reads ABI traces, optionally filters a target by either target or reverse-complement containment, maps **direct** common substrings, adds `trace` features and returns matched filenames. Stable `dseqrecord.py:730–805` verifies that behavior. The audit's stronger phrase “finds common substrings in either orientation” is **rejected as contradicted**: feature placement calls `common_sub_strings(str(self.seq), str(trace.seq), limit)` only once and never reverse-complements either argument; `common_sub_strings.py:315–342` performs literal same-orientation matching. |
| **E3 / FC-B18 — ladders and `FakeSeq`** | **applied (Table-1 Output)** | Extended the existing virtual-gel item with one clause: built-in ladder lists use `FakeSeq` length/mass/mobility records without full nucleotide sequences. Verified in stable `ladders.py`, `fakeseq.py` and `gel.py`. This directly qualifies pydna's non-sequence output. |
| **E4 / fact-check section C — balance** | **declined-by-Policy-A** | No balance or proportional-emphasis edit was made. Factual corrections and Table-1-qualified omissions were handled independently above. |
| **E5 / CA-1 consequence — `library_as_object` value basis** | **applied (evidence only)** | Kept the legacy capability value `no` and replaced the weakened absolute inventory lead with the verified distinction: `CloningStrategy.from_dseqrecords()` collects one strategy's sequences/sources/primers and `to_dseqrecords()` reconstructs and returns terminal records, so it is a real strategy/provenance collection but not a designed variant-library abstraction. `RecombinaseCollection` and `OpenDNACollectionsSource` are named separately. Stable `opencloning_models.py:1225–1404`, `recombinase.py` and the docs index verify the distinction. The adopted reworded row remains `yes` via `Assembly`, exactly as the existing supersession note says; pydna is absent from locked Table 2. No value was changed. |

### Policy C — primary-source boundary, current facts and definitions

**Outcome: applied.** Removed the four sibling prior-analysis/review rows and the related-project/live-site
rows from the source table. Removed cross-tool score arguments, competitor/referee anecdotes, the
survey-wide provenance superlative and claims about behavior in a sibling backend. OpenCloning/teemi are
mentioned only where pydna's own source or official gallery documents an export or external example; their
repositories are not evidence. The stale rate-limit hedge was replaced with current GitHub API facts.

Definitions are now explicit and consistent: *variant library* means a designed set of sequence variants;
`CloningStrategy` is called a strategy/provenance collection; reaction-product lists are not variant
libraries; “last `master` commit” and repository `pushed_at` are separate activity measures; the interface
value means the shipped Python API, not a related web interface. Unsupported live-site anecdotes were
removed rather than marked “unverified.”

### Policy D — version drift

**Outcome: no edit required.** The header's PyPI 5.5.16 remains the current stable release and its assessed
`master` commit remains the branch head. No materially different current release exists, so the single
allowed parenthetical was not added and no scoping note was renumbered or restructured.

### Rejections, escalations and value checks

- **Rejected subclaim:** FC-B17's “common substrings in either orientation,” contradicted by the direct-only
  placement call described above. The remainder of B17 was verified and applied accurately.
- **Remaining escalations:** none. E1–E5 are all resolved by the rulings and primary-source checks.
- **Capability values changed:** **0**. Corrected evidence continues to support every retained value under
  its stated definition. The legacy/adopted `library_as_object` distinction remains explicit rather than
  silently changing either value.
- **Locked-row score implications:** none. Pydna has no locked Table-2 column, and no Table-2 value was used
  as a shipping basis.

### Row-substitution report — no action

The locked file's two pre-existing risks remain only risks: **`Codon / amino-acid substitutions`** may be
near-uniform and **`Recombination and chimeras`** may be one-tool-only, but this pydna-only pass provides no
eight-tool evidence to confirm either distribution.

**Conditional candidate if a general-purpose-tool column is ever restored:** replace the flagged
`Codon / amino-acid substitutions` row with **`Executable construction provenance`** — structured
per-output construction history that can be serialized and replay-validated; `partial` if history is only
free text or cannot be validated. Pydna has unusually strong primary evidence: typed `source` trees,
`Dseqrecord.history()`, `validate_history()`, `normalize_history()`, `CloningStrategy` JSON round-trip and
SnapGene-history import (`opencloning_models.py`, `dseqrecord.py`, `snapgene_history_parser.py`). Because
pydna is not a current Table-2 column and the other eight cells have not been verified, this is not a valid
locked-table substitution today and was not applied.

### Neighbouring tension logged

1. `primer_design` still places the experimental/stock-primer limitations beside the untouched “true and
   only limit” sentence; Policy A expressly requires that result.
2. The legacy `library_as_object = no` heading remains beside its existing supersession note stating that
   the adopted reworded row is `yes`; the evidence now defines the two meanings cleanly.
3. Table 1's grouped Output cell still says only “Transformed sequences,” while verified pydna outputs also
   include protein-property values and virtual-gel images. This pass was authorized to repair the pydna
   record, not the locked table; the qualifying evidence was added to existing entries for the table authors.
4. The historical `## ESCALATED (5)` section above still records the Pass-1 state. This Pass-2 section
   supersedes those outcomes without rewriting the audit trail.
