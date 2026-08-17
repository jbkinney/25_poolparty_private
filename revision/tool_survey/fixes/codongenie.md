# CodonGenie — repair pass change log

**Record repaired:** `final/codongenie.md` (edited in place)
**Audits processed:** `citation_audit/codongenie.md` (24 findings), `factcheck/codongenie.md` (A1–A8, B1–B7, C)
**Repair date:** 2026-08-14
**Outcome counts:** 28 findings applied (30 edits, two findings occurring in two places each) · 3 rejected · 3 escalated

**No capability value was changed.** Every edit touches evidence or prose only.

## How verification was done

- Paper: `papers/Swainston2017_codongenie.pdf`, all 10 pages re-extracted with PyMuPDF (PDF page *n* = journal page *n*).
- Repo: fresh `git clone` of `neilswainston/CodonGenie`, HEAD `c26439f` (= the commit both audits froze on); all source files read locally, `git log` for authorship.
- GitHub REST: repo metadata, releases, tags, pulls, RidgeBio fork metadata and commits.
- PyPI: `pypi.org/pypi/CodonGenie/json`, plus the 1.7 wheel and sdist downloaded and unpacked (not installed).
- Live service: `https://codongenie.appspot.com` REST calls; `http://codon.synbiochem.co.uk` HTTP + `openssl s_client`; HTTP status on all 11 CDN assets.
- Nothing was installed and no container was built, per the survey constraint.

---

## APPLIED

| # | Finding | Correction made | Verified how |
|---|---|---|---|
| cit-1 | §11.3 claims `"Written in Python/Flask, ∼700 lines"` is "confirmed verbatim, p. 6" | Replaced the trailing clause: the phrase is the prior note's own wording; the two underlying passages are quoted with their true pages (p. 4 and p. 6) and the record now states the combined phrasing is not in the paper | PDF p. 4 "written in Python (using the Flask framework) and HTML/Javascript"; p. 6 "the entire application consists of ∼700 lines of code"; `prior_analyses/Swainston2017_codongenie_analysis.md` line 20 carries the composite phrase |
| cit-2 | Proline quotation silently corrects an authors' typo | Quote restored to the published text with `[sic]`: "The codon returned is **the therefore** [sic] the most frequent codon…" | PDF p. 7, verbatim |
| cit-3 | The "live JSON" block for taxid 37762 carries taxid 83333's TTT probability | `"probability":0.567` → `0.636` | Live `?aminoAcids=FILMV&organism=37762` → TTT `0.6358653093187157`; `&organism=83333` → `0.5674820143884892` |
| cit-4 / fc-A2a | "Its entire public surface is two operations" | → "Its public surface is two codon operations, plus organism listing and search" | `main.py` routes `/organisms/` and `/organisms/<term>`; `client.py` exposes 4 methods |
| cit-5 / fc-A1 | Design described as enumerating *every* ambiguous codon | Clause replaced with what the algorithm actually does, plus one sentence flagging the paper/software disagreement and the live counterexample | `codon_selector.optimise_codons`/`__analyse`/`_optimise_pos_3`; live `?aminoAcids=P&organism=37762` → only `CCG, CCT, CCA, CCC`; `?aminoAcids=M` → only `ATG`; paper p. 2 "optimized for computational efficiency, compared to a brute-force examination" |
| cit-6 / fc-A2b | "Every response is `Response(json.dumps(...))`" (§2) | → "Every response *except* `/`, which serves `static/index.html`, is …" | `main.py:40-43` `return app.send_static_file('index.html')` |
| cit-6 | Same claim restated in `vcf_vep_output` evidence | → "every data route (`/codons`, `/organisms/`, `/organisms/<term>`) returns … while `/` serves the static UI" | as above; value `no` untouched — no export/download control exists in `result.html`/`result.directive.js` |
| cit-7 | "`codon_genie/` contains five Python modules" — there are six `.py` files | §2 anchor statement → "**five** implementation modules (plus a docstring-only `__init__.py`)" | `git ls-files`; `codon_genie/__init__.py` is a 7-line docstring, no code. *Applied at the anchor only — see "Tensions left in place"* |
| cit-8 | Git-tree source URL is not a live endpoint | `github.com/neilswainston/CodonGenie/git/trees/master?recursive=1` → `api.github.com/repos/neilswainston/CodonGenie/git/trees/master?recursive=1` | non-API URL → 404; API URL → 200 |
| cit-9 / fc-A3 | "no random sampling" (`mixed_mutagenesis_one_pool`) contradicted by shipped code | Scoped: "no random sampling **in the exposed design path** (the unexposed `codon_utils.get_random_codon`/`mutate` do sample synonymous codons stochastically)"; opening clause scoped to "in the exposed tool" | `codon_utils.py:119-125` (`mutate`), `:135-157` (`get_random_codon`) — both present in the 1.7 wheel |
| cit-9 / fc-A4 | "The tool has no sequence context at all" (`combinatorial_motif_place`) | → "The **exposed** tool has no sequence context at all" | `codon_utils.get_codon_optim_seq` and `seq_utils.find_invalid` operate on whole sequences; neither is reachable from `main.py` |
| cit-9 / fc-A4 | "There are no output *sequences*" (`per_sequence_provenance`) | → "The **exposed tool emits** no output *sequences*" | `get_codon_optim_seq` returns a DNA string; `example/align.py:seq_from_alignment` returns a concatenated IUPAC string |
| fc-A5 | "it never enumerates variants" (`assay_dms`) contradicts the record's own §4.2 | → "it enumerates only the DNA variants of one codon, emits no sequences from the design API" | `codon_selector.__analyse_ambig_codon` materialises `ambiguous_codon_expansion`; live DBK → 18 variants |
| cit-10 | Analyse regex attributed to `static/index.html` | Attribution corrected to `static/codonGenie/codonGenie.ctrl.js:5` (bound in `index.html` via `data-ng-pattern`); source cell updated | `grep -rn acgtmrwsykvhdbn static/` → only `codonGenie.ctrl.js:5`; `index.html:169` has `data-ng-pattern="ctrl.codon_pattern"` |
| cit-11 | "five CDN assets" | → "11 external CDN assets (across five hosts)" in `interface` and §9 | `index.html` head lists 11 protocol-relative asset URLs across 5 hostnames; all 11 returned HTTP 200 on 2026-08-14 |
| cit-12 / fc-B1 | §4.8 "Alternative genetic codes" overstates support | Added one sentence: only the codon→amino-acid lookup follows `table_id`; candidate generation still uses the hard-coded standard-code `CODONS` map | `codon_selector.py:62-68` vs the module-level `CODONS` map `:16-36` |
| cit-13 / fc-A6 | §4.5 claims the output states "expected variant-composition skew within the designed pool" | Reframed as **encoding redundancy within the designed pool**, with the clarification that `probability` is host codon usage, not a synthesis-mixture fraction | `codon_utils.__get_codon_usage`/`_scale` (Kazusa frequencies, normalised per amino acid); result schema in `codon_selector.py:119-130` has no mixture field |
| cit-14 | "DTS 0.778, DTK 0.294 — the paper's 0.79 / 0.29 to rounding" | → "the same reversal as the paper, with DTK matching its 0.29 but DTS rounding to 0.78 against the paper's 0.79" | Live `?aminoAcids=FILMV&organism=1902` → DTS `0.7782896589104741`, DTK `0.29425093243948847`; paper p. 4 states 0.79 / 0.29 |
| cit-15 | RidgeBio timestamp labelled "last push" | → "last commit **2025-05-17T04:32:13Z**, `pushed_at` 04:32:35Z" | GitHub API: commit date 04:32:13Z, `pushed_at` 04:32:35Z |
| cit-17 | §9 attributes the wrong subject to the last 2022-06-09 commit | Subjects corrected: last is "Minor auto-reformatting." (`6ef2c3e`, 16:20:09), preceded same day by "Update to remove dependency on remote species.table flatfile." (`e38cc8e`, 13:57:08) | `git log --date=iso` |
| cit-18 / fc-A8 | "a self-hosted instance **silently** loses all codon-usage scoring" | → "not silently: `urlopen` and the `<PRE>` parse are unguarded, so failures propagate to the Flask exception handler and surface as an HTTP 500" | `codon_utils.py:159-188` (no try/except); `main.py:75-82` errorhandler → 500 |
| cit-19 | "The reviewer marked 22 of 23 verdicts 'supported'" | → "would change exactly one of the 24 verdicts (`maintained`) and let the other 23 stand" | Counted 24 capability rows in both `final/` and `extractions/codongenie.md`; `reviews/codongenie.md` says "found **one capability verdict I would change**" and never gives a 22/23 tally |
| cit-20 | clone+Docker called a "working reproduction path" although nothing was executed | "working" → "documented" in §9; `interface` and the §9 verdict now say the container was not built/run here | Record's own §10 ("Nothing in this record was verified by installing or executing the tool"); no build artefact exists and none was produced in this pass |
| cit-23 | Reconciliation table points to the wrong sections | Row 5 `(§3)` → `(§6)`; row 7 `(§6)` → `(§8)` | The record's own headings: §6 = documented examples, §8 = limitations |
| cit-24b | "only Google's `*.appspot.com` wildcard is served" | → the served certificate is Google's, CN `*.appspot-preview.com`, with `*.appspot.com` among many SANs and nothing matching this host | `openssl s_client -connect codon.synbiochem.co.uk:443`; curl exit 60 |
| fc-A7 | 2022-12-12 commit called "an automated `TrellixVulnTeam` bot PR, not author maintenance"; "last human commit 2022-06-09" | `maintained` evidence, §5 C1 and §9 now state that the 2022-12-12 **merge commit is authored by Neil Swainston** (the bot authored only the patch commit), and use "last commit carrying the author's own code changes" for 2022-06-09 | `git log`: `c26439f` author `Neil Swainston <neil@epochbiodesign.com>`, committer GitHub, message "Merge pull request #1 … Update appears to be safe to merge."; `a83d848` author `TrellixVulnTeam` |
| fc-B2 | `get_codon_optim_seq` constraint behaviour not disclosed | Added a parenthetical in `synthesis_constraints`: it re-draws synonymous codons stochastically and, with `tolerant=True`, permits otherwise invalid patterns once `max_attempts` is exhausted, so the constraint is not a hard guarantee | `codon_utils.py:62-105` |
| fc-B3 / fc-B4 | §4.9 name-drops `get_cai` and `mutate` | `get_cai` relabelled "CAI-style scoring" with the arithmetic-vs-geometric-mean point and the silent skipping of unrecognised triplets; `mutate` annotated "one sequence per call, with no check that `dna_seq` encodes `protein_seq`" | `codon_utils.py:107-117` (`sum(w_vals)/len(w_vals)`), `:119-125` |
| fc-B5 | Client library described without its configurability | Added ", and a `url=` constructor argument that can retarget it at a self-hosted instance (default: the live appspot deployment)" | `client.py:15-16` |
| fc-B6 | Flagship alignment example's selection semantics omitted | Added to §6.1: it drops gap characters, takes only the first-ranked codon at each column, and its organism lookup returns the first name-search hit | `example/align.py:13-42` |

---

## REJECTED

| # | Finding | Why rejected | Evidence |
|---|---|---|---|
| cit-16 | Claims the taxid-562 failure payload is `{"message": "KeyError: 'ATG'"}`, not the record's `'ATC'` | Did not reproduce. The record is right | Five fresh live requests to `?aminoAcids=FILMV&organism=562` on 2026-08-14 all returned HTTP 500 with `{"message": "KeyError: 'ATC'"}`. The independent reviewer (`reviews/codongenie.md` line 142) also recorded `'ATC'`. The audit's own caveat is that the key can vary with iteration order, so its `'ATG'` observation does not falsify the record |
| cit-21 | "undeployed" / "the fork is not deployed" is an unevidenced web-wide negative | Could not confirm either direction, and any softened wording would itself be unverified. Applying it would swap one unverified claim for another | I can confirm only the checkable parts, which the record already states correctly: upstream has exactly one PR ever (`TrellixVulnTeam`, closed) and 0 open issues, so there is no RidgeBio PR/issue; PyPI's `CodonGenie` project has 8 releases, all ≤2022-03-29, i.e. all predating the 2025 fork. No search can establish that no private deployment exists |
| cit-24a | "HTTP-redirects here (200)" conflates first-hop and final status | The record's wording is not false: it says the old URL redirects, and reports 200 for the destination | `curl` no-follow → 301 to `https://github.com/neilswainston/CodonGenie`; follow → 200 |

---

## ESCALATED — authors' decision required

### E1. Cross-tool claims in a single-tool record (cit-22)

The record asserts, with no source cited in it, that PoolParty has no analogous audit for user-supplied degenerate codons (§4.3), that PoolParty would be `no` on a proposed "degenerate/ambiguous codon design" row because it "enumerates explicit sequences" (§4 ROW-LIST SUGGESTION), that the redundancy output is "directly comparable to PoolParty's library-composition reporting" (§4.5), and that MPRADesignTools is 5 yr 8 mo dormant and MPRAnator ~10.6 yr (§5 C1).

What I could check: the dormancy figures **agree** with the survey's own sibling records (`final/mpradesign.md`: last commit 2020-12-17, "5 years 8 months"; `final/mpranator.md`: single commit 2015-12-27, "~10.6 years"). The PoolParty claims are the problem — `final/poolparty.md` scores `degenerate_iupac_codons` = **yes** (`base_ops/from_iupac.py` enumerates or samples any IUPAC string), with the caveat that PoolParty *expands* IUPAC rather than *emitting* a degenerate string. So "PoolParty is `no`" on a degenerate-codon row is defensible only if the proposed row is defined as *emitting* a degenerate codon, not accepting one.

**Question:** should this record cite the sibling records for its cross-tool claims, and how should the proposed new row be defined so that PoolParty's existing `degenerate_iupac_codons = yes` and a CodonGenie-favouring "PoolParty = no" do not contradict each other?
**Options:** (a) add in-record citations to `final/poolparty.md` / `final/mpradesign.md` / `final/mpranator.md` and re-word §4.3/§4.5/ROW-LIST to "emits a degenerate codon" framing; (b) drop the cross-tool sentences from this record and handle all comparisons in one place; (c) leave as is and reconcile during matrix assembly.
I made **no edit** to any of these sentences, including the "Directly comparable to PoolParty's library-composition reporting" sentence that now sits beside my corrected §4.5 text.

### E2. Missing limitation: no non-uniform nucleotide-mixture ratios (fc-B7)

I verified the fact: the entire input vocabulary is an amino-acid set, or one 3-character IUPAC codon, plus a taxonomy id (`main.py:61-72`, `static/index.html`), and nothing anywhere accepts mixture fractions, so a designed codon implies an equimolar base mixture and the tool can specify no other physical pool composition. It also computes no multi-position library arithmetic.

I did **not** add it. Whether this omission is *material* enough to add another bullet to §8 is exactly the kind of judgement reserved for the authors — and the same fact-check's section C already calls the record "somewhat limitation-heavy", so adding it trades one audit finding against another.
**Question:** add this to §8 "Unstated but material", fold it into §2's scope description, or leave it out?

### E3. Balance (fc-C)

The fact-check's section C is entirely about emphasis: that limitations are repeated across five sections while package/example capabilities are understated in the headline assessments. No factual claim is disputed. Per the repair rules this is an authors' call and nothing was changed for it. Note that the section-C complaint about "absolute denials contradicting the later acknowledgement" **is** now addressed as a side effect of the applied fc-A3/A4/A5 scoping edits.

---

## Tensions left in place (deliberately not rewritten)

1. **§0 reconciliation row 6** still summarises the added §4.5 capability as "per-amino-acid redundancy counts as pool-composition skew", which is the framing I corrected inside §4.5. Left alone: it is a historical record of what the reviewer recommended (`reviews/codongenie.md` §7 row 9 uses the same phrase), and the corrected detail is one click away.
2. **"five modules"** survives at §0 row 3, `barcode_generation`, §8 and §10; only the §2 anchor now says "five implementation modules (plus a docstring-only `__init__.py`)". Each of those four is either a complete enumerated list or a back-reference to the corrected anchor, and `codon_genie/__init__.py` holds no code, so `barcode_generation = no` is unaffected — rewriting all five sites would create more unaudited text than the imprecision warrants.
3. **`maintained` stays `partial`.** The fc-A7 correction moves the last human-attributed commit from 2022-06-09 to the author's 2022-12-12 merge — from ~4 yr 2 mo to ~3 yr 8 mo. Both fall in convention C1's 2–5-year `partial` band, so the value does not move; the same wording change was made in all four places the fact appears (`maintained` evidence, §5 C1, §9 table, §9 verdict) to avoid the record contradicting itself.
4. **§8's "`__get_codon_usage` silently builds an empty table"** (the taxid-562 bullet) was outside both audits' findings and is untouched, though it now sits next to a bullet that says the Kazusa failure mode is *not* silent. The two describe different failures (missing table vs unreachable service).
5. **`codon_optimization` cites "p. 4 Table 1"**: the scores discussed are on p. 4, but Table 1 itself is typeset on p. 6. Neither audit raised it; not touched.

## Pass 2 — policy application

**Policy date:** 2026-08-14  
**Open-item outcome counts:** 3 considered · 2 applied (`E1`, `E2`) · 1 declined (`E3`) · 0 escalated  
**Matrix impact:** 0 capability values changed · 0 value escalations · 0 Table 2 edits · 0 drift parentheticals

### Open items

| Item | Outcome | Edit | Primary-source verification |
|---|---|---|---|
| **E1 — cross-tool claims / provenance (`cit-22`, with `cit-21` follow-through)** | **APPLIED** | Chose the drop-comparisons outcome. Removed the PoolParty Analyse, library-composition, row-verdict and reproduction comparisons; removed the MPRADesignTools/MPRAnator dormancy comparison; removed all RidgeBio/fork/deployment claims, including the previously unverifiable “undeployed” statement. Replaced the fork section's evidentiary body in place with an upstream provenance boundary; its existing heading was retained under the surgical-edit rule. Internal extraction/review files are now labelled process records, not factual evidence. | PDF pp. 6–7 and upstream `c26439f`: `CodonSelector.optimise_codons` accepts one amino-acid set, and the paper/repository alignment workflow concatenates per-position calls in `example/align.py`. `git ls-remote` and GitHub API re-check on 2026-08-14 show upstream `master` still at `c26439f`; PyPI still reports 1.7. These official sources independently support the affected scope prose. |
| **E2 — no non-uniform nucleotide-mixture ratios (`fc-B7`)** | **APPLIED, TABLE-1-RELEVANT CLAUSE ONLY** | Added one brief qualifier to existing §4 item 5: the interface “accepts and reports no non-uniform nucleotide-mixture fractions.” Added no limitation bullet or new section. The multi-position-arithmetic half was not repeated because the existing single-codon Purpose/Output scope already states it and Policy B excludes non-Table-1 accretion. | Current upstream `main.py` accepts only `aminoAcids` or `codon`, plus `organism`; `static/index.html` exposes only those inputs; `codon_selector.py` returns the ambiguous codon, base sets, expansion, amino-acid/codon records and score, with no mixture-fraction field. Verified at upstream HEAD `c26439f` on 2026-08-14. |
| **E3 — balance (`fc-C`)** | **DECLINED** | No balance/emphasis edit, per ruling A. Phrases removed under E1 were removed for prohibited provenance, not to rebalance the record. | `fc-C` disputes emphasis rather than a checkable fact; the authors' ruling expressly declines it. |

### Provenance, labels, and drift follow-through

- **`cit-21` is now resolved rather than rejected.** Pass 1 could not verify the web-wide negative “not deployed”; Pass 2 removes the third-party material and the negative entirely.
- **Ambiguous labels were made consistent:** “five modules” is now “five implementation modules” everywhere relevant (with docstring-only `__init__.py` still identified at the anchor); the GitHub redirect is stated as first-hop 301 and destination 200; “GitHub releases” is distinguished from PyPI releases; `maintained` is defined by upstream commit/package-release activity and distinguished from deployment availability; the least-confident DMS/provenance summaries are scoped to full-sequence output/the exposed design operation.
- **Unverifiable anecdote dropped:** “the first id most users would try” was removed from the taxid-562 limitation. Fork-awareness, deployment and referee-risk anecdotes disappeared with E1.
- **Current drift check:** PyPI's current release remains **1.7 (2022-03-29)**, matching the record; upstream HEAD remains **`c26439f` (2022-12-12)**; GitHub still has no releases or tags; the old repository URL still returns 301 to the current upstream; the appspot home still returns 200 and the paper URL still returns 404. No material release drift exists, so ruling D requires no parenthetical. Section numbering, headings and order were preserved.

### Value basis

Removing third-party and sibling material does not weaken a value below its stated basis:

- `library_as_object = no` remains supported by the paper's user-written per-position loop and upstream `CodonSelector.optimise_codons` accepting a single amino-acid set.
- `mixed_mutagenesis_one_pool = no` remains supported by the two official `/codons` dispatch modes and the lack of a pooled multi-mode carrier in upstream source.
- `hgvs_input = no` remains supported by the three official request arguments and absence of a variant-nomenclature parser upstream.
- `maintained = partial` remains supported under the now-explicit C1 threshold by the fetched 2022-12-12 upstream commit and 2022-03-29 PyPI release.

No capability-value escalation is required; no value was changed.

### Row-substitution report (no action)

- The locked reference already flags **`Codon / amino-acid substitutions`** as a likely near-uniform row. A CodonGenie-specific candidate would be **`Degenerate/ambiguous codon design`**, replacing that row; primary evidence is the paper Abstract and pp. 2–5 plus upstream `codon_selector.py`, which designs and expands an IUPAC ambiguous codon from a requested amino-acid set.
- **Not actionable here:** CodonGenie has no Table 2 column, and the permitted CodonGenie primary sources cannot establish discrimination across the locked eight columns. The old claim that PoolParty would be `no` depended on sibling evidence and was removed. No substitution was made.
- The locked reference also flags **`Recombination and chimeras`** as a possible one-positive row. This CodonGenie record supplies no cross-tool evidence that could resolve that risk; no action.

### Neighbour tensions

1. Table 1's grouped adjacent-tools **Output** cell says “An optimised sequence or a primer set,” whereas CodonGenie's documented web output is a ranked list of ambiguous codons and expansions. Its grouped **Purpose** (“Optimise or edit individual sequences”) likewise compresses the paper's single-codon design scope. This is reported, not edited: the binding table is outside the allowed write root for this task, and CodonGenie is one member of a grouped row.
2. The Pass 1 summary and “Tensions left in place” remain a historical log. Pass 2 supersedes its E1–E3 open status and resolves its first two terminology tensions (`pool-composition skew` → encoding redundancy; `five modules` → five implementation modules) without rewriting the historical entries.

### Pass 2 rejection/escalation inventory

- **Rejected/declined:** `E3` balance/emphasis, by explicit ruling A.
- **Escalated:** none.
- **Row candidate declined for shipping:** `Degenerate/ambiguous codon design`, because CodonGenie is absent from Table 2 and cross-tool discrimination is unverified under the permitted provenance policy.
