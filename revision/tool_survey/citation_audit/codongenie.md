# Citation-integrity audit — CodonGenie

Record audited: `final/codongenie.md`

Sources checked: upstream `neilswainston/CodonGenie` at `c26439f93954d9b61b7dfccc9fd12b4936dfbcae`; `RidgeBio/CodonGenie` at `d4db59333c25c37b116360587c44aa4cebdf5c14`; all 10 pages of `papers/Swainston2017_codongenie.pdf` extracted with PyMuPDF; the official PyPI JSON plus the 1.7 wheel and sdist; GitHub repository/releases/tags/issues/pulls metadata; the official Biopython/NCBI codon-table definitions; the live CodonGenie, Kazusa, NCBI, SYNBIOCHEM, GitHub, and CDN URLs; and fresh reproductions of every API example on 2026-08-14. This audit concerns citation accuracy only, not any capability rating or conclusion.

## NOT FOUND IN ANY SOURCE

None. The two non-verbatim paper quotations below have identifiable underlying passages and are therefore classified `misquoted`, not `NOT FOUND IN ANY SOURCE`.

## Findings (every item not verified)

1. **Status: misquoted. Record line 243.** `"Written in Python/Flask, ∼700 lines"` is not verbatim on PDF p. 6, despite the explicit claim that it is. It is a composite paraphrase of two separate passages: p. 4 says CodonGenie is “written in Python (using the Flask framework) and HTML/Javascript,” while p. 6 says “the entire application consists of ∼700 lines of code.” The combined quoted wording occurs nowhere in the paper.

2. **Status: misquoted. Record line 120.** The quoted sentence `"The codon returned is therefore the most frequent codon for encoding proline in E. coli"` silently corrects the paper. PDF p. 7 actually reads, with the authors’ typo, “The codon returned is **the therefore** the most frequent codon for encoding proline in E. coli.” The substantive claim is unchanged, but the quotation is not verbatim.

3. **Status: misquoted. Record lines 64–72.** The purported top-hit JSON for `aminoAcids=FILMV&organism=37762` splices data from two different organism IDs. Fresh output for taxid 37762 has TTT `probability: 0.6358653093187157`, not `0.567`; `0.5674820143884892` is the TTT probability returned for taxid **83333**, the endpoint used in the earlier extraction. The 37762 expansion and score (`0.8777387460423176`) otherwise match. The displayed block is therefore not an accurate quotation of the stated endpoint.

4. **Status: minor-discrepancy. Record lines 49–52.** “Its entire public surface is two operations” is too broad. Design and Analyse are the two core codon operations, but the public service also exposes organism listing and organism search, and `CodonGenieClient` correspondingly exposes four methods: `get_codons`, `analyse`, `get_organisms`, and `search_organisms`. The record itself later identifies the organism routes.

5. **Status: minor-discrepancy. Record line 51.** The claim that Design enumerates and ranks *every* ambiguous codon encoding the requested amino acids is contradicted by the implementation and live API. For `aminoAcids=P`, the live service returns only the four unambiguous codons `CCG`, `CCT`, `CCA`, and `CCC`; it omits `CCN`, `CCR`, `CCY`, and other ambiguous codons that also encode only proline. The algorithm generates mergers of one selected encoding per required amino acid, not all 3,375 possible IUPAC triplets. The paper’s Abstract uses similarly broad “all ambiguous codons” wording, so this is an implementation/evidence mismatch rather than an invented quotation.

6. **Status: wrong-location. Record lines 54 and 112.** `main.py` does define four decorated routes, but not “every response” is `Response(json.dumps(...), mimetype='application/json')`, and not every Flask route returns JSON. `/` returns `app.send_static_file('index.html')` and is live as `text/html`; the exception handler uses `jsonify`. The paper’s “In all cases, results are returned in json format” occurs specifically in the web-service-results discussion, not as a claim about the home/static route. The three data-route handlers do return JSON, so the narrower API-output claim is verified.

7. **Status: number-wrong. Record lines 18, 60, 91, 199, and 221.** There are **six** `.py` files directly under `codon_genie/`, not five: the five named implementation files plus `__init__.py`. The wheel likewise contains those six top-level Python files. “Five non-`__init__` implementation modules” would be accurate; the two `test_*.py` modules and every named implementation file are otherwise correctly identified.

8. **Status: dead-link. Record line 37.** The stated tree reference `github.com/neilswainston/CodonGenie/git/trees/master?recursive=1` is not a live GitHub endpoint; with HTTPS it returns 404. The working API URL is `https://api.github.com/repos/neilswainston/CodonGenie/git/trees/master?recursive=1`, which returns 200 and supports the tree claims.

9. **Status: minor-discrepancy. Record lines 89, 90, 92, and 100.** Several absolute source-wide negatives are contradicted by shipped code that the same record later cites. “No random sampling” conflicts with `CodonOptimiser.get_random_codon()` and `mutate()`; “no sequence context at all,” “There are no output sequences,” and “never emits sequences” conflict with `get_codon_optim_seq()` and `example/align.py::seq_from_alignment()`, each of which returns a full DNA sequence. These negatives are accurate only when explicitly scoped to the user-facing single-codon `/codons` operation, not to the repository/installable package.

10. **Status: wrong-location. Record line 101.** `static/index.html` does not contain the claimed regex literal `[acgtmrwsykvhdbnACGTMRWSYKVHDBN]{3}`. It contains `data-ng-pattern="ctrl.codon_pattern"`; the literal is defined in `static/codonGenie/codonGenie.ctrl.js:5`, a file omitted from the row’s source citation and from the record’s list of raw files read. The regex claim itself is real.

11. **Status: number-wrong. Record lines 94, 127, and 209.** `static/index.html` loads **11 external asset URLs**, not five. They are distributed across five hostnames, which likely explains the record/reviewer count. Fresh requests found all 11 URLs live with HTTP 200, so link health is verified; “five CDN assets/dependencies” is the inaccurate count.

12. **Status: minor-discrepancy. Record line 142.** The evidence does not establish full “Alternative genetic codes” support. `CodonSelector(table_id=...)` does select `Bio.Data.CodonTable.unambiguous_dna_by_id[table_id]` for decoding, but candidate generation still uses the hard-coded standard-code `CODONS` mapping. For example, mitochondrial table 2 treats TGA as W and AGA/AGG as stops, while the hard-coded mapping still gives W only TGG and includes AGA/AGG under R. The parameter exists, but alternative-code design is internally inconsistent; the paper describes augmented genetic codes/alphabets as **future** adaptation.

13. **Status: uncited. Record line 139.** The JSON and UI verify per-amino-acid codon multiplicity and organism codon-usage `probability`, but they do not state “expected variant-composition skew within the designed pool.” Codon-usage probability is a host-expression statistic, not a measured or predicted synthesis abundance. A 2:1 encoding multiplicity can support a conditional inference under an equal-mixture assumption, but no cited source states that assumption or the claimed expected pool composition.

14. **Status: number-wrong. Record line 164.** The fresh taxid-1902 DTS score is `0.7782896589104741`. To two decimal places this rounds to **0.78**, not the paper’s 0.79. DTK `0.29425093243948847` does round to 0.29. The ranking reversal is verified, but “the paper’s 0.79 / 0.29 to rounding” is numerically false for the current DTS result.

15. **Status: minor-discrepancy. Record line 176.** `2025-05-17T04:32:13Z` is the RidgeBio fork’s last commit timestamp, not its last-push timestamp. GitHub reports `pushed_at: 2025-05-17T04:32:35Z`, 22 seconds later. The date, 88-versus-86 commit counts, added methods, and UI changes are verified.

16. **Status: minor-discrepancy. Record line 198.** The cited failure remains HTTP 500, but its quoted payload did not reproduce: eight fresh requests returned `{"message": "KeyError: 'ATG'"}`, not `KeyError: 'ATC'`. Which absent codon is encountered can depend on set/hash iteration order, so the exact key is unstable; the substantive unsupported-taxid failure is verified.

17. **Status: wrong-location. Record lines 129 and 211.** The last human-authored commit is indeed dated 2022-06-09, but its subject is **“Minor auto-reformatting.”** (`6ef2c3e`, 16:20:09 +01:00). **“Update to remove dependency on remote species.table flatfile.”** is the preceding same-day commit (`e38cc8e`, 13:57:08 +01:00), not the last human commit. The date/age claim is unaffected.

18. **Status: minor-discrepancy. Record line 197.** A Kazusa format change or outage would not “silently” remove scoring. `urlopen` exceptions propagate to the Flask exception handler as HTTP 500; a response with no parsable `<PRE>` content produces empty lookup tables and then a raw `KeyError`/500. The plain-HTTP live screen-scrape dependency and `<PRE>` parser are verified, but the described failure mode is not silent.

19. **Status: number-wrong. Record line 221.** The record has **24** capability cells, not 23. The review says one verdict should change (`maintained`), but does not literally mark “22 of 23” verdicts supported. If summarized arithmetically, one changed verdict out of 24 means 23 of 24 stood, not 22 of 23.

20. **Status: uncited. Record lines 127, 213, and 215.** Existence of `Dockerfile`, `start_server.sh`, and `app.yaml` supports a documented self-hosting path, but not the stronger current-runtime claims that clone+Docker is a “working reproduction path” and the tool is “reproducible via clone+Docker.” The record explicitly states at line 235 that it did not install or execute the tool, and it retains no container build/run result. The live hosted API is independently verified; this finding concerns only self-hosted reproducibility.

21. **Status: uncited. Record lines 86 and 181.** The RidgeBio code is demonstrably unmerged upstream, and the official CodonGenie PyPI artifacts predate it, but “undeployed” / “the fork is not deployed” is a web-wide negative with no deployment search, registry query, or cited host list. No deployment was found in the checked repository metadata, but the record does not supply evidence sufficient for the absolute claim.

22. **Status: uncited. Record lines 137, 146, 154, and 246.** Cross-tool claims are unsupported in this CodonGenie record: PoolParty allegedly has no analogous Analyse mode, enumerates explicit sequences, and does not use degenerate codons; MPRADesignTools and MPRAnator allegedly have the stated dormancy durations. None has a PoolParty/MPRADesignTools/MPRAnator source citation here. Line 246 itself appropriately says “if true,” but the earlier occurrences are asserted as facts.

23. **Status: wrong-location. Record lines 20 and 22.** The reconciliation’s internal section pointers are wrong. The corrected documented examples are in **§6**, not §3, and the added limitations are in **§8**, not §6.

24. **Status: minor-discrepancy. Record lines 36, 208, and 210.** The old GitHub repository URL is live through a **301** redirect whose final target is HTTP 200; wording it as an “HTTP-redirect ... (200)” conflates the first-hop and final statuses. Separately, the HTTPS certificate mismatch for `codon.synbiochem.co.uk` is verified, but the served certificate is not “only Google's `*.appspot.com` wildcard”: its CN is `*.appspot-preview.com` and it has many Google/Appspot SANs, none covering `codon.synbiochem.co.uk`. These details do not change the record’s correct dead-published-link conclusion.

## Link check

- Correctly reported dead/broken: `http://codon.synbiochem.co.uk/` and its documented API example return 404; HTTPS fails hostname validation.
- Live: the upstream and RidgeBio repositories, the redirected old repository URL, PyPI JSON, the `align.py` blob, `https://codongenie.appspot.com`, SYNBIOCHEM’s design directory, Kazusa’s working 37762 page, and the NCBI taxdump URL.
- The live UI contains 11 external asset URLs and all 11 returned 200.
- The only malformed source link found is the non-API Git-tree path in finding 8.

## Verified numerical anchors

- PDF: 10 pages; bibliographic metadata and DOI are correct.
- Upstream: 86 commits, 1 star, 3 forks, 0 open issues, 0 GitHub releases, 0 tags, not archived; HEAD 2022-12-12.
- RidgeBio fork: 88 commits; last commit 2025-05-17T04:32:13Z; `pushed_at` 2025-05-17T04:32:35Z.
- PyPI: 8 versions from 0.0.1 (2020-01-17) through 1.7 (2022-03-29); `requires_python >=3.7`; `requires_dist: null`; wheel and sdist omit `species.table` and `main.py` and contain both `test_*.py` modules.
- Live API: 35,792 organisms; FILMV@37762 gives 12 results, DTK 0.8777387460 then DTS 0.6836806630; FILMV@1902 gives DTS 0.7782896589 then DTK 0.2942509324; DBK has 18 expansions; NSS has 16 expansions and 8 amino acids; DE@4932 gives 4 results headed by GAW.
