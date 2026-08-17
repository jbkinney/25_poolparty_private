# CodonGenie — evidence memo

**Survey date:** 2026-08-10
**Slug:** `codongenie`
**Full name:** CodonGenie: optimised ambiguous codon design tools
**Citation:** Swainston N, Currin A, Green L, Breitling R, Day PJ, Kell DB (2017). *PeerJ Computer Science* 3:e120. DOI 10.7717/peerj-cs.120

---

## 1. Sources consulted

| Kind | Reference | Notes |
|---|---|---|
| pdf | `/revision/tool_survey/papers/Swainston2017_codongenie.pdf` (10 pp) | Full text extracted with PyMuPDF |
| prior_analysis | `/revision/tool_survey/prior_analyses/Swainston2017_codongenie_analysis.md` | Broadly correct; see §7 for corrections |
| repo | https://github.com/neilswainston/CodonGenie | Paper cites `synbiochem/CodonGenie`, which **now HTTP-redirects** to `neilswainston/CodonGenie` (verified, 200) |
| repo | GitHub REST API `repositories/70245462` | metadata, commits, tags, tree |
| repo files | `main.py`, `codon_genie/codon_selector.py`, `codon_genie/codon_utils.py`, `codon_genie/seq_utils.py`, `codon_genie/client.py`, `codon_genie/example/align.py`, `codon_genie/test/test_codon_genie.py`, `static/index.html`, `README.md`, `setup.py`, `requirements.txt` | read raw from `master` |
| pypi | https://pypi.org/pypi/CodonGenie/json | v1.7, 8 releases, last 2022-03-29 |
| web | http://codon.synbiochem.co.uk (paper's advertised URL) | **DEAD — Google App Engine HTTP 404** |
| web | https://codongenie.appspot.com | **LIVE — HTTP 200, REST API returns correct results** (tested 2026-08-10) |
| docs | none found | There is no documentation site. The README is 6 lines of Docker instructions. The paper *is* the documentation ("CodonGenie can provide a simple, easy-to-use interface that **requires no documentation**", p. 6) |

---

## 2. What the tool actually is

CodonGenie is a **single-codon** design microservice. Its entire public surface is two operations:

- **Design** — given a set of target amino acids + an NCBI Taxonomy ID, enumerate and rank *every* ambiguous (degenerate) codon that encodes those amino acids.
- **Analyse** — given an existing ambiguous codon + organism, report which amino acids it encodes and at what codon-usage frequency.

From `main.py`, the app has exactly **four** routes: `/`, `/organisms/`, `/organisms/<term>`, `/codons`. That is the whole tool.

> "the entire application consists of ∼700 lines of code" (p. 6)

Scoring (p. 3) is the mean, over all unambiguous codons in the ambiguous codon's expansion, of the codon's frequency relative to the most-frequent synonymous codon — with **zero assigned to any codon encoding a non-required amino acid or a stop**. Confirmed in code (`codon_selector._get_score`, which uses the `cai` field).

### Live API output shape (verified)

`GET https://codongenie.appspot.com/codons?aminoAcids=FILMV&organism=83333` returns a ranked JSON array; top hit:

```json
{ "ambiguous_codon": "DTK",
  "ambiguous_codon_expansion": ["ATG","ATT","GTG","GTT","TTG","TTT"],
  "ambiguous_codon_nucleotides": ["AGT","T","GT"],
  "amino_acids": [ {"amino_acid":"F","codons":[{"cai":1.0,"codon":"TTT","probability":0.567}],"type":1}, ... ],
  "score": 0.8151852092922643 }
```

The `type` field is the only "consequence"-like annotation: from `codon_selector._get_amino_acid_type`, `-1` = Stop, `0` = additional/off-target amino acid, `1` = required amino acid.

`GET /organisms/` returns **35,792** organisms — exactly matching the paper's "as of May 2017 provided support for 35,792 organisms" (p. 4).

---

## 3. Capability assessment

### Block A — library specification

| key | value | evidence |
|---|---|---|
| `library_as_object` | **no** | The output object is a ranked list of *ambiguous codons for one position*. There is no multi-sequence object anywhere in the codebase. The paper's own worked example requires the user to write the loop: "By iterating through the alignment, the set of amino acids required at each position can be collected... and a synthetic DNA sequence **built up** from the highest-scoring ambiguous codon returned" (p. 7). `example/align.py` implements this as `dna_seq += codon['ambiguous_codon']` inside a manual `for position in zip(...)` loop. |
| `dag_chaining` | **no** | Deliberately anti-compositional at the software level: "CodonGenie is designed to follow the concept of 'microservices'... breaking down of large, monolithic applications into simple, atomic services of limited scope" (p. 6). `/codons` is a single stateless GET. No pipeline, graph, or operation-composition primitive exists. |
| `lazy_evaluation` | **no** | `CodonSelector.optimise_codons` eagerly builds `itertools.product` over all encodings, and each candidate is fully expanded: `ambiguous_codon_expansion` materialises every unambiguous codon. *Nuance worth noting for the referee:* the ambiguous codon is itself a compressed physical representation of a variant set (NTN = 16 DNA variants synthesised as one oligo), but that is degeneracy in the wet lab, not lazy evaluation in software. |
| `mixed_mutagenesis_one_pool` | **no** | One and only one mutagenesis mode: substitution of a codon by a degenerate codon encoding a user-chosen amino acid set. No indels, no WT/replicate handling, no random sampling, no pairwise/exhaustive distinction, and no container in which such modes could be mixed. |
| `combinatorial_motif_place` | **no** | No concept of a motif, an element, a position within a larger sequence, or an orientation. Input is `aminoAcids` (a string of AA letters) and `organism`. |
| `barcode_generation` | **no** | No barcode, UMI, GC-constraint or edit-distance code anywhere in the repo. |
| `per_sequence_provenance` | **partial** | There are no output *sequences*, so no per-sequence provenance in PoolParty's sense. But the per-codon JSON is genuinely structured and explanatory: `ambiguous_codon_nucleotides`, `ambiguous_codon_expansion`, per-amino-acid `codons` with `cai` and `probability`, `type` (required / off-target / stop), and `score`. Marked partial, not no, because that record does explain *why* the design was chosen — it just annotates a codon, not a library member. |
| `automatic_naming` | **no** | Nothing is named. Results are keyed by the ambiguous codon string itself; there is no FASTA/ID emitter in `main.py` or `codon_selector.py`. |
| `design_visualization` | **partial** | The web UI renders a colour-coded interactive result table, per Fig. 1 caption (p. 5): "Variant codons are shown in grey, with their encodings shown in **green, orange and red for required amino acids, additional amino acids and stop codons**, respectively... these encodings and their codon usage frequencies may be visualised through a **tooltip**." Confirmed in `static/result/result.directive.js` + `result.css`. This visualises a *single codon design*, not a library or a design graph. |

### Block B — assay coverage

| key | value | evidence |
|---|---|---|
| `assay_dms` | **partial** | Squarely aimed at coding-sequence variant libraries: "designing ambiguous codons to support protein mutagenesis applications" (Abstract); motivated by directed evolution and by "the study of sequence-to-fitness relationships (Hietpas, Jensen & Bolon, 2011)" (p. 1) — i.e. DMS-adjacent. But it designs a codon, not a DMS library: it never enumerates variants, never emits sequences, and has no notion of scanning every position of an ORF. Partial = solves one sub-step of coding-library design. |
| `assay_mpra` | **no** | Entirely amino-acid/codon-centric. The input alphabet is the 20 amino acids (`static/index.html` renders 20 AA buttons grouped polar/non-polar/acidic/basic). No promoter, enhancer, UTR, or non-coding element concept exists. |
| `assay_insilico` | **no** | No mention of models, prediction, or in-silico probing in the paper or repo. Published 2017. |

### Block C — genomics integration

| key | value | evidence |
|---|---|---|
| `genome_coordinates` | **no** | The only "genomic" input is an NCBI **Taxonomy** id used to fetch a codon usage table (`codon_utils.CodonOptimiser.__get_codon_usage` queries `kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=<taxid>`). That is an organism selector, not a coordinate system. No chromosome/start/end anywhere. |
| `transcript_models` | **no** | No GTF/GFF/BED parsing. `requirements.txt` is two lines: `Flask`, `Biopython`. |
| `exon_intron_split_codons` | **no** | Codons are handled as isolated triplets (`CODONS` dict in `codon_selector.py` maps each AA to `[pos1, pos2, wobble-set]`). No exon, intron, or splicing concept. |
| `hgvs_input` | **no** | Input is a bare string of one-letter amino acid codes (`request.args['aminoAcids']`) or a 3-character ambiguous codon. No variant-nomenclature parser. |
| `vcf_vep_output` | **no** | "In all cases, results are returned in **json** format" (p. 3). Confirmed: `main.py` returns `Response(json.dumps(...), mimetype='application/json')` only. No VCF, no FASTA, no CSV export. |
| `consequence_annotation` | **partial** | It does classify what a design will produce: `_get_amino_acid_type` returns `-1` for Stop, `0` for an additional (off-target) amino acid, `1` for a required amino acid, and this is surfaced in both the JSON and the colour-coded UI. So stop-codon inclusion and off-target amino acids are explicitly annotated. It is **not** variant-consequence annotation against a reference transcript (no synonymous / missense / frameshift / in-frame-indel calling). Partial. |

### Block D — physical construction

| key | value | evidence |
|---|---|---|
| `primer_design` | **no** | No oligo, primer, Tm, or overlap code in the repo (`seq_utils.py` contains reagent-concentration constants but the shipped `main.py`/`codon_selector.py` never call any Tm routine). The paper defers this to the authors' *other* tool: "CodonGenie is amenable to future integration with new and existing variant library design software tools (Swainston et al., 2014)" — i.e. GeneGenie, the oligomer designer. |
| `codon_optimization` | **yes** | Two independent grounds. (1) It is the core scoring axis: codons are "ranked according to their efficiency in encoding the required amino acids while minimising the inclusion of additional amino acids and stop codons. **Organism-specific codon usage is also considered**" (Abstract); the DTK-vs-DTS *E. coli* vs *S. coelicolor* reversal (p. 4, Table 1) is a pure codon-usage result. For invariant positions the tool returns the plain best codon: "The codon returned is therefore the most frequent codon for encoding proline in E. coli" (p. 7). (2) The shipped package contains a full sequence optimiser: `codon_utils.CodonOptimiser.get_codon_optim_seq(protein_seq, excl_codons, max_repeat_nuc, restr_enzyms, ...)` plus `get_cai`, `get_best_codon`, `get_random_codon`. **Caveat:** (2) is vendored from `liv-utils` and is *not* exposed by the web UI or the REST API. |
| `synthesis_constraints` | **partial** | `seq_utils.find_invalid(seq, max_repeat_nuc, restr_enzyms)` checks homopolymer runs and restriction sites (Biopython `Restriction`), and `get_codon_optim_seq` backtracks to avoid them. Real, shipped code. But it is **not reachable from the CodonGenie web app or REST API** and is never applied to ambiguous-codon design. What the user-facing tool *does* report is library size (`len(ambiguous_codon_expansion)`, e.g. "DTK... using six variants" vs NTN's 16, p. 4) — a screening/throughput constraint, not a per-sequence synthesis constraint. Hence partial. |

### Block E — engineering

| key | value | evidence |
|---|---|---|
| `interface` | **yes** | Four surfaces, **no CLI**: (a) single-page **web GUI** (AngularJS 1.5 + Bootstrap 3.3, `static/index.html`); (b) **RESTful JSON web service** — `/codons?aminoAcids=&organism=`, `/codons?codon=&organism=`, `/organisms/`, `/organisms/<term>` (p. 3); (c) **Python client library** `CodonGenieClient` (`codon_genie/client.py`) with `get_codons`, `analyse`, `search_organisms`, `get_organisms`; (d) **pip-installable package** `CodonGenie` and a **Docker** image (`Dockerfile`, `start_server.sh`). |
| `license` | **yes** | **MIT.** `LICENSE` in repo; GitHub API `license.name = "MIT License"`; `setup.py` header: "CodonGenie is licensed under the MIT License." Paper itself is CC-BY 4.0. |
| `maintained` | **yes** | **Last commit 2022-12-12** (`Merge pull request #1 from TrellixVulnTeam` — an automated security patch). Last *substantive* author commit **2022-06-09** ("Update to remove dependency on remote species.table flatfile"). Last PyPI release **1.7, 2022-03-29**. GitHub `pushed_at = 2022-12-12T11:42:23Z`. **Zero** GitHub releases and **zero** tags. 86 commits total, 1 star, 3 forks, 0 open issues. Not archived. |

---

## 4. Documented examples / vignettes (candidates for PoolParty reproduction)

1. **MSA → degenerate synthetic gene** (the flagship example, paper pp. 7–8 + `codon_genie/example/align.py`). Alignment `['PFDMR','PIAMR','PLHLR','PMNMR','PVHMR']`, host *E. coli* → **`CCG|DTK|VMT|MTG|CGT`**. Runnable script:
   https://github.com/neilswainston/CodonGenie/blob/master/codon_genie/example/align.py
   This is the one example that produces something resembling a library specification, and it is the natural head-to-head target for PoolParty.
2. **Non-polar set F, I, L, M, V** (p. 4). Naive `NTN` = 16 DNA variants; CodonGenie finds **`DTK`** ([AGT][T][GT]) = **6** variants. Scaling claim: "when using 3 DTK codons the library size is reduced from 4,096 (16³) to 216 (6³)" *(the paper's PDF prints "213", a typesetting corruption of 6³ = 216)*.
3. **Host-dependence of the same set** (p. 4, Table 1). *E. coli*: DTK score 0.88 > DTS 0.68. *S. coelicolor*: reversed — DTS 0.79 > DTK 0.29. Driven by wobble-position T vs C preference (e.g. Phe TTT/TTC = 0.64/0.36 in *E. coli* vs 0.03/0.97 in *S. coelicolor*).
4. **Fig. 1 GUI example** (p. 5). Non-polar residues **A, F, G, I, L, M, V** → preferred codon **`DBK`** = [AGT][CGT][GT] = 18 DNA variants; V is encoded twice (GTG and GTT).
5. **Analyse mode** (p. 3). `GET /codons?codon=NSS&organism=4932` (*S. cerevisiae*) — verified live; returns 16 codons and 9 encoded amino acids, all typed `0` (off-target, since no target set was given).
6. **Unit test** (`codon_genie/test/test_codon_genie.py`): `CodonSelector().optimise_codons('FLIMV', '37762')` → asserts 12 results, best = `DTK`.
7. **Contrast case vs DYNAMCC** (p. 5). For F/I/L/M/V, DYNAMCC returns the *set* {`WTT`, `VTG`} = 5 variants with no off-target; CodonGenie returns a *single* codon DTK/DTS = 6 variants. This is CodonGenie's own statement of its scope boundary.

---

## 5. Notable capabilities NOT covered by the current row list

- **Ranked trade-off curve rather than a single answer.** "CodonGenie returns not only the most 'specific' ambiguous codons... Providing results that include less specific ambiguous codons, which may also encode additional amino acids, allows the user to perform a **trade-off between library size and codon specificity**" (pp. 2–3). Results are sorted by (n variants ascending, score descending).
- **Exact library-size arithmetic per degenerate codon** — `ambiguous_codon_expansion` gives the precise DNA-variant count, which is the screening-burden currency in directed evolution.
- **Reverse / "Analyse" direction** — decode an existing degenerate codon into encoded amino acids + usage frequencies. Useful as a validator; PoolParty has no analogous audit mode for user-supplied degenerate codons.
- **Off-target and stop-codon flagging** as a first-class output field.
- **Organism resolution service** — free-text organism search → NCBI Taxonomy ID over 35,792 organisms (verified live), backed by the Kazusa Codon Usage Database.
- **Alternative genetic codes** — `CodonSelector(table_id=1)` accepts any NCBI codon table id via `Bio.Data.CodonTable.unambiguous_dna_by_id`; not exposed in the API but present. Paper flags "augmented genetic codes and expanded genetic alphabets" as future work (p. 7).
- **CAI computation for an arbitrary DNA sequence** (`get_cai`) and stochastic codon sampling (`get_random_codon`, `mutate(protein_seq, dna_seq, mutation_rate)`) in the shipped utils.
- **Deployability as a microservice** — Docker image + `app.yaml` for Google App Engine, explicitly so it can be embedded in larger design pipelines.

**Row-list suggestion:** if the matrix adds a row it should be *degenerate/ambiguous codon design*, where CodonGenie is the reference implementation and PoolParty would be "no" (PoolParty enumerates explicit sequences). That framing is honest and strengthens rather than weakens the comparison.

---

## 6. Stated limitations (the tool's own words)

- Scope is deliberately one codon: "CodonGenie designs **single** ambiguous codons to encode a desired set of amino acids, **which may also include off-target amino acids**" (p. 5) — contrasted with DYNAMCC, which "designs **sets** of ambiguous codons... with no off-target amino acid encoding and minimal redundancy."
- Efficiency cost of that choice is acknowledged: "five DNA variants encode the five desired amino acids [DYNAMCC], while CodonGenie's solution of DTK or DTS encode six DNA variants, thus producing a **larger library**" (pp. 5–6).
- Library assembly is out of scope and deferred to sister tools: "CodonGenie is amenable to **future integration** with new and existing variant library design software tools (Swainston et al., 2014)" (p. 7).
- Minimalism is a design principle, not an oversight: "~700 lines of code"; "requires no documentation" (p. 6).
- Only 17 of 20 amino acids are cleanly handled by the fast algorithm; L, R, S require the algorithm to be run per-encoding and combined (p. 2).

---

## 7. Availability TODAY (verified 2026-08-10)

| Asset | Status |
|---|---|
| **Paper's advertised URL** `http://codon.synbiochem.co.uk` | **DEAD.** Returns Google App Engine `Error 404 (Not Found)!!1` from `142.251.45.180`; the HTTPS cert served is Google's `*.appspot.com` wildcard, so the custom domain mapping has lapsed. Both `/` and the paper's documented API example `?aminoAcids=DE&organism=4932` 404. |
| **Working deployment** `https://codongenie.appspot.com` | **LIVE, HTTP 200.** Full AngularJS UI served; REST API functional and returning correct, paper-consistent results (FILMV/*E. coli* → DTK top-ranked, score 0.815; `/organisms/` returns 35,792 entries). This is the URL hard-coded as the default in `client.py`, but it appears **nowhere in the paper** — a reader following the publication lands on a dead link. |
| **Repo** github.com/synbiochem/CodonGenie | Redirects (200) to **github.com/neilswainston/CodonGenie**. Public, not archived, MIT. 86 commits, 1 star, 3 forks, 0 open issues, **0 releases, 0 tags**. |
| **Last activity** | Commit **2022-12-12** (automated dependency/security merge); last human feature commit **2022-06-09**. ~3.5 years dormant. |
| **PyPI** `CodonGenie` | Installable, **v1.7 (2022-03-29)**, 8 releases since 2020-01-17, `python_requires>=3.7`, deps `Flask` + `Biopython`. |
| **Self-hosting** | `Dockerfile` + `start_server.sh` + `app.yaml`. Should work, but note a live external dependency: `CodonOptimiser.__get_codon_usage` scrapes `http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi` (plain HTTP, HTML `<PRE>` screen-scrape) at request time. If Kazusa changes its page format or goes offline, a self-hosted instance silently loses all codon-usage scoring. |

**Verdict:** functionally alive and reproducible, but effectively unmaintained and the published entry point is broken.

---

## 8. Corrections / notes on the prior analysis

The prior analysis is substantially correct. Corrections and additions:

1. **Repo URL is stale.** The paper (and the prior note) cite `synbiochem/CodonGenie`; the repo moved to `neilswainston/CodonGenie`. Cite the redirect-target if the manuscript names a URL.
2. **The advertised web URL is dead.** `codon.synbiochem.co.uk` 404s. The prior note describes the tool as if the service were reachable at that address. If PoolParty's manuscript or rebuttal says "freely available from...", it should point to `codongenie.appspot.com` or note the status.
3. **"Written in Python/Flask, ~700 lines"** — confirmed verbatim from p. 6.
4. **Prior note says "Operates at the single-codon level, not at the level of full sequence libraries"** — confirmed and well-supported. Strongest single citable line for the rebuttal is the paper's own microservice statement (p. 6) plus the fact that its multi-position example is a user-written `for` loop in `example/align.py`.
5. **Prior note overlooks** the shipped `CodonOptimiser` / `find_invalid` utilities. These matter for the `codon_optimization` and `synthesis_constraints` cells and a tool-author referee would raise them. Handled above with the "present in package, not exposed in the tool" caveat.
6. **Prior note's claim "PoolParty does not use degenerate codons"** — if true, this is worth stating explicitly in the comparison, because it makes the two tools complementary rather than competing, which is the safest posture given a possible author-referee.
7. Minor: the paper's "213" for 6³ is a PDF typesetting corruption of **216**. Do not quote "213".

---

## 9. Least-confident cells

- `per_sequence_provenance` (**partial**) — defensible either way. Rich structured metadata exists, but there are no sequences for it to attach to. Chose partial to avoid a referee objection that the JSON output was ignored.
- `consequence_annotation` (**partial**) — the `type` field (stop / off-target / required) is genuine consequence-flagging at the codon level, but it is not variant-consequence annotation against a transcript. Partial with an explicit scope note.
- `synthesis_constraints` (**partial**) and `codon_optimization` (**yes**) — both hinge on `codon_utils.py` / `seq_utils.py`, which ship in the pip package but are unreachable from the web UI and REST API. If the matrix's convention is "user-facing capability only", `synthesis_constraints` should drop to **no** and `codon_optimization` remains **yes** on the codon-usage-ranking ground alone.
- `assay_dms` (**partial**) — arguably **no** if the matrix requires producing a library. Chose partial because site-saturation mutagenesis for coding sequences is unambiguously the tool's declared purpose.
- `design_visualization` (**partial**) — the colour-coded table is real and is design output, but it visualises one codon, not a library.
