# CodonGenie — FINAL capability record

**Slug:** `codongenie`
**Full name:** CodonGenie: optimised ambiguous codon design tools
**Citation key:** `Swainston2017rb`
**Citation:** Swainston N, Currin A, Green L, Breitling R, Day PJ, Kell DB (2017). CodonGenie: optimised ambiguous codon design tools. *PeerJ Computer Science* 3:e120. DOI 10.7717/peerj-cs.120
**Tier:** 3
**Record status:** FINAL — extraction + adversarial review reconciled, all disputed facts re-verified from primary sources on 2026-08-10/11.

---

## 0. Reconciliation summary (what changed from the extraction)

| # | Cell / field | Extraction | FINAL | Why |
|---|---|---|---|---|
| 1 | `maintained` | yes | **partial** | Reviewer verdict "overstated" accepted. The extraction's own evidence said "effectively dormant"; value contradicted evidence. Resolved under the explicit maintenance convention (§5). |
| 2 | `primer_design` evidence | implied a Tm routine exists | reworded | Re-read `seq_utils.py` (66 lines): two functions only, no Tm/annealing/overlap code anywhere in the repo. Verdict `no` is *stronger* than the old wording. |
| 3 | `barcode_generation` evidence | "4 Python modules" | **5 implementation modules** | Verified from the git tree: `client.py`, `codon_selector.py`, `codon_utils.py`, `ncbi_tax_utils.py`, `seq_utils.py`. Verdict unaffected. |
| 4 | `interface` evidence | "pip-installable, deps Flask + Biopython" | qualified | PyPI `requires_dist` is `null`; the 1.7 wheel and sdist contain **no `species.table` and no `main.py`**. `pip install CodonGenie` cannot run the app. Verdict stays `yes` (GUI + REST + client + Docker). |
| 5 | documented examples | taxid 83333, "213 is a typesetting corruption", "9 amino acids", "only functional test" | corrected | All four re-verified and fixed (§6). |
| 6 | additional capabilities | — | +2 added | `ncbi_tax_utils.TaxonomyFactory`; per-amino-acid encoding redundancy (§4). |
| 7 | limitations | — | +2 added | HTTP 500 on non-Kazusa taxids; pip-install path is broken (§8). |
| 8 | provenance | — | official sources only | Third-party material removed; the paper and upstream source independently support the single-codon scope (§7). |

No cell moved in CodonGenie's disfavour except `maintained`, and that move is required for internal consistency of the matrix.

---

## 1. Sources consulted

| Kind | Reference | Notes |
|---|---|---|
| pdf | `revision/tool_survey/papers/Swainston2017_codongenie.pdf` (10 pp) | Text extracted with PyMuPDF; page refs below are PDF pages |
| process_record | `revision/tool_survey/prior_analyses/Swainston2017_codongenie_analysis.md` | Prior workflow context only; not used as factual evidence |
| process_record | `revision/tool_survey/extractions/codongenie.md` + `revision/tool_survey/reviews/codongenie.md` | Reconciliation inputs only; not used as factual evidence |
| repo | https://github.com/neilswainston/CodonGenie | Paper cites `synbiochem/CodonGenie`, which returns 301 to this URL; the destination returns 200 (re-verified 2026-08-14) |
| repo | `api.github.com/repos/neilswainston/CodonGenie/git/trees/master?recursive=1` | Full upstream tree; `master` remains `c26439f` (re-verified 2026-08-14) |
| repo files | `main.py`, `codon_genie/{client,codon_selector,codon_utils,ncbi_tax_utils,seq_utils}.py`, `codon_genie/example/align.py`, `codon_genie/test/{test_codon_genie,test_client}.py`, `static/index.html`, `static/result/{result.html,result.directive.js,result.css}`, `setup.py`, `requirements.txt`, `Dockerfile`, `app.yaml`, `README.md` | read raw from `master` |
| repo | https://github.com/neilswainston/CodonGenie/commit/c26439f93954d9b61b7dfccc9fd12b4936dfbcae | Current upstream HEAD, dated 2022-12-12 (re-verified 2026-08-14) |
| pypi | https://pypi.org/pypi/CodonGenie/json + the downloaded `CodonGenie-1.7-py3-none-any.whl` and `CodonGenie-1.7.tar.gz` | v1.7, 8 releases, last 2022-03-29; current release re-verified 2026-08-14; contents inspected (not installed) |
| web | http://codon.synbiochem.co.uk (the paper's advertised URL) | **DEAD — HTTP 404**, re-verified 2026-08-14 |
| web | https://codongenie.appspot.com | **LIVE — HTTP 200**, home page re-verified 2026-08-14; REST API exercised repeatedly (see §3) |
| docs | README + paper | README gives Docker run instructions; the paper describes the interface as one that "requires no documentation" (p. 6) |

---

## 2. What the tool actually is

CodonGenie is a **single-codon** design microservice. Its public surface is two codon operations, plus organism listing and search:

- **Design** — given a set of target amino acids + an NCBI Taxonomy ID, enumerate and rank the candidate ambiguous (degenerate) codons the algorithm builds from those amino acids' encodings. (The Abstract says CodonGenie "designs and analyses all ambiguous codons that encode the required amino acids", but the algorithm is "optimized for computational efficiency, compared to a brute-force examination of all possible ambiguous codons" (p. 2) and never tests all 3,375 IUPAC triplets: `?aminoAcids=P&organism=37762` returns only CCG, CCT, CCA, CCC — not CCN.)
- **Analyse** — given an existing ambiguous codon + organism, report which amino acids it encodes and at what codon-usage frequency.

`main.py` defines exactly **four** routes: `/`, `/organisms/`, `/organisms/<term>`, `/codons`. Every response *except* `/`, which serves `static/index.html`, is `Response(json.dumps(...), mimetype='application/json')`.

> "the entire application consists of ∼700 lines of code" (p. 6)

Scoring (p. 3) is the mean, over all unambiguous codons in the ambiguous codon's expansion, of each codon's frequency relative to the most-frequent synonymous codon — with **zero** assigned to any codon encoding a non-required amino acid or a stop. Confirmed in `codon_selector._get_score` (uses the `cai` field).

`codon_genie/` contains **five** implementation modules (plus a docstring-only `__init__.py`): `client.py`, `codon_selector.py`, `codon_utils.py`, `ncbi_tax_utils.py`, `seq_utils.py` — plus `example/align.py` and two test modules.

### Live API output shape (re-verified 2026-08-11)

`GET https://codongenie.appspot.com/codons?aminoAcids=FILMV&organism=37762` → 12 ranked results; top hit:

```json
{ "ambiguous_codon": "DTK",
  "ambiguous_codon_expansion": ["ATG","ATT","GTG","GTT","TTG","TTT"],
  "ambiguous_codon_nucleotides": ["AGT","T","GT"],
  "amino_acids": [ {"amino_acid":"F","codons":[{"cai":1.0,"codon":"TTT","probability":0.636}],"type":1}, ... ],
  "score": 0.8777 }
```

`type` is the only "consequence"-like field: `-1` = Stop, `0` = additional/off-target amino acid, `1` = required amino acid (`codon_selector._get_amino_acid_type`). Verified live: `?aminoAcids=WY&organism=37762` → `TRK` with `W:1, Y:1, C:0, Stop:-1`.

`GET /organisms/` returns **35,792** organisms — exactly matching the paper's "as of May 2017 provided support for 35,792 organisms" (p. 4). Re-verified live 2026-08-11.

---

## 3. Capability assessment (FINAL)

### Block A — library specification

| key | value | evidence | source |
|---|---|---|---|
| `library_as_object` | **no** | The output object is a ranked list of *ambiguous codons for one position*. `CodonSelector.optimise_codons(amino_acids, organism_id)` does `set(amino_acids.upper())`, so multi-position input is impossible by construction — `"FILMV,DE"` is one set, not two positions. There is no multi-sequence object anywhere. The paper's own multi-position example requires the user to write the loop: "By iterating through the alignment, the set of amino acids required at each position can be collected... and a synthetic DNA sequence **built up** from the highest-scoring ambiguous codon returned" (p. 7); `example/align.py` implements this as `dna_seq += codon['ambiguous_codon']` inside a hand-written `for` loop over `zip(...)`. | PDF p. 7; `codon_genie/codon_selector.py`; `codon_genie/example/align.py` |
| `dag_chaining` | **no** | Deliberately anti-compositional at the software level: "CodonGenie is designed to follow the concept of 'microservices'... breaking down of large, monolithic applications into simple, atomic services of limited scope" (p. 6, quoted verbatim and verified against the PDF). `/codons` is a single stateless GET. No pipeline, graph, or operation-composition primitive exists in any of the five implementation modules. | PDF p. 6; `main.py`; all of `codon_genie/` |
| `lazy_evaluation` | **no** | `optimise_codons` eagerly builds `itertools.product` over all per-amino-acid encodings, and `__analyse_ambig_codon` materialises `ambiguous_codon_expansion` for **every** candidate — verified in the live JSON, where all 12 FILMV results carry a full expansion. *Nuance worth keeping for the referee:* the ambiguous codon is itself a compressed physical representation of a variant set (NTN = 16 DNA variants synthesised as one oligo), but that is degeneracy in the wet lab, not lazy evaluation in software. | `codon_genie/codon_selector.py`; live `?aminoAcids=FILMV&organism=37762` |
| `mixed_mutagenesis_one_pool` | **no** | One and only one mutagenesis mode in the exposed tool: substitution of a codon by a degenerate codon encoding a user-chosen amino-acid set. API input is `aminoAcids` + `organism` (or `codon` + `organism`). No indels, no WT/replicate handling, no random sampling in the exposed design path (the unexposed `codon_utils.get_random_codon`/`mutate` do sample synonymous codons stochastically), no pairwise/exhaustive distinction, and no container in which modes could be mixed. *Reserve nuance:* a degenerate codon inherently pools WT with mutants and the tool explicitly tracks off-target amino acids (`type == 0`) — but that is one mutagenesis mode, not a mix of modes. | `main.py`; `codon_genie/codon_selector.py` |
| `combinatorial_motif_place` | **no** | The exposed tool has no sequence context at all. `CODONS` maps each amino acid to `[pos1, pos2, wobble-set]` and the algorithm operates on an isolated triplet. No motif, element, insertion position, orientation, or permutation concept anywhere in the tree. | `codon_genie/codon_selector.py`; full git tree |
| `barcode_generation` | **no** | No barcode, UMI, GC-constraint or edit-distance code in **any of the five implementation modules** (`client.py`, `codon_selector.py`, `codon_utils.py`, `ncbi_tax_utils.py`, `seq_utils.py`) or in any of the static JS. *(Corrected: the extraction said "4 modules".)* | full git tree; all five implementation modules read raw |
| `per_sequence_provenance` | **partial** | The exposed tool emits no output *sequences*, so under a strict reading this is `no`. `partial` is the deliberate generous call: each codon record genuinely explains its own derivation — `ambiguous_codon_nucleotides`, the full `ambiguous_codon_expansion`, per-amino-acid `codons` each with `cai` and `probability`, `type` (required / off-target / stop), and `score`. Verified live. The record explains *why* the design was chosen; it just annotates a codon, not a library member. **Flagged as convention-dependent** (see §5). | live JSON; `codon_genie/codon_selector.py` |
| `automatic_naming` | **no** | Output is an unnamed sorted JSON array keyed on the ambiguous codon string. No FASTA/ID/label emitter in `main.py`, `codon_selector.py`, or the UI. (The IUPAC codon string is a systematic set encoding, not a generated design name.) | `main.py`; `codon_genie/codon_selector.py`; `static/result/result.html` |
| `design_visualization` | **partial** | The web UI renders a colour-coded interactive result table, per Fig. 1 caption (p. 5): "Variant codons are shown in grey, with their encodings shown in **green, orange and red for required amino acids, additional amino acids and stop codons**, respectively... these encodings and their codon usage frequencies may be visualised through a **tooltip**." Verified in `static/result/result.html`: `label-success` / `label-warning` / `label-danger` for `type` 1 / 0 / -1, plus two `uib-tooltip`s (per-codon "codon (probability)" and the nucleotide-set string). The 2016-era CDN dependencies (bootstrapcdn, googleapis angular, typeahead.js, angular-ui) all still return 200, so the live UI genuinely renders. It visualises **one codon**, not a library or a design graph; a strict reading ("an HTML table is not a visualisation") gives `no`. `partial` is the generous-safe call. | PDF p. 5 Fig. 1; `static/result/result.html`, `result.directive.js`, `result.css`; live https://codongenie.appspot.com |

### Block B — assay coverage

| key | value | evidence | source |
|---|---|---|---|
| `assay_dms` | **partial** | Declared purpose is squarely coding-sequence variant libraries: "designing ambiguous codons to support protein mutagenesis applications" (Abstract); motivated by directed evolution and by "the study of sequence-to-fitness relationships (Hietpas, Jensen & Bolon, 2011)" (p. 1) — DMS-adjacent. But it designs a codon, not a DMS library: it enumerates only the DNA variants of one codon, emits no sequences from the design API, and has no notion of scanning every position of an ORF. `partial` = solves one sub-step of coding-library design; a stricter matrix definition ("must produce a library") gives `no`. **Flagged as convention-dependent** (§5). | PDF Abstract, p. 1 |
| `assay_mpra` | **no** | Entirely amino-acid/codon-centric. `static/index.html` renders exactly 20 amino-acid buttons (non-polar / polar / acidic / basic); the Analyse input is regex-limited to `[acgtmrwsykvhdbnACGTMRWSYKVHDBN]{3}` (the literal is in `static/codonGenie/codonGenie.ctrl.js:5`, bound in `index.html` as `data-ng-pattern="ctrl.codon_pattern"`). No promoter, enhancer, UTR, or TFBS concept in paper or repo. | `static/index.html`; `static/codonGenie/codonGenie.ctrl.js`; PDF full text |
| `assay_insilico` | **no** | The only scoring is the deterministic codon-usage formula (`_get_score` over `cai` values). No model, prediction, or ML dependency; `requirements.txt` is two lines: `Flask`, `Biopython`. Published 2017. | `codon_genie/codon_selector.py`; `requirements.txt` |

### Block C — genomics integration

| key | value | evidence | source |
|---|---|---|---|
| `genome_coordinates` | **no** | The sole "genomic" input is an NCBI **Taxonomy** id, used to fetch a codon usage table (`codon_utils.CodonOptimiser.__get_codon_usage` scrapes `http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=<taxid>&aa=1&style=GCG`). That is an organism selector, not a coordinate system. No chromosome/start/end/strand parameter in any of the four routes. | `codon_genie/codon_utils.py`; `main.py` |
| `transcript_models` | **no** | No GTF/GFF/BED reader in the tree. The only data file is `codon_genie/species.table` (name → taxid). Note `ncbi_tax_utils.py` does parse NCBI taxdump `nodes.dmp`/`names.dmp`, but that is taxonomy, not transcript annotation. | full git tree; `codon_genie/ncbi_tax_utils.py` |
| `exon_intron_split_codons` | **no** | Codons are isolated triplets: `CODONS` maps each amino acid to `[pos1, pos2, wobble-set]`, positions 1–2 are aligned and wobble sets collapsed by `_optimise_pos_3`. No exon, intron, splice-site, or reading-frame-phase concept. | `codon_genie/codon_selector.py` |
| `hgvs_input` | **no** | `request.args['aminoAcids']` is a bare one-letter amino-acid string; `request.args['codon']` is a 3-character IUPAC code; `organism` is a taxid. No variant-nomenclature parser exists. | `main.py`; `codon_genie/codon_selector.py` |
| `vcf_vep_output` | **no** | "In all cases, results are returned in **json** format" (p. 3). Confirmed by exhaustion: every data route (`/codons`, `/organisms/`, `/organisms/<term>`) returns `Response(json.dumps(...), mimetype='application/json')` while `/` serves the static UI, and `result.html` / `result.directive.js` contain no export, download, copy, or CSV control. No VCF/FASTA/GenBank/CSV writer anywhere. | PDF p. 3; `main.py`; `static/result/*` |
| `consequence_annotation` | **partial** | The tool does classify what a design will produce: `_get_amino_acid_type` returns `-1` for Stop, `0` for an additional (off-target) amino acid, `1` for a required amino acid, surfaced in both the JSON and the colour-coded UI. Verified live, not just in source: `?aminoAcids=WY&organism=37762` → `TRK` with `['W:1','Y:1','C:0','Stop:-1']`, so stop-codon and off-target flagging is real and user-visible. It is **not** variant-consequence annotation against a reference transcript (no synonymous / missense / frameshift / in-frame-indel calling). Drop to `no` if the row is defined strictly as VEP-style annotation. **Flagged as convention-dependent** (§5). | live API; `codon_genie/codon_selector.py`; `static/result/result.html` |

### Block D — physical construction

| key | value | evidence | source |
|---|---|---|---|
| `primer_design` | **no** | **No Tm, annealing, overlap, oligo, or primer routine exists anywhere in the repository.** `seq_utils.py` is 66 lines containing exactly two functions (`find_invalid`, `_get_restr_type`) plus five reagent-concentration constants (`NA`, `K`, `TRIS`, `MG`, `DNTP` and `__DEFAULT_REAG_CONC`) that are **never referenced anywhere in the repo** — vestigial from `liv-utils`. *(Corrected: the extraction's wording implied a Tm routine existed somewhere; it does not.)* The paper defers oligo design to the authors' other tool: "CodonGenie is amenable to future integration with new and existing variant library design software tools (Swainston et al., 2014)" — i.e. GeneGenie. | `codon_genie/seq_utils.py` (read in full); repo-wide grep; PDF p. 7 |
| `codon_optimization` | **yes** | Two independent grounds. **(1) User-facing:** organism-specific codon usage is the headline ranking axis — codons are "ranked according to their efficiency in encoding the required amino acids while minimising the inclusion of additional amino acids and stop codons. **Organism-specific codon usage is also considered**" (Abstract). Reproduced live: taxid 37762 (*E. coli*) → DTK 0.878 > DTS 0.684; taxid 1902 (*S. coelicolor*) → DTS 0.778 > DTK 0.294 — the host-dependence reversal of Table 1, p. 4. For invariant positions the tool returns the plain best codon: "The codon returned is the therefore [sic] the most frequent codon for encoding proline in E. coli" (p. 7). **(2) Package-level:** `codon_utils.CodonOptimiser` provides a full sequence optimiser — `get_codon_optim_seq(protein_seq, excl_codons, max_repeat_nuc, restr_enzyms, ...)`, `get_cai`, `get_best_codon`, `get_random_codon`, `get_all_codons`, `mutate` — and `codon_utils.py` is confirmed present inside `CodonGenie-1.7-py3-none-any.whl`. **Caveat:** (2) is vendored from `liv-utils` and is *not* exposed by the web UI or the REST API; the `yes` stands on (1) alone. | PDF Abstract, p. 4 Table 1, p. 7; live API (taxids 37762 and 1902); `codon_genie/codon_utils.py`; PyPI wheel contents |
| `synthesis_constraints` | **partial** | `seq_utils.find_invalid(seq, max_repeat_nuc, restr_enzyms)` performs a homopolymer-run regex check and a Biopython `Restriction` site search, and `get_codon_optim_seq` backtracks on its result (it re-draws synonymous codons stochastically, and with `tolerant=True` progressively permits otherwise invalid patterns once `max_attempts` is exhausted, so the constraint is not a hard guarantee). Real, shipped code — both modules confirmed present in the 1.7 wheel. **But there is no reachable call path from `main.py` or any REST route**, and it is never applied to ambiguous-codon design. What the user-facing tool *does* report is exact library size (`len(ambiguous_codon_expansion)`, e.g. DTK's six variants vs NTN's 16, p. 4) — a screening/throughput constraint, not a per-sequence synthesis constraint. **Convention-dependent and asymmetric:** unlike `codon_optimization`, this cell rests *entirely* on unexposed package code. Under a strict user-facing convention it is `no`; under a package-level convention arguably `yes`. `partial` is correct only under the convention written down in §5. | `codon_genie/seq_utils.py`; `codon_genie/codon_utils.py`; PyPI wheel contents; PDF p. 4 |

### Block E — engineering

| key | value | evidence | source |
|---|---|---|---|
| `interface` | **yes** | Four surfaces, **no CLI** (`setup.py` has no `entry_points`/`console_scripts`; `main(argv)` takes only a Flask port; no argparse anywhere): (a) single-page **web GUI** (AngularJS 1.5 + Bootstrap 3.3, `static/index.html`) — live at https://codongenie.appspot.com, HTTP 200, all 11 external CDN assets (across five hosts) still resolve; (b) **RESTful JSON web service** — `/codons?aminoAcids=&organism=`, `/codons?codon=&organism=`, `/organisms/`, `/organisms/<term>` (p. 3); (c) **Python client library** `CodonGenieClient` (`codon_genie/client.py`, present in the wheel) with `get_codons`, `analyse`, `search_organisms`, `get_organisms`, and a `url=` constructor argument that can retarget it at a self-hosted instance (default: the live appspot deployment); (d) **Docker** (`Dockerfile` + `start_server.sh`) and Google App Engine (`app.yaml`). **Caveat added on review:** PyPI declares **no** dependencies (`info.requires_dist` is `null`; `setup.py` has no `install_requires`), and neither `species.table` nor `main.py` is present in the 1.7 wheel or sdist (verified by unzipping/untarring both). So `pip install CodonGenie` alone cannot run the app or list organisms — the working paths are `git clone` + Docker (documented, not build-tested here), or the live REST API. | `main.py`; `setup.py`; `static/index.html`; `codon_genie/client.py`; `Dockerfile`; `app.yaml`; PyPI JSON + downloaded wheel/sdist; live HTTP probes |
| `license` | **yes** | **MIT**, confirmed three ways: `LICENSE` file in the repo tree; GitHub reports `license = MIT` (`spdx_id: MIT`); `setup.py`/`main.py` headers state "CodonGenie is licensed under the MIT License"; PyPI classifier `License :: OSI Approved :: MIT License`. (The wheel `METADATA` literally reads `License: UNKNOWN` because `setup.py` omits the `license=` kwarg — packaging noise, not a real ambiguity.) Paper itself is CC-BY 4.0. | repo `LICENSE`; GitHub repo page; `setup.py`; PyPI JSON; wheel METADATA |
| `maintained` | **partial** ⚠️ *(changed from `yes`)* | **Public, installable and one live deployment — but no development activity in ~4 years and no GitHub releases or tags ever.** Verified: last commit of any kind **2022-12-12** — the author's own merge of an external `TrellixVulnTeam` tarfile-sanitisation PR (the merge commit is authored by Neil Swainston; the patch commit is the bot's), not development work; last commit carrying the author's own code changes **2022-06-09** (4 yr 2 mo before survey); last PyPI release **1.7, 2022-03-29** (4 yr 4 mo); **zero** GitHub releases and **zero** tags, ever; 86 commits, 1 star, 3 forks, 0 open issues; not archived. The URL published in the paper (`http://codon.synbiochem.co.uk`) **404s**, and SYNBIOCHEM's own tool directory still advertises the dead link. The corresponding author's paper contact (`neil.swainston@manchester.ac.uk`) no longer matches `setup.py` (`neil.swainston@liverpool.ac.uk`). The tool still *works* — that is **availability**, which `availability_status` carries, not maintenance. Scored under the explicit convention in §5. | GitHub repository/commit metadata (re-verified 2026-08-14); PyPI JSON release dates; live 404 on `codon.synbiochem.co.uk`; `setup.py` |

---

## 4. Additional capabilities NOT covered by the current row list

1. **Ranked trade-off curve rather than a single answer.** "CodonGenie returns not only the most 'specific' ambiguous codons... Providing results that include less specific ambiguous codons, which may also encode additional amino acids, allows the user to perform a **trade-off between library size and codon specificity**" (pp. 2–3). Results are sorted by `(len(expansion) ascending, -score)` (`_format_results`).
2. **Exact library-size arithmetic per degenerate codon** — `ambiguous_codon_expansion` gives the precise DNA-variant count, which is the screening-burden currency in directed evolution.
3. **Reverse / "Analyse" direction** — decode an existing degenerate codon into encoded amino acids + host usage frequencies. A validator/audit mode.
4. **Off-target and stop-codon flagging as a first-class output field** (`type` ∈ {1, 0, −1}), surfaced in JSON and colour-coded in the UI.
5. **Per-amino-acid redundancy quantified and surfaced (ADDED on review).** Each output amino-acid record carries the *list* of expansion codons that encode it, and `result.html` renders the **count** next to the letter (`{{aa.amino_acid}} {{aa.codons.length}}`). This is an explicit statement of **encoding redundancy within the designed pool** — how many of the codon's DNA variants give each amino acid. (Each codon's `probability` is the host's codon-usage frequency, not a synthesis-mixture fraction; the interface accepts and reports no non-uniform nucleotide-mixture fractions.) Verified live: `?aminoAcids=AFGILMV&organism=37762` → `DBK`, 18 variants, with A×2, G×2, V×2, S×3, T×2 — i.e. valine is over-represented 2:1 relative to singly-encoded residues.
6. **Organism resolution service** — free-text organism search → NCBI Taxonomy ID over **35,792** organisms (re-verified live), backed by the Kazusa Codon Usage Database.
7. **NCBI Taxonomy tree descent and synonym expansion (ADDED on review).** `ncbi_tax_utils.TaxonomyFactory` downloads `ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`, parses `nodes.dmp`/`names.dmp`, and exposes `get_child_ids(parent_id)` (recursive descent of the taxonomy tree) and `get_names(tax_id)` (all synonyms); `get_codon_usage_organisms(expand=True)` uses it to expand the Kazusa species table with every child taxon and every name synonym. So organism resolution is not a flat name→taxid lookup — it can walk NCBI Taxonomy. Ships in the pip wheel; `expand` defaults to `False` and it is **not** exposed via the REST API.
8. **Alternative genetic codes** — `CodonSelector(table_id=1)` accepts any NCBI codon table id via `Bio.Data.CodonTable.unambiguous_dna_by_id`; not exposed in the API but present. Only the codon→amino-acid lookup follows the chosen table: candidate generation still uses the hard-coded standard-code `CODONS` map, so a non-standard table is not fully supported. The paper flags "augmented genetic codes and expanded genetic alphabets" as future work (p. 7).
9. **CAI-style scoring of an arbitrary DNA sequence** (`get_cai` — an *arithmetic* mean of per-codon relative adaptiveness that silently skips unrecognised triplets, not the conventional geometric-mean CAI) and stochastic codon sampling (`get_random_codon`, `mutate(protein_seq, dna_seq, mutation_rate)` — one sequence per call, with no check that `dna_seq` encodes `protein_seq`) in the shipped utils.
10. **Deployability as a microservice** — Docker image + `app.yaml` for Google App Engine, explicitly so it can be embedded in larger design pipelines (p. 7).

**ROW-LIST SUGGESTION (retained from the extraction):** **degenerate/ambiguous codon design** is a distinctive CodonGenie capability, supported by the paper and upstream `codon_selector.py`. CodonGenie has no Table 2 column, so this record does not propose or ship a Table 2 row; the non-actionable substitution candidate is reported in the Pass 2 log.

---

## 5. Scoring conventions this record relies on (must be applied uniformly across all tools)

These are written down because three cells hinge on them. If the orchestrator adopts different thresholds, the affected cells must be recomputed for **every** tool, not just this one.

- **C1 — `maintained` threshold.** Adopted here: `yes` = an upstream commit **or** package release within ~2 years of the survey date; `partial` = 2–5 years since the more recent of those events while the official repository and package remain available; `no` = >5 years since both, archived, or unavailable. Under C1, CodonGenie (latest upstream commit 2022-12-12; latest PyPI release 2022-03-29) = **partial**. Here `maintained` measures upstream project activity, not whether a deployment remains available.
- **C2 — "ships in the package" vs "user-facing".** Adopted here: a capability present in the installable package but unreachable from the tool's own UI/API scores **`partial`**, not `yes`. Consequences: `synthesis_constraints` = partial (rests *entirely* on the unexposed `seq_utils.find_invalid`); `codon_optimization` = yes anyway, because it has independent user-facing grounds (organism-specific ranking) and does not need the package argument.
- **C3 — generous reading on borderline scope cells.** Where a tool performs a genuine, user-visible sub-step of a row's capability but not the full capability, score **`partial`** with an explicit scope note rather than `no`. Affects `per_sequence_provenance`, `design_visualization`, `assay_dms`, `consequence_annotation`. Under a strict reading all four would be `no`; the extractor and the reviewer independently agreed `partial` is the referee-safe posture. Applying C3 to CodonGenie but not to other tools would be inconsistent.

---

## 6. Documented examples / vignettes (candidates for PoolParty reproduction)

1. **MSA → degenerate synthetic gene** (flagship example, paper pp. 7–8 + `codon_genie/example/align.py`). Alignment `['PFDMR','PIAMR','PLHLR','PMNMR','PVHMR']`, host *E. coli* → **`CCG|DTK|VMT|MTG|CGT`**. Runnable script: https://github.com/neilswainston/CodonGenie/blob/master/codon_genie/example/align.py — the only example producing something resembling a library specification. Note the per-position `for` loop is user code, not tool functionality: it drops gap characters, takes only the first-ranked codon at each column, and its organism lookup returns the first name-search hit.
2. **Non-polar set F, I, L, M, V** (p. 4). Naive `NTN` = 16 DNA variants; CodonGenie finds **`DTK`** = [AGT][T][GT] = **6** variants. Scaling claim: "when using 3 DTK codons the library size is reduced from 4,096 (16³) to 213 (6³)". **The paper prints 213; 6³ = 216.** *(Corrected: the extraction called 213 a PDF typesetting corruption. I re-extracted the page — both superscripts survive as trailing digits ("163", "63"), so 213 is the literal published value, i.e. an arithmetic typo in the paper. Use 216; never quote 213; and do not characterise it as a rendering artefact to a referee who co-authored the paper.)*
3. **Host-dependence of the same set** (p. 4, Table 1) — **reproduction recipe corrected**. Use **taxid 37762** (*E. coli*, the id used in the repo's own unit test): `?aminoAcids=FILMV&organism=37762` → DTK **0.878**, DTS **0.684** — exactly the paper's 0.88 / 0.68. And **taxid 1902** (*S. coelicolor*): DTS **0.778**, DTK **0.294** — the same reversal as the paper, with DTK matching its 0.29 but DTS rounding to 0.78 against the paper's 0.79. *(The extraction used taxid 83333 (E. coli K-12), which returns DTK 0.82 and does not reproduce the paper's 0.88.)* Driven by wobble-position T vs C preference (Phe TTT/TTC = 0.64/0.36 in *E. coli* vs 0.03/0.97 in *S. coelicolor*).
4. **Fig. 1 GUI example** (p. 5). Non-polar residues **A, F, G, I, L, M, V** → preferred codon **`DBK`** = [AGT][CGT][GT] = 18 DNA variants. Verified live at taxid 37762: DBK first, 18-codon expansion, with A, G, V each encoded twice and S three times (see additional capability 5).
5. **Analyse mode** (p. 3). `GET /codons?codon=NSS&organism=4932` (*S. cerevisiae*) — verified live: 16 codons in the expansion and **8** encoded amino acids (A, C, G, P, R, S, T, W), all typed `0` (off-target, since no target set was given). *(Corrected: the extraction memo said 9 amino acids.)*
6. **Unit tests — two modules, not one** *(corrected)*.
   - `codon_genie/test/test_codon_genie.py`: `CodonSelector().optimise_codons('FLIMV', '37762')` → asserts 12 results, best = `DTK`. (Verified live: 12 results, DTK first.)
   - `codon_genie/test/test_client.py` (**added on review**): five live-service integration fixtures, one of which is **the paper's own documented API example** — `client.get_codons('DE', 4932)` → asserts 4 results with best = **`GAW`**, i.e. `?aminoAcids=DE&organism=4932` from PDF p. 3 (verified live: 4 results, GAW first). Also `client.analyse('NTT', '4932')` → 4 amino acids; `search_organisms('escherich')`; `search_organisms('escherichia co')` (space handling); `get_organisms()` non-empty.
7. **Contrast case vs DYNAMCC** (p. 5). For F/I/L/M/V, DYNAMCC returns the *set* {`WTT`, `VTG`} = 5 variants with no off-target encoding; CodonGenie returns a *single* codon DTK/DTS = 6 variants. This is CodonGenie's own statement of its scope boundary.

---

## 7. Third-party fork — record for pre-emption (changes no cell)

The paper and upstream repository are sufficient to establish the tool's scope. Upstream `master` (`c26439f`, re-verified 2026-08-14) has no full-sequence or multi-position design method: `CodonSelector.optimise_codons` accepts one amino-acid set, while the paper's multi-position workflow concatenates per-position results in `example/align.py`.

Two consequences:

1. **The official evidence supports `library_as_object = no` and `mixed_mutagenesis_one_pool = no`.** No third-party corroboration is needed.
2. **No fork, sibling record, or deployment inference is used as evidence.** Capability conclusions rest on the paper, upstream source and package index.

---

## 8. Stated and material limitations

**The tool's own words:**

- Scope is deliberately one codon: "CodonGenie designs **single** ambiguous codons to encode a desired set of amino acids, **which may also include off-target amino acids**" (p. 5) — contrasted with DYNAMCC, which "designs **sets** of ambiguous codons... with no off-target amino acid encoding and minimal redundancy."
- The efficiency cost of that choice is acknowledged: "five DNA variants encode the five desired amino acids [DYNAMCC], while CodonGenie's solution of DTK or DTS encode six DNA variants, thus producing a **larger library**" (pp. 5–6).
- Library assembly is out of scope and deferred to sister tools: "CodonGenie is amenable to **future integration** with new and existing variant library design software tools (Swainston et al., 2014)" (p. 7).
- Minimalism is a design principle, not an oversight: "∼700 lines of code"; "requires no documentation" (p. 6).
- Only 17 of 20 amino acids are cleanly handled by the fast algorithm; L, R and S require the algorithm to be run per-encoding and the results combined (p. 2).

**Unstated but material (verified, added on review):**

- **Live plain-HTTP screen-scrape dependency.** `CodonOptimiser.__get_codon_usage` fetches `http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=<taxid>&aa=1&style=GCG` at request time and parses the `<PRE>` block. If Kazusa changes its page format or goes offline, a self-hosted instance loses all codon-usage scoring — not silently: `urlopen` and the `<PRE>` parse are unguarded, so failures propagate to the Flask exception handler and surface as an HTTP 500.
- **Only Kazusa-listed taxids work, and failures are raw 500s.** `?aminoAcids=FILMV&organism=562` returns **HTTP 500** with `{"message": "KeyError: 'ATC'"}`, because Kazusa has no table for that id, `__get_codon_usage` silently builds an empty table, and the lookup explodes downstream. No input validation, no useful error message, no fallback. (Re-verified 2026-08-11.)
- **The pip path is broken for the application.** The 1.7 wheel and sdist contain the five `codon_genie` implementation modules and both test modules but **no `species.table` and no `main.py`** — so a pip install cannot serve the app or run `get_codon_usage_organisms()`. `CodonSelector.optimise_codons()` works from a pip install only if the user separately installs Biopython, since PyPI declares no dependencies at all (`requires_dist: null`).
- **The paper's advertised entry point is dead** (see §9).

---

## 9. Availability TODAY (verified 2026-08-10/11; URL and release status re-verified 2026-08-14)

| Asset | Status |
|---|---|
| **Paper's advertised URL** `http://codon.synbiochem.co.uk` | **DEAD.** HTTP **404** (Google App Engine `Error 404 (Not Found)!!1`); HTTPS fails cert-name validation (the served certificate is Google's, CN `*.appspot-preview.com`, with `*.appspot.com` among its many SANs and nothing matching this host), so the custom domain mapping has lapsed. Both `/` and the paper's documented example `?aminoAcids=DE&organism=4932` 404. SYNBIOCHEM's own tool directory (synbiochem.co.uk/pipeline/design/) still advertises this dead link. |
| **Working deployment** `https://codongenie.appspot.com` | **LIVE, HTTP 200.** Full AngularJS UI served (all 11 external 2016-era CDN assets, across five hosts, still return 200); REST API functional and paper-consistent — FILMV @ 37762 → DTK 0.878 top-ranked; `/organisms/` returns exactly 35,792 entries. This URL is hard-coded as the default in `client.py` but appears **nowhere in the paper**, so a reader following the publication lands on a dead link. **Cite this URL (or the repo), never the paper's.** |
| **Repo** | `github.com/synbiochem/CodonGenie` returns **301** to **`github.com/neilswainston/CodonGenie`**, whose response is 200. Public, not archived, MIT. 86 commits, 1 star, 3 forks, 0 open issues, **0 GitHub releases, 0 tags**. |
| **Last activity** | Last commit **2022-12-12** — the author's merge of the `TrellixVulnTeam` bot's PR; last commit carrying the author's own code changes **2022-06-09** ("Minor auto-reformatting."; the preceding same-day commit is "Update to remove dependency on remote species.table flatfile."). ~3 yr 8 mo since any upstream commit; ~4 yr 2 mo since the author's own code changes. |
| **PyPI** `CodonGenie` | **v1.7 (2022-03-29)**, 8 releases since 2020-01-17, `python_requires >=3.7`, **no declared dependencies** (`requires_dist: null`). Wheel/sdist omit `species.table` and `main.py`, so pip alone cannot run the app. |
| **Self-hosting** | `Dockerfile` + `start_server.sh` + `app.yaml`. This (or the live REST API) is the documented reproduction path — **not** `pip install`. The container was not built or run here (see §10). Subject to the Kazusa screen-scrape dependency above. |

**Verdict:** functionally alive — verified live on the App Engine instance, with a documented but here-untested clone+Docker path — but **effectively unmaintained** (no development commit in ~4 years, no releases or tags ever) and the **published entry point is broken**.

---

## 10. Unresolved disagreements and least-confident cells

**No genuine unresolvable extractor-vs-reviewer disagreements remain.** The reviewer would change exactly one of the 24 verdicts (`maintained`) and let the other 23 stand, and I re-derived the disputed facts independently; every relevant factual correction was confirmed (taxid 37762 vs 83333, five implementation modules, 8 not 9 amino acids for NSS, `requires_dist: null`, missing `species.table`, the literal "213" in the PDF, two test modules, and the HTTP 500 on taxid 562). No cell is set to `unknown`.

The one changed cell, `maintained` (yes → partial), is not a factual dispute — extractor and reviewer agree on every underlying fact — but a **threshold** dispute, resolved by writing down convention C1 in §5. Flagged for the orchestrator: if C1 is changed, this cell must move with it (a strict "activity within N years" rule gives `no`).

**Least-confident / convention-dependent cells, in order of exposure:**

1. `maintained` (**partial**) — depends entirely on convention C1. `no` under a strict activity threshold; `yes` only if "maintained" is read as "still available", which duplicates `availability_status`.
2. `synthesis_constraints` (**partial**) — depends entirely on convention C2. The capability is real code shipped in the wheel but with **no reachable call path** from the web app or REST API. `no` under a strict user-facing convention.
3. `per_sequence_provenance` (**partial**) — convention C3. Rich structured metadata exists, but the exposed design operation emits no output sequences for it to attach to. `no` under a strict reading; `partial` chosen so a referee cannot say the JSON output was ignored.
4. `consequence_annotation` (**partial**) — convention C3. The `type` field (stop / off-target / required) is genuine, user-visible consequence flagging at the codon level, verified live; it is not VEP-style annotation against a reference transcript.
5. `assay_dms` (**partial**) — convention C3. Site-saturation mutagenesis of coding sequences is unambiguously the declared purpose, but the design API emits no full-sequence variants and has no ORF-scanning concept. `no` if the row requires producing a library.
6. `design_visualization` (**partial**) — convention C3. The colour-coded table and tooltips are real, render live, and are design output — but they visualise one codon, not a library or a design graph.
7. `codon_optimization` (**yes**) — highest confidence of the judgement calls; doubly grounded (user-facing ranking axis reproduced live, plus a full sequence optimiser in the package) and survives even if convention C2 is tightened.

**Nothing in this record was verified by installing or executing the tool**, per the survey constraints. All statements come from the PDF, the `master` source read via `raw.githubusercontent.com`, the GitHub commit atom feeds, the PyPI JSON and the downloaded (not installed) wheel and sdist, and live HTTP probes of both deployments.

---

## 11. Corrections to the prior analysis

1. **Repo URL is stale.** The paper (and the prior note) cite `synbiochem/CodonGenie`; the repo is now `neilswainston/CodonGenie`. Cite the redirect target.
2. **The advertised web URL is dead.** `codon.synbiochem.co.uk` 404s. If the manuscript or rebuttal names a CodonGenie URL it must be `https://codongenie.appspot.com` or the GitHub repo.
3. **"Written in Python/Flask, ∼700 lines"** — the prior note's own wording; both facts confirmed, but as two separate passages: "written in Python (using the Flask framework) and HTML/Javascript" (p. 4) and "the entire application consists of ∼700 lines of code" (p. 6). The combined phrasing appears nowhere in the paper.
4. **"Operates at the single-codon level, not at the level of full sequence libraries"** — confirmed and well-supported. The strongest citable line is the paper's own microservice statement (p. 6), backed by the fact that its multi-position example is a user-written `for` loop in `example/align.py` (§7).
5. **The prior note overlooks the shipped `CodonOptimiser` / `find_invalid` utilities.** These drive the `codon_optimization` and `synthesis_constraints` cells and an author-referee would raise them. Handled above with the explicit "present in package, not exposed in the tool" convention (C2).
6. **Ambiguous-codon output.** CodonGenie's Design result is a ranked list of IUPAC ambiguous codons and their concrete expansions; comparisons with other tools require their own primary evidence.
7. **Do not repeat the extraction's claim that "213" is a typesetting corruption.** The paper prints 213; 6³ = 216. State it that way (§6 item 2).
