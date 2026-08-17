# Adversarial review — CodonGenie extraction

**Reviewed:** `extractions/codongenie.md` + the structured extraction JSON
**Review date:** 2026-08-10
**Reviewer posture:** attempt to falsify every cell; assume a referee may be Swainston/Currin/Kell.

---

## 0. Bottom line

The extraction is **unusually good**. I independently re-derived essentially every claim from
primary sources (raw GitHub source, GitHub REST API, PyPI JSON + the actual wheel/sdist, live HTTP
probes of both deployments, the local PDF) and found **one capability verdict I would change
(`maintained`)**, plus a handful of evidence-level inaccuracies that a tool-author referee could
use to chip at credibility, and **five things the extractor missed** — one of which
(a 2025 third-party fork that adds full-sequence design) is the kind of thing a referee would
enjoy pointing at.

No cell is *understated* in the dangerous direction. The extractor consistently chose the generous
reading (`partial` rather than `no`) on the four borderline cells, which is the correct
referee-safe posture. Nothing is inflated in CodonGenie's disfavour.

---

## 1. What I verified independently (not by re-reading the extractor's quotes)

| Check | Method | Result |
|---|---|---|
| Repo identity, license, activity | `api.github.com/repos/neilswainston/CodonGenie` | MIT; `pushed_at` 2022-12-12T11:42:23Z; created 2016-10-07; 1 star; 3 forks; 0 open issues; not archived; **0 releases, 0 tags** (both endpoints return `[]`) — all as claimed |
| Commit log | `/commits?per_page=20` | Last = 2022-12-12 `Merge pull request #1 from TrellixVulnTeam/master`; previous = 2022-11-28 bot commit; **last human commit 2022-06-09** — as claimed |
| Full file tree | `/git/trees/master?recursive=1` | As claimed **except**: `codon_genie/` holds **5** Python modules, not 4 — the extractor never lists `ncbi_tax_utils.py`, and never read `test/test_client.py` |
| `main.py` routes | raw source | Exactly 4: `/`, `/organisms/`, `/organisms/<term>`, `/codons`. Every response is `json.dumps(...)`, `mimetype='application/json'`. Confirmed |
| No CLI | `setup.py` | No `entry_points`, no `console_scripts`, no argparse anywhere. `main(argv)` only takes a port for the Flask server. **"no CLI" confirmed** |
| `codon_selector.py` | raw source | `CODONS` dict, `_get_amino_acid_type` (-1/0/1), `_get_score`, `_format_results` sort key `(len(expansion), -score)` — all exactly as described |
| `seq_utils.py` | raw source | 70 lines: `find_invalid` (homopolymer regex + `Bio.Restriction` search) and `_get_restr_type`. **There is no Tm function at all** (see §3.1) |
| `codon_utils.py` | raw source | `CodonOptimiser.get_codon_optim_seq`, `get_cai`, `get_best_codon`, `get_random_codon`, `get_all_codons`, `mutate`; `__get_codon_usage` screen-scrapes `http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=<taxid>&aa=1&style=GCG` parsing the `<PRE>` block. Confirmed verbatim |
| Paper URL dead | `curl http://codon.synbiochem.co.uk/` | **HTTP 404**, Google `Error 404 (Not Found)!!1`, from 142.250.191.20. HTTPS fails cert-name validation (`no alternative certificate subject name matches`). Confirmed |
| SYNBIOCHEM's own directory | synbiochem.co.uk/pipeline/design/ | Still advertises the dead `http://codon.synbiochem.co.uk`. Corroborates |
| Live deployment | `https://codongenie.appspot.com` | HTTP 200. `/codons?aminoAcids=FILMV&organism=83333` → 12 results, DTK first, score 0.8152. `/organisms/` → **exactly 35,792** entries (3.0 MB). `/codons?codon=NSS&organism=4932` → 16-codon expansion. All confirmed |
| Stop-codon flagging live | `?aminoAcids=WY&organism=37762` | `TRK` → `['W:1','Y:1','C:0','Stop:-1']`. `type=-1` confirmed live, not just in source |
| UI colour coding | `static/result/result.html` | `label-success` (green) for `type==1`, `label-warning` (orange) for `0`, `label-danger` (red) for `-1`; two `uib-tooltip`s (codon+probability; nucleotide sets). Confirmed. **No export/download control anywhere in the UI** |
| UI CDN deps still resolve | HEAD on all 5 externals | bootstrapcdn 200, googleapis angular 200, typeahead.js 200, angular-ui 200. The 2016-era UI still loads. (Good news the extractor didn't check but should have) |
| PyPI | `pypi.org/pypi/CodonGenie/json` + downloaded wheel + sdist | v1.7, 8 releases, latest 2022-03-29, `requires_python >=3.7`. **See §3.3 — `requires_dist` is `null`, and `species.table` is NOT in the wheel or sdist** |

---

## 2. The one verdict I would change

### `maintained` = **yes** → **OVERSTATED**

This is the only cell where the extractor's own evidence contradicts the value it assigns. The
evidence string ends "**Effectively dormant for ~3.5 years as of 2026-08-10**" and the
availability field says "**effectively unmaintained**" — and the cell still reads `yes`.

Verified facts:

- Last commit of **any** kind: 2022-12-12 — and it is an automated `TrellixVulnTeam` tarfile
  sanitisation bot PR, not maintenance by the authors.
- Last commit by a human author: **2022-06-09** — i.e. **4 years 2 months** ago, not 3.5.
- **Zero** GitHub releases, **zero** tags, ever. 1 star. 0 open issues (and 0 closed issues to
  speak of — nobody is filing any).
- Last PyPI release 2022-03-29 (**4 years 4 months**).
- The URL published in the paper, and still advertised on SYNBIOCHEM's own tool directory,
  **404s** — which is itself the clearest possible evidence that nobody is minding the tool.
- Corresponding author's affiliation has moved twice since publication (Manchester → Liverpool →
  "GeneGenie Bioinformatics Ltd." per the source headers), and the paper's contact address
  (`neil.swainston@manchester.ac.uk`) no longer matches `setup.py`
  (`neil.swainston@liverpool.ac.uk`).

The only argument for `yes` is that the artefact still *works*: repo public and not archived,
PyPI installable, one deployment live. That is **availability**, not maintenance, and the matrix
already has an `availability_status` free-text field carrying it.

**Suggested value: `partial`** (public, installable, one live instance — but no development
activity in ~4 years, no releases/tags ever, and the published entry point is broken).
Under a strict "active development within N years" convention the honest value is **`no`**.

⚠️ **This must be resolved as a cross-tool convention, not per-tool.** If `maintained` is scored
generously here, the same generosity has to be extended to every other surveyed tool whose last
commit is 2022 — otherwise the matrix is internally inconsistent and a referee will find it.
Note the direction of risk: an inflated `yes` here is the *safe* error w.r.t. an author-referee,
but the *unsafe* error w.r.t. a referee who checks the matrix for even-handedness.

---

## 3. Evidence-level errors that do not change a verdict (but should be fixed before this text goes near a referee)

### 3.1 `primer_design` evidence overstates what's in `seq_utils.py`
The extraction says "*`seq_utils.py` holds reagent-concentration constants but `main.py` and
`codon_selector.py` never invoke any Tm routine*". This implies a Tm routine exists somewhere. It
does not. `seq_utils.py` is 70 lines containing exactly two functions (`find_invalid`,
`_get_restr_type`) plus five dead constants (`NA/K/TRIS/MG/DNTP` and `__DEFAULT_REAG_CONC`) that
are **never referenced anywhere in the repository**. The verdict `no` is *stronger* than the
evidence sentence suggests. Reword to: "no Tm, annealing, overlap or primer routine exists
anywhere in the repository; `seq_utils.py` retains five unused reagent-concentration constants
vestigial from `liv-utils`."

### 3.2 `barcode_generation` evidence miscounts the tree
"*4 Python modules in `codon_genie/`*" — there are **5**: `client.py`, `codon_selector.py`,
`codon_utils.py`, `ncbi_tax_utils.py`, `seq_utils.py`. Verdict unaffected (none of them contain
barcode/UMI/edit-distance code — I checked all five), but the miscount reveals that
`ncbi_tax_utils.py` was never opened (see §4.1).

### 3.3 `availability_status` misattributes the dependency declaration
"*PYPI: … dependencies just Flask + Biopython*". PyPI's metadata declares **no dependencies at
all**: `info.requires_dist` is `null`, and the wheel's `METADATA` has no `Requires-Dist` lines
(`setup.py` has no `install_requires`). Flask + Biopython come only from `requirements.txt`,
which is not packaged. Further, **`codon_genie/species.table` is absent from both the 1.7 wheel
and the 1.7 sdist** (verified by unzipping/untarring both), so a pip-installed CodonGenie cannot
run `get_codon_usage_organisms()` at all — `main.py` would fail at import time, and the organism
list is unavailable. `CodonSelector.optimise_codons()` still works from a pip install *if* the
user separately installs Biopython. Practical effect: **the working reproduction path is
`git clone` + Docker, or the live REST API — not `pip install CodonGenie`.** This qualifies
`interface` bullet (d) but does not change the `yes`.

### 3.4 "213 is a typesetting corruption of 216"
The published PDF genuinely prints "reduced from 4,096 (16³) to **213** (6³)". I re-extracted the
surrounding text; the superscripts survive as trailing digits on both numbers ("163", "63"), so
"213" is the literal published value, not a glyph-loss artifact. It is an **arithmetic/typo error
in the paper**, not a PDF-extraction corruption. The practical advice (use **216**, never quote
213) is right; the *characterisation* is a guess and should not be asserted to a referee who
co-authored the paper. Say instead: "the paper prints 213; 6³ = 216."

### 3.5 memo §4 item 5 miscounts the NSS analysis
The memo says NSS returns "16 codons and **9** encoded amino acids". Live check: 16 codons, **8**
amino acids (A, C, G, P, R, S, T, W). The structured extraction's own `documented_examples` entry
lists the 16 codons correctly and does not repeat the "9", so only the memo needs fixing.

### 3.6 "the repository's only functional test" — **false**
There are **two** test modules. `codon_genie/test/test_client.py` contains five live-service
integration tests that the extractor never read (see §4.4).

### 3.7 The FILMV verification used a taxid that does not reproduce the paper's numbers
The extraction verifies against `organism=83333` and reports score **0.8152**, while quoting the
paper's **0.88** two lines earlier, without reconciling them. I resolved this:

| taxid | organism | top result |
|---|---|---|
| **37762** | *E. coli* (the id used in the repo's own unit test) | **DTK 0.88**, DTS 0.68 — **exactly the paper** |
| 83333 | *E. coli* K-12 | DTK 0.82, DTS 0.68 — close but ≠ paper |
| **1902** | *S. coelicolor* | **DTS 0.78, DTK 0.29** — matches the paper's 0.79/0.29 to rounding |
| 562 | *E. coli* (species-level) | **HTTP 500**, `{"message": "KeyError: 'ATC'"}` |

**Recommendation:** the reproduction recipe in the extraction should specify **37762** and
**1902**, not 83333, so that the numbers quoted in any rebuttal match the paper exactly.

---

## 4. Capabilities / facts the extractor MISSED

### 4.1 `ncbi_tax_utils.TaxonomyFactory` — taxonomic expansion and synonym resolution
An entire unread module. `TaxonomyFactory` downloads
`ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`, parses `nodes.dmp`/`names.dmp`, and exposes
`get_child_ids(parent_id)` (recursive descent of the taxonomy tree) and `get_names(tax_id)` (all
synonyms). `get_codon_usage_organisms(expand=True)` uses it to expand the Kazusa species table
with every child taxon and every name synonym. This makes the "organism resolution" additional
capability materially richer than the extraction describes — it is not just a flat name→taxid
lookup, it can walk NCBI Taxonomy. Not exposed via the REST API (`expand` defaults to `False`),
but it ships in the pip wheel. **Add to `additional_capabilities`.**

### 4.2 Per-amino-acid redundancy is quantified and surfaced
Every output amino-acid record carries a *list* of the codons in the expansion that encode it, and
the UI renders the **count** next to the letter (`{{aa.amino_acid}} {{aa.codons.length}}` in
`result.html`). Combined with the per-codon `probability` (organism codon-usage frequency), this
is an explicit statement of **expected variant-composition skew within the designed pool** — e.g.
DBK encodes V twice (GTG, GTT), so V is over-represented 2:1 relative to singly-encoded residues.
The extraction mentions "valine is encoded twice" only as a Fig. 1 detail and never registers it
as a capability. This is directly comparable to whatever PoolParty reports about library
composition, and a referee would notice its omission. **Add to `additional_capabilities`.**

### 4.3 A third-party fork added exactly the capability the extraction scores as absent
`RidgeBio/CodonGenie` (fork, pushed **2025-05-17**) adds
`CodonSelector.optimise_codons_seq(aa_seq, edits, organism_id)`: it parses an edit string
(`"A12V,A12L,..."` via a `_parse_edit` helper), validates each edit against the wild-type protein
sequence, groups multiple edits at the same position into one amino-acid set, and emits a
best-scoring codon for **every position of the full sequence** — plus matching UI changes.

Two consequences:

1. **It corroborates `library_as_object = no` and `mixed_mutagenesis_one_pool = no` for upstream
   CodonGenie**, and does so in the most persuasive possible way: a third party had to fork the
   tool to get multi-position, WT-anchored, edit-list design. That is a *good* fact for the
   rebuttal — but it must be framed as "a third-party fork", never as "CodonGenie can't do X".
2. **It is a live referee risk if unmentioned.** Upstream is unaware (no PR, no issue), it is not
   deployed, not on PyPI, and not the published tool — so it correctly does not change any cell.
   But the extraction should record its existence so the authors aren't blindsided.

### 4.4 `test/test_client.py` — five ready-made, paper-aligned reproduction fixtures
Unread. It contains assertions that are better reproduction targets than the one test the
extraction cites, because one of them is **the paper's own documented API example**:

- `client.get_codons('DE', 4932)` → 4 results, best = **`GAW`** ← this is literally
  `?aminoAcids=DE&organism=4932` from PDF p. 3
- `client.analyse('NTT', '4932')` → 4 amino acids
- `search_organisms('escherich')` and `search_organisms('escherichia co')` (space handling)
- `get_organisms()` non-empty

**Add to `documented_examples`.** Also correct "the repository's only functional test".

### 4.5 Robustness limitation: only Kazusa-listed taxids work, and failures are raw 500s
`?organism=562` (the canonical *E. coli* species taxid, and the first id most users would try)
returns **HTTP 500** with `{"message": "KeyError: 'ATC'"}` — because Kazusa has no usage table for
that id, `__get_codon_usage` silently builds an empty table, and the lookup explodes downstream.
There is no validation, no useful error, and no fallback. Combined with the plain-HTTP live
screen-scrape of `kazusa.or.jp` at every request (which the extraction *did* catch), this is a
substantive robustness limitation. **Add to `stated_limitations` as an "unstated but material"
item alongside the Kazusa dependency.**

---

## 5. Cells I tried hard to break and could not

- **`library_as_object` = no.** Attacked three ways: (i) the `/codons` response *is* a
  structured object — but it is a ranked list of candidate codons for **one** position, not a set
  of designed sequences; (ii) `align.py` builds a sequence — but the loop is 100% user code, the
  client just calls `get_codons` once per position and the caller does `dna_seq += ...`;
  (iii) maybe the deployed service is ahead of master — it is not, `optimise_codons` takes a single
  `amino_acids` string and does `set(amino_acids.upper())`, so `"FILMV,DE"` is one set, not two
  positions. Verdict stands.
- **`vcf_vep_output` = no.** Checked for a hidden export path: `result.html` and
  `result.directive.js` have no download/CSV/copy control, and every Flask route returns
  `application/json`. Nothing else emits any format. Stands.
- **`design_visualization` = partial.** Verified the colour-coded table and both tooltips exist in
  source and render live. A strict reading ("a table is not a visualisation") gives `no`; `partial`
  is the generous, referee-safe call. Stands.
- **`codon_optimization` = yes.** Doubly grounded and unassailable: organism-specific usage is the
  headline ranking axis (reproduced live: DTK/DTS reversal between taxids 37762 and 1902), *and*
  `CodonOptimiser.get_codon_optim_seq` is a real full-sequence optimiser that ships in the wheel
  (verified present in `CodonGenie-1.7-py3-none-any.whl`). Stands.
- **All Block C cells = no.** The entire input vocabulary is
  `{20 AA letters} ∪ {one 3-char IUPAC codon} ∪ {NCBI taxonomy id}`; the UI regex is
  `[acgtmrwsykvhdbnACGTMRWSYKVHDBN]{3}`; `requirements.txt` is two lines. Nothing coordinate-,
  annotation-, splice-, or HGVS-shaped exists. Stands, high confidence.

---

## 6. Convention flags for the orchestrator (must be applied uniformly across all tools)

1. **`maintained` threshold.** Define it explicitly (last commit? last release? archived flag?) and
   apply it identically. CodonGenie's `yes` will not survive contact with a matrix in which other
   2022-dormant tools are scored `no`.
2. **"Ships in the package" vs "user-facing".** The extraction already flags this. My verification
   confirms the tension is real and *asymmetric*: `codon_optimization` = yes survives on
   user-facing grounds alone, but `synthesis_constraints` = partial rests **entirely** on
   `seq_utils.find_invalid`, which no code path reachable from the web app or REST API ever calls.
   Under a strict user-facing convention it is `no`; under a package convention it is arguably
   `yes`. `partial` is the correct hedge only if the convention is explicitly "package-level
   capability, discounted for non-exposure" — write that convention down.
3. **Which URL to cite.** If the manuscript or rebuttal names a CodonGenie URL, it must **not** be
   the paper's `http://codon.synbiochem.co.uk` (404, verified twice). Use
   `https://codongenie.appspot.com` (live, and hard-coded as the default in `client.py`) or the
   GitHub repo. Note the repo URL in the paper (`synbiochem/CodonGenie`) redirects to
   `neilswainston/CodonGenie`.
4. **Tone.** The `ROW-LIST SUGGESTION` in `additional_capabilities` (add a
   "degenerate/ambiguous codon design" row where CodonGenie is the reference implementation and
   PoolParty is `no`) is the single best piece of judgement in the extraction. Keep it. It converts
   the most likely author-referee objection into a compliment.

---

## 7. Corrections to make in `extractions/codongenie.md`

| # | Location | Change |
|---|---|---|
| 1 | `maintained` cell | `yes` → `partial` (or `no`), per the cross-tool convention; delete the internal contradiction with "effectively dormant" |
| 2 | `primer_design` evidence | Drop the implication that a Tm routine exists; state that none exists anywhere |
| 3 | `barcode_generation` evidence | "4 Python modules" → "5"; name `ncbi_tax_utils.py` |
| 4 | `availability_status` | PyPI declares **no** dependencies (`requires_dist: null`); `species.table` is missing from wheel and sdist; the working install path is clone+Docker, not pip |
| 5 | `documented_examples` #2 | "typesetting corruption of 216" → "the paper prints 213; 6³ = 216" |
| 6 | `documented_examples` #6 | Remove "the repository's only functional test"; add `test_client.py` and its `get_codons('DE', 4932) → GAW` fixture |
| 7 | `documented_examples` #2/#3 | Use taxid **37762** (and **1902**) for reproduction — 83333 gives 0.82, not the paper's 0.88 |
| 8 | memo §4 item 5 | NSS → 8 encoded amino acids, not 9 |
| 9 | `additional_capabilities` | Add: NCBI taxonomy expansion/synonyms (`ncbi_tax_utils`); per-amino-acid redundancy counts + per-codon usage probability as expected pool-composition skew |
| 10 | `stated_limitations` | Add: unsupported taxids (e.g. 562) return HTTP 500 with a raw `KeyError`; no input validation |
| 11 | new note | Record the `RidgeBio/CodonGenie` fork (2025-05-17) adding `optimise_codons_seq` — as corroboration of `library_as_object = no`, and as pre-emption of a referee pointing at it |
