# Controlled vocabulary for the comparison table

Every term below is used with exactly one meaning throughout Table 1, its caption,
the supplementary matrix, and the Results prose that cites them. Chosen by
**measured usage**, not preference: term frequency counted across the 11 surveyed
tool papers and across `26.05.18_bmc_submission/latex/main.tex`.

Counts are `total occurrences (papers using it / 11)`, then our manuscript.

## Artifact nouns

| Use this | Not this | Evidence |
|---|---|---|
| **library** (or *DNA sequence library*) | *oligo pool*, *sequence pool*, *pool* | library 202 (9/11), ours 66. **"oligo pool" is 1 (1/11)** — it is our house term (4-6 uses in our ms), not field convention. Keep it in prose if you like; the table says *library*. |
| **variant** — a member of the library | *mutant* | variant 217 (9/11), ours 23. mutant only 6 (3/11). |
| **mutation** — the event (substitution / insertion / deletion) | — | mutation 183 (7/11). Use *substitution*, *deletion*, *insertion* for specific kinds. |
| **construct** — the full synthesised oligo including adapters and barcodes | — | 36 (8/11), ours 7. Use ONLY where the assembled-oligo distinction matters; otherwise say *variant* or *sequence*. |

## Capability terms — one term per concept, identical in every row

| Term | Means exactly | Evidence / note |
|---|---|---|
| **saturation mutagenesis** | Every single substitution at every position of a region, enumerated | 20 (5/11); **our ms uses it 0 times — adopt it.** Replaces the three phrasings in the first draft: *exhaustive scans*, *exhaustive single-event scans*, *exhaustive in-silico saturation mutagenesis*. |
| **pairwise and higher-order variants** | Two or more mutations co-occurring in one sequence | pairwise 22 (4/11), ours 4 — well attested. *higher-order* alone is weak (2, 1/11) but is our ms's term (4); anchoring it to *pairwise* makes it legible. |
| **randomly sampled variants** | Variants drawn by a rate, count or RNG parameter rather than enumerated | No single convention exists in the corpus, so define it in the caption. Do not conflate with *combinatorial*. |
| **combinatorial placement** | Elements (motifs, sites) placed at multiple positions and orientations | combinatorial 28 (4/11). |
| **mixed variant types in one library** | Structurally different component types in a single specification, pooled | Plain English; no field term exists. |
| **directed acyclic graph (DAG)** | PoolParty's architecture | ours: DAG 28, *directed acyclic graph* 4, **computational graph 0**. Spell out at first use, then DAG. Do **not** write *computational graph*. |

## Banned from the table

| Banned | Why | Say instead |
|---|---|---|
| **in silico** | Every tool in the table is computational, so it carries no information here. 14 uses (4/11) but always to contrast with wet-lab. | nothing — drop it |
| **exhaustive** | Weakly attested (8, 2/11) and vaguer than the alternative | *saturation mutagenesis*, or *enumerates all …* |
| **comprehensive** | 9 (6/11) but evaluative rather than descriptive | state what is enumerated |
| **single mutant** | **0 occurrences in all 11 papers**; 4 in our ms | *single substitution*, or *saturation mutagenesis* |
| **tiling / tile** | 0 in papers, 1 in ours | *scanning across positions* |
| **variant content** | my coinage, meaningless outside this project | *variants* |
| **design card** | 0 in papers, 18 in ours — our contribution, must be defined where introduced | not needed in this table |

## Vendor jargon: gloss or replace

| Their term | In the table write |
|---|---|
| *targeton* (VaLiAnT) | **target region** |
| *parametric deletion* (VaLiAnT) | **deletions of user-specified length** |
| *mutator*, *mutator type* (VaLiAnT) | **mutation type** |
| *Pool*, *Operation* (PoolParty, capitalised API nouns) | lowercase *library* and *operation* in the table; reserve capitals for the API description in the text |

## Consistency rule

If two rows describe the same capability, they must use the identical phrase. A
reader comparing two rows should never have to decide whether two different
wordings mean the same thing. The first draft failed this in four places; that
failure is what this file exists to prevent.
