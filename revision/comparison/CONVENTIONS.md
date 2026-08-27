# Convention survey: how bioinformatics tool-comparison tables are built

Evidence base for the design of Tables 1 and 2. Four real comparison tables
examined, plus the 11 tool papers tracked at `../../lit_review/analyzed/`.

## Sources

| Source | Scope | Shape |
|---|---|---|
| Variant scoring tools for deep mutational scanning (2025), PMC12494760 | 12 tools x 10 columns | tools = rows |
| CRISPR-Cas9 and Its Bioinformatics Tools: A Systematic Review, PMC12109910 | **7 separate tables**, ~9 tools each | tools = rows |
| Comparative Analysis of Codon Optimization Tools (2025), PMC12010093 | 10 tools x design criteria | tools = rows |
| The 11 tool papers we cite | — | **none contains a tool-vs-tool comparison table** |

## Column headers, verbatim

**DMS review:** Name · Variant populations *(Two-population | Time-series)* ·
Input sequencing data types *(Direct | Tiles | Barcodes | Count table)* ·
Scoring approach · Scoring method · Scoring statistic · Language · Interface ·
Documentation Link · Reference

**CRISPR review**, across its seven tables: Tool Name · Purpose · Key Features ·
Availability · Scope · Pros/Cons · Data Coverage

**Codon optimization:** design criteria as columns, e.g. CAI (Codon Adaptation
Index), Individual Codon Usage (ICU), Codon Context (CC)

## Patterns

1. **Tools are rows.** 4/4.
2. **Headers are full words.** The CRISPR review is explicit: *"Full words/phrases
   only — no abbreviations or acronyms in column titles."* Acronyms appear only
   where they are established domain terms, and are expanded in the header itself
   (CAI, ICU, CC). Invented shorthand does not occur.
3. **`Purpose` + `Key Features` is the dominant descriptive pair** — all seven
   CRISPR tables use it.
4. **Two distinct styles**, which do not mix comfortably:
   - *descriptive* — text cells of 15–80 words (CRISPR review)
   - *matrix* — tools x criteria, checkmarks (codon optimization)
   - the DMS review hybridises and is the least tidy of the three
5. **Grouped sub-columns** are how matrix tables stay readable, with full-word
   sub-headers.
6. **Absent from all three:** license, version examined, installation method.
   Documentation link appears once (DMS review); Reference appears once.
7. **Separate tables per category** is the conventional answer to a heterogeneous
   tool set — the CRISPR review uses seven.

## How our tables follow and deviate

| Choice | Convention | Ours |
|---|---|---|
| Tools as rows | 4/4 | follows |
| Full-word headers | 4/4 | follows |
| Purpose + Key features | 7/7 CRISPR tables | follows (Table 1) |
| Descriptive vs matrix | kept separate | follows — Table 1 descriptive, Table 2 matrix |
| Reference column | 1/3 | follows (adds it) |
| Output column | not standard, but DMS review has *Input data types* | **adds** — the editor explicitly asked about program outputs |
| Version / license / docs columns | absent 3/3 | follows (omits; version in caption) |
| Pros/Cons | 2/7 CRISPR tables | **declines** — editorialising about tools whose authors may referee |
| Grouped rows for adjacent categories | none; CRISPR uses 7 tables instead | **deviates** — two grouped rows, members named and cited |
| Two symbol levels | 2/2 matrix tables | **deviates** — three levels (✓ / ~ / —), because the audits produced genuine partials |
| Tools as rows in the feature matrix | 4/4, incl. the codon-optimization matrix | **deviates** — Table 2 is pivoted, capabilities as rows. The codon matrix affords criteria-as-columns only because its criteria are short acronyms (CAI, ICU, CC); our full-word capability names would need three-line stacked headers. Pivoting is what makes the full-word decision affordable, and it creates room for a per-row gloss carrying the operational definition. |

Both deviations are recorded here so they can be defended in the response letter
rather than discovered by a referee.
