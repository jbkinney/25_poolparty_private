# Table 1 draft v3 — tool overview (main text)

**Shape follows surveyed convention.** Four real comparison tables were examined
(see `CONVENTIONS.md`): tools are rows in 4/4; column headers are full words
with no invented acronyms; the dominant descriptive pattern is
*Tool · Purpose · Key features · Availability*. Two columns are added beyond that
pattern — *Output*, because the editor explicitly asked for "similarities and
differences in program outputs", and *Reference*, which is near-universal.

Terminology follows `TERMINOLOGY.md`.

---

## The table

| Tool | Purpose | Key features | Output | Availability | Reference |
|---|---|---|---|---|---|
| MPRA Design Tools | Design MPRA libraries and estimate statistical power | Barcode assignment with error-correcting sets; power analysis over effect size and sequencing depth | Barcoded constructs and statistical power estimates | R package (GitHub); Shiny app | [Ghazi2018aa] |
| MPRAnator | Design MPRA libraries of motif arrangements and sequence variants | Combinatorial placement of motifs; randomly sampled and combinatorial substitution sets; scrambled and reverse-complement controls | Oligo library (FASTA) | Web service | [GeorgakopoulosSoares2017gb] |
| **PoolParty** | Design DNA sequence libraries of any type | Libraries specified as a directed acyclic graph of operations; accepts sequences, FASTA, degenerate templates and motifs, or generates k-mers with no template; saturation mutagenesis, randomly sampled variants, and pairwise and higher-order variants can be mixed in one library; per-variant design cards | Library of sequences, each with a design card (CSV, TSV, FASTA, JSONL) | Python package (PyPI) | This work |
| VaLiAnT | Design and annotate libraries for saturation genome editing and cDNA deep mutational scanning | Saturation mutagenesis from genomic coordinates; several mutation types per target region; transcript-aware, including codons split across exons | Per-region oligo libraries with metadata tables and VCF carrying variant consequence annotation | Command-line tool (source, Docker) | [Barbon2022th] |
| Oligopool Calculator | Design pool infrastructure and analyse sequencing readout | Constraint-aware barcodes, primers and spacers; degenerate compression of a supplied variant set; off-target screening against a host genome | Synthesis-ready pool; variant counts from sequencing reads | Python package (PyPI); command line | [Hossain2024oc] |
| tangermeme | Probe cis-regulatory logic with deep learning models | Saturation mutagenesis of input sequences; insertion and removal of motifs across sets of sequences | Perturbed sequences and model predictions | Python package (PyPI) | [Schreiber2025nd] |
| General-purpose toolkits<br>*(Biopython, pydna, SeqPro)* | Manipulate sequences, simulate cloning, prepare arrays for modelling | Parsing, transformation and file-format support; no library abstraction | Transformed sequences | Python packages (PyPI) | [Cock2009df, Pereira2015wj, Klie2023kg] |
| Adjacent design tools<br>*(CodonGenie, DNA Chisel, ledidi, Mutation Maker)* | Optimise or edit individual sequences | Codon choice, synthesis-constraint satisfaction, model-guided editing, mutagenic primer design | An optimised sequence or a primer set | Web services and Python packages | [Swainston2017rb, Zulkower2020jk, Schreiber2025ledidi, Hiraga2021yg] |

**Footnote:** *Tools examined August 2026. Two rows group tools that address
problems adjacent to library design; members are named and cited individually.
Em dash (—) denotes a capability the tool does not provide.*

---

## Design decisions

**Row order is alphabetical within the comparable set, with PoolParty in place
rather than first.** The DMS review orders chronologically so its own subject is
not privileged. Bolding identifies PoolParty without promoting it. Trivially
changed if you prefer it first.

**Output earns its column.** It is the most direct response to the editor's ask
about program outputs, and reading down it makes the complementarity argument
structurally: VaLiAnT emits VCF with consequence annotation, Oligopool Calculator
emits a synthesis-ready pool from variants a user already has, adjacent tools emit
one sequence. Nobody has to be told these are different problems.

**No Pros/Cons column**, though two of the seven CRISPR-review tables use one.
Writing "Cons" about tools whose authors may referee this paper is precisely the
overstatement risk the editor flagged. Purpose and Output let the reader conclude.

**No genomic-context column.** It dissolves into VaLiAnT's Purpose
("from genomic coordinates", "transcript-aware") and Output ("VCF carrying variant
consequence annotation"), and into PoolParty's Key features, which names its input
types without claiming genomic integration. R3's concession lands without a
bespoke column a referee would notice we invented.

**No version, license or documentation-link columns.** Absent from all four
surveyed tables. Version is handled by the caption footnote; license belongs in
the supplement.

## On the two grouped rows

Convention lists every tool individually (4/4 surveyed). We deviate because this
table spans tool *categories*, which none of the surveyed tables does — the CRISPR
systematic review, which faces the same problem, answers it with seven separate
tables, which is worse here.

Three mitigations: members are named in the Tool cell, all members are cited in
Reference, and the footnote states the grouping rationale.

**Known instability:** tangermeme was moved out of the general-purpose group after
the audit found it satisfies saturation mutagenesis and composable design
operations. That grouping is a judgement call, not a fact, and the caption should
own it. It also means the Background sentence claiming general-purpose toolkits
"have no notion of a variant library as a structured object" needs narrowing.

## Open decision, unchanged

**tangermeme's row.** Given its own row here. Alternative: return it to the group
and narrow the Background sentence instead.
