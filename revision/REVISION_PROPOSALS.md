# Consolidated manuscript and response proposals

**Status:** Drafting underway; runtime and memory benchmarks drafted first.

**Writer:** Codex editorial orchestrator only.

**Baseline:** `main.tex` at commit `da98f615`; SHA-256
`cbdca45f1d67280e14c0256100d22f3a2a90cc531cfceb2bbdabafbe29dfc9e6`.

## Status and evidence labels

- Proposal: `Draft`, `Needs author decision`, `Approved`, `Applied`, `Verified`
- Evidence: `Merged`, `Branch/PR`, `Analysis`, `Incomplete`

## Coverage ledger

| Batch | Referee item | Stable IDs | Principal evidence | Status |
|---:|---|---|---|---|
| 1 | Editor: benchmarks, comparisons, and calibrated claims | E1 | `benchmarks/`, `comparison/` | Not drafted |
| 1 | R1.1 and R1.4: runtime and memory | LL1.1, LL4.1 | `benchmarks/` | Draft |
| 1 | R1.2 and R1.7: incremental value of design cards | LL2.1, LL11.1 | manuscript, package behavior | Not drafted |
| 1 | R1.3 and R1.5: systematic feature/code comparison | LL3.1, LL3.2, LL5.1 | `comparison/`, `code_comparison/` | Not drafted |
| 2 | R2.1a: library statistics | LL12.1, LL12.2 | `stats/` | Not drafted |
| 2 | R2.1b: LLM use and test ownership | LL13.1, LL13.2 | repository history and tests | Not drafted |
| 2 | R2.2a: backward/forward passes and supplemental figure | LL14.1 | `mechanism_figure/` | Not drafted |
| 2 | R2.2b: output equivalence with other tools | LL15.1 | `comparison/`, `code_comparison/` | Not drafted |
| 2 | R2.2c: filters and content constraints | LL16.1, LL16.2 | `constraint_filters/`, `stats/` | Not drafted |
| 2 | R2.2d: earlier PoolParty software name | LL17.1 | literature and manuscript | Not drafted |
| 3 | R3: scope and genome/transcript awareness | LL18.1 | `comparison/`, `vcf_input/` | Not drafted |
| 3 | R3: molecular consequence and VEP-compatible output | LL19.1 | `vcf_input/` | Not drafted |
| 3 | R3: LLM disclosure and journal policy | LL23.1, LL23.2 | manuscript and repository history | Not drafted |
| 4 | R1.6: suggested references | LL6.1–LL10.1 | cited papers and relevance assessment | Not drafted |
| 4 | R3: “lack” versus “scarcity” | LL20.1 | manuscript | Not drafted |
| 4 | R3: code blocks in Figures 2B and 3B | LL21.1 | manuscript and figures | Needs author decision |
| 4 | R3: definition of “most abundant codon” | LL22.1 | manuscript and implementation | Not drafted |
| 4 | R3: remaining line edits and positive comments | LL24.1, R3-M1–R3-M8 | manuscript and relevant evidence | Not drafted |

## Supplementary asset registry

Displayed numbers are provisional. Manuscript and response drafts should use
the stable semantic labels below; final Figure S/Table S numbers will be taken
from the compiled supplement after its contents and order settle.

| Stable label | Provisional number | Asset and source | Supports | Readiness |
|---|---|---|---|---|
| `fig:benchmark-scaling` | Figure S1 | `benchmarks/figS_scaling.pdf`; wrapper `comparison/figure_s1.tex`; caption `benchmarks/SUPPLEMENTARY_TEXT.md` | R1.1, R1.4 | Figure and caption available; final placement pending |
| `fig:sequence-generation` | Figure S2 | `mechanism_figure/fig_s2.pdf`; wrapper `comparison/figure_s2.tex`; caption `mechanism_figure/fig_s2_caption.md` | R2.2a | Content available; final content and print-legibility review pending |
| `tab:benchmark-summary` | Table S1 | benchmark results and caption in `benchmarks/SUPPLEMENTARY_TEXT.md`; current wrapper `comparison/tableS1.tex` | R1.1, R1.4 | Data available; replace the image wrapper with a native LaTeX table |
| `tab:library-statistics` | Table S2 | analysis in `stats/`; final table not yet assembled | R2.1a | Analysis available; table and caption pending |

`comparison/table1.tex` and `comparison/table2.tex` are currently treated as
main-text table proposals, not supplementary assets. Add or reorder registry
rows only when a concrete asset is proposed; changing displayed order must not
require changing semantic labels.

## Cross-comment author decisions

Decisions affecting multiple responses will be recorded here, including design-
card validation, figure code blocks, suggested references, scope language, and
LLM disclosure.

- **Approved workflow decision:** When multiple comments request the same
  substantive change, provide one complete primary response and a concise
  acknowledgment and cross-reference under each repeated comment.

## Drafted proposals

### R1.1 and R1.4 — Runtime and peak memory

**Decision:** Draft

**Primary comment (LL1.1):**

> In the GB1 example, the library exceeds 540,000 sequences, yet the paper
> provides no runtime or memory consumption benchmarks. Readers cannot assess
> the tool's practical usability across different library scales. The authors
> should supplement benchmarking data, including runtime and peak memory usage
> across libraries of increasing sizes (e.g., 1K, 10K, 100K, 500K sequences),
> and specify the test environment configuration.

**Repeated comment (LL4.1):**

> Provide benchmarking data on runtime and memory usage across libraries of
> varying sizes.

**Evidence — Analysis:**

- `benchmarks/paper_benchmark.json` contains 20 benchmark points, each based on
  five replicate runs. Requested sizes range from one sequence to the full
  example-specific generation size: 547,230 for GB1, 24,000 for MPRA, and
  200,000 for each of the two separately generated matched SpliceAI Pools.
- `benchmarks/figS_scaling.pdf` shows generation time, throughput, and peak
  resident memory across library size. Above the fixed-cost regime at small
  sizes, generation time and peak memory increase approximately linearly; the
  throughput differs among examples because their DAGs contain different
  Operations.
- At the reported generation sizes, mean generation time was 132.24 s for GB1,
  9.22 s for MPRA, and 11.93 s for one 200,000-sequence SpliceAI Pool. Peak
  memory was 487, 127, and 208 MB, respectively. The SpliceAI example generates
  two matched Pools separately, so the total generation time is approximately
  twice the per-Pool value.
- The hardware, software versions, timing method, peak-memory method, and
  replicate design are recorded in `benchmarks/SUPPLEMENTARY_TEXT.md` and the
  benchmark script.

**Manuscript edit 1 — add the benchmark setup to Methods:**

- Location: after baseline line 148, at the end of `Implementation and
  extensibility` and before `SpliceAI surrogate analysis`.
- Current text: the manuscript does not describe the benchmark setup.
- Proposed insertion:

```latex
We benchmarked generation of the three worked-example libraries at their
reported sizes and at several smaller sizes. Benchmarks were run on a single
workstation under Ubuntu 22.04.3 LTS in Windows Subsystem for Linux 2 (kernel
5.10.16.3-microsoft-standard-WSL2), with an Intel Core Ultra 9 185H processor
and 31~GB RAM allocated to WSL2. We used CPython 3.12.2, PoolParty 0.1.0,
StateTracker 0.1.0, NumPy 2.2.5, pandas 2.3.3, and beartype 0.22.9. Timings are
single-threaded wall-clock measurements of \texttt{generate\_library}, reported
as mean $\pm$ SD over five replicate runs. Peak memory is the maximum resident
set size reported by \texttt{getrusage}, measured in a fresh process for each
data point.
```

**Manuscript edit 2 — cite the supplementary benchmarks:**

- Location: append to the existing GB1 Results paragraph at baseline line 179.
- Current anchor: `Users can also inspect the DAG structure by calling the
  \texttt{print\_dag} method of the root Pool (Fig.~\ref{fig:figure2}D).`
- Proposed appended sentence:

```latex
Runtime and peak memory benchmarks for this and the other two examples are
provided in Supplementary Fig.~\ref{fig:benchmark-scaling} and Supplementary
Table~\ref{tab:benchmark-summary}.
```

**Supplementary additions:**

- Include `benchmarks/figS_scaling.pdf` with the caption in
  `benchmarks/SUPPLEMENTARY_TEXT.md` and label it
  `fig:benchmark-scaling` (provisionally Figure S1).
- Set the benchmark values in `benchmarks/SUPPLEMENTARY_TEXT.md` as a native
  LaTeX table with label `tab:benchmark-summary` (provisionally Table S1). The
  existing `benchmarks/tableS1.pdf` is a draft preview, not the final table.

**Draft response to LL1.1:**

We thank the reviewer for pointing out this omission. We have now benchmarked
runtime and peak memory for the three main-text examples across a range of
library sizes. For example, the full GB1 library (547,230 sequences) took
$132 \pm 6$~s to generate and used 487~MB of peak memory on our test machine.
Runtime and memory both increase approximately linearly with library size,
apart from some initial overhead, and the per-sequence cost depends on the
operations used. We have added these results to Supplementary
Figure~\ref{fig:benchmark-scaling} and Table~\ref{tab:benchmark-summary}, with
details of the benchmarking setup in the Methods.

**Draft response to LL4.1:**

This request is addressed together with Comment 1 above.
