# Supplementary performance text

Paste-ready captions and Methods text. Numbers come from
`paper_benchmark.json` (20 points, 5 replicate runs each).

Regenerate the data with `python paper_benchmark.py`; redraw the figure and the
table image with `python make_figure.py`. The table image is for circulating a
draft — set the final version in LaTeX.

---

## Table S1 — caption

> **Runtime and memory for the three worked examples.** Each library was
> generated at its published size on a single core. Times are means ± SD over 5
> replicate runs; peak memory is maximum resident set size, including a ~106 MB
> baseline for the Python interpreter, NumPy, and pandas. DAG construction does
> not depend on library size and is reported separately. The SpliceAI
> library comprises two matched 200,000-sequence pools (cryptic GT and disrupted
> GA) generated separately; values are for one pool, and the full library takes
> twice as long.

```latex
\caption{\textbf{Runtime and memory for the three worked examples.} Each library
was generated at its published size on a single core. Times are means $\pm$ SD
over 5 replicate runs; peak memory is maximum resident set size, including a
$\sim$106~MB baseline for the Python interpreter, NumPy, and pandas. DAG
construction does not depend on library size and is reported separately.
The SpliceAI library comprises two matched 200,000-sequence pools (cryptic GT
and disrupted GA) generated separately; values are for one pool, and the full
library takes twice as long.}
```

| Example | Library size | Time (s) | Throughput (seq s⁻¹) | Peak memory (MB) | DAG (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| GB1 deep mutational scan | 547,230 | 132.24 ± 5.76 | 4,138 | 487 | 0.22 |
| MPRA regulatory grammar | 24,000 | 9.22 ± 0.22 | 2,602 | 127 | 0.06 |
| SpliceAI surrogate | 200,000 | 11.93 ± 0.27 | 16,758 | 208 | 0.05 |

---

## Figure S1 — caption

> **Library generation scales linearly in time and memory.** **(A)** Wall-clock
> generation time, **(B)** throughput, and **(C)** peak resident memory against
> the number of sequences generated, for the three worked examples. Points are
> means of 5 replicate runs; error bars show ±1 SD and are smaller than the
> plotting symbols at most points. Below ~10³ sequences a fixed setup cost
> dominates, giving the plateau in **(A)** and the rise in **(B)**; above 10⁴
> sequences throughput is constant, so generation time is linear in library
> size. Peak memory increases linearly, by 0.52–0.66 kB per sequence above a
> ~106 MB baseline.

```latex
\caption{\textbf{Library generation scales linearly in time and memory.}
\textbf{(A)} Wall-clock generation time, \textbf{(B)} throughput, and
\textbf{(C)} peak resident memory against the number of sequences generated, for
the three worked examples. Points are means of 5 replicate runs; error bars show
$\pm$1~SD and are smaller than the plotting symbols at most points. Below
$\sim$10$^3$ sequences a fixed setup cost dominates, giving the plateau in
\textbf{(A)} and the rise in \textbf{(B)}; above 10$^4$ sequences throughput is
constant, so generation time is linear in library size. Peak memory increases
linearly, by 0.52--0.66~kB per sequence above a $\sim$106~MB baseline.}
```

---

## Methods — benchmarking environment

> Benchmarks were run on a single workstation: Ubuntu 22.04.3 LTS under Windows
> Subsystem for Linux 2 (WSL2; kernel 5.10.16.3-microsoft-standard-WSL2), Intel
> Core Ultra 9 185H (22 logical processors), 31 GB RAM allocated to the WSL2
> virtual machine, CPython 3.12.2, PoolParty 0.1.0, StateTracker 0.1.0, NumPy
> 2.2.5, pandas 2.3.3, beartype 0.22.9. All timings are single-threaded,
> in-memory wall-clock measurements of `generate_library`, reported as mean ± SD
> over 5 replicate runs. Peak memory is the maximum resident set size reported
> by `getrusage`, measured in a fresh process for each data point.

```latex
Benchmarks were run on a single workstation: Ubuntu 22.04.3 LTS under Windows
Subsystem for Linux 2 (WSL2; kernel 5.10.16.3-microsoft-standard-WSL2), Intel
Core Ultra 9 185H (22 logical processors), 31~GB RAM allocated to the WSL2
virtual machine, CPython 3.12.2, PoolParty 0.1.0, StateTracker 0.1.0, NumPy
2.2.5, pandas 2.3.3, beartype 0.22.9. All timings are single-threaded, in-memory
wall-clock measurements of \texttt{generate\_library}, reported as mean $\pm$ SD
over 5 replicate runs. Peak memory is the maximum resident set size reported by
\texttt{getrusage}, measured in a fresh process for each data point.
```

**Note, not manuscript text.** "In-memory" is the one claim here a reviewer could
challenge, given the benchmark ran under WSL2 with the package on the `/mnt/c`
DrvFs mount. It was checked: the same workload against a copy of the package on
the native ext4 filesystem gave 10.24 ± 0.05 s, against 10.13 ± 0.09 s on
`/mnt/c`. Only module import is filesystem-sensitive (0.56 s vs 2.12 s), and
import is not part of any reported measurement.

---

## Availability and requirements — suggested edit

The current block gives no memory guidance. Suggested addition in bold:

> **Other requirements:** NumPy (≥ 1.20), pandas (≥ 1.3), beartype (≥ 0.22.9),
> pyfaidx (≥ 0.8.1), typing_extensions (≥ 4.0). The companion `statetracker`
> package is installed automatically alongside PoolParty. **Generating a library
> requires approximately 0.6 kB of memory per sequence, i.e. under 1 GB for
> libraries of up to 10⁶ sequences.**

---

## Not reported, and why

- **No comparison against other tools.** None of the seven design tools cited in
  the manuscript report runtime benchmarks to compare against. Single-tool
  scaling curves are normal in this literature.
- **No fitted scaling exponent.** The throughput plateau in Fig. S1B shows
  linearity directly and is easier to read than a power-law exponent estimated
  from a handful of points.
- **No toggle micro-optimisations.** `pp.toggle_styles(False)` and
  `pp.toggle_cards(False)` save a modest, workload-dependent amount of time.
  Implementation detail.
- **GB1 is reported as published**, including the cost of `pp.stack` evaluating
  every branch for each generated sequence. Its four components generated
  separately sum to about a third of the stacked time. See
  `dev/pre_publish_audit.md` Phase 8 (H3) for the known optimisation that has
  not been applied.
