#!/usr/bin/env python
"""Publication benchmark for the PoolParty manuscript (supplementary table + figure).

Measures wall-clock time and peak resident memory for the three worked examples
across a sweep of library sizes.

    python paper_benchmark.py            # run the full sweep -> paper_benchmark.json
    python paper_benchmark.py --quick    # small sizes only, for a smoke test
    python paper_benchmark.py --point mpra 10000    # internal: one measurement

Design notes (chosen to match reporting conventions in the field):

* **Peak memory is maxRSS**, read from `resource.getrusage(RUSAGE_CHILDREN)`.
  maxRSS is a high-water mark that never decreases within a process, so every
  measurement runs in its OWN subprocess -- otherwise the first large run would
  contaminate every later one. This is the same quantity `/usr/bin/time -v`
  reports. tracemalloc is deliberately NOT used: it tracks only Python-level
  allocations and misses the NumPy/pandas buffers that dominate here.
* **5 repetitions, reported as mean +- SD.** Not best-of-N, which is a
  software-engineering idiom rather than a publication one.
* Timing brackets `generate_library()` only -- DAG construction and imports are
  measured separately, since they are size-independent fixed costs.
"""
from __future__ import annotations

import argparse
import json
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "paper_benchmark.json"

# --------------------------------------------------------------------------
# The three worked examples, exactly as published.
# --------------------------------------------------------------------------
ORF = ("ATGCAGTACAAGCTGATCCTGAACGGTAAGACGCTGAAAGGTGAGACGACCACCGAAGCTGTAG"
       "ACGCTGCTACTGCAGAGAAGGTGTTCAAGCAGTACGCTAACGACAACGGCGTCGACGGTGAATG"
       "GACCTACGACGACGCTACCAAAACCTTCACGGTTACCGAA")
CARDS = {"codon_positions": "position", "wt_aas": "wt_aa", "mut_aas": "mut_aa"}
BG1_100 = ("GCAAGTCTGCCATCGTGTTCAGAAGGGCCAGAAATGCCAAGGACTCAGGGGAGG"
           "AGAATTAAGTCAGAGAGTTTCATTACTGAGTGTTGTTTGACTTTGT")
MPRA_TEMPLATE = ("ACTGGCCGCTTCACTG<cre>" + BG1_100 + "</cre>GGTACCTCTAGA<bc>"
                 + "N" * 8 + "</bc>AGATCGGAAGAGCGTCG")
SPLICEAI_POSITIONS = list(range(51, 90)) + list(range(107, 168))


def build_dms(pp):
    """Example 1: GB1 deep mutational scan. 547,230 sequences."""
    orf = pp.from_seq(ORF).named("orf_pool")
    pos = slice(1, 56)
    single = orf.mutagenize_orf(num_mutations=1, mutation_type="missense_only_first",
        codon_positions=pos, prefix="single", style="red", mode="sequential",
        cards=CARDS).named("single_pool")
    double = orf.mutagenize_orf(num_mutations=2, mutation_type="missense_only_first",
        codon_positions=pos, prefix="double", style="red", mode="sequential").named("double_pool")
    rand = orf.mutagenize_orf(mutation_rate=0.1, mutation_type="missense_only_first",
        codon_positions=pos, prefix="random", style="red", mode="random",
        num_states=10000).named("random_pool")
    wt = orf.repeat(times=100, prefix="wt").named("wt_pool")
    return pp.stack([single, double, rand, wt])


def build_mpra(pp):
    """Example 2: MPRA regulatory grammar. 24,000 sequences."""
    t = pp.from_seq(MPRA_TEMPLATE)
    h = pp.from_seq("GGGGCAAAGGTCA", style="blue").flip(mode="sequential", cards={"flip": "hnf4a_strand"})
    p = pp.from_seq("CCGGGTCATTGGGGTCAGG", style="purple").flip(mode="sequential", cards={"flip": "ppara_strand"})
    x = pp.from_seq("GTGATGACGTGTCCCAT", style="orange").flip(mode="sequential", cards={"flip": "xbp1_strand"})
    cre = t.insertion_multiscan(region="cre", insertion_pools=[h, p, x],
        insertion_mode="unordered", replace=True, min_spacing=0, num_insertions=3,
        mode="random", num_states=1000, names=["hnf4a", "ppara", "xbp1"],
        cards={"starts": "positions", "names": "tfbs"}).repeat(times=3)
    bc = pp.get_barcodes(num_barcodes=cre.num_states, length=8, gc_range=(0.3, 0.6),
                         min_edit_distance=1, style="bold", seed=42)
    return cre.replace_region(region_name="bc", content_pool=bc)


def build_spliceai(pp):
    """Example 3: SpliceAI surrogate, one arm. 200,000 sequences."""
    import pandas as pd
    bg = (HERE / "spliceai_background_201bp.txt").read_text().strip()
    seqs = pd.read_csv(HERE / "spliceai_cryptic_seqs.csv").cryptic_sequence.tolist()
    cryptic = pp.from_seqs(seqs, mode="sequential",
                           cards={"seq": "cryptic_sequence"}).named("cryptic")
    return pp.from_seq(bg).replacement_scan(
        cryptic, positions=SPLICEAI_POSITIONS, mode="sequential",
        cards={"start": "cryptic_position"}).named("library")


EXAMPLES = {"mpra": build_mpra, "spliceai": build_spliceai, "dms": build_dms}
NATIVE_SIZE = {"mpra": 24_000, "spliceai": 200_000, "dms": 547_230}

# Sweep DOWN from each library's own size rather than past it. Two benefits:
#   * n=1 measures the fixed cost directly -- peak RSS there IS the baseline
#     (interpreter + NumPy + pandas + the built DAG), so the figure shows the
#     floor rather than hiding it in an intercept.
#   * every point is a real library the design can actually produce; nothing
#     re-enumerates states, so no point needs an asterisk.
# Integer decades only. An earlier 1-3-10 ladder was tried; the extra points
# changed nothing that mattered (they were originally added to stabilise a
# power-law fit that is no longer reported -- linearity is shown by the
# throughput plateau instead), so the denser ladder was dropped.
LADDER = [1, 10, 100, 1_000, 10_000, 100_000]


def sizes_for(example: str) -> list[int]:
    n = NATIVE_SIZE[example]
    return [s for s in LADDER if s < n] + [n]


# --------------------------------------------------------------------------
# child: one measurement, in its own process so maxRSS is uncontaminated
# --------------------------------------------------------------------------
def run_point(example: str, num_seqs: int) -> dict:
    import poolparty as pp
    pp.init()
    t0 = time.perf_counter()
    pool = EXAMPLES[example](pp)
    build_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    out = pool.generate_library(num_cycles=1, num_seqs=num_seqs, seed=42)
    gen_s = time.perf_counter() - t0
    n = len(out)
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {"example": example, "requested": num_seqs, "n": n,
            "build_s": build_s, "gen_s": gen_s, "peak_mb": peak_kb / 1024}


# --------------------------------------------------------------------------
# parent: sweep
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--point", nargs=2, metavar=("EXAMPLE", "N"))
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()

    if a.point:
        print(json.dumps(run_point(a.point[0], int(a.point[1]))))
        return

    records = []
    for example in ["mpra", "spliceai", "dms"]:
        sizes = sizes_for(example)
        if a.quick:
            sizes = [s for s in sizes if s <= 10_000]
        n_reps = a.reps
        for size in sizes:
            reps = []
            for _ in range(n_reps):
                r = subprocess.run(
                    [sys.executable, __file__, "--point", example, str(size)],
                    capture_output=True, text=True, cwd=str(HERE))
                if r.returncode != 0:
                    print(f"  !! {example} n={size} FAILED\n{r.stderr[-600:]}", flush=True)
                    break
                reps.append(json.loads(r.stdout.strip().split("\n")[-1]))
            if not reps:
                continue
            gen = [x["gen_s"] for x in reps]
            mem = [x["peak_mb"] for x in reps]
            bld = [x["build_s"] for x in reps]
            rec = {
                "example": example, "n": reps[0]["n"], "reps": len(reps),
                "gen_mean": statistics.mean(gen),
                "gen_sd": statistics.stdev(gen) if len(gen) > 1 else 0.0,
                "peak_mb_mean": statistics.mean(mem),
                "peak_mb_sd": statistics.stdev(mem) if len(mem) > 1 else 0.0,
                "build_mean": statistics.mean(bld),
                "native": reps[0]["n"] == NATIVE_SIZE[example],
                "repeats_states": reps[0]["n"] > NATIVE_SIZE[example],
            }
            records.append(rec)
            print(f"  {example:9s} n={rec['n']:>7,}  "
                  f"gen {rec['gen_mean']:8.2f} +- {rec['gen_sd']:5.2f} s  "
                  f"peak {rec['peak_mb_mean']:7.1f} +- {rec['peak_mb_sd']:5.1f} MB"
                  f"{'   <- published size' if rec['native'] else ''}", flush=True)

    OUT.write_text(json.dumps(records, indent=1))
    print(f"\nwrote {OUT}  ({len(records)} points)")
    print("BENCHDONE")


if __name__ == "__main__":
    main()
