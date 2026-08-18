"""Reproduce MPRAnator's "MPRA Motif Use Case" library in PoolParty and diff it.

Two independent implementations are compared:

  reference()  -- MPRAnator's documented algorithm, re-expressed here from its
                  published description and source reading. Their code is
                  all-rights-reserved (no LICENCE file) and is NOT vendored;
                  running their own modules reproduces these same 976 arrangements
                  and the documented 5,856 total. See mpranator_motif_usecase.md.
  build()      -- the PoolParty expression of the same design.

    python3 repro_mpranator_motifs.py
"""

import itertools
import poolparty as pp

BG = {
 "bg1_chr6:77195320":"tgtgtcttaaaaaaacaaacaaacaaacaaaatcccgaaataaaacacaacaaaaaaaaccccaccccataatcttcaggacagtctgtc",
 "bg2_chr9:37271330":"gtatctactctctgcccttacaacctcctcccagaaagaataaaatgtttctcatcctggaagctacagtgtgtcacacagtatactctt",
}
MOT   = {"AP1":"TGACTCA", "ELK1":"ACCGGAAGT", "RFX":"CGTTGCTAGGCAACG"}
MIN, MAX, EDGE, INTERVAL, L = 6, 24, 15, 12, 90
RIGHT = L - EDGE                      # motif end <= 74  =>  start <= 75 - n

def build(bg_seq):
    """One Pool per (motif subset, left-to-right ordering).

    MPRAnator constrains only the LEFTMOST motif to a multiple of the interval of
    substitution. In PoolParty that is insertion_mode='ordered' (positions strictly
    increasing) with the first insert's candidate positions restricted; the motif
    permutations MPRAnator gets from itertools.product over ordered position tuples
    become the k! separate orderings stacked here.
    """
    pools = []
    bg = pp.from_seq(bg_seq.upper())
    for k in (1, 2, 3):
        for subset in itertools.combinations(MOT, k):
            for order in itertools.permutations(subset):
                seqs = [MOT[n] for n in order]
                positions = []
                for i, s in enumerate(seqs):
                    hi = RIGHT + 1 - len(s)                    # inclusive upper bound
                    cand = range(EDGE, hi)
                    positions.append([p for p in cand if p % INTERVAL == 0] if i == 0
                                     else list(cand))
                if any(len(p) == 0 for p in positions):
                    continue
                tag = "_".join(order)
                kw = dict(positions=positions, replace=True,
                          insertion_mode="ordered", mode="sequential",
                          prefix="+".join(order),
                          # unique region names: motif lengths differ between calls
                          # and region lengths are Party-global (cf. FINDINGS B3)
                          names=[f"{tag}_{i}" for i in range(k)])
                if k > 1:
                    kw.update(min_spacing=MIN, max_spacing=MAX)
                pools.append(pp.insertion_multiscan(bg, k, [pp.from_seq(s) for s in seqs], **kw))
    return pools


def reference():
    """MPRAnator's four filters, applied to ordered position tuples.

    1  placed motif must not run past the end
    2  gap between consecutive motifs in [MIN, MAX]
    3  rightmost motif end <= L-1-EDGE
    4  the LEFTMOST motif start is a multiple of INTERVAL
    """
    out = set()
    for bg_seq in BG.values():
        bg = bg_seq.upper()
        for k in (1, 2, 3):
            for subset in itertools.combinations(MOT, k):
                motifs = [MOT[n] for n in subset]
                for pos in itertools.product(range(EDGE, L), repeat=k):
                    if any(p + len(m) > L for p, m in zip(pos, motifs)):
                        continue                                          # 1
                    spans = sorted((p, p + len(m) - 1) for p, m in zip(pos, motifs))
                    gaps = [spans[i + 1][0] - spans[i][1] - 1 for i in range(k - 1)]
                    if any(g < MIN or g > MAX for g in gaps):
                        continue                                          # 2
                    if spans[-1][1] > L - 1 - EDGE:
                        continue                                          # 3
                    if min(pos) % INTERVAL:
                        continue                                          # 4
                    s = list(bg)
                    for p, m in zip(pos, motifs):
                        s[p:p + len(m)] = m
                    out.add("".join(s))
    return out


def main():
    total = set()
    per_bg = {}
    for name, seq in BG.items():
        pp.init()
        pools = build(seq)
        lib = pp.stack(pools)
        s = set(lib.to_df(num_cycles=1, show_progress=False)["seq"])
        per_bg[name] = (lib.num_states, len(s))
        total |= s
        print(f"  {name:22s} pools={len(pools):2d} states={lib.num_states:4d} distinct={len(s):4d}")
    print(f"\n  PoolParty combined distinct: {len(total)}")
    ref = reference()
    print(f"  MPRAnator reference        : {len(ref)}")
    print(f"  EXACT overlap              : {len(total & ref)} "
          f"({100 * len(total & ref) / len(ref):.1f}%)")
    print(f"  only PoolParty / only ref  : {len(total - ref)} / {len(ref - total)}")
    return total

if __name__ == "__main__":
    main()
