"""Express the SpliceAI GT-arm library in PoolParty and in tangermeme, and diff.

The library: 2,000 cryptic 5' splice sites of graded MaxEntScan strength,
substituted at each of 100 positions across a 201 bp region carrying the
canonical 5' splice site of SMN2 exon 7. 2,000 x 100 = 200,000 sequences.

This is the GT arm of the manuscript's Example 3 -- not the whole example. The
full example also builds the 9-mer ladder from a PWM, slices the region and two
5 kb flanks from GRCh38, and adds a matched GA control arm in which the cryptic
GT is disrupted. Those stages are out of scope here.

Ground truth is examples/spliceai_design_cards_published.csv.gz in the poolparty
repository, filtered to condition == "GT".

The two tools need different environments, so run each side separately:

    conda activate poolparty_dev   # or any env with poolparty on main
    python repro.py poolparty

    conda activate tangermeme
    python repro.py tangermeme

Each side is compared independently against the published cards, so if both
report an identical set and order they are also identical to each other. No
intermediate file is written.
"""

import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", "..", ".."))

# Inputs, both committed in this repository.
TARGET_REGION_FILE = os.path.join(_REPO, "25_poolparty_private", "revision",
                                  "benchmarks", "spliceai_background_201bp.txt")
CRYPTIC_SITES_FILE = os.path.join(_REPO, "25_poolparty_private", "revision",
                                  "benchmarks", "spliceai_cryptic_seqs.csv")
# Ground truth, from the poolparty repository's examples directory.
PUBLISHED_CARDS = os.path.join(_REPO, "poolparty-statecounter", "poolparty",
                               "examples", "spliceai_design_cards_published.csv.gz")

# Positions 90-106 are omitted: the canonical 5' splice site sits at 100.
SCAN_POSITIONS = list(range(51, 90)) + list(range(107, 168))

KEY = ["cryptic_position", "cryptic_sequence", "full_sequence"]


def inputs():
    target_region = open(TARGET_REGION_FILE).read().strip()
    cryptic_sites = pd.read_csv(CRYPTIC_SITES_FILE).cryptic_sequence.tolist()
    return target_region, cryptic_sites


# --------------------------------------------------------------------------
# PoolParty. Three statements. The design card columns are named by the
# operations that produce them: the pool contributes the 9-mer, the scan
# contributes the position.
# --------------------------------------------------------------------------
def build_poolparty():
    import poolparty as pp

    target_region, cryptic_sites = inputs()

    pp.init()
    cryptic_pool = pp.from_seqs(cryptic_sites, mode="sequential",
                                cards={"seq": "cryptic_sequence"})
    library = pp.from_seq(target_region).replacement_scan(
        cryptic_pool, positions=SCAN_POSITIONS, mode="sequential",
        cards={"start": "cryptic_position"})

    print(f"  num_states = {library.num_states:,}  (nothing generated yet)")
    df = library.generate_library()
    return df.rename(columns={"seq": "full_sequence"})[KEY]


# --------------------------------------------------------------------------
# tangermeme. Three statements, written as a tangermeme user would: one
# batched substitution per position, sequences left as one-hot tensors in the
# shape a model wants. `expand` rather than `repeat` because `substitute`
# clones internally, so no copy is needed here.
#
# Everything after `library_ohe` exists only to make the diff possible, and is
# not how the library would be used.
# --------------------------------------------------------------------------
def build_tangermeme():
    import torch
    from tangermeme.ersatz import substitute
    from tangermeme.utils import characters, one_hot_encode

    target_region, cryptic_sites = inputs()

    target_ohe = one_hot_encode(target_region).unsqueeze(0).expand(
        len(cryptic_sites), -1, -1)
    cryptic_ohe = torch.stack([one_hot_encode(site) for site in cryptic_sites])
    library_ohe = torch.stack([substitute(target_ohe, cryptic_ohe, start=p)
                               for p in SCAN_POSITIONS])

    mb = library_ohe.element_size() * library_ohe.nelement() / 1e6
    print(f"  library_ohe {tuple(library_ohe.shape)}  {mb:.0f} MB, "
          f"materialised, position-major")

    # --- comparison scaffolding, not part of using tangermeme ---
    # Reorder position-major -> site-major to match PoolParty's enumeration,
    # then decode to strings. `characters` takes one sequence at a time.
    flat = library_ohe.transpose(0, 1).flatten(0, 1)
    return pd.DataFrame({
        "cryptic_position": SCAN_POSITIONS * len(cryptic_sites),
        "cryptic_sequence": np.repeat(cryptic_sites, len(SCAN_POSITIONS)),
        "full_sequence": [characters(seq) for seq in flat],
    })[KEY]


def published():
    df = pd.read_csv(PUBLISHED_CARDS, usecols=KEY + ["condition"])
    return df[df.condition == "GT"][KEY].reset_index(drop=True)


def compare(name, df):
    ref = published()
    df = df.reset_index(drop=True)
    as_set = set(map(tuple, ref.itertuples(index=False)))
    mine = set(map(tuple, df.itertuples(index=False)))
    print(f"  rows            : {len(df):,}  (published GT arm: {len(ref):,})")
    print(f"  identical set   : {as_set == mine}")
    print(f"  identical order : {ref.equals(df)}")
    print(f"  only in {name:<9}: {len(mine - as_set):,}")
    print(f"  only published  : {len(as_set - mine):,}")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else ""
    if which == "poolparty":
        compare("poolparty", build_poolparty())
    elif which == "tangermeme":
        compare("tangermeme", build_tangermeme())
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
