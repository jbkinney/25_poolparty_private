"""Reproduce VaLiAnT's brca1_pep library in PoolParty and diff it, exon by exon.

Ground truth is VaLiAnT's own committed expected output, which its repository
validates by md5sum -- VaLiAnT itself is never installed. The four _meta.csv files
are downloaded on first run.

See valiant_brca1_pep.md for the full write-up.

    python3 repro_valiant_brca1_pep.py
"""

import copy
import csv
import os
import urllib.request

import poolparty as pp
from poolparty.codon_table import STANDARD_GENETIC_CODE as G

RAW = ("https://raw.githubusercontent.com/cancerit/VaLiAnT/develop/"
       "examples/sge/brca1/brca1_pep_output_exp")
META = [
    "chr17_43115634_43115878_minus_sgRNA_ex2_meta.csv",
    "chr17_43106355_43106599_minus_sgRNA_ex3_meta.csv",
    "chr17_43104794_43105038_minus_sgRNA_ex4_meta.csv",
    "chr17_43104080_43104330_minus_sgRNA_ex5_meta.csv",
]


def fetch():
    for f in META:
        if not os.path.exists(f):
            print(f"  downloading {f}")
            urllib.request.urlretrieve(f"{RAW}/{f}", f)


def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def valiant_codon_preference():
    """VaLiAnT picks CGG for arginine; PoolParty's table puts AGA first.

    That single entry is the entire difference between 95% and 100% agreement.
    """
    code = copy.deepcopy(dict(G))
    code["R"] = ["CGG"] + [c for c in code["R"] if c != "CGG"]
    return code


def reproduce(meta_csv, use_valiant_code):
    rows = list(csv.DictReader(open(meta_csv)))
    pam = rows[0]["pam_seq"]          # variants build on the PAM-protected reference
    length = len(pam)
    ref_start = int(rows[0]["ref_start"])
    sense = rc(pam)                   # minus strand, so the CDS reads 5'->3' here

    # plus-strand codon start i  ->  sense-strand start  length - 3 - i
    starts = sorted({length - 3 - (int(x["mut_position"]) - ref_start)
                     for x in rows if x["mutator"] == "aa"})
    a, b = starts[0], starts[-1] + 3

    pp.init()
    if use_valiant_code:
        pp.set_genetic_code(valiant_codon_preference())

    p = pp.annotate_orf(pp.from_seq(sense), region_name="cds", extent=(a, b), frame=1)

    # VaLiAnT `aa`: 19 substitutions per codon. Frame comes from the OrfRegion.
    aa = pp.mutagenize_orf(p, region="cds", num_mutations=1,
                           mutation_type="missense_only_first",
                           mode="sequential", prefix="aa")

    # VaLiAnT `inframe`: delete each complete codon. The slice step of 3 supplies
    # codon alignment -- deletion_scan does not read the frame from the OrfRegion.
    inframe = pp.deletion_scan(p, deletion_length=3, deletion_marker=None,
                               region="cds", positions=slice(0, None, 3),
                               mode="sequential", prefix="del")

    lib = pp.stack([aa, inframe])
    seqs = list(lib.to_df(num_cycles=1, show_progress=False)["seq"])
    ours, theirs = set(seqs), {x["mseq_no_adapt"] for x in rows}

    return dict(exon=meta_csv.split("_sgRNA_")[1][:3], tgt_len=length,
                codons=len(starts), extent=(a, b), aa=aa.num_states,
                dele=inframe.num_states, states=lib.num_states,
                distinct=len(ours), theirs=len(theirs),
                match=len(ours & theirs), only_ours=len(ours - theirs),
                only_theirs=len(theirs - ours))


def main():
    fetch()
    files = sorted(META, key=lambda f: f.split("_sgRNA_")[1])
    for label, flag in (("PoolParty default codon table", False),
                        ("VaLiAnT codon preference", True)):
        print(f"\n=== {label} ===")
        print(f"{'exon':5}{'tgt':>5}{'codons':>7}{'extent':>12}{'aa':>6}{'del':>5}"
              f"{'states':>7}{'distinct':>9}{'theirs':>7}{'exact':>7}{'%':>7}")
        tot_match = tot_theirs = 0
        for f in files:
            r = reproduce(f, flag)
            tot_match += r["match"]
            tot_theirs += r["theirs"]
            print(f"{r['exon']:5}{r['tgt_len']:>5}{r['codons']:>7}{str(r['extent']):>12}"
                  f"{r['aa']:>6}{r['dele']:>5}{r['states']:>7}{r['distinct']:>9}"
                  f"{r['theirs']:>7}{r['match']:>7}{100*r['match']/r['theirs']:>6.1f}%")
        print(f"{'ALL':5}{'':5}{'':7}{'':12}{'':6}{'':5}{'':7}{'':9}"
              f"{tot_theirs:>7}{tot_match:>7}{100*tot_match/tot_theirs:>6.1f}%")


if __name__ == "__main__":
    main()
