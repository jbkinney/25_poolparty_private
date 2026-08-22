"""Reproduce VaLiAnT's brca1_pep library in PoolParty, from published inputs only.

Nothing is read from VaLiAnT's output except the 2,339 sequences being checked
against. The four inputs used are all published in VaLiAnT's own repository:
the targeton TSV, chr17.fa, the transcript GTF, and the PAM protection VCF.

Division of labour, stated plainly because it is the point:

  PoolParty  reference extraction (from_fasta, strand-aware) and variant
             generation (mutagenize_orf + deletion_scan_orf + stack)
  harness    applies the PAM protection VCF, and reads the GTF's CDS phase to
             work out which region and reading frame to mutate

PoolParty is not annotation-aware. Choosing the region and frame is the user's
job, and that is the transcript_annotation_aware concession in Table 2.

Requires poolparty >= the merge of feature/orf-indel-scans (deletion_scan_orf).
See valiant_brca1_pep.md.

    python3 repro_valiant_brca1_pep.py
"""

import copy
import csv
import gzip
import os
import shutil
import urllib.request

import poolparty as pp
from poolparty.codon_table import STANDARD_GENETIC_CODE as G

RAW = "https://raw.githubusercontent.com/cancerit/VaLiAnT/develop/examples/sge"
FASTA = "chr17.fa"

# From parameter_input_files/brca1_pep_targeton_input.txt (1-based inclusive) and
# the GTF CDS phase. All four rows carry action_vector "(),(aa,inframe),()".
#            ref_start  ref_end   r2_start  r2_end   gtf_phase
TARGETON = {"ex2": (43115634, 43115878, 43115726, 43115779, 1),
            "ex3": (43106355, 43106599, 43106456, 43106533, 1),
            "ex4": (43104794, 43105038, 43104868, 43104956, 1),
            "ex5": (43104080, 43104330, 43104122, 43104261, 2)}

# From parameter_input_files/brca1_protection_edits.vcf. Plus-strand REF/ALT.
PAM = {"ex2": {43115764: ("C", "T"), 43115770: ("C", "T")},
       "ex3": {43106509: ("G", "A"), 43106506: ("C", "T")},
       "ex4": {43104953: ("G", "A"), 43104950: ("T", "C")},
       "ex5": {43104167: ("G", "A"), 43104164: ("A", "G")}}

META = {ex: f"chr17_{a}_{b}_minus_sgRNA_{ex}_meta.csv"
        for ex, (a, b, _, _, _) in TARGETON.items()}

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def fetch():
    """Ground truth and the reference. chr17.fa unpacks to 85 MB."""
    for ex, f in META.items():
        if not os.path.exists(f):
            print(f"  downloading {f}")
            urllib.request.urlretrieve(f"{RAW}/brca1/brca1_pep_output_exp/{f}", f)
    if not os.path.exists(FASTA):
        print("  downloading chr17.fa.gz (26 MB, unpacks to 85 MB)")
        urllib.request.urlretrieve(f"{RAW}/ref/chr17.fa.gz", "chr17.fa.gz")
        with gzip.open("chr17.fa.gz", "rb") as fi, open(FASTA, "wb") as fo:
            shutil.copyfileobj(fi, fo)


def valiant_codon_preference():
    """VaLiAnT picks CGG for arginine where PoolParty's table puts AGA first.

    That single entry is the whole difference between 95% and 100% agreement.
    """
    code = copy.deepcopy(dict(G))
    code["R"] = ["CGG"] + [c for c in code["R"] if c != "CGG"]
    return code


def orf_extent(ref_end, r2_start, r2_end, phase):
    """Sense-strand ORF extent, from the targeton row and the GTF phase alone.

    These are minus-strand targetons, so from_fasta(strand='-') returns the
    sense strand and sense index i corresponds to genomic (ref_end - i).

    GTF phase counts from the transcript 5' end, which on the minus strand is
    r2's GENOMIC RIGHT edge. Those `phase` bases belong to a codon that began in
    the previous exon; the leftover at the other end belongs to a codon that
    continues into the next one. Both are skipped -- VaLiAnT substitutes an amino
    acid only where it holds the whole codon.
    """
    n_codons = ((r2_end - r2_start + 1) - phase) // 3
    start = (ref_end - r2_end) + phase
    return start, start + 3 * n_codons, n_codons


def reproduce(ex):
    ref_start, ref_end, r2_start, r2_end, phase = TARGETON[ex]
    a, b, n_codons = orf_extent(ref_end, r2_start, r2_end, phase)

    pp.init()
    pp.set_genetic_code(valiant_codon_preference())

    # PoolParty extracts the targeton and reverse-complements it (strand='-').
    tgt = pp.from_fasta(FASTA, ("chr17", ref_start - 1, ref_end, "-"))
    sense = tgt.to_df(num_cycles=1, show_progress=False)["seq"][0].upper()

    # PAM protection edits, mapped into sense coordinates. The assertion makes a
    # coordinate or strand error fail loudly instead of producing a plausible
    # but wrong library.
    chars = list(sense)
    for pos, (ref, alt) in PAM[ex].items():
        i = ref_end - pos
        assert chars[i] == ref.translate(COMP), f"{ex}: PAM REF mismatch at {pos}"
        chars[i] = alt.translate(COMP)
    sense = "".join(chars)

    # PoolParty generates the variants.
    p = pp.annotate_orf(pp.from_seq(sense), region_name="cds", extent=(a, b), frame=1)
    aa = pp.mutagenize_orf(p, region="cds", num_mutations=1,
                           mutation_type="missense_only_first",
                           mode="sequential", prefix="aa")        # 19 per codon
    inframe = pp.deletion_scan_orf(p, deletion_codons=1, deletion_marker=None,
                                   region="cds", mode="sequential", prefix="del")
    lib = pp.stack([aa, inframe])                                 # 1 per codon

    ours = set(lib.to_df(num_cycles=1, show_progress=False)["seq"])
    theirs = {r["mseq_no_adapt"] for r in csv.DictReader(open(META[ex]))}

    return dict(ex=ex, phase=phase, codons=n_codons, extent=(a, b),
                states=lib.num_states, distinct=len(ours), theirs=len(theirs),
                match=len(ours & theirs), only_ours=len(ours - theirs),
                only_theirs=len(theirs - ours))


def main():
    fetch()
    print(f"\n{'exon':6}{'phase':>6}{'codons':>7}{'extent':>13}{'states':>8}"
          f"{'distinct':>9}{'theirs':>8}{'exact':>7}{'%':>8}")
    tot_match = tot_theirs = 0
    for ex in TARGETON:
        r = reproduce(ex)
        tot_match += r["match"]
        tot_theirs += r["theirs"]
        print(f"  {r['ex']:4}{r['phase']:>6}{r['codons']:>7}{str(r['extent']):>13}"
              f"{r['states']:>8}{r['distinct']:>9}{r['theirs']:>8}{r['match']:>7}"
              f"{100 * r['match'] / r['theirs']:>7.1f}%")
    print(f"  {'ALL':4}{'':6}{'':7}{'':13}{'':8}{'':9}{tot_theirs:>8}{tot_match:>7}"
          f"{100 * tot_match / tot_theirs:>7.1f}%")


if __name__ == "__main__":
    main()
