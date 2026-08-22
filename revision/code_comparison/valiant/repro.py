"""Express one library in both VaLiAnT and PoolParty, and diff the oligos.

This is NOT a reproduction of a VaLiAnT shipped example. The library is one we
defined -- their brca1_pep targeton with `stop` added to the action vector and PAM
protection dropped -- generated with both tools and compared.

Provenance: the same VaLiAnT installation reproduces their shipped brca1_pep
library exactly (2,339 of 2,339 oligos) and passes their own check_brca1_pep.sh
md5 validation, so the ground truth here is trustworthy despite being self-generated.

Requires poolparty >= 1a29b22 on main (deletion_scan_orf, insertion_scan_orf).
See comparison.md.

    python3 repro_valiant_brca1_modified.py
"""

import copy
import csv
import os
import sys

import poolparty as pp
from poolparty.codon_table import STANDARD_GENETIC_CODE as G

# VaLiAnT's clone is a sibling of this repository. Large inputs and VaLiAnT's own
# output live there, never in this repository. Override with VALIANT_EXAMPLE_DATA.
_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT = os.path.normpath(os.path.join(_HERE, "..", "..", "..", "..",
                                         "VaLiAnT", "examples", "sge"))
DATA_DIR = os.path.abspath(os.environ.get("VALIANT_EXAMPLE_DATA", _DEFAULT))
FASTA = os.path.join(DATA_DIR, "ref", "chr17.fa")
VALIANT_OUT = os.path.join(DATA_DIR, "brca1", "brca1_modified_output",
                           "chr17_43115634_43115878_minus_meta.csv")

# The single targeton, from inputs/targeton.tsv, plus the GTF CDS phase for that
# span. Both are published inputs; nothing here comes from VaLiAnT's output.
REF_START, REF_END = 43115634, 43115878      # 1-based inclusive
R2_START, R2_END = 43115726, 43115779
PHASE = 1                                    # GTF CDS phase for chr17:43115726-43115779

ADAPTOR_5 = "AATGATACGGCGACCACCGA"
ADAPTOR_3 = "TCGTATGCCGTCTTCTGCTTG"

VALIANT_CMD = """\
cd {clone}/examples/sge/brca1
valiant sge parameter_input_files/targeton.tsv \\
    ../ref/chr17.fa  brca1_modified_output  'homo sapiens' 'GRCh38' \\
    --revcomp-minus-strand \\
    --adaptor-5 {a5} --adaptor-3 {a3} \\
    --gff reference_input_files/ENST00000357654.9.gtf"""


def valiant_codon_preference():
    """VaLiAnT picks CGG for arginine where PoolParty's table puts AGA first.

    That single entry is the whole difference between 95% and 100% agreement.
    """
    code = copy.deepcopy(dict(G))
    code["R"] = ["CGG"] + [c for c in code["R"] if c != "CGG"]
    return code


def orf_extent():
    """Sense-strand ORF extent, from the targeton row and the GTF phase alone.

    from_fasta(strand='-') returns the sense strand, so sense index i corresponds
    to genomic (REF_END - i). GTF phase counts from the transcript 5' end, which on
    the minus strand is r2's GENOMIC RIGHT edge. Those `phase` bases belong to a
    codon that began in the previous exon; the remainder at the other end belongs to
    one that continues into the next. Both are skipped -- an amino acid can only be
    substituted where the whole codon is present.
    """
    n_codons = ((R2_END - R2_START + 1) - PHASE) // 3
    start = (REF_END - R2_END) + PHASE
    return start, start + 3 * n_codons, n_codons


def build_library():
    """The figure code. Five statements, one fan-out, no loops."""
    orf_start, orf_stop, _ = orf_extent()

    orf = pp.from_fasta(FASTA, ("chr17", REF_START - 1, REF_END, "-")) \
            .annotate_orf(region_name="cds", extent=(orf_start, orf_stop), frame=1)

    aa = orf.mutagenize_orf(region="cds", num_mutations=1,
                            mutation_type="missense_only_first", mode="sequential")
    inframe = orf.deletion_scan_orf(deletion_codons=1, deletion_marker=None,
                                    region="cds", mode="sequential")
    stop = orf.insertion_scan_orf("TGA", region="cds", replace=True, mode="sequential")

    return pp.join([ADAPTOR_5, pp.stack([aa, inframe, stop]), ADAPTOR_3]), aa, inframe, stop


def main():
    if not os.path.exists(FASTA):
        sys.exit(f"missing reference: {FASTA}\n"
                 f"run {DATA_DIR}/unpack_reference.sh in the VaLiAnT clone")
    if not os.path.exists(VALIANT_OUT):
        sys.exit("missing VaLiAnT output. In the `valiant` conda env, run:\n\n"
                 + VALIANT_CMD.format(clone=os.path.dirname(os.path.dirname(DATA_DIR)),
                                      a5=ADAPTOR_5, a3=ADAPTOR_3) + "\n")

    pp.init()
    pp.set_genetic_code(valiant_codon_preference())

    orf_start, orf_stop, n_codons = orf_extent()
    library, aa, inframe, stop = build_library()

    ours = {s.upper() for s in library.to_df(num_cycles=1, show_progress=False)["seq"]}
    theirs = {r["mseq"].upper() for r in csv.DictReader(open(VALIANT_OUT))}

    print(f"  data dir      : {DATA_DIR}")
    print(f"  codons        : {n_codons}   ORF extent (sense): ({orf_start}, {orf_stop})")
    print(f"  aa            : {aa.num_states}  (19 x {n_codons})")
    print(f"  inframe       : {inframe.num_states}")
    print(f"  stop          : {stop.num_states}")
    print(f"  library states: {library.num_states}   distinct: {len(ours)}")
    print(f"  VaLiAnT oligos: {len(theirs)}")
    print(f"  exact match   : {len(ours & theirs)}/{len(theirs)} "
          f"({100 * len(ours & theirs) / len(theirs):.1f}%)")
    print(f"  only ours     : {len(ours - theirs)}     only theirs: {len(theirs - ours)}")
    print(f"  oligo lengths : {sorted({len(s) for s in ours})}  "
          f"(245 + {len(ADAPTOR_5)} + {len(ADAPTOR_3)} = 286; deletions 283)")


if __name__ == "__main__":
    main()
