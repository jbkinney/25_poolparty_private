# VaLiAnT `brca1_pep` — exact reproduction in PoolParty

**Result: 2,339 of 2,339 oligos reproduced byte-for-byte across all four targetons,
once PoolParty is given VaLiAnT's codon preference.** With PoolParty's default codon
table the agreement is 95%; the entire 5% is one amino acid (see
[The one discrepancy](#the-one-discrepancy)).

Answers referee **R2 2b** (overlap of pool elements, and investigation of the unique
elements) and contributes to **R1 #5** (a common use case expressed in both tools).

**VaLiAnT was never installed.** Its expected outputs are committed to its own
repository and validated there by `md5sum`, so the comparison is against the
authors' published ground truth rather than against our own run of their tool.

---

## 1. The example, and why this one

**`brca1_pep`** — an amino-acid-level saturation genome editing library over BRCA1
exons 2–5, on the minus strand of chr17.

VaLiAnT ships **five** example libraries. This one was chosen on three grounds:

| Ground | Detail |
|---|---|
| **Only two mutators** | `(aa, inframe)`. Both are expressible in PoolParty; four of the other examples additionally use `snvre`, which PoolParty cannot enumerate exhaustively (see `../comparison/FINDINGS.md` C1) |
| **One deletion length** | `deletion_scan` hard-codes its internal region name, so two different `deletion_length` values collide in one Party (`FINDINGS.md` B3). The other examples combine `1del` with `2del0`/`2del1` |
| **Complete ground truth** | Committed expected output including the reference sequences, so no reference genome, GTF or install is required |

It is a **fair** target rather than a flattering one: `aa` and `inframe` are
VaLiAnT's own primitives, chosen by its authors to demonstrate their tool.

---

## 2. Exact source

Repository: `https://github.com/cancerit/VaLiAnT`, branch **`develop`**
Directory: **`examples/sge/brca1/`**

| Path | Role |
|---|---|
| `run_brca1_pep.sh` | the invocation (reproduced below) |
| `parameter_input_files/brca1_pep_targeton_input.txt` | the four targeton rows |
| `parameter_input_files/brca1_protection_edits.vcf` | PAM/protospacer protection, passed as `--pam` |
| `reference_input_files/ENST00000357654.9.gtf` | transcript model, passed as `--gff` |
| **`brca1_pep_output_exp/`** | **expected output, v4.0.0 — the ground truth used here** |
| `brca1_pep_output_exp_v3/` | expected output, v3 era |
| `check_brca1_pep.sh` | validation, via `examples/compare_exp_results.sh`, by `md5sum` |

Files consumed for this comparison, one per targeton:

```
brca1_pep_output_exp/chr17_43115634_43115878_minus_sgRNA_ex2_meta.csv
brca1_pep_output_exp/chr17_43106355_43106599_minus_sgRNA_ex3_meta.csv
brca1_pep_output_exp/chr17_43104794_43105038_minus_sgRNA_ex4_meta.csv
brca1_pep_output_exp/chr17_43104080_43104330_minus_sgRNA_ex5_meta.csv
brca1_pep_output_exp/ref_sequences.csv
```

There is **no wiki page** for this example. The `develop` repository is its only
documentation. (The wiki documents *Input files*, *Advanced usage*, *cDNA example*,
*Saturation-prime-editing-example* and *Docker-usage-example* only.)

VaLiAnT's own command, verbatim from `run_brca1_pep.sh`:

```bash
valiant sge \
    parameter_input_files/brca1_pep_targeton_input.txt \
    ../ref/chr17.fa  brca1_pep_output  'homo sapiens' 'GRCh38' \
    --pam parameter_input_files/brca1_protection_edits.vcf \
    --revcomp-minus-strand \
    --adaptor-5 AATGATACGGCGACCACCGA \
    --adaptor-3 TCGTATGCCGTCTTCTGCTTG \
    --gff reference_input_files/ENST00000357654.9.gtf
```

The targeton input, verbatim — all four rows carry the **same** action vector:

```
ref_chr ref_strand ref_start ref_end   r2_start  r2_end    ext_vector action_vector         sgrna_vector
chr17   -          43115634  43115878  43115726  43115779  "25,25"    "(),(aa,inframe),()"  sgRNA_ex2
chr17   -          43106355  43106599  43106456  43106533  "25,25"    "(),(aa,inframe),()"  sgRNA_ex3
chr17   -          43104794  43105038  43104868  43104956  "25,25"    "(),(aa,inframe),()"  sgRNA_ex4
chr17   -          43104080  43104330  43104122  43104261  "20,41"    "(),(aa,inframe),()"  sgRNA_ex5
```

---

## 3. Mutation types

VaLiAnT's mutators are a **closed enum of seven** (`src/valiant/mutator_type.py`).
This example uses two, both of which require coding-sequence annotation:

| Mutator | Definition | Variants per codon |
|---|---|---|
| **`aa`** | Substitute the codon to **every other amino acid**: 20 non-stop residues minus the current one | **19** |
| **`inframe`** | Delete the **complete reading-frame triplet**, preserving frame. Output is 3 nt shorter | **1** |
| | | **20 per codon** |

The action vector `"(),(aa,inframe),()"` has three slots — **r1, r2, r3**. Only the
central target region **r2** is mutated; the flanking intronic extensions (25/25 nt,
or 20/41 for exon 5) and the P5/P7 adaptors are carried through unchanged.

### Computation of the totals

Every oligo carries exactly one mutation, so the count is
`19 x codons` substitutions plus `1 x codons` deletions:

| Exon | Targeton | r2 span | r2 nt | Frame offset in r2 | In-frame codons | `aa` = 19n | `inframe` = n | **Total 20n** | Distinct |
|---|---|---|---|---|---|---|---|---|---|
| 2 | chr17:43115634-43115878 (245 nt) | 43115726-43115779 | 54 | 2 | **17** | 323 | 17 | **340** | 340 |
| 3 | chr17:43106355-43106599 (245 nt) | 43106456-43106533 | 78 | 2 | **25** | 475 | 25 | **500** | 500 |
| 4 | chr17:43104794-43105038 (245 nt) | 43104868-43104956 | 89 | 1 | **29** | 551 | 29 | **580** | 580 |
| 5 | chr17:43104080-43104330 (251 nt) | 43104122-43104261 | 140 | 0 | **46** | 874 | 46 | **920** | **919** |
| | | | | | **117** | 2,223 | 117 | **2,340** | **2,339** |

Three things this table settles, none of them guessable from the inputs alone:

- **r2 is not codon-phased.** Exons 2 and 3 have r2 lengths that are exact multiples
  of 3 (54, 78) yet contain only 17 and 25 complete codons, because the reading frame
  starts 2 nt into r2. The frame offsets are 2, 2, 1, 0 and come from the transcript
  model, not from r2's boundaries.
- **Exon 5 generates 920 variants but 919 distinct sequences.** Two `inframe`
  deletions of adjacent identical codons yield the same oligo. VaLiAnT's `_meta.csv`
  has 920 rows and 919 distinct `mseq_no_adapt`; its `_unique.csv` reports 919.
  PoolParty independently gives `num_states = 920` and 919 distinct sequences — exact
  agreement at **both** levels. (Compare `FINDINGS.md` C6: `num_states` counts state
  slots, not distinct sequences.)
- **The published paper and the shipped output disagree by one.** Supp. Table 3 gives
  918 for exon 5; the committed v4.0.0 output gives 919. Validate against the shipped
  file, which is md5-checkable.

---

## 4. Oligo geometry

Verified from the metadata, not assumed. Getting any of this wrong silently produces
a plausible but wrong library.

| Fact | Value |
|---|---|
| Strand | **minus** (`ref_strand = -`, `revc = 1`) |
| Template variants are built on | **`pam_seq`**, *not* `ref_seq` |
| PAM protection (exon 2) | 2 edits: 43115764 C→T, 43115770 C→T, both annotated `syn;syn` |
| Sense strand (CDS reads 5'->3') | **`reverse_complement(pam_seq)`** |
| Coordinate map, plus index `i` -> sense index | **`L - 3 - i`**, where `L` = targeton length |
| Final oligo | `adaptor5 + reverse_complement(mutated targeton) + adaptor3` |
| `ref`/`new` columns | **plus-strand** codons; `ref_aa`/`alt_aa` are the minus-strand translation |

Sanity check performed: the exon 2 sense-strand ORF at `[100, 151)` translates to
`LELIKEPVSTKCDHIFC`, which equals VaLiAnT's `ref_aa` column read in minus-strand
order.

Resulting sense-strand ORF extents: exon 2 `[100,151)`, exon 3 `[67,142)`,
exon 4 `[83,170)`, exon 5 `[71,209)`.

---

## 5. PoolParty reproduction

Run against `poolparty 0.1.0`. `repro_valiant_brca1_pep.py` in this directory is the
executable version; the core is:

```python
import copy, csv
import poolparty as pp
from poolparty.codon_table import STANDARD_GENETIC_CODE as G

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

rows = list(csv.DictReader(open(meta_csv)))       # VaLiAnT's committed _meta.csv
pam  = rows[0]["pam_seq"]                         # variants build on the PAM-protected ref
L    = len(pam)
sense = rc(pam)                                   # minus strand -> CDS reads 5'->3'

# plus-strand codon start i  ->  sense-strand start  L - 3 - i
starts = sorted({L - 3 - (int(x["mut_position"]) - int(rows[0]["ref_start"]))
                 for x in rows if x["mutator"] == "aa"})
a, b = starts[0], starts[-1] + 3                  # the in-frame ORF extent

pp.init()
pp.set_genetic_code(valiant_codon_preference())   # see section 6

p = pp.annotate_orf(pp.from_seq(sense), region_name="cds", extent=(a, b), frame=1)

# VaLiAnT `aa` -- 19 substitutions per codon.
# mutagenize_orf reads the reading frame from the OrfRegion registered above.
aa = pp.mutagenize_orf(p, region="cds", num_mutations=1,
                       mutation_type="missense_only_first",
                       mode="sequential", prefix="aa")

# VaLiAnT `inframe` -- delete each complete codon.
# The slice step of 3 gives codon-aligned starts; positions are relative to the
# region, which begins at the ORF start. deletion_scan does NOT read the frame
# from the OrfRegion -- the slice is doing that work. See FINDINGS.md B1.
inframe = pp.deletion_scan(p, deletion_length=3, deletion_marker=None,
                           region="cds", positions=slice(0, None, 3),
                           mode="sequential", prefix="del")

lib = pp.stack([aa, inframe])
df  = lib.to_df(num_cycles=1, show_progress=False)   # write_tags=False by default
```

Two undocumented details this depends on:

- **`positions` accepts a `slice`, and the step is honoured.** `deletion_scan.rst`
  describes it as an *"Explicit list of window start positions"* and never mentions
  slices. The shared `validate_positions` resolves them via `positions.indices(...)`.
  See `FINDINGS.md` B1.
- **`generate_library()` leaves region tags in the sequence** (`<cds>...</cds>`). Use
  `to_df(num_cycles=1)`, whose `write_tags` defaults to `False`.

---

## 6. Results

| Exon | VaLiAnT distinct | PoolParty states | Default codon table | **VaLiAnT codon preference** |
|---|---|---|---|---|
| 2 | 340 | 340 | 323 (95.0%) | **340 / 340 (100%)** |
| 3 | 500 | 500 | 475 (95.0%) | **500 / 500 (100%)** |
| 4 | 580 | 580 | 552 (95.2%) | **580 / 580 (100%)** |
| 5 | 919 | 920 (919 distinct) | 876 (95.3%) | **919 / 919 (100%)** |
| **All** | **2,339** | **2,340 (2,339 distinct)** | 95.2% | **2,339 / 2,339 (100%)** |

Comparison is set-wise on `mseq_no_adapt` — the mutated targeton without adaptors.
Zero sequences unique to either tool.

### The one discrepancy

Every mismatch under PoolParty's default table was the **arginine** substitution —
exactly one per codon, 117 of them.

| | Arg codon chosen |
|---|---|
| PoolParty, `missense_only_first` | **AGA** (first in its Kazusa-ordered list `['AGA','AGG','CGG','CGC','CGA','CGT']`) |
| VaLiAnT | **CGG** |

`missense_only_first` means "the most frequent codon", and the two tools disagree on
which codon that is for exactly one residue. Nothing about the *design space*
differed — only the DNA realisation of one amino acid. Reordering that one entry via
`pp.set_genetic_code()` gives 100%:

```python
def valiant_codon_preference():
    code = copy.deepcopy(dict(G))
    code["R"] = ["CGG"] + [c for c in code["R"] if c != "CGG"]
    return code
```

This is the substantive answer to R2 2b's "investigate the unique elements": the
unique elements were a codon-usage parameter, not a capability difference, and
PoolParty exposes the parameter.

---

## 7. What this does and does not show

**Shows.** Given the same target sequence and the same codon preference, PoolParty
reproduces a published VaLiAnT library exactly — all 2,339 oligos across four
targetons, including the flanking extensions, the PAM-protection edits, the
minus-strand orientation and the deletion-length bookkeeping.

**Does not show.** PoolParty did not resolve the targetons from genomic coordinates.
VaLiAnT reads chr17 plus a GTF and derives the reference sequence, reading frame and
codon phase; here those came from VaLiAnT's own `_meta.csv`/`ref_sequences.csv`.
That is the genomic-integration gap recorded in Table 2 (`genome_coordinates`
partial, `transcript_annotation_aware` no) and it is real. The claim is about
**variant generation given a target sequence**, not about coordinate resolution.

**Not attempted here.** `snvre` (exhaustive synonymous re-encoding) and mixed
deletion spans, for the reasons in section 1.

---

## 8. Planned additions

Same directory, same shape — source, mutation accounting, code, results:

| Tool | Target | Status |
|---|---|---|
| **MPRAnator** | "MPRA Motif Use Case Example" — AP1/ELK1/RFX over two mm9 backgrounds, **5,856 sequences** stated | not started. The live service returns **HTTP 500** on barcode-enabled requests, so the comparison must be against the documented count, not a run |
| **tangermeme** | Tutorial **B3 — Spacing**: motif-pair spacing sweeps, which its own docs call "a small combinatorial motif-placement library" | not started. Pip-installable and actively maintained, so this one can be executed head-to-head |

Those three tools are mutually near-disjoint in Table 2 — none can express another's
flagship example — and they map onto the paper's three declared applications: DMS,
MPRA, and in-silico probing of genomic models.

The mirror direction belongs with them: **VaLiAnT expresses one of the four
components of our GB1 library** (all single amino-acid substitutions, 1,045 of
547,230 sequences — 0.19%). It has no replication operation (its pipeline ends in
order-invariant dedup), no combination mutator, and no stochastic mutator.
