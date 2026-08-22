# VaLiAnT `brca1_pep` — exact reproduction in PoolParty

**All 2,339 oligos across four targetons reproduced byte-for-byte, from VaLiAnT's
published inputs.** Nothing is taken from its output except the sequences being
checked against.

Answers referee **R2 2b** (overlap of pool elements, and investigation of the
unique elements) and contributes to **R1 #5**.

**VaLiAnT was never installed.** Its expected outputs are committed to its own
repository and validated there by `md5sum`, so the comparison is against the
authors' published ground truth.

Verified against `poolparty` at `7e8c508` (merged `main`), which is the first
version carrying `deletion_scan_orf`.

---

## 1. The example, and why this one

**`brca1_pep`** — an amino-acid-level saturation genome editing library over BRCA1
exons 2–5, on the minus strand of chr17.

VaLiAnT ships **five** example libraries. This one uses only **two mutators**, both
of which PoolParty expresses as named operations, and a **single deletion length**.
It is a fair target rather than a flattering one: `aa` and `inframe` are VaLiAnT's
own primitives, chosen by its authors to demonstrate their tool.

> **Note on the other four.** They were originally out of reach partly because
> `deletion_scan` could not be used twice with different deletion lengths in one
> Party, and four of the five combine `1del` with `2del0`/`2del1`. That collision
> was fixed in the same merge that added `deletion_scan_orf`, so `brca1_nuc` and
> `cdna` may now be reachable. `cdna` additionally needs `snvre`, which PoolParty
> cannot enumerate exhaustively — see `../comparison/FINDINGS.md` C1.

---

## 2. Exact source

Repository: `https://github.com/cancerit/VaLiAnT`, branch **`develop`**
Directory: **`examples/sge/brca1/`**

| Path | Role |
|---|---|
| `run_brca1_pep.sh` | the invocation |
| `parameter_input_files/brca1_pep_targeton_input.txt` | four targeton rows |
| `parameter_input_files/brca1_protection_edits.vcf` | PAM protection, `--pam` |
| `reference_input_files/ENST00000357654.9.gtf` | transcript model, `--gff` |
| `../ref/chr17.fa.gz` | reference, 26 MB gzipped → 85 MB (`unpack_reference.sh`) |
| **`brca1_pep_output_exp/`** | **expected output, v4.0.0 — the ground truth** |
| `check_brca1_pep.sh` | validation by `md5sum`, via `examples/compare_exp_results.sh` |

`brca1_clinvar.vcf`, `brca1_gnomad.vcf` and the custom-variant manifest sit in the
same directory but belong to `brca1_nuc`. **`brca1_pep` passes no `--vcf` and no
`--bg`**, which is why its counts are clean multiples of 20 with no remainder.

**There is no wiki page** for this example; the repository is its only documentation.

VaLiAnT's own command, verbatim:

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

The targeton file — all four rows carry the **same** action vector:

```
ref_chr ref_strand ref_start ref_end   r2_start  r2_end    ext_vector action_vector         sgrna_vector
chr17   -          43115634  43115878  43115726  43115779  "25,25"    "(),(aa,inframe),()"  sgRNA_ex2
chr17   -          43106355  43106599  43106456  43106533  "25,25"    "(),(aa,inframe),()"  sgRNA_ex3
chr17   -          43104794  43105038  43104868  43104956  "25,25"    "(),(aa,inframe),()"  sgRNA_ex4
chr17   -          43104080  43104330  43104122  43104261  "20,41"    "(),(aa,inframe),()"  sgRNA_ex5
```

`ref_start`/`ref_end` bound the whole targeton, which is the oligo.
`r2_start`/`r2_end` locate the exon inside it. `ext_vector` sizes the two flanking
regions r1 and r3 — **in genomic order, not transcript order**, settled by exon 5
where the two values differ (`"20,41"`, and the genomically-left flank is 20 nt).
The action vector's three slots map to r1, r2, r3; only the middle is filled, so
only the exon is mutated and the flanks ride along unchanged.

---

## 3. Mutation types, and how the total is computed

VaLiAnT's mutators are a **closed enum of seven** (`src/valiant/mutator_type.py`).
This example uses two, both requiring coding-sequence annotation:

| Mutator | Definition | Variants per codon |
|---|---|---|
| **`aa`** | substitute the codon to **every other amino acid** (20 non-stop − 1) | **19** |
| **`inframe`** | delete the **complete reading-frame triplet**, preserving frame | **1** |
| | | **20 per codon** |

### The GTF is what fixes the codon count

`r2_start`/`r2_end` are **user-chosen** and need not be codon-aligned — and here
**three of the four are not, at either end**. The reading frame comes from the GTF's
CDS `phase` column, independently.

For a minus-strand targeton, phase counts from the transcript 5′ end, which is r2's
**genomic right** edge:

```
r2 = 43115726 ─────────────────────────────────────────── 43115779   (54 nt)

 726 727 │ 728 ─────────────────────────────────── 778 │ 779
 └──┬──┘   └───────────────┬──────────────────────────┘   └┬┘
    │                      │                               │
 2 orphan bases      17 complete codons (51 nt)      1 orphan base
                                                     = phase

            54  =  1  +  51  +  2
```

Both orphan groups are halves of codons **split across exon boundaries** — 43115779
completes a codon that began in exon 1; 43115726–43115727 begin one that finishes in
exon 3. VaLiAnT substitutes an amino acid only where it holds the whole codon, so
both are skipped. Verified: zero orphan bases carry an `aa` mutation, in all four
targetons.

```
n_codons = (r2_length − phase) // 3
orphans  = (r2_length − phase) %  3      at the genomic left
variants = 20 × n_codons
```

| Exon | Targeton | r2 span | r2 nt | Phase | Codons | `aa` 19n | `inframe` n | **Total** | Distinct |
|---|---|---|---|---|---|---|---|---|---|
| 2 | chr17:43115634-43115878 (245 nt) | 43115726-43115779 | 54 | 1 | **17** | 323 | 17 | **340** | 340 |
| 3 | chr17:43106355-43106599 (245 nt) | 43106456-43106533 | 78 | 1 | **25** | 475 | 25 | **500** | 500 |
| 4 | chr17:43104794-43105038 (245 nt) | 43104868-43104956 | 89 | 1 | **29** | 551 | 29 | **580** | 580 |
| 5 | chr17:43104080-43104330 (251 nt) | 43104122-43104261 | 140 | 2 | **46** | 874 | 46 | **920** | **919** |
| | | | | | **117** | 2,223 | 117 | **2,340** | **2,339** |

Two things this settles that the inputs alone do not:

- **Exon 5 generates 920 variants but 919 distinct sequences.** Two `inframe`
  deletions of adjacent identical codons yield the same oligo. VaLiAnT's `_meta.csv`
  has 920 rows and 919 distinct `mseq_no_adapt`; its `_unique.csv` reports 919.
  PoolParty independently gives `num_states = 920` and 919 distinct — exact
  agreement at **both** levels. (Compare `FINDINGS.md` C6: `num_states` counts state
  slots, not distinct sequences.)
- **The paper and the shipped output disagree by one.** Supp. Table 3 gives 918 for
  exon 5; the committed v4.0.0 output gives 919. Validate against the shipped file,
  which is md5-checkable.

---

## 4. Oligo geometry

Verified from the metadata, not assumed. Getting any of this wrong silently produces
a plausible but wrong library.

| Fact | Value |
|---|---|
| Strand | **minus** (`ref_strand = -`, `revc = 1`) |
| Variants are built on | **`pam_seq`**, *not* `ref_seq` |
| PAM protection (exon 2) | 43115764 C→T, 43115770 C→T, both annotated `syn;syn` |
| Sense strand | `from_fasta(..., strand='-')` — reverse-complemented on extraction |
| Coordinate map | sense index `i` ↔ genomic `ref_end − i` |
| Final oligo | `adaptor5 + sense-strand mutated targeton + adaptor3` |
| `ref`/`new` columns | **plus-strand** codons; `ref_aa`/`alt_aa` are the minus-strand translation |

The strand trap, concretely: exon 2's second PAM edit reads `CAA → TAA` on the plus
strand, which looks like a stop codon being created. Read in the correct orientation
it is `TTG → TTA`, Leu → Leu. Both edits are wobble-position changes.

Sanity check performed: the exon 2 sense-strand ORF at `[100, 151)` translates to
`LELIKEPVSTKCDHIFC`, equal to VaLiAnT's `ref_aa` column read in minus-strand order.

---

## 5. PoolParty reproduction

`repro_valiant_brca1_pep.py` in this directory is the executable version. Core:

```python
# PoolParty extracts the targeton and reverse-complements it
tgt   = pp.from_fasta("chr17.fa", ("chr17", ref_start - 1, ref_end, "-"))
sense = tgt.to_df(num_cycles=1, show_progress=False)["seq"][0].upper()

# PAM protection edits from the published VCF, mapped into sense coordinates
for pos, (ref, alt) in PAM[ex].items():
    i = ref_end - pos
    assert chars[i] == ref.translate(COMP)      # fail loudly on a coordinate error
    chars[i] = alt.translate(COMP)

# Region and frame from the GTF's CDS phase -- the user's job, see section 6
a, b, n_codons = orf_extent(ref_end, r2_start, r2_end, phase)

p       = pp.annotate_orf(pp.from_seq(sense), region_name="cds", extent=(a, b), frame=1)
aa      = pp.mutagenize_orf(p, region="cds", num_mutations=1,
                            mutation_type="missense_only_first",
                            mode="sequential", prefix="aa")          # VaLiAnT `aa`
inframe = pp.deletion_scan_orf(p, deletion_codons=1, deletion_marker=None,
                               region="cds", mode="sequential", prefix="del")
                                                                     # VaLiAnT `inframe`
lib     = pp.stack([aa, inframe])
df      = lib.to_df(num_cycles=1, show_progress=False)
```

Four operations: one Pool, fan-out to two mutagenesis operations, merged by `stack`
— the same shape as the manuscript's Figure 2.

`deletion_scan_orf` takes `deletion_codons` and `codon_positions` in **codon units**,
so no offset arithmetic is needed and the scan cannot run past the ORF. It arrived in
the merge at `1a29b22`; before that, the same result required
`deletion_scan(deletion_length=3, region=..., positions=slice(0, None, 3))`, which is
correct only with `region=` supplied — see `../comparison/FINDINGS.md` B1.

`to_df` is used rather than `generate_library`, because `generate_library` leaves
region tags in the sequence; `to_df`'s `write_tags` defaults to `False`.

### Results

| Exon | Phase | Codons | ORF extent (sense) | PoolParty states | Distinct | VaLiAnT distinct | Exact match |
|---|---|---|---|---|---|---|---|
| 2 | 1 | 17 | `(100, 151)` | 340 | 340 | 340 | **340 / 340** |
| 3 | 1 | 25 | `(67, 142)` | 500 | 500 | 500 | **500 / 500** |
| 4 | 1 | 29 | `(83, 170)` | 580 | 580 | 580 | **580 / 580** |
| 5 | 2 | 46 | `(71, 209)` | 920 | 919 | 919 | **919 / 919** |
| **All** | | **117** | | 2,340 | **2,339** | **2,339** | **2,339 / 2,339 (100%)** |

The extents are computed by `orf_extent()` from the targeton row and the GTF phase
alone — no read of VaLiAnT's output. They agree with the codon positions recorded in
its `_meta.csv`, which is the check that the phase handling is right.

Comparison is set-wise on `mseq_no_adapt` — the mutated targeton without adaptors.
Zero sequences unique to either tool.

### The one discrepancy

Under PoolParty's default codon table, agreement is 95.2%, and **every** mismatch is
the **arginine** substitution — exactly one per codon, 117 of them.

| | Arg codon chosen |
|---|---|
| PoolParty `missense_only_first` | **AGA** — first in its Kazusa-ordered list |
| VaLiAnT | **CGG** |

`missense_only_first` means "the most frequent codon", and the two tools disagree on
which codon that is, for exactly one residue. Nothing about the *design space*
differed — only the DNA realisation of one amino acid. Reordering that one entry via
`pp.set_genetic_code()` gives 100%.

This is the substantive answer to R2 2b's "investigate the unique elements": the
unique elements were a codon-usage parameter, and PoolParty exposes the parameter.

---

## 6. What this does and does not show

**Division of labour, stated plainly.**

| Step | Who |
|---|---|
| Extract the targeton from `chr17.fa`, reverse-complemented | **PoolParty** (`from_fasta`, `strand='-'`) |
| Apply the PAM protection VCF | harness |
| Read the GTF's CDS phase; choose region and frame | **harness** |
| Generate the variants | **PoolParty** (`mutagenize_orf`, `deletion_scan_orf`, `stack`) |

**PoolParty is not annotation-aware.** It does not read GTF/GFF files and has no
transcript model, so the target region and reading frame are supplied by the user.
Here they were derived from the GTF's CDS `phase` column by the harness. That is the
`transcript_annotation_aware` = `no` entry in Table 2, made concrete — and it is the
only step carried by anything other than PoolParty or VaLiAnT.

**Shows.** Given the same reference and the same codon partition, PoolParty
reproduces a published VaLiAnT library exactly — all 2,339 oligos, including the
flanking extensions, the PAM-protection edits, the minus-strand orientation and the
deletion-length bookkeeping. Reference extraction is PoolParty's own.

**Does not show.** Resolution of the reading frame from an annotation file. VaLiAnT
reads the GTF and derives phase, strand and codon partition; the harness does that
here.

**Does not exercise VaLiAnT's split-codon handling.** Its documented
"codons split across exons" capability is not tested by this example: every split
codon at an r2 boundary is **skipped**, not resolved. Its cDNA mode, where the CDS is
contiguous, is where that feature would bite.

**Not compared.** Adaptors and the min/max length filter — the comparison is on
`mseq_no_adapt`. Dedup *is* compared, and agrees (exon 5, 920 → 919).

**Not attempted.** `snvre`, exhaustive synonymous re-encoding, which PoolParty cannot
enumerate because the per-residue alternative count is non-uniform and the lazy state
arithmetic requires a fixed branching factor (`FINDINGS.md` C1).

---

## 7. Status of this series

| Tool | Target | Result |
|---|---|---|
| **VaLiAnT** | `brca1_pep`, BRCA1 exons 2–5 SGE | **2,339 / 2,339 exact** |
| **MPRAnator** | MPRA Motif Use Case | **976 / 976 exact** at the arrangement level — `mpranator_motif_usecase.md` |
| **tangermeme** | Tutorial B3 — motif-pair spacing sweeps | not started; pip-installable, so runnable head-to-head |

The mirror direction: **VaLiAnT expresses one of the four components of our GB1
library** — all single amino-acid substitutions, 1,045 of 547,230 sequences, 0.19%.
It has no replication operation (its pipeline ends in order-invariant dedup), no
combination mutator, and no stochastic mutator. **MPRAnator** cannot express the GB1
library at all, having no codon or amino-acid concept.
