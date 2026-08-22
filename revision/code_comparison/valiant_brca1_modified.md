# One library, two tools: VaLiAnT and PoolParty

**357 of 357 oligos byte-identical, adaptors included.**

Source for a supplementary figure comparing how the same library is specified in each
tool. Contributes to referee **R1 #5** and **R2 2b**.

## Scope — read this first

**This is not a reproduction of a VaLiAnT shipped example.** The library is one we
defined: their `brca1_pep` targeton for BRCA1 exon 2, with `stop` added to the
action vector and PAM protection dropped. Both tools were given that design and their
outputs compared.

The ground truth is therefore **our own VaLiAnT run**, not a committed artifact. What
makes it trustworthy: the same installation reproduces VaLiAnT's shipped `brca1_pep`
library exactly — **2,339 of 2,339 oligos** — and passes their own
`check_brca1_pep.sh`, which md5-compares unique CSVs, metadata CSVs and VCFs against
committed expected output. All three validated.

Environment: VaLiAnT 4.0.0 (`develop`, conda env `valiant`, Python 3.11);
`poolparty` at `7e8c508` on `main`, the first version carrying `deletion_scan_orf`.

---

## 1. The design

One targeton: BRCA1 exon 2 on the minus strand of chr17, with 25 nt of intronic
flank either side, three mutation types applied to the exon, and P5/P7 adaptors.

| | |
|---|---|
| Targeton | chr17:43115634-43115878 (245 nt), minus strand |
| Target region r2 | chr17:43115726-43115779 (54 nt) — the exon |
| Flanks | 25 nt each side, unmutated |
| Mutators | `aa`, `inframe`, `stop` |
| Adaptors | P5 `AATGATACGGCGACCACCGA` / P7, appended to every oligo |

Deliberately excluded:

- **`ala`** — VaLiAnT's alanine scan is a strict subset of `aa`, which already
  substitutes every other amino acid including alanine, and both tools pick `GCC`.
  Including it adds 17 duplicates and nothing else.
- **PAM protection** — needs no explanation here, and removing it also removes an
  input file and the whole CRISPR digression.
- **`snvre`** — exhaustive synonymous re-encoding, which PoolParty cannot enumerate;
  the per-residue alternative count is non-uniform and the lazy state arithmetic
  needs a fixed branching factor (`../comparison/FINDINGS.md` C1).

---

## 2. Input files — three, and one is 170 bytes

| File | Size | Role |
|---|---|---|
| `inputs/targeton.tsv` | **170 B** | the design. We authored it; versioned here |
| `ENST00000357654.9.gtf` | 25.7 KB | transcript model, `--gff`. Theirs |
| `chr17.fa` | 81 MB | reference. Theirs, via `unpack_reference.sh` |

The targeton file complete:

```
ref_chr	ref_strand	ref_start	ref_end	r2_start	r2_end	ext_vector	action_vector	sgrna_vector
chr17	-	43115634	43115878	43115726	43115779	"25,25"	"(),(aa,inframe,stop),()"	
```

`ref_start`/`ref_end` bound the targeton, which is the oligo. `r2_start`/`r2_end`
locate the exon inside it. `ext_vector` sizes the two flanks, **in genomic order**.
The action vector's three slots map to r1, r2, r3 — only the middle is filled.
`sgrna_vector` is empty because no PAM VCF is supplied; VaLiAnT accepts that.

The GTF is 51 lines of long attribute strings, but only one record matters — the CDS
covering r2, whose **phase** column fixes the reading frame:

```
chr17	HAVANA	CDS	43115726	43115779	.	-	1	gene_id ENSG00000012048.23; transcript_id ENST00000357654.9; ...
```

`brca1_clinvar.vcf`, `brca1_gnomad.vcf` and the custom-variant manifest sit in the
same directory but belong to `brca1_nuc`, which passes `--vcf`. Nothing here reads them.

---

## 3. Mutators and the oligo count

| Mutator | What it does | Per codon |
|---|---|---|
| `aa` | substitute the codon to every other amino acid (20 non-stop − 1) | **19** |
| `inframe` | delete the complete reading-frame triplet | **1** |
| `stop` | substitute the codon to a stop codon | **1** |
| | | **21** |

The codon count comes from the GTF, not from the targeton file — `r2_start`/`r2_end`
are user-chosen and here are **not codon-aligned at either end**:

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
completes a codon begun in exon 1; 43115726-43115727 begin one finished in exon 3.
Neither tool substitutes an amino acid without the whole codon, so both are skipped.

```
n_codons = (r2_length − phase) // 3       = (54 − 1) // 3 = 17
oligos   = 21 × n_codons                  = 357
```

For a minus-strand targeton the sense-strand ORF extent is
`(ref_end − r2_end + phase, + 3·n_codons)` = **(100, 151)**.

---

## 4. VaLiAnT

```bash
valiant sge parameter_input_files/targeton.tsv \
    ../ref/chr17.fa  brca1_modified_output  'homo sapiens' 'GRCh38' \
    --revcomp-minus-strand \
    --adaptor-5 AATGATACGGCGACCACCGA --adaptor-3 TCGTATGCCGTCTTCTGCTTG \
    --gff reference_input_files/ENST00000357654.9.gtf
```

Emits `_meta.csv` (32 columns), `_unique.csv` (357), two VCFs, `ref_sequences.csv`
and `config.json`.

---

## 5. PoolParty

```python
ADAPTOR_5 = "AATGATACGGCGACCACCGA"
ADAPTOR_3 = "TCGTATGCCGTCTTCTGCTTG"

orf = pp.from_fasta("chr17.fa", ("chr17", 43115633, 43115878, "-")) \
        .annotate_orf(region_name="cds", extent=(100, 151), frame=1)

aa      = orf.mutagenize_orf(region="cds", num_mutations=1,
                             mutation_type="missense_only_first", mode="sequential")
inframe = orf.deletion_scan_orf(deletion_codons=1, deletion_marker=None,
                                region="cds", mode="sequential")
stop    = orf.insertion_scan_orf("TGA", region="cds", replace=True, mode="sequential")

library = pp.join([ADAPTOR_5, pp.stack([aa, inframe, stop]), ADAPTOR_3])
```

Five statements. `orf` fans out to three operations, `stack` merges them, `join`
adds the adaptors — the shape of Figure 2. Method form throughout, matching the GB1
tutorial; `stack` and `join` are module-level because they take several pools.

Every argument is load-bearing. `mode="sequential"` enumerates rather than samples;
`deletion_marker=None` shortens rather than writing gap characters;
`replace=True` substitutes rather than inserts; `region="cds"` confines the scan to
the exon. `mutation_type="missense_only_first"` is the default but is written out,
because it is the line that corresponds to VaLiAnT's `aa`.

Two mappings worth naming:

- **`stop` is a fixed-codon replacement, not a mutagenesis call.** PoolParty has no
  residue-targeted mutation type, so a stop scan is expressed as
  `insertion_scan_orf("TGA", replace=True)`. Same output, different model.
- **`deletion_scan_orf` takes codon units.** `deletion_codons=1`, no offset
  arithmetic, and the scan cannot run past the ORF. Before this operation shipped, the
  same result needed `deletion_scan(deletion_length=3, region=..., positions=slice(0, None, 3))`.

---

## 6. Result

| | aa | inframe | stop | Total |
|---|---|---|---|---|
| PoolParty states | 323 | 17 | 17 | **357** |
| distinct | | | | **357** |
| VaLiAnT oligos | | | | **357** |
| **exact match** | | | | **357 / 357 (100%)** |

Zero sequences unique to either tool. Oligo lengths 286 nt, and 283 for the
`inframe` deletions — so `join` and the adaptors are verified, not just the
mutagenesis.

**One parameter had to be aligned.** Under PoolParty's default codon table, agreement
is 95%, and every mismatch is the arginine substitution — one per codon. PoolParty's
Kazusa-ordered table puts `AGA` first; VaLiAnT picks `CGG`. Reordering that single
entry through `pp.set_genetic_code()` gives 100%. The design spaces never differed —
only the DNA realisation of one amino acid. That is the substantive answer to R2 2b's
"investigate the unique elements".

---

## 7. What this shows, and what it does not

| Step | Who |
|---|---|
| Extract the targeton from `chr17.fa`, reverse-complemented | **PoolParty** (`from_fasta`, `strand='-'`) |
| Read the GTF's CDS phase; choose region and frame | **the user** |
| Generate the variants and add adaptors | **PoolParty** |

**PoolParty is not annotation-aware.** It does not read GTF/GFF files and has no
transcript model, so the target region and reading frame are supplied by the user —
here derived from the GTF's CDS phase. That is `transcript_annotation_aware` = no in
Table 2, made concrete, and it is the only step neither tool performed for us.

**One targeton, not four.** `from_fasta` accepts a list of coordinates, but
`annotate_orf` takes a single `extent`, and the four `brca1_pep` exons have different
ORF windows — (100,151), (67,142), (83,170), (71,209) — because phase and r2 length
differ. VaLiAnT handles four targetons in one invocation; PoolParty needs one DAG each.
Real asymmetry, and it is about annotation handling rather than library design.

**Split codons are skipped, not resolved.** VaLiAnT's documented
codons-across-exons capability is not exercised here.

**Not compared:** the metadata table, the VCFs, and the min/max length filter.

---

## 8. For the figure

- **their input** — `targeton.tsv` complete, two lines
- **their annotation** — the one CDS line above, showing `phase`
- **their command** — the `valiant sge` invocation from §4
- **our code** — the five statements from §5
- **caption** — the 357/357 byte-identity, the scope sentence from the top, and that
  the reading frame is user-supplied because PoolParty is not annotation-aware

The byte-identity does not need to appear in the figure itself.
