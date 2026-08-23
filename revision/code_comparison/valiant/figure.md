# Supplementary figure — VaLiAnT / PoolParty code comparison

Design record. The rendered figure is not in this repository yet; this file is the
specification it should be drawn from, and the source of the caption text.

**Shape:** upright, full width, wide aspect — see `../README.md` for the width
budget and the shared drawing spec. Full width holds ≈ 130 characters, so a
half-width column holds ≈ 64.

An earlier version of this record specified a landscape page at 552 pt. The
manuscript has no sideways float, and the character budget is the same either way,
so the layout below is unaffected — only the justification changed.

The PoolParty code is reformatted to a narrow column — one or two arguments per
line — because its natural formatting runs to 86 characters and will not fit
beside panel B. That is a deliberate cost, paid to keep the two specifications in
view at once, which is the figure's whole purpose.

---

## Layout

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│ (A)  One library: BRCA1 exon 2, three mutation types, P5/P7 adaptors            │
│                                                                                 │
│      chr17:43115634-43115878   minus strand, 245 nt                             │
│      ┌──────────┬────────────────────────────────────┬──────────┐               │
│      │   r1     │        r2 = exon 2, 54 nt          │   r3     │               │
│      │  25 nt   │  1 │◄──── 17 codons, 51 nt ────►│ 2 │  25 nt   │              │
│      └──────────┴────────────────────────────────────┴──────────┘               │
│                   ▲                                  ▲                          │
│              phase = 1                        split codon                       │
│           (split codon, exon 1)              (continues, exon 3)                │
│                                                                                 │
│      17 codons × (19 aa + 1 inframe + 1 stop)  =  357 oligos                    │
│      [P5 20 nt]────────── 245 nt ──────────[P7 21 nt]   =  286 nt each          │
└─────────────────────────────────────────────────────────────────────────────────┘
┌────────────────────────────────────┬────────────────────────────────────────────┐
│ (B)  VaLiAnT 4.0.0                 │ (C)  PoolParty                             │
│                                    │                                            │
│  targeton.tsv                      │  ADAPTOR_5 = "AATGATACGGCGACCACCGA"        │
│   ref_chr      chr17               │  ADAPTOR_3 = "TCGTATGCCGTCTTCTGCTTG"       │
│   ref_strand   -                   │                                            │
│   ref_start    43115634            │  orf = pp.from_fasta(                      │
│   ref_end      43115878            │      "chr17.fa",                           │
│   r2_start     43115726            │      ("chr17", 43115633, 43115878, "-")    │
│   r2_end       43115779            │  ).annotate_orf(                           │
│   ext_vector   "25,25"             │      region_name="cds",                    │
│   action_vec   "(),(aa,inframe,    │      extent=(100, 151), frame=1)  ◄──┐     │
│                  stop),()"         │                                      │     │
│                                    │  aa = orf.mutagenize_orf(            │     │
│  ENST00000357654.9.gtf             │      region="cds", num_mutations=1,  │     │
│   CDS 43115726 43115779 . - 1      │      mutation_type=                  │     │
│                          ▲         │        "missense_only_first",        │     │
│                        phase ──────┼──── transcribed by hand ─────────────┘     │
│                                    │      mode="sequential")                    │
│  chr17.fa   (81 MB)                │                                            │
│                                    │  inframe = orf.deletion_scan_orf(          │
│  $ valiant sge targeton.tsv \      │      deletion_codons=1,                    │
│      ../ref/chr17.fa out \         │      deletion_marker=None,                 │
│      'homo sapiens' 'GRCh38' \     │      region="cds", mode="sequential")      │
│      --revcomp-minus-strand \      │                                            │
│      --adaptor-5 ... \             │  stop = orf.insertion_scan_orf(            │
│      --adaptor-3 ... \             │      "TGA", region="cds",                  │
│      --gff ENST...gtf              │      replace=True, mode="sequential")      │
│                                    │                                            │
│  3 input files + CLI               │  library = pp.join([ADAPTOR_5,             │
│                                    │      pp.stack([aa, inframe, stop]),        │
│                                    │      ADAPTOR_3])                           │
└────────────────────────────────────┴────────────────────────────────────────────┘
```

---

## Panel notes

**(A)** Establishes what the action vector's three slots mean. The `1` and `2`
inside the r2 bar are the orphaned bases of codons split across exon boundaries —
43115779 completes a codon begun in exon 1, and 43115726-43115727 begin one
finished in exon 3. Neither tool substitutes an amino acid without the whole codon,
so both are skipped, which is why 54 nt yields 17 codons and not 18.

**(B)** `targeton.tsv` is shown **transposed**; the file itself is two tab-separated
lines, 94 characters wide, which will not fit the column. Say so in the caption.
The GTF is 51 lines of long attribute strings, so only the CDS record covering r2
appears — it is the one that carries `phase`. Flags on the `valiant sge` line are
elided; the count of input files is what matters.

**(C)** Method form throughout, matching `docs/tutorials/dms_gb1.rst` and Figure 2.
`stack` and `join` stay module-level because they take several pools. `.named()` and
`prefix=` are omitted: they feed name and design-card output, and this comparison is
on sequences.

Every argument shown is load-bearing. `mode="sequential"` enumerates rather than
samples; `deletion_marker=None` shortens rather than writing gap characters;
`replace=True` substitutes rather than inserts; `region="cds"` confines each scan to
the exon. `mutation_type="missense_only_first"` **is** the default but is written
out, because it is the line that corresponds to VaLiAnT's `aa`.

**The arrow** links the GTF's `phase` to `extent=(100, 151)`. It is the only element
of the figure that carries an argument rather than a fact, and it is there
deliberately: PoolParty does not read annotation files, so that value is supplied by
the user. Omitting the arrow would let the figure imply an equivalence that does not
hold.

---

## Caption

Written for a reader, not a referee.

**The line is tool coinage, not technical vocabulary.** Field-conventional terms are
fine unglossed — *codon*, *exon*, *intron*, *transcript*, *reading frame*, *in
frame*, *phase*, *CDS*, *adaptor*, *minus strand*, *stop codon*, *codon usage
table*. A reader of this journal knows them, and paraphrasing them makes the caption
longer and vaguer, not clearer.

What must not appear unexplained is vocabulary invented by the tool being described:

| VaLiAnT coinage | Handling |
|---|---|
| **targeton** | avoided. Its own term for the designed region; carries no meaning outside the tool |
| **r1 / r2 / r3** | unavoidable — they label panel A and appear as columns in panel B — so tied to plain words once: *"exon 2 in the middle, flanked by 25 bp of intron on either side (labelled r2, r1 and r3 …)"* |
| **action vector**, **ext_vector** | avoided; described by what they do |
| mutator names **aa**, **stop**, **inframe** | shown in panel B, described in the caption by their effect rather than their names |

Avoiding `targeton` is preferred but not absolute — one use with a gloss would be
acceptable if the sentence needed it. It does not.

> **Supplementary Figure N. The same variant library, specified in VaLiAnT and in
> PoolParty.**
>
> **(A)** The library. Within exon 2 of BRCA1, each codon in turn is replaced by
> every other amino acid, replaced by a stop codon, or deleted in frame — 357
> variants, each carrying sequencing adaptors. The exon contains 17 complete codons;
> its terminal codons are split across the neighbouring exons and are left unchanged.
> VaLiAnT's input table labels the exon r2, and the flanking intron r1 and r3.
>
> **(B)** VaLiAnT runs from the command line, with the design written into three
> files: a table of coordinates and the changes to make in them, the gene's
> transcript annotation, and the reference chromosome (table transposed, annotation
> abridged).
>
> **(C)** The same design in PoolParty.
>
> The two libraries are identical base for base, adaptors included, once both tools
> use the same codon usage table. The arrow marks the one value PoolParty cannot
> derive: VaLiAnT reads the annotation's phase field to locate the first complete
> codon, whereas PoolParty is given that position. This library is not one of
> VaLiAnT's worked examples; we defined it using their coordinates for BRCA1 exon 2.

**193 words.** Trimmed from 267 by cutting what the figure already shows (the 25 /
54 / 25 structure), compressing figure mechanics to a parenthetical, dropping the
statement count from panel C, and folding the codon-table caveat into the result
sentence.

**If R2 2b is answered in the paper rather than the response letter**, restore the
specifics: *VaLiAnT chooses CGG for arginine where PoolParty's default chooses AGA.*
Naming them makes clear the disagreement is a single parameter for a single residue,
which is more defensible than an unqualified "same codon usage table".

Deliberately not in the caption: that the same installation reproduces VaLiAnT's
shipped `brca1_pep` library exactly (2,339 of 2,339) and passes its `md5`
validation. That licenses the ground truth, but it is provenance for a referee, not
information for a reader — it belongs in the supplementary text or the response
letter.

## Cut from an earlier draft, and why

A trimmed variant moved the split-codon annotations, the adaptor lengths, the file
size and the elided flags out of the figure and into the caption, on the grounds that
a figure shows *what* and a caption explains *why*. It was rejected: the arithmetic
`54 = 1 + 51 + 2` is what makes the 357 checkable by eye, and a reader who has to
reconstruct it from prose loses that. The fuller panel A is the point of the figure,
not decoration.
