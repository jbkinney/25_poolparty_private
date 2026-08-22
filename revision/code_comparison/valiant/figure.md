# Supplementary figure — VaLiAnT / PoolParty code comparison

Design record. The rendered figure is not in this repository yet; this file is the
specification it should be drawn from, and the source of the caption text.

**Orientation:** landscape. **Width available:** 552 pt at `\textwidth` in
`sn-jnl`. **Code type size:** 7 pt monospace (≈ 4.2 pt/char), so a half-width
column holds ≈ 66 characters and full width ≈ 131.

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

> **Supplementary Figure N. One library specified in VaLiAnT and in PoolParty.**
> **(A)** The library: BRCA1 exon 2 on the minus strand of chr17, with 25 nt of
> intronic flank either side, every codon substituted to each other amino acid
> (`aa`, 19 per codon) and to a stop codon (`stop`), and each codon deleted in frame
> (`inframe`) — 357 oligos, each carrying P5 and P7 adaptors. The exon is not
> codon-aligned: one base at its 3′ end and two at its 5′ end belong to codons split
> across the neighbouring exons and are not mutated by either tool. **(B)** VaLiAnT
> 4.0.0 takes three input files and a command-line invocation; `targeton.tsv` is
> shown transposed for legibility, and only the CDS record of the transcript
> annotation is shown. **(C)** The same library in PoolParty. The arrow marks the
> reading frame, which VaLiAnT derives from the transcript annotation and PoolParty
> requires as an argument. **The two libraries are byte-identical: all 357 oligos,
> including adaptors.** Agreement requires one shared parameter, the codon usage
> table — VaLiAnT selects CGG for arginine where PoolParty's default selects AGA;
> without aligning it, 340 of 357 oligos match. This library is not a VaLiAnT shipped
> example: it is their BRCA1 exon 2 targeton with `stop` added and PAM protection
> removed. The same VaLiAnT installation reproduces their shipped `brca1_pep`
> library exactly, 2,339 of 2,339 oligos, and passes its own `md5` validation.

---

## Cut from an earlier draft, and why

A trimmed variant moved the split-codon annotations, the adaptor lengths, the file
size and the elided flags out of the figure and into the caption, on the grounds that
a figure shows *what* and a caption explains *why*. It was rejected: the arithmetic
`54 = 1 + 51 + 2` is what makes the 357 checkable by eye, and a reader who has to
reconstruct it from prose loses that. The fuller panel A is the point of the figure,
not decoration.
