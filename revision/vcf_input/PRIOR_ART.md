# How existing models turn a VCF into model input

Read from source on 2026-08-20, fetched from GitHub. Not from memory, and not from
the local copies under `KinneyLab/VEP_DNA/` — those were checked afterwards and
agree.

## Sources

| Tool | Repository | File read |
|---|---|---|
| SpliceAI | `Illumina/SpliceAI` @ `master` | `spliceai/utils.py` |
| AlphaGenome | `google-deepmind/alphagenome` @ `main` | `src/alphagenome/data/genome.py`, `docs/source/variant_scoring.md` |
| AlphaMissense | `google-deepmind/alphamissense` @ `main` | `alphamissense/data/pipeline_missense.py` |
| EVE | `OATML-Markslab/EVE` @ `master` | `utils/data_utils.py` |

## The headline: R3's model list is two different things

**DNA-sequence models** take a genomic window and one-hot encode it. A VCF plus a
reference genome is enough.

**Protein/MSA models** never see a DNA window. AlphaMissense's `ProteinVariant`
carries a protein `sequence` read from FASTA, a 1-based `position`, a
`reference_aa` and an `alternate_aa`, and builds
`alternate_sequence = sequence[:position-1] + alternate_aa + sequence[position:]`.
EVE's `MSA_processing(MSA_location=...)` consumes a multiple sequence alignment; its
`create_all_singles` builds mutants as `letter + str(pos) + mut`, so it has
protein-level positions but no genomic coordinate.

Getting from a VCF to a protein variant requires transcript models, reading frames
and split codons — i.e. VEP. That is the boundary of what a source operation can
reach.

## SpliceAI — the reference implementation

`get_delta_scores` in `spliceai/utils.py`:

```python
seq   = ann.ref_fasta[chrom][record.pos - wid//2 - 1 : record.pos + wid//2].seq
x_ref = 'N'*pad_size[0] + seq[pad_size[0]:wid-pad_size[1]] + 'N'*pad_size[1]
x_alt = x_ref[:wid//2] + str(record.alts[j]) + x_ref[wid//2 + ref_len:]
x_ref, x_alt = one_hot_encode(x_ref)[None, :], one_hot_encode(x_alt)[None, :]
```

Points that matter for our design:

1. **One row does not give one sequence.** The loop is
   `for j in range(len(record.alts)): for i in range(len(idxs)):` — once per ALT
   allele, and again per *overlapping gene*. Both the N-padding (`pad_size`, from
   the gene's boundaries) and the reverse-complement depend on which gene, so the
   same variant yields different sequences for different genes.
2. **It is always a REF/ALT pair**, never a single sequence.
3. **ALT and REF differ in length for indels.**
   `len(x_alt) = wid + len(ALT) - len(REF)`. They are fed to the model separately
   and the *predictions* re-aligned afterwards.
4. **The VCF REF is verified against the reference genome**, and the record skipped
   on mismatch:
   `if seq[wid//2 : wid//2+len(record.ref)].upper() != record.ref: skip`.
   This is the check that catches an hg19 VCF against an hg38 FASTA.
5. **Records skipped:** symbolic alts (`<`, `>`), alts containing `.`/`-`/`*`, REF
   longer than **`2*dist_var`** — the *scored span*, 100 bp by default, not the
   10,101 bp window — and windows running off a contig end (`len(seq) != wid`).
   **MNPs and delins are not skipped**: `utils.py:138-140` appends a null-score
   record and continues, because SpliceAI's output re-alignment cannot handle them.
   That is a constraint on predictions, not on sequence construction.
6. **Contig names are normalised** (`normalise_chrom`) so `chr1` and `1` interoperate.
7. It reads the reference with **`pyfaidx`** — the same library PoolParty's
   `from_fasta` already depends on.

## AlphaGenome — the modern abstraction

`Variant` is a plain dataclass: `chromosome`, `position`, `reference_bases`,
`alternate_bases`, `name`, `info`. Its docstring is explicit that this is *not* a
VCF row:

> "Differs from the Variant definition in a VCF file, which allows for multiple
> alternative bases and contains sample information. This `Variant` class does not
> include sample information or variant call quality information."

So multi-allelic rows are normalised apart and sample genotypes are discarded.

AlphaGenome's `docs/source/variant_scoring.md` (in its repository): *"the variant is treated as a pair of sequences:
reference (REF) and alternate (ALT)"*, and for indels *"the ALT allele's prediction
profile is aligned to the REF allele's coordinate space"* — the same
re-align-the-output strategy SpliceAI uses.

### The coordinate convention, which we adopt

AlphaGenome uses **both** conventions deliberately, and lets the *name of the
field* carry which one applies:

```python
Variant.position    # "The 1-based position of the variant"      <- VCF convention
Variant.start       # "the 0-based start position" = position-1  <- BED convention
Interval.start/end  # "0-based", width = end - start             <- BED convention
```

This is the field convention, not an inconsistency: `position` means 1-based (VCF,
GFF, Ensembl); `start`/`end` mean 0-based half-open (BED, pyranges, Python slices).

## What we take from this

- One state per **(record, ALT allele)**, sample data discarded — AlphaGenome.
- REF/ALT as a **pair**, unequal length for indels — both.
- **Verify REF against the FASTA** and fail — SpliceAI.
- **Normalise contig names** — SpliceAI.
- **Skip** symbolic alts, non-DNA alts (`.`/`-`/`*`), oversized alleles and
  off-contig windows — SpliceAI. **Not** MNPs.
- **1-based `pos`, 0-based half-open window** with the names carrying the
  convention — AlphaGenome.

## What we deliberately do not take

**Gene-aware padding and strand.** SpliceAI's N-padding to gene boundaries and its
per-gene reverse-complement both require an annotation file. Without one there is
no gene, so neither is meaningful. This is why `strand` is a single argument
applying to the whole call rather than a per-variant property, and why N-padding is
absent.
