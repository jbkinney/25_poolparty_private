# `from_vcf` — design

Status of every item is marked. **DECIDED** items were settled in discussion;
**OPEN** items are not, and two of them block implementation.

## Proposed signature

```python
def from_vcf(
    vcf_path,
    fasta_path,
    *,
    alleles="both",        # "ref" | "alt" | "both"
    flank_left=5000,
    flank_right=5000,
    strand="+",
    # inherited from the existing source ops, unchanged
    prefix=None, style=None, cards=None,
    mode=..., num_states=None, iter_order=None,
) -> DnaPool
```

Six new arguments. Everything else is the argument surface `from_fasta` and
`from_seqs` already expose.

## Why this is shaped like `from_fasta`

`from_fasta` already reads a reference with `pyfaidx`, takes `(chrom, start, stop,
strand)` coordinates, turns a batch into N states via `from_seqs`, and supports
`prefix`, `style`, `cards`, `iter_order`. `from_vcf` differs only in where the
coordinates come from and in emitting a REF/ALT pair per record. Most of the design
is inherited rather than invented.

## Decisions taken

### One pool, not two — **DECIDED**
REF and ALT members live in the same pool, distinguished by the `allele` design
card key. `alleles=` selects which are emitted.

### Indel length: match SpliceAI, no flag — **DECIDED**
`len(sequence) = flank_left + flank_right + len(ALT)`, so indels produce sequences
of differing length within one pool.

An `equal_length=True` option was considered and rejected. There is no coherent way
to trim: for an indel you can hold at most two of {fixed width, variant centred,
REF and ALT covering the same genomic span}. The only trim that keeps the variant
centred refills from the reference on one side, which makes the ALT window cover a
*different stretch of chromosome* than the REF window. Since the measurement is
ALT − REF, part of that difference would then be caused by unequal reference
context rather than by the variant. That is a scientific confound, not a cosmetic
one.

This is comfortable because PoolParty never emits arrays — output is sequences and
a DataFrame — there is no `pad`/`resize`/`one_hot` operation in the package, and
variable-length pools already occur (`deletion_scan(deletion_marker=None)`).
One-hot encoding happens in user code, which is where padding belongs and where the
choice is visible.

Documentation owes one sentence: *the variant begins at offset `flank_left` in the
emitted sequence; indels make sequences differ in length, so pad at encoding time
if your model needs uniform input.*

### Multi-allelic: 1 REF + N ALT — **DECIDED**
A row with two ALT alleles and `alleles="both"` yields three states, not two pairs.
The REF window is identical for both alleles and duplicating it would misstate the
library size. Pairing for ALT − REF is recoverable through the shared
`variant_id`.

### Coordinate conventions — **DECIDED**

| Card key | Convention | Round-trips to |
|---|---|---|
| `pos` | 1-based, exactly the VCF POS | the user's VCF |
| `window_start`, `window_stop` | 0-based, half-open | `from_fasta`, BED, pyranges |

Both round-trips work, which is the point. `pos` must equal the VCF POS or the card
lies about its own input. `from_fasta` documents *"Coordinates are 0-based [start,
stop)"*, so window coordinates feed straight back into it with no arithmetic.

Mixing conventions is deliberate and follows AlphaGenome (`PRIOR_ART.md`): the
field name carries the convention — `position` is 1-based, `start`/`end` are
0-based half-open.

The residual trap is that `pos - window_start` is off by one as a within-sequence
offset; the correct value is `(pos - 1) - window_start`. No card key is added for
it, because by construction the variant always begins at offset `flank_left`.

### Design card keys — **DECIDED**
`chrom` · `pos` · `ref` · `alt` · `allele` · `variant_id` · `window_start` ·
`window_stop`, plus `info:`-prefixed keys on request, e.g.
`cards=['chrom', 'pos', 'info:AF', 'info:CLNSIG']`.

Carrying the VCF `ID` column and selected `INFO` fields through to design cards is
what puts ClinVar and gnomAD identifiers on every sequence — part of what R3 wanted
from annotation, obtained on the input side.

Cut as redundant: `alt_index` (`variant_id` + `alt` already disambiguate) and
`strand` (a single argument, so identical on every row).

### Sequence naming — **DECIDED**
gnomAD/GTEx style, `.`-separated to match PoolParty's existing `wt.mut_03.bc_01`:

```
chr1_12345_A_G          alleles="alt"
chr1_12345_A_G.ref      alleles="both", REF member
chr1_12345_A_G.alt      alleles="both", ALT member
```

Underscores rather than AlphaGenome's display form `chr1:12345:A>G`, because
PoolParty exports FASTA and a `>` inside a header is hostile to parsers. Also
avoids resembling `from_fasta`'s `{chrom}:{start}-{stop}({strand})` names, which
use **0-based** coordinates — two source ops emitting similar-looking names with
different conventions would be a trap.

The rsID goes in the `variant_id` card, not the name: VCF IDs are frequently `.`,
so names must be systematic.

### Fixed behaviour, not arguments — **DECIDED**
Each of these was considered as an argument and rejected as configuration nobody
would tune:

- **REF is verified against the FASTA; a mismatch raises.** One string comparison
  on data already in hand, and it catches the catastrophic silent failure — an
  hg19 VCF against an hg38 FASTA yields a complete, plausible, entirely wrong
  library. Offering `warn`/`skip` would invite someone to ship that library.
- **Contig names are normalised** (`chr1` ↔ `1`).
- **Skipped, with a reported count:** symbolic alts, MNPs, off-contig windows.
- **Zero valid records raises** rather than returning an empty pool.
- **Sample genotypes are ignored entirely**, as AlphaGenome does.

### Out of scope for v1 — **DECIDED**
- `pool`/`region` splicing — no clear meaning for 10^4 variant windows.
- N-padding to gene boundaries — requires annotation.
- Per-variant strand — requires annotation. `strand` applies to the whole call,
  which covers the real case of a single gene on the minus strand.
- VCF **output**. See `README.md`.

## Open — blocking implementation

### Lazy or eager state generation — **OPEN**
`from_fasta` batch mode materialises every sequence before handing the list to
`from_seqs`. Copied here, a ClinVar-scale VCF (~3M records) at 10 kb windows is
~30 GB of strings.

PoolParty's operations are already lazy per row and `pyfaidx` does cheap random
access through the `.fai` index, so a lazy `from_vcf` could store only
`(chrom, pos, ref, alt)` per record — tens of bytes — and slice the reference when
a row is requested. That is ~100 MB of metadata for 3M records instead of 30 GB.

| | Eager | Lazy |
|---|---|---|
| Work | Delegate to `from_seqs` | Custom Operation |
| Ceiling | ~10^4–10^5 records | whole-VCF |
| Laziness claim | broken at the source | preserved |

MAVE-scale libraries are thousands of variants, so eager would work for R3's case.
But "libraries are declared, sequences generated only on request" is the paper's
central architectural claim, and a source operation that violates it is exactly
what a referee would notice. **Recommendation: lazy.**

### VCF parsing — dependency or stdlib — **OPEN**
Current dependencies are `numpy`, `pandas`, `beartype`, `statetracker`, `pyfaidx`,
`typing_extensions`. No `pysam`, no `cyvcf2`. `main.tex:148` claims minimal
dependencies.

`pysam` and `cyvcf2` both require htslib and are a common source of install
failure. The thing that makes them unnecessary: **`from_vcf` never needs random
access.** It iterates every record, and bgzipped VCFs are valid gzip streams, so
stdlib `gzip` reads `.vcf.gz` sequentially without tabix. A ~50-line parser
covering `CHROM POS ID REF ALT FILTER INFO` adds nothing to the dependency list.

**Recommendation: stdlib parser**, documenting that indexed random access is
unsupported because it is not needed.

## Open — not blocking

| Item | Note |
|---|---|
| **`FILTER` default** | PASS-only is what most pipelines assume, but silently dropping records from a user's VCF should be loud. A scientific default, not an engineering one. No recommendation. |
| **`mode` default** | Every other source op defaults to `"random"`. A VCF is an ordered list, which argues for `"sequential"` — but that adds a special case to an API whose value is not having special cases. No confident recommendation. |
| **Module placement** | `base_ops/` rather than `fixed_ops/`, since it is inherently state-multiplying. `from_fasta` sits in `fixed_ops/` but delegates to `from_seqs` in batch mode. Low stakes. |
| **Test fixtures** | Synthetic, tiny — a few-kb FASTA and a ten-line VCF. Generate the `.fai` rather than committing it. **No reference genomes or real VCFs in this repository**: hg38 is ~3 GB, a ClinVar VCF ~100 MB, and `.git` is already 95 MB. |
