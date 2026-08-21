# `from_vcf` — design

Status of every item is marked. **DECIDED** items were settled in discussion;
**OPEN** items are not. Nothing is blocking implementation.

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
`chrom` · `pos` · `ref` · `alt` · `allele` · `variant_id` · `filter` ·
`window_start` · `window_stop`, plus `info:`-prefixed keys on request, e.g.
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

### `FILTER` column: every record is kept — **DECIDED 2026-08-21**

`from_vcf` does not consult the VCF `FILTER` column. A record marked `LowQual`,
or anything else, produces a sequence like any other.

PASS-only would have been the other reasonable default, but silently discarding
rows from a file the user handed us is the wrong failure: the library comes back
smaller than the VCF and nothing says why.

This is consistent with the skip list rather than in tension with it. Symbolic
alts, MNPs and off-contig windows are skipped because **no sequence can be built**
for them — there are no bases to substitute, or no window to cut. `FILTER` is a
quality judgement made by whoever ran the variant caller, not a structural
impossibility. The rule is: skip what cannot be represented, keep what we merely
might not like.

**Consequence:** `filter` becomes a design card key, so a user who does want
PASS-only can select on it downstream and see exactly what they dropped.

### Out of scope for v1 — **DECIDED**
- `pool`/`region` splicing — no clear meaning for 10^4 variant windows.
- N-padding to gene boundaries — requires annotation.
- Per-variant strand — requires annotation. `strand` applies to the whole call,
  which covers the real case of a single gene on the minus strand.
- VCF **output**. See `README.md`.

### State generation: eager — **DECIDED 2026-08-21**

`from_vcf` parses the VCF, extracts every window, and hands the list to the
existing `from_seqs` machinery, exactly as `from_fasta` batch mode does.

A lazy variant was argued for and rejected. The argument for it was that PoolParty
declares libraries without materialising them, so a source op that reads everything
into memory violates the architecture. That conflates two different things:
**laziness exists for state spaces that are *generated*** — `mutagenize_orf(
num_mutations=3)` is 576,156 states from one short specification, and materialising
that is what the DAG exists to avoid. A VCF is a finite list already sitting on
disk. Reading it is not a combinatorial explosion; it is input.

Measured, on GRCh38 with ClinVar coordinates over `/mnt/c`:

| | |
|---|---|
| `Fasta()` open, 194 contigs | 0.03 s |
| One random 10,001 bp window via `.fai` | 0.18–0.25 ms |
| 20,000 variants (saturation scale), eager | ~0.4 GB, ~4 s |
| 200,000 variants | ~4 GB |
| All of ClinVar (3,375,801 records) | ~68 GB — will not run |

Lazy would not have reduced the I/O, only deferred it: each window is read once
either way. It would also have cost a `Fasta` handle whose lifetime must survive
operation copies, for a benefit that only appears at a scale nobody designs a
library at.

**Consequence to document, not to guard against:** memory scales with record count
times window width. A user pointing this at an unfiltered ClinVar VCF will exhaust
memory. That is a loud failure, not a silent wrong answer, and the docstring should
say plainly that the whole VCF is read into memory and large files should be
pre-filtered.

### VCF parsing: stdlib — **DECIDED 2026-08-21**

No new dependency. `gzip`, `io` and `re` ship with Python, so nothing is added to
`pyproject.toml` — which matters because `main.tex:148` claims minimal
dependencies, and `pysam`/`cyvcf2` both require htslib.

This works because **`from_vcf` never needs random access.** It iterates every
record. Verified against a real bgzipped ClinVar file (149 MB, BGZF confirmed by
its `BC` extra field): stdlib `gzip` read 200,000 records in 0.6 s, and all
3,375,801 in 9 s.

Not supported, and documented as such: BCF (binary) input, and tabix random access.
Neither is needed.

Risks accepted: INFO percent-encoding must be decoded by hand (~5 lines, and only
matters for `info:` card keys). Everything else in the VCF spec that is hard —
structural variants, breakends, multi-sample genotypes — we skip or discard anyway.

## Open — low stakes, awaiting confirmation

| Item | Recommendation |
|---|---|
| **Module placement** | `base_ops/`. `from_fasta` sits in `fixed_ops/` because it has a single-coordinate mode that yields one fixed sequence; `from_vcf` has no such mode — a VCF is inherently a list and every record is a library member — so it is purely state-multiplying. |
| **`mode`** | Hard-code `"sequential"`, not exposed. This is not a special case: `from_fasta`'s batch path already hard-codes `mode="sequential"` when it calls `from_seqs`. `from_vcf` is batch-shaped only, so the precedent applies directly. |
| **Test fixtures** | Synthetic and tiny — a few-kb FASTA and a ten-line VCF, `.fai` generated rather than committed. **No reference genomes or real VCFs in this repository**: hg38 is ~3 GB, a ClinVar VCF ~150 MB, and `.git` is already 95 MB. |

## Verification available at implementation time

Both files needed to test against real data already exist locally and must **not**
be copied into this repository:

- `KinneyLab/xin_collab/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa` (+ `.fai`)
- `KinneyLab/VEP_DNA/clinvar_hg38/clinvar.vcf.gz`

Note that the ClinVar file uses `1`-style contig names while UCSC-derived FASTAs
use `chr1`, so this pair is a live test of contig normalisation.

The strongest correctness check is that SpliceAI's own window construction is
reproducible: with `flank_left = flank_right = 5000`, `from_vcf` should emit
byte-identical sequences to `spliceai/utils.py` for the same records. That is also
what R2 2b asks for -- a comparison of pool outputs against another tool.
