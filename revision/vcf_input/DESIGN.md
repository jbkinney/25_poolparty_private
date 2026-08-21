# `from_vcf` — design

Status of every item is marked. **DECIDED** items were settled in discussion;
Every design question is now settled; there are no OPEN items. What remains is
implementation.

## Proposed signature

```python
def from_vcf(
    vcf_path,
    fasta_path,
    *,
    alleles="both",        # "ref" | "alt" | "both"
    flank_left,               # required, no default
    flank_right,              # required, no default
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
library size.

Worked example. One VCF record with two ALT alleles:

```
chr1  1000  rs123  A  G,T
```

With `alleles="both"` this yields three states:

| # | `allele` | `ref` | `alt` | sequence |
|---|---|---|---|---|
| 1 | `ref` | `A` | *(empty)* | window with `A` at the centre |
| 2 | `alt` | `A` | `G` | window with `G` |
| 3 | `alt` | `A` | `T` | window with `T` |

**`alt` is empty on a REF row.** There are two ALT alleles but only one REF
sequence, so no single ALT belongs to it. Filling it with `G,T` would imply a
sequence containing both, and repeating `A` would imply a substitution that did not
happen.

**Pairing is on `(chrom, pos, ref)`, not on `variant_id`.** An earlier version of
this document said `variant_id` — that was wrong. The VCF `ID` column is optional
and is `.` in most variant-caller output, so every record would share one value and
pairing would collapse. `(chrom, pos, ref)` identifies the record, and all three are
already card keys. To score ALT − REF for state 2, take the row with the same
`(chrom, pos, ref)` and `allele == "ref"`.

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

### Flanks are required, not defaulted — **DECIDED 2026-08-21**

`flank_left` and `flank_right` have no default. The caller must state them.

A default of 5000 would silently import SpliceAI's window into a general-purpose
tool. It is right for SpliceAI and absurd for a 200 bp MPRA tile, and window size
is the single most consequential parameter here — the amount of sequence context a
model sees. Making it explicit costs one line of typing and turns an inherited
assumption into a deliberate choice.

They remain separate rather than one symmetric `flank=`, because an even-width
window cannot be centred and SpliceAI's `wid//2` silently favours the left. Stating
both makes that choice visible.

### `alleles` — default `"both"`, unrecognised values raise — **DECIDED 2026-08-21**

`alleles="both"` emits REF and ALT members; `"ref"` and `"alt"` emit one side.
Anything else raises an error that names all three valid options, rather than a
generic failure or a silent fallback.

### Indel anchor: centre on POS, as SpliceAI does — **DECIDED 2026-08-21**

VCF represents a deletion with an anchor base. Deleting the `G` from `AG` is
written:

```
chr1  1000  .  AG  A
```

so `POS` points at the **anchor `A`, not the deleted `G`**. The window is centred
on `POS`, matching `spliceai/utils.py`, which slices
`[pos - wid//2 - 1 : pos + wid//2]` without adjusting for the anchor.

This must be documented explicitly. Anyone reasoning about a deletion naturally
assumes `POS` marks the deleted base, and for a long deletion the difference moves
the window by the length of the deletion.

### Malformed input: one boundary check — **DECIDED 2026-08-21**

A data line that does not split into at least 8 tab-separated fields raises,
reporting the line number, the truncated line, and what was expected.

That is the whole validation. `from_vcf` is not a VCF linter: it does not check
INFO syntax, base alphabets, POS ordering, or header consistency. One check at the
boundary catches the malformations that actually occur — space-delimited files,
truncated lines, and files that are not VCFs at all.

**No whitespace fallback.** A space-delimited file exists in this project's own
data (`VEP_DNA/annotation/test.vcf`), and splitting on arbitrary whitespace to
accommodate it would guess at intent and could mask real corruption. It is not a
VCF; say so clearly.

### `num_states` — **DECIDED 2026-08-21**

`num_states` is the number of sequences the pool yields, with no separate concept
behind it. For a VCF with **N records** carrying **M ALT alleles in total**:

| `alleles` | states |
|---|---|
| `"alt"` | M |
| `"ref"` | N |
| `"both"` | N + M |

Most VCFs have no multi-allelic rows, so M = N and `"both"` gives 2N.

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

### Module placement: `base_ops/` — **DECIDED 2026-08-21**

Not a judgement call — it is what the directories mean. `fixed_operation` takes
`seq_from_seqs_fn: Callable[[list[str]], str]`: **one output sequence per input
state**. So `fixed_ops/` holds deterministic transforms (`rc`, `upper`, `lower`,
`slice_seq`, `add_prefix`, `join`, `clear_gaps`) and `base_ops/` holds operations
that create or multiply states (`mutagenize`, `recombine`, `shuffle_seq`,
`get_kmers`, `get_barcodes`, `from_seqs`, `from_iupac`, `from_motif`).

`from_fasta` sits in `fixed_ops/` because its single-coordinate mode yields one
fixed sequence; its batch mode delegates out to `from_seqs` in `base_ops/`.
`from_vcf` has no single-record mode and turns a file into N states, so it belongs
in `base_ops/`.

### `mode="sequential"`, hard-coded — **DECIDED 2026-08-21**

Not exposed as an argument. `from_fasta`'s batch path already hard-codes
`mode="sequential"` when calling `from_seqs`, so this follows precedent rather than
inventing a special case. A VCF is an ordered list and the natural traversal is
file order.

### Test fixtures — **DECIDED 2026-08-21**

Synthetic and tiny: a few-kb FASTA and a ten-line VCF, with the `.fai` generated at
test time rather than committed.

**No reference genomes or real VCFs in this repository.** hg38 is ~3 GB, a ClinVar
VCF ~150 MB, and `.git` is already 95 MB.

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
