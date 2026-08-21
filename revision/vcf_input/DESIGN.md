# `from_vcf` — design

Reads a VCF plus a reference FASTA and returns a `DnaPool` of windows around each
variant. Answers part of Reviewer 3's comment; see `README.md` for what it does and
does not claim.

**Revised 2026-08-21** after three independent reviews (correctness, adversarial,
simplicity). Seven decisions were reversed and four factual claims corrected; the
superseded reasoning has been removed rather than annotated. Nothing is
implemented.

## Signature

```python
def from_vcf(
    vcf_path,
    fasta_path,
    flank_left,
    flank_right,
    *,
    alleles="both",        # Literal["ref","alt","both"]
    prefix=None, style=None, cards=None, iter_order=None,
) -> DnaPool
```

Five arguments, four of them positional. `mode`, `num_states` and `strand` are not
exposed; see the decisions below.

## Decisions

### Sequence length varies; no padding, no trimming

`len(sequence) = flank_left + len(allele) + flank_right`. A deletion produces a
shorter ALT sequence than its REF; an insertion a longer one. Nothing is padded or
trimmed to compensate.

Padding was considered. Gap markers (PoolParty's own convention — `deletion_scan`
emits `---` by default) fix deletions but cannot fix insertions: there is nowhere
to put the extra bases. Gapping the *reference* instead makes REF and ALT equal
**within a pair**, but lengths still differ **across** the pool because the allele
width differs per variant, so `seq_length` is `None` either way. Trimming the
flanks to compensate makes the ALT window cover a different stretch of chromosome
than the REF window, so part of any ALT − REF difference would come from unequal
reference context rather than from the variant.

**The cost, stated rather than engineered around.** A pool whose sequences differ
in length has `seq_length is None`, and these operations refuse to run on it —
measured, each against a variable-length pool and an otherwise identical
uniform-length control:

```
subseq_scan   deletion_scan   insertion_scan   deletion_multiscan
recombine     shuffle_scan    mutagenize_orf
```

`shuffle_seq`, `rc`, `upper`, `filter` and the export operations are unaffected.
`mutagenize` is affected only in `mode="sequential"`, which needs a known
`seq_length`; with `mutation_rate` it works. `insertion_multiscan` and
`annotate_orf` are unresolved — both failed on the uniform control too, for
unrelated reasons, so they were not isolated.

A pool containing only SNVs has uniform length and none of this applies — the
`variant_type` card exists to make that filter trivial.

*Correction:* an earlier draft justified this by claiming variable-length pools
already occur via `deletion_scan(deletion_marker=None)`. Measured false —
`deletion_length` is a scalar, so every state is shortened identically and
`seq_length` stays defined. `from_vcf` is the first operation in PoolParty to
produce a genuinely variable-length pool.

### One reference sequence per site, always

All records at the same `(chrom, pos, ref)` share one REF window, because it *is*
the same window. De-duplication is unconditional; there is no switch.

Without it, ClinVar yields **319,648 redundant REF sequences** — identical DNA,
emitted up to 182 times at one position. Multi-allelic VCF rows (`A → G,T`) are the
case an earlier draft designed for, and ClinVar contains **zero** of them: it splits
alleles across separate records at the same position, which is the case that
actually occurs.

State counts over ClinVar (3,375,801 records, 3,056,153 distinct sites):

| `alleles` | states |
|---|---|
| `"alt"` | 3,375,801 |
| `"ref"` | 3,056,153 |
| `"both"` | 6,431,954 |

`num_states` is exactly the number of sequences the pool yields.

### Sequence naming, and no allele suffix

```
1_939398_G          the reference at this site
1_939398_G_INS      a variant at this site
```

GTEx form, underscore-separated — `chr22_1024_A_C` minus the optional build suffix.
This is a published convention with a parser in AlphaGenome, whose `VariantFormat`
enum recognises it alongside gnomAD (`22-1024-A-C`) and its own default
(`chr22:1024:A>C`). The colon form is unusable here because PoolParty exports FASTA
and a `>` inside a header is hostile to parsers. This also avoids resembling `from_fasta`'s
`{chrom}:{start}-{stop}({strand})` names, which use 0-based coordinates.

**No `.ref`/`.alt` suffix.** Once each reference appears once per site, reference and
variant names cannot collide, so the suffix identifies nothing that position does not:
a reference name is always a strict prefix of its variants' names, and by field
position the second-from-right underscore field is digits for a reference and DNA
for a variant.

The rsID goes in the `variant_id` card, not the name — VCF IDs are frequently `.`,
so names must be systematic.

**Dotted contig names are fully supported.** RefSeq and GenBank references use
versioned accessions throughout — `GCF_000001405.39_GRCh38.p13_genomic.fna` names
chromosome 1 `NC_000001.11` — so a rule that skipped or mishandled dotted contigs
would discard entire files. The field-position test above works on them; only the
combination of a dotted contig, downstream operations appended, *and* no design
cards requested leaves a name ambiguous, because `.` is both the accession
separator and PoolParty's DAG token separator. The `allele` card is authoritative in
that case. This affects 28 of 3,375,801 ClinVar records and none of the
`alleles="alt"` workflow.

### Design cards

`chrom` · `pos` · `ref` · `alt` · `allele` · `variant_type` · `variant_id` ·
`filter` · `window_start` · `window_stop`, plus `info_`-prefixed keys on request.

Carrying the VCF `ID` column and selected `INFO` fields is what puts ClinVar and
gnomAD identifiers on every sequence — the part that actually answers R3, obtained
on the input side.

`alt` is `None` on a reference row. With one reference per site it pairs with
several variants, so no single ALT belongs to it. `None` rather than `""` because
only `None` is found by `isna()` and survives a CSV round-trip.

`variant_type` is `"snv"`, `"insertion"` or `"deletion"`, following AlphaGenome's
classification (`genome.py`): `snv` when both alleles are one base, otherwise
`insertion` or `deletion` by which is longer.

It earns a key because it makes the variable-length limitation navigable rather than
merely documented. The seven operations that refuse a variable-length pool all work
on a SNV-only pool, and `variant_type` turns that into one filter instead of the
user re-deriving it from `len(ref)` and `len(alt)`:

```python
df[df.variant_type == "snv"]     # uniform length; all seven operations available
```

AlphaGenome also exposes `is_frameshift` and `is_structural`. Neither is adopted:
frameshift is a property of the coding context rather than the variant alone, and
`structural` is an arbitrary 50 bp threshold.

Two mechanics to respect. With list-style `cards=[...]`, `generate_library` prefixes
every column with the operation name, giving `op[0]:from_vcf.allele`; bare column
names require dict-style `cards={'alt': 'alt', ...}`. And `info_AF` rather than
`info:AF`, because a colon makes the column unusable with `df.query()`.

This is ten static keys against a package maximum of five. The deviation is
deliberate: the identity fields are what make the pool traceable to its input, which
is the feature, and `variant_type` is what makes the variable-length limitation
workable.

### Coordinates: `pos` 1-based, windows 0-based half-open

| key | convention | round-trips to |
|---|---|---|
| `pos` | 1-based, exactly the VCF POS | the user's VCF |
| `window_start`, `window_stop` | 0-based, half-open | `from_fasta`, BED, pyranges |

Mixing them is deliberate and follows AlphaGenome, where the field name carries the
convention: `Variant.position` is 1-based, `Interval.start/end` are 0-based
half-open. `pos` must equal the VCF POS or the card lies about its own input, and
`from_fasta` documents "Coordinates are 0-based [start, stop)".

`pos - window_start` is off by one as a within-sequence offset; the correct value is
`(pos - 1) - window_start`. On a reference row the window describes the sequence; on
a variant row it describes the reference span the variant replaces, so
`window_stop - window_start != len(seq)` for indels.

### Indel anchor: centre on POS

VCF writes a deletion with an anchor base — deleting `G` from `AG` is
`chr1 1000 . AG A` — so POS marks the **anchor `A`, not the deleted `G`**. The
window is centred on POS, matching `spliceai/utils.py`. Documented because anyone
reasoning about a deletion assumes the opposite, and for a long deletion the
difference moves the window by its length.

### Always plus strand; no `strand` argument

Sequences are the reference plus strand. `pp.rc()` is exported and composes, so
`pp.rc(from_vcf(...))` covers the reverse-complement case without an argument.

*Correction, 2026-08-21:* an earlier draft justified dropping the argument by
saying it would leave the name reporting plus-strand bases while the sequence
carried their complements, with nothing recording which. That argument does not
hold — `pp.rc(from_vcf(...))` produces exactly the same mismatch, since `rc`
contributes nothing to the name. The real reason to omit `strand` is narrower: it
would duplicate an operation that already exists and composes.

### Guards — fixed behaviour, each reporting a count

Skipped records are counted and reported through `warnings.warn`, not dropped
silently: this operation's output is a library the user keeps, unlike SpliceAI's,
which is a score column.

| Guard | Reason | In ClinVar |
|---|---|---|
| ALT containing `.`, `-`, `*`, `<`, `>` | Not DNA. Pasted in, they produce a normal-length sequence with a literal `.` in it, and `.`/`*` are in PoolParty's `IGNORE_CHARS`, so `clear_gaps` would silently delete the variant position | 1,017 records |
| `len(REF)` or `len(ALT)` over a stated cap | A 9,983-base REF turns a requested 201 bp window into ~10 kb with no error. SpliceAI caps at `2*dist_var` (100 bp default) | max REF 9,983, max ALT 9,837 |
| Window running off a contig end | No window to cut | — |
| Contig absent from the FASTA | `pyfaidx` raises `KeyError`; a patch scaffold should not abort a 3M-record run | 28 records on scaffolds |

**REF verification, case-insensitively, log and skip.** The VCF's REF is compared
against the FASTA with `.upper()` on both sides, following
`spliceai/utils.py:117` — soft-masked references (UCSC `hg38.fa`, Ensembl `dna_sm`)
are lowercase over roughly half the genome, and a case-sensitive comparison would
reject every variant in a repeat region. A mismatch logs and skips the record rather
than raising; an hg19 VCF against an hg38 FASTA then fails visibly on essentially
every record, which is the diagnosis, while a handful of oddities do not abort the
run.

*Not a guard:* MNPs and delins are **kept**. An earlier draft skipped them, citing
SpliceAI. That was a misreading — `utils.py:138-140` *appends* a null-score record
and continues, because its output re-alignment cannot handle them. That is a
constraint on predictions, which PoolParty does not produce. `1 930222 GAACTC
TTCTTCTG` is perfectly representable. 14,738 ClinVar records affected.

### `FILTER` is not consulted

Every record produces a sequence regardless of `FILTER`. Silently discarding rows
from a supplied file is the wrong failure — the library comes back smaller than the
VCF with nothing saying why. `filter` is a card key so anyone wanting PASS-only can
select on it. Note it is `.` for all 3,375,801 ClinVar records, so it is inert on
that file and meaningful on caller output.

The rule this follows: **skip what cannot be represented, keep what we merely might
not like.**

### Flanks are required

No default. A default of 5000 would import SpliceAI's window into a general-purpose
tool — right for SpliceAI, absurd for a 200 bp MPRA tile — and window size is the
parameter that most changes the answer.

Separate rather than one symmetric `flank=` because asymmetric context is real for
splice and promoter work.

*Correction:* an earlier draft justified splitting them by claiming SpliceAI's
`wid//2` silently favours the left on an even window. False — `wid = 10000 + cov`
with `cov = 2*dist_var + 1` is always odd, and the slice is exactly symmetric.

### `alleles`

`Literal["ref","alt","both"]`, default `"both"`. `@beartype` already raises naming
all three valid options on a bad value, as it does for `get_kmers(case=)` and
`get_barcodes(padding_side=)`; no hand-rolled check.

### Eager, and its own `Operation` subclass

The VCF is parsed and every window extracted up front.

Laziness in PoolParty exists for state spaces that are *generated* —
`mutagenize_orf(num_mutations=3)` is 576,156 states from one short specification.
A VCF is a finite list already on disk. Measured: opening the FASTA costs 0.03 s, a
random 10,001 bp window 0.18–0.34 ms, and 20,000 variants take ~3.3 s and ~0.55 GiB
peak. Lazy would defer the same I/O, not reduce it. Memory scales with record count
times window width, so an unfiltered ClinVar VCF (~68 GB) exhausts memory — a loud
failure, documented in the docstring rather than guarded against.

**It cannot delegate to `from_seqs`.** `FromSeqsOp.design_card_keys` is fixed at
`["seq_name", "seq_index"]` and `Operation._validate_cards` is a hard whitelist:

```
ValueError: Invalid card key(s) ['allele','alt','chrom','pos','ref'] for from_seqs.
            Valid keys: ['seq','seq_index','seq_name','state']
```

So `from_vcf` needs a `FromVcfOp(Operation)` setting `design_card_keys` before
`super().__init__` — the `ScoreOp` pattern at `fixed_ops/score.py:95`, which is also
the precedent for the dynamic `info_` keys. `_get_copy_params` must carry the
extracted windows so `copy()` does not re-read the files.

*Correction:* an earlier draft said the design "hands the list to the existing
`from_seqs` machinery, exactly as `from_fasta` batch mode does", and that most of
the design was inherited. Both were wrong.

`prefix` is folded into the generated names and `prefix=None` passed downward;
`from_seqs` raises `ValueError: Cannot specify both seq_names and prefix`, and
`from_fasta` handles it the same way.

### VCF parsing: standard library only

`gzip`, `io` and `re`. No `pysam`, no `cyvcf2` — both need htslib, and
`main.tex` claims minimal dependencies.

This works because `from_vcf` never needs random access. Verified against the real
bgzipped ClinVar file (149 MB, BGZF confirmed by its `BC` extra field): stdlib
`gzip` read 200,000 records in 0.53 s and all 3,375,801 in 8.5 s.

Unsupported and documented: BCF input, tabix random access. INFO percent-decoding
is done with `urllib.parse.unquote` — 131 escaped values per 400k ClinVar records
and no bare `%`.

**Malformed input:** a data line not splitting into at least 8 tab-separated fields
raises, reporting the line number and what was expected. That is the whole
validation — not a VCF linter. **No whitespace fallback**, even though a
space-delimited file exists in this project's data
(`VEP_DNA/annotation/test.vcf`, which is also CRLF-terminated): accommodating it
would guess at intent. *Correction, 2026-08-21:* an earlier draft said `\r` must be stripped because a
well-formed CRLF VCF would otherwise corrupt every `info_` card. Measured false —
both `open()` and `gzip.open(..., "rt")` use universal newlines and translate
`\r\n` to `\n` before the parser sees the line. The stripping was dead code and the
test asserting it was vacuous; both were removed.

### Module placement, traversal, fixtures

`base_ops/`. `fixed_ops/` holds deterministic one-in-one-out transforms; `from_vcf`
turns a file into N states. `from_fasta` sits in `fixed_ops/` only because its
single-coordinate mode yields one fixed sequence, and `from_vcf` has no such mode.

`mode="sequential"`, hard-coded, matching `from_fasta`'s batch path. Not exposed.

Fixtures are synthetic and tiny — a few-kb FASTA and a ten-line VCF, `.fai`
generated at test time. **No reference genomes or real VCFs in this repository.**

### Out of scope for v1

`pool`/`region` splicing; N-padding to gene boundaries and per-variant strand (both
need annotation); VCF **output** (see `README.md`); `docs/operations/from_vcf.rst`
must be written, with entries in `source_operations.rst`.

## Verification at implementation time

Read-only, not copied into this repository:

- `KinneyLab/xin_collab/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa` (+ `.fai`)
- `KinneyLab/VEP_DNA/clinvar_hg38/clinvar.vcf.gz` — 3,375,801 records
- `KinneyLab/VEP_DNA/splicevardb_x_clinvar_snv.vcf` — 8,490 records, `chr1`-style

**Contig normalisation** is exercised by the *splicevardb* file against the Ensembl
FASTA — `chr1` against `1`. An earlier draft named the ClinVar/Ensembl pair for
this; both are `1`-style, so that pair tests nothing.

**SpliceAI reproduction** is a partial check only. SpliceAI's default `-D 50` gives
`wid = 10101`, so matching it needs flanks of **5050**, not 5000. Even then, byte
identity holds only for `len(REF) == 1` — SpliceAI holds the window fixed at `wid`,
so its ALT is `wid - len(REF) + len(ALT)`, while this design keeps `flank_right`
fixed. Verified: for `AGG>A` at flank 5050, SpliceAI emits 10101/10099 and this
design 10103/10101. Restrict the comparison to SNVs and insertions, and expect
further divergence from SpliceAI's gene-boundary N-padding and per-gene
reverse-complement, neither of which is adopted.
