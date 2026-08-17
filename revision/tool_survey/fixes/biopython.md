# Biopython — repair change log

**Record repaired:** `final/biopython.md` (edited in place)
**Audits processed:** `citation_audit/biopython.md` (18 findings), `factcheck/biopython.md` (A: 11, B: 15, C: 1)
**Repair date:** 2026-08-14
**Verification baseline:** tag `biopython-188` on github.com/biopython/biopython, read over `raw.githubusercontent.com`; PyPI JSON API; live biopython.org docs. Nothing installed or executed.

**Outcome tally:** 20 applied (across both audits; several findings are duplicates of each other and share one edit), 3 rejected, 20 escalated (17 of them the section-B/C coverage block).

**No capability value was changed.** All 24 `**\`key\` = value**` lines are byte-identical to the pre-repair record; only evidence and prose were touched.

---

## Applied

### CA-1 / FC-A4 — "all round-trip to GenBank" (`per_sequence_provenance`)

- **Verified:** `Bio/SeqIO/InsdcIO.py` (1556 lines) contains **zero** occurrences of `letter_annotations`. `record.annotations` is read only by specific INSDC keys (`date`, `data_file_division`, `references`, `structured_comment`, `comment`, `contig`, `keywords`, `segment`, `taxonomy`, `db_source`, `molecule_type`, …); `record.features` is written (lines 1166, 1521). Arbitrary annotation keys and per-letter data are therefore not serialised.
- **Correction:** `…plus \`letter_annotations\` for per-position data; all round-trip to GenBank.` → `…; \`features\` and the recognised INSDC \`annotations\` keys round-trip to GenBank, but arbitrary \`annotations\` keys and \`letter_annotations\` are not written by the GenBank writer (\`Bio/SeqIO/InsdcIO.py\` never reads \`letter_annotations\`).` `Bio/SeqIO/InsdcIO.py` appended to that entry's `*Source:*` line.

### CA-2 / FC-A7 (part) — bigMaf listed as an exon-block transcript model (`transcript_models`)

- **Verified:** `Bio/Align/bigmaf.py:32–57` declares `AutoSQLTable("bedMaf", "Bed3 with MAF block", [chrom, chromStart, chromEnd, mafBlock])`. No `blockCount`/`blockSizes`/`blockStarts`, no `thickStart`/`thickEnd`.
- **Correction:** `bigmaf` removed from the exon-block list at §3 `transcript_models`. bigMaf read/write support is real and is left asserted elsewhere in the record (§5 item 3); only the transcript-model attribution was removed.

### CA-3 / FC-A7 (part) — PSL credited with `thickStart`/`thickEnd`

- **Verified:** `Bio/Align/psl.py` — `grep -c -i thick` = **0**. PSL carries `blockSizes`/`qStarts`/`tStarts` (lines 173–193) but no CDS bounds. `Bio/Align/bigpsl.py:72,77,87,92` *does* carry `thickStart`, `thickEnd`, `blockCount`, `blockSizes`; `Bio/Align/bed.py:82–135` carries block arrays plus `thickStart`/`thickEnd`.
- **Correction:** three places.
  - §0 change table: `BED12/bigBed/PSL/bigPsl exon-block transcript models with thickStart/thickEnd CDS bounds` → `BED12/bigBed/bigPsl … — and PSL exon blocks without CDS bounds —`.
  - §3 `transcript_models`: exon-block list now reads `bigbed` and `bigpsl` — and, without `thickStart`/`thickEnd`, for `psl`.
  - §3 `exon_intron_split_codons`: `BED12/PSL exon blocks with thickStart/thickEnd` → `BED12 exon blocks with thickStart/thickEnd`.

### CA-4 / FC-A9 — "read *and write*" for the full `Bio.Align.formats` tuple (§5 item 1)

- **Verified:** `grep -c "^class AlignmentWriter"` over each of the 20 format modules at tag `biopython-188`: `emboss`, `hhr`, `msf`, `tabular` = 0; the other 16 = 1. `Bio/Align/__init__.py:4870–4877` — `write()` looks up `module.AlignmentWriter` and raises `ValueError(f"File writing has not yet been implemented for the {fmt} format")` on `AttributeError`. Tuple confirmed at lines 4812–4833.
- **Correction:** heading changed to `read for the full \`formats\` tuple, write for 16 of its 20 members`, with the four reader-only formats named and the `ValueError` behaviour stated.

### CA-5 — "dinucleotide-preserving shuffling" (§4 assessment)

- **Verified:** `Doc/Tutorial/chapter_cookbook.rst:119–203`. The recipe is `nuc_list = list(original_rec.seq)` then `random.shuffle(nuc_list)` then `Seq("".join(nuc_list))`. Mononucleotide composition only; no dinucleotide algorithm anywhere in the section.
- **Correction:** `(dinucleotide-preserving shuffling of an existing genome …)` → `(\`random.shuffle\` over the list of individual nucleotides of an existing genome, preserving mononucleotide composition only …)`.

### CA-6 — "Biopython no longer drives external command-line tools" (§6)

- **Verified:** in 1.88 source — `Bio/PDB/DSSP.py:176,189` `subprocess.Popen`; `Bio/PDB/NACCESS.py:65` `subprocess.Popen`; `Bio/PDB/PSEA.py:38` `subprocess.run`; `Bio/PDB/ResidueDepth.py:535` `subprocess.call(make_surface, shell=True)`; `Bio/Phylo/PAML/_paml.py:118,121` `subprocess.call`. The `Bio.Application` removal itself is confirmed verbatim at `DEPRECATED.rst:228–232`.
- **Correction:** `— Biopython no longer drives external command-line tools.` → `— Biopython no longer ships a generic framework for driving external command-line tools, though individual modules still launch executables directly via \`subprocess\` (\`Bio.PDB.DSSP\`, \`NACCESS\`, \`PSEA\`, \`ResidueDepth\`/MSMS, \`Bio.Phylo.PAML\`).`

### CA-7 / FC-A3 (part) — truncated `SeqRecord.__getitem__` quotation (`per_sequence_provenance`)

- **Verified:** `Bio/SeqRecord.py:433–437` reads *"With the exception of any molecule type, the annotations dictionary and the dbxrefs list are not used for the new SeqRecord, as in general they may not apply to the subsequence."* Lines 572–574 implement exactly that: `if "molecule_type" in self.annotations: answer.annotations["molecule_type"] = …`. The record had dropped the leading qualification while marking the quote verified verbatim.
- **Correction:** leading clause restored to the quotation; `so slicing silently drops annotations` → `so slicing silently drops annotations apart from \`molecule_type\``.

### CA-8 — Tutorial `.2bit` code block was a paraphrase (`genome_coordinates`)

- **Verified:** `Doc/Tutorial/chapter_align.rst:1591–1599`. The Tutorial assigns `filename = f"{name}.2bit"` separately, parses `filename`, asserts `len(genome_alignment.sequences[i]) == len(genome[chromosome])`, assigns the genome sequence, then sets `genome_alignment.sequences[i].id = f"{name}.{chromosome}"`.
- **Correction:** the fenced block now reproduces the Tutorial's five statements exactly (doctest prompts stripped). Nothing else in that entry was touched.

### CA-9 / FC-A1 — "Since 1.78" for undefined/partially-defined sequences (§2)

- **Verified:** `NEWS.rst` — 1.79 section (3 June 2021) at lines 554/614/627–644 introduces `_UndefinedSequenceData` and `Seq(None, length=6)`; the 1.78 section begins at line 689 and does not. `.defined` was added in **1.80** (`NEWS.rst:482–483`, "Sequences now have a `defined` attribute…").
- **Correction:** `Since 1.78` → `Since 1.79`, and `.defined` annotated `(added in 1.80)`.
- **Note:** `.defined_ranges` appears in `Bio/Seq.py` (lines 348, 2008, 2451) but has no `NEWS.rst` entry; no version was asserted for it and none was added.

### CA-10 — "19-item Tutorial cookbook" (§2 and `assay_dms`)

- **Verified:** `Doc/Tutorial/chapter_cookbook.rst` heading scan — 13 `~~~` subsections under "Working with sequence files", 4 under "Sequence parsing plus simple plots", plus the top-level BioSQL section = **18**. The record's own enumeration in §4 also totals 18.
- **Correction:** `19-item` → `18-item` in both places (§2 load-bearing paragraph, `assay_dms`).

### CA-11 — PyPI first-upload timestamp (`maintained`)

- **Verified:** PyPI JSON, release 1.88 — earliest `upload_time_iso_8601` is `2026-08-06T12:13:36.003082Z` for `biopython-1.88.tar.gz`; `2026-08-06T12:30:38.131147Z` is `biopython-1.88-cp310-cp310-macosx_11_0_arm64.whl`.
- **Correction:** `PyPI first file upload \`2026-08-06T12:30:38Z\`` → `\`2026-08-06T12:13:36Z\` (the sdist)`.
- **Not changed:** the "31 files" figure. Re-counted from the JSON: exactly 31 files had `upload_time_iso_8601 < 2026-08-11`, so it was correct on the record's 2026-08-10 finalization date. (PyPI now lists 41 after ten cp315/cp315t wheels on 2026-08-12 — matching the audit's own time-scoped note, which it flagged as "not a finding".)

### CA-12 — `Bio.motifs` "complete documented surface" (`combinatorial_motif_place`)

- **Verified:** `Bio/motifs/__init__.py` — beyond the listed members, the public `Motif` surface includes `__getitem__` (slicing, line 296), `mask` (221/224), `pseudocounts` (253/256), `background` (269/272), and attributes `name` (194), `alignment`, `counts`. The completeness claim is inaccurate.
- **Correction:** `Its complete documented surface is` → `Its documented surface includes`. The negative claims (no PWM sampler, no motif placer) are untouched and were re-confirmed against the full `def`/`@property` listing.

### CA-13 — package list omits four `[tool.setuptools].packages` entries (§2)

- **Verified:** `pyproject.toml:45–108`. `Bio.Align.substitution_matrices.data`, `Bio.Entrez.DTDs`, `Bio.Entrez.XSDs`, `Bio.SearchIO._model` are literal entries absent from the record's list.
- **Correction:** the list's lead-in now states the collapsing rule and names the four: `(\`[tool.setuptools] packages\`; data and internal subpackages — … — collapsed)`.

### CA-14 — the SeqIO format correction was itself incomplete (§0 and `transcript_models`)

- **Verified:** `Bio/SeqIO/__init__.py:414–453`, 38 direct `_FormatToIterator` keys. Beyond `twobit`, the original memo's list also omitted `abi-trim`, `fasta-blast`, `fasta-pearson`, `embl-cds`, `genbank-cds`, `sff-trim`.
- **Correction:** the six keys named in the `transcript_models` correction paragraph; §0 summary line now reads `(omitted \`twobit\`, six other direct \`_FormatToIterator\` keys, and the auto-registered \`AlignIO\` formats)`.

### CA-15 / FC-A10 — "~30 `Bio.SeqIO` formats read/write" listing read-only formats (§5 item 7)

- **Verified:** `_FormatToIterator` 38 direct keys, `_FormatToWriter` 18 direct keys (`Bio/SeqIO/__init__.py:414–475`); `Bio/AlignIO/__init__.py:153–176` adds 11 iterators and 8 writers, auto-registered into `SeqIO` at `Bio/SeqIO/__init__.py:509–519`. Totals: **49 readable / 26 writable**. `snapgene`, `gck`, `gfa1`, `gfa2`, `twobit` are absent from `_FormatToWriter` (read-only); `xdna` is present in both.
- **Correction:** item 7 now reads `49 readable \`Bio.SeqIO\` format identifiers and 26 writable (GenBank, EMBL, FASTQ variants and xDNA read/write; SnapGene, GCK, GFA1/2 and **twobit** read-only; plus auto-registered AlignIO formats)`.

### CA-16 — "`pip install biopython` therefore draws nothing" (`design_visualization`)

- **Verified:** `Bio/Graphics/__init__.py:9–23` — docstring and the `MissingPythonDependencyError` raise on missing ReportLab are exactly as quoted. But `Bio/motifs/__init__.py:467–566` `weblogo()` POSTs to `https://weblogo.threeplusone.com/create.cgi` via `urlopen` with no ReportLab involvement, and the record itself names that path twice elsewhere.
- **Correction:** `\`pip install biopython\` therefore draws nothing.` → `\`pip install biopython\` therefore renders no \`Bio.Graphics\` output.` The ReportLab point that cuts against the cell is preserved intact.

### FC-A2 — `Alignment` and `SeqIO.to_dict()` mischaracterised (`library_as_object`)

- **Verified:** `Bio/Align/__init__.py:999` `class Alignment`, constructed as `Alignment(sequences, coordinates)` where the docstring's own example passes three ungapped sequences of lengths 12/9/11 plus a coordinate array (lines 1040–1055) — not equal-length gapped rows. `Alignments` (line 4055) is a `list` subclass of such objects. `Bio/SeqIO/__init__.py:729` — `to_dict` is documented as "Turn a sequence iterator or list into a dictionary"; it materialises an ordinary in-memory dict, unlike `index()`/`index_db()`.
- **Correction:** clause (b) now distinguishes legacy `MultipleSeqAlignment` (equal-length gapped rows) from modern `Alignment` (original sequences + coordinate array) and `Alignments` (container); clause (c) now reads `the in-memory dict eagerly built by \`SeqIO.to_dict()\` together with the dict-like \`index()\` / \`index_db()\` views over a file`.
- **Left alone (tension noted):** the embedded quotation attributed to `Bio.SeqIO` — *"interprets multiple sequence alignment file formats as collections of equal length (gapped) sequences"* — is **not present** in 1.88's `Bio/SeqIO/__init__.py`, its module docstring, or the live `Bio.SeqIO` API page. Neither audit flagged it, so it was not touched. See escalation ESC-5.

### FC-A3 (part) — "Insertion … with features maintained" (`combinatorial_motif_place`)

- **Verified:** `Bio/SeqRecord.py:593` — a slice keeps a feature only when `start <= f.location.start and f.location.end <= stop`. A feature straddling the insertion point is dropped from both slices and cannot be recovered by `__add__`.
- **Correction:** `Insertion is \`rec[:i] + insert + rec[i:]\` with features maintained.` → `… with features maintained — but only those falling entirely within a slice survive, so a feature spanning position \`i\` is lost.`

### FC-A5 — obsolete ML/statistical inventory (`assay_insilico`)

- **Verified:** `Bio/HMM/` at tag `biopython-188` contains a single 283-byte `__init__.py` holding only a docstring. `DEPRECATED.rst:210–214` records that `Bio.HMM.DynamicProgramming`, `Trainer`, `MarkovModel`, `Utilities` were deprecated in 1.82 and removed in 1.86; lines 185–208 record the same for `Bio.kNN`, `Bio.LogisticRegression`, `Bio.NaiveBayes`, `Bio.MaxEntropy`, `Bio.MarkovModel`.
- **Correction:** `Biopython's only ML/statistical content is the legacy \`Bio.HMM\`, \`Bio.Cluster\`, and \`Bio.PopGen\`` → `… is \`Bio.Cluster\` and \`Bio.PopGen\` … ; the \`Bio.HMM\` implementation modules (…), together with \`Bio.kNN\`, \`Bio.LogisticRegression\`, \`Bio.NaiveBayes\` and \`Bio.MaxEntropy\`, were removed in 1.86, leaving \`Bio.HMM\` an empty package.` The cell's actual ground (no sequence-to-function model interface) is unchanged.

### FC-B1 (first half) — `MutableSeq` reduced to per-position assignment (`mixed_mutagenesis_one_pool`)

- **Verified:** `Bio/Seq.py` `MutableSeq` — `__setitem__` (2220, item **and** slice), `__delitem__` (2242), `append` (2253), `insert` (2265), `pop` (2280), `remove` (2299), `reverse` (2318), `extend` (2325); inherited `replace` (1951), `join` (1913), `__add__` (526), `__mul__` (560).
- **Correction:** `The closest primitive is per-base assignment on a \`MutableSeq\`` → `The closest primitives are \`MutableSeq\`'s in-place edits — item and slice assignment/deletion, \`insert\`, \`append\`, \`extend\`, \`pop\`, \`remove\`, plus the shared \`replace\`, \`join\`, \`+\` and \`*\``. The verbatim docstring quote and the cell's `no` reasoning are unchanged.
- **Second half deferred:** the finding also asks for the "no molecule/alphabet validation" hazard to be added. That is new coverage rather than a correction to an existing statement, so it rolls into ESC-3.

---

## Rejected

### CA-18 — license quotation not character-for-character verbatim (`license`)

- **Why rejected:** the discrepancy is real but immaterial and the "fix" would degrade the text. `LICENSE.rst:1–8` does use double quotation marks around "Biopython License Agreement" and "BSD 3-Clause License"; the record renders them as single marks because the whole passage is already inside a double-quoted italic span, where nested double quotes are ambiguous. **The record makes no "verbatim" claim at that location** (unlike the entries in §3 that do), and the audit itself states the meaning is unchanged. Changing it would introduce nested `"…"…"…"` for zero factual gain.
- **Evidence:** `LICENSE.rst:1–8` read directly; record line for `license` contains no "verbatim" marker (the §6 method caveat scopes that claim to `Quoted strings marked "verbatim"`).

### FC-A8 — `CodonAligner` misattributed split-codon functionality (`exon_intron_split_codons`)

- **Why rejected:** the audit's reading is stronger than what the record says. The record's clause is `plus \`Bio.Align.CodonAligner\` (confirmed present in \`Bio/Align/__init__.py\`) **for codon-aware alignment**` — it attributes codon-aware alignment to `CodonAligner`, not exon-model interpretation or intron removal. That attribution is correct. The genuinely wrong part of the same sentence (PSL credited with `thickStart`/`thickEnd`) was corrected under CA-3.
- **Evidence:** `Bio/Align/__init__.py:4604` — `class CodonAligner(_codonaligner.CodonAligner)`, docstring *"Aligns a nucleotide sequence to an amino acid sequence. This class implements a dynamic programming algorithm to align a nucleotide sequence to an amino acid sequence."* `score(seqA, seqB)` takes a protein and a nucleotide sequence. The record's description matches; the audit's paraphrase of the record does not.

### FC-A11 — "the linked 'latest' documentation was not 1.88" (§1 sources table)

- **Why rejected:** the finding is wrong about the URL the record actually cites. The record labels `https://biopython.org/docs/latest/Tutorial/index.html` as "Tutorial & Cookbook, 1.88" — and that page, fetched live on 2026-08-14, contains the string **"Biopython 1.88"**. What still shows 1.87 is a *different* page, `https://biopython.org/wiki/Documentation` (and the homepage), which offers a `1.87.zip` docs download. The record already states that lag correctly under `maintained` ("biopython.org's homepage still advertises 1.87").
- **Evidence:** `curl https://biopython.org/docs/latest/Tutorial/index.html` → matches "Biopython 1.88"; `curl https://biopython.org/wiki/Documentation` → matches "1.87.zip", "1.87"; `curl https://biopython.org/` → "1.87".

---

## Escalated (record left untouched at these points)

### ESC-1 (CA-17) — the load-bearing absence claim vs. `MeltingTemp.Tm_GC(valueset=2)`

- **Fact verified:** `Bio/SeqUtils/MeltingTemp.py:725–727` documents valueset 2 as the *"'QuikChange' formula. Recommended (by the manufacturer) for the design of primers for QuikChange mutagenesis."*
- **Why escalated:** the sentence at issue is flagged in the record itself as **"the load-bearing claim for the comparison table"** — that nothing in Biopython is "concerned with designing … a mutagenesis scheme". Whether a Tm formula tuned for QuikChange primers falsifies "concerned with", or is merely a calculation the claim already tolerates, decides how much of Block A/B rests on that sentence. The audit itself concedes "the narrower functional claim elsewhere in the record remains supported".
- **Question for authors:** does the load-bearing sentence need a carve-out for `Tm_GC(valueset=2)`, or is it left absolute?
- **Options:** (a) leave as is; (b) narrow to "no … operation that *generates* a variant library, mutagenesis scheme, oligo pool or combinatorial construct"; (c) append a one-clause exception naming `Tm_GC(valueset=2)` as a Tm formula only.

### ESC-2 (FC-A6) — `Bio.Alphabet` stub and `Bio.PopGen` removals left unflagged

- **Facts verified:** `Bio/Alphabet/__init__.py` is a docstring plus an unconditional `raise ImportError("Bio.Alphabet has been removed from Biopython…")`, yet `Bio.Alphabet` is a live entry in `[tool.setuptools].packages` (`pyproject.toml:52`). `pyproject.toml` lists only `Bio.PopGen` and `Bio.PopGen.GenePop` — no FDist, no SimCoal.
- **Why escalated:** nothing the record *states* is false — the package list is literally accurate and it is presented as a package list, not a capability list. The finding is that an omission is misleading, which is exactly the materiality judgment reserved for the authors. It also touches §5 item 12's "population genetics (`Bio.PopGen`)" line.
- **Question for authors:** should §2's package list and §5 item 12 carry a "these are stubs / reduced since the paper" annotation?
- **Options:** (a) leave; (b) annotate `Alphabet` in the package list as an import-failing stub removed in 1.78; (c) also note in §5 item 12 that `Bio.PopGen` is GenePop parsing only since the 1.70 FDist/SimCoal removal.

### ESC-3 (FC-B2 … FC-B15) — the fourteen remaining "incomplete" findings

Not applied. Each asks the record to cover a capability area it currently does not discuss at all, rather than to correct a statement it makes. Every one of them would add positive capability evidence to a cell already scored `no` or `partial`, which is precisely "a correction that would change what a capability value rests on".

| # | Finding | Cell it would land in |
|---|---|---|
| B2 | CDS validation (`translate(cds=True)`), `CodonTable` alt codes, naive `back_table` | `codon_optimization` = yes (already the 2nd-most-debatable cell, §6.3) |
| B3 | `Seq.search` multi-pattern, `nt_search` IUPAC-aware | `synthesis_constraints` = partial |
| B4 | per-letter qualities; ABI/`abi-trim`/PHD/ACE/SFF/FASTQ/QUAL interfaces | §4 / `lazy_evaluation` |
| B5 | `Alignment` slicing, sorting, concat, frequencies, substitutions, `from_alignments_with_same_reference` | §5 / `library_as_object` = no |
| B6 | `Bio.codonalign.build` (experimental) | §5 |
| B7 | `Bio.CAPS.CAPSMap(...).dcuts` detail | `consequence_annotation` = no |
| B8 | GenePop split/remove/stream operations | §5 item 12 |
| B9 | `Bio.Phylo.Consensus.bootstrap` collection generation | §5 |
| B10 | motif thresholds, pseudocounts, score distributions, comparison | `combinatorial_motif_place` = no (its completeness claim was fixed under CA-12) |
| B11 | digest/topology/overhang/isoschizomer modes; REBASE 404 (2024), 1,088 enzymes | `synthesis_constraints` = partial (3rd-most-debatable cell, §6.4) |
| B12 | Tm model breadth and caller-supplied context | `primer_design` = no |
| B13 | `ProteinAnalysis` / GC / molecular-weight descriptor detail | §5 item 12 |
| B14 | `Bio.PDB.PICIO.read_PIC_seq` + `internal_coords` | §5 item 11 |
| B15 | `Bio.Blast.qblast`, `ExPASy.ScanProsite.scan` | §5 item 10 |

- **Question for authors:** how much capability breadth should a *design-tool* survey record carry for an infrastructure library? Adding B11/B12 in particular would strengthen the evidence under two cells the record already flags as debatable.
- **Options:** (a) leave — the record's scope is design capability, not an API census; (b) add a single "screening and validation primitives" paragraph in §5 covering B2–B4, B10–B13 without touching any Block A–D entry; (c) work each into its cell (this would put `synthesis_constraints` and `primer_design` back in play and must precede, not follow, the ratings-reconciliation pass).

### ESC-4 (FC-C) — balance

- The fact-check's section C is a single judgment: the record "materially understates what current Biopython can do as infrastructure". Explicitly reserved for the authors; no edit made.
- **Question for authors:** accept the balance criticism and rebalance (which means resolving ESC-3 first), or record that the imbalance is deliberate given the survey's design-tool scope?
- **Options:** (a) declare the asymmetry intentional in §6 and close; (b) rebalance per ESC-3 option (b); (c) rebalance per ESC-3 option (c) and re-run the ratings pass.

### ESC-5 (incidental — from neither audit) — an unlocatable quotation in `library_as_object`

- Found while verifying FC-A2. The record attributes to `Bio.SeqIO` the quotation *"interprets multiple sequence alignment file formats as collections of equal length (gapped) sequences"*. It is **not** in `Bio/SeqIO/__init__.py` at tag `biopython-188` (no "equal length" match anywhere in the file), not in `Bio/AlignIO/__init__.py`, not in `Doc/Tutorial/chapter_seqio.rst`, and not on the live `Bio.SeqIO` API page. It reads like wording from an older release.
- Left untouched: it is outside both audits' findings, and removing it would be an unaudited edit. Flagged here because the surrounding sentence *was* edited under FC-A2, so a reviewer will encounter it.
- **Question for authors:** locate the source release for the quote and version-qualify it, or drop it?

---

## Verification notes

- Every code claim was checked against tag `biopython-188`, not master, since the record's assessed version is 1.88. Where the record cites master line numbers, spot checks matched the tag: `Bio/Align/__init__.py:3396` (liftover docstring), `:4818` (`"chain"` in `formats`), `Bio/Align/bigbed.py:1042` (`search`, including the stray double period the record notes), `Bio/SeqIO/__init__.py:450` (`twobit`), `:509–519` (AlignIO auto-registration), `Bio/Seq.py:2173` (`MutableSeq` docstring), `Bio/SeqRecord.py:433–437`, `DEPRECATED.rst:228–232`, `Bio/Graphics/__init__.py:9–23`, `Bio/Align/bed.py:7–12`.
- No file outside `revision/tool_survey/` was written. Nothing was installed or executed; all source was read over HTTPS.

## Pass 2 — policy application

**Date:** 2026-08-14

**Baseline and counts.** Every historical open item was rechecked against Cock et al. 2009 and official Biopython sources at tag `biopython-188`; current-release and availability checks used PyPI plus the official homepage and Tutorial. The PDF was read directly with PyMuPDF. Counts use the 21 atomic policy items below: **4 applied · 16 declined-by-policy · 1 verified/no edit required · 0 rejected-unverifiable · 0 remaining record escalations.** No capability value changed.

1. **ESC-1 / CA-17 — QuikChange wording: applied.** `Bio/SeqUtils/MeltingTemp.py:719–727` verifies that `Tm_GC(valueset=2)` is labelled the QuikChange formula and recommended for mutagenesis-primer design. The implementation accepts a supplied sequence and returns only its calculated Tm; it does not select a primer or generate a mutagenesis scheme. **Edit:** narrowed the load-bearing absence sentence from anything “concerned with designing” to anything that **generates** a library/scheme/pool/construct, and added the one-clause QuikChange distinction to the existing `primer_design` entry. The existing `mixed_mutagenesis_one_pool = no` and `primer_design = no` values remain supported.

2. **ESC-2 / FC-A6 — `Bio.Alphabet` and historical `Bio.PopGen`: applied.** `Bio/Alphabet/__init__.py` documents removal in 1.78 and raises `ImportError` unconditionally. `DEPRECATED.rst:619–630` documents removal of `Bio.PopGen.FDist` and `Bio.PopGen.SimCoal` in 1.70; `pyproject.toml` now packages only `Bio.PopGen` and `Bio.PopGen.GenePop`. **Edit:** brief annotations were added to the existing §2 package inventory and §5 item 12. No capability value changed.

3. **ESC-3 / FC-B1 remainder — absent molecule/alphabet validation: declined-by-policy.** Verified from the `Bio.Alphabet` stub and `Bio.Seq.Seq.translate` documentation, which says translating a protein sequence is accepted but biologically meaningless. Table 1 already describes general sequence manipulation and does not claim molecule-type validation; no Purpose, Key features, Output, or Availability cell changes.

4. **ESC-3 / FC-B2 — CDS validation and codon tables: declined-by-policy.** `Bio/Seq.py:1520–1635` verifies that `translate(cds=True)` checks a valid start, length divisible by three, and a single terminal in-frame stop; `Bio/Data/CodonTable.py:149–161` explicitly calls `make_back_table` a naive single-codon mapping and exposes alternative/ambiguous tables. These are details within Table 1's existing “transformation” description. The `codon_optimization = yes` value continues to rest independently on `CodonAdaptationIndex.optimize()`.

5. **ESC-3 / FC-B3 — exact/IUPAC sequence search: declined-by-policy.** `Bio.Seq.Seq.search` verifies multi-subsequence exact search; `Bio/SeqUtils/__init__.py:277–290` verifies forward-strand IUPAC query expansion by `nt_search`. This does not make the generic Table 1 key-feature cell inaccurate and does not change `synthesis_constraints = partial`.

6. **ESC-3 / FC-B4 — quality annotations and read formats: declined-by-policy.** `Bio/SeqRecord.py:262–310` verifies per-letter quality storage and slicing; the 1.88 `Bio.SeqIO` registry/docs verify ABI/`abi-trim` (Mott), ACE, PHD, SFF/`sff-trim`, FASTQ variants, and QUAL. Table 1 already says “Parsing, transformation and file-format support,” so no cell needs another format list.

7. **ESC-3 / FC-B5 — alignment-collection manipulation: declined-by-policy.** `Bio/Align/__init__.py` verifies row/column slicing, sorting, reverse complement, concatenation, frequencies, substitutions, counts, and `Alignment.from_alignments_with_same_reference`. These operations are already encompassed by Table 1's transformation wording; none creates the library abstraction that the same cell explicitly denies.

8. **ESC-3 / FC-B6 — experimental codon-alignment builder: declined-by-policy.** `Bio/codonalign/__init__.py:21–49` emits `BiopythonExperimentalWarning` and defines `build(pro_align, nucl_seqs, ...)` to return a `CodonAlignment`. This niche, experimental transformation does not alter a Table 1 cell.

9. **ESC-3 / FC-B7 — CAPS discrimination: declined-by-policy.** `Bio/CAPS/__init__.py:20–129` verifies equal-length alignment checking and `CAPSMap.dcuts` entries containing position, enzyme, sequences cut, and sequences blocked. CAPS detects restriction differences; it is not the molecular-consequence classifier tested by `consequence_annotation`, and it does not alter the four Table 1 cells.

10. **ESC-3 / FC-B8 — GenePop collection operations: declined-by-policy as a separate omission.** `Bio/PopGen/GenePop/__init__.py:97–223` verifies reconstruction, split-by-population/locus, and removal operations; `FileParser.py` verifies sequential large-file access. The short “GenePop parsing/manipulation only” description added under ESC-2 supplies the required current-package context; the omitted method census would not change Table 1.

11. **ESC-3 / FC-B9 — bootstrap collections: declined-by-policy.** `Bio/Phylo/Consensus.py:541–614` verifies the obsolete alignment-replicate generator plus current bootstrap-tree/consensus functions; `Bio/Nexus/Nexus.py:1932` verifies its bootstrapped-matrix method. Statistical resampling is not library design and changes none of the Table 1 cells.

12. **ESC-3 / FC-B10 — motif scoring controls: declined-by-policy.** `Bio/motifs/__init__.py` verifies configurable pseudocounts/background and motif slicing; `Bio/motifs/matrix.py:218–285,400–417` verifies consensus controls, the unimplemented substitution-matrix branch, and DNA-only PSSM calculation; `Bio/motifs/thresholds.py` verifies FPR/FNR/balanced/Patser thresholds. These are analysis controls, not motif placement, and do not alter Table 1.

13. **ESC-3 / FC-B11 — restriction digest/topology breadth: declined-by-policy.** `Bio/Restriction/Restriction.py` verifies linear/circular search, `catalyse`/`catalyze`, end types, compatible ends, isoschizomers, suppliers, and `Analysis` filters. The official 1.88 dictionary identifies REBASE release 404 (2024) and contains 1,088 enzyme records. This strengthens the already documented restriction-analysis primitive but changes neither `synthesis_constraints = partial` nor a Table 1 cell.

14. **ESC-3 / FC-B12 — Tm model breadth: declined-by-policy.** `Bio/SeqUtils/MeltingTemp.py` verifies `Tm_Wallace`, `Tm_GC`, and `Tm_NN`, custom tables, salt/ion/dNTP handling, mismatch/dangling-end inputs through caller-supplied `c_seq`/`shift`, `selfcomp`, and separate DMSO/formamide correction. Only the QuikChange wording needed a factual carve-out, applied under ESC-1; the broader API census does not change `primer_design = no` or Table 1.

15. **ESC-3 / FC-B13 — DNA/protein descriptors: declined-by-policy.** `Bio/SeqUtils/ProtParam.py` verifies composition, molecular weight, aromaticity, instability, flexibility, GRAVY scales, pI/charge, extinction coefficient, and coarse secondary-structure fraction. `Bio/SeqUtils/__init__.py` verifies GC ambiguity policies, `GC123`, GC skew, and molecular-weight strand/topology/monoisotopic modes. These analysis details do not change the generic Table 1 cells.

16. **ESC-3 / FC-B14 — sequence-to-default protein coordinates: declined-by-policy.** `Bio/PDB/PICIO.py:818–860` verifies `read_PIC_seq` returning a `Structure` with default internal coordinates, while the same file warns that defaults degrade quickly away from true coordinates; `Bio/PDB/internal_coords.py` verifies internal↔atom coordinate rebuilding and geometry edits. This specialist representation does not change Table 1's purpose, key-feature summary, output, or availability.

17. **ESC-3 / FC-B15 — remote screening entry points: declined-by-policy.** `Bio.Blast.qblast` is documented in the 1.88 source as an NCBI QBLAST call; `Bio.ExPASy.ScanProsite.scan` calls the ScanProsite endpoint; `Bio.SearchIO` officially parses BLAST, HMMER, Infernal, Exonerate, and InterProScan outputs. These network/parser details do not alter the existing Table 1 cells.

18. **ESC-4 / FC-C — balance and emphasis: declined-by-policy.** The authors' ruling expressly declines rebalancing. No edit.

19. **ESC-5 — unlocated `Bio.SeqIO` quotation: applied.** The prescribed PDF extraction located the sentence in Cock et al. 2009 p.1423: “Bio.SeqIO interprets multiple sequence alignment file formats as collections of equal length (gapped) sequences.” **Edit:** explicitly attributed and version-qualified it as the 2009 paper's description and added the paper to that entry's source line.

20. **Policy C — provenance and ambiguous labels: applied.** Removed the three sibling-analysis rows from §1, replaced master-branch references with tag `biopython-188`, and defined “current” to mean that tag absent a retrieval date. Removed reliance on the separate HGVS package and named GFF packages, replaced the unverifiable “never merged” Mapper history with its verified absence from the 1.88 tree, removed the code-author/referee anecdote and GitHub-rate-limit hedge, and replaced the unsupported “single most common” constraint claim. Current PyPI facts were fetched rather than hedged: 1.88 still governs, with 41 files as of 2026-08-14; the official homepage/Tutorial version split is now date-scoped.

21. **Policy D — version drift: verified/no edit required.** PyPI still reports 1.88 as current, and the official latest Tutorial identifies itself as 1.88. There is no materially different current release, so no version-drift parenthetical was added. The header version continues to govern.

**Value-basis result.** Every current value remains supported by admitted primary evidence. In particular, a Tm formula does not design primers; CAPS does not call molecular consequences; motif scoring does not place motifs; alignment/GenePop/bootstrap containers do not supply a library-design abstraction. **No value change or downstream value escalation was made.**

**Escalations:** none. **Rejected-unverifiable findings:** none. Every factual audit item above was independently verifiable; policy-excluded omissions were declined rather than rejected.

**Row-substitution candidates:** none established. The locked draft already flags `Codon / amino-acid substitutions` (possible near-uniform support) and `Recombination and chimeras` (possible one-tool support), but Biopython has no Table 2 column and therefore cannot resolve either distribution. Its strongest alternative axis, sequence/genomic-interval/alignment interoperability (49 readable SeqIO identifiers, 26 writable identifiers, 20 readable `Bio.Align` formats, liftover, and indexed region query), would make Biopython a ceiling only in a table that included Biopython; no verified eight-tool evidence shows it beats a named incumbent under the locked substitution protocol. No replacement row is nominated or applied.

**Neighbouring tension logged.** The historical `## Escalated` section above still says these points were untouched; this Pass 2 section supersedes those outcomes and preserves the earlier text as an audit trail. §0 retains extraction/review history but those sibling records are no longer listed or used as factual sources. PyPI's 41-file count differs from the historically correct 31 files at the record's 2026-08-10 finalization because ten cp315/cp315t wheels were uploaded on 2026-08-12.
