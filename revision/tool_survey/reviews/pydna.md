# Adversarial review — pydna extraction

Reviewer pass date: 2026-08-10 (independent re-verification of repo, docs and PyPI).
Target: `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/extractions/pydna.md` + the structured extraction.

## Method (independent, not re-reading the extractor's quotes)

- Re-fetched from `raw.githubusercontent.com/pydna-group/pydna/master`: `README.md`, `pyproject.toml`, `LICENSE.txt`, and `src/pydna/{design,utils,codon,parsers,amplify,dseqrecord,assembly2,opencloning_models,genbank,seq,seqrecord,tm,contig,primer_screen,sequence_regex,sequence_picker,snapgene_history_parser}.py`.
- GitHub REST API was rate-limited (HTTP 403) from this host, so repo-liveness was verified independently via the commit Atom feed `https://github.com/pydna-group/pydna/commits/master.atom` (latest commit 2026-08-01T09:24:09Z) and PyPI JSON (`5.5.16`, uploaded 2026-06-09; classifiers `Development Status :: 4 - Beta`, `License :: OSI Approved :: BSD License`; `requires_python >=3.10,<4.0`).
- **Whole-docs keyword audit** (stronger than the extractor's code search): downloaded the Sphinx search index `https://pydna-group.github.io/pydna/searchindex.js` (258 kB), which stems every word of every docs page *including the full API reference*, and enumerated all 820 documented objects. Keyword counts: `barcode` 0, `gff` 0, `gtf` 0, `hgvs` 0, `vcf` 0, `degenerate` 0, `optimi*` 0, `mutagen*` 1 (example_gallery only — the external teemi link), `combinatori*` 2 (example_gallery + index), `exon` 1 page, `intron` 3 pages.
- Re-extracted the paper PDF with PyMuPDF and re-checked every quoted passage.

Result: the extraction is careful, unusually well sourced, and **factually accurate on every literal statement I could check**. Three cells are nonetheless scored too harshly, and all three fail *internal consistency with the Biopython extraction in this same survey* — which is the single most exploitable weakness for a referee, because pydna is built on Biopython.

---

## The three findings that matter

### 1. `exon_intron_split_codons = no` → should be `partial` (UNDERSTATED)

pydna's own getting-started documentation demonstrates a spliced CDS end-to-end.
`docs/notebooks/Dseq_Features.ipynb`, cell 12 (markdown), verbatim:

> "### Adding a Feature with Parts — To add a feature with parts, like a CDS with introns, we need to use a `CompoundLocation` object when creating a `SeqFeature`. The example code below adds a CDS with two parts, between 3-9bp and 12-15bp … In a real-world scenario this would represent a CDS with an intron that skips the `ACG` codon: ATGCGT~~ACG~~TGA"

Cell 13 output shows the GenBank feature `CDS join(4..9,13..15)`; cell 14: *"We can even extract a protein record as follows (see how the protein sequence is `MR`, skipping the intron)"*; cell 15 executes `dummy_record.features[-1].extract(dummy_record).translate()` → `ProteinSeq('MR')`.

And it is a first-class API, not just a notebook trick: `SeqRecord.extract_feature(n)` is literally `return self.features[n].extract(self)` (`src/pydna/seqrecord.py` ~line 383), inherited by `Dseqrecord.extract_feature()`; `utils.shift_location()` / `shift_feature()` explicitly handle `CompoundLocation`; GenBank `join(...)` locations round-trip through `parsers.parse()`.

Because exons are concatenated in transcript order before translation, a codon split across an exon junction *does* translate correctly. The extraction's sentence "nothing assembles a CDS across exons or handles codons split across an exon junction" is demonstrably false as written.

**Consistency problem:** `extractions/biopython.md` scores this key `partial` on exactly this mechanism ("CompoundLocation represents joined exons … so a codon split across an exon boundary translates correctly via `feature.extract(genome).translate()`"). pydna inherits that mechanism *and documents it in its own tutorial*. Scoring Biopython `partial` and pydna `no` is indefensible in a table a pydna author will read.

Suggested value: `partial`, with the same honest caveat used for Biopython — the splicing machinery is correct and documented, but pydna never places an edit into a junction-straddling codon because it performs no variant design.

### 2. `transcript_models = no` → should be `partial` (UNDERSTATED, same root cause)

The extraction's literal claims are all true: no GFF/GTF parser (confirmed — `parsers.py` handles GenBank/EMBL/FASTA + SnapGene only, and `gff`/`gtf` return 0 hits across the entire docs index), no transcript object, no isoform selection. But "annotation comes only from GenBank feature tables" *is* the transcript model: INSDC `mRNA`/`CDS` `join()` features parse into `CompoundLocation`, and pydna documents constructing, exporting and extracting them.

Again `extractions/biopython.md` scores this `partial (GTF/GFF specifically: no)` with the reasoning "which is a transcript model in all but file format". The same reasoning applies verbatim to pydna.

Suggested value: `partial (GFF/GTF: no)`.

### 3. `synthesis_constraints = no` → should be `partial` (UNDERSTATED, and the most likely author rebuttal)

`extractions/biopython.md` scores `synthesis_constraints = partial` on the grounds that "`Bio.Restriction` is a first-class restriction-site analyser … which covers the single most common synthesis/cloning constraint check (unwanted restriction sites). Add `gc_fraction()`, `MeltingTemp`, `molecular_weight()` …".

pydna exposes strictly *more* of that than Biopython does, as documented methods on its core objects (verified in `src/pydna/dseqrecord.py` lines 925-953 and in the docs object index):

- `Dseq.no_cutters() / unique_cutters() / once_cutters() / twice_cutters() / n_cutters(n) / cutters()` — defaulting to the `CommOnly` (commercially available) REBASE enzyme set;
- `Dseqrecord.number_of_cuts(*enzymes)`;
- `Seq.gc()` / `SeqRecord.gc()`;
- `pydna.tm` — `tm_default`, `tm_dbd`, `tm_product`, `ta_default`, `Q5`, `tmbresluc`, `tm_neb` (queries the NEB API), plus `program()` / `dbd_program()` which emit a suggested PCR program;
- `Seq.rarecodons()` returning the *positions* of rare codons.

"Does my designed construct contain an unwanted site for enzyme X, and which enzymes cut it zero/once/twice?" is answerable in one call. Nothing *enforces* a constraint and there is no homopolymer/repeat/hairpin/GC-window/oligo-length rule — so `partial`, not `yes`.

Suggested value: `partial`, with the sentence "site-content and Tm queries only; nothing is enforced and there is no repeat/homopolymer/GC-window/length rule" carried with the cell. If the authors prefer the strict reading (`no`), they must also change Biopython to `no`, or the table is self-contradicting.

---

## Smaller corrections (value stands, evidence needs softening)

- **`mixed_mutagenesis_one_pool = no`** — value correct. But the evidence sentence "pydna has no mutagenesis operation at all" overreaches and is rebuttable:
  - `assembly2.pcr_assembly(template, fwd, rvs, limit=14, mismatches=0)` and `PCRAssembly(..., mismatches=N)` accept **substitution mismatches between primer and template** (`primer_template_overlap`: "Maximum number of mismatches (only substitutions, no deletion or insertion)"; regex `"(" + query + "){s<=" + str(mismatches) + "}"`), i.e. pydna *simulates* mutagenic-primer PCR and the primer's mismatched bases end up in the product;
  - `assembly2.crispr_integration` / `homologous_recombination_integration` + the documented `Example_CRISPR` notebook ("repair it with an oligo") implement a designed single edit.
  Reword to: "pydna can *simulate* an individual edit (mismatched-primer PCR, CRISPR + oligo repair) but has no mutagenesis *design* operation — no scanning, no exhaustive singles/pairs, no random sampling, no WT-replicate concept, and no way to mix mutagenesis types into one pool."

- **`primer_design = yes`** — value correct but the stated scope limit ("assembly/amplification primers only") is too narrow. `src/pydna/primer_screen.py` also provides `forward_primers`, `reverse_primers`, `primer_pairs`, `flanking_primer_pairs`, `diff_primer_pairs`, `diff_primer_triplets`: given *a set of sequences* and a primer stock list, it selects primer pairs yielding **distinguishable product sizes per sequence** ("primers are selected that result in unique product sizes from each of the input sequences … could be used to verify genetic modifications such as cloning an insert into a plasmid vector"), with `short`/`long` size windows and a gel-resolvability `callback`. That is verification/genotyping primer selection across a collection, and it is one of the four headline docs examples. Keep `yes`; widen the scope sentence. The true limit is only that there is **no mutagenic-primer designer**.

- **`codon_optimization = no`** — value defensible (pydna exposes no optimizer; `codon.py` is pure data; `utils.express` is "Not yet implemented"). Flag for the manuscript: this same survey scores **Biopython `yes (basic)`** for `CodonAdaptationIndex.optimize()`, and biopython is a *hard dependency* of pydna (`biopython = "^1.87"`). Phrase the cell as "pydna itself provides codon-usage diagnostics (CAI, rare-codon localisation, start/stop frequency, N-end rule) but no codon-optimizing function", not as "pydna cannot codon-optimize".

- **`design_visualization = yes`** — supported, and if anything under-evidenced. The extraction misses `Contig.figure_mpl()` ("Graphic representation of the assembly … Returns matplotlib.figure.Figure", `src/pydna/contig.py` line 272) and `Contig.detailed_figure()`. So pydna has a genuine *graphical* assembly view, not only ASCII + gel PNG.

- **`per_sequence_provenance = yes`** — supported and, again, under-evidenced: `snapgene_history_parser.parse_snapgene_history()` "reads a `.dna` file, resolves every step recorded in its SnapGene history, and returns a `Dseqrecord` whose `source` attribute tree mirrors the full cloning provenance", i.e. pydna *imports* provenance from a commercial GUI as well as exporting it. `normalize_history()` is also present alongside `validate_history()`. The extraction's caveat (cloning provenance ≠ design-vs-WT provenance) is the right thing to print in the table.

- **`dag_chaining = partial`** — defensible, and I could not falsify it. I checked whether `CloningStrategy` can be *replayed on new inputs* (which would force `yes`): it cannot. `CloningStrategy` has only `from_dseqrecords()`, `add_dseqrecord()`, `add_primer()`, `reassign_ids()`, `model_dump[_json]()` — export only; re-execution lives in the OpenCloning backend, a separate project. `validate_history(recursive=True)` replays operations but only against the original inputs, for verification. Note for the table: Biopython scores `no` on this key, so pydna's `partial` is already a meaningful, defensible step up.

- **`combinatorial_motif_place = partial`** — leans *generous* to pydna (Biopython gets `no`), which is the safe direction. The `Assembly` graph genuinely enumerates all fragment orders and orientations, so `partial` is fine; keep the "homology-driven simulation of possible reaction products, no user-facing placement specification" clause attached.

- **`lazy_evaluation = no`** — supported, and consistent: Biopython earned `partial` for `SeqIO.parse`/`index`, whereas pydna's `parsers.parse()` is annotated `-> list[Dseqrecord | SeqRecord]` and is documented as "a greedy function, use carefully". Verified.

- **`genome_coordinates = partial`** — supported and consistent with Biopython's `partial`. One extra data point the extraction missed, which strengthens `partial` rather than changing it: `pydna.sequence_picker.genbank_accession()` BLASTs a sequence against NCBI `nt` and returns a record described as `"{accession} REGION: {start}..{stop}"` — i.e. sequence → genomic coordinates, the reverse direction.

- **`barcode_generation = no`** — supported. Closest thing in the codebase is `primer_screen.expand_iupac_to_dna()` (`expand_iupac_to_dna("ATNG") -> ['ATGG','ATAG','ATTG','ATCG']`), a cartesian expansion of a degenerate IUPAC oligo. It enumerates a degenerate pool but applies no GC or edit-distance constraint and has no attachment machinery, so `no` holds. Worth one clause in the evidence so an author cannot say "you missed our degenerate-oligo expansion".

- **`interface`, `license`, `maintained`, `assay_dms`, `assay_mpra`, `assay_insilico`, `hgvs_input`, `vcf_vep_output`, `consequence_annotation`, `library_as_object`, `automatic_naming`** — all independently re-verified, all supported. `pyproject.toml` declares no `[project.scripts]` (confirmed by reading the whole file); `LICENSE.txt` is BSD-3-Clause; commit feed and PyPI dates match the extraction exactly; the docs index shows zero occurrences of hgvs/vcf/gff/gtf/barcode across the entire rendered documentation including the API reference.

---

## Capabilities the extractor missed entirely

1. **`parse_snapgene_history()`** — reconstructs a full cloning-provenance tree from a SnapGene `.dna` file into pydna `source` objects (provenance *import*, via the `sgffp` dependency). No other tool in this survey does provenance interop with a commercial GUI.
2. **`Contig.figure_mpl()` / `Contig.detailed_figure()`** — matplotlib graphical rendering of a linear or circular assembly (the extraction credits only ASCII figures + gel images).
3. **Diagnostic / verification primer-pair selection over a set of sequences** — `primer_screen.diff_primer_pairs`, `diff_primer_triplets`, `flanking_primer_pairs`, `primer_pairs` with product-size windows and a gel-distinguishability callback.
4. **`expand_iupac_to_dna()`** — cartesian expansion of an extended-IUPAC (degenerate) oligo into every concrete DNA sequence it encodes.
5. **PCR simulation with mismatched primers** (`mismatches=N`) — substitution mismatches between primer and template, i.e. simulation of site-directed mutagenesis by primer.
6. **Restriction-site content queries as first-class methods** — `no_cutters/unique_cutters/once_cutters/twice_cutters/n_cutters/cutters/number_of_cuts` over the `CommOnly` enzyme set (see finding 3).
7. **Multi-part (spliced) CDS features** via `CompoundLocation`, documented in `Dseq_Features.ipynb`, with correct extraction and translation (see findings 1-2).
8. **`pydna.tm`** — several Tm models (`tm_default`, `tm_dbd`, `tm_product`, `Q5`, `tmbresluc`, `tm_neb` via the NEB web API) plus `program()`/`dbd_program()` emitting a suggested PCR cycling program.
9. **`utils.shift_location()` / `shift_feature()` / `Dseqrecord.shifted()` / `synced()`** — location arithmetic that preserves features (including `CompoundLocation`s) across the origin of a circular molecule; relevant because it is how design metadata survives circular permutation.
10. **`Dseqrecord.normalize_history()`** — alongside the already-credited `validate_history()`.

None of these changes the headline conclusion: pydna is a cloning-strategy simulation and documentation tool with no notion of a variant library. But items 1-3 and 6-7 are exactly the kind of thing a pydna author would list in a referee report as "the authors did not look at our software".

## Bottom line

Change three cells to `partial` (`exon_intron_split_codons`, `transcript_models`, `synthesis_constraints`) — all three are currently inconsistent with how the same mechanisms were scored for Biopython in this survey. Soften the evidence sentences for `mixed_mutagenesis_one_pool`, `primer_design` and `codon_optimization`. Everything else survives adversarial checking.
