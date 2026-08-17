# pydna — v2 capability record

**Slug:** `pydna`
**Full name:** pydna — "a simulation and documentation tool for DNA assembly strategies using Python"
**Citation key:** `Pereira2015wj`
**Version assessed:** PyPI 5.5.16 (2026-06-09); `master` at `6.0.0-a.24.post` (head commit 2026-08-01)
**Basis:** re-scored from `final/pydna.md` (FINAL, adversarially reviewed) against the v2 row set in
`ROWS_v2.md`. No new extraction was performed; two supporting files (`extractions/pydna.md`,
`reviews/pydna.md`) were re-read to resolve the four new Block D rows.

**Scope note (important, and not a criticism).** pydna is a general-purpose molecular-cloning simulation
and documentation toolkit. Only functionality that participates in *specifying a DNA sequence library* is
scored here. Gel rendering, format I/O, GenBank download/fixing, `seguid` checksums, Tm models, CRISPR
protospacer search, restriction bookkeeping and reaction simulation are described where they touch a row,
but they are not allowed to inflate any cell. pydna's own framing of "combinatorial" means *multi-fragment
assembly of one intended construct*, not variant libraries — that is the tool's own wording, and the
record follows it.

---

## Block A — library specification

**`library_first_class_object` — no.**
There is no library / pool / collection class anywhere in the documented API. Verified against the Sphinx
`searchindex.js`, which enumerates all **820 documented objects** — no collection type appears. The unit of
work is a single `Dseqrecord` (one dsDNA molecule). Operations that can yield more than one molecule return
**plain Python containers**: `Assembly.assemble_linear(...) -> list[Dseqrecord]`, `assemble_circular(...) ->
list[Dseqrecord]`, `Dseqrecord.cut(*enzymes) -> Tuple[Dseqrecord, ...]`, and those elements are *alternative
products of one reaction*, not a designed set.
**Fairness clause that must travel with this cell:** pydna is emphatically *not* a file-writing tool — the
individual molecule (`Dseqrecord`) is a rich first-class object the user holds, inspects, transforms and
passes onward, which is the better half of what this row asks. What is absent is the *library* level: no
object represents a set-as-a-designed-thing, so building many variants means writing your own loop over
plain `list`s. Scored `no` on the row as written (the object must be the library).
*Source:* `src/pydna/assembly2.py`, `src/pydna/dseqrecord.py`, docs `searchindex.js` object index, paper abstract.
*(v1 `library_as_object` = no; this half of the split is unchanged.)*

**`composable_operations` — yes. [CHANGED from `partial` under v1 `dag_chaining`]**
The v1 `partial` was driven by the *DAG / declarative-pipeline* criterion in the old row wording ("can a
serialized strategy be re-executed on new inputs?" — it cannot; `CloningStrategy` is export-only). The v2
row asks a different and narrower question: do design steps compose/chain/nest arbitrarily, or is the
pipeline fixed by the tool? On that question pydna is unambiguous: **there is no pipeline at all**, only
ordinary Python methods, each of which returns a `Dseqrecord` that every other operation accepts —
`pcr(...)` → `.cut(...)` → `.looped()` → `gibson_assembly([...])` → `.shifted()` → `.cut(...)`, nested and
reordered freely, across ~15 assembly back-ends (`gibson_assembly, in_fusion_assembly, fusion_pcr_assembly,
in_vivo_assembly, restriction_ligation_assembly, golden_gate_assembly, ligation_assembly, gateway_assembly,
homologous_recombination_integration/excision/inversion, cre_lox_integration/excision,
recombinase_assembly`, plus `PCRAssembly`/`SingleFragmentAssembly`). Composition is also *recorded*: every
derived record carries a typed `source` pointing at its inputs and `Dseqrecord.history()`
(dseqrecord.py L1521) prints the tree:

```
>>> product, *_ = gibson_assembly(fragments, limit=4)
>>> print(product.history())
╙── product_name (Dseqrecord(o34))
    └─╼ GibsonAssemblySource
        ├─╼ fragment1 (Dseqrecord(-21))
        ├─╼ fragment2 (Dseqrecord(-12))
        └─╼ fragment3 (Dseqrecord(-13))
```

**Clause that must travel with the cell:** *what composes is construct-building operations on one molecule
at a time, executed eagerly; there is no parameterised specification that fans out over a variant set.* The
v1 caveats survive and are captured by `lazy_generation = no`, `library_algebra = no` and
`library_first_class_object = no` rather than by discounting this row.
*Source:* `src/pydna/assembly2.py`, `src/pydna/dseqrecord.py` L1521/L1554/L1580, `src/pydna/opencloning_models.py`.

**`lazy_generation` — no.** *(rename of `lazy_evaluation`, value carried unchanged)*
Execution is eager: "Execution of the script automatically yields the sequence of the final molecule(s) and
that of any intermediate constructs" (abstract). `parsers.parse()` is annotated `-> list[Dseqrecord |
SeqRecord]` and documented as "a greedy function, use carefully";
`Assembly.assemble_linear/assemble_circular/assemble_insertion` return `list[Dseqrecord]` — all products
materialised. Internal generators (`utils.limit_iterator`, utils.py L920, wrapping
`nx.cycles.simple_cycles(self.G)` capped at 10,000 in `assembly2.py`) are **combinatorial-explosion guards
on graph path enumeration**, not user-facing lazy sequence generation.
*Consistency:* Biopython earns `partial` here for `SeqIO.parse`/`index`/`index_db` lazy *parsing*; pydna
deliberately does not expose that, so `no` is correct and the two do not conflict.
*Source:* `src/pydna/parsers.py`, `src/pydna/assembly2.py`, `src/pydna/utils.py` L920, paper abstract.

**`library_algebra` — no.** *(new half of the `library_as_object` split)*
There is no library object, therefore no stack / concat / sample / repeat operation over libraries. Nothing
in the 820-object API index combines two sets of designed sequences, samples n from a set, or replicates a
set. Multi-molecule results are plain `list`/`tuple`, and any pooling, concatenation or sampling is done by
the user in their own Python — i.e. outside the tool, which the row defines as `no`. The one operation that
*consumes* many sequences at once, `Assembly(fragments)`, joins parts into single constructs by homology; it
is fragment assembly, not library algebra.
*Source:* `src/pydna/assembly2.py`, `src/pydna/dseqrecord.py`, docs `searchindex.js`.

**`exhaustive_single_scans` — no.** *(from the `mixed_mutagenesis_one_pool` split)*
No substitution/deletion/insertion scanning of any kind: no operation walks a sequence and emits every
variant at every position. `mutagen*` occurs on exactly **one** page of the entire rendered documentation
(`example_gallery`, and there only as an outbound link to the third-party `teemi` repo notebook "02 Primer
directed mutagenesis with pydna"). pydna *can simulate* an individual designed edit — `assembly2.pcr_assembly
(template, fwd, rvs, mismatches=N)` (L3244) / `PCRAssembly(frags, mismatches=N)` (L1951), whose
`primer_template_overlap()` (L697) documents `mismatches : int — "Maximum number of mismatches (only
substitutions, no deletion or insertion)"` and implements it as the fuzzy regex `"(" + query + "){s<=" +
str(mismatches) + "}"`; and `crispr_integration` / `homologous_recombination_integration` with the
`Example_CRISPR` notebook ("delete genes by making two cuts … repair it with an oligo"). That is simulation
of one user-specified edit, not enumeration of a scan.
*Source:* `src/pydna/assembly2.py` L697-L790, L1947-L1988, L3244-L3280; `docs/example_gallery.rst`; docs `searchindex.js`.

**`sampled_random_mutagenesis` — no.** *(from the split)*
No random or rate-based variant sampling. `utils.randomDNA/randomRNA/randomORF/randomprot` (utils.py
L444-L538) generate *de novo* random sequences from `length` (and optional `maxlength`) only — they do not
take a parent sequence and do not apply a per-base mutation rate, so they are not variant sampling. No error-
prone-PCR, doped-oligo or mutation-rate parameter exists anywhere in the API index.
*Source:* `src/pydna/utils.py` L444-L538, docs `searchindex.js`.

**`higher_order_combinatorial` — no.** *(from the split)*
Nothing enumerates pairwise or higher-order combinations of mutations in one sequence. The one place where
several changes land in one molecule is mismatched-primer PCR (`mismatches=N` permits up to N substitutions
in the annealing region), but the mismatches are wherever the *user's own primer* happens to differ from the
template — pydna neither designs them nor enumerates their combinations. The `Assembly` graph enumerates
combinations of **fragments**, not of mutations.
*Source:* `src/pydna/assembly2.py` L697-L790, L1947-L1988.

**`heterogeneous_components_one_library` — no.** *(from the split)*
There is no library specification into which anything could be mixed. `Assembly` accepts a heterogeneous list
of fragments, but that is the parts list for one construct (with all mechanistically possible products
enumerated as a strategy-error check), not a library in which structurally different component types coexist
as designed sub-populations. No WT-replicate concept, no per-component typing, no pool.
*Source:* `src/pydna/assembly2.py` (`Assembly` docstring), `src/pydna/dseqrecord.py`, docs `searchindex.js`.

**`combinatorial_motif_place` — partial.** *(carried from v1, unchanged)*
`Assembly` builds a directed graph over fragments in which "The sign of the node key represents the
orientation of the fragment, positive for forward orientation, negative for reverse orientation"
(`Assembly` docstring), then enumerates **all** linear / circular / insertion paths through it. From the
`new_assemblies` notebook: "there are two possible circular assemblies, one with all fragments oriented as
they were passed to the class constructor, and one with the second fragment inverted." So a set of parts
sharing homology genuinely is combined in **all valid orders and orientations** — a real combinatorial
enumeration.
Why not `yes`: (i) enumeration is driven purely by sequence homology / cut-site compatibility, never by a
user-specified placement rule; (ii) there is no API to say "place motif M at positions p₁…pₙ, in both
orientations, over all subsets"; (iii) the output is a *simulation of possible reaction products*, not a
design specification. The docs send users to the external `teemi` package for combinatorial library work.
**Keep the clause "homology-driven simulation of reaction products, no user-facing placement API" attached
to this cell.** `partial` leans generous (Biopython scores `no`), which is the safe direction.
*Source:* `src/pydna/assembly2.py` (`Assembly` docstring), `docs/notebooks/new_assemblies.ipynb`, `docs/example_gallery.rst`.

**`barcode_generation` — no.** *(carried from v1)*
Zero occurrences of `barcode` across the entire rendered documentation **including the full API reference**
(Sphinx `searchindex.js` audit). No GC- or edit-distance-constrained oligo generation and no attachment
machinery. Named defensively so no author can say it was missed: the closest item is
`primer_screen.expand_iupac_to_dna()` (primer_screen.py L144), which cartesian-expands an extended-IUPAC
oligo into every concrete sequence it encodes (`expand_iupac_to_dna("ATNG")` → `['ATGG','ATAG','ATTG','ATCG']`;
`len(expand_iupac_to_dna("ACGTURYSWKMBDHVN")) == 20736`) — a degenerate pool, but with no GC constraint, no
minimum edit distance and nothing that attaches it to a target. `utils.randomDNA` is unconstrained.
*Source:* `src/pydna/primer_screen.py` L144-L162, `src/pydna/utils.py` L444-L538, docs `searchindex.js`.

**`per_sequence_provenance` — yes.** *(carried from v1)*
The strongest cell in the record. `Dseqrecord.source` is a **typed pydantic model**; subclasses in
`opencloning_models.py` include `PCRSource, SequenceCutSource, RestrictionAndLigationSource,
GibsonAssemblySource, InFusionSource, OverlapExtensionPCRLigationSource, InVivoAssemblySource,
LigationSource, GatewaySource, HomologousRecombinationSource, CRISPRSource, CreLoxRecombinationSource,
RecombinaseSource, OligoHybridizationSource, PolymeraseExtensionSource, AnnotationSource,
ReverseComplementSource, GenomeCoordinatesSource, AddgeneIdSource, BenchlingUrlSource,
SnapGenePlasmidSource, …`, each with operation-specific typed fields. Three verified methods:
`history()` (dseqrecord.py L1521), `validate_history(recursive=True)` (L1554) — "Validate the cloning history
of this sequence by replaying operations" — and `normalize_history()` (L1580). `CloningStrategy` serialises
the whole tree to JSON for the OpenCloning web app. Provenance is also **imported**:
`snapgene_history_parser.parse_snapgene_history()` (L521) reads a SnapGene `.dna` file, re-executes every
recorded step through the corresponding pydna assembly function, and returns a `Dseqrecord` whose `source`
tree mirrors the full cloning provenance — **no other tool in this survey does provenance interop with a
commercial GUI.**
**Caveat that must travel with the cell:** this is provenance of *how one molecule was cloned*, not of *how a
designed variant relates to a wild-type or to a design specification*. It is nevertheless structured, typed,
machine-readable, per-sequence provenance, so `yes` is the honest value.
*Source:* `src/pydna/opencloning_models.py`, `src/pydna/dseqrecord.py` L1521/L1554/L1580, `src/pydna/snapgene_history_parser.py` L1-13 & L521-546, `docs/notebooks/history.ipynb`.

**`automatic_naming` — partial.** *(carried from v1)*
Derived records do get auto-generated names/ids, but they are **operation-tagged, not design-descriptive**,
and are truncated to 16 characters (the GenBank LOCUS limit). From `src/pydna/amplify.py`:

```python
prd.name = (
    identifier_from_string(new_identifier)[:16]
    or self.kwargs.get("name")
    or f"{len(prd)}bp_PCR_prod"[:16]
)
```

and from `dseqrecord.py`: `answer.id = "{name}_lin".format(name=self.name)` with `answer.name =
answer.id[:16]` (linearize, L920-922) and `answer.name = "{}_rc".format(self.name[:13])` (reverse
complement). So `1011bp_PCR_prod`, `pUC19_lin`, `pUC19_rc` appear automatically, and PCR products inherit a
feature label from the template when a feature fully spans the amplicon. No naming scheme encodes a variant's
design (position / mutation / replicate) — because there are no variants. `partial` is already generous
relative to Biopython's `no`.
*Source:* `src/pydna/amplify.py`, `src/pydna/dseqrecord.py` L920-922, `src/pydna/utils.py` L349.

**`design_visualization` — yes.** *(carried from v1)*
Five distinct documented mechanisms: `Dseqrecord.figure()` / `Amplicon.figure()` (ASCII/HTML double-stranded
and primer-annealing figures; `figure(fig_type=...)` renders a compact assembly figure from an `Assembly`
source and a PCR-alignment figure from a `PCRSource`); `Contig.figure()` (contig.py L134, "Compact ascii
representation of the assembled fragments") and `Contig.detailed_figure()` (L80); **`Contig.figure_mpl()`
(L272) — "Graphic representation of the assembly. Returns — matplotlib.figure.Figure"**, drawing fragments on
concentric arcs with per-edge `tab20` colours for circular assemblies, i.e. a genuine graphical assembly view;
`pydna.gel.gel(samples=..., gel_length=600, ...)` rendering a real gel image with a MW ladder; and the
`history()` tree plus `CloningStrategy.model_dump_json()` → loadable in the OpenCloning web GUI. Paper:
homologies "are added to each sequence as metadata in the form of Genbank features (Figure 1A) which can be
inspected graphically using a sequence editor."
**Caveat for the table:** everything above visualises *one construct / one assembly*. Nothing visualises a
library.
*Source:* `src/pydna/contig.py` L80/L134/L272-L300, `src/pydna/dseqrecord.py`, `src/pydna/gel.py`, `README.md`, paper Fig. 1A.

---

## Block B — assay coverage

**`assay_dms` — no.**
No deep-mutational-scanning, saturation-mutagenesis or codon-scanning functionality. `mutagen*` appears on
exactly one docs page in the whole index (the external teemi link). pydna does understand coding sequences —
`seqrecord.isorf()`, `translate()`, `seq.orfs()`/`orfs2()`, `Seq.cai()/rarecodons()/startcodon()/stopcodon()/
express()` (seq.py L142-L172) — but these **diagnose one ORF** and generate **no variants**.
*Source:* `src/pydna/seq.py` L142-L228, `src/pydna/seqrecord.py`, docs `searchindex.js`, paper.

**`assay_mpra` — no.**
No MPRA / STARR-seq / oligo-pool / promoter-enhancer-tiling module, function or docs page, confirmed against
the full docs object index and `example_gallery.rst`. The only "library"-flavoured gallery items are
**external teemi notebooks** (e.g. "04 Promoter library design in *S. cerevisiae*") — a different package,
and metabolic-engineering promoter swapping rather than MPRA.
*Source:* `docs/example_gallery.rst`, docs `searchindex.js`.

**`assay_insilico` — no.**
No machine-learning or model-probing functionality. `pyproject.toml` (129 lines, read in full): required deps
are `biopython ^1.87`, `networkx >=2.8.8`, `numpy`, `prettytable`, `pydivsufsort`, `pyfiglet`, `seguid`,
`regex`, `opencloning-linkml ^1`, `sgffp`, `appdirs`; extras are `pyperclip`, `pyparsing`+`requests`, `cai2`,
`scipy`+`matplotlib`+`pillow`, `pyahocorasick`. **No torch, tensorflow, jax, sklearn or model API.**
*Source:* `pyproject.toml` (full, re-read 2026-08-11).

---

## Block C — genomics integration

**`genome_coordinates` — partial.**
`pydna.genbank.Genbank.nucleotide(item, seq_start=None, seq_stop=None, strand=1)` (genbank.py L69) accepts an
accession-plus-region string; documented forms include `BK006936.2 REGION: complement(613900..615202)`,
`NM_005546 REGION: 1..100`, `21614549:1-100`, `21614549:c100-1`. `BK006936.2` is an *S. cerevisiae*
chromosome record, so chromosome-coordinate retrieval is genuinely supported. **But** only via NCBI
accessions over the network: no genome-build concept, no local genome FASTA/2bit, no `chr7:1000-2000`
interface. `GenomeCoordinatesSource(assembly_accession, locus_tag, gene_id, coordinates)` exists as a
provenance *record type*, but retrieval by assembly accession is implemented in the separate OpenCloning
backend. Corroborating: `sequence_picker.genbank_accession()` (L17) BLASTs a sequence against NCBI `nt` and
returns a `Dseqrecord` described as `f"{accession} REGION: {start}..{stop}"` (L47) — the reverse mapping.
*Consistency:* Biopython also scores `partial` (Entrez `seq_start`/`seq_stop`).
*Source:* `src/pydna/genbank.py` L69-L99, `src/pydna/sequence_picker.py` L17-L47, `src/pydna/opencloning_models.py`.

**`transcript_models` — partial (GFF/GTF specifically: no).**
There is no GFF/GTF parser and nothing imports one (`gff` and `gtf` both return 0 hits across the entire
rendered docs index, API reference included); there is no `Transcript` object and no isoform-selection logic.
`parsers.py` handles GenBank, EMBL, FASTA (`embl_gb_fasta()`), SnapGene and primers. **But annotation via
GenBank feature tables is a transcript model:** INSDC `mRNA`/`CDS` features with `join(...)` locations parse
into `Bio.SeqFeature.CompoundLocation`, and pydna documents constructing, exporting and extracting them in
its own getting-started tutorial `docs/notebooks/Dseq_Features.ipynb` — cell 9 reproduces a real *S. pombe*
EMBL record with `FT   CDS   join(1878362..1878785,1878833..1880604)`; cells 12-13 build a two-exon CDS with
`CompoundLocation` whose executed output contains `CDS   join(4..9,13..15)`; cells 18-19 handle
origin-spanning features on circular molecules.
*Consistency:* `extractions/biopython.md` scores this key `partial (GTF/GFF specifically: no)` for exactly
this mechanism; pydna inherits it from Biopython (`^1.87`, a hard dependency) and demonstrates it in its own
tutorial, so scoring pydna `no` would be self-contradicting.
*Source:* `docs/notebooks/Dseq_Features.ipynb` cells 0/9/12/13/18/19, `src/pydna/parsers.py`, docs `searchindex.js`, `extractions/biopython.md` L88-89.

**`exon_intron_split_codons` — partial.**
`Dseq_Features.ipynb` cell 12 (markdown, verbatim): "To add a feature with parts, like a CDS with introns, we
need to use a `CompoundLocation` object … In a real-world scenario this would represent a CDS with an intron
that skips the `ACG` codon". Cell 13 builds `CompoundLocation([FeatureLocation(3, 9), FeatureLocation(12,
15)])` into a `CDS` `SeqFeature`; executed output shows `location: join{[3:9], [12:15]}` and the GenBank line
`CDS   join(4..9,13..15)`. Cell 15 executes `dummy_record.features[-1].extract(dummy_record).translate()` with
executed output `ProteinSeq('MR')`. This is a first-class API path, not a notebook trick:
`SeqRecord.extract_feature(n)` is `return self.features[n].extract(self)` (seqrecord.py L366), inherited by
`Dseqrecord.extract_feature()` (dseqrecord.py L275); `utils.shift_location()` (L70), `utils.shift_feature()`
(L132) and `utils.location_boundaries()` (L749) all handle `CompoundLocation`. Because exons are concatenated
in transcript order before translation, a codon straddling an exon junction translates correctly.
Why not `yes`: (a) pydna performs no variant design, so it never has to *place* an edit into a
junction-straddling codon — the hard part in a design tool is never exercised; (b) the tutorial's exon lengths
are multiples of three, so it demonstrates intron skipping without exhibiting an actual split codon.
*Consistency:* `extractions/biopython.md` scores `partial` on the identical mechanism.
*Source:* `docs/notebooks/Dseq_Features.ipynb` cells 12-15, `src/pydna/seqrecord.py` L366, `src/pydna/dseqrecord.py` L275, `src/pydna/utils.py` L70/L132/L749, `extractions/biopython.md` L91-92.

**`vcf_vep_output` — no.**
Zero occurrences of `vcf` across the entire docs index. Output formats are GenBank / FASTA / EMBL
(`Dseqrecord.format(format="gb")`, `.write()`), gel PNG images, and the OpenCloning `CloningStrategy` JSON.
**There is no variant representation in the data model at all**, so there is nothing to emit — the absence is
structural, not an omission. No VEP interoperability.
*Source:* docs `searchindex.js`, `src/pydna/seqrecord.py`, `src/pydna/opencloning_models.py`.

**`consequence_annotation` — no.**
No molecular-consequence vocabulary (stop-gained / synonymous / missense / frameshift / in-frame indel)
among the 820 documented objects. Nearest features are `seqrecord.isorf()` and `translate()`, which leave all
interpretation to the user.
*Source:* docs `searchindex.js`, `src/pydna/seqrecord.py`, `src/pydna/seq.py`.

*(v1 `hgvs_input = no` is dropped from the row set per ROWS_v2; the underlying finding — 0 occurrences of
`hgvs` in the whole docs index, no `hgvs` dependency — moves to limitations prose.)*

---

## Block D — adjacent / complementary

**`primer_design` — yes.**
`src/pydna/design.py`: `primer_design()` — "Designs forward and reverse primers for PCR amplification of a
template sequence", able to match a target Tm or an existing primer's Tm; `assembly_fragments()` — "Adds
tails to primers for a linear assembly through homologous recombination or Gibson assembly";
`circular_assembly_fragments()`; `user_assembly_design()` (experimental, USER cloning). Paper: "Pydna
provides the pydna assembly_primers function in order to automatically design tailed primers … The algorithm
tries to create primers with balanced melting temperature for the annealing region." `Amplicon` objects store
the annealing region, per-primer Tm and a suggested PCR program; `src/pydna/tm.py` ships `tm_default`,
`tm_dbd`, `tm_product`, `ta_default`, `ta_dbd`, `program`, `dbd_program`, `Q5`, `tmbresluc`, `tm_neb`.
Additionally `src/pydna/primer_screen.py` performs diagnostic/verification primer-pair **selection across a
set of sequences** (`forward_primers` L290, `reverse_primers` L375, `primer_pairs` L461,
`flanking_primer_pairs` L573, `diff_primer_pairs` L648, `diff_primer_triplets` L760): "Given an iterable of
sequences and a primer list, primers are selected that result in unique product sizes from each of the input
sequences … meant to filter for PCR products that are likely to migrate to sufficiently distinct locations to
be distinguishable on a typical agarose gel."
**The true and only limit, which must travel with the cell:** there is **no mutagenic-primer designer** —
nothing generates QuikChange-style overlapping mutagenic primer pairs from a desired mutation set. pydna can
*simulate* mismatched-primer PCR (`mismatches=N`); it does not *design* the mismatched primer. The v2 row is
worded "Mutagenic primer / oligo design for wet-lab protocols"; `yes` is scored on the oligo/primer-design
limb, which pydna covers unusually well, with the mutagenic gap stated explicitly.
*Source:* `src/pydna/design.py`, `src/pydna/tm.py`, `src/pydna/primer_screen.py` L290-L790, paper (Implementation), `docs/example_gallery.rst`.

**`codon_optimization` — no.**
`src/pydna/codon.py` is **data only** — per-organism codon `weights`, start/stop codon frequency
distributions, `rare_codons` lists, N-end-rule half-life data (yeast `sce` shipped). `utils.cai()` (L225),
`utils.rarecodons()` (L235, returning *positions* as `List[slice]`) and `utils.express()` (L248) are surfaced
as `Seq.cai()`, `Seq.rarecodons()`, `Seq.startcodon()`, `Seq.stopcodon()`, `Seq.express()` (seq.py L146-L172).
**These diagnose codon usage. No function back-translates a protein or rewrites a CDS**; `utils.express` is
documented "Not yet implemented"; `optimi*` returns 0 hits across the docs index.
**Required wording caution:** the same survey scores Biopython `yes (basic)` for
`Bio.SeqUtils.CodonAdaptationIndex.optimize()`, and biopython `^1.87` is a hard dependency of pydna — so that
function is importable in any pydna environment, just not part of pydna's API. Phrase the cell as *"pydna
provides codon-usage diagnostics (CAI, rare-codon localisation, start/stop frequency, N-end rule) but no
codon-optimizing function"* — never "pydna cannot codon-optimize".
*Source:* `src/pydna/codon.py`, `src/pydna/utils.py` L225-L275, `src/pydna/seq.py` L142-L172, docs `searchindex.js`, `extractions/biopython.md` L108-109.

**`synthesis_constraints` — partial.**
Restriction-site content is a **one-call query** on pydna's core objects: `Dseq.no_cutters()` (dseq.py L1856,
"Enzymes in a RestrictionBatch **not** cutting sequence"), `unique_cutters()` (L1866, `n_cutters(n=1)`, with
`once_cutters = unique_cutters` aliased at L1874), `twice_cutters()` (L1876), `n_cutters(n=3, batch=None)`
(L1884), `cutters()` (L1894) — all five defaulting to `CommOnly`, the commercially-available REBASE set — all
six re-exported on `Dseqrecord` (L925-L947, each delegating with `batch or CommOnly`), plus
`Dseqrecord.number_of_cuts(*enzymes)` (L949). Also `Seq.gc()` / `SeqRecord.gc()`, `Seq.rarecodons()`
returning rare-codon *positions*, and the ten-function `pydna.tm` module.
Why not `yes`: **nothing is enforced.** No homopolymer rule, no repeat-content filter, no
hairpin/secondary-structure check, no GC-window scan (only a whole-sequence `Seq.gc()` float), no oligo-length
or vendor/pool constraint, no "reject or repair this sequence" API. `sequence_regex.py` compiles IUPAC
motifs into regexes for recognising recombinase/enzyme sites, not for constraint filtering. **Carry this
sentence with the cell:** *"site-content and Tm queries only; nothing is enforced and there is no repeat /
homopolymer / GC-window / length rule."*
*Consistency:* `extractions/biopython.md` scores `partial` for `Bio.Restriction` + `gc_fraction()` +
`MeltingTemp`; pydna wraps and extends exactly that machinery, so it cannot sit below Biopython here.
*Source:* `src/pydna/dseq.py` L1856-L1900, `src/pydna/dseqrecord.py` L925-L953, `src/pydna/seq.py` L142/L152, `src/pydna/tm.py`, `src/pydna/sequence_regex.py`, `extractions/biopython.md` L111-112.

**`degenerate_iupac_codons` — partial. [NEW ROW]**
Real IUPAC machinery exists, at the nucleotide/oligo level. `primer_screen.expand_iupac_to_dna()`
(primer_screen.py L144) cartesian-expands an extended-IUPAC degenerate oligo into every concrete DNA sequence
it encodes — docstring examples verbatim: `expand_iupac_to_dna("ATNG")` → `['ATGG', 'ATAG', 'ATTG', 'ATCG']`
and `len(expand_iupac_to_dna("ACGTURYSWKMBDHVN")) == 20736`. `sequence_regex.py` compiles IUPAC-degenerate
motifs into regexes (used for recombinase/enzyme site recognition), and `pydna.alphabet` carries the extended
alphabet.
Why `partial` and not `yes`: (i) expansion is **not codon-aware** — there is no NNK/NNS/NNB scheme, no
codon-table-driven degenerate design, and no restriction of the expansion to an amino-acid set; (ii) there is
no **compression** (nothing collapses a set of desired codons into a minimal degenerate codon); (iii) the
function lives in the *primer-screening* module and exists to enumerate a degenerate primer stock for
Aho-Corasick annealing search, not to specify a library — the word `degenerate` returns **0 hits** across the
whole rendered docs index including the API reference, so this is a raw primitive the user would assemble
into a design themselves. Under the "raw primitives → partial at most" rule, `partial`.
*Source:* `src/pydna/primer_screen.py` L144-L162, `src/pydna/sequence_regex.py`, `src/pydna/alphabet.py`, docs `searchindex.js` (`degenerate` 0 hits), `reviews/pydna.md` (keyword audit).

**`negative_control_generation` — no. [NEW ROW]**
No scramble, shuffle, dinucleotide-preserving shuffle, or control-set generator anywhere in the 820
documented objects. Named defensively: `utils.randomDNA/randomRNA/randomORF/randomprot` (utils.py L444-L538)
produce *de novo* random sequences from a length argument — they are not composition- or
dinucleotide-matched to an input and are not framed or documented as controls; and
`Dseqrecord.reverse_complement()` / `ReverseComplementSource` exist as a molecular-biology operation on a
construct (with its own provenance record), not as a control-generation feature attached to a design. Nothing
generates a matched control *set* for a designed sequence.
*Source:* `src/pydna/utils.py` L444-L538, `src/pydna/dseqrecord.py`, `src/pydna/opencloning_models.py`, docs `searchindex.js`.

**`ml_model_in_loop` — no. [NEW ROW]**
No predictive model participates anywhere. `pyproject.toml` (read in full) declares no torch, tensorflow, jax,
sklearn or model API among required deps (`biopython`, `networkx`, `numpy`, `prettytable`, `pydivsufsort`,
`pyfiglet`, `seguid`, `regex`, `opencloning-linkml`, `sgffp`, `appdirs`) or extras (`pyperclip`,
`pyparsing`+`requests`, `cai2`, `scipy`+`matplotlib`+`pillow`, `pyahocorasick`). No scoring, ranking,
optimisation-against-a-predictor, or model-guided edit selection exists. (The `tm_neb` function queries the
NEB web API for a melting temperature — a thermodynamic calculator, not a predictive model in a design loop.)
*Source:* `pyproject.toml` (full, 129 lines), docs `searchindex.js`.

**`readout_analysis` — no. [NEW ROW]**
Nothing analyses sequencing readout from a built library: no FASTQ/BAM reading, no demultiplexing, no barcode
counting, no variant-frequency or enrichment computation. The closest items are wet-lab *verification design
and simulation*, not readout analysis: `pydna.gel.gel()` renders a simulated gel image;
`primer_screen.diff_primer_pairs()` selects diagnostic primer pairs giving gel-distinguishable product sizes;
`Dseqrecord.validate_history()` replays a cloning history against its original inputs; and
`sequence_picker.genbank_accession()` BLASTs a single sequence against NCBI `nt`. Since there is no library
representation, there is also no readout to map back to one.
*Source:* `src/pydna/gel.py`, `src/pydna/primer_screen.py` L648-L790, `src/pydna/dseqrecord.py` L1554, `src/pydna/sequence_picker.py`, `pyproject.toml`, docs `searchindex.js`.

---

## Block E — engineering and availability

**`interface` — yes → Python API / importable library only (no CLI, no GUI, no web service).**
`import pydna`, intended for scripts and Jupyter notebooks. `pyproject.toml` (129 lines, full read) declares
**no `[project.scripts]` and no console entry points**, so `pip install pydna` installs no command-line tool;
the classifier `"Environment :: Console"` refers to REPL use. The paper's "pydna live" web console
(`http://pydna-shell.appspot.com/`) is **dead (HTTP 404)**. The paper's phrase "designed purely as a command
line tool" means *scripting / non-GUI*, not an installed executable — gloss it that way if quoted.
**OpenCloning** (opencloning.org) is a **separate project** providing a web GUI over pydna's `CloningStrategy`
JSON and must not be credited to this row.
*Source:* `pyproject.toml` (full), `README.md`, paper (Comparison with existing tools; Availability), live HTTP checks.

**`license` — yes → BSD 3-Clause (OSI-approved, permissive).**
`LICENSE.txt`: "Copyright (c) 2013-2026 Björn Johansson, Centro de Biologia Molecular e Ambiental (CBMA)/
Department of Biology, University of Minho, Braga, Portugal" plus the three standard BSD clauses. Every
source file carries `SPDX-License-Identifier: BSD-3-Clause`; `pyproject.toml` has `license-files =
["LICENSE.txt"]` and `[tool.poetry] license = "BSD"`; PyPI classifier "License :: OSI Approved :: BSD
License". Footnotes: the 2015 paper says "License: FreeBSD" (2-clause) whereas the shipped text is 3-clause,
and GitHub's detector reports `NOASSERTION` because the file is `LICENSE.txt` with a custom preamble —
neither affects the value.
*Source:* `LICENSE.txt`, `pyproject.toml`, PyPI JSON classifiers, paper (Availability).

**`installable_today` — yes → `pip install pydna` → 5.5.16 (uploaded 2026-06-09), pure Python.**
`requires_python = ">=3.10,<4.0"`; classifiers cover Python 3.10, 3.11, 3.12, 3.13, 3.14. Optional extras:
`clipboard`, `download`, `express`, `gel`, `primer_screen`. No compiled components; no manual steps. Release
cadence is tight: 5.5.9 (2026-03-31) → 5.5.10 (04-10) → 5.5.11 (04-20) → 5.5.12 (06-02) → 5.5.13 (06-03) →
5.5.14 (06-04) → 5.5.15 (06-05) → 5.5.16 (06-09). Docs live at https://pydna-group.github.io/pydna (HTTP 200).
*Source:* PyPI JSON (re-verified 2026-08-11), `pyproject.toml`, docs site HTTP check.

**`last_activity` — yes → `master` head commit 2026-08-01; repo `pushed_at` 2026-08-10; PyPI 5.5.16 2026-06-09.**
Verified without the GitHub REST API (rate-limited) via the commit Atom feed
`https://github.com/pydna-group/pydna/commits/master.atom`: head **2026-08-01T09:24:09Z**, with further
commits 2026-07-16, 07-15, 07-12. Repo not archived, 228 stars, 67 open issues, CI green. A **v6 major is in
development** (`master` version string `6.0.0-a.24.post.17+b7b559bd66`). **Two active maintainers** (Björn F.
Johansson, Manuel Lera-Ramirez). Last activity is **10 days before the survey date** — pydna is among the most
actively maintained tools in this survey.
*Source:* commits Atom feed, PyPI JSON, `pyproject.toml`, GitHub repo page — all 2026-08-11.

**`documented_examples` — yes → ~19 first-party examples: 12 tutorial notebooks + 4 gallery examples + 3 paper cookbook cases.**
Tutorial notebooks (`docs/notebooks/`, rendered on the docs site): `Dseq`, `Dseq_Features`, `Importing_Seqs`,
`Restrict_Ligate_Cloning`, `PCR`, `primer_design`, `Gibson`, `CRISPR`, `history`, `new_assemblies`,
`readme_example`, `hackathon_copenhagen_2025` (12). Example gallery (`docs/example_gallery.rst`, pydna's own,
no external package required): `Example_Restriction`, `Example_Gibson`, `Example_CRISPR`,
`Example_primer_screen` (4). Paper cookbook, still shipped in `docs/cookbook/`: YEp24PGK_XK (12 lines, 11,452
bp product), pGUP1 (11 lines, 9,981 bp), and the Lactose pathway ("nineteen short pydna scripts",
experimentally validated in *S. cerevisiae*) (3). Plus a full 820-object API reference.
**Not counted:** the gallery also links out to third-party `teemi` notebooks (HT CAZyme primer design, Kozak
library, primer-directed mutagenesis, promoter library) — these are a different package and are the only
"library"-flavoured items anywhere near pydna; a referee may cite them, so name them and exclude them
explicitly.
*Source:* `docs/example_gallery.rst`, `docs/notebooks/`, `docs/cookbook/`, paper, docs `searchindex.js`.

---

## v2 summary table

| Key | Value | One-line basis |
|---|---|---|
| library_first_class_object | **no** | No collection/pool class among 820 documented objects; the first-class object is one molecule (`Dseqrecord`), multi-product results are plain `list` |
| composable_operations | **yes** | No fixed pipeline: every operation returns a `Dseqrecord` any other accepts; arbitrary chaining/nesting across ~15 assembly back-ends, recorded in a typed `source` DAG |
| lazy_generation | **no** | `parse()` "a greedy function"; all products materialised; internal generators are explosion guards |
| library_algebra | **no** | No stack/concat/sample/repeat over libraries; pooling is user Python outside the tool |
| exhaustive_single_scans | **no** | No scanning operation; can only *simulate* one user-specified edit (`mismatches=N` PCR, CRISPR+oligo) |
| sampled_random_mutagenesis | **no** | `randomDNA` makes de-novo random sequences from a length, not rate-based variants of a parent |
| higher_order_combinatorial | **no** | `Assembly` enumerates fragment combinations, never mutation combinations |
| heterogeneous_components_one_library | **no** | No library specification exists; a heterogeneous fragment list is the parts list for one construct |
| combinatorial_motif_place | **partial** | `Assembly` graph enumerates all fragment orders/orientations; homology-driven, no user placement API |
| barcode_generation | **no** | `barcode` 0 hits in whole docs index; `expand_iupac_to_dna()` has no GC/edit-distance constraint and no attachment |
| per_sequence_provenance | **yes** | Typed pydantic `source` tree, `history()`, `validate_history()`, `normalize_history()`, SnapGene history *import*; cloning provenance, not design-vs-WT |
| automatic_naming | **partial** | `1011bp_PCR_prod`, `_lin`, `_rc`; 16-char truncated, operation-tagged not design-descriptive |
| design_visualization | **yes** | `figure()`, `detailed_figure()`, `figure_mpl()` → matplotlib Figure, `gel.gel()` image, `history()` tree, OpenCloning JSON — all per construct |
| assay_dms | **no** | No DMS/saturation/codon scanning; `mutagen*` on one docs page (external link) |
| assay_mpra | **no** | No MPRA/STARR/oligo-pool/tiling module or docs page |
| assay_insilico | **no** | No ML dependency or model-probing API in `pyproject.toml` |
| genome_coordinates | **partial** | `Genbank.nucleotide("BK006936.2 REGION: complement(613900..615202)")`; NCBI accessions only, no build/local genome |
| transcript_models | **partial** *(GFF/GTF: no)* | GenBank/EMBL `join()` → `CompoundLocation`, documented in own tutorial; no GFF/GTF parser, no transcript object |
| exon_intron_split_codons | **partial** | Tutorial builds spliced CDS, exports `join(4..9,13..15)`, extracts and translates to `MR`; never used for edit placement |
| vcf_vep_output | **no** | 0 occurrences of `vcf`; no variant representation exists to emit |
| consequence_annotation | **no** | No consequence vocabulary among 820 documented objects |
| primer_design | **yes** | `primer_design`, `assembly_fragments`, `user_assembly_design`, 10-function `tm` module, diagnostic primer-pair selection across sequence sets; **no mutagenic-primer designer** |
| codon_optimization | **no** | CAI / rare-codon diagnostics only; `codon.py` is data; `express` "Not yet implemented"; `optimi*` 0 hits |
| synthesis_constraints | **partial** | `no_cutters/unique_cutters/once_cutters/twice_cutters/n_cutters/cutters/number_of_cuts` over `CommOnly` + `gc()` + `tm`; nothing enforced |
| degenerate_iupac_codons | **partial** | `expand_iupac_to_dna()` expands IUPAC oligos (20,736 from `ACGTURYSWKMBDHVN`); nucleotide-level, no codon scheme, no compression, `degenerate` 0 hits in docs |
| negative_control_generation | **no** | No scramble/shuffle/matched-control generator; `randomDNA` is de-novo and unconstrained; `rc()` is a molecular operation |
| ml_model_in_loop | **no** | No torch/tf/jax/sklearn/model API in required deps or extras |
| readout_analysis | **no** | No FASTQ/counting/enrichment; gel simulation and diagnostic-PCR selection are verification *design*, not readout analysis |
| interface | **yes** | Python library API; no `[project.scripts]`, no CLI, no GUI, no web service (OpenCloning is a separate project) |
| license | **yes** | BSD 3-Clause, SPDX headers in every file, OSI classifier |
| installable_today | **yes** | `pip install pydna` → 5.5.16 (2026-06-09), pure Python, Python 3.10-3.14, 5 extras |
| last_activity | **yes** | `master` commit 2026-08-01; push 2026-08-10; PyPI 5.5.16 2026-06-09; v6 alpha in dev; two maintainers |
| documented_examples | **yes** | ~19 first-party: 12 tutorial notebooks + 4 gallery examples + 3 paper cookbook cases, plus 820-object API reference |

---

## Changes from the v1 record

| Row | v1 | v2 | Why |
|---|---|---|---|
| `dag_chaining` → `composable_operations` | partial | **yes** | Rename that genuinely changes the question. v1's `partial` rested on "the DAG is a record, not a re-executable program" — that criterion is not what v2 asks. v2 asks composition vs fixed pipeline; pydna has no pipeline and composes arbitrarily. The v1 caveats are now carried by `lazy_generation`, `library_algebra` and `library_first_class_object`. |
| `lazy_evaluation` → `lazy_generation` | no | no | Pure rename; value and evidence carried unchanged. |
| `library_as_object` → `library_first_class_object` | no | no | Split; pydna was already `no` on the combined row. |
| `library_as_object` → `library_algebra` (new) | — | **no** | No library object to combine; pooling/sampling is user code outside the tool. |
| `mixed_mutagenesis_one_pool` → 4 rows | no | **no / no / no / no** | All four resolve `no` independently. pydna simulates a single user-specified edit; it designs none, samples none, combines none, and has no pool. |
| `hgvs_input` | no | *dropped* | Removed from the row set; finding moves to limitations prose. |
| `degenerate_iupac_codons` (new) | — | **partial** | `expand_iupac_to_dna()` is genuine IUPAC expansion, but nucleotide-level, not codon-aware, no compression, undocumented as a design feature. |
| `negative_control_generation` (new) | — | **no** | Nothing framed or shipped as control generation. |
| `ml_model_in_loop` (new) | — | **no** | No ML dependency of any kind. |
| `readout_analysis` (new) | — | **no** | No sequencing-readout analysis; gel/diagnostic-PCR are verification design. |
| `installable_today` (new) | — | **yes** | From `availability_status`: pip, 5.5.16, pure Python. |
| `last_activity` (new) | — | **yes** | From `availability_status`: commit 2026-08-01, push 2026-08-10. |
| `documented_examples` (new) | — | **yes** | Counted from the memo's notebook/gallery/cookbook inventory. |
| `maintained` | yes | *folded into* `last_activity` | v2 has no `maintained` row; the evidence is carried there. |

## Cells flagged for human review

- `composable_operations = yes` — the one upgrade in this record; a defensible reading of the new wording, but
  a reviewer may prefer to hold it at `partial` for continuity with v1. Whichever is chosen, the same rule must
  be applied to every tool whose v1 `dag_chaining` was discounted for the DAG/declarative criterion rather than
  for composability.
- `library_first_class_object = no` — the counterargument (individual molecules *are* first-class,
  transformable, passable objects; the tool never merely writes files) is real and is stated in the cell. `no`
  is scored because the row asks about a *library* object.
- `combinatorial_motif_place = partial` — a v1 judgement call in the generous direction; the enumeration is
  over fragments by homology, not motifs by rule.
- `degenerate_iupac_codons = partial` — new row; hinges on whether nucleotide-level IUPAC expansion inside a
  primer-screening module counts. Could be argued to `no` (not codon-aware, `degenerate` 0 hits in docs).
- `primer_design = yes` — the v2 row wording leads with "mutagenic primer", and pydna has no mutagenic-primer
  designer. `yes` is scored on the oligo/primer-design limb.
- `synthesis_constraints = partial`, `transcript_models = partial`, `exon_intron_split_codons = partial` — all
  three are coupled to Biopython's scores in this survey (pydna hard-depends on biopython `^1.87` and wraps the
  same machinery). If any is hardened to `no`, Biopython must be hardened identically in the same pass.
- `codon_optimization = no` — a referee could argue `partial` from the CAI / rare-codon diagnostics; the cell
  must read "no codon-optimizing function", never "cannot codon-optimize".
- `per_sequence_provenance = yes` — value is right, but the caveat (cloning provenance, not design-vs-WT
  provenance) is load-bearing and must be printed.
- `documented_examples` — the count (~19) depends on excluding the linked third-party `teemi` notebooks; if
  other tools' counts include externally linked material, the rule must be applied consistently.
- `negative_control_generation = no` and `readout_analysis = no` — resolved from the memo's
  `additional_capabilities` / limitations sections rather than from a fresh targeted search of the repo; both
  are consistent with the whole-docs `searchindex.js` audit, but neither keyword was audited by name.
