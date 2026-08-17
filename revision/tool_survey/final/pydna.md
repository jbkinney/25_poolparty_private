# pydna — FINAL capability record

**Slug:** `pydna`
**Full name:** pydna — "a simulation and documentation tool for DNA assembly strategies using python"
**Citation key:** `Pereira2015wj`
**Tier:** 1
**Paper:** Pereira F, Azevedo F, Carvalho Â, Ribeiro GF, Budde MW, Johansson B. *Pydna: a simulation and
documentation tool for DNA assembly strategies using python.* **BMC Bioinformatics** 16:142 (2015).
doi:10.1186/s12859-015-0544-x
**Authors / maintainers:** Björn F. Johansson (Univ. Minho, CBMA), Manuel Lera-Ramirez
**Version assessed:** PyPI 5.5.16 (uploaded 2026-06-09); `master` at in-development `6.0.0-a.24.post`
(latest commit 2026-08-01T09:24:09Z)
**Record status:** FINAL. Adversarial review returned **21 of 24 keys "supported"** and **3 keys
"understated"**. All three understated cells were changed (`transcript_models`,
`exon_intron_split_codons`, `synthesis_constraints`: `no` → `partial`), three evidence sentences were
softened without changing their values (`mixed_mutagenesis_one_pool`, `primer_design`,
`codon_optimization`), and all ten missed capabilities were incorporated. **Every changed cell and
every missed capability was independently re-verified from `master` source and the rendered docs for
this memo** (see "Re-verification log"). No unresolved extractor/reviewer disagreements remain — no
cell is set to `unknown`.

---

## Sources consulted

| Kind | Reference | Retrieved |
|---|---|---|
| pdf | `revision/tool_survey/papers/Pereira2015_pydna.pdf` — BMC Bioinformatics 16:142 (2015) | 2026-08-10 |
| repo | https://github.com/pydna-group/pydna (the paper's `github.com/BjornFJohansson/pydna` 301-redirects here) | 2026-08-11 |
| repo (raw, `master`) | `README.md`, `LICENSE.txt`, `pyproject.toml` (full, 129 lines) | 2026-08-11 (re-verified) |
| repo (raw, `master`) | `src/pydna/dseq.py` (3,173 L), `dseqrecord.py` (1,636 L), `assembly2.py` (3,289 L), `contig.py` (462 L), `primer_screen.py` (917 L), `snapgene_history_parser.py` (570 L), `seqrecord.py` (731 L), `utils.py` (1,021 L), `seq.py` (340 L), `tm.py` (421 L), `genbank.py`, `sequence_picker.py` | 2026-08-11 (re-verified) |
| repo (raw, `master`) | `src/pydna/{design,codon,parsers,amplify,opencloning_models,gel,crispr,sequence_regex,fusionpcr,oligonucleotide_hybridization}.py` | 2026-08-10 |
| repo (raw, `master`) | `docs/notebooks/Dseq_Features.ipynb` — **downloaded and parsed cell-by-cell with executed outputs** for this memo | 2026-08-11 |
| docs | https://pydna-group.github.io/pydna (HTTP 200); `docs/example_gallery.rst`; Sphinx `searchindex.js` (258 kB — stems every word of every docs page *including the full API reference*; 820 documented objects enumerated) | 2026-08-10 / 2026-08-11 |
| docs (legacy) | https://pydna.readthedocs.io/latest/ (still HTTP 200; not the maintained site) | 2026-08-10 |
| pypi | https://pypi.org/pypi/pydna/json — version, upload timestamps, classifiers, `requires_python` | 2026-08-14 (re-verified) |
| repo metadata | GitHub repository and `master`-commit APIs | 2026-08-14 (re-verified) |

No install, no clone, no execution — document/web research only, per the standing constraint.

---

## What pydna actually is

From the abstract (verbatim):

> "We have developed pydna, an extensible, free and open source Python library for simulating basic
> molecular biology DNA unit operations such as restriction digestion, ligation, PCR, primer design,
> Gibson assembly and homologous recombination. A cloning strategy expressed as a pydna script provides
> a description that is complete, unambiguous and stable. Execution of the script automatically yields
> the sequence of the final molecule(s) and that of any intermediate constructs."

> "Most of the pydna functionality is implemented as methods for the Dseqrecord class, which was
> designed to hold sequence information necessary for describing a double-stranded DNA molecule."
> (Implementation)

> "Pydna simplifies both the planning and sharing of cloning strategies and is especially useful for
> complex or combinatorial DNA molecule construction." (Abstract)

**Read "combinatorial" carefully.** In pydna it means *multi-fragment assemblies* — many parts, one
intended product, with all mechanistically possible products enumerated to catch strategy errors. It
does **not** mean variant libraries. Throughout this record, *variant library* means a designed set of
sequence variants; reaction-product lists and the sequences collected in one `CloningStrategy` are named
explicitly rather than treated as variant libraries.

**The package has grown enormously since the 2015 paper**, which must not be quoted as the current API.
Added since: Golden Gate, Gateway (BP/LR), Cre-lox, generic
recombinases, CRISPR cut + HR repair, fusion PCR, oligo hybridization, a typed provenance / "cloning
history" data model, gel-image rendering, a new `assembly2` engine, SnapGene-history import, a
matplotlib assembly view, Aho-Corasick primer screening, and a second maintainer. Pydna now stores typed,
machine-readable per-sequence provenance (see `per_sequence_provenance`), but still has no variant-library
abstraction.

Modules in `src/pydna/`:
`all, alphabet, amplicon, amplify, assembly, assembly2, codon, common_sub_strings, contig, cre_lox,
crispr, design, dseq, dseqrecord, fakeseq, fusionpcr, gateway, gel, genbank, genbankfixer, ladders,
oligonucleotide_hybridization, opencloning_models, parsers, primer, primer_screen, readers, recombinase,
seq, seqrecord, sequence_picker, sequence_regex, snapgene_history_parser,
threading_timer_decorator_exit, tm, types, utils`

Assembly back-ends in `assembly2.py`: `gibson_assembly, in_fusion_assembly, fusion_pcr_assembly,
in_vivo_assembly, restriction_ligation_assembly, golden_gate_assembly, ligation_assembly,
gateway_assembly, homologous_recombination_integration/excision_or_inversion,
cre_lox_integration/excision_or_inversion, recombinase_integration/excision_or_inversion,
recombinase_assembly, crispr_integration` (the plain `*_excision` names are deprecated wrappers), plus
`PCRAssembly` and `SingleFragmentAssembly`.

---

## Capability assessment

Legend: **value** ∈ {yes, partial, no, unknown}. Cells whose value or evidence changed from the original
extraction are marked **[CHANGED]** or **[EVIDENCE SOFTENED]**.

### Block A — library specification

**`library_as_object` — no.**

> **Superseded — do not cite this value.** This entry scores the *original* row wording. The row was
> subsequently reworded (`ROWS_v3.md`) to drop the "passed to another operation" and "size without
> materialising" criteria, which are now separate rows. Under the original three-criterion rubric the
> verified audit scored pydna **partial** (criteria (i) and (iii) pass via `Assembly` and
> `assembly2.get_possible_assembly_number`; only (ii) fails). Under the adopted reworded rubric it is
> **yes** — `Assembly` is a purpose-built type, exported in `pydna.all.__all__`. The shipping value is
> in `verified/library_first_class_object.md`; the paragraph below is retained as the original
> reasoning only. Its absolute claim about the 820-object index was withdrawn by the citation audit
> (the index does contain `RecombinaseCollection`).

There is no **variant-library-specific** class. `CloningStrategy.from_dseqrecords()` accepts terminal
`Dseqrecord`s and collects one strategy's sequences, sources and primers; `to_dseqrecords()` reconstructs
the source trees and returns that strategy's terminal records. It is a purpose-built strategy/provenance
collection, not a designed variant set. The other documented collection-named types are
`recombinase.RecombinaseCollection` (recombinase enzymes) and
`opencloning_models.OpenDNACollectionsSource` (a repository-ID provenance record). The molecular unit passed
through operations is a single `Dseqrecord` (one dsDNA molecule).
Operations that can yield more than one molecule return **plain Python containers**:
`Assembly.assemble_linear(...) -> list[Dseqrecord]`, `assemble_circular(...) -> list[Dseqrecord]`,
`Dseqrecord.cut(*enzymes) -> Tuple[Dseqrecord, ...]`. Those multiple elements are *alternative products of
one reaction*, not a designed set. Building many variants means writing your own loop. The paper's stated
goal is "the sequence of the final molecule(s)".
*Source:* `src/pydna/assembly2.py`, `src/pydna/dseqrecord.py`, `src/pydna/opencloning_models.py`, docs `searchindex.js` object index, paper abstract.

**`dag_chaining` — partial.**
Steps compose naturally (each operation returns a `Dseqrecord` that feeds the next:
`pcr(...)` → `.cut(...)` → `.looped()` → `gibson_assembly([...])`), and pydna
**records** the resulting dependency graph — already in the current stable release, PyPI 5.5.16, not only on
the v6 development line: a `Dseqrecord` produced by a modelled operation carries a typed `source` object
pointing at its inputs (pydna's own caveat, "Not all methods generate sources so far", is real —
`Dseqrecord.__add__`, for instance, assigns none), and `Dseqrecord.history()` (dseqrecord.py L1521) prints the
tree. Verbatim from the docstring:

```
>>> product, *_ = gibson_assembly(fragments, limit=4)
>>> product.name = "product_name"
>>> print(product.history())
╙── product_name (Dseqrecord(o34))
    └─╼ GibsonAssemblySource
        ├─╼ fragment1 (Dseqrecord(-21))
        ├─╼ fragment2 (Dseqrecord(-12))
        └─╼ fragment3 (Dseqrecord(-13))
```

Why not `yes`: `CloningStrategy` does round-trip —
`from_dseqrecords()`, `add_dseqrecord()`, `add_primer()`, `reassign_ids()`, `model_dump()/model_dump_json()`
on the way out, and `to_dseqrecords()`, `normalize()`, `validate()` rebuilding the `Dseqrecord`s and their
source trees on the way back in — but each of those reads the sequences embedded in the strategy itself.
There is **no re-execution of a serialized strategy on new inputs inside
pydna**. `validate_history(recursive=True)` replays
operations *only against the original inputs*, for verification. So the DAG is a record of eagerly executed
steps on one molecule lineage, not a declarative pipeline that is built first, evaluated later, and
parameterised to fan out over a variant set.
Note the 2015 paper's Figure 5 dependency graph ("A dependency graph produced from the Lactose pathway pydna
source files") was produced by an *external software-dependency tool*, not by pydna.
*Source:* `src/pydna/dseqrecord.py` L1521/L1554/L1580, `src/pydna/opencloning_models.py`, paper Fig. 5.

**`lazy_evaluation` — no.**
Execution is eager: "Execution of the script automatically yields the sequence of the final molecule(s) and
that of any intermediate constructs" (abstract). `parsers.parse()` is annotated
`-> list[Dseqrecord | SeqRecord]` and documented as "a greedy function, use carefully";
`Assembly.assemble_linear/assemble_circular/assemble_insertion` return `list[Dseqrecord]` — all products
materialised. Internal generators (`utils.limit_iterator` at utils.py L920, wrapping
`nx.cycles.simple_cycles(self.G)` capped at 10,000 in `assembly2.py`) are **combinatorial-explosion guards on
graph path enumeration**, not user-facing lazy sequence generation.
*Source:* `src/pydna/parsers.py`, `src/pydna/assembly2.py`, `src/pydna/utils.py` L920, paper abstract.

**`mixed_mutagenesis_one_pool` — no. [EVIDENCE SOFTENED]**
Value unchanged and confirmed: no mutagenesis *design* operation — no substitution/deletion/insertion
scanning, no exhaustive singles or pairs, no random variant sampling, no WT-replicate concept, and no pool
into which different mutagenesis types could be mixed. `mutagen*` occurs on exactly **one** page of the
entire rendered documentation (`example_gallery`, and there only as an outbound link to the third-party
`teemi` repo's "02 Primer directed mutagenesis with pydna" notebook).

**The original evidence sentence "pydna has no mutagenesis operation at all" was overreaching and is
withdrawn.** pydna *can simulate* an individual designed edit, two ways, both re-verified in source:
- `assembly2.pcr_assembly(template, fwd_primer, rvs_primer, add_primer_features=False, limit=14,
  mismatches=0)` (assembly2.py L3244) and `PCRAssembly(frags, limit=25, mismatches=0)` (L1951) accept
  substitution mismatches between primer and template. `primer_template_overlap()` (L697) documents
  `mismatches : int — "Maximum number of mismatches (only substitutions, no deletion or insertion)"` and
  implements it as a fuzzy regex, verbatim `"(" + query + "){s<=" + str(mismatches) + "}"` (L765/L770).
  That is simulation of mutagenic-primer (QuikChange-style) PCR: the primer's mismatched bases end up in
  the product.
- `assembly2.crispr_integration` / `homologous_recombination_integration`, plus the documented
  `Example_CRISPR` notebook — "Using CRISPR with homologous recombination to delete genes by making two cuts
  in the genome, and repair it with an oligo."

**Correct wording for the table:** *"pydna can simulate an individual edit (mismatched-primer PCR, CRISPR +
oligo repair) but has no mutagenesis design operation: no scanning, no exhaustive singles/pairs, no random
sampling, no WT replicates, and no pool to mix them into."*
*Source:* `src/pydna/assembly2.py` L697-L790, L1947-L1988, L3244-L3280; `docs/example_gallery.rst`; docs `searchindex.js`.

**`combinatorial_motif_place` — partial.**
`Assembly` builds a directed graph over fragments in which "The sign of the node key represents the
orientation of the fragment, positive for forward orientation, negative for reverse orientation"
(`Assembly` docstring), then enumerates linear / circular / insertion paths through it. Circular enumeration
ignores fragment order and so is exhaustive over orders and orientations; linear enumeration defaults to
`use_fragment_order=True` (paths must start with the first fragment and end with the last), every getter
defaults to `max_assemblies=50`, and cycle enumeration is capped at 10,000 paths. From the
`new_assemblies` notebook: "there are two possible circular assemblies, one with all fragments oriented as
they were passed to the class constructor, and one with the second fragment inverted." So a set of parts
sharing homology genuinely is combined in **all valid orders and orientations** (circularly, or linearly with
`use_fragment_order=False`) — a real combinatorial enumeration.
Why not `yes`: (i) enumeration is driven purely by sequence homology / cut-site compatibility, never by a
user-specified placement rule; (ii) there is no API to say "place motif M at positions p₁…pₙ, in both
orientations, over all subsets"; (iii) the output is a *simulation of possible reaction products* (a
strategy-error check), not a design specification. The official gallery places its combinatorial-library
workflow under "External examples" and says that example uses the extra `teemi` dependency.
**Keep the clause "homology-driven simulation of reaction products, no user-facing placement API" attached to
this cell.**
*Source:* `src/pydna/assembly2.py` (`Assembly` docstring), `docs/notebooks/new_assemblies.ipynb`, `docs/example_gallery.rst`.

**`barcode_generation` — no.**
Zero occurrences of `barcode` across the **entire rendered documentation including the full API reference**
(verified in the Sphinx `searchindex.js` object/text index). No GC- or edit-distance-constrained oligo
generation exists.
`utils.randomDNA/randomRNA/randomORF/randomprot` (utils.py L444-L538) take only `length` (and optional
`maxlength`) — unconstrained.

The closest documented function is
`primer_screen.expand_iupac_to_dna()` (primer_screen.py L144), which cartesian-expands an extended-IUPAC
degenerate oligo into every concrete DNA sequence it encodes — verbatim docstring examples
`expand_iupac_to_dna("ATNG")` → `['ATGG', 'ATAG', 'ATTG', 'ATCG']`, and
`len(expand_iupac_to_dna("ACGTURYSWKMBDHVN")) == 20736`. That enumerates a degenerate pool but applies
**no GC constraint, no minimum edit distance, and has no attachment machinery**, so `no` still holds.
*Source:* `src/pydna/primer_screen.py` L144-L162, `src/pydna/utils.py` L444-L538, docs `searchindex.js`.

**`per_sequence_provenance` — yes.**
The supporting evidence is more specific than the original extraction claimed. From
`src/pydna/opencloning_models.py` module docstring (verbatim):

> "When using pydna to plan cloning, it stores the provenance of ``Dseqrecord`` objects in their ``source``
> attribute. Not all methods generate sources so far, so refer to the documentation notebooks for examples on
> how to use this feature. The ``history`` method of ``Dseqrecord`` objects can be used to get a string
> representation of the provenance of the sequence. You can also use the ``CloningStrategy`` class to create a
> JSON representation of the cloning strategy. That ``CloningStrategy`` can be loaded in the OpenCloning web
> interface to see a representation of the cloning strategy."

`Dseqrecord.source` is a **typed pydantic model**. Source subclasses in `opencloning_models.py`:
`PCRSource, SequenceCutSource, RestrictionAndLigationSource, GibsonAssemblySource, InFusionSource,
OverlapExtensionPCRLigationSource, InVivoAssemblySource, LigationSource, GatewaySource,
HomologousRecombinationSource, CRISPRSource, CreLoxRecombinationSource, RecombinaseSource,
OligoHybridizationSource, PolymeraseExtensionSource, AnnotationSource, ReverseComplementSource,
GenomeCoordinatesSource, AddgeneIdSource, BenchlingUrlSource, SnapGenePlasmidSource, …` — each with
operation-specific typed fields (e.g. `RestrictionAndLigationSource.restriction_enzymes: list[AbstractCut]`).
Three verified methods on `Dseqrecord`: `history()` (L1521), `validate_history(recursive=True)` (L1554)
— "Validate the cloning history of this sequence by replaying operations" — and
**`normalize_history()` (L1580)**, which the original extraction missed.

**Provenance is also IMPORTED, not only exported** (missed by the original extraction, re-verified here).
`src/pydna/snapgene_history_parser.py` module docstring, verbatim:

> "Parse SnapGene .dna files and reconstruct their cloning history as `Dseqrecord` objects. The single public
> entry point is `parse_snapgene_history`, which reads a ``.dna`` file, resolves every step recorded in its
> SnapGene history, and returns a ``Dseqrecord`` whose ``source`` attribute tree mirrors the full cloning
> provenance."

`parse_snapgene_history(data, file_name="")` (L521) re-executes each recorded SnapGene step through
`gibson_assembly`, `pcr_assembly`, `restriction_ligation_assembly`, `gateway_assembly`,
`in_fusion_assembly`, `ligation_assembly`, `fusion_pcr_assembly` and `oligonucleotide_hybridization`,
raising `ValueError` "If a recorded operation cannot be reproduced (no matching product found)". Depends on
the `sgffp` package.

**Caveat that must travel with the cell:** this is provenance of *how one molecule was cloned*, not of *how a
designed variant relates to a wild-type / design specification*. It is nevertheless structured, typed,
machine-readable, per-sequence provenance, so `yes` follows the stated definition.
*Source:* `src/pydna/opencloning_models.py`, `src/pydna/dseqrecord.py` L1521/L1554/L1580, `src/pydna/snapgene_history_parser.py` L1-13 & L521-546, `docs/notebooks/history.ipynb`.

**`automatic_naming` — partial.**
Derived records do get auto-generated names/ids, but they are **operation-tagged, not design-descriptive**,
and are truncated to 16 characters (the GenBank LOCUS limit). Verbatim from `src/pydna/amplify.py`:

```python
prd.name = (
    identifier_from_string(new_identifier)[:16]
    or self.kwargs.get("name")
    or f"{len(prd)}bp_PCR_prod"[:16]
)
```

and from `src/pydna/dseqrecord.py`: `answer.id = "{name}_lin".format(name=self.name)` with
`answer.name = answer.id[:16]` (linearize, L920-922) and `answer.name = "{}_rc".format(self.name[:13])`
(reverse complement). So `1011bp_PCR_prod`, `pUC19_lin`, `pUC19_rc` appear automatically, and PCR products
inherit a feature label from the template when a feature spans the entire template (`amplify.py`:
`f.location.start == 0 and f.location.end == len(self.template)`). There is no naming scheme
encoding a variant's design (position / mutation / replicate) — because there are no variants.
*Source:* `src/pydna/amplify.py`, `src/pydna/dseqrecord.py` L920-922, `src/pydna/utils.py` L349 (`identifier_from_string`).

**`design_visualization` — yes. [EVIDENCE STRENGTHENED]**
Five distinct mechanisms, all documented:
- `Dseqrecord.figure()` / `Amplicon.figure()` — ASCII/HTML-marked double-stranded figures and
  primer-annealing figures. `figure(fig_type=...)`: "If the source is an `Assembly` object, it returns a
  compact figure of the assembly … If the source is a `PCRSource` object, it returns a figure of the PCR
  alignment."
- `Contig.figure()` (contig.py L134) — "Compact ascii representation of the assembled fragments"; and
  `Contig.detailed_figure()` (L80) — "Returns a text representation of the assembled fragments", with worked
  linear and circular examples. Paper: "The figure method gives an text based figure outlining how the
  sequences were assembled for rapid inspection."
- **`Contig.figure_mpl()` (contig.py L272) — MISSED BY THE ORIGINAL EXTRACTION, re-verified here.** Docstring
  verbatim: *"Graphic representation of the assembly. Returns — matplotlib.figure.Figure — A representation
  of a linear or culrcular assembly."* [sic]. Implementation lazily imports `matplotlib.pyplot` and
  `matplotlib.patches`, picks per-edge colours from `tab20`, and draws fragments on concentric arcs for
  circular assemblies. **So pydna has a genuine graphical assembly view, not only ASCII plus a gel PNG.**
  (A `figure_plotly` is stubbed out at L405 but commented, i.e. not shipped.)
- `pydna.gel.gel(samples=..., gel_length=600, ...)` renders an actual gel image with a molecular-weight
  ladder via NumPy and PIL, with SciPy for the migration interpolation (README: "Gel electrophoresis of DNA
  with generation of gel images"). Matplotlib is used by `Contig.figure_mpl()`, not by `gel()`.
- `Dseqrecord.history()` tree (above) and `CloningStrategy.model_dump_json()` → loadable in the **OpenCloning**
  web interface for a graphical cloning-strategy view.
Paper additionally: homologies "are added to each sequence as metadata in the form of Genbank features
(Figure 1A) which can be inspected graphically using a sequence editor."
**Caveat for the table:** everything above visualises *one construct / one assembly*. Nothing visualises a
library.
*Source:* `src/pydna/contig.py` L80/L134/L272-L300, `src/pydna/dseqrecord.py`, `src/pydna/gel.py`, `README.md`, paper Fig. 1A.

### Block B — assay coverage

**`assay_dms` — no.**
No deep-mutational-scanning, saturation-mutagenesis or codon-scanning functionality anywhere. `mutagen*`
appears on exactly one docs page in the whole index (`example_gallery`, the external teemi link). pydna does
understand coding sequences — `seqrecord.isorf()`, `translate()`, `seq.orfs()`/`orfs2()`, and
`Seq.cai()/rarecodons()/startcodon()/stopcodon()/express()` (seq.py L142-L172) — but these **diagnose one
ORF** and generate **no variants**. Searched paper, README, docs gallery, all module names and the full
820-object API index: absent.
*Source:* `src/pydna/seq.py` L142-L228, `src/pydna/seqrecord.py`, docs `searchindex.js`, paper.

**`assay_mpra` — no.**
No MPRA / STARR-seq / oligo-pool / promoter-enhancer-tiling module, function, or docs page. Confirmed
against the full docs object index and the `example_gallery` source. The gallery lists "04 Promoter library
design in *Saccharomyces cerevisiae*" only under "Bioengineering/SynBio projects where you can apply pydna"
as an external link; it is not a pydna module or a documented MPRA workflow.
*Source:* `docs/example_gallery.rst`, docs `searchindex.js`.

**`assay_insilico` — no.**
No machine-learning or model-probing functionality of any kind. `pyproject.toml` re-read in full for this
memo: required dependencies are `biopython ^1.87`, `networkx >=2.8.8`, `numpy`, `prettytable >=3.5.0`,
`pydivsufsort ^0.0.20`, `pyfiglet`, `seguid >=0.0.5`, `regex`, `opencloning-linkml ^1`, `sgffp >=0.18.0`,
`appdirs >=1.4.4`; optional extras are `pyperclip` (clipboard), `pyparsing`+`requests` (download), `cai2`
(express), `scipy`+`matplotlib`+`pillow` (gel), `pyahocorasick` (primer_screen). **No torch, no tensorflow,
no jax, no sklearn, no model API.**
*Source:* `pyproject.toml` (full, 129 lines, re-read 2026-08-11).

### Block C — genomics integration

**`genome_coordinates` — partial.**
`pydna.genbank.Genbank.nucleotide(item, seq_start=None, seq_stop=None, strand=1)` (genbank.py L69) accepts an
accession-plus-region string. Docstring examples, verbatim:

```
| BK006936.2 REGION: complement(613900..615202)
| NM_005546 REGION: 1..100
| NM_005546 REGION: complement(1..100)
| 21614549:1-100
| 21614549:c100-1
| 21614549 1-100
| 21614549 c100-1
```

`BK006936.2` is an *S. cerevisiae* chromosome record, so chromosome-coordinate retrieval is genuinely
supported. **But** only via NCBI accessions over the network: there is no genome-build concept, no local
genome FASTA/2bit, and no chromosome-name (`chr7:1000-2000`) interface. The provenance model has
`GenomeCoordinatesSource(assembly_accession, locus_tag, gene_id, coordinates)`, but that is a *record type*
only; pydna implements no retrieval by `assembly_accession` through that model.

**Corroborating item missed by the original extraction (re-verified here):**
`pydna.sequence_picker.genbank_accession(s)` (sequence_picker.py L17) BLASTs a sequence against NCBI `nt` and
returns a `Dseqrecord` whose description is built as
`f"{best_alignment.accession} REGION: {start}..{stop}"` (L47) — i.e. the **reverse mapping**, sequence →
genomic coordinates. Strengthens `partial`; does not reach `yes`.
*Source:* `src/pydna/genbank.py` L69-L99, `src/pydna/sequence_picker.py` L17-L47, `src/pydna/opencloning_models.py`.

**`transcript_models` — partial (GFF/GTF specifically: no). [CHANGED from `no`]**

> **Changed on adversarial review; the reviewer's verdict was "understated" and I re-verified it from the
> notebook source with executed outputs before accepting it.**

The original extraction's *literal* claims all hold and are retained: `src/pydna/parsers.py` handles
GenBank, EMBL and FASTA (`embl_gb_fasta()`) plus SnapGene (`parse_snapgene()`) and primers
(`parse_primers()`); **there is no GFF/GTF parser and nothing imports one** (`gff` and `gtf` both return 0
hits across the entire rendered docs index, API reference included); there is no `Transcript` object and no
isoform-selection logic.

**But "annotation comes from GenBank feature tables only" *is* a transcript model.** INSDC `mRNA` and `CDS`
features with `join(...)` locations parse into `Bio.SeqFeature.CompoundLocation`, and pydna **documents
constructing, exporting and extracting them in its own getting-started tutorial**. From
`docs/notebooks/Dseq_Features.ipynb` (downloaded and parsed cell-by-cell for this memo):

- Cell 0 (markdown): "Examples include coding sequences (CDS), introns, promoters, etc."
- Cell 9 (markdown) reproduces a real *S. pombe* EMBL record with
  `FT   CDS             join(1878362..1878785,1878833..1880604)` as the worked motivation.
- Cell 12-13 build a two-exon CDS with `CompoundLocation` and export it; executed output contains
  `CDS             join(4..9,13..15)` and `location: join{[3:9], [12:15]}`.
- Cells 18-19 handle **origin-spanning features on circular molecules** as a `CompoundLocation`
  (`join{[19:25](+), [0:6](+)}` → GenBank `join(20..25,1..6)`), with the note "This feature will be displayed
  as a single feature in SnapGene viewer and Benchling, since they support this convention."

The value is `partial` because pydna represents and operates on joined `mRNA`/`CDS` feature locations but
has no GFF/GTF input, transcript object or isoform-selection logic.
*Source:* `docs/notebooks/Dseq_Features.ipynb` cells 0/9/12/13/18/19 (executed outputs), `src/pydna/parsers.py`, docs `searchindex.js` (`gff` 0, `gtf` 0).

**`exon_intron_split_codons` — partial. [CHANGED from `no`]**

> **Changed on adversarial review; verdict "understated". The original sentence "nothing assembles a CDS
> across exons or handles codons split across an exon junction" is demonstrably false as written and is
> withdrawn. I re-downloaded `Dseq_Features.ipynb` and read the executed cell outputs myself.**

`docs/notebooks/Dseq_Features.ipynb`, cell 12 (markdown), verbatim:

> "### Adding a Feature with Parts — To add a feature with parts, like a CDS with introns, we need to use a
> `CompoundLocation` object when creating a `SeqFeature`. The example code below adds a CDS with two parts,
> between 3-9bp and 12-15bp, to my features list. In a real-world scenario this would represent a CDS with an
> intron that skips the `ACG` codon: ATGCGT~~ACG~~TGA"

Cell 13 code builds `CompoundLocation([FeatureLocation(3, 9), FeatureLocation(12, 15)])` into a `CDS`
`SeqFeature`, appends it to a `Dseqrecord`, and its **executed output** shows both
`location: join{[3:9], [12:15]}` and the GenBank line `CDS             join(4..9,13..15)`.
Cell 14 (markdown): *"We can even extract a protein record as follows (see how the protein sequence is `MR`,
skipping the intron)"*. Cell 15 executes
`sub_record = dummy_record.features[-1].extract(dummy_record); print(sub_record.translate())` with executed
output `ProteinSeq('MR')`.

This is a **first-class API path, not a notebook trick**: `SeqRecord.extract_feature(n)` is literally
`return self.features[n].extract(self)` (`src/pydna/seqrecord.py` L366), inherited by
`Dseqrecord.extract_feature()` (`dseqrecord.py` L275); `utils.shift_location()` (L70) and
`utils.shift_feature()` (L132) explicitly handle `CompoundLocation`, as does `utils.location_boundaries()`
(L749); GenBank `join(...)` locations round-trip through `parsers.parse()`. Because exons are concatenated in
transcript order (with strand handled) *before* translation, **a codon straddling an exon junction translates
correctly**.

Why `partial` and not `yes`: (a) pydna
performs no variant design, so it never has to *place* an edit into a junction-straddling codon — the hard
part of this capability in a design tool is simply never exercised; (b) the tutorial's own example uses exon
lengths that are multiples of three (6 bp + 3 bp), so the notebook demonstrates intron skipping but does not
itself exhibit a split codon — correctness for split codons follows from the concatenate-then-translate
mechanism rather than from a shipped worked example.
*Source:* `docs/notebooks/Dseq_Features.ipynb` cells 12-15 (executed outputs), `src/pydna/seqrecord.py` L366, `src/pydna/dseqrecord.py` L275, `src/pydna/utils.py` L70/L132/L749.

**`hgvs_input` — no.**
Independently confirmed twice: **zero occurrences of `hgvs`** anywhere in the rendered documentation
including the full API reference (`searchindex.js` audit), and no `hgvs` dependency in `pyproject.toml`. The
repo-wide code-search hits — 31, in two tracked files: `tests/pUC_LAC4_correct_rotation.gb` (1) and
`docs/notebooks/gibson_eg.gb` (30) — are all incidental `HGVS` substrings inside protein `/translation=`
qualifiers, not a parser. No mention in the paper or docs.
*Source:* docs `searchindex.js`, `pyproject.toml`, `tests/pUC_LAC4_correct_rotation.gb`.

**`vcf_vep_output` — no.**
Zero occurrences of `vcf` across the entire docs index. Output formats are GenBank / FASTA / EMBL
(`Dseqrecord.format(format="gb")`, `.write()`), gel PNG images, and the OpenCloning `CloningStrategy` JSON.
**There is no variant representation in the data model at all**, so there is nothing to emit — the absence is
structural, not an omission. No VEP interoperability.
*Source:* docs `searchindex.js`, `src/pydna/seqrecord.py` (`format`/`write`), `src/pydna/opencloning_models.py`.

**`consequence_annotation` — no.**
No molecular-consequence vocabulary (stop-gained / synonymous / missense / frameshift / in-frame indel)
anywhere in the 820 documented objects. The nearest features are `seqrecord.isorf()` ("Detect if sequence is
an open reading frame (orf) in the 5'-3' direction") and
`translate()`, which leave all interpretation to the user.
*Source:* docs `searchindex.js` (full object index), `src/pydna/seqrecord.py`, `src/pydna/seq.py`.

### Block D — physical construction

**`primer_design` — yes. [EVIDENCE SOFTENED — scope limit widened]**
Value unchanged. `src/pydna/design.py` provides:
- `primer_design()` — "This function designs a forward primer and a reverse primer for PCR amplification of
  a given template sequence"; can match a target Tm or an existing primer's Tm;
- `assembly_fragments()` — *"Adds tails to primers for a linear assembly through homologous recombination or
  Gibson assembly"*;
- `circular_assembly_fragments()` (deprecated alias for `assembly_fragments(circular=True)`);
- `user_assembly_design()` (experimental, USER cloning).

Paper: "Pydna provides the pydna assembly_primers function in order to automatically design tailed primers
for a series of DNA fragments … The algorithm tries to create primers with balanced melting temperature for
the annealing region." `Amplicon` objects "store rich information about the PCR simulation, such as the DNA
region where the primer anneals, melting temperature of each primer and also a suggested PCR program"
(`src/pydna/tm.py` — verified to contain `tm_default`, `tm_dbd`, `tm_product`, `ta_default`, `ta_dbd`,
`program`, `dbd_program`, `Q5`, `tmbresluc`, `tm_neb`).

**The original scope limit "assembly/amplification primers only" was too narrow and is withdrawn.**
`src/pydna/primer_screen.py` (re-verified: `forward_primers` L290, `reverse_primers` L375, `primer_pairs`
L461, `flanking_primer_pairs` L573, `diff_primer_pairs` L648, `diff_primer_triplets` L760) performs
**diagnostic / verification primer-pair SELECTION across a set of sequences**. `diff_primer_pairs` docstring,
verbatim:

> "Primer pairs for diagnostic PCR. Given an iterable of sequences and a primer list, primers are selected
> that result in unique product sizes from each of the input sequences. … Primers 1 and 2 could be used to
> verify genetic modifications such as cloning an insert into a plasmid vector. … The callback function is
> used to return true or false for the PCR products. This score is meant to filter for PCR products that are
> likely to migrate to sufficiently distinct locations to be distinguishable on a typical agarose gel. Only
> products larger than `short` and smaller than `long` are returned."

Signature: `diff_primer_pairs(sequences, primer_list, short=500, long=1500, limit=16, automaton=None,
callback=callback)`. This is genotyping/verification primer selection over a *collection*, and it is one of
the four headline docs examples (`Example_primer_screen`). It ships behind the optional
`primer_screen`/`pyahocorasick` extra and warns on import that it "is experimental and not yet extensively
tested. api may change in future versions."

**The true and only limit:** there is **no mutagenic-primer designer** — nothing generates QuikChange-style
overlapping mutagenic primer pairs from a desired mutation set. (pydna can *simulate* mismatched-primer PCR
via `mismatches=N`; it does not *design* the mismatched primer.)
*Source:* `src/pydna/design.py`, `src/pydna/tm.py`, `src/pydna/primer_screen.py` L290-L790, paper (Implementation), `docs/example_gallery.rst`.

**`codon_optimization` — no. [EVIDENCE SOFTENED]**
Value defensible and unchanged. `src/pydna/codon.py` is **data only** — per-organism codon `weights`,
`start`/`stop` codon frequency distributions, `rare_codons` lists, and N-end-rule half-life data (yeast `sce`
is the only fully populated organism; `rare_codons` also ships a populated `eco` list, while `weights` has no
`eco` key and `start`/`stop`/`n_end` have empty `eco` entries). `src/pydna/utils.py` provides `cai()` (L225,
computes the codon adaptation index against the per-organism weights table), `rarecodons()` (L235, returning
rare-codon **positions** as `List[slice]`), and `express()` (L248) — all three carry only the placeholder
docstring `"""docstring."""`. The same diagnostics are exposed as `Seq.cai()`,
`Seq.rarecodons()`, `Seq.startcodon()`, `Seq.stopcodon()`, `Seq.express()` (seq.py L146-L172; `express`
requires the optional `cai2` extra), though `Seq.express()` is a separately implemented method returning a
`PrettyTable`, not a wrapper around `utils.express`. **These diagnose codon usage. No function anywhere
back-translates a protein or rewrites a CDS**, and `utils.express` is marked "**NOT IMPLEMENTED YET**" and
raises `NotImplementedError`. `optimi*` returns 0 hits across the entire docs index.

**Required wording caution.** Phrase the cell as
*"pydna provides codon-usage diagnostics (CAI, rare-codon localisation, start/stop frequency, N-end rule) but
no codon-optimizing function"* — **not** "pydna cannot codon-optimize".
*Source:* `src/pydna/codon.py`, `src/pydna/utils.py` L225-L275, `src/pydna/seq.py` L142-L172, `pyproject.toml`, docs `searchindex.js`.

**`synthesis_constraints` — partial. [CHANGED from `no`]**

> **Changed on adversarial review after `dseq.py` L1856-L1900 and `dseqrecord.py` L925-L953 were
> re-verified in full.**

**Pydna exposes restriction-site and sequence-property queries as documented methods on its core objects**
(all re-verified in source for this memo):

- `Dseq.no_cutters()` (dseq.py L1856) — "Enzymes in a RestrictionBatch **not** cutting sequence";
  `Dseq.unique_cutters()` (L1866) — "cutting sequence once", implemented as `n_cutters(n=1)`, with
  `once_cutters = unique_cutters` declared as an explicit alias (L1874); `Dseq.twice_cutters()` (L1876);
  `Dseq.n_cutters(n=3, batch=None)` (L1884) — "Enzymes in a RestrictionBatch cutting n times";
  `Dseq.cutters()` (L1894) — "cutting sequence at least once". **All five default to
  `CommOnly`**, the commercially-available REBASE enzyme set, when no batch is supplied.
- The same six are re-exported on the record type: `Dseqrecord.no_cutters/unique_cutters/once_cutters/
  twice_cutters/n_cutters/cutters` (dseqrecord.py L925-L947), each delegating with `batch or CommOnly`.
- `Dseqrecord.number_of_cuts(*enzymes)` (L949) — "The number of cuts by digestion with the Restriction
  enzymes contained in the iterable", verbatim
  `return sum([len(enzyme.search(self.seq)) for enzyme in flatten(enzymes)])`.
- `Seq.gc()` (seq.py L142) and `SeqRecord.gc()`.
- `Seq.rarecodons()` (seq.py L152) returning `List[slice]` — the *positions* of rare codons.
- `pydna.tm` — `tm_default` (L30), `tm_dbd` (L73), `tm_product` (L116), `ta_default` (L129), `ta_dbd` (L144),
  `Q5` (L313), `tmbresluc` (L321), `tm_neb` (L371, which queries the NEB web API), plus `program()` (L148)
  and `dbd_program()` (L212) which emit a suggested PCR cycling program.

"Does my designed construct contain an unwanted site for enzyme X, and which enzymes cut it zero / once /
twice?" is a **one-call query** in pydna.

Why `partial` and not `yes`: **nothing is enforced**. There is no homopolymer rule, no repeat-content
filter, no hairpin / secondary-structure check, no GC-window scan (only a single whole-sequence
`Seq.gc()` float), no oligo-length or vendor/pool constraint, and no "reject or repair this sequence" API.
`sequence_regex.py` compiles IUPAC-degenerate motifs into regexes for recognising recombinase/enzyme sites,
not for constraint filtering. **Carry this sentence with the cell:** *"site-content and Tm queries only;
nothing is enforced and there is no repeat / homopolymer / GC-window / length rule."*

*Source:* `src/pydna/dseq.py` L1856-L1900, `src/pydna/dseqrecord.py` L925-L953, `src/pydna/seq.py` L142/L152, `src/pydna/tm.py` (all functions), `src/pydna/sequence_regex.py`.

### Block E — engineering

**`interface` — yes (Python API / library).**
Python API only (`import pydna`), intended for scripts and Jupyter notebooks. `pyproject.toml` was re-read in
full (129 lines) for this memo: it declares **no `[project.scripts]` and no console entry points**, so
`pip install pydna` installs **no command-line tool**. The classifier `"Environment :: Console"` refers to
REPL use, not to a CLI. The current repository ships no GUI or web server.
The paper's phrase "designed purely as a command line tool" (Comparison with existing tools) is quoted
accurately but means *scripting / non-GUI*, not an installed executable — gloss it that way if it appears in
the manuscript.
Pydna's own `opencloning_models.py` documentation says its `CloningStrategy` JSON can be loaded in an
external web interface; that does not make the interface part of pydna.
*Source:* `pyproject.toml` (full), `README.md`, `src/pydna/opencloning_models.py`, paper (Comparison with existing tools; Availability).

**`license` — yes (BSD 3-Clause, OSI-approved, permissive).**
`LICENSE.txt` header: "Copyright (c) 2013-2026 Björn Johansson, Centro de Biologia Molecular e Ambiental
(CBMA)/Department of Biology, University of Minho, Braga, Portugal", followed by the three standard BSD
clauses including "Neither the name of the organizations … may be used to endorse or promote products".
All 40 `src/pydna/*.py` files on `master` carry an SPDX header: 39 are `BSD-3-Clause`, and the vendored
`threading_timer_decorator_exit.py` (imported by `assembly.py`) is `SPDX-License-Identifier: MIT` — also
permissive. The released 5.5.16 sources carry no SPDX headers at all. `pyproject.toml`:
`license-files = ["LICENSE.txt"]` and `[tool.poetry] license = "BSD"`. PyPI classifier:
`"License :: OSI Approved :: BSD License"`.
Two footnotes: the 2015 paper says "License: FreeBSD" (2-clause) whereas the shipped text is 3-clause; and
GitHub's licence detector reports `NOASSERTION` because the file is named `LICENSE.txt` with a custom
preamble. Neither affects the value.
*Source:* `LICENSE.txt`, `pyproject.toml`, PyPI JSON classifiers, paper (Availability).

**`maintained` — yes (actively, dual maintainership).**
Re-verified from the GitHub repository API and PyPI on 2026-08-14:
- The `master` head remains
  **2026-08-01T09:24:09Z**, with further commits 2026-07-16, 2026-07-15, 2026-07-12.
- PyPI JSON — current version **5.5.16, uploaded 2026-06-09T06:29:56Z**; the preceding releases
  5.5.9 (2026-03-31), 5.5.10 (2026-04-10), 5.5.11 (2026-04-20), 5.5.12 (2026-06-02), 5.5.13 (2026-06-03),
  5.5.14 (2026-06-04), 5.5.15 (2026-06-05) show a tight cadence.
- `master` `pyproject.toml` version string `6.0.0-a.24.post.17+b7b559bd66` — a v6 major is in development.
- `requires_python = ">=3.10,<4.0"`; classifiers list Python 3.10, 3.11, 3.12, 3.13, **3.14**;
  `Development Status :: 4 - Beta`.
- Two authors declared in `pyproject.toml`: Björn F. Johansson and Manuel Lera-Ramirez.
- Repo `pushed_at` **2026-08-11T14:09:17Z**, 228 stars, 67 open issues **and** PRs (GitHub's aggregate
  `open_issues_count`; the live split is 61 issues / 6 PRs), not archived, CI (tests + coverage +
  docs publish) workflows present.
*Source:* GitHub repository and commit APIs, PyPI JSON, `pyproject.toml`, repository workflows — rechecked 2026-08-14.

---

## Capability summary table

| Key | Value | One-line basis |
|---|---|---|
| library_as_object | **no** | No variant-library abstraction; `CloningStrategy` collects one strategy's records/provenance, while assembly getters return plain `list[Dseqrecord]` |
| dag_chaining | **partial** | Typed `source` DAG + `history()` recorded eagerly; `CloningStrategy` round-trips only over its own embedded inputs, no replay on new inputs |
| lazy_evaluation | **no** | `parse()` "a greedy function"; all products materialised; generators are explosion guards only |
| mixed_mutagenesis_one_pool | **no** | Can *simulate* one edit (`mismatches=N` PCR, CRISPR+oligo); no mutagenesis *design*, no pool |
| combinatorial_motif_place | **partial** | `Assembly` graph enumerates compatible fragment orders/orientations under configured bounds; homology-driven, no placement API |
| barcode_generation | **no** | `barcode` 0 hits in whole docs index; `expand_iupac_to_dna()` expands degenerate oligos but with no GC/distance constraint |
| per_sequence_provenance | **yes** | Typed pydantic `source` tree, `history()`, `validate_history()`, `normalize_history()`, SnapGene history *import* |
| automatic_naming | **partial** | `1011bp_PCR_prod`, `_lin`, `_rc` auto-names, 16-char truncated, operation-tagged not design-descriptive |
| design_visualization | **yes** | `figure()`, `detailed_figure()`, **`figure_mpl()` → matplotlib Figure**, `gel.gel()` image, `history()` tree, OpenCloning JSON |
| assay_dms | **no** | No DMS/saturation/codon-scanning; `mutagen*` on one docs page (external link) |
| assay_mpra | **no** | No MPRA/STARR/oligo-pool/tiling module or docs page |
| assay_insilico | **no** | No ML dependency or model-probing API anywhere in `pyproject.toml` |
| genome_coordinates | **partial** | `Genbank.nucleotide("BK006936.2 REGION: complement(613900..615202)")`; NCBI accessions only, no build/local genome |
| transcript_models | **partial** *(GFF/GTF: no)* | GenBank/EMBL `join()` → `CompoundLocation`, documented; no GFF/GTF parser, no transcript object |
| exon_intron_split_codons | **partial** | Tutorial builds spliced CDS, exports `join(4..9,13..15)`, extracts and translates to `MR`; never used for edit placement |
| hgvs_input | **no** | 0 occurrences of `hgvs` in the whole docs index; no `hgvs` dependency |
| vcf_vep_output | **no** | 0 occurrences of `vcf`; no variant representation exists to emit |
| consequence_annotation | **no** | No consequence vocabulary among 820 documented objects |
| primer_design | **yes** | `primer_design`, `assembly_fragments`, `user_assembly_design`, full `tm` module, **plus diagnostic primer-pair selection across sequence sets**; no mutagenic-primer designer |
| codon_optimization | **no** | CAI / rare-codon diagnostics only; `codon.py` is data; `utils.express` "**NOT IMPLEMENTED YET**"; `optimi*` 0 hits |
| synthesis_constraints | **partial** | `no_cutters/unique_cutters/once_cutters/twice_cutters/n_cutters/cutters/number_of_cuts` over `CommOnly` + `gc()` + `tm` module; nothing enforced |
| interface | **yes** | Python library API; no `[project.scripts]`, CLI, GUI or web server is shipped |
| license | **yes** | BSD 3-Clause, SPDX headers on every `master` source file (one vendored file is MIT), OSI classifier |
| maintained | **yes** | master commit 2026-08-01; PyPI 5.5.16 (2026-06-09); v6 alpha in dev; Python 3.10-3.14; two maintainers |

---

## pydna's own documented examples and vignettes

**Getting-started tutorial notebooks** (`docs/notebooks/`, rendered at https://pydna-group.github.io/pydna):

| Notebook | Content |
|---|---|
| `Dseq.ipynb` | "Representing sequences in pydna" — the `Dseq` double-stranded/circular data model |
| `Dseq_Features.ipynb` | "Working with Features using the Dseqrecord class" — includes the spliced-CDS `CompoundLocation` walkthrough (`join(4..9,13..15)` → `ProteinSeq('MR')`) and origin-spanning features on circular molecules |
| `Importing_Seqs.ipynb` | "Importing and viewing sequence files in pydna" |
| `Restrict_Ligate_Cloning.ipynb` | "Restriction and Ligation" |
| `PCR.ipynb` | PCR simulation |
| `primer_design.ipynb` | Primer design |
| `Gibson.ipynb` | Gibson assembly |
| `CRISPR.ipynb` | CRISPR-Cas9 |
| `history.ipynb` | "Cloning history / Cloning strategy" — provenance, `history()`, `CloningStrategy` JSON → OpenCloning |
| `new_assemblies.ipynb` | How the `assembly2` engine represents Gibson / restriction-ligation / Golden Gate / PCR / HR assemblies, incl. the "two possible circular assemblies … one with the second fragment inverted" enumeration |
| `readme_example.ipynb`, `hackathon_copenhagen_2025.ipynb` | Walkthroughs |

**Example gallery** (`docs/example_gallery.rst`, verbatim bullets — pydna's own, no external package required):
- **Example_Restriction** — "PCRing a gene out of the genome, and cloning into a vector using restriction and ligation."
- **Example_Gibson** — "Gibson assembly of *R. cellulolyticum* genomic fragments into a plasmid, from the original Gibson assembly paper (doi: 10.1038/nmeth.1318)."
- **Example_CRISPR** — "Using CRISPR with homologous recombination to delete genes by making two cuts in the genome, and repair it with an oligo. Used in the industrially relevant *K. phaffi*."
- **Example_primer_screen** — "Using the Aho-Corasick algorithm to quickly screen a primer list for annealing positions in a sequence."

**From the 2015 paper** (files still shipped in `docs/cookbook/`):
- **YEp24PGK_XK** — BamHI/BglII digestion + ligation, 12 lines of code, 11,452 bp product.
- **pGUP1** — construction by homologous recombination, 11 lines, 9,981 bp product.
- **Lactose pathway** — two-gene assembly (LAC4 + LAC12 with PDC1/PGI1/TPI1 promoters) described by "nineteen short pydna scripts", experimentally validated in *S. cerevisiae*. These scripts are in the paper's **Additional file 2**, not in `docs/cookbook/`.

**External examples and projects linked by the official docs are NOT pydna core functionality.** The docs
explicitly say that the first two examples use the extra `teemi` dependency:
- "Example_HT_cazyme_primer_design: Design primers for a high-throughput CAZyme library."
- "Example_designing_primers_for_kozak_library: We explore the combinatorial space of the most abundant kozak sequences and make repair-primers for the experiments."
- "02 Primer directed mutagenesis with pydna"
- "03 Investigate plastic-degrading enzymes with pydna"
- "04 Promoter library design in *Saccharomyces cerevisiae*"

**Reproducibility note.** None of pydna's own examples is a variant-library design; a pydna user must supply
the variant enumeration and orchestration outside the package.

---

## Additional capabilities not covered by the current row list

Each item below was re-verified in pydna's source or documentation. Items marked **[NEW]** were added on
adversarial review.

1. **Molecular-biology reaction simulation** — restriction digestion with correct sticky/blunt-end
   bookkeeping, ligation compatibility, PCR (incl. tailed primers and inverse PCR on circular templates),
   Gibson, In-Fusion, Golden Gate, Gateway (BP/LR), Cre-lox, generic recombinases, fusion PCR, oligo
   hybridization, polymerase extension, CRISPR cut + HR repair.
2. **Double-stranded / circular topology as a first-class data model** — `Dseq` tracks watson/crick strands,
   5′/3′ overhangs and circularity (`Dseqrecord(..., circular=True)`), with `looped()`, `synced()`,
   `shifted()`, `fill_in()`, `T4()`, `nibble_*()`, `terminal_transferase()`.
3. **Structured provenance / cloning-strategy serialisation** — typed `source` models + `history()` tree +
   `CloningStrategy` JSON interoperable with the OpenCloning web app + `validate_history()` replay
   verification **+ `normalize_history()` [NEW]**.
4. **[NEW] Provenance IMPORT from a commercial GUI** — `snapgene_history_parser.parse_snapgene_history()`
   reads a SnapGene `.dna` file, resolves every recorded step (re-executing it through the corresponding
   pydna assembly function), and returns a `Dseqrecord` whose `source` tree mirrors the full cloning
   provenance; raises `ValueError` if a recorded operation cannot be reproduced. Uses the `sgffp` dependency.
   Unsupported histories include GC/TA/TOPO cloning, destroyed restriction fragments, blunting-dependent
   removal and manual insertion/replacement/reverse-translation/codon/site edits.
5. **Round-trip with lab file formats and repositories** — GenBank/EMBL/FASTA/SnapGene parsing, GenBank
   download by accession+region, `genbankfixer` for malformed records, Addgene/Benchling/SEVA/iGEM/Euroscarf
   source records, `seguid` checksums, ApE/SnapGene-friendly feature colouring, clipboard copy;
   `utils.eq()` compares linear/reverse-complement and optional circular-rotation equivalence,
   `utils.deduplicate()` removes repeated items, and `Dseqrecord.map_trace_files()` reads ABI traces, applies
   an optional orientation-aware target filter, maps direct common substrings and adds trace features.
6. **Virtual gel electrophoresis** — `pydna.gel.gel()` renders a gel image with a molecular-weight ladder;
   built-in ladder lists use length/mass/mobility `FakeSeq` records without full nucleotide sequences.
7. **[NEW] Graphical assembly rendering** — `Contig.figure_mpl()` returns a `matplotlib.figure.Figure`
   drawing a linear or circular assembly (per-fragment colours from `tab20`, concentric arcs when circular);
   `Contig.detailed_figure()` gives the aligned text view.
8. **CRISPR protospacer/PAM search** — `pydna.crispr` with a `_cas` abstract base and enzyme classes.
9. **[NEW] Diagnostic / verification primer-pair SELECTION across a set of sequences** —
   `primer_screen.diff_primer_pairs`, `diff_primer_triplets`, `flanking_primer_pairs`, `primer_pairs`,
   `forward_primers`, `reverse_primers`: given a sequence set and a primer stock, "primers are selected that
   result in unique product sizes from each of the input sequences", with `short`/`long` product-size windows
   and a gel-distinguishability `callback`. Aho-Corasick automaton, savable/loadable.
10. **[NEW] Degenerate-oligo expansion** — `primer_screen.expand_iupac_to_dna()` cartesian-expands an
    extended-IUPAC oligo into every concrete DNA sequence it encodes
    (`"ATNG"` → `['ATGG','ATAG','ATTG','ATCG']`; `"ACGTURYSWKMBDHVN"` → 20,736 sequences).
11. **[NEW] PCR simulation with mismatched primers** — `assembly2.pcr_assembly(..., mismatches=N)` /
    `PCRAssembly(..., mismatches=N)` via `primer_template_overlap`'s fuzzy regex `"(...){s<=N}"`, i.e.
    simulation of site-directed mutagenesis by primer mismatch (substitutions only, no indels).
12. **[NEW] Restriction-site content queries as first-class methods** —
    `no_cutters()`, `unique_cutters()`/`once_cutters()`, `twice_cutters()`, `n_cutters(n)`, `cutters()` over
    the `CommOnly` REBASE set, plus `Dseqrecord.number_of_cuts(*enzymes)`.
13. **[NEW] Multi-part (spliced) CDS features** via `Bio.SeqFeature.CompoundLocation`, documented end-to-end
    in `Dseq_Features.ipynb` with GenBank `join()` export, `extract_feature()` splicing and correct
    translation — including origin-spanning features on circular molecules.
14. **[NEW] `pydna.tm` — multiple Tm models** (`tm_default`, `tm_dbd`, `tm_product`, `ta_default`, `ta_dbd`,
    `Q5`, `tmbresluc`, `tm_neb` which queries the NEB web API) plus `program()` / `dbd_program()` emitting a
    suggested PCR cycling program.
15. **[NEW] Feature/location arithmetic across a circular origin** — `utils.shift_location()`,
    `utils.shift_feature()`, `utils.location_boundaries()`, `Dseqrecord.shifted()`, `Dseqrecord.synced()`;
    handles `CompoundLocation`. Relevant because it is how design metadata survives circular permutation.
16. **Codon-usage and protein diagnostics** — CAI, rare-codon *positions*, start/stop codon frequency,
    N-end-rule half-life (`Seq.express()`), yeast tables shipped; translation returns
    `ProteinSeq`/`ProteinSeqRecord`, and `ProteinSeq` reports molecular weight, pI and instability index.

---

## Stated limitations (pydna's own scope, plus what it does not do)

pydna's authors scope it as a **simulation and documentation tool for cloning strategies**, and the record
below is a description of that scope, not a criticism of it.

- **No variant-library abstraction.** `CloningStrategy` collects the records and provenance for one cloning
  strategy, but it is not a designed variant set; there is no mutagenesis design operation (no scanning, no
  exhaustive singles or pairs, no random sampling, no WT-replicate concept, no way
  to mix mutagenesis types into one pool); no barcode generation with GC or edit-distance constraints; no
  user-facing motif-placement specification. Multi-molecule results are lists of *alternative reaction
  products*, not designed sets.
- **Simulation, not design, of edits.** `mismatches=N` PCR and CRISPR + oligo repair simulate an individual
  designed edit; nothing generates the mutagenic primer or enumerates the edits to make.
- **Eager execution.** `parsers.parse()` is documented as "a greedy function, use carefully"; assemblies
  materialise every product. Internal generators are explosion guards, not user-facing laziness.
- **Provenance is cloning provenance, not design-vs-WT provenance.** The `source` tree records *how a molecule
  was built from other molecules*; it has no concept of a reference/wild-type sequence and a designed
  deviation from it.
- **The DAG is a record, not a program.** `CloningStrategy` serialises; it cannot be re-executed on new inputs
  inside pydna. `validate_history()` replays only against the
  original inputs.
- **No variant file formats.** No HGVS input, no VCF/VEP output, no consequence vocabulary. There is no
  variant representation in the data model at all.
- **No GFF/GTF.** Annotation enters only through GenBank/EMBL feature tables (and SnapGene). No transcript
  object, no isoform selection.
- **Genome access is network-bound and accession-based.** No genome build concept, no local FASTA/2bit, no
  `chr7:1000-2000` interface.
- **Nothing is enforced.** Restriction-site, GC and Tm queries are available, but there is no
  repeat/homopolymer/hairpin/GC-window/oligo-length rule and no "reject or repair" API.
- **No codon-optimizing function.** Codon usage is diagnosed, never rewritten; `utils.express` is marked
  "**NOT IMPLEMENTED YET**" and raises `NotImplementedError`.
- **No machine learning.** No ML dependency and no model-probing API.
- **No CLI, no GUI, no web server is shipped.** Library only; the 2015 paper separately described a hosted
  interactive environment.
- **Naming is operation-tagged and 16-char truncated**, not design-descriptive.

---

## Availability status (as of 2026-08-14)

- **Installable today: YES.** `pip install pydna` → **5.5.16**, uploaded **2026-06-09**. Pure Python.
  `requires_python = ">=3.10,<4.0"`; classifiers cover Python 3.10-3.14. Optional extras: `clipboard`,
  `download`, `express`, `gel`, `primer_screen`.
- **Repo alive: YES.** github.com/pydna-group/pydna — `master` head commit **2026-08-01T09:24:09Z**, further
  commits 2026-07-16 / 07-15 / 07-12; repo `pushed_at` **2026-08-11T14:09:17Z**; not archived; CI workflows
  present; 228 stars;
  61 open issues + 6 open PRs (GitHub's `open_issues_count` aggregate of 67 counts both);
  **two active maintainers** (Björn F. Johansson, Manuel Lera-Ramirez). A **v6 major is in
  development** (`master` version string `6.0.0-a.24.post.17+b7b559bd66`).
- **Release cadence: healthy.** 5.5.9 (2026-03-31) → 5.5.10 (04-10) → 5.5.11 (04-20) → 5.5.12 (06-02) →
  5.5.13 (06-03) → 5.5.14 (06-04) → 5.5.15 (06-05) → 5.5.16 (06-09).
- **Docs alive: YES.** https://pydna-group.github.io/pydna (HTTP 200), full API reference + notebook gallery.
  Legacy https://pydna.readthedocs.io/latest/ still resolves but is not the maintained site.
- **Web service shipped by current pydna: NO.** The repository contains a Python library, not a web server;
  the 2015 paper described a separately hosted interactive environment.
- **Last `master` commit: 10 days before this survey** (2026-08-01). **Last repository push: the survey
  date** (2026-08-11T14:09:17Z). These labels are kept separate because `pushed_at` need not identify a new
  `master` commit.

---

## Unresolved disagreements

**None.** Every one of the 24 keys is settled:
- 21 keys retained their extracted values after primary-source verification.
- 3 keys (`transcript_models`, `exon_intron_split_codons`, `synthesis_constraints`) changed
  **`no` → `partial`** after verification against the executed feature-notebook outputs and the
  `dseq.py`/`dseqrecord.py` method definitions.
- 0 keys set to `unknown`. No cell required splitting the difference.

**Three carry-forward cautions rather than disagreements:**
1. `synthesis_constraints = partial` and `transcript_models = partial` and `exon_intron_split_codons = partial`
   rest on joined-feature handling and restriction/Tm query APIs, not on GFF/GTF import, variant placement or
   constraint enforcement; those qualifiers must remain attached.
2. `codon_optimization = no` must read "pydna provides codon-usage diagnostics but no codon-optimizing
   function", never the broader "pydna cannot codon-optimize".
3. `dag_chaining = partial` and `combinatorial_motif_place = partial` are judgement calls. Keep the clauses
   "record, not a program" and "homology-driven simulation, no placement API" attached.

## Confidence

**Highest confidence** (multiply verified, including a whole-documentation keyword audit of the Sphinx
`searchindex.js` covering the full API reference): `library_as_object`, `barcode_generation`, `hgvs_input`,
`vcf_vep_output`, `consequence_annotation`, `assay_dms`, `assay_mpra`, `assay_insilico`, `interface`,
`license`, `maintained`, `primer_design`.

**High confidence, newly corrected**: `exon_intron_split_codons`, `transcript_models` (both rest on executed
notebook outputs I re-downloaded and parsed), `synthesis_constraints` (rests on method definitions I re-read
in `dseq.py` L1856-L1900 and `dseqrecord.py` L925-L953).

**Remaining judgement calls** (state the qualifying clause alongside the value in the paper):
`dag_chaining` and `combinatorial_motif_place` — the machinery genuinely enumerates combinations and records
a DAG, but neither is a library-design abstraction; `per_sequence_provenance` — `yes`, but it is cloning
provenance, not design-variant provenance; `genome_coordinates` — `partial` hinges on accepting
NCBI accession + REGION as genome coordinates; `codon_optimization` — `no` because diagnostics do not rewrite
a coding sequence.

**Not verified by execution:** nothing was installed or run, per the standing constraint. All statements come
from the paper PDF, the rendered documentation (including notebook cells with their *stored executed
outputs*), the `master` source read over `raw.githubusercontent.com`, PyPI JSON, and GitHub repository/commit
metadata, with release and repository facts refreshed through 2026-08-14.

---

## Re-verification log for this final memo (2026-08-11)

Everything below was checked by me directly, not taken on the reviewer's word:

| Claim | How verified | Result |
|---|---|---|
| `Contig.figure_mpl()` returns a matplotlib Figure | fetched `src/pydna/contig.py`, read L272-L300 | Confirmed; docstring "Graphic representation of the assembly … Returns matplotlib.figure.Figure". `detailed_figure()` at L80, `figure()` at L134, `figure_plotly` commented out at L405 |
| Cutter methods over `CommOnly` | fetched `dseq.py` L1856-L1900 and `dseqrecord.py` L925-L953 | Confirmed all six + `number_of_cuts`. **Nuance:** `once_cutters` is an explicit alias `once_cutters = unique_cutters` on `Dseq`, and `unique_cutters` is `n_cutters(n=1)` |
| `primer_screen` diagnostic primer selection | fetched `primer_screen.py`, read `diff_primer_pairs` docstring L648-L700 | Confirmed verbatim, incl. "unique product sizes from each of the input sequences" and the gel-distinguishability callback |
| `expand_iupac_to_dna` | `primer_screen.py` L144-L162 | Confirmed, incl. both doctest examples |
| `parse_snapgene_history` | fetched `snapgene_history_parser.py`, read module docstring L1-13 and function L521-L546 | Confirmed verbatim; re-executes each step through the pydna assembly functions |
| `mismatches` in PCR assembly | `assembly2.py` L697-L790, L1947-L1988, L3244-L3280 | Confirmed, incl. the literal fuzzy regex `"(" + query + "){s<=" + str(mismatches) + "}"` |
| Spliced CDS tutorial | downloaded `docs/notebooks/Dseq_Features.ipynb`, parsed cells with stored outputs | Confirmed cells 12-15: `join{[3:9], [12:15]}`, GenBank `CDS join(4..9,13..15)`, `ProteinSeq('MR')`. Also found cells 18-19 (origin-spanning `CompoundLocation`) which neither extractor nor reviewer mentioned |
| `extract_feature` | `seqrecord.py` L366, `dseqrecord.py` L275 | Confirmed |
| `shift_location`/`shift_feature` handle CompoundLocation | `utils.py` L70, L132, L749 | Confirmed |
| `pydna.tm` model list | `tm.py` full function list | Confirmed all ten functions incl. `tm_neb` (NEB API) and `program`/`dbd_program` |
| `Seq.gc/cai/rarecodons/express` | `seq.py` L142-L172 | Confirmed; `rarecodons` returns `List[slice]` |
| No `[project.scripts]`; dependency list; license | `pyproject.toml` full read (129 lines) | Confirmed — no scripts table at all; deps and extras as listed; `license = "BSD"`, `license-files = ["LICENSE.txt"]`, classifier "License :: OSI Approved :: BSD License", Python 3.10-3.14 |
| Maintenance dates | PyPI JSON + commits Atom feed | Confirmed: 5.5.16 @ 2026-06-09; master head 2026-08-01T09:24:09Z; master version `6.0.0-a.24.post.17+b7b559bd66` |
| `Genbank.nucleotide` region parsing | `genbank.py` L69-L99 | Confirmed all seven documented region forms |
| `sequence_picker.genbank_accession` | `sequence_picker.py` L17-L47 | Confirmed; builds `f"{accession} REGION: {start}..{stop}"` |
