# pydna — evidence memo for PoolParty feature-comparison table

Survey date: 2026-08-10 (checks re-run 2026-08-11 UTC).
Surveyed by: automated literature/repo review for BMC Bioinformatics revision.

---

## 0. Identity and sources consulted

| Kind | Reference |
|---|---|
| Paper (PDF, local) | `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/papers/Pereira2015_pydna.pdf` — Pereira F, Azevedo F, Carvalho Â, Ribeiro GF, Budde MW, Johansson B. "Pydna: a simulation and documentation tool for DNA assembly strategies using python." *BMC Bioinformatics* 16:142 (2015). doi:10.1186/s12859-015-0544-x |
| Prior analysis (local) | `/mnt/c/Users/zhliu/Desktop/KinneyLab/25_poolparty_private/revision/tool_survey/prior_analyses/Pereira2015_pydna_analysis.md` |
| Current repo | https://github.com/pydna-group/pydna (the paper's URL `github.com/BjornFJohansson/pydna` now 301-redirects here; verified HTTP 200 → `https://github.com/pydna-group/pydna`) |
| Docs site | https://pydna-group.github.io/pydna (HTTP 200) |
| Old docs | https://pydna.readthedocs.io/latest/ (still HTTP 200, but the maintained docs are the GitHub Pages site linked from README) |
| PyPI | https://pypi.org/project/pydna/ |
| Source files read (raw.githubusercontent, `master`) | `README.md`, `LICENSE.txt`, `pyproject.toml`, `docs/example_gallery.rst`, `src/pydna/{design,codon,utils,parsers,amplify,dseqrecord,dseq,seq,seqrecord,assembly2,opencloning_models,gel,crispr,primer_screen,sequence_regex,fusionpcr,oligonucleotide_hybridization,sequence_picker,genbank}.py` |
| GitHub API | `repos/pydna-group/pydna`, `/commits`, `/releases`, `/contents/...`, code search |

**Important correction to the prior analysis:** the prior notes describe pydna as of the 2015 paper. The package has been very substantially extended since (Golden Gate, Gateway, Cre-lox, generic recombinases, CRISPR, fusion PCR, oligo hybridization, a structured provenance/"cloning history" data model, gel image rendering, a new `assembly2` engine, a second maintainer). Several prior-note claims need qualification — in particular **"Design cards: Not available"** is no longer accurate in spirit: pydna now stores structured per-sequence provenance in a `source` attribute and can serialise it to JSON (see §3, `per_sequence_provenance`). The prior note's core conclusion (pydna is a *cloning-simulation* tool with no notion of a variant library) is still correct.

---

## 1. What pydna is (from the paper)

> "We have developed pydna, an extensible, free and open source Python library for simulating basic molecular biology DNA unit operations such as restriction digestion, ligation, PCR, primer design, Gibson assembly and homologous recombination. A cloning strategy expressed as a pydna script provides a description that is complete, unambiguous and stable. Execution of the script automatically yields the sequence of the final molecule(s) and that of any intermediate constructs." (Abstract)

> "Most of the pydna functionality is implemented as methods for the Dseqrecord class, which was designed to hold sequence information necessary for describing a double-stranded DNA molecule." (Implementation)

> "Pydna simplifies both the planning and sharing of cloning strategies and is especially useful for complex or combinatorial DNA molecule construction." (Abstract) — note "combinatorial" here means multi-fragment assemblies, not variant libraries.

Availability section of the paper: "Project home page: https://pypi.python.org/pypi/pydna … Other requirements: Python 2.7 … License: FreeBSD".

## 2. What shipped (current repo/docs)

README feature list, verbatim:

> Pydna provides simulation of:
> - Primer design
> - PCR
> - Restriction digestion
> - Ligation
> - Gel electrophoresis of DNA with generation of gel images
> - Homologous recombination
> - Gibson assembly
> - Golden gate assembly (in progress)

Modules present in `src/pydna/` (GitHub contents API):
`all, alphabet, amplicon, amplify, assembly, assembly2, codon, common_sub_strings, contig, cre_lox, crispr, design, dseq, dseqrecord, fakeseq, fusionpcr, gateway, gel, genbank, genbankfixer, ladders, oligonucleotide_hybridization, opencloning_models, parsers, primer, primer_screen, readers, recombinase, seq, seqrecord, sequence_picker, sequence_regex, snapgene_history_parser, threading_timer_decorator_exit, tm, types, utils`

Assembly back-ends implemented in `assembly2.py` (function names read from source): `gibson_assembly, in_fusion_assembly, fusion_pcr_assembly, in_vivo_assembly, restriction_ligation_assembly, golden_gate_assembly, ligation_assembly, gateway_assembly, homologous_recombination_integration/excision/inversion, cre_lox_integration/excision, recombinase_assembly`, plus `PCRAssembly` and `SingleFragmentAssembly`.

---

## 3. Capability-by-capability assessment

### Block A — library specification

**`library_as_object` = no.**
There is no library/pool/collection class anywhere in `src/pydna/`. The unit of work is a single `Dseqrecord` (one dsDNA molecule). Functions that can return more than one molecule return plain Python lists: `Assembly.assemble_linear(...) -> list[Dseqrecord]`, `assemble_circular(...) -> list[Dseqrecord]`, `Dseqrecord.cut(*enzymes) -> Tuple[Dseqrecord, ...]`. The multiple products are *alternative outcomes of one assembly reaction*, not a designed library. The paper is explicit that the goal is "the sequence of the final molecule(s)". To build many variants the user writes their own loop. Repo-wide code search for `barcode` → 0 hits; for `mutagenesis` → 0 hits.

**`dag_chaining` = partial.**
Steps compose naturally because each operation returns a `Dseqrecord` that is the input to the next (`pcr(...)` → `.cut(...)` → `.looped()` → `gibson_assembly([...])`), and since v6-dev pydna *records* the resulting dependency graph: every derived `Dseqrecord` carries a `source` object pointing at its inputs, and `Dseqrecord.history()` prints the tree. Docstring of `Dseqrecord.history()` (verbatim from `src/pydna/dseqrecord.py`):

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

The 2015 paper also showed a script-level dependency graph (Figure 5: "A dependency graph produced from the Lactose pathway pydna source files"), but that graph was produced by an *external* software-dependency tool, not by pydna.
Why only *partial*: the DAG is a **record of eagerly executed steps on one molecule lineage**, not a declarative pipeline object that is constructed first and evaluated later, and it cannot be parameterised to fan out over a set of variants.

**`lazy_evaluation` = no.**
Execution is eager: "Execution of the script automatically yields the sequence of the final molecule(s) and that of any intermediate constructs" (Abstract). `Assembly.assemble_linear/assemble_circular/assemble_insertion` return `list[Dseqrecord]`; all products are materialised. Internally the *graph path search* is generator-based and capped (`utils.limit_iterator`, e.g. `limit_iterator(nx.cycles.simple_cycles(self.G), 10000)` in `assembly2.py`), but that is a combinatorial-explosion guard on path enumeration, not lazy sequence generation exposed to the user.

**`mixed_mutagenesis_one_pool` = no.**
pydna has no mutagenesis operation at all. GitHub code search `mutagenesis repo:pydna-group/pydna` → `total_count: 0`. There is no substitution/deletion/insertion scanning, no exhaustive-singles, no random sampling of variants, no WT-replicate concept. The docs' "Bioengineering/SynBio projects" section links out to a *third-party* notebook ("02 Primer directed mutagenesis with pydna", in the `teemi` repo) — i.e. mutagenesis is done by the user on top of pydna, one edit at a time.

**`combinatorial_motif_place` = partial.**
pydna's `Assembly` builds a directed graph over fragments where "The sign of the node key represents the orientation of the fragment, positive for forward orientation, negative for reverse orientation" (`Assembly` docstring) and then enumerates **all** linear/circular/insertion paths through it. From the `new_assemblies` notebook: "there are two possible circular assemblies, one with all fragments oriented as they were passed to the class constructor, and one with the second fragment inverted." So a set of parts that share homology does get combined in all valid orders and orientations, which is a genuine combinatorial enumeration. However: (i) it is driven purely by sequence homology/cut-site compatibility, not by a user-specified placement rule; (ii) there is no API to say "place motif M at positions p1…pn, in both orientations, in all subsets"; (iii) the enumeration is a *simulation of possible reaction products* (used to catch strategy errors), not a design specification. The docs point users to the external `teemi` package for actual combinatorial library work (see §4).

**`barcode_generation` = no.**
Repo-wide code search for `barcode` → 0 hits. No GC/edit-distance-constrained oligo generation anywhere. `utils.randomDNA/randomRNA/randomORF/randomprot` generate unconstrained random sequences (no GC or distance constraints, no attachment machinery).

**`per_sequence_provenance` = yes.**
From `src/pydna/opencloning_models.py` module docstring (verbatim):

> "When using pydna to plan cloning, it stores the provenance of ``Dseqrecord`` objects in their ``source`` attribute. Not all methods generate sources so far, so refer to the documentation notebooks for examples on how to use this feature. The ``history`` method of ``Dseqrecord`` objects can be used to get a string representation of the provenance of the sequence. You can also use the ``CloningStrategy`` class to create a JSON representation of the cloning strategy. That ``CloningStrategy`` can be loaded in the OpenCloning web interface to see a representation of the cloning strategy."

`Dseqrecord.source` is a typed pydantic model. Source classes present in `opencloning_models.py` include `PCRSource, SequenceCutSource, RestrictionAndLigationSource, GibsonAssemblySource, InFusionSource, OverlapExtensionPCRLigationSource, InVivoAssemblySource, LigationSource, GatewaySource, HomologousRecombinationSource, CRISPRSource, CreLoxRecombinationSource, RecombinaseSource, OligoHybridizationSource, PolymeraseExtensionSource, AnnotationSource, ReverseComplementSource, GenomeCoordinatesSource, AddgeneIdSource, BenchlingUrlSource, SnapGenePlasmidSource, …`, and each carries operation-specific fields (e.g. `RestrictionAndLigationSource.restriction_enzymes: list[AbstractCut]`). There is also `Dseqrecord.validate_history(recursive=True)` which "Validate[s] the cloning history of this sequence by replaying operations."
Caveat for the table: this is provenance of **how one molecule was cloned**, not of **how a designed variant relates to a wild-type/design specification**. It is nevertheless structured, typed, machine-readable per-sequence provenance, so "yes" is the honest value.

**`automatic_naming` = partial.**
Derived records get auto-generated names/ids, but they are operation-tagged rather than design-descriptive, and are truncated to 16 chars (GenBank LOCUS limit). From `src/pydna/amplify.py`:

```python
prd.name = (
    identifier_from_string(new_identifier)[:16]
    or self.kwargs.get("name")
    or f"{len(prd)}bp_PCR_prod"[:16]
)
```

From `src/pydna/dseqrecord.py`: `answer.id = "{name}_lin".format(name=self.name)`, `answer.name = answer.id[:16]` (linearize) and `answer.name = "{}_rc".format(self.name[:13])` (reverse complement). So names like `1011bp_PCR_prod`, `pUC19_lin`, `pUC19_rc` are produced automatically, and PCR products inherit a feature label from the template when one fully spans the amplicon. There is no systematic naming scheme encoding a variant's design (position/mutation/replicate), because there are no variants.

**`design_visualization` = yes.**
Multiple mechanisms, all documented:
- `Dseqrecord.figure()` / `Amplicon.figure()` — ASCII/HTML-marked double-stranded figures and primer-annealing figures (see README example output, and `figure(fig_type=...)`: "If the source is an ``Assembly`` object, it returns a compact figure of the assembly … If the source is a ``PCRSource`` object, it returns a figure of the PCR alignment.").
- `Contig.figure()` from the paper: "The figure method gives an text based figure outlining how the sequences were assembled for rapid inspection."
- `pydna.gel.gel(samples=..., gel_length=600, ...)` renders an actual gel image via PIL/matplotlib (README: "Gel electrophoresis of DNA with generation of gel images").
- `Dseqrecord.history()` tree (above) and `CloningStrategy.model_dump_json()` → loadable in the **OpenCloning** web interface for a graphical view of the cloning strategy.
- Paper: homologies "are added to each sequence as metadata in the form of Genbank features (Figure 1A) which can be inspected graphically using a sequence editor."
Caveat: all of this visualises *one construct / one assembly*, never a library.

### Block B — assay coverage

**`assay_dms` = no.** No coding-variant / codon-scanning / saturation-mutagenesis functionality (see `mixed_mutagenesis_one_pool`). pydna does know about coding sequences (`seqrecord.isorf()`, `translate()`, `orfs()`), and `seq.py` has `cai()`, `rarecodons()`, `startcodon()`, `stopcodon()`, `express()`, but these analyse a single ORF for expression prospects — they do not design a variant library. Searched paper, README, docs gallery, and all module names for DMS/deep mutational scanning: absent.

**`assay_mpra` = no.** Nothing about regulatory element libraries, oligo pools, MPRA/STARR, or promoter/enhancer tiling in pydna itself. The docs' external-examples section links to a *teemi* notebook on "Promoter library design in *Saccharomyces cerevisiae*" — that is an external package, and the yeast promoter-library context is metabolic engineering, not MPRA.

**`assay_insilico` = no.** No machine-learning/model-probing functionality of any kind; no torch/tensorflow/model dependency in `pyproject.toml` (dependencies are biopython, networkx, numpy, prettytable, pydivsufsort, pyfiglet, seguid, regex, opencloning-linkml, sgffp, plus optional pyperclip/requests/pyparsing/cai2/scipy/matplotlib/pillow/pyahocorasick).

### Block C — genomics integration

**`genome_coordinates` = partial.**
`pydna.genbank.Genbank.nucleotide(item, seq_start=None, seq_stop=None, strand=1)` accepts an accession-plus-region string and parses it with regexes; docstring examples include:

```
| BK006936.2 REGION: complement(613900..615202)
| NM_005546 REGION: 1..100
| NM_005546 REGION: complement(1..100)
```

and the code also accepts `accession:start-stop` / `accession:cstart-stop` forms (`re.search(r"(:|\s)(?P<start>\d+)-(?P<stop>\d+)", item)`). `BK006936.2` is a *S. cerevisiae* chromosome record, so chromosome-coordinate retrieval is genuinely supported — but only via NCBI accessions over the network, with no genome-build concept, no local genome FASTA/2bit, and no chromosome-name (`chr7:1000-2000`) interface. The provenance model additionally has `GenomeCoordinatesSource(assembly_accession, locus_tag, gene_id, coordinates)` — but that is a record type for the OpenCloning data model (retrieval by assembly accession is implemented in the OpenCloning backend), not a pydna design entry point.

**`transcript_models` = no.**
`src/pydna/parsers.py` handles GenBank, EMBL, FASTA (`embl_gb_fasta()`) and SnapGene (`parse_snapgene()`); it also has `parse_primers()`. No GFF/GTF parser exists and no module imports one. Annotation comes from GenBank feature tables only. There is no transcript/exon model object, no CDS-splicing logic across features.

**`exon_intron_split_codons` = no.**
Consequence: nothing that assembles a CDS across exons or handles codons split across an exon junction. `seqrecord.isorf()` / `seq.orfs()` operate on a contiguous stretch. Checked all of `seq.py`, `seqrecord.py`, `dseqrecord.py` method lists — no exon/intron handling.

**`hgvs_input` = no.**
Code search `hgvs repo:pydna-group/pydna` → 1 hit, and it is inside a test GenBank fixture (`tests/pUC_LAC4_correct_rotation.gb`), i.e. an incidental substring in sequence/annotation text, not an HGVS parser. No `hgvs` dependency in `pyproject.toml`. No mention in paper or docs.

**`vcf_vep_output` = no.**
Code search `vcf repo:pydna-group/pydna` → `total_count: 0`. Output formats are GenBank/FASTA/EMBL (`Dseqrecord.format(format="gb")`, `.write()`), gel PNGs, and the OpenCloning `CloningStrategy` JSON. No variant-file emission.

**`consequence_annotation` = no.**
No stop-gained/synonymous/frameshift classification anywhere. The closest features are `seqrecord.isorf()` (is this an ORF?) and `translate()`, which the user would have to interpret themselves. No variant-consequence vocabulary in the codebase.

### Block D — physical construction

**`primer_design` = yes**, but only assembly/amplification primers — not mutagenic primers.
`src/pydna/design.py` provides `primer_design()` ("Designs forward and reverse primers for PCR amplification of a template sequence"; can match an existing primer's Tm), `assembly_fragments()` — docstring: *"Adds tails to primers for a linear assembly through homologous recombination or Gibson assembly"* — `circular_assembly_fragments()` (deprecated alias for `assembly_fragments(circular=True)`), and an experimental `user_assembly_design()` for USER cloning. Paper: "Pydna provides the pydna assembly_primers function in order to automatically design tailed primers for a series of DNA fragments … The algorithm tries to create primers with balanced melting temperature for the annealing region." Amplicon objects "store rich information about the PCR simulation, such as the DNA region where the primer anneals, melting temperature of each primer and also a suggested PCR program" (`src/pydna/tm.py` holds the Tm models). There is **no** mutagenic-primer designer (no QuikChange-style overlapping mutagenic primer pairs generated from a desired mutation set).

**`codon_optimization` = no** (with a caveat worth stating).
`src/pydna/codon.py` contains only data — per-organism codon `weights`, `start`/`stop` codon frequency distributions, `rare_codons` lists, and `n_end` half-life data (yeast `sce` is the shipped organism). `src/pydna/utils.py` provides `cai()` ("Calculates codon adaptation index for a sequence relative to an organism") and `rarecodons()` ("Identifies rare codons within a coding sequence"), surfaced as `Seq.cai()`, `Seq.rarecodons()`, `Seq.startcodon()`, `Seq.stopcodon()`, `Seq.express()` (optional `express` extra → `cai2`). These **diagnose** codon usage; there is no function anywhere that back-translates a protein or rewrites a CDS to improve codon usage. `utils.express` is documented as "Not yet implemented".

**`synthesis_constraints` = no.**
No synthesis-feasibility checking: no GC-window scan, no homopolymer/repeat detection, no hairpin/secondary-structure or repeat-content filter, no oligo-length/pool constraints. `Seq.gc()` returns a single whole-sequence GC fraction and nothing enforces anything with it. `primer_screen.py` (Aho-Corasick) screens a *primer list* for annealing positions — a diagnostic PCR utility, not a synthesis constraint checker. `sequence_regex.py` compiles IUPAC-degenerate motifs into regexes, used for recognising recombinase/enzyme sites, not for constraint filtering.

### Block E — engineering

**`interface`**: Python API only (importable library, `import pydna`), intended for scripts and Jupyter notebooks. `pyproject.toml` declares **no** `[project.scripts]`/console entry points → no CLI. No GUI and no web service of its own: the paper's "pydna live" web console (ref. 9, `http://pydna-shell.appspot.com/`) is **dead (HTTP 404)**. A separate project, **OpenCloning** (https://opencloning.org, github.com/OpenCloning/OpenCloning), provides a web GUI built on pydna and consumes pydna's `CloningStrategy` JSON — but that is a different tool.

**`license`**: BSD 3-Clause. `LICENSE.txt` header: "Copyright (c) 2013-2026 Björn Johansson, Centro de Biologia Molecular e Ambiental (CBMA)/Department of Biology, University of Minho, Braga, Portugal" with the three standard clauses including "Neither the name of the organizations … may be used to endorse or promote products". Every source file carries `SPDX-License-Identifier: BSD-3-Clause`. PyPI classifier: "License :: OSI Approved :: BSD License". (The 2015 paper stated "License: FreeBSD"; the shipped licence text is 3-clause. GitHub's licence detector reports `NOASSERTION` because the file is named `LICENSE.txt` with a custom preamble.)

**`maintained`**: Yes, actively. Latest commit on `master`: **2026-08-01** (`4e02f81d "bandwith argument (#671)"`); repository `pushed_at` **2026-08-10T08:46:14Z**. Latest PyPI release **5.5.16 on 2026-06-09** (releases v5.5.12–v5.5.16 all in June 2026); `pyproject.toml` on master shows an in-development `6.0.0-a.24`. 228 stars, 67 open issues, not archived, CI (tests + coverage + docs publish) green badges. Two active maintainers (Björn F. Johansson, Manuel Lera-Ramirez). Python 3.10–3.14 supported.

---

## 4. pydna's own documented examples / vignettes

Getting-started tutorial notebooks (`docs/notebooks/`, rendered at https://pydna-group.github.io/pydna):
- `Dseq.ipynb` — "Representing sequences in pydna"
- `Dseq_Features.ipynb` — "Working with Features using the Dseqrecord class"
- `Importing_Seqs.ipynb` — "Importing and viewing sequence files in pydna"
- `Restrict_Ligate_Cloning.ipynb` — "Restriction and Ligation"
- `PCR.ipynb` — PCR simulation
- `primer_design.ipynb` — Primer design
- `Gibson.ipynb` — Gibson assembly
- `CRISPR.ipynb` — CRISPR-Cas9
- `history.ipynb` — "Cloning history / Cloning strategy" (provenance, `history()`, `CloningStrategy` JSON → OpenCloning)
- `new_assemblies.ipynb` — how the `assembly2` engine represents Gibson/restriction-ligation/Golden Gate/PCR/HR assemblies
- `readme_example.ipynb`, `hackathon_copenhagen_2025.ipynb`

Example gallery (`docs/example_gallery.rst`, verbatim bullets):
- **Example_Restriction** — "PCRing a gene out of the genome, and cloning into a vector using restriction and ligation."
- **Example_Gibson** — "Gibson assembly of *R. cellulolyticum* genomic fragments into a plasmid, from the original Gibson assembly paper (doi: 10.1038/nmeth.1318)."
- **Example_CRISPR** — "Using CRISPR with homologous recombination to delete genes by making two cuts in the genome, and repair it with an oligo. Used in the industrially relevant *K. phaffi*."
- **Example_primer_screen** — "Using the Aho-Corasick algorithm to quickly screen a primer list for annealing positions in a sequence."

External examples the docs point to (all require the third-party `teemi` package — **not** pydna functionality):
- "Example_HT_cazyme_primer_design: Design primers for a high-throughput CAZyme library."
- "Example_designing_primers_for_kozak_library: We explore the combinatorial space of the most abundant kozak sequences and make repair-primers for the experiments."
- "02 Primer directed mutagenesis with pydna"
- "03 Investigate plastic-degrading enzymes with pydna"
- "04 Promoter library design in *Saccharomyces cerevisiae*"

From the 2015 paper (files still shipped in `docs/cookbook/`): construction of **YEp24PGK_XK** by BamHI/BglII digestion + ligation (12 lines of code, 11452 bp product); construction of **pGUP1** by homologous recombination (11 lines, 9981 bp product); assembly of a two-gene **lactose pathway** (LAC4 + LAC12 with PDC1/PGI1/TPI1 promoters) described by "nineteen short pydna scripts", experimentally validated in *S. cerevisiae*.

**Reproducibility note for PoolParty:** none of pydna's own examples is a variant-library design. Attempting to reproduce them in PoolParty would be a category error; the honest framing is that pydna's use cases are cloning-strategy simulation/documentation and are orthogonal to PoolParty's.

---

## 5. Notable capabilities NOT covered by the current row list

Candidate extra rows (pydna would score "yes" on all of these; most other surveyed tools would score "no"):

1. **Molecular-biology reaction simulation** — restriction digestion with correct sticky/blunt end bookkeeping, ligation compatibility, PCR (incl. tailed primers and inverse PCR on circular templates), Gibson, In-Fusion, Golden Gate, Gateway (BP/LR), Cre-lox, generic recombinases, fusion PCR, oligo hybridization, polymerase extension, CRISPR cut + HR repair.
2. **Double-stranded / circular topology as a first-class data model** — `Dseq` tracks watson/crick strands, 5′/3′ overhangs, and circularity (`Dseqrecord(..., circular=True)`), with `looped()`, `synced()`, `shifted()`, `fill_in()`, `T4()`, `nibble_*()`, `terminal_transferase()`.
3. **Structured provenance / cloning-strategy serialisation** — typed `source` models + `history()` tree + `CloningStrategy` JSON interoperable with the OpenCloning web app + `validate_history()` replay verification. (This is the row that overlaps PoolParty's "design cards" most directly, and it is the one place where the prior analysis is out of date.)
4. **Round-trip with lab file formats and repositories** — GenBank/EMBL/FASTA/SnapGene parsing, GenBank download by accession+region, `genbankfixer` for malformed records, Addgene/Benchling/SEVA/iGEM/Euroscarf source records, `seguid` checksums, ApE/SnapGene-friendly feature colouring, clipboard copy.
5. **Virtual gel electrophoresis** — `pydna.gel.gel()` renders a gel image with a molecular-weight ladder.
6. **CRISPR protospacer/PAM search** — `pydna.crispr` with a `_cas` abstract base and enzyme classes.
7. **Fast primer screening over many sequences** — `primer_screen.py`, Aho-Corasick automaton, savable/loadable.
8. **Codon-usage diagnostics** — CAI, rare-codon location, start/stop codon frequency, N-end-rule half-life (`Seq.express()`), yeast tables shipped.

---

## 6. Availability today

- **Installable:** yes. `pip install pydna` → 5.5.16 (2026-06-09) on PyPI; Python 3.10–3.14; pure Python with optional extras `clipboard, download, express, gel, primer_screen`.
- **Repo alive:** yes. github.com/pydna-group/pydna, last commit 2026-08-01, pushed 2026-08-10, not archived, CI green, 228 stars, 67 open issues, active dual maintainership.
- **Docs alive:** yes, https://pydna-group.github.io/pydna (200). Legacy readthedocs still resolves.
- **Dead links from the paper:** the "pydna live" interactive web console (http://pydna-shell.appspot.com/) is gone (404); the Google Code source link is obsolete; the Dropbox-hosted cookbook nbviewer link is obsolete. The old GitHub URL redirects correctly.

## 7. Confidence

Highest confidence: the Block B/C "no"s (verified by repo-wide code search and by reading every module name plus the parser and utils APIs), `primer_design`, `license`, `maintained`, `interface`, `barcode_generation`, `mixed_mutagenesis_one_pool`.
Least certain judgement calls: `dag_chaining` and `combinatorial_motif_place` (both scored "partial" — the underlying machinery genuinely enumerates combinations and records a DAG, but neither is a library-design abstraction, and a pydna author could reasonably argue for "yes" on either); `per_sequence_provenance` ("yes", but it is cloning provenance rather than design-variant provenance); `codon_optimization` ("no", though CAI/rare-codon analysis exists and someone could argue "partial"); `genome_coordinates` ("partial" hinges on accepting NCBI accession+REGION as genome coordinates).
Not verified by execution: nothing was installed or run, per instructions — all statements come from the paper, the rendered docs, and the `master` source read via raw.githubusercontent.com on 2026-08-11.
