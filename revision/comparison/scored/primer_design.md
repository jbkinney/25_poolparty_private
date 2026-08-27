# Primer design (`primer_design`, row 20)

Audit date: 2026-08-17. Each repository was inspected at the immutable commit
named below. The Mutation Maker repository linked by its paper was unavailable,
so its own open-access methods paper is the primary evidence for that tool.

## Binding operational test

Apply the global rule literally: credit only a tool-provided operation, parameter,
or mode, not a workflow reconstructed by the user from generic sequence
primitives. This row asks whether the tool designs the PCR or mutagenic primers
used to construct a library in the laboratory (including their sequences, melting
temperatures, overlaps, or mixing ratios).

- `yes`: primer sequences are a documented tool output.
- `partial`: the tool provides primer-specific handling, such as checking a
  proposed primer-binding segment against supplied primers for
  heterodimerisation, but does not itself expose construction primers as a
  documented output.
- `no`: neither primer output nor primer-specific handling is provided.
- `unknown`: the public primary evidence genuinely cannot determine the result.

This test does **not** credit an oligo merely because it is synthesised, an adaptor
merely because it could be used as a priming site, or a user-supplied primer that
the tool only carries into an assembled construct. It also does not credit generic
sequence design merely because a user labels the resulting sequence a primer.
Checking/filtering physical defects in a designed library sequence belongs to row
19; repairing such a sequence belongs to row 11.

## Results

| Tool | Value | Confidence | Operational finding |
|---|---|---:|---|
| PoolParty | `no` | high | Its documented operation catalogue and public API contain no primer output or primer-specific operation. |
| VaLiAnT | `no` | high | It emits variant oligos and metadata and accepts optional adaptors; it neither designs nor handles primers. |
| MPRAnator | `no` | high | Its three design surfaces generate motif-placement, SNP, and transmutation sequences, with no primer API or output. |
| MPRA Design Tools (`mpradesigntools` / `designMPRA`) | `no` | high | Forward and reverse PCR primers are user inputs that are concatenated into constructs, not tool-designed outputs. |
| Oligopool Calculator (`oligopool`) | `yes` | high | Its `primer` operation designs forward/reverse, paired, and per-set primers and returns their sequences in a named output column plus thermodynamic statistics. |
| Mutation Maker | `yes` | high | Its documented SSSM, MSDM, and PAS workflows output mutagenic oligos/primers and associated Tm/structure information; PAS also calculates oligo mixing ratios. |
| DNA Chisel | `partial` | high | It provides primer-specific Tm, homology, repeat, and heterodimer constraints, but the documented “primer collection” is user orchestration of generic random-sequence optimization rather than a primer-design output API. |
| tangermeme | `no` | high | Its complete source/API surface concerns model interrogation and model-guided sequence edits and has no primer output or primer-specific handling. |

## Per-tool evidence and disconfirmation

### 1. PoolParty — `no` (high confidence)

**Evidence.** At commit
[`1bb0179e1c3720b1fffd471802b3040f9336de28`](https://github.com/jbkinney/poolparty-statetracker/tree/1bb0179e1c3720b1fffd471802b3040f9336de28),
the documentation calls its operations index the catalogue of “over 50 built-in
operations” and partitions it into sources, mutagenesis, scanning, regions, ORF,
utilities, composition, and state operations
([operations index, lines 1–20 and 43–75](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/docs/operations/index.rst#L1-L75)).
The exported public surface enumerates the operations and utilities; it ends with
sequence-property and restriction-site helpers and contains no primer operation
([`__all__`, lines 146–262](https://github.com/jbkinney/poolparty-statetracker/blob/1bb0179e1c3720b1fffd471802b3040f9336de28/poolparty/src/poolparty/__init__.py#L146-L262)).

**Reasoning.** A designed `Pool` can be materialised as oligos, but no documented
operation emits primers or applies primer-specific handling. Generic prefixing,
barcode generation, restriction-site tests, or filtering are not primer design
under the global rule.

**Disconfirmation attempt.** I searched the complete tracked package source and
documentation (excluding generated images/notebooks) for `primer`, `PCR`, `Tm`,
`melting`, `hairpin`, `homodimer`, `heterodimer`, and `mixing ratio`, then checked
the complete operations index and public exports. The sole relevant-looking source
hit was a comment describing restriction enzymes to avoid in Gibson *primers*; it
does not define an operation. A documented primer sequence/output or a
primer-specific constraint would have changed this score.

### 2. VaLiAnT — `no` (high confidence)

**Evidence.** At commit
[`8796cc112dafd4919fec59913f58cd2be87c45eb`](https://github.com/cancerit/VaLiAnT/tree/8796cc112dafd4919fec59913f58cd2be87c45eb),
VaLiAnT's exhaustive usage summary lists optional 5′/3′ adaptors as inputs and
lists its outputs as a retrieval-QC file, oligonucleotide metadata, variants,
unique oligonucleotides, and configuration—no primer output
([README lines 61–102](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/README.md#L61-L102)).
The shared CLI likewise describes `adaptor-5` and `adaptor-3` as sequences “to be
added” to the oligo
([README command reference](https://github.com/cancerit/VaLiAnT/blob/8796cc112dafd4919fec59913f58cd2be87c45eb/README.md#L196-L218)).

**Reasoning.** Supplied adaptors and variant-library oligos are not designed PCR
or mutagenic primers. There is no primer-specific operation or output.

**Disconfirmation attempt.** I searched every tracked source, test, example, and
README/wiki-facing file at the commit using the same primer/Tm/PCR/structure terms.
The only `primer` hits were `primer_bind` annotations already present in an example
GenBank plasmid—input annotations, not generated primers. I also inspected the CLI
inputs and complete output-file list. A generated forward/reverse or mutagenic
primer field, or a named primer constraint, would have changed the score.

### 3. MPRAnator — `no` (high confidence)

**Evidence.** The tool's own paper describes exactly three components: MPRA Motif
Design, MPRA SNP Design, and Transmutation (scrambling, reversing,
complementing, or random mutation)
([MPRAnator paper, Abstract/Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC5198521/)).
The corresponding public source snapshot is commit
[`9969790d62410138d4281b7955da6d085f07b1bc`](https://github.com/hemberg-lab/MPRAnator/tree/9969790d62410138d4281b7955da6d085f07b1bc),
whose downloadable interfaces are `MpraMotifs_script.py` and
`MpraSnps_script.py`; the remaining sequence transformation surface is implemented
in `part3.py`.

**Reasoning.** Those operations generate experimental insert sequences or
controls. None emits construction primers or provides primer-specific handling.

**Disconfirmation attempt.** I searched the full repository—application forms,
views/tasks, downloadable scripts, core modules, and tests—for `primer`, `PCR`,
`Tm`, `melting`, `hairpin`, `homodimer`, `heterodimer`, and mixing/ratio terms;
there were no primer-related hits. I then checked the paper's feature and output
descriptions for a wet-lab primer stage. A documented primer sequence/output or
primer-specific parameter would have changed the score.

### 4. MPRA Design Tools — `no` (high confidence)

**Evidence.** In `mpradesigntools` commit
[`afd386ef12051bb0a37ad63a6f92acd555246584`](https://github.com/andrewGhazi/mpradesigntools/tree/afd386ef12051bb0a37ad63a6f92acd555246584),
the documented example supplies literal `fwprimer` and `revprimer` strings to
`processVCF`
([README lines 106–120](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/README.md#L106-L120)).
The function documentation is unambiguous: both arguments are strings containing
the PCR primers “to be used”
([`processVCFfast.R` lines 193–218](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L193-L218)).
The implementation simply concatenates those supplied strings around context,
restriction sites, and barcodes
([lines 398–418](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L398-L418)).
The companion Shiny app at commit
[`0cf56ee602fc86dde705906d071a39cbdf6e99a8`](https://github.com/andrewGhazi/designMPRA/tree/0cf56ee602fc86dde705906d071a39cbdf6e99a8)
likewise presents the forward and reverse primers as text inputs
([`ui.R` lines 50–55](https://github.com/andrewGhazi/designMPRA/blob/0cf56ee602fc86dde705906d071a39cbdf6e99a8/ui.R#L50-L55))
and passes them unchanged into `processVCF`
([`server.R` lines 113–123](https://github.com/andrewGhazi/designMPRA/blob/0cf56ee602fc86dde705906d071a39cbdf6e99a8/server.R#L113-L123)).

**Reasoning.** Carrying user-authored primers into the final construct is not
designing or emitting primers. No Tm, overlap, primer-pair, dimer, or other
primer-specific mechanism is applied, so even `partial` is unwarranted.

**Disconfirmation attempt.** I searched all R source, generated Rd pages, README,
Shiny UI/server code, and scripts in both repositories. Every meaningful primer hit
resolved to an input parameter or literal default followed by concatenation. I
also searched the thermodynamic/structure terms and found no primer-design
routine. A tool-selected primer sequence or primer-specific validation/constraint
would have changed the score.

### 5. Oligopool Calculator — `yes` (high confidence)

**Evidence.** At tagged commit
[`b88fa394ca67ed4c48ec41127b5707694ee7cf0a` (`v2026.02.22.1`)](https://github.com/ayaanhossain/oligopool/tree/b88fa394ca67ed4c48ec41127b5707694ee7cf0a),
the first-party API says `primer` designs amplification primers with Tm,
hairpin/homodimer/heterodimer/cross-dimer, repeat, and optional background
screening; its signature includes direction, Tm range, pairing, per-set grouping,
and a caller-named `primer_column`
([API lines 145–200](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/api.md#L145-L200)).
The same reference states that paired primers are chained for matched Tm and that
per-set primers are designed with cross-set dimer screening
([lines 202–209](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/api.md#L202-L209)).
The implementation's return contract is “A pandas DataFrame of designed primers”
([`primer.py` lines 17–93](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/primer.py#L17-L93));
on success it inserts the generated primer into `primer_column`, optionally writes
CSV, and reports Tm and structural metrics
([lines 1280–1353](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/primer.py#L1280-L1353)).

**Reasoning.** The tool directly emits primer sequences and primer thermodynamic
properties, including paired and multiplexed designs. This meets `yes`; it is not
merely synthesis-constraint checking because the checks operate inside a named
primer-generating API and the generated primer is returned.

**Disconfirmation attempt.** I tested the contrary hypothesis that `primer` only
checks a caller-supplied sequence. Its required `primer_sequence_constraint` is an
IUPAC design space, not a completed primer; source invokes `primer_engine`, then
inserts the selected `primer`/per-set primers into the output column. I inspected
the API, implementation, CLI/example pipeline, paired-primer path, and success
output. The score would fall to `partial` if the API only screened an existing
primer or returned no selected sequence, but both are directly disproved.

### 6. Mutation Maker — `yes` (high confidence)

**Evidence.** Mutation Maker's own open-access methods paper describes it as a
“mutagenic oligo design software” for SSSM, multisite directed mutagenesis, and
PCR-based accurate synthesis, and says the workflows are designed around
overlap-extension PCR, QuikChange multisite mutagenesis, and PCR-based gene
synthesis
([Hiraga et al. 2021, Introduction/Methods](https://pubs.acs.org/doi/10.1021/acssynbio.0c00542)).
The same paper's “Melting Temperature Calculations, Primer Dimer and Hairpin
Checks” section documents Primer3-based Tm, hairpin, homodimer, and heterodimer
calculations; it explicitly compares SSSM flanking primers with candidate
mutagenic forward/reverse primers and all MSDM primers. It further states in the
author-contribution record that PAS oligo generation and ratio calculation were
implemented. These are direct descriptions of the tool's outputs and scoring
logic, not a downstream use invented by a user.

**Reasoning.** Candidate mutagenic forward/reverse primers and PAS assembly oligos
are designed outputs used to construct the library, with documented Tm/structure
properties and, for PAS, ratios. That meets `yes`.

**Disconfirmation attempt.** I checked whether “oligo” meant only a final library
member or whether primers were merely supplied. The methods distinguish flanking
primers from *candidate mutagenic forward or reverse primers*, describe their
selection penalties, and describe multiple PCR construction workflows. I searched
the paper for output, primer, Tm, dimer, hairpin, overlap, and ratio evidence. The
linked GitHub repository was unavailable at audit time, so no unverifiable source
claim is used. The score would fall if the paper only described checking supplied
primers, but its candidate-design and workflow language directly rules that out.

### 7. DNA Chisel — `partial` (high confidence)

**Evidence.** At version 3.2.16, commit
[`68c09304341c3656f3dfe63eda37757d6a7b3917`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/tree/68c09304341c3656f3dfe63eda37757d6a7b3917),
the built-in `AllowPrimer` specification makes a sequence location
“primer-friendly” for PCR/Sanger sequencing and groups named melting-temperature,
homology/repeat, and heterodimer constraints
([`AllowPrimer.py` lines 1–24 and 29–52](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/builtin_specifications/AllowPrimer.py#L1-L52)).
Its implementation registers `EnforceMeltingTemperature` and optional
`AvoidHeterodimerization` against supplied external primer sequences
([lines 56–90](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/builtin_specifications/AllowPrimer.py#L56-L90)).

The documented “collection of compatible primers” example also supplies the
boundary counterexample. It says DNA Chisel is “not originally meant for creating
collections,” defines its own `create_new_primer`, starts from a generic random
20-mer, applies generic constraints, returns `problem.sequence`, and performs its
own loop/list bookkeeping
([example lines 1–46](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/examples/common_scenarios/primers_collection.py#L1-L46)).

**Reasoning.** `AllowPrimer` and `AvoidHeterodimerization` are genuine,
tool-provided primer-adjacent handling, so `no` would be too low. But they make a
subsequence compatible with priming; they do not choose PCR amplicon endpoints,
emit forward/reverse construction primers, or expose a primer-design result type.
The example's user-defined loop relabels generic optimized random sequences as
primers, exactly the reconstruction excluded by the global rule and by the
row-19/20 boundary note. Therefore `partial`, not `yes`.

**Disconfirmation attempt.** I inspected the full built-in specification catalogue,
`AllowPrimer`, `AvoidHeterodimerization`, their tests, the GenBank annotation API,
and the primer-collection example. The apparent positive was challenged by tracing
where the “primer” sequence originates and who performs collection construction:
the user script does. Conversely, the dedicated `AllowPrimer` API disproves `no`.
A documented operation returning construction primer sequences (rather than an
optimized parent sequence/location) would have raised the score to `yes`.

### 8. tangermeme — `no` (high confidence)

**Evidence.** At commit
[`2006b310cd72a28c56c3ea4ba67f738fff74bb89`](https://github.com/jmschrei/tangermeme/tree/2006b310cd72a28c56c3ea4ba67f738fff74bb89),
the project defines its scope as atomic sequence operations, predictive-model
application/interpretation, and model-guided sequence design
([README lines 7–12](https://github.com/jmschrei/tangermeme/blob/2006b310cd72a28c56c3ea4ba67f738fff74bb89/README.md#L7-L12)).
Its usage explains that “design” consumes arbitrary PyTorch models and returns
tensor-based sequence designs
([README lines 80–86](https://github.com/jmschrei/tangermeme/blob/2006b310cd72a28c56c3ea4ba67f738fff74bb89/README.md#L80-L86));
the first-party design API likewise returns ranked designed sequence tensors, not
laboratory primers
([design API](https://tangermeme.readthedocs.io/en/latest/api/design.html)).

**Reasoning.** Model-guided sequence edits are biological sequence design, not
physical construction-primer design. No primer output or primer-specific handling
is provided.

**Disconfirmation attempt.** I searched every Python module, API `.rst`, test,
README, and notebook source text for the full primer/PCR/Tm/melting/structure term
set and inspected the complete first-party API module list (ersatz, prediction,
attribution, marginalization/ablation, spacing, saturation mutagenesis, variant
effect, product, and design). There were no primer-specific hits. A function that
returned PCR/mutagenic primers or applied a named primer constraint would have
changed the score.

## Final consistency and boundary check

- **One threshold across tools:** only Oligopool Calculator and Mutation Maker
  have primary evidence that the tool emits construction-primer sequences.
- **Input is not output:** MPRA Design Tools' `fwprimer`/`revprimer` arguments and
  VaLiAnT's adaptors are user-supplied sequence components, so they receive no
  credit.
- **Primer-adjacent is only partial:** DNA Chisel receives `partial` because its
  named APIs handle Tm and primer heterodimerisation, but its “primer collection”
  depends on user-defined generic-sequence orchestration. Avoiding
  heterodimerisation while designing a sequence is not itself primer design.
- **No row-19 double credit:** Oligopool Calculator's physical checks are discussed
  here only because they are internal to a primer-generating operation; the score
  is earned by the returned primer sequence. DNA Chisel's generic synthesis
  constraints are not re-labelled as primer output. Checks on final library oligos
  are outside this row.
- **No row-11 double credit:** sequence repair/optimization alone is insufficient.
  DNA Chisel's optimisation is capped at `partial` here because only its dedicated
  primer-adjacent specifications count; PoolParty/tangermeme sequence operations
  do not become primer design through user composition.
- **No `unknown` hedge:** every tool had enough primary evidence for a determinate
  score. Mutation Maker's deleted/unavailable code repository reduces source-level
  reproducibility but not determinacy because the tool's own permanent,
  open-access methods paper directly documents primer outputs and calculations.
