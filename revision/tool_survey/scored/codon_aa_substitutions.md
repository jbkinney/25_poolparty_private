# Audit: `codon_aa_substitutions`

Date assessed: 2026-08-15

## Measurement instrument and operational test

I read the preamble and the binding global rule in `revision/tables/ROW_DEFINITIONS.md:1-34`, and the row definition at lines 117-124, before examining any tool. In particular, the global rule says that a capability counts only when the tool supplies an operation, parameter, or mode, that `partial` cannot be used for user-assembled logic, and that an example-only capability is at most `partial`. The row requires a codon- or amino-acid-level substitution operation, tool-handled codon-to-residue mapping, and a choice of replacement policy; its special rule says that a user-supplied residue list is `no`, while a fixed tool-provided codon-aware substitution set is `partial`.

**Operational test stated before tool inspection:** For each tool, I checked whether one documented operation lets the user specify substitutions at the codon or amino-acid level, has the tool perform codon↔residue mapping, and offers a tool-defined replacement-policy choice; I scored `yes` if all hold, `partial` if codon-aware substitution exists but its substitution set is fixed, and `no` if substitutions are only nucleotide-level or require the user to supply residue targets/lists or assemble the behavior themselves.

This test does not require the operation to enumerate a library: the row itself does not impose that requirement, and names "most-frequent codon" as an example policy. Thus a documented codon-optimization operation can pass if it maps residues to synonymous codons and offers multiple replacement policies. Conversely, a script that merely loops over triplet-sized nucleotide strings without a genetic-code mapping does not become codon-aware in the row's sense.

## Results

| Tool | Score | Confidence |
|---|---|---|
| poolparty | **yes** | high |
| VaLiAnT | **yes** | high |
| MPRAnator | **no** | high |
| MPRADesignTools / designMPRA | **no** | high |
| Oligopool Calculator | **no** | high |
| Mutation Maker | **partial** | high |
| DnaChisel | **yes** | high |
| tangermeme | **no** | high |

## Per-tool audits

### poolparty — **yes** (high confidence)

**Positive evidence.** The documented operation is explicitly codon-level: `poolparty-statecounter/poolparty/docs/operations/mutagenize_orf.rst:4-7` says, "Introduce codon-level mutations into an ORF sequence" and documents restriction by eligible codons. The `mutation_type` parameter offers seven policies: `any_codon`, `nonsynonymous_first`, `nonsynonymous_random`, `missense_only_first`, `missense_only_random`, `synonymous`, and `nonsense` (`mutagenize_orf.rst:41-47`). The output cards include `wt_codons`, `mut_codons`, `wt_aas`, and `mut_aas` (`mutagenize_orf.rst:88-92`).

The mapping and exact policy semantics are implemented in `poolparty-statecounter/poolparty/src/poolparty/codon_table.py`. The class exposes both "amino acid -> list of codons" and "codon -> amino acid" maps (`:58-68`); builds the inverse codon-to-AA table (`:108-126`); and defines, among other policies, "First codon of each different amino acid including stop," all codons for different amino acids, missense-only equivalents, synonymous codons, and stop codons (`:148-155`, implemented at `:162-206`). The built-in codon lists are frequency ordered, so the `_first` policies choose the most frequent codon (`:26-32`). This directly satisfies codon-level specification, mapping, and a tool-defined policy choice.

**Behavioral check.** As required for this local tool, I ran the read-only venv with bytecode disabled:

```text
PYTHONDONTWRITEBYTECODE=1 poolparty-statecounter/.venv/bin/python
```

After adding the repository `src` directory to `sys.path` and calling `pp.init()`, `mutagenize_orf('ATG', num_mutations=1, mutation_type='missense_only_first', mode='sequential', cards=['mut_codons','mut_aas']).generate_library()` returned 19 variants: one representative codon for each non-methionine, non-stop amino acid. A second run on `GCT` with `mutation_type='synonymous'`, random mode, and five states returned only alanine-encoding codons (the sampled rows contained GCA/GCG). A preliminary invocation using the nonexistent `pp.init(seed=1)` argument failed with a signature error; I corrected it to the documented `pp.init()` before making the behavioral observations. No repository file was written.

**Disconfirmation attempt / what would change the score.** `partial` would have resulted if only one fixed codon-aware set existed; `no` would have resulted if the API only accepted nucleotide alternatives or a caller-authored residue list. I searched the operation documentation and API source for `mutation_type`, `codon`, `amino`, `wt_aas`, `mut_aas`, `synonymous`, `missense`, and `nonsense`, and executed two policies. The multiple built-in policies and explicit maps disconfirm those alternatives.

### VaLiAnT — **yes** (high confidence)

Primary repository revision: `cancerit/VaLiAnT@8796cc112dafd4919fec59913f58cd2be87c45eb` (`develop`).

**Positive evidence.** The README separates nucleotide mutation from four CDS-only codon-aware modes: alanine (`ala`), stop (`stop`), all-amino-acid (`aa`), and `snvre` substitutions (`README.md:255-268`). Their policies are concrete:

* `ala`: replaces each codon with the highest-ranked alanine codon (`README.md:324-336`).
* `stop`: replaces each codon with the highest-ranked stop codon (`README.md:338-350`).
* `aa`: chooses the highest-ranked codon for every other amino acid, producing 19 sequences per amino-acid codon and 20 per stop (`README.md:352-380`).
* `snvre`: use all synonymous triplets for a synonymous change, the top-ranked synonymous triplet for a missense change, and the top-ranked stop for a nonsense change (`README.md:383-393`).

The source performs the mapping rather than asking the caller for residue lists. `src/valiant/codon_table.py:53-71` builds codon→amino-acid and rank-sorted amino-acid→codon maps, while `:95-117` selects top codons and translates codons. `src/valiant/mutators/codon.py:76-103` implements dedicated alanine, stop, and all-amino-acid mutators; the all-AA mutator translates the reference codon and asks the table for top codons excluding the reference residue and stop. The authors' paper independently describes the same operation: protein-level saturation's `aa` function "exchanges each wild-type codon for the most frequent triplet code of each other amino acid" using a default or user-defined frequency table (`Barbon2022_VaLiAnT_all.pdf`, p. 2).

**Disconfirmation attempt / what would change the score.** `partial` would have resulted if `aa` were the only fixed codon set, and `no` if the only operation were nucleotide `snv` plus post-hoc amino-acid annotation. I searched the README mutation-type and codon-table sections, `src/valiant/mutator_type.py`, all codon mutators, `src/valiant/codon_table.py`, and the paper for `codon`, `amino`, `alanine`, `stop`, `aa`, `snvre`, `translate`, `top`, and `rank`. The four documented built-in modes and implemented bidirectional mappings support `yes`.

### MPRAnator — **no** (high confidence)

Primary repository revision: `hemberg-lab/MPRAnator@9969790d62410138d4281b7955da6d085f07b1bc` (`master`).

**Evidence of the available behavior.** The official documentation says the SNP page designs oligos to study SNPs and lets the user include or exclude SNP combinations (`iliasApp/templates/iliasApp/docs.html:40-43`). Its Transmutation page offers "Mutate random positions in the input sequence," accepts nucleic-acid IUPAC letters, and reports `Mutated_nucleotides` (`docs.html:138-159`). The SNP implementation reads `REF` and `ALT` nucleotide fields and substitutes the supplied alternate nucleotide sequence; the relevant primary file is repository-root `parseSNPs.py` (inputs and fields at `:5-6`, `:43-46`; nucleotide substitution at `:382-392`, and caller-supplied alternate targets at `:431-451`). The companion downloadable client likewise illustrates `REF` and `ALT` nucleotide SNP input (`downloadables/MpraSnps_script.py:12-14`). The paper describes MPRA SNP Design as acting on single or combinations of SNPs and Transmutation as introducing random mutations, not coding substitutions (`Georgakopoulos-Soares2017_MPRAnator.pdf`, p. 1).

**Absence search.** I listed the complete 50-path repository tree, then streamed all 29 relevant `.py`, `.html`, and `.md` files from the official archive and searched case-insensitively for:

```text
codon|amino.?acid|protein|residue|translat|NNK|NNS|degenerate
```

There were zero matches. I also inspected `docs.html`, `parseSNPs.py`, `downloadables/MpraSnps_script.py`, `myfunctions.py`, `oligo.py`, `part1.py`, `part3.py`, `viewsCore.py`, and the paper around every mutation/SNP description.

**Disconfirmation attempt / what would change the score.** A fixed documented codon/AA replacement operation with an internal genetic-code map would change this to `partial`; multiple policies would change it to `yes`. The full-tree term search and focused inspection found neither. Because the negative result is based on an exhaustive relevant-source search plus positive evidence that the exposed alternatives are nucleotide REF/ALT and random nucleotides, this is `no`, not `unknown`.

### MPRADesignTools / designMPRA — **no** (high confidence)

Primary revisions: `andrewGhazi/mpradesigntools@afd386ef12051bb0a37ad63a6f92acd555246584` and `andrewGhazi/designMPRA@0cf56ee602fc86dde705906d071a39cbdf6e99a8` (`master`).

**Evidence of the available behavior.** The package README says its main design function is `processVCF`, calculates output from "ref/alt alleles," and uses only the VCF `CHROM`, `POS`, `REF`, and `ALT` columns (`mpradesigntools/README.md:43-58`). It further says multiple alternatives must be supplied in `ALT` and some alleles must be oriented manually (`:56-61`). In `R/processVCFfast.R:5-15`, the function splits a caller's comma-separated `ALT` values into rows. The implementation classifies a substitution from single A/C/G/T `REF` and `ALT` values and constructs the alternate sequence by replacing the reference letter with `snp$ALT` (the same code is in `designMPRA/scripts/processVCFfast.R:74-76,107-119`). The paper likewise describes oligos containing reference and alternate alleles of variants and says variants are input by the user (`Ghazi2018_MPRADesignTools_all.pdf`, pp. 1 and 4).

**Absence search.** I listed all 72 paths in `mpradesigntools` and all 26 paths in `designMPRA`. I streamed and searched all R source, R Markdown, Markdown, and application source for:

```text
codon|amino.?acid|protein|residue|translat|NNK|NNS|degenerate
```

No coding-aware operation appeared. I separately inspected `README.md`, all three package R files, `designMPRA`'s server/UI/scripts, and both copies of `processVCFfast.R`. The compiled `Rmd/Supplement.html` contains embedded/base64 asset noise, so it was not treated as a meaningful text hit; its corresponding source R Markdown was searched instead.

**Disconfirmation attempt / what would change the score.** A package or Shiny parameter that accepted a coding sequence/codon position and internally generated a fixed residue/codon replacement set would produce `partial`; multiple built-in policies would produce `yes`. I searched the full meaningful source set, README/API surface, and paper for exactly those concepts and found only caller-authored nucleotide alleles. This supports `no`.

### Oligopool Calculator — **no** (high confidence)

Primary repository revision: `ayaanhossain/oligopool@b88fa394ca67ed4c48ec41127b5707694ee7cf0a` (`master`).

**Evidence of the documented API.** The README enumerates design modules for barcodes, primers, motifs/anchors, and spacers, and says Degenerate Mode "compress[es] variant libraries" (`README.md:34-42`). The `compress` API's purpose is to "Compress concrete DNA sequences into IUPAC-degenerate oligos" and its required `input_data` must already contain concrete A/T/G/C sequences (`docs/api.md:770-799`). It promises no invented sequences and returns a 1:1 mapping when inputs are too diverse (`:810-819`). `expand` merely expands a caller-supplied IUPAC sequence into concrete A/T/G/C strings (`:834-872`). Neither interface accepts a CDS, genetic code, residue, or replacement policy.

**Why the codon example does not count.** The only codon generator is an example-side helper, not the documented `oligopool` API. `examples/library-compressor/mutant_generator.py:20-24` defines all 64 triplets by the Cartesian product of DNA bases; `generate_codon_variants` then loops over those strings and splices each into a nucleotide offset (`:82-126`). It contains no genetic-code translation or residue mapping. The demo imports this separate `mutant_generator`, generates all 64 input sequences first, and only then calls `op.compress` (`run_degenerate_demo.py:26-30,54-79`). Thus even the example is triplet-sized nucleotide enumeration, followed by user/example composition, rather than the row's tool-handled codon↔residue substitution capability.

**Absence search.** I listed all 85 repository paths and searched the code/docs for `codon`, `amino`, `protein`, `residue`, `translate`, `genetic code`, `NNK`, and `NNS`. All meaningful codon hits were confined to the library-compressor example files just described; none occurred in a documented package substitution interface. I inspected the README feature/API lists, complete Degenerate Mode API, `compress`, `expand`, the core degenerate implementation, and every codon-bearing example.

**Disconfirmation attempt / what would change the score.** A documented `oligopool.*` call that accepted a parent CDS or amino-acid target and internally produced a fixed codon-aware replacement set would produce `partial`; policy choices would produce `yes`. If the example helper had at least performed a genetic-code mapping, it could have supported `partial` at most under the global example-only rule. The searches confirmed none of these; `no` is therefore warranted.

### Mutation Maker — **partial** (high confidence)

Primary repository revision: `ra100/Mutation_Maker@396c1c0ede7529f3dedf65e821e8c1f20c9a7043` (`master`). The paper's printed `github.com/Merck/Mutation_Maker` location is unavailable; the supplied authorship-bearing mirror is the repository assessed here.

**Positive evidence for the restricted capability.** The authors define SSM as substituting a residue with all possible amino acids (`Hiraga2021_MutationMaker.pdf`, p. 2). For their SSSM algorithm, the paper says the minimal input includes mutation sites with their degenerate codon, with default `NNK`; a degenerate codon may be supplied directly, or generated "from a user-specified lists of amino acids" (`Hiraga2021_MutationMaker.pdf`, p. 3). The current interface is explicitly amino-acid/codon aware: mutation syntax is `[Amino Acid Codon][Location][X]` (`frontend/src/scenes/SSM/components/SSMForm.tsx:217-234`), and the parser turns the one-based amino-acid position into a three-base offset and represents the change as `AminoMutation` (`backend/mutation_maker/mutation.py:28-69`). The form provides a **tool default** `NNK` degenerate codon (`SSMForm.tsx:304-328`). This is a usable fixed, tool-provided codon-aware substitution set without requiring a caller-authored residue list, so it reaches `partial`.

**Why it is not `yes`.** The alternate amino-acid route does not provide a tool-defined replacement-policy choice: it asks the user to decide the replacement set. `SSMForm.tsx:97-106` filters amino acids that the user marked `include` or `avoid` and sends those lists to the degenerate-codon generator. The component tells the user, "Select which amino acids you want to include" and "exclude" (`frontend/src/shared/components/AminoAcidInput/index.tsx:58-84`). The mapping implementation imports an amino-acid-to-codon table and expands precisely the supplied amino list (`frontend/src/scenes/SSM/components/amino.ts:19-21,70-79`). That path is explicitly excluded by the row's "If the user must supply the residue list" rule and cannot upgrade the fixed default to `yes`. Directly typing a different degenerate codon likewise makes the caller specify the replacement set; the only tool-chosen set established by the sources is default NNK.

The repository also contains a separate `MutationSite` abstraction that maps its supplied target amino acids to concrete/minimal codons (`backend/mutation_maker/mutation.py:119-163`), but its constructor derives `new_aminos` from caller-provided mutations. It therefore confirms mapping while also confirming that the non-NNK target set is caller supplied.

**Disconfirmation attempt / what would change the score.** `yes` would require at least two tool-defined policy options chosen without the user listing residues or directly spelling the replacement set (for example all residues versus alanine-only or stop-only). `no` would result if even the default NNK path were absent and every target had to be listed. I searched the full tree and paper for `SSM`, `SSSM`, `degenerate codon`, `NNK`, `amino`, `include`, `avoid`, `codon usage`, `minimal triplets`, and mutation input; inspected the SSM form, AA selector and mapper, backend mutation types, degeneracy code, tests, and the paper's SSSM/MSDM/gene-synthesis sections. The fixed default supports `partial`; the searched alternatives are user-authored sets and do not support `yes`.

### DnaChisel — **yes** (high confidence)

Primary repository revision: `Edinburgh-Genome-Foundry/DnaChisel@68c09304341c3656f3dfe63eda37757d6a7b3917` (`master`).

**Positive evidence.** `CodonOptimize` is a documented codon-level operation with three tool-defined replacement policies. Its docstring says it "Codon-optimize[s] a coding sequence using a user-selected method" (`dnachisel/builtin_specifications/codon_optimization/CodonOptimize.py:6-18`) and defines:

* `use_best_codon`: replace every codon by the most frequent synonymous codon;
* `match_codon_usage`: match the target organism's codon-usage profile; and
* `harmonize_rca`: choose a synonymous codon whose target-organism usage matches the original codon's host-organism usage (`CodonOptimize.py:20-29,47-49`).

The same documentation explicitly requires `EnforceTranslation` to preserve the amino-acid sequence (`:31-32`). `EnforceTranslation` is a `CodonSpecification`; it enforces a protein translation and accepts a named genetic table (`dnachisel/builtin_specifications/EnforceTranslation.py:13-18,28-38,51-54`). `CodonOptimize`'s codon-usage table is keyed by amino acid and gives the relative use of each of its codons (`CodonOptimize.py:59-64`). Together these are explicit tool-handled residue↔codon mapping and multiple replacement policies. The operation is synonymous rather than residue-changing, but the row expressly includes substitution "at the codon or amino-acid level" and names most-frequent codon as a qualifying policy; it does not require a changed residue or multi-sequence enumeration.

**Boundary / non-evidence.** The general `MutationSpace` does not independently earn this score. Its public example requires the caller to provide segment variants (`dnachisel/MutationSpace/MutationSpace.py:9-30`), and its default space is individual nucleotide choices (`:166-190`). Similarly, `EnforceChanges` is explicitly measured in nucleotides (`dnachisel/builtin_specifications/EnforceChanges.py:16-31`). Those would be user reconstruction or nucleotide-only under the row. The score rests only on the documented `CodonOptimize`/`EnforceTranslation` interfaces.

**Disconfirmation attempt / what would change the score.** `partial` would have resulted if DnaChisel exposed only one fixed synonymous replacement rule; `no` if codon awareness existed only as a preservation constraint around nucleotide edits or required the caller to enumerate residue targets. I searched the 323-path tree and official paper/docs for `codon`, `amino acid substitution`, `mutagenesis`, `saturation`, `translation`, `mutation space`, `reverse translate`, and inspected README examples, built-in specification docs/source, all codon-optimization methods, `EnforceTranslation`, `EnforceChanges`, `MutationSpace`, and reverse-translation utilities. The three codon policies and genetic-table mapping meet the operational test.

### tangermeme — **no** (high confidence)

Primary repository revision: `jmschrei/tangermeme@2006b310cd72a28c56c3ea4ba67f738fff74bb89` (`main`).

**Evidence of the available behavior.** Tangermeme describes assumption-free atomic sequence operations extensible to any alphabet (`README.md:10-12`). Its saturation-mutagenesis API associates substitutions with entries of the input sequence tensor and interprets the result as the influence of each "observed nucleotide" (`README.md:224-234`). Its variant interface accepts a tensor of coordinate/channel substitutions (`README.md:238-252`). These are alphabet-character or nucleotide operations; there is no reading-frame or genetic-code input.

**Absence search.** I listed the complete 160-path tree (93 source/docs/test files) and searched all meaningful files for:

```text
codon|amino|residue|ORF|translation|genetic_code|NNK|NNS
```

There were zero relevant hits (generic prose uses of "protein" referred to model/readout context, not a protein alphabet or genetic code). I also inspected `saturation_mutagenesis.py`, `variant_effect.py`, `ersatz.py`, the greedy design/substitution modules, the skill documentation for saturation mutagenesis and design, the README API examples, tutorials indicated by the repository, and `Schreiber2025_tangermeme_all.pdf`. The design documentation specifies a default DNA alphabet of A/C/G/T and single-character edits, not codons.

**Disconfirmation attempt / what would change the score.** A documented call accepting a frame/genetic code and mapping codons to a fixed residue replacement set would change this to `partial`; a built-in replacement-policy choice would change it to `yes`. The complete-tree vocabulary search and focused API/paper inspection found neither. This is consequently `no`, not `unknown`.

## Search and provenance log

Only primary sources were used as evidence: each official repository at the pinned revision above, its README/documentation/source, and the supplied author papers. Per-tool survey records were treated as locator leads only; no statement or rating from them was used or cited.

Searches performed:

1. Read `ROW_DEFINITIONS.md:1-34` and `:117-124`, then stated the operational test above before opening tool evidence.
2. Searched the local PoolParty repository with `rg -n -i` over Python/RST/Markdown for `codon`, `amino`, `residue`, `translation`, `mutation_type`, `NNK`, and related policy terms; an initial broad search also encountered generated SVG text, which was excluded as non-source noise. Inspected `mutagenize_orf.rst`, the API implementation, `codon_table.py`, and relevant tests. Ran the two read-only behavioral checks described above.
3. Extracted text from every supplied PDF with PyMuPDF (`fitz`) and searched it for the same coding/policy terms; page-numbered evidence was rechecked per PDF. This was done for VaLiAnT, MPRAnator, MPRADesignTools, Mutation Maker, DnaChisel, and tangermeme. PoolParty's author PDF was used only as a pointer, per the task instruction; the local repository supplied its evidence.
4. Queried the official GitHub trees and pinned branch heads for all seven remote tools. For each `no`, searched the entire meaningful tree, then manually inspected the documented entry points and mutation implementations. The exact exhaustive patterns and file sets are given in each tool section.
5. A request for `hemberg-lab/MPRAnator/iliasApp/parseSNPs.py` returned 404; the tree showed and I inspected the actual repository-root `parseSNPs.py`. Requests for DnaChisel's old flat `CodonOptimize.py` path and tangermeme's nonexistent `tangermeme/design.py` likewise returned 404; I used the tree to locate DnaChisel's current `builtin_specifications/codon_optimization/CodonOptimize.py` and inspected tangermeme's actual design modules. These failed paths contributed no evidence.
6. GitHub access initially failed inside the restricted network sandbox; the same read-only official-repository queries were rerun with approved network access. No package was installed and no source repository was modified.

## Row-level finding

The row discriminates on a single scale: three tools pass with multiple built-in codon replacement policies (PoolParty, VaLiAnT, DnaChisel), Mutation Maker is restricted to a fixed tool-selected NNK set unless the caller supplies the replacement set, and four tools lack a codon↔residue substitution operation. The important boundary is not simply whether the word *codon* appears: Oligopool Calculator's example enumerates 64 triplet strings but never maps them to residues, whereas DnaChisel's single-output codon optimization qualifies because the row explicitly admits codon-level most-frequent-codon policies and does not require residue-changing or library-enumerating output.
