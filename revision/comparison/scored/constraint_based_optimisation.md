# Audit: constraint-based sequence optimisation

> Row renamed after scoring: `constraint_repaired_variants` -> `constraint_based_optimisation`. Label only — the operational test was unchanged (the tool *modifies* a sequence to satisfy declared constraints; rejection-only filtering does not count), so the scores stand without rescoring. Renamed because "repair" is not this field's vocabulary: "optimisation" appears 110 times across 11/11 tool papers, "repair" 7 times across 3, and never in the sequence-design sense.

Date: 2026-08-15

## Measurement instrument and operational test

I read `revision/tables/ROW_DEFINITIONS.md` before examining any tool. The binding global rule says that a capability counts only when the tool supplies an operation, parameter, or mode for it; user reconstruction from unrelated primitives is `no`, and `partial` is a restricted tool-provided capability, not user assembly. Row 11 says:

> The tool **modifies** a sequence so it satisfies declared constraints, rather than only rejecting sequences that violate them. **Rejection-only filtering does not count** — this row is specifically about repair. `partial` = repair for a narrow constraint class only.

Operational test, stated before tool inspection: **For each tool, determine whether a documented tool-provided operation can take a sequence that violates a user-declared constraint and return a modified sequence that satisfies it—`yes` for general or multiple constraint classes, `partial` for repair restricted to a narrow class, and `no` when the tool only rejects/filters or exposes no repair operation after repository, documentation, and paper searches.** An in-operation candidate is a sequence for this test, but it must actually be edited; drawing a new candidate or discarding the old one is rejection/generation, not repair. A model objective is not silently reinterpreted as a sequence constraint, because doing so would require user reconstruction prohibited by the global rule.

## Results

| Tool | Value | Confidence |
|---|---|---|
| PoolParty | **no** | high |
| VaLiAnT | **no** | high |
| MPRAnator | **partial** | high |
| MPRA Design Tools | **partial** | high |
| Oligopool Calculator | **yes** | high |
| Mutation Maker | **yes** | high |
| DNA Chisel | **yes** | high |
| tangermeme | **no** | high |

### PoolParty — no

**Positive-looking evidence tested.** `get_barcodes` exposes several constraints, but its documented operation begins without an input sequence and uses candidate rejection. The documentation says:

> “Generate a pool of DNA barcodes satisfying distance and quality constraints. All barcodes are pre-generated at construction time using a greedy random algorithm.”

Source: canonical repository commit `1bb0179e1c37`, `poolparty/docs/operations/get_barcodes.rst:4-7`; its constraint list is at `:9-12`.

The source makes the rejection behavior explicit:

> `candidates = rng.integers(0, 4, size=(batch_size, length), dtype=np.uint8)`
>
> `candidates = candidates[_gc_filter_batch(candidates, gc_range)]`
>
> `candidates = candidates[_homopolymer_filter_batch(candidates, max_homopolymer)]`

and later rejects distance failures with `continue` before accepting the candidate. Source: `poolparty/src/poolparty/base_ops/get_barcodes.py:73-79,157-200`. This is candidate filtering, not editing a failing candidate.

The general predicate API is also rejection-only:

> “If the predicate returns False for a sequence, the operation returns `NullSeq`, which propagates through downstream operations.”

Source: `poolparty/src/poolparty/base_ops/filter_seq.py:13-18`; implementation returns the original parent on pass and `NullSeq()` on failure at `:60-65`.

**Behavioral check required for this repository.** With `PYTHONDONTWRITEBYTECODE=1` and the repository's read-only Python 3.12 environment, I ran `get_barcodes(num_barcodes=4, length=6, gc_range=(0.5,0.5), max_homopolymer=2, seed=7)` and obtained `GTGTCT, CACTGA, GACCAA, GGTCAA`; these were newly generated candidates. I then filtered `['AAAA', 'ACGT']` with a predicate rejecting `AAAA`; output contained only `ACGT`, not a repaired version of `AAAA`. This agrees with the quoted source.

**Disconfirmation attempt.** I searched all 177 relevant `.py`, `.rst`, and `.md` files under `poolparty/src`, `poolparty/docs`, and the README for `repair|constraint|resolve_constraints|satisf|violat|filter|reject|optimi[sz]|forbid|avoid|restriction|homopolymer|hairpin|hamming|gc|melting|temperature|length`. I read the complete relevant APIs and docs for `get_barcodes`, `filter_seq`, `shuffle_seq`, sequence-property checkers, and `reverse_translate`. `shuffle_seq` only preserves mono-/dinucleotide composition (`docs/operations/shuffle_seq.rst:35-40`), and reverse translation creates a DNA sequence from protein rather than repairing a violating DNA input. A tool-provided operation that edits a failing input/candidate (rather than discarding it) under any declared constraint would have changed this to at least `partial`; a broad constraint solver would have changed it to `yes`. I searched for both and found neither.

### VaLiAnT — no

VaLiAnT's synthesis-related declared constraint is a length filter. The README says:

> “Oligonucleotides exceeding a given length (`max-length` option) will not be included in the unique oligonucleotide files and their metadata will be stored in separate files marked as 'excluded'.”

Source: canonical repository commit `8796cc112daf`, `README.md:112`.

The implementation returns a Boolean and increments exclusion counters; it does not shorten an oligo:

> `elif oligo_length > opt.oligo_max_length:`
>
> `    self.too_long += 1`
>
> `    return False`

Source: `src/valiant/oligo_generation_info.py:54-65`. The paper independently describes the same behavior: sequences over the user-defined length are “excluded from standard output files,” with separate metadata (Barbon et al. 2022, main text p. 5 of the article, lines corresponding to the local PDF extraction around 305-307).

Potential PAM protection was checked separately. It does not select a repair edit: the documented input is a user-authored “VCF file containing single-nucleotide substitution variants linked to sgRNA identifiers via the `SGRNA` tag,” followed by explicit REF/ALT examples. Source: `README.md:720-735`. Applying pre-specified edits is not a tool-provided constraint-repair decision under the global rule.

**Disconfirmation attempt.** I searched the full `src/valiant` tree, all examples, README, CHANGELOG, configuration JSON, and the VaLiAnT paper for the general term set above and separately for `PAM|protect|protection`, `repair`, `constraint`, `optimize`, and all common synthesis constraints. I examined `oligo_generation_info.py`, `meta_table.py`, `sge_proc.py`, `pam_variant.py`, CLI/config code, and the README's length/PAM/background sections. Other behavior is also discard/error based: background variants may be filtered (`README.md:116`), invalid backgrounds raise errors (`:118`), and overlapping mutations are discarded (`:120`). A mode that automatically chooses an edit to shorten an overlong sequence, remove a site, repair GC/repeats/structure, or choose a PAM-protection edit would have produced `partial` or `yes`; I explicitly searched for such a mode and found only exclusion and user-supplied variants.

### MPRAnator — partial

MPRAnator supplies narrow repair inside its barcode generator. It starts with a random barcode candidate, tests the user-declared GC bounds, and, on failure, changes bases until the bounds pass:

> `while short_term in barcode_storage or gc_cont(short_term) * 100 > maxgc ... or gc_cont(short_term) * 100 < mingc ...:`
>
> `    position = rd.choice(range(barCodeLength))`
>
> `    short_term = short_term.replace(short_term[position], rd.choice(l1))`

Source: canonical repository commit `9969790d6241`, `part1.py:52-88`, particularly `:82-87`. The user-facing path passes `minimumGCContent` and `maximumGCContent` into this generator (`part1.py:295-315`; form labels in `iliasApp/forms.py:392-393,439-447`). This is actual sequence modification rather than simply dropping the candidate, so it passes the repair test.

The restriction-site feature is not repair. It reports a duplicate in the output header:

> `outputHeader += "|DUPLICATE_RESTRICTION_SITES - RESTRICTION1"`

Source: `part1.py:277-286`; the application's documentation says the label “Report[s] the restriction site which has multiple copies present in the oligo” (`iliasApp/views.py:232-243`).

**Why partial.** Repair is restricted to the internally generated barcode and essentially its GC/duplicate loop; restriction sites are only flagged, and motif spacing is handled by enumerating admissible placements rather than repairing a violating sequence.

**Disconfirmation attempt.** I searched all 28 repository `.py`/`.md` files plus the paper for the general term set, then read `part1.py`, `part3.py`, `myfunctions.py`, `oligo.py`, forms/views, both downloadable scripts, and the restriction/GC matches. Broader repair for assembled oligos, input backgrounds, restriction sites, repeats, length, or structure would have changed the score to `yes`; zero such operation was found. If the GC loop had redrawn/discarded candidates rather than assigning a changed `short_term`, the value would have been `no`; I checked that exact branch and confirmed the quoted edit.

### MPRA Design Tools — partial

The package exposes two narrow repairs. First, `randomly_fix` is a documented API whose title is “Randomly correct aberrant digestion sites” and whose description is:

> “For a SNP with aberrant digestion sites in the context, randomly change bases in the site across barcodes.”

Source: canonical `mpradesigntools` repository commit `afd386ef1205`, `man/randomly_fix.Rd:5-24`; corresponding source `R/processVCFfast.R:112-127`. The implementation selects a nonmatching allele and writes it back into the construct:

> `fixed_pattern = rep(altered_patterns$altered_pattern, times = 2)`
>
> `constrseq_fixed = map2_chr(constrseq, fixed_pattern, reassign_pattern, ...)`

Source: `R/processVCFfast.R:150-172`. `processVCF` exposes this through `alter_aberrant` (`:1019-1020,1099-1115`) and calls `randomly_fix` before rebuilding and recounting the sequence (`:452-481`).

Second, overlong constructs can be shortened:

> “If provided, constructs that end up longer than this have sequence context evenly removed from both sides until sufficiently short.”

Source: `R/processVCFfast.R:1024-1027`; implementation adjusts both context ranges at `:317-324`.

The same API also documents rejection-only neighbors: `filterPatterns` removes barcodes from the pool (`man/processVCF.Rd:89-92`), and some complex digestion-site cases return failure (`R/processVCFfast.R:421-449`). Those do not add credit.

**Why partial.** Repair is limited to aberrant restriction/digestion sites and a maximum construct length within this fixed MPRA construction workflow; `alter_aberrant` and `max_construct_size` are both labeled under development. There is no general constraint framework.

**Disconfirmation attempt.** I searched all 36 `.R`, `.Rd`, and `.md` files in `mpradesigntools`, all 16 `.R`, `.Rmd`, and `.md` files in the companion `designMPRA` repository at commit `0cf56ee602fc`, and the paper/supplement PDF for the general term set plus `randomly_fix|alter_aberrant|aberrant|digestion|restriction|max_construct_size|shorten`. I read the full relevant sections of `processVCFfast.R`, generated Rd pages, README, and the older companion scripts. A tool-provided solver spanning independently declared GC, repeat, structure, homology, motif, and other constraints would have changed the score to `yes`; searches found filtering/precomputed barcode pools but no broader solver. If `alter_aberrant` only resampled barcodes or failed, it would be `no`; the current package source explicitly changes the genomic-context bases, so `partial` is warranted.

### Oligopool Calculator — yes

The documented API directly declares multiple classes of constraints. For example, `primer`:

> “Design constrained primers under an IUPAC sequence constraint with Tm/repeat/hairpin/dimer screening.”

and exposes a sequence constraint, minimum/maximum melting temperature, maximum repeat length, excluded motifs, background databases, paired-primer matching, and context (`oligopool/primer.py:17-39,41-64,78-90`). `motif` similarly exposes IUPAC, repeats, excluded motifs, background, and junction-aware context (`oligopool/motif.py:17-40,51-68`). The README summarizes “constraint-based design” and “Rich constraints” covering IUPAC sequence constraints, motif exclusion, repeats, Hamming distance, and primer thermodynamics (`README.md:39,46`).

This is a solver, not only a verifier. For motifs, the engine defines an objective over repeat, excluded-motif, edge, and background constraints (`oligopool/base/core_motif.py:995-1016`), passes it together with the sequence/structure constraint into `NRPMaker.nrp_maker` (`:1018-1034`), and records the designed sequence with status `solved` (`:1043-1051,1092-1098`). Primer design does the analogous operation for Tm, repeats, paired/cross-primer compatibility, excluded motifs, junctions, background, and structure (`oligopool/base/core_primer.py:2093-2141`), returning the designed primer and computed Tm/GC/structure statistics (`:2143-2185`).

The paper corroborates breadth and search behavior. It states that primer pairs satisfy seven design criteria, including excluded sequences, Tm, off-target binding, secondary structure, terminal bases, homo- and heterodimers (Hossain et al. 2024, pp. 4221-4222), and that sequence, structure, and function constraints are applied “during path generation”; the procedure uses “path generation and tracebacks” (p. 4222). It also says motif/spacer pathfinding “evaluates the input specifications and identifies an oligonucleotide-specific sequence solution” avoiding researcher-defined motifs and homopolymers (p. 4222).

**Disconfirmation attempt.** I searched all 46 `.py`/`.md` package and documentation files and the paper for the general term set plus `candidate|path|traceback|maker|objective|check_constraints|patch_mode`. I read the public `barcode`, `primer`, `motif`, `spacer`, `pad`, `verify`, and background interfaces and the three core design engines. Barcode candidates are indeed filtered (`core_barcode.py:1073-1173`; paper p. 4221), which alone would be `no`, but primer/motif/spacer paths synthesize/edit sequences under directly declared constraints and are not user reconstruction. If every design module only discarded whole candidates, the score would have been `no`; if repair existed for only one narrow class, `partial`. The multiple direct sequence, structural, thermodynamic, repeat, motif, distance, background, and edge constraints support `yes`.

### Mutation Maker — yes

Mutation Maker supplies multiple explicit constraint-satisfaction and backtracking workflows. The paper describes it as enabled by “optimization, constraint-satisfaction and backtracking algorithms” (Hiraga et al. 2021, p. 357) and gives the concrete repair behavior:

> “the fast-approximation algorithm [works] by which mutagenic primers are dynamically extended until they satisfy the design criteria.”

Source: paper p. 359. It lists user-adjustable hairpin/dimer checks, GC, primer length, 5′/3′ lengths, 3′ Tm, overlap size/Tm, and cross-primer Tm consistency on the same page.

The repository implements that edit directly:

> “Grows a forward primer from a given overlap, and returns the shortest primer which has 3' Tm above the temperature threshold.”

Source: canonical repository commit `396c1c0ede75`, `backend/mutation_maker/ssm_fast_approximation.py:72-88`; the reverse-primer equivalent is at `:91-109`. The public configuration exposes min/opt/max primer size, GC, end/overlap size and Tm, GC clamp, hairpin and dimer scoring/checks (`backend/mutation_maker/ssm_types.py:132-175`).

The PAS workflow adds fragment/overlap lengths and Tm, GC, organism/codon usage, avoided motifs, and hairpin/homodimer settings (`pas_types.py:29-66`) and uses a backtracking/dynamic-programming `PASOptimizer` (`pas_back_track.py:32-64,83-95,149-215`). Thus this is not one narrow predicate or only candidate rejection.

**Caveat considered.** The SSM paper says primer/3′-end lengths are hard constraints while remaining SSM parameters are soft and may be violated (p. 359); it reports violations rather than pretending all soft preferences passed. That does not negate the capability because tool-provided hard repair exists, and the overall platform supplies multiple constraint-solving workflows, including PAS backtracking and motif-aware reverse translation. It does mean this score is for the presence of broad repair operations, not a guarantee that every optional preference is hard.

**Disconfirmation attempt.** I searched all 72 relevant backend `.py`/`.md` files, the root README, and the paper for the general term set plus `grow|extend|backtrack|PASOptimizer|hard constraint|soft|avoided_motifs|reverse translation`. I read all SSM/QCLM/PAS config types, the fast-approximation grow/score functions, PAS backtracking, reverse translation, and paper algorithm sections. If the program only enumerated and ranked fixed primers, or merely reported violations, it would have been `no`; the explicit primer-growth repair disproves that. If only this one Tm/length repair existed, it would have been `partial`; the independent SSM, QCLM, and PAS constraints/backtracking across geometry, thermodynamics, GC, structure, motifs, and codon choice support `yes`.

### DNA Chisel — yes

DNA Chisel is the clearest positive control. Its README says it optimizes DNA with respect to constraints/objectives and supplies more than 15 composable specification classes, including provider constraints, homology avoidance, and GC tuning (`README.rst:19-31`). Its runnable example declares BsaI avoidance, windowed GC, and translation preservation, calls:

> `problem.resolve_constraints()`
>
> `problem.optimize()`

and retrieves `problem.sequence` as the final sequence (`README.rst:47-82`). The implementation stores the original sequence and exposes sequence edits (`dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py:167-180,191-200`).

The documentation states the repair mechanism directly:

> “DNA Chisel hunts down every constraint breach and suboptimal region by recreating local version of the problem around these regions. Each type of constraint can be locally reduced and solved in its own way.”

Source: `README.rst:150-156`. The paper likewise defines a starting sequence plus hard constraints that “must be satisfied in the final sequence,” reports modified nucleotides, and says unsatisfactory regions are optimized by random mutations or exhaustive variants (Zulkower & Rosser 2020, pp. 4508-4509).

**Disconfirmation attempt.** I searched all 97 `.py`/`.rst` source and documentation files and the paper for the general term set plus `resolve_constraints|breach|mutation_space|local resolution|NoSolutionError`. I read the README example/how-it-works section, core problem class, constraint solver mixin, built-in specification index/classes, CLI, GenBank API, and the paper implementation section. I specifically checked the negative boundary: `SequenceLengthBounds` says it “can't really be solved or optimized” (`builtin_specifications/SequenceLengthBounds.py:9-10`), so that one checker does not inflate the score. If the framework only evaluated specifications or only one narrow class were solvable, the score would be `no`/`partial`; broad, composable, directly invoked repair is explicit, so `yes`.

### tangermeme — no

tangermeme does edit sequences, but its documented design operation is model-objective optimization, not declared sequence-constraint repair. `greedy_substitution` takes `model`, desired model output `y`, motifs, a model-output `loss`, masks, tolerance, and iteration count (`tangermeme/design/greedy_substitution.py:24-40`). Its description is:

> “This design function will greedily add motifs to achieve a desired output from the model.”

Source: canonical repository commit `8d732b8c0876`, `tangermeme/design/greedy_substitution.py:41-53`. The README says the algorithm tries motif/position combinations and takes the one whose prediction is closest to the desired model output (`README.md:268-280`). This belongs to model-guided design, not a tool-provided constraint declaration.

The GC routine is only a calculator:

> “Compute the GC content of one-hot encoded sequences.”

Source: `tangermeme/utils.py:896-922`. It is not wired into a repair operation. Locus exclusion and GC matching elsewhere in `io.py`/`match.py` filter or sample genomic regions, not modify sequences.

**Disconfirmation attempt.** I searched all 66 `.py`, `.rst`, and `.md` source/docs files, all design modules and tutorials/API pages, bundled skill reference docs, and the paper for the general term set plus `loss|objective|input_mask|GC|restriction|homopolymer|hairpin|melting|sequence constraint`. There were no tool-provided sequence-constraint parameters or modes for repair. What would change the score: a documented constraint argument/mode that directly repairs GC, motifs, repeats, structure, distance, or similar would yield at least `partial`; a general sequence-constraint interface would yield `yes`. I searched for both. A user could write a predictive model/loss that represents such a constraint, but that is precisely user reconstruction from an interface provided for model targets and therefore does not count under the binding global rule.

## Search and source audit trail

The per-tool `final/*.md` records were used only once for repository/page leads. I ran a lead-only `rg` over those files for repository, documentation, and likely source-file names; no statement, quote, score, or conclusion from them appears as evidence above.

Primary snapshots inspected:

| Tool | Repository snapshot |
|---|---|
| PoolParty | `jbkinney/poolparty-statetracker` `1bb0179e1c37` (local canonical read-only checkout) |
| VaLiAnT | `cancerit/VaLiAnT` `8796cc112daf` (`develop`) |
| MPRAnator | `hemberg-lab/MPRAnator` `9969790d6241` |
| MPRA Design Tools | `andrewGhazi/mpradesigntools` `afd386ef1205`; companion `andrewGhazi/designMPRA` `0cf56ee602fc` |
| Oligopool Calculator | `ayaanhossain/oligopool` `b88fa394ca67` |
| Mutation Maker | `ra100/Mutation_Maker` `396c1c0ede75` |
| DNA Chisel | `Edinburgh-Genome-Foundry/DnaChisel` `68c09304341c` |
| tangermeme | `jmschrei/tangermeme` `8d732b8c0876` |

For each snapshot I first enumerated files with `find`/`rg --files`, recorded the Git origin and commit with `git rev-parse`/`git remote get-url`, then ran recursive case-insensitive searches. The common search expression, sometimes split across invocations to keep output readable, was:

```text
constraint|repair|resolve_constraints|satisf|violate|breach|forbid|avoid|
filter|reject|optimi[sz]|homopolymer|gc(_content|_range)?|restriction|
hairpin|secondary structure|homodimer|heterodimer|melting|temperature|
hamming|repeat|length|candidate|modify|change|backtrack|traceback
```

Tool-specific follow-ups were: PoolParty `get_barcodes|filter_seq|shuffle_seq|reverse_translate`; VaLiAnT `PAM|protect|protection|excluded|eval_in_range`; MPRAnator `barcode_generator|minimumGCContent|duplicate restriction`; MPRA Design Tools `randomly_fix|alter_aberrant|digestion|max_construct_size|shorten`; Oligopool Calculator `maker|objective|path|traceback|patch_mode|check_constraints`; Mutation Maker `grow|extend|PASOptimizer|hard|soft|avoided_motifs|reverse_translation`; DNA Chisel `resolve_constraints|mutation_space|local resolution|NoSolutionError`; tangermeme `loss|model|input_mask|sequence constraint|GC`.

PDFs were text-extracted with PyMuPDF as instructed and searched with `rg -n -i -C`: `Barbon2022_VaLiAnT_all.pdf`, `Georgakopoulos-Soares2017_MPRAnator.pdf`, `Ghazi2018_MPRADesignTools_all.pdf`, `Hossain2024_OligopoolCalculator.pdf`, `Hiraga2021_MutationMaker.pdf`, `Zulkower2020_dnachisel.pdf`, and `Schreiber2025_tangermeme_all.pdf`. Search terms were the common expression narrowed to each paper's vocabulary (`constraint`, `repair`, `optimize`, `filter`, `length`, `restriction`, `GC`, `barcode`, `primer`, `motif`, `structure`, `greedy`, `design`). PoolParty's developer-authored PDF was treated only as a pointer as required; repository source/docs and the live read-only execution determined its result.

## Row-level finding

The row discriminates on one consistent scale: three tools provide broad, directly declared constraint-solving repair (`yes`); two actually edit sequences but only for narrow, workflow-specific constraint classes (`partial`); and three either reject/filter or offer edits only for a different declared purpose (`no`). The most important boundary is observable in closely comparable barcode code: PoolParty filters out failing random candidates (`no`), while MPRAnator mutates its failing barcode candidate until the declared GC bounds pass (`partial`). No `unknown` was needed.
