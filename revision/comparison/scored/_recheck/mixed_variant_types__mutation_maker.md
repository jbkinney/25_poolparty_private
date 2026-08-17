# Mutation Maker — mixed variant types in one library (row 6)

- **Value:** `partial`
- **Confidence:** High (0.95)
- **Version inspected:** maintained Mutation Maker repository owned by tool coauthor Rastislav Švarba, commit `396c1c0ede7529f3dedf65e821e8c1f20c9a7043` (2026-02-14), whose history contains the original Merck publication commits; tool paper: Hiraga et al., *ACS Synthetic Biology* 2021, DOI `10.1021/acssynbio.0c00542`

## Binding operational test

Row 6 requires two or more structurally different component types to be declared in one specification and emitted as one pooled output. Separate tool runs that the caller merges are `no`. Under the global rule, a user cannot reconstruct the capability by combining unrelated workflows.

## Primary-source evidence

Mutation Maker exposes its three design modalities as separate request endpoints:

> `@hug.post('/ssm', versions=1)`  
> `Return task if type Site Saturation Mutagenesis (SSM)`

> `@hug.post('/qclm', versions=1)`  
> `Return task if type QuikChange Lightning Multi Site-Directed Mutagenesis (QCLM)`

> `@hug.post('/pas', versions=1)`  
> `Return task if type PCR-based Accurate Synthesis task (PAS)`

Source: [`api/api.py`, lines 40–61, commit `396c1c0...`](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/api/api.py#L40-L61).

Each workflow's first-party description is substitution-only: SSSM performs “random amino acid substitutions,” MSDM combines “specific amino acid changes,” and PAS combines “specific amino acids changes.” Source: [`frontend/src/App.tsx`, lines 149–189, same commit](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/frontend/src/App.tsx#L149-L189).

PAS does provide a one-request, one-output mixed pool within a narrow substitution family. Its request declares mutation choices and frequencies by position:

> `class PASMutationFormattedInput(JsonObject):`  
> `    mutants = ListProperty(required=True)`  
> `    position = IntegerProperty(required=True)`  
> `    frequency = FloatProperty(required=True)`

and each emitted oligo carries an explicit pool mixing ratio:

> `class PASOligoOutput(JsonObject):`  
> `    sequence = StringProperty(required=True)`  
> `    mix_ratio = FloatProperty(required=True)`

Source: [`backend/mutation_maker/pas_types.py`, lines 145–173, same commit](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/pas_types.py#L145-L173).

The tool itself adds the residual wild-type component, generates cross-site combinations, and emits them together:

> `if sum_of_probabilities != 1:`  
> `    wild_type_prob = 1 - sum_of_probabilities`  
> `    ... # creating wild type record with corresponding codon`  
> `mutations_combinations_with_probabilitites = Mutations.generate_mutation_combinations(...)`  
> `return generate_oligos_from_combinations(...)`

Source: [`backend/mutation_maker/generate_oligos.py`, lines 238–309, same commit](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/generate_oligos.py#L238-L309). `generate_oligos_from_combinations` says it generates DNA sequence “for all combinations,” appends each to one `oligos` list, normalizes their ratios to sum to one, and returns that list: [lines 122–142](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/generate_oligos.py#L122-L142).

The shared mutation parser confirms the restriction:

> `Represents a mutation at a single location from one AA to another.`  
> `self.length = 3`

It accepts an original amino acid, a position, and one target amino acid (or `X`), and all mutation objects remain length three. Source: [`backend/mutation_maker/mutation.py`, lines 28–71, same commit](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/mutation.py#L28-L71).

Even the most flexible PAS request has one global input-type switch rather than mixed types:

> `# List of mutations by position in the gene`  
> `mutations = ListProperty(PASMutationFormattedInput, required=False)`  
> `# Are mutations given as codons?`  
> `is_mutations_as_codons = BooleanProperty(required=True)`

Source: [`backend/mutation_maker/pas_types.py`, lines 199–210, same commit](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/pas_types.py#L199-L210). The parser branches the entire request on that one Boolean: [`generate_oligos.py`, lines 46–57](https://github.com/ra100/Mutation_Maker/blob/396c1c0ede7529f3dedf65e821e8c1f20c9a7043/backend/mutation_maker/generate_oligos.py#L46-L57).

The tool paper likewise says that Mutation Maker workflows were designed primarily to support SSSM, MSDM/QCLM, and PAS as distinct modalities, and describes their output mutations as codon/amino-acid substitutions. Source: [Hiraga et al. 2021, first-party tool paper, Abstract/Results and Figure 1](https://pubs.acs.org/doi/10.1021/acssynbio.0c00542).

## Reasoning

A single PAS specification declares multiple target sites/alternatives and their frequencies. Mutation Maker supplies the residual wild-type state when a site's declared mutant frequencies sum below one, computes the Cartesian product across sites, and returns all resulting oligos together with normalized mixing ratios. Thus one output can include parental, single-substitution, and multi-substitution components; the caller does not run separate workflows or concatenate their results.

This is `partial`, not `yes`, because the mixed components are restricted to one structural mutation family: fixed-length codon/amino-acid substitutions and their combinations. The common mutation model has invariant `length = 3`, PAS globally chooses amino-acid versus codon interpretation, and no indel or other structurally broader component class can be co-declared. The separate SSSM, QCLM, and PAS endpoints cannot themselves be combined for extra credit.

## Disconfirmation attempt

I searched the complete maintained repository at the pinned commit across backend schemas, parsers, algorithms, API routes, frontend forms, exporters, and tests for `insertion`, `deletion`, `indel`, `mutation type`, pooling/merge behavior, wild-type/parental handling, degeneracy, mutation ratios, and random/directed workflows. I inspected all three input schemas (`SSMInput`, `QCLMInput`, `PASInput`), the common mutation parser, PAS Cartesian-product and ratio generation, output schemas, and the three API endpoints. I also checked the tool paper's descriptions of SSSM, MSDM, and PAS. This disconfirmed `no` by locating the PAS one-output mixture, but found no supported insertion, deletion, shuffling, recombination, or other broader structural class, and no endpoint accepts multiple workflow types or pools their outputs.

The score would change to `yes` if one supported specification could co-declare broader structurally distinct classes—for example substitutions plus insertions/deletions—and Mutation Maker emitted them as one pool. It would change to `no` if PAS merely generated one combination per run or required the caller to append combinations, but the source explicitly performs the Cartesian product and returns one normalized oligo list. Neither condition for changing the final `partial` value was found.
