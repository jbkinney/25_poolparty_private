# Oligopool Calculator — saturation mutagenesis (row 3)

- **Value:** `partial`
- **Confidence:** High (0.98)
- **Version inspected:** Oligopool Calculator `2026.02.22.1`, canonical repository commit `b88fa394ca67ed4c48ec41127b5707694ee7cf0a` (2026-02-22)

## Binding operational test

The row requires one specification to make the tool enumerate every single substitution at every eligible position. It assigns `partial` when that exhaustive capability exists only in an example script rather than the documented API. The global rule additionally requires a tool-provided operation, parameter, or mode rather than user reconstruction.

## Primary-source evidence

The canonical repository contains a first-party example helper with a single parent-sequence input:

> `def generate_single_mutants(`  
> `    sequence: str,`  
> `    include_wildtype: bool = True,`  
> `)`  
> `Generate all single-nucleotide mutants of a DNA sequence.`  
> `For a sequence of length L, this generates up to 3*L mutants.`

The implementation exhausts both dimensions required by the row: all positions and all non-reference substitutions:

> `for pos in range(len(sequence)):`  
> `    ref_base = sequence[pos]`  
> `    for alt_base in DNA_BASES:`  
> `        if alt_base != ref_base:`  
> `            ...`  
> `            yield (variant_id, mutant_seq, pos + 1, ref_base, alt_base)`

The same example exposes this through one strategy selection, without requiring the caller to provide positions or substitutions:

> `if strategy == 'single':`  
> `    data = [(vid, seq) for vid, seq, *_ in generate_single_mutants(sequence, include_wildtype)]`

Sources: [`examples/library-compressor/mutant_generator.py`, lines 43–79, commit `b88fa394...`](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/examples/library-compressor/mutant_generator.py#L43-L79) and [lines 200–233](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/examples/library-compressor/mutant_generator.py#L200-L233). The repository labels this file as an example helper: “`mutant_generator.py` - Helper utilities for generating test variants.” Source: [`examples/library-compressor/README.md`, lines 1–8, same commit](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/examples/library-compressor/README.md#L1-L8).

By contrast, the first-party user guide's documented saturation-mutagenesis recipe begins with variants already enumerated:

> **The decomposition**: Generate all substitution variants, then `compress` into IUPAC-degenerate oligos.

and its input is explicitly user-supplied:

> `# Your substitution variants (e.g., alanine scanning, full saturation)`  
> `'Sequence': [...]  # 1000 single-substitution variants, strict ATGC`

Source: [`docs/docs.md`, lines 2331–2353, same commit](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/docs.md#L2331-L2353).

## Reasoning

`generate_single_mutants(sequence)` is a tool-authored operation that takes one parent sequence, derives all eligible positions with `range(len(sequence))`, derives the substitution set internally from `DNA_BASES`, excludes only the reference base, and emits every resulting single-nucleotide substitution. It therefore meets the substantive exhaustive-enumeration test, and the user does not supply a position or substitution list.

However, the function lives under `examples/library-compressor/`, is described as a helper for test variants, and is absent from the documented package API. The documented saturation workflow instead expects the caller to arrive with the substitution variants and only offers `compress` afterward. The row definition explicitly maps an exhaustive implementation available only in an example script to `partial`; therefore the correct value is `partial`, not `yes`. It is not `no` because the first-party example really does implement the enumeration rather than merely showing user-written pseudocode or requiring a supplied mutation list.

## Disconfirmation attempt

I actively searched the complete canonical repository at the pinned commit for `saturat`, `mutagen`, `mutation`, `substitut`, `variant`, `generate_single`, and `single mutant`; inspected the complete API table of contents in `docs/api.md`; inspected the package export registry (`oligopool/__init__.py`, `__api__` and `_LAZY_ATTRS`); and enumerated the public top-level module functions. This found no documented/package-exported saturation-mutagenesis operation or parameter. The only exhaustive single-substitution generator found is the example helper above. The user guide and agent guide instead say to “generate all substitutions” before using the package's documented `compress` workflow.

Evidence that would change the score to `yes` would be a documented package/API/CLI operation accepting a parent sequence or region specification and internally enumerating every position and alternative substitution. Evidence that the example helper required caller-supplied positions/substitutions, skipped eligible positions/alternatives, or was not first-party tool material would instead change the score to `no`. None was found.

Disconfirmation locations: [`docs/api.md`, lines 14–48](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/api.md#L14-L48); [`oligopool/__init__.py`, lines 6–31 and 42–69](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/__init__.py#L6-L31); [`docs/agent-skills.md`, lines 371–378](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/agent-skills.md#L371-L378).
