# DNA Chisel — saturation mutagenesis (row 3)

- **Value:** `no`
- **Confidence:** high
- **Tool version / source revision checked:** DNA Chisel 3.2.16, canonical repository commit `68c09304341c3656f3dfe63eda37757d6a7b3917` (current `master` checked 2026-08-17).

## Primary-source evidence

1. The documented `MutationSpace.all_variants(sequence)` API says: **“Iterate through all sequence variants in this mutation space.”** Its implementation builds one slot for every mutable choice and then iterates `itertools.product(*variants_slots)`, applying every selected choice to the same output sequence. Thus, over a multi-position region it emits the Cartesian product (including multi-position combinations), not the library consisting of each single substitution at each position.

   Source: canonical source, [`dnachisel/MutationSpace/MutationSpace.py`, lines 132–164 at commit `68c0930`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/MutationSpace/MutationSpace.py#L132-L164).

2. The one-mutation path is random rather than exhaustive. `pick_random_mutations` is documented as **“Draw N random mutations.”** When `n_mutations == 1`, it executes `index = np.random.randint(len(self.multichoices))`, selects that one choice, and returns one random alternative. `apply_random_mutations` merely applies the returned random choices. This cannot enumerate every substitution at every position from one call/specification.

   Source: canonical source, [`dnachisel/MutationSpace/MutationSpace.py`, lines 106–130 at commit `68c0930`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/MutationSpace/MutationSpace.py#L106-L130).

3. A default mutation space does make every nucleotide position mutable: `from_optimization_problem` defines the alternatives for each original base and constructs one `MutationChoice((i, i + 1), ...)` for every sequence position. However, the only provided restriction operation is `localized(location)`, documented as **“Return a new version with only mutations overlapping the location,”** and it requires the caller to supply that location. Combining repeated caller-chosen locations with `all_variants` and aggregating the results would therefore be caller-authored enumeration/bookkeeping, which the binding global rule excludes.

   Sources: canonical source, [`dnachisel/MutationSpace/MutationSpace.py`, lines 88–94](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/MutationSpace/MutationSpace.py#L88-L94) and [lines 166–180](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/MutationSpace/MutationSpace.py#L166-L180), both at commit `68c0930`.

## Reasoning against the operational definition

Row 3 requires one tool-provided specification/operation to enumerate every **single** substitution at every eligible position. DNA Chisel provides (a) Cartesian enumeration of the entire mutation space and (b) random selection of a requested number of mutation choices. It does not provide a saturation-mutagenesis operation that emits the exhaustive single-substitution set. A caller could loop over positions, localize the mutation space, invoke `all_variants`, remove the reference allele, and concatenate outputs, but the row explicitly says that supplying positions/substitutions is not tool enumeration, and the global rule classifies such caller reconstruction as `no`, not `partial`.

## Disconfirmation attempt

I searched the complete first-party repository at the pinned commit—including `dnachisel/`, `docs/`, `examples/`, and `tests/`—case-insensitively for `saturat`, `mutagen`, `single mutation`, `single_mut`, `all_mut`, `all_variants`, `pick_random_mutations`, `apply_random_mutations`, and all function/class names containing `variant` or `mutat`. I also inspected the generated first-party Core Classes documentation for every public `MutationSpace` method. The only relevant APIs were `all_variants`, `localized`, `pick_random_mutations`, and `apply_random_mutations`, with the behaviors documented above; no saturation or exhaustive-single-substitution API or example script was present.

Evidence that would have changed the score: a documented method/mode (or, for `partial`, a first-party example script) accepting a parent sequence plus a region and emitting all three non-reference nucleotide substitutions at every eligible position without caller-supplied position/substitution lists or caller aggregation. No such evidence was found.
