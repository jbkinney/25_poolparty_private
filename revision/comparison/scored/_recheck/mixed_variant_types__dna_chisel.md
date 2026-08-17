# DNA Chisel — mixed variant types in one library (row 6)

- **Value:** `no`
- **Confidence:** high
- **Tool version / source revision checked:** DNA Chisel 3.2.16, canonical repository commit `68c09304341c3656f3dfe63eda37757d6a7b3917` (current `master` checked 2026-08-17).

## Primary-source evidence

1. The core documented object is singular: **“Problem specifications: sequence, constraints, optimization objectives.”** Its `sequence` parameter is **“A string of ATGC characters,”** and initialization stores it as `self.sequence`. Multiple constraints/objectives therefore act on one sequence-design problem; they do not declare different library component types.

   Sources: canonical source, [`dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py`, lines 20–71 at commit `68c0930`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py#L20-L71) and [lines 115–139](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/dnachisel/DnaOptimizationProblem/DnaOptimizationProblem.py#L115-L139).

2. The first-party usage documentation labels retrieval **“GET THE FINAL SEQUENCE”** and exposes `final_sequence = problem.sequence` or one annotated `final_record`. It documents no pooled multi-sequence result from a specification.

   Source: canonical README, [`README.rst`, lines 69–82 at commit `68c0930`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/README.rst#L69-L82).

3. The repository's own collection example explicitly states: **“DNA Chisel is not originally meant for creating collections of sequences”** and **“Here we create 20 short sequences one after the other.”** Its main loop calls a single-sequence design function repeatedly and performs caller-side pooling with `existing_primers.append(new_primer)`.

   Sources: first-party example, [`examples/common_scenarios/primers_collection.py`, lines 1–8 at commit `68c0930`](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/examples/common_scenarios/primers_collection.py#L1-L8) and [lines 39–46](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/68c09304341c3656f3dfe63eda37757d6a7b3917/examples/common_scenarios/primers_collection.py#L39-L46).

## Reasoning against the operational definition

Row 6 requires two or more structurally different component types declared in one specification and emitted as one tool-pooled output. DNA Chisel's documented design unit and result are one sequence. Even its first-party multi-sequence example performs separate design calls and appends each result to a caller-owned list. Under the binding global rule and row-6 statement that separate runs merged by the user are `no`, this is `no`, not `partial`.

## Disconfirmation attempt

I searched the complete first-party repository at the pinned commit—including `dnachisel/`, `docs/`, `examples/`, and tests—for `mixed`, `pool`, `pooled`, `library`, `collection`, `variant`, `sequences`, `output`, and multi-sequence design interfaces. I inspected the complete `DnaOptimizationProblem` constructor/output path, `MutationSpace`, report APIs, CLI documentation, and the one example explicitly concerning a collection. Multi-sequence inputs found in BLAST/Bowtie utilities and constraint-report helpers are comparison/reporting inputs, not a design specification that generates and pools component types. No tool-provided mixed-library operation, parameter, or mode was found.

Evidence that would have changed the value: a documented API/mode accepting two or more structural design components (for example substitutions plus indels or designed variants plus controls) in one specification and returning/writing them as one pooled output. No such evidence was found.
