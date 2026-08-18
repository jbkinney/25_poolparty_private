# Capability row definitions

Operational tests for scoring eight DNA library-design tools. These are tests, not
descriptions: apply them as written rather than substituting your own reading of
the terms.

**Scoring vocabulary.** `yes` = supported · `partial` = supported with a stated
restriction · `no` = not supported · `unknown` = genuinely undeterminable from
available sources, never a hedge.

Apply one threshold across all tools. If two tools present equivalent evidence they
receive the same value.

## Global rule — binding on every row below

**A capability counts only where the tool provides an operation, parameter or mode
for it.** It does not count where the user reconstructs it by composing primitives
the tool provides for other purposes, plus their own bookkeeping.

`partial` describes a **tool-provided capability that is restricted** — restricted
in scope, in the types it accepts, in requiring one component to be supplied
pre-built, or in existing only in an example script rather than the documented API.
**`partial` never describes a capability the user assembles.** If the user must
write the logic that delivers the capability, the value is `no`.

Two clarifications, because the line is not "did the user write code" — these are
libraries, so the user always writes code:

- Feeding a tool's **own documented output back into its own documented input** is
  tool-provided composition, and may be `partial`. Both ends are interfaces the
  tool ships for that purpose.
- Calling several unrelated functions and reconciling their outputs yourself — for
  example normalising incompatible shapes, or appending results to a list across
  repeated calls — is user reconstruction, and is `no`.

---

### 1. Library object
A documented API call returns a value representing the whole multi-sequence library
as an instance of a type the tool defines for that purpose (carrying member
identity, design provenance, or membership operations), which the user can bind to
a variable and inspect — as opposed to the tool only writing files, or handing back
a general-purpose container (ndarray, DataFrame, ragged array, torch tensor,
list/dict of sequences) whose library meaning lives only in the caller's head.
`partial` = a durable, documented, in-memory value representing the **whole
library** is returned by the tool's own design API, but its type is general-purpose
(DataFrame, list, tensor).
An object representing a single sequence or a single optimization problem is `no`,
however purpose-built its type.

### 2. Composable operations
Two or more distinct **library-design** operations share a carrier type and compose
in either order, with fan-out (one intermediate feeding two downstream operations)
expressible. General sequence utilities (parse, revcomp, translate, pad, one-hot,
write) and cloning-simulation primitives are **not** design operations.
`partial` = design operations exist and chain, but only in a tool-fixed order, only
for a small documented subset, or only by feeding the tool's own documented output
back into its own documented input (which counts, per the global rule, because both
ends are tool interfaces).
`no` = the tool fixes the order of its design operations in source, **or** the tool
exposes no library-design operations at all; distinguish which in your evidence.
Internal caching or memoisation of shared intermediates is **excluded** from the
threshold — this row measures expressibility only.

### 3. Saturation mutagenesis
From one specification, the tool enumerates every single substitution at every
eligible position of a region. `partial` = exhaustive over positions but only for
one narrow event type, or the capability exists only in an example script rather
than the documented API.
**If the user must supply the list of positions or substitutions, the tool is not
performing the enumeration — score `no`.**

### 4. Randomly sampled variants
A documented rate, count or RNG parameter produces a *sample* of the variant space.
**Random sequence generation is not sampled mutagenesis** — generating random
barcodes, shuffled controls, or stochastic search during optimization does not
count.

### 5. Pairwise and higher-order variants
Two or more **independent** mutation events co-occurring in a single output
sequence, with the combinations enumerated or sampled **by the tool**.

Explicitly excluded: hand-authored multi-variant input (a user-written VCF);
uniform context edits applied to every oligo (e.g. PAM-protection edits); one
multi-base single event (a 3 bp deletion is one event); combinations the user
produces by running the tool repeatedly and merging.

**Boundary rules — binding.** This row sits beside two others describing adjacent
capabilities. Without these rules the same capability is credited twice, once here
and once in its own row.

- **vs. row 9 (Combinatorial multi-motif placement): inserting multiple motifs is
  *placement*, not mutagenesis.** A tool that positions two or more motifs across a
  background sequence scores in row 9 and **not** here, however many combinations
  it produces. Score here only if the tool combines *mutation events* —
  substitutions, insertions, deletions applied to a parent sequence.
- **vs. degenerate templates: IUPAC expansion does not count.** Expanding
  `ATG NNK GGC` into concrete sequences yields multi-position variation, but the
  user authored the combinatorial specification and the tool merely decompressed
  it. The tool must generate the combinations from a design specification, not
  enumerate a template the user already wrote.
- **vs. row 12 (Model-guided variants): objective-driven optimisers score `partial`
  at most.** A tool that emits a sequence carrying several edits because an
  optimiser converged there has not enumerated or sampled a *declarable*
  combinatorial space. Score such tools `partial` here, and give full credit in
  row 12.

### 6. Mixed variant types in one library
Two or more structurally different component types declared in **one**
specification and emitted as **one** pooled output. `partial` = the tool emits one
pooled output containing two or more component types, but the types are restricted
(limited to two, or to variants of the same kind), or one component must be
supplied pre-built.
**Separate runs that the user merges are `no`** — the pooling must be performed by
the tool.

### 7. Codon / amino-acid substitutions
Substitution specified at the codon or amino-acid level, with the tool handling the
codon-to-residue mapping and offering a choice of replacement policy (e.g. all
residues, most-frequent codon, degenerate codon). Nucleotide substitution that
merely happens to fall inside a coding region does **not** count — the tool must be
codon-aware. `partial` = codon awareness exists but the substitution set is fixed.
**If the user must supply the residue list, the tool is not designing the
substitutions — score `no`.**

### 8. Insertions and deletions
Generation of variants that insert or remove bases as a **designed variant type**,
via a supported operation. Supplying pre-made indel sequences as input does not
count. Gap-marked deletions count provided the tool either shortens the sequence or
explicitly marks the deletion in the output. `partial` = one of insertion or
deletion only, or fixed length only.

### 9. Combinatorial multi-motif placement
Two or more distinct motifs placed across multiple positions and/or orientations,
with the combinations enumerated or sampled by the tool. `partial` = single-motif
placement, or positions enumerated but orientations fixed.

### 10. Recombination
Construction of variants by joining segments from two or more parent sequences at
breakpoints the tool chooses or enumerates. Concatenating user-supplied fragments
does not count; the tool must generate the breakpoint combinations. This covers
both the wet-lab sense (chimeric constructs) and the in-silico sense (swapping
segments between parents to ask which region drives a model's prediction).

### 11. Constraint-based sequence optimization
The tool **modifies** a sequence so it satisfies declared constraints, rather than
only rejecting sequences that violate them. **Rejection-only filtering does not
count** — this row is about the tool altering the sequence to comply. `partial` =
optimization against a narrow constraint class only.

*Displayed in the **Physical construction** block.*

### 12. Model-guided variants
Design driven by the output of a predictive model, where the model's prediction
determines the edit. `partial` = the tool accepts an arbitrary scoring callable so
a model can be attached, but performs no optimization against it. `yes` = an
optimization loop against the model's output.

*Displayed in the variant-generation block. This row is about **biological
function**, not manufacturability — it does not belong with the physical-construction
rows even though all three involve an objective.*

### 13. Per-sequence design cards
Structured, machine-readable metadata attached to each output sequence recording
**how it was constructed**, beyond naming the mutation it carries. `partial` = some
provenance recorded, but not per sequence, not structured, or limited to the
variant identity.

### 14. Automatic naming
Informative identifiers generated **by the tool** for each output sequence,
composed from the construction history. **Carrying a user-supplied identifier
through to the output does not count.** `partial` = names generated but not
composed from construction history, or generated for only some operations.

### 15. Sequence styling
Visual marking of **what changed inside the sequence itself** — highlighting,
casing, colour — so a user can audit a variant by eye. Plots, graph views and
report documents are **not** this row.

### 16. Genome coordinates
Accepts reference-genome coordinates as input, and/or emits coordinates that locate
each designed sequence in a reference. `partial` = coordinates accepted for
extraction only, with no coordinate output for designed variants.

### 17. Transcript / annotation aware
Represents transcript models, exon/intron structure, or annotation files (GTF/GFF),
such that design respects them. `partial` = annotation consumed for one purpose
only, or exon structure handled without transcript models.

### 18. Shuffling
The tool generates shuffled variants of an input sequence, performing the shuffle
itself. A user shuffling sequences with their own code does not count.

`yes` = **dinucleotide-preserving** shuffling is available.
`partial` = mononucleotide shuffling or scrambling only.
`no` = no shuffling operation.

The dinucleotide threshold is deliberate: dinucleotide composition alone drives much
of a sequence model's output, so a mononucleotide-shuffled control systematically
overstates motif effects. A tool offering only a scramble is genuinely weaker for
this purpose and the row should record that.

**Boundary — vs. row 9 (Combinatorial multi-motif placement):** reverse-complement
and orientation flipping are **not** shuffling. Orientation is scored in row 9,
which tests positions *and orientations*. Score here only for operations that
permute the residues of a sequence.

*Displayed in the variant-generation block of the table, alongside insertions and
deletions, recombination and codon substitutions.*

### 19. Synthesis-constraint checking
The tool checks designed sequences against physical synthesis or cloning
constraints — oligo length caps, homopolymer runs, restriction-site occurrences,
extreme GC, melting temperature, secondary structure, repeat content, or
cross-hybridisation and off-target matching against a supplied background reference
— and reports or filters the violations.

Screening designed sequences against an external reference (a host genome, a
vector, or another oligo set) counts here; it is the relational case of the same
capability, not a separate row.

`yes` = several such constraint types are modelled by the tool as named,
documented checks.
`partial` = a few constraint types are provided as named checks, **or** the tool
supplies only a generic predicate mechanism through which a user can express such a
check without the constraint types being modelled.
`no` = no constraint checking against physical synthesis or cloning limits.

*Displayed in the **Physical construction** block.*

Count only constraints on **physical realisability**. Checks on biological content —
"does this contain a stop codon", "is this in frame" — are not this row.

**Boundary — vs. row 11 (Constraint-based sequence optimization):** these two are
complementary and must not both be credited for the same mechanism. Row 19 is
**detecting** violations — checking, reporting, rejecting, filtering. Row 11 is
**fixing** them — altering the sequence so it complies. A tool that only rejects
scores in row 19 and `no` in row 11. A tool that repairs scores in row 11, and in
row 19 only if it also exposes the checks independently.

### 20. Primer design
Design of PCR or mutagenic primers for constructing the library in the laboratory —
primer sequences, melting temperatures, overlaps, or mixing ratios.

`yes` = primers are a documented output of the tool.
`partial` = primer-adjacent constraints are handled (e.g. avoiding heterodimerisation
against a supplied primer set) without primers being designed or emitted.
`no` = no primer output and no primer-specific handling.

**Boundary — vs. rows 11 and 19:** this row is about producing the oligos used to
*build* the library, not about checking or optimising the designed sequences
themselves.

*Displayed in the **Physical construction** block.*
