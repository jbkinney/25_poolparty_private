# `from_vcf` — a source operation reading variants from a VCF

Answers part of **Reviewer 3**'s comment in the BMC Bioinformatics report of
26.07.29 (`../REFEREE_REPORT_26.07.29.md`).

Sibling of `../benchmarks/` and `../comparison/`. No shared artifacts.

## What R3 asked, and what this does and does not answer

R3 raised two distinct things. They point in **opposite directions** and must not
be conflated in the response letter.

| R3's point | Direction | Answered here |
|---|---|---|
| *"Base sequences … are represented/input as sequence string (pasted by user) rather than standardized co-ordinate or HGVS notation (this limits interpretation/ingress of variants from clinical and genomics resources)"* | VCF **in** | **Yes** |
| *"if this is possible for PoolParty (the production of a VCF or a format that VEP will accept) … provide an easy means to generate VEP input format from PoolParty outputs"* | VCF **out** | **No** |

**Decision, 2026-08-20: input only.** VCF output requires composing a per-variant
genomic coordinate through arbitrary operation chains, which is `../comparison/FINDINGS.md` B2
and is genuinely undefined for `recombine`, `stack`, `join` and `shuffle`. Reading
a VCF sidesteps that entirely — coordinates flow in, not out.

The response letter must therefore concede the VEP-output gap separately rather
than presenting `from_vcf` as answering it.

## A second scoping limit

R3 names **SpliceAI, EVE, AlphaMissense and PolyPhen**. Those are not one class of
tool, and only one of the four takes a DNA window:

| Model | Named by R3 | Input | Served by `from_vcf` |
|---|---|---|---|
| SpliceAI | yes | DNA window, one-hot | **partly** — see below |
| AlphaMissense | yes | protein FASTA + MSA | no |
| EVE | yes | multiple sequence alignment | no |
| PolyPhen | yes | protein sequence + alignment | no |
| AlphaGenome | no (added for reference) | DNA interval, one-hot | **partly** |

The three protein tools need transcript-aware translation from a genomic variant to
an amino-acid change — which is exactly what VEP does, and is the part `from_vcf`
cannot reach.

**"Partly", not "yes", for the two DNA models.** SpliceAI's input depends on
gene-boundary N-padding and per-gene reverse-complement, so the same variant yields
different sequences for different overlapping genes. Both are annotation-dependent
and both are out of scope here, so `from_vcf` windows are close approximations, not
identical inputs, for minus-strand genes and variants near transcript boundaries.
Verified against each repository; see `PRIOR_ART.md`.

## Files

| Path | What it is |
|---|---|
| `PRIOR_ART.md` | How SpliceAI, AlphaGenome, AlphaMissense and EVE actually handle VCF, read from source |
| `DESIGN.md` | The proposed API, decisions taken, and decisions still open |

## Status

**Design complete, revised 2026-08-21** after three independent reviews
(correctness, adversarial, simplicity). Seven decisions were reversed and four
factual claims corrected; `DESIGN.md` records the current reasoning and the
corrections, not the superseded arguments.

Nothing implemented; nothing in `poolparty-statecounter/` modified.
