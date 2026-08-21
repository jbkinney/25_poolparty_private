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

R3 names SpliceAI, EVE, AlphaMissense, PolyPhen. Those are not one class of tool:

| Model | Input | Served by `from_vcf` |
|---|---|---|
| SpliceAI | DNA window, one-hot | **yes** |
| AlphaGenome | DNA interval, one-hot | **yes** |
| AlphaMissense | protein FASTA + MSA | no |
| EVE | multiple sequence alignment | no |

The protein/MSA half needs transcript-aware translation — which is what VEP does,
and is the part `from_vcf` cannot reach. Verified against each repository; see
`PRIOR_ART.md`.

## Files

| Path | What it is |
|---|---|
| `PRIOR_ART.md` | How SpliceAI, AlphaGenome, AlphaMissense and EVE actually handle VCF, read from source |
| `DESIGN.md` | The proposed API, decisions taken, and decisions still open |

## Status

Nothing implemented. Nothing in `poolparty-statecounter/` has been modified. This
directory holds a design under discussion.

The two architectural questions are settled: state generation is **eager**, and VCF
parsing uses the **standard library**, adding no dependency. Both are argued with
measurements in `DESIGN.md`.

All design questions are settled. What remains is implementation, and confirming
three low-stakes placement choices noted at the end of `DESIGN.md`.
