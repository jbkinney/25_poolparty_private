# MPRA Design Tools — mixed variant types in one library (row 6)

- **Value:** `yes`
- **Confidence:** high
- **Tool version / source revision checked:** `mpradesigntools` 0.2.0, canonical repository commit `afd386ef12051bb0a37ad63a6f92acd555246584` (current `master` checked 2026-08-17).

## Primary-source evidence

1. The first-party README identifies one library-design entry point and one file specification: **“Currently the main function of MPRA Design Tools package is to design a set of barcoded sequences for MPRA experiments ... This is done with the `processVCF` function.”** Under the input constraints for that same VCF, it says: **“Insertions and deletions must encode the reference and alternate alleles (respectively) as a dash character '-'.”** Thus the documented input format for the single `processVCF` operation includes indels, rather than requiring a separate indel-design run.

   Source: canonical repository README, [`README.md`, lines 45–57 at commit `afd386e`](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/README.md#L45-L57).

2. The implementation establishes three structurally different supported component types on every input row: **`isSNV = ...`**, **`isINS = snp$REF == '-'`**, and **`isDEL = snp$ALT == '-'`**. The following comment states that the function has separate generation blocks for **“each of these cases”** and that there are **“subtle differences in each of the types of variants.”** The dispatch implements an SNV branch, an insertion branch (`else if (isINS)`), and a deletion branch (`else if (isDEL)`); inputs outside those three are reported as **“Not identifiable as SNV, insertion, or deletion.”**

   Sources: canonical source, [`R/processVCFfast.R`, lines 232–242](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L232-L242), [lines 544–548](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L544-L548), and [lines 764–775 and 967–975](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L764-L775) (the final quoted failure text is at [lines 967–975](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L967-L975)), all at commit `afd386e`.

3. Pooling is tool-provided within that one call. `processVCF` reads one VCF, invokes `processSnp` **rowwise** over it, then collects all successful per-row outputs using **`Reduce('rbind', .)`** into `successes`. With `outPath`, it performs **`write_tsv(successes, path = outPath)`**. This is one combined output produced by the tool, not separate runs or caller-side merging.

   Source: canonical source, [`R/processVCFfast.R`, lines 1099–1138](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L1099-L1138) and [lines 1221–1258](https://github.com/andrewGhazi/mpradesigntools/blob/afd386ef12051bb0a37ad63a6f92acd555246584/R/processVCFfast.R#L1221-L1258), both at commit `afd386e`.

## Reasoning against the operational definition

One VCF is one specification to the tool's documented `processVCF` operation. That VCF can declare SNVs, insertions, and deletions in different rows. The function determines the type row by row, generates constructs through three distinct type-specific branches, and then itself row-binds the constructs into one result and, when requested, one TSV. Therefore MPRA Design Tools supplies mixed-type pooling directly and satisfies `yes`. It is not the `partial` case described by row 6 because support is not limited to two types or to variants of the same structural kind: the implementation supports three structurally distinct kinds (substitution, insertion, deletion).

## Disconfirmation attempt

I searched the complete first-party `mpradesigntools` 0.2.0 repository at the pinned commit (all `R/`, `man/`, and README files) for `variant`, `allele`, `SNV`, `SNP`, `insert`, `delet`, `indel`, `processVCF`, `outPath`, `write_tsv`, and `rbind`. I specifically checked whether the VCF was validated as a uniform variant type, whether SNVs and indels required separate public operations or separate calls, whether any type branch wrote its own output, and whether combining results was left to the caller. No uniform-type restriction or separate-run workflow exists: `processVCF` dispatches each row independently among SNV/insertion/deletion and performs the final `rbind` and optional TSV write itself.

Evidence that would have changed the score: a first-party restriction requiring all records in one VCF to share a type, separate entry points/output files that users must merge, or support for only two component types (which would make the score `partial`). None was found; the source instead directly demonstrates three-type pooling in one invocation and output.
