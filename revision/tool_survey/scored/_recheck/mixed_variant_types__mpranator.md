# MPRAnator — mixed variant types in one library (row 6)

- **Value:** `yes`
- **Confidence:** high
- **Tool source revision checked:** canonical repository commit `9969790d62410138d4281b7955da6d085f07b1bc` (current `master` checked 2026-08-17); tool paper: Georgakopoulos-Soares et al., *Bioinformatics* 33(1), 2017, pp. 137–138, DOI `10.1093/bioinformatics/btw584`.

## Primary-source evidence

1. One SNP-design form contains one input field labelled **“Enter your SNPs (VCF Format)”**. `SnpFile.getSnps()` parses all lines of that one field into one `snps` collection, accepts REF and ALT strings of differing lengths (including `.`), and does not impose a same-event-type rule across records.

   Sources: canonical source, [`iliasApp/forms.py`, lines 501–505 at commit `9969790`](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/iliasApp/forms.py#L501-L505), and [`mycustom.py`, lines 241–290](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/mycustom.py#L241-L290).

2. The first-party variant generator explicitly documents mixed allele-length handling: **“if the SNP contains a small deletion or insertion (smaller than 10nt) we either remove part of the sequence or we insert adenines in one edge.”** For every VCF allele it computes `difference = len(REF) - len(change)` and has separate branches for `difference < 0` (insertion), `difference == 0` (substitution), and `difference > 0` (deletion). All generated records are accumulated into the same `ExtraSequences`/`ExtraNames` result returned by the call.

   Source: canonical source, [`parseSNPs.py`, lines 29–31 and 134–206 at commit `9969790`](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/parseSNPs.py#L29-L31) and [`parseSNPs.py`, lines 134–206](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/parseSNPs.py#L134-L206).

3. The request handler sends that one parsed SNP collection through one `make_sequence_copies(...)` call, adds every returned sequence to one `finalOutput`, passes that pool once to `createMPRAResultOutput(...)`, and writes one `mpraOutput` response / `MPRA_SNP_result.txt`. The first-party documentation describes the result as **“the synthesized oligonucleotides in FASTA format.”**

   Sources: canonical source, [`iliasApp/views.py`, lines 351–386 at commit `9969790`](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/iliasApp/views.py#L351-L386), and first-party documentation template, [`iliasApp/templates/iliasApp/docs.html`, lines 40–50](https://github.com/hemberg-lab/MPRAnator/blob/9969790d62410138d4281b7955da6d085f07b1bc/iliasApp/templates/iliasApp/docs.html#L40-L50).

## Reasoning against the operational definition

MPRAnator's SNP-design workflow supports three structurally different variant components—substitutions, insertions, and deletions—declared together as records in one VCF-format specification and emitted together as one FASTA pool. This exceeds row 6's `partial` boundary of only two types or variants of the same kind. Indels are restricted to less than 10 nt (`abs(difference) < 10`), but all three structural types remain supported in the same run/output, so the operational value is `yes`.

## Disconfirmation attempt

I searched and inspected the complete first-party repository at the pinned commit, independently of DNA Chisel: all forms, request handlers, SNP parsing/generation code, output formatting, documentation templates, downloadable API scripts, Transmutation, motif design, and tests. I specifically checked whether insertion/deletion support required a separate endpoint or separate result, whether validation forced homogeneous REF/ALT lengths, and whether the outputs were left for the caller to merge. It does not: one `SnpS` field is parsed together, the generator handles all three `difference` cases in the same call, and the view builds one `finalOutput` and one FASTA response. I also checked Transmutation's apparently mixed checkboxes; those transformations are sequentially applied to every input sequence rather than pooled as distinct component classes, so they were not used as evidence for the `yes` score.

Evidence that would have changed the value: source showing that VCF records must all have the same REF/ALT length class, that insertions/deletions are processed through separate runs/files, or that the caller must concatenate results. The pinned implementation shows the opposite.
