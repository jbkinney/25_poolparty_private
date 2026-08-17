# Ledidi citation-integrity audit

Audited record: `final/ledidi.md`

Repository checked at master commit `adbca708d45340fb7f375e4d0e2438d3cfa7c852`; the published PyPI 2.1.0 wheel was also checked against release commit `7adfcd6453a6851e2145bec328d5063bd19f4c4f`. The local 40-page paper PDF, live documentation and package URLs, GitHub/PyPI metadata, and Zenodo record 14604495 were checked independently.

## NOT FOUND IN ANY SOURCE

None.

## Non-verified items

1. **misquoted — Tutorial 2 quotation (record lines 136–139).** The record quotes, “we could force that edit to be included in the design.” The primary notebook actually says, “we could force that edit ot be included in the design.” The record silently corrects the source's `ot` typo while presenting the result as a verbatim quotation. The substantive meaning is unchanged.

2. **misquoted — README modality quotation (record lines 216–218).** The record ends the quotation with “inputs)”. The README says “inputs such as small molecules).” The quoted text deletes `such as small molecules` and supplies a closing parenthesis without an ellipsis.

3. **wrong-location — affinity-catalog paper quotation (record lines 240–243).** The exact sentence “an affinity catalog for GATA2 binding reveals a sophisticated usage of a learned cis-regulatory code” occurs once in the paper, in the abstract on physical PDF page 4. It is not a quotation from §3.3/Fig. 4E as the record's citation placement claims. Section 3.3 discusses the same subject, but does not contain that exact sentence.

4. **wrong-location — ReadTheDocs source inventory (record line 40).** The cited live `en/latest/` site identifies itself as the Ledidi v2.0.0 documentation. At audit time, `getting_started.html`, `input_output.html`, `parameters.html`, and `faq.html` returned HTTP 404; Tutorials 0, 7, and 8 also returned HTTP 404. Those files exist on repository master, but they are not available at the cited live documentation location. The ReadTheDocs root itself is live, so this is not a dead-link finding.

5. **number-wrong — “Nine tutorial notebooks … all rendered on readthedocs” (record lines 393–394).** Nine notebooks exist in repository master, but only Tutorials 1–6 are rendered at the cited ReadTheDocs site. Tutorials 0, 7, and 8 return HTTP 404. The live rendered count is six, not nine.

6. **wrong-location — Zenodo provenance for tutorial oracle models (record lines 43, 393–394, and 523–525).** Zenodo record 14604495 contains BPNet model/configuration files, including GATA2, E2F3, and ATAC models. Its file inventory contains no Enformer, Malinois, ChromBPNet, or Beluga model. Therefore it does not support either “the rest use BPNet/Enformer/Malinois oracles from Zenodo record 14604495” or the broader statement “Oracle models used by the tutorials are on Zenodo.” The notebooks refer to local model paths, but that does not establish this Zenodo record as the source of the absent model families.

7. **number-wrong — `.rst` page count (record lines 154–155 and 546–547).** Repository master contains 11 `.rst` files: six top-level documentation pages and five API stubs. The claimed count of seven is not the master-tree count. The release-era documentation tree has seven `.rst` files, so the evidence combines counts from one version with the master documentation/notebook inventory from another. The reported zero keyword hits remain reproducible; only the stated search-universe count is wrong.

8. **wrong-location — exhaustive characterization of repository/paper readouts (record lines 306–309).** The cited repository and paper do not support “Every readout in the repo and the paper is either an oracle prediction or an edit count/position.” They also report attribution scores/logos, FIMO motif-hit counts, correlations, and runtime benchmarks. This finding concerns only the exhaustive factual description used as evidence, not the record's consequence-annotation rating.

9. **minor-discrepancy — E2F model identity (record lines 420–422).** The paper's main text/Fig. 2 identifies E2F6, while Supplementary Table 1 identifies E2F3. Zenodo record 14604495 contains `E2F3.torch`; the tutorials request `E2F6.torch`, which is absent from that record. The sources conflict, so the asserted E2F6 identity cannot be cleanly verified across the cited evidence.

10. **minor-discrepancy — “17 curated motifs” (record lines 459–463).** The paper says the baseline used “a set of 17 motifs thought to be relevant.” It does not call them `curated`; that word does not occur in the paper text. The number 17 and the baseline procedure are verified, but `curated` is an unsupported characterization.

11. **misquoted — FAQ quotation (record lines 485–487).** The source sentence is “If a motif could have been inserted at several locations, Ledidi commits to one and all sampled sequences reflect that choice.” The record terminates the quotation after “one” and adds a period, without an ellipsis. The omitted clause materially explains the correlation claim the quotation is being used to support.

12. **wrong-location — `greedy_pruning(threshold=0)` described as unreleased (record lines 493–495 and 515–517).** The cited master changelog says threshold 0 is accepted “again,” and commit history describes the intervening behavior as a backward-incompatible regression from 2.1.0. In the published v2.1.0 implementation, the comparison `if best_score < threshold` already permits `threshold=0` behaviorally. Thus the cited source and release code do not support classifying this capability as unreleased 2.2.0 work. The other capabilities grouped with it in those sentences—new input validation, `plot_loss`, and the tangermeme-based `plot_edits`—are unreleased relative to v2.1.0.

## Uncited factual claims

None identified after excluding evaluative language and conclusions, as instructed.
