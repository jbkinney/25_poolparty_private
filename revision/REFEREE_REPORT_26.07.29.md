# Referee report — BMC Bioinformatics, 26.07.29

**Verbatim, as supplied by the author.** Pasted into a working session on
2026-08-15 and recovered from that session's transcript on 2026-08-17. Nothing
below is paraphrased, reordered or abridged.

This file exists because the report was not in the repository: every planning
document in `comparison/` cites it (R1 #3/#5, R2 2b, R3, editor) and none of them
could be checked against it. One consequence had already reached
`MAIN_TEX_CHANGES.md` -- see the note at the end of this file.

**Do not edit.** Corrections and readings belong in `comparison/PLAN.md` and
`comparison/MAIN_TEX_CHANGES.md`, not here.

---

25_poolparty referee report
26.07.29
Feedback from the Editor
Dear Dr Kinney,

Your manuscript, "PoolParty: streamlined design of DNA sequence libraries in Python", has now been assessed.

We invite you to revise your paper, carefully addressing the comments from the reviewers and the editor. Please ensure the results are accurately reported, any overstated conclusions are rewritten and the limitations of the work fully explained. When your revision is ready, please submit the updated manuscript and a point-by-point response. This will help us move to a swift decision.

Editor Comments

"In general, the reviewers had positive impressions of PoolParty, and it can likely be revised to a state that will be acceptable. Please carefully address each point raised by the reviewers. It's especially important to provide benchmarking data, comparisons to existing tools, and not to overstate the benefits the tool may provide or fail to discuss limitations. Reviewer 1 recommends comparing PoolParty to 2-3 other comparable tools. However, rather than arbitrarily picking two or three other tools for comparisons, you should use a well-defined process for picking which tools, and how many, to include in comparisons. Perhaps 2-3 is adequate; however, that may be too few tools for comparisons. The comparisons should include, at a minimum, benchmarking (e.g., runtime, memory usage), program features, and similarities/differences in program outputs."
-Perry Ridge

In reviewing the referee reports, we noticed that you were encouraged to add specific references to your manuscript. While it is not required to include these specific references, you may also choose not to include any if you feel it won’t enhance your manuscript. You are also welcome to seek alternative publications that are relevant to your manuscript's content if you feel that those proposed are not suitable.

We recommend submitting all revisions within the mentioned deadline.

If you need more time, please contact us and include your submission ID.

Kind regards,

Payel Karmakar
Editor
BMC Bioinformatics

Feedback from the Reviewer

Reviewer 1

This study presents PoolParty, a Python package for designing DNA sequence libraries. Its core design employs a directed acyclic graph (DAG) architecture, representing libraries as Pool objects and sequence generation steps as Operation objects. The manuscript could be strengthened in benchmarking data, comparison with existing tools, and empirical support for design card utility. These are supplementary rather than fundamental issues, and should be readily addressable in revision.
1. In the GB1 example, the library exceeds 540,000 sequences, yet the paper provides no runtime or memory consumption benchmarks. Readers cannot assess the tool's practical usability across different library scales. The authors should supplement benchmarking data, including runtime and peak memory usage across libraries of increasing sizes (e.g., 1K, 10K, 100K, 500K sequences), and specify the test environment configuration.[LL1.1]
2. The paper highlights design cards as a core feature, emphasizing their use as covariates in downstream analyses. However, it does not evaluate whether design card variables provide significant incremental predictive value beyond sequence-derived features. The authors should either supplement a comparative experiment quantifying the added information gain of design cards, or explicitly state in the Discussion that this validation is planned as future work.[LL2.1]
3. The introduction mentions VaLiAnT and MPRAnator, noting their limitations, but does not provide a systematic comparison between PoolParty and these tools. The authors should provide a table comparing PoolParty with 2-3 existing tools on core features, and demonstrate code conciseness in at least one common use case, to help readers objectively assess PoolParty's advantages.[LL3.1][LL3.2]
4. Provide benchmarking data on runtime and memory usage across libraries of varying sizes.[LL4.1]
5. Include a table comparing features with 2-3 existing tools.[LL5.1]
6. To further contextualize this research within the existing academic landscape and strengthen its theoretical support, it is recommended that the authors supplement several relevant representative references.[LL6.1][LL6.2]
[1] https://doi.org/10.1093/nsr/nwaf495[LL7.1]
[2] https://doi.org/10.1038/s41467-026-72882-y[LL8.1]
[3] https://doi.org/10.1186/s13059-025-03606-6[LL9.1]
[4] https://doi.org/10.1016/j.patcog.2025.112078[LL10.1]
7. Either provide empirical analysis of design cards' incremental contribution to surrogate modeling performance, or explicitly discuss this as future work.[LL11.1]

Reviewer 2

In “PoolParty: streamlined design of DNA sequence libraries in Python”, Liu et al. describe a tool to design complex oligo pool sequences to ease the design of DMS, MPRA, and other multiplex libraries. Users specify what permutations are to be made, and the tool operates in two steps to generate mutated sequences: first, a backward pass is to assign each pool a state that describes which sequence to generate, and second a forward pass generates each sequence. The resulting pool can be visualized in a variety of ways, including being exported to a pandas dataframe, a file, or visualized using highlighted characters. Over 50 operations are implemented, as well as other methods to combine and sample sequence pools. The toolkit is implanted in python and available on github and PyPI. The accompanying documentation is thorough, and provided examples demonstrate its utility. Overall, the manuscript is well-written. In my opinion, this tool will improve reproducibility and the ability to design complex sequences of pools.
Major comments:
1a) Additional statistical readouts describing pool generation would be useful. For example, how many unique sequences vs duplicates were produced in each pool (i.e. out of the universe of all permutations that meet the specifications)? How unique are the sequences (min/max/average hamming distance)? How frequent are homopolymers? Etc.[LL12.1][LL12.2]
1b) The authors disclose that LLMs were used to develop, debug, and document the code. Did the authors also write a test harness, and what edge cases were nominated by the LLMs for testing vs what edge cases were nominated by the human developers for testing? [LL13.1][LL13.2]A
Minor comments:
2a) The backward/forward pass steps are not well described in the manuscript, and figure 1b is hard to interpret. I suggest adding a supplemental figure showing how these passes are performed. For example, the text explains that the backward pass happens as information is passed to the root node – does this mean that the state is passed from ‘Mutagenize’ to ‘Pool A’ in or first from ‘Pool Final’ to ‘Stack’? And if so, what information is passed? Perhaps showing each step in an expanded figure or labeling arrows with the order would be more useful.[LL14.1]
2b) A comparison to the pools produced by other tools such as VaLiAnT or MPRAnator would provide additional confidence in the ability of PoolParty to produce sequence pools as expected. Metrics such as the number of pool elements that overlap and an investigation of any unique pool elements would ensure that the generate pools are accurate.[LL15.1]
2c) additional operations to improve library content would be useful. E.g. to maximize hamming distance, remove sequences containing restriction enzyme sites, remove homopolymers, etc. would make the tool more useful in generating ready-to-use pools.[LL16.1][LL16.2]
2d) The name PoolParty may be confused with the existing PoolParty analysis softwre (https://doi.org/10.1111/1755-0998.12784)[LL17.1]

Reviewer 3

This paper is well written, well organised/clear and the figures beautifully presented. I was able to clone the Github repo, install both packages and run the tests, and tutorial examples and some additional operations with ease. The documentation that accompanies the software is excellent and follows standards for the field. The highlighting feature is useful for sanity checking sequence generation on example sequences.

A systematic and more intuitive approach to multiplexed oligo library generation is needed and PoolParty fills a gap not available with current tools, which are either bespoke or require extensive upstream thought/planning rather than a more cut/paste approach.

Overall, PoolParty is a powerful sequence library design tool, but it is not a genome-/transcript-aware or clinically oriented multiplexed variant library design platform;[LL18.1] the authors don’t claim that it is, but there is tacit, and sometimes explicit claims that this tool might replace other tools with an equivalent capability of operations; this should be addressed. The lack of a clear path (as far as I can see) to generating more extensive variant annotations (e.g VEP) is a limitation that may limit the tool being built into workflows with other tools[LL19.1], but I do think the novel features of PoolParty are strong (especially for use in in silico experiments [as the authors demonstrate], and ML); these could be brought more to the fore in the publication text.

I recommend accepting this manuscript for publication in Bioinformatics. I have some suggestions, mostly around language/emphasis.

Comments:

• I think that the emphasis made on PoolParty being a universal tool that can ‘replace’ (line 256) existing code/tools needs to be qualified more. The benefits of PoolParty are that it is straightforward to use, interpretable, modular, with a nice ‘design card’ tracking and graphed interpretation – I think this is particularly useful for ML and in silico work.

However, whist it does streamline design processes, the package operates on input sequences rather than genomic loci co-ordinates/transcript so in terms of wet-lab experiments it is most useful for fundamental biology cDNA DMS and MPRA assays that are episomal (i.e. from plasmid / using reporters) rather than MAVEs that rely upon genome editing (which is becoming the gold standard for clinical variant work). It does not intuitively and natively support reference genome coordinates, transcript models, exon/intron structures, alternative isoforms, or genome assembly awareness (as far as I can tell); VaLiAnT does do this. Base sequences into which variants are installed are represented/input as sequence string (pasted by user) rather than standardized co-ordinate or HGVS notation (this limits interpretation/ingress of variants from clinical and genomics resources). If PoolParty is presented as intuitive and adaptable rather than a universal replacement to alternative tools and code this is less of a concern and may allow researchers to more readily select which software to invest time into learning/using.

• Molecular consequence beyond codon change for missense is not reported as far as I can tell, VaLiAnT does this to some extent (e.g. stop-gained, synonymous, inframe deletion, etc, PoolParty may do this also, but I couldn’t see this after running a ‘replacement_scan’?) but generally researchers use VEP on their tool outputs to achieve extensive annotation. VEP requires standard annotation or format (e.g. VCF) for input rather than sequence strings – if this is possible for PoolParty (the production of a VCF or a format that VEP will accept) then one could obtain extensive annotations (e.g. SpliceAI predictions, EVE, AlphaMissense, PolyPhen, Consequence, genomic coordinates, Clinvar and gnomAD identifiers etc etc) which I think would greatly increase the uptake of the tool by researchers in the MAVE community that are more clinically/genomically minded – if this is possible I would state this/provide an easy means to generate VEP input format from PoolParty outputs.

• The authors note limitations at the end of the discussion. I think updating this section to address that PoolParty is mostly of benefit for DMS and MPRA rather than genome editing MAVE assays (e.g. not as useful for Saturation Genome Editing HDR libraries and Saturation Prime Editing pegRNA libraries) which require extensive and specific annotation/standardized nomenclature and consideration of transcript (generally through ingress of GTF/GFF files) and exon/intron boundaries (including split codons between exons, again through GTF/GFF boundary annotations).

Alternatively/additionally the authors might wish to highlight that the tool is agnostic to species/genome/transcript and present it as a more intuitive sequence generation tool which is particularly useful for ML testing/benchmarking which is well tracked with the design cards and interpretable/iteratively mutable through graph interrogation. This builds on the comment that generally annotations are variant specific (line 56-58), with the emphasis being that PoolParty serves as a useful tool to manifest oligo libraries as structured objects (less of a concern for wet-lab investigations in my opinion, but much more salient for ML based studies on raw sequence out of genomic context).

• Line 12 – there is not a “lack” of purpose-built software, replace with “scarcity” perhaps.[LL20.1]

• Line 51-53 - stated that VaLinAnT is assay specific, then say used for DMS and saturation genome editing two different assays (SGE is a type of DMS, but VaLinAnT is used for SGE and cDNA DMS which are quite different assays).

• Line 66 – “expensive computation” this is only true for truly massive sequence generation – most full gene DMS/saturation libraries of >20,000 variants take <5 minutes on a standard laptop CPU (true of Valiant) – this statement would be true for many permutations of sequence for use in ML (e.g. pairwise, full genome CDS on GPU) – I would qualify this point on cost.

• The use cases of protein coding, MPRA promoter and splice are good examples (line 68).

• Figure 2 B and Figure 3 B – is it necessary helpful to put these code blocks here? I’m not sure how useful this is for readership.[LL21.1]

• Line 194 – ‘most abundant codon’ – this needs to be defined, is this codon usage (I think not) or instance in the input sequence?[LL22.1]

• Line 235 – SpliceAI is defined here for the first time but is introduced before this in section 2.8 – I would move the ‘a deep learning model of mRNA splicing’ to the earlier section.

• Lien 238 “varying strength” – maxentscan in legend, I would define this is the amin body text too.

• Lien 248 – Fig.4G is introduced before 4E-F, perhaps generally introduce Fig.4. then this could stay as is.

• Line 256 – “replaces the ad hoc” is too bold and a little presumptive as for major points above– needs clarification/qualification. Also, Line 274 – “in contrast”- I think these claims are too strong – valiant can do DMS and MPRA as well as SGE/SPE – I would tone these statements down and/or emphasise what is unique to PoolParty.

• Line 271 – “in the spliceAI example” – the one in section 2.8 or 3.3? (3.3 I imagine, but best to state).

• The code is very well written and documented, and it is noted by the authors that Anthropic models have been used to do this/help which I think is sensible, editorially this may be something to confirm is acceptable at Bioinformatics (as reviewers were requested not to use these tools, but this is likely due to confidentiality).[LL23.1][LL23.2]

• The name ‘PoolParty’ is great.[LL24.1]

---

## Note on one reading already recorded elsewhere

`comparison/MAIN_TEX_CHANGES.md` D3 currently instructs: *"Do not concede HGVS.
R3's premise is wrong -- `hgvs_input` is `no` for all 13 tools surveyed, VaLiAnT
included."*

Against the text above, R3's "VaLiAnT does do this" attaches to the preceding
list -- reference genome coordinates, transcript models, exon/intron structures,
alternative isoforms, genome assembly awareness -- all of which VaLiAnT does
provide. The HGVS clause that follows is a separate statement about **PoolParty
only**, and it is correct: PoolParty accepts sequence strings, not coordinates or
HGVS notation.

R3 does not claim VaLiAnT accepts HGVS input. That instruction rebuts a claim the
report does not make, and is flagged here for resolution in
`MAIN_TEX_CHANGES.md`. This file records the observation; the fix belongs there.
