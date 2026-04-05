# Plan: Web Search for Additional Papers Relevant to PoolParty's Novelty

## Goal

Identify papers not currently cited in the manuscript that could be relevant to PoolParty's novelty claims. Focus on the Introduction's three categories of related work, plus any tools that might blur the boundaries.

## PoolParty's Core Novelty Claims (to search against)

1. First unified, declarative, composable framework for DNA sequence library design
2. DAG-based design with lazy evaluation
3. Design cards providing structured metadata for each sequence
4. Automatic sequence naming and styling
5. Applicable across assay types (MPRA, DMS, in silico)
6. Over 20 built-in Operations; extensible via subclassing

## Search Strategy

### Round 1: Direct competitors — oligo pool / variant library design software

These are the most important to find. Any tool that designs DNA sequence libraries programmatically could be a direct competitor.

**Searches:**
- `"oligo pool design" software tool Python`
- `"variant library design" software tool`
- `"oligonucleotide library design" software`
- `"sequence library design" Python package`
- `"DNA library design" computational tool`
- `"mutant library design" software`

### Round 2: DMS / saturation mutagenesis library design tools

VaLiAnT was the closest prior art. Are there others?

**Searches:**
- `"deep mutational scanning" library design software tool`
- `"saturation mutagenesis" library design software`
- `"saturation genome editing" oligo design tool` (to find VaLiAnT competitors)
- `"DMS library" design software`
- `Mutation Maker oligo design protein engineering` (flagged from CodonGenie analysis)

### Round 3: MPRA library design tools

MPRAnator and MPRA Design Tools were the only two found. Are there others?

**Searches:**
- `"MPRA" library design software tool`
- `"massively parallel reporter assay" oligo design`
- `"reporter assay" sequence library design software`

### Round 4: In silico experiment / genomic AI library design

PoolParty claims utility for designing libraries for in silico probing of genomic AI models. Are there tools specifically for this?

**Searches:**
- `"in silico mutagenesis" library design genomic deep learning`
- `"sequence design" probing genomic AI model`
- `Schreiber 2025 "programmatic design" cis-regulatory elements` (flagged from tangermeme analysis)

### Round 5: Composable / declarative / DAG-based sequence design

PoolParty's architectural novelty (DAG of composable Operations) is distinctive. Search for anything with a similar design philosophy applied to sequences.

**Searches:**
- `"composable" DNA sequence design framework`
- `"directed acyclic graph" sequence library`
- `"declarative" DNA library design`

### Round 6: General review articles on MAVE / library design

Review articles may cite tools we've missed and help ensure comprehensive coverage.

**Searches:**
- `"multiplex assay of variant effect" review 2023 2024 2025`
- `"oligo pool" design review`
- `"variant effect" library design computational review`
- `"deep mutational scanning" review 2023 2024 2025`
- `"deep mutational scanning" library construction methods review`

## Evaluation Criteria

For each paper found, I will ask:
1. Does it design DNA sequence libraries (not just analyze them)?
2. Is it a software tool (not just a wet-lab protocol)?
3. Does it overlap with any of PoolParty's novelty claims?
4. Is it already cited in the manuscript?
5. Should it be cited, and if so, where?

## Output

A summary document listing all papers found, with a brief assessment of each and a recommendation (cite / investigate further / not relevant). Any papers that warrant full analysis will be flagged for the same treatment as the 8 PDFs already analyzed.

## Potential Limitations

- Web searches may miss very recent preprints (last few weeks)
- Some tools may exist only as GitHub repos without a publication
- Non-English publications could be missed
- The search terms above are biased toward the terminology used in our manuscript; tools from adjacent fields (e.g., synthetic biology, protein engineering) may use different vocabulary — Rounds 2, 5, and 6 attempt to address this
