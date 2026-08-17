# Adversarial review — Mutation Maker extraction

**Reviewed:** `extractions/mutationmaker.md` + the structured capability block
**Reviewer method:** independent re-read of `papers/Hiraga2021_MutationMaker.pdf` (PyMuPDF), independent
GitHub API / raw-file inspection of `ra100/Mutation_Maker@master`, GitHub code search, PyPI checks.
**Date:** 2026-08-10

---

## Headline

The value column is **almost entirely correct**. I could not falsify a single yes/no/partial
*value*. But there are **three things a Mutation Maker author would attack**, and one of them is
serious enough that it must be fixed before the referee response goes out.

### 1. SERIOUS — "community copy" is wrong. `ra100` is a paper co-author.

The extraction says, repeatedly, that the surviving repository is a *"community-maintained copy"* /
*"community copy"* / *"surviving copy"* and frames `maintained: yes` as charity toward a third-party fork.

```
$ gh api users/ra100 --jq '{login,name,company,blog}'
{"login":"ra100","name":"Rastislav Švarba","company":"MSD IT","blog":"https://rast.io"}
```

**Rastislav Švarba** is the 8th listed author of Hiraga et al. 2021 ("Rastislav Svarba"), and his
GitHub company field is **MSD IT** — Merck Sharp & Dohme, the copyright holder. The other
contributors to the repo are `spalemartin` (Martin Spale, co-author) and `prihoda` (David Prihoda,
co-author). The fork network also contains `prihoda/Mutation_Maker` and `swiewiora/Mutation_Maker`
(Sebastian Wiewiora, co-author).

So `ra100/Mutation_Maker` is **an author-and-MSD-employee-maintained repository**, not a community
salvage. Calling it a community fork in a referee response, when the referee may be Bitton or Švarba
himself, is the single most damaging line in the document. Rewrite to something like: *"the repository
URL printed in the paper (github.com/Merck/Mutation_Maker) now 404s; the code is maintained by
co-author Rastislav Švarba (MSD IT) at github.com/ra100/Mutation_Maker."*

Maintenance is also stronger than presented — it is not a single 2026 burst:

| Month | Commits |
|---|---|
| 2020-07 | 9 |
| 2020-11 | 3 |
| 2021-01 | 4 |
| 2021-02 | 10 |
| 2022-08/09/10 | 20 |
| 2025-07 | 7 |
| 2026-02 | 32 |

CI + Primer3 workflows both `success` on 2026-02-14; not archived; 0 open issues; 19 stars; 13 forks;
`/tags` and `/releases` both empty (that part of the extraction is correct).

### 2. SERIOUS — the absolute "never emits the enumerated library members" is falsifiable.

This sentence appears in **four** places (`library_as_object`, `per_sequence_provenance`, `assay_dms`,
`stated_limitations`). It is **wrong for the PAS workflow**.

`backend/mutation_maker/generate_oligos.py`:

```python
def generate_oligos_from_combinations(mutations_combinations_with_probabilitites, chosen_codons_on_sites,
                                      dna, mutations_sites, start, goi_offset) -> List[PASOligo]:
    """ From combination of mutations generates a dna sequence with respective codons on mutation sites
    for all combinations. """
    oligos = []
    for combination in mutations_combinations_with_probabilitites:
        concentration = mutations_combinations_with_probabilitites[combination]
        for position, mutation in enumerate(combination):
            codon = chosen_codons_on_sites[position][mutation]
            dna = Codons.replace_codon(dna, mutations_sites[position], codon, start, goi_offset)
        oligos.append(PASOligo(sequence=dna, ratio=concentration))
```

with `cartesian_dictionary(...)` = `itertools.product` over the per-site codon dictionaries, folding
the per-site frequencies by multiplication, then normalising the ratios to sum to 1. That is
**explicit combinatorial enumeration of variant sequences with per-variant frequency** — it is just
scoped to a fragment rather than a full-length gene. Each enumerated sequence is emitted as
`PASOligoOutput{sequence, mix_ratio, mutations, reds, blues}`.

Two more emissions the memo denies:
- **PAS returns the reverse-translated, codon-optimised full-length gene.**
  `pas.py` calls `sequences.translate_goi_sequences(translator)`, which mutates
  `PASSequences.gene_of_interest` **in place**, and then returns
  `PASOutput(input_data=workflow_input, results=results)`. So the designed DNA sequence for the whole
  gene comes back in the payload.
- `PASResult.fragment` is itself a designed sequence with `start/end/length/overlap/overlap_Tm/overlap_GC`.

The correct, defensible framing is: *"Mutation Maker enumerates variant sequences only at the oligo /
fragment level (PAS), never as full-length assembled library members, and never for SSM or QCLM — those
two emit primers only."* That is still a real gap vs. PoolParty and it cannot be rebutted.

### 3. MODERATE — the "raw sequence strings only" evidence for Block C has a hole.

`genome_coordinates` and `transcript_models` are both correctly **no**, but the stated reason
("Inputs are raw sequence strings only") invites the rebuttal *"Mutation Maker reads GenBank."*
`frontend/src/shared/components/FileUploadInput/index.tsx`:

```ts
const FILE_EXTENSIONS = '.fa, .fasta, .gb';
...
if (text.indexOf('ORIGIN') > -1) {
  return text.split('ORIGIN')[1].split('')
    .filter(letter => /[ATCGatcg...]/.test(letter)).join('').toUpperCase()
}
```

It accepts `.fa/.fasta/.gb` — but the "parser" splits on `ORIGIN` and keeps only the letters, i.e. it
**discards every FEATURES key, every coordinate, every annotation**. Say so explicitly; that turns a
potential rebuttal into a stronger point.

Likewise PAS mutations can be uploaded as an **XLSX table with columns `Site` / `MT` / `MT%`**
(`PASForm/components/InputMutations/FileUploadMutations.tsx`). That is a tabular library specification
document, which is relevant to `library_as_object` and is not mentioned anywhere in the memo.

---

## Verified-correct claims (spot checks I ran independently)

- `Merck/Mutation_Maker` → `gh api` returns `{"message":"Not Found","status":"404"}`. **Confirmed.**
- `SSMInput.degenerate_codon = StringProperty(required=True, default="NNS")` — **confirmed**, one
  job-wide field. The extraction's open question ("does the UI support per-site degenerate codons?")
  is now **settled: no.** `frontend/src/scenes/SSM/components/SSMForm.tsx` has a single
  `degenerateCodon` form field plus an "Amino acids" tab that calls `getDegenerateCodons(include, avoid)`
  and writes the result back into that **one** field. Note the paper says default **"NNK"** (p. 359)
  while the code defaults to **"NNS"** — a paper/code discrepancy worth a footnote.
  Also: the field accepts a **comma-separated list** of codons (`degenerate_codon.split(",")[0]` in
  `ssm.py`; `joinGrammars` in `amino.ts` emits `"NDT,VHG"`-style multi-codon grammars), so it is a
  job-wide codon *set*, not strictly a single codon.
- `api/api.py` schedules `tasks.export_ssm` / `tasks.export_qclm` while `backend/tasks.py` defines only
  `tasks.ssm`, `tasks.qclm`, `tasks.pas`, `tasks.species_table` — **confirmed, the export endpoints are
  dead.** *Correction:* there is **no `/v1/export_pas` endpoint at all**; the three that exist are
  `export_qclm`, `export_ssm`, `export_species_table`. Fix the `interface` evidence string.
- No CLI: `argparse` code-search returns 2 hits, both in `LICENSES_THIRD_PARTY` and
  `frontend/package-lock.json`. **Confirmed.**
- No PyPI, and stronger than stated: there is **no `setup.py` and no `pyproject.toml` anywhere in the
  repository tree** (367 paths). So it is not pip-installable even from a git URL; `mutation_maker` is
  importable only as a source directory on `PYTHONPATH`. Tighten the `interface` wording.
- Code search over `ra100/Mutation_Maker`: `hgvs` 0, `vcf` 0, `gtf` 0, `exon` 0, `intron` 0,
  `chromosome` 0, `barcode` 1 (bootstrap.css). **All confirmed.**
- 35,799 species / 33 genetic-code tables: paper p. 366 confirms; `codon_usage_data/species.table`
  has 35,793 lines. **Confirmed.**
- GPL-3.0 with per-file MSD headers: **confirmed** on every file I fetched.
- Supporting Information is *only* "Figures S1−S4 illustrate Mutation Maker algorithms (PDF)" —
  there is **no user guide**, which supports the "no documentation site" finding.

---

## Things the extraction missed

1. **PAS combinatorial variant enumeration with normalised mixing probabilities** (see §2 above) —
   this is the closest thing in Mutation Maker to a PoolParty-style enumerated pool and it should be
   named explicitly rather than denied.
2. **Reverse-translated full-length gene returned in the output payload** (`translate_goi_sequences`
   mutates `PASSequences.gene_of_interest` in place before `PASOutput` is built). Mutation Maker
   therefore *does* emit a designed full-length DNA sequence — just a single codon-optimised parent,
   not a library.
3. **Automatic parental/WT inclusion at every QCLM site.** `mutation.py`:
   `self.new_aminos = frozenset(old_aminos.union({m.new_amino for m in mutations}))` — the wild-type
   amino acid is unioned into every site's target set *by construction*. WT retention is not an opt-in
   trick you get by typing `P269P`; it is unconditional. This strengthens `mixed_mutagenesis_one_pool`.
4. **Mutations may be specified as raw DNA codons, not only amino acids.**
   `PASInput.is_mutations_as_codons: BooleanProperty(required=True)`; `PASMutationSite` docstring:
   *"Each mutation is a amino acid IUPAC code, or a DNA codon."* This is nucleotide-level control over
   the exact codon placed — a capability axis the row list does not have.
5. **Tabular / file-based library specification input**: PAS mutations from `.xlsx` (`Site`, `MT`, `MT%`),
   sequences from `.fa/.fasta/.gb`.
6. **The SSSM degeneracy solver lives in the browser, in TypeScript**
   (`frontend/src/scenes/SSM/components/amino.ts` — `modifyGrammar`, `joinGrammars`,
   `generateCombinationsFromDegenerateCodon`, `MAX_NUMBER_OF_COMBINATIONS = 100`,
   `MAX_DURATION = 10 min`), *not* in the Python backend. Anyone using the REST API directly does not
   get it. Worth a footnote if you claim the API is a usable programmatic entry point.
7. **No indel support of any kind.** `AminoMutation.__init__` hardcodes `self.length = 3` and every
   mutation is a codon-for-codon substitution; code search `insertion` → 1 hit (irrelevant),
   `deletion` → 0. The memo never states this, and "no insertions or deletions" is a clean,
   unattackable limitation to put in the matrix.
8. **BDD/Gherkin executable specs** — `backend/features/*.feature` (`degenerate_codon.feature`,
   `mutation.feature`, `primer.feature`, `temperature_calculator.feature`,
   `possible_mutagenic_primer_position.feature`, `extract_sequence_from_plasmid.feature`) with
   `features/steps/*.py`. A second documented-example source alongside `tests/test_support.py`.
9. **`api/server_fastapi.py`** — an alternative ASGI server present in the current tree, i.e. the
   `hug`/`falcon`-only description is now slightly dated.
10. **Paper/code discrepancy on the default degenerate codon**: paper p. 359 says default `"NNK"`,
    `ssm_types.py` says `default="NNS"`.

---

## Per-claim verdicts

| key | claimed | verdict | note |
|---|---|---|---|
| library_as_object | partial | supported (evidence must be fixed) | value fine; the "never emits library members" justification is false for PAS |
| dag_chaining | no | supported | three independent endpoints, no composition primitive; verified |
| lazy_evaluation | no | supported | eager; no generator/iterator API; Redis-materialised JSON |
| mixed_mutagenesis_one_pool | partial | supported | strengthened: WT union is automatic at every QCLM site; PAS accepts codon-level specs |
| combinatorial_motif_place | no | supported | motifs = restriction sites to avoid; no insertion, and no indels at all |
| barcode_generation | no | supported | verified |
| per_sequence_provenance | partial | understated evidence | every emitted PAS oligo carries a full `PASMutationFormatted` record; "no records exist" is false |
| automatic_naming | yes | supported | verified in `SaveFile/index.tsx`; SSM header is `['Well Position','Name','Sequence','Notes']` |
| design_visualization | yes | supported | verified in paper + repo |
| assay_dms | yes | supported | correct call; a "no" here would be indefensible |
| assay_mpra | no | supported | note PAS can synthesise arbitrary DNA, but positions are still codon-indexed |
| assay_insilico | no | supported | verified |
| genome_coordinates | no | supported (evidence gap) | must pre-empt "it reads .gb" — the GenBank reader discards all annotation |
| transcript_models | no | supported (evidence gap) | same GenBank caveat |
| exon_intron_split_codons | no | supported | `self.length = 3` hardcoded |
| hgvs_input | no | supported | verified |
| vcf_vep_output | no | supported | verified |
| consequence_annotation | no | supported | `PASMutationFormatted.wild_type` is display-only, not a consequence call |
| primer_design | yes | supported | core of the tool |
| codon_optimization | yes | supported | plus: the optimised gene is actually returned in the output |
| synthesis_constraints | yes | supported | verified in `PASConfig` |
| interface | yes | supported (two factual errors) | no `/v1/export_pas`; no setup.py/pyproject anywhere |
| license | yes | supported | GPL-3.0 verified |
| maintained | yes | **understated** | maintainer is co-author Rastislav Švarba (MSD IT), not "the community" |

## Bottom line

Zero value changes required. Three text fixes are mandatory before this goes anywhere near a referee:
(a) drop "community copy", (b) delete the absolute "never emits the enumerated library members",
(c) pre-empt the FASTA/GenBank upload. Everything else is solid and the Block C "no"s are bulletproof.
