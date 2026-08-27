# Oligopool Calculator — mixed variant types in one library (row 6)

- **Value:** `no`
- **Confidence:** High (0.99)
- **Version inspected:** Oligopool Calculator `2026.02.22.1`, canonical repository commit `b88fa394ca67ed4c48ec41127b5707694ee7cf0a` (2026-02-22)

## Binding operational test

Row 6 requires two or more structurally different component types to be declared in one specification and emitted as one pooled output. Separate runs that a user merges are `no`; under the global rule, user-side row appending or bookkeeping cannot establish the capability.

## Primary-source evidence

The first-party user guide says Design Mode begins with already supplied biological content:

> `1. Start with your core sequences (variants, promoters, genes, etc.)`  
> `2. Add functional elements: primer → motif → barcode → spacer`

Source: [`docs/docs.md`, lines 303–315, commit `b88fa394...`](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/docs.md#L303-L315). Its basic workflow is still more explicit:

> `Start with your variants in a CSV`  
> `df = pd.read_csv('variants.csv')`

Source: [`docs/docs.md`, lines 1300–1313, same commit](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/docs.md#L1300-L1313).

The documented recipe for extending a pool assigns the actual row pooling to the caller:

> `Extend pool with new variants | append rows → barcode(patch_mode=True) → spacer(patch_mode=True)`

Source: [`docs/agent-skills.md`, lines 347–358, same commit](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/docs/agent-skills.md#L347-L358).

The potentially relevant `join` operation does not pool rows or variant classes. It joins columns for the same members:

> `Join two oligo tables on ID and insert new columns from other_data into the input_data column order.`  
> `input_data and other_data must contain exactly the same ID set`

Source: [`oligopool/join.py`, lines 13–45, same commit](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/join.py#L13-L45).

Likewise, `compress` accepts concrete variants with no structural type declaration:

> `Compress concrete DNA sequences into IUPAC-degenerate oligos.`  
> `input_data must contain a unique 'ID' column; all other columns must be non-empty strict-ATGC DNA strings.`  
> `Sequences of different lengths are compressed independently by length group.`

Source: [`oligopool/compress.py`, lines 15–56, same commit](https://github.com/ayaanhossain/oligopool/blob/b88fa394ca67ed4c48ec41127b5707694ee7cf0a/oligopool/compress.py#L15-L56).

## Reasoning

Oligopool Calculator builds synthesis architecture around a caller-supplied row set. It can accept any strict-ATGC rows that the caller has already mixed and can add primers, motifs, barcodes, and spacers to those rows, but it neither declares nor generates structurally distinct variant classes in one specification. `join` is a horizontal same-ID architecture join, not a union of component rows; `merge` similarly concatenates columns within each row. `compress` emits one synthesis DataFrame but merely compresses a pre-existing concrete sequence set and has no operation/parameter representing variant types. Therefore those operations do not satisfy the row's declaration-and-pooling requirement.

The value is `no`, not `partial`. The partial cases still require the tool to emit a pool containing at least two component types. Here the entire variant row set must be supplied pre-built or appended by the user; the tool does not pool even one generated variant type with one pre-built type.

## Disconfirmation attempt

I searched the complete pinned repository, including README, user guide, API reference, AI-agent guide, examples, all package modules, CLI, and pipeline configuration code, for `mixed`, `variant type`, `component type`, `pool`, `combine`, `union`, `append`, `join`, `merge`, `stack`, `insertion`, `deletion`, and related terms. I inspected all design/generation entry points, with particular attention to `join`, `merge`, `compress`, `expand`, motif modes, `patch_mode`, parallel DAG fan-out/recombination, and `final`. No operation declares multiple structural variant classes or unions their members into one row set. `join` requires the same IDs, `merge` concatenates columns, DAG recombination uses `join`, and the documented extension recipe instructs the caller to append rows.

Evidence that would change the score would be a documented one-request mode that declares at least two structurally different library components and returns their pooled rows, or a documented tool-provided row-union operation intended to combine distinct generated component outputs. A restricted two-class mode, or a mode that generated one class while accepting another pre-built class, would support `partial`; broader native mixing would support `yes`. None was found. Merely accepting a CSV in which the user has already mixed substitutions, indels, controls, or other sequences does not change the score because the structural types are not declared to or pooled by Oligopool Calculator.
