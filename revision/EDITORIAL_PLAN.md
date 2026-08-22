# Editorial plan for the referee response

## Objective

Produce three document products within the revision deadline:

1. a point-by-point response to reviewers;
2. a clean revised manuscript; and
3. a tracked-changes manuscript.

Codex is the editorial orchestrator and primary writer. Existing technical
plans, findings, code, tests, benchmarks, branches, and PRs are evidence inputs.
Other agents are consulted only when those records leave a factual gap. They do
not independently draft or edit the manuscript or central ledger.

## Authoritative sources

- Referee report: `revision/REFEREE_REPORT_26.07.29.md` (verbatim; do not edit).
- Submitted manuscript: `26.05.18_bmc_submission/latex/main.tex`.
- Do not use the older `manuscript/main.tex`.
- Frozen manuscript baseline:
  - last manuscript commit: `da98f615d67aa3e4a76860ad3796f243ba8209c4`
  - SHA-256: `cbdca45f1d67280e14c0256100d22f3a2a90cc531cfceb2bbdabafbe29dfc9e6`
  - source length: 294 lines
- Response style:
  - `revision/templates/response_to_reviewers.tex`
  - `revision/templates/26.08.07_gaugefixer_response_to_reviewers.pdf`
  - `revision/templates/13059_2022_2661_MOESM2_ESM (4).docx`
- Tracked-change example:
  `revision/templates/BIOINF-2026-0541.R1_Proof_hi.pdf`.

If the manuscript hash changes before approved edits are applied, re-anchor all
affected proposals before continuing.

## Working files and ownership

- This file records the workflow and rules.
- `revision/REVISION_PROPOSALS.md` is the only active editorial ledger.
- Codex alone drafts and updates the ledger.
- The author approves, revises, or rejects each proposal.
- Codex later applies only approved changes to the authoritative manuscript and
  constructs the response and tracked version.

No topic-specific `EDITORIAL_PROPOSAL.md` files are required. Existing topic
files remain unchanged and serve as evidence.

## Evidence rule

Every factual proposal records one evidence state:

- `Merged`: present and verified on the final package branch.
- `Branch/PR`: implemented and tested but not yet verified on final `main`.
- `Analysis`: supported by completed analysis, data, or literature review only.
- `Incomplete`: not ready to support a final manuscript claim.

Drafting may begin from any state. Before submission, every claim that a feature
was added or implemented must be verified against final package code and docs.

## Proposal format

Each central-ledger entry contains only:

1. **Comment:** exact reviewer text and stable ID.
2. **Evidence:** verified facts, evidence state, and source paths.
3. **Manuscript edit:** section, baseline source line, exact text anchor, current
   text, and minimal proposed text; or `No manuscript change`.
4. **Draft response:** direct response in the authors' voice.
5. **Decision:** `Draft`, `Needs author decision`, `Approved`, `Applied`, or
   `Verified`.

Baseline source lines locate proposals. Final response letters use rendered line
numbers from the compiled revised manuscript, calculated only after edits settle.

## Supplementary assets

- `REVISION_PROPOSALS.md` contains the single supplementary-asset registry.
- Give every figure and table a stable semantic LaTeX label that describes its
  content, such as `fig:benchmark-scaling`; do not encode `S1` or `S2` in the
  label.
- Treat displayed numbers such as Figure S1 and Table S1 as provisional until
  all approved assets are assembled in their final order.
- During drafting, cite the semantic `\ref{...}` label rather than typing a
  supplementary number into manuscript or response text.
- After the supplement is compiled, update the response letter with the final
  rendered numbers and verify every cross-reference.
- Do not create a separate planning file for supplementary material. The asset
  registry tracks ordering and readiness; the central proposals still contain
  the substantive manuscript and response wording.

## Writing rules

- Match the vocabulary, tense, and plain style of the surrounding manuscript.
- Make the smallest change that fully addresses the comment.
- Treat author-provided drafts as the source text. Change them only to correct a
  factual or grammatical problem or to insert final references; flag any larger
  proposed rewrite for author approval.
- Prefer direct, slightly abrupt, naturally uneven prose over polished
  transitions or symmetrical AI-style phrasing. Do not manufacture typos or
  grammatical errors.
- Preserve LaTeX commands, citations, labels, and terminology unless necessary.
- Narrow unsupported claims and state limitations directly.
- Do not add completed software work unless it materially answers a comment.
- In responses, answer the point first, state what changed, give evidence where
  needed, and quote the exact revised passage.
- When comments substantially overlap, designate one complete primary response
  and answer each repeated comment with a brief acknowledgment and explicit
  cross-reference. Do not duplicate the full response or leave any comment
  unanswered.
- Disagree respectfully where the evidence does not support the premise.
- Avoid promotional or stock AI language.
- Follow the concise, factual voice and point-by-point structure of the
  GaugeFixer and MAVE-NN examples.

## Fast workflow

1. **Survey and draft:** read the manuscript, report, templates, and existing
   evidence; draft directly in `REVISION_PROPOSALS.md` in four batches.
2. **Author decisions:** present each batch for approval while drafting the next.
3. **Apply:** make one controlled manuscript pass containing approved edits,
   assemble approved supplementary assets in their final order, then compile
   and inspect the clean manuscript package.
4. **Finalize:** write the response from applied text, add final rendered line
   numbers and supplementary numbers, generate the tracked version with
   `latexdiff`, compile all outputs, and verify every
   comment/change/quotation/cross-reference.

## Drafting batches

1. Editor and Reviewer 1: benchmarks, comparisons, and design cards.
2. Reviewer 2: statistics, test ownership, mechanism figure, comparisons,
   filters, and naming.
3. Reviewer 3: scope, VCF/VEP, limitations, calibrated claims, and disclosure.
4. References, minor wording, presentation choices, and cross-comment cleanup.

The author may review completed batches immediately; consolidation does not wait
for every topic to be drafted.
