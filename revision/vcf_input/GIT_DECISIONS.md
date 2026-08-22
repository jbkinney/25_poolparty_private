# `from_vcf`: git status and pending decisions

**Written** 2026-08-22 · **Audience** a decision-maker or a fresh agent with no prior context
**Scope** everything needed to decide how to get the `from_vcf` branch onto GitHub. Self-contained;
no earlier conversation required.

---

## 1. What this is about

`from_vcf` is a new PoolParty source operation: it reads a VCF plus a reference FASTA and returns a
pool of reference-genome windows around each variant. It was written to answer part of Reviewer 3's
comment on the BMC Bioinformatics submission — how variants from clinical resources reach
sequence-based models (SpliceAI, AlphaGenome, AlphaMissense, EVE).

The code is finished and verified. **Nothing has been pushed.** What remains is a set of decisions
about *how* to publish it, plus a handful of defects found in review. This document holds the facts
and the open questions.

---

## 2. Repository facts (all independently verified)

| Fact | Value |
|---|---|
| Remote | `jbkinney/poolparty-statetracker` |
| **Visibility** | **PUBLIC** (`isPrivate: false`) — matters for §4.1 |
| Layout | monorepo: `poolparty/` and `statetracker/` under the git root |
| Tooling | git 2.34.1, gh 2.4.0, authenticated as `Zhihan-Leo-Liu`, token `push: true, admin: false` |
| Working branch | `from-vcf` @ `f1044d8`, in worktree `/mnt/c/Users/zhliu/Desktop/KinneyLab/poolparty-from-vcf` |
| Merge base | `4f125f2` (the PR #16 merge commit) |
| `origin/main` | `7b990ae` — **7 commits ahead** of the merge base (PR #17, ORF geometry). Confirmed current by `ls-remote` |
| Branch size | **6 commits**, 9 files, **+1688 / −0**, 70 objects, ~70 KB of text |
| Remote branch | `from-vcf` and `feature/from-vcf` **do not exist** on origin — the push creates a new ref, no force needed |
| Local `main` | `1bb0179`, **19 commits behind** origin/main, strict ancestor of the branch. Checked out in the `poolparty-statecounter` worktree, which is **dirty** (`M poolparty/README.md`) |
| Worktrees | **6**, sharing one object + ref store (common dir in `poolparty-statecounter/.git`) |
| Shared stash | `stash@{0}` holds **another session's** WIP on `orf-geometry-consistency` |
| Tags | `v0.1.0` (`f6ce1b0`), `v0.1.1` (`e98fb64`) are ancestors of HEAD and exist on the remote at the same SHAs |
| Branch protection | **NONE.** `branches/main/protection` → 404, rulesets → `[]`. No required checks, no required reviews |

### The 9 files

New (created by this branch):
- `poolparty/src/poolparty/base_ops/from_vcf.py` (+538)
- `poolparty/tests/test_from_vcf.py` (+843)
- `poolparty/docs/operations/from_vcf.rst` (+275)

Pre-existing, insertions only:
- `poolparty/CHANGELOG.md` (+14), `poolparty/docs/api.rst` (+2),
  `poolparty/docs/metadata/design_cards.rst` (+6),
  `poolparty/docs/operations/source_operations.rst` (+3),
  `poolparty/src/poolparty/__init__.py` (+4), `poolparty/src/poolparty/base_ops/__init__.py` (+3)

**Precision note:** "−0" is true of the diff against `origin/main`. Per-commit, `8d79b5c` deletes 7
lines from `CHANGELOG.md` — but those 7 were added by `3c48cd8` earlier on the same branch. **No
pre-existing content is removed at any point.** Both statements matter; use the right one.

### Clean (checked, nothing to worry about)

No binaries, no `__pycache__`/`.fai`/`.pyc`/editor droppings, no secrets or absolute paths in
content, no CRLF in tracked blobs, no submodules, no LFS, no `.gitattributes`, everything under
`poolparty/` (nothing touches `statetracker/` or repo-root files), no in-progress rebase/merge in
any of the 6 worktrees, branch checked out in exactly one worktree, working tree pristine.

---

## 3. Verification status

| Check | Result |
|---|---|
| Full suite on the **merged** state (branch + origin/main) | **3072 passed, 14 xfailed** — reproduced twice, independently |
| Full suite, branch alone | 3037 passed, 14 xfailed (Δ35 = main's new ORF tests) |
| `test_from_vcf.py` alone | 63 passed |
| Merge with `origin/main` | **Conflict-free.** `CHANGELOG.md` is the only file both sides touch; it auto-merges |
| `ruff check` + `ruff format --check`, new files | Clean at the **pinned** v0.9.1 (local ruff is exactly 0.9.1) |
| pre-commit hooks (all 8) | New files pass every one |
| Sphinx build | Succeeds. `from_vcf` is in the toctree, `api.html`, `source_operations.html`. 11 warnings, **none** attributable to the new files |
| Dependencies | `pyfaidx>=0.8.1` already declared; **no new dependency**. `gzip` is stdlib |
| Python floor | `requires-python = ">=3.10"`; CI matrix is 3.10–3.14 × ubuntu/macos/windows |
| Windows portability | Simulated Windows CRLF fixture writes: still 3072 passed. `pyfaidx` handle is closed in a `finally` |

### How CI actually behaves (settled with three independent proofs)

- `test.yml` triggers on `push: [main, master, develop]`, `pull_request: [main, master, develop]`,
  `workflow_dispatch`.
- **Pushing the feature branch runs no CI at all.** Opening the PR is what triggers it.
  (Empirical: every non-`main` branch in this repo's history has exactly one run, event
  `pull_request`.)
- **`pull_request` CI tests the *merged* state.** `actions/checkout` with no `ref:` checks out
  `GITHUB_SHA` at `GITHUB_REF` = `refs/pull/N/merge`. Proven by this repo's own runner log for
  PR #17: `HEAD is now at 518876c Merge 975c2ef... into 4f125f2...`.
  **Therefore merging `origin/main` into the branch first is unnecessary.**
- The `lint` job is `continue-on-error: true`, so style cannot fail the PR.
- **No Sphinx build in CI, and Read the Docs posts no PR check.** The 275 new RST lines are gated by
  nothing automated.
- **Because there is no branch protection, CI is advisory** — a red PR is mergeable — and CI does
  **not** recompute when `main` advances, so a green check can go stale.

---

## 4. Decisions required

### 4.1 — BLOCKER: what may the commit messages say in public?

The repo is public. Commit `31e4a78`'s message body contains:

```
line 4:  ...Answers the ingress half of Reviewer 3's comment: variants from clinical
line 8:  Design is recorded in 25_poolparty_private/revision/vcf_input/DESIGN.md; this
```

`25_poolparty_private` is tracked nowhere in the repo — it exists only as a string in that message.
Pushing therefore publishes the path of a private design document **and** the fact that a manuscript
is under peer review together with what a specific anonymous referee asked. BMC Bioinformatics peer
review is confidential until publication.

Secondary, same theme: `3c48cd8` and `6e821e5` say "three independent reviews" and "all three
reviewers converging" — these were AI review agents, not journal referees, but a public reader
cannot tell, and next to `31e4a78`'s "Reviewer 3" the ambiguity reads the wrong way. All six commits
also quote internal suite counts (`Full suite 3037 passed, 14 xfailed`), which nothing upstream does.

**This is fixable only before the first push.** A later force-push does not remove a message from
GitHub's events API or from anything that already fetched it.

| Option | Effect |
|---|---|
| **A (recommended)** — strip the private path and the referee reference; describe the feature on its own technical terms | Nothing confidential published. Loses the "why now" from the permanent history, which belongs in the PR body anyway |
| B — strip the path, keep a generic "requested by a reviewer" | Path protected, but still discloses that the work is review-driven |
| C — push as-is | Publishes a private tree's layout and referee content on a public repo, permanently |

### 4.2 — BLOCKER: which rewrite tool?

`git filter-branch` is **ruled out**, not a choice: un-ranged it rewrites **429** commits rather
than 6 (measured), destroying the merge base so the PR would render the whole repository as added,
with `poolparty: poolparty:` on the 247 commits that already carry the prefix. It also writes to the
shared `refs/original/`, creates an un-gitignored `.git-rewrite` inside the working tree, and its
`-f` retry destroys its own backup.

The two live options both keep the branch ref pristine until success:

| | `git rebase 4f125f2 --exec` | `git commit-tree` replay |
|---|---|---|
| Commits touched | exactly 6 | exactly 6 |
| Working-tree writes | 6 checkouts, on DrvFs | **none** |
| Result lands on | the branch, at the end | a **new** ref; original untouched |
| Recovery | `git rebase --abort`; then `reset --hard ORIG_HEAD` | delete the new ref |
| Validation done | **byte-identical round-trip of all 6 real messages**, `Co-Authored-By` trailers preserved, plus a negative test for a concurrent commit | tree-hash equality gate designed, not yet executed |

Recommendation: **`rebase --exec`**, on the strength of the empirical validation. The counter-argument
is real — this filesystem (WSL2 DrvFs over `/mnt/c`) has a documented corruption incident in this
project, and `commit-tree` touches no working tree at all. Either is defensible.

### 4.3 — What should the reworded messages actually say?

Measured, so the two halves can be decided separately:

- **The `poolparty:` prefix is a real convention.** `CONTRIBUTING.md` mandates it; of the last 20
  non-merge commits on `main`, 14 carry it, including all 14 from PRs #15/#16/#17. **All 6 commits
  on this branch lack it.**
- **It is contested only on unmerged branches.** `feature/orf-indel-scans` (0/6) and
  `feature/pool-stats` (0/7) use Conventional Commits (`feat(orf):`, `fix:`). None of that is
  merged, so it is not evidence about what the maintainer accepts.
- **Subject length is *not* a house rule.** 68 of 247 `poolparty:` subjects on `main` exceed 72
  characters; the longest is 592. Shortening is optional polish.
- **Bodies are the real outlier.** Upstream bodies run 3–12 lines. These run 17–69 and narrate the
  review process ("Three fix rounds each introduced regressions in the same place"). Three of the
  six subjects also describe the review rather than the change.

Candidate subjects (≤72, prefixed, lowercased to match the most recent 14):

```
poolparty: add from_vcf, a source operation reading variants from a VCF
poolparty: fix from_vcf defects found by three independent reviews
poolparty: fix from_vcf ordering bugs and three untrue doc claims
poolparty: rebuild the from_vcf guard sequence from an explicit table
poolparty: rewrite the from_vcf docstrings in the package's voice
poolparty: align from_vcf docs with the house page conventions
```

Note subject 2 still says "three independent reviews" — decide it together with §4.1.

### 4.4 — The ~32 hour clock skew

The WSL clock is genuinely **~32 h behind** real time (local `2026-08-22T10:32Z` vs origin
`pushed_at: 2026-08-23T18:25Z`). Windows is correct; Linux drifted. A reword resets committer dates
to "now", i.e. **~32 h earlier than the base commit they sit on**, and GitHub orders PR commit lists
by date.

Options: resync (`sudo hwclock -s`) before rewriting · accept the backdating · preserve the original
committer dates explicitly during the rewrite.

### 4.5 — Branch name

Every human branch on the remote and in the other worktrees uses a `type/` prefix
(`feature/orf-indel-scans`, `fix/orf-geometry-consistency`, `docs/constraint-filters`), and
`CONTRIBUTING.md` says `feature/your-feature`. `from-vcf` is the only bare name.

Recommendation: rename to **`feature/from-vcf`**. Free right now (exists nowhere but locally); after
a push it becomes a delete-and-recreate on the remote. **Must happen after §5 step 1** — see below.

### 4.6 — Fork or same-repo branch?

`CONTRIBUTING.md` says to fork, but observed practice contradicts it: all 12 PR merges read
`Merge pull request #N from jbkinney/<branch>`, i.e. same-repo branches; there are zero fork PRs.
Three `claude/*` branches already live on origin. Decisive technical point: `test.yml` passes
`CODECOV_TOKEN` from secrets, which GitHub withholds from fork PRs.

Recommendation: **same-repo branch.**

### 4.7 — Which review findings to fix before the PR?

| # | Finding | Severity | Note |
|---|---|---|---|
| a | `CHANGELOG.md` now has **two** `### Added` headings under `[Unreleased]` (the merge base had one at line 65; this branch added a second at line 10) | Objective defect | Introduced by this branch |
| b | `docs/operations/library_size.rst:139` lists all seven source ops **exhaustively** and `from_vcf` is missing | Real gap | Missed in the docs audit |
| c | The reference-mismatch `ValueError` asserts a conclusion the data cannot support. At one record: `1 of 1 records ... disagree ... The two files are most likely different genome builds or assemblies.` Also `_MISMATCH_LIMIT = 0.2` is hard-coded while its neighbour `max_allele_length` is a parameter | Blocking, per review | Reproduced |
| d | A BOM-prefixed VCF dies with `line 1: expected at least 8 tab-separated fields`. `encoding="utf-8"` should be `utf-8-sig`. Windows tooling writes BOMs and Windows is in the CI matrix | Real bug | Reproduced |
| e | 3 of 8 test classes are named after the review (`TestReviewFindings`, `TestGuardsAtTheirDefaults`, `TestGuardTable`); siblings use `Test<Op><Aspect>`. Several duplicate tests | Reviewability | Tests are 843 of 1688 lines — the first thing a maintainer reads |
| f | `warnings.simplefilter("error")` in two tests turns any future `pyfaidx`/`pandas` deprecation into a failure; floors are unpinned | Latent | Not a today-risk |
| g | `README.md:117` "Create pools" row omits `from_vcf` | **Do not touch** | Another session has `README.md` uncommitted in the `main` worktree |

### 4.8 — PR body disclosures

- **`INFO` fields with `Number=A` are not split per allele** (a multi-allelic `AF` gives every
  alternate the whole comma-separated value). Correct splitting needs the header's `Number=`.
  Review verdict: **disclosure is sufficient** — nothing fails silently, and the remedy
  (`bcftools norm -m-`) is documented.
- **The 20% reference-mismatch threshold.** Review verdict: **disclosure is not sufficient** — see
  4.7c. Fix the message wording at minimum.

---

## 5. Mandatory config fix — not a decision, and it must go first

```
branch.from-vcf.remote = origin
branch.from-vcf.merge  = refs/heads/main      <-- wrong
push.default           = unset in all 4 scopes -> "simple"
```

Two independent consequences:

1. **`gh` already resolves this branch to the wrong head.** `gh pr status` prints
   `There is no pull request associated with [main]`. `gh pr create` would target head `main`.
2. **A bare `git push` is refused, and git's own first suggested remedy is
   `git push origin HEAD:main`** — which would land 6 unreviewed commits on an unprotected default
   branch with a push-capable token.

`git branch -m` **carries `branch.<name>.merge` through the rename**, so renaming first bakes the
problem into the new name. Unset it before anything else:

```bash
git config --unset branch.from-vcf.merge
git config --unset branch.from-vcf.remote
```

(The same misconfiguration exists on `feature/pool-stats` and `docs/constraint-filters` — not in
scope here.)

---

## 6. Execution order, once the decisions are made

```
0. Preflight, read-only. Abort on any mismatch:
     HEAD == f1044d8 · origin/main == 7b990ae · merge-base == 4f125f2
     rev-list --count origin/main..HEAD == 6 · status --porcelain empty
     for-each-ref refs/original empty · no .git-rewrite · stash list == 1 (DO NOT TOUCH)
     no *.lock under the common dir

1. Unset branch.from-vcf.{merge,remote}          <-- FIRST, see §5
     verify: gh pr status no longer says [main]

2. Reword the 6 messages (tool per §4.2, content per §4.1 + §4.3)

3. GATE. Nothing leaves the machine until all pass:
     git diff --stat ORIG_HEAD HEAD                    -> MUST BE EMPTY  (proves messages only)
     git rev-list --count origin/main..HEAD            -> 6
     git merge-base origin/main HEAD                   -> 4f125f2
     git log --format='%s' origin/main..HEAD | grep -c '^poolparty: '  -> 6  (not 253, not 429)
     git log --format='%B' origin/main..HEAD | grep -c 'Co-Authored-By' -> 6
     git for-each-ref refs/original                    -> empty
     git show-ref --tags                               -> v0.1.0/v0.1.1 unchanged

4. git branch -m from-vcf feature/from-vcf
     verify worktree HEAD and `git worktree list`

5. git push --dry-run --verbose origin \
        refs/heads/feature/from-vcf:refs/heads/feature/from-vcf
     expect exactly one "* [new branch]" line, then push for real.
     Never -u, never --tags, never --follow-tags, never HEAD:main.

6. gh pr create --base main --head feature/from-vcf --draft --title ... --body-file ...

7. Verify the PR shows 6 commits. If it says 429, the rewrite went wrong:
     close the PR and delete the branch immediately.
```

### Irreversibility

| Step | Reversibility | Recovery |
|---|---|---|
| 1 config unset | trivial | re-set the two keys |
| 4 rename | trivial | `git branch -m feature/from-vcf from-vcf` |
| 2 reword | easy | `git rebase --abort` mid-flight; `git reset --hard ORIG_HEAD` after |
| 5 push | semi | `git push origin --delete feature/from-vcf` — but SHAs are public and CI runs are permanent |
| 6 PR | **point of no return** | closing does not undo it: `refs/pull/N/head` persists on GitHub permanently, and the notification to the maintainer cannot be unsent |

### Standing prohibitions for the whole run

- **No `git stash`** — `refs/stash` is one shared ref and holds another session's orf-geometry WIP.
- **No `git add -A`** — has already staged another session's work once in this project.
- **No `git gc`, `reflog expire`, or `worktree prune`** — after a rewrite the reflog is the only backup.
- **Never `HEAD:main`, never a bare `git push`, never `-f`.**
- **Never `--tag-name-filter`** — `v0.1.0`/`v0.1.1` are ancestors of HEAD and exist on the remote.
- **Do not touch** local `main`, the merged-but-undeleted `fix/*` remote branches, or
  `poolparty/README.md` (dirty in another worktree).
- **Do not merge `origin/main` into the branch** — verified unnecessary (§3).

---

## 7. After the PR

- **Re-check that `main` has not advanced before merging.** There is no branch protection and no
  up-to-date requirement, and CI does not recompute on base-only movement, so a green check can be
  stale. If `main` moved, re-trigger before merging.
- `CHANGELOG.md` `[Unreleased]` is the one file both sides edit; expect that to be where any future
  conflict lands.
- No version bump is owed — both packages sit at 0.1.1 and features accumulate under `[Unreleased]`.
- Consider a separate PR adding a Sphinx build to CI; docs are currently ungated.

## 8. Explicitly out of scope

Local `main` being 19 commits behind · the merged-but-undeleted `fix/*` branches on origin · the
stale `poolparty-orf-geometry-fixes` worktree (its commit is already in `origin/main`) · the two
stale dependabot PRs (#13, #14) · the same upstream misconfiguration on `feature/pool-stats` and
`docs/constraint-filters` · the repo-wide pre-existing ruff baseline (108 errors, 65 files would
reformat; CI's lint job cannot fail).

---

## 9. Decision checklist

| # | Decision | Recommendation | Blocking? |
|---|---|---|---|
| 1 | What commit messages may say publicly (§4.1) | A — strip the private path and referee reference | **Yes** |
| 2 | Rewrite tool (§4.2) | `rebase --exec` | **Yes** |
| 3 | Reword scope (§4.3) | Add the prefix (required); rewrite the review-narrating subjects and cut bodies; shortening is optional | **Yes** |
| 4 | Clock skew (§4.4) | Resync before rewriting — cheap | **Yes** |
| 5 | Branch name (§4.5) | `feature/from-vcf` | **Yes** |
| 6 | Fork vs same-repo (§4.6) | Same-repo | **Yes** |
| 7 | Which findings to fix first (§4.7) | a, b, c, d before the PR; e if the maintainer's first read matters; f later; **skip g** | **Yes** |
| 8 | PR disclosures (§4.8) | Disclose `Number=A`; **fix** the mismatch message | **Yes** |
