# Prompt: Polish Text

When asked to critically review manuscript text (e.g., a LaTeX paper), apply the following levels of scrutiny. Be thorough and do not hold back — the user wants to hear every issue, not just the major ones.

## What to check

### Typos and spelling
- Misspelled words (e.g., "nonlinar" → "nonlinear", "recieved" → "received")
- Hyphenation errors (e.g., "upper-most" → "uppermost")

### Grammar
- Subject-verb agreement (e.g., "the departure ... have been" → "has been")
- Dangling or misattached modifiers (e.g., "When modeled in the weak-mutation limit, the distribution..." — the distribution isn't what's being modeled)
- Broken parallel structure (e.g., "show either X or can be Y" — "either...or" must connect parallel grammatical forms)
- Correlative conjunctions: "both...and" (not "both...as well as"), "either...or", "neither...nor"
- Missing commas after introductory phrases (e.g., "In particular, ...")

### Consistency
- Hyphenation consistency: e.g., "power law" as a noun vs. "power-law" as an adjective — flag if the same form is used inconsistently
- Terminology: flag if the same concept is described with different words in nearby passages without clear reason

### Precision
- Singular vs. plural when the referent is unambiguous (e.g., "highest-fitness sequences" when there is one unique optimum)
- Vague phrasing where the abstract or another section is more specific (e.g., "though with caveats" vs. "though with a reduced range of validity")
- Claims that don't quite match the cited source — when analysis notes are available for cited works, verify each citation is used correctly

### Style and flow
- Colloquial or informal phrasing that doesn't match the register of the rest of the paper (e.g., "methods of doing this" → "methods for doing so"). However, respect intentional informality — if the paper uses playful language deliberately (e.g., "stubby"), don't flatten it.
- Filler words that add nothing ("actually", "really", "very" when not earning its keep)
- Repetition of the same word in close proximity (e.g., "evolution ... evolution" in the same sentence)
- Repetition of the same claim across multiple paragraphs — flag when a point is made too many times and suggest which instances could be cut or reframed

### Structural
- Whether each paragraph earns its place and advances the argument
- Whether the ordering of ideas is logical
- Whether transitions between paragraphs are smooth or abrupt

## How to present findings

- Order findings from most to least significant.
- Group by category (Typos, Grammar, Style, etc.).
- For each issue, quote the problematic text and explain the problem concisely.
- When possible, suggest a specific fix.
- Do not suggest changes that are purely matters of taste unless the current phrasing is genuinely awkward.

## Citation verification

When analysis notes (e.g., `*_analysis.md` files) are available for cited works:
- Read the analysis notes for every work cited in the passage under review.
- Verify that each citation supports the claim it is attached to.
- Flag citations that are correct but weak (e.g., a general review cited for a specific claim that another paper makes more directly).
- Suggest additional citations only if there is a clear gap and a strong candidate exists in the analyzed literature.
