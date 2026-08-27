# Expressiveness experiment — rebuilding other tools' own examples in PoolParty

Answers **R2 2b** (compare pool outputs against existing tools) and **R1 #5** (code
conciseness on a common use case). These cannot be answered by a table; they need
code that runs and output that can be counted.

**Method.** Take a tool's *own published example library* — the design its authors
chose to demonstrate it with — and express it as a PoolParty DAG. Report per
attempt: can express (Y/N) · invocations needed · lines of code · any capability
that blocked it.

Using each tool's own example, rather than one we invent, removes the objection
that the comparison was rigged by choosing a task that suits PoolParty.

## Layout

```
examples/
├── <tool>_<example>/
│   ├── source.md      their spec, cited to paper/repo, with expected counts
│   ├── poolparty.py   our reconstruction
│   └── output/        what we produced, beside what they state
```

## Candidate targets

Inventoried during the tool survey — see each record's "documented examples /
vignettes" heading in `../records/`. Several records label these explicitly as
*candidates for PoolParty reproduction*.

| Tool | Examples inventoried | Strongest target | Why |
|---|---|---|---|
| VaLiAnT | 3 shipped in `examples/`, each with `run*.sh`, `check*.sh` and committed expected output | `brca1_nuc` — BRCA1 exons 2–5 SGE, 4 targetons in one invocation | States unique-variant counts per exon: **583 / 740 / 825 / 1185** |
| MPRAnator | 5 | "MPRA Motif Use Case Example" — AP1/ELK1/RFX over two mm9 backgrounds | States a total: **"This results in a total of 5856 sequences."** |
| Oligopool Calculator | 13 | not yet chosen | |
| tangermeme | 20 | not yet chosen | |
| DNA Chisel | several | not yet chosen | |

## Two constraints on how the comparison can be run

**MPRAnator cannot be executed.** Its web tool returns **HTTP 500 on its own sample
data** — verified during the survey on 3 of its 5 examples, including the Motifs
and SNPs sample inputs. The comparison against MPRAnator must therefore be against
its *documented* expected count (5856), not against a live run. State this in the
response letter rather than letting a referee discover it.

**VaLiAnT validates by md5sum against committed output files.** A PoolParty
reconstruction will not be byte-identical — different tool, different output
format. "Comparison of pool outputs" has to mean library **content**: unique
variant counts and variant identities. Say so explicitly, so the work is not
measured against a bar no reimplementation could clear.

## The mirror direction

The strongest single result needs no installation: **VaLiAnT cannot express the GB1
library.** Its `MutatorType` is a closed 7-member enum with no stochastic member and
no combination mutator, so sampled and higher-order mutagenesis are structurally
unavailable. This follows from the source and can be stated without running
anything — but it belongs here, beside the forward direction, so the comparison is
symmetric rather than a list of things we can do.
