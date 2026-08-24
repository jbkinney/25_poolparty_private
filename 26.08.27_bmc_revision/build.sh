#!/usr/bin/env bash
set -euo pipefail

assembly_dir="$(cd "$(dirname "$0")" && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf -- "$build_dir"' EXIT

cp "$assembly_dir/manuscript_clean.tex" "$build_dir/"
cp "$assembly_dir/baseline_submission.tex" "$build_dir/"
cp "$assembly_dir/response_to_reviewers.tex" "$build_dir/"
cp "$assembly_dir/supplementary.tex" "$build_dir/"
cp "$assembly_dir/references.bib" "$build_dir/"
cp "$assembly_dir/sn-jnl.cls" "$build_dir/"
cp "$assembly_dir/sn-vancouver-num.bst" "$build_dir/"
cp -R "$assembly_dir/assets" "$build_dir/"

(
  cd "$build_dir"
  latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript_clean.tex
  latexdiff baseline_submission.tex manuscript_clean.tex \
    > manuscript_tracked.tex
  latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript_tracked.tex
  latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary.tex
  latexmk -pdf -interaction=nonstopmode -halt-on-error response_to_reviewers.tex
)

if rg -n \
  'There were undefined (references|citations)|Reference .* undefined|Citation .* undefined|Float too large|multiply defined (references|citations)' \
  "$build_dir"/*.log; then
  echo "Build failed: unresolved LaTeX references or citations." >&2
  exit 1
fi

cp "$build_dir/manuscript_tracked.tex" "$assembly_dir/manuscript_tracked.tex"
cp "$build_dir/manuscript_clean.pdf" "$assembly_dir/outputs/manuscript_clean.pdf"
cp "$build_dir/manuscript_tracked.pdf" "$assembly_dir/outputs/manuscript_tracked.pdf"
cp "$build_dir/supplementary.pdf" "$assembly_dir/outputs/Additional_file_1.pdf"
cp "$build_dir/response_to_reviewers.pdf" \
  "$assembly_dir/outputs/response_to_reviewers.pdf"

echo "Built clean and tracked manuscripts, Additional file 1, and response letter in outputs/."
