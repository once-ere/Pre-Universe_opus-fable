#!/usr/bin/env bash

set -euo pipefail

repository_root="$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")"
cd "$repository_root"

export SOURCE_DATE_EPOCH=1789603200
export FORCE_SOURCE_DATE=1

mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory build docs/gpt5_6_cosmology.tex
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory build docs/gpt5_6_cosmology.tex
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory build docs/gpt5_6_cosmology.tex
cp build/gpt5_6_cosmology.pdf docs/gpt5_6_cosmology.pdf

printf 'report: docs/gpt5_6_cosmology.pdf\n'
pdfinfo docs/gpt5_6_cosmology.pdf | sed -n '/^Pages:/p;/^File size:/p'
