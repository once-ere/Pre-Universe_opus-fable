#!/usr/bin/env bash

set -euo pipefail

repository_root="$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")"
cd "$repository_root"

export SOURCE_DATE_EPOCH=1789603200
export FORCE_SOURCE_DATE=1

mkdir -p build
reports=(
  gpt5_6_cosmology
  gpt5_6_dark_sector_relationships
)

for report in "${reports[@]}"; do
  for _ in 1 2 3; do
    pdflatex -interaction=nonstopmode -halt-on-error \
      -output-directory build "docs/${report}.tex"
  done
  cp "build/${report}.pdf" "docs/${report}.pdf"

  printf 'report: docs/%s.pdf\n' "$report"
  pdfinfo "docs/${report}.pdf" | sed -n '/^Pages:/p;/^File size:/p'
done
