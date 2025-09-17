#!/usr/bin/env bash
set -euo pipefail

# Simple builder for the JOSS-style paper.
# Requires: pandoc (and xelatex for PDF output).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p build

if ! command -v pandoc >/dev/null 2>&1; then
  echo "pandoc not found. Install pandoc to build the paper." >&2
  echo "macOS: brew install pandoc" >&2
  echo "Linux: sudo apt-get install pandoc" >&2
  exit 1
fi

TARGET="${1:-all}"

build_html() {
  pandoc paper.md \
    --from=markdown+tex_math_dollars+raw_tex+smart \
    --citeproc \
    --metadata=link-citations:true \
    --bibliography=paper.bib \
    --standalone \
    --toc \
    -o build/paper.html
  echo "Built build/paper.html"
}

build_pdf() {
  if command -v xelatex >/dev/null 2>&1; then
    pandoc paper.md \
      --from=markdown+tex_math_dollars+raw_tex+smart \
      --citeproc \
      --metadata=link-citations:true \
      --bibliography=paper.bib \
      --standalone \
      -V geometry:margin=1in \
      --pdf-engine=xelatex \
      -o build/paper.pdf
    echo "Built build/paper.pdf"
  else
    echo "Skipping PDF: xelatex not found (install TeX Live)." >&2
  fi
}

case "$TARGET" in
  html) build_html ;;
  pdf) build_pdf ;;
  all) build_html; build_pdf ;;
  clean) rm -rf build ;;
  *) echo "Usage: $0 [html|pdf|all|clean]" >&2; exit 2 ;;
esac

