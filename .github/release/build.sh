#!/usr/bin/env bash
# Build the release assets of "The Mirror Clause of the Geometric Selection Conjecture Is Empty":
# the PDF and the LaTeX source (tex.zip).
#
# Usage: build.sh <version> <output directory>
# Needs pdflatex with amsart, xurl and microtype, pdfinfo (poppler) and zip.
set -euo pipefail

version=$1
mkdir -p "$2"
out=$(cd "$2" && pwd)
root=$(cd "$(dirname "$0")/../.." && pwd)
src=$root/paper
pdf=$root/paper/paper1/MirrorClauseGSC.pdf
name=MirrorClauseGSC-v$version
stamp=$(git -C "$root" log -1 --format=%cI)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

mkdir -p "$work/tex"
cp "$src/main.tex" "$work/tex/"

# The source has to compile from this file alone, with no BibTeX pass and no shell escape, to the
# same number of pages as the committed PDF, with no undefined reference and no bad box.
cp -r "$work/tex" "$work/build"
(
  cd "$work/build"
  for pass in 1 2 3; do
    pdflatex -interaction=nonstopmode -halt-on-error -no-shell-escape main.tex > /dev/null
  done
  pages() { pdfinfo "$1" | awk '/^Pages:/ {print $2}'; }
  test "$(pages main.pdf)" = "$(pages "$pdf")"
  if grep -iE 'undefined|overfull|underfull' main.log; then exit 1; fi
)

find "$work/tex" -exec touch -h -d "$stamp" {} +
cp "$pdf" "$out/$name.pdf"
(cd "$work/tex" && find . -type f | LC_ALL=C sort | sed 's|^\./||' | TZ=UTC zip -X -q -9 "$out/$name-tex.zip" -@)
ls -l "$out"
