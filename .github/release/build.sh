#!/usr/bin/env bash
# Build the release assets of "Explicit String Vacua and the Geometric Selection Conjecture":
# the PDF, the LaTeX source (tex.zip) and the arXiv tarball, from final-revision/.
#
# Usage: build.sh <version> <output directory>
# Needs pdflatex with amsart, xurl, adjustbox, orcidlink and microtype, pdfinfo (poppler), zip and tar.
set -euo pipefail

version=$1
mkdir -p "$2"
out=$(cd "$2" && pwd)
root=$(cd "$(dirname "$0")/../.." && pwd)
src=$root/final-revision
pdf=$src/ExplicitStringVacuaGSC.pdf
name=ExplicitStringVacuaGSC-v$version
stamp=$(git -C "$root" log -1 --format=%cI)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

mkdir -p "$work/tex"
cp "$src"/main.tex "$src"/paper-body.tex "$src"/split-branch.tex "$src"/matrix-certificates.tex "$work/tex/"
cp -r "$src/figures" "$src/anc" "$work/tex/"
find "$work/tex" -name '__pycache__' -prune -exec rm -rf {} +

# The source has to compile from these files alone, with no BibTeX pass and no shell escape, to the
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

if [ "${CHECK_ONLY:-0}" = 1 ]; then
  echo "final-revision compiles to $(pdfinfo "$pdf" | awk '/^Pages:/ {print $2}') pages with no bad box."
  exit 0
fi

find "$work/tex" -exec touch -h -d "$stamp" {} +
cp "$pdf" "$out/$name.pdf"
(cd "$work/tex" && find . -type f | LC_ALL=C sort | sed 's|^\./||' | TZ=UTC zip -X -q -9 "$out/$name-tex.zip" -@)
(cd "$work/tex" && find . -type f | LC_ALL=C sort | sed 's|^\./||' | tar --mtime="$stamp" --owner=0 --group=0 --numeric-owner -czf "$out/$name-arxiv.tar.gz" -T -)
ls -l "$out"
