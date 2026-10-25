#!/usr/bin/env bash
# Build the release assets of "The Mirror Clause of the Geometric Selection Conjecture Is Empty":
# the PDF and the LaTeX source (tex.zip).
#
# Usage: build.sh <version> <output directory>
# Needs pdflatex with amsart, xurl and microtype, pdfinfo (poppler) and zip.
set -euo pipefail

version=$1
root=$(cd "$(dirname "$0")/../.." && pwd)
if [ -d "$root/[Paper_10] package" ]; then
  mkdir -p "$2"
  out=$(cd "$2" && pwd)
  for edition in Annals-review-draft Hexagon-v2-review-draft; do
    src="$root/[Paper_10] package/$edition"
    work=$(mktemp -d)
    cp -r "$src/tex-source/." "$work/"
    (
      cd "$work"
      for pass in 1 2 3; do
        pdflatex -interaction=nonstopmode -halt-on-error -no-shell-escape main.tex > /dev/null
      done
      pages() { pdfinfo "$1" | awk '/^Pages:/ {print $2}'; }
      test "$(pages main.pdf)" = "$(pages "$src/ExplicitStringVacuaGSC.pdf")"
      if grep -iE 'undefined|overfull' main.log; then exit 1; fi
    )
    if [ "${CHECK_ONLY:-0}" != 1 ]; then
      test -f "$root/[Paper_10] package/DOI-RECEIPT.json" || {
        echo 'Obtain and record the Zenodo version DOI before creating final ZIPs.' >&2
        exit 1
      }
      python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["version"]==sys.argv[2]' "$root/[Paper_10] package/DOI-RECEIPT.json" "$version"
      cp "$src/ExplicitStringVacuaGSC.pdf" "$out/ExplicitStringVacuaGSC-$edition-v$version.pdf"
      (cd "$src/tex-source" && find . -type f | LC_ALL=C sort | sed 's|^./||' | zip -X -q -9 "$out/ExplicitStringVacuaGSC-$edition-v$version-tex.zip" -@)
    fi
  done
  exit 0
fi
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
