# The Mirror Clause of the Geometric Selection Conjecture Is Empty

Deep Bhattacharjee and Ushashi Bhattacharya

This repository holds the note, its LaTeX source and the check for *The Mirror
Clause of the Geometric Selection Conjecture Is Empty*. The compiled note is
[`paper/paper1/MirrorClauseGSC.pdf`](paper/paper1/MirrorClauseGSC.pdf) and its
source, as one zip, is
[`paper/paper2/MirrorClauseGSC-tex.zip`](paper/paper2/MirrorClauseGSC-tex.zip);
[`paper/main.tex`](paper/main.tex) is the same file unzipped.

## What is proved

The Geometric Selection Conjecture of Bhattacharjee, Bhattacharya and Ghosh
asks for a functional on the heterotic landscape, built from the
Weil–Petersson metric and invariant under mirror symmetry, whose unique
minimum is a phenomenologically viable vacuum (conditions (C1)–(C9)). The
paper *Which String Vibrations Become Which Generations on the Dwork Quintic*
([hexagonmath.org 2610.00446v1](https://hexagonmath.org/2610.00446v1), code at
[creelie/strings v1.2.0](https://github.com/creelie/strings/releases/tag/v1.2.0))
proved that the conjecture holds if and only if there is a viable vacuum fixed
by the mirror involution μ wherever μ is defined.

The note proves that the fixed-point clause is never met.

- A vacuum fixed by μ is isomorphic to its mirror, so its threefold has
  h¹¹ = h²¹ and Euler number χ = 0.
- Its component of the landscape contains a standard embedding (X₀, T_{X₀}),
  and Chern classes and Euler numbers are constant on a component, so
  c₃(V) = χ = 0.
- By Hirzebruch–Riemann–Roch, χ(X, V) = c₃(V)/2 for c₁(V) = 0 on a
  Calabi–Yau threefold, so the vacuum has zero net families and fails (C2).

Hence the conjecture holds if and only if a viable vacuum exists on a
component of the landscape that carries no standard embedding, with no
condition on mirror symmetry; the same holds when invariance under a compact
group of isometries containing μ is demanded; and no vacuum built on the
tangent bundle of any Calabi–Yau threefold, or a deformation of it, is ever
the selected vacuum, three families or not (the Tian–Yau quotient has three).
If mirror pairs are identified as one point, the conjecture is equivalent to
the existence of a viable vacuum of any type.

What the note does not do is exhibit a viable vacuum. Conditions (C7)–(C9),
stabilised moduli, the Yukawa hierarchy and supersymmetry breaking, are fixed
by comparison with measurement, and no published model establishes all three.
The conjecture stays open there.

## The check

`code/mirror_clause.py` carries the lemma and the corollary as functions on the
discrete data of a vacuum (Hodge numbers, c₃, order of a free quotient) and
checks them on self-mirror pairs such as (19, 19), (11, 11), (3, 3) and
three-family pairs such as (6, 9), with refuters; `bash code/run_checks.sh`
prints `511 checks passed, 0 failed`. The same file is the module
`construct/gsc_mirror_fixed.py` of the Anaxia engine,
[creelie/eBLLM](https://github.com/creelie/eBLLM).

## Releases

The note cites this repository by its release page,
https://github.com/creelie/geometric-selection-conjecture/releases/tag/v1.0.0.
Merging a change to `.github/release/notes.md` into `main` publishes the
release for the version in `CITATION.cff`, with the PDF and the source zip
attached.

## License

MIT, see [`LICENSE`](LICENSE).
