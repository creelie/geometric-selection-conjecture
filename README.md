# Paper 10 — Explicit String Vacua and the Geometric Selection Conjecture

**Unpublished author-review revision. The full GSC is not proved.**

Author of this revision: Deep Bhattacharjee. Ushashi Bhattacharya is credited in the manuscript footnote and the original coauthored references.

This is one coherent manuscript consolidating the two supplied Hexagon TeX
archives and Paper 10, with a new fixed-complex-structure construction, an analytic singlet obstruction, and
an explicit audit of the nine conditions. No outside paper PDF is embedded.
There is one title page and one bibliography, and no margin line numbering.

## What was added

- An explicit invariant pencil `F_theta = F_0 + theta G`, fixed at `theta = pi`.
- A stratified proof that `F_0 = 0` is smooth, with rational Bezout certificates.
- Smoothness of `X_pi` and freeness of the degree-four action, using a rational
  discriminant and transcendence of pi, rather than numerical sampling.
- Exact equivariant Koszul cohomology for the published tetraquadric bundle.
  Four six-dimensional integer character blocks have determinants
  `792, -792, 792, -792` at theta=1. Their nonzero product certifies the
  nonzero determinant polynomial and hence invertibility at theta=pi.
- Three-family charged spectrum, one Higgs doublet pair and no massless
  colored Higgs triplets for the specified Wilson line.
- An explicit residual Abelian generator, perturbative Yukawa rank obstruction,
  classical flat volume direction, and a targeted holomorphic extension route.
- A condition-by-condition ledger for this fixed datum, not a combination
  of results at different vacua.

## Actual result of the nine-condition audit

| Condition | Result for the fixed split-bundle datum |
| --- | --- |
| C1 | Polystability and cover anomaly effectivity proved; integral quotient anomaly and hidden completion not certified. |
| C2 | Exact descended index -3. |
| C3 | Extra U(1) survives in the specified perturbative model. |
| C4 | Three charged families and no additional chiral charged exotics. |
| C5 | No additional massless charged vectorlike exotics to lift; neutral/hidden masses not asserted. |
| C6 | One Higgs doublet pair, no massless Higgs triplets in the projection. |
| C7 | Classical volume direction is flat; controlled corrected stabilization absent. |
| C8 | Perturbatively rank(Y_u) <= 2 and Y_d = 0 at zero singlet expectation values. |
| C9 | No verified supersymmetry-breaking sector and phenomenological certificate. |

This candidate **does not pass all nine conditions**. The non-Abelian extension
construction is a different branch requiring a fresh stability, operator and
potential audit. It is not silently substituted into the fixed-point ledger.

## Files and reproduction

The PDF is the review copy. In either edition folder's `tex-source`, compile `main.tex` with pdfLaTeX
three times. The class, body, four matrix certificates and all figures are
included. No shell escape, external article PDF, network access, bibliography
download or Lean installation is required for TeX compilation.

The `anc` directory contains executable controls and their exact JSON results.
Use Python with Sympy 1.14.0. Run `smoothness.py`, `spectrum.py`,
`matrix_certificate.py`, and `controls.py`. The matrix generator refreshes
`matrix-certificates.tex` from integer input. These are finite algebra checks;
they are not a complete formalization or a proof of physical viability.

The consolidated spectral proof text is retained from the supplied Dwork
source. This session rechecked the stated finite counts and arithmetic, not
every inherited analytic assertion or every old repository toolchain.
Independent expert review and novelty assessment are pending.

## Sources and publication status

- https://hexagonmath.org/2610.00446v1
- https://hexagonmath.org/2610.00493v1
- Paper 10 research revision: https://doi.org/10.5281/zenodo.23646340
- Bundle data credited to Buchbinder, Constantin and Lukas,
  JHEP 03 (2014) 025, https://doi.org/10.1007/JHEP03(2014)025

The existing Paper 10 DOI is a source citation, **not a DOI for Paper 10**.
No new release, DOI, public Hexagon version or journal submission has been made.
No final ZIP is created before a new DOI, as previously requested by the author.

Mathematical and physical closure is not implied by use of the Annals class.


Current files: [`[Paper_10] package`](%5BPaper_10%5D%20package).
Earlier immutable release: [v2.0.0](https://github.com/creelie/geometric-selection-conjecture/releases/tag/v2.0.0).
