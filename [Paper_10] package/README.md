# Explicit String Vacua and the Geometric Selection Conjecture

Paper 10 - research revision v2.0.0, 25 October 2026.

**The full GSC remains unproved. No single heterotic vacuum satisfying
all nine conditions is supplied. These files must not be described as an
unconditional resolution or an accepted Annals paper.**

The two editions use the same mathematical body:

- `Annals-review-draft/`: a PDF and a source ZIP using the `aomart` class.
- `Hexagon-v2-review-draft/`: a PDF and a source ZIP with root `main.tex`
  using `amsart`, suitable for a compilation review.
- `verification/`: exact and symbolic checks, proof-status ledger and
  compilation/visual-review record.
- `AI_PROVENANCE.md`, `hexagonsubmissionform.md` and
  `journalsubmissiondetails.md`: scope and proposed submission information.

There are no margin line numbers. Ordinary page, theorem and equation
numbers remain. No new Hexagon version is asserted to exist.

The proposed revision target is the existing GSC note,
https://hexagonmath.org/2610.00493v1, with the Dwork paper,
https://hexagonmath.org/2610.00446v1, cited separately. This package does
not overwrite the Dwork paper or its previous submission files.

What is established in the draft:

1. A marked compatible point gives a smooth globally proper invariant
   selection function under the stated manifold/involution hypotheses.
2. The specified tetraquadric bundle has the stated slopes, Chern classes,
   effective anomaly class on the cover and a degree-four quotient index
   of three in absolute value. Quotient anomaly torsion remains unaudited.
3. A positive matter metric can erase a holomorphic Yukawa hierarchy;
   uniform positive upper/lower bounds preserve its asymptotic orders.
4. Explicit Hessian and residual bounds imply a unique nearby stable
   critical point and preserve strictly separated mass inequalities.
5. A specified classical global-supersymmetry effective model has positive
   scalar Hessian and broken supersymmetry. Its string realization is not
   established. A no-scale calculation records a tree-level obstruction.

These are elementary construction criteria and independent checks of a
published example, not a newly established breakthrough. Definitions are
credited to the existing preprints; no unrefereed existence theorem is
used to certify a viable vacuum.

The source ZIPs compile with pdfLaTeX. Run it three times on `main.tex`.
The optional `anc/` scripts do not enter TeX compilation. Exact controls
need Python 3. Symbolic controls additionally need Sympy 1.14.0.

Author review of all new material is pending. This research revision is
published through the existing GSC repository as v2.0.0. Its Zenodo version
DOI is recorded in the release and receipt once issued. No journal
submission, acceptance or automatic Hexagon upload is claimed.
