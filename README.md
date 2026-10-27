# Explicit String Vacua and the Geometric Selection Conjecture

Deep Bhattacharjee and Ushashi Bhattacharya. Release v3.0.0 of the paper, its
LaTeX source and the programs that check its finite arithmetic, archived at
https://doi.org/10.5281/zenodo.23682773.

**The Geometric Selection Conjecture is not proved.** The paper exhibits one
explicit heterotic datum on a free quotient of a tetraquadric, audits it
against the nine conditions C1 to C9, and proves that the two failures it finds
at the split locus, a surviving Abelian gauge factor and a vanishing Yukawa
matrix, are forced.

## What the paper proves

- The character labels, cubic couplings and invariant-direction mass maps of
  the Dwork and weighted Fermat pencils, and the exclusion of three
  standard-embedding families on free quotients of thirteen hypergeometric
  complete-intersection types.
- The mirror clause of the conjecture is empty, so the conjecture is the
  existence of one viable vacuum on a component without a standard embedding.
- One tetraquadric equation at a transcendental parameter is smooth with a
  free action of order four; the published rank-five line-bundle sum of
  Buchbinder, Constantin and Lukas on it has index -3 and, after the Wilson
  line, three charged families and one Higgs doublet pair, by an exact
  determinant certificate.
- Section 16 (new in this version). On every Calabi-Yau threefold with at
  most four Kahler parameters, every polystable sum of five line bundles with
  c_1 = 0 has charge matrix of rank at most three and keeps at least one
  Abelian gauge boson massless at its split locus (Theorem 16.1): condition C3
  fails for every rank-five split datum there. For the fixed bundle the
  polystable Kahler classes are exactly the volume ray. The perturbative
  selection rules, decided by a flow criterion (Lemma 16.3), exclude six of
  eleven three-family split data on the tetraquadric from C8 to all orders in
  the singlets and leave five (Proposition 16.4); the five are candidates for
  the extension branch, not vacua.

| Condition | Fixed datum |
| --- | --- |
| C1 | Polystability, c_1 = 0 and cover anomaly effectivity proved; quotient anomaly and hidden completion not certified. |
| C2 | Descended index exactly -3. |
| C3 | Fails: an Abelian generator survives, and Theorem 16.1 forces this for every rank-five split datum. |
| C4 | Three charged families, no additional chiral charged exotics. |
| C5 | No additional massless charged vectorlike exotics; neutral and hidden masses not asserted. |
| C6 | One Higgs doublet pair, no massless Higgs triplets. |
| C7 | The volume ray is flat at tree level. |
| C8 | Fails: rank Y_u <= 2 and Y_d = 0, and no singlet insertion restores Y_d (Proposition 16.4). |
| C9 | No supersymmetry-breaking sector constructed. |

## Layout

- `final-revision/`: the paper (`main.tex`, `paper-body.tex`, `split-branch.tex`,
  `matrix-certificates.tex`, `figures/`), its compiled PDF
  `ExplicitStringVacuaGSC.pdf` (amsart, 44 pages), and `anc/`, the exact
  controls of the tetraquadric construction with their JSON results.
- `code/mirror_clause.py`: the mirror clause on Hodge data (511 checks).
- `code/yukawa_charges.py`: the selection rules of Section 16 on the eleven
  split data, by a maximal flow, with the closed sets or singlet insertions
  that certify each entry of Table 6 (45 checks; output
  `yukawa_charges_results.json`). `code/fig_charges.py` draws Figure 5.
- `code/run_checks.sh` runs both checks; `bash .github/release/build.sh <version> dist`
  rebuilds the release assets.
- `paper/`: the mirror-clause note of releases 1.x. `[Paper_10] package/` and
  `archive/`: the review drafts of releases 2.x.

Compile the paper with `pdflatex main.tex` three times in `final-revision/`;
no BibTeX pass and no shell escape are needed. The ten further line-bundle
sums of Section 16 come from the module `construct/lbs_tetraquadric.py` of
the Anaxia engine, https://github.com/creelie/eBLLM, at commit `a85ae17`.

## Earlier versions

- Hexagon 2610.00446v1 (spectral calculations) and 2610.00493v1 (mirror clause).
- Releases: v1.0.1 mirror-clause note, 10.5281/zenodo.23567115; v2.0.0,
  10.5281/zenodo.23646340; v2.1.0, 10.5281/zenodo.23647986.
- Release v3.0.0 (this version): https://doi.org/10.5281/zenodo.23682773. Concept DOI of all
  versions: 10.5281/zenodo.23567114.
