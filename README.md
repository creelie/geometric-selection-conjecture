# Explicit String Vacua and the Geometric Selection Conjecture

Deep Bhattacharjee and Ushashi Bhattacharya. Release v3.1.0 of the paper, its
LaTeX source and the programs that check its finite arithmetic, archived at
https://doi.org/10.5281/zenodo.23688897.

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
- Section 16. On every Calabi-Yau threefold with at
  most four Kahler parameters, every polystable sum of five line bundles with
  c_1 = 0 has charge matrix of rank at most three and keeps at least one
  Abelian gauge boson massless at its split locus (Theorem 16.1): condition C3
  fails for every rank-five split datum there. For the fixed bundle the
  polystable Kahler classes are exactly the volume ray. The perturbative
  selection rules, decided by a flow criterion (Lemma 16.3), exclude six of
  eleven three-family split data on the tetraquadric from C8 to all orders in
  the singlets and leave five (Proposition 16.4); the five are candidates for
  the extension branch, not vacua.
- Section 16.3 (new in this version). The extension branch of the eleven data
  is decided for its two smallest shapes. A single rank-two extension never
  removes the forced Abelian factor, at any Kahler class (Proposition 16.8).
  Two rank-two extensions remove it for six of the eleven data (the fixed
  datum and classes 2, 6, 7, 8, 10), each realised on one side of a stability
  wall through an explicit polystable class (Theorem 16.10, with the wall
  lemma 16.7 and exact extension groups), and are excluded for classes 4, 5, 9
  at every Kahler class and for class 1 on the general tetraquadric (Lemma
  16.9). Descent to the quotient, spectrum and Yukawa couplings of the
  extensions are not recomputed.
- Definition 9.1 (new in this version) fixes acceptance windows W for the
  three phenomenological conditions C7, C8, C9, with the reason for each, and
  Conjecture 9.2 states GSC(W). Every negative result of the paper holds for
  every choice of windows (Remark 9.3).

| Condition | Fixed datum |
| --- | --- |
| C1 | Polystability, c_1 = 0 and cover anomaly effectivity proved; quotient anomaly and hidden completion not certified. |
| C2 | Descended index exactly -3. |
| C3 | Fails at the split locus, and Theorem 16.1 forces this for every rank-five split datum; on the two-block extension branch the generator is removed on the cover (Theorem 16.10), descent not certified. |
| C4 | Three charged families, no additional chiral charged exotics. |
| C5 | No additional massless charged vectorlike exotics; neutral and hidden masses not asserted. |
| C6 | One Higgs doublet pair, no massless Higgs triplets. |
| C7 | The volume ray is flat at tree level: a massless modulus, so C7_W fails for every window. |
| C8 | Fails: rank Y_u <= 2 and Y_d = 0, and no singlet insertion restores Y_d (Proposition 16.4); C8_W fails for every window. |
| C9 | No supersymmetry-breaking sector constructed; the F- and D-terms vanish, so C9_W fails for every window. |

## Layout

- `final-revision/`: the paper (`main.tex`, `paper-body.tex`, `split-branch.tex`,
  `matrix-certificates.tex`, `figures/`), its compiled PDF
  `ExplicitStringVacuaGSC.pdf` (amsart, 49 pages), and `anc/`, the exact
  controls of the tetraquadric construction with their JSON results.
- `code/mirror_clause.py`: the mirror clause on Hodge data (511 checks).
- `code/yukawa_charges.py`: the selection rules of Section 16 on the eleven
  split data, by a maximal flow, with the closed sets or singlet insertions
  that certify each entry of Table 6 (45 checks; output
  `yukawa_charges_results.json`). `code/fig_charges.py` draws Figure 5.
- `code/extension_branch.py`: the extension branch of Section 16.3 in exact
  arithmetic (sympy): slope spaces, orientations, Koszul extension groups,
  polystable points and wall certificates, and the permutation certificate of
  Lemma 16.9 (61 checks; output `extension_branch_results.json`).
- `code/run_checks.sh` runs the three checks; `bash .github/release/build.sh <version> dist`
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
- Release v3.0.0 (split branch): https://doi.org/10.5281/zenodo.23682773.
- Release v3.1.0 (this version): https://doi.org/10.5281/zenodo.23688897. Concept DOI of all
  versions: 10.5281/zenodo.23567114.
