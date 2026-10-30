# Explicit String Vacua and the Geometric Selection Conjecture

Deep Bhattacharjee and Ushashi Bhattacharya. Release v3.3.0 of the paper, its
LaTeX source and the programs that check its finite arithmetic. Contact for
this repository: egspl@proton.me (GitHub account creelie).

**The Geometric Selection Conjecture is not proved.** The paper exhibits one
explicit heterotic datum on a free quotient of a tetraquadric, audits it
against the nine conditions C1 to C9, and proves that the two failures it finds
at the split locus, a surviving Abelian gauge factor and a vanishing Yukawa
matrix, are forced. On the two-block extension branch over the quotient the
first six conditions hold together, while C8 fails for every bundle on the
quotient built from the eleven data by extensions, with blocks of any size.
No result depends on an unrefereed source.

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
- Section 16.3. The extension branch of the eleven data
  is decided for its two smallest shapes. A single rank-two extension never
  removes the forced Abelian factor, at any Kahler class (Proposition 16.8).
  Two rank-two extensions remove it for six of the eleven data (the fixed
  datum and classes 2, 6, 7, 8, 10), each realised on one side of a stability
  wall through an explicit polystable class (Theorem 16.10, with the wall
  lemma 16.7 and exact extension groups), and are excluded for classes 4, 5, 9
  at every Kahler class and for class 1 on the general tetraquadric (Lemma
  16.9).
- Section 16.4. The branch on the quotient. The Gamma-
  compatible data are the fixed datum and classes 6, 8, 9, 10; class 9 has no
  polystable two-block bundle. By holomorphic Lefschetz the cohomology of each
  summand is a multiple of the regular representation (Lemma 16.12), and
  selection rules from the block torus and from closed subbundles (Lemma
  16.14) bound the cup-product Yukawas for every extension class and every
  equivariant structure (Theorem 16.15): rank Y_u <= 2 and Y_d = Y_e = 0 for
  the fixed datum, rank Y_u <= 1 for class 6, rank Y_u <= 2 for class 8, and
  Y_u = Y_d = Y_e = 0 for class 10. Hence C8_W fails for every polystable
  Gamma-bundle with two blocks of size two from the eleven data (Corollary
  16.16). For the fixed datum on the shape L2->L1, L5->L3, L4 the bundle E_Y
  is polystable on the quotient, has no massless extra U(1), and has exactly
  three families, one Higgs pair and no triplets (Theorem 16.17): C1 to C6 hold
  together, C7 is flat and C8 fails. Worldsheet instantons are not treated.
- Section 16.5 (new in this version). Bundles filtered by the line bundles,
  with blocks of any size. An extension of L_b by L_a can be reordered when
  H^1(L_a L_b^-1) = 0 (Lemma 16.20), so every such bundle has a subbundle on
  which the selection rules act. Theorem 16.21: Y_d = Y_e = 0 for the fixed
  datum; rank Y_u <= 2, or rank Y_d, Y_e <= 2, for class 6; rank Y_u <= 2 for
  class 8; Y_u = Y_d = Y_e = 0 for classes 9 and 10. Corollary 16.22: C8_W
  fails for every bundle on X_pi/Gamma filtered by the line bundles of the
  eleven data that has three families.
- The eleven data are checked independently (first Chern class, index, c_2
  bound, H^2(V) = 0 and an exact polystable Kahler class for each), and the
  five unrefereed preprints cited before are no longer used.
- Definition 9.1 fixes acceptance windows W for the
  three phenomenological conditions C7, C8, C9, with the reason for each, and
  Conjecture 9.2 states GSC(W). Every negative result of the paper holds for
  every choice of windows (Remark 9.3).

| Condition | Fixed datum |
| --- | --- |
| C1 | Polystability, c_1 = 0 and cover anomaly effectivity proved, at the split point and for E_Y; quotient anomaly and hidden completion not certified. |
| C2 | Descended index exactly -3, for E_Y as well. |
| C3 | Fails at the split locus, as Theorem 16.1 forces for every rank-five split datum; holds for E_Y on the quotient (Theorem 16.17). |
| C4 | Three charged families, no additional chiral charged exotics, for E_Y as well. |
| C5 | No additional massless charged vectorlike exotics; neutral and hidden masses not asserted. |
| C6 | One Higgs doublet pair, no massless Higgs triplets. |
| C7 | The volume ray is flat at tree level, also on the branch: a massless modulus, so C7_W fails for every window. |
| C8 | Fails: rank Y_u <= 2 and Y_d = 0, and no singlet insertion restores Y_d (Proposition 16.4); on the branch rank Y_u <= 2 and Y_d = Y_e = 0 for every extension class (Theorem 16.15), and Y_d = Y_e = 0 for every bundle filtered by the same line bundles (Theorem 16.21); C8_W fails for every window. |
| C9 | No supersymmetry-breaking sector constructed; the F- and D-terms vanish, so C9_W fails for every window. |

## Layout

- `final-revision/`: the paper (`main.tex`, `paper-body.tex`, `split-branch.tex`,
  `matrix-certificates.tex`, `figures/`), its compiled PDF
  `ExplicitStringVacuaGSC.pdf` (amsart, 56 pages), and `anc/`, the exact
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
- `code/gamma_branch.py`: Section 16.4 and the independent check of the eleven
  data: first Chern class, index, c_2 bound, H^2(V) = 0 and an exact
  polystable Kahler class for each; equivariant cohomology on X_pi over
  Q(theta), the Gamma-compatible data and shapes, and the term-rank bounds of
  Theorem 16.15 over all equivariant structures (152 checks; output
  `gamma_branch_results.json`).
- `code/filtered_branch.py`: Section 16.5: the links H^1(L_a L_b^-1) != 0 and
  the pieces on X_pi for the five Gamma-data, the facts used in the proof of
  Theorem 16.21, and its bounds recomputed for all 120 flag orders and every
  down-set (644 checks; output `filtered_branch_results.json`).
- `code/run_checks.sh` runs the five checks; `bash .github/release/build.sh <version> dist`
  rebuilds the release assets.
- `paper/`: the mirror-clause note of releases 1.x. `[Paper_10] package/` and
  `archive/`: the review drafts of releases 2.x.

Compile the paper with `pdflatex main.tex` three times in `final-revision/`;
no BibTeX pass and no shell escape are needed. The ten further line-bundle
sums of Section 16 were found by a computer search; `code/gamma_branch.py`
checks every property the paper uses, so nothing rests on that search.

## Earlier versions

- Hexagon 2610.00446v1 (spectral calculations), 2610.00493v1 (mirror clause)
  and 2610.00493v2; this version is 2610.00493v3.
- Releases: v1.0.1 mirror-clause note, 10.5281/zenodo.23567115; v2.0.0,
  10.5281/zenodo.23646340; v2.1.0, 10.5281/zenodo.23647986.
- Release v3.0.0 (split branch): https://doi.org/10.5281/zenodo.23682773.
- Release v3.1.0 (extension branch): https://doi.org/10.5281/zenodo.23688897.
- Release v3.2.0 (branch on the quotient): https://doi.org/10.5281/zenodo.23788190.
- Release v3.3.0 (this version): the version DOI is recorded after Zenodo
  publishes it. Concept DOI of all versions: 10.5281/zenodo.23567114.
