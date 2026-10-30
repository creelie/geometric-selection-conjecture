# Explicit String Vacua and the Geometric Selection Conjecture, release 3.2.0

The paper by Deep Bhattacharjee and Ushashi Bhattacharya in its amsart form (54 pages), with its LaTeX source, the arXiv tarball and the checking programs. Contact for this repository: egspl@proton.me (GitHub account creelie).

- New Section 16.4, "The branch on the quotient". Lemma 16.12: by holomorphic Lefschetz, cohomology of a summand in a single degree is a multiple of the regular representation of Γ = (Z/2)², whatever the equivariant structure. Lemma 16.13 and Tables 10, 11: the Γ-compatible data are the fixed datum and classes 6, 8, 9, 10, and class 9 has no polystable two-block bundle. Lemma 16.14: selection rules from the block torus and from closed subbundles. Theorem 16.15: for every extension class and every equivariant structure, rank Y_u ≤ 2 and Y_d = Y_e = 0 for the fixed datum, rank Y_u ≤ 1 for class 6, rank Y_u ≤ 2 for class 8, and Y_u = Y_d = Y_e = 0 for class 10. Corollary 16.16: C8_W fails for every polystable Γ-bundle with two blocks of size two built from the eleven data.
- Theorem 16.17: for the fixed datum on the shape L2→L1, L5→L3, L4 the descended bundle E_Y is polystable, has no massless Abelian gauge boson beyond SU(5), and has exactly three families, one Higgs pair and no Higgs triplets. So C1 to C6 hold together on the branch, C7 is flat under J ↦ λJ and C8 fails. Worldsheet instantons are not treated (Remark 16.18).
- The eleven three-family data are checked independently (first Chern class, index, c₂ bound, H²(V) = 0 and an exact polystable Kähler class for each). The five unrefereed preprints cited before are no longer used; every remaining reference is refereed.
- Corrections: condition C5 is read as the absence of charged vectorlike exotics, in line with GSC(W); a sign in a figure caption, the ledger (Table 12) and the conclusion are brought in line with the results.
- `code/gamma_branch.py` (152 checks, output `gamma_branch_results.json`) checks the eleven data and Section 16.4: equivariant cohomology on X_π over Q(θ), the Γ-compatible shapes and the term-rank bounds over all equivariant structures. `code/extension_branch.py` (61), `code/yukawa_charges.py` (45) and `code/mirror_clause.py` (511) are unchanged.
- The Geometric Selection Conjecture is not proved.

Assets: the PDF, the LaTeX source (tex.zip) and the arXiv tarball.
