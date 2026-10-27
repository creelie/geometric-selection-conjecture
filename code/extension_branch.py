"""The extension branch of the eleven split data on a tetraquadric X in (P^1)^4 (Section 16.3).

For V = L_1 + ... + L_5 with charge matrix K (rows k_a in Z^4) the slope vector
mu(J) = (mu_J(L_a))_a = K kappa(J), kappa_i(J) = J_i J^2, lies in the column space Im K of K, which has
dimension rk K (three for all data but class 3) and is contained in the hyperplane sum_a mu_a = 0.

Two rank-two blocks.  A bundle E_1 + E_2 + L_o with E_1 a non-split extension 0 -> L_p -> E_1 -> L_q -> 0
(class in H^1(L_p L_q^-1)) and E_2 one of L_s by L_r is polystable of slope zero only if mu_o = 0,
mu_p + mu_q = 0 (hence mu_r + mu_s = 0), mu_p < 0 < mu_q and mu_r < 0 < mu_s.  The slopes then lie in
W = {mu in Im K : mu_o = 0, mu_p + mu_q = 0}.  When W is a line Q w the signs are those of f w for a
scalar f, so an orientation (which member of each pair is the sub-line-bundle) is forced by the sign of
f; when W is a plane the achievable orientations are read off from the two functionals mu_q, mu_s on W.
A block structure is REALISABLE when for some achievable orientation both extension groups
H^1(L_sub L_quot^-1) are nonzero on every smooth tetraquadric, and EXCLUDED when for every achievable
orientation one of them vanishes on every smooth tetraquadric (the Koszul sequence gives the exact value
whenever h^1(A, k - 2) = 0 = h^2(A, k) on A = (P^1)^4).  Otherwise it is UNDECIDED (only upper bounds).
The forced Abelian generator u (K u = 0, sum u = 0) is Higgsed by the extension classes exactly when u is
not constant on both pairs; for a line W this holds for one pair iff it holds for the other, by u . w = 0.

One rank-two block.  E_1 + L_r + L_s + L_o is polystable only if mu_r = mu_s = mu_o = 0.  When the rows
k_r, k_s, k_o are independent this forces kappa into the line K^T kappa = 0, so mu = 0 and E_1 cannot be
stable: no single extension appears at any Kaehler class.

For every datum with a realisable structure the script certifies the polystable point J_0 of the split
bundle (all slopes zero, J_0 in the open Kaehler cone), that the two constraint gradients are independent
there (the surface S is smooth at J_0), and that mu_q has nonzero differential along S (f takes both signs
on S near J_0).  Exact arithmetic throughout (sympy).  Output: extension_branch_results.json.
"""
import itertools, json, os, sys
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # yukawa_charges.py beside this file, also under python3 -I
from yukawa_charges import BCL, ENGINE, hA, sub as vsub, rank

t = sp.symbols('t1:5')
r_sym = sp.Symbol('r')

def kappa(i):  # J_i J^2 / 4 on the tetraquadric: e_2 of the other three coordinates
    idx = [j for j in range(4) if j != i]
    return sum(t[a] * t[b] for a, b in itertools.combinations(idx, 2))

def slope(k):
    return sp.expand(sum(k[i] * kappa(i) for i in range(4)))

def h1(k):
    """Exact h^1(X, O_X(k)) on every smooth tetraquadric when the Koszul connecting maps vanish for
    degree reasons (h^1(A, k-2) = 0 = h^2(A, k)); otherwise None (only an upper bound is known)."""
    a = hA(k); b = hA(tuple(x - 2 for x in k))
    if a[1] + b[2] == 0: return 0          # the Koszul upper bound vanishes
    return a[1] + b[2] if (b[1] == 0 and a[2] == 0) else None

# Cech bases on P^1: H^0(O(m)) has the monomials x^a y^b with a, b >= 0, a + b = m; H^1(O(m)) those with
# a, b <= -1, a + b = m.  Multiplication by a polynomial acts on these representatives, and a product whose
# exponent in x or y leaves the range is a coboundary.  The Kuenneth basis of H^d(A, O(k)) is the product.
def cech_basis(m, i):
    if i == 0: return [(a, m - a) for a in range(0, m + 1)] if m >= 0 else []
    return [(a, m - a) for a in range(-1, m, -1) if m - a <= -1] if m <= -2 else []

def koszul_matrix(k, degs, quadric):
    """Matrix of multiplication by the (2,2,2,2) form `quadric` (list of (coefficient, exponent 4-tuple of pairs))
    from the Kuenneth component H^d(A, O(k - 2)) with factor degrees `degs` to the same component of H^d(A, O(k))."""
    src = list(itertools.product(*[cech_basis(k[j] - 2, degs[j]) for j in range(4)]))
    tgt = list(itertools.product(*[cech_basis(k[j], degs[j]) for j in range(4)]))
    idx = {m: i for i, m in enumerate(tgt)}
    M = sp.zeros(len(tgt), len(src))
    for j, mono in enumerate(src):
        for coeff, exps in quadric:
            new = tuple((mono[f][0] + exps[f][0], mono[f][1] + exps[f][1]) for f in range(4))
            if all((new[f][0] >= 0 and new[f][1] >= 0) if degs[f] == 0 else (new[f][0] <= -1 and new[f][1] <= -1) for f in range(4)):
                M[idx[new], j] += coeff
    return M

# The quadric s_0 of Section 16.3: six monomials, written as exponent pairs (power of x_i, power of y_i).
S0 = [(1, ((2, 0), (2, 0), (0, 2), (2, 0))), (1, ((0, 2), (0, 2), (2, 0), (2, 0))), (1, ((0, 2), (2, 0), (2, 0), (0, 2))),
      (1, ((2, 0), (0, 2), (0, 2), (0, 2))), (1, ((0, 2), (1, 1), (1, 1), (0, 2))), (1, ((2, 0), (1, 1), (1, 1), (2, 0)))]

def traceless_kernel(K):
    M = sp.Matrix([[K[a][i] for a in range(5)] for i in range(4)] + [[1] * 5])
    out = []
    for v in M.nullspace():
        v = v * sp.ilcm(*[sp.fraction(x)[1] for x in v])
        g = sp.igcd(*[int(x) for x in v]); v = v / g
        if v[0] < 0 or (v[0] == 0 and sum(v) < 0): v = -v
        out.append([int(x) for x in v])
    return out

def primitive(v):
    v = sp.Matrix(v); v = v * sp.ilcm(*[sp.fraction(x)[1] for x in v]); g = sp.igcd(*[int(x) for x in v])
    return [int(x) for x in v / g]

def space_W(K, p, q, o):
    """Basis of {mu in Im K : mu_o = 0, mu_p + mu_q = 0} as a list of integer vectors (possibly empty)."""
    Kt = sp.Matrix(K)  # 5 x 4, columns span Im K (slopes = K kappa)
    C = sp.Matrix([[1 if a == o else 0 for a in range(5)], [1 if a in (p, q) else 0 for a in range(5)]])
    vecs = [Kt * n for n in (C * Kt).nullspace()]
    vecs = [v for v in vecs if any(x != 0 for x in v)]
    if not vecs: return []
    M = sp.Matrix.hstack(*vecs)
    cols = M.columnspace()
    return [primitive(c) for c in cols]

def orientations(W, p, q, r, s):
    """Achievable (sub_1, quot_1, sub_2, quot_2) with mu_sub < 0 < mu_quot on both pairs for some mu in W."""
    out = []
    if len(W) == 1:
        w = W[0]
        if w[q] == 0 or w[s] == 0: return out
        for f in (1, -1):
            s1, q1 = (p, q) if f * w[p] < 0 else (q, p)
            s2, q2 = (r, s) if f * w[r] < 0 else (s, r)
            out.append((s1, q1, s2, q2))
        return out
    if len(W) == 2:
        # functionals mu_q and mu_s on W in the basis W = [w1, w2]
        fq = (W[0][q], W[1][q]); fs = (W[0][s], W[1][s])
        for sq_sign in (1, -1):
            for ss_sign in (1, -1):
                # is there (x, y) with sign(fq . (x,y)) = sq_sign and sign(fs . (x,y)) = ss_sign ?
                ok = False
                if fq != (0, 0) and fs != (0, 0):
                    det = fq[0] * fs[1] - fq[1] * fs[0]
                    if det != 0: ok = True
                    else:  # proportional: fs = c fq, c = sign-determining scalar
                        c = sp.Rational(fs[0], fq[0]) if fq[0] != 0 else sp.Rational(fs[1], fq[1])
                        ok = (c > 0) == (sq_sign == ss_sign)
                if ok:
                    s1, q1 = (p, q) if sq_sign > 0 else (q, p)
                    s2, q2 = (r, s) if ss_sign > 0 else (s, r)
                    out.append((s1, q1, s2, q2))
        return out
    return out

POLYSTABLE = {  # exact Kaehler classes with all five slopes zero (t4 = 1), checked below
    "BCL fixed datum": [(1, 1, 1, 1)],
    "class 1": [(1, (sp.sqrt(7) - 1) / 3, (sp.sqrt(7) - 1) / 3, 1)],
    "class 2": [(sp.Rational(-28, 5) + (-39 * r_sym + 18 * r_sym**3 + 27 * r_sym**2) / 5,
                 sp.Rational(-2, 5) + (-22 * r_sym + 6 * r_sym**2 + 9 * r_sym**3) / 10, r_sym, 1)],
    "class 4": [(1, 1, sp.Rational(1, 4), 1)],
    "class 5": [(4, 4, 4, 1)],
    "class 6": [(1, 1, 1, 1)],
    "class 7": [(1, (5 + sp.sqrt(105)) / 8, (7 + sp.sqrt(105)) / 4, 1)],
    "class 8": [(1, 1, 1, 1)],
    "class 9": [(2 + sp.sqrt(7), 1, 2 + sp.sqrt(7), 1)],
    "class 10": [(2 + sp.sqrt(7), 1, 2 + sp.sqrt(7), 1)],
}
R_POLY = 9 * r_sym**4 - 36 * r_sym**2 + 14 * r_sym + 16   # class 2: r is its root in (1.42, 1.43)
R_VALUE = sp.CRootOf(R_POLY, 3)

def at_point(expr, t0):
    e = expr.subs(dict(zip(t, t0)))
    if e.has(r_sym):
        e = sp.rem(sp.Poly(sp.expand(e), r_sym), sp.Poly(R_POLY, r_sym)).as_expr()
    return sp.simplify(e)

def main():
    data = [("BCL fixed datum", [tuple(k) for k in BCL])] + [(f"class {i}", [tuple(k) for k in K]) for i, (g, K) in enumerate(ENGINE, 1)]
    results = []; passed = failed = 0
    def check(ok, what):
        nonlocal passed, failed
        if ok: passed += 1
        else: failed += 1; print("FAILED:", what)
    for name, K in data:
        mus = [slope(k) for k in K]
        U = traceless_kernel(K)
        rec = {"name": name, "rank": rank(K), "u": U, "partitions": [], "polystable": [], "single_block": []}
        for t0 in POLYSTABLE.get(name, []):
            vals = [at_point(m, t0) for m in mus]
            check(all(v == 0 for v in vals), f"{name}: slopes vanish at {t0}")
            num = [float(sp.N(sp.sympify(x).subs(r_sym, R_VALUE))) for x in t0]
            check(all(v > 0 for v in num), f"{name}: polystable point in the open Kaehler cone")
            rec["polystable"].append([str(x) for x in t0])
        if name == "class 3":
            rec["note"] = "rank 2: slopes vanish on a curve, two Abelian generators; excluded by C8 (Proposition 16.4)"
        # one rank-two block: the other three rows independent => no single extension at any Kaehler class
        single_open = []
        for p, q in itertools.combinations(range(5), 2):
            others = [K[a] for a in range(5) if a not in (p, q)]
            N = sp.Matrix(others).nullspace()          # kappa with mu_a = 0 for the three other summands
            mu_p_on_N = [sum(K[p][i] * n[i] for i in range(4)) for n in N]
            e = {"pair": [p, q], "rank_of_others": rank(others), "dim_N": len(N),
                 "mu_p_vanishes_on_N": all(x == 0 for x in mu_p_on_N),
                 "ext_dims": {f"{p}{q}": h1(vsub(K[p], K[q])), f"{q}{p}": h1(vsub(K[q], K[p]))},
                 "u_higgsed": any(u[p] != u[q] for u in U), "isomorphic": tuple(K[p]) == tuple(K[q])}
            if e["dim_N"] == 1 or e["mu_p_vanishes_on_N"]: e["verdict"] = "no slope-zero point"
            elif e["isomorphic"]: e["verdict"] = "excluded (isomorphic summands)"
            elif all(v == 0 for v in e["ext_dims"].values()): e["verdict"] = "excluded"
            elif any(v is None for v in e["ext_dims"].values()) and not any(v for v in e["ext_dims"].values() if v): e["verdict"] = "undecided"; single_open.append((p, q))
            else: e["verdict"] = "slope condition satisfiable" + (", Abelian factor kept" if not e["u_higgsed"] else ""); single_open.append((p, q))
            rec["single_block"].append(e)
        rec["single_block_excluded"] = not single_open
        # two rank-two blocks
        realisable = []; undecided = []
        seen = set()
        for o in range(5):
            rest = [a for a in range(5) if a != o]
            for p, q in itertools.combinations(rest, 2):
                r, s = [a for a in rest if a not in (p, q)]
                key = (o, frozenset([frozenset([p, q]), frozenset([r, s])]))
                if key in seen: continue
                seen.add(key)
                W = space_W(K, p, q, o)
                entry = {"pairs": [[p, q], [r, s]], "o": o, "dim_W": len(W), "W": W,
                         "same_determinant": vsub(vsub(K[p], K[r]), vsub(K[s], K[q])) == (0, 0, 0, 0)}
                ors = orientations(W, p, q, r, s)
                entry["orientations"] = []
                verdict = "excluded" if ors else "no slope-zero point"
                for (s1, q1, s2, q2) in ors:
                    h_1 = h1(vsub(K[s1], K[q1])); h_2 = h1(vsub(K[s2], K[q2]))
                    higgs = any(u[s1] != u[q1] or u[s2] != u[q2] for u in U)
                    orec = {"sub": [s1, s2], "quotient": [q1, q2], "ext_dims": [h_1, h_2], "u_higgsed": higgs}
                    entry["orientations"].append(orec)
                    if h_1 and h_2:
                        if higgs and not entry["same_determinant"]:
                            realisable.append((s1, q1, s2, q2, o, h_1, h_2, len(W)))
                            verdict = "realisable"
                        else:
                            verdict = "realisable, Abelian factor kept" if verdict != "realisable" else verdict
                    elif h_1 is None or h_2 is None:
                        if h_1 == 0 or h_2 == 0: pass
                        elif verdict not in ("realisable",): verdict = "undecided"; undecided.append((s1, q1, s2, q2, o))
                entry["verdict"] = verdict
                rec["partitions"].append(entry)
        rec["realisable"] = [list(x) for x in realisable]
        rec["undecided"] = [list(x) for x in undecided]
        rec["two_blocks_excluded"] = not realisable and not undecided and not any(
            e["verdict"].startswith("realisable") for e in rec["partitions"])
        # certificate at the polystable point for the first realisable structure with a line W
        lines = [x for x in realisable if x[7] == 1]
        if lines and name in POLYSTABLE:
            s1, q1, s2, q2, o, h_1, h_2, _ = lines[0]
            t0 = POLYSTABLE[name][0]
            cons = [mus[o], mus[s1] + mus[q1]]
            G = sp.Matrix([[at_point(sp.diff(c, v), t0) for v in t] for c in cons])
            check(G.rank() == 2, f"{name}: S is a smooth surface at J_0")
            T = G.nullspace()
            dq = sp.Matrix([[at_point(sp.diff(mus[q1], v), t0) for v in t]])
            check(any(sp.simplify((dq * n)[0]) != 0 for n in T), f"{name}: mu_q has nonzero differential along S")
            rec["certificate"] = {"sub": [s1, s2], "quotient": [q1, q2], "o": o, "ext_dims": [h_1, h_2], "J0": [str(x) for x in t0]}
        results.append(rec)
    # class 1: the two groups H^1(L_2 L_5^-1), H^1(L_5 L_2^-1) are bounded by 24 and vanish on every member for which
    # multiplication by the quadric is injective on H^2(A, k - 2) -> H^2(A, k); for s_0 that map is a bijection of
    # the monomial bases (a permutation matrix), so the groups vanish on the general tetraquadric.
    K1 = [tuple(k) for k in ENGINE[0][1]]
    for k, degs in ((vsub(K1[1], K1[4]), (1, 1, 0, 0)), (vsub(K1[4], K1[1]), (0, 0, 1, 1))):
        a = hA(k); b = hA(tuple(x - 2 for x in k))
        check(a[1] == 0 and a[2] == 24 and b[2] == 24, f"class 1: Koszul data of {k}: h^1(A,k) = 0, h^2(A,k) = h^2(A,k-2) = 24")
        M = koszul_matrix(k, degs, S0)
        check(M.shape == (24, 24) and all(sum(1 for x in M.row(i) if x != 0) == 1 for i in range(24))
              and all(sum(1 for x in M.col(j) if x != 0) == 1 for j in range(24)), f"class 1: s_0 permutes the 24 monomial classes for {k}")
        check(abs(M.det()) == 1, f"class 1: the Koszul map for {k} is an isomorphism at s_0")
    for rec in results:
        if rec["name"] == "class 1": rec["generic_vanishing"] = {"groups": ["L_2 L_5^-1", "L_5 L_2^-1"], "quadric": S0, "note": "zero on the general tetraquadric"}
    EXPECT = {"BCL fixed datum": "realisable", "class 1": "undecided", "class 2": "realisable", "class 3": "excluded",
              "class 4": "excluded", "class 5": "excluded", "class 6": "realisable", "class 7": "realisable",
              "class 8": "realisable", "class 9": "excluded", "class 10": "realisable"}
    for rec in results:
        status = "realisable" if rec["realisable"] else ("excluded" if rec["two_blocks_excluded"] else "undecided")
        rec["two_blocks"] = status
        check(status == EXPECT[rec["name"]], f"{rec['name']}: two rank-two blocks {status}, expected {EXPECT[rec['name']]}")
        # a single rank-two block never removes the generator u: either the other three rows are independent
        # (no slope-zero point) or u is constant on the pair
        check(all(e["verdict"] == "no slope-zero point" or not e["u_higgsed"] for e in rec["single_block"]),
              f"{rec['name']}: a single extension keeps the Abelian factor")
        if rec["name"] == "class 1":
            # every undecided orientation has one exact nonzero group and one group equal to H^1(L_2 L_5^-1) or its
            # inverse (L_2 and L_3 are isomorphic for this class), which vanish on the general tetraquadric
            K1 = [tuple(k) for k in ENGINE[0][1]]; target = {vsub(K1[1], K1[4]), vsub(K1[4], K1[1])}
            und = [o for e in rec["partitions"] for o in e["orientations"] if None in o["ext_dims"] and 0 not in o["ext_dims"]]
            check(len(und) == 4 and all(vsub(K1[o["sub"][i]], K1[o["quotient"][i]]) in target
                                        for o in und for i in range(2) if o["ext_dims"][i] is None),
                  "class 1: every undecided orientation needs H^1(L_2 L_5^-1) or H^1(L_5 L_2^-1)")
            rec["two_blocks"] = status = "excluded on the general tetraquadric"
        print(rec["name"], "rank", rec["rank"], "u =", rec["u"], "| two blocks:", status,
              "| single block excluded:", rec["single_block_excluded"], "| polystable:", rec["polystable"])
        for e in rec["partitions"]:
            if e["verdict"] != "no slope-zero point":
                print("   ", e["pairs"], "o", e["o"], "dim W", e["dim_W"], "W", e["W"], e["verdict"],
                      [(x["sub"], x["quotient"], x["ext_dims"], x["u_higgsed"]) for x in e["orientations"]])
        if "certificate" in rec: print("    certificate", rec["certificate"])
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(results, open(os.path.join(here, "extension_branch_results.json"), "w"), indent=1, default=str)
    print(f"{passed} checks passed, {failed} failed")
    sys.exit(1 if failed else 0)

if __name__ == "__main__": main()
