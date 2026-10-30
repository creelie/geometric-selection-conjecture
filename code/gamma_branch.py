"""The eleven data and the extension branch on the quotient Y = X_pi / Gamma (Sections 16.2 and 16.4).

Part A checks, for the fixed datum and the ten further sums of Table 7, the properties stated in
Section 16.2: c_1(V) = 0, the index chi(X, V) = -3|G| for the stated group order, c_2(V).J_i <= 24 =
c_2(T_X).J_i, and H^2(X, V) = 0 on every smooth tetraquadric (Koszul bound).

Part B works on the fixed tetraquadric X_pi of Section 13 with the free group Gamma = <a, b>,
a: [x:y] -> [x:-y], b: [x:y] -> [y:x] on every factor.  A line bundle O(k) carries the untwisted lift when
sum(k) is even (the lifts of a and b then commute); a structure twisted by a character chi multiplies every
character of its cohomology by chi.  Characters are ordered (++), (+-), (-+), (--) as values on (a, b).

Cohomology of O_X(k) by characters, exactly: Koszul 0 -> O_A(k-2) -> O_A(k) -> O_X(k) -> 0 on
A = (P^1)^4; every O_A(m) has at most one nonzero group (Kuenneth), with the Cech monomial basis (a acts
by (-1)^(sum of y-exponents), b sends v to m - v with the sign (-1)^q).  When source and target sit in one
degree the map is multiplication by F_theta = F_0 + theta G of Section 13, whose rank on each character is
computed over Q(theta); at theta = pi, which is transcendental, the rank is the same.

A two-block bundle is E = E_1 + E_2 + L_o with 0 -> L_s(i) -> E_i -> L_q(i) -> 0 non-split; its shapes and
orientations (those allowed by the slopes) are read from extension_branch_results.json.  It descends to Y
when the line bundles carry equivariant structures with trivial product and both extension classes are
Gamma-invariant.  A Wilson line rho(a), rho(b) in {1, D}, D = diag(1,1,1,-1,-1), not both 1, gives the
character eps with eps(g) = -1 iff rho(g) = D; Q lies in H^1(V)^eps, u^c in H^1(V)^1, d^c in
H^1(wedge^2 V)^1, H_d in H^1(wedge^2 V)^eps, and H_u is dual to H^2(wedge^2 V)^eps.

Selection rules (Lemma 16.14).  (T) The block torus (t_1, t_2, t_3), scalar on E_1, E_2, L_o, acts on
H^3(det V) by t_1^2 t_2^2 t_3; a cubic coupling of weight vectors vanishes unless the weights add up to
(2, 2, 1).  (F) For a set S of indices containing s(i) whenever it contains q(i), V_S = sum of the L_a with
a in S, completed to the blocks, is a subbundle; with the filtration by the number of factors from V_S,
a coupling of classes of levels l_1, l_2, l_3 vanishes when l_1 + l_2 + l_3 > |S|.  A class of exact level l
and weight w maps to a nonzero class of the graded piece of level l and weight w, so the multiplicity of
(weight, level) is at most the cohomology of that graded piece.  The rank of a Yukawa matrix is at most its
term rank (a maximum matching in the bipartite graph of the allowed pairs, with multiplicities).

The same bounds are also recomputed with the Koszul upper bounds (16.2) in place of the exact dimensions and
over every orientation allowed by the slopes, which makes them valid on every smooth Gamma-invariant
tetraquadric on which Gamma acts freely (Remark 16.18).

Output: gamma_branch_results.json.
"""
import itertools, json, os, sys
from collections import Counter
from itertools import combinations, product
import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from yukawa_charges import BCL, ENGINE, hX_upper

HERE = os.path.dirname(os.path.abspath(__file__))
CHARS = list(product((1, -1), repeat=2))
TRIV = (1, 1)
EPS = [(-1, 1), (1, -1), (-1, -1)]
TH = sp.Symbol('theta')
DET = (2, 2, 1)


def cmul(x, y):
    return (x[0] * y[0], x[1] * y[1])


def add(*ks):
    return tuple(sum(x) for x in zip(*ks))


def neg(k):
    return tuple(-x for x in k)


# ---------------------------------------------------------------- the equation F_theta of Section 13
def deformation_terms():
    terms = {}
    for coefficient, I in zip((1, 2, 3, 5, 7, 11), combinations(range(4), 2)):
        J = [j for j in range(4) if j not in I]
        for bits in product((0, 2), repeat=2):
            ex = [0] * 4
            for i in I:
                ex[i] = 1
            for j, bb in zip(J, bits):
                ex[j] = bb
            terms[tuple(ex)] = terms.get(tuple(ex), 0) + coefficient
    terms[(1, 1, 1, 1)] = 13
    return terms


G_TERMS = deformation_terms()
F0_TERMS = {bits: (1, 2, 5, 2, 1)[sum(b // 2 for b in bits)] for bits in product((0, 2), repeat=4)}


def ambient(k):
    """(degree, basis of y-exponent vectors) of the unique nonzero group of O_A(k), or (None, [])."""
    if -1 in k:
        return None, []
    q = sum(d < 0 for d in k)
    return q, list(product(*(range(d + 1) if d >= 0 else range(d + 1, 0) for d in k)))


def eigenbasis(k, ch):
    q, bs = ambient(k)
    if q is None:
        return []
    look = {v: i for i, v in enumerate(bs)}
    seen, out = set(), []
    for i, v in enumerate(bs):
        if i in seen:
            continue
        j = look[tuple(d - e for d, e in zip(k, v))]
        seen.update([i, j])
        if (-1) ** sum(v) != ch[0]:
            continue
        sb = (-1) ** q
        if i == j:
            if sb == ch[1]:
                out.append({i: 1})
        else:
            out.append({i: 1, j: ch[1] * sb})
    return out


def mult_matrix(k):
    km = tuple(d - 2 for d in k)
    qs, src = ambient(km)
    qt, tgt = ambient(k)
    if qs is None or qt is None or qs != qt:
        return None, src, tgt
    look = {v: i for i, v in enumerate(tgt)}
    M = sp.zeros(len(tgt), len(src))
    for j, v in enumerate(src):
        for terms, c in ((F0_TERMS, 1), (G_TERMS, TH)):
            for ex, a in terms.items():
                w = tuple(x + e for x, e in zip(v, ex))
                if w in look:
                    M[look[w], j] += a * c
    return M, src, tgt


def generic_rank(M):
    if M.rows == 0 or M.cols == 0:
        return 0
    return DomainMatrix.from_Matrix(M).convert_to(sp.QQ.frac_field(TH)).rank()


_CACHE = {}


def coh(k):
    """(dims h^0..h^3, {q: {character: multiplicity}}) of O_{X_pi}(k) with the untwisted structure."""
    k = tuple(k)
    if k in _CACHE:
        return _CACHE[k]
    km = tuple(d - 2 for d in k)
    qs, _ = ambient(km)
    qt, _ = ambient(k)
    out = {q: {ch: 0 for ch in CHARS} for q in range(4)}
    M, src, tgt = mult_matrix(k)
    for ch in CHARS:
        se, te = eigenbasis(km, ch), eigenbasis(k, ch)
        r = 0
        if M is not None and se and te:
            cols = []
            for vec in se:
                col = sp.zeros(len(src), 1)
                for i, c in vec.items():
                    col[i] = c
                image = M * col
                coords = sp.Matrix([image[min(vv)] for vv in te])
                recon = sp.zeros(len(tgt), 1)
                for c, vv in zip(coords, te):
                    for i, a in vv.items():
                        recon[i] += c * a
                assert sp.expand(recon - image) == sp.zeros(len(tgt), 1), ("equivariance", k, ch)
                cols.append(coords)
            r = generic_rank(sp.Matrix.hstack(*cols))
        if qt is not None:
            out[qt][ch] += len(te) - (r if M is not None else 0)
        if qs is not None:
            out[qs - 1][ch] += len(se) - (r if M is not None else 0)
    dims = [sum(out[q].values()) for q in range(4)]
    _CACHE[k] = (dims, out)
    return _CACHE[k]


def twisted(k, chi):
    _, o = coh(k)
    return {q: {cmul(c, chi): o[q][c] for c in CHARS} for q in range(4)}


def index(k):
    return 2 * sum(k[i] * k[j] * k[l] for i, j, l in combinations(range(4), 3)) + 2 * sum(k)


def term_rank(rows, cols, allowed):
    """Maximum matching between row slots and column slots with multiplicities (augmenting paths)."""
    R = [(r, i) for r, m in rows.items() for i in range(m)]
    Cs = [(c, j) for c, m in cols.items() for j in range(m)]
    match = {}

    def aug(u, seen):
        for v in Cs:
            if (u[0], v[0]) in allowed and v not in seen:
                seen.add(v)
                if v not in match or aug(match[v], seen):
                    match[v] = u
                    return True
        return False
    return sum(1 for u in R if aug(u, set()))


R_SYM = sp.Symbol('r')
R_POLY = 9 * R_SYM ** 4 - 36 * R_SYM ** 2 + 14 * R_SYM + 16
R_ROOT = sp.CRootOf(sp.Poly(R_POLY, R_SYM), 3)   # the root 1.42764... in (1.4, 1.5)
_s7, _s105 = sp.sqrt(7), sp.sqrt(105)
SLOPE_ZERO = {
    "fixed": (1, 1, 1, 1), "class 1": (1, (_s7 - 1) / 3, (_s7 - 1) / 3, 1),
    "class 2": ((18 * R_SYM ** 3 + 27 * R_SYM ** 2 - 39 * R_SYM - 28) / 5, (9 * R_SYM ** 3 + 6 * R_SYM ** 2 - 22 * R_SYM - 4) / 10, R_SYM, 1),
    "class 3": (1, sp.Rational(5, 2), 1, 1), "class 4": (1, 1, sp.Rational(1, 4), 1), "class 5": (4, 4, 4, 1),
    "class 6": (1, 1, 1, 1), "class 7": (1, (5 + _s105) / 8, (7 + _s105) / 4, 1), "class 8": (1, 1, 1, 1),
    "class 9": (2 + _s7, 1, 2 + _s7, 1), "class 10": (2 + _s7, 1, 2 + _s7, 1)}

DATA = [("fixed", 4, [tuple(r) for r in BCL])] + [(f"class {i}", g, [tuple(r) for r in K]) for i, (g, K) in enumerate(ENGINE, 1)]
RES = json.load(open(os.path.join(HERE, "extension_branch_results.json")))


def structures(name):
    entry = next(e for e in RES if (e["name"] == "BCL fixed datum" and name == "fixed") or e["name"] == name)
    return [(o["sub"][0], o["quotient"][0], o["sub"][1], o["quotient"][1], p["o"])
            for p in entry["partitions"] for o in p["orientations"]]


def live_structures(name, K):
    """Shapes and orientations allowed by the slopes whose two extension groups are nonzero on X_pi."""
    return [st for st in structures(name)
            if coh(add(K[st[0]], neg(K[st[1]])))[0][1] and coh(add(K[st[2]], neg(K[st[3]])))[0][1]]


def blocks(st):
    s1, q1, s2, q2, o = st
    return {s1: 0, q1: 0, s2: 1, q2: 1, o: 2}


def weight(st, *idx):
    b = blocks(st)
    v = [0, 0, 0]
    for a in idx:
        v[b[a]] += 1
    return tuple(v)


def bounds(K, st, mult1, mult_pair1, hu_pairs, hd_pairs, rules):
    """Term-rank bounds of Y_u and Y_d under the torus rule and the filtrations of `rules` (sets S).

    mult1[a] (dict char -> mult) for Q (index 0) and u^c (index 1) from H^1(L_a); mult_pair1[(a,b)] the
    d^c multiplicity from H^1(L_a L_b); hu_pairs / hd_pairs the pairs whose H^2 / H^1 can carry H_u / H_d."""
    out = {}
    for name, S in rules.items():
        lev = {a: (1 if a in S else 0) for a in range(5)}
        Ld = len(S)
        Q, U, D = Counter(), Counter(), Counter()
        for a in range(5):
            Q[(weight(st, a), lev[a])] += mult1[a][0]
            U[(weight(st, a), lev[a])] += mult1[a][1]
        for (a, b), m in mult_pair1.items():
            D[(weight(st, a, b), lev[a] + lev[b])] += m
        hu = {(weight(st, a, b), -(lev[a] + lev[b])) for (a, b) in hu_pairs}   # the dual piece of wedge^3 V
        hd = {(weight(st, a, b), lev[a] + lev[b]) for (a, b) in hd_pairs}
        au = {(x, y) for x in Q for y in U for (h, c) in hu
              if add(x[0], y[0]) == h and x[1] + y[1] + c <= 0}
        ad = {(x, y) for x in Q for y in D for (h, c) in hd
              if add(x[0], y[0], h) == DET and x[1] + y[1] + c <= Ld}
        Qp = {k: m for k, m in Q.items() if m}
        out[name] = (term_rank(Qp, {k: m for k, m in U.items() if m}, au),
                     term_rank(Qp, {k: m for k, m in D.items() if m}, ad))
    return out


def admissible(st, S):
    s1, q1, s2, q2, o = st
    return (q1 not in S or s1 in S) and (q2 not in S or s2 in S)


def main():
    passed = failed = 0

    def check(ok, what):
        nonlocal passed, failed
        if ok:
            passed += 1
        else:
            failed += 1
            print("FAILED:", what)

    results = {"data": [], "gamma": {}}

    # ---------------------------------------------------------------- part A: the eleven data
    for name, g, K in DATA:
        c1 = [sum(k[i] for k in K) for i in range(4)]
        ind = sum(index(k) for k in K)
        c2 = [-2 * sum(sum(k[x] * k[y] for x, y in combinations([j for j in range(4) if j != i], 2)) for k in K)
              for i in range(4)]
        h2 = sum(hX_upper(k)[2] for k in K)
        check(c1 == [0, 0, 0, 0], f"{name}: c_1(V) = 0")
        check(ind == -3 * g, f"{name}: chi(X, V) = -3|G| with |G| = {g}")
        check(all(c <= 24 for c in c2), f"{name}: c_2(V).J_i <= 24 = c_2(T_X).J_i")
        check(h2 == 0, f"{name}: H^2(X, V) = 0 on every smooth tetraquadric")
        lo = -3 if name == "fixed" else -2
        check(all(lo <= x <= 2 for k in K for x in k), f"{name}: entries in [{lo}, 2]")
        # a Kaehler class with all slopes zero, at an explicit algebraic point (exact arithmetic)
        tt = SLOPE_ZERO[name]
        kap = [sum(tt[x] * tt[y] for x, y in combinations([j for j in range(4) if j != i], 2)) for i in range(4)]
        mu = [sum(k[i] * kap[i] for i in range(4)) for k in K]
        if name == "class 2":
            mu = [sp.rem(sp.expand(m), R_POLY, R_SYM) for m in mu]
            positive = all(sp.sympify(x).subs(R_SYM, R_ROOT).evalf(30) > 0 for x in tt)
        else:
            positive = all(sp.simplify(x) > 0 for x in tt)
        check(all(sp.simplify(m) == 0 for m in mu) and positive, f"{name}: all slopes vanish at t = {tt}")
        results["data"].append(dict(name=name, order=g, c2J=c2, index=ind,
                                    gamma_compatible=all(sum(k) % 2 == 0 for k in K)))
    gam = [name for name, g, K in DATA if all(sum(k) % 2 == 0 for k in K)]
    check(gam == ["fixed", "class 6", "class 8", "class 9", "class 10"],
          "the data whose rows all have even sum are the fixed datum and the classes 6, 8, 9, 10")

    # ---------------------------------------------------------------- part B: the Gamma-branch
    # the class 9 has no two-block bundle (Theorem 16.10(i))
    expect_families = {"fixed": [0, 2, 0, 0, 1], "class 6": [0, 0, 1, 1, 1], "class 8": [0, 0, 1, 1, 1],
                       "class 9": [0, 0, 0, 0, 3], "class 10": [0, 0, 0, 0, 3]}
    expect_live = {"fixed": [(2, 4, 5, 3, 1), (2, 1, 5, 3, 4), (5, 1, 2, 3, 4)],
                   "class 6": [(3, 2, 5, 4, 1), (2, 5, 4, 3, 1), (4, 1, 5, 2, 3), (3, 1, 5, 2, 4), (4, 1, 3, 2, 5)],
                   "class 8": [(3, 2, 4, 5, 1), (3, 2, 5, 4, 1), (4, 1, 5, 2, 3), (5, 1, 4, 2, 3)],
                   "class 9": [],
                   "class 10": [(5, 1, 4, 2, 3), (2, 1, 5, 3, 4)]}
    expect_pieces = {"fixed": {(1,): [0, 8, 0, 0], (4,): [0, 4, 0, 0], (1, 3): [0, 4, 0, 0],
                               (1, 4): [0, 3, 3, 0], (3, 4): [0, 8, 0, 0]},
                     "class 6": {(2,): [0, 4, 0, 0], (3,): [0, 4, 0, 0], (4,): [0, 4, 0, 0], (0, 1): [0, 4, 0, 0],
                                 (0, 4): [0, 4, 0, 0], (1, 3): [0, 4, 0, 0], (1, 4): [0, 3, 3, 0]},
                     "class 8": {(2,): [0, 4, 0, 0], (3,): [0, 4, 0, 0], (4,): [0, 4, 0, 0], (0, 1): [0, 4, 0, 0],
                                 (0, 2): [0, 8, 0, 0], (1, 3): [0, 3, 3, 0], (1, 4): [0, 3, 3, 0], (3, 4): [0, 3, 3, 0]},
                     "class 9": {(4,): [0, 12, 0, 0], (0, 4): [0, 12, 0, 0], (2, 3): [0, 3, 3, 0]},
                     "class 10": {(4,): [0, 12, 0, 0], (0, 4): [0, 8, 0, 0], (1, 3): [0, 4, 0, 0], (1, 4): [0, 3, 3, 0]}}
    expect_blind = {"fixed": (2, 0), "class 6": (1, 3), "class 8": (2, 3), "class 10": (0, 0)}
    for name, g, K in DATA:
        if name not in expect_families:
            continue
        rec = {}
        # every line bundle and every pair with nonzero cohomology on X_pi, exactly
        pieces = {}
        for a in range(5):
            d, _ = coh(K[a])
            if any(d):
                pieces[(a,)] = d
        for a, b in combinations(range(5), 2):
            d, _ = coh(add(K[a], K[b]))
            if any(d):
                pieces[(a, b)] = d
        check(pieces == expect_pieces[name], f"{name}: the pieces with nonzero cohomology on X_pi")
        rec["pieces"] = {"".join(f"L{x + 1}" for x in key): d for key, d in pieces.items()}
        # H^1 of every summand is a multiple of the regular representation (Gamma acts freely)
        fam = []
        for a in range(5):
            _, o = coh(K[a])
            vals = set(o[1].values())
            check(len(vals) == 1 and all(o[q][c] == 0 for q in (0, 2, 3) for c in CHARS),
                  f"{name}: H^*(L_{a + 1}) is H^1, a multiple of the regular representation")
            fam.append(vals.pop())
        check(fam == expect_families[name], f"{name}: families per summand {fam}")
        rec["families"] = fam
        live = live_structures(name, K)
        check(sorted(tuple(x + 1 for x in st) for st in live) == sorted(expect_live[name]),
              f"{name}: shapes with both extension groups nonzero on X_pi")
        rec["structures"] = []
        for st in live:
            s1, q1, s2, q2, o = st
            srec = dict(structure=[x + 1 for x in st],
                        ext=[coh(add(K[s1], neg(K[q1])))[0][1], coh(add(K[s2], neg(K[q2])))[0][1]])
            # ---- descent: equivariant structures (characters chi_a with trivial product) with invariant classes
            worst_u = worst_d = 0
            ncases = 0
            for chis in product(CHARS, repeat=4):
                last = TRIV
                for c in chis:
                    last = cmul(last, c)
                chi = list(chis) + [last]
                inv1 = twisted(add(K[s1], neg(K[q1])), cmul(chi[s1], chi[q1]))[1][TRIV]
                inv2 = twisted(add(K[s2], neg(K[q2])), cmul(chi[s2], chi[q2]))[1][TRIV]
                if not inv1 or not inv2:
                    continue
                for eps in EPS:
                    mult1 = {a: (twisted(K[a], chi[a])[1][eps], twisted(K[a], chi[a])[1][TRIV]) for a in range(5)}
                    mp, hu, hd = {}, [], []
                    for a, b in combinations(range(5), 2):
                        tw = twisted(add(K[a], K[b]), cmul(chi[a], chi[b]))
                        if tw[1][TRIV]:
                            mp[(a, b)] = tw[1][TRIV]
                        if tw[2][eps]:
                            hu.append((a, b))
                        if tw[1][eps]:
                            hd.append((a, b))
                    rules = {"T": set(), "F_s1": {s1}, "F_s2": {s2}}
                    b = bounds(K, st, mult1, mp, hu, hd, rules)
                    ru = min(v[0] for v in b.values())
                    rd = min(v[1] for v in b.values())
                    worst_u, worst_d = max(worst_u, ru), max(worst_d, rd)
                    ncases += 1
            srec["descending_cases"] = ncases
            srec["twisted_bounds"] = [worst_u, worst_d]
            # ---- character-blind bounds: every structure and Wilson line (multiplicity of every character
            # of H^1(L_a) is fam[a]; every pair with h^2 != 0 may carry H_u, every pair with h^1 != 0 H_d)
            mult1 = {a: (fam[a], fam[a]) for a in range(5)}
            mp = {(a, b): coh(add(K[a], K[b]))[0][1] for a, b in combinations(range(5), 2) if coh(add(K[a], K[b]))[0][1]}
            hu = [(a, b) for a, b in combinations(range(5), 2) if coh(add(K[a], K[b]))[0][2]]
            hd = list(mp)
            rules = {"T": set(), "F_s1": {s1}, "F_s2": {s2}}
            b = bounds(K, st, mult1, mp, hu, hd, rules)
            srec["blind_bounds"] = [min(v[0] for v in b.values()), min(v[1] for v in b.values())]
            check(srec["twisted_bounds"][0] <= srec["blind_bounds"][0] and srec["twisted_bounds"][1] <= srec["blind_bounds"][1],
                  f"{name} {srec['structure']}: the bounds for each twist do not exceed the character-blind bound")
            rec["structures"].append(srec)
        if rec["structures"]:
            bu = max(s["blind_bounds"][0] for s in rec["structures"])
            bd = max(s["blind_bounds"][1] for s in rec["structures"])
            check((bu, bd) == expect_blind[name], f"{name}: rank Y_u <= {bu}, rank Y_d <= {bd} on every two-block bundle")
            check(bu <= 2, f"{name}: rank Y_u <= 2, so C8 fails")
        results["gamma"][name] = rec

    # ---------------------------------------------------------------- the same bounds from the Koszul upper bounds
    for name, g, K in DATA:
        if name not in expect_blind:
            continue
        up = lambda *ks: hX_upper(add(*ks))
        check(all(hX_upper(K[a])[0] == hX_upper(K[a])[2] == hX_upper(K[a])[3] == 0 for a in range(5)),
              f"{name}: every summand has only H^1 on every smooth tetraquadric")
        fam = [-index(K[a]) // 4 for a in range(5)]
        mult1 = {a: (fam[a], fam[a]) for a in range(5)}
        mp = {(a, b): up(K[a], K[b])[1] for a, b in combinations(range(5), 2) if up(K[a], K[b])[1]}
        hu = [(a, b) for a, b in combinations(range(5), 2) if up(K[a], K[b])[2]]
        worst = [0, 0]
        for st in structures(name):
            s1, q1, s2, q2, o = st
            if not up(K[s1], neg(K[q1]))[1] or not up(K[s2], neg(K[q2]))[1]:
                continue
            b = bounds(K, st, mult1, mp, hu, list(mp), {"T": set(), "F_s1": {s1}, "F_s2": {s2}})
            worst = [max(worst[0], min(v[0] for v in b.values())), max(worst[1], min(v[1] for v in b.values()))]
        check(tuple(worst) == expect_blind[name], f"{name}: the bounds {expect_blind[name]} from the Koszul upper bounds alone")
        results["gamma"][name]["koszul_bounds"] = worst

    # ---------------------------------------------------------------- the hand proofs of Theorem 16.15
    def hbit(K, *idx):
        return coh(add(*[K[a] for a in idx]))[0]

    K = DATA[0][2]
    # fixed datum: S = {2, 4, 5} is admissible in all three shapes, L_1, L_3 and every pair meeting {1, 3} are
    # acyclic, so H^*(V) = H^*(V_S), H^*(wedge^2 V) = H^*(wedge^2 V_S) and Y_d factors through wedge^5 of rank 3
    S = {1, 3, 4}
    check(all(admissible(st, S) for st in live_structures("fixed", K)), "fixed: V_S, S = {2,4,5}, is a subbundle in every shape")
    check(all(not any(hbit(K, a)) for a in (0, 2)) and all(not any(hbit(K, a, b)) for a, b in combinations(range(5), 2)
                                                           if {a, b} & {0, 2}),
          "fixed: L_1, L_3 and every pair meeting them are acyclic")
    check(all({st[0], st[2]} == {1, 4} for st in live_structures("fixed", K)), "fixed: L_2 and L_5 are the two sub-line-bundles")
    # spectrum: in S1 = (2,1,5,3,4) and S2 = (5,1,2,3,4) the bundle V_S is split, so the equivariant cohomology of V and
    # wedge^2 V is that of the split sum; the invariant parts of the extension groups are nonzero (descent)
    for st, dims in (((1, 0, 4, 2, 3), (12, 4, 3, 1)), ((4, 0, 1, 2, 3), (12, 20, 3, 5))):
        s1, q1, s2, q2, o = st
        e1, e2 = coh(add(K[s1], neg(K[q1]))), coh(add(K[s2], neg(K[q2])))
        got = (e1[0][1], e2[0][1], e1[1][1][TRIV], e2[1][1][TRIV])
        check(got == dims, f"fixed {[x + 1 for x in st]}: extension groups {got[:2]}, invariant parts {got[2:]}")
        check(q1 in (0, 2) and q2 in (0, 2) and o == 3, f"fixed {[x + 1 for x in st]}: V_S = L_2 + L_5 + L_4 is split")
    # the third shape (2,4,5,3,1) contains the extension of L_4 by L_2, and H^1(L_4 L_5) -> H^2(L_2 L_5) may be nonzero
    check(hbit(K, 3, 4)[1] and hbit(K, 1, 4)[2], "fixed (2,4,5,3,1): a connecting map H^1(L4L5) -> H^2(L2L5) is possible")
    # class 10: S = {1, 2, 4, 5} admissible, L_3 and every pair with L_3 acyclic; up: H_u from L_2 L_5 only, L_2, L_5 in
    # distinct blocks and every 10 from L_5
    K = DATA[10][2]
    check(all(admissible(st, {0, 1, 3, 4}) for st in live_structures("class 10", K)), "class 10: V_S, S = {1,2,4,5}, is a subbundle")
    check(not any(hbit(K, 2)) and all(not any(hbit(K, 2, a)) for a in (0, 1, 3, 4)), "class 10: L_3 and its pairs are acyclic")
    check(all(blocks(st)[1] != blocks(st)[4] for st in live_structures("class 10", K)), "class 10: L_2, L_5 in distinct blocks")
    # class 8: the family of L_3 never meets H_u: L_3 is a sub-line-bundle (filtration by L_3) or the lone summand
    K = DATA[8][2]
    for st in live_structures("class 8", K):
        if st[4] == 2:
            check(all(weight(st, 2)[2] == 1 and weight(st, a, b)[2] == 0 for a, b in [(1, 3), (1, 4), (3, 4)]),
                  f"class 8 {[x + 1 for x in st]}: L_3 is alone and no H_u piece contains it")
        else:
            check(st[0] == 2 and all(2 not in (a, b) for a, b in [(1, 3), (1, 4), (3, 4)]),
                  f"class 8 {[x + 1 for x in st]}: L_3 is a sub-line-bundle and every H_u dual piece contains it")
    print("Gamma-compatible data:", gam)
    for name, rec in results["gamma"].items():
        print(name, "families", rec["families"], "pieces", rec["pieces"])
        for s in rec["structures"]:
            print("   ", s)
    json.dump(results, open(os.path.join(HERE, "gamma_branch_results.json"), "w"), indent=1, default=str)
    print(f"{passed} checks passed, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
