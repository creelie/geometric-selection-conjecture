"""Perturbative U(1)^5 selection rules for SU(5) line-bundle sums on a tetraquadric X in (P^1)^4.

For each class K (five line bundles k_a in Z^4, rows), every smooth X in the class:
  10's   from H^1(L_a)            charge e_a
  5bar's from H^1(L_a L_b)        charge e_a+e_b          (families or H_d)
  5's    from H^2(L_a L_b)        charge -e_a-e_b         (H_u)
  singlets from H^1(L_a L_b^-1)   charge e_a-e_b
Upper bounds valid for every smooth member: h^i(O_X(k)) <= h^i(O_A(k)) + h^{i+1}(O_A(k-2)), Koszul.
A coupling 10_a 5bar_bc 5bar_de with singlet insertions is allowed only if
  (1,1,1,1,1) - (e_a+e_b+e_c+e_d+e_e) is a nonnegative integer combination of available singlet charges,
i.e. a flow from the positive to the negative entries along available singlet edges p->q (charge e_p-e_q).
Output: for each class, whether some choice of (10, 5bar family, 5bar Higgs) admits a down coupling.
"""
import itertools, json, sys
from collections import deque

def hP1(n):
    if n >= 0: return (n + 1, 0)
    if n <= -2: return (0, -n - 1)
    return (0, 0)

def hA(k):
    out = [0] * 5
    for degs in itertools.product(range(2), repeat=4):
        v = 1
        for n, d in zip(k, degs):
            v *= hP1(n)[d]
        out[sum(degs)] += v
    return out

def hX_upper(k):
    a = hA(k); b = hA(tuple(x - 2 for x in k))
    return [a[i] + (b[i + 1] if i + 1 < 5 else 0) for i in range(4)]

def hX_koszul_exact_if_degenerate(k):
    """Exact value when the Koszul connecting map is forced to be zero or an iso by degree reasons: we only use upper bounds."""
    return hX_upper(k)

def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))

def feasible(v, edges):
    """Is v (sum zero) = sum n_pq (e_p - e_q), n_pq >= 0, over edges (p,q)? Max flow with infinite inner capacities."""
    n = len(v)
    S, T = n, n + 1
    cap = {}
    def addc(u, w, c):
        cap[(u, w)] = cap.get((u, w), 0) + c
        cap.setdefault((w, u), 0)
    supply = 0
    for p in range(n):
        if v[p] > 0: addc(S, p, v[p]); supply += v[p]
        if v[p] < 0: addc(p, T, -v[p])
    for (p, q) in edges: addc(p, q, 10 ** 6)
    flow = 0
    while True:
        par = {S: None}; dq = deque([S])
        while dq and T not in par:
            u = dq.popleft()
            for w in range(n + 2):
                if (u, w) in cap and cap[(u, w)] > 0 and w not in par:
                    par[w] = u; dq.append(w)
        if T not in par: break
        path = []; w = T
        while par[w] is not None: path.append((par[w], w)); w = par[w]
        f = min(cap[e] for e in path)
        for (u, w) in path: cap[(u, w)] -= f; cap[(w, u)] += f
        flow += f
    return flow == supply

def analyse(K, exact_h=None):
    K = [tuple(k) for k in K]
    h = exact_h if exact_h else hX_upper
    tens = [a for a in range(5) if h(K[a])[1] > 0]
    fivebars = [(a, b) for a in range(5) for b in range(a + 1, 5) if h(add(K[a], K[b]))[1] > 0]
    fives = [(a, b) for a in range(5) for b in range(a + 1, 5) if h(add(K[a], K[b]))[2] > 0]
    singlets = [(p, q) for p in range(5) for q in range(5) if p != q and h(sub(K[p], K[q]))[1] > 0]
    one = (1,) * 5
    def e(*idx):
        v = [0] * 5
        for i in idx: v[i] += 1
        return tuple(v)
    down_allowed = []
    for a in tens:
        for (b, c) in fivebars:
            for (d, f) in fivebars:
                Q = e(a, b, c, d, f)
                v = sub(one, Q)
                if feasible(v, singlets):
                    down_allowed.append((a, (b, c), (d, f)))
    up_allowed = []
    for a in tens:
        for b in tens:
            if b < a: continue
            for (d, f) in fives:
                v = tuple(-x for x in sub(e(a, b), e(d, f)))  # Q = e_a+e_b-e_d-e_f ; need Q + sum = 0
                if feasible(v, singlets):
                    up_allowed.append((a, b, (d, f)))
    return dict(tens=tens, fivebars=fivebars, fives=fives, singlets=singlets,
                down_allowed=down_allowed, up_allowed=up_allowed)

BCL = [(-1, 0, 0, 1), (-1, -3, 2, 2), (0, 1, -1, 0), (1, 1, -1, -1), (1, 1, 0, -2)]

ENGINE = [[2, [[-2, 1, 1, -1], [-1, -1, 1, 1], [-1, -1, 1, 1], [2, 0, -2, 1], [2, 1, -1, -2]]], [4, [[-2, 1, 0, -1], [-1, -1, 1, 1], [-1, -1, 1, 1], [2, 0, -2, 1], [2, 1, 0, -2]]], [2, [[-1, 0, 0, 1], [-1, 0, 0, 1], [-1, 0, 0, 1], [1, -2, 1, -1], [2, 2, -1, -2]]], [2, [[-1, -1, 1, 0], [-1, 0, 0, 1], [-1, 0, 0, 1], [1, -1, 0, 0], [2, 2, -1, -2]]], [2, [[-1, 0, 1, 0], [-1, 0, 1, 0], [0, -2, 0, 1], [1, 1, -2, 0], [1, 1, 0, -1]]], [4, [[-1, 0, 1, 0], [-1, 1, 0, 0], [0, -2, 1, 1], [1, 0, -2, 1], [1, 1, 0, -2]]], [4, [[-1, 0, 0, 1], [-1, 0, 0, 1], [0, -2, 1, 1], [1, 1, -2, -1], [1, 1, 1, -2]]], [4, [[-1, -1, 1, 1], [-1, 1, 0, 0], [0, -2, 1, 1], [1, 1, -2, 0], [1, 1, 0, -2]]], [4, [[-1, 0, 1, 0], [-1, 1, -2, 0], [0, -1, 0, 1], [0, -1, 0, 1], [2, 1, 1, -2]]], [4, [[-1, 0, 1, 0], [-1, 1, -2, 0], [0, -1, 0, 1], [1, -1, -1, 1], [1, 1, 2, -2]]]]

if __name__ == "__main__" and False:
    print("BCL fixed datum, upper-bound (every smooth member) singlet set:")
    r = analyse(BCL); print(json.dumps(r))
    for g, K in ENGINE:
        r = analyse(K)
        print(g, K, "down:", len(r["down_allowed"]), "up:", len(r["up_allowed"]), "singlet edges:", len(r["singlets"]))

# ---- certificates -------------------------------------------------------------------------------
def closed_sets_violation(v, edges):
    """A set B with no outgoing edge and sum_B v > 0 (then v is not a nonnegative combination), or None."""
    n = len(v)
    for r in range(1, n):
        for B in itertools.combinations(range(n), r):
            Bs = set(B)
            if any(p in Bs and q not in Bs for (p, q) in edges): continue
            if sum(v[p] for p in B) > 0: return list(B)
    return None

def decomposition(v, edges):
    """Explicit n_pq >= 0 with v = sum n_pq (e_p - e_q), by greedy augmenting paths; None if infeasible."""
    n = len(v); v = list(v); used = {}
    adj = {p: [q for (pp, q) in edges if pp == p] for p in range(n)}
    while any(x > 0 for x in v):
        s = next(p for p in range(n) if v[p] > 0)
        # BFS to a sink
        par = {s: None}; dq = deque([s]); t = None
        while dq:
            u = dq.popleft()
            if v[u] < 0: t = u; break
            for w in adj[u]:
                if w not in par: par[w] = u; dq.append(w)
        if t is None: return None
        w = t
        while par[w] is not None:
            used[(par[w], w)] = used.get((par[w], w), 0) + 1; w = par[w]
        v[s] -= 1; v[t] += 1
    return used

def report(K, name):
    r = analyse(K)
    one = (1,) * 5
    def e(*idx):
        v = [0] * 5
        for i in idx: v[i] += 1
        return tuple(v)
    out = {"name": name, "K": [list(k) for k in K], "tens": r["tens"], "fivebars": r["fivebars"], "fives": r["fives"],
           "singlets": r["singlets"], "down_allowed": len(r["down_allowed"]), "up_allowed": len(r["up_allowed"])}
    if r["down_allowed"]:
        a, (b, c), (d, f) = r["down_allowed"][0]
        out["down_example"] = {"ten": a, "fivebar": [b, c], "higgs": [d, f], "singlets": {f"{p}->{q}": m for (p, q), m in decomposition(sub(one, e(a, b, c, d, f)), r["singlets"]).items()}}
    else:
        certs = set()
        for a in r["tens"]:
            for (b, c) in r["fivebars"]:
                for (d, f) in r["fivebars"]:
                    B = closed_sets_violation(sub(one, e(a, b, c, d, f)), r["singlets"]); certs.add(tuple(B) if B else ("??",))
        out["down_certificates"] = sorted(certs)
    if r["up_allowed"]:
        a, b, (d, f) = r["up_allowed"][0]
        v = sub(e(d, f), e(a, b))
        out["up_example"] = {"tens": [a, b], "higgs": [d, f], "singlets": {f"{p}->{q}": m for (p, q), m in decomposition(v, r["singlets"]).items()}}
    else:
        certs = set()
        for a in r["tens"]:
            for b in r["tens"]:
                if b < a: continue
                for (d, f) in r["fives"]:
                    B = closed_sets_violation(sub(e(d, f), e(a, b)), r["singlets"]); certs.add(tuple(B) if B else ("??",))
        out["up_certificates"] = sorted(certs) if r["fives"] else ["no 5 at all"]
    return out

EXPECTED = {  # (down allowed, up allowed) as in Table 6 and Proposition 16.4 of the paper
    "BCL fixed datum": (False, True),
    "class 1": (True, True), "class 2": (True, True), "class 3": (True, False), "class 4": (True, False),
    "class 5": (True, False), "class 6": (True, True), "class 7": (True, True), "class 8": (True, True),
    "class 9": (False, False), "class 10": (False, False),
}

def rank(K):
    """Rank of the integer matrix with rows K over the rationals."""
    from fractions import Fraction
    M = [[Fraction(x) for x in row] for row in K]
    r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r

def main():
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    data = [("BCL fixed datum", BCL)] + [(f"class {i}", K) for i, (g, K) in enumerate(ENGINE, 1)]
    results = [report(K, name if name == "BCL fixed datum" else f"{name} (group order {ENGINE[int(name.split()[1]) - 1][0]})") for name, K in data]
    json.dump(results, open(os.path.join(here, "yukawa_charges_results.json"), "w"), indent=1)
    passed = failed = 0
    def check(ok, what):
        nonlocal passed, failed
        if ok: passed += 1
        else:
            failed += 1; print("FAILED:", what)
    for (name, K), o in zip(data, results):
        print(o["name"], "| 10s", o["tens"], "| down", o["down_allowed"], o.get("down_certificates", o.get("down_example")), "| up", o["up_allowed"], o.get("up_certificates", o.get("up_example")))
        down, up = EXPECTED[name]
        check((o["down_allowed"] > 0) == down, f"{name}: down-type coupling")
        check((o["up_allowed"] > 0) == up, f"{name}: up-type coupling")
        check(sum(sum(k) for k in K) == 0 and all(sum(k[i] for k in K) == 0 for i in range(4)), f"{name}: c_1(V) = 0")
        check(rank(K) <= 3, f"{name}: rank K <= h^{{1,1}} - 1 = 3 (Theorem 16.1)")
    check(rank(BCL) == 3, "fixed datum: rank K = 3")
    print(f"{passed} checks passed, {failed} failed")
    sys.exit(1 if failed else 0)
if __name__ == "__main__": main()
