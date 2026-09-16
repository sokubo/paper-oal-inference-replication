#!/usr/bin/env python3
"""example_A1.py -- Example A1 of Web Appendix B (paper 3, review round 2, 2026-09-16).

A bounded dependent design in which the two reduced logistic fits (on x = (1, X) and on
x~ = (1, X, T)) have exactly the same pseudo-true propensity and both reduced AIPW targets
are unbiased, yet the corrected influence functions phi*_{P,{X}} and phi*_{P,{X,T}} of Lemma A4
differ: stability of the reduced propensity is not sufficient for (A3''). The mechanism is
E[T | X] = X^2 - 5/2 != 0, which makes the T-block of D_{b_a} (the derivative of the AIPW map
with respect to the added outcome coefficient) nonzero. Replacing T by eta, for which
E[eta | X, I] = 0, restores equality, as Lemma A5 predicts.

All population quantities are computed with exact rational arithmetic (fractions.Fraction);
no data, no simulation, no dependencies beyond the standard library.  The design was
communicated to us by a reviewer; this script is our own implementation of the calculation.

Design.  X uniform on {-2,-1,1,2}; eta independent Rademacher; T = X^2 - 5/2 + eta;
  I in {-1,1} with P(I=1|X=1) = 799/2250, P(I=1|X=2) = 13/480, P(I=1|X=-x) = 1 - P(I=1|X=x),
  independent of eta given X;  e = expit(log(4)(X+I));  Y(a) = a + X + beta_T T + eps, eps ~ N(0,1).
  Reduced pseudo-true propensity of both models: e* = expit(log(2) X).  Influence functions are
  evaluated at beta_T = 0 (the limit along the sequence beta_T = 2 kappa_n -> 0); the AIPW bias is
  linear in beta_T and is checked at beta_T = 1.

Usage:  python3 example_A1.py
"""
from fractions import Fraction as F
from itertools import product


def solve(a, b):
    n = len(b)
    m = [[F(v) for v in row] + [F(r)] for row, r in zip(a, b)]
    for j in range(n):
        ix = next(i for i in range(j, n) if m[i][j])
        m[j], m[ix] = m[ix], m[j]
        pv = m[j][j]
        m[j] = [v / pv for v in m[j]]
        for i in range(n):
            if i != j:
                c = m[i][j]
                m[i] = [u - c * v for u, v in zip(m[i], m[j])]
    return [r[-1] for r in m]


def dot(a, b):
    return sum((u * v for u, v in zip(a, b)), F(0))


def support(t_mode):
    """Support points with probabilities, true and reduced pseudo-true propensities."""
    rows = []
    for x, ins, eta in product((-2, -1, 1, 2), (-1, 1), (-1, 1)):
        p_plus = F(799, 2250) if abs(x) == 1 else F(13, 480)
        if x < 0:
            p_plus = 1 - p_plus
        prob = F(1, 8) * (p_plus if ins == 1 else 1 - p_plus)
        odds = F(4) ** (x + ins)
        e = odds / (1 + odds)                       # expit(log(4)(X+I))
        q = F(2) ** x / (1 + F(2) ** x)             # expit(log(2) X)
        t = (F(x * x) - F(5, 2) + eta) if t_mode == "T" else F(eta)
        rows.append(dict(x=F(x), i=F(ins), t=t, p=prob, e=e, q=q))
    assert sum(r["p"] for r in rows) == 1
    return rows


def corrected_if(rows, augmented, beta_T=F(0)):
    """Corrected influence function phi*_{P,S} of Lemma A4 for S = {X} or {X,T}, at beta_T.

    Sign convention of the paper: J = -E[dU/dtheta'], IF_theta = J^{-1} U, correction +D' IF_theta.
    Since the reduced outcome model is correct at beta_T = 0 (and for the augmented set at any
    beta_T), D_gamma = 0 there; the outcome blocks D_{b_1}, D_{b_0} are computed explicitly and
    the normal error is integrated analytically (E eps = 0, E eps^2 = 1).
    Returns (h1, h0, D-blocks, variance): phi* = h_A(z) * eps on {A = a}, so Var phi* =
    E[e h1^2 + (1-e) h0^2] (the mean-difference part vanishes because the outcome model is exact).
    """
    design = [[F(1), r["x"]] + ([r["t"]] if augmented else []) for r in rows]
    d = len(design[0])
    scores = [sum(r["p"] * xx[j] * (r["e"] - r["q"]) for r, xx in zip(rows, design)) for j in range(d)]
    assert all(v == 0 for v in scores), "reduced logistic score must vanish at e* for both models"
    hs, Ds = [], []
    for arm in (1, 0):
        pa = [r["e"] if arm else 1 - r["e"] for r in rows]
        qa = [r["q"] if arm else 1 - r["q"] for r in rows]
        J = [[sum(r["p"] * ea * xx[j] * xx[k] for r, ea, xx in zip(rows, pa, design)) for k in range(d)] for j in range(d)]
        D = [(1 if arm else -1) * sum(r["p"] * xx[j] * (1 - ea / qi) for r, ea, qi, xx in zip(rows, pa, qa, design)) for j in range(d)]
        Ds.append(D)
        coef = solve(J, D)                          # J^{-1} D
        hs.append([(F(1) if arm else F(-1)) / qi + dot(xx, coef) for qi, xx in zip(qa, design)])
    var = sum(r["p"] * (r["e"] * u * u + (1 - r["e"]) * v * v) for r, u, v in zip(rows, hs[0], hs[1]))
    return hs, Ds, var


def aipw_bias(rows, augmented, beta_T):
    """Population AIPW bias tau(e*_S, m*_S) - tau with the arm-specific projections."""
    design = [[F(1), r["x"]] + ([r["t"]] if augmented else []) for r in rows]
    d = len(design[0])
    proj = []
    for arm in (1, 0):
        ea = [r["e"] if arm else 1 - r["e"] for r in rows]
        mean = [F(arm) + r["x"] + beta_T * r["t"] for r in rows]
        J = [[sum(r["p"] * e_i * xx[j] * xx[k] for r, e_i, xx in zip(rows, ea, design)) for k in range(d)] for j in range(d)]
        rhs = [sum(r["p"] * e_i * yy * xx[j] for r, e_i, yy, xx in zip(rows, ea, mean, design)) for j in range(d)]
        c = solve(J, rhs)
        proj.append([dot(xx, c) for xx in design])
    m1, m0 = proj
    return sum(r["p"] * (u - v + r["e"] * (1 + r["x"] + beta_T * r["t"] - u) / r["q"]
                         - (1 - r["e"]) * (r["x"] + beta_T * r["t"] - v) / (1 - r["q"]))
               for r, u, v in zip(rows, m1, m0)) - 1


def l2_distance_sq(rows, a, b):
    return sum(r["p"] * (r["e"] * (u - v) ** 2 + (1 - r["e"]) * (w - z) ** 2)
               for r, u, v, w, z in zip(rows, a[0], b[0], a[1], b[1]))


if __name__ == "__main__":
    for mode, label in (("T", "Example A1: T = X^2 - 5/2 + eta  (E[T|X] != 0)"),
                        ("eta", "Modified design: T replaced by eta  (E[eta|X,I] = 0)")):
        rows = support(mode)
        print("=" * 78)
        print(label)
        for xv in (-2, -1, 1, 2):
            num = sum(r["p"] * r["t"] for r in rows if r["x"] == xv)
            den = sum(r["p"] for r in rows if r["x"] == xv)
            print(f"  E[T | X={xv:2d}] = {num / den}", end="")
        print()
        print("  AIPW bias at beta_T = 1, base / augmented:", aipw_bias(rows, False, F(1)), aipw_bias(rows, True, F(1)))
        hb, Db, vb = corrected_if(rows, False)
        ha, Da, va = corrected_if(rows, True)
        print(f"  Var phi*_{{X}}   = {vb} = {float(vb):.10f}")
        print(f"  Var phi*_{{X,T}} = {va} = {float(va):.10f}")
        print("  D_{b_1} augmented (intercept, X, T blocks):", [float(v) for v in Da[0]])
        print("  D_{b_0} augmented (intercept, X, T blocks):", [float(v) for v in Da[1]])
        d2 = l2_distance_sq(rows, hb, ha)
        print(f"  ||phi*_{{X,T}} - phi*_{{X}}||_2^2 = {d2} = {float(d2):.6e};  ||.||_2 = {float(d2) ** 0.5:.6f}")
