#!/usr/bin/env python3
"""Solve the public P2 v3 MILP; all problem data comes from the public JSON."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


ROOT = Path(__file__).resolve().parent
PROBLEM = ROOT / "shared" / "problema_p2_public.json"


class Model:
    def __init__(self):
        self.names, self.lb, self.ub, self.integer = [], [], [], []
        self.rows, self.row_lb, self.row_ub = [], [], []

    def var(self, name, lo=-math.inf, hi=math.inf, binary=False):
        k = len(self.names)
        self.names.append(name)
        self.lb.append(0.0 if binary else lo)
        self.ub.append(1.0 if binary else hi)
        self.integer.append(int(binary))
        return k

    def row(self, terms, lo=-math.inf, hi=math.inf):
        self.rows.append({k: v for k, v in terms.items() if v})
        self.row_lb.append(lo)
        self.row_ub.append(hi)

    def eq(self, terms, rhs):
        self.row(terms, rhs, rhs)

    def solve(self, objective):
        ii, jj, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, v in row.items():
                ii.append(i)
                jj.append(j)
                vv.append(v)
        A = coo_matrix((vv, (ii, jj)), shape=(len(self.rows), len(self.names))).tocsr()
        c = np.zeros(len(self.names))
        for k, v in objective.items():
            c[k] = v
        return milp(c, integrality=np.array(self.integer),
                    bounds=Bounds(self.lb, self.ub),
                    constraints=LinearConstraint(A, self.row_lb, self.row_ub),
                    options={"mip_rel_gap": 1e-8})


def pwl(model, label, x, spec):
    """Convex combination with a selected adjacent segment (no SOS2 API needed)."""
    xs, ys = spec["nodes"], spec["values"]
    assert len(xs) == len(ys) and len(xs) >= 2
    assert all(a < b for a, b in zip(xs, xs[1:]))
    lam = [model.var(f"{label}_lambda_{j}", 0, 1) for j in range(len(xs))]
    seg = [model.var(f"{label}_segment_{j}", binary=True) for j in range(len(xs)-1)]
    y = model.var(label)
    model.eq({j: 1 for j in lam}, 1)
    model.eq({j: 1 for j in seg}, 1)
    model.eq({x: 1, **{j: -v for j, v in zip(lam, xs)}}, 0)
    model.eq({y: 1, **{j: -v for j, v in zip(lam, ys)}}, 0)
    for j, lj in enumerate(lam):
        terms = {lj: 1}
        if j > 0:
            terms[seg[j-1]] = -1
        if j < len(seg):
            terms[seg[j]] = -1
        model.row(terms, hi=0)
    return y


def clip(model, well, q, d, r, qlo, qhi, dlo, dhi):
    """Exact union: d<=0,r=q; 0<=d<=q,r=q-d; d>=q,r=0.

    M values are valid over the other two regimes. The supplied d bounds are
    used only here; d and r receive no auxiliary bound constraints.
    """
    a, b, c = [model.var(f"{well}_regime_{s}", binary=True)
               for s in ("nonpositive", "interior", "saturated")]
    model.eq({a: 1, b: 1, c: 1}, 1)
    # Each inequality is relaxed by a worst-case value over the public bounds.
    m_da = max(0, dhi)
    m_db = max(0, -dlo)
    m_dqb = max(0, dhi-qlo)
    m_dqc = max(0, qhi-dlo)
    model.row({d: 1, a: m_da}, hi=m_da)                # a => d<=0
    model.row({d: 1, b: -m_db}, lo=-m_db)             # b => d>=0
    model.row({d: 1, q: -1, b: m_dqb}, hi=m_dqb)     # b => d<=q
    model.row({d: 1, q: -1, c: -m_dqc}, lo=-m_dqc)   # c => d>=q
    # r-q belongs to [-qhi,0] in the union; r-q+d belongs to
    # [dlo,max(0,dhi-qlo)]; r belongs to [0,qhi].
    model.row({r: 1, q: -1, a: qhi}, hi=qhi)
    model.row({r: 1, q: -1, a: -qhi}, lo=-qhi)
    m_mid_pos = max(0, dhi-qlo)
    m_mid_neg = max(0, -dlo)
    model.row({r: 1, q: -1, d: 1, b: m_mid_pos}, hi=m_mid_pos)
    model.row({r: 1, q: -1, d: 1, b: -m_mid_neg}, lo=-m_mid_neg)
    model.row({r: 1, c: qhi}, hi=qhi)
    model.row({r: 1, c: -qhi}, lo=-qhi)
    return a, b, c


def build(p, epsilon):
    m = Model()
    wells = p["wells"]
    domains = p["decision_variables"]["domains"]
    fixed = p["decision_variables"]["fixed"]
    q = {w: m.var(f"q_{w}", *(fixed[w], fixed[w]) if w in fixed else domains[w])
         for w in wells}
    r = {w: m.var(f"r_{w}") for w in wells}
    d = {w: m.var(f"deficit_{w}") for w in wells}
    Q = m.var("Q")
    own = {w: m.var(f"own_{w}") for w in wells}
    real = {w: m.var(f"real_{w}") for w in wells}
    R = m.var("R")
    Z = m.var("Z")
    FOPT = m.var("FOPT")
    FWIT = m.var("FWIT")
    m.row({v: 1 for v in q.values()}, hi=p["decision_variables"]["capacity_requested"])
    m.eq({Q: 20, **{v: -1 for v in q.values()}}, -600)
    for w in wells:
        m.eq({own[w]: 20, q[w]: -1}, -80)
        m.eq({real[w]: 20, r[w]: -1}, -80)
    m.eq({R: 1, **{v: -1 for v in r.values()}}, 0)
    m.eq({Z: 20, R: -1}, -600)
    pw = p["pwl_functions"]
    own2 = pwl(m, "own2_INJECT3", own["INJECT3"], pw["own2_INJECT3"])
    Q2 = pwl(m, "Q2", Q, pw["Q2"])
    Z2 = pwl(m, "Z2", Z, pw["Z2"])
    regimes = {}
    for w in wells:
        coeff = p["layer1_deficit_models"][w]
        terms = {d[w]: 1}
        for key, var in (("own", own[w]), ("Q", Q), ("own2", own2), ("Q2", Q2)):
            if key in coeff:
                terms[var] = -coeff[key]
        m.eq(terms, coeff.get("intercept", 0))
        qlo, qhi = (fixed[w], fixed[w]) if w in fixed else domains[w]
        dlo, dhi = p["deficit_bounds"][w]
        regimes[w] = clip(m, w, q[w], d[w], r[w], qlo, qhi, dlo, dhi)
    fc = p["layer2"]["FOPT"]["coefficients"]
    fterms = {FOPT: 1, Z2: -fc["Z2"]}
    for i, w in enumerate(wells, 1):
        fterms[real[w]] = -fc[f"r{i}"]
    m.eq(fterms, fc["intercept"])
    wc = p["layer2"]["FWIT"]["coefficients"]
    m.eq({FWIT: 1, Z: -wc["Z"]}, wc["intercept"])
    m.row({FWIT: 1}, hi=epsilon)
    return m, {"q": q, "r": r, "d": d, "Q": Q, "own": own,
               "real": real, "R": R, "Z": Z, "FOPT": FOPT,
               "FWIT": FWIT, "own2": own2, "Q2": Q2, "Z2": Z2,
               "regimes": regimes}


def interpolate(x, spec):
    xs, ys = spec["nodes"], spec["values"]
    for j in range(len(xs)-1):
        if x <= xs[j+1] + 1e-8:
            t = (x-xs[j])/(xs[j+1]-xs[j])
            return ys[j] + t*(ys[j+1]-ys[j])
    raise ValueError(f"PWL argument out of range: {x}")


def audit(p, v):
    wells = p["wells"]
    q, r = v["requested"], v["realized"]
    Q = (sum(q.values())-600)/20
    for w in wells:
        c = p["layer1_deficit_models"][w]
        own = (q[w]-80)/20
        deficit = c.get("intercept", 0)+c.get("own", 0)*own+c.get("Q", 0)*Q
        if w == "INJECT3":
            deficit += c["own2"]*interpolate(own, p["pwl_functions"]["own2_INJECT3"])
            deficit += c["Q2"]*interpolate(Q, p["pwl_functions"]["Q2"])
        v["deficit"][w] = deficit
        assert abs(r[w] - (q[w]-min(max(deficit, 0), q[w]))) < 1e-4
    R = sum(r.values())
    Z = (R-600)/20
    fc = p["layer2"]["FOPT"]["coefficients"]
    fopt = fc["intercept"] + sum(fc[f"r{i}"]*(r[w]-80)/20
                                   for i, w in enumerate(wells, 1))
    fopt += fc["Z2"]*interpolate(Z, p["pwl_functions"]["Z2"])
    wc = p["layer2"]["FWIT"]["coefficients"]
    fwit = wc["intercept"] + wc["Z"]*Z
    assert abs(R-v["R"]) < 1e-4
    assert abs(fopt-v["FOPT"]) < 1e-3
    assert abs(fwit-v["FWIT"]) < 1e-3
    assert sum(q.values()) <= p["decision_variables"]["capacity_requested"]+1e-5
    assert fwit <= v["epsilon_FWIT"]+1e-4
    v["audit"] = {"Q": Q, "Z": Z, "max_clip_residual": max(
        abs(r[w]-(q[w]-min(max(v["deficit"][w], 0), q[w]))) for w in wells),
        "FOPT_residual": v["FOPT"]-fopt, "FWIT_residual": v["FWIT"]-fwit}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epsilon", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if not math.isfinite(args.epsilon):
        ap.error("epsilon must be finite")
    p = json.loads(PROBLEM.read_text())
    m, idx = build(p, args.epsilon)
    result = m.solve({idx["FOPT"]: -1})
    if result.x is None or result.status not in (0, 1):
        raise RuntimeError(f"MILP failed: status={result.status}, message={result.message}")
    x = result.x
    val = lambda k: float(x[k])
    out = {"status": "optimal" if result.status == 0 else "limit_reached",
           "epsilon_FWIT": args.epsilon,
           "requested": {w: val(j) for w, j in idx["q"].items()},
           "realized": {w: val(j) for w, j in idx["r"].items()},
           "deficit": {}, "R": val(idx["R"]), "FOPT": val(idx["FOPT"]),
           "FWIT": val(idx["FWIT"]), "solver_gap": float(result.mip_gap),
           "solver_message": result.message,
           "regimes": {w: [round(val(j)) for j in js]
                       for w, js in idx["regimes"].items()}}
    audit(p, out)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
