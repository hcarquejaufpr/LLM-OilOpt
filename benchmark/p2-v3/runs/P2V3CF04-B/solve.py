#!/usr/bin/env python3
"""MILP exato para a instância pública P2 v3."""

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
        self.names, self.lo, self.hi, self.integer, self.cost = [], [], [], [], []
        self.rows, self.lower, self.upper = [], [], []

    def var(self, name, lo=-math.inf, hi=math.inf, binary=False, cost=0.0):
        k = len(self.names)
        self.names.append(name)
        self.lo.append(0.0 if binary else lo)
        self.hi.append(1.0 if binary else hi)
        self.integer.append(1 if binary else 0)
        self.cost.append(cost)
        return k

    def row(self, terms, lo=-math.inf, hi=math.inf):
        terms = {k: v for k, v in terms.items() if v != 0}
        self.rows.append(terms)
        self.lower.append(lo)
        self.upper.append(hi)

    def eq(self, terms, rhs):
        self.row(terms, rhs, rhs)

    def solve(self):
        rr, cc, vv = [], [], []
        for i, terms in enumerate(self.rows):
            for k, v in terms.items():
                rr.append(i)
                cc.append(k)
                vv.append(v)
        a = coo_matrix((vv, (rr, cc)), shape=(len(self.rows), len(self.names))).tocsr()
        return milp(np.array(self.cost), integrality=np.array(self.integer),
                    bounds=Bounds(self.lo, self.hi),
                    constraints=LinearConstraint(a, self.lower, self.upper),
                    options={"mip_rel_gap": 1e-9})


def pwl(model, name, arg, spec):
    xs, ys = spec["nodes"], spec["values"]
    assert len(xs) == len(ys) and len(xs) >= 2
    assert all(a < b for a, b in zip(xs, xs[1:]))
    lam = [model.var(f"lambda_{name}_{j}", lo=0) for j in range(len(xs))]
    seg = [model.var(f"segment_{name}_{j}", binary=True) for j in range(len(xs)-1)]
    value = model.var(name)
    model.eq({k: 1 for k in lam}, 1)
    model.eq({k: 1 for k in seg}, 1)
    model.eq({arg: 1, **{k: -v for k, v in zip(lam, xs)}}, 0)
    model.eq({value: 1, **{k: -v for k, v in zip(lam, ys)}}, 0)
    for j, k in enumerate(lam):
        row = {k: 1}
        if j > 0:
            row[seg[j-1]] = -1
        if j < len(seg):
            row[seg[j]] = -1
        model.row(row, hi=0)
    return value


def clip(model, name, q, d, qlo, qhi, dlo, dhi):
    """Three disjoint interiors; shared boundary points have the same r."""
    r = model.var(f"r_{name}")
    z0, z1, z2 = (model.var(f"clip_{name}_{j}", binary=True) for j in range(3))
    model.eq({z0: 1, z1: 1, z2: 1}, 1)

    # Each inactive inequality is redundant by d in [dlo,dhi] and q in [qlo,qhi].
    model.row({d: 1, z0: max(0, dhi)}, hi=max(0, dhi))  # z0 => d <= 0
    model.row({d: 1, z1: min(0, dlo)}, lo=min(0, dlo))  # z1 => d >= 0
    m_up = max(0, dhi-qlo)
    model.row({d: 1, q: -1, z1: m_up}, hi=m_up)  # z1 => d <= q
    m_low = max(0, qhi-dlo)
    model.row({d: 1, q: -1, z2: -m_low}, lo=-m_low)  # z2 => d >= q

    # z0 => r=q; inactive regimes imply |r-q| <= qhi.
    model.row({r: 1, q: -1, z0: qhi}, hi=qhi)
    model.row({r: 1, q: -1, z0: -qhi}, lo=-qhi)
    # z1 => r=q-d; inactive z0 gives r-q+d=d, z2 gives d-q.
    m_mid = max(abs(dlo), abs(dhi), qhi+max(0, -dlo))
    model.row({r: 1, q: -1, d: 1, z1: m_mid}, hi=m_mid)
    model.row({r: 1, q: -1, d: 1, z1: -m_mid}, lo=-m_mid)
    # z2 => r=0; other regimes imply 0 <= r <= qhi.
    model.row({r: 1, z2: qhi}, hi=qhi)
    model.row({r: 1, z2: -qhi}, lo=-qhi)
    return r


def construct(p, epsilon):
    m = Model()
    wells = p["wells"]
    domains = {**{w: [v, v] for w, v in p["decision_variables"]["fixed"].items()},
               **p["decision_variables"]["domains"]}
    q = {w: m.var(f"q_{w}", *domains[w]) for w in wells}
    m.row({q[w]: 1 for w in wells}, hi=p["decision_variables"]["capacity_requested"])
    Q = m.var("Q")
    m.eq({Q: 20, **{q[w]: -1 for w in wells}}, -600)
    own = {w: m.var(f"own_{w}") for w in wells}
    for w in wells:
        m.eq({own[w]: 20, q[w]: -1}, -80)
    own2 = pwl(m, "own2_INJECT3", own["INJECT3"], p["pwl_functions"]["own2_INJECT3"])
    Q2 = pwl(m, "Q2", Q, p["pwl_functions"]["Q2"])

    deficit, r = {}, {}
    for w in wells:
        coeff = p["layer1_deficit_models"][w]
        deficit[w] = m.var(f"deficit_{w}")
        terms = {deficit[w]: 1, own[w]: -coeff.get("own", 0), Q: -coeff.get("Q", 0)}
        if w == "INJECT3":
            terms[own2] = -coeff["own2"]
            terms[Q2] = -coeff["Q2"]
        m.eq(terms, coeff.get("intercept", 0))
        if not coeff:
            r[w] = m.var(f"r_{w}")
            m.eq({r[w]: 1, q[w]: -1}, 0)
        else:
            r[w] = clip(m, w, q[w], deficit[w], *domains[w],
                        *p["deficit_bounds"][w])

    R, Z = m.var("R"), m.var("Z")
    m.eq({R: 1, **{r[w]: -1 for w in wells}}, 0)
    m.eq({Z: 20, R: -1}, -600)
    Z2 = pwl(m, "Z2", Z, p["pwl_functions"]["Z2"])
    real = {w: m.var(f"real_{w}") for w in wells}
    for w in wells:
        m.eq({real[w]: 20, r[w]: -1}, -80)
    fc = p["layer2"]["FOPT"]["coefficients"]
    FOPT = m.var("FOPT", cost=-1)
    m.eq({FOPT: 1, **{real[w]: -fc[f"r{i}"] for i, w in enumerate(wells, 1)},
          Z2: -fc["Z2"]}, fc["intercept"])
    wc = p["layer2"]["FWIT"]["coefficients"]
    FWIT = m.var("FWIT")
    m.eq({FWIT: 1, Z: -wc["Z"]}, wc["intercept"])
    m.row({FWIT: 1}, hi=epsilon)
    return m, {"q": q, "r": r, "deficit": deficit, "R": R, "Z": Z,
               "FOPT": FOPT, "FWIT": FWIT, "Q": Q, "own2": own2,
               "Q2": Q2, "Z2": Z2}


def interp(x, spec):
    xs, ys = spec["nodes"], spec["values"]
    for j in range(len(xs)-1):
        if xs[j]-1e-7 <= x <= xs[j+1]+1e-7:
            t = (x-xs[j])/(xs[j+1]-xs[j])
            return ys[j]+t*(ys[j+1]-ys[j])
    raise ValueError(f"PWL argument {x} outside public nodes")


def audit(p, indexes, x):
    wells = p["wells"]
    q = {w: float(x[indexes["q"][w]]) for w in wells}
    r = {w: float(x[indexes["r"][w]]) for w in wells}
    Q = (sum(q.values())-600)/20
    deficits = {}
    for w in wells:
        c = p["layer1_deficit_models"][w]
        d = c.get("intercept", 0) + c.get("own", 0)*(q[w]-80)/20 + c.get("Q", 0)*Q
        if w == "INJECT3":
            d += c["own2"]*interp((q[w]-80)/20, p["pwl_functions"]["own2_INJECT3"])
            d += c["Q2"]*interp(Q, p["pwl_functions"]["Q2"])
        deficits[w] = d
        if abs(r[w]-(q[w]-min(max(d, 0), q[w]))) > 1e-5:
            raise RuntimeError(f"clipping validation failed: {w}")
    R = sum(r.values())
    Z = (R-600)/20
    fc = p["layer2"]["FOPT"]["coefficients"]
    fopt = fc["intercept"] + sum(fc[f"r{i}"]*(r[w]-80)/20 for i,w in enumerate(wells,1))
    fopt += fc["Z2"]*interp(Z, p["pwl_functions"]["Z2"])
    wc = p["layer2"]["FWIT"]["coefficients"]
    fwit = wc["intercept"] + wc["Z"]*Z
    if abs(fopt-x[indexes["FOPT"]]) > 1e-4 or abs(fwit-x[indexes["FWIT"]]) > 1e-4:
        raise RuntimeError("layer 2 validation failed")
    return q, r, deficits, Q, R, Z, fopt, fwit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    p = json.loads(PROBLEM.read_text())
    model, idx = construct(p, args.epsilon)
    result = model.solve()
    if result.status != 0:
        raise RuntimeError(f"MILP did not reach optimality: {result.message}")
    q, r, deficits, Q, R, Z, fopt, fwit = audit(p, idx, result.x)
    if fwit > args.epsilon + 1e-5:
        raise RuntimeError("epsilon validation failed")
    payload = {"status": "optimal", "epsilon_FWIT": args.epsilon,
               "requested": q, "realized": r, "R": R, "FOPT": fopt,
               "FWIT": fwit, "solver_gap": float(result.mip_gap),
               "Q": Q, "Z": Z, "deficits": deficits,
               "pwl_values": {key: float(result.x[idx[key]]) for key in ("own2", "Q2", "Z2")},
               "solver": "scipy.optimize.milp (HiGHS)",
               "mip_dual_bound": float(result.mip_dual_bound),
               "mip_node_count": int(result.mip_node_count)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
