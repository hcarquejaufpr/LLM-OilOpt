#!/usr/bin/env python3
"""Solve the public P2 v3 water-injection MILP with SciPy/HiGHS."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


ROOT = Path(__file__).resolve().parent
PUBLIC = ROOT / "shared" / "problema_p2_public.json"


class Model:
    def __init__(self):
        self.names, self.lb, self.ub, self.integer = [], [], [], []
        self.rows, self.clb, self.cub = [], [], []

    def var(self, name, lo=-math.inf, hi=math.inf, binary=False):
        j = len(self.names)
        self.names.append(name)
        self.lb.append(0.0 if binary else lo)
        self.ub.append(1.0 if binary else hi)
        self.integer.append(1 if binary else 0)
        return j

    def row(self, terms, lo=-math.inf, hi=math.inf):
        self.rows.append({k: v for k, v in terms.items() if v != 0})
        self.clb.append(lo)
        self.cub.append(hi)

    def eq(self, terms, rhs):
        self.row(terms, rhs, rhs)

    def pwl(self, name, x, nodes, values):
        if len(nodes) != len(values) or len(nodes) < 2 or any(
            a >= b for a, b in zip(nodes, nodes[1:])
        ):
            raise ValueError(f"Invalid PWL nodes for {name}")
        y = self.var(name, min(values), max(values))
        lam = [self.var(f"{name}.lambda[{j}]", 0, 1) for j in range(len(nodes))]
        seg = [self.var(f"{name}.segment[{j}]", binary=True)
               for j in range(len(nodes) - 1)]
        self.eq({j: 1 for j in lam}, 1)
        self.eq({x: 1, **{j: -v for j, v in zip(lam, nodes)}}, 0)
        self.eq({y: 1, **{j: -v for j, v in zip(lam, values)}}, 0)
        self.eq({j: 1 for j in seg}, 1)
        for j, v in enumerate(lam):
            allowed = ([seg[j - 1]] if j else []) + ([seg[j]] if j < len(seg) else [])
            self.row({v: 1, **{s: -1 for s in allowed}}, hi=0)
        return y

    def clip(self, name, q, r, d, qlo, qhi, dlo, dhi):
        """Exact disjunction; auxiliary deficit bounds set M only, never d bounds."""
        low, mid, high = [self.var(f"{name}.regime.{s}", binary=True)
                          for s in ("low", "middle", "high")]
        self.eq({low: 1, mid: 1, high: 1}, 1)
        # Global 0 <= r <= q is an implication of clipping, not an extra condition.
        self.row({r: 1, q: -1}, hi=0)

        # d <= 0 and r = q when low=1.
        self.row({d: 1, low: max(dhi, 0)}, hi=max(dhi, 0))
        self.row({q: 1, r: -1, low: qhi}, hi=qhi)

        # 0 <= d <= q and r = q-d when mid=1.
        self.row({d: 1, mid: min(dlo, 0)}, lo=min(dlo, 0))
        m = max(dhi - qlo, 0)
        self.row({d: 1, q: -1, mid: m}, hi=m)
        mpos = max(dhi, 0)
        mneg = max(qhi - dlo, 0)
        self.row({r: 1, q: -1, d: 1, mid: mpos}, hi=mpos)
        self.row({r: 1, q: -1, d: 1, mid: -mneg}, lo=-mneg)

        # d >= q and r = 0 when high=1.
        m = max(qhi - dlo, 0)
        self.row({d: 1, q: -1, high: -m}, lo=-m)
        self.row({r: 1, high: qhi}, hi=qhi)

    def solve(self, objective):
        ii, jj, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, v in row.items():
                ii.append(i); jj.append(j); vv.append(v)
        matrix = coo_matrix((vv, (ii, jj)), shape=(len(self.rows), len(self.names))).tocsr()
        c = np.zeros(len(self.names))
        c[objective] = -1 / 10000  # Scale objective for HiGHS numerics.
        return milp(c, integrality=np.array(self.integer),
                    bounds=Bounds(self.lb, self.ub),
                    constraints=LinearConstraint(matrix, self.clb, self.cub),
                    options={"mip_rel_gap": 1e-8})


def build(data, epsilon):
    m = Model()
    wells = data["wells"]
    fixed = data["decision_variables"]["fixed"]
    domains = data["decision_variables"]["domains"]
    qbounds = {w: [fixed[w], fixed[w]] if w in fixed else domains[w] for w in wells}
    q = {w: m.var(f"q.{w}", *qbounds[w]) for w in wells}
    r = {w: m.var(f"r.{w}", 0, qbounds[w][1]) for w in wells}
    d = {w: m.var(f"deficit.{w}") for w in wells}
    m.row({v: 1 for v in q.values()}, hi=data["decision_variables"]["capacity_requested"])

    Q = m.var("Q")
    m.eq({Q: 20, **{v: -1 for v in q.values()}}, -600)
    own = {}
    real = {}
    for w in wells:
        own[w] = m.var(f"own.{w}")
        real[w] = m.var(f"real.{w}")
        m.eq({own[w]: 20, q[w]: -1}, -80)
        m.eq({real[w]: 20, r[w]: -1}, -80)
    Z = m.var("Z")
    m.eq({Z: 20, **{v: -1 for v in r.values()}}, -600)

    pwl_data = data["pwl_functions"]
    pwl = {}
    for key, x in (("own2_INJECT3", own["INJECT3"]), ("Q2", Q), ("Z2", Z)):
        p = pwl_data[key]
        pwl[key] = m.pwl(key, x, p["nodes"], p["values"])

    for w in wells:
        co = data["layer1_deficit_models"][w]
        terms = {d[w]: 1}
        for key, v in co.items():
            if key == "intercept":
                continue
            feature = {"Q": Q, "own": own[w], "own2": pwl["own2_INJECT3"],
                       "Q2": pwl["Q2"]}[key]
            terms[feature] = terms.get(feature, 0) - v
        m.eq(terms, co.get("intercept", 0))
        m.clip(w, q[w], r[w], d[w], *qbounds[w], *data["deficit_bounds"][w])

    fopt = m.var("FOPT")
    co = data["layer2"]["FOPT"]["coefficients"]
    terms = {fopt: 1, pwl["Z2"]: -co["Z2"]}
    for i, w in enumerate(wells, 1):
        terms[real[w]] = -co[f"r{i}"]
    m.eq(terms, co["intercept"])
    fwit = m.var("FWIT")
    co = data["layer2"]["FWIT"]["coefficients"]
    m.eq({fwit: 1, Z: -co["Z"]}, co["intercept"])
    m.row({fwit: 1}, hi=epsilon)
    return m, q, r, d, Q, Z, fopt, fwit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    data = json.loads(PUBLIC.read_text())
    m, q, r, d, Q, Z, fopt, fwit = build(data, args.epsilon)
    result = m.solve(fopt)
    status = {0: "optimal", 1: "limit", 2: "infeasible", 3: "unbounded", 4: "error"}.get(result.status, "error")
    output = {"status": status, "epsilon_FWIT": args.epsilon,
              "requested": None, "realized": None, "R": None,
              "FOPT": None, "FWIT": None,
              "solver_gap": float(result.mip_gap) if getattr(result, "mip_gap", None) is not None else None,
              "solver_message": result.message}
    if result.x is not None:
        x = result.x
        output.update(requested={w: float(x[j]) for w, j in q.items()},
                      realized={w: float(x[j]) for w, j in r.items()},
                      R=float(sum(x[j] for j in r.values())),
                      FOPT=float(x[fopt]), FWIT=float(x[fwit]),
                      Q=float(x[Q]), Z=float(x[Z]),
                      deficits={w: float(x[j]) for w, j in d.items()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    if status != "optimal":
        raise SystemExit(f"Solver status: {status}: {result.message}")


if __name__ == "__main__":
    main()
