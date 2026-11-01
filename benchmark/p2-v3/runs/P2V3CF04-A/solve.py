#!/usr/bin/env python3
"""Solve the public P2 v3 water-injection MILP using SciPy/HiGHS."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


HERE = Path(__file__).resolve().parent
DATA = HERE / "shared" / "problema_p2_public.json"


class Model:
    def __init__(self):
        self.names, self.lb, self.ub, self.integrality = [], [], [], []
        self.rows, self.row_lb, self.row_ub = [], [], []

    def var(self, name, lb=-math.inf, ub=math.inf, binary=False):
        k = len(self.names)
        self.names.append(name)
        self.lb.append(0.0 if binary else lb)
        self.ub.append(1.0 if binary else ub)
        self.integrality.append(1 if binary else 0)
        return k

    def row(self, terms, lb=-math.inf, ub=math.inf):
        merged = {}
        for k, a in terms:
            merged[k] = merged.get(k, 0.0) + a
        self.rows.append({k: a for k, a in merged.items() if a})
        self.row_lb.append(lb)
        self.row_ub.append(ub)

    def eq(self, terms, rhs):
        self.row(terms, rhs, rhs)

    def pwl(self, label, x, function):
        nodes, values = function["nodes"], function["values"]
        assert len(nodes) == len(values) and len(nodes) >= 2
        assert all(a < b for a, b in zip(nodes, nodes[1:]))
        y = self.var(label)
        lambdas = [self.var(f"lambda_{label}_{j}", lb=0) for j in range(len(nodes))]
        segments = [self.var(f"segment_{label}_{j}", binary=True)
                    for j in range(len(nodes) - 1)]
        self.eq([(v, 1) for v in lambdas], 1)
        self.eq([(v, 1) for v in segments], 1)
        self.eq([(x, 1)] + [(v, -a) for v, a in zip(lambdas, nodes)], 0)
        self.eq([(y, 1)] + [(v, -a) for v, a in zip(lambdas, values)], 0)
        for j, lam in enumerate(lambdas):
            adjacent = []
            if j:
                adjacent.append((segments[j - 1], -1))
            if j < len(segments):
                adjacent.append((segments[j], -1))
            self.row([(lam, 1)] + adjacent, ub=0)
        return y

    def solve(self, fopt):
        ii, jj, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, val in row.items():
                ii.append(i)
                jj.append(j)
                vv.append(val)
        matrix = coo_matrix((vv, (ii, jj)), shape=(len(self.rows), len(self.names))).tocsr()
        objective = np.zeros(len(self.names))
        objective[fopt] = -1
        return milp(objective, integrality=np.array(self.integrality),
                    bounds=Bounds(self.lb, self.ub),
                    constraints=LinearConstraint(matrix, self.row_lb, self.row_ub),
                    options={"mip_rel_gap": 1e-9})


def build(p, epsilon):
    m = Model()
    wells = p["wells"]
    dv = p["decision_variables"]
    q = {}
    for w in wells:
        lo, hi = (dv["fixed"][w],) * 2 if w in dv["fixed"] else dv["domains"][w]
        q[w] = m.var("q_" + w, lo, hi)
    m.row([(q[w], 1) for w in wells], ub=dv["capacity_requested"])

    Q = m.var("Q")
    m.eq([(Q, 20)] + [(q[w], -1) for w in wells], -600)
    own = {}
    for w in wells:
        own[w] = m.var("own_" + w)
        m.eq([(own[w], 20), (q[w], -1)], -80)
    functions = p["pwl_functions"]
    own2 = m.pwl("own2_INJECT3", own["INJECT3"], functions["own2_INJECT3"])
    Q2 = m.pwl("Q2", Q, functions["Q2"])

    d, r = {}, {}
    for w in wells:
        d[w] = m.var("deficit_" + w)  # no auxiliary deficit bound imposed
        c = p["layer1_deficit_models"][w]
        terms = [(d[w], 1), (Q, -c.get("Q", 0)),
                 (own[w], -c.get("own", 0))]
        if w == "INJECT3":
            terms += [(own2, -c.get("own2", 0)), (Q2, -c.get("Q2", 0))]
        m.eq(terms, c.get("intercept", 0))

        # h=max(d,0); g=max(d-q,0); r=q-h+g.
        # The two hinge binaries encode all three clipping regimes exactly.
        dl, du = p["deficit_bounds"][w]
        ql, qu = (dv["fixed"][w],) * 2 if w in dv["fixed"] else dv["domains"][w]
        el, eu = dl - qu, du - ql
        h = m.var("positive_deficit_" + w)
        g = m.var("excess_deficit_" + w)
        b = m.var("positive_switch_" + w, binary=True)
        a = m.var("excess_switch_" + w, binary=True)
        # h >= d, h >= 0; h <= d - dl*(1-b); h <= du*b
        m.row([(h, 1), (d[w], -1)], lb=0)
        m.row([(h, 1)], lb=0)
        m.row([(h, 1), (d[w], -1), (b, -dl)], ub=-dl)
        m.row([(h, 1), (b, -du)], ub=0)
        # g >= d-q, g >= 0; g <= d-q-el*(1-a); g <= eu*a
        m.row([(g, 1), (d[w], -1), (q[w], 1)], lb=0)
        m.row([(g, 1)], lb=0)
        m.row([(g, 1), (d[w], -1), (q[w], 1), (a, -el)], ub=-el)
        m.row([(g, 1), (a, -eu)], ub=0)
        r[w] = m.var("r_" + w)  # no auxiliary realized bound imposed
        m.eq([(r[w], 1), (q[w], -1), (h, 1), (g, -1)], 0)

    Z = m.var("Z")
    m.eq([(Z, 20)] + [(r[w], -1) for w in wells], -600)
    real = {}
    for w in wells:
        real[w] = m.var("real_" + w)
        m.eq([(real[w], 20), (r[w], -1)], -80)
    Z2 = m.pwl("Z2", Z, functions["Z2"])

    fopt = m.var("FOPT")
    fc = p["layer2"]["FOPT"]["coefficients"]
    m.eq([(fopt, 1), (Z2, -fc["Z2"])] +
         [(real[w], -fc[f"r{i}"]) for i, w in enumerate(wells, 1)],
         fc["intercept"])
    fwit = m.var("FWIT")
    wc = p["layer2"]["FWIT"]["coefficients"]
    m.eq([(fwit, 1), (Z, -wc["Z"])], wc["intercept"])
    m.row([(fwit, 1)], ub=epsilon)
    return m, q, r, d, Q, Z, fopt, fwit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("--epsilon must be finite")
    p = json.loads(DATA.read_text())
    m, q, r, d, Q, Z, fopt, fwit = build(p, args.epsilon)
    result = m.solve(fopt)
    status = {0: "optimal", 1: "limit_reached", 2: "infeasible",
              3: "unbounded", 4: "solver_error"}.get(result.status, "solver_error")
    output = {"status": status, "epsilon_FWIT": args.epsilon,
              "solver_gap": float(result.mip_gap) if result.mip_gap is not None else None,
              "solver_message": result.message}
    if result.x is not None:
        x = result.x
        output.update({
            "requested": {w: float(x[q[w]]) for w in p["wells"]},
            "realized": {w: float(x[r[w]]) for w in p["wells"]},
            "deficit": {w: float(x[d[w]]) for w in p["wells"]},
            "R": float(sum(x[r[w]] for w in p["wells"])),
            "Q": float(x[Q]), "Z": float(x[Z]),
            "FOPT": float(x[fopt]), "FWIT": float(x[fwit]),
            "objective_bound": float(-result.mip_dual_bound)
            if result.mip_dual_bound is not None else None})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    if status != "optimal":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
