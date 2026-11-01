#!/usr/bin/env python3
"""Solve the public P2 v3 mixed integer linear model."""

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
        self.names, self.lower, self.upper, self.integer = [], [], [], []
        self.rows, self.row_lb, self.row_ub = [], [], []

    def var(self, name, lb=-math.inf, ub=math.inf, integer=False):
        j = len(self.names)
        self.names.append(name)
        self.lower.append(lb)
        self.upper.append(ub)
        self.integer.append(int(integer))
        return j

    def row(self, coefficients, lb=-math.inf, ub=math.inf):
        self.rows.append({j: v for j, v in coefficients.items() if v})
        self.row_lb.append(lb)
        self.row_ub.append(ub)

    def eq(self, coefficients, rhs):
        self.row(coefficients, rhs, rhs)

    def le(self, coefficients, rhs):
        self.row(coefficients, ub=rhs)

    def ge(self, coefficients, rhs):
        self.row(coefficients, lb=rhs)

    def solve(self, objective):
        rr, cc, data = [], [], []
        for i, row in enumerate(self.rows):
            for j, value in row.items():
                rr.append(i)
                cc.append(j)
                data.append(value)
        matrix = coo_matrix((data, (rr, cc)), shape=(len(self.rows), len(self.names))).tocsr()
        c = np.zeros(len(self.names))
        c[objective] = -0.001  # Scale the large FOPT objective for HiGHS.
        result = milp(c, integrality=np.array(self.integer),
                      bounds=Bounds(self.lower, self.upper),
                      constraints=LinearConstraint(matrix, self.row_lb, self.row_ub),
                      options={"mip_rel_gap": 1e-9})
        return result


def pwl(model, name, argument, spec):
    """Continuous lambda weights with exactly one adjacent segment selected."""
    nodes, values = spec["nodes"], spec["values"]
    n = len(nodes)
    assert n == len(values) and n >= 2
    assert all(nodes[k] < nodes[k + 1] for k in range(n - 1))
    y = model.var(name)
    lam = [model.var(f"{name}.lambda.{j}", lb=0, ub=1) for j in range(n)]
    seg = [model.var(f"{name}.segment.{j}", lb=0, ub=1, integer=True)
           for j in range(n - 1)]
    model.eq({j: 1 for j in lam}, 1)
    model.eq({j: 1 for j in seg}, 1)
    model.eq({argument: 1, **{j: -x for j, x in zip(lam, nodes)}}, 0)
    model.eq({y: 1, **{j: -v for j, v in zip(lam, values)}}, 0)
    for j, weight in enumerate(lam):
        row = {weight: 1}
        if j > 0:
            row[seg[j - 1]] = -1
        if j < n - 1:
            row[seg[j]] = -1
        model.le(row, 0)
    return y


def big_m_values(deficit_bounds, q_bounds):
    """Sufficient implication constants, using bounds only as certificates."""
    L, U = deficit_bounds
    qmin, qmax = q_bounds
    assert L <= U and 0 <= qmin <= qmax
    return {
        "d_le_zero": max(U, 0),
        "d_ge_zero": max(-L, 0),
        "d_le_q": max(U - qmin, 0),
        "d_ge_q": max(qmax - L, 0),
        "t_zero": qmax,
        "t_minus_d": max(qmax - L, 0),
        "d_minus_t": max(U, 0),
        "t_equals_q": qmax,
    }


def add_clip(model, well, q, d, r, q_bounds, deficit_bounds):
    """t=clip(d,0,q), r=q-t; auxiliary bounds never constrain d or r."""
    m = big_m_values(deficit_bounds, q_bounds)
    t = model.var(f"clip.{well}", lb=0)  # t >= 0 is intrinsic to clipping.
    a, b, c = [model.var(f"regime.{well}.{s}", lb=0, ub=1, integer=True)
               for s in ("lower", "middle", "upper")]
    model.eq({a: 1, b: 1, c: 1}, 1)
    model.eq({r: 1, q: -1, t: 1}, 0)
    model.le({t: 1, q: -1}, 0)  # t <= q is intrinsic to clipping.

    # lower: d <= 0, t = 0
    model.le({d: 1, a: m["d_le_zero"]}, m["d_le_zero"])
    model.le({t: 1, a: m["t_zero"]}, m["t_zero"])
    # middle: 0 <= d <= q, t = d
    model.ge({d: 1, b: -m["d_ge_zero"]}, -m["d_ge_zero"])
    model.le({d: 1, q: -1, b: m["d_le_q"]}, m["d_le_q"])
    model.le({t: 1, d: -1, b: m["t_minus_d"]}, m["t_minus_d"])
    model.le({d: 1, t: -1, b: m["d_minus_t"]}, m["d_minus_t"])
    # upper: d >= q, t = q
    model.ge({d: 1, q: -1, c: -m["d_ge_q"]}, -m["d_ge_q"])
    model.le({q: 1, t: -1, c: m["t_equals_q"]}, m["t_equals_q"])
    return m


def build(problem, epsilon):
    model = Model()
    wells = problem["wells"]
    decisions = problem["decision_variables"]
    q, r, own, real, deficit = {}, {}, {}, {}, {}
    q_bounds = {}
    for well in wells:
        bounds = ([decisions["fixed"][well]] * 2 if well in decisions["fixed"]
                  else decisions["domains"][well])
        q_bounds[well] = bounds
        q[well] = model.var(f"q.{well}", lb=bounds[0], ub=bounds[1])
        r[well] = model.var(f"r.{well}")
        own[well] = model.var(f"own.{well}")
        real[well] = model.var(f"real.{well}")
        deficit[well] = model.var(f"deficit.{well}")
        model.eq({own[well]: 20, q[well]: -1}, -80)
        model.eq({real[well]: 20, r[well]: -1}, -80)

    model.le({j: 1 for j in q.values()}, decisions["capacity_requested"])
    Q = model.var("Q")
    model.eq({Q: 20, **{j: -1 for j in q.values()}}, -600)
    R = model.var("R")
    model.eq({R: 1, **{j: -1 for j in r.values()}}, 0)
    Z = model.var("Z")
    model.eq({Z: 20, R: -1}, -600)
    funcs = problem["pwl_functions"]
    own2 = pwl(model, "own2_INJECT3", own["INJECT3"], funcs["own2_INJECT3"])
    Q2 = pwl(model, "Q2", Q, funcs["Q2"])
    Z2 = pwl(model, "Z2", Z, funcs["Z2"])

    big_m = {}
    for well in wells:
        coeffs = problem["layer1_deficit_models"][well]
        row = {deficit[well]: 1}
        for key, feature in (("own", own[well]), ("Q", Q),
                             ("own2", own2), ("Q2", Q2)):
            if key in coeffs:
                row[feature] = -coeffs[key]
        model.eq(row, coeffs.get("intercept", 0))
        big_m[well] = add_clip(model, well, q[well], deficit[well], r[well],
                               q_bounds[well], problem["deficit_bounds"][well])

    fopt = model.var("FOPT")
    fc = problem["layer2"]["FOPT"]["coefficients"]
    row = {fopt: 1, Z2: -fc["Z2"]}
    for i, well in enumerate(wells, 1):
        row[real[well]] = -fc[f"r{i}"]
    model.eq(row, fc["intercept"])
    fwit = model.var("FWIT")
    wc = problem["layer2"]["FWIT"]["coefficients"]
    model.eq({fwit: 1, Z: -wc["Z"]}, wc["intercept"])
    model.le({fwit: 1}, epsilon)
    return model, {"q": q, "r": r, "own": own, "real": real,
                   "deficit": deficit, "Q": Q, "R": R, "Z": Z,
                   "own2_INJECT3": own2, "Q2": Q2, "Z2": Z2,
                   "FOPT": fopt, "FWIT": fwit}, big_m


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    problem = json.loads(PROBLEM.read_text())
    model, index, big_m = build(problem, args.epsilon)
    result = model.solve(index["FOPT"])
    output = {"status": "optimal" if result.status == 0 else "not_optimal",
              "epsilon_FWIT": args.epsilon,
              "solver_status": int(result.status), "solver_message": result.message,
              "solver_gap": None if result.mip_gap is None else float(result.mip_gap),
              "big_m": big_m,
              "model_size": {"variables": len(model.names),
                             "constraints": len(model.rows)}}
    if result.x is not None:
        x = result.x
        output.update({"requested": {w: float(x[j]) for w, j in index["q"].items()},
                       "realized": {w: float(x[j]) for w, j in index["r"].items()},
                       "deficit": {w: float(x[j]) for w, j in index["deficit"].items()},
                       "features": {key: float(x[index[key]]) for key in
                                    ("Q", "Z", "own2_INJECT3", "Q2", "Z2")},
                       "R": float(x[index["R"]]),
                       "FOPT": float(x[index["FOPT"]]),
                       "FWIT": float(x[index["FWIT"]])})
        output["max_clip_residual"] = max(abs(output["realized"][w] -
            (output["requested"][w] - min(max(output["deficit"][w], 0),
                                           output["requested"][w])))
            for w in problem["wells"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    if result.status != 0:
        raise SystemExit(f"Solver did not prove optimality: {result.message}")


if __name__ == "__main__":
    main()
