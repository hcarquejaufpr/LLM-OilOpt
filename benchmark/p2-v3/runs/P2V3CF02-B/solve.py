#!/usr/bin/env python3
"""Solve the public P2 v3 MILP; no reference solutions or OPM data are used."""

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
        self.rows, self.row_lower, self.row_upper = [], [], []

    def var(self, name, lower=-math.inf, upper=math.inf, binary=False):
        index = len(self.names)
        self.names.append(name)
        self.lower.append(0.0 if binary else lower)
        self.upper.append(1.0 if binary else upper)
        self.integer.append(1 if binary else 0)
        return index

    def row(self, terms, lower=-math.inf, upper=math.inf):
        self.rows.append(terms)
        self.row_lower.append(lower)
        self.row_upper.append(upper)

    def eq(self, terms, rhs):
        self.row(terms, rhs, rhs)

    def le(self, terms, rhs):
        self.row(terms, upper=rhs)

    def solve(self, objective):
        ri, ci, data = [], [], []
        for i, row in enumerate(self.rows):
            for j, coefficient in row.items():
                if coefficient:
                    ri.append(i)
                    ci.append(j)
                    data.append(coefficient)
        matrix = coo_matrix((data, (ri, ci)), shape=(len(self.rows), len(self.names))).tocsr()
        c = np.zeros(len(self.names))
        c[objective] = -1.0
        return milp(
            c,
            integrality=np.asarray(self.integer),
            bounds=Bounds(self.lower, self.upper),
            constraints=LinearConstraint(matrix, self.row_lower, self.row_upper),
            options={"mip_rel_gap": 1e-9},
        )


def add_pwl(model, name, argument, specification):
    """Select one adjacent public segment and interpolate within that segment."""
    nodes, values = specification["nodes"], specification["values"]
    assert len(nodes) == len(values) and len(nodes) >= 2
    assert all(b > a for a, b in zip(nodes, nodes[1:]))
    y = [model.var(f"{name}.segment.{j}", binary=True) for j in range(len(nodes) - 1)]
    t = [model.var(f"{name}.fraction.{j}", lower=0, upper=1) for j in range(len(y))]
    value = model.var(name)
    model.eq({j: 1 for j in y}, 1)
    for tj, yj in zip(t, y):
        model.le({tj: 1, yj: -1}, 0)
    xrow, vrow = {argument: 1}, {value: 1}
    for j, (yj, tj) in enumerate(zip(y, t)):
        xrow[yj] = -nodes[j]
        xrow[tj] = -(nodes[j + 1] - nodes[j])
        vrow[yj] = -values[j]
        vrow[tj] = -(values[j + 1] - values[j])
    model.eq(xrow, 0)
    model.eq(vrow, 0)
    return value


def add_clip(model, well, q, d, r, q_bounds, deficit_bounds):
    """Three exact regimes. Auxiliary bounds appear only in valid M coefficients."""
    qmin, qmax = q_bounds
    dmin, dmax = deficit_bounds
    assert 0 <= qmin <= qmax and dmin <= dmax
    # With 0 <= r <= q, r-q lies in [-qmax, 0]. The supplied deficit
    # interval bounds d only for calculating these implication constants.
    m_d_upper = max(0.0, dmax)
    m_d_lower = max(0.0, -dmin)
    m_d_minus_q_upper = max(0.0, dmax - qmin)
    m_d_minus_q_lower = max(0.0, qmax - dmin)
    m_middle_upper = max(0.0, dmax)  # r-q+d <= dmax
    m_middle_lower = max(0.0, qmax - dmin)  # r-q+d >= dmin-qmax
    z0, z1, z2 = [model.var(f"{well}.regime.{j}", binary=True) for j in range(3)]
    model.eq({z0: 1, z1: 1, z2: 1}, 1)
    model.le({r: 1, q: -1}, 0)  # implied by clipping, not an extra restriction
    # Regime 0: d <= 0, r = q.
    model.le({d: 1, z0: m_d_upper}, m_d_upper)
    model.le({r: -1, q: 1, z0: qmax}, qmax)
    # Regime 1: 0 <= d <= q, r = q-d.
    model.le({d: -1, z1: m_d_lower}, m_d_lower)
    model.le({d: 1, q: -1, z1: m_d_minus_q_upper}, m_d_minus_q_upper)
    model.le({r: 1, q: -1, d: 1, z1: m_middle_upper}, m_middle_upper)
    model.le({r: -1, q: 1, d: -1, z1: m_middle_lower}, m_middle_lower)
    # Regime 2: d >= q, r = 0.
    model.le({d: -1, q: 1, z2: m_d_minus_q_lower}, m_d_minus_q_lower)
    model.le({r: 1, z2: qmax}, qmax)
    return {
        "d_upper": m_d_upper, "d_lower": m_d_lower,
        "d_minus_q_upper": m_d_minus_q_upper,
        "d_minus_q_lower": m_d_minus_q_lower,
        "r_minus_q_lower": qmax,
        "middle_upper": m_middle_upper, "middle_lower": m_middle_lower,
        "r_zero_upper": qmax,
    }


def build(problem, epsilon):
    model = Model()
    wells = problem["wells"]
    fixed = problem["decision_variables"]["fixed"]
    domains = problem["decision_variables"]["domains"]
    q_bounds = {w: [fixed[w], fixed[w]] if w in fixed else domains[w] for w in wells}
    q = {w: model.var(f"q.{w}", *q_bounds[w]) for w in wells}
    r = {w: model.var(f"r.{w}", lower=0, upper=q_bounds[w][1]) for w in wells}
    own = {w: model.var(f"own.{w}") for w in wells}
    real = {w: model.var(f"real.{w}") for w in wells}
    Q = model.var("Q")
    R = model.var("R")
    Z = model.var("Z")
    model.le({j: 1 for j in q.values()}, problem["decision_variables"]["capacity_requested"])
    model.eq({Q: 20, **{j: -1 for j in q.values()}}, -600)
    model.eq({R: 1, **{j: -1 for j in r.values()}}, 0)
    model.eq({Z: 20, R: -1}, -600)
    for w in wells:
        model.eq({own[w]: 20, q[w]: -1}, -80)
        model.eq({real[w]: 20, r[w]: -1}, -80)
    pwl = problem["pwl_functions"]
    own2 = add_pwl(model, "own2_INJECT3", own["INJECT3"], pwl["own2_INJECT3"])
    Q2 = add_pwl(model, "Q2", Q, pwl["Q2"])
    Z2 = add_pwl(model, "Z2", Z, pwl["Z2"])
    deficits, big_m = {}, {}
    for w in wells:
        coefficients = problem["layer1_deficit_models"][w]
        d = model.var(f"deficit.{w}")  # deliberately unbounded
        deficits[w] = d
        row = {d: 1}
        for key, j in (("own", own[w]), ("Q", Q), ("own2", own2), ("Q2", Q2)):
            if key in coefficients:
                row[j] = -coefficients[key]
        model.eq(row, coefficients.get("intercept", 0))
        if coefficients:
            big_m[w] = add_clip(model, w, q[w], d, r[w], q_bounds[w], problem["deficit_bounds"][w])
        else:
            model.eq({r[w]: 1, q[w]: -1}, 0)
    fc = problem["layer2"]["FOPT"]["coefficients"]
    wc = problem["layer2"]["FWIT"]["coefficients"]
    fopt = model.var("FOPT")
    fwit = model.var("FWIT")
    model.eq({fopt: 1, **{real[w]: -fc[f"r{i}"] for i, w in enumerate(wells, 1)}, Z2: -fc["Z2"]}, fc["intercept"])
    model.eq({fwit: 1, Z: -wc["Z"]}, wc["intercept"])
    model.le({fwit: 1}, epsilon)
    return model, {"q": q, "r": r, "d": deficits, "Q": Q, "R": R, "Z": Z,
                   "FOPT": fopt, "FWIT": fwit, "own2_INJECT3": own2, "Q2": Q2, "Z2": Z2,
                   "big_m": big_m}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    problem = json.loads(PROBLEM.read_text())
    model, v = build(problem, args.epsilon)
    result = model.solve(v["FOPT"])
    if result.x is None or result.status != 0:
        raise RuntimeError(f"MILP was not solved to optimality: {result.message}")
    x = result.x
    wells = problem["wells"]
    output = {
        "status": "optimal", "epsilon_FWIT": args.epsilon,
        "requested": {w: float(x[v["q"][w]]) for w in wells},
        "realized": {w: float(x[v["r"][w]]) for w in wells},
        "deficit": {w: float(x[v["d"][w]]) for w in wells},
        "Q": float(x[v["Q"]]), "R": float(x[v["R"]]), "Z": float(x[v["Z"]]),
        "own2_INJECT3": float(x[v["own2_INJECT3"]]),
        "Q2": float(x[v["Q2"]]), "Z2": float(x[v["Z2"]]),
        "FOPT": float(x[v["FOPT"]]), "FWIT": float(x[v["FWIT"]]),
        "solver_gap": float(result.mip_gap),
        "solver": {"name": "scipy.optimize.milp/HiGHS", "status_code": int(result.status),
                   "message": result.message, "dual_bound": float(-result.mip_dual_bound),
                   "node_count": int(result.mip_node_count)},
        "big_m": v["big_m"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
