#!/usr/bin/env python3
"""Auditable MILP for the public P2 v3 water-injection problem."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "shared" / "problema_p2_public.json"


class LinearMILP:
    def __init__(self):
        self.names, self.lower, self.upper, self.integrality, self.cost = [], [], [], [], []
        self.rows, self.row_lb, self.row_ub = [], [], []

    def var(self, name, lb=-np.inf, ub=np.inf, binary=False, cost=0.0):
        index = len(self.names)
        self.names.append(name)
        self.lower.append(0.0 if binary else lb)
        self.upper.append(1.0 if binary else ub)
        self.integrality.append(1 if binary else 0)
        self.cost.append(cost)
        return index

    def add(self, terms, lb=-np.inf, ub=np.inf):
        row = {}
        for index, coefficient in terms.items():
            row[index] = row.get(index, 0.0) + coefficient
        self.rows.append(row)
        self.row_lb.append(lb)
        self.row_ub.append(ub)

    def eq(self, terms, rhs=0.0):
        self.add(terms, rhs, rhs)

    def solve(self):
        row, col, value = [], [], []
        for i, terms in enumerate(self.rows):
            for j, coefficient in terms.items():
                if coefficient:
                    row.append(i)
                    col.append(j)
                    value.append(coefficient)
        matrix = coo_matrix((value, (row, col)), shape=(len(self.rows), len(self.names))).tocsr()
        return milp(
            np.asarray(self.cost),
            integrality=np.asarray(self.integrality),
            bounds=Bounds(self.lower, self.upper),
            constraints=LinearConstraint(matrix, self.row_lb, self.row_ub),
            options={"mip_rel_gap": 1e-9},
        )


def pwl(model, name, x, spec):
    nodes, values = spec["nodes"], spec["values"]
    assert len(nodes) == len(values) and len(nodes) >= 2
    assert all(a < b for a, b in zip(nodes, nodes[1:]))
    lam = [model.var(f"lambda_{name}_{j}", lb=0.0, ub=1.0) for j in range(len(nodes))]
    seg = [model.var(f"segment_{name}_{j}", binary=True) for j in range(len(nodes) - 1)]
    y = model.var(f"pwl_{name}")
    model.eq({j: 1.0 for j in lam}, 1.0)
    model.eq({j: 1.0 for j in seg}, 1.0)
    model.eq({x: 1.0, **{j: -v for j, v in zip(lam, nodes)}})
    model.eq({y: 1.0, **{j: -v for j, v in zip(lam, values)}})
    for j, weight in enumerate(lam):
        adjacent = ([seg[j - 1]] if j > 0 else []) + ([seg[j]] if j < len(seg) else [])
        model.add({weight: 1.0, **{k: -1.0 for k in adjacent}}, ub=0.0)
    return y


def add_clip(model, name, q, r, d, q_lo, q_hi, d_lo, d_hi):
    """Three exact implications, with constants from auxiliary deficit intervals.

    r has only the logical bounds 0 <= r <= q (q <= q_hi), not realized_bounds.
    Every M covers the corresponding inactive inequality over those bounds.
    """
    neg = model.var(f"clip_negative_{name}", binary=True)
    mid = model.var(f"clip_middle_{name}", binary=True)
    sat = model.var(f"clip_saturated_{name}", binary=True)
    model.eq({neg: 1, mid: 1, sat: 1}, 1)
    model.add({r: 1, q: -1}, ub=0)  # consequence of the clip identity

    # d <= 0, r = q when neg = 1.
    model.add({d: 1, neg: max(d_hi, 0)}, ub=max(d_hi, 0))
    model.add({r: 1, q: -1, neg: -q_hi}, lb=-q_hi)

    # 0 <= d <= q, r = q-d when mid = 1.
    m = max(-d_lo, 0)
    model.add({d: -1, mid: m}, ub=m)
    m = max(d_hi - q_lo, 0)
    model.add({d: 1, q: -1, mid: m}, ub=m)
    m_low = max(q_hi - d_lo, 0)
    m_high = max(d_hi, 0)
    model.add({r: 1, q: -1, d: 1, mid: -m_low}, lb=-m_low)
    model.add({r: 1, q: -1, d: 1, mid: m_high}, ub=m_high)

    # d >= q, r = 0 when sat = 1.
    m = max(q_hi - d_lo, 0)
    model.add({q: 1, d: -1, sat: m}, ub=m)
    model.add({r: 1, sat: q_hi}, ub=q_hi)


def build(data, epsilon):
    wells = data["wells"]
    domains = data["decision_variables"]["domains"]
    fixed = data["decision_variables"]["fixed"]
    model = LinearMILP()
    q, r, d, q_domain = {}, {}, {}, {}
    for well in wells:
        lo, hi = ([fixed[well]] * 2) if well in fixed else domains[well]
        q_domain[well] = (lo, hi)
        q[well] = model.var(f"q_{well}", lo, hi)
        r[well] = model.var(f"r_{well}", 0, hi)
        d[well] = model.var(f"deficit_{well}")

    capacity = data["decision_variables"]["capacity_requested"]
    model.add({q[w]: 1 for w in wells}, ub=capacity)
    Q = model.var("Q")
    Z = model.var("Z")
    model.eq({Q: 20, **{q[w]: -1 for w in wells}}, -600)
    model.eq({Z: 20, **{r[w]: -1 for w in wells}}, -600)
    own3 = model.var("own_INJECT3")
    model.eq({own3: 20, q["INJECT3"]: -1}, -80)
    own2 = pwl(model, "own2_INJECT3", own3, data["pwl_functions"]["own2_INJECT3"])
    Q2 = pwl(model, "Q2", Q, data["pwl_functions"]["Q2"])
    Z2 = pwl(model, "Z2", Z, data["pwl_functions"]["Z2"])

    for well in wells:
        co = data["layer1_deficit_models"][well]
        terms = {d[well]: 1.0}
        rhs = co.get("intercept", 0.0)
        for key, v in (("own", own3), ("Q", Q), ("own2", own2), ("Q2", Q2)):
            if key in co:
                terms[v] = terms.get(v, 0) - co[key]
        model.eq(terms, rhs)
        d_lo, d_hi = data["deficit_bounds"][well]
        q_lo, q_hi = q_domain[well]
        add_clip(model, well, q[well], r[well], d[well], q_lo, q_hi, d_lo, d_hi)

    fopt = model.var("FOPT", cost=-1)
    fwit = model.var("FWIT")
    f_co = data["layer2"]["FOPT"]["coefficients"]
    # real_i = (r_i - 80)/20; keep the constant shift in the equality.
    f_terms = {fopt: 1, Z2: -f_co["Z2"]}
    for i, well in enumerate(wells, 1):
        f_terms[r[well]] = -f_co[f"r{i}"] / 20
    model.eq(f_terms, f_co["intercept"] - 4 * sum(f_co[f"r{i}"] for i in range(1, 9)))
    w_co = data["layer2"]["FWIT"]["coefficients"]
    model.eq({fwit: 1, Z: -w_co["Z"]}, w_co["intercept"])
    model.add({fwit: 1}, ub=epsilon)
    return model, q, r, d, Q, Z, Z2, fopt, fwit


def interpolate(x, spec):
    nodes, values = spec["nodes"], spec["values"]
    for j in range(len(nodes) - 1):
        if nodes[j] - 1e-7 <= x <= nodes[j + 1] + 1e-7:
            t = (x - nodes[j]) / (nodes[j + 1] - nodes[j])
            return values[j] * (1 - t) + values[j + 1] * t
    raise ValueError(f"PWL argument outside public nodes: {x}")


def audit(data, epsilon, requested, realized, deficits, FOPT, FWIT):
    wells = data["wells"]
    qsum, rsum = sum(requested.values()), sum(realized.values())
    Q, Z = (qsum - 600) / 20, (rsum - 600) / 20
    own = (requested["INJECT3"] - 80) / 20
    funcs = data["pwl_functions"]
    own2 = interpolate(own, funcs["own2_INJECT3"])
    Q2 = interpolate(Q, funcs["Q2"])
    Z2 = interpolate(Z, funcs["Z2"])
    errors = {}
    for well in wells:
        co = data["layer1_deficit_models"][well]
        expected_d = sum(co.get(k, 0) * v for k, v in (("intercept", 1), ("own", own), ("Q", Q), ("own2", own2), ("Q2", Q2)))
        errors[f"deficit_{well}"] = abs(deficits[well] - expected_d)
        expected_r = requested[well] - min(max(expected_d, 0), requested[well])
        errors[f"clip_{well}"] = abs(realized[well] - expected_r)
    fco = data["layer2"]["FOPT"]["coefficients"]
    expected_fopt = fco["intercept"] + sum(fco[f"r{i}"] * (realized[w] - 80) / 20 for i, w in enumerate(wells, 1)) + fco["Z2"] * Z2
    wco = data["layer2"]["FWIT"]["coefficients"]
    expected_fwit = wco["intercept"] + wco["Z"] * Z
    errors["FOPT"] = abs(FOPT - expected_fopt)
    errors["FWIT"] = abs(FWIT - expected_fwit)
    return {"max_equation_error": max(errors.values()), "equation_errors": errors,
            "requested_capacity_slack": data["decision_variables"]["capacity_requested"] - qsum,
            "epsilon_slack": epsilon - FWIT,
            "Q": Q, "Z": Z, "PWL": {"own2_INJECT3": own2, "Q2": Q2, "Z2": Z2}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    data = json.loads(DATA.read_text())
    model, q, r, d, Q, Z, Z2, fopt, fwit = build(data, args.epsilon)
    solution = model.solve()
    status = {0: "optimal", 1: "limit_reached", 2: "infeasible", 3: "unbounded", 4: "solver_error"}.get(solution.status, "solver_error")
    result = {"status": status, "epsilon_FWIT": args.epsilon,
              "solver_gap": float(solution.mip_gap) if getattr(solution, "mip_gap", None) is not None else None,
              "solver_message": solution.message,
              "model_size": {"variables": len(model.names), "constraints": len(model.rows)}}
    if solution.x is not None:
        x = solution.x
        wells = data["wells"]
        requested = {w: float(x[q[w]]) for w in wells}
        realized = {w: float(x[r[w]]) for w in wells}
        deficits = {w: float(x[d[w]]) for w in wells}
        result.update({"requested": requested, "realized": realized, "R": float(sum(realized.values())),
                       "FOPT": float(x[fopt]), "FWIT": float(x[fwit]), "deficits": deficits,
                       "audit": audit(data, args.epsilon, requested, realized, deficits, float(x[fopt]), float(x[fwit]))})
    else:
        result.update({"requested": None, "realized": None, "R": None, "FOPT": None, "FWIT": None})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if status not in ("optimal", "infeasible", "unbounded"):
        raise SystemExit(f"Solver did not finish normally: {solution.message}")


if __name__ == "__main__":
    main()
