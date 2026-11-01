#!/usr/bin/env python3
"""Solve the public P2 v3 MILP, without external physical constraints."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


ROOT = Path(__file__).resolve().parent
PROBLEM = ROOT / "shared" / "problema_p2_public.json"
INF = math.inf


class Model:
    def __init__(self):
        self.names, self.lo, self.hi, self.integer = [], [], [], []
        self.rows, self.cols, self.data = [], [], []
        self.clo, self.chi, self.labels = [], [], []

    def var(self, name, lo=-INF, hi=INF, binary=False):
        j = len(self.names)
        self.names.append(name)
        self.lo.append(0.0 if binary else lo)
        self.hi.append(1.0 if binary else hi)
        self.integer.append(1 if binary else 0)
        return j

    def con(self, name, terms, lo=-INF, hi=INF):
        row = len(self.clo)
        for col, value in terms.items():
            if value:
                self.rows.append(row)
                self.cols.append(col)
                self.data.append(value)
        self.clo.append(lo)
        self.chi.append(hi)
        self.labels.append(name)

    def pwl(self, name, arg, spec):
        xs, ys = spec["nodes"], spec["values"]
        assert len(xs) == len(ys) and len(xs) > 1
        assert all(a < b for a, b in zip(xs, xs[1:]))
        value = self.var(name)
        weights = [self.var(f"{name}.lambda.{j}", lo=0) for j in range(len(xs))]
        segments = [self.var(f"{name}.segment.{j}", binary=True) for j in range(len(xs) - 1)]
        self.con(f"{name}.sum_lambda", {j: 1 for j in weights}, 1, 1)
        self.con(f"{name}.one_segment", {j: 1 for j in segments}, 1, 1)
        self.con(f"{name}.argument", {arg: 1, **{j: -v for j, v in zip(weights, xs)}}, 0, 0)
        self.con(f"{name}.value", {value: 1, **{j: -v for j, v in zip(weights, ys)}}, 0, 0)
        for j, weight in enumerate(weights):
            near = ([segments[j - 1]] if j else []) + ([segments[j]] if j < len(segments) else [])
            self.con(f"{name}.adjacency.{j}", {weight: 1, **{k: -1 for k in near}}, hi=0)
        return value

    def solve(self, objective):
        matrix = coo_matrix((self.data, (self.rows, self.cols)), shape=(len(self.clo), len(self.names))).tocsr()
        c = np.zeros(len(self.names))
        for j, v in objective.items():
            c[j] = -v
        result = milp(c, integrality=np.array(self.integer), bounds=Bounds(self.lo, self.hi),
                      constraints=LinearConstraint(matrix, self.clo, self.chi),
                      options={"mip_rel_gap": 1e-9})
        return result


def build(problem, epsilon):
    m = Model()
    wells = problem["wells"]
    decision = problem["decision_variables"]
    qb = {w: [float(decision["fixed"][w])] * 2 if w in decision["fixed"]
          else list(map(float, decision["domains"][w])) for w in wells}
    q = {w: m.var(f"q.{w}", *qb[w]) for w in wells}
    r = {w: m.var(f"r.{w}", lo=0, hi=qb[w][1]) for w in wells}
    # These broad r bounds follow directly from 0 <= clip(d,0,q) <= q;
    # the supplied realized_bounds are never used as variable bounds.
    capacity = float(decision["capacity_requested"])
    m.con("capacity_requested", {q[w]: 1 for w in wells}, hi=capacity)
    Q = m.var("Q")
    m.con("Q.definition", {Q: 20, **{q[w]: -1 for w in wells}}, -600, -600)
    own = {w: m.var(f"own.{w}") for w in wells}
    for w in wells:
        m.con(f"own.{w}.definition", {own[w]: 20, q[w]: -1}, -80, -80)
    own2 = m.pwl("own2_INJECT3", own["INJECT3"], problem["pwl_functions"]["own2_INJECT3"])
    Q2 = m.pwl("Q2", Q, problem["pwl_functions"]["Q2"])
    d = {}
    big_m = {}
    for w in wells:
        d[w] = m.var(f"deficit.{w}")  # no deficit_bounds as feasibility bounds
        coeff = problem["layer1_deficit_models"][w]
        terms = {d[w]: 1}
        for key, col in [("own", own[w]), ("Q", Q), ("own2", own2), ("Q2", Q2)]:
            if key in coeff:
                terms[col] = -float(coeff[key])
        intercept = float(coeff.get("intercept", 0))
        m.con(f"deficit.{w}.definition", terms, intercept, intercept)
        if not coeff:
            m.con(f"realized.{w}.zero_deficit", {r[w]: 1, q[w]: -1}, 0, 0)
            continue
        L, U = map(float, problem["deficit_bounds"][w])
        qlo, qhi = qb[w]
        pad = 1e-8
        M = {
            "negative_d": max(0, U) + pad,
            "negative_rq": qhi + pad,
            "middle_d": max(0, -L) + pad,
            "middle_dq": max(0, U - qlo) + pad,
            "middle_rqd": max(abs(L - qhi), abs(qhi - qlo + U)) + pad,
            "saturated_dq": max(0, qhi - L) + pad,
            "saturated_r": qhi + pad,
        }
        big_m[w] = M
        b0, b1, b2 = [m.var(f"regime.{w}.{s}", binary=True) for s in ("negative", "middle", "saturated")]
        m.con(f"regime.{w}.one", {b0: 1, b1: 1, b2: 1}, 1, 1)
        # b0: d <= 0, r = q.
        m.con(f"clip.{w}.negative.d", {d[w]: 1, b0: M["negative_d"]}, hi=M["negative_d"])
        m.con(f"clip.{w}.negative.rq.upper", {r[w]: 1, q[w]: -1, b0: M["negative_rq"]}, hi=M["negative_rq"])
        m.con(f"clip.{w}.negative.rq.lower", {r[w]: 1, q[w]: -1, b0: -M["negative_rq"]}, lo=-M["negative_rq"])
        # b1: 0 <= d <= q, r = q-d.
        m.con(f"clip.{w}.middle.d.lower", {d[w]: 1, b1: -M["middle_d"]}, lo=-M["middle_d"])
        m.con(f"clip.{w}.middle.dq.upper", {d[w]: 1, q[w]: -1, b1: M["middle_dq"]}, hi=M["middle_dq"])
        m.con(f"clip.{w}.middle.rqd.upper", {r[w]: 1, q[w]: -1, d[w]: 1, b1: M["middle_rqd"]}, hi=M["middle_rqd"])
        m.con(f"clip.{w}.middle.rqd.lower", {r[w]: 1, q[w]: -1, d[w]: 1, b1: -M["middle_rqd"]}, lo=-M["middle_rqd"])
        # b2: d >= q, r = 0.
        m.con(f"clip.{w}.saturated.dq", {d[w]: 1, q[w]: -1, b2: -M["saturated_dq"]}, lo=-M["saturated_dq"])
        m.con(f"clip.{w}.saturated.r", {r[w]: 1, b2: M["saturated_r"]}, hi=M["saturated_r"])
    R = m.var("R")
    m.con("R.definition", {R: 1, **{r[w]: -1 for w in wells}}, 0, 0)
    Z = m.var("Z")
    m.con("Z.definition", {Z: 20, R: -1}, -600, -600)
    real = {w: m.var(f"real.{w}") for w in wells}
    for w in wells:
        m.con(f"real.{w}.definition", {real[w]: 20, r[w]: -1}, -80, -80)
    Z2 = m.pwl("Z2", Z, problem["pwl_functions"]["Z2"])
    fopt = problem["layer2"]["FOPT"]["coefficients"]
    objective = {real[w]: float(fopt[f"r{i}"]) for i, w in enumerate(wells, 1)}
    objective[Z2] = float(fopt["Z2"])
    fwit = problem["layer2"]["FWIT"]["coefficients"]
    m.con("FWIT.epsilon", {Z: float(fwit["Z"])}, hi=epsilon - float(fwit["intercept"]))
    return m, objective, {"q": q, "r": r, "d": d, "Q": Q, "R": R, "Z": Z,
                          "own": own, "real": real, "own2": own2, "Q2": Q2, "Z2": Z2,
                          "big_m": big_m}


def interpolate(x, spec):
    xs, ys = spec["nodes"], spec["values"]
    if x < xs[0] - 1e-6 or x > xs[-1] + 1e-6:
        raise ValueError(f"PWL argument outside public nodes: {x}")
    for j in range(len(xs) - 1):
        if x <= xs[j + 1] + 1e-10:
            t = min(1, max(0, (x - xs[j]) / (xs[j + 1] - xs[j])))
            return ys[j] + t * (ys[j + 1] - ys[j])
    return ys[-1]


def audit(problem, epsilon, m, ids, result):
    x = result.x
    wells = problem["wells"]
    q = {w: float(x[ids["q"][w]]) for w in wells}
    r = {w: float(x[ids["r"][w]]) for w in wells}
    deficits = {w: float(x[ids["d"][w]]) for w in wells}
    Q = (sum(q.values()) - 600) / 20
    R = sum(r.values())
    Z = (R - 600) / 20
    fc = problem["layer2"]["FOPT"]["coefficients"]
    wc = problem["layer2"]["FWIT"]["coefficients"]
    fopt = float(fc["intercept"]) + sum(float(fc[f"r{i}"]) * (r[w] - 80) / 20 for i, w in enumerate(wells, 1))
    z2 = interpolate(Z, problem["pwl_functions"]["Z2"])
    fopt += float(fc["Z2"]) * z2
    fwit = float(wc["intercept"]) + float(wc["Z"]) * Z
    errors = {}
    for w in wells:
        c = problem["layer1_deficit_models"][w]
        own = (q[w] - 80) / 20
        d = float(c.get("intercept", 0)) + float(c.get("own", 0)) * own + float(c.get("Q", 0)) * Q
        if "own2" in c:
            d += float(c["own2"]) * interpolate(own, problem["pwl_functions"]["own2_INJECT3"])
        if "Q2" in c:
            d += float(c["Q2"]) * interpolate(Q, problem["pwl_functions"]["Q2"])
        errors[w] = {"deficit_model": abs(deficits[w] - d),
                     "clipping": abs(r[w] - (q[w] - min(max(d, 0), q[w])))}
    max_error = max(v for e in errors.values() for v in e.values())
    if max_error > 1e-5 or fwit > epsilon + 1e-4 or sum(q.values()) > problem["decision_variables"]["capacity_requested"] + 1e-5:
        raise RuntimeError(f"solution audit failed: {max_error=}, {fwit=}, {epsilon=}")
    return {"status": "optimal", "epsilon_FWIT": epsilon, "requested": q, "realized": r,
            "R": R, "FOPT": fopt, "FWIT": fwit, "solver_gap": float(result.mip_gap),
            "Q": Q, "Z": Z, "Z2": z2, "deficits": deficits, "audit_errors": errors,
            "big_m": ids["big_m"], "solver": "SciPy milp / HiGHS",
            "solver_message": result.message}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    problem = json.loads(PROBLEM.read_text())
    m, objective, ids = build(problem, args.epsilon)
    result = m.solve(objective)
    if result.status != 0 or result.x is None:
        raise RuntimeError(f"MILP failed (status {result.status}): {result.message}")
    output = audit(problem, args.epsilon, m, ids, result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
