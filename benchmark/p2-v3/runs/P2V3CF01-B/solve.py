#!/usr/bin/env python3
"""Exact public P2 v3 MILP, using SciPy/HiGHS."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "shared" / "problema_p2_public.json"


class Builder:
    def __init__(self):
        self.names, self.lower, self.upper, self.integer = [], [], [], []
        self.rows, self.cols, self.data, self.cl, self.cu = [], [], [], [], []

    def var(self, name, lo=-np.inf, hi=np.inf, integer=False):
        i = len(self.names)
        self.names.append(name)
        self.lower.append(lo)
        self.upper.append(hi)
        self.integer.append(int(integer))
        return i

    def row(self, terms, lo=-np.inf, hi=np.inf):
        k = len(self.cl)
        for col, value in terms.items():
            if value:
                self.rows.append(k)
                self.cols.append(col)
                self.data.append(value)
        self.cl.append(lo)
        self.cu.append(hi)

    def equal(self, terms, value=0.0):
        self.row(terms, value, value)

    def pwl(self, name, x, nodes, values):
        assert len(nodes) == len(values) and len(nodes) >= 2
        assert all(a < b for a, b in zip(nodes, nodes[1:]))
        y = self.var(name)
        lam = [self.var(f"{name}.lambda.{j}", 0, 1) for j in range(len(nodes))]
        seg = [self.var(f"{name}.segment.{j}", 0, 1, True) for j in range(len(nodes)-1)]
        self.equal({v: 1 for v in lam}, 1)
        self.equal({v: 1 for v in seg}, 1)
        self.equal({x: 1, **{v: -node for v, node in zip(lam, nodes)}})
        self.equal({y: 1, **{v: -value for v, value in zip(lam, values)}})
        for j, v in enumerate(lam):
            terms = {v: 1}
            if j > 0:
                terms[seg[j-1]] = -1
            if j < len(seg):
                terms[seg[j]] = -1
            self.row(terms, hi=0)
        return y

    def matrices(self):
        a = coo_matrix((self.data, (self.rows, self.cols)),
                       shape=(len(self.cl), len(self.names))).tocsr()
        return Bounds(self.lower, self.upper), LinearConstraint(a, self.cl, self.cu)


def interval_bounds(problem):
    """Check that published deficit bounds cover box interval arithmetic."""
    wells = problem["wells"]
    domains = {**{w: [v, v] for w, v in problem["decision_variables"]["fixed"].items()},
               **problem["decision_variables"]["domains"]}
    qlo = sum(domains[w][0] for w in wells)
    qhi = min(sum(domains[w][1] for w in wells),
              problem["decision_variables"]["capacity_requested"])
    q_interval = ((qlo-600)/20, (qhi-600)/20)
    def add_term(box, coefficient, limits):
        vals = [coefficient * limits[0], coefficient * limits[1]]
        return box[0] + min(vals), box[1] + max(vals)
    checked = {}
    for w in wells:
        model = problem["layer1_deficit_models"][w]
        box = (model.get("intercept", 0.0),) * 2
        own = tuple((v-80)/20 for v in domains[w])
        for key, limits in (("own", own), ("Q", q_interval)):
            box = add_term(box, model.get(key, 0.0), limits)
        if w == "INJECT3":
            for key, pwl_name in (("own2", "own2_INJECT3"), ("Q2", "Q2")):
                values = problem["pwl_functions"][pwl_name]["values"]
                box = add_term(box, model[key], (min(values), max(values)))
        supplied = problem["deficit_bounds"][w]
        if supplied[0] > box[0] + 1e-8 or supplied[1] < box[1] - 1e-8:
            raise ValueError(f"Auxiliary deficit bounds do not cover interval for {w}: {box}")
        checked[w] = {"interval_arithmetic": list(box), "published": supplied}
    return domains, q_interval, checked


def big_m(problem, domains):
    result = {}
    for w in problem["wells"]:
        dl, du = problem["deficit_bounds"][w]
        ql, qu = domains[w]
        result[w] = {
            "d_le_0": max(0.0, du),
            "d_ge_0": max(0.0, -dl),
            "d_le_q": max(0.0, du-ql),
            "d_ge_q": max(0.0, qu-dl),
            "r_eq_q": max(0.0, du),
            "r_eq_q_minus_d": max(0.0, -dl, du-ql),
            "r_eq_0": qu,
        }
    return result


def activated_equality(b, terms, binary, m):
    # -M(1-b) <= expression <= M(1-b)
    b.row({**terms, binary: m}, hi=m)
    b.row({**terms, binary: -m}, lo=-m)


def build(problem, epsilon):
    domains, q_interval, checked = interval_bounds(problem)
    ms = big_m(problem, domains)
    b = Builder()
    wells = problem["wells"]
    q = {w: b.var(f"q.{w}", *domains[w]) for w in wells}
    r = {w: b.var(f"r.{w}") for w in wells}
    d = {w: b.var(f"deficit.{w}") for w in wells}
    own = {w: b.var(f"own.{w}") for w in wells}
    real = {w: b.var(f"real.{w}") for w in wells}
    Q, R, Z = (b.var(v) for v in ("Q", "R", "Z"))
    b.row({v: 1 for v in q.values()}, hi=problem["decision_variables"]["capacity_requested"])
    b.equal({Q: 1, **{v: -1/20 for v in q.values()}}, -30)
    b.equal({R: 1, **{v: -1 for v in r.values()}})
    b.equal({Z: 1, R: -1/20}, -30)
    for w in wells:
        b.equal({own[w]: 1, q[w]: -1/20}, -4)
        b.equal({real[w]: 1, r[w]: -1/20}, -4)
    pwls = problem["pwl_functions"]
    own2 = b.pwl("own2_INJECT3", own["INJECT3"], **{k: pwls["own2_INJECT3"][k] for k in ("nodes", "values")})
    Q2 = b.pwl("Q2", Q, **{k: pwls["Q2"][k] for k in ("nodes", "values")})
    Z2 = b.pwl("Z2", Z, **{k: pwls["Z2"][k] for k in ("nodes", "values")})
    for w in wells:
        model = problem["layer1_deficit_models"][w]
        terms = {d[w]: 1, own[w]: -model.get("own", 0), Q: -model.get("Q", 0)}
        if w == "INJECT3":
            terms[own2] = -model["own2"]
            terms[Q2] = -model["Q2"]
        b.equal(terms, model.get("intercept", 0))
        low = b.var(f"regime.low.{w}", 0, 1, True)
        mid = b.var(f"regime.mid.{w}", 0, 1, True)
        high = b.var(f"regime.high.{w}", 0, 1, True)
        b.equal({low: 1, mid: 1, high: 1}, 1)
        m = ms[w]
        # low: d<=0, r=q
        b.row({d[w]: 1, low: m["d_le_0"]}, hi=m["d_le_0"])
        activated_equality(b, {r[w]: 1, q[w]: -1}, low, m["r_eq_q"])
        # mid: 0<=d<=q, r=q-d
        b.row({d[w]: -1, mid: m["d_ge_0"]}, hi=m["d_ge_0"])
        b.row({d[w]: 1, q[w]: -1, mid: m["d_le_q"]}, hi=m["d_le_q"])
        activated_equality(b, {r[w]: 1, q[w]: -1, d[w]: 1}, mid, m["r_eq_q_minus_d"])
        # high: d>=q, r=0
        b.row({q[w]: 1, d[w]: -1, high: m["d_ge_q"]}, hi=m["d_ge_q"])
        activated_equality(b, {r[w]: 1}, high, m["r_eq_0"])
    fc = problem["layer2"]["FOPT"]["coefficients"]
    wc = problem["layer2"]["FWIT"]["coefficients"]
    fopt, fwit = b.var("FOPT"), b.var("FWIT")
    b.equal({fopt: 1, **{real[w]: -fc[f"r{i}"] for i, w in enumerate(wells, 1)},
             Z2: -fc["Z2"]}, fc["intercept"])
    b.equal({fwit: 1, Z: -wc["Z"]}, wc["intercept"])
    b.row({fwit: 1}, hi=epsilon)
    indices = {"q": q, "r": r, "d": d, "own": own, "real": real,
               "Q": Q, "R": R, "Z": Z, "own2_INJECT3": own2,
               "Q2": Q2, "Z2": Z2, "FOPT": fopt, "FWIT": fwit}
    return b, indices, {"big_m": ms, "deficit_interval_validation": checked,
                        "Q_auxiliary_interval": q_interval}


def interpolate(x, spec):
    nodes, values = spec["nodes"], spec["values"]
    if x < nodes[0]-1e-6 or x > nodes[-1]+1e-6:
        raise ValueError(f"PWL argument outside public nodes: {x}")
    j = min(max(int(np.searchsorted(nodes, x, side="right"))-1, 0), len(nodes)-2)
    t = (x-nodes[j])/(nodes[j+1]-nodes[j])
    return values[j]+t*(values[j+1]-values[j])


def audit(problem, epsilon, x, idx, raw):
    wells = problem["wells"]
    q = {w: float(x[idx["q"][w]]) for w in wells}
    Q = (sum(q.values())-600)/20
    pwl = problem["pwl_functions"]
    deficits = {}
    r = {}
    for w in wells:
        m = problem["layer1_deficit_models"][w]
        own = (q[w]-80)/20
        deficit = m.get("intercept", 0)+m.get("own", 0)*own+m.get("Q", 0)*Q
        if w == "INJECT3":
            deficit += m["own2"]*interpolate(own, pwl["own2_INJECT3"])
            deficit += m["Q2"]*interpolate(Q, pwl["Q2"])
        deficits[w] = deficit
        r[w] = q[w]-min(max(deficit, 0), q[w])
    R = sum(r.values())
    Z = (R-600)/20
    Z2 = interpolate(Z, pwl["Z2"])
    fc = problem["layer2"]["FOPT"]["coefficients"]
    wc = problem["layer2"]["FWIT"]["coefficients"]
    fopt = fc["intercept"]+sum(fc[f"r{i}"]*(r[w]-80)/20 for i, w in enumerate(wells, 1))+fc["Z2"]*Z2
    fwit = wc["intercept"]+wc["Z"]*Z
    max_r_diff = max(abs(r[w]-x[idx["r"][w]]) for w in wells)
    if max_r_diff > 1e-4 or abs(fopt-x[idx["FOPT"]]) > 1e-3 or fwit > epsilon+1e-3:
        raise RuntimeError("MILP solution failed independent public-equation audit")
    return {"status": "optimal" if raw.status == 0 else "feasible_limit",
            "epsilon_FWIT": epsilon, "requested": q, "realized": r,
            "R": R, "FOPT": fopt, "FWIT": fwit,
            "solver_gap": float(raw.mip_gap),
            "Q": Q, "Z": Z, "deficits": deficits,
            "features": {"own": {w: (q[w]-80)/20 for w in wells},
                         "real": {w: (r[w]-80)/20 for w in wells},
                         "own2_INJECT3": interpolate((q["INJECT3"]-80)/20, pwl["own2_INJECT3"]),
                         "Q2": interpolate(Q, pwl["Q2"]), "Z2": Z2},
            "audit": {"max_realized_milp_difference": max_r_diff,
                      "FOPT_milp_difference": abs(fopt-x[idx["FOPT"]]),
                      "FWIT_milp_difference": abs(fwit-x[idx["FWIT"]])},
            "solver_message": raw.message}


def solve(epsilon):
    problem = json.loads(SOURCE.read_text())
    b, idx, metadata = build(problem, epsilon)
    bounds, constraint = b.matrices()
    c = np.zeros(len(b.names))
    c[idx["FOPT"]] = -1
    result = milp(c, integrality=b.integer, bounds=bounds, constraints=constraint,
                  options={"mip_rel_gap": 1e-8, "time_limit": 100.0})
    if result.x is None or result.status not in (0, 1):
        raise RuntimeError(f"HiGHS failed: {result.message}")
    solution = audit(problem, epsilon, result.x, idx, result)
    solution["model_sha256"] = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    solution["milp_size"] = {"variables": len(b.names), "binary": sum(b.integer),
                             "constraints": len(b.cl)}
    return solution


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    solution = solve(args.epsilon)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(solution, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
