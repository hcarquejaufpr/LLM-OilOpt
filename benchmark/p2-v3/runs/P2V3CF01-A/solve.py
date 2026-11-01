#!/usr/bin/env python3
"""MILP for the public P2 v3 water-injection problem.

Only the JSON in shared/ is model data.  scipy.optimize.milp uses HiGHS.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "shared" / "problema_p2_public.json"


class LinearMILP:
    def __init__(self):
        self.names = []
        self.lower = []
        self.upper = []
        self.integer = []
        self.rows = []
        self.rlower = []
        self.rupper = []

    def var(self, name, low=-math.inf, high=math.inf, binary=False):
        k = len(self.names)
        self.names.append(name)
        self.lower.append(0.0 if binary else low)
        self.upper.append(1.0 if binary else high)
        self.integer.append(1 if binary else 0)
        return k

    def row(self, terms, low=-math.inf, high=math.inf):
        self.rows.append({k: v for k, v in terms.items() if v != 0})
        self.rlower.append(low)
        self.rupper.append(high)

    def equal(self, terms, value):
        self.row(terms, value, value)

    def matrices(self):
        rr, cc, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, a in row.items():
                rr.append(i)
                cc.append(j)
                vv.append(a)
        a = coo_matrix((vv, (rr, cc)), shape=(len(self.rows), len(self.names))).tocsr()
        return Bounds(self.lower, self.upper), LinearConstraint(a, self.rlower, self.rupper)


def pwl(m, label, x, nodes, values):
    """Convex combination of only the two vertices of one chosen segment."""
    assert len(nodes) == len(values) and len(nodes) >= 2
    assert all(b > a for a, b in zip(nodes, nodes[1:]))
    lam = [m.var(f"{label}.lambda[{j}]", 0, 1) for j in range(len(nodes))]
    seg = [m.var(f"{label}.segment[{j}]", binary=True) for j in range(len(nodes) - 1)]
    y = m.var(f"{label}.value")
    m.equal({j: 1 for j in lam}, 1)
    m.equal({j: 1 for j in seg}, 1)
    m.equal({x: 1, **{j: -float(v) for j, v in zip(lam, nodes)}}, 0)
    m.equal({y: 1, **{j: -float(v) for j, v in zip(lam, values)}}, 0)
    for j, lj in enumerate(lam):
        row = {lj: 1}
        if j > 0:
            row[seg[j - 1]] = -1
        if j < len(seg):
            row[seg[j]] = -1
        m.row(row, high=0)
    return y


def clipping(m, name, q, d, qlo, qhi, dlo, dhi):
    """Three exact regimes; deficit bounds enter only the Big-M coefficients."""
    # A small outward allowance covers decimal serialization and roundoff.
    L, U = dlo - 1e-7, dhi + 1e-7
    t = m.var(f"{name}.clipped_deficit", 0, qhi)
    r = m.var(f"{name}.realized")
    z0 = m.var(f"{name}.regime_below_zero", binary=True)
    z1 = m.var(f"{name}.regime_inside", binary=True)
    z2 = m.var(f"{name}.regime_above_request", binary=True)
    m.equal({z0: 1, z1: 1, z2: 1}, 1)

    # z0: d <= 0 and t = 0.
    m.row({d: 1, z0: max(0, U)}, high=max(0, U))
    m.row({t: 1, z0: qhi}, high=qhi)

    # z1: 0 <= d <= q and t = d.
    m.row({d: 1, z1: min(0, L)}, low=min(0, L))
    mid_upper = max(0, U - qlo)
    m.row({d: 1, q: -1, z1: mid_upper}, high=mid_upper)
    td_upper = qhi - L
    dt_upper = max(0, U)
    m.row({t: 1, d: -1, z1: td_upper}, high=td_upper)
    m.row({d: 1, t: -1, z1: dt_upper}, high=dt_upper)

    # z2: d >= q and t = q.
    high_lower = max(0, qhi - L)
    m.row({d: 1, q: -1, z2: -high_lower}, low=-high_lower)
    tq_upper = qhi - qlo
    m.row({t: 1, q: -1, z2: tq_upper}, high=tq_upper)
    m.row({q: 1, t: -1, z2: qhi}, high=qhi)

    m.equal({r: 1, q: -1, t: 1}, 0)
    big_m = dict(d_low=L, d_high=U, d_nonpositive=max(0, U),
                 d_nonnegative=max(0, -L), d_at_most_q=mid_upper,
                 d_at_least_q=high_lower, t_zero=qhi,
                 t_minus_d=td_upper, d_minus_t=dt_upper,
                 t_minus_q=tq_upper, q_minus_t=qhi)
    return r, t, (z0, z1, z2), big_m


def build(data):
    m = LinearMILP()
    wells = data["wells"]
    decision = data["decision_variables"]
    q, d, r, t, regimes, big_m = {}, {}, {}, {}, {}, {}
    for w in wells:
        domain = decision["domains"].get(w)
        lo, hi = domain if domain is not None else [decision["fixed"][w]] * 2
        q[w] = m.var(f"{w}.requested", lo, hi)
    m.row({q[w]: 1 for w in wells}, high=decision["capacity_requested"])

    Q = m.var("Q")
    m.equal({Q: 20, **{q[w]: -1 for w in wells}}, -600)
    own = {}
    for w in wells:
        own[w] = m.var(f"{w}.own")
        m.equal({own[w]: 20, q[w]: -1}, -80)

    p = data["pwl_functions"]
    own2 = pwl(m, "own2_INJECT3", own["INJECT3"],
               p["own2_INJECT3"]["nodes"], p["own2_INJECT3"]["values"])
    Q2 = pwl(m, "Q2", Q, p["Q2"]["nodes"], p["Q2"]["values"])
    for w in wells:
        d[w] = m.var(f"{w}.deficit")
        co = data["layer1_deficit_models"][w]
        row = {d[w]: 1}
        for key, idx in (("own", own[w]), ("Q", Q),
                         ("own2", own2), ("Q2", Q2)):
            if key in co:
                row[idx] = -co[key]
        m.equal(row, co.get("intercept", 0))
        domain = decision["domains"].get(w)
        lo, hi = domain if domain is not None else [decision["fixed"][w]] * 2
        r[w], t[w], regimes[w], big_m[w] = clipping(
            m, w, q[w], d[w], lo, hi, *data["deficit_bounds"][w])

    Z = m.var("Z")
    m.equal({Z: 20, **{r[w]: -1 for w in wells}}, -600)
    Z2 = pwl(m, "Z2", Z, p["Z2"]["nodes"], p["Z2"]["values"])
    real = {}
    for w in wells:
        real[w] = m.var(f"{w}.real_feature")
        m.equal({real[w]: 20, r[w]: -1}, -80)

    FOPT = m.var("FOPT")
    co = data["layer2"]["FOPT"]["coefficients"]
    row = {FOPT: 1, Z2: -co["Z2"]}
    for i, w in enumerate(wells, 1):
        row[real[w]] = -co[f"r{i}"]
    m.equal(row, co["intercept"])
    FWIT = m.var("FWIT")
    co = data["layer2"]["FWIT"]["coefficients"]
    m.equal({FWIT: 1, Z: -co["Z"]}, co["intercept"])
    indices = dict(q=q, d=d, r=r, t=t, regimes=regimes, Q=Q, Z=Z,
                   own=own, real=real, own2=own2, Q2=Q2, Z2=Z2,
                   FOPT=FOPT, FWIT=FWIT)
    return m, indices, big_m


def interpolate(x, nodes, values):
    if x < nodes[0] - 1e-6 or x > nodes[-1] + 1e-6:
        raise ValueError(f"PWL argument {x} outside [{nodes[0]}, {nodes[-1]}]")
    return float(np.interp(x, nodes, values))


def validate_auxiliary_bounds(data):
    """Check the supplied deficit intervals by independent interval arithmetic.

    This does not add any row or variable bound to the optimization model.
    """
    wells = data["wells"]
    decision = data["decision_variables"]
    qrange = {w: decision["domains"].get(w, [decision["fixed"].get(w)] * 2)
              for w in wells}
    q_total = [sum(qrange[w][0] for w in wells),
               min(sum(qrange[w][1] for w in wells), decision["capacity_requested"])]
    Qrange = [(x - 600) / 20 for x in q_total]
    p = data["pwl_functions"]
    own3 = [(x - 80) / 20 for x in qrange["INJECT3"]]
    if own3[0] < p["own2_INJECT3"]["nodes"][0] - 1e-10 or \
       own3[1] > p["own2_INJECT3"]["nodes"][-1] + 1e-10:
        raise ValueError("own3 domain is not covered by its PWL grid")
    if Qrange[0] < p["Q2"]["nodes"][0] - 1e-10 or \
       Qrange[1] > p["Q2"]["nodes"][-1] + 1e-10:
        raise ValueError("Q domain is not covered by its PWL grid")

    # A PWL extremum on a closed interval occurs at an endpoint or a node.
    def pwl_range(name, domain):
        nodes, values = p[name]["nodes"], p[name]["values"]
        samples = [interpolate(x, nodes, values) for x in domain]
        samples += [y for x, y in zip(nodes, values) if domain[0] <= x <= domain[1]]
        return min(samples), max(samples)

    common = {"Q": Qrange, "own2": pwl_range("own2_INJECT3", own3),
              "Q2": pwl_range("Q2", Qrange), "intercept": [1, 1]}
    for w in wells:
        ranges = dict(common, own=[(x - 80) / 20 for x in qrange[w]])
        low = high = 0.0
        for key, coefficient in data["layer1_deficit_models"][w].items():
            a, b = ranges[key]
            low += min(coefficient * a, coefficient * b)
            high += max(coefficient * a, coefficient * b)
        supplied_low, supplied_high = data["deficit_bounds"][w]
        if supplied_low > low + 1e-10 or supplied_high < high - 1e-10:
            raise ValueError(f"deficit_bounds for {w} do not contain interval bound")

    # This auxiliary interval proves coverage only; it is never imposed on r.
    z_low = (sum(data["realized_bounds"][w][0] for w in wells) - 600) / 20
    z_high = (min(sum(data["realized_bounds"][w][1] for w in wells),
                  decision["capacity_requested"]) - 600) / 20
    if z_low < p["Z2"]["nodes"][0] - 1e-10 or \
       z_high > p["Z2"]["nodes"][-1] + 1e-10:
        raise ValueError("Z interval is not covered by its PWL grid")


def audit(data, idx, values, epsilon):
    wells = data["wells"]
    get = lambda j: float(values[j])
    q = {w: get(idx["q"][w]) for w in wells}
    r = {w: get(idx["r"][w]) for w in wells}
    d = {w: get(idx["d"][w]) for w in wells}
    p = data["pwl_functions"]
    Q = (sum(q.values()) - 600) / 20
    Z = (sum(r.values()) - 600) / 20
    own = {w: (q[w] - 80) / 20 for w in wells}
    real = {w: (r[w] - 80) / 20 for w in wells}
    own2 = interpolate(own["INJECT3"], **{k: p["own2_INJECT3"][k] for k in ("nodes", "values")})
    Q2 = interpolate(Q, **{k: p["Q2"][k] for k in ("nodes", "values")})
    Z2 = interpolate(Z, **{k: p["Z2"][k] for k in ("nodes", "values")})
    deficit_expected = {}
    for w in wells:
        co = data["layer1_deficit_models"][w]
        deficit_expected[w] = sum(co.get(k, 0) * v for k, v in (
            ("intercept", 1), ("own", own[w]), ("Q", Q),
            ("own2", own2), ("Q2", Q2)))
    co = data["layer2"]["FOPT"]["coefficients"]
    fopt_expected = co["intercept"] + sum(
        co[f"r{i}"] * real[w] for i, w in enumerate(wells, 1)) + co["Z2"] * Z2
    co = data["layer2"]["FWIT"]["coefficients"]
    fwit_expected = co["intercept"] + co["Z"] * Z
    residuals = {
        "deficit_max_abs": max(abs(d[w] - deficit_expected[w]) for w in wells),
        "clipping_max_abs": max(abs(r[w] - (q[w] - min(max(d[w], 0), q[w]))) for w in wells),
        "FOPT_abs": abs(get(idx["FOPT"]) - fopt_expected),
        "FWIT_abs": abs(get(idx["FWIT"]) - fwit_expected),
        "capacity_excess": max(0, sum(q.values()) - data["decision_variables"]["capacity_requested"]),
        "epsilon_excess": max(0, fwit_expected - epsilon),
    }
    if any(v > 1e-4 for v in residuals.values()):
        raise RuntimeError(f"MILP result failed independent formula audit: {residuals}")
    return q, r, d, Q, Z, fopt_expected, fwit_expected, residuals


def solve(epsilon, data, time_limit=None):
    validate_auxiliary_bounds(data)
    m, idx, big_m = build(data)
    m.row({idx["FWIT"]: 1}, high=epsilon)
    bounds, constraints = m.matrices()
    objective = np.zeros(len(m.names))
    objective[idx["FOPT"]] = -1
    options = {"mip_rel_gap": 1e-8}
    if time_limit is not None:
        options["time_limit"] = time_limit
    result = milp(objective, integrality=m.integer, bounds=bounds,
                  constraints=constraints, options=options)
    status = {0: "optimal", 1: "limit_reached", 2: "infeasible",
              3: "unbounded", 4: "solver_error"}.get(result.status, "solver_error")
    output = {"status": status, "epsilon_FWIT": epsilon,
              "requested": None, "realized": None, "R": None,
              "FOPT": None, "FWIT": None, "solver_gap": None,
              "solver_message": result.message}
    if result.x is not None:
        q, r, d, Q, Z, fopt, fwit, residuals = audit(data, idx, result.x, epsilon)
        output.update(requested=q, realized=r, R=sum(r.values()),
                      FOPT=fopt, FWIT=fwit,
                      solver_gap=float(result.mip_gap),
                      Q=Q, Z=Z, deficit=d, audit_residuals=residuals,
                      binary_regimes={w: [float(result.x[k]) for k in idx["regimes"][w]]
                                      for w in data["wells"]})
    return output, m, big_m


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--time-limit", type=float)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("--epsilon must be finite")
    data = json.loads(SOURCE.read_text())
    output, _, _ = solve(args.epsilon, data, args.time_limit)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    if output["status"] != "optimal":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
