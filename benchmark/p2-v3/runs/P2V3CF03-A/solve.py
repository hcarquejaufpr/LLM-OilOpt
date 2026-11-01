#!/usr/bin/env python3
"""Public P2 v3 MILP, using SciPy/HiGHS and the supplied JSON coefficients."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "shared" / "problema_p2_public.json"


class MILP:
    def __init__(self):
        self.names, self.lb, self.ub, self.integer = [], [], [], []
        self.rows, self.lo, self.hi = [], [], []

    def var(self, name, lo=-np.inf, hi=np.inf, binary=False):
        j = len(self.names)
        self.names.append(name)
        self.lb.append(0.0 if binary else lo)
        self.ub.append(1.0 if binary else hi)
        self.integer.append(1 if binary else 0)
        return j

    def row(self, terms, lo=-np.inf, hi=np.inf):
        merged = {}
        for j, a in terms:
            merged[j] = merged.get(j, 0.0) + a
        self.rows.append({j: a for j, a in merged.items() if a})
        self.lo.append(lo)
        self.hi.append(hi)

    def eq(self, terms, rhs):
        self.row(terms, rhs, rhs)

    def le(self, terms, rhs):
        self.row(terms, hi=rhs)

    def ge(self, terms, rhs):
        self.row(terms, lo=rhs)

    def matrices(self):
        ii, jj, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, a in row.items():
                ii.append(i)
                jj.append(j)
                vv.append(a)
        a = coo_matrix((vv, (ii, jj)), shape=(len(self.rows), len(self.names))).tocsr()
        return Bounds(self.lb, self.ub), LinearConstraint(a, self.lo, self.hi)


def pwl(m, label, x, spec):
    """Exactly one segment; two nonnegative weights on its adjacent endpoints."""
    nodes, vals = spec["nodes"], spec["values"]
    if len(nodes) != len(vals) or len(nodes) < 2 or any(b <= a for a, b in zip(nodes, nodes[1:])):
        raise ValueError(f"Invalid PWL nodes for {label}")
    value = m.var(label + ".value")
    selections, weights = [], []
    for k in range(len(nodes) - 1):
        y = m.var(f"{label}.segment{k}", binary=True)
        a = m.var(f"{label}.left{k}", lo=0)
        b = m.var(f"{label}.right{k}", lo=0)
        m.eq([(a, 1), (b, 1), (y, -1)], 0)
        selections.append(y)
        weights.append((a, b))
    m.eq([(y, 1) for y in selections], 1)
    m.eq([(x, 1)] + [(v, -nodes[k + side]) for k, pair in enumerate(weights)
                       for side, v in enumerate(pair)], 0)
    m.eq([(value, 1)] + [(v, -vals[k + side]) for k, pair in enumerate(weights)
                           for side, v in enumerate(pair)], 0)
    return value


def clip(m, label, q, d, r, qlo, qhi, dlo, dhi):
    """Exact disjunction. Auxiliary deficit bounds are used only for valid M."""
    # Enlarge rounded auxiliary bounds for floating-point and JSON rounding.
    margin = 1e-7
    L, U = dlo - margin, dhi + margin
    low, middle, high = [m.var(f"{label}.clip.{s}", binary=True)
                         for s in ("low", "middle", "high")]
    m.eq([(low, 1), (middle, 1), (high, 1)], 1)
    # 0 <= r <= q follows directly from the clipping identity.
    m.le([(r, 1), (q, -1)], 0)
    # Low: d <= 0 and r = q.
    m.le([(d, 1), (low, U)], U)
    m.ge([(r, 1), (q, -1), (low, -qhi)], -qhi)
    # Middle: 0 <= d <= q and r = q-d.
    m.ge([(d, 1), (middle, L)], L)
    m.le([(d, 1), (q, -1), (middle, max(0.0, U - qlo))], max(0.0, U - qlo))
    m.ge([(r, 1), (q, -1), (d, 1), (middle, -(qhi - L))], -(qhi - L))
    m.le([(r, 1), (q, -1), (d, 1), (middle, U)], U)
    # High: d >= q and r = 0.
    m.ge([(d, 1), (q, -1), (high, -(qhi - L))], -(qhi - L))
    m.le([(r, 1), (high, qhi)], qhi)


def build(data, epsilon):
    m = MILP()
    wells = data["wells"]
    q, r, d = {}, {}, {}
    domains = data["decision_variables"]
    for w in wells:
        lo, hi = ([domains["fixed"][w]] * 2 if w in domains["fixed"]
                  else domains["domains"][w])
        q[w] = m.var("q." + w, lo, hi)
        # 0 <= r <= q is a direct consequence of clipping, independent of auxiliary bounds.
        r[w] = m.var("r." + w, 0, hi)
        d[w] = m.var("deficit." + w)
    m.le([(q[w], 1) for w in wells], domains["capacity_requested"])
    Q = m.var("Q")
    Z = m.var("Z")
    m.eq([(Q, 20)] + [(q[w], -1) for w in wells], -600)
    m.eq([(Z, 20)] + [(r[w], -1) for w in wells], -600)
    own = {}
    for w in wells:
        own[w] = m.var("own." + w)
        m.eq([(own[w], 20), (q[w], -1)], -80)
    pw = data["pwl_functions"]
    own2 = pwl(m, "own2_INJECT3", own["INJECT3"], pw["own2_INJECT3"])
    Q2 = pwl(m, "Q2", Q, pw["Q2"])
    Z2 = pwl(m, "Z2", Z, pw["Z2"])
    for w in wells:
        c = data["layer1_deficit_models"][w]
        terms = [(d[w], 1), (Q, -c.get("Q", 0)),
                 (own[w], -c.get("own", 0))]
        if w == "INJECT3":
            terms += [(own2, -c.get("own2", 0)), (Q2, -c.get("Q2", 0))]
        m.eq(terms, c.get("intercept", 0))
        lo, hi = ([domains["fixed"][w]] * 2 if w in domains["fixed"]
                  else domains["domains"][w])
        clip(m, w, q[w], d[w], r[w], lo, hi, *data["deficit_bounds"][w])
    FOPT = m.var("FOPT")
    FWIT = m.var("FWIT")
    fc = data["layer2"]["FOPT"]["coefficients"]
    # Each real_i is (r_i-80)/20. Shift the intercept exactly as this implies.
    m.eq([(FOPT, 1), (Z2, -fc["Z2"])] +
         [(r[w], -fc[f"r{k}"] / 20) for k, w in enumerate(wells, 1)],
         fc["intercept"] - 4 * sum(fc[f"r{k}"] for k in range(1, 9)))
    wc = data["layer2"]["FWIT"]["coefficients"]
    m.eq([(FWIT, 1), (Z, -wc["Z"])], wc["intercept"])
    m.le([(FWIT, 1)], epsilon)
    return m, q, r, d, FOPT, FWIT


def calculate(data, requested):
    """Independent forward evaluation of the published algebra and PWL curves."""
    wells = data["wells"]
    Q = (sum(requested.values()) - 600) / 20
    def interp(name, x):
        p = data["pwl_functions"][name]
        nodes, values = p["nodes"], p["values"]
        if x < nodes[0] - 1e-7 or x > nodes[-1] + 1e-7:
            raise ValueError(f"{name} argument outside public nodes: {x}")
        k = max(0, min(len(nodes)-2, int(np.searchsorted(nodes, x, side="right")-1)))
        t = (x-nodes[k])/(nodes[k+1]-nodes[k])
        return values[k]*(1-t)+values[k+1]*t
    deficit, realized = {}, {}
    for w in wells:
        c = data["layer1_deficit_models"][w]
        own = (requested[w]-80)/20
        x = c.get("intercept", 0) + c.get("Q", 0)*Q + c.get("own", 0)*own
        if w == "INJECT3":
            x += c.get("own2", 0)*interp("own2_INJECT3", own)
            x += c.get("Q2", 0)*interp("Q2", Q)
        deficit[w] = x
        realized[w] = requested[w] - min(max(x, 0), requested[w])
    Z = (sum(realized.values())-600)/20
    fc = data["layer2"]["FOPT"]["coefficients"]
    FOPT = fc["intercept"] + sum(fc[f"r{k}"]*(realized[w]-80)/20
                                    for k, w in enumerate(wells, 1)) + fc["Z2"]*interp("Z2", Z)
    wc = data["layer2"]["FWIT"]["coefficients"]
    FWIT = wc["intercept"] + wc["Z"]*Z
    return deficit, realized, FOPT, FWIT


def solve(epsilon):
    data = json.loads(DATA.read_text())
    m, q, r, d, f, fw = build(data, epsilon)
    bounds, constraint = m.matrices()
    obj = np.zeros(len(m.names))
    obj[f] = -1
    result = milp(obj, integrality=m.integer, bounds=bounds,
                  constraints=constraint, options={"mip_rel_gap": 1e-8})
    output = {"status": "optimal" if result.status == 0 else
              ("limit" if result.status == 1 else "infeasible" if result.status == 2 else "error"),
              "epsilon_FWIT": epsilon,
              "solver_gap": float(result.mip_gap) if getattr(result, "mip_gap", None) is not None else None,
              "solver_message": result.message}
    if result.x is not None:
        # HiGHS may return a bound with a roughly 1e-7 feasibility residual.
        # Snap only such numerical residuals to the published domain endpoints.
        domains = data["decision_variables"]
        requested = {}
        for w, j in q.items():
            lo, hi = ([domains["fixed"][w]] * 2 if w in domains["fixed"]
                      else domains["domains"][w])
            value = float(result.x[j])
            if abs(value - lo) <= 2e-7:
                value = lo
            elif abs(value - hi) <= 2e-7:
                value = hi
            requested[w] = value
        if sum(requested.values()) > domains["capacity_requested"]:
            excess = sum(requested.values()) - domains["capacity_requested"] + 1e-7
            movable = [w for w in q if w not in domains["fixed"]
                       and requested[w] - domains["domains"][w][0] >= excess]
            if not movable:
                raise RuntimeError("Cannot remove numerical capacity residual")
            requested[max(movable, key=lambda w: requested[w]-domains["domains"][w][0])] -= excess
        deficit, realized, FOPT, FWIT = calculate(data, requested)
        if FWIT > epsilon:
            # Re-evaluate the published equations after boundary snapping, then
            # remove only the resulting numerical budget residual.
            candidates = []
            for w in q:
                if w in domains["fixed"]:
                    continue
                lo = domains["domains"][w][0]
                if requested[w] <= lo + 1e-6:
                    continue
                trial = dict(requested)
                trial[w] = lo
                if calculate(data, trial)[3] <= epsilon - 1e-5:
                    candidates.append(w)
            if not candidates:
                raise RuntimeError("Cannot remove numerical FWIT residual")
            w = max(candidates, key=lambda w: requested[w]-domains["domains"][w][0])
            left, right = domains["domains"][w][0], requested[w]
            for _ in range(50):
                mid = (left+right)/2
                trial = dict(requested)
                trial[w] = mid
                if calculate(data, trial)[3] <= epsilon - 1e-5:
                    left = mid
                else:
                    right = mid
            requested[w] = left
            deficit, realized, FOPT, FWIT = calculate(data, requested)
        output.update(requested=requested, realized=realized,
                      R=float(sum(realized.values())), FOPT=float(FOPT), FWIT=float(FWIT),
                      deficit=deficit,
                      max_realized_residual=max(abs(realized[w]-result.x[r[w]]) for w in q),
                      FOPT_milp=float(result.x[f]), FWIT_milp=float(result.x[fw]),
                      requested_total=float(sum(requested.values())))
        if getattr(result, "mip_dual_bound", None) is not None:
            # A feasible numerical polish may move the reported objective by
            # a few micro-units; include that move in the reported gap.
            upper_bound = -float(result.mip_dual_bound)
            output["solver_gap"] = max(output["solver_gap"] or 0.0,
                                       (upper_bound - FOPT) / max(1.0, abs(FOPT)))
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("epsilon must be finite")
    output = solve(args.epsilon)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    if output["status"] == "error":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
