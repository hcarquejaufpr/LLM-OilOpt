#!/usr/bin/env python3
"""Solve the public P2 v3 water-injection MILP from its JSON data."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "shared" / "problema_p2_public.json"


class MILP:
    def __init__(self):
        self.names = []
        self.lower = []
        self.upper = []
        self.integer = []
        self.rows = []
        self.row_lower = []
        self.row_upper = []

    def var(self, name, lo=-np.inf, hi=np.inf, binary=False):
        idx = len(self.names)
        self.names.append(name)
        self.lower.append(0 if binary else lo)
        self.upper.append(1 if binary else hi)
        self.integer.append(1 if binary else 0)
        return idx

    def row(self, terms, lo=-np.inf, hi=np.inf):
        combined = {}
        for idx, coef in terms:
            combined[idx] = combined.get(idx, 0.0) + coef
        self.rows.append(combined)
        self.row_lower.append(lo)
        self.row_upper.append(hi)

    def eq(self, terms, rhs=0.0):
        self.row(terms, rhs, rhs)

    def solve(self, objective):
        matrix = lil_matrix((len(self.rows), len(self.names)), dtype=float)
        for j, row in enumerate(self.rows):
            for i, value in row.items():
                matrix[j, i] = value
        c = np.zeros(len(self.names))
        c[objective] = -1.0
        return milp(
            c=c,
            integrality=np.asarray(self.integer),
            bounds=Bounds(self.lower, self.upper),
            constraints=LinearConstraint(matrix.tocsr(), self.row_lower, self.row_upper),
            options={"mip_rel_gap": 1e-9},
        )


def add_pwl(model, name, x, nodes, values):
    if len(nodes) != len(values) or len(nodes) < 2 or any(
        nodes[j] >= nodes[j + 1] for j in range(len(nodes) - 1)
    ):
        raise ValueError(f"Invalid PWL table: {name}")
    y = model.var(f"pwl_{name}")
    lam = [model.var(f"lambda_{name}_{j}", 0, 1) for j in range(len(nodes))]
    seg = [model.var(f"segment_{name}_{j}", binary=True) for j in range(len(nodes) - 1)]
    model.eq([(v, 1) for v in lam], 1)
    model.eq([(v, 1) for v in seg], 1)
    model.eq([(x, 1)] + [(v, -a) for v, a in zip(lam, nodes)])
    model.eq([(y, 1)] + [(v, -a) for v, a in zip(lam, values)])
    for j, v in enumerate(lam):
        adjacent = []
        if j > 0:
            adjacent.append((seg[j - 1], -1))
        if j < len(seg):
            adjacent.append((seg[j], -1))
        model.row([(v, 1)] + adjacent, hi=0)
    return y


def add_clip(model, name, q, r, d, qlo, qhi, dlo, dhi):
    """Exactly encode r=q-clip(d,0,q); auxiliary d bounds enter M only."""
    a, b, c = [model.var(f"regime_{name}_{s}", binary=True) for s in ("low", "mid", "high")]
    model.eq([(a, 1), (b, 1), (c, 1)], 1)

    # Low: d<=0, r=q. Each inactive side uses a valid global expression bound.
    model.row([(d, 1), (a, max(0, dhi))], hi=max(0, dhi))
    model.row([(r, 1), (q, -1), (a, max(0, qhi - qlo))], hi=max(0, qhi - qlo))
    model.row([(r, 1), (q, -1), (a, -qhi)], lo=-qhi)

    # Middle: 0<=d<=q, r=q-d.
    model.row([(d, -1), (b, max(0, -dlo))], hi=max(0, -dlo))
    model.row([(d, 1), (q, -1), (b, max(0, dhi - qlo))], hi=max(0, dhi - qlo))
    mid_upper = max(0, qhi - qlo + dhi)
    mid_lower = max(0, qhi - dlo)
    model.row([(r, 1), (q, -1), (d, 1), (b, mid_upper)], hi=mid_upper)
    model.row([(r, 1), (q, -1), (d, 1), (b, -mid_lower)], lo=-mid_lower)

    # High: d>=q, r=0. This regime may be impossible for some wells.
    model.row([(q, 1), (d, -1), (c, max(0, qhi - dlo))], hi=max(0, qhi - dlo))
    model.row([(r, 1), (c, qhi)], hi=qhi)


def build(data, epsilon):
    wells = data["wells"]
    fixed = data["decision_variables"]["fixed"]
    domains = data["decision_variables"]["domains"]
    m = MILP()
    q, r, d, own, real = {}, {}, {}, {}, {}
    for well in wells:
        lo, hi = (fixed[well], fixed[well]) if well in fixed else domains[well]
        q[well] = m.var(f"q_{well}", lo, hi)
        # 0<=r<=q_hi follows from the clipping equation and q>=0.
        r[well] = m.var(f"r_{well}", 0, hi)
        d[well] = m.var(f"deficit_{well}")
        own[well] = m.var(f"own_{well}")
        real[well] = m.var(f"real_{well}")
        m.eq([(own[well], 20), (q[well], -1)], -80)
        m.eq([(real[well], 20), (r[well], -1)], -80)

    Q, Z = m.var("Q"), m.var("Z")
    m.eq([(Q, 20)] + [(v, -1) for v in q.values()], -600)
    m.eq([(Z, 20)] + [(v, -1) for v in r.values()], -600)
    m.row([(v, 1) for v in q.values()], hi=data["decision_variables"]["capacity_requested"])
    pwl = {}
    for name, table in data["pwl_functions"].items():
        arg = table["argument"]
        x = Q if arg == "Q" else Z if arg == "Z" else own[arg.removeprefix("own_")]
        pwl[name] = add_pwl(m, name, x, table["nodes"], table["values"])

    for well in wells:
        coeff = data["layer1_deficit_models"][well]
        terms = [(d[well], 1)]
        feature = {"Q": Q, "own": own[well], "own2": pwl.get("own2_" + well), "Q2": pwl["Q2"]}
        for key, variable in feature.items():
            if key in coeff:
                terms.append((variable, -coeff[key]))
        m.eq(terms, coeff.get("intercept", 0))
        lo, hi = (fixed[well], fixed[well]) if well in fixed else domains[well]
        dlo, dhi = data["deficit_bounds"][well]
        if dlo > dhi:
            raise ValueError(f"Invalid auxiliary deficit bounds for {well}")
        add_clip(m, well, q[well], r[well], d[well], lo, hi, dlo, dhi)

    fopt, fwit = m.var("FOPT"), m.var("FWIT")
    fc = data["layer2"]["FOPT"]["coefficients"]
    m.eq([(fopt, 1), (pwl["Z2"], -fc["Z2"])] +
         [(real[well], -fc[f"r{j}"]) for j, well in enumerate(wells, 1)], fc["intercept"])
    wc = data["layer2"]["FWIT"]["coefficients"]
    m.eq([(fwit, 1), (Z, -wc["Z"])], wc["intercept"])
    m.row([(fwit, 1)], hi=epsilon)
    return m, fopt, fwit, q, r, d, Q, Z, pwl


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error("--epsilon must be finite")
    data = json.loads(DATA_PATH.read_text())
    m, fopt, fwit, q, r, d, Q, Z, pwl = build(data, args.epsilon)
    result = m.solve(fopt)
    status = {0: "optimal", 1: "limit", 2: "infeasible", 3: "unbounded"}.get(result.status, "error")
    out = {
        "status": status,
        "epsilon_FWIT": args.epsilon,
        "requested": None,
        "realized": None,
        "R": None,
        "FOPT": None,
        "FWIT": None,
        "solver_gap": None if getattr(result, "mip_gap", None) is None else float(result.mip_gap),
        "solver_message": result.message,
    }
    if result.x is not None:
        x = result.x
        out.update({
            "requested": {w: float(x[i]) for w, i in q.items()},
            "realized": {w: float(x[i]) for w, i in r.items()},
            "R": float(sum(x[i] for i in r.values())),
            "FOPT": float(x[fopt]),
            "FWIT": float(x[fwit]),
            "deficit": {w: float(x[i]) for w, i in d.items()},
            "Q": float(x[Q]),
            "Z": float(x[Z]),
            "pwl": {name: float(x[i]) for name, i in pwl.items()},
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    if result.status != 0:
        raise SystemExit(f"Solver status {status}: {result.message}")


if __name__ == "__main__":
    main()
