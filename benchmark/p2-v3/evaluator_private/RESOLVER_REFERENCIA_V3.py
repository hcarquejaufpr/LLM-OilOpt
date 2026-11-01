#!/usr/bin/env python3
"""Solver de referencia P2 v3 (MILP/SOS2 exato via scipy.optimize.milp -> HiGHS).

Le somente shared/problema_p2_public.json (v3). Nao usa OPM, sessoes LLM, v1/v2 nem valores P1-P5 pre-conhecidos.

Modos:
  bigM_only               deficit_i LIVRE; deficit_bounds/realized_bounds apenas derivam Big-M (formulacao correta, v3).
  loose_bigM              deficit_i livre com Big-M largo e independente dos bounds (controle cruzado).
  impose_deficit_bounds   INCORRETO: impoe bounds do arquivo recebido como limites de deficit_i (reproduz a ambiguidade v2).
O modo e escolhido pelo chamador; a referencia oficial usa bigM_only e e verificada contra loose_bigM e a
avaliacao exata (sem MILP).
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

HERE = Path(__file__).resolve().parent
PUBLIC = HERE.parent / "shared" / "problema_p2_public.json"
MIP_REL_GAP = 1e-9
TIME_LIMIT = 900


def load(path=PUBLIC):
    return json.loads(Path(path).read_text())


def domains(P):
    W = P["wells"]
    dom = {w: tuple(P["decision_variables"]["domains"][w]) for w in W if w != "INJECT1"}
    dom["INJECT1"] = (79.5, 79.5)
    return W, dom


def pw(f, x, tol=1e-7):
    assert f["nodes"][0] - tol <= x <= f["nodes"][-1] + tol, ("extrapolacao PWL", x)
    return float(np.interp(x, f["nodes"], f["values"]))


def evaluate(P, q):
    """Avaliacao exata (sem MILP) do modelo publico para um vetor q (ordem de P['wells'])."""
    W, _ = domains(P)
    L, PW = P["layer1_deficit_models"], P["pwl_functions"]
    co, fw = P["layer2"]["FOPT"]["coefficients"], P["layer2"]["FWIT"]["coefficients"]
    q = np.asarray(q, float)
    Q = (q.sum() - 600) / 20
    d = np.zeros(8)
    for i, w in enumerate(W):
        m = L[w]
        if m:
            o = (q[i] - 80) / 20
            d[i] = m.get("intercept", 0) + m.get("own", 0) * o + m.get("Q", 0) * Q
            if w == "INJECT3":
                d[i] += m["own2"] * pw(PW["own2_INJECT3"], o) + m["Q2"] * pw(PW["Q2"], Q)
    r = q - np.clip(d, 0, q)
    Z = (r.sum() - 600) / 20
    fopt = co["intercept"] + sum(co["r%d" % (i + 1)] * (r[i] - 80) / 20 for i in range(8)) + co["Z2"] * pw(PW["Z2"], Z)
    fwit = fw["intercept"] + fw["Z"] * Z
    return {"FOPT": float(fopt), "FWIT": float(fwit), "r": r, "deficit": d, "R": float(r.sum()), "Q": float(Q), "Z": float(Z)}


class Model:
    def __init__(self):
        self.lb, self.ub, self.ig, self.rows = [], [], [], []

    def var(self, lo, hi, integer=0):
        self.lb.append(lo)
        self.ub.append(hi)
        self.ig.append(integer)
        return len(self.lb) - 1

    @staticmethod
    def expr(*terms):
        out = {}
        for c, x in terms:
            for k, v in x.items():
                out[k] = out.get(k, 0.0) + c * v
        return out

    def le(self, a, b):
        x = self.expr((1, a), (-1, b))
        c = x.pop(-1, 0.0)
        self.rows.append((x, -np.inf, -c))

    def ge(self, a, b):
        self.le(b, a)

    def eq(self, a, b):
        x = self.expr((1, a), (-1, b))
        c = x.pop(-1, 0.0)
        self.rows.append((x, -c, -c))


V = lambda i: {i: 1.0}
C = lambda c: {-1: c}


def build_and_solve(P, epsilon, mode="bigM_only", imposed_bounds=None, log=None):
    W, dom = domains(P)
    L, PW = P["layer1_deficit_models"], P["pwl_functions"]
    co, fw = P["layer2"]["FOPT"]["coefficients"], P["layer2"]["FWIT"]["coefficients"]
    M = Model()
    E = M.expr

    def pwl(x, f):
        n = len(f["nodes"])
        lam = [M.var(0, 1) for _ in range(n)]
        seg = [M.var(0, 1, 1) for _ in range(n - 1)]
        M.eq(E(*[(1, V(l)) for l in lam]), C(1))
        M.eq(E(*[(1, V(s)) for s in seg]), C(1))
        for k in range(n):
            adj = [seg[j] for j in (k - 1, k) if 0 <= j < n - 1]
            M.le(V(lam[k]), E(*[(1, V(j)) for j in adj]))
        M.eq(x, E(*[(f["nodes"][k], V(lam[k])) for k in range(n)]))
        return E(*[(f["values"][k], V(lam[k])) for k in range(n)])

    lo = [dom[w][0] for w in W]
    hi = [dom[w][1] for w in W]
    q = [C(79.5) if i == 0 else V(M.var(lo[i], hi[i])) for i in range(8)]
    S = E(*[(1, x) for x in q])
    Qx = E((1 / 20, S), (1, C(-30)))
    r = []
    for i, w in enumerate(W):
        m = L[w]
        if not m:
            r.append(q[i])
            continue
        own = E((1 / 20, q[i]), (1, C(-4)))
        de = E((1, C(m.get("intercept", 0))), (m.get("own", 0), own), (m.get("Q", 0), Qx))
        if w == "INJECT3":
            de = E((1, de), (m["own2"], pwl(own, PW["own2_INJECT3"])), (m["Q2"], pwl(Qx, PW["Q2"])))
        if mode == "impose_deficit_bounds":
            dl, dh = imposed_bounds[w]
            d = M.var(dl, dh)
        else:
            d = M.var(-1e3, 1e3)
        M.eq(V(d), de)
        if mode == "loose_bigM":
            dl, dh, BM = -50.0, 50.0, 400.0
        else:
            src = P["deficit_bounds"] if mode == "bigM_only" else imposed_bounds
            dl, dh = src[w]
            BM = hi[i] + max(abs(dl), abs(dh)) + 1.0
        ri = M.var(0, hi[i])
        A, B, Cc = (M.var(0, 1, 1) for _ in range(3))
        M.eq(E((1, V(A)), (1, V(B)), (1, V(Cc))), C(1))
        # regime A: d<=0, r=q | B: 0<=d<=q, r=q-d | C: d>=q, r=0
        M.le(V(d), E((dh, C(1)), (-dh, V(A))))
        M.ge(V(ri), E((1, q[i]), (-BM, C(1)), (BM, V(A))))
        M.ge(V(d), E((dl, C(1)), (-dl, V(B))))
        M2 = dh - lo[i]
        M.le(E((1, V(d)), (-1, q[i])), E((M2, C(1)), (-M2, V(B))))
        M.ge(V(ri), E((1, q[i]), (-1, V(d)), (-BM, C(1)), (BM, V(B))))
        M.le(V(ri), E((1, q[i]), (-1, V(d)), (BM, C(1)), (-BM, V(B))))
        M3 = dl - hi[i]
        M.ge(V(d), E((1, q[i]), (M3, C(1)), (-M3, V(Cc))))
        M.le(V(ri), E((BM, C(1)), (-BM, V(Cc))))
        M.le(V(ri), q[i])
        r.append(V(ri))
    R = E(*[(1, x) for x in r])
    Zx = E((1 / 20, R), (1, C(-30)))
    z2 = pwl(Zx, PW["Z2"])
    M.le(S, C(P["decision_variables"]["capacity_requested"]))
    M.le(Zx, C((epsilon - fw["intercept"]) / fw["Z"]))
    obj = E(*[(co["r%d" % (i + 1)] / 20, r[i]) for i in range(8)], (co["Z2"], z2))
    const = co["intercept"] - sum(co["r%d" % (i + 1)] * 4 for i in range(8))
    n = len(M.lb)
    c = np.zeros(n)
    for k, v in obj.items():
        if k >= 0:
            c[k] = -v
    A_ = lil_matrix((len(M.rows), n))
    lo_, up_ = [], []
    for j, (x, a, b) in enumerate(M.rows):
        for k, v in x.items():
            A_[j, k] = v
        lo_.append(a)
        up_.append(b)
    t0 = time.time()
    res = milp(c, constraints=LinearConstraint(A_.tocsr(), lo_, up_), integrality=M.ig, bounds=Bounds(M.lb, M.ub),
               options={"mip_rel_gap": MIP_REL_GAP, "time_limit": TIME_LIMIT, "presolve": True, "disp": False})
    out = {"mode": mode, "status_code": int(res.status), "message": res.message, "success": bool(res.success),
           "mip_gap": getattr(res, "mip_gap", None), "dual_bound": getattr(res, "mip_dual_bound", None),
           "time_s": round(time.time() - t0, 3), "n_vars": n, "n_integer": int(sum(M.ig)), "n_rows": len(M.rows)}
    if res.x is not None:
        qv = np.array([79.5] + list(res.x[:7]))
        out["q"] = qv.tolist()
        out["FOPT_milp"] = float(-res.fun + const)
    return out


def reference_for(P, epsilon, mode="bigM_only", imposed_bounds=None):
    o = build_and_solve(P, epsilon, mode, imposed_bounds)
    if "q" in o:
        ev = evaluate(P, o["q"])
        o.update({"requested": o["q"], "realized": ev["r"].tolist(), "R": ev["R"], "FOPT": ev["FOPT"], "FWIT": ev["FWIT"],
                  "deficit": ev["deficit"].tolist(), "FWIT_slack": float(epsilon - ev["FWIT"]),
                  "capacity_slack": float(P["decision_variables"]["capacity_requested"] - sum(o["q"])),
                  "milp_vs_exact_FOPT_abs": abs(o["FOPT_milp"] - ev["FOPT"])})
    return o


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=str(HERE / "REFERENCIA_P2_PRIVADA_V3.json"))
    a = ap.parse_args()
    P = load()
    pol = {}
    for e in P["optimization"]["epsilon_candidates"]:
        main_ = reference_for(P, e["epsilon_FWIT"], "bigM_only")
        ctrl = reference_for(P, e["epsilon_FWIT"], "loose_bigM")
        main_["control_loose_bigM"] = {k: ctrl.get(k) for k in ("status_code", "mip_gap", "FOPT", "FWIT", "requested", "time_s")}
        main_["control_FOPT_diff_abs"] = abs(main_["FOPT"] - ctrl["FOPT"])
        main_["epsilon_FWIT"] = e["epsilon_FWIT"]
        pol[e["policy"]] = main_
    ref = {
        "schema": "p2-private-reference-v3",
        "formulation": "exact 3-regime clipping, SOS2 via lambda+adjacent binaries; deficit_i free",
        "deficit_bounds_imposed_as_constraints": False,
        "deficit_bounds_used_only_for": "Big-M constants",
        "solver": {"name": "HiGHS via scipy.optimize.milp", "scipy": scipy.__version__, "numpy": np.__version__,
                   "python": sys.version.split()[0], "platform": platform.platform(),
                   "mip_rel_gap_requested": MIP_REL_GAP, "time_limit_s": TIME_LIMIT, "presolve": True},
        "inputs_sha256": {"shared/problema_p2_public.json": sha256(PUBLIC),
                          "evaluator_private/RESOLVER_REFERENCIA_V3.py": sha256(__file__)},
        "policies": pol,
    }
    Path(a.output).write_text(json.dumps(ref, indent=2) + "\n")
    for k, v in pol.items():
        print(k, v["status_code"], v["mip_gap"], round(v["FOPT"], 6), round(v["FWIT"], 3), round(v["control_FOPT_diff_abs"], 6),
              [round(x, 4) for x in v["requested"]])


if __name__ == "__main__":
    main()
