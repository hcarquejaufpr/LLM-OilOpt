#!/usr/bin/env python3
"""MILP for the public P2 v3 problem. Requires scipy and numpy."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
PROBLEM = ROOT / 'shared' / 'problema_p2_public.json'


class LP:
    def __init__(self):
        self.names, self.lb, self.ub, self.integer, self.rows = [], [], [], [], []

    def var(self, name, lo=-math.inf, hi=math.inf, binary=False):
        j = len(self.names)
        self.names.append(name)
        self.lb.append(0.0 if binary else lo)
        self.ub.append(1.0 if binary else hi)
        self.integer.append(1 if binary else 0)
        return j

    def row(self, terms, lo=-math.inf, hi=math.inf):
        self.rows.append((terms, lo, hi))

    def eq(self, terms, rhs):
        self.row(terms, rhs, rhs)

    def solve(self, objective):
        rr, cc, vv, lower, upper = [], [], [], [], []
        for k, (terms, lo, hi) in enumerate(self.rows):
            for j, val in terms.items():
                if val:
                    rr.append(k); cc.append(j); vv.append(val)
            lower.append(lo); upper.append(hi)
        matrix = coo_matrix((vv, (rr, cc)), shape=(len(self.rows), len(self.names))).tocsr()
        c = np.zeros(len(self.names))
        for j, val in objective.items():
            c[j] = val
        return milp(c, integrality=np.array(self.integer),
                    bounds=Bounds(self.lb, self.ub),
                    constraints=LinearConstraint(matrix, lower, upper),
                    options={'mip_rel_gap': 1e-9})


def pwl(lp, name, argument, spec):
    """Adjacent segment formulation: one segment and convex node weights."""
    xs, vs = spec['nodes'], spec['values']
    n = len(xs)
    assert n == len(vs) and n >= 2 and all(xs[i] < xs[i+1] for i in range(n-1))
    value = lp.var(name, min(vs), max(vs))
    lam = [lp.var(f'{name}_lambda_{i}', 0, 1) for i in range(n)]
    seg = [lp.var(f'{name}_segment_{i}', binary=True) for i in range(n-1)]
    lp.eq({j: 1 for j in lam}, 1)
    lp.eq({j: 1 for j in seg}, 1)
    lp.eq({argument: 1, **{lam[i]: -xs[i] for i in range(n)}}, 0)
    lp.eq({value: 1, **{lam[i]: -vs[i] for i in range(n)}}, 0)
    for i, j in enumerate(lam):
        terms = {j: 1}
        if i: terms[seg[i-1]] = -1
        if i < n-1: terms[seg[i]] = -1
        lp.row(terms, hi=0)
    return value


def interp(x, spec):
    xs, vs = spec['nodes'], spec['values']
    if x < xs[0] - 1e-6 or x > xs[-1] + 1e-6:
        raise ValueError(f'PWL argument {x} outside [{xs[0]}, {xs[-1]}]')
    for i in range(len(xs)-1):
        if x <= xs[i+1] + 1e-10:
            t = (x-xs[i])/(xs[i+1]-xs[i])
            return vs[i] + t*(vs[i+1]-vs[i])
    return vs[-1]


def build(p, epsilon):
    lp = LP()
    wells = p['wells']
    domains = p['decision_variables']['domains']
    fixed = p['decision_variables']['fixed']
    cap = p['decision_variables']['capacity_requested']
    q = {w: lp.var('q_'+w, *(fixed[w], fixed[w]) if w in fixed else domains[w]) for w in wells}
    lp.row({q[w]: 1 for w in wells}, hi=cap)
    Q = lp.var('Q')
    lp.eq({Q: 20, **{q[w]: -1 for w in wells}}, -600)
    own = {}
    for w in wells:
        own[w] = lp.var('own_'+w)
        lp.eq({own[w]: 20, q[w]: -1}, -80)
    pf = p['pwl_functions']
    own2 = pwl(lp, 'own2_INJECT3', own['INJECT3'], pf['own2_INJECT3'])
    Q2 = pwl(lp, 'Q2', Q, pf['Q2'])
    d, r, real, clipping_M = {}, {}, {}, {}
    for w in wells:
        qlo, qhi = (fixed[w], fixed[w]) if w in fixed else domains[w]
        dlo, dhi = p['deficit_bounds'][w]
        if qlo <= 0 or dlo > dhi:
            raise ValueError('Invalid public q or deficit bounds')
        d[w] = lp.var('deficit_'+w)  # auxiliary bounds are deliberately not variable bounds
        r[w] = lp.var('r_'+w, 0, qhi)  # generic consequence of clipping; no realized_bounds used
        real[w] = lp.var('real_'+w)
        lp.eq({real[w]: 20, r[w]: -1}, -80)
        m = p['layer1_deficit_models'][w]
        terms = {d[w]: 1}
        for key, idx in [('own', own[w]), ('Q', Q), ('own2', own2), ('Q2', Q2)]:
            if key in m:
                terms[idx] = -m[key]
        lp.eq(terms, m.get('intercept', 0))
        # These are logical consequences of clipping, useful to tighten the relaxation.
        lp.row({r[w]: 1, q[w]: -1}, hi=0)
        neg, mid, high = [lp.var(f'{w}_{s}', binary=True) for s in ('negative', 'middle', 'high')]
        lp.eq({neg: 1, mid: 1, high: 1}, 1)
        M = {
            'd_le_0': max(0, dhi),
            'd_ge_0': max(0, -dlo),
            'd_le_q': max(0, dhi-qlo),
            'd_ge_q': max(0, qhi-dlo),
            'r_eq_q': qhi,
            'r_eq_q_minus_d_lower': max(0, qhi-dlo),
            'r_eq_q_minus_d_upper': max(0, dhi),
            'r_eq_0': qhi,
        }
        clipping_M[w] = M
        # neg: d <= 0, r = q
        lp.row({d[w]: 1, neg: M['d_le_0']}, hi=M['d_le_0'])
        lp.row({r[w]: 1, q[w]: -1, neg: -M['r_eq_q']}, lo=-M['r_eq_q'])
        # mid: 0 <= d <= q, r = q-d
        lp.row({d[w]: 1, mid: -M['d_ge_0']}, lo=-M['d_ge_0'])
        lp.row({d[w]: 1, q[w]: -1, mid: M['d_le_q']}, hi=M['d_le_q'])
        lp.row({r[w]: 1, q[w]: -1, d[w]: 1, mid: -M['r_eq_q_minus_d_lower']}, lo=-M['r_eq_q_minus_d_lower'])
        lp.row({r[w]: 1, q[w]: -1, d[w]: 1, mid: M['r_eq_q_minus_d_upper']}, hi=M['r_eq_q_minus_d_upper'])
        # high: d >= q, r = 0
        lp.row({d[w]: 1, q[w]: -1, high: -M['d_ge_q']}, lo=-M['d_ge_q'])
        lp.row({r[w]: 1, high: M['r_eq_0']}, hi=M['r_eq_0'])
    R = lp.var('R')
    lp.eq({R: 1, **{r[w]: -1 for w in wells}}, 0)
    Z = lp.var('Z')
    lp.eq({Z: 20, R: -1}, -600)
    Z2 = pwl(lp, 'Z2', Z, pf['Z2'])
    fc = p['layer2']['FOPT']['coefficients']
    FOPT = lp.var('FOPT')
    lp.eq({FOPT: 1, **{real[w]: -fc[f'r{i}'] for i,w in enumerate(wells,1)}, Z2: -fc['Z2']}, fc['intercept'])
    wc = p['layer2']['FWIT']['coefficients']
    FWIT = lp.var('FWIT')
    lp.eq({FWIT: 1, Z: -wc['Z']}, wc['intercept'])
    lp.row({FWIT: 1}, hi=epsilon)
    return lp, {'q': q, 'r': r, 'd': d, 'R': R, 'Q': Q, 'Z': Z, 'FOPT': FOPT, 'FWIT': FWIT,
                'own2': own2, 'Q2': Q2, 'Z2': Z2, 'M': clipping_M}


def audit(p, x, idx, epsilon):
    wells = p['wells']
    q = {w: float(x[idx['q'][w]]) for w in wells}
    r = {w: float(x[idx['r'][w]]) for w in wells}
    Q = (sum(q.values())-600)/20
    deficits = {}
    for w in wells:
        m = p['layer1_deficit_models'][w]
        ow = (q[w]-80)/20
        d = m.get('intercept',0) + m.get('own',0)*ow + m.get('Q',0)*Q
        if w == 'INJECT3':
            d += m['own2']*interp(ow,p['pwl_functions']['own2_INJECT3'])
            d += m['Q2']*interp(Q,p['pwl_functions']['Q2'])
        deficits[w] = d
        if abs(r[w] - (q[w]-min(max(d,0),q[w]))) > 2e-5:
            raise RuntimeError(f'Clipping audit failed at {w}')
        if abs(d-x[idx['d'][w]]) > 2e-5:
            raise RuntimeError(f'Deficit audit failed at {w}')
    R = sum(r.values()); Z = (R-600)/20
    z2 = interp(Z,p['pwl_functions']['Z2'])
    fc = p['layer2']['FOPT']['coefficients']
    fopt = fc['intercept'] + sum(fc[f'r{i}']*(r[w]-80)/20 for i,w in enumerate(wells,1)) + fc['Z2']*z2
    wc = p['layer2']['FWIT']['coefficients']
    fwit = wc['intercept']+wc['Z']*Z
    if abs(fopt-x[idx['FOPT']]) > 2e-4 or abs(fwit-x[idx['FWIT']]) > 2e-4:
        raise RuntimeError('Layer 2 audit failed')
    if fwit > epsilon + 1e-4 or sum(q.values()) > p['decision_variables']['capacity_requested']+1e-5:
        raise RuntimeError('Feasibility audit failed')
    return {'requested': q, 'realized': r, 'deficit': deficits, 'Q': Q, 'R': R, 'Z': Z,
            'Z2': z2, 'FOPT': fopt, 'FWIT': fwit}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--epsilon', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon):
        parser.error('--epsilon must be finite')
    p = json.loads(PROBLEM.read_text())
    lp, idx = build(p,args.epsilon)
    result = lp.solve({idx['FOPT']: -1})
    output = {'status': 'optimal' if result.status == 0 else 'solver_status_'+str(result.status),
              'epsilon_FWIT': args.epsilon, 'solver_gap': float(result.mip_gap) if result.x is not None else None,
              'solver_message': result.message}
    if result.x is not None:
        output.update(audit(p,result.x,idx,args.epsilon))
        output['clipping_big_m'] = idx['M']
        output['solver_objective'] = float(-result.fun)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    if result.status != 0:
        raise SystemExit(f'MILP did not solve optimally: {result.message}')


if __name__ == '__main__':
    main()
