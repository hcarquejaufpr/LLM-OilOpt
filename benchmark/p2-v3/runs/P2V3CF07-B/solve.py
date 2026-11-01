#!/usr/bin/env python3
"""Public P2 v3 MILP; scipy.optimize.milp (HiGHS)."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
PROBLEM = ROOT / 'shared' / 'problema_p2_public.json'


class MILP:
    def __init__(self):
        self.names, self.lower, self.upper, self.integer = [], [], [], []
        self.rows, self.lhs, self.rhs = [], [], []

    def var(self, name, lo=-np.inf, hi=np.inf, binary=False):
        k = len(self.names)
        self.names.append(name)
        self.lower.append(0 if binary else lo)
        self.upper.append(1 if binary else hi)
        self.integer.append(int(binary))
        return k

    def row(self, terms, lo=-np.inf, hi=np.inf):
        self.rows.append({k: v for k, v in terms.items() if v})
        self.lhs.append(lo)
        self.rhs.append(hi)

    def eq(self, terms, rhs=0):
        self.row(terms, rhs, rhs)

    def solve(self, objective):
        ii, jj, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, v in row.items():
                ii.append(i); jj.append(j); vv.append(v)
        a = coo_matrix((vv, (ii, jj)), shape=(len(self.rows), len(self.names))).tocsr()
        c = np.zeros(len(self.names)); c[objective] = -1
        result = milp(c, integrality=np.array(self.integer),
                      bounds=Bounds(self.lower, self.upper),
                      constraints=LinearConstraint(a, self.lhs, self.rhs),
                      options={'mip_rel_gap': 1e-9, 'time_limit': 110})
        return result


def pwl(m, label, argument, spec):
    nodes, values = spec['nodes'], spec['values']
    n = len(nodes)
    y = m.var(label)
    lam = [m.var(f'{label}.lambda.{j}', 0, 1) for j in range(n)]
    seg = [m.var(f'{label}.segment.{j}', binary=True) for j in range(n - 1)]
    m.eq({k: 1 for k in lam}, 1)
    m.eq({k: 1 for k in seg}, 1)
    m.eq({argument: 1, **{k: -float(x) for k, x in zip(lam, nodes)}})
    m.eq({y: 1, **{k: -float(v) for k, v in zip(lam, values)}})
    for j, k in enumerate(lam):
        t = {k: 1}
        if j > 0: t[seg[j - 1]] = -1
        if j < n - 1: t[seg[j]] = -1
        m.row(t, hi=0)
    return y


def clipping(m, label, q, d, r, q_bounds, d_bounds):
    """Three disjoint regimes, with ties allowed at transition points.

    M values are valid over the *union* of the other two regimes. The
    auxiliary deficit bounds and q domain only derive M; neither is added
    as a constraint on d or r.
    """
    qlo, qhi = q_bounds
    dlo, dhi = d_bounds
    a, b, c = [m.var(f'{label}.regime.{x}', binary=True)
               for x in ('no_deficit', 'partial', 'full')]
    m.eq({a: 1, b: 1, c: 1}, 1)
    # a=1 => d<=0, r=q. For b/c inactive: d<=dhi, |r-q|<=qhi.
    m.row({d: 1, a: max(0, dhi)}, hi=max(0, dhi))
    m.row({r: 1, q: -1, a: qhi}, hi=qhi)
    m.row({q: 1, r: -1, a: qhi}, hi=qhi)
    # b=1 => 0<=d<=q, r=q-d.
    neg_dlo = max(0, -dlo)
    over_q = max(0, dhi - qlo)
    m.row({d: 1, b: over_q, q: -1}, hi=over_q)
    m.row({d: 1, b: -neg_dlo}, lo=-neg_dlo)
    mb = max(neg_dlo, over_q)
    m.row({r: 1, q: -1, d: 1, b: mb}, hi=mb)
    m.row({q: 1, r: -1, d: -1, b: mb}, hi=mb)
    # c=1 => d>=q, r=0. In a/b, q-d <= qhi-dlo, r<=qhi.
    under_q = max(0, qhi - dlo)
    m.row({d: 1, q: -1, c: -under_q}, lo=-under_q)
    m.row({r: 1, c: qhi}, hi=qhi)
    m.row({r: -1, c: qhi}, hi=qhi)
    return {'d_upper': max(0, dhi), 'd_lower': neg_dlo,
            'd_minus_q_upper': over_q, 'd_minus_q_lower': under_q,
            'r_minus_q': qhi, 'r_minus_q_plus_d': mb, 'r_zero': qhi}


def build(data, epsilon):
    m = MILP()
    wells = data['wells']
    dv = data['decision_variables']
    qb = {w: dv['domains'][w] if w in dv['domains'] else
          [dv['fixed'][w], dv['fixed'][w]] for w in wells}
    q = {w: m.var(f'q.{w}', *qb[w]) for w in wells}
    r = {w: m.var(f'r.{w}') for w in wells}
    own = {w: m.var(f'own.{w}') for w in wells}
    real = {w: m.var(f'real.{w}') for w in wells}
    d = {w: m.var(f'deficit.{w}') for w in wells}
    Q, R, Z = [m.var(x) for x in ('Q', 'R', 'Z')]
    m.row({k: 1 for k in q.values()}, hi=dv['capacity_requested'])
    m.eq({Q: 20, **{k: -1 for k in q.values()}}, -600)
    for w in wells:
        m.eq({own[w]: 20, q[w]: -1}, -80)
        m.eq({real[w]: 20, r[w]: -1}, -80)
    P = data['pwl_functions']
    own2 = pwl(m, 'own2_INJECT3', own['INJECT3'], P['own2_INJECT3'])
    Q2 = pwl(m, 'Q2', Q, P['Q2'])
    m.eq({R: 1, **{k: -1 for k in r.values()}})
    m.eq({Z: 20, R: -1}, -600)
    Z2 = pwl(m, 'Z2', Z, P['Z2'])
    constants = {}
    for w in wells:
        spec = data['layer1_deficit_models'][w]
        t = {d[w]: 1}
        if 'own' in spec: t[own[w]] = -spec['own']
        if 'Q' in spec: t[Q] = -spec['Q']
        if 'own2' in spec: t[own2] = -spec['own2']
        if 'Q2' in spec: t[Q2] = -spec['Q2']
        m.eq(t, spec.get('intercept', 0))
        constants[w] = clipping(m, w, q[w], d[w], r[w], qb[w],
                                data['deficit_bounds'][w])
    fopt, fwit = m.var('FOPT'), m.var('FWIT')
    fc = data['layer2']['FOPT']['coefficients']
    terms = {fopt: 1, Z2: -fc['Z2']}
    for i, w in enumerate(wells, 1):
        terms[real[w]] = -fc[f'r{i}']
    m.eq(terms, fc['intercept'])
    wc = data['layer2']['FWIT']['coefficients']
    m.eq({fwit: 1, Z: -wc['Z']}, wc['intercept'])
    m.row({fwit: 1}, hi=epsilon)
    return m, {'q': q, 'r': r, 'd': d, 'Q': Q, 'R': R, 'Z': Z,
               'own': own, 'real': real, 'own2': own2, 'Q2': Q2,
               'Z2': Z2, 'FOPT': fopt, 'FWIT': fwit}, constants


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--epsilon', required=True, type=float)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.epsilon): parser.error('epsilon must be finite')
    data = json.loads(PROBLEM.read_text())
    m, v, constants = build(data, args.epsilon)
    result = m.solve(v['FOPT'])
    if result.x is None or result.status != 0:
        raise RuntimeError(f'HiGHS failed: status={result.status}; {result.message}')
    x = result.x
    wells = data['wells']
    out = {
        'status': 'optimal', 'epsilon_FWIT': args.epsilon,
        'requested': {w: float(x[v['q'][w]]) for w in wells},
        'realized': {w: float(x[v['r'][w]]) for w in wells},
        'deficits': {w: float(x[v['d'][w]]) for w in wells},
        'Q': float(x[v['Q']]), 'R': float(x[v['R']]), 'Z': float(x[v['Z']]),
        'own': {w: float(x[v['own'][w]]) for w in wells},
        'real': {w: float(x[v['real'][w]]) for w in wells},
        'own2_INJECT3': float(x[v['own2']]), 'Q2': float(x[v['Q2']]),
        'Z2': float(x[v['Z2']]),
        'FOPT': float(x[v['FOPT']]), 'FWIT': float(x[v['FWIT']]),
        'solver_gap': float(result.mip_gap),
        'solver_message': result.message,
        'regimes': {w: max(('no_deficit', 'partial', 'full'),
                           key=lambda name: x[m.names.index(f'{w}.regime.{name}')])
                    for w in wells},
        'model_size': {'variables': len(m.names), 'constraints': len(m.rows)},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
