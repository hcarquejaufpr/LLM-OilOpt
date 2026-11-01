#!/usr/bin/env python3
"""Auditable MILP for the public P2 v3 water-injection problem.

Requires scipy and numpy. All numerical model data come from the public JSON.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_array

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'shared' / 'problema_p2_public.json'


class Model:
    def __init__(self):
        self.variables = []
        self.rows = []
        self.objective = {}

    def var(self, name, lower=-math.inf, upper=math.inf, binary=False):
        i = len(self.variables)
        self.variables.append(dict(name=name, lower=lower, upper=upper, binary=binary))
        return i

    def row(self, name, coefficients, lower=-math.inf, upper=math.inf):
        c = {int(i): float(v) for i, v in coefficients.items() if v != 0}
        self.rows.append(dict(name=name, coefficients=c, lower=lower, upper=upper))

    def eq(self, name, coefficients, rhs=0):
        self.row(name, coefficients, rhs, rhs)

    def le(self, name, coefficients, rhs):
        self.row(name, coefficients, upper=rhs)

    def ge(self, name, coefficients, rhs):
        self.row(name, coefficients, lower=rhs)

    def pwl(self, name, x, y, nodes, values):
        assert len(nodes) == len(values) and len(nodes) >= 2
        assert all(a < b for a, b in zip(nodes, nodes[1:]))
        lam = [self.var(f'{name}.lambda[{j}]', 0, 1) for j in range(len(nodes))]
        seg = [self.var(f'{name}.segment[{j}]', 0, 1, True) for j in range(len(nodes)-1)]
        self.eq(f'{name}.lambda_sum', {v: 1 for v in lam}, 1)
        self.eq(f'{name}.segment_sum', {v: 1 for v in seg}, 1)
        self.eq(f'{name}.argument', {x: 1, **{v: -nodes[j] for j, v in enumerate(lam)}})
        self.eq(f'{name}.value', {y: 1, **{v: -values[j] for j, v in enumerate(lam)}})
        for j, v in enumerate(lam):
            c = {v: 1}
            if j > 0:
                c[seg[j-1]] = -1
            if j < len(seg):
                c[seg[j]] = -1
            self.le(f'{name}.adjacency[{j}]', c, 0)

    def gated_le(self, name, expression, M, z):
        assert M >= 0 and math.isfinite(M)
        c = dict(expression)
        c[z] = c.get(z, 0) + M
        self.le(name, c, M)

    def gated_eq(self, name, expression, M, z):
        self.gated_le(name+'.upper', expression, M, z)
        self.gated_le(name+'.lower', {i: -v for i, v in expression.items()}, M, z)

    def matrix(self):
        ii, jj, vv = [], [], []
        for i, row in enumerate(self.rows):
            for j, v in row['coefficients'].items():
                ii.append(i); jj.append(j); vv.append(v)
        A = coo_array((np.asarray(vv, dtype=float),
                       (np.asarray(ii, dtype=np.int32), np.asarray(jj, dtype=np.int32))),
                      shape=(len(self.rows), len(self.variables))).tocsr()
        return A

    def export(self, source_sha256, auxiliary_M):
        epsilon_row = next(r for r in self.rows if r['name'] == 'epsilon_FWIT')
        return dict(schema='p2-v3-sparse-milp-v1', source='shared/problema_p2_public.json',
                    source_sha256=source_sha256, sense='maximize',
                    parameter=dict(name='epsilon_FWIT', row='epsilon_FWIT',
                                   field='upper', instance_value=epsilon_row['upper']),
                    objective={str(i): v for i, v in self.objective.items()},
                    variables=[dict(v, lower=None if not math.isfinite(v['lower']) else v['lower'],
                                    upper=None if not math.isfinite(v['upper']) else v['upper'])
                               for v in self.variables],
                    constraints=[dict(name=r['name'], coefficients={str(i): v for i, v in r['coefficients'].items()},
                                      lower=None if not math.isfinite(r['lower']) else r['lower'],
                                      upper=None if not math.isfinite(r['upper']) else r['upper']) for r in self.rows],
                    clipping_big_m=auxiliary_M,
                    note='No direct bound is imposed on deficit or realized variables; auxiliary bounds only set Big-M.')


def build(data, epsilon):
    m = Model()
    wells = data['wells']
    domains = {w: tuple(data['decision_variables']['domains'].get(w,
                        [data['decision_variables']['fixed'][w]]*2 if w in data['decision_variables']['fixed'] else []))
               for w in wells}
    q = {w: m.var(f'q.{w}', *domains[w]) for w in wells}
    r = {w: m.var(f'r.{w}') for w in wells}
    d = {w: m.var(f'deficit.{w}') for w in wells}
    Q = m.var('Q'); Z = m.var('Z'); own3 = m.var('own.INJECT3')
    own2 = m.var('PWL.own2_INJECT3'); Q2 = m.var('PWL.Q2'); Z2 = m.var('PWL.Z2')
    fopt = m.var('FOPT'); fwit = m.var('FWIT')
    m.eq('feature.Q', {Q: 20, **{q[w]: -1 for w in wells}}, -600)
    m.eq('feature.Z', {Z: 20, **{r[w]: -1 for w in wells}}, -600)
    m.eq('feature.own.INJECT3', {own3: 20, q['INJECT3']: -1}, -80)
    m.le('capacity_requested', {q[w]: 1 for w in wells}, data['decision_variables']['capacity_requested'])
    for name, x, y in [('own2_INJECT3', own3, own2), ('Q2', Q, Q2), ('Z2', Z, Z2)]:
        p = data['pwl_functions'][name]
        m.pwl(name, x, y, p['nodes'], p['values'])

    big_m = {}
    for w in wells:
        c = data['layer1_deficit_models'][w]
        expr = {d[w]: 1}
        rhs = c.get('intercept', 0)
        for key, idx, scale in [('Q', Q, 1), ('Q2', Q2, 1), ('own2', own2, 1),
                                ('own', q[w], 1/20)]:
            if key in c:
                expr[idx] = expr.get(idx, 0) - c[key]*scale
                if key == 'own':
                    rhs -= 4*c[key]
        m.eq(f'deficit_model.{w}', expr, rhs)
        L, U = map(float, data['deficit_bounds'][w])
        qL, qU = domains[w]
        # Valid for every intended clipping solution: 0 <= r <= q <= qU.
        # No auxiliary deficit/realized bounds are added as model constraints.
        M = dict(d_le_zero=max(0, U), d_ge_zero=max(0, -L),
                 d_le_q=max(0, U-qL), d_ge_q=max(0, qU-L),
                 r_eq_q=qU, r_eq_q_minus_d=max(qU-L, qU-qL+U), r_eq_zero=qU)
        big_m[w] = M
        a = [m.var(f'clip.{w}.regime[{k}]', 0, 1, True) for k in range(3)]
        m.eq(f'clip.{w}.one_regime', {v: 1 for v in a}, 1)
        # 0: d<=0, r=q. 1: 0<=d<=q, r=q-d. 2: d>=q, r=0.
        m.gated_le(f'clip.{w}.zero.d', {d[w]: 1}, M['d_le_zero'], a[0])
        m.gated_eq(f'clip.{w}.zero.r', {r[w]: 1, q[w]: -1}, M['r_eq_q'], a[0])
        m.gated_le(f'clip.{w}.middle.d_lower', {d[w]: -1}, M['d_ge_zero'], a[1])
        m.gated_le(f'clip.{w}.middle.d_upper', {d[w]: 1, q[w]: -1}, M['d_le_q'], a[1])
        m.gated_eq(f'clip.{w}.middle.r', {r[w]: 1, q[w]: -1, d[w]: 1}, M['r_eq_q_minus_d'], a[1])
        m.gated_le(f'clip.{w}.upper.d', {q[w]: 1, d[w]: -1}, M['d_ge_q'], a[2])
        m.gated_eq(f'clip.{w}.upper.r', {r[w]: 1}, M['r_eq_zero'], a[2])

    fc = data['layer2']['FOPT']['coefficients']
    fexpr = {fopt: 1, Z2: -fc['Z2']}
    rhs = fc['intercept']
    for j, w in enumerate(wells, 1):
        beta = fc[f'r{j}']
        fexpr[r[w]] = -beta/20
        rhs -= 4*beta
    m.eq('FOPT_definition', fexpr, rhs)
    wc = data['layer2']['FWIT']['coefficients']
    m.eq('FWIT_definition', {fwit: 1, Z: -wc['Z']}, wc['intercept'])
    m.le('epsilon_FWIT', {fwit: 1}, epsilon)
    m.objective[fopt] = 1
    return m, dict(q=q, r=r, d=d, Q=Q, Z=Z, FOPT=fopt, FWIT=fwit), big_m


def interp(nodes, values, x):
    if x < nodes[0]-1e-7 or x > nodes[-1]+1e-7:
        raise ValueError(f'PWL argument {x} outside [{nodes[0]}, {nodes[-1]}]')
    j = min(max(int(np.searchsorted(nodes, x, side='right')-1), 0), len(nodes)-2)
    t = (x-nodes[j])/(nodes[j+1]-nodes[j])
    return values[j]*(1-t)+values[j+1]*t


def independent_check(data, result):
    wells = data['wells']; q = result['requested']; r = result['realized']
    Q = (sum(q.values())-600)/20; Z = (sum(r.values())-600)/20
    pwls = data['pwl_functions']
    own2 = interp(**{k: pwls['own2_INJECT3'][k] for k in ('nodes', 'values')}, x=(q['INJECT3']-80)/20)
    Q2 = interp(**{k: pwls['Q2'][k] for k in ('nodes', 'values')}, x=Q)
    Z2 = interp(**{k: pwls['Z2'][k] for k in ('nodes', 'values')}, x=Z)
    dr = {}
    for w in wells:
        c = data['layer1_deficit_models'][w]
        own = (q[w]-80)/20
        deficit = c.get('intercept', 0)+c.get('own', 0)*own+c.get('Q', 0)*Q+c.get('own2', 0)*own2+c.get('Q2', 0)*Q2
        dr[w] = abs(r[w]-(q[w]-min(max(deficit, 0), q[w])))
    fc = data['layer2']['FOPT']['coefficients']
    f = fc['intercept']+sum(fc[f'r{i}']*(r[w]-80)/20 for i,w in enumerate(wells,1))+fc['Z2']*Z2
    wc = data['layer2']['FWIT']['coefficients']; fw = wc['intercept']+wc['Z']*Z
    return dict(max_clipping_residual=max(dr.values()), clipping_residuals=dr,
                FOPT_residual=abs(result['FOPT']-f), FWIT_residual=abs(result['FWIT']-fw),
                capacity_slack=data['decision_variables']['capacity_requested']-sum(q.values()),
                epsilon_slack=result['epsilon_FWIT']-result['FWIT'])


def solve(data, epsilon, export_model=None):
    m, ids, big_m = build(data, epsilon)
    if export_model is not None:
        sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        Path(export_model).write_text(json.dumps(m.export(sha, big_m), indent=2, allow_nan=False)+'\n')
    n = len(m.variables)
    obj = np.zeros(n)
    for i,v in m.objective.items(): obj[i] = -v
    bounds = Bounds([v['lower'] for v in m.variables], [v['upper'] for v in m.variables])
    integrality = np.array([int(v['binary']) for v in m.variables], dtype=np.int32)
    A = m.matrix()
    con = LinearConstraint(A, [x['lower'] for x in m.rows], [x['upper'] for x in m.rows])
    opt = milp(obj, integrality=integrality, bounds=bounds, constraints=con,
               options={'mip_rel_gap': 1e-8})
    status = {0:'optimal', 1:'limit_reached', 2:'infeasible', 3:'unbounded', 4:'solver_error'}.get(opt.status, 'unknown')
    result = dict(status=status, epsilon_FWIT=epsilon, requested=None, realized=None,
                  R=None, FOPT=None, FWIT=None, solver_gap=None,
                  solver_message=str(opt.message))
    if opt.x is not None:
        x = opt.x
        result.update(requested={w: float(x[i]) for w,i in ids['q'].items()},
                      realized={w: float(x[i]) for w,i in ids['r'].items()},
                      R=float(sum(x[i] for i in ids['r'].values())),
                      FOPT=float(x[ids['FOPT']]), FWIT=float(x[ids['FWIT']]),
                      solver_gap=float(opt.mip_gap) if opt.mip_gap is not None else None)
        result['independent_check'] = independent_check(data, result)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--epsilon', type=float, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--export-model', type=Path, help='Optional sparse model JSON')
    args = ap.parse_args()
    if not math.isfinite(args.epsilon): ap.error('epsilon must be finite')
    data = json.loads(SOURCE.read_text())
    result = solve(data, args.epsilon, args.export_model)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if result['status'] != 'optimal': raise SystemExit(f"Solver status: {result['status']}: {result['solver_message']}")


if __name__ == '__main__':
    main()
