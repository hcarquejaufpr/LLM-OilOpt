#!/usr/bin/env python3
"""MILP/PWL model of the public P2 v3 water-injection problem."""
import argparse
import hashlib
import json
import math
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'shared' / 'problema_p2_public.json'


class Model:
    def __init__(self):
        self.names, self.lo, self.hi, self.integer, self.obj = [], [], [], [], []
        self.rows = []

    def var(self, name, lo=-math.inf, hi=math.inf, integer=False, objective=0.0):
        j = len(self.names)
        self.names.append(name)
        self.lo.append(lo)
        self.hi.append(hi)
        self.integer.append(int(integer))
        self.obj.append(objective)
        return j

    def row(self, name, terms, lo=-math.inf, hi=math.inf):
        terms = {j: v for j, v in terms.items() if v != 0}
        self.rows.append((name, terms, lo, hi))

    def eq(self, name, terms, rhs):
        self.row(name, terms, rhs, rhs)

    def le(self, name, terms, rhs):
        self.row(name, terms, -math.inf, rhs)

    def pwl(self, name, x, nodes, values):
        assert len(nodes) == len(values) and len(nodes) >= 2
        assert all(nodes[k] < nodes[k + 1] for k in range(len(nodes) - 1))
        y = self.var(name)
        lambdas = [self.var(f'{name}.lambda[{k}]', 0, 1) for k in range(len(nodes))]
        segs = [self.var(f'{name}.segment[{k}]', 0, 1, True) for k in range(len(nodes) - 1)]
        self.eq(f'{name}.lambda_sum', {j: 1 for j in lambdas}, 1)
        self.eq(f'{name}.segment_sum', {j: 1 for j in segs}, 1)
        self.eq(f'{name}.argument', {x: 1, **{j: -v for j, v in zip(lambdas, nodes)}}, 0)
        self.eq(f'{name}.value', {y: 1, **{j: -v for j, v in zip(lambdas, values)}}, 0)
        for k, lam in enumerate(lambdas):
            terms = {lam: 1}
            if k > 0:
                terms[segs[k - 1]] = -1
            if k < len(segs):
                terms[segs[k]] = -1
            self.le(f'{name}.adjacency[{k}]', terms, 0)
        return y

    def export(self, source_hash, data):
        return {
            'schema': 'p2-v3-milp-matrix-1',
            'source_sha256': source_hash,
            'sense': 'maximize',
            'objective_note': 'FOPT affine expression, including real_i=(r_i-80)/20; constant in objective_constant',
            'objective_constant': data['layer2']['FOPT']['coefficients']['intercept'] - 4 * sum(data['layer2']['FOPT']['coefficients'][f'r{i}'] for i in range(1, 9)),
            'variables': [{'name': n, 'lower': finite(l), 'upper': finite(u), 'integer': bool(z), 'objective': c}
                          for n, l, u, z, c in zip(self.names, self.lo, self.hi, self.integer, self.obj)],
            'constraints': [{'name': n, 'coefficients': {self.names[j]: v for j, v in t.items()},
                             'lower': finite(l), 'upper': finite(u)} for n, t, l, u in self.rows],
            'epsilon_row': 'FWIT_limit',
            'note': 'Null bounds mean infinity. The epsilon row is set at solve time. Auxiliary deficit_bounds and realized_bounds are not variable bounds or explicit feasibility rows.'
        }

    def solve(self):
        ri, ci, vv = [], [], []
        for k, (_, terms, _, _) in enumerate(self.rows):
            for j, value in terms.items():
                ri.append(k); ci.append(j); vv.append(value)
        A = coo_matrix((vv, (ri, ci)), shape=(len(self.rows), len(self.names))).tocsr()
        constraints = LinearConstraint(A, [r[2] for r in self.rows], [r[3] for r in self.rows])
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Unrecognized options detected:.*')
            return milp(c=-np.asarray(self.obj), integrality=np.asarray(self.integer),
                        bounds=Bounds(self.lo, self.hi), constraints=constraints,
                        options={'mip_rel_gap': 1e-9, 'mip_feasibility_tolerance': 1e-9,
                                 'primal_feasibility_tolerance': 1e-9})


def finite(x):
    return x if math.isfinite(x) else None


def build(data, epsilon):
    m = Model()
    wells = data['wells']
    domains = data['decision_variables']['domains']
    fixed = data['decision_variables']['fixed']
    capacity = data['decision_variables']['capacity_requested']
    q, r, d = {}, {}, {}
    qlo, qhi = {}, {}
    beta = data['layer2']['FOPT']['coefficients']
    for i, w in enumerate(wells, 1):
        a, b = ([fixed[w]] * 2) if w in fixed else domains[w]
        qlo[w], qhi[w] = a, b
        q[w] = m.var(f'q.{w}', a, b)
        # 0 <= r <= max(q) follows from the clipping identity, and is not from realized_bounds.
        r[w] = m.var(f'r.{w}', 0, b, objective=beta[f'r{i}'] / 20)
        d[w] = m.var(f'd.{w}')
    m.le('requested_capacity', {q[w]: 1 for w in wells}, capacity)
    Q = m.var('Q')
    m.eq('Q_definition', {Q: 20, **{q[w]: -1 for w in wells}}, -600)
    Z = m.var('Z')
    m.eq('Z_definition', {Z: 20, **{r[w]: -1 for w in wells}}, -600)
    pwl = data['pwl_functions']
    Q2 = m.pwl('PWL.Q2', Q, pwl['Q2']['nodes'], pwl['Q2']['values'])
    own3 = m.var('own.INJECT3')
    m.eq('own.INJECT3_definition', {own3: 20, q['INJECT3']: -1}, -80)
    own3sq = m.pwl('PWL.own2_INJECT3', own3, pwl['own2_INJECT3']['nodes'], pwl['own2_INJECT3']['values'])
    Z2 = m.pwl('PWL.Z2', Z, pwl['Z2']['nodes'], pwl['Z2']['values'])
    m.obj[Z2] = beta['Z2']
    for w in wells:
        co = data['layer1_deficit_models'][w]
        terms = {d[w]: 1}
        rhs = co.get('intercept', 0)
        for key, j in [('own', q[w]), ('Q', Q), ('own2', own3sq), ('Q2', Q2)]:
            if key not in co:
                continue
            if key == 'own':
                terms[j] = terms.get(j, 0) - co[key] / 20
                rhs -= 4 * co[key]
            else:
                terms[j] = terms.get(j, 0) - co[key]
        m.eq(f'deficit.{w}', terms, rhs)
        if not co:
            # A zero deficit makes clipping identically r=q.
            m.eq(f'clip.{w}.zero_deficit', {r[w]: 1, q[w]: -1}, 0)
            continue
        L, U = data['deficit_bounds'][w]
        a, b = qlo[w], qhi[w]
        # Each M bounds the inactive implication, based on semantic intervals only.
        pad = 1e-7
        M_dlow = max(0, U) + pad
        M_dmidlow = max(0, -L) + pad
        M_dmidup = max(0, U - a) + pad
        M_dup = max(0, b - L) + pad
        M_rq = min(b, max(0, U)) + pad
        M_rqd = max(0, -L, U - a) + pad
        M_r0 = b + pad
        low = m.var(f'clip.{w}.low', 0, 1, True)
        mid = m.var(f'clip.{w}.middle', 0, 1, True)
        up = m.var(f'clip.{w}.upper', 0, 1, True)
        m.eq(f'clip.{w}.choose', {low: 1, mid: 1, up: 1}, 1)
        m.le(f'clip.{w}.low_sign', {d[w]: 1, low: M_dlow}, M_dlow)
        m.le(f'clip.{w}.mid_lower', {d[w]: -1, mid: M_dmidlow}, M_dmidlow)
        m.le(f'clip.{w}.mid_upper', {d[w]: 1, q[w]: -1, mid: M_dmidup}, M_dmidup)
        m.le(f'clip.{w}.upper_sign', {q[w]: 1, d[w]: -1, up: M_dup}, M_dup)
        for sign in (1, -1):
            m.le(f'clip.{w}.low_value.{sign}', {r[w]: sign, q[w]: -sign, low: M_rq}, M_rq)
            m.le(f'clip.{w}.mid_value.{sign}', {r[w]: sign, q[w]: -sign, d[w]: sign, mid: M_rqd}, M_rqd)
            m.le(f'clip.{w}.upper_value.{sign}', {r[w]: sign, up: M_r0}, M_r0)
    fwit = data['layer2']['FWIT']['coefficients']
    m.le('FWIT_limit', {Z: fwit['Z']}, epsilon - fwit['intercept'])
    return m, q, r, d, Q, Z, Q2, own3sq, Z2


def pwl_eval(x, nodes, values):
    if x < nodes[0] - 1e-6 or x > nodes[-1] + 1e-6:
        raise ValueError('PWL argument outside public node domain')
    for k in range(len(nodes) - 1):
        if x <= nodes[k + 1] + 1e-10:
            t = (x - nodes[k]) / (nodes[k + 1] - nodes[k])
            return values[k] + t * (values[k + 1] - values[k])
    return values[-1]


def audit(data, requested, realized, deficit, FOPT, FWIT):
    wells = data['wells']
    Q = (sum(requested.values()) - 600) / 20
    Z = (sum(realized.values()) - 600) / 20
    pwl = data['pwl_functions']
    q2 = pwl_eval(Q, pwl['Q2']['nodes'], pwl['Q2']['values'])
    own3 = (requested['INJECT3'] - 80) / 20
    own2 = pwl_eval(own3, pwl['own2_INJECT3']['nodes'], pwl['own2_INJECT3']['values'])
    z2 = pwl_eval(Z, pwl['Z2']['nodes'], pwl['Z2']['values'])
    dres, rres = {}, {}
    aux_d_violations, aux_r_violations = {}, {}
    for w in wells:
        co = data['layer1_deficit_models'][w]
        computed = (co.get('intercept', 0) + co.get('own', 0) * (requested[w] - 80) / 20
                    + co.get('Q', 0) * Q + co.get('own2', 0) * own2 + co.get('Q2', 0) * q2)
        dres[w] = deficit[w] - computed
        rres[w] = realized[w] - (requested[w] - min(max(computed, 0), requested[w]))
        dl, du = data['deficit_bounds'][w]
        rl, ru = data['realized_bounds'][w]
        aux_d_violations[w] = max(0, dl - deficit[w], deficit[w] - du)
        aux_r_violations[w] = max(0, rl - realized[w], realized[w] - ru)
    co = data['layer2']['FOPT']['coefficients']
    computed_fopt = co['intercept'] + sum(co[f'r{i}'] * (realized[w] - 80) / 20 for i, w in enumerate(wells, 1)) + co['Z2'] * z2
    fwc = data['layer2']['FWIT']['coefficients']
    computed_fwit = fwc['intercept'] + fwc['Z'] * Z
    return {'max_deficit_residual': max(abs(v) for v in dres.values()),
            'max_clipping_residual': max(abs(v) for v in rres.values()),
            'FOPT_residual': FOPT - computed_fopt, 'FWIT_residual': FWIT - computed_fwit,
            'max_auxiliary_deficit_bound_violation': max(aux_d_violations.values()),
            'max_auxiliary_realized_bound_violation': max(aux_r_violations.values()),
            'requested_capacity_slack': data['decision_variables']['capacity_requested'] - sum(requested.values())}


def run(epsilon, export_model=None):
    data = json.loads(SOURCE.read_text())
    m, q, r, d, Q, Z, Q2, own3sq, Z2 = build(data, epsilon)
    if export_model:
        export_model.write_text(json.dumps(m.export(hashlib.sha256(SOURCE.read_bytes()).hexdigest(), data), indent=2) + '\n')
    result = m.solve()
    status = {0: 'optimal', 1: 'limit_reached', 2: 'infeasible', 3: 'unbounded', 4: 'solver_error'}.get(result.status, 'solver_error')
    output = {'status': status, 'epsilon_FWIT': epsilon, 'requested': None, 'realized': None,
              'R': None, 'FOPT': None, 'FWIT': None,
              'solver_gap': finite(float(result.mip_gap)) if getattr(result, 'mip_gap', None) is not None else None,
              'solver_message': result.message}
    if result.x is None:
        return output
    x = result.x
    requested = {w: float(x[j]) for w, j in q.items()}
    realized = {w: float(x[j]) for w, j in r.items()}
    deficits = {w: float(x[j]) for w, j in d.items()}
    co = data['layer2']['FOPT']['coefficients']
    fopt = co['intercept'] + sum(co[f'r{i}'] * (realized[w] - 80) / 20 for i, w in enumerate(data['wells'], 1)) + co['Z2'] * float(x[Z2])
    fwc = data['layer2']['FWIT']['coefficients']
    fwit = fwc['intercept'] + fwc['Z'] * float(x[Z])
    output.update({'requested': requested, 'realized': realized, 'R': sum(realized.values()),
                   'FOPT': fopt, 'FWIT': fwit, 'deficit': deficits,
                   'Q': float(x[Q]), 'Z': float(x[Z]),
                   'PWL_Q2': float(x[Q2]), 'PWL_own2_INJECT3': float(x[own3sq]), 'PWL_Z2': float(x[Z2]),
                   'audit': audit(data, requested, realized, deficits, fopt, fwit)})
    return output


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--epsilon', type=float, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if not math.isfinite(args.epsilon):
        p.error('--epsilon must be finite')
    output = run(args.epsilon)
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    if output['status'] not in ('optimal', 'infeasible'):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
