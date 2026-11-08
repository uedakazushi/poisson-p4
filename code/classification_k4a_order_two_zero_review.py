#!/usr/bin/env python3
"""Independent Fraction reconstruction of the full K4a order-two-zero operator.

No producer, arithmetic helper or CAS is imported. The existing same-orbit
curve value e=1=t/174 is read back from its full 150-coordinate graph.
"""
import ast
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[1]
JOB = 'research/classification/jobs/k4a_order_two_zero_separate_review.json'
OUTPUT = 'reports/CLASSIFICATION_K4A_ORDER_TWO_ZERO_SEPARATE_REVIEW_20261108.json'


def add(*terms):
    out = {}
    for term in terms:
        for exp, c in term.items():
            out[exp] = out.get(exp, F(0)) + c
            if not out[exp]:
                del out[exp]
    return out


def mul(a, b):
    out = {}
    for e, c in a.items():
        for g, d in b.items():
            h = tuple(x+y for x, y in zip(e, g))
            out[h] = out.get(h, F(0)) + c*d
            if not out[h]:
                del out[h]
    return out


def scale(p, c):
    return {e: c*v for e, v in p.items() if c*v}


def diff(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            h = list(e)
            h[i] -= 1
            out[tuple(h)] = c*e[i]
    return out


def parse(text, names):
    n = len(names)
    def walk(a):
        if isinstance(a, ast.Constant) and isinstance(a.value, int):
            return {(0,)*n: F(a.value)} if a.value else {}
        if isinstance(a, ast.Name):
            e = [0]*n
            e[names.index(a.id)] = 1
            return {tuple(e): F(1)}
        if isinstance(a, ast.UnaryOp):
            p = walk(a.operand)
            return scale(p, -1) if isinstance(a.op, ast.USub) else p
        if isinstance(a, ast.BinOp):
            p = walk(a.left)
            if isinstance(a.op, ast.Pow):
                k = int(ast.literal_eval(a.right))
                if k < 0:
                    assert len(p) == 1
                    e, c = next(iter(p.items()))
                    return {tuple(k*x for x in e): c**k}
                out = {(0,)*n: F(1)}
                for _ in range(k):
                    out = mul(out, p)
                return out
            q = walk(a.right)
            if isinstance(a.op, ast.Add):
                return add(p, q)
            if isinstance(a.op, ast.Sub):
                return add(p, scale(q, -1))
            if isinstance(a.op, ast.Mult):
                return mul(p, q)
            if isinstance(a.op, ast.Div):
                assert len(q) == 1
                e, c = next(iter(q.items()))
                return mul(p, {tuple(-x for x in e): 1/c})
        raise ValueError(ast.dump(a))
    return walk(ast.parse(str(text), mode='eval').body)


def det(a):
    b = [list(row) for row in a]
    value = F(1)
    for i in range(len(b)):
        j = next(j for j in range(i, len(b)) if b[j][i])
        if i != j:
            b[i], b[j] = b[j], b[i]
            value = -value
        pivot = b[i][i]
        value *= pivot
        for k in range(i+1, len(b)):
            if b[k][i]:
                ratio = b[k][i]/pivot
                for ell in range(i+1, len(b)):
                    b[k][ell] -= ratio*b[i][ell]
                b[k][i] = F(0)
    return value


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (20, 21))
    resource.setrlimit(resource.RLIMIT_AS, (512*1024**2,)*2)
    begun = time.monotonic()
    job = json.loads((ROOT/JOB).read_bytes())
    for r in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest() == r['sha256'], r['path']
    def load(path):
        return json.loads((ROOT/path).read_bytes())
    graph = load('reports/CLASSIFICATION_K4A_COMPLETE_TORUS_SLICE_20261108.json')
    original = load('reports/CLASSIFICATION_K4A_CANONICAL_FULL_AMBIENT_20261108.json')
    producer = load('reports/CLASSIFICATION_K4A_ORDER_TWO_ZERO_LOCUS_20261108.json')
    pairs = list(itertools.combinations(range(5), 2))
    coords = list(itertools.combinations_with_replacement(range(5), 2))
    exps = []
    for i, j in coords:
        e = [0]*5
        e[i] += 1
        e[j] += 1
        exps.append(tuple(e))
    coeff_order = original['full_original_coefficient_order']
    assert len(coeff_order) == 150
    assert {(tuple(r['pair']), tuple(r['exponent'])) for r in coeff_order} == {(p,e) for p in pairs for e in exps}
    p = {ij: {} for ij in pairs}
    for key, value in graph['graph_coefficients_sparse'].items():
        i, degree = map(int, key.split(','))
        r = coeff_order[i]
        ij, e = tuple(r['pair']), tuple(r['exponent'])
        p[ij] = add(p[ij], {e: F(value)*174**degree})
    assert all(c.denominator == 1 for q in p.values() for c in q.values())
    for ij in pairs:
        assert p[ij] == parse(producer['generic_actual_brackets'].get(','.join(map(str,ij)), '0'), ['c','a','b','u','v'])
    x = [{tuple(int(k == i) for k in range(5)): F(1)} for i in range(5)]
    rows, polys = [], []
    for i, j, k in itertools.combinations(range(5), 3):
        cubic = add(mul(x[i], p[j,k]), scale(mul(x[j], p[i,k]), -1), mul(x[k], p[i,j]))
        for ell in range(5):
            rows.append({'triple':[i,j,k], 'derivative':ell})
            polys.append(diff(cubic, ell))
    assert rows == producer['operator_equation_order']
    assert [list(q) for q in coords] == producer['operator_monomial_order']
    a = [[q.get(e, F(0)) for e in exps] for q in polys]
    sparse = {f'{i},{j}':str(c) for i, row in enumerate(a) for j,c in enumerate(row) if c}
    assert sparse == producer['full_operator_sparse']
    assert all(parse(text, ['c','a','b','u','v']) == q for text,q in zip(producer['all50_quadratic_order2_zero_equations'],polys))
    ir, ic = producer['independent_rows'], producer['independent_columns']
    assert len(ir) == len(set(ir)) == len(ic) == 15 and ic == list(range(15))
    d = det([[a[i][j] for j in ic] for i in ir])
    assert d and d == F(producer['nonzero_minor_determinant'])
    completed = []
    names = ['x0','x1','x2','x3','x4','p0','t0','t1','t2','t3']
    for name in ['T36','T38']:
        family = load(f'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/{name}.json')
        brackets = {tuple(map(int,k.split(','))):parse(v,names) for k,v in family['homogeneous_bivector'].items()}
        for q in brackets.values():
            assert all(e[0] == 0 and sum(e[:5]) == 2 for e in q)
        xx = [{tuple(int(k == i) for k in range(10)): F(1)} for i in range(5)]
        zerojet_count = 0
        for i,j,k in itertools.combinations(range(5),3):
            cubic = add(mul(xx[i],brackets.get((j,k),{})),scale(mul(xx[j],brackets.get((i,k),{})),-1),mul(xx[k],brackets.get((i,j),{})))
            for ell in range(5):
                derivative = diff(cubic,ell)
                # At q=[1:0:0:0:0], coefficients in all p0,t0..t3 remain symbolic.
                assert not {e:c for e,c in derivative.items() if all(e[j] == 0 for j in range(1,5))}
                zerojet_count += 1
        completed.append({'component':name,'all_full_source_coefficients_independent_x0':True,'all_base_quadratic':True,'all50_first_derivative_jets_identically_zero_at_vertex':zerojet_count == 50,'all_period_and_smoothing_parameters_retained':True})
    source = Path(__file__).read_bytes()
    result = {'task_id':'E_reduced_nonnormal','all_independent_checks_passed':True,
        'method':'Standard-library Fraction and independent sparse Laurent-polynomial parser; no producer/helper/CAS imports.',
        'same_existing_dense_orbit_parameter':'e=1, original t=174; accepted complete-slice one-orbit theorem is reused, no new anchor/rank/profile.',
        'full_actual_graph_evaluation_matches_all10_brackets':True,'full50_equations_and_50x15_operator_reconstructed':True,
        'full_operator_shape':[50,15],'full_operator_sparse':sparse,'independent_rows':ir,'independent_columns':ic,
        'nonzero_minor_determinant':str(d),'rank_QQ':15,'full_kernel_dimension':0,
        'order_two_zero_projective_scheme_empty_over_C':True,'completed_source_vertex_jets':completed,
        'pure_mathematical_bridge':'Intrinsic first-jet/order-two vanishing of E wedge P is representative/PGL invariant; proper P4 flag makes existence closed. Completed T36/T38 whole sources satisfy it. Same dense24 orbit violates it; accepted canonical smooth unique24 point excludes canonical K4a source from those components.',
        'T36_T38_whole_parents_excluded':True,'canonical_K4a_source_outside_T36_T38_by_unique_component':True,
        'T37_membership_or_new_component_proved':False,'new_component_count':0,
        'inputs':job['exact_inputs_and_hashes'],'reviewer_source_bytes':source.decode(),'reviewer_sha256':hashlib.sha256(source).hexdigest(),
        'seconds':time.monotonic()-begun,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    with (ROOT/OUTPUT).open('x') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
        f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_independent_checks_passed','rank_QQ','full_kernel_dimension','nonzero_minor_determinant','seconds','maximum_RSS_KiB']}))


if __name__ == '__main__':
    main()
