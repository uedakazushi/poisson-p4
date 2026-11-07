#!/usr/bin/env python3
"""Independent whole-parameter Fraction audit of cusp tangent/secant parents.

Uses only a pre-existing independent sparse Laurent/Fraction ring.  No producer,
S61 arithmetic engine, symbolic algebra system, sampling or field extension.
"""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

from classification_minimal_quartic_review import LaurentRing, determinant

ROOT = Path(__file__).resolve().parents[1]
PAIRS = list(combinations(range(5), 2))
QUADS = list(combinations(range(5), 4))


def matrix(R, entries):
    return [[R.parse(str(v)) for v in row] for row in entries]


def eye(R):
    return [[R.constant(i == j) for j in range(5)] for i in range(5)]


def mm(R, A, B):
    return [[R.add(*(R.mul(a, b) for a, b in zip(row, col)))
             for col in zip(*B)] for row in A]


def skew(R, entries):
    P = [[{} for _ in range(5)] for _ in range(5)]
    for (i, j), p in entries.items():
        P[i][j], P[j][i] = p, R.scale(p, -1)
    return P


def decode(R, data, tuple_keys=False):
    return skew(R, {tuple(map(int, k.strip('()').split(','))): R.parse(v)
                    for k, v in data.items()})


def sub(R, p, values):
    out = {}
    for e, c in p.items():
        term = R.constant(c)
        for name, power in zip(R.names, e):
            a = values.get(name, R.variable(name))
            if power < 0:
                a = R.inverse_monomial(a)
            for _ in range(abs(power)):
                term = R.mul(term, a)
        out = R.add(out, term)
    return out


def push(R, P, G, Gi):
    x = [R.variable(n) for n in R.names[R.n:]]
    values = {n: R.add(*(R.mul(a, b) for a, b in zip(row, x)))
              for n, row in zip(R.names[R.n:], Gi)}
    S = [[sub(R, p, values) for p in row] for row in P]
    return mm(R, mm(R, G, S), list(map(list, zip(*G))))


def jacobi(R, P):
    strict = {}
    for i, j, k in combinations(range(5), 3):
        strict[i, j, k] = R.add(*(R.mul(P[a][l], R.derivative(P[b][c], l))
            for a, b, c in [(i, j, k), (j, k, i), (k, i, j)] for l in range(5)))
    projective = {q: R.add(*(R.scale(R.mul(R.variable(R.names[R.n + q[k]]),
                      strict[q[:k] + q[k+1:]]), (-1)**k) for k in range(4))) for q in QUADS}
    return strict, projective


def actual_D(R, P):
    out = {}
    for i in range(5):
        a, b, c, d = [j for j in range(5) if j != i]
        pf = R.add(R.mul(P[a][b], P[c][d]), R.scale(R.mul(P[a][c], P[b][d]), -1),
                   R.mul(P[a][d], P[b][c]))
        out = R.add(out, R.scale(R.mul(R.variable(R.names[R.n+i]), pf), (-1)**i))
    return out


def same_projective(R, P, Q):
    for i, j, k in combinations(range(5), 3):
        value = R.add(*(R.scale(R.mul(R.variable(R.names[R.n+a]),
                         R.add(P[b][c], R.scale(Q[b][c], -1))), sign)
            for a, b, c, sign in [(i, j, k, 1), (j, i, k, -1), (k, i, j, 1)]))
        assert not value


def normalizer(source, saved):
    R = LaurentRing(tuple(f'r{i}' for i in range(6)), ('h', 'x', 'y', 'z', 'w'), invert=('r5',))
    P = decode(R, source)
    assert P == skew(R, {(1,3):R.parse('-h*x'),(1,4):R.parse('3*h*x-3*x**2'),
        (2,3):R.parse('3*(h-x)**2-3*h*y'),(2,4):R.parse('-9*x*y'),
        (3,4):R.parse('3*x*z+3*(r0*h**2+r1*h*x+r2*h*y+r3*x**2+r4*x*y+r5*y**2)')})
    assert not any(jacobi(R, P)[0].values())
    assert not any(jacobi(R, P)[1].values())
    assert actual_D(R, P) == R.parse('9*h*x*((h-x)**3-h**2*y)')
    G = eye(R)
    ell = R.parse('(r3+r2)*x/2+r4*y/4')
    kap = R.parse('(-3*r1-3*r3/2+9*r2/2)*x-r2*y')
    for i, p in [(3, ell), (4, kap)]:
        for j in range(5):
            G[i][j] = R.add(G[i][j], R.derivative(p, j))
    Gi = [[R.add(R.scale(G[i][j], -1), R.constant(2*(i == j))) for j in range(5)] for i in range(5)]
    assert mm(R, G, Gi) == eye(R) == mm(R, Gi, G)
    assert determinant(R, G) == R.constant(1)
    assert G == matrix(R, saved['magnetic_translation_matrix'])
    assert Gi == matrix(R, saved['magnetic_translation_inverse'])
    translated = push(R, P, G, Gi)
    vals = {f'r{i}': R.parse(v) for i, v in enumerate(['r0+r2', '0', '0', '0', '0', 'r5'])}
    assert translated == [[sub(R, p, vals) for p in row] for row in P]
    scale = eye(R); inv = eye(R)
    scale[3][3] = scale[4][4] = R.parse('1/r5')
    inv[3][3] = inv[4][4] = R.parse('r5')
    assert mm(R, scale, inv) == eye(R)
    normalized = [[R.mul(R.parse('r5'), p) for p in row] for row in push(R, translated, scale, inv)]
    vals = {f'r{i}': R.parse(v) for i, v in enumerate(['(r0+r2)/r5', '0', '0', '0', '0', '1'])}
    assert normalized == [[sub(R, p, vals) for p in row] for row in P]
    # All six, including I=0, retain the same actual modular pair and failure.
    Z = [R.add(*(R.derivative(P[i][j], i) for i in range(1, 5))) for j in range(5)]
    assert Z == [R.parse(v) for v in ['0', '0', '0', '-4*h', '3*h-12*x']]
    X = []
    for p in [P[i][1] for i in range(5)]:
        out = {}
        for e, c in p.items():
            assert e[R.n+1] > 0
            d = list(e); d[R.n+1] -= 1; out[tuple(d)] = c
        X.append(out)
    assert X == [R.parse(v) for v in ['0', '0', '0', 'h', '-3*h+3*x']]
    B = [R.add(a, b) for a, b in zip(Z, X)]
    fail = R.add(R.mul(R.parse('h'), B[3]), R.scale(R.mul(R.parse('z'), B[0]), -1))
    assert sub(R, fail, {'x': {}}) == R.parse('-3*h**2')
    return dict(entire_six_magnetic_translation_and_inverse=True,
                only_inverted_parameter='I=r5', full_nonzero_I_scaling_inverse=True,
                actual_D='9*h*x*((h-x)^3-h^2*y)', all_six_T21_field_identity=True,
                T21_condition_a_on_x_zero_failure='-3*h^2', I_zero_not_excluded_from_C39_C40=True)


def whole_graphs(saved, cert, model, equations, sections):
    R = LaurentRing(('lam',), ('h', 'x', 'y', 'z', 'w'))
    P = decode(R, saved['whole_lambda_brackets'])
    assert P == skew(R, {(1,3):R.parse('-h*x'),(1,4):R.parse('3*h*x-3*x**2'),
        (2,3):R.parse('3*(h-x)**2-3*h*y'),(2,4):R.parse('-9*x*y'),
        (3,4):R.parse('3*x*z+3*lam*h**2+3*y**2')})
    strict, J = jacobi(R, P)
    assert not any(strict.values()) and not any(J.values())
    move, frame = matrix(R, saved['affine_translation']), matrix(R, saved['Darboux_frame'])
    point = dict(zip(('h','x','y','z','w'), map(R.constant, (1,1,1,0,0))))
    K = [[sub(R, P[i][j], point) for j in range(1,5)] for i in range(1,5)]
    pf = R.add(R.mul(K[0][1], K[2][3]), R.scale(R.mul(K[0][2], K[1][3]), -1), R.mul(K[0][3], K[1][2]))
    assert pf == R.constant(-9)
    F4 = [row[1:] for row in frame[1:]]
    standard = [[R.constant(((i,j) in [(0,1),(2,3)])-((i,j) in [(1,0),(3,2)])) for j in range(4)] for i in range(4)]
    assert mm(R, mm(R, F4, K), list(map(list, zip(*F4)))) == standard
    perspective = eye(R);perspective[0][1:] = [R.parse(v) for v in saved['perspective']]
    M, Mi = matrix(R, saved['complete_PGL_matrix']), matrix(R, saved['complete_PGL_inverse'])
    assert M == mm(R, perspective, mm(R, frame, move))
    assert mm(R, M, Mi) == eye(R) == mm(R, Mi, M)
    assert determinant(R, M) == R.constant(F(-1,9))
    pushed = push(R, P, M, Mi)
    params = [R.parse(v) for v in saved['all61_parameters']]
    assert len(params) == 61 and all(not any(e[1:]) for p in params for e in p)
    def ev(data):
        return R.add(*(R.mul(R.constant(F(num,den)), *(params[k] for k in mon)) for mon,num,den in data))
    mons=[]
    for terms in combinations_with_replacement(range(5),2):
        e=[0]*5
        for j in terms:e[j]+=1
        mons.append(tuple(e))
    cols=[(ij,e) for ij in PAIRS for e in mons]
    v=[{} for _ in range(150)]
    for i,(ij,e) in enumerate(cols):
        if ij in [(1,2),(3,4)] and e == (2,0,0,0,0):v[i]=R.constant(1)
    for i,p in zip(cert['raw_free_indices'],params):v[i]=p
    for a in model['eliminated_raw_variables']:v[a['raw_index']]=ev(a['polynomial'])
    entries={}
    for (ij,e),p in zip(cols,v):
        mon={((0,)+e):F(1)}
        entries[ij]=R.add(entries.get(ij,{}),R.mul(p,mon))
    represented=skew(R,entries)
    same_projective(R,pushed,represented)
    same_projective(R,push(R,represented,Mi,M),P)
    assert not any(jacobi(R,represented)[1].values())
    assert len(equations['equations'])==134 and not any(ev(p) for p in equations['equations'])
    values={f'u{i+1}':p for i,p in enumerate(params)}
    # Extend names solely for exact AST polynomial substitution in graph formulas.
    T=LaurentRing(tuple(values),())
    out={}
    for label, section in sections.items():
        def graph_ev(s):
            ans={}
            for e,c in T.parse(s).items():
                term=R.constant(c)
                for name,power in zip(T.names,e):
                    for _ in range(power):term=R.mul(term,values[name])
                ans=R.add(ans,term)
            return ans
        residual=[R.add(params[k],R.scale(graph_ev(s),-1)) for k,s in zip(section['weight2_u_indices_zero_based'],section['quadratic_section_q'])]
        assert len(residual)==45
        assert residual == [R.parse(s) for s in saved['graph_membership'][label]['all45_residuals']]
        constants=[p[R.zero] for p in residual if p and set(p)=={R.zero}]
        witness=F(1,4) if label=='T39' else F(-2,3)
        assert witness in constants
        assert sum(bool(p) for p in residual)==(11 if label=='T39' else 13)
        out[label]=dict(all45_residuals=[R.encode(p) for p in residual],
                        nonzero_count=sum(bool(p) for p in residual),
                        exact_unit_ideal_witness=str(witness), all_finite_lambda_excluded=True)
    return dict(entire_QQ_lambda_not_sampled=True, uniform_Pf='-9', PGL_det='-1/9',
                all150_forward_and_inverse_Euler_quotients=True, strict10_zero=True,
                all350_projective_Jacobi_zero=True, all134_equations_zero=True,
                constant_parameter_denominators_only=True, graphs=out)


def T21_fields(raw):
    R=LaurentRing(('p0','p1','t0','t1','t2'),('A','B','C','H','W'))
    P=decode(R,{k:v.replace('x0','A').replace('x1','B').replace('x2','C').replace('x3','H').replace('x4','W') for k,v in raw['homogeneous_bivector'].items()})
    L=[R.parse(v) for v in ['0','-p1*B/3','p1*C/3','0','p0*W']]
    E=[R.variable(v) for v in ('A','B','C','H','W')]
    Q=[[R.add(P[i][j],R.scale(R.mul(E[i],L[j]),-1),R.mul(E[j],L[i])) for j in range(5)] for i in range(5)]
    f=R.parse('(p1*A*B*C+t0*C**3-t1*B**3+t2*A**3)/3')
    entries={(0,1):R.derivative(f,2),(0,2):R.scale(R.derivative(f,1),-1),(1,2):R.derivative(f,0)}
    entries.update({(i,4):R.parse('-p0*W*'+name) for i,name in enumerate(('A','B','C'))})
    assert Q==skew(R,entries)
    assert not any(jacobi(R,P)[1].values())
    Z=[R.add(*(R.derivative(Q[i][j],i) for i in (0,1,2,4))) for j in range(5)]
    assert Z==[R.parse(v) for v in ['p0*A','p0*B','p0*C','0','-3*p0*W']]
    X=[R.parse(v) for v in ['-p0*A','-p0*B','-p0*C','0','0']]
    assert [R.mul(R.parse('W'),v) for v in X]==[Q[i][4] for i in range(5)]
    assert [R.add(a,b) for a,b in zip(Z,X)]==[R.parse(v) for v in ['0','0','0','0','-3*p0*W']]
    assert [R.add(a,R.scale(b,4)) for a,b in zip(Z,X)]==[R.parse(v) for v in ['-3*p0*A','-3*p0*B','-3*p0*C','0','-3*p0*W']]
    return dict(entire_five_parameter_raw_source_and_Euler_gauge=True,
                no_period_or_edge_inversions=True, both_intrinsic_flag_conditions_identically=True,
                proper_flag_and_swap_scope='separately checked by geometric proof in review note')


def main():
    signal.alarm(110)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begin=time.monotonic()
    job=json.loads((ROOT/'research/classification/jobs/two_lines_global_parent_separate_review.json').read_text())
    for a in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/a['path']).read_bytes()).hexdigest()==a['sha256'],a['path']
    def load(p):return json.loads((ROOT/p).read_text())
    source=next(c for c in load('reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json')['cases'] if c['case']=='cusp_cusp_tangent_cusp_secant')
    saved=load('reports/CLASSIFICATION_TWO_LINES_CUSP_TANGENT_SECANT_WHOLE_GRAPHS_20261107.json')
    normal=normalizer(source['full_six_magnetic_source'],saved)
    graph=whole_graphs(saved,load('legacy/v02/data/symplectic_slice_certificate.json'),load('legacy/v02/data/symplectic_slice_equations.json'),load('legacy/v02/data/symplectic_slice_reduced_equations.json'),{a:load(p) for a,p in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]})
    t21=T21_fields(load('reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T21.json'))
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,method='Independent standard-library Fraction sparse Laurent polynomials; no producer, S61 engine or CAS imports.',normalization=normal,whole_graphs=graph,T21=t21,input_sha256=job['exact_inputs_and_hashes'],reviewer_source=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-begin,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,no_new_component_or_other_parent_exclusion=True,geometric_scope='C39/C40 use the separately accepted global evaluation fiber theorem; T21 uses all five-parameter source plus independent proper flag and order swap proof.')
    out=ROOT/'reports/CLASSIFICATION_TWO_LINES_GLOBAL_PARENT_SEPARATE_REVIEW_20261107.json'
    with out.open('x') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print('PASS whole lambda + whole T21:',round(result['seconds'],3),'seconds',flush=True)


if __name__=='__main__':main()
