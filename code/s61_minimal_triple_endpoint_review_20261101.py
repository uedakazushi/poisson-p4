#!/usr/bin/env python3
"""Independent raw135 audit of the triple-tail binary endpoint."""
import argparse
import hashlib
import json
import time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61, X, J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, sp_basis, infinitesimal_action

ROOT = Path(__file__).resolve().parents[1]
M = lambda a: s.Matrix(a).applyfunc(s.sympify)


def run():
    start = time.time()
    rp = ROOT/'reports/S61_MINIMAL_TRIPLE_ENDPOINT_20261101.json'
    r = json.loads(rp.read_text())
    assert r['all_passed']
    for name, digest in r['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    model = S61()
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 2]
    qi, ui = {i:j for j,i in enumerate(two)}, {i:j for j,i in enumerate(one)}
    equations = model.reduced['equations'] + [json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    encoded = [[(tuple(qi[i] for i in mon if i in qi), tuple(ui[i] for i in mon if i in ui), s.Rational(a,b)) for mon,a,b in eq] for eq in equations]
    ut, geo, local = r['universal_thickening'], r['source_geometry'], r['point_local_data']
    K,L,V,U,d,z = s.symbols('K L V U d s')
    nuisance = s.symbols('n0:6')
    beta, delta, q0 = M(ut['beta']), M(ut['delta']), M(ut['q0'])
    nn, ww = M(ut['n']), M(ut['w'])
    dn, dw = s.sympify(ut['n_den']), s.sympify(ut['w_den'])
    assert dn == 80*(-23*L+10*V*z)
    assert s.expand(dw-19200*V*(-23*L+10*V*z)**3) == 0
    R, *_ = ring([K,L,V]+list(nuisance)+[z], QQ)
    q = [R.from_expr(v) for v in q0+z*delta]
    n = [R.from_expr(s.cancel(v*dn)) for v in nn]
    w = [R.from_expr(s.cancel(v*dw)) for v in ww]
    DN, DW = R.from_expr(dn), R.from_expr(dw)
    coefficients = {}
    for index, eq in enumerate(encoded):
        residual = [R.zero, R.zero, R.zero]
        for qm,bm,c in eq:
            if len(qm) != 2:
                continue
            if bm not in coefficients:
                coefficients[bm] = R.from_expr(s.prod(beta[j] for j in bm))
            a,b = qm
            C = c*coefficients[bm]
            residual[0] += C*q[a]*q[b]
            residual[1] += C*(q[a]*n[b]+n[a]*q[b])
            residual[2] += C*((q[a]*w[b]+w[a]*q[b])*DN**2+n[a]*n[b]*DW)
        assert all(value == 0 for value in residual), index
    fixed = two.index(29)
    assert q0[fixed] == 1 and nn[fixed] == ww[fixed] == delta[fixed] == 0
    n0 = nn.subs(z,0)
    assert (n0-M(ut['n0'])).applyfunc(s.cancel) == s.zeros(45,1)
    assert s.factor(delta.row_join(n0)[[0,24],:].det()) == -31*L/920

    # Rebuild the local derivative and quadratic obstruction from all raw equations.
    point = M(local['point61'])
    b, qq = point[one,:], point[two,:]
    assert qq == q0
    chart = [i for i in range(45) if i != fixed]
    derivative, tcol, polys, cache = s.zeros(135,45), s.zeros(135,1), [], {}
    for row, eq in enumerate(encoded):
        f = {}
        for qm,bm,c in eq:
            if bm not in cache:
                cache[bm] = s.prod(b[k] for k in bm)
            co = c*cache[bm]
            f[qm] = f.get(qm,0)+co
            if len(qm) == 2:
                i,j = qm
                derivative[row,i] += co*qq[j]
                derivative[row,j] += co*qq[i]
            elif len(qm) == 1:
                tcol[row] += co*qq[qm[0]]
        assert sum(c*qq[m[0]]*qq[m[1]] for m,c in f.items() if len(m)==2) == 0
        polys.append(f)
    A = derivative[:,chart].row_join(tcol)
    T, annihilator = M(local['tangent_basis']), M(local['cokernel_rows'])
    assert A == M(local['jacobian']) and A*T == s.zeros(135,10)
    assert annihilator*A == s.zeros(annihilator.rows,45) and T[-1,:] == s.zeros(1,10)
    forms = {i:([0]*10 if i==fixed else list(T[chart.index(i),:])) for i in range(45)}
    monos = list(combinations_with_replacement(range(10),2))
    indices = {m:i for i,m in enumerate(monos)}
    Q = s.zeros(135,len(monos))
    for row,f in enumerate(polys):
        for qm,c in f.items():
            fs = [forms[i] for i in qm]+[list(T[-1,:])]*(2-len(qm))
            for i in range(10):
                for j in range(10):
                    Q[row,indices[tuple(sorted((i,j)))]] += c*fs[0][i]*fs[1][j]
    obstruction = annihilator*Q
    assert obstruction == M(local['quadratic_obstruction_basis'])
    hessians = []
    for coeff in obstruction.tolist():
        hs = s.zeros(10)
        for c,(i,j) in zip(coeff,monos):
            hs[i,j] += c
            hs[j,i] += c
        hessians += hs.tolist()
    ranks = (derivative.rank(), A.rank(), obstruction.rank(), s.Matrix(hessians).rank())
    assert ranks == (34,35,5,7)
    free = [i for i in range(45) if i not in A.rref()[1]]
    assert T[free,:] == s.eye(10)
    source_beta = M(geo['beta'])
    source_variables = list(map(s.Symbol,geo['parameters']))
    sample = {s.Symbol(k):s.Rational(v) for k,v in geo['beta_sample'].items()}
    assert source_beta.subs(sample) == b
    recovery = {s.Symbol(k):s.sympify(v) for k,v in geo['source_recovery'].items()}
    scale = s.sympify(geo['source_beta_scale'])
    assert (scale*beta.subs(recovery,simultaneous=True)-source_beta).applyfunc(s.cancel) == s.zeros(16,1)
    env = {k:v.subs(sample) for k,v in recovery.items()}
    dt = delta.subs(env)[chart,:].col_join(s.zeros(1,1))
    nt = n0.subs(env)[chart,:].col_join(s.zeros(1,1))
    tangents = dt.row_join(nt)[free,:]
    assert T*tangents == dt.row_join(nt)
    assert tangents == M(r['actual_triple_line_tangent'])
    assert tangents[3,0] == 0 and tangents[3,1] != 0
    assert M(r['actual_triple_line_initial_annihilator'])*tangents == s.zeros(8,2)
    assert M(r['actual_triple_line_initial_annihilator']).rank() == 8
    for vector in [nn,ww]:
        for value in vector:
            numerator, denominator = s.fraction(s.cancel(value.subs(env)))
            assert int(s.Rational(denominator.subs(z,0)).p) % 32003
            assert all(int(s.Rational(c).q) % 32003 for c in s.Poly(numerator,z).all_coeffs()+s.Poly(denominator,z).all_coeffs())

    # Independently identify the same smooth binary threefold.
    zz,u,bf = primitive_beta(model)
    a,c,h = s.symbols('a c h')
    v1, v2 = s.Matrix([a,1,0,c]), s.Matrix([-c,0,1,h])
    bp = bf.subs(dict(zip(u,b)))
    section = s.Matrix([(v1.T*bp.subs(dict(zip(zz,v)))*v2)[0] for v in [v1,v2]])
    CJ = section.jacobian([a,c,h]).subs({a:0,c:0,h:0})
    assert CJ == M(r['curve_jacobian']) and CJ.rank() == 2
    family = M(r['actual_binary_family_tangent'])
    assert s.Matrix(hessians)*family == s.zeros(50,3)
    assert family[r['binary_parameters'],:].det() != 0

    # Reconstruct beta action and the full original U-dependent endpoint.
    vec = s.Matrix([bf[i,j].diff(zz[k]) for (i,j),k in COORDS])
    basis = vec.jacobian(u)
    rows = list(basis.T.rref()[1])
    iv = basis[rows,:].inv()
    action = s.Matrix.hstack(*[iv*infinitesimal_action(vec,B)[rows,:] for B in sp_basis()])
    spatial = s.Matrix(X[1:]); x,y,z0,w0 = spatial
    baseQ = s.zeros(4); baseQ[1,2] = x*w0; baseQ -= baseQ.T
    def coordinates(Q):
        P = X[0]**2*J
        P[1:,1:] += Q
        return s.Matrix(model.normalize(P)['parameters'])[two,:]
    orbit = []
    for B in sp_basis():
        dq = coordinates(sum(((B*spatial)[i]*baseQ.diff(spatial[i]) for i in range(4)),s.zeros(4))-B*baseQ-baseQ*B.T)
        orbit.append(dq-q0*dq[fixed])
    O = s.Matrix.hstack(*orbit)
    stabilizer = M(geo['stabilizer'])
    assert O.rank() == 5 and stabilizer.rank() == 5 and O*stabilizer == s.zeros(45,5)
    H = source_beta.jacobian(source_variables)
    ann = M(geo['source_annihilator'])
    assert ann.rank() == 6 and ann*H == s.zeros(6,10)
    ar = action.subs(dict(zip(u,source_beta)))
    normal = (ann*ar*stabilizer).applyfunc(s.expand)
    assert (normal-M(geo['normal_stabilizer_matrix'])).applyfunc(s.expand) == s.zeros(6,5)
    assert normal.to_DM().rank() == 1
    differential = H.row_join(ar)
    joint = differential.col_join(s.zeros(45,10).row_join(O))
    assert differential.subs(sample).rank() == joint.subs(sample).rank() == 16
    assert differential.subs(sample)[:,geo['beta_minor_columns']].det() == s.Rational(geo['beta_minor'])
    old = json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text())
    fam = next(f for f in old['families'] if f['label']=='M1')
    tail = next(t for t in fam['rational_tails'] if t['axis']==2 and t['k']=='-1/3')
    gg = s.symbols('g0:9')
    monomials = [x*x*z0,x*x*w0,x*z0*z0,x*z0*w0,x*w0*w0,z0**3,z0*z0*w0,z0*w0*w0,w0**3]
    potential = sum(s.sympify(tail['coordinate_formulas'].get(str(g),str(g)))*m for g,m in zip(gg,monomials)).subs(gg[2],d)
    highest = s.diff(potential,d)
    assert s.expand(highest-x*z0*z0-(-3*L+14*U)*x*x*z0) == 0
    Pinf = s.zeros(4); Pinf[1,2] = -s.diff(highest,w0); Pinf[1,3] = s.diff(highest,z0); Pinf -= Pinf.T
    transform = M(geo['full_source_normalization'])
    assert transform*J[1:,1:]*transform.T == J[1:,1:]
    assert (transform.inv()*Pinf.subs(dict(zip(spatial,transform*spatial)),simultaneous=True)*transform.inv().T+2*baseQ).applyfunc(s.expand) == s.zeros(4)
    fullsource = s.Matrix([0,0,0,0,-s.Rational(1,3),nuisance[0],U,nuisance[1],nuisance[2],V,nuisance[3],nuisance[4],nuisance[5],s.Rational(1,3),K-5*V,5*U-L])
    Bb = bf.subs(dict(zip(u,fullsource)))
    pullback = transform.T*Bb.subs(dict(zip(zz,transform*s.Matrix(zz))),simultaneous=True)*transform
    pulled = (iv*s.Matrix([pullback[i,j].diff(zz[k]) for (i,j),k in COORDS])[rows,:]).applyfunc(s.expand)
    assert (pulled-M(geo['full_source_normalized_beta'])).applyfunc(s.expand) == s.zeros(16,1)
    assert ann*pulled == s.zeros(6,1)
    paths = [Path(__file__),rp,ROOT/'code/s61_minimal_triple_endpoint_20261101.py',ROOT/'research/notes/S61_MINIMAL_TRIPLE_ENDPOINT_20261101.md',ROOT/'research/notes/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.md',ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'code/truncated_local_ideal_20261101.cpp',ROOT/'code/truncated_local_ideal_io_20261101.py']
    return dict(id='S61_MINIMAL_TRIPLE_ENDPOINT_REVIEW_20261101',all_passed=True,raw135_all405_triple_line_coefficients=True,regular_flat_triple_line_and_embedding=True,raw135_point_ranks=list(ranks),same_smooth_binary_threefold=True,full_original_U_endpoint_coverage=True,projective_orbit_rank=5,normal_stabilizer_rank=1,beta_rank=16,joint_rank=16,conceptual_review=dict(strict_filtered_kernel_accepted=True,flat_nonreduced_quotient_sufficient=True,DVR_equal_rank_torsion_free_accepted=True,origin_factor_completion_accepted=True,relative_ambient_t2_and_generic16D_removal_accepted=True),scope=['Independent QQ raw135 universal triple-line identities, local derivatives and quadratic obstructions, same smooth binary chart, source action and full U endpoint coverage.','Primary degree7 finite-field closure and Macaulay2 decomposition are separately reproduced; no independent implementation of that elimination.','The kernel is taken over the flat triple line, not over its reduction.','No free66 assertion, no remaining flag/binary exhaustion or affine exhaustion.'],seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
