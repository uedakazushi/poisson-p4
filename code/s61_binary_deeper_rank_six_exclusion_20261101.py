#!/usr/bin/env python3
"""B-rank at most six on K=Qstar=0 has proper beta image, by divisor minors."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sympy as s
ROOT = Path(__file__).resolve().parents[1]


def verify():
    start = time.time()
    cover_path = ROOT/'reports/S61_BINARY_RANK_DROP_COVER_20261101.json'
    rank_path = ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    previous_path = ROOT/'reports/S61_BINARY_DEEPER_REDUCTION_20261101.json'
    cover = json.loads(cover_path.read_text())
    rank = json.loads(rank_path.read_text())
    previous = json.loads(previous_path.read_text())
    assert cover['all_passed'] and rank['all_passed'] and previous['all_passed']
    assert previous['excluded_components'][0]['ideal_index'] == 0
    assert previous['excluded_components'][0]['beta_image_dimension_upper_bound'] == 15
    B = s.Matrix(cover['block_matrices'][1]).applyfunc(s.sympify)
    u = s.symbols('u1:62')
    a,b,c,d,e,f = [u[i] for i in [3,13,34,35,48,49]]
    core = [a,b,c,d,e,f]
    K,L = 5*c+f,5*b-e
    Q = 600*b*b-255*b*e-63*e*e+100*a*f
    param = {c:-f/5, a:-(600*b*b-255*b*e-63*e*e)/(100*f)}
    restricted = B.subs(param, simultaneous=True).applyfunc(s.cancel)
    rows = [0,1,2,3,5,6,10]
    cols = [0,1,2,3,4,5,7]
    minor = s.factor(restricted[rows,cols].det(method='domain-ge'))
    A = 100*b*b-505*b*e+187*e*e
    R = d*A-40*f*f*(2*b-e)
    expected = s.Rational(27,12500000)*(5*b-11*e)*(5*b-3*e)*L**2*R
    assert s.expand(minor-expected) == 0
    # R is primitive linear in d even after base extension: gcd(A,2b-e)=1.
    assert s.expand(A.subs(b,e/2)+s.Rational(81,2)*e**2) == 0
    assert A.subs(e,0) == 100*b*b
    assert s.expand(Q.subs(f,0)-3*(5*b-3*e)*(40*b+7*e)) == 0

    test_parameters = [
        ('5b-11e',s.Rational(11,5),s.Integer(1),s.Integer(3),s.Integer(2)),
        ('5b-3e',s.Rational(3,5),s.Integer(1),s.Integer(3),s.Integer(2)),
        ('R',s.Integer(2),s.Integer(3),s.Integer(5),-s.Rational(1000,947)),
    ]
    results = []
    for label,bv,ev,fv,dv in test_parameters:
        av = -(600*bv*bv-255*bv*ev-63*ev*ev)/(100*fv)
        sample = {a:av,b:bv,c:-fv/5,d:dv,e:ev,f:fv}
        assert K.subs(sample)==0 and Q.subs(sample)==0 and L.subs(sample)!=0
        factor = {'5b-11e':5*b-11*e,'5b-3e':5*b-3*e,'R':R}[label]
        assert factor.subs(sample)==0
        matrix = B.subs(sample)
        cp = list(matrix.to_DM().rref()[1])
        rp = list(matrix[:,cp].T.to_DM().rref()[1])
        assert len(cp)==7
        determinant=matrix[rp,cp].det()
        assert determinant!=0
        results.append(dict(divisor=label,sample={str(k):str(v) for k,v in sample.items()},
            rows=rp,columns=cp,minor=str(determinant)))
    for label,bv in [('f=0,5b-3e=0',s.Rational(3,5)),('f=0,40b+7e=0',-s.Rational(7,40))]:
        sample={a:s.Integer(2),b:bv,c:s.Integer(0),d:s.Integer(3),e:s.Integer(1),f:s.Integer(0)}
        assert K.subs(sample)==0 and Q.subs(sample)==0 and L.subs(sample)!=0
        matrix=B.subs(sample)
        cp=list(matrix.to_DM().rref()[1]);rp=list(matrix[:,cp].T.to_DM().rref()[1])
        assert len(cp)==7
        determinant=matrix[rp,cp].det();assert determinant!=0
        results.append(dict(divisor=label,sample={str(k):str(v) for k,v in sample.items()},
            rows=rp,columns=cp,minor=str(determinant)))

    # All five fixed-Q stabilizer directions preserve the rank condition.
    one=[3,4,11,12,13,24,25,26,33,34,35,36,47,48,49,50]
    action=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    restricted_action=action.subs({u[11]:0,u[12]:0})
    stabilizer=s.Matrix.hstack(*[s.Matrix(H).applyfunc(s.sympify).reshape(16,1) for H in rank['stabilizer_generators']])
    assert stabilizer.rank()==5
    assert all(restricted_action[one.index(i),j]==0 for i in [11,12] for j in range(5))
    core_action=restricted_action[[one.index(i) for i in [3,13,34,35,48,49]],:]
    assert core_action[:,[0,1,3]]==s.zeros(6,3)
    weights=[(-3,-1,1,3,-1,1),(2,1,0,-1,1,0)]
    assert core_action[:,2]==s.Matrix([w*v for w,v in zip(weights[0],core)])
    assert 2*core_action[:,4]==s.Matrix([w*v for w,v in zip(weights[1],core)])
    rw=s.symbols('r0:12');cw=s.symbols('c0:9');weight_certificates=[]
    for w in weights:
        equations=[]
        for i in range(12):
            for j in range(9):
                if B[i,j]:
                    observed={sum(mi*wi for mi,wi in zip(mon,w)) for mon,_ in s.Poly(B[i,j],*core).terms()}
                    assert len(observed)==1
                    equations.append(rw[i]+cw[j]-observed.pop())
        sol=list(s.linsolve(equations,rw+cw))[0]
        free=set().union(*(v.free_symbols for v in sol))
        sol=[v.subs(dict.fromkeys(free,s.S.Zero)) for v in sol]
        row_weights,column_weights=sol[:12],sol[12:]
        for i in range(12):
            for j in range(9):
                deriv=sum(wi*vi*s.diff(B[i,j],vi) for wi,vi in zip(w,core))
                assert s.expand(deriv-(row_weights[i]+column_weights[j])*B[i,j])==0
        for polynomial in [K,L,Q]:
            deriv=sum(wi*vi*s.diff(polynomial,vi) for wi,vi in zip(w,core))
            assert not s.cancel(deriv/polynomial).free_symbols
        weight_certificates.append(dict(core_weights=w,row_weights=list(map(str,row_weights)),column_weights=list(map(str,column_weights))))
    paths=[Path(__file__),cover_path,rank_path,previous_path,
        ROOT/'research/notes/S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_20261101.md']
    return dict(id='S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_20261101',all_passed=True,field='QQ',
        source_equations=list(map(str,[K,Q])),rational_chart_parameterization={str(k):str(v) for k,v in param.items()},
        generic_minor_rows=rows,generic_minor_columns=cols,generic_seven_minor=str(minor),
        candidate_divisors_away_from_f_zero=list(map(str,[5*b-11*e,5*b-3*e,L,R])),
        nonzero_minor_on_each_remaining_divisor=results,
        rank_condition_stabilizer_invariant=True,stabilizer_core_torus_weights=weight_certificates,
        core_rank_defect_dimension_bound_away_from_L=2,
        unrestricted_other_beta_coordinates=8,stabilizer_dimension=5,
        beta_image_dimension_upper_bound=15,
        all_K_zero_B_rank_at_most_six_avoided_by_general_beta=True,
        proof_note='research/notes/S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_20261101.md',
        limitations=['This is a proper-beta-image exclusion for the B-block rank<=6 subset.',
            'It does not determine all C-block defects or the full free65 open complement.',
            'The surviving K=Qstar=0 stratum with B-rank7 is not exhausted.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','beta_image_dimension_upper_bound','seconds']}))
