#!/usr/bin/env python3
"""Exact zero-first-jet CA/CB fibers and their established Q45 labels.

All computations are over QQ. The separate proof uses the already accepted
proper line incidences; no full S61 zero-fiber decomposition is asserted.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

from weight1_structure import J,PAIRS,TRIPLES,sp_basis,minor_witness

ROOT=Path(__file__).resolve().parents[1]
Y=s.Matrix(s.symbols('x y z w'));x,y,z,w=Y
MON=[s.prod(Y[i] for i in ids) for ids in itertools.combinations_with_replacement(range(4),2)]

def flatten(P):
    return s.Matrix([s.Poly(P[i,j],*Y).coeff_monomial(m)
                     for i,j in PAIRS for m in MON])

def skew(entries):
    Q=s.zeros(4)
    for a,(i,j) in zip(entries,PAIRS):Q[i,j]=a;Q[j,i]=-a
    return Q

def divergence(Q):
    return s.Matrix([s.expand(sum(s.diff(Q[i,j],Y[i]) for i in range(4))) for j in range(4)])

def schouten(P,Q):
    return [s.expand(sum(P[i,h]*s.diff(Q[j,k],Y[h])+P[j,h]*s.diff(Q[k,i],Y[h])+P[k,h]*s.diff(Q[i,j],Y[h])
        +Q[i,h]*s.diff(P[j,k],Y[h])+Q[j,h]*s.diff(P[k,i],Y[h])+Q[k,h]*s.diff(P[i,j],Y[h]) for h in range(4)))
        for i,j,k in TRIPLES]

def wedge(v,Q):
    return [s.expand(v[i]*Q[j,k]-v[j]*Q[i,k]+v[k]*Q[i,j]) for i,j,k in TRIPLES]

def equations_coefficients(polynomials):
    result=[]
    for f in polynomials:result.extend(s.Poly(f,*Y).coeffs())
    return result

def factor_matrix(Q):
    v=s.Matrix(s.symbols('v0:4'))
    rows=[s.Poly(f,*Y).coeff_monomial(m) for f in wedge(v,Q) for m in MON]
    return s.linear_eq_to_matrix(rows,list(v))[0]

def linear_kernel_matrix(rows,variables):
    matrix,_=s.linear_eq_to_matrix(equations_coefficients(rows),variables)
    return matrix,s.Matrix.hstack(*matrix.nullspace())

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    r=s.Symbol('r');k=s.symbols('k0:9')
    cubics=[s.prod((x,z,w)[i] for i in ids) for ids in itertools.combinations_with_replacement(range(3),3) if ids!=(0,0,0)]
    G=sum(a*b for a,b in zip(k,cubics))
    QA=skew([r*x*x,0,0,-s.diff(G,w),s.diff(G,z),0])
    QB=skew([r*x*x,0,0,-r*x*z-s.diff(G,w),-r*x*w+s.diff(G,z),r*x*x])
    QR=QA.subs(r,0)
    families={'M0':QA,'M1':QB,'R0':QR};records={}
    for name,Q in families.items():
        assert schouten(J,Q)==[0]*4 and schouten(Q,Q)==[0]*4
        div=divergence(Q);expected=s.Matrix([0,(2 if name=='M0' else 4)*r*x,0,0]) if name!='R0' else s.zeros(4,1)
        assert div==expected
        pars=[*k] if name=='R0' else [r,*k]
        embedding=flatten(Q).jacobian(pars);assert embedding.rank()==len(pars)
        sample={r:1}|dict(zip(k,[2,3,5,7,11,13,17,19,23]))
        if name=='R0':sample={**sample,**dict.fromkeys(k,0),k[5]:1,k[8]:1}
        value=Q.subs(sample)
        # The three infinitesimal chart generators move e_y to e_x,e_z,e_w.
        chart_dirs=[]
        for dr,ds,dt in [(1,0,0),(0,1,0),(0,0,1)]:
            A=s.Matrix([[0,dr,-dt,ds],[0,0,0,0],[0,ds,0,0],[0,dt,0,0]])
            field=A*Y
            chart_dirs.append(flatten(A*value+value*A.T-sum((field[i]*value.diff(Y[i]) for i in range(4)),s.zeros(4))))
        differential=embedding.row_join(s.Matrix.hstack(*chart_dirs))
        assert differential.rank()==len(pars)+3
        F=factor_matrix(value)
        records[name]=dict(raw_coefficients=list(map(str,flatten(Q))),divergence=list(map(str,div)),
            fixed_line_dimension=len(pars),total_source_dimension=len(pars)+3,
            differential_witness=minor_witness(differential),generic_factor_dimension=4-F.rank())
        print(name,'EXACT DIMENSION',len(pars)+3,flush=True)

    # Independent fixed-line recognition in the full 60-coefficient raw space.
    unknown=s.symbols('q0:60');Q=skew([sum(unknown[10*a+i]*m for i,m in enumerate(MON)) for a in range(6)])
    compatibility=schouten(J,Q);Cmat,W=linear_kernel_matrix(compatibility,unknown)
    assert W.cols==45
    div=divergence(Q);trace=Q[0,1]+Q[2,3]
    # This identity follows from compatibility, not from the Poisson equations.
    relation=div+J*s.Matrix([s.diff(trace,v) for v in Y])
    relmat,_=s.linear_eq_to_matrix(equations_coefficients(list(relation)),unknown)
    assert relmat*W==s.zeros(relmat.rows,45)
    line=s.Matrix([0,1,0,0]);image_line=[div[i] for i in (0,2,3)]
    checks={}
    for name,target in [('M0',Q),('M1',Q-trace*J/2),('R0',Q)]:
        rows=compatibility+wedge(line,target)+(list(div) if name=='R0' else image_line)
        constraints,basis=linear_kernel_matrix(rows,unknown)
        model=families[name];pars=list(k) if name=='R0' else [r,*k]
        embedding=flatten(model).jacobian(pars)
        assert constraints*embedding==s.zeros(constraints.rows,len(pars))
        assert basis.cols==len(pars)
        checks[name]=dict(ambient_raw_dimension=60,linear_constraint_rank=constraints.rank(),kernel_dimension=basis.cols,
            full_matrix=constraints.tolist(),family_embedding=embedding.tolist())
    # Exact old M0/M1 table crosswalk; normalization is trace=x^2/2.
    t=s.symbols('t1:10');crosswalk={}
    for name,c,rv in [('M0',s.S.Zero,s.Rational(1,2)),('M1',s.Rational(1,4),s.Rational(1,4))]:
        old=skew([(s.Rational(1,2)-c)*x*x,0,0,
            t[0]*x*x+(-2*c-t[5])*x*z-t[7]*z*z/2-2*t[8]*z*w+t[1]*x*w+t[2]*w*w,
            t[3]*x*x+t[4]*x*z+t[5]*x*w+t[6]*z*z+t[7]*z*w+t[8]*w*w,c*x*x])
        substitution=dict(zip(t,[-k[1],-2*k[4],-3*k[8],k[0],2*k[2],k[3]-(rv if name=='M1' else 0),3*k[5],2*k[6],k[7]]))
        assert (old.subs(substitution,simultaneous=True)-families[name].subs(r,rv)).applyfunc(s.expand)==s.zeros(4)
        crosswalk[name]=dict(fixed_r=str(rv),fixed_trace='x^2/2',old_c=str(c),t_in_k={str(v):str(substitution[v]) for v in t},
            linear_part_determinant=str(s.Matrix(list(substitution.values())).jacobian(k).det()))

    # All possible constant-factor dimensions on R0: point, P1, origin P3.
    examples={
        'unique':QR.subs({**dict.fromkeys(k,0),k[5]:1,k[8]:1}),
        'exceptional_P1':QR.subs({**dict.fromkeys(k,0),k[5]:1}),
        'origin_P3':s.zeros(4)}
    factor_examples={name:dict(matrix=factor_matrix(Qv).tolist(),kernel=[list(v) for v in factor_matrix(Qv).nullspace()],
        projective_fiber_dimension=3-factor_matrix(Qv).rank()) for name,Qv in examples.items()}
    assert [v['projective_fiber_dimension'] for v in factor_examples.values()]==[0,1,3]

    # Distinguish coefficient rank from the true base flag incidence.
    from ambient_weight1_image import common_plane
    from weight1_structure import coefficient_map,form_matrices,infinitesimal_action
    p,b=common_plane();bs=b.subs(dict(zip(p,[0,3,5,7,11,13])));K=coefficient_map(bs)
    assert K.rank()==2
    diff=b.jacobian(p).row_join(s.Matrix.hstack(*(infinitesimal_action(bs,A) for A in sp_basis())))
    assert diff.rank()==9
    a=s.Symbol('a');v=s.Matrix([243*a,1,-23*a,-37*a]);mats=form_matrices(bs)
    contraction=s.Matrix.hstack(*[(v.T*M).T for M in mats]);ell=(v.T*J).T
    left=list(contraction);right=list(ell*ell.T)
    equations=[s.expand(left[i]*right[j]-left[j]*right[i]) for i in range(16) for j in range(i)]
    gb=s.groebner(equations,a,domain=s.QQ);assert list(gb)==[a]
    voff=s.Matrix([243,0,-23,-37]);coff=list(s.Matrix.hstack(*[(voff.T*M).T for M in mats]));loff=(voff.T*J).T;roff=list(loff*loff.T)
    assert any(coff[i]*roff[j]!=coff[j]*roff[i] for i in range(16) for j in range(i))
    diagnostic=dict(base_point=[0,3,5,7,11,13],K=K.tolist(),coefficient_rank=2,
        base_incidence_differential_rank=9,base_flag_chart='v=e_y+a*(243,0,-23,-37)',
        contraction_groebner_basis=list(map(str,gb)),off_chart_kernel_flag_excluded=True,
        warning='The base incidence has a contraction condition in addition to K_beta*v=0; it is not P(ker K_beta).')
    inputs=[Path(__file__),ROOT/'code/weight1_structure.py',ROOT/'code/ambient_weight1_image.py',
        ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json',
        ROOT/'research/notes/AMBIENT_FULL_FIRSTJET_FIBERS_20261103_v63.md',
        ROOT/'legacy/v07/quadratic_classification.tex',ROOT/'legacy/v07/zero_boundary_section.tex']
    result=dict(id='AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_v63',all_passed=True,field='QQ',
        coordinate_order=list(map(str,Y)),quadratic_monomials=list(map(str,MON)),cubic_monomials=list(map(str,cubics)),
        families=records,full_fixed_line_recognition=checks,legacy_crosswalk=crosswalk,factor_fiber_examples=factor_examples,
        rankdrop_diagnostic=diagnostic,
        scope='Full REDUCED CA/CB zero-first-jet fibers equal existing M0/M1 (dimensions13/13). Common-flag N zero fiber has reduced image R0 dimension12. This does not assert reducedness of scheme-theoretic target zero fibers or all intersection fibers.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seconds=time.monotonic()-started)
    def serial(v):
        if isinstance(v,s.Rational):return str(v)
        raise TypeError(type(v))
    (output/'report.json').write_text(json.dumps(result,indent=2,default=serial)+'\n')
    print('ZERO FIRSTJET CERTIFICATES SAVED',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
