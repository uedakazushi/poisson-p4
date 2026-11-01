#!/usr/bin/env python3
"""Independent pointwise Jacobian and source-differential audit of four axes."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
ROOT=Path(__file__).resolve().parents[1]

def run():
    start=time.time();m=S61();main_path=ROOT/'reports/S61_MINIMAL_ALL_AXIS_REMOVAL_20261101.json'
    main=json.loads(main_path.read_text());assert main['all_passed']
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    encoded=m.reduced['equations']+[json.loads(extra_path.read_text())['universal_new_equation']]
    weights=m.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2];qi={i:j for j,i in enumerate(two)}
    z,u,beta=primitive_beta(m);pairs=list(combinations(range(4),2));coords=[(i,j,h) for i,j in pairs for h in range(4)]
    def vec(B):return s.Matrix([s.diff(B[i,j],z[h]) for i,j,h in coords])
    basis=vec(beta).jacobian(u);selector=list(basis.T.to_DM().rref()[1]);left=basis[selector,:].inv()
    J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]);eps=s.Symbol('eps');Z=s.Matrix(z);generators=[]
    for i in range(4):
        for j in range(i,4):
            C=s.zeros(4);C[i,j]=C[j,i]=1;H=J*C;assert H.T*J+J*H==s.zeros(4);generators.append(H)
    results=[]
    for case in main['cases']:
        pt=s.Matrix(list(map(s.Rational,case['sample_point61'])));q=pt[two,:];b=pt[one,:];A=s.zeros(135,45);tc=s.zeros(135,1);values=[]
        # Reconstruct derivatives directly from the immutable encoded61-variable
        # monomials; do not call the main script's specialization helper.
        for row,f in enumerate(encoded):
            value=0
            for mon,n,d in f:
                coefficient=s.Rational(n,d);qm=[]
                for i in mon:
                    if weights[i]==1:coefficient*=pt[i]
                    else:qm.append(qi[i])
                if len(qm)==2:
                    i,j=qm;value+=coefficient*q[i]*q[j];A[row,i]+=coefficient*q[j];A[row,j]+=coefficient*q[i]
                elif len(qm)==1:tc[row]+=coefficient*q[qm[0]]
            values.append(value)
        assert values==[0]*135 and A*q==s.zeros(135,1)
        assert len(A.to_DM().rref()[1])==44
        fixed=case['chart_fixed_q'];cols=[i for i in range(45) if i!=fixed];full=A[:,cols].row_join(tc)
        assert full==s.Matrix(case['homogenized_jacobian']).applyfunc(s.sympify)
        rows=case['original134_minor_rows'];M=full[rows,:];D=M.to_DM();det=D.domain.to_sympy(D.det())
        assert det==s.Rational(case['original134_minor_determinant']) and det!=0 and max(rows)<134
        # Actual finite pullback at this beta, differentiated in eps. This
        # independently avoids infinitesimal_action and its coordinate flattening.
        B=beta.subs(dict(zip(u,b)));action=[]
        for H in generators:
            G=s.eye(4)+eps*H;pull=G.T*B.subs(dict(zip(z,G*Z)),simultaneous=True)*G
            derivative=pull.diff(eps).subs(eps,0).applyfunc(s.expand);w=vec(derivative)
            coefficients=left*w[selector,:];assert basis*coefficients==w;action.append(coefficients)
        bv=s.Matrix(list(map(s.sympify,case['source_beta'])));params=list(map(s.Symbol,case['source_variables']))
        sample={s.Symbol('k'):1,s.Symbol('U'):2,s.Symbol('V'):3,s.Symbol('L'):5,s.Symbol('K'):19}
        sample.update(dict(zip(s.symbols('b0:6'),[2,3,5,7,11,13])))
        assert bv.subs(sample)==b
        differential=bv.jacobian(params).subs(sample).row_join(s.Matrix.hstack(*action))
        assert differential==s.Matrix(case['source_beta_differential']).applyfunc(s.sympify)
        pivot=list(differential.to_DM().rref()[1]);assert len(pivot)==16
        determinant=differential[:,pivot].det();assert determinant!=0
        expected=s.Rational(569) if case['axis']==6 else s.Rational(38082,25);assert determinant==expected
        results.append(dict(kind=case['kind'],axis=case['axis'],leading_rank=44,original134_full_minor=str(det),
            source_differential_columns=pivot,source_differential_minor=str(determinant),fixed_beta_local_ring='C',
            relative_local_ring_regular_dimension=16,whole_generic_axis_fiber_removed=True))
    paths=[Path(__file__),main_path,ROOT/'code/s61_minimal_all_axis_removal_20261101.py',
        ROOT/'research/notes/S61_MINIMAL_ALL_AXIS_REMOVAL_20261101.md',ROOT/'code/s61_model.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',extra_path]
    return dict(id='S61_MINIMAL_ALL_AXIS_REMOVAL_REVIEW_20261101',all_passed=True,cases=results,
        independently_checked=['All135 point values and full45 original134 minor, directly from61-variable monomials.',
            'Finite first-order pullback and actual sampled source differential at all four points.',
            'Relative regularity and proper16-dimensional image argument reviewed; no positive-dimensional-family extrapolation.'],
        dependence='Universal family equations and minimal reduced exhaustion remain dependencies of the reviewed main theorem.',
        proof_note='research/notes/S61_MINIMAL_ALL_AXIS_REMOVAL_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'source_minors':[c['source_differential_minor'] for c in r['cases']],'seconds':r['seconds']}))

if __name__=='__main__':main()
