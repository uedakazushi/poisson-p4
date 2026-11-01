#!/usr/bin/env python3
"""Independent full-beta leading equations and finite-pullback slice audit."""
import argparse,hashlib,json
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
ROOT=Path(__file__).resolve().parents[1]

def verify():
    model=S61();path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json';report=json.loads(path.read_text())
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    u=s.symbols('u1:62');index={v:i for i,v in enumerate(two)};fixed=index[29];D=s.zeros(135,45);tcol=s.zeros(135,1);values=[]
    extra=json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']
    for row,eq in enumerate(model.reduced['equations']+[extra]):
        val=0
        for mon,num,den in eq:
            qs=[index[i] for i in mon if i in index];b=s.Rational(num,den)*s.prod(u[i] for i in mon if i in one)
            if len(qs)==2:
                i,j=qs
                if i==fixed and j==fixed:val+=b
                if i==fixed:D[row,j]+=b
                if j==fixed:D[row,i]+=b
            if len(qs)==1 and qs[0]==fixed:tcol[row,0]+=b
        values.append(s.expand(val))
    gb=s.groebner([v for v in values if v],*[u[i] for i in one]);assert list(gb)==[u[11],u[12]]
    sub={u[11]:0,u[12]:0};D=D.subs(sub);tcol=tcol.subs(sub)
    parse=lambda v:s.Matrix(v).applyfunc(s.sympify)
    K=parse(report['constant_kernel']);N=parse(report['residual_kernel'])
    assert D[:70,:]*K==s.zeros(70,29) and K.rank()==29
    residual=D[70:,:]*K;assert (residual*N).applyfunc(s.cancel)==s.zeros(65,11)
    assert N[report['residual_free_columns'],:]==s.eye(11)
    sample={s.Symbol(k):s.Rational(v) for k,v in report['sample'].items()}
    leading=D.subs(sample).extract(report['leading_minor_rows'],report['leading_minor_columns'])
    assert leading.det(method='domain-ge')==s.Rational(report['leading_minor_at_sample'])!=0
    full=D[:,[i for i in range(45) if i!=fixed]].row_join(tcol).subs(sample)
    assert full.extract(report['full_minor_rows'],report['full_minor_columns']).det(method='domain-ge')==s.Rational(report['full_minor_at_sample'])!=0
    z,bu,beta=primitive_beta(model);V=s.Matrix(z);eps=s.Symbol('eps');pairs=list(combinations(range(4),2))
    def vec(P):return s.Matrix([s.expand(P[i,j]).coeff(z[k]) for i,j in pairs for k in range(4)])
    basis=vec(beta).jacobian(bu);_,piv=basis.T.to_DM().rref();piv=list(piv);inverse=basis[piv,:].inv()
    Q=s.zeros(4);Q[1,2]=z[0]*z[3];Q[2,1]=-Q[1,2]
    omega=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]));actions=[];generators=[]
    for raw,stored in zip(report['stabilizer_generators'],report['primitive_actions']):
        A=parse(raw);generators.append(A);assert A*omega+omega*A.T==s.zeros(4)
        dq=A*Q+Q*A.T-sum(((A*V)[i]*Q.diff(z[i]) for i in range(4)),s.zeros(4))
        scalar=s.Poly(s.expand(dq[1,2]),*z).coeff_monomial(z[0]*z[3]);assert (dq-scalar*Q).applyfunc(s.expand)==s.zeros(4)
        M=s.eye(4)+eps*A;changed=M.T*beta.subs(dict(zip(V,M*V)),simultaneous=True)*M
        action=(inverse*vec(changed.diff(eps).subs(eps,0))[piv,:]).applyfunc(s.expand)
        assert action==parse(stored)
        assert action[bu.index(u[11])].subs(sub)==0 and action[bu.index(u[12])].subs(sub)==0
        actions.append(action)
    assert s.Matrix.hstack(*[A.reshape(16,1) for A in generators]).rank()==5
    fixed_beta={s.Symbol(k):s.Rational(v) for k,v in report['slice_fixed_coordinates'].items()};rows=[bu.index(v) for v in fixed_beta]
    determinant=s.Matrix.hstack(*actions)[rows,:].subs(sample|fixed_beta).det();assert determinant==553
    paths=[Path(__file__),path,ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_BINARY_UNIFORM_RANK_REVIEW_20261101',all_passed=True,
        all_beta_leading_ideal=['u12','u13'],leading_rank=34,full_chart_rank=35,
        universal_kernel_multiplications_verified=True,finite_pullback_primitive_actions_verified=5,
        slice_dominance_minor=str(determinant),proof_note='research/notes/S61_BINARY_UNIFORM_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,slice_minor=553)))
if __name__=='__main__':main()
