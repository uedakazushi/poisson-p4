#!/usr/bin/env python3
"""Independent exact-field audit of the nonbinary double-root minor open."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement,permutations
import sympy as s
from s61_model import S61,X,J
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def run():
    start=time.time();path=ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json';report=json.loads(path.read_text());assert report['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};u=s.symbols('u1:62');beta={i:u[i] for i in one};beta.update({3:0,4:0,11:0,12:0,49:-5*u[34],50:5*u[25]})
    source=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];F=[]
    for eq in source:
        f={}
        for mon,a,b in eq:
            mm=tuple(qi[i] for i in mon if i in qi);cc=s.Rational(a,b)*s.prod(beta[i] for i in mon if i in beta);f[mm]=f.get(mm,s.S.Zero)+cc
        F.append({mm:s.expand(cc) for mm,cc in f.items() if cc})
    fixed=qi[45];cols=[i for i in range(45) if i!=fixed];q=s.zeros(45,1);q[fixed]=-2;Jac=s.zeros(135,45)
    for row,f in enumerate(F):
        assert sum(cc*s.prod(q[i] for i in mm) for mm,cc in f.items() if len(mm)==2)==0
        for mm,cc in f.items():
            if len(mm)==2:
                for j,i in enumerate(cols):Jac[row,j]+=cc*sum(q[mm[1-h]] for h in range(2) if mm[h]==i)
            elif len(mm)==1:Jac[row,-1]+=cc*q[mm[0]]
    assert Jac==M(report['jacobian']);T=M(report['tangent']);L=M(report['left_kernel']);assert (Jac.to_DM()*T.to_DM()).is_zero_matrix and (L.to_DM()*Jac.to_DM()).is_zero_matrix
    free=[i for i in range(45) if i not in report['pivot_columns']];leftfree=[i for i in range(135) if i not in report['pivot_rows']];assert T[free,:]==s.eye(3) and L[:,leftfree]==s.eye(93) and T[-1,:]==s.zeros(1,3)
    dm=Jac.extract(report['pivot_rows'],report['pivot_columns']).to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));assert s.expand(det-s.sympify(report['jacobian_minor']))==0
    # Three independent q-tangent vectors give q-rank<=41; full-rank42 gives >=41.
    assert all(i<44 for i in free)
    h=s.sympify(report['normal_form_linear_shift']);denoms=set()
    for value in set(T).union(set(L)).union({h}):
        for fac,power in s.factor_list(s.denom(s.cancel(value)))[1]:
            assert s.rem(s.Poly(det,*[u[i] for i in one]),s.Poly(fac,*[u[i] for i in one])).is_zero;denoms.add(str(fac))
    field=s.QQ.frac_field(*[u[i] for i in one]);zero=field.zero;cache={}
    def cv(value):
        if value not in cache:cache[value]=field.from_sympy(value)
        return cache[value]
    ET=[[zero]*3 for _ in range(46)]
    for j,i in enumerate(cols):ET[i]=[cv(T[j,c]) for c in range(3)]
    ET[45]=[cv(T[-1,c]) for c in range(3)];mons=list(combinations_with_replacement(range(3),2));H=[[zero]*6 for _ in range(135)]
    for row,f in enumerate(F):
        for mm,c in f.items():
            if fixed in mm:continue
            i,j=(*mm,*([45]*(2-len(mm))))
            if not any(ET[i]) or not any(ET[j]):continue
            cc=cv(c)
            for col,(a,b) in enumerate(mons):H[row][col]+=cc*(ET[i][a]*ET[j][b]+(ET[i][b]*ET[j][a] if a!=b else zero))
    obs=[]
    for row in L.tolist():
        val=[zero]*6
        for j,c in enumerate(row):
            if c:
                cc=cv(c)
                for col in range(6):val[col]+=cc*H[j][col]
        obs.append(val)
    assert obs==[[cv(s.sympify(c)) for c in row] for row in report['quadratic']]
    assert [i for i,row in enumerate(obs) if any(row)]==[41,46,57,73]
    hh=cv(h);normal=[[cv(1),zero,2*hh,zero,zero,hh*hh],[zero,cv(1),zero,zero,hh,zero],[zero,zero,zero,cv(1),zero,zero]]
    for row in obs:
        assert row==[sum((row[j]*normal[k][col] for k,j in enumerate([0,1,3])),zero) for col in range(6)]
    C=[[obs[i][j] for j in [0,1,3]] for i in [41,46,57]];nd=zero
    for p in permutations(range(3)):
        inversions=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3));term=cv((-1)**inversions)
        for i in range(3):term*=C[i][p[i]]
        nd+=term
    expected=u[25]*u[34]*(u[13]-2*u[48])**2*(5*u[13]-4*u[48])**2/(4*u[48]*(2*u[13]-u[48])**3)
    assert nd==cv(expected)==cv(s.sympify(report['normal_quadratic_coefficient_minor']))
    for fac,power in s.factor_list(s.fraction(s.cancel(expected))[0])[1]:assert s.rem(s.Poly(det,*[u[i] for i in one]),s.Poly(fac,*[u[i] for i in one])).is_zero
    # Finite symmetry and identification of the actual q=-2 coordinate convention.
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=z*z;Q[1,3]=-2*z*w;Q-=Q.T;Biv=X[0]**2*J;Biv[1:,1:]+=Q;pp=m.normalize(Biv)['parameters'];assert s.Matrix([pp[i] for i in two])==q
    S=M(report['symmetry_matrix']);assert S.T*J[1:,1:]*S==J[1:,1:];target=s.zeros(4);target[1,2]=2*z*w;target[1,3]=-w*w;target-=target.T
    assert (S.inv()*Q.subs(dict(zip(Y,S*Y)),simultaneous=True)*S.inv().T-target).applyfunc(s.expand)==s.zeros(4)
    paths=[Path(__file__),path,ROOT/'code/s61_nonbinary_double_root_open_20261101.py',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.md',m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'research/notes/S61_R0_UNIFORM_VERTEX_20261101.md',ROOT/'research/notes/S61_R0_GENERIC_EXHAUSTION_20261101.md',ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_NONBINARY_DOUBLE_ROOT_REVIEW_20261101',all_passed=True,field='QQ(beta)',full_rank=42,leading_rank=41,universal_minor=str(det),denominator_factors=sorted(denoms),normal_shift=str(h),normal_quadratic_minor=str(expected),finite_symmetry_verified=True,accepted_completed_ring='C[[u,e1,e2]]/(e1,e2)^2 with t=0',scope='The explicit minor open over sufficiently general beta only; no exceptional-divisor or whole-orbit removal.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
