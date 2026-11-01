#!/usr/bin/env python3
"""Independent full-open and exceptional-image audit of four mixed curve orbits."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def factors(expr,b):return {str(s.Poly(f,*b).monic().as_expr()) for f,e in s.factor_list(expr,*b)[1]}
def run():
    started=time.time();main=json.loads((ROOT/'reports/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.json').read_text());orbit=json.loads((ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json').read_text())['orbits'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)};eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];enc=[[(tuple(qi[i] for i in mm if i in qi),tuple(ui[i] for i in mm if i in ui),s.Rational(nu,de)) for mm,nu,de in eq] for eq in eqs];b=s.symbols('b0:9');field=s.QQ.frac_field(*b);zero=field.zero
    z,bu,be=primitive_beta(m);v=s.Matrix([be[i,j].diff(z[k]) for (i,j),k in COORDS]);base=v.jacobian(bu);rows=list(base.T.rref()[1]);left=base[rows,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(v,A)[rows,:] for A in sp_basis()]);out=[]
    for label,r in main['cases'].items():
        st=time.time();H=M(r['beta_kernel']);beta=H*s.Matrix(b);q=s.Matrix([s.Rational(orbit[label]['point61'][i]) for i in two]);fixed=next(i for i,c in enumerate(q) if c);cols=[i for i in range(45) if i!=fixed];qk=list(map(field.from_sympy,q));coeff=[];products={};rawJ=[[zero]*45 for _ in range(135)]
        for row,eq in enumerate(enc):
            ff={}
            for mm,bm,c in eq:
                if bm not in products:
                    p=field.one
                    for i in bm:p*=field.from_sympy(beta[i])
                    products[bm]=p
                c=field.from_sympy(c)*products[bm]
                if c:ff[mm]=ff.get(mm,zero)+c
            coeff.append(ff);point=zero
            for mm,c in ff.items():
                if len(mm)==2:
                    i,j=mm;point+=c*qk[i]*qk[j]
                    if i!=fixed:rawJ[row][cols.index(i)]+=c*qk[j]
                    if j!=fixed:rawJ[row][cols.index(j)]+=c*qk[i]
                elif len(mm)==1:rawJ[row][-1]+=c*qk[mm[0]]
            assert point==zero
        JD=DM(rawJ,(135,45),field);J=M(r['jacobian']);T=M(r['tangent']);L=M(r['left_kernel']);A=M(r['quadratic']);assert (JD-J.to_DM().convert_to(field)).is_zero_matrix;TD=T.to_DM().convert_to(field);LD=L.to_DM().convert_to(field);assert (JD*TD).is_zero_matrix and (LD*JD).is_zero_matrix
        freecols=[i for i in range(45) if i not in r['pivot_columns']];freerows=[i for i in range(135) if i not in r['pivot_rows']];assert T[freecols,:]==s.eye(2) and T[-1,:]==s.zeros(1,2) and L[:,freerows]==s.eye(92)
        det=JD.extract(r['pivot_rows'],r['pivot_columns']).det();assert det==field.from_sympy(s.sympify(r['jacobian_minor']));allowed=factors(s.sympify(r['jacobian_minor']),b)
        for matrix in [T,L,A]:
            for value in set(matrix):assert factors(s.denom(s.cancel(value)),b)<=allowed
        ET=[[zero]*2 for _ in range(46)];tt=TD.to_list()
        for j,i in enumerate(cols):ET[i]=tt[j]
        ET[45]=tt[-1];raw=[[zero]*3 for _ in range(135)]
        for row,ff in enumerate(coeff):
            for mm,c in ff.items():
                i,j=(*mm,*([45]*(2-len(mm))))
                if fixed in (i,j):continue
                raw[row][0]+=c*ET[i][0]*ET[j][0];raw[row][1]+=c*(ET[i][0]*ET[j][1]+ET[i][1]*ET[j][0]);raw[row][2]+=c*ET[i][1]*ET[j][1]
        assert (LD*DM(raw,(135,3),field)-A.to_DM().convert_to(field)).is_zero_matrix and A[:,:2]==s.zeros(92,2)
        nonzero=[i for i in range(92) if A[i,2]!=0];assert nonzero==r['quadratic_nonzero_rows'];coef=s.cancel(A[nonzero[0],2]);bad=allowed|factors(s.numer(coef),b);assert bad==set(r['bad_factors'])
        ar=acts.xreplace(dict(zip(bu,beta)));exceptions=[]
        expected=set()
        for item in main['exceptional_factors'][label]:
            f=s.sympify(item['factor']);coord=s.sympify(item['solve_variable']);sol=s.sympify(item['solve_value']);expected|=factors(f,b);assert s.degree(f,coord)==1 and s.cancel(f.subs(coord,sol))==0 and s.gcd(s.diff(f,coord),f.subs(coord,0))==1
            deg=s.Poly(f,*b).total_degree();assert deg in [1,2]
            if deg==2:assert s.hessian(f,b).rank()==4
            # Full16x18 source differential, independently of the saved8x10 normal matrix.
            bp=beta.subs(coord,sol);pars=[bb for bb in b if bb!=coord];D=bp.jacobian(pars).row_join(ar.subs(coord,sol));DD=D.to_DM().convert_to(field);rr,pv=DD.transpose().rref(method='GJ');rank=len(pv);assert rank==item['beta_image_dimension']<16
            # Exact left kernel verifies the upper bound, and a rational minor the reverse bound.
            rl=rr.to_list();fr=[i for i in range(16) if i not in pv];nr=[]
            for ff in fr:
                row=[zero]*16;row[ff]=field.one
                for j,pp in enumerate(pv):row[pp]=-rl[j][ff]
                nr.append(row)
            NN=DM(nr,(len(fr),16),field);assert (NN.extract(range(len(fr)),fr)-DM.eye(len(fr),field)).is_zero_matrix and (NN*DD).is_zero_matrix;sample=dict(zip(b,range(2,11)));sample_rank=D.subs(sample).rank()
            if sample_rank!=rank:
                sample=dict(zip(b,range(3,12)));sample_rank=D.subs(sample).rank()
            assert sample_rank==rank
            exceptions.append(dict(factor=str(f),degree=deg,hessian_rank=4 if deg==2 else None,geometrically_irreducible=True,full_source_beta_differential_rank=rank,exact_left_kernel_dimension=16-rank))
        assert expected==bad
        out.append(dict(label=label,raw135_jacobian_and_quadratic_verified=True,full_rank_on_minor_open=43,tangent_identity_rows=freecols,left_kernel_identity_columns=freerows,minor=str(s.sympify(r['jacobian_minor'])),square_coefficient=str(coef),all_certificate_denominators_supported_on_minor=True,all_bad_factors_exhausted=sorted(bad),exceptions=exceptions,seconds=time.time()-st));print(label,'PASS',time.time()-st,flush=True)
    assert sum(len(r['exceptions']) for r in out)==11
    paths=[Path(__file__),ROOT/'code/s61_r0_mixed_curve_exhaustion_20261101.py',ROOT/'reports/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.json',ROOT/'research/notes/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.md',ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json',ROOT/'reports/S61_R0_MIXED_ORBIT_REVIEW_20261101.json',ROOT/'reports/S61_R0_MIXED_CURVE_UNIFORM_RANK_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_R0_MIXED_CURVE_EXHAUSTION_REVIEW_20261101',all_passed=True,cases=out,finite_R_argument_accepted=True,no_smooth_actual_curve_assumption_required=True,all_four_orbits_removed_over_sufficiently_general_beta=True,completed_rings='k[[u,e]]/(e^2), t=0',scope='Allpoints of T1,T2,T4,T6 over sufficiently general beta; T5/T7 surfaces, binary/flag endpoints and full affine exhaustion remain separate.',seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
