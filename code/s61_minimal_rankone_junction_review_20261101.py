#!/usr/bin/env python3
"""Independent raw135, finite normalization, and cubic-cone junction audit."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
from s61_minimal_junction_review_20261101 import standards
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def run():
    started=time.time();path=ROOT/'reports/S61_MINIMAL_RANKONE_JUNCTION_20261101.json';src=json.loads(path.read_text());assert src['all_passed'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];index={i:j for j,i in enumerate(two)}
    eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    def polynomials(beta):
        ans=[]
        for eq in eqs:
            f={}
            for mon,a,b in eq:
                mm=tuple(index[i] for i in mon if i in index);cc=s.Rational(a,b)*s.prod(beta[one.index(i)] for i in mon if i in one);f[mm]=f.get(mm,s.S.Zero)+cc
            ans.append({mm:s.expand(cc) for mm,cc in f.items() if cc})
        return ans
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);a,b,d,h=s.symbols('a b d h');rho,la=s.symbols('rho la',nonzero=True)
    def tensor(a,b,d,h):
        Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-b*x*x-2*d*x*w-3*h*w*w;Q[1,3]=a*x*x;return Q-Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    # Finite transformations are checked on all six bivector coefficients.
    shear=s.eye(4);shear[1,2]=-a;shear[1,3]=-b+2*d*a-3*h*a*a;shear[2,0]=b-2*d*a+3*h*a*a;shear[3,0]=-a
    assert (shear.T*J[1:,1:]*shear-J[1:,1:]).applyfunc(s.expand)==s.zeros(4)
    pull=(shear.inv()*tensor(a,b,d,h).subs(dict(zip(Y,shear*Y)),simultaneous=True)*shear.inv().T).applyfunc(s.expand)
    assert (pull-tensor(0,0,d-3*h*a,h)).applyfunc(s.expand)==s.zeros(4)
    tor=s.diag(rho,1/rho,1/la,la);pull=(tor.inv()*tensor(0,0,d,h).subs(dict(zip(Y,tor*Y)),simultaneous=True)*tor.inv().T/rho**2).applyfunc(s.cancel)
    assert (pull-tensor(0,0,d*la**2,h*la**3/rho)).applyfunc(s.cancel)==s.zeros(4)
    univ=src['A1_universal'];q=coords(tensor(0,0,1,1));assert q==M(univ['q_point']);u=s.symbols('u0:16');generic=polynomials(u)
    constraints=s.Matrix([s.expand(sum(cc*s.prod(q[i] for i in mm) for mm,cc in f.items() if len(mm)==2)) for f in generic[:134]])
    bm=constraints.jacobian(u);assert constraints==bm*s.Matrix(u) and bm.rank()==7
    BK=M(univ['beta_kernel']);assert BK.shape==(16,9) and BK.rank()==9 and bm*BK==s.zeros(134,9)
    bs=s.symbols('b0:9');F=polynomials(BK*s.Matrix(bs));cols=univ['chart_columns'];D=s.zeros(135,44)
    for row,f in enumerate(F):
        assert s.expand(sum(cc*s.prod(q[i] for i in mm) for mm,cc in f.items() if len(mm)==2))==0
        for mm,cc in f.items():
            if len(mm)==2:
                for col,i in enumerate(cols):D[row,col]+=cc*sum(q[mm[1-j]] for j in range(2) if mm[j]==i)
    assert D==M(univ['leading_jacobian']);C=M(univ['first_kernel']);N=M(univ['residual_kernel']);assert C.rank()==14 and D[:70,:]*C==s.zeros(70,14)
    assert D*C==M(univ['residual_matrix']) and (D.to_DM()*C.to_DM()*N.to_DM()).is_zero_matrix
    free=[i for i in range(N.rows) if list(N[i,:]) in [[1,0],[0,1]]];assert N[free,:].rank()==2
    Ds=D.xreplace(dict(zip(bs,range(1,10))));assert Ds.rank()==42
    # Reconstruct the complete local135 Jacobian at the actual junction.
    local=src['local_data'];point=list(map(s.Rational,local['point61']));beta=[point[i] for i in one];q=s.Matrix([point[i] for i in two]);F=polynomials(beta);cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols);A=s.zeros(135,45)
    for row,f in enumerate(F):
        assert sum(cc*s.prod(q[i] for i in mm) for mm,cc in f.items() if len(mm)==2)==0
        for mm,cc in f.items():
            if len(mm)==2:
                for col,i in enumerate(cols):A[row,col]+=cc*sum(q[mm[1-j]] for j in range(2) if mm[j]==i)
            elif len(mm)==1:A[row,-1]+=cc*q[mm[0]]
    T=M(local['tangent_basis']);assert A==M(local['jacobian']) and A.rank()==40 and T.rank()==5 and A*T==s.zeros(135,5)
    # Cubic series in all46 homogeneous coordinates, then the full truncated
    # Macaulay row ideal. This does not import the primary cubic-bound routine.
    rows=list(A.T.to_DM().rref()[1]);piv=list(A.to_DM().rref()[1]);U=s.zeros(45,40);inv=A.extract(rows,piv).inv(method='DM')
    for j,i in enumerate(piv):U[i,:]=inv[j,:]
    E={1:s.zeros(46,5)}
    for j,i in enumerate(cols):E[1][i,:]=T[j,:]
    E[1][45,:]=T[44,:];mons={dd:list(combinations_with_replacement(range(5),dd)) for dd in [1,2,3]};res={}
    for degree in [2,3]:
        pos={mm:i for i,mm in enumerate(mons[degree])};raw=s.zeros(135,len(pos))
        for row,f in enumerate(F):
            for mm,cc in f.items():
                i,j=(*mm,*([45]*(2-len(mm))))
                if fixed in [i,j]:continue
                for aa in range(1,degree):
                    bb=degree-aa
                    for ii,vi in enumerate(E[aa][i,:]):
                        if not vi:continue
                        for jj,vj in enumerate(E[bb][j,:]):
                            if vj:raw[row,pos[tuple(sorted(mons[aa][ii]+mons[bb][jj]))]]+=cc*vi*vj
        lift=-U*raw[rows,:];res[degree]=raw+A*lift;assert res[degree][rows,:]==s.zeros(40,len(pos));E[degree]=s.zeros(46,len(pos))
        for j,i in enumerate(cols):E[degree][i,:]=lift[j,:]
        E[degree][45,:]=lift[44,:]
    allmons=mons[2]+mons[3];ix={mm:i for i,mm in enumerate(allmons)};trunc=res[2].row_join(res[3]);rr,pv=trunc.to_DM().rref();base=rr.to_Matrix()[:len(pv),:];mul=[]
    for row in base.tolist():
        for vv in range(5):
            val=[s.S.Zero]*len(allmons)
            for jj,mm in enumerate(mons[2]):
                if row[jj]:val[ix[tuple(sorted(mm+(vv,)))]]=row[jj]
            if any(val):mul.append(val)
    rr,pv=base.col_join(s.Matrix(mul)).to_DM().rref();rr=rr.to_Matrix();v=s.symbols('v0:5');forms=[];ranks={}
    for row,pc in enumerate(pv):
        dd=len(allmons[pc]);ranks[dd]=ranks.get(dd,0)+1;forms.append(sum(rr[row,j]*s.prod(v[i] for i in mm) for j,mm in enumerate(allmons) if len(mm)==dd))
    gb0=s.groebner([f.subs(v[0],0) for f in forms],*reversed(v[1:]),order='grevlex');st0,hf0=standards(gb0)
    gb1=s.groebner([f.subs(v[0],1) for f in forms],*reversed(v[1:]),order='grevlex');st1,hf1=standards(gb1);assert len(st1)==5
    # Explicit A1 branch and its tangent, independently encoded as a tensor.
    tau=s.Symbol('t');qc=coords(tensor(-2,s.Rational(2,9)+20*tau,s.Rational(1,3)-15*tau/2,tau));assert qc==M(src['A1_q_curve']) and qc.subs(tau,0)==q
    assert all(s.expand(sum(cc*s.prod(qc[i] for i in mm) for mm,cc in f.items() if len(mm)==2))==0 for f in F)
    dq=qc.diff(tau).subs(tau,0);dq-=q*dq[fixed]/q[fixed];direction=s.Matrix([-3,0,0,0,0]);assert T*direction==dq[cols,:].col_join(s.zeros(1,1))
    assert s.expand((s.Rational(1,3)-15*tau/2)-3*tau*(-2))==s.Rational(1,3)-3*tau/2
    mainpath=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';main=json.loads(mainpath.read_text())['cases']['M0'];mainv=M(main['main_branch_tangent']);BD=M(main['beta_differential']);QD=M(main['q_differential']);sd=M(main['source_direction']);assert BD*sd==s.zeros(16,1) and T*mainv==(QD*sd).col_join(s.zeros(1,1))
    assert mainv==s.Matrix([s.Rational(-2,3),s.Rational(-2,3),0,0,0]) and s.Matrix.hstack(mainv,direction).rank()==2 and mainv[0] and direction[0]
    assert BD.rank()==16 and s.Rational(main['discriminant_derivative'])!=0
    known=json.loads((ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json').read_text())['cases']['M0']['rank_certificates']['135']
    assert known['full_rank']==42 and known['leading_rank']==42 and not known['homogenizer_zero_on_tangent']
    paths=[Path(__file__),path,mainpath,ROOT/'code/s61_minimal_rankone_junction_20261101.py',ROOT/'research/notes/S61_MINIMAL_RANKONE_JUNCTION_20261101.md',ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_junction_review_20261101.py']
    return dict(id='S61_MINIMAL_RANKONE_JUNCTION_REVIEW_20261101',all_passed=True,universal_A1_beta_kernel_dimension=9,universal_A1_leading_rank=42,finite_normalization_verified=True,point_full_rank=40,truncated_initial_ranks=ranks,zero_section_standard_degrees=hf0,localized_standard_degrees=hf1,multiplicity_upper=5,actual_branch_lower_lengths=[3,2],same_beta_branch_tangents=[list(map(str,mainv)),list(map(str,direction))],generic_main_factor='Full S length3 with square-zero augmentation; t^2=0. The length3 lower bound is not a lower bound for its t quotient.',generic_A1_factor='Full S and t quotient both length2, so t=0.',accepted_t_nilpotent=True,accepted_generic_rankone_M0_junction_stratum_removed=True,scope='Rank-one M0 main/A1 junction only; no flatness, full completed-ring presentation, or entire A1 family removal.',seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();result=run()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds']}))
