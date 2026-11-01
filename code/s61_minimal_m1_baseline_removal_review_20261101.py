#!/usr/bin/env python3
"""Independent normal-minor, geometry, and finite-cover review for M1 baselines."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61
from itertools import combinations_with_replacement
from s61_model import X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from s61_minimal_junction_review_20261101 import standards
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
def normal_minor():
    st=time.time();m=S61();bs=s.symbols('b0:10');field=s.QQ.frac_field(*bs);R,*gens=ring(','.join(map(str,bs)),QQ);old=json.loads((ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json').read_text())['cases']['M1']['universal'];main=json.loads((ROOT/'reports/S61_MINIMAL_M1_BASELINE_REMOVAL_20261101.json').read_text());r=main['uniform'];M=lambda a:s.Matrix(a).applyfunc(s.sympify);A=M(old['leading_jacobian']);tc=M(old['t_column']);C=M(old['first_kernel']);pr=old['first_pivot_rows'];pc=old['first_pivot_columns'];U=s.zeros(44,len(pc));iv=A[pr,pc].inv()
    for j,i in enumerate(pc):U[i,:]=iv[j,:]
    E=M(old['residual_matrix']).row_join(M(old['residual_t_column']));rp=r['jacobian_minor_rows'];cp=r['jacobian_minor_columns'];chosen=r['quadratic_rows'];ED=E.to_DM().convert_to(field);Cr=ED.extract(chosen,cp)*ED.extract(rp,cp).inv();L=s.zeros(3,135)
    for i,j in enumerate(chosen):L[i,j]=1
    LD=L.to_DM().convert_to(field);proj=s.zeros(len(rp),135)
    for i,j in enumerate(rp):proj[i,j]=1
    LD=LD-Cr*proj.to_DM().convert_to(field);sel=s.zeros(len(pr),135)
    for i,j in enumerate(pr):sel[i,j]=1
    LD=LD-LD*A.to_DM().convert_to(field)*U.to_DM().convert_to(field)*sel.to_DM().convert_to(field)
    assert (LD*A.row_join(tc).to_DM().convert_to(field)).is_zero_matrix
    print('left coefficients',time.time()-st,flush=True)
    def poly(v):return R.from_expr(v.as_expr())
    rows=[];den=[]
    for row in LD.to_list():
     d=R.one
     for c in row:d=d.lcm(poly(c.denom))
     out=[]
     for c in row:
      q,rem=d.div(poly(c.denom));assert not rem;out.append(poly(c.numer)*q)
     rows.append(out);den.append(d)
    used=[j for j in range(135) if any(row[j] for row in rows)];print('rows used',len(used),time.time()-st,flush=True)
    T=M(r['normal_frame_polynomial']);cols=old['chart_columns'];fixed=next(i for i in range(45) if i not in cols);ET=[[R.zero,R.zero] for _ in range(45)]
    for j,i in enumerate(cols):ET[i]=[R.from_expr(v) for v in T[j,:]]
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};betap=[R.from_expr(v) for v in M(old['beta_kernel'])*s.Matrix(bs)];eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];prodc={};pairs={};HQ={}
    for row in used:
     ff={}
     for mon,nu,de in eqs[row]:
      mm=tuple(qi[i] for i in mon if i in qi)
      if len(mm)!=2 or fixed in mm:continue
      bm=tuple(one.index(i) for i in mon if i in one)
      if bm not in prodc:
       v=R.one
       for i in bm:v*=betap[i]
       prodc[bm]=v
      cc=prodc[bm]*QQ(nu,de)
      if cc:ff[mm]=ff.get(mm,R.zero)+cc
     val=[R.zero]*3
     for mm,c in ff.items():
      if mm not in pairs:
       i,j=mm;pairs[mm]=[ET[i][0]*ET[j][0],ET[i][0]*ET[j][1]+ET[i][1]*ET[j][0],ET[i][1]*ET[j][1]]
      for k in range(3):val[k]+=c*pairs[mm][k]
     HQ[row]=val
    print('quadratic polynomials',time.time()-st,flush=True)
    expected=M(r['quadratic_minor_matrix']);ds=list(map(s.sympify,r['normal_denominators']));norm=[(0,0),(0,1),(1,1)]
    for i in range(3):
     for j,(a,b) in enumerate(norm):
      lhs=sum((rows[i][row]*HQ[row][j] for row in used),R.zero);target=field.from_sympy(expected[i,j]*ds[a]*ds[b]);assert lhs*poly(target.denom)==den[i]*poly(target.numer),(i,j)
    print('all9 selected coefficients verified',time.time()-st,flush=True)

    assert ED.extract(rp,cp).det()==field.from_sympy(s.sympify(r['jacobian_minor']))
    # A three-term expansion after clearing row denominators avoids costly
    # multivariate gcds in fraction-field Gaussian elimination.
    pm=[];pd=[]
    for row in expected.to_DM().convert_to(field).to_list():
        dd=R.one
        for c in row:dd=dd.lcm(poly(c.denom))
        prow=[]
        for c in row:
            quot,rem=dd.div(poly(c.denom));assert not rem;prow.append(poly(c.numer)*quot)
        pm.append(prow);pd.append(dd)
    determinant=pm[0][0]*(pm[1][1]*pm[2][2]-pm[1][2]*pm[2][1])-pm[0][1]*(pm[1][0]*pm[2][2]-pm[1][2]*pm[2][0])+pm[0][2]*(pm[1][0]*pm[2][1]-pm[1][1]*pm[2][0])
    target=field.from_sympy(s.sympify(r['quadratic_minor']));assert determinant*poly(target.denom)==pd[0]*pd[1]*pd[2]*poly(target.numer)
    frame=M(r['full_kernel_frame']);assert frame[[24,36,42],:]==s.eye(3) and (A.row_join(tc).to_DM()*frame.to_DM()).is_zero_matrix
    assert old['first_rank']==25 and old['leading_rank']==41 and old['full_rank']==42
    factors=main['factor_polynomials'];assert set(factors)=={'main','P1','P2','P3','W2','W3','b7','b9'};allowed=s.prod(s.sympify(v) for v in factors.values());denfacts=set()
    for expression in [s.sympify(r['jacobian_minor']),s.sympify(r['quadratic_minor'])]:
        for part in s.fraction(expression):
            for fac,power in s.factor_list(part)[1]:assert s.rem(s.Poly(allowed,*bs),s.Poly(fac,*bs)).is_zero
    for v in set(frame):
     for fac,power in s.factor_list(s.denom(s.cancel(v)))[1]:
      assert s.rem(s.Poly(allowed,*bs),s.Poly(fac,*bs)).is_zero;denfacts.add(str(fac))
    return dict(all9_selected_quadratic_coefficients_verified=True,selected_raw_equation_rows=used,normal_minor=str(s.sympify(r['quadratic_minor'])),jacobian_minor=str(s.sympify(r['jacobian_minor'])),frame_identity_rows=[24,36,42],frame_denominator_factors=sorted(denfacts),seconds=time.time()-st)

def geometry_and_points():
    started=time.time();m=S61();main=json.loads((ROOT/'reports/S61_MINIMAL_M1_BASELINE_REMOVAL_20261101.json').read_text());old=json.loads((ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json').read_text())['cases']['M1']['universal'];one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};bs=s.symbols('b0:10');K=M(old['beta_kernel']);bv=K*s.Matrix(bs);factors={k:s.sympify(v) for k,v in main['factor_polynomials'].items()}
    eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];encoded=[[(tuple(qi[i] for i in mon if i in qi),tuple(one.index(i) for i in mon if i in one),s.Rational(nu,de)) for mon,nu,de in eq] for eq in eqs]
    def polys(beta,top=False):
        fs=[];products={}
        for eq in encoded:
            f={}
            for mm,bm,c in eq:
                if top and len(mm)!=2:continue
                if bm not in products:products[bm]=s.S.Zero if any(beta[i]==0 for i in bm) else s.prod(beta[i] for i in bm)
                if products[bm]:f.setdefault(mm,[]).append(c*products[bm])
            fs.append({mm:s.expand(s.Add(*terms)) for mm,terms in f.items()})
        return fs
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);a=s.Symbol('a');Q=s.zeros(4);Q[0,1]=Q[2,3]=x*x/2;Q[1,2]=-2*x*z;Q[1,3]=x*w;Q-=Q.T
    def coords(q):
        P=X[0]**2*J;P[1:,1:]+=q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed];Qcurve=Q.copy();Qcurve[1,3]+=a*x*x;Qcurve[3,1]-=a*x*x;qc=coords(Qcurve);assert qc==M(main['actual_curve_q'])
    derivative=qc.diff(a);derivative-=q0*derivative[fixed]/q0[fixed];assert M(main['uniform']['full_kernel_frame'])[:,0]==derivative[cols,:].col_join(s.zeros(1,1))
    FG=polys(bv,top=True);assert all(s.expand(sum(c*s.prod(qc[i] for i in mm) for mm,c in f.items()))==0 for f in FG)
    orbit=[]
    for H in sp_basis():
        dq=coords(sum(((H*Y)[i]*Q.diff(Y[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq)
    OD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*OD.nullspace());assert OD.rank()==7 and stab.cols==3
    zz,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);kr=list(K.T.rref()[1]);ki=K[kr,:].inv();ba=(acts*stab).xreplace(dict(zip(bu,bv)));bactions=(ki*ba[kr,:]).applyfunc(s.expand);assert (K*bactions-ba).applyfunc(s.expand)==s.zeros(16,3)
    preservation={}
    for label,f in factors.items():
        derivatives=[s.expand(sum(s.diff(f,b)*bactions[j,i] for j,b in enumerate(bs))) for i in range(3)];quot=[]
        for value in derivatives:
            qq,rem=s.div(s.Poly(value,*bs),s.Poly(f,*bs));assert rem.is_zero;quot.append(str(qq.as_expr()))
        preservation[label]=quot
    out=[]
    for rec in main['factor_cases']:
        label=rec['label'];f=factors[label];drop=s.sympify(rec['solved_variable']);sol=s.sympify(rec['solution']);assert s.degree(f,drop)==1 and s.cancel(f.subs(drop,sol))==0;assert s.gcd(s.diff(f,drop),f.subs(drop,0)) in [1,-1]
        pars=[b for b in bs if b!=drop];beta=bv.subs(drop,sol);BD=beta.jacobian(pars).row_join(acts.xreplace(dict(zip(bu,beta))))
        if label in ['b7','b9','W2']:
            N=M(rec['left_kernel']);assert (BD.to_DM()-M(rec['beta_differential']).to_DM()).is_zero_matrix and (BD.T.to_DM()*N.to_DM()).is_zero_matrix
            free=[i for i in range(16) if any(N[i,j]==1 and all(N[i,k]==0 for k in range(N.cols) if k!=j) for j in range(N.cols))];assert N[free,:].rank()==N.cols
            sample=dict(zip(bs,range(1,11)));rank=BD.subs(sample).rank();assert rank==16-N.cols==rec['generic_beta_rank'];out.append(dict(label=label,primitive_linear_irreducible=True,universal_beta_image_dimension=rank));continue
        sample={s.sympify(k):s.Rational(v) for k,v in rec['sample'].items()};assert f.xreplace(sample)==0;BD0=BD.subs(sample);assert BD0.rank()==16
        ld=rec['local_data'];point=list(map(s.Rational,ld['point61']));q=s.Matrix([point[i] for i in two]);bet=s.Matrix([point[i] for i in one]);assert q==q0 and bet==bv.xreplace(sample);F=polys(bet);A=s.zeros(135,45)
        for row,ff in enumerate(F):
            for mm,c in ff.items():
                if len(mm)==2:
                    for col,i in enumerate(cols):A[row,col]+=c*sum(q[mm[1-j]] for j in range(2) if mm[j]==i)
                elif len(mm)==1:A[row,-1]+=c*q[mm[0]]
        T=M(ld['tangent_basis']);assert A==M(ld['jacobian']) and A.rank()==42 and A[:,:44].rank()==41 and T.cols==3 and A*T==s.zeros(135,3)
        rows=list(A.T.to_DM().rref()[1]);piv=list(A.to_DM().rref()[1]);U=s.zeros(45,42);iv=A.extract(rows,piv).inv(method='DM')
        for j,i in enumerate(piv):U[i,:]=iv[j,:]
        E={1:s.zeros(46,3)}
        for j,i in enumerate(cols):E[1][i,:]=T[j,:]
        E[1][45,:]=T[-1,:];mons={d:list(combinations_with_replacement(range(3),d)) for d in [1,2,3]};res={}
        for degree in [2,3]:
            ix={mm:i for i,mm in enumerate(mons[degree])};raw=s.zeros(135,len(ix))
            for row,ff in enumerate(F):
                for mm,c in ff.items():
                    i,j=(*mm,*([45]*(2-len(mm))))
                    if fixed in [i,j]:continue
                    for aa in range(1,degree):
                        bb=degree-aa
                        for ii,vi in enumerate(E[aa][i,:]):
                            if not vi:continue
                            for jj,vj in enumerate(E[bb][j,:]):
                                if vj:raw[row,ix[tuple(sorted(mons[aa][ii]+mons[bb][jj]))]]+=c*vi*vj
            lift=-U*raw[rows,:];res[degree]=raw+A*lift;E[degree]=s.zeros(46,len(ix))
            for j,i in enumerate(cols):E[degree][i,:]=lift[j,:]
            E[degree][45,:]=lift[-1,:]
        allmons=mons[2]+mons[3];ix={mm:i for i,mm in enumerate(allmons)};rr,pv=res[2].row_join(res[3]).to_DM().rref();base=rr.to_Matrix()[:len(pv),:];multiples=[]
        for row in base.tolist():
            for vv in range(3):
                val=[s.S.Zero]*len(allmons)
                for j,mm in enumerate(mons[2]):
                    if row[j]:val[ix[tuple(sorted(mm+(vv,)))]]=row[j]
                if any(val):multiples.append(val)
        rr,pv=base.col_join(s.Matrix(multiples)).to_DM().rref();rr=rr.to_Matrix();v=s.symbols('v0:3');forms=[];ranks={}
        for row,pc in enumerate(pv):
            d=len(allmons[pc]);ranks[d]=ranks.get(d,0)+1;forms.append(sum(rr[row,j]*s.prod(v[i] for i in mm) for j,mm in enumerate(allmons) if len(mm)==d))
        gb0=s.groebner([g.subs(v[0],0) for g in forms],*reversed(v[1:]),order='grevlex');st0,hf0=standards(gb0);gb1=s.groebner([g.subs(v[0],1) for g in forms],*reversed(v[1:]),order='grevlex');st1,hf1=standards(gb1);assert len(st1)==3
        dq=qc.diff(a);dq-=q*dq[fixed]/q[fixed];tv=M(rec['actual_curve_tangent']);assert T*tv==dq[cols,:].col_join(s.zeros(1,1)) and tv[0]!=0
        assert all(s.expand(g.subs(dict(zip(v,tv))))==0 for g in forms)
        out.append(dict(label=label,primitive_linear_irreducible=True,beta_rank=16,raw135_full_rank=42,leading_rank=41,initial_ranks=ranks,zero_section_standard_degrees=hf0,localized_standard_degrees=hf1,generic_rank_upper=3,actual_same_beta_curve_verified=True));print(json.dumps(out[-1]),flush=True)
    mr=json.loads((ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json').read_text())['cases']['M1'];PP,RR,AA,BB,UU,VV,KK,LL=s.symbols('P R A B U V K L');sub={PP:bv[1],RR:bv[3],AA:bv[4],BB:bv[13],UU:bv[6],VV:bv[9],KK:bv[14]+5*bv[9],LL:5*bv[6]-bv[15]}
    assert s.cancel(s.sympify(mr['solution_a']).subs(sub,simultaneous=True)-factors['main']/(27*bs[7]))==0 and s.cancel(s.sympify(mr['solution_b']).subs(sub,simultaneous=True))==0 and s.cancel(s.sympify(mr['H']).subs(sub,simultaneous=True)-z*w)==0
    return dict(projective_orbit_dimension=7,stabilizer_dimension=3,factor_semi_invariant_quotients=preservation,cases=out,main_intersection_verified=True,all_eight_open_factors_covered=True,seconds=time.time()-started)

def run():
    start=time.time();normal=normal_minor();geo=geometry_and_points();paths=[Path(__file__),ROOT/'reports/S61_MINIMAL_M1_BASELINE_REMOVAL_20261101.json',ROOT/'code/s61_minimal_M1_baseline_removal_20261101.py',ROOT/'research/notes/S61_MINIMAL_M1_BASELINE_REMOVAL_20261101.md',ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json',ROOT/'reports/S61_MINIMAL_BASELINE_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'code/s61_minimal_junction_review_20261101.py']
    return dict(id='S61_MINIMAL_M1_BASELINE_REMOVAL_REVIEW_20261101',all_passed=True,normal_minor=normal,geometry_and_points=geo,accepted_all_finite_M1_baselines_removed_over_general_beta=True,scope='Both finite M1 free baselines; parameter-infinity flag endpoints and M0 baselines remain separate.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
