#!/usr/bin/env python3
"""Independent normal-minor, geometry, and finite-cover review for M0 baselines."""
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
    st=time.time();m=S61();bs=s.symbols('b0:10');field=s.QQ.frac_field(*bs);R,*gens=ring(','.join(map(str,bs)),QQ);old=json.loads((ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json').read_text())['cases']['M0']['universal'];main=json.loads((ROOT/'reports/S61_MINIMAL_M0_BASELINE_REMOVAL_20261101.json').read_text());r=main['uniform'];M=lambda a:s.Matrix(a).applyfunc(s.sympify);A=M(old['leading_jacobian']);tc=M(old['t_column']);C=M(old['first_kernel']);pr=old['first_pivot_rows'];pc=old['first_pivot_columns'];U=s.zeros(44,len(pc));iv=A[pr,pc].inv()
    for j,i in enumerate(pc):U[i,:]=iv[j,:]
    E=M(old['residual_matrix']).row_join(M(old['residual_t_column']));rp=r['jacobian_minor_rows'];cp=r['jacobian_minor_columns'];chosen=r['quadratic_rows'];ED=E.to_DM().convert_to(field);Cr=ED.extract(chosen,cp)*ED.extract(rp,cp).inv();L=s.zeros(6,135)
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
    T=M(r['normal_frame_polynomial']);cols=old['chart_columns'];fixed=next(i for i in range(45) if i not in cols);ET=[[R.zero]*3 for _ in range(46)]
    for j,i in enumerate(cols):ET[i]=[R.from_expr(v) for v in T[j,:]]
    ET[45]=[R.from_expr(v) for v in T[-1,:]]
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};betap=[R.from_expr(v) for v in M(old['beta_kernel'])*s.Matrix(bs)];eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']];prodc={};pairs={};HQ={}
    for row in used:
     ff={}
     for mon,nu,de in eqs[row]:
      mm=tuple(qi[i] for i in mon if i in qi)
      if fixed in mm:continue
      mm=mm+(45,)*(2-len(mm))
      bm=tuple(one.index(i) for i in mon if i in one)
      if bm not in prodc:
       v=R.one
       for i in bm:v*=betap[i]
       prodc[bm]=v
      cc=prodc[bm]*QQ(nu,de)
      if cc:ff[mm]=ff.get(mm,R.zero)+cc
     val=[R.zero]*6
     for mm,c in ff.items():
      if mm not in pairs:
       i,j=mm;pairs[mm]=[ET[i][a]*ET[j][b]+(ET[i][b]*ET[j][a] if a!=b else R.zero) for a,b in combinations_with_replacement(range(3),2)]
      for k in range(6):val[k]+=c*pairs[mm][k]
     HQ[row]=val
    print('quadratic polynomials',time.time()-st,flush=True)
    expected=M(r['quadratic_minor_matrix']);ds=list(map(s.sympify,r['normal_denominators']));norm=list(combinations_with_replacement(range(3),2))
    for i in range(6):
     for j,(a,b) in enumerate(norm):
      lhs=sum((rows[i][row]*HQ[row][j] for row in used),R.zero);target=field.from_sympy(expected[i,j]*ds[a]*ds[b]);assert lhs*poly(target.denom)==den[i]*poly(target.numer),(i,j)
    print('all36 selected coefficients verified',time.time()-st,flush=True)

    assert ED.extract(rp,cp).det()==field.from_sympy(s.sympify(r['jacobian_minor']))
    # Clear row denominators before a polynomial determinant.
    pm=[];pd=[]
    for row in expected.to_DM().convert_to(field).to_list():
        dd=R.one
        for c in row:dd=dd.lcm(poly(c.denom))
        prow=[]
        for c in row:
            quot,rem=dd.div(poly(c.denom));assert not rem;prow.append(poly(c.numer)*quot)
        pm.append(prow);pd.append(dd)
    from sympy.polys.matrices import DomainMatrix
    ddm=s.Matrix([[v.as_expr() for v in row] for row in pm]).to_DM();determinant=R.from_expr(ddm.domain.to_sympy(ddm.det()))
    totalden=R.one
    for d in pd:totalden*=d
    target=field.from_sympy(s.sympify(r['quadratic_minor']));assert determinant*poly(target.denom)==totalden*poly(target.numer)
    frame=M(r['full_kernel_frame']);assert frame[[18,37,42,44],:]==s.eye(4) and (A.row_join(tc).to_DM()*frame.to_DM()).is_zero_matrix
    assert old['first_rank']==25 and old['leading_rank']==41 and old['full_rank']==41
    factors=[bs[1],bs[7],bs[8],8*bs[1]-bs[9],11*bs[1]-bs[9],194*bs[1]-3*bs[9],44*bs[1]-3*bs[9],3*bs[1]-bs[9],4938*bs[1]**2-727*bs[1]*bs[9]+27*bs[9]**2];allowed=s.prod(factors);denfacts=set()
    for expression in [s.sympify(r['jacobian_minor']),s.sympify(r['quadratic_minor'])]:
        for part in s.fraction(expression):
            for fac,power in s.factor_list(part)[1]:assert s.rem(s.Poly(allowed,*bs),s.Poly(fac,*bs)).is_zero
    for v in set(frame):
        for fac,power in s.factor_list(s.denom(s.cancel(v)))[1]:
            assert s.rem(s.Poly(allowed,*bs),s.Poly(fac,*bs)).is_zero;denfacts.add(str(fac))
    return dict(all36_selected_quadratic_coefficients_verified=True,selected_raw_equation_rows=used,normal_minor=str(s.sympify(r['quadratic_minor'])),jacobian_minor=str(s.sympify(r['jacobian_minor'])),frame_identity_rows=[18,37,42,44],t_tangent_row_is_constant=[0,0,0,1],frame_denominator_factors=sorted(denfacts),all_open_factors_list=list(map(str,factors)),seconds=time.time()-st)

def review_geometry_and_points():
    started=time.time();m=S61();main=json.loads((ROOT/'reports/S61_MINIMAL_M0_BASELINE_REMOVAL_20261101.json').read_text());old=json.loads((ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json').read_text())['cases']['M0']['universal'];one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};bs=s.symbols('b0:10');K=M(old['beta_kernel']);bv=K*s.Matrix(bs)
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
    x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);gamma=s.Symbol('gamma');Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-2*x*z;Q[1,3]=2*x*w;Q-=Q.T
    def coords(q):
        P=X[0]**2*J;P[1:,1:]+=q;return s.Matrix(m.normalize(P)['parameters'])[two,:]
    q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed];assert cols==old['chart_columns'];Qcurve=Q.copy();Qcurve[1,2]-=gamma*x*x;Qcurve[2,1]+=gamma*x*x;qc=coords(Qcurve);assert qc==s.Matrix([s.sympify(v,locals={'gamma':gamma}) for v in main['actual_curve_q']])
    dq=qc.diff(gamma);dq-=q0*dq[fixed]/q0[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));assert M(main['uniform']['full_kernel_frame'])[:,0]==-tv
    FG=polys(bv,top=True);assert all(s.expand(sum(c*s.prod(qc[i] for i in mm) for mm,c in f.items()))==0 for f in FG)
    orbit=[]
    for H in sp_basis():
        d=coords(sum(((H*Y)[i]*Q.diff(Y[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T);d-=q0*d[fixed]/q0[fixed];orbit.append(d)
    OD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*OD.nullspace());assert OD.rank()==7 and stab.cols==3
    zz,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);kr=list(K.T.rref()[1]);ki=K[kr,:].inv();ba=(acts*stab).xreplace(dict(zip(bu,bv)));bactions=(ki*ba[kr,:]).applyfunc(s.expand);assert (K*bactions-ba).applyfunc(s.expand)==s.zeros(16,3)
    alpha=s.Symbol('a');factors={'b1':bs[1],'b7':bs[7],'b8':bs[8],'all_ratio_hyperplanes':alpha*bs[1]-bs[9]};preservation={}
    for label,f in factors.items():
        derivatives=[s.expand(sum(s.diff(f,b)*bactions[j,i] for j,b in enumerate(bs))) for i in range(3)];quot=[]
        for value in derivatives:
            qq,rem=s.div(s.Poly(value,*bs),s.Poly(f,*bs));assert rem.is_zero;quot.append(str(qq.as_expr()))
        preservation[label]=quot
    proper=[]
    for rec in main['geometry']['proper_beta_images']:
        drop=s.sympify(rec['factor']);pars=[b for b in bs if b!=drop];beta=bv.subs(drop,0);BD=beta.jacobian(pars).row_join(acts.xreplace(dict(zip(bu,beta))));N=M(rec['left_kernel']);assert (BD.to_DM()-M(rec['beta_differential']).to_DM()).is_zero_matrix and (BD.T.to_DM()*N.to_DM()).is_zero_matrix;assert N.rank()==N.cols
        rank=BD.subs(dict(zip(bs,range(1,11)))).rank();assert rank==16-N.cols==rec['image_dimension'];proper.append(dict(factor=str(drop),universal_beta_image_dimension=rank))
    def raw_jacobian(F):
        A=s.zeros(135,45)
        for row,ff in enumerate(F):
            for mm,c in ff.items():
                if len(mm)==2:
                    for col,i in enumerate(cols):A[row,col]+=c*sum(q0[mm[1-j]] for j in range(2) if mm[j]==i)
                elif len(mm)==1:A[row,-1]+=c*q0[mm[0]]
        return A.applyfunc(s.expand)
    rational=[]
    for rec in main['rational_exceptional_strata']:
        ratio=s.Rational(rec['ratio']);sample={s.sympify(k):s.Rational(v) for k,v in rec['sample'].items()};assert sample[bs[9]]==ratio*sample[bs[1]];beta=bv.xreplace(sample);bsource=bv.subs(bs[9],ratio*bs[1]);BD=bsource.jacobian(bs[:9]).row_join(acts.xreplace(dict(zip(bu,beta))));assert BD.rank()==16
        ld=rec['local_data'];point=list(map(s.Rational,ld['point61']));assert s.Matrix([point[i] for i in two])==q0 and s.Matrix([point[i] for i in one])==beta;F=polys(beta);A=raw_jacobian(F);T=M(ld['tangent_basis']);assert A==M(ld['jacobian']) and A.rank()==41 and A[:,:44].rank()==ld['leading_jacobian_rank'] and T.cols==4 and A*T==s.zeros(135,4)
        vt=M(rec['actual_curve_tangent']);assert T*vt==tv and vt[0]!=0
        depth=4 if ratio==3 else 2;initial=initial_bound(F,A,T,cols,fixed,vt,depth);assert initial['localized_length']==4
        if ratio!=3:assert initial['initial_ranks'][2]==6 and initial['zero_section_length']==4
        rational.append(dict(ratio=str(ratio),raw135_full_rank=41,leading_rank=A[:,:44].rank(),restricted_beta_rank=16,initial=initial,actual_same_beta_curve_verified=True));print('rational',ratio,'PASS',flush=True)
    field=field_point(main['algebraic_exceptional_strata'],bv,bs,bu,acts,polys,raw_jacobian,tv,cols,fixed)
    mr=json.loads((ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json').read_text())['cases']['M0'];PP,RR,AA,BB,UU,VV,KK,LL=s.symbols('P R A B U V K L');sub={PP:bv[1],RR:bv[3],AA:bv[4],BB:bv[13],UU:bv[6],VV:bv[9],KK:bv[14]+5*bv[9],LL:5*bv[6]-bv[15]}
    assert s.cancel(s.sympify(mr['solution_a']).subs(sub,simultaneous=True))==0 and s.cancel(s.sympify(mr['solution_b']).subs(sub,simultaneous=True)-(8*bs[1]-bs[9])/(6*bs[7]))==0
    return dict(projective_orbit_dimension=7,stabilizer_dimension=3,factor_semi_invariant_quotients=preservation,proper_beta_images=proper,rational_strata=rational,algebraic_strata=field,main_ratio8_verified=True,all_open_factors_covered=True,seconds=time.time()-started)

def initial_bound(F,A,T,cols,fixed,vt,depth):
    """Independent graded46-coordinate recursion and full truncated multiples."""
    rows=list(A.T.to_DM().rref()[1]);piv=list(A.to_DM().rref()[1]);U=s.zeros(45,len(piv));iv=A.extract(rows,piv).inv(method='DM')
    for j,i in enumerate(piv):U[i,:]=iv[j,:]
    n=T.cols;E={1:s.zeros(46,n)}
    for j,i in enumerate(cols):E[1][i,:]=T[j,:]
    E[1][45,:]=T[-1,:];mons={d:list(combinations_with_replacement(range(n),d)) for d in range(1,depth+1)};res={}
    for degree in range(2,depth+1):
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
    allmons=sum([mons[d] for d in range(2,depth+1)],[]);ix={mm:i for i,mm in enumerate(allmons)};rr,pv=s.Matrix.hstack(*[res[d] for d in range(2,depth+1)]).to_DM().rref();base=rr.to_Matrix()[:len(pv),:];multiples=[]
    for row in base.tolist():
        for factor in sum([mons[d] for d in range(1,depth-1)],[]):
            val=[s.S.Zero]*len(allmons)
            for j,mm in enumerate(allmons):
                if len(mm)+len(factor)<=depth and row[j]:val[ix[tuple(sorted(mm+factor))]]+=row[j]
            if any(val):multiples.append(val)
    mac=base.col_join(s.Matrix(multiples)) if multiples else base;rr,pv=mac.to_DM().rref();rr=rr.to_Matrix();v=s.symbols(f'v0:{n}');forms=[];ranks={}
    for row,pc in enumerate(pv):
        d=len(allmons[pc]);ranks[d]=ranks.get(d,0)+1;forms.append(sum(rr[row,j]*s.prod(v[i] for i in mm) for j,mm in enumerate(allmons) if len(mm)==d))
    assert all(s.expand(g.xreplace(dict(zip(v,vt))))==0 for g in forms)
    gb0=s.groebner([g.subs(v[0],0) for g in forms],*reversed(v[1:]),order='grevlex');st0,hf0=standards(gb0);gb1=s.groebner([g.subs(v[0],1) for g in forms],*reversed(v[1:]),order='grevlex');st1,hf1=standards(gb1)
    return dict(independent_graded_series_depth=depth,full_truncated_macaulay_shape=list(mac.shape),initial_ranks=ranks,zero_section_standard_degrees=hf0,zero_section_length=len(st0),localized_standard_degrees=hf1,localized_length=len(st1),reversed_variable_order=True)

def field_point(rec,bv,bs,bu,acts,polys,raw_jacobian,tv,cols,fixed):
    from sympy.polys.matrices import DomainMatrix as DM
    a=s.Symbol('a');f=27*a*a-727*a+4938;assert s.discriminant(f,a)==-25*191 and s.Poly(f,a).is_irreducible;root=s.CRootOf(f,1);field=QQ.algebraic_field(root);alpha=field.unit
    def cv(e):return sum((field.from_sympy(c)*alpha**mm[0] for mm,c in s.Poly(e,a).terms()),field.zero)
    def dm(MM):return DM([[cv(e) for e in row] for row in MM.tolist()],MM.shape,field)
    sample=dict(zip(bs,range(1,11)));sample[bs[9]]=2*a;beta=bv.xreplace(sample);assert list(map(str,beta))==rec['beta_sample'];F=polys(beta);A=dm(raw_jacobian(F));T=dm(M(rec['tangent_basis']));assert A.rank()==41 and A.extract(range(135),range(44)).rank()==41 and (A*T).is_zero_matrix and T.rank()==4
    # Derive our own left nullspace, instead of importing stored cokernel rows.
    L=A.transpose().nullspace();assert (L*A).is_zero_matrix
    tt=T.to_list();ET=[[field.zero]*4 for _ in range(46)]
    for j,i in enumerate(cols):ET[i]=tt[j]
    ET[45]=tt[-1];mons=list(combinations_with_replacement(range(4),2));raw=[[field.zero]*len(mons) for _ in range(135)]
    for row,ff in enumerate(F):
        for mm,coef in ff.items():
            i,j=(*mm,*([45]*(2-len(mm))))
            if fixed in (i,j):continue
            c=cv(coef)
            for k,(u,v) in enumerate(mons):raw[row][k]+=c*(ET[i][u]*ET[j][v]+(ET[i][v]*ET[j][u] if u!=v else field.zero))
    quad=L*DM(raw,(135,10),field);assert quad.rank()==6;H=[]
    for row in quad.to_list():
        if not any(row):continue
        h=[[field.zero]*4 for _ in range(4)]
        for c,(i,j) in zip(row,mons):
            if i==j:h[i][i]+=2*c
            else:h[i][j]+=c;h[j][i]+=c
        H+=h
    Hess=DM(H,(len(H),4),field);assert Hess.rank()==3;v=dm(M(rec['actual_curve_tangent']));assert T*v==dm(tv) and (Hess*v).is_zero_matrix and v.rank()==1
    bsource=bv.subs(bs[9],a*bs[1]);BD=dm(bsource.jacobian(bs[:9]).row_join(acts.xreplace(dict(zip(bu,beta)))));wit=rec['restricted_beta_minor'];det=BD.extract(wit['rows'],wit['columns']).det();assert det==field.from_sympy(s.Integer(1289945088)) and BD.rank()==16
    return dict(field_polynomial=str(f),chosen_embedding='CRootOf(f,1), opposite to main implementation',both_embeddings_certified=True,raw135_full_rank=41,leading_rank=41,independent_cokernel_quadratic_rank=6,hessian_rank=3,actual_curve_is_common_radical=True,transverse_length=4,restricted_beta_rank=16,restricted_beta_minor='1289945088')

def run():
    start=time.time();normal=normal_minor();geo=review_geometry_and_points();paths=[Path(__file__),ROOT/'reports/S61_MINIMAL_M0_BASELINE_REMOVAL_20261101.json',ROOT/'code/s61_minimal_M0_baseline_removal_20261101.py',ROOT/'research/notes/S61_MINIMAL_M0_BASELINE_REMOVAL_20261101.md',ROOT/'code/s61_minimal_local_quartic_bound_20261101.py',ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json',ROOT/'reports/S61_MINIMAL_BASELINE_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'code/s61_minimal_junction_review_20261101.py']
    return dict(id='S61_MINIMAL_M0_BASELINE_REMOVAL_REVIEW_20261101',all_passed=True,normal_minor=normal,geometry_and_points=geo,accepted_all_finite_M0_baselines_removed_over_general_beta=True,scope='Both finite M0 baselines; parameter-infinity flag endpoints and other minimal tails remain separate. Length4 is for full S, not S/(t).',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
