#!/usr/bin/env python3
"""A nonflat triple-branch infinity endpoint: rank7, transverse length12, t nilpotent."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_nil22_infinity_20261101 import beta_matrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import sp_basis,COORDS,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(x) for x in row] for row in M.tolist()]
def standards(gb):
    assert gb.is_zero_dimensional
    leads=[p.LM(order=gb.order).exponents for p in gb.polys];out=[];counts=[];n=len(gb.gens)
    for degree in range(20):
        current=[]
        for inds in combinations_with_replacement(range(n),degree):
            e=tuple(inds.count(i) for i in range(n))
            if not any(all(a<=b for a,b in zip(l,e)) for l in leads):current.append(e)
        out+=current;counts.append(len(current))
        if not current:return out,counts
    raise AssertionError('standard-monomial bound exceeded')
def local_certificate():
    m=S61();point=s.zeros(61,1);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    for j,i in enumerate(one):point[i]=j+2
    for i in [3,4,11,12]:point[i]=0
    point[13]=1;point[48]=0;point[49]=-5*point[34];point[50]=5*point[25];point[40]=-1
    r=calculate(m,True,point);assert r['homogenized_jacobian_rank']==39 and r['leading_jacobian_rank']==39
    beta={i:point[i] for i in one};_,polys=specialize(m,beta,True)
    cols=r['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
    J=s.Matrix(r['jacobian']).applyfunc(s.Rational);T=s.Matrix(r['tangent_basis']).applyfunc(s.Rational)[:,:]
    rows=list(J.T.to_DM().rref()[1]);pivots=list(J.to_DM().rref()[1]);assert len(rows)==len(pivots)==39
    U=s.zeros(45,39);iv=J.extract(rows,pivots).inv(method='DM')
    for i,c in enumerate(pivots):U[c,:]=iv[i,:]
    assert J[rows,:]*U==s.eye(39) and J*T==s.zeros(135,6)
    JU=J*U;mon2=list(combinations_with_replacement(range(6),2));mon3=list(combinations_with_replacement(range(6),3));i3={x:i for i,x in enumerate(mon3)}
    Q=s.zeros(135,21);terms=[]
    for row,f in enumerate(polys):
     ts=[]
     for mon,c in f.items():
      if fixed in mon:continue
      i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
      for n,(a,b) in enumerate(mon2):Q[row,n]+=c*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
     terms.append(ts)
    w2=-U*Q[rows,:];A2=Q-JU*Q[rows,:];Q3=s.zeros(135,56)
    for row,ts in enumerate(terms):
     for i,j,c in ts:
      for a in range(6):
       for n,mon in enumerate(mon2):Q3[row,i3[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w2[j,n]+T[j,a]*w2[i,n])
    w3=-U*Q3[rows,:];A3=Q3-JU*Q3[rows,:]

    a2rows=list(A2.T.to_DM().rref()[1]);a2cols=list(A2.to_DM().rref()[1]);assert len(a2rows)==9
    C=A2[:,a2cols]*A2.extract(a2rows,a2cols).inv(method='DM');assert C*A2[a2rows,:]==A2
    B=A3-C*A3[a2rows,:];brows=list(B.T.to_DM().rref()[1]);assert len(brows)==39
    v=s.symbols('v0:6');quadrics=[sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,mon2)) for row in A2[a2rows,:].tolist()]
    cubics=[sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,mon3)) for row in B[brows,:].tolist()]
    transverse_gb=s.groebner([f.subs(v[0],0) for f in quadrics],*v[1:],order='grevlex');basis,counts=standards(transverse_gb);assert counts==[1,5,6,2,0]
    cone_gb=s.groebner([f.subs(v[0],1) for f in quadrics+cubics],*v[1:],order='grevlex');cone_basis,cone_counts=standards(cone_gb);assert len(cone_basis)==7
    # Compute the full nonhomogeneous transverse quotient, retaining A2+A3.
    i2=[i for i,mon in enumerate(mon2) if 0 not in mon];i3=[i for i,mon in enumerate(mon3) if 0 not in mon];m2=[mon2[i] for i in i2];m3=[mon3[i] for i in i3];mons=[(i,) for i in range(1,6)]+m2+m3;ix={x:i for i,x in enumerate(mons)}
    relations=[]
    for row in range(135):
        relations.append([s.S.Zero]*5+list(A2[row,i2])+list(A3[row,i3]))
        for a in range(1,6):
            val=[s.S.Zero]*55
            for i,mon in enumerate(m2):val[ix[tuple(sorted((a,)+mon))]]+=A2[row,i2[i]]
            if any(val):relations.append(val)
    Rel=s.Matrix(relations);rank=len(Rel.to_DM().rref()[1]);assert rank==44
    t=list(T[-1,1:])+list(w2[-1,i2])+list(w3[-1,i3])
    def product(a,b):
        out=[s.S.Zero]*55
        for i,x in enumerate(a):
            if not x:continue
            for j,y in enumerate(b):
                if y and len(mons[i])+len(mons[j])<=3:out[ix[tuple(sorted(mons[i]+mons[j]))]]+=x*y
        return out
    t2=product(t,t);assert len(Rel.col_join(s.Matrix([t2])).to_DM().rref()[1])==rank
    assert len(Rel.col_join(s.Matrix([t])).to_DM().rref()[1])==rank+1
    extra=[t]
    for mon in [(i,) for i in range(1,6)]+m2:
        a=[s.S.Zero]*55;a[ix[mon]]=1;extra.append(product(t,a))
    assert len(Rel.col_join(s.Matrix(extra)).to_DM().rref()[1])==46
    return dict(local_data=r,jacobian_rows=rows,jacobian_columns=pivots,quadratic_monomials=mon2,cubic_monomials=mon3,A2=enc(A2),A3=enc(A3),quadratic_rows=a2rows,cubic_rows=brows,initial_quadrics=list(map(str,quadrics)),initial_cubics=list(map(str,cubics)),transverse_quadratic_Hilbert=counts,localized_cone_basis=[str(p.as_expr()) for p in cone_gb.polys],localized_cone_standard_monomials=cone_basis,localized_cone_length=7,transverse_relation_rank=rank,transverse_length=12,transverse_infinity_length=10,transverse_t_nonzero=True,transverse_t_square_zero=True,transverse_t_class=list(map(str,t)))

def universal_rank():
    m=S61();u=s.symbols('u1:62');beta={i:u[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};beta.update({3:0,4:0,11:0,12:0,48:0,49:-5*u[34],50:5*u[25]})
    two,F=specialize(m,beta,True);fixed=two.index(40);cols=[i for i in range(45) if i!=fixed];D=s.zeros(135,44);tc=s.zeros(135,1)
    for row,f in enumerate(F):
     for j,col in enumerate(cols):D[row,j]=-f.get(tuple(sorted((fixed,col))),0)
     tc[row]=-f.get((fixed,),0)
    D0=D[:70,:];rows=list(D0.T.rref()[1]);piv=list(D0.rref()[1]);assert len(rows)==len(piv)==24
    K=s.Matrix.hstack(*D0.nullspace());U=s.zeros(44,24);iv=D0.extract(rows,piv).inv()
    for j,i in enumerate(piv):U[i,:]=iv[j,:]
    shift=U*tc[rows,:];assert (D0*shift-tc[:70,:]).applyfunc(s.expand)==s.zeros(70,1)
    R=(D[70:,:]*K).row_join(tc[70:,:]-D[70:,:]*shift).applyfunc(s.expand);rr,pp=R.to_DM().rref();rr=rr.to_Matrix();free=[i for i in range(R.cols) if i not in pp];N=s.zeros(R.cols,len(free))
    for j,i in enumerate(free):
     N[i,j]=1
     for a,c in enumerate(pp):N[c,j]=-rr[a,i]
    assert (R.to_DM()*N.to_DM()).is_zero_matrix
    assert len(pp)==15
    return dict(beta={str(u[i]):str(v) for i,v in beta.items()},first_rank=24,residual_rank=15,full_rank_upper=39,constant_kernel=enc(K),t_lift=enc(shift),residual=enc(R),residual_kernel=enc(N),free=free)

def branch_geometry(local):
    m=S61();x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);one=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==1];two=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==2]
    Q=s.zeros(4);Q[1,2]=x*z;Q[1,3]=-x*w;Q=Q-Q.T
    M=beta_matrix(m,Q);H=s.Matrix.hstack(*M.nullspace());assert H.cols==10
    zz,bu,B=primitive_beta(m);vector=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vector.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();actions=s.Matrix.hstack(*[(left*infinitesimal_action(vector,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()])
    def coords(Q):
        P=X[0]**2*J0;P[1:,1:]+=Q;pp=m.normalize(P)['parameters'];return s.Matrix([pp[i] for i in two])
    q0=coords(Q);fixed=next(i for i,v in enumerate(q0) if v);chart=[i for i in range(45) if i!=fixed];qd=[]
    for A in sp_basis():
        dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T
        dq=coords(dQ);dq-=q0*dq[fixed]/q0[fixed];qd.append(dq[chart,:])
    orbit=s.Matrix.hstack(*qd);QD=s.zeros(44,10).row_join(orbit);point=s.Matrix(local['local_data']['point61']).applyfunc(s.Rational);sub=dict(zip(bu,[point[i] for i in one]));BD=H.row_join(actions.subs(sub));assert BD.rank()==16
    tangent=QD*s.Matrix.hstack(*BD.nullspace());piv=list(tangent.to_DM().rref()[1]);assert len(piv)==1
    T=s.Matrix(local['local_data']['tangent_basis']).applyfunc(s.Rational);rr=list(T.T.rref()[1]);vm=T[rr,:].inv()*tangent[:,piv].col_join(s.zeros(1,1))[rr,:]
    assert vm==s.Matrix([-s.Rational(3596,2605),-s.Rational(32,521),0,s.Rational(32,521),0,0])
    constraint=s.zeros(1,16);constraint[0,one.index(48)]=1;HH=H*s.Matrix.hstack(*(constraint*H).nullspace());assert HH.cols==9
    restricted=HH.row_join(actions.subs(sub));qrestricted=s.zeros(44,9).row_join(orbit);assert restricted.rank()==16 and restricted.col_join(qrestricted).rank()==16
    stabilizer=s.Matrix.hstack(*orbit.nullspace());assert orbit.rank()==7 and stabilizer.cols==3
    b=s.symbols('b0:9');stab_action=(actions*stabilizer).subs(dict(zip(bu,HH*s.Matrix(b))),simultaneous=True)
    assert (constraint*stab_action).applyfunc(s.expand)==s.zeros(1,3) and (M*stab_action).applyfunc(s.expand)==s.zeros(134,3)
    cols=list(restricted.to_DM().rref()[1]);minor=restricted[:,cols].det(method='domain-ge');assert minor!=0
    return dict(beta_kernel=enc(H),beta_differential=enc(BD),main_tangent=enc(vm),main_generic_length_lower_bound=2,restricted_beta_kernel=enc(HH),restricted_beta_differential=enc(restricted),restricted_joint_rank=16,restricted_beta_minor_columns=cols,restricted_beta_minor=str(minor),projective_orbit_rank=7,stabilizer_kernel=enc(stabilizer),stratum_dimension=16)

def verify():
    start=time.time();local=local_certificate();universal=universal_rank();geometry=branch_geometry(local)
    path=ROOT/'reports/S61_MINIMAL_ZERO49_BRANCHES_20261101.json';branches=json.loads(path.read_text());assert branches['all_passed']
    inputs=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_NONBINARY_FORMAL_RING_20261101.json',path,ROOT/'research/notes/S61_MINIMAL_ZERO49_BRANCHES_20261101.md']
    return dict(id='S61_NONBINARY_ZERO49_ENDPOINT_20261101',all_passed=True,field='QQ and QQ(nine beta parameters)',local=local,universal=universal,geometry=geometry,generic_branch_lengths=[2,3,2],completed_ring_rank_over_parameter=7,transverse_length=12,torsion_minimal_generator_count=5,not_Cohen_Macaulay=True,homogenizer_nilpotent=True,generic16D_stratum_removed=True,proof_note='research/notes/S61_NONBINARY_ZERO49_ENDPOINT_20261101.md',limitations=['Full t nilpotence exponent and complete multiplication are not determined.','Do not transfer transverse t-square-zero to the full local ring.','Generic16D stratum removal uses the finite-jet specialization and descent argument in section6; surrounding curve strata are not exhausted.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seconds=time.time()-start)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','completed_ring_rank_over_parameter','transverse_length','homogenizer_nilpotent','seconds']}))
if __name__=='__main__':main()
