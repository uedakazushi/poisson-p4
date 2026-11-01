#!/usr/bin/env python3
"""A dominant eighteen-dimensional binary rank-drop stratum and its local open.

The two displayed linear hypersurface charts are symplectically equivalent.
No assertion removes the entire generic-beta fiber of this stratum.
"""
import sys,json,time,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import sympy as s,numpy as np
from itertools import combinations_with_replacement
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta


def verify_endpoint():
    extra_zero="u12"
    start=time.time();p=32003;m=S61()
    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';r=json.loads(family_path.read_text());assert r['all_passed'];r['sample'][extra_zero]='0';r['sample']['u35']='11';r['sample']['u50']='-55'
    b={int(k[1:])-1:s.Rational(v) for k,v in r['sample'].items()};two,F=specialize(m,b,True);fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];q0=s.zeros(45,1);q0[fixed]=1
    J=s.zeros(135,45)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if len(mon)==2:
       i,j=mon
       if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
       if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
      elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
    T=s.Matrix.hstack(*J.nullspace());n=T.cols;assert n==11 and J*T==s.zeros(135,11) and T[-1,:]==s.zeros(1,11)
    _,piv=J.to_DM().rref();free=[i for i in range(45) if i not in piv];assert T[free,:]==s.eye(n)
    z,uu,bf=primitive_beta(m);bf=bf.subs({s.Symbol(k):s.Rational(v) for k,v in r['sample'].items()});a,c,d,h,k=s.symbols('a c d h k');v1=s.Matrix([a,1,0,c]);v2=s.Matrix([-c,0,1,d]);origin={a:0,c:0,d:0,h:0,k:0}
    sec=s.Matrix([(v1.T*bf.subs(dict(zip(z,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(z,v2)))*v2)[0]])
    assert sec.subs(origin)==s.zeros(2,1)
    curve_jac=sec.jacobian([a,c,d]).subs(origin);assert curve_jac.rank()==2
    direction=curve_jac.nullspace()[0]
    x,y,zz,w=X[1:];ell1=x-a*y+c*zz;ell2=w-c*y-d*zz;fam=(ell1*ell2+h*ell1**2+k*ell2**2)*(v1*v2.T-v2*v1.T)
    ds=[]
    for v in [a,c,d,h,k]:
     P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin);qs=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);qs-=q0*qs[fixed];ds.append(qs[cols,:].col_join(s.zeros(1,1)))
    family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4]);coordinates=family[free,:];assert T*coordinates==family and J*family==s.zeros(135,3)
    _,params=coordinates.T.to_DM().rref();params=list(params);normal=[i for i in range(n) if i not in params];assert params==[0,1,3] and len(normal)==8
    def mod(v):v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,p)%p
    def arr(M):return np.array([[mod(v) for v in row] for row in M.tolist()],dtype=np.int64)
    def rref(A):
     A=A.copy()%p;row=0;piv=[]
     for col in range(A.shape[1]):
      choices=np.flatnonzero(A[row:,col])
      if len(choices)==0:continue
      j=row+int(choices[0]);A[[row,j]]=A[[j,row]];A[row]=A[row]*pow(int(A[row,col]),-1,p)%p
      mask=np.flatnonzero(A[:,col]);mask=mask[mask!=row]
      if len(mask):A[mask]=(A[mask]-A[mask,col,None]*A[row,None,:])%p
      piv.append(col);row+=1
      if row==A.shape[0]:break
     return A,piv
    j=arr(J);t=arr(T);_,C=rref(j);_,R=rref(j[:,C].T);rank=len(C);rr,_=rref(np.concatenate([j[np.ix_(R,C)],np.eye(rank,dtype=np.int64)],axis=1));U=np.zeros((45,rank),dtype=np.int64);U[C,:]=rr[:,rank:];JU=j@U%p
    assert len(C)==len(R)==34 and np.array_equal(j[R,:]@U%p,np.eye(34,dtype=np.int64))
    H=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      i,l=[cols.index(v) for v in mon]+[44]*(2-len(mon));v=mod(c)
      if i==l:H[row,i,l]+=v
      else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p

    # Exact QQ quadratic rank in the eight transverse variables.
    Tn=T[:,normal];exact_mons=list(combinations_with_replacement(range(8),2));em={mon:i for i,mon in enumerate(exact_mons)};ecache={}
    def eprod(i,j):
     key=tuple(sorted((i,j)))
     if key not in ecache:
      out=[s.S.Zero]*36
      for a,v in enumerate(Tn[i,:]):
       if not v:continue
       for bb,w in enumerate(Tn[j,:]):
        if w:out[em[tuple(sorted((a,bb)))]]+=v*w
      ecache[key]=out
     return ecache[key]
    Qe=[]
    for f in F:
     out=[s.S.Zero]*36
     for mon,c in f.items():
      if fixed in mon:continue
      i,jj=[cols.index(v) for v in mon]+[44]*(2-len(mon))
      for k,v in enumerate(eprod(i,jj)):out[k]+=c*v
     Qe.append(out)
    coe=s.Matrix.hstack(*J.T.nullspace()).T;Ae=(coe.to_DM()*s.Matrix(Qe).to_DM()).to_Matrix();assert len(Ae.to_DM().rref()[1])==14

    H%=p;t=t[:,normal];n=len(normal)
    m2=list(combinations_with_replacement(range(n),2));m3=list(combinations_with_replacement(range(n),3));m4=list(combinations_with_replacement(range(n),4));i3={mon:i for i,mon in enumerate(m3)};i4={mon:i for i,mon in enumerate(m4)}
    Q=np.zeros((135,len(m2)),dtype=np.int64)
    for row in range(135):
     q=t.T@(H[row]@t%p)%p
     for k,(a,c) in enumerate(m2):Q[row,k]=q[a,c]*(1 if a==c else 2)%p
    w2=-U@Q[R,:]%p;A2=(Q-JU@Q[R,:])%p
    C3=np.zeros((135,len(m3)),dtype=np.int64)
    for row in range(135):
     cross=2*t.T@(H[row]@w2%p)%p
     for i in range(n):
      for k,mon in enumerate(m2):C3[row,i3[tuple(sorted((i,)+mon))]]+=cross[i,k]
    C3%=p;w3=-U@C3[R,:]%p;A3=(C3-JU@C3[R,:])%p
    C4=np.zeros((135,len(m4)),dtype=np.int64)
    for row in range(135):
     cross=2*t.T@(H[row]@w3%p)%p;square=w2.T@(H[row]@w2%p)%p
     for i in range(n):
      for k,mon in enumerate(m3):C4[row,i4[tuple(sorted((i,)+mon))]]+=cross[i,k]
     for i,mi in enumerate(m2):
      for k,mk in enumerate(m2):C4[row,i4[tuple(sorted(mi+mk))]]+=square[i,k]
    C4%=p;w4=-U@C4[R,:]%p;A4=(C4-JU@C4[R,:])%p
    base=np.concatenate([A2,A3,A4],axis=1);ncols=base.shape[1];rows=list(base)
    for row2,row3 in zip(A2,A3):
     for a in range(n):
      out=np.zeros(ncols,dtype=np.int64)
      for mon,v in zip(m2,row2):out[len(m2)+i3[tuple(sorted((a,)+mon))]]+=v
      for mon,v in zip(m3,row3):out[len(m2)+len(m3)+i4[tuple(sorted((a,)+mon))]]+=v
      if np.any(out):rows.append(out%p)
     for ma in m2:
      out=np.zeros(ncols,dtype=np.int64)
      for mon,v in zip(m2,row2):out[len(m2)+len(m3)+i4[tuple(sorted(ma+mon))]]+=v
      if np.any(out):rows.append(out%p)
    rels=np.array(rows);rr,piv=rref(rels);rk=len(piv);length=1+n+ncols-rk
    trow=np.concatenate([w2[-1,:],w3[-1,:],w4[-1,:]])%p
    rankt=len(rref(np.concatenate([rr[:rk,:],trow[None,:]],axis=0))[1])
    tmult=[trow]
    for a in range(n):
     out=np.zeros(ncols,dtype=np.int64)
     for mon,v in zip(m2,w2[-1,:]):out[len(m2)+i3[tuple(sorted((a,)+mon))]]+=v
     for mon,v in zip(m3,w3[-1,:]):out[len(m2)+len(m3)+i4[tuple(sorted((a,)+mon))]]+=v
     tmult.append(out%p)
    for ma in m2:
     out=np.zeros(ncols,dtype=np.int64)
     for mon,v in zip(m2,w2[-1,:]):out[len(m2)+len(m3)+i4[tuple(sorted(ma+mon))]]+=v
     tmult.append(out%p)
    rankideal=len(rref(np.concatenate([rr[:rk,:],np.array(tmult)],axis=0))[1])
    length3=1+n+len(m2)-len(rref(rels[:,:len(m2)])[1]);length4=1+n+len(m2)+len(m3)-len(rref(rels[:,:len(m2)+len(m3)])[1])

    assert length==65 and rk==430 and rankt==431 and rankideal==431
    assert [1,n,length3-1-n,length4-length3,length-length4]==[1,8,22,34,0]
    # Save nonzero minors of the actual nonhomogeneous relation matrix, and of
    # that matrix with t appended. Their characteristic-zero use is in the note.
    _,selected_rows=rref(rels.T);assert len(selected_rows)==430
    selected_columns=piv
    minor430=rels[np.ix_(selected_rows,selected_columns)].copy()
    def determinant_mod(M):
     M=M.copy();det=1
     for col in range(len(M)):
      ii=col+int(np.flatnonzero(M[col:,col])[0])
      if ii!=col:M[[col,ii]]=M[[ii,col]];det=-det
      value=int(M[col,col]);det=det*value%p;M[col]=M[col]*pow(value,-1,p)%p
      mask=np.flatnonzero(M[col+1:,col])+col+1
      if len(mask):M[mask]=(M[mask]-M[mask,col,None]*M[col,None,:])%p
     return det%p
    det430=determinant_mod(minor430);assert det430

    # Initial cubics lift to QQ because the exact and modular quadratic ranks agree.
    initial,ipiv=rref(np.concatenate([A2,A3],axis=1));assert sum(i<36 for i in ipiv)==14
    initialq=initial[:14,:36];initialc=initial[14:len(ipiv),36:]
    ms={d:list(combinations_with_replacement(range(8),d)) for d in range(6)};idx5={mon:i for i,mon in enumerate(ms[5])}
    mat5=[]
    for d,gens in [(2,initialq),(3,initialc)]:
     for row in gens:
      for mul in ms[5-d]:
       out=np.zeros(792,dtype=np.int64)
       for mon,c in zip(ms[d],row):out[idx5[tuple(sorted(mon+mul))]]+=c
       mat5.append(out%p)
    mat5=np.array(mat5);_,piv5=rref(mat5);_,rows5=rref(mat5.T);assert len(piv5)==len(rows5)==792
    quintic_minor=determinant_mod(mat5[rows5,:]);assert quintic_minor

    aug=np.concatenate([rels[selected_rows,:],trow[None,:]],axis=0)
    _,augcols=rref(aug);assert len(augcols)==431
    det431=determinant_mod(aug[:,augcols]);assert det431
    # A regular binary point on the SAME fixed-beta threefold: Q1=xw+x^2+w^2.
    P=X[0]**2*J0;xx,yy,zz,ww=X[1:];P[2,3]+=xx*ww+xx*xx+ww*ww;P[3,2]-=xx*ww+xx*xx+ww*ww
    q1=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);assert q1[fixed]==1
    J1=s.zeros(135,45)
    for row,f in enumerate(F):
     value=0
     for mon,c in f.items():
      if len(mon)==2:
       i,jj=mon;value+=c*q1[i]*q1[jj]
       if i!=fixed:J1[row,cols.index(i)]+=c*q1[jj]
       if jj!=fixed:J1[row,cols.index(jj)]+=c*q1[i]
      elif len(mon)==1:J1[row,44]+=c*q1[mon[0]]
     assert value==0
    T1=s.Matrix.hstack(*J1.nullspace());assert T1.shape==(45,10) and T1[-1,:]==s.zeros(1,10)
    mons1=list(combinations_with_replacement(range(10),2));idx1={mon:i for i,mon in enumerate(mons1)};cache={}
    def product1(i,j):
     key=tuple(sorted((i,j)))
     if key not in cache:
      vals=[s.S.Zero]*55
      for a,v in enumerate(T1[i,:]):
       if not v:continue
       for bb,w in enumerate(T1[j,:]):
        if w:vals[idx1[tuple(sorted((a,bb)))]]+=v*w
      cache[key]=vals
     return cache[key]
    Q1=[]
    for f in F:
     row=[s.S.Zero]*55
     for mon,c in f.items():
      if fixed in mon:continue
      i,jj=[cols.index(v) for v in mon]+[44]*(2-len(mon))
      for k,v in enumerate(product1(i,jj)):row[k]+=c*v
     Q1.append(row)
    co1=s.Matrix.hstack(*J1.T.nullspace()).T;A1=(co1.to_DM()*s.Matrix(Q1).to_DM()).to_Matrix();_,goodrows=A1.T.to_DM().rref();assert len(goodrows)==5
    hess=[]
    for row in A1[list(goodrows),:].tolist():
     M=s.zeros(10)
     for c,(i,jj) in zip(row,mons1):M[i,jj]+=c;M[jj,i]+=c
     hess+=M.tolist()
    assert len(s.Matrix(hess).to_DM().rref()[1])==7
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    return dict(extra_zero=extra_zero,sample=r['sample'],linear_parameters=params,normal_variables=normal,
     actual_binary_family_tangent=enc(coordinates),parameter_minor=str(coordinates[params,:].det()),curve_derivative=enc(curve_jac),
     jacobian_pivot_rows=R,jacobian_pivot_columns=C,exact_transverse_quadratic_rank=14,
     quintic_initial_minor_rows=rows5,quintic_initial_minor_mod_prime=quintic_minor,
     residual_quadratic_mod_prime=A2.tolist(),residual_cubic_mod_prime=A3.tolist(),residual_quartic_mod_prime=A4.tolist(),
     homogenizer_quadratic_mod_prime=w2[-1,:].tolist(),homogenizer_cubic_mod_prime=w3[-1,:].tolist(),homogenizer_quartic_mod_prime=w4[-1,:].tolist(),
     relation_matrix_shape=list(rels.shape),relation_minor_rows=selected_rows,relation_minor_columns=selected_columns,relation_minor_mod_prime=det430,
     augmented_homogenizer_minor_columns=augcols,augmented_homogenizer_minor_mod_prime=det431,
     same_beta_regular_point_q=list(map(str,q1)),same_beta_regular_full_rank=35,same_beta_regular_quadratic_rank=5,same_beta_regular_hessian_radical=3,
     transverse_length=65,transverse_hilbert_function=[1,8,22,34],completed_ring_rank=65,nilradical_fourth_power_zero=True,
     homogenizer_nonzero=True,homogenizer_square_zero=True,homogenizer_annihilates_nilradical=True,completed_infinity_quotient_rank=64,
     seconds=time.time()-start)


def verify():
    from weight1_structure import J as J4
    from itertools import combinations
    start=time.time();point=verify_endpoint();m=S61();u=s.symbols('u1:62')
    weights=m.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    beta={i:(s.S.Zero if i in [11,12] else u[i]) for i in one};two,F=specialize(m,beta,True);fixed=two.index(29);q0=s.zeros(45,1);q0[fixed]=1
    D=s.zeros(135,45);tc=s.zeros(135,1)
    for row,f in enumerate(F):
        val=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;val+=c*q0[i]*q0[j];D[row,i]+=c*q0[j];D[row,j]+=c*q0[i]
            elif len(mon)==1:tc[row]+=c*q0[mon[0]]
        assert val==0
    K=s.Matrix.hstack(*D[:70,:].nullspace());assert K.shape==(45,29)
    E=(D[70:,:]*K).applyfunc(s.expand)
    rr=[[0,1,3,5],[2,4,11,13,14,17,23,33,44,47,53,63],[9,10,12,15,16,20,22,29,32,43,45,46,50,52,59,62,64]]
    cc=[[0,1,2,5,7,10,12,14,18,19,22,28],[4,9,13,15,17,20,21,25,27],[3,16,23,24,26]]
    rr.append([i for i in range(65) if all(i not in r for r in rr)]);cc.append([i for i in range(29) if all(i not in c for c in cc)])
    for i in range(4):
        for j in range(4):
            if i>j or i==3 or j==3 or (i,j)==(0,2):assert E[rr[i],cc[j]]==s.zeros(len(rr[i]),len(cc[j]))
    core={u[3],u[13],u[34],u[35],u[48],u[49]}
    for i in range(3):assert E[rr[i],cc[i]].free_symbols<=core
    # Upper triangular blocks have sizes4x12,12x9,17x5. On ell=0, B has
    # a universal nonzero polynomial kernel vector, hence rank at most8.
    a,b,c,d,e=u[3],u[13],u[34],u[35],u[48];sub={u[49]:-5*c};B=E[rr[1],cc[1]].subs(sub)
    ker=s.Matrix([0,-6*b**3-3*b**2*e+8*b*c*a+3*b*e**2+2*c*a*e-d*a**2,
        -2*(2*b*c*d-16*c**3-4*c*d*e+d**2*a),0,0,
        -2*(12*b**2*c-6*b*c*e+2*b*d*a-8*c**2*a-d*a*e),
        2*(-2*b**2*d-8*b*c**2-b*d*e+4*c**2*e+2*c*d*a+d*e**2),0,0])
    assert (B*ker).applyfunc(s.expand)==s.zeros(12,1) and ker!=s.zeros(9,1)
    sample={s.Symbol(k):s.Rational(v) for k,v in point['sample'].items()}
    assert len(B.subs(sample).to_DM().rref()[1])==8
    Ds=D.subs(sample);chart=[i for i in range(45) if i!=fixed];full=Ds[:,chart].row_join(tc.subs(sample))
    assert len(Ds.to_DM().rref()[1])==33 and len(full.to_DM().rref()[1])==34

    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';rank_path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    family=json.loads(family_path.read_text());rank=json.loads(rank_path.read_text());z,uu,form=primitive_beta(m)
    actions=s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    stabilizer=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    constraints=s.Matrix([u[11],u[12],5*u[34]+u[49]]);constraint_matrix=constraints.jacobian(uu)
    hsub={u[11]:0,u[12]:0,u[49]:-5*u[34]}
    assert (constraint_matrix*stabilizer).subs(hsub).applyfunc(s.expand)==s.zeros(3,5)
    normal=(constraint_matrix*actions).subs(sample);columns=list(normal.to_DM().rref()[1]);assert len(columns)==3
    det=normal[:,columns].det();assert det==-18300

    # A finite symplectic involution exchanges the two hypersurface charts.
    S=s.Matrix([[0,0,0,1],[0,0,-1,0],[0,-1,0,0],[1,0,0,0]])
    assert S.T*J4*S==J4 and S*S==s.eye(4)
    def vector(M):return s.Matrix([s.expand(M[i,j]).coeff(z[k]) for i,j in combinations(range(4),2) for k in range(4)])
    basis=vector(form).jacobian(uu);rows=list(basis.T.to_DM().rref()[1]);changed=S.T*form.subs(dict(zip(z,S*s.Matrix(z))),simultaneous=True)*S
    new=(basis[rows,:].inv()*vector(changed)[rows,:]).applyfunc(s.expand)
    assert new[uu.index(u[11])]==u[12] and new[uu.index(u[12])]==u[11]
    assert s.expand(5*new[uu.index(u[34])]+new[uu.index(u[49])])==5*u[13]-u[48]
    Q=s.zeros(4);Q[1,2]=z[0]*z[3];Q[2,1]=-Q[1,2]
    assert (S*Q.subs(dict(zip(z,S*s.Matrix(z))),simultaneous=True)*S.T+Q).applyfunc(s.expand)==s.zeros(4)
    point.update(joint_image_dimension=18,generic_beta_stratum_dimension=2,first_jet_dominance_minor=str(det),
        first_jet_dominance_columns=columns,whole_generic_beta_stratum_removed=False)
    paths=[Path(__file__),family_path,rank_path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
        ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',ROOT/'research/notes/S61_BINARY_UNIFORM_CUBIC_20261101.md',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    return dict(id='S61_BINARY_HYPERSURFACE_RANK_DROP_20261101',all_passed=True,prime=32003,point=point,
        hypersurface_equation=str(constraints[2]),equivalent_hypersurface_equation=str(5*u[13]-u[48]),
        constant_first70_kernel=enc(K),residual_block_rows=rr,residual_block_columns=cc,
        diagonal_blocks=[enc(E[rr[i],cc[i]]) for i in range(3)],B_kernel=list(map(str,ker)),
        universal_residual_rank_upper_bound=17,generic_leading_rank=33,generic_full_rank=34,
        beta_normal_action=enc(normal),projective_stabilizer_preserves_hypersurface=True,
        symmetry_matrix=enc(S),symmetry_hypersurface_image=str(5*u[13]-u[48]),
        old_seventeen_dimensional_substratum='u35=0 inside 5*u35+u50=0',
        proof_note='research/notes/S61_BINARY_HYPERSURFACE_RANK_DROP_20261101.md',
        limitations=['The local free65 open does not remove the whole generic-beta surface of this18-dimensional stratum.',
                    'Its complement may have dimension17 and dominate W16; no full binary rank-drop exhaustion.',
                    'No full affine-fiber finiteness or five-point exhaustion.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    q=r['point'];print(json.dumps(dict(all_passed=True,joint_image_dimension=18,generic_fixed_beta_dimension=2,
        point={k:q[k] for k in ['first_jet_dominance_minor','quintic_initial_minor_mod_prime','relation_minor_mod_prime','augmented_homogenizer_minor_mod_prime']},seconds=r['seconds'])))


if __name__=='__main__':main()
