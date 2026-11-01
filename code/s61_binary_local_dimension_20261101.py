#!/usr/bin/env python3
"""Exact characteristic-zero local dimension via cubic obstructions and a quartic minor."""
import sys,json,time,hashlib
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
import numpy as np
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

def verify():
    start=time.time();p=32003
    root=Path(__file__).resolve().parents[1];source_path=root/'reports/S61_BINARY_QUADRATIC_OBSTRUCTION_20261101.json';source=json.loads(source_path.read_text());assert source['all_passed'] and source['quadratic_span_rank']==5;pt=list(map(s.Rational,source['point']));m=S61();beta={i:pt[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};two,polys=specialize(m,beta,True)
    k=source['chart_fixed_q'];cols=[i for i in range(45) if i!=k];J=s.Matrix(source['jacobian']).applyfunc(s.Rational);T=s.Matrix(source['tangent']).applyfunc(s.Rational)
    def mod(v):v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,p)%p
    def arr(M):return np.array([[mod(x) for x in row] for row in M.tolist()],dtype=np.int64)
    def rref(A):
     A=A.copy()%p;row=0;piv=[]
     for col in range(A.shape[1]):
      choices=np.flatnonzero(A[row:,col])
      if len(choices)==0:continue
      j=row+int(choices[0]);A[[row,j]]=A[[j,row]];A[row]=A[row]*pow(int(A[row,col]),-1,p)%p
      for i in range(A.shape[0]):
       if i!=row and A[i,col]:A[i]=(A[i]-A[i,col]*A[row])%p
      piv.append(col);row+=1
      if row==A.shape[0]:break
     return A,piv
    j=arr(J);t=arr(T);_,C=rref(j);_,R=rref(j[:,C].T);assert len(C)==len(R)==35
    mi=np.concatenate([j[np.ix_(R,C)],np.eye(35,dtype=np.int64)],axis=1);rr,_=rref(mi);inv=rr[:,35:];assert np.array_equal(j[np.ix_(R,C)]@inv%p,np.eye(35,dtype=np.int64));U=np.zeros((45,35),dtype=np.int64);U[C,:]=inv;JU=j@U%p
    H=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(polys):
     for mon,c in f.items():
      if k in mon:continue
      ind=[cols.index(i) for i in mon]+[44]*(2-len(mon));i,l=ind;v=mod(c)
      if i==l:H[row,i,l]+=v
      else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p
    H%=p;mons2=list(combinations_with_replacement(range(10),2));mons3=list(combinations_with_replacement(range(10),3));idx3={v:i for i,v in enumerate(mons3)}
    Q=np.zeros((135,55),dtype=np.int64)
    for row in range(135):
     a=t.T@(H[row]@t%p)%p
     for l,(i,jj) in enumerate(mons2):Q[row,l]=a[i,jj]*(1 if i==jj else 2)%p
    w2=(-U@Q[R,:])%p;R2=(Q-JU@Q[R,:])%p;assert len(rref(R2)[1])==5
    Cubic=np.zeros((135,220),dtype=np.int64)
    for row in range(135):
     cross=2*t.T@(H[row]@w2%p)%p
     for i in range(10):
      for j2,mon in enumerate(mons2):Cubic[row,idx3[tuple(sorted((i,)+mon))]]+=cross[i,j2]
    Cubic%=p;R3=(Cubic-JU@Cubic[R,:])%p
    rr,piv=rref(np.concatenate([R2,R3],axis=1));assert sum(i<55 for i in piv)==5
    quadrics=rr[:5,:55];cubics=rr[5:len(piv),55:]
    # Add all variable multiples of quadrics to distinguish genuinely new cubic initial forms.
    mult=[]
    for q in quadrics:
     for i in range(10):
      row=np.zeros(220,dtype=np.int64)
      for mon,c in zip(mons2,q):row[idx3[tuple(sorted((i,)+mon))]]+=c
      mult.append(row%p)
    rankold=len(rref(np.array(mult))[1]);ranknew=len(rref(np.concatenate([np.array(mult),cubics]))[1]);assert (rankold,ranknew)==(50,65)
    # Pick three linear parameters and certify finite normal colength by quartic rank.
    hess=[]
    for q in quadrics:
     a=np.zeros((10,10),dtype=np.int64)
     for (i,jj),c in zip(mons2,q):
      a[i,jj]+=c*(2 if i==jj else 1)
      if i!=jj:a[jj,i]+=c
     hess.extend(a%p)
    rrh,ph=rref(np.array(hess));freeh=[i for i in range(10) if i not in ph]
    krad=np.zeros((10,len(freeh)),dtype=np.int64)
    for l,i in enumerate(freeh):
     krad[i,l]=1
     for row,jj in enumerate(ph):krad[jj,l]=-rrh[row,i]%p
    _,parameters=rref(krad.T);normal=[i for i in range(10) if i not in parameters]
    assert len(parameters)==3 and len(normal)==7
    normalmons={d:list(combinations_with_replacement(normal,d)) for d in (1,2,3,4)}
    normalindices={d:{m:i for i,m in enumerate(normalmons[d])} for d in (1,2,3,4)}
    normalq=np.array([[q[mons2.index(m)] for m in normalmons[2]] for q in quadrics],dtype=np.int64)
    normalc=np.array([[c[mons3.index(m)] for m in normalmons[3]] for c in cubics],dtype=np.int64)
    mat4=[]
    for d,gens in [(2,normalq),(3,normalc)]:
     for row in gens:
      for mult in normalmons[4-d]:
       rr=np.zeros(len(normalmons[4]),dtype=np.int64)
       for mon,c in zip(normalmons[d],row):rr[normalindices[4][tuple(sorted(mon+mult))]]+=c
       mat4.append(rr%p)
    mat4=np.array(mat4);_,p4=rref(mat4);_,rows4=rref(mat4.T);assert len(p4)==len(rows4)==210


    mat3=[]
    for row in normalq:
     for i in normal:
      rr=np.zeros(84,dtype=np.int64)
      for mon,c in zip(normalmons[2],row):rr[normalindices[3][tuple(sorted(mon+(i,)))]]+=c
      mat3.append(rr%p)
    mat3=np.concatenate([np.array(mat3),normalc]);assert len(rref(mat3)[1])==50
    assert len(rref(normalq)[1])==5
    minor=mat4[rows4,:].copy();det=1
    for col in range(210):
     i=next(i for i in range(col,210) if minor[i,col])
     if i!=col:minor[[i,col]]=minor[[col,i]];det=-det
     value=int(minor[col,col]);det=det*value%p;minor[col]=minor[col]*pow(value,-1,p)%p
     for i in range(col+1,210):
      if minor[i,col]:minor[i]=(minor[i]-minor[i,col]*minor[col])%p
    assert det
    return dict(id='S61_BINARY_LOCAL_DIMENSION_20261101',all_passed=True,coefficient_field='QQ with exact modular nonzero-minor certificates',prime=p,
     jacobian_pivot_columns=C,jacobian_pivot_rows=R,quadratic_monomials=mons2,cubic_monomials=mons3,
     residual_quadratic_matrix=R2.tolist(),residual_cubic_matrix=R3.tolist(),initial_quadrics=quadrics.tolist(),initial_cubics=cubics.tolist(),
     quadratic_rank=5,cubic_span_before=rankold,cubic_span_after=ranknew,
     linear_parameters=parameters,remaining_normal_variables=normal,normal_quadrics=normalq.tolist(),normal_cubics=normalc.tolist(),
     quartic_monomials=normalmons[4],quartic_matrix=mat4.tolist(),quartic_minor_rows=rows4,quartic_minor_determinant_mod_prime=det,
     quartic_rank=210,normal_hilbert_function_upper_bound=[1,7,23,34],normal_colength_upper_bound=65,local_dimension=3,
     proof_note='research/notes/S61_BINARY_LOCAL_DIMENSION_20261101.md',
     limitations=['The 135-equation unsaturated homogenization at one rational binary point.','The full local ring and t-saturation are not determined.','No full affine-fiber finiteness or exhaustion.'],
     source_sha256={str(q.relative_to(root)):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),source_path,root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']})


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,local_dimension=3,quartic_rank=210,normal_colength_upper_bound=65)))
if __name__=='__main__':main()
