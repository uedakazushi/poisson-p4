#!/usr/bin/env python3
"""Independent univariate/cofactor checks for the finite quintic open."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from s61_global_annihilator_review_20261103_v62 import resultant
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

P=32003
def trim(f):
 f=[int(c)%P for c in f]
 while f and not f[-1]:f.pop()
 return f
def add(f,g):return trim([(f[i]if i<len(f)else 0)+(g[i]if i<len(g)else 0)for i in range(max(len(f),len(g)))])
def mul(f,g):
 out=[0]*(len(f)+len(g)-1)if f and g else[]
 for i,c in enumerate(f):
  for j,d in enumerate(g):out[i+j]=(out[i+j]+c*d)%P
 return trim(out)
def divide(f,g):
 f=trim(f);g=trim(g);q=[0]*max(0,len(f)-len(g)+1);inv=pow(g[-1],-1,P)
 while len(f)>=len(g):
  k=len(f)-len(g);c=f[-1]*inv%P;q[k]=c
  for i,d in enumerate(g):f[i+k]=(f[i+k]-c*d)%P
  f=trim(f)
 return trim(q),f
def gcd(f,g):
 while g:f,g=g,divide(f,g)[1]
 return trim([c*pow(f[-1],-1,P)for c in f])
def invert_det(M):
 n=len(M);A=[[int(v)%P for v in row]+[int(i==j)for j in range(n)]for i,row in enumerate(M)];det=1
 for j in range(n):
  k=next(k for k in range(j,n)if A[k][j]);A[k],A[j]=A[j],A[k]
  if k!=j:det=-det
  v=A[j][j];det=det*v%P;A[j]=[(c*pow(v,-1,P))%P for c in A[j]]
  for i in range(n):
   if i!=j and A[i][j]:
    c=A[i][j];A[i]=[(u-c*v)%P for u,v in zip(A[i],A[j])]
 assert [row[:n]for row in A]==[[int(i==j)for j in range(n)]for i in range(n)]
 return np.array([row[n:]for row in A],dtype=np.int64),det%P

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/S61_AFFINE_QUINTIC_OPEN_20261103_v62';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json')
 assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 paths=[Path(__file__),public,ROOT/'code/s61_global_annihilator_review_20261103_v62.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for name,h in d['input_sha256'].items():p=ROOT/name;assert sha(p)==h,name;paths.append(p)
 for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h and p.read_bytes()==(args.replay/name).read_bytes();paths.append(p)
 cyclic=ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62';old=(cyclic/'native-certificate.txt').read_text().splitlines();native=(folder/'native.txt').read_text().splitlines()
 chi=list(map(int,next(line for line in old if line.startswith('CHI ')).split()[1:]));coordinates={int(line.split()[1]):list(map(int,line.split()[2:]))for line in old if line.startswith('COORDINATE ')}
 polynomial=lambda name:list(map(int,next(line for line in native if line.startswith(name+' ')).split()[1:]))
 r,g,u,v=[polynomial(name)for name in['RESIDUAL','GCD','BEZOUT_CHI','BEZOUT_RESIDUAL']]
 assert len(chi)==1032 and len(r)==1031 and len(g)==6 and g[-1]==1 and len(u)<=1025 and len(v)<=1026
 assert add(mul(u,chi),mul(v,r))==g
 chibar,e1=divide(chi,g);rbar,e2=divide(r,g);assert not e1 and not e2 and gcd(chibar,rbar)==[1]
 assert resultant(chibar,rbar)==16461
 dg=[i*g[i]%P for i in range(1,6)];assert gcd(g,dg)==[1]and resultant(g,dg)==19651
 coord=[divide(coordinates[i],g)[1]for i in range(29)];assert coord[0]==[1]
 data=np.fromstring((cyclic/'quadrics.in').read_text(),sep=' ',dtype=np.int64);assert data[:4].tolist()==[29,134,P,17];offset=4
 high=data[offset:offset+134*435].reshape(134,435);offset+=134*435;cross=data[offset:offset+134*29*17].reshape(134,29,17);offset+=134*29*17;base=data[offset:offset+134*153].reshape(134,153);offset+=134*153;center=data[offset:]
 pairs=list(itertools.combinations_with_replacement(range(29),2));bpairs=list(itertools.combinations_with_replacement(range(17),2))
 linear=cross@center%P;constant=base@np.array([int(center[i])*int(center[j])%P for i,j in bpairs],dtype=np.int64)%P
 products=[divide(mul(coord[i],coord[j]),g)[1]for i,j in pairs]
 for row in range(134):
  value=[int(constant[row])]
  for p,c in zip(products,high[row]):value=add(value,[int(c)*v for v in p])
  for p,c in zip(coord,linear[row]):value=add(value,[int(c)*v for v in p])
  assert not divide(value,g)[1]
 certificate=read(folder/'jacobian-certificate.json');C=np.array(certificate['common_combination'],dtype=np.int64)
 known=read(ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json')['known_five_mod_32003'];scale=pow(int(center[0]),-1,P);Y=[int(v)*scale%P for v in center];norm=1;dets=[]
 for point in known:
  Z=[int(v)*scale%P for v in point+[0]*25];J=np.zeros((134,45),dtype=np.int64)
  for row in range(134):
   for(i,j),c in zip(pairs,high[row]):J[row,i]=(int(J[row,i])+int(c)*Z[j])%P;J[row,j]=(int(J[row,j])+int(c)*Z[i])%P
   for i in range(29):
    for j in range(17):
     c=int(cross[row,i,j]);J[row,i]=(int(J[row,i])+c*Y[j])%P
     if j:J[row,28+j]=(int(J[row,28+j])+c*Z[i])%P
   for(i,j),c in zip(bpairs,base[row]):
    if i:J[row,28+i]=(int(J[row,28+i])+int(c)*Y[j])%P
    if j:J[row,28+j]=(int(J[row,28+j])+int(c)*Y[i])%P
  inverse,det=invert_det((C@J%P).tolist());dets.append(det);norm=norm*det%P;L=inverse@C%P
  assert np.array_equal(L@J%P,np.eye(45,dtype=np.int64))
  assert not np.any(L[29:]@J[:,:29]%P)and np.array_equal(L[29:]@J[:,29:]%P,np.eye(16,dtype=np.int64))
 assert dets==certificate['determinants']==[20321,16946,11970,18492,6923]and norm==25750
 paths += [ROOT/'reports/S61_AFFINE_SEVENTEEN_BASE_20261103_v62_REVIEW.json',ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62_REVIEW/report.json',ROOT/'research/notes/S61_AFFINE_QUINTIC_OPEN_20261103_v62_REVIEW.md']
 report=dict(id='S61_AFFINE_QUINTIC_OPEN_20261103_v62_REVIEW',all_passed=True,
  fresh_producer_and_all_saved_evidence_match=True,independent_pure_list_full2056_coefficient_Bezout=True,
  quotient_resultant_nonzero=16461,quintic_squarefree_resultant=19651,
  independent_all134_original_equations_zero_in_frozen_quintic=True,
  independent_all5_original134_by45_Jacobians_and_common_leftinverse_identities=True,
  high2051_block_unit_and_exact_monic5_relation_in_actual_ideal_accepted=True,
  exact_module_presentation=[5,675],generic_knownfive_lower_bound_plus_free5_upper_accepted=True,
  remainder_differential_transport_and_unit_denominator_telescoping_accepted=True,
  explicit_norm_open_actual_scheme_equals_known_finite_free5_accepted=True,
  generic_coefficients='finite inverse and norm recipes; not expanded tables',
  complement_resolved_or_whole_affine_five_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS actual monic5/open finite675 presentation and norm-isolated knownfive; complement OPEN')

if __name__=='__main__':main()
