#!/usr/bin/env python3
"""One bounded modular degree18 probe on ambient and actual saved QQ points.

No sampled kernel is called an equation. Reuses accepted universal theta/H
conjugation and the complete Inv14 basis; all evaluations here are new.
"""
import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
import random
import resource
import signal
import time
import numpy as np
import sympy as s
from q_new_family_degree16_invariant_probe_20261103_v63 import ROOT,PRIME,read,sha,matrix,mod,modmat,necklaces,echelon
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import sp_basis

STEM='Q_NEW_FAMILY_DEGREE18_INVARIANT_PROBE_20261103_v63'
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 d8path=ROOT/'reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
 d14path=ROOT/'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
 samplepath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63/samples.json'
 mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
 oldreview=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json'
 charpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
 info,d8,d14,samples,payload=map(read,[invpath,d8path,d14path,samplepath,mappath])
 assert read(oldreview)['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
 chars=read(charpath)['character_records'];assert chars[18]['invariant_dimension']==9 and chars[14]['invariant_dimension']==3
 O=matrix(info['symplectic_form16']);reps=[modmat(matrix(v)) for v in info['all_ten_representation_matrices']]
 gs=sp_basis();gram=s.Matrix([[s.trace(A*B) for B in gs] for A in gs])
 om,oi,gi=modmat(O),modmat(O.inv()),modmat(gram.inv())
 w4,w9=necklaces(4),necklaces(9);assert len(w9)==60
 b8=d8['independent_five_columns'];pairs=d14['degree14_bracket_pairs'];b14=d14['degree14_basis_columns'];assert len(b14)==3
 defs=[dict(kind='trace9',word=list(w)) for w in w9]
 defs += [dict(kind='F4_times_basis14_bracket',degree14_basis_column=j,degree8_pair=pairs[j],degree8_basis_columns=b8) for j in b14]
 dv=np.stack([np.hstack([om@R[:,[k]]%PRIME for R in reps]) for k in range(16)])
 def invariants(beta):
  b=np.array(beta,dtype=np.int64).reshape(16,1)
  vv=np.hstack([om@R@b%PRIME for R in reps]);mu=(b.T@vv%PRIME).reshape(10,1)*pow(2,-1,PRIME)%PRIME
  co=gi@mu%PRIME;theta=sum((int(c)*R for c,R in zip(co[:,0],reps)),np.zeros((16,16),dtype=np.int64))%PRIME
  U=vv@gi%PRIME;H=(theta+oi@U%PRIME@vv.T)%PRIME;mats=[theta,H]
  table={():np.eye(16,dtype=np.int64)}
  for length in range(1,6):
   for w in itertools.product(range(2),repeat=length):table[w]=table[w[:-1]]@mats[w[-1]]%PRIME
  traces=[]
  for w in w9:traces.append(int(np.sum(table[w[:4]]*table[w[4:]].T,dtype=np.int64)%PRIME))
  f4=int((mu.T@co)[0,0]%PRIME);dc=gi@vv.T%PRIME
  dtheta=np.stack([sum((int(dc[j,k])*reps[j] for j in range(10)),np.zeros((16,16),dtype=np.int64))%PRIME for k in range(16)])
  dh=np.stack([(dtheta[k]+oi@dv[k]%PRIME@U.T+oi@U%PRIME@dv[k].T)%PRIME for k in range(16)])
  gradients=[]
  for idx in b8:
   if idx==6:gradients.append((4*f4*(vv@co%PRIME)).reshape(16)%PRIME);continue
   word=w4[idx];kk=[np.zeros((16,16),dtype=np.int64) for _ in range(2)]
   for r,letter in enumerate(word):kk[letter]=(kk[letter]+table[word[r+1:]+word[:r]])%PRIME
   grad=(np.einsum('ij,kji->k',kk[0],dtheta)+np.einsum('ij,kji->k',kk[1],dh))%PRIME
   gradients.append(grad)
  products=[]
  for j in b14:
   a,b=pairs[j];bracket=int((gradients[a]@oi%PRIME)@gradients[b]%PRIME)
   products.append(f4*bracket%PRIME)
  return traces+products
 # Full actual source evaluation, compared with already exact QQ beta16 values.
 basis=matrix(info['primitive24_by16_basis']);br=list(basis.T.rref()[1]);left=basis[br,:].inv()
 model=S61();z,u,beta=primitive_beta(model);one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
 coeff=[beta[i,j].coeff(t) for i,j in itertools.combinations(range(4),2) for t in z]
 betaU=s.Matrix([[s.expand(c).coeff(v) for v in u] for c in coeff]);convert=modmat(left*betaU[br,:])
 actual=[];points=[]
 for sample in samples['actual_QQ_points_and_betas']:
  point=sample['point'];values=[]
  for node in payload['nodes']:
   op=node[0]
   if op=='q':v=mod(node[1])
   elif op=='v':v=point[payload['input_names'][node[1]]]%PRIME
   elif op=='+':v=sum(values[j] for j in node[1])%PRIME
   elif op=='*':
    v=1
    for j in node[1]:v=v*values[j]%PRIME
   elif op=='inv':v=pow(values[node[1]],-1,PRIME)
   else:raise AssertionError(op)
   values.append(v)
  assert values[payload['normalization']['value_pivot_node']] and values[payload['normalization']['value_pfaffian_node']]
  lin=np.array([values[payload['output_u1_through_u61'][j]] for j in one],dtype=np.int64)
  b=list(map(int,convert@lin%PRIME));assert b==[mod(v) for v in sample['beta16']]
  actual.append(b);points.append(point)
 rng=random.Random(2026110318);ambient=[[rng.randrange(PRIME) for _ in range(16)] for _ in range(24)]
 ea=[invariants(b) for b in ambient];er=[invariants(b) for b in actual]
 ar,ac,_=echelon(ea);rr,rc,_=echelon(er);atrace=len(echelon([r[:60] for r in ea])[0]);rtrace=len(echelon([r[:60] for r in er])[0])
 assert len(ar)<=9 and len(rr)<=9
 # A possible restricted kernel is retained only in the ambient-independent
 # nine-dimensional basis; no sample relation is promoted to a polynomial.
 restricted=s.Matrix([[row[j] for j in ac] for row in er]);kernel=[]
 if len(ar)==9:
  from sympy.polys.matrices import DomainMatrix
  kernel=DomainMatrix.from_Matrix(restricted).convert_to(s.GF(PRIME)).nullspace().to_Matrix().tolist()
  kernel=[[int(v)%PRIME for v in row] for row in kernel]
 result=dict(id=STEM,all_passed=True,prime=PRIME,time_cap_seconds=60,complete_character_dimension18=9,
  number_invariants=len(defs),definitions=defs,ambient_count=len(ambient),actual_QQ_source_count=len(actual),
  ambient_rank=len(ar),actual_restriction_rank=len(rr),ambient_trace9_rank=atrace,actual_trace9_rank=rtrace,
  ambient_independent_rows=ar,ambient_independent_columns=ac,actual_independent_rows=rr,actual_independent_columns=rc,
  ambient_beta16=ambient,actual_QQ_source_points=points,actual_beta16_modp=actual,ambient_evaluation=ea,actual_evaluation=er,
  restricted_modular_kernel_in_ambient_basis=kernel,
  ambient_span_complete_over_QQ_if_rank9=len(ar)==9,
  candidate_identity_or_degree18_equation_claimed=False,exact_QQ_restriction_minor_not_yet_computed=True,
  no_higher_degree_search=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'code/q_new_family_degree16_invariant_probe_20261103_v63.py',invpath,d8path,d14path,samplepath,mappath,oldreview,charpath]},
  seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:result[k] for k in ['all_passed','ambient_rank','actual_restriction_rank','ambient_trace9_rank','actual_trace9_rank','seconds']},indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
