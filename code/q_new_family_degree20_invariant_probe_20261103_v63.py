#!/usr/bin/env python3
"""One60s modular degree20 probe; full ambient span is checked by Weyl36."""
import argparse
import itertools
import json
from pathlib import Path
import random
import resource
import signal
import time
import numpy as np
import sympy as s
from q_new_family_degree16_invariant_probe_20261103_v63 import ROOT,PRIME,read,sha,matrix,mod,modmat,necklaces,echelon,character
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import sp_basis

STEM='Q_NEW_FAMILY_DEGREE20_INVARIANT_PROBE_20261103_v63'
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 d8path=ROOT/'reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
 d12path=ROOT/'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
 d16path=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63/samples.json'
 mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
 oldreview=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json'
 info,d8,d12,d16,payload=map(read,[invpath,d8path,d12path,d16path,mappath])
 assert read(oldreview)['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
 chars=character(list(map(tuple,info['Cartan_weights'])),20);dimension=chars[20]['invariant_dimension'];assert dimension==36
 O=matrix(info['symplectic_form16']);reps=[modmat(matrix(R)) for R in info['all_ten_representation_matrices']]
 gs=sp_basis();gram=s.Matrix([[s.trace(A*B) for B in gs] for A in gs]);om,oi,gi=modmat(O),modmat(O.inv()),modmat(gram.inv())
 w4,w6,w8,w10=[necklaces(n) for n in [4,6,8,10]];assert len(w10)==108
 b8=d8['independent_five_columns'];b12=d12['degree12_basis_columns'];b16=d16['selected_columns'];assert len(b8)==5 and len(b12)==9 and len(b16)==20 and max(b16)<36
 products=list(itertools.product(range(5),range(9)))
 defs=[dict(kind='trace10',word=list(w)) for w in w10]
 defs += [dict(kind='F4_times_basis16_trace8',basis_column=j,word=list(w8[j])) for j in b16]
 defs += [dict(kind='basis8_times_basis12',indices=[i,j],degree8_original_column=b8[i],degree12_trace_word=list(w6[b12[j]])) for i,j in products]
 assert len(defs)==173
 def invariants(beta):
  b=np.array(beta,dtype=np.int64).reshape(16,1);vv=np.hstack([om@R@b%PRIME for R in reps]);mu=(b.T@vv%PRIME).reshape(10,1)*pow(2,-1,PRIME)%PRIME
  co=gi@mu%PRIME;theta=sum((int(c)*R for c,R in zip(co[:,0],reps)),np.zeros((16,16),dtype=np.int64))%PRIME
  H=(theta+oi@vv%PRIME@gi%PRIME@vv.T)%PRIME;tables={():np.eye(16,dtype=np.int64)}
  for length in range(1,6):
   for w in itertools.product(range(2),repeat=length):tables[w]=tables[w[:-1]]@[theta,H][w[-1]]%PRIME
  def trace(w):
   k=len(w)//2;return int(np.sum(tables[w[:k]]*tables[w[k:]].T,dtype=np.int64)%PRIME)
  f4=int((mu.T@co)[0,0]%PRIME);v8=[trace(w) for w in w4]+[f4*f4%PRIME];v8=[v8[i] for i in b8];v12=[trace(w6[i]) for i in b12]
  return [trace(w) for w in w10]+[f4*trace(w8[i])%PRIME for i in b16]+[v8[i]*v12[j]%PRIME for i,j in products]
 basis=matrix(info['primitive24_by16_basis']);br=list(basis.T.rref()[1]);left=basis[br,:].inv()
 model=S61();z,u,beta=primitive_beta(model);one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
 coeff=[beta[i,j].coeff(t) for i,j in itertools.combinations(range(4),2) for t in z];betaU=s.Matrix([[s.expand(c).coeff(v) for v in u] for c in coeff]);convert=modmat(left*betaU[br,:])
 def actual(point):
  values=[]
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
  pivot=values[payload['normalization']['value_pivot_node']];pf=values[payload['normalization']['value_pfaffian_node']];assert pivot and pf
  lin=np.array([values[payload['output_u1_through_u61'][j]] for j in one],dtype=np.int64)
  return list(map(int,convert@lin%PRIME)),pivot,pf
 points=[];actual_beta=[];domain=[]
 for sample in d16['actual_QQ_points_and_betas']:
  point=sample['point'];b,pivot,pf=actual(point);assert b==[mod(v) for v in sample['beta16']]
  points.append(point);actual_beta.append(b);domain.append([pivot,pf])
 rng=random.Random(2026110320);seen={tuple(p[k] for k in payload['input_names']) for p in points}
 while len(points)<64:
  point={f'n{i}':rng.randint(1,13) for i in range(3)};point.update(scale=1);point.update({f'base{i}':rng.randint(-4,4) for i in range(1,5)})
  point.update({f'g{i}':v for i,v in enumerate([1,0,0,1,0,0,0,0,0,0])});key=tuple(point[k] for k in payload['input_names'])
  if key in seen:continue
  seen.add(key)
  try:b,pivot,pf=actual(point)
  except(ZeroDivisionError,ValueError,AssertionError):continue
  points.append(point);actual_beta.append(b);domain.append([pivot,pf])
 ambient=[[rng.randrange(PRIME) for _ in range(16)] for _ in range(64)]
 ae=[invariants(b) for b in ambient];re=[invariants(b) for b in actual_beta]
 ar,ac,_=echelon(ae);rr,rc,_=echelon(re);assert len(ar)<=dimension and len(rr)<=dimension
 from sympy.polys.matrices import DomainMatrix
 kernel=[]
 if len(ar)==dimension:
  M=s.Matrix([[row[j] for j in ac] for row in re]);kernel=[[int(v)%PRIME for v in row] for row in DomainMatrix.from_Matrix(M).convert_to(s.GF(PRIME)).nullspace().to_Matrix().tolist()]
 result=dict(id=STEM,all_passed=True,prime=PRIME,full_character_dimension20=dimension,character_records=chars,
  definitions=defs,number_invariants=173,ambient_count=64,actual_count=64,previous_QQ_point_count=20,
  ambient_rank=len(ar),actual_restriction_rank=len(rr),ambient_trace10_rank=len(echelon([row[:108] for row in ae])[0]),actual_trace10_rank=len(echelon([row[:108] for row in re])[0]),
  ambient_independent_rows=ar,ambient_independent_columns=ac,actual_independent_rows=rr,actual_independent_columns=rc,
  ambient_beta16=ambient,actual_source_points=points,actual_beta16=actual_beta,value_pivots_and_Pfaffians_modp=domain,
  ambient_evaluation=ae,actual_evaluation=re,restricted_modular_kernel_in_ambient_basis=kernel,
  complete_ambient_span_verified=len(ar)==dimension,
  actual_symbolic_identity_or_degree20_equation_claimed=False,exact_QQ_restriction_minor_computed=False,
  no_degree22_probe_or_Lie_GB=True,input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'code/q_new_family_degree16_invariant_probe_20261103_v63.py',invpath,d8path,d12path,d16path,mappath,oldreview]},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:result[k] for k in ['all_passed','full_character_dimension20','ambient_rank','actual_restriction_rank','ambient_trace10_rank','actual_trace10_rank','seconds']},indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
