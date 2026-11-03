#!/usr/bin/env python3
"""Exact QQ degree18 restriction minor after the single bounded modular probe."""
import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
import resource
import shlex
import subprocess
import sys
import time
import sympy as s
from q_new_family_degree16_invariant_probe_20261103_v63 import ROOT,PRIME,read,sha,matrix,mod
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import sp_basis

STEM='Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63'
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
def main(out,probe):
 sys.set_int_max_str_digits(0);started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 pp=probe/'report.json';data=read(pp);assert data['all_passed'] and data['ambient_trace9_rank']==data['actual_trace9_rank']==9
 for n,h in data['input_sha256'].items():assert sha(ROOT/n)==h,n
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
 samplepath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63/samples.json'
 info,payload,samples=map(read,[invpath,mappath,samplepath]);B=matrix(info['primitive24_by16_basis']);br=list(B.T.rref()[1]);left=B[br,:].inv()
 O=matrix(info['symplectic_form16']);reps=[matrix(M) for M in info['all_ten_representation_matrices']]
 gs=sp_basis();gram=s.Matrix([[s.trace(A*C) for C in gs] for A in gs])
 rows=data['actual_independent_rows'];cols=data['actual_independent_columns'];assert len(rows)==len(cols)==9 and max(cols)<60
 model=S61();z,u,beta=primitive_beta(model);one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
 coeff=[beta[i,j].coeff(t) for i,j in itertools.combinations(range(4),2) for t in z]
 betaU=s.Matrix([[s.expand(c).coeff(v) for v in u] for c in coeff]);points=[];betas=[]
 for row in rows:
  point=data['actual_QQ_source_points'][row];values=[]
  for node in payload['nodes']:
   op=node[0]
   if op=='q':v=Fraction(node[1])
   elif op=='v':v=Fraction(point[payload['input_names'][node[1]]])
   elif op=='+':v=sum((values[j] for j in node[1]),Fraction(0))
   elif op=='*':
    v=Fraction(1)
    for j in node[1]:v*=values[j]
   elif op=='inv':v=1/values[node[1]]
   else:raise AssertionError(op)
   values.append(v)
  b24=betaU*s.Matrix([s.Rational(values[payload['output_u1_through_u61'][i]]) for i in one]);b=left*b24[br,:]
  assert B*b==b24 and list(map(str,b))==samples['actual_QQ_points_and_betas'][row]['beta16']
  assert [mod(v) for v in b]==data['actual_beta16_modp'][row]
  vp=values[payload['normalization']['value_pivot_node']];pf=values[payload['normalization']['value_pfaffian_node']];assert vp and pf
  points.append(dict(point=point,beta16=list(map(str,b)),value_pivot=str(vp),value_pfaffian=str(pf)));betas.append(b)
 native=ROOT/'code/q_new_family_degree18_invariant_exact_20261103_v63.cpp';np=out/'native-input.txt'
 with np.open('w') as f:
  def write(M):
   f.write(f'{M.rows} {M.cols}\n')
   for row in M.tolist():f.write(' '.join(map(str,row))+'\n')
  f.write('9 60\n')
  for M in [O,O.inv(),gram.inv(),*reps]:write(M)
  for d in data['definitions'][:60]:assert d['kind']=='trace9';f.write(''.join(map(str,d['word']))+'\n')
  f.write(' '.join(map(str,cols))+'\n')
  for b in betas:write(b)
 binary=ROOT/'work/new-family-degree18-exact-native-v63';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
 subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 evidence=out/'exact-evaluation.json'
 with (out/'native.log').open('w') as f:subprocess.run([str(binary),str(np),str(evidence)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=65,preexec_fn=cap)
 exact=read(evidence);assert exact['all_passed'] and len(exact['evaluation'])==9
 for i,row in enumerate(rows):assert [mod(v) for v in exact['evaluation'][i]]==data['actual_evaluation'][row][:60]
 det=Fraction(exact['exact_selected9_minor']);assert det
 from sympy.polys.matrices import DomainMatrix
 modular=s.Matrix([[data['actual_evaluation'][r][c] for c in cols] for r in rows]);md=int(DomainMatrix.from_Matrix(modular).convert_to(s.GF(PRIME)).det())%PRIME
 assert md==mod(str(det))!=0
 (out/'samples.json').write_text(json.dumps(dict(probe_rows=rows,selected_columns=cols,selected_words=[data['definitions'][j]['word'] for j in cols],actual_QQ_source_points=points),indent=2)+'\n')
 invariant_review=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json';old=read(invariant_review)
 assert old['all_passed'] and old['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
 dimreview=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63_REVIEW.json';assert read(dimreview)['all_passed']
 files=[Path(__file__),native,pp,invpath,mappath,samplepath,invariant_review,dimreview,ROOT/'research/notes'/f'{STEM}.md']
 result=dict(id=STEM,all_passed=True,full_degree18_invariant_dimension=9,
  candidate_trace9_count=60,complete_ambient_Inv18_spanned_by_traces=True,
  actual_new_image_QQ_degree18_restriction_rank=9,exact_selected9_minor_nonzero=True,
  exact9_minor_file='exact-evaluation.json',exact9_minor_numerator_bits=abs(det.numerator).bit_length(),exact9_minor_denominator_bits=det.denominator.bit_length(),
  exact9_minor_mod32003=md,all540_QQ_trace9_values_verified_modp=True,
  all9_QQ_source_points_freshly_evaluated_on_actual_dense_map=True,
  complete_degree18_restriction_injective=True,defining_homogeneous_invariant_degree_lower_bound=20,
  no_equation_or_degree20_probe_or_global_fiber_exhaustion_claimed=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},
  evidence_sha256={p.name:sha(p) for p in [np,evidence,out/'samples.json',out/'native.log']},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS',STEM,result['seconds'],'minor bits',result['exact9_minor_numerator_bits'],result['exact9_minor_denominator_bits'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.probe.resolve())
