#!/usr/bin/env python3
"""One exact QQ36minor; no22+probe after the single degree20modular probe."""
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
STEM='Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63'
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
def main(out,probe):
 sys.set_int_max_str_digits(0);started=time.monotonic();out.mkdir(parents=True,exist_ok=False);pp=probe/'report.json';data=read(pp)
 assert data['all_passed'] and data['full_character_dimension20']==data['ambient_trace10_rank']==data['actual_trace10_rank']==36
 for n,h in data['input_sha256'].items():assert sha(ROOT/n)==h,n
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
 info,payload=map(read,[invpath,mappath]);basis=matrix(info['primitive24_by16_basis']);br=list(basis.T.rref()[1]);left=basis[br,:].inv()
 O=matrix(info['symplectic_form16']);reps=[matrix(M) for M in info['all_ten_representation_matrices']];gs=sp_basis();gram=s.Matrix([[s.trace(A*B) for B in gs] for A in gs])
 rows=data['actual_independent_rows'];cols=data['actual_independent_columns'];assert len(rows)==len(cols)==36 and max(cols)<108
 model=S61();z,u,beta=primitive_beta(model);one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
 coeff=[beta[i,j].coeff(t) for i,j in itertools.combinations(range(4),2) for t in z];betaU=s.Matrix([[s.expand(c).coeff(v) for v in u] for c in coeff])
 points=[];betas=[]
 for row in rows:
  point=data['actual_source_points'][row];values=[]
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
  b24=betaU*s.Matrix([s.Rational(values[payload['output_u1_through_u61'][i]]) for i in one]);b=left*b24[br,:];assert basis*b==b24
  assert [mod(v) for v in b]==data['actual_beta16'][row]
  vp=values[payload['normalization']['value_pivot_node']];pf=values[payload['normalization']['value_pfaffian_node']];assert vp and pf
  points.append(dict(point=point,beta16=list(map(str,b)),value_pivot=str(vp),value_pfaffian=str(pf)));betas.append(b)
 native=ROOT/'code/q_new_family_degree20_invariant_exact_20261103_v63.cpp';np=out/'native-input.txt'
 with np.open('w') as f:
  def write(M):
   f.write(f'{M.rows} {M.cols}\n')
   for row in M.tolist():f.write(' '.join(map(str,row))+'\n')
  f.write('36 108\n')
  for M in [O,O.inv(),gram.inv(),*reps]:write(M)
  for d in data['definitions'][:108]:assert d['kind']=='trace10';f.write(''.join(map(str,d['word']))+'\n')
  f.write(' '.join(map(str,cols))+'\n')
  for b in betas:write(b)
 binary=ROOT/'work/new-family-degree20-exact-native-v63';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
 subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 evidence=out/'exact-evaluation.json'
 with (out/'native.log').open('w') as f:subprocess.run([str(binary),str(np),str(evidence)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=65,preexec_fn=cap)
 exact=read(evidence);assert exact['all_passed'] and len(exact['evaluation'])==36
 for i,row in enumerate(rows):assert [mod(v) for v in exact['evaluation'][i]]==data['actual_evaluation'][row][:108]
 det=Fraction(exact['exact_selected36_minor']);assert det
 from sympy.polys.matrices import DomainMatrix
 modular=s.Matrix([[data['actual_evaluation'][r][c] for c in cols] for r in rows]);md=int(DomainMatrix.from_Matrix(modular).convert_to(s.GF(PRIME)).det())%PRIME
 assert md==mod(str(det))!=0
 (out/'samples.json').write_text(json.dumps(dict(probe_rows=rows,selected_columns=cols,selected_words=[data['definitions'][j]['word'] for j in cols],actual_QQ_source_points=points),indent=2)+'\n')
 invreview=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json';assert read(invreview)['all_passed']
 prev18=ROOT/'reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63/report.json';assert read(prev18)['all_passed']
 files=[Path(__file__),native,pp,invpath,mappath,invreview,prev18,ROOT/'code/q_new_family_degree16_invariant_probe_20261103_v63.py',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes'/f'{STEM}.md']
 result=dict(id=STEM,all_passed=True,full_degree20_invariant_dimension=36,all108_trace10_span_complete_Inv20=True,
  actual_new_image_QQ_degree20_restriction_rank=36,exact_selected36_minor_nonzero=True,exact36_minor_file='exact-evaluation.json',
  exact36_minor_numerator_bits=abs(det.numerator).bit_length(),exact36_minor_denominator_bits=det.denominator.bit_length(),exact36_minor_mod32003=md,
  all3888_QQ_trace10_values_verified_modp=True,all36_QQ_source_points_freshly_evaluated_on_actual_dense_map=True,
  complete_degree20_restriction_injective=True,defining_homogeneous_invariant_degree_lower_bound=22,
  exact_native_run_count=1,no_equation_or_degree22_probe_or_global_fiber_exhaustion_claimed=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},evidence_sha256={p.name:sha(p) for p in [np,evidence,out/'samples.json',out/'native.log']},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',STEM,result['seconds'],'minorbits',result['exact36_minor_numerator_bits'],result['exact36_minor_denominator_bits'])
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.probe.resolve())
