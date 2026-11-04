#!/usr/bin/env python3
"""ONE universal Poisson-generated Inv26 basis job on the SAME formal plane."""
import argparse,gzip,itertools,json,resource,shlex,subprocess,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from q_new_family_degree16_invariant_probe_20261103_v63 import character
from q_new_family_full_degree24_universal_basis_20261104_v63 import read,sha,echelon,det_mod
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63';PRIME=32003
def mod(z):
 q=Fraction(str(z));assert q.denominator%PRIME;return q.numerator*pow(q.denominator,-1,PRIME)%PRIME
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def main(out,work):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 ambientpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
 prepath=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_PREFLIGHT_20261104_v63/report.json'
 basisdir=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63'
 gatepath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW/report.json'
 oldprepath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63/report.json'
 cdpath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63/native-certificate.json.gz'
 inv,spin,ambient,pre,gate,oldpre,oldcd=map(read,[invpath,spinpath,ambientpath,prepath,gatepath,oldprepath,cdpath]);assert pre['all_passed']and gate['all_passed']
 for n,h in gate['input_sha256'].items():assert sha(ROOT/n)==h,n
 vectors=ambient['ambient_beta16'][:3];assert vectors==oldpre['fixed_ambient_vectors_from_saved_probe']
 old=read(basisdir/'native-certificate.json.gz');assert old['status']=='complete_Inv24_basis_certified';defs12=read(basisdir/'definitions.json')['degree12_basis']
 Inv12=[old['degree12_columns'][i]for i in old['degree12_selected_columns']];assert len(Inv12)==9
 chars=character(list(map(tuple,inv['Cartan_weights'])),26);dims={r['degree']:r['invariant_dimension']for r in chars}
 assert[dims[d]for d in[8,12,14,18,24,25,26]]==[5,9,3,9,67,0,45]
 (out/'character-records.json').write_text(json.dumps(chars,indent=2)+'\n')
 Oi=s.Matrix(inv['symplectic_form16']).inv();coeff=[]
 for key in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  for row in spin[key]:
   v=[4*Fraction(str(z))for z in row];assert all(z.denominator==1 for z in v);coeff.append(list(map(int,v)))
 inp=work/'formal-input.txt'
 with inp.open('w')as f:
  f.write('16 136 32003\n')
  for r in vectors:f.write(' '.join(map(str,r))+'\n')
  for r in spin['quadratic_monomials']+coeff:f.write(' '.join(map(str,r))+'\n')
  f.write(' '.join(str(mod(z))for z in Oi)+'\n');f.write(' '.join(map(str,oldcd['F4']))+'\n')
  for r in Inv12:f.write(' '.join(map(str,r))+'\n')
 with gzip.open(out/'formal-input.txt.gz','wb')as f:f.write(inp.read_bytes())
 native=ROOT/'code/q_new_family_degree26_poisson_basis_20261104_v63.cpp';helper=ROOT/'code/q_new_family_full_degree24_universal_basis_20261104_v63.cpp';binary=work/'poisson26-native'
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 raw=work/'native-certificate.json';error=None;data=None
 try:
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(raw)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=60,preexec_fn=cap)
  assert raw.stat().st_size<=32*1024*1024;data=read(raw);assert data['all_passed']and data['ambient16_derivatives_before_substitution']and data['all32_quadratic_and_all8_trace_Euler_identities']
 except Exception as e:error=repr(e)
 witnesses={};definitions={};complete=False
 w4=['0000','0001','0011','0101','0111'];w6=['000101','000111','011011'];pairs=list(itertools.combinations(range(5),2))
 definitions['Inv8_trace_words']=w4;definitions['CD_trace12_words']=w6;definitions['full_Inv12_basis']=defs12
 definitions['degree14']=[dict(kind='Poisson_Inv8_Inv8',left=i,right=j)for i,j in pairs]
 if data:
  assert len(data['ambient16_gradient_columns'])==8
  for k,grad in enumerate(data['ambient16_gradient_columns']):assert len(grad)==16 and all(len(g)==(36 if k<5 else 78)for g in grad)
  assert len(echelon(data['Inv8_trace4_columns'])[0])==5
  for degree in[14,18,26]:
   columns=data[f'degree{degree}_columns']
   if not columns:continue
   assert all(len(c)==(degree+1)*(degree+2)//2 and all(isinstance(x,int)and 0<=x<PRIME for x in c)for c in columns)
   cols,rows=echelon(columns);assert cols==data[f'degree{degree}_selected_columns']and rows==data[f'degree{degree}_selected_coefficient_rows']
   assert len(cols)==data[f'degree{degree}_rank']<=dims[degree]
   det=det_mod([[columns[j][i]for j in cols]for i in rows]);assert det
   witnesses[str(degree)]=dict(rank=len(cols),coefficient_rows=rows,columns=cols,minor_mod32003=det)
  definitions['Inv14_basis']=[definitions['degree14'][i]for i in data['degree14_selected_columns']]
  if data['degree18_columns']:
   assert len(echelon(data['degree18_columns'][:3])[0])==data['degree18_F4_product_rank']
   definitions['degree18']=[dict(kind='F4_Inv14',Inv14_basis_index=i)for i in range(3)]+[dict(kind='Poisson_Inv8_CD_trace12',Inv8_index=i,CD_trace12_index=j)for i in range(5)for j in range(3)]
   definitions['Inv18_basis']=[definitions['degree18'][i]for i in data['degree18_selected_columns']]
  complete=data['status']=='complete_QQ_Inv26_basis_certified'
  if complete:
   assert[data[f'degree{d}_rank']for d in[14,18,26]]==[3,9,45]and data['degree18_F4_product_rank']==3
   p18=data['primitive18_selected_candidate_columns'];assert len(p18)==6 and all(i>=3 for i in p18)
   definitions['primitive18_basis']=[definitions['degree18'][i]for i in p18]
   definitions['degree26']=[dict(kind='Inv14_times_Inv12',Inv14_basis_index=i,Inv12_basis_index=j)for i in range(3)for j in range(9)]+[dict(kind='primitive18_times_Inv8',primitive18_basis_index=i,Inv8_index=j)for i in range(6)for j in range(5)]
   definitions['Inv26_basis']=[definitions['degree26'][i]for i in data['degree26_selected_columns']];assert len(definitions['Inv26_basis'])==45
  with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
 (out/'definitions.json').write_text(json.dumps(definitions,indent=2)+'\n');(out/'minor-witnesses.json').write_text(json.dumps(witnesses,indent=2)+'\n')
 inputs=[Path(__file__),native,helper,ROOT/'code/q_new_family_degree16_invariant_probe_20261103_v63.py',invpath,spinpath,ambientpath,prepath,gatepath,basisdir/'native-certificate.json.gz',basisdir/'definitions.json',oldprepath,cdpath]
 report=dict(id=STEM,all_passed=error is None,status=data['status']if data else 'bounded_failure',complete_QQ_Inv26_basis=complete,
  character_dimensions={str(d):dims[d]for d in[8,12,14,18,24,25,26]},ranks={k:v['rank']for k,v in witnesses.items()},
  ambient16_derivatives_of_quadratic_CD_before_restriction=True,Poisson_convention='gradient(f)^T Omega^-1 gradient(g)',
  SAME_first3_saved_ambient_vectors=True,one_native_formal_coefficient_job_only=True,no_H16_derivative_or_new_point_or_source26_job=True,
  no_automatic_candidate_growth_on_failure=True,NEW25_defining_degree_lower_bound_unchanged=26,NEW25_equation_or_degree26_restriction_claimed=False,
  limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,certificate_bytes=32*1024*1024),error=error,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','status','complete_QQ_Inv26_basis','ranks','seconds','error']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();main(args.output.resolve(),args.work.resolve())
