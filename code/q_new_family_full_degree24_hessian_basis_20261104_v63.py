#!/usr/bin/env python3
"""ONE corrected same-plane Hessian16 coefficient basis job; no source expansion."""
import argparse,gzip,itertools,json,resource,shlex,subprocess,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from weight1_structure import sp_basis
from q_new_family_full_degree24_universal_basis_20261104_v63 import read,sha,echelon,det_mod,necklaces
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63'
P=32003
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def mod(z):
 q=Fraction(str(z));assert q.denominator%P;return q.numerator*pow(q.denominator,-1,P)%P
def main(out,work):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 faildir=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63'
 failpath=faildir/'report.json';fail=read(failpath);assert fail['all_passed']and fail['status']=='degree12_rank_gate_failed'
 for n,h in fail['evidence_sha256'].items():assert sha(faildir/n)==h
 cdpath=faildir/'native-certificate.json.gz';cd=read(cdpath)
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 ambientpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
 old12path=ROOT/'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
 oldcode=ROOT/'code/q_new_family_degree12_14_invariant_exclusion_20261103_v63.py'
 prepath=ROOT/'reports/Q_NEW_FAMILY_DEGREE24_THETA_H_REPAIR_PREFLIGHT_20261104_v63/report.json'
 charpath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
 info,ambient,old12,pre,char=map(read,[invpath,ambientpath,old12path,prepath,charpath])
 assert old12['full_degree12_dimension']==9 and pre['preflight_only_no_new_invariant_or_rank_computation']
 dimensions={int(r['degree']):int(r['invariant_dimension'])for r in char['character_records_through24']};assert [dimensions[i]for i in[12,16,24]]==[9,20,67]
 O=s.Matrix(info['symplectic_form16']);Oi=O.inv();rep=[s.Matrix(r)for r in info['all_ten_representation_matrices']]
 gs=sp_basis();Gi=s.Matrix([[s.trace(a*b)for b in gs]for a in gs]).inv();vectors=ambient['ambient_beta16'][:3]
 oldpre=read(ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63/report.json');assert vectors==oldpre['fixed_ambient_vectors_from_saved_probe']
 inp=work/'formal-input.txt'
 with inp.open('w')as f:
  f.write('16 32003\n')
  for v in vectors:f.write(' '.join(map(str,v))+'\n')
  for mat in[O,Oi,Gi,*rep]:f.write(' '.join(str(mod(z))for z in mat)+'\n')
  f.write(' '.join(map(str,cd['F4']))+'\n')
  for row in cd['E5']+cd['degree12_columns']+cd['trace8_columns']:f.write(' '.join(map(str,row))+'\n')
 with gzip.open(out/'formal-input.txt.gz','wb')as f:f.write(inp.read_bytes())
 native=ROOT/'code/q_new_family_full_degree24_hessian_basis_20261104_v63.cpp';header=ROOT/'code/q_new_family_full_degree24_universal_basis_20261104_v63.cpp'
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));binary=work/'hessian-basis-native'
 subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 raw=work/'native-certificate.json';error=None;data=None
 try:
  with(out/'native.log').open('w')as log:subprocess.run([str(binary),str(inp),str(raw)],stdout=log,stderr=subprocess.STDOUT,timeout=60,preexec_fn=cap,check=True)
  assert raw.stat().st_size<=32*1024*1024;data=read(raw);assert data['all_passed']and data['H_is16x16_Hessian_not4x4_D']and data['old_H_symplectic_matrix_identity_verified']
 except Exception as e:error=repr(e)
 witnesses={};definitions={};complete=False;new_source=[]
 words6,words8=necklaces(6),necklaces(8)
 if data:
  assert data['degree12_columns'][:20]==cd['degree12_columns']
  assert data['old_trace6_columns']==data['degree12_columns'][20:]
  assert len(data['old_trace6_columns'])==14
  for degree in[12,16,24]:
   columns=data[f'degree{degree}_columns']
   if not columns:continue
   assert all(len(row)==(degree+1)*(degree+2)//2 and all(isinstance(z,int)and 0<=z<P for z in row)for row in columns)
   cols,rows=echelon(columns);assert cols==data[f'degree{degree}_selected_columns']and rows==data[f'degree{degree}_selected_coefficient_rows']
   assert len(cols)==data[f'degree{degree}_rank']<=dimensions[degree]
   square=[[columns[j][i]for j in cols]for i in rows];det=det_mod(square);assert det
   witnesses[str(degree)]=dict(rank=len(cols),coefficient_rows=rows,columns=cols,minor_mod32003=det)
  cddefs=read(faildir/'definitions.json');defs12=cddefs['degree12']+[dict(kind='theta_H_trace6',word=list(w))for w in words6]
  definitions['degree12']=defs12;definitions['degree12_basis']=[defs12[i]for i in data['degree12_selected_columns']]
  if data['degree16_columns']:
   assert data['degree16_columns'][24:60]==cd['trace8_columns']
   assert len(echelon(data['degree16_columns'][:24])[0])==data['degree16_product_rank']
   if data['old_trace8_attempted']:assert len(data['old_trace8_columns'])==36 and data['degree16_columns'][60:]==data['old_trace8_columns']
  complete=data['status']=='complete_Inv24_basis_certified'
  if complete:
   assert [data[f'degree{d}_rank']for d in[12,16,24]]==[9,20,67]and data['degree16_product_rank']==19
   selected=data['degree12_selected_columns'];assert selected[:8]==[0,1,2,3,4,5,9,10]and selected[-1]>=20
   prim12=[defs12[i]for i in selected[5:]];assert prim12[0]['kind']=='N12'
   definitions['primitive12']=prim12
   index=data['primitive16_CD_trace8_index'];oldindex=data['primitive16_old_trace8_index']
   assert (index>=0)!=(oldindex>=0)
   prim16=dict(kind='trace8',word=list(words8[index]))if index>=0 else dict(kind='theta_H_trace8',word=list(words8[oldindex]))
   definitions['primitive16']=prim16
   spectral=[]
   for r in range(4):
    for ij in itertools.combinations_with_replacement(range(4),r):
     powers=[6-2*r,0,0,0,0]
     for j in ij:powers[1+j]+=1
     spectral.append(dict(kind='spectral',F4_p1_p2_p3_p4_exponents=powers))
   definitions['degree24']=spectral+[dict(kind='primitive12_times_spectral12',primitive=i,spectral=j)for i in range(4)for j in range(5)]+[dict(kind='primitive12_product',indices=list(ij))for ij in itertools.combinations_with_replacement(range(4),2)]+[dict(kind='primitive16_times_Inv8',index=i)for i in range(5)]
   assert len(definitions['degree24'])==70
   definitions['degree24_basis']=[definitions['degree24'][i]for i in data['degree24_selected_columns']]
   new_source=prim12[1:]+[prim16]
  with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
 (out/'definitions.json').write_text(json.dumps(dict(Inv8_order=['F4^2','p1','p2','p3','p4'],**definitions),indent=2)+'\n')
 (out/'minor-witnesses.json').write_text(json.dumps(witnesses,indent=2)+'\n')
 inputs=[Path(__file__),native,header,failpath,cdpath,invpath,ambientpath,old12path,oldcode,prepath,charpath,ROOT/'code/weight1_structure.py']
 report=dict(id=STEM,all_passed=error is None,status=data['status']if data else 'bounded_failure',universal_complete_QQ_Inv24_basis=complete,
  character_dimensions={str(k):dimensions[k]for k in[12,16,24]},ranks={k:v['rank']for k,v in witnesses.items()},degree16_product_rank=data['degree16_product_rank']if data else None,
  structural_repair='old theta/H act on W16, H is sp16 Hessian; CD-only packet preserved',
  formal_plane='SAME first3 saved ambient beta vectors: beta0+s*beta1+t*beta2',prime=P,
  one_corrected_native_job_only=True,no_actual_NEW25_point_evaluations=True,no_source_invariant_expansion=True,
  old_trace8_attempted=data['old_trace8_attempted']if data else None,
  QQ_independence_reason='nonzero mod32003 polynomial-coefficient minors, p-unit constant matrices and character upper bounds',
  selected_new_future_source_words=new_source,no_future_source_job_authorized_or_started=True,
  limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,certificate_bytes=32*1024*1024),error=error,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','status','universal_complete_QQ_Inv24_basis','ranks','degree16_product_rank','old_trace8_attempted','selected_new_future_source_words','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();main(args.output.resolve(),args.work.resolve())
