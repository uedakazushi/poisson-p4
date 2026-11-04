#!/usr/bin/env python3
"""ONE formal ambient two-plane coefficient test, with explicit rank failure gates."""
import argparse,gzip,hashlib,itertools,json,resource,shlex,subprocess,time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63'
P=32003
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def necklaces(n):return [w for w in itertools.product(range(2),repeat=n)if w==min(w[i:]+w[:i]for i in range(n))]
def echelon(columns):
 basis={};chosen=[];rows=[]
 for k,c in enumerate(columns):
  r=c[:]
  for i,b in sorted(basis.items()):
   q=r[i]
   if q:r=[(a-q*z)%P for a,z in zip(r,b)]
  i=next((j for j,z in enumerate(r)if z),None)
  if i is None:continue
  inv=pow(r[i],-1,P);basis[i]=[z*inv%P for z in r];chosen.append(k);rows.append(i)
 return chosen,rows
def det_mod(A):
 A=[r[:]for r in A];det=1
 for i in range(len(A)):
  j=next(k for k in range(i,len(A))if A[k][i])
  if j!=i:A[i],A[j]=A[j],A[i];det=-det
  pivot=A[i][i];det=det*pivot%P;inv=pow(pivot,-1,P)
  for k in range(i+1,len(A)):
   q=A[k][i]*inv%P
   if q:A[k]=[(a-q*b)%P for a,b in zip(A[k],A[i])]
 return det%P
def main(out,work):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 prepath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63/report.json'
 ambientpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 charpath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
 pre,ambient,spin,char=map(read,[prepath,ambientpath,spinpath,charpath]);assert pre['all_passed']
 vectors=ambient['ambient_beta16'][:3];assert vectors==pre['fixed_ambient_vectors_from_saved_probe']
 dimensions={int(r['degree']):int(r['invariant_dimension'])for r in char['character_records_through24']}
 assert [dimensions[i]for i in[12,16,24]]==[9,20,67]
 mon=spin['quadratic_monomials'];assert len(mon)==136
 coeff=[]
 for key in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  for row in spin[key]:
   rr=[4*Fraction(str(c))for c in row];assert all(x.denominator==1 for x in rr);coeff.append([int(x)for x in rr])
 inp=work/'formal-input.txt'
 with inp.open('w')as f:
  f.write('16 136 32003\n')
  for row in vectors:f.write(' '.join(map(str,row))+'\n')
  for pair in mon:f.write(' '.join(map(str,pair))+'\n')
  for row in coeff:f.write(' '.join(map(str,row))+'\n')
 with gzip.open(out/'formal-input.txt.gz','wb')as f:f.write(inp.read_bytes())
 native=ROOT/'code/q_new_family_full_degree24_universal_basis_20261104_v63.cpp';binary=work/'formal-basis-native'
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
 subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 raw=work/'native-certificate.json';error=None;data=None
 try:
  with(out/'native.log').open('w')as log:subprocess.run([str(binary),str(inp),str(raw)],stdout=log,stderr=subprocess.STDOUT,timeout=60,preexec_fn=cap,check=True)
  assert raw.stat().st_size<=32*1024*1024;data=read(raw);assert data['all_passed']
 except Exception as e:error=repr(e)
 witnesses={};defs={};complete=False
 if data:
  assert data['prime']==P and data['max_intermediate_terms']<=325
  for degree in[12,16,24]:
   columns=data[f'degree{degree}_columns']
   if not columns:
    assert data[f'degree{degree}_rank']==0 and not data[f'degree{degree}_selected_columns']
    continue
   for c in columns:assert len(c)==(degree+1)*(degree+2)//2 and all(isinstance(x,int)and 0<=x<P for x in c)
   cols,rows=echelon(columns)
   assert cols==data[f'degree{degree}_selected_columns']and rows==data[f'degree{degree}_selected_coefficient_rows']
   assert len(cols)==data[f'degree{degree}_rank']<=dimensions[degree]
   matrix=[[columns[j][i]for j in cols]for i in rows];det=det_mod(matrix)
   assert det
   witnesses[str(degree)]=dict(rank=len(cols),coefficient_rows=rows,columns=cols,minor_mod32003=det)
  complete=data['status']=='complete_Inv24_basis_certified'
  words6,words8=necklaces(6),necklaces(8)
  defs['degree12']=[dict(kind='F4_times_Inv8',index=i)for i in range(5)]+[dict(kind='N12')]+[dict(kind='trace6',word=list(w))for w in words6]
  selected12=data['degree12_selected_columns']
  defs['degree12_basis']=[defs['degree12'][i]for i in selected12]
  if data['degree16_columns']:
   assert len(data['degree16_columns'])==60
   rankprod=len(echelon(data['degree16_columns'][:24])[0]);assert rankprod==data['degree16_product_rank']
   defs['degree16']=[dict(kind='F4_times_basis12',index=i)for i in range(9)]+[dict(kind='Inv8_product',indices=list(ij))for ij in itertools.combinations_with_replacement(range(5),2)]+[dict(kind='trace8',word=list(w))for w in words8]
  if complete:
   assert data['degree12_lower_rank']==5 and selected12[5]==5 and len(selected12)==9
   assert data['degree16_product_rank']==19 and data['degree16_rank']==20 and data['degree24_rank']==67
   defs['primitive12']=[defs['degree12'][i]for i in selected12[5:]]
   index=data['primitive16_trace8_index'];assert 0<=index<36;defs['primitive16']=dict(kind='trace8',word=list(words8[index]))
   spectral=[]
   for r in range(4):
    for ij in itertools.combinations_with_replacement(range(4),r):
     powers=[6-2*r,0,0,0,0]
     for j in ij:powers[1+j]+=1
     spectral.append(dict(kind='spectral',F4_p1_p2_p3_p4_exponents=powers))
   defs['degree24']=spectral+[dict(kind='primitive12_times_spectral12',primitive=i,spectral=j)for i in range(4)for j in range(5)]+[dict(kind='primitive12_product',indices=list(ij))for ij in itertools.combinations_with_replacement(range(4),2)]+[dict(kind='primitive16_times_Inv8',index=i)for i in range(5)]
   assert len(defs['degree24'])==70
   defs['degree24_basis']=[defs['degree24'][i]for i in data['degree24_selected_columns']]
  with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
 (out/'definitions.json').write_text(json.dumps(dict(Inv8_order=['F4^2','p1','p2','p3','p4'],**defs),indent=2)+'\n')
 (out/'minor-witnesses.json').write_text(json.dumps(witnesses,indent=2)+'\n')
 inputs=[Path(__file__),native,prepath,ambientpath,spinpath,charpath]
 report=dict(id=STEM,all_passed=error is None,status=data['status']if data else 'bounded_failure',
  universal_complete_QQ_Inv24_basis=complete,character_dimensions={str(k):dimensions[k]for k in[12,16,24]},
  ranks={k:v['rank']for k,v in witnesses.items()},degree16_product_rank=data['degree16_product_rank']if data else None,
  formal_plane='first three saved ambient beta vectors: beta0+s*beta1+t*beta2',prime=P,
  one_native_coefficient_job_only=True,no_actual_NEW25_point_evaluations=True,no_source_invariant_expansion=True,
  trace8_coefficient_columns_precomputed_but_not_rank_tested_on_degree12_failure=bool(data and data['status']=='degree12_rank_gate_failed'),
  degree16_and_degree24_rank_tests_attempted=bool(data and data['degree16_columns']),
  QQ_independence_reason='nonzero mod32003 polynomial-coefficient minor and 4-power unit denominators',
  rank_failure_stops_without_alternative_plane_or_candidate_growth=True,
  limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,certificate_bytes=32*1024*1024),
  error=error,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','status','universal_complete_QQ_Inv24_basis','ranks','degree16_product_rank','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
 args=p.parse_args();main(args.output.resolve(),args.work.resolve())
