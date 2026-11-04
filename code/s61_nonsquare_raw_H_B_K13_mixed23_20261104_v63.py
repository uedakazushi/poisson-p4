#!/usr/bin/env python3
"""ONE m2/m3 evaluation, optional saved B21 coefficients iff bothzero."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_MIXED23_SCALARS_20261104_v63';sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(a):int(v)for a,v in p}
def add(a,b,c=1):
 z=dict(a)
 for i,v in b.items():
  z[i]=z.get(i,0)+c*v
  if not z[i]:del z[i]
 return z
def mul(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():z[i+j]=z.get(i+j,0)+c*d
 return{k:v for k,v in z.items()if v}
def powers(p,n):
 out=[{0:1}]
 for j in range(n):out.append(mul(out[-1],p))
 return out
def exact_division(p,g):
 z=dict(p);q={};dg=max(g);lc=g[dg]
 while z and max(z)>=dg:
  k=max(z)-dg;v=z[max(z)];assert v%lc==0;c=v//lc;q[k]=c;z=add(z,{a+k:c*v for a,v in g.items()},-1)
 assert not z;return q
def common(ps):return math.lcm(*(Q(v).denominator for p in ps for m,v in p))
def raw_eval(p,tp,td,ee,D,ddegree=None):
 ans={}
 for m,v in p:
  if ddegree is not None and m[2]!=ddegree:continue
  assert m[0]>=0 and m[1]<=ee and(m[2]==0 or ddegree is not None)
  ans=add(ans,{a+m[0]:c for a,c in tp[m[1]].items()},int(Q(v)*D)*td**(ee-m[1]))
 return ans
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63';eb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_E_RECONSTRUCTION_20261104_v63';pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63';bb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63';assert sha(hb/'report.json')=='3d54591f9a13846089d3dfa8f1afba3665f41c2078b794893f323079768d24cd'
 f=read(hb/'H4-monic-field-coefficients.json');g={a:int(v)for a,v in f['g13']};tn=ip(f['e_numerator']);td=int(f['e_denominator']);mixed=read(pb/'native-certificate.json')['mixed3_quotients_by_final'][1:3];B=read(bb/'certificate.json')['prepared_integer7x3'];DB=common([p for row in B for p in row]);Ds=[common([p])for p in mixed];ees=[max(m[1]for m,v in p)for p in mixed]
 inp=out/'native-input.txt'
 with inp.open('w')as s:
  for p in[g,tn]:
   s.write(str(len(p))+'\n')
   for a,v in sorted(p.items(),reverse=True):s.write(str(a)+' '+str(v)+'\n')
  s.write(str(td)+'\n')
  for p,D,ee in zip(mixed,Ds,ees):
   s.write(str(len(p))+' '+str(ee)+' '+str(D)+'\n')
   for m,v in p:s.write(str(m[0])+' '+str(m[1])+' '+str(int(Q(v)*D))+'\n')
  s.write(str(DB)+'\n')
  for row in B:
   for p in row:
    s.write(str(len(p))+'\n')
    for m,v in p:s.write(' '.join(map(str,m))+' '+str(int(Q(v)*DB))+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_mixed23_20261104_v63.cpp';exe=work/'mixed23'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','mixed23',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 n=read(out/'native-certificate.json');tp=powers(tn,23);scalar_stats=[]
 for row in n['scalars']:
  j=row['mixed_index']-2;ee=ees[j];raw=raw_eval(mixed[j],tp,td,ee,Ds[j]);assert raw==ip(row['raw_numerator']);cd=Ds[j]*td**ee;assert cd==int(row['common_source_denominator']);r=ip(row['remainder']);lp=row['L_power'];assert {a:v*g[13]**lp for a,v in raw.items()}==add(mul(g,ip(row['g13_quotient'])),r)
  num=ip(row['numerator']);den=int(row['denominator']);assert {a:v*cd*g[13]**lp for a,v in num.items()}=={a:v*den for a,v in r.items()}
  if row['unit_modg13']:assert add(mul(ip(row['unit_S']),ip(row['primitive_numerator'])),mul(ip(row['unit_T']),g))=={0:int(row['unit_constant'])}
  scalar_stats.append(dict(mixed_index=j+2,zero=not num,terms=len(num),unit=row['unit_modg13']))
 B_stats=None
 if n['B21_coefficients_prepared']:
  bc=read(out/'native-certificate.json.B21.json');assert len(bc['entries'])==21;bitsn=bitsd=terms=nonzero=0;cd=DB*td**4
  for entry,p in enumerate([p for row in B for p in row]):
   for d,r in enumerate(bc['entries'][entry]):
    raw=raw_eval(p,tp,td,4,DB,d);num=ip(r['numerator']);den=int(r['denominator']);diff=add({a:v*den for a,v in raw.items()},{a:v*cd for a,v in num.items()},-1);exact_division(diff,g);terms+=len(num);nonzero+=bool(num);bitsn=max(bitsn,max((abs(v).bit_length()for v in num.values()),default=0));bitsd=max(bitsd,den.bit_length())
  B_stats=dict(Ka13_coefficient_terms=terms,nonzero_field_coefficients=nonzero,max_numerator_bits=bitsn,max_denominator_bits=bitsd,all84_coefficient_source_images_independent_integer_checked=True,no_H4reduction_matrix_Gaussian_or_rightinverse=True)
 classification='additional_scalar_unit_old_DG_bad_EMPTY_global_root_gate'if n['any_scalar_unit']else'both_mixed23_zero_candidate_retained_B21_input_only'if n['both_scalars_zero']else'proper_scalar_gcd_candidate_retained'
 proof=dict(classification=classification,scalar_stats=scalar_stats,all_raw_powers_pseudodivisions_canonical_outputs_replayed=True,B21_actual=B_stats,source_membership_scope='Mi=G*mi, rawMi frompenciladjugate inactualminorideal afterknownFULLUunits. DivideG ONLYonD(G); G0separate.',no_global_Bbad_empty_or_existence_claim=True,no_e_Euclid_newnorm_matrix_Gaussian=True)
 (out/'source-product-checks.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',hb/'report.json',hb/'next-scalar-and-B21-preflight.json',hb/'H4-monic-field-coefficients.json',eb/'report.json',pb/'report.json',pb/'native-certificate.json',bb/'report.json',bb/'certificate.json']
 report=dict(id=ID,all_passed=True,classification=classification,scalar_stats=scalar_stats,B21_prepared=n['B21_coefficients_prepared'],B21_actual=B_stats,actual_old_DG_empty=bool(n['any_scalar_unit']),no_global_Bbad_claim=True,no_H4reduction_eEuclid_newnorm_matrix_Gaussian_point_prime=True,native_seconds=n['native_seconds'],input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output,a.work)
 except Exception as e:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'all_passed':False,'error':str(e),'completed_native_certificate_preserved':(a.output/'native-certificate.json').exists(),'no_auto_retry':True},indent=2)+'\n')
  raise
