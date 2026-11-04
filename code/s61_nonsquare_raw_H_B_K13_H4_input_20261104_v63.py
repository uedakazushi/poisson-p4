#!/usr/bin/env python3
"""ONE H4-only prep; source COST only for m2/m3 and prepared B21."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63';sys.set_int_max_str_digits(0)
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
 z=[{0:1}]
 for i in range(n):z.append(mul(z[-1],p))
 return z
def integer_sum(polys):
 vals=[Q(v)for p in polys for m,v in p];D=math.lcm(*(v.denominator for v in vals));return D,sum(abs(int(v*D))for v in vals)
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 eb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_E_RECONSTRUCTION_20261104_v63';base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_LENGTH_K13_PREFLIGHT_20261104_v63';hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63';bb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63';assert sha(eb/'report.json')=='27a59e1018afe49f266285bfaac9662b5ab1e5f76718226fb9aa90c55b2bc2b5'
 e=read(eb/'native-certificate.json');g={a:int(v)for a,v in read(base/'K13-source-input.json')['g13']};tn=ip(e['e_numerator']);td=int(e['e_denominator']);H=read(hb/'certificate.json')['monic_H_laurent'];DH,_=integer_sum([H]);assert len(H)==74
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for p in[g,tn]:
   f.write(str(len(p))+'\n')
   for a,v in sorted(p.items(),reverse=True):f.write(str(a)+' '+str(v)+'\n')
  f.write(str(td)+'\n'+str(DH)+'\n74\n')
  for m,v in H:
   q=Q(v)*DH;assert q.denominator==1;f.write(str(-m[0])+' '+str(m[1])+' '+str(m[2])+' '+str(q.numerator)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_H4_input_20261104_v63.cpp';exe=work/'H4-input'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','H4-input',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 n=read(out/'native-certificate.json');ain={i-1:-v for i,v in g.items()if i};ap=powers(ain,5);tp=powers(tn,8);common=DH*g[0]**5*td**8;assert common==int(n['common_source_denominator']);raw=[{}for _ in range(5)]
 for m,v in H:raw[m[2]]=add(raw[m[2]],mul(ap[-m[0]],tp[m[1]]),int(Q(v)*DH)*g[0]**(5+m[0])*td**(8-m[1]))
 fields=[]
 for i,r in enumerate(n['H4_coefficients']):
  assert raw[i]==ip(r['raw_common_numerator']);R=ip(r['remainder']);assert {a:v*g[13]**r['L_power']for a,v in raw[i].items()}==add(mul(g,ip(r['g13_quotient'])),R)
  N=ip(r['numerator']);D=int(r['denominator']);assert {a:v*common*g[13]**r['L_power']for a,v in N.items()}=={a:v*D for a,v in R.items()};fields.append([[a,str(Q(v,D))]for a,v in sorted(N.items(),reverse=True)])
 assert fields[4]==[[0,'1']]
 actual=dict(H4_QQ_coefficient_terms=sum(len(p)for p in fields),coefficient_num_bits=max(max(abs(ip(r['numerator'])[a]).bit_length()for a in ip(r['numerator']))for r in n['H4_coefficients']),coefficient_den_bits=max(int(r['denominator']).bit_length()for r in n['H4_coefficients']),source_raw_terms=sum(len(p)for p in raw),full_22_NF_or_matrix_or_Gaussian_NOT_computed=True)
 (out/'H4-monic-field-coefficients.json').write_text(json.dumps(dict(g13=[[a,str(v)]for a,v in sorted(g.items(),reverse=True)],e_numerator=e['e_numerator'],e_denominator=e['e_denominator'],H4_coefficients_d0_to4=fields,actual=actual),indent=2)+'\n')
 # COST ONLY: no scalar mixed-polynomial evaluation.
 pencil=read(pb/'native-certificate.json');TM=max(sum(abs(v)for v in tn.values()),td);growth=abs(g[13])+sum(abs(v)for v in g.values());scalar_cost=[]
 for j in[1,2]:
  p=pencil['mixed3_quotients_by_final'][j];D,S=integer_sum([p]);ea=max(m[0]for m,v in p);ee=max(m[1]for m,v in p);deg=ea+12*ee;red=max(0,deg-12);rawN=S*TM**ee;NN=rawN*growth**red;DD=D*td**ee*abs(g[13])**red
  bound=lambda count,B:count*(150+(B*30103+99999)//100000)
  output=13*(150+((NN.bit_length()+DD.bit_length())*30103+99999)//100000);trace=bound(deg+1,rawN.bit_length())+bound(max(0,deg-12),NN.bit_length())+bound(13,NN.bit_length())+output+10000
  scalar_cost.append(dict(mixed_index=j+1,source_terms=len(p),source_a_e_degree=[ea,ee],raw_a_degree_upper=deg,g13_reduction_steps_upper=red,output_QQ_slots13=13,output_num_bits_upper=NN.bit_length(),output_den_bits_upper=DD.bit_length(),output_RAW_bytes_upper=output,whole_raw_reduction_certificate_bytes_upper=trace,no_evaluation_performed=True))
 prepared=read(bb/'certificate.json')['prepared_integer7x3'];DB,SB=integer_sum([p for row in prepared for p in row]);ba=max(m[0]for row in prepared for p in row for m,v in p);be=max(m[1]for row in prepared for p in row for m,v in p);bd=max(m[2]for row in prepared for p in row for m,v in p);Bdeg=ba+12*be;Bred=max(0,Bdeg-12);NB=SB*TM**be*growth**Bred;BD=DB*td**be*abs(g[13])**Bred;Bbytes=21*4*13*(150+((NB.bit_length()+BD.bit_length())*30103+99999)//100000)
 costs=dict(m2_m3=scalar_cost,two_scalar_whole_trace_RAW_upper=sum(c['whole_raw_reduction_certificate_bytes_upper']for c in scalar_cost),two_scalar_within32MiB=sum(c['whole_raw_reduction_certificate_bytes_upper']for c in scalar_cost)<32*1024*1024,B21_source_terms=sum(len(p)for row in prepared for p in row),B21_a_e_d_degrees=[ba,be,bd],B21_d_degree3_already_below_H4=True,B21_output_scalar_slots_upper=21*4*13,B21_num_bits_upper=NB.bit_length(),B21_den_bits_upper=BD.bit_length(),B21_output_RAW_bytes_upper=Bbytes,B21_entries_NF_NOT_computed=True,B21_localized_source_scaling_gate='MINORS_COSTpreparedinteger7x3 hasacceptedrow/columnknownFULLUunitandQQscalings; retainFULLU.',no_new_job_authorized_by_this_metadata=True)
 (out/'next-scalar-and-B21-preflight.json').write_text(json.dumps(costs,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',eb/'report.json',eb/'native-certificate.json',base/'K13-source-input.json',hb/'report.json',hb/'certificate.json',pb/'report.json',pb/'native-certificate.json',bb/'report.json',bb/'certificate.json']
 report=dict(id=ID,all_passed=True,H4_only_all5_coefficients_source_bound_exact=True,all_source_substitution_reductions_replayed_integer=True,monic4=True,actual=actual,no_m2_m3_evaluation_or_B21_NF_or_fullmatrix_Gaussian=True,native_seconds=n['native_seconds'],input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(json.dumps(costs,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output,a.work)
 except Exception as e:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'all_passed':False,'error':str(e),'native_certificate_preserved':(a.output/'native-certificate.json').exists(),'no_auto_retry':True},indent=2)+'\n')
  raise
