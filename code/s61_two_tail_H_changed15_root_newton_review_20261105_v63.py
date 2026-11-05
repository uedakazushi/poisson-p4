"""Independent fixed changed15 postcheck: Horner plus four Newton/Faddeev steps."""
from pathlib import Path
import json,hashlib,resource,signal,time
R=Path(__file__).resolve().parents[1];P=32003
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=R/'reports/S61_TWO_TAIL_H_CHANGED15_ROOT_NEWTON_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);start=time.monotonic()
pins={};products=0;bits=0;state={'stage':'SOURCE'}
def load(s,h):
 p=R/s;assert sha(p)==h;pins[s]=h;return json.loads(p.read_text())
c=load('reports/S61_TWO_TAIL_H_CHANGED15_CUSTODY_20261105_v63/report.json','38216fcd90439f1efcdbf2acbc9e236f63092dce7c01e49c8c767b7f6e833755')
p=load('reports/S61_TWO_TAIL_H_CHANGED15_AXIS_ENTRY_20261105_v63/report.json','2b93dce6522f3f08300270b278bbb860f0cb29e7cbfcc31c9a777ed4bed275ee')
review=load('reports/S61_TWO_TAIL_H_CHANGED15_ROOT_SOURCE_REVIEW_20261105_v63/report.json','95d94bfb53418cc6094e8d2bdacdb704496cacf627d98f7389830007306cae74')
q=load('reports/S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63/report.json','041d9234ccf704286eacd0da6e0f6161d8b268b27dfaa678e2080657f5f32193')
(out/'phase.json').write_text(json.dumps(dict(SOURCE_prior=True,only_same_axis_prime_quotient=True,no_new_physical_source=True,input_sha256=pins,field_temp_bits=30,existing_ROM_bits=188,hard_products=10000,owned_upper=4096,RAW_cap=4194304,hard_seconds=60,independent_method='Horner scalar polys and Newton4; no subsetDP or producer arithmetic functions'))+'\n')
H=[17430,11013,30453,27809]
def mul(x,y):
 global products,bits
 z=x*y;products+=1;bits=max(bits,z.bit_length());assert products<=10000 and bits<=30
 return z%P
def inv(x):
 n=P-2;a=x%P;z=1;assert a
 while n:
  if n&1:z=mul(z,a)
  a=mul(a,a);n>>=1
 assert mul(x%P,z)==1;return z
def add(a,b):return [(x+y)%P for x,y in zip(a,b)]
def rm(a,b):
 z=[0]*7
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]=(z[i+j]+mul(x,y))%P
 for k in [6,5,4]:
  for j,h in enumerate(H):z[k-4+j]=(z[k-4+j]-mul(z[k],h))%P
 return z[:4]
def by_a(v):return [(-mul(v[3],H[0]))%P]+[(v[i-1]-mul(v[3],H[i]))%P for i in range(1,4)]
try:
 M=[];state['stage']='Horner16'
 for cell in c['axis_input_polynomials']:
  co={a:int(v)%P for a,v in cell['coefficients']};v=[0]*4
  for a in range(9,-1,-1):v=by_a(v);v[0]=(v[0]+co.get(a,0))%P
  den=inv(int(cell['D'])%P);M.append([mul(x,den)for x in v])
 state['matrix']=M;B=[[int(i==j),0,0,0]for i in range(4)for j in range(4)];stages=[]
 for k in range(1,5):
  state.update(stage='Newton4',k=k,B=B)
  C=[[0]*4 for _ in range(16)]
  for i in range(4):
   for j in range(4):
    for a in range(4):C[4*i+j]=add(C[4*i+j],rm(M[4*i+a],B[4*a+j]))
  trace=[0]*4
  for i in range(4):trace=add(trace,C[5*i])
  ik=inv(k);ck=[(-mul(x,ik))%P for x in trace]
  B=[v[:]for v in C]
  for i in range(4):B[5*i]=add(B[5*i],ck)
  stages.append(dict(k=k,trace=trace,c=ck,B=B))
 assert all(v==[0]*4 for v in B)
 det=stages[-1]['c'];combined=rm(q['D12_remainder'],det)
 assert det==p['detPair_remainder'] and combined==p['combined_remainder']==[0]*4
 evidence=dict(matrix16=M,Newton_stages4=stages,CayleyHamilton_B4_zero=True,detPair=det,combined=combined)
 ep=out/'independent-Horner-Newton4-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 r=dict(id=out.name,all_passed=True,independent_same_fixed_quotient_review=True,input_sha256=pins,detPair_remainder=det,combined_remainder=combined,
 zero_combined_inconclusive=True,no_global_H_divisibility_or_rank_bad_inference=True,all_four_Newton_stages_and_CayleyHamilton_checked=True,
 field_products=products,max_field_temp_bits=bits,owned_upper=4096,existing_ROM_input_bits=188,code_sha256=sha(Path(__file__)),evidence_sha256={ep.name:sha(ep),'phase.json':sha(out/'phase.json')},seconds=time.monotonic()-start)
 rp=out/'report.json';rp.write_text(json.dumps(r,indent=2)+'\n');print(sha(rp),products,det,combined)
except BaseException as e:
 (out/'failure-prefix.json').write_text(json.dumps(dict(error=repr(e),state=state,no_retry=True),indent=2)+'\n');raise
