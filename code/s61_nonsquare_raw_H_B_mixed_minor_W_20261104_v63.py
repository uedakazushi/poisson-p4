#!/usr/bin/env python3
"""ONE source-bound saved W product/gcd, no field inverse/e-Euclid/new norm."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_MIXED_MINOR_W_20261104_v63';sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return {int(a):Q(v)for a,v in p}
def integers(p):return {int(a):int(v)for a,v in p}
def mul(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():z[i+j]=z.get(i+j,0)+c*d
 return{k:v for k,v in z.items()if v}
def add(a,b,c=1):
 z=dict(a)
 for i,v in b.items():
  z[i]=z.get(i,0)+c*v
  if not z[i]:del z[i]
 return z
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_LENGTH_K13_PREFLIGHT_20261104_v63';assert sha(base/'report.json')=='5f87e4e632fa9975a6936bca3f34cf01dfb4912f56569e3a8b29abbe2c3a9555';pref=read(base/'mixed-minor-W-source-preflight.json');source=read(base/'K13-source-input.json')
 g=poly(source['g13']);assert all(v.denominator==1 for v in g.values());g={a:int(v)for a,v in g.items()};A=[[poly(p)for p in row]for row in pref['Nm1_modg13']];V=[poly(p)for p in pref['normalized_V_modg13']]
 DA=math.lcm(*(v.denominator for row in A for p in row for v in p.values()));DV=math.lcm(*(v.denominator for p in V for v in p.values()));N=[[{a:int(v*DA)for a,v in p.items()}for p in row]for row in A];vi=[{a:int(v*DV)for a,v in p.items()}for p in V];assert [DA.bit_length(),DV.bit_length()]==pref['common_numeric_denominator_bits']
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for p in[g]+[p for row in N for p in row]+vi:
   f.write(str(len(p))+'\n')
   for a,v in sorted(p.items(),reverse=True):f.write(str(a)+' '+str(v)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_mixed_minor_W_20261104_v63.cpp';exe=work/'mixed-minor-W'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','W-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 c=read(out/'native-certificate.json');w=read(out/'native-certificate.json.W-checkpoint.json');assert len(w['rows'])==16;remainders=[];lc=g[13]
 for i,row in enumerate(w['rows']):
  raw={}
  for j in range(16):raw=add(raw,mul(N[i][j],vi[j]))
  assert raw==integers(row['raw_integer_product']);qq=integers(row['g13_quotient']);rr=integers(row['remainder']);power=row['leading_coefficient_power'];assert 0<=power<=12
  assert {a:v*lc**power for a,v in raw.items()}==add(mul(g,qq),rr);assert max(rr,default=-1)<13;remainders.append(rr)
 previous=g
 for step in c['gcd_chain']:
  i=step['component'];new=integers(step['new_gcd']);assert previous==integers(step['previous_gcd']);assert mul(new,integers(step['previous_quotient']))==previous;assert mul(new,integers(step['component_quotient']))==remainders[i]
  lhs=add(mul(integers(step['Bezout_S']),previous),mul(integers(step['Bezout_T']),remainders[i]));assert lhs=={a:int(step['Bezout_constant'])*v for a,v in new.items()};previous=new
 assert previous==integers(c['final_gcd'])
 if c['gcd_is_QQunit']:classification='triple_T_ZERO_old_DG_bad_EMPTY_global_combine_root_only'
 elif c['all_W_components_zero']:classification='all16_W_zero_proper_g13_candidate_retained_no_next_method'
 else:classification='proper_gcd_candidate_retained_no_next_method'
 proof=dict(classification=classification,actual_W_modg13='remainder_i/(DA*DV*leading(g13)^power_i)',DA=str(DA),DV=str(DV),V_common_content=pref['V_common_QQ_content'],V_common_known_a_power=pref['V_common_known_a_power'],normalization_only_ONE_common_QQ_a_unit=True,all16_raw_products_and_pseudoreduction_identities_independent_integer_replayed=True,all_gcd_division_and_Bezout_identities_replayed=True,no_field_inverse=True,unit_criterion='W_i isQQ/aunit timesmixedmaximalminorof[Nk|Nm1]. Fitt0(T) containsg13andallWi. Unitgcd=>T=0.',actual_bad_scope='OnlyD(G); G0andwholeglobalcombinationareparentsseparategates.',all_zero_not_nonempty_proof=True,no_e_Euclid_newnorm_or_method_switch=True)
 (out/'proof-and-product-checks.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',base/'report.json',base/'mixed-minor-W-source-preflight.json',base/'K13-source-input.json',base/'fitting-length-lemma.json']
 report=dict(id=ID,all_passed=True,classification=classification,all_W_components_zero=c['all_W_components_zero'],nonzero_W_components=c['nonzero_W_components'],gcd_degree=max(integers(c['final_gcd'])),old_DG_empty_proved=bool(c['gcd_is_QQunit']),no_global_Bbad_empty_claim=True,no_field_inverse_e_Euclid_newnorm_point_prime_GB=True,native_seconds=c['native_seconds'],input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output,a.work)
 except Exception as e:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'all_passed':False,'error':str(e),'partial_checkpoints_preserved':[p.name for p in a.output.glob('*checkpoint.json')],'no_auto_retry':True},indent=2)+'\n')
  raise
