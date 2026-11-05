#!/usr/bin/env python3
"""Independent p4 source, corrected-zero3 witnesses and prepared fixed3x8 gate."""
import ast,hashlib,itertools,json,resource,signal,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def field(p,key):
 PINS[p]=sha(ROOT/p);s=(ROOT/p).read_text();token=json.dumps(key);assert s.count(token)==1;i=s.index(':',s.index(token))+1
 while s[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(s,i)[0]
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(enc,v))
 return v
def add(p,k,c):
 p[k]=p.get(k,F(0))+c
 if not p[k]:del p[k]
def count(v):
 if isinstance(v,(int,F)):return 1
 if isinstance(v,str):
  try:F(v);return 1
  except ValueError:return 0
 if isinstance(v,dict):return sum(count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 gp='reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_PREFLIGHT_20261105_v63/report.json';zp='reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_ZERO3_SOURCE_20261105_v63/report.json';g=read(gp);z=read(zp)
 for r in[g,z]:
  assert r['all_passed']and not r['actual_rank_dual_or_cube_computed'] if 'actual_rank_dual_or_cube_computed'in r else r['all_passed']
  for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 D=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_PREFLIGHT_20261105_v63/tiny-five-relations-and-constructive-quartic-dual-source.json')
 Z=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_ZERO3_SOURCE_20261105_v63/three-actual-zero-heads-and-corrected-first-source-relations.json')
 cp=[p for p in g['input_sha256']if p.endswith('independent-homogeneous600-closure-and-r-reduction-proof.json')][0]
 q=[{int(k):F(c)for k,c in row.items()}for row in field(cp,'known_cubic_quotient3_by600')]
 assert field(cp,'r_multiplication_degree3_representative_rows600')==[119,479,599]
 assert q==[{int(k):F(c)for k,c in row}for row in D['source_known_cubic_quotient3_by600_sparse']]
 ms3=list(itertools.combinations_with_replacement(range(8),3));ms4=list(itertools.combinations_with_replacement(range(8),4));idx3={m:i for i,m in enumerate(ms3)};heads=[0,3,4];pure4=[330*h+329 for h in heads]
 assert pure4==D['known_pure_r4_basis_rows1650']
 p4=[]
 for target in range(5):
  for mon in ms4:
   k=target*120+idx3[mon[:3]];v=mon[-1];coeff=[row.get(k,F(0))for row in q];outcol=[F(0)]*3;witness=[];relation={}
   add(relation,(target,mon),F(1))
   # First lifted actual degree3 relation times v.
   for h,c in zip(heads,coeff):add(relation,(h,tuple(sorted((7,7,7,v)))),-c)
   for f,c in enumerate(coeff):
    if not c:continue
    second=heads[f]*120+idx3[tuple(sorted((7,7,v)))];b=[row.get(second,F(0))for row in q];witness.append([f,c,second])
    add(relation,(heads[f],tuple(sorted((7,7,v,7)))),c)
    for j,d in enumerate(b):
     outcol[j]+=c*d;add(relation,(heads[j],(7,7,7,7)),-c*d)
   assert [F(dict(D['all1650_constructive_quartic_reduction_columns3_sparse'][len(p4)]).get(i,0))for i in range(3)]==outcol
   saved=D['all1650_two_cubic_relation_witnesses'][len(p4)];assert saved[0:2]==[k,v]and [[int(i),F(c),int(j)]for i,c,j in saved[2]]==witness
   wanted={(target,mon):F(1)}
   for j,c in enumerate(outcol):add(wanted,(heads[j],(7,7,7,7)),-c)
   assert relation==wanted;p4.append(outcol)
 assert all(p4[pure4[j]][i]==int(i==j)for i in range(3)for j in range(3))
 np=[p for p in g['input_sha256']if p.endswith('actual-natural-Told-to-Tnil-source-map-and-origin.json')][0]
 M=field(np,'actual_origin_data')['actual_origin_map_columns5'];assert M==D['source_actual_natural_origin_columns5']
 assert all(all(F(c)==0 for c in M[j])for j in[1,2,3])
 pure=[[row.get(120*h+119,F(0))for h in range(5)]for row in q];assert pure==[[F(c)for c in row]for row in D['source_q3_restricted_to_five_pure_r3_heads']]
 # Verify each old degree3 zero witness independently by multiplication only.
 oldp=[p for p in z['input_sha256']if p.endswith('known-initial40-and-length-upper17.json')][0]
 old=read(oldp);gens=old['homogeneous_module_generators40'];checks=[]
 for proof in Z['three_actual_old_head_r3_zero_source_proofs']:
  back={}
  for j,shift,c in proof['homogeneous_degree3_zero_reduction_trace']:
   for target,mon,a in gens[j]:add(back,(target,tuple(x+y for x,y in zip(mon,shift))),F(c)*F(a))
  assert back=={(proof['oldhead'],(0,0,0,0,0,0,0,3)):F(1)};checks.append(proof['oldhead'])
 assert checks==[1,2,3]
 firstp=[p for p in z['input_sha256']if p.endswith('actual-three-first-columns-and-fixed-quotient-proof.json')][0]
 first=read(firstp);assert first['all576_actual_first_Koszul_checks']and first['old_incorrect_cache_and_rank_diagnostic_withdrawn']
 assert first['actual_old_free56_labels']==[23,29,30]and Z['three_corrected_actual_M1_columns40']==first['actual_three_minimal_first_columns40']
 for p in sorted(set(g['input_sha256'])|set(z['input_sha256'])):
  if p.endswith('/report.json'):assert read(p)['all_passed']
 custody=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_PREPARED_WORKER_20261105_v63/report.json');driver=ROOT/custody['prepared_code'];PINS[custody['prepared_code']]=sha(driver)
 assert sha(driver)==custody['prepared_code_sha256']=='d53306dac45e2ec131e9e776348ee5ac6b2c91c0acd6f67fca1a2d9ef7802606'
 tree=ast.parse(driver.read_text());assert not any(isinstance(n,ast.ImportFrom)and n.module not in['fractions','pathlib']for n in ast.walk(tree))
 helpcheck=subprocess.run([sys.executable,str(driver),'--help'],capture_output=True,text=True,timeout=5);assert helpcheck.returncode==0
 assert not (ROOT/'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63').exists()
 owned=count([g,z,D,Z,q,M,p4,old,first,custody]);assert owned<60000,owned
 data={'all1650_source_columns_and_two_stage_multiplybacks':True,'all9_pure_r4_embedding_products':True,'actual_relation_sources_M0_and_correctedM1_pinned':True,'three_old_r3_zero_multiplybacks':checks,
 'actual_zero_upgrade':'known initial relation gives oldhead*r3 in oldm4; independently accepted actual oldm4=0 makes it exactly zero. M0=0 for those3, so M1*r3 gives true degree4 relation; no withdrawn V1 firsts.',
 'p4_generation_proof':'For every degree4 monomial, first q3 relation times last variable, then each head*r2*last q3 relation times r. These explicit lifted true initial relations express it modulo actual degree4 relations as a linear combination of3 pure-r4 classes.',
 'rank2_completion_proof':'actual gr4 dimension1 and surjection k3→gr4 imply a2dim kernel. If the8 fixed actual relations span rank2, their annihilator composed with p4 is the actual complete dual. Rank0/1 inconclusive; rank3 contradictory and worker asserts.',
 'numeric_job_not_executed_by_reviewer':True,'fixed_future_shape':[3,8],'prepared_worker_code_checked':True,'review_owned_numeric_slots':owned,'worker_static18000_hard24000':True,'source_bit265_hard2048':True,'no_cube_values_authorized_by_this_source_report':True}
 ep=out/'independent-p4-cofactor-zero3-and-fixed-rank-worker-source.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'SOURCE_PASS_constructive_actual_quartic_dual_fixed3x8_gate','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'prepared_code_sha256':sha(driver),'rank_values_or_quartic_dual_evaluated':False,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
