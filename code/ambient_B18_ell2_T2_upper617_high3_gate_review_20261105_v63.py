#!/usr/bin/env python3
"""Fresh V2 source/positive/custody audit; no signed617 fold or high values."""
import hashlib,json,math,resource,signal,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='2c00ddfe0a10dad52925b490f16787c31db03c8d4b032a443c43492469841ef0'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/617-direct-old-highestthird-coherent-source-positive-input.json';d=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 err=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_ERRATUM_20261105_v63/report.json');assert err['all_passed'] and err['old_bytes_preserved'] and err['no_signed_formal_gradient_or_actual_source_values_executed']
 assert sha(ROOT/'code/ambient_B18_ell2_T2_upper617_high3_input_20261105_v63.py')==err['old_code_sha256'];assert sha(ROOT/'reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_20261105_v63/report.json')==err['withdrawn_positive_gate_sha256']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_BOUNDARY_PLANE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 assert read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 sup=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_SUPPORT_PREFLIGHT_20261105_v63/factored617-lower-kernel-direct-third-support-and-source-frontier.json');p8=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_P8_ANCESTOR_SOURCE_20261105_v63/literal86-P8-composite-program-and-exact-boundary-catalog.json')
 assert d['typed_source_program']==sup['pruned_typed_source_closure'] and d['literal_C0_program']==p8['selected_original_C0_arithmetic_nodes'];typed=dict(d['typed_source_program']);c0=dict(d['literal_C0_program']);assert len(typed)==629 and len(c0)==532
 proof=read('reports/AMBIENT_B18_ELL2_T2_LOWER_RANK2_HF_TAIL_SOURCE_20261105_v63/complete1012-kernel617-and-filtered-r-tail-source-proof.json',['all617_kernel_lower_constant_source_recipes'])['all617_kernel_lower_constant_source_recipes'];assert len(proof)==617
 recs=[]
 for rec in proof:
  factors=[{'case619':rec['free_case619'],'scale_operator':['literal','1']}]+[{'case619':j,'scale_operator':['negate_literal',x]} for j,x in zip([563,609],rec['subtract_actual_basis563609'],strict=True) if Q(x)]
  recs.append({'free_case619':rec['free_case619'],'constant_case619_source_factors':factors})
 assert recs==d['case_source_recipes617']
 original=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations'];assert [v['original768_combination'] for v in original]==d['original619_source768_combinations'];C=d['original619_source768_combinations']
 mug=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes','Mnew0_columns']);B=list(map(dict,mug['old640_source_columns']));allops=dict(mug['old640_typed_closure_nodes']);M=list(map(sp,mug['Mnew0_columns']))
 assert all(allops[i]==n for i,n in typed.items())
 lam=list(map(sp,read('reports/AMBIENT_B18_ELL2_T2_UPPER2_CUBIC_OUTPUT_SOURCE_20261105_v63/fixed-upper2-one-cubic-output-identity-and-source-frontier.json',['lambda_two'])['lambda_two']));gamma=[{i:z for i in range(56) if (z:=sum((a*M[i].get(j,Q(0)) for j,a in row.items()),Q(0)))} for row in lam];assert gamma==list(map(sp,d['Gamma_two']))
 P=read('reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',['honest_D0_relation_columns'])['honest_D0_relation_columns'];need=sorted({j for rec in recs for f in rec['constant_case619_source_factors'] for j,x in C[f['case619']] if j<640});pruned={}
 for j in need:
  root={}
  for row in gamma:
   for i in row:
    if i<48:
     k=B[j].get(i)
     if k is None:continue
     n=allops[k]
     if n[0]=='c':tag=['constant',n[1],n[2]]
     elif n[0]=='actual_C0_scalar_function':tag=['C0',n[1]]
     else:assert k in typed;tag=['typed',k]
    else:tag=['C0',P[j][i-48+1]]
    root[i]=tag
  pruned[j]=sorted(root.items())
 assert sorted(pruned.items())==[(j,[(i,v) for i,v in row]) for j,row in d['pruned_original_old_columns_roots']] and len(need)==618
 a79={i:sp(v).get(0,Q(0)) for i,v in read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json',['first640'])['first640']['all79_atom_values_full10']}
 a164={i:sp(v).get(0,Q(0)) for i,v in read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',['all_requested164_root_values'])['all_requested164_root_values']}
 leaves=({n[1] for n in typed.values() if n[0]=='actual_C0_scalar_function'}|{i for i,kind,n in p8['exact_literal_cache_boundary_nodes']})-set(c0);assert len(leaves)==96
 origins={};classes=[]
 for i in sorted(leaves):
  assert i in a79 or i in a164
  if i in a79 and i in a164:assert a79[i]==a164[i]
  origins[i]=a79[i] if i in a79 else a164[i];classes.append([i,'source164' if i in a164 else 'C0_79_only',str(origins[i])])
 assert classes==d['literal_boundary_origins_and_membership'] and sum(k=='source164' for i,k,x in classes)==48 and sum(k=='C0_79_only' for i,k,x in classes)==48
 maximum=0
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());return v
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def add(a,b,c=(1,1)):
  z=dict(a)
  for i,x in b.items():
   v=plus(z.get(i,(1,0)),mul(c,x))
   if v[1]:z[i]=v
  return z
 def scale(a,c):return {i:mul(x,c) for i,x in a.items() if x[1] and c[1]}
 co={};cg={}
 def CC(i):
  if i in co:return co[i],cg[i]
  assert i in origins;return ex(origins[i]),{i:ex(1)}
 for i,n in sorted(c0.items()):
  if n[0]=='c':o=ex(Q(n[1],n[2]));g={}
  elif n[0]=='n':o,g=CC(n[1]);g=dict(g)
  else:
   a,u=CC(n[1]);b,v=CC(n[2])
   if n[0]=='a':o=plus(a,b);g=add(u,v)
   else:assert n[0]=='m';o=mul(a,b);g=add(scale(u,b),v,a)
  co[i]=o;cg[i]=g
 to={};tg={}
 for i,n in sorted(typed.items()):
  if n[0]=='actual_C0_scalar_function':o,g=CC(n[1]);g=dict(g)
  elif n[0]=='c':o=ex(Q(n[1],n[2]));g={}
  elif n[0]=='sum':
   o=(1,0);g={}
   for j,a,b in n[1]:s=ex(Q(a,b));o=plus(o,mul(s,to[j]));g=add(g,tg[j],s)
  else:
   assert n[0]=='product';o=ex(Q(n[1],n[2]));g={}
   for j in n[3]:g=add(scale(g,to[j]),tg[j],o);o=mul(o,to[j])
  to[i]=o;tg[i]=g
 col={}
 for j,roots in pruned.items():
  roots=dict(roots);pair=[]
  for row in gamma:
   z={}
   for i,x in row.items():
    if i not in roots:continue
    tag=roots[i]
    if tag[0]=='constant':v={}
    elif tag[0]=='C0':v=CC(tag[1])[1]
    else:assert tag[0]=='typed' and tag[1] in tg;v=tg[tag[1]]
    z=add(z,v,ex(x))
   pair.append(z)
  col[j]=pair
 def case(j,ch):
  z={}
  for i,x in C[j]:
   if i<640:z=add(z,col[i][ch],ex(x))
  return z
 bas={j:[case(j,ch) for ch in range(2)] for j in [563,609]};global2=[(1,0)]*2;output_support=0
 for rec in recs:
  for ch in range(2):
   z={}
   for f in rec['constant_case619_source_factors']:
    j=f['case619'];z=add(z,bas[j][ch] if j in bas else case(j,ch),ex(f['scale_operator'][1]))
   output_support+=len(z)
   for D,U in z.values():
    a,b=global2[ch];dd=math.lcm(a,D);global2[ch]=mark((dd,max(b*(dd//a),U*(dd//D))))
 ngrad=sum(map(len,cg.values()))+sum(map(len,tg.values()));cgrad=sum(sum(map(len,pair)) for pair in col.values())
 assert ngrad==d['positive_sparse_node_gradient_entries']==218 and cgrad==d['positive_sparse_old_column_gradient_entries']==136 and output_support==d['positive617_two_row_gradient_entries']==322
 assert [[str(a),str(b)] for a,b in global2]==d['positive_common_two_row_envelopes'] and maximum==76 and 3*maximum+1==d['temporary_prior_bits']==229
 # Explicit bound on the actual worker method, using FINAL serialized input.
 owned=2*count(d)+(3*1161+2*218)+(618+2*136)+(2+8*96)+2*(617+2*322+2*96+30)+(3*96+40)+1000+4000;assert owned==38097<d['future_formal_numeric_bound']==40920 and d['future_formal_RAW_bound']==5519120
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_PREPARED_CUSTODY_20261105_v63/report.json';cust=read(cp);assert cust['all_passed'];wp=ROOT/cust['worker'];assert sha(wp)==cust['worker_sha256']=='050071b528172309b5cf5bcd2ca36fbac1173b49a1c7518b2b2c8bb2292e0bf8'
 for p,h in cust['literal_input_sha256'].items():assert sha(ROOT/p)==h
 z=subprocess.run([sys.executable,str(wp),'--help'],capture_output=True,text=True,check=True);assert not z.stderr and hashlib.sha256(z.stdout.encode()).hexdigest()==cust['safe_help_stdout_sha256']
 code=wp.read_text();assert 'tg.get(k,{})' not in code and "assert i in origins" in code and "assert tag[0]=='typed' and tag[1] in tg" in code
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent-V2-617-source-classification-positive-and-worker-proof.json';ve.write_text(json.dumps({'all617_lower_kernel_factors_exact':True,'Gamma_lambda_M0_same_frame':True,'all618_original_old_column_root_tags_exact':True,'direct_C0_roots_share_same96_symbols':True,'source164_boundary_ids':[i for i,k,x in classes if k=='source164'],'only79_boundary_ids':[i for i,k,x in classes if k=='C0_79_only'],'all_positive_node_column617_stats_and_envelopes_rebuilt':[218,136,322],'positive_committed_bits':76,'temporary_bits':229,'actual_worker_final_input_owned_bound':owned,'whole_cap':40920,'literal_worker_and_help_and_withdrawn_V1_preserved':True,'no_signed617_fold_or_higher_coefficient_values':True,'scope':'OnlyformalindependenthighestthirdsymbolweightsintheDIRECTOLDterm. Actualthirds/native/knowncross/b1/A2/M1/upperresponseareunassigned. Source164plane9availabilitydoesnotextendto48only79symbols. No necessityclaim.'},separators=(',',':'))+'\n');report={'id':STEM,'all_passed':True,'source_gate_sha256':pins[rp]['sha256'],'formal_only_source_code_custody_condition_closed':True,'higher_scalar_or_upper_values_not_queried':True,'positive_bits':76,'temporary_bits':229,'actual_worker_owned_bound':owned,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)}};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owned}))
if __name__=='__main__':main()
