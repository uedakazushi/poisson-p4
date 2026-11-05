#!/usr/bin/env python3
"""Source/positive/custody review only; no signed third-symbol fold."""
import hashlib,json,math,resource,signal,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_INPUT_GATE_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='4fdfae0b28a4a16a3153d4ad27728aa7cb612fb371a1fed32d896153c6293018'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/direct-old-highestthird-origin-dual-source-and-positive-input.json';d=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 support=read('reports/AMBIENT_B18_ELL2_T2_UPPER563_SUPPORT_PREFLIGHT_20261105_v63/one563-coefficient-free-cubic-second-and-unit-carry-ancestor-inventory.json');assert support['case619']==563
 up=read('reports/AMBIENT_B18_ELL2_T2_UPPER2_CUBIC_OUTPUT_SOURCE_20261105_v63/fixed-upper2-one-cubic-output-identity-and-source-frontier.json',['lambda_two']);lam=list(map(sp,up['lambda_two']))
 mug=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes','Mnew0_columns']);B=list(map(dict,mug['old640_source_columns']));ops=dict(mug['old640_typed_closure_nodes']);M=list(map(sp,mug['Mnew0_columns']))
 gamma=[{i:sum((a*M[i].get(j,Q(0)) for j,a in row.items()),Q(0)) for i in range(56)} for row in lam];gamma=[{i:x for i,x in row.items() if x} for row in gamma];assert gamma==list(map(sp,support['Gamma_two']))
 c0=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63/actual-native70-lower-prefix-and-witness-if-any.json')['nonzero_witness']['original768_constant_combination'];assert c0==support['original768_constant_combination']
 P=read('reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',['honest_D0_relation_columns'])['honest_D0_relation_columns'];expected=[];roots=set()
 for channel,row in enumerate(gamma):
  for j,c in c0:
   if j>=640:continue
   for i,g in row.items():
    if i<48:
     k=B[j].get(i)
     if k is None or ops[k][0]=='c':continue
     kind='typed';roots.add(k)
    else:kind='literal_C0';k=P[j][i-48+1]
    expected.append({'channel':channel,'original_source':j,'source_row56':i,'root_kind':kind,'root':k,'coefficient_factors':[str(g),c]})
 assert expected==d['original_factor_terms'];closure=set()
 def visit(k):
  if k in closure:return
  closure.add(k);n=ops[k]
  for child in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3] if n[0]=='product' else []):visit(child)
 for k in roots:visit(k)
 assert [[k,ops[k]] for k in sorted(closure)]==d['typed_source_program'] and len(closure)==45
 atoms={n[1] for i,n in d['typed_source_program'] if n[0]=='actual_C0_scalar_function'};assert len(atoms)==24
 old=dict(read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json',['first640'])['first640']['all79_atom_values_full10']);assert all(i in old and not sp(old[i]).get(0) for i in atoms)
 assert dict(d['actual_scalar_leaf_origins'])=={i:'0' for i in atoms}
 X=read('reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',['full_X0_56_by768'])['full_X0_56_by768']
 for root,row,col,value in d['literalP8_point_origin_source_cells']:assert value=='0' and Q(X[row][col])==0 and any(t['root']==root and t['source_row56']==row and t['original_source']==col for t in expected)
 del X,P
 maximum=0
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());return v
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def add(a,b,c):
  v=dict(a)
  for i,x in b.items():v[i]=plus(v.get(i,(1,0)),mul(c,x))
  return v
 origin={};dual={}
 for i,n in d['typed_source_program']:
  if n[0]=='actual_C0_scalar_function':o=(1,0);v={n[1]:ex(1)}
  elif n[0]=='c':o=ex(Q(n[1],n[2]));v={}
  elif n[0]=='sum':
   o=(1,0);v={}
   for j,a,b in n[1]:c=ex(Q(a,b));o=plus(o,mul(c,origin[j]));v=add(v,dual[j],c)
  else:
   assert n[0]=='product';o=ex(Q(n[1],n[2]));v={}
   for j in n[3]:v=add({k:mul(x,origin[j]) for k,x in v.items()},dual[j],o);o=mul(o,origin[j])
  origin[i]=o;dual[i]=v
 target=[{},{}]
 for t in expected:target[t['channel']]=add(target[t['channel']],dual[t['root']] if t['root_kind']=='typed' else {t['root']:ex(1)},mul(ex(t['coefficient_factors'][0]),ex(t['coefficient_factors'][1])))
 assert [[i,list(map(str,v))] for i,v in origin.items()]==d['positive_node_origin_envelopes']
 assert [[i,[[j,list(map(str,v))] for j,v in sorted(row.items())]] for i,row in dual.items()]==d['positive_node_high3_sensitivity_envelopes']
 assert [[[i,list(map(str,v))] for i,v in sorted(row.items())] for row in target]==d['positive_two_weighted_high3_envelopes'] and maximum==57 and 3*maximum+1==172
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_PREPARED_CUSTODY_20261105_v63/report.json';c=read(cp);assert c['all_passed'];wp=ROOT/c['worker'];assert sha(wp)==c['worker_sha256']=='25f8746b261be0573fc5f57f21a33a881f51b319453424114661f85aa8779355'
 for p,h in c['literal_input_sha256'].items():assert sha(ROOT/p)==h
 z=subprocess.run([sys.executable,str(wp),'--help'],capture_output=True,text=True,check=True);assert not z.stderr and hashlib.sha256(z.stdout.encode()).hexdigest()==c['safe_help_stdout_sha256']
 assert 2*count(d)+6*45+8*24+4000<d['future_whole_numeric_bound']==7618 and d['future_RAW_bound']==1578968
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent45-node24-origin20-factor-formal-only-source-proof.json';ve.write_text(json.dumps({'Gamma_lambda_M0_literal':True,'all20_original_factor_terms_verified':True,'all45_original_nodes_and24_cached_zero_origins_verified':True,'literal_two_P8_cells_verified':True,'all_positive_bounds57_172_rebuilt':True,'worker_and_safe_help_verified':True,'scope':'OnlyformalindependenthighestthirdsymbolsintheDIRECTOLDGammaXc3term. No signedfold/scalarresponse/native/firstsecondcross/b1/A2/M1/uppervalue. Nonzeroformalgradientdoesnotproveanewphysicaljetisneeded.'},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'source_gate_sha256':pins[rp]['sha256'],'no_signed_formal_or_response_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)}};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
