#!/usr/bin/env python3
"""Independent aliases/origin/formal-factor input review; no signed fold."""
import argparse,collections,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as OWN
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_GATE_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(v):return {int(i):Q(c) for i,c in v if Q(c)}
def count(v):
 if isinstance(v,(int,Q)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(k)+count(w) for k,w in v.items())
 if isinstance(v,(tuple,list,set)):return sum(map(count,v))
 return 0
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=sha(ROOT/p);return OWN.selected(p,k)
 gp='reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_GATE_20261105_v63/report.json';g=read(gp);assert hs[gp]=='c591e50d54e1d536d0ef2f25b79e7cdea8dd5dedb7c5efae889bd2b31f665c6a' and g['all_passed']
 dp=str(Path(gp).parent/'factored619-old-honest-second-linear-input-and-positive-bound.json');d=read(dp);assert hs[dp]==g['evidence_sha256'][Path(dp).name]
 for r in g['input_sha256'].values():assert sha(ROOT/r['path'])==r['sha256']
 custody=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_PREPARED_CUSTODY_20261105_v63/report.json');assert custody['worker_sha256']=='d5f39093aa493f756b2fb5d2a7569d830a7b135bca842b9d692a06b1101fc68a' and sha(ROOT/custody['worker'])==custody['worker_sha256']
 for p,h in custody['literal_inputs'].items():assert sha(ROOT/p)==h
 assert count(d)==9458 and 2*count(d)+10000==28916<g['future_formal_numeric_bound']==31822
 ap='reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_ALIAS_ORIGIN_20261105_v63/twelve-arithmetic-alias-closure-and-zero-origin-bindings.json';aliases=read(ap)
 prof=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_CALLBACK_PROFILE_20261105_v63/literal-oldB-native-X2-callback-and-firstcarry-profile.json')
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';known=field(mug,'all_exact_cached_leaf_bindings_origin_and_ds');K={r['original_node']:r for r in known};ops=dict(field(mug,'old640_typed_closure_nodes'));B=[dict(v) for v in field(mug,'old640_source_columns')]
 C=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations');assert len(C)==619 and all(r['full56_Xgamma1_equals_Xp0a1'] for r in C)
 G=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json');Gamma=list(map(sp,G['Gamma4_on_source56_times8']))
 originalM=list(map(sp,field(mug,'Mnew0_columns')));PB=list(map(sp,field('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json','PiBottom3_from56_columns')))
 adj=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63/actual-two-mixed-output-adjoint-certificates-or-obstruction.json');assert read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63_REVIEW/report.json')['all_passed'];rho=[{},{}]+list(map(sp,adj['solutions_on136_literal_outputs']))
 fixed=count([g,d,aliases,prof,known,ops,B,C,G,Gamma,originalM,PB,adj])+20000;assert fixed<150000
 (out/'phase.json').write_text(json.dumps({'phase':'independent_known_source_and_formal_factor_checks_only','owned_fixed_plus_scratch':fixed,'no_signed_weighted_gradient_or_function_values':True})+'\n')
 # Recover all four Gamma rows from original M0/Pi3 and the independently
 # checked adjoints. In particular retain the two DIRECT chiA rows.
 duals=[[Q(0),Q(1),Q(4)],[Q(11),Q(0),Q(270)]];ownGamma=[]
 for channel in range(4):
  ell=duals[channel//2];pair=channel%2;row={}
  for j in range(56):
   direct=sum((x*PB[j].get(k,Q(0)) for k,x in enumerate(ell)),Q(0))
   for m in range(8):
    val=(direct if m==pair else Q(0))-sum((x*originalM[j].get(t,Q(0)) for flat,x in rho[channel].items() for t,p in [divmod(flat,8)] if p==m),Q(0))
    if val:row[8*j+m]=val
  ownGamma.append(row)
 assert ownGamma==Gamma and all(k//8<48 for row in Gamma for k in row)
 # Literal original source roots, all619/fourchannels. This is support only,
 # so no old root coefficient or weighted gradient is formed.
 roots=set()
 for rec in C:
  for j,c in rec['original768_combination']:
   if j>=640:continue
   for row in Gamma:
    for k in row:
     root=B[j].get(k//8)
     if root is not None:roots.add(root)
 assert roots==set(prof['old_distinct_roots']) and len(roots)==576
 reachable=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in reachable:continue
  reachable.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0] in ['c','actual_C0_scalar_function']
 assert reachable==set(prof['reachable_typed_node_ids']) and len(reachable)==668
 leaves={ops[i][1] for i in reachable if ops[i][0]=='actual_C0_scalar_function'};assert len(leaves)==62 and leaves==set(prof['reachable_actual_C0_scalar_labels'])
 missing=leaves-set(K);nodes=dict(aliases['selected_arithmetic_nodes']);assert missing==set(aliases['missing12_literal_nodes']) and len(missing)==12 and len(nodes)==36
 # Independently verify each of the36 literal alias AST nodes by streaming
 # the immutable original DAG; retain neither a whole AST nor fresh jets.
 cp=OWN.PATHS['C0'];hs[cp]=sha(ROOT/cp);text=(ROOT/cp).read_text();assert len(text.encode())<12*1024*1024;pos=text.index('[',text.index('"DAG_scalar_nodes"'))+1;dec=json.JSONDecoder();wanted=set(nodes)
 for i in range(max(wanted)+1):
  while text[pos].isspace() or text[pos]==',':pos+=1
  node,pos=dec.raw_decode(text,pos)
  if i in wanted:assert node==nodes[i]
  del node
 del text
 vals={i:sp(r['origin_and_ds']).get(0,Q(0)) for i,r in K.items()};aliasops=set()
 for i,n in sorted(nodes.items()):
  op=n[0];aliasops.add(op)
  if op=='c':v=Q(n[1],n[2])
  elif op=='a':assert n[1] in vals and n[2] in vals;v=vals[n[1]]+vals[n[2]]
  elif op=='m':assert n[1] in vals and n[2] in vals;v=vals[n[1]]*vals[n[2]]
  elif op=='n':v=-vals[n[1]]
  else:assert op=='i' and vals[n[1]];v=1/vals[n[1]]
  assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=384;vals[i]=v
 assert all(i in vals and vals[i]==0 for i in leaves)
 for rec in aliases['cached_boundary_records_used']:assert rec==K[rec['original_node']]
 # Formal derivative paths are products of literal rational constants.
 # Zero-origin nonconstant factors remove higher products' honest-second
 # sensitivities, but their first cross-products are explicitly not touched.
 origin={};paths={}
 for i in sorted(reachable):
  n=ops[i];op=n[0]
  if op=='actual_C0_scalar_function':origin[i]=vals[n[1]];terms=[(n[1],())]
  elif op=='c':origin[i]=Q(n[1],n[2]);terms=[]
  elif op=='sum':
   origin[i]=sum((Q(a,b)*origin[j] for j,a,b in n[1]),Q(0));terms=[]
   for j,a,b in n[1]:terms.extend((label,fac+(str(Q(a,b)),)) for label,fac in paths[j])
  else:
   assert op=='product';origin[i]=Q(n[1],n[2])*math.prod(origin[j] for j in n[3]);terms=[]
   for slot,j in enumerate(n[3]):
    others=[u for k,u in enumerate(n[3]) if k!=slot]
    if any(origin[u]==0 for u in others):continue
    assert all(ops[u][0]=='c' for u in others);f=(str(Q(n[1],n[2])),)+tuple(str(origin[u]) for u in others);terms.extend((label,fac+f) for label,fac in paths[j])
  assert op=='c' or origin[i]==0;paths[i]=terms
 sensitive=sorted(i for i in roots if paths[i]);assert len(sensitive)==13 and sensitive==d['nonzero_linear_parts_only13_raw_roots']
 def canon(case,ch,pair,label,source,root,factors):return (case,ch,pair,label,source,root,tuple(sorted(str(Q(x)) for x in factors)))
 expected=collections.Counter()
 for case,rec in enumerate(C):
  for source,c in rec['original768_combination']:
   if source>=640:continue
   for ch,row in enumerate(Gamma):
    for k,a in row.items():
     target,pair=divmod(k,8);root=B[source].get(target)
     if root is None:continue
     for label,f in paths[root]:expected[canon(case,ch,pair,label,source,root,(str(c),str(a))+f)]+=1
 actual=collections.Counter(canon(r['case619'],r['channel4'],r['ordinary_pair8'],r['original_C0_symbol'],r['literal_old640_source'],r['literal_typed_root'],r['rational_factors']) for r in d['all619_four_channel_old_second_factor_records']);assert expected==actual and sum(expected.values())==702
 # Independent positive envelopes only. There is NO signed weighted fold.
 bounds={};prior=0
 for r in d['all619_four_channel_old_second_factor_records']:
  D,U=1,1
  for x in r['rational_factors']:
   q=Q(x);D*=q.denominator;U*=abs(q.numerator);prior=max(prior,D.bit_length(),U.bit_length())
  key=(r['case619'],r['channel4'],r['ordinary_pair8'],r['original_C0_symbol']);dd,uu=bounds.get(key,(1,0));nd=math.lcm(D,dd);U=uu*(nd//dd)+U*(nd//D);bounds[key]=(nd,U);prior=max(prior,nd.bit_length(),U.bit_length())
 saved={tuple(k):(int(D),int(U)) for k,D,U in d['all_positive_symbol_coefficient_envelopes']};assert bounds==saved and len(bounds)==576 and prior==46
 labels=sorted({k[3] for k in bounds});assert labels==d['all_independent_C0_second_symbols']==g['surviving_original_symbols']
 worker=ROOT/custody['worker'];code=worker.read_text();assert 'grad.get(key,Q(0))+q' in code and 'no_actual_honestsecond_native_gauge_firstcross_or_pivotcarry_response' in code
 help_result=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,timeout=10);assert help_result.returncode==0
 owned=fixed+count([vals,nodes,paths,origin,expected,bounds])+10000;assert owned<150000
 ep=out/'independent668-origin-source-factors-and-positive-prior-checks.json';ep.write_text(json.dumps({'all_four_Gamma_rebuilt_from_original_M0_Pi3_and_verified_rho':True,'direct_chiA_channels_retained':True,'roots576_and_reachable668_and_leaves62_verified':True,'literal36_alias_nodes_stream_verified':True,'all62_original_leaf_origins_zero':True,'alias_operation_types':sorted(aliasops),'nonconstant_typed_origins_allzero':True,'thirteen_linear_sensitive_raw_roots':sensitive,'all702_factored_terms_verified_including_all619_four_channels':True,'all576_positive_envelopes_verified':True,'positive_prior_bits':prior,'serialized_input_including3metadata':count(d),'two_input_copies_plus10000':2*count(d)+10000,'unchanged_future_numeric':31822,'all13_literal_second_labels_retained':labels,'no_signed_weighted_gradient_or_unknown_response_evaluated':True,'product_rule_scope':'Products of two positive-order factors have no linear sensitivity to a single honest-second symbol; their known first cross-products remain untouched.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all619_four_channel_source_factors_verified':True,'all62_zero_origins_and_alias_literals_verified':True,'formal_factor_terms':702,'weighted_symbol_slots':576,'positive_prior_bits':46,'temporary_prior_bits':139,'whole_future_numeric':31822,'whole_future_RAW':4163896,'worker_sha256':sha(worker),'worker_safe_help_exit':help_result.returncode,'no_signed_gradient_or_function_values_run':True,'input_sha256':hs,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'owned_numeric_bound':owned,'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','formal_factor_terms','weighted_symbol_slots','positive_prior_bits','owned_numeric_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
