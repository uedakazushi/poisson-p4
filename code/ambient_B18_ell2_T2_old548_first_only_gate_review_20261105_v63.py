#!/usr/bin/env python3
"""Independent first-only548 SOURCE/worker review; no image values."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as OWN
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_INPUT_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(v):return {int(i):F(x) for i,x in v if F(x)}
def count(v):
 if isinstance(v,(int,F)):return 1
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(tuple,list,set)):return sum(map(count,v))
 return 0
def main(out):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=sha(ROOT/p);return OWN.selected(p,k)
 gp='reports/AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_INPUT_20261105_v63/report.json';g=read(gp);assert hs[gp]=='b74e5722306e454fec3b2b8d199d7df359631075c340634f694a785ac5051679' and g['all_passed'];ip=str(Path(gp).parent/'pruned548-first-only-cache-source-and-positive-bound.json');d=read(ip);assert hs[ip]==g['evidence_sha256'][Path(ip).name]
 for r in g['input_sha256'].values():assert sha(ROOT/r['path'])==r['sha256']
 cu=read('reports/AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_PREPARED_CUSTODY_20261105_v63/report.json');assert cu['worker_sha256'].startswith('1c003578') and sha(ROOT/cu['worker'])==cu['worker_sha256']
 formal=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63/actual-formal619-oldsecond-weights-and-symbol-support.json');assert read('reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63_REVIEW/report.json')['all_passed'];nonzero=set(formal['nonzero_case619_indices'])
 C=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations');selected=[]
 for i,r in enumerate(C):
  cc=sp(r['original768_combination'])
  if i in nonzero or any(j>=640 for j in cc):continue
  assert len(cc)==1 and next(iter(cc.values()))==1 and all(k%9<7 for k in sp(r['actual_a1_27_times9']));selected.append({'case619':i,'original_old640_source':next(iter(cc)),'actual_a1_27_times9':r['actual_a1_27_times9']})
 assert selected==d['selected548_source_records'] and len(selected)==548 and sum(bool(r['actual_a1_27_times9']) for r in selected)==32 and sum(len(r['actual_a1_27_times9']) for r in selected)==103
 gg='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json';Gamma=list(map(sp,field(gg,'Gamma4_on_source56_times8')));assert d['Gamma4_source56_times8']==field(gg,'Gamma4_on_source56_times8');assert read('reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_GATE_20261105_v63_REVIEW/report.json')['all_passed'];rows=sorted({k//8 for v in Gamma for k in v});assert max(rows)<48
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';B=list(map(dict,field(mug,'old640_source_columns')));ops=dict(field(mug,'old640_typed_closure_nodes'));pb={r['original_old640_source']:{i:B[r['original_old640_source']][i] for i in rows if i in B[r['original_old640_source']]} for r in selected};assert pb=={j:dict(v) for j,v in d['oldB_source_rows_pruned']}
 closure=set();todo=[k for v in pb.values() for k in v.values()]
 while todo:
  k=todo.pop()
  if k in closure:continue
  closure.add(k);n=ops[k]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0] in ['c','actual_C0_scalar_function']
 assert len(closure)==636 and {k:ops[k] for k in closure}==dict(d['typed_closure_nodes']);leaves={ops[k][1] for k in closure if ops[k][0]=='actual_C0_scalar_function'};assert len(leaves)==62
 fp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';cache=field(fp,'first640');allfirst=dict(cache['all79_atom_values_full10']);first={i:sp(allfirst[i]) for i in leaves};assert first=={i:sp(v) for i,v in d['actual_C0_leaf_first10']};assert all(not v.get(0) for v in first.values())
 # Literal source-interface check: the full10 cache ALREADY uses new r/ell.
 # Old bar8 matchesindices0..8; independent source(ds)/5 givesindex9.
 oldatoms=dict(field('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json','all79_atom_values'));boundary={r['original_node']:sp(r['origin_and_ds']) for r in field(mug,'all_exact_cached_leaf_bindings_origin_and_ds')}
 def pairmul(a,b):
  v={};q=a.get(0,F(0))*b.get(0,F(0));r=a.get(0,F(0))*b.get(1,F(0))+a.get(1,F(0))*b.get(0,F(0))
  if q:v[0]=q
  if r:v[1]=r
  return v
 for i,n in field(mug,'old_exact84_original_arithmetic_nodes'):
  if n[0]=='c':v={0:F(n[1],n[2])} if n[1] else {}
  elif n[0]=='n':v={k:-x for k,x in boundary[n[1]].items()}
  elif n[0]=='a':v={k:x for k in set(boundary[n[1]])|set(boundary[n[2]]) if (x:=boundary[n[1]].get(k,F(0))+boundary[n[2]].get(k,F(0)))}
  else:assert n[0]=='m';v=pairmul(boundary[n[1]],boundary[n[2]])
  boundary[i]=v
 for i,v in first.items():assert {k:x for k,x in v.items() if k<9}==sp(oldatoms[i]) and v.get(9,F(0))==boundary[i].get(1,F(0))/5
 piv=field('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json','original27_X_source_columns');assert piv==d['point_pivot27_original_source_ids'];needed=sorted({k//9 for r in selected for k,x in r['actual_a1_27_times9']});assert needed==d['needed_pivot_indices']
 R={}
 for off,path,key in [(0,OWN.PATHS['old640'],'source_columns640_rows56_sparse_dual9'),(640,OWN.PATHS['old128'],'source_columns128_rows48_sparse_dual9')]:
  raw=field(path,key)
  for j in needed:
   idx=piv[j]-off
   if 0<=idx<len(raw):R[j]={row:sp(v)[8] for row,v in raw[idx] if row in rows and sp(v).get(8)}
 assert R=={p:sp(v) for p,v in d['actual_pivot_r_first_pruned']}
 olde=field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json','actual_first640_ell_columns');E={p:{i:x for i,x in sp(olde[piv[p]]).items() if i in rows} for p in needed if piv[p]<640};assert E=={p:sp(v) for p,v in d['actual_old_pivot_ell_first_pruned']}
 native={p:piv[p]-640 for p in needed if piv[p]>=640};assert native==dict(d['native_pivot_local128_ids']);L0=[{i:x for i,x in sp(v).items() if i in rows} for v in field(mug,'L0_Res4_Ntranspose_columns')];Le=[{i:x for i,x in sp(v).items() if i in rows} for v in field('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json','actual_L1ell_columns')];Se0=field(mug,'native_S1ell_columns');S00=field(mug,'native_S0_columns');Se={p:sp(Se0[j]) for p,j in native.items()};S0={p:sp(S00[j]) for p,j in native.items()}
 assert L0==list(map(sp,d['native_L0_rows_pruned'])) and Le==list(map(sp,d['native_L1ell_rows_pruned'])) and Se=={p:sp(v) for p,v in d['native_pivot_S1ell_columns']} and S0=={p:sp(v) for p,v in d['native_pivot_S0_columns']}
 fixed=count([d,selected,Gamma,ops,B,cache,first,boundary,R,E,L0,Le,Se,S0])+10000;assert fixed<150000
 (out/'phase.json').write_text(json.dumps({'phase':'SOURCE_positive_polynomial_check_only_no_actual548_images','fixed_plus_scratch':fixed})+'\n')
 pairs=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)];assert d['ordinary_pair_order8']==[list(v) for v in pairs];allowed={()}|{(i,) for i in range(9)}|set(pairs);maximum=0
 def mark(v):
  nonlocal maximum
  maximum=max(maximum,v[0].bit_length(),v[1].bit_length());return v
 def ex(v):v=F(v);return mark((v.denominator,abs(v.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def accum(v,k,b):
  if b[1]:v[k]=plus(v.get(k,(1,0)),b)
 def product(a,b):
  v={}
  for u,x in a.items():
   for w,y in b.items():
    q=tuple(sorted(u+w))
    if q in allowed:accum(v,q,mul(x,y))
  return v
 jets={}
 for i in sorted(closure):
  n=ops[i]
  if n[0]=='actual_C0_scalar_function':v={(k-1,) if k else ():ex(x) for k,x in first[n[1]].items()}
  elif n[0]=='c':v={():ex(F(n[1],n[2]))} if n[1] else {}
  elif n[0]=='sum':
   v={}
   for j,a,b in n[1]:
    for word,x in jets[j].items():accum(v,word,mul(ex(F(a,b)),x))
  else:
   assert n[0]=='product';v={():ex(F(n[1],n[2]))} if n[1] else {}
   for j in n[3]:v=product(v,jets[j])
  jets[i]=v
 EP={p:{i:ex(x) for i,x in v.items()} for p,v in E.items()}
 for p in native:
  v={}
  for A,b in [(L0,Se[p]),(Le,S0[p])]:
   for j,x in b.items():
    for i,y in A[j].items():accum(v,i,mul(ex(x),ex(y)))
  EP[p]=v
 globalDU=[(1,0)]*4
 for rec in selected:
  a1=sp(rec['actual_a1_27_times9']);cb=pb[rec['original_old640_source']]
  for ch,row in enumerate(Gamma):
   val=(1,0)
   for k,x in row.items():
    target,p=divmod(k,8);u,v=pairs[p];root=cb.get(target)
    if root is not None:val=plus(val,mul(ex(x),jets[root].get((u,v),(1,0))))
    if u<7:
     for j in needed:
      a=a1.get(9*j+u,F(0))
      if a:val=plus(val,mul(ex(x),mul(ex(a),ex(R[j].get(target,F(0))) if v==7 else EP[j].get(target,(1,0)))))
   D,U=globalDU[ch];dd,uu=val;nd=math.lcm(D,dd);globalDU[ch]=mark((nd,max(U*(nd//D),uu*(nd//dd))))
 assert maximum==60 and [[str(D),str(U)] for D,U in globalDU]==d['all548_four_images_global_D_U'] and d['positive_temporary_bit_bound']==181
 assert count(d)==7974 and 2*count(d)+18*len(closure)+8*548+6000==37780<g['future_numeric_bound']==40368
 worker=ROOT/cu['worker'];code=worker.read_text();assert 'for A,b in [(L0,Se[p]),(Le,S0[p])]' in code and 'if any(values)' in code and 'actual_function_second_parts_not_set_zero' in code and 'jets[child]' in code;helpres=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,timeout=10);assert helpres.returncode==0
 ep=out/'independent548-source-firsts-and-positive-polynomial-proof.json';ep.write_text(json.dumps({'all548_single_old_sources_complete':True,'all_four_formal_oldsecond_weights_zero_pinned':True,'no_native_X2_sources':True,'all103_transverse_a1_terms32cases_verified':True,'all636_original_typed_nodes_and62literal_firsts_verified':True,'full9_first_frame':'indices0,z0..6,rA,ell; oldbarfirst0..8 literalcompare andell=ds/5 independentlyreplayed','all_individual_pivot_r_and_ell_first_inputs_verified':True,'individual_native_pivot_ell_keeps_L0Se_plus_LeS0':True,'ordinary_selected_poly_multiplication_no_square_factor2':True,'positive_prior60_reproduced_no_signed_image_evaluated':True,'global_four_D_U':[[str(D),str(U)] for D,U in globalDU],'serialized_numeric_including_metadata':7974,'live_two_copies_jets_prefix_scratch':37780,'unchanged_future_numeric':40368,'zero_extension_scope':'Onlythewhole548functional,afterall4honestsecondweightscancel. Noactualleafsecondcoefficientassignedzero. Knownfirstcrossproductsretained.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'selected548_source_and_zero_weight_eligibility_verified':True,'all62_cached_full9_firsts_verified':True,'all_native_pivot_LeS0_terms_retained':True,'positive_prior_bits':60,'temporary_prior_bits':181,'future_whole_numeric':40368,'future_whole_RAW':5148704,'worker_sha256':sha(worker),'safe_help_exit':helpres.returncode,'no_actual548_image_or_unknown_coefficient_evaluated':True,'input_sha256':hs,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','selected548_source_and_zero_weight_eligibility_verified','positive_prior_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
