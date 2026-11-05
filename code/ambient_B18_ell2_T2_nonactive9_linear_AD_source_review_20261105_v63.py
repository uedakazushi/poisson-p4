#!/usr/bin/env python3
"""Independent18psi factor/source/prior review; no signed sensitivity fold."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,add,image
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_NONACTIVE9_LINEAR_AD_INPUT_GATE_V2_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):
  inputs[p]={'path':p,'sha256':sha(ROOT/p)}
  return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and not r['signed_sensitivities_or_values_evaluated']
 assert r['code_sha256']==sha(ROOT/'code/ambient_B18_ell2_T2_nonactive9_linear_ad_gate_v2_20261105_v63.py')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/eighteen-factored-old-linear-second-and-native-source-inputs.json';d=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 up='reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_UPPER_SUPPORT_PREFLIGHT_20261105_v63/eighteen-formal-weights-and-coefficient-free-source-ancestor-inventory.json';u=read(up)
 low=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63/exact171-decomposition-162-zero-9-lower-first-source-input.json')
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 m=read(mug,['all741_zero_origin_source_combinations','old640_source_columns','old640_typed_closure_nodes','native_S0_columns'])
 mu=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json',['actual_selected_Y0_source_columns'])
 C=[sp(m['all741_zero_origin_source_combinations'][j]['original768']) for j in mu['actual_selected_Y0_source_columns']]
 B=[dict(v) for v in m['old640_source_columns']];catalog=dict(m['old640_typed_closure_nodes']);S0=list(map(sp,m['native_S0_columns']))
 assert len(catalog)==1302
 assert all(not image(S0,{j-640:x for j,x in c.items() if j>=640}) for c in C)
 P=read('reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',['honest_D0_relation_columns'])['honest_D0_relation_columns']
 gam=read('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json',['Gamma56_two_channels'])
 Gamma=list(map(sp,gam['Gamma56_two_channels']))
 ops=dict(u['typed_closure_nodes']);assert len(ops)==8 and ops==dict(d['typed_linear_root_provenance'])
 assert all(catalog[k]==node for k,node in ops.items())
 leaf={}
 for k,node in sorted(ops.items()):
  if node[0]=='actual_C0_scalar_function':v={node[1]:Q(1)}
  elif node[0]=='c':v={}
  elif node[0]=='product':
   assert len(node[3])==1;v={j:Q(node[1],node[2])*x for j,x in leaf[node[3][0]].items()}
  else:
   assert node[0]=='sum';v={}
   for child,a,b in node[1]:v=add(v,leaf[child],Q(a,b))
  leaf[k]=v
 assert len({label for node in ops.values() if node[0]=='actual_C0_scalar_function' for label in [node[1]]})==4
 expected_weights=[]
 for rec in low['remaining9_sources']:
  br,be=sp(rec['b1_r_selected99']),sp(rec['b1_ell_selected99']);wa={2*j:x for j,x in br.items()};wb={2*j:x for j,x in be.items()}
  for j,x in br.items():wb=add(wb,{2*j+1:x})
  expected_weights.extend([(rec['source_index'],'A',wa),(rec['source_index'],'B',wb)])
 rows=d['all18_factored_linear_second_sensitivity_inputs'];assert len(rows)==18
 maximum=0;oldterms=nativeterms=0;labels=set()
 def mark(b):
  nonlocal maximum
  maximum=max(maximum,b[0].bit_length(),b[1].bit_length());return b
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1]))
 def plus(a,b):
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 for rec,inventory,(idx,ch,weights) in zip(rows,u['eighteen_formal_contractions'],expected_weights,strict=True):
  assert(rec['original418_source'],rec['contraction'])==(inventory['original418_source'],inventory['contraction'])==(idx,ch)
  assert sp(rec['formal_weights_on99psi2'])==weights==sp(inventory['formal_weights_on99psi2'])
  terms=[];native=[];bounds={};nbounds={}
  for flat,w in rec['formal_weights_on99psi2']:
   j,channel=divmod(flat,2)
   for source,c in C[j].items():
    if source>=640:
     f=[w,str(c)];native.append([channel,source-640,f]);b=mul(ex(w),ex(c));key=(channel,source-640);nbounds[key]=plus(nbounds.get(key,(1,0)),b);continue
    for row,a in Gamma[channel].items():
     if row>=48:leaves={P[source][row-47]:Q(1)}
     else:
      root=B[source].get(row)
      if root is None:continue
      node=catalog[root]
      if node[0]=='c':continue
      # Direct leaves are resolved in the FULL1302 source catalog. They
      # need not be among the eight composite nodes.
      if node[0]=='actual_C0_scalar_function':leaves={node[1]:Q(1)}
      else:assert root in leaf;leaves=leaf[root]
     for label,a0 in leaves.items():
      f=[w,str(c),str(a),str(a0)];terms.append([label,f]);b=(1,1)
      for x in f:b=mul(b,ex(x))
      bounds[label]=plus(bounds.get(label,(1,0)),b);labels.add(label)
  assert terms==rec['old_literal_second_terms_as_four_factor_products'] and native==rec['native128_terms_as_two_factor_products']
  assert [[j,str(a),str(b)] for j,(a,b) in sorted(bounds.items())]==rec['positive_old_symbol_coefficient_D_U']
  assert [[list(j),str(a),str(b)] for j,(a,b) in sorted(nbounds.items())]==rec['positive_native_coefficient_D_U']
  oldterms+=len(terms);nativeterms+=len(native)
 assert (oldterms,nativeterms,len(labels),maximum)==(192,126,10,68)
 assert d['committed_prior_bit_bound']==68 and d['temporary_prior_bit_bound']==205
 assert d['opaque_second_symbols_independent_with_literal_original_C0_labels'] and not d['signed_source_fold_AD_sensitivity_or_values_evaluated']
 assert count(d)==2572 and 2*count(d)+4000==9144<r['future_whole_numeric_bound']==11726
 assert r['future_RAW_bound']==2961532<32*1024*1024
 custody=read('work/a-zero-proof-v63/nonactive9-linear-ad-prepared-custody.json');wp=ROOT/custody['worker']
 assert sha(wp)==custody['worker_sha256']=='eb100fc7006534eb21fc169e73a7caf01b154291b908b7afa3573a070ecf4756'
 assert custody['gate_sha256']==sha(ROOT/rp) and custody['signed_sensitivity_or_values_not_executed'] and custody['safe_help_exit']==0
 worker=wp.read_text();assert 'old_second_function_values_not_queried' in worker and 'M1_pivot_G2_remaining_formal' in worker and 'count([g,p])<600' in worker
 inputs['prepared_formal_worker']={'path':custody['worker'],'sha256':sha(wp)}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent18-factor-linear-source-prior-and-custody-proof.json'
 e.write_text(json.dumps({'all18_weights_from_exact_b1_verified':True,'all8_nodes_match_full1302_catalog_and_are_linear':True,'all192_old_and126_native_literal_factor_records_verified':True,'all10_independent_honest_second_labels':sorted(labels),'all_positive_envelopes68bits_rebuilt':True,'all99_intrinsic_nativeS0_combinations_zero':True,'decoded_input_includes_three_metadata_fields':2572,'whole_live_count_with4000_scratch':9144,'no_signed_gradient_or_callback_values':True,
 'scope':'Only formal old-second/source fold. M1, pivot, G2 and lower terms remain formal; nativeS2/L1 and all honestsecond functionvalues remain unknown. The later W2 covariance theorem is not a premise for this restricted formal fold. The withdrawn composite-eight-catalog shortcut is not used.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_gate_sha256':sha(ROOT/rp),'formal_fold_worker_source_condition_closed':True,'no_signed_gradient_or_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'metadata_count':count([r,report]),'labels':sorted(labels)}))
if __name__=='__main__':main()
