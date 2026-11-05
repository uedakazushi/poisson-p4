#!/usr/bin/env python3
"""Positive/factored gate only for eighteen old linear-second sensitivities."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_LINEAR_AD_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 d=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_UPPER_SUPPORT_PREFLIGHT_20261105_v63/eighteen-formal-weights-and-coefficient-free-source-ancestor-inventory.json')
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
 allc=field(mug,'all741_zero_origin_source_combinations');sel=field(mu,'actual_selected_Y0_source_columns');C=[H.parse(allc[j]['original768']) for j in sel];del allc
 B=[dict(v) for v in field(mug,'old640_source_columns')];P=field(H.PATHS['C0'],'honest_D0_relation_columns');Gamma=list(map(H.parse,field('reports/AMBIENT_B18_ELL2_T2_G3_ACTIVE4_SOURCE8_INPUT_GATE_20261105_v63/actual-eight-source-root-folds-M1-operators-and-positive-AD-budget.json','Gamma56_two_channels')))
 ops=dict(d['typed_closure_nodes']);leaf={}
 for i,n in sorted(ops.items()):
  if n[0]=='actual_C0_scalar_function':v={n[1]:F(1)}
  elif n[0]=='c':v={}
  elif n[0]=='product':
   assert len(n[3])==1;v={j:F(n[1],n[2])*c for j,c in leaf[n[3][0]].items()}
  else:
   assert n[0]=='sum';v={}
   for j,a,b in n[1]:
    for k,c in leaf[j].items():v[k]=v.get(k,F(0))+F(a,b)*c
   v={k:c for k,c in v.items() if c}
  leaf[i]=v
 fixed=count([d,C,B,P,Gamma,leaf])+5000;assert fixed<60000
 H.dump(out/'phase.json',{'phase':'positive_linear_AD_input_only_before_any_signed_source_fold','fixed_source_ownership':fixed})
 prior=0
 def ck(v):
  nonlocal prior
  prior=max(prior,v[0].bit_length(),v[1].bit_length());return v
 def exact(v):v=F(v);return ck((v.denominator,abs(v.numerator)))
 def prod(a,b):return ck((a[0]*b[0],a[1]*b[1]))
 def summ(a,b):D=math.lcm(a[0],b[0]);return ck((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 rows=[]
 for r in d['eighteen_formal_contractions']:
  terms=[];native=[];bound={};nbound={}
  for k,ws in r['formal_weights_on99psi2']:
   j,ch=divmod(k,2)
   for source,c in C[j].items():
    if source>=640:
     factors=[ws,str(c)];native.append([ch,source-640,factors]);b=prod(exact(ws),exact(c));key=(ch,source-640);nbound[key]=summ(nbound.get(key,(1,0)),b);continue
    for row,a in Gamma[ch].items():
     if row<48:
      root=B[source].get(row)
      if root is None:continue
      assert root in ops or root not in {i for kind,i in r['possible_old_typed_or_literal_root_labels'] if kind==0}
      # Constant typed entries have zero second coefficient. The needed
      # nonconstant closure is exactly the support inventory's eight nodes.
      if root not in ops:continue
      leaves=leaf[root]
     else:leaves={P[source][row-48+1]:F(1)}
     for label,s in leaves.items():
      factors=[ws,str(c),str(a),str(s)];terms.append([label,factors]);b=(1,1)
      for x in factors:b=prod(b,exact(x))
      bound[label]=summ(bound.get(label,(1,0)),b)
  assert set(bound)<=set(r['uncancelled_structural_original_second_labels'])
  rows.append({'original418_source':r['original418_source'],'contraction':r['contraction'],'formal_weights_on99psi2':r['formal_weights_on99psi2'],
   'old_literal_second_terms_as_four_factor_products':terms,'native128_terms_as_two_factor_products':native,
   'positive_old_symbol_coefficient_D_U':[[j,str(v[0]),str(v[1])] for j,v in sorted(bound.items())],
   'positive_native_coefficient_D_U':[[list(j),str(v[0]),str(v[1])] for j,v in sorted(nbound.items())]})
 assert prior<=384,('source_positive_height',prior)
 data={'all18_factored_linear_second_sensitivity_inputs':rows,'typed_linear_root_provenance':d['typed_closure_nodes'],
  'opaque_second_symbols_independent_with_literal_original_C0_labels':True,'signed_source_fold_AD_sensitivity_or_values_evaluated':False,
  'native_source0_eligibility':'each selected99 source has separatelyzeroS0cN; every futureweightednativecombination inheritszero bylinearity. ActualS2/L1responses remain future callbacks, not evaluated.',
  'remaining_M1_pivot_G2_terms':'remain EXACT18psi source identity; this old-only sensitivity test does not discard anyM1, pivot, lower or G2 term.',
  'input_source_gauge_proof':'full56W2planezero source4b1c pendingindependentpeer; never setindividualM1/psi/mu tozero',
  'committed_prior_bit_bound':prior,'temporary_prior_bit_bound':3*prior+1}
 # Data plus raw/parsed input overlap, all180 possible old coefficients and
 # all actual114 native supports, one factored term stream and metadata.
 future=count(data)*2+2*(18*10+114)+6000;source=fixed+count(data)+4000;assert future<60000 and source<60000
 raw=2_000_000+future*(2*math.ceil(prior*math.log10(2))+40);assert raw<32*1024*1024
 data.update({'future_whole_numeric_bound':future,'future_RAW_bound':raw,'source_gate_ownership_bound':source})
 ep=out/'eighteen-factored-old-linear-second-and-native-source-inputs.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'signed_sensitivities_or_values_evaluated':False,'old_linear_AST_nodes':len(ops),'old_possible_symbols':len(set(j for r in rows for j,f in r['old_literal_second_terms_as_four_factor_products'])),
  'old_factored_term_count':sum(len(r['old_literal_second_terms_as_four_factor_products']) for r in rows),'native_factored_term_count':sum(len(r['native128_terms_as_two_factor_products']) for r in rows),
  'future_whole_numeric_bound':future,'future_RAW_bound':raw,'source_gate_ownership_bound':source,'positive_prior_bits':prior,'temporary_prior_bits':3*prior+1,
  'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','old_possible_symbols','old_factored_term_count','native_factored_term_count','future_whole_numeric_bound','future_RAW_bound','source_gate_ownership_bound','positive_prior_bits']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
