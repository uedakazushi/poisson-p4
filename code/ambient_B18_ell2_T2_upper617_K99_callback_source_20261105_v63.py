#!/usr/bin/env python3
"""99 Kfull callback ancestry/source/count only; no formal weights or values."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63'

def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);pins={}
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 gate='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 allc=field(gate,'all741_zero_origin_source_combinations')
 yp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json'
 selected=field(yp,'actual_selected_Y0_source_columns');combs=[allc[i]['original768'] for i in selected];del allc
 assert len(combs)==99
 qpath='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json'
 gamma=field(qpath,'Gamma_two');rows={i for row in gamma for i,x in row};assert rows=={1,2,4,17,18,21,22,28,55}
 allops=dict(field(gate,'old640_typed_closure_nodes'));B=[dict(v) for v in field(gate,'old640_source_columns')]
 cp=H.PATHS['C0'];P=field(cp,'honest_D0_relation_columns')
 old=sorted({int(j) for c in combs for j,x in c if int(j)<640})
 roots=[];typedroots=set();c0roots=set()
 for j in old:
  v=[]
  for row in sorted(rows):
   if row<48 and row in B[j]:v.append([row,['typed',B[j][row]]]);typedroots.add(B[j][row])
   elif row==55:v.append([row,['C0',P[j][8]]]);c0roots.add(P[j][8])
  roots.append([j,v])
 del B,P
 closure=H.closure(allops,typedroots);typed={i:allops[i] for i in closure};del allops
 c0roots.update(n[1] for n in typed.values() if n[0]=='actual_C0_scalar_function')
 known79=set(field(gate,'old79_atom_ids'))
 raw='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json'
 known164={i for i,v in field(raw,'all_requested164_root_values')};boundary=known79|known164
 # At most1500 literal arithmetic nodes and400 boundaries are retained;
 # metadata/source IDs only, not any cached point/first/second response.
 basecount=count([combs,selected,gamma,roots,typed,c0roots,known79,known164,pins])
 assert basecount+20*(1500+400)+10000<150000
 pins[cp]=H.sha(ROOT/cp);text=(ROOT/cp).read_text();assert len(text.encode())<12*1024*1024
 begin=text.index('[',text.index('"DAG_scalar_nodes"'))+1;decoder=json.JSONDecoder()
 nodes={};hits={};unknown={};front=set(c0roots);layers=[]
 while front:
  assert len(nodes)<1500 and len(hits)<400 and len(layers)<100
  wanted=set(front);front=set();pos=begin;found={}
  for i in range(max(wanted)+1):
   while text[pos].isspace() or text[pos]==',':pos+=1
   n,pos=decoder.raw_decode(text,pos)
   if i in wanted:found[i]=n
  assert set(found)==wanted
  for i,n in found.items():
   if i in boundary and n[0]!='c':hits[i]=['source164' if i in known164 else 'C0_79_only',n];continue
   nodes[i]=n;op=n[0]
   if op in ('a','m'):kids=n[1:3]
   elif op in ('n','i'):kids=n[1:2]
   elif op=='c':kids=[]
   else:unknown[i]=n;kids=[]
   for k in kids:
    assert isinstance(k,int) and k<i
    if k not in nodes and k not in hits:front.add(k)
  layers.append(sorted(wanted))
  H.dump(out/'phase.json',{'phase':'source_metadata_ancestry_no_gradient_or_secondresponse','layers_complete':len(layers),'typed_nodes':len(typed),'literal_nodes':len(nodes),'cache_boundaries':len(hits),'unknown_boundaries':len(unknown),'next_frontier':sorted(front)})
 del text
 ep=out/'actual99-selected-source-K-callback-roots-and-literal-ancestry.json'
 data={'selected99_kernel741_indices':selected,'selected99_original768_constant_combinations':combs,'Gamma_two':gamma,'needed_original_old_columns':roots,'typed_program':sorted(typed.items()),'literal_C0_program':sorted(nodes.items()),'literal_cache_boundaries':[[i,*v] for i,v in sorted(hits.items())],'unrecognized_boundaries_exact_frontier':sorted(unknown.items()),'literal_inverse_nodes_retained':[i for i,n in sorted(nodes.items()) if n[0]=='i'],'metadata_dependency_layers':layers,'Cfull_source':'reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/report.json','Cfull_source_sha256':H.sha(ROOT/'reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/report.json'),'formula':'Kfull_j_rr=Gamma Xstar_j_rr+Cfull_r Xstar_j_r;Kfull_j_rell=Gamma Xstar_j_rell+Cfull_r Xstar_j_ell+Cfull_ell Xstar_j_r','native_required_source_condition':'S0*cN=0 separately for all99; THENrr=L0S2rr+LrSr,re=L0S2re+LrSe+LeSr. Rebindbeforeanyresponse,neveruseindividualpivoteffectiveformula.','no_signed_rootweights_formal_AD_or_secondresponses':True,'scope':'ONLYliteralcallbacksourceancestry;noh2cacheavailabilityoractualKvalueclaim.'}
 H.dump(ep,data);owned=basecount+count([nodes,hits,unknown,layers,data])+10000;assert owned<150000
 report={'id':STEM,'all_passed':True,'metadata_only_no_response_values':True,'selected99':99,'original_old_columns':len(old),'typed_nodes':len(typed),'C0_arithmetic_nodes':len(nodes),'boundaries':len(hits),'only79_boundaries':sum(v[0]=='C0_79_only' for v in hits.values()),'unknown_boundaries':len(unknown),'inverse_nodes':len(data['literal_inverse_nodes_retained']),'future_formal_sensitivity_or_K_value_GO_not_implied':True,'source_owned_numeric_bound':owned,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','original_old_columns','typed_nodes','C0_arithmetic_nodes','boundaries','only79_boundaries','unknown_boundaries','inverse_nodes','source_owned_numeric_bound','seconds']}))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
