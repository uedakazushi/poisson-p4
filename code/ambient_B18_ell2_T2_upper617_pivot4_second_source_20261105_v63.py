#!/usr/bin/env python3
"""Literal ancestry for four OLD27 pivots, no second response evaluation."""
import argparse,json,resource,signal,time
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER617_PIVOT4_SECOND_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';allops=dict(field(mug,'old640_typed_closure_nodes'));B=[dict(v) for v in field(mug,'old640_source_columns')];P=field(H.PATHS['C0'],'honest_D0_relation_columns')
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63/pruned-known-Cfull224-input-and-positive-source-envelopes.json';piv=[j for j in field(cp,'pivot27_source_indices') if j<640];gamma=field(cp,'Gamma_two');rows={i for v in gamma for i,x in v};assert piv==[604,613,622,631]
 roots=[];tr=set();cr=set()
 for j in piv:
  col=[]
  for row in sorted(rows):
   if row<48 and row in B[j]:col.append([row,['typed',B[j][row]]]);tr.add(B[j][row])
   elif row==55:col.append([row,['C0',P[j][8]]]);cr.add(P[j][8])
  roots.append([j,col])
 del B,P
 keep=H.closure(allops,tr);ops={i:allops[i] for i in keep};del allops
 cr.update(n[1] for n in ops.values() if n[0]=='actual_C0_scalar_function')
 known79=set(field(mug,'old79_atom_ids'));raw='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';known164={i for i,v in field(raw,'all_requested164_root_values')};boundary=known79|known164
 sourcebase=count([ops,roots,gamma,piv,known79,known164]);assert sourcebase+20*1900+10000<150000
 src=H.PATHS['C0'];pins[src]=H.sha(ROOT/src);text=(ROOT/src).read_text();assert len(text.encode())<12*1024*1024;begin=text.index('[',text.index('"DAG_scalar_nodes"'))+1;dec=json.JSONDecoder();nodes={};hits={};unknown={};layers=[];front=set(cr)
 while front:
  assert len(nodes)<1500 and len(hits)<400 and len(layers)<100
  wanted=set(front);front=set();found={};pos=begin
  for i in range(max(wanted)+1):
   while text[pos].isspace() or text[pos]==',':pos+=1
   n,pos=dec.raw_decode(text,pos)
   if i in wanted:found[i]=n
  assert set(found)==wanted
  for i,n in found.items():
   if i in boundary and n[0]!='c':hits[i]=['source164' if i in known164 else 'C0_79_only',n];continue
   nodes[i]=n
   if n[0] in ['a','m']:kids=n[1:3]
   elif n[0] in ['n','i']:kids=n[1:2]
   elif n[0]=='c':kids=[]
   else:unknown[i]=n;kids=[]
   for k in kids:
    assert k<i
    if k not in nodes and k not in hits:front.add(k)
  layers.append(sorted(wanted));H.dump(out/'phase.json',{'phase':'OLD4pivot-literal-ancestry-only','typed_nodes':len(ops),'C0_nodes':len(nodes),'boundaries':len(hits),'unknown':len(unknown),'layers':len(layers),'no_secondresponse_or_gradient':True})
 del text
 ep=out/'literal-four-OLDpivot-Gamma-plane2-callbacks-and-frontier.json';ev={'original4_old_pivots':piv,'Gamma_two':gamma,'four_source_root_columns':roots,'typed_program':sorted(ops.items()),'literal_C0_program':sorted(nodes.items()),'cached_boundaries':[[i,*v] for i,v in sorted(hits.items())],'unrecognized_opaque_frontier':sorted(unknown.items()),'inverse_nodes':[i for i,n in nodes.items() if n[0]=='i'],'metadata_layers':layers,'K99old_order_proof_not_reused':True,'native_other23_pivots_full_L2S0_retained_separately':True,'future_required_pairs':['rr','rell'],'no_sourcevalue_or_opaquegradient_evaluated':True,'scope':'Onlysourceofouter-GammaXp2*a1term;actuala1/full617weightfoldandD*Xc2sensitivitymustcombinebeforehigherrequest.'};H.dump(ep,ev)
 owned=sourcebase+count([nodes,hits,unknown,layers,ev])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'metadata_only_no_secondresponses':True,'typed_nodes':len(ops),'C0_nodes':len(nodes),'boundaries':len(hits),'only79_boundaries':sum(v[0]=='C0_79_only' for v in hits.values()),'unknown_boundaries':len(unknown),'inverse_nodes':len(ev['inverse_nodes']),'source_owned_bound':owned,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','typed_nodes','C0_nodes','boundaries','only79_boundaries','unknown_boundaries','inverse_nodes','source_owned_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
