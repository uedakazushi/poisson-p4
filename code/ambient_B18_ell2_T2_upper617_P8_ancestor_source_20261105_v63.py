#!/usr/bin/env python3
"""Stream literal two-P8 ancestor program; metadata only, no source values."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_P8_ANCESTOR_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def field(p,k):pins[p]=sha(ROOT/p);return fields(ROOT/p,[k])[k]
 cache='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';known79={i for i,v in field(cache,'first640')['all79_atom_values_full10']}
 physical='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';known164={i for i,v in field(physical,'all_requested164_root_values')}
 cp='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';pins[cp]=sha(ROOT/cp);text=(ROOT/cp).read_text();assert len(text.encode())<12*1024*1024;begin=text.index('[',text.index('"DAG_scalar_nodes"'))+1;dec=json.JSONDecoder()
 roots=field('reports/AMBIENT_B18_ELL2_T2_UPPER617_SUPPORT_PREFLIGHT_20261105_v63/factored617-lower-kernel-direct-third-support-and-source-frontier.json','additional_literal_C0_third_symbols');boundary=known79|known164;nodes={};hits={};unknown={};front=set(roots);layers=[]
 while front:
  assert len(nodes)<1500 and len(layers)<100
  wanted=set(front);front=set();pos=begin;found={}
  for i in range(max(wanted)+1):
   while text[pos].isspace() or text[pos]==',':pos+=1
   n,pos=dec.raw_decode(text,pos)
   if i in wanted:found[i]=n
  assert set(found)==wanted
  for i,n in found.items():
   if i in boundary and n[0]!='c':hits[i]=['C0_cached79' if i in known79 else 'physical164',n];continue
   nodes[i]=n;op=n[0]
   if op in ['a','m']:kids=n[1:3]
   elif op in ['n','i']:kids=n[1:2]
   elif op=='c':kids=[]
   else:unknown[i]=n;kids=[]
   for k in kids:
    assert isinstance(k,int) and k<i
    if k not in nodes and k not in hits:front.add(k)
  layers.append(sorted(wanted))
  (out/'phase.json').write_text(json.dumps({'phase':'literal_ancestor_metadata_only','layers_complete':len(layers),'retained_arithmetic_nodes':len(nodes),'boundaries':len(hits),'unknown_boundaries':len(unknown),'next_frontier':sorted(front)})+'\n')
 del text
 ep=out/'literal86-P8-composite-program-and-exact-boundary-catalog.json';ev={'root_nodes':roots,'selected_original_C0_arithmetic_nodes':sorted(nodes.items()),'exact_literal_cache_boundary_nodes':[[i,*v] for i,v in sorted(hits.items())],'uncached_or_unrecognized_opaque_boundaries':sorted(unknown.items()),'streamed_dependency_layers':layers,'all_parent_indices_strictly_lower':True,'source_inverse_nodes_remain_symbolic':sorted(i for i,n in nodes.items() if n[0]=='i'),'no_scalar_values_gradients_or_highestsource_jet_coefficients_evaluated':True,'scope':'ONLYliteralfunctionidentity;noavailabilityorfulluppercancellationclaim'};ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'literal_metadata_only':True,'arithmetic_nodes':len(nodes),'cached_boundaries':len(hits),'unknown_boundaries':len(unknown),'inverse_nodes':len(ev['source_inverse_nodes_remain_symbolic']),'layers':len(layers),'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'retained_numeric_upper':20*(len(nodes)+len(hits))+10000,'seconds':time.monotonic()-start};assert report['retained_numeric_upper']<60000;(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','arithmetic_nodes','cached_boundaries','unknown_boundaries','inverse_nodes','layers','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
