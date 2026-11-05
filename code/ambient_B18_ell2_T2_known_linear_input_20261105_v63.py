#!/usr/bin/env python3
"""Actual complete L48 source and bounded degree2 homogeneous closure input."""
import argparse,itertools,json,math,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 vp=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63_REVIEW/report.json');assert vp['all_passed']
 actual=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63/actual-two-invariants-A23-and-normalized-linear-graph.json')
 G=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json');B=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json')
 cols=[H.parse(G['actual_projected_G3linear_columns'][j]) for j in G['selected_columns']];labels=[['G3',j] for j in G['selected_columns']];assert len(cols)==25
 for k,j in enumerate(B['selected_bottom_columns']):
  v={27+i:F(c) for i,c in B['all642_bottom27_columns'][j]};a,b=map(F,actual['actual23_normalized_graph_values'][k]);v=add(v,{i:c for i,c in [(25,a),(26,b)] if c});cols.append(v);labels.append(['normalized_graph',k,j])
 piv={};receipts=[];maximum=0
 def guard(v):
  nonlocal maximum
  maximum=max([maximum]+[max(abs(c.numerator).bit_length(),c.denominator.bit_length()) for c in v.values()]);assert maximum<=2048;return v
 for j,orig in enumerate(cols):
  v=dict(orig);steps=[]
  while v and min(v) in piv:p=min(v);c=v[p];steps.append([p,str(c)]);v=guard(add(v,piv[p],-c))
  assert v,'all48 independent by exact first-layer source proof';p=min(v);c=v[p];piv[p]=guard(scale(v,1/c));receipts.append({'source_column':j,'pivot':p,'divide_by':str(c),'subtract':steps,'normalized':srec(piv[p])})
 assert len(piv)==48;free=[i for i in range(54) if i not in piv];allowed=[[i for i in range(9) if 9*t+i not in piv] for t in range(6)];assert sum(map(len,allowed))==6
 mons=[tuple(w.count(i) for i in range(9)) for w in itertools.combinations_with_replacement(range(9),2) if w!=(8,8)];assert len(mons)==44
 standards=[(t,m) for t in range(6) for m in mons if all(not m[v] or v in allowed[t] for v in range(9))];assert len(standards)<=21
 si={x:i for i,x in enumerate(standards)};memo={};nodes=[];depth={};height={}
 def dependency(t,m):
  key=t,m
  if key in memo:return memo[key]
  idx=len(nodes);memo[key]=idx;nodes.append(None)
  if m[8]>=2:n=['zero_ell2'];dep=0;D,U=1,0
  elif key in si:n=['standard',si[key]];dep=0;D,U=1,1
  else:
   v=next(v for v in range(9) if m[v] and 9*t+v in piv);row=piv[9*t+v];edges=[];bounds=[]
   for k,c in row.items():
    if k==9*t+v:continue
    target,w=divmod(k,9);q=list(m);q[v]-=1;q[w]+=1;q=tuple(q)
    # Position-over-lex order is multiplicative; all tails are smaller.
    assert (target,tuple(-x for x in q))>(t,tuple(-x for x in m))
    child=dependency(target,q);edges.append([child,str(-c)]);cd,cu=height[child];bounds.append((c.denominator*cd,abs(c.numerator)*cu))
   n=['rewrite',9*t+v,edges];dep=1+max([depth[j] for j,c in edges]+[0]);D=math.lcm(*(d for d,u in bounds)) if bounds else 1;U=sum(u*(D//d) for d,u in bounds)
  nodes[idx]=n;depth[idx]=dep;height[idx]=(D,U);return idx
 roots=[[dependency(t,m) for m in mons] for t in range(6)];assert all(n is not None for n in nodes)
 nfheight=max(max(d.bit_length(),u.bit_length()) for d,u in height.values());slots={'all48_source_and_normalized_linear_dense':2*48*54,
  'all264_linear_NF_degree2_dense':len(nodes)*len(standards),
  'all432_true_source_multiple_columns_dense':48*9*len(standards),
  'rank_basis_representations_and_factorizations_upper':3*len(standards)**2+48*9*len(standards),
  'source_dependency_edges':sum(len(n[2]) if n[0]=='rewrite' else 1 for n in nodes),'metadata_and_scratch':10000}
 total=sum(slots.values());assert total<60000 and nfheight<=2048
 ep=out/'actual-complete-L48-and-degree2-rewrite-source-input.json';H.dump(ep,{'ring':'k[zA0,...,zA6,rA,ell]/ell²',
  'all48_true_initial_linear_columns54':[srec(v) for v in cols],'actual_source_labels48':labels,
  'fixed48_independence_and_normalization_receipts':receipts,
  'normalized_linear_rewrite_rows':[[p,srec(v)] for p,v in sorted(piv.items())],
  'degree1_free_coordinates':free,'allowed_variables_per_target':allowed,
  'degree2_monomial_order44':[list(m) for m in mons],'degree2_linear_standard_module_monomials':[[t,list(m)] for t,m in standards],
  'degree2_monomial_rewrite_nodes':nodes,'all264_NF_root_nodes':roots,
  'rewrite_termination':'homogeneous position-over-lex; smaller targetindex leads, thenlex z0>...>r>ell; each linear tail isstrictlysmaller aftermultiplication',
  'all432_source_relations_to_reduce':'all48 normalizedtrueL rows timesall9variables, includingell evenifleadingell²vanishes. DoNOTassumeL isGroebner.',
  'source_completeness':'knownG3top25 plus23actualnormalizedgraph spanALL1060linearinitialrelations. Existenceofcoherentconstantunitframeisproved, canonicalM1other2notzero-filled.',
  'degree2_NF_positive_bound_bits':nfheight,'degree2_rewrite_max_depth':max(depth.values()),
  'whole_degree2_closure_slot_bound':total,'slot_breakdown':slots,
  'future_exact_known_cover_prefix':'h0+h1+h2_known=12+(standard2count-rankall432). If>20stop: thiscovercannotimproveacceptedactualupper20, notactualHF2lowerbound. Otherwisefreshdegree3overlap/slotgate beforeextension.',
  'actual_new_source_or_higher_T2_relations_not_queried':True,'degree2_source_rank_started':False,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':2048}})
 report={'id':STEM,'all_passed':True,'actual_known_linear_rank':48,'first_standard_count':len(free),'degree2_standard_count':len(standards),
  'degree2_NF_height_bound_bits':nfheight,'fixed_linear_max_actual_height_bits':maximum,'whole_degree2_closure_slot_bound':total,
  'degree2_rank_started':False,'new_source_or_actual_quadratic_relations_computed':False,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),
  'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','degree2_standard_count','degree2_NF_height_bound_bits','whole_degree2_closure_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
