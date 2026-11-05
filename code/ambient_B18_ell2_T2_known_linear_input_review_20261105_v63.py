#!/usr/bin/env python3
"""Independent true-linear-source and finite degree2 rewrite-DAG audit."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,math,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual-complete-L48-and-degree2-rewrite-source-input.json',
'graph':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63/actual-two-invariants-A23-and-normalized-linear-graph.json',
'graph_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63_REVIEW/report.json',
'HF_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63_REVIEW/report.json',
'G3':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json',
'bottom':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['graph_peer']['all_passed'] and D['HF_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['gate'];cols=[sp(D['G3']['actual_projected_G3linear_columns'][j]) for j in D['G3']['selected_columns']];labels=[['G3',j] for j in D['G3']['selected_columns']]
 for k,j in enumerate(D['bottom']['selected_bottom_columns']):
  p={27+i:Q(c) for i,c in D['bottom']['all642_bottom27_columns'][j]};a,b=map(Q,D['graph']['actual23_normalized_graph_values'][k]);p=add(p,{i:c for i,c in [(25,a),(26,b)] if c});cols.append(p);labels.append(['normalized_graph',k,j])
 assert len(cols)==48 and list(map(enc,cols))==E['all48_true_initial_linear_columns54'] and labels==E['actual_source_labels48']
 piv={}
 for r in E['fixed48_independence_and_normalization_receipts']:
  p=dict(cols[r['source_column']])
  for k,c in r['subtract']:p=add(p,scale(piv[k],-Q(c)))
  k=r['pivot'];assert p and min(p)==k and k not in piv and p[k]==Q(r['divide_by']);piv[k]=scale(p,1/Q(r['divide_by']));assert enc(piv[k])==r['normalized']
 assert len(piv)==48 and [[k,enc(p)] for k,p in sorted(piv.items())]==E['normalized_linear_rewrite_rows'];free=[i for i in range(54) if i not in piv];assert free==E['degree1_free_coordinates'] and len(free)==6
 allowed=[[v for v in range(9) if 9*i+v not in piv] for i in range(6)];assert allowed==E['allowed_variables_per_target']
 mons=[tuple(w.count(i) for i in range(9)) for w in itertools.combinations_with_replacement(range(9),2) if w!=(8,8)];assert list(map(list,mons))==E['degree2_monomial_order44']
 standards=[(i,m) for i in range(6) for m in mons if all(not e or v in allowed[i] for v,e in enumerate(m))];assert [[i,list(m)] for i,m in standards]==E['degree2_linear_standard_module_monomials'] and len(standards)==6
 nodes=E['degree2_monomial_rewrite_nodes'];labels_nodes={};bounds={};depth={};standardindex={p:i for i,p in enumerate(standards)}
 def verify(node,target,m):
  label=target,m
  if node in labels_nodes:assert labels_nodes[node]==label;return
  labels_nodes[node]=label;n=nodes[node]
  if m[8]>=2:assert n==['zero_ell2'];bounds[node]=(1,0);depth[node]=0;return
  if label in standardindex:assert n==['standard',standardindex[label]];bounds[node]=(1,1);depth[node]=0;return
  v=next(v for v,e in enumerate(m) if e and 9*target+v in piv);pivot=9*target+v;assert n[0]=='rewrite' and n[1]==pivot;tail=[(k,c) for k,c in piv[pivot].items() if k!=pivot];assert len(n[2])==len(tail);terms=[];childdepth=[]
  for(child,q),(k,c) in zip(n[2],tail):
   i,w=divmod(k,9);new=list(m);new[v]-=1;new[w]+=1;new=tuple(new);assert Q(q)==-c and (i,tuple(-x for x in new))>(target,tuple(-x for x in m));verify(child,i,new);d,u=bounds[child];terms.append((c.denominator*d,abs(c.numerator)*u));childdepth.append(depth[child])
  d=math.lcm(*(a for a,b in terms)) if terms else 1;u=sum(b*(d//a) for a,b in terms);bounds[node]=(d,u);depth[node]=1+max(childdepth+[0])
 roots=E['all264_NF_root_nodes'];assert len(roots)==6 and all(len(row)==44 for row in roots)
 for i,row in enumerate(roots):
  for node,m in zip(row,mons):verify(node,i,m)
 assert len(labels_nodes)==len(nodes)==270 and len(set(labels_nodes.values()))==270
 h=max(max(d.bit_length(),u.bit_length()) for d,u in bounds.values());assert h==139==E['degree2_NF_positive_bound_bits'] and max(depth.values())==E['degree2_rewrite_max_depth']
 assert sum(E['slot_breakdown'].values())==22375==E['whole_degree2_closure_slot_bound']<150000 and E['degree2_source_rank_started']==False
 ep=o/'independent-trueL48-and-all264-rewrite-source-DAG.json';dump(o,ep.name,{'inputs':inputs,'all48_true_linear_generators_and_fixed_independence_receipts_verified':True,'actual_first_standard_count':6,'all264_degree2_roots_and270_dependency_nodes_verified':True,'strict_homogeneous_position_over_lex_termination':True,'degree2_candidate_standard_count':6,'NF_positive_height_bound_bits':139,'whole_degree2_closure_slot_bound':22375,
'source_representative_scope':'Graphuppercomponentsarecanonicalrepresentativesmodtheknown25G3linearrelations. EqualityofPi2classesprovestheyaretrueinitialrelationsinthereviewedexistingcoherentconstantunitframe; itdoesnotprovideunknownhigherliftsorassigncanonicalM1values0.',
'next_complete_degree2_step':'ReduceALL48trueLrows timesALL9variables(432columns), includingellwhenleadingell²vanishes. Thisspansentiredegree2linearinitialsubmodule,sofixedrankandcompletefactorizationgiveEXACTdegree2oftheknownlinearcover. NoGroebnerassumptionoroverlapdiscardingispermitted.',
'safe_HF_scope':'KnownlinearcovermapsontoactualgrT2; itsdegree2dimensionisanUPPERboundforactualHF2,neverlowerbound. Ifcoverprefixexceeds20itmerelycannotimprovetheacceptedactualupper20. NoactualhigherHF/totalT2lengthcomputed.',
'no_new_source_implicit_gauge_or_actual_quadratic_relation_values':True})
 report={'id':STEM,'all_passed':True,'actualL48_source_and_complete_degree2_rewrite_input_accepted':True,'actual_degree2_relation_rank_or_HF2_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
