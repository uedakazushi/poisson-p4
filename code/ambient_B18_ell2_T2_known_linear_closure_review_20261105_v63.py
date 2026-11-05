#!/usr/bin/env python3
"""Independent432 reductions and global linear-module Buchberger certificate."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63';STEM=A+'_REVIEW';B='AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63'
P={'primary':f'reports/{A}/report.json','gb':f'reports/{A}/all-critical-pairs-and-known-cover-stop-proof.json',
'degree2_report':f'reports/{B}/report.json','degree2':f'reports/{B}/actual432-overlap-rank-and-known-linear-cover-degree2.json',
'failure':f'reports/{B}/FIRST_REPORT_PATH_SHADOW.json','executed_snapshot':f'reports/{B}/FIRST_REPORT_PATH_SHADOW_CODE.py',
'input':'reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json',
'input_peer':'reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63_REVIEW/report.json',
'length_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k!='executed_snapshot'};assert all(D[k]['all_passed'] for k in ['primary','degree2_report','input_peer','length_peer'])
 for key,stem in [('primary',A),('degree2_report',B)]:
  for p,h in D[key]['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in D[key]['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/p)==h
 assert inputs['executed_snapshot']['sha256']==D['degree2_report']['executed_math_code_sha256'] and inputs['degree2']['sha256']==D['failure']['completed_math_evidence_sha256'];assert D['failure']['report_only_repair'] and not D['failure']['arithmetic_rerun']
 G=D['input'];E=D['degree2'];nodes=G['degree2_monomial_rewrite_nodes'];memo={}
 def nf(i):
  if i in memo:return memo[i]
  n=nodes[i]
  if n[0]=='standard':p={n[1]:Q(1)}
  elif n[0]=='zero_ell2':p={}
  else:
   assert n[0]=='rewrite';p=add(*(scale(nf(j),Q(c)) for j,c in n[2]))
  memo[i]=p;return p
 roots=G['all264_NF_root_nodes']
 for row in roots:
  for i in row:nf(i)
 assert [[i,enc(p)] for i,p in sorted(memo.items())]==E['all264_linear_NF_columns'] and len(memo)==270
 mons=list(map(tuple,G['degree2_monomial_order44']));mi={m:i for i,m in enumerate(mons)};rows={k:sp(p) for k,p in G['normalized_linear_rewrite_rows']};multiple={};labels=[];remainders=[]
 for pivot,row in rows.items():
  for var in range(9):
   p={}
   for k,c in row.items():
    target,w=divmod(k,9);m=[0]*9;m[var]+=1;m[w]+=1
    if m[8]<2:p=add(p,scale(memo[roots[target][mi[tuple(m)]]],c))
   assert not p;multiple[pivot,var]=len(labels);labels.append([pivot,var]);remainders.append(enc(p))
 assert len(labels)==432 and labels==E['all432_true_source_multiple_labels'] and remainders==E['all432_remainder_columns']
 assert E['rank']==0 and E['fixed_rank_receipts']==[] and E['fixed_original_minor_source_columns']==[] and all(not p for p in E['all432_original_selected_source_factorizations'])
 assert E['known_degree2_quotient_projection_columns']==[[[i,'1']] for i in range(6)] and E['degree2_quotient_free_standard_indices']==list(range(6))
 critical=[]
 for p,q in itertools.combinations(rows,2):
  if p//9==q//9:critical.append({'linear_leads':[p,q],'S_pair_zero_NF_source_multiple_columns':[multiple[p,q%9],multiple[q,p%9]],'combination':[1,-1]})
 assert len(critical)==171 and critical==D['gb']['all_linear_linear_critical_pairs']
 ell=[]
 for p in rows:
  if p%9==8:ell.append({'linear_lead':p,'kind':'leading_ell_overlap_degree2','zero_NF_source_column':multiple[p,8]})
  else:ell.append({'linear_lead':p,'kind':'lcm_degree3_all_terms_have_ell_squared','source':'ell²*L -leadvar*e_target*ell²'})
 assert len(ell)==48 and ell==D['gb']['all_linear_ell_squared_critical_pairs']
 allowed=G['allowed_variables_per_target'];assert allowed==[[],[],[],[7,8],[7,8],[7,8]]==D['gb']['allowed_variables_per_target']
 hf=[]
 for degree in range(4):
  count=0
  for target in range(6):
   for word in itertools.combinations_with_replacement(range(9),degree):
    if word.count(8)<2 and all(v in allowed[target] for v in word):count+=1
  hf.append(count)
 assert hf==[6,6,6,6]==D['gb']['known_cover_HF_through3'] and sum(hf)==24>20
 ep=o/'independent432-zero-NFs-219-critical-pairs-and-cover-stop.json';dump(o,ep.name,{'inputs':inputs,'all270_NF_nodes_and432_true_multiples_reconstructed_zero':True,'degree2_known_relation_rank':0,'all171_linear_linear_Spairs_and48_baseell2_pairs_complete':True,
'global_GB_proof':'Take48normalizedlinearrowsand6ell²e_targetgeneratorsinthefullpolynomialmodule. Same-targetlinearS-pairsarethe171listeddegree2differencesandreduce0. Ell²S-pairsare48listedcases:leadingellusesverifiedell*L NF0; othercasesareentirelyell²multiplesandreduce0. DifferenttargetcomponentpairsdonotformmoduleS-pairs. Thehomogeneousmultiplicativewellorderwasindependentlycheckedbytheinputpeer, soBuchbergercriterionprovesglobalGB.',
'known_cover_Hilbert_series':'6/(1-t)','known_cover_HF_through3':hf,'known_cover_prefix_length':24,
'safe_I2_and_stop_scope':'IA,IB arecontainedinm, soI²=(IAIB)²iscontainedinm4. AddingONLYthisannihilatorcannotchangecoverdegrees0..3; its24prefixalreadyexceedsacceptedactualupper20. Stopthiscoarsecover. ActualgrT2mayhaveadditionalquadraticandhigherinitialrelations, soNOactualHF2/3lowerboundorequalityisclaimed.',
'metadata_failure_scope':'Executedcodeandcompletedc130evidencehasharepreserved; failurewasthereportpathshadowafterevidence. Thispeerreconstructsallmathematicalcontentsindependentlyanddoesnotreruntheproduceroritsranksearch.',
'new_source_or_actual_quadratic_coefficients_computed':False})
 report={'id':STEM,'all_passed':True,'all432_zero_reductions_and_global_knownlinear_GB_accepted':True,'known_cover_HF_through3':hf,'known_cover_prefix_length':24,'actual_T2_higher_HF_length_or_H0_not_claimed':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
