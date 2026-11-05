#!/usr/bin/env python3
"""Literal degree2 matrix + singleton-unit graph profile; no general rank."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_G2_STRAND_SOURCE_PROFILE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 pp='reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json';var=field(pp,'variable_order6');G=field(pp,'actual_linear_G2_active102_by42_columns');assert var==['zA0','zA1','zA4','zA5','rA','ell'] and len(G)==42
 mon=[(i,j) for i in range(6) for j in range(i,6) if (i,j)!=(5,5)];ix={m:i for i,m in enumerate(mon)};assert len(mon)==20
 A=[];byrow=[{} for _ in range(340)]
 for j,col in enumerate(G):
  for v in range(6):
   q={}
   for slot,x in col:
    row,u=divmod(slot,6);pair=tuple(sorted((u,v)))
    if pair==(5,5):continue
    k=20*row+ix[pair];assert k not in q and F(x);q[k]=x
   A.append(q)
 for j,col in enumerate(A):
  for i,x in col.items():byrow[i][j]=x
 # Graph-only leaf deletion. Remaining coefficients never change.
 liveR=set(range(340));liveC=set(range(252));rs=[set(v) for v in byrow];cs=[set(v) for v in A];piv=[]
 while True:
  cand=next(((i,next(iter(rs[i])),'row_singleton') for i in sorted(liveR) if len(rs[i])==1),None)
  if cand is None:
   cand=next(((next(iter(cs[j])),j,'column_singleton') for j in sorted(liveC) if len(cs[j])==1),None)
  if cand is None:break
  i,j,kind=cand;x=A[j][i];piv.append({'row340':i,'column252':j,'literal_nonzero_unit':x,'singleton_kind':kind})
  for k in tuple(rs[i]):cs[k].discard(i)
  for k in tuple(cs[j]):rs[k].discard(j)
  rs[i].clear();cs[j].clear();liveR.remove(i);liveC.remove(j)
 R=sorted(i for i in liveR if rs[i]);C=sorted(j for j in liveC if cs[j]);res=[[j,[[i,A[j][i]] for i in sorted(cs[j])]] for j in C];assert all(i in R for j,col in res for i,x in col)
 # Integer row envelopes only: no determinant or residual rank is evaluated.
 clear=[]
 for i in R:
  D=1
  for j in rs[i]:D=math.lcm(D,F(A[j][i]).denominator)
  U=sum(abs(F(A[j][i]).numerator)*(D//F(A[j][i]).denominator) for j in rs[i]);clear.append([i,str(D),str(U),U.bit_length()])
 upper=190-len(piv);assert upper>=0
 dbits=sum(sorted((v[3] for v in clear),reverse=True)[:upper]) if upper else 0
 src=read('reports/AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63/actual-full-G3linear-active-fixed-rank-and-kernel.json');assert src['active_rank']==62 and src['full_G3linear_rank']==242
 sourceproof={'active_domain':'F2_active42 times6linearcoefficients =252','active_codomain':'F1_active17 times20ordinaryquadratics(ell²omitted)=340','rank_upper190':'exactactiveallorderG2G3=0 andbothorigins0 =>linearproduct0; acceptedG3linear62independentintheSAME42basis',
 'full918_reduction_pending_independent_source_proof':'linearcomplexistheexactzero3Koszultensorofactivecomplexintheacceptedunittransport. IfactiveH2degree3=0,zero3KresolvesregularnewvariablesandfullH2degree3=0. Thisisaconditionalreduction,notagradedexactnessassumption.',
 'singleton_minor_proof':'recursivelyexpandselectedminorbythecertifiedsingletonroworcolumn;det=sign*productliteralnonzeropivots;nofillandnolargeevaluatedproduct','residual_rank_NOT_queried':True}
 ep=out/'literal-active-G2-linear-strand-matrix-singleton-units-and-residual-source.json';data={'active_variable_order6':var,'ordinary_quadratic_pairs20':mon,'matrix_shape':[340,252],'literal_matrix_columns':[sorted(v.items()) for v in A],'original_linear_G2_columns102_by42':G,'fixed_singleton_unit_pivots':piv,'remaining_nonzero_row_indices':R,'remaining_nonzero_column_indices':C,'remaining_literal_matrix_columns':res,'remaining_row_integer_D_L1_bounds':clear,'known_total_rank_upper190':True,'remaining_rank_upper':upper,'all_possible_remaining_rank_minor_integer_bit_bound_from_rowL1':dbits,'source_proof':sourceproof,'general_rank_or_kernel_computed':False,'new_physical_or_source_coefficients':False};H.dump(ep,data)
 owned=count([G,A,byrow,piv,R,C,res,clear,data])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'SOURCE_literal_pattern_and_singleton_units_only':True,'general_rank_query_run':False,'matrix_nonzeros':sum(map(len,A)),'singleton_unit_count':len(piv),'remaining_nonzero_shape':[len(R),len(C)],'remaining_rank_upper':upper,'residual_minor_integer_bit_bound':dbits,'source_owned_bound':owned,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','matrix_nonzeros','singleton_unit_count','remaining_nonzero_shape','remaining_rank_upper','residual_minor_integer_bit_bound','source_owned_bound']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
