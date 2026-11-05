#!/usr/bin/env python3
"""Postcheck fixed active source bases and coherent full T2 frame transport."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv,mp,eye
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json',
 'constant':f'reports/{A}/actual-active-unit-bases-and26-frame-transport.json',
 'op':f'reports/{A}/actual-honest-T2-132-by1186-factored-source-operator.json',
 'old':'reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63/actual-honest-T2-459-by2070-factored-source-operator.json',
 'old_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json',
 'active':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
 'rank':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_RANK_20261105_v63/active-d3-rank-and-all-source-factorizations.json',
 'rank_peer':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_RANK_20261105_v63_REVIEW/report.json',
 'first':'reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json',
 'res':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
 'adj':'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json',
 'conormal':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/actual-V-conormal26-by56-unit-minor-and-source-factorization.json',
 'conormal_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def fixed_postcheck(n,cols,receipts,free,Pi):
 basis={}
 for r in receipts:
  p=dict(cols[r['source_column']])
  for i,c in r['subtract']:
   assert p.get(i)==Q(c);p=add(p,scale(basis[i],-Q(c)))
  i=r['pivot'];c=Q(r['divide_by']);assert c and min(p)==i and p[i]==c and i not in basis
  p=scale(p,1/c);assert enc(p)==r['normalized'];basis[i]=p
 assert free==[i for i in range(n) if i not in basis]
 assert all(mv(Pi,p)=={} for p in cols)
 assert all(Pi[i]=={j:Q(1)} for j,i in enumerate(free))
 # Fixed constraints plus the free unit values determine the quotient.
 for i,p in basis.items():assert not mv(Pi,p)
 return len(basis)
def determinant(M):
 a=[list(row) for row in M];n=len(a);out=Q(1)
 for j in range(n):
  k=next(i for i in range(j,n) if a[i][j])
  if k!=j:a[k],a[j]=a[j],a[k];out=-out
  c=a[j][j];out*=c
  for i in range(j+1,n):
   d=a[i][j]/c
   for q in range(j+1,n):a[i][q]-=d*a[j][q]
   a[i][j]=Q(0)
 return out
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'FIXED_ACTIVE_UNITS_AND_COHERENT26_FRAME','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in P if k=='primary' or k.endswith('_peer'))
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 C=D['constant'];S=D['active'];res=D['res'];d1=list(map(sp,S['active_d1_columns']));d2=list(map(sp,S['active_d2_columns']));d3=list(map(sp,S['active_d3_columns']))
 split1=S['signed_zero_wedge_decomposition']['1'];split2=S['signed_zero_wedge_decomposition']['2']
 sel56=[split2[j][1] for j in D['first']['fixed_original_minor_columns'] if split2[j][0]==[]]
 sel126=D['rank']['fixed_original_minor_columns'];assert len(set(sel56))==56 and len(sel126)==126
 assert C['active56_source_columns']==sel56 and C['active126_source_columns']==sel126
 a0=[];ap=[]
 for vec in res['source_a11_columns']:
  p={};pol=[]
  for row,coeff in vec:
   z,i,s=split1[row];assert not z
   pol.append([i,[[m,str(Q(c)*s)] for m,c in coeff]])
   c=sum((Q(c) for m,c in coeff if m==[0,0]),Q(0))*s
   if c:p[i]=c
  a0.append(p);ap.append(pol)
 assert ap==C['active11_source_columns_polynomial_r_ell']
 fixed1=a0+[d2[j] for j in sel56];fixed2=[{j:Q(1)} for j in sel56]+[d3[j] for j in sel126]
 assert list(map(enc,fixed1))==C['active_F1_fixed67_columns'] and list(map(enc,fixed2))==C['active_F2_fixed182_columns']
 free17=C['active_F1_free17'];free42=C['active_F2_free42'];Pi17=list(map(sp,C['active_Pi17_columns']));Pi42=list(map(sp,C['active_Pi42_columns']))
 assert fixed_postcheck(84,fixed1,C['active_F1_fixed67_receipts'],free17,Pi17)==67
 assert fixed_postcheck(224,fixed2,C['active_F2_fixed182_receipts'],free42,Pi42)==182
 assert len(free17)==17 and len(free42)==42
 assert all(not mv(Pi17,p) for p in d2) and all(not mv(Pi42,p) for p in d3)
 eta=res['physical_target3_basis_positions'];assert eta==[0,5,10]
 nonphys=[i for i in range(14) if i not in eta];ia=[]
 for k in free17:
  p={k:Q(1)}
  for i,e in enumerate(nonphys):p=add(p,scale(a0[i],-d1[k].get(e,Q(0))))
  assert not mv(d1,p) and mv(Pi17,p)=={len(ia):Q(1)};ia.append(p)
 inv={(tuple(z),i):(j,s) for j,(z,i,s) in enumerate(split1)}
 inc=[{inv[(),i][0]:c*inv[(),i][1] for i,c in p.items()} for p in ia]
 for z in S['zero_actions']:
  for i in eta:
   j,s=inv[(z,),i];inc.append({j:Q(s)})
 assert len(inc)==26 and list(map(enc,inc))==C['new_full_I1_origin26_columns']
 oldPi=list(map(sp,D['adj']['Pi26_source126_columns']));T=[mv(oldPi,p) for p in inc];Ti=list(map(sp,C['T0_inverse_columns']))
 assert list(map(enc,T))==C['T0_columns'] and mp(T,Ti)==mp(Ti,T)==eye(26)
 det=determinant([[T[j].get(i,Q(0)) for j in range(26)] for i in range(26)])
 assert str(det)==C['T0_determinant']=='-1'
 Mold=list(map(sp,D['conormal']['actual26_by56_columns']));Mnew=[mv(Ti,p) for p in Mold]
 assert list(map(enc,Mnew))==C['actual_Mnew0_columns']
 selected=D['conormal']['fixed_selected26_columns'];assert selected==C['same_selected26_V_columns']
 newdet=Q(D['conormal']['fixed26_minor_determinant'])/det
 assert str(newdet)==C['new_selected26_minor_determinant']=='-38773173828125/2404631929946112'
 factors=list(map(sp,D['conormal']['all56_columns_in_selected26_basis']))
 assert list(map(enc,factors))==C['all56_new_column_factorizations_in_same_selectedJ_basis']
 assert all(mv([Mnew[j] for j in selected],p)==v for p,v in zip(factors,Mnew))
 # Check all source-typed transport nodes and all768 unchanged cycles.
 E=D['op'];N=E['operator_nodes'];old=D['old'];ON=old['operator_nodes'];oldpath=P['old']
 def one(tag):
  matches=[i for i,n in enumerate(N) if n[0]==tag];assert len(matches)==1;(i,)=matches;return i
 def check(i,tag,*args):assert N[i]==[tag,*args]
 active=one('active5_periodic_Koszul_C2_resolution');assert N[active][1]==S['active_variables'] and N[active][2]['initial_ranks']==[14,84,224,364,434]
 b1=one('active_F1_unit_basis');assert N[b1][1]==[['a11',ap],['d2_columns',sel56],['Efree17_minus_a11_nonphysical_d1',free17]] and N[b1][2]['block_sizes']==[11,56,17]
 b2=one('active_F2_unit_basis');assert N[b2][1]==[['Eselected56',sel56],['d3_columns',sel126],['Efree42_minus_Eselected56_middle56_B1inv_d2',free42]] and N[b2][2]['block_sizes']==[56,126,42]
 comp238=[i for i in range(364) if i not in sel126]
 b3=one('active_F3_source_unit_clear');check(b3,'active_F3_source_unit_clear',sel126,comp238,'middle126_B2inv_d3',{'block_sizes':[126,238]})
 ag=one('active_reduced_complete_resolution');assert N[ag][1:5]==[active,b1,b2,b3] and N[ag][5]['first_ranks']==[3,17,42,238]
 tensor=one('tensor_active_resolution_with_Koszul_zero3');assert N[tensor][1:3]==[ag,[2,3,6]] and N[tensor][3]['first_ranks']==[3,26,102,418]
 i1=one('tensor_active_unit_chain_inclusion_degree1');assert N[i1][1]==tensor
 p2=one('tensor_active_unit_chain_projection_degree2');assert N[p2][1]==tensor and N[p2][2]['target_rank']==102
 op1=next(i for i,n in enumerate(ON) if n[0]=='last26_coordinates_B1inv');hp=next(i for i,n in enumerate(ON) if n[0]=='exact_periodic_Koszul_k_linear_homotopy')
 tt=next(i for i,n in enumerate(N) if n==['compose',['external_operator',oldpath,op1],i1])
 ti=one('formal_unit_inverse');assert N[ti][1]==tt and N[ti][2]['constant_determinant']=='-1'
 check(E['actual_Mnew_node'],'compose',ti,['external_operator',oldpath,old['actual_M26_by56_node']])
 check(E['actual_Nnew_node'],'compose',['external_operator',oldpath,old['actual_N56_by26_node']],tt)
 check(E['actual_Hfree_node'],'external_operator',oldpath,old['actual_Hfree56_by30_node'])
 h2=next(i for i,n in enumerate(N) if n==['compose',p2,['H_degree1_to2_external',oldpath,hp],i1])
 yy=E['all768_cycle_lift_top_nodes'];assert len(yy)==768
 for x,y in zip(old['source_X768_nodes'],yy):
  assert N[y][0:2]==['apply_vector',h2]
  check(N[y][2],'apply_vector',E['actual_Mnew_node'],['external_operator',oldpath,x])
 final=E['complete_T2_presentation_node'];check(final,'complete_honest_T2_presentation',{'shape':[132,1186],'first418_columns':['G3_of_tensor_resolution',tensor,'bottomzero30'],'last768_columns':['top102vectornodes',yy,'bottom30sourcevectors',oldpath,old['all768_actual_cycle_lift_bottom_nodes']]})
 # The known active contractible pairs tensor to these full unit counts.
 paircounts=[11,3*11+56,3*11+3*56+126,11+3*56+3*126]
 assert paircounts==[11,89,327,557]
 assert [3,17+3*3,42+3*17+3*3,238+3*42+3*17+3]==[3,26,102,418]
 assert E['source_proof']['rank418_not_claimed_minimal_beta3']
 ep=o/'independent-active-unit-and-132-by1186-source-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all_fixed67_and182_source_receipts_and_quotients_verified':True,'all364_active_d3_annihilations_checked':True,'common_physical_eta3_and_new26_inclusion_verified':True,'T0_two_sided_inverse_and_independent_determinant':str(det),'new_selected26_V_minor':str(newdet),'all56_factorizations_transferred_exactly':True,'full_tensor_unit_pair_counts':paircounts,'no_numeric_d4_rank_inferred':True,'complete_C2_initial_ranks':[3,26,102,418],'all768_original_source_cycles_and_typed_lifts_retained':True,'complete_honest_T2_source_shape':[132,1186],'new_physical_coefficient_or_T2_length_computed':False})
 report={'id':STEM,'all_passed':True,'complete_actual_T2_typed132_by1186_accepted':True,'active_unit_and_tensor557_transport_accepted':True,'coherent26_frame_det_minus1_and_new_minor_accepted':True,'new_source_coefficient_or_T2_length_or_H0_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
