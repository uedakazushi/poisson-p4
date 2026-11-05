#!/usr/bin/env python3
"""Independent typed-source and exact-homotopy review; no coefficient query."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,math,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,add,scale,mv,mp,eye,enc
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json',
 'op':f'reports/{A}/actual-honest-T2-459-by2070-factored-source-operator.json',
 'old':'reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63/phase1-actual-X768-source.json',
 'old_report':'reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63/report.json',
 'old_peer':'reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63_REVIEW/report.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'chain_peer':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'SyzQ':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'SyzQ_peer':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW/report.json',
 'conormal':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/actual-V-conormal26-by56-unit-minor-and-source-factorization.json',
 'conormal_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'res':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
 'res_peer':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json',
 'betti':'reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json',
 'betti_peer':'reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 'adj':'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json',
 'source_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'global_tensor_peer':'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
 'global_conductor_peer':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW/report.json',
 'closed348_peer':'reports/AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63_REVIEW/report.json',
 'raw33_peer':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):
 (o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in o.iterdir())<32*1024*1024
def put(d,k,c):
 c=d.get(k,Q(0))+c
 if c:d[k]=c
 else:d.pop(k,None)
def plus(*ps):
 out={}
 for p in ps:
  for k,c in p.items():put(out,k,c)
 return out
class HomotopyFixtures:
 """Own polynomial total-complex implementation for two actual actions."""
 def __init__(self,A,E):self.A=A;self.E=E
 def power(self,vec,a,b):
  for _ in range(a):vec=mv(self.A[0],vec)
  for _ in range(b):vec=mv(self.A[1],vec)
  return vec
 def translation(self,f,sign):
  out={}
  for (p,w,m,j),c in f.items():
   a,b,e=m
   for u in range(a+1):
    for v in range(b+1):
     factor=c*math.comb(a,u)*math.comb(b,v)*sign**(a+b-u-v)
     for i,d in self.power({j:Q(1)},a-u,b-v).items():put(out,(p,w,(u,v,e),i),factor*d)
  return out
 def differential(self,f):
  out={}
  for (p,w,m,j),c in f.items():
   if p:
    if not m[2]:put(out,(p-1,w,(m[0],m[1],1),j),c)
    for i,d in self.E[j].items():put(out,(p-1,w,m,i),c*(-1)**p*d)
   for t,k in enumerate(w):
    w2=w[:t]+w[t+1:];sgn=(-1)**(p+t);m2=list(m);m2[k]+=1
    put(out,(p,w2,tuple(m2),j),c*sgn)
    for i,d in self.A[k][j].items():put(out,(p,w2,m,i),-c*sgn*d)
  return out
 def homotopy(self,f):
  shifted=self.translation(f,1);hk={}
  for (p,w,m,j),c in shifted.items():
   for k in range(min(w) if w else 2):
    if m[k] and not any(m[:k]):
     n=list(m);n[k]-=1;put(hk,(p,(k,)+w,tuple(n),j),c*(-1)**p)
  out=self.translation(hk,-1)
  for (p,w,m,j),c in f.items():
   if not w and m[2]:
    for i,d in self.power({j:Q(1)},m[0],m[1]).items():put(out,(p+1,(),(0,0,0),i),c*d)
  return out
 def projection(self,f):
  out={}
  for (p,w,m,j),c in f.items():
   if p or w:continue
   v=self.power({j:Q(1)},m[0],m[1])
   if m[2]:v=mv(self.E,v)
   for i,d in v.items():put(out,(0,(),(0,0,0),i),c*d)
  return out
 def run(self):
  count=0;mons=[(a,b,e) for a in range(3) for b in range(3-a) for e in range(2)]
  for p in range(3):
   for q in range(3-p):
    for w in itertools.combinations(range(2),q):
     for m in mons:
      for j in range(14):
       f={(p,w,m,j):Q(1)}
       actual=plus(self.differential(self.homotopy(f)),self.homotopy(self.differential(f)))
       expected=plus(f,{k:-v for k,v in self.projection(f).items()})
       assert actual==expected,('TOTAL_HOMOTOPY_SIGN',p,w,m,j)
       count+=1
  return count
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'TYPED_ALLORDER_SOURCE_AND_HOMOTOPY_REVIEW','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in P if k=='primary' or k.endswith('_peer'))
 R=D['primary']
 for p,h in R['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in R['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 assert sha(ROOT/'code/ambient_B18_ell2_honest_T2_source_operator_20261105_v63.py')==R['code_sha256']
 assert sha(ROOT/'research/notes'/f'{A}.md')==R['note_sha256']
 old=D['old'];O=old['operator_nodes'];manifest=old['source_manifest'];E=D['op'];N=E['operator_nodes']
 assert len(O)==8107 and manifest['actual_X_shape']==[56,768]
 # Verify source ancestry directly, rather than using the old ell0 values.
 original_sources={k:old['source_hashes'][k] for k in ['B','B_peer','P','SyzQ','SyzQ_peer']}
 for v in original_sources.values():assert sha(ROOT/v['path'])==v['sha256']
 assert E['source_operator_proof']['original_full_R0_callback_artifacts']==original_sources
 C0=D['C0'];B=D['chain']
 for i,n in enumerate(O):
  assert N[i]==[n[0].replace('_restrict_ell0','_restrict_ell_squared')]+n[1:]
  if n[0]=='actual_C0_scalar_function_restrict_ell0':assert 0<=n[1]<len(C0['DAG_scalar_nodes'])
  elif n[0]=='actual_chain_scalar_function_restrict_ell0':assert 0<=n[1]<len(B['operator_scalar_nodes'])
  elif n[0]=='actual_Res4_SyzQ_scalar_entry_restrict_ell0':assert 0<=n[1]<48 and 0<=n[2]<128
  else:assert n==['zero_function']
 X=E['source_X768_nodes'];assert len(X)==768
 for j,(node,col) in enumerate(zip(X,manifest['source_X_sparse_columns'])):
  assert N[node]==['actual_X_column_restrict_ell_squared',j,col]
  actual={i:O[a] for i,a in col}
  if j<640:
   expected={i:['actual_chain_scalar_function_restrict_ell0',v] for i,v in B['B48_by640_sparse_columns'][j]}
   expected.update({47+i:['actual_C0_scalar_function_restrict_ell0',v] for i,v in enumerate(C0['honest_D0_relation_columns'][j]) if i and v})
  else:expected={i:['actual_Res4_SyzQ_scalar_entry_restrict_ell0',i,j-640] for i in range(48)}
  assert actual==expected,('ORIGINAL_SOURCE_COLUMN',j)
 assert B['all640_UP_equals_QB_verified'] and B['all28_intrinsic_C0_columns_have_zero_lift']
 assert D['SyzQ']['actual_d2']['shape']==[48,128] and D['SyzQ']['actual_d2']['value']=='Res4(N^t)*Res4(S12_by32)'
 assert E['source_operator_proof']['native_source_definition']==manifest['actual_SyzQ_entry_definition']
 V=E['actual_V56_nodes'];assert len(V)==56
 for j,v in enumerate(V):assert N[v]==['actual_V_column_from_same_allorder_chain_restrict_ell_squared',j,D['adj']['actual_V56_source_columns'][j]]
 def one(tag):
  a=[i for i,n in enumerate(N) if n[0]==tag];assert len(a)==1,(tag,a);return a[0]
 def node(i,tag,*args):assert N[i]==[tag,*args],(i,tag)
 res=D['res'];betti=D['betti'];free26=D['adj']['conormal_standard_source_indices26'];sel89=betti['fixed_original_minor_columns'];free429=[i for i in range(518) if i not in sel89]
 raw=one('periodic_Koszul_C2_resolution');assert N[raw][1]['first_ranks']==[14,126,518,1302]
 assert N[raw][1]['periodic_d_p']=='ell I+(-1)^p L' and N[raw][1]['horizontal_d']=='(-1)^p sum_i(x_iI-A_i) exterior_contraction_i'
 a11=one('literal_a11_source_columns');node(a11,'literal_a11_source_columns',P['res'],res['source_a11_columns'])
 ds=one('selected_columns_of_d2');node(ds,'selected_columns_of_d2',raw,sel89)
 ef=one('source_standard_columns_F1');node(ef,'source_standard_columns_F1',126,free26)
 rr=one('nonphysical_target_projection_of_d1_Efree');node(rr,'nonphysical_target_projection_of_d1_Efree',raw,ef,res['physical_target3_basis_positions'])
 ec=one('Efree_minus_a11_R');node(ec,'Efree_minus_a11_R',ef,a11,rr)
 b1=one('F1_basis_a11_d2sel89_Ecorrect26');node(b1,'F1_basis_a11_d2sel89_Ecorrect26',a11,ds,ec)
 b1i=[i for i,n in enumerate(N) if n[0]=='formal_unit_matrix_inverse' and n[1]==b1];assert len(b1i)==1;b1i=b1i[0]
 z=one('middle89_coordinates_B1inv_d2_on_F2_complement');node(z,'middle89_coordinates_B1inv_d2_on_F2_complement',raw,b1i,free429)
 i2=one('F2_complement_inclusion_minus_selected89_Z');node(i2,'F2_complement_inclusion_minus_selected89_Z',518,free429,sel89,z)
 p2=one('F2_complement_projection');node(p2,'F2_complement_projection',518,free429)
 p1=one('last26_coordinates_B1inv');node(p1,'last26_coordinates_B1inv',b1i)
 i0=one('eta_physical3_in_C2_basis14');node(i0,'eta_physical3_in_C2_basis14',res['physical_target3_basis_positions'])
 p0=one('physical_m0_polynomial_projection');node(p0,'physical_m0_polynomial_projection',P['res'],res['physical_lift_m0_columns'])
 g1,g2,g3=E['actual_G1_G2_G3_nodes'];node(g1,'physical_d1_on_Ecorrect26',raw,ec,P['res']);node(g2,'compose',p1,['d2',raw],i2);node(g3,'compose',p2,['d3',raw])
 hp=one('exact_periodic_Koszul_k_linear_homotopy');assert N[hp][1]==raw and N[hp][2]['total_H']=='(-1)^p H_Kdelta + i_K h_periodic p_K'
 hg1=[i for i,n in enumerate(N) if n==['compose',p1,['H_degree0_to1',hp],i0]];hg2=[i for i,n in enumerate(N) if n==['compose',p2,['H_degree1_to2',hp],ec]];assert len(hg1)==len(hg2)==1
 m=E['actual_M26_by56_node'];assert N[m][:2]==['matrix_from_columns',[26,56]]
 assert N[m][2]==[next(i for i,n in enumerate(N) if n==['apply_vector',hg1[0],v]) for v in V]
 selected=D['conormal']['fixed_selected26_columns'];complement=[i for i in range(56) if i not in selected]
 mu=one('selected_source_columns');node(mu,'selected_source_columns',m,selected)
 mui=[i for i,n in enumerate(N) if n[0]=='formal_unit_matrix_inverse' and n[1]==mu];assert len(mui)==1
 assert N[mui[0]][2]['constant_determinant']==D['conormal']['fixed26_minor_determinant']
 n=E['actual_N56_by26_node'];node(n,'selected_source_inclusion_times_inverse',56,selected,mui[0])
 hf=E['actual_Hfree56_by30_node'];node(hf,'E_Cfree_minus_N_M_Cfree',56,complement,n,m)
 s=E['actual_S56_by459_node'];node(s,'kerV_source',[g2,hf],n,{'formula':'[N G2,Hfree]','shape':[56,459]})
 yy=E['all768_actual_cycle_lift_top_nodes'];bb=E['all768_actual_cycle_lift_bottom_nodes'];assert len(yy)==len(bb)==768
 for x,y,b in zip(X,yy,bb):
  assert N[y][0]=='apply_vector' and N[y][1]==hg2[0]
  node(N[y][2],'apply_vector',m,x);node(b,'select_vector_components',x,complement)
 final=E['complete_actual_T2_presentation_node'];node(final,'complete_T2_presentation_block',{'shape':[459,2070],'first1302_columns':['top',g3,'bottom_zero30'],'last768_columns':['top_vectors',yy,'bottom_vectors',bb]})
 # Known source identities and low constants justify ALL formal inverses.
 assert len(free26)==26 and len(sel89)==89 and len(free429)==429 and len(complement)==30
 assert betti['rank']==89 and D['conormal']['rank']==26
 assert D['source_peer']['allorder_gA505938_prefix0_source_and_primitiveV_gate_accepted']
 assert D['conormal_peer']['formal_unit_section_n_and_h0_source_gate_accepted']
 actions=[[sp(c) for c in M] for M in D['C2']['C2_coordinate_actions9_sparse_columns']];Eell=actions[8]
 assert all(not p for p in mp(Eell,Eell))
 assert all(mp(A0,B0)==mp(B0,A0) for A0 in actions for B0 in actions)
 count=HomotopyFixtures([actions[0],actions[7]],Eell).run();assert count==1344
 proof=E['source_operator_proof'];assert proof['old_ell0_m_resolution_or_X_coefficient_arrays_reused'] is False
 assert proof['new_scalar_or_T2_length_or_H0_query_performed'] is False
 ep=o/'independent-768source-100units-homotopy-and-completeT2-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all8107_original_fullR0_callback_restrictions_and768_source_columns_verified':True,'all640_B_P_and128_native_source_bindings_rebuilt':True,'ell0_coefficient_arrays_not_reused':True,'all11318_typed_node_interfaces_and459_by2070_blocks_verified':True,'allorder100_unit_cancellations_and_chain_projection_proved':True,'periodic_total_homotopy_sign_fixtures':count,'allorder_ordered_tensor_contraction_and_finite_translation_proved':True,'M_N_h0_free30_and768_actual_cycle_lifts_accepted':True,'global_I_annihilates_K_and_C_premise':'accepted original staged tensor normalization conductor theorem, not inferred from Ann(C2)','derived_reduction_sequence':'0 -> K/ell²K -> T2 -> C[ell²] -> 0 from original free comparison complex','I_squared_annihilates_T2':True,'finite_base_dimension_upper':696,'maximal_ideal_power10_annihilates_T2':True,'actual_finite_base_or_T2_coefficients_or_length_computed':False})
 report={'id':STEM,'all_passed':True,'actual_complete_T2_typed_presentation_459_by2070_accepted':True,'fullR0_ell_squared_source_and768_lifts_accepted':True,'M_N_h0_and100_unit_cancellations_accepted':True,'global_conductor_I_squared_and_finite696_m10_bounds_accepted':True,'new_source_coefficient_or_T2_length_or_H0_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
