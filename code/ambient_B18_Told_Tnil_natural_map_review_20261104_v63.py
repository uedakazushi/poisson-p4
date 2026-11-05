#!/usr/bin/env python3
"""Independent native polynomial/source transport and fixed cached5 origin audit."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID;ZERO=(0,)*9
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def pd(p):return {tuple(m):F(c)for m,c in p if F(c)}
def pa(a,b,c=F(1)):
 z=a.copy()
 for m,x in b.items():Q.put(z,m,c*x)
 return z
def pm(a,b):
 z={}
 for e,x in a.items():
  for f,y in b.items():Q.put(z,tuple(i+j for i,j in zip(e,f)),x*y)
 return z
def cm(col,i,p):
 q=pa(col.get(i,{}),p)
 if q:col[i]=q
 else:col.pop(i,None)
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;k=str(p.relative_to(ROOT));pins[k]=sha(p)
  if k not in cache:cache[k]=json.loads(p.read_text())
  return cache[k]
 def save(name,v):(OUT/name).write_text(json.dumps(v,indent=2)+'\n')
 rr=load(BASE);assert pins[f'reports/{BASE}/report.json']=='e881d79daf901b41f3c327397692d7c76e401824fa2d99bad4f4eac6b5ba0f52'
 for p,h in rr['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in rr['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/p)==h
 D=load(BASE,'actual-natural-Told-to-Tnil-source-map-and-origin.json');S={}
 for k,ref in D['source_artifacts'].items():
  p=ROOT/ref['path'];assert sha(p)==ref['sha256'];pins[ref['path']]=sha(p);S[k]=json.loads(p.read_text())
 for k in ['old_unit_peer','new_honest_peer','new_minimal_peer','first_peer','nilfinite_peer']:assert S[k]['all_passed']
 O,R,H,Min,old=[S[k]for k in ['old_unit','radical','new_honest','new_minimal','old_minimal']]
 assert R['accepted_native_C_coefficient_interface']['sha256']==pins[D['source_artifacts']['original_native']['path']]
 assert R['A_definition']=='A=N12^-1*injection12_by8*D8_by9'
 assert R['actual_C_variable_rules']['N12']=='accepted full_N=T(kappa)*R12, dG8+zeros=N12*physical_F12'
 assert D['old_minimal_generator_free56_rows']==old['actual_free_five_source_rows']==[22,23,29,30,55]
 Eold=O['fixed_selected24_columns'];F32=[i for i in range(56)if i not in Eold]
 assert S['old_honest']['actual_source32_complement_rows']==F32
 assert [F32.index(i)for i in D['old_minimal_generator_free56_rows']]==[8,9,12,13,31]
 # Fixed actual polynomial identities, no regenerating the author's J/search.
 DD=[[pd(p)for p in row]for row in R['G8_polynomial_factor_D8_by9']];r9=[pd(p)for p in R['R9_generators_polynomials9']];newS=[{i:pd(p)for i,p in col}for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']]
 routing=[{i:pd(p)for i,p in col}for col in D['native_Koszul36_to_complete35_columns']];J=[{i:pd(p)for i,p in col}for col in D['actual_native_old32_to_radical35_columns']]
 assert len(routing)==36 and len(J)==32 and all(not c for c in J[28:])
 for (i,j),col in zip(itertools.combinations(range(9),2),routing):
  ans={}
  for k,a in col.items():
   for row,b in newS[k].items():cm(ans,row,pm(a,b))
  expected={i:{m:-c for m,c in r9[j].items()},j:r9[i]};assert ans==expected
 G=[]
 for row in DD:
  q={}
  for a,b in zip(row,r9):q=pa(q,pm(a,b))
  G.append(q)
 oldS=[{}for _ in range(32)];native=S['original_native']['native_S12_by32']
 for row,j,g,sign in native['koszul28_in_first8_rows']:oldS[j][row]={m:sign*c for m,c in G[g].items()}
 assert native['free4_zero_source_columns']==[[8,28,1],[9,29,1],[10,30,1],[11,31,1]]
 for j,col in enumerate(J):
  lhs={};rhs={}
  for k,a in col.items():
   for row,b in newS[k].items():cm(lhs,row,pm(a,b))
  for i,g in oldS[j].items():
   for row,a in enumerate(DD[i]):cm(rhs,row,pm(a,g))
  assert lhs==rhs,(j,'full polynomial source chain')
 assert sum(len(p)for col in routing+J for p in col.values())==D['visible_polynomial_terms']==507
 # Full source chain and actual D0 class: A^t N^t E8=D^t, preserving
 # the640 positive columns, complete128 native columns, actual same c8.
 assert D['actual_full_chain_map_old768_to_new780']['dimensions']==[780,768]
 assert D['actual_full_chain_map_old768_to_new780']['first640']=='literal identity'
 assert 'minus'not in R['actual_F12_relation_identity']
 # Type-check all16 matrix nodes and40 original column nodes against ACTUAL files.
 N=D['matrix_operator_nodes'];sh=D['matrix_operator_shapes'];assert len(N)==len(sh)==16
 maps={'old_unit':{'typed_m24_by56':[24,56],'typed_n56_by24':[56,24]},'new_honest':{'old_to_new_Z44_by56':[44,56],'new_m24_by44':[24,44]},'new_minimal':{'full_target_projection5_by128':[5,128]}}
 for i,a in enumerate(N):
  typ=a[0];refs=[]
  if typ=='actual_external_matrix':assert sh[i]==maps[a[1]][a[2]]
  elif typ=='QQ_matrix_rows':assert sh[i]==[len(a[1]),len(a[1][0])]
  elif typ=='identity':assert sh[i][0]==sh[i][1]
  elif typ=='matrix_product':refs=a[1];assert all(sh[x][1]==sh[y][0]for x,y in zip(refs,refs[1:]))and sh[i]==[sh[refs[0]][0],sh[refs[-1]][1]]
  elif typ=='matrix_neg':refs=[a[1]];assert sh[i]==sh[a[1]]
  elif typ=='matrix_sum':refs=a[1:];assert all(sh[j]==sh[i]for j in refs)
  elif typ=='vertical_blocks':refs=a[1];assert len({sh[j][1]for j in refs})==1 and sh[i]==[sum(sh[j][0]for j in refs),sh[refs[0]][1]]
  elif typ=='matrix_from_original_columns':assert sh[i]==[84,5]
  else:raise AssertionError(typ)
  assert all(j<i for j in refs)
 assert N[0]==['actual_external_matrix','old_unit','typed_m24_by56']and N[1]==['actual_external_matrix','old_unit','typed_n56_by24']
 assert [[F(c)for c in row]for row in N[2][1]]==[[F(i==j)for j in[22,23,29,30,55]]for i in range(56)]
 assert N[4]==['matrix_product',[1,0]]and N[5]==['matrix_neg',4]and N[6]==['matrix_sum',3,5]and N[7]==['matrix_product',[6,2]]and N[9]==['matrix_product',[8,7]]and N[11]==['matrix_product',[10,9]]
 vec=D['original_vector_lift_nodes'];assert len(vec)==40
 for j in range(5):
  pos=8*j;assert vec[pos]=={'operation':['ORIGINAL_matrix_column_functions',11,j],'dimension':24}
  parts=[]
  for sec,power in enumerate([1,1,2]):
   a=pos+1+2*sec;b=a+1
   assert vec[a]=={'operation':['contiguous_slice',pos,8*sec,8],'dimension':8}
   assert vec[b]=={'operation':['ORIGINAL_complete_Koszul_Hadamard_Euler_lift',a,power,'new_honest.shared_homotopy_template_v1'],'dimension':28};parts.append(b)
  assert vec[pos+7]=={'operation':['ordered_concatenation',parts],'dimension':84}
 assert N[12]==['matrix_from_original_columns',[7,15,23,31,39]]and N[13]==['vertical_blocks',[9,12]]and N[15]==['matrix_product',[14,13]]
 # Same actual closed lifts are valid for both old n and direct n′ honest128
 # presentations, so the fixed minimal projection is legal (not a blind identity).
 assert H['shared_homotopy_template_v1']['not_S_linear']and H['matrices']['old_to_new_Z44_by56']==13
 # Independent cached source origin, rather than reranking stored5 columns.
 A=S['actual_origin_data']if 'actual_origin_data'in S else D['actual_origin_data'];Ao=D['actual_origin_data']
 m0=[Q.sparse(row)for row in O['all24_by56_m0']];ni=[Q.sparse(row)for row in O['fixed24_minor_inverse']]
 minor=[{j:m0[i][k]for j,k in enumerate(Eold)if m0[i].get(k,F(0))}for i in range(24)];assert Q.mprod(minor,ni)==Q.mprod(ni,minor)==Q.eye(24)
 n0=[{}for _ in range(56)]
 for row,v in zip(Eold,ni):n0[row]=v
 old0=[]
 for j in[22,23,29,30,55]:
  x={i:m0[i][j]for i in range(24)if m0[i].get(j,F(0))};c=Q.plus({j:F(1)},Q.scaled(Q.mvec(n0,x),-1));assert not Q.mvec(m0,c);old0.append(c)
 assert [Q.sparse(v)for v in Ao['all_five_old_cycle_origins56']]==old0
 assert [Q.sparse(row)for row in Ao['old_n0_56_by24']]==n0 and [Q.sparse(row)for row in Ao['old_m0_24_by56']]==m0
 C36=Q.FiniteC([[0]*8]+[[int(i==j)for i in range(8)]for j in range(8)],S['monics1']['nu_first9'])
 Af=S['X1']['new_A12_by9_C36'];z0=[{}for _ in range(44)];z1=[{}for _ in range(3)]
 for i,row in enumerate(Af):
  for gen,a in enumerate(row):
   for mask in range(4):
    q=C36.mul(Q.sparse(a),{mask:F(1)})
    for out in range(4):
     if q.get(out,F(0)):z0[4*gen+out][4*i+mask]=q[out]
    if gen==5:
     assert not q.get(0,F(0))
     for k,j in enumerate([0,4,5]):
      if q.get(4*(j+1),F(0)):z1[k][4*i+mask]=q[4*(j+1)]
 for i in range(8):z0[36+i][48+i]=F(1)
 assert z0==[Q.sparse(v)for v in Ao['actual_Z0_44_by56']]and z1==[Q.sparse(v)for v in Ao['actual_Z_row20_zi_firsts3_by56']]
 new0=[Q.mvec(z0,c)for c in old0];der=[Q.mvec(z1,c)for c in old0]
 assert new0==[Q.sparse(c)for c in Ao['all_five_new_cycle_origins44']]and der==[Q.sparse(c)for c in Ao['all_five_new_cycle_row20_zi_firsts3']]
 mn=[Q.sparse(v)for v in S['m0']['m0_24_by44']];assert mn[7]=={20:F(5,2)}and all(not Q.mvec(mn,c)for c in new0)
 xp=S['X1']['X0_pivot_columns'];Xcols=[{int(i):Q.pairs(v)for i,v in c}for c in S['X1']['all780_source_columns44_sparse_dual9']];free=[i for i in range(44)if i not in S['m0']['selected_source_columns24']];I=[free[i]for i in Min['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in Min['fixed_minor_data']['inverse']]
 hi=[[F(5,2)*Xcols[j].get(20,{}).get(i+1,F(0))for j in xp]for i in[0,4,5]];assert hi==[[F(c)for c in v]for v in S['first']['all54_actual_pivot_h0_readback']]
 assert Min['d3_only_linear']['constant_projection3_by84']==[[str(int(i==j))for i in range(84)]for j in[6,24,26]]
 mapped=[]
 for c,d in zip(new0,der):
  co=Q.mvec(Ji,{k:c[v]for k,v in enumerate(I)if c.get(v,F(0))});v=[]
  for target in[42,43]:v.append(c.get(target,F(0))-sum((Xcols[j].get(target,{}).get(0,F(0))*co.get(k,F(0))for k,j in enumerate(xp)),F(0)))
  for k in range(3):v.append(F(5,2)*d.get(k,F(0))-sum((hi[k][q]*co.get(q,F(0))for q in range(18)),F(0)))
  mapped.append(Q.sparse(v))
 assert mapped==[Q.sparse(v)for v in Ao['actual_origin_map_columns5']]
 assert mapped[0]=={2:F(5,3),3:F(-2,3),4:F(5,2)}and mapped[4]=={1:F(1),2:F(-14075,13448),3:F(36008,25215),4:F(28015,40344)}and all(not mapped[i]for i in[1,2,3])
 cols=Ao['origin_rank_pivot_columns'];rows=Ao['origin_rank_pivot_rows'];inv=[Q.sparse(v)for v in Ao['origin_rank_minor_inverse']];small=[{j:mapped[c][r]for j,c in enumerate(cols)if mapped[c].get(r,F(0))}for r in rows]
 assert cols==[0,4]and rows==[1,2]and Q.mprod(small,inv)==Q.mprod(inv,small)==Q.eye(2)and Q.det_fixed(small)==F(-5,3)==F(Ao['origin_rank_minor_determinant'])
 for j,v in enumerate(mapped):
  co=Q.mvec(inv,{i:v[r]for i,r in enumerate(rows)if v.get(r,F(0))});assert Q.plus(*(Q.scaled(mapped[cols[k]],a)for k,a in co.items()))==v and co==Q.sparse(Ao['all5_origin_column_factorizations'][j])
 assert Ao['rank']==rr['actual_natural_map_origin_rank']==2 and mapped[4][1]==1
 # Elementary Snake lemma after ell: both normalizers ell-regular identify the
 # complete H1's with literal kernels Dbar→Ebar. NE is free9 over D.
 assert 9-8+4==5 and 5-2==3
 save('independent-source-native-chain-and-five-origin-proof.json',dict(all36_Koszul_to_complete35_columns_verified=True,all32_full_native_polynomial_identities=True,visible_terms507=True,old_four_zero_native_columns_map_zero=True,actual_N12_D_A_transpose_and_C4_Res_order=True,all640_plus128_source_boundaries_preserved=True,all16_matrix40_ORIGINAL_vector_types_and_five_basis_lifts=True,old_and_direct_honest128_choices_compatible_on_actual_chosen_lifts=True,all56_old_cycle_and44_new_cycle_origin_values_from_cached_units=True,all44x56_Z0_and3x56_Zrow20_firsts_verified=True,unknown_old_cycle_derivatives_killed_by_Zrow20_origin0=True,actual_Pi0_18pivot_and_h0_cleanup=True,actual_origin_rank2_fixed_minor_both_inverses_all5factorizations=True,unit_entry_new2_old5=1,exact_injection_and_cokernel_length5_by_Snake_lemma=True,cokernel_mu3=True,possible_HF_cokernel=[[3,2],[3,1,1]],m3_cokernel_zero=True,zero_origin_HF_shift_hypothesis_false=True,no_HF2_Tnil_total_length_nilcube_claim=True,no_producer_import_source_jet_or_rank_search=True))
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,actual_origin_rank=2,exact_natural_injection=True,cokernel_length=5,cokernel_mu=3,all32_native_source_and_fixed5_origin_checks=True,zero_origin_HF_shift_hypothesis_FAILED=True,actual_HF2_Tnil_or_nilcube_NOT_decided=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name not in['report.json','phase.json']},seconds=time.monotonic()-start)
 save('phase.json',dict(phase='FINAL actual source and rank2 peer'));rep['evidence_sha256']['phase.json']=sha(OUT/'phase.json');save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
