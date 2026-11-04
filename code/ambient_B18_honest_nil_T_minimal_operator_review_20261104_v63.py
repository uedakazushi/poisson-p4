#!/usr/bin/env python3
"""Independent three-unit-Schur source/rank audit, no new source values."""
import hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
import s61_nonsquare_prePOWER_kappa11_schur_review_20261104_v63 as NC
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63';ID=PRIMARY+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(path):
  p=ROOT/path;pins[path]=sha(p);return json.loads(p.read_text())
 base='reports/'+PRIMARY+'/'
 r=load(base+'report.json');assert pins[base+'report.json']=='5e28bc2af334d0dd340a2300a520e2bd746482f4d25e6add2acc77c25c2b5dba'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 D=load(base+'actual-minimal-Tnil5-849-and-d3-linear-operator.json')
 S={k:load(v['path'])for k,v in D['source_artifacts'].items()}
 for k,v in D['source_artifacts'].items():assert pins[v['path']]==v['sha256']
 assert S['honest_peer']['all_passed']and S['values_peer']['all_passed']and S['values_report']['mu_R0_Tnil']==5
 H,M,X,W,K=[S[k]for k in ['honest','m0','X1','W0','Koszul']]
 N=D['operator_nodes'];sh=D['operator_shapes'];assert len(N)==len(sh)==60
 for i,n in enumerate(N):
  k=n[0];refs=[]
  if k=='actual_external_matrix':
   assert n[1]=='honest'and n[2]in H['matrices'];assert sh[i]==H['matrix_operator_shapes'][H['matrices'][n[2]]]
  elif k=='QQ_matrix_rows':assert sh[i]==[len(n[1]),len(n[1][0])]and all(len(r)==sh[i][1]for r in n[1])
  elif k=='identity':assert sh[i][0]==sh[i][1]
  elif k=='zero':pass
  elif k=='row_column_selection':
   refs=[n[1]];assert sh[i]==[len(n[2]),len(n[3])]and all(x in range(sh[n[1]][0])for x in n[2])and all(x in range(sh[n[1]][1])for x in n[3])
  elif k=='matrix_product':
   refs=n[1];assert all(sh[a][1]==sh[b][0]for a,b in zip(refs,refs[1:]))and sh[i]==[sh[refs[0]][0],sh[refs[-1]][1]]
  elif k=='matrix_sum':refs=n[1:];assert all(sh[x]==sh[i]for x in refs)
  elif k=='matrix_neg':refs=[n[1]];assert sh[n[1]]==sh[i]
  elif k=='origin_unit_matrix_inverse':
   refs=[n[1]];assert sh[i]==sh[n[1]]and sh[i][0]==sh[i][1]
   C=[Q.sparse(r)for r in n[2]['literal_origin_matrix']];assert len(C)==sh[i][0]and Q.det_fixed(C)==F(n[2]['determinant'])!=0
  elif k in ['horizontal_blocks','vertical_blocks']:
   refs=n[1];axis=0 if k=='horizontal_blocks'else 1;assert len({sh[j][axis]for j in refs})==1
   shape=sh[refs[0]].copy();shape[1-axis]=sum(sh[j][1-axis]for j in refs);assert sh[i]==shape
  else:raise AssertionError(k)
  assert all(j<i for j in refs)
 def selection(i,p,I,J):assert N[i]==['row_column_selection',p,I,J]
 m0=[[F(c)for c in row]for row in M['m0_24_by44']];E=M['selected_source_columns24'];C=[i for i in range(44)if i not in E]
 assert N[5][1]==[[str(int(i==j))for j in E]for i in range(44)]
 assert N[6]==['matrix_product',[0,5]]and N[7][1]==6 and N[8]==['matrix_product',[5,7]]
 assert [[F(c)for c in row]for row in N[7][2]['literal_origin_matrix']]==[[m0[i][j]for j in E]for i in range(24)]
 assert N[13]==['vertical_blocks',[10,12]]and N[10]==['horizontal_blocks',[1,8,9]]and N[12]==['horizontal_blocks',[2,11,4]]
 assert N[9]==N[11]==['zero']and sh[9]==[44,168]and sh[11]==[84,24]
 # Every stage's exact selection/Schur/projection wiring, not merely dimensions.
 original=list(range(972));targets=list(range(128));cancelled=[]
 for stage in D['three_exact_unit_Schur_steps']:
  p=stage['input'];I=stage['unit_rows'];J=stage['unit_columns'];Ir=[i for i in range(sh[p][0])if i not in I];Jr=[j for j in range(sh[p][1])if j not in J]
  assert Ir==stage['retained_rows']and Jr==stage['retained_columns']
  u,b,c,d=[stage[k]for k in ['unit_block','B','C','D']]
  selection(u,p,I,J);selection(b,p,I,Jr);selection(c,p,Ir,J);selection(d,p,Ir,Jr)
  ui=stage['unit_inverse'];assert N[ui]==['origin_unit_matrix_inverse',u,{'literal_origin_matrix':stage['unit_origin_matrix'],'determinant':stage['unit_origin_determinant']}]
  s=stage['Schur'];assert N[s][0]=='matrix_sum'and N[s][1]==d
  minus=N[s][2];assert N[minus][0]=='matrix_neg'and N[N[minus][1]]==['matrix_product',[c,ui,b]]
  proj=stage['target_projection'];assert N[proj][0]=='matrix_sum';pr,minus=N[proj][1:];selection(pr,N[pr][1],Ir,list(range(sh[p][0])));assert N[N[pr][1]]==['identity']
  assert N[minus][0]=='matrix_neg';mul=N[N[minus][1]];assert mul[0]=='matrix_product'and mul[1][:2]==[c,ui]
  pi=mul[1][2];selection(pi,N[pi][1],I,list(range(sh[p][0])));assert N[N[pi][1]]==['identity']
  cancelled.append([original[j]for j in J]);original=[original[j]for j in Jr];targets=[targets[i]for i in Ir]
 assert [len(x)for x in cancelled]==[24,18,81]
 assert [x['cancelled_original_relation_columns']for x in D['all972_original_relation_routes']]==cancelled
 assert original==D['retained_original_relation_columns849']and len(original)==849
 assert sorted(original+sum(cancelled,[]))==list(range(972))and sum(i<780 for i in original)==681 and [i for i in original if i>=804]==list(range(804,972))
 S1,S2,S3=D['three_exact_unit_Schur_steps']
 assert S1['input']==13 and S1['unit_rows']==E and S1['unit_columns']==list(range(780,804))and S1['Schur']==21
 # n is supported only on E: lower section block is zero to ALL orders.
 assert S1['retained_rows']==C+list(range(44,128))
 mi=[[F(c)for c in row]for row in M['selected_minor_inverse']]
 assert [[F(c)for c in row]for row in S1['unit_origin_matrix']]==mi and F(S1['unit_origin_determinant'])==1/F(M['selected_minor_determinant'])
 xp=X['X0_pivot_columns'];const=[{int(i):dict((int(k),F(v))for k,v in sp).get(0,F(0))for i,sp in col}for col in X['all780_source_columns44_sparse_dual9']]
 Jrows=D['fixed_minor_data']['compatible_X0_top18_rows_in_C20'];J0=[[const[j].get(C[i],F(0))for j in xp]for i in Jrows]
 assert S2['input']==21 and S2['unit_rows']==Jrows and S2['unit_columns']==xp and S2['Schur']==36
 assert J0==[[F(c)for c in row]for row in S2['unit_origin_matrix']]
 Ji=[Q.sparse(r)for r in D['fixed_minor_data']['inverse']];Jsp=[Q.sparse(row)for row in J0]
 assert Q.mprod(Jsp,Ji)==Q.mprod(Ji,Jsp)==Q.eye(18)
 for j,c in enumerate(X['all780_X0_factorization_coefficients']):
  co=Q.sparse(c);assert Q.mvec(Jsp,co)=={i:const[j].get(C[k],F(0))for i,k in enumerate(Jrows)if const[j].get(C[k],F(0))}
 assert S2['retained_columns'][:762]==W['residual_original_X_columns']
 wcols=[Q.pairs(c)for c in W['all_W84_columns_sparse']];wr=W['W_pivot_rows'];wp=W['W_pivot_residual_columns'];wi=[Q.sparse(r)for r in W['W_fixed_minor_inverse']]
 wj=[[wcols[j].get(i,F(0))for j in wp]for i in wr]
 assert S3['input']==36 and S3['unit_rows']==[2+i for i in wr]and S3['unit_columns']==wp and S3['Schur']==51
 assert wj==[[F(c)for c in row]for row in S3['unit_origin_matrix']]
 assert N[59]==['matrix_product',[S3['target_projection'],S2['target_projection'],S1['target_projection']]]and sh[51]==[5,849]and sh[59]==[5,128]
 NC.nc_block_proof() # Noncommuting two-sided Schur proof; applies at each origin-unit block.
 # Constant d3 projection onto the surviving THREE wedge coordinates.
 small=D['d3_only_linear'];wc=[i for i in range(84)if i not in wr];assert wc==small['remaining_wedge_indices3']==[6,24,26]
 projection=[]
 for t in wc:
  row={t:F(1)}
  for j,r in enumerate(wr):
   c=sum(wcols[wp[k]].get(t,F(0))*wi[k].get(j,F(0))for k in range(81))
   if c:Q.put(row,r,-c)
  projection.append(row)
 assert projection==[Q.sparse(r)for r in small['constant_projection3_by84']]
 assert all(not Q.mvec(projection,w)for w in wcols)
 known=[]
 for col in K['all8_differentials_sparse_columns'][2]:
  out={}
  for row,mon,sgn in col:
   if sum(mon)==1:
    j=mon.index(1)
    for t,p in enumerate(projection):
     if row in p:Q.put(out,8*(t+2)+j,p[row]*sgn)
  known.append(out)
 assert known==[Q.pairs(c)for c in small['all168_known_minimal_linear_columns40_sparse']]
 kp=small['pivot_columns'];kr=small['pivot_rows'];ki=[Q.sparse(row)for row in small['minor_inverse']]
 minor=[{j:known[c].get(r,F(0))for j,c in enumerate(kp)if known[c].get(r,F(0))}for r in kr]
 assert len(kp)==len(kr)==15 and Q.mprod(minor,ki)==Q.mprod(ki,minor)==Q.eye(15)and Q.det_fixed(minor)==F(small['minor_determinant'])
 for col in known:
  co=Q.mvec(ki,{i:col[r]for i,r in enumerate(kr)if r in col})
  assert Q.plus(*(Q.scaled(known[kp[i]],a)for i,a in co.items()))==col
 assert small['is_full_linear_relation_space']is False and small['known_first_Loewy_dimension_upper']==25
 # Audit the missing-adapter2 obstruction, keeping nonzero old Bold0 residuals.
 b0=[]
 for col in S['oldBfirst']['source_columns640_rows56_sparse_dual9']:
  d={}
  for i,v in col:
   if i<48:
    a=Q.pairs(v).get(0,F(0))
    if a:d[i]=a
  b0.append(d)
 b0 += [{}for _ in range(140)];survive=[]
 for j in W['residual_original_X_columns']:
  c=Q.sparse(X['all780_X0_factorization_coefficients'][j]);v=Q.plus(b0[j],*(Q.scaled(b0[xp[k]],-a)for k,a in c.items()))
  if v:survive.append([j,[[i,str(a)]for i,a in sorted(v.items())]])
 prof=D['full_first_layer_missing_data_profile'];bs=prof['X2_data'];assert survive==bs['constant_residual_Bold0_nonzero_list']and len(survive)==bs['constant_residual_Bold0_nonzero_count']==9
 support=sorted({i for j,v in survive for i,a in v});assert support==bs['needed_ResAT2_old_columns_only']==[32,36,40,44]and bs['ResAT2_Bold0_drop_claimed']is False
 bm=[Q.sparse(r)for r in bs['surviving_Bold0_span_minor']];bi=[Q.sparse(r)for r in bs['surviving_Bold0_span_inverse']]
 assert Q.mprod(bm,bi)==Q.mprod(bi,bm)==Q.eye(4)and Q.det_fixed(bm)==F(bs['surviving_Bold0_span_determinant'])
 vectors=[dict((i,F(a))for i,a in v)for j,v in survive]
 for col in range(4):assert {support[i]:bm[i][col]for i in range(4)if col in bm[i]}in vectors
 assert prof['no_X2_m1_values_computed']and prof['full_first_Loewy_layer_or_nil_generators_or_cube_decided']is False
 labels=D['minimal_target_labels'];assert [x.get('index')for x in labels]==[42,43,6,24,26]
 save('independent-three-Schur-and-d3-proof.json',{'all60_matrix_nodes_actual_honest_namespaces_and_source_paths_checked':True,'all3_origin_unit_blocks_and_24_18_81_fixed_minors_verified':True,'section_supported_on_E_lower_block_zero_to_all_orders':True,'all972_relations_routed_once_681X_plus_all168d3_retained':True,'complete_minimal_presentation_shape':[5,849],'all_target_projection_and_noncommuting_Schur_identities':True,'minimal_target_indices':[42,43,6,24,26],'all168_literal_d3_firsts_and_fixed_rank15_inverse_determinant_factorizations':True,'d3_only_not_full_linear_relation_space':True,'all9_nonzero_Bold0_residuals_support4_rank4_verified':True,'m1_X1_and_adapter2_Bold0_not_discarded':True,'no_new_source_jet_X2_m1_length_nilgens_or_cube_claim':True})
 for helper in [Path(Q.__file__),Path(NC.__file__)]:pins[str(helper.relative_to(ROOT))]=sha(helper)
 rep={'id':ID,'all_passed':True,'frozen':True,'complete_minimal_target_shape':[5,849],'all_order_unit_pair_counts':[24,18,81],'typed_matrix_nodes':60,'all972_source_relations_and_actual_honest_source_maps_accepted':True,'known_d3_only_linear_rank':15,'full_linear_layer_NOT_computed':True,'nonzero_Bold0_residual_count':9,'no_source_values_jets_X2_m1_scalar_minimal849_or_nil_cube_computation':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{ID}.md'),'evidence_sha256':{'independent-three-Schur-and-d3-proof.json':sha(OUT/'independent-three-Schur-and-d3-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in OUT.iterdir())<32*1024*1024;save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':main()
