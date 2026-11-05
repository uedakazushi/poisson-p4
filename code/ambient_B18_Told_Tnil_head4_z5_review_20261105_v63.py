#!/usr/bin/env python3
"""Independent projected natural-map source/value and filtration proof.
Only existing physical source2 is composed. No producer, solve or new jet.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
from ambient_B18_honest_nil_T_gauge120_review_20261105_v63 import Base
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
M2=list(itertools.combinations_with_replacement(range(8),2));M3=list(itertools.combinations_with_replacement(range(8),3));P8=list(itertools.combinations(range(8),2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 sid='AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_SOURCE_PREFLIGHT_20261105_v63';vid='AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63'
 sr=load(sid);vr=load(vid);assert sr['all_passed']and vr['all_passed']
 assert pins[f'reports/{sid}/report.json']=='b291d472cf9a47d706ae4972a79fd8a96d9a8d420a55978fcc5a92a2c26fac23'
 assert pins[f'reports/{vid}/report.json']=='b25e9cf95c649d6a789d6f05d027caf8cced24a8b6c40e5852850dffe172fd4c'
 for stem,rep in[(sid,sr),(vid,vr)]:
  for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
  for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/n)==h
 H=load(sid,'head4-z5-projected-source-and-cache-gate.json');V=load(vid,'actual-head4-z5-map-and-F2-kernel-coefficient.json')
 names={
 'nat':('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json'),
 'oldunit':('AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63','actual-fullV56-m0-and-typed-comparison.json'),
 'X':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json'),
 'm0':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json'),
 'min':('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json'),
 'first':('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json'),
 'W':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase3-actual-W84-rank-and-mu.json'),
 'C1':('AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63','actual81-pivot-first-carries-and-six-H1-source-adjoints.json'),
 'q':('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63','literal-known-cubic-input-and-exact24-row-reduction.json'),
 'cub':('AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63','complete-cubic917-typed-source-and-three-NF-adjoints.json'),
 'rg':('AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63','right106-full-first-span-and-cache-footprint.json'),
 'rv':('AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63','phase1-all35-nine-convolution-projections.json'),
 'physics':('AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63','phase1-common105-and-physical-source-base25.json'),
 'source':('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json'),
 'source2':('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63','full9-firstjet-certificate.json'),
 'fil':('AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63','actual-F2-source-products-and-m5-zero-proof.json'),
 'zero3':('AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63','actual-three-first-columns-and-fixed-quotient-proof.json')}
 D={k:load(*p)for k,p in names.items()}
 for p in['AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW','AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_GAUGE120_20261105_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_RIGHT106_20261105_v63_REVIEW','AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW']:assert load(p)['all_passed']
 exact=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63');assert exact['actual_T_length']==11 and exact['actual_Tnil_length']==16
 save('phase.json',dict(phase='all actual source/filtration/cache peers pinned; fixed source proof before value'))
 N=D['nat']['actual_origin_data'];cold=Q.sparse(N['all_five_old_cycle_origins56'][4]);cnew=Q.sparse(N['all_five_new_cycle_origins44'][4]);origin=Q.sparse(N['actual_origin_map_columns5'][4]);n0=[Q.sparse(v)for v in N['old_n0_56_by24']];Z0=[Q.sparse(v)for v in N['actual_Z0_44_by56']]
 assert set(cold)=={1,2,17,18,21,25,28,55};assert F(N['old_m0_24_by56'][7][55])==F(-125,24)
 O=D['oldunit'];ops=O['typed_regular_function_operators'];sel=O['fixed_selected24_columns'];un=O['typed_unit24_inverse_node'];assert ops[un][0]=='formal_unit_matrix_inverse'
 assert ops[un][1]==[[row[j]for j in sel]for row in O['typed_m24_by56']]
 for i,row in enumerate(O['typed_n56_by24']):
  for j,k in enumerate(row):assert ops[k]==(['unit_matrix_inverse_entry',un,sel.index(i),j]if i in sel else['zero'])
 # n=E(mE)^-1 gives n1=-n0*m1*n0, hence FULL c1=-n0*m1*c0.
 # There is no use of ZERO3's m0e=0 shortcut for this nonzero head.
 X=D['X'];raw=[{int(i):Q.pairs(v)for i,v in col}for col in X['all780_source_columns44_sparse_dual9']];xp=X['X0_pivot_columns'];newm0=[Q.sparse(v)for v in D['m0']['m0_24_by44']];free=[i for i in range(44)if i not in D['m0']['selected_source_columns24']];I=[free[i]for i in D['min']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in D['min']['fixed_minor_data']['inverse']]
 a0=Q.mvec(Ji,{j:cnew.get(i,F(0))for j,i in enumerate(I)});assert [a0.get(i,F(0))for i in range(18)]==list(map(F,H['all18_first_origin_cleanup_a0']))
 assert all(not raw[p].get(42,{}).get(0,F(0))for p in xp)
 oldadj={j:-sum((a*n0[k].get(j,F(0))for k,a in Z0[42].items()),F(0))for j in range(24)};oldadj={j:c for j,c in oldadj.items()if c}
 assert oldadj=={6:F(1),8:F(-5,9),12:F(-125,36),13:F(-5,9),14:F(-4)}
 assert [oldadj.get(i,F(0))for i in range(24)]==list(map(F,H['actual_old_m1_r_adjoint24']))
 W=D['W'];wi=[Q.sparse(v)for v in W['W_fixed_minor_inverse']];wr=W['W_pivot_rows'];C1=D['C1']['all81_actual_visible_C1_top2_by8_and_wedge3_r'];hadj={h:sum((F(C1[k][0][7])*wi[k].get(j,F(0))for k in range(81)),F(0))for j,h in enumerate(wr)};hadj={h:c for h,c in hadj.items()if c}
 assert hadj=={27:F(-8,3)}and P8[27]==(6,7)and newm0[7]=={20:F(5,2)}
 pi=[Q.pairs(v)for v in D['q']['quadratic_projection3_by180_sparse']];im2={m:i for i,m in enumerate(M2)};factors=[]
 for row in pi:
  adj={8*t+v:row[36*t+im2[tuple(sorted((v,5)))]]for t in range(5)for v in range(8)if row.get(36*t+im2[tuple(sorted((v,5)))],F(0))};assert set(adj)<={7};factors.append(adj.get(7,F(0)))
 assert factors==list(map(F,H['quadratic_projection_factors_times_Jr']))
 F1=D['first'];R35=F1['rank_minor_columns'];L=[Q.pairs(F1['all849_actual_minimal_first_columns40_sparse'][j])for j in R35];rr=F1['rank_minor_rows'];iv=[Q.sparse(v)for v in F1['rank_minor_inverse']];rhs={8*i+5:c for i,c in origin.items()};lift=Q.mvec(iv,{j:rhs[r]for j,r in enumerate(rr)if rhs.get(r,F(0))});rebuilt={}
 for j,c in lift.items():rebuilt=Q.plus(rebuilt,Q.scaled(L[j],c))
 assert rebuilt==rhs and [lift.get(i,F(0))for i in range(35)]==list(map(F,H['exact_R35_constant_lift35']))
 delta=[Q.pairs(v)for v in D['cub']['cubic_NF_dual3_by600_sparse']];im3={m:i for i,m in enumerate(M3)}
 for f in range(3):assert {36*t+j:delta[f][120*t+im3[tuple(sorted(mon+(7,)))]]for t in range(5)for j,mon in enumerate(M2)if delta[f].get(120*t+im3[tuple(sorted(mon+(7,)))],F(0))}==pi[f]
 G=D['rg'];res=[list(map(F,v))for v in G['all35_by8_k42_residual_coefficients']];eff={(f,v):{(h,w):F(c)for h,w,c in e}for f,v,e in G['all_actual_effective_C1_H2_callbacks']};us=[];vs=[]
 for f in range(3):
  us.append(sum((c*res[j][zi]*delta[f].get(120*(n+2)+119,F(0))for j,c in lift.items()for n,zi in enumerate([0,4,5])),F(0)))
  vs.append(-sum((c*res[j][i]*eff[f,7].get((56+P8.index((i,7)),7),F(0))for j,c in lift.items()for i in range(7)),F(0)))
 assert us==vs==[F(0)]*3
 save('phase1-independent-projected-source-and-formal-cancellation.json',dict(full_old_n1_retained=True,old_cycle_first_formula='-n0*m1*c0',old_m0e55_nonzero=F(-125,24),old_m1_adjoint=oldadj,H0_carry_adjoint=hadj,pi3_z5_factors=factors,R35_lift=lift,all40_actual_first_source_equalities=True,U_sensitivities=us,V_sensitivities=vs,all_last_r_division_adjoints_zero=True))
 # Independently recompose the48 actual physical coefficient inputs from
 # the immutable scalar DAG and accepted105 actual source2 leaf values.
 P=D['physics'];B=Base(P['base25_monomials']);MI=B.idx;rows=[18,19,28,29,60,67,72,63,64,68,73,74];SRC=D['source'];nodes=SRC['DAG_scalar_nodes'];leaves={i:Q.sparse(v)for i,v in P['all105_actual_source2_values']};family=SRC['families'];roots=[i for f in['CA','CB']for row in rows for i in family[f]['complete_normalized_physical_w_coefficients'][0][row]]+[family[f]['actual_tracefree_nu']for f in['CA','CB']]+[7790];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i)
  if i in leaves:continue
  n=nodes[i];assert n[0]in['c','a','m','n','i']
  if n[0]!='c':todo.extend(n[1:])
 uses=collections.Counter(roots)
 for i in need:
  if i not in leaves and nodes[i][0]!='c':uses.update(nodes[i][1:])
 val={};peak=0
 for i in sorted(need):
  n=nodes[i];o=n[0]
  if i in leaves:z=leaves[i]
  elif o=='c':z={0:F(n[1],n[2])}if n[1]else{}
  elif o=='a':z=Q.plus(val[n[1]],val[n[2]])
  elif o=='m':z=B.mul(val[n[1]],val[n[2]])
  elif o=='n':z=Q.scaled(val[n[1]],-1)
  else:z=B.inv(val[n[1]])
  val[i]=z;peak=max(peak,sum(map(len,val.values())));assert peak+35700<150000
  if i not in leaves and o!='c':
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del val[j]
 pp={f:[[val[i]for i in family[f]['complete_normalized_physical_w_coefficients'][0][row]]for row in rows]for f in['CA','CB']};nus=[val[family[f]['actual_tracefree_nu']]for f in['CA','CB']]
 for f in['CA','CB']:
  for j,row in enumerate(rows):assert pp[f][j]==[Q.sparse(v)for v in P['physical150_tracepairs_base25'][f][row]]
 assert nus==[Q.sparse(v)for v in P['nu_A_B_base25']]
 jet=D['source2'];gj={}
 for mon,c in zip(jet['source2_monomials'],jet['both57_source2']['CA']['pure_function2jets']['7790']):
  k=MI.get(tuple(sorted(7 if i==8 else i for i in mon)))
  if k is not None:Q.put(gj,k,F(c))
 assert gj==val[7790]
 C=Q.FiniteC([tuple(m.count(i)for i in range(8))for m in B.mon],[[v.get(i,F(0))for i in range(25)]for v in nus]);cols=[]
 assert [x['physical_row']for x in O['V56_column_order'][:48:4]]==rows
 for j in range(12):
  parts=[Q.plus(pp['CA'][j][0],Q.scaled(pp['CB'][j][0],-1)),pp['CA'][j][1],Q.scaled(pp['CB'][j][1],-1),{}];f={4*i+k:c for k,v in enumerate(parts)for i,c in v.items()}
  cols.extend(C.mul(f,{mask:F(1)})for mask in range(4))
 for i in range(8):cols.append(C.mul({4*MI[(i,)]:F(-1)}if i<7 and(i,)in MI else Q.scaled({4*j:c for j,c in gj.items()},-1)if i==7 else{},{1:F(1)}))
 assert len(cols)==56;body={}
 for j,c in cold.items():body=Q.plus(body,Q.scaled(cols[j],c))
 oldm0=[Q.sparse(v)for v in N['old_m0_24_by56']]
 for j,col in enumerate(cols):
  for s,mask in enumerate([1,2,3]):
   for i in range(7):
    if(i,)in MI:assert col.get(4*MI[(i,)]+mask,F(0))==oldm0[8*s+i].get(j,F(0))
 # A shift zi -> zi+c_i*r² cannot alter any zi*r coefficient; these
 # are all the canonical prefix derivatives selected by the exact dual.
 prefixes=[];oldterm=F(0)
 for row,c in oldadj.items():
  sec,i=divmod(row,8);assert i<7 and(i,7)in MI;v=body.get(4*MI[(i,7)]+[1,2,3][sec],F(0));prefixes.append([row,v,c]);oldterm+=c*v
 assert enc(prefixes)==V['all5_actual_old_prefix_r_values_with_adjoint']
 # Exact Res4(A^t) first derivatives, keeping ALL monic multiplication
 # carries. Its lower8 identity makes Z42_r literally zero.
 Z206=F(0)
 for j,c in cold.items():
  if j>=48:continue
  gen,mask=divmod(j,4);av=X['new_A12_by9_C36'][gen][5];z={k:F(av[k])for k in range(4)if F(av[k])}
  for v in range(8):
   if(v,)in MI:
    for k in range(4):Q.put(z,4*MI[(v,)]+k,F(av[4*(v+1)+k]))
  Z206+=c*C.mul(z,{mask:F(1)}).get(4*MI[(6,)],F(0))
 xp42=sum((a0.get(k,F(0))*raw[p].get(42,{}).get(8,F(0))for k,p in enumerate(xp)),F(0));xp20=sum((a0.get(k,F(0))*raw[p].get(20,{}).get(7,F(0))for k,p in enumerate(xp)),F(0));Jr=oldterm-xp42+F(20,3)*(Z206-xp20)
 assert Jr==F(V['actual_projected_natural_head4_first_target0_r']) and oldterm==F(V['actual_old_projected_first_contribution']) and Z206==F(V['actual_Z20_z6_cold0']) and xp42==F(V['actual_Xp42_r_a0']) and xp20==F(V['actual_Xp20_z6_a0'])
 labels=[tuple(v)for v in D['rv']['convolution_labels']];lower=[list(map(F,v))for v in D['rv']['all35_projected_M2_lower_representatives']];cor=[sum((c*lower[j][labels.index((f,7))]for j,c in lift.items()),F(0))for f in range(3)];wimage=[factors[f]*Jr-cor[f]for f in range(3)]
 assert cor==list(map(F,V['actual_R35_quadratic_lift_correction3'])) and wimage==list(map(F,V['actual_head4_z5_quadratic_image3']))
 # Rebind the two source images rather than treating saved coordinates
 # as anonymous vectors: head3 has zero origin, head4 uses its true origin.
 first3=[Q.pairs(v.items())for v in D['zero3']['actual_three_minimal_first_columns40']];rfirst={}
 for j,c in first3[2].items():t,i=divmod(j,8);Q.put(rfirst,36*t+im2[tuple(sorted((i,7)))],c)
 u=Q.mvec(pi,rfirst);q2=Q.mvec(pi,{36*t+im2[(7,7)]:c for t,c in origin.items()});assert [u.get(i,F(0))for i in range(3)]==list(map(F,V['known_actual_u_equals_e3r_image3'])) and [q2.get(i,F(0))for i in range(3)]==list(map(F,V['known_actual_q2_equals_e4r2_image3']))
 d=u[1]*q2[2]-u[2]*q2[1];assert d and d==F(V['fixed_pair_determinant']);inv=[[q2[2]/d,-q2[1]/d],[-u[2]/d,u[1]/d]];assert inv==[list(map(F,v))for v in V['fixed_pair_inverse2']]
 aa=inv[0][0]*wimage[1]+inv[0][1]*wimage[2];bb=inv[1][0]*wimage[1]+inv[1][1]*wimage[2];assert [aa,bb]==list(map(F,V['actual_w_equals_a_u_plus_b_q2_coefficients']))
 assert all(aa*u.get(i,F(0))+bb*q2.get(i,F(0))==wimage[i]for i in range(3))
 test=F(1)+aa*F(11767,3652);assert test and test==F(V['actual_qC_coefficient_of_r_kernel_vector'])
 assert D['fil']['remaining_actual_HF_options']==[[5,5,3,3],[5,5,3,2,1]]
 proof=dict(all48_actual_physical_root_values_and_nu_recomposed_from_scalar_DAG=True,source_DAG_ancestor_count=len(need),source_DAG_live_QQ_peak=peak,all56_old_C4_columns_and_all_available_m0_readbacks=True,actual_five_old_prefixes=prefixes,actual_Jr=Jr,actual_R35_M2_correction3=cor,actual_w_image3=wimage,actual_u_image3=u,actual_q2_image3=q2,source_pair_determinant=d,actual_a=aa,actual_b=bb,actual_nonzero_qC_coefficient=test,
  filtration_orientation='u=e3*r,w=e4*z5,q2=e4*r²,qC=e4*z5*r,g=e4*r³. t=w-a*u-b*q2 lies in F3, and r*t=(1+a*c)qC-(a*alpha+b)g. qC is nonzero modulo oldm3, hence rt is nonzero without evaluating alpha.',
  exact_conclusion='Natural injection sends rt to a nonzero newm4 class. Accepted newm5=0, r-reduction and exact length16/HF553 force completeHF55321.',actual_full_HF_Tnil=[5,5,3,2,1],actual_m_adic_lengths=[16,11,6,3,1,0],actual_m5_zero=True,actual_m4_dimension1=True,
  actual_deep_generator_old_expression='i(e4*z5*r - a*e3*r² - b*e4*r³)',deep_generator_nonzero_initial_old_qC_coefficient=test,old_alpha_not_evaluated_or_assumed=True,rA_not_DVRell=True,no_original_nilideal_cube_or_generators_claim=True,no_new_source_gauge_jet_or_solve=True)
 save('phase2-independent-physical-source2-map-and-complete-HF-proof.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/'code/ambient_B18_honest_nil_T_gauge120_review_20261105_v63.py',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,actual_full_HF_Tnil=[5,5,3,2,1],actual_total_length16=True,actual_m5_zero=True,actual_m4_dimension1=True,nonzero_qC_coefficient=str(test),full_source_old_n1_C1_R35_correction_and_all3_NF_coordinates=True,no_new_source_jet_or_original_nilcube_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20;save('report.json',rep);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
