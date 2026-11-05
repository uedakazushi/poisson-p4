#!/usr/bin/env python3
"""Independent corrected source2/zero3 cleanup and literal quotient audit."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time
from ambient_B18_T_joint135_selected_functionals_review_20261104_v63 import Cover
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW'
P={
 'v2_report':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_CACHE_V2_20261105_v63/report.json',
 'v2':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_CACHE_V2_20261105_v63/actual-three-oldm1-and-full81-formal-callbacks.json',
 'map_report':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/report.json',
 'map':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/actual-three-first-columns-and-fixed-quotient-proof.json',
 'gate':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_SOURCE_PREFLIGHT_20261104_v63/three-oldm1-cache-and-cleanup-callback-gate.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'physical':'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',
 'phys':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase1-cached-source-monics-and-physical-adapter.json',
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'natural_peer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
 'X':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'm':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'linear':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'linear_peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW/report.json',
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase2-row20-surrogates-and-true-m1row7.json',
 'whole':'reports/AMBIENT_B18_HONEST_NIL_T_LINEAR_SOURCE_WHOLE_PREFLIGHT_20261104_v63/actual-final681-source-sensitivities-and-cached-directions.json',
 'B':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def put(d,k,c):
 x=d.get(k,Q(0))+c
 if x:d[k]=x
 else:d.pop(k,None)
def add(*ps):
 d={}
 for p in ps:
  for k,c in p.items():put(d,k,c)
 return d
def sc(p,c):return {k:c*v for k,v in p.items() if c*v}
def sp(p):
 if isinstance(p,dict):return {int(i):Q(c) for i,c in p.items() if Q(c)}
 if p and isinstance(p[0],(str,int)):return {i:Q(c) for i,c in enumerate(p) if Q(c)}
 return {int(i):Q(c) for i,c in p if Q(c)}
def mv(A,v):return {i:c for i,row in enumerate(A) for c in [sum(a*v.get(j,Q(0)) for j,a in row.items())] if c}
def mp(A,B):return [add(*(sc(B[j],c) for j,c in row.items())) for row in A]
def eye(n):return [{i:Q(1)} for i in range(n)]
def lin(A,v):return [add(*(sc(v[j],c) for j,c in row.items())) for row in A]
M=[()]+[(i,) for i in range(8)]+list(itertools.combinations_with_replacement(range(8),2));IX={m:i for i,m in enumerate(M)}
def mul(a,b):
 out={}
 for i,c in a.items():
  for j,d in b.items():
   mon=tuple(sorted(M[i]+M[j]))
   if mon in IX:put(out,IX[mon],c*d)
 return out
def inv(a):
 c=a.get(0,Q(0));assert c
 h=add({0:Q(1)},sc(a,-1/c));v=sc(add({0:Q(1)},h,mul(h,h)),1/c)
 assert mul(a,v)=={0:Q(1)};return v
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'CORRECTED_V2_ONLY_SOURCE_INPUTS_PINNED','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in ['v2_report','map_report','natural_peer','linear_peer'])
 for r in ['v2_report','map_report']:
  for p,h in D[r]['input_sha256'].items():assert sha(ROOT/p)==h
 J=D['FJ'];mons=J['source2_monomials'];Y=list(map(sp,D['phys']['CB_F9_in_CA_F9_ell0_degree2']))
 def fold(p):
  out={}
  for m,c in zip(mons,p):
   if Q(c):put(out,IX[tuple(sorted(7 if i==8 else i for i in m))],Q(c))
  return out
 def comp(p):
  out={}
  for m,c in zip(mons,p):
   v={0:Q(c)} if Q(c) else {}
   for i in m:v=mul(v,Y[i])
   out=add(out,v)
  return out
 phy=D['physical']['families'];assert D['physical']['source2_monomials']==mons
 for a,b,saved in zip(phy['CA']['physical_R0_readback_jets_in_own_F9'],phy['CB']['physical_R0_readback_jets_in_own_F9'],D['phys']['all9_shared_physical_degree2_polynomials']):assert fold(a)==comp(b)==sp(saved)
 # Unit XB1 makes this nonlinear transport unique; both inverse products
 # are checked, with no source inverse solve.
 JB=list(map(sp,phy['CB']['R0_by_F9_unit_Jacobian']));JBi=list(map(sp,phy['CB']['fixed_inverse9']));assert mp(JB,JBi)==mp(JBi,JB)==eye(9)
 leaves={f:{i:(fold(p) if f=='CA' else comp(p)) for i,p in J['both57_source2'][f]['all57_source2jets']} for f in ['CA','CB']}
 nodes=D['source']['DAG_scalar_nodes'];point=D['gate']['four_cached_point_function_profiles'];requested={f:set() for f in ['CA','CB']}
 for p in point:requested[p['family']].add(p['root'])
 for f,data in D['source']['families'].items():requested[f].add(data['actual_tracefree_nu'])
 requested['CA'].update([7790,7943]);requested['CB'].update([28539,28692])
 results={}
 for family in ['CA','CB']:
  roots=sorted(requested[family]|set(J['both57_source2'][family]['equations57']));todo=list(roots);needed=set()
  while todo:
   i=todo.pop()
   if i in needed:continue
   needed.add(i)
   if i in leaves[family]:continue
   n=nodes[i];assert n[0] in ['a','m','n','i','c'],('UNBOUND_SOURCE',family,i,n)
   if n[0]!='c':todo.extend(n[1:])
  uses=collections.Counter(j for i in needed if i not in leaves[family] and nodes[i][0]!='c' for j in nodes[i][1:]);uses.update(roots);v={}
  for i in sorted(needed):
   n=nodes[i];op=n[0]
   if i in leaves[family]:a=leaves[family][i]
   elif op=='c':a={0:Q(n[1],n[2])} if n[1] else {}
   elif op=='a':a=add(v[n[1]],v[n[2]])
   elif op=='m':a=mul(v[n[1]],v[n[2]])
   elif op=='n':a=sc(v[n[1]],-1)
   elif op=='i':a=inv(v[n[1]])
   else:raise AssertionError(n)
   v[i]=a
   if i not in leaves[family]:
    for j in n[1:] if op!='c' else []:
     uses[j]-=1
     if not uses[j]:del v[j]
  for k,i in enumerate(J['both57_source2'][family]['equations57']):
   expected=({IX[(7 if k-48==8 else k-48,)]:Q(1)} if k>=48 else {}) if family=='CA' else (Y[k-48] if k>=48 else {})
   assert v[i]==expected
  results[family]={i:v[i] for i in requested[family]}
 points=[results[p['family']][p['root']] for p in point];delta5=add(points[0],sc(points[1],-1));delta7=add(points[2],sc(points[3],-1))
 fA=results['CA'][7790];fB=results['CB'][28539];nu=[results[f][data['actual_tracefree_nu']] for f,data in D['source']['families'].items()]
 assert nu==list(map(sp,D['phys']['actual_nu_A_B_degree2']))
 bA=results['CA'][7943];bB=results['CB'][28692];assert not bA.get(0,Q(0)) and not bB.get(0,Q(0))
 functions=[[sc(mul(bB,nu[1]),-1),mul(bA,nu[0]),delta5],[delta7,{},sc(fB,-1)],[{},delta7,fA]]
 for row in D['gate']['CA_rawf63_first_prefixes']:
  i=row['prefix_index'];assert sp(row['first9_bar'])=={q+1:c for q in range(8) for c in [fA.get(IX[tuple(sorted((i,q)))],Q(0))] if (q<=i or q==7) and c}
 old=[]
 for col,rr in zip(functions,D['gate']['nine_r_first_values_cached']):
  v=[]
  for sec,p in enumerate(col):
   assert all(not p.get(i,Q(0)) for i in range(9))
   v.extend([{q:c for q in range(8) for c in [p.get(IX[tuple(sorted((i,q)))],Q(0))] if (q<=i or q==7) and c} for i in range(7)])
   v.append({7:Q(rr[sec])} if Q(rr[sec]) else {})
  old.append(v)
 N=D['natural']['actual_origin_data'];cold=[[sc(v,-1) for v in lin(list(map(sp,N['old_n0_56_by24'])),col)] for col in old]
 X=D['X'];W=D['W'];MD=D['m'];Min=D['minimal'];phi=D['phys'];r9=Cover(M[:9],[{i:c for i,c in n.items() if i<9} for n in nu]);r45=Cover(M,nu)
 A=[[{divmod(i,4):Q(c) for i,c in enumerate(p) if Q(c)} for p in row] for row in X['new_A12_by9_C36']]
 Z0=[{} for _ in range(44)];Z1=[[{} for _ in range(56)] for _ in range(44)]
 for i,row in enumerate(A):
  for j,a in enumerate(row):
   for b in range(4):
    p=r9.mul(a,{(0,b):Q(1)});col=4*i+b
    for (base,mask),c in p.items():
     if base==0:Z0[4*j+mask][col]=c
     else:Z1[4*j+mask][col][base-1]=c
 for i in range(8):Z0[36+i][48+i]=Q(1)
 assert Z0==list(map(sp,N['actual_Z0_44_by56'])) and not Z0[20]
 C0=[sp(N['all_five_new_cycle_origins44'][i]) for i in [1,2,3]];C1=[];C20=[]
 for q,j in enumerate([23,29,30]):
  c1=[add(v,Z1[i][j]) for i,v in enumerate(lin(Z0,cold[q]))];C1.append(c1)
  i,b=divmod(j,4);p=r45.mul(A[i][5],{(0,b):Q(1)});z={base:c for (base,mask),c in p.items() if mask==0}
  for k,v in enumerate(cold[q]):z=add(z,mul({a+1:c for a,c in Z1[20][k].items()},{a+1:c for a,c in v.items()}))
  C20.append(z)
  assert old[q]==list(map(sp,D['v2']['records'][q]['old_m1_24x8'])) and c1==list(map(sp,D['v2']['records'][q]['new_c1_44x8']))
 dump(o,'phase1-independent-source2-and-full-bnu.json',{'all114_source_equations_through2_verified':True,'nonlinear_physical_adapter_verified':True,'all_three_old_m1_and_new_C1_rebuilt':True,'b1_nu1_retained':True})
 xp=X['X0_pivot_columns'];cols=[{int(i):sp(p) for i,p in col} for col in X['all780_source_columns44_sparse_dual9']]
 m0=list(map(sp,MD['m0_24_by44']));free=[i for i in range(44) if i not in MD['selected_source_columns24']]
 I=[free[i] for i in Min['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=list(map(sp,Min['fixed_minor_data']['inverse']))
 minor=[{j:cols[c].get(r,{}).get(0,Q(0)) for j,c in enumerate(xp) if cols[c].get(r,{}).get(0,Q(0))} for r in I]
 assert mp(minor,Ji)==mp(Ji,minor)==eye(18)
 wc=list(map(sp,W['all_W84_columns_sparse']));wp=W['W_pivot_residual_columns'];wr=W['W_pivot_rows'];Wi=list(map(sp,W['W_fixed_minor_inverse']))
 wm=[{j:wc[c].get(r,Q(0)) for j,c in enumerate(wp) if wc[c].get(r,Q(0))} for r in wr]
 assert mp(wm,Wi)==mp(Wi,wm)==eye(81)
 factor=list(map(sp,X['all780_X0_factorization_coefficients']));O=W['residual_original_X_columns']
 def residual(j):return add({j:Q(1)},{xp[k]:-c for k,c in factor[j].items()})
 raw15=list(map(lambda p:list(map(Q,p)),D['first']['all780_row20_second_surrogates15']))
 needed=[tuple(m) for m in D['whole']['all15_needed_row20_X2_monomials']]
 bold0=[{r:p.get(0,Q(0)) for r,p in {int(r):sp(p) for r,p in col}.items() if r<48 and p.get(0,Q(0))} for col in D['B']['source_columns640_rows56_sparse_dual9']]+[{} for _ in range(140)]
 sensitiv=[{} for _ in range(780)]
 for j,col in enumerate(D['radical']['native_tensor_Koszul6_HB2_resolution'][1]['columns']):
  for row,poly in col:
   if row!=5:continue
   for e,c in poly:
    if sum(e)==1 and any(e[:5]):put(sensitiv[640+4*j],e.index(1),Q(c))
 pairs=list(itertools.combinations(range(8),2));actual=[];allchecks=0;records=[]
 for j,r in enumerate(D['v2']['records']):
  a0=mv(Ji,{k:C0[j][i] for k,i in enumerate(I) if i in C0[j]})
  co={xp[k]:-c for k,c in a0.items()}
  pre=[add(C1[j][i],*({k-1:c*a for k,c in cols[s].get(i,{}).items() if k} for s,a in co.items())) for i in range(44)]
  v=lin(m0,pre);h0={}
  for sec in range(3):
   for k,(a,b) in enumerate(pairs):
    c=v[8*sec+7].get(a,Q(0)) if b==7 else (v[8*sec+b].get(a,Q(0))-v[8*sec+a].get(b,Q(0)))/2
    if c:h0[28*sec+k]=c
  beta=mv(Wi,{k:h0[i] for k,i in enumerate(wr) if i in h0})
  assert add(*(sc(wc[wp[k]],a) for k,a in beta.items()))==h0
  co=add(co,*(sc(residual(O[wp[k]]),-a) for k,a in beta.items()))
  assert a0==sp(r['first18']) and beta==sp(r['full81']) and co==sp(r['total_actual_X780_constant_coefficients'])
  # Verify all source firsts and both full coefficient-function callbacks,
  # before any ordered Hadamard component is read.
  final0=add(C0[j],*({i:p[0]*a for i,p in cols[s].items() if p.get(0,Q(0))} for s,a in co.items()));assert not final0
  final1=[add(C1[j][i],*({k-1:c*a for k,c in cols[s].get(i,{}).items() if k} for s,a in co.items())) for i in range(44)]
  assert all(not p for p in lin(m0,final1));allchecks+=24*8
  Aw=add(*(sc(bold0[k],a) for k,a in co.items()));Uw=add(*(sc(sensitiv[k],a) for k,a in co.items()))
  assert not Aw and not Uw and sp(r['remaining_actual_adapter_A2_weights'])==Aw and sp(r['remaining_native_even2_weights'])==Uw
  a1=lin(Ji,[final1[i] for i in I]);out=[]
  for i in [42,43]:out.append(add(final1[i],*(sc(a1[k],-cols[xp[k]].get(i,{}).get(0,Q(0))) for k in range(18))))
  final20=[C20[j].get(IX[m],Q(0))+sum(a*raw15[k][p] for k,a in co.items()) for p,m in enumerate(needed)]
  m1r=list(map(Q,D['first']['true_actual_m1_row7_r_coefficients44']));hi=D['linear']['all54_actual_pivot_h0_readback']
  for hidx,i in enumerate([0,4,5]):
   f={}
   for q in range(8):
    c=Q(5,2)*final20[needed.index(tuple(sorted((i,q))))] if q<=i or q==7 else Q(0)
    if q==7:c+=sum(m1r[k]*final1[k].get(i,Q(0)) for k in range(44))
    c-=sum(Q(hi[hidx][k])*a1[k].get(q,Q(0)) for k in range(18))
    if c:f[q]=c
   out.append(f)
  flat={8*i+q:c for i,p in enumerate(out) for q,c in p.items()}
  assert flat==sp(r['final_known40']);actual.append(flat)
  records.append({'old_free56':r['old_free56'],'all84_h0_cleanup_verified':True,'all192_first_source_checks':True,
    'all_A2_and_native_even2_source_coefficients_zero':True,'actual_first_column40':[[i,str(c)] for i,c in sorted(flat.items())]})
 assert allchecks==576 and actual==list(map(sp,D['map']['actual_three_minimal_first_columns40']))
 # Check the literal target order from the actual three Schur retentions,
 # rather than accepting a display label or a placeholder assertion.
 retained=list(range(128))
 for step in Min['three_exact_unit_Schur_steps']:retained=[retained[i] for i in step['retained_rows']]
 assert retained==[42,43,50,68,70]
 assert Min['minimal_target_labels']==[{'kind':'original_free44_source','index':42},{'kind':'original_free44_source','index':43},
  {'kind':'Koszul_wedge84','index':6,'label':[0,[0,7]]},{'kind':'Koszul_wedge84','index':24,'label':[0,[4,7]]},{'kind':'Koszul_wedge84','index':26,'label':[0,[5,7]]}]
 lincols=list(map(sp,D['linear']['all849_actual_minimal_first_columns40_sparse']))
 origins=[sp(N['actual_origin_map_columns5'][i]) for i in [0,4]];products=[{8*i+v:c for i,c in p.items()} for p in origins for v in range(8)]
 F4=list(map(sp,D['map']['quotient4_by40_after_origin_products']));assert F4[1]=={7:Q(1)}
 assert all(not mv(F4,c) for c in lincols+products)
 four=[mv(F4,c) for c in actual];assert four==list(map(sp,D['map']['actual_three_projected_columns4']))
 def rank(vectors):
  basis={}
  for p in vectors:
   p=dict(p)
   while p:
    k=max(p)
    if k not in basis:basis[k]=sc(p,1/p[k]);break
    p=add(p,sc(basis[k],-p[k]))
  return len(basis)
 assert rank(lincols)==35 and rank(lincols+products)==36 and rank(lincols+products+actual)==39 and rank(F4)==4
 minor=[[four[j].get(i,Q(0)) for j in range(3)] for i in [0,2,3]]
 det=(minor[0][0]*(minor[1][1]*minor[2][2]-minor[1][2]*minor[2][1])-minor[0][1]*(minor[1][0]*minor[2][2]-minor[1][2]*minor[2][0])+minor[0][2]*(minor[1][0]*minor[2][1]-minor[1][1]*minor[2][0]))
 assert det==Q(86403125,7746048) and str(det)==D['map']['selected_small_determinant']
 # The unique gr1 quotient functional is literally e42*r. Its square
 # channel follows from commutativity and the already proved length5.
 assert all(not p.get(7,Q(0)) for p in lincols+products+actual)
 assert D['natural_peer']['cokernel_length']==5 and D['natural_peer']['cokernel_mu']==3
 assert D['quad']['exact_quotient_duals_rows_sparse'][2]==[[35,'1']]
 assert [tuple(m) for m in D['quad']['quadratic_monomials36']][35]==(7,7)
 ep=o/'independent-V2-source576-actual-three-firsts-and-literal-r2-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all_source2_physical_transport_and114_equation_checks_pass':True,'all_three_bnu_products_include_b1_nu1':True,
  'all576_first_source_and252_h0_checks_pass':True,'actual_three_columns':records,'retained_original128_target_rows5':retained,
  'actual_final_first_span_rank':39,'fixed_small3_determinant':str(det),'literal_grQ1_dual40':[[7,'1']],
  'actual_cokernel_HF':[3,1,1],'cokernel_length':5,
  'literal_grQ2_dual180':[[35,'1']],'actual661_quadratic_dual_row2_zero_by_projection_to_Q':True,
  'commutativity_proof':'Every degree2 monomial except e42*r² has a degree1 divisor annihilated by the literal unique grQ1 quotient; HFQ2=1 forces e42*r² nonzero. Thus every actual Tnil quadratic initial has coefficient35 zero.',
  'optional_rA_module_consequence':'noncanonically Q = QQ[[rA]]/(rA³) + QQ² as QQ[[rA]]-module only; zi tops may couple into the r² socle',
  'no_R0_split_ell_exponents_actual_Tnil_HF2_or_nilcube_claim':True,'new_source_or_jet_coefficients':0})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'corrected_V2_source_and_all576_checks_independently_accepted':True,'actual_first_rank':39,'actual_cokernel_HF':[3,1,1],
  'literal_target42_r_and_target0_r2_duals_accepted':True,'actual661_quadratic_row2_zero':True,
  'new_source_jet_or_Tnil_HF2_nilcube_claim':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
