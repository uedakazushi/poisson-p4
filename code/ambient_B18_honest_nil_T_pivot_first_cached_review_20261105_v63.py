#!/usr/bin/env python3
"""Independent PROVIDED A via N*A=E8*D and all81 source carries.
No inverse solve, gauge solve, source jet or producer imports.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 r=load(BASE);assert pins[f'reports/{BASE}/report.json']=='00342a7baf3e0b81f68e0e1d4d4e5c47ee5103ae02a27ecc2f0b229178c3fc0b'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 supplied=load(BASE,'actual81-pivot-first-carries-and-six-H1-source-adjoints.json')
 assert load('AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW')['all_passed']
 names={
 'N2':('AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63','phase2-actual-N-B11.json'),
 'B11':('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase3-actual-proper-G8-B11.json'),
 'nu':('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase1-source-B11-values.json'),
 'basis':('AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63','actual-B11-source-and-gauge-preflight.json'),
 'radical':('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json'),
 'GL':('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63','actual-E4-nil-two-generators-and-relations.json'),
 'X1':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json'),
 'm0':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json'),
 'W':('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase3-actual-W84-rank-and-mu.json'),
 'minimal':('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json'),
 'first2':('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase2-row20-surrogates-and-true-m1row7.json'),
 'first3':('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json'),
 'oldB1':('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','actual-first640-X-mod-m2.json'),
 'd3':('AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_VALUES_20261104_v63','actual153-d3-quadratic-quotient-values-and-known-rank.json'),
 'd3source':('AMBIENT_B18_HONEST_NIL_T_D3_QUADRATIC_SOURCE_PREFLIGHT_20261104_v63','actual153-d3-quadratic-adjoints-and12-callback-sensitivities.json'),
 'quad':('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63','known-quadratic-closure-and-complete814-source-projection.json'),
 'Koszul':('AMBIENT_B18_REDUCED_CBAR_KOSZUL_INPUT_20261104_v63','actual-reduced-Cbar-three-sector-Koszul.json')}
 D={k:load(*p)for k,p in names.items()};save('phase.json',dict(phase='all source/frame inputs pinned; no solve'))
 C=Q.FiniteC(D['basis']['B11_basis_monomials8'],D['nu']['nu_A_B_B11']);R=D['radical'];P=D['B11']
 core=[Q.sparse(P['kappa9_C_B44'][i])for i in[0,1,2,6]];vs=[Q.sparse(v)for v in P['actualG8_C_B44'][:5]]+[C.polynomial(D['GL']['forward_GL4'][s],core)for s in['h','c','r','l']]
 N=[[Q.sparse(v)for v in row]for row in D['N2']['N_C_B44']];A=[[Q.sparse(v)for v in row]for row in supplied['actual_A12_by9_C44']]
 DD=[[C.polynomial(p,vs)for p in row]for row in R['G8_polynomial_factor_D8_by9']]
 assert Q.det_fixed([{j:p.get(0,F(0))for j,p in enumerate(row)if p.get(0,F(0))}for row in N])
 for i in range(12):
  for j in range(9):assert Q.plus(*(C.mul(N[i][k],A[k][j])for k in range(12)))==(DD[i][j]if i<8 else{})
 AF=[C.component(C.mul(row[5],{mask:F(1)}),0)for row in A for mask in range(4)]
 assert [[v.get(i,F(0))for i in range(11)]for v in AF]==[[F(c)for c in row]for row in supplied['all48_actual_ResAT_row20_B11']]
 assert all(not v.get(0,F(0))for v in AF)
 mons15=[tuple(m)for m in D['first2']['monomials15']];raw=[]
 for j,col in enumerate(D['oldB1']['source_columns640_rows56_sparse_dual9']):
  b0={row:Q.pairs(sp).get(0,F(0))for row,sp in col if row<48}
  raw.append([F(D['first2']['all780_row20_second_surrogates15'][j][mons15.index((i,7))])+sum((a*AF[row].get(bi,F(0))for row,a in b0.items()),F(0))for i,bi in zip([0,4,5],[7,9,10])])
 native=[]
 for col in R['native_tensor_Koszul6_HB2_resolution'][1]['columns']:
  poly=next((p for row,p in col if row==5),[]);v=C.polynomial(poly,vs)
  for mask in range(4):
   p=C.component(C.mul(v,{mask:F(1)}),0);native.append([p.get(b,F(0))for b in[7,9,10]])
 assert native==[[F(c)for c in row]for row in D['d3']['all140_actual_native_row20_z0r_z4r_z5r']]
 raw+=native
 X,W,M=D['X1'],D['W'],D['minimal'];xs=[{int(i):Q.pairs(v)for i,v in col}for col in X['all780_source_columns44_sparse_dual9']];xp=X['X0_pivot_columns'];xc=[Q.sparse(v)for v in X['all780_X0_factorization_coefficients']]
 Crows=[i for i in range(44)if i not in D['m0']['selected_source_columns24']];I=[Crows[i]for i in M['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in M['fixed_minor_data']['inverse']]
 hi=[[F(c)for c in row]for row in D['first3']['all54_actual_pivot_h0_readback']];mu1=list(map(F,D['first2']['true_actual_m1_row7_r_coefficients44']));O=W['residual_original_X_columns'];Wp=W['W_pivot_residual_columns'];Wr=W['W_pivot_rows'];Wi=[Q.sparse(v)for v in W['W_fixed_minor_inverse']]
 C1=[]
 for k in Wp:
  j=O[k];co=xc[j]
  x1=[Q.plus({v-1:a for v,a in xs[j].get(i,{}).items()if v},*(Q.scaled({v-1:a for v,a in xs[xp[p]].get(i,{}).items()if v},-c)for p,c in co.items()))for i in range(44)]
  a1=Q.mprod(Ji,[x1[i]for i in I]);top=[Q.plus(x1[i],*(Q.scaled(a1[p],-xs[xp[p]].get(i,{}).get(0,F(0)))for p in range(18)))for i in[42,43]];bottom=[]
  for n,i in enumerate([0,4,5]):
   r=raw[j][n]-sum((c*raw[xp[p]][n]for p,c in co.items()),F(0))
   bottom.append(F(5,2)*r+sum((mu1[p]*x1[p].get(i,F(0))for p in range(44)),F(0))-sum((hi[n][p]*a1[p].get(7,F(0))for p in range(18)),F(0)))
  C1.append([[v.get(q,F(0))for q in range(8)]for v in top]+[bottom])
 assert C1==[[[F(c)for c in p]for p in row]for row in supplied['all81_actual_visible_C1_top2_by8_and_wedge3_r']]
 Bg=[[F(c)for c in row]for row in D['d3source']['actual_Bold0_residual_graph4_by81']]
 for k in range(81):
  for n,bi in enumerate([7,9,10]):assert C1[k][2][n]-F(D['d3']['all81_pivot_C1_visible_surrogates5_by8'][k][2+n][7])==F(5,2)*sum((Bg[g][k]*AF[row].get(bi,F(0))for g,row in enumerate([32,36,40,44])),F(0))
 adj=[[{}for _ in range(81)]for f in range(6)];mons=list(itertools.combinations_with_replacement(range(8),2))
 for f,sp in enumerate(D['quad']['exact_quotient_duals_rows_sparse']):
  for index,c in sp:
   target,m=divmod(index,36);p,q=mons[m];c=F(c)
   for k in range(81):
    if target>=2:assert p==q==7;Q.put(adj[f][k],7,c*C1[k][2][target-2])
    else:
     Q.put(adj[f][k],q,c*C1[k][target][p])
     if p!=q:Q.put(adj[f][k],p,c*C1[k][target][q])
 assert [[[v.get(i,F(0))for i in range(8)]for v in row]for row in adj]==[[[F(c)for c in p]for p in row]for row in supplied['all6_actual_C1_adjoint_81_by8']]
 eff=[]
 for f in range(6):
  z=[]
  for h,row in enumerate(Wr):
   for var in range(8):
    c=sum((adj[f][k].get(var,F(0))*Wi[k].get(h,F(0))for k in range(81)),F(0))
    if c:z.append([row,var,c])
  eff.append(z)
 assert eff==[[[h,v,F(c)]for h,v,c in row]for row in supplied['all6_actual_effective_H1_callbacks_row84_var8_sparse']]
 projected=[]
 for col in D['Koszul']['all8_differentials_sparse_columns'][2]:
  u={}
  for row,mon,c in col:
   if sum(mon)==1 and row in Wr:
    v=mon.index(1);h=Wr.index(row)
    for k in range(81):Q.put(u,(k,v),F(c)*Wi[k].get(h,F(0)))
  projected.append([-sum((adj[f][k].get(v,F(0))*c for(k,v),c in u.items()),F(0))for f in range(6)])
 for comb in D['d3source']['all153_constant_d3_source_combinations']:assert all(sum((F(c)*projected[j][f]for j,c in comb),F(0))==0 for f in range(6))
 save('independent-NA-constraint-all81-carries-and-six-source-adjoints.json',dict(provided_A_unique_by_NA_equals_E8D=True,all1296_C44_equations=True,N_origin_unit=True,
  all48_actual_ResAT_row20_B11=True,own_direct_native140_three_callbacks=True,all81_actual_C1_top2_all8_bottom3_r=True,
  all243_graph_delta_checks=True,all6_actual_C1_adjoint81x8_and_effectiveH1_callbacks=True,all153_d3_six_projections_zero=True,
  no_unseen_C1_components_assigned_zero=True,no_new_source_jet_or_inverse_solve=True))
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,actual_A_and_all81_visible_C1_source_verified=True,all6_H1_adjoints_and153_d3_zero_verified=True,no_remaining661_scalar_values_or_rank_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20
 save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
