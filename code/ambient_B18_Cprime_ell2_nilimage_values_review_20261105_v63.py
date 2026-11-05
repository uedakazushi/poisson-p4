#!/usr/bin/env python3
"""Independent PROVIDED finite33 gauge and actual nilimage quotient.
Sparse table arithmetic; no source, gauge or coefficient solve.
"""
import collections,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
from ambient_B18_honest_nil_T_gauge120_review_20261105_v63 import Cover,V
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def bits(p):
 assert all(max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=2048 for c in p.values());return p
class Base33:
 def __init__(self,A):self.table=[[Q.pairs(c)for c in row]for row in A['table_columns_for_each_left_basis']];self.deg=A['basis_degrees']
 def mul(self,x,y):
  z={}
  for i,c in x.items():
   for j,d in y.items():
    for k,a in self.table[i][j].items():Q.put(z,k,c*d*a)
  return bits(z)
 def inv(self,p):
  c=p[0];nil=Q.plus({0:F(1)},Q.scaled(p,-1/c));ans={0:F(1)};v=ans
  for _ in range(4):v=self.mul(v,nil);ans=Q.plus(ans,v)
  assert not self.mul(v,nil);ans=Q.scaled(ans,1/c);assert self.mul(p,ans)=={0:F(1)};return bits(ans)
class Ring132(Cover):
 def inv(self,x):
  c=x[(0,0)];h=self.scale(self.add(self.c(c),self.neg(x)),1/c);v=self.ONE;z=self.ONE
  for _ in range(10):v=self.mul(v,h);z=self.add(z,v)
  assert not self.mul(v,h);z=self.scale(z,1/c);assert self.mul(x,z)==self.ONE;self.inverse_count+=1;return z
 def install(self):
  super().install()
  def mi(A):
   n=len(A);C=V.constmat(V.rational_inverse([[x.get((0,0),F(0))for x in row]for row in A]));E=V.matmul(C,[[self.add(x,self.c(-x.get((0,0),F(0))))for x in row]for row in A]);term=V.eye(n);ans=V.eye(n)
   for k in range(1,11):term=V.matmul(term,E);ans=[[self.add(ans[i][j],self.scale(term[i][j],(-1)**k))for j in range(n)]for i in range(n)]
   assert all(not x for row in V.matmul(term,E)for x in row);ans=V.matmul(ans,C);assert V.matmul(A,ans)==V.matmul(ans,A)==V.eye(n);return ans
  V.matrix_inverse=mi
def apply(cols,v):
 z={}
 for j,c in v.items():z=Q.plus(z,Q.scaled(cols[j],c))
 return z
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63';rep=load(primary);assert pins[f'reports/{primary}/report.json']=='cffac50257d126e9f0b4cae908228ec7d1b8ea2f547d5765cf331bdfd6633782'and rep['all_passed']
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 assert load('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_SOURCE_20261105_v63_REVIEW')['source_frame_budget_gate_passes']
 phase1=load(primary,'phase1-source150-readback-fingerprints.json');g=load(primary,'phase2-canonical30-exact-values.json');proper=load(primary,'phase3-actual-kappa9-and-normalized-nilframe.json');result=load(primary,'phase4-actual-Cprime2-quotient-and-actions.json')
 names={
 'source':('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json'),
 'T':('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-source105-ell2-raw33-values.json'),
 'table':('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63','exact-existing-raw33-multiplication-table.json'),
 'interface':('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63','actual-source-physical-firstjet-interface.json'),
 'proper':('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63','literal-proper-source14-preflight.json'),
 'first':('AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63','actual-native128-first-input-and-mixed-index-profile.json'),
 'C2':('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-physical12-and-C2-right-exact-quotient.json'),
 'raw':('AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63','actual-ell2-product33-and-raw45-actions.json'),
 'B11':('AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63','actual-B11-source-and-gauge-preflight.json'),
 'B11gauge':('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase2-gauge30-B11-correction.json'),
 'B11proper':('AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63','phase3-actual-proper-G8-B11.json'),
 'bar':('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63','actual-physical12-actions-and-C-special.json'),
 'primebar':('AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63','actual-eight-nil-images-and-reduced-Cbar-cokernel.json')}
 D={k:load(*p)for k,p in names.items()};save('phase.json',dict(phase='immutable PROVIDED990 and actualsource33 inputs pinned; no solve'))
 A=D['table'];B=Base33(A);T=D['T'];leaves={i:Q.pairs(v)for i,v in zip(T['joint105_unknown_order'],T['all105_unknown_final_values'])};S=D['source'];nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values()for pair in f['complete_normalized_physical_w_coefficients'][0]for i in pair]+[f['actual_tracefree_nu']for f in S['families'].values()];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i)
  if i in leaves:continue
  n=nodes[i];assert n[0]in['a','m','n','i','c']
  if n[0]!='c':todo.extend(n[1:])
 uses=collections.Counter(roots)
 for i in need:
  if i not in leaves and nodes[i][0]!='c':uses.update(nodes[i][1:])
 vals={};peak=0
 for i in sorted(need):
  n=nodes[i];o=n[0]
  if i in leaves:p=leaves[i]
  elif o=='c':p={0:F(n[1],n[2])}if n[1]else{}
  elif o=='a':p=Q.plus(vals[n[1]],vals[n[2]])
  elif o=='m':p=B.mul(vals[n[1]],vals[n[2]])
  elif o=='n':p=Q.scaled(vals[n[1]],-1)
  else:p=B.inv(vals[n[1]])
  vals[i]=bits(p);peak=max(peak,sum(map(len,vals.values())));assert peak+100000<150000
  if i not in leaves and o!='c':
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j]
 pairs={f:[[vals[i]for i in pair]for pair in d['complete_normalized_physical_w_coefficients'][0]]for f,d in S['families'].items()};nu=[vals[d['actual_tracefree_nu']]for d in S['families'].values()]
 for f,p in pairs.items():
  dense=[[[x.get(k,F(0))for k in range(33)]for x in pair]for pair in p];assert hashlib.sha256(json.dumps(enc(dense),separators=(',',':')).encode()).hexdigest()==phase1['all600_physical_root_readbacks_sha256'][f]
 assert [[p.get(k,F(0))for k in range(33)]for p in nu]==[list(map(F,p))for p in phase1['nu_A_B33']]
 overlap=0
 for i,v in T['all_requested164_root_values']:
  if i in vals:assert vals[i]==Q.pairs(v);overlap+=1
 assert overlap==phase1['source164_exact_overlap_count'];del vals
 ring=Ring132(B,nu);ring.install();I,P=D['interface'],D['proper'];GL=D['first']['GL18_source_parameter_positions']
 def frame(fam):
  d=S['families'][fam];p=ring.add(ring.tau(0 if fam=='CA'else 1),ring.scale(ring.base(leaves[d['p_quadratic_relation']['s']]),F(1,2)));M=V.eye(5)
  for name,(x,y)in d['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':
    a,b=divmod(GL[int(name[5:])-1],5);M[a][b]=ring.add(M[a][b],ring.base(leaves[x]),ring.mul(ring.base(leaves[y]),p))
  C=V.eye(5)
  if fam=='CA':C[0][3]=ring.scale(p,-2);C[1][0]=p;C[1][3]=ring.neg(ring.mul(p,p))
  else:C[0][3]=ring.scale(p,F(1,2));C[1][0]=p;C[1][3]=ring.scale(ring.mul(p,p),F(1,4))
  return V.matmul(M,C)
 HA,HB=frame('CA'),frame('CB');PA=[ring.add(ring.base(a),ring.mul(ring.base(b),ring.tau(0)))for a,b in pairs['CA']];ev=[Q.sparse(v)for v in g['actual_even15_B33']];od=[Q.sparse(v)for v in g['actual_odd15_B33']];G=V.eye(5)
 for(i,j),x,y in zip(I['actual_group15_positions'],ev,od):G[i][j]=ring.add(G[i][j],ring.base(x),ring.mul(ring.base(y),ring.tau(0)))
 Fm=V.normalized_flag_frame(V.matmul(G,HA));assert Fm==V.eye(5)
 rows=P['zero_coordinate_rows36'];rr=P['CA_constant_left_inverse_rows19'];center=list(map(F,P['actual_CA_source19_center']));tr=V.matmul(Fm,G);moved=V.tensor_push(PA,tr);co=[r[0]for r in V.matmul(V.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in rr])]
 assert all(co[j]==ring.c(center[j])for j in P['CA_fiber10_fixed_coefficient_indices'])
 FB=V.normalized_flag_frame(V.matmul(G,HB));actual=[ring.add(co[j],ring.c(-center[j]))for j in P['CA_free_fiber9_coefficient_indices']]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j]for j in I['retained_source9_indices']]
 def dense(p):return [p.get((i,m),F(0))for i in range(33)for m in range(4)]
 assert [dense(p)for p in actual]==[list(map(F,p))for p in proper['actual14_C132']] and [dense(p)for p in kappa]==[list(map(F,p))for p in proper['kappa9_C132']]
 # All supplied990 values are tested via the ORIGINAL nonlinear30
 # canonical equations plus the same unit origin; no correction solve.
 base=A['basis33'];bm={tuple(m):i for i,m in enumerate(D['B11']['B11_basis_monomials8'])};maps=[]
 for z,r,e in base:
  mon=[0]*8;mon[7]=r
  if z:mon[z-1]+=1
  maps.append({bm[tuple(mon)]:F(1)}if not e and tuple(mon)in bm else{})
 def project11(p):return apply(maps,p)
 for given,old in[(ev,D['B11gauge']['even15_B11']),(od,D['B11gauge']['odd15_B11'])]:assert [project11(p)for p in given]==[Q.sparse(v)for v in old]
 for p,old in zip(kappa,D['B11proper']['kappa9_C_B44']):assert [project11({i:c for(i,m),c in p.items()if m==mask}).get(k,F(0))for k in range(11)for mask in range(4)]==list(map(F,old))
 av,d,x,t=[kappa[i]for i in[0,1,2,6]];h=ring.add(av,ring.neg(d),x,ring.scale(t,-2));c=ring.scale(ring.add(d,x),F(1,2));l=ring.scale(ring.add(av,d,ring.neg(x)),F(1,3));nil2=[ring.mul(d,h),ring.tm(d,c,l)];nil8=[ring.mul(p,{(0,mask):F(1)})for p in nil2 for mask in range(4)]
 assert [dense(p)for p in nil2]==[list(map(F,p))for p in proper['actual_nil2_C132']] and [dense(p)for p in nil8]==[list(map(F,p))for p in proper['actual_nil8_C132']]
 save('phase1-independent-gauge990-constraints-and-actual-nilframe.json',dict(original600physical_roots_and_source164_readbacks=True,all30_nonlinear_constraints=True,all_actual14_kappa9_and_nil8=True,allB11gauge_and_kappa_overlaps=True,trace_unit_inverse_uses0through10=True,inverse_products_verified=ring.inverse_count,no_canonical_or_source_solve=True))
 C2=D['C2'];q6=[[Q.pairs(v)for v in a['quotient_map6_by33_columns']]for a in D['raw']['canonical_branch_quotients6']];proj=[Q.pairs(v)for v in C2['all45_projection_columns']]
 def raw45(p):
  a=apply(q6[0],{i:c for(i,m),c in p.items()if m==1});b=apply(q6[1],{i:c for(i,m),c in p.items()if m==2});return Q.plus(a,{i+6:c for i,c in b.items()},{i+12:c for(i,m),c in p.items()if m==3})
 images=[apply(proj,raw45(p))for p in nil8];assert images==[Q.pairs(v)for v in result['all8_nil_images_C2']];acts=[[Q.pairs(v)for v in a]for a in C2['C2_coordinate_actions9_sparse_columns']]
 closure=[]
 for p in images:
  for z,r,e in base:
   v=apply(acts[z-1],p)if z else p
   for _ in range(r):v=apply(acts[7],v)
   for _ in range(e):v=apply(acts[8],v)
   closure.append(v)
 assert closure==[Q.pairs(v)for v in result['all264_R0_coefficient_closure_columns']]
 piv={}
 def reduce(v):
  for k,p in sorted(piv.items()):v=Q.plus(v,Q.scaled(p,-v.get(k,F(0))))
  return v
 for j in result['fixed_image_pivot_columns']:
  p=reduce(closure[j]);assert p;k=min(p);piv[k]=Q.scaled(p,1/p[k])
 assert len(piv)==7 and piv=={p:Q.pairs(v)for p,v in result['fixed_image_normalized_basis']}
 assert all(not reduce(p)for p in closure)and all(not reduce(apply(a,p))for a in acts for p in piv.values())
 free=[i for i in range(14)if i not in piv];assert free==result['quotient_original_C2_basis_indices']and len(free)==7
 def qp(p):
  p=reduce(p);return {i:p[k]for i,k in enumerate(free)if p.get(k,F(0))}
 qa=[[qp(a[j])for j in free]for a in acts];assert qa==[[Q.pairs(v)for v in a]for a in result['actual_Cprime2_coordinate_actions9']]
 assert all(apply(a,apply(b,{j:F(1)}))==apply(b,apply(a,{j:F(1)}))for a in qa for b in qa for j in range(7))and all(not apply(qa[8],v)for v in qa[8])
 old26=[Q.pairs(v)for v in D['bar']['all26_projection_columns']];to26=[Q.pairs(v)for v in C2['raw45_to_raw26_ell0_columns']];pb=D['primebar'];oldnil=[Q.sparse(v)for v in pb['all8_actual_nil_image_vectors_in_Cbar8']]
 assert [apply(old26,apply(to26,raw45(p)))for p in nil8]==oldnil
 bp={p:Q.sparse(v)for p,v in pb['normalized_triangular_image_basis']};bf=pb['quotient_original_basis_indices']
 def bproj(v):
  for k,p in sorted(bp.items()):v=Q.plus(v,Q.scaled(p,-v.get(k,F(0))))
  return {i:v[k]for i,k in enumerate(bf)if v.get(k,F(0))}
 qbar=[bproj(apply(old26,to26[C2['C2_basis_raw45_indices'][j]]))for j in free];assert qbar==[Q.pairs(v)for v in result['actual_Cprime2_to_Cprimebar4_columns']]
 ba=[[Q.sparse(v)for v in a]for a in pb['all8_quotient_coordinate_actions_columns']]+[[{}for _ in range(4)]];assert all(apply(qbar,apply(a,{j:F(1)}))==apply(b,apply(qbar,{j:F(1)}))for a,b in zip(qa,ba)for j in range(7))
 # Small literal QQ checks only; the actual264 image rank above used
 # the PROVIDED seven pivot columns, not a new source-column search.
 def rr(cols):
  ps={};selected=[]
  for j,v in enumerate(cols):
   p=v.copy()
   for k,a in sorted(ps.items()):p=Q.plus(p,Q.scaled(a,-p.get(k,F(0))))
   if p:k=min(p);ps[k]=Q.scaled(p,1/p[k]);selected.append(j)
  return ps,selected
 ellbasis,ellcols=rr(qa[8]);barbasis,barcols=rr(qbar);assert len(ellbasis)==3 and len(barbasis)==4
 proof=dict(all990_provided_gauge_values_verified_without_solving=True,all600physical_DAG_root_readbacks=True,all30_original_canonical_constraints=True,all14_kappa9_and8_nilframe_values=True,all_B11_and_ell0_nil8_overlaps=True,all264_closure_columns_and_fixed7_basis_verified=True,all9_action_invariances_and_Cprimebar_intertwiners=True,
  actual_Cprime2_length=7,actual_image_rank7=True,actual_ell_rank3=True,ell_fixed_independent_columns=ellcols,ell_normalized_basis=ellbasis,bar_fixed_independent_columns=barcols,
  actual_DVR_torsion_consequence='rank2 andtorsionmu2 imply exponent multiset{1,e},e>=2. The second exact exponent remains open; kerell/imell of Cprime2 sees only one genuine visible short-torsion direction.',
  genuine_target_not_identified=True,no_actual_connecting_values_or_Nbar14_lifts=True,no_ell_division_or_finite_H1_shortcut=True,no_new_source_or_gauge_solve=True)
 save('phase2-independent-nilimage-quotient-actions-and-scope.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/'code/ambient_B18_honest_nil_T_gauge120_review_20261105_v63.py',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 final=dict(id=ID,all_passed=True,frozen=True,actual_Cprime2_length=7,actual_image_rank=7,actual_ell_action_rank=3,all_source_gauge_nilimage_and_action_gates_verified=True,genuine_ell_torsion_target_still_open=True,no_new_source_jet_or_solve=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20;save('report.json',final);print(json.dumps({k:v for k,v in final.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
