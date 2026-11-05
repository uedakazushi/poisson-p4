#!/usr/bin/env python3
"""Independent fixed Cprime quadratic source/value peer; no coefficient solve."""
import collections,gc,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
from ambient_B18_honest_nil_T_gauge120_review_20261105_v63 import Base,V
from ambient_B18_Cprime_ell2_nilimage_values_review_20261105_v63 import Ring132
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def bits(p):
 assert all(max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=2048 for c in p.values());return p
class B27(Base):
 def mul(self,a,b):return bits(super().mul(a,b))
class R108(Ring132):
 def mul(self,a,b):return bits(super().mul(a,b))
def source_inputs(D,phase):
 B=B27(phase['base27_monomials']);J=D['subset'];assert len(B.mon)==27
 y=[{B.idx[(i,)]:F(1)}for i in range(9)];YA=y[:8]+[Q.plus(y[7],Q.scaled(y[8],F(1,5)))];YB=[]
 def compose(p,args):
  v={}
  for mon,c in zip(J['source2_monomials'],p):
   if not F(c):continue
   z={0:F(c)}
   for i in mon:z=B.mul(z,args[i])
   v=Q.plus(v,z)
  return bits(v)
 for lin,quad in zip(D['adapter']['YB_in_YA_linear'],D['adapter']['YB_in_YA_quadratic']):
  v={}
  for i,c in enumerate(lin):
   if F(c):v=Q.plus(v,Q.scaled(YA[i],F(c)))
  for exp,c in quad:
   if not F(c):continue
   z={0:F(c)}
   for i,n in enumerate(exp):
    for _ in range(n):z=B.mul(z,YA[i])
   v=Q.plus(v,z)
  YB.append(v)
 leaves={i:compose(p,YA)for i,p in J['both57_source2']['CA']['all57_source2jets']};shared=[]
 for i,p in J['both57_source2']['CB']['all57_source2jets']:
  z=compose(p,YB)
  if i in leaves:assert leaves[i]==z;shared.append(i)
  else:leaves[i]=z
 assert sorted(shared)==list(range(2,11))and len(leaves)==105
 assert leaves=={i:Q.sparse(v)for i,v in phase['all105_cached_composite_values']}
 maps=[Q.pairs(v)for v in D['gate']['actual33_to26_columns']]
 def old26(v):
  z={}
  for i,c in v.items():z=Q.plus(z,Q.scaled(maps[i],c))
  return z
 old=D['oldsource'];old105={i:Q.pairs(v)for i,v in zip(old['joint105_unknown_order'],old['all105_unknown_final_values'])}
 assert all({j:c for j,c in v.items()if j<26}==old26(old105[i])for i,v in leaves.items())
 S=D['source'];nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values()for pair in f['complete_normalized_physical_w_coefficients'][0]for i in pair]+[f['actual_tracefree_nu']for f in S['families'].values()]+[i for f in ['CA','CB']for i in J['both57_source2'][f]['equations57']]+[7790,28539];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i)
  if i in leaves:continue
  n=nodes[i];assert n[0]in ['c','a','m','n','i']
  if n[0]!='c':todo.extend(n[1:])
 assert len(need)==20383;uses=collections.Counter(roots)
 for i in need:
  if i not in leaves and nodes[i][0]!='c':uses.update(nodes[i][1:])
 vals={};peak=0
 for i in sorted(need):
  n=nodes[i];op=n[0]
  if i in leaves:v=leaves[i]
  elif op=='c':v={0:F(n[1],n[2])}if n[1]else{}
  elif op=='a':v=Q.plus(vals[n[1]],vals[n[2]])
  elif op=='m':v=B.mul(vals[n[1]],vals[n[2]])
  elif op=='n':v=Q.scaled(vals[n[1]],-1)
  else:v=B.inv(vals[n[1]])
  vals[i]=bits(v);peak=max(peak,sum(map(len,vals.values())));assert peak+100000<150000
  if i not in leaves and op!='c':
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j]
 for f,args in [('CA',YA),('CB',YB)]:
  ev=[vals[i]for i in J['both57_source2'][f]['equations57']];assert ev[:48]==[{}]*48 and ev[48:]==args
 assert all({j:c for j,c in vals[i].items()if j<26}==old26(Q.pairs(v))for i,v in old['all_requested164_root_values'])
 pairs={f:[[vals[a],vals[b]]for a,b in data['complete_normalized_physical_w_coefficients'][0]]for f,data in S['families'].items()};nu=[vals[data['actual_tracefree_nu']]for data in S['families'].values()]
 for f,p in pairs.items():
  dense=[[[z.get(i,F(0))for i in range(27)]for z in pair]for pair in p];assert hashlib.sha256(json.dumps(enc(dense),separators=(',',':')).encode()).hexdigest()==phase['all600_physical_root_fingerprints'][f]
 assert [[p.get(i,F(0))for i in range(27)]for p in nu]==[list(map(F,p))for p in phase['nu_A_B27']]
 ga=vals[7790];assert [ga.get(i,F(0))for i in range(27)]==list(map(F,phase['gA_source7790_B27']))
 assert all(not p.get(0)for p in nu)
 return B,leaves,pairs,nu,ga,old26,peak
def actual_proper(D,B,leaves,pairs,nu,gauge,given):
 ring=R108(B,nu);ring.install();I,P=D['interface'],D['proper'];GL=D['first']['GL18_source_parameter_positions'];S=D['source']
 def frame(family):
  f=S['families'][family];p=ring.add(ring.tau(0 if family=='CA'else 1),ring.scale(ring.base(leaves[f['p_quadratic_relation']['s']]),F(1,2)));M=V.eye(5)
  for name,(x,z)in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);M[i][j]=ring.add(M[i][j],ring.base(leaves[x]),ring.mul(ring.base(leaves[z]),p))
  C=V.eye(5)
  if family=='CA':C[0][3]=ring.scale(p,-2);C[1][0]=p;C[1][3]=ring.neg(ring.mul(p,p))
  else:C[0][3]=ring.scale(p,F(1,2));C[1][0]=p;C[1][3]=ring.scale(ring.mul(p,p),F(1,4))
  return V.matmul(M,C)
 HA,HB=frame('CA'),frame('CB');PA=[ring.add(ring.base(x),ring.mul(ring.base(z),ring.tau(0)))for x,z in pairs['CA']];ev=[Q.sparse(v)for v in gauge['actual_even15_B27']];od=[Q.sparse(v)for v in gauge['actual_odd15_B27']];G=V.eye(5)
 for(i,j),x,z in zip(I['actual_group15_positions'],ev,od):G[i][j]=ring.add(G[i][j],ring.base(x),ring.mul(ring.base(z),ring.tau(0)))
 Fm=V.normalized_flag_frame(V.matmul(G,HA));assert Fm==V.eye(5)
 moved=V.tensor_push(PA,V.matmul(Fm,G));rows=P['zero_coordinate_rows36'];rr=P['CA_constant_left_inverse_rows19'];centers=list(map(F,P['actual_CA_source19_center']));co=[v[0]for v in V.matmul(V.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in rr])];assert all(co[j]==ring.c(centers[j])for j in P['CA_fiber10_fixed_coefficient_indices'])
 FB=V.normalized_flag_frame(V.matmul(G,HB));actual=[ring.add(co[j],ring.c(-centers[j]))for j in P['CA_free_fiber9_coefficient_indices']]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j]for j in I['retained_source9_indices']]
 dense=lambda p:[p.get((i,m),F(0))for i in range(27)for m in range(4)]
 assert list(map(dense,actual))==[list(map(F,p))for p in given['actual14_C108']]and list(map(dense,kappa))==[list(map(F,p))for p in given['kappa9_C108']]
 return ring,kappa,ev,od
def poly(r,p,vs):
 z={}
 for term in p:
  e=term[0];c=F(term[1],term[2])if len(term)==3 else F(term[1]);v=r.c(c)
  for x,n in zip(vs,e):
   if n:v=r.mul(v,r.power(x,n))
  z=r.add(z,v)
 return z
def Vcolumns(D,ring,kappa,ga):
 a,d,x,t=[kappa[i]for i in [0,1,2,6]];du=ring.add(ring.c(1),ring.neg(ring.mul(d,d)));pu=ring.add(ring.c(5),ring.mul(d,d));u=[];K=D['core5']
 assert du.get((0,0),F(0))==1 and pu.get((0,0),F(0))==5
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];v=ring.tm(poly(ring,p['numerator'],[a,d,x,t]),ring.power(du,-p['D_power']),ring.power(pu,-p['P_power']));u.append(ring.add(kappa[K['source9_variable_order'].index(name)],ring.neg(v)))
 h=ring.add(a,ring.neg(d),x,ring.scale(t,-2));c=ring.scale(ring.add(d,x),F(1,2));r=ring.scale(ring.add(ring.scale(a,-2),d,ring.neg(x)),F(1,6));l=ring.scale(ring.add(a,d,ring.neg(x)),F(1,3));vs=u+[h,c,r,l];R=D['rad'];dd=poly(ring,R['d_polynomial'],vs);assert dd==d;cols=[]
 for p in R['R9_generators_polynomials9']:
  v=ring.mul(dd,poly(ring,p,vs))
  for mask in range(4):cols.append(ring.mul(v,{(0,mask):F(1)}))
 for i in range(8):cols.append({(ring.b.idx[(i,)],1):F(-1)}if i<7 else{(j,1):-c for j,c in ga.items()})
 assert len(cols)==44;return cols
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(stem,name='report.json',keys=None):
  p=ROOT/'reports'/stem/name;pins[str(p.relative_to(ROOT))]=sha(p)
  if keys is None:return json.loads(p.read_text())
  text=p.read_text();decoder=json.JSONDecoder();d={}
  for key in keys:
   token=json.dumps(key);assert text.count(token)==1;pos=text.index(':',text.index(token))+1
   while text[pos].isspace():pos+=1
   d[key]=decoder.raw_decode(text,pos)[0]
  return d
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63';rep=load(primary);assert rep['all_passed']
 assert rep['code_sha256']=='16e83e2bc7d50236eb65440465d76eff09388d17d5fb02bf7e1783dcc6541a14'
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 assert load('AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE_20261105_v63_REVIEW')['source_frame_budget_and_prepared_V2_gate_passes']
 assert load('AMBIENT_B18_CPRIME_QUADRATIC20_IMPORT_CORRECTION_20261105_v63_REVIEW')['authorized_one_V3_value_gate_ready']
 p1=load(primary,'phase1-common105-and-114-source-equations.json');gauge=load(primary,'phase2-thirty-ell-square-gauge-values.json');proper=load(primary,'phase3-actual-proper-B27-values.json');values=load(primary,'phase4-all-twenty-actual-quadratic-values.json')
 names={
  'gate':('AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE_PREFLIGHT_20261105_v63','complete-four-by20-and-ell-square-cache-source-gate.json',None),
  'subset':('AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE2_SUBSET_20261105_v63','minimal-two57-source2-and-labels.json',None),
  'source':('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json',['DAG_scalar_nodes','families']),
  'adapter':('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63','actual-special-fiber20-coefficients.json',['YB_in_YA_linear','YB_in_YA_quadratic']),
  'oldsource':('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-source105-ell2-raw33-values.json',['joint105_unknown_order','all105_unknown_final_values','all_requested164_root_values']),
  'interface':('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63','actual-source-physical-firstjet-interface.json',['actual_group15_positions','retained_source9_indices']),
  'proper':('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63','literal-proper-source14-preflight.json',['zero_coordinate_rows36','CA_constant_left_inverse_rows19','CA_fiber10_fixed_coefficient_indices','CA_free_fiber9_coefficient_indices','actual_CA_source19_center','CA_left_inverse19']),
  'first':('AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63','actual-native128-first-input-and-mixed-index-profile.json',['GL18_source_parameter_positions']),
  'core5':('AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63','actual4variable7core-polynomials.json',['source9_variable_order','selected5_pivot_order','five_actual_rational_roots']),
  'rad':('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json',['R9_generators_polynomials9','d_polynomial']),
  'oldgauge':('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63','phase2-canonical30-exact-values.json',None),
  'oldproper':('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63','phase3-actual-kappa9-and-normalized-nilframe.json',['actual14_C132','kappa9_C132']),
  'Vfirst':('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase1-actual-fullnine-first-V-columns.json',None),
  'J30':('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase2-fixed30-Jacobian.json',None)}
 D={k:load(*args)for k,args in names.items()};save('phase.json',dict(phase='immutable primary/source and PROVIDED30 gauge/4x20 certificates pinned; no solve'))
 B,leaves,pairs,nu,ga,old26,streampeak=source_inputs(D,p1)
 ring,kappa,ev,od=actual_proper(D,B,leaves,pairs,nu,gauge,proper)
 assert all({j:c for j,c in v.items()if j<26}==old26(Q.sparse(old))for current,oldrows in [(ev,D['oldgauge']['actual_even15_B33']),(od,D['oldgauge']['actual_odd15_B33'])]for v,old in zip(current,oldrows))
 cor=[p.get(26,F(0))for p in ev+od];assert cor==list(map(F,gauge['correction30']))
 J=D['J30'];a0=[[F(c)for c in row]for row in J['A0']];a1=[[F(c)for c in row]for row in J['A1_eta_delta']];jj=[a0[i]+[F(0)]*15 for i in range(15)]+[a1[i]+a0[i]for i in range(15)];assert [-sum((c*v for c,v in zip(row,cor)),F(0))for row in jj]==list(map(F,gauge['RHS30']))
 for key,new in [('actual14_C132',[Q.sparse(p)for p in proper['actual14_C108']]),('kappa9_C132',kappa)]:
  # Convert dense actual14 using its explicit four-mask indexing.
  if key=='actual14_C132':new=[{(i//4,i%4):c for i,c in p.items()}for p in new]
  for p,old in zip(new,D['oldproper'][key]):
   for mask in range(4):assert {i:c for(i,m),c in p.items()if m==mask and i<26}==old26(Q.sparse(old[mask::4]))
 cols=Vcolumns(D,ring,kappa,ga);first=[Q.pairs(c)for c in D['Vfirst']['all44_actual_linear_columns27_sparse']];mons2=list(itertools.combinations_with_replacement(range(9),2));plane=[(7,7),(7,8),(8,8)];pi=[Q.pairs(v)for v in D['gate']['actual_NF4_rows_by135_sparse']];raw=[];nf=[]
 def gamma(v,i):return [v.get((i,1),F(0)),v.get((i,3),F(0)),v.get((i,2),F(0))-F(11,12)*v.get((i,1),F(0))]
 for j,v in enumerate(cols):
  assert gamma(v,0)==[F(0)]*3
  assert {9*t+i:c for i in range(9)for t,c in enumerate(gamma(v,B.idx[(i,)]))if c}==first[j]
  row=[gamma(v,B.idx[m])[t]for t in range(3)for m in plane];raw.append(row);pol={45*t+mons2.index(m):gamma(v,B.idx[m])[t]for t in range(3)for m in plane};nf.append([sum((c*pol.get(i,F(0))for i,c in p.items()),F(0))for p in pi])
 assert raw==[list(map(F,v))for v in values['all44_actual_V2_pureplane_nine_values']]and nf==[list(map(F,v))for v in values['all44_projected_NF4_values']]
 kernel=[Q.pairs(v)for v in D['gate']['all20_constant_kernel_columns44_sparse']];M=[[sum((c*nf[j][i]for j,c in k.items()),F(0))for i in range(4)]for k in kernel];assert M==[list(map(F,v))for v in values['all20_actual_quadratic_values4']]
 cp,rr=values['fixed_rank_columns'],values['fixed_rank_rows'];rank=len(cp);iv=[[F(c)for c in row]for row in values['fixed_minor_inverse']];assert len(rr)==len(iv)==rank and len(set(cp))==len(set(rr))==rank
 A=[[M[j][i]for j in cp]for i in rr]
 assert all(sum((A[i][k]*iv[k][j]for k in range(rank)),F(0))==F(i==j)and sum((iv[i][k]*A[k][j]for k in range(rank)),F(0))==F(i==j)for i in range(rank)for j in range(rank))
 det=F(0)
 for perm in itertools.permutations(range(rank)):
  p=F((-1)**sum(perm[i]>perm[j]for i in range(rank)for j in range(i+1,rank)))
  for i,j in enumerate(perm):p*=A[i][j]
  det+=p
 assert det==F(values['fixed_minor_determinant'])and det
 digest=hashlib.sha256()
 for col in M:
  coef=[sum((c*col[rr[j]]for j,c in enumerate(row)),F(0))for row in iv];assert [sum((c*M[cp[j]][i]for j,c in enumerate(coef)),F(0))for i in range(4)]==col;digest.update(json.dumps(enc(coef),separators=(',',':')).encode()+b'\n')
 assert digest.hexdigest()==values['all20_factorization_sha256']and rank==values['actual_complete_rank']and 4-rank==values['actual_Cprime_HF2']
 save('independent-actual114-source-gauge-and-complete4by20-certificates.json',dict(all105_full9_cache_compositions=True,all9_shared_physical_coordinates=True,all114_original_equations=True,all164_B26_root_overlaps=True,all600_physical_readback_hashes=True,provided30_gauge_coefficients_verified_without_solving=True,all30_original_constraints=True,all_B26_gauge_actual14_kappa9_overlaps=True,all44_fullnine_V1_and_origins_match=True,fullproper_rational_functions_beforetrace=True,all396_pureplane_V2_coefficients_and176_projected_values=True,all80_complete_quadratic_entries=M,provided_fixed_rank=rank,provided_fixed_minor_det=det,both_inverse_products=True,all20_factorizations_sha256=digest.hexdigest(),actual_Cprime_HF_through2=[3,3,4-rank],no_long_exponent_or_nilcube_conclusion=True,no_source_gauge_jet_or_coefficient_solve=True,sparse_source_peak=streampeak,trace_inverse_products=ring.inverse_count))
 for p in [Path(__file__),Path(Q.__file__),ROOT/'code/ambient_B18_honest_nil_T_gauge120_review_20261105_v63.py',ROOT/'code/ambient_B18_Cprime_ell2_nilimage_values_review_20261105_v63.py',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 assert time.monotonic()-st<60 and sum(p.stat().st_size for p in OUT.iterdir()if p.is_file())<32<<20
 save('report.json',dict(id=ID,all_passed=True,frozen=True,actual_complete_quadratic_rank=rank,actual_Cprime_HF_through2=[3,3,4-rank],all_source_gauge_proper_and_complete_rank_gates_verified=True,no_new_source_jet_or_solve=True,cubic_or_long_exponent_not_decided=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
 print(json.dumps(dict(all_passed=True,actual_rank=rank,actual_HF=[3,3,4-rank],seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
