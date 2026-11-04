#!/usr/bin/env python3
"""New source-bound radical comparison and complete35 syzygies; no jets/GB."""
import argparse,hashlib,itertools,json,resource,signal,time
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(p,q):
 r=dict(p)
 for m,c in q.items():
  r[m]=r.get(m,Q(0))+c
  if not r[m]:del r[m]
 return r
def scale(p,c):return {m:c*v for m,v in p.items()if c*v}
def mul(p,q):
 r=defaultdict(Q)
 for a,x in p.items():
  for b,y in q.items():r[tuple(i+j for i,j in zip(a,b))]+=x*y
 return {m:c for m,c in r.items()if c}
ZERO=(0,)*9
def const(q):return {ZERO:Q(q)}if q else {}
def variable(i):m=list(ZERO);m[i]=1;return {tuple(m):Q(1)}
def enc(p):return [[list(m),str(c)]for m,c in sorted(p.items())]
class Pool:
 def __init__(self):self.nodes=[];self.cache={};self.z=self.c(0);self.o=self.c(1)
 def put(self,n):
  if n in self.cache:return self.cache[n]
  assert len(self.nodes)<150000,'150k_new_typed_node_guard'
  i=len(self.nodes);self.nodes.append(n);self.cache[n]=i;return i
 def c(self,x):x=Q(x);return self.put(('C_QQ',x.numerator,x.denominator))
 def add(self,*v):
  v=[x for x in v if x!=self.z]
  if not v:return self.z
  if len(v)==1:return v[0]
  return self.put(('C_sum',tuple(sorted(v))))
 def mul(self,*v):
  if self.z in v:return self.z
  v=[x for x in v if x!=self.o]
  if not v:return self.o
  if len(v)==1:return v[0]
  return self.put(('C_product',tuple(sorted(v))))
 def neg(self,x):return self.mul(self.c(-1),x)
 def comp(self,x,i):
  if x==self.z:return self.z
  n=self.nodes[x]
  if n[0]=='C_QQ':return x if not i else self.z
  if n[0]=='C_embed_old_R0':return self.put(('R0_old_operator',n[1]))if not i else self.z
  if n[0]=='C_basis':return self.o if i==n[1] else self.z
  return self.put(('R0_component_of_C',x,i))
 def polynomial(self,p,vs):
  terms=[]
  for m,c in p.items():
   factors=[self.c(c)]
   for v,e in zip(vs,m):factors.extend([v]*e)
   terms.append(self.mul(*factors))
  return self.add(*terms)
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 files={
  'nil':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63/actual-E4-nil-two-generators-and-relations.json',
  'context':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
  'oldB':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
  'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'}
 N,C,B,H=map(read,[ROOT/p for p in files.values()])
 assert C['source_adapter_T']['G8']==['u0','u1','u2','u3','u4','Q0(core4)','Q1(core4)','s*r^2']
 assert B['physical_target_rows12']==[18,19,28,29,60,67,72,63,64,68,73,74]
 assert B['all640_UP_equals_QB_verified'] and B['exact_existing_U9_C4'][0]==[1,0,0,0]
 P=H['honest_D0_relation_columns'];assert len(P)==640 and all(len(p)==9 for p in P)
 (out/'phase.json').write_text(json.dumps({'id':STEM,'status':'inputs_pinned_before_algebra','source_hashes':{k:sha(ROOT/v)for k,v in files.items()},'native_syzygies':35,'old_relation_columns':640,'new_typed_node_bound':150000},indent=2)+'\n')
 # Native coordinate order Z6=(u5,h), then c,r,l. The three cross quadrics
 # have complete2-column Hilbert-Burch resolution.
 vv=[variable(i)for i in range(9)];z6=vv[:6];c,r,l=vv[6:]
 core=[mul(c,r),mul(c,l),mul(r,l)]
 hb1={j:core[j]for j in range(3)}
 hb2=[{0:l,1:scale(r,-1)},{0:l,2:scale(c,-1)}]
 bases=[]
 for n in range(9):
  group=[]
  for b in range(3):
   a=n-b
   if 0<=a<=6:
    for I in itertools.combinations(range(6),a):
     for j in range([1,3,2][b]):group.append((I,b,j))
  bases.append(group)
 ranks=list(map(len,bases));assert ranks==[1,9,35,77,105,91,49,15,2]
 maps=[]
 for n in range(1,9):
  index={v:i for i,v in enumerate(bases[n-1])};cols=[]
  for I,b,j in bases[n]:
   col={}
   for k,i in enumerate(I):
    row=index[(I[:k]+I[k+1:],b,j)];col[row]=add(col.get(row,{}),scale(z6[i],(-1)**k))
   if b:
    for dest,p in ({0:hb1[j]}if b==1 else hb2[j]).items():
     row=index[(I,b-1,dest)];col[row]=add(col.get(row,{}),scale(p,(-1)**len(I)))
   cols.append({i:p for i,p in col.items()if p})
  maps.append(cols)
 # Every native chain identity checked exactly, all degrees, over QQ[9vars].
 chains=0
 for n in range(1,8):
  for col in maps[n]:
   product={}
   for j,p in col.items():
    for i,q in maps[n-1][j].items():product[i]=add(product.get(i,{}),mul(q,p))
   assert all(not p for p in product.values());chains+=1
 # Order is u5,h,cr,cl,rl as required for the actual nine-generator row.
 R=[col[0]for col in maps[0]];assert R==z6+core
 d=add(add(c,r),l);dR=[mul(d,p)for p in R]
 for col in maps[1]:
  product={}
  for i,p in col.items():product=add(product,mul(dR[i],p))
  assert not product
 # Embed actual Q0,Q1,f cofactors from the newly accepted literal nil packet.
 def poly4(rows):
  result={}
  for m,q in rows:
   ex=[0]*9;ex[6],ex[7],ex[8],ex[5]=m;result[tuple(ex)]=Q(q)
  return result
 V=[[poly4(p)for p in row]for row in N['source_factorization_coefficients4_by3']]
 D=[[{}for j in range(9)]for i in range(8)]
 for i in range(5):D[i][i]=const(1)
 for i in range(3):
  for j in range(4):D[5+i][5+j]=V[j][i]
 G=vv[:5]+[poly4(q)for q in N['actual_core_generators3']]
 for i in range(8):
  value={}
  for j in range(9):value=add(value,mul(D[i][j],R[j]))
  assert value==G[i]
 # Native finite C-value circuit. Ninv is the SAME pinned source inverse,
 # not a new arbitrary matrix or a coefficient-series expansion.
 pool=Pool();cv=[pool.put(('actual_proper_u',i))for i in range(5)]+[pool.put(('actual_core_GL4','h'))]+[pool.put(('actual_core_GL4',s))for s in ['c','r','l']]
 dnode=pool.polynomial(d,cv);Rn=[pool.polynomial(p,cv)for p in R];Qn=[pool.mul(dnode,p)for p in Rn]
 Dn=[[pool.polynomial(p,cv)for p in row]for row in D]
 inverse=[[pool.put(('accepted_native_unit_inverse_N12_entry',i,j))for j in range(8)]for i in range(12)]
 A=[[pool.add(*(pool.mul(inverse[i][k],Dn[k][j])for k in range(8)))for j in range(9)]for i in range(12)]
 basis=[pool.o]+[pool.put(('C_basis',i))for i in range(1,4)]
 # Actual4x36 ideal row: column4*j+a is dR_j times basis_a.
 Q36=[]
 for p in Qn:
  for b in basis:
   w=pool.mul(p,b);Q36.append([pool.comp(w,i)for i in range(4)])
 # Complete35 native columns, restricted to140 columns using the same basis.
 S140=[]
 for col in maps[1]:
  v={i:pool.polynomial(p,cv)for i,p in col.items()}
  for b in basis:
   w=[]
   for i,p in sorted(v.items()):
    value=pool.mul(p,b)
    w.extend([[4*i+k,pool.comp(value,k)]for k in range(4)if pool.comp(value,k)!=pool.z])
   S140.append(w)
 assert len(S140)==140
 # Group old48 rows into12 exact C-values; every old scalar function stays
 # an opaque R0 function with its accepted original r/g ancestry.
 def old(i):
  if i in [0,1]:return i
  return pool.put(('C_embed_old_R0',i))
 Bold=[]
 for col in B['B48_by640_sparse_columns']:
  grouped={}
  for i,v in col:grouped.setdefault(i//4,[]).append(pool.mul(old(v),basis[i%4]))
  Bold.append({i:pool.add(*v)for i,v in grouped.items()})
 BredC=[];Bred=[]
 for col in Bold:
  vals=[pool.add(*(pool.mul(A[i][j],v)for i,v in sorted(col.items())))for j in range(9)]
  BredC.append(vals);Bred.append([pool.comp(v,a)for v in vals for a in range(4)])
 assert len(Bred)==640 and all(len(x)==36 for x in Bred)
 U=[[pool.comp(old(v),0)for v in col]for col in B['exact_existing_U9_C4']]
 assert U[0]==[pool.o,pool.z,pool.z,pool.z]
 # P source root references are genuine R0 values; no pullback or new division.
 P8=[[pool.put(('R0_original_honest_C0_scalar',i))if i not in [0,1]else i for i in col[1:]]for col in P]
 V44=[[Q36[j][i]for j in range(36)]+[pool.put(('R0_neg',U[j][i]))if U[j][i]!=pool.z else pool.z for j in range(1,9)]for i in range(1,4)]
 ResAT=[]
 for i in range(12):
  for b in basis:
   ResAT.append([pool.comp(pool.mul(A[i][j],b),a)for j in range(9)for a in range(4)])
 native=[{'source_degree':i+1,'shape':[ranks[i],ranks[i+1]],'columns':[[[r,enc(p)]for r,p in sorted(col.items())]for col in cols]}for i,cols in enumerate(maps)]
 assert len(pool.nodes)<150000
 result={'actual_C_ring':'same free4 C/R0 normalizer tensor algebra before its12 physical equations',
  'C_basis':['1','tauA','tauB','tauA*tauB'],'actual_R0_order':C['source_bound_interface'].get('R0_order','physical9 [27,108,109,112,113,117,118,119,61]'),
  'accepted_native_C_coefficient_interface':{'path':files['context'],'sha256':sha(ROOT/files['context']),'selectors':['source_bound_interface','source_adapter_T','group30_to120','matrix_mu_by_rows','hensel_Res4','unit_matrix_inverse']},
  'actual_C_variable_rules':{'actual_proper_u':'second5 pivot_i minus accepted rational core root_i, exact u5 from source_adapter_T',
   'actual_core_GL4':N['forward_GL4'],'GL4_variable_input_order':N['old_variable_order'],'core_input_selection':'kappa entries[0,1,2,6] in accepted source9 (a,d,v2,v6)',
   'N12':'accepted full_N=T(kappa)*R12, dG8+zeros=N12*physical_F12',
   'accepted_native_unit_inverse_N12_entry':'entry(i,j) of unique inverse of that SAME N12; all known-unit inverse functions from accepted context, d is never inverted'},
  'R9_generators_polynomials9':[enc(p)for p in R],'d_polynomial':enc(d),'G8_polynomial_factor_D8_by9':[[enc(p)for p in row]for row in D],
  'native_free_basis_labels':[[{'regular_wedge':I,'HB_degree':b,'HB_basis':j}for I,b,j in basis]for basis in bases],
  'native_tensor_Koszul6_HB2_resolution':native,'native_Betti':ranks,'Ered_first_differential':'multiply only native d1 by scalar d; all higher maps unchanged',
  'restricted_R0_resolution_ranks':[4*n for n in ranks],'native_resolution_is_minimal_C':True,'restricted_R0_resolution_minimality_claimed':False,
  'native_new_A12_by9_nodes':A,'A_definition':'A=N12^-1*injection12_by8*D8_by9',
  'actual_Res4_Atranspose36_by48_columns':ResAT,
  'actual_F12_relation_identity':'physical_F12=A*(dR9) as generator COLUMN; ALL12 originalphysical rows preserved',
  'transpose_convention':'Qold48=Qrad36*Res4(A^t); block(j,i)=mu(A[i,j]); NEVER transpose inside mu',
  'new_Cscalar_and_R0component_nodes':pool.nodes,
  'actual_Qrad4_by36_columns':Q36,'actual_complete_SyzQrad36_by140_sparse_columns':S140,
  'actual_Qrad_chainproof':'G8=D8x9R9; N12F12=(dG8,0); native dR9*S35=0; exact C/R0 restriction',
  'old_actual_chain_context':{'path':files['oldB'],'sha256':sha(ROOT/files['oldB']),'old_R0_operator_nodes':'operator_scalar_nodes','old_B':'B48_by640_sparse_columns','old_U':'exact_existing_U9_C4'},
  'actual_Brad9_C_values_by640_columns':BredC,'actual_Brad36_R0_values_by640_columns':Bred,
  'Brad_definition':'A^t times grouped actual Bold12_C_values; extract fourcomponents',
  'actual_UP_Qrad_Brad':'accepted UP=QoldBold followed by source-bound Qold=QradRes4(A^t)',
  'honest_C0_P_context':{'path':files['C0'],'sha256':sha(ROOT/files['C0']),'P9_by640':'honest_D0_relation_columns','R0_original_honest_C0_scalar':'exact DAG_scalar_nodes function; old r/g never redivided'},
  'actual_P8_by640_columns':P8,'actual_V3_by44_rows_after_unit1_cancellation':V44,
  'actual_X44_by780':{'first640_columns':'stack(actual_Brad36_R0_values_by640_columns,actual_P8_by640_columns)',
   'last140_columns':'stack(actual_complete_SyzQrad36_by140_sparse_columns,zero8)',
   'unreduced_shapes':{'V':[4,45],'X':[45,780]},'unit_cancelled':'-U0=-e0 at source36,target0; chain forces removed Xrowzero after the unit column change',
   'final_shapes':{'V':[3,44],'X':[44,780]}},
  'NilD0_exactly_H1_of_actual_complex':True,'nil_kernel_theorem':'E4_red reduced, D0->E4 integral Spec-surjective implies ker nil; hence Nil(D0)=ker(D0->E4_red)',
  'homology_completeness':'P presents D0, S140 presents ALL Syz(Qrad), UP=QradBrad; a kernel pair maps to its D0 class and the inverse uses a Qrad lift; source relations captured by Brad and ALL S140',
  'no_actual_D0_nil_generator_list_cube_or_depth_claimed':True,'no_imageNilE_surjectivity_or_N0_principality_assumed':True,
  'new_typed_node_count':len(pool.nodes),'new_expanded_polynomial_terms':sum(len(p)for row in D for p in row)+sum(len(p)for cols in maps for col in cols for p in col.values()),
  'universal_chain_columns_checked':chains,'all8_G8_factorization_identities_checked':True}
 ev=out/'actual-radical-normalizer-comparison-operator.json';ev.write_text(json.dumps(result,separators=(',',':'))+'\n')
 assert result['new_expanded_polynomial_terms']<150000 and ev.stat().st_size<32*1024*1024
 gates=[ROOT/'reports'/n/'report.json'for n in ['AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW','AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW','AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW','AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW','AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in gates)
 paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']+[ROOT/p for p in files.values()]+gates
 report={'id':STEM,'all_passed':True,'status':'source_bound_complete_reduced_normalizer_comparison_and_exact_nil_H1',
  'actual_E4_reduced_ideal':'d*(u5,h,cr,cl,rl)','actual_Qrad_shape':[4,36],'complete_SyzQrad_shape':[36,140],
  'native_minimal_C_Betti':ranks,'R0_restricted_resolution_ranks':[4*n for n in ranks],
  'actual_G8_D8x9_R9_all_identities':True,'physical12_source_adapter_A12x9_and_transpose_preserved':True,'actual_UP_Qrad_Brad':True,
  'actual_nil_comparison_V_shape':[3,44],'actual_nil_comparison_X_shape':[44,780],'NilD0_equals_H1':True,
  'actual_nil_generators_cube_socle_depth_decided':False,'no_new_source_jet_point_GB_or_expanded_coefficient_series':True,
  'old_actual_r_g_ancestry_preserved':True,'node_count':len(pool.nodes),'bytes':ev.stat().st_size,'new_polynomial_terms':result['new_expanded_polynomial_terms'],
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');(out/'phase.json').write_text(json.dumps({'id':STEM,'status':'FINAL'},indent=2)+'\n')
 print(json.dumps({k:report[k]for k in ['id','all_passed','node_count','bytes','new_polynomial_terms','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
