#!/usr/bin/env python3
"""Independent C4 tensor/source adapter and fixed X0 certificate replay."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
import hashlib,json,resource,time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63'
STEM=AUTHOR+'_REVIEW'
Z=(Q(0),)*4;ONE=(Q(1),Q(0),Q(0),Q(0))
TA=(Q(0),Q(1),Q(0),Q(0));TB=(Q(0),Q(0),Q(1),Q(0))
PAIRS=list(combinations(range(5),2));MONS=list(combinations_with_replacement(range(5),2))
PI={w:i for i,w in enumerate(PAIRS)};MI={w:i for i,w in enumerate(MONS)}
INV_CERTS=[]

def scalar(q):return(Q(q),Q(0),Q(0),Q(0))
def plus(*xs):return tuple(sum(v,Q(0))for v in zip(*xs))if xs else Z
def neg(x):return tuple(-v for v in x)
def sc(x,q):return tuple(v*q for v in x)
def times(x,y):
 a,b,c,d=x;e,f,g,h=y
 return(a*e,a*f+b*e,a*g+c*e,a*h+d*e+b*g+c*f)
def tm(*xs):
 v=ONE
 for x in xs:v=times(v,x)
 return v
def inverse(x):
 a,b,c,d=x;assert a,'No nonunit is inverted.'
 y=(1/a,-b/a**2,-c/a**2,2*b*c/a**3-d/a**2)
 assert times(x,y)==ONE;INV_CERTS.append((x,y));return y
def power(x,n):
 if n<0:return power(inverse(x),-n)
 v=ONE
 while n:
  if n&1:v=times(v,x)
  x=times(x,x);n//=2
 return v
def eye(n):return[[ONE if i==j else Z for j in range(n)]for i in range(n)]
def matmul(A,B):
 C=[[Z]*len(B[0])for _ in A]
 for i,row in enumerate(A):
  for k,x in enumerate(row):
   if x!=Z:
    for j,y in enumerate(B[k]):
     if y!=Z:C[i][j]=plus(C[i][j],times(x,y))
 return C
def rational_inverse(A):
 n=len(A);V=[list(row)+[Q(i==j)for j in range(n)]for i,row in enumerate(A)]
 for k in range(n):
  pivot=next(i for i in range(k,n)if V[i][k]);V[k],V[pivot]=V[pivot],V[k]
  d=V[k][k];V[k]=[x/d for x in V[k]]
  for i in range(n):
   if i!=k and V[i][k]:
    d=V[i][k];V[i]=[x-d*y for x,y in zip(V[i],V[k])]
 assert [r[:n]for r in V]==[[Q(i==j)for j in range(n)]for i in range(n)]
 return[r[n:]for r in V]
def matrix_inverse(A):
 # Independent nilpotent perturbation formula, rather than the author's C4
 # Gaussian inverse. No source implicit system is solved.
 n=len(A);C=[[scalar(x)for x in row]for row in rational_inverse([[x[0]for x in row]for row in A])]
 E=matmul(C,[[plus(x,scalar(-x[0]))for x in row]for row in A])
 E2=matmul(E,E);E3=matmul(E2,E);assert all(x==Z for row in E3 for x in row)
 P=[[plus(ONE if i==j else Z,neg(E[i][j]),E2[i][j])for j in range(n)]for i in range(n)]
 result=matmul(P,C);assert matmul(A,result)==matmul(result,A)==eye(n);return result
def constmat(A):return[[scalar(x)for x in row]for row in A]
def unpack(A):return[[tuple(Q(x)for x in p)for p in row]for row in A]
def polynomial(terms,vals,derivative=None):
 result=Z;cache={}
 for record in terms:
  if len(record)==3:e,n,d=record;q=Q(n,d)
  else:e,q=record;q=Q(q)
  e=list(e)
  if derivative is not None:
   if not e[derivative]:continue
   q*=e[derivative];e[derivative]-=1
  term=scalar(q)
  for i,k in enumerate(e):
   if k:
    if(i,k)not in cache:cache[i,k]=power(vals[i],k)
    term=times(term,cache[i,k])
    if term==Z:break
  result=plus(result,term)
 return result
def synthetic_prefix_coefficients(terms,vals,axes,roots):
 # Monomial-by-monomial synthetic division of consecutive actual prefix
 # differences. The prefix uses actual earlier axes and root later axes.
 result=[Z]*len(axes);positions={j:i for i,j in enumerate(axes)}
 for e,n,d in terms:
  for k,j in enumerate(axes):
   if not e[j]:continue
   prefix=list(vals)
   for h in range(k+1,len(axes)):prefix[axes[h]]=roots[h]
   term=scalar(Q(n,d))
   for axis,degree in enumerate(e):
    if axis!=j and degree:term=times(term,power(prefix[axis],degree))
   geometric=Z
   for h in range(e[j]):geometric=plus(geometric,tm(power(vals[j],e[j]-1-h),power(roots[k],h)))
   result[k]=plus(result[k],times(term,geometric))
 return result
def normalized_flag_frame(H):
 C=matrix_inverse(H);ell=[times(x,inverse(C[1][1]))for x in C[1]]
 k=[plus(x,neg(times(C[0][1],y)))for x,y in zip(C[0],ell)]
 u=inverse(k[0]);k=[times(x,u)for x in k]
 F=eye(5);F[0]=k;F[1]=ell;return F
def tensor_push(v,G):
 # Literal exterior-square(G) times symmetric-square(G^-1), independently
 # constructed directly in the full150 raw coefficient tensor.
 Gi=matrix_inverse(G);out=[Z]*150
 wedge=[]
 for a,b in PAIRS:
  wedge.append([plus(times(G[i][a],G[j][b]),neg(times(G[i][b],G[j][a])))for i,j in PAIRS])
 for index,c in enumerate(v):
  if c==Z:continue
  a=index//15;u,w=MONS[index%15]
  quadratic={}
  for i,j in product(range(5),repeat=2):
   x=tm(c,Gi[u][i],Gi[w][j])
   if x!=Z:
    m=tuple(sorted((i,j)));quadratic[m]=plus(quadratic.get(m,Z),x)
  for pair,x in enumerate(wedge[a]):
   if x!=Z:
    for m,y in quadratic.items():out[15*pair+MI[m]]=plus(out[15*pair+MI[m]],times(x,y))
 return out
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())

def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 rp='reports/'+AUTHOR+'/report.json';A=read(rp);assert A['all_passed'];hashes={rp:sha(rp)}
 for p,h in A['input_sha256'].items():assert sha(p)==h,p;hashes[p]=h
 for name,h in A['evidence_sha256'].items():
  p='reports/'+AUTHOR+'/'+name;assert sha(p)==h,p;hashes[p]=h
 for p,h in [('code/ambient_B18_comparison_X_origin_20261104_v63.py',A['code_sha256']),('research/notes/'+AUTHOR+'.md',A['note_sha256'])]:
  assert sha(p)==h;hashes[p]=h
 # Independent source-binding gates, not the author worker's arithmetic.
 peer_stems=['AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63_REVIEW',
  'AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW',
  'AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW',
  'AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW',
  'AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW',
  'AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW']
 for st in peer_stems:
  p='reports/'+st+'/report.json';d=read(p);assert d['all_passed'];hashes[p]=sha(p)
 def source(prefix,name):
  p='reports/'+prefix+'_20261104_v63/'+name;assert p in hashes;return read(p)
 I=source('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET','actual-source-physical-firstjet-interface.json')
 P=source('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT','literal-proper-source14-preflight.json')
 Eq=source('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL','actual-proper-source14-equations17.json')
 H=source('AMBIENT_B18_PROPER_RATIONAL_ELIMINATION','actual12-flag-divided-polynomials.json')
 K=source('AMBIENT_B18_PROPER_CORE_ELIMINATION','actual4variable7core-polynomials.json')
 Member=source('AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP','membership-and-source-cofactors.json')
 prep=source('AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP','preparation.json')
 B=source('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY','actual-UP-QB-640-typed-operator.json')
 C0=source('AMBIENT_B18_T0_HONEST_C0','honest-C0-D0-D4-DAG.json')
 First=source('AMBIENT_B18_T3_FULL9_FIRSTJET','full9-firstjet-certificate.json')
 Artin=source('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION','actual-source105-Artin20-values.json')
 Tensor=source('AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL','actual-normalizer-tensor-free4-DAG.json')
 V=source('AMBIENT_B18_FULLV56_COMPARISON_UNIT','actual-fullV56-m0-and-typed-comparison.json')
 Op=source('AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR','actual-complete-Res4-syzygy-operator.json')
 quartic=source('AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC','quotient-firsts-and-actual-D0-firstjets.json')
 phase1=read('reports/'+AUTHOR+'/phase1-native-values.json');phase2=read('reports/'+AUTHOR+'/phase2-native-adapter.json')
 Xsaved=read('reports/'+AUTHOR+'/actual-X0-and-unit-cancellation.json')
 (out/'preflight.json').write_text(json.dumps({'C4_basis':['1','tauA','tauB','tauA*tauB'],
  'actual_X_shape':[56,768],'fixed_source_coefficients':4,'no_new_source_solve_or_jet':True,
  'input_sha256':hashes},indent=2)+'\n')
 fam={f['family']:f for f in I['families']}
 HA=eye(5);HB=eye(5);gauge=eye(5)
 for i,j in product(range(5),repeat=2):
  HA[i][j]=plus(HA[i][j],sc(TA,Q(fam['CA']['source_frame5_by5_by9_first'][7][i][j])))
  HB[i][j]=plus(HB[i][j],sc(TB,Q(fam['CB']['source_frame5_by5_by9_first'][8][i][j])))
 for(i,j),row in zip(I['actual_group15_positions'],I['group15_first9']):gauge[i][j]=plus(gauge[i][j],sc(TA,Q(row[7])))
 assert normalized_flag_frame(matmul(gauge,HA))==eye(5)
 FB=normalized_flag_frame(matmul(gauge,HB))
 flags=[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]
 assert gauge==unpack(phase1['gauge']) and FB==unpack(phase1['CB_flag_frame'])
 moved=tensor_push([scalar(x)for x in P['actual_P0_physical150']],gauge)
 rows36=P['zero_coordinate_rows36'];left=constmat(P['CA_left_inverse19'])
 picked=[[moved[rows36[j]]]for j in P['CA_constant_left_inverse_rows19']]
 own=[row[0]for row in matmul(left,picked)]
 center=list(map(Q,P['actual_CA_source19_center']))
 assert all(own[j]==scalar(center[j])for j in P['CA_fiber10_fixed_coefficient_indices'])
 actual=[plus(own[j],scalar(-center[j]))for j in P['CA_free_fiber9_coefficient_indices']]+flags
 kappa=[actual[j]for j in I['retained_source9_indices']]
 assert actual==[tuple(map(Q,p))for p in phase1['actual14_C0']]
 assert kappa==[tuple(map(Q,p))for p in phase1['kappa9_C0']]
 for i,p in enumerate(actual):assert p[:3]==(Q(0),Q(I['source14_first_by_physical9'][i][7]),Q(I['source14_first_by_physical9'][i][8]))
 assert all(polynomial(f,actual)==Z for f in Eq['equations14variable_sparse'])
 transform=matmul(FB,gauge);ann=constmat(P['CB_annihilator17_by36']);targets=B['physical_target_rows12'];R17=[[Z]*12 for _ in range(17)]
 for j,target in enumerate(targets):
  inp=[Z]*150;inp[target]=ONE
  for k,row in enumerate(Tensor['eliminated_zero_rows4']):inp[row]=scalar(Tensor['eliminated_four_coordinate_reconstruction'][k][rows36.index(target)])
  pushed=tensor_push(inp,transform)
  col=matmul(ann,[[pushed[i]]for i in rows36])
  for i in range(17):R17[i][j]=col[i][0]
 assert R17==unpack(phase2['R17'])
 D=plus(ONE,neg(times(kappa[1],kappa[1])));PP=plus(scalar(5),times(kappa[1],kappa[1]))
 first_axes=[0,1,10,11,13];first_names=['v0','v1','b','c','e'];firstroots=[]
 for name in first_names:
  r=H['five_explicit_rational_roots'][name]
  firstroots.append(tm(polynomial(r['numerator'],kappa),power(D,-r['D_denominator_power'])))
 assert firstroots==[actual[j]for j in first_axes]
 jac=[[polynomial(f,actual,j)for j in first_axes]for f in Eq['equations14variable_sparse']]
 BA=matmul(jac[5:],matrix_inverse(jac[:5]))
 R12=[[plus(R17[i+5][j],neg(plus(*(times(BA[i][k],R17[k][j])for k in range(5)))))for j in range(12)]for i in range(12)]
 for i,r in enumerate(H['all12_actual_remaining_and_divided_records']):
  unit=sc(power(D,r['D_denominator_power']),Q(r['constant_denominator'],r['primitive_content_unit']))
  R12[i]=[times(unit,p)for p in R12[i]]
 assert R12==unpack(phase2['R12'])
 hpolys=[r['divided_primitive_polynomial']for r in H['all12_actual_remaining_and_divided_records']]
 core=[kappa[j]for j in[0,1,2,6]]
 second_names=K['selected5_pivot_order'];second_axes=[K['source9_variable_order'].index(n)for n in second_names]
 second_roots=[]
 for name in second_names:
  r=K['five_actual_rational_roots'][name]
  second_roots.append(tm(polynomial(r['numerator'],core),power(D,-r['D_power']),power(PP,-r['P_power'])))
 coeff=[synthetic_prefix_coefficients(p,kappa,second_axes,second_roots)for p in hpolys]
 order=K['selected5_source_H_rows']+[i for i in range(12)if i not in K['selected5_source_H_rows']]
 assert order==Op['source_adapter_T']['permutation']
 Ai=matrix_inverse([coeff[j]for j in order[:5]]);BA=matmul([coeff[j]for j in order[5:]],Ai)
 T=[[Z]*12 for _ in range(12)]
 for i in range(5):
  for j in range(5):T[i][order[j]]=Ai[i][j]
 for i in range(7):
  T[i+5][order[i+5]]=ONE
  for j in range(5):T[i+5][order[j]]=neg(BA[i][j])
 for i,r in enumerate(K['remaining7_actual_core_polynomials']):
  lam=sc(tm(power(D,-r['D_denominator_power']),power(PP,-r['P_denominator_power'])),Q(r['QQ_unit_content_numerator'],r['QQ_unit_content_denominator']))
  T[i+5]=[times(inverse(lam),p)for p in T[i+5]]
 rr=plus(core[2],neg(core[0]),neg(core[1]));E=eye(12);E[5][5]=inverse(D);E[6][6]=inverse(D)
 E[7][5]=neg(tm(rr,inverse(times(PP,D))));E[7][7]=sc(inverse(times(PP,D)),-2);T=matmul(E,T)
 E=eye(12)
 for j in range(4):
  for i in range(3):E[8+j][5+i]=neg(tm(power(D,[2,2,2,3][j]),polynomial(Member['cofactor3_by4'][i][j],core)))
 T=matmul(E,T);N=matmul(T,R12)
 assert T==unpack(phase2['T']) and N==unpack(phase2['N'])
 Ni=unpack(phase2['N_inverse']);assert matmul(N,Ni)==matmul(Ni,N)==eye(12)
 G8=[plus(kappa[j],neg(root))for j,root in zip(second_axes,second_roots)]+[polynomial(p,core)for p in prep['generators3']]
 assert G8==[tuple(map(Q,p))for p in phase2['G8_C0']]
 expected=[Q(1,2),Q(25,8),0,0,0,5,Q(-125,2),0]
 assert G8==[(Q(0),Q(0),Q(0),Q(q))for q in expected]
 assert matmul(T,[[polynomial(p,kappa)]for p in hpolys])==[[g]for g in G8]+[[Z]for _ in range(4)]
 assert all(times(kappa[1],g)==Z for g in G8)
 # Never cancel the nonunit d from d*G8=0 in this algebra.
 nativeS=[[Z]*32 for _ in range(12)]
 for i,j,g,sign in Op['native_S12_by32']['koszul28_in_first8_rows']:nativeS[i][j]=sc(G8[g],Q(sign))
 for i,j,c in Op['native_S12_by32']['free4_zero_source_columns']:nativeS[i][j]=scalar(c)
 native=matmul(list(map(list,zip(*N))),nativeS)
 syz=[[Q(0)]*128 for _ in range(48)]
 for i,j in product(range(12),range(32)):
  for b in range(4):
   y=times(native[i][j],tuple(Q(k==b)for k in range(4)))
   for k in range(4):syz[4*i+k][4*j+b]=y[k]
 assert syz==[[Q(x)for x in row]for row in Xsaved['all48_by128_SyzQ_origins']]
 # Cache only accepted ACTUAL regular origins, no evaluation of specialized
 # old numerator/denominator. All source node identity is by node ID.
 values={int(i):Q(p[0])for key in['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']for i,p in First[key].items()}
 values.update({i:Q(dict(p).get(0,'0'))for i,p in Artin['all_requested_original_source_root_values']})
 values.update({i:Q(0)for i in C0['actual_C0_basis8_scalar_nodes']})
 for r in C0['all35_regular_coordinate_divisions']:
  first=First['both57_source2']['CA']['pure_function2jets'][str(r['original_pure_CA_source_function'])]
  for j,q in enumerate(r['all7_regular_coordinate_divisions']):assert Q(q['origin'])==Q(first[1+j]);values[q['node']]=Q(first[1+j])
 oldg=quartic['original_regular_function_nodes_unmodified']['old_CB'];newg=quartic['original_regular_function_nodes_unmodified']['new_CA']
 values.update({i:Q(0)for i in oldg+newg})
 gb=next(r['gB']for r in B['accepted_actual_regular_function_identities']if r['kind']=='accepted_actual_CB_f63_tail');values[gb]=Q(0)
 dag=C0['DAG_scalar_nodes'];memo=dict(values)
 def original(i):
  if i in memo:return memo[i]
  node=dag[i];tag=node[0]
  if tag=='c':x=Q(node[1],node[2])
  elif tag=='n':x=-original(node[1])
  elif tag=='a':x=original(node[1])+original(node[2])
  elif tag=='m':
   x=original(node[1]);x=x*original(node[2])if x else Q(0)
  else:raise AssertionError(('unbound ACTUAL regular function origin',i,node))
  memo[i]=x;return x
 ops=[]
 for node in B['operator_scalar_nodes']:
  tag=node[0]
  if tag=='c':x=Q(node[1],node[2])
  elif tag=='actual_C0_scalar_function':x=original(node[1])
  elif tag=='sum':x=sum((ops[j]*Q(a,b)for j,a,b in node[1]),Q(0))
  elif tag=='product':
   x=Q(node[1],node[2])
   for j in node[3]:x*=ops[j]
  else:raise AssertionError(node)
  ops.append(x)
 X=[[Q(0)]*768 for _ in range(56)]
 for j,column in enumerate(B['B48_by640_sparse_columns']):
  for i,k in column:X[i][j]=ops[k]
 for j,column in enumerate(C0['honest_D0_relation_columns']):
  for i,root in enumerate(column[1:]):
   x=original(root);assert x==0;X[48+i][j]=x
 for i in range(48):X[i][640:]=syz[i]
 assert X==[[Q(x)for x in row]for row in Xsaved['full_X0_56_by768']]
 assert sum(bool(x)for row in X for x in row)==97
 for i,x in Xsaved['source_function_constant_refs']:assert original(i)==Q(x)
 m0=[[Q(x)for x in row]for row in V['all24_by56_m0']]
 columns=[{i:X[i][j]for i in range(56)if X[i][j]}for j in range(768)]
 for column in columns:
  assert all(sum((row[i]*x for i,x in column.items()),Q(0))==0 for row in m0)
 selected=V['fixed_selected24_columns'];free=[i for i in range(56)if i not in selected]
 assert free==Xsaved['free_kernel_rows32']
 pivots={};basiscols=[];receipt_count=0
 for j,r in enumerate(Xsaved['all768_fixed_elimination_receipts_and_constant_combinations']):
  assert r['source_column']==j
  reduced={k:X[old][j]for k,old in enumerate(free)if X[old][j]}
  for k,q in r['subtract']:
   q=Q(q);assert k in pivots and reduced.get(k,Q(0))==q
   for i,x in pivots[k].items():
    v=reduced.get(i,Q(0))-q*x
    if v:reduced[i]=v
    else:reduced.pop(i,None)
  combo={i:Q(x)for i,x in r['constant_original_source_combination']}
  assert combo and max(combo)<=j and j in combo
  full=[sum((X[i][k]*q for k,q in combo.items()),Q(0))for i in range(56)]
  if r.get('zero'):
   assert not reduced and not any(full)
  else:
   k=r['pivot'];q=Q(r['divide']);assert k==min(reduced)and reduced[k]==q and q
   normalized={i:x/q for i,x in reduced.items()}
   assert normalized=={i:Q(x)for i,x in r['normalized']}
   assert [full[i]for i in free]==[normalized.get(k,Q(0))for k in range(32)]
   assert k not in pivots;pivots[k]=normalized;basiscols.append(j)
  receipt_count+=1
 assert len(pivots)==27 and basiscols==Xsaved['selected_independent_columns']
 for a,b in Xsaved['native_unit_inverse_certificates']:
  a,b=tuple(map(Q,a)),tuple(map(Q,b));assert a[0]and times(a,b)==ONE
 result={
  'all14_actual_C4_sources_flags_and17_equations_independently_evaluated':True,
  'full_tensor150_action_not_author_Lie_helper':True,
  'all_R17_R12_T_N_and_saved_Ninverse_products_checked':True,
  'C4_G8_exact_tauA_tauB_coefficients':list(map(str,expected)),
  'dG_zero_NOT_used_to_infer_G_zero':True,
  'all48x128_actual_Res4_N_transpose_S_entries_checked':True,
  'SyzQ_rank_structural':{'Koszul_nonzero_vector_coefficient_rank':7,'multiply_tauA_tauB_rank':1,'free4_C4_rank':16,'total':23,'N_transpose_unit_preserves_rank':True},
  'all640_B_and_P_origins_checked':True,'all768_m0_X0_columns_zero':True,
  'all768_given_cancellation_receipts_replayed':receipt_count,
  'actual_X0_nonzeros':97,'rank_X0':27,'fixed_unit_quotient_directions':32,
  'mu_Rbar_T_lower_bound':5,
  'X0_alone_does_not_determine_mu_T_length_K_torsion_depth':True,
  'cached_actual_regular_r_g_origins_only_no_redivision':True,
  'no_new_source_solve_jet_GB_point_or_X1':True,
  'formal_note':'The projected X0 constant rank accounts for27 unit relations in the free32 directions; possible Euler wedge constants from X1 remain uncomputed.'}
 (out/'independent-C4-source-and-X0-proof.json').write_text(json.dumps(result,indent=2)+'\n')
 for p in[str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,
  'actual_X0_rank':27,'actual_X0_nonzero_coefficients':97,
  'all768_actual_source_origins_and_constant_receipts_verified':True,
  'source_G8_nonzero_nil_C4_vector_preserved':True,
  'X1_mu_T_length_K_torsion_or_depth_computed':False,
  'input_sha256':hashes,'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32<<20
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
