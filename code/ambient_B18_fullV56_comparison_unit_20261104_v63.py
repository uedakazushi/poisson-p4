#!/usr/bin/env python3
"""One fixed Artin20 extraction of m(0), with regular typed comparison maps."""
from fractions import Fraction as Q
import json,resource,time
from pathlib import Path
from ambient_B18_raw_cokernel_Artin20_evaluation_20261104_v63 import Algebra20,add,scale,rec,read_rec,sha
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63'
PATHS={
 'H20':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'source105':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
 'chainlift':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json'}
def mm(A,B):return [[sum((a*b for a,b in zip(row,col)),Q(0)) for col in zip(*B)] for row in A]
def inv(A):
 n=len(A);M=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(A)];det=Q(1)
 for j in range(n):
  k=next(k for k in range(j,n) if M[k][j])
  if k!=j:M[k],M[j]=M[j],M[k];det=-det
  c=M[j][j];det*=c;M[j]=[x/c for x in M[j]]
  for k in range(n):
   if k!=j and M[k][j]:
    c=M[k][j];M[k]=[x-c*y for x,y in zip(M[k],M[j])]
 return [r[n:] for r in M],det
def p5mul(a,b):
 z={}
 for i,x in a.items():
  for j,y in b.items():
   if i+j<5:
    v=z.get(i+j,Q(0))+x*y
    if v:z[i+j]=v
    else:z.pop(i+j,None)
 return z
def main():
 t=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 d={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()};H,S,B,R=(d[k] for k in ['H20','source105','chainlift','resolution'])
 A=Algebra20(H);assert S['all105_equations_exactly_zero'] and B['all640_UP_equals_QB_verified']
 cs=list(map(Q,R['CB_shift_coefficients7']));maps=[]
 for shifted in [False,True]:
  cols=[]
  for m in H['quotient_basis_monomials']:
   degree=m[7]+2*sum(m[:7]);c=Q(1)
   for x,k in zip(cs,m[:7]):c*=(-x if shifted else Q(0))**k
   cols.append({degree:c} if c and degree<5 else {})
  def pi(v):
   z={}
   for i,c in v.items():z=add(z,cols[i],c)
   return z
  assert cols[0]=={0:Q(1)}
  for i in range(20):
   for j in range(20):assert pi(A.mul({i:Q(1)},{j:Q(1)}))==p5mul(cols[i],cols[j])
  maps.append(cols)
 def project(k,v):
  z={}
  for i,c in v.items():z=add(z,maps[k][i],c)
  return z
 nuA,nuB=[read_rec(x) for x in S['nu_CA_CB_values']]
 def cmul(v,w):
  z=[{} for _ in range(4)]
  for i,x in enumerate(v):
   for j,y in enumerate(w):
    q=A.mul(x,y)
    if i&j&1:q=A.mul(q,nuA)
    if i&j&2:q=A.mul(q,nuB)
    z[i^j]=add(z[i^j],q)
  return z
 V=[];labels=[]
 for j,r in enumerate(S['all12_physical_coefficients']):
  f=[read_rec(r['delta']),read_rec(r['bA']),scale(read_rec(r['bB']),-1),{}]
  for k in range(4):
   e=[{} for _ in range(4)];e[k]={0:Q(1)}
   V.append(cmul(f,e)[1:]);labels.append({'Q48_column':4*j+k,'physical_row':r['physical_row'],'tensor_basis':k,'exact_H20_value':True})
 for i in range(7):
  V.append([{1+i:Q(-1)},{},{}]);labels.append({'negative_U8_column':i,'actual_function':'-zA'+str(i),'exact_H20_value':True})
 # Only its cubic coefficient can affect m(0).  Do not pretend the unknown
 # quartic tail is zero in the actual function.
 r3=H['quotient_basis_monomials'].index([0,0,0,0,0,0,0,3])
 V.append([{r3:-Q(125,24)},{},{}]);labels.append({'negative_U8_column':7,'actual_function':'-gA_plane(r,r)',
  'exact_H20_value':False,'m0_relevant_representative':'-(125/24)r3','discarded_order4_tail_does_not_affect_any_m0_slot':True})
 columns=[];curve_values=[]
 for col in V:
  ans=[];cv=[]
  for s,(v,n) in enumerate(zip(col,[3,3,2])):
   p=project(int(s==1),v);assert all(p.get(k,Q(0))==0 for k in range(n))
   ans.extend(v.get(i,Q(0)) for i in range(1,8));ans.append(p.get(n,Q(0)));cv.append(rec(p))
  columns.append(ans);curve_values.append(cv)
 assert len(columns)==56 and all(len(c)==24 for c in columns)
 piv={};chosen=[];receipt=[]
 for j,v in enumerate(columns):
  v={i:x for i,x in enumerate(v) if x};steps=[]
  while v:
   k=min(v)
   if k not in piv:
    c=v[k];piv[k]={i:x/c for i,x in v.items()};chosen.append(j);receipt.append({'source_column':j,'pivot':k,'divide':str(c),'subtract':steps});break
   c=v[k];steps.append([k,str(c)]);v=add(v,piv[k],-c)
 assert len(piv)==24
 m0=[list(row) for row in zip(*columns)];M0=[[row[j] for j in chosen] for row in m0]
 iM0,det=inv(M0);factor=mm(iM0,m0);assert mm(M0,factor)==m0
 # Finite regular-function operators, defined before any singular quotient.
 ops=[];cache={}
 def put(n):
  key=json.dumps(n,separators=(',',':'))
  if key in cache:return cache[key]
  i=len(ops);ops.append(n);cache[key]=i;return i
 zero=put(['zero']);functions=[]
 for c in B['exact_existing_Q48_C4']:
  functions.append([put(['actual_chain_scalar_restrict_ell0',op]) if op else zero for op in c[1:]])
 for c in B['exact_existing_U9_C4'][1:]:
  functions.append([put(['negative_actual_chain_scalar_restrict_ell0',op]) if op else zero for op in c[1:]])
 assert len(functions)==56
 m=[[zero]*56 for _ in range(24)]
 for j,f in enumerate(functions):
  for s,root in enumerate(f):
   if root==zero:continue
   prefixes=[put(['adapted_sector_prefix',root,s,k]) for k in range(8)]
   for i in range(7):m[8*s+i][j]=put(['regular_coordinate_Hadamard',prefixes[i+1],prefixes[i],s,i])
   m[8*s+7][j]=put(['exact_monomial_r_division',prefixes[0],[3,3,2][s],
                      'actual_sector_annihilator_membership_not_numerical_zeroing'])
 selected=[[row[j] for j in chosen] for row in m]
 unit=put(['formal_unit_matrix_inverse',selected,'actual_center_matrix_inverse',[[str(x) for x in row] for row in iM0]])
 n=[[zero]*24 for _ in range(56)]
 for i,j in enumerate(chosen):
  for k in range(24):n[j][k]=put(['unit_matrix_inverse_entry',unit,i,k])
 e={'actual_variable_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
  'V56_target_order':['tauA','tauB','tauA*tauB'],'V56_column_order':labels,
  'H20_to_r5_projections20_column_sparse':[[rec(x) for x in M] for M in maps],
  'all800_projection_product_identities':True,
  'V56_H20_values_or_explicit_m0_representative':[[rec(v) for v in c] for c in V],
  'all56_sector_curve_restrictions':curve_values,
  'm0_row_order':'sector A,B,H each zA0..zA6 then r^n, n=3,3,2',
  'all24_by56_m0':[[str(x) for x in row] for row in m0],
  'fixed_selected24_columns':chosen,'greedy_rank24_receipts':receipt,
  'fixed24_minor_determinant':str(det),'fixed24_minor_matrix':[[str(x) for x in row] for row in M0],
  'fixed24_minor_inverse':[[str(x) for x in row] for row in iM0],
  'all56_column_factorization':[[str(x) for x in row] for row in factor],
  'M0_times_factor_equals_m0':True,
  'typed_regular_function_operators':ops,'actual_V56_function_operator_refs':functions,
  'typed_m24_by56':m,'typed_n56_by24':n,'typed_unit24_inverse_node':unit,
  'h_is_zero':'m*n=I24 by the selected unit block; hence V*n=d1 exactly',
  'typed_map_proof':'ordered coordinate Hadamards plus exact univariate r^n division reconstruct each actual sector coefficient; the latter is divisible by its indicated r^n since actual Cbar augmentation is zero',
  'actual_Xfull_source':'(B48,P8) 640 columns plus(SyzQ128,0)128; constant target1 and P0 are removed by its exact first-row identity',
  'conditional_T_source_shape_before_freezero16_cancellation':[116,936],
  'remaining_cycle_gate':'construct actual Koszul2 lifts of m*Xfull; coefficient constant/linear queries require actual V second and X first data, not arbitrary quotient zeros',
  'scalar_m_n_coefficients_or_348_T_or_K_computed':False,
  'input_hashes':{k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}}
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True)
 ep=out/'actual-fullV56-m0-and-typed-comparison.json';data=json.dumps(e,separators=(',',':'))+'\n'
 assert len(data.encode())<32<<20 and len(ops)<150000;ep.write_text(data)
 rep={'all_passed':True,'actual_fullV56_m0_rank':24,'fixed24_minor':str(det),
  'typed_actual_m_n_unit_inverse_constructed':True,'h_zero_exact':True,
  'source105_or_new_source_jet_solve':False,'full_T_cycle_lifts_or_K_torsion_computed':False,
  'operator_nodes':len(ops),'bytes':len(data.encode()),'input_hashes':e['input_hashes'],
  'evidence':{'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-t,'rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep),flush=True)
if __name__=='__main__':main()
