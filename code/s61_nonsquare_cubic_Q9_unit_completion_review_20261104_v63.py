#!/usr/bin/env python3
"""Independent sparse free-algebra review of the actual GL10 target completion."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE='S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63';ID=BASE+'_REVIEW'
ZERO={};ONE={():Q(1)};DEADLINE=0;MAXMONS=0
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def constant(c):return{}if not c else{():Q(c)}
def variable(name,i=0):return{((name,i),):Q(1)}
def add(a,b):
 global MAXMONS
 z=a.copy()
 for m,c in b.items():
  v=z.get(m,Q())+c
  if v:z[m]=v
  else:z.pop(m,None)
 MAXMONS=max(MAXMONS,len(z));assert len(z)<20000 and time.monotonic()<DEADLINE;return z
def neg(a):return{m:-c for m,c in a.items()}
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 global MAXMONS
 if not a or not b:return{}
 z={}
 for m,c in a.items():
  for n,d in b.items():
   k=tuple(sorted(m+n));v=z.get(k,Q())+c*d
   if v:z[k]=v
   else:z.pop(k,None)
 MAXMONS=max(MAXMONS,len(z));assert len(z)<20000 and time.monotonic()<DEADLINE;return z
def eye(n):return[[ONE if i==j else ZERO for j in range(n)]for i in range(n)]
def mm(A,B):
 n,k,m=len(A),len(B),len(B[0]);assert all(len(row)==k for row in A);z=[[{}for j in range(m)]for i in range(n)]
 for i in range(n):
  for t in range(k):
   if A[i][t]:
    for j in range(m):
     if B[t][j]:z[i][j]=add(z[i][j],mul(A[i][t],B[t][j]))
 return z
def matsub(A,B):return[[sub(x,y)for x,y in zip(a,b)]for a,b in zip(A,B)]
def scaled(c,A):return[[mul(c,x)for x in row]for row in A]
def summands(values):
 z={}
 for q in values:z=add(z,q)
 return z
def universal():
 h=[variable('h',i)for i in range(9)];c=[variable('c',i)for i in range(9)];t=summands(mul(x,y)for x,y in zip(h,c));B=[[sub(ONE if i==j else ZERO,mul(c[i],h[j]))for j in range(9)]for i in range(9)]
 U=[B[i]+[c[i]]for i in range(9)]+[h+[ZERO]];W=[[mul(c[i],h[j])for j in range(9)]+[neg(c[i])]for i in range(9)]+[[neg(x)for x in h]+[ONE]]
 assert matsub(mm(U,U),eye(10))==scaled(sub(t,ONE),W)
 assert mm([h+[ZERO]],U)==[[mul(sub(ONE,t),x)for x in h]+[t]]
 Adj=[[add(sub(ONE,t)if i==j else ZERO,mul(c[i],h[j]))for j in range(9)]for i in range(9)]
 diagonal=scaled(sub(ONE,t),eye(9));assert mm(B,Adj)==mm(Adj,B)==diagonal
 # Rank-one minors make every determinant term with>=2 changed columns zero.
 rankone_checks=0
 for i,k in itertools.combinations(range(9),2):
  for j,l in itertools.combinations(range(9),2):assert mul(mul(c[i],h[j]),mul(c[k],h[l]))==mul(mul(c[i],h[l]),mul(c[k],h[j]));rankone_checks+=1
 bordered=neg(summands(mul(mul(h[i],Adj[i][j]),c[j])for i in range(9)for j in range(9)));assert bordered==neg(t)
 return{'all100_involution_remainders_checked':True,'all10_row_coefficients_checked':True,'all162_adjugate_products_checked':True,'all1296_rankone2minors_checked':rankone_checks==1296,
  'det_I_minus_ch':'1-t by column multilinearity; >=2 changed columns have rank-one minors zero','adjugate':'(1-t)I+ch; verified products and det identify it as a polynomial matrix over the universal domain','bordered_det_U':'-h*adj(B)*c=-t','actual_t1_det_U':'-1'}
def main(out):
 global DEADLINE
 start=time.monotonic();DEADLINE=start+50;signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports'/BASE;pr=read(base/'report.json');cert=read(base/'certificate.json');assert sha(base/'certificate.json')==pr['certificate_sha256'];hashes={}
 for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
 powerpath=ROOT/cert['scalar_DAG']['external_power_certificate'];matrixpath=ROOT/cert['scalar_DAG']['external_matrix_certificate'];power=read(powerpath);parent=read(matrixpath)
 gate=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63_REVIEW/report.json';peer=read(gate);assert peer['all_passed']and peer['frozen'];assert peer['input_sha256'][str(powerpath.relative_to(ROOT))]==sha(powerpath)and peer['input_sha256'][str(matrixpath.relative_to(ROOT))]==sha(matrixpath)
 pg=cert['accepted_power_inner_product_gate'];assert pg['power_certificate_sha256']==sha(powerpath)and pg['independent_review_sha256']==sha(gate)
 oldrows=parent['new15_source_original25_labels'];oldcols=parent['new22_target_original32_columns'];Aold=parent['remaining15x22_scalar_roots'];assert len(oldrows)==15 and oldrows[0]==9 and len(oldcols)==22 and len(Aold)==15
 selected=cert['selected_original_target_labels'];assert selected==power['actual_h9_original_target_labels']==[11,13,14,15,17,18,19,22,25];pos=[oldcols.index(j)for j in selected];spare=cert['pure_spare_original_target_label'];tails=[j for j in oldcols if j not in selected];assert tails==cert['all13_tail_original_labels']and len(tails)==13 and spare==1 and spare in tails
 heavy=[27,28,29,30,31];assert spare not in heavy and all(Aold[0][oldcols.index(j)]==0 for j in heavy)
 mallowed={r for row in Aold for r in row if r};pallowed={r for r in power['actual_h9_cofactor_roots']if r};nodes=cert['scalar_DAG']['nodes'];assert len(nodes)==pr['new_scalar_nodes']==6880
 values=[];matrix_ids={};power_ids={};inverse_ids=[]
 for i,n in enumerate(nodes):
  tag,*a=n
  if tag=='z':v=ZERO
  elif tag=='one':v=ONE
  elif tag=='known_inverse_a29':assert not a;v=variable('u29');inverse_ids.append(i)
  elif tag=='matrix_parent':assert len(a)==1 and a[0]in mallowed and a[0]not in matrix_ids;matrix_ids[a[0]]=i;v=variable('M',a[0])
  elif tag=='power_parent':assert len(a)==1 and a[0]in pallowed and a[0]not in power_ids;power_ids[a[0]]=i;v=variable('H',a[0])
  else:
   assert tag in['a','s','m']and len(a)==2 and all(0<=j<i for j in a);v=add(values[a[0]],values[a[1]])if tag=='a'else sub(values[a[0]],values[a[1]])if tag=='s'else mul(values[a[0]],values[a[1]])
  values.append(v)
 assert len(inverse_ids)==1 and set(power_ids)==pallowed and set(matrix_ids)==mallowed
 get=lambda A:[[values[q]for q in row]for row in A]
 A=[[ZERO if not r else variable('M',r)for r in row]for row in Aold];h=[A[0][j]for j in pos];cc=[ZERO if not r else mul(variable('H',r),variable('u29'))for r in power['actual_h9_cofactor_roots']]
 assert[values[r]for r in cert['actual_h9_roots']]==h and[values[r]for r in cert['actual_c9_dual_roots']]==cc
 t=summands(mul(x,y)for x,y in zip(h,cc));assert values[pg['inner_root']]==t
 univ=universal();T=eye(22);Ti=eye(22)
 for label in tails:
  j=oldcols.index(label);b=A[0][j]
  for i,ci in zip(pos,cc):T[i][j]=neg(mul(ci,b));Ti[i][j]=mul(ci,b)
 assert mm(T,Ti)==mm(Ti,T)==eye(22);assert get(cert['tail_clear22_roots'])==T and get(cert['tail_clear22_inverse_roots'])==Ti
 aftertail=mm(A,T)
 for j,label in enumerate(oldcols):assert aftertail[0][j]==(A[0][j]if label in selected else mul(sub(ONE,t),A[0][j]))
 U=[[sub(ONE if i==j else ZERO,mul(cc[i],h[j]))for j in range(9)]+[cc[i]]for i in range(9)]+[h+[ZERO]];assert get(cert['GL10_roots'])==U
 block=selected+[spare];assert block==cert['GL10_original_target_block_order'];positions=[oldcols.index(j)for j in block];U22=eye(22)
 for i,ri in enumerate(positions):
  for j,cj in enumerate(positions):U22[ri][cj]=U[i][j]
 target=mm(T,U22);inverse=mm(U22,Ti);assert get(cert['target22_roots'])==target and get(cert['target22_inverse_roots'])==inverse
 AT=mm(A,target);pivot=oldcols.index(spare);assert[values[r]for r in cert['pivot_row_before_exact_simplification_roots']]==AT[0]
 # Exact first-row remainder factor, before invoking the accepted actual t=1 gate.
 remainder=[]
 for j,label in enumerate(oldcols):
  if label in selected:g=mul(add(ONE,A[0][pivot]),A[0][j])
  elif label==spare:g=neg(ONE)
  else:g=A[0][j]
  assert sub(AT[0][j],ONE if j==pivot else ZERO)==mul(sub(ONE,t),g);remainder.append(g)
 assert cert['pivot_row_exact_simplified']==[int(j==pivot)for j in range(22)]
 source=eye(15);sourceinverse=eye(15)
 for i in range(1,15):source[i][0]=neg(AT[i][pivot]);sourceinverse[i][0]=AT[i][pivot]
 assert get(cert['source15_roots'])==source and mm(source,sourceinverse)==mm(sourceinverse,source)==eye(15)
 actual=mm(source,AT)
 for i in range(15):
  for j in range(22):
   wanted=(ONE if j==pivot else ZERO)if i==0 else ZERO if j==pivot else AT[i][j]
   g=remainder[j]if i==0 else AT[i][pivot]if j==pivot else neg(mul(AT[i][pivot],remainder[j]))
   assert sub(actual[i][j],wanted)==mul(sub(ONE,t),g)
 rows=oldrows[1:];cols=[j for j in oldcols if j!=spare];F=[[AT[i][oldcols.index(j)]for j in cols]for i in range(1,15)]
 assert get(cert['actual14x21_scalar_roots'])==F and cert['old15_source_labels']==oldrows and cert['old22_target_labels']==oldcols and cert['new14_source_labels']==rows and cert['new21_target_labels']==cols
 assert len(F)==14 and all(len(row)==21 for row in F)
 for i in range(1,15):
  for j in heavy:
   k=oldcols.index(j);assert AT[i][k]==A[i][k];r=Aold[i][k];expected=0 if not r else matrix_ids[r];assert cert['actual14x21_scalar_roots'][i-1][cols.index(j)]==expected
 assert cert['all5_heavy_quotient_labels_unchanged']==heavy and cert['source_change_determinant']=='1'and cert['target_change_determinant']=='-1'and cert['total_global_unit_pivots_in_augmented25_presentation']==11
 proof={'own_sparse_free_commutative_Fraction_algebra':True,'all6880_nodes_and_external_root_labels_replayed':True,'actual_POWER_gate_bound_to_both_actual9_and_current_parent15x22':True,
  'new_inverse_only_certified_a29':True,'universal_GL10':univ,'all13_tail_changes_and_both_inverse_products_checked':True,'spare_original_pure_target1_verified':True,
  'embedded_GL10_target22_and_explicit_inverse_factorization_checked':True,'all22_first_row_remainders_explicitly_divisible_by1_minus_t':True,'all330_full_source_clearing_identities_mod_actual_dot1_checked':True,
  'all294_final_entries_and_all_source_target_labels_checked':True,'all70_heavy_entries_preserved_at_literal_root_level':True,'actual14x21_same_cokernel_accepted':True,
  'no_power207_coefficient_expansion_authorbuilder_point_GB_or_extra_completion':True,'no_full57_full239_original_fixed11minor_or_kernel_count_claim':True,'max_free_monomials_in_one_scalar':MAXMONS}
 (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',gate]:hashes[str(p.relative_to(ROOT))]=sha(p)
 report={'id':ID,'all_passed':True,'actual_global15x22_to14x21_same_cokernel_independently_accepted':True,'all294_entries_and70_heavy_literalroots_verified':True,'input_sha256':hashes,'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},'seconds':time.monotonic()-start,'max_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
 signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'all_passed':True,'seconds':report['seconds'],'max_RSS_KiB':report['max_RSS_KiB'],'report_sha256':sha(out/'report.json')},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
