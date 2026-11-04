#!/usr/bin/env python3
"""Independent saved-point matrix/Schur/minor and degree3-colon proof audit."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def inv(A):
 n=len(A);a=[row[:]+[Q(int(i==j))for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  i=next(i for i in range(j,n)if a[i][j]);a[j],a[i]=a[i],a[j];p=a[j][j];a[j]=[c/p for c in a[j]]
  for i in range(n):
   if i!=j and a[i][j]:
    c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[j])]
 return[row[n:]for row in a]
def mm(A,B):return[[sum((a*b for a,b in zip(row,col)),Q())for col in zip(*B)]for row in A]
def sub(A,rows,cols):return[[A[i][j]for j in cols]for i in rows]
def det_integer(A):
 n=len(A);den=[];a=[]
 for row in A:
  d=math.lcm(*(v.denominator for v in row));den.append(d);a.append([int(v*d)for v in row])
 prior=1;sign=1
 for k in range(n-1):
  p=next(i for i in range(k,n)if a[i][k])
  if p!=k:a[p],a[k]=a[k],a[p];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    val=a[i][j]*pivot-a[i][k]*a[k][j];assert val%prior==0;a[i][j]=val//prior
   a[i][k]=0
  prior=pivot
 return Q(sign*a[-1][-1],math.prod(den))
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63';c=ROOT/'reports/S61_NONSQUARE_ELL_CUBIC_RANK10_20261104_v63'
 r1,r2=read(b/'report.json'),read(c/'report.json')
 assert r1['cubic11_rank_at_saved_point']==10 and not r1['symbolic_QQ_a_e_d_rank11_proved']
 assert r2['generic_actual_cubic11_rank_atleast10_proved']and r2['generic_rank11_or_upper10_not_claimed']
 for rec in[r1,r2]:
  for p,h in rec['input_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in r1['evidence_sha256'].items():assert sha(b/p)==h,p
 assert sha(c/'certificate.json')==r2['certificate_sha256']
 n=read(b/'native.json');cir=read(b/'circuit.json');w=read(c/'certificate.json')
 old=ROOT/'reports/S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63';oc=read(old/'circuit.json');assert oc['point']==[1,2,3]
 assert cir['point']==[1,2,3] and cir['parametric_source']==oc['source']
 assert cir['genuine_basis_rows']==list(range(41,57)) and cir['quad_source_rows']==[0,1,3,5]
 assert cir['restriction_v8']==['x0','x1','x2','0','x3','H','x1-3*x0','x4']
 assert cir['quad19_pivot_rows']==n['quad_pivot_rows'] and cir['quad19_pivot_columns']==n['quad_pivot_columns']
 assert cir['genuine_heavy5_pivot_rows']==n['restricted_heavy5_pivot_rows']
 assert cir['remaining11_genuine_rows']==n['remaining11_genuine_rows'] and cir['remaining16_pure_columns']==n['quotient16_columns']
 assert cir['symbolic_rank_upper']==11 and not cir['symbolic_rank11_proved_by_saved_oldpoint_nonzero35minor']
 assert w['same_already_computed_QQ_point']==[1,2,3] and w['generic_rank_interval_for_actual_R11']==[10,11]
 sourcep=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';source=read(sourcep);mon8=source['normal_monomials']
 lines=(old/'point-matrix.txt').read_text().splitlines();assert lines[0]=='57 91';M=[[Q(v)for v in line.split()]for line in lines[2:]];assert len(M)==57
 # Independent polynomial substitution: repeated convolution, not producer's
 # binomial formula. Six target exponents are x0,x1,x2,x3,x4,H.
 zero=(0,)*6;basis=[tuple(int(i==j)for i in range(6))for j in range(6)]
 images=[{basis[0]:Q(1)},{basis[1]:Q(1)},{basis[2]:Q(1)},{},{basis[3]:Q(1)},{basis[5]:Q(1)},{basis[1]:Q(1),basis[0]:Q(-3)},{basis[4]:Q(1)}]
 def times(p,q):
  ans={}
  for m,c0 in p.items():
   for v,d0 in q.items():
    k=tuple(a+b for a,b in zip(m,v));ans[k]=ans.get(k,Q())+c0*d0
  return{k:v for k,v in ans.items()if v}
 cache={}
 def monmap(m):
  m=tuple(m)
  if m not in cache:
   z={zero:Q(1)}
   for power,image in zip(m,images):
    for _ in range(power):z=times(z,image)
   cache[m]=z
  return cache[m]
 def restrict(entries):
  p={}
  for mon,c0 in entries:
   if not c0:continue
   for ex,d0 in monmap(mon).items():p[ex]=p.get(ex,Q())+c0*d0
  return{k:v for k,v in p.items()if v}
 pm=[tuple(m)for m in cir['restricted40_column_order']['pure35_exponents']];assert len(pm)==len(set(pm))==35
 assert all(len(m)==5 and sum(m)==3 for m in pm);pi={m:i for i,m in enumerate(pm)}
 qpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';qd=read(qpath);qrows=[]
 for qi in[0,1,3,5]:
  q=restrict([(mon,sum((Q(cc)*2**ex[1]*3**ex[2]for ex,cc in terms),Q()))for mon,terms in qd['quadratic_rows_sparse'][qi]])
  assert all(not m[5]and sum(m)==2 for m in q)
  for j in range(5):
   row=[Q()]*35
   for mon,co in times(q,{basis[j]:Q(1)}).items():row[pi[mon[:5]]]+=co
   qrows.append(row)
 G=[]
 for row in M[41:]:
  p=restrict(zip(mon8,row));z=[Q()]*40
  for mon,co in p.items():
   if not mon[5]:z[pi[mon[:5]]]+=co
   else:assert mon[5]==1 and sum(mon[:5])==1;z[35+mon[:5].index(1)]+=co
  G.append(z)
 light=[0,1,2,3,4,6,7];hcols=[]
 for j in light:
  m=[0]*8;m[5]=m[j]=1;hcols.append(mon8.index(m))
 H7=[[row[j]for j in hcols]for row in M[41:]]
 toks=iter((b/'oldpoint-native-input.txt').read_text().split());assert next(toks)=='ELLCOLON11_FROZEN_QQ_POINT'
 for actual in[qrows,G,H7]:
  rr,cc=int(next(toks)),int(next(toks));assert[rr,cc]==[len(actual),len(actual[0])]
  assert[[Q(next(toks))for _ in range(cc)]for _ in range(rr)]==actual
 assert next(toks,None)is None
 full=[row+[Q()]*5 for row in qrows]+G
 qr,qc=n['quad_pivot_rows'],n['quad_pivot_columns'];gr=n['restricted_heavy5_pivot_rows'];rg=n['remaining11_genuine_rows'];j16=n['quotient16_columns']
 pivotrows=qr+[20+i for i in gr];pivotcols=qc+list(range(35,40));assert len(pivotrows)==len(pivotcols)==24
 A24=sub(full,pivotrows,pivotcols);Ainv=inv(A24);target=sub(full,pivotrows,j16)
 def schur(rows):
  correction=mm(mm(sub(full,rows,pivotcols),Ainv),target)
  return[[x-y for x,y in zip(row,co)]for row,co in zip(sub(full,rows,j16),correction)]
 assert all(not v for row in schur(list(range(20)))for v in row)
 S=schur([20+i for i in rg]);assert S==[[Q(v)for v in row]for row in n['actual11x16']]
 tenrows=w['extracted10_rows'];tencols=w['fixed_columns_from_completed_RREF'];minor=sub(S,tenrows,tencols)
 assert minor==[[Q(v)for v in row]for row in w['exact10x10_minor']]
 delta10=det_integer(minor);assert delta10==Q(w['exact_nonzero10_minor_determinant'])and delta10
 # Upper10 at THIS point only: all11 rows reconstruct from the saved10.
 coeff=mm(sub(S,range(11),tencols),inv(minor));assert mm(coeff,sub(S,tenrows,range(16)))==S
 delta19=det_integer(sub(qrows,qr,qc));delta5=det_integer(sub(G,gr,range(35,40)));delta7=det_integer(sub(H7,n['heavy7_pivot_rows'],range(7)))
 assert[delta19,delta5,delta7]==[Q(n[k])for k in['quad19_minor','heavy5_minor','heavy7_minor']]
 p34=w['parametric_nonzero34minor'];delta34=det_integer(sub(full,p34['rows'],p34['columns']))
 assert delta34==Q(p34['oldpoint_value'])==delta19*delta5*delta10
 proof=dict(actual1452_QQ_source_entries_independently_reconstructed=True,
  independent_restriction='repeated sparse convolution v6=x1-3x0,v3=0, preservesell; not author binomial helper',
  actual11x16=[[str(v)for v in row]for row in S],
  fixed10_minor_determinant=str(delta10),fixed34_minor_determinant=str(delta34),
  exact_all11_rows_in_saved10_rowspan_coefficients=[[str(v)for v in row]for row in coeff],
  saved_point_rank10_only=True,generic_actual_rank_interval=[10,11],
  conceptual_colon='OnU78, Koszul6isWHOLEkernel ofraw84->126; lift/subtract6, allcofactors becomeell-divisible; cancelell inpolynomialring ->(J:ell)3=J3',
  exact_extension='0->(S/(J,Z,W))3 rank16 ->ell (S/I)4 rank64 ->(S/(I,ell))4 rank48 ->0',
  corrected_source_lifts=['ell*r9','H*(ellZ)-ell*fZ=ell*pZ','H*(ellW)-ell*fW=ell*pW'],
  canonical11_in_C122_source_split=True,target_rank11_not_assumed=True,
  reduction='OnU78,heavy7,Delta10 andontoV48chart, quotientkernel15/sourceimage10, then48projection gives sameactualcoker E49->V6',
  fullCI_regular_sequence_or_ell_nonzerodivisor_modJ_not_used=True,
  projection48_not_proved_by_this_packet=True,no_allpoint57_full239_or_global_unit_claim=True)
 (out/'independent-certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 deps=[ROOT/'reports'/z/'report.json'for z in['S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW','S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v63_REVIEW','S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63_REVIEW','S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63_REVIEW','S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW','S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in deps)
 note=ROOT/f'research/notes/{ID}.md'
 result=dict(id=ID,all_passed=True,degree3_colon_and_exact16_64_48_bundle_sequence_accepted=True,
  actual_canonical11_source_and_required_heavy_commutator_corrections_accepted=True,
  independent_saved_point1452_entry_binding_and_full11x16_Schur_match=True,
  saved_QQ_point_exact_cubic_rank10_not_generic_upper10=True,
  integer_Bareiss_19_5_7_10_34_minors_independently_checked=True,
  actual_parametric_Schur_minor_gives_generic_rank_atleast10=True,
  conditional49_to6_same_cokernel_recipe_accepted=True,projection48_or_generic_rank11_not_proved=True,
  no_new_point_prime_GB_field_solver_fullquartic_expansion_or_author_rerun=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,b/'report.json',b/'native.json',b/'circuit.json',b/'oldpoint-native-input.txt',c/'report.json',c/'certificate.json',old/'point-matrix.txt',sourcep,qpath,*deps]},
  independent_certificate_sha256=sha(out/'independent-certificate.json'),seconds=time.monotonic()-started,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
