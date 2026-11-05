#!/usr/bin/env python3
"""Independent fixed linear elimination, no producer Groebner division."""
import json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63';ID=PRI+'_REVIEW';PRE=0
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def op(a,b,m=False):
 global PRE
 if a and b:
  na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
  h=max(na+nb,da+db)if m else max(na+db,nb+da)+1;PRE=max(PRE,h);assert h<=1153
 q=a*b if m else a+b;assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=384;return q
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text())if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json','09091109926dd0f4cdcfb6a519df11588209924f1d0dec1063d1ad4d2ebe3c32')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_Kbar_r2_initial_source_20261105_v63.py',pr['code_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 ownbasis=pin(ROOT/'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63_REVIEW/report.json','a5ec34be9b391a8dc9b6a34ce242c962046249301faab7c43deea42cc5cfe71a');assert ownbasis['all_passed']
 def src(s):return ROOT/next(p for p in pins if p.endswith(s))
 old=src('known-initial40-and-length-upper17.json');g=field(old,'homogeneous_module_generators40')
 linear=field(src('linear-module-rewrite-certificate.json'),'linear_basis37')
 exact=src('exact-Loewy11-and-source-basis-corollary.json');assert field(exact,'exact_Loewy_HF')==[5,3,2,1]
 assert field(old,'standard_degree1')==['e3*r','e4*z5','e4*r']and field(old,'standard_degree2')==['e4*z5*r','e4*r²']
 def unit(i):return tuple(int(j==i)for j in range(8))
 for a,row in zip(g[:37],linear):assert {(t,tuple(m)):Q(q)for t,m,q in a}=={(i//8,unit(i%8)):Q(q)for i,q in row}
 R=[[Q(0)]*40 for _ in range(37)]
 for j,row in enumerate(linear):
  for i,q in row:R[j][i]=Q(q)
 standard=[31,37,39];pivots=[i for i in range(40)if i not in standard]
 A=[[R[j][i]for j in range(37)]+[Q(i==15),Q(i==23)]for i in pivots]
 K=[[Q(q)for q in row]for row in field(src('actual-seven-barbasis-typed56-cycle-operator-and-full11-unit-complement.json'),'actual_Kbar11_by7')]
 owned=count([pr,D,ownbasis,g,linear,R,A,K,pins])+15000;assert owned<150000 and raw+200000<32<<20
 save(out/'positive-before-fixed37_linear_elimination.json',{'whole_SOURCE_upper':owned,'linear37_transpose_system_two_fixed_RHS':True,'no_GB_or_new_source_values':True,'known_degree1_standards':standard,'committed384_temporary1153':True})
 # Solve source cofactor coordinates with bottom-most row pivots, unlike
 # the producer's monomial-leading division algorithm.
 for j in range(37):
  i=next(i for i in reversed(range(j,37))if A[i][j])
  if i!=j:A[j],A[i]=A[i],A[j]
  p=A[j][j];A[j]=[op(q,1/p,True)for q in A[j]]
  for i in reversed(range(37)):
   if i==j or not A[i][j]:continue
   ratio=A[i][j];A[i]=[op(a,-op(ratio,b,True))for a,b in zip(A[i],A[j])]
 cofs=[[A[j][37+k]for j in range(37)]for k in range(2)];remainders=[]
 for target,cf in zip([15,23],cofs):
  row=[Q(i==target)for i in range(40)]
  for j in reversed(range(37)):
   for i in reversed(range(40)):row[i]=op(row[i],-op(cf[j],R[j][i],True))
  assert all(not q for i,q in enumerate(row)if i not in standard);remainders.append(row)
 assert remainders[0][31]==remainders[0][39]==0 and remainders[1][31]==Q(1,4)and remainders[1][39]==0
 alpha1=remainders[0][37];alpha2=remainders[1][37]
 assert [alpha1,alpha2]==list(map(Q,D['actual_alpha1_alpha2']))==[Q(-1681,1992),Q(75645,14608)]
 # Independent multiplyback of ALL supplied polynomial certificates.
 gens=[{(t,tuple(m)):Q(q)for t,m,q in row}for row in g]
 def decode(v):return{(t,tuple(m)):Q(q)for t,m,q in v}
 def reconstruct(receipt):
  z=decode(receipt['normal_form'])
  key='exact_linear_source_cofactors'if'exact_linear_source_cofactors'in receipt else'exact_old40_source_cofactors'
  for j,shift,q in reversed(receipt[key]):
   for (t,m),a in gens[j].items():
    k=(t,tuple(x+y for x,y in zip(m,shift)));z[k]=op(z.get(k,Q(0)),op(Q(q),a,True))
    if not z[k]:del z[k]
  assert z==decode(receipt.get('input',receipt.get('degree2_input')))
 for rec in D['two_exact_linear_reductions_and_membership_cofactors']+D['all7_exact_degree2_reductions_and_membership_cofactors']:reconstruct(rec)
 c=Q(11767,3652);rr=tuple(2*int(i==7)for i in range(8));z5r=tuple(int(i in[5,7])for i in range(8))
 assert g[38]==[[3,list(rr),'1'],[4,list(z5r),'11767/3652']]
 # Actual gr2 map uses r applied to the two independently solved gr1
 # equations and G38; all six higher oldbasis rows map into oldm3.
 matrix=[[op(op(alpha1,K[1][j],True),op(op(alpha2,-c/4),K[2][j],True))for j in range(7)],K[4].copy()]
 for j in range(7):matrix[0][j]=op(matrix[0][j],-op(c,K[3][j],True))
 assert matrix==[[Q(q)for q in row]for row in D['actual_projected_r2_matrix2_by7']]
 assert matrix[0][:5]==[0]*5 and matrix[1]==[0]*7
 assert matrix[0][5]==Q(273412940423,2952943290)and matrix[0][6]==Q(-6957653957,295294329)
 assert D['projected_rank']==1 and D['fixed_nonzero_witness_column']==5
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'same37_linear_source_generators_rechecked':True,'independent_fixed_linear_elimination_not_producer_GBdivision':True,
  'two_independent_source_cofactor_vectors':list(map(lambda r:list(map(str,r)),cofs)),
  'actual_gr1_alpha1_alpha2':[str(alpha1),str(alpha2)],'all_two_plus_seven_saved_cofactor_multiplybacks_exact':True,
  'actual_projected_r2_matrix2by7':[[str(q)for q in row]for row in matrix],
  'exact_HF2_completeness':'The two q1,q2 standard classes are actual gr2 Told because the source initial relations give their span and accepted exactHF2=2 gives independence. Linear reductions are gr1 identities; multiplying by r is valid modulo oldm3. Full actual r2 images may carry additional oldm3=gdeep components, kept uncomputed.',
  'projection_rank1_fixed_nonzero_U_witness':True,
  'one_valid_r_block_consequence':'Actual r2 Kbar is nonzero by its nonzero actual gr2 projection. Independently proved r3 Kbar=0 then forces at least one length3 nilpotent r block. No full r-Jordan partition or ell-DVR exponent follows.',
  'not_full_r2_action_matrix_or_higher_tail_values':True,'no_new_source_jet_or_H0_nil_claim':True,'owned_upper':owned,'actual_prebits':PRE,'RAW':raw}
 ep=out/'independent-linear-system-alphas-cofactor-multiplybacks-and-projected-r2.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'actual_projected_r2_rank1':True,'full_r2_or_Jordan_type_claimed':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
