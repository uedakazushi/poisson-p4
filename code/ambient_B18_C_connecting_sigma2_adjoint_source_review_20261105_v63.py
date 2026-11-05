#!/usr/bin/env python3
"""Independent constant adjoint, original-z Hadamards and formal frontier."""
import hashlib,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha,mons,mon
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRE=0;LIMIT=2048
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def check(a,b,mul=False):
 global PRE
 if not a or not b:return
 na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
 h=max(na+nb,da+db) if mul else max(max(na+db,nb+da)+1,da+db)
 PRE=max(PRE,h);assert h<=LIMIT<=2048
def addq(a,b):check(a,b);return a+b
def mulq(a,b):check(a,b,True);return a*b
def total(it):
 v=Q(0)
 for x in it:v=addq(v,x)
 return v
def mat(A,B):return [[total(mulq(x,B[k][j]) for k,x in enumerate(row)) for j in range(len(B[0]))] for row in A]
def plus(a,b):
 d=dict(a)
 for k,q in b.items():
  d[k]=addq(d.get(k,Q(0)),q)
  if not d[k]:del d[k]
 return d
def scale(a,c):return {k:mulq(q,c) for k,q in a.items() if c and q}
def decode(v):return {k:Q(q) for k,q in v}
def main():
 global LIMIT
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);h0=sha(p)
  if h:assert h0==h,(str(p),h0,h)
  key=str(p.relative_to(ROOT))
  if key not in pins:RAW+=p.stat().st_size
  pins[key]=h0
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_connecting_sigma2_adjoint_source_20261105_v63.py',pr['code_sha256'],False)
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def find(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 u=find('actual-fullV56-m0-and-typed-comparison.json');p=find('all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json')
 f=find('actual-five-head-formal-RegDivell-operator-and-source-binding.json');cb=find('actual-Cbar-minimal-resolution-and-three-targets.json')
 m0=[[Q(x) for x in row] for row in field(u,'all24_by56_m0')]
 sel=field(u,'fixed_selected24_columns');iv=[[Q(x) for x in row] for row in field(u,'fixed24_minor_inverse')]
 A=[[Q(x) for x in row] for row in field(u,'fixed24_minor_matrix')]
 assert A==[[row[i] for i in sel] for row in m0]
 ops=field(u,'typed_regular_function_operators');mtyped=field(u,'typed_m24_by56')
 shifts=[Q(x) for x in field(cb,'CB_shift_coefficients7')];assert len(shifts)==7
 C=[Q(0)]*56
 for i,x in field(p,'sigma2_surviving_B_high_coordinate_adjoint'):C[i]=Q(x)
 heads=field(f,'old_five_heads_free56_rows');assert heads==[22,23,29,30,55]
 literals=[q for row in m0+iv+A for q in row]+C+shifts+[1/c for c in shifts if c]
 denominator=math.lcm(*(x.denominator for x in literals));K=max(1,max((abs(x.numerator)+x.denominator-1)//x.denominator for x in literals))
 # All computations below are degree<=4 polynomials in these literal
 # rationals (including reciprocals of KNOWN shifts), with <=57^4 terms.
 H=max(((57*K*denominator)**4).bit_length(),(denominator**4).bit_length())
 LIMIT=2*H+1;assert LIMIT<2048
 owned=count([pr,D,m0,sel,iv,A,ops,mtyped,shifts,C,heads,literals,pins])+40000
 assert owned<150000 and RAW+1000000<32<<20
 save(out/'positive-before-constant-source-checks.json',{'source_only':True,'literal_common_denominator':str(denominator),'literal_absolute_integer_bound':K,'degree_upper':4,'terms_upper':57**4,'coefficient_height_upper':H,'operation_prebound':LIMIT,'whole_numeric_upper':owned,'RAW_input_plus_reserve':RAW+1000000})
 I=[[Q(i==j) for j in range(24)] for i in range(24)]
 assert mat(A,iv)==I==mat(iv,A)
 n0=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(sel):n0[j]=iv[i]
 assert mat(m0,n0)==I
 E=[[Q(i==h) for h in heads] for i in range(56)]
 N=mat(n0,[[row[h] for h in heads] for row in m0])
 sigma=[[addq(x,-y) for x,y in zip(a,b)] for a,b in zip(E,N)]
 assert mat(m0,sigma)==[[Q(0)]*5 for _ in range(24)]
 assert D['actual_sigma0_56_by5_sparse']==[[[i,str(sigma[i][j])] for i in range(56) if sigma[i][j]] for j in range(5)]
 L=mat([C],n0)[0];sparseL=[[i,str(x)] for i,x in enumerate(L) if x]
 assert sparseL==D['fixed_L_Cn0_row24_sparse']==[[8,'-43303/59400'],[12,'-149887/237600'],[13,'30713/59400']]
 # Check EVERY surviving path back to original source-column function.
 checked=[]
 for i in [8,12,13]:
  for j in range(56):
   op=ops[mtyped[i][j]]
   if op==['zero']:continue
   assert op[0]=='regular_coordinate_Hadamard' and op[-2:]==[1,i-8]
   upper,lower=ops[op[1]],ops[op[2]]
   assert upper[:3]==lower[:3] and upper[0]=='adapted_sector_prefix' and upper[2]==1
   assert upper[3]==i-7 and lower[3]==i-8
   checked.append([i,j,mtyped[i][j]])
 assert checked==D['all_168_nonzero_m_operator_paths_checked']
 # Independent formal original-z evaluation, 3*165 monomials. Prefix k
 # retains first k y coordinates; ORIGINAL z=0 has y=c*r², not y=0.
 checks=0
 for i in [0,4,5]:
  ci=shifts[i];assert ci
  for e in mons(3):
   if e[8]:continue
   value=[]
   for k in [i+1,i]:
    v=Q(1)
    for j,n in enumerate(e[:7]):
     for _ in range(n):v=mulq(v,Q(0) if j<k else -shifts[j])
    value.append(v)
   actual=mulq(addq(value[0],-value[1]),1/ci) if e[7]+2*sum(e[:7])-2==2 else Q(0)
   expected=Q(e==mon(r=2,z=i))
   if e[i]==2 and sum(e)==2:expected=addq(expected,-ci)
   if e[i]==1 and sum(e)==2:
    for j in range(i+1,7):
     if e[j]==1:expected=addq(expected,-shifts[j])
   assert actual==expected,(i,e);checks+=1
 assert checks==495
 # Source aliases for the FULL plane4 frontier, not only its low jets.
 fdoc=pin(f,pr['input_sha256'][str(f.relative_to(ROOT))]);chaininfo=fdoc['source_namespace_bindings']['chain']
 chain=ROOT/chaininfo['path'];pin(chain,chaininfo['sha256'],False)
 identities=field(chain,'accepted_actual_regular_function_identities')
 aliases={r['physical_odd_coefficient']:r['odd_generator'] for r in identities if r['kind']=='exact_QQ_linear_source_alias'}
 tail=next(r for r in identities if r['kind']=='accepted_actual_CA_f63_tail')
 assert tail['function']==7790 and tail['gA']==505938 and len(tail['coordinate_quotients'])==7
 c2=ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
 pin(c2,'c29b3794da3410228b7620102ec1906af68d5a12df3222d0ffa31c6ca8f439a6',False)
 physical=field(c2,'physical12_coefficients');ca=[7148,7199,7454,7505,7688,7943,8045]
 for i,rec in enumerate(physical[:7]):assert aliases[rec['roots']['bA_node']]==ca[i]
 assert aliases[next(r for r in physical if r['physical_row']==63)['roots']['bA_node']]==7790
 # Universal C4 multiplication gives exact Bfull and raw plane4 primitives.
 B=[];VA=[];VB=[]
 for rec in physical:
  row=rec['physical_row']
  B.extend([{f'bB{row}':Q(-1)},{},{f'delta{row}':Q(1)},{f'nuA_times_bA{row}':Q(1)}])
  ba={} if row in [r['physical_row'] for r in physical[:7]] else {('gA4' if row==63 else f'bA4_{row}'):Q(1)}
  VA.extend([ba,{f'delta4_{row}':Q(1)},{},{f'KNOWN_nuB2_times_bB2_{row}':Q(-1)}])
  VB.extend([{f'bB4_{row}':Q(-1)},{},{f'delta4_{row}':Q(1)},{}])
 B.extend([{} for _ in range(8)]);VA.extend([{} for _ in range(7)]+[{'gA4':Q(-1)}]);VB.extend([{} for _ in range(8)])
 combos=[];fourth=[]
 for h in range(5):
  fs=[]
  for table in [B,VA,VB]:
   ans={}
   for j,c in enumerate(sigma):ans=plus(ans,scale(table[j],c[h]))
   fs.append(ans)
  assert fs[0]==decode(D['five_actual_B_function_sigma0_combinations'][h]);combos.append(fs[0])
  for k,label in [(1,'A4_sigma0_primitive'),(2,'B4_sigma0_primitive')]:assert fs[k]==decode(D['coherent_V4_sigma0_primitive_columns'][h][label])
  fourth.extend(fs[1:])
 primitive=sorted({key for poly in fourth for key in poly if not key.startswith('KNOWN_')})
 assert primitive==D['seven_remaining_primitive_plane_fourth_functions_for_V4_sigma0']
 assert len(primitive)==7 and D['fourteen_primitive_slots_upper_for_V4sigma0_only']==14
 actual_owned=count([pr,D,m0,sel,iv,A,ops,mtyped,shifts,C,heads,literals,pins,n0,E,N,sigma,L,checked,fdoc,identities,aliases,tail,physical,B,VA,VB,combos,fourth,primitive])+10000
 assert actual_owned<=owned
 assert all(sha(ROOT/p)==h for p,h in pins.items()) and RAW+1000000<32<<20
 proof={'all_constant_unit_inverse_cycle_and_L_identities':True,'L':sparseL,
  'all168_candidate_B_sector_paths_zero_or_coordinate_Hadamard_checked':True,'nonzero_Hadamard_paths':len(checked),'all495_original_z_prefix_monomials':True,
  'ordinary_recurrence_with_n1_retained':[
   'n=Esel*(mEsel)^-1 has support on the selected24 rows. Every positive-order sigma coefficient also has that support, and m*sigma=0. Its order1 and2 equations uniquely give sigma1=-n0*m1*sigma0 and sigma2=-n0*(m2*sigma0+m1*sigma1).',
   'Equivalently n1=-n0*m1*n0 and n2=n0*m1*n0*m1*n0-n0*m2*n0. Substitution gives the same formulas without assuming n1zero or m0headzero. These are ordinary coefficients, not second derivatives.'
  ],
  'original_z_Hadamard_proof':'At originalz=0, y=z+c*r2 equals c*r2. Prefixk retains firstk y and sets later yzero, so originalz_j=0 forj<k and-c_j*r2 forj>=k. Difference ofprefixi+1 andprefixi dividedbyc_i*r2 gives [r2]Hi(f)=[zi*r2]f-ci[zi2]f-sum(j>i)cj[zi*zj]f. Any monomialtotaldegree>=4 involvingzi has weighted rdegree>=3 afterdivision, soonly originaldegree<=3 is needed.',
  'L_no_curve_division':'Rows8,12,13 are B-sector z0,z4,z5 regular coordinate Hadamards. None is curve division row7/15/23. This does not remove last-row source4 from sigma1_r elsewhere.',
  'all5_exact_B_function_sigma0_combinations_and7_primitive_V4_functions':True,
  'exact_plane_frontier_source':'First7bAi are actualzi by sourcealiases; bA63 is rawf637790 and equalsgA505938 onz0 by full7telescopes. V4_A(tauAB)=-nuB2*bB2 is retained; V4_B(tauAB)=0 becausebAplane2=0 andnu0/first=0. No primitive fourth coefficient is assignedzero.',
  'fourteen_slots_scope':'Only rawV4*sigma0 uses these7functions times(r3ell,r2ell2). The FULL coherent sum also containsV3sigma1,mixedV2sigma1,Lm2sigma0,Lm1sigma1. No fullhighest cancellation, directional cache coverage, edge value or valueGO is inferred.',
  'no_new_source3_4_sigma_jet_or_connecting_values':True,
  'positive_degree4_height_bound':H,'positive_preoperation_bound':LIMIT,'actual_preoperation_bits':PRE,
  'whole_owned_upper':owned,'RAW_unique_inputs':RAW}
 ep=out/'independent-sigma2-recurrence-three-row-Hadamards-and-full-frontier.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'actual_L_paths495_prefix_checks_and14slot_SOURCE_PASS':True,
  'connecting_or_new_higher_values_computed':False,'seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'positive_pre':LIMIT,'actual_pre':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
