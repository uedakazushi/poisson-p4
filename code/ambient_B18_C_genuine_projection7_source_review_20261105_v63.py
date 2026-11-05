#!/usr/bin/env python3
"""Independent local-reduction audit of the actual D-only sparse q7."""
import hashlib,itertools,json,math,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63'
ID=PRI+'_REVIEW'
PIN='de02c3e5fc8b542298127b208099e9513a91def5ade9e5f3db79359416b01341'
ZERO=(0,)*9; C=(); PRE=0
def sha(p):
 h=hashlib.sha256()
 with pathlib.Path(p).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def count(x):
 if isinstance(x,bool) or x is None:return 0
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(tuple,list)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.keys()))+sum(map(count,x.values()))
 return 0
def check(a,b,mul=False):
 global PRE
 if not a or not b:return
 na,da=abs(a.numerator).bit_length(),a.denominator.bit_length()
 nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
 h=max(na+nb,da+db) if mul else max(max(na+db,nb+da)+1,da+db)
 PRE=max(PRE,h);assert h<=2048
def qa(a,b):check(a,b);return a+b
def qm(a,b):check(a,b,True);return a*b
def qsum(it):
 v=Q(0)
 for q in it:v=qa(v,q)
 return v
def plus(a,b):
 d=dict(a)
 for k,q in b.items():
  d[k]=qa(d.get(k,Q(0)),q)
  if not d[k]:del d[k]
 return d
def scale(a,q):return {k:qm(v,q) for k,v in a.items() if v and q}
def scalar(q):return {C:Q(q)} if q else {}
def decode(a):return {tuple(k) if k!=['QQ'] else C:Q(v) for k,v in a}
def leaf(row,k):return {('cache',row,k):Q(1)}
def times(a,b):
 if not a or not b:return {}
 if set(a)<={C}:return scale(b,a.get(C,Q(0)))
 assert set(b)<={C},'No multiplication of two unpaired cache leaves'
 return scale(a,b.get(C,Q(0)))
def vadd(a,b):
 d={i:dict(v) for i,v in a.items()}
 for i,v in b.items():
  d[i]=plus(d.get(i,{}),v)
  if not d[i]:del d[i]
 return d
def vmul(a,q):return {i:times(v,q) for i,v in a.items() if times(v,q)}
def unit(i,q=1):return {i:scalar(q)} if q else {}
def mon(r=0,ell=0,z=None):
 v=[0]*9;v[7]=r;v[8]=ell
 if z is not None:v[z]+=1
 return tuple(v)
def mm(a,b):return tuple(x+y for x,y in zip(a,b))
def mons(d):
 for n in range(d+1):
  for ids in itertools.combinations_with_replacement(range(9),n):
   v=[0]*9
   for i in ids:v[i]+=1
   yield tuple(v)
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID
 if out.exists():assert {p.name for p in out.iterdir()}=={'ATTEMPT1_HISTORY'}
 else:out.mkdir()
 pins={};raw=0
 def pin(p,h=None,decode_json=True):
  nonlocal raw
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if decode_json else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PIN)
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 primary_code=ROOT/'code/ambient_B18_C_genuine_projection7_source_v2_20261105_v63.py'
 pin(primary_code,pr['code_sha256'],False)
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'],False)
 evname=next(iter(pr['evidence_sha256']));d=pin(ROOT/'reports'/PRI/evname,pr['evidence_sha256'][evname])
 # Own top-level lexer, previously frozen; no author producer imports.
 paths=list(pr['input_sha256'])
 def path(fragment):return ROOT/next(p for p in paths if fragment in p and not p.endswith('/report.json'))
 t7p=path('actual-filtered23-');bp=path('actual-special-fiber20-');fp=path('full-HB6-');hp=path('two-coherent-HB-')
 cp=path('actual-physical12-');rp=path('actual-ell2-product33-');frp=path('actual-first21-')
 P=[[[decode(x) for x in row] for row in field(fp,k)] for k in ['P3_A','P3_B']]
 basis=field(t7p,'literal_filtered23_basis_F16_then_t7');q23=field(t7p,'actual_D_linear_quotient_q7_by23')
 assert q23==[[int(j==16+i) for j in range(23)] for i in range(7)]
 A0=[[Q(x) for x in row] for row in field(bp,'A0')];Ai=[[Q(x) for x in row] for row in field(bp,'A0_inverse')]
 assert Q(field(bp,'A0_determinant'))==Q(-1,1024)
 for i in range(7):
  for j in range(7):
   assert qsum(qm(A0[i][k],Ai[k][j]) for k in range(7))==int(i==j)
 norm=field(bp,'actual_generator_normalization')
 assert 'A(z,r,s)^-1' in norm['IB'] and 'A0^-1*(b3-A1*A0^-1*b2)'==norm['beta3']
 B=[]
 for degree,key in [(2,'beta2_before_ell'),(3,'beta3_before_ell')]:
  nr=[]
  for row in field(bp,key):
   v=[Q(0)]*(degree+1)
   for e,c in row:
    assert sum(e)==degree
    if any(e[:7]):
     assert degree==3 # beta2 is plane; beta3 nonplane terms vanish by order.
     continue
    for j in range(e[8]+1):v[j]=qa(v[j],qm(Q(c),Q(math.comb(e[8],j),5**j)))
   nr.append(v)
  B.append(nr)
 beta,beta3=B
 assert beta==[[Q(x) for x in row] for row in d['beta2_plane_after_actual_s_equals_r_plus_ell_over5']]
 assert beta3==[[Q(x) for x in row] for row in d['beta3_plane_after_actual_s_equals_r_plus_ell_over5']]
 assert [i for i,row in enumerate(beta) if any(row)]==[0,4,5]
 raw33=field(rp,'actual33_basis');phys={x['physical_row']:dict(x['delta']) for x in field(cp,'physical12_coefficients')}
 recs=field(frp,'all12_verified_source_records')[:7];rows=[x['physical_row'] for x in recs]
 assert rows==[18,19,28,29,60,67,72]
 # Independently build FULL relevant RH mixed and cubic coefficients.
 mixed=[];D=[];literal=[]
 for i,row in enumerate(rows):
  rec=recs[i];assert rec['delta_first9_zero'] and rec['delta_and_bA_pureplane_second_zero']
  assert [Q(x) for x in rec['bA_first_z7']]==[Q(i==j) for j in range(7)]
  use=[]
  for z in [0,4,5]:
   for s in range(2):use.append((z,s,raw33.index([z+1,1-s,s])))
  if i not in [0,4,5]:
   for ix in [29,30]+[ix for z,s,ix in use]:
    assert Q(phys[row].get(ix,'0'))==0;literal.append(('cache',row,ix))
   mixed.append([]);D.append([{},{}])
  else:mixed.append([(z,s,leaf(row,ix)) for z,s,ix in use]);D.append([leaf(row,29),leaf(row,30)])
 a,b=Q(8,25),Q(96,275)
 kap=[[plus(scale(P[side][k][0],ab),scale(P[side][k][1],-Q(1))) for k in range(3)] for side,ab in enumerate([a,b])]
 savedkap=[[[decode(x) for x in row] for row in d['kappa_k_A_B_UNPAIRED']]][0]
 assert kap==savedkap
 for item in field(t7p,'four_ell_coefficients_UNPAIRED'):
  side=0 if item['label'].endswith('_A') else 1;k=1 if item['label'].startswith('alpha') else 2
  assert decode(item['operator_old_C2_only'])==kap[side][k]
 # Casewise actual reductions, independent of author's sparse-table builder.
 # z_i B uses full beta2/beta3. z_i H uses FULL mixed RH followed by RB,
 # so the negative mixed*beta2 carry arises from recursion, not a copied E3.
 def reduce(t,e):
  n=sum(e)
  if n>=4:return {}
  if any(e[:7]):
   z=next(i for i in range(7) if e[i]);rest=list(e);rest[z]-=1;rest=tuple(rest)
   if t==0:return {}
   if t==1:
    ans={}
    for deg,br in [(2,beta[z]),(3,beta3[z])]:
     for j,c in enumerate(br):ans=vadd(ans,vmul(reduce(1,mm(rest,mon(deg-j,j))),scalar(-c)))
    return ans
   if n>=2:return {} # z_i H lies in actual m^3C.
   ans={}
   for zz,s,c in mixed[z]:ans=vadd(ans,vmul(reduce(1,mon(1-s,s,zz)),scale(c,-Q(1))))
   for j,c in enumerate(D[z]):ans=vadd(ans,vmul(reduce(1,mon(3-j,j)),scale(c,-Q(1))))
   return ans
  r,l=e[7:]
  if t<2:
   if (r,l)==(2,0):return unit(3+t)
   if (r,l)==(2,1):return unit(5+t)
   if (r,l)==(3,0):return unit(5+t,-[a,b][t])
   return {} # all other pure-plane A/B monomials are in the actual F.
  if n<=1:return unit({(0,0):0,(1,0):1,(0,1):2}[(r,l)])
  if n==3:return {}
  k=l;ans={}
  for side in range(2):
   for j in range(2):ans=vadd(ans,vmul(reduce(side,mon(3-j,j)),scale(P[side][k][j],-Q(1))))
  return ans
 saved={(x['raw_target'],tuple(x['monomial'])):{i:decode(v) for i,v in x['target7_affine']} for x in d['projector7_by660_sparse_columns']}
 allm=list(mons(3));assert len(allm)==220
 rebuilt={}
 for t in range(3):
  for e in allm:
   v=reduce(t,e);assert v==saved.get((t,e),{}),(t,e)
   if v:rebuilt[t,e]=v
 assert len(rebuilt)==25
 # Source normalization P=T^t(D,H*E), retaining mixed RH*RB carry.
 E=[]
 for i in range(7):
  es=[dict(x) for x in D[i]]
  for z,s,c in mixed[i]:
   for j in range(2-s):es[j+s]=plus(es[j+s],scale(c,-beta[z][j]))
  E.append(es)
 assert E==[[decode(x) for x in row] for row in d['first7_effective_E3_first2_UNPAIRED']]
 T=[[Q(x) for x in row] for row in field(hp,'constant_source_normalization_T')]
 H=[[Q(x) for x in row] for row in field(hp,'source_first7_h_rows')]
 for k in range(3):
  for j in range(2):
   pA={};pB={}
   for u,i in enumerate([0,4,5]):
    pA=plus(pA,scale(D[i][j],T[u][k]))
    h=qsum(qm(T[w][k],H[w][i]) for w in range(3));pB=plus(pB,scale(E[i][j],h))
   assert pA==P[0][k][j] and pB==P[1][k][j]
 # Build original normalized relation families and stream ALL allowed multiples.
 rel=[]
 for i in range(7):
  rel.append(('RA'+str(i),1,[(0,mon(z=i),scalar(1))]))
  rb=[(1,mon(z=i),scalar(1))]
  for deg,br in [(2,beta[i]),(3,beta3[i])]:rb += [(1,mon(deg-j,j),scalar(c)) for j,c in enumerate(br)]
  rel.append(('RB'+str(i),1,rb))
  rh=[(2,mon(z=i),scalar(1))]+[(1,mon(1-s,s,z),c) for z,s,c in mixed[i]]
  rh += [(1,mon(3-j,j),c) for j,c in enumerate(D[i])]
  rel.append(('RH'+str(i),1,rh))
 for k in range(3):
  terms=[(2,mon(2-k,k),scalar(1))]
  terms += [(side,mon(3-j,j),P[side][k][j]) for side in range(2) for j in range(2)]
  rel.append(('P'+str(k),2,terms))
 for side,ab in enumerate([a,b]):rel.append(('curve'+str(side),3,[(side,mon(3),scalar(1)),(side,mon(2,1),scalar(ab))]))
 labels=[];hasher=hashlib.sha256()
 for name,n,terms in rel:
  done=0
  for mul in mons(3-n):
   v={}
   for t,e,c in terms:v=vadd(v,vmul(reduce(t,mm(e,mul)),c))
   assert not v,(name,mul);done+=1
   hasher.update(json.dumps([name,list(mul)],separators=(',',':')).encode())
  labels.append([name,n,done])
 assert labels==d['normalized26_source_relation_multiplier_checks'] and sum(x[2] for x in labels)==1187
 assert hasher.hexdigest()==d['source_multiplier_labels_hash']
 def Eact(v):
  w={}
  for i,c in v.items():
   if i in [0,3,4]:w=vadd(w,{ {0:2,3:5,4:6}[i]:c})
   if i in [1,2]:
    assert set(c)<={C}
    for side in range(2):w=vadd(w,vmul({5+side:kap[side][i]},c))
  return w
 for t in range(3):
  for e in allm:assert reduce(t,mm(e,mon(ell=1)))==Eact(reduce(t,e))
 basis_checks=0
 for j,item in enumerate(basis):
  if item['m_degree']>=4:continue
  assert reduce(item['raw_target'],mon(item['r_power'],item['ell_power']))==(unit(j-16) if j>=16 else {})
  basis_checks+=1
 assert basis_checks==21
 # Exact V1 history: one field label only. Do not mutate either producer.
 hist=ROOT/'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_HISTORY_20261105_v63'
 h=pin(hist/'v1-field-label-failure.json');old=hist/'source-v1-before-field-label-correction.py'
 pin(old,h['old_code_sha256'],False)
 assert old.read_text().replace("'actual_four_coefficients'","'all4_actual_cached_coefficients'")==primary_code.read_text()
 # Recheck every frozen author/dependency byte after the independent work.
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 owned=count([pr,d,P,B,A0,Ai,norm,basis,q23,phys,raw33,recs,mixed,D,E,kap,T,H,saved,rebuilt,rel,labels,literal,pins])+10000
 assert owned<60000<150000 and raw+1000000<32<<20 and PRE<=2048
 proof={
  'all660_columns_rebuilt_by_actual_local_reduction':True,'nonzero_sparse_columns':25,
  'all26_actual_normalized_source_families_and1187_multiples':True,
  'all660_D_linearity_and21_actual_basis_images':True,
  'negative_RH_mixed_RB_beta2_carry_retained':True,
  'four_E7_operator_bindings_checked_without_pairing':True,
  'actual_source_order_and_completeness_proof':[
   'The accepted all-order ordered-prefix frame zB=A*z+b has detA0=-1/1024. Replacing its first7 generators by A^-1*zB is an actual unit change; beta2=A0^-1*b2 and beta3=A0^-1*(b3-A1*A0^-1*b2). Unknown A2*b2 starts at order4.',
   'IA gives exact zi*A=0. RB gives zi*B=-beta2_i*B-beta3_i*B+O(m4C). beta2 is proved plane; beta3 nonplane terms contain z and are cubic, hence in m4C because z*target is in actual m2C.',
   'RH_i=zi*H+delta_i*B+(bAi-zi)*H is actual. The accepted zero plane quadratic coefficients imply zi*H lies in m3C. Quadratic H terms containing z and cubic H terms therefore vanish modulo m4. RH delta2 mixed*B must instead be reduced through RB, producing exactly -delta2_mixed*beta2.',
   'Two z factors times B lie in m4 by RB; one z in any cubic target lies in m4. All plane H cubics lie in m4 by the actual quadratic Pk. None of these statements sets a missing coefficient to zero.',
   'The actual first-frame Pk have plane quadratic H plus the peered FULL P3_A/P3_B; the latter contains the negative beta carry. Their other cubic pieces and order4 terms are q-zero by the preceding order argument. Curve r3 coefficients are known units; their accepted ratios are 8/25 and96/275. Their remaining pure cubic terms belong to F.',
   'These local rules reduce EVERY raw degree<=3 column to a linear combination of the literal actual21 basis modulo m4C. Thus the formula is the existing quotient C4/F16 on generators, rather than a claim based only on26 relations annihilated in an abstract larger cover. Exact length21/HF3666 and the accepted literal basis ensure the reductions are actual and complete.'
  ],
  'D_only_and_genuine_cycle_restriction':'F=span_D(A,B,rA,rB) is primitive D-free4 but not R0-stable. q is D-linear only. For a genuine eta inC[ell], q(eta) is inkerE7 and rows1,2,5,6 give its four genuine coordinates; the accepted Jordan inverse may then be applied. Arbitrary vectors must not be called genuine cycles.',
  'formal_edge_frontier':'Form V_full times an ell-independent lift of the old cycle, prove ell0 zero and divide formally byell BEFORE finiteq. Highest remaining needed numerator monomials are A/B r3ell andr2ell2. V4sigma0,V3sigma1,V2sigma2,V1sigma3 remain UNKNOWN; no highest cancellation or edge value is asserted.',
  'V1_exact_single_field_label_history_verified':True,
  'edge_or_new_affine_responses_evaluated':False,
  'owned_numeric_including_all_labels_and10000_reserve':owned,
  'all_QQ_preoperation_bits':PRE,'unique_RAW_inputs':raw,
  'hard_caps':{'seconds':60,'memory':2<<30,'numeric':150000,'bits':2048,'RAW':32<<20}}
 ep=out/'independent-local-source-reductions-and-Dlinear-q7-checks.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
  'actual_D_only_projector7_source_and1187_multiples_PASS':True,
  'new_affine_or_connecting_values_computed':False,'seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'all_passed':True,'owned':owned,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
