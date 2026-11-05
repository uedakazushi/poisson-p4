#!/usr/bin/env python3
"""Sparse, source-bound D-linear C/m^4 -> T7 projector; no edge values."""
import hashlib,itertools,json,math,re,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63'
PATHS={
 'T7':'reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63/actual-filtered23-Dfree16-quotient7-and-genuine4-source.json',
 'T7_peer':'reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63_REVIEW/report.json',
 'E7':'reports/AMBIENT_B18_C_GENUINE_T7_CACHED_VALUES_20261105_v63/four-cached-coefficients-actual-E7-and-genuine-four-target-coordinates.json',
 'E7_peer':'reports/AMBIENT_B18_C_GENUINE_T7_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'R0':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63/genuine-four-R0-action-transport-and-cached-HB-operator-identities.json',
 'R0_peer':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63_REVIEW/report.json',
 'formal':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/full-HB6-by2-affine-source-operator-and-old27-leaves.json',
 'formal_peer':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63_REVIEW/report.json',
 'HB':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json',
 'HB_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json',
 'beta':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'beta_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'exact':'reports/AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63/report.json',
 'exact_peer':'reports/AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63_REVIEW/report.json',
}
Z=(0,)*9; CONST=('QQ',)
MAX_PRE=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.keys()))+sum(map(count,x.values()))
 return 0
def height(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def addq(a,b):
 global MAX_PRE
 h=2*max(height(a),height(b))+1;MAX_PRE=max(MAX_PRE,h);assert h<=2048
 return a+b
def mulq(a,b):
 global MAX_PRE
 h=height(a)+height(b);MAX_PRE=max(MAX_PRE,h);assert h<=2048
 return a*b
def affine_add(a,b):
 v=a.copy()
 for k,q in b.items():
  v[k]=addq(v.get(k,Q(0)),q)
  if not v[k]:del v[k]
 return v
def scale(a,q):return{k:mulq(v,q)for k,v in a.items()if v and q}
def const(q):return{CONST:Q(q)}if q else{}
def leaf(row,ix):return{('cache',row,ix):Q(1)}
def decode(a):return{tuple(k):Q(v)for k,v in a}
def encode(a):return[[list(k),str(v)]for k,v in sorted(a.items())]
def mon(**kw):
 v=[0]*9
 for k,p in kw.items():v[int(k[1:])if k.startswith('z')else {'r':7,'ell':8}[k]]=p
 return tuple(v)
def plusmon(a,b):return tuple(x+y for x,y in zip(a,b))
def mons(d):
 out=[]
 for degree in range(d+1):
  for ids in itertools.combinations_with_replacement(range(9),degree):
   v=[0]*9
   for i in ids:v[i]+=1
   out.append(tuple(v))
 return out
def selective(path,fields,receipts):
 b=path.read_bytes();s=b.decode();dec=json.JSONDecoder();out={}
 for key in fields:
  matches=list(re.finditer(r'"'+re.escape(key)+r'"\s*:',s));assert len(matches)==1,(path,key)
  start=matches[0].end()
  while s[start].isspace():start+=1
  value,end=dec.raw_decode(s,start);lo=len(s[:start].encode());hi=len(s[:end].encode())
  receipts.append({'path':str(path.relative_to(ROOT)),'field':key,'byte_slice':[lo,hi],
    'slice_sha256':hashlib.sha256(b[lo:hi]).hexdigest(),'numeric_leaves':count(value)})
  out[key]=value
 return out
def main():
 started=time.monotonic();signal.alarm(60)
 pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[];d={}
 fields={
  'T7':['literal_filtered23_basis_F16_then_t7','actual_D_linear_quotient_q7_by23'],
  'E7':['all4_actual_cached_coefficients','actual_E7','actual_kernel7_by4'],
  'R0':['basis_change_inverse4'],
  'formal':['P3_A','P3_B'],
  'HB':['constant_source_normalization_T','source_first7_h_rows'],
  'beta':['beta2_before_ell','beta3_before_ell','gA_plane3','gB_own_plane3'],
  'frame':['all12_verified_source_records'],
  'C2':['physical12_coefficients'],
  'raw33':['actual33_basis'],
 }
 # Only named top-level fields are decoded: large unused action/source tables
 # are never parsed, including their rational-string copies.
 for k,p in PATHS.items():
  d[k]=selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 assert d['exact']['actual_comparison_C_torsion_exponents']==[1,1,2,3]
 source_count=count(d)+count(receipts);assert source_count+40000<150000
 P=[[list(map(decode,row))for row in d['formal'][side]]for side in ['P3_A','P3_B']]
 a1,b1=Q(8,25),Q(96,275)
 kap=[[affine_add(scale(P[side][k][0],ab),scale(P[side][k][1],Q(-1)))for k in range(3)]for side,ab in enumerate([a1,b1])]
 beta=[];beta3=[]
 for order,key,dest in[(2,'beta2_before_ell',beta),(3,'beta3_before_ell',beta3)]:
  for row in d['beta'][key]:
   arr=[Q(0)]*(order+1)
   for exp,s in row:
    if any(exp[:7]):continue
    nr,ns=exp[7:];assert nr+ns==order
    for j in range(ns+1):arr[j]=addq(arr[j],mulq(Q(s),Q(math.comb(ns,j),5**j)))
   dest.append(arr)
 assert[i for i,v in enumerate(beta)if any(v)]==[0,4,5]
 raw=d['raw33']['actual33_basis'];phys={r['physical_row']:dict(r['delta'])for r in d['C2']['physical12_coefficients']}
 records=d['frame']['all12_verified_source_records'][:7];rows=[x['physical_row']for x in records]
 assert rows==[18,19,28,29,60,67,72]
 D=[];Gamma=[];E=[];literal_zero=[]
 for i,row in enumerate(rows):
  rec=records[i];assert rec['delta_first9_zero']and rec['delta_and_bA_pureplane_second_zero']
  assert list(map(Q,rec['bA_first_z7']))==[Q(i==j)for j in range(7)]
  di=[leaf(row,29),leaf(row,30)];ga=[{},{}]
  for z in[0,4,5]:
   for shift in[0,1]:
    ix=raw.index([z+1,1-shift,shift]);u=leaf(row,ix)
    for j in range(2-shift):ga[j+shift]=affine_add(ga[j+shift],scale(u,beta[z][j]))
  if i not in[0,4,5]:
   used={k for a in di+ga for k in a}
   assert all(Q(phys[row].get(k[2],'0'))==0 for k in used)
   literal_zero.extend([list(k)for k in sorted(used)]);di=[{},{}];ga=[{},{}]
  D.append(di);Gamma.append(ga);E.append([affine_add(x,scale(y,Q(-1)))for x,y in zip(di,ga)])
 T=[[Q(v)for v in r]for r in d['HB']['constant_source_normalization_T']]
 H=[[Q(v)for v in r]for r in d['HB']['source_first7_h_rows']]
 for k in range(3):
  for m in range(2):
   PA={};PB={}
   for j,i in enumerate([0,4,5]):
    PA=affine_add(PA,scale(D[i][m],T[j][k]))
    coeff=sum((mulq(T[p][k],H[p][i])for p in range(3)),Q(0))
    PB=affine_add(PB,scale(E[i][m],coeff))
   assert PA==P[0][k][m]and PB==P[1][k][m]
 # Sparse vector entries are affine OLD-cache operators, never paired.
 def vec(*pairs):return{i:a for i,a in pairs if a}
 def va(a,b):
  v={i:x.copy()for i,x in a.items()}
  for i,x in b.items():
   v[i]=affine_add(v.get(i,{}),x)
   if not v[i]:del v[i]
  return v
 def vs(a,q):return{i:scale(x,q)for i,x in a.items()if q}
 table={}
 def put(target,m,v):
  if v:assert(target,m)not in table;table[(target,m)]=v
 put(2,Z,vec((0,const(1))));put(2,mon(r=1),vec((1,const(1))));put(2,mon(ell=1),vec((2,const(1))))
 for i in range(7):
  zi=mon(**{'z'+str(i):1})
  put(1,zi,vec((4,const(-beta[i][0])),(6,const(addq(mulq(b1,beta3[i][0]),-beta3[i][1])))))
  put(2,zi,vec((6,affine_add(scale(E[i][0],b1),scale(E[i][1],Q(-1))))))
  put(1,plusmon(zi,mon(r=1)),vec((6,const(addq(mulq(b1,beta[i][0]),-beta[i][1])))))
  put(1,plusmon(zi,mon(ell=1)),vec((6,const(-beta[i][0]))))
 for side,ab in enumerate([a1,b1]):
  put(side,mon(r=2),vec((3+side,const(1))))
  put(side,mon(r=3),vec((5+side,const(-ab))))
  put(side,mon(r=2,ell=1),vec((5+side,const(1))))
 for k,m in enumerate([mon(r=2),mon(r=1,ell=1),mon(ell=2)]):put(2,m,vec((5,kap[0][k]),(6,kap[1][k])))
 def proj(poly):
  v={}
  for (target,m),a in poly.items():
   for i,b in table.get((target,m),{}).items():
    # Source coefficients a below are scalars or affine leaves, never both
    # nonlinear: only constants multiply affine projector entries.
    if set(a)<= {CONST}:term=scale(b,a.get(CONST,Q(0)))
    else:
     assert set(b)<={CONST};term=scale(a,b.get(CONST,Q(0)))
    v=va(v,vec((i,term)))
  return v
 def relation(terms):
  p={}
  for target,m,a in terms:
   key=(target,m);p[key]=affine_add(p.get(key,{}),a)
  return{k:v for k,v in p.items()if v}
 def multiply(poly,m):return{(t,plusmon(e,m)):a for(t,e),a in poly.items()if sum(plusmon(e,m))<=3}
 rels=[]
 for i in range(7):
  zi=mon(**{'z'+str(i):1});rels.append(('RA'+str(i),1,{(0,zi):const(1)}))
  rb=[(1,zi,const(1))]+[(1,mon(r=2-k,ell=k),const(q))for k,q in enumerate(beta[i])]
  rb +=[(1,mon(r=3-k,ell=k),const(q))for k,q in enumerate(beta3[i])]
  rels.append(('RB'+str(i),1,relation(rb)))
  rh=[(2,zi,const(1))]
  for z in[0,4,5]:
   for shift in[0,1]:
    ix=raw.index([z+1,1-shift,shift]);a=leaf(rows[i],ix)
    if i not in[0,4,5]:assert Q(phys[rows[i]].get(ix,'0'))==0;a={}
    rh.append((1,plusmon(mon(**{'z'+str(z):1}),mon(r=1-shift,ell=shift)),a))
  rh +=[(1,mon(r=3-m,ell=m),D[i][m])for m in range(2)]
  rels.append(('RH'+str(i),1,relation(rh)))
 for k,m in enumerate([mon(r=2),mon(r=1,ell=1),mon(ell=2)]):
  terms=[(2,m,const(1))]
  for side in range(2):terms +=[(side,mon(r=3-j,ell=j),P[side][k][j])for j in range(2)]
  rels.append(('P'+str(k),2,relation(terms)))
 for side,key in enumerate(['gA_plane3','gB_own_plane3']):
  # Use only the accepted leading curve RATIO, independent of multiplying
  # the actual curve by a scalar unit or adding z*RA/RB.
  ab=[a1,b1][side];rels.append(('curve'+str(side),3,{(side,mon(r=3)):const(1),(side,mon(r=2,ell=1)):const(ab)}))
 checks=[];digest=hashlib.sha256()
 for name,order,p in rels:
  n=0
  for m in mons(3-order):
   v=proj(multiply(p,m));assert not v,(name,m,v);n+=1
   digest.update(json.dumps([name,list(m)],separators=(',',':')).encode())
  checks.append([name,order,n])
 assert len(rels)==26 and sum(n for _,_,n in checks)==1187
 # D-linearity, without pretending q is R0-linear. Symbolic E7 is affine;
 # a scalar or a high-socle vector is sufficient on all sparse columns.
 def Eaction(v):
  w={}
  for i,a in v.items():
   if i in[0,3,4]:w=va(w,vec(({0:2,3:5,4:6}[i],a)))
   elif i in[1,2]:
    assert set(a)<={CONST}
    for side in range(2):w=va(w,vec((5+side,scale(kap[side][i],a.get(CONST,Q(0))))))
  return w
 allmons=mons(3);assert len(allmons)==220
 for t in range(3):
  for m in allmons:
   assert table.get((t,plusmon(m,mon(ell=1))),{})==Eaction(table.get((t,m),{}))
 basis=d['T7']['literal_filtered23_basis_F16_then_t7']
 for j,b in enumerate(basis):
  if b['m_degree']>=4:continue
  expected=vec((j-16,const(1)))if j>=16 else{}
  assert table.get((b['raw_target'],mon(r=b['r_power'],ell=b['ell_power'])),{})==expected
 assert all(not table.get((t,m),{})for t in range(3)for m in allmons if sum(m)==3 and(any(m[:7])or t==2))
 sparse=[{'raw_target':t,'monomial':list(m),'target7_affine':[[i,encode(a)]for i,a in sorted(v.items())]}for(t,m),v in sorted(table.items())]
 symbolic_slots=count(sparse)+count(checks)+count(P)+count(E)+count(beta)+count(beta3)+count(table)+count(d)+count(receipts)+20000
 assert symbolic_slots<100000
 data={
  'scope':'Actual D-linear q:C/m^4 C→T7. NO R0 quotient action, edge value, ell division in an Artin ring, or new physical/gauge/IFT coefficient.',
  'same_actual_R0_variables':['z0','z1','z2','z3','z4','z5','z6','r','ell'],
  'raw_targets':['A=tauA','B=tauB','H=tauA*tauB'],
  'coefficient_labels':'cache(row,index)=actual OLD C2 delta coefficient. Leaves are symbolic even when an already accepted value exists; no new affine response is paired.',
  'all220_monomials_degree_at_most3':[list(m)for m in allmons],
  'projector7_by660_sparse_columns':sparse,
  'unlisted_columns_zero_by_proof_not_unknown_filling':True,
  'beta2_plane_after_actual_s_equals_r_plus_ell_over5':[[str(q)for q in row]for row in beta],
  'beta3_plane_after_actual_s_equals_r_plus_ell_over5':[[str(q)for q in row]for row in beta3],
  'first7_source_rows':rows,'first7_effective_E3_first2_UNPAIRED':[[encode(a)for a in row]for row in E],
  'literal_existing_zero_leaf_occurrences':literal_zero,
  'kappa_k_A_B_UNPAIRED':[[encode(a)for a in row]for row in kap],
  'normalized26_source_relation_multiplier_checks':checks,
  'all1187_projected_source_multiples_zero':True,'source_multiplier_labels_hash':digest.hexdigest(),
  'all660_D_linearity_identities':True,'all21_literal_basis_images_checked':True,
  'normalization_proof':'RA_i=z_i A exactly. RB_i=(z_i+beta_i2+beta_i3+O4)B in the accepted regular7 unit frame. RH_i=z_iH+delta_i2 B+delta_i3 B+(bAi-zi)H. Nonplane beta3, delta3 have z and project0. For RH, delta2_zi*r/ell B must use the nonzero RB quadratic carries: E_i3=delta_i3−delta_i2_mixed*beta2. First7 bAi2_plane=delta_i2_plane=0; all omitted terms contain a z in degree3 or are plane-cubic H, and q kills them. Normalize Pk by the actual first-frame operations; their negative RB beta2 carry is exactly the peered P3_B. Thus Pk=planequad_k H+P3_A A+P3_B B+omitted q-zero terms. Curve gA/gB leading r³ coefficients are units, and q sees only their r³/r²ell ratio a1/b1. Higher terms vanish modulo m4.',
  'completeness_proof':'The accepted actual HF(C) through3 is(3,6,6,6); its literal21 basis is the degree≤3 part of filtered23. The normalized26 actual source columns and all permitted scalar multiples supply the checked q-zero identities. q agrees with the already proved quotient C4/F16 on all21 basis vectors, and both kill m4. Hence this is the actual q, not a new cover NF or a presumed R0-linear quotient.',
  'kernel_projection':'For genuine eta∈C[ell], q(eta)∈kerE7 and its four coordinates are rows[1,2,5,6]. Apply the accepted genuine basis_change_inverse4 for J2+k² coordinates. Do NOT apply those selectors to an arbitrary noncycle and call it genuine.',
  'edge_SOURCE_order_frontier':{
    'numerator':'V_full * ell-independent lift((I−nbar*mbar) E56×5), divided formally by ell BEFORE this finite q',
    'highest_needed_raw_eta':'A/B plane degree3 r³,r²ell; numerator4 r³ell,r²ell²',
    'lower_needed':'B zi degree1, B zi*r/zi*ell degree2, H zi degree1, H plane degree≤2',
    'unknowns_retained':'V4,sigma3 and all coherent V3*sigma1/V2*sigma2 carry; no highest response cancellation assumed',
    'no_edge_values_computed':True,
  },
  'selective_source_extraction_receipts':receipts,
  'whole_source_numeric_upper':symbolic_slots,'input_numeric_count':source_count,
  'maximum_QQ_coefficient_preoperation_bound':MAX_PRE,
  'hard_caps':{'owned_numeric':150000,'QQ_bits':2048,'seconds':60,'memory_bytes':2147483648,'RAW_bytes':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json'
 ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ACTUAL_SPARSE_DLINEAR_Q7_NO_CONNECTING_VALUES',
  'actual_q7_source_proved':True,'source_multiplier_checks':1187,'D_linearity_checks':660,
  'edge_values_or_new_affine_responses_computed':False,'input_sha256':pins,
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':symbolic_slots,'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()
