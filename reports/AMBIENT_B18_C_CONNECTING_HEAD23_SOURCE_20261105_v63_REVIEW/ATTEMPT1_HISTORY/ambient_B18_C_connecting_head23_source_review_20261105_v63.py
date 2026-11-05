#!/usr/bin/env python3
"""Independent q7 coefficient convolution for coherent original head23."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha,mon
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRE=0
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def qa(a,b,mul=False):
 global PRE
 if a and b:
  na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
  h=max(na+nb,da+db) if mul else max(max(na+db,nb+da)+1,da+db)
  PRE=max(PRE,h);assert h<=2048
 return a*b if mul else a+b
def add(a,b):
 d=dict(a)
 for k,x in b.items():
  d[k]=qa(d.get(k,Q(0)),x)
  if not d[k]:del d[k]
 return d
def sc(p,q):return {k:qa(x,q,True) for k,x in p.items() if x and q}
def mul(p,q):
 d={}
 for a,x in p.items():
  for b,y in q.items():
   k=tuple(sorted(a+b));assert len(k)<=2
   d[k]=qa(d.get(k,Q(0)),qa(x,y,True))
 return {k:x for k,x in d.items() if x}
def sum_poly(xs):
 p={}
 for x in xs:p=add(p,x)
 return p
def const(x):return {():Q(x)} if x else {}
def var(*key):return {(tuple(key),):Q(1)}
def dec(v):return {tuple(tuple(label) for label in labels):Q(q) for labels,q in v}
def dot(row,v):return sum_poly(sc(x,q) for q,x in zip(row,v) if q and x)
def sub(e,f):return tuple(x-y for x,y in zip(e,f)) if all(x>=y for x,y in zip(e,f)) else None
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_connecting_head23_source_20261105_v63.py',pr['code_sha256'],False)
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def find(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 up=find('actual-fullV56-m0-and-typed-comparison.json');sel=field(up,'fixed_selected24_columns')
 iv=[[Q(x) for x in row] for row in field(up,'fixed24_minor_inverse')]
 n=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(sel):n[j]=iv[i]
 adj=find('actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json')
 assert field(adj,'actual_sigma0_56_by5_sparse')[1]==[[23,'1']]
 L=dict((i,Q(x)) for i,x in field(adj,'fixed_L_Cn0_row24_sparse'));assert set(L)=={8,12,13}
 qp=find('sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json')
 qcols=field(qp,'projector7_by660_sparse_columns');beta=[[Q(x) for x in row] for row in field(qp,'beta2_plane_after_actual_s_equals_r_plus_ell_over5')]
 plane=find('all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json')
 B2=[[Q(x) for x in row] for row in field(plane,'actual_plane2_B_three_monomials_rows3_columns56')]
 H2=[[Q(x) for x in row] for row in field(plane,'actual_plane2_H_three_monomials_rows3_columns56')]
 cb=find('actual-Cbar-minimal-resolution-and-three-targets.json');shift=[Q(x) for x in field(cb,'CB_shift_coefficients7')]
 assert [i for i,c in enumerate(shift) if c]==[0,4,5]
 raw=field(find('actual-ell2-product33-and-raw45-actions.json'),'actual33_basis')
 physical=field(find('actual-physical12-and-C2-right-exact-quotient.json'),'physical12_coefficients')
 phys={x['physical_row']:x for x in physical};pindex={r:i for i,r in enumerate(phys)}
 q20=find('phase1-common105-and-114-source-equations.json');base=field(q20,'base27_monomials');nus=field(q20,'nu_A_B27')
 lp=find('phase1-six-root-existing-cache-readbacks.json');mp=find('phase1-six-cached-pivot-pure-ell4-values.json')
 lambdas={x['physical_row']:x['delta_r_ell2'] for x in field(lp,'all3_lambda_readbacks')}
 mus={x['physical_row']:x['delta_ell3'] for x in field(mp,'all3_pivot_mu03_readbacks')}
 # Explicit accepted lambda-free zeros, including the omitted v25 source.
 azstem='AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63';azr=pin(ROOT/'reports'/azstem/'report.json')
 azpeer=pin(ROOT/'reports'/(azstem+'_REVIEW')/'report.json');assert azpeer['all_passed']
 azname='eight-cached-even-root-values-and-four-complete-alpha-readbacks.json'
 az=pin(ROOT/'reports'/azstem/azname,azr['evidence_sha256'][azname])
 zeros={x['physical_row'] for x in az['all4_actual_alpha_delta_B6'] if x['delta_r_ell2']=='0'}
 assert zeros=={19,28,29,72}
 # Source-positive ownership gate, before constructing output polynomials.
 owned=count([pr,D,n,iv,sel,qcols,beta,B2,H2,shift,raw,physical,base,nus,lambdas,mus,azr,azpeer,az,pins])+70000
 assert owned<150000 and RAW+1000000<32<<20
 save(out/'SOURCE-positive-before-symbolic-convolution.json',{'no_output_pairing':True,'all_polynomials_degree_at_most2':True,'symbolic_source_reserve':70000,'whole_numeric_upper':owned,'hard_bits':2048,'RAW_plus_output_reserve':RAW+1000000})
 def delta(row,e):
  if e==mon(r=1,ell=2):return {} if row in zeros else var('lambda12',row)
  if e==mon(ell=3):return var('mu03',row)
  if e[8]>=2:return {} # only lower planeell² is requested here: actualzero.
  ix=raw.index([next((i+1 for i,x in enumerate(e[:7]) if x),0),e[7],e[8]])
  return var('cache',row,ix)
 def nu(f,slots):return var('nu',f,base.index(sorted(slots)))
 Az={i:nu(0,[i]) for i in [0,4,5]};Arr=nu(0,[7,7]);Brr=nu(1,[7,7])
 def qaff(expr):return sum_poly(sc(var(*label),Q(q)) if label!=['QQ'] else const(Q(q)) for label,q in expr)
 # Actual m1r*e23, including all three distinct curve rows.
 w=[{} for _ in range(24)];bb67=[-H2[k][21] for k in range(3)]
 w[7]=sc(Brr,-bb67[0]);w[15]=sc(add(Arr,sc(sum_poly(sc(Az[i],shift[i]) for i in Az),-Q(1))),-shift[5])
 for i in range(7):w[16+i]=delta(67,mon(r=1,z=i))
 w[23]=delta(67,mon(r=3));v=[sc(dot(row,w),-Q(1)) for row in n]
 assert [[i,p] for i,p in enumerate(w) if p]==[[i,dec(p)] for i,p in D['sigma1_r_source24']]
 assert [[i,p] for i,p in enumerate(v) if p]==[[i,dec(p)] for i,p in D['actual_sigma1_r56_as_cache_polynomials']]
 assert all(not v[i] for i in [0,4,8,12,16,20,24])
 def Bmixed(i,ell,j):
  if j>=48:return {}
  row=physical[j//4]['physical_row'];part=j%4;ix=raw.index([i+1,int(not ell),int(ell)])
  if part==0:return sc(var('bB',row,ix),-Q(1))
  if part==2:return var('cache',row,ix)
  return {}
 # Formal sigma1_zi outside the proven three-row adjoint remains unknown.
 vz={};rzproof=[]
 for i in [0,4,5]:
  aa=Q(96,275)*beta[i][0]-beta[i][1];bb=-beta[i][0]
  row=[aa*B2[1][j]+bb*B2[2][j] for j in range(56)]
  rz=[sum((row[j]*n[j][k] for j in range(56)),Q(0)) for k in range(24)]
  assert all(not q for k,q in enumerate(rz) if k not in [8,12,13])
  wz=[var('UNKNOWN_wz',i,k) for k in range(24)]
  wz[8]={};wz[12]={};wz[13]=Az[i]
  vz[i]=[sc(dot(row,wz),-Q(1)) for row in n]
  rzproof.append([i,[[k,str(q)] for k,q in enumerate(rz) if q]])
 assert rzproof==D['sigma1_z_complete_leftadjoints_Rz_n0']
 # Formal sigma2 preserves ALL outside rows until q contraction kills them.
 U=[var('UNKNOWN_U',k) for k in range(24)]
 for i in [0,4,5]:
  m2=sc(Az[i],-shift[5])
  if i==5:m2=add(m2,Arr)
  U[8+i]=add(m2,sum_poly(mul(Bmixed(i,False,j),v[j]) for j in range(56) if v[j]))
 v2=[sc(dot(row,U),-Q(1)) for row in n]
 # Generic actual V coefficient accessor, restricted ONLY to q-required
 # slots. All other sigma1_z/sigma2 dependencies are formal, not zero-filled.
 def V(t,e,j):
  deg=sum(e)
  if deg<=1:return {} # every required slot here has positiveell / plane.
  if deg==2:
   if not any(e[:7]):
    ix=e[8]
    return const((B2 if t==1 else H2)[ix][j]) if t in [1,2] else {}
   assert t==1 and sum(e[:7])==1 and e[7]+e[8]==1
   return Bmixed(next(i for i,x in enumerate(e[:7]) if x),bool(e[8]),j)
  assert deg==3 and not any(e[:7]) and e[8] in [1,2]
  if j>=48:return {} # v55=0; negativez pureplane=0.
  row=physical[j//4]['physical_row'];part=j%4
  if t==0:
   if part==1:return delta(row,e)
   if part==0:
    if row in [r['physical_row'] for r in physical[:7]]:return {}
    if row==63:return const(Q(5,3) if e[8]==1 else Q(1,6))
    return var('UNUSED_graph_A3',row,e[8])
   return {}
  assert t==1
  if part==0:
   assert row==63 # other bB3 columns have sigma1r=0 before accessor.
   return sc(var('bB',63,30) if e[8]==1 else var('NEW_bB63_lambda12'),-Q(1))
  if part==2:return delta(row,e)
  return {}
 def native(t,M):
  deg=sum(M)
  if t==0:
   assert deg==4 and not any(M[:7]);r,l=M[7:];ans={}
   for a in range(3):
    b=l-a
    if 0<=b<=2:ans=add(ans,sc(nu(1,[7]*(2-a)+[8]*a),-bb67[b]))
   return ans
  if t==1:
   if M[5]!=1 or sum(M[:7])!=1:return {}
   rem=list(M);rem[5]-=1
   if sum(rem)!=2:return {} # nu pureplanefirst=0.
   return nu(0,[7]*rem[7]+[8]*rem[8])
  # Onlydelta67 plane2/3 andziell enters the selected q columns.
  if deg==2 and not any(M[:7]):return {}
  return delta(67,M)
 # Coherent numerator convolution sigma0+sigma1_r+sum sigma1_zi+sigma2_rr.
 qresult={i:{} for i in [1,2,5,6]}
 for column in qcols:
  selected=[(i,qaff(a)) for i,a in column['target7_affine'] if i in qresult]
  if not selected:continue
  M=list(column['monomial']);M[8]+=1;M=tuple(M);t=column['raw_target']
  value=native(t,M)
  for sm,vector in [(mon(r=1),v)]+[(mon(z=i),vz[i]) for i in [0,4,5]]+[(mon(r=2),v2)]:
   rem=sub(M,sm)
   if rem is None:continue
   if sum(rem)<=1:continue # positiveell V1 and V0 areexactzero.
   value=add(value,sum_poly(mul(V(t,rem,j),x) for j,x in enumerate(vector) if x))
  for i,weight in selected:qresult[i]=add(qresult[i],mul(weight,value))
 assert not qresult[1] and not qresult[2]
 highs=[dec(p) for p in D['head1_genuine_high_A_B_UNPAIRED']]
 assert [qresult[5],qresult[6]]==highs
 canonical=[sc(qresult[5],Q(5,3)),add(qresult[6],sc(qresult[5],Q(-121,48)))]
 assert canonical==[dec(p) for p in D['head1_canonical_socle_coordinates_rg_k2_UNPAIRED']]
 labels={label for p in qresult.values() for monomial in p for label in monomial}
 assert not any(label[0].startswith('UNKNOWN') or label[0].startswith('UNUSED') for label in labels)
 # Two lower coordinates onALLfiveheads, kept as actualunpairedcachelabels.
 c0=field(adj,'actual_sigma0_56_by5_sparse');low=[]
 for h,vec in enumerate(c0):
  ds=dict((i,Q(q)) for i,q in vec);row=[]
  for order in [1,2]:row.append(sum_poly(sc(var('bB2',r,order),-ds.get(j,Q(0))) for r,j in [(18,1),(60,17),(67,21)]))
  low.append(row)
 assert low==[[dec(p) for p in row] for row in D['five_head_low_rows_rH_ellH_UNPAIRED']]
 # Check every cached leaf's ORIGINAL locator without evaluatingoutputs.
 for label,value,locator in D['all_input_values_already_cached']:
  key=tuple(label)
  if key[0] in ['cache','bB']:expected=dict(phys[key[1]]['delta' if key[0]=='cache' else 'bB']).get(key[2],'0')
  elif key[0]=='nu':expected=nus[key[1]][key[2]]
  elif key[0]=='lambda12':expected=lambdas[key[1]]
  elif key[0]=='mu03':expected=mus[key[1]]
  else:assert key[0]=='bB2';expected=-H2[key[2]][4*pindex[key[1]]+1]
  assert Q(value)==Q(expected)
 actual_owned=count([pr,D,n,iv,qcols,beta,B2,H2,shift,raw,physical,base,nus,lambdas,mus,az,labels,U,v,w,v2,vz,qresult,canonical,low,pins])+20000
 assert actual_owned<=owned and all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'full_q_required_coefficient_convolution_rebuilt':True,'two_highs_and_two_canonical_polynomials_match':True,
  'all14_high_canonical_and_low_fivehead_polynomials_checked':True,
  'all_formal_wzi_outside_B_rows_and_Uoutside_L_cancel_before_values':True,
  'sigma1_source_proof':[
   'Headfree56row23 isphysical67*tauAB andsigma0=e23. A source is-nuB*bB67: itsunshiftedr3quotient firstcoefficient=-nuBrr*bB67rr. Bsource isnuA*z5 ontheBcurve z=-c*r2, yielding-c5*(nuArr-sum ci*nuAzi). HsectorisUNshifted(z7,r2),soitscurvefirst isdelta67r3.',
   'Remainingm1r entries areH zHadamards delta67_zi*r. Allfirst7Funitcolumn components ofsigma1r vanish exactly. No primitivefourthcoefficient enters THIShead. This uses literal bA67=z5 and distinctB/Hsectorcurves.',
   'For sigma1_zi, retainformalall24sourcecomponents. The FULLqhigh leftadjoint kills everyrow outsideBz0/4/5; theBsource isnuA*z5, soits survivingcomponent is onlynuA_zi inrow13. Allotherformal unknowns cancel before any delta2zz read.',
   'The fullsigma2 correction is-n0*(m2sigma0+m1sigma1). OriginalzHadamard gives m2 rows8/12=-c5nuAz0/4, row13=nuArr-c5nuAz5. m1r onALL56sourcecolumns is the actualBmixedzi*r matrix; bothterms survive in the independentconvolution.'
  ],
  'source_carry_and_ancestry':'AllqHquadratic/zH andBmixedcoefficients, native-nuB2*bB67second carry, completeV3sigma1 andsigma2 terms were included by q coefficient extraction. Free-row lambda12zeros19/28/29/72 were explicitly pinned to acceptedALPHA4 values/peer, so omittedv25terms are provedzero, not missing-input defaults.',
  'remaining_new_label':['NEW_bB63_lambda12'],'original_root':38648,
  'no_polynomial_or_missing_root_values_paired':True,
  'kernel_type_conditional_proof':'M7=J2²+k³ surjectsG4=J2+k². rankrM=2,rankrG=1. Ifsimpletarget projections ofh0,h1 areindependent, imageSocMcontainsrG plusbothsimpletargets, hencehasdimension>=3. Asemisimplekernel3 wouldinsteadgiveimageSocMdimension2. Thereforekernelhasrankr1 andtypeJ2+k. Thetwo-targetminor islow_ellH(h0)*k2(h1) becausehead1low0; it is UNPAIRED.',
  'no_original_H0_or_nil_inference':True,'owned_upper':owned,'actual_preoperation_bits':PRE}
 ep=out/'independent-full-q-convolution-and-one-root-head23-SOURCE.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'coherent_head23_full_SOURCE_PASS':True,'cached_root_or_edge_values_evaluated':False,'seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
