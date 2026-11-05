#!/usr/bin/env python3
"""Independent coefficient convolution for head30, safe63/trusted scopes."""
import itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha,mon
from ambient_B18_C_connecting_head23_source_review_20261105_v63 import add,sc,mul,sum_poly,const,var,dec,dot,sub
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_C_CONNECTING_HEAD30_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=True);assert not (out/'report.json').exists();pins={};rawbytes=0
 def pin(p,h=None,parse=True):
  nonlocal rawbytes
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:rawbytes+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text())if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_connecting_head30_source_20261105_v63.py',pr['code_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def src(s):return ROOT/next(p for p in pins if p.endswith(s))
 up=src('actual-fullV56-m0-and-typed-comparison.json');selected=field(up,'fixed_selected24_columns');iv=[[Q(q)for q in row]for row in field(up,'fixed24_minor_inverse')]
 n=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(selected):n[j]=iv[i]
 adj=src('actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json')
 assert field(adj,'actual_sigma0_56_by5_sparse')[3]==[[30,'1']]
 assert field(src('actual-five-head-formal-RegDivell-operator-and-source-binding.json'),'full_V56_reference_pairs_columns56_targets3')[30][0]==[0,1]
 L={i:Q(q)for i,q in field(adj,'fixed_L_Cn0_row24_sparse')};shifts=list(map(Q,field(adj,'actual_CB_shift_c7')));active=[i for i,q in enumerate(shifts)if q]
 assert set(L)=={8,12,13}and active==[0,4,5]
 qp=src('sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json');cols=field(qp,'projector7_by660_sparse_columns')
 plane=src('all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json')
 B2=[[Q(q)for q in row]for row in field(plane,'actual_plane2_B_three_monomials_rows3_columns56')]
 H2=[[Q(q)for q in row]for row in field(plane,'actual_plane2_H_three_monomials_rows3_columns56')]
 raw=field(src('actual-ell2-product33-and-raw45-actions.json'),'actual33_basis');physical=field(src('actual-physical12-and-C2-right-exact-quotient.json'),'physical12_coefficients')
 phys={r['physical_row']:r for r in physical};safe63=set(map(tuple,field(src('finite63-socle-safe-output-theorem-and-actual-lower-input-overlaps.json'),'basis63_monomials')))
 ga=[Q(0)]*4
 for (r,s),q in field(src('actual-special-fiber20-coefficients.json'),'gA_plane3'):
  for j in range(s+1):ga[j]+=Q(q)*Q(math.comb(s,j),5**j)
 assert ga==[Q(125,24),Q(5,3),Q(1,6),Q(2,375)]
 gA=[var('gA_plane3',i)for i in range(4)]
 # Bind the7 oldnegativeC0 entries to the EXACT same CA functions, rather
 # than silently dropping alltheir plane3 coefficients.
 chain=ROOT/'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json'
 pin(chain,'d57cdec83150cde5f2f342a83fa817c71b60dfa36d158c9e8636b04d96fef379',False)
 aliases=field(chain,'accepted_actual_regular_function_identities')
 expected=[(7148,16376),(7199,16378),(7454,16396),(7505,16398),(7688,17862),(7943,17876),(8045,17886)]
 assert [(x['odd_generator'],x['physical_odd_coefficient'])for x in aliases if x['kind']=='exact_QQ_linear_source_alias'][:7]==expected
 assert [phys[x['physical_row']]['roots']['bA_node']for x in physical[:7]]==[b for a,b in expected]
 owned=count([pr,D,n,iv,selected,L,shifts,cols,B2,H2,raw,physical,ga,safe63,aliases,pins])+75000
 assert owned<150000 and rawbytes+1000000<32<<20
 save(out/'positive-before-formal-head30-convolution.json',{'whole_SOURCE_upper':owned,'all_polynomial_values_UNPAIRED':True,'safe63_and_trusted_plane3_separate':True,'all24_first_and_second_source_slots_formal_until_contraction':True})
 def safe(e):
  m=tuple(i for i,q in enumerate(e)for _ in range(q));assert m in safe63;return var('SAFE_delta63',*m)
 def lower(side,row,e):
  ix=raw.index([next((i+1 for i,q in enumerate(e[:7])if q),0),e[7],e[8]])
  return var(side,row,ix)
 def dl(row):return var('delta_lambda12',row)
 def odd3(side,row,e):
  if side=='bA'and row==63:return gA[e[8]]
  if e[8]==2:return var('TRUSTED_plane3',side,row,1,2)
  if e[8]in[0,1]:return lower(side,row,e)
  return var('UNUSED_odd03',side,row)
 def delta3(row,e):
  if e[8]==2:return dl(row)
  if e[8]in[0,1]:return lower('delta',row,e)
  return var('UNUSED_delta03',row)
 def qaff(expr):return sum_poly(const(Q(q))if label==['QQ']else sc(var('delta',label[1],label[2]),Q(q))for label,q in expr)
 w=[{}for _ in range(24)]
 curve=safe(mon(r=4))
 for i in active:curve=add(curve,sc(safe(mon(r=2,z=i)),-shifts[i]))
 for i,j in itertools.combinations_with_replacement(active,2):
  e=[0]*9;e[i]+=1;e[j]+=1;curve=add(curve,sc(safe(e),shifts[i]*shifts[j]))
 w[15]=curve;w[23]=gA[0]
 for i in range(7):w[8+i]=lower('delta',63,mon(r=1,z=i));w[16+i]=lower('bA',63,mon(r=1,z=i))
 v=[sc(dot(row,w),-Q(1))for row in n]
 assert [[i,p]for i,p in enumerate(w)if p]==[[i,dec(p)]for i,p in D['sigma1_r_source24']]
 assert [[i,p]for i,p in enumerate(v)if p]==[[i,dec(p)]for i,p in D['sigma1_r56']]
 assert not v[55]and all(not v[j]for j in[29,30,31])
 def mixedB(i,ell,j):
  if j>=48:return {}
  row=physical[j//4]['physical_row'];part=j%4;e=mon(z=i,**({'ell':1}if ell else{'r':1}))
  return sc(lower('bB',row,e),-Q(1))if part==0 else(lower('delta',row,e)if part==2 else{})
 vz={}
 for i in active:
  wz=[var('UNKNOWN_wz',i,k)for k in range(24)]
  for k in active:
   e=[0]*9;e[k]+=1;e[i]+=1;wz[8+k]=safe(e)if k<=i else{}
  vz[i]=[sc(dot(row,wz),-Q(1))for row in n]
 U=[var('UNKNOWN_U',i)for i in range(24)]
 for i in active:
  hi=safe(mon(z=i,r=2));e=[0]*9;e[i]=2;hi=add(hi,sc(safe(e),-shifts[i]))
  for j in active:
   if j>i:e=[0]*9;e[i]+=1;e[j]+=1;hi=add(hi,sc(safe(e),-shifts[j]))
  U[8+i]=add(hi,sum_poly(mul(mixedB(i,False,j),x)for j,x in enumerate(v)if x))
 v2=[sc(dot(row,U),-Q(1))for row in n]
 def V(t,e,j):
  degree=sum(e)
  if degree<=1:return {}
  if degree==2:
   if not any(e[:7]):return const((B2 if t==1 else H2)[e[8]][j])if t in[1,2]else{}
   assert t==1;return mixedB(next(i for i,q in enumerate(e[:7])if q),bool(e[8]),j)
  assert degree==3 and not any(e[:7])and e[8]in[1,2]
  if j>=48:
   if t==0 and j<55:return sc(odd3('bA',physical[j-48]['physical_row'],e),-Q(1))
   return {}
  row=physical[j//4]['physical_row'];part=j%4
  if t==0:
   if part==0:return odd3('bA',row,e)
   if part==1:return delta3(row,e)
   return {}
  if t==1:
   if part==0:return sc(odd3('bB',row,e),-Q(1))
   if part==2:return delta3(row,e)
   return {}
  return var('UNUSED_H3',j,e[8])
 def native(t,e):
  degree=sum(e)
  if t==0:return {}
  if t==1:
   if degree==2:return lower('delta',63,e)if any(e[:7])else const(B2[e[8]][30])
   return safe(e)
  if degree==2:return lower('bA',63,e)if any(e[:7])else const(H2[e[8]][30])
  assert degree==3 and not any(e[:7]);return gA[e[8]]
 result={i:{}for i in[1,2,5,6]}
 for column in cols:
  selected=[(i,qaff(a))for i,a in column['target7_affine']if i in result]
  if not selected:continue
  e=list(column['monomial']);e[8]+=1;e=tuple(e);t0=column['raw_target'];p=native(t0,e)
  for exponent,vector in[(mon(r=1),v)]+[(mon(z=i),vz[i])for i in active]+[(mon(r=2),v2)]:
   rem=sub(e,exponent)
   if rem is None or sum(rem)<=1:continue
   p=add(p,sum_poly(mul(V(t0,rem,j),x)for j,x in enumerate(vector)if x))
  for i,q in selected:result[i]=add(result[i],mul(q,p))
 assert not result[1]and not result[2]
 polys=[result[5],result[6],sc(result[5],Q(5,3)),add(result[6],sc(result[5],Q(-121,48)))]
 expected_polys=list(map(dec,D['four_genuine_and_canonical_high_polynomials_UNPAIRED']))
 diffs=[add(a,sc(b,-Q(1)))for a,b in zip(polys,expected_polys)]
 zero_label=('delta',63,21)
 assert raw[21]==[2,0,1]and 21 not in dict(phys[63]['delta'])
 beta3=field(qp,'beta3_plane_after_actual_s_equals_r_plus_ell_over5')
 gamma=Q(96,275)*Q(beta3[1][0])-Q(beta3[1][1]);assert gamma==Q(-827,63360)
 assert diffs==[{},sc(var(*zero_label),gamma),{},sc(var(*zero_label),gamma)]
 def normalize_zero_leaf(p):return{m:q for m,q in p.items()if zero_label not in m}
 polys=list(map(normalize_zero_leaf,polys));expected_polys=list(map(normalize_zero_leaf,expected_polys))
 assert polys==expected_polys
 labels={label for p in polys for m in p for label in m}
 assert not any(label[0].startswith(('UNKNOWN','UNUSED'))for label in labels)
 assert not any(label[:2]==('TRUSTED_plane3','bA')for label in labels)
 assert {x[2]for x in labels if x[:2]==('TRUSTED_plane3','bB')}=={18,19,28,29,60,63,67,72}
 for k,q,locator in D['known_cached_inputs']:
  if k[0]in['delta','bA','bB']:assert Q(q)==Q(dict(phys[k[1]][k[0]]).get(k[2],'0'))
  elif k[0]=='gA_plane3':assert Q(q)==ga[k[1]]
  elif k[0]=='delta_lambda12':
   cache={r['physical_row']:Q(r['delta_r_ell2'])for r in field(src('phase1-six-root-existing-cache-readbacks.json'),'all3_lambda_readbacks')}
   assert Q(q)==cache.get(k[1],Q(0))
  else:raise AssertionError(k)
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'independent_generic_q7_convolution_matches_four_operators_MOD_one_explicit_proved_cached_zero':True,'head30_actual_sigma0_e30_native_V_0_delta63_bA63':True,
  'full_adapted_Bcurve4_with_all_six_diagonal_and_offdiagonal_zz_terms_retained':True,'original_z_Hsector_and_Hadamard_second_prefix_corrected':True,
  'all_formal_sigma1z_and_sigma2_components_retained_until_full_projection':True,'full_minus_Lm1sigma1_and_all_Bmixed_qH_qB_carry_kept':True,
  'seven_actual_CA_bA_and_negative_C0_function_aliases_verified':expected,
  'all_first7_bA_plane3_callbacks_cancel_algebraically_NOT_set_zero':True,
  'inactive_beta3_carry_clarification':{'unconstrained_polynomial_difference_B_and_k2':'-827/63360 * delta63_z1_ell','actual_delta63_raw33_slot':21,'raw33_slot21':[2,0,1],'canonical_sparse_delta63_slot21_absent_exactzero':True,'all_other_formal_coefficients_match_before_normalization':True,'not_unrestricted_free_coefficient_equality':True},
  'remaining_eight_bB_plane3_callbacks':[18,19,28,29,60,63,67,72],
  'safe_delta_callback_precision':'Every SAFE_delta63 label lies in the accepted63monomial ring. Other new callback labels are ONLY trusted pureplane3 bB functions, not their fullB63 images.',
  'no_root_composition_or_head_response_value_or_new_source_solve_evaluated':True,'no_original_H0_or_nil_claim':True,'whole_SOURCE_upper':owned,'RAW':rawbytes}
 ep=out/'independent-full-coherent-head30-q-convolution-and-literal-seven-CA-aliases.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'full_head30_source_and_callback_scopes_PASS':True,'actual_head_response_evaluated':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':main()
