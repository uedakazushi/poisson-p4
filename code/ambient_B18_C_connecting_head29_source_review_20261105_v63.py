#!/usr/bin/env python3
"""Independent full q7 convolution and all297 raw33-plane9 intertwiners."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha,mon
from ambient_B18_C_connecting_head23_source_review_20261105_v63 import add,sc,mul,sum_poly,const,var,dec,dot,sub
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63';AV='AMBIENT_B18_C_CONNECTING_HEAD29_CACHE_AVAILABILITY_20261105_v63';ID=PRI+'_REVIEW'
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
  return json.loads(p.read_text()) if parse else None
 def packet(st,pinh):
  r=pin(ROOT/'reports'/st/'report.json',pinh)
  for p,h in r['input_sha256'].items():
   pin(ROOT/p,h,False)
   if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
  en=next(iter(r['evidence_sha256']));v=pin(ROOT/'reports'/st/en,r['evidence_sha256'][en]);return r,v
 pr,D=packet(PRI,'f18d6c07e1ea37a460840ea9ae45581354d17d099d8ed09278ace01d4feb0dc8')
 ar,AD=packet(AV,'fc451d7f99198ef962d565221ec707255f42d4a477f581366075f9a0cc54e7ef')
 pin(ROOT/'code/ambient_B18_C_connecting_head29_source_20261105_v63.py',pr['code_sha256'],False)
 pin(ROOT/'code/ambient_B18_C_connecting_head29_cache_availability_source_v2_20261105_v63.py',ar['code_sha256'],False)
 def src(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 up=src('actual-fullV56-m0-and-typed-comparison.json');sel=field(up,'fixed_selected24_columns');iv=[[Q(x)for x in a]for a in field(up,'fixed24_minor_inverse')]
 n=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(sel):n[j]=iv[i]
 adj=src('actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json')
 assert field(adj,'actual_sigma0_56_by5_sparse')[2]==[[29,'1']]
 five=src('actual-five-head-formal-RegDivell-operator-and-source-binding.json');assert field(five,'full_V56_reference_pairs_columns56_targets3')[29][1]==[0,1]
 L=dict((i,Q(q))for i,q in field(adj,'fixed_L_Cn0_row24_sparse'));assert set(L)=={8,12,13}
 qp=src('sparse-actual-Dlinear-projector7-source-and-five-head-order-frontier.json');cols=field(qp,'projector7_by660_sparse_columns')
 beta=[[Q(q)for q in row]for row in field(qp,'beta2_plane_after_actual_s_equals_r_plus_ell_over5')]
 plane=src('all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json')
 B2=[[Q(q)for q in row]for row in field(plane,'actual_plane2_B_three_monomials_rows3_columns56')]
 H2=[[Q(q)for q in row]for row in field(plane,'actual_plane2_H_three_monomials_rows3_columns56')]
 assert all(row[29]==0 for row in B2+H2)
 rp=src('actual-ell2-product33-and-raw45-actions.json');basis=field(rp,'actual33_basis')
 physical=field(src('actual-physical12-and-C2-right-exact-quotient.json'),'physical12_coefficients');phys={r['physical_row']:r for r in physical}
 zeros={19,28,29,72};prior=pin(ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63_REVIEW/report.json','25cb1e164eb54b85f6835a51781df9cfa0f5adae5279193e22c1a48e34d2b281');assert prior['all_passed']
 # Positive gate includes formal-polynomial construction, not their values.
 owned=count([pr,D,ar,AD,n,iv,cols,beta,B2,H2,basis,physical,prior,pins])+70000
 assert owned<150000 and rawbytes+1000000<32<<20
 save(out/'positive-before-unpaired-head29-convolution.json',{'whole_SOURCE_upper':owned,'formal_degree_at_most2':True,'no_high_polynomial_pairing':True,'all_unknown_second_source_components_retained_until_projection':True})
 def delta(row,e):
  if sum(e)==3 and not any(e[:7]) and e[8]==2:return {} if row in zeros else var('lambda12',row)
  if e[8]>=2:return {}
  ix=basis.index([next((i+1 for i,x in enumerate(e[:7])if x),0),e[7],e[8]])
  if sum(e)==3 and not any(e[:7]) and row in zeros:assert Q(dict(phys[row]['delta']).get(ix,'0'))==0;return {}
  return var('cache',row,ix)
 def bb(ix):return var('bB',63,ix)
 def qaff(expr):return sum_poly(sc(var(*label),Q(q))if label!=['QQ']else const(Q(q))for label,q in expr)
 w=[{}for _ in range(24)];w[7]=var('NEW_delta63_plane4','r4');w[23]=sc(bb(29),-Q(1))
 for i in range(7):
  ix=basis.index([i+1,1,0]);w[i]=var('cache',63,ix);w[16+i]=sc(bb(ix),-Q(1))
 v=[sc(dot(row,w),-Q(1))for row in n]
 assert [[i,p]for i,p in enumerate(w)if p]==[[i,dec(p)]for i,p in D['sigma1_r_w24']]
 assert [[i,p]for i,p in enumerate(v)if p]==[[i,dec(p)]for i,p in D['sigma1_r_v56']]
 assert all(not v[i]for i in[0,4,8,12,16,20,24])
 def bmixed(i,ell,j):
  if j>=48:return {}
  row=physical[j//4]['physical_row'];part=j%4;ix=basis.index([i+1,int(not ell),int(ell)])
  return sc(var('bB',row,ix),-Q(1))if part==0 else(var('cache',row,ix)if part==2 else{})
 # All24 unprojected first-z and second-rr source slots stay independent.
 vz={}
 for i in[0,4,5]:
  wz=[var('UNKNOWN_wz',i,k)for k in range(24)]
  for k in[8,12,13]:wz[k]={}
  vz[i]=[sc(dot(row,wz),-Q(1))for row in n]
 U=[var('UNKNOWN_U',k)for k in range(24)]
 for i in[0,4,5]:U[8+i]=sum_poly(mul(bmixed(i,False,j),v[j])for j in range(56)if v[j])
 v2=[sc(dot(row,U),-Q(1))for row in n]
 def V(t,e,j):
  degree=sum(e)
  if degree<=1:return {}
  if degree==2:
   if not any(e[:7]):return const((B2 if t==1 else H2)[e[8]][j])if t in[1,2]else{}
   assert t==1;return bmixed(next(i for i,x in enumerate(e[:7])if x),bool(e[8]),j)
  assert degree==3 and not any(e[:7]) and e[8]in[1,2]
  if j>=48:return {}
  row=physical[j//4]['physical_row'];part=j%4
  if t==0:
   if part==1:return delta(row,e)
   if part==0:
    if row in[r['physical_row']for r in physical[:7]]:return {}
    if row==63:return const(Q(5,3)if e[8]==1 else Q(1,6))
    return var('UNUSED_A3',row,e[8])
   return {}
  assert t==1
  if part==0:
   assert row==63;return sc(bb(30)if e[8]==1 else var('bB63_lambda12'),-Q(1))
  return delta(row,e)if part==2 else{}
 def native(t,e):
  degree=sum(e)
  if t==0:
   assert degree==4 and not any(e[:7])
   return var('NEW_delta63_plane4',{0:'r4',1:'r3ell',2:'r2ell2'}[e[8]])
  if t==1:return {}
  if degree==2:
   if any(e[:7]):return sc(bb(basis.index([next(i+1 for i,x in enumerate(e[:7])if x),e[7],e[8]])),-Q(1))
   return const(H2[e[8]][29])
  assert degree==3
  return sc(bb(30)if e[8]==1 else(var('bB63_lambda12')if e[8]==2 else var('NEW_CACHED_bB63_ell3')),-Q(1))
 result={i:{}for i in[1,2,5,6]}
 for col in cols:
  selected=[(i,qaff(a))for i,a in col['target7_affine']if i in result]
  if not selected:continue
  M=list(col['monomial']);M[8]+=1;M=tuple(M);target=col['raw_target'];value=native(target,M)
  for exponent,vec in[(mon(r=1),v)]+[(mon(z=i),vz[i])for i in[0,4,5]]+[(mon(r=2),v2)]:
   rem=sub(M,exponent)
   if rem is None or sum(rem)<=1:continue
   value=add(value,sum_poly(mul(V(target,rem,j),x)for j,x in enumerate(vec)if x))
  for i,q in selected:result[i]=add(result[i],mul(q,value))
 assert not result[1]and not result[2]
 polynomials=[result[5],result[6],sc(result[5],Q(5,3)),add(result[6],sc(result[5],Q(-121,48)))]
 assert polynomials==list(map(dec,D['four_genuine_and_canonical_high_polynomials_UNPAIRED']))
 labels={x for p in polynomials for m in p for x in m};assert not any(x[0].startswith(('UNKNOWN','UNUSED'))for x in labels)
 known={tuple(k):(Q(q),loc)for k,q,loc in D['known_cached_inputs']}
 for k,(q,loc)in known.items():
  if k[0]in['cache','bB']:expect=dict(phys[k[1]]['delta'if k[0]=='cache'else'bB']).get(k[2],'0')
  elif k[0]=='lambda12':expect=next(r['delta_r_ell2']for r in field(src('phase1-six-root-existing-cache-readbacks.json'),'all3_lambda_readbacks')if r['physical_row']==k[1])
  else:assert k==('bB63_lambda12',);expect=field(src('phase1-single-root-existing-cache-readback.json'),'requested_bB63_r_ell2')
  assert q==Q(expect)
 # Separate completeness of the existing raw33 projection: all9 actions on
 # all33 basis classes, not a plain-monomial-NF assumption.
 actions=field(rp,'coordinate_actions9_sparse_columns');plane9=[(0,0),(1,0),(0,1),(2,0),(1,1),(3,0),(2,1),(4,0),(3,1)];pi={e:i for i,e in enumerate(plane9)}
 P=[({pi[(r,e)]:Q(1)}if z==0 else{})for z,r,e in basis]
 digest=hashlib.sha256()
 for i in range(9):
  for j in range(33):
   a={}
   for k,q in actions[i][j]:
    for v,c in P[k].items():a[v]=a.get(v,Q(0))+Q(q)*c
   a={k:q for k,q in a.items()if q};b={}
   if i>=7:
    for k,q in P[j].items():
     e=tuple(plane9[k][l]+int(l==i-7)for l in range(2))
     if e in pi:b[pi[e]]=b.get(pi[e],Q(0))+q
   assert a==b
   digest.update(json.dumps([i,j,sorted(a.items())],default=str,separators=(',',':')).encode())
 assert digest.hexdigest()==AD['all297_intertwiner_hash']
 assert [[[k,str(v)]for k,v in p.items()]for p in P]==AD['actual_projection9_by33_sparse_columns']
 assert basis[31:]==[[0,4,0],[0,3,1]]
 values=[Q(dict(phys[63]['delta'])[i])for i in[31,32]];assert values==[Q(-1236125,1728),Q(-91325,288)]
 replacement={('NEW_delta63_plane4','r4'):('CACHED_delta63_plane4','r4'),('NEW_delta63_plane4','r3ell'):('CACHED_delta63_plane4','r3ell')}
 refined=[]
 for p in polynomials:
  q={}
  for m,c in p.items():
   k=tuple(sorted(replacement.get(x,x)for x in m));q[k]=q.get(k,Q(0))+c
  refined.append({k:v for k,v in q.items()if v})
 assert refined==list(map(dec,AD['refined_four_head29_polynomials_UNPAIRED']))
 missing={x for p in refined for m in p for x in m if x[0].startswith('NEW_')}
 assert missing=={('NEW_delta63_plane4','r2ell2'),('NEW_CACHED_bB63_ell3',)}
 assert [p.get((('NEW_delta63_plane4','r2ell2'),),Q(0))for p in refined]==[1,0,Q(5,3),Q(-121,48)]
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'full_independent_q7_coefficient_convolution_matches_all_four_unpaired_high_operators':True,'head29_sigma0_e29_physical63_tauA_full_B_identically_zero_verified':True,
  'all_unknown_sigma1_z_and_sigma2_slots_project_out_before_values':True,'full_m1r_sigma1r_term_retained':True,
  'low_two_head29_coordinates_source_zero':True,'all_cached_literal_ancestry_verified_without_high_pairing':True,
  'all297_raw33_plane9_coordinate_basis_intertwiners':True,'intertwiner_hash':digest.hexdigest(),
  'two_existing_cache_fourth_coefficients':list(map(str,values)),
  'remaining_formal_callbacks':[list(k)for k in sorted(missing)],'selected_delta22_sensitivity':['1','0','5/3','-121/48'],
  'selectedIFT4_value_gate_or_code_NOT_authorized':True,'no_new_root_composition_or_head29_response_evaluated':True,
  'honest_projection_proof':'P(1)=1 and P*x_i=plane_x_i*P for all nine generators on all33 cyclic basis classes implies P(f)=f(z=0,r,ell) for every function in this finite quotient. Thus raw33 slots31/32 are literal Taylor coefficients even with nonhomogeneous carry.',
  'scope':'Source-only coherent operators and cache availability. Delta63 r2ell2 and bB63 ell3 remain formal here. No head29 edge, original H0 or nil result.',
  'whole_source_numeric_upper':owned,'RAW_bytes':rawbytes}
 ep=out/'independent-coherent-head29-four-operators-and-all297-cache-descent-identities.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'coherent_source_and_cache_availability_PASS':True,'actual_head29_response_evaluated':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'source_owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':main()
