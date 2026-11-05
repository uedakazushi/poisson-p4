#!/usr/bin/env python3
"""Independent original connecting source/types/plane orders; no edge values."""
import hashlib,itertools,json,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha,mon,mm,mons
ROOT=pathlib.Path(__file__).resolve().parents[1]
FIVE='AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63'
PLANE='AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def fastfield(p,k):
 # Only a unique literal field is decoded; large scalar arrays remain text.
 text=p.read_text();m=list(re.finditer(r'"'+re.escape(k)+r'"\s*:',text));assert len(m)==1
 pos=m[0].end()
 while text[pos].isspace():pos+=1
 return json.JSONDecoder().raw_decode(text,pos)[0]
def prefix(p,k,n):
 text=p.read_text();matches=list(re.finditer(r'"'+k+r'"\s*:',text));assert len(matches)==1
 pos=matches[0].end()
 while text[pos].isspace():pos+=1
 assert text[pos]=='[';pos+=1;dec=json.JSONDecoder();rows=[]
 for i in range(n):
  while text[pos].isspace() or text[pos]==',':pos+=1
  row,pos=dec.raw_decode(text,pos);rows.append(row)
 return rows
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 outputs=[ROOT/'reports'/(x+'_REVIEW') for x in [FIVE,PLANE]]
 for out in outputs:out.mkdir(exist_ok=False)
 pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  name=str(p.relative_to(ROOT))
  if name not in pins:RAW+=p.stat().st_size
  pins[name]=s
  return json.loads(p.read_text()) if parse else None
 reports=[];docs=[]
 for stem in [FIVE,PLANE]:
  pr=pin(ROOT/'reports'/stem/'report.json')
  for p,h in pr['input_sha256'].items():
   pin(ROOT/p,h,False)
   if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
  code=ROOT/'code'/('ambient_B18_C_connecting_five_formal_source_20261105_v63.py' if stem==FIVE else 'ambient_B18_C_connecting_plane_order_source_20261105_v63.py')
  pin(code,pr['code_sha256'],False);pin(ROOT/'research/notes'/f'{stem}.md',pr['note_sha256'],False)
  en=next(iter(pr['evidence_sha256']));docs.append(pin(ROOT/'reports'/stem/en,pr['evidence_sha256'][en]));reports.append(pr)
 F,P=docs;rf,rp=reports
 assert sha(ROOT/'reports'/FIVE/'report.json')=='8d11f0e8f70a023faba8f58ecb98b2959cc90ead92e2c8560e42ba31c64e2b95'
 assert rp['sigma2_rr_adjoint_retained'] and not rp['actual_connecting_values_computed']
 qp=ROOT/'reports/AMBIENT_B18_C_GENUINE_PROJECTION7_SOURCE_20261105_v63_REVIEW/report.json'
 qr=pin(qp,'36e17392dfd4abba9aa69132065a2809e2cc6f73e1bdc7d7ba2f225ff349db10');assert qr['all_passed']
 assert F['source_namespace_bindings']['q7']['sha256']=='4894cf94c936546773826182c859ccfba7323f882ce381db8c17e8e59a9f6a05'
 def fpath(key):return ROOT/F['source_namespace_bindings'][key]['path']
 chain=fpath('chain');unit=fpath('unit');nat=fpath('natural');clp=fpath('clones')
 Q48=field(chain,'exact_existing_Q48_C4');U9=field(chain,'exact_existing_U9_C4')
 assert len(Q48)==48 and all(len(row)==4 for row in Q48)
 assert U9==[[1,0,0,0]]+[[0,j+2,0,0] for j in range(8)]
 full=[[[i,1] for i in row[1:]] for row in Q48]
 full += [[[i,-1 if i else 1] for i in row[1:]] for row in U9[1:]]
 assert full==F['full_V56_reference_pairs_columns56_targets3']
 assert len(full)==56 and sum(len(x) for x in full)==168
 ops=field(unit,'typed_regular_function_operators');old=field(unit,'actual_V56_function_operator_refs')
 for col in range(56):
  for t in range(3):
   op=ops[old[col][t]]
   if op==['zero']:pair=[0,1]
   else:
    assert op[0] in ['actual_chain_scalar_restrict_ell0','negative_actual_chain_scalar_restrict_ell0']
    pair=[op[1],-1 if op[0].startswith('negative') else 1]
   assert pair==full[col][t]
 roots={i for row in full for i,s in row};assert max(roots)==2160
 original=prefix(chain,'operator_scalar_nodes',2161);peakprefix=count(original)
 needed=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);o=original[i]
  if o[0]=='c':assert len(o)==3 and o[2]!=0;children=[]
  elif o[0]=='actual_C0_scalar_function':assert len(o)==2;children=[]
  elif o[0]=='sum':
   assert all(d!=0 for j,n,d in o[1]);children=[j for j,n,d in o[1]]
  else:
   assert o[0]=='product' and len(o)==4 and o[2]!=0;children=o[3]
  assert all(0<=j<i for j in children);todo.extend(children)
 selected=[[i,original[i]] for i in sorted(needed)]
 assert len(selected)==269 and selected==F['selected_original_chain_scalar_nodes']
 first10=original[:10];del original
 # Exact m,n unit-section source, rather than a dense inverse expansion.
 m=field(unit,'typed_m24_by56');n=field(unit,'typed_n56_by24');sel=field(unit,'fixed_selected24_columns')
 inv=field(unit,'typed_unit24_inverse_node');det=Q(field(unit,'fixed24_minor_determinant'))
 assert det==Q(-397267578125,50096498540544) and len(sel)==24 and len(set(sel))==24
 assert len(m)==24 and all(len(row)==56 for row in m)
 assert len(n)==56 and all(len(row)==24 for row in n)
 assert ops[inv][0]=='formal_unit_matrix_inverse'
 assert ops[inv][1]==[[m[i][j] for j in sel] for i in range(24)]
 for i in range(56):
  for j in range(24):
   if i in sel:assert ops[n[i][j]]==['unit_matrix_inverse_entry',inv,sel.index(i),j]
   else:assert ops[n[i][j]]==['zero']
 heads=field(nat,'old_minimal_generator_free56_rows');assert heads==[22,23,29,30,55]
 native=field(nat,'matrix_operator_nodes');shapes=field(nat,'matrix_operator_shapes')
 E=[[str(int(i==h)) for h in heads] for i in range(56)]
 assert native[:3]==[['actual_external_matrix','old_unit','typed_m24_by56'],['actual_external_matrix','old_unit','typed_n56_by24'],['QQ_matrix_rows',E]]
 assert native[3:8]==[['identity'],['matrix_product',[1,0]],['matrix_neg',4],['matrix_sum',3,5],['matrix_product',[6,2]]]
 assert shapes[:8]==[[24,56],[56,24],[56,5],[56,56],[56,56],[56,56],[56,56],[56,5]]
 # Distinct same105 clones define a literal coefficient-ring section.
 clones=field(clp,'two_distinct_original105_source_clones')
 assert [x['namespace'] for x in clones]==['FULL_R0_F9','ELL0_S_F9']
 assert clones[0]['target9'][:8]==clones[1]['target9'][:8]
 assert clones[0]['target9'][-1]==['sum',['r'],['scale','1/5',['ell']]]
 assert clones[1]['target9'][-1]==['r']
 assert all(x['equation_target_indices']==list(range(48,57)) for x in clones)
 assert all(x['origin_Jacobian_determinant']=='1202315964973056/390625' for x in clones)
 assert F['same_original105_two_separate_source_clones']==clones
 # Independent type inference for every 17 constructor nodes.
 nodes=F['matrix_operator_nodes'];declared=F['matrix_operator_shapes'];got=[]
 assert len(nodes)==17
 for ix,o in enumerate(nodes):
  op=o[0]
  if op=='actual_external_ELL0_matrix':
   assert o[-1]=='ELL0_S_F9';s=[24,56] if o[2]=='typed_m24_by56' else [56,24]
  elif op=='QQ_matrix_rows':
   assert len({len(row) for row in o[1]})==1;s=[len(o[1]),len(o[1][0])]
  elif op=='identity_matrix':s=[o[1],o[1]]
  elif op=='matrix_product':
   a,b=got[o[1]],got[o[2]];assert a[1]==b[0];s=[a[0],b[1]]
  elif op=='matrix_difference':assert got[o[1]]==got[o[2]];s=got[o[1]]
  elif op=='ell_independent_ring_section_Rbar_into_R0':assert o[1]==6;s=got[o[1]]
  elif op=='actual_full_V3_by56_from_original_chain_Q48_negative_U8':assert o[1]==full and o[-1]=='FULL_R0_F9';s=[3,56]
  elif op=='entrywise_formal_RegDiv_ell_before_finite_quotient':assert o[1]==9 and o[2:] == ['ell','Vbar_sigma=0'];s=got[o[1]]
  elif op=='raw_coefficient_extraction_through_degree3':assert got[o[1]]==[3,5];s=[660,5]
  elif op=='actual_external_Dlinear_q7_sparse':assert got[o[3]]==[660,5];s=[7,5]
  elif op=='select_rows_of_genuine_ker_E7':assert o[2]==[1,2,5,6] and got[o[1]]==[7,5];s=[4,5]
  else:assert op=='actual_three_chi_constraints_from_source7' and got[o[1]]==[4,5];s=[3,1]
  assert s==declared[ix];got.append(s)
 assert nodes[2]==['QQ_matrix_rows',E] and nodes[4:8]==[['matrix_product',1,0],['matrix_difference',3,4],['matrix_product',5,2],['ell_independent_ring_section_Rbar_into_R0',6,'DO_NOT_REPLACE_BY_FULL_SOURCE_m_n']]
 assert nodes[9]==['matrix_product',8,7] and nodes[15]==['matrix_product',14,13]
 # Full boundary descent is inherited from source-verified identities, not
 # from tensoring homology or assuming a finite kernel equals torsion.
 assert field(chain,'all640_UP_equals_QB_verified')
 assert field(chain,'all28_intrinsic_C0_columns_have_zero_lift')
 assert field(chain,'all612_source_pair_and_honest_lift_bindings') # actual references retained
 assert '128' in field(unit,'actual_Xfull_source') and '640' in field(unit,'actual_Xfull_source')
 # Correct U9 versus OLD640 namespace from literal immutable records.
 c0=ROOT/next(p for p in rp['input_sha256'] if p.endswith('honest-C0-D0-D4-DAG.json'))
 c0basis=fastfield(c0,'actual_C0_basis8_scalar_nodes')
 assert c0basis==[7148,7199,7454,7505,7688,7943,8045,505938]
 assert first10[2:10]==[['actual_C0_scalar_function',i] for i in c0basis]
 reg=fastfield(c0,'all35_regular_coordinate_divisions')[0]
 ids=field(chain,'accepted_actual_regular_function_identities')
 tail=next(x for x in ids if x['kind']=='accepted_actual_CA_f63_tail')
 assert reg['original_pure_CA_source_function']==tail['function']==7790
 assert reg['restriction_r_s_node']==tail['gA']==505938
 assert tail['coordinate_quotients']==[x['node'] for x in reg['all7_regular_coordinate_divisions']]
 correction=P['namespace_correction_to_frozen_five_formal_prose']
 assert correction['actual_C0_basis8_original_nodes']==c0basis and correction['exact_tail_identity_record']==tail
 gA=[Q(x) for x in fastfield(c0,'cached_actual_CA_gA_all_degree3_coefficients')]
 assert len(gA)==10 and not any(gA[:6]) and any(gA[6:])
 alias={x['physical_odd_coefficient']:x['odd_generator'] for x in ids if x['kind']=='exact_QQ_linear_source_alias'}
 # Rebuild all56 plane low-order entries by multiplication in the free C4
 # algebra: tauA²=nuA,tauB²=nuB, Fj=delta+bA*tauA-bB*tauB.
 c2=ROOT/next(p for p in rp['input_sha256'] if p.endswith('actual-physical12-and-C2-right-exact-quotient.json'))
 rawp=ROOT/next(p for p in rp['input_sha256'] if p.endswith('actual-ell2-product33-and-raw45-actions.json'))
 ellp=ROOT/next(p for p in rp['input_sha256'] if p.endswith('all164-cached-pure-ell-readbacks-and-physical12-second-values.json'))
 framep=ROOT/next(p for p in rp['input_sha256'] if p.endswith('actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json'))
 raw33=field(rawp,'actual33_basis');physical=field(c2,'physical12_coefficients')
 ell={x['physical_row']:x for x in field(ellp,'all12_actual_physical_coefficients_pure_ell3')}
 assert field(framep,'both_nu_constant_and_pureplane_first_zero')
 B2=[[Q(0)]*56 for _ in range(3)];H2=[[Q(0)]*56 for _ in range(3)];checks=[]
 first=[raw33.index(x) for x in [[0,1,0],[0,0,1]]]
 second=[raw33.index(x) for x in [[0,2,0],[0,1,1]]]
 for j,rec in enumerate(physical):
  row=rec['physical_row'];v={k:{int(i):Q(x) for i,x in rec[k]} for k in ['delta','bA','bB']}
  for x in v.values():assert all(not x.get(i,Q(0)) for i in [0]+first)
  for k in ['delta','bA']:assert all(not v[k].get(i,Q(0)) for i in second) and Q(ell[row][k+'_jet3'][2])==0
  bb=[v['bB'].get(i,Q(0)) for i in second]+[Q(ell[row]['bB_jet3'][2])]
  for k in range(3):B2[k][4*j]=-bb[k];H2[k][4*j+1]=-bb[k]
  if j<7:assert alias[rec['roots']['bA_node']]==c0basis[j]
  checks.append([row,list(map(str,bb))])
 assert B2==[[Q(x) for x in row] for row in P['actual_plane2_B_three_monomials_rows3_columns56']]
 assert H2==[[Q(x) for x in row] for row in P['actual_plane2_H_three_monomials_rows3_columns56']]
 adj=[[j,str(B2[2][j]-Q(96,275)*B2[1][j])] for j in range(56) if B2[2][j]-Q(96,275)*B2[1][j]]
 assert adj==P['sigma2_surviving_B_high_coordinate_adjoint'] and [x[0] for x in adj]==[0,16,20]
 # Reconstruct the19 required numerator selectors from independent q7.
 qpath=fpath('q7');qcols=field(qpath,'projector7_by660_sparse_columns');needed_raw=[]
 for col in qcols:
  weights=[[i,a] for i,a in col['target7_affine'] if i in [1,2,5,6]]
  if weights:needed_raw.append({'raw_target':col['raw_target'],'numerator_monomial':list(mm(tuple(col['monomial']),mon(ell=1))),'q7_free_rows_affine':weights})
 assert needed_raw==F['required_raw_numerator_coefficient_functionals'] and len(needed_raw)==19
 route=[];survive=set();drops=[]
 for req in needed_raw:
  M=tuple(req['numerator_monomial']);t=req['raw_target'];terms=[]
  # Enumerate by componentwise divisors, rather than scanning a generic jet.
  for sm in itertools.product(*(range(x+1) for x in M[:8]),[0]):
   vm=tuple(x-y for x,y in zip(M,sm));n=sum(vm)
   if n<=1:reason='V0=0 or V1_ell=0'
   elif n==2 and not any(vm[:7]) and t==0:reason='V2_A,plane=0'
   else:reason=None
   if reason:drops.append([t,list(M),list(sm),reason])
   else:
    survive.add(sm);terms.append({'sigma_monomial':list(sm),'V_monomial':list(vm),'V_order':n})
  route.append((t,M,terms))
 def canonical(paths):return sorted(json.dumps(x,sort_keys=True) for x in paths)
 for (t,M,terms),saved in zip(route,P['all_selected_q7_product_paths']):
  assert t==saved['raw_target'] and list(M)==saved['numerator_monomial'] and canonical(terms)==canonical(saved['retained_product_paths'])
 assert canonical(drops)==canonical(P['formally_dropped_paths_with_literal_reason'])
 assert survive=={mon(),mon(r=1),mon(z=0),mon(z=4),mon(z=5),mon(r=2)}
 assert not any(sum(x)>=3 for x in survive)
 high=[x for x in needed_raw if sum(x['numerator_monomial'])==4]
 assert high==F['highest_four_raw_numerator_monomials'] and len(high)==4
 assert all(x['raw_target'] in [0,1] and tuple(x['numerator_monomial'][7:]) in [(3,1),(2,2)] for x in high)
 owned=count([reports,docs,ops,old,m,n,native,shapes,clones,full,selected,heads,E,nodes,declared,c0basis,reg,ids,raw33,physical,ell,B2,H2,qcols,needed_raw,route,drops,checks])+peakprefix+15000
 assert owned<150000 and RAW+1000000<32<<20
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={
  'all168_full_V_to_ell0_literal_references_and_signs':True,'selected_original_scalar_DAG_nodes':269,
  'typed_unit_m_n_and_same_natural_five_cycles_checked':True,'all17_operator_nodes_shapes_and_foreign_namespaces':True,
  'same105_full_and_ell0_clones_with_nonzero_Jacobian':True,
  'regular_division_and_genuine_target_proof':[
   'The accepted mbar is actual: d1*mbar=Vbar. The literal selected24 unit block gives nbar=E_S*(mbar E_S)^-1 and mbar*nbar=I24. Thus Vbar*(I-nbar*mbar)E=0 for the actual five heads.',
   'Lift the FINAL ELL0 coefficients through the ring section Rbar→R0, without replacing mbar/nbar by FULL counterparts. Each of the15 FULL numerators restricts to0 and is uniquely divisible byell in the formal regular ring. This division precedes finite q.',
   'For eta=numerator/ell, ell[eta]=0 in actual C. The independent D-only q7 therefore yields kerE7; rows1,2,5,6 and the accepted genuine Jordan inverse give C[ell], not the finiteker alone.',
   'For any mod-ell boundary lift the FULL original X and coefficients. FULL V X=0 follows from the accepted actual640 UP=QB columns plus128 complete native normalizer syzygies. An ell-independent boundary lift differs byell*t, so its eta differs byV*t and is zero inC. This proves descent for all768 paths and their825 Schur/Euler images, without a flatness or homology-base-change shortcut.',
   'The formal coefficient section and this boundary argument give Rbar-linearity of the CONNECTING MAP. Rbar target actions use the accepted genuine transported module, never a presumed R0-linear action onD-only T7.'
  ],
  'precise_namespace_correction':'Actual U9 uses honest C0 zA0..6,gA505938; source scalar9 is gA. OLD640 CA multipliers use f63=7790, related by the accepted full seven-coordinate telescope f63=gA+sum(e_i*zAi). Frozen FIVE_FORMAL old_eighth prose is overbroad and is superseded ONLY there by PLANE_ORDER correction. All actual168 references, sigma and RegDiv remain unchanged.',
  'all56_origin_and_pureplane_first_zero_Aplane_second_zero':True,
  'all56_B_H_plane_second_matrices_from_original_C4_products':True,
  'surviving_sigma2_B_high_adjoint':adj,
  'all19_selected_numerator_functionals_and_coherent_convolution_routes':True,
  'only_sigma_jets':['0','r','z0','z4','z5','r2'],
  'sigma3_drops_but_sigma2_rr_does_not':True,
  'open_highest_frontier':'V4_A/B*sigma0 at r3ell/r2ell2, coherent V3*sigma1 and mixed V2*sigma1, and the NONZERO three-support B2*sigma2_rr remain. No general source4 availability or sigma2 cancellation is claimed.',
  'no_connecting_new_source_affine_or_jet_values_computed':True,
  'whole_owned_labels_and15000reserve_and_prefix':owned,'RAW_unique_inputs':RAW,
  'caps':{'numeric':150000,'RAW':32<<20,'seconds':60,'memory':2<<30}}
 result=[]
 for stem,out in zip([FIVE,PLANE],outputs):
  ep=out/'independent-full-source-formal-division-and-plane-order-proof.json';save(ep,proof)
  r={'id':stem+'_REVIEW','all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
   'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
   'actual168_references_five_cycles_formal_division_descent_and_plane_routing_PASS':True,
   'connecting_values_computed':False,'seconds':time.monotonic()-start}
  save(out/'report.json',r);result.append([stem,sha(out/'report.json')])
 signal.alarm(0);print(json.dumps({'all_passed':True,'reports':result,'owned':owned,'RAW':RAW,'seconds':time.monotonic()-start}))
if __name__=='__main__':main()
