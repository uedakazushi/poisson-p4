#!/usr/bin/env python3
"""Independent literal zero7/unit57 and pointed physical-R0 order proof."""
import json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRE=0
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def ck(a,b,m=False):
 global PRE
 if a and b:
  na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
  h=max(na+nb,da+db)if m else max(na+db,nb+da)+1;PRE=max(PRE,h);assert h<=2048
 return a*b if m else a+b
def sparse(a):return[{i:Q(q)for i,q in enumerate(row)if Q(q)}for row in a]
def products(A,B):
 for i,row in enumerate(A):
  v={}
  for k,q in row.items():
   for j,r in B[k].items():
    a=ck(v.get(j,Q(0)),ck(q,r,True))
    if a:v[j]=a
    else:v.pop(j,None)
  assert v=={i:Q(1)}
def determinant9(a):
 a=[[Q(q)for q in row]for row in a];d=Q(1)
 for j in range(9):
  i=next(k for k in range(j,9)if a[k][j])
  if i!=j:a[i],a[j]=a[j],a[i];d=-d
  p=a[j][j];d=ck(d,p,True)
  for k in range(j+1,9):
   if not a[k][j]:continue
   ratio=ck(a[k][j],1/p,True)
   for v in range(j,9):a[k][v]=ck(a[k][v],-ck(ratio,a[j][v],True))
 return d
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
 pr=pin(ROOT/'reports'/PRI/'report.json','107e8c6b4a06b6ea479bbcf6c297a75985e3455e263d0ddfea87b6009e6aa414')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_ell2_conductors_global_m3_source_20261105_v63.py',pr['code_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 add=pin(ROOT/'reports/AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_OWNERSHIP_ADDENDUM_V2_20261105_v63/report.json')
 assert add['all_passed']and add['preserved_primary_sha256']==sha(ROOT/'reports'/PRI/'report.json')
 assert sum(add['explicit_whole_parsed_objects_and_all_temporaries_even_if_retained_together'].values())==add['whole_numeric_upper']==55138<60000
 assert add['whole_RAW_from_primary_plus_addendum_max']==17358586<32<<20
 def src(s):return ROOT/next(p for p in pins if p.endswith(s))
 fj=src('full9-firstjet-certificate.json');monomials=field(fj,'source2_monomials');families=field(fj,'both57_source2')
 Z=field(src('cached-source2-and-exact-E4-quadrics.json'),'families')
 rel=src('full-relative186-DAG.json');instances=field(rel,'all16_pointed_unit57_inverse_instances')
 cut=[k for k,m in enumerate(monomials)if set(m)<={7,8}]
 assert {tuple(monomials[k])for k in cut}=={(),(7,),(8,),(7,7),(7,8),(8,8)}
 actualrows=[27,108,109,112,113,117,118,119,61]
 owned=count([pr,D,add,monomials,families,Z,instances,pins])+11000
 assert owned<60000 and raw+1000000<32<<20
 save(out/'positive-before-fixed-unit-matrix-checks.json',{'all_parsed_parent_blocks_and_aliases_charged':True,'source_private_cap60000':True,'whole_numeric_upper':owned,'no_scalar_response_evaluation':True,'source2_slot_multisets_not_exponents':True})
 records=[]
 for name,root in [('CA',7790),('CB',28539)]:
  f=families[name];z=Z[name]
  assert f['f63']==root and len(f['unknown57'])==len(f['equations57'])==57
  assert len(set(f['unknown57']))==len(set(f['equations57']))==57
  assert f['all57_equations_all55_coefficients_verified']and len(f['pure_function2jets'][str(root)])==55
  assert all(Q(f['pure_function2jets'][str(root)][k])==0 for k in cut)
  assert z['own_F9_coordinate_order']==f['F9_nodes']and len(f['F9_nodes'])==9
  assert z['physical_R0_rows9']==actualrows and z['all9_inverse_composition_identities_through2_checked']
  # The constants prove locality. Their actual first coefficients prove the
  # stored9matrix is the Jacobian of this SAME physical R0 readback.
  readbacks=z['physical_R0_readback_jets_in_own_F9'];inverses=z['inverse_F9_in_common_R0_through2']
  assert all(Q(v[0])==0 for v in readbacks+inverses)
  J9=z['R0_by_F9_unit_Jacobian'];Inv9=z['fixed_inverse9']
  assert [[Q(v[1+i])for i in range(9)]for v in readbacks]==[[Q(q)for q in row]for row in J9]
  assert [[Q(v[1+i])for i in range(9)]for v in inverses]==[[Q(q)for q in row]for row in Inv9]
  products(sparse(J9),sparse(Inv9));products(sparse(Inv9),sparse(J9))
  det=determinant9(J9);assert det==Q(z['unit9_determinant'])and det
  # Unit57 is checked through both saved products. No inverse/IFT solve.
  J57=sparse(f['fixed_J57']);I57=sparse(f['fixed_J57_inverse'])
  products(J57,I57);products(I57,J57);assert Q(f['determinant'])!=0
  records.append({'family':name,'unrestricted_f63':root,'own_F9_nodes':f['F9_nodes'],'zero7_all6_exact_coefficient_slots':cut,'all9_map_and_inverse_constants0':True,'same_physical_R0_rows':actualrows,'own9_determinant':str(det),'both57_stored_inverse_products_PASS':True,'unit57_determinant':f['determinant']})
  assert owned+count([J57,I57,J9,Inv9])+500<60000
 # Literal allorder zero7 names, not unrestricted source f63 labels.
 ca_peer=pin(ROOT/'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63_REVIEW/report.json','273577c4778771f8d467e89ab5a7fd0b1c15ab7afde063f2bd38ee37e239dbda')
 assert ca_peer['allorder_gA505938_prefix0_source_and_primitiveV_gate_accepted']
 c0=ROOT/ca_peer['input_sha256']['C0']['path'];pin(c0,ca_peer['input_sha256']['C0']['sha256'],False)
 prefix=field(c0,'all7_pure_CA_tail_prefix57_systems')[0];tel=field(c0,'all35_regular_coordinate_divisions')[0]
 assert prefix['retained_zero_prefix_length']==0 and Q(prefix['unit57_determinant'])==Q(families['CA']['determinant'])
 assert len(prefix['bindings'])==len(prefix['all57_unit_equations'])==57
 assert tel['original_pure_CA_source_function']==7790 and tel['restriction_r_s_node']==505938
 assert field(c0,'actual_C0_basis8_scalar_nodes')[-1]==505938
 assert [x['coordinate']for x in tel['all7_regular_coordinate_divisions']]==families['CA']['F9_nodes'][:7]
 assert field(rel,'CB_coordinate_Y9_nodes')==families['CB']['F9_nodes']and field(rel,'CB_f63_node')==28539 and field(rel,'CB_g_restriction_node')==299650
 cb=next(x for x in instances if x['retained_Y_indices']==[7,8]);assert len(cb['bindings'])==len(cb['all57_pointed_equations'])==57
 assert Q(cb['unit_Jacobian'])==Q(families['CB']['determinant'])
 assert D['literal_conductor_ids']==[505938,299650]and D['original_unrestricted_f63_ids_NOT_assigned_order3']==[7790,28539]
 assert raw+1000000<32<<20 and all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'both_own57_six_exact_zero7_coefficients_rechecked':True,'all9_physical_and_inverse_constants0_for_both_maps':True,
  'both_actual9_Jacobians_from_readback_firsts_and_both_inverse_products_verified':True,'both_fixed57_inverse_products_verified_without_solve':True,
  'same_physical_R0_rows_for_CA_CB':actualrows,'source_records':records,
  'literal_CA_prefix0_7790_to505938_and_CB_7_8_28539_to299650':True,
  'allorder_prefix_AST_source_inheritance':'CA accepted independent273577 checked all57 prefix equations plus clone root505938 symbolically; CB accepted fullrelative peer checked the literal pointed57 constructor and retainedY7/8 root299650. Current literal systems/root/coordinate locators agree with those exact inputs.',
  'global_order_proof':'In each own F9 chart the zero7 function psi(u,v)=f63(0,...,0,u,v) is defined to all orders by the same unit57 IFT. Its constant, two first and three ordinary second coefficients are all0, hence psi is in(u,v)^3. Its literal pullback is gA505938 or gB299650. Both own-F9 pullbacks are pointed local homomorphisms into the SAME physical R0 because all9 constants0 and the actual unit9 Jacobian agrees with the physical readback. Thus gA,gB lie in m_R0^3 in ALL directions, including zi*r/zi*ell. Ell=5(s-r) is an invertible pointed linear change, so ell² base change preserves this.',
  'only_restricted_literal_gA_gB_not_unrestricted_f63':True,'not_other_g_h_or_nonunit_quotients_or_upper_HF_values':True,
  'producer_curated18882_not_whole_parse_peak_and_addendum55138_adopted':True,
  'own_parsed_parent_blocks_prior':owned,'actual_matrix_check_prebits':PRE,'RAW_unique':raw,
  'no_new_response_or_physical_sourcejet':True}
 ep=out/'independent-unit57-zero7-and-same-pointed-R0-global-m3-proof.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'literal_gA_gB_global_mR0_cubed_proved':True,'unrestricted_f63_order3_claimed':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
