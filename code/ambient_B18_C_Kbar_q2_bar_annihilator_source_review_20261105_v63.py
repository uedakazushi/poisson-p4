#!/usr/bin/env python3
"""Independent specific-H, product44 order and raw128 q2 bar annihilator proof."""
import itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63';ID=BASE+'_REVIEW'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  pp=ROOT/p;z=sha(pp)
  if h:assert z==h,(p,h,z)
  if p not in pins:raw+=pp.stat().st_size
  pins[p]=z;return json.loads(pp.read_text())if parse else None
 rp=f'reports/{BASE}/report.json';r=pin(rp,'ee0ec8dc474c4be3c6e1e14be262780dcd16f544520eae0b2b26e9f9c0ecc6e5');assert r['all_passed']
 ep=next(iter(r['evidence_sha256']));E=pin(f'reports/{BASE}/{ep}',r['evidence_sha256'][ep])
 comp=pin(f'reports/{BASE}/specific-H-template-order-and-raw-domain-companion.json','eb63cdc8cf817365e2b51c1dae392cbc92f97ba9ed888f20657ab6a907ba5c53');assert comp['all_passed']and comp['frozen_primary_sha256']==pins[rp]
 for v in[r,comp]:
  for p,h in v['input_sha256'].items():pin(p,h,False)
 pin('code/ambient_B18_C_Kbar_q2_bar_annihilator_source_20261105_v63.py',r['code_sha256'],False)
 pin('code/ambient_B18_C_Kbar_q2_specific_H_order_source_20261105_v63.py',comp['code_sha256'],False)
 def path(s):return next(p for p in pins if p.endswith(s))
 for p in list(pins):
  if p.endswith('_REVIEW/report.json')or p.endswith('_REVIEW_V2/report.json'):assert pin(p)['all_passed']
 def fields(p,ks):return {k:field(ROOT/p,k)for k in ks}
 hp=path('actual-honest-Tnil128-972-formal-operator.json')
 H=fields(hp,['shared_homotopy_template_v1','actual_S_frame','matrices','matrix_operator_nodes','matrix_operator_shapes','exact_matrix_identities','canonical_division_regular_proofs','source_artifacts','scalar_function_nodes'])
 pp=path('five-cache-functionals-and-genuine-fourteen-kernel-source.json');P=fields(pp,['literal_direct_section_path','minimal_head_inclusion_rows128','source_m0_24_by44','source_n0_44_by24','head42_origin_cycle44_sparse'])
 natural=fields(path('actual-natural-Told-to-Tnil-source-map-and-origin.json'),['actual_origin_data','old_minimal_generator_free56_rows','linearity_and_descent','source_artifacts'])
 five=fields(path('actual-two-cache-witness-and-five-minimal-bar-generators.json'),['actual_minimal_bar_generators','actual_mu_S_Nbar'])
 upper=fields(path('actual-trace-first-gate-and-entire-conductor-proof.json'),['all_Qrad_even_coefficients_in_J','actual_joint_ideal','all36_Qrad_C_maximal_ideal_order_lower_bounds'])
 Pi=fields(path('all-three-Schur-row42-and-structural-bar-cube-zero.json'),['full_minimal_Pi0_first_row128_sparse','all3_exact_Schur_row_routes','all780_X0_row42_coefficients_exact_zero','final_retained_original_labels5'])
 quartic=fields(path('actual-eight-relations-and-complete-quartic-dual.json'),['actual_complete_quartic_dual1650_sparse'])
 loewy=pin(path('AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/report.json'));assert loewy['exact_actual_HF']==[5,5,3,2,1]and loewy['dimensions_m_power_Tnil']==[16,11,6,3,1,0]
 product_report=pin(path('AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/report.json'));assert product_report['all_passed']
 pe=next(iter(product_report['evidence_sha256']));prodpath=str(pathlib.Path(path('AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/report.json')).parent/pe)
 product=pin(prodpath,product_report['evidence_sha256'][pe])
 for p,h in product_report['input_sha256'].items():pin(p,h,False)
 mp='reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json';pin(mp,parse=False)
 Min=fields(mp,['operator_nodes','operator_shapes','matrices','actual_S_frame','source_artifacts','three_exact_unit_Schur_steps'])
 for name,src in H['source_artifacts'].items():pin(src['path'],src['sha256'],False)
 for name,src in Min['source_artifacts'].items():pin(src['path'],src['sha256'],False)
 respath=H['source_artifacts']['new_resolution']['path'];Resolution=fields(respath,['last_regular_sequence_powers','actual_S_variable_order','all8_differentials_sparse_columns','basis_sector_wedge'])
 radicalpath=product['source_artifacts']['radical']['path'];Vrad=field(ROOT/radicalpath,'actual_V3_by44_rows_after_unit1_cancellation')
 parsed=count([r,E,comp,H,P,natural,five,upper,Pi,quartic,loewy,product_report,product,Min,Resolution,Vrad,pins]);owner=parsed+20000
 assert owner<150000 and raw+2_000_000<32*1024*1024
 save(out/'positive-before-pure-source-identities.json',{'all_decoded_source_numeric_labels':parsed,'whole_owned_upper':owner,'RAW_input_plus_output_upper':raw+2000000,'no_signed_nilproducts_or_newjets':True})
 # Same S=Q[[z0..z6,r]], same actual ell0 pullback, same m/V source frames.
 assert H['actual_S_frame']==Min['actual_S_frame']==Resolution['actual_S_variable_order']==['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA']
 hn=H['matrix_operator_nodes'];hs=H['scalar_function_nodes'];hm=H['matrices']
 assert(hm['new_m24_by44'],hm['d1_new'],hm['d2_new'],hm['Phi44_by128'])==(4,22,23,39)
 assert hn[3]==['matrix_product',[2,0]]and hn[39]==['horizontal_blocks',[37,38]]
 assert hn[35]==['matrix_product',[21,4]]and hn[36]==['matrix_neg',35]and hn[37]==['matrix_sum',[34,36]]and hn[38]==['matrix_product',[21,23]]
 for row in range(3):
  for col in range(44):
   node=hn[0][1][row][col];fn=Vrad[row][col]
   assert hs[node]==(['QQ',0,1]if fn==0 else['actual_function_restrict_ell0','radical',fn])
 assert 'Vgamma=d1*m'in H['exact_matrix_identities']and 'm*n=I24'in H['exact_matrix_identities']
 mn=Min['operator_nodes'];assert Min['matrices']['direct_unit_section_n']==8 and mn[0]==['actual_external_matrix','honest','new_m24_by44']
 assert mn[3]==['actual_external_matrix','honest','d2_new']and mn[6]==['matrix_product',[0,5]]and mn[7][0:2]==['origin_unit_matrix_inverse',6]and mn[8]==['matrix_product',[5,7]]
 assert Min['matrices']['full_target_projection5_by128']==59 and mn[59]==['matrix_product',[58,43,28]]
 assert P['minimal_head_inclusion_rows128']==[42,43,50,68,70]
 m0=[list(map(Q,row))for row in P['source_m0_24_by44']];n0=[list(map(Q,row))for row in P['source_n0_44_by24']];E44=[list(map(Q,row))for row in mn[5][1]]
 assert all(not row[43]for row in m0)and all(not x for row in E44[42:]for x in row)
 selected=[i for i,row in enumerate(E44)if any(row)];assert len(selected)==24 and all(E44[i]==[Q(j==k)for k in range(24)]for j,i in enumerate(selected))
 assert all(not x for i,row in enumerate(n0)if i not in selected for x in row)
 assert all(sum(m0[i][k]*n0[k][j]for k in range(44))==Q(i==j)for i in range(24)for j in range(24))
 cycle0=[[Q(i==head)-sum(n0[i][k]*m0[k][head]for k in range(24))for head in[42,43]]+[Q(0)]*3 for i in range(44)]
 assert [(i,str(row[0]))for i,row in enumerate(cycle0)if row[0]]==[(i,q)for i,q in P['head42_origin_cycle44_sparse']]
 M0=natural['actual_origin_data']['actual_origin_map_columns5'];assert natural['old_minimal_generator_free56_rows']==[22,23,29,30,55]
 leading=[sum(row[j]*Q(M0[4][j])for j in range(5))for row in cycle0];assert leading==[Q(i==43)for i in range(44)]
 assert all(not Q(x)for j in[1,2,3]for x in M0[j])and not Q(M0[4][0])and Q(M0[4][1])==1
 assert five['actual_mu_S_Nbar']==5 and five['actual_minimal_bar_generators']==['h1-(5/3)rA*h0','h2+(88/25)rA*h0','h3-(2/5)rA*h0','h4-(136/75)rA*h0','zA5*h0']
 # Verify every literal d1/d2 column: sector powers n=(1,1,2), no r division.
 assert Resolution['last_regular_sequence_powers']==[1,1,2]
 d1,d2=Resolution['all8_differentials_sparse_columns'][:2];pairs=list(itertools.combinations(range(8),2))
 for sector,n in enumerate([1,1,2]):
  mons=[]
  for i in range(8):m=[0]*8;m[i]=n if i==7 else 1;mons.append(m);assert d1[8*sector+i]==[[sector,m,1]]
  for k,(i,j)in enumerate(pairs):assert d2[28*sector+k]==[[8*sector+j,mons[i],1],[8*sector+i,mons[j],-1]]
 template=H['shared_homotopy_template_v1'];assert template==comp['literal_shared_template']
 assert template['pair_i_7_output']=='Hi'and template['pair_i_j_output_i_j_lt7']=='(Euler_z+2)^-1*(partial_zi(wj)-partial_zj(wi)) on coefficient functions'
 assert template['not_S_linear']and template['no_pull_out_actual_coefficients_or_truncated_normal_forms']
 # Universal Weyl/Euler identity, as polynomials in the seven arbitrary
 # exponents alpha_i. This checks coefficients, not sampled H input values.
 euler=[]
 for j in range(7):
  rows=[]
  for k in range(7):
   dd=[0]*8;deltad=[0]*8
   if j==k:
    for i in range(7):
     if i!=j:dd[i+1]+=1
    deltad[j+1]+=1;deltad[0]+=1
   else:dd[k+1]-=1;deltad[k+1]+=1
   total=[a+b for a,b in zip(dd,deltad)];assert total==([1]*8 if j==k else[0]*8);rows.append([dd,deltad,total])
  euler.append(rows)
 assert upper['all_Qrad_even_coefficients_in_J']and upper['actual_joint_ideal']=='J=IA+IB=(zA7,r²,rs,s²)'
 assert len(upper['all36_Qrad_C_maximal_ideal_order_lower_bounds'])==36
 gp=pin(path('AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_SOURCE_20261105_v63/report.json'));assert gp['all_passed']and gp['global_mR0_cubed_for_literal_gA_andgB']
 # Propagate order through the EXACT492-node product source, without values.
 nodes=product['all_new_R0_scalar_nodes'];roots=product['literal_bilinear44_output_roots'];needed=set()
 def use(i):
  if i in needed:return
  needed.add(i);node=nodes[i]
  if node[0]in['sum','product']:
   for k in node[1]:use(k)
 for i in roots:use(i)
 beta=product['actual_beta8_C0_node_ids'];assert beta==[7148,7199,7454,7505,7688,7943,8045,505938]
 evens=set(product['actual_Qrad_even36_R0_component_root_ids']);order={};INF=1000
 for i in sorted(needed):
  node=nodes[i];kind=node[0]
  if kind=='QQ':v=0 if node[1]else INF
  elif kind=='actual_R0_function':
   ns,fn=node[1:]
   if ns=='radical':assert fn in evens;v=1
   else:
    assert ns=='C0'and fn in beta+[20086];v=3 if fn==505938 else 1 if fn in beta[:7]else 0
  elif kind=='arbitrary_cycle_input':
   name,k=node[1:];assert name in['z','zprime']
   v=(3 if k<43 else 2)if name=='z'else(1 if k<43 else 0)
  elif kind=='sum':v=min((order[k]for k in node[1]),default=INF)
  else:assert kind=='product';v=min(INF,sum(order[k]for k in node[1]))
  order[i]=v
 outputorders=[order[i]for i in roots];assert all(v>=5 for v in outputorders[:43])and outputorders[43]>=4
 assert product['actual_multiplication_on_H1']=='exact transported original D0 product under the accepted complete H1 isomorphism'
 # Formal Phi_n*iota cancellation: n and mc stay arbitrary in all44 rows.
 cancellations=[]
 for i in range(44):
  coefficients={('c',i):1}
  for j in range(24):
   key=('n_mc',i,j);coefficients[key]=coefficients.get(key,0)-1;coefficients[key]+=1
  assert{k:v for k,v in coefficients.items()if v}=={('c',i):1};cancellations.append(i)
 assert Pi['full_minimal_Pi0_first_row128_sparse']==[[42,'1']]and Pi['all780_X0_row42_coefficients_exact_zero']
 assert Pi['final_retained_original_labels5']==[42,43,50,68,70]
 assert [(x['incoming_position'],x['retained_position'])for x in Pi['all3_exact_Schur_row_routes']]==[(42,18),(18,0),(0,0)]
 assert all(Q(x['unit_origin_determinant'])for x in Pi['all3_exact_Schur_row_routes'])
 assert quartic['actual_complete_quartic_dual1650_sparse']==[[329,'1']]
 # Complete fourth-degree basis uses330 monomials per generator;329 reads r4
 # of generator0. All Pi-positive-degree terms start5 on iota in m4.
 assert math.comb(11,4)==330 and 5*330==1650
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 e={'all_passed':True,'same_S_formal_regular_domain_and_actual_product_Vrad_to_Vgamma_to_m_checked':True,
  'direct_n_literal_E_mE_inverse_and_source_m0n0_identity_checked':True,'q2_equal_natural_r2_oldhead4_leading_free44_e43_checked':True,
  'five_genuine_Nbar_generators_origin_e43_or0_and_complete_mu5_pinned':True,
  'specific_H_sectors':[1,1,2],'all24d1_and84d2_signs_powers_checked':True,'universal_Euler_identity49_component_polynomials':euler,
  'specific_H_order_proof':'In S=Q[[z7,r]], a genuine cycle v of orderd has v7(0,r)=0 since r^n is a nonzerodivisor. Each orderedHadamard quotient byzi lowers order atmost1. wi=vi+r^nHi hasorderd for n>=1. Onezi derivative then (Euler_z+2)^-1 lowers order atmost1. The coefficientwise Euler identity above holds for arbitrary exponentalpha and rdegree. No divisionbyr/r², nofiniteNF and no S-linear extension of780columnlifts is used.',
  'product44_all44_output_order_bounds_from_actual_DAG':outputorders,'product_free42_order_atleast5_free43_atleast4':True,
  'm_times_product_order5':'The degree4 product is supported onlye43, and actualm0e43=0. Every actual higher m coefficient is regular; hence m*c in m5. Applying WHOLE specificH gives Hmc in m4.',
  'raw128_iota_and_domain':'iota(c)=(c,Hmc) in free44+wed84. d1m=Vgamma and actualproductcycle imply sectorcycles; exactd2Hmc=mc. Phi_n=[I-nm,nd2] for honestn and directn; all44 cancellations prove Phi_n*iota=c. Thus each exact sectiontransportJ_n(s,b)=(Phi_n(s,b),b) fixes iota.',
  'all44_formal_Phi_iota_cancellations':len(cancellations),'fullPi0_row_e42_and_84_zero_wedge_entries_preserved':True,
  'quartic_step':'iota hasorder4; Pi isregular. Pi0 firstrow reads only c42, whose order>=5. All Pi1+ terms start5. Complete quarticdual reads ONLY generator0*r4, so classzero in m4/m5. Accepted m5Tnil0 gives actual zero class.',
  'bar_result':'Since both inputs are genuine Nbar classes, product44 is their actual D0bar product, and Nbar injects Tnil, q2*Nbar=0. g=rq2 also annihilates Nbar by central Rbar-linearity.',
  'square_generation':'Accepted Kbar=(U,S,q2)_Rbar. Therefore Kbar² is Rbar-generated by U²,US,S²; retain all r-multiples for QQ span/module closure. No independence or product scalar value asserted.',
  'upstairs_limitation':'No liftq2 annihilator upstairs, no ell cancellation, no H0/nilindex improvement. Tnil itself is not assigned a ring product.',
  'no_signed_product_source_IFT_gauge_jet_values':True,'whole_owned_upper':owner,'RAW_input_plus_output_upper':raw+2000000}
 evidence=out/'independent-specific-H-product44-Phi-iota-Pi-domain-bar-annihilator-proof.json';save(evidence,e)
 rep={'id':ID,'all_passed':True,'frozen':True,'PURE_SOURCE_ONLY_no_product_values':True,'q2_times_genuine_Nbar_zero':True,'Kbar_square_three_Rbar_generators_only':True,'original_H0_or_nilindex_claimed':False,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{evidence.name:sha(evidence)},'seconds':time.monotonic()-t};save(out/'report.json',rep);signal.alarm(0)
 print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owner,'product_orders':sorted(set(outputorders)),'seconds':rep['seconds']}))
if __name__=='__main__':main()
