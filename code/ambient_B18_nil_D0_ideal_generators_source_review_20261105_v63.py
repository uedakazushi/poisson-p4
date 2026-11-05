#!/usr/bin/env python3
"""Independent actual odd action, nil-fiber coordinates and D0 ideal minimality."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NIL_D0_IDEAL_GENERATORS_SOURCE_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def field(p,key):
 PINS[p]=sha(ROOT/p);s=(ROOT/p).read_text();tok=json.dumps(key);assert s.count(tok)==1;i=s.index(':',s.index(tok))+1
 while s[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(s,i)[0]
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return list(map(enc,x))
 return x
def pv(i):e=[0]*7;e[i]=1;return {tuple(e):1}
def add(*ps):
 z={}
 for p in ps:
  for m,c in p.items():z[m]=z.get(m,0)+c
 return {m:c for m,c in z.items()if c}
def neg(p):return {m:-c for m,c in p.items()}
def mul(*ps):
 z={(0,)*7:1}
 for p in ps:
  t={}
  for a,c in z.items():
   for b,d in p.items():m=tuple(i+j for i,j in zip(a,b));t[m]=t.get(m,0)+c*d
  z={m:c for m,c in t.items()if c}
 return z
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_NIL_D0_IDEAL_GENERATORS_SOURCE_20261105_v63/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p in r['input_sha256']:
  if p.endswith('/report.json'):assert read(p)['all_passed']
 S=read('reports/AMBIENT_B18_NIL_D0_IDEAL_GENERATORS_SOURCE_20261105_v63/literal-odd-eight-action-and-five-dimensional-nil-fiber-proof.json')
 binding=read('reports/AMBIENT_B18_NIL_D0_ACTION_INPUT_BINDING_20261105_v63/report.json');chain=read('reports/AMBIENT_B18_NIL_D0_ODD_ACTION_CHAIN_IDENTITY_20261105_v63/report.json')
 for receipt in[binding,chain]:
  assert receipt['all_passed']
  for p,h in receipt['input_sha256'].items():assert sha(ROOT/p)==h
 P=read('reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/actual-nil-product44-and-cube-DAG.json')
 cp=P['source_artifacts']['C0']['path'];assert sha(ROOT/cp)==P['source_artifacts']['C0']['sha256']
 # Select only source-interface fields, NEVER decode the567193-node C0 DAG here.
 beta=field(cp,'actual_C0_basis8_scalar_nodes');images=field(cp,'actual_C0_free9_images_in_CA_trace_pairs');relations=field(cp,'honest_D0_relation_columns')
 assert beta==[7148,7199,7454,7505,7688,7943,8045,505938]
 assert images==[[1,0]]+[[0,b]for b in beta]and len(relations)==640 and all(len(col)==9 for col in relations)
 assert field(cp,'D0_relation_origin_rank')==0
 assert field(cp,'cached_actual_CA_gA_all_degree3_coefficients')[:6]==['0']*6
 cpeer=read('reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json');assert cpeer['actual_D0_free9_relations640_and_D4_map_checked']
 assert cpeer['input_sha256'][cp]==sha(ROOT/cp)
 assert P['actual_beta8_C0_node_ids']==beta and P['actual_nuA_C0_node_id']==20086
 # Bind b0 to the SAME a36 inputs as the new odd-action operator.
 pool=P['all_new_R0_scalar_nodes'];broot=P['cycle_representative_even_b0'];assert broot==171
 assert pool[broot][0]=='sum'and len(pool[broot][1])==36
 for j,k in enumerate(pool[broot][1]):
  assert pool[k][0]=='product';f,z=pool[k][1]
  assert pool[f]==['actual_R0_function','radical',P['actual_Qrad_even36_R0_component_root_ids'][j]]
  assert pool[z]==['arbitrary_cycle_input','z',j]
 nodes=S['actual_odd_action_nodes'];roots=S['all8_odd_action44_roots']
 assert nodes[:3]==[['QQ',0],['actual_R0_function','C0',20086],['actual_nil_even_b0','product',broot]]
 assert nodes[3:39]==[['nil_a36_input',j]for j in range(36)]
 # Tiny commutative monomial normal form of EVERY352 output, with b0 symbolic
 # but bound above to the actual product171 source expression.
 pol=[]
 for node in nodes:
  tag=node[0]
  if tag=='QQ':assert node[1]==0;p={}
  elif tag=='actual_R0_function':p={(('nu',)if node[2]==20086 else ('beta'+str(beta.index(node[2])),)):1}
  elif tag=='actual_nil_even_b0':p={('b0',):1}
  elif tag=='nil_a36_input':p={('a'+str(node[1]),):1}
  elif tag=='product':
   p={():1}
   for k in node[1]:
    p={tuple(sorted(a+b)):c*d for a,c in p.items()for b,d in pol[k].items()}
  else:raise AssertionError(node)
  pol.append(p)
 for i,outroots in enumerate(roots):
  expected=[]
  for g in range(9):
   for src,hasnu in[(4*g+1,True),(4*g,False),(4*g+3,True),(4*g+2,False)]:expected.append({tuple(sorted((('nu',)if hasnu else ())+('beta'+str(i),'a'+str(src)))):1})
  expected +=[{('b0',):1}if j==i else{}for j in range(8)]
  assert [pol[k]for k in outroots]==expected
 # Distinct free seven-variable engine verifies the off-cycle equations.
 q0,qA,qB,qAB,B,nu,beta_poly=[pv(i)for i in range(7)]
 qa=[mul(nu,beta_poly,qA),mul(beta_poly,q0),mul(nu,beta_poly,qAB),mul(beta_poly,qB)]
 v=[add(qA,neg(B)),qB,qAB];newV=[add(qa[1],neg(mul(beta_poly,q0))),qa[2],qa[3]]
 assert newV==[{},mul(nu,beta_poly,v[2]),mul(beta_poly,v[1])]
 assert add(qa[0],neg(mul(nu,beta_poly,B)))==mul(nu,beta_poly,v[0])
 assert newV==[{tuple(m):int(c)for m,c in row}for row in chain['V_action3']]
 fp='reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json'
 first=field(fp,'all849_actual_minimal_first_columns40_sparse');L={0:Q(1),4:Q(25,4),5:Q(1)}
 assert S['fifth_functional_Lz']==[[i,str(c)]for i,c in L.items()]
 for col in first:
  col={int(i):Q(c)for i,c in col};assert sum(c*col.get(i,Q(0))for i,c in L.items())==0
 mu=read('reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63_REVIEW/report.json');assert mu['actual_mu_R0_N']==mu['actual_mu_S_Nbar']==5
 deep=read('reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_REDUCTION_20261105_v63/actual-bar-cube-order-four-and-leading-source-proof.json');assert deep['all14_a36_and_zero7_odd_coordinate_origins_zero']
 Pi=read('reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63_REVIEW/independent-Pi0-row42-and-all-genuine-bar-cubes-zero-proof.json');assert Pi['independent_full_Pi0_firstrow128_sparse']==[[42,'1']]
 data={'all352_literal_odd_action_outputs_independently_rebuilt':True,'input_namespaces_bound_to_same_actual_product171_and_a36':True,'all36_b0_source_terms_bound':True,'offcycle_chain_and_even_identity_free_polynomial_checks':True,
 'C0_actual_same_base_surjection':'Accepted actual C0 is R0+(z7,gA)*tauA, not a free9 algebra; its free9 module cover and complete640 relations give the SAME zero-degree D0 quotient. Exact source beta last is curve505938, not rawf637790.',
 'localmax_generation_proof':'All beta_i*tauA vanish at the actual closed point. More intrinsically their pairwise products nuA*beta_i*beta_j lie m_R0, so the finite R0-fiber is a quotient of k plus squarezero8. It is local with residue k; hence m_D0=m_R0D0+sum R0*(beta_i*tauA). No positive-graded extension is asserted.',
 'all849_Lz_annihilations':True,'Lz':[[i,str(c)]for i,c in L.items()],
 'fiber_I5_proof':'Genuine Nbar has v0(0)=0. Four head origins and Lz(v0_first) descend through all849 relations and kill m_SNbar; the four correctedheads and z5h0 give I5 because Lz(r)=0,Lz(z5)=1. Actual mu_SNbar5 makes them complete fiber coordinates.',
 'canonical_nil_action_order_proof':'The five canonical barcycle lifts have a36∈m, Qrad_even∈m gives b0∈m², beta_i∈m. All44 action coordinates therefore lie m². Membership m is regular, whole-function H loses≤1, so (cprime,Hmcprime) lies m. The minimal origin is0, and head0 first part is Pi0row0*cprime_first=e42*cprime_first=0; Pi1 cannot contribute because the honest origin is0. Thus all5 fiber functionals vanish on every8 odd-generator action.',
 'genuine_bar_action_descent':'Lift each canonical Nbar class to a full nilcycle and correct its literal reduction by a lifted boundary. The polynomial odd action represents actual C0 multiplication and preserves cycles by the checked chain identity. Full boundary images represent0 in D0, so it induces the actual action on Nbar, without multiplying arbitrary Tnil classes.',
 'ideal_minimality_proof':'All8 odd generators kill N/m_R0N on the5 minimal R0 generators; R0-linearity gives odd_i*N⊂m_R0N. The actual maximalideal decomposition gives m_D0N=m_R0N. Both residue fields are k, so mu_D0idealN=mu_R0N=5.',
 'actual_minimal_D0_ideal_generator_count':5,'scope':'zero-degree D0 only; arbitrary true nil lifts of the5 bar generators are minimal, but full coefficient lifts are not built; whole positive-graded nil ideal, products and index are not computed','no_scalar8_by5_action_or_new_source_gauge_jet':True,'bounded_worker_uses_selective_C0_fields_not_full_DAG':True}
 ep=out/'independent-eight-actions-correct-fiber-and-D0-ideal-five-proof.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_independent_actual_D0_nilideal_minimal_generator_count5','actual_mu_D0_ideal_N':5,'actual_mu_R0_N':5,'whole_positive_graded_nil_count_claimed':False,'input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
