#!/usr/bin/env python3
"""Independent conormal adjoints and allorder primitive-V source gate."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','profile':f'reports/{A}/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json',
 'betti':'reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json',
 'betti_peer':'reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 'resolution':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
 'resolution_peer':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'chain_peer':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'actual20':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'actual20_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
 'physical_catalog':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
class Terms:
 """Canonical rational expression DAG; variables use node IDs, not names."""
 def __init__(self):self.nodes=[];self.cache={};self.zero=self.c(0);self.one=self.c(1)
 def put(self,key):
  if key not in self.cache:self.cache[key]=len(self.nodes);self.nodes.append(key)
  return self.cache[key]
 def c(self,c):return self.put(('c',Q(c)))
 def neg(self,a):
  n=self.nodes[a]
  return self.c(-n[1]) if n[0]=='c' else n[1] if n[0]=='n' else self.put(('n',a))
 def plus(self,a,b):
  if a==self.zero:return b
  if b==self.zero:return a
  if self.nodes[a]==('n',b) or self.nodes[b]==('n',a):return self.zero
  if self.nodes[a][0]==self.nodes[b][0]=='c':return self.c(self.nodes[a][1]+self.nodes[b][1])
  return self.put(('a',*sorted((a,b))))
 def mul(self,a,b):
  if self.zero in [a,b]:return self.zero
  if a==self.one:return b
  if b==self.one:return a
  if self.nodes[a][0]==self.nodes[b][0]=='c':return self.c(self.nodes[a][1]*self.nodes[b][1])
  return self.put(('m',*sorted((a,b))))
 def inv(self,a):return self.c(1/self.nodes[a][1]) if self.nodes[a][0]=='c' else self.put(('i',a))
 def eval(self,nodes,roots,binding):
  needed=set();todo=list(roots)
  while todo:
   i=todo.pop()
   if i in needed:continue
   needed.add(i)
   if i in binding:continue
   n=nodes[i]
   if n[0] in ['a','m','n','i']:todo.extend(n[1:])
   else:assert n[0] in ['c','v']
  vals={}
  for i in sorted(needed):
   n=nodes[i];op=n[0]
   if i in binding:v=binding[i]
   elif op=='v':v=self.put(('v',i))
   elif op=='c':v=self.c(Q(n[1],n[2]))
   elif op=='a':v=self.plus(vals[n[1]],vals[n[2]])
   elif op=='m':v=self.mul(vals[n[1]],vals[n[2]])
   elif op=='n':v=self.neg(vals[n[1]])
   elif op=='i':v=self.inv(vals[n[1]])
   else:raise AssertionError(n)
   vals[i]=v
  return {i:vals[i] for i in roots},len(needed)
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'CONORMAL_AND_LITERAL_GA_PLANE_SOURCE_GATES','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert all(D[k]['all_passed'] for k in ['primary','betti_peer','resolution_peer','chain_peer','C0_peer','actual20_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 profile=D['profile'];d2=list(map(sp,D['betti']['all_d2_constant_columns']));anchors=D['resolution']['source_anchors11'];selected=D['betti']['fixed_original_minor_columns']
 fixed=[{i:Q(-1)} for i in anchors]+[d2[i] for i in selected];basis={}
 for receipt in profile['fixed100_elimination_receipts']:
  p=fixed[receipt['fixed_basis_column']]
  for k,c in receipt['subtract']:p=add(p,scale(basis[k],-Q(c)))
  k=receipt['pivot'];assert k not in basis and p.get(k)==Q(receipt['divide_by'])
  p=scale(p,1/Q(receipt['divide_by']));assert enc(p)==receipt['normalized'];basis[k]=p
 assert len(basis)==100
 Pi=list(map(sp,profile['Pi26_source126_columns']));free=profile['conormal_standard_source_indices26'];assert free==[i for i in range(126) if i not in basis] and len(free)==26
 assert all(Pi[i]=={j:Q(1)} for j,i in enumerate(free))
 assert all(not mv(Pi,p) for p in fixed+d2)
 actions=[[sp(p) for p in M] for M in D['C2']['C2_coordinate_actions9_sparse_columns']];eta=[0,5,10]
 monomials=[]
 for total in range(1,5):
  for ell in range(2):
   for indices in itertools.combinations_with_replacement(range(8),total-ell):
    m=[0]*9
    for i in indices:m[i]+=1
    m[8]=ell;monomials.append(tuple(m))
 assert len(monomials)==659
 def mon_action(m,v):
  for i,n in enumerate(m):
   for _ in range(n):v=mv(actions[i],v)
  return v
 library=[];counts={}
 for physical,pos in enumerate(eta):
  for mon in monomials:
   if mon[8]:gamma=mon_action(mon[:8],{pos:Q(1)})
   else:
    i=max(i for i,n in enumerate(mon[:8]) if n);rest=list(mon[:8]);rest[i]-=1
    gamma={14+14*i+j:c for j,c in mon_action(rest,{pos:Q(1)}).items()}
   p=mv(Pi,gamma)
   if p:library.append({'physical_target':physical,'monomial9':list(mon),'conormal26_value':enc(p)});counts[sum(mon)]=counts.get(sum(mon),0)+1
 assert library==profile['scalar_monomial_adjoint_library'] and counts=={1:21,2:7,3:4}
 assert all(sum(x['monomial9'])<4 for x in library)
 chain=D['chain'];V=[[[i,'1'] for i in c[1:]] for c in chain['exact_existing_Q48_C4']]+[[[i,'-1'] for i in c[1:]] for c in chain['exact_existing_U9_C4'][1:]]
 assert V==profile['actual_V56_source_columns'] and len(V)==56
 request={}
 for j,col in enumerate(V):
  for entry in library:
   root,sign=col[entry['physical_target']]
   if not root:continue
   key=(j,root,tuple(entry['monomial9']));request[key]=add(request.get(key,{}),scale(sp(entry['conormal26_value']),Q(sign)))
 actual=[{'V_column':j,'chain_scalar_root':r,'monomial9':list(m),'conormal26_weight':enc(p)} for (j,r,m),p in sorted(request.items()) if p]
 assert actual==profile['requested_actual_V_source_coefficients']
 # Independent allorder prefix0 AST binding, including every57 equation.
 S=D['source'];C=D['C0'];ca=S['families']['CA'];FJ=D['FJ']['both57_source2']['CA'];prefix=C['all7_pure_CA_tail_prefix57_systems'][0]
 assert prefix['retained_zero_prefix_length']==0 and Q(prefix['unit57_determinant'])==Q(FJ['determinant'])!=0
 term=Terms();binding={old:term.put(('v',new)) for old,new in prefix['bindings']}
 f63=ca['actual_zero_odd_I0_generators'][7][1];assert f63==7790
 Y=FJ['F9_nodes'];original=ca['all48_zero_IFT_quadratic_coefficient_equations']+Y+[f63]
 expected,nold=term.eval(S['DAG_scalar_nodes'],original,binding)
 oldactual,_=term.eval(C['DAG_scalar_nodes'],Y,{})
 actualroots=prefix['all57_unit_equations']+[505938]
 have,nnew=term.eval(C['DAG_scalar_nodes'],actualroots,{})
 equations=[expected[i] for i in ca['all48_zero_IFT_quadratic_coefficient_equations']]+[term.plus(expected[y],term.neg(oldactual[y]) if j>=7 else term.zero) for j,y in enumerate(Y)]
 assert [have[i] for i in prefix['all57_unit_equations']]==equations
 assert have[505938]==expected[f63]
 assert C['actual_C0_basis8_scalar_nodes'][-1]==505938
 telescope=C['all35_regular_coordinate_divisions'][0];assert telescope['original_pure_CA_source_function']==7790 and telescope['restriction_r_s_node']==505938
 assert chain['operator_scalar_nodes'][9]==['actual_C0_scalar_function',505938]
 assert chain['exact_existing_U9_C4']==[[1,0,0,0]]+[[0,i,0,0] for i in range(2,10)]
 assert chain['all48_existing_Q_frame_columns_verified'] and chain['coefficient_C_basis']==['1','tauA','tauB','tauA*tauB']
 assert D['joint']['g_CA_all10_coefficients']==C['cached_actual_CA_gA_all_degree3_coefficients']
 assert all(not Q(c) for mon,c in zip(D['joint']['monomials'],D['joint']['g_CA_all10_coefficients']) if sum(mon)<3)
 # The primitive physical catalog is exactly the original source; accepted
 # exactQ48 frame identities permit bypassing expanded membership quotients.
 catalog=D['physical_catalog']['physical12_coefficient_roots'];assert [r['physical_row'] for r in catalog]==chain['physical_target_rows12']
 for r in catalog:
  j=r['physical_row'];aA,bA=S['families']['CA']['complete_normalized_physical_w_coefficients'][0][j];aB,bB=S['families']['CB']['complete_normalized_physical_w_coefficients'][0][j]
  assert [r[k] for k in ['aA_node','aB_node','bA_node','bB_node']]==[aA,aB,bA,bB]
 assert len(D['FJ']['both57_source2']['CA']['all57_source2jets'])==len(D['FJ']['both57_source2']['CB']['all57_source2jets'])==57
 assert len(D['actual20']['YB_in_YA_quadratic'])==9 and len(D['joint']['all105_unknown_jets'])==105
 assert all(len(p)==10 for i,p in D['joint']['all105_unknown_jets'])
 ep=o/'independent-conormal26-adjoints-primitiveV-and-plane-gA-source-gate.json'
 dump(o,ep.name,{'inputs':inputs,'all100_fixed_source_basis_receipts_and26quotient_verified':True,
  'all1977_scalar_monomial_adjoint_checks_complete':True,'actual_nonzero_patterns':32,'active_orders':counts,
  'all_order4_callbacks_annihilated_after_actual_conormal_projection':True,'all56_literal_chainV_sources_and_requested_weights_verified':True,
  'gA505938_is_original7790_in_prefix0_zero7_unit57_system':True,'all57_prefix_equations_and_clone_function_AST_verified':True,
  'prefix_original_AST_nodes':nold,'prefix_actual_AST_nodes':nnew,'actual_prefix_unit_determinant':prefix['unit57_determinant'],
  'primitive_Q48_V_recipe_per_physical_row':['(bA,-bB,0)','(delta,0,-bB)','(0,delta,bA)','(-nuB*bB,nuA*bA,delta)'],
  'primitive_U8_V_recipe':['7*(-zAi,0,0)','(-gA_plane,0,0)'],
  'allorder_primitive_replacement_hypothesis':'exact existing Q48 multiplication/frame identities independently accepted in pinned chain peer; no quotient value is set zero',
  'source_callback_availability':'mixedquadratics from fullnine own57 source2 and nonlinearCB transport; cubicpure(r,ell) from COMMON105 bivariate3 and gAplane3 under s=r+ell/5',
  'no_independent_ownCB_zero7_copy_or_honest_r_g_higher_coefficient_required_for_primitiveV':True,
  'actual_conormal_V26_by56_values_not_evaluated':True,'new_physical_or_gauge_or_source_coefficients':0})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'conormal26_Hadamard_source_profile_and32patterns_accepted':True,'allorder_gA505938_prefix0_source_and_primitiveV_gate_accepted':True,
  'actual_V_values_or_new_source_coefficients_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
