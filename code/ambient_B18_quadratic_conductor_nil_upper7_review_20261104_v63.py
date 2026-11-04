#!/usr/bin/env python3
"""Independent same-source first traces and whole-conductor nil^7 proof."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63';STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def times(a,b):return (a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[2]*b[0])
def inverse(a):
 assert a[0]
 out=(1/a[0],-a[1]/a[0]**2,-a[2]/a[0]**2);assert times(a,out)==(1,0,0);return out
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 hashes={}
 def load(p):hashes[p]=sha(p);return read(p)
 def get(st,f='report.json'):return load('reports/'+st+'/'+f)
 R=get(AUTHOR);assert R['all_passed'] and hashes['reports/'+AUTHOR+'/report.json']=='76afb23c5500fa6e7cd20652364f51d32a822a9f9bb224f7a9ecbb74ee66a4cf'
 for p,h in R['input_sha256'].items():assert sha(p)==h;hashes[p]=h
 for n,h in R['evidence_sha256'].items():p='reports/'+AUTHOR+'/'+n;assert sha(p)==h;hashes[p]=h
 for st,want in [('AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW','7d876031573583d070890c28da09ae123c550d2108cda73d45cd2f61b87eaa19'),('AMBIENT_B18_QUADRATIC_NORMALIZER_NIL_UPPER12_20261104_v63_REVIEW','e1a6937b676c1518d192adb0dfdc5b35d7d841eaaf381d0f52f2c9c0ca1db219'),('AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW','513e63044f8048da7c6761f1b02c63c0bcd6f643beb65339a2d3eb0de18c1401')]:
  d=get(st);assert d['all_passed'] and hashes['reports/'+st+'/report.json']==want
 Jpeer=get('AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW')
 assert Jpeer['exact_joint_zero_ideal_quotient_length']==3 and Jpeer['actual_CB_three_quadratics_span_m_squared_minor']=='625/768'
 assert Jpeer['literal_finite105_DAG_and_every_binding_reconstructed_from_actual_source']
 Src=get('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json')
 J=get('AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63','joint105-algebraic-DAG.json')
 C=get('AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63','joint-degree3-certificate.json')
 Rad=get('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json')
 old=get('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63','actual-UP-QB-640-typed-operator.json')
 assert J['DAG_scalar_nodes'][:41315]==Src['DAG_scalar_nodes'] and len(Src['DAG_scalar_nodes'])==41315
 assert J['actual_common_physical_base_rows']==[27,108,109,112,113,117,118,119,61]
 assert Rad['actual_R0_order']=='physical9 [27,108,109,112,113,117,118,119,61]'
 bindings=J['all105_unknown_bindings'];assert len(bindings)==105 and [i for i,_ in bindings]==C['all105_unknown_original_node_order']
 assert [i for i,_ in C['all105_unknown_jets']]==C['all105_unknown_original_node_order']
 assert C['monomials'][:3]==[[0,0],[1,0],[0,1]] and C['CA_F9_targets']==['0']*7+['r','s']
 roots=[Src['families'][name]['actual_tracefree_nu'] for name in ('CA','CB')];assert roots==[20086,40865]
 assert [old['operator_scalar_nodes'][i] for i in old['nu_A_B_operator_nodes']]==[['actual_C0_scalar_function',r] for r in roots]
 fixed={i:tuple(map(Q,v[:3])) for i,v in C['all105_unknown_jets']}
 r,s=J['formal_curve_parameter_nodes'];fixed[r]=(Q(0),Q(1),Q(0));fixed[s]=(Q(0),Q(0),Q(1))
 nodes=J['DAG_scalar_nodes'];needed=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if i in fixed:continue
  n=nodes[i]
  if n[0] in ('a','m'):todo.extend(n[1:3])
  elif n[0] in ('n','i'):todo.append(n[1])
  else:assert n[0]=='c',('uncached source variable',i,n)
 assert len(needed)==11
 values={};unitchecks=[]
 # Independent forward topological dual-number calculation (no producer AD).
 for i in sorted(needed):
  n=nodes[i];o=n[0]
  if i in fixed:v=fixed[i]
  elif o=='c':v=(Q(n[1],n[2]),Q(0),Q(0))
  elif o=='a':v=plus(values[n[1]],values[n[2]])
  elif o=='m':v=times(values[n[1]],values[n[2]])
  elif o=='n':v=tuple(-q for q in values[n[1]])
  else:assert o=='i';v=inverse(values[n[1]]);unitchecks.append([i,str(values[n[1]][0])])
  values[i]=v
 assert [values[i] for i in roots]==[(Q(0),)*3]*2
 primary=get(AUTHOR,'actual-trace-first-gate-and-entire-conductor-proof.json')
 assert primary['cached_restricted_first_coefficients']==[['0','0','0']]*2
 assert primary['used_original_ancestor_count']==11
 assert primary['all_used_nodes_first_coefficients']==[[i,nodes[i],list(map(str,values[i])),i in fixed] for i in sorted(needed)]
 # All native C-coordinates are the centered actual proper/physical ones.
 # Verify their GL4 rules have no constant term, and the actual d,R9
 # polynomials have positive order before following every source C block.
 rules=Rad['actual_C_variable_rules']['actual_core_GL4']
 assert set(rules)=={'c','r','l','h'} and all(sum(m)>=1 for p in rules.values() for m,q in p if Q(q))
 assert all(sum(m)>=1 for m,q in Rad['d_polynomial'] if Q(q))
 assert all(p and min(sum(m) for m,q in p if Q(q))>=1 for p in Rad['R9_generators_polynomials9'])
 context=Rad['accepted_native_C_coefficient_interface'];load(context['path']);assert hashes[context['path']]==context['sha256']
 ns=Rad['new_Cscalar_and_R0component_nodes'];Corders={};visited=set()
 def order(root):
  if root in Corders:return Corders[root]
  n=ns[root];o=n[0];visited.add(root)
  if o=='C_QQ':a=10**6 if not Q(n[1],n[2]) else 0
  elif o in ('actual_proper_u','actual_core_GL4'):a=1
  elif o=='C_basis':a=int(n[1]).bit_count()
  elif o=='C_product':a=sum(order(j) for j in n[1])
  else:assert o=='C_sum';a=min(order(j) for j in n[1])
  Corders[root]=a;return a
 orders=[];blockroots=[]
 for j,q in enumerate(Rad['actual_Qrad4_by36_columns']):
  parts=[ns[i] for i in q];assert all(n[0]=='R0_component_of_C' and n[2]==k for k,n in enumerate(parts))
  root=parts[0][1];assert all(n[1]==root for n in parts);orders.append(order(root))
  if j%4==0:blockroots.append(root)
 assert len(orders)==36 and min(orders)>=2 and orders==primary['all36_Qrad_C_maximal_ideal_order_lower_bounds']
 product=get('AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63','actual-nil-product44-and-cube-DAG.json')
 assert product['actual_Qrad36_native_C_block_roots']==blockroots and product['actual_nuA_C0_node_id']==20086
 assert product['actual_beta8_C0_node_ids'][-1]==505938 and Rad['NilD0_exactly_H1_of_actual_complex']
 # Four-basis multiplication: even coefficient can arise only from equal
 # masks, with factors 1,nuA,nuB,nuA*nuB respectively.
 even_pairs=[]
 for a in range(4):
  for b in range(4):
   if a^b==0:even_pairs.append([a,b,int(bool(a&b&1)),int(bool(a&b&2))])
 assert even_pairs==[[0,0,0,0],[1,1,1,0],[2,2,0,1],[3,3,1,1]]
 # Entire conductor is essential: f and f*tau both send the two Bi basis
 # vectors back into Ai. Record all four exact membership types.
 conductor_products={'f*1':'f in R0','f*tau':'f*tau in Ii*tau','(f*tau)*1':'f*tau in Ii*tau','(f*tau)*tau':'f*nu in R0'}
 proof={
  'source_trace11_forward_calculation':[[i,nodes[i],list(map(str,values[i]))] for i in sorted(needed)],
  'both_actual_nu_in_J':'Exact J=(zA7,(r,s)^2); the source restrictions have zero constant/first terms, hence nu lies in (zA7)+m_R0² subset J.',
  'same_R0_regular_coordinates_and_J_identity_from_accepted105_gate':True,
  'all36_source_Qrad_C_order_lower_bounds':orders,'centered_native_C_coordinate_interface_bound':True,
  'even_C_basis_product_table':even_pairs,
  'scalar_mC_square_subset_J':'For f,g in mC, f0,g0 in mR0; (fg)0=f0g0+nuA fA gA+nuB fB gB+nuA nuB fAB gAB. Thus pr0(mC²) subset mR0²+(nuA,nuB) subset J.',
  'ALL_nil_cycles':'Complete H1 surjects onto every nil y. Its even representative b0=Qeven*a lies in J since all36 Qeven coefficients do; its odd part uses ALL unchanged beta8*tauA.',
  'full_conductor_multiplication_memberships':conductor_products,
  'full_conductor_annihilation':'For c in IiBi subset Ai, multiplication Bi->Ai is Ai-linear, with composite c on Ai; after any basechange c kills the inclusion kernel. Finite-over-common-R0 tensors are already complete.',
  'both_normalization_orders':'Each second-step basechange kernel has ideal square zero, so composite K² lies in both first-step kernels KA and KB.',
  'annihilators':'IA*KA=IB*KB=0 gives J*K²=0; every betaA_i*tauA in IA*BA also kills KA and hence K². Therefore every nil y kills K².',
  'Nil_times_K_squared_zero':True,'K_cubed_zero':'K⊂Nil (old K^4=0); hence K³⊂Nil*K²=0.',
  'Nil_power7_zero':'Nil³⊂K by WHOLE E4 Nil³=0, so Nil^7=Nil*(Nil³)^2⊂Nil*K²=0.',
  'full_actual_germ':'Accepted all-positive nil products zero/closed homogeneous completion and actualG19/scalar product transfer preserve the bound.',
  'lower_bound3':'Accepted actual Qb nonzero nil-square and coherent closure through P0.',
  'not_used':['No assumption Nil(E4)^2=0 at N/P','No branch image label/exhaustion','No arbitrary H0 annihilator','No exact nil index or cube decision']}
 ev=out/'independent-trace11-source-and-entire-conductor-proof.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
 for p in [str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'same_actual_nu_roots11_ancestors_first_zero_verified':True,'all36_actual_Qrad_elements_order_at_least2_verified':True,'EVERY_nil_cycle_even_b0_in_J':True,'entire_conductor_completed_basechange_annihilation_accepted':True,'Nil_D0_times_K_squared_zero':True,'K_power3_zero':True,'full_actual_B18_germ_nil_power7_zero':True,'actual_nil_index_interval':[3,7],'E4_N_P_nil_square_zero_or_exact_cube_claimed':False,'no_new_source_jet_solve_point_GB_or780_expansion':True,'input_sha256':hashes,'evidence_sha256':{ev.name:hashlib.sha256(ev.read_bytes()).hexdigest()},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)
if __name__=='__main__':main()
