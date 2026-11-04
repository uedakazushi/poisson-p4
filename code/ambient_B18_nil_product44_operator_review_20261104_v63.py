#!/usr/bin/env python3
"""Independent finite-source reader and universal polynomial nil products."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63';STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def add(*ps):
 out={}
 for p in ps:
  for m,c in p.items():
   v=out.get(m,Q(0))+c
   if v:out[m]=v
   else:out.pop(m,None)
 return out
def scale(p,q):return {m:c*q for m,c in p.items() if c*q}
def mul(*ps):
 out={():Q(1)}
 for p in ps:
  dest={}
  for m,c in out.items():
   for n,d in p.items():
    k=tuple(sorted(m+n));v=dest.get(k,Q(0))+c*d
    if v:dest[k]=v
    else:dest.pop(k,None)
  out=dest
 return out
def cmul(a,b,na,nb):
 out=[{} for _ in range(4)]
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   v=mul(x,y)
   if i&j&1:v=mul(v,na)
   if i&j&2:v=mul(v,nb)
   out[i^j]=add(out[i^j],v)
 return out
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 hashes={}
 def load(p):hashes[p]=sha(p);return read(p)
 def get(st,n='report.json'):return load('reports/'+st+'/'+n)
 R=get(AUTHOR);assert R['all_passed'] and hashes['reports/'+AUTHOR+'/report.json']=='d093004d8fdad173825cc4b8eb5d5b70ea77615b2c1dfbc8aaa3d1116105eb4f'
 for p,h in R['input_sha256'].items():assert sha(p)==h;hashes[p]=h
 for n,h in R['evidence_sha256'].items():p='reports/'+AUTHOR+'/'+n;assert sha(p)==h;hashes[p]=h
 radical=get('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json')
 radpeer=get('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW');assert radpeer['all_passed'] and hashes['reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json']=='07b8ec2ea5d6dc419c15fad7cb3ea9526d817037b48c8d96b9966882c1dc851f'
 old=get('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63','actual-UP-QB-640-typed-operator.json')
 C0=get('AMBIENT_B18_T0_HONEST_C0_20261104_v63','honest-C0-D0-D4-DAG.json')
 D=get(AUTHOR,'actual-nil-product44-and-cube-DAG.json')
 assert radical['NilD0_exactly_H1_of_actual_complex'] and radical['actual_X44_by780']['final_shapes']=={'V':[3,44],'X':[44,780]}
 beta=C0['actual_C0_basis8_scalar_nodes'];assert beta==[7148,7199,7454,7505,7688,7943,8045,505938]
 assert C0['actual_C0_free9_images_in_CA_trace_pairs']==[[1,0]]+[[0,i] for i in beta]
 assert old['exact_existing_U9_C4']==[[1,0,0,0]]+[[0,i+2,0,0] for i in range(8)]
 ons=old['operator_scalar_nodes'];assert [ons[i+2] for i in range(8)]==[['actual_C0_scalar_function',i] for i in beta]
 assert ons[old['nu_A_B_operator_nodes'][0]]==['actual_C0_scalar_function',20086]
 qcols=radical['actual_Qrad4_by36_columns'];ns=radical['new_Cscalar_and_R0component_nodes'];assert len(qcols)==36
 blocks=[]
 for g in range(9):
  base=None
  for a in range(4):
   nodes=[ns[i] for i in qcols[4*g+a]]
   assert all(n[0]=='R0_component_of_C' and n[2]==i for i,n in enumerate(nodes))
   root=nodes[0][1];assert all(n[1]==root for n in nodes)
   if a==0:base=root
   else:
    assert ns[root][0]=='C_product' and len(ns[root][1])==2 and base in ns[root][1]
    other=[j for j in ns[root][1] if j!=base];assert len(other)==1 and ns[other[0]]==['C_basis',a]
  blocks.append(base)
 assert blocks==D['actual_Qrad36_native_C_block_roots']
 assert D['actual_Qrad_even36_R0_component_root_ids']==[q[0] for q in qcols]
 V=radical['actual_V3_by44_rows_after_unit1_cancellation']
 assert all(V[k][:36]==[q[k+1] for q in qcols] for k in range(3))
 for i in range(8):
  assert ns[V[0][36+i]][0]=='R0_neg'
  assert ns[ns[V[0][36+i]][1]]==['R0_old_operator',i+2]
  assert ns[V[1][36+i]]==ns[V[2][36+i]]==['C_QQ',0,1]
 assert D['actual_beta8_C0_node_ids']==beta and D['actual_nuA_C0_node_id']==20086
 registry={}
 def variable(key):
  if key not in registry:registry[key]=len(registry)
  return {(registry[key],):Q(1)}
 def c(q):return {():Q(q)} if q else {}
 coeff={('C0',i) for i in beta+[20086]}|{('radical',q[0]) for q in qcols}
 pool=D['all_new_R0_scalar_nodes'];assert len(pool)==492
 roots=set(D['literal_bilinear44_output_roots']+D['literal_one_input_cubic44_output_roots']+[D[k] for k in ['cycle_representative_even_b0','cycle_representative_odd_scalar_B','second_cycle_representative_even','second_cycle_odd_scalar','actual_even_product_on_cycles','actual_odd_product_scalar_on_cycles','actual_even_cube_on_cycles','actual_odd_cube_coefficient_factor']])
 uses=Counter(roots)
 for i,n in enumerate(pool):
  if n[0] in ('sum','product'):assert all(0<=j<i for j in n[1]);uses.update(n[1])
 values={};peak=0
 for i,n in enumerate(pool):
  o=n[0]
  if o=='QQ':v=c(Q(n[1],n[2]))
  elif o=='actual_R0_function':assert tuple(n[1:]) in coeff;v=variable(('actual',n[1],n[2]))
  elif o=='arbitrary_cycle_input':assert n[1] in ('z','zprime') and 0<=n[2]<44;v=variable(('input',n[1],n[2]))
  elif o=='sum':v=add(*(values[j] for j in n[1]))
  else:assert o=='product';v=mul(*(values[j] for j in n[1]))
  values[i]=v;peak=max(peak,sum(map(len,values.values())));assert peak<150000
  if o in ('sum','product'):
   for j in n[1]:
    uses[j]-=1
    if uses[j]==0:values.pop(j,None)
  if uses[i]==0:values.pop(i,None)
 # Independently form products by the four-basis multiplication law, rather
 # than reproducing the producer's explicit four-row mu implementation.
 nu=variable(('actual','C0',20086));nb=variable(('universal','nuB'))
 bet=[variable(('actual','C0',i)) for i in beta];q0=[variable(('actual','radical',q[0])) for q in qcols]
 z=[variable(('input','z',i)) for i in range(44)];zp=[variable(('input','zprime',i)) for i in range(44)]
 a0=add(*(mul(q,x) for q,x in zip(q0,z[:36])));ap0=add(*(mul(q,x) for q,x in zip(q0,zp[:36])))
 B=add(*(mul(q,x) for q,x in zip(bet,z[36:])));Bp=add(*(mul(q,x) for q,x in zip(bet,zp[36:])))
 assert values[D['cycle_representative_even_b0']]==a0 and values[D['cycle_representative_odd_scalar_B']]==B
 assert values[D['second_cycle_representative_even']]==ap0 and values[D['second_cycle_odd_scalar']]==Bp
 pair=[a0,B,{},{}];pairprime=[ap0,Bp,{},{}];pairproduct=cmul(pair,pairprime,nu,nb)
 expected=[]
 for g in range(9):expected.extend(cmul(pair,zp[4*g:4*g+4],nu,nb))
 expected.extend(add(mul(a0,y),mul(ap0,x)) for x,y in zip(z[36:],zp[36:]))
 assert len(expected)==44 and all(values[i]==p for i,p in zip(D['literal_bilinear44_output_roots'],expected))
 assert values[D['actual_even_product_on_cycles']]==pairproduct[0] and values[D['actual_odd_product_scalar_on_cycles']]==pairproduct[1]
 square=cmul(pair,pair,nu,nb);triple=cmul(square,pair,nu,nb);cube=[]
 for g in range(9):cube.extend(cmul(square,z[4*g:4*g+4],nu,nb))
 factor=add(scale(mul(a0,a0),3),mul(nu,B,B));cube.extend(mul(factor,x) for x in z[36:])
 assert all(values[i]==p for i,p in zip(D['literal_one_input_cubic44_output_roots'],cube))
 assert values[D['actual_even_cube_on_cycles']]==triple[0] and values[D['actual_odd_cube_coefficient_factor']]==factor
 assert mul(B,factor)==triple[1]
 # Universal identities with arbitrary off-cycle odd residues.
 g=[variable(('universal','q',i)) for i in range(4)];a=[variable(('universal','a',i)) for i in range(4)]
 b0=variable(('universal','b0'));b=variable(('universal','B'));bprime=variable(('universal','Bp'))
 p=[b0,b,{},{}];qa=cmul(g,a,nu,nb);qprod=cmul(g,cmul(p,a,nu,nb),nu,nb)
 assert qprod==cmul(p,qa,nu,nb)
 odd=add(mul(b0,bprime),mul(qa[0],b));res=[add(qa[1],scale(bprime,-1)),qa[2],qa[3]]
 got=[add(qprod[1],scale(odd,-1)),qprod[2],qprod[3]]
 want=[mul(b0,res[0]),add(mul(b0,res[1]),mul(nu,b,res[2])),add(mul(b,res[1]),mul(b0,res[2]))]
 assert got==want
 assert add(qprod[0],scale(add(mul(b0,qa[0]),mul(nu,b,bprime)),-1))==mul(nu,b,res[0])
 # Characteristic-zero diagonal cubes detect all triple products.
 x,y,t=[variable(('polarization',k)) for k in range(3)]
 cubeof=lambda v:mul(v,v,v)
 lhs=add(cubeof(add(x,y,t)),scale(cubeof(add(x,y)),-1),scale(cubeof(add(x,t)),-1),scale(cubeof(add(y,t)),-1),cubeof(x),cubeof(y),cubeof(t))
 assert lhs==scale(mul(x,y,t),6)
 evidence={'all492_DAG_nodes_and44_bilinear44_cubic_roots_verified':True,'all36_actual_Qrad_Clinear_blocks_and_same_C0_U_beta_nu_bound':True,'off_cycle_three_residual_and_even_correction_identities':True,'all_cycle_even_recovery_and_actual_D0_multiplication':True,'boundary_descent':'Complete H1 isomorphism sends boundary to zero D0 class; products have zero class and hence are complete X boundaries. No new780-column homotopy claimed.','chain_level_commutative_or_associative_free9_ring_assumed':False,'arbitrary_cycle_cube_gate_equivalent_by_polarization':True,'all_cycles_generated_or_cube_zero_nonzero_decided':False,'live_DAG_polynomial_coefficients_peak':peak,'independent_monomial_representation':'sorted variable-index tuples, generic actual coefficient symbols'}
 ev=out/'independent-product44-source-and-universal-proof.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
 for p in [str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'actual_44_bilinear_and44_cubic_operator_source_and_coefficients_verified':True,'off_cycle_and_even_correction_proof_verified':True,'product_descends_to_actual_NilD0_H1':True,'cubic_gate_quantifies_ALL_exact_cycles':True,'nil_generator_list_cube_decision_or780_homotopy_claimed':False,'no_new_source_jet_or_scalar_coefficients':True,'input_sha256':hashes,'evidence_sha256':{ev.name:hashlib.sha256(ev.read_bytes()).hexdigest()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)
if __name__=='__main__':main()
