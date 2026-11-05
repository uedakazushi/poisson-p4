#!/usr/bin/env python3
"""Independent99-source/typedC0 ancestry and separate intrinsic S0 identity."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,fields,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,d):(o/n).write_text(json.dumps(d,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='b22cf94701dc41112609322c47c87d651e28bb7b71c934814153560459415faf'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 d=rd(b+'actual99-selected-source-K-callback-roots-and-literal-ancestry.json')
 put('phase.json',{'phase':'original99_source_and_literal_callback_ancestry','no_formal_gradient_or_second_response':True})
 gp='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'
 z=fld(gp,['all741_zero_origin_source_combinations','old640_typed_closure_nodes','old640_source_columns','old79_atom_ids','native_S0_columns'])
 yp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';sel=fld(yp,['actual_selected_Y0_source_columns'])['actual_selected_Y0_source_columns']
 comb=[z['all741_zero_origin_source_combinations'][j]['original768'] for j in sel]
 assert len(sel)==99 and sel==d['selected99_kernel741_indices'] and comb==d['selected99_original768_constant_combinations']
 cols=[dict(v) for v in z['old640_source_columns']];ops=dict(z['old640_typed_closure_nodes']);known79=set(z['old79_atom_ids'])
 cp='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';P=fld(cp,['honest_D0_relation_columns'])['honest_D0_relation_columns']
 gamma=fld('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json',['Gamma_two'])['Gamma_two'];rows={i for row in gamma for i,x in row};assert rows=={1,2,4,17,18,21,22,28,55} and gamma==d['Gamma_two']
 old=sorted({int(j) for c in comb for j,x in c if int(j)<640});roots=[];seed=set();c0seed=set()
 for j in old:
  rec=[]
  for i in sorted(rows):
   if i<48 and i in cols[j]:rec.append([i,['typed',cols[j][i]]]);seed.add(cols[j][i])
   elif i==55:rec.append([i,['C0',P[j][8]]]);c0seed.add(P[j][8])
  roots.append([j,rec])
 assert len(old)==80 and roots==d['needed_original_old_columns']
 seen={};front=set(seed)
 while front:
  j=front.pop()
  if j in seen:continue
  n=ops[j];seen[j]=n
  if n[0]=='actual_C0_scalar_function':c0seed.add(n[1]);kids=[]
  elif n[0]=='product':kids=n[3]
  elif n[0]=='sum':kids=[v[0] for v in n[1]]
  else:raise AssertionError(('unrecognized_typed_op',j,n))
  for k in kids:
   assert k<j;front.add(k)
 assert sorted(seen.items())==[(i,n) for i,n in d['typed_program']] and len(seen)==50
 rawp='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json'
 known164={i for i,v in fld(rawp,['all_requested164_root_values'])['all_requested164_root_values']}
 wanted=c0seed;assert len(wanted)==27
 # Scan the immutable original scalar AST, retaining ONLY wanted metadata.
 text=(ROOT/cp).read_text();assert len(text.encode())<12*1024*1024;dec=json.JSONDecoder();pos=text.index('[',text.index('"DAG_scalar_nodes"'))+1;literal={}
 for i in range(max(wanted)+1):
  while text[pos].isspace() or text[pos]==',':pos+=1
  node,pos=dec.raw_decode(text,pos)
  if i in wanted:literal[i]=node
 del text
 boundary=[];arith=[]
 for i,n in sorted(literal.items()):
  if i in known79|known164 and n[0]!='c':boundary.append([i,'source164' if i in known164 else 'C0_79_only',n])
  else:arith.append([i,n])
 assert arith==d['literal_C0_program']==[[0,['c',0,1]]] and boundary==d['literal_cache_boundaries']
 assert len(boundary)==26 and sum(k=='source164' for i,k,n in boundary)==24
 assert [i for i,k,n in boundary if k=='C0_79_only']==[299650,327554]
 assert not d['literal_inverse_nodes_retained'] and not d['unrecognized_boundaries_exact_frontier']
 before=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_HIGH3_INPUT_GATE_V2_20261105_v63/617-direct-old-highestthird-coherent-source-positive-input.json',['pruned_original_old_columns_roots'])['pruned_original_old_columns_roots'];oldbefore={int(v[0]) for v in before};new=sorted(set(old)-oldbefore);assert len(new)==18
 # SAME selected99 original768 source coefficients, with ALL128 native columns
 # including the16 scalar free generators, checked intrinsically beforeL0.
 S0=list(map(sp,z['native_S0_columns']));assert len(S0)==128
 prior=0;native_combs=[]
 for c in comb:
  nc={int(j)-640:Q(x) for j,x in c if int(j)>=640};native_combs.append(nc);terms={}
  for j,x in nc.items():
   for i,y in S0[j].items():terms.setdefault(i,[]).append((x,y))
  for ts in terms.values():
   den=math.lcm(*(x.denominator*y.denominator for x,y in ts));U=sum(abs(x.numerator*y.numerator)*(den//(x.denominator*y.denominator)) for x,y in ts);prior=max(prior,den.bit_length(),U.bit_length())
 assert prior<=384
 put('phase2-positive-point-native-prior.json',{'phase':'positive_bound_before_99_intrinsic_nativeS0_checks','positive_prior_bits':prior,'not_a_second_response':True})
 bits=0
 for c in native_combs:
  q={}
  for j,x in c.items():
   for i,y in S0[j].items():
    q[i]=q.get(i,Q(0))+x*y;bits=max(bits,abs(q[i].numerator).bit_length(),q[i].denominator.bit_length())
    if not q[i]:del q[i]
  assert not q
 assert bits<=prior
 # Pure arithmetic/order source only; this does not evaluate sensitivities.
 uses=[j for j,n in seen.items() if n[0]=='product' and 94 in n[3]]
 assert sorted(uses)==[516,602,663,724,785,846,907]
 assert all(seen[j][3]==[94,100+k] for k,j in enumerate(sorted(uses)))
 assert seen[94]==['actual_C0_scalar_function',327554]
 orderpeer='reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63_REVIEW/report.json';order=rd(orderpeer);assert order['all_passed'] and pins[orderpeer]=='dee271e97f4904d80c812f046666a8101fda1805adddffe17e8fe7e59f2d51f1'
 # No coefficient values are imported from the order theorem here.
 cfullpeer='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63_REVIEW/report.json';cf=rd(cfullpeer);assert cf['all_passed'] and pins[cfullpeer]=='76e1aa8842efe0d6d5af78572641c23437ebc902d612c8424876d0058624a48e'
 proof={'selected99_original741_and768_indices_exact':True,'all80_original_old_columns_and50_typed_nodes_literal':True,'18_new_original_old_columns_not_reused_from617_pruning':new,'all27_C0_nodes_scanned_at_original_indices':True,'boundary24_literal164_and_gB299650_h327554':True,'no_inverse_or_unrecognized_opaque_defaults':True,'native99_full_intrinsic_S0_combinations_zero_separately':True,'all128_intrinsic_native_columns_include_free16':True,'native_source_product_prior_bits':prior,'native_source_actual_point_bits':bits,'native_effective_second_formula':'For the99constantcombinations S0*cN=0, [LS]rr=L0Srr+LrSr and [LS]rell=L0Srell+LrSell+LellSr in ordinary coefficients; no factorial2. These are whole99 formulas, never individualpivot replacements. Res4 monic carries are included inside the actual L/S component series.','honest327554_control_flow_only':'Leaf94 occurs solely in seven productnodes with zB100..106. The accepted source order theorem has h0=0 and zB point/r/ell0, so each product has planeorder>=3. This gives a structural route to zero second response; it is not a catalog-based h2=0 assertion or a new numeric response.','Kfull_scope':'CfullXstar1 uses peered76e1 coefficients; selected99 callback seconds and any future fullK values require separate cache/formal source controls and wholejobgate.','no_K_gradient_or_second_response_value_evaluated':True}
 ep='independent99-source-ancestry-intrinsic-nativeS0-and-control-proof.json';put(ep,proof)
 owned=count([r,d,z,sel,comb,cols,ops,P,gamma,rows,old,roots,seed,c0seed,seen,known79,known164,literal,boundary,arith,before,oldbefore,new,S0,native_combs,order,cf,proof,pins])+12000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_ancestry_point_native_identity_only':True,'no_K_or_second_response_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()
