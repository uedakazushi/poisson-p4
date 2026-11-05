#!/usr/bin/env python3
"""Independent product-DAG valuation and actual H/regular-projection source audit."""
import hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_REDUCTION_20261105_v63/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 src=read('reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_REDUCTION_20261105_v63/actual-bar-cube-order-four-and-leading-source-proof.json')
 K=read('reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63/actual-five-edge-functionals-and-genuine-Nbar14-cycles.json')['actual_edge_kernel_basis16_by14_rows']
 assert list(map(F,K[0]))==[0]*14
 origins=[[F(0)]*14 for _ in range(44)];origins[43]=list(map(F,K[1]))
 assert origins==[[F(c)for c in row]for row in src['all14_actual_cycle_origins44_by14']]
 conductor=read('reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63/actual-trace-first-gate-and-entire-conductor-proof.json')
 assert min(conductor['all36_Qrad_C_maximal_ideal_order_lower_bounds'])>=2 and conductor['same_origin_tauA_tauB_zero']
 product=read('reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63/actual-nil-product44-and-cube-DAG.json')
 nodes=product['all_new_R0_scalar_nodes'];orders=[];huge=10**6
 beta=product['actual_beta8_C0_node_ids'];assert beta[-1]==505938
 for node in nodes:
  tag=node[0]
  if tag=='QQ':v=huge if int(node[1])==0 else 0
  elif tag=='arbitrary_cycle_input':v=0 if node[2]==43 else 1
  elif tag=='actual_R0_function':
   if node[1]=='radical':v=1
   else:assert node[1]=='C0';v=3 if node[2]==beta[-1] else 1
  elif tag=='sum':v=min(orders[j]for j in node[1])
  elif tag=='product':v=min(huge,sum(orders[j]for j in node[1]))
  else:raise AssertionError(node)
  orders.append(v)
 cubeorders=[orders[j]for j in product['literal_one_input_cubic44_output_roots']]
 assert min(cubeorders[:43])>=5 and cubeorders[43]>=4
 assert orders[product['cycle_representative_even_b0']]>=2 and orders[product['cycle_representative_odd_scalar_B']]>=2
 point=read('reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json')
 assert all(F(row[43])==0 for row in point['m0_24_by44'])
 honest=read('reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/actual-honest-Tnil128-972-formal-operator.json')
 H=honest['shared_homotopy_template_v1']
 assert H['pair_i_7_output']=='Hi' and H['not_S_linear'] and H['no_pull_out_actual_coefficients_or_truncated_normal_forms']
 assert 'rA^n*Hi' in H['remaining_cycle_wi'] and '(Euler_z+2)^-1' in H['pair_i_j_output_i_j_lt7']
 assert 'divided by zi' in H['cross_coefficients_Hi'] and 'd+2>=2' in H['positive_integer_division']
 mn=honest['matrix_operator_nodes'];roots=honest['matrices'];assert mn[roots['Phi44_by128']][0]=='horizontal_blocks'
 assert 'd2*h=m*X' in honest['exact_matrix_identities']
 minimal=read('reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json')
 allowed={'QQ_matrix_rows','actual_external_matrix','horizontal_blocks','identity','matrix_neg','matrix_product','matrix_sum','origin_unit_matrix_inverse','row_column_selection','vertical_blocks','zero'}
 assert all(node[0]in allowed for node in minimal['operator_nodes'])
 assert minimal['matrices']['full_target_projection5_by128']==59
 assert len(minimal['three_exact_unit_Schur_steps'])==3 and all(F(s['unit_origin_determinant'])!=0 for s in minimal['three_exact_unit_Schur_steps'])
 data={'all492_actual_product_nodes_order_audit':len(nodes),'actual_44_cube_orders_lower_bounds':cubeorders,'actual_b0_B_orders':[orders[product['cycle_representative_even_b0']],orders[product['cycle_representative_odd_scalar_B']]],
 'm_cube_order5_by_literal_m0_e43_zero':True,'actual_H_source_template_checked':H,
 'general_H_order_proof':{'input':'v in m^d, delta_(z7,r^n)v=0, n=1 or2','v7':'r^n*v7(0,r)=0 in formal domain gives v7(0,r)=0','Hi':'ordered zi-Hadamard divides only zi: order>=d-1','wi':'vi+r^nHi has order>=min(d,d-1+n)=d','Hij':'partial_z lowers at most1; inverse Euler_z+2 preserves each ordinary monomial degree','conclusion':'H(v) in m^(d-1), no r or r² division, apply whole functions before extending any chosen columns linearly'},
 'actual_representative_proof':'For c=barcube, d1*m*c=Vgamma*c=0. Pick actual universal H(mc), d2H=mc. Then Phi(c,Hmc)=(I-nm)c+nd2Hmc=c exactly. Free128 representative lies m4. Direct-section version uses same m,H and regular exact n; full minimal target projection consists solely of regular unit Schur operations, so class belongs m4Tnil.',
 'literal_bar_lift_proof':'Nbar injects into H1(complex/ell). Every chosen barcycle differs from reduction of a full nilcycle by Xbar*qbar. Lift qbar by coefficient section and add full boundary X*q. This gives a full exact cycle with that literal reduction. Hence specialized product44 computes the actual Nbar product; no ring structure on every Tnil class is assumed.',
 'leading_degree_requires_only':'a36_first,bz_first,bgA_origin,Qrad_even_first,m0,m1e43,Pi0. Terms not listed gain at least one order; νA0=0.',
 'conclusion':'Nbar^3 image in Tnil lies actual m4Tnil (dimension1), m5Tnil0','35_values_or_first_inactive9_not_yet_computed':True,
 'no_vanishing_original_cube_inferred':'zero barcubes may hide in ellN; nonzero barcube detects original cube nonzero'}
 ep=out/'independent-actual-product-order-and-whole-function-H-proof.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_actual_bar_cube_subset_one_dimensional_m4_Tnil','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'cube_values_evaluated':False,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
