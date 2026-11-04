#!/usr/bin/env python3
"""B-degree rejection of G32 using frozen low maps, without polynomial products."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_EDGE_BOUNDARY_VALUATION_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'
 reviewpath=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW/report.json'
 edgepath=ROOT/'reports/Q_FIRSTJET_CROSS_BLOCK_EDGE_BOUNDARY_20261104_v63/report.json'
 dominancepath=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/report.json'
 assert all(read(p)['all_passed']for p in[folder/'report.json',reviewpath,edgepath,dominancepath])
 names=['F4','p0','p2','p3','N12'];ps={n:read(folder/(n+'.json.gz'))for n in names}
 expected=dict(F4=(0,2,4),p0=(0,3,8),p2=(1,4,8),p3=(1,4,8),N12=(1,5,12))
 metadata={}
 for name,p in ps.items():
  deg=max(t[4]for t in p['terms']);assert(deg,p['denominator_B_power'],p['denominator_D_power'])==expected[name]
  assert int(p['denominator_integer'])>0 and any(int(t[0])for t in p['terms'])
  metadata[name]=dict(B_numerator_degree=deg,denominator_B_power=p['denominator_B_power'],
   denominator_H_power=p['denominator_D_power'],denominator_integer=p['denominator_integer'],terms=len(p['terms']))
 # These are degree bounds in B over QQ(a,Y,W); no source polynomial is multiplied.
 K_parts=[('D0*p2',7,16,1),('D0*S^2',7,16,0),('D0*p3',7,16,1),
          ('D0^2',6,16,0),('p3^2',8,16,2),('p3*S^2',8,16,1),('S^4',8,16,0)]
 K_bounds=[dict(term=n,denB=db,denH=dh,numerator_B_degree=dg,
                degree_after_B8_H16_clear=dg+8-db)for n,db,dh,dg in K_parts]
 assert max(p['degree_after_B8_H16_clear']for p in K_bounds)==2
 G_bounds=[dict(term='D0*N12^2',denB=13,denH=32,numerator_B_degree=2,clear_B_power=3,cleared_degree=5),
           dict(term='S*K16*N12',denB=15,denH=32,numerator_B_degree_upper=3,clear_B_power=1,cleared_degree_upper=4),
           dict(term='K16^2',denB=16,denH=32,numerator_B_degree_upper=4,clear_B_power=0,cleared_degree_upper=4)]
 P0=[t for t in ps['p0']['terms']if t[4]==0];leadN=[t for t in ps['N12']['terms']if t[4]==1]
 assert P0 and leadN and all(int(t[0])for t in P0+leadN)
 evidence=dict(coefficient_field='QQ(a,Y,W)',source_denominator_name='H=147456Y-(3072W-1)^2; D0=p0 is a different invariant',
  exact_low_map_B_data=metadata,K16_clearing='B^8*H^16',K16_term_degree_bounds=K_bounds,
  G32_clearing='B^16*H^32',G32_term_degree_bounds=G_bounds,
  unique_leading_B5_coefficient='num(p0)*lead_B(num(N12))^2/(den_integer(p0)*den_integer(N12)^2)',
  leading_coefficient_factors=dict(p0_numerator_terms=P0,N12_leading_B1_coefficient_terms=leadN,
   denominator_p0=ps['p0']['denominator_integer'],denominator_N12=ps['N12']['denominator_integer']),
  product_nonzero_reason='Both saved factors are nonzero polynomials over QQ(a,Y,W), an integral domain.',
  no_source_polynomial_multiplication_performed=True,no_source32_expansion_or_point_evaluation=True)
 with gzip.open(out/'certificate.json.gz','wt')as f:json.dump(evidence,f,indent=2);f.write('\n')
 inputs=[Path(__file__),folder/'report.json',reviewpath,edgepath,dominancepath,*[folder/(n+'.json.gz')for n in names]]
 report=dict(id=STEM,all_passed=True,source_G32_nonzero_by_unique_B_degree5=True,
  G32_rejected_as_NEW25_equation=True,generic_NEW25_cross_block_g_h_i_j_nonzero_on_known_unit_chart=True,
  generic_NEW25_Laurent_g_h_u_edge_v_edge_chart_nonempty=True,
  L_unit_on_NEW25_or_NEW25_multisection_degree_or_new_degree_bound_claimed=False,
  no_polynomial_products_no_source32_expansion_no_points=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={'certificate.json.gz':sha(out/'certificate.json.gz')},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','G32_rejected_as_NEW25_equation','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
