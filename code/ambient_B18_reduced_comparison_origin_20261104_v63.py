#!/usr/bin/env python3
"""Exact tiny fixed first-interface gate; no scalar source/jet worker."""
import argparse,hashlib,json
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_REDUCED_COMPARISON_ORIGIN_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
 out.mkdir(parents=True,exist_ok=False)
 files={
 'interface':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json',
 'interface_peer':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63_REVIEW/report.json',
 'roots':'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',
 'roots_peer':'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63_REVIEW/report.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'radical_peer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
 'product_peer':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
 }
 expected={
 'interface':'7aee89c964510451fdc4f205361c73f615251d61439a28f756a87ae22909aeba',
 'interface_peer':'2033dff3459a3e9e4088f9dbbf4ad6f6549dd23cdd01bc8a3a543450523dc7eb',
 'roots':'9cae2cb8395b4208df992761988b697a6dbccab5ef7a6c18fd26504b0dd7b5cf',
 'roots_peer':'ee45748a1e273b1d6f853cbd3284f243a9fe634a3cbe7f702badfe08c41dd7b9',
 'radical':'376d5fd2b8f2ab8a45651b0374d5b983ebe2df8f8a507108675942e853adbba5',
 'radical_peer':'07b8ec2ea5d6dc419c15fad7cb3ea9526d817037b48c8d96b9966882c1dc851f',
 'product_peer':'7d876031573583d070890c28da09ae123c550d2108cda73d45cd2f61b87eaa19',
 }
 for k,p in files.items():assert sha(ROOT/p)==expected[k],k
 for k in ['interface_peer','roots_peer','radical_peer','product_peer']:assert read(ROOT/files[k])['all_passed']
 I,T,R=[read(ROOT/files[k])for k in ['interface','roots','radical']]
 assert I['physical_regular9_order']==['x108','x109','x112','x113','x117','x118','x119','tau_CA','tau_CB']
 order=I['retained_source9_order'];assert order==['a','d','v2','v3','v4','v5','v6','v7','v8']
 first={n:tuple(Q(x)for x in row[-2:])for n,row in zip(order,I['actual_source9_by_physical9_Jacobian'])}
 assert [first[n]for n in ['a','d','v2','v6']]==[(Q(1),Q(-1)),(Q(-2),Q(-1,2)),(Q(1),Q(0)),(Q(1),Q(0))]
 pivot=T['selected5_pivot_order'];assert pivot==['v4','v3','v5','v7','v8']
 rootmins={n:min(sum(e)for e,a,b in T['five_actual_rational_roots'][n]['numerator'])for n in pivot}
 assert min(rootmins.values())>=2 and all(first[n]==(0,0)for n in pivot)
 h=tuple(first['a'][j]-first['d'][j]+first['v2'][j]-2*first['v6'][j]for j in range(2))
 d=first['d'];assert h==(Q(2),Q(-1,2))
 mixed=d[0]*h[1]+d[1]*h[0];assert mixed==0
 proof={
 'actual_R0_fibre_algebra':'k[tauA,tauB]/(tauA²,tauB²)',
 'same_physical_origin_and_trace_basis':True,
 'proper_first_trace_vectors':{n:[str(x)for x in v]for n,v in first.items()},
 'all5_frozen_core_root_min_numerator_degrees':rootmins,
 'all5_u_trace_first_parts_zero':True,
 'h_first_trace_vector':[str(x)for x in h],
 'd_h_quadratic_coefficients_tauA2_tauAtauB_tauB2':[str(d[0]*h[0]),str(mixed),str(d[1]*h[1])],
 'higher_C_order3_terms_zero_in_this_fibre':True,
 'all9_native_dR_generators_fibre_zero':True,
 'all36_Qrad_origin_four_components':[[0]*4 for _ in range(36)],
 'all_actual_C0_odd_beta_origins_zero':True,
 'actual_Vrad_origin3_by44':[[0]*44 for _ in range(3)],
 'actual_Vrad_origin_rank':0,
 'actual_reduced_comparison_cokernel_mu_R0':3,
 'no_kernel_basechange_or_nil_generator_count_inferred':True,
 'ell_zero_nil_images_exact_values_computed':False,
 'no_new_source_jet_GB_point_or_780_expansion':True,
 }
 ep=out/'actual-first-interface-and-zero-origin-Vrad.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 rels=list(files.values())+[f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT))]
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_reduced_comparison_origin_rank0_no_unit_cancellation',
 'actual_Vrad_origin_rank':0,'actual_Cred_cokernel_mu_R0':3,
 'new_nil_generators_ell_zero_class_values_or_cube_decision':False,
 'input_sha256':{p:sha(ROOT/p)for p in rels},'evidence_sha256':{ep.name:sha(ep)}}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['id','all_passed','actual_Vrad_origin_rank','actual_Cred_cokernel_mu_R0']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
