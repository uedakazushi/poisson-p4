#!/usr/bin/env python3
"""Independent five existing products and exact natural-diagram dimension proof."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_TOLD_GENUINE_G3_VALUES_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return list(map(enc,x))
 return x
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_TOLD_GENUINE_G3_VALUES_20261105_v63/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 sr=read('reports/AMBIENT_B18_TOLD_GENUINE_G3_SOURCE_20261105_v63_REVIEW/report.json');assert sr['all_passed']and not sr['witness_paired']
 V=read('reports/AMBIENT_B18_TOLD_GENUINE_G3_VALUES_20261105_v63/actual-five-cached-products-and-normalizer-image-rank-four.json')
 N=read('reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json')
 M=[[Q(v)for v in col]for col in N['actual_origin_data']['actual_origin_map_columns5']];assert all(col[0]==0 for col in M)
 E=read('reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63/actual-five-edge-functionals-and-genuine-Nbar14-cycles.json')
 a,*b=list(map(Q,E['five_actual_cached_functional_values_a_b1_b2_b3_b4']));xi=Q(E['xi_equals_minus11a_by16']);assert xi==-Q(11,16)*a
 J=read('reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63/actual-head4-z5-map-and-F2-kernel-coefficient.json')
 Jr=Q(J['actual_projected_natural_head4_first_target0_r']);assert Jr==Q(-43373,15129)
 terms=[b[i]*M[4][i+1]for i in range(4)]+[xi*Jr];w=sum(terms,Q(0));assert w==Q(2991593,10086000)and w
 assert terms[:4]==list(map(Q,V['all_four_origin_products_bj_Mj4']))and terms[4]==Q(V['xi_times_Jr'])
 assert V['actual_old_head4_genuine_G3_image']==['0',str(w)]and V['canonical_head0_beta_not_used']
 old=read('reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json');assert old['actual_T_length']==11 and old['actual_Tnil_length']==16
 finite=read('reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/actual-ell-line-values-and-nil-DVR-proof.json');assert finite['image_NilD0_in_NilE4_free_DVR_rank']==7 and finite['normalized_nil_DVR_finite_free_rank']==9
 assert N['exact_sequence_proof']['sequence']=='0→Told→Tnil→NE/ellNE→C/ellC→Cprime/ellCprime→0'
 assert E['genuine_Nbar_dimension']==14
 oldedge=1;intersection=11-oldedge;normalimage=14-intersection;Irank=7;NErank=9;Lrank=NErank-Irank;tormu=Irank-normalimage;Lbar=NErank-normalimage
 assert (intersection,normalimage,Lrank,tormu,Lbar)==(10,4,2,3,5)
 assert V['actual_intersection_Told_Nbar_length']==intersection and V['actual_normalized_nil_image_Nbar_dimension']==normalimage and V['actual_L_DVR_rank']==Lrank and V['actual_L_DVR_torsion_min_generator_count']==tormu
 assert r['actual_edge_rank']==1 and r['actual_L_generic_DVR_rank']==2 and r['actual_L_torsion_mu']==3
 data={'independent_five_actual_cache_products':terms,'actual_oldhead4_edge_Q1_coefficient':w,'actual_oldhead4_K_coefficient':0,'canonical_beta_not_computed_or_used':True,
 'actual_old_edge_rank_k':1,'actual_Nbar_intersection_Told_length':intersection,'actual_normalizer_nil_fiber_image_dimension':normalimage,
 'I_definition':'R0/D-submodule image(N→NE), equivalently N/K; not asserted to be an E4 ideal','actual_I_Dfree_rank':Irank,'actual_NE_Dfree_rank':NErank,'L_definition':'R0/D-module coker(I→NE), not an algebra quotient claim','actual_L_generic_Drank':Lrank,'actual_L_torsion_minimal_Dgenerators':tormu,'actual_Lbar_dimension':Lbar,
 'exact_source_dimension_proof':'epsilonold image⊂kQ1 and one positive witness⇒rank1. Exact compatible sameD0 triangle hasker(theta)=Told, so Nbar∩Told length11-1=10, theta(Nbar) dimension14-10=4. Nbar surjects Ibar, giving rank(Ibar→NEbar)=4. Tor_D exact0L[ell]→k7→k9→Lbar0 gives torsion-generator dimension3 andLbar5; fraction-field rankL=9-7=2.',
 'Q_vs_Lbar':'Q=Tnil/Told=im(theta) is a canonical length5 submodule ofNEbar. Lbar=NEbar/im(Ibar) is a canonical length5 quotient. They are not identified and no S-splitting follows from r-Jordan311.',
 'no_normalizer_image_basis_or_full_L_actions_computed':True,'L_torsion_lengths_and_exponents':'OPEN','original_H0_lengths_exponents_cube_and_exact_nilindex':'not decided by this calculation','no_new_source_gauge_jet_or_derivatives':True}
 ep=out/'independent-five-products-and-exact-normalizer-image-diagram.json';ep.write_text(json.dumps(enc(data),indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_independent_Told_edge_rank1_normalizer_image4_and_L_torsion_mu3','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'actual_edge_rank':1,'actual_normalizer_nil_fiber_image':4,'actual_L_Drank':2,'actual_L_torsion_mu':3,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
