#!/usr/bin/env python3
"""Formal source-bound exact Loewy/length corollary after final witness peer."""
import argparse,json,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63'
def main(out):
    start=time.monotonic();out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'witness':'reports/AMBIENT_B18_T_COMBINED619_QUARTIC_CLOSURE_20261105_v63/actual-nonzero-quartic-witness.json',
      'witness_primary':'reports/AMBIENT_B18_T_COMBINED619_QUARTIC_CLOSURE_20261105_v63/report.json',
      'witness_peer':'reports/AMBIENT_B18_T_COMBINED619_QUARTIC_WITNESS_20261105_v63_REVIEW/report.json',
      'initial40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
      'initial40_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
      'HF3':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/report.json',
      'HF3_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
      'complete5':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
      'complete5_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
      'rank4':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json',
      'Nfinite':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63/report.json',
      'Nfinite_peer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
      'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/report.json',
      'natural_peer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};D={k:read(ROOT/p) for k,p in paths.items()}
    assert all(D[k]['all_passed'] for k in ['witness_primary','witness_peer','initial40_peer','HF3','HF3_peer','complete5_peer','rank4','Nfinite','Nfinite_peer','natural','natural_peer'])
    assert sha(ROOT/paths['witness_peer'])=='a73a13d105b05d5dd24af5f22ce069a340d5ad1ba000a3bc9e6059258c571c32'
    w=D['witness'];initial=D['initial40'];assert w['kernel619_index']==580 and F(w['actual_corrected_pure_quartic_if_bound'])==F(40553074375,178412544)
    assert w['honest_second_linear_and_quadratic_terms_cancel'] and w['all_source_inputs_fixed']
    assert initial['standard_degree0_count']==5 and initial['standard_degree1']==['e3*r','e4*z5','e4*r']
    assert initial['standard_degree2']==['e4*z5*r','e4*r²'] and initial['standard_degree_at_least3']==['e4*r^d']
    # The literal pure-last-target functional kills every multiple of40gens.
    for g in initial['homogeneous_module_generators40']:
        assert all(not(t==4 and not any(m[:7]) and F(c)) for t,m,c in g)
    assert D['HF3']['actual_HF_through_degree3']==[5,3,2,1]
    assert D['complete5']['actual_minimal_generator_presentation_shape']==[5,825]
    source_rows=D['complete5']['actual_free_five_source_rows'];assert source_rows==[22,23,29,30,55]
    basis=[{'target':i,'monomial_original8':[0]*8,'degree':0} for i in range(5)]
    for target,variables in [(3,[7]),(4,[5]),(4,[7]),(4,[5,7]),(4,[7,7]),(4,[7,7,7])]:
        m=[0]*8
        for v in variables:m[v]+=1
        basis.append({'target':target,'source_V56_row':source_rows[target],'monomial_original8':m,'degree':len(variables)})
    hf=[sum(b['degree']==i for b in basis) for i in range(4)];assert hf==[5,3,2,1] and len(basis)==11
    r=D['rank4'];assert r['rank_D_K']==4 and r['dim_k_C_ell_torsion']==4
    muK=11-r['dim_k_C_ell_torsion'];tors=muK-r['rank_D_K'];assert (muK,tors)==(7,3)
    n=D['Nfinite'];assert n['actual_NilD0_DVR_finite_rank']==11
    assert D['natural']['actual_cokernel_length']==5
    proof={
      'actual_base':'Rbar=QQ[[zA0,...,zA6,rA]]=R0/(ell), ell=5(sA-rA)',
      'actual_module':'complete actual five-generator825-relation comparison T',
      'known_initial_module40':'all40 have actual source lifts; accepted Groebner proof gives onlye4*r^d standard in degree>=3',
      'pure_quartic_functional':'literal coefficient ofe4*rA^4; annihilates allknown40 initialgenerators andalltheirdegree4multiples',
      'new_relation':'actual correctedcase580 startsdegree4 andhaspurelastcoefficient40553074375/178412544!=0',
      'degree4_consequence':'known initial quotientdegree4 is1-dimensional, soactualgr4T=0',
      'Nakayama':'m4T=m5T=m*(m4T); m4T is finitelygenerated overNoetherianlocalRbar, hence m4T=0',
      'lower_exactness':'independentactualHF0..3=(5,3,2,1)',
      'exact_length':11,'exact_Loewy_HF':[5,3,2,1],
      'dimensions_m_power_T':[11,6,3,1,0],
      'actual_filtered_monomial_basis11':basis,
      'basis_proof':'theirleadingclasses areexactgradedbases andm4T=0; filteredlifting provesbothspanning andlinearindependence',
      'm3T':'one-dimensional k span ofclass(source55*rA^3), nonzero andannihilatedbyallm',
      'maximal_ideal_Loewy_length_T':4,
      'full_socle_dimension_not_computed':True,
      'actual_multiplication_matrices_not_claimed_from_graded_data':True,
      'rA_is_NOT_the_DVR_parameter':'DVRparameter isell; do notidentifythese',
      'derived_DVR_values':{'rank_D_K':4,'mu_D_K':muK,'mu_D_torsion_K':tors,'rank_D_N':11,'mu_D_N':muK+7,'mu_D_torsion_N':tors,'length_Tnil':11+5},
      'DVR_derivation':'acceptedexactsequence0->K/ellK->T->C[ell]->0 withdimC[ell]=4; NnoncanonicallyK+free7; naturalT->Tnilinjectivewithcokerlength5',
      'not_concluded':['exactH0length','DVRtorsionexponents','originalnilradicalcube','nilradicalnilpotencyindex','all619sourcevalues','fullTnilHilbertfunction'],
      'no_new_coefficient_rank_GB_or_matrix_job':True}
    ep=out/'exact-Loewy11-and-source-basis-corollary.json';ep.write_text(json.dumps(proof,separators=(',',':'))+'\n')
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'status':'formal_corollary_from_FINAL_independent_source_peers','actual_T_length':11,
      'actual_T_Loewy_HF':[5,3,2,1],'actual_T_Loewy_length':4,'actual_filtered_basis_size':11,
      'actual_mu_D_K':7,'actual_mu_D_torsion_K':3,'actual_mu_D_N':14,'actual_mu_D_torsion_N':3,'actual_Tnil_length':16,
      'exact_H0_length_torsion_exponents_or_nil_cube_computed':False,
      'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','actual_T_length','actual_T_Loewy_HF','actual_mu_D_K','actual_mu_D_N','actual_Tnil_length']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
