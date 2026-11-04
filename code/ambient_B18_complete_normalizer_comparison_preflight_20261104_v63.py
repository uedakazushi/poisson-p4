#!/usr/bin/env python3
"""Finite-block/proof preflight only; no syzygy solve or coefficient expansion."""
import argparse,hashlib,json,math,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_COMPLETE_NORMALIZER_COMPARISON_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    stems=['AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63','AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW','AMBIENT_B18_PROPER_SYZYGY_ADAPTER_PREFLIGHT_20261104_v63','AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW','AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW']
    files=[ROOT/'reports'/s/'report.json' for s in stems];reports=list(map(read,files));assert all(r['all_passed'] for r in reports)
    checked={}
    for r in reports:
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p;checked[p]=h
    primary,peer,adapter,normalizer,dvr=reports
    assert primary['normalized_E4_H0_zero_depth1_proved'] and primary['actual_local_core_CI_2_2_3_CM1_multiplicity12_proved']
    assert adapter['complete_source_syzygy_recipe_count']==32 and adapter['regular_generators8'] and not adapter['physical_free4_coefficients_or128_syzygies_claimed']
    assert dvr['product_base_DVR_special_fiber']==20 and dvr['raw_coker_DVR_min_generators']==26 and dvr['conditional_actual_coker_relation_bound']==326
    np=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json';N=read(np)
    assert N['actual_R0_rows9']==[27,108,109,112,113,117,118,119,61] and N['free4_basis']==['1','tau_CA','tau_CB','tau_CA*tau_CB']
    assert len(N['all12_zero_physical_difference_columns'])==12 and len(N['normalizer_tensor_relation_columns4_by48'])==48 and len(N['actual_C0_free9_to_E4_columns'])==9
    Eoriginal=[4,48,128,224,280,224,112,32,4];Ereduced=[4*choose(8,i) for i in range(9)]
    assert Eoriginal==[4,48,128,*Ereduced[3:]] and Ereduced[:4]==[4,32,112,224]
    Cbounds=[26*choose(9,i) for i in range(10)];assert Cbounds[:4]==[26,234,936,2184]
    # Cancel the23 obvious extra DVR generators when physical3 are included.
    F1bound=234-23;Lgen=936+40;Lrel=2184+F1bound;Krel=Lrel+752
    assert (F1bound,Lgen,Lrel,Krel)==(211,976,2395,3147)
    naive_Dgens=20*Lgen;naive_Drels=20*Krel+14*Lgen
    assert (naive_Dgens,naive_Drels)==(19520,76604)
    recipe={'source_status':{'CI8_and_E4_H0_zero':'accepted','physical_T12_Cfree4_interface':'MS separate typed-interface/Res4 producer, not computed in this packet or inferred from source type','actual_DVR20_basis_actions_and26by326_entries':'not yet serialized by prior preflight'},
      'actual_normalizer_basis':N['free4_basis'],'R0_physical_rows':N['actual_R0_rows9'],
      'source_adapter_orientation':'Hnew=T Hold as columns; Qnew=Qold Res(T^t). Res(T^t) hasblock(j,i)=mu(T[i,j]), not transpose(mu(T)).',
      'complete_Q_syzygies':{'over_normalizer_C':'transported Koszul28 plus4 free redundant-source relations','R0_d2_shape':[48,128],'R0_d3_shape':[128,224],'E_resolution_sizes_nonminimal':Eoriginal,'after_redundant16_source_columns_removed':Ereduced,'no_d_inverse':'cancel the nonzero common d in the regular domain, never call it a unit'},
      'literal_low_comparison':{'P':[9,640],'Q_reduced':[4,32],'U':[4,9],'B_lift':[32,640],'identity':'U P=Q B','after_eliminating_unit1_shapes':[[40,752],[3,40]],'H1':'K','H0':'actual coker M','higher_homology':'unwanted; truncated H2 is not K[ell]'},
      'complete_DVR_M_input':{'D':'QQ[[ell]], actualell=5(s_CA-r_CA)','R':'D[[x1,...,x8]]','generators_max':26,'presentation_relations_max':326,'image_basis_rank_max':26,'required_two_inclusions':['P_D=A B','A=P_D H'],'required_injectivity':'nonzero r-minor ofA, or exact DVR image-basis proof','actions':['T_i A=A S_i','[T_i,T_j]=A H_ij'],'consequences':['[S_i,S_j]=H_ij A','cyclicsum(H_jk T_i-S_i H_jk)=0 byJacobi+injectivityA']},
      'complete_R0_M_resolution':{'modules':'F_n=(R^m tensor exterior^n D8) directsum (R^r tensor exterior^(n-1) D8)','worst_sizes':Cbounds,'differential':'A oninternal1; sum(x_i-T_i)iota_i oninternal0; minus sum(x_i-S_i)iota_i oninternal1; minus sum H_ij iota_i iota_j frominternal0 to1','exactness':'filter byexterior degree: A injective leavesM; commuting x_i-actions give completed monicKoszul resolution. Actual actions topologicallynilpotent mod ell.','no_automatic_higher_homotopies':'the2-row amplitude killsH*H; Jacobi forces theonlycubic identity'},
      'Schanuel_complete_L':{'existing_M_shape':[3,40],'free_resolution_of_M':'G3->G2->G1->R3->M','required_lifts':['d1 m=V','V n=d1','d2 h=I-m n'],'L_generators':'S=[n d2, I-n m]','L_relations':'[(d3,0),(-h,n)]','completeness':'ifS(w,x)=0,set y=m x-d2w; thenx=n y andd2(w+h y)=0','L_generator_upper':Lgen,'L_relation_upper':Lrel},
      'complete_K_presentation':{'input_cycles':752,'lift_formula':'for v inkerV choose w via completeC Koszul/Weierstrass contraction with d2w=m v, thenlift v as(w,v)','generator_upper':Lgen,'relation_upper':Krel,'no_SyzP_required':True,'source_ancestry_required':'all752cycles andtheirlifts checked; oldr/gstayexistingregularR0functions, nofresh0/0 orinventedderivatives'},
      'DVR_torsion_finish':{'after_complete_K_cokernel_presentation_only':'tensor byH=R/(IAIB), thenbyD/ell; right exactness is enough becauseIAIB K=0','naive_upper_D_generators':naive_Dgens,'naive_upper_D_relations':naive_Drels,'naive_dense_query_not_authorized':True,'test':'H0(K)!=0 iff dimQQ K/ellK >dimFracD K[1/ell]','torsion_cyclic_count':'difference equals number, not length, ofDVRtorsioncyclicsummands','H0D0_equals_H0K':'follows fromH0E4=0; fullDdepth is not computed'},
      'derived_special_fiber_warning':'H2=K[ell] holds only for a genuine comparison modelwithhomologyKdegree1,Mdegree0,nootherhomology. Not fromthe752->40->3truncation.'}
    ep=out/'complete-syzygy-and-DVR-resolution-operator-preflight.json';ep.write_text(json.dumps(recipe,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',np,*files]
    report={'id':STEM,'all_passed':True,'frozen':True,'proof_and_operator_size_preflight_only':True,'CI8_supplies_complete_normalizer_source_syzygy_recipe':True,'DVR_coker_resolution_and_complete_K_cokernel_algorithm_proved':True,'four_R0_coefficient_extraction_or_DVR_basis_actions_computed_here':False,'no_truncated_specialization_kernel_or_Tor_shortcut':True,'actual_K_torsion_or_D0_depth_computed':False,'naive_large_presentation_not_approved_for_materialization':True,'no_new_solve_GB_physical_jets_or_oldcofactor_derivatives':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'checked_prior_source_hashes':checked,'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
