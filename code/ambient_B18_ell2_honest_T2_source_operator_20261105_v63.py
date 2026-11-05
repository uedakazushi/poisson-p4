#!/usr/bin/env python3
"""Exact factored allorder honest T2 source operator; no scalar queries."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();out=out.resolve();out.mkdir(parents=True,exist_ok=False);hs={}
    def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    sourcepath='reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63/phase1-actual-X768-source.json'
    old=read(sourcepath)
    original_sources={k:old['source_hashes'][k] for k in ['B','B_peer','P','SyzQ','SyzQ_peer']}
    for value in original_sources.values():
        assert sha(ROOT/value['path'])==value['sha256'];hs[value['path']]=value['sha256']
    conpath='reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/actual-V-conormal26-by56-unit-minor-and-source-factorization.json'
    con=read(conpath);assert con['rank']==26
    peer=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json');assert peer['all_passed']
    rp='reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json'
    res=read(rp)
    betti=read('reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json')
    adj=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json')
    c2path='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
    c2=read(c2path);assert c2['actual_C2_dimension']==14
    assert res['after11_unit_cancellation']==[3,115,518,1302]
    assert betti['rank']==89 and betti['actual_beta1_C2_over_B']==26
    manifest=old['source_manifest'];assert manifest['actual_X_shape']==[56,768] and manifest['actual_VX_zero']
    C=con['fixed_selected26_columns'];assert len(set(C))==26;Cfree=[i for i in range(56) if i not in C];assert len(Cfree)==30
    sel89=betti['fixed_original_minor_columns'];free429=[i for i in range(518) if i not in sel89];assert len(free429)==429
    free26=adj['conormal_standard_source_indices26'];assert len(free26)==26
    nodes=[]
    def put(n):nodes.append(n);return len(nodes)-1
    translated=[]
    allowed={'zero_function','actual_C0_scalar_function_restrict_ell0','actual_chain_scalar_function_restrict_ell0','actual_Res4_SyzQ_scalar_entry_restrict_ell0'}
    for n in old['operator_nodes']:
        assert n[0] in allowed
        translated.append(put([n[0].replace('_restrict_ell0','_restrict_ell_squared')]+n[1:]))
    X=[]
    for j,col in enumerate(manifest['source_X_sparse_columns']):
        assert all(0<=i<56 and 0<=a<len(translated) for i,a in col)
        X.append(put(['actual_X_column_restrict_ell_squared',j,[[i,translated[a]] for i,a in col]]))
    assert len(X)==768
    V=[]
    for j,col in enumerate(adj['actual_V56_source_columns']):
        V.append(put(['actual_V_column_from_same_allorder_chain_restrict_ell_squared',j,col]))
    P=put(['periodic_Koszul_C2_resolution',{'actions':c2path,'basis_size':14,'variables8':['zA'+str(i) for i in range(7)]+['rA'],
      'periodic_d_p':'ell I+(-1)^p L','horizontal_d':'(-1)^p sum_i(x_iI-A_i) exterior_contraction_i',
      'first_ranks':[14,126,518,1302]}])
    a11=put(['literal_a11_source_columns',rp,res['source_a11_columns']])
    d2sel=put(['selected_columns_of_d2',P,sel89])
    Efree=put(['source_standard_columns_F1',126,free26])
    R=put(['nonphysical_target_projection_of_d1_Efree',P,Efree,res['physical_target3_basis_positions']])
    Ecorrect=put(['Efree_minus_a11_R',Efree,a11,R])
    B1=put(['F1_basis_a11_d2sel89_Ecorrect26',a11,d2sel,Ecorrect])
    B1inv=put(['formal_unit_matrix_inverse',B1,{'constant_unit_certificate':'fixed100independent+Pi26(Efree)=I26',
      'certificate_path':'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json'}])
    G1=put(['physical_d1_on_Ecorrect26',P,Ecorrect,rp])
    Z=put(['middle89_coordinates_B1inv_d2_on_F2_complement',P,B1inv,free429])
    I2=put(['F2_complement_inclusion_minus_selected89_Z',518,free429,sel89,Z])
    P2=put(['F2_complement_projection',518,free429])
    P1=put(['last26_coordinates_B1inv',B1inv])
    I1=Ecorrect
    I0=put(['eta_physical3_in_C2_basis14',res['physical_target3_basis_positions']])
    P0=put(['physical_m0_polynomial_projection',rp,res['physical_lift_m0_columns']])
    G2=put(['compose',P1,['d2',P],I2])
    G3=put(['compose',P2,['d3',P]])
    HP=put(['exact_periodic_Koszul_k_linear_homotopy',P,{
      'translation':'T f(x)=f(x+A), inverse T^-1 f=f(x-A); finite matrix derivative degree<=3',
      'Koszul_H':'T^-1 H_ordered_K(x) T; prefixes0,...,7',
      'periodic_h_p':'h_p(u+ell w)=w as a constant vector, for everyp>=0',
      'total_H':'(-1)^p H_Kdelta + i_K h_periodic p_K',
      'p_K':'matrixevaluation f(A,ell), zero on positive exteriordegree',
      'linearity':'continuous k-linear, not claimed B-linear'}])
    HG1=put(['compose',P1,['H_degree0_to1',HP],I0])
    HG2=put(['compose',P2,['H_degree1_to2',HP],I1])
    Mcols=[put(['apply_vector',HG1,v]) for v in V]
    M=put(['matrix_from_columns',[26,56],Mcols])
    MU=put(['selected_source_columns',M,C])
    MUinv=put(['formal_unit_matrix_inverse',MU,{'constant_matrix_and_minor':conpath,'constant_determinant':con['fixed26_minor_determinant'],
      'peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json'}])
    N=put(['selected_source_inclusion_times_inverse',56,C,MUinv])
    Hfree=put(['E_Cfree_minus_N_M_Cfree',56,Cfree,N,M])
    S=put(['kerV_source',[G2,Hfree],N,{'formula':'[N G2,Hfree]','shape':[56,459]}])
    Y=[];bottom=[]
    for j,x in enumerate(X):
        Y.append(put(['apply_vector',HG2,put(['apply_vector',M,x])]))
        bottom.append(put(['select_vector_components',x,Cfree]))
    matrix=put(['complete_T2_presentation_block',{'shape':[459,2070],
      'first1302_columns':['top',G3,'bottom_zero30'],'last768_columns':['top_vectors',Y,'bottom_vectors',bottom]}])
    proof={'base':'B=(k[ell]/ell²)[[zA0,...,zA6,rA]], actualCA-F9 pullback before quotient',
      'restriction_semantics':'each translated scalar source refers to the SAME ORIGINAL regular R0function beforeell, then sA=rA+ell/5 andell²=0; never reuse its ell0value or recompute0/0',
      'native_source_definition':manifest['actual_SyzQ_entry_definition'],
      'original_full_R0_callback_artifacts':original_sources,
      'old_ell0_m_resolution_or_X_coefficient_arrays_reused':False,
      'all_original768_cycles_retained':{'intrinsic28':list(range(28)),'CB612':list(range(28,640)),'complete_native128':list(range(640,768))},
      'source_VX_zero':'restrict exact original UP=QB and complete actual Res4Syzygy identity; targetandrowframes unchanged',
      'unit100':'a11 andactuald2 selected89 are independent atorigin withPi26freeidentity; completeB1is a formalunit matrix',
      'unit89_source':'afterB1change,d2 selected89 areexactmiddleidentity, so F2clear is unipotent',
      'reduced_C2_resolution_first_ranks':[3,26,429,1302],
      'p_i_chainmaps':'explicit target/source basis changes and100contractibleunitpairs; p*i=identity on reducedresolution',
      'homotopy_exactness':'periodic h byepsilon coefficient andorderedKoszulcontraction give dH+Hd=I onpositivecycles; reducedH=pHi',
      'comparison_M':'columnj=HG1(V_j); G1*M=V by actual C2augmentationzero',
      'comparison_M0':'the SAME orderedconstantHadamardprojection as actual26x56peeredmatrix',
      'unit_section_N':'Eselected*(Mselected)^-1, justifiedby peered nonzero constant26minor',
      'MN_identity':True,'VN_equals_G1':True,'comparison_h':'zero becauseMN=I26',
      'free30_split':'B56=N B26 directsum Hfree B30; inverse coordinates (M v,v_Cfree)',
      'kernel_generation':'S(w,c)=N G2w+Hfree c; all kerV because Mv inkerG1=imG2',
      'kernel_relations':'S(w,c)=0 iff c=0 andw inimageG3',
      'actual_X_lift':'Y_j=HG2(MX_j), soG2Y_j=MX_j; S(Y_j,(X_j)_Cfree)=X_j',
      'complete_T2_cokernel':'block [G3,Y768;0,X_Cfree] with459generators2070relations',
      'coefficient_linearity':'H applied only to actualcolumns then extended B-linearly; no false B-linearity ofH is used',
      'finite_base_safety':'onlyAFTERcompletecokernel may tensor by(IAIB)^2; IannK,C implies I²annT2',
      'finite_base_dimension_upper':696,'finite_base_actual_dimension_not_computed':True,
      'annihilator_order_upper':'m5 R/(I,ell²)=0 ->m10 R/(I²,ell²)=0 ->m10T2=0',
      'source_query_order_profile':{
        'ordinary_H_input_order_upper_for_output_degree_q':'q+4, oneordereddivision plusatmost3matrixAi powers; thisis conservative, notnecessity',
        'M_output_degree_q_ordinary_V_input_upper':'q+4',
        'Y_output_degree_q_ordinary_inputs_upper':{'V':'q+8','X':'q+4'},
        'coefficient_degree9_blanket_query_upper':{'ordinary_V':17,'ordinary_X':13},
        'not_authorized_or_claimed_needed':'nohighjetquery follows; useactualsourceadjoint/low-unitcancellation tolowerandbound first',
        'honest_quotients':'ordinaryregularX input canitselfneed highernumeratororders; preserveactualsourcefunctions',
        'existing_low_cache_not_complete':'cachedB56conormalconstantdoesnot supplyfullhigherM,N,Ycoefficients'},
      'new_scalar_or_T2_length_or_H0_query_performed':False}
    evidence={'operator_nodes':nodes,'source_X768_nodes':X,'actual_V56_nodes':V,
      'actual_M26_by56_node':M,'actual_N56_by26_node':N,'actual_Hfree56_by30_node':Hfree,
      'actual_G1_G2_G3_nodes':[G1,G2,G3],'actual_S56_by459_node':S,
      'all768_actual_cycle_lift_top_nodes':Y,'all768_actual_cycle_lift_bottom_nodes':bottom,
      'complete_actual_T2_presentation_node':matrix,'source_operator_proof':proof}
    ep=out/'actual-honest-T2-459-by2070-factored-source-operator.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    assert len(nodes)+10271+len(Y)*2+459+2070<150000 and ep.stat().st_size<32*1024*1024
    report={'id':STEM,'all_passed':True,'scope':'complete factored formal-source operator; scalar coefficients and T2length not queried',
      'complete_actual_T2_presentation_shape':[459,2070],'reduced_C2_first_ranks':[3,26,429,1302],
      'all_original768_cycles_retained':True,'comparison_m_n_h_source_defined':True,'new_physical_coefficients_evaluated':False,
      'actual_T2_or_H0_length_computed':False,'operator_nodes':len(nodes),'bytes':ep.stat().st_size,
      'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','complete_actual_T2_presentation_shape','operator_nodes','bytes','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
