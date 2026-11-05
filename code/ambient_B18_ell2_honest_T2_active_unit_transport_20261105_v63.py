#!/usr/bin/env python3
"""Known-source active unit transport/tensor reduction of honest T2.

No new rank search, physical coefficient, or scalar homotopy query.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def guard(v):
    assert all(max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=384 for x in v.values());return v
def fixed_quotient(n,cols):
    piv={};receipts=[]
    for j,orig in enumerate(cols):
        v=dict(orig);steps=[]
        while v and min(v) in piv:
            p=min(v);c=v[p];steps.append([p,str(c)]);v=guard(add(v,piv[p],-c))
        assert v,'known fixed source columns must be independent'
        p=min(v);c=v[p];piv[p]=guard(scale(v,1/c));receipts.append({'source_column':j,'pivot':p,'divide_by':str(c),'subtract':steps,'normalized':srec(piv[p])})
    free=[i for i in range(n) if i not in piv];fi={i:j for j,i in enumerate(free)};red={}
    for i in reversed(range(n)):
        if i in fi:red[i]={fi[i]:F(1)}
        else:
            q={}
            for j,c in piv[i].items():
                if j!=i:q=add(q,red[j],-c)
            red[i]=guard(q)
    Pi=[red[i] for i in range(n)];assert all(not mv(Pi,v) for v in cols)
    return free,Pi,receipts
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hs={}
    def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    oldpath='reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63/actual-honest-T2-459-by2070-factored-source-operator.json'
    old=read(oldpath)
    oldpeer=read('reports/AMBIENT_B18_ELL2_HONEST_T2_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json');assert oldpeer['all_passed']
    pref=read('reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json')
    rank=read('reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_RANK_20261105_v63/active-d3-rank-and-all-source-factorizations.json')
    rpeer=read('reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_RANK_20261105_v63_REVIEW/report.json');assert rpeer['all_passed']
    first=read('reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json')
    respath='reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json';res=read(respath)
    adj=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json')
    con=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/actual-V-conormal26-by56-unit-minor-and-source-factorization.json')
    assert rank['rank']==126 and rank['actual_beta2_C2_over_B']==102
    D1=list(map(parse,pref['active_d1_columns']));D2=list(map(parse,pref['active_d2_columns']));D3=list(map(parse,pref['active_d3_columns']))
    split1=pref['signed_zero_wedge_decomposition']['1'];split2=pref['signed_zero_wedge_decomposition']['2']
    selected56=[split2[j][1] for j in first['fixed_original_minor_columns'] if not split2[j][0]]
    assert len(selected56)==len(set(selected56))==56
    selected126=rank['fixed_original_minor_columns'];assert len(selected126)==126
    a0=[];a_poly=[]
    for v in res['source_a11_columns']:
        p={};q=[]
        for row,pol in v:
            z,a,sign=split1[row];assert not z
            q.append([a,[[m,str(F(c)*sign)] for m,c in pol]])
            c=sum((F(c) for m,c in pol if m==[0,0]),F(0))*sign
            if c:p[a]=c
        a0.append(p);a_poly.append(q)
    fixed1=a0+[D2[j] for j in selected56]
    free17,Pi17,rec1=fixed_quotient(84,fixed1);assert len(free17)==17 and all(not mv(Pi17,v) for v in D2)
    fixed2=[{j:F(1)} for j in selected56]+[D3[j] for j in selected126]
    free42,Pi42,rec2=fixed_quotient(224,fixed2);assert len(free42)==42
    assert all(not mv(Pi42,v) for v in D3)
    eta=res['physical_target3_basis_positions'];nonphys=[i for i in range(14) if i not in eta]
    I1active=[]
    for j in free17:
        v={j:F(1)}
        for k,i in enumerate(nonphys):v=add(v,a0[k],-D1[j].get(i,F(0)))
        assert not mv(D1,v) and mv(Pi17,v)=={len(I1active):F(1)}
        I1active.append(v)
    inverse1={(tuple(z),i):(j,sign) for j,(z,i,sign) in enumerate(split1)}
    I1full=[]
    for v in I1active:
        I1full.append({inverse1[((),i)][0]:c*inverse1[((),i)][1] for i,c in v.items()})
    for z in pref['zero_actions']:
        for i in eta:
            row,sign=inverse1[((z,),i)];I1full.append({row:F(sign)})
    assert len(I1full)==26
    Piold=list(map(parse,adj['Pi26_source126_columns']));Tcols=[mv(Piold,v) for v in I1full]
    # The frame isomorphism is known invertible; save an exact two-sided inverse.
    aug=[[Tcols[j].get(i,F(0)) for j in range(26)]+[F(i==j) for j in range(26)] for i in range(26)]
    det=F(1);invreceipts=[]
    for j in range(26):
        k=next(i for i in range(j,26) if aug[i][j])
        if k!=j:aug[k],aug[j]=aug[j],aug[k];det=-det
        p=aug[j][j];det*=p;aug[j]=[c/p for c in aug[j]];ops=[]
        for i in range(26):
            if i!=j and aug[i][j]:
                c=aug[i][j];ops.append([i,str(c)]);aug[i]=[a-c*b for a,b in zip(aug[i],aug[j])]
        for row in aug:guard({i:c for i,c in enumerate(row) if c})
        invreceipts.append({'column':j,'swap_with':k,'normalize_by':str(p),'clear':ops})
    inversecols=[{i:aug[i][26+j] for i in range(26) if aug[i][26+j]} for j in range(26)]
    assert all(mv(Tcols,v)=={j:F(1)} for j,v in enumerate(inversecols))
    assert all(mv(inversecols,v)=={j:F(1)} for j,v in enumerate(Tcols))
    Mold=list(map(parse,con['actual26_by56_columns']));Mnew=[guard(mv(inversecols,v)) for v in Mold]
    selectedJ=con['fixed_selected26_columns'];factors=list(map(parse,con['all56_columns_in_selected26_basis']))
    assert all(mv([Mnew[j] for j in selectedJ],v)==m for v,m in zip(factors,Mnew))
    newdet=F(con['fixed26_minor_determinant'])/det;guard({0:newdet});assert newdet
    constant={'active11_source_columns_polynomial_r_ell':a_poly,'active56_source_columns':selected56,'active126_source_columns':selected126,
      'active_F1_fixed67_columns':[srec(v) for v in fixed1],'active_F1_fixed67_receipts':rec1,
      'active_F1_free17':free17,'active_Pi17_columns':[srec(v) for v in Pi17],
      'active_F2_fixed182_columns':[srec(v) for v in fixed2],'active_F2_fixed182_receipts':rec2,
      'active_F2_free42':free42,'active_Pi42_columns':[srec(v) for v in Pi42],
      'new_full_I1_origin26_columns':[srec(v) for v in I1full],
      'frame_T_old_to_new_definition':'T=p_old1*i_new1, soG1new=G1old*T',
      'T0_columns':[srec(v) for v in Tcols],'T0_inverse_columns':[srec(v) for v in inversecols],
      'T0_determinant':str(det),'T0_inverse_receipts':invreceipts,
      'actual_Mnew0_columns':[srec(v) for v in Mnew],'same_selected26_V_columns':selectedJ,
      'new_selected26_minor_determinant':str(newdet),'new_to_old_determinant_ratio':'newminor=oldminor/detT0',
      'all56_new_column_factorizations_in_same_selectedJ_basis':[srec(v) for v in factors]}
    cp=out/'actual-active-unit-bases-and26-frame-transport.json';cp.write_text(json.dumps(constant,separators=(',',':'))+'\n')
    nodes=[]
    def put(n):nodes.append(n);return len(nodes)-1
    AP=put(['active5_periodic_Koszul_C2_resolution',pref['active_variables'],{'initial_ranks':[14,84,224,364,434]}])
    B1=put(['active_F1_unit_basis',[['a11',a_poly],['d2_columns',selected56],['Efree17_minus_a11_nonphysical_d1',free17]],{'constant_certificate':cp.name,'block_sizes':[11,56,17]}])
    B2=put(['active_F2_unit_basis',[['Eselected56',selected56],['d3_columns',selected126],['Efree42_minus_Eselected56_middle56_B1inv_d2',free42]],{'constant_certificate':cp.name,'block_sizes':[56,126,42]}])
    comp238=[j for j in range(364) if j not in selected126]
    B3=put(['active_F3_source_unit_clear',selected126,comp238,'middle126_B2inv_d3',{'block_sizes':[126,238]}])
    AG=put(['active_reduced_complete_resolution',AP,B1,B2,B3,{'first_ranks':[3,17,42,238],
      'd4_and_higher':'retainoriginaluppercomplexwithsameunit-chaintransports; no newd4rank'}])
    tensor=put(['tensor_active_resolution_with_Koszul_zero3',AG,pref['zero_actions'],{
      'basis_order':'activefactor first, zero-wedge factor last','sign':'activehomologicaldegree in zeroKoszul differential',
      'first_ranks':[3,26,102,418],'higher_completeness':'tensorofcompleteunit-equivalentactive resolutionwithregularzero3Koszul'}])
    i_new1=put(['tensor_active_unit_chain_inclusion_degree1',tensor,{'origin':cp.name}])
    p_new2=put(['tensor_active_unit_chain_projection_degree2',tensor,{'target_rank':102}])
    oldnodes=old['operator_nodes'];oldp1=next(i for i,n in enumerate(oldnodes) if n[0]=='last26_coordinates_B1inv')
    oldHP=next(i for i,n in enumerate(oldnodes) if n[0]=='exact_periodic_Koszul_k_linear_homotopy')
    T=put(['compose',['external_operator',oldpath,oldp1],i_new1])
    Tinv=put(['formal_unit_inverse',T,{'constant_certificate':cp.name,'constant_determinant':str(det)}])
    M=put(['compose',Tinv,['external_operator',oldpath,old['actual_M26_by56_node']]])
    N=put(['compose',['external_operator',oldpath,old['actual_N56_by26_node']],T])
    Hfree=put(['external_operator',oldpath,old['actual_Hfree56_by30_node']])
    H2=put(['compose',p_new2,['H_degree1_to2_external',oldpath,oldHP],i_new1])
    Y=[]
    for x in old['source_X768_nodes']:Y.append(put(['apply_vector',H2,put(['apply_vector',M,['external_operator',oldpath,x]])]))
    pres=put(['complete_honest_T2_presentation',{'shape':[132,1186],
      'first418_columns':['G3_of_tensor_resolution',tensor,'bottomzero30'],
      'last768_columns':['top102vectornodes',Y,'bottom30sourcevectors',oldpath,old['all768_actual_cycle_lift_bottom_nodes']]}])
    proof={'complete_source_base':'B=(k[ell]/ell²)[[zA0,...,zA6,rA]]',
      'active_unit_pair_counts':[11,56,126],
      'all_active_unit_bases_at_origin_certified_without_new_ranksearch':True,
      'formal_unit_basis_definitions':'B1/B2 chosenfromactuald2/d3sourcecolumns; inverseunitscertifiedbyfixedconstantbases; B3clearunipotent',
      'active_resolution_first_ranks':[3,17,42,238],
      'tensor_zero3_removes_additional_degree3_4_unit_pairs':11+3*56+3*126,
      'global_resolution_first_ranks':[3,26,102,418],
      'rank418_not_claimed_minimal_beta3':True,
      'frame_change':'T=p_old1*i_new1 is a formalunit; G1new=G1old*T; exactconstantT0/twosidedinversechecked',
      'comparison_change':'Mnew=T^-1 Mold, Nnew=Nold*T, Hfreeunchanged; MnewNnew=I26 andVNnew=G1new',
      'cycle_lifts':'Ynew_j=p_new2 H_original i_new1(Mnew X_j); all768sameoriginalfullR0cycles',
      'exact_complete_T2_cokernel_shape':[132,1186],
      'finite_coefficient_or_T2_length_queries':False,
      'unchanged_source_order_caution':'unittransportdoesnotmakehighercoefficientsknown; ordinaryinputordersmustbeadjoint-profiledbeforequeries'}
    ep=out/'actual-honest-T2-132-by1186-factored-source-operator.json';ep.write_text(json.dumps({'operator_nodes':nodes,'actual_Mnew_node':M,'actual_Nnew_node':N,'actual_Hfree_node':Hfree,'all768_cycle_lift_top_nodes':Y,'complete_T2_presentation_node':pres,'source_proof':proof},separators=(',',':'))+'\n')
    slots=sum(len(v) for v in fixed1+fixed2+Pi17+Pi42+I1full+Tcols+inversecols+Mnew)+sum(len(q['subtract']) for q in rec1+rec2)
    assert slots+len(nodes)+10000<150000 and sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    report={'id':STEM,'all_passed':True,'complete_actual_T2_presentation_shape':[132,1186],
      'full_C2_initial_resolution_ranks':[3,26,102,418],'active_unit_pairs':[11,56,126],
      'constant26_frame_change_determinant':str(det),'new_V_selected26_minor':str(newdet),
      'all_original768_cycles_retained':True,'actual_T2_or_H0_length_computed':False,'new_physical_coefficients_queried':False,
      'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in [cp,ep]},'whole_named_slots':slots+len(nodes)+10000,'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','complete_actual_T2_presentation_shape','constant26_frame_change_determinant','new_V_selected26_minor','whole_named_slots','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
