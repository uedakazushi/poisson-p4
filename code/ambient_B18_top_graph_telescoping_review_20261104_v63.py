#!/usr/bin/env python3
"""Read-only exact circuit audit of actual graph lifts; no source rerun."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview

ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_TOP_GRAPH_TELESCOPING_20261104_v63'
STEM=BASE+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def solve(M,B):
    n=len(M);a=[list(map(F,row))+list(map(F,b)) for row,b in zip(M,B)];det=F(1)
    for j in range(n):
        i=next(i for i in range(j,n) if a[i][j])
        if i!=j:a[i],a[j]=a[j],a[i];det=-det
        p=a[j][j];det*=p;a[j]=[x/p for x in a[j]]
        for i in range(n):
            if i!=j and a[i][j]:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return [r[n:] for r in a],det


def closure(nodes,roots):
    wanted=set();stack=list(roots)
    while stack:
        i=stack.pop()
        if i in wanted:continue
        wanted.add(i);n=nodes[i]
        if n[0] in ['n','i']:stack.append(n[1])
        elif n[0] in ['a','m']:stack+=list(n[1:])
    return wanted


def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json')
    assert primary['all_passed'] and primary['actual_CB_Hodd_10_generators_Jacobian_rank']==9
    assert primary['nonunit_regular_divisions_only_with_exact_coordinate_restriction_zero_proof']==54
    assert not primary['depth_nil_T3_membership_or_socle_answer_claimed']
    cert=read(folder/'actual-top12-DAG.json');nodes=[tuple(n) for n in cert['DAG_scalar_nodes']]
    oldfolder=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63'
    oldpath=oldfolder/'actual-source-coefficients.json';old=read(oldpath)
    peerpath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW/report.json'
    peer=read(peerpath);assert peer['all_passed']
    assert peer['independent_sparse_physical_source_GL_M24_and_full48_Jacobian_blocks']
    assert peer['all150_w0to3_finite_source_DAG_binding_replay_byte_identical']
    assert cert['old_DAG_input_sha256']==sha(oldpath)
    nold=len(old['DAG_scalar_nodes']);assert nodes[:nold]==[tuple(n) for n in old['DAG_scalar_nodes']]
    D=DAGReview(old['DAG_scalar_nodes'])
    cb=old['families']['CB'];ca=old['families']['CA']
    xids=[D.variables[f'x{row}'] for row in old['physical_zero_base_rows9']]
    zids=[i for p in cb['coefficient_unknown_pairs'].values() for i in p]
    zids+=[cb['p_quadratic_relation'][k] for k in ['s','t']]
    assert len(zids)==48 and len(set(xids+zids))==57
    eq=cb['all48_zero_IFT_quadratic_coefficient_equations']
    J=[[D.grad[i].get(v,F(0)) for v in zids] for i in eq]
    rhs=[[-D.grad[i].get(v,F(0)) for v in xids] for i in eq]
    zx,detJ=solve(J,rhs)
    gate=cert['coordinate_gate']
    assert J==[list(map(F,r)) for r in gate['source48_actual_AD_Jacobian']]
    assert zx==[list(map(F,r)) for r in gate['source48_dz_dx']]
    assert detJ==F(gate['source48_Jacobian_determinant'])==F(-331776,625)
    def implicit_gradient(node):
        answer=[D.grad[node].get(x,F(0)) for x in xids]
        for i,z in enumerate(zids):
            c=D.grad[node].get(z,F(0))
            if c:answer=[a+c*b for a,b in zip(answer,zx[i])]
        return answer
    generators=cb['actual_zero_odd_I0_generators']+cb['actual_positive_odd_J0_generators']
    matrix=[implicit_gradient(node) for row,node in generators]
    assert matrix==[list(map(F,r)) for r in gate['actual_CB_10_by9_odd_Jacobian']]
    assert gate['generator_physical_rows']==[18,19,28,29,60,67,72,63,46,95]
    selected=gate['selected_generator_indices'];assert selected==[0,1,2,3,4,5,6,8,9]
    assert not any(matrix[7])
    inv,detF=solve([matrix[i] for i in selected],[[F(i==j) for j in range(9)] for i in range(9)])
    assert detF==F(1,96) and detJ*detF==F(-3456,625)
    gids=[node for row,node in generators];Y=[gids[i] for i in selected]
    assert Y==cert['selected_Y_coordinate_nodes'] and all(D.values[i]==0 for i in gids)
    bids=[pair[1] for row,pair in cb['actual_individual_physical_cubic_trace_pairs']]
    centers=[]
    for node in bids:
        assert D.values[node]==0
        g=implicit_gradient(node)
        centers.append([sum(g[k]*inv[k][j] for k in range(9)) for j in range(9)])
    # Extend only scalar origin values; no large AD on the nine cloned systems.
    for i,n in enumerate(nodes[nold:],nold):
        op=n[0]
        if op=='c':value=F(n[1],n[2])
        elif op=='v':value=F(0)
        elif op=='r':
            assert len(n)==5 and 0<=n[3]<6 and 0<=n[4]<9
            assert n[1]<i and n[2]<i
            value=centers[n[3]][n[4]]
        else:
            assert all(isinstance(j,int) and 0<=j<i for j in n[1:])
            a=n[1]
            if op=='n':value=-D.values[a]
            elif op=='a':value=D.values[a]+D.values[n[2]]
            elif op=='m':value=D.values[a]*D.values[n[2]]
            elif op=='i':assert D.values[a];value=1/D.values[a]
            else:raise ValueError(n)
        D.values.append(value)
    D.nodes=nodes;D.cache={n:i for i,n in enumerate(nodes)}
    assert len(D.cache)==len(nodes)==176910
    assert sum(n[0]=='r' for n in nodes)==54
    needed=closure(nodes[:nold],eq+[gids[i] for i in selected]+bids)
    assert all(i in xids+zids for i in needed if nodes[i][0]=='v')
    partial=[]
    for level,instance in enumerate(cert['all9_pointed_unit57_inverse_instances'],1):
        assert instance['level']==level
        bindings=dict(instance['variable_bindings']);assert set(bindings)==set(xids+zids)
        assert all(nodes[new]==('v',f'INV{level}_{nodes[oldid][1]}') for oldid,new in bindings.items())
        mapping={}
        for i in sorted(needed):
            n=nodes[i];op=n[0]
            if op=='c':new=D.const(F(n[1],n[2]))
            elif op=='v':new=bindings[i]
            elif op=='n':new=D.neg(mapping[n[1]])
            elif op=='a':new=D.add(mapping[n[1]],mapping[n[2]])
            elif op=='m':new=D.mul(mapping[n[1]],mapping[n[2]])
            elif op=='i':new=D.inverse(mapping[n[1]])
            else:raise ValueError(n)
            mapping[i]=new
        expected=[mapping[i] for i in eq]
        expected += [D.add(mapping[gids[index]],D.neg(Y[k] if k<level else D.zero))
                     for k,index in enumerate(selected)]
        assert expected==instance['all57_equations'] and all(D.values[i]==0 for i in expected)
        assert F(instance['unit57_Jacobian_determinant_in_z_then_x_order'])==detJ*detF
        partial.append([mapping[i] for i in bids])
    assert len(partial)==9
    ca0=ca['actual_zero_odd_I0_generators'];cap=ca['actual_positive_odd_J0_generators']
    assert [r for r,node in ca0]==[r for r,node in cb['actual_zero_odd_I0_generators']]
    assert [r for r,node in cap]==[46,95]
    assert cert['all_lower_CB_kernel_and_graph_relations_must_still_be_kept']
    for alpha,hrec in enumerate(cert['all6_actual_CB_odd_membership_lifts']):
        row=cb['actual_individual_physical_cubic_trace_pairs'][alpha][0]
        assert hrec['raw_cubic_row']==row
        hs=[D.zero]*10;previous=D.zero
        for i,index in enumerate(selected):
            numerator=D.add(partial[i][alpha],D.neg(previous))
            hid=D.lookup(('r',numerator,Y[i],alpha,i))
            h=hrec['regular_coordinate_divisions9'][i]
            assert h['h_node']==hid and h['numerator']==numerator and h['coordinate']==Y[i]
            assert F(h['origin_value'])==centers[alpha][i]
            hs[index]=hid;previous=partial[i][alpha]
        assert hs==hrec['membership_coefficients10'] and hs[7]==D.zero
        aA,bA=ca['actual_individual_physical_cubic_trace_pairs'][alpha][1]
        aB,bB=cb['actual_individual_physical_cubic_trace_pairs'][alpha][1]
        even=D.add(aB,D.neg(aA));odd=D.neg(bA)
        for j,(physical,fA) in enumerate(ca0):
            even=D.add(even,D.mul(hs[j],D.add(ca['complete_normalized_physical_w_coefficients'][0][physical][0],
                                             D.neg(cb['complete_normalized_physical_w_coefficients'][0][physical][0]))))
            odd=D.add(odd,D.mul(hs[j],fA))
        for j,(physical,fA) in enumerate(cap):
            even=D.add(even,D.mul(hs[8+j],D.add(ca['complete_normalized_physical_w_coefficients'][1][physical][0],
                                               D.neg(cb['complete_normalized_physical_w_coefficients'][1][physical][0]))))
            odd=D.add(odd,D.mul(hs[8+j],fA))
        rec=cert['all6_actual_relative_cubic_pairs_in_same_CA_trace_basis'][alpha]
        assert rec=={'raw_cubic_row':row,'same_CA_tau_basis_CB_graph_minus_CA_pair':[even,odd]}
        for j,(physical,fA) in enumerate(cap):
            product=[D.mul(ca['actual_tracefree_nu'],D.mul(fA,odd)),D.mul(fA,even)]
            assert cert['all12_actual_top_T3_pairs'][2*alpha+j]=={
                'positive_odd_raw_row':physical,'cubic_raw_row':row,'pair_nodes':product}
    checked={}
    for p in [folder/'report.json',peerpath]:
        data=read(p)
        for name,want in data.get('input_sha256',{}).items():
            target=ROOT/name;assert sha(target)==want,str(target);checked[name]=want
        for name,want in data.get('evidence_sha256',{}).items():
            target=p.parent/name;assert sha(target)==want,str(target);checked[str(target.relative_to(ROOT))]=want
    proof={'implicit_CB_10by9_Jacobian':[[str(x) for x in r] for r in matrix],
        'det48':str(detJ),'selected9_indices':selected,'det9':str(detF),'det57':str(detJ*detF),
        'six_by9_regular_quotient_origin_values':[[str(x) for x in row] for row in centers],
        'all9_inverse_clones_and_all513_equations_exactly_bound':True,
        'all54_regular_coordinate_division_nodes_bound':True,
        'restriction_zero_proof':'consecutive pointed unit57 systems have identical targets after Yi=0; uniqueness, regular coordinate and domain give unique algebraic quotient',
        'all6_relative_graph_pairs_and12_products_reconstructed_in_CA_trace_basis':True,
        'T3_membership_socle_depth_nil_answer':False,'author_source_or_followup_Q186_rerun':False}
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),
           folder/'report.json',folder/'actual-top12-DAG.json',folder/'rank-gate.json',oldpath,peerpath]
    for p in paths:checked[str(p.relative_to(ROOT))]=sha(p)
    result={'id':STEM,'all_passed':True,
        'independent_original_DAG_AD_and_physical_implicit_10by9_gate':True,
        'selected9_determinant':'1/96','pointed_inverse57_determinant':'-3456/625',
        'all9_inverse_instances_and_all57_equations_each_exact_DAG_binding':True,
        'regularity_of54_coordinate_quotients_without_unit_Y_or_p_inverse':True,
        'all6_actual_memberships_and_relative_graphs_in_same_CA_tau_basis':True,
        'all12_top_T3_pair_finite_algebraic_circuits_verified':True,
        'full_lower_CB_kernel_and_graph_relations_retained':True,
        'T3_membership_depth_socle_or_full_nil_claim':False,
        'author_full_source_newpoint_GB_or_Q186_followup_job':False,
        'input_sha256':checked,'evidence_sha256':{'independent-proof.json':sha(ep)},
        'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
