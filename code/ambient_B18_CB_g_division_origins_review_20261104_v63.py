#!/usr/bin/env python3
"""Independent sparse-dictionary two-variable microjet of the SAME unit57 system."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure,solve

ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63';STEM=BASE+'_REVIEW'
MON=[(d-i,i) for d in range(4) for i in range(d+1)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def plus(a,b):
    c=dict(a)
    for e,q in b.items():c[e]=c.get(e,F(0))+q
    return {e:q for e,q in c.items() if q}
def scale(a,c):return {e:c*q for e,q in a.items() if c*q}
def times(a,b):
    c={}
    for e,p in a.items():
        for f,q in b.items():
            g=(e[0]+f[0],e[1]+f[1])
            if sum(g)<=3:c[g]=c.get(g,F(0))+p*q
    return {e:q for e,q in c.items() if q}
def invert(a):
    c=a.get((0,0),F(0));assert c
    h=scale(a,1/c);h=plus(h,{(0,0):F(-1)})
    answer={(0,0):F(1)};power={(0,0):F(1)}
    for degree in range(1,4):
        power=times(power,h);answer=plus(answer,scale(power,(-1)**degree))
    return scale(answer,1/c)
def vector(a):return [a.get(e,F(0)) for e in MON]


def evaluate(nodes,needed,bindings):
    values={}
    for i in needed:
        n=nodes[i];op=n[0]
        if op=='c':p={(0,0):F(n[1],n[2])} if n[1] else {}
        elif op=='v':p=bindings[i]
        elif op=='n':p=scale(values[n[1]],-1)
        elif op=='a':p=plus(values[n[1]],values[n[2]])
        elif op=='m':p=times(values[n[1]],values[n[2]])
        elif op=='i':p=invert(values[n[1]])
        else:raise ValueError(n)
        values[i]=p
    return values


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json');cert=read(folder/'degree3-certificate.json')
    assert primary['all_passed'] and primary['single_pointed57_system_only']
    assert primary['homogeneous_orders_computed']==[1,2,3]
    assert primary['origin_values_by_physical_row']=={'64':'0','68':'0','73':'0','74':'0'}
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    data=read(sourcepath);nodes=data['DAG_scalar_nodes'];cb=data['families']['CB']
    assert cert['single_actual_original_source_DAG_hash']==sha(sourcepath)
    fullpath=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json';full=read(fullpath)
    fullpeer=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json'
    assert read(fullpeer)['all_passed']
    divisions=full['all4_known_regular_g_divisions']
    assert all(q['origin_value']=='not evaluated; no unit use' for q in divisions)
    D=DAGReview(nodes)
    xids=[D.variables[f'x{r}'] for r in data['physical_zero_base_rows9']]
    zids=[i for p in cb['coefficient_unknown_pairs'].values() for i in p]+[cb['p_quadratic_relation'][k] for k in ['s','t']]
    variables=zids+xids
    assert variables==cert['single57_unknown_source_node_ids'] and len(set(variables))==57
    assert [nodes[i][1] for i in variables]==cert['single57_unknown_source_variable_names']
    eq48=cb['all48_zero_IFT_quadratic_coefficient_equations']
    gens=cb['actual_zero_odd_I0_generators']+cb['actual_positive_odd_J0_generators']
    selected=[0,1,2,3,4,5,6,8,9];F9=[gens[i][1] for i in selected]
    g=gens[7][1];brows=[64,68,73,74]
    bnodes=[cb['complete_normalized_physical_w_coefficients'][0][r][1] for r in brows]
    assert F9==cert['selected_F9_original_nodes'] and eq48==cert['original48_equation_nodes']
    assert g==cert['g_CB_f63_original_node']
    equations=eq48+F9
    needed=sorted(closure(nodes,equations+[g]+bnodes))
    assert needed==cert['original_source_needed_node_ids'] and len(needed)==9978
    assert {i for i in needed if nodes[i][0]=='v'}==set(variables)
    J=[[D.grad[i].get(v,F(0)) for v in variables] for i in equations]
    assert J==[list(map(F,r)) for r in cert['fixed_J57']]
    _,det=solve(J,[[] for _ in range(57)]);assert det==F(-3456,625)
    inverse=[list(map(F,r)) for r in cert['fixed_J57_inverse']]
    for A,B in [(J,inverse),(inverse,J)]:
        assert all(sum(A[i][k]*B[k][j] for k in range(57) if A[i][k] and B[k][j])==int(i==j)
                   for i in range(57) for j in range(57))
    # Independent polynomial representation and inversion by finite geometric series.
    targets=[{} for _ in range(55)]+[{(1,0):F(1)},{(0,1):F(1)}]
    solution={i:{} for i in variables};correction_receipts=[]
    for degree in range(1,4):
        v=evaluate(nodes,needed,solution)
        exponents=[e for e in MON if sum(e)==degree]
        residual=[plus(v[row],scale(target,-1)) for row,target in zip(equations,targets)]
        correction=[[-sum(inverse[i][k]*residual[k].get(e,F(0)) for k in range(57) if inverse[i][k])
                     for e in exponents] for i in range(57)]
        saved=cert['homogeneous_Newton_corrections'][degree-1]
        assert saved['order']==degree and saved['monomials']==list(map(list,exponents))
        assert correction==[list(map(F,r)) for r in saved['correction57_by_order_plus1']]
        for i,var in enumerate(variables):solution[var]=plus(solution[var],dict(zip(exponents,correction[i])))
        correction_receipts.append({'degree':degree,'independent_sparse_correction_matches':True})
    values=evaluate(nodes,needed,solution)
    assert all(values[row]==target for row,target in zip(equations,targets))
    saved_jets={i:list(map(F,v)) for i,v in cert['all57_unknown_jets']}
    assert all(vector(solution[i])==saved_jets[i] for i in variables)
    rr={(1,0):F(1)};ss={(0,1):F(1)}
    cubic=scale(times(plus(rr,scale(ss,-1)),times(plus(scale(rr,2),scale(ss,3)),plus(scale(rr,2),scale(ss,3)))),-6)
    assert values[g]==cubic
    assert vector(cubic)==list(map(F,cert['g_CB_all10_coefficients']))
    assert vector(cubic)==[F(0)]*6+[F(-24),F(-48),F(18),F(54)]
    derived=[]
    for row,node,division,saved in zip(brows,bnodes,divisions,cert['all4_regular_quotient_origin_values']):
        assert not values[node] and saved['all10_bbar_coefficients']==['0']*10
        assert division['node']==saved['known_regular_full_DAG_g_quotient_node']
        assert full['DAG_scalar_nodes'][division['node']][3]==f'graph{row}:known_I0_membership'
        assert saved['physical_graph_row']==row and saved['old_source_odd_coefficient_node']==node
        h0=values[node].get((3,0),F(0))/cubic[(3,0)]
        assert h0==0==F(saved['quotient_origin_value'])
        derived.append({'physical_row':row,'full_DAG_g_node':division['node'],'derived_exact_origin':str(h0)})
    # Confirm the frozen primary full-DAG origin metadata was not overwritten.
    assert all(q['origin_value']=='not evaluated; no unit use' for q in read(fullpath)['all4_known_regular_g_divisions'])
    checked={}
    for p in [folder/'report.json',fullpeer]:
        d=read(p)
        for name,want in d.get('input_sha256',{}).items():
            t=ROOT/name;assert sha(t)==want,str(t);checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():
            t=p.parent/name;assert sha(t)==want,str(t);checked[str(t.relative_to(ROOT))]=want
    manifest=read(folder/'manifest.json')
    for name,want in manifest['files_sha256'].items():
        p=ROOT/name;assert sha(p)==want,str(p);checked[name]=want
    proof={'independent_sparse_dictionary_jets_with_geometric_series_unit_inverse':True,
        'same57_system_Jacobian_determinant':str(det),'J_and_Jinv_both_products_checked':True,
        'all9978_original_ancestor_nodes_bound':True,'all57_equations_all10_monomials_zero':True,
        'homogeneous_correction_checks':correction_receipts,
        'g_CB_cubic_coefficients':list(map(str,vector(cubic))),
        'g_CB_cubic_factorization':'-6*(r-s)*(2*r+3*s)^2',
        'exact_origin_overrides_separate_from_frozen_primary':derived,
        'known_global_regular_divisibility_used_before_jet_constant_extraction':True,
        'full_quotient_zero_or_depth_claim':False,'higher_jet_or_author_rerun':False,
        'next_actual_colon_work':'coordinate with manuscript for full186 residue/pivot and symmetric annihilator preflight; not duplicated here'}
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',
              folder/'degree3-certificate.json',folder/'manifest.json',sourcepath,fullpath,fullpeer]:
        checked[str(p.relative_to(ROOT))]=sha(p)
    result={'id':STEM,'all_passed':True,
        'independent_same57_Jacobian_and_sparse_bivariate_inverse_through3':True,
        'all57_equations_and_saved_unknown_jets_independently_verified':True,
        'g_CB_cubic':'-6*(r-s)*(2*r+3*s)^2',
        'all4_bbar_jets_through3_zero_and_exact_regular_quotient_origins_zero':True,
        'global_membership_not_inferred_from_jets':True,
        'frozen_full_DAG_None_metadata_preserved':True,
        'full_quotient_zero_higher_order_colon_depth_nil_claim':False,
        'newpoint_source_rerun_GB_or_higherjet':False,
        'input_sha256':checked,'evidence_sha256':{'independent-proof.json':sha(ep)},
        'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
