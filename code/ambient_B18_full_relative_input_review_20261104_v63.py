#!/usr/bin/env python3
"""Read-only all99+87/module circuit audit, with explicit nonunit-division scopes."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure,solve

ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63';STEM=BASE+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json')
    assert primary['all_passed'] and primary['actual_Q_generator_count']==186
    assert not primary['g_division_origin_constants_evaluated']
    assert not primary['full_nil_depth_socle_T3_membership_or_primary_answer_claimed']
    cert=read(folder/'full-relative186-DAG.json');nodes=[tuple(n) for n in cert['DAG_scalar_nodes']]
    oldpath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    old=read(oldpath);ca,cb=old['families']['CA'],old['families']['CB']
    nold=len(old['DAG_scalar_nodes']);assert nodes[:nold]==[tuple(n) for n in old['DAG_scalar_nodes']]
    peerpath=ROOT/'reports/AMBIENT_B18_TOP_GRAPH_TELESCOPING_20261104_v63_REVIEW/report.json'
    assert read(peerpath)['all_passed']
    mupath=ROOT/'reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json'
    mu=read(mupath);assert mu['all_passed'] and mu['CB_actual_I_minimal10_and_image_depth21_edim39_e2_accepted']
    relativepeer=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW/report.json'
    assert read(relativepeer)['all_passed']
    gatepath=ROOT/'reports/AMBIENT_B18_TOP_GRAPH_TELESCOPING_20261104_v63/rank-gate.json';gate=read(gatepath)
    D=DAGReview(old['DAG_scalar_nodes'])
    xids=[D.variables[f'x{i}'] for i in old['physical_zero_base_rows9']]
    zids=[i for p in cb['coefficient_unknown_pairs'].values() for i in p]+[cb['p_quadratic_relation'][k] for k in ['s','t']]
    zx=[list(map(F,row)) for row in gate['source48_dz_dx']]
    def deriv(node):
        answer=[D.grad[node].get(x,F(0)) for x in xids]
        for k,z in enumerate(zids):
            c=D.grad[node].get(z,F(0))
            if c:answer=[a+c*b for a,b in zip(answer,zx[k])]
        return answer
    gens=cb['actual_zero_odd_I0_generators']+cb['actual_positive_odd_J0_generators']
    gids=[n for r,n in gens];selected=gate['selected_generator_indices']
    assert selected==[0,1,2,3,4,5,6,8,9]
    inv,detF=solve([deriv(gids[i]) for i in selected],[[F(i==j) for j in range(9)] for i in range(9)])
    assert detF==F(1,96)
    def dy(node):
        g=deriv(node);return [sum(g[k]*inv[k][j] for k in range(9)) for j in range(9)]
    Y=[gids[i] for i in selected];r,s=Y[7:];f63=gids[7]
    assert not any(dy(f63))
    graphpath=ROOT/'reports/AMBIENT_B18_NOETHER_FIBER_EVEN_QUADRATIC_20261103_v63/certificate.json.gz'
    graph=read(graphpath)
    weights=read(ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json')['physical150_torus_weights']
    rows=graph['independent87_physical_graph_rows']
    assert len(rows)==len(set(rows))==87
    assert Counter(weights[i] for i in rows)=={0:4,1:53,2:24,3:6}
    assert [i for i in rows if weights[i]==0]==[64,68,73,74]
    Bnodes={row:cb['complete_normalized_physical_w_coefficients'][weights[row]][row][1] for row in rows}
    assert all(D.values[n]==0 for n in Bnodes.values())
    rrecords={v['node']:v for v in cert['all_coordinate_divisions']}
    grecords={v['node']:v for v in cert['all4_known_regular_g_divisions']}
    assert len(rrecords)==784 and len(grecords)==4
    has_g=[False]*nold
    for i,n in enumerate(nodes[nold:],nold):
        op=n[0]
        if op=='c':value=F(n[1],n[2]);gflag=False
        elif op=='v':value=F(0);gflag=False
        elif op in ['r','g']:
            assert len(n)==4 and n[1]<i and n[2]<i
            if op=='r':
                rec=rrecords[i];assert n[1]==rec['numerator'] and n[2]==rec['coordinate']
                value=F(rec['origin_value']);gflag=has_g[n[1]] or has_g[n[2]]
            else:
                rec=grecords[i];assert n[1]==rec['numerator'] and n[2]==rec['denominator_g']
                assert rec['origin_value']=='not evaluated; no unit use'
                value=None;gflag=True
        else:
            assert all(isinstance(j,int) and 0<=j<i for j in n[1:])
            gflag=any(has_g[j] for j in n[1:]);a=n[1]
            if op=='i':
                assert not gflag and D.values[a] is not None and D.values[a]!=0
                value=1/D.values[a]
            elif op=='n':value=None if D.values[a] is None else -D.values[a]
            else:
                b=n[2]
                value=None if D.values[a] is None or D.values[b] is None else (
                    D.values[a]+D.values[b] if op=='a' else D.values[a]*D.values[b])
                assert op in ['a','m']
        D.values.append(value);has_g.append(gflag)
    D.nodes=nodes;D.cache={n:i for i,n in enumerate(nodes)}
    assert len(D.cache)==len(nodes)==339669
    assert {i for i,n in enumerate(nodes) if n[0]=='g'}==set(grecords)
    assert {i for i,n in enumerate(nodes) if n[0]=='r'}==set(rrecords)
    needed=closure(nodes[:nold],cb['all48_zero_IFT_quadratic_coefficient_equations']+[gids[i] for i in selected]+list(Bnodes.values())+[f63])
    full=tuple(range(9));tail=(7,8)
    subsets={tuple(range(i)) for i in range(1,9)}
    subsets.update(tuple(list(range(i))+[7,8]) for i in range(7));subsets.add((8,))
    instances={full:None};records=cert['all16_pointed_unit57_inverse_instances']
    assert [tuple(v['retained_Y_indices']) for v in records]==sorted(subsets)
    for number,rec in enumerate(records,1):
        sub=tuple(rec['retained_Y_indices']);bindings=dict(rec['bindings'])
        assert set(bindings)==set(xids+zids)
        assert all(nodes[new]==('v',f'FULLINV{number}_{nodes[oldid][1]}') for oldid,new in bindings.items())
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
        expected=[mapping[i] for i in cb['all48_zero_IFT_quadratic_coefficient_equations']]
        expected += [D.add(mapping[gids[index]],D.neg(Y[j] if j in sub else D.zero)) for j,index in enumerate(selected)]
        assert expected==rec['all57_pointed_equations']
        assert F(rec['unit_Jacobian'])==F(-3456,625)
        instances[sub]=mapping
    def at(node,sub):
        sub=tuple(sub)
        if not sub:assert D.values[node]==0;return D.zero
        if sub==full:return node
        return instances[sub][node]
    seen_r=set()
    def cdiv(num,coord,label,value):
        node=D.lookup(('r',num,coord,label));rec=rrecords[node]
        assert F(rec['origin_value'])==value
        seen_r.add(node);return node
    def zero7(node,label):
        previous=at(node,tail);out=[];d=dy(node)
        for i in range(7):
            current=at(node,tuple(list(range(i+1))+[7,8]))
            out.append(cdiv(D.add(current,D.neg(previous)),Y[i],f'{label}:z{i}',d[i]));previous=current
        return out
    g=at(f63,tail);e=zero7(f63,'f63');gs=at(f63,(8,))
    A=cdiv(D.add(g,D.neg(gs)),r,'g_split_r',F(0));B=cdiv(gs,s,'g_split_s',F(0))
    assert cert['CB_g_restriction_node']==g and cert['CB_f63_zero7_telescoping_coefficients']==e
    assert cert['CB_g_split_rA_plus_sB_nodes']==[A,B]
    assert cert['CB_coordinate_Y9_nodes']==Y and cert['CB_f63_node']==f63
    graph_pairs=[];seen_g=set()
    stored_members=dict(cert['all87_actual_odd_membership_lifts'])
    for row in rows:
        node=Bnodes[row];degree=weights[row];hs=[D.zero]*10
        if degree:
            previous=D.zero;d=dy(node)
            for i,index in enumerate(selected):
                current=at(node,tuple(range(i+1)))
                hs[index]=cdiv(D.add(current,D.neg(previous)),Y[i],f'graph{row}:Y{i}',d[i]);previous=current
        else:
            h8=D.lookup(('g',at(node,tail),g,f'graph{row}:known_I0_membership'));seen_g.add(h8)
            hi=zero7(node,f'zero_graph{row}')
            for i in range(7):hs[i]=D.add(hi[i],D.neg(D.mul(h8,e[i])))
            hs[7]=h8
        assert hs==stored_members[row]
        aB,bB=cb['complete_normalized_physical_w_coefficients'][degree][row]
        aA,bA=ca['complete_normalized_physical_w_coefficients'][degree][row]
        even=D.add(aB,D.neg(aA));odd=D.neg(bA)
        for i,(raw,fA) in enumerate(ca['actual_zero_odd_I0_generators']):
            even=D.add(even,D.mul(hs[i],D.add(ca['complete_normalized_physical_w_coefficients'][0][raw][0],D.neg(cb['complete_normalized_physical_w_coefficients'][0][raw][0]))))
            odd=D.add(odd,D.mul(hs[i],fA))
        for j,(raw,jA) in enumerate(ca['actual_positive_odd_J0_generators']):
            even=D.add(even,D.mul(hs[8+j],D.add(ca['complete_normalized_physical_w_coefficients'][1][raw][0],D.neg(cb['complete_normalized_physical_w_coefficients'][1][raw][0]))))
            odd=D.add(odd,D.mul(hs[8+j],jA))
        graph_pairs.append({'kind':'graph_CB_minus_CA','raw_row':row,'weight':degree,'normalized_CA_trace_pair':[even,odd]})
    assert seen_r==set(rrecords) and seen_g==set(grecords)
    def pa(a,b):return (D.add(a[0],b[0]),D.add(a[1],b[1]))
    def pn(a):return tuple(D.neg(i) for i in a)
    def ps(a,c):return tuple(D.mul(i,c) for i in a)
    nuA=ca['actual_tracefree_nu'];nuB=cb['actual_tracefree_nu']
    def pm(a,b):return (D.add(D.mul(a[0],b[0]),D.mul(nuA,D.mul(a[1],b[1]))),D.add(D.mul(a[0],b[1]),D.mul(a[1],b[0])))
    rawodd=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators'];Z=[]
    for i,(row,f) in enumerate(rawodd):
        degree=int(i>=8);aA,bA=ca['complete_normalized_physical_w_coefficients'][degree][row]
        aB,bB=cb['complete_normalized_physical_w_coefficients'][degree][row]
        Z.append((D.add(aA,D.neg(aB)),bA))
    originalZ=Z[:]
    for i in range(7):Z[7]=pa(Z[7],pn(ps(originalZ[i],e[i])))
    bn=Y[:7]+[g,r,s]
    assert [list(z) for z in Z]==cert['CB_new_odd_generators10_evaluated_on_CA_normalized_pairs']
    assert bn==cert['CB_scalar_odd_generators_normalized_coefficients']
    kernel=[]
    def save(label,degree,p):kernel.append({'kind':'CB_kernel','label':label,'weight':degree,'normalized_CA_trace_pair':list(p)})
    for i,j in it.combinations(range(7),2):save(f'regular_Koszul_{i}_{j}',0,pa(ps(Z[j],bn[i]),pn(ps(Z[i],bn[j]))))
    for i in range(7):
        for j in [7,8,9]:save(f'cross_{i}_{j}',int(j>=8),pa(ps(Z[j],bn[i]),pn(ps(Z[i],bn[j]))))
    save('HB_wg_minus_A_wr_minus_B_ws',1,pa(pa(Z[7],pn(ps(Z[8],A))),pn(ps(Z[9],B))))
    save('HB_s_wr_minus_r_ws',1,pa(ps(Z[8],s),pn(ps(Z[9],r))))
    assert len(kernel)==44
    for i in range(10):
        for j in range(i,10):
            p=pm(Z[i],Z[j]);p=(D.add(p[0],D.neg(D.mul(nuB,D.mul(bn[i],bn[j])))),p[1])
            save(f'product_{i}_{j}',int(i>=8)+int(j>=8),p)
    assert Counter(v['weight'] for v in kernel)=={0:64,1:32,2:3}
    Q=kernel+graph_pairs
    assert Q==cert['all99_CB_kernel_and87_graph_generators_evaluated_on_CA']
    assert Counter(v['weight'] for v in Q)=={0:68,1:85,2:27,3:6}
    H=[(D.one,D.zero)]+[(D.zero,n) for row,n in rawodd];C0=H[:9]
    assert [list(h) for h in H]==cert['actual_CA_H_generators']
    assert [list(h) for h in C0]==cert['actual_CA_C0_generators']
    def span(indices,multipliers):return [list(pm(tuple(Q[i]['normalized_CA_trace_pair']),h)) for i in indices for h in multipliers]
    idx={d:[i for i,q in enumerate(Q) if q['weight']==d] for d in range(4)}
    Jinf=span(range(186),H)
    Qn={str(n):span([i for d in range(n) for i in idx[d]],H)+span(idx[n],C0) for n in range(4)}
    assert Jinf==cert['actual_Jinf_R0_squared_generator_matrix']
    assert Qn==cert['actual_Qn_R0_squared_generator_matrices']
    top=[list(pm(tuple(q['normalized_CA_trace_pair']),h)) for q in graph_pairs if q['weight']==3 for h in H[-2:]]
    assert top==cert['actual12_top_T3_additional_columns'] and len(top)==12
    counts=[len(Qn[str(i)]) for i in range(4)]
    assert counts==[9*68,11*68+9*85,11*153+9*27,11*180+9*6]==[612,1513,1926,2034]
    assert len(Jinf)==11*186==2046
    checked={}
    for p in [folder/'report.json',peerpath,mupath,relativepeer]:
        d=read(p)
        for name,want in d.get('input_sha256',{}).items():
            t=ROOT/name;assert sha(t)==want,str(t);checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():
            t=p.parent/name;assert sha(t)==want,str(t);checked[str(t.relative_to(ROOT))]=want
    proof={'all16_unit57_systems_all912_equations_bound':True,'coordinate_divisions_bound':784,
        'regular_g_divisions_bound':4,'g_origins_remain_unknown_and_never_inverted':True,
        'g_nonzero_gate':'actual mu(I_CB)=10 excludes g=0',
        'regular_g_membership_gate':'complete actual minimal image algebra generation; nonnegative weights; closed I0; restrict z=0',
        'full_syzygy_count':44,'full_product_count':55,'full_graph_count':87,
        'Q_weight_counts':dict(Counter(v['weight'] for v in Q)),
        'literal_Qn_column_counts':counts,'literal_Jinf_column_count':len(Jinf),
        'all8131_module_columns_and12_extra_columns_DAG_bound':True,
        'membership_colon_socle_depth_nil_answer':False,'author_source_or_full_builder_rerun':False}
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',folder/'full-relative186-DAG.json',peerpath,mupath,relativepeer]:
        checked[str(p.relative_to(ROOT))]=sha(p)
    result={'id':STEM,'all_passed':True,
        'complete_actual44_syzygies_plus55_products_kernel_proof_accepted':True,
        'g_nonzero_mu10_and_four_regular_divisions_without_circular_graph_assumption':True,
        'all16_inverse_instances_and784_coordinate_quotients_bound':True,
        'unknown_g_division_origins_never_assigned_zero_or_used_as_units':True,
        'all99_CB_kernel_and87_actual_graph_pairs_same_CA_trace_basis_verified':True,
        'nuB_product_RHS_and_nuA_pair_multiplication_distinguished':True,
        'all_four_Qn_and_Jinf_literal_module_columns_verified':True,
        'Qn_column_counts':counts,'Jinf_column_count':2046,
        'source_rerun_newpoint_GB_actual_module_rank_depth_nil_claim':False,
        'input_sha256':checked,'evidence_sha256':{'independent-proof.json':sha(ep)},
        'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
