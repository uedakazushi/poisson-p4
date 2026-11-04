#!/usr/bin/env python3
"""Independent opaque-atom weighted algebra for the saved LOW3 negative test."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_CB_g_division_origins_review_20261104_v63 import evaluate as source_evaluate,MON
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_T3_SPECIALIZED_LOW3_OBSTRUCTION_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/AMBIENT_B18_T3_SPECIALIZED_LOW3_OBSTRUCTION_20261104_v63'
    primary=read(base/'report.json');cert=read(base/'specialized-low3-certificate.json')
    assert primary['all_passed']and not primary['T3_nonzero_and_depth19_proved']
    for rel,h in primary['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in primary['evidence_sha256'].items():assert sha(base/rel)==h,rel
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fullpath=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    ownpath=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/independent-certificate.json'
    ownreport=ownpath.parent/'report.json';assert read(ownreport)['all_passed']
    gpeer=ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63_REVIEW/report.json'
    assert read(gpeer)['all_passed']
    old=read(sourcepath);full=read(fullpath);nodes=full['DAG_scalar_nodes'];oldnodes=old['DAG_scalar_nodes']
    assert nodes[:len(oldnodes)]==oldnodes
    solution={i:{m:F(q)for m,q in zip(MON,v)if F(q)}for i,v in read(ownpath)['all105_unknown_sparse_jets']}
    originpath=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    origins=read(originpath)['all339669_scalar_origin_values_evaluated']
    atomic={q['node']for q in full['all_coordinate_divisions']}|{q['node']for q in full['all4_known_regular_g_divisions']}
    order={i:1 for i in atomic}
    for i in full['CB_g_split_rA_plus_sB_nodes']:assert i in atomic;order[i]=2
    for rec in full['all_coordinate_divisions']:assert F(origins[rec['node']])==F(rec['origin_value'])
    for rec in full['all4_known_regular_g_divisions']:assert F(origins[rec['node']])==0
    qcols=full['actual_Qn_R0_squared_generator_matrices']['3'];top=full['actual12_top_T3_additional_columns']
    assert len(qcols)==2034 and len(top)==12
    gtail=full['CB_g_restriction_node'];needed=set();terminal=set();todo=[i for q in qcols+top for i in q]
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i)
        if i<len(oldnodes):terminal.add(i);continue
        if i in atomic or i==gtail:continue
        op,*args=nodes[i]
        if op in('a','m'):todo+=args
        elif op in('n','i'):todo.append(args[0])
        elif op!='c':raise ValueError(nodes[i])
    cb_rs=[i for row,i in old['families']['CB']['actual_positive_odd_J0_generators']]
    assert cb_rs==cert['CB_actual_positive_r_s_original_nodes']
    usedold=sorted(closure(oldnodes,list(terminal)+cb_rs))
    values=source_evaluate(oldnodes,usedold,solution)
    def actual(i):return {(m[0],m[1],()):q for m,q in values[i].items()if q}
    def degree(m):return m[0]+m[1]+sum(order[i]for i in m[2])
    def add(a,b):
        q=dict(a)
        for m,c in b.items():q[m]=q.get(m,F(0))+c
        return {m:c for m,c in q.items()if c}
    def scale(a,c):return {m:c*v for m,v in a.items()if c*v}
    def mul(a,b):
        q={}
        for m,c in a.items():
            for n,d in b.items():
                e=(m[0]+n[0],m[1]+n[1],tuple(sorted(m[2]+n[2])))
                if degree(e)<=3:q[e]=q.get(e,F(0))+c*d
        return {m:c for m,c in q.items()if c}
    one={(0,0,()):F(1)}
    # Recompose the known CB cubic using actual CB-coordinate LINEAR terms
    # on the CA plane; no CA/CB coordinate names are identified.
    rB={m:c for m,c in actual(cb_rs[0]).items()if degree(m)==1}
    sB={m:c for m,c in actual(cb_rs[1]).items()if degree(m)==1}
    tail=scale(mul(add(rB,scale(sB,-1)),mul(add(scale(rB,2),scale(sB,3)),add(scale(rB,2),scale(sB,3)))),-6)
    tail_expected={(mon.count(0),mon.count(1),()):F(q)for mon,q in cert['CB_g_tail_exact_cubic_polynomial']}
    assert tail==tail_expected
    pol={}
    for i in sorted(needed):
        if i<len(oldnodes):p=actual(i)
        elif i in atomic:
            p={(0,0,(i,)):F(1)};c=F(origins[i])
            if c:p[(0,0,())]=c
        elif i==gtail:p=tail
        else:
            op,*args=nodes[i]
            if op=='c':p={(0,0,()):F(args[0],args[1])}if args[0]else{}
            elif op=='n':p=scale(pol[args[0]],-1)
            elif op=='a':p=add(pol[args[0]],pol[args[1]])
            elif op=='m':p=mul(pol[args[0]],pol[args[1]])
            elif op=='i':
                c=pol[args[0]].get((0,0,()),F(0));assert c
                h=scale(add(pol[args[0]],{(0,0,()):-c}),1/c)
                p=dict(one);power=dict(one)
                for k in range(1,4):power=mul(power,h);p=add(p,scale(power,(-1)**k))
                p=scale(p,1/c)
            else:raise ValueError(nodes[i])
        pol[i]=p
    position={m:i for i,m in enumerate(MON)}
    def vectors(pair):
        both={}
        for row,node in enumerate(pair):
            for m,c in pol[node].items():both.setdefault(m,[F(0),F(0)])[row]+=c
        known=[F(0)]*20;unc=[]
        for(r,s,atoms),co in both.items():
            if not atoms:
                for row in range(2):known[10*row+position[r,s]]+=co[row]
            else:
                lo=sum(order[i]for i in atoms)
                for mr,ms in MON:
                    if lo<=mr+ms<=3-r-s:
                        v=[F(0)]*20
                        for row in range(2):v[10*row+position[r+mr,s+ms]]=co[row]
                        if any(v):unc.append(v)
        return known,unc
    def shift(v,m):
        w=[F(0)]*20
        for row in range(2):
            for n,j in position.items():
                dest=(m[0]+n[0],m[1]+n[1])
                if dest in position:w[10*row+position[dest]]+=v[10*row+j]
        return w
    basis={}
    def remainder(v):
        v=v[:]
        for j,row in sorted(basis.items()):
            if v[j]:v=[a-v[j]*b for a,b in zip(v,row)]
        return v
    def insert(v):
        v=remainder(v);j=next((i for i,c in enumerate(v)if c),None)
        if j is not None:
            c=v[j];basis[j]=[a/c for a in v]
    records=cert['all2034_Q_column_known_and_enlarged_uncertain_vectors']
    for index,pair in enumerate(qcols):
        known,unc=vectors(pair);saved=records[index]
        assert saved['column']==index and known==list(map(F,saved['known_part']))
        assert Counter(map(tuple,unc))==Counter(tuple(map(F,q))for q in saved['uncertainty_vectors'])
        for v in [known,*unc]:
            for m in MON:insert(shift(v,m))
    assert len(basis)==11
    for _,v in cert['Q_enlarged_span_basis']:assert not any(remainder(list(map(F,v))))
    assert all(not pol[node]for pair in top for node in pair)
    assert not cert['negative_nonmembership_columns'] and cert['literal_formal_top3_all_zero']
    assert all(not any(map(F,q['known_part']))and not q['uncertainty_vectors']for q in cert['all12_top_results'])
    evidence=dict(independent_source_jet_namespace='own frozen joint105 source jets; no Newton solve',
        independent_formal_encoding='(r exponent,s exponent,opaque atomic-node multiset), minimum atom orders retained',
        original_source_ancestor_count=len(usedold),post_source_node_count=len(needed-terminal),
        all2034_Q_columns_and_uncertainty_vector_multisets_match=True,
        rank_enlarged_Q=11,all12_top_pairs_formally_zero_through3=True,
        CB_tail_composed_from_actual_CB_linear_r_s_and_independent_reviewed_cubic=True,
        unknown_division_derivatives_never_zeroed=True,
        scope='No negative certificate was found. This does not prove membership, T3=0, depth20, or nil structure. Method stops at3.',
        exact_presentation_next='Construct an abstract presentation of Jinf, then reduce its syzygies and the lifted Q3 relations over A=R0/(I0A+I0B). Do not replace Jinf tensor A by its ambient image in A^2.')
    ep=out/'independent-proof.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'specialized-low3-certificate.json',
            ROOT/'code/ambient_B18_T3_specialized_low3_obstruction_20261104_v63.py',
            ROOT/'research/notes/AMBIENT_B18_T3_SPECIALIZED_LOW3_OBSTRUCTION_20261104_v63.md',sourcepath,fullpath,ownpath,ownreport,gpeer,originpath,
            ROOT/'code/ambient_B18_CB_g_division_origins_review_20261104_v63.py']
    report=dict(id=ID,all_passed=True,independent_all2034_enlarged_Q_vectors_and_rank11=True,
        independent_formal_atomic_algebra_all_top12_zero_through3=True,
        actual_CB_tail_binding_and_all_unknown_division_deviations_retained=True,
        no_negative_certificate_and_no_actual_membership_T3_zero_or_depth_conclusion=True,
        no_source_Newton_degree4_newpoint_GB_or_author_rerun=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={ep.name:sha(ep)},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','seconds','max_rss_kib')},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
