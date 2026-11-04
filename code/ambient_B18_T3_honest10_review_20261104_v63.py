#!/usr/bin/env python3
"""Independent pure-unit telescopes, exact lifts and complete H10 kernel audit."""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip,hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure,solve
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63_REVIEW'
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    if out.exists():
        assert not any(out.iterdir()), 'refuse to overwrite evidence'
    else:out.mkdir(parents=True)
    base=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63'
    primary=read(base/'report.json');cert=read(base/'honest10-presentation-DAG.json');assert primary['all_passed']
    for rel,h in primary['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in primary['evidence_sha256'].items():assert sha(base/rel)==h,rel
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fullpath=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    originpath=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    graphpath=ROOT/'reports/AMBIENT_B18_NOETHER_FIBER_EVEN_QUADRATIC_20261103_v63/certificate.json.gz'
    weightpath=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json'
    peers=[ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW/report.json',
           ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json',
           ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json']
    assert all(read(p)['all_passed']for p in peers)
    source,full,orig,graph,weight=map(read,(sourcepath,fullpath,originpath,graphpath,weightpath))
    oldnodes=source['DAG_scalar_nodes'];fullnodes=full['DAG_scalar_nodes'];nodes=[tuple(n)for n in cert['DAG_scalar_nodes']]
    assert nodes[:len(fullnodes)]==[tuple(n)for n in fullnodes]
    assert fullnodes[:len(oldnodes)]==oldnodes
    D=DAGReview(oldnodes);ca=source['families']['CA']
    xids=[D.variables[f'x{r}']for r in source['physical_zero_base_rows9']]
    az=[i for pair in ca['coefficient_unknown_pairs'].values()for i in pair]+[ca['p_quadratic_relation'][k]for k in('s','t')]
    gens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']
    selected=[0,1,2,3,4,5,6,8,9];Y=[gens[i][1]for i in selected];f63=gens[7][1]
    assert Y==cert['actual_R0_coordinate_nodes_Y9']
    rows=graph['independent87_physical_graph_rows'];assert len(rows)==87
    roots=list(dict.fromkeys([f63]+[ca['complete_normalized_physical_w_coefficients'][weight['physical150_torus_weights'][r]][r][1]for r in rows]))
    assert len(roots)==88 and all(D.values[i]==0 for i in roots+Y)
    eq48=ca['all48_zero_IFT_quadratic_coefficient_equations']
    J=[[D.grad[i].get(z,F(0))for z in az]for i in eq48]
    rhs=[[-D.grad[i].get(x,F(0))for x in xids]for i in eq48]
    zx,det48=solve(J,rhs)
    def dx(node):
        return [D.grad[node].get(x,F(0))+sum((D.grad[node].get(z,F(0))*zx[k][j]for k,z in enumerate(az)),F(0))for j,x in enumerate(xids)]
    F9=[dx(i)for i in Y];invF,det9=solve(F9,[[F(i==j)for j in range(9)]for i in range(9)])
    assert det9==F(128,3)and det48*det9==F(-3623878656,625)
    derivatives={i:[sum((dx(i)[k]*invF[k][j]for k in range(9)),F(0))for j in range(9)]for i in roots}
    needed=sorted(closure(oldnodes,eq48+Y+roots))
    assert {i for i in needed if oldnodes[i][0]=='v'}<=set(az+xids)
    assert all(oldnodes[i][0]in('c','v','n','a','m','i')for i in needed)
    # All added coordinate divisions have independently calculated origin
    # derivatives. No old r/g function is cloned or divided again.
    telescopes=cert['all_pure_source_Hadamard_membership_telescopes']
    assert [r['original_source_function']for r in telescopes]==roots
    new_divisions={r['node']:(root,j,r)for root,tel in zip(roots,telescopes)for j,r in enumerate(tel['all9_regular_coordinate_divisions'])}
    assert len(new_divisions)==792
    D.values=list(map(F,orig['all339669_scalar_origin_values_evaluated']))
    assert len(D.values)==len(fullnodes)
    for i,n in enumerate(nodes[len(fullnodes):],len(fullnodes)):
        op=n[0]
        if op=='c':v=F(n[1],n[2])
        elif op=='v':v=F(0)
        elif op=='r':
            assert i in new_divisions
            root,j,rec=new_divisions[i]
            assert n==('r',rec['numerator'],Y[j],f'H10_CA_source{root}:coordinate{j}')
            v=derivatives[root][j];assert v==F(rec['origin'])
        elif op=='n':v=-D.values[n[1]]
        elif op=='a':v=D.values[n[1]]+D.values[n[2]]
        elif op=='m':v=D.values[n[1]]*D.values[n[2]]
        elif op=='i':assert D.values[n[1]];v=1/D.values[n[1]]
        else:raise ValueError(n)
        D.values.append(v)
    assert {i for i,n in enumerate(nodes)if i>=len(fullnodes)and n[0]=='r'}==set(new_divisions)
    D.nodes=nodes;D.cache={n:i for i,n in enumerate(nodes)};assert len(nodes)==486437
    # The producer seeds .memo whereas RegularDAG uses .cache. Some old
    # scalar/variable keys therefore reappear as aliases; no global
    # hash-cons uniqueness is assumed. Explicit binding checks follow.
    duplicate_keys=len(nodes)-len(D.cache);assert duplicate_keys==96
    prefixes=[]
    for level,system in enumerate(cert['all8_pure_CA_prefix57_unit_systems'],1):
        assert system['level']==level
        binding=dict(system['bindings']);assert set(binding)==set(az+xids)
        assert all(nodes[new]==('v',f'H10_CA_PREFIX{level}_{oldnodes[old][1]}')for old,new in binding.items())
        mapped={}
        for i in needed:
            n=oldnodes[i];op=n[0]
            if op=='c':new=D.const(F(n[1],n[2]))
            elif op=='v':new=binding[i]
            elif op=='n':new=D.neg(mapped[n[1]])
            elif op=='a':new=D.add(mapped[n[1]],mapped[n[2]])
            elif op=='m':new=D.mul(mapped[n[1]],mapped[n[2]])
            elif op=='i':new=D.inverse(mapped[n[1]])
            else:raise ValueError(n)
            mapped[i]=new
        expected=[mapped[i]for i in eq48]+[D.add(mapped[i],D.neg(Y[j]if j<level else D.zero))for j,i in enumerate(Y)]
        assert expected==system['all57_unit_equations']
        assert all(D.values[i]==0 for i in expected)
        assert F(system['determinant_z48_then_x9'])==det48*det9
        prefixes.append(mapped)
    seed={D.zero:(D.zero,)*9}
    for j,i in enumerate(Y):seed[i]=tuple(D.one if k==j else D.zero for k in range(9))
    for root,tel in zip(roots,telescopes):
        previous=D.zero;vector=[]
        for j,rec in enumerate(tel['all9_regular_coordinate_divisions']):
            current=prefixes[j][root]if j<8 else root
            assert rec['current_source_restriction']==current and rec['previous_source_restriction']==previous
            assert rec['numerator']==D.add(current,D.neg(previous))and rec['coordinate']==Y[j]
            vector.append(rec['node']);previous=current
        assert previous==root and vector==tel['membership_vector9'];seed[root]=tuple(vector)
    # Independent proof-rule interpreter, not the author's ideal-factor search.
    rules={int(i):q for i,q in cert['odd_membership_propagation_rules'].items()};memo={}
    def vector(i):
        if i in memo:return memo[i]
        rule=rules[i];kind=rule[0]
        if kind in('zero','coordinate','pure_CA_unit57_telescope'):ans=seed[i]
        elif kind=='negative_of_certified_function':
            assert nodes[rule[1]]==('n',i);ans=tuple(D.neg(q)for q in vector(rule[1]))
        elif kind=='neg':
            assert nodes[i]==('n',rule[1]);ans=tuple(D.neg(q)for q in vector(rule[1]))
        elif kind=='sum':
            assert nodes[i]==('a',*sorted(rule[1:]));ans=tuple(D.add(a,b)for a,b in zip(vector(rule[1]),vector(rule[2])))
        elif kind=='product_left_known_ideal_factor':
            assert nodes[i]==('m',*sorted(rule[1:]));ans=tuple(D.mul(rule[2],q)for q in vector(rule[1]))
        elif kind=='product_right_known_ideal_factor':
            assert nodes[i]==('m',*sorted(rule[1:]));ans=tuple(D.mul(rule[1],q)for q in vector(rule[2]))
        else:raise ValueError(rule)
        memo[i]=ans;return ans
    for i in rules:vector(i)
    source_columns=[full['actual_Qn_R0_squared_generator_matrices']['3'],full['actual12_top_T3_additional_columns'],full['actual_Jinf_R0_squared_generator_matrix']]
    lift_keys=['all2034_Q3_columns_lifted_R0_10','all12_top_columns_lifted_R0_10','all2046_Jinf_columns_lifted_R0_10']
    for original,key in zip(source_columns,lift_keys):
        assert [[a,*vector(b)]for a,b in original]==cert[key]
    assert Counter(map(tuple,source_columns[2]))==Counter(map(tuple,source_columns[0]+source_columns[1]))
    koszul=[]
    for i in range(9):
        for j in range(i+1,9):
            col=[D.zero]*10;col[i+1]=D.neg(Y[j]);col[j+1]=Y[i];koszul.append(col)
    assert koszul==cert['all36_Koszul_kernel_columns_R0_10']
    assert cert['H10_generator_images_in_actual_CA_trace_pairs']==[[D.one,D.zero]]+[[D.zero,i]for i in Y]
    assert koszul+cert[lift_keys[0]]==cert['honest_D3_presentation_relation_columns_R0_10']
    assert koszul+cert[lift_keys[0]]+cert[lift_keys[1]]==cert['honest_D4_presentation_relation_columns_R0_10']
    proof=dict(source_binding='original full physical CA48 Hensel system and the actual nine regular odd-coordinate functions',
        det_CA48=str(det48),det_Y9=str(det9),det_CA57=str(det48*det9),
        pure_original_source_functions=len(roots),prefix57_systems=8,regular_coordinate_divisions=792,
        all_prefix_equations_and_telescoping_numerators_rebuilt=True,
        all_original_r_g_functions_left_as_scalars_no_pullback=True,
        independently_checked_propagation_rules=len(rules),exact_lifted_column_counts=[2034,12,2046],
        whole_kernel='H=R0 directsum m9*tau; map R0^10 sends E0 to1 and Ei toYi*tau. Normalizer basis1,tau is free2. Y9 is a regular sequence, so all first syzygies are the36 Koszul columns.',
        exact_D3='coker(R0^2070->R0^10)',exact_D4='coker(R0^2082->R0^10)',
        actual_T3='image of the12 added relation columns inside D3, equivalently ker(D3->D4)',
        not_complete_T3_presentation='Relations among those12 images still require syzygies against all2070 D3 columns. Reduction of these matrices modulo m^2 alone loses Tor.',
        no_T3_nonzero_depth_or_nil_claim=True)
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'honest10-presentation-DAG.json',
            ROOT/'code/ambient_B18_T3_honest10_presentation_20261104_v63.py',
            ROOT/'research/notes/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63.md',
            sourcepath,fullpath,originpath,graphpath,weightpath,*peers,
            ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',
            ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    report=dict(id=ID,all_passed=True,independent_source48_AD_and_actual_Y9_regular_coordinate_binding=True,
        all8_pure_CA57_systems_and792_regular_quotients_independently_bound=True,
        all3993_exact_propagation_rules_and4092_lifts_independently_verified=True,
        complete_H10_kernel_exactly36_Koszul_by_regular_sequence=True,
        exact_D3_free10_relation_columns=2070,exact_D4_free10_relation_columns=2082,
        actual_T3_image12_only_relations_still_uncomputed=True,
        no_old_nonunit_function_pullback_new_jets_points_GB_membership_or_depth_claim=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={ep.name:sha(ep)},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','seconds','max_rss_kib')},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
