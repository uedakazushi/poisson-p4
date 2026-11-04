#!/usr/bin/env python3
"""ONE literal honest H=R0+ m9*tau presentation, retaining the Koszul kernel.

Only pure-unit CA source functions are pulled back.  Old regular r/g coefficient
functions stay on the original base and are never redivided after restriction.
"""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_full_relative_input_20261104_v63 import RegularDAG
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors,forward,linear_solve,replay
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fullp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    originp=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    graphp=ROOT/'reports/AMBIENT_B18_NOETHER_FIBER_EVEN_QUADRATIC_20261103_v63/certificate.json.gz'
    lengths=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63/report.json'
    fullpeer=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json'
    old=read(source);full=read(fullp);orig=read(originp);assert read(lengths)['all_passed'] and read(fullpeer)['all_passed']
    import gzip
    graph=json.loads(gzip.decompress(graphp.read_bytes())); rows=graph['independent87_physical_graph_rows']
    oldnodes=old['DAG_scalar_nodes'];nodes=full['DAG_scalar_nodes'];assert nodes[:len(oldnodes)]==oldnodes
    ca=old['families']['CA'];varids={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'}
    xids=[varids[f'x{i}'] for i in old['physical_zero_base_rows9']]
    az=[i for pair in ca['coefficient_unknown_pairs'].values() for i in pair]+[ca['p_quadratic_relation']['s'],ca['p_quadratic_relation']['t']]
    gens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']
    selected=[0,1,2,3,4,5,6,8,9];Y=[gens[i][1] for i in selected]; f63=gens[7][1]
    # All graph odd functions are original CA source functions, not composed r/g functions.
    weightpath=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json'
    weights=read(weightpath)['physical150_torus_weights']
    graphB={r:ca['complete_normalized_physical_w_coefficients'][weights[r]][r][1] for r in rows}
    roots=list(dict.fromkeys([f63]+list(graphB.values())))
    assert len(rows)==87 and len(roots)<=88
    vals,grad=forward(oldnodes,set(xids+az));assert all(vals[i]==0 for i in roots+Y)
    eq48=ca['all48_zero_IFT_quadratic_coefficient_equations']
    J=[[grad[i].get(v,F(0)) for v in az] for i in eq48]
    B=[[-grad[i].get(v,F(0)) for v in xids] for i in eq48]
    zx,det48=linear_solve(J,B)
    def dx(i):
        a=[grad[i].get(v,F(0)) for v in xids]
        for k,z in enumerate(az):
            c=grad[i].get(z,F(0))
            if c:a=[u+c*v for u,v in zip(a,zx[k])]
        return a
    Fj=[dx(i) for i in Y];invF,det9=linear_solve(Fj,[[F(int(i==j)) for j in range(9)] for i in range(9)])
    assert det9==F(128,3) and det48*det9==F(-3623878656,625)
    derivatives={i:[sum(dx(i)[k]*invF[k][j] for k in range(9)) for j in range(9)] for i in roots}
    D=RegularDAG();assert D.zero==0 and D.one==1
    D.nodes=[tuple(n) for n in nodes];D.values=list(map(F,orig['all339669_scalar_origin_values_evaluated']))
    D.memo={n:i for i,n in enumerate(D.nodes)}
    assert len(D.nodes)==len(D.values) and len(D.memo)==len(D.nodes)
    needed=ancestors(oldnodes,eq48+Y+roots)
    assert {i for i in needed if oldnodes[i][0]=='v'}<=set(xids+az)
    prefix=[];systems=[]
    for level in range(1,9):
        binding={i:D.var(f'H10_CA_PREFIX{level}_{oldnodes[i][1]}') for i in az+xids}
        m=replay(D,oldnodes,binding,needed);prefix.append(m)
        equations=[m[i] for i in eq48]+[D.add(m[i],D.neg(Y[j] if j<level else D.zero)) for j,i in enumerate(Y)]
        assert len(equations)==57 and all(D.values[i]==0 for i in equations)
        systems.append({'level':level,'bindings':[[i,j] for i,j in binding.items()],
                        'all57_unit_equations':equations,'determinant_z48_then_x9':str(det48*det9)})
    zero9=(D.zero,)*9;membership={D.zero:zero9};proof={D.zero:['zero']}
    for j,i in enumerate(Y):
        v=[D.zero]*9;v[j]=D.one;membership[i]=tuple(v);proof[i]=['coordinate',j]
    telescopes=[]
    for root in roots:
        previous=D.zero;coeff=[];division=[]
        for j in range(9):
            current=prefix[j][root] if j<8 else root
            numerator=D.add(current,D.neg(previous))
            h=D.put(('r',numerator,Y[j],f'H10_CA_source{root}:coordinate{j}'),derivatives[root][j])
            coeff.append(h);division.append({'node':h,'numerator':numerator,'coordinate':Y[j],
                                            'current_source_restriction':current,'previous_source_restriction':previous,
                                            'origin':str(derivatives[root][j])})
            previous=current
        assert previous==root
        membership[root]=tuple(coeff);proof[root]=['pure_CA_unit57_telescope',len(telescopes)]
        telescopes.append({'original_source_function':root,'membership_vector9':coeff,'all9_regular_coordinate_divisions':division,
                           'exact_identity':'sum Yj*hj telescopesfromconstant0tooriginalroot; prefix9=original bypointedunituniqueness'})
    # The original DAG simplifies -(-b) to its child, so recognize both signs
    # of every certified source function before traversing the relative inputs.
    for root in list(membership):
        node=D.nodes[root]
        if node[0]=='n' and node[1] not in membership:
            membership[node[1]]=tuple(D.neg(v) for v in membership[root])
            proof[node[1]]=['negative_of_certified_function',root]
    # Propagate ideal membership, treating all old r/g functions as scalar cofactors.
    # No pullback or derivative of those functions is needed anywhere in this step.
    def get(i):
        if i in membership:return membership[i]
        n=D.nodes[i];op=n[0];ans=None;why=None
        if op=='n':
            a=get(n[1])
            if a is not None:ans=tuple(D.neg(v) for v in a);why=['neg',n[1]]
        elif op=='a':
            a,b=get(n[1]),get(n[2])
            if a is not None and b is not None:ans=tuple(D.add(u,v) for u,v in zip(a,b));why=['sum',n[1],n[2]]
        elif op=='m':
            a=get(n[1])
            if a is not None:ans=tuple(D.mul(n[2],v) for v in a);why=['product_left_known_ideal_factor',n[1],n[2]]
            else:
                b=get(n[2])
                if b is not None:ans=tuple(D.mul(n[1],v) for v in b);why=['product_right_known_ideal_factor',n[1],n[2]]
        membership[i]=ans
        if why:proof[i]=why
        return ans
    def lift(pair):
        a,b=pair;v=get(b)
        assert v is not None,('unrecognized_odd_membership_factor',b,D.nodes[b])
        return [a]+list(v)
    Q3=full['actual_Qn_R0_squared_generator_matrices']['3'];tops=full['actual12_top_T3_additional_columns'];Jinf=full['actual_Jinf_R0_squared_generator_matrix']
    qlift=[lift(p) for p in Q3];tlift=[lift(p) for p in tops];jlift=[lift(p) for p in Jinf]
    assert len(qlift)==2034 and len(tlift)==12 and len(jlift)==2046
    koszul=[]
    for i in range(9):
        for j in range(i+1,9):
            v=[D.zero]*10;v[1+i]=D.neg(Y[j]);v[1+j]=Y[i];koszul.append(v)
    assert len(koszul)==36
    cert={'DAG_scalar_nodes':D.nodes,
          'actual_R0_coordinate_nodes_Y9':Y,'H10_generator_images_in_actual_CA_trace_pairs':[[D.one,D.zero]]+[[D.zero,i] for i in Y],
          'all8_pure_CA_prefix57_unit_systems':systems,'all_pure_source_Hadamard_membership_telescopes':telescopes,
          'old_r_g_function_pullbacks_performed':False,
          'odd_membership_propagation_rules':{str(i):v for i,v in proof.items()},
          'all36_Koszul_kernel_columns_R0_10':koszul,
          'all2034_Q3_columns_lifted_R0_10':qlift,'all12_top_columns_lifted_R0_10':tlift,
          'all2046_Jinf_columns_lifted_R0_10':jlift,
          'honest_D3_presentation_relation_columns_R0_10':koszul+qlift,
          'honest_D4_presentation_relation_columns_R0_10':koszul+qlift+tlift,
          'presentation_proof':'ker(R0^10->R0+ m9 tau)=36Koszul; lift mapsrecoveroriginalpair byrecursiveexactidentities',
          'actual_T3':'imageof12liftedtopcolumns incoker(honest_D3_relation_matrix)',
          'T3_presentation_not_yet_computed':'needsrelationsamong12images, i.e syzygiesagainst2070D3relationcolumns',
          'no_mod_m_squared_of_ambient_matrix_or_module_membership_claim':True}
    encoded=json.dumps(cert,separators=(',',':'))+'\n';assert len(encoded.encode())<=32*1024*1024
    (out/'honest10-presentation-DAG.json').write_text(encoded)
    paths=[Path(__file__),source,fullp,originp,graphp,lengths,fullpeer,weightpath]
    report={'id':STEM,'all_passed':True,'status':'literal_honest10_D3_D4_input_not_T3_presentation_or_depth_solution',
            'actual_H10_presentation_kernel_exact_Koszul36':True,
            'all_Q3_top_and_Jinf_odd_coefficients_lifted_by_exact_rules':True,
            'D3_free10_relation_columns':len(koszul+qlift),'D4_free10_relation_columns':len(koszul+qlift+tlift),
            'new_pure_CA_unit57_prefix_systems':8,'new_regular_coordinate_divisions':9*len(roots),
            'old_nonunit_r_g_functions_never_restricted_or_redivided':True,
            'actual_T3_relations_nonzero_colon_depth_or_nil_solved':False,
            'no_GB_new_jets_points_or_ungauged_core':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{'honest10-presentation-DAG.json':sha(out/'honest10-presentation-DAG.json')},
            'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','D3_free10_relation_columns','D4_free10_relation_columns','new_regular_coordinate_divisions','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
