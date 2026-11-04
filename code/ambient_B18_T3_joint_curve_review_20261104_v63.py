#!/usr/bin/env python3
"""Independent actual joint105 sparse microjets and symmetric annihilator audit."""
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
from ambient_B18_CB_g_division_origins_review_20261104_v63 import MON,plus,scale,evaluate,vector

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW'


def read(p):return json.loads(p.read_text())


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63'
    author=read(base/'report.json');cert=read(base/'joint-degree3-certificate.json')
    literal=read(base/'joint105-algebraic-DAG.json');assert author['all_passed']
    for rel,h in author['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in author['evidence_sha256'].items():assert sha(base/rel)==h,rel
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    sourcepeer=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW/report.json'
    capeer=ROOT/'reports/AMBIENT_B18_T3_CURVE_ANNIHILATOR_20261104_v63_REVIEW/report.json'
    weightedpeer=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW/report.json'
    fullpeer=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json'
    for p in(sourcepeer,capeer,weightedpeer,fullpeer):assert read(p)['all_passed']
    source=read(sourcepath);nodes=source['DAG_scalar_nodes'];ca,cb=[source['families'][f]for f in('CA','CB')]
    D=DAGReview(nodes)
    xids=[D.variables[f'x{r}']for r in source['physical_zero_base_rows9']]
    def coefficient_ids(f):return [i for pair in f['coefficient_unknown_pairs'].values()for i in pair]+[f['p_quadratic_relation'][k]for k in('s','t')]
    az,bz=coefficient_ids(ca),coefficient_ids(cb);variables=az+xids+bz
    assert len(variables)==len(set(variables))==105
    assert variables==cert['all105_unknown_original_node_order']
    agen=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']
    assert [r for r,i in agen]==[18,19,28,29,60,67,72,63,46,95]
    selected=[0,1,2,3,4,5,6,8,9]
    F9=[agen[i][1]for i in selected];g=agen[7][1]
    bnodes=[i for r,i in cb['actual_zero_odd_I0_generators']]
    brows=[r for r,i in cb['actual_zero_odd_I0_generators']]
    eq=ca['all48_zero_IFT_quadratic_coefficient_equations']+F9+cb['all48_zero_IFT_quadratic_coefficient_equations']
    assert eq==cert['all105_equation_original_node_order']
    needed=sorted(closure(nodes,eq+[g]+bnodes));assert len(needed)==17956
    assert {i for i in needed if nodes[i][0]=='v'}<=set(variables)
    J=[[D.grad[row].get(var,F(0))for var in variables]for row in eq]
    assert J==[list(map(F,row))for row in cert['fixed105_Jacobian']]
    inverse=[list(map(F,row))for row in cert['fixed105_Jacobian_inverse']]
    for A,B in((J,inverse),(inverse,J)):
        assert all(sum((A[i][k]*B[k][j]for k in range(105)if A[i][k]and B[k][j]),F(0))==int(i==j)
                   for i in range(105)for j in range(105))
    assert all(J[i][j]==0 for i in range(57)for j in range(57,105))
    _,detA=solve([row[:57]for row in J[:57]],[[]for _ in range(57)])
    _,detB=solve([row[57:]for row in J[57:]],[[]for _ in range(48)])
    assert detA==F(-3623878656,625)and detB==F(-331776,625)
    assert detA*detB==F(cert['det_joint105'])
    # Independent sparse dictionary jets/geometric-series unit inverses.
    targets=[{}for _ in range(55)]+[{(1,0):F(1)},{(0,1):F(1)}]+[{}for _ in range(48)]
    solution={i:{}for i in variables};stages=[]
    for degree in range(1,4):
        values=evaluate(nodes,needed,solution);mons=[m for m in MON if sum(m)==degree]
        residual=[plus(values[row],scale(target,-1))for row,target in zip(eq,targets)]
        corrections=[[-sum((inverse[i][k]*residual[k].get(mon,F(0))for k in range(105)if inverse[i][k]),F(0))
                      for mon in mons]for i in range(105)]
        saved=cert['homogeneous_corrections'][degree-1]
        assert saved['degree']==degree and saved['monomials']==list(map(list,mons))
        assert corrections==[list(map(F,row))for row in saved['all105_corrections']]
        for i,var in enumerate(variables):solution[var]=plus(solution[var],dict(zip(mons,corrections[i])))
        stages.append(dict(degree=degree,all105_corrections_match=True))
    values=evaluate(nodes,needed,solution)
    assert all(values[row]==target for row,target in zip(eq,targets))
    saved={i:list(map(F,v))for i,v in cert['all105_unknown_jets']}
    assert all(vector(solution[i])==saved[i]for i in variables)
    assert vector(values[g])==list(map(F,cert['g_CA_all10_coefficients']))
    for row,node,saved in zip(brows,bnodes,cert['all8_CB_zero_restrictions']):
        assert saved['CB_physical_zero_row']==row and saved['old_source_node']==node
        assert vector(values[node])==list(map(F,saved['all_degree_le3_coefficients']))
        assert all(not co for mon,co in values[node].items()if sum(mon)<2)
    assert all(not co for mon,co in values[g].items()if sum(mon)<3)
    # Rebuild the entire finite105 algebraic clone from original-source
    # nodes and supplied variable binding; do not trust only its jet output.
    ln=[tuple(n)for n in literal['DAG_scalar_nodes']]
    assert ln[:len(nodes)]==[tuple(n)for n in nodes]
    for n in ln[len(nodes):]:
        op=n[0]
        if op=='c':q=F(n[1],n[2])
        elif op=='v':q=F(0)
        elif op=='n':q=-D.values[n[1]]
        elif op=='a':q=D.values[n[1]]+D.values[n[2]]
        elif op=='m':q=D.values[n[1]]*D.values[n[2]]
        elif op=='i':assert D.values[n[1]];q=1/D.values[n[1]]
        else:raise ValueError(n)
        D.values.append(q)
    D.nodes=ln;D.cache={q:i for i,q in enumerate(ln)};assert len(D.cache)==len(ln)
    binding=dict(literal['all105_unknown_bindings']);assert set(binding)==set(variables)
    assert all(ln[new]==('v','JOINT_INV_'+nodes[old][1])for old,new in binding.items())
    mapping={}
    for i in needed:
        n=nodes[i];op=n[0]
        if op=='c':new=D.const(F(n[1],n[2]))
        elif op=='v':new=binding[i]
        elif op=='n':new=D.neg(mapping[n[1]])
        elif op=='a':new=D.add(mapping[n[1]],mapping[n[2]])
        elif op=='m':new=D.mul(mapping[n[1]],mapping[n[2]])
        elif op=='i':new=D.inverse(mapping[n[1]])
        else:raise ValueError(n)
        mapping[i]=new
    rp,sp=literal['formal_curve_parameter_nodes']
    cloned=[mapping[i]for i in eq]
    cloned[55]=D.add(cloned[55],D.neg(rp));cloned[56]=D.add(cloned[56],D.neg(sp))
    assert cloned==literal['all105_pointed_unit_equations']
    assert mapping[g]==literal['actual_CA_curve_g_node']
    assert [[r,mapping[i]]for r,i in zip(brows,bnodes)]==literal['actual_CB_zero8_restriction_nodes']
    actual_roots=cloned+[mapping[g]]+[mapping[i]for i in bnodes]
    needclone=sorted(closure(ln,actual_roots))
    clonebindings={new:solution[old]for old,new in binding.items()}
    clonebindings[rp]={(1,0):F(1)};clonebindings[sp]={(0,1):F(1)}
    assert {i for i in needclone if ln[i][0]=='v'}<=set(clonebindings)
    clonevalues=evaluate(ln,needclone,clonebindings)
    assert all(not clonevalues[i]for i in cloned)
    assert clonevalues[mapping[g]]==values[g]
    assert all(clonevalues[mapping[i]]==values[i]for i in bnodes)
    quadratic=[]
    for row in(18,60,67):
        node=bnodes[brows.index(row)]
        quadratic.append([values[node].get(mon,F(0))for mon in((2,0),(1,1),(0,2))])
    _,detQ=solve(quadratic,[[]for _ in range(3)]);assert detQ==F(625,768)
    proof=dict(actual105_source_DAG_prefix_and_every_cloned_equation_bound=True,
        actual_common_physical_zero_base_rows=source['physical_zero_base_rows9'],
        det_CA57=str(detA),det_CB48=str(detB),det_joint105=str(detA*detB),
        all105_unknown_sparse_jets=[[i,list(map(str,vector(solution[i])))]for i in variables],
        homogeneous_stages=stages,g_CA_10_coefficients=list(map(str,vector(values[g]))),
        actual_CB_zero8_restrictions=[dict(row=r,coefficients=list(map(str,vector(values[i]))))for r,i in zip(brows,bnodes)],
        quadratic_rows=[18,60,67],quadratic_matrix=[[str(q)for q in row]for row in quadratic],quadratic_minor=str(detQ),
        whole_ideal_Nakayama='J=(g_CA,B_CB,1,...,B_CB,8) is contained in m^2. Three actual restrictions span m^2/m^3. Thus m^2=J+m^3; Nakayama on m^2/J gives J=m^2 exactly.',
        joint_zero_ideal_quotient='R0/(I0_CA+I0_CB)=k[[r,s]]/(r,s)^2',length=3,Hilbert_function=[1,2],
        symmetric_annihilator='For either ambient component A_d=w^d H (d>=1), A_0=C0. Actual full99+87 generators have relative weights<=3. Thus Q_d stabilizes at d>=4, and ker(w:D3->D4)=Qinf/Q3 is generated by cubic relations times H/C0. For every zero odd f_i and positive coefficient j_e, f_i*(j_e*tau*c)=j_e*((f_i*tau)*c) is a degree3 C0 relation. Therefore both I0_CA and I0_CB kill the SAME intrinsic T3.',
        no_Tor_or_trace_basis_identification='This is a quotient-module identity on each actual ambient presentation, not recomputing an image after base change or equating tau_A with tau_B.',
        T3_finiteness='T3 is a finite module over the actual length3 joint zero-ideal quotient. With its12 known generators, length(T3)<=36.',
        no_T3_nonzero_membership_socle_depth_or_full_nil_claim=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'joint105-algebraic-DAG.json',base/'joint-degree3-certificate.json',
            ROOT/'code/ambient_B18_T3_joint_curve_restrictions_20261104_v63.py',
            ROOT/'research/notes/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63.md',
            ROOT/'research/notes/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63.md',
            ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/report.json',
            sourcepath,sourcepeer,capeer,weightedpeer,fullpeer,
            ROOT/'code/ambient_B18_CB_g_division_origins_review_20261104_v63.py',
            ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',
            ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    result=dict(id=ID,all_passed=True,
        independently_rebuilt_actual105_AD_Jacobian_block_determinants_and_inverse_products=True,
        independent_sparse_dict_jets_all105_equations_and_unknowns_through3=True,
        literal_finite105_DAG_and_every_binding_reconstructed_from_actual_source=True,
        actual_CB_three_quadratics_span_m_squared_minor='625/768',
        exact_joint_zero_ideal_quotient_length=3,
        symmetric_intrinsic_T3_annihilation_and_reverse_weight3_stability_accepted=True,
        T3_finite_length_proved=True,T3_length_upper=36,
        no_T3_nonzero_colon_membership_depth19_full_nil_or_matrix_fiber_inference=True,
        no_author_rerun_new_point_GB_or_degree_above3=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,frozen=True)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k]for k in('all_passed','exact_joint_zero_ideal_quotient_length','T3_finite_length_proved','seconds','max_rss_kib')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
