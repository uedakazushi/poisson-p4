#!/usr/bin/env python3
"""ONE exact global-free4 tangency test and literal77/linear77 nil extractor."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_top_graph_telescoping_20261104_v63 import forward
from ambient_W1_S61_intersection_recognition_review_20261103_v63 import inverse

ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    import hashlib
    return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(A):return [[str(x) for x in row] for row in A]

def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('ONE nil-extractor60s cap; no automatic larger retry')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    pref=ROOT/'reports/AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_PREFLIGHT_20261104_v63';pre=read(pref/'preflight.json')
    bp=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json';oldp=bp.parent/'report.json'
    builder,old=map(read,[bp,oldp]);names=builder['ambient_variables'];vs=s.symbols(' '.join(names));index={n:i for i,n in enumerate(names)};local=dict(zip(names,vs))
    free=['v1','v4','v5','v7'];freeidx=[index[n] for n in free];annidx=[index[n] for n in old['annihilator_prime_generators']]
    assert len(annidx)==45 and set(annidx)==set(range(49))-set(freeidx)
    assert pre['candidate_nonzero_coordinates']=={'v0':'1'} and pre['dF_on_candidate']=='1'
    # Literal polynomial restrictions over QQ[free4], never just a center test.
    identities=[];equations=builder['component_presentations']['CA']+builder['component_presentations']['CB']
    for row,text in enumerate(equations):
        polynomial=s.Poly(s.sympify(text,locals=local),*vs);zero={};tangent={};contributions=0
        for mon,c in polynomial.terms():
            ex=tuple(mon[i] for i in freeidx)
            if not any(mon[i] for i in annidx):zero[ex]=zero.get(ex,s.S.Zero)+c
            if mon[index['v0']]==1 and not any(mon[i] for i in annidx if i!=index['v0']):
                tangent[ex]=tangent.get(ex,s.S.Zero)+c;contributions+=1
        zero={m:str(c) for m,c in zero.items() if c};tangent={m:str(c) for m,c in tangent.items() if c}
        if zero or tangent:
            counter={'all_passed':False,'counterexample_row':row,'equation':text,'zero_restriction':[[list(m),c] for m,c in zero.items()],
                'tangent_restriction':[[list(m),c] for m,c in tangent.items()],'stopped_without_changing_ansatz':True}
            (out/'counterexample.json').write_text(json.dumps(counter,indent=2)+'\n');print(json.dumps(counter));return
        identities.append({'row':row,'source_component':'CA' if row<len(builder['component_presentations']['CA']) else 'CB',
            'polynomial_terms':len(polynomial.terms()),'selected_linear_v0_terms':contributions,'constant_and_dual_coefficients_exactly_zero':True})
    assert len(identities)==218 and old['nilradical_generator']=='2*aa + v0 + z25'
    Fpoly=s.Poly(s.sympify(old['nilradical_generator'],locals=local),*vs)
    assert s.diff(Fpoly.as_expr(),vs[index['v0']])==1
    # Fresh finite input; immutable parentDAGs are read, never edited.
    tailp=ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63/literal-tail-binding-DAG.json';tail=read(tailp)
    np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json';nb=read(np)
    cp=ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json';nc=read(cp)
    gates=[pref/'report.json',ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json']
    for rp in gates:
        r=read(rp);assert r['all_passed']
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
    D=DAG();g=[D.var(f'NIL_CURVE_GROUP_{i}') for i in range(17)];f=[D.var('NIL_CURVE_'+n) for n in free]
    eta=[D.var(f'NIL_CURVE_NCOEF_{i}') for i in range(56)];unknown=g+f+eta;assert len(unknown)==77
    epsilon=D.var('NIL_AUX_V0');Z=D.var('NIL_CURVE_Z_Z1')
    bind={i:D.zero for i in tail['old49_centered_variable_nodes'].values()}
    bind.update(zip([tail['old49_centered_variable_nodes'][n] for n in free],f));bind[tail['old49_centered_variable_nodes']['v0']]=epsilon
    bind.update(zip(tail['actual_group17_parameter_nodes'],g));tcache={}
    def ts(i):
        if i in tcache:return tcache[i]
        n=tail['DAG_scalar_nodes'][i]
        if n[0]=='v':v=bind[i]
        elif n[0]=='c':v=D.const(F(n[1],n[2]))
        elif n[0]=='n':v=D.neg(ts(n[1]))
        elif n[0]=='i':v=D.inverse(ts(n[1]))
        else:v={'a':D.add,'m':D.mul}[n[0]](ts(n[1]),ts(n[2]))
        tcache[i]=v;return v
    cut20=[ts(i) for i in tail['all20_common_cut_equation_nodes']]
    rho=[ts(i) for i in tail['same_physical_R0_readback_nodes9']];physical=[ts(i) for i in tail['actual_full_physical150_nodes']]
    Rrows=tail['same_physical_R0_rows9'];assert Rrows==nb['actual_same_R0_rows9']==nc['same_actual_R0_rows9']
    cohenrows=nc['actual_R8_physical_rows8'];assert cohenrows==[108,109,112,113,117,118,119,27]
    nbind=dict(zip(nb['all56_coefficient_unknown_nodes'],eta))
    nbind.update((node,rho[Rrows.index(row)]) for node,row in zip(nb['actual_Cohen_R8_parameter_nodes'],cohenrows))
    nbind[nc['actual_independent_R0_rho61_node']]=rho[Rrows.index(61)]
    ncache={}
    def ns(i):
        if i in ncache:return ncache[i]
        n=nc['DAG_scalar_nodes'][i]
        if n[0]=='v':assert i in nbind,(i,n);v=nbind[i]
        elif n[0]=='c':v=D.const(F(n[1],n[2]))
        elif n[0]=='n':v=D.neg(ns(n[1]))
        elif n[0]=='i':v=D.inverse(ns(n[1]))
        else:v={'a':D.add,'m':D.mul}[n[0]](ns(n[1]),ns(n[2]))
        ncache[i]=v;return v
    n56=[ns(i) for i in nb['all56_actual_coefficient_Hensel_equations']]
    z=[ns(i) for i in nc['z1_through_z8_actual_odd_coefficient_nodes']];nu=ns(nc['nu_actual_node']);a61=ns(nc['actual_a61_node']);vnode=ns(nc['actual_v_node'])
    equations77=cut20+n56+[D.add(z[0],D.neg(Z))]
    assert len(equations77)==77 and all(D.values[i]==0 for i in equations77)
    values,grad=forward(D.nodes,set(unknown+[epsilon]))
    A=[[grad[i].get(v,F(0)) for v in unknown] for i in equations77]
    forcing=[grad[i].get(epsilon,F(0)) for i in equations77]
    Ai,adet=inverse(A);assert adet==F(-47775744000)==-80*298598400*2
    center=[-sum(a*b for a,b in zip(row,forcing)) for row in Ai]
    assert all(sum(a*b for a,b in zip(row,center))+b==0 for row,b in zip(A,forcing))
    # Actual full Nmodule lifts and coefficients, all in the SAME physical R0.
    n_lifts=[D.one]+[D.add(physical[row],D.neg(ns(rec['actual_even_part_node']))) for row,rec in zip(nc['z1_through_z8_physical_target_rows'][1:],nc['eight_N_module_generators'][1:])]
    assert len(n_lifts)==8
    d1=[{int(i):ns(node) for i,node in col.items()} for col in nc['R0_sparse_d1_8_by36_columns']]
    d2=[{int(i):ns(node) for i,node in col.items()} for col in nc['R0_sparse_d2_36_by84_columns']]
    relations36=[D.total(D.mul(node,n_lifts[i]) for i,node in col.items()) for col in d1]
    # Original scalar77-system is used ONLY to take its first-order epsilon map.
    # It is not claimed to parametrize actual D for finite epsilon.
    original_nodes=D.nodes[:]
    dualvars=[D.var(f'NIL_DUAL_CENTERED_{i}') for i in range(77)]
    seed={v:(v,D.add(du,D.const(c))) for v,du,c in zip(unknown,dualvars,center)}
    seed[epsilon]=(D.zero,D.one);seed[Z]=(Z,D.zero)
    paircache={}
    def dual(i):
        if i in paircache:return paircache[i]
        n=original_nodes[i];op=n[0]
        if op=='v':assert i in seed,(i,n);p=seed[i]
        elif op=='c':p=(D.const(F(n[1],n[2])),D.zero)
        elif op=='a':x,y=dual(n[1]),dual(n[2]);p=(D.add(x[0],y[0]),D.add(x[1],y[1]))
        elif op=='n':x=dual(n[1]);p=(D.neg(x[0]),D.neg(x[1]))
        elif op=='m':
            x,y=dual(n[1]),dual(n[2]);p=(D.mul(x[0],y[0]),D.add(D.mul(x[0],y[1]),D.mul(x[1],y[0])))
        else:
            assert op=='i';x=dual(n[1]);iv=D.inverse(x[0]);p=(iv,D.neg(D.mul(x[1],D.mul(iv,iv))))
        paircache[i]=p;return p
    e77=[dual(i) for i in equations77];zero77=[p[0] for p in e77];linear77=[p[1] for p in e77]
    physical_pairs=[dual(i) for i in physical];R0pairs=[dual(i) for i in rho]
    Npairs=[dual(i) for i in n_lifts];z_pairs=[dual(i) for i in z];nu_pair=dual(nu);v_pair=dual(vnode)
    relation_pairs=[dual(i) for i in relations36];c36=[p[1] for p in relation_pairs]
    psi8=[p[0] for p in z_pairs[1:]]+[v_pair[0]]
    d1_pairs=[{str(i):dual(node) for i,node in col.items()} for col in d1]
    d2_curve=[{str(i):dual(node)[0] for i,node in col.items()} for col in d2]
    # Exactly linear in dual77, all inverses independent of them.
    degree=[];dualset=set(dualvars)
    for i,n in enumerate(D.nodes):
        if n[0]=='v':v=1 if i in dualset else 0
        elif n[0]=='c':v=0
        elif n[0] in ['n','i']:
            v=degree[n[1]]
            if n[0]=='i':assert v==0
        elif n[0]=='a':v=max(degree[n[1]],degree[n[2]])
        else:assert n[0]=='m';v=degree[n[1]]+degree[n[2]]
        assert v<=1;degree.append(v)
    assert all(degree[i]==0 for i in zero77+psi8) and all(degree[i]<=1 for i in linear77+c36)
    assert all(D.values[i]==0 for i in zero77+linear77)
    _,dg=forward(D.nodes,set(unknown+dualvars))
    A0=[[dg[i].get(v,F(0)) for v in unknown] for i in zero77]
    Ad=[[dg[i].get(v,F(0)) for v in dualvars] for i in linear77]
    assert A0==Ad==A
    # The original old49 source map has F0=0 and Fdelta=1 literally.
    oldpairs={}
    for n,oldnode in tail['old49_centered_variable_nodes'].items():
        p=dual(bind[oldnode])
        if n in ['v4','v7']:p=(D.add(p[0],D.one),p[1])
        oldpairs[n]=p
    assert oldpairs['v0']==(D.zero,D.one) and oldpairs['aa']==oldpairs['z25']==(D.zero,D.zero)
    assert all(D.values[p[0]]==0 for p in relation_pairs)
    coefficient={'DAG_scalar_nodes':D.nodes,
        'actual_old49_unshifted_dual_map':oldpairs,'actual_group17_dual_map':[dual(i) for i in g],
        'old_nil_F_expression':old['nilradical_generator'],'old_nil_F_dual_pair':[D.zero,D.one],
        'curve_parameter_Z_node':Z,'curve_parameter_is_actual_z1_b61':True,
        'curve77_unknowns_group17_free4_N56':unknown,'curve77_equations':zero77,
        'curve77_equation_order':'oldphysical20cuts,actualN56coefficients,z1-Z',
        'curve77_origin_Jacobian':enc(A),'curve77_unit_determinant':str(adet),
        'linear_dual77_centered_unknowns':dualvars,'linear_dual77_unknown_center_shifts':list(map(str,center)),
        'linear_dual77_equations':linear77,'dual77_origin_Jacobian_equals_curve77':True,
        'linear_dual_seed_is_actual_old_v0':True,'linear_dual_all_division_inputs_have_dual_degree0':True,
        'actual_physical150_dual_pairs':physical_pairs,
        'same_actual_R0_rows9':Rrows,'same_actual_R0_dual_pairs':R0pairs,
        'actual_R8_physical_rows8':cohenrows,
        'actual_N_z1_through_z8_dual_pairs':z_pairs,'actual_N_nu_dual_pair':nu_pair,'actual_R0_v_dual_pair':v_pair,
        'actual_N_module_generator_lifts8_dual_pairs':Npairs,
        'N_d1_actual36_dual_columns':d1_pairs,'N_d2_actual84_curve_columns':d2_curve,
        'actual36_N_relation_on_D_lifts_dual_pairs':relation_pairs,
        'positive_nil_residue_coefficients_c36':c36,
        'positive_sign_convention':'Nrelation(lift)=c_l(Z)*F_old;horseshoe_d1_upperright=-c',
        'nilcurve_graph_ideal_P8_variables':['z2','z3','z4','z5','z6','z7','z8','v'],
        'nilcurve_graph_ideal_P8_right_sides_psi':psi8,
        'nilcurve_graph_ideal_formula':'P=(z2-psi2(z1),...,z8-psi8(z1),v-psiv(z1)),z1=Z',
        'all77_and_linear77_are_pointed_unit_inputs_not_expanded_series':True,
        'ringmap':'actualDtail -> QQ[[Z]][delta]/delta^2, with curve77 and linear-dual77 giving unique coefficient functions',
        'nil_module_extractor':'delta(a*F_old)=a|curve; R0-linear on Nil only',
        'whole_ring_delta_not_assumed_R0_linear':True,
        'all9_R0_delta_components_retained_for_honest_coefficient_corrections':True,
        'Hadamard_or_regular_coordinate_quotient_nodes_not_composed':True,
        'no_36_coefficient_series_expansion_or_depth_T3_solution':True,
        'all_unit_inverse_requests':D.divisions}
    text=json.dumps(coefficient,separators=(',',':'))+'\n';assert len(text.encode())<30*1024*1024
    cert=out/'nil-delta-extractor-DAG.json';cert.write_text(text)
    proof={'all218_QQfree4_identity_records':identities,'free4_unshifted_coordinate_order':free,
        'source_preflight_derivation_is_partial_v0':True,'actual_dF_equals1':True,
        'curve77_determinant_block_factors':['-80','298598400','2'],'curve77_determinant':str(adet),
        'all_curve77_and_linear77_center_values_zero':True,
        'actual36_nil_residue_constant_terms':[str(D.values[i]) for i in c36],
        'constant_terms_not_promoted_to_full_curve_function_vanishing':True}
    proofpath=out/'tangency-and-unit-proof.json';proofpath.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),pref/'report.json',pref/'preflight.json',bp,oldp,tailp,np,cp,*gates,
        ROOT/'code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py',ROOT/'code/ambient_B18_top_graph_telescoping_20261104_v63.py']
    report={'id':ID,'all_passed':True,'all218_partial_v0_restrictions_identically_zero_over_QQfree4':True,
        'actual_old_nil_F_maps_exactly_delta':True,'actual_same_R0_curve77_and_linear_dual77_constructed':True,
        'curve77_unit_determinant':str(adet),'all36_positive_nil_residue_functions_specified_implicitly':True,
        'all9_actual_R0_dual_readbacks_and_psi8_kept':True,
        'no_full_curve_series_Hadamard_composition_huge_inverse_GB_or_jet':True,
        'T3_zero_or_B18_origin_depth_or_extension_split_claim':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'evidence_sha256':{cert.name:sha(cert),proofpath.name:sha(proofpath)},
        'DAG_nodes':len(D.nodes),'artifact_bytes':sum(p.stat().st_size for p in out.iterdir()),
        'seconds':time.monotonic()-began,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True,
        'independent_review_pending':True}
    assert report['artifact_bytes']<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','curve77_unit_determinant','DAG_nodes','artifact_bytes','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except (TimeoutError,MemoryError) as e:
        a.output.mkdir(parents=True,exist_ok=True)
        (a.output/'bounded-failure.json').write_text(json.dumps({'all_passed':False,'reason':str(e),'no_automatic_larger_retry':True})+'\n');raise
