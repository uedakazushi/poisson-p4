#!/usr/bin/env python3
"""Independent exact tangency, fixed circuit routing and dual77 unit audit."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview

ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_20261104_v63'
ID=BASE+'_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def solve(A,b):
    """Own exact augmented elimination; no author matrix routine."""
    n=len(A);a=[list(row)+[v] for row,v in zip(A,b)];d=F(1)
    for j in range(n):
        k=next(k for k in range(j,n) if a[k][j])
        if k!=j:a[k],a[j]=a[j],a[k];d=-d
        p=a[j][j];d*=p;a[j]=[v/p for v in a[j]]
        for i in range(j+1,n):
            q=a[i][j]
            if q:a[i]=[v-q*w for v,w in zip(a[i],a[j])]
    x=[F(0)]*n
    for i in reversed(range(n)):x[i]=a[i][-1]-sum(a[i][j]*x[j] for j in range(i+1,n))
    assert all(sum(a*b for a,b in zip(row,x))==v for row,v in zip(A,b))
    return x,d

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('Independent extractor review60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;rp=folder/'report.json';primary=read(rp);assert primary['all_passed']
    checked={}
    for key in ['input_sha256','evidence_sha256']:
        for name,want in primary[key].items():
            p=ROOT/name if key=='input_sha256' else folder/name
            assert sha(p)==want,name;checked[str(p.relative_to(ROOT))]=want
    bp=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json'
    oldp=bp.parent/'report.json';builder,old=read(bp),read(oldp)
    names=builder['ambient_variables'];variables=s.symbols(' '.join(names));loc=dict(zip(names,variables))
    free=['v1','v4','v5','v7'];ann=set(old['annihilator_prime_generators'])
    assert set(names)-ann==set(free) and len(ann)==45
    # Restrict using symbolic substitution and differentiation, independently
    # of the producer's selected-monomial implementation.
    zero_sub={loc[n]:0 for n in ann};records=[]
    equations=builder['component_presentations']['CA']+builder['component_presentations']['CB']
    for j,text in enumerate(equations):
        f=s.sympify(text,locals=loc)
        a=s.Poly(f.subs(zero_sub,simultaneous=True),*(loc[n] for n in free))
        b=s.Poly(s.diff(f,loc['v0']).subs(zero_sub,simultaneous=True),*(loc[n] for n in free))
        assert a.is_zero and b.is_zero,j
        records.append({'row':j,'value_and_partial_v0_identically_zero':True})
    assert len(records)==218
    f=s.sympify(old['nilradical_generator'],locals=loc)
    assert f.subs(zero_sub,simultaneous=True)==0
    assert s.diff(f,loc['v0']).subs(zero_sub,simultaneous=True)==1
    cp=folder/'nil-delta-extractor-DAG.json';c=read(cp);D=DAGReview(c['DAG_scalar_nodes'])
    vn=D.variables;unk=c['curve77_unknowns_group17_free4_N56'];du=c['linear_dual77_centered_unknowns']
    assert unk==[vn[f'NIL_CURVE_GROUP_{i}'] for i in range(17)]+[vn['NIL_CURVE_'+n] for n in free]+[vn[f'NIL_CURVE_NCOEF_{i}'] for i in range(56)]
    assert du==[vn[f'NIL_DUAL_CENTERED_{i}'] for i in range(77)]
    epsilon=vn['NIL_AUX_V0'];Z=vn['NIL_CURVE_Z_Z1'];assert c['curve_parameter_Z_node']==Z

    def transform(nodes,roots,seeds,pair=False):
        memo=dict(seeds)
        def visit(i):
            if i in memo:return memo[i]
            n=nodes[i];op=n[0]
            if op=='c':v=D.const(F(n[1],n[2]));v=(v,D.zero) if pair else v
            else:
                assert op!='v',('unbound actual input',i,n)
                a=visit(n[1]);b=visit(n[2]) if op in ['a','m'] else None
                if not pair:
                    v={'a':lambda:D.add(a,b),'m':lambda:D.mul(a,b),'n':lambda:D.neg(a),'i':lambda:D.inverse(a)}[op]()
                elif op=='a':v=(D.add(a[0],b[0]),D.add(a[1],b[1]))
                elif op=='n':v=(D.neg(a[0]),D.neg(a[1]))
                elif op=='m':v=(D.mul(a[0],b[0]),D.add(D.mul(a[0],b[1]),D.mul(a[1],b[0])))
                else:
                    iv=D.inverse(a[0]);v=(iv,D.neg(D.mul(a[1],D.mul(iv,iv))))
            memo[i]=v;return v
        return [visit(i) for i in roots]

    tp=ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63/literal-tail-binding-DAG.json';tail=read(tp)
    np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json';nb=read(np)
    ncpath=ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json';nc=read(ncpath)
    rows=tail['same_physical_R0_rows9'];assert rows==c['same_actual_R0_rows9']==nc['same_actual_R0_rows9']
    assert nc['actual_R8_physical_rows8']==[108,109,112,113,117,118,119,27]
    seed={v:D.zero for v in tail['old49_centered_variable_nodes'].values()}
    seed.update((tail['old49_centered_variable_nodes'][n],v) for n,v in zip(free,unk[17:21]))
    seed[tail['old49_centered_variable_nodes']['v0']]=epsilon
    seed.update(zip(tail['actual_group17_parameter_nodes'],unk[:17]))
    tailroots=tail['all20_common_cut_equation_nodes']+tail['same_physical_R0_readback_nodes9']+tail['actual_full_physical150_nodes']
    vals=transform(tail['DAG_scalar_nodes'],tailroots,seed);cut,rho,physical=vals[:20],vals[20:29],vals[29:]
    nseed=dict(zip(nb['all56_coefficient_unknown_nodes'],unk[21:]))
    nseed.update((v,rho[rows.index(row)]) for v,row in zip(nb['actual_Cohen_R8_parameter_nodes'],nc['actual_R8_physical_rows8']))
    nseed[nc['actual_independent_R0_rho61_node']]=rho[rows.index(61)]
    ev=[q['actual_even_part_node'] for q in nc['eight_N_module_generators'][1:]]
    nroots=nb['all56_actual_coefficient_Hensel_equations']+nc['z1_through_z8_actual_odd_coefficient_nodes']+[nc['nu_actual_node'],nc['actual_v_node']]+ev
    nv=transform(nc['DAG_scalar_nodes'],nroots,nseed)
    eq77=cut+nv[:56]+[D.add(nv[56],D.neg(Z))]
    A=[[D.grad[i].get(v,F(0)) for v in unk] for i in eq77]
    rhs=[-D.grad[i].get(epsilon,F(0)) for i in eq77]
    shift,det=solve(A,rhs);assert det==F(-47775744000)
    assert list(map(F,c['linear_dual77_unknown_center_shifts']))==shift
    assert [[str(v) for v in row] for row in A]==c['curve77_origin_Jacobian']
    # Route the independently bound actual scalar functions through dual
    # arithmetic. Every equality below is an exact circuit-node equality,
    # using ring/dual identities, not equality at an evaluation point.
    dualseed={v:(v,D.add(w,D.const(q))) for v,w,q in zip(unk,du,shift)}
    dualseed[epsilon]=(D.zero,D.one);dualseed[Z]=(Z,D.zero)
    original_roots=eq77+physical+rho+nv[56:]
    pp=transform(D.nodes,original_roots,dualseed,True)
    e77=pp[:77];p150=pp[77:227];r9=pp[227:236];nz=pp[236:244];nnu,nvpair=pp[244:246];even=pp[246:]
    assert [p[0] for p in e77]==c['curve77_equations']
    assert [p[1] for p in e77]==c['linear_dual77_equations']
    assert p150==list(map(tuple,c['actual_physical150_dual_pairs']))
    assert r9==list(map(tuple,c['same_actual_R0_dual_pairs']))
    assert nz==list(map(tuple,c['actual_N_z1_through_z8_dual_pairs']))
    assert nnu==tuple(c['actual_N_nu_dual_pair']) and nvpair==tuple(c['actual_R0_v_dual_pair'])
    assert [p[0] for p in nz[1:]]+[nvpair[0]]==c['nilcurve_graph_ideal_P8_right_sides_psi']
    n_lifts=[(D.one,D.zero)]+[(D.add(p150[row][0],D.neg(a[0])),D.add(p150[row][1],D.neg(a[1]))) for row,a in zip(nc['z1_through_z8_physical_target_rows'][1:],even)]
    assert n_lifts==list(map(tuple,c['actual_N_module_generator_lifts8_dual_pairs']))
    d1roots=[v for col in nc['R0_sparse_d1_8_by36_columns'] for v in col.values()]
    nscalar=transform(nc['DAG_scalar_nodes'],d1roots,nseed)
    d1pairs=iter(transform(D.nodes,nscalar,dualseed,True));relations=[]
    for col,want in zip(nc['R0_sparse_d1_8_by36_columns'],c['N_d1_actual36_dual_columns']):
        terms=[]
        for key in col:
            a=next(d1pairs);assert a==tuple(want[key]);b=n_lifts[int(key)]
            terms.append((D.mul(a[0],b[0]),D.add(D.mul(a[0],b[1]),D.mul(a[1],b[0]))))
        r=(D.total(p[0] for p in terms),D.total(p[1] for p in terms));relations.append(r)
    assert relations==list(map(tuple,c['actual36_N_relation_on_D_lifts_dual_pairs']))
    assert [p[1] for p in relations]==c['positive_nil_residue_coefficients_c36']
    for i in c['curve77_equations']+c['linear_dual77_equations']:assert D.values[i]==0
    A0=[[D.grad[i].get(v,F(0)) for v in unk] for i in c['curve77_equations']]
    A1=[[D.grad[i].get(v,F(0)) for v in du] for i in c['linear_dual77_equations']]
    assert A0==A1==A
    degree=[]
    for i,n in enumerate(D.nodes):
        op=n[0]
        if op=='v':v=int(i in du)
        elif op=='c':v=0
        elif op in ['n','i']:
            v=degree[n[1]]
            if op=='i':assert not v
        else:v=max(degree[n[1]],degree[n[2]]) if op=='a' else degree[n[1]]+degree[n[2]]
        assert v<=1;degree.append(v)
    proof={'all218_exact_free4_tangency_records':records,'actual_same_R0_rows':rows,
           'curve77_and_dual77_origin_matrix':[[str(v) for v in row] for row in A],
           'unit_determinant':str(det),'independent_dual_center_shifts':list(map(str,shift)),
           'all150_plus9_physical_dual_functions_and36_residues_exactly_routed':True,
           'whole_ring_delta_R0_linearity_assumed':False,
           'only_nil_module_extractor_is_R0_linear':True,
           'all_inverse_inputs_have_dual_degree_zero':True,
           'c36_series_vanishing_depth_or_T0_answer_claimed':False}
    ep=out/'independent-tangency-unit-and-routing.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    sources=[Path(__file__),ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',
             ROOT/'research/notes'/f'{ID}.md',rp,cp,bp,oldp,tp,np,ncpath]
    report={'id':ID,'all_passed':True,'all218_free4_tangencies_exact':True,
            'all77_plus_dual77_actual_circuit_bindings_independently_rebuilt':True,
            'all150_physical_and9_R0_dual_pairs_and36_residues_rebuilt':True,
            'both77_Jacobians_equal_and_unit_determinant':str(det),
            'not_author_solver_rerun_or_new_source_jet':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sources},
            'checked_dependency_sha256':checked,'evidence_sha256':{ep.name:sha(ep)},
            'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items() if 'sha256' not in k}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
