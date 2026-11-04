#!/usr/bin/env python3
"""Independent full9 first-jet audit and structural T3=0 application.

Re-evaluates SAVED source2 coefficients; no new Newton solve or higher jets.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure, solve

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63'
MON = [()] + [(i,) for i in range(9)] + list(itertools.combinations_with_replacement(range(9), 2))
POS = {m: i for i, m in enumerate(MON)}
ZERO55 = (F(0),)*55


def read(p): return json.loads(p.read_text())


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for c in iter(lambda: f.read(1048576), b''): h.update(c)
    return h.hexdigest()


def qmul(a, b):
    v = [a[0]*b[0]]
    v += [a[0]*b[1+i]+b[0]*a[1+i] for i in range(9)]
    for i, j in MON[10:]:
        c = a[0]*b[POS[(i,j)]]+b[0]*a[POS[(i,j)]]+a[1+i]*b[1+j]
        if i != j: c += a[1+j]*b[1+i]
        v.append(c)
    return tuple(v)


def qinv(a):
    assert a[0]
    v = [1/a[0]]+[-a[1+i]/a[0]**2 for i in range(9)]
    for i, j in MON[10:]:
        v.append(-a[POS[(i,j)]]/a[0]**2+(1 if i == j else 2)*a[1+i]*a[1+j]/a[0]**3)
    return tuple(v)


def eval2(nodes, needed, bindings):
    v = {}
    for i in needed:
        n = nodes[i];o = n[0]
        if o == 'c': q = (F(n[1], n[2]),)+(F(0),)*54
        elif o == 'v': q = bindings[i]
        elif o == 'a': q = tuple(a+b for a,b in zip(v[n[1]],v[n[2]]))
        elif o == 'n': q = tuple(-a for a in v[n[1]])
        elif o == 'm': q = qmul(v[n[1]],v[n[2]])
        elif o == 'i': q = qinv(v[n[1]])
        else: raise ValueError(o)
        v[i] = q
    return v


def linear_mul(a, b):
    return [a[0]*b[0]]+[a[0]*y+b[0]*x for x,y in zip(a[1:],b[1:])]


def insert_basis(basis, vector):
    v = {i:a for i,a in enumerate(vector) if a}
    for pivot, b in basis:
        if pivot in v:
            c = v[pivot]
            for i,a in b.items():
                q = v.get(i,F(0))-c*a
                if q:v[i]=q
                else:v.pop(i,None)
    if not v:return None
    p=min(v);c=v[p];basis.append((p,{i:a/c for i,a in v.items()}));return p


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63'
    fp=folder/'full9-firstjet-certificate.json';fr=folder/'report.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    rp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    op=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    proof_peers=[
        ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',
    ]
    assert read(fr)['all_passed'] and all(read(p)['all_passed'] for p in proof_peers)
    cert,source,full,honest,unit,origin=map(read,[fp,sp,rp,hp,up,op])
    assert cert['source2_monomials']==list(map(list,MON))
    oldnodes=source['DAG_scalar_nodes'];D=DAGReview(oldnodes)
    xids=[D.variables[f'x{i}'] for i in source['physical_zero_base_rows9']]
    source_jets={};function_jets={};metadata={}
    unit_determinants={}
    for name in ['CA','CB']:
        family=source['families'][name];saved=cert['both57_source2'][name]
        variables=[i for p in family['coefficient_unknown_pairs'].values() for i in p]
        variables += [family['p_quadratic_relation'][k] for k in ['s','t']]
        variables += xids
        gen=family['actual_zero_odd_I0_generators']+family['actual_positive_odd_J0_generators']
        Y=[gen[i][1] for i in [0,1,2,3,4,5,6,8,9]]
        equations=family['all48_zero_IFT_quadratic_coefficient_equations']+Y
        assert variables==saved['unknown57'] and equations==saved['equations57'] and Y==saved['F9_nodes']
        J=[[D.grad[i].get(j,F(0)) for j in variables] for i in equations]
        assert J==[list(map(F,row)) for row in saved['fixed_J57']]
        Ji,det=solve(J,[[F(i==j) for j in range(57)] for i in range(57)])
        assert det==F(saved['determinant'])
        assert Ji==[list(map(F,row)) for row in saved['fixed_J57_inverse']]
        sol={int(i):tuple(map(F,a)) for i,a in saved['all57_source2jets']}
        assert set(sol)==set(variables) and all(len(a)==55 and not a[0] for a in sol.values())
        roots=list(map(int,saved['pure_function2jets']))
        needed=sorted(closure(oldnodes,equations+roots))
        v=eval2(oldnodes,needed,sol)
        for j,node in enumerate(equations):
            expected=[F(0)]*55
            if j>=48:expected[1+j-48]=1
            assert v[node]==tuple(expected)
        assert all(v[i]==tuple(map(F,saved['pure_function2jets'][str(i)])) for i in roots)
        source_jets[name]=sol;function_jets[name]=v
        metadata[name]=(family,Y,gen[7][1])
        unit_determinants[name]=str(det)
    # Independent coordinate adapter from the same physical x9.
    XA=[[source_jets['CA'][i][1+j] for j in range(9)] for i in xids]
    XB=[[source_jets['CB'][i][1+j] for j in range(9)] for i in xids]
    XBi,_=solve(XB,[[F(i==j) for j in range(9)] for i in range(9)])
    adapter=[[sum(XBi[i][k]*XA[k][j] for k in range(9)) for j in range(9)] for i in range(9)]
    assert adapter==[list(map(F,row)) for row in cert['YB_to_YA_linear']]
    _,adet=solve(adapter,[[F(i==j) for j in range(9)] for i in range(9)])
    assert adet==F(1,4096)

    def hadamard(p,j,tail):
        return [p[1+j]]+[p[POS[tuple(sorted((j,k)))]] if k<=j or (tail and k>=7) else F(0) for k in range(9)]
    known={};r_records={}
    for rec in honest['all_pure_source_Hadamard_membership_telescopes']:
        p=function_jets['CA'][rec['original_source_function']]
        for j,q in enumerate(rec['all9_regular_coordinate_divisions']):
            a=hadamard(p,j,False);assert a[0]==F(q['origin'])
            assert a==list(map(F,cert['all792_CA_h_firstjets'][str(q['node'])]))
            known[q['node']]=a+[F(0)]*4
            r_records[q['node']]=a
    graph_weights={q['raw_row']:q['weight'] for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['kind']=='graph_CB_minus_CA'}
    cb,_,cbf63=metadata['CB']
    p_g=function_jets['CB'][cbf63]
    assert all(p_g[i]==0 for i,m in enumerate(MON) if all(j>=7 for j in m))
    for rec in full['all_coordinate_divisions']:
        node=rec['node'];label=full['DAG_scalar_nodes'][node][3]
        if label.startswith('g_split_'):a=[F(0)]*10
        else:
            head,coordinate=label.split(':');j=int(coordinate[1:])
            if head=='f63':root=cbf63
            else:
                row=int(head.replace('zero_graph','').replace('graph',''))
                root=cb['complete_normalized_physical_w_coefficients'][graph_weights[row]][row][1]
            a=hadamard(function_jets['CB'][root],j,coordinate.startswith('z'))
        a=[a[0]]+[sum(a[1+k]*adapter[k][j] for k in range(9)) for j in range(9)]
        assert a[0]==F(rec['origin_value'])
        assert a==list(map(F,cert['all784_CB_h_firstjets_in_YA'][str(node)]))
        known[node]=a+[F(0)]*4
    gnodes=[q['node'] for q in full['all4_known_regular_g_divisions']]
    for j,node in enumerate(gnodes):
        a=[F(0)]*14;a[10+j]=1;known[node]=a
    known[full['CB_g_restriction_node']]=[F(0)]*14

    binding={}
    for name in ['CA','CB']:
        for i,q in source_jets[name].items():
            lin=list(q[1:10])
            if name=='CB':lin=[sum(lin[k]*adapter[k][j] for k in range(9)) for j in range(9)]
            a=[F(0)]+lin+[F(0)]*4
            if i in binding:assert binding[i]==a
            binding[i]=a
    original=[]
    for n in oldnodes:
        o=n[0]
        if o=='c':a=[F(n[1],n[2])]+[F(0)]*13
        elif o=='v':a=binding[D.variables[n[1]]]
        elif o=='a':a=[x+y for x,y in zip(original[n[1]],original[n[2]])]
        elif o=='n':a=[-x for x in original[n[1]]]
        elif o=='m':a=linear_mul(original[n[1]],original[n[2]])
        elif o=='i':
            b=original[n[1]];assert b[0];a=[1/b[0]]+[-x/b[0]**2 for x in b[1:]]
        else:raise ValueError(o)
        original.append(a)
    nodes=unit['DAG_scalar_nodes'];relations=unit['all_reduced_D3_relation_columns'];top=unit['all12_transformed_top_columns']
    needed=set();stack=[i for c in relations+top for i in c]
    while stack:
        i=stack.pop()
        if i in needed:continue
        needed.add(i)
        if i<len(oldnodes) or i in known:continue
        o=nodes[i][0]
        if o in ['a','m']:stack+=nodes[i][1:]
        elif o in ['n','i']:stack.append(nodes[i][1])
        elif o!='c':raise ValueError((i,nodes[i]))
    values={}
    for i in sorted(needed):
        if i<len(original):a=original[i]
        elif i in known:a=known[i]
        else:
            n=nodes[i];o=n[0]
            if o=='c':a=[F(n[1],n[2])]+[F(0)]*13
            elif o=='a':a=[x+y for x,y in zip(values[n[1]],values[n[2]])]
            elif o=='n':a=[-x for x in values[n[1]]]
            elif o=='m':a=linear_mul(values[n[1]],values[n[2]])
            elif o=='i':
                b=values[n[1]];assert b[0];a=[1/b[0]]+[-x/b[0]**2 for x in b[1:]]
            else:raise ValueError(o)
        values[i]=a
    basis=[];pivots=[];matrix=[]
    for j,col in enumerate(relations):
        a=[x for i in col for x in values[i][:10]]
        assert a==list(map(F,cert['all2069_known_firstjet90_columns'][j]))
        assert all(values[i][10:]==[F(0)]*4 for i in col)
        assert all(a[10*k]==0 for k in range(9)) and a[:10]==[F(0)]*10
        if insert_basis(basis,a) is not None:pivots.append(j)
        matrix.append(a)
    assert len(basis)==43
    for j,col in enumerate(top):
        assert all(values[i]==[F(0)]*14 for i in col)
        assert not any(map(F,cert['all12_top_results'][j]['known_firstjet90']))
    # All opaque g derivative sensitivities have disappeared; no derivatives
    # of these four regular functions have been set to zero.
    assert all(not any(map(F,row)) for c in cert['all2069_four_g_sensitivities'] for row in c)

    # Independent all-order order bound, using only proved scalar origins.
    origins=list(map(F,origin['all339669_scalar_origin_values_evaluated']))
    order=[]
    for i,n in enumerate(honest['DAG_scalar_nodes'][:len(origins)]):
        o=n[0]
        if origins[i]:a=0
        elif o=='c':a=999
        elif o=='m':a=order[n[1]]+order[n[2]]
        elif o=='a':a=max(1,min(order[n[1]],order[n[2]]))
        elif o=='n':a=order[n[1]]
        elif o in ['v','r','g']:a=1
        else:raise ValueError((i,n))
        order.append(a)
    order_records=[]
    for old,new in zip(honest['all12_top_columns_lifted_R0_10'],top):
        assert old[5]==0 and new[0]==old[0]
        assert all(order[i]>=2 for i in old)
        assert order[old[0]]>=3
        order_records.append({'top_E0_node':old[0],'E0_order_lower_bound':order[old[0]],
                              'all_entries_order_lower_bound':min(order[i] for i in old),
                              'unit_pivot_component5_literal_zero':True,'Schur_E0_pointer_unchanged':True})

    # Complete algebraic gates: unit cone and actual cyclic nilcurve, not jets.
    ncone=read(ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json')
    assert ncone['N_full_minimal_free_resolution_ranks'][:3]==[8,36,84]
    assert unit['remaining_original_free_generator_indices']==[0,1,2,3,4,6,7,8,9]
    nbasis=[];exception=[]
    for alpha,col in enumerate(ncone['d1_universal_sparse_polynomial_columns']):
        a=[F(0)]*72
        n0terms=[]
        for row,poly in col:
            for c,mon in poly:
                if row==0:n0terms.append((int(c),tuple(mon)))
                if len(mon)==1:
                    assert mon[0]!=8  # nu occurs only in higher-order products.
                    j=mon[0] if mon[0]<8 else 8
                    a[9*row+j]+=int(c)
        assert a[:9]==[F(0)]*9
        if any(a):assert insert_basis(nbasis,a) is not None
        else:
            exception.append(alpha)
            assert [(c,m) for c,m in n0terms if len(m)==2]==[(1,(9,9))]
    assert len(nbasis)==35 and exception==[28]
    ev={'both57_unit_determinants':unit_determinants,'all114_saved_source2_series_all55_equations_independently_verified':True,
        'all2069_rebuilt_exact_firstjet90_columns':[[str(x) for x in a] for a in matrix],
        'independent43_pivot_original_columns':pivots,
        'independent43_span_basis':[[i,[[j,str(a)] for j,a in sorted(b.items())]] for i,b in basis],
        'all1576_coordinate_quotient_firstjets_independently_rebuilt':True,
        'four_unknown_g_derivative_sensitivities_identically_zero':True,
        'all12_top_full_firstjets_zero':True,'all12_top_order_records':order_records,
        'literal_N_linear_relation_rank_independently_rechecked':35,
        'literal_N_exceptional_relation_column':28,
        'structural_kernel_generators':{'total_upper_bound':44,'independent_linear_generators':43,
          'exceptional_relation':'(v^2-nu*z1^2)*n0-c0(z1)*F','exceptional_c0_order_at_least2_by_rank43':True,
          'all43_linear_generators_n0_component_order_at_least2':True,
          'exceptional_n0_quadratic_initial':'v^2'},
        'proof_conclusion':'alltop subset m*Jinf, hence Jinf=Q3 by Nakayama and T3=0',
        'D_depth19_or20_not_determined':True}
    ep=out/'independent-firstjet-and-structural-proof.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',fp,fr,sp,rp,hp,up,op,*proof_peers,
           ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63/actual-N-literal-cone-DAG.json',
           ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',
           ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'actual_T3_zero_by_complete_44_generator_kernel_and_independent_finite_gates',
        'independent_full9_relation_firstjet_rank':43,'all12_top_full_firstjets_zero':True,
        'all12_top_E0_in_m_cubed_after_exact_Schur':True,
        'actual_Jinf_equals_Q3_proved':True,'actual_T3_is_zero_proved':True,
        'requires_no_actual36_extension_residue_evaluation':True,
        'actual_T0_T1_T2_or_full_D_depth_solved':False,
        'no_new_source_solve_jets_above_saved2_GB_points_or_primes':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
        'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'certificate_bytes':ep.stat().st_size}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','independent_full9_relation_firstjet_rank','actual_T3_is_zero_proved','seconds','max_rss_kib','certificate_bytes']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
