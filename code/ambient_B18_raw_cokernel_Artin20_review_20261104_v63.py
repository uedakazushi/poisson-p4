#!/usr/bin/env python3
"""Independent substitution of saved Artin105 roots and fixed image certificate.

No Hensel corrections or source solver are executed by this checker.
"""
from collections import Counter
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

from ambient_B18_raw_coker26_closed_actions_20261104_v63 import (
    zeros,eye,mm,mv,add,scale,dense_columns,blockdiag)

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dense(record,n=20):
    value=[F(0)]*n
    for i,c in record:value[i]=F(c)
    return value
def vecadd(a,b):return tuple(x+y for x,y in zip(a,b))
def vecscale(a,c):return tuple(x*c for x in a)
def col(A,j):return [row[j] for row in A]


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):
        hashes[str(p.relative_to(ROOT))]=sha(p)
        return json.loads(p.read_text())
    def gate(stem):
        p=ROOT/'reports'/stem/'report.json';r=read(p);assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            q=p.parent/name;assert sha(q)==h,name;hashes[str(q.relative_to(ROOT))]=h
        for rec in r.get('input_hashes',{}).values():
            assert sha(ROOT/rec['path'])==rec['sha256'];hashes[rec['path']]=rec['sha256']
        return r
    b=ROOT/'reports'
    for stem in ['AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW',
        'AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63_REVIEW',
        'AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW',
        'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW']:
        gate(stem)
    primary=gate('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63')
    phase=read(ROOT/primary['phase1_evidence']['path']);assert hashes[primary['phase1_evidence']['path']]==primary['phase1_evidence']['sha256']
    acts=read(ROOT/primary['action_evidence']['path']);assert hashes[primary['action_evidence']['path']]==primary['action_evidence']['sha256']
    raw=read(ROOT/primary['raw26_evidence']['path']);assert hashes[primary['raw26_evidence']['path']]==primary['raw26_evidence']['sha256']
    A=read(b/'AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json')
    S=read(b/'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    J=read(b/'AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    P=read(b/'AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json')
    E=read(b/'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json')
    assert phase['A20_basis']==A['quotient_basis_monomials']==raw['H20_basis_monomials']
    nodes=S['DAG_scalar_nodes'];ids=J['all105_unknown_original_node_order'];eq=J['all105_equation_original_node_order']
    assert ids==phase['joint105_unknown_order']==P['unknown105_order']
    assert eq==phase['joint105_equation_order']==P['equation105_order']
    assert len(ids)==len(set(ids))==len(eq)==105
    requested=eq[:]
    for rec in P['physical12_coefficient_roots']:
        requested.extend(rec[k] for k in ['aA_node','aB_node','bA_node','bB_node'])
    nus=[S['families'][f]['actual_tracefree_nu'] for f in ['CA','CB']]
    requested+=nus
    for f in ['CA','CB']:requested.extend(i for _,i in S['families'][f]['actual_zero_odd_I0_generators'])
    requested=sorted(set(requested));assert requested==P['requested_root_nodes'] and len(requested)==164
    needed=set();todo=requested[:]
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i)
        if nodes[i][0] not in ['c','v']:todo.extend(nodes[i][1:])
    order=sorted(needed);assert order==P['needed_closure_order'] and len(order)==18584
    assert {i for i in order if nodes[i][0]=='v'}==set(ids)
    # Bind the whole original/pruned interface structurally, including source
    # root names. This adds no coefficient or source-function query.
    pruning=dict(E['original_node_to_pruned_node']);pruning[0]=0;pnodes=E['DAG_scalar_nodes']
    shared_prefix={old:new for old,new in pruning.items() if old<len(nodes)}
    for old,new in shared_prefix.items():
        op,*q=nodes[old];p=pnodes[new];assert op==p[0]
        if op=='c':assert F(q[0],q[1])==F(p[1],p[2])
        elif op=='v':assert p[1]==f'ACTUAL_NODE{old}_{q[0]}'
        elif op in ['a','m']:assert sorted(pruning[k] for k in q)==sorted(p[1:])
        else:assert [pruning[k] for k in q]==p[1:]
    for fam in ['CA','CB']:
        sys=next(s for s in E['normalizer_sources48_each'] if s['family']==fam)
        source_eq=S['families'][fam]['all48_zero_IFT_quadratic_coefficient_equations']
        assert [pruning[i] for i in source_eq]==sys['equations48']
    assert S['physical_zero_base_rows9']==E['actual_R0_rows9']==phase['actual_physical_R0_rows']
    assert [pruning[i] for i in range(2,11)]==E['R0_variable_nodes']

    # Multiplication is independently read from the accepted H20 coordinate
    # matrices. The producer used495 monomial reductions instead.
    basis=[tuple(m) for m in A['quotient_basis_monomials']]
    coordinate=[dense_columns(v,20) for v in A['multiplication8_column_sparse']]
    unit=tuple([F(1)]+[F(0)]*19);zero=(F(0),)*20
    table=[]
    for i,mon in enumerate(basis):
        B=eye(20)
        for k,e in enumerate(mon):
            for _ in range(e):B=mm(B,coordinate[k])
        assert mv(B,unit)==[F(i==j) for j in range(20)]
        for j in range(20):
            for k,c in enumerate(col(B,j)):
                if c:table.append((i,j,k,c))
    assert len(table)==80
    def multiply(a,b):
        if not any(a) or not any(b):return zero
        if not any(a[1:]):return vecscale(b,a[0])
        if not any(b[1:]):return vecscale(a,b[0])
        c=[F(0)]*20
        for i,j,k,q in table:
            if a[i] and b[j]:c[k]+=q*a[i]*b[j]
        return tuple(c)
    def invert(a):
        assert a[0]
        n=vecscale(a,-1/a[0]);n=list(n);n[0]+=1;n=tuple(n)
        answer=unit;term=unit
        for _ in range(4):term=multiply(term,n);answer=vecadd(answer,term)
        answer=vecscale(answer,1/a[0]);assert multiply(a,answer)==unit
        return answer
    unknown={i:tuple(dense(r)) for i,r in zip(ids,phase['all105_unknown_final_values'])}
    assert all(v[0]==0 for v in unknown.values())
    targets=[tuple(F(i==j) for i in range(20)) for j in range(1,8)]+[tuple(F(i==8) for i in range(20))]*2
    assert [list(v) for v in targets]==[dense(v) for v in phase['actual_CA_F9_targets']]
    uses=Counter()
    for i in order:
        if nodes[i][0] not in ['c','v']:uses.update(nodes[i][1:])
    uses.update(requested);values={};inverse_cert={q[0]:q for q in phase['unit_inverse_exact_checks']}
    inverse_count=0;peak=0
    (out/'preflight.json').write_text(json.dumps(dict(saved_roots_only=True,source_nodes18584=True,
        unknown_vectors105=True,A20_products=80,source_constants='Fraction(num,den)',
        no_J0_update_Hensel_sourcejet_or_author_worker=True,
        caps=dict(seconds=60,address_bytes=2_000_000_000,output_bytes=32*1024*1024,live_coefficients=150000)),indent=2)+'\n')
    for i in order:
        op,*q=nodes[i]
        if op=='c':c=F(q[0],q[1]);v=tuple([c]+[F(0)]*19)
        elif op=='v':v=unknown[i]
        elif op=='a':v=vecadd(values[q[0]],values[q[1]])
        elif op=='n':v=vecscale(values[q[0]],-1)
        elif op=='m':v=multiply(values[q[0]],values[q[1]])
        elif op=='i':
            v=invert(values[q[0]]);cert=inverse_cert[i]
            assert cert[1]==q[0] and tuple(dense(cert[2]))==values[q[0]] and tuple(dense(cert[3]))==v
            inverse_count+=1
        else:raise AssertionError(op)
        values[i]=v;peak=max(peak,len(values));assert 20*(peak+105)<150000
        if op not in ['c','v']:
            for j in q:
                uses[j]-=1
                if not uses[j]:values.pop(j)
    assert inverse_count==6 and set(values)==set(requested)
    assert {i:tuple(dense(v)) for i,v in phase['all_requested_original_source_root_values']}==values
    for k,i in enumerate(eq):assert values[i]==(targets[k-48] if 48<=k<57 else zero)
    assert all(not r for r in phase['all105_final_equation_residuals'])
    assert [list(unknown[i]) for i in range(2,11)]==[dense(v) for v in phase['actual_physical_R0_readback']]
    assert nus==phase['nu_roots_CA_CB'] and [list(values[i]) for i in nus]==[dense(v) for v in phase['nu_CA_CB_values']]
    # A positive exact root check and the accepted unit105 Jacobian prove these
    # are actual pointed Hensel values; correction history is not re-solved.
    assert len(phase['stages'])==4 and [q['stage'] for q in phase['stages']]==[1,2,3,4]

    pruned_cache={pruning[i]:v for i,v in unknown.items()}
    pruned_cache.update({pruning[i]:v for i,v in values.items() if i in pruning})
    def pvalue(i):
        if i in pruned_cache:return pruned_cache[i]
        op,*a=pnodes[i]
        if op=='c':v=tuple([F(a[0],a[1])]+[F(0)]*19)
        elif op=='a':v=vecadd(pvalue(a[0]),pvalue(a[1]))
        elif op=='n':v=vecscale(pvalue(a[0]),-1)
        elif op=='m':v=multiply(pvalue(a[0]),pvalue(a[1]))
        elif op=='i':v=invert(pvalue(a[0]))
        else:raise AssertionError(('unbound_actual_E4_coefficient',i,op))
        pruned_cache[i]=v;return v
    coefficient_values=[]
    for recipe,cert,e4 in zip(P['physical12_coefficient_roots'],phase['all12_physical_coefficients'],E['all12_zero_physical_difference_columns']):
        assert recipe['physical_row']==cert['physical_row']==e4['physical_row']
        v={k:values[recipe[k+'_node']] for k in ['aA','aB','bA','bB']}
        v['delta']=vecadd(v['aA'],vecscale(v['aB'],-1))
        for k,c in v.items():assert list(c)==dense(cert[k])
        assert [pvalue(i) for i in e4['column4']]==[v['delta'],v['bA'],vecscale(v['bB'],-1),zero]
        coefficient_values.append(v)
    # These same actual F coefficients, not just their quadratic approximation,
    # underlie the raw26 physical action formula.
    QA=[[F(c) for c in row] for row in raw['quotient_maps3_by20']['CA']]
    QB=[[F(c) for c in row] for row in raw['quotient_maps3_by20']['CB']]
    lift=[[F(c) for c in row] for row in raw['chosen_common_lift20_by3']]
    assert mm(QA,lift)==mm(QB,lift)==eye(3)
    assert raw['raw26_basis']==acts['raw26_basis'] and acts['sector_ranges']==[[0,3],[3,6],[6,26]]
    IA=[values[i] for _,i in raw['source_conductor_generators_rows_and_ids']['CA']]
    IB=[values[i] for _,i in raw['source_conductor_generators_rows_and_ids']['CB']]
    for v in IA:assert not any(mv(QA,v))
    for v in IB:assert not any(mv(QB,v))
    assert all(multiply(a,b)==zero for a in IA for b in IB)
    def mu(v):return [list(row) for row in zip(*(multiply(v,tuple(F(k==j) for k in range(20))) for j in range(20)))]
    physical=[]
    for i,v in enumerate(coefficient_values):
        delta,bA,bB=[v[k] for k in ['delta','bA','bB']]
        assert not any(mv(QA,bA)) and not any(mv(QB,bB))
        BAmu,BBmu,Dmu=mu(bA),mu(bB),mu(delta)
        # Prove independence of ALL possible quotient lifts, not just a chosen
        # representative: I-lift*pi spans the entire kernel of each quotient.
        assert mm(BBmu,add(eye(20),scale(mm(lift,QA),-1)))==zeros(20,20)
        assert mm(BAmu,add(eye(20),scale(mm(lift,QB),-1)))==zeros(20,20)
        matrix=zeros(26,26)
        blocks=[(0,0,mm(mm(QA,Dmu),lift)),(3,3,mm(mm(QB,Dmu),lift)),
            (6,0,scale(mm(BBmu,lift),-1)),(6,3,mm(BAmu,lift)),
            (0,6,mm(QA,mu(vecscale(multiply(values[nus[1]],bB),-1)))),
            (3,6,mm(QB,mu(multiply(values[nus[0]],bA)))),(6,6,Dmu)]
        for r,c,block in blocks:
            for k,row in enumerate(block):
                for j,a in enumerate(row):matrix[r+k][c+j]=a
        assert matrix==dense_columns(acts['physical_actions12_column_sparse'][i],26)
        physical.append(matrix)
    X=[dense_columns(q,26) for q in raw['raw26_coordinate_actions8_column_sparse']]
    for A0 in physical:
        for B0 in X:assert mm(A0,B0)==mm(B0,A0)
    for i,A0 in enumerate(physical):
        for B0 in physical[i+1:]:assert mm(A0,B0)==mm(B0,A0)
    columns=[col(M,j) for M in physical for j in range(26)]
    assert columns==[dense(r,26) for r in acts['all312_relation_columns_sparse']]
    # Replay the saved18 pivots; select no new pivot and perform no rank search.
    piv={}
    for receipt in acts['fixed_pivot_receipts']:
        p=columns[receipt['source_column']][:]
        for i,c in receipt['subtract']:
            assert p[i]==F(c)
            p=[a-F(c)*b for a,b in zip(p,piv[i])]
        i=receipt['pivot'];c=F(receipt['divide'])
        assert p[i]==c and c and all(not p[j] for j in range(i)) and i not in piv
        piv[i]=[a/c for a in p]
    assert len(piv)==acts['relation_rank']==18
    assert piv=={i:dense(v,26) for i,v in acts['echelon_rows']}
    for p in columns:
        p=p[:]
        for i in sorted(piv):p=[a-p[i]*b for a,b in zip(p,piv[i])]
        assert not any(p)
    free=acts['C_basis_raw26_indices'];assert free==[i for i in range(26) if i not in piv]
    proj=dense_columns(acts['all26_projection_columns'],8)
    assert len(free)==acts['C_special_dimension']==primary['actual_C_special_dimension']==8
    assert [[row[i] for i in free] for row in proj]==eye(8)
    for p in columns:assert not any(mv(proj,p))
    for i,M0 in enumerate(X):
        C0=dense_columns(acts['C_coordinate_actions8_column_sparse'][i],8)
        assert mm(proj,M0)==mm(C0,proj)
    proof=dict(all_passed=True,source_root_substitution='Own dense Fraction20 vectors, multiplication from accepted8 coordinate matrices; no producer evaluator or Hensel iteration.',
        shared_original_source_prefix_nodes_structurally_bound=len(shared_prefix),
        all18584_source_nodes_checked=True,all105_equations_exact_with_actual_CA_F9_targets=True,
        all164_saved_source_roots_equal=True,all6_unit_inverses_and_constant_denominators_checked=True,
        actual_physical_R0_readback9_checked=True,all12_aA_aB_delta_bA_bB_and_nu2_checked=True,
        all12_E4_F_columns_bound_to_same_actual_roots=True,
        all64_actual_source_conductor_products_zero=True,all_quotient_lift_kernel_matrices_zero=True,
        all12_actual26_actions_reconstructed=True,coordinate_intertwiners=96,physical_commutators=66,
        all312_columns_checked_and_fixed18_pivot_receipts_replayed=True,
        relation_rank=18,actual_C_mod_ell_dimension=8,all8_quotient_coordinate_actions_checked=True,
        actual_cokernel_proof='C_raw is quotiented by the image of all twelve actual physical differences. Their R0-linear actions are reconstructed on the whole26 basis; hence the sum of their312 image columns is the entire closed relation image. Right exactness computes C/(ell); no kernel is specialized.',
        no_source_solve_corrections_new_jets_GB_point_prime_generic_DVR_rank_K_torsion_or_depth=True,
        live_vector_peak=peak,live_source_coefficient_bound=20*(peak+105))
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'code/ambient_B18_raw_coker26_closed_actions_20261104_v63.py',ROOT/f'research/notes/{ID}.md']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,all105_actual_source_equations_independently_verified=True,
        all12_actual_raw26_actions_independently_verified=True,actual_C_mod_ell_dimension=8,
        no_generic_DVR_rank_K_torsion_socle_or_D0_depth_claim=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
