#!/usr/bin/env python3
"""Independent exact Schur equivalence of the honest10 local presentation."""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def rank(rows):
    basis={}
    for row in rows:
        a=list(row)
        for j,b in sorted(basis.items()):
            if a[j]:c=a[j];a=[x-c*y for x,y in zip(a,b)]
        j=next((i for i,x in enumerate(a)if x),None)
        if j is not None:c=a[j];basis[j]=[x/c for x in a]
    return len(basis)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63'
    author=read(base/'report.json');cert=read(base/'unit-reduced-presentation-DAG.json');assert author['all_passed']
    for rel,h in author['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in author['evidence_sha256'].items():assert sha(base/rel)==h,rel
    honestpath=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    peer=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63_REVIEW/report.json'
    assert read(peer)['all_passed'];honest=read(honestpath)
    sourcepath=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    oldorigin=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    nodes=[tuple(n)for n in cert['DAG_scalar_nodes']]
    assert nodes[:len(honest['DAG_scalar_nodes'])]==[tuple(n)for n in honest['DAG_scalar_nodes']]
    vals=list(map(F,read(oldorigin)['all339669_scalar_origin_values_evaluated']))
    quot={q['node']:F(q['origin'])for rec in honest['all_pure_source_Hadamard_membership_telescopes']for q in rec['all9_regular_coordinate_divisions']}
    assert len(quot)==792
    for i,n in enumerate(nodes[len(vals):],len(vals)):
        op=n[0]
        if op=='c':v=F(n[1],n[2])
        elif op=='v':v=F(0)
        elif op=='r':v=quot[i]
        elif op=='n':v=-vals[n[1]]
        elif op=='a':v=vals[n[1]]+vals[n[2]]
        elif op=='m':v=vals[n[1]]*vals[n[2]]
        elif op=='i':assert vals[n[1]];v=1/vals[n[1]]
        else:raise ValueError(n)
        vals.append(v)
    assert [[i,str(v)]for i,v in sorted(quot.items())]==cert['all792_coordinate_quotient_origin_values']
    oldA=honest['honest_D3_presentation_relation_columns_R0_10'];oldT=honest['all12_top_columns_lifted_R0_10']
    A0=[[vals[i]for i in col]for col in oldA];T0=[[vals[i]for i in col]for col in oldT]
    assert A0==[list(map(F,col))for col in cert['honest10_origin_relation_matrix']]
    assert T0==[list(map(F,col))for col in cert['honest10_origin_top12_matrix']]
    assert rank(A0)==rank(A0+T0)==1 and rank(T0)==0
    step,=cert['all_constant_unit_pivot_steps'];pivot=oldA[1389]
    assert step['eliminated_relation_original_column']==1389 and step['eliminated_free_generator_original_index']==5
    assert step['relation_before_projection']==pivot and step['pivot_node']==pivot[5]
    unit=pivot[5];inverse=step['inverse_unit_node']
    assert vals[unit]==F(10,3) and nodes[inverse]==('i',unit)
    # Use only the independent DAG operator lookup rules; no author Schur code.
    D=DAGReview(read(sourcepath)['DAG_scalar_nodes']);D.nodes=nodes;D.values=vals;D.cache={n:i for i,n in enumerate(nodes)}
    keep=[0,1,2,3,4,6,7,8,9];assert keep==cert['remaining_original_free_generator_indices']
    def projection(v):return [D.add(v[k],D.neg(D.mul(pivot[k],D.mul(inverse,v[5]))))for k in keep]
    # Multiplication has different association from the primary. Verify the
    # stored coefficient as a scalar product multiset instead of pointer form.
    def factors(i):
        n=nodes[i]
        if n[0]=='m':return factors(n[1])+factors(n[2])
        return [i]
    def same_scalar(a,b):
        if a==b:return True
        return sorted(factors(a))==sorted(factors(b))
    def expected_entry(v,k,actual):
        if v[5]==D.zero:return actual==v[k]
        n=nodes[actual]
        if v[k]==D.zero:
            assert n[0]=='n';prod=n[1]
        else:
            assert n[0]=='a' and v[k]in n[1:]
            other=n[2]if n[1]==v[k]else n[1];assert nodes[other][0]=='n';prod=nodes[other][1]
        expected=[q for q in(pivot[k],inverse,v[5])if q!=D.one]
        if D.zero in expected:return actual==v[k]
        return sorted(factors(prod))==sorted(q for i in expected for q in factors(i))
    newA=cert['all_reduced_D3_relation_columns'];newT=cert['all12_transformed_top_columns']
    remaining=[i for i in range(2070)if i!=1389]
    assert remaining==cert['remaining_D3_relation_original_indices']
    for index,col in zip(remaining,newA):
        assert all(expected_entry(oldA[index],k,col[j])for j,k in enumerate(keep))
    for old,col in zip(oldT,newT):assert all(expected_entry(old,k,col[j])for j,k in enumerate(keep))
    P=cert['exact_projection_matrix_remaining_by_original10']
    for j,k in enumerate(keep):
        assert all(P[j][l]==(D.one if l==k else D.zero)for l in keep)
        assert nodes[P[j][5]][0]=='n'
        assert sorted(factors(nodes[P[j][5]][1]))==sorted(factors(pivot[k])+factors(inverse))
    assert all(vals[i]==0 for col in newA+newT for i in col)
    assert not cert['literally_zero_reduced_D3_column_positions']
    assert newA==cert['reduced_D3_relations_after_dropping_only_literal_zero_columns']
    assert newA+newT==cert['honest_reduced_D4_relations']
    assert len(newA)==2069 and len(newA+newT)==2081
    proof=dict(honest10_peer_bound=True,all_scalar_origins_evaluated_by_typed_semantics=True,
        original_relation_rank=1,original_relation_plus_top_rank=1,top_origin_rank=0,
        pivot_relation=1389,pivot_free_generator=5,unit_origin='10/3',
        exact_Schur='v_k -> v_k-c_k*u^-1*v_5 for k!=5; this is a local module isomorphism, not a jet operation',
        all2069_relations_and12_top_columns_verified_by_independent_product_multisets=True,
        complete_Koszul_kernel_retained_via_the_prior_exact_presentation=True,
        D3='coker(R0^2069->R0^9)',D4='coker(R0^2081->R0^9)',
        exact_minimal_module_generator_counts=[9,9],
        reason='Every remaining relation coefficient lies in the local maximal ideal, so tensoring the cokernel with the residue field has dimension9.',
        T3='image12 in D3 only; its syzygies/nonzero/depth remain uncomputed',
        no_global_unique_DAG_key_assumption=True,no_new_jets_GB_points_or_nonunit_inverse=True)
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'unit-reduced-presentation-DAG.json',
           ROOT/'code/ambient_B18_T3_honest10_unit_reduction_20261104_v63.py',
           ROOT/'research/notes/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63.md',honestpath,peer,sourcepath,oldorigin]
    report=dict(id=ID,all_passed=True,independent_scalar_origins_rank1_and_unit10_over3=True,
        exact_free9_Schur_all_relations_and_top12_accepted=True,mu_D3=9,mu_D4=9,
        D3_relation_columns=2069,D4_relation_columns=2081,
        T3_relation_kernel_nonzero_membership_depth_or_nil_not_claimed=True,
        no_author_rerun_new_jets_points_GB_or_nonunit_inverses=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={ep.name:sha(ep)},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in('all_passed','mu_D3','seconds','max_rss_kib')},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
