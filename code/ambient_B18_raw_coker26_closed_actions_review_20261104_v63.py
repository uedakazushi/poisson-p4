#!/usr/bin/env python3
"""Independent sparse-vector quotient/action review; no producer invocation."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def vec(x):return {int(i):F(c) for i,c in x if F(c)}
def plus(a,b,c=F(1)):
    q=dict(a)
    for j,v in b.items():
        x=q.get(j,F(0))+c*v
        if x:q[j]=x
        else:q.pop(j,None)
    return q
def apply(A,v):
    q={}
    for j,c in v.items():q=plus(q,A[j],c)
    return q
def three_product(a,b):
    q={}
    for i,x in a.items():
        for j,y in b.items():
            if i+j<3:q=plus(q,{i+j:x*y})
    return q
def image_polynomial(records,values):
    result={}
    for mon,c in records:
        q={0:F(c)}
        for i,e in enumerate(mon):
            for _ in range(e):q=three_product(q,values[i])
        result=plus(result,q)
    return result
def poly_apply(records,actions,v):
    result={}
    for mon,c in records:
        q={j:F(c)*x for j,x in v.items()}
        for i,e in enumerate(mon):
            for _ in range(e):q=apply(actions[i],q)
        result=plus(result,q)
    return result
def block(A,B,C):
    result=[]
    for offset,M in [(0,A),(3,B),(6,C)]:
        result.extend([{i+offset:c for i,c in v.items()} for v in M])
    assert len(result)==26
    return result
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    b=ROOT/'reports/AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63';p=b/'report.json';r=read(p)
    assert sha(p)=='2a1633574d301cbebd39e75285111a763e33b5afb5c2cddfdf84b740870acacc' and r['all_passed']
    for q,h in r['input_sha256'].items():assert sha(ROOT/q)==h,q
    for q,h in r['evidence_sha256'].items():assert sha(b/q)==h,q
    ep=b/'actual-raw26-actions-and-quotients.json';E=read(ep)
    basepath=ROOT/E['source_H20_file'];assert sha(basepath)==E['source_H20_sha256'];H=read(basepath)
    peer=ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json';assert read(peer)['all_passed']
    assert E['actual_CA_F9_coordinate_order9']==[f'zA_{i}' for i in range(7)]+['r_CA','s_CA']
    assert E['source_conductor_generators_rows_and_ids']==H['source_conductor_generators_rows_and_ids']
    assert E['ell']==H['ell_exact']=='5*(s_CA-r_CA)'
    basis=list(map(tuple,H['quotient_basis_monomials']));assert list(map(list,basis))==E['H20_basis_monomials']
    actions=[[vec(q) for q in A] for A in H['multiplication8_column_sparse']]
    assert E['H20_coordinate_actions8_column_sparse']==H['multiplication8_column_sparse']
    beta=[F(rec[0][1]) if rec else F(0) for rec in H['beta2_after_ell']]
    assert beta==list(map(F,E['beta2_coefficients7']))
    images={'CA':[],'CB':[]}
    for m in basis:
        for side in ['CA','CB']:
            degree=m[7]+2*sum(m[:7]);c=F(1)
            for j,e in enumerate(m[:7]):
                if e:c*=F(0) if side=='CA' else (-beta[j])**e
            images[side].append({degree:c} if degree<3 and c else {})
    for side in ['CA','CB']:
        saved=E['quotient_maps3_by20'][side]
        columns=[{i:F(saved[i][j]) for i in range(3) if F(saved[i][j])} for j in range(20)]
        assert columns==images[side]
    indices=E['chosen_lift_basis_indices'];assert indices==[0,8,16]
    lift=[[F(v) for v in row] for row in E['chosen_common_lift20_by3']]
    assert all(lift[i][j]==F(i==indices[j]) for i in range(20) for j in range(3))
    for side in images:assert [images[side][j] for j in indices]==[{0:F(1)},{1:F(1)},{2:F(1)}]
    coords={'CA':[{} for _ in range(7)]+[{1:F(1)}],
            'CB':[{2:-c} if c else {} for c in beta]+[{1:F(1)}]}
    for side in images:
        for i,j in itertools.product(range(8),range(20)):
            assert apply(images[side],actions[i][j])==three_product(coords[side][i],images[side][j])
    # Recompute each basis scalar by sparse coordinate actions, independently
    # of the producer's dense-matrix multiplication implementation.
    basis_actions=[]
    for m in basis:
        A=[]
        for j in range(20):
            v={j:F(1)}
            for i,e in enumerate(m):
                for _ in range(e):v=apply(actions[i],v)
            A.append(v)
        basis_actions.append(A)
    assert basis_actions==[[vec(q) for q in A] for A in E['H20_basis_multiplication_column_sparse']]
    assert all(A[0]=={j:F(1)} for j,A in enumerate(basis_actions))
    checks=0
    for side in images:
        for i,j in itertools.product(range(20),repeat=2):
            assert apply(images[side],basis_actions[i][j])==three_product(images[side][i],images[side][j]);checks+=1
    assert checks==800
    for q in H['IA_truncated_generators']:assert not image_polynomial(q,coords['CA'])
    for q in H['IB_truncated_generators']:assert not image_polynomial(q,coords['CB'])
    assert E['source_IA_truncated_generators']==H['IA_truncated_generators'] and E['source_IB_truncated_generators']==H['IB_truncated_generators']
    rawactions=[]
    for i in range(8):
        A=[three_product(coords['CA'][i],{j:F(1)}) for j in range(3)]
        B=[three_product(coords['CB'][i],{j:F(1)}) for j in range(3)]
        rawactions.append(block(A,B,actions[i]))
    assert rawactions==[[vec(q) for q in A] for A in E['raw26_coordinate_actions8_column_sparse']]
    for i,j,k in itertools.product(range(8),range(8),range(26)):
        assert apply(rawactions[i],rawactions[j][k])==apply(rawactions[j],rawactions[i][k])
    # Every20 basis scalar must act through its two canonical quotients,
    # never by silently relabelling a raw physical coefficient as zA/r.
    for i,raw in enumerate(E['raw26_scalar_action_by_H20_basis_column_sparse']):
        A=[three_product(images['CA'][i],{j:F(1)}) for j in range(3)]
        B=[three_product(images['CB'][i],{j:F(1)}) for j in range(3)]
        assert list(map(vec,raw))==block(A,B,basis_actions[i])
    for q in H['all64_actual_pair_products_mod_m5']:
        for j in range(20):assert not poly_apply(q,actions,{j:F(1)})
        for j in range(26):assert not poly_apply(q,rawactions,{j:F(1)})
    proof={'own_engine':'sparse Fraction vector actions and truncated univariate convolution; no author import/run',
      'actual_CA_F9_frame_and_source20_bound':True,'both_canonical_quotient_maps_and_common_r_verified':True,
      'surjective_common_lifts':[0,8,16],'all16_coordinate_intertwiners':True,
      'all800_basis_multiplicativity_identities':checks,'IA_and_IB_target_annihilations':True,
      'all64_products_annihilate_both_H20_and_raw26':True,
      'all8_raw26_coordinate_actions_and20_basis_scalar_callbacks_verified':True,
      'all26_actions_commute':True,'raw_module_is_3_plus3_plus20':True,
      'right_exact_specialization_not_a_kernel_basechange':True,
      'coefficients_must_first_be_actual_evaluated_H20_vectors':True,
      'no_nu_delta_b_physicalrho_or_comparison_torsion_evaluation':True}
    ev=out/'independent-quotient-and-sparse-action-checks.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',p,ep,basepath,peer,
      ROOT/'code/ambient_B18_raw_coker26_closed_actions_20261104_v63.py',ROOT/'research/notes/AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63.md']
    rep={'id':STEM,'all_passed':True,'canonical_actual_maps20_to3_and_closed_raw26_actions_accepted':True,
      'all800_products16_intertwiners64annihilators_independently_checked':True,
      'no_producer_source_jet_GB_point_prime_or_new_inverse_run':True,
      'no_rawphysical_coefficient_or_generic_DVR_K_torsion_D0_depth_claim':True,
      'input_sha256':{str(q.relative_to(ROOT)):sha(q) for q in paths},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['id','all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
