#!/usr/bin/env python3
"""Fixed exact-certificate audit of finiteTor63 targets; no second-carry solve."""
import argparse, hashlib, itertools, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_zero_normalizer_quadratic_review_20261104_v63 import decode,composition,vector,rank,determinant
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_FINITE_TOR63_IMAGE_20261104_v63_REVIEW'
MON=[()]+[(i,) for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2));POS={m:i for i,m in enumerate(MON)}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def qq(v):return list(map(F,v))
def sparse(v):return {i:F(c) for i,c in enumerate(v) if F(c)}
def linear_combination(cols,coeff,residual=()):
    out={int(i):F(c) for i,c in residual if F(c)}
    for k,c in coeff:
        k=int(k);c=F(c)
        for i,v in sparse(cols[k]).items():
            out[i]=out.get(i,F(0))+c*v
            if not out[i]:del out[i]
    return out
def multiples(first,selected):
    out=[]
    for k in selected:
        for variable in range(9):
            a=[F(0)]*405
            for i,c in sparse(first[k]).items():
                gen,j=divmod(i,9);a[45*gen+POS[tuple(sorted((j,variable)))]-10]+=c
            out.append(a)
    return out
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);assert d['all_passed'];bind(p)
        for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;bind(p.parent/n)
        return d
    base=ROOT/'reports';names=['AMBIENT_B18_FINITE_TOR63_IMAGE_PREFLIGHT_20261104_v63','AMBIENT_B18_TOR63_QUADRATIC_STOP_GATE_20261104_v63','AMBIENT_B18_TOR63_FULL_FIRST_MULTIPLE_STOP_20261104_v63']
    for n in names:audit(base/n/'report.json')
    for n in ['AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63_REVIEW','AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW','AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW','AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW']:
        audit(base/n/'report.json')
    T=read(base/names[0]/'finite63-proof-and-input-preflight.json')
    Q=read(base/names[1]/'Tor63-quadratic-targets-and-intrinsic-stop.json')
    U=read(base/names[2]/'full-first-multiples-and-Tor63-targets.json')
    E=read(base/'AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json')
    C=read(base/'AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json')
    J=read(base/'AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    joint=read(base/'AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    source=read(base/'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    rows=[18,19,28,29,60,67,72,63]
    assert T['actual_ideal_generator_physical_rows8']==rows
    for fam,key in [('CA','IA_generator_source_nodes8'),('CB','IB_generator_source_nodes8')]:
        assert T[key]==[p[1] for p in source['families'][fam]['actual_zero_odd_I0_generators']]
    X=[decode(p) for p in E['families']['CA']['physical_R0_readback_jets_in_own_F9']]
    pairs={f:{r:[decode(p) for p in pp] for r,pp in E['families'][f]['all12_trace_pairs_in_common_R0_through2']} for f in ('CA','CB')}
    ca=[composition(pairs['CA'][r][1],X) for r in rows]
    cb=[composition(pairs['CB'][r][1],X) for r in rows]
    # Direct common-physical-root compositions bind the constant16x10 matrix;
    # no unexpanded all-order cycle coefficients are treated as known numbers.
    A=[]
    for p in ca+cb:
        lin=[p.get((i,),F(0)) for i in range(7)]
        assert not p.get((7,),F(0)) and not p.get((8,),F(0))
        A.append(lin+[p.get(m,F(0)) for m in ((7,7),(7,8),(8,8))])
    assert A==[qq(row) for row in T['actual_constant16_by10_generator_matrix']]
    assert all(A[i]==[F(i==j) for j in range(7)]+[F(0)]*3 for i in range(7))
    restricted={p['CB_physical_zero_row']:qq(p['all_degree_le3_coefficients'][3:6]) for p in joint['all8_CB_zero_restrictions']}
    assert all(A[8+i][7:]==restricted[r] for i,r in enumerate(rows))
    selected=[0,1,2,3,4,5,6,8,12,13];redundant=[i for i in range(16) if i not in selected]
    assert selected==T['selected10_combined_indices'] and redundant==T['redundant6_combined_indices']
    V=[A[i] for i in selected];Vi=[qq(row) for row in T['selected10_unit_inverse']]
    assert determinant(V)==F(625,768)==F(T['selected10_unit_minor'])
    for M,N in ((V,Vi),(Vi,V)):
        assert all(sum(M[i][k]*N[k][j] for k in range(10))==int(i==j) for i in range(10) for j in range(10))
    q0=[]
    for k in redundant:
        c=[sum(A[k][i]*Vi[i][j] for i in range(10)) for j in range(10)]
        q=[F(i==k) for i in range(16)]
        for i,v in zip(selected,c):q[i]-=v
        assert all(sum(q[i]*A[i][j] for i in range(16))==0 for j in range(10));q0.append(q)
    assert q0==[qq(q) for q in T['all6_zero_column_cycle_origins16']]
    # C0 origin images, including f63's seven Hadamard coefficients and gA.
    Xi=[[F(0)]+[F(i==j) for j in range(7)]+[F(0)] for i in range(7)]
    Xi.append([F(0)]+A[7][:7]+[F(1)])
    assert Xi==[qq(r) for r in Q['actual_C0_basis_origin_images_of_CA8_odd_generators']]
    from ambient_B18_zero_normalizer_quadratic_review_20261104_v63 import plus,scale
    delta=[composition(plus(pairs['CA'][r][0],scale(pairs['CB'][r][0],-1)),X) for r in rows]
    assert all(all(len(m)==2 for m in p) for p in delta)
    assert [vector(p) for p in delta]==Q['all8_even_differences_in_actual_CA_F9_through2']
    targets=[]
    assert len(Q['all15_Tor2_quadratic_targets_and_reductions'])==15
    assert len(U['all15_saved_Tor2_targets_and_exact_reductions'])==15
    for target_index,(alpha,beta) in enumerate(itertools.combinations(range(6),2)):
        h=[[q0[alpha][i]*q0[beta][8+j]-q0[beta][i]*q0[alpha][8+j] for j in range(8)] for i in range(8)]
        assert h==[qq(row) for row in Q['all15_Tor2_quadratic_targets_and_reductions'][target_index]['constant_h8_by8']]
        v=[F(0)]*405
        for i,j,g in itertools.product(range(8),range(8),range(9)):
            c=h[i][j]*Xi[i][g]
            if c:
                for m,a in delta[j].items():v[45*g+POS[m]-10]+=c*a
        targets.append(v)
    assert len(targets)==15
    first=[qq(v) for v in C['all640_exact_D0_firstjet81_columns']]
    assert first==[qq(v) for v in U['all640_actual_D0_first_columns81']] and len(first)==640
    intrinsic=[]
    for i,j in itertools.combinations(range(8),2):
        c=[F(0)]*81
        if j<7:c[9*(1+i)+j]=-1
        if i<7:c[9*(1+j)+i]=1
        intrinsic.append(c)
    assert first[:28]==intrinsic==[qq(v) for v in Q['all28_intrinsic_C0_Koszul_first_columns81']]
    c252=multiples(intrinsic,range(28));assert c252==[qq(v) for v in Q['all252_intrinsic_maxideal_times_first_columns405']]
    assert rank(c252)==196
    outside=[]
    for pair,v,rec in zip(itertools.combinations(range(6),2),targets,Q['all15_Tor2_quadratic_targets_and_reductions']):
        assert rec['q_cycle_pair']==list(pair) and v==qq(rec['actual_C0_quadratic_target405'])
        assert linear_combination(c252,rec['intrinsic252_membership_coefficients'],rec['exact_residual'])==sparse(v)
        assert rec['in_intrinsic_m_times_Koszul_first_span']==(not rec['exact_residual'])
        if rec['exact_residual']:outside.append(pair)
    assert outside==[(1,2),(1,3),(1,4),(1,5),(2,5),(3,5),(4,5)]
    # Verify all640 expressions in the original first columns and an exact
    # rank38 lower bound, then reconstruct every one of the342 multiples.
    sel=U['independent38_original_relation_indices'];assert len(sel)==38 and rank([first[i] for i in sel])==38
    assert len(U['all640_first_span_membership_coefficients'])==640
    for v,coeff in zip(first,U['all640_first_span_membership_coefficients']):
        assert all(int(k) in sel for k,c in coeff)
        assert linear_combination(first,coeff)==sparse(v)
    c342=multiples(first,sel);assert len(c342)==342
    assert c342==[qq(v) for v in U['all342_maxideal_times_independent38_columns405']]
    labels=[dict(actual_D0_original_relation_index=i,F9_multiplier_index=j) for i in sel for j in range(9)]
    assert labels==U['all342_column_provenance'] and rank(c342)==261
    for v,rec in zip(targets,U['all15_saved_Tor2_targets_and_exact_reductions']):
        assert v==qq(rec['saved_exact_quadratic_target405']) and not rec['exact_residual']
        assert linear_combination(c342,rec['full_first342_membership_coefficients'])==sparse(v)
    # Tiny explicit A3 and omega modules: squared maximal ideal acts zero.
    # For omega basis(u,v,su), r*v=-su, s*u=su, otherszero.
    ra=[[0,0,0],[1,0,0],[0,0,0]];sa=[[0,0,0],[0,0,0],[1,0,0]]
    rw=[[0,0,0],[0,0,0],[0,-1,0]];sw=[[0,0,0],[0,0,0],[1,0,0]]
    socles=[];minimal=[]
    for r,s in ((ra,sa),(rw,sw)):
        assert all(sum(P[i][k]*Q[k][j] for k in range(3))==0 for P,Q in ((r,r),(r,s),(s,r),(s,s)) for i in range(3) for j in range(3))
        socles.append(3-rank(r+s));minimal.append(3-rank([list(c) for c in zip(*r)]+[list(c) for c in zip(*s)]))
    assert socles==[2,1] and minimal==[1,2]
    assert 15*3+6*3==63 and 15*socles[0]+6*socles[1]==36 and 15*minimal[0]+6*minimal[1]==27
    assert len(T['all27_minimal_Tor2_generator_specs'])==27 and len(T['all63_QQ_basis_specs'])==63
    assert T['all27_minimal_Tor2_generator_specs']==[
        *[dict(type='q_wedge_q',q_indices=[i,j]) for i,j in itertools.combinations(range(6),2)],
        *[dict(type='eta_wedge_q',eta=e,q_index=i) for e in ('eta_u','eta_v') for i in range(6)]]
    assert T['all63_QQ_basis_specs']==[
        *[dict(type='q_wedge_q',q_indices=[i,j],scalar=s) for i,j in itertools.combinations(range(6),2) for s in ('1','r','s')],
        *[dict(type='eta_wedge_q',eta=e,q_index=i,scalar=s) for i in range(6) for e,s in (('eta_u','1'),('eta_v','1'),('eta_u','s'))]]
    proof=dict(raw_kernel_is_Tor2_not_Tor1=True,raw_module='A3^15 + omega_A3^6',length=63,socle=36,minimal_generators=27,
        image_in_D0_is_finite_m_R0_squared_annihilated=True,full_comparison_kernel_not_identified_with_image=True,
        actual_constant16x10_source_matrix_rebuilt_from_common_physical_E4_pairs=True,selected_unit_minor='625/768',
        all6_cycle_origins_source_bound=True,allorder_cycle_transform_exists_by_regular_Hadamard_and_unit3_but_not_expanded=True,
        all15_actual_quadratic_targets_reconstructed=True,intrinsic252_rank196=True,intrinsic_residual_pairs=list(map(list,outside)),
        all640_first_memberships_checked=True,first_rank=38,all342_source_multiples_checked=True,multiple_rank=261,
        all15_exact_target_memberships_checked=True,remaining12_eta_images_have_module_order_atleast3=True,
        entire_27_generated_finite_image_is_contained_in_m_R0_cubed_D0=True,
        finite_image_zero_nonzero_full_kernel_socle_or_depth_not_decided=True,
        no_full640_second_carry_newjets_author_worker_or_GB=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/'code/ambient_B18_zero_normalizer_quadratic_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Finite rawTor63 image, source-bound15 quadratic targets and exact full-first multiple membership; all27/image in m_R0^3D0 only, no vanishing or depth result.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
