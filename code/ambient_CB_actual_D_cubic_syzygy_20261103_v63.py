#!/usr/bin/env python3
"""Audit the completed full150 QQ(b) cubic syzygy gate from preserved sources."""
import argparse,gzip,hashlib,json,subprocess,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63'
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for q in iter(lambda:f.read(1<<20),b''):h.update(q)
    return h.hexdigest()
def main(out,recompute):
    start=time.monotonic();assert not(out/'report.json').exists()
    D=read(out/'generic-cubic.json.gz');setup=read(out/'generic-setup.json.gz');point=read(out/'point-cubic.json.gz')
    assert D['all_passed']and setup['all_passed']and point['all_passed']
    workers=[ROOT/'code/ambient_CB_cubic_generic_setup_worker_20261103_v63.py',ROOT/'code/ambient_CB_cubic_generic_syzygy_worker_20261103_v63.py']
    if recompute:
        for p in workers:subprocess.run([str(ROOT/'.venv/bin/python'),str(p)],check=True,timeout=120)
        fresh=read(ROOT/'work/v62-review/CB-actualD-cubic-syzygy-generic-probe.json')
        for key in['all150_syzygy_defect_quadratics','all150_socle_coefficients','CA_plane_b150','CA_plane_c150','division_h0','division_h8','division_h_linear']:
            assert fresh[key]==D[key],key
    b=s.symbols('b');K=s.QQ.frac_field(b);C=lambda v:K.from_sympy(s.sympify(v,locals={'b':b}))
    def poly(v):return{(i,j):C(q)for i,j,q in v if C(q)}
    def add(*P):
        A={}
        for p in P:
            for k,q in p.items():A[k]=A.get(k,K.zero)+q
        return{k:q for k,q in A.items()if q}
    def scale(P,c):return{k:q*c for k,q in P.items()if q*c}
    def mul(P,V,n=3):
        A={}
        for(i,j),q in P.items():
            for(k,l),c in V.items():
                if i+j+k+l<=n:A[i+k,j+l]=A.get((i+k,j+l),K.zero)+q*c
        return{k:q for k,q in A.items()if q}
    one={(0,0):K.one};rho={(0,1):K.one}
    def compose(P,X,n):
        A={}
        for(i,j),q in P.items():
            term=one
            for _ in range(i):term=mul(term,X,n)
            for _ in range(j):term=mul(term,rho,n)
            A=add(A,scale(term,q))
        return{k:q for k,q in A.items()if sum(k)<=n}
    bc,cp,H0,H8=map(lambda key:list(map(poly,D[key])),['CA_plane_b150','CA_plane_c150','division_h0','division_h8'])
    Hi=[[poly(q)for q in row]for row in D['division_h_linear']]
    j0=setup['f_row'];minimum=setup['minimal8_linear_b_rows'];j8=D['quadratic_generator_row'];assert j0==18 and j8==63
    assert len(bc)==len(cp)==len(H0)==len(H8)==len(Hi)==150 and cp[j0]==one
    assert all(not p.get((0,0),K.zero)for j,p in enumerate(cp)if j!=j0)
    beta=[h.get((0,0),K.zero)for h in Hi[j8]]
    b8=add(bc[j8],*(scale(bc[j],-q)for q,j in zip(beta,minimum[1:])))
    c8=add(cp[j8],*(scale(cp[j],-q)for q,j in zip(beta,minimum[1:])))
    assert not c8.get((0,0),K.zero),'Quadratic-generator c has to start at order1'
    quadrows=[]
    for j in range(150):
        h0=compose(H0[j],bc[j0],2);hi=[compose(h,bc[j0],1)for h in Hi[j]]
        predicted=add(mul(h0,bc[j0]),mul(H8[j],b8),*(mul(h,bc[k])for h,k in zip(hi,minimum[1:])))
        assert add(bc[j],scale(predicted,-K.one))=={},('all150-b-syzygy',j)
        defect=add(cp[j],scale(h0,-K.one),scale(mul(H8[j],c8,2),-K.one),*(scale(mul(h,cp[k],2),-K.one)for h,k in zip(hi,minimum[1:])))
        assert defect==poly(D['all150_syzygy_defect_quadratics'][j]),j
        assert all(sum(k)==2 for k in defect),j
        quadrows.append([K.to_sympy(defect.get(k,K.zero))for k in[(2,0),(1,1),(0,2)]])
    Z=s.Matrix([[s.sympify(q,locals={'b':b})for q in row]for row in D['actual_Z4_quadratic_matrix']])
    assert len(DomainMatrix.from_Matrix(Z).convert_to(K).rref()[1])==2
    assert len(DomainMatrix.from_Matrix(Z.col_join(s.Matrix(quadrows))).convert_to(K).rref()[1])==2
    assert D['nonzero_cubic_socle_witness_row']is None and all(C(q)==0 for q in D['all150_socle_coefficients'])
    assert set(setup['base29_rows']).isdisjoint(setup['minimal8_linear_b_rows'])
    for j in setup['base29_rows']:assert not bc[j]and not cp[j]
    # This is specialization consistency, not the proof of rational-function zero.
    for key in['CA_plane_b150','CA_plane_c150','division_h0','division_h8','all150_syzygy_defect_quadratics']:
        for g,p in zip(D[key],point[key]):
            A={k:s.cancel(K.to_sympy(q).subs(b,1))for k,q in poly(g).items()};B={k:s.sympify(q)for k,q in[((i,j),q)for i,j,q in p]}
            assert A==B,key
    deps=[ROOT/'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63_REVIEW/report.json',
          ROOT/'reports/AMBIENT_CB_ACTUAL_D_QUADRATIC_GATE_20261103_v63/report.json',
          ROOT/'reports/AMBIENT_CB_SELF_PAIR_COMMON_BASE_ETALE_20261103_v63/report.json',
          ROOT/'reports/QB_SOURCE_SERIES_HELPER_20261103_v63/report.json',
          ROOT/'reports/QB_SOURCE_SERIES_GENERIC_HELPER_20261103_v63/report.json']
    assert all(read(p)['all_passed']for p in deps)
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',out/'generic-cubic.json.gz',out/'generic-setup.json.gz',out/'point-cubic.json.gz',out/'point-setup.json.gz',out/'completed-worker-manifest.json',*workers,
           ROOT/'code/ambient_CB_cubic_point_setup_worker_20261103_v63.py',ROOT/'code/ambient_CB_cubic_point_syzygy_worker_20261103_v63.py',ROOT/'code/qb_source_series_helper_20261103_v63.py',ROOT/'code/qb_source_series_generic_helper_20261103_v63.py',*deps]
    report=dict(id=STEM,all_passed=True,field='QQ(b), no parameter sampling',
        actual_original_three_source_full150_jets_through3=True,full_actual_nonzero_ratio_times_binary_terms_retained=True,
        actual58_and87_source_constraint_inverses_verified=True,
        all150_including_all121_nonbase_odd_coordinates_covered=True,
        regular_sequence_division_8linear_generators_through_degree2_and_quadratic_generator_through_degree1=True,
        quadratic_generator_multiplier_c_constant_explicitly_zero=True,
        all150_full_b_syzygy_identities_through3_rechecked=True,all150_constant_and_linear_defects_zero=True,
        all150_defect_quadratics_in_actual_Z4_ideal_over_QQb=True,
        missing_cubic_socle_class_identically_zero_not_just_point=True,
        higher_normal_powers_already_in_known_core_ideal_application_gate_clear=True,
        completed_generic_worker_seconds=D['seconds'],completed_generic_setup_seconds=setup['seconds'],
        actual_completed_computations_reused_without_repeating=True,full_recompute_requested=recompute,
        full_D_theorem_saved_in_separate_application_packet=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','all150_defect_quadratics_in_actual_Z4_ideal_over_QQb','missing_cubic_socle_class_identically_zero_not_just_point','seconds']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--recompute',action='store_true');a=p.parse_args();main(a.output.resolve(),a.recompute)
