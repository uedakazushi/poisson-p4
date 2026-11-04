#!/usr/bin/env python3
"""Closed actual raw-cokernel26 actions from the accepted conductor-product20.

No source function, source jet, comparison map or torsion is evaluated here.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zeros(n,m):return [[F(0) for _ in range(m)] for _ in range(n)]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def mm(A,B):
    C=zeros(len(A),len(B[0]))
    for i,row in enumerate(A):
        for k,a in enumerate(row):
            if a:
                for j,b in enumerate(B[k]):
                    if b:C[i][j]+=a*b
    return C
def mv(A,v):return [sum((a*b for a,b in zip(row,v)),F(0)) for row in A]
def add(A,B):return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def scale(A,c):return [[a*c for a in row] for row in A]
def dense_columns(columns,n):
    A=zeros(n,len(columns))
    for j,column in enumerate(columns):
        for i,c in column:A[i][j]=F(c)
    return A
def sparse_columns(A):
    return [[[i,str(A[i][j])] for i in range(len(A)) if A[i][j]] for j in range(len(A[0]))]
def stringify(A):return [[str(c) for c in row] for row in A]
def blockdiag(A,B,C):
    out=zeros(len(A)+len(B)+len(C),len(A[0])+len(B[0])+len(C[0]))
    r=c=0
    for M in [A,B,C]:
        for i,row in enumerate(M):
            for j,a in enumerate(row):out[r+i][c+j]=a
        r+=len(M);c+=len(M[0])
    return out
def mul3(v):
    return [[v[i-j] if i>=j else F(0) for j in range(3)] for i in range(3)]


def multiplication20(coefficients, evidence):
    """Return multiplication by an ALREADY EVALUATED actual H20 coefficient.

    This callback does not turn a physical-coordinate label into a value.
    """
    assert len(coefficients)==20
    out=zeros(20,20)
    for c,columns in zip(coefficients,evidence['H20_basis_multiplication_column_sparse']):
        if F(c):out=add(out,scale(dense_columns(columns,20),F(c)))
    return out


def scalar_action26(coefficients,evidence):
    """Actual scalar action, given its previously evaluated H20 vector."""
    assert len(coefficients)==20
    A=[[F(c) for c in row] for row in evidence['quotient_maps3_by20']['CA']]
    B=[[F(c) for c in row] for row in evidence['quotient_maps3_by20']['CB']]
    v=[F(c) for c in coefficients]
    return blockdiag(mul3(mv(A,v)),mul3(mv(B,v)),multiplication20(v,evidence))


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
        for record in r.get('input_hashes',{}).values():
            q=ROOT/record['path'];assert sha(q)==record['sha256'];hashes[record['path']]=record['sha256']
        return r
    gate('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW')
    gate('AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW')
    primary=gate('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63')
    path=ROOT/'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json'
    x=read(path);assert sha(path)==primary['evidence']['sha256']
    assert x['variables']==[f'zA_{i}' for i in range(7)]+['r_CA']
    basis=[tuple(m) for m in x['quotient_basis_monomials']]
    assert len(basis)==20 and basis[0]==(0,)*8
    idx={m:i for i,m in enumerate(basis)}
    coord=[]
    for i in range(8):m=[0]*8;m[i]=1;coord.append(idx[tuple(m)])
    assert coord==list(range(1,9))
    raw=[dense_columns(c,20) for c in x['multiplication8_column_sparse']]
    assert len(raw)==8 and all(len(a)==20 and len(a[0])==20 for a in raw)
    (out/'preflight.json').write_text(json.dumps(dict(source_H20_bytes=path.stat().st_size,
        source_actions8=20,basis_actions20=20,raw_coordinate_actions8=26,
        quotient_maps2=[3,20],source_ideal_generator_pairs=64,
        no_new_source_inverse_jet_point_prime_GB_or_rank_search=True,
        resource_caps=dict(seconds=60,address_bytes=2_000_000_000,output_bytes=32*1024*1024,polynomial_terms=150000)),indent=2)+'\n')
    for A in raw:
        for B in raw:assert mm(A,B)==mm(B,A)

    # The two canonical H20 quotient maps use actual CA-F9 r on BOTH sides.
    # beta3 and all later terms vanish modulo r^3 after z_i=-c_i*r^2.
    beta=[]
    for rec in x['beta2_after_ell']:
        assert not rec or len(rec)==1 and rec[0][0]==[0]*7+[2]
        beta.append(F(rec[0][1]) if rec else F(0))
    assert beta==list(map(F,['9835/864',0,0,0,'29875/3456','-4565/864',0]))
    maps=[]
    for side in ['CA','CB']:
        Q=zeros(3,20)
        for j,m in enumerate(basis):
            degree=m[7]+2*sum(m[:7]);c=F(1)
            if side=='CA' and any(m[:7]):c=F(0)
            elif side=='CB':
                for k,e in enumerate(m[:7]):c*=(-beta[k])**e
            if degree<3:Q[degree][j]=c
        maps.append(Q)
    PA,PB=maps
    lifts=zeros(20,3)
    lift_indices=[idx[(0,)*7+(j,)] for j in range(3)]
    assert lift_indices==[0,8,16]
    for j,i in enumerate(lift_indices):lifts[i][j]=1
    assert mm(PA,lifts)==mm(PB,lifts)==eye(3)
    JR=zeros(3,3);JR[1][0]=1;JR[2][1]=1
    JA=[zeros(3,3) for _ in range(7)]+[JR]
    JB=[scale(mm(JR,JR),-c) for c in beta]+[JR]
    for i in range(8):
        assert mm(PA,raw[i])==mm(JA[i],PA)
        assert mm(PB,raw[i])==mm(JB[i],PB)

    def poly_value(records,actions):
        ans=zeros(len(actions[0]),len(actions[0]))
        for mon,c in records:
            term=eye(len(actions[0]))
            for k,e in enumerate(mon):
                for _ in range(e):term=mm(term,actions[k])
            ans=add(ans,scale(term,F(c)))
        return ans
    IA=x['IA_truncated_generators'];IB=x['IB_truncated_generators']
    assert len(IA)==len(IB)==8
    for gen in IA:assert poly_value(gen,JA)==zeros(3,3)
    for gen in IB:assert poly_value(gen,JB)==zeros(3,3)
    assert poly_value(x['gA_after_ell'],JA)==zeros(3,3)
    assert poly_value(x['gB_after_ell'],JB)==zeros(3,3)
    assert x['gA_after_ell']==[[[0]*7+[3],'125/24']]
    assert x['gB_after_ell']==[[[0]*7+[3],'-15125/144']]
    # Re-evaluate every literal IA*IB polynomial on the accepted exact H20.
    for gen in x['all64_actual_pair_products_mod_m5']:
        assert poly_value(gen,raw)==zeros(20,20)

    basis20=[]
    for m in basis:
        A=eye(20)
        for k,e in enumerate(m):
            for _ in range(e):A=mm(A,raw[k])
        basis20.append(A)
    e0=[F(i==0) for i in range(20)]
    for j,A in enumerate(basis20):assert mv(A,e0)==[F(i==j) for i in range(20)]
    for j,A in enumerate(basis20):
        assert mm(PA,A)==mm(mul3([row[j] for row in PA]),PA)
        assert mm(PB,A)==mm(mul3([row[j] for row in PB]),PB)
    # These are canonical algebra quotients, not merely arbitrary vector maps.
    for i,A in enumerate(basis20):
        for j in range(20):
            product=[row[j] for row in A]
            assert mv(PA,product)==mv(mul3([row[i] for row in PA]),[row[j] for row in PA])
            assert mv(PB,product)==mv(mul3([row[i] for row in PB]),[row[j] for row in PB])

    action26=[blockdiag(A,B,H) for A,B,H in zip(JA,JB,raw)]
    source_basis26=[]
    for side in ['tauA','tauB']:
        source_basis26.extend([dict(sector=side,CA_F9_r_power=j) for j in range(3)])
    source_basis26.extend([dict(sector='tauA*tauB',H20_basis_index=j,CA_F9_monomial=list(m)) for j,m in enumerate(basis)])
    for A in action26:
        for B in action26:assert mm(A,B)==mm(B,A)
    for a in IA:
        aa=poly_value(a,action26)
        for b in IB:assert mm(aa,poly_value(b,action26))==zeros(26,26)
    evidence=dict(source_H20_file=str(path.relative_to(ROOT)),source_H20_sha256=sha(path),
        exact_common_base='Actual CA-F9 coordinates, not raw physical rho',
        actual_CA_F9_coordinate_order9=[f'zA_{i}' for i in range(7)]+['r_CA','s_CA'],
        source_conductor_generators_rows_and_ids=x['source_conductor_generators_rows_and_ids'],
        ell=x['ell_exact'],ell_zero_is_exact_s_equals_r=True,
        CA_quotient='R0/(IA,ell) = QQ[r_CA]/r_CA^3',CB_quotient='R0/(IB,ell) = QQ[r_CA]/r_CA^3',
        common3_basis=['1','r_CA','r_CA^2'],H20_basis_monomials=[list(m) for m in basis],
        beta2_coefficients7=list(map(str,beta)),
        quotient_maps3_by20=dict(CA=stringify(PA),CB=stringify(PB)),
        quotient_map_definition=dict(CA='zA_i=0,r_CA=r',CB='zA_i=-beta2_i*r^2,r_CA=r',both='r^3=0'),
        chosen_common_lift20_by3=stringify(lifts),chosen_lift_basis_indices=lift_indices,
        H20_coordinate_actions8_column_sparse=x['multiplication8_column_sparse'],
        H20_basis_multiplication_column_sparse=[sparse_columns(A) for A in basis20],
        raw26_basis=source_basis26,raw26_sector_ranges=dict(tauA=[0,3],tauB=[3,6],tauA_tauB=[6,26]),
        raw26_coordinate_actions8_column_sparse=[sparse_columns(A) for A in action26],
        ninth_s_CA_action_equals_coordinate7_r_CA=True,
        raw26_scalar_action_by_H20_basis_column_sparse=[sparse_columns(blockdiag(mul3([row[j] for row in PA]),mul3([row[j] for row in PB]),basis20[j])) for j in range(20)],
        source_IA_truncated_generators=IA,source_IB_truncated_generators=IB,
        actual_raw_module_decomposition='(R0/(IA,ell))*tauA directsum (R0/(IB,ell))*tauB directsum (R0/(IA*IB,ell))*tauA*tauB',
        quotient_exactness_proof='The source-bound ideal generators are zA7,gA and the legal unit-normalized zA+beta2+beta3+O(m4),gB. On either QQ[r]/r3 target all omitted terms vanish. gA/gB have nonzero cubic leading units, so these are surjections from the actual length3 conductor quotients and therefore isomorphisms. The H20 quotient maps are induced by IA*IB subset IA and IB. Intertwining with every actual H20 generator verifies the same canonical maps on its fixed basis.',
        raw_module_no_new_ring_structure_claim=True,
        no_physical_rho_nu_delta_a_b_function_values_or_comparison_torsion_evaluated=True)
    # Exercise the public callbacks on every accepted generator and every basis
    # scalar, still without evaluating any previously unknown source function.
    for i in range(8):assert scalar_action26([F(j==coord[i]) for j in range(20)],evidence)==action26[i]
    for j in range(20):assert multiplication20([F(i==j) for i in range(20)],evidence)==basis20[j]
    (out/'actual-raw26-actions-and-quotients.json').write_text(json.dumps(evidence,indent=2)+'\n')
    proof=dict(all_passed=True,canonical_actual_quotient_maps20_to3=True,
        both_maps_surjective_via_lift_indices=[0,8,16],all16_generator_intertwining_identities=True,
        all800_basis_quotient_multiplicativity_checks=True,all8_IA_and8_IB_target_annihilators=True,
        all64_actual_products_annihilate_H20_and_raw26=True,
        raw26_coordinates_commute=True,closed_raw_module_basis_dimension=26,
        actual_raw26_scalar_callback_complete_for_supplied_H20_vectors=True,
        frame='Actual CA-F9; r used on both A3 sectors; s=r only after ell=0',
        no_nu_a_b_delta_evaluation_rawphysical_substitution_DVR_action_or_comparison_claim=True)
    (out/'finite-identities.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{STEM}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=STEM,all_passed=True,actual_closed_raw_cokernel_dimension=26,
        actual_coordinate_actions8_by26_by26=True,actual_H20_quotient_maps_both3_by20=True,
        no_new_source_jet_inverse_solve_GB_point_prime_or_rank_search=True,
        no_actual_comparison_kernel_cokernel_torsion_or_D0_depth_claim=True,
        input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        frozen=True,independent_review_pending=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
