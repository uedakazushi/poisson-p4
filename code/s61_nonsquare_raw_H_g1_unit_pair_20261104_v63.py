#!/usr/bin/env python3
"""One fixed g1-branch pair/inverse check. No minor search or rank solver."""
import argparse, hashlib, json, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_NONSQUARE_RAW_H_G1_UNIT_PAIR_20261104_v63'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()
def add(p,q,c=Q(1)):
    r=p.copy()
    for m,v in q.items():
        r[m]=r.get(m,Q())+c*v
        if not r[m]: del r[m]
    return r
def mul(p,q):
    r={}
    for m,v in p.items():
        for n,w in q.items():
            s=tuple(x+y for x,y in zip(m,n));r[s]=r.get(s,Q())+v*w
    return {m:v for m,v in r.items() if v}
def sm(p,c): return {m:c*v for m,v in p.items() if c*v}
def shift(p,s,c=Q(1)):
    return {tuple(x+y for x,y in zip(m,s)):c*v for m,v in p.items() if c*v}
def power(p,k):
    out={(0,0,0):Q(1)}
    for _ in range(k): out=mul(out,p)
    return out
def dec(p): return {tuple(m):Q(v) for m,v in p if Q(v)}
def enc(p): return [[list(m),str(v)] for m,v in sorted(p.items())]
def matmul(A,B):
    return [[add(mul(A[i][0],B[0][j]),mul(A[i][1],B[1][j]))
             for j in range(2)] for i in range(2)]

def main(out):
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63'
    raw=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63'
    paths=[Path(__file__),base/'report.json',base/'certificate.json',raw/'report.json',raw/'checkpoint.json',
           ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63_REVIEW/report.json']
    rd=lambda p:json.loads(p.read_text())
    report,cert,rawreport,checkpoint=[rd(p) for p in paths[1:5]]
    assert sha(paths[2])==report['certificate_sha256']
    assert sha(paths[4])==rawreport['evidence_sha256']['checkpoint.json']
    assert rd(paths[5])['all_passed']
    rawM=[[dec(p) for p in row] for row in checkpoint['matrix15x6']]
    branch=next(b for b in cert['branches'] if b['branch']=='g1')
    records={(r['raw_checkpoint_row'],r['column']):r for r in branch['full14x2_records']}
    a={(1,0,0):Q(1)};e={(0,1,0):Q(1)};d={(0,0,1):Q(1)};one={(0,0,0):Q(1)}
    T=add(add(mul(a,power(d,2)),sm(d,-2)),sm(one,8))
    g1=dec(cert['g1'])
    B=add(mul(add(add(add(sm(a,50),sm(power(e,2),18)),sm(e,-9)),sm(one,-2)),d),add(sm(e,60),sm(one,-40)))
    # Literal remainder formulas, not interpolation or a field-rank test.
    M=[[shift(power(T,3),(0,0,-3),Q(1,100)),shift(power(T,4),(0,0,-4),Q(1,10))],
       [sm(power(T,2),Q(1,300)),shift(power(T,3),(0,0,-1),Q(1,6))]]
    selected=[]
    for i,row in enumerate([5,8]):
        assert row in checkpoint['remaining_rows']
        assert all(not rawM[row][j] for j in [2,3,4])
        for j,col in enumerate([0,5]):
            rec=records[row,col];quot=dec(rec['quotient_by_original_branch_polynomial'])
            assert M[i][j]==dec(rec['remainder'])
            assert add(mul(quot,g1),M[i][j])==rawM[row][col]
            selected.append({'row':row,'column':col,'original':enc(rawM[row][col]),
                             'quotient_by_g1':enc(quot),'remainder':enc(M[i][j])})
    # B = 50T/d modulo g1, with explicit original-ring quotient.
    bq=add(add(sm(e,6),sm(one,1)),shift(one,(0,0,-1),Q(-20)))
    assert add(B,shift(T,(0,0,-1),Q(-50)))==mul(g1,bq)
    det=add(mul(M[0][0],M[1][1]),mul(M[0][1],M[1][0]),Q(-1))
    assert det==shift(power(T,6),(0,0,-4),Q(1,750))
    # det = B^6*d^2/(750*50^6) modulo g1: difference-of-powers witness.
    V=shift(T,(0,0,-1),Q(50));diffq={}
    for k in range(6): diffq=add(diffq,mul(power(B,5-k),power(V,k)))
    dq=shift(mul(bq,diffq),(0,0,2),Q(-1,750*50**6))
    unitdet=shift(power(B,6),(0,0,2),Q(1,750*50**6))
    assert add(det,unitdet,Q(-1))==mul(g1,dq)
    adj=[[M[1][1],sm(M[0][1],Q(-1))],[sm(M[1][0],Q(-1)),M[0][0]]]
    for product in [matmul(M,adj),matmul(adj,M)]:
        assert all(product[i][j]==(det if i==j else {}) for i in range(2) for j in range(2))
    # Inverse entries are saved as common numerator/denominator. Its denominator
    # is a unit by the preceding quotient identity, not an unknown localization.
    inverse=[[{'numerator':enc(adj[i][j]),'common_denominator':enc(det)} for j in range(2)] for i in range(2)]
    proof={'T':enc(T),'g1':enc(g1),'B':enc(B),'selected_source_rows':[5,8],'selected_heavy_columns':[0,5],
           'other_heavy_columns':[2,3,4],'all_six_other_heavy_entries_literally_zero':True,
           'selected_full_source_quotient_identities':selected,
           'B_minus_50T_over_d_quotient_by_g1':enc(bq),
           'pair_matrix':[[enc(p) for p in row] for row in M],
           'determinant':enc(det),'unit_determinant_representative':enc(unitdet),
           'det_minus_unitrepresentative_quotient_by_g1':enc(dq),
           'unit_formula':'det=B^6*d^2/(750*50^6) in R/(g1)',
           'inverse':inverse,'both_inverse_products_verified_after_exact_denominator_clearing':True,
           'block_reduction':'On FULLU/(g1), this pair occupies only first/last columns and is invertible; eliminating these columns from the other rows leaves exactly the unchanged middle B7x3 block. Thus rank(heavy14x5)=2+rank(B7x3) at every residue field on this branch.',
           'conditional_scope':'Any independently proved finite-or-empty Bbad cover restricts the possible heavy-rank defect on the g1 branch. No Bbad existence, no global B rank3, no Theta result, no full cubic/quartic rank claimed.'}
    (out/'certificate.json').write_text(json.dumps(proof,separators=(',',':'))+'\n')
    result={'id':ID,'all_four_source_remainders_verified':True,'all_four_source_quotient_identities_verified':True,
            'middle_three_coefficients_of_both_pair_rows_literally_zero':True,
            'determinant_and_known_unit_representative_verified':True,'both_inverse_products_verified':True,
            'no_minor_search_GB_solver_new_point_prime_or_other_branch_operation':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-started,
            'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())
