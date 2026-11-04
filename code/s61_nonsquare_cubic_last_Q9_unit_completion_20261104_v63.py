#!/usr/bin/env python3
"""Final fixed GL10 completion; requires the independent POWER review gate."""
import argparse
import hashlib
import json
import resource
import signal
import time
from pathlib import Path

import sympy as sy
from sympy.polys.rings import ring
from s61_nonsquare_cubic_zw25_pair_unit_completion_20261104_v63 import Lazy, matmul

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def universal():
    tup=ring(','.join([f'a{i}'for i in range(9)]+[f'c{i}'for i in range(9)]),sy.ZZ)
    R=tup[0];aa=list(tup[1:10]);cc=list(tup[10:]);z=R.zero;one=R.one
    mm=lambda A,B:matmul(A,B,z,lambda x,y:x+y,lambda x,y:x*y)
    t=sum((x*y for x,y in zip(aa,cc)),z)
    B=[[(one if i==j else z)-cc[i]*aa[j]for j in range(9)]for i in range(9)]
    U=[B[i]+[cc[i]]for i in range(9)]+[aa+[z]]
    W=[[cc[i]*aa[j]for j in range(9)]+[-cc[i]]for i in range(9)]+[[-x for x in aa]+[one]]
    UU=mm(U,U)
    for i in range(10):
        for j in range(10):
            assert UU[i][j]-(one if i==j else z)==(t-one)*W[i][j]
    assert mm([aa+[z]],U)[0]==[(one-t)*x for x in aa]+[t]
    # Universal adjugate/rank-one determinant proof, no 10! expansion.
    adj=[[(one-t if i==j else z)+cc[i]*aa[j]for j in range(9)]for i in range(9)]
    for product in[mm(B,adj),mm(adj,B)]:
        assert product==[[(one-t if i==j else z)for j in range(9)]for i in range(9)]
    bordered=-sum((aa[i]*adj[i][j]*cc[j]for i in range(9)for j in range(9)),z)
    assert bordered==-t
    return {'integer_variables':[str(x)for x in aa+cc],
            'dot_t':'sum(a_i*c_i)',
            'all100_U_square_remainders':'U^2-I=(t-1)*[[c*a,-c],[-a,1]]',
            'row_identity':'[a,0]*U=[(1-t)*a,t]',
            'all162_rankone_adjugate_products_verified':True,
            'det_B_proof':'Column multilinearity of I-c*a: all terms changing two columns have proportional c columns and vanish; the base plus nine singleton terms equals1-t.',
            'adj_B':'(1-t)*I+c*a',
            'bordered_determinant_polynomial':'-a*adj(B)*c=-t',
            'at_actual_dot1':'U is an involution with determinant-1; no free basis of ker(a) is assumed.'}


def main(out,power_review):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text())
    p=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63'
    m=ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63'
    q=ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_PREFLIGHT_20261104_v63'
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_cubic_zw25_pair_unit_completion_20261104_v63.py',
           p/'report.json',p/'certificate.json',m/'report.json',m/'certificate.json',
           q/'report.json',q/'certificate.json',power_review]
    pr,pc,mr,mc,qr,qc,peer=[read(x)for x in paths[2:]]
    assert peer['all_passed'] and peer['frozen']
    assert peer['input_sha256'][str(paths[3].relative_to(ROOT))]==sha(paths[3])
    assert sha(paths[3])==pr['certificate_sha256'] and sha(paths[5])==mr['certificate_sha256']
    assert sha(paths[7])==qr['certificate_sha256']
    assert pr['literal_finite_DAG_Bezout_target']=='a^29'
    original_rows=mc['new15_source_original25_labels'];original_cols=mc['new22_target_original32_columns']
    assert original_rows[0]==9 and len(original_rows)==15 and len(original_cols)==22
    selected=pc['actual_h9_original_target_labels'];assert selected==qc['query_only_original32_target_columns']
    assert selected==[11,13,14,15,17,18,19,22,25]
    def asdict(poly):return{tuple(m):sy.QQ(c)for m,c in poly}
    for f,v in zip(pc['actual_h9'],qc['actual_nine_coefficients']):
        n,d=v['actual_coefficient'];assert len(d)==1 and d[0][0]==[0,0,0]
        assert asdict(f)=={tuple(m):sy.QQ(c)/sy.QQ(d[0][1])for m,c in n}
    Aold=mc['remaining15x22_scalar_roots']
    assert [Aold[0][original_cols.index(j)]for j in selected]==qc['queried_roots']
    # All previously expanded core arithmetic stays external and immutable.
    # The new circuit contains only finite references and scalar operations.
    lazy=Lazy()
    matrix=lambda root:0 if not root else lazy.put(['matrix_parent',root])
    power=lambda root:0 if not root else lazy.put(['power_parent',root])
    inversea=lazy.put(['known_inverse_a29'])
    A=[[matrix(root)for root in row]for row in Aold]
    h=[A[0][original_cols.index(j)]for j in selected]
    c=[lazy.mul(power(root),inversea)for root in pc['actual_h9_cofactor_roots']]
    inner=0
    for x,y in zip(h,c):inner=lazy.add(inner,lazy.mul(x,y))
    proof=universal()
    # Clear all thirteen target tails with the exact known row functional.
    selected_positions=[original_cols.index(j)for j in selected]
    tails=[j for j in original_cols if j not in selected]
    assert len(tails)==13
    heavy=[27,28,29,30,31]
    assert all(Aold[0][original_cols.index(j)]==0 for j in heavy)
    T,Ti=lazy.eye(22),lazy.eye(22)
    for j in tails:
        bj=A[0][original_cols.index(j)]
        for i,ci in enumerate(c):
            value=lazy.mul(ci,bj)
            T[selected_positions[i]][original_cols.index(j)]=lazy.sub(0,value)
            Ti[selected_positions[i]][original_cols.index(j)]=value
    assert set(selected).isdisjoint(tails)
    spare=1;assert spare in tails and spare not in heavy
    block=selected+[spare];positions=[original_cols.index(j)for j in block]
    U=[ [lazy.sub(int(i==j),lazy.mul(c[i],h[j]))for j in range(9)]+[c[i]]for i in range(9)]
    U.append(h+[0])
    U22=lazy.eye(22)
    for i,ri in enumerate(positions):
        for j,cj in enumerate(positions):U22[ri][cj]=U[i][j]
    target=lazy.mm(T,U22);targeti=lazy.mm(U22,Ti)
    transformed=lazy.mm(A,target);pivot=original_cols.index(spare)
    source=lazy.eye(15)
    for i in range(1,15):source[i][0]=lazy.sub(0,transformed[i][pivot])
    rows=original_rows[1:];cols=[j for j in original_cols if j!=spare]
    final=[[transformed[i][original_cols.index(j)]for j in cols]for i in range(1,15)]
    for i in range(1,15):
        for j in heavy:
            assert transformed[i][original_cols.index(j)]==A[i][original_cols.index(j)]
    cert={'accepted_power_inner_product_gate':{'power_certificate_sha256':sha(paths[3]),
             'independent_review_sha256':sha(power_review),'inner_root':inner,
             'identity':'sum h_i*(H_i/a29)=1 by pinned actual-coordinate POWER identity'},
          'universal_integer_proof':proof,
          'scalar_DAG':{'nodes':lazy.nodes,
            'external_matrix_certificate':str(paths[5].relative_to(ROOT)),
            'matrix_parent_node_definition':'parent.scalar_DAG.nodes[root]',
            'external_power_certificate':str(paths[3].relative_to(ROOT)),
            'power_parent_node_definition':'parent.scalar_DAG.nodes[root]',
            'known_inverse_a29_definition':'1/a^29, certified FULLU unit'},
          'actual_h9_roots':h,'actual_c9_dual_roots':c,
          'selected_original_target_labels':selected,'all13_tail_original_labels':tails,
          'tail_clear22_roots':T,'tail_clear22_inverse_roots':Ti,
          'tail_inverse_proof':'T=I-N,Ti=I+N; support(N) has rows in h9 and columns in disjoint tails13, henceN^2=0 identically.',
          'pure_spare_original_target_label':spare,
          'GL10_original_target_block_order':block,'GL10_roots':U,
          'target22_roots':target,'target22_inverse_roots':targeti,
          'target_inverse_proof':'U^2=I by the pinned inner=1 and universal identity; tail products exact by disjoint support. (T*U)^-1=U*Ti.',
          'source15_roots':source,'source_change_determinant':'1','target_change_determinant':'-1',
          'old15_source_labels':original_rows,'old22_target_labels':original_cols,
          'new14_source_labels':rows,'new21_target_labels':cols,
          'pivot_row_before_exact_simplification_roots':transformed[0],
          'pivot_row_exact_simplified':[int(i==pivot)for i in range(22)],
          'actual14x21_scalar_roots':final,'all5_heavy_quotient_labels_unchanged':heavy,
          'new_pivot_value':'1','new_pivot_source_label':9,'new_pivot_target_label':spare,
          'total_global_unit_pivots_in_augmented25_presentation':11,
          'same_cokernel_proof':'After tail clearing row9=[h9,0^13]. Embedded U sends it to a unit vector in the spare pure target. The explicit source15 clearing deletes that row/column; all other full row coefficients remain in the saved14x21 roots.',
          'no_original_fixed11minor_no_full57_or_quartic239_claim':True}
    (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
    report={'id':ID,'independent_POWER_gate_bound_to_actual9_matrix_entries':True,
            'all100_universal_involution_and_row_polynomial_checks_passed':True,
            'universal_determinant_minus1_rankone_proof_verified':True,
            'all13_tails_and_all15_rows_transported_finite_DAG':True,
            'actual_shape_same_cokernel':[14,21],
            'new_global_unit_pivots':1,'total_global_pivots_in_transformed_presentation':11,
            'all5_heavy_quotient_coefficients_unchanged_literal_roots':True,
            'core_25k_power_DAG_reused_external_never_expanded':True,
            'new_inverse_only_a29_knownunit':True,
            'no_Gaussian_point_GB_minorsearch_or_automatic_nextstep':True,
            'no_full57_quartic239_or_kernel38_gap42_change_claim':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'certificate_sha256':sha(out/'certificate.json'),'new_scalar_nodes':len(lazy.nodes),
            'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--power-review',type=Path,required=True)
    args=ap.parse_args();main(args.output.resolve(),args.power_review.resolve())
