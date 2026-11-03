#!/usr/bin/env python3
"""Independent corrected shear, full28 polynomial kernel and physical Hessian."""
import argparse
from pathlib import Path
import json
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
from s61_model import X,bivector_from_coefficients,pushforward
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha
from q_positive_weight6_small_model_review_20261103_v63 import physical_coefficients,own_action


def main(replay,output):
    directory=ROOT/'reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63'
    d=read(directory/'report.json');fresh=read(replay/'report.json')
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
    assert (directory/'certificate.json.gz').read_bytes()==(replay/'certificate.json.gz').read_bytes()
    for filename,digest in d['input_sha256'].items():assert sha(ROOT/filename)==digest,filename
    c=read(directory/'certificate.json.gz');assert sha(directory/'certificate.json.gz')==d['certificate_sha256']
    family=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz')
    base=read(ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz')
    n0,n1,n2,r=s.symbols('n0 n1 n2 r')
    full=s.Matrix(list(map(s.sympify,family['canonical150'])))
    P=bivector_from_coefficients(full)
    # Independently test the ENTIRE canonical homogeneous tensor, rather
    # than a diagonal firstjet or one affine bracket.
    for epsilon in [1,-1]:
        M=s.diag(1,r,r**2,r**3,epsilon*r);M[3,0]=(r**3-epsilon)/(6*n0)
        assert s.factor(M.det())==epsilon*r**7
        changed=epsilon/r**2*pushforward(P,M)
        target=P.subs({n0:epsilon*n0,n1:epsilon*r*n1,n2:n2/r},simultaneous=True)
        assert (changed-target).applyfunc(s.cancel)==s.zeros(5)
    q=full.subs(n1,1);C=s.Matrix(base['canonical_basis']);free=base['canonical_free_raw_indices']
    T=s.Matrix(c['extra_normal126_by2']);physical=C*T
    assert physical==s.Matrix(c['extra_normal_physical150_by2'])
    actual=own_action(q).row_join(q.diff(n0)).row_join(q.diff(n2))[free,:]
    selected=c['actual_family_selected26_columns']
    K=actual[:,selected].row_join(T)
    assert K==s.Matrix(c['full126_by28_polynomial_kernel'])
    R,*variables=ring('n0,n2',s.QQ)
    qr=[R.from_expr(v)for v in q];crows=[[(j,s.QQ.convert(C[i,j]))for j in range(126)if C[i,j]]for i in range(150)]
    constant,_=physical_coefficients();J=[]
    for terms in constant:
        raw={}
        for(a,b),value in terms.items():
            raw[a]=raw.get(a,R.zero)+value*qr[b];raw[b]=raw.get(b,R.zero)+value*qr[a]
        row={}
        for i,value in raw.items():
            if value:
                for j,coefficient in crows[i]:row[j]=row.get(j,R.zero)+coefficient*value
        J.append({j:v for j,v in row.items()if v})
    Kr=[[R.from_expr(v)for v in row]for row in K.tolist()]
    for row in J:
        for j in range(28):assert sum((value*Kr[i][j]for i,value in row.items()),R.zero)==R.zero
    sample={n0:1,n2:1}
    Kp=K.subs(sample)
    assert Kp[c['kernel28_rows'],:].det()==s.Rational(c['kernel28_minor'])!=0
    Jp=s.zeros(350,126)
    for i,row in enumerate(J):
        for j,value in row.items():Jp[i,j]=sum(s.QQ.to_sympy(v)for v in value.values())
    determinant=DomainMatrix.from_Matrix(Jp[c['Jacobian98_rows'],c['Jacobian98_columns']]).convert_to(s.QQ).det()
    assert s.Rational(determinant)==s.Rational(c['Jacobian98_minor'])!=0
    Tr=[[R.from_expr(v)for v in row]for row in physical.tolist()]
    H=[]
    for terms in constant:
        row=[R.zero]*3
        for(a,b),value in terms.items():
            row[0]+=value*Tr[a][0]*Tr[b][0]
            row[1]+=value*(Tr[a][0]*Tr[b][1]+Tr[a][1]*Tr[b][0])
            row[2]+=value*Tr[a][1]*Tr[b][1]
        H.append(row)
    for row,saved in zip(H,c['full350_by3_Hessian_columns']):
        assert all(value==R.from_expr(s.sympify(text))for value,text in zip(row,saved))
    paths=[Path(__file__),directory/'report.json',directory/'certificate.json.gz',
        ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.py',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json',
        ROOT/'code/s61_model.py',ROOT/'code/q_positive_weight6_quadratic_model_20261103_v63.py',
        ROOT/'research/notes/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63_REVIEW.md']
    result=dict(id='Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63_REVIEW',all_passed=True,
        fresh_short_producer_and_complete_certificate_match=True,
        independent_full150_canonical_finite_sheared_transform_both_signs=True,
        dense_n1_equals1_normalization_on_D_n0_n1_accepted=True,
        independent_all26_actual_group_scalar_two_family_tangent_columns=True,
        independent_physical350_by126_QQn0n2_Jacobian=True,
        independent_all350_by28_universal_polynomial_kernel_products_zero=True,
        exact28kernel_minor_and98Jacobian_minor_checked=True,
        exact_generic_full_canonical_Jacobian_rank=98,native192second_rank_not_used=True,
        independent_all350_by3_physical_quadratic_Hessian_coefficients=True,
        generic_two_normal_tangent_quotient_basis_accepted=True,
        generic_quadratic_cokernel_rank_or_exact_Artin_length_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
    with output.open('x')as f:
        json.dump(result,f,indent=2);f.write('\n')
    print('PASS independent sheared normalization+28polynomial kernel+98minor; no nativeRF needed')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();main(args.replay,args.output)
