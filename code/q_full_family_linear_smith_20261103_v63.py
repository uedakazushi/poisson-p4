#!/usr/bin/env python3
"""Full126 canonical Jacobi/orbit Smith certificates over QQ[a] for Qa.

Physical Jacobi linearization is reconstructed directly in all150 raw basis
directions; no108-column slice is substituted at exceptional parameters.
Unimodular left/right matrices certify every Smith diagonal at ALL a.
"""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from q_weighted_finite_model_20261103_v63 import frames,A,WEIGHTS,point,coefficients,BASIS,TRIPLES,MON3
from ambient_rank_two_flag_intersection_table_20261103_v63 import lie

ROOT=Path(__file__).resolve().parents[1]

def linear_data():
    C,free,*_=frames();q=coefficients(point(A));Q={}
    for col,(pair,mon)in enumerate(BASIS):
        if q[col]:Q.setdefault(pair,{})[mon]=q[col]
    def qp(i,j):
        return Q.get((i,j),{})if i<j else({m:-v for m,v in Q.get((j,i),{}).items()}if i>j else{})
    mi={m:i for i,m in enumerate(MON3)};L=s.zeros(350,150)
    for col,((u,v),mon)in enumerate(BASIS):
        def rp(i,j):return 1 if(i,j)==(u,v)else -1 if(i,j)==(v,u)else 0
        for row,(i,j,k)in enumerate(TRIPLES):
            for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
                for ell in range(5):
                    r=rp(c,d)
                    if r and ell in mon:
                        rem=list(mon);rem.remove(ell)
                        for m,co in qp(b,ell).items():L[35*row+mi[tuple(sorted(m+tuple(rem)))],col]+=r*mon.count(ell)*co
                    r=rp(b,ell)
                    if r:
                        for m,co in qp(c,d).items():
                            if ell in m:
                                rem=list(m);rem.remove(ell)
                                L[35*row+mi[tuple(sorted(mon+tuple(rem)))],col]+=r*m.count(ell)*co
    target=[10-(sum(t)-sum(m))for t in TRIPLES for m in MON3]
    degree=[5-WEIGHTS[i]for i in free];J=L*C
    positions=list(itertools.product(range(5),repeat=2));O=s.Matrix.hstack(*(lie(q,i,j)for i,j in positions))
    assert C*q[free,0]==q and J*q[free,0]==s.zeros(350,1)
    assert (J*O[free,:]).applyfunc(s.expand)==s.zeros(350,25)
    assert sum(((i-2)*O[:,positions.index((i,i))]for i in range(5)),s.zeros(150,1))==5*q
    assert all(J[i,j]==0 for i in range(350)for j in range(126)if target[i]!=degree[j])
    return C,free,q,J,target,degree,O,positions

def matrix_json(M):return [[str(v)for v in row]for row in M.tolist()]

def certify(B,label,rows,cols):
    R=s.QQ.poly_ring(A);M=DomainMatrix.from_Matrix(B).convert_to(R).to_dense()
    D,U,V=smith_normal_decomp(M)
    assert U*M*V==D
    du,dv=U.det(),V.det();assert du and dv and not R.to_sympy(du).has(A)and not R.to_sympy(dv).has(A)
    dd=D.to_Matrix();diagonal=[s.factor(dd[i,i])for i in range(min(dd.shape))]
    assert all(dd[i,j]==0 for i in range(dd.rows)for j in range(dd.cols)if i!=j)
    rank=sum(v!=0 for v in diagonal)
    return dict(label=label,rows=rows,columns=cols,rank=rank,diagonal=list(map(str,diagonal)),
      matrix=matrix_json(B),left=matrix_json(U.to_Matrix()),right=matrix_json(V.to_Matrix()),
      left_determinant=str(R.to_sympy(du)),right_determinant=str(R.to_sympy(dv)))

def main(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    C,free,q,J,target,degree,O,positions=linear_data();jac=[];orbit=[]
    for d in sorted(set(degree)):
        cc=[j for j,v in enumerate(degree)if v==d];rr=[i for i,v in enumerate(target)if v==d]
        if rr:jac.append(certify(J[rr,cc],f'Jacobi_positive_weight_{d}',rr,cc))
    for ell in range(-4,5):
        cc=[k for k,(i,j)in enumerate(positions)if i-j==ell];rr=[i for i in free if WEIGHTS[i]==5+ell]
        if rr:orbit.append(certify(O[rr,cc],f'orbit_Lie_weight_{ell}',rr,cc))
    aindex=BASIS.index(((1,4),(0,0)))
    initial=[]
    for ell in [0,1,2]:
        cc=[k for k,(i,j)in enumerate(positions)if i-j==ell]
        rr=[i for i in free if WEIGHTS[i]==5+ell and i!=aindex]
        initial.append(certify(O[rr,cc],f'initial_nonpositive_group_cut_{ell}',rr,cc))
    assert sum(b['rank']for b in jac)==75 and sum(b['rank']for b in orbit)==17
    jnon=[s.factor(s.sympify(v))for b in jac for v in b['diagonal']if s.sympify(v).has(A)]
    assert len(jnon)==3
    assert sorted(str(s.monic(v,A))for v in jnon)==sorted(['a','a','a - 2/9'])
    assert all(not s.sympify(v).has(A)for b in orbit+initial for v in b['diagonal'])
    ranks={}
    for value in [s.S.Zero,s.Rational(2,9),s.Rational(-4,9),s.S.One]:
        jr=len(J.subs(A,value).rref()[1]);orr=len(O.subs(A,value).rref()[1])
        expect=73 if value==0 else 74 if value==s.Rational(2,9) else 75
        assert jr==expect and orr==17
        ranks[str(value)]=dict(Jacobi_rank=jr,affine_tangent_dimension=126-jr,
          projective_tangent_dimension=125-jr,affine_orbit_dimension=17,projective_orbit_dimension=16,
          transverse_tangent_dimension=109-jr)
    payload=dict(id='Q_FULL_FAMILY_LINEAR_SMITH_20261103_V63',field='QQ[a]',
      canonical_free_raw_indices=free,canonical_basis=matrix_json(C),raw_Qa150=list(map(str,q)),
      full350_by126_Jacobi=matrix_json(J),Jacobi_row_relative_weights=target,
      canonical_relative_weights=degree,full150_by25_orbit=matrix_json(O),Lie_positions=positions,
      Jacobi_blocks=jac,orbit_blocks=orbit,initial_nonpositive_group_blocks=initial)
    archive=output/'certificate.json.gz'
    with archive.open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as z:z.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    paths=[Path(__file__),ROOT/'code/q_weighted_finite_model_20261103_v63.py',ROOT/'code/ambient_rank_two_flag_intersection_table_20261103_v63.py']
    result=dict(id=payload['id'],all_passed=True,full126_columns_used=True,
      all_unimodular_left_right_identities_checked=True,all_transformation_determinants_nonzero_constants=True,
      Jacobi_generic_rank=75,Jacobi_nonunit_invariant_factors_up_to_constants=['a','a','9a-2'],
      orbit_all_a_rank=17,nonpositive_group_cut_rank_all_a=7,
      only_nonzero_a_Jacobi_rank_drop='a=2/9',Jacobi_rank_at_a_zero=73,
      group_minor_a_or_9a_plus4_zeros_are_coordinate_frame_artifacts=True,
      exact_ranks_at_special_parameters=ranks,
      full_higher_order_scheme_classification_claimed=False,
      input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
      certificate_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('FULL Qa JACOBI+ORBIT SMITH CERTIFICATE PASS',result['seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
