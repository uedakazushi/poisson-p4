#!/usr/bin/env python3
"""Necessary ordinary CA paired-flag3jets atB18, all9 flags andCA20/free10."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from ambient_CA_V1_independent_review_20261103_v63 import tensor_factor_operator,act,rational

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63'


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def mm(A,B):return(DomainMatrix.from_Matrix(A).convert_to(s.QQ)*DomainMatrix.from_Matrix(B).convert_to(s.QQ)).to_Matrix()
def piv(M):return list(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(60)
    source=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';data=read(source)
    accepted=ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json';assert read(accepted)['all_passed']
    C=s.Matrix(data['fixed_fiber_linear_basis150_by20']);free=data['affine_slice_free_indices'];assert free==[1,2,3,9,10,14,15,17,18,19]
    P0=list(map(rational,C*s.Matrix(data['center_fixed20'])));F=[list(map(rational,C[:,j]))for j in free]
    fullflag=list(map(tuple,data['flag_positions']));assert len(fullflag)==9 and (1,0)in fullflag
    flag=[v for v in fullflag if v!=(1,0)];assert len(flag)==8
    ops=[tensor_factor_operator(*divmod(i,5))for i in range(25)]
    M=s.zeros(5);M[1,0]=1;M[0,3]=-2
    E=[]
    for i,j in flag:
        B=s.zeros(5);B[i,j]=1;E.append(B)
    bs=s.symbols('b0:8');cs=s.symbols('c0:8')
    genericB=sum((v*B for v,B in zip(bs,E)),s.zeros(5));genericC=sum((v*B for v,B in zip(cs,E)),s.zeros(5))
    assert(M+genericB)**4==s.zeros(5)
    log2=genericB-M*M/2;log3=genericC-(M*genericB+genericB*M)/2+M**3/3
    assert(log2+M*M/2-genericB).applyfunc(s.expand)==s.zeros(5)
    assert(log3+(M*log2+log2*M)/2+M**3/6-genericC).applyfunc(s.expand)==s.zeros(5)
    def A(matrix,vector):
        answer=[Q(0)]*150
        for i,j in itertools.product(range(5),repeat=2):
            if matrix[i,j]:
                coefficient=rational(matrix[i,j]);value=act(ops[5*i+j],vector)
                for k,c in enumerate(value):
                    if c:answer[k]+=coefficient*c
        return answer
    LM=lambda v:A(M,v);LM2=lambda v:A(M*M,v)
    def K2(v):return[(a-b)/2 for a,b in zip(LM(LM(v)),LM2(v))]
    def K30(v):
        return[a/3-(b+c)/4+d/6 for a,b,c,d in zip(A(M**3,v),LM(LM2(v)),LM2(LM(v)),LM(LM(LM(v))))]
    def KB(B,v):
        return[-a/2+(b+c)/2 for a,b,c in zip(A(M*B+B*M,v),LM(A(B,v)),A(B,LM(v)))]
    selected=piv(C.T);left=C[selected,:].inv()
    def project(V):return V-mm(C,mm(left,V[selected,:]))
    assert project(s.Matrix(LM(P0)))==s.zeros(150,1)
    D=project(s.Matrix.hstack(*[s.Matrix(LM(v))for v in F],*[s.Matrix(A(B,P0))for B in E]))
    constant=project(s.Matrix(K2(P0)));columns=piv(D);rows=piv(D[:,columns].T)
    inverse=D[rows,columns].inv();assert mm(inverse,D[rows,columns])==s.eye(len(columns))
    coker_constant=constant-mm(D[:,columns],mm(inverse,constant[rows,:]))
    consistent=coker_constant==s.zeros(150,1)
    free_columns=[i for i in range(18)if i not in columns]
    setup=dict(CA_fixed_fiber_dimension=20,CA_parabolic_free10=free,full_relative_flag9=fullflag,
        normalized_flag_coordinate=[1,0],tangent_M=mat(M),second_and_third_order_flag8=flag,
        all150_outside_fiber_order2_linear_matrix=mat(D),all150_order2_constant=list(map(str,constant)),
        rank_order2=len(columns),pivot_columns=columns,pivot_rows=rows,inverse_pivot=mat(inverse),
        ordinary_order2_consistent=consistent,constant_coker=list(map(str,coker_constant)))
    (out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n')
    if not consistent:
        result=dict(id=STEM,all_passed=True,ordinary_order2_consistent=False,
            no_any_order_or_nonordinary_pair_exclusion_claimed=True,seconds=time.monotonic()-start)
    else:
        d=s.zeros(18,1);T=s.zeros(18,len(free_columns))
        for j,col in enumerate(free_columns):T[col,j]=1
        solved=-mm(inverse,constant[rows,:]);mixed=-mm(inverse,D[rows,free_columns])
        for j,col in enumerate(columns):d[col]=solved[j];T[col,:]=mixed[j,:]
        assert mm(D,d)+constant==s.zeros(150,1)and mm(D,T)==s.zeros(150,len(free_columns))
        # Order3 unknownhigherP2/free10 andhigherflagC8 span EXACTLY the sameD.
        labels=[(),*[(i,)for i in range(18)],*[(i,10+j)for i in range(10)for j in range(8)]]
        coefficient_vectors=[s.Matrix(K30(P0)),*[s.Matrix(K2(v))for v in F],*[s.Matrix(KB(B,P0))for B in E],*[s.Matrix(A(B,F[i]))for i in range(10)for B in E]]
        R3=project(s.Matrix.hstack(*coefficient_vectors))
        R3c=R3-mm(D[:,columns],mm(inverse,R3[rows,:]))
        f=len(free_columns);names=[f'v{i}'for i in range(f)]
        if f:
            R,*vv=ring(','.join(names),s.QQ)
            variables=[R(s.Rational(d[i]))+sum((R(s.Rational(T[i,j]))*vv[j]for j in range(f)),R.zero)for i in range(18)]
            words=[m for degree in range(3)for m in itertools.combinations_with_replacement(range(f),degree)]
            expansion=s.zeros(len(labels),len(words))
            for i,label in enumerate(labels):
                polynomial=R.one
                for j in label:polynomial*=variables[j]
                for k,word in enumerate(words):expansion[i,k]=s.QQ.to_sympy(polynomial.get(tuple(word.count(j)for j in range(f)),s.QQ.zero))
        else:
            vv=[];words=[()];expansion=s.Matrix([s.prod(d[i]for i in label)for label in labels])
        residual=mm(R3c,expansion);eqrows=piv(residual.T);basis=residual[eqrows,:]
        syms=s.symbols(' '.join(names))if f else[]
        if f==1:syms=(syms,)
        polynomials=[s.factor(sum(c*s.prod(syms[i]for i in word)for c,word in zip(row,words)if c))for row in basis.tolist()]
        gb=[];status='not_run_more_than4_free_variables'
        if f<=4:
            if f and polynomials:gb=list(s.groebner(polynomials,*syms,domain=s.QQ,order='grevlex'));status='tiny_QQ_complete'
            else:gb=polynomials;status='constant_or_zero'
        extra={}
        if f==3:
            v0,v1,v2=syms;conic=4*v0*v1+4*v0+7*v1*v1+18*v1+8;H=4*v0+7*v1+11
            expected=[2*v2*(2*v0+3*v1+5)/5,-conic/2,-v2*(v0+4*v1+5)/5]
            assert all(s.expand(a-b)==0 for a,b in zip(polynomials,expected))and len(polynomials)==3
            cof=[-H/6,2*v2/3,-2*H/3]
            assert s.expand(sum(a*b for a,b in zip(cof,polynomials))-v2)==0
            assert s.expand(-2*polynomials[1]-conic)==0
            # One exact finite3jet witness; it is NOT an all-orders pair.
            witness=s.Matrix([-2,0,0]);low=d+T*witness
            products=s.Matrix([s.prod(low[i]for i in label)for label in labels])
            intrinsic=mm(R3,products);high=s.zeros(18,1);hsol=-mm(inverse,intrinsic[rows,:])
            for i,col in enumerate(columns):high[col]=hsol[i]
            assert mm(D,low)+constant==s.zeros(150,1)and mm(D,high)+intrinsic==s.zeros(150,1)
            assert T[1,:]==s.zeros(1,3)and 1 in columns
            functional=-mm(inverse[columns.index(1),:],R3[rows,:]);affine=d+T*s.Matrix(syms)
            second_extra=s.factor(sum(c*s.prod(affine[i]for i in label)for c,label in zip(functional,labels)))
            assert s.expand(second_extra-(v1+2)*(3*v0+7*v1+5)/2)==0
            assert high[1]==-1 #the extraCA coefficient is genuinely used in this witness
            extra=dict(exact_coker_ideal_generators=[str(v2),str(conic)],
                v2_as_cofactor_combination=list(map(str,cof)),conic_equals_minus2_second_generator=True,
                smooth_necessary_jet_curve_parameterization='v1=t-1,v0=-7*t/4-1+3/(4*t),v2=0; t!=0',
                unique_higher_first_source_CA_extra_C2_polynomial=str(second_extra),
                finite3jet_witness_free_values=['-2','0','0'],finite3jet_first_source_rho10=list(map(str,low[:10,:])),
                finite3jet_second_order_flag_B8=list(map(str,low[10:,:])),
                finite3jet_higher_first_source_P2_free10=list(map(str,high[:10,:])),
                finite3jet_third_order_flag_C8=list(map(str,high[10:,:])),
                all150_membership_order2_and3_witness_verified=True,
                witness_first_source_extra_CA_C2_a2_coefficient='-1')
        allsub={f'rho{i}':str(d[i]+sum(T[i,j]*syms[j]for j in range(f)))for i in range(10)}
        allsub.update({f'B{i}':str(d[10+i]+sum(T[10+i,j]*syms[j]for j in range(f)))for i in range(8)})
        payload=dict(order2_affine_origin=list(map(str,d)),order2_free_original_columns=free_columns,
            order2_affine_embedding=mat(T),free_parameters=names,full_order2_solution=allsub,
            order3_low_monomials=[list(v)for v in labels],order3_raw_quotient_coefficient_matrix=mat(R3),
            order3_coker_coefficient_matrix_before_substitution=mat(R3c),free_monomials=[list(v)for v in words],
            all150_order3_coker_after_full_higher_variable_elimination=mat(residual),independent_coker_rows=eqrows,
            exact_independent_coker_polynomials=list(map(str,polynomials)),tiny_QQ_GB=list(map(str,gb)),tiny_QQ_status=status,
            ordinary_order3_necessary_jet_nonempty=not(any(p==1 for p in gb)),**extra)
        with (out/'certificate.json.gz').open('wb')as file:
            with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as z:z.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
        result=dict(id=STEM,all_passed=True,ordinary_order2_consistent=True,order2_rank=len(columns),order2_free_variables=f,
            all9_flag_directions_and_CA20_free10_retained=True,all150_membership_rows_retained=True,
            order3_full_higher_P2_and_flag8_elimination_uses_same18column_image=True,
            independent_order3_coker_rank=len(eqrows),necessary_order3_equations=list(map(str,polynomials)),
            tiny_QQ_GB=list(map(str,gb)),ordinary_order3_necessary_jet_nonempty=payload['ordinary_order3_necessary_jet_nonempty'],
            actual_formal_or_algebraic_pair_branch_proved=False,all_order_or_nonordinary_pair_exclusion_claimed=False,
            individual_conductor_or_R19_conductor_identification_claimed=False,seconds=time.monotonic()-start)
        if extra:result.update(exact_necessary_coker_ideal='(v2,4*v0*v1+4*v0+7*v1^2+18*v1+8)',
            exact_nonempty_finite3jet_witness_all150_rows_verified=True,
            CA_extra_C2_used_at_order2_in_witness=True,group_logI_plusN_order3_identity_checked=True)
    paths=[Path(__file__),source,accepted,ROOT/'code/ambient_CA_V1_independent_review_20261103_v63.py']
    note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths}
    result['evidence_sha256']={p.name:sha(p)for p in out.iterdir()if p.is_file()}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
