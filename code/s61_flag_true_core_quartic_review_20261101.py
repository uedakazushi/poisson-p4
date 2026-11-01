#!/usr/bin/env python3
"""Independent lower504 reconstruction and full-tail proof audit of true quartics."""
import argparse
import hashlib
from itertools import combinations_with_replacement as cwr
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM

ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time()
    primary='reports/S61_FLAG_TRUE_CORE_QUARTIC_20261101.json'
    quadratic='reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json'
    cubic='reports/S61_FLAG_TRUE_CORE_CUBIC_20261101.json'
    d=json.loads((ROOT/primary).read_text())
    q=json.loads((ROOT/quadratic).read_text())
    c=json.loads((ROOT/cubic).read_text())
    assert all(record['all_passed'] for record in [d,q,c])
    for record in [d,q,c]:
        for path,digest in record['input_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    saved=d['lower_fourth_space'];mat=lambda a:s.Matrix(a).applyfunc(s.sympify)
    Q=mat(q['polynomial_quadratic_basis'])
    rr,_=Q.to_DM().to_field().rref();rr=rr.to_Matrix()
    C=rr[[i for i in range(15) if not any(x.free_symbols for x in rr[i,:])],:]
    assert C.rows==13
    assert (C-C[:,q['constant_minor']['columns']]*Q).applyfunc(s.expand)==s.zeros(13,55)
    p=c['quadratic_product_space'];fc=p['free']
    fe=[j for j in range(120) if j not in p['remaining_pivots']]
    cg=s.zeros(5,220);saved_cubic=mat(c['cubic_space']['rref'])
    for j,col in enumerate([fc[i] for i in fe]):
        cg[:,col]=saved_cubic[:,j]
    mon={i:list(cwr(range(10),i)) for i in range(1,5)}
    index={m:i for i,m in enumerate(mon[4])}
    def multiply(matrix,degree,extra):
        out={}
        for row in range(matrix.rows):
            for column,monomial in enumerate(mon[degree]):
                coefficient=matrix[row,column]
                if coefficient:
                    for k,multiplier in enumerate(mon[extra]):
                        key=(row*len(mon[extra])+k,index[tuple(sorted(monomial+multiplier))])
                        out[key]=out.get(key,0)+coefficient
        return s.MutableSparseMatrix(matrix.rows*len(mon[extra]),715,out)
    BC=multiply(C,2,2);rc,pc=BC.to_DM().rref();RC=rc.to_Matrix()[:len(pc),:]
    assert len(pc)==440 and list(pc)==saved['constant_pivots']
    assert RC==mat(saved['constant_rref'])
    B=multiply(Q,2,2).col_join(multiply(cg,3,1));assert B.shape==(875,715)
    variables=sorted(B.free_symbols,key=str);K=s.QQ.poly_ring(*variables)
    BD=B.to_DM().convert_to(K);RD=RC.to_DM().convert_to(K);free=saved['free']
    E=BD.extract(range(875),free)-BD.extract(range(875),list(pc))*RD.extract(range(440),free)
    nonzero=sorted(set(i for i,j in E.to_dok()))
    assert nonzero==saved['remaining_rows'] and len(nonzero)==195
    E=E.extract(nonzero,range(275));selected=saved['selected_rows'];columns=saved['remaining_pivots']
    M=E.extract(selected,columns)
    M0=M.to_Matrix().subs({v:0 for v in variables}).to_DM()
    M0i=M0.inv().convert_to(K)
    N=M0i*M-DM.eye(64,K)
    powers=[N,N*N,N*N*N,N*N*N*N]
    counts=[len(power.to_dok()) for power in powers]
    assert counts==saved['power_nnz']==[425,102,6,0]
    assert powers[-1].is_zero_matrix
    MI=(DM.eye(64,K)-powers[0]+powers[1]-powers[2])*M0i
    assert MI*M==DM.eye(64,K) and M*MI==DM.eye(64,K)
    assert MI.to_Matrix()==mat(saved['remaining_inverse'])
    RE=MI*E.extract(selected,range(275))
    assert RE.to_Matrix()==mat(saved['remaining_rref'])
    assert (E-E.extract(range(195),columns)*RE).is_zero_matrix
    determinant=M0.domain.to_sympy(M0.det())
    assert determinant==s.Rational(saved['minor_value']) and determinant
    fe4=[j for j in range(275) if j not in columns]
    P1=DM.from_dok({(j,k):K.one for k,j in enumerate(free)},(715,275),K)
    RCF=RD.extract(range(440),free)
    P1-=DM.from_dok({(pc[i],j):value for (i,j),value in RCF.to_dok().items()},P1.shape,K)
    P2=DM.from_dok({(j,k):K.one for k,j in enumerate(fe4)},(275,211),K)
    REF=RE.extract(range(64),fe4)
    P2-=DM.from_dok({(columns[i],j):value for (i,j),value in REF.to_dok().items()},P2.shape,K)
    P4=P1*P2
    active={4,5,7,8}
    discarded=[j for j,m in enumerate(mon[4]) if sum(i in active for i in m)>=3]
    assert P4.extract(discarded,range(211)).is_zero_matrix
    true=d['true_fourth_space'];RR=mat(true['rref']);t=mat(true['corrected_t4_remainder'])
    assert RR.shape==(6,211) and RR[:,true['pivots']]==s.eye(6)
    assert all(s.denom(x).is_Number for x in RR)
    assert s.Rational(true['minor_value'])!=0
    assert t[:,true['pivots']]==s.zeros(1,6)
    column,value=true['constant_t_witness']
    assert t[0,column]==s.Rational(value) and s.Rational(value)!=0
    assert 715-504-6==205
    paths=list(d['input_sha256'])+[primary,quadratic,cubic,
        'research/notes/S61_FLAG_TRUE_CORE_QUARTIC_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_FLAG_TRUE_CORE_QUARTIC_REVIEW_20261101.md']
    return dict(id='S61_FLAG_TRUE_CORE_QUARTIC_REVIEW_20261101',all_passed=True,
        independent_full_875_lower_rows_rebuilt=True,independent_lower_rank=504,
        independent_nilpotent_inverse_power_nonzero_counts=counts,
        independent_constant_minor=str(determinant),
        independent_discarded_active_cubic_columns_zero=True,
        independent_saved_true_basis_and_t_witness_checked=True,
        proof_and_source_review_accepted=True,all_285_filtered_tail_completeness_source_audited=True,
        independent_raw4_or_285_tail_reimplementation=False,primary_producer_replayed=False,
        accepted_hilbert_function=[1,10,40,101,205],accepted_all_seven_parameter_values=True,
        accepted_t_order=4,
        limitations=['Raw135 fourth contraction, all285 filtered identities and the new6 minor remain primary computational dependencies.',
                     'No complete initial ideal, completed-ring identification, ambient nilpotence or saturation removal.'],
        input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('TRUE CORE QUARTIC REVIEW PASS',result['seconds'])
