#!/usr/bin/env python3
"""Exact local linear and quadratic data on the two minimal infinity curves."""
import argparse,hashlib,json
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]

def calculate(model,extra,point):
    beta={i:point[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    two,polynomials=specialize(model,beta,extra)
    q0=s.Matrix([point[i] for i in two]);fixed=next(i for i,v in enumerate(q0) if v)
    columns=[i for i in range(45) if i!=fixed];count=len(polynomials)
    derivative=s.zeros(count,45);inhomogeneous=s.zeros(count,1)
    for row,f in enumerate(polynomials):
        value=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q0[i]*q0[j]
                derivative[row,i]+=c*q0[j];derivative[row,j]+=c*q0[i]
            elif len(mon)==1:inhomogeneous[row,0]+=c*q0[mon[0]]
        assert value==0
    jacobian=derivative[:,columns].row_join(inhomogeneous)
    tangent=s.Matrix.hstack(*jacobian.nullspace())
    cokernel=s.Matrix.hstack(*jacobian.T.nullspace()).T
    n=tangent.cols
    linear_forms={i:list(tangent[columns.index(i),:]) if i!=fixed else [0]*n for i in range(45)}
    t=list(tangent[-1,:]);monomials=list(combinations_with_replacement(range(n),2));coefficients=[]
    for polynomial in polynomials:
        out={}
        for mon,c in polynomial.items():
            forms=[linear_forms[i] for i in mon]+[t]*(2-len(mon))
            for i in range(n):
                for j in range(n):
                    key=tuple(sorted((i,j)));out[key]=out.get(key,0)+c*forms[0][i]*forms[1][j]
        coefficients.append([out.get(mon,0) for mon in monomials])
    projected=(cokernel.to_DM()*s.Matrix(coefficients).to_DM()).to_Matrix()
    _,pivots=projected.T.to_DM().rref();selected=list(pivots);basis=projected[selected,:]
    hessians=[]
    for row in basis.tolist():
        f=dict(zip(monomials,row))
        hessians.extend([[f[tuple(sorted((i,j)))]*(2 if i==j else 1) for j in range(n)] for i in range(n)])
    hessian=s.Matrix(hessians);radical=s.Matrix.hstack(*hessian.nullspace())
    assert jacobian*tangent==s.zeros(count,n)
    assert cokernel*jacobian==s.zeros(cokernel.rows,45)
    def rows(matrix):return [[str(v) for v in row] for row in matrix.tolist()]
    return dict(equation_count=count,point61=list(map(str,point)),q_to_u_indices_zero_based=two,
        chart_fixed_u=two[fixed]+1,chart_fixed_value=str(q0[fixed]),chart_columns=columns,
        leading_jacobian_rank=derivative.rank(),homogenized_jacobian_rank=jacobian.rank(),
        jacobian=rows(jacobian),tangent_basis=rows(tangent),homogenizer_zero_on_tangent=tangent[-1,:]==s.zeros(1,n),
        cokernel_rows=rows(cokernel[selected,:]),quadratic_monomials=monomials,
        quadratic_obstruction_basis=rows(basis),quadratic_common_radical=rows(radical),
        quadratic_span_rank=basis.rows,common_hessian_rank=hessian.rank())



def verify():
    source_path=ROOT/'reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json'
    source=json.loads(source_path.read_text());model=S61();out={}
    expected={('M0',134):(41,41,4,6,3,False),('M0',135):(42,42,3,3,2,False),('M1',134):(42,42,3,3,2,False),('M1',135):(42,43,2,1,1,True)}
    for label,record in [('M0',source),('M1',source['M1_family'])]:
        assert record['all_passed']
        assert record['dominance_minor']['rank']==16 and record['beta_line_incidence_minor']['rank']==17
        point=s.Matrix(list(map(s.Rational,record['sample_S61_coordinates'])));out[label]={}
        for extra in [False,True]:
            result=calculate(model,extra,point);n=s.Matrix(result['tangent_basis']).cols
            got=(result['leading_jacobian_rank'],result['homogenized_jacobian_rank'],n,result['quadratic_span_rank'],result['common_hessian_rank'],result['homogenizer_zero_on_tangent'])
            assert got==expected[label,134+extra]
            radical=s.Matrix(result['quadratic_common_radical']);T=s.Matrix(result['tangent_basis'])
            assert radical.cols==1 and T[-1,:]*radical==s.zeros(1,1)
            assert result['quadratic_span_rank']==n*(n-1)//2
            result['local_dimension']=1;result['transverse_length_upper_bound']=n
            out[label][str(134+extra)]=result
    paths=[Path(__file__),source_path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_INFINITY_LOCAL_DATA_20261101',all_passed=True,field='QQ',cases=out,
        proof_note='research/notes/S61_MINIMAL_INFINITY_FORMAL_RINGS_20261101.md',
        limitations=['Exact pointwise local dimension and quadratic data; completed rings require the separate universal rank argument.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();result=verify()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,cases=4,local_dimension=1)))
if __name__=='__main__':main()
