#!/usr/bin/env python3
"""Determine the completed local rings at the previously certified infinity flag.

The new input is a universal leading-Jacobian rank bound along flag incidence.
Combined with the smooth flag threefold and quadratic length bound, finite
module arguments remove t and determine a square-zero thickening.
"""
import argparse,hashlib,json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize,chart,EXTRA
from weight1_structure import minor_witness
ROOT=Path(__file__).resolve().parents[1]

def rows(matrix):return [[str(v) for v in row] for row in matrix.tolist()]

def generic_ranks(model):
    u=s.symbols('u1:62');beta={i:u[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1};beta[3]=0
    two,polynomials=specialize(model,beta,True);k=two.index(37)
    derivative=s.Matrix([[f.get(tuple(sorted((k,j))),0)*(2 if j==k else 1) for j in range(45)] for f in polynomials])
    assert all(not f.get((k,k),0) for f in polynomials)
    constant=derivative[:70,:];kernel=s.Matrix.hstack(*constant.nullspace());assert constant.rank()==15
    residual=derivative[70:,:]*kernel;out=[]
    for extra in (False,True):
        matrix=residual if extra else residual[:64,:]
        reduced,pivots=matrix.to_DM().rref();reduced=reduced.to_Matrix()
        free=[j for j in range(30) if j not in pivots];relation=s.zeros(30,len(free))
        for i,j in enumerate(free):
            relation[j,i]=1
            for row,col in enumerate(pivots):relation[col,i]=-reduced[row,j]
        assert all(s.cancel(f)==0 for f in matrix*relation)
        assert relation[free,:]==s.eye(len(free))
        variables=sorted(matrix.free_symbols,key=str)
        sample={v:s.Integer((-1)**i*(i*i+3*i+5)) for i,v in enumerate(variables)}
        minor=minor_witness(matrix.subs(sample));expected=22 if extra else 21
        assert minor['rank']==len(pivots)==expected
        out.append(dict(equation_count=134+int(extra),constant_rank=15,constant_kernel=rows(kernel),
            residual=rows(matrix),rational_kernel=rows(relation),sample={str(k):str(v) for k,v in sample.items()},
            residual_minor=minor,generic_leading_jacobian_rank=15+expected,
            universal_rank_upper_bound=True))
    return out

def local134(model):
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    beta={i:s.Integer(k+2) for k,i in enumerate(one)};beta[3]=0
    two,polynomials=specialize(model,beta,False);k,columns,a,b=chart(polynomials,two);jacobian=a.row_join(b)
    tangent=s.Matrix.hstack(*jacobian.nullspace());cokernel=s.Matrix.hstack(*jacobian.T.nullspace()).T
    assert tangent.shape==(45,8) and tangent[-1,:]==s.zeros(1,8)
    forms={j:list(tangent[i,:]) for i,j in enumerate(columns)};forms[k]=[0]*8
    monomials=list(combinations_with_replacement(range(8),2));coefficients=[]
    for f in polynomials:
        out={}
        for mon,c in f.items():
            if len(mon)!=2:continue
            for i in range(8):
                for j in range(8):
                    key=tuple(sorted((i,j)));out[key]=out.get(key,0)+c*forms[mon[0]][i]*forms[mon[1]][j]
        coefficients.append([out.get(mon,0) for mon in monomials])
    projected=(cokernel.to_DM()*s.Matrix(coefficients).to_DM()).to_Matrix();_,pivots=projected.T.to_DM().rref()
    selected=list(pivots);basis=projected[selected,:];hessians=[]
    for row in basis.tolist():
        f=dict(zip(monomials,row));hessians.extend([[f[tuple(sorted((i,j)))]*(2 if i==j else 1) for j in range(8)] for i in range(8)])
    hessian=s.Matrix(hessians);radical=s.Matrix.hstack(*hessian.nullspace())
    assert a.rank()==36 and jacobian.rank()==37 and basis.rows==15 and hessian.rank()==5 and radical.cols==3
    assert jacobian*tangent==s.zeros(134,8) and cokernel*jacobian==s.zeros(cokernel.rows,45)
    return dict(beta={str(i+1):str(v) for i,v in beta.items()},jacobian=rows(jacobian),tangent=rows(tangent),
        cokernel_rows=rows(cokernel[selected,:]),monomials=monomials,quadratic_obstructions=rows(basis),
        common_radical=rows(radical),jacobian_rank=37,tangent_dimension=8,quadratic_span_dimension=15,
        normal_dimension=5,smooth_flag_dimension=3)

def verify():
    model=S61();generic=generic_ranks(model);primary=local134(model)
    old=ROOT/'reports/S61_FLAG_LOCAL_DIMENSION_20261101.json';augmented=json.loads(old.read_text())
    assert augmented['all_passed'] and augmented['jacobian_rank']==38
    assert augmented['quadratic_span_rank']==10 and augmented['common_hessian_rank']==4
    inputs=[Path(__file__),ROOT/'code/s61_flag_transverse_obstruction_20261101.py',EXTRA,old,
        ROOT/'research/notes/S61_NONBINARY_FORMAL_RING_20261101.md',
        ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md']
    return dict(id='S61_FLAG_FORMAL_RING_20261101',all_passed=True,field='QQ',generic_ranks=generic,local134=primary,
        local135_source=str(old.relative_to(ROOT)),
        complete_ring134='C[[u1,u2,u3,e1,e2,e3,e4,e5]]/(e1,e2,e3,e4,e5)^2',
        complete_ring135='C[[u1,u2,u3,e1,e2,e3,e4]]/(e1,e2,e3,e4)^2',
        homogenizer_zero=True,proof_note='research/notes/S61_FLAG_FORMAL_RING_20261101.md',
        limitations=['The specified flag and the open incidence family with the same rank conditions.','No assertion that every infinity flag or every infinity stratum is removed.','No full affine-fiber finiteness or exhaustion.'],
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    data=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(data,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,generic_leading_ranks=[36,37],squarezero_ranks=[5,4],homogenizer_zero=True)))
if __name__=='__main__':main()
