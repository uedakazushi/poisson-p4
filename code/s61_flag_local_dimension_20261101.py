#!/usr/bin/env python3
"""Certify local dimension three at one excess infinity flag over QQ.

Ten quadratic tangent obstructions span Sym² of four normal linear forms.
The included smooth three-dimensional flag incidence supplies the lower bound.
This is the unsaturated 135-generator homogenization, not its true closure.
"""
import argparse,hashlib,json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize,chart,EXTRA
ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    beta={i:s.Integer(k+2) for k,i in enumerate(one)};beta[3]=s.S.Zero
    two,polynomials=specialize(model,beta,True);k,columns,a,b=chart(polynomials,two);jacobian=a.row_join(b)
    tangent=s.Matrix.hstack(*jacobian.nullspace());assert tangent.shape==(45,7) and tangent[44,:]==s.zeros(1,7)
    cokernel=s.Matrix.hstack(*jacobian.T.nullspace()).T;assert cokernel.shape==(97,135)
    rows={j:list(tangent[i,:]) for i,j in enumerate(columns)};rows[k]=[0]*7
    monomials=list(combinations_with_replacement(range(7),2));coefficients=[]
    for polynomial in polynomials:
        out={}
        for mon,c in polynomial.items():
            if len(mon)!=2:continue
            for i in range(7):
                for j in range(7):
                    key=tuple(sorted((i,j)))
                    out[key]=out.get(key,0)+c*rows[mon[0]][i]*rows[mon[1]][j]
        coefficients.append([s.cancel(out.get(mon,0)) for mon in monomials])
    projected=(cokernel.to_DM()*s.Matrix(coefficients).to_DM()).to_Matrix()
    _,pivots=projected.T.to_DM().rref();selected=list(pivots);basis=projected[selected,:];assert basis.rows==10
    hessians=[]
    for row in basis.tolist():
        f=dict(zip(monomials,row))
        hessians.extend([[f[tuple(sorted((i,j)))]*(2 if i==j else 1) for j in range(7)] for i in range(7)])
    hessian=s.Matrix(hessians);assert hessian.rank()==4
    radical=s.Matrix.hstack(*hessian.nullspace());assert radical.shape==(7,3)
    assert jacobian*tangent==s.zeros(135,7) and cokernel*jacobian==s.zeros(97,45)
    report={'id':'S61_FLAG_LOCAL_DIMENSION_20261101','field':'QQ','all_passed':True,
        'beta':{str(i+1):str(v) for i,v in beta.items()},
        'homogenized_chart_variables':'44 weight2 directions except u38, followed by homogenizing t',
        'q_to_u_indices_zero_based':two,'chart_columns':columns,
        'homogenized_jacobian':[[str(v) for v in row] for row in jacobian.tolist()],
        'jacobian_rank':38,'tangent_basis':[[str(v) for v in row] for row in tangent.tolist()],
        'cokernel_rows':[[str(v) for v in row] for row in cokernel[selected,:].tolist()],
        'quadratic_monomials':[list(mon) for mon in monomials],
        'quadratic_obstruction_basis':[[str(v) for v in row] for row in basis.tolist()],
        'quadratic_common_radical':[[str(v) for v in row] for row in radical.tolist()],
        'quadratic_span_rank':10,'common_hessian_rank':4,'local_dimension':3,
        'lower_bound':'The included flag incidence is smooth dimension3 here: varying w=-e4+a*e3 at v=-e2 differentiates its incidence equation by a nonzero multiple of u12=4.',
        'upper_bound':'Ten independent quadrics with common radical3 in a7-dimensional tangent space span Sym² of the four normal forms, so the tangent cone has dimension at most3.',
        'limitations':['No assertion about the t-saturated closure.','No full affine fiber finiteness or degree result.','Quadratic tangent obstructions do not describe the entire completed local ring.'],
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (Path(__file__),ROOT/'code/s61_flag_transverse_obstruction_20261101.py',EXTRA,ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'local_dimension':3,'jacobian_rank':38,'quadratic_obstructions':10,'output':str(args.output)}))

if __name__=='__main__':main()
