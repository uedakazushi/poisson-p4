#!/usr/bin/env python3
"""Exact whole linear kernels for the ten regular quadratic pencil types.

Four distinct eigenvalues and the zero/factor cases have separate uniform
proofs.  These are orbit representatives of complete Segre strata, not
sample points in one family.  No Groebner or jet lifting is used.
"""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import sys
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix


def main():
    started=time.monotonic();signal.alarm(30)
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    xx=s.symbols('x1:5');pairs=tuple(itertools.combinations(range(4),2))
    types=[('single_4',[(0,4)]),('single_31',[(0,3),(0,1)]),
           ('single_22',[(0,2),(0,2)]),('eigen31_block3',[(0,3),(1,1)]),
           ('eigen31_block21',[(0,2),(0,1),(1,1)]),
           ('eigen22_blocks2_2',[(0,2),(1,2)]),
           ('eigen22_blocks2_11',[(0,2),(1,1),(1,1)]),
           ('eigen22_blocks11_11',[(0,1),(0,1),(1,1),(1,1)]),
           ('eigen211_block2',[(0,2),(1,1),(-1,1)]),
           ('eigen211_block11',[(0,1),(0,1),(1,1),(-1,1)])]
    singular='--singular' in sys.argv
    if singular:
        types=[('singular_no_common_kernel',[])]
    def wedge(a,b):
        return s.expand(a[0]*b[5]+a[5]*b[0]-a[1]*b[4]-a[4]*b[1]+a[2]*b[3]+a[3]*b[2])
    results=[]
    for name,blocks in types:
        F=s.zeros(4);G=s.zeros(4);offset=0
        for eigen,size in blocks:
            symmetric=s.zeros(size);jordan=s.eye(size)*eigen
            for i in range(size):symmetric[i,size-1-i]=1
            for i in range(size-1):jordan[i,i+1]=1
            F[offset:offset+size,offset:offset+size]=symmetric
            G[offset:offset+size,offset:offset+size]=symmetric*jordan
            offset+=size
        if singular:
            F[0,1]=F[1,0]=1
            G[0,2]=G[2,0]=1
            G[3,3]=1
        xv=s.Matrix(xx);f=(xv.T*F*xv)[0]/2;g=(xv.T*G*xv)[0]/2
        df=[s.diff(f,x) for x in xx];dg=[s.diff(g,x) for x in xx]
        gamma=[s.expand(df[i]*dg[j]-df[j]*dg[i]) for i,j in pairs]
        assert wedge(gamma,gamma)==0
        cc=s.symbols('b0:24')
        beta=[sum(cc[4*p+i]*x for i,x in enumerate(xx)) for p in range(6)]
        beta_matrix=s.zeros(4)
        for (i,j),value in zip(pairs,beta):beta_matrix[i,j]=value;beta_matrix[j,i]=-value
        equations=[s.diff(beta_matrix[j,k],xx[i])-s.diff(beta_matrix[i,k],xx[j])+s.diff(beta_matrix[i,j],xx[k])
                   for i,j,k in itertools.combinations(range(4),3)]
        equations+=s.Poly(wedge(beta,gamma),*xx).coeffs()
        matrix,_=s.linear_eq_to_matrix(equations,cc)
        domain=DomainMatrix.from_Matrix(matrix).convert_to(s.QQ)
        rank=domain.rank()
        uu=s.symbols('u0:4');vv=s.symbols('v0:4')
        expected=[s.expand(df[i]*vv[j]-df[j]*vv[i]+uu[i]*dg[j]-uu[j]*dg[i]) for i,j in pairs]
        expected_coeff=s.Matrix([s.diff(value,x) for value in expected for x in xx])
        basis=expected_coeff.jacobian(uu+vv)
        assert matrix*basis==s.zeros(matrix.rows,8)
        basis_domain=DomainMatrix.from_Matrix(basis).convert_to(s.QQ)
        basis_rank=basis_domain.rank()
        omega=s.symbols('o0:6')
        omega_equations=s.Poly(wedge(omega,gamma),*xx).coeffs()
        constant_matrix,_=s.linear_eq_to_matrix(omega_equations,omega)
        constant_domain=DomainMatrix.from_Matrix(constant_matrix).convert_to(s.QQ)
        constant_kernel=constant_domain.nullspace().to_Matrix()
        rr=s.symbols('r0:'+str(constant_kernel.rows))
        omega_source=constant_kernel.T*s.Matrix(rr) if rr else s.zeros(6,1)
        pf=s.factor(wedge(list(omega_source),list(omega_source))/2)
        results.append({'name':name,'self_adjoint_Jordan_blocks':blocks,
                        'nondegenerate_quadric':str(s.expand(f)),'second_quadric':str(s.expand(g)),
                        'highest_six_coefficients':list(map(str,gamma)),
                        'full_linear_two_form_coefficient_count':24,
                        'closed_and_orthogonal_matrix':matrix.tolist(),
                        'closed_and_orthogonal_rank':rank,'whole_kernel_dimension':24-rank,
                        'expected_eight_parameter_basis':basis.tolist(),'expected_basis_rank':basis_rank,
                        'expected_basis_annihilated_exactly':True,
                        'constant_orthogonal_matrix':constant_matrix.tolist(),
                        'constant_orthogonal_kernel':constant_kernel.tolist(),
                        'constant_orthogonal_six_form':list(map(str,omega_source)),
                        'constant_orthogonal_Pfaffian':str(pf),
                        'constant_nondegenerate_possible':pf!=0,
                        'whole_linear_kernel_is_expected':24-rank==8 and basis_rank==8})
    here=Path(__file__)
    report={'status':('SINGULAR_LINEAR_KERNEL_PROVED' if singular else 'ALL_TEN_LINEAR_KERNELS_PROVED') if all(x['whole_linear_kernel_is_expected'] for x in results) else 'EXTRA_KERNELS_REMAIN',
            'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),
            'coefficient_field':'QQ','cases':results,
            'scope':('The full symmetric Kronecker3 plus regular1 pencil, f=x1*x2,g=x1*x3+x4^2/2. Completeness of the no-common-kernel singular case and the constant-Pfaffian contradiction are separate geometric dependencies.' if singular else 'All ten factor-free Segre strata with a nondegenerate pencil member and at most three distinct eigenvalues. Geometry of pair changes and the four-distinct, factor, and no-nondegenerate-member complements is separate. Exact whole linear kernels do not themselves assert exhaustion of all highest parts.'),
            'seconds':round(time.monotonic()-started,4)}
    output='reports/CLASSIFICATION_UNIMODULAR_SINGULAR_PENCIL_20261107.json' if singular else 'reports/CLASSIFICATION_UNIMODULAR_QUADRATIC_PENCIL_20261107.json'
    Path(output).write_text(json.dumps(report,indent=2,default=str)+'\n')
    signal.alarm(0)
    print(json.dumps({'status':report['status'],'cases':[{k:r[k] for k in ('name','whole_kernel_dimension','constant_nondegenerate_possible')} for r in results],'seconds':report['seconds']}))


if __name__=='__main__':
    main()
