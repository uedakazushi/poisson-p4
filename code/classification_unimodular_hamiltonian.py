#!/usr/bin/env python3
"""Complete unimodular quadratic sources with a fixed affine Hamiltonian.

The affine Hamiltonian has no zero.  Its six nilpotent pair types are kept.
The three unspecified brackets retain all45 coefficients.  With div P=0,
constant Pfaffian is equivalent to Jacobi; it is linear after fixing P#dl.
"""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix


def main():
    started=time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    x,y,z,w=xx=s.symbols('x y z w')
    pairs=tuple(itertools.combinations(range(4),2))
    monomials=[s.Integer(1),*xx]+[xx[i]*xx[j] for i,j in itertools.combinations_with_replacement(range(4),2)]
    cases=[('nil31_long',s.Matrix([y,z,1,0]),3),
           ('nil31_short',s.Matrix([y,z,0,1]),2),
           ('nil22_top',s.Matrix([y,1,w,0]),3),
           ('nil211_long',s.Matrix([y,1,0,0]),3),
           ('nil211_short_chain_covector',s.Matrix([y,0,1,0]),1),
           ('nil211_short_other_covector',s.Matrix([y,0,1,0]),3)]
    results=[]
    for name,V,ell in cases:
        begun=time.monotonic()
        assert V[ell]==0
        cc=s.symbols('c0:45')
        P=s.zeros(4);cursor=0
        for i,j in pairs:
            if i==ell:
                value=V[j]
            elif j==ell:
                value=-V[i]
            else:
                value=sum(cc[cursor+k]*m for k,m in enumerate(monomials));cursor+=15
            P[i,j]=value;P[j,i]=-value
        assert cursor==45
        div=[s.expand(sum(s.diff(P[i,j],xx[i]) for i in range(4))) for j in range(4)]
        pf=s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
        equations=[]
        for p in div:
            equations += s.Poly(p,*xx).coeffs() if p!=0 else []
        equations += [coefficient for exponent,coefficient in s.Poly(pf,*xx).terms() if any(exponent)]
        matrix,rhs=s.linear_eq_to_matrix(equations,cc)
        augmented=DomainMatrix.from_Matrix(matrix.row_join(rhs)).convert_to(s.QQ)
        reduced,pivots=augmented.rref();R=reduced.to_Matrix()
        if 45 in pivots:
            results.append({'name':name,'affine_Hamiltonian':list(map(str,V)),
                            'linear_Hamiltonian_coordinate':str(xx[ell]),'full_unknown_count':45,
                            'linear_system_shape':list(matrix.shape),'augmented_rank':len(pivots),
                            'inconsistent_augmented_pivot':45,'source_dimension':None,
                            'status':'EMPTY_LINEAR_SOURCE','seconds':round(time.monotonic()-begun,4)})
            continue
        free=[i for i in range(45) if i not in pivots]
        nn=s.symbols('n0:'+str(len(free)))
        solutions=dict(zip((cc[i] for i in free),nn))
        for row,col in enumerate(pivots):
            solutions[cc[col]]=R[row,45]-sum(R[row,j]*solutions[cc[j]] for j in free)
        source=P.subs(solutions,simultaneous=True).applyfunc(s.expand)
        assert all(s.expand(e.subs(solutions,simultaneous=True))==0 for e in equations)
        actual_pf=s.expand(pf.subs(solutions,simultaneous=True))
        assert not any(actual_pf.has(q) for q in xx)
        assert [s.expand(sum(s.diff(source[i,j],xx[i]) for i in range(4))) for j in range(4)]==[0]*4
        assert source[ell,:].T==V
        jacobi=[s.expand(sum(source[i,l]*s.diff(source[j,k],xx[l])
                              +source[j,l]*s.diff(source[k,i],xx[l])
                              +source[k,l]*s.diff(source[i,j],xx[l]) for l in range(4)))
                for i,j,k in itertools.combinations(range(4),3)]
        assert jacobi==[0]*4
        atzero=dict.fromkeys(xx,0)
        linear=s.Matrix([[s.diff(source[i,j],q).subs(atzero) for q in xx] for i,j in pairs])
        minors={','.join(map(str,rows)):s.factor(linear[list(rows),:].det())
                for rows in itertools.combinations(range(6),4)}
        nonzero={key:str(value) for key,value in minors.items() if value!=0}
        status='EMPTY_RANK4_SOURCE' if actual_pf==0 else 'ALL_PRIMITIVE_EVALUATIONS_NONINJECTIVE' if not nonzero else 'INJECTIVE_SOURCE_REMAINS'
        results.append({'name':name,'affine_Hamiltonian':list(map(str,V)),
                        'linear_Hamiltonian_coordinate':str(xx[ell]),'full_unknown_count':45,
                        'linear_system_shape':list(matrix.shape),'linear_system_rank':len(pivots),
                        'source_dimension':len(free),'free_column_indices':free,
                        'full_affine_solution':{str(c):str(solutions[c]) for c in cc},
                        'full_six_brackets':{str(i)+','+str(j):str(source[i,j]) for i,j in pairs},
                        'constant_Pfaffian':str(actual_pf),'all_universal_Jacobi_zero':True,
                        'all_divergence_zero':True,'exact_Hamiltonian_identity':True,
                        'linear_evaluation_matrix':[[str(t) for t in row] for row in linear.tolist()],
                        'all_fifteen_evaluation_minors_zero':not nonzero,'nonzero_evaluation_minors':nonzero,
                        'status':status,'seconds':round(time.monotonic()-begun,4)})
    here=Path(__file__)
    report={'status':'ALL_SIX_TYPES_PROPER' if all(q['status']!='INJECTIVE_SOURCE_REMAINS' for q in results) else 'FULL_SOURCES_WITH_REMAINDER',
            'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),
            'coefficient_field':'QQ','all_six_nilpotent_affine_Hamiltonian_pairs_retained':True,
            'cases':results,'scope':'All unimodular quadratic P with the displayed nowhere-zero affine Hamiltonian and linear Hamiltonian function. Exhaustion of pair types and the semisimple/constant complement are separate elementary geometry dependencies. No sampled point is used.',
            'seconds':round(time.monotonic()-started,4)}
    Path('reports/CLASSIFICATION_UNIMODULAR_AFFINE_HAMILTONIAN_20261107.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({'status':report['status'],'cases':[{q:r[q] for q in ('name','source_dimension','status')} for r in results],'seconds':report['seconds']}))


if __name__=='__main__':
    main()
