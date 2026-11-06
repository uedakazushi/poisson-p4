#!/usr/bin/env python3
"""Whole linear extension sources of the five degenerate L112 quadric types.

The generic residue field is QQ(u). All exceptional residues are recorded
by explicit nonzero pivot minors; no specialization proves the generic case.
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
    started=time.monotonic(); signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    x,y,z,w=s.symbols('x y z w'); xx=(x,y,z,w); u=s.symbols('u')
    pairs=tuple(itertools.combinations(range(4),2)); field=s.QQ.frac_field(u)
    types=[('rank1_aligned_zero',z*z/2+x*w),
           ('rank1_aligned_nonzero',z*z/2+x*w+y*y),
           ('rank1_transverse_zero',z*z/2+(x+y)*w),
           ('rank1_transverse_nonzero',z*z/2+(x+y)*w+x*x),
           ('rank0_full_mixed',x*z+y*w)]
    def dm(m): return DomainMatrix.from_Matrix(m).convert_to(field)
    def pf(p): return s.expand(p[0]*p[5]-p[1]*p[4]+p[2]*p[3])
    def clean(p): return s.factor(s.cancel(p))
    boundary='--boundary' in sys.argv
    double_lines='--double-two-lines' in sys.argv
    complete='--complete' in sys.argv or boundary or double_lines
    if boundary:
        types=[('x3y_aligned_zero',z*z/2+x*w),
               ('x3y_aligned_nonzero',z*z/2+x*w+y*y),
               ('x3y_transverse',z*z/2+(x+y)*w),
               ('x3y_full_mixed',x*z+y*w),
               ('x2y2_aligned_zero',z*z/2+x*w),
               ('x2y2_aligned_nonzero',z*z/2+x*w+y*y),
               ('x2y2_transverse',z*z/2+(x+y)*w+u*x*x),
               ('x2y2_full_mixed',x*z+y*w)]
    if double_lines: types=[('x2yz_b_nonzero',s.Integer(0)),('x2yz_b_zero',s.Integer(0))]
    records=[]
    queue=[(name,q,None) for name,q in types]
    def null_basis(matrix):
        rref,pivots=dm(matrix).rref()
        rref=rref.to_Matrix()
        free=[j for j in range(matrix.cols) if j not in pivots]
        result=s.zeros(matrix.cols,len(free))
        for j,k in enumerate(free):
            result[k,j]=1
            for i,p in enumerate(pivots): result[p,j]=-rref[i,k]
        return result
    for name,q,value in queue:
        field=s.QQ.frac_field(u) if value is None else s.QQ
        residue=s.Integer(3) if boundary and name.startswith('x3y') else s.Integer(1) if boundary else u
        dq=[s.diff(q,v) for v in xx]
        gamma=[]
        for i,j in pairs:
            gamma.append(s.expand((residue-1)*q*int(i==0 and j==1)
                         +residue*y*(int(i==0)*dq[j]-int(j==0)*dq[i])
                         +x*(int(i==1)*dq[j]-int(j==1)*dq[i])))
        if boundary: gamma[0]+=x*x if name.startswith('x3y') else x*y
        Q=s.Matrix([gamma[5],-gamma[4],gamma[3],gamma[2],-gamma[1],gamma[0]])
        if double_lines:
            aa,bb=(u,s.Integer(1)) if name=='x2yz_b_nonzero' else (s.Integer(1),s.Integer(0))
            Q=s.Matrix([x*y,-x*z,bb*x*x,2*y*z,-y*w-(2*bb+aa)*x*y,z*w+aa*x*z])
        if value is not None: Q=Q.subs(u,value)
        assert pf(Q)==0
        cc=s.symbols('l0:24'); L=s.Matrix([sum(cc[4*i+j]*xx[j] for j in range(4)) for i in range(6)])
        Pmat=s.zeros(4)
        for (i,j),v in zip(pairs,L): Pmat[i,j]=v;Pmat[j,i]=-v
        equations=[sum(s.diff(Pmat[i,j],xx[i]) for i in range(4)) for j in range(4)]
        equations+=s.Poly(pf(Q+L)-pf(Q)-pf(L),*xx).coeffs()
        M,_=s.linear_eq_to_matrix(equations,cc)
        _,cols=dm(M).rref(); _,rows=dm(M.T).rref(); rank=len(cols)
        pivot=clean(M.extract(rows,cols).det(method='domain-ge'))
        K=null_basis(M) if complete else dm(M).nullspace().to_Matrix().T
        assert all(clean(v)==0 for v in M*K)
        trans=s.Matrix([s.diff(v,xx[j]) for v in Q for j in range(4)]).jacobian(xx)
        assert all(clean(v)==0 for v in M*trans)
        _,trrows=dm(trans.T).rref(); _,trcols=dm(trans).rref(); trrank=len(trrows)
        trminor=clean(trans.extract(trrows,trcols).det(method='domain-ge'))
        # Coordinate translations preserve the full Pfaffian condition. Fix
        # exactly four source coordinates using a certified translation minor.
        cuts=s.zeros(trrank,24)
        for i,r in enumerate(trrows): cuts[i,r]=1
        sliced=M.col_join(cuts)
        NK=null_basis(sliced) if complete else dm(sliced).nullspace().to_Matrix().T
        normal=s.symbols('n0:'+str(NK.cols))
        coeff=NK*s.Matrix(normal) if normal else s.zeros(24,1)
        Ln=s.Matrix([sum(coeff[4*i+j]*xx[j] for j in range(4)) for i in range(6)])
        constants=s.symbols('c0:6'); C=s.Matrix(constants)
        degree2=[s.Poly(pf(Q+C)-pf(Q)-pf(C),*xx).coeff_monomial(mon) for mon in [xx[i]*xx[j] for i in range(4) for j in range(i,4)]]
        G,_=s.linear_eq_to_matrix(degree2,constants)
        _,gc=dm(G).rref(); _,gr=dm(G.T).rref()
        gminor=clean(G.extract(gr,gc).det(method='domain-ge'))
        free=[i for i in range(6) if i not in gc]
        cf=s.symbols('a0:'+str(len(free)))
        Csource=s.zeros(6,1)
        for i,a in zip(free,cf): Csource[i]=a
        lpf=s.Poly(pf(Ln),*xx)
        mon2=[xx[i]*xx[j] for i in range(4) for j in range(i,4)]
        rhs=s.Matrix([-lpf.coeff_monomial(mon2[i]) for i in gr])-G.extract(gr,free)*s.Matrix(len(cf),1,list(cf))
        sol=G.extract(gr,gc).inv()*rhs
        for i,v in zip(gc,sol): Csource[i]=clean(v)
        polynomial=s.Poly(pf(Q+Ln+Csource),*xx)
        residual=sorted(set(str(clean(v)) for mon,v in polynomial.terms() if sum(mon)>0 and clean(v)!=0))
        finalpf=clean(polynomial.coeff_monomial(1))
        bad=clean(pivot*trminor*gminor)
        if complete:
            denominators=[s.denom(s.cancel(e)) for e in list(K)+list(NK)+list(Csource)]
            bad=clean(bad*s.prod(set(denominators)))
        if complete and value is None:
            numerator,denominator=s.fraction(bad)
            factors=s.factor_list(numerator*denominator,u)[1]
            exceptional=[]
            for fac,power in factors:
                assert s.degree(fac,u)==1
                root=s.solve(fac,u)[0]
                if boundary or double_lines or root not in (0,1): exceptional.append(root)
            for root in sorted(set(exceptional)): queue.append((name,q,root))
        records.append({'name':name,'residue':str(value) if value is not None else 'generic u','quadric':str(q),'highest_bivector':list(map(str,Q)),
                        'full_linear_matrix':M.tolist(),'rank':rank,'pivot_rows':list(rows),'pivot_columns':list(cols),'pivot_determinant':str(pivot),
                        'full_kernel':K.tolist(),'translation_basis':trans.tolist(),'translation_rank':trrank,'translation_cut_rows':list(trrows),'translation_minor':str(trminor),
                        'normal_dimension':NK.cols,'normal_linear_bivector':list(map(str,Ln)),
                        'constant_equation_matrix':G.tolist(),'constant_pivot_rows':list(gr),'constant_pivot_columns':list(gc),'constant_minor':str(gminor),
                        'constant_source':list(map(str,Csource)),'remaining_full_Pfaffian_equations':residual,'constant_Pfaffian':str(finalpf),
                        'exceptional_parameter_polynomial':str(bad),'already_covered_residues':['zero magnetic coefficient'] if boundary else ['u=0','u=1','t=0']})
    here=Path(__file__)
    report={'status':'DOUBLE_AND_TWO_DISTINCT_LINEAR_FACTORS_ALL_SOURCES' if double_lines else 'TWO_QUARTIC_POLE_TYPES_ALL_DEGENERATE_QUADRICS' if boundary else 'ALL_FIVE_TYPES_AND_EVERY_EXCEPTIONAL_RESIDUE' if complete else 'WHOLE_GENERIC_SOURCES_WITH_ALL_EXCEPTIONAL_MINORS','coefficient_field':'QQ(u), and QQ at every exceptional residue' if complete else 'QQ(u)','cases':records,
            'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),
            'scope':'For D=x^2 yz and independent linear numerator f=w, Q=(X+aN) wedge (Y+bN), X=diag(1,-2,0,1),Y=diag(0,1,-1,0),N=x partial_w. b nonzero normalizes to1, otherwise a nonzero normalizes to1. The zero pair is the solved five-edge log type; dependent linears are the proved pullback type.' if double_lines else 'For D=x^3 y or x^2 y^2, the complete degenerate-restriction quadric normal forms with nonzero magnetic coefficient normalized to one; all parameter-minor zeros are included. Full nondegenerate-restriction proofs and rational partial-fraction completeness are separate geometric dependencies.' if boundary else 'All five canonical quadrics with degenerate restriction to ker(x,y), at residues outside the displayed pivot minors. Exceptional residues, source equations and component-boundary limits are not discarded.',
            'seconds':round(time.monotonic()-started,4)}
    output='reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_DOUBLE_LINES_20261107.json' if double_lines else 'reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_BOUNDARY_20261107.json' if boundary else 'reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_COMPLETE_20261107.json' if complete else 'reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_20261107.json'
    Path(output).write_text(json.dumps(report,indent=2,default=str)+'\n')
    print(json.dumps({'status':report['status'],'seconds':report['seconds'],'cases':[{k:v[k] for k in ('name','residue','rank','normal_dimension','exceptional_parameter_polynomial','remaining_full_Pfaffian_equations','constant_Pfaffian')} for v in records]}))


if __name__=='__main__': main()
