#!/usr/bin/env python3
"""Three full compatible-quadratic components with square-zero normal plane.

Generic ranks are computed over QQ(p3,p4,p5); obstruction ranks are checked
at a point where all vector bundles have their generic ranks. Thus the
certificate supplies a family argument, not merely sample tangent ranks.
"""
from v7_paths import ROOT, DATA, REFERENCE
from pathlib import Path
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix
from independent_minimal_thickness_audit import Calculus,rank,constant_compatibility,symplectic_orbit
HERE=Path(__file__).resolve().parent

def field_rank(M,field):
    rr,piv=DomainMatrix.from_Matrix(M).convert_to(field).rref()
    return len(piv),list(piv)

def core(C,name,p3,p4,p5):
    x,y,z,w=C.x;p=p3*x*x+p4*x*z+p5*z*z
    if name=='Cminus':
        q=[x*x/16+x*z/2-z*z/16,0,3*x*x/16-x*z/4+z*z/16,
           -x*x/16+x*z/4-3*z*z/16,p+(x+z)*(y-w)/8,-x*x/16+x*z/2+z*z/16]
    elif name=='U1':
        q=[-z*z,0,-x*x/2-4*x*z-4*z*z,-z*z/2,p+x*y+2*y*z,z*z+x*z]
    elif name=='AB':
        q=[x*z/2,0,z*z+x*x/4,-x*x/16-z*z/4,p,x*z/2]
    else:raise ValueError(name)
    return C.tensor(dict(zip(C.pairs,q)))

def main():
    C=Calculus(4);C.x=s.symbols('x y z w');x,y,z,w=C.x
    # Allow polynomial parameters; point computations still use exact QQ.
    C.pol=lambda f:{m:a for m,a in s.Poly(f,*C.x).terms() if a}
    p3,p4,p5=s.symbols('p3 p4 p5');params=[p3,p4,p5];field=s.QQ.frac_field(*params)
    J=C.tensor({(0,1):1,(2,3):1});CC,B=constant_compatibility(C,J)
    assert rank(CC)==15 and rank(B)==45
    report={}
    for name in ['Cminus','U1','AB']:
        F=core(C,name,*params);FM=C.as_matrix(F)
        assert not any(C.cubic_vector(C.jacobi(F)))
        assert CC*C.vector(F)==s.zeros(16,1)
        assert s.expand(FM[0,1]+FM[2,3]-x*z)==0
        DG=C.derivative(F)*B
        gr,gp=field_rank(DG,field);assert gr==33
        OG=symplectic_orbit(C,J,F)
        family=[C.tensor({I:s.diff(FM[I],p) for I in C.pairs}) for p in params]
        AG=s.Matrix.hstack(*[C.vector(Z) for Z in OG+family])
        agr,agp=field_rank(AG,field);assert agr==10
        at={p3:1,p4:2,p5:3};Q=core(C,name,1,2,3);D=DG.subs(at)
        A=AG.subs(at);assert rank(D)==33 and rank(A)==10
        O=s.Matrix.hstack(*[C.vector(Z) for Z in symplectic_orbit(C,J,Q)])
        orbit_rank=rank(O);orbit_scale_rank=rank(O.row_join(C.vector(Q)))
        K=s.Matrix.hstack(*D.nullspace());KV=B*K
        piv=A.row_join(KV).rref()[1]
        normal=[C.from_vector(KV[:,i-A.cols]) for i in piv if i>=A.cols]
        assert len(normal)==2
        R,S=normal
        obs=s.Matrix.hstack(C.cubic_vector(C.jacobi(R)),C.cubic_vector(C.polar(R,S)),C.cubic_vector(C.jacobi(S)))
        assert rank(D.row_join(obs))==36
        assert not any(C.cubic_vector(C.polar(Q,R))) and not any(C.cubic_vector(C.polar(Q,S)))
        assert CC*C.vector(R)==s.zeros(16,1) and CC*C.vector(S)==s.zeros(16,1)
        # An explicit residual3-by3 determinant certifies all normal quadrics.
        dp=D.rref()[1];D0=D[:,list(dp)];rows=list(D0.T.rref()[1])
        residual=obs-D0*D0[rows,:].inv()*obs[rows,:]
        rr=list(residual.T.rref()[1]);minor=residual[rr,:]
        assert len(rr)==3 and minor.det()!=0
        entry={'generic_Jacobi_rank':gr,'generic_actual_family_rank':agr,
               'generic_rank_field':'QQ(p3,p4,p5)','generic_Jacobi_pivots':gp,'generic_family_pivots':agp,
               'sample':{'p3':1,'p4':2,'p5':3},'sample_Jacobi_rank':33,'sample_actual_rank':10,
               'Sp4_orbit_rank':orbit_rank,'orbit_times_scale_rank':orbit_scale_rank,
               'normal_dimension':2,'normal_quadratic_obstruction_rank':3,
               'normal_quadratic_residual_rows':rr,'normal_quadratic_residual_minor':[[str(v) for v in minor.row(i)] for i in range(3)],
               'normal_quadratic_residual_determinant':str(s.factor(minor.det())),
               'Q':C.encode(Q),'R':C.encode(R),'S':C.encode(S),
               'component_dimension':10,'generic_transverse_ring':'K[u,v]/(u,v)^2','generic_length':3,
               'all_passed':True}
        report[name]=entry
        print(name,{k:entry[k] for k in ['generic_Jacobi_rank','generic_actual_family_rank','Sp4_orbit_rank','orbit_times_scale_rank','normal_quadratic_residual_determinant']},flush=True)
    report['all_passed']=True
    (DATA/'closed_squarezero_certificate.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
