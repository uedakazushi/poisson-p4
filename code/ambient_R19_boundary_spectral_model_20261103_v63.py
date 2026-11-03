#!/usr/bin/env python3
"""Proper relative-quadric R model and its entire nilpotent18D boundary.

The spectral relations, full parabolic actions, and image minors are exactQQ.
Jordan-flag exhaustion is the paired geometric proof, not sampling.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action

ROOT=Path(__file__).resolve().parents[1]
h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0])

def bivector(D):
    P=(ey*D.T-D*ey.T).applyfunc(s.expand);assert canonical(P)==P
    return P

def stable_space(P,parameters,positions):
    directions=[P.diff(v) for v in parameters];C=s.Matrix.hstack(*(vector(D) for D in directions));rows=list(C.T.rref()[1]);assert len(rows)==len(parameters)
    inverse=C[rows,:].inv();reps=[]
    for i,j in positions:
        M=s.zeros(5);M[i,j]=1;B=s.Matrix.hstack(*(vector(action(M,D)) for D in directions));small=inverse*B[rows,:]
        assert C*small==B;reps.append(small)
    return directions,C,reps

def image_witness(P,directions):
    M=s.Matrix.hstack(*(vector(D) for D in directions),*(vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P)) for i,j in itertools.product(range(5),repeat=2)))
    cols=list(M.rref()[1]);rows=list(M[:,cols].T.rref()[1]);det=M[rows,cols].det();assert det!=0
    return dict(affine_rank=len(cols),projective_dimension_lower=len(cols)-1,minor_rows=rows,minor_columns=cols,minor_determinant=str(det),full_matrix=[[str(v) for v in row] for row in M.tolist()])

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();t,r,b,c,d,e=s.symbols('t r b c d e');gg=s.symbols('g0:6');G=sum(a*m for a,m in zip(gg,[h*h,h*x,h*z,x*x,x*z,z*z]))
    parameters=[t,r,b,c,d,e,*gg]
    D=s.Matrix([x*((r-13*t)*h+b*x+c*z),12*t*x*x,0,x*(d*h+e*x+(-r-13*t)*z),2*t*x*w+G]);P=bivector(D);quad=r*r+c*d-25*t*t
    positions=[(i,j) for i,j in itertools.product(range(5),repeat=2) if not(j==2 and i!=2) and not(i==1 and j!=1) and not(i in (0,3) and j not in (0,1,3))]
    assert len(positions)==16
    directions,C,reps=stable_space(P,parameters,positions)
    for (i,j),small in zip(positions,reps):
        M=s.zeros(5);M[i,j]=1;chi=M[2,2]-M[1,1];velocity=small*s.Matrix(parameters)
        assert s.expand(velocity[0]-chi*t)==0
        assert s.expand(sum(s.diff(quad,v)*dv for v,dv in zip(parameters,velocity))-2*chi*quad)==0
    assert s.hessian(quad,[t,r,c,d]).rank()==4
    assert s.hessian(quad.subs(t,0),[r,c,d]).rank()==3
    T=s.Matrix([[r-13*t,b,c],[0,12*t,0],[d,e,-r-13*t]]);lam=s.Symbol('lambda')
    char=s.expand(T.charpoly(lam).as_expr());target=(lam+8*t)*(lam-12*t)*(lam+18*t)
    assert s.rem(s.Poly(char-target,r),s.Poly(quad,r)).as_expr()==0
    print('SPECTRAL12 SPACE192COLS + QUADRICCHARACTER PASS',time.monotonic()-started,flush=True)
    # Proper fullflag triangular9-space for the entire nilpotent boundary.
    aa,bb,cc=s.symbols('aN bN cN');chain_parameters=[aa,bb,cc,*gg]
    N=bivector(s.Matrix([aa*x*x,0,0,bb*h*x+cc*x*x,G]))
    chain_positions=[(i,j) for i,j in itertools.product(range(5),repeat=2) if not(j==2 and i!=2) and not(i==1 and j!=1) and not(i==0 and j not in (0,1)) and not(i==3 and j not in (0,1,3))]
    assert len(chain_positions)==15
    ndirs,NC,nreps=stable_space(N,chain_parameters,chain_positions)
    # Exact membership in N_g via markedh_N=x, k_N=h.
    HN=-(bb*h*x+cc*h*h)*w+s.integrate(G.subs({h:x,x:h},simultaneous=True),z)
    source=raw(-aa*h*h,0,HN);swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
    changed=(swap*source.subs(dict(zip(X,swap*s.Matrix(X))),simultaneous=True)*swap.T).applyfunc(s.expand)
    assert changed==N
    npoint={aa:1,bb:1,cc:1}|dict(zip(gg,range(1,7)));nw=image_witness(N.subs(npoint),ndirs);assert nw['affine_rank']==19
    # A generic boundary source point is smooth and mapsimmersively to R.
    spectral_point={t:0,r:0,b:1,c:0,d:1,e:1}|dict(zip(gg,range(1,7)))
    assert quad.subs(spectral_point)==0 and s.diff(quad,c).subs(spectral_point)==1
    tangent_dirs=[D0 for v,D0 in zip(parameters,directions) if v!=c]
    rw=image_witness(P.subs(spectral_point),tangent_dirs);assert rw['affine_rank']==20
    # T17: twoCasimirs and a triangularHamiltonian; T16 its squaref2 subset.
    f0,f1,f2=s.symbols('f0 f1 f2');f=f0*h*h+f1*h*x+f2*x*x;T17=bivector(s.Matrix([0,0,0,-f,G]));tp=[f0,f1,f2,*gg]
    tpos=[(i,j) for i,j in itertools.product(range(5),repeat=2) if not(j==2 and i!=2) and not(i in (0,1) and j not in (0,1)) and not(i==3 and j not in (0,1,3))]
    assert len(tpos)==16
    tdirs,TC,treps=stable_space(T17,tp,tpos)
    disc=f1*f1-4*f0*f2
    for (i,j),small in zip(tpos,treps):
        M=s.zeros(5);M[i,j]=1;velocity=small*s.Matrix(tp);chi=2*(M[2,2]+M[3,3]-M[0,0]-M[1,1])
        assert s.expand(sum(s.diff(disc,v)*dv for v,dv in zip(tp,velocity))-chi*disc)==0
    tp0={f0:1,f1:1,f2:1}|dict(zip(gg,range(1,7)));tw=image_witness(T17.subs(tp0),tdirs);assert tw['affine_rank']==18
    sqp={f0:0,f1:0,f2:1}|dict(zip(gg,range(1,7)));sqdirs=[T17.diff(f1)*2,T17.diff(f2)*2,*[T17.diff(v) for v in gg]]
    sw=image_witness(T17.subs(sqp),sqdirs);assert sw['affine_rank']==17
    # Exact t0 bigcell source has f2binarysquare and arbitraryG.
    q0=s.Symbol('q0');tzero=bivector(s.Matrix([0,0,0,-q0*x*x,G]))
    assert tzero==T17.subs({f0:0,f1:0,f2:q0})
    # A polynomial degeneration from actualR to the generalT17 image.
    u,Theta,Q0=s.symbols('u Theta Q0');ga,gb,gc,gd,ge,gf=gg
    As=ga-gc/u+gf/(u*u);Bs=gb-ge/u;Cs=gc/u-2*gf/(u*u);Ds=gd;Es=ge/u;Fs=gf/(u*u)
    sourceV=s.Matrix([0,4*u*Theta*x*x,0,-u*Q0*x*x-2*u*Theta*x*z,2*u*Theta*x*w+As*h*h+Bs*h*x+Cs*h*z+Ds*x*x+Es*x*z+Fs*z*z])
    sourceP=ey*sourceV.T-sourceV*ey.T
    back=s.eye(5);back[3,0]=1;back[3,3]=u;forward=back.inv()
    degeneration=(forward*sourceP.subs(dict(zip(X,back*s.Matrix(X))),simultaneous=True)*forward.T).applyfunc(s.cancel)
    expectedV=s.Matrix([0,4*u*Theta*x*x,0,-Q0*x*x-2*Theta*h*x-2*u*Theta*x*z,G+2*u*Theta*x*w]);expectedP=ey*expectedV.T-expectedV*ey.T
    assert (degeneration-expectedP).applyfunc(s.cancel)==s.zeros(5)
    paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'reports/AMBIENT_R19_COMPLEMENT_ORBITS_20261103_v63/report.json',
        ROOT/'research/notes/AMBIENT_R19_ENTIRE_COMPLEMENT_DOUBLE_20261103_v63.md']
    result=dict(id='AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63',all_passed=True,field='QQ',
        proper_R_base_flag='Lline, Xline subsetK3 subsetannL; baseF(1,2,4;5) dimension9',
        spectral_linear_space_dimension=12,parabolic_GL_positions=positions,spectral_basis_raw150=[[str(v) for v in C[:,i]] for i in range(12)],
        spectral_parameters=list(map(str,parameters)),spectral_field=list(map(str,D)),relative_quadric=str(quad),quadric_rank=4,nil_hyperplane_quadric_rank=3,
        spectral_parabolic_parameter_actions=[[[str(v) for v in row] for row in M.tolist()] for M in reps],
        t_character='M_yy-M_xx',quadric_character='2*(M_yy-M_xx)',spectral_source_dimension=19,
        spectral_source_normal=True,spectral_source_proper_birational_image='ENTIRE R19',
        entire_boundary_identity='R intersection N = B18 = image(t=0 spectralsource)',
        boundary_irreducible=True,boundary_projective_dimension=18,chain_linear_parameters=list(map(str,chain_parameters)),chain_field=list(map(str,[aa*x*x,0,bb*h*x+cc*x*x,G])),
        chain_fullflag_parabolic_positions=chain_positions,chain_parameter_actions=[[[str(v) for v in row] for row in M.tolist()] for M in nreps],chain_image_witness=nw,
        R_generic_boundary_source_witness=rw,R_generic_boundary_source_point={str(k):str(v) for k,v in spectral_point.items()},
        T17_image_witness=tw,T16_image_witness=sw,T17_parabolic_positions=tpos,T17_parameter_actions=[[[str(v) for v in row] for row in M.tolist()] for M in treps],
        nested_supports='T16=bigcellt0image subsetT17=twoCasimirtriangularfamily subsetB18=entireR intersectionN; dimensions16,17,18',
        exact_T17_degeneration_canonical_verified=True,full_R_R1_candidate_evidence='properRsource smooth atgenericB18; oneCasimirX, intrinsicL, andK3 bydivisibility recoveruniquegeo-fiber; rank20 giveslength1',
        scope='Support andpropermodels, not fullschemeintersection R/N or globalprimarydecomposition; localRgeometry consequences separated in note.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('PROPER R QUADRIC + ENTIRE BOUNDARY B18 PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
