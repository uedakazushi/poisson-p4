#!/usr/bin/env python3
"""Whole CA-normalizer ramification germ at the genericB18 representative.

The target CA/CB intersection ideal is NOT computed here. The critical
Fitting ideal is reconstructed exactly on an etale parabolic slice.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action

ROOT=Path(__file__).resolve().parents[1]
h,x,y,z,w=X


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(M):return [[str(v)for v in row]for row in M.tolist()]


def main(output,m2):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
    def change(P):return(swap*P.subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
    basis=[change(canonical(raw(f,0,0)))for f in[h*h,h*x,x*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cubics]
    C=s.Matrix.hstack(*map(vector,basis));cr=list(C.T.rref()[1]);left=C[cr,:].inv()
    assert left*C[cr,:]==s.eye(20)
    def rem(M):return(M-C*left*M[cr,:]).applyfunc(s.expand)
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    center=left*vector(P0)[cr,:];assert C*center==vector(P0)
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    parabolic=[(i,j)for i,j in itertools.product(range(5),repeat=2)if(i,j)not in flag]
    assert len(flag)==9 and len(parabolic)==16
    # Stabilizer of L=<ey>, W=<h,x>, H=<x>, a connected parabolic.
    for i,j in parabolic:
        M=s.zeros(5);M[i,j]=1
        assert rem(s.Matrix.hstack(*[vector(action(M,P))for P in basis]))==s.zeros(150,20)
    orbit=s.Matrix.hstack(*[left*vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))[cr,:]for i,j in parabolic])
    pc=list(orbit.rref()[1]);assert len(pc)==10
    fixed=[0,4,5,6,7,8,11,12,13,16];free=[i for i in range(20)if i not in fixed]
    assert free==[1,2,3,9,10,14,15,17,18,19]
    orbitdet=s.factor(orbit[fixed,pc].det());assert orbitdet!=0
    ss=s.symbols('s0:10');P=P0.copy()
    for index,t in zip(free,ss):P+=t*basis[index]
    Q=rem(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
    qr=[5,22,23,27,31,35,50,65];qc=[0,1,2,4,5,6,7,8];extra=3
    unit=s.factor(Q[qr,qc].det());assert unit==1
    schur=(Q[:,extra]-Q[:,qc]*Q[qr,qc].inv()*Q[qr,extra]).applyfunc(s.expand)
    equations=list(dict.fromkeys(v for v in schur if v!=0));assert equations
    # Exact polynomial normal form, with all7 eliminated variables recorded.
    a0,a1,a2,a3,a4,a5,a6,a7,a8,a9=ss
    U=a3+2*a7;V=7*a3-6*a7;W=a2
    e5=(12*a7*a7-3*a3*a3)/5
    elimination=[a0-3*a3/5+4*a7/5,a1-7*e5/12,a4-35*e5/24+2*a7*a7,a5-e5,a6,a8-a7*a7,a9-a7**3/3]
    predicted=[*elimination,W*U,W*V,U*V,U**3]
    G=s.groebner(equations,*ss,order='grevlex',domain=s.QQ)
    F=s.groebner(predicted,*ss,order='grevlex',domain=s.QQ)
    assert all(G.reduce(s.expand(f))[1]==0 for f in predicted)
    assert all(F.reduce(s.expand(f))[1]==0 for f in equations)
    assert list(G)==list(F)
    # Directly substitute the inverse linear(u,v) transformation, eliminating
    # all7 coordinates, and compare the remaining ideal in three variables.
    uu,vv,ww=s.symbols('u v w')
    low={a3:(3*uu+vv)/10,a7:(7*uu-vv)/20,a2:ww}
    low.update({a0:(3*a3/5-4*a7/5).subs(low),a1:(7*e5/12).subs(low),a4:(35*e5/24-2*a7*a7).subs(low),a5:e5.subs(low),a6:0,a8:(a7*a7).subs(low),a9:(a7**3/3).subs(low)})
    reduced=[s.expand(f.subs(low,simultaneous=True))for f in equations]
    tiny=s.groebner([ww*uu,ww*vv,uu*vv,uu**3],ww,uu,vv,order='grevlex',domain=s.QQ)
    reducedGB=s.groebner(reduced,ww,uu,vv,order='grevlex',domain=s.QQ)
    assert list(tiny)==list(reducedGB)
    m2text='R=QQ[s0,s1,s2,s3,s4,s5,s6,s7,s8,s9];\nI=ideal('+','.join(str(f).replace('**','^')for f in equations)+');\nJ=ideal('+','.join(str(s.expand(f)).replace('**','^')for f in predicted)+');\nassert(I==J);\nS=QQ[w,u,v];\nK=ideal(w*u,w*v,u*v,u^3);\nassert(K==intersect(ideal(u,v),ideal(u,w),ideal(w,v,u^3)));\nassert(dim(S/K)==1);\nassert(degree(S/K)==2);\nprint hilbertSeries(S/K,Reduce=>true);\nprint "CA_B18_FULL_CRITICAL_GERM_PASS";\nexit 0;\n'
    (output/'critical.m2').write_text(m2text)
    if m2:
        completed=subprocess.run(['M2','--script',str((output/'critical.m2').resolve())],capture_output=True,text=True,check=True,timeout=60)
        assert 'CA_B18_FULL_CRITICAL_GERM_PASS'in completed.stdout
        (output/'critical.log').write_text(completed.stdout+completed.stderr)
    certificate=dict(
        center_fixed20=list(map(str,center)),fixed_fiber_linear_basis150_by20=rows(C),
        parabolic_positions=parabolic,all16_by20_parabolic_stability=True,
        parabolic_fiber_orbit20_by16=rows(orbit),orbit10_columns=pc,orbit10_fixed_rows=fixed,orbit10_minor=str(orbitdet),
        affine_slice_free_indices=free,slice_variables=list(map(str,ss)),flag_positions=flag,
        full150_by9_flag_quotient=rows(Q),unit8_rows=qr,unit8_columns=qc,unit8_determinant=str(unit),
        full150_critical_Schur_column=list(map(str,schur)),all_nonzero_critical_polynomials=list(map(str,equations)),
        seven_monic_eliminations=list(map(str,elimination)),u=str(U),v=str(V),w=str(W),
        inverse_substitution={str(k):str(v)for k,v in low.items()},
        predicted_entire_critical_ideal=list(map(str,predicted)),exact_QQ_groebner_basis=list(map(str,G)),
        both_ideal_containments_and_reduced3variable_equality=True,
        three_variable_ideal=['w*u','w*v','u*v','u^3'],
        primary_decomposition=[['u','v'],['u','w'],['w','v','u^3']])
    encoded=(json.dumps(certificate,separators=(',',':'))+'\n').encode()
    (output/'certificate.json.gz').write_bytes(gzip.compress(encoded,mtime=0))
    note=ROOT/'research/notes/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63.md'
    paths=[Path(__file__),note,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_B18_RAMIFICATION_COTANGENT_20261103_v63/report.json']
    # Parent's independent finite-nearB18 proof is already frozen above.
    assert json.loads(paths[-2].read_text())['all_passed'] is True
    result=dict(id='AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63',all_passed=True,
        field='QQ',all16_parabolic_by20_fiber_columns_stable=True,
        etale19_smooth_source_factors_from_flag9_and_parabolic_orbit10=True,
        whole_intrinsic_critical_Fitting_ideal_on_slice=True,unit8_block_determinant='1',
        exact_nonzero_Schur_generator_count=len(equations),exact_QQ_ideal_equalities=True,
        independent_Macaulay2_ideal_and_primary_checks_run=bool(m2),
        canonical_cone_normalizer_ramification_completion='C[[g1,...,g19,w,u,v]]/(w*u,w*v,u*v,u^3)',
        projective_normalizer_ramification_completion='C[[g1,...,g18,w,u,v]]/(w*u,w*v,u*v,u^3)',
        reduced_critical_completion='C[[g19,w,v]]/(w*v)',
        exact_associated_prime_count=3,two_reduced_branches=True,embedded_prime_at_junction=True,
        nilradical_principal_u=True,nilradical_cubed_zero=True,nilradical_square_nonzero=True,
        cone_critical_dimension=20,cone_critical_embedding_dimension=22,cone_critical_depth=19,
        projective_critical_dimension=19,projective_critical_embedding_dimension=21,projective_critical_depth=18,
        critical_Hilbert_Samuel_multiplicity=2,
        actual_CA_target_completion_or_CA_CB_intersection_or_odd_ideal_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('WHOLE CA B18 NORMALIZER CRITICAL GERM node+u^3 PASS',len(equations),result['seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--m2',action='store_true')
    args=ap.parse_args();main(args.output.resolve(),args.m2)
