#!/usr/bin/env python3
"""Full projective reduced CA/CB intersection: N27 plus one rank2 R19.

Finite double-flag linear algebra covers the rank2 boundary. One explicit
point excludes the remaining19-dimensional family from the global N closure.
No new full Poisson component or global scheme-theoretic splitting is claimed.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import raw,canonical,vector,X,h,x,y,z,w
from ambient_rank_two_flag_intersection_table_20261103_v63 import main as table,permutation,lie

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    table(output/'double-flag-table')
    data=json.loads((output/'double-flag-table/report.json').read_text())
    assert data['endpoint_cases']==66 and data['entire_N_same_flag_cases']==64
    exceptions=[r for r in data['records']if not r['entire_linear_family_lies_in_N_with_the_CA_flag']]
    assert len(exceptions)==2
    ca_path=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    ca=json.loads(ca_path.read_text());cubics=list(map(s.sympify,ca['cubic_monomials']))
    Aall=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,v)))for v in ca['basis_canonical150']))
    A=Aall[:,[0,1,2,*range(4,20)]]
    def rb(f,H):
        P=raw(f,0,H);P[2,3]-=s.diff(f,x)*z/2;P[3,2]=-P[2,3]
        P[2,4]-=s.diff(f,x)*w/2;P[4,2]=-P[2,4]
        return canonical(P)
    B=s.Matrix.hstack(*(vector(rb(f,0))for f in [h*h,h*x,x*x]),*(vector(rb(0,H))for H in cubics))
    small=[r for r in exceptions if r['relative_permutation']==[3,1,2,0,4]][0]
    perm=small['relative_permutation'];Bp=s.Matrix.hstack(*(permutation(B[:,j],perm)for j in range(19)))
    kernel=(A.row_join(-Bp)).nullspace();KB=s.Matrix.hstack(*(v[19:,0]for v in kernel))
    assert len(kernel)==7 and all(KB[i,j]==0 for i in [0,1]for j in range(7))
    assert any(KB[2,j]for j in range(7))
    # This whole exceptional family has CB binary F proportional to x_B²,
    # hence belongs to N by the separately certified square recognition.
    big=[r for r in exceptions if r['relative_permutation']==[1,3,2,0,4]][0]
    K=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,v)))for v in big['intersection_basis_in_CA_coordinates']))
    L=A*K;assert K.cols==8 and big['full_flag_pair_orbit_dimension']==12
    q=s.symbols('q0:8');parameters=K*s.Matrix(q)
    f=s.expand(sum(parameters[i]*m for i,m in enumerate([h*h,h*x,x*x])))
    H=s.expand(sum(parameters[3+i]*v for i,v in enumerate(cubics)))
    assert f==-4*q[4]*x*x
    expected=3*h*h*q[7]*z+2*h*q[3]*x*z+h*q[6]*z*z+q[0]*w*x*x+q[1]*x*x*z+q[2]*x*z*z/2+2*q[4]*w*x*z+q[5]*z**3/3
    assert s.expand(H-expected)==0
    generic=L*s.Matrix([j%7+1 for j in range(8)])
    D=L.row_join(s.Matrix.hstack(*(lie(generic,a,b)for a,b in itertools.product(range(5),repeat=2))))
    pc=list(D.rref()[1]);pr=list(D[:,pc].T.rref()[1]);assert len(pc)==20
    det=D[pr,pc].det();assert det!=0
    # Explicit rank2 witness: q4=q5=q7=1, all other q=0.
    V=s.Matrix([0,4*x*x,0,-2*x*z,3*h*h+z*z+2*x*w])
    P=s.zeros(5);P[2,:]=V.T;P[:,2]=-V;Q=canonical(P)
    assert vector(Q)==L*s.Matrix([0,0,0,0,1,1,0,1])
    C=s.Matrix([-8*h*x,12*x*x,0,-18*x*z,15*h*h+5*z*z+2*x*w])
    assert Q==(s.eye(5)[:,2]*C.T-C*s.eye(5)[:,2].T)/5
    assert sum(s.diff(C[i],X[i])for i in [0,1,3,4])==0
    # Unique intrinsic constant vector factor: E wedge v wedge Q=0.
    mon3=list(itertools.combinations_with_replacement(range(5),3));factor=[]
    for p in range(5):
        v=[int(i==p)for i in range(5)];row=[]
        for i,j,k,l in itertools.combinations(range(5),4):
            value=s.expand(X[i]*(v[j]*Q[k,l]-v[k]*Q[j,l]+v[l]*Q[j,k])
              -X[j]*(v[i]*Q[k,l]-v[k]*Q[i,l]+v[l]*Q[i,k])
              +X[k]*(v[i]*Q[j,l]-v[j]*Q[i,l]+v[l]*Q[i,j])
              -X[l]*(v[i]*Q[j,k]-v[j]*Q[i,k]+v[k]*Q[i,j]))
            poly=s.Poly(value,*X);row.extend(poly.coeff_monomial(s.prod(X[t]for t in m))for m in mon3)
        factor.append(s.Matrix(row))
    factor=s.Matrix.hstack(*factor);assert factor.rank()==4
    assert factor.nullspace()==[s.eye(5)[:,2]]
    mon2=list(itertools.combinations_with_replacement(range(5),2))
    casimirs=s.Matrix.hstack(*(s.Matrix([s.Poly(Q[j,i],*X).coeff_monomial(s.prod(X[t]for t in m))for i in range(5)for m in mon2])for j in range(5)))
    assert casimirs.rank()==5
    # N_g0 requires C(ell)=lambda ell² and independent k with
    # C(k)=-alpha ell²-4lambda ell*k. lambda0 is excluded by rank4.
    coeff=s.symbols('A B C D');aa,bb,cc,dd=coeff;lam=s.Symbol('lam')
    ell=aa*h+bb*x+cc*z+dd*w;cell=sum(coeff[j]*C[i]for j,i in enumerate([0,1,3,4]))
    poly=s.Poly(s.expand(cell-lam*ell**2),h,x,z,w)
    assert poly.coeff_monomial(w*w)==-lam*dd*dd
    assert poly.coeff_monomial(h*h).subs(dd,0)==-lam*aa*aa
    assert poly.coeff_monomial(z*z).subs(dd,0)==-lam*cc*cc
    assert poly.coeff_monomial(x*x).subs({dd:0,aa:0,cc:0})==12*bb-lam*bb*bb
    # Thus nonzero ell is a multiple of x and lambda*ell=12x.
    kk=aa*h+bb*x+cc*z+dd*w;ck=sum(coeff[j]*C[i]for j,i in enumerate([0,1,3,4]))
    test=s.expand(ck+48*x*kk)
    assert s.expand(test-(40*aa*h*x+60*bb*x*x+30*cc*x*z+dd*(15*h*h+5*z*z+50*x*w)))==0
    # Requiring test proportional to x² forces A=C=D=0, contradicting
    # independence of k from ell. No Darboux classification is assumed.
    input_paths=[Path(__file__),ROOT/'code/ambient_rank_two_flag_intersection_table_20261103_v63.py',
      ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ca_path,
      ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json',
      ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63/report.json',
      ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_INTERSECTION_20261103_v63_REVIEW.json']
    result=dict(id='AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63',all_passed=True,
      full_projective_reduced_intersection_components=['global_N27','rank2_R19'],
      global_reduced_intersection_exhaustion_proved=True,
      new_full_Poisson_component_claimed=False,
      flag_double_cosets=33,endpoint_cases=66,N_with_same_CA_flag_cases=64,
      remaining_N_case_CB_square=True,remaining_R_case_dimension=19,
      R_raw_f=str(f),R_raw_H=str(H),R_constant_factor='d_y',
      R_relative_flag_permutation=[1,3,2,0,4],R_linear_parameter_count=8,
      R_flag_pair_orbit_dimension=12,R_source_projective_dimension_upper=19,
      R_actual_affine_image_rank=20,R_dimension_minor_rows=pr,R_dimension_minor_columns=pc,
      R_dimension_minor_determinant=str(det),
      R_witness_parameters=[0,0,0,0,1,1,0,1],R_witness_C_times_five=list(map(str,C)),
      R_witness_unique_intrinsic_constant_vector_factor=True,
      R_witness_constant_Casimir_covector_rank=5,
      R_witness_excluded_from_N_f_zero=True,R_witness_excluded_from_N_g_zero=True,
      global_scheme_intersection_primary_decomposition_claimed=False,
      local_at_Q_results_unchanged=True,
      input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in input_paths},
      seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('FULL REDUCED GLOBAL INTERSECTION=N27 UNION R19 PASS',result['seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
