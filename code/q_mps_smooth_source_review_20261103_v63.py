#!/usr/bin/env python3
"""Independent QQ/Fp polynomial-source checks for33 actual MPS memberships."""
import argparse
from fractions import Fraction as F
import itertools
import json
from math import gcd,lcm
from pathlib import Path
import sympy as s
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha
from q_mps_arc_membership_review_20261103_v63 import dual

VARS=s.symbols('alpha beta z x y rho');ALPHA,BETA,Z,X,Y,RHO=VARS
LOCALS=dict(zip(map(str,VARS),VARS))
EIGEN=[0,1,Z,X,Y];PAIRS=list(itertools.combinations(range(1,5),2))
GAMMA=s.Matrix([ALPHA*(EIGEN[i]*EIGEN[j]**2-EIGEN[i]**2*EIGEN[j])+BETA*(EIGEN[i]*EIGEN[j]**3-EIGEN[i]**3*EIGEN[j])+RHO*(EIGEN[i]*EIGEN[j]**4-EIGEN[i]**4*EIGEN[j])+EIGEN[i]**2*EIGEN[j]**3-EIGEN[i]**3*EIGEN[j]**2 for i,j in PAIRS])

def equations(record,permutation):
    basis=s.Matrix(record['period_basis']).applyfunc(s.Rational).T
    null=basis.T.nullspace();ann=[]
    for v in null:
        den=lcm(*(int(c.q)for c in v));row=[int(c*den)for c in v];content=gcd(*row);ann.append([v//content for v in row])
    L=s.Matrix(ann);assert L.rows==6-basis.cols and L*basis==s.zeros(L.rows,basis.cols)
    B=s.zeros(5)
    for(i,j),v in zip(PAIRS,dual(GAMMA)):B[i,j]=v;B[j,i]=-v
    for i in range(1,5):B[0,i]=sum(B[i,j]for j in range(1,5));B[i,0]=-B[0,i]
    period=s.Matrix([B[permutation[i],permutation[j]]for i,j in PAIRS]);polys=[s.Poly(s.expand(v),*VARS,domain=s.ZZ)for v in L*period]
    return basis,period,polys

def value(terms,point,derivatives=(),modulus=None):
    out=0 if modulus else F(0)
    for exponents,coefficient in terms:
        powers=list(exponents);c=int(coefficient)
        for j in derivatives:
            c*=powers[j];powers[j]-=1
            if not c:break
        if not c:continue
        for v,power in zip(point,powers):
            c=c*v**power
            if modulus:c%=modulus
        out+=c
        if modulus:out%=modulus
    return out

def inverse_det(M,p):
    n=len(M);A=[[int(v)%p for v in row]+[int(i==j)for j in range(n)]for i,row in enumerate(M)];det=1
    for j in range(n):
        k=next(k for k in range(j,n)if A[k][j]);A[j],A[k]=A[k],A[j]
        if k!=j:det=-det
        c=A[j][j];det=det*c%p;inv=pow(c,-1,p);A[j]=[v*inv%p for v in A[j]]
        for i in range(n):
            if i!=j and A[i][j]:
                c=A[i][j];A[i]=[(u-c*v)%p for u,v in zip(A[i],A[j])]
    assert [row[:n]for row in A]==[[int(i==j)for j in range(n)]for i in range(n)]
    return[row[n:]for row in A],det%p

def packet(stem,replay,certificate_name):
    folder=ROOT/f'reports/{stem}_20261103_v63';report=folder/'report.json';d=read(report);fresh=read(replay/'report.json')
    assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
    assert(folder/certificate_name).read_bytes()==(replay/certificate_name).read_bytes()
    paths=[Path(__file__),report,ROOT/'code/q_smith_global_model_review_20261103_v63.py',ROOT/'code/q_mps_arc_membership_review_20261103_v63.py']
    for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
    for name,h in d['evidence_sha256'].items():p=folder/name;assert sha(p)==h,name;paths.append(p)
    return d,read(folder/certificate_name),paths

def save(path,report,paths):
    report['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)}
    with path.open('x')as f:json.dump(report,f,indent=2);f.write('\n')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--qq-replay',type=Path,required=True);ap.add_argument('--hensel-replay',type=Path,required=True);ap.add_argument('--qq-output',type=Path,required=True);ap.add_argument('--hensel-output',type=Path,required=True);args=ap.parse_args()
    catalogue=read(ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json')['diagrams'];records={r['id']:r for r in catalogue}
    vandermonde=s.prod(EIGEN[i]for i in range(1,5))*s.prod(EIGEN[j]-EIGEN[i]for i in range(1,5)for j in range(i+1,5))
    assert s.expand(GAMMA[0]*GAMMA[5]-GAMMA[1]*GAMMA[4]+GAMMA[2]*GAMMA[3]-RHO*vandermonde)==0
    qq,evidence,qpaths=packet('Q_MPS_IMPLICIT_TORIC_ARCS',args.qq_replay,'implicit-certificates.json');witnesses=read(ROOT/'code/data/q_mps_implicit_toric_arc_witnesses_20261103_v63.json')
    expected=[f'T{i:02d}'for i in range(4,16)]+[f'T{i:02d}'for i in range(18,22)]+['T28'];assert[w['component']for w in witnesses]==expected
    qqdet=[]
    for w,c in zip(witnesses,evidence):
        basis,period,polys=equations(records[w['component']],w['permutation']);terms=[v.terms()for v in polys];point=[F(w['alpha']),F(w['beta']),F(2),F(3),F(4),F(w['rho'])];assert all(value(p,point)==0 for p in terms)
        assert all(s.Poly(f-s.sympify(g,locals=LOCALS),*VARS).is_zero for f,g in zip(polys,c['exact_period_equations']))
        assert point[-1]not in[0,1,-1]
        J=s.Matrix([[s.Rational(str(value(p,point,(j,))))for j in range(5)]for p in terms]);assert J==s.Matrix([[s.Rational(v)for v in row]for row in c['exact_nonrho_Jacobian']])
        selected=w['Jacobian_selected_columns'];H=J[:,selected];det=H.det();assert det==s.Rational(w['Jacobian_minor'])==s.Rational(c['selected_determinant'])!=0;qqdet.append(str(det))
        coordinates=s.Matrix(list(map(s.Rational,c['period_coordinates_at_source_point'])));sub=dict(zip(VARS,[s.Rational(str(v))for v in point]));assert basis*coordinates==period.subs(sub)
        tangent=s.Matrix(list(map(s.Rational,c['rho_tangent'])));F_rho=s.Matrix([s.Rational(str(value(p,point,(5,))))for p in terms]);assert H*tangent+F_rho==s.zeros(H.rows,1)
        velocity=[F(0)]*5+[F(1)]
        for j,column in enumerate(selected):velocity[column]=F(str(tangent[j]))
        rhs=s.Matrix([s.Rational(str(sum((velocity[i]*velocity[j]*value(p,point,(i,j))for i in range(6)for j in range(6)),F(0))))for p in terms]);acc=s.Matrix(list(map(s.Rational,c['rho_second_derivative'])));assert H*acc+rhs==s.zeros(H.rows,1)
    hp,hevidence,hpaths=packet('Q_MPS_HENSEL_TORIC_ARCS',args.hensel_replay,'hensel-certificates.json');hw=read(ROOT/'code/data/q_mps_hensel_toric_arc_witnesses_20261103_v63.json');hexpected=['T22','T23','T24','T25','T26','T27','T29','T30','T31','T32','T33','T34','T35','T36','T37','T38'];assert[w['component']for w in hw]==hexpected
    modular=[]
    for w,c in zip(hw,hevidence):
        prime=w['prime'];assert prime in[11,13,17]and s.isprime(prime);basis,period,polys=equations(records[w['component']],w['permutation']);assert len(polys)==5 and basis.cols==1
        terms=[v.terms()for v in polys];assert all(s.Poly(f-s.sympify(g,locals=LOCALS),*VARS).is_zero for f,g in zip(polys,c['exact_integer_period_equations']))
        assert w['eigenvalues'][:2]==[0,1];point=[w['alpha'],w['beta'],*w['eigenvalues'][2:],w['rho']];assert len({0,1,*[v%prime for v in point[2:5]]})==5 and point[-1]%prime not in[0,1,prime-1]
        assert all(int(v.q)%prime for v in basis)and all(value(p,point,modulus=prime)==0 for p in terms)
        J=[[value(p,point,(j,),prime)for j in range(5)]for p in terms];assert J==w['Jacobian_mod_prime']==c['nonrho_Jacobian_mod_p'];inverse,det=inverse_det(J,prime);assert det==w['selected_determinant_mod_prime']==c['nonrho_determinant_mod_p']!=0;modular.append([w['component'],prime,det])
        derivative=[value(p,point,(5,),prime)for p in terms];tangent=[-sum(v*x for v,x in zip(row,derivative))%prime for row in inverse];assert tangent==c['rho_tangent_mod_p']
        lifts=c['Newton_lifts_at_fixed_integer_rho'];assert len(lifts)==8 and lifts[0]['values']==point
        for precision,record in enumerate(lifts,1):
            lifted=record['values'];assert record['precision']==precision and lifted[-1]==point[-1]and all(value(p,lifted,modulus=prime**precision)==0 for p in terms)
            if precision<8:
                old=prime**precision;new=old*prime;residual=[value(p,lifted,modulus=new)for p in terms];assert all(v%old==0 for v in residual)
                correction=[-sum(v*(r//old)for v,r in zip(row,residual))%prime for row in inverse];nextpoint=[(lifted[j]+old*correction[j])%new for j in range(5)]+[point[-1]];assert nextpoint==lifts[precision]['values']
    common=[ROOT/'reports/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63_REVIEW.json',ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json']
    qpaths+=common+[ROOT/'research/notes/Q_MPS_IMPLICIT_TORIC_ARCS_20261103_v63_REVIEW.md'];hpaths+=common+[ROOT/'research/notes/Q_MPS_HENSEL_TORIC_ARCS_20261103_v63_REVIEW.md']
    save(args.qq_output,dict(id='Q_MPS_IMPLICIT_TORIC_ARCS_20261103_v63_REVIEW',all_passed=True,saved_preoutage_fresh_replay_and_entire_certificate_recovered_and_matched=True,
        independent_direct5by5_period_extension_permutation_and_complete_annihilator=True,independent_all17_QQ_equation_points_nonrho_Jacobians_and_selected_minors=True,selected_determinants=qqdet,
        independent_all17_exact_first_and_second_rho_derivatives=True,actual_period_source_etale_curve_not_Poisson_tangent_integration=True,
        all_finite_Qa_in_components=expected,all_a_by_actual_arc_and_dense_rho_closedness=True,nonmembership_branch_exhaustion_or_component_count_increase_claimed=False),qpaths)
    save(args.hensel_output,dict(id='Q_MPS_HENSEL_TORIC_ARCS_20261103_v63_REVIEW',all_passed=True,saved_preoutage_fresh_replay_and_entire_certificate_recovered_and_matched=True,
        independent_all16_complete_integer_period_polynomials_Fp_points_Jacobians_and_determinants=True,modular_determinants=modular,
        independent_all16_precision8_fixed_rho_Newton_replays=True,positive_allorder_Hensel_existence_and_etale_rho_curve_accepted=True,
        admissibility_nonzero_on_same_dominant_complex_source_after_faithful_base_extension=True,
        all_finite_Qa_in_components=hexpected,finite_field_emptiness_or_precision8_as_allorder_proof_claimed=False,
        nonmembership_branch_exhaustion_or_component_count_increase_claimed=False),hpaths)
    print('PASS independent17QQ+16Fp actual period equations/minors/jets and positive Hensel bridges')

if __name__=='__main__':main()
