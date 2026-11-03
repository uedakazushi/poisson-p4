#!/usr/bin/env python3
"""Independent small coefficient and one-old-point residual-curve review.

No author module, new beta point, high-degree probe, or Groebner calculation.
"""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

from q_firstjet_spin4_mixed_discriminant_review_20261103_v63 import mat,mm,det,inv,add

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW'

def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zero(n,m):return [[Q(0)for _ in range(m)]for _ in range(n)]
def ident(n):return [[Q(i==j)for j in range(n)]for i in range(n)]
def trace(M):return sum(M[i][i]for i in range(len(M)))

def main(output):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63/report.json'
    certpath=author.parent/'certificate.json.gz';ar,cert=read(author),read(certpath);assert ar['all_passed']
    pairpath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    pairreview=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/report.json'
    normalpath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63/certificate.json.gz'
    normalreview=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63_REVIEW/report.json'
    dompath=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/certificate.json.gz'
    domreview=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63_REVIEW/report.json'
    independent_dom=domreview.parent/'independent-certificate.json.gz'
    invariantpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    checked={}
    for p in [author,pairreview,normalreview,domreview]:
        data=read(p);assert data['all_passed']
        for key in ['input_sha256','evidence_sha256']:
            for name,want in data.get(key,{}).items():
                found=next((q for q in [ROOT/name,p.parent/name]if q.is_file()and sha(q)==want),None)
                assert found is not None,(p,name);checked[str(found.relative_to(ROOT))]=want
    pair,normal,dom,own_dom,information=map(read,[pairpath,normalpath,dompath,independent_dom,invariantpath])
    x,y,a,g,h,i,j,t,z,k1,k2,m,R,T=s.symbols('x y a g h i j t z k1 k2 m R T')
    delta=x*x-y*y;S=x*x+y*y;ell=(a+x)**2-y*y
    d=-(a*x+S)/y;c=m-ell;f=m-x*x*ell/y**2
    C=s.Matrix([[a,1,g,h],[c,-a,i,j],[-j,h,d,1],[i,-g,f,-d]])
    D=s.diag(x,-x,y,-y)
    #24-term Leibniz expansion, independent of the author's matrix determinant backend.
    def determinant4(M):
        answer=0
        for permutation in itertools.permutations(range(4)):
            parity=sum(permutation[ii]>permutation[jj]for ii in range(4)for jj in range(ii+1,4))%2
            answer+=(-1)**parity*s.prod(M[row,permutation[row]]for row in range(4))
        return s.expand(answer)
    original=C.subs(m,g*j-h*i)
    pencil=s.Poly(determinant4(t*original+z*D),t,z)
    p=[s.cancel(pencil.coeff_monomial(t**(4-ii)*z**ii))for ii in range(5)]
    assert all(s.cancel(q-s.sympify(saved))==0 for q,saved in zip(p,cert['full_pencil_coefficients']))
    B2=-4*a*a*x*x-2*a*x**3-6*a*x*y*y+3*x**4+y**4
    assert s.cancel(p[2]-B2-(x+y)**2*g*j+(x-y)**2*h*i)==0
    assert s.cancel(p[3]+2*a*x*delta+2*x*x*S)==0 and p[4]==x*x*y*y
    kappa=k2-B2;u=(kappa-(x-y)**2*m)/(4*x*y);v=(kappa-(x+y)**2*m)/(4*x*y)
    assert s.cancel(u-v-m)==0
    A0=a*d*(a*y+d*x)+(a*y+d*x)*(u+v)+c*d*y+f*a*x
    K1=k1/2-A0
    B0=(a*a+c)*(d*d+f)+m*m+2*a*d*(u+v)-delta**2
    Q1=-c*y*R*T*T-x*v*R*R*T+f*x*u*T+y*u*v*R-K1*R*T
    Q2=c*R*R*T*T-c*f*T*T-v*v*R*R+f*u*u+B0*R*T-2*c*d*R*T*T+2*a*f*u*T-2*a*v*R*R*T+2*d*u*v*R
    assert s.cancel(Q1-s.sympify(cert['polynomial22_equation1']))==0
    assert s.cancel(Q2-s.sympify(cert['polynomial22_equation2']))==0
    #Pull back to the genuine double cover. g,h remain arbitrary invertible coordinates.
    pull={R:g/h,T:g*h}
    source={i:v/h,j:u/g}
    assert s.cancel((p[1]-k1).subs({m:g*j-h*i}).subs(source,simultaneous=True)-2*(Q1/(R*T)).subs(pull,simultaneous=True))==0
    assert s.cancel((determinant4(C)-delta**2).subs(source,simultaneous=True)-(Q2/(R*T)).subs(pull,simultaneous=True))==0
    assert s.expand(y*y*(p[0]-delta**2)-s.sympify(normal['final_polynomial']))==0
    #Check bidegrees and every homogeneous boundary coefficient, not only corners.
    R0,R1,T0,T1=s.symbols('R0 R1 T0 T1')
    def homogeneous(q):
        poly=s.Poly(s.expand(q),R,T);assert poly.degree(R)==poly.degree(T)==2
        return s.expand(sum(co*R1**er*R0**(2-er)*T1**et*T0**(2-et)for (er,et),co in poly.terms()))
    H1,H2=map(homogeneous,[Q1,Q2])
    edges=[({R0:1,R1:0},f*x*u*T0*T1),({R0:0,R1:1},-x*v*T0*T1),
           ({T0:1,T1:0},y*u*v*R0*R1),({T0:0,T1:1},-c*y*R0*R1)]
    for substitution,expected in edges:assert s.cancel(H1.subs(substitution)-expected)==0
    corners=[{R0:1,R1:0,T0:1,T1:0},{R0:1,R1:0,T0:0,T1:1},
             {R0:0,R1:1,T0:1,T1:0},{R0:0,R1:1,T0:0,T1:1}]
    for substitution,expected in zip(corners,[f*u*u,-c*f,-v*v,c]):
        assert H1.subs(substitution)==0 and s.cancel(H2.subs(substitution)-expected)==0
    assert s.cancel(s.diff(Q1,k1)+R*T/2)==0 and s.diff(Q2,k1)==0
    #m identity in the FULL diagonal-D chart, with a,c,d,f independently free.
    cc,dd,ff=s.symbols('cc dd ff');general=C.subs({c:cc,d:dd,f:ff},simultaneous=True)
    #The graph expressions overlap; construct directly to retain true independent entries.
    general=s.Matrix([[a,1,g,h],[cc,-a,i,j],[-j,h,dd,1],[i,-g,ff,-dd]])
    intrinsic=(s.trace(general*D*D*general*D*D)-S*s.trace(general*general*D*D)+x*x*y*y*s.trace(general*general))/(2*delta*delta)
    assert s.cancel(intrinsic-g*j+h*i)==0
    #One old QQ point: own forward matrix-dual arithmetic, not term replacement.
    C0,D0=mat(pair['actual_C4']),mat(pair['actual_D4']);basis=list(map(mat,pair['pair_coordinate_basis4']))
    def dualtrace(word,E=None,letter=None):
        value=ident(4);variation=zero(4,4)
        for ch in word:
            M=C0 if ch=='C'else D0;change=E if ch==letter else zero(4,4)
            variation=add(mm(variation,M),mm(value,change));value=mm(value,M)
        return trace(value),trace(variation)
    words=['CC','DD','CCDD','CDDCDD','DDDD']
    qC,qD,q4,q6,qD4=[dualtrace(w)[0]for w in words]
    SD=qD/2;PD=qD*qD/8-qD4/4;dis=SD*SD-4*PD;num=q6-SD*q4+PD*qC
    assert PD==det(D0)and dis and num/(2*dis)==Q(cert['actual_same_point_intrinsic_m'])
    grad=[]
    for letter in ['C','D']:
        for E in basis:
            dqC,dqD,dq4,dq6,dqD4=[dualtrace(w,E,letter)[1]for w in words]
            dS=dqD/2;dP=qD*dqD/4-dqD4/4;ddis=2*SD*dS-4*dP
            dnum=dq6-dS*q4-SD*dq4+dP*qC+PD*dqC
            grad.append((dnum*dis-num*ddis)/(2*dis*dis))
    assert [grad]==mat(cert['intrinsic_m_gradient1x20'])
    beta=mm([grad],mat(pair['jacobian20x16']));assert beta==mat(cert['intrinsic_m_gradient1x16_beta'])
    five=mat(own_dom['independent_trace8_direct_beta_gradient5x16']);six=five+beta
    assert six==mat(cert['universal_six_gradient6x16'])
    cols=cert['selected_beta_columns'];minor=[[row[col]for col in cols]for row in six];dminor=det(minor)
    assert len(cols)==6 and minor==mat(cert['exact6_minor'])and dminor==Q(cert['exact6_minor_determinant'])and dminor
    assert abs(dminor.numerator).bit_length()==544 and dminor.denominator.bit_length()==374
    b0=mat([[q]for q in dom['beta16']]);orb=[]
    orbitcols=[mm(mat(A),b0)for A in information['all_ten_representation_matrices']]
    orbit=[[column[row][0]for column in orbitcols]for row in range(16)]
    assert mm(six,orbit)==zero(6,10)
    new=mm(six,mat(dom['full_beta16x18_source_Jacobian']));assert new==mat(cert['same_NEW_source_six_gradient6x18'])
    selected=dom['selected_source_columns'];alpha=mm([[new[5][j]for j in selected]],inv([[row[j]for j in selected]for row in new[:5]]))
    assert mm(alpha,new[:5])==[new[5]]and alpha==mat(cert['NEW_source_m_differential_linear_combination'])
    evidence={
        'all_five_pencil_coefficients_own24_term_determinant':True,
        'both_exact_quotient_equations_and_all_homogeneous_edges_corners':True,
        'generic_coprimeness_coefficient_minus_RT_over2':True,
        'full_diagonal_chart_intrinsic_m_identity':True,
        'own_forward_dual_Fraction_gradient20_and_beta16':True,
        'same_existing_point_index':0,'rank6_minor_columns':cols,'exact_QQ_rank6_minor':str(dminor),
        'all10_orbit_derivatives_zero':True,'NEW_source_rank_stays5_not6':True,
        'proof_of_degree_at_most8_in_review_note':True,
    }
    ev=output/'independent-certificate.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,author,certpath,pairpath,pairreview,normalpath,normalreview,dompath,domreview,independent_dom,invariantpath,
           ROOT/'code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py',
           ROOT/'research/notes/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63.md']
    report={
        'id':ID,'all_passed':True,
        'full_three_actual_fixed_base_equations_and_relative_sign_quotient_accepted':True,
        'generic_P1xP1_proper_intersection_scheme_total_length8_accepted':True,
        'intrinsic_rational_m_and_universal_six_invariants_independent_accepted':True,
        'selected_actual_component_m_degree_at_most8_after_ordered_finite_base_and_relative_sign_quotient':True,
        'exact_degree8_reduced_eight_points_irreducibility_genus_or_Prym_identification_claim':False,
        'NEW25_multisection_degree_or_equation_claim':False,
        'one_existing_point_only_no_GB_high_degree_probe_or_author_rerun':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{ev.name:sha(ev)},'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-start,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'m_degree_upper':8,'seconds':report['seconds']}))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
