#!/usr/bin/env python3
"""Whole long-top B source deformation splitting its repeated hyperplane.
The construction uses a closed rational one-form on P3. All identities are
checked as polynomial identities in every source parameter and epsilon.
"""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main():
    begun=time.monotonic();h,x,y,z,w=s.symbols('h x1 x2 x3 x4');xx=(h,x,y,z)
    a,c,f,ep=s.symbols('a c f epsilon');E=s.Matrix(xx)
    R=-h*y+z*z/2;H=h+ep*z;Q=a*H*H+f*R
    dh=s.Matrix([1,0,0,0]);dz=s.Matrix([0,0,0,1]);dH=s.Matrix([s.diff(H,v) for v in xx]);dR=s.Matrix([s.diff(R,v) for v in xx])
    omega=(Q*(h*dz-z*dh)+c*h*(H*dR-2*R*dH)).applyfunc(s.expand)
    B=s.Matrix(4,4,lambda i,j:s.expand((s.diff(omega[j],xx[i])-s.diff(omega[i],xx[j]))/4))
    Pi=s.zeros(4)
    for i,j,k,l,sign in [(0,1,2,3,1),(0,2,1,3,-1),(0,3,1,2,1),(1,2,0,3,1),(1,3,0,2,-1),(2,3,0,1,1)]:
        Pi[i,j]=sign*B[k,l];Pi[j,i]=-Pi[i,j]
    N=s.Matrix([0,y,z,h])+ep*s.Matrix([0,0,2*y,z]);K=s.Matrix([0,h,0,0])
    A=-N+3*ep*E/4+5*f*K/4
    P=(Pi+E*A.T-A*E.T).applyfunc(s.expand)
    def jac(M):return [s.expand(sum(M[i,l]*s.diff(M[j,k],xx[l])+M[j,l]*s.diff(M[k,i],xx[l])+M[k,l]*s.diff(M[i,j],xx[l]) for l in range(4))) for i,j,k in itertools.combinations(range(4),3)]
    assert jac(P)==[0]*4 and jac(Pi)==[0]*4
    assert all(s.expand(sum(s.diff(Pi[i,j],xx[i]) for i in range(4)))==0 for j in range(4))
    assert sum(s.diff(A[i],xx[i]) for i in range(4))==0
    lie=Pi*s.zeros(4)
    for i,j in itertools.combinations(range(4),2):
        value=sum(A[k]*s.diff(Pi[i,j],xx[k])-Pi[k,j]*s.diff(A[i],xx[k])-Pi[i,k]*s.diff(A[j],xx[k]) for k in range(4))
        assert s.expand(value)==0
    pf=s.factor(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
    assert s.expand(pf+h*H*Q)==0
    assert all(s.Poly(P[i,j],*xx).total_degree()==2 for i,j in itertools.combinations(range(4),2) if P[i,j]!=0)
    source_path=ROOT/'reports/CLASSIFICATION_NIL31_LONG_B_CANONICAL_20261107.json';source=json.loads(source_path.read_text())
    env={'h':h,'x1':x,'x2':y,'x3':z,'x4':w,'n0':a,'n2':c,'n3':0,'n5':f,'n9':0}
    C=s.zeros(5)
    for key,value in source['canonical_brackets'].items():
        i,j=map(int,key.split(','));C[i,j]=s.expand(s.sympify(value,locals=env)/2)
    C-=C.T
    v=s.Matrix([0,f*h-y,-z,-h,0]);EE=s.Matrix([h,x,y,z,w]);lift=s.zeros(5)
    lift[:4,:4]=P.subs(ep,0)
    assert all(s.expand(v)==0 for v in lift-C-(EE*v.T-v*EE.T))
    Qmatrix=s.hessian(Q,(h,y,z))/2;det=s.factor(Qmatrix.det())
    assert s.expand(det+f*f*(2*a*ep*ep+f)/8)==0
    # The rational expression is only an explanatory identity on its open;
    # the displayed omega and P themselves are polynomial everywhere.
    numerator=h*H*Q
    for i,var in enumerate(xx):
        rational=(s.diff(H,var)/H-s.diff(h,var)/h)/ep+(c/f)*(s.diff(Q,var)/Q-2*s.diff(H,var)/H)
        assert s.cancel(omega[i]-numerator*rational)==0
    assert s.expand((omega.T*N)[0]-numerator)==0
    here=Path(__file__)
    result={'task_id':'E_multiple_hyperplane','all_passed':True,'coefficient_ring':'QQ[a,c,f,epsilon,h,x1,x2,x3]',
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source_path,here]},
            'producer_bytes':here.read_text(),'quadratic_R':str(R),'split_hyperplane_H':str(H),'quadric_Q':str(s.expand(Q)),
            'omega_coefficients':list(map(str,omega)),'canonical_divergence_free_bivector':{str(i)+','+str(j):str(Pi[i,j]) for i,j in itertools.combinations(range(4),2) if Pi[i,j]!=0},
            'linear_symmetry_A':list(map(str,A)),
            'full_homogeneous_affine_bivector':{str(i)+','+str(j):str(P[i,j]) for i,j in itertools.combinations(range(4),2) if P[i,j]!=0},
            'whole_source_universal_strict_Jacobi_zero':True,'canonical_part_strict_Jacobi_and_divergence_zero':True,
            'linear_symmetry_preserves_canonical_part':True,'epsilon_zero_equals_whole_original_source_mod_Euler':True,
            'epsilon_zero_Euler_gauge_vector':list(map(str,v)),
            'actual_affine_pfaffian':str(pf),'projective_divisor':'{x4=0}+{h=0}+{h+epsilon*x3=0}+{Q=0}',
            'quadric_rank3_minor_determinant':str(det),
            'covered_open':'f*epsilon*(2*a*epsilon**2+f) !=0',
            'complements':'For every f!=0 choose a dense epsilon open and take epsilon->0. The f=0 rank-four source boundary is included by closure of the entire irreducible source. a=f=0 is rank at most2.',
            'coverage_dependency':'THREE_HYPERPLANES_QUADRIC_SUPPORT.md and its complete independent/dependent incidence extensions: all three distinct hyperplanes plus an irreducible quadric.',
            'whole_long_top_B_source_covered':True,'seconds':time.monotonic()-begun}
    out=ROOT/'reports/CLASSIFICATION_NIL31_LONG_B_DEFORMATION_20261107.json'
    if out.exists():raise FileExistsError(str(out))
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'pfaffian':str(pf),'quadric_minor':str(det)}))


def half_extra():
    begun=time.monotonic();x,y,z,w=s.symbols('x1 x2 x3 x4');xx=(x,y,z,w);h=s.Symbol('h')
    a,t,ep=s.symbols('a s epsilon')
    entries=[2*a*(1+ep*y),-y*z+4*t*y*w-2*t*z,w*(2*t-y),-z*(1+ep*y),-w*(1+ep*y),-w*w-ep*z*w]
    P=s.zeros(4)
    for (i,j),v in zip(itertools.combinations(range(4),2),entries):P[i,j]=v;P[j,i]=-v
    jac=[s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l]) for l in range(4))) for i,j,k in itertools.combinations(range(4),3)]
    assert jac==[0]*4
    homogeneous=P.applyfunc(lambda v:s.expand(h*h*v.subs({q:q/h for q in xx},simultaneous=True)))
    assert all(s.Poly(v,h,*xx).total_degree()==2 for v in homogeneous if v!=0)
    Q=a*h*w+(a*ep+2*t)*h*z-2*t*y*w
    pf=s.factor(homogeneous[0,1]*homogeneous[2,3]-homogeneous[0,2]*homogeneous[1,3]+homogeneous[0,3]*homogeneous[1,2])
    assert s.expand(pf+2*w*(h+ep*y)*Q)==0
    det=s.factor((s.hessian(Q,(h,y,z,w))/2).det())
    assert s.expand(det-t*t*(a*ep+2*t)**2/4)==0
    source_path=ROOT/'reports/CLASSIFICATION_NIL22_HALF_EXTRA_CANONICAL_20261107.json'
    source=json.loads(source_path.read_text())
    env={str(q):q for q in (h,)+xx};env.update({'n'+str(i):s.Integer(0) for i in range(15)});env['n0']=a;env['n5']=t
    for key,value in source['canonical_brackets'].items():
        i,j=map(int,key.split(','));assert s.expand(homogeneous[i-1,j-1].subs(ep,0)-s.sympify(value,locals=env))==0
    # Independently check the elementary commuting-field description on w!=0.
    r=s.Symbol('r');vars3=(y,r,w)
    U=s.Matrix([1+ep*y,1+ep*r,0]);W=s.Matrix([0,0,w])
    Z=s.Matrix([2*a*(1+ep*y),4*t*(y-r),(2*t-y)*w])
    R=-U*W.T+W*U.T
    for i,j in itertools.combinations(range(3),2):
        lie=sum(Z[k]*s.diff(R[i,j],vars3[k])-R[k,j]*s.diff(Z[i],vars3[k])-R[i,k]*s.diff(Z[j],vars3[k]) for k in range(3))
        assert s.expand(lie)==0
    result={'task_id':'E_multiple_hyperplane','all_passed':True,
        'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source_path,Path(__file__)]},
        'producer_bytes':Path(__file__).read_text(),'coefficient_ring':'QQ[a,s,epsilon,x1,x2,x3,x4]',
        'full_affine_brackets':list(map(str,entries)),
        'full_homogeneous_h_Casimir_brackets':{str(i+1)+','+str(j+1):str(homogeneous[i,j]) for i,j in itertools.combinations(range(4),2)},
        'all_four_universal_strict_Jacobi_zero':True,'epsilon_zero_is_whole_normalized_source':True,
        'three_variable_rank_two_part_preserved_by_extension_vector_field':True,
        'quadric_Q':str(Q),'actual_affine_homogeneous_pfaffian':str(pf),'rank4_quadric_matrix_determinant':str(det),
        'covered_open':'s*epsilon*(a*epsilon+2*s) !=0',
        'projective_actual_divisor':'{h=0}+{x4=0}+{h+epsilon*x2=0}+{Q=0}',
        'complements':'For each s!=0, choose a dense epsilon open and take epsilon->0. The s=0 original source is the already T01-covered half n5zero source.',
        'coverage_dependency':'Complete three-distinct-hyperplane plus irreducible-quadric support theorem',
        'whole_nil22_half_extra_source_covered':True,'seconds':time.monotonic()-begun}
    out=ROOT/'reports/CLASSIFICATION_NIL22_HALF_EXTRA_DEFORMATION_20261107.json'
    if out.exists():raise FileExistsError(str(out))
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'pfaffian':str(pf),'quadric_minor':str(det)}))

def short_third_main():
    """Independent universal affine check of the semisimple modular pencil."""
    begun=time.monotonic();x,y,z,w=s.symbols('x1 x2 x3 x4');xx=(x,y,z,w)
    h,v,g,ep,lam=s.symbols('h v g epsilon lambda')
    q=2*x*z-y*y+ep*x*x
    entries=[(2*v+1)*x*z-(v+1)*y*y/2+(v+1)*ep*x*x/2,
             v*y*z,x/3+g*y,v*z*(z+ep*x),y/3+g*(z+ep*x),z/3]
    P=s.zeros(4)
    for (i,j),value in zip(itertools.combinations(range(4),2),entries):P[i,j]=value;P[j,i]=-value
    J=[s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l]) for l in range(4))) for i,j,k in itertools.combinations(range(4),3)]
    assert J==[0]*4
    modular=s.Matrix([s.expand(sum(s.diff(P[i,j],xx[i]) for i in range(4))) for j in range(4)])
    assert modular==s.Matrix([y,z+ep*x,0,1])
    N=s.Matrix([y,z+ep*x,0]);xyz=(x,y,z)
    assert s.expand(sum(N[i]*s.diff(q,xyz[i]) for i in range(3)))==0
    assert sum(s.diff(N[i],xyz[i]) for i in range(3))==0
    for i,j in itertools.combinations(range(3),2):
        lie=sum(N[k]*s.diff(P[i,j],xyz[k])-P[k,j]*s.diff(N[i],xyz[k])-P[i,k]*s.diff(N[j],xyz[k]) for k in range(3))
        assert s.expand(lie)==0
    linear=modular.jacobian(xx)
    assert linear.rank()==2
    char=s.factor(linear.charpoly(lam).as_expr());assert s.expand(char-lam*lam*(lam*lam-ep))==0
    assert linear*(linear*linear-ep*s.eye(4))==s.zeros(4)
    homogeneous=P.applyfunc(lambda value:s.expand(h*h*value.subs({var:var/h for var in xx},simultaneous=True)))
    pf=s.factor(homogeneous[0,1]*homogeneous[2,3]-homogeneous[0,2]*homogeneous[1,3]+homogeneous[0,3]*homogeneous[1,2])
    assert s.expand(pf-(3*v+1)*h*z*q/6)==0
    source_path=ROOT/'reports/CLASSIFICATION_NIL31_SHORT_THIRD_MAIN_CANONICAL_20261107.json'
    source=json.loads(source_path.read_text());env={str(var):var for var in (h,)+xx};env.update({'n9':v,'n11':g})
    for key,value in source['canonical_brackets'].items():
        i,j=map(int,key.split(','));assert s.expand(6*homogeneous[i-1,j-1].subs(ep,0)-s.sympify(value,locals=env))==0
    out={'task_id':'E_multiple_hyperplane','all_passed':True,
         'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [source_path,Path(__file__)]},
         'producer_bytes':Path(__file__).read_text(),'coefficient_ring':'QQ[v,g,epsilon,x1,x2,x3,x4]',
         'full_affine_brackets':list(map(str,entries)),
         'full_homogeneous_brackets':{str(i+1)+','+str(j+1):str(homogeneous[i,j]) for i,j in itertools.combinations(range(4),2)},
         'all_four_universal_strict_Jacobi_zero':True,'epsilon_zero_is_whole_canonical_source_scaled_by_one_sixth':True,
         'modular_vector':list(map(str,modular)),'modular_linear_rank':2,
         'modular_characteristic_polynomial':str(char),'modular_minimal_polynomial_divides':'lambda*(lambda**2-epsilon)',
         'modular_semisimple_for':'epsilon!=0','divergence_free_linear_symmetry_preserves_q_and_pi3':True,
         'quadric_q':str(q),'actual_homogeneous_affine_pfaffian':str(pf),
         'covered_open':'epsilon*(3*v+1)!=0',
         'coverage_dependency':'Complete rank-two semisimple modular support theorem in MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md, including a nonzero constant component in the zero eigenspace.',
         'complements':'For every normalized rank-four point (3*v+1)!=0, epsilon=0 follows by closedness of the finite union of known components. Transfer across normalization denominators uses the original irreducible main-source closure, as recorded in the parent source proof.',
         'scope':'An actual polynomial family for every v,g,epsilon; no numerical-point rank or finite-jet integration is used.',
         'seconds':time.monotonic()-begun}
    path=ROOT/'reports/CLASSIFICATION_NIL31_SHORT_THIRD_MAIN_DEFORMATION_20261107.json'
    if path.exists():raise FileExistsError(str(path))
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds'],'Pf':str(pf),'modular_charpoly':str(char)}))

if __name__=='__main__':
    import sys
    if '--short-main' in sys.argv:short_third_main()
    elif '--half-extra' in sys.argv:half_extra()
    else:main()
