#!/usr/bin/env python3
"""One actual quadratic-curve attempt on the unique surviving tangent line.
This does not infer a formal solution or increase a jet when the polynomial
curve fails. The point is an explicit diagonal image of the frozen witness.
"""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import classification_rank3_mixed_component_review as engine
ROOT=Path(__file__).resolve().parents[1]

def main():
    begun=time.monotonic();xx=s.symbols('h x1 x2 x3 x4');h,x,y,z,w=xx;u=s.Symbol('u')
    path=ROOT/'reports/CLASSIFICATION_NIL22_ONE_QUADRATIC_SUPPORT_20261107.json'
    packet=json.loads(path.read_text());env={str(q):q for q in xx}
    def decode(data):
        P=s.zeros(5)
        for key,value in data.items():
            i,j=map(int,key.split(','));P[i,j]=s.sympify(value,locals=env)
        return P-P.T
    old=decode(packet['canonical_brackets']);oldR=decode(packet['normal_brackets'][0])+318*decode(packet['normal_brackets'][1])
    mu=s.Rational(1,159);M=s.diag(1,mu**2,mu,mu,1);inv=M.inv()*s.Matrix(xx)
    def transform(P):return (M*P*M.T).subs(dict(zip(xx,inv)),simultaneous=True).applyfunc(s.expand)
    P=transform(old)/mu;R=transform(oldR)
    expected=decode({'0,3':'-x2*x4/2','1,3':'2*x2*x3','1,4':'h*x2+2*x2*x4','2,3':'2*h*x3-x1*x4/2','2,4':'h**2+2*h*x4','3,4':'h**2-x2*x4/2+x3*x4+2*x4**2'})
    assert P[2,3]==24*(h*h+y*w) and R==expected
    def sparse(A):return {ij:dict(s.Poly(A[ij],*xx).terms()) for ij in engine.PAIRS if A[ij]!=0}
    assert engine.jacobi(sparse(P))==s.zeros(350,1)
    J=engine.jacobi(sparse(P+u*R)).applyfunc(s.expand)
    assert all(s.Poly(q,u).coeff_monomial(u)==0 for q in J)
    Rjac=engine.jacobi(sparse(R));Q=s.zeros(5);minor=None
    if Rjac!=s.zeros(350,1):
        Dcols=[];ps=sparse(P)
        for ij,e in engine.COLUMNS:
            v={k:dict(f) for k,f in ps.items()};v.setdefault(ij,{})[e]=v.get(ij,{}).get(e,0)+1
            Dcols.append(engine.jacobi(v))
        D=s.Matrix.hstack(*Dcols);indices=packet['implicit_coordinate_indices'];rows=packet['selected_equation_rows']
        block=DomainMatrix.from_Matrix(D[rows,indices]).convert_to(s.QQ)
        minor=block.det();assert minor
        correction=-block.inv().to_Matrix()*Rjac[rows,:]
        qv=s.zeros(150,1)
        for i,index in enumerate(indices):qv[index]=correction[i]
        assert D*qv+Rjac==s.zeros(350,1)
        for ((i,j),e),a in zip(engine.COLUMNS,qv):Q[i,j]+=a*s.prod(q**v for q,v in zip(xx,e))
        Q-=Q.T
    curve=P+u*R+u*u*Q;actual=engine.jacobi(sparse(curve)).applyfunc(s.expand)
    success=actual==s.zeros(350,1)
    def encode(A):return {str(i)+','+str(j):str(s.factor(A[i,j])) for i,j in engine.PAIRS if A[i,j]!=0}
    out={'task_id':'E_multiple_hyperplane','status':'EXACT_ACTUAL_QUADRATIC_POISSON_CURVE' if success else 'QUADRATIC_POLYNOMIAL_CURVE_NOT_POISSON',
         'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,Path(__file__),Path(engine.__file__)]},
         'producer_bytes':Path(__file__).read_text(),'diagonal_transport':list(map(str,M.diagonal())),
         'point_scalar_factor':'159','normal_scalar_after_transport':'1/159','same_frozen_witness_exact_transport':True,
         'point_brackets':encode(P),'normal_brackets':encode(R),'normal_self_Jacobi_zero':Rjac==s.zeros(350,1),
         'quadratic_correction_brackets':encode(Q),'implicit_minor_determinant':str(minor),
         'all350_quadratic_coefficient_zero':True,'full_curve_brackets':encode(curve),
         'whole350_actual_Jacobi_zero':success,
         'all_nonzero_actual_Jacobi_coefficients':{str(i):str(s.factor(v)) for i,v in enumerate(actual) if v},
         'scope':'One genuine degree-two polynomial curve was tested. If it fails, no formal integration or new componenthood is inferred.',
         'seconds':time.monotonic()-begun}
    target=ROOT/'reports/CLASSIFICATION_NIL22_ONE_ACTUAL_CURVE_20261107.json'
    if target.exists():raise FileExistsError(str(target))
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'seconds':out['seconds'],'nonzero_residual_rows':len(out['all_nonzero_actual_Jacobi_coefficients'])}))


def cubic_support():
    begun=time.monotonic();xx=s.symbols('h x1 x2 x3 x4');env={str(v):v for v in xx}
    source_path=ROOT/'reports/CLASSIFICATION_NIL22_ONE_ACTUAL_CURVE_20261107.json'
    packet=json.loads(source_path.read_text())
    def decode(data):
        M=s.zeros(5)
        for key,value in data.items():
            i,j=map(int,key.split(','));M[i,j]=s.sympify(value,locals=env)
        return M-M.T
    def sparse(M):return {ij:dict(s.Poly(s.expand(M[ij]),*xx).terms()) for ij in engine.PAIRS if M[ij]!=0}
    def encode(M):return {str(i)+','+str(j):str(s.factor(M[i,j])) for i,j in engine.PAIRS if M[i,j]!=0}
    def tensor(v):
        M=s.zeros(5)
        for ((i,j),e),a in zip(engine.COLUMNS,v):M[i,j]+=a*s.prod(x**k for x,k in zip(xx,e))
        return M-M.T
    def dm(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ)
    P=decode(packet['point_brackets']);R=decode(packet['normal_brackets']);Q=decode(packet['quadratic_correction_brackets'])
    J=lambda M:engine.jacobi(sparse(M))
    assert J(P)==s.zeros(350,1)
    Dcols=[];ps=sparse(P)
    for ij,e in engine.COLUMNS:
        value={k:dict(p) for k,p in ps.items()};value.setdefault(ij,{})[e]=value.get(ij,{}).get(e,0)+1
        Dcols.append(engine.jacobi(value))
    D=s.Matrix.hstack(*Dcols);K=dm(D).nullspace().to_Matrix().T
    assert K.shape==(150,49) and D*K==s.zeros(350,49)
    assert D*engine.vector(sparse(Q))+J(R)==s.zeros(350,1)
    JR=J(R);C=J(R+Q)-JR-J(Q)
    kernel_tensors=[tensor(K[:,i]) for i in range(49)]
    mixed=s.Matrix.hstack(*[J(R+T)-JR-J(T) for T in kernel_tensors])
    M=D.row_join(mixed)
    def project(matrix,C):
        cols=list(dm(matrix).rref()[1]);basis=matrix[:,cols]
        rows=list(dm(basis.T).rref()[1]);block=dm(basis[rows,:]);inv=block.inv().to_Matrix()
        rem=C-basis*inv*C[rows,:]
        if rem==s.zeros(350,1):return dict(rank=len(cols),nonzero=False),None
        nz=next(i for i,value in enumerate(rem) if value)
        lam=s.zeros(1,350);lam[0,nz]=1;selected=basis[nz,:]*inv
        for j,row in enumerate(rows):lam[0,row]-=selected[0,j]
        assert lam*matrix==s.zeros(1,matrix.cols)
        value=(lam*C)[0];assert value
        return dict(rank=len(cols),nonzero=True,left_functional={str(i):str(lam[0,i]) for i in range(350) if lam[0,i]},
                    nonzero_cubic_pairing=str(value),minor_rows=rows,minor_columns=cols,minor_determinant=str(block.det())),lam
    result,lam=project(M,C)
    direct=None
    if result['nonzero']:
        oldpath=ROOT/'reports/CLASSIFICATION_NIL22_ONE_QUADRATIC_SUPPORT_20261107.json';old=json.loads(oldpath.read_text())
        T=decode(old['normal_brackets'][1]);mu=s.Rational(1,159);A=s.diag(1,mu**2,mu,mu,1);inv=A.inv()*s.Matrix(xx)
        T=(A*T*A.T).subs(dict(zip(xx,inv)),simultaneous=True).applyfunc(s.expand)
        assert D*engine.vector(sparse(T))==s.zeros(350,1)
        direct,_=project(M.row_join(J(T)),C)
    out={'task_id':'E_multiple_hyperplane','status':'CUBIC_SUPPORT_OBSTRUCTION_NONZERO' if result['nonzero'] else 'CUBIC_AMBIGUITY_QUOTIENT_ZERO',
         'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source_path,Path(__file__),Path(engine.__file__)]},
         'producer_bytes':Path(__file__).read_text(),'full_D_rank':101,'full_kernel_dimension':49,
         'all_kernel_vectors_killed_by_full_D':True,'second_order_equation_exact':True,
         'kernel_basis_brackets':[encode(T) for T in kernel_tensors],
         'cubic_B_R_Q':{str(i):str(v) for i,v in enumerate(C) if v},
         'all_second_order_ambiguities_quotient':result,
         'including_remaining_full_normal_quadratic_image_quotient':direct,
         'scope':'Nonzero quotient kills every second-order correction Q+kerD. Together with the quadratic initial support (a²,ab), the Puiseux valuation argument proves reduced support is a point. A literal cubic initial generator is asserted only if the strengthened quotient remains nonzero.',
         'exact_nilpotence_thickness_not_claimed':True,'seconds':time.monotonic()-begun}
    path=ROOT/'reports/CLASSIFICATION_NIL22_ONE_CUBIC_SUPPORT_20261107.json'
    if path.exists():raise FileExistsError(str(path))
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'quotient_rank':result['rank'],'direct_cubic':None if direct is None else direct['nonzero'],'seconds':out['seconds']}))

def aff32_limit():
    """An actual Laurent PGL degeneration of the full v14 Aff32 family."""
    import zipfile
    begun=time.monotonic()
    xx=s.symbols('h x1 x2 x3 x4');h,x,y,z,w=xx
    eps,B,t,b=s.symbols('eps B t b')
    # The old coordinates are (h,y,x,w,z); matrices below use (h,x,y,z,w).
    N=s.Matrix([0,y,h,w,0]);H=s.Matrix([0,2*x,y,(b+1)*z,b*w])
    original=N*H.T-H*N.T
    eta=s.zeros(5);theta=s.zeros(5)
    def entry(M,i,j,value):M[i,j]=value;M[j,i]=-value
    entry(eta,0,2,-(2*b*b+b-1)*h*y)
    entry(eta,0,1,-(b*b+b-1)*h*x-b*y*y)
    entry(eta,0,4,-(2*b-1)*h*w)
    entry(eta,0,3,-(b-1)*y*w)
    entry(eta,2,1,(b*b-4*b+2)*y*x)
    entry(eta,1,3,(b-2)*x*z)
    entry(theta,0,2,b*(b-1)*(b+1)*y*y)
    entry(theta,0,1,b*(b-1)**2*y*x)
    entry(theta,0,4,b*(b-1)*y*w)
    original+=t*eta+t*t*theta
    def sparse(M):return {ij:dict(s.Poly(s.expand(M[ij]),*xx).terms()) for ij in engine.PAIRS if M[ij]!=0}
    def encode(M):return {str(i)+','+str(j):str(s.factor(M[i,j])) for i,j in engine.PAIRS if M[i,j]!=0}
    assert all(s.expand(q)==0 for q in engine.jacobi(sparse(original)))
    # old h=H, old x=eps^2 X, old y=eps Y, old z=H+eps Z, old w=W.
    A=s.diag(1,eps**2,eps,eps,1);A[3,0]=1
    inverse=A.inv();old_values=A*s.Matrix(xx)
    transported=eps*inverse*original.subs(dict(zip(xx,old_values)),simultaneous=True)*inverse.T
    transported=transported.subs(b,-1+eps*B).applyfunc(s.cancel)
    polynomial=all(s.denom(q)==1 for q in transported)
    assert polynomial
    limit=transported.subs(eps,0).applyfunc(s.expand)
    # New Z'=Z+tY and sign -1 give precisely A1(t_A=4t,e_A=-B).
    C=s.eye(5);C[3,2]=t
    correction=C.inv()*s.Matrix(xx)
    normalized=-(C*limit*C.T).subs(dict(zip(xx,correction)),simultaneous=True).applyfunc(s.expand)
    expected=s.zeros(5)
    entry(expected,1,2,2*h*x-y*y)
    entry(expected,1,3,4*t*(h*x-y*y/2)+2*x*w-B*h*y)
    entry(expected,1,4,y*w)
    entry(expected,2,3,-B*h*h+y*w)
    entry(expected,2,4,h*w)
    entry(expected,3,4,w*w+4*t*h*w)
    assert all(s.expand(q)==0 for q in normalized-expected), normalized-expected
    assert all(s.expand(q)==0 for q in engine.jacobi(sparse(limit)))
    def pfaffian(M):
        ans=0
        for i in range(5):
            j,k,l,m=[v for v in range(5) if v!=i]
            ans+=(-1)**i*xx[i]*(M[j,k]*M[l,m]-M[j,l]*M[k,m]+M[j,m]*M[k,l])
        return s.factor(ans)
    actual_pf=pfaffian(normalized)
    assert s.expand(actual_pf-2*t*h*h*w*(2*h*x-y*y))==0
    source=ROOT/'sources/v14.zip'
    with zipfile.ZipFile(source) as bundle:
        oldcode=bundle.read('code/rank2_lower_bounds/verify_aff32_deformation.py')
    canonical=ROOT/'reports/CLASSIFICATION_NIL22_ONE_CANONICAL_20261107.json'
    cp=json.loads(canonical.read_text());env={str(v):v for v in xx}
    env.update({'n7':s.Integer(159),'n13':s.Integer(1),'e':s.Integer(159),'t':s.Integer(1)})
    # The point comparison is independent of source symbol naming.
    frozen=24*normalized.subs({t:s.Rational(1,4),B:-159})
    qp=ROOT/'reports/CLASSIFICATION_NIL22_ONE_QUADRATIC_SUPPORT_20261107.json'
    qdata=json.loads(qp.read_text());oldpoint=s.zeros(5)
    for key,value in qdata['canonical_brackets'].items():
        i,j=map(int,key.split(','));entry(oldpoint,i,j,s.sympify(value,locals=env))
    assert all(s.expand(q)==0 for q in frozen-oldpoint), frozen-oldpoint
    out={'task_id':'E_multiple_hyperplane','status':'ENTIRE_A_ONE_SOURCE_IS_AFF32_COMPONENT_BOUNDARY',
         'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,canonical,qp,Path(__file__),Path(engine.__file__)]},
         'producer_bytes':Path(__file__).read_text(),
         'historical_aff32_producer_entry':'code/rank2_lower_bounds/verify_aff32_deformation.py',
         'historical_aff32_producer_sha256':hashlib.sha256(oldcode).hexdigest(),
         'aff32_original_brackets':encode(original),'whole_aff32_projective_Jacobi_zero':True,
         'old_coordinates_in_new_order_h_x_y_z_w':list(map(str,old_values)),
         'coordinate_determinant':str(A.det()),'aff32_b':'-1+eps*B','aff32_t':'t','tensor_scale':'eps',
         'all150_transformed_coefficients_polynomial_in_eps':polynomial,
         'whole_transformed_family_brackets':encode(transported),
         'limit_brackets':encode(limit),'all350_limit_projective_Jacobi_zero':True,
         'final_coordinate_change':{'new_z':'old_z+t*y'},'final_tensor_scale':'-1',
         'normalized_limit_brackets':encode(normalized),'whole_A_one_parameter_identity':True,
         'parameter_map':{'A_one_t':'4*t','A_one_e':'-B'},'actual_projective_pfaffian':str(actual_pf),
         'frozen_witness':{'aff32_t':'1/4','B':'-159','scale_to_original_canonical':'24','all150_identity':True},
         'known_family_dependency':'The full v14 Aff32 family lies in its 23-dimensional rank-four component T19; its component proof and SNC degeneration are prior accepted inputs, not reproved by this coefficient identity.',
         'coverage_logic':'For every t_A!=0 and every e_A choose t=t_A/4,B=-e_A. For eps!=0 the displayed actual family is a PGL transform and scalar of Aff32; its polynomial eps=0 endpoint is the full canonical A1 source. Closedness gives all endpoints, including e_A=0. No local smoothness or finite-jet integration is used.',
         'seconds':time.monotonic()-begun}
    path=ROOT/'reports/CLASSIFICATION_NIL22_ONE_AFF32_DEGENERATION_20261107.json'
    if path.exists():raise FileExistsError(str(path))
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'seconds':out['seconds'],'Pf':out['actual_projective_pfaffian']}))

if __name__=='__main__':
    import sys
    if '--aff32-limit' in sys.argv:aff32_limit()
    elif '--cubic-support' in sys.argv:cubic_support()
    else:main()
