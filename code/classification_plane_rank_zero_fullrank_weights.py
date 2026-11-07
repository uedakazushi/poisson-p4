#!/usr/bin/env python3
"""Full 30/9 linear normal-weight operators for the two complete k=1 rank3 symbols.

All tau walls come from an exact polynomial maximal minor; every irreducible
factor is specialized in characteristic zero. No nonlinear Groebner basis.
"""
from pathlib import Path
from itertools import combinations
import hashlib, importlib.util, json, resource, signal, time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / 'code/classification_plane_smooth_cubic_weights.py'
SPEC = importlib.util.spec_from_file_location('plane_weight_arithmetic', HELPER)
H = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(H)
TAU = s.Symbol('tau')
POLY = s.QQ.poly_ring(TAU)
FIELD = POLY.get_field()


def serial_matrix(mat):
    return [[str(s.factor(value)) for value in row] for row in mat.tolist()]


def field_dm(mat):
    return DomainMatrix.from_Matrix(mat).convert_to(POLY).convert_to(FIELD)


def rank_certificate(mat, variable=True):
    dm = field_dm(mat) if variable else DomainMatrix.from_Matrix(mat).to_field()
    rr, columns = dm.rref()
    rows = dm.transpose().rref()[1]
    rank = len(columns)
    assert len(rows) == rank
    reduced = rr.to_Matrix()
    free = [j for j in range(mat.cols) if j not in columns]
    kernel = s.zeros(mat.cols, len(free))
    for k, j in enumerate(free):
        kernel[j, k] = 1
        for i, pivot in enumerate(columns):
            kernel[pivot, k] = -reduced[i, j]
    assert all(s.cancel(value) == 0 for value in mat * kernel)
    minor = mat.extract(list(rows), list(columns))
    if variable:
        determinant = DomainMatrix.from_Matrix(minor).convert_to(POLY).det().as_expr()
    else:
        determinant = minor.det(method='domain-ge')
    assert determinant != 0 or rank == 0
    return {'rank': rank, 'nullity': len(free),
            'minor_rows': list(rows), 'minor_columns': list(columns),
            'minor_determinant': str(s.factor(determinant)),
            'complete_kernel': serial_matrix(kernel),
            'kernel_substitution_zero': True}, determinant, kernel


def operator(P0, P, source, x):
    jac = H.projective(H.bilinear(P0, P, x), x)
    rows = [s.sympify(value) for value in H.all350(jac, x)]
    matrix = s.Matrix([[s.expand(row).coeff(v) for v in source] for row in rows])
    assert all(s.expand(row - sum(matrix[i,j]*v for j,v in enumerate(source))) == 0
               for i,row in enumerate(rows))
    return matrix


def sources(x):
    z,w=x[3:]
    P1=s.zeros(5); variables=[]
    for name,normal in [('l',z),('n',w)]:
        part=s.zeros(5); cc=s.symbols(name+'0:9'); variables.extend(cc)
        for k,(i,j) in enumerate(combinations(range(3),2)):
            part[i,j]=sum(cc[3*k+r]*x[r] for r in range(3)); part[j,i]=-part[i,j]
        P1 += normal*part
    for name,normal,j in [('d',z*z,4),('b',z*w,3),('e',z*w,4),('c',w*w,3)]:
        cc=s.symbols(name+'0:3'); variables.extend(cc)
        for i in range(3):P1[i,j]+=normal*cc[i];P1[j,i]-=normal*cc[i]
    P2=s.zeros(5); cc=s.symbols('r0:9')
    for k,(i,j) in enumerate(combinations(range(3),2)):
        P2[i,j]=cc[3*k]*z*z+cc[3*k+1]*z*w+cc[3*k+2]*w*w;P2[j,i]=-P2[i,j]
    return P1,variables,P2,list(cc)


def shear(P0, x):
    cc=s.symbols('y0:6');z,w=x[3:]
    Y=s.Matrix([z*cc[i]+w*cc[i+3] for i in range(3)]+[0,0])
    orbit=s.zeros(5)
    for i,j in H.PAIRS:
        orbit[i,j]=s.expand(-sum(Y[k]*s.diff(P0[i,j],x[k]) for k in range(5))
                            +sum(P0[k,j]*s.diff(Y[i],x[k])+P0[i,k]*s.diff(Y[j],x[k])for k in range(5)))
        orbit[j,i]=-orbit[i,j]
    first=s.Matrix([s.Poly(orbit[i,3],z,w).coeff_monomial(z*z)for i in range(3)]+[0,0])
    second=s.Matrix([s.Poly(orbit[i,4],z,w).coeff_monomial(w*w)for i in range(3)]+[0,0])
    canonical=orbit+H.wedge(x,z*first+w*second)
    coeff=[]
    for normal in [z,w]:
        for i,j in combinations(range(3),2):
            linear=s.Poly(canonical[i,j],z,w).coeff_monomial(normal)
            coeff.extend([s.expand(linear).coeff(x[k])for k in range(3)])
    for normal,j in [(z*z,4),(z*w,3),(z*w,4),(w*w,3)]:
        coeff.extend([s.Poly(canonical[i,j],z,w).coeff_monomial(normal)for i in range(3)])
    mat=s.Matrix([[s.expand(value).coeff(c)for c in cc]for value in coeff])
    return mat


def main():
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite symbolic linear-operator budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    En=s.Matrix([0,0,0,z,w]); V=s.Matrix([0,0,0,b*z+c*w,a*z-b*w])
    ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1])
    P1,source1,P2,source2=sources(x)
    records=[]
    for name,U,T,K,offset in [
        ('semisimple',s.Matrix([0,1,0,0,0]),s.Matrix([2*TAU*a,0,-2*TAU*c,0,0]),(1-2*TAU)*z*w,17),
        ('nilpotent',s.Matrix([1,0,0,0,0]),s.Matrix([-2*TAU*b,TAU*c,0,0,0]),(TAU-s.Rational(1,2))*z*z,16)]:
        P0=H.wedge(T,En)+H.wedge(U,V)+K*H.wedge(ez,ew)
        assert all(s.expand(v)==0 for v in H.projective(H.bilinear(P0,P0,x),x).values())
        D1=operator(P0,P1,source1,x);D2=operator(P0,P2,source2,x)
        S=shear(P0,x);assert all(s.expand(v)==0 for v in D1*S)
        c1,m1,k1=rank_certificate(D1);c2,m2,k2=rank_certificate(D2);cs,ms,ks=rank_certificate(S)
        wall=s.Poly(m1*m2*ms,TAU)
        factors=s.factor_list(wall.as_expr())[1]
        walls=[]
        for factor,exponent in factors:
            polynomial=s.Poly(factor,TAU)
            if polynomial.degree()!=1:
                raise ValueError('Nonlinear wall needs exact algebraic-field certificate: '+str(factor))
            value=-polynomial.nth(0)/polynomial.nth(1)
            r1,_,_=rank_certificate(D1.subs(TAU,value),False)
            r2,_,_=rank_certificate(D2.subs(TAU,value),False)
            rs,_,_=rank_certificate(S.subs(TAU,value),False)
            higher=r1['nullity']+r2['nullity']
            walls.append({'minimal_polynomial':str(factor),'tau':str(value),'excluded_tau_zero':value==0,
                          'D1':r1,'D2':r2,'actual_shear':rs,'higher_fiber_upper_bound':higher,
                          'PGL_image_upper_bound':offset-1+higher,
                          'G_dimension_criterion':value==0 or offset-1+higher<=24})
        generic_higher=c1['nullity']+c2['nullity']
        generic_shear=(cs['rank']==6 and c1['nullity']==6 and c2['nullity']==0)
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])
                            for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
        records.append({'type':name,'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},
                        'P0_full_projective_Jacobi_zero':True,'actual_D':str(s.factor(actual)),
                        'complete_D1_350_by30':serial_matrix(D1),'complete_D2_350_by9':serial_matrix(D2),
                        'complete_actual_shear_30_by6':serial_matrix(S),
                        'generic_D1':c1,'generic_D2':c2,'generic_actual_shear':cs,
                        'exception_polynomial':str(s.factor(wall.as_expr())),'all_irreducible_factors':[[str(f),int(e)]for f,e in factors],
                        'all_exceptional_roots':walls,'generic_higher_fiber_upper_bound':generic_higher,
                        'generic_PGL_image_upper_bound':offset+generic_higher,
                        'generic_kernel_equals_actual_shear6_and_D2_injective':generic_shear,
                        'G_dimension_criterion':offset+generic_higher<=24 and all(r['G_dimension_criterion']for r in walls)})
        print(name,'generic nullities',c1['nullity'],c2['nullity'],'walls',[(r['tau'],r['higher_fiber_upper_bound'])for r in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md'
    result={'task_id':'G','all_passed':all(r['G_dimension_criterion']for r in records),
            'scope':'Entire mD2 rankzero-plane k1 coefficient-map rank3 stratum, only tau!=0. Other coefficient ranks,k0,m>=3 remain outside. Geometric normalization,source completeness,finite shears and group restoration are separate proof inputs.',
            'complete_source1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},'complete_source2':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},
            'cases':records,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,proof]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')




def rank2_kernel_main():
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite rank2-kernel operator budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old);after=functions(Path(__file__).read_text())
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    En=s.Matrix([0,0,0,z,w]);ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1]);U=s.Matrix([1,0,0,0,0])
    VB=s.Matrix([0,0,0,b*z+c*w,-b*w]);VN=s.Matrix([0,0,0,c*w,b*z])
    P1,source1,P2,source2=sources(x)
    models=[('ND',VN,s.Matrix([0,2*b,-2*c,0,0]),-2*z*w,16,False),
            ('BSS_rho',VB,s.Matrix([2*a+TAU*b,2*b,-4*c,0,0]),-6*z*w,16,True),
            ('BN1',VB,s.Matrix([c,0,2*b,0,0]),-w*w,15,False),
            ('BN0',VB,s.Matrix([0,0,2*b,0,0]),-w*w,14,False)]
    records=[]
    for name,V,T,K,offset,variable in models:
        P0=H.wedge(T,En)+H.wedge(U,V)+K*H.wedge(ez,ew)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,P0,x),x).values())
        D1=operator(P0,P1,source1,x);D2=operator(P0,P2,source2,x);S=shear(P0,x)
        assert all(s.expand(value)==0 for value in D1*S)
        c1,m1,k1=rank_certificate(D1,variable);c2,m2,k2=rank_certificate(D2,variable);cs,ms,ks=rank_certificate(S,variable)
        factors=s.factor_list(s.Poly(m1*m2*ms,TAU).as_expr())[1]if variable else []
        walls=[]
        for factor,exponent in factors:
            polynomial=s.Poly(factor,TAU)
            if polynomial.degree()!=1:raise ValueError('Nonlinear rho-wall needs exact algebraic-field certificate: '+str(factor))
            value=-polynomial.nth(0)/polynomial.nth(1)
            r1,_,_=rank_certificate(D1.subs(TAU,value),False);r2,_,_=rank_certificate(D2.subs(TAU,value),False);rs,_,_=rank_certificate(S.subs(TAU,value),False)
            higher=r1['nullity']+r2['nullity']
            walls.append({'minimal_polynomial':str(factor),'rho':str(value),'D1':r1,'D2':r2,'actual_shear':rs,
                          'higher_fiber_upper_bound':higher,'PGL_image_upper_bound':offset-1+higher,'G_dimension_criterion':offset-1+higher<=24})
        generic_higher=c1['nullity']+c2['nullity']
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
        record={'type':name,'rho_represented_by_tau':variable,'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},
                'P0_full_projective_Jacobi_zero':True,'actual_D':str(s.factor(actual)),
                'complete_D1_350_by30':serial_matrix(D1),'complete_D2_350_by9':serial_matrix(D2),'complete_actual_shear_30_by6':serial_matrix(S),
                'generic_D1':c1,'generic_D2':c2,'generic_actual_shear':cs,'generic_higher_fiber_upper_bound':generic_higher,
                'generic_PGL_image_upper_bound':offset+generic_higher,'exception_polynomial':str(s.factor(m1*m2*ms))if variable else None,
                'all_irreducible_factors':[[str(f),int(e)]for f,e in factors],'all_exceptional_roots':walls,
                'G_dimension_criterion':offset+generic_higher<=24 and all(row['G_dimension_criterion']for row in walls)}
        records.append(record)
        print(name,'nullities',c1['nullity'],c2['nullity'],'PGL upper',offset+generic_higher,'walls',[(row['rho'],row['higher_fiber_upper_bound'])for row in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_RANK_ZERO_RANK_TWO_KERNEL_FIRST_JET.md'
    result={'task_id':'G','all_passed':all(row['G_dimension_criterion']for row in records),
            'scope':'Entire mD2 rankzero-plane k1 rankM2/M(p)0 first-jet stratum: ND,allBSSrho,BN1,BN0. No conclusion if declared kernel-sum target fails; no nonlinear GB.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'complete_source1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},'complete_source2':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},
            'cases':records,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,proof,previous]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')



def multivariate_rank_certificate(mat, parameters):
    domain=s.QQ.poly_ring(*parameters) if parameters else s.QQ
    dm=DomainMatrix.from_Matrix(mat).convert_to(domain).to_field()
    rr,columns=dm.rref();rows=dm.transpose().rref()[1];rank=len(columns)
    reduced=rr.to_Matrix();free=[j for j in range(mat.cols)if j not in columns]
    kernel=s.zeros(mat.cols,len(free))
    for k,j in enumerate(free):
        kernel[j,k]=1
        for i,pivot in enumerate(columns):kernel[pivot,k]=-reduced[i,j]
    assert all(s.cancel(value)==0 for value in mat*kernel)
    minor=mat.extract(list(rows),list(columns))
    determinant=DomainMatrix.from_Matrix(minor).convert_to(domain).det()
    determinant=domain.to_sympy(determinant)
    assert determinant!=0 or rank==0
    return {'rank':rank,'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
            'minor_determinant':str(s.factor(determinant)),'complete_kernel':serial_matrix(kernel),'kernel_substitution_zero':True},determinant


def nonkernel_matrix_certificate(D1,D2,S,parameters,offset):
    c1,m1=multivariate_rank_certificate(D1,parameters);c2,m2=multivariate_rank_certificate(D2,parameters);cs,ms=multivariate_rank_certificate(S,parameters)
    higher=c1['nullity']+c2['nullity']
    return {'D1':c1,'D2':c2,'actual_shear':cs,'higher_fiber_upper_bound':higher,
            'PGL_image_upper_bound':offset+higher,'G_dimension_criterion':offset+higher<=24},s.factor(m1*m2*ms)


def algebraic_rank_certificate(mat,domain,alpha):
    values=[[domain.from_sympy(value.subs(TAU,alpha))for value in row]for row in mat.tolist()]
    dm=DomainMatrix(values,mat.shape,domain);rr,columns=dm.rref();rows=dm.transpose().rref()[1]
    reduced=rr.to_list();free=[j for j in range(mat.cols)if j not in columns]
    kernel=[[domain.zero for j in free]for i in range(mat.cols)]
    for k,j in enumerate(free):
        kernel[j][k]=domain.one
        for i,pivot in enumerate(columns):kernel[pivot][k]=-reduced[i][j]
    product=dm*DomainMatrix(kernel,(mat.cols,len(free)),domain)
    assert not any(any(row)for row in product.to_list())
    determinant=dm.extract(list(rows),list(columns)).det()
    assert determinant!=domain.zero or len(columns)==0
    return {'rank':len(columns),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
            'minor_determinant':str(domain.to_sympy(determinant)),
            'complete_kernel':[[str(domain.to_sympy(value))for value in row]for row in kernel],
            'kernel_substitution_zero':True}


def nonkernel_univariate(D1,D2,S,offset,excluded):
    generic,minor=nonkernel_matrix_certificate(D1,D2,S,(TAU,),offset)
    factors=s.factor_list(minor)[1];walls=[]
    for factor,multiplicity in factors:
        poly=s.Poly(factor,TAU)
        if poly.degree()==1:
            value=-poly.nth(0)/poly.nth(1)
            record,_=nonkernel_matrix_certificate(D1.subs(TAU,value),D2.subs(TAU,value),S.subs(TAU,value),(),offset-1)
            is_excluded=s.cancel(excluded.subs(TAU,value))==0
        else:
            value=s.CRootOf(poly.as_expr(),0);domain=s.QQ.algebraic_field(value)
            r1=algebraic_rank_certificate(D1,domain,value);r2=algebraic_rank_certificate(D2,domain,value);rs=algebraic_rank_certificate(S,domain,value)
            higher=r1['nullity']+r2['nullity']
            record={'D1':r1,'D2':r2,'actual_shear':rs,'higher_fiber_upper_bound':higher,
                    'PGL_image_upper_bound':offset-1+higher,'G_dimension_criterion':offset-1+higher<=24,
                    'exact_number_field':str(domain),'minimal_polynomial':str(factor),'all_conjugate_roots_covered':True}
            is_excluded=domain.from_sympy(excluded.subs(TAU,value))==domain.zero
        record.update({'parameter':str(value),'excluded_multiplicity_higher':is_excluded})
        if is_excluded:record['G_dimension_criterion']=True
        walls.append(record)
    generic.update({'exception_polynomial':str(minor),'all_irreducible_factors':[[str(f),int(e)]for f,e in factors],
                    'all_exceptional_roots':walls,'all_strata_dimension_criterion':generic['G_dimension_criterion']and all(w['G_dimension_criterion']for w in walls)})
    return generic


def rank2_nonkernel_main():
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite all-nonkernel operator budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes'];parse=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;ETA=s.Symbol('eta')
    En=s.Matrix([0,0,0,z,w]);ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1])
    ua=s.Matrix([1,0,0,0,0]);ub=s.Matrix([0,1,0,0,0]);uc=s.Matrix([0,0,1,0,0])
    VN=s.Matrix([0,0,0,c*w,b*z]);VB=s.Matrix([0,0,0,b*z+c*w,-b*w])
    diagonal=s.Matrix([2*a,-b,-c,0,0])
    P1,source1,P2,source2=sources(x)
    models=[('NDSS',ub+uc,VN,TAU*diagonal,-(z*z-w*w)/2,17,(TAU,),TAU),
            ('NDSSN',ub+uc,VN,s.Matrix([b-c,0,0,0,0]),-(z*z-w*w)/2,16,(),s.Integer(1)),
            ('NDNSS',ub,VN,diagonal,-z*z/2,16,(),s.Integer(1)),
            ('NDN',ub,VN,s.Matrix([c,0,0,0,0]),-z*z/2,15,(),s.Integer(1)),
            ('BSS',ub,VB,s.Matrix([(2*TAU-2*ETA)*a,ETA*b,(ETA-2*TAU)*c,0,0]),(1-2*TAU)*z*w,17,(TAU,ETA),4*TAU-3*ETA),
            ('BSSR',ub,VB,s.Matrix([-2*TAU*a/3+c,4*TAU*b/3,-2*TAU*c/3,0,0]),(1-2*TAU)*z*w,16,(TAU,),s.Integer(1)),
            ('BNS',uc,VB,s.Matrix([-2*a,b,c+2*TAU*b,0,0]),(s.Rational(1,2)-TAU)*w*w,16,(TAU,),s.Integer(1)),
            ('BNR',uc,VB,s.Matrix([b,0,2*TAU*b,0,0]),(s.Rational(1,2)-TAU)*w*w,15,(TAU,),s.Integer(1))]
    records=[]
    for name,U,V,T,K,offset,parameters,excluded in models:
        P0=H.wedge(T,En)+H.wedge(U,V)+K*H.wedge(ez,ew)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,P0,x),x).values())
        D1=operator(P0,P1,source1,x);D2=operator(P0,P2,source2,x);S=shear(P0,x)
        assert all(s.expand(value)==0 for value in D1*S)
        record={'type':name,'parameters':[str(v)for v in parameters],'excluded_polynomial':str(excluded),
                'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},
                'complete_D1_350_by30':serial_matrix(D1),'complete_D2_350_by9':serial_matrix(D2),'complete_actual_shear_30_by6':serial_matrix(S)}
        if len(parameters)==0:
            analysis,_=nonkernel_matrix_certificate(D1,D2,S,(),offset)
            analysis['all_strata_dimension_criterion']=analysis['G_dimension_criterion']
        elif len(parameters)==1:
            analysis=nonkernel_univariate(D1,D2,S,offset,excluded)
        else:
            analysis,minor=nonkernel_matrix_certificate(D1,D2,S,parameters,offset)
            factors=s.factor_list(minor)[1];lines=[]
            for factor,multiplicity in factors:
                poly=s.Poly(factor,TAU,ETA)
                if poly.total_degree()!=1:raise ValueError('Nonlinear twoparameter wall: '+str(factor))
                ca=poly.coeff_monomial(TAU);cb=poly.coeff_monomial(ETA);constant=poly.coeff_monomial(1)
                substitution={TAU:TAU,ETA:-(ca*TAU+constant)/cb}if cb!=0 else {TAU:-constant/ca,ETA:TAU}
                ex=s.expand(excluded.subs(substitution,simultaneous=True))
                line={'equation':str(factor),'parameterization':{str(k):str(v)for k,v in substitution.items()},'excluded_multiplicity_higher':ex==0}
                if ex==0:line['all_strata_dimension_criterion']=True
                else:
                    matrices=[M.subs(substitution,simultaneous=True)for M in [D1,D2,S]]
                    line['analysis']=nonkernel_univariate(*matrices,offset-1,ex)
                    line['all_strata_dimension_criterion']=line['analysis']['all_strata_dimension_criterion']
                lines.append(line)
            analysis.update({'exception_polynomial':str(minor),'all_irreducible_factors':[[str(f),int(e)]for f,e in factors],'all_exceptional_lines':lines,
                             'all_strata_dimension_criterion':analysis['G_dimension_criterion']and all(line['all_strata_dimension_criterion']for line in lines)})
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
        record.update({'analysis':analysis,'actual_D':str(s.factor(actual)),'P0_full_projective_Jacobi_zero':True})
        records.append(record)
        print(name,'generic higher',analysis['higher_fiber_upper_bound'],'genericPGL',analysis['PGL_image_upper_bound'],'allGcriterion',analysis['all_strata_dimension_criterion'],flush=True)
    proof=ROOT/'research/classification/PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_FIRST_JET.md'
    result={'task_id':'G','all_passed':all(r['analysis']['all_strata_dimension_criterion']for r in records),'cases':records,
            'scope':'Entire mD2/k1/rankM2/M(p)!=0 eight wholefirstjet families; excludesonlylistedhighermultiplicityloci; recordsanynonclosinglinearbound withoutnonlinearGB.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'complete_source1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},'complete_source2':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,proof,previous]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    try:
        if '--rank2-nonkernel' in __import__('sys').argv:rank2_nonkernel_main()
        elif '--rank2-kernel' in __import__('sys').argv:rank2_kernel_main()
        else:main()
    except FileExistsError:raise
    except Exception as error:
        failure=ROOT/('reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-nonkernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-kernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_RUNTIME_20261107.json')
        if not failure.exists():
            failure.write_text(json.dumps({'task_id':'G','all_passed':False,'failure':repr(error),
                'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'producer_bytes':Path(__file__).read_text(),'scope':'Stopped bounded symbolic linear-operator job; no mathematical G conclusion.'},indent=2)+'\n')
        raise
