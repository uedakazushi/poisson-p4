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

def global_lowparameter_main():
    """Whole at-most-one-parameter global-vector symbols; no anchor sampling."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole-lowparameter operator budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text())
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    En=s.Matrix([0,0,0,z,w]);ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1])
    Hbase=s.Matrix([a,b,-2*c,0,0]);N=s.Matrix([b,0,0,0,0]);B=s.Matrix([c,0,0,0,0]);C=s.Matrix([0,0,b,0,0])
    J=s.Matrix([b,c,0,0,0]);J2=s.Matrix([c,0,0,0,0]);G=s.Matrix([a,0,-c,0,0])
    I0=s.Matrix([c,0,0,0,0]);I1=s.Matrix([0,c,0,0,0])
    normal_ss=s.Matrix([0,0,0,z,-w]);normal_nil=s.Matrix([0,0,0,0,z])
    Drho=s.Matrix([(TAU+1)*a,TAU*b,(-1-2*TAU)*c,0,0])
    # tuple: whole source, normal type, T,V,K, affine-symbol upper, varying parameter, excluded polynomial
    models=[
        ('SS_NILSQUARE','semisimple',J,J2,TAU*z*w,10,True,s.Integer(1)),
        ('SS_COMMON_K1','semisimple',I0,I1,z*w,9,False,s.Integer(1)),
        ('SS_COMMON_K0','semisimple',I0,I1,s.Integer(0),8,False,s.Integer(1)),
        ('NIL_B1_REGULAR','nilpotent',G,J,z*w,11,False,s.Integer(1)),
        ('NIL_B1_DIAGONAL','nilpotent',Drho,N,z*w,11,True,1+3*TAU),
        ('NIL_B1_RESONANT_B','nilpotent',Drho.subs(TAU,-s.Rational(2,3))+B,N,z*w,10,False,s.Integer(1)),
        ('NIL_B1_RESONANT_C','nilpotent',Drho.subs(TAU,-s.Rational(1,3))+C,N,z*w,10,False,s.Integer(1)),
        ('NIL_B0_NILREG','nilpotent',J+J2,J,TAU*z*z,10,True,s.Integer(1)),
        ('NIL_B0_TSQUARE_K1','nilpotent',J2,J,z*z,9,False,s.Integer(1)),
        ('NIL_B0_TSQUARE_K0','nilpotent',J2,J,s.Integer(0),8,False,s.Integer(1)),
        ('NIL_B0_VSQUARE_K1','nilpotent',J,J2,z*z,9,False,s.Integer(1)),
        ('NIL_B0_VSQUARE_K0','nilpotent',J,J2,s.Integer(0),8,False,s.Integer(1)),
        ('NIL_B0_COMMON_K1','nilpotent',I0,I1,z*z,8,False,s.Integer(1)),
        ('NIL_B0_COMMON_K0','nilpotent',I0,I1,s.Integer(0),7,False,s.Integer(1)),
        ('NIL_B0_MIX_TN','nilpotent',N,Hbase+N,TAU*z*z,10,True,s.Integer(1)),
        ('NIL_B0_MIX_VH','nilpotent',Hbase+N,Hbase,TAU*z*z,10,True,s.Integer(1)),
        ('NIL_B0_MIX_TN_VH','nilpotent',N,Hbase,TAU*z*z,9,True,s.Integer(1)),
        ('NIL_B0_MIX_VN','nilpotent',Hbase+N,N,TAU*z*z,10,True,s.Integer(1)),
        ('NIL_B0_MIX_TH','nilpotent',Hbase,Hbase+N,TAU*z*z,10,True,s.Integer(1)),
        ('NIL_B0_MIX_TH_VN_K1','nilpotent',Hbase,N,z*z,9,False,s.Integer(1)),
        ('NIL_B0_MIX_TH_VN_K0','nilpotent',Hbase,N,s.Integer(0),8,False,s.Integer(1))]
    P1,source1,P2,source2=sources(x);records=[]
    for name,normal_type,T,V,K,base_dimension,variable,excluded in models:
        normal=normal_ss if normal_type=='semisimple' else normal_nil
        P0=H.wedge(T,En)+H.wedge(V,normal)+K*H.wedge(ez,ew)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,P0,x),x).values())
        determinant=s.expand(s.det(s.Matrix.hstack(x[:3,0],T[:3,0],V[:3,0])))
        assert determinant!=0
        D1=operator(P0,P1,source1,x);D2=operator(P0,P2,source2,x);S=shear(P0,x)
        assert all(s.expand(value)==0 for value in D1*S)
        if variable:analysis=nonkernel_univariate(D1,D2,S,base_dimension+5,excluded)
        else:
            analysis,_=nonkernel_matrix_certificate(D1,D2,S,(),base_dimension+5)
            analysis['all_strata_dimension_criterion']=analysis['G_dimension_criterion']
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
        record={'type':name,'normal_type':normal_type,'varying_parameter':str(TAU)if variable else None,
                'affine_first_symbol_dimension_upper':base_dimension,'excluded_polynomial':str(excluded),
                'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},
                'P0_full_projective_Jacobi_zero':True,'base_determinant':str(s.factor(determinant)),
                'actual_D':str(s.factor(actual)),
                'complete_D1_350_by30':serial_matrix(D1),'complete_D2_350_by9':serial_matrix(D2),
                'complete_actual_shear_30_by6':serial_matrix(S),'analysis':analysis}
        records.append(record)
        print(name,'higher',analysis['higher_fiber_upper_bound'],'image',analysis['PGL_image_upper_bound'],
              'allGcriterion',analysis['all_strata_dimension_criterion'],flush=True)
    proofs=[ROOT/'research/classification/PLANE_RANK_ZERO_GLOBAL_VECTOR_SEMISIMPLE_FIRST_JETS.md',
            ROOT/'research/classification/PLANE_RANK_ZERO_GLOBAL_VECTOR_NILPOTENT_FIRST_JETS.md']
    result={'task_id':'G','all_passed':all(record['analysis']['all_strata_dimension_criterion']for record in records),
            'scope':'Twentyone complete at-most-one-parameter mD2 globalV first-symbol strata. All parameters and all exact irreducible maximal-minor root complements retained. Generic two/three-parameter mixed/Cartan/nil sources remain outside this finite packet.',
            'cases':records,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'complete_source1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},
            'complete_source2':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,previous]+proofs],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_onecut_main():
    """One genuine higher-Jacobi cut suffices for each of two complete K=0 fibers."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ONECUT_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite two-fiber onecut budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json'
    packet=json.loads(previous.read_text());old=packet['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));environment={str(value):value for value in x};environment['tau']=TAU
    P1,source1,P2,source2=sources(x);records=[]
    for name in ['SS_NILSQUARE','NIL_B0_NILREG']:
        original=next(case for case in packet['cases']if case['type']==name)
        wall=next(row for row in original['analysis']['all_exceptional_roots']if row['parameter']=='0')
        K1=s.Matrix([[s.sympify(value)for value in row]for row in wall['D1']['complete_kernel']])
        K2=s.Matrix([[s.sympify(value)for value in row]for row in wall['D2']['complete_kernel']])
        assert K1.cols==8 and K2.cols==3
        D2=s.Matrix([[s.sympify(value,locals=environment).subs(TAU,0)for value in row]for row in original['complete_D2_350_by9']])
        P0=s.zeros(5)
        for i,j in H.PAIRS:
            P0[i,j]=s.sympify(original['normal_form_P0'][str((i,j))],locals=environment).subs(TAU,0);P0[j,i]=-P0[i,j]
        u=s.symbols('u0:8');v=s.symbols('v0:3')
        first=P1.subs(dict(zip(source1,K1*s.Matrix(u))),simultaneous=True)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,first,x),x).values())
        Q=s.Matrix([s.sympify(value)/2 for value in H.all350(H.projective(H.bilinear(first,first,x),x),x)])
        rows=wall['D2']['minor_rows'];columns=wall['D2']['minor_columns'];minor=D2.extract(rows,columns)
        assert len(rows)==len(columns)==6 and minor.det()!=0
        inverse=minor.inv();particular=-inverse*Q.extract(rows,[0])
        residual=Q+D2[:,columns]*particular
        nonzero=next((i for i,value in enumerate(residual)if s.expand(value)!=0),None)
        record={'type':name,'full_weight1_kernel':serial_matrix(K1),'full_weight2_kernel':serial_matrix(K2),
                'weight2_pivot_rows':rows,'weight2_pivot_columns':columns,'weight2_pivot_inverse':serial_matrix(inverse),
                'full_Q2_350':[str(s.expand(value))for value in Q],
                'affine_first_symbol_dimension':9,'initial_higher_linear_upper':11}
        if nonzero is not None:
            phi=s.zeros(1,350);phi[0,nonzero]=1
            coefficients=D2[nonzero:nonzero+1,columns]*inverse
            for j,row in enumerate(rows):phi[0,row]-=coefficients[j]
            assert phi*D2==s.zeros(1,9)
            polynomial=s.expand((phi*Q)[0]);assert polynomial!=0
            record.update({'certificate_type':'nonzero_weight2_cokernel_obstruction',
                           'cokernel_functional':[[i,str(value)]for i,value in enumerate(list(phi))if value!=0],
                           'functional_times_full_D2_zero':True,'obstruction_polynomial':str(s.factor(polynomial)),
                           'full_weight1_zero_locus_dimension_upper':7,'weight2_fiber_dimension':3,
                           'full_higher_dimension_upper':10,'PGL_image_upper':24,'G_dimension_criterion':True})
        else:
            coefficients=s.zeros(9,1)
            for i,column in enumerate(columns):coefficients[column]=particular[i]
            coefficients+=K2*s.Matrix(v)
            second=P2.subs(dict(zip(source2,coefficients)),simultaneous=True)
            assert all(s.expand(value)==0 for value in D2*coefficients+Q)
            J3=[s.sympify(value)for value in H.all350(H.projective(H.bilinear(first,second,x),x),x)]
            J4=[s.sympify(value)/2 for value in H.all350(H.projective(H.bilinear(second,second,x),x),x)]
            witness=next(((weight,i,s.expand(value))for weight,values in [(3,J3),(4,J4)]for i,value in enumerate(values)if s.expand(value)!=0),None)
            record.update({'full_weight2_parametrization':[str(s.expand(value))for value in coefficients],
                           'weight2_equations_identically_solved':True,'higher_weight3_350':[str(s.expand(value))for value in J3],
                           'higher_weight4_350':[str(s.expand(value))for value in J4]})
            if witness is not None:
                weight,index,polynomial=witness
                record.update({'certificate_type':'nonzero_higher_Jacobi_on_full_weight2_solution',
                               'witness_weight':weight,'witness_row':index,'obstruction_polynomial':str(s.factor(polynomial)),
                               'full_higher_dimension_upper':10,'PGL_image_upper':24,'G_dimension_criterion':True})
            else:
                record.update({'certificate_type':'all_higher_Jacobi_zero_on_full_11parameter_space',
                               'full_higher_dimension_upper':11,'PGL_image_upper':25,'G_dimension_criterion':False})
        records.append(record)
        print(name,record['certificate_type'],'image',record['PGL_image_upper'],flush=True)
    result={'task_id':'G','all_passed':all(case['G_dimension_criterion']for case in records),
            'scope':'Two entire K0 higher fibers, with base9 each; a single genuine Jacobi cut suffices. No nonlinear Groebner basis or finite sample membership is used.',
            'cases':records,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_compatibility_main():
    """Three small whole shear slices, using only necessary weight2 quadrics."""
    import ast,os,subprocess,tempfile
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMPATIBILITY_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite three compatibility cores budget210s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(210)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json'
    oldpath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ONECUT_20261107.json'
    packet=json.loads(previous.read_text());old=json.loads(oldpath.read_text())['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    prefix=ROOT/'work/nonsquare-v62/singular-runtime/prefix';native=prefix/'usr/bin/Singular'
    nativeenv=dict(os.environ)
    nativeenv['LD_LIBRARY_PATH']='/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu:'+str(prefix/'usr/lib/x86_64-linux-gnu')+':'+nativeenv.get('LD_LIBRARY_PATH','')
    nativeenv['SINGULAR_PROCS_DIR']=str(prefix/'usr/lib/x86_64-linux-gnu/singular/MOD')
    nativeenv['SINGULARPATH']='/workspace/cloud-setup/prefix/usr/share/singular/LIB'
    x=s.Matrix(s.symbols('a b c z w'));environment={str(value):value for value in x};environment['tau']=TAU
    P1,source1,_,_=sources(x);records=[]
    for name in ['SS_COMMON_K0','NIL_B0_TSQUARE_K0','NIL_B0_COMMON_K0']:
        original=next(case for case in packet['cases']if case['type']==name);analysis=original['analysis']
        K1=s.Matrix([[s.sympify(value)for value in row]for row in analysis['D1']['complete_kernel']])
        D1=s.Matrix([[s.sympify(value)for value in row]for row in original['complete_D1_350_by30']])
        D2=s.Matrix([[s.sympify(value)for value in row]for row in original['complete_D2_350_by9']])
        S=s.Matrix([[s.sympify(value)for value in row]for row in original['complete_actual_shear_30_by6']])
        shearpivots=analysis['actual_shear']['minor_columns'];B=S[:,shearpivots];r=len(shearpivots)
        combined=B.row_join(K1);pivots=list(DomainMatrix.from_Matrix(combined).to_field().rref()[1])
        assert pivots[:r]==list(range(r)) and len(pivots)==K1.cols
        complement=[column-r for column in pivots[r:]];R=K1[:,complement];basis=B.row_join(R)
        assert D1*basis==s.zeros(350,K1.cols)
        basis_rows=list(DomainMatrix.from_Matrix(basis.T).to_field().rref()[1]);basis_minor=basis[basis_rows,:]
        assert basis_minor.det()!=0
        u=s.symbols('u0:'+str(R.cols));uenv={str(value):value for value in u}
        first=P1.subs(dict(zip(source1,R*s.Matrix(u))),simultaneous=True)
        Q=s.Matrix([s.sympify(value)/2 for value in H.all350(H.projective(H.bilinear(first,first,x),x),x)])
        rows=analysis['D2']['minor_rows'];columns=analysis['D2']['minor_columns'];inverse=D2.extract(rows,columns).inv()
        residual=Q-D2[:,columns]*inverse*Q.extract(rows,[0])
        monomials=[s.prod(u[i]**exponent for i,exponent in enumerate(powers))for powers in H.monomials(len(u),2)]
        coefficient_matrix=s.Matrix([[s.expand(value).coeff(monomial)for monomial in monomials]for value in residual])
        assert all(s.expand(value-sum(coefficient_matrix[i,j]*m for j,m in enumerate(monomials)))==0 for i,value in enumerate(residual))
        selected=list(DomainMatrix.from_Matrix(coefficient_matrix.T).to_field().rref()[1])
        generators=[s.Poly(s.expand(residual[row]),*u,domain=s.QQ)for row in selected]
        functionals=[]
        for row in selected:
            phi=s.zeros(1,350);phi[0,row]=1;coefficients=D2[row:row+1,columns]*inverse
            for j,pivotrow in enumerate(rows):phi[0,pivotrow]-=coefficients[j]
            assert phi*D2==s.zeros(1,9) and s.expand((phi*Q)[0]-residual[row])==0
            functionals.append([[i,str(value)]for i,value in enumerate(list(phi))if value!=0])
        record={'type':name,'affine_first_symbol_dimension_upper':original['affine_first_symbol_dimension_upper'],
                'full_first_weight_kernel':serial_matrix(K1),'actual_shear_rank':r,'actual_shear_selected_columns':shearpivots,
                'complete_transverse_weight1_30_by_n':serial_matrix(R),'shear_plus_slice_basis_rows':basis_rows,
                'shear_plus_slice_basis_minor':serial_matrix(basis_minor),'shear_plus_slice_basis_inverse':serial_matrix(basis_minor.inv()),
                'slice_variables':list(map(str,u)),'full_Q2_350':[str(s.expand(value))for value in Q],
                'D2_pivot_rows':rows,'D2_pivot_columns':columns,'D2_pivot_inverse':serial_matrix(inverse),
                'selected_obstruction_rows':selected,'obstruction_cokernel_functionals':functionals,
                'necessary_compatibility_generators':[str(value.as_expr())for value in generators],
                'weight2_fiber_dimension':analysis['D2']['nullity']}
        if not generators:
            higher=r+len(u)+analysis['D2']['nullity'];ambient=original['affine_first_symbol_dimension_upper']+higher+5
            record.update(status='all_weight2_compatibility_zero_no_improvement',compatibility_slice_dimension_upper=len(u),
                          full_higher_dimension_upper=higher,PGL_image_upper=ambient,G_dimension_criterion=ambient<=24)
            records.append(record)
            print(name,record['status'],'slice',len(u),'image',ambient,flush=True)
            continue
        script='ring R=0,('+','.join(map(str,u))+'),dp;\nshort=0;\n'
        script+='ideal I='+','.join(str(value.as_expr()).replace('**','^')for value in generators)+';\nmatrix T;ideal G=liftstd(I,T);\n'
        script+='print("DIM:"+string(dim(G)));\nint i;int j;for(i=1;i<=size(G);i++){print("G:"+string(i)+":"+string(G[i]));}\n'
        script+='for(i=1;i<=nrows(T);i++){for(j=1;j<=ncols(T);j++){if(T[i,j]!=0){print("L:"+string(i)+":"+string(j)+":"+string(T[i,j]));}}}\nquit;\n'
        record['exact_native_input']=script
        try:
            with tempfile.TemporaryDirectory(prefix='plane-compatibility-')as temp:
                task=Path(temp)/'source.sing';task.write_text(script)
                run=subprocess.run([str(native),'-q',str(task)],env=nativeenv,capture_output=True,text=True,timeout=45,check=False)
            if run.returncode or '?'in run.stdout:
                record.update(status='native_failure_no_dimension_proof',returncode=run.returncode,stdout=run.stdout[-20000:],stderr=run.stderr[-20000:],G_dimension_criterion=False)
            else:
                lifted={};lifts=[];native_dim=None
                for line in run.stdout.splitlines():
                    if line.startswith('G:'):
                        _,j,polynomial=line.split(':',2);lifted[int(j)-1]=s.Poly(s.sympify(polynomial.replace('^','**'),locals=uenv),*u,domain=s.QQ)
                    elif line.startswith('L:'):
                        _,i,j,polynomial=line.split(':',3);lifts.append((int(i)-1,int(j)-1,s.Poly(s.sympify(polynomial.replace('^','**'),locals=uenv),*u,domain=s.QQ)))
                    elif line.startswith('DIM:'):native_dim=int(line.split(':',1)[1])
                assert lifted and native_dim is not None
                sums={j:s.Poly(0,*u,domain=s.QQ)for j in lifted}
                for i,j,value in lifts:sums[j]+=generators[i]*value
                assert lifted==sums
                leading=[lifted[j].monoms(order='grevlex')[0]for j in sorted(lifted)]
                edges=[sum(1<<i for i,exponent in enumerate(powers)if exponent)for powers in leading]
                cover=next(indices for k in range(len(u)+1)for indices in combinations(range(len(u)),k)if all(edge&sum(1<<i for i in indices)for edge in edges))
                height=len(cover);bound=len(u)-height;assert bound==native_dim
                counterexamples=[]
                if height:
                    for indices in combinations(range(len(u)),height-1):
                        mask=sum(1<<i for i in indices);counterexamples.append(next(i for i,edge in enumerate(edges)if not edge&mask))
                higher=r+bound+analysis['D2']['nullity'];ambient=original['affine_first_symbol_dimension_upper']+higher+5
                record.update(status='exact_lifted_necessary_compatibility_dimension',all_ideal_lifts_verified=True,
                    lifted_polynomials=[str(lifted[j].as_expr())for j in sorted(lifted)],
                    ideal_membership_lifts=[{'generator':i,'polynomial':j,'coefficient':str(value.as_expr())}for i,j,value in lifts],
                    leading_monomial_exponents=leading,leading_support_minimum_cover=height,minimum_cover_variables=list(cover),
                    all_size_one_less_cover_counterexample_indices=counterexamples,
                    compatibility_slice_dimension_upper=bound,full_higher_dimension_upper=higher,
                    PGL_image_upper=ambient,G_dimension_criterion=ambient<=24,native_stdout=run.stdout)
        except subprocess.TimeoutExpired as error:
            record.update(status='finite45second_native_limit_no_dimension_proof',stdout=str(error.stdout or '')[-20000:],stderr=str(error.stderr or '')[-20000:],G_dimension_criterion=False)
        records.append(record)
        print(name,record['status'],'slice',record.get('compatibility_slice_dimension_upper'),'image',record.get('PGL_image_upper'),flush=True)
    result={'task_id':'G','all_passed':all(case['G_dimension_criterion']for case in records),
            'scope':'Three entire fixed K0 globalV symbol fibers. Actual finite shears cover full first-weight source; only full-weight2 necessary compatibility quadrics on8/6/8variable transverse slices are used. No full Poisson ideal or80variable GB.',
            'cases':records,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous,oldpath,native]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_residual_main():
    """Complete original higher Jacobi on seven whole linear compatibility branches."""
    import ast,os,subprocess,tempfile
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_RESIDUAL_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite original-Jacobi residual budget210s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(210)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    prior=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMPATIBILITY_20261107.json'
    sourcepath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json'
    cores=json.loads(prior.read_text());packet=json.loads(sourcepath.read_text());old=cores['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    records=[];x=s.Matrix(s.symbols('a b c z w'));env={str(value):value for value in x};env['tau']=TAU
    template1,vars1,template2,vars2=sources(x)
    for name in ['SS_COMMON_K0','NIL_B0_COMMON_K0']:
        core=next(case for case in cores['cases']if case['type']==name)
        original=next(case for case in packet['cases']if case['type']==name)
        u=s.symbols('u0:8');v=s.symbols('v0:3');R=s.Matrix([[s.sympify(value)for value in row]for row in core['complete_transverse_weight1_30_by_n']])
        K2=s.Matrix([[s.sympify(value)for value in row]for row in original['analysis']['D2']['complete_kernel']])
        D2=s.Matrix([[s.sympify(value)for value in row]for row in original['complete_D2_350_by9']])
        Q=s.Matrix([s.sympify(value)for value in core['full_Q2_350']]);inverse=s.Matrix([[s.sympify(value)for value in row]for row in core['D2_pivot_inverse']])
        particular=-inverse*Q.extract(core['D2_pivot_rows'],[0]);coefficient=s.zeros(9,1)
        for i,column in enumerate(core['D2_pivot_columns']):coefficient[column]=particular[i]
        coefficient+=K2*s.Matrix(v)
        first=template1.subs(dict(zip(vars1,R*s.Matrix(u))),simultaneous=True)
        second=template2.subs(dict(zip(vars2,coefficient)),simultaneous=True)
        if name=='SS_COMMON_K0':
            branches=[{u[7]:0,u[4]:0},{u[7]:0,u[1]:2*u[0]+2*u[6]},
                      {u[3]:-2*u[2]-2*u[5],u[4]:0},{u[3]:-2*u[2]-2*u[5],u[1]:2*u[0]+2*u[6]}]
            target_height=2
        else:
            branches=[{u[7]:0,u[3]:0},{u[7]:0,u[5]:0},{u[3]:-u[7],u[1]:2*u[2]+4*u[5]}]
            target_height=1
        for number,substitution in enumerate(branches):
            free=tuple(value for value in u if value not in substitution)+tuple(v);assert len(free)==9
            first_branch=first.subs(substitution,simultaneous=True);second_branch=second.subs(substitution,simultaneous=True)
            assert all(s.expand(value)==0 for value in (D2*coefficient+Q).subs(substitution,simultaneous=True))
            J3=[s.sympify(value)for value in H.all350(H.projective(H.bilinear(first_branch,second_branch,x),x),x)]
            J4=[s.sympify(value)/2 for value in H.all350(H.projective(H.bilinear(second_branch,second_branch,x),x),x)]
            nonzero=[(weight,i,s.Poly(value,*free,domain=s.QQ))for weight,values in [(3,J3),(4,J4)]for i,value in enumerate(values)if value!=0]
            witnesses=[];gcd=None
            for weight,i,polynomial in nonzero:
                new_gcd=polynomial.monic()if gcd is None else s.polys.polytools.gcd(gcd,polynomial).monic()
                if gcd is None or new_gcd!=gcd:
                    witnesses.append({'weight':weight,'row':i,'polynomial':str(polynomial.as_expr()),'cumulative_gcd':str(new_gcd.as_expr())})
                    gcd=new_gcd
                if target_height==1 or gcd.total_degree()==0:break
            height=0 if gcd is None else 1 if gcd.total_degree()!=0 else 2
            higher=core['actual_shear_rank']+9-height;ambient=original['affine_first_symbol_dimension_upper']+higher+5
            record={'type':name,'linear_branch':number,'substitution':{str(key):str(value)for key,value in substitution.items()},
                    'all_weight2_solutions_parametrized':True,'free_variables':list(map(str,free)),
                    'complete_weight1_coefficients':{str((i,j)):str(s.expand(first_branch[i,j]))for i,j in H.PAIRS},
                    'complete_weight2_coefficients':{str((i,j)):str(s.expand(second_branch[i,j]))for i,j in H.PAIRS},
                    'all350_weight3_Jacobi':[str(s.expand(value))for value in J3],
                    'all350_weight4_Jacobi':[str(s.expand(value))for value in J4],
                    'exact_nonzero_polynomial_or_gcd_chain':witnesses,'proven_height_lower_bound':height,
                    'dimension_argument':'A nonzero ideal has height>=1; gcd1 in the polynomial UFD excludes every height1 prime, so height>=2, not unit ideal.',
                    'normalized_full_higher_dimension_upper':9-height,'actual_shear_rank':core['actual_shear_rank'],
                    'full_higher_dimension_upper':higher,'PGL_image_upper':ambient,'G_dimension_criterion':ambient<=24}
            records.append(record)
            print(name,'branch',number,'height',height,'image',ambient,flush=True)
    # Same exact six-variable core as before; fix native rational-power serialization only.
    core=next(case for case in cores['cases']if case['type']=='NIL_B0_TSQUARE_K0')
    u=s.symbols('u0:6');uenv={str(value):value for value in u}
    generators=[s.Poly(s.sympify(value),*u,domain=s.QQ)for value in core['necessary_compatibility_generators']]
    def native_polynomial(polynomial):
        terms=[]
        for powers,coefficient in polynomial.terms():
            monomial='*'.join(str(variable)+'^'+str(exponent)for variable,exponent in zip(u,powers)if exponent)
            terms.append('('+str(coefficient)+')'+('*'+monomial if monomial else ''))
        return '+'.join(terms)or'0'
    script='ring R=0,('+','.join(map(str,u))+'),dp;\nshort=0;\n'
    script+='ideal I='+','.join(native_polynomial(value)for value in generators)+';\nmatrix T;ideal G=liftstd(I,T);\n'
    script+='print("DIM:"+string(dim(G)));\nint i;int j;for(i=1;i<=size(G);i++){print("G:"+string(i)+":"+string(G[i]));}\n'
    script+='for(i=1;i<=nrows(T);i++){for(j=1;j<=ncols(T);j++){if(T[i,j]!=0){print("L:"+string(i)+":"+string(j)+":"+string(T[i,j]));}}}\nquit;\n'
    # Independently parse the corrected polynomial strings back to QQ before running native.
    assert [s.Poly(s.sympify(native_polynomial(value).replace('^','**'),locals=uenv),*u,domain=s.QQ)for value in generators]==generators
    prefix=ROOT/'work/nonsquare-v62/singular-runtime/prefix';native=prefix/'usr/bin/Singular';nativeenv=dict(os.environ)
    nativeenv['LD_LIBRARY_PATH']='/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu:'+str(prefix/'usr/lib/x86_64-linux-gnu')+':'+nativeenv.get('LD_LIBRARY_PATH','')
    nativeenv['SINGULAR_PROCS_DIR']=str(prefix/'usr/lib/x86_64-linux-gnu/singular/MOD');nativeenv['SINGULARPATH']='/workspace/cloud-setup/prefix/usr/share/singular/LIB'
    repaired={'type':'NIL_B0_TSQUARE_K0','same_original_core_all_coefficients_exact':True,'exact_native_input':script,
              'necessary_compatibility_generators':core['necessary_compatibility_generators']}
    try:
        with tempfile.TemporaryDirectory(prefix='plane-compatibility-repair-')as temp:
            task=Path(temp)/'source.sing';task.write_text(script)
            run=subprocess.run([str(native),'-q',str(task)],env=nativeenv,capture_output=True,text=True,timeout=45,check=False)
        if run.returncode or '?'in run.stdout:
            repaired.update(status='native_failure_no_dimension_proof',stdout=run.stdout[-20000:],stderr=run.stderr[-20000:],G_dimension_criterion=False)
        else:
            lifted={};lifts=[]
            for line in run.stdout.splitlines():
                if line.startswith('G:'):
                    _,j,value=line.split(':',2);lifted[int(j)-1]=s.Poly(s.sympify(value.replace('^','**'),locals=uenv),*u,domain=s.QQ)
                elif line.startswith('L:'):
                    _,i,j,value=line.split(':',3);lifts.append((int(i)-1,int(j)-1,s.Poly(s.sympify(value.replace('^','**'),locals=uenv),*u,domain=s.QQ)))
            sums={j:s.Poly(0,*u,domain=s.QQ)for j in lifted}
            for i,j,value in lifts:sums[j]+=generators[i]*value
            assert lifted and lifted==sums
            leading=[lifted[j].monoms(order='grevlex')[0]for j in sorted(lifted)]
            edges=[sum(1<<i for i,exponent in enumerate(powers)if exponent)for powers in leading]
            cover=next(indices for k in range(7)for indices in combinations(range(6),k)if all(edge&sum(1<<i for i in indices)for edge in edges))
            height=len(cover);bound=6-height;higher=5+bound+3;ambient=8+higher+5
            repaired.update(status='exact_lifted_same_core_dimension',all_ideal_lifts_verified=True,
                lifted_polynomials=[str(lifted[j].as_expr())for j in sorted(lifted)],
                ideal_membership_lifts=[{'generator':i,'polynomial':j,'coefficient':str(value.as_expr())}for i,j,value in lifts],
                leading_monomial_exponents=leading,leading_support_minimum_cover=height,minimum_cover_variables=list(cover),
                compatibility_slice_dimension_upper=bound,full_higher_dimension_upper=higher,PGL_image_upper=ambient,
                G_dimension_criterion=ambient<=24,native_stdout=run.stdout)
    except subprocess.TimeoutExpired as error:
        repaired.update(status='finite45second_native_limit_no_dimension_proof',stdout=str(error.stdout or '')[-20000:],G_dimension_criterion=False)
    records.append(repaired);print('NIL_B0_TSQUARE_K0 repair',repaired['status'],'image',repaired.get('PGL_image_upper'),flush=True)
    result={'task_id':'G','all_passed':all(case['G_dimension_criterion']for case in records),'cases':records,
            'scope':'Complete original weights3/4 on all four SS and three nil common compatibility branches, each with all three second-weight kernel variables, plus exact same TSQUARE compatibility core serialization repair. No new jet or changed source.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,prior,sourcepath,native]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_tsquare_main():
    """The last compatibility support is one linear space, not a sampled chart."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TSQUARE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole TSQUARE original-Jacobi budget60s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    corepath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMPATIBILITY_20261107.json'
    sourcepath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json'
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_RESIDUAL_20261107.json'
    cores=json.loads(corepath.read_text());packet=json.loads(sourcepath.read_text());old=json.loads(previous.read_text())['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    core=next(case for case in cores['cases']if case['type']=='NIL_B0_TSQUARE_K0')
    original=next(case for case in packet['cases']if case['type']=='NIL_B0_TSQUARE_K0')
    u=s.symbols('u0:6');v=s.symbols('v0:3');L=-u[2]+3*u[3]+u[4];difference=u[2]-u[4]
    generators=[s.sympify(value)for value in core['necessary_compatibility_generators']]
    assert s.expand(generators[1]+12*u[5]**2)==0
    assert s.expand(generators[3]-u[3]*L/3)==0
    assert s.expand(generators[5]-14*u[0]*u[5]/3-L*(2*difference+9*u[3])/6)==0
    # At any reduced point u5=0; then g3=0 and g5=0 imply L^2=0, hence L=0.
    assert s.expand(L**2-((s.Rational(15,2)*u[3]*L)-(s.Rational(1,2)*L*(2*difference+9*u[3]))))==0
    substitution={u[5]:0,u[2]:3*u[3]+u[4]}
    assert all(s.expand(value.subs(substitution,simultaneous=True))==0 for value in generators)
    x=s.Matrix(s.symbols('a b c z w'));P1,source1,P2,source2=sources(x)
    R=s.Matrix([[s.sympify(value)for value in row]for row in core['complete_transverse_weight1_30_by_n']])
    K2=s.Matrix([[s.sympify(value)for value in row]for row in original['analysis']['D2']['complete_kernel']])
    D2=s.Matrix([[s.sympify(value)for value in row]for row in original['complete_D2_350_by9']])
    Q=s.Matrix([s.sympify(value)for value in core['full_Q2_350']]);inverse=s.Matrix([[s.sympify(value)for value in row]for row in core['D2_pivot_inverse']])
    particular=-inverse*Q.extract(core['D2_pivot_rows'],[0]);coefficient=s.zeros(9,1)
    for i,column in enumerate(core['D2_pivot_columns']):coefficient[column]=particular[i]
    coefficient+=K2*s.Matrix(v)
    assert all(s.expand(value)==0 for value in (D2*coefficient+Q).subs(substitution,simultaneous=True))
    first=P1.subs(dict(zip(source1,R*s.Matrix(u))),simultaneous=True).subs(substitution,simultaneous=True)
    second=P2.subs(dict(zip(source2,coefficient)),simultaneous=True).subs(substitution,simultaneous=True)
    J3=[s.sympify(value)for value in H.all350(H.projective(H.bilinear(first,second,x),x),x)]
    J4=[s.sympify(value)/2 for value in H.all350(H.projective(H.bilinear(second,second,x),x),x)]
    witness=next(((weight,i,s.expand(value))for weight,values in [(3,J3),(4,J4)]for i,value in enumerate(values)if value!=0),None)
    height=1 if witness is not None else 0;ambient=8+5+7-height+5
    result={'task_id':'G','all_passed':ambient<=24,'type':'NIL_B0_TSQUARE_K0',
            'scope':'Entire reduced compatibility support, one whole linear4space, all3 secondweight kernel variables; remaining original Jacobi gives one nonzero cut if present.',
            'reduced_support_equal_linear_space':True,'support_linear_equations':['u5=0',str(L)+'=0'],
            'support_identities':['g1=-12*u5^2','g3=u3*L/3','g5=14*u0*u5/3+L*(2*(u2-u4)+9*u3)/6','L^2=(15/2)*u3*L-(1/2)*L*(2*(u2-u4)+9*u3)'],
            'linear_substitution':{str(key):str(value)for key,value in substitution.items()},
            'all_secondweight_solutions_parametrized':True,'firstweight_slice_dimension':4,'secondweight_kernel_dimension':3,
            'complete_weight1_coefficients':{str((i,j)):str(s.expand(first[i,j]))for i,j in H.PAIRS},
            'complete_weight2_coefficients':{str((i,j)):str(s.expand(second[i,j]))for i,j in H.PAIRS},
            'all350_weight3_Jacobi':[str(s.expand(value))for value in J3],
            'all350_weight4_Jacobi':[str(s.expand(value))for value in J4],
            'nonzero_witness':None if witness is None else {'weight':witness[0],'row':witness[1],'polynomial':str(s.factor(witness[2]))},
            'proven_higher_Jacobi_height_lower_bound':height,'actual_shear_rank':5,'affine_symbol_base_dimension':8,
            'full_higher_dimension_upper':5+7-height,'PGL_image_upper':ambient,
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,corepath,sourcepath,previous]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n');print('NIL_B0_TSQUARE_K0 whole support image',ambient,flush=True)


def global_twoparameter_main():
    """Four complete two-parameter families and every exact linear wall."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TWOPARAMETER_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite four whole twoparameter families budget240s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(240)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TSQUARE_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;ETA=s.Symbol('eta')
    En=s.Matrix([0,0,0,z,w]);ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1])
    Hbase=s.Matrix([a,b,-2*c,0,0]);N=s.Matrix([b,0,0,0,0]);J=s.Matrix([b,c,0,0,0]);J2=s.Matrix([c,0,0,0,0])
    normal_ss=s.Matrix([0,0,0,z,-w]);normal_nil=s.Matrix([0,0,0,0,z]);P1,source1,P2,source2=sources(x)
    models=[('SS_MIXH',TAU*Hbase+N,Hbase,normal_ss,ETA*z*w,s.Integer(1)),
            ('SS_MIXN',Hbase+TAU*N,N,normal_ss,ETA*z*w,s.Integer(1)),
            ('SS_NILREG',TAU*J+J2,J,normal_ss,ETA*z*w,s.Integer(1)),
            ('NIL_B0_MIX',Hbase+TAU*N,Hbase+N,normal_nil,ETA*z*z,TAU-1)]
    records=[]
    for name,T,V,normal,K,excluded in models:
        P0=H.wedge(T,En)+H.wedge(V,normal)+K*H.wedge(ez,ew)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,P0,x),x).values())
        D1=operator(P0,P1,source1,x);D2=operator(P0,P2,source2,x);S=shear(P0,x)
        assert all(s.expand(value)==0 for value in D1*S)
        analysis,minor=nonkernel_matrix_certificate(D1,D2,S,(TAU,ETA),16)
        factors=s.factor_list(minor)[1];lines=[];unsupported=[]
        for factor,multiplicity in factors:
            polynomial=s.Poly(factor,TAU,ETA)
            if polynomial.total_degree()!=1:
                unsupported.append({'whole_curve_equation':str(factor),'multiplicity':int(multiplicity),
                                    'curve_dimension_upper':1,'uniform_trivial_higher_fiber_upper':39,
                                    'coarse_restored_image_upper':54,'coarse_bound_sufficient':False,
                                    'status':'whole_nonlinear_curve_requires_new_exact_functionfield_or_geometric_proof'})
                continue
            ca=polynomial.coeff_monomial(TAU);cb=polynomial.coeff_monomial(ETA);constant=polynomial.coeff_monomial(1)
            substitution={TAU:TAU,ETA:-(ca*TAU+constant)/cb}if cb!=0 else{TAU:-constant/ca,ETA:TAU}
            ex=s.expand(excluded.subs(substitution,simultaneous=True))
            line={'equation':str(factor),'parameterization':{str(key):str(value)for key,value in substitution.items()},'excluded_multiplicity_higher':ex==0}
            if ex==0:line['all_strata_dimension_criterion']=True
            else:
                matrices=[matrix.subs(substitution,simultaneous=True)for matrix in [D1,D2,S]]
                line['analysis']=nonkernel_univariate(*matrices,15,ex)
                line['all_strata_dimension_criterion']=line['analysis']['all_strata_dimension_criterion']
            lines.append(line)
        analysis.update({'exception_polynomial':str(minor),'all_irreducible_factors':[[str(factor),int(multiplicity)]for factor,multiplicity in factors],
                         'all_exceptional_lines':lines,'unsupported_whole_curves':unsupported,
                         'all_strata_dimension_criterion':analysis['G_dimension_criterion']and not unsupported and all(line['all_strata_dimension_criterion']for line in lines)})
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
        record={'type':name,'parameters':['tau','eta'],'affine_first_symbol_dimension_upper':11,
                'excluded_polynomial':str(excluded),'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},
                'P0_full_projective_Jacobi_zero':True,'actual_D':str(s.factor(actual)),
                'complete_D1_350_by30':serial_matrix(D1),'complete_D2_350_by9':serial_matrix(D2),
                'complete_actual_shear_30_by6':serial_matrix(S),'analysis':analysis}
        records.append(record)
        destination.write_text(json.dumps({'task_id':'G','all_passed':False,'status':'bounded_partial_completed_whole_families',
            'cases':records,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-started},indent=2)+'\n')
        print(name,'generichigher',analysis['higher_fiber_upper_bound'],'lines',len(lines),'nonlinearcurves',len(unsupported),'allGcriterion',analysis['all_strata_dimension_criterion'],flush=True)
    proofs=[ROOT/'research/classification/PLANE_RANK_ZERO_GLOBAL_VECTOR_SEMISIMPLE_FIRST_JETS.md',ROOT/'research/classification/PLANE_RANK_ZERO_GLOBAL_VECTOR_NILPOTENT_FIRST_JETS.md']
    result={'task_id':'G','all_passed':all(record['analysis']['all_strata_dimension_criterion']for record in records),'cases':records,
            'scope':'Four entire mD2/globalV two-parameter source families. Every linear maximal-minor wall is retained as a wholeline, with all exact point roots. Unsupported nonlinear curves retain explicit insufficient coarse bound and are not excluded. Generic threeparameter/wholehighmultiplicity/scalarfirstjets remain separate.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous]+proofs],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_hyperbola_main():
    """One whole rational maximal-minor curve; tau zero is empty on it."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_HYPERBOLA_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic();signal.alarm(90)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TWOPARAMETER_20261107.json'
    packet=json.loads(previous.read_text());old=packet['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    ETA=s.Symbol('eta');environment={'tau':TAU,'eta':ETA}
    original=next(case for case in packet['cases']if case['type']=='NIL_B0_MIX')
    matrices=[s.Matrix([[s.sympify(value,locals=environment)for value in row]for row in original[key]])for key in ['complete_D1_350_by30','complete_D2_350_by9','complete_actual_shear_30_by6']]
    cleared=[matrix.subs(ETA,1-3/TAU).applyfunc(lambda value:s.cancel(TAU*value))for matrix in matrices]
    assert all(s.denom(value)==1 for matrix in cleared for value in matrix)
    analysis=nonkernel_univariate(*cleared,15,TAU*(TAU-1))
    result={'task_id':'G','all_passed':analysis['all_strata_dimension_criterion'],
            'scope':'Entire NILb0MIX rational hyperbola eta*tau-tau+3=0. Its tau0 complement is empty; tau1 is excluded exactly as m>=3, all other roots and conjugates retained.',
            'whole_parameterization':'eta=1-3/tau, tau!=0','tau_zero_curve_complement_empty':True,
            'rank_equivalent_denominator_clearing_scalar':'tau',
            'full_cleared_D1_350_by30':serial_matrix(cleared[0]),'full_cleared_D2_350_by9':serial_matrix(cleared[1]),
            'full_cleared_actual_shear_30_by6':serial_matrix(cleared[2]),'analysis':analysis,
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print('whole NILb0MIX hyperbola image',analysis['PGL_image_upper_bound'],'allGcriterion',analysis['all_strata_dimension_criterion'],flush=True)


def global_nil_line_main():
    """Entire SS nil K0 line: exact shear slices, binary resultant and every root."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_NIL_LINE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole SS nil line compatibility budget240s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(240)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_HYPERBOLA_20261107.json'
    sourcepath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TWOPARAMETER_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes'];packet=json.loads(sourcepath.read_text())
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    original=next(case for case in packet['cases']if case['type']=='SS_NILREG');ETA=s.Symbol('eta')
    matrices=[s.Matrix([[s.sympify(value,locals={'tau':TAU,'eta':ETA}).subs(ETA,0)for value in row]for row in original[key]])for key in ['complete_D1_350_by30','complete_D2_350_by9','complete_actual_shear_30_by6']]
    x=s.Matrix(s.symbols('a b c z w'));template,vars1,_,_=sources(x)
    def exact_certificate(matrix,domain):
        dm=DomainMatrix.from_Matrix(matrix).convert_to(domain);rr,columns=dm.rref();rows=dm.transpose().rref()[1]
        reduced=rr.to_Matrix();free=[j for j in range(matrix.cols)if j not in columns];kernel=s.zeros(matrix.cols,len(free))
        for k,j in enumerate(free):
            kernel[j,k]=1
            for i,column in enumerate(columns):kernel[column,k]=-reduced[i,j]
        product=dm*DomainMatrix.from_Matrix(kernel).convert_to(domain)
        assert not any(entry for row in product.to_list()for entry in row)
        determinant=domain.to_sympy(dm.extract(list(rows),list(columns)).det())
        return {'rank':len(columns),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
                'minor_determinant':str(determinant),'complete_kernel':serial_matrix(kernel)},kernel
    def one(value=None):
        domain=s.QQ.frac_field(TAU)if value is None else s.QQ if value.is_Rational else s.QQ.algebraic_field(value)
        current=matrices if value is None else[matrix.subs(TAU,value)for matrix in matrices]
        D1,D2,S=[DomainMatrix.from_Matrix(matrix).convert_to(domain).to_Matrix()for matrix in current]
        c1,K1=exact_certificate(D1,domain);c2,K2=exact_certificate(D2,domain);cs,_=exact_certificate(S,domain)
        B=S[:,cs['minor_columns']];r=B.cols;combined=B.row_join(K1)
        pivots=DomainMatrix.from_Matrix(combined).convert_to(domain).rref()[1]
        assert list(pivots[:r])==list(range(r))and len(pivots)==K1.cols
        R=K1[:,[column-r for column in pivots[r:]]];basis=B.row_join(R)
        basisrows=DomainMatrix.from_Matrix(basis.T).convert_to(domain).rref()[1]
        basis_minor=basis[list(basisrows),:]
        basisdet=domain.to_sympy(DomainMatrix.from_Matrix(basis_minor).convert_to(domain).det())
        u=s.symbols('u0:'+str(R.cols));first=template.subs(dict(zip(vars1,R*s.Matrix(u))),simultaneous=True)
        Q=s.Matrix([s.Poly(s.sympify(entry)/2,*u,domain=domain).as_expr()for entry in H.all350(H.projective(H.bilinear(first,first,x),x),x)])
        inverse=DomainMatrix.from_Matrix(D2.extract(c2['minor_rows'],c2['minor_columns'])).convert_to(domain).inv().to_Matrix()
        residual=Q-D2[:,c2['minor_columns']]*inverse*Q.extract(c2['minor_rows'],[0])
        polynomials=[s.Poly(entry,*u,domain=domain)for entry in residual]
        gcd=None;chain=[]
        for row,polynomial in enumerate(polynomials):
            if polynomial.is_zero:continue
            updated=polynomial.monic()if gcd is None else s.polys.polytools.gcd(gcd,polynomial).monic()
            if gcd is None or updated!=gcd:
                gcd=updated;chain.append({'row':row,'polynomial':str(polynomial.as_expr()),'cumulative_gcd':str(gcd.as_expr())})
            if gcd.total_degree()==0:break
        height=0 if gcd is None else 1 if gcd.total_degree()!=0 else 2
        base=10 if value is None else 9;higher=r+len(u)-height+c2['nullity'];ambient=base+higher+5
        record={'parameter':None if value is None else str(value),'coefficient_field':str(domain),'D1':c1,'D2':c2,'actual_shear':cs,
                'full_firstweight_slice_30_by_n':serial_matrix(R),'slice_variables':list(map(str,u)),
                'shear_plus_slice_basis_rows':list(basisrows),'shear_plus_slice_basis_minor':serial_matrix(basis_minor),
                'full_Q2_350':[str(entry)for entry in Q],'D2_pivot_inverse':serial_matrix(inverse),
                'full_compatibility_350':[str(polynomial.as_expr())for polynomial in polynomials],
                'exact_nonzero_or_gcd_chain':chain,'proven_compatibility_height_lower_bound':height,
                'affine_symbol_base_dimension':base,'full_higher_dimension_upper':higher,'PGL_image_upper':ambient,
                'G_dimension_criterion':ambient<=24}
        exception=None
        if value is None:
            assert len(u)==2
            monomials=[u[0]**2,u[0]*u[1],u[1]**2]
            coefficient_matrix=s.Matrix([[polynomial.coeff_monomial(monomial)for monomial in monomials]for polynomial in polynomials])
            selected=DomainMatrix.from_Matrix(coefficient_matrix.T).convert_to(domain).rref()[1]
            generators=[polynomials[row].as_expr()for row in selected]
            candidates=generators+[sum(integer**i*g for i,g in enumerate(generators))for integer in range(5)]
            resultant_certificate=None
            for p,q in combinations(candidates,2):
                pp=s.Poly(p,*u,domain=domain);qq=s.Poly(q,*u,domain=domain)
                a,b,c=[pp.coeff_monomial(monomial)for monomial in monomials];d,e,f=[qq.coeff_monomial(monomial)for monomial in monomials]
                sylvester=s.Matrix([[a,b,c,0],[0,a,b,c],[d,e,f,0],[0,d,e,f]])
                resultant=domain.to_sympy(DomainMatrix.from_Matrix(sylvester).convert_to(domain).det())
                if resultant!=0:
                    resultant_certificate={'first_binary_quadric':str(p),'second_binary_quadric':str(q),'Sylvester_matrix':serial_matrix(sylvester),'nonzero_binary_resultant':str(s.factor(resultant))}
                    break
            expressions=[s.sympify(c1['minor_determinant']),s.sympify(c2['minor_determinant']),s.sympify(cs['minor_determinant']),basisdet]
            expressions.extend(s.denom(entry)for entry in list(R)+list(inverse))
            if resultant_certificate:expressions.append(s.sympify(resultant_certificate['nonzero_binary_resultant']))
            factors={}
            for expression in expressions:
                for part in s.fraction(s.cancel(expression)):
                    for factor,multiplicity in s.factor_list(part,TAU)[1]:
                        canonical=s.Poly(factor,TAU).monic().as_expr();factors[str(canonical)]=canonical
            exception=s.prod(factors.values())
            record.update({'binary_resultant_certificate':resultant_certificate,'generic_exception_polynomial':str(s.factor(exception)),
                           'all_irreducible_exception_factors':sorted(factors),'selected_independent_compatibility_rows':list(selected)})
        return record,exception
    generic,exception=one();walls=[]
    destination.write_text(json.dumps({'task_id':'G','all_passed':False,'status':'bounded_partial_generic_line_complete',
        'generic':generic,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'producer_bytes':Path(__file__).read_text()},indent=2)+'\n')
    print('SS_NILREG K0 generic height',generic['proven_compatibility_height_lower_bound'],'image',generic['PGL_image_upper'],flush=True)
    for factor,multiplicity in s.factor_list(exception,TAU)[1]:
        polynomial=s.Poly(factor,TAU)
        value=-polynomial.nth(0)/polynomial.nth(1)if polynomial.degree()==1 else s.CRootOf(polynomial.as_expr(),0)
        record,_=one(value);record['minimal_polynomial']=str(factor);record['all_conjugate_roots_covered']=True
        walls.append(record)
        print('SS_NILREG K0 root',value,'height',record['proven_compatibility_height_lower_bound'],'image',record['PGL_image_upper'],flush=True)
    result={'task_id':'G','all_passed':generic['G_dimension_criterion']and all(record['G_dimension_criterion']for record in walls),
            'generic':generic,'all_exceptional_points':walls,
            'scope':'Entire SSNILREG K0 line via complete actual shear slice and full weight2 compatibility. Binary resultant controls the generic two-variable locus; every exact basis/rank/resultant root restored with its fullsource and P2kernel3. No missing point or chosen higherpoint.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous,sourcepath]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_threeparameter_main():
    """Two whole three-parameter sources; every linear plane and rational curve wall."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_THREEPARAMETER_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole threeparameter walls budget300s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(300)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_NIL_LINE_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    ETA,KAPPA=s.symbols('eta kappa');parameters=(TAU,ETA,KAPPA)
    def plane_analysis(matrices,excluded):
        D1,D2,S=matrices;analysis,minor=nonkernel_matrix_certificate(D1,D2,S,(TAU,ETA),16)
        factors=s.factor_list(minor)[1];curves=[];unsupported=[]
        for factor,multiplicity in factors:
            polynomial=s.Poly(factor,TAU,ETA)
            if polynomial.degree(ETA)==1:
                A=s.Poly(factor,ETA).coeff_monomial(ETA);B=s.Poly(factor,ETA).coeff_monomial(1)
                substitution={TAU:TAU,ETA:-B/A};univariate_A=A;univariate_B=B
            elif polynomial.degree(TAU)==1:
                A=s.Poly(factor,TAU).coeff_monomial(TAU);B=s.Poly(factor,TAU).coeff_monomial(1)
                substitution={TAU:(-B/A).subs(ETA,TAU),ETA:TAU};univariate_A=A.subs(ETA,TAU);univariate_B=B.subs(ETA,TAU)
            else:
                unsupported.append({'whole_curve_equation':str(factor),'dimension_upper':1,'coarse_higher_upper':39,'coarse_ambient_upper':54,'criterion_sufficient':False})
                continue
            assert s.Poly(s.gcd(univariate_A,univariate_B),TAU).degree()==0
            pulled=[matrix.subs(substitution,simultaneous=True).applyfunc(s.cancel)for matrix in matrices]
            denominator=s.Integer(1)
            for matrix in pulled:
                for value in matrix:denominator=s.lcm(denominator,s.denom(value))
            denominator=s.factor(denominator)
            cleared=[matrix.applyfunc(lambda value:s.cancel(denominator*value))for matrix in pulled]
            assert all(s.denom(value)==1 for matrix in cleared for value in matrix)
            ex=s.cancel(excluded.subs(substitution,simultaneous=True));exnumerator=s.fraction(ex)[0]
            curve={'equation':str(factor),'parameterization':{str(key):str(value)for key,value in substitution.items()},
                   'graph_A':str(univariate_A),'graph_B':str(univariate_B),'coprime_A_B':True,
                   'A_zero_curve_complement_empty':True,'rank_equivalent_denominator_clearing':str(denominator),
                   'excluded_multiplicity_higher':exnumerator==0}
            if exnumerator==0:curve['all_strata_dimension_criterion']=True
            else:
                curve['full_cleared_D1_350_by30']=serial_matrix(cleared[0]);curve['full_cleared_D2_350_by9']=serial_matrix(cleared[1]);curve['full_cleared_actual_shear_30_by6']=serial_matrix(cleared[2])
                curve['analysis']=nonkernel_univariate(*cleared,15,s.expand(exnumerator*denominator*univariate_A))
                curve['all_strata_dimension_criterion']=curve['analysis']['all_strata_dimension_criterion']
            curves.append(curve)
        analysis.update({'exception_polynomial':str(minor),'all_irreducible_factors':[[str(factor),int(multiplicity)]for factor,multiplicity in factors],
                         'all_whole_curve_walls':curves,'unsupported_whole_curves':unsupported,
                         'all_strata_dimension_criterion':analysis['G_dimension_criterion']and not unsupported and all(curve['all_strata_dimension_criterion']for curve in curves)})
        return analysis
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    En=s.Matrix([0,0,0,z,w]);ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1])
    Hbase=s.Matrix([a,b,-2*c,0,0]);N=s.Matrix([b,0,0,0,0])
    normal_ss=s.Matrix([0,0,0,z,-w]);normal_nil=s.Matrix([0,0,0,0,z])
    Tcartan=s.Matrix([-(1+TAU)*a/3,(2-TAU)*b/3,(2*TAU-1)*c/3,0,0])
    Vcartan=s.Matrix([-(1+ETA)*a/3,(2-ETA)*b/3,(2*ETA-1)*c/3,0,0])
    models=[('SS_MIX',TAU*Hbase+ETA*N,Hbase+N,normal_ss,KAPPA*z*w),
            ('NIL_B0_CARTAN',Tcartan,Vcartan,normal_nil,KAPPA*z*z)]
    P1,source1,P2,source2=sources(x);records=[]
    for name,T,V,normal,K in models:
        P0=H.wedge(T,En)+H.wedge(V,normal)+K*H.wedge(ez,ew)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,P0,x),x).values())
        matrices=[operator(P0,P1,source1,x),operator(P0,P2,source2,x),shear(P0,x)]
        assert all(s.expand(value)==0 for value in matrices[0]*matrices[2])
        analysis,minor=nonkernel_matrix_certificate(*matrices,parameters,17)
        factors=s.factor_list(minor)[1];planes=[];unsupported=[]
        for factor,multiplicity in factors:
            polynomial=s.Poly(factor,*parameters)
            if polynomial.total_degree()!=1:
                graph=[]
                for variable in parameters:
                    if polynomial.degree(variable)==1:
                        univariate=s.Poly(factor,variable);graph.append({'solved_variable':str(variable),'A':str(univariate.coeff_monomial(variable)),'B':str(univariate.coeff_monomial(1))})
                unsupported.append({'whole_surface_equation':str(factor),'dimension_upper':2,'possible_rational_graphs':graph,
                                    'vertical_complement_requires_all_A_B_zero_fibers':True,'coarse_higher_upper':39,'coarse_ambient_upper':55,'criterion_sufficient':False})
                continue
            coefficients=[polynomial.coeff_monomial(variable)for variable in parameters];constant=polynomial.coeff_monomial(1)
            pivot=max(i for i,coefficient in enumerate(coefficients)if coefficient!=0);free=[i for i in range(3)if i!=pivot]
            p,q=s.symbols('plane_p plane_q');temporary={parameters[free[0]]:p,parameters[free[1]]:q}
            expression=-(constant+sum(coefficients[i]*parameters[i]for i in free))/coefficients[pivot]
            temporary[parameters[pivot]]=expression.subs(temporary,simultaneous=True)
            substitution={variable:expression.subs({p:TAU,q:ETA},simultaneous=True)for variable,expression in temporary.items()}
            assert s.expand(factor.subs(substitution,simultaneous=True))==0
            excluded=s.expand((TAU-ETA).subs(substitution,simultaneous=True))
            plane={'equation':str(factor),'parameterization':{str(key):str(value)for key,value in substitution.items()},'excluded_multiplicity_higher':excluded==0}
            if excluded==0:plane['all_strata_dimension_criterion']=True
            else:
                pulled=[matrix.subs(substitution,simultaneous=True)for matrix in matrices]
                plane['full_D1_350_by30']=serial_matrix(pulled[0]);plane['full_D2_350_by9']=serial_matrix(pulled[1]);plane['full_actual_shear_30_by6']=serial_matrix(pulled[2])
                plane['analysis']=plane_analysis(pulled,excluded)
                plane['all_strata_dimension_criterion']=plane['analysis']['all_strata_dimension_criterion']
            planes.append(plane)
            destination.write_text(json.dumps({'task_id':'G','all_passed':False,'status':'bounded_partial_whole_plane_walls',
                'current_type':name,'completed_whole_planes':planes,'prior_completed_sources':records,
                'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text()},indent=2)+'\n')
            print(name,'whole plane',factor,'allGcriterion',plane['all_strata_dimension_criterion'],flush=True)
        analysis.update({'exception_polynomial':str(minor),'all_irreducible_factors':[[str(factor),int(multiplicity)]for factor,multiplicity in factors],
                         'all_whole_plane_walls':planes,'unsupported_whole_surfaces':unsupported,
                         'all_strata_dimension_criterion':analysis['G_dimension_criterion']and not unsupported and all(plane['all_strata_dimension_criterion']for plane in planes)})
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
        generic_empty=analysis['D1']['nullity']==6 and analysis['D2']['nullity']==0 and analysis['actual_shear']['rank']==6
        record={'type':name,'parameters':['tau','eta','kappa'],'affine_symbol_base_dimension':12,'excluded_polynomial':'tau-eta',
                'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},'actual_D':str(s.factor(actual)),
                'P0_full_projective_Jacobi_zero':True,'full_D1_350_by30':serial_matrix(matrices[0]),'full_D2_350_by9':serial_matrix(matrices[1]),
                'full_actual_shear_30_by6':serial_matrix(matrices[2]),'generic_actual_finite_shear_eliminates_all_higher':generic_empty,
                'generic_reduced_irreducible_quintic_source_empty':generic_empty,'analysis':analysis}
        records.append(record)
        destination.write_text(json.dumps({'task_id':'G','all_passed':False,'status':'bounded_partial_completed_whole_sources','cases':records,
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text()},indent=2)+'\n')
        print(name,'genericimage',analysis['PGL_image_upper_bound'],'planes',len(planes),'nonlinearsurfaces',len(unsupported),'allGcriterion',analysis['all_strata_dimension_criterion'],flush=True)
    proofs=[ROOT/'research/classification/PLANE_RANK_ZERO_GLOBAL_VECTOR_SEMISIMPLE_FIRST_JETS.md',ROOT/'research/classification/PLANE_RANK_ZERO_GLOBAL_VECTOR_NILPOTENT_FIRST_JETS.md']
    result={'task_id':'G','all_passed':all(record['analysis']['all_strata_dimension_criterion']for record in records),'cases':records,
            'scope':'Two entire remaining mD2/globalV threeparameter sources. Generic fullsource and every linear plane wall, every linear-in-one-variable rational curve wall and every exact finite root retained. Nonlinear surfaces and unsupported curves remain entire explicit frontiers, never discarded by codimension alone.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous]+proofs],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_cartan_surface_main():
    """The two remaining Cartan surfaces are one actual orbit-equivalent graph."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_CARTAN_SURFACE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole Cartan graph-surface budget240s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(240)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_THREEPARAMETER_20261107.json'
    packet=json.loads(previous.read_text());old=packet['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    ETA,KAPPA=s.symbols('eta kappa');x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    environment={str(value):value for value in list(x)+[TAU,ETA,KAPPA]}
    original=next(case for case in packet['cases']if case['type']=='NIL_B0_CARTAN')
    P0=s.zeros(5)
    for i,j in H.PAIRS:P0[i,j]=s.sympify(original['normal_form_P0'][str((i,j))],locals=environment);P0[j,i]=-P0[i,j]
    g=s.eye(5);g[0,0]=g[1,1]=0;g[0,1]=g[1,0]=1
    parameter_involution={TAU:1-TAU,ETA:1-ETA,KAPPA:-KAPPA}
    transformed=-g*P0.subs({a:b,b:a},simultaneous=True)*g.T
    assert all(s.expand(value)==0 for value in transformed-P0.subs(parameter_involution,simultaneous=True))
    F1=ETA*TAU-3*ETA+3*KAPPA*TAU+TAU
    F2=ETA*TAU+2*ETA+3*KAPPA*TAU-3*KAPPA-2*TAU-1
    assert s.expand(F1.subs(parameter_involution,simultaneous=True)-F2)==0
    substitution={KAPPA:-ETA/3+ETA/TAU-s.Rational(1,3)}
    assert s.cancel(F1.subs(substitution))==0
    assert F1.subs(TAU,0)==-3*ETA and F2.subs(TAU,1)==3*ETA-3
    matrices=[s.Matrix([[s.sympify(value,locals=environment)for value in row]for row in original[key]])for key in ['full_D1_350_by30','full_D2_350_by9','full_actual_shear_30_by6']]
    cleared=[matrix.subs(substitution,simultaneous=True).applyfunc(lambda value:s.cancel(TAU*value))for matrix in matrices]
    assert all(s.denom(value)==1 for matrix in cleared for value in matrix)
    excluded=TAU*(TAU-ETA)
    analysis,minor=nonkernel_matrix_certificate(*cleared,(TAU,ETA),16)
    factors=s.factor_list(minor)[1];curves=[];unsupported=[]
    for factor,multiplicity in factors:
        polynomial=s.Poly(factor,TAU,ETA)
        if polynomial.degree(ETA)==1:
            A=s.Poly(factor,ETA).coeff_monomial(ETA);B=s.Poly(factor,ETA).coeff_monomial(1)
            curve_substitution={TAU:TAU,ETA:-B/A};univariate_A=A;univariate_B=B
        elif polynomial.degree(TAU)==1:
            A=s.Poly(factor,TAU).coeff_monomial(TAU);B=s.Poly(factor,TAU).coeff_monomial(1)
            curve_substitution={TAU:(-B/A).subs(ETA,TAU),ETA:TAU};univariate_A=A.subs(ETA,TAU);univariate_B=B.subs(ETA,TAU)
        else:
            unsupported.append({'whole_curve_equation':str(factor),'curve_dimension_upper':1,'coarse_ambient_upper':54,'coarse_bound_sufficient':False})
            continue
        assert s.Poly(s.gcd(univariate_A,univariate_B),TAU).degree()==0
        pulled=[matrix.subs(curve_substitution,simultaneous=True).applyfunc(s.cancel)for matrix in cleared]
        denominator=s.Integer(1)
        for matrix in pulled:
            for value in matrix:denominator=s.lcm(denominator,s.denom(value))
        denominator=s.factor(denominator);curve_matrices=[matrix.applyfunc(lambda value:s.cancel(denominator*value))for matrix in pulled]
        assert all(s.denom(value)==1 for matrix in curve_matrices for value in matrix)
        ex=s.fraction(s.cancel(excluded.subs(curve_substitution,simultaneous=True)))[0]
        curve={'equation':str(factor),'parameterization':{str(key):str(value)for key,value in curve_substitution.items()},
               'graph_A':str(univariate_A),'graph_B':str(univariate_B),'coprime_A_B':True,'A_zero_curve_complement_empty':True,
               'rank_equivalent_denominator_clearing':str(denominator),'excluded_by_graph_domain_or_higher_multiplicity':ex==0}
        if ex==0:curve['all_strata_dimension_criterion']=True
        else:
            curve['full_cleared_D1_350_by30']=serial_matrix(curve_matrices[0]);curve['full_cleared_D2_350_by9']=serial_matrix(curve_matrices[1]);curve['full_cleared_actual_shear_30_by6']=serial_matrix(curve_matrices[2])
            curve['analysis']=nonkernel_univariate(*curve_matrices,15,s.expand(ex*denominator*univariate_A))
            curve['all_strata_dimension_criterion']=curve['analysis']['all_strata_dimension_criterion']
        curves.append(curve)
        destination.write_text(json.dumps({'task_id':'G','all_passed':False,'status':'bounded_partial_completed_whole_graph_curves',
            'completed_curves':curves,'unsupported_whole_curves':unsupported,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'producer_bytes':Path(__file__).read_text()},indent=2)+'\n')
        print('Cartan whole surface curve',factor,'allGcriterion',curve['all_strata_dimension_criterion'],flush=True)
    analysis.update({'exception_polynomial':str(minor),'all_irreducible_factors':[[str(factor),int(multiplicity)]for factor,multiplicity in factors],
                     'all_whole_curve_walls':curves,'unsupported_whole_curves':unsupported,
                     'all_strata_dimension_criterion':analysis['G_dimension_criterion']and not unsupported and all(curve['all_strata_dimension_criterion']for curve in curves)})
    result={'task_id':'G','all_passed':analysis['all_strata_dimension_criterion'],
            'scope':'Both entire remaining NILb0Cartan mD2 surfaces. Actualbase permutation plus tensor-1 identifies them. One full rational graph is covered by allQQ2curves/roots; denominatorzero vertical lines have tau=eta so lie entirely outside m2.',
            'first_surface':str(F1),'second_surface':str(F2),'actual_parameter_involution':{str(key):str(value)for key,value in parameter_involution.items()},
            'actual_linear_permutation_plus_tensor_minus_one_verified':True,'full_P0_involution_identity_zero':True,
            'first_surface_graph':str(substitution[KAPPA]),'graph_domain':'tau!=0',
            'first_vertical_complement':'tau=eta=0, kappa arbitrary, higher multiplicity',
            'second_vertical_complement':'tau=eta=1, kappa arbitrary, higher multiplicity',
            'source_graph_exclusion':str(excluded),'rank_equivalent_denominator_clearing':'tau',
            'full_cleared_D1_350_by30':serial_matrix(cleared[0]),'full_cleared_D2_350_by9':serial_matrix(cleared[1]),
            'full_cleared_actual_shear_30_by6':serial_matrix(cleared[2]),'analysis':analysis,
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_commonimage_main():
    """Four whole independent pointwise-parallel first-symbol orbits, all multiplicities."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMMONIMAGE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_CARTAN_SURFACE_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    T=s.Matrix([b,0,0,0,0]);V=s.Matrix([c,0,0,0,0]);En=s.Matrix([0,0,0,z,w])
    normal_ss=s.Matrix([0,0,0,z,-w]);normal_nil=s.Matrix([0,0,0,0,z])
    ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1]);P1,source1,P2,source2=sources(x)
    models=[('SS_COMMON_IMAGE_K1',normal_ss,z*w,9),('SS_COMMON_IMAGE_K0',normal_ss,s.Integer(0),8),
            ('NIL_COMMON_IMAGE_K1',normal_nil,z*z,8),('NIL_COMMON_IMAGE_K0',normal_nil,s.Integer(0),7)]
    records=[]
    for name,normal,K,base in models:
        P0=H.wedge(T,En)+H.wedge(V,normal)+K*H.wedge(ez,ew)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,P0,x),x).values())
        D1=operator(P0,P1,source1,x);D2=operator(P0,P2,source2,x);S=shear(P0,x)
        assert D1*S==s.zeros(350,6)
        analysis,_=nonkernel_matrix_certificate(D1,D2,S,(),base+5)
        records.append({'type':name,'affine_first_symbol_dimension':base,'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},
                        'full_projective_P0_Jacobi_zero':True,'full_D1_350_by30':serial_matrix(D1),'full_D2_350_by9':serial_matrix(D2),
                        'full_actual_shear_30_by6':serial_matrix(S),'analysis':analysis})
        print(name,'higher',analysis['higher_fiber_upper_bound'],'image',analysis['PGL_image_upper_bound'],'Gcriterion',analysis['G_dimension_criterion'],flush=True)
    proof=ROOT/'research/classification/PLANE_RANK_ZERO_COMMON_IMAGE_FIRST_JETS.md'
    result={'task_id':'G','all_passed':all(record['analysis']['G_dimension_criterion']for record in records),'cases':records,
            'scope':'Four entire independent pointwise-parallel common-image globalV first-symbol orbits; all multiplicities via reduced-D modular lemma. DependentT=tV/scalar/zeroP0 branches remain separate.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous,proof]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_dependent_main():
    """Five fixed whole symbols and one whole line after hyperplane recognition."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_DEPENDENT_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite dependent six-whole-source budget90s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(90)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMMONIMAGE_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    Hbase=s.Matrix([a,b,-2*c,0,0]);N=s.Matrix([b,0,0,0,0]);zero=s.zeros(5,1)
    En=s.Matrix([0,0,0,z,w]);normal_ss=s.Matrix([0,0,0,z,-w]);normal_nil=s.Matrix([0,0,0,0,z])
    ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1]);P1,source1,P2,source2=sources(x)
    models=[('SS_H_RESONANCE',zero,Hbase,normal_ss,-5*z*w,7,False),
            ('SS_N_LINE',TAU*N,N,normal_ss,s.Integer(0),7,True),
            ('NIL_H_RESONANCE',zero,Hbase,normal_nil,s.Rational(5,2)*z*z,6,False),
            ('NIL_N_TRACE_ZERO',zero,N,normal_nil,s.Integer(0),5,False),
            ('NIL_N_TRACE_ONE_K1',N,N,normal_nil,z*z,7,False),
            ('NIL_N_TRACE_ONE_K0',N,N,normal_nil,s.Integer(0),6,False)]
    def hyperplane_certificate(analysis,semisimple):
        kernel=s.Matrix([[s.sympify(value,locals={'tau':TAU})for value in row]for row in analysis['D1']['complete_kernel']])
        C=kernel[27:30,:];D=kernel[18:21,:]
        Czero=all(s.cancel(value)==0 for value in C);Dzero=all(s.cancel(value)==0 for value in D)
        empty=Czero or semisimple and Dzero
        analysis.update({'full_C_projection':serial_matrix(C),'full_D_projection':serial_matrix(D),
                         'required_C_nonzero':True,'required_D_nonzero':semisimple,
                         'all_C_zero':Czero,'all_D_zero':Dzero,'G_empty_required_normal_hyperplane':empty})
        if empty:analysis['G_dimension_criterion']=True
    records=[]
    for name,T,V,normal,K,base,variable in models:
        P0=H.wedge(T,En)+H.wedge(V,normal)+K*H.wedge(ez,ew)
        assert all(s.expand(value)==0 for value in H.projective(H.bilinear(P0,P0,x),x).values())
        D1=operator(P0,P1,source1,x);D2=operator(P0,P2,source2,x);S=shear(P0,x)
        assert all(s.expand(value)==0 for value in D1*S)
        if variable:
            analysis=nonkernel_univariate(D1,D2,S,base+5,s.Integer(1))
            for wall in analysis['all_exceptional_roots']:hyperplane_certificate(wall,True)
        else:analysis,_=nonkernel_matrix_certificate(D1,D2,S,(),base+5)
        hyperplane_certificate(analysis,normal==normal_ss)
        analysis['all_strata_dimension_criterion']=analysis['G_dimension_criterion']and all(wall['G_dimension_criterion']for wall in analysis.get('all_exceptional_roots',[]))
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]));assert actual==0
        record={'type':name,'affine_first_symbol_dimension_upper':base,'varying_parameter':'tau'if variable else None,
                'normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},'full_P0_projective_Jacobi_zero':True,'actual_D0_zero':True,
                'full_D1_350_by30':serial_matrix(D1),'full_D2_350_by9':serial_matrix(D2),'full_actual_shear_30_by6':serial_matrix(S),'analysis':analysis}
        records.append(record)
        print(name,'higher',analysis['higher_fiber_upper_bound'],'image',analysis['PGL_image_upper_bound'],'allGcriterion',analysis['all_strata_dimension_criterion'],flush=True)
    proof=ROOT/'research/classification/PLANE_RANK_ZERO_DEPENDENT_GLOBAL_VECTOR_HYPERPLANE_REDUCTION.md'
    result={'task_id':'G','all_passed':all(record['analysis']['all_strata_dimension_criterion']for record in records),'cases':records,
            'scope':'Five fixed whole dependent symbols plus entire SSNline, after full-Jordan hyperplane reduction. C/D nonzero are necessary for irreducible actualD; zero projection proves Gempty. Any nonclosing higher bound remains entire exact source, no newjet orsample.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous,proof]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def global_ore_main():
    """Whole linear-stage Ore identities, anchor support, and eight fixed affine modules."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ORE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole Ore recognition and8module budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_DEPENDENT_20261107.json'
    packet=json.loads(previous.read_text());old=packet['producer_bytes']
    parse=lambda value:{node.name:ast.dump(node,include_attributes=False)for node in ast.parse(value).body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=parse(old);after=parse(Path(__file__).read_text());assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;y=s.Matrix([b,c,z,w]);zero=s.Integer(0)
    d0,d1,b0,b1,e0,e1,f=s.symbols('d0 d1 b0 b1 e0 e1 f');q0,q1,q2=s.symbols('q0 q1 q2')
    generic_quad=q0*z*z+q1*z*w+q2*w*w
    def makeQ(kind,D,B,E,F,quad,t=TAU):
        if kind in ['ss','ss1','ss3']:
            p=t+1;qq=t-1;XB=s.Matrix([t*z*w,qq*w*w]);XE=s.Matrix([-p*z*z,-t*z*w])
            XF=s.Matrix([w*w,0])if kind=='ss3'else s.Matrix([z*w,w*w])if kind=='ss1'else s.zeros(2,1)
            xb=B[0]*XB+E[0]*XE;xc=B[1]*XB+E[1]*XE+F*XF
            hb=-qq*w*B[1]/2+p*z*E[1]/2;hc=qq*w*B[0]/2-p*z*E[0]/2
        elif kind=='n1':
            XD=s.Matrix([2*z*z,z*(z+2*w)]);XE=s.Matrix([z*(z+2*w),2*w*(z+w)])
            xb=D[0]*XD+E[0]*XE;xc=D[1]*XD+E[1]*XE
            hb=-z*D[1]-(z+w)*E[1];hc=z*D[0]+(z+w)*E[0]
        else:
            XD=s.Matrix([0,z*z]);XB=s.Matrix([-2*z*w,-3*w*w]);XE=s.Matrix([z*z,2*z*w]);XF=s.Matrix([z*z,z*w])
            xb=D[0]*XD+B[0]*XB+E[0]*XE;xc=D[1]*XD+B[1]*XB+E[1]*XE+F*XF
            hb=2*w*B[1]-z*E[1];hc=-2*w*B[0]+z*E[0]
        Q=s.zeros(4);Q[0,1]=s.expand(hb*b+hc*c+quad)
        Q[0,2]=s.expand(xb[0]);Q[0,3]=s.expand(xb[1]);Q[1,2]=s.expand(xc[0]);Q[1,3]=s.expand(xc[1])
        for i,j in combinations(range(4),2):Q[j,i]=-Q[i,j]
        return Q,hb,hc,xb,xc
    def pad(Q):
        P=s.zeros(5)
        for i,j in combinations(range(4),2):P[i+1,j+1]=Q[i,j];P[j+1,i+1]=-Q[i,j]
        return P
    def Wbasis(u,W0):
        return s.Matrix([z*(u[0]*b+u[1]*c)+w*(u[2]*b+u[3]*c)+u[4]*z*z+u[5]*z*w+u[6]*w*w,
                         z*(u[7]*b+u[8]*c)+w*(u[9]*b+u[10]*c)+u[11]*z*z+u[12]*z*w+u[13]*w*w,
                         W0[2]+u[14]*z*w+u[15]*w*w,W0[3]+u[16]*z*z+u[17]*z*w])
    def actualD(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    template,vars1,template2,vars2=sources(x);coverage=[]
    source_models=[('ss','SS_N_LINE',None),('ss1','SS_N_LINE',s.Integer(1)),('ss3','SS_N_LINE',s.Integer(3)),
                   ('n1','NIL_N_TRACE_ONE_K0',None),('n0','NIL_N_TRACE_ZERO',None)]
    for kind,name,value in source_models:
        original=next(case for case in packet['cases']if case['type']==name);analysis=original['analysis']
        if value is not None:analysis=next(wall for wall in analysis['all_exceptional_roots']if wall['parameter']==str(value))
        K1=s.Matrix([[s.sympify(entry,locals={'tau':TAU})for entry in row]for row in analysis['D1']['complete_kernel']]);k=s.symbols('k0:'+str(K1.cols))
        source=K1*s.Matrix(k);P1=template.subs(dict(zip(vars1,source)),simultaneous=True)
        P0=s.zeros(5)
        for i,j in H.PAIRS:
            P0[i,j]=s.sympify(original['normal_form_P0'][str((i,j))],locals={str(v):v for v in list(x)+[TAU]});P0[j,i]=-P0[i,j]
        if value is not None:P0=P0.subs(TAU,value)
        P=P0+P1+template2;L=s.Matrix([0]+[s.diff(P[0,j],a)for j in range(1,5)])
        ore=P-H.wedge(x,L);assert all(s.expand(s.diff(entry,a))==0 for entry in ore)
        D=[source[19],source[20]];B=[source[22],source[23]];E=[source[25],source[26]]
        F=source[12]if kind=='ss1'else source[29]if kind=='ss3'else source[7]if kind=='n0'else zero
        t=TAU if kind=='ss'else value if value is not None else zero
        Q,*_=makeQ(kind,D,B,E,F,template2[1,2],t)
        assert all(s.expand(entry)==0 for entry in ore.extract([1,2,3,4],[1,2,3,4])-Q)
        W=s.Matrix([ore[0,j]for j in range(1,5)]);W0=s.Matrix([P0[0,j]for j in range(1,5)])
        allowed=[[b*z,c*z,b*w,c*w,z*z,z*w,w*w],[b*z,c*z,b*w,c*w,z*z,z*w,w*w],[z*w,w*w],[z*z,z*w]]
        coefficients=[]
        for i,monomials in enumerate(allowed):
            poly=s.Poly(s.expand(W[i]-W0[i]),*y)
            coeff=[poly.coeff_monomial(monomial)for monomial in monomials];coefficients.extend(coeff)
            assert s.expand(poly.as_expr()-sum(v*m for v,m in zip(coeff,monomials)))==0
        parameter_list=list(k)+list(vars2)
        qcoords=B+E if kind.startswith('ss')else D+E if kind=='n1'else D+B+E
        if kind in ['ss1','ss3','n0']:qcoords=qcoords+[F]
        outputs=qcoords+list(vars2[6:])+coefficients
        assert len(outputs)==len(parameter_list)
        change=s.Matrix([[s.expand(output).coeff(v)for v in parameter_list]for output in outputs])
        determinant=DomainMatrix.from_Matrix(change).to_field().det();assert determinant!=0
        detexpr=DomainMatrix.from_Matrix(change).to_field().domain.to_sympy(determinant)
        assert not s.sympify(detexpr).has(TAU)
        coverage.append({'type':kind,'full_D1_kernel':serial_matrix(K1),'complete_Euler_gauge_L':list(map(str,L)),
                         'full_source_a_independent_after_Euler_gauge':True,'Q_expected_identity_zero':True,
                         'W_in_complete18source':True,'change_source_to_Q_and_W':serial_matrix(change),
                         'nonzero_constant_change_determinant':str(detexpr),
                         'Q_master_source':{str((i,j)):str(s.expand(Q[i,j]))for i,j in combinations(range(4),2)},
                         'complete_W_source':list(map(str,W))})
    D=(d0,d1);B=(b0,b1);E=(e0,e1);anchor_certificates=[]
    specs=[('ss',D,B,E,zero,TAU,[b0*e1-b1*e0],[b0,b1,e0,e1]),
           ('n1',D,B,E,zero,zero,[d0*e1-d1*e0],[d0,d1,e0,e1]),
           ('n0',D,B,E,f,zero,[2*(d0*b1-d1*b0)+f*e0,b0*e1-b1*e0+f*b0,d0*e1-d1*e0-f*d0],[d0,d1,b0,b1,e0,e1,f]),
           ('ss3',D,B,E,f,s.Integer(3),[b0*e1-b1*e0,f*e0],[b0,b1,e0,e1,f]),
           ('ss1',D,B,E,f,s.Integer(1),[b0*e1-b1*e0-f*e0,f*b0],[b0,b1,e0,e1,f])]
    for kind,DD,BB,EE,FF,t,gens,variables in specs:
        Q,*_=makeQ(kind,DD,BB,EE,FF,generic_quad,t);strict=H.bilinear(pad(Q),pad(Q),x)
        rows=[]
        for tri in combinations(range(1,5),3):
            polynomial=s.Poly(strict[tri]/2,*y)
            rows.extend([s.expand(polynomial.coeff_monomial(powers))for powers in H.monomials(4,3)])
        assert not any(entry.has(q0,q1,q2)for entry in rows)
        if kind=='ss':
            delta=gens[0];multipliers=[s.expand(entry).coeff(b0*e1)for entry in rows]
            assert all(s.expand(entry-coefficient*delta)==0 for entry,coefficient in zip(rows,multipliers))
            gcd=zero;bezout=[zero]*len(rows)
            for i,coefficient in enumerate(multipliers):
                if coefficient==0:continue
                aa,bb,newgcd=s.gcdex(gcd,coefficient,TAU)
                bezout=[s.cancel(aa*v)for v in bezout];bezout[i]+=bb;gcd=newgcd
                if s.Poly(gcd,TAU).degree()==0:break
            bezout=[s.cancel(entry/gcd)for entry in bezout]
            assert s.expand(sum(entry*row for entry,row in zip(bezout,rows))-delta)==0
            certificate={'all_rows_scalar_times_Delta':list(map(str,multipliers)),
                         'uniform_polynomial_Bezout_to_Delta':[[i,str(entry)]for i,entry in enumerate(bezout)if entry!=0]}
        else:
            monomials=[s.prod(variables[i]**exponent for i,exponent in enumerate(powers))for powers in H.monomials(len(variables),2)]
            M=s.Matrix([[s.Poly(entry,*variables).coeff_monomial(monomial)for monomial in monomials]for entry in rows])
            G=s.Matrix([[s.Poly(entry,*variables).coeff_monomial(monomial)for monomial in monomials]for entry in gens])
            selected=list(DomainMatrix.from_Matrix(M.T).to_field().rref()[1]);assert len(selected)==len(gens)
            basis=M[selected,:];columns=list(DomainMatrix.from_Matrix(basis).to_field().rref()[1])
            forward=G[:,columns]*basis[:,columns].inv();assert forward*basis==G
            gcolumns=list(DomainMatrix.from_Matrix(G).to_field().rref()[1]);back=M[:,gcolumns]*G[:,gcolumns].inv();assert back*G==M
            certificate={'selected_Jacobi_rows':selected,'anchor_generators_from_selected_Jacobi':serial_matrix(forward),
                         'all_Jacobi_rows_from_anchor_generators':serial_matrix(back)}
        anchor_certificates.append({'type':kind,'anchor_generators':list(map(str,gens)),
                                    'full80_strict_Q_Jacobi':list(map(str,rows)),'exact_ideal_equal_certificate':certificate})
    Qnil,*_=makeQ('n0',D,B,E,f,generic_quad);h,ss=s.symbols('h shear')
    g=s.eye(4);g[3,2]=h
    pushed=g*Qnil.subs(w,w-h*z)*g.T
    predicted,*_=makeQ('n0',(d0-h*e0-h*h*b0,d1-h*e1-h*h*b1),B,(e0+2*h*b0,e1+2*h*b1),f,generic_quad.subs(w,w-h*z))
    assert all(s.expand(entry)==0 for entry in pushed-predicted)
    g=s.eye(4);g[1,0]=ss
    pushed=g*Qnil.subs(c,c-ss*b)*g.T
    predicted,*_=makeQ('n0',(d0,d1+ss*d0),(b0,b1+ss*b0),(e0,e1+ss*e0),f,generic_quad)
    assert all(s.expand(entry)==0 for entry in pushed-predicted)
    gg=specs[2][6]
    assert s.expand(f*(4*d0*b0+e0*e0)-(2*d0*gg[1]-2*b0*gg[2]+e0*gg[0]))==0
    v0,v1,dd,bb,ee=s.symbols('v0 v1 dd bb ee');v=s.Matrix([v0,v1,0,0]);constant_factor_checks=[]
    for kind,DD,BB,EE,FF,t in [('ss',(0,0),(bb*v0,bb*v1),(ee*v0,ee*v1),zero,TAU),
                              ('n1',(dd*v0,dd*v1),(0,0),(ee*v0,ee*v1),zero,zero),
                              ('n0',(dd*v0,dd*v1),(bb*v0,bb*v1),(ee*v0,ee*v1),zero,zero),
                              ('ss1',(0,0),(0,bb),(0,ee),f,s.Integer(1)),
                              ('ss3',(0,0),(0,bb),(0,ee),f,s.Integer(3)),
                              ('n0',(0,dd),(0,bb),(0,ee),f,zero)]:
        Q,*_=makeQ(kind,DD,BB,EE,FF,generic_quad,t);vv=v if FF==0 else s.Matrix([0,1,0,0])
        assert all(s.expand(vv[i]*Q[j,k]-vv[j]*Q[i,k]+vv[k]*Q[i,j])==0 for i,j,k in combinations(range(4),3))
        assert all(s.expand(sum(vv[k]*s.diff(Q[i,j],y[k])for k in range(4)))==0 for i,j in combinations(range(4),2))
        constant_factor_checks.append({'type':kind,'constant_factor_and_invariance_exact':True,'first_column_zero_branch':FF!=0})
    u=s.symbols('u0:18');ua,ub,va,vb,lam=s.symbols('uz uw vz vw lam');shift_u=ua*z+ub*w;shift_v=va*z+vb*w
    canonical=[('A','n0',(1,0),(0,0),(0,1),1,zero,w*w,5,5,6,s.Matrix([0,0,0,b*z])),
               ('B','n0',(0,0),(1,0),(0,-1),1,zero,z*w,5,6,7,s.Matrix([0,0,0,b*z])),
               ('S3','ss3',(0,0),(1,0),(0,0),1,s.Integer(3),z*z,6,5,6,s.Matrix([0,0,4*b*z,2*b*w])),
               ('S1','ss1',(0,0),(0,-1),(1,0),1,s.Integer(1),z*w,6,5,6,s.Matrix([0,0,2*b*z,0]))]
    cases=[];cocycle_checks=[]
    for name,kind,DD,BB,EE,FF,t,residue,base,orbit0,orbit1,W0 in canonical:
        Q,hb,hc,xb,xc=makeQ(kind,DD,BB,EE,FF,generic_quad,t)
        coboundary=s.expand(-hb*shift_u-hc*shift_v+sum(xb[i]*s.diff(shift_v,[z,w][i])-xc[i]*s.diff(shift_u,[z,w][i])for i in range(2)))
        g=s.eye(4);g[0,2]=ua;g[0,3]=ub;g[1,2]=va;g[1,3]=vb
        pushed=g*Q.subs({b:b-shift_u,c:c-shift_v},simultaneous=True)*g.T
        predicted,*_=makeQ(kind,DD,BB,EE,FF,generic_quad+coboundary,t)
        assert all(s.expand(entry)==0 for entry in pushed-predicted)
        g=s.diag(1/lam,1/lam,1,1)
        pushed=lam*g*Q.subs({b:lam*b,c:lam*c},simultaneous=True)*g.T
        predicted,*_=makeQ(kind,DD,BB,EE,FF,generic_quad/lam,t)
        assert all(s.cancel(entry)==0 for entry in pushed-predicted)
        cocycle_checks.append({'type':name,'full_shear_identity_zero':True,'coboundary':str(coboundary),
                               'residual_quadratic':str(residue),'actual_diagonal_tensor_scales_residue_by_inverse_lam':True})
        for delta in [0,1]:
            Q,*_=makeQ(kind,DD,BB,EE,FF,delta*residue,t);Q5=pad(Q)
            assert all(s.expand(entry)==0 for entry in H.bilinear(Q5,Q5,x).values())
            W=Wbasis(u,W0);P=Q5+H.wedge(s.Matrix([1,0,0,0,0]),s.Matrix([0,*W]))
            rows=[s.sympify(entry)for entry in H.all350(H.projective(H.bilinear(P,P,x),x),x)]
            constants=s.Matrix([entry.subs(dict.fromkeys(u,0))for entry in rows])
            M=s.Matrix([[s.expand(entry).coeff(variable)for variable in u]for entry in rows])
            assert all(s.expand(entry-constants[i]-sum(M[i,j]*variable for j,variable in enumerate(u)))==0 for i,entry in enumerate(rows))
            cert,_=multivariate_rank_certificate(M,());rank=cert['rank'];kernel=s.Matrix([[s.sympify(entry)for entry in row]for row in cert['complete_kernel']])
            augmented=M.row_join(constants);arank=DomainMatrix.from_Matrix(augmented).to_field().rank()
            LW=s.zeros(4)
            for i,j in combinations(range(4),2):
                LW[i,j]=s.expand(sum(W[k]*s.diff(Q[i,j],y[k])-Q[k,j]*s.diff(W[i],y[k])-Q[i,k]*s.diff(W[j],y[k])for k in range(4)));LW[j,i]=-LW[i,j]
            module=[]
            for i,j,k in combinations(range(4),3):
                polynomial=s.Poly(s.expand(y[i]*LW[j,k]-y[j]*LW[i,k]+y[k]*LW[i,j]),*y)
                module.extend([polynomial.coeff_monomial(powers)for powers in H.monomials(4,4)])
            MC=s.Matrix([entry.subs(dict.fromkeys(u,0))for entry in module]);MM=s.Matrix([[s.expand(entry).coeff(variable)for variable in u]for entry in module])
            MA=MM.row_join(MC)
            assert DomainMatrix.from_Matrix(MM).to_field().rank()==rank
            assert DomainMatrix.from_Matrix(MA).to_field().rank()==arank
            assert DomainMatrix.from_Matrix(augmented.col_join(MA)).to_field().rank()==arank
            orbit=orbit0 if delta==0 else orbit1
            record={'type':name,'central_delta':delta,'fixed_symbol_base_dimension':base,'Q_orbit_dimension_upper':orbit,
                    'complete_Q':{str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},'W0':list(map(str,W0)),
                    'full350_affine_matrix':serial_matrix(M),'full350_affine_constant':list(map(str,constants)),
                    'full140_module_matrix':serial_matrix(MM),'full140_module_constant':list(map(str,MC)),
                    'full350_and_module140_same_affine_rowspace':True,'linear_certificate':cert,'augmented_rank':arank}
            if arank>rank:
                ri=cert['minor_rows'];ci=cert['minor_columns'];inverse=M.extract(ri,ci).inv();particular=s.zeros(18,1)
                partial=-inverse*constants.extract(ri,[0])
                for i,column in enumerate(ci):particular[column]=partial[i]
                residual=M*particular+constants;index=next(i for i,entry in enumerate(residual)if entry!=0)
                phi=s.zeros(1,350);phi[0,index]=1;weights=M[index:index+1,ci]*inverse
                for j,row in enumerate(ri):phi[0,row]-=weights[j]
                assert phi*M==s.zeros(1,18)and(phi*constants)[0]!=0
                record.update(status='whole_affine_module_inconsistent',cokernel_witness=[[i,str(entry)]for i,entry in enumerate(list(phi))if entry!=0],
                              nonzero_constant_contradiction=str((phi*constants)[0]),G_empty=True,G_criterion=True)
            else:
                ri=cert['minor_rows'];ci=cert['minor_columns'];particular=s.zeros(18,1)
                partial=-M.extract(ri,ci).inv()*constants.extract(ri,[0])
                for i,column in enumerate(ci):particular[column]=partial[i]
                assert M*particular+constants==s.zeros(350,1)
                Dfull=actualD(P);Dconstant=Dfull.subs(dict.fromkeys(u,0));Dlinear=s.Matrix([[s.expand(Dfull).coeff(variable)for variable in u]])
                assert s.expand(Dfull-Dconstant-(Dlinear*s.Matrix(u))[0])==0
                Dpoint=s.expand(Dconstant+(Dlinear*particular)[0]);directions=(Dlinear*kernel).applyfunc(s.expand)
                common=s.Poly(Dpoint,*x,domain=s.QQ)
                for entry in directions:
                    if entry!=0:common=s.polys.polytools.gcd(common,s.Poly(entry,*x,domain=s.QQ)).monic()
                factor=None
                if common.total_degree()>0:
                    for value,multiplicity in s.factor_list(common.as_expr(),*x)[1]:
                        if 0<s.Poly(value,*x).total_degree()<5:factor=value;break
                bound=base+orbit+kernel.cols+5
                record.update(status='complete_affine_module',particular_solution=list(map(str,particular)),full_W_fiber_dimension=kernel.cols,
                              whole_restored_image_upper=bound,actual_D_particular=str(Dpoint),actual_D_kernel_coefficients=list(map(str,directions)),
                              common_D_polynomial_gcd=str(common.as_expr()),proper_fixed_D_factor=None if factor is None else str(factor),
                              G_empty=factor is not None,G_criterion=factor is not None or bound<=24)
            cases.append(record)
            print(name,delta,'rank',rank,'aug',arank,'k',18-rank,'Gempty',record['G_empty'],'Gcriterion',record['G_criterion'],flush=True)
    proof=ROOT/'research/classification/PLANE_RANK_ZERO_RANK_ONE_NILPOTENT_ORE_REDUCTION.md'
    saito=ROOT/'research/classification/ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md'
    result={'task_id':'G','all_passed':all(record['G_criterion']for record in cases),'full_source_Ore_coverage':coverage,
            'strict_anchor_ideal_certificates':anchor_certificates,'nil_normal_and_base_shear_identities_zero':True,
            'nil_Pluckers_lift_to_conic_exact':True,'constant_factor_and_invariance_checks':constant_factor_checks,
            'central_coboundary_and_diagonal_checks':cocycle_checks,'cases':cases,
            'scope':'Entire three rank1nil K0 branches: constant-factor Q gives actualD essential3; all other Q lie in8 actual canonical orbits. Fullsource Euler maps, strictQJac ideals, normalizer/cocycle identities and all8 full350/140 affine modules retained. No GB/newjet/anchor sampling.',
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}for path in [HELPER,previous,proof,saito]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_cubic_main():
    """Two complete triple-line/m3/squarefree first connections, no point sampling."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite nonreduced triple-line two-model budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ORE_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old);after=functions(Path(__file__).read_text());assert len(before)==26
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x
    P1,v1,P2,v2=sources(x)
    def actualD(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    monbase=[tuple(v.count(i)for i in range(3))for v in __import__('itertools').combinations_with_replacement(range(3),2)]
    mon3=[base+(j,3-j)for base in monbase for j in range(4)]
    def leading3(poly):return s.Matrix([s.Poly(poly,*x).coeff_monomial(mon)for mon in mon3])
    results=[]
    for name,unipotent,B in [('scalar',False,3),('unipotent',True,6)]:
        P0=s.zeros(5);P0[1,2]=3*a*a;P0[2,1]=-P0[1,2]
        if unipotent:P0[2,4]=3*a*z;P0[4,2]=-P0[2,4]
        assert all(s.expand(value)==0 for value in H.all350(H.projective(H.bilinear(P0,P0,x),x),x))
        assert actualD(P0)==0
        allowed=[a*a*z**j*w**(3-j)for j in range(4)]if not unipotent else [a*a*z**3,a*z*z*(a*w+b*z),z*(a*w+b*z)**2]
        A=s.Matrix.hstack(*(leading3(poly)for poly in allowed));ca,_,_=rank_certificate(A,False)
        assert ca['rank']==len(allowed)
        Ann=s.Matrix.hstack(*A.T.nullspace()).T;assert Ann.cols==24 and Ann*A==s.zeros(Ann.rows,A.cols)
        assert rank_certificate(Ann,False)[0]['rank']==24-len(allowed)
        D3=leading3(actualD(P0+P1));D3M=D3.jacobian(v1)
        assert D3-D3M*s.Matrix(v1)==s.zeros(24,1)
        D1=operator(P0,P1,v1,x);D2=operator(P0,P2,v2,x)
        D1complete=D1.col_join(Ann*D3M)
        c1,_,K1=rank_certificate(D1complete,False);c2,_,K2=rank_certificate(D2,False)
        S=shear(P0,x);assert all(s.expand(value)==0 for value in D1complete*S)
        cs,_,_=rank_certificate(S,False)
        d=c1['nullity']+c2['nullity'];upper=B+5+d
        results.append({'type':name,'full_normal_form_P0':{str((i,j)):str(s.expand(P0[i,j]))for i,j in H.PAIRS},
            'P0_full350_Jacobi_zero_and_actualD_zero':True,'whole_fixed_plane_symbol_base_dimension_upper':B,
            'all24_leading_monomials':mon3,'full_parallel_leading_basis':list(map(str,allowed)),
            'parallel_basis_24byN':serial_matrix(A),'parallel_membership_full_annihilator':serial_matrix(Ann),
            'full_D3_linear_map_24by30':serial_matrix(D3M),'full350_weight1_Jacobi':serial_matrix(D1),
            'full_weight1_with_parallel_D3':serial_matrix(D1complete),'full350_weight2_Jacobi':serial_matrix(D2),
            'complete_actual_shear_30by6':serial_matrix(S),'weight1_certificate':c1,'weight2_certificate':c2,'actual_shear_certificate':cs,
            'all_inhomogeneous_P2_fibers_upper':c2['nullity'],'whole_higher_source_dimension_upper':d,
            'whole_projective_image_upper':upper,'G_dimension_criterion':upper<=24})
        print(name,'fullP1/P2 kernels',c1['nullity'],c2['nullity'],'whole image upper',upper,flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_FIRST_JETS.md'
    entrance=ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md'
    archive=ROOT/'research/classification/inputs/plane_rank_zero_pre_nonreduced_cubic_frozen_inputs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_G24_closed':all(row['G_dimension_criterion']for row in results),
            'scope':'Rank2 plane with zero cubic3L, multiplicity3 and generic squarefree normal cubic. Fullroot polynomial-degree proof gives exactly two first connections. Allfull350/30/9 and full24leading coefficients retained; m2,m3repeated,m4,2L+M outside.',
            'full_source1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},'full_source2':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},
            'cases':results,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,previous,proof,entrance,archive]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_cubic_leading_main():
    """Actual leading image closes both full connections, despite coarse bounds."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_LEADING_EXCLUSION_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite two-model actual-leading image budget60s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_20261107.json'
    packet=json.loads(previous.read_text());old=packet['producer_bytes'];assert packet['all_G24_closed']is False
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old);after=functions(Path(__file__).read_text());assert len(before)==27
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;loc={str(v):v for v in x}
    P1,v1,P2,v2=sources(x)
    def actualD(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    results=[]
    for row in packet['cases']:
        P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=s.sympify(row['full_normal_form_P0'][str((i,j))],locals=loc);P0[j,i]=-P0[i,j]
        D1=operator(P0,P1,v1,x);assert serial_matrix(D1)==row['full350_weight1_Jacobi']
        fullD=actualD(P0+P1+P2)
        mon3=[tuple(mon)for mon in row['all24_leading_monomials']]
        D3=s.Matrix([s.Poly(fullD,*x).coeff_monomial(mon)for mon in mon3]);L=D3.jacobian(v1)
        assert D3==L*s.Matrix(v1)and serial_matrix(L)==row['full_D3_linear_map_24by30']
        Ann=s.Matrix(row['parallel_membership_full_annihilator']);M=D1.col_join(Ann*L)
        cert,_,K=rank_certificate(M,False);assert cert==row['weight1_certificate']
        image=L*K;ci,_,_=rank_certificate(image,False)
        allowed=(2,0,0,3,0);at=mon3.index(allowed)
        assert all(image[i,j]==0 for i in range(24)if i!=at for j in range(image.cols))
        if row['type']=='scalar':assert ci['rank']==0
        else:assert ci['rank']==1
        u=s.symbols('u0:'+str(K.cols));G3=s.expand(sum(image[i,j]*u[j]*s.prod(v**e for v,e in zip(x,mon3[i]))for i in range(24)for j in range(K.cols)))
        assert s.expand(s.discriminant(G3.subs(w,1),z))==0
        results.append({'type':row['type'],'full150_P0':row['full_normal_form_P0'],'all24_leading_monomials':mon3,
            'full350_weight1_operator':serial_matrix(D1),'full_actual_D3_linear_map':serial_matrix(L),
            'complete_kernel_certificate':cert,'full_leading_image_matrix':serial_matrix(image),'image_rank_certificate':ci,
            'whole_leading_polynomial':str(s.factor(G3)),'only_possible_image_monomial':list(allowed),
            'full_actualD3_independent_of_all9_P2':True,'normal_cubic_discriminant_identically_zero':True,
            'entire_squarefree_m3_stratum_G_empty':True})
    proof=ROOT/'research/classification/PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_FIRST_JETS.md'
    archive=ROOT/'research/classification/inputs/plane_nonreduced_triple_line_pre_leading_frozen_inputs_20261107.json'
    result={'task_id':'G','all_passed':True,'all_two_whole_sources_G_empty':True,'cases':results,
            'scope':'Entire rank-two-plane R=3L, multiplicity3 and generic squarefree normal cubic. Fullsource root-degree proof gives two connections; full linear Jacobi and actual leading map make G3 zero or a cube. No higher jet, nonlinear elimination or generic-only bound.',
            'old_coarse_dimension_report_false_unchanged':True,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,previous,proof,archive]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def m4_two_line_casimir_main():
    """Whole nonreduced-R generic211 actual-D families, necessary canonical source."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite two-fixed/one-whole-line fullCasimir budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_LEADING_EXCLUSION_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old);after=functions(Path(__file__).read_text());assert len(before)==28
    assert all(after[name]==value for name,value in before.items())
    x=s.symbols('A B u v w');A,B,u,v,w=x
    mons=lambda d:[tuple(z.count(i)for i in range(5))for z in __import__('itertools').combinations_with_replacement(range(5),d)]
    mon2=mons(2);mon6=mons(6);indices={mon:i for i,mon in enumerate(mon6)}
    basis=[(i,j,mon)for i,j in combinations(range(5),2)for mon in mon2]
    def matrix_for(F):
        grad=[s.Poly(s.diff(F,t),*x)for t in x];M=s.MutableSparseMatrix(1075,150,{})
        for k,(i,j,mon)in enumerate(basis):
            for component,g,sgn in [(i,grad[j],1),(j,grad[i],-1)]:
                for power,coefficient in g.terms():
                    if coefficient==0:continue
                    at=component*210+indices[tuple(a+b for a,b in zip(mon,power))];M[at,k]+=sgn*coefficient
            for coordinate,component,sgn in [(i,j,1),(j,i,-1)]:
                if mon[coordinate]:
                    power=list(mon);power[coordinate]-=1;at=power.index(1)
                    M[1050+component*5+at,k]+=sgn*mon[coordinate]
        return s.Matrix(M)
    def algebraic_certificate(M,factor):
        poly=s.Poly(factor,TAU)
        if poly.degree()==1:
            theta=-poly.nth(0)/poly.nth(1);cert,_,_=rank_certificate(M.subs(TAU,theta),False)
            return {'minimal_polynomial':str(factor),'tau':str(theta),'all_complex_embeddings_covered':True,'source_certificate':cert}
        theta=s.CRootOf(poly.as_expr(),0);domain=s.QQ.alg_field(theta)
        dm=DomainMatrix.from_Matrix(M.subs(TAU,theta)).convert_to(domain);rr,columns=dm.rref();rows=dm.transpose().rref()[1]
        R=rr.to_Matrix();free=[j for j in range(M.cols)if j not in columns];K=s.zeros(M.cols,len(free))
        for k,j in enumerate(free):
            K[j,k]=1
            for i,pivot in enumerate(columns):K[pivot,k]=-R[i,j]
        kdm=DomainMatrix.from_Matrix(K).convert_to(domain);assert (dm*kdm).is_zero_matrix
        det=dm.extract(list(rows),list(columns)).det();assert det!=domain.zero
        cert={'rank':len(columns),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
              'minor_determinant':str(domain.to_sympy(det)),'complete_kernel':serial_matrix(K),'kernel_substitution_zero':True}
        return {'minimal_polynomial':str(factor),'primitive_root':str(theta),'all_complex_embeddings_covered':True,
                'number_field_degree':poly.degree(),'source_certificate':cert}
    records=[]
    for name,F,variable,Ddim,target in [
        ('sem_zero',u**4*A+u**2*v**2*B+v**5,False,18,7),
        ('sem_nonzero',u**4*A+u**2*v**2*B+u*v**4+v**5,False,19,6),
        ('Jordan_whole',u**2*v**2*A+u**3*v*B+u**5+TAU*u*v**4+v**5,True,19,6)]:
        M=matrix_for(F);cert,minor,K=rank_certificate(M,variable)
        assert all(s.cancel(value)==0 for value in M*K)
        linearZ=[0,0,0,0,0]if name=='sem_nonzero'else [-4*A,-2*B,u,0,0]if name=='sem_zero'else [u,-v,0,0,0]
        assert s.expand(sum(z*s.diff(F,t)for z,t in zip(linearZ,x)))==0
        factors=s.factor_list(minor)[1]if variable else [];walls=[]
        for factor,power in factors:
            wall=algebraic_certificate(M,factor);wall['minor_factor_multiplicity']=int(power)
            wall['whole_D_family_dimension_upper']=18;wall['whole_projective_image_upper']=17+wall['source_certificate']['nullity']
            wall['G_dimension_criterion']=wall['source_certificate']['nullity']<=7;walls.append(wall)
        upper=Ddim+cert['nullity']-1
        records.append({'type':name,'actual_D':str(F),'coefficient_parameter_tau':variable,
            'full1075by150_Casimir_and_divergence_matrix':serial_matrix(M),'full_canonical_source_certificate':cert,
            'generic_minor_all_factors':[[str(f),int(power)]for f,power in factors],
            'all_exceptional_fields':walls,'whole_D_family_dimension_upper':Ddim,'whole_projective_image_upper':upper,
            'generic_target_linear_source_dimension':target,'actual_linear_stabilizer_field':list(map(str,linearZ)),
            'G_dimension_criterion':cert['nullity']<=target and all(row['G_dimension_criterion']for row in walls)})
        print(name,'source dimension',cert['nullity'],'generic upper',upper,'wall dims',[row['source_certificate']['nullity']for row in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_M4_TWO_LINE_COLLISION_QUINTICS.md'
    geometry=ROOT/'research/classification/PLANE_M4_REPEATED_LEADING_GEOMETRY.md'
    identity=ROOT/'research/classification/ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md'
    archive=ROOT/'research/classification/inputs/plane_nonreduced_triple_line_m3_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_G24_closed':all(row['G_dimension_criterion']for row in records),
            'scope':'All highrank generic211 m4 sources whose rank2-plane zero cubic has support atmost2lines. Root collisions exclude rank3/unalignedsem and all3L highrank. Actual D reduces to two fixed sem forms or one whole Jordan line, allparameter walls covered by exact fields. Necessary div0/Casimir only; fullJacobi/actualPf can only restrict.',
            'all150_quadratic_basis':[[i,j,list(mon)]for i,j,mon in basis],'all210_degree6_monomials':mon6,'cases':records,
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,previous,proof,geometry,identity,archive]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def full_quadratic_casimir_matrix(F,x):
    """Full five-variable150 source; no vertex coefficients removed."""
    mons=lambda d:[tuple(z.count(i)for i in range(5))for z in __import__('itertools').combinations_with_replacement(range(5),d)]
    mon2=mons(2);mon6=mons(6);indices={mon:i for i,mon in enumerate(mon6)}
    basis=[(i,j,mon)for i,j in combinations(range(5),2)for mon in mon2]
    grad=[s.Poly(s.diff(F,t),*x)for t in x];M=s.MutableSparseMatrix(1075,150,{})
    for k,(i,j,mon)in enumerate(basis):
        for component,g,sgn in [(i,grad[j],1),(j,grad[i],-1)]:
            for power,coefficient in g.terms():
                if coefficient==0:continue
                at=component*210+indices[tuple(a+b for a,b in zip(mon,power))];M[at,k]+=sgn*coefficient
        for coordinate,component,sgn in [(i,j,1),(j,i,-1)]:
            if mon[coordinate]:
                power=list(mon);power[coordinate]-=1;at=power.index(1)
                M[1050+component*5+at,k]+=sgn*mon[coordinate]
    return s.Matrix(M),basis,mon6


def m4_nonreduced_squarefree_main():
    """Whole nonreduced-R squarefree m4 sources: six complete actual quintics."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_CASIMIR_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole six-squarefree-actualD source budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_20261107.json';old=json.loads(previous.read_text())['producer_bytes']
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old);after=functions(Path(__file__).read_text());assert len(before)==29
    assert all(after[name]==value for name,value in before.items())
    x=s.symbols('A B u v w');A,B,u,v,w=x
    models=[('J0_secant_zero',A*u**4+B*v**4,17,[[-4*A,0,u,0,0],[0,-4*B,0,v,0]]),
            ('J0_secant_one',A*u**4+B*v**4+u**3*v**2,18,[[-8*A,12*B,2*u,-3*v,0]]),
            ('J0_secant_both',A*u**4+B*v**4+u**3*v**2+u**2*v**3,19,[]),
            ('J0_compression',A*u**3*v+B*u*v**3+u**5+v**5,19,[]),
            ('I0_mixed_zero',A*u**4+B*u*v**3+v**5,18,[[-4*A,-B,u,0,0]]),
            ('I0_mixed_nonzero',A*u**4+B*u*v**3+u**3*v**2+v**5,19,[])]
    records=[]
    for name,F,Ddim,stabilizers in models:
        M,basis,mon6=full_quadratic_casimir_matrix(F,x);cert,_,K=rank_certificate(M,False)
        for Z in stabilizers:assert s.expand(sum(z*s.diff(F,t)for z,t in zip(Z,x)))==0
        cc=s.symbols('k0:'+str(K.cols));values=K*s.Matrix(cc)if K.cols else s.zeros(150,1);P=s.zeros(5)
        for r,(i,j,mon)in enumerate(basis):
            if values[r]!=0:P[i,j]+=values[r]*s.prod(t**e for t,e in zip(x,mon));P[j,i]=-P[i,j]
        V=[]
        for i in range(5):
            j,k,l,m=[r for r in range(5)if r!=i]
            V.append(s.expand((-1)**i*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])))
        gradient=[s.diff(F,t)for t in x];mon,coefficient=s.Poly(gradient[0],*x).terms()[0]
        q=s.expand(s.Poly(V[0],*x).coeff_monomial(mon)/coefficient)
        assert not set(x)&q.free_symbols
        assert all(s.expand(V[i]-q*gradient[i])==0 for i in range(5))
        empty=q==0;upper=Ddim+cert['nullity']-1
        records.append({'type':name,'actual_D':str(F),'full1075by150_Casimir_and_divergence_matrix':serial_matrix(M),
            'full_canonical_source_certificate':cert,'actual_linear_stabilizer_fields':[list(map(str,Z))for Z in stabilizers],
            'universal_source_Pi':{str((i,j)):str(s.expand(P[i,j]))for i,j in H.PAIRS},
            'full_five_Pfaffian_components':list(map(str,V)),'Pfaffian_gradient_scalar_quadratic':str(s.factor(q)),
            'full_Pfaffian_vector_equals_q_gradientF':True,'actual_Euler_Pfaffian_equals_5qF':True,
            'entire_necessary_source_rank_at_most2':empty,'whole_D_family_dimension_upper':Ddim,
            'whole_projective_image_upper':upper,'G_dimension_criterion':empty or upper<=24})
        print(name,'canonicalsource',cert['nullity'],'whole upper',upper,'rank4empty',empty,flush=True)
    proof=ROOT/'research/classification/PLANE_M4_NONREDUCED_SQUAREFREE_QUINTICS.md'
    pencils=ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_PENCILS.md'
    identity=ROOT/'research/classification/ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md'
    archive=ROOT/'research/classification/inputs/plane_m4_two_line_casimir_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_G24_closed':all(row['G_dimension_criterion']for row in records),
            'scope':'Entire squarefree normal-leading m4 rank2-plane source with nonred zero cubic. Nonexceptionalj has rank1 by rootcover overP2minusdiscriminantline, independent ofother connection poles. Special stablepencils reduce to sixactualD, allfullcanonical modules/Pfaffians computed. Rank1 goes to accepted essential3; all3L rank2 excluded byusedlines.',
            'all150_quadratic_basis':[[i,j,list(mon)]for i,j,mon in basis],'all210_degree6_monomials':mon6,'cases':records,
            'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [HELPER,previous,proof,pencils,identity,archive]],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def double_line_m3_squarefree_main():
    """One whole finite classification batch for all split/1+2/3 root covers."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole nonred m3squarefree batch budget300s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(300)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_CASIMIR_20261107.json'
    old=json.loads(previous.read_text())['producer_bytes']
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old);after=functions(Path(__file__).read_text());assert len(before)==31
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;f=a*a*b;normals=[z,w]
    mons=lambda n,d:[tuple(t.count(i)for i in range(n))for t in __import__('itertools').combinations_with_replacement(range(n),d)]
    monG=[base+(j,3-j)for base in mons(3,2)for j in range(4)]
    monPDE=[base+(j,3-j)for base in mons(3,4)for j in range(4)]
    monD2=[base+(j,2-j)for base in mons(3,3)for j in range(3)]
    spatial=lambda mon:s.prod(t**e for t,e in zip(x,mon))
    P1,v1,P2,v2=sources(x)
    def actualD(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def baseP():
        P=s.zeros(5);P[0,2]=-a*a;P[2,0]=a*a;P[1,2]=2*a*b;P[2,1]=-2*a*b;return P
    def setbr(P,i,j,value):P[i,j]=value;P[j,i]=-value
    pairs=set();intersections=[]
    for k,l in [(0,2),(0,3),(1,2),(1,3)]:
        det=(k-1)*(2-l)-(2-k)*(l-1)
        for p in range(3):
            for q in range(3):
                r0=4-3*p;r1=4-3*q
                lam=s.Rational(r0*(2-l)-(2-k)*r1,det);mu=s.Rational((k-1)*r1-r0*(l-1),det)
                assert (k-1)*lam+(2-k)*mu==r0 and(l-1)*lam+(2-l)*mu==r1
                pairs.add(tuple(sorted([lam,mu])));intersections.append({'normal_support':[k,l],'base_powers':[p,q],'lambda_mu':[str(lam),str(mu)]})
    models=[]
    for lam,mu in sorted(pairs):
        P=baseP();setbr(P,2,3,lam*a*z);setbr(P,2,4,mu*a*w)
        models.append(('DIAGONAL_'+str(lam)+'_'+str(mu),P,5 if lam==mu else 7,'diagonal',{'lambda':str(lam),'mu':str(mu)}))
    for name,lam,mu,Bdim in [('NIL_SPLIT',s.Integer(-2),s.Integer(-2),7),('NIL_ONE_PLUS_TWO',s.Rational(-1,2),s.Integer(-2),8)]:
        P=baseP();setbr(P,2,3,lam*a*z);setbr(P,2,4,mu*a*w+b*z);models.append((name,P,Bdim,'nil',{}))
    for g in [2,3,4]:
        h=s.Rational(g,2);P=baseP();setbr(P,2,3,(4-3*h)*a*z);setbr(P,0,4,(3-h)*a*z)
        setbr(P,1,4,-h*b*z);setbr(P,2,4,(g-3)*c*z-2*a*w)
        models.append(('C_DEPENDENT_A2_G'+str(g),P,10,'cdependent',{'g':g}))
    P=baseP();setbr(P,2,3,a*z);setbr(P,1,4,-3*a*z);setbr(P,2,4,4*a*w)
    models.append(('C_DEPENDENT_B2_CANDIDATE',P,10,'cdependent_inconsistent',{}))
    records=[];kk=s.symbols('qK0:3');gg=s.symbols('G0:24');G=sum(t*spatial(mon)for t,mon in zip(gg,monG))
    for name,P0,Bdim,kind,metadata in models:
        Pk=P0.copy();setbr(Pk,3,4,kk[0]*z*z+kk[1]*z*w+kk[2]*w*w)
        Dk=actualD(Pk);poly=s.Poly(Dk,*x);rows=s.Matrix([poly.coeff_monomial(mon)for mon in monD2])
        assert s.expand(Dk-sum(rows[i]*spatial(mon)for i,mon in enumerate(monD2)))==0
        KM=rows.jacobian(kk);constant=rows.subs(dict.fromkeys(kk,0));assert rows==KM*s.Matrix(kk)+constant
        kc,_,_=rank_certificate(KM,False);ka,_,_=rank_certificate(KM.row_join(-constant),False)
        record={'type':name,'kind':kind,'metadata':metadata,'whole_fixed_plane_symbol_dimension_upper':Bdim,
            'full30_K_global_consistency_matrix':serial_matrix(KM),'full30_K_global_consistency_constant':serial_matrix(constant),
            'K_matrix_certificate':kc,'K_augmented_certificate':ka}
        if ka['rank']>kc['rank']:
            record.update(global_quadratic_P0_impossible=True,entire_squarefree_m3_source_empty=True,G_dimension_criterion=True)
            records.append(record);print(name,'whole globalK inconsistent',flush=True);continue
        solution=list(s.linsolve((KM,-constant),kk));assert len(solution)==1
        Kvalue=s.expand(sum(t*m for t,m in zip(solution[0],[z*z,z*w,w*w])));setbr(P0,3,4,Kvalue)
        assert actualD(P0)==0
        Jac0=[s.expand(t)for t in H.all350(H.projective(H.bilinear(P0,P0,x),x),x)]
        if any(t!=0 for t in Jac0):
            record.update(global_quadratic_P0_impossible=True,unique_K_solution=list(map(str,solution[0])),
                nonzero_full350_P0_Jacobi=[(i,str(t))for i,t in enumerate(Jac0)if t!=0],
                entire_squarefree_m3_source_empty=True,G_dimension_criterion=True)
            records.append(record);print(name,'whole P0Jacobi inconsistent',flush=True);continue
        S=[[[0 for i in range(3)]for j in range(2)]for k in range(2)]
        for i in range(2):
            for j in range(2):
                V=s.Matrix([-s.diff(P0[k,3+i],normals[j])for k in range(3)])
                assert s.expand(sum(s.diff(V[k],x[k])for k in range(3)))==0
                numerator=[c*V[1]-b*V[2],a*V[2]-c*V[0],b*V[0]-a*V[1]]
                S[i][j]=[s.expand(numerator[k]+(s.diff(f,x[k])if i==j else 0))for k in range(3)]
        expressions=[]
        for i in range(3):
            expr=3*f*s.diff(G,x[i])-3*s.diff(f,x[i])*G
            expr+=sum(S[j][k][i]*normals[k]*s.diff(G,normals[j])for j in range(2)for k in range(2))
            expr-=(S[0][0][i]+S[1][1][i])*G;expr=s.expand(expr);p=s.Poly(expr,*x)
            expressions.extend(p.coeff_monomial(mon)for mon in monPDE)
        PM=s.Matrix(expressions).jacobian(gg);assert s.Matrix(expressions)==PM*s.Matrix(gg)
        pc,_,KG=rank_certificate(PM,False);assert KG.cols>0
        Ann=s.Matrix.hstack(*KG.T.nullspace()).T;assert Ann*KG==s.zeros(Ann.rows,KG.cols)
        D1=operator(P0,P1,v1,x);D2=operator(P0,P2,v2,x)
        full=actualD(P0+P1+P2);p=s.Poly(full,*x);D3=s.Matrix([p.coeff_monomial(mon)for mon in monG]);L=D3.jacobian(v1)
        assert D3==L*s.Matrix(v1)
        M1=D1.col_join(Ann*L);c1,_,K1=rank_certificate(M1,False);c2,_,K2=rank_certificate(D2,False)
        Sh=shear(P0,x);assert all(s.expand(t)==0 for t in M1*Sh)
        cs,_,_=rank_certificate(Sh,False);image=L*K1;ic,_,_=rank_certificate(image,False)
        cols=ic['minor_columns'];tt=s.symbols('I0:'+str(len(cols)))
        gim=s.expand(sum(image[i,j]*tt[k]*spatial(monG[i])for k,j in enumerate(cols)for i in range(24)))
        gp=s.Poly(gim,z,w);A3=gp.coeff_monomial(z**3);B3=gp.coeff_monomial(z*z*w);C3=gp.coeff_monomial(z*w*w);D3n=gp.coeff_monomial(w**3)
        disc=s.expand(B3**2*C3**2-4*A3*C3**3-4*B3**3*D3n-27*A3**2*D3n**2+18*A3*B3*C3*D3n)
        empty=disc==0;higher=c1['nullity']+c2['nullity'];upper=Bdim+5+higher
        record.update(unique_global_K=str(Kvalue),full150_P0={str((i,j)):str(P0[i,j])for i,j in H.PAIRS},P0_full350_Jacobi_zero_and_actualD_zero=True,
            full_3f_connection_numerators=[[[str(v)for v in t]for t in row]for row in S],
            full180by24_parallel_matrix=serial_matrix(PM),parallel_certificate=pc,parallel_membership_annihilator=serial_matrix(Ann),
            full350_weight1_operator=serial_matrix(D1),full350_weight2_operator=serial_matrix(D2),full24by30_actualD3_map=serial_matrix(L),
            full_weight1_plus_parallel_matrix=serial_matrix(M1),weight1_certificate=c1,weight2_certificate=c2,
            complete_actual_shear=serial_matrix(Sh),actual_shear_certificate=cs,full_leading_image_matrix=serial_matrix(image),leading_image_certificate=ic,
            whole_leading_image_polynomial=str(gim),whole_leading_discriminant_polynomial=str(disc),
            entire_squarefree_m3_source_empty=empty,whole_higher_source_dimension_upper=higher,whole_projective_image_upper=upper,
            G_dimension_criterion=empty or upper<=24)
        records.append(record);print(name,'P1/P2',c1['nullity'],c2['nullity'],'image',ic['rank'],'SFempty',empty,'upper',upper,flush=True)
    proofs=[ROOT/'research/classification'/name for name in ['PLANE_DOUBLE_LINE_M3_SQUAREFREE_FINITE_CONNECTIONS.md','PLANE_DOUBLE_LINE_M3_CYCLIC_CUBIC_FIRST_JETS.md','PLANE_DOUBLE_LINE_M3_CANDIDATE_LIFTS.md']]
    archive=ROOT/'research/classification/inputs/plane_double_line_m3_pre_batch_frozen_inputs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_G24_closed':all(t['G_dimension_criterion']for t in records),
            'all_squarefree_m3_sources_empty':all(t['entire_squarefree_m3_source_empty']for t in records),
            'scope':'Entire R=2L+M rank2-plane multiplicity3 squarefree normal cubic, allsplit/1+2/connected3 covers. Wholefinite root/Riccati/degree classification, globalK consistency, allparallel spaces, complete30/9 and actualD3 images. R3L separate alreadyempty; m2/repeatedm3/reducedR remain outside.',
            'all36_eigenweight_intersections':intersections,'all_normal_swap_distinct_diagonal_pairs':[[str(t)for t in pair]for pair in sorted(pairs)],
            'all24_leading_monomials':monG,'all30_K_monomials':monD2,'all60_parallel_component_monomials':monPDE,
            'full_source1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},'full_source2':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},
            'cases':records,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
            'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
            'input_bindings':[{'path':str(t.relative_to(ROOT)),'sha256':hashlib.sha256(t.read_bytes()).hexdigest()}for t in [HELPER,previous,archive]+proofs],
            'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def double_line_m3_ore_main():
    """Whole remaining two sources: exact Euler bridge and four affine modules."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_ORE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole4 Ore-module budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_REPLAY_20261107.json';packet=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(packet['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==32
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;y=s.Matrix([a,b,z,w]);indices=[0,1,3,4]
    P1,v1,P2,v2=sources(x);ww=s.symbols('W0:18');qv=s.symbols('q0:3');F=sum(t*m for t,m in zip(qv,[z*z,z*w,w*w]))
    def Wbasis(W0):
        return W0+s.Matrix([z*(ww[0]*a+ww[1]*b)+w*(ww[2]*a+ww[3]*b)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*a+ww[8]*b)+w*(ww[9]*a+ww[10]*b)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w])
    Wd=s.Matrix([a*a,-2*a*b,-2*a*z,-2*a*w]);Wn=Wd+s.Matrix([0,0,0,b*z])
    def QD(A,Xz,Xw,F):
        Q=s.zeros(4)
        for i,j,value in [(0,1,-b*A+F),(0,2,-A*z),(0,3,-A*w),(1,2,Xz),(1,3,Xw)]:Q[i,j]=value;Q[j,i]=-value
        return Q
    def QN(U,V,R,F):
        Q=s.zeros(4)
        for i,j,value in [(0,1,b*U*z+F),(0,2,2*U*z*z),(0,3,3*U*z*w-V*z*z/6),
                          (1,2,V*z*z-6*U*z*w),(1,3,R*z*z+2*V*z*w-9*U*w*w)]:Q[i,j]=value;Q[j,i]=-value
        return Q
    def pad(Q,W):
        P=s.zeros(5)
        for i,j in combinations(range(4),2):P[indices[i],indices[j]]=Q[i,j];P[indices[j],indices[i]]=-Q[i,j]
        for i in range(4):P[2,indices[i]]=W[i];P[indices[i],2]=-W[i]
        return P
    def actualD(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def jac4(Q):
        return {(i,j,k):s.expand(sum(Q[i,l]*s.diff(Q[j,k],y[l])+Q[j,l]*s.diff(Q[k,i],y[l])+Q[k,l]*s.diff(Q[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
    coverage=[]
    for kind,name,W0 in [('D','DIAGONAL_-2_-2',Wd),('N','NIL_SPLIT',Wn)]:
        old=next(row for row in packet['cases']if row['type']==name);K=s.Matrix(old['weight1_certificate']['complete_kernel'])
        kk=s.symbols('source'+kind+'0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=s.sympify(old['full150_P0'][str((i,j))],locals={str(t):t for t in x});P0[j,i]=-P0[i,j]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(kk))),simultaneous=True)+P2
        L=s.Matrix([-s.diff(P[0,2],c),-s.diff(P[1,2],c),0,0,0]);gauged=(P-H.wedge(x,L)).applyfunc(s.expand)
        assert all(s.diff(gauged[i,j],c)==0 for i,j in H.PAIRS)
        Q=gauged.extract(indices,indices);W=s.Matrix([gauged[2,i]for i in indices]);central=P2[0,1]
        if kind=='D':
            A=s.cancel(-Q[0,2]/z);Xz=Q[1,2];Xw=Q[1,3];expected=QD(A,Xz,Xw,central)
            coordinates=[s.expand(A).coeff(t)for t in [z,w]]+[s.Poly(part,z,w).coeff_monomial(mon)for part in [Xz,Xw]for mon in [z*z,z*w,w*w]]
        else:
            U=s.Poly(Q[0,2],z,w).coeff_monomial(z*z)/2;V=-6*s.Poly(Q[0,3],z,w).coeff_monomial(z*z);R=s.Poly(Q[1,3],z,w).coeff_monomial(z*z)
            expected=QN(U,V,R,central);coordinates=[U,V,R]
        assert all(s.expand(t)==0 for t in Q-expected)
        coordinates+=[s.Poly(central,z,w).coeff_monomial(mon)for mon in [z*z,z*w,w*w]]
        wd=W-W0;wcoords=[]
        for i in range(2):wcoords.extend([s.Poly(wd[i],a,b,z,w).coeff_monomial(mon)for mon in [a*z,b*z,a*w,b*w,z*z,z*w,w*w]])
        wcoords+=[s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),
                  s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)]
        assert all(s.expand(t)==0 for t in W-Wbasis(W0).subs(dict(zip(ww,wcoords)),simultaneous=True))
        variables=list(kk)+v2;J=s.Matrix(coordinates+wcoords).jacobian(variables);det=J.det(method='domain-ge')
        assert J.rows==J.cols and det!=0 and not s.sympify(det).free_symbols
        coverage.append({'type':kind,'full_D1_kernel':serial_matrix(K),'literal_Euler_gauge_L':list(map(str,L)),
            'all_full_source_c_independent':True,'full_Q':{str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},
            'full_W':list(map(str,W)),'complete_joint_source_change_matrix':serial_matrix(J),'nonzero_constant_determinant':str(det)})
    ad=s.symbols('A0:2');xd=s.symbols('X0:6');A=ad[0]*z+ad[1]*w;Xz=xd[0]*z*z+xd[1]*z*w+xd[2]*w*w;Xw=xd[3]*z*z+xd[4]*z*w+xd[5]*w*w
    JD=jac4(QD(A,Xz,Xw,F));XA=s.expand(Xz*s.diff(A,z)+Xw*s.diff(A,w))
    assert JD[(0,1,2)]==s.expand(z*XA)and JD[(0,1,3)]==s.expand(w*XA)and JD[(0,2,3)]==0 and JD[(1,2,3)]==0
    nu,nv,nr=s.symbols('nil_u nil_v nil_r');QNgen=QN(nu,nv,nr,F);assert all(t==0 for t in jac4(QNgen).values())
    gp=s.Poly(actualD(pad(QN(nu,nv,nr,0),Wn)),z,w);G3=s.expand(sum(coef*z**power[0]*w**power[1]for power,coef in gp.terms()if sum(power)==3))
    hand=z*(3*nr*a*a*z*z+nv*a*z*(3*a*w+b*z)-nu*(3*a*w+b*z)**2)
    assert s.expand(G3-hand)==0 or s.expand(G3+hand)==0
    coefficients=[s.Poly(G3,z,w).coeff_monomial(mon)for mon in [z**3,z*z*w,z*w*w,w**3]];AA,BB,CC,DD=coefficients
    disc=s.expand(BB*BB*CC*CC-4*AA*CC**3-4*BB**3*DD-27*AA*AA*DD*DD+18*AA*BB*CC*DD)
    assert s.expand(disc-729*nu*nu*(nv*nv+12*nu*nr)*a**8)==0
    h=s.Symbol('normal_shear');M=s.eye(4);M[3,2]=h;new= M*QNgen.subs(w,w-h*z,simultaneous=True)*M.T
    predicted=QN(nu,nv+6*h*nu,nr-h*nv-3*h*h*nu,F.subs(w,w-h*z))
    assert all(s.expand(t)==0 for t in new-predicted)
    assert s.expand((nv+6*h*nu)**2+12*nu*(nr-h*nv-3*h*h*nu)-nv*nv-12*nu*nr)==0
    uz,uw,vz,vw=s.symbols('shift_uz shift_uw shift_vz shift_vw');shiftU=uz*z+uw*w;shiftV=vz*z+vw*w;M=s.eye(4);M[0,2]=uz;M[0,3]=uw;M[1,2]=vz;M[1,3]=vw
    R2=w*(w-z);QDcan=QD(z,0,R2,F);QNcan=QN(1,0,1,F)
    Dpush=M*QDcan.subs({a:a-shiftU,b:b-shiftV},simultaneous=True)*M.T
    Npush=M*QNcan.subs({a:a-shiftU,b:b-shiftV},simultaneous=True)*M.T
    Nchange=(vz-uw)*z*z+(2*vw+6*uz)*z*w+9*uw*w*w
    assert all(s.expand(t)==0 for t in Dpush-QD(z,0,R2,F-uw*R2))
    assert all(s.expand(t)==0 for t in Npush-QN(1,0,1,F+Nchange))
    scaleS,scaleR,scaleK,scaleL=s.symbols('scale_s scale_r scale_k scale_l',nonzero=True)
    def push(P,coords,M,scalar):
        old=M.inv()*s.Matrix(coords);return (scalar*M*P.subs(dict(zip(coords,old)),simultaneous=True)*M.T).applyfunc(s.cancel)
    kscale=scaleS**2/scaleR;scalar=scaleS/scaleR;M5=s.diag(scaleS,scaleS,scaleR,kscale,kscale);M4=s.diag(scaleS,scaleS,kscale,kscale)
    assert all(t==0 for t in push(pad(s.zeros(4),Wd),x,M5,scalar)-pad(s.zeros(4),Wd))
    assert all(s.cancel(t)==0 for t in push(QDcan,y,M4,scalar)-QD(z,0,R2,scaleR*F/scaleS))
    bscale=scaleS*scaleL/scaleK;M5=s.diag(scaleS,bscale,scaleR,scaleK,scaleL);M4=s.diag(scaleS,bscale,scaleK,scaleL)
    assert all(t==0 for t in push(pad(s.zeros(4),Wn),x,M5,scalar)-pad(s.zeros(4),Wn))
    up=scaleS**2*nu/(scaleR*scaleK);vp=scaleS**2*scaleL*nv/(scaleR*scaleK**2);rp=scaleS**2*scaleL**2*nr/(scaleR*scaleK**3)
    fp=scaleS*bscale*scalar*F.subs({z:z/scaleK,w:w/scaleL},simultaneous=True)
    assert all(s.cancel(t)==0 for t in push(QNgen,y,M4,scalar)-QN(up,vp,rp,fp))
    def cert(M,domain):
        dm=DomainMatrix.from_Matrix(M).convert_to(domain);rr,pivots=dm.rref();rows=dm.transpose().rref()[1];R=rr.to_Matrix()
        free=[i for i in range(M.cols)if i not in pivots];K=s.zeros(M.cols,len(free))
        for j,column in enumerate(free):
            K[column,j]=1
            for i,pivot in enumerate(pivots):K[pivot,j]=-R[i,column]
        assert (dm*DomainMatrix.from_Matrix(K).convert_to(domain)).is_zero_matrix
        determinant=dm.extract(list(rows),list(pivots)).det();assert determinant!=domain.zero or not pivots
        C={'rank':len(pivots),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(pivots),
           'minor_determinant':str(domain.to_sympy(determinant)),'complete_kernel':serial_matrix(K),'kernel_substitution_zero':True}
        return C,dm,K,domain.to_sympy(determinant)
    def analyze(M,C,MM,MC,Dconstant,Dlinear,domain,offset):
        cc,dm,K,minor=cert(M,domain);ac,ad,_,aminor=cert(M.row_join(C),domain)
        m140=DomainMatrix.from_Matrix(MM.row_join(MC)).convert_to(domain);rr1=ad.rref()[0].to_Matrix();rr2=m140.rref()[0].to_Matrix()
        assert len(m140.rref()[1])==ac['rank']and all(domain.from_sympy(t)==domain.zero for t in rr1[:ac['rank'],:]-rr2[:ac['rank'],:])
        ri=cc['minor_rows'];ci=cc['minor_columns'];inverse=dm.extract(ri,ci).inv();cdm=DomainMatrix.from_Matrix(C).convert_to(domain)
        solved=(-inverse*cdm.extract(ri,[0])).to_Matrix();part=s.zeros(18,1)
        for i,j in enumerate(ci):part[j]=solved[i]
        residual=(dm*DomainMatrix.from_Matrix(part).convert_to(domain)+cdm).to_Matrix()
        record={'full350_affine_matrix':serial_matrix(dm.to_Matrix()),'full350_affine_constant':list(map(str,cdm.to_Matrix())),
            'full140_module_matrix':serial_matrix(DomainMatrix.from_Matrix(MM).convert_to(domain).to_Matrix()),
            'full140_module_constant':list(map(str,DomainMatrix.from_Matrix(MC).convert_to(domain).to_Matrix())),
            'full350_and140_same_affine_rowspace':True,'linear_certificate':cc,'augmented_certificate':ac}
        denominators=[]
        if ac['rank']>cc['rank']:
            at=next(i for i,t in enumerate(residual)if t!=0);weights=(dm.extract([at],ci)*inverse).to_Matrix();phi=s.zeros(1,M.rows);phi[0,at]=1
            for j,i in enumerate(ri):phi[0,i]-=weights[j]
            pd=DomainMatrix.from_Matrix(phi).convert_to(domain);assert(pd*dm).is_zero_matrix
            contradiction=(pd*cdm).to_Matrix()[0];assert contradiction!=0
            record.update(status='whole_affine_module_inconsistent',cokernel_witness=[[i,str(t)]for i,t in enumerate(phi)if t!=0],
                nonzero_constant_contradiction=str(contradiction),G_empty=True,G_dimension_criterion=True)
        else:
            assert all(t==0 for t in residual)
            D0=s.Poly(s.expand(Dconstant+sum(Dlinear[i]*part[i]for i in range(18))),*x,domain=domain)
            directions=[s.Poly(s.expand(sum(Dlinear[i]*K[i,j]for i in range(18))),*x,domain=domain)for j in range(K.cols)]
            common=D0
            for poly in directions:
                if not poly.is_zero:common=s.polys.polytools.gcd(common,poly).monic()
            factor=None
            if 0<common.total_degree()<5:factor=common.monic()
            elif common.total_degree()==5:
                for poly,power in common.factor_list()[1]:
                    if 0<poly.total_degree()<5:factor=poly.monic();break
            if factor is not None and domain==FIELD:
                for poly in [factor]+[poly.exquo(factor)for poly in [D0]+directions]:
                    for value in poly.coeffs():denominators.append(s.denom(s.cancel(value)))
            upper=offset+K.cols
            record.update(status='complete_affine_module',particular_solution=list(map(str,part)),full_W_fiber_dimension=K.cols,
                actual_D_particular=str(D0.as_expr()),actual_D_kernel_coefficients=[str(poly.as_expr())for poly in directions],
                common_D_polynomial_gcd=str(common.as_expr()),proper_fixed_D_factor=None if factor is None else str(factor.as_expr()),
                whole_projective_image_upper=upper,G_empty=factor is not None,G_dimension_criterion=factor is not None or upper<=24)
        return record,minor,aminor,denominators
    cases=[]
    for name,Q,W0,offset,variable in [('D_ZERO',QD(z,0,R2,0),Wd,16,False),('D_AXIS',QD(z,0,R2,z*w),Wd,17,False),
                                     ('D_LINE',QD(z,0,R2,z*z+TAU*z*w),Wd,18,True),('N_FIXED',QN(1,0,1,0),Wn,18,False)]:
        assert all(t==0 for t in jac4(Q).values());W=Wbasis(W0);P=pad(Q,W)
        rows=s.Matrix([s.expand(t)for t in H.all350(H.projective(H.bilinear(P,P,x),x),x)]);C=rows.subs(dict.fromkeys(ww,0));M=rows.jacobian(ww)
        assert rows==C+M*s.Matrix(ww)
        LW=s.zeros(4)
        for i,j in combinations(range(4),2):
            LW[i,j]=s.expand(sum(W[k]*s.diff(Q[i,j],y[k])-Q[k,j]*s.diff(W[i],y[k])-Q[i,k]*s.diff(W[j],y[k])for k in range(4)));LW[j,i]=-LW[i,j]
        module=[]
        for i,j,k in combinations(range(4),3):
            poly=s.Poly(s.expand(y[i]*LW[j,k]-y[j]*LW[i,k]+y[k]*LW[i,j]),*y)
            module.extend(poly.coeff_monomial(mon)for mon in H.monomials(4,4))
        module=s.Matrix(module);MC=module.subs(dict.fromkeys(ww,0));MM=module.jacobian(ww);assert module==MC+MM*s.Matrix(ww)
        D=actualD(P);Dconstant=D.subs(dict.fromkeys(ww,0));Dlinear=[s.expand(D).coeff(t)for t in ww]
        assert s.expand(D-Dconstant-sum(value*t for value,t in zip(Dlinear,ww)))==0
        record,minor,aminor,denominators=analyze(M,C,MM,MC,Dconstant,Dlinear,FIELD if variable else s.QQ,offset)
        record.update(type=name,whole_Q={str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},W0=list(map(str,W0)),whole_parameter_tau=variable)
        if variable:
            exceptional=s.Poly(s.cancel(minor*aminor),TAU)
            for value in denominators:exceptional=s.lcm(exceptional,s.Poly(value,TAU))
            factors=s.factor_list(exceptional.as_expr())[1];walls=[]
            for polynomial,power in factors:
                p=s.Poly(polynomial,TAU)
                if p.degree()==1:theta=-p.nth(0)/p.nth(1);domain=s.QQ
                else:theta=s.CRootOf(p.as_expr(),0);domain=s.QQ.alg_field(theta)
                sub=lambda expr:expr.subs(TAU,theta)
                result,_,_,_=analyze(sub(M),sub(C),sub(MM),sub(MC),sub(Dconstant),[sub(t)for t in Dlinear],domain,17)
                result.update(minimal_polynomial=str(polynomial),tau=str(theta),all_complex_embeddings_covered=True,minor_factor_multiplicity=int(power))
                walls.append(result)
            record.update(all_exception_polynomial=str(s.factor(exceptional.as_expr())),all_exception_factors=[[str(t),int(e)]for t,e in factors],
                factor_specialization_denominators=list(map(str,denominators)),all_exceptional_fields=walls,
                whole_line_G_closed=record['G_dimension_criterion']and all(t['G_dimension_criterion']for t in walls))
        else:record['whole_line_G_closed']=record['G_dimension_criterion']
        cases.append(record);print(name,record['status'],'k',record['linear_certificate']['nullity'],'Gempty',record['G_empty'],'Gclosed',record['whole_line_G_closed'],flush=True)
    proof=ROOT/'research/classification/PLANE_DOUBLE_LINE_M3_ORE_REDUCTION.md';saito=ROOT/'research/classification/ESSENTIAL_THREE_QUINTIC_EIGHT_COMPONENTS.md'
    archive=ROOT/'research/classification/inputs/plane_double_line_m3_ore_preappend_inputs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_G24_closed':all(t['whole_line_G_closed']for t in cases),
        'full_two_source_Euler_Ore_coverage':coverage,'strict_D_Q_Jacobi_equals_z_w_times_XA':{'full_Jacobi':{str(k):str(v)for k,v in JD.items()},'three_generators':[str(s.Poly(XA,z,w).coeff_monomial(mon))for mon in [z*z,z*w,w*w]]},
        'strict_N_Q_Jacobi_identically_zero':True,'N_squarefree_discriminant_exact':str(s.factor(disc)),
        'all_normal_shear_coboundary_diagonal_identities_zero':True,'N_central_coboundary':str(Nchange),'D_central_coboundary':str(-uw*R2),
        'cases':cases,'scope':'Entire remaining DIAGscalar/NIL_SPLIT m3squarefree source. DIAG A0 actualDessential3; A!=0 andN SFopen reduce byactualnormalizers to twofixedD,wholeDline andonefixedN. Fulloriginal350/140 affineW18, alltauwalls andfactor-specialization denominators retained. No jet/GB/pointfarm.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(t.relative_to(ROOT)),'sha256':hashlib.sha256(t.read_bytes()).hexdigest()}for t in [HELPER,previous,proof,saito,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    try:
        if '--double-line-m3-ore' in __import__('sys').argv:double_line_m3_ore_main()
        elif '--double-line-m3-squarefree' in __import__('sys').argv:double_line_m3_squarefree_main()
        elif '--m4-nonreduced-squarefree' in __import__('sys').argv:m4_nonreduced_squarefree_main()
        elif '--m4-two-line-casimir' in __import__('sys').argv:m4_two_line_casimir_main()
        elif '--nonreduced-cubic-leading' in __import__('sys').argv:nonreduced_cubic_leading_main()
        elif '--nonreduced-cubic' in __import__('sys').argv:nonreduced_cubic_main()
        elif '--global-ore' in __import__('sys').argv:global_ore_main()
        elif '--global-dependent' in __import__('sys').argv:global_dependent_main()
        elif '--global-commonimage' in __import__('sys').argv:global_commonimage_main()
        elif '--global-cartan-surface' in __import__('sys').argv:global_cartan_surface_main()
        elif '--global-threeparameter' in __import__('sys').argv:global_threeparameter_main()
        elif '--global-nil-line' in __import__('sys').argv:global_nil_line_main()
        elif '--global-hyperbola' in __import__('sys').argv:global_hyperbola_main()
        elif '--global-twoparameter' in __import__('sys').argv:global_twoparameter_main()
        elif '--global-tsquare' in __import__('sys').argv:global_tsquare_main()
        elif '--global-residual' in __import__('sys').argv:global_residual_main()
        elif '--global-compatibility' in __import__('sys').argv:global_compatibility_main()
        elif '--global-onecut' in __import__('sys').argv:global_onecut_main()
        elif '--global-lowparameter' in __import__('sys').argv:global_lowparameter_main()
        elif '--rank2-nonkernel' in __import__('sys').argv:rank2_nonkernel_main()
        elif '--rank2-kernel' in __import__('sys').argv:rank2_kernel_main()
        else:main()
    except FileExistsError:raise
    except Exception as error:
        failure=ROOT/('reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_ORE_RUNTIME_20261107.json' if '--double-line-m3-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_RUNTIME_20261107.json' if '--double-line-m3-squarefree' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_RUNTIME_20261107.json' if '--m4-nonreduced-squarefree' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_RUNTIME_20261107.json' if '--m4-two-line-casimir' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_LEADING_RUNTIME_20261107.json' if '--nonreduced-cubic-leading' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_RUNTIME_20261107.json' if '--nonreduced-cubic' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ORE_RUNTIME_20261107.json' if '--global-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_DEPENDENT_RUNTIME_20261107.json' if '--global-dependent' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMMONIMAGE_RUNTIME_20261107.json' if '--global-commonimage' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_CARTAN_SURFACE_RUNTIME_20261107.json' if '--global-cartan-surface' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_THREEPARAMETER_RUNTIME_20261107.json' if '--global-threeparameter' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_NIL_LINE_RUNTIME_20261107.json' if '--global-nil-line' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_HYPERBOLA_RUNTIME_20261107.json' if '--global-hyperbola' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TWOPARAMETER_RUNTIME_20261107.json' if '--global-twoparameter' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TSQUARE_RUNTIME_20261107.json' if '--global-tsquare' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_RESIDUAL_RUNTIME_20261107.json' if '--global-residual' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMPATIBILITY_RUNTIME_20261107.json' if '--global-compatibility' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ONECUT_RUNTIME_20261107.json' if '--global-onecut' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_RUNTIME_20261107.json' if '--global-lowparameter' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-nonkernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-kernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_RUNTIME_20261107.json')
        if not failure.exists():
            failure.write_text(json.dumps({'task_id':'G','all_passed':False,'failure':repr(error),
                'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'producer_bytes':Path(__file__).read_text(),'scope':'Stopped bounded symbolic linear-operator job; no mathematical G conclusion.'},indent=2)+'\n')
        raise
