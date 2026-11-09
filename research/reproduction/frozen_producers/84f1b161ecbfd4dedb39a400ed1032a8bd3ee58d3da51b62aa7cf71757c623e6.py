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


if __name__=='__main__':
    try:
        if '--global-tsquare' in __import__('sys').argv:global_tsquare_main()
        elif '--global-residual' in __import__('sys').argv:global_residual_main()
        elif '--global-compatibility' in __import__('sys').argv:global_compatibility_main()
        elif '--global-onecut' in __import__('sys').argv:global_onecut_main()
        elif '--global-lowparameter' in __import__('sys').argv:global_lowparameter_main()
        elif '--rank2-nonkernel' in __import__('sys').argv:rank2_nonkernel_main()
        elif '--rank2-kernel' in __import__('sys').argv:rank2_kernel_main()
        else:main()
    except FileExistsError:raise
    except Exception as error:
        failure=ROOT/('reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TSQUARE_RUNTIME_20261107.json' if '--global-tsquare' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_RESIDUAL_RUNTIME_20261107.json' if '--global-residual' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMPATIBILITY_RUNTIME_20261107.json' if '--global-compatibility' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ONECUT_RUNTIME_20261107.json' if '--global-onecut' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_RUNTIME_20261107.json' if '--global-lowparameter' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-nonkernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-kernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_RUNTIME_20261107.json')
        if not failure.exists():
            failure.write_text(json.dumps({'task_id':'G','all_passed':False,'failure':repr(error),
                'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'producer_bytes':Path(__file__).read_text(),'scope':'Stopped bounded symbolic linear-operator job; no mathematical G conclusion.'},indent=2)+'\n')
        raise
