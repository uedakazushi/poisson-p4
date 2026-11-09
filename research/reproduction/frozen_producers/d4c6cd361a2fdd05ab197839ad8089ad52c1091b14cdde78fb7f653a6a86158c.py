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


def nonreduced_m3_constant_roots_main():
    """All14 whole constant double/simple-root symbols, exact leading line first."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_ROOTS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole14 constantroot finite budget300s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(300)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_ORE_20261107.json';old=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==33
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;P1,v1,P2,v2=sources(x)
    monG=[base+(j,3-j)for base in H.monomials(3,2)for j in range(4)]
    monD2=[base+(j,2-j)for base in H.monomials(3,3)for j in range(3)]
    kk=s.symbols('K0:3');spatial=lambda mon:s.prod(t**e for t,e in zip(x,mon))
    def setbr(P,i,j,value):P[i,j]=value;P[j,i]=-value
    def actualD(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(t)for t in row]for row in M.tolist()],M.shape,domain)
        rr,columns=dm.rref();rows=dm.transpose().rref()[1];free=[i for i in range(M.cols)if i not in columns]
        kernel=[[domain.zero for j in free]for i in range(M.cols)];reduced=rr.to_list()
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(columns):kernel[pivot][j]=-reduced[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in (dm*kd).to_list())
        minor=dm.extract(list(rows),list(columns)).det();assert minor!=domain.zero or not columns
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(columns),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
                'minor_determinant':str(domain.to_sympy(minor)),'complete_kernel':serial_matrix(K),'kernel_substitution_zero':True},domain.to_sympy(minor),K
    def analyze(M1,D2,Sh,L,target_index,Bdim,domain):
        c1,m1,K1=certify(M1,domain);c2,m2,K2=certify(D2,domain);cs,ms,KS=certify(Sh,domain)
        image=(L*K1).applyfunc(s.cancel)
        assert all(domain.from_sympy(image[i,j])==domain.zero for i in range(image.rows)if i!=target_index for j in range(image.cols))
        coefficients=[domain.from_sympy(image[target_index,j])for j in range(image.cols)]
        nonzero=[i for i,value in enumerate(coefficients)if value!=domain.zero]
        image_rank=int(bool(nonzero));empty=not nonzero;upper=Bdim+c1['nullity']+c2['nullity']+5
        record={'full_weight1_plus_line_matrix':serial_matrix(M1),'full350_weight2_matrix':serial_matrix(D2),'full_shear_matrix':serial_matrix(Sh),
            'weight1_certificate':c1,'weight2_certificate':c2,'actual_shear_certificate':cs,
            'actual_D3_complete_kernel_image':serial_matrix(image),'actual_D3_line_image_rank':image_rank,
            'actual_D3_nonzero_certificate':None if empty else {'kernel_column':nonzero[0],'coefficient':str(domain.to_sympy(coefficients[nonzero[0]]))},
            'whole_m3_source_empty':empty,'whole_fixed_plane_symbol_dimension_upper':Bdim,'whole_projective_image_upper':upper,
            'G_dimension_criterion':empty or upper<=24}
        denominators=[]
        for K in [K1,K2,KS]:denominators.extend(s.denom(s.cancel(t))for t in K)
        return record,s.factor(m1*m2*ms),denominators
    models=[]
    for e in range(3):
        for name,gamma,delta,Bdim in [('EXACT_C',0,1,10),('EXACT_B',1,0,9),('LOG',0,0,8)]:
            P=s.zeros(5)
            for i,j,value in [(0,2,-a*a),(1,2,2*a*b),(2,3,(4-3*e)*a*z),(0,4,-delta*a*w),
                              (1,4,2*delta*b*w),(2,4,-delta*c*w+(3*TAU*a-3*gamma*b)*w)]:setbr(P,i,j,value)
            models.append(('DOUBLE_LINE_E'+str(e)+'_'+name,P,a**e*b**(2-e)*z*z*w,Bdim,True,delta*(3*e-4)*z*w))
    for name,V,Bdim in [('ZERO',s.zeros(3,1),5),('LINEAR',s.Matrix([0,0,3*a]),7),('SQUARE',s.Matrix([0,0,6*b]),8),
                        ('SQUARE_PLUS_LINEAR',s.Matrix([0,-3*a,6*b]),9),('RANK_TWO',s.Matrix([0,-3*b,3*c]),10)]:
        P=s.zeros(5);setbr(P,1,2,3*a*a)
        for i in range(3):setbr(P,i,4,-w*V[i])
        models.append(('TRIPLE_LINE_'+name,P,a*a*z*z*w,Bdim,False,s.Integer(0)))
    records=[]
    for name,P0,target,Bdim,variable,expected_K in models:
        domain=FIELD if variable else s.QQ;Pk=P0.copy();setbr(Pk,3,4,kk[0]*z*z+kk[1]*z*w+kk[2]*w*w)
        Dk=actualD(Pk);p=s.Poly(Dk,*x);D2rows=s.Matrix([p.coeff_monomial(mon)for mon in monD2])
        assert s.expand(Dk-sum(D2rows[i]*spatial(mon)for i,mon in enumerate(monD2)))==0
        KM=D2rows.jacobian(kk);C=D2rows.subs(dict.fromkeys(kk,0));kc,_,_=certify(KM,domain);ka,_,_=certify(KM.row_join(-C),domain)
        assert kc['rank']==ka['rank']==3
        sol=list(s.linsolve((KM,-C),kk));assert len(sol)==1
        Kvalue=s.expand(sum(t*m for t,m in zip(sol[0],[z*z,z*w,w*w])));assert s.expand(Kvalue-expected_K)==0
        setbr(P0,3,4,Kvalue);assert actualD(P0)==0
        Jac0=[s.expand(t)for t in H.all350(H.projective(H.bilinear(P0,P0,x),x),x)]
        assert all(t==0 for t in Jac0)
        D1=operator(P0,P1,v1,x);D2=operator(P0,P2,v2,x);full=actualD(P0+P1+P2);p=s.Poly(full,*x)
        leading=s.Matrix([p.coeff_monomial(mon)for mon in monG]);L=leading.jacobian(v1);assert all(s.expand(t)==0 for t in leading-L*s.Matrix(v1))
        target_vector=s.Matrix([s.Poly(target,*x).coeff_monomial(mon)for mon in monG]);nonzero=[i for i,t in enumerate(target_vector)if t!=0]
        assert len(nonzero)==1 and target_vector[nonzero[0]]==1;idx=nonzero[0]
        Ann=s.eye(24).extract([i for i in range(24)if i!=idx],list(range(24)));assert Ann*target_vector==s.zeros(23,1)
        M1=D1.col_join(Ann*L);Sh=shear(P0,x);assert all(s.expand(t)==0 for t in M1*Sh)
        record,minor,denominators=analyze(M1,D2,Sh,L,idx,Bdim,domain)
        record.update(type=name,whole_alpha_line=variable,full150_P0={str((i,j)):str(P0[i,j])for i,j in H.PAIRS},
            full30_K_consistency_matrix=serial_matrix(KM),full30_K_consistency_constant=serial_matrix(C),K_certificate=kc,K_augmented_certificate=ka,
            unique_global_K=str(Kvalue),P0_full350_Jacobi_zero_and_actual_D_zero=True,full350_weight1_operator=serial_matrix(D1),
            full24by30_actualD3_map=serial_matrix(L),exact_nonzero_leading_target=str(target),full_line_annihilator=serial_matrix(Ann))
        if variable:
            exception=s.Poly(minor,TAU)
            for value in denominators:exception=s.lcm(exception,s.Poly(value,TAU))
            factors=s.factor_list(exception.as_expr())[1];walls=[]
            for factor,multiplicity in factors:
                poly=s.Poly(factor,TAU)
                theta=-poly.nth(0)/poly.nth(1)if poly.degree()==1 else s.CRootOf(poly.as_expr(),0)
                dom=s.QQ if poly.degree()==1 else s.QQ.alg_field(theta)
                sub=lambda M:M.subs(TAU,theta)
                wall,_,_=analyze(sub(M1),sub(D2),sub(Sh),sub(L),idx,Bdim-1,dom)
                wall.update(alpha=str(theta),minimal_polynomial=str(factor),all_complex_embeddings_covered=True,multiplicity=int(multiplicity))
                walls.append(wall)
            record.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_irreducible_factors=[[str(t),int(e)]for t,e in factors],
                all_exceptional_roots=walls,generic_kernel_specialization_denominators=list(map(str,denominators)),
                all_strata_G_closed=record['G_dimension_criterion']and all(t['G_dimension_criterion']for t in walls))
        else:record['all_strata_G_closed']=record['G_dimension_criterion']
        records.append(record);print(name,'k',record['weight1_certificate']['nullity'],record['weight2_certificate']['nullity'],
            'D3image',record['actual_D3_line_image_rank'],'upper',record['whole_projective_image_upper'],
            'walls',[(t['alpha'],t['weight1_certificate']['nullity'],t['weight2_certificate']['nullity'],t['actual_D3_line_image_rank'],t['whole_projective_image_upper'])for t in record.get('all_exceptional_roots',[])],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_DOUBLE_SIMPLE_CONSTANT_ROOTS.md'
    archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_constant_roots_preappend_inputs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_G24_closed':all(t['all_strata_G_closed']for t in records),
        'scope':'AllnonreducedR rank2plane multiplicity3 generic double+simple normal roots, bothprimitive rows constant. Whole9QQalpha+5QQ sources, fullnonlog closedforms retained, fullglobalK/full350/30/9/actualD3 exact nonzero line. Other rowdegrees/generictriple/m2 outside.',
        'cases':records,'all24_leading_monomials':monG,'all30_K_monomials':monD2,'full_source1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},
        'full_source2':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(t.relative_to(ROOT)),'sha256':hashlib.sha256(t.read_bytes()).hexdigest()}for t in [HELPER,previous,proof,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_constant_ore_main():
    """Whole six residual sources: source Euler bridge and full quadratic row test."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_ORE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole6 constantroot Ore finite budget300s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_ROOTS_20261107.json';packet=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(packet['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==34
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;y=s.Matrix([a,b,z,w]);indices=[0,1,3,4]
    P1,v1,P2,v2=sources(x);ww=s.symbols('W0:18');ymons=H.monomials(4,2)
    def parse(t):return s.sympify(t,locals={str(u):u for u in [*x,TAU]})
    def matrix(raw):return s.Matrix([[parse(t)for t in row]for row in raw])
    # Keep the literal five Pfaffians explicit; no source quotient in the divisor.
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def Wbasis(W0):
        return W0+s.Matrix([z*(ww[0]*a+ww[1]*b)+w*(ww[2]*a+ww[3]*b)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*a+ww[8]*b)+w*(ww[9]*a+ww[10]*b)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w])
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(t)for t in row]for row in M.tolist()],M.shape,domain)
        rr,columns=dm.rref();rows=dm.transpose().rref()[1];free=[i for i in range(M.cols)if i not in columns]
        kernel=[[domain.zero for j in free]for i in range(M.cols)];reduced=rr.to_list()
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(columns):kernel[pivot][j]=-reduced[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in (dm*kd).to_list())
        minor=dm.extract(list(rows),list(columns)).det();assert minor!=domain.zero or not columns
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(columns),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
            'minor_determinant':str(s.factor(domain.to_sympy(minor))),'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    def sparse(M):return [{str(j):str(s.factor(M[i,j]))for j in range(M.cols)if M[i,j]!=0}for i in range(M.rows)]
    def analyze(old,alpha=None,domain=FIELD):
        variable=alpha is None and old['whole_alpha_line'];sub=lambda expr:expr if alpha is None else expr.subs(TAU,alpha)
        M1=sub(matrix(old['full_weight1_plus_line_matrix']));cc,_,K=certify(M1,domain)
        tt=s.symbols('p0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=sub(parse(old['full150_P0'][str((i,j))]));P0[j,i]=-P0[i,j]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(tt))),simultaneous=True)+P2
        gauge=s.Matrix([-s.diff(P[0,2],c),-s.diff(P[1,2],c),0,0,0]);gauged=(P-H.wedge(x,gauge)).applyfunc(s.expand)
        bad=[(str((i,j)),str(s.diff(gauged[i,j],c)))for i,j in H.PAIRS if s.diff(gauged[i,j],c)!=0]
        record={'alpha':None if alpha is None else str(alpha),'whole_alpha_line':variable,'original_full_M1_certificate':cc,
            'literal_first_Euler_gauge':list(map(str,gauge)),'c_independent_gauge':not bad}
        if bad:record.update(nonzero_c_derivatives=bad,G_closed=False);return record,[]
        W=s.Matrix([gauged[2,i]for i in indices])
        ell=-s.Poly(W[2],z,w).coeff_monomial(z*z)*z-s.Poly(W[3],z,w).coeff_monomial(w*w)*w
        W=(W+ell*y).applyfunc(s.expand);Q=gauged.extract(indices,indices);Q0=Q.subs(dict.fromkeys(list(tt)+v2,0));W0=W.subs(dict.fromkeys(list(tt)+v2,0))
        coords=s.Matrix([s.Poly(Q[i,j]-Q0[i,j],*y).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in ymons]);variables=list(tt)+v2
        J=coords.jacobian(variables);assert all(s.expand(t)==0 for t in coords-J*s.Matrix(variables))
        qc,qminor,_=certify(J,domain);columns=qc['minor_columns'];uu=s.symbols('Q0:'+str(len(columns)));Qcan=Q0.copy()
        for k,(i,j)in enumerate(combinations(range(4),2)):
            Qcan[i,j]+=sum(J[10*k+r,jj]*spatial*svar for jj,svar in zip(columns,uu)for r,mon in enumerate(ymons)for spatial in [s.prod(t**e for t,e in zip(y,mon))]);Qcan[j,i]=-Qcan[i,j]
        wd=W-W0;wcoords=[]
        for i in range(2):wcoords.extend([s.Poly(wd[i],*y).coeff_monomial(mon)for mon in [a*z,b*z,a*w,b*w,z*z,z*w,w*w]])
        wcoords+=[s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)]
        assert all(s.expand(t)==0 for t in W-Wbasis(W0).subs(dict(zip(ww,wcoords)),simultaneous=True))
        JQ={(i,j,k):s.expand(sum(Qcan[i,l]*s.diff(Qcan[j,k],y[l])+Qcan[j,l]*s.diff(Qcan[k,i],y[l])+Qcan[k,l]*s.diff(Qcan[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        Wcan=Wbasis(W0);Pi=s.zeros(5)
        for i,j in combinations(range(4),2):Pi[indices[i],indices[j]]=Qcan[i,j];Pi[indices[j],indices[i]]=-Qcan[i,j]
        for i in range(4):Pi[2,indices[i]]=Wcan[i];Pi[indices[i],2]=-Wcan[i]
        full=[s.expand(t)for t in H.all350(H.projective(H.bilinear(Pi,Pi,x),x),x)];D=divisor(Pi);dp=s.Poly(D,*x)
        allvars=list(uu)+list(ww);jp=[s.Poly(t,*allvars)for t in full];dterms=[(mon,s.Poly(value,*allvars))for mon,value in dp.terms()]
        monomials=sorted(set(mon for poly in jp+[poly for _,poly in dterms]for mon in poly.monoms()),reverse=True)
        M=s.Matrix([[poly.coeff_monomial(mon)for mon in monomials]for poly in jp]);jc,jminor,JK=certify(M,domain)
        Drows=s.Matrix([[poly.coeff_monomial(mon)for mon in monomials]for _,poly in dterms]);factors=[]
        for spatial_index in [0,1,3,4]:
            rows=[i for i,(mon,_)in enumerate(dterms)if mon[spatial_index]==0];T=Drows.extract(rows,list(range(M.cols)))
            residual=(T*JK).applyfunc(s.cancel);zero=all(domain.from_sympy(t)==domain.zero for t in residual)
            factors.append({'candidate_factor':str(x[spatial_index]),'full_D_restriction_rows':rows,'Jacobi_quadratic_rowspan_membership':zero,
                'nonzero_annihilator_residual':None if zero else next(str(t)for t in residual if domain.from_sympy(t)!=domain.zero)})
        record.update(all_full_source_Euler_gauge_covered=True,second_Euler_c_normal_component=str(ell),Q_source_dimension=qc['rank'],
            full60by_source_Q_coefficient_matrix=serial_matrix(J),Q_space_certificate=qc,full_affine_Q={str((i,j)):str(Qcan[i,j])for i,j in combinations(range(4),2)},
            Q0={str((i,j)):str(Q0[i,j])for i,j in combinations(range(4),2)},full_W0=list(map(str,W0)),full_W18_covers_actual_source=True,
            full_strict_Q_Jacobi={str(k):str(v)for k,v in JQ.items()},full_original350_quadratic_coefficient_matrix=sparse(M),quadratic_source_monomials=monomials,
            quadratic_source_variables=list(map(str,allvars)),full350_quadratic_rowspace_certificate=jc,all_actual_D_spatial_monomials=[mon for mon,_ in dterms],
            all_actual_D_coefficient_matrix=sparse(Drows),candidate_linear_factors=factors,proper_fixed_D_factor_proved=next((t['candidate_factor']for t in factors if t['Jacobi_quadratic_rowspan_membership']),None))
        record['G_closed']=record['proper_fixed_D_factor_proved']is not None
        den=[qminor,jminor]
        if variable:
            den.extend(s.denom(s.cancel(t))for t in list(K)+list(J)+list(JK))
        return record,den
    results=[]
    for old in packet['cases']:
        if old['all_strata_G_closed']:continue
        variable=old['whole_alpha_line'];generic,den=analyze(old,domain=FIELD if variable else s.QQ);generic['type']=old['type'];walls=[]
        if variable:
            exception=s.Poly(s.sympify(old['all_exception_polynomial']),TAU)
            if generic['G_closed']:
                for val in den:
                    num,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(num,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            for factor,power in s.factor_list(exception.as_expr())[1]:
                p=s.Poly(factor,TAU);alpha=-p.nth(0)/p.nth(1)if p.degree()==1 else s.CRootOf(p.as_expr(),0)
                domain=s.QQ if p.degree()==1 else s.QQ.alg_field(alpha);wall,_=analyze(old,alpha,domain)
                wall.update(minimal_polynomial=str(factor),all_complex_embeddings_covered=True);walls.append(wall)
            generic.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,
                whole_source_G_closed=generic['G_closed']and all(t['G_closed']for t in walls))
        else:generic['whole_source_G_closed']=generic['G_closed']
        results.append(generic);print(old['type'],'Qdim',generic.get('Q_source_dimension'),'factor',generic.get('proper_fixed_D_factor_proved'),
            'walls',[(t['alpha'],t.get('Q_source_dimension'),t.get('proper_fixed_D_factor_proved'))for t in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_CONSTANT_ROOTS_PARTIAL_BOUND.md'
    archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_constant_roots_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_six_G_closed':all(t['whole_source_G_closed']for t in results),'cases':results,
        'scope':'Whole6residualconstantrootm3sources. EulerfullP1kernel+FULLP2nine; Qwholelinearimage andW18superset; strictQJac andall350quadraticrowspans testactualDproperfactors. Negativecriterionisnotabsenceoffactor/Gclosure.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(t.relative_to(ROOT)),'sha256':hashlib.sha256(t.read_bytes()).hexdigest()}for t in [HELPER,previous,proof,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_constant_modules_main():
    """The complete22 canonical modules; every scalar/resonant complement retained."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_MODULES_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole22 affineOre finite budget300s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_ORE_20261107.json';old=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==35
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;y=s.Matrix([a,b,z,w]);indices=[0,1,3,4];ww=s.symbols('W0:18')
    def makeQ(ab,az,aw,bz,bw):
        Q=s.zeros(4)
        for i,j,v in [(0,1,ab),(0,2,az),(0,3,aw),(1,2,bz),(1,3,bw)]:Q[i,j]=v;Q[j,i]=-v
        return Q
    def Wbasis(W0):
        return W0+s.Matrix([z*(ww[0]*a+ww[1]*b)+w*(ww[2]*a+ww[3]*b)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*a+ww[8]*b)+w*(ww[9]*a+ww[10]*b)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w])
    def pad(Q,W):
        P=s.zeros(5)
        for i,j in combinations(range(4),2):P[indices[i],indices[j]]=Q[i,j];P[indices[j],indices[i]]=-Q[i,j]
        for i in range(4):P[2,indices[i]]=W[i];P[indices[i],2]=-W[i]
        return P
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(t)for t in row]for row in M.tolist()],M.shape,domain);rr,columns=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in columns];kernel=[[domain.zero for j in free]for i in range(M.cols)];reduced=rr.to_list()
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(columns):kernel[pivot][j]=-reduced[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list())
        minor=dm.extract(list(rows),list(columns)).det();assert minor!=domain.zero or not columns
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(columns),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
            'minor_determinant':str(s.factor(domain.to_sympy(minor))),'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    def analyze(M,C,Dc,Dl,domain,offset):
        cm,minor,K=certify(M,domain);ca,aminor,KA=certify(M.row_join(-C),domain)
        record={'full350_affine_matrix':serial_matrix(M),'full350_affine_constant':serial_matrix(C),'linear_certificate':cm,
            'augmented_certificate':ca,'image_dimension_offset':offset};denominators=[]
        if cm['rank']<ca['rank']:
            record.update(status='inconsistent',G_closed=True,G_empty=True);return record,[minor,aminor]
        dm=DomainMatrix([[domain.from_sympy(t)for t in row]for row in M.row_join(-C).tolist()],(M.rows,19),domain);rr,pivots=dm.rref();raw=rr.to_list()
        part=s.zeros(18,1)
        for i,p in enumerate(pivots):part[p]=domain.to_sympy(raw[i][18])
        assert all(domain.from_sympy(t)==domain.zero for t in M*part+C)
        D0=s.Poly(s.expand(Dc+sum(Dl[i]*part[i]for i in range(18))),*x,domain=domain)
        directions=[s.Poly(s.expand(sum(Dl[i]*K[i,j]for i in range(18))),*x,domain=domain)for j in range(K.cols)]
        common=D0
        for poly in directions:
            if not poly.is_zero:common=s.polys.polytools.gcd(common,poly).monic()
        factor=None
        if common.is_zero:zero=True
        else:
            zero=False
            if 0<common.total_degree()<5:factor=common.monic()
            elif common.total_degree()==5:
                for poly,power in common.factor_list()[1]:
                    if 0<poly.total_degree()<5:factor=poly.monic();break
        upper=offset+K.cols
        record.update(status='whole_affine_module',particular_solution=list(map(str,part)),W_kernel_dimension=K.cols,
            full_actual_D_particular=str(D0.as_expr()),full_actual_D_kernel_coefficients=[str(poly.as_expr())for poly in directions],
            common_D_polynomial_gcd=str(common.as_expr()),proper_D_factor=None if factor is None else str(factor.as_expr()),
            whole_actual_D_zero=zero,whole_projective_image_upper=upper,G_empty=zero or factor is not None,G_closed=zero or factor is not None or upper<=24)
        if domain==FIELD:
            denominators.extend(s.denom(s.cancel(t))for t in list(part)+list(K))
            if factor is not None:
                for poly in [factor]+[poly.exquo(factor)for poly in [D0]+directions]:denominators.extend(s.denom(s.cancel(t))for t in poly.coeffs())
        return record,[minor,aminor]+denominators
    models=[];t=TAU;Wlog=s.Matrix([a*a,-2*a*b,-2*a*z,3*t*a*w]);We1=s.Matrix([a*a,-2*a*b,a*z,3*t*a*w-3*b*w]);Wtri=s.Matrix([0,-3*a*a,0,-3*a*w])
    for eps in [0,1]:
        for eta in [0,1]:
            F=eps*z*z+eta*w*w;zeros=2-eps-eta
            models.append(('EXACT_B_E1_'+str(eps)+str(eta),makeQ(-a*z+F,0,-z*w,z*z,(t+s.Rational(2,3))*z*w),We1,9,4-zeros,True,[]))
            models.append(('LOG_D_'+str(eps)+str(eta),makeQ(b*z+F,z*z,(4-3*t)*z*w/6,0,z*w),Wlog,8,5-zeros,True,[]))
            models.append(('SCALAR_B_'+str(eps)+str(eta),makeQ(b*w+eps*z*w+eta*w*w,z*w,w*w,z*z,0),Wlog.subs(t,-s.Rational(2,3)),7,5-zeros,False,[]))
            models.append(('TRIPLE_LINEAR_'+str(eps)+str(eta),makeQ(2*a*z+F,0,z*w,-2*z*z,-2*z*w),Wtri,7,5-zeros,False,[]))
    for eps in [0,1]:
        Lam=(4-3*t)/(2-6*t)
        models.append(('LOG_B_'+str(eps),makeQ(Lam*b*w+eps*z*w,Lam*z*w,w*w,z*z,0),Wlog,8,4+eps,True,[s.Rational(1,3),-s.Rational(2,3)]))
        models.append(('CRITICAL_CD_'+str(eps),makeQ(a*w+2*b*z+eps*z*z,2*z*z,3*z*w,-2*z*w,-3*w*w+z*w),Wlog.subs(t,-s.Rational(5,3)),7,5+eps,False,[]))
        models.append(('CRITICAL_A_'+str(eps),makeQ(eps*w*w,0,z*z,-z*z,-2*z*w),Wlog.subs(t,-s.Rational(8,3)),7,4+eps,False,[]))
    cases=[]
    for name,Q,W0,Bdim,Qdim,variable,excluded in models:
        strict={(i,j,k):s.cancel(sum(Q[i,l]*s.diff(Q[j,k],y[l])+Q[j,l]*s.diff(Q[k,i],y[l])+Q[k,l]*s.diff(Q[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        assert all(v==0 for v in strict.values())
        Pi=pad(Q,Wbasis(W0));full=s.Matrix([s.expand(v)for v in H.all350(H.projective(H.bilinear(Pi,Pi,x),x),x)])
        C=full.subs(dict.fromkeys(ww,0));M=full.jacobian(ww);assert all(s.cancel(v)==0 for v in full-C-M*s.Matrix(ww))
        D=divisor(Pi);Dc=D.subs(dict.fromkeys(ww,0));Dl=[s.expand(D).coeff(v)for v in ww];assert s.expand(D-Dc-sum(v*q for v,q in zip(Dl,ww)))==0
        offset=Bdim+Qdim+5;record,den=analyze(M,C,Dc,Dl,FIELD if variable else s.QQ,offset)
        record.update(type=name,whole_alpha_line=variable,full_Q={str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},W0=list(map(str,W0)),
            full_strict_Q_Jacobi_identically_zero=True,whole_fixed_plane_symbol_dimension_upper=Bdim,whole_Q_source_dimension_upper=Qdim)
        if variable:
            exception=s.Poly(1,t)
            for value in den+list(M)+list(C)+Dl+[Dc]:
                if value in den:num,de=s.fraction(s.cancel(value));exception=s.lcm(exception,s.Poly(num,t));exception=s.lcm(exception,s.Poly(de,t))
                else:exception=s.lcm(exception,s.Poly(s.denom(s.cancel(value)),t))
            for value in excluded:exception=s.lcm(exception,s.Poly(t-value,t))
            walls=[]
            for factor,power in s.factor_list(exception.as_expr())[1]:
                p=s.Poly(factor,t);alpha=-p.nth(0)/p.nth(1)if p.degree()==1 else s.CRootOf(p.as_expr(),0)
                if alpha in excluded:
                    walls.append({'alpha':str(alpha),'source_excluded_here':True,'reason':'leadingzero'if alpha==s.Rational(1,3)else'coveredbyallfourSCALAR_Bfixedmodules','G_closed':True});continue
                dom=s.QQ if p.degree()==1 else s.QQ.algebraic_field(alpha);sub=lambda value:value.subs(t,alpha)
                wall,_=analyze(sub(M),sub(C),sub(Dc),[sub(v)for v in Dl],dom,offset-1)
                wall.update(alpha=str(alpha),minimal_polynomial=str(factor),all_complex_embeddings_covered=True);walls.append(wall)
            record.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_irreducible_factors=[[str(v),int(e)]for v,e in s.factor_list(exception.as_expr())[1]],
                all_exceptional_fields=walls,all_strata_G_closed=record['G_closed']and all(v['G_closed']for v in walls))
        else:record['all_strata_G_closed']=record['G_closed']
        cases.append(record);print(name,record['status'],'k',record['linear_certificate']['nullity'],'factor',record.get('proper_D_factor'),
            'upper',record.get('whole_projective_image_upper'),'allG',record['all_strata_G_closed'],
            'walls',[(v['alpha'],v.get('W_kernel_dimension'),v.get('proper_D_factor'),v['G_closed'])for v in record.get('all_exceptional_fields',[])],flush=True)
    proofs=[ROOT/'research/classification'/name for name in ['PLANE_NONREDUCED_M3_CONSTANT_ROOTS_ORE_GEOMETRY.md','PLANE_NONREDUCED_M3_LOG_ORE_NORMALFORMS.md']]
    archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_constant_ore_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all22_G_closed':all(v['all_strata_G_closed']for v in cases),'cases':cases,
        'scope':'Entire22canonical modulescover allnonconstantfactor branches of6constantrootm3 residualsources, withallresonances/zeroresidualcentralclasses. Allfull350affineW18 andallactualD affinecoefficients, entirealphalines/allexceptionfields. Allsupportnormalizerproofs areexplicitseparateinputs.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,archive]+proofs],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_moving_double_main():
    """Whole11 first connections with moving double row degree1/simple row degree0."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole11 movingdouble finite budget300s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_MODULES_20261107.json';old=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==36
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;base=[a,b,c];normals=[z,w];P1,v1,P2,v2=sources(x);kk=s.symbols('K0:3')
    monG=[v+(j,3-j)for v in H.monomials(3,2)for j in range(4)];monK=[v+(j,2-j)for v in H.monomials(3,3)for j in range(3)]
    def dlog(g):return s.Matrix([s.cancel(s.diff(g,t)/g)for t in base])
    def differential(g):return s.Matrix([s.diff(g,t)for t in base])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list();kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());minor=dm.extract(list(rows),list(cols)).det()
        assert minor!=domain.zero or not cols;K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(minor))),
            'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    models=[];f=a*a*b
    for name,C,A,Bdim in [('REPEATED_LINE_INTERSECTION',a,b,9),('REPEATED_LINE_OTHER',a,c,10)]:
        p=dlog(f)/3+TAU*dlog(a/b);models.append((name,f,A,C,p,Bdim,True))
    for geom,C,A,offset in [('SIMPLE_LINE_INTERSECTION',b,a,8),('SIMPLE_LINE_OTHER',b,c,9)]:
        for name,gamma,delta,dim in [('EXACT_C',0,1,2),('EXACT_B',1,0,1),('ZERO',0,0,0)]:
            p=dlog(f)/3+dlog(a/b)/3+gamma*differential(b/a)+delta*differential(c/a)
            models.append((geom+'_'+name,f,A,C,p,offset+dim,False))
    f=a**3
    for name,gamma,delta,dim in [('EXACT_C',0,1,2),('EXACT_B',1,0,1),('ZERO',0,0,0)]:
        models.append(('TRIPLE_LINE_'+name,f,b,a,dlog(a)+gamma*differential(b/a)+delta*differential(c/a),7+dim,False))
    records=[]
    for name,f,A,C,p,Bdim,variable in models:
        q=dlog(f/C**2);r=(-differential(A/C)+(A/C)*(q-p)).applyfunc(s.cancel);beta=[[p,s.zeros(3,1)],[r,q]]
        Pbase=s.zeros(5)
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:Pbase[i,j]=s.diff(f,base[k]);Pbase[j,i]=-Pbase[i,j]
        omegas=[]
        for i in range(2):
            for j in range(2):
                om=(3*f*beta[i][j]-(differential(f)if i==j else s.zeros(3,1))).applyfunc(s.cancel)
                assert s.expand(sum(base[k]*om[k]for k in range(3)))==0
                assert all(s.denom(v)==1 for v in om)
                V=s.Matrix([s.diff(om[2],b)-s.diff(om[1],c),s.diff(om[0],c)-s.diff(om[2],a),s.diff(om[1],a)-s.diff(om[0],b)])/3
                assert s.expand(sum(s.diff(V[k],base[k])for k in range(3)))==0
                inverse=s.Matrix([c*V[1]-b*V[2],a*V[2]-c*V[0],b*V[0]-a*V[1]])
                assert all(s.expand(v)==0 for v in inverse-om)
                for k in range(3):Pbase[k,3+i]-=V[k]*normals[j];Pbase[3+i,k]+=V[k]*normals[j]
                omegas.append({'i':i,'j':j,'beta':list(map(str,beta[i][j])),'horizontal_polynomial_form':list(map(str,om)),'inverse_tracefree_field':list(map(str,V))})
        target=s.expand(z*(A*z+C*w)**2);target_vector=s.Matrix([s.Poly(target,*x).coeff_monomial(mon)for mon in monG]);idx=next(i for i,v in enumerate(target_vector)if v!=0)
        Ann=s.zeros(23,24)
        for row,i in enumerate(j for j in range(24)if j!=idx):Ann[row,i]=1;Ann[row,idx]=-target_vector[i]/target_vector[idx]
        assert Ann*target_vector==s.zeros(23,1)
        def analyze(alpha=None,domain=FIELD):
            sub=lambda v:v if alpha is None else v.subs(TAU,alpha)
            Pk=sub(Pbase);Pk[3,4]=kk[0]*z*z+kk[1]*z*w+kk[2]*w*w;Pk[4,3]=-Pk[3,4]
            Dk=divisor(Pk);poly=s.Poly(Dk,*x);rows=s.Matrix([poly.coeff_monomial(mon)for mon in monK]);KM=rows.jacobian(kk);const=rows.subs(dict.fromkeys(kk,0))
            kc,km,_=certify(KM,domain);ka,am,_=certify(KM.row_join(-const),domain);den=[km,am]
            record={'full30_global_K_matrix':serial_matrix(KM),'full30_global_K_constant':serial_matrix(const),'K_certificate':kc,'K_augmented_certificate':ka}
            if ka['rank']>kc['rank']:record.update(global_K_inconsistent=True,whole_m3_source_empty=True,G_closed=True);return record,den
            sol=list(s.linsolve((KM,-const),kk));assert len(sol)==1 and kc['rank']==3
            P0=sub(Pbase);P0[3,4]=s.expand(sum(v*m for v,m in zip(sol[0],[z*z,z*w,w*w])));P0[4,3]=-P0[3,4]
            assert all(domain.from_sympy(v)==domain.zero for v in s.Poly(divisor(P0),*x).coeffs())
            Jac0=[s.expand(v)for v in H.all350(H.projective(H.bilinear(P0,P0,x),x),x)];bad=[(i,v)for i,v in enumerate(Jac0)if domain.from_sympy(v)!=domain.zero]
            record.update(unique_global_K=str(P0[3,4]),full150_P0={str((i,j)):str(P0[i,j])for i,j in H.PAIRS})
            if bad:
                record.update(P0_nonzero_Jacobi=[(i,str(v))for i,v in bad],whole_m3_source_empty=True,G_closed=True);den.append(bad[0][1]);return record,den
            D1=operator(P0,P1,v1,x);D2=operator(P0,P2,v2,x);full=divisor(P0+P1+P2);pp=s.Poly(full,*x)
            D3=s.Matrix([pp.coeff_monomial(mon)for mon in monG]);L=D3.jacobian(v1);assert all(s.expand(v)==0 for v in D3-L*s.Matrix(v1))
            M1=D1.col_join(Ann*L);Sh=shear(P0,x);assert all(s.expand(v)==0 for v in M1*Sh)
            c1,m1,K1=certify(M1,domain);c2,m2,K2=certify(D2,domain);cs,ms,KS=certify(Sh,domain);image=(L*K1).applyfunc(s.cancel)
            assert all(domain.from_sympy(v)==domain.zero for v in Ann*image)
            leading=[domain.from_sympy(image[idx,j]/target_vector[idx])for j in range(image.cols)];nonzero=[j for j,v in enumerate(leading)if v!=domain.zero]
            offset=Bdim+5-(1 if alpha is not None and variable else 0);upper=offset+K1.cols+K2.cols
            record.update(P0_full350_Jacobi_zero_and_actualD_zero=True,full350_weight1_operator=serial_matrix(D1),full350_weight2_operator=serial_matrix(D2),
                full24by30_actual_D3_map=serial_matrix(L),full_weight1_plus_line_matrix=serial_matrix(M1),complete_actual_shear=serial_matrix(Sh),
                weight1_certificate=c1,weight2_certificate=c2,actual_shear_certificate=cs,whole_D3_kernel_image=serial_matrix(image),
                nonzero_D3_line_certificate=None if not nonzero else {'kernel_column':nonzero[0],'coefficient':str(domain.to_sympy(leading[nonzero[0]]))},
                D3_image_rank=int(bool(nonzero)),whole_m3_source_empty=not nonzero,whole_projective_image_upper=upper,G_closed=not nonzero or upper<=24)
            den.extend([m1,m2,ms]);den.extend(s.denom(s.cancel(v))for v in list(K1)+list(K2)+list(KS)+list(sol[0]))
            return record,den
        result,den=analyze(domain=FIELD if variable else s.QQ);result.update(type=name,whole_alpha_line=variable,full_connection_inverse_forms=omegas,
            whole_fixed_plane_symbol_dimension_upper=Bdim,exact_nonzero_leading_target=str(target),full_line_annihilator=serial_matrix(Ann))
        if variable:
            exception=s.Poly(1,TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            walls=[]
            for factor,power in s.factor_list(exception.as_expr())[1]:
                poly=s.Poly(factor,TAU);alpha=-poly.nth(0)/poly.nth(1)if poly.degree()==1 else s.CRootOf(poly.as_expr(),0)
                domain=s.QQ if poly.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(factor),all_complex_embeddings_covered=True);walls.append(wall)
            result.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,
                all_strata_G_closed=result['G_closed']and all(v['G_closed']for v in walls))
        else:result['all_strata_G_closed']=result['G_closed']
        records.append(result);print(name,'empty',result['whole_m3_source_empty'],'k',result.get('weight1_certificate',{}).get('nullity'),
            result.get('weight2_certificate',{}).get('nullity'),'upper',result.get('whole_projective_image_upper'),
            'walls',[(v['alpha'],v['whole_m3_source_empty'],v.get('whole_projective_image_upper'))for v in result.get('all_exceptional_sources',[])],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_MOVING_DOUBLE_FIRST_CONNECTIONS.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_constant_modules_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all11_G_closed':all(v['all_strata_G_closed']for v in records),'cases':records,
        'scope':'Whole nonredR m3generic21 primitive degrees(1,0): allglobalclosedfirstconnections/twoalphalines+ninefixed, inverseforms/globalK/P0Jac/full30/9/actualD3exactline,allalphawalls. Constantroot/otherdegrees/generictriple/m2/reducedR outside.',
        'all24_D3_monomials':monG,'all30_K_monomials':monK,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_moving_double_ore_main():
    """Whole two residual sources, exactQ/W bridge and three complete modules."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole2source3module Ore finite budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_20261107.json';packet=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(packet['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==37
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;y=s.Matrix([a,b,z,w]);indices=[0,1,3,4];P1,v1,P2,v2=sources(x)
    ww=s.symbols('W0:18');qf=s.symbols('F0:3');u=s.Symbol('u');F=qf[0]*z*z+qf[1]*z*w+qf[2]*w*w;A=(2-6*TAU)/3;B=A+1
    def makeQ(ab,az,aw,bz,bw):
        Q=s.zeros(4)
        for i,j,v in [(0,1,ab),(0,2,az),(0,3,aw),(1,2,bz),(1,3,bw)]:Q[i,j]=v;Q[j,i]=-v
        return Q
    def Qdouble(u,F):return makeQ(u*b*z+F,u*A*z*z,u*B*z*w,-u*A*z*w,-u*B*w*w)
    def Qtriple(u,F):return makeQ(u*a*z+F,2*u*z*z,u*(-z*z+2*z*w),-2*u*z*w,u*(z*w-2*w*w))
    Wd=s.Matrix([a*a,-2*a*b,3*TAU*a*z,-2*a*w+(1-3*TAU)*b*z]);Wt=s.Matrix([0,-3*a*a,-3*a*z,3*(a+b)*z])
    def Wbasis(W0):
        return W0+s.Matrix([z*(ww[0]*a+ww[1]*b)+w*(ww[2]*a+ww[3]*b)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*a+ww[8]*b)+w*(ww[9]*a+ww[10]*b)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w])
    def pad(Q,W):
        P=s.zeros(5)
        for i,j in combinations(range(4),2):P[indices[i],indices[j]]=Q[i,j];P[indices[j],indices[i]]=-Q[i,j]
        for i in range(4):P[2,indices[i]]=W[i];P[indices[i],2]=-W[i]
        return P
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    coverage=[]
    for name,Qfun,W0,variable in [('REPEATED_LINE_INTERSECTION',Qdouble,Wd,True),('TRIPLE_LINE_EXACT_B',Qtriple,Wt,False)]:
        old=next(v for v in packet['cases']if v['type']==name);K=s.Matrix([[parse(v)for v in row]for row in old['weight1_certificate']['complete_kernel']]);tt=s.symbols('p0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=parse(old['full150_P0'][str((i,j))]);P0[j,i]=-P0[i,j]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(tt))),simultaneous=True)+P2
        gauge=s.Matrix([-s.diff(P[0,2],c),-s.diff(P[1,2],c),0,0,0]);gauged=(P-H.wedge(x,gauge)).applyfunc(s.expand)
        assert all(s.diff(gauged[i,j],c)==0 for i,j in H.PAIRS)
        Q=gauged.extract(indices,indices);W=s.Matrix([gauged[2,i]for i in indices]);ell=-s.Poly(W[2],z,w).coeff_monomial(z*z)*z-s.Poly(W[3],z,w).coeff_monomial(w*w)*w
        W=(W+ell*y).applyfunc(s.expand);noncommon=s.Poly(Q[0,1],a,b,z,w).coeff_monomial(b*z if variable else a*z);central=s.expand(Q[0,1]-noncommon*(b*z if variable else a*z))
        expected=Qfun(noncommon,central);assert all(s.expand(v)==0 for v in Q-expected)
        wd=W-W0;wc=[]
        for i in range(2):wc.extend([s.Poly(wd[i],*y).coeff_monomial(mon)for mon in[a*z,b*z,a*w,b*w,z*z,z*w,w*w]])
        wc.extend([s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)])
        assert all(s.expand(v)==0 for v in W-Wbasis(W0).subs(dict(zip(ww,wc)),simultaneous=True))
        coords=[noncommon]+[s.Poly(central,z,w).coeff_monomial(mon)for mon in[z*z,z*w,w*w]]+wc;J=s.Matrix(coords).jacobian(list(tt)+v2);det=s.factor(J.det(method='domain-ge'))
        assert J.rows==J.cols==22 and det!=0 and not s.sympify(det).free_symbols
        sourceQ=Qfun(u,F);strict={(i,j,k):s.expand(sum(sourceQ[i,l]*s.diff(sourceQ[j,k],y[l])+sourceQ[j,l]*s.diff(sourceQ[k,i],y[l])+sourceQ[k,l]*s.diff(sourceQ[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        assert all(v==0 for v in strict.values());D=divisor(pad(sourceQ,Wbasis(W0)));poly=s.Poly(D,*x);D3=sum(value*s.prod(v**e for v,e in zip(x,mon))for mon,value in poly.terms()if mon[3]+mon[4]==3)
        if variable:
            expectedD3=s.Rational(3,2)*u*A*(A-1)*z*(b*z+a*w)**2;assert s.expand(D3-expectedD3)==0
            for wall in old['all_exceptional_sources']:
                alpha=parse(wall['alpha'])
                if wall['whole_m3_source_empty']:assert alpha in[-s.Rational(1,6),s.Rational(1,3)]
                else:assert wall['weight1_certificate']['nullity']==13
        else:assert s.expand(D3-6*u*z*(b*z+a*w)**2)==0
        coverage.append({'type':name,'full_source_first_Euler_gauge':list(map(str,gauge)),'second_gauge_c_component':str(ell),
            'full_c_independence_and_QW18_identity':True,'joint22by22_source_coordinate_matrix':serial_matrix(J),'joint_nonzero_constant_determinant':str(det),
            'all_strict_Q_Jacobi_zero':True,'full_D3':str(s.factor(D3)),'all_other_old_walls_have_same_complete13_source':True})
    aa,bb,cc,dd=s.symbols('sa sb sc sd');normalshift=s.Matrix([aa*z+bb*w,cc*z+dd*w,0,0]);S=y+normalshift;Sinv=y-normalshift;Jac=S.jacobian(y)
    shiftdata=[]
    for name,Qfun in [('DOUBLE',Qdouble),('TRIPLE',Qtriple)]:
        Q=Qfun(u,F);transport=(Jac*Q.subs(dict(zip(y,Sinv)),simultaneous=True)*Jac.T).applyfunc(s.expand)
        delta=u*((A-1)*cc*z*z+A*(aa+dd)*z*w+B*bb*w*w)if name=='DOUBLE'else u*((2*cc-dd-aa)*z*z+(2*dd+2*aa-2*bb)*z*w+2*bb*w*w)
        assert all(s.expand(v)==0 for v in transport-Qfun(u,F+delta));shiftdata.append({'type':name,'full_actual_central_shear_verified':True,'central_change':str(delta)})
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list();kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());minor=dm.extract(list(rows),list(cols)).det();assert minor!=domain.zero or not cols
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(minor))),'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    def analyze(M,C,Dc,Dl,domain,offset):
        cm,mi,K=certify(M,domain);ca,am,KA=certify(M.row_join(-C),domain);rec={'full350_matrix':serial_matrix(M),'full350_constant':serial_matrix(C),'linear_certificate':cm,'augmented_certificate':ca,'image_offset':offset};den=[mi,am]
        if ca['rank']>cm['rank']:rec.update(status='inconsistent',G_closed=True,G_empty=True);return rec,den
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.row_join(-C).tolist()],(350,19),domain);rr,pivots=dm.rref();raw=rr.to_list();part=s.zeros(18,1)
        for i,j in enumerate(pivots):part[j]=domain.to_sympy(raw[i][18])
        assert all(domain.from_sympy(v)==domain.zero for v in M*part+C)
        D0=s.Poly(s.expand(Dc+sum(Dl[i]*part[i]for i in range(18))),*x,domain=domain);directions=[s.Poly(s.expand(sum(Dl[i]*K[i,j]for i in range(18))),*x,domain=domain)for j in range(K.cols)]
        common=D0
        for poly in directions:
            if not poly.is_zero:common=s.polys.polytools.gcd(common,poly).monic()
        factor=None;zero=common.is_zero
        if not zero and 0<common.total_degree()<5:factor=common.monic()
        elif not zero and common.total_degree()==5:
            factor=next((poly.monic()for poly,power in common.factor_list()[1]if 0<poly.total_degree()<5),None)
        rec.update(status='complete_affine_module',particular_solution=list(map(str,part)),full_D_particular=str(D0.as_expr()),
            full_D_kernel_coefficients=[str(poly.as_expr())for poly in directions],common_D_gcd=str(common.as_expr()),proper_D_factor=None if factor is None else str(factor.as_expr()),
            whole_D_zero=zero,W_kernel_dimension=K.cols,whole_projective_image_upper=offset+K.cols,G_empty=zero or factor is not None,G_closed=zero or factor is not None or offset+K.cols<=24)
        if domain==FIELD:
            den.extend(s.denom(s.cancel(v))for v in list(K)+list(part))
            if factor is not None:
                for poly in[factor]+[poly.exquo(factor)for poly in[D0]+directions]:den.extend(s.denom(s.cancel(v))for v in poly.coeffs())
        return rec,den
    cases=[]
    for name,Q,W0,variable,offset in [('DOUBLE_WHOLE',Qdouble(1,0),Wd,True,18),('DOUBLE_CRITICAL_F',Qdouble(1,w*w).subs(TAU,s.Rational(5,6)),Wd.subs(TAU,s.Rational(5,6)),False,17),('TRIPLE_FIXED',Qtriple(1,0),Wt,False,17)]:
        P=pad(Q,Wbasis(W0));rows=s.Matrix([s.expand(v)for v in H.all350(H.projective(H.bilinear(P,P,x),x),x)]);C=rows.subs(dict.fromkeys(ww,0));M=rows.jacobian(ww)
        assert all(s.expand(v)==0 for v in rows-C-M*s.Matrix(ww));D=divisor(P);Dc=D.subs(dict.fromkeys(ww,0));Dl=[s.expand(D).coeff(v)for v in ww];assert s.expand(D-Dc-sum(v*q for v,q in zip(Dl,ww)))==0
        rec,den=analyze(M,C,Dc,Dl,FIELD if variable else s.QQ,offset);rec.update(type=name,full_Q={str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},W0=list(map(str,W0)))
        if variable:
            exception=s.Poly((TAU+s.Rational(1,6))*(TAU-s.Rational(1,3)),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            walls=[]
            for pol,power in s.factor_list(exception.as_expr())[1]:
                p=s.Poly(pol,TAU);alpha=-p.nth(0)/p.nth(1)if p.degree()==1 else s.CRootOf(p.as_expr(),0)
                if alpha in[-s.Rational(1,6),s.Rational(1,3)]:walls.append({'alpha':str(alpha),'source_empty_by_original_full_D3_image':True,'G_closed':True});continue
                domain=s.QQ if p.degree()==1 else s.QQ.algebraic_field(alpha);sub=lambda v:v.subs(TAU,alpha);wall,_=analyze(sub(M),sub(C),sub(Dc),[sub(v)for v in Dl],domain,offset-1)
                wall.update(alpha=str(alpha),minimal_polynomial=str(pol),all_complex_embeddings_covered=True);walls.append(wall)
            rec.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,all_strata_G_closed=rec['G_closed']and all(v['G_closed']for v in walls))
        else:rec['all_strata_G_closed']=rec['G_closed']
        cases.append(rec);print(name,rec['status'],'k',rec['linear_certificate']['nullity'],'factor',rec.get('proper_D_factor'),'G',rec['all_strata_G_closed'],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_REDUCTION.md';erratum=ROOT/'research/classification/PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_LEADING_ERRATUM.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_moving_double_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_three_modules_G_closed':all(v['all_strata_G_closed']for v in cases),'full_two_source_coverage':coverage,'actual_full_central_shear_identities':shiftdata,'cases':cases,
        'scope':'Whole2movingdouble residualsources; full22source Euler/Q/W18bijection, allcohomcentralinverse/leadingcriteria, onewholeline+criticalF+fixedmodule withfull350/actualD andallparametercomplements. Otherprimitiveflags/generictriple/m2/reducedRoutside.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,erratum,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_moving_simple_main():
    """Whole37 first connections forconstantdouble/movingsimple degrees1and2."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole37 movingsimple finite budget300s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_20261107.json';old=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==38
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;base=[a,b,c];normals=[z,w];P1,v1,P2,v2=sources(x);kk=s.symbols('K0:3')
    monG=[v+(j,3-j)for v in H.monomials(3,2)for j in range(4)];monK=[v+(j,2-j)for v in H.monomials(3,3)for j in range(3)]
    def dlog(g):return s.Matrix([s.cancel(s.diff(g,t)/g)for t in base])
    def differential(g):return s.Matrix([s.diff(g,t)for t in base])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list();kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());minor=dm.extract(list(rows),list(cols)).det()
        assert minor!=domain.zero or not cols;K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(minor))),
            'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    models=[];f=a*a*b
    for h in[a,b]:
        hs='A'if h==a else'B'
        for name,C,A,Bdim in [('INTERSECTION',a,b,9),('OTHER',a,c,10)]:
            q=dlog(f)/3+TAU*dlog(a/b);models.append(('DEG1_H'+hs+'_CA_'+name,f,A,C,h,q,Bdim,True))
        for geom,A,offset in [('INTERSECTION',a,8),('OTHER',c,9)]:
            for nm,gamma,dim in [('EXACT_B',1,1),('ZERO',0,0)]:
                alpha=s.Rational(4,3)if h==a else s.Rational(7,3)
                q=dlog(f)/3+alpha*dlog(a/b)+gamma*differential(b/a)
                models.append(('DEG1_H'+hs+'_CB_'+geom+'_'+nm,f,A,b,h,q,offset+dim,False))
    f=a**3
    for nm,gamma,delta,dim in [('EXACT_C',0,1,2),('EXACT_B',1,0,1),('ZERO',0,0,0)]:
        q=dlog(a)+gamma*differential(b/a)+delta*differential(c/a)
        models.append(('DEG1_TRIPLE_'+nm,f,b,a,a,q,7+dim,False))
    f=a*a*b
    Asquare=[c*c+eps*a*b+eta*b*b for eps in[0,1]for eta in[0,1]]+[b*c,c*(a+b),a*c+b*b,b*b,a*b+b*b]
    Across=[c*c+eps*a*a+eta*b*b for eps in[0,1]for eta in[0,1]]+[c*(a+b),a*c+b*b,b*c+a*a,a*a+b*b]
    for j,Arow in enumerate(Asquare):models.append(('DEG2_CA2_'+str(j),f,Arow,a*a,s.Integer(1),-2*dlog(a)+3*dlog(b),12,False))
    for j,Arow in enumerate(Across):models.append(('DEG2_CAB_'+str(j),f,Arow,a*b,s.Integer(1),2*dlog(a)-dlog(b),12,False))
    for gamma in[0,1]:
        q=4*dlog(a)-3*dlog(b)+gamma*differential(b/a)
        models.append(('DEG2_CB2_'+str(gamma),f,a*a+gamma*a*b,b*b,s.Integer(1),q,8+gamma,False))
    f=a**3
    for j,Arow,Bdim in[(0,b*c,10),(1,b*b,8),(2,b*b+a*c,9)]:
        models.append(('DEG2_TRIPLE_'+str(j),f,Arow,a*a,s.Integer(1),dlog(a),Bdim,False))
    assert len(models)==37
    records=[]
    for name,f,A,C,h,q,Bdim,variable in models:
        p=dlog(f/(h*C));r=(-differential(A/C)+(A/C)*(q-p)).applyfunc(s.cancel);beta=[[p,s.zeros(3,1)],[r,q]]
        Pbase=s.zeros(5)
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:Pbase[i,j]=s.diff(f,base[k]);Pbase[j,i]=-Pbase[i,j]
        omegas=[]
        for i in range(2):
            for j in range(2):
                om=(3*f*beta[i][j]-(differential(f)if i==j else s.zeros(3,1))).applyfunc(s.cancel)
                assert s.expand(sum(base[k]*om[k]for k in range(3)))==0
                assert all(s.denom(v)==1 for v in om)
                V=s.Matrix([s.diff(om[2],b)-s.diff(om[1],c),s.diff(om[0],c)-s.diff(om[2],a),s.diff(om[1],a)-s.diff(om[0],b)])/3
                assert s.expand(sum(s.diff(V[k],base[k])for k in range(3)))==0
                inverse=s.Matrix([c*V[1]-b*V[2],a*V[2]-c*V[0],b*V[0]-a*V[1]])
                assert all(s.expand(v)==0 for v in inverse-om)
                for k in range(3):Pbase[k,3+i]-=V[k]*normals[j];Pbase[3+i,k]+=V[k]*normals[j]
                omegas.append({'i':i,'j':j,'beta':list(map(str,beta[i][j])),'horizontal_polynomial_form':list(map(str,om)),'inverse_tracefree_field':list(map(str,V))})
        target=s.expand(h*z*z*(A*z+C*w));target_vector=s.Matrix([s.Poly(target,*x).coeff_monomial(mon)for mon in monG]);idx=next(i for i,v in enumerate(target_vector)if v!=0)
        Ann=s.zeros(23,24)
        for row,i in enumerate(j for j in range(24)if j!=idx):Ann[row,i]=1;Ann[row,idx]=-target_vector[i]/target_vector[idx]
        assert Ann*target_vector==s.zeros(23,1)
        def analyze(alpha=None,domain=FIELD):
            sub=lambda v:v if alpha is None else v.subs(TAU,alpha)
            Pk=sub(Pbase);Pk[3,4]=kk[0]*z*z+kk[1]*z*w+kk[2]*w*w;Pk[4,3]=-Pk[3,4]
            Dk=divisor(Pk);poly=s.Poly(Dk,*x);rows=s.Matrix([poly.coeff_monomial(mon)for mon in monK]);KM=rows.jacobian(kk);const=rows.subs(dict.fromkeys(kk,0))
            kc,km,_=certify(KM,domain);ka,am,_=certify(KM.row_join(-const),domain);den=[km,am]
            record={'full30_global_K_matrix':serial_matrix(KM),'full30_global_K_constant':serial_matrix(const),'K_certificate':kc,'K_augmented_certificate':ka}
            if ka['rank']>kc['rank']:record.update(global_K_inconsistent=True,whole_m3_source_empty=True,G_closed=True);return record,den
            sol=list(s.linsolve((KM,-const),kk));assert len(sol)==1 and kc['rank']==3
            P0=sub(Pbase);P0[3,4]=s.expand(sum(v*m for v,m in zip(sol[0],[z*z,z*w,w*w])));P0[4,3]=-P0[3,4]
            assert all(domain.from_sympy(v)==domain.zero for v in s.Poly(divisor(P0),*x).coeffs())
            Jac0=[s.expand(v)for v in H.all350(H.projective(H.bilinear(P0,P0,x),x),x)];bad=[(i,v)for i,v in enumerate(Jac0)if domain.from_sympy(v)!=domain.zero]
            record.update(unique_global_K=str(P0[3,4]),full150_P0={str((i,j)):str(P0[i,j])for i,j in H.PAIRS})
            if bad:
                record.update(P0_nonzero_Jacobi=[(i,str(v))for i,v in bad],whole_m3_source_empty=True,G_closed=True);den.append(bad[0][1]);return record,den
            D1=operator(P0,P1,v1,x);D2=operator(P0,P2,v2,x);full=divisor(P0+P1+P2);pp=s.Poly(full,*x)
            D3=s.Matrix([pp.coeff_monomial(mon)for mon in monG]);L=D3.jacobian(v1);assert all(s.expand(v)==0 for v in D3-L*s.Matrix(v1))
            M1=D1.col_join(Ann*L);Sh=shear(P0,x);assert all(s.expand(v)==0 for v in M1*Sh)
            c1,m1,K1=certify(M1,domain);c2,m2,K2=certify(D2,domain);cs,ms,KS=certify(Sh,domain);image=(L*K1).applyfunc(s.cancel)
            qbad=s.cancel(3*f*sum(sub(q[i])*v1[-3+i]for i in range(3)));qbpoly=s.Poly(qbad,*base)
            QM=s.Matrix([qbpoly.coeff_monomial(mon)for mon in H.monomials(3,2)]).jacobian(v1)
            assert all(domain.from_sympy(v)==domain.zero for v in QM*K1)
            assert all(domain.from_sympy(v)==domain.zero for v in Ann*image)
            leading=[domain.from_sympy(image[idx,j]/target_vector[idx])for j in range(image.cols)];nonzero=[j for j,v in enumerate(leading)if v!=domain.zero]
            offset=Bdim+5-(1 if alpha is not None and variable else 0);upper=offset+K1.cols+K2.cols
            record.update(P0_full350_Jacobi_zero_and_actualD_zero=True,full350_weight1_operator=serial_matrix(D1),full350_weight2_operator=serial_matrix(D2),
                full24by30_actual_D3_map=serial_matrix(L),full_weight1_plus_line_matrix=serial_matrix(M1),complete_actual_shear=serial_matrix(Sh),
                weight1_certificate=c1,weight2_certificate=c2,actual_shear_certificate=cs,whole_D3_kernel_image=serial_matrix(image),
                full6by30_constant_root_character_matrix=serial_matrix(QM),all_full_kernel_q_Cbad_zero=True,
                nonzero_D3_line_certificate=None if not nonzero else {'kernel_column':nonzero[0],'coefficient':str(domain.to_sympy(leading[nonzero[0]]))},
                D3_image_rank=int(bool(nonzero)),whole_m3_source_empty=not nonzero,whole_projective_image_upper=upper,G_closed=not nonzero or upper<=24)
            den.extend([m1,m2,ms]);den.extend(s.denom(s.cancel(v))for v in list(K1)+list(K2)+list(KS)+list(sol[0]))
            return record,den
        result,den=analyze(domain=FIELD if variable else s.QQ);result.update(type=name,whole_alpha_line=variable,full_connection_inverse_forms=omegas,
            whole_fixed_plane_symbol_dimension_upper=Bdim,exact_nonzero_leading_target=str(target),full_line_annihilator=serial_matrix(Ann))
        if variable:
            exception=s.Poly(1,TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            walls=[]
            for factor,power in s.factor_list(exception.as_expr())[1]:
                poly=s.Poly(factor,TAU);alpha=-poly.nth(0)/poly.nth(1)if poly.degree()==1 else s.CRootOf(poly.as_expr(),0)
                domain=s.QQ if poly.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(factor),all_complex_embeddings_covered=True);walls.append(wall)
            result.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,
                all_strata_G_closed=result['G_closed']and all(v['G_closed']for v in walls))
        else:result['all_strata_G_closed']=result['G_closed']
        records.append(result);print(name,'empty',result['whole_m3_source_empty'],'k',result.get('weight1_certificate',{}).get('nullity'),
            result.get('weight2_certificate',{}).get('nullity'),'upper',result.get('whole_projective_image_upper'),
            'walls',[(v['alpha'],v['whole_m3_source_empty'],v.get('whole_projective_image_upper'))for v in result.get('all_exceptional_sources',[])],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_MOVING_SIMPLE_FIRST_CONNECTIONS.md';obstruction=ROOT/'research/classification/PLANE_NONREDUCED_M3_CONSTANT_ROOT_OBSTRUCTION_APPLICATIONS.md';rootlemma=ROOT/'research/classification/PLANE_M3_CONSTANT_ROOT_OBSTRUCTION.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_moving_double_ore_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all37_G_closed':all(v['all_strata_G_closed']for v in records),'cases':records,
        'scope':'Whole nonredR m3generic21 primitive degrees(0,1)/(0,2): fourwholealpha+33fixed, allglobalclosedq/nonlog/quadraticroworbits, fullK/P0/full30/9/actualD3exactline/allqCbadkernels. Fourdelta!=0degree1sourcesalreadyemptybyuniversalcontraction. Generictriple/m2/reducedRoutside.',
        'all24_D3_monomials':monG,'all30_K_monomials':monK,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,obstruction,rootlemma,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')



def nonreduced_m3_moving_simple_ore_main():
    """All seven residual first-symbol strata, their full vertex Euler sources."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_ORE_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole7vertex sources finite budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_20261107.json';packet=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(packet['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==39
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;P1,v1,P2,v2=sources(x);ww=s.symbols('W0:18');mons=H.monomials(4,2)
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    def matrix(raw):return s.Matrix([[parse(v)for v in row]for row in raw])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1];free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list()
        kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());mi=dm.extract(list(rows),list(cols)).det();assert mi!=domain.zero or not cols
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(mi))),'complete_kernel':serial_matrix(K)},domain.to_sympy(mi),K
    selected=[]
    for old in packet['cases']:
        if old['all_strata_G_closed']:continue
        if old['type'].endswith('_OTHER'):
            bad=[v for v in old['all_exceptional_sources']if not v['G_closed']];assert len(bad)==1
            fixed=dict(bad[0]);fixed.update(type=old['type']+'_WALL',whole_alpha_line=False,exact_nonzero_leading_target=old['exact_nonzero_leading_target'],
                full_line_annihilator=old['full_line_annihilator'],whole_fixed_plane_symbol_dimension_upper=old['whole_fixed_plane_symbol_dimension_upper']-1)
            selected.append(fixed)
        else:selected.append(old)
    assert len(selected)==7
    def analyze(old,alpha=None,domain=FIELD):
        sub=lambda v:v if alpha is None else v.subs(TAU,alpha);vertex=1 if old['type']=='DEG2_CA2_0'else 2;indices=[i for i in range(5)if i!=vertex];y=s.Matrix([x[i]for i in indices]);base0,base1=y[:2]
        L0=s.zeros(5,1)
        if old['type']=='DEG1_HA_CA_OTHER_WALL':L0[4]=-3*z
        elif old['type']=='DEG1_HB_CA_OTHER_WALL':L0[4]=3*z
        elif vertex==1:L0[2]=2*a
        M1=sub(matrix(old['full_weight1_plus_line_matrix']));c1,_,K=certify(M1,domain);tt=s.symbols('p0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=sub(parse(old['full150_P0'][str((i,j))]));P0[j,i]=-P0[i,j]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(tt))),simultaneous=True)+P2-H.wedge(x,L0)
        L=s.Matrix([-s.diff(P[i,vertex],x[vertex])if i!=vertex else 0 for i in range(5)])
        gauged=(P-H.wedge(x,L)).applyfunc(s.expand);bad=[(str((i,j)),str(s.diff(gauged[i,j],x[vertex])))for i,j in H.PAIRS if s.diff(gauged[i,j],x[vertex])!=0]
        rec={'alpha':None if alpha is None else str(alpha),'vertex_index':vertex,'weight_zero_Euler_gauge':list(map(str,L0)),
            'full_positive_weight_Euler_gauge':list(map(str,L)),'full_original_M1_certificate':c1,'all_vertex_independence':not bad}
        if bad:rec.update(nonzero_vertex_derivatives=bad,source_covered=False,G_closed=False);return rec,[]
        Q=gauged.extract(indices,indices);W=s.Matrix([gauged[vertex,i]for i in indices]);ell=-s.Poly(W[2],z,w).coeff_monomial(z*z)*z-s.Poly(W[3],z,w).coeff_monomial(w*w)*w
        W=(W+ell*y).applyfunc(s.expand);variables=list(tt)+v2;Q0=Q.subs(dict.fromkeys(variables,0));W0=W.subs(dict.fromkeys(variables,0))
        qcoords=s.Matrix([s.Poly(Q[i,j]-Q0[i,j],*y).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in mons]);J=qcoords.jacobian(variables)
        assert all(s.expand(v)==0 for v in qcoords-J*s.Matrix(variables));qc,qminor,_=certify(J,domain);cols=qc['minor_columns'];uu=s.symbols('Q0:'+str(len(cols)));Qcan=Q0.copy()
        for k,(i,j)in enumerate(combinations(range(4),2)):
            Qcan[i,j]+=sum(J[10*k+r,jj]*s.prod(v**e for v,e in zip(y,mon))*q for jj,q in zip(cols,uu)for r,mon in enumerate(mons));Qcan[j,i]=-Qcan[i,j]
        Wcan=W0+s.Matrix([z*(ww[0]*base0+ww[1]*base1)+w*(ww[2]*base0+ww[3]*base1)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*base0+ww[8]*base1)+w*(ww[9]*base0+ww[10]*base1)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w]);wc=[];wd=W-W0
        for i in range(2):wc.extend([s.Poly(wd[i],*y).coeff_monomial(mon)for mon in[base0*z,base1*z,base0*w,base1*w,z*z,z*w,w*w]])
        wc.extend([s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)])
        assert all(s.expand(v)==0 for v in W-Wcan.subs(dict(zip(ww,wc)),simultaneous=True));Wmap=s.Matrix(wc).jacobian(variables)
        strict={(i,j,k):s.expand(sum(Qcan[i,l]*s.diff(Qcan[j,k],y[l])+Qcan[j,l]*s.diff(Qcan[k,i],y[l])+Qcan[k,l]*s.diff(Qcan[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        Pi=s.zeros(5)
        for i,j in combinations(range(4),2):Pi[indices[i],indices[j]]=Qcan[i,j];Pi[indices[j],indices[i]]=-Qcan[i,j]
        for i in range(4):Pi[vertex,indices[i]]=Wcan[i];Pi[indices[i],vertex]=-Wcan[i]
        D=divisor(Pi);dp=s.Poly(D,*x);D3=s.expand(sum(v*s.prod(t**e for t,e in zip(x,mon))for mon,v in dp.terms()if mon[3]+mon[4]==3))
        Ann=matrix(old['full_line_annihilator']);leading=s.Matrix([dp.coeff_monomial(mon)for mon in packet['all24_D3_monomials']]);line=(Ann*leading).applyfunc(s.expand)
        constant=[]
        for v in range(4):
            if all(s.expand(Qcan[i,j])==0 for i,j in combinations(range(4),2)if v not in[i,j])and all(s.diff(Qcan[i,j],y[v])==0 for i,j in combinations(range(4),2)):constant.append(str(y[v]))
        rec.update(source_covered=True,remaining_Euler_c_component=str(ell),all_Q_coordinates=list(map(str,y)),full_Q0={str((i,j)):str(Q0[i,j])for i,j in combinations(range(4),2)},
            full_W0=list(map(str,W0)),full60_Q_source_matrix=serial_matrix(J),full_W18_source_matrix=serial_matrix(Wmap),Q_space_certificate=qc,
            full_affine_Q={str((i,j)):str(Qcan[i,j])for i,j in combinations(range(4),2)},full_strict_Q_Jacobi={str(k):str(v)for k,v in strict.items()},
            full_W18_cover=True,full_actual_D_polynomial=str(D),full_actual_D3_polynomial=str(s.factor(D3)),
            exact_nonzero_leading_target=old['exact_nonzero_leading_target'],full23_D3_line_equations=list(map(str,line)),
            entire_Q_constant_factor_directions=constant,G_closed=bool(constant))
        den=[qminor];den.extend(s.denom(s.cancel(v))for v in list(K)+list(J));return rec,den
    results=[]
    for old in selected:
        variable=old['whole_alpha_line'];rec,den=analyze(old,domain=FIELD if variable else s.QQ);rec.update(type=old['type'],whole_alpha_line=variable,whole_fixed_plane_symbol_dimension_upper=old['whole_fixed_plane_symbol_dimension_upper']);walls=[]
        if variable:
            exception=s.Poly(parse(old['all_exception_polynomial']),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            for pol,power in s.factor_list(exception.as_expr())[1]:
                pp=s.Poly(pol,TAU);alpha=-pp.nth(0)/pp.nth(1)if pp.degree()==1 else s.CRootOf(pp.as_expr(),0)
                known=next((v for v in old['all_exceptional_sources']if parse(v['alpha'])==alpha),None)
                if known is not None and known['whole_m3_source_empty']:walls.append({'alpha':str(alpha),'whole_m3_source_empty_by_original_certificate':True,'G_closed':True});continue
                domain=s.QQ if pp.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(old,alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(pol),all_complex_embeddings_covered=True);walls.append(wall)
            rec.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,all_strata_G_closed=rec['G_closed']and all(v['G_closed']for v in walls))
        else:rec['all_strata_G_closed']=rec['G_closed']
        results.append(rec);print(old['type'],'covered',rec['source_covered'],'Qdim',rec.get('Q_space_certificate',{}).get('rank'),'factor',rec.get('entire_Q_constant_factor_directions'),
            'walls',[(v['alpha'],v.get('Q_space_certificate',{}).get('rank'))for v in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_MOVING_SIMPLE_ORE_ENTRANCE.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_moving_simple_frozen_jobs_20261107.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all7_G_closed':all(v['all_strata_G_closed']for v in results),'cases':results,
        'scope':'Whole2alpha+5fixed residualmoving-simplem3 sources: fullweightzeroEuler vertexandpositivefullkernel+FULLP2, exactQspace/W18cover/allstrictQJac/allactualD3line. No repeatednegativefactor-rowspan/jet/GB. FurthernonlinearGsupportmayremain.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_moving_simple_modules_main():
    """Three bad-column-zero exclusions, four complete source bridges, ten modules."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_MODULES_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole10affine modules finite budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_ORE_20261108.json';whole=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_20261107.json';packet=json.loads(whole.read_text());prior=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(prior['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==40
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;y=s.Matrix([a,b,z,w]);indices=[0,1,3,4];P1,v1,P2,v2=sources(x)
    ww=s.symbols('W0:18');qf=s.symbols('F0:3');u=s.Symbol('u');F=qf[0]*z*z+qf[1]*z*w+qf[2]*w*w
    k=-(3*TAU+5)/(3*TAU+8);R=(3*TAU+2)/(3*TAU-4);S=-(6*TAU+1)/(3*TAU-4)
    def makeQ(ab,az,aw,bz,bw):
        Q=s.zeros(4)
        for i,j,v in[(0,1,ab),(0,2,az),(0,3,aw),(1,2,bz),(1,3,bw)]:Q[i,j]=v;Q[j,i]=-v
        return Q
    def QHA(u,F):return makeQ(F,0,u*k*z*z,u*z*z,2*u*z*w)
    def QHB(u,F):return makeQ(u*(R*a*w+b*z)+F,0,u*z*w,-u*R*z*w,u*S*w*w)
    def QTR(u,F):return makeQ(F,0,-u*z*z,u*z*z,u*(2*z*w+z*z))
    def QCB(u,F):return makeQ(F,u*z*z,2*u*z*w,0,-u*z*z)
    WHA=s.Matrix([a*a,-2*a*b,-2*a*z,3*TAU*a*w+(3*TAU+5)*b*z]);WHB=s.Matrix([a*a,-2*a*b,a*z,3*TAU*a*w+(3*TAU+2)*b*z])
    WTR=s.Matrix([0,-3*a*a,0,-3*a*w+3*(a-b)*z]);WCB=s.Matrix([a*a,-2*a*b,4*a*z,10*a*w-3*a*z-3*b*w])
    def Wbasis(W0):
        return W0+s.Matrix([z*(ww[0]*a+ww[1]*b)+w*(ww[2]*a+ww[3]*b)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*a+ww[8]*b)+w*(ww[9]*a+ww[10]*b)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w])
    def pad(Q,W):
        P=s.zeros(5)
        for i,j in combinations(range(4),2):P[indices[i],indices[j]]=Q[i,j];P[indices[j],indices[i]]=-Q[i,j]
        for i in range(4):P[2,indices[i]]=W[i];P[indices[i],2]=-W[i]
        return P
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    exclusions=[]
    for name,alpha in[('DEG1_HA_CA_OTHER',-s.Rational(8,3)),('DEG1_HB_CA_OTHER',s.Rational(4,3)),('DEG2_CA2_0',None)]:
        old=next(v for v in packet['cases']if v['type']==name)
        if alpha is not None:old=next(v for v in old['all_exceptional_sources']if parse(v['alpha'])==alpha)
        K=s.Matrix([[parse(v)for v in row]for row in old['weight1_certificate']['complete_kernel']]);assert K[-3:,:]==s.zeros(3,K.cols)
        assert parse(old['full150_P0']['(3, 4)'])==0
        for i in range(3):assert s.rem(s.Poly(parse(old['full150_P0'][str((i,3))]),z),s.Poly(z,z)).is_zero
        exclusions.append({'type':name,'alpha':None if alpha is None else str(alpha),'complete_P1_kernel':serial_matrix(K),'all_three_bad_normal_rows_zero':True,'proper_D_factor':'z','whole_irreducible_source_empty':True})
    coverage=[]
    for name,Qfun,W0,variable,extract in[('DEG1_HA_CA_INTERSECTION',QHA,WHA,True,'bz'),('DEG1_HB_CA_INTERSECTION',QHB,WHB,True,'aw'),
                                      ('DEG1_TRIPLE_EXACT_B',QTR,WTR,False,'bz'),('DEG2_CB2_1',QCB,WCB,False,'az')]:
        old=next(v for v in packet['cases']if v['type']==name);K=s.Matrix([[parse(v)for v in row]for row in old['weight1_certificate']['complete_kernel']]);tt=s.symbols('p0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=parse(old['full150_P0'][str((i,j))]);P0[j,i]=-P0[i,j]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(tt))),simultaneous=True)+P2;L=s.Matrix([-s.diff(P[0,2],c),-s.diff(P[1,2],c),0,0,0]);gauged=(P-H.wedge(x,L)).applyfunc(s.expand)
        assert all(s.diff(gauged[i,j],c)==0 for i,j in H.PAIRS);Q=gauged.extract(indices,indices);W=s.Matrix([gauged[2,i]for i in indices])
        ell=-s.Poly(W[2],z,w).coeff_monomial(z*z)*z-s.Poly(W[3],z,w).coeff_monomial(w*w)*w;W=(W+ell*y).applyfunc(s.expand)
        noncommon=s.cancel(Q[1,2]/(z*z))if extract=='bz'else s.cancel(Q[0,3]/(z*w))if extract=='aw'else s.cancel(Q[0,2]/(z*z))
        central=s.expand(Q[0,1]-(noncommon*(R*a*w+b*z)if Qfun==QHB else 0))
        assert all(s.cancel(v)==0 for v in Q-Qfun(noncommon,central));wd=W-W0;wc=[]
        for i in range(2):wc.extend([s.Poly(wd[i],*y).coeff_monomial(mon)for mon in[a*z,b*z,a*w,b*w,z*z,z*w,w*w]])
        wc.extend([s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)])
        assert all(s.expand(v)==0 for v in W-Wbasis(W0).subs(dict(zip(ww,wc)),simultaneous=True));coords=[noncommon]+[s.Poly(central,z,w).coeff_monomial(mon)for mon in[z*z,z*w,w*w]]+wc
        J=s.Matrix(coords).jacobian(list(tt)+v2);det=s.factor(J.det(method='domain-ge'));assert J.rows==J.cols==22 and det!=0 and not s.sympify(det).free_symbols
        Qs=Qfun(u,F);strict={(i,j,l):s.cancel(sum(Qs[i,r]*s.diff(Qs[j,l],y[r])+Qs[j,r]*s.diff(Qs[l,i],y[r])+Qs[l,r]*s.diff(Qs[i,j],y[r])for r in range(4)))for i,j,l in combinations(range(4),3)};assert all(v==0 for v in strict.values())
        D=divisor(pad(Qs,Wbasis(W0)));dp=s.Poly(D,*x);D3=s.Matrix([dp.coeff_monomial(mon)for mon in packet['all24_D3_monomials']]);Ann=s.Matrix([[parse(v)for v in row]for row in old['full_line_annihilator']]);assert all(s.cancel(v)==0 for v in Ann*D3)
        target=s.Poly(parse(old['exact_nonzero_leading_target']),*x);idx=next(i for i,mon in enumerate(packet['all24_D3_monomials'])if target.coeff_monomial(mon)!=0);amplitude=s.factor(D3[idx]/target.coeff_monomial(packet['all24_D3_monomials'][idx]));assert amplitude!=0 and amplitude.free_symbols<={u,TAU}
        if variable:assert all(wall['weight1_certificate']['nullity']==13 for wall in old['all_exceptional_sources']if not wall['whole_m3_source_empty'])
        coverage.append({'type':name,'full_Euler_gauge':list(map(str,L)),'full_joint22_coordinate_matrix':serial_matrix(J),'nonzero_constant_joint_determinant':str(det),
            'whole_Q_W18_bijection':True,'full_strict_Q_Jacobi_zero':True,'exact_D3_target_amplitude':str(amplitude),'all_nonempty_old_walls_same_complete13_source':True})
    aa,bb,cc,dd=s.symbols('sa sb sc sd');shift=s.Matrix([aa*z+bb*w,cc*z+dd*w,0,0]);T=y+shift;Tinv=y-shift;Jac=T.jacobian(y);shears=[]
    for name,Qfun,delta in[('HA',QHA,u*((k*dd-aa)*z*z-2*bb*z*w)),('HB',QHB,u*(-cc*z*z-(S+R)*bb*w*w)),
                            ('TR',QTR,u*(-(dd+aa+bb)*z*z-2*bb*z*w)),('CB',QCB,u*((cc+bb)*z*z+2*dd*z*w))]:
        Q=Qfun(u,F);transport=Jac*Q.subs(dict(zip(y,Tinv)),simultaneous=True)*Jac.T
        assert all(s.cancel(v)==0 for v in transport-Qfun(u,F+delta));shears.append({'type':name,'full_actual_central_shear_identity_zero':True,'central_change':str(s.factor(delta))})
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list();kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());minor=dm.extract(list(rows),list(cols)).det();assert minor!=domain.zero or not cols
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(minor))),'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    def analyze(M,C,Dc,Dl,domain,offset):
        cm,mi,K=certify(M,domain);ca,am,KA=certify(M.row_join(-C),domain);rec={'full350_matrix':serial_matrix(M),'full350_constant':serial_matrix(C),'linear_certificate':cm,'augmented_certificate':ca,'image_offset':offset};den=[mi,am]
        if ca['rank']>cm['rank']:rec.update(status='inconsistent',G_closed=True,G_empty=True);return rec,den
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.row_join(-C).tolist()],(350,19),domain);rr,pivots=dm.rref();raw=rr.to_list();part=s.zeros(18,1)
        for i,j in enumerate(pivots):part[j]=domain.to_sympy(raw[i][18])
        assert all(domain.from_sympy(v)==domain.zero for v in M*part+C)
        D0=s.Poly(s.expand(Dc+sum(Dl[i]*part[i]for i in range(18))),*x,domain=domain);directions=[s.Poly(s.expand(sum(Dl[i]*K[i,j]for i in range(18))),*x,domain=domain)for j in range(K.cols)]
        common=D0
        for poly in directions:
            if not poly.is_zero:common=s.polys.polytools.gcd(common,poly).monic()
        factor=None;zero=common.is_zero
        if not zero and 0<common.total_degree()<5:factor=common.monic()
        elif not zero and common.total_degree()==5:
            factor=next((poly.monic()for poly,power in common.factor_list()[1]if 0<poly.total_degree()<5),None)
        rec.update(status='complete_affine_module',particular_solution=list(map(str,part)),full_D_particular=str(D0.as_expr()),
            full_D_kernel_coefficients=[str(poly.as_expr())for poly in directions],common_D_gcd=str(common.as_expr()),proper_D_factor=None if factor is None else str(factor.as_expr()),
            whole_D_zero=zero,W_kernel_dimension=K.cols,whole_projective_image_upper=offset+K.cols,G_empty=zero or factor is not None,G_closed=zero or factor is not None or offset+K.cols<=24)
        if domain==FIELD:
            den.extend(s.denom(s.cancel(v))for v in list(K)+list(part))
            if factor is not None:
                for poly in[factor]+[poly.exquo(factor)for poly in[D0]+directions]:den.extend(s.denom(s.cancel(v))for v in poly.coeffs())
        return rec,den
    cases=[];models=[]
    for eps in[0,1]:
        models.extend([('HA_'+str(eps),QHA(1,eps*w*w),WHA,True,17+eps,[-s.Rational(5,3),-s.Rational(8,3)]),
            ('HB_'+str(eps),QHB(1,eps*z*w),WHB,True,17+eps,[s.Rational(4,3),-s.Rational(2,3)]),
            ('HB_EXTRA_'+str(eps),QHB(1,eps*z*w+w*w).subs(TAU,s.Rational(1,3)),WHB.subs(TAU,s.Rational(1,3)),False,16+eps,[]),
            ('TRIPLE_'+str(eps),QTR(1,eps*w*w),WTR,False,16+eps,[]),('CB2_'+str(eps),QCB(1,eps*w*w),WCB,False,17+eps,[])])
    for name,Q,W0,variable,offset,excluded in models:
        P=pad(Q,Wbasis(W0));rows=s.Matrix([s.expand(v)for v in H.all350(H.projective(H.bilinear(P,P,x),x),x)]);C=rows.subs(dict.fromkeys(ww,0));M=rows.jacobian(ww)
        assert all(s.cancel(v)==0 for v in rows-C-M*s.Matrix(ww));D=divisor(P);Dc=D.subs(dict.fromkeys(ww,0));Dl=[s.expand(D).coeff(v)for v in ww];assert s.expand(D-Dc-sum(v*q for v,q in zip(Dl,ww)))==0
        rec,den=analyze(M,C,Dc,Dl,FIELD if variable else s.QQ,offset);rec.update(type=name,full_Q={str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},W0=list(map(str,W0)))
        if variable:
            exception=s.Poly(s.prod(TAU-v for v in excluded),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            walls=[]
            for pol,power in s.factor_list(exception.as_expr())[1]:
                p=s.Poly(pol,TAU);alpha=-p.nth(0)/p.nth(1)if p.degree()==1 else s.CRootOf(p.as_expr(),0)
                if alpha in excluded:walls.append({'alpha':str(alpha),'source_empty_by_original_full_D3_image':True,'G_closed':True});continue
                domain=s.QQ if p.degree()==1 else s.QQ.algebraic_field(alpha);sub=lambda v:v.subs(TAU,alpha);wall,_=analyze(sub(M),sub(C),sub(Dc),[sub(v)for v in Dl],domain,offset-1)
                wall.update(alpha=str(alpha),minimal_polynomial=str(pol),all_complex_embeddings_covered=True);walls.append(wall)
            rec.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,all_strata_G_closed=rec['G_closed']and all(v['G_closed']for v in walls))
        else:rec['all_strata_G_closed']=rec['G_closed']
        cases.append(rec);print(name,rec['status'],'k',rec['linear_certificate']['nullity'],'factor',rec.get('proper_D_factor'),'G',rec['all_strata_G_closed'],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_MOVING_SIMPLE_TEN_MODULE_REDUCTION.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_moving_simple_ore_frozen_jobs_20261108.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all_ten_modules_G_closed':all(v['all_strata_G_closed']for v in cases),'three_full_bad_column_zero_sources':exclusions,
        'full_four_source_coverage':coverage,'actual_full_central_shear_identities':shears,'cases':cases,
        'scope':'Wholemoving-simplem3 residualseven: threefullbadcolumnzero zfactor, fourfullQ/W22bijections and10actualcanonicalW18modules withallalphawalls. Iftrue closesbothprimitive01/02 G; generictriple/m2/reducedRoutside. OldCA2vertexbfalseobstructionkeptbutnotneededforG.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,whole,proof,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_triple_main():
    """Whole90 generictriple firstconnections, badnormalcolumn factor first."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole90 triple firstconnections finite budget600s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_MODULES_20261108.json';old=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(old['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==41
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;base=[a,b,c];normals=[z,w];P1,v1,P2,v2=sources(x);kk=s.symbols('K0:3')
    monG=[v+(j,3-j)for v in H.monomials(3,2)for j in range(4)];monK=[v+(j,2-j)for v in H.monomials(3,3)for j in range(3)]
    def dlog(g):return s.Matrix([s.cancel(s.diff(g,t)/g)for t in base])
    def differential(g):return s.Matrix([s.diff(g,t)for t in base])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list();kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());minor=dm.extract(list(rows),list(cols)).det()
        assert minor!=domain.zero or not cols;K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(minor))),
            'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    models=[];xx=b/a;yy=c/a;dx=differential(xx);dy=differential(yy);lx=dlog(xx)
    for e in range(3):
        f=a*a*b;h=a**e*b**(2-e);eta=s.Rational(4,3)-e-TAU;p=dlog(f)/3+TAU*dlog(a/b)+dx;q=2*p-dlog(f/h)
        for nm,rr in[('ZERO',s.zeros(3,1)),('LOG',lx)]:models.append(('E'+str(e)+'_EXACT_B_'+nm,f,h,p,q,rr,10,True,[]))
        alpha=s.Rational(1,3)-e;pv=p.subs(TAU,alpha);qv=q.subs(TAU,alpha)
        models.append(('E'+str(e)+'_EXACT_B_RESONANT',f,h,pv,qv,(yy*dx+dy)/xx,10,False,[]))
        p=dlog(f)/3+TAU*dlog(a/b);q=2*p-dlog(f/h);excluded=[s.Rational(4,3)-e-v for v in[-2,-1,0,1]]
        for nm,rr in[('ZERO',s.zeros(3,1)),('DX',dx),('Y',eta*yy*lx+dy)]:models.append(('E'+str(e)+'_LOG_GENERIC_'+nm,f,h,p,q,rr,10,True,excluded))
        special=[(-1,[s.zeros(3,1),dx,-yy*lx+dy,-yy*lx+dy+dx],9),
            (0,[s.zeros(3,1),dx,lx,lx+dx,dy,lx+dy],9),
            (1,[differential(psi)/xx for psi in[xx*yy,xx*yy+yy,s.Integer(0),xx*xx,yy,xx*xx+yy]],11),
            (-2,[xx*xx*differential(psi/(xx*xx))for psi in[xx*yy,xx*yy+yy,yy,s.Integer(0),xx]],11)]
        for ev,rows,Bdim in special:
            alpha=s.Rational(4,3)-e-ev;pv=p.subs(TAU,alpha);qv=q.subs(TAU,alpha)
            for j,rr in enumerate(rows):models.append(('E'+str(e)+'_LOG_ETA'+str(ev)+'_'+str(j),f,h,pv,qv,rr,Bdim,False,[]))
    f=a**3;h=a*a
    for j,(psi,Bdim)in enumerate([(s.Integer(0),4),(xx,6),(xx*xx,7),(xx*xx+yy,8)]):models.append(('TRIPLE_PHI_ZERO_'+str(j),f,h,dlog(a),dlog(a),differential(psi),Bdim,False,[]))
    for j,(rr,Bdim)in enumerate([(s.zeros(3,1),7),(xx*dx,8),(yy*dx+dy,9)]):models.append(('TRIPLE_PHI_LINEAR_'+str(j),f,h,dlog(a)+dx,dlog(a)+2*dx,rr,Bdim,False,[]))
    for j,(rr,Bdim)in enumerate([(s.zeros(3,1),8),(dx,9)]):models.append(('TRIPLE_PHI_SQUARE_'+str(j),f,h,dlog(a)+2*xx*dx,dlog(a)+4*xx*dx,rr,Bdim,False,[]))
    assert len(models)==90 and sum(row[-2]for row in models)==15
    records=[]
    checkpoint=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CHECKPOINT_20261108.json'
    if checkpoint.exists():raise FileExistsError(checkpoint)
    for name,f,h,p,q,r,Bdim,variable,excluded in models:
        beta=[[p,s.zeros(3,1)],[r,q]]
        for i,j in combinations(range(3),2):
            assert s.cancel(s.diff(p[j],base[i])-s.diff(p[i],base[j]))==0
            assert s.cancel(s.diff(q[j],base[i])-s.diff(q[i],base[j]))==0
            assert s.cancel(s.diff(r[j],base[i])-s.diff(r[i],base[j])+(q[i]-p[i])*r[j]-(q[j]-p[j])*r[i])==0
        Pbase=s.zeros(5)
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:Pbase[i,j]=s.diff(f,base[k]);Pbase[j,i]=-Pbase[i,j]
        omegas=[]
        for i in range(2):
            for j in range(2):
                om=(3*f*beta[i][j]-(differential(f)if i==j else s.zeros(3,1))).applyfunc(s.cancel)
                assert s.expand(sum(base[k]*om[k]for k in range(3)))==0
                assert all(s.denom(v)==1 for v in om)
                V=s.Matrix([s.diff(om[2],b)-s.diff(om[1],c),s.diff(om[0],c)-s.diff(om[2],a),s.diff(om[1],a)-s.diff(om[0],b)])/3
                assert s.expand(sum(s.diff(V[k],base[k])for k in range(3)))==0
                inverse=s.Matrix([c*V[1]-b*V[2],a*V[2]-c*V[0],b*V[0]-a*V[1]])
                assert all(s.expand(v)==0 for v in inverse-om)
                for k in range(3):Pbase[k,3+i]-=V[k]*normals[j];Pbase[3+i,k]+=V[k]*normals[j]
                omegas.append({'i':i,'j':j,'beta':list(map(str,beta[i][j])),'horizontal_polynomial_form':list(map(str,om)),'inverse_tracefree_field':list(map(str,V))})
        target=s.expand(h*z**3);target_vector=s.Matrix([s.Poly(target,*x).coeff_monomial(mon)for mon in monG]);idx=next(i for i,v in enumerate(target_vector)if v!=0)
        Ann=s.zeros(23,24)
        for row,i in enumerate(j for j in range(24)if j!=idx):Ann[row,i]=1;Ann[row,idx]=-target_vector[i]/target_vector[idx]
        assert Ann*target_vector==s.zeros(23,1)
        def analyze(alpha=None,domain=FIELD):
            sub=lambda v:v if alpha is None else v.subs(TAU,alpha)
            Pk=sub(Pbase);Pk[3,4]=kk[0]*z*z+kk[1]*z*w+kk[2]*w*w;Pk[4,3]=-Pk[3,4]
            Dk=divisor(Pk);poly=s.Poly(Dk,*x);rows=s.Matrix([poly.coeff_monomial(mon)for mon in monK]);KM=rows.jacobian(kk);const=rows.subs(dict.fromkeys(kk,0))
            kc,km,_=certify(KM,domain);ka,am,_=certify(KM.row_join(-const),domain);den=[km,am]
            record={'full30_global_K_matrix':serial_matrix(KM),'full30_global_K_constant':serial_matrix(const),'K_certificate':kc,'K_augmented_certificate':ka}
            if ka['rank']>kc['rank']:record.update(global_K_inconsistent=True,whole_m3_source_empty=True,G_closed=True);return record,den
            sol=list(s.linsolve((KM,-const),kk));assert len(sol)==1 and kc['rank']==3
            P0=sub(Pbase);P0[3,4]=s.expand(sum(v*m for v,m in zip(sol[0],[z*z,z*w,w*w])));P0[4,3]=-P0[3,4]
            assert all(domain.from_sympy(v)==domain.zero for v in s.Poly(divisor(P0),*x).coeffs())
            Jac0=[s.expand(v)for v in H.all350(H.projective(H.bilinear(P0,P0,x),x),x)];bad=[(i,v)for i,v in enumerate(Jac0)if domain.from_sympy(v)!=domain.zero]
            record.update(unique_global_K=str(P0[3,4]),full150_P0={str((i,j)):str(P0[i,j])for i,j in H.PAIRS})
            if bad:
                record.update(P0_nonzero_Jacobi=[(i,str(v))for i,v in bad],whole_m3_source_empty=True,G_closed=True);den.append(bad[0][1]);return record,den
            D1=operator(P0,P1,v1,x);D2=operator(P0,P2,v2,x);full=divisor(P0+P1+P2);pp=s.Poly(full,*x)
            D3=s.Matrix([pp.coeff_monomial(mon)for mon in monG]);L=D3.jacobian(v1);assert all(s.expand(v)==0 for v in D3-L*s.Matrix(v1))
            M1=D1.col_join(Ann*L)
            c1,m1,K1=certify(M1,domain);c2,m2,K2=certify(D2,domain);image=(L*K1).applyfunc(s.cancel)
            qbad=s.cancel(3*f*sum(sub(q[i])*v1[-3+i]for i in range(3)));qbpoly=s.Poly(qbad,*base)
            QM=s.Matrix([qbpoly.coeff_monomial(mon)for mon in H.monomials(3,2)]).jacobian(v1)
            assert all(domain.from_sympy(v)==domain.zero for v in QM*K1)
            assert all(domain.from_sympy(v)==domain.zero for v in Ann*image)
            leading=[domain.from_sympy(image[idx,j]/target_vector[idx])for j in range(image.cols)];nonzero=[j for j,v in enumerate(leading)if v!=domain.zero]
            offset=Bdim+5-(1 if alpha is not None and variable else 0);upper=offset+K1.cols+K2.cols
            badzero=all(domain.from_sympy(v)==domain.zero for v in K1[-3:,:])
            assert all(domain.from_sympy(v)==domain.zero for i in range(5)for v in s.Poly(P0[i,3].subs(z,0),*x).coeffs())
            record.update(P0_full350_Jacobi_zero_and_actualD_zero=True,full350_weight1_operator=serial_matrix(D1),full350_weight2_operator=serial_matrix(D2),
                full24by30_actual_D3_map=serial_matrix(L),full_weight1_plus_line_matrix=serial_matrix(M1),
                weight1_certificate=c1,weight2_certificate=c2,whole_D3_kernel_image=serial_matrix(image),
                all_three_bad_normal_kernel_rows_zero=badzero,proper_D_factor_by_bad_normal_rows='z'if badzero else None,
                full6by30_constant_root_character_matrix=serial_matrix(QM),all_full_kernel_q_Cbad_zero=True,
                nonzero_D3_line_certificate=None if not nonzero else {'kernel_column':nonzero[0],'coefficient':str(domain.to_sympy(leading[nonzero[0]]))},
                D3_image_rank=int(bool(nonzero)),whole_m3_source_empty=not nonzero,whole_irreducible_source_empty=not nonzero or badzero,whole_projective_image_upper=upper,G_closed=not nonzero or badzero or upper<=24)
            den.extend([m1,m2]);den.extend(s.denom(s.cancel(v))for v in list(K1)+list(K2)+list(sol[0]))
            return record,den
        result,den=analyze(domain=FIELD if variable else s.QQ);result.update(type=name,whole_alpha_line=variable,full_connection_inverse_forms=omegas,
            whole_fixed_plane_symbol_dimension_upper=Bdim,exact_nonzero_leading_target=str(target),full_line_annihilator=serial_matrix(Ann))
        if variable:
            exception=s.Poly(s.prod(TAU-v for v in excluded),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            walls=[]
            for factor,power in s.factor_list(exception.as_expr())[1]:
                poly=s.Poly(factor,TAU);alpha=-poly.nth(0)/poly.nth(1)if poly.degree()==1 else s.CRootOf(poly.as_expr(),0)
                if alpha in excluded:
                    walls.append({'alpha':str(alpha),'whole_m3_source_empty':False,'full_special_source_elsewhere':True,'G_closed':True});continue
                domain=s.QQ if poly.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(factor),all_complex_embeddings_covered=True);walls.append(wall)
            result.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,
                all_strata_G_closed=result['G_closed']and all(v['G_closed']for v in walls))
        else:result['all_strata_G_closed']=result['G_closed']
        records.append(result)
        if len(records)%10==0:checkpoint.write_text(json.dumps({'task_id':'G','complete':False,'completed_source_count':len(records),'target_source_count':90,'cases':records,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text()},indent=2)+'\n')
        print(name,'empty',result.get('whole_irreducible_source_empty',result['whole_m3_source_empty']),'k',result.get('weight1_certificate',{}).get('nullity'),
            result.get('weight2_certificate',{}).get('nullity'),'upper',result.get('whole_projective_image_upper'),
            'walls',[(v['alpha'],v['whole_m3_source_empty'],v.get('whole_projective_image_upper'))for v in result.get('all_exceptional_sources',[])],flush=True)
    proofs=[ROOT/'research/classification'/name for name in['PLANE_NONREDUCED_M3_TRIPLE_ROOT_FIRST_CONNECTIONS.md','PLANE_NONREDUCED_M3_TRIPLE_FINITE_NORMALFORMS.md','PLANE_NONREDUCED_M3_TRIPLE_CANONICAL_REDUCTION.md','PLANE_NONREDUCED_M3_TRIPLE_IRREGULAR_C_EXCLUSION.md','PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md','PLANE_M3_CONSTANT_ROOT_OBSTRUCTION.md']]
    archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_triple_preappend_inputs_20261108.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all90_G_closed':all(v['all_strata_G_closed']for v in records),'cases':records,
        'scope':'Entiregenerictripleroot nonredR firstconnectionsafterpurefullsourceobstructions:15wholealpha+75fixed. Allglobalclosedflatconnections/inverseforms/K/P0Jac/full30/9/D3line/allbadnormalrows; fullspecialfibers retained. Badcolumnzero impliesproperzfactor. Higher/Gonlyifexplicittrue, m2/reducedRoutside.',
        'all24_D3_monomials':monG,'all30_K_monomials':monK,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,archive]+proofs],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')
    if checkpoint.exists():checkpoint.unlink()


def nonreduced_m3_triple_ore_main():
    """All14 remaininggenerictriple sources: fullEuler/Q/W, universalconstant-v test."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_ORE_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('whole14triple Ore-source finite budget240s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_20261108.json';packet=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(packet['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==42
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;P1,v1,P2,v2=sources(x);ww=s.symbols('W0:18');mons=H.monomials(4,2)
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    def matrix(raw):return s.Matrix([[parse(v)for v in row]for row in raw])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1];free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list()
        kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());mi=dm.extract(list(rows),list(cols)).det();assert mi!=domain.zero or not cols
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(mi))),'complete_kernel':serial_matrix(K)},domain.to_sympy(mi),K
    selected=[]
    for e in range(3):
        row=dict(next(v for v in packet['cases']if v['type']=='E'+str(e)+'_LOG_GENERIC_DX'));row['whole_fixed_plane_symbol_dimension_upper']=9;selected.append(row)
    row=dict(next(v for v in packet['cases']if v['type']=='E2_LOG_GENERIC_ZERO'));row['whole_fixed_plane_symbol_dimension_upper']=8;selected.append(row)
    for suffix,Bdim in[('ZERO',8),('LOG',9)]:
        outer=next(v for v in packet['cases']if v['type']=='E0_EXACT_B_'+suffix);row=dict(next(v for v in outer['all_exceptional_sources']if parse(v['alpha'])==-s.Rational(2,3)))
        row.update(type='EXACT_B_FIXED_'+suffix,whole_alpha_line=False,whole_fixed_plane_symbol_dimension_upper=Bdim,exact_nonzero_leading_target=outer['exact_nonzero_leading_target'],full_line_annihilator=outer['full_line_annihilator']);selected.append(row)
    for e in range(3):
        row=dict(next(v for v in packet['cases']if v['type']=='E'+str(e)+'_LOG_ETA0_3'));row['whole_fixed_plane_symbol_dimension_upper']=8;selected.append(row)
    row=dict(next(v for v in packet['cases']if v['type']=='E2_LOG_ETA0_2'));row['whole_fixed_plane_symbol_dimension_upper']=7;selected.append(row)
    for nm in['TRIPLE_PHI_ZERO_1','TRIPLE_PHI_LINEAR_0','TRIPLE_PHI_LINEAR_1','TRIPLE_PHI_SQUARE_1']:selected.append(dict(next(v for v in packet['cases']if v['type']==nm)))
    assert len(selected)==14 and sum(v['whole_alpha_line']for v in selected)==4
    def analyze(old,alpha=None,domain=FIELD):
        sub=lambda v:v if alpha is None else v.subs(TAU,alpha);vertex=2;indices=[0,1,3,4];y=s.Matrix([a,b,z,w]);base0,base1=y[:2];L0=s.zeros(5,1)
        M1=sub(matrix(old['full_weight1_plus_line_matrix']));c1,_,K=certify(M1,domain);tt=s.symbols('p0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=sub(parse(old['full150_P0'][str((i,j))]));P0[j,i]=-P0[i,j]
        assert P0[3,4]==0 and all(s.diff(P0[i,j],c)==0 for i,j in H.PAIRS)
        if all(domain.from_sympy(v)==domain.zero for v in K[-3:,:]):
            return {'source_covered':True,'full_original_M1_certificate':c1,'all_bad_normal_rows_zero':True,'Ore_not_needed_because_irreducible_source_empty':True,'proper_D_factor':'z','G_closed':True},[]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(tt))),simultaneous=True)+P2-H.wedge(x,L0)
        L=s.Matrix([-s.diff(P[i,vertex],x[vertex])if i!=vertex else 0 for i in range(5)])
        gauged=(P-H.wedge(x,L)).applyfunc(s.expand);bad=[(str((i,j)),str(s.diff(gauged[i,j],x[vertex])))for i,j in H.PAIRS if s.diff(gauged[i,j],x[vertex])!=0]
        rec={'alpha':None if alpha is None else str(alpha),'vertex_index':vertex,'weight_zero_Euler_gauge':list(map(str,L0)),
            'full_positive_weight_Euler_gauge':list(map(str,L)),'full_original_M1_certificate':c1,'all_vertex_independence':not bad}
        if bad:rec.update(nonzero_vertex_derivatives=bad,source_covered=False,G_closed=False);return rec,[]
        Q=gauged.extract(indices,indices);W=s.Matrix([gauged[vertex,i]for i in indices]);ell=-s.Poly(W[2],z,w).coeff_monomial(z*z)*z-s.Poly(W[3],z,w).coeff_monomial(w*w)*w
        W=(W+ell*y).applyfunc(s.expand);variables=list(tt)+v2;Q0=Q.subs(dict.fromkeys(variables,0));W0=W.subs(dict.fromkeys(variables,0))
        qcoords=s.Matrix([s.Poly(Q[i,j]-Q0[i,j],*y).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in mons]);J=qcoords.jacobian(variables)
        assert all(s.expand(v)==0 for v in qcoords-J*s.Matrix(variables));qc,qminor,_=certify(J,domain);cols=qc['minor_columns'];uu=s.symbols('Q0:'+str(len(cols)));Qcan=Q0.copy()
        for k,(i,j)in enumerate(combinations(range(4),2)):
            Qcan[i,j]+=sum(J[10*k+r,jj]*s.prod(v**e for v,e in zip(y,mon))*q for jj,q in zip(cols,uu)for r,mon in enumerate(mons));Qcan[j,i]=-Qcan[i,j]
        Wcan=W0+s.Matrix([z*(ww[0]*base0+ww[1]*base1)+w*(ww[2]*base0+ww[3]*base1)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*base0+ww[8]*base1)+w*(ww[9]*base0+ww[10]*base1)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w]);wc=[];wd=W-W0
        for i in range(2):wc.extend([s.Poly(wd[i],*y).coeff_monomial(mon)for mon in[base0*z,base1*z,base0*w,base1*w,z*z,z*w,w*w]])
        wc.extend([s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)])
        assert all(s.expand(v)==0 for v in W-Wcan.subs(dict(zip(ww,wc)),simultaneous=True));Wmap=s.Matrix(wc).jacobian(variables)
        joint=J.extract(qc['minor_rows'],list(range(J.cols))).col_join(Wmap);jc,jminor,JK=certify(joint,domain)
        vv=s.symbols('V0:4');V=s.Matrix(vv);cv=[]
        for i,j,k in combinations(range(4),3):cv.append(V[i]*Qcan[j,k]-V[j]*Qcan[i,k]+V[k]*Qcan[i,j])
        cv.extend(sum(V[k]*s.diff(Qcan[i,j],y[k])for k in range(4))for i,j in combinations(range(4),2))
        cvcoeff=[val for expr in cv for val in s.Poly(s.expand(expr),*y,*uu).coeffs()if val!=0]
        Cvm=s.Matrix(cvcoeff).jacobian(vv)if cvcoeff else s.zeros(0,4);vc,vminor,VK=certify(Cvm,domain)
        strict={(i,j,k):s.expand(sum(Qcan[i,l]*s.diff(Qcan[j,k],y[l])+Qcan[j,l]*s.diff(Qcan[k,i],y[l])+Qcan[k,l]*s.diff(Qcan[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        Pi=s.zeros(5)
        for i,j in combinations(range(4),2):Pi[indices[i],indices[j]]=Qcan[i,j];Pi[indices[j],indices[i]]=-Qcan[i,j]
        for i in range(4):Pi[vertex,indices[i]]=Wcan[i];Pi[indices[i],vertex]=-Wcan[i]
        D=divisor(Pi);dp=s.Poly(D,*x);D3=s.expand(sum(v*s.prod(t**e for t,e in zip(x,mon))for mon,v in dp.terms()if mon[3]+mon[4]==3))
        Ann=matrix(old['full_line_annihilator']);leading=s.Matrix([dp.coeff_monomial(mon)for mon in packet['all24_D3_monomials']]);line=(Ann*leading).applyfunc(s.expand)
        constant=[]
        for v in range(4):
            if all(s.expand(Qcan[i,j])==0 for i,j in combinations(range(4),2)if v not in[i,j])and all(s.diff(Qcan[i,j],y[v])==0 for i,j in combinations(range(4),2)):constant.append(str(y[v]))
        rec.update(source_covered=True,remaining_Euler_c_component=str(ell),all_Q_coordinates=list(map(str,y)),full_Q0={str((i,j)):str(Q0[i,j])for i,j in combinations(range(4),2)},
            full_W0=list(map(str,W0)),full60_Q_source_matrix=serial_matrix(J),full_W18_source_matrix=serial_matrix(Wmap),Q_space_certificate=qc,
            full_affine_Q={str((i,j)):str(Qcan[i,j])for i,j in combinations(range(4),2)},full_strict_Q_Jacobi={str(k):str(v)for k,v in strict.items()},
            full_W18_cover=True,full_actual_D_polynomial=str(D),full_actual_D3_polynomial=str(s.factor(D3)),
            exact_nonzero_leading_target=old['exact_nonzero_leading_target'],full23_D3_line_equations=list(map(str,line)),
            entire_Q_constant_factor_directions=constant,whole_QW_source_joint_matrix=serial_matrix(joint),joint_QW_certificate=jc,
            joint_QW_bijection=joint.rows==joint.cols==jc['rank'],full_common_constant_v_matrix=serial_matrix(Cvm),common_constant_v_certificate=vc,
            whole_Q_common_constant_v=None if VK.cols==0 else list(map(str,VK[:,0])),G_closed=VK.cols>0)
        den=[qminor,jminor,vminor];den.extend(s.denom(s.cancel(v))for v in list(K)+list(J)+list(JK)+list(VK));return rec,den
    results=[]
    for old in selected:
        variable=old['whole_alpha_line'];rec,den=analyze(old,domain=FIELD if variable else s.QQ);rec.update(type=old['type'],whole_alpha_line=variable,whole_fixed_plane_symbol_dimension_upper=old['whole_fixed_plane_symbol_dimension_upper']);walls=[]
        if variable:
            exception=s.Poly(parse(old['all_exception_polynomial']),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            for pol,power in s.factor_list(exception.as_expr())[1]:
                pp=s.Poly(pol,TAU);alpha=-pp.nth(0)/pp.nth(1)if pp.degree()==1 else s.CRootOf(pp.as_expr(),0)
                known=next((v for v in old['all_exceptional_sources']if parse(v['alpha'])==alpha),None)
                if known is not None and known.get('whole_m3_source_empty',False):walls.append({'alpha':str(alpha),'whole_m3_source_empty_by_original_certificate':True,'G_closed':True});continue
                domain=s.QQ if pp.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(old,alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(pol),all_complex_embeddings_covered=True);walls.append(wall)
            rec.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,all_strata_G_closed=rec['G_closed']and all(v['G_closed']for v in walls))
        else:rec['all_strata_G_closed']=rec['G_closed']
        results.append(rec);print(old['type'],'covered',rec['source_covered'],'Qdim',rec.get('Q_space_certificate',{}).get('rank'),'factor',rec.get('entire_Q_constant_factor_directions'),
            'walls',[(v['alpha'],v.get('Q_space_certificate',{}).get('rank'))for v in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_TRIPLE_ORE_ENTRANCE.md';partial=ROOT/'research/classification/PLANE_NONREDUCED_M3_TRIPLE_PARTIAL_BOUND.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_triple_ore_preappend_inputs_20261108.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all14_G_closed':all(v['all_strata_G_closed']for v in results),'cases':results,
        'scope':'Whole4alpha+10fixed remaininggenerictriple sources. Existingfull90M1/P0 datareused, nooperator/profile/jet rerun. Allspecialeta fibers restoreddirectoriginalmatrix, wholepositiveEuler/Q/Wsource/jointbijection/strictQJac/D3line anduniversalconstant-v criterion. FalsewholeGsources remainexplicit.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,partial,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')



def nonreduced_m3_triple_corrected_v2_main():
    """Correct-minus full r8 classification; reuse only exact old P0/operator certificates."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('corrected90 minus-curvature finite budget600s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(600);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_20261108.json';old=json.loads(previous.read_text())
    archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_triple_corrected_v2_preappend_inputs_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==43
    assert all(after[name]==value for name,value in before.items())
    helper_old=json.loads((ROOT/'research/classification/inputs/plane_nonreduced_m3_triple_ore_preappend_inputs_20261108.json').read_text())['helper_source']
    assert all(functions(HELPER.read_text())[name]==value for name,value in functions(helper_old).items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;base=[a,b,c];normals=[z,w];P1,v1,P2,v2=sources(x);kk=s.symbols('K0:3')
    monG=[v+(j,3-j)for v in H.monomials(3,2)for j in range(4)];monK=[v+(j,2-j)for v in H.monomials(3,3)for j in range(3)]
    def dlog(g):return s.Matrix([s.cancel(s.diff(g,t)/g)for t in base])
    def differential(g):return s.Matrix([s.diff(g,t)for t in base])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list();kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());minor=dm.extract(list(rows),list(cols)).det()
        assert minor!=domain.zero or not cols;K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(minor))),
            'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    def matrix(raw):return s.Matrix([[parse(v)for v in row]for row in raw])
    def p0_key(P,target):return(str(s.expand(target)),tuple(str(s.cancel(P[i,j]))for i,j in H.PAIRS))
    old_cache={}
    for outer in old['cases']:
        for record,wall in[(outer,None)]+[(v,v.get('alpha'))for v in outer.get('all_exceptional_sources',[])]:
            if 'full350_weight1_operator'not in record:continue
            P=s.zeros(5)
            for i,j in H.PAIRS:P[i,j]=parse(record['full150_P0'][str((i,j))]);P[j,i]=-P[i,j]
            target=parse(outer['exact_nonzero_leading_target'])
            old_cache[p0_key(P,target)]=(record,{'path':str(previous.relative_to(ROOT)),'source_type':outer['type'],'alpha':wall})
    reuse_fields=['full350_weight1_operator','full350_weight2_operator','full24by30_actual_D3_map','full_weight1_plus_line_matrix','weight1_certificate','weight2_certificate','whole_D3_kernel_image','all_three_bad_normal_kernel_rows_zero','proper_D_factor_by_bad_normal_rows','full6by30_constant_root_character_matrix','all_full_kernel_q_Cbad_zero','nonzero_D3_line_certificate','D3_image_rank','whole_m3_source_empty','whole_irreducible_source_empty','whole_projective_image_upper','G_closed']
    models=[];xx=b/a;yy=c/a;dx=differential(xx);dy=differential(yy);lx=dlog(xx)
    for e in range(3):
        f=a*a*b;h=a**e*b**(2-e);eta=s.Rational(4,3)-e-TAU;p=dlog(f)/3+TAU*dlog(a/b)+dx;q=2*p-dlog(f/h)
        for nm,rr in[('ZERO',s.zeros(3,1)),('LOG',lx)]:models.append(('E'+str(e)+'_EXACT_B_'+nm,f,h,p,q,rr,10,True,[]))
        alpha=s.Rational(7,3)-e;pv=p.subs(TAU,alpha);qv=q.subs(TAU,alpha)
        models.append(('E'+str(e)+'_EXACT_B_RESONANT',f,h,pv,qv,(-yy*dx+dy)/xx,10,False,[]))
        p=dlog(f)/3+TAU*dlog(a/b);q=2*p-dlog(f/h);excluded=[s.Rational(4,3)-e-v for v in[-1,0,1,2]]
        for nm,rr in[('ZERO',s.zeros(3,1)),('DX',dx),('Y',-eta*yy*lx+dy)]:models.append(('E'+str(e)+'_LOG_GENERIC_'+nm,f,h,p,q,rr,10,True,excluded))
        special=[(1,[s.zeros(3,1),dx,-yy*lx+dy,-yy*lx+dy+dx],9),
            (0,[s.zeros(3,1),dx,lx,lx+dx,dy,lx+dy],9),
            (-1,[differential(psi)/xx for psi in[xx*yy,xx*yy+yy,s.Integer(0),xx*xx,yy,xx*xx+yy]],11),
            (2,[xx*xx*differential(psi/(xx*xx))for psi in[xx*yy,xx*yy+yy,yy,s.Integer(0),xx]],11)]
        for ev,rows,Bdim in special:
            alpha=s.Rational(4,3)-e-ev;pv=p.subs(TAU,alpha);qv=q.subs(TAU,alpha)
            for j,rr in enumerate(rows):models.append(('E'+str(e)+'_LOG_ETA'+str(ev)+'_'+str(j),f,h,pv,qv,rr,Bdim,False,[]))
    f=a**3;h=a*a
    for j,(psi,Bdim)in enumerate([(s.Integer(0),4),(xx,6),(xx*xx,7),(xx*xx+yy,8)]):models.append(('TRIPLE_PHI_ZERO_'+str(j),f,h,dlog(a),dlog(a),differential(psi),Bdim,False,[]))
    for j,(rr,Bdim)in enumerate([(s.zeros(3,1),7),(xx*dx,8),(yy*dx-dy,9)]):models.append(('TRIPLE_PHI_LINEAR_'+str(j),f,h,dlog(a)+dx,dlog(a)+2*dx,rr,Bdim,False,[]))
    for j,(rr,Bdim)in enumerate([(s.zeros(3,1),8),(dx,9)]):models.append(('TRIPLE_PHI_SQUARE_'+str(j),f,h,dlog(a)+2*xx*dx,dlog(a)+4*xx*dx,rr,Bdim,False,[]))
    assert len(models)==90 and sum(row[-2]for row in models)==15
    records=[]
    checkpoint=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_CHECKPOINT_20261108.json'
    if checkpoint.exists():raise FileExistsError(checkpoint)
    for name,f,h,p,q,r,Bdim,variable,excluded in models:
        beta=[[p,s.zeros(3,1)],[r,q]]
        for i,j in combinations(range(3),2):
            assert s.cancel(s.diff(p[j],base[i])-s.diff(p[i],base[j]))==0
            assert s.cancel(s.diff(q[j],base[i])-s.diff(q[i],base[j]))==0
            assert s.cancel(s.diff(r[j],base[i])-s.diff(r[i],base[j])-(q[i]-p[i])*r[j]+(q[j]-p[j])*r[i])==0
        Pbase=s.zeros(5)
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:Pbase[i,j]=s.diff(f,base[k]);Pbase[j,i]=-Pbase[i,j]
        omegas=[]
        for i in range(2):
            for j in range(2):
                om=(3*f*beta[i][j]-(differential(f)if i==j else s.zeros(3,1))).applyfunc(s.cancel)
                assert s.expand(sum(base[k]*om[k]for k in range(3)))==0
                assert all(s.denom(v)==1 for v in om)
                V=s.Matrix([s.diff(om[2],b)-s.diff(om[1],c),s.diff(om[0],c)-s.diff(om[2],a),s.diff(om[1],a)-s.diff(om[0],b)])/3
                assert s.expand(sum(s.diff(V[k],base[k])for k in range(3)))==0
                inverse=s.Matrix([c*V[1]-b*V[2],a*V[2]-c*V[0],b*V[0]-a*V[1]])
                assert all(s.expand(v)==0 for v in inverse-om)
                for k in range(3):Pbase[k,3+i]-=V[k]*normals[j];Pbase[3+i,k]+=V[k]*normals[j]
                omegas.append({'i':i,'j':j,'beta':list(map(str,beta[i][j])),'horizontal_polynomial_form':list(map(str,om)),'inverse_tracefree_field':list(map(str,V))})
        target=s.expand(h*z**3);target_vector=s.Matrix([s.Poly(target,*x).coeff_monomial(mon)for mon in monG]);idx=next(i for i,v in enumerate(target_vector)if v!=0)
        Ann=s.zeros(23,24)
        for row,i in enumerate(j for j in range(24)if j!=idx):Ann[row,i]=1;Ann[row,idx]=-target_vector[i]/target_vector[idx]
        assert Ann*target_vector==s.zeros(23,1)
        def analyze(alpha=None,domain=FIELD):
            sub=lambda v:v if alpha is None else v.subs(TAU,alpha)
            Pk=sub(Pbase);Pk[3,4]=kk[0]*z*z+kk[1]*z*w+kk[2]*w*w;Pk[4,3]=-Pk[3,4]
            Dk=divisor(Pk);poly=s.Poly(Dk,*x);rows=s.Matrix([poly.coeff_monomial(mon)for mon in monK]);KM=rows.jacobian(kk);const=rows.subs(dict.fromkeys(kk,0))
            kc,km,_=certify(KM,domain);ka,am,_=certify(KM.row_join(-const),domain);den=[km,am]
            record={'full30_global_K_matrix':serial_matrix(KM),'full30_global_K_constant':serial_matrix(const),'K_certificate':kc,'K_augmented_certificate':ka}
            if ka['rank']>kc['rank']:record.update(global_K_inconsistent=True,whole_m3_source_empty=True,G_closed=True);return record,den
            sol=list(s.linsolve((KM,-const),kk));assert len(sol)==1 and kc['rank']==3
            P0=sub(Pbase);P0[3,4]=s.expand(sum(v*m for v,m in zip(sol[0],[z*z,z*w,w*w])));P0[4,3]=-P0[3,4]
            assert all(domain.from_sympy(v)==domain.zero for v in s.Poly(divisor(P0),*x).coeffs())
            Jac0=[s.expand(v)for v in H.all350(H.projective(H.bilinear(P0,P0,x),x),x)];bad=[(i,v)for i,v in enumerate(Jac0)if domain.from_sympy(v)!=domain.zero]
            record.update(unique_global_K=str(P0[3,4]),full150_P0={str((i,j)):str(P0[i,j])for i,j in H.PAIRS})
            if bad:
                record.update(P0_nonzero_Jacobi=[(i,str(v))for i,v in bad],whole_m3_source_empty=True,G_closed=True);den.append(bad[0][1]);return record,den
            cache=old_cache.get(p0_key(P0,sub(target)))
            if cache is not None:
                original,reference=cache
                copied=json.loads(json.dumps({key:original[key]for key in reuse_fields}))
                copied['reused_exact_P0_operator_certificate']=reference
                record.update(copied)
                record['whole_projective_image_upper']=Bdim+5-(1 if alpha is not None and variable else 0)+record['weight1_certificate']['nullity']+record['weight2_certificate']['nullity']
                record['G_closed']=record['whole_irreducible_source_empty']or record['whole_projective_image_upper']<=24
                den.extend(parse(record[key]['minor_determinant'])for key in['weight1_certificate','weight2_certificate'])
                den.extend(s.denom(s.cancel(v))for key in['weight1_certificate','weight2_certificate']for v in matrix(record[key]['complete_kernel']))
                den.extend(s.denom(s.cancel(v))for v in sol[0])
                return record,den
            D1=operator(P0,P1,v1,x);D2=operator(P0,P2,v2,x);full=divisor(P0+P1+P2);pp=s.Poly(full,*x)
            D3=s.Matrix([pp.coeff_monomial(mon)for mon in monG]);L=D3.jacobian(v1);assert all(s.expand(v)==0 for v in D3-L*s.Matrix(v1))
            M1=D1.col_join(Ann*L)
            c1,m1,K1=certify(M1,domain);c2,m2,K2=certify(D2,domain);image=(L*K1).applyfunc(s.cancel)
            qbad=s.cancel(3*f*sum(sub(q[i])*v1[-3+i]for i in range(3)));qbpoly=s.Poly(qbad,*base)
            QM=s.Matrix([qbpoly.coeff_monomial(mon)for mon in H.monomials(3,2)]).jacobian(v1)
            assert all(domain.from_sympy(v)==domain.zero for v in QM*K1)
            assert all(domain.from_sympy(v)==domain.zero for v in Ann*image)
            leading=[domain.from_sympy(image[idx,j]/target_vector[idx])for j in range(image.cols)];nonzero=[j for j,v in enumerate(leading)if v!=domain.zero]
            offset=Bdim+5-(1 if alpha is not None and variable else 0);upper=offset+K1.cols+K2.cols
            badzero=all(domain.from_sympy(v)==domain.zero for v in K1[-3:,:])
            assert all(domain.from_sympy(v)==domain.zero for i in range(5)for v in s.Poly(P0[i,3].subs(z,0),*x).coeffs())
            record.update(P0_full350_Jacobi_zero_and_actualD_zero=True,full350_weight1_operator=serial_matrix(D1),full350_weight2_operator=serial_matrix(D2),
                full24by30_actual_D3_map=serial_matrix(L),full_weight1_plus_line_matrix=serial_matrix(M1),
                weight1_certificate=c1,weight2_certificate=c2,whole_D3_kernel_image=serial_matrix(image),
                all_three_bad_normal_kernel_rows_zero=badzero,proper_D_factor_by_bad_normal_rows='z'if badzero else None,
                full6by30_constant_root_character_matrix=serial_matrix(QM),all_full_kernel_q_Cbad_zero=True,
                nonzero_D3_line_certificate=None if not nonzero else {'kernel_column':nonzero[0],'coefficient':str(domain.to_sympy(leading[nonzero[0]]))},
                D3_image_rank=int(bool(nonzero)),whole_m3_source_empty=not nonzero,whole_irreducible_source_empty=not nonzero or badzero,whole_projective_image_upper=upper,G_closed=not nonzero or badzero or upper<=24)
            den.extend([m1,m2]);den.extend(s.denom(s.cancel(v))for v in list(K1)+list(K2)+list(sol[0]))
            return record,den
        result,den=analyze(domain=FIELD if variable else s.QQ);result.update(type=name,whole_alpha_line=variable,full_connection_inverse_forms=omegas,
            whole_fixed_plane_symbol_dimension_upper=Bdim,exact_nonzero_leading_target=str(target),full_line_annihilator=serial_matrix(Ann))
        if variable:
            exception=s.Poly(s.prod(TAU-v for v in excluded),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            walls=[]
            for factor,power in s.factor_list(exception.as_expr())[1]:
                poly=s.Poly(factor,TAU);alpha=-poly.nth(0)/poly.nth(1)if poly.degree()==1 else s.CRootOf(poly.as_expr(),0)
                if alpha in excluded:
                    walls.append({'alpha':str(alpha),'whole_m3_source_empty':False,'full_special_source_elsewhere':True,'G_closed':True});continue
                domain=s.QQ if poly.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(factor),all_complex_embeddings_covered=True);walls.append(wall)
            result.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,
                all_strata_G_closed=result['G_closed']and all(v['G_closed']for v in walls))
        else:result['all_strata_G_closed']=result['G_closed']
        records.append(result)
        if len(records)%10==0:checkpoint.write_text(json.dumps({'task_id':'G','complete':False,'completed_source_count':len(records),'target_source_count':90,'cases':records,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text()},indent=2)+'\n')
        print(name,'reused',bool(result.get('reused_exact_P0_operator_certificate')),'empty',result.get('whole_irreducible_source_empty',result['whole_m3_source_empty']),'k',result.get('weight1_certificate',{}).get('nullity'),
            result.get('weight2_certificate',{}).get('nullity'),'upper',result.get('whole_projective_image_upper'),
            'walls',[(v['alpha'],v['whole_m3_source_empty'],v.get('whole_projective_image_upper'))for v in result.get('all_exceptional_sources',[])],flush=True)
    proofs=[ROOT/'research/classification'/name for name in['PLANE_M3_FREE_CONNECTION_SIGN_ERRATUM.md','PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_CONNECTIONS_V2.md','PLANE_NONREDUCED_M3_TRIPLE_IRREGULAR_C_EXCLUSION.md','PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md','PLANE_M3_CONSTANT_ROOT_OBSTRUCTION.md']]
    rx,ry,et,ga,de=s.symbols('chart_x chart_y eta gamma delta');coef=s.symbols('A B C D E F G H');A,B,C,D,E,F,G,Hh=coef
    PP=A+B*rx+C*ry+D*rx*ry+E*ry**2;QQ=F+G*rx+Hh*ry-D*rx**2-E*rx*ry
    covariant=s.expand(rx*(s.diff(QQ,rx)-s.diff(PP,ry))-(et+1)*QQ-ga*rx*QQ+de*rx*PP)
    kernel_receipt={'all_eight_coefficients':list(map(str,coef)),'P':str(PP),'Q':str(QQ),'full_covariant_identity':str(covariant),
        'log_full_coefficient_equations':list(map(str,s.Poly(covariant.subs({ga:0,de:0}),rx,ry).coeffs())),
        'gamma_one_full_coefficient_equations':list(map(str,s.Poly(covariant.subs({ga:1,de:0}),rx,ry).coeffs()))}
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all90_G_closed':all(v['all_strata_G_closed']for v in records),'cases':records,'corrected_covariant_r8':kernel_receipt,'whole_old90_coverage_withdrawn':True,'exact_main_analyses_reused':sum('reused_exact_P0_operator_certificate'in v for v in records),
        'scope':'Correct minus curvature whole generictriple firstconnections after fullsource obstructions:15 wholealpha+75 fixed. Old90 coverage is explicitly withdrawn. Reuse only exact full rawP0+target equality and unchanged helper FunctionDefs, retaining copied explicit fullmatrix/kernel/minor certificates. All new beta/globalK/P0/full30/9/D3/badrows and every wholeline wall checked; G only where explicit true.',
        'all24_D3_monomials':monG,'all30_K_monomials':monK,'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,archive]+proofs],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')
    if checkpoint.exists():checkpoint.unlink()


def nonreduced_m3_triple_corrected_ore_main():
    """New correct-minus r=dy/x fixed source: complete affine Q/W, original Jacobi/Pf structure."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_ORE_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('one corrected-r affine Ore-source finite budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json';packet=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(packet['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==44
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;P1,v1,P2,v2=sources(x);ww=s.symbols('W0:18');mons=H.monomials(4,2)
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    def matrix(raw):return s.Matrix([[parse(v)for v in row]for row in raw])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1];free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list()
        kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());mi=dm.extract(list(rows),list(cols)).det();assert mi!=domain.zero or not cols
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(mi))),'complete_kernel':serial_matrix(K)},domain.to_sympy(mi),K
    selected=[dict(next(v for v in packet['cases']if v['type']=='E2_LOG_ETA-1_4'))]
    selected[0]['whole_fixed_plane_symbol_dimension_upper']=8
    def analyze(old,alpha=None,domain=FIELD):
        sub=lambda v:v if alpha is None else v.subs(TAU,alpha);vertex=2;indices=[0,1,3,4];y=s.Matrix([a,b,z,w]);base0,base1=y[:2];L0=s.zeros(5,1)
        M1=sub(matrix(old['full_weight1_plus_line_matrix']));c1,_,K=certify(M1,domain);tt=s.symbols('p0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=sub(parse(old['full150_P0'][str((i,j))]));P0[j,i]=-P0[i,j]
        assert P0[3,4]==0 and all(s.diff(P0[i,j],c)==0 for i,j in H.PAIRS)
        if all(domain.from_sympy(v)==domain.zero for v in K[-3:,:]):
            return {'source_covered':True,'full_original_M1_certificate':c1,'all_bad_normal_rows_zero':True,'Ore_not_needed_because_irreducible_source_empty':True,'proper_D_factor':'z','G_closed':True},[]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(tt))),simultaneous=True)+P2-H.wedge(x,L0)
        L=s.Matrix([-s.diff(P[i,vertex],x[vertex])if i!=vertex else 0 for i in range(5)])
        gauged=(P-H.wedge(x,L)).applyfunc(s.expand);bad=[(str((i,j)),str(s.diff(gauged[i,j],x[vertex])))for i,j in H.PAIRS if s.diff(gauged[i,j],x[vertex])!=0]
        rec={'alpha':None if alpha is None else str(alpha),'vertex_index':vertex,'weight_zero_Euler_gauge':list(map(str,L0)),
            'full_positive_weight_Euler_gauge':list(map(str,L)),'full_original_M1_certificate':c1,'all_vertex_independence':not bad}
        if bad:rec.update(nonzero_vertex_derivatives=bad,source_covered=False,G_closed=False);return rec,[]
        Q=gauged.extract(indices,indices);W=s.Matrix([gauged[vertex,i]for i in indices]);ell=-s.Poly(W[2],z,w).coeff_monomial(z*z)*z-s.Poly(W[3],z,w).coeff_monomial(w*w)*w
        W=(W+ell*y).applyfunc(s.expand);variables=list(tt)+v2;Q0=Q.subs(dict.fromkeys(variables,0));W0=W.subs(dict.fromkeys(variables,0))
        qcoords=s.Matrix([s.Poly(Q[i,j]-Q0[i,j],*y).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in mons]);J=qcoords.jacobian(variables)
        assert all(s.expand(v)==0 for v in qcoords-J*s.Matrix(variables));qc,qminor,_=certify(J,domain);cols=qc['minor_columns'];uu=s.symbols('Q0:'+str(len(cols)));Qcan=Q0.copy()
        for k,(i,j)in enumerate(combinations(range(4),2)):
            Qcan[i,j]+=sum(J[10*k+r,jj]*s.prod(v**e for v,e in zip(y,mon))*q for jj,q in zip(cols,uu)for r,mon in enumerate(mons));Qcan[j,i]=-Qcan[i,j]
        Wcan=W0+s.Matrix([z*(ww[0]*base0+ww[1]*base1)+w*(ww[2]*base0+ww[3]*base1)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*base0+ww[8]*base1)+w*(ww[9]*base0+ww[10]*base1)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w]);wc=[];wd=W-W0
        for i in range(2):wc.extend([s.Poly(wd[i],*y).coeff_monomial(mon)for mon in[base0*z,base1*z,base0*w,base1*w,z*z,z*w,w*w]])
        wc.extend([s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)])
        assert all(s.expand(v)==0 for v in W-Wcan.subs(dict(zip(ww,wc)),simultaneous=True));Wmap=s.Matrix(wc).jacobian(variables)
        joint=J.extract(qc['minor_rows'],list(range(J.cols))).col_join(Wmap);jc,jminor,JK=certify(joint,domain)
        vv=s.symbols('V0:4');V=s.Matrix(vv);cv=[]
        for i,j,k in combinations(range(4),3):cv.append(V[i]*Qcan[j,k]-V[j]*Qcan[i,k]+V[k]*Qcan[i,j])
        cv.extend(sum(V[k]*s.diff(Qcan[i,j],y[k])for k in range(4))for i,j in combinations(range(4),2))
        cvcoeff=[val for expr in cv for val in s.Poly(s.expand(expr),*y,*uu).coeffs()if val!=0]
        Cvm=s.Matrix(cvcoeff).jacobian(vv)if cvcoeff else s.zeros(0,4);vc,vminor,VK=certify(Cvm,domain)
        strict={(i,j,k):s.expand(sum(Qcan[i,l]*s.diff(Qcan[j,k],y[l])+Qcan[j,l]*s.diff(Qcan[k,i],y[l])+Qcan[k,l]*s.diff(Qcan[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        Pi=s.zeros(5)
        for i,j in combinations(range(4),2):Pi[indices[i],indices[j]]=Qcan[i,j];Pi[indices[j],indices[i]]=-Qcan[i,j]
        for i in range(4):Pi[vertex,indices[i]]=Wcan[i];Pi[indices[i],vertex]=-Wcan[i]
        D=divisor(Pi);dp=s.Poly(D,*x);D3=s.expand(sum(v*s.prod(t**e for t,e in zip(x,mon))for mon,v in dp.terms()if mon[3]+mon[4]==3))
        Ann=matrix(old['full_line_annihilator']);leading=s.Matrix([dp.coeff_monomial(mon)for mon in packet['all24_D3_monomials']]);line=(Ann*leading).applyfunc(s.expand)
        constant=[]
        for v in range(4):
            if all(s.expand(Qcan[i,j])==0 for i,j in combinations(range(4),2)if v not in[i,j])and all(s.diff(Qcan[i,j],y[v])==0 for i,j in combinations(range(4),2)):constant.append(str(y[v]))
        rec.update(source_covered=True,remaining_Euler_c_component=str(ell),all_Q_coordinates=list(map(str,y)),full_Q0={str((i,j)):str(Q0[i,j])for i,j in combinations(range(4),2)},
            full_W0=list(map(str,W0)),full60_Q_source_matrix=serial_matrix(J),full_W18_source_matrix=serial_matrix(Wmap),Q_space_certificate=qc,
            full_affine_Q={str((i,j)):str(Qcan[i,j])for i,j in combinations(range(4),2)},full_strict_Q_Jacobi={str(k):str(v)for k,v in strict.items()},
            full_W18_cover=True,full_actual_D_polynomial=str(D),full_actual_D3_polynomial=str(s.factor(D3)),
            exact_nonzero_leading_target=old['exact_nonzero_leading_target'],full23_D3_line_equations=list(map(str,line)),
            entire_Q_constant_factor_directions=constant,whole_QW_source_joint_matrix=serial_matrix(joint),joint_QW_certificate=jc,
            joint_QW_bijection=joint.rows==joint.cols==jc['rank'],full_common_constant_v_matrix=serial_matrix(Cvm),common_constant_v_certificate=vc,
            whole_Q_common_constant_v=None if VK.cols==0 else list(map(str,VK[:,0])),G_closed=VK.cols>0)
        den=[qminor,jminor,vminor];den.extend(s.denom(s.cancel(v))for v in list(K)+list(J)+list(JK)+list(VK));return rec,den
    results=[]
    for old in selected:
        variable=old['whole_alpha_line'];rec,den=analyze(old,domain=FIELD if variable else s.QQ);rec.update(type=old['type'],whole_alpha_line=variable,whole_fixed_plane_symbol_dimension_upper=old['whole_fixed_plane_symbol_dimension_upper']);walls=[]
        if variable:
            exception=s.Poly(parse(old['all_exception_polynomial']),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            for pol,power in s.factor_list(exception.as_expr())[1]:
                pp=s.Poly(pol,TAU);alpha=-pp.nth(0)/pp.nth(1)if pp.degree()==1 else s.CRootOf(pp.as_expr(),0)
                known=next((v for v in old['all_exceptional_sources']if parse(v['alpha'])==alpha),None)
                if known is not None and known.get('whole_m3_source_empty',False):walls.append({'alpha':str(alpha),'whole_m3_source_empty_by_original_certificate':True,'G_closed':True});continue
                domain=s.QQ if pp.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(old,alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(pol),all_complex_embeddings_covered=True);walls.append(wall)
            rec.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,all_strata_G_closed=rec['G_closed']and all(v['G_closed']for v in walls))
        else:rec['all_strata_G_closed']=rec['G_closed']
        results.append(rec);print(old['type'],'covered',rec['source_covered'],'Qdim',rec.get('Q_space_certificate',{}).get('rank'),'factor',rec.get('entire_Q_constant_factor_directions'),
            'walls',[(v['alpha'],v.get('Q_space_certificate',{}).get('rank'))for v in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_ORE_ENTRANCE.md';partial=ROOT/'research/classification/PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_CONNECTIONS_V2.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_triple_corrected_ore_preappend_20261108.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'new_fixed_G_closed':all(v['all_strata_G_closed']for v in results),'cases':results,
        'scope':'Single genuinely new fixed correct-minus connection E2 eta=-1 r=dy/x. Entire original M1 kernel15+fullP2nine, affine Q0 retained, fullEuler/QWjointsource/strictJac/D3/commonconstant-v. Old coverage explicitly withdrawn; corrected whole proof separately. Failed criterion remains unresolved.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,partial,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def nonreduced_m3_triple_vertex_b_main():
    """New correct-minus r=dy/x fixed source: alternate actual vertex-b Euler gauge, complete Q/W source."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_VERTEX_B_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('one corrected-r affine Ore-source finite budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json';packet=json.loads(previous.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(packet['producer_bytes']);after=functions(Path(__file__).read_text());assert len(before)==44
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;P1,v1,P2,v2=sources(x);ww=s.symbols('W0:18');mons=H.monomials(4,2)
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    def matrix(raw):return s.Matrix([[parse(v)for v in row]for row in raw])
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(v)for v in row]for row in M.tolist()],M.shape,domain);rr,cols=dm.rref();rows=dm.transpose().rref()[1];free=[i for i in range(M.cols)if i not in cols];raw=rr.to_list()
        kernel=[[domain.zero for j in free]for i in range(M.cols)]
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(cols):kernel[pivot][j]=-raw[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list());mi=dm.extract(list(rows),list(cols)).det();assert mi!=domain.zero or not cols
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(cols),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(cols),'minor_determinant':str(s.factor(domain.to_sympy(mi))),'complete_kernel':serial_matrix(K)},domain.to_sympy(mi),K
    selected=[dict(next(v for v in packet['cases']if v['type']=='E2_LOG_ETA-1_4'))]
    selected[0]['whole_fixed_plane_symbol_dimension_upper']=8
    def analyze(old,alpha=None,domain=FIELD):
        sub=lambda v:v if alpha is None else v.subs(TAU,alpha);vertex=1;indices=[0,2,3,4];y=s.Matrix([a,c,z,w]);base0,base1=y[:2];L0=s.Matrix([0,0,2*a,0,0])
        M1=sub(matrix(old['full_weight1_plus_line_matrix']));c1,_,K=certify(M1,domain);tt=s.symbols('p0:'+str(K.cols));P0=s.zeros(5)
        for i,j in H.PAIRS:P0[i,j]=sub(parse(old['full150_P0'][str((i,j))]));P0[j,i]=-P0[i,j]
        assert all(s.diff((P0-H.wedge(x,L0))[i,j],b)==0 for i,j in H.PAIRS)
        if all(domain.from_sympy(v)==domain.zero for v in K[-3:,:]):
            return {'source_covered':True,'full_original_M1_certificate':c1,'all_bad_normal_rows_zero':True,'Ore_not_needed_because_irreducible_source_empty':True,'proper_D_factor':'z','G_closed':True},[]
        P=P0+P1.subs(dict(zip(v1,K*s.Matrix(tt))),simultaneous=True)+P2-H.wedge(x,L0)
        L=s.Matrix([-s.diff(P[i,vertex],x[vertex])if i!=vertex else 0 for i in range(5)])
        gauged=(P-H.wedge(x,L)).applyfunc(s.expand);bad=[(str((i,j)),str(s.diff(gauged[i,j],x[vertex])))for i,j in H.PAIRS if s.diff(gauged[i,j],x[vertex])!=0]
        rec={'alpha':None if alpha is None else str(alpha),'vertex_index':vertex,'weight_zero_Euler_gauge':list(map(str,L0)),
            'full_positive_weight_Euler_gauge':list(map(str,L)),'full_original_M1_certificate':c1,'all_vertex_independence':not bad}
        if bad:rec.update(nonzero_vertex_derivatives=bad,source_covered=False,G_closed=False);return rec,[]
        Q=gauged.extract(indices,indices);W=s.Matrix([gauged[vertex,i]for i in indices]);ell=-s.Poly(W[2],z,w).coeff_monomial(z*z)*z-s.Poly(W[3],z,w).coeff_monomial(w*w)*w
        W=(W+ell*y).applyfunc(s.expand);variables=list(tt)+v2;Q0=Q.subs(dict.fromkeys(variables,0));W0=W.subs(dict.fromkeys(variables,0))
        qcoords=s.Matrix([s.Poly(Q[i,j]-Q0[i,j],*y).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in mons]);J=qcoords.jacobian(variables)
        assert all(s.expand(v)==0 for v in qcoords-J*s.Matrix(variables));qc,qminor,_=certify(J,domain);cols=qc['minor_columns'];uu=s.symbols('Q0:'+str(len(cols)));Qcan=Q0.copy()
        for k,(i,j)in enumerate(combinations(range(4),2)):
            Qcan[i,j]+=sum(J[10*k+r,jj]*s.prod(v**e for v,e in zip(y,mon))*q for jj,q in zip(cols,uu)for r,mon in enumerate(mons));Qcan[j,i]=-Qcan[i,j]
        Wcan=W0+s.Matrix([z*(ww[0]*base0+ww[1]*base1)+w*(ww[2]*base0+ww[3]*base1)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*base0+ww[8]*base1)+w*(ww[9]*base0+ww[10]*base1)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w]);wc=[];wd=W-W0
        for i in range(2):wc.extend([s.Poly(wd[i],*y).coeff_monomial(mon)for mon in[base0*z,base1*z,base0*w,base1*w,z*z,z*w,w*w]])
        wc.extend([s.Poly(wd[2],z,w).coeff_monomial(z*w),s.Poly(wd[2],z,w).coeff_monomial(w*w),s.Poly(wd[3],z,w).coeff_monomial(z*z),s.Poly(wd[3],z,w).coeff_monomial(z*w)])
        assert all(s.expand(v)==0 for v in W-Wcan.subs(dict(zip(ww,wc)),simultaneous=True));Wmap=s.Matrix(wc).jacobian(variables)
        joint=J.extract(qc['minor_rows'],list(range(J.cols))).col_join(Wmap);jc,jminor,JK=certify(joint,domain)
        vv=s.symbols('V0:4');V=s.Matrix(vv);cv=[]
        for i,j,k in combinations(range(4),3):cv.append(V[i]*Qcan[j,k]-V[j]*Qcan[i,k]+V[k]*Qcan[i,j])
        cv.extend(sum(V[k]*s.diff(Qcan[i,j],y[k])for k in range(4))for i,j in combinations(range(4),2))
        cvcoeff=[val for expr in cv for val in s.Poly(s.expand(expr),*y,*uu).coeffs()if val!=0]
        Cvm=s.Matrix(cvcoeff).jacobian(vv)if cvcoeff else s.zeros(0,4);vc,vminor,VK=certify(Cvm,domain)
        strict={(i,j,k):s.expand(sum(Qcan[i,l]*s.diff(Qcan[j,k],y[l])+Qcan[j,l]*s.diff(Qcan[k,i],y[l])+Qcan[k,l]*s.diff(Qcan[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        Pi=s.zeros(5)
        for i,j in combinations(range(4),2):Pi[indices[i],indices[j]]=Qcan[i,j];Pi[indices[j],indices[i]]=-Qcan[i,j]
        for i in range(4):Pi[vertex,indices[i]]=Wcan[i];Pi[indices[i],vertex]=-Wcan[i]
        D=divisor(Pi);dp=s.Poly(D,*x);D3=s.expand(sum(v*s.prod(t**e for t,e in zip(x,mon))for mon,v in dp.terms()if mon[3]+mon[4]==3))
        Ann=matrix(old['full_line_annihilator']);leading=s.Matrix([dp.coeff_monomial(mon)for mon in packet['all24_D3_monomials']]);line=(Ann*leading).applyfunc(s.expand)
        constant=[]
        for v in range(4):
            if all(s.expand(Qcan[i,j])==0 for i,j in combinations(range(4),2)if v not in[i,j])and all(s.diff(Qcan[i,j],y[v])==0 for i,j in combinations(range(4),2)):constant.append(str(y[v]))
        rec.update(source_covered=True,remaining_Euler_c_component=str(ell),all_Q_coordinates=list(map(str,y)),full_Q0={str((i,j)):str(Q0[i,j])for i,j in combinations(range(4),2)},
            full_W0=list(map(str,W0)),full60_Q_source_matrix=serial_matrix(J),full_W18_source_matrix=serial_matrix(Wmap),Q_space_certificate=qc,
            full_affine_Q={str((i,j)):str(Qcan[i,j])for i,j in combinations(range(4),2)},full_strict_Q_Jacobi={str(k):str(v)for k,v in strict.items()},
            full_W18_cover=True,full_actual_D_polynomial=str(D),full_actual_D3_polynomial=str(s.factor(D3)),
            exact_nonzero_leading_target=old['exact_nonzero_leading_target'],full23_D3_line_equations=list(map(str,line)),
            entire_Q_constant_factor_directions=constant,whole_QW_source_joint_matrix=serial_matrix(joint),joint_QW_certificate=jc,
            joint_QW_bijection=joint.rows==joint.cols==jc['rank'],full_common_constant_v_matrix=serial_matrix(Cvm),common_constant_v_certificate=vc,
            whole_Q_common_constant_v=None if VK.cols==0 else list(map(str,VK[:,0])),G_closed=VK.cols>0)
        den=[qminor,jminor,vminor];den.extend(s.denom(s.cancel(v))for v in list(K)+list(J)+list(JK)+list(VK));return rec,den
    results=[]
    for old in selected:
        variable=old['whole_alpha_line'];rec,den=analyze(old,domain=FIELD if variable else s.QQ);rec.update(type=old['type'],whole_alpha_line=variable,whole_fixed_plane_symbol_dimension_upper=old['whole_fixed_plane_symbol_dimension_upper']);walls=[]
        if variable:
            exception=s.Poly(parse(old['all_exception_polynomial']),TAU)
            for val in den:
                nu,de=s.fraction(s.cancel(val));exception=s.lcm(exception,s.Poly(nu,TAU));exception=s.lcm(exception,s.Poly(de,TAU))
            for pol,power in s.factor_list(exception.as_expr())[1]:
                pp=s.Poly(pol,TAU);alpha=-pp.nth(0)/pp.nth(1)if pp.degree()==1 else s.CRootOf(pp.as_expr(),0)
                known=next((v for v in old['all_exceptional_sources']if parse(v['alpha'])==alpha),None)
                if known is not None and known.get('whole_m3_source_empty',False):walls.append({'alpha':str(alpha),'whole_m3_source_empty_by_original_certificate':True,'G_closed':True});continue
                domain=s.QQ if pp.degree()==1 else s.QQ.algebraic_field(alpha);wall,_=analyze(old,alpha,domain)
                wall.update(alpha=str(alpha),minimal_polynomial=str(pol),all_complex_embeddings_covered=True);walls.append(wall)
            rec.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_exceptional_sources=walls,all_strata_G_closed=rec['G_closed']and all(v['G_closed']for v in walls))
        else:rec['all_strata_G_closed']=rec['G_closed']
        results.append(rec);print(old['type'],'covered',rec['source_covered'],'Qdim',rec.get('Q_space_certificate',{}).get('rank'),'factor',rec.get('entire_Q_constant_factor_directions'),
            'walls',[(v['alpha'],v.get('Q_space_certificate',{}).get('rank'))for v in walls],flush=True)
    proof=ROOT/'research/classification/PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_ORE_ENTRANCE.md';partial=ROOT/'research/classification/PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_CONNECTIONS_V2.md';archive=ROOT/'research/classification/inputs/plane_nonreduced_m3_triple_vertex_b_preappend_20261108.json'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'new_fixed_G_closed':all(v['all_strata_G_closed']for v in results),'cases':results,
        'scope':'Alternate geometric vertex-b of genuinely new fixed correct-minus connection E2 eta=-1 r=dy/x; zero-weight Euler gauge 2a partial_c is retained. Entire original M1 kernel15+fullP2nine, affine Q0 retained, fullEuler/QWjointsource/strictJac/D3/commonconstant-v. Old coverage explicitly withdrawn; corrected whole proof separately. Failed criterion remains unresolved.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,partial,archive]],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def m3_triple_valid_ore_modules_main():
    """Six valid literal triple Ore sources: full strict support reduced to18 affine modules, all walls retained."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_M3_TRIPLE_VALID_ORE_MODULES_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('valid18triple affineOre finite budget240s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_ORE_20261108.json';old=json.loads(previous.read_text())
    archive=ROOT/'research/classification/inputs/plane_m3_triple_valid_modules_preappend_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==46
    assert all(after[name]==value for name,value in before.items())
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;y=s.Matrix([a,b,z,w]);indices=[0,1,3,4];ww=s.symbols('W0:18')
    def makeQ(ab,az,aw,bz,bw):
        Q=s.zeros(4)
        for i,j,v in [(0,1,ab),(0,2,az),(0,3,aw),(1,2,bz),(1,3,bw)]:Q[i,j]=v;Q[j,i]=-v
        return Q
    def Wbasis(W0):
        return W0+s.Matrix([z*(ww[0]*a+ww[1]*b)+w*(ww[2]*a+ww[3]*b)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*a+ww[8]*b)+w*(ww[9]*a+ww[10]*b)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w])
    def pad(Q,W):
        P=s.zeros(5)
        for i,j in combinations(range(4),2):P[indices[i],indices[j]]=Q[i,j];P[indices[j],indices[i]]=-Q[i,j]
        for i in range(4):P[2,indices[i]]=W[i];P[indices[i],2]=-W[i]
        return P
    def divisor(P):
        return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    def certify(M,domain):
        dm=DomainMatrix([[domain.from_sympy(t)for t in row]for row in M.tolist()],M.shape,domain);rr,columns=dm.rref();rows=dm.transpose().rref()[1]
        free=[i for i in range(M.cols)if i not in columns];kernel=[[domain.zero for j in free]for i in range(M.cols)];reduced=rr.to_list()
        for j,k in enumerate(free):
            kernel[k][j]=domain.one
            for i,pivot in enumerate(columns):kernel[pivot][j]=-reduced[i][k]
        kd=DomainMatrix(kernel,(M.cols,len(free)),domain);assert not any(any(row)for row in(dm*kd).to_list())
        minor=dm.extract(list(rows),list(columns)).det();assert minor!=domain.zero or not columns
        K=s.Matrix(M.cols,len(free),lambda i,j:domain.to_sympy(kernel[i][j]))
        return {'rank':len(columns),'nullity':len(free),'minor_rows':list(rows),'minor_columns':list(columns),
            'minor_determinant':str(s.factor(domain.to_sympy(minor))),'complete_kernel':serial_matrix(K)},domain.to_sympy(minor),K
    def analyze(M,C,Dc,Dl,domain,offset):
        cm,minor,K=certify(M,domain);ca,aminor,KA=certify(M.row_join(-C),domain)
        record={'full350_affine_matrix':serial_matrix(M),'full350_affine_constant':serial_matrix(C),'linear_certificate':cm,
            'augmented_certificate':ca,'image_dimension_offset':offset};denominators=[]
        if cm['rank']<ca['rank']:
            record.update(status='inconsistent',G_closed=True,G_empty=True);return record,[minor,aminor]
        dm=DomainMatrix([[domain.from_sympy(t)for t in row]for row in M.row_join(-C).tolist()],(M.rows,19),domain);rr,pivots=dm.rref();raw=rr.to_list()
        part=s.zeros(18,1)
        for i,p in enumerate(pivots):part[p]=domain.to_sympy(raw[i][18])
        assert all(domain.from_sympy(t)==domain.zero for t in M*part+C)
        D0=s.Poly(s.expand(Dc+sum(Dl[i]*part[i]for i in range(18))),*x,domain=domain)
        directions=[s.Poly(s.expand(sum(Dl[i]*K[i,j]for i in range(18))),*x,domain=domain)for j in range(K.cols)]
        common=D0
        for poly in directions:
            if not poly.is_zero:common=s.polys.polytools.gcd(common,poly).monic()
        factor=None
        if common.is_zero:zero=True
        else:
            zero=False
            if 0<common.total_degree()<5:factor=common.monic()
            elif common.total_degree()==5:
                for poly,power in common.factor_list()[1]:
                    if 0<poly.total_degree()<5:factor=poly.monic();break
        upper=offset+K.cols
        record.update(status='whole_affine_module',particular_solution=list(map(str,part)),W_kernel_dimension=K.cols,
            full_actual_D_particular=str(D0.as_expr()),full_actual_D_kernel_coefficients=[str(poly.as_expr())for poly in directions],
            common_D_polynomial_gcd=str(common.as_expr()),proper_D_factor=None if factor is None else str(factor.as_expr()),
            whole_actual_D_zero=zero,whole_projective_image_upper=upper,G_empty=zero or factor is not None,G_closed=zero or factor is not None or upper<=24)
        if domain==FIELD:
            denominators.extend(s.denom(s.cancel(t))for t in list(part)+list(K))
            if factor is not None:
                for poly in [factor]+[poly.exquo(factor)for poly in [D0]+directions]:denominators.extend(s.denom(s.cancel(t))for t in poly.coeffs())
        return record,[minor,aminor]+denominators
    models=[];t=TAU;ell=(6*t-2)/(3*t+2)
    Wdx=s.Matrix([a*a,-2*a*b,3*t*a*z,(6*t-1)*a*w-3*b*z])
    Wzero=s.Matrix([a*a,-2*a*b,3*t*a*z,(6*t+2)*a*w])
    Wphi=s.Matrix([0,-3*a*a,0,-3*a*z])
    for eps in[0,1]:
        models.append(('E1DX_'+str(eps),makeQ(-(3*t+2)*a*z+eps*w*w,0,-3*z*z,(6*t+1)*z*z,9*t*z*w),Wdx,9,3+eps,True,[-s.Rational(1,6)]))
        models.append(('E2R_'+str(eps),makeQ((ell-1)*b*z+eps*w*w,ell*z*z,(ell+1)*z*w,0,z*z),Wzero,8,4+eps,True,[s.Rational(1,3),-s.Rational(2,3)]))
        models.append(('PHI0_RHO_'+str(eps),makeQ((t+1)*a*z+eps*w*w,0,t*z*z,-2*z*z,(t-3)*z*w),Wphi,6,5+eps,True,[s.Integer(0)]))
    models.append(('E1DX_EXTRA_ZW',makeQ(-3*a*z+z*w,0,-3*z*z,3*z*z,3*z*w),Wdx.subs(t,s.Rational(1,3)),7,3,False,[]))
    models.append(('E2_CRITICAL_V',makeQ(2*a*w+3*b*z,4*z*z,5*z*w,-4*z*w,z*z-5*w*w),Wzero.subs(t,-s.Rational(5,3)),7,6,False,[]))
    for j,F in enumerate([s.Integer(0),z*w,w*w]):
        models.append(('E2_D0_'+str(j),makeQ(b*z+F,z*z,z*w,0,z*z),Wzero.subs(t,-s.Rational(2,3)),6,3+j,False,[]))
        models.append(('E1_SCALAR_MIX_'+str(j),makeQ(a*z+F,0,z*z,-z*z,z*z-z*w),s.Matrix([a*a,-2*a*b,a*z,a*w-3*(a+b)*z]),8,2+j,False,[]))
        models.append(('E2_SCALAR_LOG_'+str(j),makeQ(2*b*z+F,2*z*z,z*z+2*z*w,0,z*z),s.Matrix([a*a,-2*a*b,-2*a*z,-2*a*w-3*a*z]),7,3+j,False,[]))
    models.append(('PHI0_RHO1_EXTRA_ZW',makeQ(2*a*z+z*w,0,z*z,-2*z*z,-2*z*w),Wphi,6,4,False,[]))
    assert len(models)==18 and sum(v[-2]for v in models)==6
    cases=[]
    for name,Q,W0,Bdim,Qdim,variable,excluded in models:
        strict={(i,j,k):s.cancel(sum(Q[i,l]*s.diff(Q[j,k],y[l])+Q[j,l]*s.diff(Q[k,i],y[l])+Q[k,l]*s.diff(Q[i,j],y[l])for l in range(4)))for i,j,k in combinations(range(4),3)}
        assert all(v==0 for v in strict.values())
        Pi=pad(Q,Wbasis(W0));full=s.Matrix([s.expand(v)for v in H.all350(H.projective(H.bilinear(Pi,Pi,x),x),x)])
        C=full.subs(dict.fromkeys(ww,0));M=full.jacobian(ww);assert all(s.cancel(v)==0 for v in full-C-M*s.Matrix(ww))
        D=divisor(Pi);Dc=D.subs(dict.fromkeys(ww,0));Dl=[s.expand(D).coeff(v)for v in ww];assert s.expand(D-Dc-sum(v*q for v,q in zip(Dl,ww)))==0
        offset=Bdim+Qdim+5;record,den=analyze(M,C,Dc,Dl,FIELD if variable else s.QQ,offset)
        record.update(type=name,whole_alpha_line=variable,full_Q={str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},W0=list(map(str,W0)),
            full_strict_Q_Jacobi_identically_zero=True,whole_fixed_plane_symbol_dimension_upper=Bdim,whole_Q_source_dimension_upper=Qdim)
        if variable:
            exception=s.Poly(1,t)
            for value in den+list(M)+list(C)+Dl+[Dc]:
                if value in den:num,de=s.fraction(s.cancel(value));exception=s.lcm(exception,s.Poly(num,t));exception=s.lcm(exception,s.Poly(de,t))
                else:exception=s.lcm(exception,s.Poly(s.denom(s.cancel(value)),t))
            for value in excluded:exception=s.lcm(exception,s.Poly(t-value,t))
            walls=[]
            for factor,power in s.factor_list(exception.as_expr())[1]:
                p=s.Poly(factor,t);alpha=-p.nth(0)/p.nth(1)if p.degree()==1 else s.CRootOf(p.as_expr(),0)
                if alpha in excluded:
                    walls.append({'alpha':str(alpha),'source_excluded_here':True,'reason':'E1DX leadingzero'if name.startswith('E1DX')else('E2 leadingzero'if alpha==s.Rational(1,3)else'E2_D0fullthreefixed')if name.startswith('E2R')else'PHI0_Bzero_commonconstantfactor' ,'G_closed':True});continue
                dom=s.QQ if p.degree()==1 else s.QQ.algebraic_field(alpha);sub=lambda value:value.subs(t,alpha)
                wall,_=analyze(sub(M),sub(C),sub(Dc),[sub(v)for v in Dl],dom,offset-1)
                wall.update(alpha=str(alpha),minimal_polynomial=str(factor),all_complex_embeddings_covered=True);walls.append(wall)
            record.update(all_exception_polynomial=str(s.factor(exception.as_expr())),all_irreducible_factors=[[str(v),int(e)]for v,e in s.factor_list(exception.as_expr())[1]],
                all_exceptional_fields=walls,all_strata_G_closed=record['G_closed']and all(v['G_closed']for v in walls))
        else:record['all_strata_G_closed']=record['G_closed']
        cases.append(record);print(name,record['status'],'k',record['linear_certificate']['nullity'],'factor',record.get('proper_D_factor'),
            'upper',record.get('whole_projective_image_upper'),'allG',record['all_strata_G_closed'],
            'walls',[(v['alpha'],v.get('W_kernel_dimension'),v.get('proper_D_factor'),v['G_closed'])for v in record.get('all_exceptional_fields',[])],flush=True)
    proofs=[ROOT/'research/classification'/name for name in['PLANE_M3_FREE_CONNECTION_SIGN_ERRATUM.md','PLANE_M3_TRIPLE_VALID_ORE_STRICT_SUPPORT_REDUCTION.md','PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_ORE_ENTRANCE.md']]
    result={'task_id':'G','all_exact_arithmetic_passed':True,'all18_G_closed':all(v['all_strata_G_closed']for v in cases),'cases':cases,
        'scope':'Six valid literal triple Ore source strict supports:18 canonical modules (6whole+12fixed), all full350 affineW18/all actualD coefficients and every parameter wall. This does not restore wrong old coverage or cover new corrected r=dy/x source. Allzero/constant-v branches and special central strata returned by separate pureproof.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,archive]+proofs],
        'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def m3_triple_new_fixed_original_support_main():
    """One entire corrected fixed source: original Jacobi/Pf support and its two residual directions."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_M3_TRIPLE_NEW_FIXED_ORIGINAL_SUPPORT_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('one entire corrected fixed originalsupport finite budget240s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    archive=ROOT/'research/classification/inputs/plane_m3_triple_new_fixed_support_preappend_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==47 and all(after[name]==value for name,value in before.items())
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json';packet=json.loads(previous.read_text());row=next(v for v in packet['cases']if v['type']=='E2_LOG_ETA-1_4')
    x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;indices=[0,1,3,4];y=s.Matrix([a,b,z,w]);P1,v1,P2,v2=sources(x);pp=s.symbols('p0:15');qq=s.symbols('q0:9');variables=list(pp)+list(qq)
    def parse(v):return s.sympify(v,locals={str(t):t for t in[*x,TAU]})
    def matrix(raw):return s.Matrix([[parse(v)for v in rr]for rr in raw])
    K=matrix(row['weight1_certificate']['complete_kernel']);M1=matrix(row['full_weight1_plus_line_matrix']);assert K.shape==(30,15)and M1*K==s.zeros(M1.rows,15)
    rr=row['weight1_certificate']['minor_rows'];cc=row['weight1_certificate']['minor_columns'];assert M1.extract(rr,cc).det()!=0
    P0=s.zeros(5)
    for i,j in H.PAIRS:P0[i,j]=parse(row['full150_P0'][str((i,j))]);P0[j,i]=-P0[i,j]
    P=P0+P1.subs(dict(zip(v1,K*s.Matrix(pp))),simultaneous=True)+P2.subs(dict(zip(v2,qq)),simultaneous=True)
    L=s.Matrix([-s.diff(P[i,2],c)if i!=2 else 0 for i in range(5)]);Pi=(P-H.wedge(x,L)).applyfunc(s.expand)
    W=s.Matrix([Pi[2,i]for i in indices]);R=Pi.extract(indices,indices).diff(c);Q=Pi.extract(indices,indices).subs(c,0)
    assert all(s.diff(v,c)==0 for v in W)and R[0,1]==pp[8]*z and all(R[i,j]==0 for i,j in combinations(range(4),2)if(i,j)!=(0,1))
    strict=H.bilinear(Pi,Pi,x);full=[s.expand(v)for v in H.all350(H.projective(strict,x),x)]
    polys=[s.Poly(v,*variables,domain=s.QQ)for v in full];assert all(poly.total_degree()<=2 for poly in polys)
    targets=[pp[5]*pp[8],pp[5]**2,pp[8]**2];targetpolys=[s.Poly(v,*variables,domain=s.QQ)for v in targets]
    monomials=sorted(set(mon for poly in polys+targetpolys for mon,coef in poly.terms()),key=lambda mon:(sum(mon),mon))
    DM=DomainMatrix([[s.QQ.from_sympy(poly.coeff_monomial(mon))for mon in monomials]for poly in polys],(350,len(monomials)),s.QQ)
    reduced,pivots=DM.rref();rows=DM.transpose().rref()[1];basis=DM.extract(list(rows),list(range(len(monomials))));minor=basis.extract(list(range(len(rows))),list(pivots));det=minor.det();assert det!=0 or not rows
    reducedrows=reduced.to_list();rrefpolys=[s.factor(sum(s.QQ.to_sympy(coef)*s.prod(t**e for t,e in zip(variables,mon))for coef,mon in zip(reducedrows[i],monomials)))for i in range(len(pivots))]
    receipts=[];inverse=None
    for expr,poly in zip(targets,targetpolys):
        vec=[s.QQ.from_sympy(poly.coeff_monomial(mon))for mon in monomials];residual=list(vec)
        for i,pivot in enumerate(pivots):
            scalar=vec[pivot]
            if scalar!=0:residual=[a0-scalar*b0 for a0,b0 in zip(residual,reducedrows[i])]
        member=all(v==0 for v in residual);receipt={'target':str(expr),'in_original_equation_rowspan':member}
        if member:
            if inverse is None:inverse=minor.inv()
            coeff=DomainMatrix([[vec[j]for j in pivots]],(1,len(pivots)),s.QQ)*inverse
            targetdm=DomainMatrix([vec],(1,len(monomials)),s.QQ);assert coeff*basis==targetdm
            receipt['exact_original_row_lift']=[{'row':int(rows[j]),'coefficient':str(s.QQ.to_sympy(v))}for j,v in enumerate(coeff.to_list()[0])if v!=0]
        receipts.append(receipt)
    D=s.expand(sum((-1)**i*x[i]*(Pi[j,k]*Pi[l,m]-Pi[j,l]*Pi[k,m]+Pi[j,m]*Pi[k,l])for i in range(5)for j,k,l,m in [[t for t in range(5)if t!=i]]))
    result={'task_id':'G','all_exact_arithmetic_passed':True,'G_closed':False,'scope':'One entire new correct-minus fixed source E2eta-1/r=dy/x. Original M1kernel15+fullP2nine; nohomogeneousP2 restriction, noextra jet/minor profile. Original350quadratic support andactualPf retained. Product-rowspan target only structural necessary relation, notautomaticclosure.',
        'source_variables':list(map(str,variables)),'full150_original_Pi':{str((i,j)):str(P[i,j])for i,j in H.PAIRS},'full_positive_Euler_gauge':list(map(str,L)),
        'full150_c_gauged_Pi':{str((i,j)):str(Pi[i,j])for i,j in H.PAIRS},'full_affine_Qbar':{str((i,j)):str(Q[i,j])for i,j in combinations(range(4),2)},'full_c_linear_R':{str((i,j)):str(R[i,j])for i,j in combinations(range(4),2)},'full_affine_W':list(map(str,W)),
        'full350_original_Jacobi_equations':list(map(str,full)),'parameter_monomials':monomials,'equation_rowspace_rank':len(pivots),'equation_rowspace_basis_original_rows':list(rows),'equation_rowspace_pivot_columns':list(pivots),'equation_rowspace_minor_determinant':str(s.QQ.to_sympy(det)),
        'full_row_reduced_equations':list(map(str,rrefpolys)),'two_residual_direction_receipts':receipts,'full_actual_D':str(D),'all_c_degree_actual_D_coefficients':list(map(str,s.Poly(D,c).all_coeffs())),
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,archive]],'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n');print('fullsource_rank',len(pivots),'products',[(v['target'],v['in_original_equation_rowspan'])for v in receipts],flush=True)


def m3_triple_new_fixed_actual_factor_main():
    """Exact original-Jacobi ideal lift of actual divisor restriction: no cohom localization."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_M3_TRIPLE_NEW_FIXED_ACTUAL_FACTOR_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('three exact original-row factor lifts finite budget60s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    archive=ROOT/'research/classification/inputs/plane_m3_triple_new_fixed_factor_preappend_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==48 and all(after[name]==value for name,value in before.items())
    previous=ROOT/'reports/CLASSIFICATION_PLANE_M3_TRIPLE_NEW_FIXED_ORIGINAL_SUPPORT_20261108.json';old=json.loads(previous.read_text())
    x=s.symbols('a b c z w');a,b,c,z,w=x;variables=s.symbols(' '.join(old['source_variables']));local={str(v):v for v in list(x)+list(variables)};parse=lambda v:s.sympify(v,locals=local)
    pp=variables[:15];qq=variables[15:];equations=list(map(parse,old['full350_original_Jacobi_equations']));monomials=old['parameter_monomials'];rows=old['equation_rowspace_basis_original_rows'];pivots=old['equation_rowspace_pivot_columns']
    basis=DomainMatrix([[s.QQ.from_sympy(s.Poly(equations[i],*variables).coeff_monomial(mon))for mon in monomials]for i in rows],(len(rows),len(monomials)),s.QQ)
    minor=basis.extract(list(range(len(rows))),pivots);inverse=minor.inv();targets=[pp[14]*(pp[14]-pp[9]),pp[14]*pp[8],qq[2]-(pp[14]-pp[9])*(pp[10]+4*pp[2]-pp[4])/15-pp[3]*pp[8]/3]
    receipts=[];lifts=[]
    for target in targets:
        poly=s.Poly(target,*variables);vec=[s.QQ.from_sympy(poly.coeff_monomial(mon))for mon in monomials]
        coeff=DomainMatrix([[vec[j]for j in pivots]],(1,len(pivots)),s.QQ)*inverse;assert coeff*basis==DomainMatrix([vec],(1,len(monomials)),s.QQ)
        lift={int(rows[j]):s.QQ.to_sympy(value)for j,value in enumerate(coeff.to_list()[0])if value!=0}
        assert s.expand(sum(value*equations[i]for i,value in lift.items())-target)==0
        lifts.append(lift);receipts.append({'target':str(target),'exact_original350_row_lift':[{'row':i,'coefficient':str(value)}for i,value in lift.items()]})
    A,B,C=targets;D=parse(old['full_actual_D']);restriction=s.expand(D.subs(z,0));claimed=a*A*w**4+2*b*B*w**4+pp[14]*qq[2]*w**5
    assert s.expand(restriction-claimed)==0
    coefficients=[a*w**4+(pp[10]+4*pp[2]-pp[4])*w**5/15,2*b*w**4+pp[3]*w**5/3,pp[14]*w**5]
    assert s.expand(restriction-sum(v*t for v,t in zip(coefficients,targets)))==0
    total={}
    for scalar,lift in zip(coefficients,lifts):
        for i,value in lift.items():total[i]=s.expand(total.get(i,0)+scalar*value)
    total={i:value for i,value in total.items()if value!=0};assert s.expand(restriction-sum(value*equations[i]for i,value in total.items()))==0
    proof=ROOT/'research/classification/PLANE_M3_TRIPLE_NEW_FIXED_ACTUAL_FACTOR.md'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'whole_fixed_source_irreducible_quintic_empty':True,'G_closed':True,'scope':'Entire genuinely new corrected fixed source, full original15P1+9P2. Exact polynomial ideal identity D|z0 in original350Jacobi, no source coordinate localization; actualD has zfactor at every Poisson point. Othergenerictriple sources depend on separate corrected90/18modules proofs.',
        'all_source_variables':list(map(str,variables)),'three_equations_and_original_row_lifts':receipts,'full_actual_D_restriction':str(restriction),'three_equation_multiplier_coefficients':list(map(str,coefficients)),
        'exact_full_original350_D_restriction_lift':[{'row':i,'coefficient':str(value)}for i,value in total.items()],'all_three_row_identities_exact':True,'full_actual_D_ideal_identity_exact':True,'inverted_source_coordinates':[],
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in[HELPER,previous,proof,archive]],'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n');print('whole_fixed_actual_z_factor_PASS','row_lift_terms',len(total),flush=True)


def double_line_m2_nil_resonance_main():
    """Whole necessary p=-2/7 nil fiber: exact old kernel, actual two-torus, original support."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M2_NIL_RESONANCE_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('wholep=-2/7 originalsource finite budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    archive=ROOT/'research/classification/inputs/plane_double_line_m2_nil_resonance_preappend_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==49 and all(after[name]==value for name,value in before.items())
    previous=ROOT/'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M2_CONSTANT_DOUBLE_THIRTYFOUR_20261108.json';packet=json.loads(previous.read_text());row=next(v for v in packet['cases']if v['name']=='HA3_NIL1_WHOLE_P');pvalue=-s.Rational(2,7)
    def decode(v):
        if isinstance(v,str):return s.Rational(v)
        assert set(v)=={'numerator','denominator'}
        nu=sum(s.Rational(c)*pvalue**i for i,c in enumerate(v['numerator']));de=sum(s.Rational(c)*pvalue**i for i,c in enumerate(v['denominator']));assert de!=0;return nu/de
    x=s.Matrix(s.symbols('a b c u v'));a,b,c,u,v=x;mons=H.monomials(5,2)
    def fromvector(raw):
        assert len(raw)==150;P=s.zeros(5)
        for k,(i,j)in enumerate(H.PAIRS):P[i,j]=s.expand(sum(decode(raw[15*k+t])*s.prod(z**e for z,e in zip(x,mon))for t,mon in enumerate(mons)));P[j,i]=-P[i,j]
        return P
    P0=fromvector(row['full150_P0']);expected=s.zeros(5)
    for i,j,val in[(0,2,-a*a),(1,2,2*a*b),(2,3,(pvalue-3)*a*u),(1,4,a*u),(2,4,b*u+pvalue*a*v)]:expected[i,j]=val;expected[j,i]=-val
    assert P0==expected
    K=s.Matrix([[decode(value)for value in rr]for rr in row['weight1_certificate']['complete_kernel']]).T;M1=s.Matrix([[decode(value)for value in rr]for rr in row['complete350x30']]);assert K.shape==(30,11)and M1*K==s.zeros(350,11)
    free=[j for j in range(30)if j not in row['weight1_certificate']['columns']];assert K.extract(free,list(range(11)))==s.eye(11)
    determinant=decode(row['weight1_certificate']['minor_determinant']);assert determinant!=0
    B1=list(map(fromvector,packet['full30_weight1_basis']));B2=list(map(fromvector,packet['full9_weight2_basis']));kk=s.symbols('k0:11');qq=s.symbols('q0:9')
    cohom=K*s.Matrix(kk);P1=sum((B*val for B,val in zip(B1,cohom)),s.zeros(5));P2=sum((B*val for B,val in zip(B2,qq)),s.zeros(5))
    m=s.expand(P1[2,3]).coeff(v,2);T=s.expand(P1[0,4]).coeff(u,2);assert m==kk[5]and T==kk[6]
    raw=(P0+P1+P2).applyfunc(s.expand);unit={kk[5]:s.Integer(1),kk[6]:s.Integer(1)};Pi=raw.subs(unit,simultaneous=True);variables=[z for z in list(kk)+list(qq)if z not in unit]
    mm,tt=s.symbols('normal_m normal_T');diag=[1,mm/tt,mm**2/tt**2,tt**2/mm,tt];coordsub=dict(zip(x,[xx/dd for xx,dd in zip(x,diag)]));scale=tt**2/mm**2
    pushed=s.Matrix(5,5,lambda i,j:s.cancel(scale*diag[i]*diag[j]*P0[i,j].subs(coordsub,simultaneous=True)));assert pushed==P0
    assert s.cancel(diag[2]*diag[3]*scale/diag[4]**2*mm)==1 and s.cancel(diag[0]*diag[4]*scale/diag[3]**2*tt)==1
    def divisor(P):return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,n]-P[j,l]*P[k,n]+P[j,n]*P[k,l])for i in range(5)for j,k,l,n in [[r for r in range(5)if r!=i]]))
    assert s.expand(divisor(P0)-(pvalue-4)*a**3*u*u)==0
    original=[s.expand(z)for z in H.all350(H.projective(H.bilinear(Pi,Pi,x),x),x)];actual=divisor(Pi);current=list(original);D=actual;remaining=list(variables);rules={};rounds=[];inconsistent=False;final=[]
    for turn in range(19):
        polys=[s.Poly(z,*remaining,domain=s.QQ)for z in current]if remaining else[s.Poly(z,domain=s.QQ)for z in current]
        monomials=sorted(set(mon for poly in polys for mon,val in poly.terms()),key=lambda mon:(sum(mon),mon))
        if not monomials:break
        DM=DomainMatrix([[s.QQ.from_sympy(poly.coeff_monomial(mon))for mon in monomials]for poly in polys],(350,len(monomials)),s.QQ)
        reduced,pivots=DM.rref();rows=DM.transpose().rref()[1];reducedrows=reduced.to_list();rank=len(pivots);basis=DM.extract(list(rows),list(range(len(monomials))));minor=basis.extract(list(range(rank)),list(pivots));det=minor.det();assert det!=0 or rank==0
        final=[s.factor(sum(s.QQ.to_sympy(coef)*s.prod(t**e for t,e in zip(remaining,mon))for coef,mon in zip(reducedrows[i],monomials)))for i in range(rank)]
        linear=[z for z in final if s.Poly(z,*remaining).total_degree()<=1]if remaining else final
        receipt={'round':turn,'free_variables':list(map(str,remaining)),'original_equation_rowspace_rank':rank,'rowspace_basis_input350_rows':list(rows),'rowspace_pivot_columns':list(pivots),'rowspace_minor_determinant':str(s.QQ.to_sympy(det)),'all_reduced_equations':list(map(str,final)),'linear_consequences':[]}
        if not linear:rounds.append(receipt);break
        inverse=minor.inv()
        for consequence in linear:
            pol=s.Poly(consequence,*remaining)if remaining else s.Poly(consequence);vec=[s.QQ.from_sympy(pol.coeff_monomial(mon))for mon in monomials];lift=DomainMatrix([[vec[j]for j in pivots]],(1,rank),s.QQ)*inverse
            assert lift*basis==DomainMatrix([vec],(1,len(monomials)),s.QQ)
            sparse={int(rows[j]):s.QQ.to_sympy(val)for j,val in enumerate(lift.to_list()[0])if val!=0};assert s.expand(sum(val*current[i]for i,val in sparse.items())-consequence)==0
            receipt['linear_consequences'].append({'equation':str(consequence),'exact_input350_row_lift':[{'row':i,'coefficient':str(val)}for i,val in sparse.items()]})
        if any(not z.free_symbols for z in linear):inconsistent=True;receipt['nonzero_constant_consequence']=True;rounds.append(receipt);break
        LM=s.Matrix([[s.expand(z).coeff(t)for t in remaining]+[s.expand(z).subs(dict.fromkeys(remaining,0))]for z in linear]);rr,piv=LM.rref();assert len(remaining)not in piv
        sub={remaining[j]:s.expand(-rr[i,len(remaining)]-sum(rr[i,k]*remaining[k]for k in range(len(remaining))if k not in piv))for i,j in enumerate(piv)}
        assert sub and all(s.expand(z.subs(sub,simultaneous=True))==0 for z in linear)
        receipt['exact_affine_substitution']={str(k):str(v)for k,v in sub.items()};rounds.append(receipt);current=[s.expand(z.subs(sub,simultaneous=True))for z in current];D=s.expand(D.subs(sub,simultaneous=True));rules={key:s.expand(value.subs(sub,simultaneous=True))for key,value in rules.items()};rules.update(sub);remaining=[z for z in remaining if z not in sub]
    else:raise AssertionError('finite Gaussian elimination did not decrease source-variable count')
    proper=None
    if not inconsistent:
        for factor,power in s.factor_list(D)[1]:
            degree=s.Poly(factor,*x).total_degree()
            if 0<degree<5:proper=factor;assert s.denom(s.cancel(D/factor))==1;break
    upper=None if inconsistent else 16+len(remaining);closed=inconsistent or proper is not None or upper<=24
    proof=ROOT/'research/classification/PLANE_DOUBLE_LINE_M2_NIL_RESONANCE_ENTRANCE.md'
    result={'task_id':'G','all_exact_arithmetic_passed':True,'whole_p_fixed_value':str(pvalue),'whole_Tm_unit_G_closed':closed,'scope':'Only entire p=-2/7 nilsource T*m!=0, actualtwo-torus normalizesboth1 with inverseunits. Original complete11P1 evaluatedfromfrozen34kernel+FULL9P2; original350/Pf and exactlinear Gaussian quotient only. T0/m0 complements andwholep proof separate.',
        'source_first_symbol':{str((i,j)):str(P0[i,j])for i,j in H.PAIRS},'old_generic_kernel_specialized':serial_matrix(K),'old_generic_minor_value':str(determinant),'old_generic_kernel_free_identity_rows':free,
        'raw_m':str(m),'raw_T':str(T),'actual_normalizer_diagonal':list(map(str,diag)),'actual_tensor_scale':str(scale),'normalizer_P0_and_two_unit_identities_exact':True,
        'all150_unlocalized_original_source':{str((i,j)):str(raw[i,j])for i,j in H.PAIRS},'all150_normalized_original_source':{str((i,j)):str(Pi[i,j])for i,j in H.PAIRS},'normalized_source_variables':list(map(str,variables)),
        'full350_normalized_original_Jacobi':list(map(str,original)),'full_actual_D_before_elimination':str(actual),'exact_linear_elimination_rounds':rounds,'all_affine_elimination_rules':{str(k):str(v)for k,v in rules.items()},'remaining_source_variables':list(map(str,remaining)),
        'remaining_original_Jacobi_equations':list(map(str,current)),'full_actual_D_after_elimination':str(D),'whole_unit_source_inconsistent':inconsistent,'proper_polynomial_actual_D_factor':None if proper is None else str(proper),'whole_unit_image_upper':upper,'normalization_residual_group_fiber_lower':8,
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_bindings':[{'path':str(z.relative_to(ROOT)),'sha256':hashlib.sha256(z.read_bytes()).hexdigest()}for z in[HELPER,previous,proof,archive]],'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n');print('resonance',str(pvalue),'unitclosed',closed,'inconsistent',inconsistent,'remaining',len(remaining),'factor',proper,'upper',upper,flush=True)


def triple_weak_wall_global_graphs_main():
    """One faithful rank-four chart per entire weak-wall orbit; global T39/T40 graphs."""
    import ast
    from s61_model import S61,X,pushforward
    destination=ROOT/'reports/CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_GLOBAL_GRAPHS_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('two whole weak-orbit globalgraph evaluations budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    archive=ROOT/'research/classification/inputs/plane_triple_weak_wall_graphs_preappend_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==50 and all(after[name]==value for name,value in before.items())
    previous=ROOT/'reports/CLASSIFICATION_PLANE_M3_TRIPLE_VALID_ORE_MODULES_20261108.json';packet=json.loads(previous.read_text());family=next(v for v in packet['cases']if v['type']=='E1DX_1')
    coords=s.symbols('a b c z w');a,b,c,z,w=coords;ww=s.symbols('W0:18');y=[a,b,z,w];indices=[0,1,3,4];locals_map={str(v):v for v in list(coords)+[TAU]};parse=lambda v:s.sympify(v,locals=locals_map)
    def divisor(P):return s.expand(sum((-1)**i*coords[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]))
    model=S61();results=[]
    for value in[s.Rational(4,21),s.Rational(1,12)]:
        wall=next(v for v in family['all_exceptional_fields']if s.Rational(v['alpha'])==value);P=s.zeros(5)
        for i,j in combinations(range(4),2):P[indices[i],indices[j]]=parse(family['full_Q'][str((i,j))]).subs(TAU,value);P[indices[j],indices[i]]=-P[indices[i],indices[j]]
        W0=s.Matrix([parse(v).subs(TAU,value)for v in family['W0']]);W=W0+s.Matrix([z*(ww[0]*a+ww[1]*b)+w*(ww[2]*a+ww[3]*b)+ww[4]*z*z+ww[5]*z*w+ww[6]*w*w,
            z*(ww[7]*a+ww[8]*b)+w*(ww[9]*a+ww[10]*b)+ww[11]*z*z+ww[12]*z*w+ww[13]*w*w,
            ww[14]*z*w+ww[15]*w*w,ww[16]*z*z+ww[17]*z*w])
        W=W.subs(dict(zip(ww,map(s.Rational,wall['particular_solution']))),simultaneous=True)
        for i in range(4):P[2,indices[i]]=W[i];P[indices[i],2]=-W[i]
        D=divisor(P);assert s.expand(D-parse(wall['full_actual_D_particular']))==0
        evaluation=[0,0,0,0,1];E=D.subs(dict(zip(coords,evaluation)));assert E!=0
        original=P.subs(dict(zip(coords,X)),simultaneous=True);permutation=s.zeros(5)
        for i,j in enumerate([4,0,1,2,3]):permutation[i,j]=1
        moved=pushforward(original,permutation);K=moved.subs(dict(zip(X,[1,0,0,0,0])))[1:,1:];assert K[0,1]==1 and K.det()!=0
        unit=[s.eye(4)[:,i]for i in range(4)];r0=unit[0];r1=unit[1]/K[0,1]
        r2=unit[2]-(unit[2].T*K*r1)[0]*r0+(unit[2].T*K*r0)[0]*r1
        r3=unit[3]-(unit[3].T*K*r1)[0]*r0+(unit[3].T*K*r0)[0]*r1;r3=r3/s.cancel((r2.T*K*r3)[0])
        frame=s.eye(5);frame[1:,1:]=s.Matrix.vstack(*(v.T for v in[r0,r1,r2,r3]));J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1
        assert (frame[1:,1:]*K*frame[1:,1:].T-J).applyfunc(s.cancel)==s.zeros(4)
        G=frame*permutation;GI=G.inv();assert G*GI==s.eye(5);framed=pushforward(original,G);assert (pushforward(framed,GI)-original).applyfunc(s.cancel)==s.zeros(5)
        normalized=model.normalize(framed);verification=model.verify_point(normalized,check_equations=True);point=normalized['parameters'];us=s.symbols('u1:62');graphs={}
        for name,path in[('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
            data=json.loads((ROOT/path).read_text());ii=data['weight2_u_indices_zero_based'];residual=[s.factor(s.cancel(point[k]-s.sympify(v).subs(dict(zip(us,point)),simultaneous=True)))for k,v in zip(ii,data['quadratic_section_q'])];assert len(residual)==45
            graphs[name]={'in_full_global_component':all(v==0 for v in residual),'all45_graph_residuals':list(map(str,residual)),'nonzero_residuals':{str(k):str(v)for k,v in zip(ii,residual)if v!=0}}
        total=normalized['matrix']*G;assert total.det()!=0
        rec={'tau':str(value),'full_original_Pi':{str((i,j)):str(P[i,j])for i,j in H.PAIRS},'actual_D':str(D),'structural_rank4_point':evaluation,'actual_D_at_point':str(E),'constant_value_matrix':serial_matrix(K),'actual_permutation':serial_matrix(permutation),'Darboux_frame':serial_matrix(frame),'full_framing_matrix':serial_matrix(G),'full_inverse_framing_matrix':serial_matrix(GI),'perspective_matrix':serial_matrix(normalized['matrix']),'full_actual_PGL_matrix':serial_matrix(total),'full_inverse_actual_PGL_matrix':serial_matrix(total.inv()),'perspective':list(map(str,normalized['perspective'])),'all61_S61_parameters':list(map(str,point)),'faithful150_full134_verification':verification,'graph_membership':graphs,'whole_source_implication':'EntirecanonicalW6isoneactualPGLorbit by separatepureproof. Positiveclosedgraph giveswholefamilycontainment; negativeexcludeswholeglobalT39/T40 by acceptedconnectedSp4/Ga4evaluation-fiber theorem. Otherparents/componenthoodseparate.'}
        results.append(rec);print('weakwall',str(value),'globalgraphs',{name:v['in_full_global_component']for name,v in graphs.items()},flush=True)
    proofs=[ROOT/'research/classification'/name for name in['PLANE_M3_TRIPLE_WEAK_WALLS_TWO_ACTUAL_ORBITS.md','PLANE_M3_TRIPLE_WEAK_WALL_ORBIT_DIMENSIONS.md','S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md','S61_GLOBAL_COMPONENT_FIBER_SEPARATE_REVIEW.md']]
    paths=[HELPER,ROOT/'code/s61_model.py',previous,archive,ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json']+proofs
    paths.extend(ROOT/'legacy/v02/data'/name for name in['symplectic_slice_certificate.json','symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    result={'task_id':'E_reduced_nonnormal','all_exact_arithmetic_passed':True,'cases':results,'scope':'Exactly2 entire weak-orbit knownparent tests, one structurallyguaranteed symplecticevaluation each, all150/PGL inverse/full134 andall45graphrows. No additionalpoint/Jacobian/operator/jet/profile batch. Global T39/T40 result uses separatelyacceptedfullcomponentfiber theorem, notgenericchartguess. Otherparent/newcomponentunresolved.',
        'compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in paths],'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def triple_weak_wall_ambient_main():
    """One fullambient matrix per alreadyproved whole actual orbit, no new anchors."""
    import ast
    from collections import Counter
    from classification_three_hyperplanes_quadric_independent import ambient_matrices
    from classification_rank3_mixed_component_review import COLUMNS,ROWS
    destination=ROOT/'reports/CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_FULL_AMBIENT_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('twoentireorbits fullambient exactbudget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(3072*1024**2,)*2)
    archive=ROOT/'research/classification/inputs/plane_triple_weak_wall_ambient_preappend_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==51 and all(after[name]==value for name,value in before.items())
    previous=ROOT/'reports/CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_GLOBAL_GRAPHS_20261108.json';packet=json.loads(previous.read_text());x=s.symbols('a b c z w');mapping=dict(zip(map(str,x),x));weights=[3,2,0,5,4]
    wc=[sum(weights[i]for i in ij)-sum(v*q for v,q in zip(e,weights))+3 for ij,e in COLUMNS]
    wr=[sum(weights[i]for i in ij)-sum(v*q for v,q in zip(e,weights))+6 for ij,e in ROWS]
    def enc(A):return {'shape':list(A.shape),'entries':{str(i)+','+str(j):str(z)for(i,j),z in A.todok().items()}}
    records=[]
    for case in packet['cases']:
        P=s.zeros(5);scale=14 if case['tau']=='4/21'else 4;assert case['tau']in['4/21','1/12']
        for key,text in case['full_original_Pi'].items():
            i,j=ast.literal_eval(key);P[i,j]=s.expand(scale*s.sympify(text,locals=mapping));P[j,i]=-P[i,j]
        assert all(sum(weights[i]for i in ij)-sum(v*q for v,q in zip(e,weights))+3==0 for ij in H.PAIRS for e,z in s.Poly(P[ij],*x).terms()if z)
        D,V,M,unused=ambient_matrices(P,x,[P]);assert unused==M.row_join(M[:,-1]);assert all(wr[i]==wc[j]for(i,j),z in D.todok().items())
        certs={};kernels={}
        for name,A in [('full350_D',D),('Euler25',V),('Euler_GL5_scalar51',M)]:
            cert,det,K=rank_certificate(A,variable=False);certs[name]=cert;kernels[name]=K
        assert certs['Euler25']['rank']==24 and D*M==s.zeros(350,51)
        kernel=kernels['full350_D'];free=[j for j in range(150)if j not in certs['full350_D']['minor_columns']]
        assert kernel[free,:]==s.eye(len(free))
        assert all(wc[i]==wc[free[j]]for(i,j),z in kernel.todok().items())
        group_weights=[]
        for j in range(51):
            support={wc[i]for i in range(150)if M[i,j]};assert len(support)<=1;group_weights.append(next(iter(support))if support else None)
        kcounts=Counter(wc[j]for j in free);gcounts=Counter(group_weights[j]for j in certs['Euler_GL5_scalar51']['minor_columns']);normal={str(k):kcounts[k]-gcounts[k]for k in sorted(set(kcounts)|set(gcounts))if kcounts[k]!=gcounts[k]};assert all(v>=0 for v in normal.values())
        ranks={name:cert['rank']for name,cert in certs.items()};orbit=ranks['Euler_GL5_scalar51']-25;tangent=125-ranks['full350_D'];assert orbit==23
        determinant=s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[r for r in range(5)if r!=i]]));assert s.expand(determinant-scale**2*s.sympify(case['actual_D'],locals=mapping))==0
        records.append({'tau':case['tau'],'integer_scale':scale,'coordinate_order':list(map(str,x)),'complete_actual_brackets':{str((i,j)):str(P[i,j])for i,j in H.PAIRS},'actual_signed_Euler_Pfaffian':str(determinant),'ranks_QQ':ranks,'all350_Jacobi_zero':True,'all_D_times_actual51_zero':True,'complete_sparse_matrices':{name:enc(A)for name,A in [('full350_D',D),('Euler25',V),('Euler_GL5_scalar51',M)]},'certificates':certs,'full_D_kernel_free_columns':free,'full_D_kernel_sparse':enc(kernel),'kernel_all350_substitution_and_free_identity_verified':True,'actual_projective_orbit_dimension':orbit,'projective_tangent_dimension':tangent,'smooth_unique_component_proved':orbit==tangent,'canonical_orbitclosure_full_component':orbit==tangent,'normalized_actual_torus':{'coordinate_weights':weights,'tensor_scale_exponent':3,'all150_point_identity':True},'all150_tangent_weights':wc,'all350_equation_weights':wr,'all51_group_column_weights':group_weights,'kernel_weight_dimensions':dict(sorted(kcounts.items())),'group_weight_dimensions':dict(sorted(gcounts.items())),'normal_quotient_weight_dimensions':normal,'next_mathematical_branch':'D102 smooth23candidate; D101 unique49d kernel/48actualgroup, finiteWHOLEequivariantslice possibleonlyifnormalcharacter nonzero; D<=100 stoplocalrouteforactualparentgeometry.'})
        print('weakambient',case['tau'],ranks,'orbit',orbit,'tangent',tangent,'normal',normal,flush=True)
    paths=[HELPER,ROOT/'code/classification_three_hyperplanes_quadric_independent.py',ROOT/'code/classification_rank3_mixed_component_review.py',previous,archive,ROOT/'research/classification/PLANE_M3_TRIPLE_WEAK_WALLS_TWO_ACTUAL_ORBITS.md',ROOT/'research/classification/PLANE_M3_TRIPLE_WEAK_WALL_ORBIT_DIMENSIONS.md',ROOT/'research/classification/PLANE_TRIPLE_WEAK_WALL_GLOBAL_PARENT_EXCLUSIONS.md']
    result={'task_id':'E_reduced_nonnormal','all_exact_arithmetic_passed':True,'cases':records,'full_original_coefficient_order':[{'pair':list(ij),'exponent':list(e)}for ij,e in COLUMNS],'full_original_equation_order':[{'quad':list(ij),'exponent':list(e)}for ij,e in ROWS],'scope':'ONE fullambientD/Euler/full51 per TWOentirealreadyprovedactualorbits, integral scale only, fullrationalkernel/nonzerominors. No alternativepoint/profile/normal/jet. Componentdistinctness or24parentrecognition separate. CompletedGjobsreadonly.','compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'input_bindings':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in paths],'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


def triple_y_walls_ore_main():
    """ExactEuler/Ore/inhomogeneousP2 from the TWO existing9dim kernels only."""
    import ast
    destination=ROOT/'reports/CLASSIFICATION_PLANE_TRIPLE_Y_WALLS_ORE_20261108.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('two saved9kernelwholeYsource exactOre budget180s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(3072*1024**2,)*2)
    archive=ROOT/'research/classification/inputs/plane_triple_y_walls_ore_preappend_20261108.json';snapshot=json.loads(archive.read_text())
    functions=lambda text:{n.name:ast.dump(n,include_attributes=False)for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
    before=functions(snapshot['runner_source']);after=functions(Path(__file__).read_text());assert len(before)==52 and all(after[name]==value for name,value in before.items())
    previous=ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json';packet=json.loads(previous.read_text());x=s.Matrix(s.symbols('a b c z w'));a,b,c,z,w=x;parse=lambda v:s.sympify(v,locals={str(z):z for z in[*x,TAU]});p=s.symbols('p0:9');r=s.symbols('r0:9');P1raw,src1,P2raw,src2=sources(x);ys=[0,1,3,4]
    records=[]
    def jac(P):return {q:s.expand(v/2)for q,v in H.projective(H.bilinear(P,P,x),x).items()}
    def actualD(P):return s.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[t for t in range(5)if t!=i]]))
    for name,alpha in [('E1_LOG_GENERIC_Y',s.Integer(-1)),('E2_LOG_GENERIC_Y',-s.Rational(4,3))]:
        family=next(c for c in packet['cases']if c['type']==name);wall=next(c for c in family['all_exceptional_sources']if parse(c['alpha'])==alpha);assert wall['weight1_certificate']['nullity']==9 and wall['weight2_certificate']['nullity']==0
        K=s.Matrix([[parse(v)for v in row]for row in wall['weight1_certificate']['complete_kernel']]);B=s.Matrix([[parse(v)for v in row]for row in wall['full350_weight2_operator']]);assert K.shape==(30,9)and B.shape==(350,9)and K[2,:]==s.zeros(1,9)and K[11,:]==s.zeros(1,9)
        P0=s.zeros(5)
        for key,text in wall['full150_P0'].items():i,j=ast.literal_eval(key);P0[i,j]=parse(text);P0[j,i]=-P0[i,j]
        P1=P1raw.subs(dict(zip(src1,K*s.Matrix(p))),simultaneous=True).applyfunc(s.expand);P2=P2raw.subs(dict(zip(src2,r)),simultaneous=True);raw=P0+P1+P2
        R=s.Matrix([s.diff(raw[2,j],c)if j!=2 else 0 for j in range(5)]);assert all(s.diff(v,c)==0 for v in R)
        ore=(raw-H.wedge(x,R)).applyfunc(s.expand);assert all(s.diff(v,c)==0 for v in ore)
        Q=ore.extract(ys,ys);W=s.Matrix([ore[2,j]for j in ys]);yf=[x[i]for i in ys]
        strictQ={tri:s.expand(sum(Q[i,l]*s.diff(Q[j,k],yf[l])for i,j,k in[tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(4)))for tri in combinations(range(4),3)}
        original=jac(raw);oreJ=jac(ore);assert all(s.expand(original[q]-oreJ[q])==0 for q in original)
        for tri in combinations(range(4),3):
            quad=tuple(sorted([2]+[ys[i]for i in tri]));coefficient=s.expand(oreJ[quad]).coeff(c)
            sign=(-1)**quad.index(2);assert s.expand(coefficient-sign*strictQ[tri])==0
        # The full nine P2 variables are retained up to an exact constant invertible nine-row elimination.
        q11=[parse(v)/2 for v in H.all350(H.projective(H.bilinear(P1,P1,x),x),x)]
        rr=wall['weight2_certificate']['minor_rows'];cc=wall['weight2_certificate']['minor_columns'];assert cc==list(range(9))
        block=B.extract(rr,cc);det=block.det(method='domain-ge');assert det==parse(wall['weight2_certificate']['minor_determinant'])and det!=0
        particular=(-block.inv()*s.Matrix([q11[i]for i in rr])).applyfunc(s.expand);compat=(B*particular+s.Matrix(q11)).applyfunc(s.expand)
        replacement=dict(zip(r,particular));full=raw.subs(replacement,simultaneous=True).applyfunc(s.expand);Qf=Q.subs(replacement,simultaneous=True).applyfunc(s.expand);Wf=W.subs(replacement,simultaneous=True).applyfunc(s.expand)
        fullJ=H.all350(jac(full),x);remaining=list(dict.fromkeys(str(s.factor(parse(v)))for v in fullJ if parse(v)!=0));D=actualD(full)
        pf=s.factor(Qf[0,1]*Qf[2,3]-Qf[0,2]*Qf[1,3]+Qf[0,3]*Qf[1,2])
        qstrict=[s.factor(v.subs(replacement,simultaneous=True))for v in strictQ.values()]
        rec={'type':name,'alpha':str(alpha),'all9_P1_parameters':list(map(str,p)),'FULL9_P2_parameters':list(map(str,r)),'saved_full30by9_kernel':serial_matrix(K),'saved_P0':{str((i,j)):str(P0[i,j])for i,j in H.PAIRS},'full_P1':{str((i,j)):str(P1[i,j])for i,j in H.PAIRS},'full_P2_before_elimination':{str((i,j)):str(P2[i,j])for i,j in H.PAIRS},'actual_Euler_gauge_R':list(map(str,R)),'complete_literal_Ore_Q_before_P2_elimination':serial_matrix(Q),'complete_literal_Ore_W_before_P2_elimination':list(map(str,W)),'whole_c_independence_and_all350_Euler_identity_verified':True,'full_strict_Q_Jacobi_before_elimination':{str(tri):str(v)for tri,v in strictQ.items()},'all_strictQ_equations_are_original_c_coefficients':True,'same_existing_P2_minor_rows':rr,'same_existing_P2_minor_determinant':str(det),'full9_inhomogeneous_P2_unique_particular':list(map(str,particular)),'all350_weight2_compatibility_polynomials':list(map(str,compat)),'all350_original_Jacobi_after_elimination':fullJ,'all_distinct_original_support_polynomials':remaining,'complete_Q_after_unique_P2_elimination':serial_matrix(Qf),'complete_W_after_unique_P2_elimination':list(map(str,Wf)),'strictQ_Jacobi_after_elimination':list(map(str,qstrict)),'full_Q_Pfaffian_after_elimination':str(pf),'full_actual_D_after_elimination':str(D),'whole_support_closed':False,'scope':'Exactwhole9P1/FULL9P2 input/equivalence andoriginalsupport, no newM1/point/Jacobian/minorprofile; support recognition requires nextpureproof.'}
        records.append(rec);print('YwholeOre',name,str(alpha),'support',remaining,'Qpf',pf,'D',s.factor(D),flush=True)
    paths=[HELPER,previous,archive,ROOT/'research/classification/PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_CONNECTIONS_V2.md',ROOT/'research/classification/PLANE_M3_FREE_CONNECTION_SIGN_ERRATUM.md',ROOT/'research/classification/PLANE_DEPENDENT_ORE_SEPARATE_REVIEW.md']
    result={'task_id':'E_reduced_nonnormal','all_exact_arithmetic_passed':True,'cases':records,'scope':'Twoexistingdimension-onlyYwalls, full9P1andFULL9inhomogeneousP2, actualEulerwholeOre inverse/constantP2elimination/originalall350. No newM1/normal/jet/point/minorprofiles. NoEclosureclaimeduntilwholepolynomialsupport analyzed.','compatible_extension_old_function_AST_unchanged':True,'old_function_names':sorted(before),'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'input_bindings':[{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()}for v in paths],'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    try:
        if '--triple-y-walls-ore' in __import__('sys').argv:triple_y_walls_ore_main()
        elif '--triple-weak-wall-ambient' in __import__('sys').argv:triple_weak_wall_ambient_main()
        elif '--triple-weak-wall-global-graphs' in __import__('sys').argv:triple_weak_wall_global_graphs_main()
        elif '--double-line-m2-nil-resonance' in __import__('sys').argv:double_line_m2_nil_resonance_main()
        elif '--m3-triple-new-fixed-factor' in __import__('sys').argv:m3_triple_new_fixed_actual_factor_main()
        elif '--m3-triple-new-fixed-support' in __import__('sys').argv:m3_triple_new_fixed_original_support_main()
        elif '--m3-triple-valid-ore-modules' in __import__('sys').argv:m3_triple_valid_ore_modules_main()
        elif '--nonreduced-m3-triple-vertex-b' in __import__('sys').argv:nonreduced_m3_triple_vertex_b_main()
        elif '--nonreduced-m3-triple-corrected-ore' in __import__('sys').argv:nonreduced_m3_triple_corrected_ore_main()
        elif '--nonreduced-m3-triple-corrected-v2' in __import__('sys').argv:nonreduced_m3_triple_corrected_v2_main()
        elif '--nonreduced-m3-triple-ore' in __import__('sys').argv:nonreduced_m3_triple_ore_main()
        elif '--nonreduced-m3-triple' in __import__('sys').argv:nonreduced_m3_triple_main()
        elif '--nonreduced-m3-moving-simple-modules' in __import__('sys').argv:nonreduced_m3_moving_simple_modules_main()
        elif '--nonreduced-m3-moving-simple-ore' in __import__('sys').argv:nonreduced_m3_moving_simple_ore_main()
        elif '--nonreduced-m3-moving-simple' in __import__('sys').argv:nonreduced_m3_moving_simple_main()
        elif '--nonreduced-m3-moving-double-ore' in __import__('sys').argv:nonreduced_m3_moving_double_ore_main()
        elif '--nonreduced-m3-moving-double' in __import__('sys').argv:nonreduced_m3_moving_double_main()
        elif '--nonreduced-m3-constant-modules' in __import__('sys').argv:nonreduced_m3_constant_modules_main()
        elif '--nonreduced-m3-constant-ore' in __import__('sys').argv:nonreduced_m3_constant_ore_main()
        elif '--nonreduced-m3-constant-roots' in __import__('sys').argv:nonreduced_m3_constant_roots_main()
        elif '--double-line-m3-ore' in __import__('sys').argv:double_line_m3_ore_main()
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
        failure=ROOT/('reports/CLASSIFICATION_PLANE_TRIPLE_Y_WALLS_ORE_RUNTIME_20261108.json' if '--triple-y-walls-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_FULL_AMBIENT_RUNTIME_20261108.json' if '--triple-weak-wall-ambient' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_GLOBAL_GRAPHS_RUNTIME_20261108.json' if '--triple-weak-wall-global-graphs' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M2_NIL_RESONANCE_RUNTIME_20261108.json' if '--double-line-m2-nil-resonance' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_M3_TRIPLE_NEW_FIXED_ACTUAL_FACTOR_RUNTIME_20261108.json' if '--m3-triple-new-fixed-factor' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_M3_TRIPLE_NEW_FIXED_ORIGINAL_SUPPORT_RUNTIME_20261108.json' if '--m3-triple-new-fixed-support' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_M3_TRIPLE_VALID_ORE_MODULES_RUNTIME_20261108.json' if '--m3-triple-valid-ore-modules' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_VERTEX_B_RUNTIME_20261108.json' if '--nonreduced-m3-triple-vertex-b' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_ORE_RUNTIME_20261108.json' if '--nonreduced-m3-triple-corrected-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_RUNTIME_20261108.json' if '--nonreduced-m3-triple-corrected-v2' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_ORE_RUNTIME_20261108.json' if '--nonreduced-m3-triple-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_RUNTIME_20261108.json' if '--nonreduced-m3-triple' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_MODULES_RUNTIME_20261108.json' if '--nonreduced-m3-moving-simple-modules' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_ORE_RUNTIME_20261108.json' if '--nonreduced-m3-moving-simple-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_RUNTIME_20261107.json' if '--nonreduced-m3-moving-simple' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_RUNTIME_20261107.json' if '--nonreduced-m3-moving-double-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_RUNTIME_20261107.json' if '--nonreduced-m3-moving-double' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_MODULES_RUNTIME_20261107.json' if '--nonreduced-m3-constant-modules' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_ORE_RUNTIME_20261107.json' if '--nonreduced-m3-constant-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_ROOTS_RUNTIME_20261107.json' if '--nonreduced-m3-constant-roots' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_ORE_RUNTIME_20261107.json' if '--double-line-m3-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_RUNTIME_20261107.json' if '--double-line-m3-squarefree' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_RUNTIME_20261107.json' if '--m4-nonreduced-squarefree' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_RUNTIME_20261107.json' if '--m4-two-line-casimir' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_LEADING_RUNTIME_20261107.json' if '--nonreduced-cubic-leading' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_RUNTIME_20261107.json' if '--nonreduced-cubic' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ORE_RUNTIME_20261107.json' if '--global-ore' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_DEPENDENT_RUNTIME_20261107.json' if '--global-dependent' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMMONIMAGE_RUNTIME_20261107.json' if '--global-commonimage' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_CARTAN_SURFACE_RUNTIME_20261107.json' if '--global-cartan-surface' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_THREEPARAMETER_RUNTIME_20261107.json' if '--global-threeparameter' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_NIL_LINE_RUNTIME_20261107.json' if '--global-nil-line' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_HYPERBOLA_RUNTIME_20261107.json' if '--global-hyperbola' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TWOPARAMETER_RUNTIME_20261107.json' if '--global-twoparameter' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TSQUARE_RUNTIME_20261107.json' if '--global-tsquare' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_RESIDUAL_RUNTIME_20261107.json' if '--global-residual' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMPATIBILITY_RUNTIME_20261107.json' if '--global-compatibility' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ONECUT_RUNTIME_20261107.json' if '--global-onecut' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_RUNTIME_20261107.json' if '--global-lowparameter' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-nonkernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_RUNTIME_20261107.json' if '--rank2-kernel' in __import__('sys').argv else 'reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_RUNTIME_20261107.json')
        if not failure.exists():
            failure.write_text(json.dumps({'task_id':'G','all_passed':False,'failure':repr(error),
                'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'producer_bytes':Path(__file__).read_text(),'scope':'Stopped bounded symbolic linear-operator job; no mathematical G conclusion.'},indent=2)+'\n')
        raise
