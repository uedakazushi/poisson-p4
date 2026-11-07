#!/usr/bin/env python3
"""Whole rankM2 high-m plane symbols: seven finite operators and one reuse.

All39 higher coefficients, exact univariate minor complements and full
characteristic-zero kernels are retained. No shared producer is imported or
edited, and no nonlinear GB or point scan is performed.
"""
import argparse
from fractions import Fraction as F
import hashlib,itertools as it,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]
TAU=s.Symbol('tau');POLY=s.QQ.poly_ring(TAU)


def encode(M):return [[str(s.factor(v))for v in row]for row in M.tolist()]


def certificate(M,domain):
    dm=DomainMatrix.from_Matrix(M).convert_to(domain).to_field();rr,cols=dm.rref();rows=dm.transpose().rref()[1]
    field=dm.domain;R=rr.to_list();free=[j for j in range(M.cols)if j not in cols]
    K=[[field.zero for j in free]for i in range(M.cols)]
    for k,j in enumerate(free):
        K[j][k]=field.one
        for i,p in enumerate(cols):K[p][k]=-R[i][j]
    assert not any(any(row)for row in (dm*DomainMatrix(K,(M.cols,len(free)),field)).to_list())
    det=DomainMatrix.from_Matrix(M.extract(list(rows),list(cols))).convert_to(domain).det()
    expression=domain.to_sympy(det);assert expression!=0 or len(cols)==0
    return dict(rank=len(cols),nullity=len(free),minor_rows=list(rows),minor_columns=list(cols),minor_determinant=str(s.factor(expression)),
                complete_kernel=[[str(s.factor(field.to_sympy(v)))for v in row]for row in K],full_kernel_substitution_zero=True),expression


def algebraic_certificate(M,factor):
    alpha=s.CRootOf(factor,0);field=s.QQ.algebraic_field(alpha)
    values=[[field.from_sympy(v.subs(TAU,alpha))for v in row]for row in M.tolist()]
    dm=DomainMatrix(values,M.shape,field);rr,cols=dm.rref();rows=dm.transpose().rref()[1];R=rr.to_list();free=[j for j in range(M.cols)if j not in cols]
    K=[[field.zero for j in free]for i in range(M.cols)]
    for k,j in enumerate(free):
        K[j][k]=field.one
        for i,p in enumerate(cols):K[p][k]=-R[i][j]
    assert not any(any(row)for row in (dm*DomainMatrix(K,(M.cols,len(free)),field)).to_list())
    value=dm.extract(list(rows),list(cols)).det();assert value!=field.zero or not cols
    return dict(rank=len(cols),nullity=len(free),minor_rows=list(rows),minor_columns=list(cols),minor_determinant=str(field.to_sympy(value)),
                complete_kernel=[[str(field.to_sympy(v))for v in row]for row in K],full_kernel_substitution_zero=True,
                exact_number_field=str(field),all_conjugate_roots_covered=True)


def wedge(U,V):return q.antisymmetric({(i,j):q.add(q.mul(U[i],V[j]),q.scale(q.mul(U[j],V[i]),-1))for i,j in q.PAIRS})
def plus(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in q.PAIRS})


def actual(P):
    out={}
    for missing in range(5):
        i,j,k,l=[v for v in range(5)if v!=missing]
        pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
        out=q.add(out,q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
    return out


def coefficient(poly,spatial):
    out=0
    for e,c in poly.items():
        assert not e[6]and not e[7]
        if e[:5]==spatial:out+=s.Rational(c.numerator,c.denominator)*TAU**e[5]
    return s.expand(out)


def operator(P0,bases):
    columns=[]
    for P in bases:
        J=q.jacobi(plus(P0,P))
        columns.append([coefficient(J[quad],mon)for quad,mon in q.ROWS])
    return s.Matrix(columns).T


def analyze(D1,D2,affine_dimension,variable):
    c1,m1=certificate(D1,POLY if variable else s.QQ);c2,m2=certificate(D2,POLY if variable else s.QQ)
    d=c1['nullity']+c2['nullity'];offset=affine_dimension+5
    out=dict(D1=c1,D2=c2,higher_fiber_upper=d,whole_projective_image_upper=offset+d,G24_closed=offset+d<=24)
    if not variable:return out
    product=s.factor(m1*m2*TAU);factors=s.factor_list(product)[1];walls=[]
    for factor,multiplicity in factors:
        poly=s.Poly(factor,TAU)
        if poly.degree()==1:
            value=-poly.nth(0)/poly.nth(1)
            r1,_=certificate(D1.subs(TAU,value),s.QQ);r2,_=certificate(D2.subs(TAU,value),s.QQ)
            parameter=str(value)
        else:
            r1=algebraic_certificate(D1,factor);r2=algebraic_certificate(D2,factor);parameter='all roots of '+str(factor)
        d=r1['nullity']+r2['nullity'];image=offset-1+d
        walls.append(dict(minimal_polynomial=str(factor),multiplicity=int(multiplicity),parameter=parameter,all_complex_roots_covered=True,
                          D1=r1,D2=r2,higher_fiber_upper=d,fixed_parameter_affine_symbol_dimension=affine_dimension-1,
                          whole_projective_image_upper=image,G24_closed=image<=24))
    out.update(exception_polynomial=str(product),all_irreducible_factors=[[str(v),int(m)]for v,m in factors],all_exceptional_roots=walls,
               zero_parameter_explicitly_retained=True,all_strata_G24_closed=out['G24_closed']and all(r['G24_closed']for r in walls))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole rankM2 high-m operator budget120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    X=[q.atom(1,i)for i in range(5)];a,b,c,z,w=X;t=q.atom(1,5);zero={};En=[zero,zero,zero,z,w]
    U=lambda i:[q.atom(int(j==i))for j in range(5)]
    VN=[zero,zero,zero,q.mul(c,w),q.mul(b,z)]
    VB=[zero,zero,zero,q.add(q.mul(b,z),q.mul(c,w)),q.scale(q.mul(b,w),-1)]
    first=[]
    for normal in [3,4]:
        for ij in it.combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for mon,j in [([3,3],4),([3,4],3),([3,4],4),([4,4],3)]:
        for i in range(3):first.append(q.antisymmetric({(i,j):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in it.combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    models=[('NDN0',U(1),VN,[zero]*5,q.scale(q.mul(z,z),-F(1,2)),9,False,'M0 nonzero nilpotent, nondegenerate W'),
            ('BSS0_tau',U(1),VB,[q.scale(q.mul(t,a),-F(2,3)),q.scale(q.mul(t,b),F(4,3)),q.scale(q.mul(t,c),-F(2,3)),zero,zero],q.mul(q.add(q.atom(1),q.scale(t,-2)),q.mul(z,w)),10,True,'tau is whole Lax R/M0 parameter'),
            ('BN0_tau',U(2),VB,[zero,zero,q.scale(q.mul(t,b),2),zero,zero],q.mul(q.add(q.atom(F(1,2)),q.scale(t,-1)),q.mul(w,w)),9,True,'tau is whole Lax R/M0 parameter'),
            ('R0_ND_nonzero_row_rho',U(0),VN,[q.add(b,q.mul(t,c)),zero,zero,zero,zero],zero,10,True,'rho=rb*rc; rho0 is nonzero isotropic axis'),
            ('R0_ND_zero_row',U(0),VN,[zero]*5,zero,8,False,'both original row coefficients zero'),
            ('R0_B_rc_nonzero',U(0),VB,[c,zero,zero,zero,zero],zero,9,False,'whole open rc!=0 normalized by actual Borel action'),
            ('R0_B_rcline_rho',U(0),VB,[q.mul(t,b),zero,zero,zero,zero],zero,8,True,'rc0 with rho=rb allC including bothrow0')]
    records=[]
    # Reuse a previously excluded but fully checked exact point. Do not rerun
    # the immutable m2 producer or reinterpret its former scope as a theorem.
    oldpath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_20261107.json';old=json.loads(oldpath.read_text())
    ndss=next(r for r in old['cases']if r['type']=='NDSS');wall=next(r for r in ndss['analysis']['all_exceptional_roots']if r['parameter']=='0')
    assert wall['D1']['nullity']==8 and wall['D2']['nullity']==0 and wall['excluded_multiplicity_higher']
    records.append(dict(type='NDSS0_reused',affine_first_symbol_dimension=10,source_report=str(oldpath.relative_to(ROOT)),source_report_sha256=hashlib.sha256(oldpath.read_bytes()).hexdigest(),
                        exact_source_location='cases[type=NDSS].analysis.all_exceptional_roots[parameter=0]',old_scope_exclusion_retained_but_new_geometric_application=True,
                        old_full_certificate=wall,higher_fiber_upper=8,whole_projective_image_upper=23,G24_closed=True))
    for name,u,V,T,K,dim,variable,meaning in models:
        P0=plus(wedge(T,En),wedge(u,V),q.antisymmetric({(3,4):K}))
        assert all(not v for v in q.jacobi(P0).values())and not actual(P0)
        D1=operator(P0,first);D2=operator(P0,second);analysis=analyze(D1,D2,dim,variable)
        row=dict(type=name,affine_first_symbol_dimension=dim,variable=variable,parameter_meaning=meaning,
                 full150_P0=[str(coefficient(P0.get(ij,{}),mon))for ij,mon in q.COLS],P0_full350_Jacobi_zero=True,P0_actual_Pf_zero=True,
                 full350_by30_D1=encode(D1),full350_by9_D2=encode(D2),analysis=analysis)
        records.append(row)
        Path('/tmp/plane_rank_zero_rank_two_high_m_progress.json').write_text(json.dumps(dict(completed=records,seconds=time.monotonic()-begun),indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',oldpath,
            ROOT/'research/classification/PLANE_RANK_ZERO_RANK_TWO_HIGH_MULTIPLICITY_FIRST_JET.md',
            ROOT/'research/classification/REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md',ROOT/'research/classification/REDUCED_DEGENERACY_ZERO_SURFACE_SEPARATE_REVIEW.md']
    result=dict(task_id='G',all_arithmetic_passed=True,all_eight_whole_high_m_rank2_strata=records,
                all_G24_closed=all(r.get('G24_closed',r.get('analysis',{}).get('all_strata_G24_closed',r.get('analysis',{}).get('G24_closed',False)))for r in records),
                all30_source_basis=[[str(v)for v in q.vec(P)]for P in first],all9_source_basis=[[str(v)for v in q.vec(P)]for P in second],full350_row_labels=[dict(quad=ij,monomial=mon)for ij,mon in q.ROWS],
                method='One fixed/univariate whole-source operator calculation, all irreducible minor roots exact, one immutable old excluded certificate reused. No nonlinear GB or shear-injectivity inference needed.',
                scope='Only reduced-actualD rankzero plane k1 rankM2 high-m; rankM3, rankM<=1, k0, full E and global classification are not asserted.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_source_bytes=Path(__file__).read_text(),
                seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_G24_closed'],seconds=result['seconds'],cases=[dict(type=r['type'],analysis=r.get('analysis',{}).get('higher_fiber_upper'),image=r.get('analysis',{}).get('whole_projective_image_upper',r.get('whole_projective_image_upper')))for r in records])))


def pfaffian_flag_main():
    """One exact original-Jacobi rowspace test on the three remaining walls."""
    import ast
    ap=argparse.ArgumentParser(description=pfaffian_flag_main.__doc__);ap.add_argument('--pfaffian-flag',action='store_true');ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite threewholewall original-Jacobi flag test120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    path=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_HIGH_M_20261107.json';old=json.loads(path.read_text())
    frozen={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(old['producer_source_bytes']).body if isinstance(v,(ast.FunctionDef,ast.AsyncFunctionDef))}
    current={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(Path(__file__).read_text()).body if isinstance(v,(ast.FunctionDef,ast.AsyncFunctionDef))}
    assert all(current[k]==v for k,v in frozen.items())
    def decode(vector):
        out={}
        for value,(ij,e)in zip(vector,q.COLS):
            value=F(value)
            if value:out.setdefault(ij,{})[e+(0,0,0)]=value
        return q.antisymmetric(out)
    first=list(map(decode,old['all30_source_basis']));second=list(map(decode,old['all9_source_basis']))
    saved=old['all_eight_whole_high_m_rank2_strata'][-1]
    dmon=[e for e in q.monomials(5)if e[4]==0]
    results=[]
    for value in [F(1),F(5,3),F(3)]:
        wall=next(r for r in saved['analysis']['all_exceptional_roots']if r['parameter']==str(value))
        P0=decode([str(s.sympify(v).subs(TAU,s.Rational(value.numerator,value.denominator)))for v in saved['full150_P0']])
        K=[[F(v)for v in row]for row in wall['D1']['complete_kernel']]
        family=[q.antisymmetric({ij:q.add(*(q.scale(first[i].get(ij,{}),K[i][j])for i in range(30)if K[i][j]))for ij in q.PAIRS})for j in range(len(K[0]))]+second
        n=len(family);pairs=list(it.combinations_with_replacement(range(n),2))
        ownJ=[q.jvec(q.jacobi(P))for P in family];ownD=[actual(P)for P in family];JC=[];DC=[]
        for i,P in enumerate(family):
            JC.append([a-b for a,b in zip(q.jvec(q.jacobi(plus(P0,P))),ownJ[i])]);DC.append(q.add(actual(plus(P0,P)),q.scale(ownD[i],-1)))
        for i,j in pairs:
            if i==j:JC.append(ownJ[i]);DC.append(ownD[i])
            else:
                JC.append([a-b-c for a,b,c in zip(q.jvec(q.jacobi(plus(family[i],family[j]))),ownJ[i],ownJ[j])])
                DC.append(q.add(actual(plus(family[i],family[j])),q.scale(ownD[i],-1),q.scale(ownD[j],-1)))
        M=list(map(list,zip(*JC)));target=[[P.get(e+(0,0,0),F())for P in DC]for e in dmon]
        def vector_sum(A,B,c):
            out=dict(A)
            for k,v in B.items():out[k]=out.get(k,F())+c*v
            return {k:v for k,v in out.items()if v}
        basis={};lifts={}
        for i,row in enumerate(M):
            v={j:c for j,c in enumerate(row)if c};lift={i:F(1)}
            for p in sorted(basis):
                c=v.get(p,F())
                if c:v=vector_sum(v,basis[p],-c);lift=vector_sum(lift,lifts[p],-c)
            if v:
                p=min(v);c=v[p];basis[p]={j:a/c for j,a in v.items()};lifts[p]={j:a/c for j,a in lift.items()}
        identities=[]
        for mon,row in zip(dmon,target):
            if not any(row):continue
            v={j:c for j,c in enumerate(row)if c};lift={}
            for p in sorted(basis):
                c=v.get(p,F())
                if c:v=vector_sum(v,basis[p],-c);lift=vector_sum(lift,lifts[p],c)
            if not v:
                reconstructed={}
                for i,c in lift.items():reconstructed=vector_sum(reconstructed,{j:a for j,a in enumerate(M[i])if a},c)
                assert reconstructed=={j:a for j,a in enumerate(row)if a}
            identities.append(dict(spatial_monomial=mon,member=not v,original350_row_functional={str(i):str(c)for i,c in lift.items()},nonmembership_remainder={str(j):str(c)for j,c in v.items()}))
        all_in=all(r['member']for r in identities)
        results.append(dict(rho=str(value),whole_P1_kernel_dimension=len(K[0]),whole_P2_dimension=9,parameters=n,parameter_monomials=[dict(linear=i)for i in range(n)]+[dict(quadratic_pair=p)for p in pairs],
                            full150_P0=list(map(str,q.vec(P0))),full_P1_kernel_tensor_basis=[list(map(str,q.vec(P)))for P in family[:len(K[0])]],
                            full350_Jacobi_coefficient_matrix=[[str(v)for v in row]for row in M],Jacobi_polynomial_rowspace_rank=len(basis),
                            all_Pf_mod_w_spatial_monomials=dmon,full_Pf_mod_w_parameter_coefficient_matrix=[[str(v)for v in row]for row in target],
                            all_nonzero_actual_Pf_mod_w_membership_tests=identities,whole_Pf_divisible_by_w_from_constant_Jacobi_row_combinations=all_in,
                            scope='Success proves wholeoriginal Poisson source has w factor in actualD; failure only excludes this degree0-multiplier ideal-membership criterion.'))
        Path('/tmp/plane_rank_zero_rank_two_high_m_flag_progress.json').write_text(json.dumps(dict(cases=results,seconds=time.monotonic()-begun),indent=2)+'\n')
    inputs=[Path(__file__),path,ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_RANK_ZERO_RANK_TWO_HIGH_MULTIPLICITY_FIRST_JET.md']
    result=dict(task_id='G',all_arithmetic_passed=True,cases=results,all_three_whole_flag_exclusions_pass=all(r['whole_Pf_divisible_by_w_from_constant_Jacobi_row_combinations']for r in results),
                unchanged_original_function_AST_names=sorted(frozen),original_functions_AST_unchanged=True,old_producer_bytes_preserved_in_input_report=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_flag_exclusions=result['all_three_whole_flag_exclusions_pass'],cases=[dict(rho=r['rho'],Jacobi_rank=r['Jacobi_polynomial_rowspace_rank'],all_members=r['whole_Pf_divisible_by_w_from_constant_Jacobi_row_combinations'],failed=sum(not t['member']for t in r['all_nonzero_actual_Pf_mod_w_membership_tests']))for r in results],seconds=result['seconds'])))


if __name__=='__main__':
    import sys
    pfaffian_flag_main()if '--pfaffian-flag'in sys.argv else main()
