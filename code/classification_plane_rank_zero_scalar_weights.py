#!/usr/bin/env python3
"""One finite full30/9 operator batch for all ten fixed scalar symbols.

These are whole Jordan strata under actual transformations, not generic
parameter samples. All first-symbol dimensions are restored explicitly.
"""
import argparse
from fractions import Fraction as F
import hashlib,itertools as it,json,resource,signal,time
from pathlib import Path
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]


def plus(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in q.PAIRS})
def wedge(U,V):return q.antisymmetric({(i,j):q.add(q.mul(U[i],V[j]),q.scale(q.mul(U[j],V[i]),-1))for i,j in q.PAIRS})


def actual(P):
    out={}
    for missing in range(5):
        i,j,k,l=[v for v in range(5)if v!=missing]
        pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
        out=q.add(out,q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
    return out


def certificate(M,n):
    A=[[F(v)for v in row]for row in M];r=0;piv=[]
    for j in range(n):
        k=next((k for k in range(r,len(A))if A[k][j]),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r];c=A[r][j];A[r]=[v/c for v in A[r]]
        for k in range(len(A)):
            if k!=r and A[k][j]:
                c=A[k][j];A[k]=[v-c*w for v,w in zip(A[k],A[r])]
        piv.append(j);r+=1
    free=[j for j in range(n)if j not in piv];K=[[F()for j in free]for i in range(n)]
    for k,j in enumerate(free):
        K[j][k]=F(1)
        for i,p in enumerate(piv):K[p][k]=-A[i][j]
    assert all(not sum(row[i]*K[i][j]for i in range(n))for row in M for j in range(len(free)))
    _,rows=q.rank(q.transpose(M));det=q.determinant([[M[i][j]for j in piv]for i in rows]);assert det or not piv
    return dict(rank=len(piv),nullity=len(free),minor_rows=rows,minor_columns=piv,minor_determinant=str(det),complete_kernel=[[str(v)for v in row]for row in K],full_kernel_substitution_zero=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite tenwholefixed scalar operator120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    a,b,c,z,w=[q.atom(1,i)for i in range(5)];zero={};En=[zero,zero,zero,z,w]
    Ts={'N3':([b,c,zero,zero,zero],6),'N2':([b,zero,zero,zero,zero],4),
        'D':([a,b,q.scale(c,-2),zero,zero],5),'D_plus_N2':([q.add(a,b),b,q.scale(c,-2),zero,zero],7)}
    Ks={'zero':(zero,0),'nil':(q.scale(q.mul(w,w),-1),2),'ss':(q.scale(q.mul(z,w),-2),3)}
    first=[]
    for normal in [3,4]:
        for ij in it.combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for mon,j in [([3,3],4),([3,4],3),([3,4],4),([4,4],3)]:
        for i in range(3):first.append(q.antisymmetric({(i,j):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in it.combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    results=[]
    for tname,(T,tdim)in Ts.items():
        for sname,(K,sdim)in Ks.items():
            if tname in ['D','D_plus_N2']and sname=='ss':continue
            P0=plus(wedge(T,En),q.antisymmetric({(3,4):K}))
            assert all(not v for v in q.jacobi(P0).values())and not actual(P0)
            D1=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in first]);D2=q.transpose([q.jvec(q.jacobi(plus(P0,P)))for P in second])
            c1=certificate(D1,30);c2=certificate(D2,9);d=c1['nullity']+c2['nullity'];base=tdim+sdim;image=base+5+d
            results.append(dict(type=tname+'_'+sname,T_Jordan_type=tname,S_Jordan_type=sname,affine_T_stratum_dimension=tdim,affine_S_stratum_dimension=sdim,
                                restored_affine_first_symbol_dimension=base,full150_P0=list(map(str,q.vec(P0))),all350_P0_Jacobi_zero=True,actual_P0_Pfaffian_zero=True,
                                full350_by30_D1=[[str(v)for v in row]for row in D1],full350_by9_D2=[[str(v)for v in row]for row in D2],D1=c1,D2=c2,
                                higher_fiber_dimension_upper=d,plane_choice_dimension=6,projectivization_subtraction=1,whole_projective_image_upper=image,G24_closed=image<=24))
    assert len(results)==10
    paths=[Path(__file__),ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_RANK_ZERO_SCALAR_FIRST_SYMBOL.md',
           ROOT/'research/classification/REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md',ROOT/'research/classification/REDUCED_DEGENERACY_ZERO_SURFACE_SEPARATE_REVIEW.md']
    result=dict(task_id='G',all_arithmetic_passed=True,cases=results,all_ten_fixed_G24_closed=all(r['G24_closed']for r in results),
                all30_source_basis=[list(map(str,q.vec(P)))for P in first],all9_source_basis=[list(map(str,q.vec(P)))for P in second],full350_row_labels=[dict(quad=ij,monomial=mon)for ij,mon in q.ROWS],
                scope='Ten wholefixed Tnonzero scalar first-symbol Jordan strata only. Four oneparameter/one twoparameter regular sources, T0K0 and tracelessnonzero dependent globalV sources remain separate. No E or emptiness from symbol Pf0.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),
                seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_fixed_G24_closed=result['all_ten_fixed_G24_closed'],cases=[{k:r[k]for k in ['type','higher_fiber_dimension_upper','whole_projective_image_upper','G24_closed']}for r in results],seconds=result['seconds'])))


def one_parameter_main():
    """All four whole one-parameter scalar symbols, with every minor root."""
    import ast
    import sympy as s
    import classification_plane_rank_zero_rank_two_high_m as H
    ap=argparse.ArgumentParser(description=one_parameter_main.__doc__);ap.add_argument('--one-parameter',action='store_true');ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite fourscalar wholeparameter operator120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    fixedpath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_FIXED_WEIGHTS_20261107.json';fixed=json.loads(fixedpath.read_text())
    frozen={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(fixed['producer_source_bytes']).body if isinstance(v,ast.FunctionDef)}
    now={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(Path(__file__).read_text()).body if isinstance(v,ast.FunctionDef)}
    assert all(now[k]==v for k,v in frozen.items())
    def decode(vector):
        out={}
        for value,(ij,e)in zip(vector,q.COLS):
            if F(value):out.setdefault(ij,{})[e+(0,0,0)]=F(value)
        return q.antisymmetric(out)
    first=list(map(decode,fixed['all30_source_basis']));second=list(map(decode,fixed['all9_source_basis']))
    a,b,c,z,w=[q.atom(1,i)for i in range(5)];t=q.atom(1,5);zero={};En=[zero,zero,zero,z,w]
    regular=[a,q.mul(t,b),q.scale(q.mul(q.add(q.atom(1),t),c),-1),zero,zero]
    models=[('repeated_ss_Kss',[q.mul(t,a),q.mul(t,b),q.scale(q.mul(t,c),-2),zero,zero],q.scale(q.mul(z,w),-2),8,[F(0)]),
            ('mixed_Kss',[q.add(q.mul(t,a),b),q.mul(t,b),q.scale(q.mul(t,c),-2),zero,zero],q.scale(q.mul(z,w),-2),10,[F(0)]),
            ('regular_ss_Kzero',regular,zero,8,[F(1),F(-2),-F(1,2)]),
            ('regular_ss_Knil',regular,q.scale(q.mul(w,w),-1),10,[F(1),F(-2),-F(1,2)])]
    results=[]
    for name,T,K,dim,excluded in models:
        P0=plus(wedge(T,En),q.antisymmetric({(3,4):K}))
        assert all(not v for v in q.jacobi(P0).values())and not actual(P0)
        D1=H.operator(P0,first);D2=H.operator(P0,second);analysis=H.analyze(D1,D2,dim,True)
        for wall in analysis['all_exceptional_roots']:
            value=s.sympify(wall['parameter'])if wall['minimal_polynomial'].count('tau')==1 and s.Poly(s.sympify(wall['minimal_polynomial']),H.TAU).degree()==1 else None
            wall['outside_declared_Jordan_domain']=value is not None and any(value==s.Rational(v.numerator,v.denominator)for v in excluded)
        domain_pass=analysis['G24_closed']and all(wall['G24_closed']or wall['outside_declared_Jordan_domain']for wall in analysis['all_exceptional_roots'])
        results.append(dict(type=name,restored_affine_first_symbol_dimension=dim,excluded_parameter_values=list(map(str,excluded)),outside_domain_covered_by='Separate zero-trace or fixed repeated/nilpotent Jordan strata; arithmetic of excluded roots retained without claiming this source includes them.',
                            full150_P0=[str(H.coefficient(P0.get(ij,{}),e))for ij,e in q.COLS],P0_all350_Jacobi_and_actualPf_zero=True,
                            full350_by30_D1=H.encode(D1),full350_by9_D2=H.encode(D2),analysis=analysis,all_domain_strata_G24_closed=domain_pass))
    paths=[Path(__file__),fixedpath,ROOT/'code/classification_plane_rank_zero_rank_two_high_m.py',ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_RANK_ZERO_SCALAR_FIRST_SYMBOL.md']
    result=dict(task_id='G',all_arithmetic_passed=True,cases=results,all_four_domain_G24_closed=all(r['all_domain_strata_G24_closed']for r in results),
                all_original_fixed_functions_AST_unchanged=True,original_fixed_function_names=sorted(frozen),old_fixed_report_bytes_preserved=True,
                scope='All four complete Tnonzero scalar oneparameter Jordan sources with exact domain exclusions. Two-parameter regularSS/Kss and T0K0 remain outside; fixed ten strata handled separately.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_four_G24_closed=result['all_four_domain_G24_closed'],cases=[dict(type=r['type'],generic_d=r['analysis']['higher_fiber_upper'],generic_image=r['analysis']['whole_projective_image_upper'],max_nonexcluded_wall=max((w['whole_projective_image_upper']for w in r['analysis']['all_exceptional_roots']if not w['outside_declared_Jordan_domain']),default=None),passed=r['all_domain_strata_G24_closed'])for r in results],seconds=result['seconds'])))


def two_parameter_main():
    """Whole regularSS/H source via a finite exact line/root cover."""
    import ast
    import sympy as s
    import classification_plane_rank_zero_rank_two_high_m as H
    ap=argparse.ArgumentParser(description=two_parameter_main.__doc__);ap.add_argument('--two-parameter',action='store_true');ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite scalar2parameter wholeline/root cover120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    fixedpath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_FIXED_WEIGHTS_20261107.json';fixed=json.loads(fixedpath.read_text())
    onepath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_ONE_PARAMETER_20261107.json';one=json.loads(onepath.read_text())
    frozen={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(one['producer_source_bytes']).body if isinstance(v,ast.FunctionDef)}
    now={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(Path(__file__).read_text()).body if isinstance(v,ast.FunctionDef)};assert all(now[k]==v for k,v in frozen.items())
    def decode(vector):
        out={}
        for value,(ij,e)in zip(vector,q.COLS):
            if F(value):out.setdefault(ij,{})[e+(0,0,0)]=F(value)
        return q.antisymmetric(out)
    first=list(map(decode,fixed['all30_source_basis']));second=list(map(decode,fixed['all9_source_basis']))
    L,M=s.symbols('lambda0 lambda1');DOMAIN=s.QQ.poly_ring(L,M)
    def coeff(poly,spatial):
        value=0
        for e,c in poly.items():
            assert not e[7]
            if e[:5]==spatial:value+=s.Rational(c.numerator,c.denominator)*L**e[5]*M**e[6]
        return s.expand(value)
    def op(P0,bases):
        columns=[]
        for P in bases:
            J=q.jacobi(plus(P0,P));columns.append([coeff(J[quad],e)for quad,e in q.ROWS])
        return s.Matrix(columns).T
    a,b,c,z,w=[q.atom(1,i)for i in range(5)];l=q.atom(1,5);m=q.atom(1,6);zero={}
    T=[q.mul(l,a),q.mul(m,b),q.scale(q.mul(q.add(l,m),c),-1),zero,zero];En=[zero,zero,zero,z,w]
    P0=plus(wedge(T,En),q.antisymmetric({(3,4):q.scale(q.mul(z,w),-2)}))
    assert all(not v for v in q.jacobi(P0).values())and not actual(P0)
    D1,D2=op(P0,first),op(P0,second);c1,m1=H.certificate(D1,DOMAIN);c2,m2=H.certificate(D2,DOMAIN)
    delta=(L-M)*(2*L+M)*(L+2*M);product=s.factor(m1*m2*delta);factors=s.factor_list(product)[1]
    generic_d=c1['nullity']+c2['nullity'];generic=dict(D1=c1,D2=c2,restored_affine_symbol_dimension=11,higher_fiber_upper=generic_d,whole_image_upper=16+generic_d,G24_closed=16+generic_d<=24)
    lines=[]
    for factor,multiplicity in factors:
        poly=s.Poly(factor,L,M)
        if poly.total_degree()!=1:
            checkpoint=dict(status='nonlinear_factor_stop',exact_factor=str(factor),whole_exception_polynomial=str(product),source_bytes=Path(__file__).read_text())
            args.output.with_name(args.output.stem+'_RUNTIME.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
            raise ValueError('Nonlinear exceptional curve: preserve and stop')
        A,B,C=poly.coeff_monomial(L),poly.coeff_monomial(M),poly.coeff_monomial(1)
        substitution={L:H.TAU,M:-(A*H.TAU+C)/B}if B else{L:-C/A,M:H.TAU}
        ex=s.factor(delta.subs(substitution,simultaneous=True))
        line=dict(equation=str(factor),multiplicity=int(multiplicity),parameterization={str(k):str(v)for k,v in substitution.items()},outside_Jordan_domain_entire_line=ex==0)
        if ex!=0:
            E1=D1.subs(substitution,simultaneous=True);E2=D2.subs(substitution,simultaneous=True);analysis=H.analyze(E1,E2,10,True)
            for wall in analysis['all_exceptional_roots']:
                f=s.Poly(s.sympify(wall['minimal_polynomial']),H.TAU)
                wall['outside_Jordan_domain']=s.rem(s.Poly(ex,H.TAU),f).is_zero
            line['analysis']=analysis;line['all_domain_points_G24_closed']=analysis['G24_closed']and all(wall['G24_closed']or wall['outside_Jordan_domain']for wall in analysis['all_exceptional_roots'])
        else:line['separate_coverage']='All repeated-eigenvalue scalar models in complete fixed/oneparameter packets, plus T0/H handproof. No arithmetic rank assertion on this excluded line.'
        lines.append(line)
        Path('/tmp/plane_scalar_two_parameter_progress.json').write_text(json.dumps(dict(lines=lines,seconds=time.monotonic()-begun),indent=2)+'\n')
    paths=[Path(__file__),fixedpath,onepath,ROOT/'code/classification_plane_rank_zero_rank_two_high_m.py',ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_RANK_ZERO_SCALAR_FIRST_SYMBOL.md']
    result=dict(task_id='G',all_arithmetic_passed=True,all350_P0_Jacobi_and_actualPf_zero=True,full150_P0=[str(coeff(P0.get(ij,{}),e))for ij,e in q.COLS],full350_by30_D1=H.encode(D1),full350_by9_D2=H.encode(D2),generic=generic,
                declared_regular_Jordan_domain=str(delta),whole_exception_polynomial=str(product),all_irreducible_linear_factors=[[str(f),int(m)]for f,m in factors],all_lines=lines,
                all_regularSS_H_G24_closed=generic['G24_closed']and all(line.get('all_domain_points_G24_closed',True)for line in lines),
                whole_cover_logic='FullC2 outside selected minors and Jordan delta is generic certificate; every remaining regular point lies on a listed nonexcluded line; each line has fullQQ(tau) certificate away its own selectedminor and full exact certificates at every root. Thus all line intersections and special points are retained.',
                all_old_functions_AST_unchanged=True,old_function_names=sorted(frozen),old_report_bytes_preserved=True,
                scope='Whole two-eigenparameter regular-semisimple T and S=H only; excluded repeated/Jordan/T0 strata separately retained. No E, component count or symbolfactor inference.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    real=[r for r in lines if not r['outside_Jordan_domain_entire_line']]
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_regularSS_H_G24_closed'],generic_image=generic['whole_image_upper'],all_lines=len(lines),nonexcluded_lines=len(real),all_line_root_certificates=sum(len(r['analysis']['all_exceptional_roots'])for r in real),max_line_image=max(r['analysis']['whole_projective_image_upper']for r in real),max_nonexcluded_point_image=max(w['whole_projective_image_upper']for r in real for w in r['analysis']['all_exceptional_roots']if not w['outside_Jordan_domain']),seconds=result['seconds'])))


if __name__=='__main__':
    import sys
    two_parameter_main()if '--two-parameter'in sys.argv else one_parameter_main()if '--one-parameter'in sys.argv else main()
