#!/usr/bin/env python3
"""Geometrically necessary abelian-conormal reduction of whole cubic sources.

Rankzero+irreducible quintic+nonplane surface forces multiplicity2 and a
Frobenius first jet, hence {I,I} in symbolic I^2. All first partials of all
three generator brackets are tested on the dominant surface parametrization.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
from fractions import Fraction as Q
import hashlib,json,resource,signal,time,os,subprocess,tempfile,sys
from functools import lru_cache
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_projected_veronese_zero_review import plus,scale,sparse_rref,monomials
ROOT=Path(__file__).resolve().parents[1];PAIRS=list(combinations(range(5),2))


def compute(case,old):
    start=time.monotonic();name=case['name'];x=s.symbols('x0:5');ss,tt,u,v,c=s.symbols('ss tt u v c')
    if name=='smooth_scroll_12':param=(ss,tt,u,v);mapping=(ss*u,tt*u,ss*ss*v,ss*tt*v,tt*tt*v)
    else:param=(ss,tt,c);mapping=(ss**3,ss*ss*tt,ss*tt*tt,tt**3,c)
    env=dict(zip(map(str,x),x));quadrics=[s.sympify(t,locals=env)for t in case['quadratic_ideal_generators']];sub=dict(zip(x,mapping))
    assert all(s.expand(q.subs(sub,simultaneous=True))==0 for q in quadrics)
    derivatives=[[s.expand(s.diff(q,xx).subs(sub,simultaneous=True))for xx in x]for q in quadrics]
    rowlabels=[];matrix=[]
    for qa,qb in combinations(range(3),2):
        for normal_derivative in range(5):
            columns=[];exponents=set()
            for qg in range(3):
                for i,j in PAIRS:
                    polynomial=s.Poly(s.expand(derivatives[qg][normal_derivative]*(derivatives[qa][i]*derivatives[qb][j]-derivatives[qa][j]*derivatives[qb][i])),*param,domain=s.QQ)
                    terms={e:co for e,co in polynomial.terms()if co};columns.append(terms);exponents.update(terms)
            for e in sorted(exponents,reverse=True):
                rowlabels.append({'generator_pair':[qa,qb],'ambient_partial':normal_derivative,'parametrization_exponents':e});matrix.append([col.get(e,s.Integer(0))for col in columns])
    M=s.Matrix(matrix);dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    red,piv=dm(M).rref();R=red.to_Matrix();piv=list(piv);free=[j for j in range(30)if j not in piv];B=s.zeros(30,len(free))
    for k,j in enumerate(free):
        B[j,k]=1
        for row,p in enumerate(piv):B[p,k]=-R[row,j]
    assert M*B==s.zeros(M.rows,B.cols)
    selected=list(dm(M.T).rref()[1]);minor=M[selected,piv];determinant=minor.det();assert determinant!=0 and len(piv)+B.cols==30
    d=B.cols;newcols=list(combinations_with_replacement(range(d),2));newci={a:i for i,a in enumerate(newcols)}
    oldcols=[tuple(a)for a in old['all465_quadratic_parameter_monomials']]
    def restrict(row):
        out={}
        for k,value in row.items():
            i,j=oldcols[int(k)];value=Q(value)
            for r in range(d):
                if not B[i,r]:continue
                for t in range(d):
                    if not B[j,t]:continue
                    col=newci[tuple(sorted((r,t)))];out[col]=out.get(col,Q(0))+value*Q(B[i,r])*Q(B[j,t])
        return {j:c for j,c in out.items()if c}
    J=[restrict(row)for row in case['all350_Jacobi_rows']];Drows=[restrict(row)for row in case['all126_actual_D_rows']]
    basis,proof,chosen=sparse_rref(J,True);targetcert=[];allmember=True
    for k,target in enumerate(Drows):
        residual=dict(target);functional={}
        for p in sorted(basis):
            if p in residual:
                coefficient=residual[p];residual=plus(residual,scale(basis[p],-coefficient));functional=plus(functional,scale(proof[p],coefficient))
        allmember &=not residual
        assert plus(plus(*(scale(J[i],a)for i,a in functional.items())),residual)==target
        if target:targetcert.append({'original126_D_row':k,'full350_Jacobi_functional':functional,'residual':residual})
    answer={'name':name,'parametrization_variables':list(map(str,param)),'dominant_surface_parametrization':list(map(str,mapping)),
        'all_three_ideal_generators':case['quadratic_ideal_generators'],'all_generator_pair_all_firstpartial_equations':rowlabels,
        'full_linear_matrix':[[str(c)for c in row]for row in matrix],'linear_shape':list(M.shape),'linear_rank':len(piv),
        'rank_minor':{'rows':selected,'columns':piv,'determinant':str(determinant)},'free_original30_indices':free,
        'complete_abelian_conormal_source_dimension':d,'complete_original30_by_reduced_source_matrix':[[str(a)for a in row]for row in B.tolist()],
        'all_reduced_source_quadratic_monomials':newcols,'all350_reduced_Jacobi_rows':J,'all126_reduced_actual_D_rows':Drows,
        'reduced_Jacobi_row_rank':len(basis),'all_actual_D_rows_in_reduced_full350_Jacobi_rowspace':allmember,'actual_D_certificates':targetcert,
        'seconds':time.monotonic()-start,
        'scope':'Complete necessary abelian-conormal source for rankzero and irreducible quintic on this nonplane surface. Positive actualD membership proves this entire geometric branch empty. It does not impose this source condition on arbitrary rank2 restrictions or on plane multiplicity>=3.'}
    print(name,'full linear',M.shape,'rank',len(piv),'source',d,'all350 rank',len(basis),'actualD membership',allmember,flush=True)
    return answer


def main():
    destination=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    oldpath=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_JACOBI_20261107.json';old=json.loads(oldpath.read_text())
    cases=[compute(case,old)for case in old['cases']]
    result={'task_id':'G','all_passed':True,'input_report':str(oldpath.relative_to(ROOT)),'input_sha256':hashlib.sha256(oldpath.read_bytes()).hexdigest(),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'cases':cases,
        'both_rankzero_irreducible_quintic_cubic_surface_branches_empty':all(c['all_actual_D_rows_in_reduced_full350_Jacobi_rowspace']for c in cases),
        'seconds':time.monotonic()-start,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    destination.write_text(json.dumps(result,indent=2,default=str)+'\n')


def dimension(localized=False):
    """A dimension upper bound suffices; retain exact lifted leading monomials."""
    destination=ROOT/('reports/CLASSIFICATION_CUBIC_CONE_ABELIAN_RANK4_DIMENSION_20261107.json'if localized else'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_DIMENSION_20261107.json')
    if destination.exists():raise FileExistsError(destination)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));signal.alarm(165)
    path=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json';source=json.loads(path.read_text())
    answer={'task_id':'G','input_report':str(path.relative_to(ROOT)),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'cases':[]}
    prefix=ROOT/'work/nonsquare-v62/singular-runtime/prefix';env_run=dict(os.environ)
    env_run['LD_LIBRARY_PATH']='/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu:'+str(prefix/'usr/lib/x86_64-linux-gnu')+':'+env_run.get('LD_LIBRARY_PATH','')
    env_run['SINGULAR_PROCS_DIR']=str(prefix/'usr/lib/x86_64-linux-gnu/singular/MOD');env_run['SINGULARPATH']='/workspace/cloud-setup/prefix/usr/share/singular/LIB'
    for case in source['cases']:
        if localized and case['name']!='twisted_cubic_cone_03':continue
        start=time.monotonic();d=case['complete_abelian_conormal_source_dimension'];parameters=s.symbols('r0:'+str(d));z=(*parameters,s.Symbol('inverse'))if localized else parameters;env={str(v):v for v in z};cols=case['all_reduced_source_quadratic_monomials']
        rows=[{int(j):Q(c)for j,c in row.items()}for row in case['all350_reduced_Jacobi_rows']];_,_,selected=sparse_rref(rows)
        gens=[s.Poly(sum(s.Rational(c)*z[cols[j][0]]*z[cols[j][1]]for j,c in rows[i].items()),*z,domain=s.QQ).clear_denoms()[1].primitive()[1]for i in selected]
        packet={'name':case['name'],'all_reduced_source_parameters_retained':d,'original350_row_indices':selected,'generators':[str(g.as_expr())for g in gens],
            'target_fixed_affine_dimension':7 if case['name']=='smooth_scroll_12'else 9}
        if localized:
            point=[1,0,0,1,0];Dpoint=s.expand(sum(s.prod(s.Integer(a)**e for a,e in zip(point,exponent))*sum(s.Rational(c)*parameters[cols[int(j)][0]]*parameters[cols[int(j)][1]]for j,c in row.items())for exponent,row in zip(monomials(5,5),case['all126_reduced_actual_D_rows'])))
            assert Dpoint!=0;gens.append(s.Poly(z[-1]*Dpoint-1,*z,domain=s.QQ));packet['actual_D_at_dense_orbit_point']=str(Dpoint);packet['rank4_localization_point']=point;packet['generators']=[str(g.as_expr())for g in gens]
        script='ring R=0,('+','.join(map(str,z))+'),dp;\nshort=0;\nideal I='+','.join(str(g.as_expr()).replace('**','^')for g in gens)+';\nmatrix T;ideal G=liftstd(I,T);\n'
        output_limit='min(size(G),200)'if localized else'size(G)'
        script+='print("DIM:"+string(dim(G)));\nint i;int j;int count='+output_limit+';for(i=1;i<=count;i++){print("G:"+string(i)+":"+string(G[i]));}\n'
        script+='for(i=1;i<=nrows(T);i++){for(j=1;j<=count;j++){if(T[i,j]!=0){print("L:"+string(i)+":"+string(j)+":"+string(T[i,j]));}}}\nquit;\n';packet['exact_native_input']=script
        try:
            with tempfile.TemporaryDirectory(prefix='cubic-abelian-dim-')as temp:
                task=Path(temp)/'source.sing';task.write_text(script)
                run=subprocess.run([str(prefix/'usr/bin/Singular'),'-q',str(task)],env=env_run,capture_output=True,text=True,timeout=60,check=False)
            assert run.returncode==0 and '?'not in run.stdout,run.stdout[-20000:]+run.stderr
            G={};lifts=[];native_dim=None
            for line in run.stdout.splitlines():
                if line.startswith('G:'):
                    _,j,p=line.split(':',2);G[int(j)-1]=s.Poly(s.sympify(p.replace('^','**'),locals=env),*z,domain=s.QQ)
                elif line.startswith('L:'):
                    _,i,j,p=line.split(':',3);lifts.append((int(i)-1,int(j)-1,s.Poly(s.sympify(p.replace('^','**'),locals=env),*z,domain=s.QQ)))
                elif line.startswith('DIM:'):native_dim=int(line.split(':',1)[1])
            assert G and native_dim is not None
            sums={j:s.Poly(0,*z,domain=s.QQ)for j in G}
            for i,j,p in lifts:sums[j]+=gens[i]*p
            assert sums==G
            leading=[G[j].monoms(order='grevlex')[0]for j in sorted(G)];edges=tuple(sorted(set(sum(1<<i for i,e in enumerate(mon)if e)for mon in leading)))
            @lru_cache(None)
            def cover(edges):
                if not edges:return 0
                if 0 in edges:return len(z)+1
                edge=min(edges,key=int.bit_count)
                return 1+min(cover(tuple(e for e in edges if not e&(1<<i)))for i in range(len(z))if edge&(1<<i))
            height=cover(edges);bound=len(z)-height
            if not localized:assert bound==native_dim
            else:assert bound>=native_dim
            groupdim=6 if case['name']=='smooth_scroll_12'else 8;imagebound=24+bound-1-groupdim
            packet.update(status='exact_lifted_leading_dimension',all_ideal_lifts_exactly_checked=True,
                lifted_polynomials=[str(G[j].as_expr())for j in sorted(G)],ideal_membership_lifts=[{'generator':i,'polynomial':j,'coefficient':str(p.as_expr())}for i,j,p in lifts],
                leading_monomial_exponents=leading,leading_support_minimum_vertex_cover=height,fixed_Poisson_source_dimension_upper=bound,
                source_preserving_group_dimension=groupdim,whole_projective_PGL_image_dimension_upper=imagebound,excludes_unknown_G25_to27=imagebound<25,native_stdout=run.stdout)
        except subprocess.TimeoutExpired as error:
            packet.update(status='finite60second_native_limit_no_dimension_proof',partial_stdout=str(error.stdout or '')[-20000:],partial_stderr=str(error.stderr or '')[-20000:])
        packet['seconds']=time.monotonic()-start;answer['cases'].append(packet);destination.write_text(json.dumps(answer,indent=2)+'\n')
        print(case['name'],packet['status'],packet.get('fixed_Poisson_source_dimension_upper'),packet.get('whole_projective_PGL_image_dimension_upper'),packet['seconds'],flush=True)


if __name__=='__main__':
    if '--cone-localized'in sys.argv:dimension(True)
    elif '--dimension'in sys.argv:dimension()
    else:main()
