#!/usr/bin/env python3
"""One geometrically exhaustive rank4 patch on each cubic-surface orbit.

The source-preserving group has a dense open point orbit. Thus every nonzero
actual quintic can be moved to D(p0)!=0. Exact lifted leading monomials of
J+(inverse*D(p0)-1) suffice; no assertion requires a reduced Groebner basis.
"""
from pathlib import Path
from functools import lru_cache
from fractions import Fraction
import hashlib,json,os,resource,signal,subprocess,tempfile,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]


def one(case,source):
    start=time.monotonic();name=case['name'];z=s.symbols('z0:30');inv=s.Symbol('inverse');variables=(*z,inv);env={str(v):v for v in variables}
    columns=source['all465_quadratic_parameter_monomials'];monomials=[z[i]*z[j]for i,j in columns]
    selected=case['selected_independent_original350_rows']
    original=[s.Poly(sum(s.Rational(c)*monomials[int(j)]for j,c in case['all350_Jacobi_rows'][i].items()),*variables,domain=s.QQ).clear_denoms()[1].primitive()[1]for i in selected]
    point=[0,0,0,1,0]if name=='smooth_scroll_12'else[1,0,0,1,0]
    def mon5():
        from classification_projected_veronese_zero_review import monomials
        return monomials(5,5)
    D=s.expand(sum(s.prod(s.Integer(a)**e for a,e in zip(point,powers))*sum(s.Rational(c)*monomials[int(j)]for j,c in row.items())for powers,row in zip(mon5(),case['all126_actual_D_rows'])))
    assert D!=0
    gens=original+[s.Poly(inv*D-1,*variables,domain=s.QQ)]
    script='ring R=0,('+','.join(map(str,variables))+'),dp;\nshort=0;\n'
    script+='ideal I='+','.join(str(g.as_expr()).replace('**','^')for g in gens)+';\nmatrix T;ideal G=liftstd(I,T);\n'
    script+='print("DIM:"+string(dim(G)));\nint i;int j;for(i=1;i<=size(G);i++){print("G:"+string(i)+":"+string(G[i]));}\n'
    script+='for(i=1;i<=nrows(T);i++){for(j=1;j<=ncols(T);j++){if(T[i,j]!=0){print("L:"+string(i)+":"+string(j)+":"+string(T[i,j]));}}}\nquit;\n'
    result={'name':name,'source_point':point,'actual_D_at_point':str(D),'all30_source_parameters_retained':True,'selected_independent_original350_rows':selected,
        'localized_ideal_generators':[str(g.as_expr())for g in gens],'exact_native_input':script,'source_preserving_group_dimension':6 if name=='smooth_scroll_12'else 8,
        'fixed_localized_source_dimension_target':7 if name=='smooth_scroll_12'else 9}
    prefix=ROOT/'work/nonsquare-v62/singular-runtime/prefix';nativeenv=dict(os.environ)
    nativeenv['LD_LIBRARY_PATH']='/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu:'+str(prefix/'usr/lib/x86_64-linux-gnu')+':'+nativeenv.get('LD_LIBRARY_PATH','')
    nativeenv['SINGULAR_PROCS_DIR']=str(prefix/'usr/lib/x86_64-linux-gnu/singular/MOD');nativeenv['SINGULARPATH']='/workspace/cloud-setup/prefix/usr/share/singular/LIB'
    try:
        with tempfile.TemporaryDirectory(prefix='cubic-scroll-rank4-')as temp:
            task=Path(temp)/'source.sing';task.write_text(script)
            run=subprocess.run([str(prefix/'usr/bin/Singular'),'-q',str(task)],env=nativeenv,capture_output=True,text=True,timeout=60,check=False)
        if run.returncode or '?'in run.stdout:
            result.update(status='native_runtime_failure_no_proof',returncode=run.returncode,partial_stdout=run.stdout[-20000:],stderr=run.stderr[-20000:])
        else:
            lifted={};lifts=[];native_dim=None
            for line in run.stdout.splitlines():
                if line.startswith('G:'):
                    _,j,p=line.split(':',2);lifted[int(j)-1]=s.Poly(s.sympify(p.replace('^','**'),locals=env),*variables,domain=s.QQ)
                elif line.startswith('L:'):
                    _,i,j,p=line.split(':',3);lifts.append((int(i)-1,int(j)-1,s.Poly(s.sympify(p.replace('^','**'),locals=env),*variables,domain=s.QQ)))
                elif line.startswith('DIM:'):native_dim=int(line.split(':',1)[1])
            assert lifted and native_dim is not None
            sums={j:s.Poly(0,*variables,domain=s.QQ)for j in lifted}
            for i,j,p in lifts:sums[j]+=gens[i]*p
            assert lifted==sums
            leading=[lifted[j].monoms(order='grevlex')[0]for j in sorted(lifted)]
            edges=tuple(sorted(set(sum(1<<i for i,e in enumerate(mon)if e)for mon in leading)))
            def minimal(edges):return tuple(e for e in edges if not any(f!=e and f&e==f for f in edges))
            @lru_cache(None)
            def cover(edges):
                if not edges:return 0
                if 0 in edges:return 32
                edge=min(edges,key=int.bit_count)
                return 1+min(cover(minimal(tuple(e for e in edges if not e&(1<<i))))for i in range(31)if edge&(1<<i))
            height=cover(minimal(edges));bound=-1 if height==32 else 31-height
            assert bound==native_dim
            groupdim=result['source_preserving_group_dimension'];imagebound=-1 if bound<0 else 24+bound-1-groupdim
            result.update(status='exact_lifted_dimension_certificate',all_exact_lifts_verified=True,native_dimension=native_dim,
                lifted_polynomials=[str(lifted[j].as_expr())for j in sorted(lifted)],exact_ideal_membership_lifts=[{'generator':i,'polynomial':j,'coefficient':str(p.as_expr())}for i,j,p in lifts],
                leading_monomial_exponents=leading,leading_support_minimum_vertex_cover=height,fixed_localized_source_dimension_upper=bound,
                whole_projective_PGL_rank4_image_dimension_upper=imagebound,excludes_unknown_G25_to27=imagebound<25,native_stdout=run.stdout,stderr=run.stderr)
    except subprocess.TimeoutExpired as error:
        result.update(status='finite60second_native_limit_no_dimension_proof',partial_stdout=str(error.stdout or '')[-20000:],partial_stderr=str(error.stderr or '')[-20000:])
    result['seconds']=time.monotonic()-start
    return result


def main():
    destination=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_RANK4_DIMENSION_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));signal.alarm(170);start=time.monotonic()
    sourcepath=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_JACOBI_20261107.json';source=json.loads(sourcepath.read_text())
    result={'task_id':'G','input_report':str(sourcepath.relative_to(ROOT)),'input_sha256':hashlib.sha256(sourcepath.read_bytes()).hexdigest(),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'cases':[],
        'scope':'All rank4 fixed-surface sources meet the one D(p0)!=0 patch after the explicit dense-open source-preserving point orbit action. Only lifted leading-monomial dimension bounds are asserted. Rank2 source components are irrelevant to this localization; any finite timeout is no mathematical proof.'}
    for case in source['cases']:
        answer=one(case,source);result['cases'].append(answer);destination.write_text(json.dumps(result,indent=2)+'\n')
        print(answer['name'],answer['status'],answer.get('fixed_localized_source_dimension_upper'),answer.get('whole_projective_PGL_rank4_image_dimension_upper'),answer['seconds'],flush=True)
    result['seconds']=time.monotonic()-start;result['maximum_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;destination.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
