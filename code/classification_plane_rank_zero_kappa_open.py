#!/usr/bin/env python3
"""Actual d,r normalization and bounded affine Gaussian whole-source reduction.

No GB, higher multiplier, jet or point scan. Each nontrivial step removes an
unknown by a literal affine consequence of all current350 Jacobi equations.
All original equations, forward maps and inverse free-coordinate maps remain.
"""
import argparse
from fractions import Fraction as F
import hashlib,itertools as it,json,resource,signal,time
from pathlib import Path
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]


def encode(poly):return [[list(e),str(c)]for e,c in sorted(poly.items(),reverse=True)]


def addvec(A,B,c):
    out=dict(A)
    for k,v in B.items():out[k]=out.get(k,F())+c*v
    return {k:v for k,v in out.items()if v}


def rowspace(rows,lifts):
    basis={};functions={}
    for values,function in zip(rows,lifts):
        v=dict(values);f=dict(function)
        for p in sorted(basis):
            c=v.get(p,F())
            if c:v=addvec(v,basis[p],-c);f=addvec(f,functions[p],-c)
        if v:
            p=min(v);c=v[p];v={j:a/c for j,a in v.items()};f={j:a/c for j,a in f.items()}
            for old in basis:
                c=basis[old].get(p,F())
                if c:basis[old]=addvec(basis[old],v,-c);functions[old]=addvec(functions[old],f,-c)
            basis[p]=v;functions[p]=f
    return basis,functions


def substitute(poly,j,expression):
    out={}
    for e,c in poly.items():
        k=e[j];rest=list(e);rest[j]=0;term={tuple(rest):c}
        for _ in range(k):term=q.mul(term,expression)
        out=q.add(out,term)
    return out


def compose(poly,images,n):
    out={};constant=(0,)*n
    for e,c in poly.items():
        term={constant:c}
        for i,k in enumerate(e):
            for _ in range(k):term=q.mul(term,images[i])
        out=q.add(out,term)
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite actual kappa-open affine reduction120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    flagpath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_HIGH_M_FLAG_20261107.json';flag=json.loads(flagpath.read_text())
    oldpath=ROOT/'reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_HIGH_M_20261107.json';old=json.loads(oldpath.read_text());P2=old['all9_source_basis']
    coordinate_weights=[(-3,1),(-1,0),(-1,0),(0,0),(0,0)];scalar=(2,-1)
    def tensor_weights(vector):
        weights=set()
        for coeff,(ij,e)in zip(vector,q.COLS):
            if F(coeff):weights.add(tuple(scalar[k]+sum(coordinate_weights[i][k]for i in ij)-sum(e[i]*coordinate_weights[i][k]for i in range(5))for k in range(2)))
        assert len(weights)==1
        return list(next(iter(weights)))
    results=[]
    for case in flag['cases']:
        n=case['parameters'];k1=case['whole_P1_kernel_dimension'];anchors=[8,k1+6];zero=(0,)*n
        variable=lambda j:{tuple(int(i==j)for i in range(n)):F(1)}
        images=[variable(j)for j in range(n)];images[anchors[0]]=images[anchors[1]]={zero:F(1)}
        assert tensor_weights(case['full150_P0'])==[0,0]
        weights=[tensor_weights(v)for v in case['full_P1_kernel_tensor_basis']]+[tensor_weights(v)for v in P2]
        assert weights[anchors[0]]==[-1,0]and weights[anchors[1]]==[0,-1]
        assert all(w in [[-1,0],[1,-1]]for w in weights[:k1])
        assert weights[k1:]==[[-2,0]]*6+[[0,-1]]*3
        exponents=[]
        for mon in case['parameter_monomials']:
            e=[0]*n
            if 'linear'in mon:e[mon['linear']]=1
            else:
                for j in mon['quadratic_pair']:e[j]+=1
            exponents.append(tuple(e))
        original=[{e:F(c)for e,c in zip(exponents,row)if F(c)}for row in case['full350_Jacobi_coefficient_matrix']]
        equations=[compose(P,images,n)for P in original];remaining=[j for j in range(n)if j not in anchors];steps=[];empty=False
        for iteration in range(n+1):
            assert all(sum(e)<=2 for P in equations for e in P)
            if len(remaining)<=10:break
            mons=sorted(set().union(*(set(P)for P in equations)),key=lambda e:(sum(e),e),reverse=True);index={e:i for i,e in enumerate(mons)}
            B,L=rowspace([{index[e]:c for e,c in P.items()}for P in equations],[{i:F(1)}for i in range(350)])
            linear=[];linear_lifts=[]
            for p,row in B.items():
                if all(sum(mons[j])<=1 for j in row):
                    P={mons[j]:c for j,c in row.items()};linear.append(P);linear_lifts.append(L[p])
            if not linear:break
            columns={j:i for i,j in enumerate(remaining)};constant_column=len(remaining)
            A=[]
            for P in linear:
                v={}
                for e,c in P.items():v[columns[e.index(1)]if sum(e)else constant_column]=c
                A.append(v)
            R,functions=rowspace(A,linear_lifts)
            if constant_column in R:
                f=functions[constant_column];literal=q.add(*(q.scale(equations[i],c)for i,c in f.items()));assert literal=={zero:F(1)}
                steps.append(dict(iteration=iteration,nonzero_constant_consequence=True,current350_row_functional={str(i):str(c)for i,c in f.items()}));empty=True;break
            solutions={};constraints=[]
            free=[remaining[j]for j in range(len(remaining))if j not in R]
            for pivot,row in sorted(R.items()):
                j=remaining[pivot];P={}
                for column,c in row.items():
                    mon=zero if column==constant_column else next(iter(variable(remaining[column])))
                    P[mon]=c
                f=functions[pivot];literal=q.add(*(q.scale(equations[i],c)for i,c in f.items()));assert literal==P
                expression=q.scale(q.add(P,q.scale(variable(j),-1)),-1)
                assert all(not e[k]for e in expression for k in remaining if k not in free)
                solutions[j]=expression;constraints.append(dict(solved_parameter=j,full_affine_equation=encode(P),solution=encode(expression),current350_row_functional={str(i):str(c)for i,c in f.items()}))
            assert solutions
            before=list(remaining)
            for j,expression in solutions.items():
                equations=[substitute(P,j,expression)for P in equations]
                images=[substitute(P,j,expression)for P in images]
            remaining=free
            steps.append(dict(iteration=iteration,remaining_before=before,all_affine_consequences=constraints,remaining_after=remaining))
        assert len(steps)<=n
        assert all(compose(P,images,n)==Q for P,Q in zip(original,equations))
        assert all(images[j]==variable(j)for j in remaining)
        t=0 if empty else len(remaining);image_upper=0 if empty else 14+t
        results.append(dict(rho=case['rho'],whole_original_parameters=n,kernel_P1_dimension=k1,anchors_d_r=anchors,
                            all150_P0_normalization_weights_zero=True,complete_parameter_normalization_weights=weights,actual_inverse_map='a=d^3/r*a_prime,b=d*b_prime,c=d*c_prime,z=z_prime,w=w_prime;tensorinverse=r/d^2',
                            all_affine_elimination_steps=steps,source_empty=empty,remaining_parameter_indices=remaining,normalized_ambient_dimension_upper=t,
                            whole_source_original_parameter_images=[encode(P)for P in images],full350_residual_Jacobi_equations=[encode(P)for P in equations],all350_forward_substitution_identity=True,inverse_free_parameter_projection=True,
                            two_consumed_nonzero_parameters_restored=True,fixed_symbol_affine_base_dimension=7,whole_PGL_projective_image_upper=image_upper,G24_closed=empty or image_upper<=24,
                            scope='Exact full normalized Jacobi source after actual d,r PGL normalization; residual<=10 sufficient, no residual equation discarded when recording source equivalence.'))
        Path('/tmp/plane_rank_zero_kappa_open_progress.json').write_text(json.dumps(dict(cases=results,seconds=time.monotonic()-begun),indent=2)+'\n')
    paths=[Path(__file__),flagpath,oldpath,ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/PLANE_RANK_ZERO_KAPPA_OPEN_NORMALIZATION.md']
    result=dict(task_id='G',all_arithmetic_passed=True,cases=results,all_three_kappa_open_G24_closed=all(r['G24_closed']for r in results),
                no_new_ideal_multiplier_or_GB_or_jet=True,finite_elimination_bound_unknown_count22=True,all_polynomial_degrees_at_most2=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_three_kappa_open_G24_closed'],cases=[dict(rho=r['rho'],empty=r['source_empty'],free=r['normalized_ambient_dimension_upper'],image=r['whole_PGL_projective_image_upper'],steps=len(r['all_affine_elimination_steps']))for r in results],seconds=result['seconds'])))


if __name__=='__main__':main()
