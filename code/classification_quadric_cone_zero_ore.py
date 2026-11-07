#!/usr/bin/env python3
"""Whole quadric-cone Ore reduction and complete one-variable rank bounds.

No GB or frame samples. Every exceptional root of the chosen generic minors
is retained and evaluated exactly. The source bridge includes Euler5.
"""
import argparse,hashlib,itertools as it,json,resource,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
PAIRS5=list(it.combinations(range(5),2));PAIRS3=list(it.combinations(range(3),2))


def rank(M):return len(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])
def lie(Y,B,u):return s.Matrix(3,3,lambda i,j:s.expand(sum(Y[k]*s.diff(B[i,j],u[k])-B[k,j]*s.diff(Y[i],u[k])-B[i,k]*s.diff(Y[j],u[k])for k in range(3))))
def jac(P,u):return {tri:s.expand(sum(P[i,l]*s.diff(P[j,k],u[l])+P[j,l]*s.diff(P[k,i],u[l])+P[k,l]*s.diff(P[i,j],u[l])for l in range(len(u))))for tri in it.combinations(range(len(u)),3)for i,j,k in [tri]}
def wedge(v,w):return v*w.T-w*v.T
def coeffs(expressions,u):return [c for f in expressions for e,c in s.Poly(s.expand(f),*u).terms()if c!=0]
def encode(P):return {f'{i},{j}':str(s.factor(P[i,j]))for i,j in it.combinations(range(P.rows),2)if P[i,j]!=0}


def eliminate_constant_linear(M,rhs,variables):
    basis={};compat=[]
    for row,value in zip(M.tolist(),rhs):
        row={i:s.Rational(c)for i,c in enumerate(row)if c};value=s.expand(value)
        for p in sorted(basis):
            if p in row:
                old,oldrhs=basis[p];a=row[p];value=s.expand(value-a*oldrhs)
                for j,c in old.items():row[j]=row.get(j,0)-a*c
                row={j:c for j,c in row.items()if c}
        if row:
            p=min(row);a=row[p];basis[p]=({j:c/a for j,c in row.items()},s.expand(value/a))
        elif value!=0:compat.append(value)
    for p in sorted(basis,reverse=True):
        row,value=basis[p]
        for q in sorted(basis):
            old,oldrhs=basis[q]
            if q<p and p in old:
                a=old[p];new={j:old.get(j,0)-a*row.get(j,0)for j in set(old)|set(row)};basis[q]=({j:c for j,c in new.items()if c},s.expand(oldrhs-a*value))
    free=[j for j in range(len(variables))if j not in basis]
    solution={variables[p]:s.expand(value-sum(c*variables[j]for j,c in row.items()if j!=p))for p,(row,value)in basis.items()}
    return solution,free,compat


def principal_block(Ybasis,B0,T1,u,c,degree):
    monomials=[s.prod(u[i]for i in terms)for terms in it.combinations_with_replacement(range(3),degree)]
    columns=[]
    for Y in Ybasis:
        F=lie(Y,B0,u)-wedge(Y,T1)
        columns.append(s.Matrix([s.Poly(F[i,j],*u).coeff_monomial(mon)for i,j in PAIRS3 for mon in monomials]))
    M=s.Matrix.hstack(*columns);field=s.QQ.frac_field(c);DM=DomainMatrix.from_Matrix(M).convert_to(field);RR,piv=DM.rref();R=RR.to_Matrix();rows=list(DomainMatrix.from_Matrix(M.T).convert_to(field).rref()[1]);piv=list(piv)
    minor=s.factor(M[rows,piv].det());assert minor!=0
    free=[j for j in range(M.cols)if j not in piv];K=s.zeros(M.cols,len(free))
    for k,j in enumerate(free):
        K[j,k]=1
        for i,p in enumerate(piv):K[p,k]=-R[i,j]
    assert(M*K).applyfunc(s.cancel)==s.zeros(M.rows,K.cols)
    factors=s.factor_list(s.Poly(minor,c,domain=s.QQ))[1];roots=[]
    for factor,multiplicity in factors:
        assert factor.degree()==1,'Nonlinear exceptional polynomial retained; no incomplete root claim permitted.'
        roots.append(-factor.nth(0)/factor.nth(1))
    return M,dict(shape=list(M.shape),all_matrix_entries=[[str(v)for v in row]for row in M.tolist()],generic_rank=len(piv),generic_minor_rows=rows,generic_minor_columns=piv,generic_minor=str(minor),generic_complete_kernel=[[str(v)for v in row]for row in K.tolist()],all_exceptional_roots=list(map(str,roots))),roots


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    sourcepath=ROOT/'reports/CLASSIFICATION_QUADRIC_CONE_ZERO_SOURCE_20261107.json';source=json.loads(sourcepath.read_text());u=s.Matrix(s.symbols('x0:3'));z,w=s.symbols('z w');x=s.Matrix([*u,z,w]);q=u[0]*u[2]-u[1]**2;c=s.Symbol('c');records=[]
    for case,ij,offset in zip(source['cases'],[(0,2),(0,1)],[18,17]):
        started=time.monotonic();C=s.zeros(3);C[ij]=1;C[ij[::-1]]=-1;X=C*s.Matrix([s.diff(q,t)for t in u]);A=s.zeros(5);A[:3,:3]=C;A[4,3]=1;A[3,4]=-1
        DD=s.symbols('D0:3');aa=s.Matrix(s.symbols('a0:3'));bb=s.Matrix(s.symbols('b0:3'));cc=s.Symbol('cc');MM=s.Matrix(3,3,s.symbols('M0:9'));dd=s.Matrix(s.symbols('d0:3'));ee=s.Matrix(s.symbols('e0:3'));parameters=(*DD,*aa,*bb,cc,*MM,*dd,*ee);assert len(parameters)==25
        D=s.zeros(3)
        for value,pair in zip(DD,PAIRS3):D[pair]=value;D[pair[::-1]]=-value
        L=s.zeros(5);L[:3,:3]=w*D
        for i in range(3):L[4,i]=-X[i]+w*aa[i];L[i,4]=-L[4,i];L[3,i]=(MM*u+z*dd+w*ee)[i];L[i,3]=-L[3,i]
        L[4,3]=(bb.dot(u)+cc*z);L[3,4]=-L[4,3]
        rawcoeff=s.Matrix([s.Poly(L[i,j],*x).coeff_monomial(t)for i,j in PAIRS5 for t in x]);point=rawcoeff.subs(dict.fromkeys(parameters,0));bridge=rawcoeff.jacobian(parameters)
        E=s.Matrix.hstack(*[s.Matrix(v)for v in [[int((j==k and l==i))-int((i==k and l==j))for i,j in PAIRS5 for l in range(5)]for k in range(5)]])
        full=bridge.row_join(E);oldK=s.Matrix([[s.Rational(v)for v in row]for row in case['raw_L_homogeneous_kernel']]).T;oldpoint=s.Matrix([s.Rational(v)for v in case['raw_L_particular']]);assert rank(full)==rank(oldK)==30 and rank(full.row_join(oldK))==30 and rank(full.row_join(oldpoint-point))==30
        piv=list(DomainMatrix.from_Matrix(full.T).convert_to(s.QQ).rref()[1]);recovery=s.zeros(30,50);inverse=full[piv,:].inv()
        for k,j in enumerate(piv):recovery[:,j]=inverse[:,k]
        assert recovery*full==s.eye(30);assert full*recovery*oldK==oldK and full*recovery*(oldpoint-point)==oldpoint-point
        P=(w*L+q*A).applyfunc(s.expand);Ph=s.Matrix(4,4,lambda i,j:s.expand((P[i,j]-x[i]*P[4,j]+x[j]*P[4,i]).subs(w,1)))
        B=(q*C+wedge(u,X)+D-wedge(u,aa)).applyfunc(s.expand);T=(aa-dd-cc*u-X).applyfunc(s.expand);Y=(-u*q-u*bb.dot(u)-MM*u-ee).applyfunc(s.expand)
        expected=s.zeros(4);expected[:3,:3]=B;expected[:3,3]=Y+z*T;expected[3,:3]=-(Y+z*T).T;assert(Ph-expected).applyfunc(s.expand)==s.zeros(4)
        Jaff=jac(Ph,(*u,z));Jraw=jac(P,x);PJ={}
        for quad in it.combinations(range(5),4):PJ[quad]=s.expand(sum((-1)**k*x[i]*Jraw[quad[:k]+quad[k+1:]]for k,i in enumerate(quad)))
        for tri,value in Jaff.items():assert s.expand(value+PJ[(*tri,4)].subs(w,1))==0
        assert s.expand(sum((-1)**i*x[i]*PJ[tuple(j for j in range(5)if j!=i)]for i in range(5)))==0
        BJac=jac(B,u)[(0,1,2)];LT=lie(T,B,u);FY=lie(Y,B,u)-wedge(Y,T)
        assert s.expand(Jaff[(0,1,2)]-BJac)==0
        for i,j in PAIRS3:assert s.expand(Jaff[(i,j,3)]+FY[i,j]+z*LT[i,j])==0
        # Highest spatial degrees of BJac are linear in all six parameters.
        spatial=s.Poly(BJac,*u);linear_equations=[value for exponent,value in spatial.terms()if sum(exponent)>=1];M0,rhs0=s.linear_eq_to_matrix(linear_equations,(*DD,*aa));assert rhs0==s.zeros(M0.rows,1)and rank(M0)==3
        p0,p1,delta=s.symbols('p0 p1 delta');a=s.Matrix([p0,0,p1])if ij==(0,2)else s.Matrix([p0,p1,0]);Dspecial=delta*C
        B=(q*C+wedge(u,X)+Dspecial-wedge(u,a)).applyfunc(s.expand);assert jac(B,u)[(0,1,2)]==0
        substitute={**dict(zip(aa,a)),**dict(zip(DD,[Dspecial[i,j]for i,j in PAIRS3]))};assert M0*s.Matrix([*dict(zip(DD,[Dspecial[i,j]for i,j in PAIRS3])).values(),*a])==s.zeros(M0.rows,1)
        t=s.Matrix(s.symbols('t0:3'));T=t-c*u-X;LT=lie(T,B,u);linear_spatial=[value for f in [LT[i,j]for i,j in PAIRS3]for exponent,value in s.Poly(f,*u).terms()if sum(exponent)==1];Mt,rhst=s.linear_eq_to_matrix(linear_spatial,tuple(t));assert not any(v.free_symbols for v in Mt)
        tsolution,tfree,compat=eliminate_constant_linear(Mt,rhst,t);assert rank(Mt)==(3 if ij==(0,2)else 2)
        residual=[s.factor(v)for v in coeffs([LT[i,j].subs(tsolution,simultaneous=True)for i,j in PAIRS3],u)if v!=0];residual=list(dict.fromkeys(residual));compat=list(dict.fromkeys(map(s.factor,compat)))
        all_conditions=[*compat,*residual]
        delta_equation=next(v for v in all_conditions if s.Poly(v,delta).degree()==1 and s.cancel(s.Poly(v,delta).coeff_monomial(delta)/c).is_Rational and s.cancel(s.Poly(v,delta).coeff_monomial(delta)/c)!=0)
        zero_c_nonzero=None;nonzero_c_a_constraint=None
        if ij==(0,1):
            nonzero_c_a_constraint=next(v for v in compat if s.cancel(v/(c*p1)).is_Rational and v!=0)
            zero_c_nonzero=next(s.expand(v.subs(c,0))for v in all_conditions if s.expand(v.subs(c,0))!=0)
        # The Y-system has an affine fixed cubic term. Its homogeneous fiber
        # is upper-triangular by spatial degree with these diagonal blocks.
        B0=q*C+wedge(u,X);T1=-X-c*u;Yb=[-u*u[j]for j in range(3)];YM=[];Ye=[]
        for i,j in it.product(range(3),repeat=2):v=s.zeros(3,1);v[i]=-u[j];YM.append(v)
        for i in range(3):v=s.zeros(3,1);v[i]=-1;Ye.append(v)
        matrices=[];blocks=[];roots={s.Rational(0)}
        for label,ys,degree in [('radial_quadratic',Yb,3),('linear',YM,2),('constant',Ye,1)]:
            matrix,packet,rr=principal_block(ys,B0,T1,u,c,degree);packet['name']=label;matrices.append(matrix);blocks.append(packet);roots.update(rr)
        generic_fiber=sum(M.cols-p['generic_rank']for M,p in zip(matrices,blocks));exceptions=[]
        for root in sorted(roots):
            ranks=[];certs=[]
            for M in matrices:
                value=M.subs(c,root);DM=DomainMatrix.from_Matrix(value).convert_to(s.QQ);RR,pc=DM.rref();rs=list(DomainMatrix.from_Matrix(value.T).convert_to(s.QQ).rref()[1]);pc=list(pc);det=value[rs,pc].det();assert det!=0
                ranks.append(len(pc));certs.append(dict(rows=rs,columns=pc,minor=str(det)))
            fiber=sum(M.cols-r for M,r in zip(matrices,ranks));base_bound=3 if root==0 else 2;d=base_bound+fiber
            exceptions.append(dict(c=str(root),exact_principal_ranks=ranks,rank_minors=certs,Y_affine_fiber_dimension_upper=fiber,BT_base_dimension_upper=base_bound,normalized_source_dimension_upper=d,PGL_image_dimension_upper=offset+d))
        normalized=max([3+generic_fiber,*[v['normalized_source_dimension_upper']for v in exceptions]]);global_bound=offset+normalized
        records.append(dict(name=case['name'],complete_affine25_bridge=dict(parameter_order=list(map(str,parameters)),rawL_particular=list(map(str,point)),rawL_parameter_matrix=bridge.tolist(),Euler5_matrix=E.tolist(),all30_recovery_matrix=recovery.tolist(),all_old30_raw_kernel_and_affine_point_reconstructed=True),
            all_projective_Jacobi_to_affine_and_Ore_identities_verified=True,Ore_base_B=encode(B),Ore_T=list(map(str,T)),Ore_Y=list(map(str,Y)),
            whole_B_Poisson_linear_rank=3,whole_B_parameters=['p0','p1','delta'],full_T_linear_matrix=Mt.tolist(),T_complete_pivot_solution={str(k):str(v)for k,v in tsolution.items()},T_free_translation_indices=tfree,all_BT_compatibility_conditions=list(map(str,all_conditions)),nonzero_c_delta_linear_equation=str(delta_equation),nonzero_c_a1_equation=str(nonzero_c_a_constraint),zero_c_nonzero_hypersurface=str(zero_c_nonzero),whole_BT_base_dimension_upper=3,fixed_nonzero_c_BT_dimension_upper=2,
            all_principal_blocks=blocks,generic_Y_affine_fiber_dimension_upper=generic_fiber,all_exceptional_c_including_zero=exceptions,normalized_full_Poisson_source_dimension_upper=normalized,whole_PGL_image_dimension_upper=global_bound,closes_this_entire_branch_G=global_bound<=24,seconds=round(time.monotonic()-started,6)))
        Path('/tmp/quadric_cone_zero_ore_checkpoint.json').write_text(json.dumps(dict(cases=records),indent=2,default=str)+'\n')
    paths=[Path(__file__),sourcepath,ROOT/'research/classification/QUADRIC_CONE_ZERO_SOURCE.md']
    result=dict(task_id='G',all_passed=True,cases=records,all_quadric_cone_rankzero_G_images_at_most24=all(v['closes_this_entire_branch_G']for v in records),
        triangular_rank_logic='For fixed B,T, solve the15unknown Y coefficients from spatial degrees3,2,1 in that order. Each diagonal matrix depends only on c; all lower-parameter couplings change only affine right sides, so the sum of the three nullities is an upper bound, including inconsistent fibers. All roots of each chosen generic minor are included exactly. At fixed nonzero c the BT base is<=2, at c0 <=3.',
        no_GB_or_arbitrary_point_probes=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as f:json.dump(result,f,indent=2,default=str);f.write('\n')
    print(json.dumps(dict(all_passed=True,bounds=[v['whole_PGL_image_dimension_upper']for v in records],seconds=result['seconds'])))


if __name__=='__main__':main()
