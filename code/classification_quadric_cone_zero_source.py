#!/usr/bin/env python3
"""Complete quadric-cone source and two normalized mod-H Jacobi kernels.

Fraction arithmetic; no GB, parameter samples or Poisson-support claims.
All affine kernels impose necessary conormal and mod-H Jacobi conditions.
"""
import argparse,hashlib,itertools as it,json,resource,time
from pathlib import Path
from classification_cubic_scroll_abelian_review import Q,PAIRS,add,scale,mul,diff,parse,rref,monomials,jacobi_and_divisor,determinant

ROOT=Path(__file__).resolve().parents[1]


def kernel(rows,n):
    R=rref(rows);free=[j for j in range(n)if j not in R];K=[]
    for j in free:
        v=[Q(i==j)for i in range(n)]
        for p in R:v[p]=-R[p].get(j,Q(0))
        K.append(v)
    return R,free,K


def affine(rows,rhs,n):
    R=rref([list(row)+[b]for row,b in zip(rows,rhs)]);assert n not in R
    free=[j for j in range(n)if j not in R];point=[R[j].get(n,Q(0))if j in R else Q(0)for j in range(n)];K=[]
    for j in free:
        v=[Q(i==j)for i in range(n)]
        for p in R:v[p]=-R[p].get(j,Q(0))
        K.append(v)
    assert all(sum(c*v[j]for j,c in enumerate(row))==0 for row in rows for v in K)
    assert all(sum(c*point[j]for j,c in enumerate(row))==b for row,b in zip(rows,rhs))
    return R,free,point,K


def empty():return [[{}for _ in range(5)]for _ in range(5)]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(115,120));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    q=parse('x0*x2-x1**2',[f'x{i}'for i in range(5)]);w=parse('x4',[f'x{i}'for i in range(5)]);xs=[parse(f'x{i}',[f'x{j}'for j in range(5)])for i in range(5)]
    mons2=monomials(5,2);mons3=monomials(5,3);mons4=monomials(5,4);QUADS=list(it.combinations(range(5),4));raw=[]
    for i,j in PAIRS:
        for k in range(5):
            P=empty();P[i][j]=mul(w,xs[k]);P[j][i]=scale(P[i][j],-1);raw.append(P)
    for i,j in PAIRS:
        P=empty();P[i][j]=q;P[j][i]=scale(q,-1);raw.append(P)
    def vec(P):return [P[i][j].get(e,Q(0))for i,j in PAIRS for e in mons2]
    source=[vec(P)for P in raw];assert len(rref(source))==60
    vertical=[];eulerL=[]
    for a,b in it.product(range(5),repeat=2):
        P=empty()
        for i,j in PAIRS:P[i][j]=add(mul(xs[i],xs[b])if j==a else{},scale(mul(xs[j],xs[b]),-1)if i==a else{})
        vertical.append(vec(P))
    for k in range(5):eulerL.append([Q((j==k and l==i)-(i==k and l==j))for i,j in PAIRS for l in range(5)])
    assert len(rref(vertical))==24 and len(rref(vertical+source))==79
    assert len(rref(eulerL))==5
    H=[add(*(mul(P[4][j],diff(q,j))for j in range(5)))for P in raw]
    square_generators=[mul(mul(w,w),xs[j])for j in range(5)]+[mul(w,q)]
    full=[[p.get(e,Q(0))for p in H]+[-p.get(e,Q(0))for p in square_generators]for e in mons3]
    R,free,K=kernel(full,66);assert len(R)==17 and len(K)==49
    projected=[[v[i]for i in range(60)]for v in K];assert len(rref(projected))==49
    # Eliminate only membership coefficients; preserve the exact full raw-source condition.
    swapped=[row[60:]+row[:60]for row in full];RR=rref(swapped);constraints=[]
    for p,row in RR.items():
        if p>=6:constraints.append([row.get(6+j,Q(0))for j in range(60)])
    assert len(constraints)==11
    assert len(rref([[v[50+i]for i in range(10)]for v in projected]))==7
    assert all(all(sum(c*v[j]for j,c in enumerate(row))==0 for row in constraints)for v in projected)
    euler_pivots=sorted(rref(eulerL));minor=[[v[j]for j in euler_pivots]for v in eulerL];assert determinant(minor)
    cases=[]
    for name,base_pair in [('Delta_nonzero',(0,2)),('Delta_zero_nonzeroC',(0,1))]:
        Avec=[Q(0)]*10;Avec[PAIRS.index(base_pair)]=1;Avec[PAIRS.index((3,4))]=-1
        A=empty()
        for c,(i,j)in zip(Avec,PAIRS):A[i][j]=scale(q,c);A[j][i]=scale(q,-c)
        constant_matrix=[[Q(0)for _ in range(5)]for _ in range(5)]
        for c,(i,j)in zip(Avec,PAIRS):constant_matrix[i][j]=c;constant_matrix[j][i]=-c
        assert len(rref(constant_matrix))==4
        # Treat fixed q*A as formal parameter0, then set that parameter to1.
        J,_=jacobi_and_divisor([A,*raw[:50]]);cols=list(it.combinations_with_replacement(range(51),2));ci={pair:i for i,pair in enumerate(cols)};linear=[];rhs=[];all350=[]
        for index,row in enumerate(J):
            e=mons4[index%70]
            if e[4]:all350.append(dict(coordinates_have_w=True));continue
            assert all(a==0 for col,c in row.items()for a,b in [cols[col]])
            coefficients=[row.get(ci[(0,j+1)],Q(0))for j in range(50)];constant=row.get(ci[(0,0)],Q(0))
            linear.append(coefficients);rhs.append(-constant);all350.append(dict(coordinates_have_w=False,linear_coefficients=list(map(str,coefficients)),constant=str(constant)))
        conormal_rows=[row[:50]for row in constraints];conormal_rhs=[-sum(row[50+j]*Avec[j]for j in range(10))for row in constraints]
        combined=linear+conormal_rows;combined_rhs=rhs+conormal_rhs
        original_R,original_free,original_point,original_K=affine(combined,combined_rhs,50)
        assert all(sum(c*v[j]for j,c in enumerate(row))==0 for row in combined for v in eulerL)
        gauge=[[Q(j==i)for j in range(50)]for i in euler_pivots]
        quotient_R,quotient_free,point,basis=affine(combined+gauge,combined_rhs+[0]*5,50)
        assert len(basis)+5==len(original_K)
        cases.append(dict(name=name,constant_A_upper=list(map(str,Avec)),constant_A_rank=4,all350_modH_Jacobi_rows=all350,
            complete_modH_linear_shape=[len(linear),50],modH_Jacobi_linear_rank=len(rref(linear)),combined_conormal_Jacobi_affine_rank=len(original_R),complete_raw_L_affine_dimension=len(original_K),raw_L_particular=list(map(str,original_point)),raw_L_homogeneous_kernel=[list(map(str,v))for v in original_K],
            Euler5_kernel_verified=True,Euler_gauge_coordinate_indices=euler_pivots,Euler_gauge_minor=str(determinant(minor)),complete_projective_affine_dimension=len(basis),quotient_particular=list(map(str,point)),complete_quotient_kernel=[list(map(str,v))for v in basis],quotient_free_rawL_indices=quotient_free,
            unnormalized_A_stratum_affine_dimension=7 if name=='Delta_nonzero'else 6,
            subsequent_PGL_image_bound_formula='18+d'if name=='Delta_nonzero'else'17+d',future_sufficient_full_Poisson_fiber_dimension=6 if name=='Delta_nonzero'else 7,
            scope='Complete affine supersource imposing necessary conormal and every mod-w projective Jacobi coefficient. Full Poisson support and a G bound are not yet established.'))
    # Universal inverse shear: lambdaG=lam*I-b*e_w; no exceptional b/C pivot.
    names=['lam','b0','b1','b2','c01','c02','c12'];one=parse('1',names);lam=parse('lam',names);bs=[parse(f'b{i}',names)for i in range(3)];Cs=[parse(v,names)for v in ['c01','c02','c12']]
    AA=[[{}for _ in range(5)]for _ in range(5)]
    for c,(i,j)in zip(Cs,[(0,1),(0,2),(1,2)]):AA[i][j]=c;AA[j][i]=scale(c,-1)
    for i in range(3):AA[i][3]=bs[i];AA[3][i]=scale(bs[i],-1)
    AA[4][3]=lam;AA[3][4]=scale(lam,-1)
    GG=[[lam if i==j else{}for j in range(5)]for i in range(5)];GI=[[lam if i==j else{}for j in range(5)]for i in range(5)]
    for i in range(3):GG[i][4]=scale(bs[i],-1);GI[i][4]=bs[i]
    def pm(A,B):return [[add(*(mul(a,b)for a,b in zip(row,col)))for col in zip(*B)]for row in A]
    product=pm(GG,GI);assert product==[[mul(lam,lam)if i==j else{}for j in range(5)]for i in range(5)]
    transformed=pm(pm(GG,AA),list(map(list,zip(*GG))));normal=[row[:]for row in AA]
    for i in range(3):normal[i][3]={};normal[3][i]={}
    assert transformed==[[mul(mul(lam,lam),v)for v in row]for row in normal]
    paths=[Path(__file__),ROOT/'code/classification_cubic_scroll_abelian_review.py',ROOT/'research/classification/RANK_ZERO_SMOOTH_NONPLANE_SEPARATE_REVIEW.md',ROOT/'research/classification/RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md']
    result=dict(task_id='G',all_passed=True,quadric='x0*x2-x1^2',vertex_coordinate='x3',hyperplane_coordinate='x4',raw_source_parameter_order=dict(L='pair(i,j) lex, then x_k coefficient;50 variables',A='pair(i,j)lex;10 variables'),raw_source_dimension=60,projective_source_dimension=55,Euler_redundancy_dimension=5,
        full_conormal_matrix=[[str(c)for c in row]for row in full],full_conormal_shape=[35,66],full_conormal_rank=17,complete_conormal_raw_source_dimension=49,complete_conormal_kernel_including_membership=[list(map(str,v))for v in K],all11_eliminated_raw_source_constraints=[list(map(str,v))for v in constraints],allowed_constant_A_dimension=7,
        inverse_normal_shear_and_constant_A_identity_verified=True,shear_formula='x_i_new=x_i-b_i*x4/lam; inverse x_i=x_i_new+b_i*x4/lam; q changes only by x4*linear; x3_new=x3/lam then makes A^{w,z}=1.',
        cases=cases,full_A_cover='v=0 or C=0 gives H|D; otherwise two full SO(q) plus tensor-scalar orbits Delta=4C01*C12-C02^2 nonzero orzero with Cnonzero. No C point outside those strata is discarded.',normalization_dimension_account='Aut(S) dimension12, scalar1. Allowed A strata dimensions7/6, normalized stabilizers in Aut times scalar6/7. If full normalized Poisson fiber has dimension d, PGL image <=18+d /17+d.',
        no_G_bound_claimed=True,no_GB_run=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,raw_affine=[r['complete_raw_L_affine_dimension']for r in cases],projective_affine=[r['complete_projective_affine_dimension']for r in cases],seconds=result['seconds'])))


if __name__=='__main__':main()
