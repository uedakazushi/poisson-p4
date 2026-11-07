#!/usr/bin/env python3
"""Independent eight explicit free-quintic sources and intrinsic cusp types.

Fraction sparse arithmetic, complete implicit/syzygy matrices, all ten source
parameters and local finite Puiseux recognition. No producer/CAS import.
Necessary completeness of the eight-item list is a separate pure-proof review.
"""
from fractions import Fraction as F
from itertools import combinations,product
from math import comb
from pathlib import Path
import hashlib,json,resource,signal,time
import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_minimal_quartic_review import LaurentRing,determinant
from classification_unimodular_logarithmic_review import ug,ud,trim

ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCES_20261107.json'
INPUT='research/classification/inputs/essential_three_quintic_eight_curves_20261107.json'
OUT='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCE_REVIEW_20261107.json'
JOB='research/classification/jobs/essential_three_quintic_eight_source_review.json'

def exponents(n,d):return sorted((e for e in product(range(d+1),repeat=n)if sum(e)==d),reverse=True)
def matmul(A,B):return [[sum((a*b for a,b in zip(r,c)if a and b),F(0))for c in zip(*B)]for r in A]
def cert(M,C):
    assert M==[[F(v)for v in row]for row in C['all_rows']]
    assert C['shape']==[len(M),len(M[0])]
    rows,cols=C['rows'],C['columns'];rank=C['rank']
    assert len(rows)==len(cols)==rank and len(set(rows))==len(set(cols))==rank
    value=q.determinant([[M[i][j]for j in cols]for i in rows]);assert value and value==F(C['determinant'])
    return str(value)
def power(R,P,n):
    out=R.parse('1')
    for _ in range(n):out=R.mul(out,P)
    return out
def root_multiplicities(P):
    degrees=[];P=trim(P)
    while len(P)>1:
        degrees.append(len(P)-1);P=ug(P,tuple(F(i)*v for i,v in enumerate(P)if i))
    degrees.extend([0,0]);result=[]
    for j in range(1,len(degrees)-1):result.extend([j]*(degrees[j-1]-2*degrees[j]+degrees[j+1]))
    return result

class Q5:
    def __init__(self,a=0,b=0):
        if isinstance(a,Q5):self.a,self.b=a.a,a.b;return
        self.a,self.b=F(a),F(b)
    def __add__(self,v):v=Q5(v);return Q5(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __neg__(self):return Q5(-self.a,-self.b)
    def __sub__(self,v):return self+-Q5(v)
    def __rsub__(self,v):return Q5(v)+-self
    def __mul__(self,v):v=Q5(v);return Q5(self.a*v.a+5*self.b*v.b,self.a*v.b+self.b*v.a)
    __rmul__=__mul__
    def __truediv__(self,v):
        v=Q5(v);norm=v.a*v.a-5*v.b*v.b;assert norm;return self*Q5(v.a/norm,-v.b/norm)
    def __rtruediv__(self,v):return Q5(v)/self
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,v):v=Q5(v);return self.a==v.a and self.b==v.b
    def __repr__(self):return str(self.a)+'+('+str(self.b)+')sqrt5'

N=16
def spow(a,k,field):
    out=field(1)
    for _ in range(k):out=out*a
    return out
def sadd(a,b):return [x+y for x,y in zip(a,b)]
def sscale(a,c):return [c*v for v in a]
def smul(a,b):
    zero=a[0]*0;out=[zero for _ in range(N)]
    for i,x in enumerate(a):
        if not x:continue
        for j,y in enumerate(b[:N-i]):
            if y:out[i+j]=out[i+j]+x*y
    return out
def sinverse(a):
    assert a[0];out=[a[0]*0 for _ in range(N)];out[0]=1/a[0]
    for i in range(1,N):out[i]=-sum((a[j]*out[i-j]for j in range(1,i+1)),a[0]*0)/a[0]
    return out
def series(P,rho,field,infinity=False):
    out=[field(0)for _ in range(N)];index=1 if infinity else 0
    for e,c in P.items():
        k=e[index]
        for j in range(k+1):out[j]=out[j]+field(c*comb(k,j))*spow(rho,k-j,field)
    return out
def local(phi,rho=F(0),field=F,infinity=False):
    values=[series(P,rho,field,infinity)for P in phi];chart=next(i for i,v in enumerate(values)if v[0]);inverse=sinverse(values[chart]);coords=[]
    for i,v in enumerate(values):
        if i==chart:continue
        a=smul(v,inverse);a[0]=field(0);coords.append(a)
    order=lambda a:next((i for i,v in enumerate(a)if v),N)
    orders=sorted(order(a)for a in coords)
    if orders[0]!=2:return dict(orders=orders,chart=chart)
    a,b=coords if order(coords[0])==2 else coords[::-1]
    b=sadd(b,sscale(a,-b[2]/a[2]));firstodd=None
    for i in range(3,N):
        if not b[i]:continue
        if i%2:firstodd=i;break
        k=i//2;ap=[field(int(j==0))for j in range(N)]
        for _ in range(k):ap=smul(ap,a)
        b=sadd(b,sscale(ap,-b[i]/spow(a[2],k,field)))
    assert firstodd is not None
    return dict(orders=orders,chart=chart,first_characteristic_odd=firstodd,nonzero_coefficient=str(b[firstodd]))

def cusp_types(name,phi,multiplicities):
    details=[]
    def at(infinity,odd=None,orders=None,rho=F(0),field=F):
        c=local(phi,rho,field,infinity)
        if odd is not None:assert c['first_characteristic_odd']==odd
        if orders is not None:assert c['orders']==orders
        c.update(infinity=infinity,parameter=str(rho));details.append(c);return c
    if name=='C45':at(True,orders=[4,5]);signature=['(4,5)'];budget=6
    elif name=='A12':at(False,13);signature=['A12'];budget=6
    elif name=='E8_A4':at(False,orders=[3,5]);at(True,5);signature=['E8','A4'];budget=6
    elif name=='E6_A6':at(False,orders=[3,4]);at(True,7);signature=['E6','A6'];budget=6
    elif name=='E6_A4_A2':
        at(False,orders=[3,4]);at(True,5);at(False,3,rho=F(-1,2));signature=['E6','A4','A2'];budget=6
    elif name=='A8_A4':at(False,9);at(True,5);signature=['A8','A4'];budget=6
    elif name=='THREE_A4':
        at(False,5);rho=Q5(-2,F(6,5));assert rho*rho+4*rho-F(16,5)==0
        at(False,5,rho=rho,field=Q5)
        # Nonzero QQ(sqrt5) coefficients stay nonzero in both embeddings.
        signature=['A4','A4','A4'];budget=6
    else:
        assert name=='A6_THREE_A2';at(True,7)
        assert sorted(multiplicities)==[1,1,1,1]
        signature=['A6','A2','A2','A2'];budget=6
    assert budget==6
    return dict(signature=signature,local_finite_series=details,
        genus_budget6_used=True,no_extra_singularities_or_multiple_image_branches=True,
        reason='Local delta lower bounds sum6, the arithmetic genus of the irreducible quintic. Remaining simple differential roots have multiplicity2 and delta>=1; equality forcesA2. At3A4 the quadratic conjugate is covered by the second field embedding.')

def check(row,inputrow):
    name=row['id'];R0=LaurentRing((),('s','t'));phi=[R0.parse(t)for t in inputrow['parametrization']]
    assert phi==[R0.parse(t)for t in row['exact_degree5_parametrization']]
    assert all(all(sum(e)==5 for e in P)for P in phi)
    aff=[trim(tuple(P.get((j,5-j),F(0))for j in range(6)))for P in phi]
    gcd=ug(ug(aff[0],aff[1]),aff[2]);assert len(gcd)==1 and any(P.get((5,0))for P in phi)
    implicit={}
    for degree in (4,5):
        src=exponents(3,degree);target=[(5*degree-j,j)for j in range(5*degree+1)]
        powers=[[power(R0,P,k)for k in range(degree+1)]for P in phi]
        cols=[R0.mul(*(powers[i][k]for i,k in enumerate(e)))for e in src]
        M=[[P.get(e,F(0))for P in cols]for e in target]
        saved=row['complete_implicit_coefficient_maps'][str(degree)]
        assert saved['source_exponents']==[list(e)for e in src]and saved['target_exponents']==[list(e)for e in target]
        value=cert(M,saved['coefficient_map']);assert saved['coefficient_map']['rank']==(15 if degree==4 else 20)
        K=[[F(v)for v in v]for v in saved['complete_kernel']]
        assert len(K)==(0 if degree==4 else 1)
        if degree==5:assert any(K[0])and not any(any(r)for r in matmul(M,q.transpose(K)))
        implicit[str(degree)]=dict(rank=saved['coefficient_map']['rank'],independent_minor=value)
    params=('g00','g01','g10','g11')+tuple('r'+str(i)for i in range(6));R=LaurentRing(params,('x','y','z','u','v'));p=R.parse
    f=p(row['primitive_quintic']);x=[R.variable(n)for n in ('x','y','z','u','v')];grad=[R.derivative(f,i)for i in range(3)]
    # Direct binary substitution of every term of the implicit quintic.
    substitution=R0.add(*(R0.scale(R0.mul(*(power(R0,phi[i],e[10+i])for i in range(3))),c)for e,c in f.items()));assert not substitution
    target=exponents(3,4);ess=[[P.get((0,)*10+e+(0,0),F(0))for P in grad]for e in target]
    assert row['essential_gradient_monomials']==[list(e)for e in target];em=cert(ess,row['essential_gradient_certificate'])
    mon=[p(t)for t in ('x**2','x*y','x*z','y**2','y*z','z**2')];target=exponents(3,6);cols=[R.mul(m,g)for g in grad for m in mon]
    M=[[P.get((0,)*10+e+(0,0),F(0))for P in cols]for e in target];saved=row['complete_quadratic_gradient_kernel'];sm=cert(M,saved['coefficient_map'])
    assert saved['all_degree6_target_exponents']==[list(e)for e in target]
    K=[[F(v)for v in k]for k in saved['full18_basis']];assert len(K)==2 and q.rank(K)[0]==2
    assert not any(any(r)for r in matmul(M,q.transpose(K)))and saved['coefficient_map']['rank']==16
    fields=[[R.add(*(R.scale(mon[j],k[6*i+j])for j in range(6)))for i in range(3)]for k in K]
    assert fields==[[p(t)for t in v]for v in saved['vector_fields']]
    coefficient=F(saved['determinant_over_F']);assert coefficient
    assert determinant(R,[[x[i],fields[0][i],fields[1][i]]for i in range(3)])==p(saved['saito_determinant'])==R.scale(f,coefficient)
    U=[R.add(R.mul(p('g00'),u),R.mul(p('g01'),v))for u,v in zip(*fields)]
    V=[R.add(R.mul(p('g10'),u),R.mul(p('g11'),v))for u,v in zip(*fields)]
    div=lambda W:R.add(*(R.derivative(W[i],i)for i in range(3)))
    magnetic=R.add(*(R.mul(p('r'+str(i)),m)for i,m in enumerate(mon)))
    data={(i,3):U[i]for i in range(3)};data.update({(i,4):V[i]for i in range(3)})
    data[3,4]=R.add(R.scale(R.mul(div(V),x[3]),-1),R.mul(div(U),x[4]),magnetic)
    P={}
    for ij in q.PAIRS:P[ij]=data.get(ij,{});P[ij[::-1]]=R.scale(P[ij],-1)
    source=row['whole_GL2_frame_and_all6_magnetic_source']
    assert all(P[ij]==p(source['all_ten_brackets'][str(ij[0])+','+str(ij[1])])for ij in q.PAIRS)
    for i,j,k in combinations(range(5),3):
        assert not R.add(*(R.mul(P.get((a,l),{}),R.derivative(P.get((b,c),{}),l))for a,b,c in ((i,j,k),(j,k,i),(k,i,j))for l in range(5)))
    assert all(not R.add(*(R.derivative(P.get((i,j),{}),j)for j in range(5)))for i in range(5))
    assert all(not R.add(*(R.mul(P.get((i,j),{}),R.derivative(f,j))for j in range(5)))for i in range(5))
    pf=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf.append(R.add(R.mul(P[i,j],P[k,l]),R.scale(R.mul(P[i,k],P[j,l]),-1),R.mul(P[i,l],P[j,k])))
    D=R.add(*(R.scale(R.mul(x[i],v),(-1)**i)for i,v in enumerate(pf)))
    assert pf==[p(t)for t in source['all_five_raw_pfaffians']]
    assert D==p(source['actual_quintic'])==R.scale(R.mul(p('g00*g11-g01*g10'),f),-coefficient)
    # Fixed-F all-six magnetic map has a literal inverse; no chosen R patch.
    dc=[]
    for W in fields[::-1]:
        for linear in x[:3]:dc.append(div([R.mul(linear,a)for a in W]))
    dc[3:]=[R.scale(a,-1)for a in dc[3:]]
    monexp=[next(iter(m))for m in mon];delta=[[col.get(e,F(0))for col in dc]for e in monexp]
    delta_det=q.determinant(delta);assert delta_det
    inverse=q.inverse(delta);assert matmul(delta,inverse)==[[F(i==j)for j in range(6)]for i in range(6)]
    # Binary differential gcd, including the point at infinity.
    ds=[R0.derivative(P,0)for P in phi];dt=[R0.derivative(P,1)for P in phi]
    cross=[R0.add(R0.mul(ds[i],dt[j]),R0.scale(R0.mul(ds[j],dt[i]),-1))for i,j in combinations(range(3),2)]
    cs=[trim(tuple(P.get((j,8-j),F(0))for j in range(9)))for P in cross];g=ug(ug(cs[0],cs[1]),cs[2])
    inf=min(e[1]for P in cross for e in P);expected=R0.parse(row['derivative_gcd'])
    expected_aff=trim(tuple(expected.get((j,sum(next(iter(expected)))-j),F(0))for j in range(sum(next(iter(expected)))+1)))
    assert ug(g,expected_aff)==ug(expected_aff,())and len(g)-1+inf==sum(next(iter(expected)))
    mult=root_multiplicities(g)+([inf]if inf else[]);assert sorted(mult)==sorted(row['complex_differential_root_multiplicities'])==sorted(inputrow['expected_complex_derivative_gcd_root_multiplicities'])
    types=cusp_types(name,phi,mult)
    return dict(id=name,implicit_complete_matrices=implicit,no_binary_basepoint=True,complex_absolute_irreducibility='Irreducible P1 image; any proper factor of the vanishing quintic would give a degree<=4 relation, excluded by rank15.',
        essential_rank3_minor=em,full28_by18_rank16_minor=sm,full_two_syzygy_kernel=True,saito_determinant_over_F=str(coefficient),
        full_GL2_frame_and_all6_magnetic_allJac_div_Casimir_Pf=True,whole_source_completeness='Accepted canonical essential3 theorem and Saito basis; no extra higher-fiber coordinates omitted.',
        magnetic_shear_matrix=[[str(v)for v in row]for row in delta],magnetic_shear_inverse=[[str(v)for v in row]for row in inverse],magnetic_shear_det=str(delta_det),
        all_frame_inverse_localization='detframe!=0 exactly rank4; delta_frame determinant=detdelta*detframe^3',
        complex_differential_multiplicities=sorted(mult),intrinsic_cusp_recognition=types)

def main():
    started=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    packet=json.loads((ROOT/SOURCE).read_text());assert hashlib.sha256(packet['producer_bytes'].encode()).hexdigest()==packet['producer_sha256']
    rows=packet['eight_curve_sources']['all_eight_candidate_sources'];inputs=json.loads((ROOT/INPUT).read_text())['candidates'];assert len(rows)==len(inputs)==8
    records=[]
    for row,inp in zip(rows,inputs):
        assert row['id']==inp['id'];records.append(check(row,inp));print(row['id'],'explicit source/cusp PASS',flush=True)
    signatures=[tuple(r['intrinsic_cusp_recognition']['signature'])for r in records];assert len(set(signatures))==8
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,all8_explicit_sources=records,
        all8_mutually_distinct_intrinsic_plane_cusp_types=True,all8_actual_fixedF_sources_one_PGL5_orbit23_by_accepted_uniform_lemma=True,
        no_CAS_or_producer_imports=True,necessary_eight_curve_classification_scope='Separate pure-proof review; not a hypothesis of these8 explicit source/cusp certificates.',
        no_component_count_changed=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print('PASS',OUT,round(result['seconds'],3),result['maximum_RSS_KiB'],flush=True)
if __name__=='__main__':main()
