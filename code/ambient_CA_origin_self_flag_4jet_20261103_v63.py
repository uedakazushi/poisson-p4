#!/usr/bin/env python3
"""All150 fourth-order CA ordinary paired-flag obstructions, with all earlier freedoms."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_V1_independent_review_20261103_v63 import tensor_factor_operator,act,rational

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63'
ZERO=(0,0,0,0)


class LP(dict):
    """Sparse QQ[t,t^-1,y0,y1,y2]; no rational-function simplifier is needed."""
    def __init__(self,value=0):
        if isinstance(value,dict):super().__init__((m,Q(c))for m,c in value.items()if c)
        else:super().__init__({ZERO:Q(value)}if value else{})
    def __add__(self,other):
        other=other if isinstance(other,LP)else LP(other);out=LP(self)
        for m,c in other.items():
            out[m]=out.get(m,Q(0))+c
            if not out[m]:del out[m]
        return out
    __radd__=__add__
    def __neg__(self):return LP({m:-c for m,c in self.items()})
    def __sub__(self,other):return self+-LP(other)
    def __rsub__(self,other):return LP(other)+-self
    def __mul__(self,other):
        other=other if isinstance(other,LP)else LP(other);out={}
        for a,c in self.items():
            for b,d in other.items():
                m=tuple(x+y for x,y in zip(a,b));out[m]=out.get(m,Q(0))+c*d
        return LP(out)
    __rmul__=__mul__
    def __truediv__(self,c):return self*Q(1,c)
    def __pow__(self,n):
        assert n>=0;out=LP(1)
        for _ in range(n):out*=self
        return out
    def expr(self,variables):return sum(s.Rational(c.numerator,c.denominator)*s.prod(v**e for v,e in zip(variables,m))for m,c in self.items())
    def save(self):return[[list(m),str(c)]for m,c in sorted(self.items())]
    def evaluate(self,values):return sum((c*Q(values[0])**m[0]*Q(values[1])**m[1]*Q(values[2])**m[2]*Q(values[3])**m[3]for m,c in self.items()),Q(0))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def QQM(rows):return s.Matrix([[s.Rational(c)for c in row]for row in rows])
def vadd(*vectors):return[sum((v[i]for v in vectors),LP())for i in range(len(vectors[0]))]
def vscale(v,c):return[c*x for x in v]
def cmul(M,V):
    return[sum((rational(M[i,j])*V[j]for j in range(M.cols)if M[i,j]),LP())for i in range(M.rows)]
def pmat(M):return[[LP(rational(M[i,j]))for j in range(M.cols)]for i in range(M.rows)]
def madd(*matrices):return[[sum((M[i][j]for M in matrices),LP())for j in range(5)]for i in range(5)]
def mscale(M,c):return[[c*v for v in row]for row in M]
def mm(A,B):return[[sum((A[i][k]*B[k][j]for k in range(5)),LP())for j in range(5)]for i in range(5)]


def physical_pushforward4(gcoeff,source):
    """Independent quadratic-coordinate substitution: exterior g and inverse coordinates."""
    N=5;zero=[Q(0)]*N
    def add(*vectors):return[sum((v[i]for v in vectors),Q(0))for i in range(N)]
    def mul(a,b):return[sum((a[j]*b[i-j]for j in range(i+1)),Q(0))for i in range(N)]
    def neg(a):return[-c for c in a]
    def matmul(A,B):return[[add(*[mul(A[i][k],B[k][j])for k in range(5)])for j in range(5)]for i in range(5)]
    g=[[[gcoeff[d][i][j]for d in range(N)]for j in range(5)]for i in range(5)]
    ident=[[[Q(i==j),*([Q(0)]*4)]for j in range(5)]for i in range(5)]
    nil=[[add(g[i][j],neg(ident[i][j]))for j in range(5)]for i in range(5)];inverse=ident;power=ident
    for k in range(1,5):
        power=matmul(power,nil);inverse=[[add(inverse[i][j],power[i][j]if k%2==0 else neg(power[i][j]))for j in range(5)]for i in range(5)]
    assert matmul(g,inverse)==ident and matmul(inverse,g)==ident
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));lookup={m:i for i,m in enumerate(mons)}
    wedges=[]
    for i,j in pairs:wedges.append([add(mul(g[i][k],g[j][l]),neg(mul(g[i][l],g[j][k])))for k,l in pairs])
    substituted=[]
    for m,n in mons:
        result=[zero[:]for _ in mons]
        for a,b in itertools.product(range(5),repeat=2):
            j=lookup[tuple(sorted((a,b)))];result[j]=add(result[j],mul(inverse[m][a],inverse[n][b]))
        substituted.append(result)
    output=[zero[:]for _ in range(150)]
    for k in range(10):
        for mon in range(15):
            p=source[15*k+mon]
            if not any(p):continue
            for out in range(10):
                if not any(wedges[out][k]):continue
                w=mul(wedges[out][k],p)
                for new in range(15):
                    if any(substituted[mon][new]):output[15*out+new]=add(output[15*out+new],mul(w,substituted[mon][new]))
    return output


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    base=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63'
    oldreport=read(base/'report.json');assert oldreport['all_passed']
    accepted3=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63_REVIEW/report.json';assert read(accepted3)['all_passed']
    setup=read(base/'setup.json');old=read(base/'certificate.json.gz')
    source=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';data=read(source)
    C=QQM(data['fixed_fiber_linear_basis150_by20']);free=data['affine_slice_free_indices'];assert len(free)==10
    P0=cmul(C,[LP(Q(c))for c in data['center_fixed20']]);F=[cmul(C[:,j],[LP(1)])for j in free]
    selected=list(DomainMatrix.from_Matrix(C.T).convert_to(s.QQ).rref()[1]);left=C[selected,:].inv()
    def project(v):return vadd(v,vscale(cmul(C,cmul(left,[v[j]for j in selected])),-1))
    D=QQM(setup['all150_outside_fiber_order2_linear_matrix']);cols=setup['pivot_columns'];rows=setup['pivot_rows'];inv=QQM(setup['inverse_pivot'])
    def coker(v):return vadd(v,vscale(cmul(D[:,cols],cmul(inv,[v[j]for j in rows])),-1))
    d=[LP(Q(c))for c in old['order2_affine_origin']];T=QQM(old['order2_affine_embedding'])
    t=LP({(1,0,0,0):1});ti=LP({(-1,0,0,0):1});ys=[LP({tuple(int(j==i+1)for j in range(4)):1})for i in range(3)]
    low=vadd(d,cmul(T,[-7*t/4-1+3*ti/4,t-1,LP()]))
    products=[]
    for label in old['order3_low_monomials']:
        value=LP(1)
        for j in label:value*=low[j]
        products.append(value)
    intrinsic=cmul(QQM(old['order3_raw_quotient_coefficient_matrix']),products)
    high=cmul(T,ys);solution=cmul(inv,[intrinsic[j]for j in rows])
    for i,col in enumerate(cols):high[col]-=solution[i]
    assert not any(vadd(cmul(D,high),intrinsic))
    flag=list(map(tuple,setup['second_and_third_order_flag8']));assert len(flag)==8
    M=pmat(QQM(setup['tangent_M']));B=[[LP()for _ in range(5)]for _ in range(5)];C3=[[LP()for _ in range(5)]for _ in range(5)]
    for k,(i,j)in enumerate(flag):B[i][j]=low[10+k];C3[i][j]=high[10+k]
    M2=mm(M,M);M3=mm(M2,M);M4=mm(M3,M);assert not any(v for row in M4 for v in row)
    L2=madd(B,mscale(M2,Q(-1,2)))
    L3=madd(C3,mscale(madd(mm(M,B),mm(B,M)),Q(-1,2)),mscale(M3,Q(1,3)))
    L4=madd(mscale(madd(mm(M,C3),mm(C3,M),mm(B,B)),Q(-1,2)),mscale(madd(mm(M2,B),mm(mm(M,B),M),mm(B,M2)),Q(1,3)))
    logs=[None,M,L2,L3,L4]
    # Verify the noncommutative5x5 logarithm/exponential identity itself, through degree4.
    zmat=[[LP()for _ in range(5)]for _ in range(5)];ident=pmat(s.eye(5));term=[ident,*[zmat for _ in range(4)]];exp=[ident,*[zmat for _ in range(4)]]
    for power in range(1,5):
        nxt=[zmat]
        for degree in range(1,5):nxt.append(mscale(madd(*[mm(logs[k],term[degree-k])for k in range(1,degree+1)]),Q(1,power)))
        term=nxt;exp=[madd(a,b)for a,b in zip(exp,term)]
    assert exp==[ident,M,B,C3,zmat]
    ops=[tensor_factor_operator(*divmod(i,5))for i in range(25)]
    def action(matrix,vector):
        ans=[LP()for _ in range(150)]
        for i,j in itertools.product(range(5),repeat=2):
            if matrix[i][j]:ans=vadd(ans,vscale(act(ops[5*i+j],vector),matrix[i][j]))
        return ans
    P1=[sum((low[j]*F[j][i]for j in range(10)),LP())for i in range(150)]
    P2=[sum((high[j]*F[j][i]for j in range(10)),LP())for i in range(150)]
    zero=[LP()for _ in range(150)];term=[P0,P1,P2,zero,zero];exp=[v[:]for v in term]
    for power in range(1,5):
        nxt=[zero]
        for degree in range(1,5):nxt.append(vscale(vadd(*[action(logs[k],term[degree-k])for k in range(1,degree+1)]),Q(1,power)))
        term=nxt;exp=[vadd(a,b)for a,b in zip(exp,term)]
    assert not any(project(exp[1]))and not any(project(exp[2]))and not any(project(exp[3]))
    fourth=project(exp[4]);obstruction=coker(fourth)
    assert all(sum(m[1:])<=1 for p in obstruction for m in p)
    vars=s.symbols('t y0 y1 y2');words=sorted(set(m for p in obstruction for m in p));coeff=s.Matrix([[s.Rational(p.get(m,0))for m in words]for p in obstruction])
    eqrows=list(DomainMatrix.from_Matrix(coeff.T).convert_to(s.QQ).rref()[1]);equations=[s.factor(obstruction[i].expr(vars))for i in eqrows]
    # Save reusable exact input BEFORE the small residual solver.
    payload=dict(parent_3jet_report_sha256=sha(base/'report.json'),variables=list(map(str,vars)),laurent_polynomial_encoding='[exponents(t,y0,y1,y2), rational coefficient]',
        full_order2_solution=[p.save()for p in low],full_order3_solution_with_all3_higher_freedoms=[p.save()for p in high],
        group_log1_to4=[[[p.save()for p in row]for row in A]for A in logs[1:]],
        all150_order4_before_fresh18_elimination=[p.save()for p in fourth],all150_order4_coker=[p.save()for p in obstruction],
        independent_coker_rows=eqrows,independent_coker_equations=list(map(str,equations)),fresh_order4_source10_plus_flag8_image_rank=len(cols),
        all_earlier_order3_free3_retained=True,all150_order1_to3_membership_verified=True,group_log_exp_order4_checked=True)
    with(out/'certificate.json.gz').open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as z:z.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    print(json.dumps(dict(stage='all150_order4_exact_input_saved',seconds=time.monotonic()-start,independent_coker_equations=list(map(str,equations)))),flush=True)
    # Laurent unit t is retained exactly as q*t-1; only five small variables.
    q=s.symbols('q');cleared=[s.Poly(s.expand(e*vars[0]**max(0,-min(m[0]for m in obstruction[i]))),*vars).as_expr()for i,e in zip(eqrows,equations)]
    gb=s.groebner([*cleared,q*vars[0]-1],q,*vars,domain=s.QQ,order='lex')
    assert len(equations)==4
    tv,y0,y1,y2=vars
    target=[y2,(tv*tv-1)**2,4*tv*tv*y0+(7*tv*tv+3)*y1]
    cofactors=[[3*tv/2,0,0,-2*tv/3],[0,2*tv/3,0,0],[0,0,-2*tv,0]]
    assert all(s.expand(sum(c*e for c,e in zip(cofs,equations))-wanted)==0 for cofs,wanted in zip(cofactors,target))
    assert all(s.rem(s.together(e).as_numer_denom()[0],s.Poly(y2,y2)).as_expr()==0 for e in[equations[0],equations[3]])
    witnesses=[]
    for sign in[1,-1]:
        values=[Q(sign),Q(0),Q(0),Q(0)];fresh=[LP()for _ in range(18)];right=cmul(inv,[fourth[j]for j in rows])
        for i,col in enumerate(cols):fresh[col]=LP(-right[i].evaluate(values))
        residual=vadd(fourth,cmul(D,fresh));assert not any(p.evaluate(values)for p in residual)
        P3=[sum((fresh[j]*F[j][i]for j in range(10)),LP())for i in range(150)]
        flag4=[[LP()for _ in range(5)]for _ in range(5)]
        for k,(i,j)in enumerate(flag):flag4[i][j]=fresh[10+k]
        expected=[v[:]for v in exp];expected[3]=vadd(expected[3],P3);expected[4]=vadd(expected[4],action(M,P3),action(flag4,P0))
        gcoeff=[[[p.evaluate(values)for p in row]for row in mat]for mat in[ident,M,B,C3,flag4]]
        sourcecoeff=[[p.evaluate(values)for p in v]for v in[P0,P1,P2,P3,zero]]
        direct=physical_pushforward4(gcoeff,list(map(list,zip(*sourcecoeff))))
        assert direct==[[expected[d][i].evaluate(values)for d in range(5)]for i in range(150)]
        assert all(not any(project([LP(direct[i][degree])for i in range(150)]))for degree in range(5))
        witnesses.append(dict(t=str(sign),earlier_higher_free3=['0','0','0'],fresh_order4_source10_and_flag8=list(map(lambda p:str(p.evaluate(values)),fresh)),
            source_physical150_by5=list(map(lambda row:list(map(str,row)),zip(*sourcecoeff))),
            group5_by5_by5=[[[str(c)for c in row]for row in mat]for mat in gcoeff],
            second_source_physical150_by5=[list(map(str,row))for row in direct],all150_direct_physical_pushforward_and_membership_verified=True))
    payload.update(exact_Laurent_ideal_generators=list(map(str,target)),exact_generator_cofactor_matrix=list(map(lambda row:list(map(str,row)),cofactors)),
        quotient_isomorphism='QQ[t,t^-1,y1]/((t^2-1)^2); y2=0, y0=-(7*t^2+3)*y1/(4*t^2)',finite4jet_witnesses=witnesses)
    with(out/'certificate.json.gz').open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as z:z.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,all150_physical_membership_rows_retained=True,CA20_free10_and_all9_flags_retained=True,
        order4_fresh18_rank=len(cols),order3_curve_and_all3_earlier_higher_freedoms_retained=True,order4_coker_constant_coefficient_rank=len(eqrows),
        order4_coker_equations=list(map(str,equations)),small_Laurent_QQ_ideal_GB=list(map(str,gb)),
        exact_Laurent_ideal_generators=list(map(str,target)),Gm_base_cut_to_two_doubled_points_t_plus_minus1=True,
        full_necessary_4jet_solution_is_two_doubled_lines_with_fresh_order4_free3=True,
        both_finite4jet_witnesses_all150_direct_physical_checked=True,
        necessary_ordinary4jet_empty=any(p==1 for p in gb),actual_formal_or_algebraic_pair_branch_proved=False,
        all_order_nonordinary_exclusion_or_individual_conductor_claimed=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),source,base/'report.json',base/'setup.json',base/'certificate.json.gz',accepted3,ROOT/'code/ambient_CA_V1_independent_review_20261103_v63.py']
    note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
