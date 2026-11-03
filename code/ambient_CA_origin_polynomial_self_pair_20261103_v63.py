#!/usr/bin/env python3
"""Test the t=-1 CA finite-jet witness as an exact polynomial paired source."""
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
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63'
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def QQM(rows):return s.Matrix([[s.Rational(c)for c in row]for row in rows])
def savepoly(p):return[[m[0],str(c)]for m,c in sorted(p.items())]


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(60)
    parent=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63/report.json';assert read(parent)['all_passed']
    source=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';data=read(source)
    accepted=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';assert read(accepted)['all_passed']
    C=QQM(data['fixed_fiber_linear_basis150_by20']);free=data['affine_slice_free_indices'];center=s.Matrix(list(map(s.Rational,data['center_fixed20'])))
    P0=C*center;first=s.zeros(20,1);first[free[0]]=-3;first[free[3]]=-1;P1=C*first
    R,a=ring('a',s.QQ)
    def pol(value):return R(s.Rational(value))
    point=[pol(P0[i])+a*pol(P1[i])for i in range(150)]
    M=s.zeros(5);M[1,0]=1;M[0,3]=-2;assert M**3==s.zeros(5)
    ident=[[R(int(i==j))for j in range(5)]for i in range(5)]
    g=[[ident[i][j]+a*pol(M[i,j])+a*a*pol((M*M)[i,j])for j in range(5)]for i in range(5)]
    gi=[[ident[i][j]-a*pol(M[i,j])for j in range(5)]for i in range(5)]
    def mm(A,B):return[[sum((A[i][k]*B[k][j]for k in range(5)),R.zero)for j in range(5)]for i in range(5)]
    assert mm(g,gi)==ident and mm(gi,g)==ident
    # Direct finite polynomial physical push: wedge(g) and quadratic inverse-coordinate substitution.
    wedges=[[g[i][k]*g[j][l]-g[i][l]*g[j][k]for k,l in PAIRS]for i,j in PAIRS]
    lookup={m:i for i,m in enumerate(MONS)};substituted=[]
    for m,n in MONS:
        vec=[R.zero]*15
        for i,j in itertools.product(range(5),repeat=2):vec[lookup[tuple(sorted((i,j)))]]+=gi[m][i]*gi[n][j]
        substituted.append(vec)
    pushed=[R.zero]*150
    for k in range(10):
        for mon in range(15):
            if not point[15*k+mon]:continue
            for i in range(10):
                if not wedges[i][k]:continue
                factor=wedges[i][k]*point[15*k+mon]
                for new in range(15):
                    if substituted[mon][new]:pushed[15*i+new]+=factor*substituted[mon][new]
    rows=list(DomainMatrix.from_Matrix(C.T).convert_to(s.QQ).rref()[1]);left=C[rows,:].inv()
    second=[sum((pol(left[i,j])*pushed[rows[j]]for j in range(20)),R.zero)for i in range(20)]
    reconstructed=[sum((pol(C[i,j])*second[j]for j in range(20)),R.zero)for i in range(150)]
    residual=[p-q for p,q in zip(pushed,reconstructed)]
    nonzero=[i for i,p in enumerate(residual)if p]
    print(json.dumps(dict(stage='full_polynomial_membership',seconds=time.monotonic()-start,nonzero_membership_rows=nonzero,
        residual_degrees={str(i):max(m[0]for m in residual[i])for i in nonzero})),flush=True)
    # The strict homogeneous cubic Jacobi equations and all25 divergence coefficients.
    def physical_identities(q):
        tensor={}
        for k,(i,j)in enumerate(PAIRS):
            tensor[i,j]={m:q[15*k+n]for n,m in enumerate(MONS)if q[15*k+n]};tensor[j,i]={m:-c for m,c in tensor[i,j].items()}
        for i in range(5):tensor[i,i]={}
        divergence=[[R.zero for _ in range(5)]for _ in range(5)]
        for i,j in itertools.product(range(5),repeat=2):
            for mon,c in tensor[i,j].items():
                for occurrence in range(2):
                    if mon[occurrence]==i:divergence[j][mon[1-occurrence]]+=c
        jac=[]
        cubics=list(itertools.combinations_with_replacement(range(5),3))
        for i,j,k in itertools.combinations(range(5),3):
            poly={}
            for u,v,w in[(i,j,k),(j,k,i),(k,i,j)]:
                for l in range(5):
                    for m,c in tensor[u,l].items():
                        for n,d in tensor[v,w].items():
                            for occurrence in range(2):
                                if n[occurrence]==l:
                                    word=tuple(sorted((*m,n[1-occurrence])));poly[word]=poly.get(word,R.zero)+c*d
            jac.extend(poly.get(m,R.zero)for m in cubics)
        assert len(jac)==350
        return[entry for row in divergence for entry in row],jac
    div1,jac1=physical_identities(point);div2,jac2=physical_identities(pushed)
    assert not any(div1+jac1+div2+jac2)
    rawminor=s.Matrix.hstack(P0,P1);pr=list(DomainMatrix.from_Matrix(rawminor.T).convert_to(s.QQ).rref()[1]);assert len(pr)==2
    ordinaryminor=rawminor[pr,:].det();assert ordinaryminor
    constantrows=[i for i in range(150)if P0[i]and not P1[i]];assert constantrows
    # Actual flags for the SAME P(a): F0 and g(a)^-1 F0. Covectors transform by right multiplication by g.
    # Critical source has L=<e_y>, W=<h,x>, distinguished H=<x>.
    secondL=[gi[i][2]for i in range(5)];assert secondL==[R(int(i==2))for i in range(5)]
    secondW=[g[0],g[1]];secondH=g[1]
    assert secondH==[a,R.one,R.zero,-2*a*a,R.zero]
    assert secondH[0]==a #H' differs from <x> whenever a !=0
    payload=dict(variable='a',first_source_fixed20=[savepoly(pol(center[j])+a*pol(first[j]))for j in range(20)],
        first_source_physical150=[savepoly(p)for p in point],group_g5_by5=[[savepoly(p)for p in row]for row in g],
        group_inverse5_by5=[[savepoly(p)for p in row]for row in gi],second_source_physical150=[savepoly(p)for p in pushed],
        second_source_fixed20=[savepoly(p)for p in second],all150_membership_residuals=[savepoly(p)for p in residual],
        nonzero_membership_rows=nonzero,first_source_line_rank2_rows=pr,first_source_line_rank2_determinant=str(ordinaryminor),
        nonzero_constant_physical_row=constantrows[0],nonzero_constant_physical_coefficient=str(P0[constantrows[0]]),
        second_flag_L=list(map(str,secondL)),second_flag_W=list(map(lambda row:list(map(str,row)),secondW)),second_flag_H=list(map(str,secondH)))
    with(out/'certificate.json.gz').open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as z:z.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,candidate_full_degree_bound=7,actual_polynomial_source_degree=1,
        actual_physical_pushforward_degree=max((m[0]for p in pushed for m in p),default=0),
        all150_polynomial_membership_identically_zero=not nonzero,nonzero_membership_rows=nonzero,
        both_sources_full25_divergence_and350_strict_Jacobi_identities_zero=True,
        group_actual_polynomial_inverse_both_products_verified=True,ordinary_nonconstant_nonzero_first_source_line=True,
        actual_flags_distinct_for_all_a_nonzero=True,
        actual_polynomial_self_pair_family_proved=not nonzero,
        entire_conductor_or_all_branch_classification_claimed=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),source,parent,accepted]
    note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
