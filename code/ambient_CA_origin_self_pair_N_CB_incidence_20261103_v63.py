#!/usr/bin/env python3
"""Actual N and CB source coordinates for the CA ordinary polynomial self-pair."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63'
X=s.symbols('h x y z w');h,x,y,z,w=X;a=s.symbols('a',nonzero=True)
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[m]*X[n])for i,j in PAIRS for m,n in MONS])
def canonical(P):
    div=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)]);E=s.Matrix(X)
    return(P-(E*div.T-div*E.T)/5).applyfunc(s.expand)
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(60)
    parent=ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63';assert read(parent/'report.json')['all_passed']
    reviewed=ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW.json';assert read(reviewed)['all_passed']
    ca=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';A=read(ca)
    np=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';N=read(np)
    nd=ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json';assert read(nd)['all_passed']
    boundary=ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json';assert read(boundary)['all_passed']
    data=read(parent/'certificate.json.gz');first=s.Matrix([sum(s.Rational(c)*a**n for n,c in pol)for pol in data['first_source_fixed20']])
    second=s.Matrix([sum(s.Rational(c)*a**n for n,c in pol)for pol in data['second_source_fixed20']])
    C=s.Matrix(A['fixed_fiber_linear_basis150_by20']);NC=s.Matrix(N['fiber19_basis']);indices=[0,1,*range(3,20)]
    assert NC==C[:,indices]and first[2]==0 and second[2]==-a*a
    ncoords=first[indices,:];q=NC*ncoords
    physical=s.Matrix([sum(s.Rational(c)*a**n for n,c in pol)for pol in data['first_source_physical150']]);assert q==physical
    P=s.zeros(5)
    for k,(i,j)in enumerate(PAIRS):P[i,j]=sum(q[15*k+m]*X[u]*X[v]for m,(u,v)in enumerate(MONS));P[j,i]=-P[i,j]
    # Explicit full CB fixed-source after the square-boundary change of variables.
    # In the un-swapped CA source, f=h*(-h-3a*x), H=-a*h*w*z-h*w*x+x^2*z+z^3/3.
    G=s.eye(5);G[0,0]=-3*a;G[0,1]=-1
    Gi=G.inv();assert (G*Gi)==s.eye(5)and(Gi*G)==s.eye(5)and G.det()==-3*a
    pushed=(G*P.subs(dict(zip(X,Gi*s.Matrix(X))),simultaneous=True)*G.T).applyfunc(s.cancel)
    HB=-a*x*w*z+x*w*(h+x)/(3*a)+(h+x)**2*z/(9*a*a)+z**3/3;fB=3*a*x*x
    rawB=s.zeros(5);rawB[1,2]=fB;rawB[2,3]=-s.diff(fB,x)*z/2-s.diff(HB,w);rawB[2,4]=-s.diff(fB,x)*w/2+s.diff(HB,z);rawB=rawB-rawB.T
    expected=canonical(rawB);difference=(pushed-expected).applyfunc(s.cancel);assert difference==s.zeros(5)
    assert vector(pushed)==vector(expected)and len(vector(pushed))==150
    # The canonical equality is exact on D(a), not a truncation at a=0.
    div=s.Matrix([sum(s.diff(pushed[i,j],X[i])for i in range(5))for j in range(5)]);assert div==s.zeros(5,1)
    jac=[]
    for i,j,k in itertools.combinations(range(5),3):
        J=s.cancel(sum(pushed[i,l]*s.diff(pushed[j,k],X[l])+pushed[j,l]*s.diff(pushed[k,i],X[l])+pushed[k,l]*s.diff(pushed[i,j],X[l])for l in range(5)))
        assert J==0;jac.append(str(J))
    payload=dict(N_fixed19_CA_indices=indices,first_source_N19_coordinates=list(map(str,ncoords)),all150_N_source_coordinates_match=True,
        first_CA_lift_missing_C2='0',second_CA_lift_missing_C2=str(second[2]),
        actual_CB_change_of_coordinates=[[str(c)for c in row]for row in G.tolist()],actual_CB_inverse=[[str(c)for c in row]for row in Gi.tolist()],
        actual_CB_binary_f=str(fB),actual_CB_g='0',actual_CB_cubic_H=str(s.expand(HB)),
        all150_changed_P=list(map(str,vector(pushed))),all150_CB_source=list(map(str,vector(expected))),
        all150_difference=list(map(str,vector(difference))),strict_cubic_Jacobi_ten_polynomials=jac)
    with(out/'certificate.json.gz').open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,actual_N_source_membership_all_a=True,
        first_CA_flag_is_N_fixed_source=True,second_CA_flag_is_not_N_fixed_source_for_a_nonzero=True,
        actual_CB_source_membership_explicit_on_Da=True,actual_CB_membership_at_a0_by_existing_P0_and_closedness=True,
        actual_CB_coordinate_determinant=str(G.det()),all150_canonical_CB_source_identities_zero=True,
        two_N_flags_or_R19_membership_claimed=False,contradiction_with_N_subintegral_or_CA_smooth_RminusN_claimed=False,
        seconds=time.monotonic()-start)
    paths=[Path(__file__),parent/'report.json',parent/'certificate.json.gz',reviewed,ca,np,nd,boundary];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
