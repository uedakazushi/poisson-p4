#!/usr/bin/env python3
"""No Sp4-invariant equation of degree at most24 vanishes on the T02 image.

Exact characteristic-zero upper bounds come from Weyl characters. Nonzero
modular evaluation minors certify characteristic-zero independence of explicit
contractions on rational T02 jets. This proves a degree lower bound, not an
implicit equation or exclusion of unknown S61 components.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from datetime import datetime,timezone
import hashlib,json,random,time
from itertools import combinations,permutations,product
from pathlib import Path
import sympy as s
from weight1_structure import J,form_matrices,constraints
from elliptic_weight1_dominance import symplectic_frame,primitive_vector,projective_affine
from q51_universal_section_20261101 import character_certificate
from s61_model import X,jacobi_projective
ROOT=Path(__file__).resolve().parents[1]
PRIME=1000003
P8=[0,1,2,3,4]
P12=[0,1,2,3,4,10,11,17,23]
P16=[0,1,2,3,4,7,8,9,10,13,14,15,18,19,22,38,39,45,51,95]
P18=[0,1,2,22,23,30,56,68,86]


def necklaces(n):
    return [w for w in product(range(2),repeat=n) if w==min(w[i:]+w[:i] for i in range(n))]
WORDS={n:necklaces(n) for n in (4,6,8,9)}


def invariant_values(vector):
    T=form_matrices(vector)
    nz=[(i,k,J[i,k]) for i in range(4) for k in range(4) if J[i,k]]
    M=s.Matrix(4,4,lambda a,b:sum(T[a][i,j]*T[b][k,l]*v*w
        for i,k,v in nz for j,l,w in nz))
    N=s.Matrix(4,4,lambda a,b:sum(T[j][a,i]*T[l][b,k]*v*w
        for i,k,v in nz for j,l,w in nz))
    A=[J*M,J*N];i4=s.trace(A[0]**2);cache={():s.eye(4)}
    def mat(word):
        if word not in cache:cache[word]=mat(word[:-1])*A[word[-1]]
        return cache[word]
    raised=J*M*J.T
    v=J*s.Matrix([sum(T[a][b,i]*raised[a,b] for a in range(4) for b in range(4)) for i in range(4)])
    av,bv=A[0]*v,A[1]*v
    raw8=[i4*i4]+[s.trace(mat(w)) for w in WORDS[4]]
    raw12=[i4*f for f in raw8]+[s.trace(mat(w)) for w in WORDS[6]]+[
        (v.T*J*mat(w)*v)[0] for w in product(range(2),repeat=3)]
    raw14=[sum(xx[a]*(v.T*T[a]*zz)[0] for a in range(4))
           for xx in (av,bv) for zz in (av,bv)]
    raw16=[a*b for i,a in enumerate(raw8) for b in raw8[i:]]+[i4*f for f in raw12]+[
        s.trace(mat(w)) for w in WORDS[8]]+[(v.T*J*mat(w)*v)[0]
        for w in product(range(2),repeat=5)]
    raw18=[i4*f for f in raw14]+[s.trace(mat(w)) for w in WORDS[9]]+[
        (v.T*J*mat(w)*v)[0] for w in product(range(2),repeat=6)]
    g8=[raw8[j] for j in P8];g12=[raw12[j] for j in P12]
    g14=raw14[:3];g16=[raw16[j] for j in P16];g18=[raw18[j] for j in P18]
    return {4:[i4],8:g8,12:g12,14:g14,16:g16,18:g18,
        20:[a*b for a in g8 for b in g12]+[i4*f for f in g16],
        22:[i4*f for f in g18]+[a*b for a in g8 for b in g14],
        24:[a*b for a in g8 for b in g16]+[a*b for i,a in enumerate(g12) for b in g12[i:]]}


def invariant_multiplicities(maximum=24):
    weights=character_certificate()['W_weights_with_multiplicity']
    character=[defaultdict(int) for _ in range(maximum+1)];character[0][0,0]=1
    for a,b in weights:
        for n in range(1,maximum+1):
            for (x,y),value in list(character[n-1].items()):
                character[n][x+a,y+b]+=value
    dimensions={}
    for degree in range(maximum+1):
        dimension=0
        for permutation in permutations((0,1)):
            for signs in product((-1,1),repeat=2):
                rho=(2,1);wr=tuple(signs[i]*rho[permutation[i]] for i in range(2))
                sign=(1 if permutation==(0,1) else -1)*signs[0]*signs[1]
                dimension+=sign*character[degree][2-wr[0],1-wr[1]]
        dimensions[degree]=dimension
    return weights,dimensions


def modular(value):
    value=s.Rational(value)
    denominator=int(value.q)%PRIME
    assert denominator
    return int(value.p)%PRIME*pow(denominator,-1,PRIME)%PRIME


def pivot_minor(matrix):
    work=[row[:] for row in matrix];row_ids=list(range(len(work)))
    pivot_rows=[];pivot_columns=[];row=0
    for column in range(len(work[0])):
        found=next((j for j in range(row,len(work)) if work[j][column]),None)
        if found is None:continue
        work[row],work[found]=work[found],work[row]
        row_ids[row],row_ids[found]=row_ids[found],row_ids[row]
        pivot_rows.append(row_ids[row]);pivot_columns.append(column)
        inverse=pow(work[row][column],-1,PRIME)
        work[row]=[(v*inverse)%PRIME for v in work[row]]
        for j in range(row+1,len(work)):
            coefficient=work[j][column]
            if coefficient:work[j]=[(a-coefficient*b)%PRIME for a,b in zip(work[j],work[row],strict=True)]
        row+=1
        if row==len(work):break
    square=[[matrix[i][j] for j in pivot_columns] for i in pivot_rows]
    determinant=1
    for i in range(row):
        found=next(j for j in range(i,row) if square[j][i])
        if found!=i:square[found],square[i]=square[i],square[found];determinant=-determinant
        pivot=square[i][i];determinant=determinant*pivot%PRIME;inverse=pow(pivot,-1,PRIME)
        for j in range(i+1,row):
            coefficient=square[j][i]*inverse%PRIME
            square[j]=[(a-coefficient*b)%PRIME for a,b in zip(square[j],square[i],strict=True)]
    assert determinant
    return {'rank':row,'rows':pivot_rows,'columns':pivot_columns,'determinant_mod_prime':determinant%PRIME}


def verify():
    started=time.monotonic();assert s.isprime(PRIME)
    weights,dimensions=invariant_multiplicities()
    expected={4:1,8:5,12:9,14:3,16:20,18:9,20:36,22:21,24:67}
    assert {degree:n for degree,n in dimensions.items() if degree and n}==expected
    a,b,c,d=s.symbols('a b c d');parameters=[a,b,c,d]
    homogeneous=s.zeros(5);homogeneous[0,1]=-a*X[2]**2
    for (i,j),coefficient in zip(combinations(range(1,5),2),[a,a+2*b,a+2*c,b,c,d],strict=True):
        homogeneous[i,j]=coefficient*X[i]*X[j]
    homogeneous=homogeneous-homogeneous.T
    assert all(value==0 for value in jacobi_projective(homogeneous))
    y=s.symbols('y0:4');affine=projective_affine(homogeneous,X,y)
    rng=random.Random(20261101)
    # Fixed deterministic exploration seed; no random correctness claim is used.
    for _ in range(72*16):rng.randrange(-3,5)
    evaluations={degree:[] for degree in expected};samples=[]
    closure=s.Matrix.vstack(*constraints())
    while len(samples)<80:
        values=[rng.randrange(1,9) for _ in range(8)]
        substitution=dict(zip([*parameters,*y],values,strict=True));constant=affine.subs(substitution)
        if not constant.det() or not constant[0,3]:continue
        frame=symplectic_frame(constant);inverse=frame.inv()
        derivative=[affine.diff(z).subs(substitution) for z in y]
        linear=[frame*sum((inverse[h,j]*derivative[h] for h in range(4)),s.zeros(4))*frame.T
                for j in range(4)]
        jet=primitive_vector(linear)
        assert closure*jet==s.zeros(8,1)
        value=invariant_values(jet)
        for degree,row in value.items():evaluations[degree].append(list(map(modular,row)))
        samples.append({'parameters_a_b_c_d_and_point':values,'primitive_beta':list(map(str,jet))})
    certificates={degree:pivot_minor(matrix) for degree,matrix in evaluations.items()}
    assert all(certificates[degree]['rank']==dimension for degree,dimension in expected.items())
    inputs=[Path(__file__),ROOT/'code/weight1_structure.py',ROOT/'code/elliptic_weight1_dominance.py',
            ROOT/'code/q51_universal_section_20261101.py',ROOT/'code/s61_model.py',
            ROOT/'reports/MPS_T02_FIRST_JET_20261101.json']
    return {'id':'MPS_T02_INVARIANT_OBSTRUCTION_20261101','all_passed':True,
        'generated_utc':datetime.now(timezone.utc).isoformat(),'W_weights':weights,
        'invariant_dimensions_degrees0_through24':dimensions,
        'coefficient_field':'QQ; modular nonzero minors are used only for exact characteristic-zero rank lower bounds',
        'prime':PRIME,'sample_count':len(samples),'samples':samples,
        'evaluation_matrices_mod_prime':evaluations,'rank_certificates':certificates,
        'family_universal_projective_Jacobi_checked':True,
        'family_is_MPS_T02_dependency':'MPS_T02_FIRST_JET_20261101; this is its deformation coefficient b_smooth=-a subfamily',
        'invariant_equation_degree_lower_bound':26,
        'homogeneity_and_invariance_argument':'The T02 image is a conical Sp4-stable irreducible hypersurface; its irreducible equation is homogeneous and Sp4-invariant because Sp4 has no nontrivial characters.',
        'scope':'No nonzero invariant of positive degree<=24 vanishes on the T02 first-jet image. Its defining equation is not computed. No full S61 fiber exclusion or exhaustion follows.',
        'proof_note':'research/notes/MPS_T02_INVARIANT_OBSTRUCTION_20261101.md',
        'elapsed_seconds':round(time.monotonic()-started,3),
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as handle:json.dump(report,handle,indent=2);handle.write('\n')
    print(json.dumps({'all_passed':True,'degree_lower_bound':26,'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
