#!/usr/bin/env python3
"""Independent degree16 basis-restriction and actual-point audit."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61, X, pushforward
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63_REVIEW'
P=1000003


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes())


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def matrix(rows):return s.Matrix([[s.Rational(v)for v in row]for row in rows])
def mod(v):
    q=Fraction(str(v));return q.numerator*pow(q.denominator,-1,P)%P
def modmat(M):return np.array([[mod(v)for v in row]for row in M.tolist()],dtype=np.int64)


def own_dag(payload,point):
    values=[]
    for node in payload['nodes']:
        op=node[0]
        if op=='q':v=Fraction(node[1])
        elif op=='v':v=Fraction(point[payload['input_names'][node[1]]])
        elif op=='+':v=sum((values[j]for j in node[1]),Fraction(0))
        elif op=='*':v=math.prod(values[j]for j in node[1])
        elif op=='inv':v=1/values[node[1]]
        else:raise AssertionError(op)
        values.append(v)
    return values


def grouped_character(weights):
    # Expand each DISTINCT weight with its multiplicity using binomial coefficients.
    # This is independent of the author's in-place, individual-weight coin recursion.
    maximum=18
    chars=[defaultdict(int)for _ in range(maximum+1)];chars[0][0,0]=1
    for (a,b),m in sorted(Counter(weights).items()):
        updated=[defaultdict(int)for _ in range(maximum+1)]
        for degree,old in enumerate(chars):
            for k in range(maximum-degree+1):
                coefficient=math.comb(m+k-1,k)
                target=updated[degree+k]
                for (u,v),count in old.items():target[u+k*a,v+k*b]+=coefficient*count
        chars=updated
    denominator=Counter({(0,0):1})
    for a,b in [(2,0),(0,2),(1,1),(1,-1)]:
        updated=denominator.copy()
        for (u,v),c in denominator.items():updated[u+a,v+b]-=c
        denominator=Counter({w:c for w,c in updated.items()if c})
    for degree,char in enumerate(chars):assert sum(char.values())==math.comb(15+degree,degree)
    dims=[sum(c*char[w]for w,c in denominator.items())for char in chars]
    assert dims[:15]==[1,0,0,0,1,0,0,0,5,0,0,0,9,0,3]
    assert dims[15:]==[0,20,0,9]
    return chars,dims


def main(out,replay):
    start=time.monotonic()
    folder=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63'
    probe=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
    accepted=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json'
    mapreview=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.json'
    imagereview=ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/report.json'
    records=[read(p)for p in [folder/'report.json',probe,accepted,mapreview,imagereview]]
    checked={}
    for path,d in zip([folder/'report.json',probe,accepted,mapreview,imagereview],records):
        assert d['all_passed']
        for name,want in d.get('input_sha256',{}).items():
            source=ROOT/name;assert sha(source)==want,(path,name);checked[str(source)]=want
        for name,want in d.get('evidence_sha256',{}).items():
            source=path.parent/name;assert sha(source)==want,(path,name);checked[str(source)]=want
    author,pr,old,_,_=records
    assert old['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
    assert old['defining_hypersurface_degree_lower_bound']==16
    for name in ['native-input.txt','exact-evaluation.json','samples.json','native.log']:
        assert (folder/name).read_bytes()==(replay/name).read_bytes(),name
    lowpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    low=read(lowpath);B=matrix(low['primitive24_by16_basis']);basisrows=list(B.T.rref()[1]);left=B[basisrows,:].inv()
    O=matrix(low['symplectic_form16']);rep=list(map(matrix,low['all_ten_representation_matrices']))
    gs=sp_basis();G=s.Matrix([[s.trace(A*C)for C in gs]for A in gs]);Gi=G.inv();Oi=O.inv()
    chars,dims=grouped_character(list(map(tuple,low['Cartan_weights'])))
    for item in pr['character_records']:
        degree=item['degree'];assert dims[degree]==item['invariant_dimension']
        for term in item['Weyl_terms']:assert chars[degree][tuple(term['weight'])]==term['multiplicity']
    samples=read(folder/'samples.json');exact=read(folder/'exact-evaluation.json')
    definitions=samples['invariant_definitions'];columns=samples['selected_columns']
    words8=sorted({min(w[i:]+w[:i]for i in range(8))for w in itertools.product(range(2),repeat=8)})
    assert len(words8)==36
    assert [tuple(d['word'])for d in definitions[:36]]==words8
    assert len(columns)==20 and len(set(columns))==20
    assert all(definitions[j]['kind']=='trace8'for j in columns)
    d8=read(ROOT/'reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/certificate.json')
    basis8=d8['independent_five_columns']
    payload=read(ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz')
    family=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz')
    oldX=s.symbols('h x y z w')
    physical=s.Matrix([[s.sympify(v).subs(dict(zip(oldX,X)),simultaneous=True)for v in row]for row in family['full_canonical_homogeneous_P']])
    model=S61();z,u,beta=primitive_beta(model);ones=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    beta24=[beta[i,j].coeff(t)for i,j in itertools.combinations(range(4),2)for t in z]
    betaU=s.Matrix([[s.expand(c).coeff(v)for v in u]for c in beta24])
    betas=[]
    for index,sample in enumerate(samples['actual_QQ_points_and_betas']):
        point=sample['point'];values=own_dag(payload,point)
        U=s.Matrix([s.Rational(values[i])for i in payload['output_u1_through_u61']])
        pivot=values[payload['normalization']['value_pivot_node']];pf=values[payload['normalization']['value_pfaffian_node']]
        assert pivot and pf and str(pivot)==sample['value_pivot'] and str(pf)==sample['value_pfaffian']
        linear=s.Matrix([[s.Rational(values[j])for j in row]for row in payload['normalization']['linear_normalizing_matrix_nodes']])
        translation=s.eye(5)
        for i in range(1,5):translation[i,0]=1-point['base'+str(i)]
        actual=physical.subs({s.Symbol(name):point[name]for name in ['n0','n1','n2']})
        normalized=model.normalize(pushforward(pushforward(actual,translation),linear))
        assert list(U)==normalized['parameters']
        b24=betaU*U[ones,:];b=left*b24[basisrows,:];assert B*b==b24
        assert list(b)==list(map(s.Rational,sample['beta16']));betas.append(b)
    # Independently bind EVERY rational scalar in the native input.
    tokens=iter((folder/'native-input.txt').read_text().split())
    assert next(tokens)=='20'
    def input_matrix(expected):
        n,m=int(next(tokens)),int(next(tokens));assert (n,m)==expected.shape
        for i in range(n):
            for j in range(m):assert s.Rational(next(tokens))==expected[i,j]
    for M in [O,Oi,Gi,*rep]:input_matrix(M)
    assert [int(next(tokens))for _ in range(5)]==basis8
    assert int(next(tokens))==len(definitions)==60
    for d in definitions:
        kind=d['kind'];assert next(tokens)==('P'if kind=='basis8_product'else'T'if kind=='trace8'else'F')
        if kind=='basis8_product':assert [int(next(tokens)),int(next(tokens))]==d['indices']
        else:assert next(tokens)==''.join(map(str,d['word']))
    assert [int(next(tokens))for _ in range(20)]==columns
    for b in betas:input_matrix(b)
    assert next(tokens,None) is None
    # Different prime, different evaluator: all60 exact QQ outputs are checked.
    om,oi,gi=map(modmat,[O,Oi,Gi]);reps=list(map(modmat,rep))
    words4=sorted({min(w[i:]+w[:i]for i in range(4))for w in itertools.product(range(2),repeat=4)})
    E=[]
    for b,want in zip(betas,exact['evaluation']):
        bb=modmat(b);V=np.hstack([om@R@bb%P for R in reps])
        mu=(bb.T@V%P).T*pow(2,-1,P)%P;co=gi@mu%P
        theta=sum((int(c)*R for c,R in zip(co[:,0],reps)),np.zeros((16,16),dtype=np.int64))%P
        H=(theta+(oi@V%P@gi%P)@V.T)%P
        cache={():np.eye(16,dtype=np.int64)}
        def product(word):
            word=tuple(word)
            if word not in cache:cache[word]=product(word[:-1])@[theta,H][word[-1]]%P
            return cache[word]
        def trace(word):return int(np.trace(product(word))%P)
        f4=int((mu.T@co)[0,0]%P)
        e8=[trace(w)for w in words4]+[f4*f4%P]
        e8=[e8[j]for j in basis8];row=[]
        for d in definitions:
            if d['kind']=='trace8':v=trace(d['word'])
            elif d['kind']=='F4_times_trace6':v=f4*trace(d['word'])%P
            else:v=e8[d['indices'][0]]*e8[d['indices'][1]]%P
            row.append(v)
        assert row==list(map(mod,want));E.append(row)
    exact_minor=s.Matrix([[s.Rational(exact['evaluation'][i][j])for j in columns]for i in range(20)])
    determinant=s.Rational(DomainMatrix.from_Matrix(exact_minor).convert_to(s.QQ).det())
    assert determinant==s.Rational(exact['exact_selected20_minor'])!=0
    assert mod(determinant)!=0
    note=ROOT/f'research/notes/{STEM}.md'
    paths=[Path(__file__),note,folder/'report.json',probe,accepted,mapreview,imagereview,lowpath,
           ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz',
           ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    paths += [folder/name for name in ['samples.json','exact-evaluation.json','native-input.txt','native.log']]
    report=dict(id=STEM,all_passed=True,
        checked_unique_input_and_evidence_hashes=len(checked),
        fresh_QQ_native_all4_evidence_files_byte_equal=True,
        independent_all20_actual_QQ_source_jets_direct_full150_normalization_and_Fraction_DAG=True,
        all_native_input_matrices_definitions_and20_beta_vectors_bound=True,
        independent_all1200_values_checked_at_different_prime=P,
        independent_QQ20minor_exactly_matches_native=True,
        exact_nonzero20minor=str(determinant),
        complete_integer_characters_via_grouped_weight_binomial_factors=True,
        independent_positive_root_denominator=True,
        all_invariant_dimensions_through18=dims,
        full_degree16_invariant_dimension=20,full_degree17_invariant_dimension=0,
        degree16_restriction_rank_on_actual_new_image=20,
        universal_conjugation_is_prior_independent_proof_not_sample_inference=True,
        irreducible_invariant_cone_hypersurface_generator_invariant_by_Sp4_no_characters=True,
        actual_defining_equation_degree_lower_bound=18,
        degree18_restriction_kernel_or_equation_claimed=False,
        sample_vanishing_extrapolated_to_identity=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        seconds=time.monotonic()-start)
    with out.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print('INDEPENDENT DEGREE16 RESTRICTION INJECTIVE / EQUATION DEGREE>=18 PASS',report['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--replay',type=Path,required=True);args=p.parse_args();main(args.output,args.replay)
