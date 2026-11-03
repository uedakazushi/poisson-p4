#!/usr/bin/env python3
"""Independent full-polynomial canonical family, rank26 and component25 audit."""
import argparse
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix

from q_positive_weight6_small_model_review_20261103_v63 import physical_coefficients,own_action,read,qm,mm,rank,limits,PAIRS,MON2,BI

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    producer=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63';report=read(producer/'report.json');assert report['all_passed']
    for path,digest in report['input_sha256'].items():assert sha(ROOT/path)==digest
    assert sha(producer/'certificate.json.gz')==report['certificate_sha256']
    public=read(producer/'certificate.json.gz')
    parentpath=ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/certificate.json.gz';parent=read(parentpath)
    quadratic=ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63';old=read(quadratic/'certificate.json.gz')
    R,*ns=ring('n0,n1,n2',s.QQ);zero=R.zero
    parse=lambda value:R.from_expr(s.sympify(value))
    T=qm(parent['normal99_by3']);graph=[sum((R.from_expr(T[i,j])*ns[j]for j in range(3)),zero)for i in range(99)]
    graph[17]+=ns[1]*ns[2]/48;graph[49]-=ns[1]**2/1536;graph[54]+=ns[0]*ns[2]/2
    assert all(a==parse(b)for a,b in zip(graph,parent['exact99_polynomial_graph']))
    raw=[R.from_expr(s.Rational(v))for v in old['physical_center150']]
    affine_mons=[(),(0,),(1,),(2,),(3,),*itertools.combinations_with_replacement(range(4),2)]
    affine_pairs=list(itertools.combinations(range(4),2))
    for k,(i,j)in enumerate(affine_pairs):
        for m,mon in enumerate(affine_mons):
            homogeneous=tuple([0]*(2-len(mon))+[v+1 for v in mon])
            raw[BI[(i+1,j+1),homogeneous]]+=graph[15*k+m]
    pairtable,D=physical_coefficients();div=[]
    for row in D.tolist():div.append(sum((int(c)*v for c,v in zip(row,raw)if c),zero))
    canonical=list(raw)
    for i,j in PAIRS:
        for k in range(5):
            canonical[BI[(i,j),tuple(sorted((i,k)))]]-=div[5*j+k]/5
            canonical[BI[(i,j),tuple(sorted((j,k)))]]+=div[5*i+k]/5
    assert all(sum((int(c)*v for c,v in zip(row,canonical)if c),zero)==zero for row in D.tolist())
    assert all(poly==parse(saved)for poly,saved in zip(canonical,public['canonical150']))
    assert all(sum(exponent)<=2 for poly in canonical for exponent in poly)
    # Direct affine4 Jacobi engine, coefficient ring QQ[n0,n1,n2].
    tensor={}
    for i,j in affine_pairs:
        poly={}
        for mon in MON2:
            c=raw[BI[(i+1,j+1),mon]]
            if c:
                affine=tuple(v-1 for v in mon if v);poly[affine]=poly.get(affine,zero)+c
        tensor[i,j]={mon:c for mon,c in poly.items()if c};tensor[j,i]={mon:-c for mon,c in poly.items()if c}
    for i in range(4):tensor[i,i]={}
    for i,j,k in itertools.combinations(range(4),3):
        result={}
        for a,b,c in[(i,j,k),(j,k,i),(k,i,j)]:
            for ell in range(4):
                for mon,left in tensor[a,ell].items():
                    for other,right in tensor[b,c].items():
                        if ell not in other:continue
                        remaining=list(other);remaining.remove(ell);word=tuple(sorted(mon+tuple(remaining)))
                        result[word]=result.get(word,zero)+other.count(ell)*left*right
        assert all(not value for value in result.values())
    # Native STRICT cubic Jacobi350, full parameter degree<=4, not the
    # producer's Euler-projective quartic equations and not truncated jets.
    native=ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.cpp';binary=ROOT/'work/q-new-family-review-native'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    subprocess.run(['g++','-O2','-std=c++17',str(native),*flags,'-o',str(binary)],check=True,timeout=30)
    inputpath=args.output/'full-polynomial-physical-input.txt'
    with inputpath.open('w')as stream:
        for degree in range(5):
            for poly in canonical:
                terms=[(exponent,c)for exponent,c in poly.items()if sum(exponent)==degree]
                stream.write(str(len(terms))+'\n')
                for exponent,c in terms:stream.write(str(c)+' '+' '.join(map(str,exponent+(0,)*5))+'\n')
        for row in pairtable:
            stream.write(str(len(row))+'\n')
            for(a,b),c in sorted(row.items()):stream.write(f'{a} {b} {c}\n')
        for row in D.tolist():
            terms=[(i,int(c))for i,c in enumerate(row)if c];stream.write(str(len(terms))+'\n')
            for i,c in terms:stream.write(f'{i} {c}\n')
        for _ in range(5*350):stream.write('0\n')
    with(args.output/'full-polynomial-physical.log').open('w')as log:
        subprocess.run([str(binary),str(inputpath)],stdout=log,stderr=subprocess.STDOUT,timeout=60,preexec_fn=limits,check=True)
    assert 'ALL25_DIVERGENCE_ALL5GRADES_PASS'in(args.output/'full-polynomial-physical.log').read_text()
    def evaluate(poly,values):
        return sum((s.Rational(c)*s.prod(values[i]**e for i,e in enumerate(exponent))for exponent,c in poly.items()),s.S.Zero)
    rank_certificates={};ambient=[]
    C=qm(read(ROOT/old['canonical_basis_dependency'])['canonical_basis'])
    for point in[(0,0,0),(1,1,1),(1,2,3),(2,3,5)]:
        q=s.Matrix([evaluate(poly,point)for poly in canonical]);action=own_action(q)
        tangent=s.Matrix([[evaluate(poly.diff(ns[j]),point)for j in range(3)]for poly in canonical])
        M=action.row_join(tangent)
        if point in[(0,0,0),(1,1,1)]:
            expected=public['rank_certificates'][str(point)]
            assert M==qm(expected['full150_by29_matrix'])and rank(M)==26
            minor=M[expected['minor_rows'],expected['minor_columns']]
            det=DomainMatrix.from_Matrix(minor).convert_to(s.QQ).det()
            assert det==s.Rational(expected['determinant'])!=0
            rank_certificates[str(point)]=dict(rank=26,determinant=str(det),rows=expected['minor_rows'],columns=expected['minor_columns'])
        if point!=(0,0,0):
            L=s.zeros(350,150)
            for row,poly in enumerate(pairtable):
                for(a,b),c in poly.items():L[row,a]+=c*q[b];L[row,b]+=c*q[a]
            rr=rank(mm(L,C));assert rr==98
            ambient.append(dict(parameters=point,actual_strict350_canonical126_rank=rr,affine_tangent_dimension=126-rr,projective_tangent_dimension=125-rr))
    assert s.Matrix([evaluate(poly,[0,0,0])for poly in canonical])==s.Matrix([s.Rational(v)for v in old['physical_center150']])
    # Independent full quintic from the Euler contraction of four-Pfaffians.
    P={pair:{mon:canonical[BI[pair,mon]]for mon in MON2 if canonical[BI[pair,mon]]}for pair in PAIRS}
    def product(a,b):
        out={}
        for mon,c in a.items():
            for word,d in b.items():out[tuple(sorted(mon+word))]=out.get(tuple(sorted(mon+word)),zero)+c*d
        return out
    quintic={}
    for i in range(5):
        a,b,c,d=[j for j in range(5)if j!=i]
        for sign,left,right in[(1,(a,b),(c,d)),(-1,(a,c),(b,d)),(1,(a,d),(b,c))]:
            for mon,value in product(P[left],P[right]).items():
                word=tuple(sorted(mon+(i,)));quintic[word]=quintic.get(word,zero)+(-1)**i*sign*value
    quintic={mon:c for mon,c in quintic.items()if c};assert quintic and all(mon.count(0)>=3 for mon in quintic)
    Qmatrix=[[zero for _ in range(5)]for _ in range(5)]
    for mon,value in quintic.items():
        word=list(mon)
        for _ in range(3):word.remove(0)
        i,j=word;Qmatrix[i][j]+=value/(1 if i==j else 2)
        if i!=j:Qmatrix[j][i]+=value/2
    symbolic=s.Matrix([[value.as_expr()for value in row]for row in Qmatrix])
    assert symbolic==s.Matrix([[s.sympify(v)for v in row]for row in public['Q_symmetric_matrix']])
    assert s.expand(symbolic[1:,1:].det(method='domain-ge'))==0
    sample=symbolic.subs(dict(zip(s.symbols('n0:3'),[1,1,1])))
    assert sample.det()!=0 and sample[1:,1:].rank()==3
    cataloguepath=ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json'
    catalogue=read(cataloguepath);records=catalogue['records'];assert len(records)==51
    large={r['component_id']:r['expected_s61_dimension']+9 for r in records if r['expected_s61_dimension']+9>=25}
    assert large=={'T01':25,'T39':25,'T40':25,'R4_CA':28,'R4_CB':28}
    proof_inputs=[ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63_REVIEW.json',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW/report.json']
    for path in proof_inputs:
        saved=read(path);assert saved['all_passed']
        for p,digest in saved.get('input_sha256',{}).items():assert sha(ROOT/p)==digest
    (args.output/'component-certificate.json').write_text(json.dumps(dict(
        actual_full_three_parameter_family=True,all140_affine_Jacobi_exact_QQ=True,
        all350_STRICT_cubic_canonical_Jacobi_exact_full_polynomial=True,
        canonical150_independently_reconstructed_from_old99graph=True,
        auxiliary_N_squared_equation_not_imposed=True,rank26_certificates=rank_certificates,
        irreducible_actual_projective_image_dimension_exact=25,
        full_Poisson_irreducible_component_proved=True,
        component_dimension_upper_from_full_formal_initial4_5=True,
        generic_divisor='3H+Q, Q smooth and H tangent',
        hyperplane_section_determinant_identically_zero_over_QQn3=True,
        known51_only_components_of_dimension_atleast25=large,
        known_component_dimension_uppers_reused_from_audited_legacy_catalogue=True,
        known_six_rank2_closed_loci_cannot_contain_generic_rank4=True,
        actual_new_component_distinct_from_prior57=True,proved_global_component_lower_bound=58,
        actual_full_Jacobian_point_checks=ambient,
        pointwise_tangent_excess_not_used_as_nilpotence_or_generic_nonreducedness=True,
        firstjet_dominance_or_generic_five_point_exhaustion_claimed=False),indent=2)+'\n')
    paths=[Path(__file__),native,producer/'report.json',producer/'certificate.json.gz',parentpath,
        cataloguepath,*proof_inputs,ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.py',
        ROOT/'research/notes/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW.md']
    result=dict(id='Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW',all_passed=True,
        independent_full_polynomial_affine140_and_strict350_and_canonical150=True,
        independent_rank26_and_exact_component_dimension25=True,
        exact_irreducible_new_component_distinct_from_prior57=True,
        global_component_lower_bound58_proved=True,generic_smoothness_or_reducedness_not_proved=True,
        ambient_canonical_Jacobian_rank_at111_123_235=[98,98,98],
        firstjet_dominance_or_full_generic_fiber_exhaustion_not_proved=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('INDEPENDENT ACTUAL NEW25 COMPONENT / LOWER58 PASS',result['seconds'],flush=True)


if __name__=='__main__':main()
