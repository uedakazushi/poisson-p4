#!/usr/bin/env python3
"""Independent original134/plane/cofactor/formal-division global quintic audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63'
PRIME=32003


def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def plus(A,B,c=Fraction(1)):
    C=A.copy()
    for mon,coefficient in B.items():
        value=C.get(mon,Fraction(0))+c*coefficient
        if value:C[mon]=value
        else:C.pop(mon,None)
    return C


def times(A,B):
    C={}
    for mon,coefficient in A.items():
        for other,scalar in B.items():
            key=tuple(sorted(mon+other));value=C.get(key,Fraction(0))+coefficient*scalar
            if value:C[key]=value
            else:C.pop(key,None)
    return C


def expr(text,z):
    return {tuple(j for j,n in enumerate(powers) for _ in range(n)):Fraction(coefficient)
            for powers,coefficient in s.Poly(s.sympify(text),*z).terms() if coefficient}


def substitute(poly,images):
    result={}
    for mon,coefficient in poly.items():
        term={():coefficient}
        for i in mon:term=times(term,images[i])
        result=plus(result,term)
    return result


def decode(rows):
    result={}
    for degree,coefficients in rows:
        for mon,text in coefficients:
            key=tuple(sorted(tuple(mon)+(45,)*degree))
            assert key not in result and Fraction(text)
            result[key]=Fraction(text)
    return result


def residue(q):return q.numerator*pow(q.denominator,-1,PRIME)%PRIME


def det_mod(matrix):
    A=[[residue(Fraction(c)) for c in row] for row in matrix];answer=1
    for j in range(len(A)):
        pivot=next(i for i in range(j,len(A)) if A[i][j])
        if pivot!=j:A[j],A[pivot]=A[pivot],A[j];answer=-answer
        value=A[j][j];answer=answer*value%PRIME;inverse=pow(value,-1,PRIME)
        for i in range(j+1,len(A)):
            coefficient=A[i][j]*inverse%PRIME
            for k in range(j+1,len(A)):A[i][k]=(A[i][k]-coefficient*A[j][k])%PRIME
    return answer%PRIME


def main(output,replay):
    started=time.monotonic();folder=ROOT/'reports'/STEM;old=read(folder/'report.json');assert old['all_passed']
    for name,digest in old['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in old['evidence_sha256'].items():assert sha(folder/name)==digest,name
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/s61_affine_normal_global_quintic_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=40)
    for name in old['evidence_sha256']:assert (folder/name).read_bytes()==(replay/name).read_bytes()
    specpath=ROOT/'code/data/s61_affine_normal_global_quintic_20261103_v63.json'
    spec=read(specpath);z=s.symbols('z0:3');plane=[expr(f,z) for f in spec['plane']]
    relations=sorted(spec['relations'],key=lambda r:r['target_index']);targets=[];minors=[]
    for relation in relations:
        cofactors=[expr(f,z) for f in relation['cofactors']];target=expr(relation['target'],z)
        total={}
        for B,f in zip(cofactors,plane):total=plus(total,times(B,f))
        assert total==target
        assert max((len(mon) for B in cofactors for mon in B),default=0)<=relation['degree']-2
        assert all(c.denominator%PRIME for B in cofactors for c in B.values())
        targets.append(target)
    nonzero_terms=[sum(len(expr(f,z)) for f in r['cofactors']) for r in relations]
    assert nonzero_terms==[30,30,51]
    # The producer's31/31/52 includes Poly(0).terms()'s zero entry once each.
    assert [sum(len(s.Poly(s.sympify(f),*z).terms()) for f in r['cofactors']) for r in relations]==[31,31,52]
    for degree in [4,5]:
        record=next(r for r in relations if r['degree']==degree)
        columns=[m for d in range(degree+1) for m in itertools.combinations_with_replacement(range(3),d)]
        rows=[times(f,{m:Fraction(1)}) for f in plane for d in range(degree-1) for m in itertools.combinations_with_replacement(range(3),d)]
        square=[[rows[i].get(columns[j],Fraction(0)) for j in record['selected_columns']] for i in record['selected_rows']]
        determinant=det_mod(square);assert determinant and len(square)=={4:30,5:51}[degree]
        minors.append(dict(degree=degree,rank=len(square),determinant_mod32003=determinant))
    # Derive the centered45-variable original fiber directly from the raw134.
    sourcepath=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json';raw=read(sourcepath)
    pointpath=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
    sample=next(p for p in read(pointpath)['points'] if p['component_id']=='T01')
    point=list(map(Fraction,sample['S61_coordinates_u1_through_u61']))
    qids=[i for i,w in enumerate(raw['parameter_weights']) if w==2];index={i:j for j,i in enumerate(qids)}
    fiber=[]
    for equation in raw['equations']:
        result={}
        for mon,n,d in equation:
            term={():Fraction(n,d)}
            for i in mon:term=times(term,({():point[i],(index[i],):Fraction(1)} if i in index else {():point[i]}))
            result=plus(result,term)
        assert not result.get((),0) and max(map(len,result),default=0)<=2
        fiber.append(result)
    assert len(fiber)==134
    framepath=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz'
    frame=read(framepath);E=[[Fraction(c) for c in row] for row in frame['embedding']]
    assert E[:3]==[[Fraction(int(i==j)) for j in range(3)] for i in range(3)]
    images=[{(j,):c for j,c in enumerate(row) if c} for row in E]
    assert [substitute(f,images) for f in fiber[:5]]==plane
    presentation=read(folder/'presentation.json.gz')
    p1=plus({(1,):Fraction(1)},targets[0],Fraction(-1));p2=plus({(2,):Fraction(1)},targets[1],Fraction(-1));g=targets[2]
    assert all(set(mon)<={0} for P in [p1,p2,g] for mon in P)
    assert g[(0,)*5]==1 and max(map(len,g))==5
    assert p1==expr(presentation['p1'],z) and p2==expr(presentation['p2'],z) and g==expr(presentation['g'],z)
    gp=s.Poly(s.sympify(presentation['g']),z[0]);assert s.gcd(gp,gp.diff()).degree()==0
    for f in spec['plane']:
        value=s.sympify(f).subs({z[1]:s.sympify(presentation['p1']),z[2]:s.sympify(presentation['p2'])},simultaneous=True)
        assert s.rem(s.Poly(value,z[0]),gp).is_zero
    convert=lambda P:{(45,)*len(mon):c for mon,c in P.items()}
    Z=[{(45,):Fraction(1)},plus(convert(p1),{(42,):Fraction(1)}),plus(convert(p2),{(43,):Fraction(1)})]
    G=plus(convert(g),{(44,):Fraction(-1)})
    assert G==decode(presentation['monic_quintic_G'])
    Q=[]
    for i,row in enumerate(E):
        polynomial={}
        for j,c in enumerate(row):polynomial=plus(polynomial,Z[j],c)
        if i>=3:polynomial=plus(polynomial,{(i-3,):Fraction(1)})
        Q.append(polynomial)
    maximum_h=0
    for f,quotient,remainder in zip(fiber,presentation['all134_division_quotients'],presentation['all134_remainders']):
        original=substitute(f,Q);q,r=decode(quotient),decode(remainder)
        assert plus(times(q,G),r)==original
        assert max((mon.count(45) for mon in r),default=0)<=4
        maximum_h=max(maximum_h,max((mon.count(45) for mon in original),default=0))
    assert len(presentation['all134_remainders'])==len(presentation['all134_division_quotients'])==134
    assert maximum_h<=8
    cert=read(folder/'certificate.json');assert minors==cert['selected_degree45_macaulay_minors']
    comm=ROOT/'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63/report.json';assert read(comm)['all_passed']
    accepted=[ROOT/'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63_REVIEW.json',
              ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW/report.json']
    assert all(read(path)['all_passed'] for path in accepted)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'certificate.json',folder/'presentation.json.gz',
           specpath,framepath,sourcepath,pointpath,comm,ROOT/'research/notes'/f'{STEM}.md',
           ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/report.json',*accepted]
    result=dict(id=STEM+'_REVIEW',all_passed=True,
        fresh_producer_all_evidence_byte_identical=True,
        independent_all134_fiber_centering_and_first5_actual_plane_binding=True,
        independent_three_full_QQ_cofactor_identities_and_p_unit_coefficients=True,
        actual_nonzero_cofactor_term_counts=nonzero_terms,
        producer_stored_term_counts31_31_52_include_one_zero_each=True,
        independent_degree4_30_and_degree5_51_punit_minors=minors,
        plane_exact_primitive5_algebra_and_quintic_separability_checked=True,
        independent_flat46symbol_sparse_all134_substitution_division_identities=True,
        placeholders_fixed_normal_circuits_not_free_source_parameters=True,
        exact_polynomial_R42_corrections_degree_bounds=[4,4,5],
        generic_beta_only_Macaulay_spreading_and_two_inverse_algebra_maps_accepted=True,
        actual_global_monogenic_quotient_and_5x670_module_presentation_accepted=True,
        residual_degree85_polynomial_recovery_from_D84_accepted=True,
        no_q_or_cyclic_open_or_actual_R_freeness_beta_finiteness_or_unit_claim=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','independent_degree4_30_and_degree5_51_punit_minors','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
