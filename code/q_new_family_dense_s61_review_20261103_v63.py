#!/usr/bin/env python3
"""Independent finite-source rebuild, denseQQ AD and full61 normalizer review."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61,X,bivector_from_coefficients,pushforward,coefficients
from mps_s61_arithmetic_circuits_20261103_v62 import construct
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha


def dual(nodes,point,names):
    dimension=len(names);values=[];derivatives=[]
    for i,node in enumerate(nodes):
        op,data=node
        if op=='q':value=Fraction(data);gradient=[Fraction(0)]*dimension
        elif op=='v':
            assert 0<=data<dimension
            value=Fraction(point[names[data]]);gradient=[Fraction(int(j==data))for j in range(dimension)]
        elif op=='+':
            assert all(j<i for j in data)
            value=sum((values[j]for j in data),Fraction(0))
            gradient=[sum((derivatives[j][k]for j in data),Fraction(0))for k in range(dimension)]
        elif op=='*':
            assert all(j<i for j in data)
            value=Fraction(1);gradient=[Fraction(0)]*dimension
            for j in data:
                gradient=[value*d+values[j]*old for d,old in zip(derivatives[j],gradient)]
                value*=values[j]
        elif op=='inv':
            assert data<i and values[data]
            value=1/values[data];gradient=[-g*value*value for g in derivatives[data]]
        else:raise AssertionError(op)
        values.append(value);derivatives.append(gradient)
    return values,derivatives


def main(replay,output):
    directory=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63'
    d=read(directory/'report.json');f=read(replay/'report.json');assert d['all_passed']
    for name,digest in d['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in d['evidence_sha256'].items():assert sha(directory/name)==digest,name
    assert (directory/'map.json.gz').read_bytes()==(replay/'map.json.gz').read_bytes()
    assert (directory/'exact-sample.json').read_bytes()==(replay/'exact-sample.json').read_bytes()
    def replace(value):
        if isinstance(value,str):return value.replace(str(replay.relative_to(ROOT)),str(directory.relative_to(ROOT)))
        if isinstance(value,list):return [replace(v)for v in value]
        if isinstance(value,dict):return {k:replace(v)for k,v in value.items()}
        return value
    record=read(directory/'map-record.json');assert record==replace(read(replay/'map-record.json'))
    payload=read(directory/'map.json.gz');sample=read(directory/'exact-sample.json')
    assert payload['component_id']==record['component_id']=='R4_POSITIVE_WEIGHT6_25'
    names=payload['input_names'];nodes=payload['nodes'];outputs=payload['output_u1_through_u61'];meta=payload['normalization']
    assert len(names)==18 and len(nodes)==6079 and len(outputs)==61
    assert names==['n0','n1','n2','scale',*['base'+str(i)for i in range(1,5)],*['g'+str(i)for i in range(10)]]
    family_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz'
    family=read(family_path);n=s.symbols('n0:3');scale=s.Symbol('scale');base=s.symbols('base1:5')
    physical=bivector_from_coefficients(list(map(s.sympify,family['canonical150'])))
    translation=s.eye(5)
    for i,parameter in enumerate(base,1):translation[i,0]=1-parameter
    source=pushforward(scale*physical,translation)
    # The old universal normalizer is audited already. Reconstruct its ENTIRE
    # generic circuit from this exact family and compare all nodes and outputs.
    dag,expected,metadata=construct({'id':'R4_POSITIVE_WEIGHT6_25'},source,[*n,scale,*base])
    assert json.loads(json.dumps(dag.nodes))==nodes and expected==outputs and json.loads(json.dumps(metadata))==meta
    values,gradient=dual(nodes,sample['point'],names)
    U=[s.Rational(values[i])for i in outputs]
    J=s.Matrix([[s.Rational(v)for v in gradient[i]]for i in outputs])
    assert U==list(map(s.Rational,sample['output_u1_through_u61']))
    assert J==s.Matrix(sample['full61_by18_Jacobian'])
    assert len(DomainMatrix.from_Matrix(J).rref()[1])==16
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    for matrix,key,rank in [(J,'full_rank16_minor',16),(J[one,:],'primitive_rank15_minor',15)]:
        witness=sample[key];assert witness['rank']==rank
        determinant=DomainMatrix.from_Matrix(matrix[witness['rows'],witness['columns']]).convert_to(s.QQ).det()
        assert s.Rational(determinant)==s.Rational(witness['determinant'])!=0
        assert len(DomainMatrix.from_Matrix(matrix).rref()[1])==rank
    assert values[meta['value_pivot_node']]==Fraction(1,36)
    assert values[meta['value_pfaffian_node']]==Fraction(-16,81)
    # Independently compose the source-point translation with the evaluated
    # normalizer, then normalize the ORIGINAL physical family directly.
    point=sample['point'];M=s.Matrix([[s.Rational(values[j])for j in row]for row in meta['linear_normalizing_matrix_nodes']])
    sub={parameter:s.Integer(point[str(parameter)])for parameter in [*n,scale,*base]}
    full=M*translation.subs(sub)
    original=(scale*physical).subs(sub)
    direct=model.normalize(pushforward(original,full))
    assert direct['parameters']==U
    checks=model.verify_point(direct);assert checks==sample['all_original_S61_checks']
    assert coefficients(model.reconstruct(U))==direct['coefficients']
    z,uu,beta=primitive_beta(model)
    primitive=beta.subs(dict(zip(uu,[U[i]for i in one])))
    import itertools
    beta24=[primitive[i,j].coeff(variable)for i,j in itertools.combinations(range(4),2)for variable in z]
    assert list(map(str,beta24))==sample['primitive_beta24']
    inv=curve_invariants_from_beta(beta24)
    assert list(map(str,inv))==[sample[k]for k in ['canonical_curve_c4','canonical_curve_c6','canonical_curve_discriminant','canonical_curve_j']]
    assert inv[2]!=0 and inv[2]==(inv[0]**3-inv[1]**2)/1728
    firstjet=read(ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/report.json');assert firstjet['all_passed']
    assert firstjet['exact_primitive_image_dimension15_not_dominant'] and firstjet['new_component_generic_fiber_dimension1_over_its_image']
    component=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json')
    assert component['actual_global_projective_irreducible_component_dimension_exact']==25
    assert record['expected_s61_dimension']==25-9==16 and record['generic_firstjet_image_dimension']==15
    assert record['generic_fiber_over_image_dimension']==16-15==1
    assert record['full_source_dimension']==18 and record['full_s61_component_dominance']
    assert ROOT/record['exact_rank16_witness']==directory/'exact-sample.json'
    assert ROOT/record['map']['path']==directory/'map.json.gz' and record['map']['sha256']==sha(directory/'map.json.gz')
    decoded=gzip.decompress((directory/'map.json.gz').read_bytes())
    assert hashlib.sha256(decoded).hexdigest()==record['decoded_sha256']
    assert record['input_names']==names and record['coordinate_order']=='u1,...,u61, existing S61 model'
    assert not record['defining_prime_ideal_constructed'] and not record['firstjet_image_hypersurface_equation_constructed']
    historical=read(ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json')['records']
    assert sha(ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json')==component['input_sha256']['code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json']
    ids=[r['component_id']for r in historical]+[record['component_id']]
    assert len(ids)==len(set(ids))==52
    assert record['component_id']not in [r['component_id']for r in historical]
    paths=[Path(__file__),directory/'report.json',directory/'map.json.gz',directory/'map-record.json',directory/'exact-sample.json',
        family_path,ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json',
        ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json',
        ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',ROOT/'code/s61_model.py',
        ROOT/'code/elliptic_curve_invariants_20261101.py',
        ROOT/'research/notes/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.md']
    report=dict(id='Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW',all_passed=True,
        fresh_producer_map_and_sample_bytes_and_record_after_path_rebasing_match=True,
        entire6079node_generic_circuit_rebuilt_from_actual_family_and_audited_normalizer=True,
        independent_dense18input_Fraction_AD_all61values_and61by18matrix=True,
        exact_rank16_and_rank15_minors_checked=True,
        direct_ORIGINAL_physical_family_point_frame_all61_and134_and150_and_Jacobi_checked=True,
        all_inverse_denominators_nonzero_at_exact_regular_QQ_sample=True,
        image_dimension16_full_new_S61_component_dominance_accepted=True,
        generic_firstjet_dimension15_from_separate_exact_function_field_proof=True,
        generic_fiber_over_image_dimension=1,
        exact_curve_Delta_nonzero_proves_new_image_not_discriminant_hypersurface=True,
        supplementary_source_record_hash_order_input_and_component_identity_valid=True,
        old51_records_plus_new_record_form52_distinct_rank4_maps=True,
        historical51_catalogue_not_mutated=True,
        sample_extrapolation_defining_prime_hypersurface_equation_or_generic_five_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
    with output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print('PASS independent densefull61map16/beta15, genericconstruction+AD+physicalnormalization+52recordintegration')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();main(args.replay.resolve(),args.output)
