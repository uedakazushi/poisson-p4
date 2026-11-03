#!/usr/bin/env python3
"""Explicit dominant61-coordinate circuit for the52nd known rank4 component."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

from s61_model import S61,X,pushforward
from mps_s61_arithmetic_circuits_20261103_v62 import construct
from q_positive_weight6_small_model_review_20261103_v63 import read,rank
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def differentiate(dag,point):
    values=[];gradients=[]
    for node in dag.nodes:
        op=node[0];gradient={}
        if op=='q':value=Fraction(node[1])
        elif op=='v':value=Fraction(point[dag.names[node[1]]]);gradient={node[1]:Fraction(1)}
        elif op=='+':
            value=sum((values[j]for j in node[1]),Fraction(0))
            for j in node[1]:
                for k,c in gradients[j].items():gradient[k]=gradient.get(k,Fraction(0))+c
        elif op=='*':
            prefix=[Fraction(1)]
            for j in node[1]:prefix.append(prefix[-1]*values[j])
            value=prefix[-1];suffix=Fraction(1)
            for i in reversed(range(len(node[1]))):
                factor=prefix[i]*suffix
                if factor:
                    for k,c in gradients[node[1][i]].items():gradient[k]=gradient.get(k,Fraction(0))+factor*c
                suffix*=values[node[1][i]]
        elif op=='inv':value=1/values[node[1]];gradient={k:-value*value*c for k,c in gradients[node[1]].items()}
        else:raise ValueError(op)
        values.append(value);gradients.append({k:c for k,c in gradient.items()if c})
    return values,gradients


def minor(M):
    columns=list(DomainMatrix.from_Matrix(M).rref()[1]);rows=list(DomainMatrix.from_Matrix(M[:,columns].T).rref()[1])
    det=DomainMatrix.from_Matrix(M[rows,columns]).convert_to(s.QQ).det();assert det
    return dict(rank=len(columns),rows=rows,columns=columns,determinant=str(det))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    familypath=ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz';family=read(familypath)
    n=s.symbols('n0:3');scale=s.Symbol('scale');base=s.symbols('base1:5');parameters=[*n,scale,*base]
    physical=s.Matrix([[s.sympify(v).subs(dict(zip(s.symbols('h x y z w'),X)),simultaneous=True)for v in row]for row in family['full_canonical_homogeneous_P']])
    translation=s.eye(5)
    for i,p in enumerate(base,1):translation[i,0]=1-p
    source=pushforward(scale*physical,translation)
    dag,outputs,metadata=construct({'id':'R4_POSITIVE_WEIGHT6_25'},source,parameters)
    point=dict(n0=1,n1=1,n2=1,scale=1,base1=0,base2=0,base3=0,base4=0)
    point.update({f'g{i}':v for i,v in enumerate([1,0,0,1,0,0,0,0,0,0])})
    values,gradients=differentiate(dag,point)
    U=[s.Rational(values[j])for j in outputs]
    Jacobian=s.Matrix([[s.Rational(gradients[j].get(k,0))for k in range(len(dag.names))]for j in outputs])
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    fullminor=minor(Jacobian);betaminor=minor(Jacobian[one,:]);assert fullminor['rank']==16 and betaminor['rank']==15
    matrix=s.Matrix([[s.Rational(values[j])for j in row]for row in metadata['linear_normalizing_matrix_nodes']])
    source_sample=source.subs({p:point[str(p)]for p in parameters})
    normalized=model.normalize(pushforward(source_sample,matrix))
    assert U==normalized['parameters'];checks=model.verify_point(normalized);assert all(checks.values())
    z,u,beta=primitive_beta(model);b=beta.subs(dict(zip(u,[U[i]for i in one])))
    b24=[b[i,j].coeff(zz)for i,j in __import__('itertools').combinations(range(4),2)for zz in z]
    c4,c6,disc,j=curve_invariants_from_beta(b24);assert disc!=0
    payload=dict(component_id='R4_POSITIVE_WEIGHT6_25',input_names=dag.names,nodes=dag.nodes,
        output_u1_through_u61=outputs,normalization=metadata,
        coefficient_field='QQ',field_relation=None,residual_sp4_included=True,overall_scalar_retained=True,
        full_source_dimension=18,expected_s61_dimension=16,firstjet_image_dimension=15,
        original_family_translation='(translation(1-base))_*(scale*canonicalP(n)); fixed-one normalizer therefore evaluates originalP at(1:base).',
        generic_identity_proof='Exact full source Poisson; finite rational projective transformations and universal S61 normalization theorem. Sample verifies implementation, not generic identity by extrapolation.',
        dominance_proof='Image lies in actual16D S61 component; exact61x18 differential rank16 at the recorded regular QQ point.',
        component_identification_dependency='Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW: actualPGL-saturatedprojectivecomponent25, distinctprior57.')
    encoded=(json.dumps(payload,separators=(',',':'))+'\n').encode();map_path=args.output/'map.json.gz'
    map_path.write_bytes(gzip.compress(encoded,mtime=0))
    (args.output/'exact-sample.json').write_text(json.dumps(dict(point=point,
        value_pivot=str(values[metadata['value_pivot_node']]),value_pfaffian=str(values[metadata['value_pfaffian_node']]),
        output_u1_through_u61=list(map(str,U)),full61_by18_Jacobian=[[str(v)for v in row]for row in Jacobian.tolist()],
        full_rank16_minor=fullminor,primitive_rank15_minor=betaminor,all_original_S61_checks=checks,
        primitive_beta24=list(map(str,b24)),canonical_curve_c4=str(c4),canonical_curve_c6=str(c6),
        canonical_curve_discriminant=str(disc),canonical_curve_j=str(j),
        beta_image_is_not_canonical_curve_discriminant_hypersurface=True),indent=2)+'\n')
    record=dict(component_id='R4_POSITIVE_WEIGHT6_25',map=dict(path=str(map_path.relative_to(ROOT)),sha256=sha(map_path)),
        decoded_sha256=hashlib.sha256(encoded).hexdigest(),coordinate_order='u1,...,u61, existing S61 model',
        input_names=dag.names,map_kind='rational_family_point_scalar_Sp4_arithmetic_circuit',coefficient_field='QQ',
        field_relation=None,expected_s61_dimension=16,full_source_dimension=18,full_s61_component_dominance=True,
        exact_rank16_witness=str((args.output/'exact-sample.json').relative_to(ROOT)),
        generic_firstjet_image_dimension=15,generic_fiber_over_image_dimension=1,
        defining_prime_ideal_constructed=False,firstjet_image_hypersurface_equation_constructed=False)
    (args.output/'map-record.json').write_text(json.dumps(record,indent=2)+'\n')
    dependencies=[ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/report.json']
    for dependency in dependencies:
        data=read(dependency);assert data['all_passed']
        for p,digest in data['input_sha256'].items():assert sha(ROOT/p)==digest
        for p,digest in data.get('evidence_sha256',{}).items():assert sha(dependency.parent/p)==digest
    inputs=[Path(__file__),familypath,*dependencies,
        ROOT/'code/mps_s61_arithmetic_circuits_20261103_v62.py',ROOT/'code/s61_model.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/elliptic_curve_invariants_20261101.py',
        ROOT/'research/notes/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63.md']
    report=dict(id='Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63',all_passed=True,
        new52nd_rank4_component_explicit_full61_map=True,source_inputs=18,circuit_nodes=len(dag.nodes),
        exact_QQ_sample_rank61map=16,exact_QQ_sample_rank_firstjet=15,
        full_dense_S61_image_dimension16_proved=True,
        all61_normalized_coordinates_and_original_equations_verified_at_regular_sample=True,
        generic_identity_uses_exact_source_and_audited_coordinate_theorem_not_sampling=True,
        generic_curve_discriminant_nonzero_on_new_beta_image=True,
        new_beta_hypersurface_not_degree24_curve_discriminant=True,
        new_beta_hypersurface_equation_or_recognition_not_yet_proved=True,
        generic_five_point_exhaustion_or_component_defining_prime_not_claimed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('DENSE NEW52ND RANK4 S61 MAP RANK16/18;BETA15;DISC NONZERO PASS',report['seconds'],flush=True)


if __name__=='__main__':main()
