#!/usr/bin/env python3
"""True degree-three initial equations of the full resonant finite model.

Carry cubic parts through exact constant elimination of the quadratic parts.
The resulting cubic initials, together with all linear multiples of the 57
quadrics, span the degree-three initial-ideal component.  No full GB is used.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

from q_ordinary_quadratic_cone_20261103_v63 import TruncatedQQ
from q_weighted_finite_model_validation_20261103_v63 import evaluate,physical_equations

ROOT=Path(__file__).resolve().parents[1]
P=32003


class Truncated3(TruncatedQQ):
    inverse_constants=[]

    @staticmethod
    def mul(left,right):
        if not left or not right:return {}
        a=[[],[],[],[]];b=[[],[],[],[]]
        for mon,c in left.items():a[len(mon)].append((mon,c))
        for mon,c in right.items():b[len(mon)].append((mon,c))
        out={}
        for i in range(4):
            for j in range(4-i):
                for m,c in a[i]:
                    for n,d in b[j]:
                        mon=tuple(sorted(m+n));value=out.get(mon,Fraction(0))+c*d
                        if value:out[mon]=value
                        else:out.pop(mon,None)
        return out

    @classmethod
    def inv(cls,value):
        c=value.get((),Fraction(0));assert c and c.numerator%P and c.denominator%P
        cls.inverse_constants.append(str(c))
        positive={mon:-coefficient/c for mon,coefficient in value.items() if mon}
        out=cls.one;power=cls.one
        for _ in range(3):
            power=cls.mul(power,positive);out=cls.add(out,power)
        out={mon:coefficient/c for mon,coefficient in out.items()}
        assert cls.mul(value,out)==cls.one
        return out


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(poly):return [[list(mon),str(c)] for mon,c in sorted(poly.items())]


def subtract_scaled(row,other,scale):
    for mon,c in other.items():
        value=row.get(mon,Fraction(0))-scale*c
        if value:row[mon]=value
        else:row.pop(mon,None)


def rank_QQ(rows):
    basis={};pivots=[]
    for index,original in sorted(enumerate(rows),key=lambda pair:(len(pair[1]),pair[0])):
        row=original.copy()
        while row:
            pivot=min(row);coefficient=row[pivot]
            if pivot not in basis:
                basis[pivot]={mon:c/coefficient for mon,c in row.items()};pivots.append(index);break
            subtract_scaled(row,basis[pivot],coefficient)
    return basis,pivots


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz')
    parser.add_argument('--quadratic-cone',type=Path,default=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(3000000000,3000000000))
    model=json.load(gzip.open(args.model,'rt'))
    old=json.load(gzip.open(args.quadratic_cone/'quadrics.json.gz','rt'))
    names=old['variables'];weights=model['input_weights'];assert len(names)==35
    inputs=[{(i,):Fraction(1)} for i in range(35)];inputs[0][()]=Fraction(2,9)
    values,inverses=evaluate(model['nodes'],inputs,Truncated3)
    reconstructed=[values[node] for node in model['raw150_reconstruction_nodes']]
    residuals=[values[node] for node in model['all350_residual_nodes']]
    assert all(len(mon) in (2,3) for row in residuals for mon in row)
    _,physical=physical_equations(reconstructed,Truncated3);assert residuals==physical
    quadrics=[{mon:c for mon,c in row.items() if len(mon)==2} for row in residuals]
    cubics=[{mon:c for mon,c in row.items() if len(mon)==3} for row in residuals]
    for new,previous in zip(quadrics,old['all350_quadratic_rows']):
        assert new=={(i,j):Fraction(c) for i,j,c in previous}
    print('Exact full350 degree3 DAG and physical replay complete',time.monotonic()-start,flush=True)
    basis={};carries={};basis_sources=[];operations=[];relations=[]
    for i in sorted(range(350),key=lambda i:(len(quadrics[i]),i)):
        q=quadrics[i].copy();c=cubics[i].copy();reductions=[]
        while q:
            pivot=min(q);coefficient=q[pivot]
            if pivot not in basis:
                assert coefficient.numerator%P and coefficient.denominator%P
                basis[pivot]={mon:value/coefficient for mon,value in q.items()}
                carries[pivot]={mon:value/coefficient for mon,value in c.items()}
                basis_sources.append(i)
                operations.append(dict(source=i,new_pivot=list(pivot),normalizer=str(coefficient),reductions=reductions))
                break
            reductions.append([list(pivot),str(coefficient)])
            subtract_scaled(q,basis[pivot],coefficient)
            subtract_scaled(c,carries[pivot],coefficient)
        else:
            operations.append(dict(source=i,cubic_initial=True,reductions=reductions))
            relations.append((i,c))
    assert len(basis)==57 and basis_sources==old['selected_original_rows']
    assert all(c.denominator%P for row in [*quadrics,*cubics,*[poly for i,poly in relations]] for c in row.values())
    for i,poly in relations:
        assert all(sum(weights[j] for j in mon)==model['residual_weight_bounds'][i] for mon in poly)
    active_names=['u','z3','z5','z8','z13','z14','z15','z16','z17','z24','z25','z26','z28']
    active={names.index(name) for name in active_names}
    plane=[(i,{mon:c for mon,c in poly.items() if all(j in active for j in mon)}) for i,poly in relations]
    assert all(not any(all(j in active for j in mon) for mon in q) for q in quadrics)
    plane_basis,plane_independent=rank_QQ([poly for i,poly in plane])
    # Build the complete cubic initial-space matrix, not just surviving G_i.
    byweight=defaultdict(lambda:dict(quadratic_multiples=[],new_initials=[]))
    for source in basis_sources:
        q=quadrics[source];weight=model['residual_weight_bounds'][source]
        for variable in range(35):
            product={tuple(sorted((*mon,variable))):c for mon,c in q.items()}
            byweight[weight+weights[variable]]['quadratic_multiples'].append((source,variable,product))
    for source,poly in relations:
        if poly:byweight[model['residual_weight_bounds'][source]]['new_initials'].append((source,poly))
    assert sum(len(block['quadratic_multiples']) for block in byweight.values())==35*57
    matrix=args.output/'cubic-blocks.in'
    with matrix.open('w') as stream:
        stream.write(str(len(byweight))+'\n')
        for weight,block in sorted(byweight.items()):
            rows=[poly for source,variable,poly in block['quadratic_multiples']]+[poly for source,poly in block['new_initials']]
            columns=sorted({mon for poly in rows for mon in poly});index={mon:i for i,mon in enumerate(columns)}
            stream.write(f"{weight} {len(columns)} {len(block['quadratic_multiples'])} {len(block['new_initials'])}\n")
            for poly in rows:
                data=[(index[mon],c.numerator*pow(c.denominator,-1,P)%P) for mon,c in sorted(poly.items())]
                stream.write(str(len(data))+' '+' '.join(f'{j} {c}' for j,c in data)+'\n')
    artifact=dict(variables=names,old_weights=weights,
                  raw350_cubic_terms=[encode(row) for row in cubics],
                  exact_constant_row_elimination=operations,
                  cubic_initial_relations=[dict(source=i,terms=encode(poly)) for i,poly in relations],
                  plane_active_variables=active_names,
                  plane_restricted_relations=[dict(source=i,terms=encode(poly)) for i,poly in plane if poly],
                  plane_normalized_QQ_basis=[encode(poly) for pivot,poly in sorted(plane_basis.items())],
                  quadratic_linear_multiples_recipe=dict(quadratic_original_rows=basis_sources,variables=list(range(35))))
    encoded=(json.dumps(artifact,separators=(',',':'))+'\n').encode();archive=args.output/'cubic-initials.json.gz'
    with archive.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as z:z.write(encoded)
    blockinfo={weight:dict(quadratic_multiples=len(block['quadratic_multiples']),
                           new_cubic_rows=len(block['new_initials']),
                           distinct_cubic_monomials=len({mon for row in block['quadratic_multiples'] for mon in row[2]}|
                                                       {mon for row in block['new_initials'] for mon in row[1]}))
               for weight,block in sorted(byweight.items())}
    report=dict(id='Q_ORDINARY_CUBIC_INITIAL_20261103_V63',all_passed=True,
                coefficient_ring='QQ',ordinary_local_variables=35,ordinary_truncation_degree=3,
                all350_constant_linear_terms_zero=True,all350_quadratic_parts_match_frozen_packet=True,
                all350_degree3_jets_match_independent_physical_Jacobi=True,
                inverse_Taylor_products_verified_to_degree3=True,
                quadratic_pivot_rank=57,all_quadratic_elimination_pivots_32003_units=True,
                cubic_constant_syzygy_rows=len(relations),nonzero_true_cubic_initial_rows=sum(bool(poly) for i,poly in relations),
                all1995_linear_quadratic_multiples_included=True,
                complete_cubic_initial_space='span_QQ{constant-quadratic-eliminated cubic terms, all35 variable multiples of57 quadrics}',
                all_saved_coefficients_32003_integral=True,
                false13plane_active_variables=active_names,
                nonzero_cubic_restrictions_to_false13plane=sum(bool(poly) for i,poly in plane),
                restricted_cubic_span_rank_QQ=len(plane_basis),
                false13plane_cut_by_true_cubic_initials=bool(plane_basis),
                block_shapes=blockinfo,
                matrix_rank_mod_p_not_yet_computed=True,
                full_local_dimension_or_component_decomposition_claimed=False,
                archive_sha256=sha(archive),matrix_sha256=sha(matrix),
                peak_rss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),args.model,args.quadratic_cone/'quadrics.json.gz',
                               ROOT/'code/q_ordinary_quadratic_cone_20261103_v63.py',
                               ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py']},
                seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','block_shapes']},indent=2),flush=True)


if __name__=='__main__':main()
