#!/usr/bin/env python3
"""Separate-author fullambient check of the two-conic common-bitangent source.

Uses earlier independent Fraction/S61 coefficient-loop arithmetic, not the
two-conic producer or its author-written component checker. No CAS imports.
"""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_quadratic_cofactor_review import Ring


def actual_source(ambient,source):
    R = Ring(7) # n0=h, n1:7=all six magnetic coefficients.
    def parse(text):
        return R.parse(re.sub(r'\b(h|m[0-5])\b',lambda m:'n0' if m[1]=='h' else 'n'+str(int(m[1][1:])+1),text))
    p = {ij:{} for ij in q.PAIRS}
    for key,text in ambient['full_six_magnetic_source'].items():
        p[ast.literal_eval(key)] = parse(text)
    scale = F(ambient['integer_scale'])
    for ij,text in zip(it.combinations(range(1,5),2),source['brackets']):
        assert p[ij]==R.scale(parse(text),scale)
    full = dict(p)
    full.update({(j,i):R.scale(v,-1) for (i,j),v in p.items()})
    for i,j,k in it.combinations(range(1,5),3):
        assert not R.add(*(R.mul(full.get((a,t),{}),R.derivative(full.get((b,c),{}),t-1))
                           for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for t in range(1,5)))
    pf = R.add(R.mul(p[1,2],p[3,4]),R.scale(R.mul(p[1,3],p[2,4]),-1),R.mul(p[1,4],p[2,3]))
    conics = [parse(source[k]) for k in ('q1','q2')]
    ratio = F(source['Pfaffian_multiple_of_q1q2'])*scale**2
    assert pf==R.scale(R.mul(*conics),ratio)
    D = R.mul(R.variable('n0'),pf)
    assert all(not any(ex[1:7]) and not any(ex[-2:]) for ex in D)
    # The two smooth conics are h*x-(x+/-4*y)^2, up to the unit1/16.
    x,y,h = R.variable('x'),R.variable('y'),R.variable('n0')
    if source['name']=='contact31_common_bitangent':
        expected = [R.scale(R.add(R.mul(h,x),R.scale(R.mul(R.add(x,R.scale(y,v)),R.add(x,R.scale(y,v))),-1)),F(1,16)) for v in (4,-4)]
    elif source['name']=='contact4_other_tangent':
        q1 = R.add(R.mul(h,x),R.scale(R.mul(y,y),-1))
        expected = [q1,R.add(q1,R.mul(x,x))]
    else:
        raise ValueError(source['name'])
    assert conics==expected
    # R only differentiates x,y,z,w; differentiate magnetic parameters directly.
    def parameter_derivative(poly,index):
        out = {}
        for ex,c in poly.items():
            if ex[index]:
                e = list(ex)
                e[index] -= 1
                out[tuple(e)] = c*ex[index]
        return out
    values = [F(ambient['parameter_values'][f'm{i}']) for i in range(6)]
    def specialize(poly):
        out = {}
        for ex,c in poly.items():
            for power,value in zip(ex[1:7],values):
                c *= value**power
            exponent = (ex[0],)+ex[-4:]+(0,0,0)
            out[exponent] = out.get(exponent,F(0))+c
        return {e:c for e,c in out.items() if c}
    point = q.antisymmetric({ij:specialize(v) for ij,v in p.items()})
    directions = [q.antisymmetric({ij:specialize(parameter_derivative(v,i+1)) for ij,v in p.items()}) for i in range(6)]
    def parse_q(text):
        return q.parse(re.sub(r'\b([xyzw])\b',lambda m:dict(x='x1',y='x2',z='x3',w='x4')[m[1]],text))
    for ij in q.PAIRS:
        text = ambient['point_brackets'].get(str(ij),'0')
        assert point.get(ij,{})==parse_q(text)
        for i in range(6):
            assert directions[i].get(ij,{})==parse_q(ambient['all6_parameter_directions'][i].get(str(ij),'0'))
    return point,directions,dict(all_six_magnetic_parameters_kept=True,full_source_raw_Jacobi_zero=True,
                                 actual_Pfaffian_constant_times_q1q2_for_every_parameter=True,
                                 actual_D_essential_at_most3_on_entire_source=True,
                                 q1_q2_explicit_smooth_distinct=True,source_irreducible='A6 magnetic parameters with connected PGL5 saturation')


def ambient_check(P,directions,source,normal):
    assert len(normal['normals'])==1
    R = q.antisymmetric({tuple(map(int,k.split(','))):q.parse(v) for k,v in normal['normals'][0].items()})
    Dcols = []
    for ij,e in q.COLS:
        variation = dict(P)
        variation[ij] = q.add(P.get(ij,{}),{e+(0,0,0):F(1)})
        variation[ij[::-1]] = q.scale(variation[ij],-1)
        Dcols.append(q.jvec(q.jacobi(variation)))
    D = q.transpose(Dcols)
    Vcols,Ocols = [],[]
    for a,b in it.product(range(5),repeat=2):
        V,O = {},{}
        for i,j in q.PAIRS:
            V[i,j] = q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
            O[i,j] = q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},
                           q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
        Vcols.append(q.vec(V));Ocols.append(q.vec(O))
    Mcols = Vcols+Ocols+[q.vec(P)]
    Fcols = Mcols+[q.vec(v) for v in directions]
    Kcols = Fcols+[q.vec(R)]
    assert q.rank(q.transpose(Vcols))[0]==24
    assert q.rank(q.transpose(Mcols))[0]==47
    rankF,basisF = q.rank(q.transpose(Fcols))
    assert rankF==48 and q.rank(q.transpose(Kcols))[0]==49
    assert all(sum(v*c for v,c in zip(row,col))==0 for row in D for col in Kcols)
    certificate = source['certificates']['D']
    rows,cols = certificate['rows'],certificate['columns']
    assert len(rows)==len(cols)==101
    determinant = q.determinant([[D[i][j] for j in cols] for i in rows])
    assert determinant==F(certificate['determinant']) and determinant
    # Rank>=101 from the minor; the49 independent kernel columns give rank<=101.
    left = {int(i):F(v) for i,v in normal['one_normal_quadratic_cokernel_functional']['nonzero_coefficients'].items()}
    assert all(sum(v*D[i][j] for i,v in left.items())==0 for j in range(150))
    Q = q.jvec(q.jacobi(R))
    pairing = sum(v*Q[i] for i,v in left.items())
    assert pairing==F(normal['one_normal_quadratic_cokernel_functional']['pairing'])==F(-120295458720)
    # A literal complement uses all101 verified standard-coordinate D columns
    # plus R; together with48 actual-source basis columns they span C150.
    return dict(full_D_shape=[350,150],full_D_exact_rank=101,vertical_rank=24,
                actual_group_vertical_scale_rank=47,all_six_source_derivative_directions_reconstructed=True,
                actual_source_exact_rank=48,entire_kernel_dimension=49,
                complete_kernel_equals_actual_source_plus_R=True,all_D_times_full_kernel_zero=True,
                D_minor=dict(rows=rows,columns=cols,determinant=str(determinant)),
                actual_source_basis_columns=basisF,true_slice_standard_coordinate_columns=cols,
                true_slice_dimension=102,true_slice_cut_count=48,
                left_annihilates_all150_original_columns=True,
                nonzero_quadratic_pairing=str(pairing),normal_quadratic_support_only_origin=True,
                local_raw_dimension_upper=48,actual_raw_family_dimension_lower=48,
                full_projective_component_dimension=23,
                componenthood=True,local_uniqueness_not_claimed=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    start = time.monotonic()
    paths = ['reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json',
             'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json',
             'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_NORMALS_20261107.json',
             'reports/CLASSIFICATION_BASELINE90_COMPONENT_PROJECTION_20261107.json']
    sources,ambient,normals,baseline = [json.loads(Path(p).read_text()) for p in paths]
    name = 'contact31_common_bitangent'
    source = next(r for r in sources['cases'] if r['name']==name)
    actual = next(r for r in ambient['cases'] if r['name']==name)
    normal = next(r for r in normals['cases'] if r['case']==name)
    point,dirs,whole = actual_source(actual,source)
    support = ambient_check(point,dirs,actual,normal)
    assert baseline['count']==len(baseline['rows'])==90
    assert baseline['rank4_count']==84 and baseline['rank2_count']==6
    dispositions = []
    for r in baseline['rows']:
        divisor = r.get('generic_actual_degeneracy_divisor') or {}
        if r['generic_rank']==2:
            reason = 'different_generic_rank'
        elif r['projective_dimension']!=23:
            reason = 'different_component_dimension'
        elif r['id'] in {f'T{i:02}' for i in range(1,41)}:
            reason = 'entire_new_component_essential3_closed_condition_excludes_SNC'
        else:
            degrees = divisor.get('irreducible_factor_degrees')
            multiplicities = divisor.get('corresponding_multiplicities')
            assert (degrees,multiplicities)!=([1,2,2],[1,1,1])
            assert degrees is not None and multiplicities is not None
            reason = 'different_generic_irreducible_factor_degree_multiplicity_stratum'
        dispositions.append(dict(id=r['id'],projective_dimension=r['projective_dimension'],generic_rank=r['generic_rank'],reason=reason))
    cplus = next(r for r in baseline['rows'] if r['id']=='R4_CPLUS_EXT')
    assert cplus['projective_dimension']==24
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
                  whole_source=whole,fullambient_support=support,
                  proof='On S=U101+R complementary to actual F48, selected101 equations eliminate U by IFT. The full left cokernel pairing gives nonzero quadratic coefficient in the remaining variable, so S has zero-dimensional support. Krull dimension inequality bounds full rawlocal dimension by48, matching the actual irreducible source image. Its saturated closure is therefore a full component of projective dimension23.',
                  full_baseline90_distinctness=dispositions,
                  CPLUS24_excluded_by_dimension=True,
                  MPS_exclusion='Essential<=3 is a gradient-rank closed condition on the full degeneracy section, including zero. A nonzero quintic in three linear variables vanishes to order5 along a nonempty projective vertex; an SNC divisor on a fourfold has order<=4 at every point. Thus the entire new component has no SNC point.',
                  other_same_dimension_distinctness='Generic exact factor-degree/multiplicity strata are constructible and disjoint. If two irreducible component closures were equal, their dense constructible generic source images would intersect, contradicting the different factor stratum.',
                  full_source_componenthood_and_distinctness_accepted=True,permitted_lower_bound=91,permitted_rank4_lower_bound=85,
                  local_uniqueness_claimed=False,global_exhaustion_claimed=False,
                  input_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths+
                                [str(Path(q.__file__).relative_to(Path.cwd())),'code/classification_quadratic_cofactor_review.py']},
                  reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),
                  seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','maximum_RSS_KiB','permitted_lower_bound']}))


if __name__=='__main__':
    main()
