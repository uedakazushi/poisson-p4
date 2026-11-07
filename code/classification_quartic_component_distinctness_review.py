#!/usr/bin/env python3
"""Independent nil4 quartic distinctness audit against immutable baseline89.

Recompute full witness Jacobi/Pfaffian and the plane normal two-jet using only
the independent standard-library sparse polynomial arithmetic.
"""
import argparse
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path

from classification_quadratic_cofactor_review import Ring, PAIRS


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--baseline',type=Path,required=True)
    ap.add_argument('--candidate-audit',type=Path,required=True)
    ap.add_argument('--ambient-review',type=Path,required=True)
    ap.add_argument('--ambient-producer',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    baseline_raw = gzip.decompress(args.baseline.read_bytes())
    catalog = json.loads(baseline_raw)
    producer = json.loads(args.candidate_audit.read_bytes())
    ambient = json.loads(args.ambient_review.read_bytes())
    assert hashlib.sha256(baseline_raw).hexdigest()==producer['baseline_catalogue_original_sha256']
    rows = [r for r in catalog['components'] if r.get('counted_in_global_lower_bound')]
    assert len(rows)==89 and sum(r.get('count',1) for r in rows)==89
    assert len({r['id'] for r in rows})==89
    assert sum(r.get('count',1) for r in rows if r['generic_rank']==4)==83
    assert sum(r.get('count',1) for r in rows if r['generic_rank']==2)==6
    toric_ids = {f'T{i:02}' for i in range(1,41)}
    mps = [r for r in rows if r['id'] in toric_ids]
    assert {r['id'] for r in mps}==toric_ids
    assert {r['projective_dimension'] for r in mps}=={23,24,25}
    matching = [r for r in rows if r['generic_rank']==4 and r['projective_dimension']==22]
    assert len(matching)==14
    assert {r['id'] for r in matching}==set(producer['known_rank4_dimension22_ids'])
    nonreduced = [r for r in matching if not r['generic_actual_degeneracy_divisor']['reduced']]
    reduced = [r for r in matching if r['generic_actual_degeneracy_divisor']['reduced']]
    assert len(nonreduced)==13
    assert [r['id'] for r in reduced]==['R4_TWO_HYPERPLANES_CUBIC_22']
    for r in nonreduced:
        divisor = r['generic_actual_degeneracy_divisor']
        degrees = divisor['irreducible_factor_degrees']
        powers = divisor['corresponding_multiplicities']
        assert len(degrees)==len(powers) and sum(d*m for d,m in zip(degrees,powers))==5
        assert max(powers)>1
    divisor = reduced[0]['generic_actual_degeneracy_divisor']
    assert divisor['irreducible_factor_degrees']==[1,1,3]
    assert divisor['corresponding_multiplicities']==[1,1,1]

    # First variable n0 is h; all four differentiation coordinates remain
    # x,y,z,w. Since h is a Casimir, the other six raw Jacobi expressions vanish.
    R = Ring(1)
    h,x,y,z,w = [R.variable(n) for n in R.names]
    expressions = ['36*n0*w+6*x*z-8*y**2','6*w*x-10*y*z','-10*w*y-z**2',
                   '6*w*y-9*z**2','-12*w*z','-12*w**2']
    p = [R.parse(s) for s in expressions]
    ambient_producer = json.loads(args.ambient_producer.read_bytes())['nil_regular_component']
    assert p==[R.parse(ambient_producer['witness_brackets'][f'{i+1},{j+1}'].replace('h','n0')) for i,j in PAIRS]
    assert all(sum(ex)==2 for f in p for ex in f)
    P = [[{} for j in range(4)] for i in range(4)]
    for (i,j),v in zip(PAIRS,p):
        P[i][j],P[j][i] = v,R.scale(v,-1)
    jac = [R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t))
                   for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4)))
           for i,j,k in it.combinations(range(4),3)]
    assert all(not j for j in jac)
    pf = R.add(R.mul(p[0],p[5]),R.scale(R.mul(p[1],p[4]),-1),R.mul(p[2],p[3]))
    F = R.scale(R.add(R.mul(R.add(R.scale(R.mul(w,y),2),R.scale(R.mul(z,z),-1)),
                           R.add(R.scale(R.mul(w,y),2),R.scale(R.mul(z,z),-1))),
                     R.scale(R.mul(h,w,w,w),-48)),9)
    assert pf==F
    assert {ex[0] for ex in F}=={0,1}
    linear_h = {ex[1:]:cv for ex,cv in F.items() if ex[0]==1}
    assert linear_h=={(0,0,0,3):-432}
    constant_h_at_w_zero = {ex[1:]:cv for ex,cv in F.items() if ex[0]==0 and ex[4]==0}
    assert constant_h_at_w_zero=={(0,0,4,0):9}
    D = R.mul(h,F)
    assert all(ex[3]+ex[4]>=2 for ex in D)
    normal2 = {ex:cv for ex,cv in D.items() if ex[3]+ex[4]==2}
    assert normal2==R.scale(R.mul(h,y,y,w,w),36)
    # Ambient proof is an input, not recomputed here. Its witness must coincide
    # with these six brackets; its verified local support bound is necessary.
    assert ambient['all_passed'] and ambient['orbit_closure_is_full_component']
    assert ambient['projective_orbit_dimension']==ambient['projective_ambient_local_dimension']==22
    assert ambient['derivative_rank']==102 and ambient['actual_orbit_vertical_scale_rank']==47
    assert ambient['cubic_augmented_rank']==103
    assert ambient['left_annihilates_full_D_and_every_kernel_mixed_column']
    dispositions = []
    for r in rows:
        if r['generic_rank']==2:
            reason = 'generic_rank_2_versus_4'
        elif r['projective_dimension']!=22:
            reason = 'different_proven_projective_dimension'
        elif r['id']==reduced[0]['id']:
            reason = 'closed_1_plus_1_plus_3_factorization_image'
        else:
            assert r in nonreduced
            reason = 'closed_nonreduced_quintic_cone'
        dispositions.append(dict(id=r['id'],projective_dimension=r['projective_dimension'],
                                 generic_rank=r['generic_rank'],reason=reason))
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,
                  scope='Independent distinctness audit conditional on separately verified full ambient component proof; that input is accepted, not rerun.',
                  baseline_count=89,baseline_rank4_count=83,baseline_rank2_count=6,
                  uncompressed_baseline_sha256=hashlib.sha256(baseline_raw).hexdigest(),
                  full_baseline_dispositions=dispositions,
                  MPS40_dimensions=[23,24,25],same_rank_dimension14_verified=True,
                  nonreduced22_catalogue_facts=[dict(id=r['id'],divisor=r['generic_actual_degeneracy_divisor'],
                                                  proof_sources=r.get('proof_sources',[])) for r in nonreduced],
                  unique_reduced22_id=reduced[0]['id'],
                  all_six_witness_brackets=expressions,all_ten_raw_Jacobi_zero=True,
                  full_Pfaffian_identity_verified=True,
                  irreducibility_argument='Primitive degree-one polynomial in h over C[x,y,z,w]: h coefficient=-432w^3, constant coefficient at w=0 is9z^4, hence content1. Gauss proves absolute irreducibility.',
                  actual_D_reduced_factor_degrees=[1,4],new_generic_rank=4,
                  plane_normal_twojet=R.encode(normal2),plane_normal_generic_rank=1,
                  plane_incidence_closed_argument='Proper Grassmannian incidence: D in I_plane^2 and identically zero binary quadratic discriminant. Include D=0 in the affine cone.',
                  SNC_exclusion_argument='Any plane in Sing(D) for an SNC divisor lies generically in a transverse double intersection; triple strata have dimension at most1. Its normal quadratic therefore has rank2, incompatible with incidence.',
                  closed_factorization_argument='The nonreduced quintic cone is union of cones over proper images (G,R)->G^2 R, degG=1 or2. The 1+1+3 image is proper. Include zero; all old component closures stay inside their conditions.',
                  full_ambient_component_proof_accepted_as_independent_input=True,
                  new_component_distinct_from_every_baseline_component=True,
                  permitted_new_lower_bound=90,permitted_rank4_lower_bound=84,
                  global_exhaustion_claimed=False,
                  input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                [args.baseline,args.candidate_audit,args.ambient_review,args.ambient_producer,Path('code/classification_quadratic_cofactor_review.py')]},
                  reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  reviewer_source_bytes=Path(__file__).read_text())
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_passed=True,new_distinct_component_dimension=22,permitted_lower_bound=90,permitted_rank4_lower_bound=84)))


if __name__=='__main__':
    main()
