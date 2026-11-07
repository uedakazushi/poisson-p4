#!/usr/bin/env python3
"""Independent actual Ga9 curve/component/absolute-irred/prior98 review."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_minimal_quartic_review import LaurentRing

ROOT=Path(__file__).resolve().parents[1]

def parse(text):
    alias={'h':'h','x':'x1','y':'x2','z':'x3','w':'x4','epsilon':'tau0'}
    return q.parse(re.sub(r'\b(h|x|y|z|w|epsilon)\b',lambda m:alias[m[1]],text))

def modelparse(text):
    alias=['h','x1','x2','x3','x4']
    return q.parse(re.sub(r'\bu([0-4])\b',lambda m:alias[int(m[1])],text))

def decode(data):
    return q.antisymmetric({tuple(map(int,k.split(','))):parse(v)for k,v in data.items()})

def at_parameter(p,value):
    out={}
    for ex,c in p.items():
        co=c*value**ex[5]
        if co:
            e=list(ex);e[5]=0;e=tuple(e);out[e]=out.get(e,F())+co
    return {e:c for e,c in out.items()if c}

def quintic(P):
    out={}
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),
                 q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),
                 q.mul(P.get((i,l),{}),P.get((j,k),{})))
        out=q.add(out,q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
    return out

def actions(P):
    V=[];O=[]
    for a,b in it.product(range(5),repeat=2):
        vertical={};orbit={}
        for i,j in q.PAIRS:
            vertical[i,j]=q.add(q.atom(1,i,b)if j==a else{},q.atom(-1,j,b)if i==a else{})
            orbit[i,j]=q.add(P.get((b,j),{})if i==a else{},P.get((i,b),{})if j==a else{},
                             q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
        V.append(q.vec(vertical));O.append(q.vec(orbit))
    return V,V+O+[q.vec(P)]

def linear_substitution(M):
    return [q.add(*(q.atom(v,j)for j,v in enumerate(row)))for row in M]+[q.atom(1,i)for i in range(5,8)]

def mps_factors(packet,prior,inputs):
    records=[]
    wanted={c['id']for c in prior['components']if c['generic_rank']==4 and c['projective_dimension']==23 and c['meets_SNC_log_symplectic_locus']}
    assert len(wanted)==29 and {c['id']for c in packet['cases']}==wanted
    for case in packet['cases']:
        path=ROOT/case['input_path'];inputs.append(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest()==case['input_sha256']
        source=json.loads(path.read_text());params=case['all_general_parameters']
        R=LaurentRing(params,[f'x{i}'for i in range(5)],invert=[p for p in params if p.startswith('p')])
        # Invertible period reparametrizations make the only nonmonomial
        # source denominators single Laurent variables, without a sample.
        def read(text):
            if case['id']=='T18':text=re.sub(r'\bp1\b','(2*p0-p1)',text)
            if case['id']=='T19':text=re.sub(r'\bp1\b','(p0-p1)',text)
            return R.parse(text)
        P={}
        for ij,text in source['homogeneous_bivector'].items():
            i,j=map(int,ij.split(','));P[i,j]=read(text);P[j,i]=R.scale(P[i,j],-1)
        D={}
        for missing in range(5):
            i,j,k,l=[a for a in range(5)if a!=missing]
            pf=R.add(R.mul(P.get((i,j),{}),P.get((k,l),{})),
                R.scale(R.mul(P.get((i,k),{}),P.get((j,l),{})),-1),
                R.mul(P.get((i,l),{}),P.get((j,k),{})))
            D=R.add(D,R.scale(R.mul(R.variable(f'x{missing}'),pf),(-1)**missing))
        assert D==read(case['actual_quintic'])and D
        product=read(case['factor_coefficient']);assert product
        degrees=[]
        for item in case['factors']:
            factor=read(item['factor']);assert factor
            actual={sum(ex[-5:])for ex in factor}
            assert actual=={item['coordinate_degree']}and item['coordinate_degree']>0
            for _ in range(item['multiplicity']):product=R.mul(product,factor);degrees.append(item['coordinate_degree'])
        assert D==product and sum(degrees)==5 and len(degrees)>1
        records.append(dict(id=case['id'],whole_parameter_signed_Pfaffian_reconstructed=True,
            exact_full_factor_product=True,factor_degrees=degrees,
            period_reparametrization='p1_old=2*p0-p1_new'if case['id']=='T18'else'p1_old=p0-p1_new'if case['id']=='T19'else'none'))
    return records

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    names=['reports/CLASSIFICATION_GA9_DISTINCT_CURVE_20261107.json',
           'reports/CLASSIFICATION_GA9_DISTINCT_CURVE_ORBIT_20261107.json',
           'reports/CLASSIFICATION_GA9_DISTINCT_COMPONENT_AUDIT_20261107.json',
           'reports/CLASSIFICATION_MPS23_DIVISOR_FACTORS_20261107.json',
           'reports/CLASSIFICATION_GA9_DISTINCT_PRIOR98_20261107.json',
           'reports/CLASSIFICATION_GA9_DISTINCT_WEIGHT_BLOCKS_20261107.json',
           'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_SHORT1_DISTINCT_C_20261107.json',
           'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json',
           'reports/CLASSIFICATION_NIL211_SHORT_TRANSLATION_REVIEW_20261107.json',
           'research/classification/NIL211_SHORT_TRANSLATION_REVIEW.md',
           'research/notes/MPS_S61_CATALOGUE_20261103_v62.md']
    inputs=[ROOT/n for n in names]
    curve_data,orbit_data,audit_data,mps,prior,frame,source,original,oldreview=[json.loads(p.read_text())for p in inputs[:9]]
    assert prior['component_count']==len(prior['components'])==98
    curve=decode(curve_data['whole_actual_cubic_Poisson_curve'])
    assert all(not p for p in q.jacobi(curve).values())
    P0=q.antisymmetric({ij:at_parameter(curve.get(ij,{}),F(0))for ij in q.PAIRS})
    P1=q.antisymmetric({ij:at_parameter(curve.get(ij,{}),F(1))for ij in q.PAIRS})
    assert all(P1.get(ij,{})==decode(curve_data['nonzero_point_brackets']).get(ij,{})for ij in q.PAIRS)
    assert all(P0.get(ij,{})==decode(frame['canonical_point']).get(ij,{})for ij in q.PAIRS)
    weights=[0,1,0,1,1]
    assert all(weights[i]+weights[j]-sum(weights[k]*ex[k]for k in range(5))==1-ex[5]
               for(i,j),p in curve.items()if i<j for ex in p)
    Dcurve=quintic(curve);D0=quintic(P0);D1=quintic(P1)
    assert Dcurve==parse(curve_data['actual_anticanonical_quintic'])
    assert D1==parse(curve_data['nonzero_point_actual_quintic'])
    assert D0==parse('3*h*(3*(h*h-y*y)*(h*x-y*z)-2*h**3*w)')
    # Verify actual P0 is the zero five-parameter canonical source point.
    row=original['full_nilpotent_sources'][8]
    S=LaurentRing([f'n{i}'for i in range(28)],['x','y','z','w'])
    fixed=source['nil_short_translation_distinct_c']['canonical_five_parameter_source']
    assert {int(k[1:])for k in fixed}==set(range(28))
    values={}
    for key,value in fixed.items():
        p=S.parse(re.sub(r'\b(D[1-4]|d0)\b','0',value));assert not p or set(p)=={S.zero}
        values[int(key[1:])]=p.get(S.zero,F())
    for ij,text in zip(it.combinations(range(1,5),2),row['brackets']):
        p=S.parse(text);homogeneous={}
        for ex,c in p.items():
            coefficient=3*c
            for i,n in enumerate(ex[:28]):coefficient*=values[i]**n
            if coefficient:
                mon=(2-sum(ex[-4:]),)+ex[-4:]+(0,0,0)
                homogeneous[mon]=homogeneous.get(mon,F())+coefficient
        assert {e:c for e,c in homogeneous.items()if c}==P0.get(ij,{})
    # Entire full150/350 derivative and actual group/vertical/scalar directions.
    Dcols=[]
    for ij,ex in q.COLS:
        change=dict(P1);change[ij]=q.add(P1.get(ij,{}),{ex+(0,0,0):F(1)})
        change[ij[::-1]]=q.scale(change[ij],-1)
        Dcols.append(q.jvec(q.jacobi(change)))
    derivative=q.transpose(Dcols)
    vertical,Gcols=actions(P1);G=q.transpose(Gcols)
    assert q.rank(q.transpose(vertical))[0]==24 and q.rank(G)[0]==48
    assert all(not sum((row[i]*v for i,v in enumerate(col)if v),F())for row in derivative for col in Gcols)
    DC=orbit_data['full_D_minor'];GC=orbit_data['full_G_minor']
    detD=q.determinant([[derivative[i][j]for j in DC['columns']]for i in DC['rows']])
    detG=q.determinant([[G[i][j]for j in GC['columns']]for i in GC['rows']])
    assert len(DC['rows'])==102 and detD==F(DC['determinant'])and detD
    assert len(GC['rows'])==48 and detG==F(GC['determinant'])and detG
    V0,M0=actions(P0);assert q.rank(q.transpose(M0))[0]==47
    # Exact invertible eigen-coordinate model, canonical Euler correction.
    B=[[F(v)for v in r]for r in audit_data['old_coordinates_from_eigen_coordinates']]
    Binv=q.inverse(B);assert q.matmul(B,Binv)==[[F(i==j)for j in range(5)]for i in range(5)]
    oldD=parse(orbit_data['actual_quintic_primitive'])
    primitive=modelparse(audit_data['canonical_eigen_coordinate_quintic'])
    assert q.substitute(oldD,linear_substitution(B))==q.scale(primitive,F(audit_data['eigen_quintic_original_scalar']))
    div=[q.add(*(q.der(P1.get((i,j),{}),i)for i in range(5)))for j in range(5)]
    canonical=q.antisymmetric({(i,j):q.add(P1.get((i,j),{}),
           q.scale(q.add(q.mul(q.atom(1,i),div[j]),q.scale(q.mul(div[i],q.atom(1,j)),-1)),F(-1,5)))for i,j in q.PAIRS})
    assert all(canonical.get(ij,{})==decode(orbit_data['canonical_divergence_free_point']).get(ij,{})for ij in q.PAIRS)
    pushed=q.antisymmetric({ij:q.substitute(q.add(*(q.scale(canonical.get((a,b),{}),Binv[ij[0]][a]*Binv[ij[1]][b])
                   for a,b in it.product(range(5),repeat=2))),linear_substitution(B))for ij in q.PAIRS})
    model=q.antisymmetric({tuple(map(int,k.split(','))):modelparse(v)for k,v in audit_data['canonical_eigen_coordinate_bivector'].items()})
    scaling=F(audit_data['bivector_original_scaling'])
    assert all(model.get(ij,{})==q.scale(pushed.get(ij,{}),scaling)for ij in q.PAIRS)
    assert quintic(model)==q.scale(primitive,60)
    # All-field Gauss proof: linear in u4, only possible coefficient factors u0,u1.
    coefficient=q.der(primitive,4);assert coefficient==modelparse('-240*u0**2*u1**2')
    assert not q.der(coefficient,4)
    remainder=q.add(primitive,q.scale(q.mul(q.atom(1,4),coefficient),-1))
    lin=[q.atom(1,i)for i in range(8)];lin[0]={}
    assert q.substitute(remainder,lin)==modelparse('-u1**4*u2')
    lin=[q.atom(1,i)for i in range(8)];lin[1]={}
    assert q.substitute(remainder,lin)==modelparse('-80*u0**3*u3**2')
    # Exact reduced singular support: d/du4 forces u0*u1=0;
    # u0=0 then d/du2 forces u1=0; u1=0,u0!=0 then d/du3 forces u3=0.
    assert q.der(primitive,2)==modelparse('160*u0**2*u1*u3-u1**4')
    assert q.der(primitive,3)==modelparse('-160*u0**3*u3+160*u0**2*u1*u2+6*u0*u1**3')
    for zeros in [(0,1),(1,3)]:
        lin=[q.atom(1,i)for i in range(8)]
        for i in zeros:lin[i]={}
        assert not q.substitute(primitive,lin)
        assert all(not q.substitute(q.der(primitive,i),lin)for i in range(5))
    mps_records=mps_factors(mps,prior,inputs)
    mps_ids={r['id']for r in mps_records};prior_audit=[]
    for record in prior['components']:
        cid=record['id']
        if record['generic_rank']!=4:reason='different generic rank'
        elif record['projective_dimension']!=23:reason='different full component dimension'
        elif cid in mps_ids:reason='whole accepted dense MPS family has verified reducible actual quintic'
        else:
            geom=record['generic_actual_degeneracy_divisor']
            degrees,multiplicities=geom['irreducible_factor_degrees'],geom['corresponding_multiplicities']
            assert sum(a*b for a,b in zip(degrees,multiplicities))==5 and sum(multiplicities)>1
            assert record['proof_sources']and record['verification_records']
            reason='accepted nonMPS generic divisor factorization is reducible/nonreduced'
        prior_audit.append(dict(id=cid,reason=reason))
    assert len(prior_audit)==98 and len(mps_records)==29
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_CAS_import=True,no_producer_import=True,
        full_literal_cubic_curve_all350_Jacobi_zero=True,actual_Cstar_all_nonzero_orbit_identity=True,
        P0_exactly_the_accepted_full_Ga9_source_point=True,P0_actual_projective_orbit_dimension=22,
        full150_350_D_rank=102,vertical_rank=24,full_actual_group_vertical_scale_rank=48,
        actual_projective_orbit_dimension=23,full_projective_tangent_dimension=23,
        independent_D_minor=dict(rows=DC['rows'],columns=DC['columns'],determinant=str(detD)),
        independent_G_minor=dict(rows=GC['rows'],columns=GC['columns'],determinant=str(detG)),
        full_D_times_all_actual_orbit_directions_zero=True,smooth_unique_23_orbit_component_proved=True,
        exact_eigen_coordinate_change_and_all9_bivector_entries_checked=True,
        primitive_quintic=audit_data['canonical_eigen_coordinate_quintic'],
        absolute_irreducibility='Primitive linearu4 with coefficient -240u0^2u1^2; remainder not divisible byu0 oru1 by exact restrictions, so Gauss over C proves irreducibility.',
        reduced_singular_support=['u0=u1=0','u1=u3=0'],
        nonnormality_directly_proved_without_native_singular_ideal=True,
        all29_whole_generic_MPS_factorizations=mps_records,
        all98_distinctness_audit=prior_audit,distinct_from_all_prior98=True,
        whole_Ga9_distinct_source_contained='All nonzero actualcurve points are in the P1orbit, so its component closure contains P0. Accepted irreducible source saturation is orbitclosure(P0), hence fully contained.',
        componenthood_and_distinctness_supported=True,global_registers_or_count_not_edited=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','componenthood_and_distinctness_supported','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
