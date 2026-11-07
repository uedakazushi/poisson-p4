#!/usr/bin/env python3
"""Independent whole-source type7 curves, known99 map and new23 audit.

No producer/CAS imports; all150/350 calculations use sparse Fraction helpers.
"""
import argparse
from collections import Counter
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
    aliases={'h':'h','x':'x1','y':'x2','z':'x3','w':'x4','t':'tau0','u':'tau1'}
    return q.parse(re.sub(r'\b(h|x|y|z|w|t|u)\b',lambda m:aliases[m[1]],text))

def decode(data):return q.antisymmetric({tuple(map(int,k.split(','))):parse(v)for k,v in data.items()})
def modeldecode(data):
    aliases=['h','x','y','z','w']
    return decode({k:re.sub(r'\bu([0-4])\b',lambda m:aliases[int(m[1])],v)for k,v in data.items()})

def at_t(p,value):
    out={}
    for ex,c in p.items():
        co=c*value**ex[5]
        if co:
            e=list(ex);e[5]=0;e=tuple(e);out[e]=out.get(e,F())+co
    return {e:c for e,c in out.items()if c}

def at_zero_coordinate(p,i):return {e:c for e,c in p.items()if not e[i]}

def quintic(P):
    out={}
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
        out=q.add(out,q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
    return out

def pf4(P):return q.add(q.mul(P.get((1,2),{}),P.get((3,4),{})),q.scale(q.mul(P.get((1,3),{}),P.get((2,4),{})),-1),q.mul(P.get((1,4),{}),P.get((2,3),{})))

def div(P):return [q.add(*(q.der(P.get((i,j),{}),i)for i in range(5)))for j in range(5)]
def canonical(P):
    d=div(P)
    out=q.antisymmetric({(i,j):q.add(P.get((i,j),{}),q.scale(q.add(q.mul(q.atom(1,i),d[j]),q.scale(q.mul(q.atom(1,j),d[i]),-1)),-F(1,5)))for i,j in q.PAIRS})
    assert all(not p for p in div(out))
    return out

def poly_determinant(A):
    out={}
    for perm in it.permutations(range(len(A))):
        term=q.atom((-1)**sum(perm[i]>perm[j]for i in range(len(A))for j in range(i+1,len(A))))
        for i,j in enumerate(perm):term=q.mul(term,A[i][j])
        out=q.add(out,term)
    return out

def linear_substitution(M):return [q.add(*(q.atom(v,j)for j,v in enumerate(row)))for row in M]+[q.atom(1,i)for i in range(5,8)]

def push(P,M):
    inv=q.inverse(M);changed={ij:q.substitute(p,linear_substitution(inv))for ij,p in P.items()}
    return q.antisymmetric({(i,j):q.add(*(q.scale(changed.get((a,b),{}),M[i][a]*M[j][b])for a,b in it.product(range(5),repeat=2)))for i,j in q.PAIRS})

def action(P,a,b,positive=True):
    return q.antisymmetric({(i,j):q.scale(q.add(P.get((b,j),{})if i==a else{},P.get((i,b),{})if j==a else{},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1)),1 if positive else -1)for i,j in q.PAIRS})

def group_columns(P):
    vertical=[];orbit=[]
    for a,b in it.product(range(5),repeat=2):
        vertical.append(q.vec({(i,j):q.add(q.atom(1,i,b)if j==a else{},q.atom(-1,j,b)if i==a else{})for i,j in q.PAIRS}))
        orbit.append(q.vec(action(P,a,b,False)))
    return vertical,vertical+orbit+[q.vec(P)]

def derivative(P):
    columns=[]
    assert all(not j for j in q.jacobi(P).values())
    for ij,ex in q.COLS:
        inc=dict(P);inc[ij]=q.add(inc.get(ij,{}),{ex+(0,0,0):F(1)});inc[ij[::-1]]=q.scale(inc[ij],-1)
        columns.append(q.jvec(q.jacobi(inc)))
    return q.transpose(columns)

def check_minor(A,cert):
    rows,cols=cert['rows'],cert['columns']
    value=q.determinant([[A[i][j]for j in cols]for i in rows]);assert value==F(cert['determinant'])and value
    return dict(rank=len(rows),determinant=str(value),rows=rows,columns=cols)

def stabilizer(P,saved=None):
    P=canonical(P);columns=[q.vec(action(P,i,j))for i,j in it.product(range(5),repeat=2)]+[[-v for v in q.vec(P)]]
    A=q.transpose(columns);r,piv=q.rank(A);assert r==24
    if saved:
        kernel=[]
        for col in saved['complete_stabilizer_kernel_sparse']:
            v=[F(col.get(str(i),'0'))for i in range(26)]
            assert all(sum(a*b for a,b in zip(row,v))==0 for row in A);kernel.append(v)
        assert len(kernel)==2 and q.rank(q.transpose(kernel))[0]==2
        L=[[F(v)for v in row]for row in saved['stabilizer_matrix']];kappa=F(saved['stabilizer_tensor_character'])
        v=[x for row in L for x in row]+[kappa]
        assert all(sum(a*b for a,b in zip(row,v))==0 for row in A)
    return dict(full150by26_rank=24,complete_GL_scalar_kernel_dimension=2,complete_projective_stabilizer_dimension=1)

def anchor_from_original(row,family,scale):
    R=LaurentRing([f'n{i}'for i in range(20)],['x','y','z','w'])
    values={int(k[1:]):F(v)for k,v in family['canonical_source_values'].items()}
    assert set(values)==set(range(20));out={}
    for ij,text in zip(it.combinations(range(1,5),2),row['brackets']):
        p={}
        for ex,c in R.parse(text).items():
            co=scale*c
            for i,power in enumerate(ex[:20]):co*=values[i]**power
            if co:
                assert sum(ex[-4:])<=2
                e=(2-sum(ex[-4:]),)+ex[-4:]+(0,0,0);p[e]=p.get(e,F())+co
        out[ij]={e:c for e,c in p.items()if c}
    return q.antisymmetric(out)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    names=['reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_CURVE_RECOGNITION_20261107.json',
      'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_ORBIT_MODELS_20261107.json',
      'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES_20261107.json',
      'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json',
      'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_REVIEW_V2_20261107.json',
      'reports/CLASSIFICATION_GA9_DISTINCT_COMPONENT_AUDIT_20261107.json',
      'reports/CLASSIFICATION_GA9_DISTINCT_COMPONENT_REVIEW_20261107.json',
      'reports/CLASSIFICATION_GA9_DISTINCT_PRIOR98_20261107.json',
      'research/classification/HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT.md',
      'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT_V2_20261107.json',
      'research/classification/HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_CURVE_RECOGNITION.md']
    paths=[ROOT/n for n in names]
    producer,orbits,families,original,cover,known99,priorreview,prior98=[json.loads(p.read_text())for p in paths[:8]]
    claims=producer['nil22_long_quartic_curve_recognition'];orbits=orbits['nil22_long_quartic_orbit_models']
    families=families['nil22_long_quartic_sources'];row=original['full_nilpotent_sources'][6]
    assert cover['all_passed']and cover['projective_saturated_dimension_upper_bound']==22
    assert priorreview['all_passed']and priorreview['distinct_from_all_prior98']
    records=[];models={};points={}
    for saved in claims['all_five_actual_curve_recognitions']:
        index=saved['anchor'];anchor=orbits['all_seven_full_ambient_anchors'][index-1]
        family=families['all_seven_complete_families'][index-1]
        curve=decode(saved['literal_curve_brackets']);P0=q.antisymmetric({ij:at_t(curve.get(ij,{}),F())for ij in q.PAIRS})
        assert all(not p for p in q.jacobi(curve).values())
        expected=anchor_from_original(row,family,F(anchor['integer_scale']))
        assert all(P0.get(ij,{})==expected.get(ij,{})==decode(anchor['witness_brackets']).get(ij,{})for ij in q.PAIRS)
        weights=saved['Cstar_coordinate_weights'];mu=saved['curve_parameter_weight']
        assert weights==family['canonical_Cstar_weights']
        assert all(weights[i]+weights[j]-sum(weights[k]*ex[k]for k in range(5))==1-mu*ex[5]
                   for(i,j),poly in curve.items()if i<j for ex in poly)
        vertical,G0cols=group_columns(P0);G0=q.transpose(G0cols)
        assert q.rank(q.transpose(vertical))[0]==24 and q.rank(G0)[0]==47
        anchor_minor=check_minor(G0,anchor['orbit_vertical_scale'])
        record=dict(anchor=index,all350_literal_curve_Jacobi=True,constant_matches_full20_source_and_anchor=True,
            exact_Cstar_pushforward_all_monomials=True,nonzero_parameter_weight=mu,
            anchor_full_group_rank=47,anchor_projective_orbit_dimension=22,anchor_minor=anchor_minor,
            whole_source_orbitclosure_bridge='Accepted irreducible fullsource saturation has dimension<=22 and contains this actual orbit22; hence equals its orbitclosure.')
        if index in [1,2,7]:
            assert all(not curve.get((0,j),{})for j in range(1,5))
            Z=div(curve);N=[[q.der(z,j)for j in range(5)]for z in Z]
            assert all(z==q.add(*(q.mul(N[i][j],q.atom(1,j))for j in range(5)))for i,z in enumerate(Z))
            assert Z==list(map(parse,saved['cofactor']))
            characteristic=poly_determinant([[q.add(q.atom(1,6)if i==j else{},q.scale(N[i][j],-1))for j in range(5)]for i in range(5)])
            expression={1:'u*(t-u)*(t+u)*(2*t-u)*(2*t+u)',2:'u*(t-u**2)*(4*t-u**2)',7:'u*(4*t+u**2)*(16*t+u**2)'}[index]
            assert characteristic==parse(expression)==parse(saved['cofactor_characteristic'])
            quartic=pf4(curve);assert quartic==parse(saved['actual_quartic_pfaffian'])
            assert at_t(quartic,F())==pf4(P0)
            assert any(sum(ex[:5])==4 for ex in quartic)
            record.update(h_Casimir_and_full_cofactor_verified=True,exact_characteristic=expression,
                generic_regular_semisimple_type=12,anchor_quartic_irreducibility_dependency='All-parameter primitive/odd-power certificates independently verified in accepted full-seven-source V2.',
                whole_source_containment='finite known closed MPS union via accepted entire type12 theorem')
        else:
            point=decode(saved['nonzero_point_brackets']);scale=F(saved['nonzero_point_integer_scale'])
            assert all(point.get(ij,{})==q.scale(at_t(curve.get(ij,{}),F(1)),scale)for ij in q.PAIRS)
            D=derivative(point);V,Gcols=group_columns(point);G=q.transpose(Gcols)
            assert q.rank(q.transpose(V))[0]==24 and q.rank(G)[0]==48
            assert all(not sum((row[j]*col[j]for j in range(150)if row[j]and col[j]),F())for row in D for col in Gcols)
            dminor=check_minor(D,saved['full_derivative']);gminor=check_minor(G,saved['actual_group_vertical_scalar'])
            assert dminor['rank']==102 and gminor['rank']==48
            can=canonical(point);assert all(can.get(ij,{})==decode(saved['canonical_divergence_free_brackets']).get(ij,{})for ij in q.PAIRS)
            stab=stabilizer(can,saved)
            B=[[F(v)for v in row]for row in saved['old_coordinates_from_eigen_coordinates']];Binv=q.inverse(B)
            L=[[F(v)for v in row]for row in saved['stabilizer_matrix']]
            diagonal=q.matmul(q.matmul(Binv,L),B);eigen=list(map(F,saved['stabilizer_eigenvalues']))
            assert diagonal==[[eigen[i]if i==j else F()for j in range(5)]for i in range(5)]
            model=decode(saved['primitive_eigen_tensor']);eigenscale=F(saved['eigen_tensor_scale'])
            moved=push(can,Binv)
            assert all(model.get(ij,{})==q.scale(moved.get(ij,{}),eigenscale)for ij in q.PAIRS)
            Dmodel=quintic(model)
            assert Dmodel==q.scale(parse(saved['primitive_eigen_quintic']),F(saved['eigen_quintic_content']))
            assert all(not p for p in q.jacobi(model).values())
            models[index]=model;points[index]=point
            record.update(full150_350_D_minor=dminor,full_group_minor=gminor,full_D_all_group_directions_zero=True,
                vertical_rank=24,D_rank=102,group_rank=48,projective_orbit_and_tangent_dimension=23,
                smooth_unique_23_orbitclosure_component=True,full150_Euler_canonical_and_eigen_tensor_equal=True,
                complete_stabilizer=stab,semisimple_stabilizer_eigenvalues=list(map(str,eigen)),
                stabilizer_eigenvalue_multiplicities=sorted(Counter(eigen).values(),reverse=True))
            if index==6:
                known=modeldecode(known99['canonical_eigen_coordinate_bivector'])
                total=[[F(v)for v in row]for row in saved['exact_map_to_known99']]
                assert total==q.matmul([[F(1)if i==j and i<2 or i==j==3 else F(-1)if i==j==2 else F(-1,48)if i==j==4 else F()for j in range(5)]for i in range(5)],Binv)
                moved_point=push(point,total);moved_can=canonical(moved_point)
                scalar=F(saved['exact_scalar_to_known99'])
                assert scalar==F(-125,2)
                assert all(moved_can.get(ij,{})==q.scale(known.get(ij,{}),scalar)for ij in q.PAIRS)
                # This is a full150 identity with the explicitly reconstructed Euler vertical correction.
                d=div(moved_point)
                assert all(q.add(moved_point.get((i,j),{}),q.scale(known.get((i,j),{}),-scalar))==q.scale(q.add(q.mul(q.atom(1,i),d[j]),q.scale(q.mul(q.atom(1,j),d[i]),-1)),F(1,5))for i,j in q.PAIRS)
                record.update(exact_known99_full150_map_scalar_Euler_identity=True,known99_scalar=str(scalar),
                    whole_source_containment='R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23')
        records.append(record)
    assert {r['anchor']for r in records}=={1,2,5,6,7}
    # New case5 nine-bracket normal form and complex Gauss certificate.
    expected5=decode({'0,2':'-6*h**2','0,3':'3*h*x','0,4':'-3*h*y+162*h*z-45*x*y',
       '1,2':'9*h*x','1,3':'-2*x**2','1,4':'-5*h*z+2*x*y-243*x*z','2,3':'-9*h*z+3*x*y',
       '2,4':'6*h*w-3*y**2+567*y*z','3,4':'2*x*w+7*y*z-243*z**2'})
    assert all(models[5].get(ij,{})==expected5.get(ij,{})for ij in q.PAIRS)
    F5=parse('h**3*z**2+h**2*w*x**2+h**2*x*y*z+27*h*x**2*y*z+3*x**3*y**2')
    assert quintic(models[5])==q.scale(F5,F(75))
    assert q.der(F5,4)==parse('h**2*x**2')
    remainder=at_zero_coordinate(F5,4)
    assert at_zero_coordinate(remainder,0)==parse('3*x**3*y**2')
    assert at_zero_coordinate(remainder,1)==parse('h**3*z**2')
    singular_planes=[(0,1),(0,2),(1,3)]
    for pair in singular_planes:
        for i in range(5):assert not at_zero_coordinate(at_zero_coordinate(q.der(F5,i),pair[0]),pair[1])
    assert at_zero_coordinate(q.der(F5,2),0)==parse('6*x**3*y')
    assert at_zero_coordinate(q.der(F5,3),1)==parse('2*h**3*z')
    # The known99 dense orbit has a unique projective stabilizer with five simple eigenvalues.
    known=modeldecode(known99['canonical_eigen_coordinate_bivector']);knownstab=stabilizer(known)
    known_eigen=list(map(F,known99['eigenvalues']));assert len(set(known_eigen))==5
    knownL=[[known_eigen[i]if i==j else F()for j in range(5)]for i in range(5)]
    known_action={ij:q.add(*(q.scale(action(known,a,b).get(ij,{}),knownL[a][b])for a,b in it.product(range(5),repeat=2)))for ij in q.PAIRS}
    assert all(known_action[ij]==q.scale(known.get(ij,{}),F(50))for ij in q.PAIRS)
    assert next(r for r in records if r['anchor']==5)['stabilizer_eigenvalue_multiplicities']==[2,2,1]
    # Transfer the already accepted entire-source reducibility audit, without rerunning its29 factorization proofs.
    audit=priorreview['all98_distinctness_audit'];assert len(audit)==98
    assert {c['id']for c in audit}=={c['id']for c in prior98['components']}
    reasons=Counter(c['reason']for c in audit)
    assert sorted(reasons.values())==[6,27,29,36]
    assert set(reasons)=={'different full component dimension','different generic rank',
        'whole accepted dense MPS family has verified reducible actual quintic',
        'accepted nonMPS generic divisor factorization is reducible/nonreduced'}
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
        all_five_whole_source_curve_reviews=records,
        case5_exact_nine_bracket_normal_form=True,case5_actual_quintic='75*(h³z²+h²wx²+h²xyz+27hx²yz+3x³y²)',
        case5_complex_absolute_irreducibility='primitive linear in w; coefficient h²x²; remainder modh=3x³y² and modx=h³z². Hence gcd1 in C[h,x,y,z], and Gauss applies over C.',
        case5_reduced_singular_support=['h=x=0','h=y=0','x=z=0'],case5_nonnormality=True,
        known99_complete_projective_stabilizer=knownstab,known99_simple_semisimple_stabilizer_eigenvalues=list(map(str,known_eigen)),
        case5_distinct_from_known99='Both have dense actual23 orbits and unique projective stabilizer lines. Multiplicities221 versus11111 are invariant under conjugacy/nonzero scalar and scalar identity shifts. Different reduced singular supports3planes versus2planes also independently distinguish the divisors.',
        inherited_prior98_full_audit=audit,case5_distinct_from_all_prior99=True,
        case5_new_smooth23_component_supported=True,new_global_count_or_ID_not_edited=True,
        whole_MPS_source_indices=[1,2,7],whole_known99_source_indices=[6],whole_new23_source_indices=[5],
        proof_scope='Whole-source closure bridge uses the previously independently proved finite seven-source theorem and actual orbit22 checks; generic curve closure and fullD102/G48 suffice. No IFT, private jet, formal slice, or assumptions on generic/all-fiber equivalence.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','case5_distinct_from_all_prior99','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
