#!/usr/bin/env python3
"""Independent type4 whole-quartic MPS recognition; no CAS/producer import."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

from classification_quadratic_cofactor_review import verify_source,PAIRS as AP
from classification_minimal_quartic_review import LaurentRing,determinant
from classification_quartic_nil31_short_review import core_completeness
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]
LOWER=[3,4,8,12,13,14];MAG=[0,1,2,5,6,7]

def power(R,p,n):
    if n<0:p=R.inverse_monomial(p);n=-n
    out=R.constant(1)
    for _ in range(n):out=R.mul(out,p)
    return out

def sub(R,p,values):
    out={}
    for ex,co in p.items():
        term=R.constant(co)
        for i,n in enumerate(ex):
            if n:term=R.mul(term,power(R,values.get(i,R.variable(R.names[i])),n))
        out=R.add(out,term)
    return out

def linear_system(R,eqs,variables):
    matrix=[];rhs=[]
    for p in eqs:
        row=[{}for _ in variables];constant={}
        for ex,co in p.items():
            degree=sum(ex[i]for i in variables);assert degree<=1
            if not degree:
                constant=R.add(constant,{ex:co});continue
            j=next(j for j,i in enumerate(variables)if ex[i])
            e=list(ex);e[variables[j]]=0
            row[j]=R.add(row[j],{tuple(e):co})
        matrix.append(row);rhs.append(R.scale(constant,-1))
    return matrix,rhs

def normalize_numerator(R,p):
    if not p:return ()
    exmin=[min(ex[i]for ex in p)for i in range(len(R.names))]
    value={tuple(n-min(0,m)for n,m in zip(ex,exmin)):c for ex,c in p.items()}
    lead=max(value);den=value[lead]
    return tuple(sorted((ex,c/den)for ex,c in value.items()))

def saved_parse(text):
    aliases={'x':'x1','y':'x2','z':'x3','w':'x4','e':'tau0','u':'tau1'}
    return q.parse(re.sub(r'\b(x|y|z|w|e|u)\b',lambda m:aliases[m[1]],text))

def poly_det(M):
    out={}
    for perm in it.permutations(range(len(M))):
        value=q.atom((-1)**sum(perm[i]>perm[j]for i in range(len(M))for j in range(i+1,len(M))))
        for i,j in enumerate(perm):value=q.mul(value,M[i][j])
        out=q.add(out,value)
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    paths=[ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json',
      ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_SUPPORT_20261107.json',
      *[ROOT/f'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_{name}_20261107.json'for name in ['SOURCES','SUPPORT','AMBIENT','COMPONENT']],
      ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_RECOGNITION_20261107.json']
    full,corefile,sources,support,ambient,component,final=[json.loads(p.read_text())for p in paths]
    assert hashlib.sha256(final['source_bytes'].encode()).hexdigest()==final['source_sha256']
    source=full['full_nilpotent_sources'][3];original_source=verify_source(source,4)
    R=LaurentRing([f'n{i}'for i in range(18)],['x','y','z','w'],invert=['n17','n15','n9'])
    parse=R.parse
    J=list(map(parse,corefile['all_reduced_jacobi_equations']))
    assert {normalize_numerator(R,p)for p in J}=={normalize_numerator(R,parse(p))for p in source['distinct_nonzero_jacobi_coefficients']}
    pf=parse(source['actual_affine_pfaffian']);top=parse(source['highest_quartic_pfaffian'])
    complete_core=core_completeness()
    zero={17:{},10:{},15:{},9:{}}
    for fixed,first in [(zero|{16:parse('1/3')},'(3*n11+1)*n8/6'),
                         (zero|{11:parse('-n16')},'(1-3*n16)*n8/9')]:
        assert sub(R,J[2],fixed)==parse(first)
        assert sub(R,J[8],fixed|{8:{}})==sub(R,parse('(n11+1/4)*n7'),fixed)
        assert sub(R,top,fixed|{8:{}})==sub(R,parse('n11*n7*z**4'),fixed)
    # The core report's seven nonzero-parameter branches plus all three
    # continuous-core quartic exceptions are exactly the ten candidates.
    expected=[row['substitution']for row in corefile['whole_source_dimension_bounds']if row['parameter_nonzero']]
    expected += [{f'n{i}':str(v)for i,v in {17:0,10:0,15:0,9:0,11:a,16:b}.items()}
                 for a,b in [(F(-1,3),F(1,3)),(F(-1,4),F(1,3)),(F(-1,4),F(1,4))]]
    rows=sources['nil31_short_quartic_sources']['all_ten_candidate_sources'];assert len(rows)==10
    records=[];allfamilies=[]
    for index,(claim,expectedcore)in enumerate(zip(rows,expected)):
        core={int(k[1:]):parse(v)for k,v in claim['core_substitution'].items()}
        assert core=={int(k[1:]):parse(v)for k,v in expectedcore.items()}
        values={int(k[1:]):parse(v)for k,v in claim['source_substitution'].items()}
        assert set(values)==set(range(18))
        for i,p in core.items():assert sub(R,p,values)==values[i]
        allowed={int(v.split('!')[0][1:])for v in claim['localization']}
        assert all(not any(ex[i]<0 for i in range(18)if i not in allowed)for p in values.values()for ex in p)
        residual=[sub(R,p,core)for p in J]
        eqL=[residual[i]for i in [2,3,4]]
        ML,bL=linear_system(R,eqL,LOWER)
        LC=claim['lower_linear_certificate'];lr,lc=LC['rows'],LC['columns']
        assert ML==[[parse(v)for v in row]for row in LC['matrix']] and bL==list(map(parse,LC['rhs']))
        minorL=determinant(R,[[ML[i][j]for j in lc]for i in lr]);assert minorL==parse(LC['determinant'])
        assert len(minorL)==1 and next(iter(minorL.values()))
        assert all(not sub(R,eqL[i],values)for i in lr)
        assert all(values[LOWER[j]]==R.variable(f'n{LOWER[j]}')for j in range(6)if j not in lc)
        lower_values={i:values[i]for i in LOWER}
        eqM=[sub(R,residual[i],lower_values)for i in [8,9,10,11]]
        MM,bM=linear_system(R,eqM,MAG)
        MC=claim['magnetic_linear_certificate'];mr,mc=MC['rows'],MC['columns']
        assert MM==[[parse(v)for v in row]for row in MC['matrix']]and bM==list(map(parse,MC['rhs']))
        minorM=determinant(R,[[MM[i][j]for j in mc]for i in mr]);assert minorM==parse(MC['determinant'])
        assert len(minorM)==1 and next(iter(minorM.values()))
        assert all(not sub(R,eqM[i],values)for i in mr)
        assert all(values[MAG[j]]==R.variable(f'n{MAG[j]}')for j in range(6)if j not in mc)
        # Every pivot monomial is a unit on exactly its named localization.
        for value in [minorL,minorM]:
            ex=next(iter(value))
            assert all(not n for i,n in enumerate(ex)if i not in allowed)
        jac=[sub(R,p,values)for p in J]
        actual={normalize_numerator(R,p)for p in jac if p}
        constraints={normalize_numerator(R,parse(p))for p in claim['all_remaining_Jacobi_constraints']}
        assert actual==constraints
        brackets=[sub(R,parse(p),values)for p in source['brackets']]
        assert brackets==list(map(parse,claim['universal_brackets']))
        actualpf=sub(R,pf,values);actualtop=sub(R,top,values)
        assert actualpf==parse(claim['actual_affine_pfaffian']) and actualtop==parse(claim['highest_quartic_pfaffian'])
        allfamilies.append((values,brackets,actualpf))
        records.append(dict(candidate=index,all18_substitutions_verified=True,
            two_uniform_unit_pivots=[R.encode(minorL),R.encode(minorM)],
            all_remaining_compatibility_exact=claim['all_remaining_Jacobi_constraints'],
            all_six_brackets_and_actual_and_top_pfaffian_verified=True))
    supp=support['nil31_short_quartic_support']
    for exclusion in supp['all_nine_excluded_candidate_certificates']:
        i=exclusion['candidate_index']
        if i==6:assert not parse(rows[i]['highest_quartic_pfaffian']);continue
        factor,quotient=parse(exclusion['divisor']),parse(exclusion['quotient'])
        assert R.mul(factor,quotient)==allfamilies[i][2]
        assert any(any(ex[-4:])for ex in factor)
        assert all(sum(ex[-4:])>=0 for ex in quotient)
    values,family,familyPf=allfamilies[3]
    assert values=={int(k[1:]):parse(v)for k,v in supp['source_substitution'].items()}
    free={12,13,17,6,7,8}
    assert set(i for p in values.values()for ex in p for i,n in enumerate(ex[:18])if n)==free
    assert all(values[i]==R.variable(f'n{i}')for i in free)
    coeff=R.derivative(familyPf,0)
    assert not R.derivative(coeff,0)and coeff==parse('n17*(z-1/(6*n17))**3')
    restriction=sub(R,familyPf,{20:parse('1/(6*n17)')})
    assert restriction==parse('-(12*n12*n17+2*n13+6*n17*y-w)**2/(3456*n17**3)')
    # Source point and entire actual quadratic Poisson curve.
    comp=component['source_point_one_normal'];fin=final['nil31_short_quartic_recognition']
    for field in ['point','integer_scale','witness_brackets','derivative','orbit_vertical_scale']:
        assert ambient['source_point'][field]==comp[field]
    for dependency in fin['known_support_dependencies']:
        path=ROOT/dependency['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==dependency['sha256']
        paths.append(path)
    semisimple_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_SEMISIMPLE_SOURCES_20261107.json'
    assert json.loads(semisimple_path.read_text())['full_semisimple_sources'][11]['eigenvalues']==[-2,-1,1,2]
    paths.append(semisimple_path)
    paths.append(ROOT/'code/classification_unimodular_logarithmic_review.py')
    point_values={int(k[1:]):R.constant(F(v))for k,v in comp['point']['substitutions'].items()}
    P=q.decode({k:re.sub(r'\b(x|y|z|w)\b',lambda m:dict(x='x1',y='x2',z='x3',w='x4')[m[1]],v)for k,v in fin['P_brackets'].items()})
    scale=F(comp['integer_scale']);assert scale==432
    for (i,j),p in zip(AP,source['brackets']):
        actual=sub(R,parse(p),point_values)
        homogeneous={}
        for ex,co in actual.items():
            assert not any(ex[:18])
            mon=(2-sum(ex[-4:]),)+ex[-4:]+(0,0,0)
            homogeneous[mon]=scale*co
        assert homogeneous==P[i+1,j+1]
    RV=q.antisymmetric({tuple(map(int,k.split(','))):saved_parse(v)for k,v in fin['R_brackets'].items()})
    TV=q.antisymmetric({tuple(map(int,k.split(','))):saved_parse(v)for k,v in fin['T_brackets'].items()})
    decodevec=lambda v:[F(v.get(str(i),'0'))for i in range(150)]
    assert q.vec(RV)==decodevec(comp['independent_normal_sparse_vectors'][0])
    assert q.vec(TV)==decodevec(comp['quadratic_correction_sparse'])
    curve=q.antisymmetric({ij:q.add(P.get(ij,{}),q.mul(q.atom(1,5),RV.get(ij,{})),
              q.mul(q.atom(1,5,5),TV.get(ij,{})))for ij in q.PAIRS})
    assert all(not p for p in q.jacobi(curve).values())
    rawJ=[q.add(*(q.mul(curve.get((a,l),{}),q.der(curve.get((b,c),{}),l))
                 for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(5)))for i,j,k in it.combinations(range(5),3)]
    assert all(not p for p in rawJ)
    assert all(not curve.get((0,j),{})for j in range(5))
    Dcols=[]
    for ij,ex in q.COLS:
        variation=dict(P);variation[ij]=q.add(P.get(ij,{}),{ex+(0,0,0):F(1)})
        variation[ij[::-1]]=q.scale(variation[ij],-1)
        Dcols.append(q.jvec(q.jacobi(variation)))
    D=q.transpose(Dcols);vertical=[];orbit=[]
    for a,b in it.product(range(5),repeat=2):
        V={};O={}
        for i,j in q.PAIRS:
            V[i,j]=q.add(q.atom(1,i,b)if j==a else{},q.atom(-1,j,b)if i==a else{})
            O[i,j]=q.add(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),
                         q.scale(P.get((b,j),{}),-1)if i==a else{},
                         q.scale(P.get((i,b),{}),-1)if j==a else{})
        vertical.append(q.vec(V));orbit.append(q.vec(O))
    Mcols=vertical+orbit+[q.vec(P)];actualcols=Mcols+[q.vec(RV)]
    assert all(not sum((row[i]*v for i,v in enumerate(col)if v),F())for row in D for col in actualcols)
    M=q.transpose(Mcols);actual=q.transpose(actualcols)
    assert q.rank(q.transpose(vertical))[0]==24 and q.rank(M)[0]==47 and q.rank(actual)[0]==48
    dc=fin['derivative_minor'];detD=q.determinant([[D[i][j]for j in dc['columns']]for i in dc['rows']])
    ac=fin['actual_curve_lower_minor'];detA=q.determinant([[actual[i][j]for j in ac['columns']]for i in ac['rows']])
    assert detD==F(dc['determinant'])and detD and len(dc['rows'])==102
    assert detA==F(ac['determinant'])and detA and len(ac['rows'])==48
    Z=[q.add(*(q.der(curve.get((i,j),{}),i)for i in range(5)))for j in range(5)]
    assert Z==list(map(saved_parse,fin['cofactor']))
    N=[[q.der(Z[i],j)for j in range(5)]for i in range(5)]
    assert N==[[saved_parse(v)for v in row]for row in fin['homogeneous_cofactor_matrix']]
    char=poly_det([[q.add(q.atom(1,6)if i==j else{},q.scale(N[i][j],-1))for j in range(5)]for i in range(5)])
    assert char==saved_parse('u*(u**2+93312*e)*(u**2+373248*e)')
    weights=[0,2,1,0,1]
    assert all(weights[i]+weights[j]-sum(weights[k]*ex[k]for k in range(5))==1-2*ex[5]
               for(i,j),p in curve.items()if i<j for ex in p)
    curvePf=q.add(q.mul(curve[1,2],curve[3,4]),q.scale(q.mul(curve[1,3],curve[2,4]),-1),q.mul(curve[1,4],curve[2,3]))
    assert curvePf==saved_parse(fin['actual_quartic_pfaffian'])
    A=[[F(v)for v in row]for row in source['field_linear_matrix']]
    c=list(map(F,source['field_constant']));affinecols=[]
    for i,j in it.product(range(4),repeat=2):
        B=[[F(a==i and b==j)for b in range(4)]for a in range(4)]
        affinecols.append([sum(B[a][k]*A[k][b]-A[a][k]*B[k][b]for k in range(4))for a,b in it.product(range(4),repeat=2)]+[sum(B[a][k]*c[k]for k in range(4))for a in range(4)])
    for j in range(4):affinecols.append([F()]*16+[-A[a][j]for a in range(4)])
    assert q.rank(q.transpose(affinecols))[0]==13
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_CAS_import=True,no_producer_import=True,
        complete_original_source=original_source,complete_core=complete_core,
        all_ten_complete_quartic_candidates=records,
        all_nine_excluded_candidates_verified=True,
        remaining_Gm_times_A5_full_source_Pfaffian_Gauss_certificate=True,
        canonical432scaled_point_reconstructed_from_fullsource=True,
        literal_P_plus_eR_plus_e2T_full350_and_all_rawJacobi_zero=True,
        exact_full_D_rank=102,vertical_rank=24,orbit_vertical_scalar_rank=47,actual_curve_rank=48,
        independent_D_minor=dict(rows=dc['rows'],columns=dc['columns'],determinant=str(detD)),
        independent_actual_minor=dict(rows=ac['rows'],columns=ac['columns'],determinant=str(detA)),
        full_cofactor_and_characteristic_polynomial_verified=True,
        generic_affine_regular_semisimple_type=12,
        actual_Cstar_identity_verified=True,
        complete_source_parameter_count=6,affine_cofactor_orbit_dimension=13,
        entire_source_projective_saturation_dimension_bound=22,actual_point_projective_orbit_dimension=22,
        weakest_source_whole_MPS_bridge='Generic curve H+irreducquartic type12 lies in finite closed PGLstable MPS union, hence P0 does. Irreducible source saturation has dimension<=6+13+4-1=22 and contains the actual22dim P0orbit, so equals its orbitclosure and lies in that union.',
        additional_smooth_component_bridge='Full tangent=actual curve dimension23 also proves P0 regular unique-component. The irreducible whole source contains P0, giving the same membership. No formalmapcompleteness is needed.',
        frozen_producer_source_sha256=final['source_sha256'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
