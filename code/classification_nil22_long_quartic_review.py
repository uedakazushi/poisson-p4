#!/usr/bin/env python3
"""Independent complete type7 potential/resonance/seven-source finite review."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

from classification_quadratic_cofactor_review import verify_source,PAIRS,MONS
from classification_minimal_quartic_review import LaurentRing,determinant
from classification_degenerate_quadric_q4_nil31_I_review import rank

ROOT=Path(__file__).resolve().parents[1]
CORE=[2,8,11,12,14,17,18,19]

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

def normalize(R,p):
    return R.scale(p,1/p[max(p)])if p else {}

def unique(R,polys):
    out={}
    for p in polys:
        if p:
            p=normalize(R,p);out[tuple(sorted(p.items()))]=p
    return list(out.values())

def sparse_add(left,right,scale=1):
    out=dict(left)
    for k,v in right.items():
        value=out.get(k,F())+scale*v
        if value:out[k]=value
        else:out.pop(k,None)
    return out

def span_certificates(R,gens,targets):
    """Each target as a checked degree<=2-multiplier combination of gens.

    Quadratic multipliers are required to eliminate n19*n8**2 after the
    linear n19 specialization. This is one finite Macaulay matrix, no GB.
    """
    active=sorted({i for p in gens+targets for ex in p for i,n in enumerate(ex[:20])if n})
    names=['1']+[f'n{i}'for i in active]+[f'n{i}*n{j}'for i,j in it.combinations_with_replacement(active,2)]
    multipliers=list(map(R.parse,names))
    rowdata=[];basis={}
    for gi,g in enumerate(gens):
        for mi,m in enumerate(multipliers):
            row=R.mul(g,m)
            if not row:continue
            index=len(rowdata);rowdata.append((gi,mi,row))
            proof={index:F(1)}
            while row:
                lead=max(row);co=row[lead]
                if lead not in basis:
                    basis[lead]=({k:v/co for k,v in row.items()},{k:v/co for k,v in proof.items()})
                    break
                p,pr=basis[lead];row=sparse_add(row,p,-co);proof=sparse_add(proof,pr,-co)
    certificates=[]
    for target in targets:
        row=dict(target);proof={}
        while row:
            lead=max(row);co=row[lead]
            assert lead in basis,('degreetwo membership failed',R.encode(target))
            p,pr=basis[lead];row=sparse_add(row,p,-co);proof=sparse_add(proof,pr,co)
        check={};entries=[]
        for index,co in sorted(proof.items()):
            gi,mi,p=rowdata[index];check=R.add(check,R.scale(p,co))
            entries.append([gi,names[mi],str(co)])
        assert check==target
        certificates.append(entries)
    return certificates

def ideal_equality(R,left,right):
    left=unique(R,left);right=unique(R,right)
    assert right
    reverse=span_certificates(R,left,right)
    if len(right)==1 and right[0]==R.constant(1):
        forward='Unit ideal generates every polynomial.'
    else:forward=span_certificates(R,right,left)
    return dict(left_generators=[R.encode(p)for p in left],right_generators=[R.encode(p)for p in right],
                right_in_left_degree_two_certificates=reverse,left_in_right_degree_two_certificates=forward,
                both_ideal_inclusions_checked=True)

def matrix_from_pairs(R,entries):
    P=[[{}for _ in range(4)]for _ in range(4)]
    for(i,j),v in zip(PAIRS,entries):P[i][j],P[j][i]=v,R.scale(v,-1)
    return P

def push(R,entries,G,inverse_values):
    P=matrix_from_pairs(R,entries)
    return [sub(R,R.add(*(R.mul(G[i][a],P[a][b],G[j][b])for a,b in it.product(range(4),repeat=2))),inverse_values)
            for i,j in PAIRS]

def parameter_map_from_brackets(R,entries,free):
    coefficients=[]
    for p in entries:
        for mon in MONS:
            coefficients.append({ex[:-4]+(0,)*4:co for ex,co in p.items()if ex[-4:]==mon})
    return {i:coefficients[j]for i,j in enumerate(free)}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    paths=[ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json',
        ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES_20261107.json',
        ROOT/'research/classification/HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES.md']
    full,public=[json.loads(p.read_text())for p in paths[:2]]
    whole=public['nil22_long_quartic_sources'];source=full['full_nilpotent_sources'][6]
    geometry,action,resonance=[whole[key]for key in ['potential_geometry','actual_core_action','rank_resonance']]
    original_review=verify_source(source,7);assert len(original_review['free_indices'])==20
    R=LaurentRing([f'n{i}'for i in range(20)]+['a','r','s','rho'],['x','y','z','w'],
                  invert=['r','rho','n17','n8','n14'])
    parse=R.parse;n=[R.variable(f'n{i}')for i in range(20)]
    x,y,z,w=[R.variable(v)for v in ['x','y','z','w']]
    positions={name:i for i,name in enumerate(R.names)}
    aliases=lambda text:re.sub(r'\b(t|u|v)\b',lambda m:{'t':'x','u':'y','v':'z'}[m[1]],text)
    potparse=lambda text:parse(aliases(text))
    raw=list(map(parse,source['brackets']));J=list(map(parse,source['distinct_nonzero_jacobi_coefficients']))
    pf=parse(source['actual_affine_pfaffian'])
    # (t,u,v,w)=(oldy,oldx-oldy^2/2,oldz-oldw*oldy,oldw).
    Jac=[[parse(v)for v in row]for row in [['0','1','0','0'],['1','-y','0','0'],
                ['0','-w','1','-y'],['0','0','0','1']]]
    assert determinant(R,Jac)==parse('-1')
    inv={positions['x']:parse('y+x**2/2'),positions['y']:x,
         positions['z']:parse('z+w*x'),positions['w']:w}
    straight=push(R,raw,Jac,inv)
    assert straight==list(map(potparse,geometry['all_six_transformed_brackets']))
    phi=potparse(geometry['potential']);V=list(map(potparse,geometry['vector_field']))
    der=lambda p,i:R.derivative(p,i)
    assert straight==[V[0],V[1],V[2],der(phi,3),R.scale(der(phi,2),-1),der(phi,1)]
    div=R.add(*(der(V[i],i+1)for i in range(3)));assert div==parse('-1')
    Q=R.add(phi,*(R.mul(V[i],der(phi,i+1))for i in range(3)))
    assert Q==potparse(geometry['Q0'])
    PP=matrix_from_pairs(R,straight)
    actualJ=[R.add(*(R.mul(PP[a][l],der(PP[b][c],l))
               for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(4)))for i,j,k in it.combinations(range(4),3)]
    assert actualJ==[der(Q,3),R.scale(der(Q,2),-1),der(Q,1),{}]
    assert sub(R,pf,inv)==R.add(phi,R.scale(Q,-1))
    H=[[der(der(phi,i+1),j+1)for j in range(2)]for i in range(2)]
    M=[[der(V[i],j+1)for j in range(2)]for i in range(2)]
    C=[[R.add(H[i][j],*(R.mul(M[k][i],H[k][j])for k in range(2)),
               *(R.mul(H[i][k],M[k][j])for k in range(2)))for j in range(2)]for i in range(2)]
    adj=[[H[1][1],R.scale(H[0][1],-1)],[R.scale(H[1][0],-1),H[0][0]]]
    detH=R.add(R.mul(H[0][0],H[1][1]),R.scale(R.mul(H[0][1],H[1][0]),-1))
    assert R.add(*(R.mul(adj[i][j],C[j][i])for i,j in it.product(range(2),repeat=2)))==R.scale(R.mul(n[14],detH),-2)
    # Full actual affine centralizer and all20 transformed free parameters.
    G=[[parse(v)for v in row]for row in [['1','0','a','0'],['0','1','0','a'],
                                        ['0','0','r','0'],['0','0','0','r']]]
    inverse={positions['x']:parse('x-a*z/r'),positions['y']:parse('y-a*w/r'),
             positions['z']:parse('z/r'),positions['w']:parse('w/r')}
    moved=push(R,raw,G,inverse)
    pmap=parameter_map_from_brackets(R,moved,source['free_coefficient_indices'])
    assert pmap=={int(k[1:]):parse(v)for k,v in action['full_parameter_map'].items()}
    assert moved==[sub(R,p,pmap)for p in raw]and pmap[14]==n[14]
    assert determinant(R,G)==parse('r**2')
    Hnew=[[sub(R,p,pmap)for p in row]for row in H]
    assert Hnew[0][0]==parse('-2*n17*r**2')
    assert Hnew[0][1]==parse('r*(2*a*n17-n18)')
    assert Hnew[1][1]==parse('-3*n11/2+n19+1/2-2*a*a*n17+2*a*n18')
    # The actual w-centering centralizer, not the nonlinear auxiliary coordinates.
    S=[[parse(v)for v in row]for row in [['1','0','0','0'],['0','1','0','0'],
                                       ['0','s','1','0'],['0','0','0','1']]]
    sinv={positions['x']:x,positions['y']:y,positions['z']:parse('z-s*y'),positions['w']:parse('w-s')}
    centered=push(R,raw,S,sinv)
    smap=parameter_map_from_brackets(R,centered,source['free_coefficient_indices'])
    assert smap=={int(k[1:]):parse(re.sub(r'\bshift\b','s',v))for k,v in whole['w_centering_full_parameter_map'].items()}
    assert centered==[sub(R,p,smap)for p in raw]
    assert all(smap[i]==n[i]for i in CORE)and smap[13]==parse('n13-n14*s')
    # Exactly the full source coefficients supported on the8core parameters.
    core=unique(R,[p for p in J if all(not any(ex[i]for i in range(20)if i not in CORE)for ex in p)])
    claimedcore=unique(R,list(map(parse,resonance['all_core_equations'])))
    assert {tuple(sorted(p.items()))for p in core}=={tuple(sorted(p.items()))for p in claimedcore}
    core_certificates=[]
    for claim in resonance['rank1_U_finite_q_checks']:
        qv=F(claim['n14'])
        fixed={17:parse('1'),18:{},14:R.constant(qv),11:R.constant(2*qv+F(1,2)),
               12:parse('-3*n8'),19:R.constant(3*qv+F(1,4))}
        core_certificates.append(dict(name='rank1U_q'+str(qv),
            **ideal_equality(R,[sub(R,p,fixed)for p in core],list(map(parse,claim['core_groebner'])))))
    for claim in resonance['rank1_V_three_q_checks']:
        qv=F(claim['n14']);fixed={17:{},18:{},14:R.constant(qv),11:R.constant(qv+F(1,2))}
        core_certificates.append(dict(name='rank1V_q'+str(qv),
            **ideal_equality(R,[sub(R,p,fixed)for p in core],list(map(parse,claim['core_groebner'])))))
        expected={17:{},18:{},14:R.constant(qv),11:R.constant(qv+F(1,2)),19:R.constant(qv)}
        p=F(-3)if qv==F(1,2)else F(-6)if qv==-1 else F(2)
        t=F(4)if qv==F(1,2)else F(-8)if qv==-1 else F(24)
        expected|={12:R.scale(n[8],p),2:R.scale(R.mul(n[8],n[8]),t)}
        movedcore={i:sub(R,v,expected|{positions['a']:parse(claim['normalizing_a']),positions['r']:parse('1')})for i,v in pmap.items()if i in CORE}
        assert movedcore=={i:expected[i]if i not in [2,8,12]else{}for i in CORE}
    half={17:{},18:{},14:parse('-1/2'),11:{}}
    half_eq=[parse('(3*n19+1)*n8'),parse('4*n12*n8+(4*n19+1)*n2')]
    core_certificates.append(dict(name='rank1V_q_minus_half',
        **ideal_equality(R,[sub(R,p,half)for p in core],half_eq)))
    original_phi=sub(R,phi,{positions['x']:parse('y'),positions['y']:parse('x-y**2/2'),
                          positions['z']:parse('z-w*y'),positions['w']:w})
    top_phi={ex:co for ex,co in original_phi.items()if sum(ex[-4:])==4}
    assert sub(R,top_phi,half)==parse('w**3*(n8*y-n2*w/4)')
    for claim in resonance['rank0_primitive_core_checks']:
        qv=F(claim['n14'])
        fixed={17:{},18:{},19:parse('3*n11/2-1/2'),14:R.constant(qv)}
        eigencondition=parse('(n11-2*n14)*(-n11+n14+1)')
        equations=[sub(R,p,fixed)for p in core]+[sub(R,eigencondition,fixed)]
        core_certificates.append(dict(name='rank0primitive_q'+str(qv),
            **ideal_equality(R,equations,list(map(parse,claim['core_with_primitive_eigencondition_groebner'])))))
    # Every specialized normalizer keeps ALL8core coordinates, including zeros.
    for claim in resonance['all_special_normalizations']:
        name=claim['name']
        if name=='rank1_V_half_yw3':
            fixed=half|{19:parse('-1/3'),2:parse('12*n12*n8')}
        elif name=='rank1_V_half_w4':
            fixed=half|{19:parse('-1/4'),8:{},2:parse('rho**2')}
        else:
            fixed={17:{},18:{},11:{},14:parse('-1'),19:parse('-1/2'),
                   12:n[8],2:parse('4*n8**2')}
        fixed|={positions['a']:parse(claim['a']),positions['r']:parse(claim['r'])}
        assert {i:sub(R,pmap[i],fixed)for i in CORE}=={int(k[1:]):parse(v)for k,v in claim['core'].items()}
    # Full normalized fibers: exact ideal equality, not only substitution of examples.
    families=whole['all_seven_complete_families'];assert len(families)==7
    family_records=[];dimensions=[]
    for index,claim in enumerate(families):
        fixed={int(k[1:]):parse(v)for k,v in claim['core_substitution'].items()}|{13:{}}
        left=[sub(R,p,fixed)for p in J];right=list(map(parse,claim['expected_Jacobi_relations']))
        certificate=ideal_equality(R,left,right)
        rawmap={int(k[1:]):parse(v)for k,v in claim['raw_source_substitution'].items()}
        values={int(k[1:]):parse(v)for k,v in claim['full_source_substitution'].items()}
        assert set(rawmap)==set(values)==set(range(20))
        rawjac=[sub(R,p,rawmap)for p in J]
        rawres=unique(R,[p for p in rawjac if p])
        expectedraw=unique(R,list(map(parse,claim['raw_residual_compatibility'])))
        assert {tuple(sorted(p.items()))for p in rawres}=={tuple(sorted(p.items()))for p in expectedraw}
        assert all(not sub(R,p,values)for p in J)
        brackets=[sub(R,p,values)for p in raw]
        assert brackets==list(map(parse,claim['universal_brackets']))
        actualpf=sub(R,pf,values);assert actualpf==parse(claim['actual_affine_pfaffian'])
        assert sub(R,top_phi,values)==parse(claim['highest_quartic_phi'])
        free={int(s[1:])for s in claim['parameters']}
        assert set(i for p in values.values()for ex in p for i,n in enumerate(ex[:20])if n)==free
        assert all(values[i]==n[i]for i in free)
        irr=claim['irreducibility_certificate'];kind=irr['kind']
        zero_certificate=claim['excluded_zero_complement']
        if index in [0,4,5]:
            coord=2 if index==0 else 0
            coefficient=R.derivative(actualpf,coord)
            assert coefficient==parse(irr['coefficient'].replace('^','**'))and not R.derivative(coefficient,coord)
            assert sub(R,actualpf,{positions['w']:{}})==parse(irr['restriction_at_w_zero'])
            assert any(ex[-4+(0 if index==0 else 2)]==2 for ex in parse(irr['restriction_at_w_zero']))
        elif index in [1,6]:
            coord=0 if index==1 else 2
            assert R.derivative(actualpf,coord)==(n[15]if index==1 else R.scale(n[6],-1))
        else:assert actualpf==parse(irr['exact_form'])and irr['odd_exponent']in [1,3]
        if zero_certificate:
            if index==1:zero={15:{}}
            elif index==2:zero={3:parse('4*n6*n7')}
            elif index==3:zero={5:n[1]}
            else:zero={6:{}}
            # Keep the whole raw rank0 n10 complement before restricting n6!=0.
            zeropf=sub(R,sub(R,pf,rawmap),zero)
            assert zeropf==parse(zero_certificate['exact_pfaffian'])
        point={int(k[1:]):R.constant(F(v))for k,v in claim['canonical_source_values'].items()}
        canonical=[sub(R,p,point)for p in raw]
        weights=claim['canonical_Cstar_weights'];character=claim['canonical_tensor_character']
        for(i,j),p in zip(PAIRS,canonical):
            for ex in p:
                hpower=2-sum(ex[-4:]);assert hpower>=0
                assert hpower*weights[0]+sum(a*b for a,b in zip(ex[-4:],weights[1:]))-weights[i+1]-weights[j+1]==character
        undo=[2,1,1,1,2,2,1][index]
        total=len(free)+undo+1;assert total==7
        dimensions.append(total)
        family_records.append(dict(name=claim['name'],ideal_certificate=certificate,
            all20_substitutions_sixbrackets_Jacobi_Pfaffian_verified=True,
            raw_compatibility_retained=claim['raw_residual_compatibility'],
            zero_unit_complement_checked_before_localization=True,
            primitive_linear_or_odd_power_irreducibility_verified=True,
            canonical_tensor_weights_verified=True,normal_source_dimension=len(free),
            inverse_core_normalizer_parameter_count=undo,w_centering_parameter_count=1,
            fixed_source_dimension_bound=total))
    A=[[F(v)for v in row]for row in source['field_linear_matrix']]
    c=list(map(F,source['field_constant']));cols=[]
    for i,j in it.product(range(4),repeat=2):
        B=[[F(a==i and b==j)for b in range(4)]for a in range(4)]
        cols.append([sum(B[a][k]*A[k][b]-A[a][k]*B[k][b]for k in range(4))for a,b in it.product(range(4),repeat=2)]+[sum(B[a][k]*c[k]for k in range(4))for a in range(4)])
    for j in range(4):cols.append([F()]*16+[-A[a][j]for a in range(4)])
    assert rank(list(map(list,zip(*cols))))[0]==12
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_CAS_import=True,no_producer_import=True,
        original_full20_source=original_review,
        full_polynomial_straightening_and_all_potential_Jacobi_Pfaffian_identities=True,
        full_actual_centralizer_twenty_parameter_map_and_all90_brackets_verified=True,
        full_w_centering_map_and_all_core_complements_verified=True,
        rank_resonance_core_ideal_equalities=core_certificates,
        all_seven_complete_sources=family_records,
        full_fixed_source_dimension_upper_bound=max(dimensions),
        affine_cofactor_orbit_dimension=12,projective_saturated_dimension_upper_bound=22,
        componenthood_or_membership_claimed=False,
        public_migration='V2 loads only the complete public producer packet/public proofnote/original nilpotent table; V1 immutable reviewer bytes and private-input receipt are preserved separately. All four former private structured mathematical sections were also compared equal before migration.',
        mathematical_cover_dependency='Full rank2/rank1/rank0 potential case proof with all zero/localization complements, reviewed separately; computation validates exact normalizers/core/fiber ideals but does not replace that proof.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','full_fixed_source_dimension_upper_bound','projective_saturated_dimension_upper_bound','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
