#!/usr/bin/env python3
"""Independent complete nil22 nonzero-base support and old-component witness.
Standard-library Fraction sparse arithmetic only; no producer/CAS imported.
"""
import hashlib,itertools,json,resource,time,sys
from fractions import Fraction as F
from pathlib import Path
from classification_quadratic_cofactor_review import Ring,PAIRS
from classification_two_lines_plane_cubic_component_review import rank

def plane_even():
    resource.setrlimit(resource.RLIMIT_CPU,(25,26));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1];R=Ring(6)
    source_path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json'
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_EVEN_20261107.json'
    source=next(v for v in json.loads(source_path.read_text())['cases'] if v['name']=='contact4_secant_at_contact');saved=json.loads(path.read_text())
    def parse(p):
        for i in range(6):p=p.replace('m'+str(i),'n'+str(i))
        return R.parse(p.replace('h','1'))
    p=list(map(parse,source['brackets']));n=[R.variable('n'+str(i)) for i in range(6)]
    x,y,z,w=[R.variable(v) for v in ['x','y','z','w']]
    def coordinate_sub(poly,values):
        out={}
        for exp,c in poly.items():
            term={exp[:6]+(0,)*4:c}
            for value,power in zip(values,exp[-4:]):
                for _ in range(power):term=R.mul(term,value)
            out=R.add(out,term)
        return out
    alpha=R.scale(n[2],F(-3,5));beta=R.add(R.scale(n[1],2),R.scale(n[2],F(-4,3)))
    translated=[coordinate_sub(v,[x,y,R.add(z,R.scale(alpha,-1)),R.add(w,R.scale(beta,-1))]) for v in p]
    assert all(sum(e[-4:])!=1 for v in translated for e in v)
    # Independent inverse-coordinate substitution and explicit six-bracket
    # pushforward for (X,Y,Z,W)=(x,9w/10,y,z-9*m0*x/10).
    old_sub=[x,z,R.add(w,R.scale(R.mul(n[0],x),F(9,10))),R.scale(y,F(10,9))]
    pulled=[coordinate_sub(v,old_sub) for v in translated]
    transformed=[R.scale(pulled[2],F(9,10)),pulled[0],pulled[1],R.scale(pulled[4],F(-9,10)),
        R.scale(R.add(R.scale(pulled[5],-1),R.scale(R.mul(n[0],pulled[2]),F(9,10))),F(9,10)),
        R.add(pulled[3],R.scale(R.mul(n[0],pulled[0]),F(9,10)))]
    assert transformed==list(map(parse,saved['full_transformed_six_brackets']))
    J=[R.constant(1),{},{},{},{},R.constant(1)]
    constant=[{e:c for e,c in v.items() if not any(e[-4:])} for v in transformed]
    assert constant==J
    Q=[R.add(v,R.scale(a,-1)) for v,a in zip(transformed,J)]
    assert Q==list(map(parse,saved['quadratic_Q45_six_brackets']))
    assert all(sum(e[-4:])==2 for v in Q for e in v)
    def jac(entries):
        matrix={}
        for (i,j),p in zip(PAIRS,entries):matrix[i,j]=p;matrix[j,i]=R.scale(p,-1)
        return [R.add(*(R.mul(matrix.get((a,t),{}),R.derivative(matrix.get((b,d),{}),t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
    assert not any(jac(transformed)) and not any(jac(Q))
    paths=[source_path,path,root/'code/classification_quadratic_cofactor_review.py',root/'research/classification/BETA_ZERO_COVERAGE.md']
    result={'task_id':'E_reduced_nonnormal','name':'contact4_secant_at_contact','all_passed':True,
        'whole_six_magnetic_source_affine_translation_reconstructed':True,
        'global_Darboux_pushforward_independently_reconstructed':True,'no_parameter_value_excluded':True,
        'constant_J_and_all_linear_coefficient_zeros':True,'all_full_Jacobi_and_entire_Q45_conditions':True,
        'known_support_conclusion':'Every source point is carried to the audited beta=0 support, so the entire six-parameter source lies in its known finite closed union; local uniqueness is unnecessary.',
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'method':'Independent standard-library Fraction sparse arithmetic and explicit inverse affine-coordinate substitutions; no producer/CAS import.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_EVEN_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'whole_source_beta_zero':True,'seconds':result['seconds']}))


def main():
    if '--plane-even' in sys.argv:return plane_even()
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1];R=Ring(27)
    op=root/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    bp=root/'reports/CLASSIFICATION_TWO_QUADRICS_NIL22_BASE_20261107.json'
    rp=root/'reports/CLASSIFICATION_TWO_QUADRICS_NIL22_RANK1_SUPPORT_20261107.json'
    original=json.loads(op.read_text())['full_nilpotent_sources'][5];base=json.loads(bp.read_text())['cases'][1];saved=json.loads(rp.read_text())
    n=[R.variable('n'+str(i)) for i in range(27)];x,y,z,w=[R.variable(v) for v in ['x','y','z','w']]
    def eq(a,b):assert a==b
    def substitute(poly,index,value):
        out={}
        for exp,c in poly.items():
            ee=list(exp);power=ee[index];ee[index]=0;term={tuple(ee):c}
            for _ in range(power):term=R.mul(term,value)
            out=R.add(out,term)
        return out
    def subs(poly,values):
        for i,value in values:poly=substitute(poly,i,value)
        return poly
    def tensor(entries):
        out={}
        for (i,j),v in zip(PAIRS,entries):out[i,j]=v;out[j,i]=R.scale(v,-1)
        return out
    def jac(entries):
        p=tensor(entries)
        return [R.add(*(R.mul(p.get((a,t),{}),R.derivative(p.get((b,d),{}),t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
    def pf(p):return R.add(R.mul(p[0],p[5]),R.scale(R.mul(p[1],p[4]),-1),R.mul(p[2],p[3]))
    def coefficients(polys):
        out=[]
        for poly in polys:
            by_monomial={}
            for exp,c in poly.items():by_monomial.setdefault(exp[-4:],{})[exp[:-4]+(0,)*4]=c
            out+=list(by_monomial.values())
        return out
    def normalize(poly):
        if not poly:return ()
        a=poly[min(poly)];return tuple(sorted((e,c/a) for e,c in poly.items()))
    source=list(map(R.parse,original['brackets']));p=tensor(source)
    div=[R.add(*(R.derivative(p.get((i,j),{}),i) for i in range(4))) for j in range(4)]
    assert div==[y,{},w,{}]
    coeff=coefficients(jac(source))
    assert {normalize(v) for v in coeff}=={normalize(R.parse(v)) for v in original['distinct_nonzero_jacobi_coefficients']}
    disc=R.parse('-4*(3*n18**2+4*n19*n24)/3')
    assert disc in coeff and R.scale(disc,-1) in coeff
    eq(disc,R.parse(saved['global_base_discriminant_Jacobi_coefficient']))
    assert source[4]==R.parse('n18*w*y+n19*w*w-n24*y*y/3')
    fixed=[(18,{}),(19,{}),(24,R.constant(-3))]
    source=[subs(v,fixed) for v in source];coeff=coefficients(jac(source))
    linear=[p for p in coeff if all(sum(e[:27])<=1 for e in p)]
    def row(poly):
        out=[F(0)]*28
        for exp,c in poly.items():
            assert not any(exp[27:]) and sum(exp[:27])<=1
            index=next((i for i,v in enumerate(exp[:27]) if v),27);out[index]+=c
        return out
    expected=list(map(R.parse,['n0','n1','n2','n9+2*n14/5','n10+n16/2-n26/2','n11+n17/4']))
    assert rank([row(v) for v in linear])[0]==6 and rank([row(v) for v in linear+expected])[0]==6 and rank([row(v) for v in expected])[0]==6
    initial=[(0,{}),(1,{}),(2,{}),(9,R.parse('-2*n14/5')),(10,R.parse('(-n16+n26)/2')),(11,R.parse('-n17/4'))]
    source=[subs(v,initial) for v in source]
    assert source==list(map(R.parse,base['full_brackets']))
    coeff=coefficients(jac(source))
    assert {normalize(v) for v in coeff}=={normalize(R.parse(v)) for v in base['remaining_full_Jacobi_coefficients']}
    assert R.parse('-n17**2') in coeff
    after17=[substitute(v,17,{}) for v in coeff]
    assert R.parse('-4*n14**2/5') in after17
    radical=[(17,{}),(14,{})];source=[subs(v,radical) for v in source]
    coeff=[subs(v,radical) for v in coeff if subs(v,radical)]
    for stage in saved['global_polynomial_eliminations']:
        index=int(stage['parameter'][1:]);polynomial=R.parse(stage['actual_Jacobi_equation']);solution=R.parse(stage['exact_polynomial_solution']);a=F(stage['nonzero_constant_coefficient'])
        assert a and normalize(polynomial) in {normalize(v) for v in coeff}
        eq(polynomial,R.scale(R.add(n[index],R.scale(solution,-1)),a))
        source=[substitute(v,index,solution) for v in source]
        coeff=[substitute(v,index,solution) for v in coeff if substitute(v,index,solution)]
    assert source==list(map(R.parse,saved['full_reduced_source_before_branching']))
    expected=list(map(R.parse,['n12*(3*n16-n26-2)','n12*(3*n13-n22)']))
    assert {normalize(v) for v in coeff}=={normalize(v) for v in expected}
    branches=[]
    for branch in saved['branches']:
        fixes=[(int(key[1:]),R.parse(value)) for key,value in branch['fixed_parameters'].items()]
        entries=[subs(v,fixes) for v in source]
        assert entries==list(map(R.parse,branch['full_brackets'])) and not any(jac(entries))
        actual=pf(entries);eq(actual,R.parse(branch['actual_Pfaffian']))
        free={i for poly in entries for exp in poly for i,value in enumerate(exp[:27]) if value}
        assert free=={int(v[1:]) for v in branch['free_parameters']}
        if branch['name']=='closed_n12_zero':
            eq(actual,R.mul(y,y,R.parse(branch['Pfaffian_divided_by_y_squared'])))
        else:
            eq(R.derivative(actual,0),R.scale(R.mul(y,y,w),-1));eq(R.derivative(actual,2),R.mul(y,y,y))
            eq(substitute(actual,28,{}),R.scale(R.mul(n[12],n[12]),-1))
        branches.append({'name':branch['name'],'full_source_free_dimension':len(free),'all_full_Jacobi_and_Pfaffian_identities':True})
    # Reconstruct Q_* in ITS ORIGINAL coordinate order, then perform the
    # independent permutation/sign pushforward, not a copied transformed row.
    old_entries=list(map(R.parse,['1-x*x+x*z-x*w','-x*x','0','3*x*y-4*x*z+z*z-3*z*w',
        'x*x+x*y-2*x*z+2*x*w+z*z/2-2*z*w+w*w','1+x*x+x*w']))
    old=tensor(old_entries);old_sub=[y,R.scale(z,-1),R.scale(w,-1),x]
    def coordinate_sub(poly):
        out={}
        for exp,c in poly.items():
            assert not any(exp[:27]);term=R.constant(c)
            for value,power in zip(old_sub,exp[-4:]):
                for _ in range(power):term=R.mul(term,value)
            out=R.add(out,term)
        return out
    permutation=[3,0,1,2];sign=[1,1,-1,-1]
    pushed=[R.scale(coordinate_sub(old.get((permutation[i],permutation[j]),{})),sign[i]*sign[j]) for i,j in PAIRS]
    recognition=saved['known_quartic_component_recognition'];witness=[(int(key[1:]),R.constant(F(v))) for key,v in recognition['branch_parameter_values'].items()]
    actual=[subs(R.parse(v),witness) for v in saved['branches'][1]['full_brackets']]
    assert actual==pushed==list(map(R.parse,recognition['all_six_transformed_old_brackets']))
    assert dict(witness)[12]==R.constant(-1)
    paths=[op,bp,rp,root/'legacy/v06/nilpotent22.tex',root/'code/classification_quadratic_cofactor_review.py',root/'code/classification_two_lines_plane_cubic_component_review.py']
    result={'task_id':'E_reduced_nonnormal','all_passed':True,
        'global_discriminant_is_actual_Jacobi_coefficient':True,'complete_six_linear_constraint_rowspace_reconstructed':True,
        'actual_square_equations_and_all_six_unit_polynomial_eliminations':True,
        'remaining_product_equation_ideal_exact':True,'both_closed_branches_reconstructed':branches,
        'rank1_nonzero_Pfaffian_is_irreducible_by_primitive_linear_Gauss':True,
        'old_Qstar_original_six_brackets_independently_pushed_forward':True,'exact_old_witness_in_eight_parameter_branch':True,
        'whole_component_membership_input':'The audited old unique smooth23-dimensional H_QUARTIC_ESSENTIAL5 theorem; this arithmetic checker reconstructs the actual source inclusion, not that already audited ambient theorem.',
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'Complete reduced support of nil22 image with nonzero binary base bracket; all inverse/closed complements retained. Hence every reduced H+two-irreducible-quadrics point of this modular type is plane; irreducible-quartic open branch lies in the audited old H_QUARTIC component by irreducibility and its unique-component witness.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_NIL22_RANK1_REVIEW_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'branch_dimensions':[v['full_source_free_dimension'] for v in branches],'seconds':result['seconds']}))


if __name__=='__main__':main()
