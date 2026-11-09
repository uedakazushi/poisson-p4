#!/usr/bin/env python3
"""Exact full-magnetic normal forms for the remaining nil22 image sources.

No point sampling or jets. Fiber affine translations are solved as a complete
six-by-six action, with inverse maps and all tensor coefficients retained.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

import sympy as s
from classification_nil22_image_source_review import build_P,coords,p,F0,mag,x,y,q1,q2
import classification_translated_minimal_quartic_review as exact

ROOT=Path(__file__).resolve().parents[1]


def jacobi(P):
    return [s.expand(sum(P[i,l]*s.diff(P[j,k],coords[l])+P[j,l]*s.diff(P[k,i],coords[l])+P[k,l]*s.diff(P[i,j],coords[l])for l in range(4)))for i,j,k in it.combinations(range(4),3)]


def full_vector(P):
    upper={}
    for i,j in it.combinations(range(4),2):
        polynomial={(2-sum(e),)+e+(0,0):Fraction(c)for e,c in s.Poly(P[i,j],*coords).terms()if c}
        assert all(sum(e)<=2 for e,c in s.Poly(P[i,j],*coords).terms())
        if polynomial:upper[i+1,j+1]=polynomial
    return exact.arithmetic.antisymmetric(upper)


def twists(output,koszul=False):
    begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    P,_,_=build_P(x*x*y,s.Matrix([[s.Rational(8,3),1],[1,-s.Rational(3,64)]]),s.Matrix([-s.Rational(1,8),-s.Rational(5,3)]),s.diag(-s.Rational(1,4),-s.Rational(1,2)))
    P=P.subs(dict(zip(mag,[0,1,0,0,0,0])))
    N=s.Matrix([x,y,0,0]);z=s.symbols('z0:'+('25'if koszul else'20'));tau=s.Symbol('tau')
    if koszul:
        mons=[s.Integer(1),*coords,*[coords[i]*coords[j]for i,j in it.combinations_with_replacement(range(4),2)]]
        W=s.zeros(4);W[0,1]=sum(a*b for a,b in zip(z[:15],mons))
        for j in range(2):
            row=sum(z[15+5*j+k]*f for k,f in enumerate([1,*coords]))
            W[0,j+2]=x*row;W[1,j+2]=y*row
        for i,j in it.combinations(range(4),2):W[j,i]=-W[i,j]
        for i,j,k in it.combinations(range(4),3):assert s.expand(N[i]*W[j,k]-N[j]*W[i,k]+N[k]*W[i,j])==0
        Y=None
    else:
        Y=s.Matrix([sum(z[5*i+j]*f for j,f in enumerate([1,*coords]))for i in range(4)])
        W=N*Y.T-Y*N.T
    linear=[s.Poly(f,tau).coeff_monomial(tau)for f in jacobi(P+tau*W)]
    equations=[c for f in linear for e,c in s.Poly(f,*coords).terms()if c]
    A,rhs=s.linear_eq_to_matrix(equations,z);assert rhs==s.zeros(A.rows,1)
    K=A.nullspace();fields=[None]*len(K)if koszul else[Y.subs(dict(zip(z,v)))for v in K]
    bivectors=[W.subs(dict(zip(z,v))).applyfunc(s.expand)for v in K]
    D,G,rows=exact.matrices(full_vector(P))
    assert exact.arithmetic.rank(D)[0]==100 and exact.arithmetic.rank(G)[0]==48
    vectors=[exact.vector(full_vector(W))for W in bivectors]
    assert all(sum(value*v[j]for j,value in row)==0 for row in rows for v in vectors)
    joint=[row+[v[i]for v in vectors]for i,row in enumerate(G)]
    transverse=exact.arithmetic.rank(joint)[0]-48
    eta=s.symbols('eta0:'+str(len(K)))
    total=sum([v*coefficient for v,coefficient in zip(bivectors,eta)],s.zeros(4))
    self_equations=sorted({str(s.factor(c))for f in jacobi(total)for e,c in s.Poly(f,*coords).terms()if c})
    candidates=[]
    for index,(field,W)in enumerate(zip(fields,bivectors)):
        if W==s.zeros(4):continue
        magnetic_correction=None
        if any(jacobi(W)):
            if not koszul:continue
            mu=s.symbols('mu0:6');magnetic=s.zeros(4)
            magnetic[0,1]=sum(a*b for a,b in zip(mu,[1,x,y,x*x,x*y,y*y]));magnetic[1,0]=-magnetic[0,1]
            eq=[c for f in jacobi(W+magnetic)for e,c in s.Poly(f,*coords).terms()if c]
            solution_set=s.linsolve(eq,mu)
            if solution_set is s.EmptySet or not solution_set:continue
            solution=list(solution_set)[0];free=set().union(*(value.free_symbols for value in solution))
            values=[value.subs(dict.fromkeys(free,0))for value in solution]
            magnetic_correction=list(map(str,values))
            W=(W+magnetic.subs(dict(zip(mu,values)))).applyfunc(s.expand)
            assert not any(jacobi(W))
        family=P+tau*W;assert not any(jacobi(family))
        Z=s.Matrix([s.expand(sum(s.diff(family[i,j],coords[i])for i in range(4)))for j in range(4)])
        B=Z.jacobian(s.Matrix(coords));c=Z.subs(dict.fromkeys(coords,0))
        assert (Z-B*s.Matrix(coords)-c).applyfunc(s.expand)==s.zeros(4,1)
        scalar=s.factor(s.trace(B*B)/4)
        equal_pair=scalar!=0 and (B*B-scalar*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
        determinant=s.factor(B.det());pf_at_zero=None;zero=None;accepted=False
        if equal_pair and determinant!=0:
            zero=(-B.inv()*c).applyfunc(s.cancel)
            C=family.subs(dict(zip(coords,zero)),simultaneous=True).applyfunc(s.cancel)
            pf_at_zero=s.factor(C[0,1]*C[2,3]-C[0,2]*C[1,3]+C[0,3]*C[1,2])
            accepted=pf_at_zero!=0
        candidates.append(dict(compatibility_basis_index=index,linear_field=list(map(str,field))if field is not None else None,magnetic_selfPoisson_correction=magnetic_correction,
            entire_pencil_strictJacobi_zero=True,modular=list(map(str,Z)),cofactor_characteristic_polynomial=str(s.factor(B.charpoly(s.Symbol('lambda')).as_expr())),
            cofactor_square_is_nonzero_scalar=bool(equal_pair),cofactor_square_scalar=str(scalar),
            modular_zero=list(map(str,zero))if zero is not None else None,
            constant_Pf_at_modular_zero=str(pf_at_zero),generic_SS20_rank2_hypotheses_verified=bool(accepted),
            whole_family_brackets={f'{i},{j}':str(s.factor(family[i,j]))for i,j in it.combinations(range(4),2)if family[i,j]!=0}))
    accepted=[p for p in candidates if p['generic_SS20_rank2_hypotheses_verified']]
    paths=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'code/classification_translated_minimal_quartic_review.py',ROOT/'code/classification_nil22_weighted_slice_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_MAGNETIC_NORMALFORMS_20261107.json',ROOT/'research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,whole_linear_ansatz_dimension=25 if koszul else 20,
        ansatz='All quadraticW with NwedgeWzero: Wqq arbitraryquadratic, cross=p tensor arbitraryaffinerow, Wppzero'if koszul else'All affineY via NwedgeY',
        compatibility_matrix_rank=A.rank(),entire_compatibility_kernel_dimension=len(K),compatibility_kernel_basis=[list(map(str,v))for v in K],
        genuine_normal_image_dimension=transverse,all_full350_derivative_products_zero=True,
        full_selfPoisson_quadratic_conditions=self_equations,exact_actual_basis_pencils=candidates,
        SS20_actual_pencil_count=len(accepted),whole_double_A_source_MPS_covered=bool(accepted),
        coverage_logic='If an acceptedpencil is present, generic tau gives semisimplecofactor with two opposite multiplicity2 eigenvalues and nondegenerate constant value at its modularzero. The alreadyproved wholeSS20rank2 theorem puts that densepuncturedactualfamily in the finiteclosedMPSunion; closedness puts tau0 canonicalpoint there, and explicitwhole6source orbitnormalization gives the entire double_A source.'if accepted else'No membership conclusion from this finitekernel/basispencil test; it does not exclude all other actualfamilies or highercontact curves.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},reviewer_source_bytes=Path(__file__).read_text(),
        seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k]for k in ['all_passed','entire_compatibility_kernel_dimension','genuine_normal_image_dimension','SS20_actual_pencil_count','seconds']}))


def graphs(output):
    from s61_model import S61,X,pushforward
    begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    model=S61();I0,I1,I2=s.symbols('I0 I1 I2');records=[]
    specifications=[(0,s.Matrix([[s.Rational(8,3),1],[1,-s.Rational(3,64)]]),s.Matrix([-s.Rational(1,8),-s.Rational(5,3)]),s.diag(-s.Rational(1,4),-s.Rational(1,2)),x,[],(0,2)),
        (1,s.diag(-8,1),s.Matrix([0,1]),s.diag(-s.Rational(1,4),-s.Rational(1,2)),I0+I1*x,[I0,I1],(0,2)),
        (2,s.Matrix([[2,1],[1,0]]),s.Matrix([0,-2]),s.diag(0,-1),I0+I1*x+I2*x*x,[I0,I1,I2],(0,3))]
    for index,S,V0,V1,residual,parameters,pivot in specifications:
        affine,_,_=build_P(x*x*y,S,V0,V1)
        affine=affine.subs(dict.fromkeys(mag,0));affine[0,1]+=residual;affine[1,0]=-affine[0,1]
        P=s.zeros(5)
        for i,j in it.combinations(range(4),2):
            entry=sum(co*X[0]**(2-sum(e))*s.prod(z**a for z,a in zip(X[1:],e))for e,co in s.Poly(affine[i,j],*coords).terms())
            P[i+1,j+1],P[j+1,i+1]=entry,-entry
        K=P.subs(dict(zip(X,[1,0,0,0,0])))[1:,1:]
        pf=s.factor(K[0,1]*K[2,3]-K[0,2]*K[1,3]+K[0,3]*K[1,2]);assert not pf.free_symbols and pf!=0
        i,j=pivot;unit=[s.eye(4)[:,k]for k in range(4)];r0=unit[i];r1=unit[j]/K[i,j]
        other=[]
        for k in range(4):
            if k in pivot:continue
            r=unit[k];other.append(r-(r.T*K*r1)[0]*r0+(r.T*K*r0)[0]*r1)
        r2,r3=other;r3=r3/s.cancel((r2.T*K*r3)[0])
        frame=s.eye(5);frame[1:,1:]=s.Matrix.vstack(*(r.T for r in [r0,r1,r2,r3]))
        J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1
        assert (frame[1:,1:]*K*frame[1:,1:].T-J).applyfunc(s.cancel)==s.zeros(4)
        assert all(not s.denom(s.cancel(value)).free_symbols for value in frame)
        framed=pushforward(P,frame);normalized=model.normalize(framed)
        verification=model.verify_point(normalized,check_equations=True)
        point=normalized['parameters'];us=s.symbols('u1:62')
        assert all(not s.denom(s.cancel(value)).free_symbols for value in [*point,*normalized['perspective']])
        results={}
        for name,path in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
            packet=json.loads((ROOT/path).read_text());indices=packet['weight2_u_indices_zero_based']
            residuals=[s.factor(s.cancel(point[k]-s.sympify(expression).subs(dict(zip(us,point)),simultaneous=True)))for k,expression in zip(indices,packet['quadratic_section_q'])]
            results[name]=dict(all45_graph_residuals=list(map(str,residuals)),entire_parameter_source_in_component_graph=all(value==0 for value in residuals),
                nonzero_residuals={str(k):str(value)for k,value in zip(indices,residuals)if value!=0})
        records.append(dict(source_case_index=index,residual_source_parameters=list(map(str,parameters)),canonical_F0=str(residual),
            rank4_evaluation_Pf=str(pf),no_residual_parameter_localization=True,
            Darboux_frame=[[str(a)for a in row]for row in frame.tolist()],perspective=list(map(str,normalized['perspective'])),
            all61_parameters=list(map(str,point)),faithful150_reconstruction_and_full134_equations=verification,
            exact_graph_membership=results,
            whole_original_six_source_membership_in_known_union=any(value['entire_parameter_source_in_component_graph']for value in results.values()),
            scope='For A, explicitactualorbit normalization plus its boundarycurve transfers anypositivecanonicalpoint result to the entire6source. For B/C all2or3 magneticinvariants are retained with a globallyinvertible fixedpivotDarboux frame; positiveidentities cover allparameters. Negativeidentities exclude only whole-source containment in these two graphs, not otherMPScomponents.'))
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_MAGNETIC_NORMALFORMS_20261107.json',ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json',*[ROOT/('legacy/v02/data/'+name)for name in ['symplectic_slice_certificate.json','symplectic_slice_equations.json','symplectic_slice_reduced_equations.json']]]
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,cases=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        no_arbitrary_point_sampling=True,no_higher_jets=True,reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(all_passed=True,membership=[{k:v['entire_parameter_source_in_component_graph']for k,v in r['exact_graph_membership'].items()}for r in records],seconds=result['seconds'])))


def normal_form(index,S,V0,V1):
    P,M,D=build_P(x*x*y,S,V0,V1)
    monomials=[s.Integer(1),x,y,x*x,x*y,y*y]
    alpha=s.diff(P[0,1],q1);beta=s.diff(P[0,1],q2)
    def deriv(field,value):return sum(field[j]*s.diff(value,p[j])for j in range(2))
    effects=[s.expand(-alpha*z-deriv(M[1,:],z))for z in [1,x,y]]+[s.expand(-beta*z+deriv(M[0,:],z))for z in [1,x,y]]
    T=s.Matrix([[s.Poly(f,x,y).coeff_monomial(mon)for f in effects]for mon in monomials])
    invariants=s.Matrix.vstack(*[a.T for a in T.T.nullspace()])
    assert invariants*T==s.zeros(invariants.rows,6)
    pivots=list(invariants.rref()[1]);assert len(pivots)==invariants.rows
    R=invariants[:,pivots].inv()*invariants
    canonical=s.zeros(6,1)
    for j,position in enumerate(pivots):canonical[position]=(R*s.Matrix(mag))[j]
    solution,free=T.gauss_jordan_solve(canonical-s.Matrix(mag))
    solution=solution.subs({a:0 for a in free})
    assert (T*solution-canonical+s.Matrix(mag)).applyfunc(s.expand)==s.zeros(6,1)
    L1=sum(solution[j]*mon for j,mon in enumerate([1,x,y]))
    L2=sum(solution[j+3]*mon for j,mon in enumerate([1,x,y]))
    newcoordinates=s.Matrix([q1+L1,q2+L2,x,y]);G=newcoordinates.jacobian(s.Matrix(coords))
    assert G.det()==1
    transported=(G*P*G.T).subs({q1:q1-L1,q2:q2-L2},simultaneous=True).applyfunc(s.expand)
    canonical_F=sum(canonical[j]*mon for j,mon in enumerate(monomials))
    expected=P.copy();expected[0,1]+=canonical_F-F0;expected[1,0]=-expected[0,1]
    assert (transported-expected).applyfunc(s.expand)==s.zeros(4)
    record=dict(source_case_index=index,full_magnetic_action_matrix=[[str(a)for a in row]for row in T.tolist()],
        translation_rank=T.rank(),complete_translation_kernel=[[str(a)for a in col]for col in T.nullspace()],
        full_magnetic_invariant_map=[[str(a)for a in row]for row in R.tolist()],
        residual_monomials=list(map(str,[monomials[j]for j in pivots])),residual_invariants=list(map(str,R*s.Matrix(mag))),
        exact_fiber_translation=[str(s.factor(L1)),str(s.factor(L2))],
        all_six_brackets_and_all90_coefficients_conjugacy_checked=True,translation_determinant='1',
        no_magnetic_coefficient_inverted=True,canonical_residual_F0=str(s.factor(canonical_F)),
        actual_affine_determinant=str(s.factor(D)),ambient_component_membership_claimed=False)
    if index==0:
        assert R.rows==1 and pivots==[1]
        invariant=s.expand((R*s.Matrix(mag))[0])
        assert invariant==mag[1]-mag[3]-4*mag[4]/3
        I=s.Symbol('I',nonzero=True);scale=1/I
        canonical_P=P.subs(dict(zip(mag,[0,I,0,0,0,0])))
        scaling=s.diag(scale,scale,1,1)
        actual=(scaling*canonical_P*scaling.T).subs({q1:q1/scale,q2:q2/scale},simultaneous=True)/scale
        target=P.subs(dict(zip(mag,[0,1,0,0,0,0])))
        assert (actual-target).applyfunc(s.cancel)==s.zeros(4)
        tau=s.Symbol('tau');curve=P.subs(dict(zip(mag,[0,tau,0,0,0,0])))
        for i,j,k in it.combinations(range(4),3):
            assert s.expand(sum(curve[i,l]*s.diff(curve[j,k],coords[l])+curve[j,l]*s.diff(curve[k,i],coords[l])+curve[k,l]*s.diff(curve[i,j],coords[l])for l in range(4)))==0
        record.update(entire_nonzero_invariant_source_in_one_actual_projective_orbit=True,
            invariant_nonzero_normalizer='q_new=q/I, tensor_scale=I',
            zero_invariant_complement='All zero-invariant points translate exactly to F0=0; P_tau with F0=tau*x is an actual polynomial Poisson curve from the nonzero orbit to it.',
            entire_six_parameter_source_in_single_actual_orbit_closure=True,
            all_points_in_single_actual_orbit_claimed=False,
            generic_canonical_F0='x',boundary_canonical_F0='0',
            generic_point_torus_from_q_scaling_available=False,
            boundary_point_q_scaling_torus_available=True)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--twists',action='store_true',help='Complete affine N-wedge-Y compatibility and exact actual pencils at canonical double_A')
    parser.add_argument('--koszul-pencils',action='store_true',help='All quadratic N-wedge-Wzero corrections, including six Koszul cohomology classes')
    parser.add_argument('--graphs',action='store_true',help='Faithful S61/all134/full45 T39/T40 graph identities with all residual magneticinvariants')
    args=parser.parse_args()
    if args.graphs:return graphs(args.output)
    if args.koszul_pencils:return twists(args.output,True)
    if args.twists:return twists(args.output)
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    specifications=[(0,s.Matrix([[s.Rational(8,3),1],[1,-s.Rational(3,64)]]),s.Matrix([-s.Rational(1,8),-s.Rational(5,3)]),s.diag(-s.Rational(1,4),-s.Rational(1,2))),
        (1,s.diag(-8,1),s.Matrix([0,1]),s.diag(-s.Rational(1,4),-s.Rational(1,2))),
        (2,s.Matrix([[2,1],[1,0]]),s.Matrix([0,-2]),s.diag(0,-1))]
    records=[normal_form(*entry)for entry in specifications]
    inputs=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_import=True,
        input_sha256={str(z.relative_to(ROOT)):hashlib.sha256(z.read_bytes()).hexdigest()for z in inputs},cases=records,
        implication='double_A entire6source is reduced by explicit affine/projective maps to its F0=x actualorbit and the F0=0 boundary. This reduces recognition to one actual orbit, without proving any ambient component membership.',
        reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(all_passed=True,translation_ranks=[a['translation_rank']for a in records],seconds=result['seconds'])))


if __name__=='__main__':main()
