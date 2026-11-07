#!/usr/bin/env python3
"""Independent complete type10 source bridge and three actual orbit ranks.

This excludes type10 as the GENERIC point of a dimension25--27 component.
It does not assert that its primitive-beta image is non-dominant.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

from classification_quadratic_cofactor_review import Ring, MONS, PAIRS
import classification_degenerate_quadric_q4_nil31_I_review as q


def substitute_constant(R,p,index,value):
    out = {}
    for e,c in p.items():
        c *= value**e[index]
        v = list(e);v[index]=0
        out[tuple(v)] = out.get(tuple(v),F(0))+c
    return {e:c for e,c in out.items() if c}


def source_bridge(source,certificate):
    R = Ring(25)
    n = [R.variable(f'n{i}') for i in range(25)]
    x,y,z,w = [R.variable(s) for s in ('x','y','z','w')]
    one = R.constant(1)
    r = R.add(x,R.scale(R.mul(y,y),F(-1,2)))
    A = R.add(n[1],R.mul(n[2],z),R.mul(n[3],w),R.mul(n[4],z,z),R.mul(n[5],z,w),R.mul(n[6],w,w))
    d = R.add(n[7],R.mul(n[8],z),R.mul(n[9],w))
    e = R.add(n[10],R.mul(n[11],z),R.mul(R.add(n[0],R.scale(one,-1),R.scale(n[8],-1)),w))
    L = R.add(n[15],R.mul(n[13],w),R.scale(R.mul(n[14],z),-1))
    G = R.add(R.mul(n[16],z),R.mul(n[17],w),R.mul(n[18],z,z),R.mul(n[19],z,w),R.mul(n[20],w,w),
              R.mul(n[21],z,z,z),R.mul(n[22],z,z,w),R.mul(n[23],z,w,w),R.mul(n[24],w,w,w))
    def data(v,original):
        phi = R.add(R.scale(R.mul(n[12],v,v),F(1,2)),R.mul(v,L),G)
        a = R.add(R.mul(n[0],v),A)
        entries = [a,R.add(R.mul(n[13],v),R.derivative(G,3)),
                   R.add(R.mul(n[14],v),R.scale(R.derivative(G,2),-1)),d,e,R.add(R.mul(n[12],v),L)]
        if original:
            entries[1] = R.add(entries[1],R.mul(y,d))
            entries[2] = R.add(entries[2],R.mul(y,e))
        return entries,phi,a
    entries,phi,a = data(r,True)
    assert all(sum(ex[-4:])<=2 for p in entries for ex in p)
    P = [[{} for j in range(4)] for i in range(4)]
    for (i,j),p in zip(PAIRS,entries):P[i][j],P[j][i]=p,R.scale(p,-1)
    div = [R.add(*(R.derivative(P[i][j],i) for i in range(4))) for j in range(4)]
    assert div==[y,one,{},{}]
    for i,j in PAIRS:
        lie = R.add(R.mul(y,R.derivative(P[i][j],0)),R.derivative(P[i][j],1),
                    R.scale(P[1][j],-1) if i==0 else {},R.scale(P[i][1],-1) if j==0 else {})
        assert not lie
    vectors = []
    for p in entries:
        for mon in MONS:
            vector = [F(0)]*26
            for ex,c in p.items():
                if ex[-4:]==mon:
                    assert sum(ex[:-4])<=1
                    idx = 0 if not any(ex[:-4]) else ex[:-4].index(1)+1
                    vector[idx] += c
            vectors.append(vector)
    free = source['free_coefficient_indices']
    assert source['free_parameter_count']==len(free)==25
    matrix = [vectors[j][1:] for j in free]
    offsets = [vectors[j][0] for j in free]
    inverse = [[F(v) for v in row] for row in certificate['inverse_parameter_matrix']]
    assert q.matmul(matrix,inverse)==q.matmul(inverse,matrix)==[[F(i==j) for j in range(25)] for i in range(25)]
    assert q.determinant(matrix)==F(certificate['affine_parameter_matrix_determinant'])
    for actual,raw in zip(vectors,source['all_ninety_coefficients']):
        coeff = [sum(actual[k+1]*inverse[k][j] for k in range(25)) for j in range(25)]
        constant = actual[0]-sum(v*o for v,o in zip(coeff,offsets))
        expected = R.constant(constant)
        for j,v in enumerate(coeff):expected=R.add(expected,R.scale(n[j],v))
        assert expected==R.parse(raw)
    pf = R.add(R.mul(entries[0],entries[5]),R.scale(R.mul(entries[1],entries[4]),-1),R.mul(entries[2],entries[3]))
    Vphi = R.add(R.scale(R.mul(a,R.derivative(phi,0)),-1),R.mul(d,R.derivative(phi,2)),R.mul(e,R.derivative(phi,3)))
    assert pf==R.scale(Vphi,-1)
    pf_delta0 = substitute_constant(R,pf,12,F(0))
    assert not any(sum(ex[-4:])==4 for ex in pf_delta0)
    # In auxiliary (r,s,z,w), use x=r,y=s only as an algebraic calculation.
    pp,ph,aa = data(x,False)
    PP = [[{} for j in range(4)] for i in range(4)]
    for (i,j),p in zip(PAIRS,pp):PP[i][j],PP[j][i]=p,R.scale(p,-1)
    Vph = R.add(R.scale(R.mul(aa,R.derivative(ph,0)),-1),R.mul(d,R.derivative(ph,2)),R.mul(e,R.derivative(ph,3)))
    psi = R.add(Vph,ph)
    jac = [R.add(*(R.mul(PP[a][l],R.derivative(PP[b][c],l)) for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(4))) for i,j,k in it.combinations(range(4),3)]
    assert jac==[R.scale(R.derivative(psi,3),-1),R.derivative(psi,2),{},R.derivative(psi,0)]
    half = lambda p:substitute_constant(R,p,0,F(1,2))
    dd,ee,ps = half(d),half(e),half(psi)
    baseL = R.add(R.mul(dd,R.derivative(L,2)),R.mul(ee,R.derivative(L,3)))
    constraint = R.add(R.mul(n[12],A),R.scale(L,F(-1,2)),R.scale(baseL,-1))
    Y = R.add(R.mul(n[12],x),L)
    Htilde = R.add(R.scale(R.mul(n[12],G),2),R.scale(R.mul(L,L),-1))
    baseH = R.add(R.mul(dd,R.derivative(Htilde,2)),R.mul(ee,R.derivative(Htilde,3)))
    assert R.scale(R.mul(n[12],ps),2)==R.add(baseH,Htilde,R.scale(R.mul(Y,constraint),-2))
    assert R.add(R.scale(R.mul(n[12],half(ph)),2),R.scale(R.mul(Y,Y),-1))==Htilde
    return dict(all25_parameter_affine_inverse_checked=True,all90_original_coefficients_matched=True,
                all_original_degree_constraints_and_divergence_and_Lie_invariance=True,
                raw_pfaffian_minus_Vphi_checked=True,delta0_highestPf0=True,
                auxiliary_coordinate_Jacobi_equals_gradient_psi=True,
                complete_square_polynomial_identity_with_delta_denominator_cleared=True,
                actual_affine_change='x_new=x+L(z,w)/delta, with delta!=0; auxiliary x-y^2/2 is never used as PGL')


def weight_table():
    found = set()
    for i in range(4):
        for j in range(4-i):
            if 1<=i+j<=3 and i!=j:
                lam = (F(j,2)-1)/F(i-j)
                mu = F(-1,2)-lam
                found.add(tuple(sorted((lam,mu))))
    expected = {(F(-1),F(1,2)),(F(-1,3),F(-1,6)),(F(-1,2),F(0))}
    assert found==expected
    monomials = {}
    for lam,mu in sorted(found):
        allowed = [(i,j) for i in range(4) for j in range(4-i) if 1<=i+j<=3 and i*lam+j*mu==-1]
        monomials[f'{lam},{mu}'] = allowed
    assert monomials['-1,1/2']==[(1,0)]
    assert monomials['-1/3,-1/6']==[(3,0)]
    assert monomials['-1/2,0']==[(2,0),(2,1)]
    assert all(F(-d,4)!=-1 for d in range(4))
    return dict(all_trace_minus_half_spectra_and_all_degree_le3_monomials=monomials,
                repeated_eigenvalue_no_weight_minus1=True,
                zero_eigen_translation='H=z^2(a+bw), kappa*d_w(H)=kappa*b*z^2=0. Irreducibility forces b!=0, so kappa=0; affine w translation removes a.',
                all_normalizing_scales_invertible=True,
                cubic_root_choice='For H=b*z^3 choose t^3=b/delta !=0,u=1/(delta*t); finite algebraic root choice is allowed over C.')


def orbit_ranks(normal_forms,ambient):
    R = Ring(0)
    result = []
    for model in normal_forms['three_actual_Poisson_models']:
        saved = next(v for v in ambient['results'] if v['name']==model['name'])
        scale = F(saved['integer_scale'])
        P = {}
        for (i,j),text in zip(it.combinations(range(1,5),2),model['brackets']):
            P[i,j] = {(2-sum(ex),)+ex+(0,0,0):c*scale for ex,c in R.parse(text).items()}
        P = q.antisymmetric(P)
        def parse(text):
            return q.parse(re.sub(r'\b([xyzw])\b',lambda m:dict(x='x1',y='x2',z='x3',w='x4')[m[1]],text))
        for ij in q.PAIRS:
            assert P.get(ij,{})==parse(saved['full_homogeneous_witness'].get(f'{ij[0]},{ij[1]}','0'))
        assert all(not j for j in q.jacobi(P).values())
        Vcols,Ocols = [],[]
        for a,b in it.product(range(5),repeat=2):
            V,O = {},{}
            for i,j in q.PAIRS:
                V[i,j] = q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
                O[i,j] = q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
            Vcols.append(q.vec(V));Ocols.append(q.vec(O))
        Gcols = Vcols+Ocols+[q.vec(P)]
        assert q.rank(q.transpose(Vcols))[0]==24
        assert q.rank(q.transpose(Vcols+[q.vec(P)]))[0]==25
        assert q.rank(q.transpose(Gcols))[0]==47
        # Use the saved row/column choice as a separately evaluated lower-rank
        # certificate; sign conventions for actual GL columns may differ.
        cert = saved['certificates']['group']
        matrix = q.transpose(Gcols)
        determinant = q.determinant([[matrix[i][j] for j in cert['columns']] for i in cert['rows']])
        assert determinant
        result.append(dict(name=model['name'],vertical_rank=24,vertical_plus_scalar_rank=25,
                           full_GL5_vertical_scalar_rank=47,projective_PGL_orbit_dimension=22,
                           group_minor_rows=cert['rows'],group_minor_columns=cert['columns'],group_minor_determinant=str(determinant)))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(90,95));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    start = time.monotonic()
    paths = ['reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json',
             'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS_20261107.json',
             'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_AMBIENT_20261107.json']
    source,forms,ambient = [json.loads(Path(p).read_text()) for p in paths]
    bridge = source_bridge(source['full_nilpotent_sources'][9],forms['complete_original_source_bridge'])
    table = weight_table()
    ranks = orbit_ranks(forms,ambient)
    result = dict(task_id='G',all_passed=True,no_producer_import=True,no_CAS_import=True,
                  complete_full25_source_bridge=bridge,full_weight_completeness=table,all_three_actual_orbits=ranks,
                  normal_form_count_at_most=3,type10_PGL_saturated_stratum_dimension_upper=22,
                  type10_not_generic_in_any_dimension25_26_27_component=True,
                  primitive_beta_image_nondominance_claimed=False,
                  A_B_component_membership_claimed=False,new_component_count_claimed=False,
                  completeness_read_review='Closed rational/affine field reduction retains every trace/Jordan/zero-eigen translation case. Auxiliary r=x-y^2/2 is used only to prove identities; completion R=r+L/delta is an actual affine x shift. All subsequent diagonalizations/translations/scales are affine/PGL. Potential constant subtraction changes no tensor.',
                  dimension_implication='A finite union of at most three actual22-dimensional PGL orbits has dimension<=22. It cannot contain a dense open of a component of dimension25,26 or27. This does not rule out a beta-dominant image of the22-dimensional stratum itself.',
                  input_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths+
                                [str(Path(q.__file__).relative_to(Path.cwd())),'code/classification_quadratic_cofactor_review.py']},
                  reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),
                  seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','maximum_RSS_KiB','type10_PGL_saturated_stratum_dimension_upper']}))


if __name__=='__main__':
    main()
