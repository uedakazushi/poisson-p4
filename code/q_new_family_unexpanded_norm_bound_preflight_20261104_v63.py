#!/usr/bin/env python3
"""Support-only finite10 norm polynomial recipe; no polynomial products or norm."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_UNEXPANDED_NORM_BOUND_PREFLIGHT_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    recipepath=ROOT/'reports/Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63/recipe.json.gz'
    recipe=read(recipepath)
    gates=[ROOT/'reports/Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63_REVIEW/report.json',
           ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_PRIMITIVITY_20261104_v63_REVIEW/report.json',
           ROOT/'reports/Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63_REVIEW/report.json',
           ROOT/'reports/Q_FIRSTJET_SIX_LOW_INVARIANT_INDEPENDENCE_20261104_v63/report.json',
           ROOT/'reports/Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    G=recipe['G_terms_embedded'];C=recipe['C1_terms_constant_then_linear_embedded']
    # Three support degrees: r; R-Z0; 2*(Z0+Z1+Z2).
    gdeg=[max(q[2]for q in G),max(q[3]-q[4]for q in G),max(2*sum(q[4:7])for q in G)]
    cdeg=[[max(q[5]for q in p),max(q[4]-q[3]for q in p),max(2*sum(q[1:4])for q in p)]for p in C]
    assert gdeg==[10,0,10]and cdeg==[[31,8,10],[31,8,10]]
    small=recipe['N79_N41_Fn_J_embedded'];hom={}
    for name,poly in small.items():
        terms=poly['coefficients'];ds=max(q[3]for q in terms)
        hom[name]=[max(q[2]for q in terms)+31*ds,max(q[4]for q in terms)+8*ds,10*ds]
    assert hom=={'N79':[99,27,30],'N41':[67,18,20],'Fn':[67,18,20],'J':[32,8,10]}
    def add(*args):return [sum(p[i]for p in args)for i in range(3)]
    U1=add(hom['J'],hom['N79'],[3,-1,2])
    U2=add(hom['Fn'],hom['N41'],[0,-2,4])
    V=add(hom['Fn'],hom['J'],hom['J'],[3,0,0])
    U=[max(a,b)for a,b in zip(U1,U2)]
    q=[max(U[0],V[0]),max(U[1],V[1]),max(U[2],V[2]+3)]
    assert U==[134,34,44]and V==[134,34,40]and q==[134,34,44]
    n,m=10,q[0];A=m*gdeg[1]+n*q[1];K=m*gdeg[2]+n*q[2];degree=8*A+4*K
    assert (A,K,degree)==(340,1780,9840)
    proof=dict(coefficient_base='K0=QQ(R,Z0,Z1,Z2)',finite_field='L=K0[r]/(monic G10)',
        G_source_terms=len(G),C1_source_terms=list(map(len,C)),
        support_degree_order=['degree_r','max(exponent_R-exponent_Z0)','max(2*sum_Z_exponents)'],
        G_degrees=gdeg,C1_constant_linear_degrees=cdeg,small_homogenized_bounds=hom,
        homogenization='H_d[f]=A_lin^d*f(r,-B_const/A_lin,R), s=-B_const/A_lin from the accepted primitive C1',
        balanced_U='Z0*r*(r-1)^2*JH*N79H+25*Z0^2*FnH*N41H',
        balanced_V='20*r*(r-1)^2*FnH*JH^2',U_degree_bounds=U,V_degree_bounds=V,
        formal_eta_weight_for_support=3,eta_V_minus_U_bounds=q,
        finite_operator='P(eta)=lc_r(G)^134 det_10(eta*M_V-M_U), multiplication in the monic G10 companion algebra',
        polynomialness='The root-product formula identifies this with the fixed-degree(10,134) Sylvester determinant, including the lc(G) padding power if the second degree is smaller. No determinant is expanded here.',
        nonzero='V is nonzero in the accepted degree10 field; the eta^10 coefficient is lc(G)^134 Norm(V)!=0. Primitivity additionally makes the degree10 norm minimal over K0 up to a base factor.',
        source_vanishing='Actual eta=N12/F4^3=U/V in Kshape, hence P(eta) vanishes on actual NEW.',
        polynomial_support_bounds=dict(R_minus_Z0=A,two_Z_sum_plus_three_eta=K),
        target_substitution=dict(R='(p3+4*p4)/p0',Z0='p0/F4^2',Z1='p1/F4^2',Z2='p2/F4^2',eta='N12/F4^3',p0='F4^2/400-4*p4'),
        explicit_target_polynomial='H=p0^340*F4^1780*P((p3+4p4)/p0,p0/F4^2,p1/F4^2,p2/F4^2,N12/F4^3), then substitute p0=F4^2/400-4p4',
        clearing_proof='Each P monomial R^a Z0^b Z1^c Z2^d eta^e has a-b<=340 and 2(b+c+d)+3e<=1780. Its substituted p0 exponent340+b-a and F4 exponent1780-2(b+c+d)-3e are nonnegative.',
        nonzero_universally='The six selected low invariants are algebraically independent; the rational base change with F4,R,Z0,Z1,Z2,eta is invertible generically. Clearing by nonzero F4,p0 cannot turn P into zero.',
        homogeneous_beta_degree=degree,minimal_NEW_defining_degree_upper=degree,
        irreducible_factor_reason='The actual geometric image is an irreducible hypersurface. Its homogeneous prime generator divides the nonzero invariant H.',
        not_minimal_equation='H may have other ambient hypersurface factors because the low-six map is not finite/HSOP and has generic degree8. Base content and boundary factors need separate isolation.',
        finite_isolation_strategy='Factor H over QQ and test each factor on the accepted dense actual rational source; exactly one irreducible factor vanishes on actual NEW. This is a finite recipe, not a factorization job.',
        ordered_genus_poles_not_used_to_force_degree=True,
        no_expanded_norm_resultant_products_newpoint_prime_or_equation_search=True)
    ep=out/'preflight.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),recipepath,*gates]
    report=dict(id=ID,all_passed=True,status='unexpanded_finite10_operator_and_support_degree_preflight',
        explicit_nonzero_invariant_multiple_of_actual_NEW_defined_by_finite_operators=True,
        concrete_homogeneous_degree_upper=9840,minimal_defining_degree_upper_bound=9840,
        actual_minimal_equation_or_expanded_coefficients_not_constructed=True,
        no_norm_resultant_expansion_source_point_prime_or_new_invariant_job=True,
        independent_review_pending=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start,frozen=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'degree_upper':9840,'seconds':report['seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
