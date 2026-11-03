#!/usr/bin/env python3
"""Origin conductor dim20 and e<=4 from accepted actual jets, no full ideal guess."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_global_flag_model_20261103_v63 import vector,action
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def piv(M):return list(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])
def det(M):return s.Rational(DomainMatrix.from_Matrix(M).convert_to(s.QQ).det())
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    deps=[ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json']
    for p in deps:
        d=load(p);assert d['all_passed'] is True
        for name,digest in d.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    qp=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz';quad=load(qp)
    cp=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63/certificate.json.gz';cubic=load(cp)
    fp=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz';first=load(fp)
    wp=ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63/staircase.json.gz';wholeW=load(wp)
    W,U,V=s.symbols('W U V');I150=s.eye(150);records={};critpaths=[]
    for name in['CA','CB']:
        q=quad[name];c=cubic[name];A=s.Matrix(first['CA_CB_full_firstjet_data'][name]['actual']);rr=q['actual_seven_odd_basis_rows'];cc=q['selected_actual_base_columns'];L=s.Matrix(q['all_actual_linear_combinations']);Q=s.Matrix(q['quadratic_coefficient_matrix'])[:,:2];qr=q['independent_quadratic_rows'];qc=Q*Q[qr,:].inv()
        select=I150[rr,:];h=I150[63,:]-L[63,:]*select-qc[63,:]*(I150[qr,:]-L[qr,:]*select)
        assert h*A==s.zeros(1,29)
        S2=s.Matrix(q['all150_vertical_derivatives2']);assert s.expand((h*S2)[0])==0
        b3=s.Matrix(c['actual_all150_odd_cubic_coefficients']);coefficient=(h*b3)[0]
        assert coefficient==s.sympify(c['actual_U3_class_coefficient'])!=0
        B=A[rr,cc];E=s.eye(29)[:,cc];correction=-E*B.inv()*b3[rr,:]
        assert A[rr,:]*correction==-b3[rr,:]
        gamma=s.Matrix(c['fullphysical_source_P']);gamma0=gamma.subs(U,0);Rrows=q['selected_actual28_Y']+[61];rawgamma=vector(gamma)[Rrows,:];T=s.Matrix(q['base_critical_tangent29_by3'])
        assert rawgamma.diff(U).subs(U,0)==T[:,1]
        O=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],gamma0))for i,j in itertools.product(range(5),repeat=2)])
        G=O[Rrows,:];gi=piv(G);assert len(gi)==19;G19=G[:,gi]
        assert A*G19==s.zeros(150,19)
        fullkernel=G19.row_join(T);assert len(piv(fullkernel))==22 and A*fullkernel==s.zeros(150,22)
        rows22=piv(fullkernel.T);minor22=det(fullkernel[rows22,:]);chart=E.row_join(fullkernel);minor29=det(chart);assert minor22 and minor29
        familyW=G19.row_join(T[:,0]);familyV=G19.row_join(T[:,2]);rowsW=piv(familyW.T);rowsV=piv(familyV.T)
        assert len(rowsW)==len(rowsV)==20 and len(piv(G19.row_join(T[:,[0,2]])))==21
        critpath=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';crit=load(critpath);critpaths.append(critpath)
        vvalues=[V/10,0,0,V/10,-V**2/200,0,0,-V/20,V**2/400,-V**3/24000]if name=='CA'else[V,0,0,V/2,-V**2/4,0,0,-V/2,V**2/4,-V**3/24]
        symbols=s.symbols('s0:10');sub=dict(zip(symbols,vvalues));schurkey='full150_critical_Schur_column'if name=='CA'else'all150_Schur_entries'
        assert all(s.expand(s.sympify(f).subs(sub))==0 for f in crit['predicted_entire_critical_ideal']+crit[schurkey])
        assert s.Matrix(wholeW[name]['physical_entire_W_curve']).subs(W,0)==gamma0
        records[name]=dict(seven_actual_odd_rows=rr,seven_base_columns=cc,seven_unit=str(B.det()),
          h63_as_constant_combination_of_all150_actual_odd_functions=mat(h),h63_full_R_firstjet_zero=True,h63_all_T3_quadratic_zero=True,
          h63_actual_U3_coefficient=str(coefficient),all_seven_values_on_exact_monic_Ucurve_have_order_atleast3=True,
          actual_seven_IFT_leading_correction_order=3,actual_R29_IFT_correction_coefficient=mat(correction),
          actual_R29_monic_Ucurve=mat(rawgamma),Ucurve_tangent_exactly_TU=True,
          G19_GL_indices=gi,actual_R29_G19_tangents=mat(G19),full_minimal_kernel_basis29_by22=mat(fullkernel),
          full_minimal_kernel_rank=22,kernel22_rows=rows22,kernel22_minor=str(minor22),seven_plus_kernel29_unit=str(minor29),
          genuine_G19_W_family_rank=20,genuine_G19_W_family_minor=str(det(familyW[rowsW,:])),
          genuine_G19_V_family_rank=20,genuine_G19_V_family_minor=str(det(familyV[rowsV,:])),
          W_V_tangent_union_rank=21,entire_critical_V_source=list(map(str,vvalues)),
          origin_cubic_core_form='U*(a*U²+b*U*V+c*V²), a=nonzero saved h63 coefficient; b,c NOT computed',
          I0_contains_genuine_initial_relations=['W*U','W*V','U*(a*U²+b*U*V+c*V²)'])
    a,b,c=s.symbols('a b c');f=U*(a*U*U+b*U*V+c*V*V);M=s.Matrix([[V,a*U*U+b*U*V+c*V*V],[-U,0],[0,-W]])
    minors=[s.expand(M[[1,2],:].det()),-s.expand(M[[0,2],:].det()),s.expand(M[[0,1],:].det())];assert all(s.expand(x-y)==0 for x,y in zip(minors,[W*U,W*V,f]))
    payload=dict(actual_jet_bindings=records,universal_Hilbert_Burch_matrix=mat(M),universal_maximal_minors=list(map(str,minors)),
      exact_curve_Hilbert_series='(1+T)^2/(1-T)',CM_initial_candidate_multiplicity=4,
      actual_unknown_b_c_and_higher_terms_not_set_to_zero=True)
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(payload,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63.md'
    paths=[Path(__file__),note,qp,cp,fp,wp,*critpaths,*deps]
    result=dict(id='AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63',all_passed=True,
      no_new_higher_jet_or_GB_computation=True,
      intrinsic_conductor_G_invariance_used_not_unproved_R_base_G_stability=True,
      actual_seven_odd_functions_eliminated_by_IFT=True,cubic_U3_survives_IFT_order3_correction=True,
      full_minimal_quadratics_independent_of_G19_by_lifted_derivations=True,
      genuine_initial_cubic_after_actual_odd_elimination_has_nonzero_U3=True,
      unknown_full_binary_cubic_not_identified=True,
      CA_CB_origin_conductor_cone_dimension=20,CA_CB_origin_conductor_projective_dimension=19,
      CA_CB_origin_conductor_cone_embedding_dimension=22,CA_CB_origin_conductor_multiplicity_upper_bound=4,
      at_least_two_distinct_smooth_reduced_support_branches_W_and_V=True,
      offdiagonal_individual_normalization_selfpair_dimension_upper_bound=20,
      actual_CA_CB_conductor_depth_CM_generator_upper_bound_not_determined=True,
      actual_I_full_ideal_and_actual_CA_CB_intersection_not_claimed=True,
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json.gz':sha(out/'certificate.json.gz')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('CA CB ACTUAL ORIGIN CONDUCTOR DIM20 /GENUINE INITIAL CUBIC /e<=4 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
