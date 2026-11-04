#!/usr/bin/env python3
"""Constant-QQ source gauges for a proposed actual14-variable/17-polynomial E4.

Builds NO parameter-polynomial matrix, source jet, or Groebner basis.
"""
import argparse,gzip,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63'
h,x,y,z,w=X
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def em(M):return [[str(c) for c in r] for r in M.tolist()]
def scalar_bit(c):return max(abs(int(s.numer(c))).bit_length(),int(s.denom(c)).bit_length())
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    op=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json'
    O,E=read(op),read(ep);rows36=E['coordinate_rows36'];cuts=E['common_cut_zero_rows15']
    cubics=[s.prod((h,x,z,w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in m)]
    swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
    def change(P):return (swap*P.subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
    CA=[change(canonical(raw(f,0,0))) for f in [h*h,h*x,x*x]]+[change(canonical(raw(0,0,H))) for H in cubics]
    CB=[change(canonical(rawB(f,0,0))) for f in [h*h,h*x,x*x]]+[change(canonical(rawB(0,0,H))) for H in cubics]
    AC=s.Matrix.hstack(*map(vector,CA));BC=s.Matrix.hstack(*map(vector,CB));assert AC.shape==BC.shape==(150,19)
    P=s.Matrix([0,0,1,0,0])*s.Matrix([x*x,0,0,h*x,h*h+z*z]).T-s.Matrix([x*x,0,0,h*x,h*h+z*z])*s.Matrix([0,0,1,0,0]).T
    point=vector(P);assert list(point)==list(map(s.Rational,O['physical_point150'])) and canonical(P)==P
    A36=AC[rows36,:];B36=BC[rows36,:]
    ar=list(A36.T.rref()[1]);br=list(B36.T.rref()[1]);assert len(ar)==len(br)==19
    ai=A36[ar,:].inv();bi=B36[br,:].inv();ap=ai*point[[rows36[i] for i in ar],0];bp=bi*point[[rows36[i] for i in br],0]
    assert AC*ap==BC*bp==point
    annrows=[i for i in range(36) if i not in br];ann=s.zeros(17,36)
    for k,r in enumerate(annrows):
        ann[k,r]=1
        row=B36[r,:]*bi
        for j,c in enumerate(br):ann[k,c]-=row[j]
    assert ann.rank()==17 and ann*B36==s.zeros(17,19)
    gl4=[(i,j) for i in [0,1,3,4] for j in [0,1,3,4]]
    parabolic=[(i,j) for i,j in gl4 if not(i==1 and j!=1) and not(i==0 and j not in (0,1))]
    bad=[(1,0),(1,3),(1,4),(0,3),(0,4)]
    assert len(parabolic)==11 and sorted(parabolic+bad)==sorted(gl4)
    def Eij(i,j):M=s.zeros(5);M[i,j]=1;return M
    orbits={ij:vector(action(Eij(*ij),P)) for ij in gl4}
    orbit16=s.Matrix.hstack(*[orbits[ij] for ij in gl4]);assert orbit16.rank()==15
    para=s.Matrix.hstack(*[ai*orbits[ij][[rows36[i] for i in ar],0] for ij in parabolic])
    assert AC*para==s.Matrix.hstack(*[orbits[ij] for ij in parabolic]) and para.rank()==10
    picked=list(para.rref()[1]);p10=para[:,picked];fixed=list(p10.T.rref()[1]);free=[i for i in range(19) if i not in fixed]
    assert len(fixed)==10 and len(free)==9
    d10=s.factor(p10[fixed,:].det());assert d10
    chosen=bad+[parabolic[i] for i in picked];orbit15=s.Matrix.hstack(*[orbits[ij] for ij in chosen]);assert orbit15.rank()==15
    response=orbit15[cuts,:];d15=s.factor(response.det());assert d15
    # The source flag coordinates read by the reference frame have derivative
    # -bad(M). Composing its correcting frame cancels the bad Lie direction in
    # fiber coordinates, so the 15 source gauge derivatives are block diagonal.
    gauge=s.diag(-s.eye(5),p10[fixed,:]);assert gauge.det()!=0
    sliceC=AC[:,free];fiber_der=ann*sliceC[rows36,:]
    flag_der=ann*s.Matrix.hstack(*[orbits[ij][rows36,0] for ij in bad])
    J=fiber_der.row_join(flag_der);assert J.shape==(17,14) and J.rank()==5
    jc=list(J.rref()[1]);jr=list(J[:,jc].T.rref()[1]);d5=s.factor(J[jr,jc].det());assert d5
    # Both source sections contain scalar motion. GL4 center acts on a
    # quadratic ey-wedge tensor by -1; no projective quotient is silently taken.
    center=sum((orbits[i,i] for i in [0,1,3,4]),s.zeros(150,1));assert center==-point
    denA=s.ilcm(*[s.denom(c) for c in list(AC)+list(point)]);denN=s.ilcm(*[s.denom(c) for c in ann])
    boundA=max(abs(int(c*denA)) for c in list(AC)+list(point));boundN=max(abs(int(c*denN)) for c in ann)
    numerator_bound=40*36*4*36*boundA*boundN
    coefficient_degree=5;terms_per_coefficient=math.comb(5+coefficient_degree,5)
    terms_bound=17*10*terms_per_coefficient;assert terms_bound==42840 and terms_bound<150000
    ca,bb,cc,dd,ee=s.symbols('a b c d e');frame=s.eye(5);frame[0,3]=dd;frame[0,4]=ee;frame[1,0]=ca;frame[1,3]=bb;frame[1,4]=cc
    frameinv=frame.inv();assert frame.det()==1 and (frame*frameinv)==s.eye(5)
    assert max(s.Poly(q,ca,bb,cc,dd,ee).total_degree() for q in frameinv)<=2
    proof={'coordinates':list(map(str,X)),'all150_CA_fixed_flag19':em(AC),'all150_CB_fixed_flag19':em(BC),
      'source_cubic_monomials16':list(map(str,cubics)),'actual_P0_physical150':list(map(str,point)),
      'actual_CA_source19_center':list(map(str,ap)),'actual_CB_source19_center':list(map(str,bp)),
      'zero_coordinate_rows36':rows36,'CA_constant_left_inverse_rows19':ar,'CA_left_inverse19':em(ai),
      'CB_constant_left_inverse_rows19':br,'CB_left_inverse19':em(bi),'CB_annihilator17_by36':em(ann),
      'GL4_generator_positions16':gl4,'all150_GL4_orbit16':em(orbit16),'GL4_orbit_rank':15,
      'parabolic_positions11':parabolic,'CA_parabolic_coefficient_matrix19_by11':em(para),
      'chosen_parabolic10_column_indices':picked,'CA_fiber10_fixed_coefficient_indices':fixed,'CA_free_fiber9_coefficient_indices':free,
      'parabolic_fiber10_minor':str(d10),'source_gauge15_flag5_then_fiber10_Jacobian':em(gauge),'source_gauge15_minor':str(s.factor(gauge.det())),
      'actual_G15_generator_positions':chosen,'actual_common_cut15_rows':cuts,'actual_common_cut15_response':em(response),'actual_common_cut15_unit_minor':str(d15),
      'proposed_CA_source_slice':'CAcenter + sum nine free source basis columns * v_i; fixed10 coefficients unchanged',
      'relative_CB_flag_covectors':['ell=x+a*h+b*z+c*w','k=h+d*z+e*w'],
      'relative_CB_frame5_by5':em(frame),'relative_CB_inverse_frame5_by5':em(frameinv),'relative_frame_determinant':'1',
      'proposed_polynomial_equations':'CB_ann17 * coordinates36(frame_* (CAcenter + CAfree9*v)) = 0',
      'polynomial_Jacobian17_by14_at_origin':em(J),'Jacobian_rank':5,'Jacobian_pivot_rows':jr,'Jacobian_pivot_columns':jc,'Jacobian_unit5_minor':str(d5),
      'dimension_conventions':{'affine_CA_zero_source':24,'affine_CB_zero_source':24,'diagonal_G0_section':15,'CB_linear_fiber_eliminated':19,'variables_remaining':14,'equations_remaining':17},
      'GL4_center_derivative_equals_minus_P0':True,'no_extra_projectivization_or_full_GL4_16_quotient':True,
      'future_budget':{'coefficient_matrix_shape':[17,10],'relative_flag_variables':5,'max_flag_degree':5,'max_total_degree':6,'max_monomials_per_coefficient':252,'expanded_term_upper_bound':terms_bound,
        'constant_common_denominator_divides':str(denA*denN),'conservative_integer_numerator_abs_bound':str(numerator_bound),'numerator_bit_bound':int(numerator_bound).bit_length()},
      'parameter_polynomial_matrix_built':False,'no_GB_sourcejet_or_corehigher_coefficient_computed':True}
    p=out/'literal-proper-source14-preflight.json';p.write_text(json.dumps(proof,separators=(',',':'))+'\n')
    gates=[ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',op,ep,*gates,
      ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py']
    rep={'id':STEM,'all_passed':True,'status':'actual_QQ_source_gauges_and14_variable17_polynomial_E4_model_preflight',
      'fixed_flag_CA_CB_fiber_dimensions':[19,19],'GL4_orbit_rank':15,'parabolic_rank':10,
      'source_gauge15_minor':str(s.factor(gauge.det())),'actual_physical_cut15_minor':str(d15),
      'proposed_equation_Jacobian_rank':5,'unit5_minor':str(d5),'expanded_polynomial_term_bound':terms_bound,
      'parameter_polynomial_model_built':False,'allorder_gauge_equivalence_proof_uses_actual_unit15_matrices':True,
      'no_new_sourcejet_GB_fullcore_or_CM_claim':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
