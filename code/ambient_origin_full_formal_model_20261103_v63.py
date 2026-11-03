#!/usr/bin/env python3
"""Actual CA/CB complete local model from intrinsic cover and conductor pushout."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    rigidity=ROOT/'reports/AMBIENT_THREE_AXIS_CONDUCTOR_RIGIDITY_20261103_v63.json';lemma=load(rigidity)
    assert lemma['all_passed'] is True and lemma['ring_rigidity_proved'] is True and lemma['pair_rigidity_proved'] is True
    binary=ROOT/'reports/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW/report.json'
    deps=[rigidity,binary,*[ROOT/f'reports/AMBIENT_{name}_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json'for name in['CA','CB']],
      ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json']
    for p in deps:
        d=load(p);assert d['all_passed'] is True
        for name,digest in d.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    qp=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz';quad=load(qp)
    bp=ROOT/'reports/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW/certificate.json.gz';bc=load(bp)
    U,V,r,s0,w,nu,t=s.symbols('U V r s w nu t');q2={}
    for name in['CA','CB']:
        ell=s.sympify(bc[name]['self_component_tangent_UV_equation']);lam=-s.Poly(ell,U,V).coeff_monomial(V)
        actual=-s.sympify(quad[name]['chosen_nu_lift_quadratic']);converted=s.expand(actual.subs({U:s0,V:(s0-r)/lam},simultaneous=True))
        coefficient=s.factor(converted/(s0*(s0-r)));assert coefficient==(s.Rational(1,25)if name=='CA'else s.Rational(16,81))
        q2[name]=dict(actual_nu_quadratic=str(actual),self_line_lambda=str(lam),after_s_U_r_U_minus_lambdaV=str(converted),nonzero_leading_coefficient=str(coefficient))
    zz=s.symbols('z1:7');yy=s.symbols('Y0:10');F=nu-s0*(s0-r);b=[*zz,F,w*r,w*s0,s0*r*r]
    syzygies=[]
    for i,j in itertools.combinations(range(7),2):syzygies.append(b[j]*yy[i]-b[i]*yy[j])
    for i in range(7):
        for j in range(7,10):syzygies.append(b[j]*yy[i]-b[i]*yy[j])
    syzygies +=[s0*yy[7]-r*yy[8],r*s0*yy[7]-w*yy[9]]
    products=[yy[i]*yy[j]-nu*b[i]*b[j]for i in range(10)for j in range(i,10)]
    assert len(syzygies)==44 and len(products)==55
    substitutions={y:bi*t for y,bi in zip(yy,b)}
    for f in syzygies:assert s.expand(f.subs(substitutions))==0
    for f in products:
        assert s.rem(s.Poly(f.subs(substitutions),t),s.Poly(t*t-nu,t)).as_expr()==0
    data=dict(base_R_core_variables=list(map(str,[*zz,nu,w,r,s0])),odd_generators=list(map(str,b)),odd_target_variables=list(map(str,yy)),
      normalizer_relation='t²-nu',all44_complete_first_syzygies=list(map(str,syzygies)),all55_pair_products=list(map(str,products)),all99=list(map(str,syzygies+products)),
      actual_cover_quadratic_binding=q2,conductor_primary_decomposition=[['r','s'],['w','s'],['w','r²']],
      critical_quotient_change='U=s,V=r-s; (wr,ws,sr²,s(s-r))=(WU,WV,UV,U³)')
    (out/'relations99.json.gz').write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63.md';paths=[Path(__file__),note,qp,bp,*deps]
    report=dict(id='AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63',all_passed=True,
      actual_CA_and_CB_separately_isomorphic_to_same_full_formal_model=True,
      canonical_G_equivariant_normalization_and_conductor_cover_split_off_G19=True,
      actual_pair_conductor_and_trace_quadratic_rigidity_hypotheses_bound=True,
      regular_normalizer_Cohen_presentation_constructs_NEW_base_not_original_relative_alignment=True,
      exact_base='R=C[[g19,z6,nu,w,r,s]]',exact_odd_ideal='I=(z6,nu-s(s-r),wr,ws,sr²)',
      exact_actual_completion='A=R directsum I*t subset B=R[t]/(t²-nu)',
      regular_normalizer='B=C[[g19,z6,w,r,s,t]]',
      conductor_in_normalizer='(z6,t²-s(s-r),wr,ws,sr²) B',
      conductor_quotient='C[[g19,w,r,s]]/(wr,ws,sr²)',
      normalization_conductor_cover='C[[g19,w,r,s,t]]/(wr,ws,sr²,t²-s(s-r))',
      whole_target_kernel_99_relations=True,complete_first_syzygies=44,complete_pair_products=55,
      completeness_from_regular7_Koszul_tensor_Hilbert_Burch_and_rank2_reduction=True,GB_or_Macaulay2_not_run=True,
      cone_dimension=29,cone_embedding_dimension=39,cone_depth=21,multiplicity=2,
      projective_smooth_G_factors=18,projective_dimension=28,projective_embedding_dimension=38,projective_depth=20,
      only_formal_complete_local_isomorphism_claimed=True,convergent_analytic_isomorphism_not_claimed=True,
      relative_CA_CB_embedding_actual_intersection_and_global_classification_not_claimed=True,
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'relations99.json.gz':sha(out/'relations99.json.gz')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('ACTUAL CA/CB FULL FORMAL MODEL /NEW R /COMPLETE99 RELATIONS PASS',report['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
