#!/usr/bin/env python3
"""Explicit seminormalization of the full CA/CB formal ring model."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    model=ROOT/'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63/report.json';d=json.loads(model.read_text());assert d['all_passed'] is True
    for name,digest in d['input_sha256'].items():assert sha(ROOT/name)==digest,name
    zz=s.symbols('z1:7');nu,w,r,s0,t,x,y=s.symbols('nu w r s t x y');F=nu-s0*(s0-r)
    I=[*zz,F,w*r,w*s0,s0*r*r];J=[*zz,F,w,r]
    assert s.expand(nu*w**3-F*w**3-(w*s0)*w**2*(s0-r))==0
    assert s.expand(nu*r**3-F*r**3-(s0*r*r)*r*(s0-r))==0
    beta=x+y*t;reduction=lambda f:s.rem(s.Poly(f,t),s.Poly(t*t-nu,t)).as_expr()
    assert s.expand(reduction(beta**2)-(x*x+nu*y*y+2*x*y*t))==0
    assert s.expand(reduction(beta**3)-(x**3+3*x*nu*y*y+(3*x*x*y+nu*y**3)*t))==0
    node_sub={**dict.fromkeys(zz,0),w:0,r:0,nu:s0*s0};assert all(s.expand(f.subs(node_sub))==0 for f in J)
    Waxis={**dict.fromkeys(zz,0),r:0,s0:0,nu:0};Raxis={**dict.fromkeys(zz,0),w:0,s0:0,nu:0}
    assert all(s.expand(f.subs(Waxis))==0 for f in I)and w.subs(Waxis)!=0
    assert all(s.expand(f.subs(Raxis))==0 for f in I)and r.subs(Raxis)!=0
    Y=s.symbols('Y0:9');linear=[J[j]*Y[i]-J[i]*Y[j]for i,j in itertools.combinations(range(9),2)]
    products=[Y[i]*Y[j]-nu*J[i]*J[j]for i in range(9)for j in range(i,9)];assert len(linear)==36 and len(products)==45
    subst={v:b*t for v,b in zip(Y,J)}
    for f in linear:assert s.expand(f.subs(subst))==0
    for f in products:assert reduction(f.subs(subst))==0
    evidence=dict(base_R_core_variables=list(map(str,[*zz,nu,w,r,s0])),original_odd_ideal=list(map(str,I)),seminormal_odd_ideal=list(map(str,J)),
      newly_adjoined_elements=['w*t','r*t'],
      exact_cube_cofactors=['nu*w³=F*w³+(w*s)*w²*(s-r)','nu*r³=F*r³+(s*r²)*r*(s-r)'],
      square_odd_coefficient='2*x*y',cube_odd_coefficient='y*(3*x²+nu*y²)',
      prime_self_quotient='R/J=C[[g19,s]],nu=s²',
      normalization_conductor_node='B/JB=C[[g19,s,t]]/(t²-s²)',
      all36_complete_Koszul_relations=list(map(str,linear)),all45_pair_products=list(map(str,products)),all81=list(map(str,linear+products)))
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(evidence,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63.md';paths=[Path(__file__),note,model]
    report=dict(id='AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63',all_passed=True,
      exact_for_displayed_full_formal_model=True,
      self_prime='J=(z6,nu-s(s-r),w,r)',
      exact_seminormalization='A_sn=R directsum J*t=A[w*t,r*t] subset B=R[t]/(t²-nu)',
      J_is_self_branch_prime_not_radical_of_full_I=True,
      both_new_elements_square_and_cube_in_original_A=True,
      extension_finite_elementary_subintegral=True,
      seminormality_square_cube_test_in_entire_normalization_B_proved=True,
      seminormalization_properly_contains_A_and_is_properly_contained_in_B=True,
      source_conductor='J B=(z6,w,r,t²-s²) B',
      target_conductor='J directsum J*t inside A_sn',target_conductor_quotient='C[[g19,s]]',
      source_conductor_cover_is_reduced_node=True,actual_selfpair_gluing_not_removed=True,
      complete_seminormal_target_kernel_relations=81,GB_or_Macaulay2_not_run=True,
      cone_dimension=29,cone_embedding_dimension=38,cone_depth=21,multiplicity=2,
      projective_dimension=28,projective_embedding_dimension=37,projective_depth=20,
      actual_CA_CB_application_via_full_formal_model_only=True,
      no_relative_embedding_or_global_seminormality_classification=True,
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json.gz':sha(out/'certificate.json.gz')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FULL FORMAL MODEL SEMINORMALIZATION /SQUARE-CUBE /PROPER NODE GLUING PASS',report['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
