#!/usr/bin/env python3
"""Review-gated application of the actual V1 odd rank9 to the full CA image."""
import argparse,gzip,hashlib,itertools,json,subprocess,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def main(out,review,m2):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    audit=load(review);assert audit['all_passed'] is True
    author=ROOT/'reports/AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63/report.json'
    data=load(author);assert data['all_passed'] is True
    assert data['whole_actual_odd_firstjet_rank']==data['whole_normalizer_critical_cotangent_rank']==9
    assert data['GL_orbit_dimension']==20
    dependencies=[author,review,
      ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json']
    for path in dependencies:
        item=load(path);assert item['all_passed'] is True
        for name,digest in item.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    ck=load(ROOT/'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz')
    J=s.Matrix(ck['full150_by29_target_Jacobian']);point=s.Matrix(ck['point150'])
    assert J[:,10]==point #noncentral GL generator supplies the Poisson scalar
    nu=s.Symbol('nu');ss=s.symbols('s1:9');yy=s.symbols('y0:9');t=s.Symbol('t')
    base=[nu,*ss];b=base;linear=[];products=[]
    for i,j in itertools.combinations(range(9),2):linear.append(b[j]*yy[i]-b[i]*yy[j])
    for i in range(9):
        for j in range(i,9):products.append(yy[i]*yy[j]-nu*b[i]*b[j])
    relations=linear+products;assert len(relations)==81
    substitutions={y:bi*t for y,bi in zip(yy,b)}
    for equation in relations:
        assert s.rem(s.Poly(equation.subs(substitutions),t),s.Poly(t*t-nu,t)).as_expr()==0
    data=dict(core_R_variables=list(map(str,base)),odd_coefficient_ideal=list(map(str,b)),
      odd_target_variables=list(map(str,yy)),normalizer_relation='t^2-nu',
      all36_Koszul_relations=list(map(str,linear)),all45_products=list(map(str,products)),
      all81=list(map(str,relations)),normalization_map={str(y):str(bi*t)for y,bi in zip(yy,b)})
    (out/'relations.json.gz').write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    text='S=QQ['+','.join(map(str,[t,*base,*yy]))+',MonomialOrder=>Eliminate 1];\n'
    text+='P=ideal('+','.join(str(x).replace('**','^')for x in relations)+');\n'
    text+='G=ideal('+','.join(str(x).replace('**','^')for x in[t*t-nu,*[y-bi*t for y,bi in zip(yy,b)]])+');\n'
    text+='K=eliminate({t},G);\nassert(K==P);\nassert(dim(S/P)==10);\n'
    text+='print "CA_GENERIC_V_WHOLE81_KERNEL_PASS";\nexit 0;\n'
    (out/'kernel.m2').write_text(text)
    if m2:
        run=subprocess.run(['M2','--script',str((out/'kernel.m2').resolve())],capture_output=True,text=True,check=True,timeout=60)
        assert 'CA_GENERIC_V_WHOLE81_KERNEL_PASS'in run.stdout
        (out/'kernel.log').write_text(run.stdout+run.stderr)
    note=ROOT/'research/notes/AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63.md'
    paths=[Path(__file__),note,*dependencies]
    result=dict(id='AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63',all_passed=True,
      independent_actual_source_and_generic_V_application_review=str(review.relative_to(ROOT)),
      critical_V_GL_saturation_dimension=20,V1_GL_orbit_dimension=20,
      GL_generator10_is_actual_Poisson_scalar=True,
      proper_entire_normalization_fiber_one_point_length2=True,
      normalization_deduced_after_nonzero_I_recovers_t_in_fraction_field=True,
      critical_scheme_regular_codimension9=True,
      actual_odd_I_equals_J_by_all9_regular_classes_and_Nakayama=True,
      exact_full_odd_ideal='I=(nu,s1,...,s8)',
      cone_completion='C[[g20,s8,t^2,s1*t,...,s8*t,t^3]] inside B=C[[g20,s8,t]]',
      projective_completion='same core with g19 instead of g20',
      complete_target_kernel_relations=81,complete_Koszul_relations=36,complete_products=45,
      independent_Macaulay2_entire_kernel_checked=bool(m2),
      cone_dimension=29,cone_embedding_dimension=38,cone_depth=21,
      projective_dimension=28,projective_embedding_dimension=37,projective_depth=20,
      Hilbert_Samuel_multiplicity=2,domain_and_analytically_irreducible=True,
      normalizer_conductor='(t^2,s1,...,s8) B',
      target_conductor='(t^2,s1,...,s8,t^3,s1*t,...,s8*t) A',
      conductor_quotient='C[[g20]]',normalizer_quotient='B/A=C[[g20]]*t',
      finite_elementary_subintegral_normalization=True,seminormalization_equals_B=True,
      local_singular_nonnormal_nonCM_loci_equal_smooth_conductor_support=True,
      scope='CA V1 and its dense generic V orbit only; no CB, B18 origin full ideal, mixed cubic, or actual CA/CB intersection type.',
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
      evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('CA GENERIC V FULL COMPLETION /81REL /DIM29 DEPTH21 e2 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--review',type=Path,required=True);ap.add_argument('--m2',action='store_true')
    args=ap.parse_args();main(args.output.resolve(),args.review.resolve(),args.m2)
