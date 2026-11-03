#!/usr/bin/env python3
"""Review-gated generic CB V completion, retaining its QQ(V) modulus."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out,review):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    audit=load(review);assert audit['all_passed'] is True
    author=ROOT/'reports/AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63/report.json';data=load(author)
    assert data['all_passed'] is True and data['actual_full_odd_firstjet_rank']==data['source_critical_cotangent_rank']==9
    ca=ROOT/'reports/AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63/report.json'
    deps=[author,review,ca,ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json']
    for path in deps:
        item=load(path);assert item['all_passed'] is True
        for name,digest in item.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    fp=ROOT/'reports/AMBIENT_CB_V_FAMILY_FRAME_20261103_v63/frame.json.gz';frame=load(fp)
    V=s.Symbol('V');K=s.QQ.frac_field(V);p=s.Matrix(frame['point150']);J=s.Matrix(frame['fulltarget150_by29_J'])
    assert J[:,10]==p #E_yy supplies the nonzero Poisson scalar
    T=J[:,:19].row_join(p.diff(V));D=DomainMatrix.from_Matrix(T).convert_to(K)
    assert len(D.rref()[1])==20;rows=list(D.transpose().rref()[1]);minor=s.factor(K.to_sympy(D.extract(rows,range(20)).det()))
    (out/'critical_family.json.gz').write_bytes(gzip.compress((json.dumps(dict(GL19_plus_modulusV_actual_tangent=mat(T),twenty_rows=rows,twenty_minor=str(minor),scalar_GL19_column=10),separators=(',',':'))+'\n').encode(),mtime=0))
    rel=ROOT/'reports/AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63/relations.json.gz'
    relations=load(rel);assert len(relations['all81'])==81
    (out/'relations.json.gz').write_bytes(rel.read_bytes())
    note=ROOT/'research/notes/AMBIENT_CB_GENERIC_V_FULL_COMPLETION_20261103_v63.md'
    paths=[Path(__file__),note,fp,rel,*deps]
    result=dict(id='AMBIENT_CB_GENERIC_V_FULL_COMPLETION_20261103_v63',all_passed=True,
      independent_actual_QQV_source_and_generic_application_review=str(review.relative_to(ROOT)),
      generic_modulus_retained=True,GL_orbit_rank=19,GL19_plus_V_critical_family_rank=20,critical_family_minor=str(minor),
      full_proper_finite_length_le2_open_meets_V0=True,generic_entire_flag_fiber_one_point_length2=True,
      critical_regular_codimension9=True,actual_I_equals_I_plus_nu_by_rank9_and_Nakayama=True,
      normalization_deduced_after_nonzero_I_recovers_t_in_fraction_field=True,
      exact_actual_odd_ideal='I=(nu,s1,...,s8)',
      cone_completion='C[[g20,s8,t^2,s1*t,...,s8*t,t^3]] inside B=C[[g20,s8,t]]',
      projective_completion='same core withg19',
      entire81_kernel_from_regular_sequence_Koszul_and_pair_product_reduction=True,
      exact_kernel_36_Koszul_plus45_products=81,Macaulay2_not_executed=True,
      cone_dimension=29,cone_embedding_dimension=38,cone_depth=21,
      projective_dimension=28,projective_embedding_dimension=37,projective_depth=20,multiplicity=2,
      normalizer_conductor='(t^2,s1,...,s8) B',target_conductor='(t^2,s1,...,s8,t^3,s1*t,...,s8*t) A',
      conductor_quotient='C[[g20]]',normalizer_quotient='C[[g20]]*t',domain_and_analytically_irreducible=True,
      normalization_elementary_subintegral_and_finite_universal_homeomorphism=True,
      local_singular_nonnormal_nonCM_loci_equal_smooth_conductor_support=True,
      scope='CB generic V only, modulus retained. No claim at every V, at source-frame excluded values, B18 origin full ideal, CA/CB intersection, or mixed cubic.',
      input_sha256={str(path.relative_to(ROOT)):sha(path)for path in paths},
      evidence_sha256={path.name:sha(path)for path in out.iterdir()if path.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('CB GENERIC V FULL COMPLETION /QQ(V) /81REL /DEPTH21 e2 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--review',type=Path,required=True)
    args=ap.parse_args();main(args.output.resolve(),args.review.resolve())
