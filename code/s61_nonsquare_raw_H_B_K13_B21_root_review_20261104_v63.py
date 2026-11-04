#!/usr/bin/env python3
"""Independent replay of all21 fixed evaluations and all35 fixed minors."""
import hashlib,importlib.util,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_B21_ROOT_TEST_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    helper=ROOT/'code/s61_nonsquare_raw_H_B_K13_d_Euclid_review_20261104_v63.py';spec=importlib.util.spec_from_file_location('own_fixed_table_algebra',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
    pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_B21_ROOT_TEST_20261104_v63';r=read(pb/'report.json');assert r['all_passed'] and r['status']=='all35_minors_zero_candidate_ONLY'
    paths={Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in r['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    dp=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63_REVIEW/report.json';assert read(dp)['all_passed'];paths.add(dp)
    cp=read(pb/'checkpoint.json');values=read(pb/'B21-values.json');g=ar.ip(cp['g13']);assert ar.ip(values['g13'])==g
    A=ar.OwnAlgebra(g);s=A.dec(cp['root_s']);assert s==A.dec(values['root_s'])
    dr=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63/result.json';assert s==A.dec(read(dr)['root_s']);paths.add(dr)
    bb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_MIXED23_SCALARS_20261104_v63/native-certificate.json.B21.json'
    bridge=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW/report.json';assert read(bridge)['all_passed'];paths.update([bb,bridge])
    original=[A.poly(row)for row in read(bb)['entries']];assert len(original)==21 and sum(len([v for v in x[0]if v])for p in original for x in p)==789
    checked=[]
    for i,row in enumerate(values['entries']):
        assert row['entry']==i and [row['row'],row['column']]==[i//3,i%3]
        p=A.poly(row['input_coefficients']);assert p==original[i]
        p=p+[A.z]*(4-len(p));current=p[3];steps=row['Horner_steps'];assert len(steps)==3
        for j,step in zip([2,1,0],steps):
            assert step['d_coefficient']==j and A.dec(step['previous'])==current
            product=A.mul(current,s);assert product==A.dec(step['times_s']);current=A.add(product,p[j]);assert current==A.dec(step['after_add'])
        value=A.dec(row['value']);assert value==current
        q=A.poly(row['division_by_d_minus_s_quotient']);assert A.pa(A.pm([A.neg(s),A.o],q),[value])==A.trim(p)
        checked.append(value)
    assert len(checked)==21 and sum(sum(bool(v)for v in x[0])for x in checked)==273
    matrix=[checked[3*i:3*i+3]for i in range(7)]
    triples=list(itertools.combinations(range(7),3));assert len(cp['minors'])==len(triples)==35 and cp['unit']is None
    for row,ids in zip(cp['minors'],triples):
        assert row['rows7']==list(ids);b=[matrix[i]for i in ids]
        co=[A.sub(A.mul(b[1][1],b[2][2]),A.mul(b[1][2],b[2][1])),
          A.sub(A.mul(b[1][0],b[2][2]),A.mul(b[1][2],b[2][0])),
          A.sub(A.mul(b[1][0],b[2][1]),A.mul(b[1][1],b[2][0]))]
        assert co==[A.dec(v)for v in row['Laplace2x2_cofactors']]
        det=A.add(A.sub(A.mul(b[0][0],co[0]),A.mul(b[0][1],co[1])),A.mul(b[0][2],co[2]))
        assert det==A.dec(row['value'])==A.z and row['zero']
    proof={'all_passed':True,'arithmetic':'Own fixed-table dense13 QQ algebra, no author Ring or inverse/minor search.',
      'source_all84_B_coefficients_bound_to_independent_source_bridge':True,'all21_three_stage_Horner_and_synthetic_divisions':True,
      'all21_values_QQ_terms':273,'all35_fixed_rowtriples_Laplace_cofactors_and_zero_values':True,
      'whole_Ka13_no_field_irreducibility_or_reducedness_assumption':True,
      'scope':'The prepared7x3 matrix at d=s has all maximal minors zero in the whole Ka13 algebra. Actual Bbad existence or the exact actual bad scheme still requires explicit acceptance of the original minor/open source map. No such existence/emptiness conclusion, full cubic57 or quartic239 is inferred here.',
      'quotient_arithmetic_operations':A.ops,'new_solve_inverse_GCD_minor_search_GB_point_prime_jobs':0}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'all21_evaluations_all35_zero_minors_source_verified':True,'candidate_scope_only_not_actual_existence_or_empty':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
