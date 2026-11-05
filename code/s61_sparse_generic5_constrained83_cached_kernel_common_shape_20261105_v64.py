#!/usr/bin/env python3
"""SOURCE only: fixed scalar/monomial equality of14 cached denominators."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_CACHED_KERNEL_COMMON_SHAPE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for z in iter(lambda:f.read(65536),b''):h.update(z)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(map(numeric,v))
    if isinstance(v,dict):return sum(map(numeric,v.values()))
    return 0
def main():
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    sd=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FACTORED_CONSTRAINT_SCHUR_PREFLIGHT_20261105_v64';rp=sd/'report.json'
    assert sha(rp)=='8d0c175945ac170d8275da8fb2ecd1b7673d40dfd6c3bdc035f54d58b75912bf';r=read(rp);assert r['all_passed']and r['source_only']
    op=sd/'actual-factored177x24-and-top111x24.json';assert sha(op)==r['evidence_sha256'][op.name];o=read(op)
    paths={rp,op,Path(__file__),ROOT/f'research/notes/{ID}.md'};base=None;records=[];owned=0;checks=0;total_num=0
    for d in o['old_kernel_payload_index7']:
        p=ROOT/d['path'];paths.add(p)
        if len(records)==0:assert sha(p)==r['input_sha256'][d['path']]
        with p.open('rb')as f:f.seek(d['offset']);blob=f.read(d['bytes'])
        assert hashlib.sha256(blob).hexdigest()==d['payload_sha256'];v=json.loads(blob)
        for j,a in enumerate(v):
            den=a['denominator'];assert len(den)==290;terms={tuple(t[1:]):int(t[0])for t in den};assert len(terms)==290 and all(terms.values())
            shift=[min(e[k]for e in terms)for k in range(5)];relative={tuple(e[k]-shift[k]for k in range(5)):c for e,c in terms.items()}
            if base is None:base=relative;first=max(base,key=lambda e:(sum(e),e));baselead=base[first]
            assert set(relative)==set(base);scale=Q(relative[first],baselead)
            # Fixed integer coefficient cross-products only. No polynomial
            # multiplication, division, GCD, factorization or field job.
            for e,c in relative.items():assert c*scale.denominator==base[e]*scale.numerator;checks+=1
            records.append({'row_original0':d['row_original0'],'component':j,'kernel_node':d['components'][j]['node'],
                'denominator_QQ_scalar':str(scale),'denominator_parameter_exponent':shift,
                'literal_denominator_coefficient_crosschecks':len(relative),
                'common_numerator_recipe':{'original_numerator_component':j,'source_row_payload':d,'QQ_scalar':str(1/scale),'parameter_exponent':[-z for z in shift]},
                'numerator_terms_unchanged':len(a['numerator'])})
            total_num+=len(a['numerator']);owned=max(owned,numeric(o)+numeric(v)+5*len(base)+5*len(relative)+numeric(records)+8192);assert owned<150000
    assert len(records)==14 and total_num==4313 and checks==4060
    P={'content':'1','shift':[0]*5,'terms':[[str(base[e]),*e]for e in sorted(base,key=lambda e:(sum(e),e),reverse=True)]}
    pp=out/'actual-common-denominator-P290.json';pp.write_bytes(enc(P))
    result={'all_passed':True,'source_only':True,'variable_order':['H','L','s0','q0','t0'],'actual_common_polynomial':pp.name,'actual_common_polynomial_sha256':sha(pp),
        'common_polynomial_nonzero':True,'common_polynomial_terms':290,'all14_fixed_denominator_equalities':records,
        'all4060_literal_integer_crossproducts_verified':checks,'all4313_numerator_terms_unchanged':total_num,
        'exact_formula':'d_i=lambda_i*parameter^e_i*P290; cached T_i=(lambda_i^-1*parameter^-e_i*N_i)/P290. New numerators remain literal source/scale recipes, not re-expanded coefficients.',
        'factored_top_formula':'For each of four free columns: P290*old_free_entry + sum(old_pivot_entry*rescaled_cached_numerator). Divide by P290; all20 other columns stay unchanged.',
        'scope':'K5 equality; the cached alias chart requires P290 and the recorded parameter-monomial units nonzero, together with the actual constraint determinant chart. No complementary divisor or complete inverse claim.',
        'no_polynomial_GCD_or_coefficient_product_worker':True,'no_top_products_evaluated':True}
    ep=out/'fourteen-source-denominators-and-factored-P290.json';ep.write_bytes(enc(result))
    raw=sum(p.stat().st_size for p in paths)+pp.stat().st_size+ep.stat().st_size+1024*1024;assert raw<32*1024*1024
    report={'id':ID,'all_passed':True,'source_only':True,'actual_common_denominator_terms':290,'literal_crosschecks':4060,'cached_numerator_terms':4313,
        'source_owned_upper':owned,'RAW_upper':raw,'no_new_scalar_polynomial_product_or_kernel_query':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{pp.name:sha(pp),ep.name:sha(ep)},'seconds':time.monotonic()-begin}
    reportp=out/'report.json';reportp.write_bytes(enc(report));assert sum(p.stat().st_size for p in paths)+sum(p.stat().st_size for p in out.iterdir())<=raw
    signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(reportp))
if __name__=='__main__':main()
