#!/usr/bin/env python3
"""ONE exact named-unit strip of seven fixed normalizers; no factor search."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PURE_KERNEL9_NORMALIZER_STRIP_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(line):
    t=line.split();n=int(t[0]);assert len(t)==1+4*n
    return {tuple(map(int,t[2+4*k:5+4*k])):Q(t[1+4*k])for k in range(n)}
def enc(p):return [[list(m),str(v)]for m,v in p.terms()]
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    pb=ROOT/'reports/S61_NONSQUARE_PURE_KERNEL9_CHARTS_20261104_v63';assert sha(pb/'report.json')=='7186322f73ff3be7c31bb8f7ed43d28f29c8fbc9c117e8e6e1cdf954abdb0e6d'
    c=read(pb/'localizer-inventory.json');occ=[x for x in c['all136_inverse_occurrences']if x['kind']=='NORMALIZE_g'and x['classification'].startswith('UNKNOWN')];assert len(occ)==26
    a,e,d=s.symbols('a e d');gens=(a,e,d);pol=lambda q:s.Poly.from_dict({m:s.Rational(v.numerator,v.denominator)for m,v in q.items()},gens,domain=s.QQ)
    source=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63/operations.txt';unique={};records=[];terms=0
    with source.open()as f:
        for x in occ:
            f.seek(x['byte_offset']);raw=parse(f.readline());p=pol(raw);terms+=len(raw);monic=p.monic();key=tuple(monic.terms())
            if key not in unique:unique[key]=(len(unique),monic)
            records.append(dict(source_byte_offset=x['byte_offset'],source_operation=x['operation'],source_phase=x['introduced_phase'],
                projective_polynomial_index=unique[key][0],source_QQ_unit=str(p.LC()),source_polynomial=enc(p)))
    assert terms==590 and len(unique)==7 and sum(len(p.terms())for i,p in unique.values())==374
    fixed=[('a',s.Poly(a,gens,domain=s.QQ)),('d',s.Poly(d,gens,domain=s.QQ)),
        ('F',s.Poly(500*a+90*e**2+3*e-24,gens,domain=s.QQ)),
        ('B',s.Poly((50*a+18*e**2-9*e-2)*d+60*e-40,gens,domain=s.QQ)),
        ('C',s.Poly(100*a+18*e**2-9*e-2,gens,domain=s.QQ)),
        ('J',s.Poly(5*a*d-6*e-4,gens,domain=s.QQ))]
    result=[];explicit=0
    for idx,p in unique.values():
        cur=p;steps=[];powers=[]
        for name,f in fixed:
            n=0
            while cur.total_degree()>=f.total_degree():
                q,r=cur.div(f)
                if not r.is_zero:break
                assert f*q==cur
                steps.append(dict(named_unit=name,before=enc(cur),quotient=enc(q),exact_product_all_coefficients_checked=True));cur=q;n+=1
            powers.append(n)
        constant=cur.LC();rem=cur.monic();reconstructed=rem.mul_ground(constant)
        for(name,f),n in zip(fixed,powers):
            for _ in range(n):reconstructed*=f
        assert reconstructed==p
        explicit+=len(p.terms())+len(rem.terms())+sum(len(z['before'])+len(z['quotient'])for z in steps);assert explicit<150000
        result.append(dict(index=idx,source_monic_polynomial=enc(p),named_units_order=[n for n,f in fixed],named_unit_powers=powers,
            QQ_unit=str(constant),remaining_monic_polynomial=enc(rem),remaining_is_QQ_unit=rem.total_degree()==0,
            exact_division_steps=steps,whole_product_reconstruction_checked=True,no_unknown_factor_removed=True))
        (out/'checkpoint.json').write_text(json.dumps(dict(completed=len(result),active_explicit_terms=explicit,no_auto_other_factorization=True),indent=2)+'\n')
    leftovers={tuple(tuple(m)+tuple([v])for m,v in r['remaining_monic_polynomial'])for r in result if not r['remaining_is_QQ_unit']}
    cert=dict(source26_occurrences=records,seven_distinct_normalizers=result,named_unit_polynomials={n:enc(f)for n,f in fixed},
        named_only_exactdivisions=True,general_factorization_GCD_rank_minor_point_prime_NOT_run=True,
        unknown_remaining_polynomials_preserved=True,distinct_nonconstant_leftovers=len(leftovers))
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    proof=dict(domain='FULLU=adFBCJ !=0',source_forward_chart='D(product38NORMALIZE_g), includingthe26formerlyUNKNOWNoccurrences',
        literal_each_occurrence='source_g=source_QQ_unit*remaining_QQ_unit*product_named_units*remaining_monic_polynomial',
        consequence='OnFULLU,exactlythe remaining_polynomial localizers must be retained for this forward source recipe. No claim those zero sets are actualcubicdefect or fullycovered.',
        old_equivalence_chart_136_plus15_UNCHANGED=True,no_automatic_chart_extension_or_rank_claim=True)
    (out/'chart-consequence.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',pb/'report.json',pb/'localizer-inventory.json',source]
    report=dict(id=ID,all_passed=True,source_unknown_occurrences=26,source_distinct_polynomials=7,source_distinct_terms=374,
        active_explicit_terms=explicit,named_unit_division_steps=sum(len(r['exact_division_steps'])for r in result),
        distinct_nonconstant_leftovers=len(leftovers),no_unknown_factorization_coefficient_rank_or_other_job=True,
        independent_review_pending=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),no_auto_retry=True),indent=2)+'\n')
        raise
