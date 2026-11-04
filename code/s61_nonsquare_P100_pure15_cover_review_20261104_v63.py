#!/usr/bin/env python3
"""Source-bound fixed fifteen-minor/Bezout cover, independent FLINT replay."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_PURE15_COVER_20261104_v63_REVIEW'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);sys.set_int_max_str_digits(0)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);work=ROOT/'work/P100-pure15-independent-review';work.mkdir(parents=True,exist_ok=True);hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def gate(stem,passed=True):
        p=ROOT/'reports'/stem/'report.json';r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            q=p.parent/name;assert sha(q)==h,name;hashes[str(q.relative_to(ROOT))]=h
        return r
    b=ROOT/'reports'
    oldpeer=gate('S61_NONSQUARE_P100_G18_PURE8_20261104_v63_REVIEW')
    newpeer=gate('S61_NONSQUARE_P100_K2_PURE7_20261104_v63_REVIEW')
    assert hashes['reports/S61_NONSQUARE_P100_G18_PURE8_20261104_v63_REVIEW/report.json']=='8ef758e0872b7e8b54b478fe82e275daffc24605f1b92a80ff17d472e373ae39'
    assert hashes['reports/S61_NONSQUARE_P100_K2_PURE7_20261104_v63_REVIEW/report.json']=='ea4e71d95e47a71069392aff8ac64170ab6e5681e7de2cbc8485b41aa21e77b3'
    base=b/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63';r=gate(base.name,False)
    assert r['complete'] and r['whole_P100_quartic_pointwise_lower149']
    assert sha(base/'report.json')=='3a7e735027b75f7e9c1ac6d7c77615710da0db387f2129c0f4ffe77c8e8cf6d2'
    native=read(base/'native.json');binding=read(base/'source-binding.json');assert sha(base/'native.json')=='b2811cdb0a0651339313ac21b98e88bd7a64b1170455a986c7db724de18d126c'
    oldb=b/'S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63'
    old=read(oldb/'native.json');oldbinding=read(oldb/'selected-source-binding.json')
    newb=b/'S61_NONSQUARE_P100_K2_PURE7_20261104_v63'
    newbinding=read(newb/'source-binding.json');rows=[read(newb/f'pure98-row{i}.json') for i in range(7)]
    pf=b/'S61_NONSQUARE_P100_PURE15_COVER_PREFLIGHT_20261104_v63';pr=gate(pf.name,False);pc=read(pf/'certificate.json')
    assert sha(pf/'certificate.json')==pr['certificate_sha256']
    assert binding['old_selected_source_binding']==oldbinding and binding['new_corrected_K2_source_binding']==newbinding
    columns=old['selected_columns'];assert columns==pc['fixed_target_local98_columns8']==[0,1,2,3,4,5,6,9]
    assert binding['fixed_old_selected_target98_columns']==columns
    assert binding['fixed_original122_target_labels']==[0,1,2,4,5,6,7,10]
    assert binding['fixed_original122_target_labels']==[rows[0]['target122_columns'][i] for i in columns]
    powers=[rec['new_known_unit_LCD_powers5'] for rec in pc['profiles7']]
    for i,(row,rec,ex) in enumerate(zip(rows,binding['new_sources_in_order'],powers)):
        assert rec==dict(kernel=2,actual_nonpivot135=23,physical_light=row['physical_light'],known_unit_rowLCD_powers5=ex)
        assert row['target122_columns']==rows[0]['target122_columns']
    assert len(native['candidates'])==native['computed_fixed_replacements']==15
    assert [(rec['K2_light_position'],rec['replaced_old_selected_row_position']) for rec in native['candidates']]==[(0,i) for i in range(8)]+[(1,i) for i in range(7)]
    assert [c['basis_index'] for c in native['literal_primitive_Bezout_contributors']]==[0,7,9,10,15]
    assert [c['basis_index'] for c in native['actual_minor_polynomial_Bezout_contributors']]==[0,7,9,10,15]
    cpp=Path(__file__).with_suffix('.cpp')
    for p in [Path(__file__),cpp,ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    info=dict(source_native_bytes=(base/'native.json').stat().st_size,
        fixed_adjugate_products=128,selected_new_RF_entries=56,fixed_replacements=15,
        primitive_and_actual_Bezout_contributors=5,no_gcd_rank_search_new_source_or_producer=True,
        large_QQ_strings_forwarded_directly_to_FLINT=True,
        resource_caps=dict(seconds=60,address_bytes=2_000_000_000,output_bytes=32*1024*1024,live_with_literal_reserve=150000))
    (out/'preflight-and-source-binding.json').write_text(json.dumps(info,indent=2)+'\n')
    inputfile=work/'fixed-minors-and-Bezout-input.txt'
    with inputfile.open('w') as f:
        def poly(p):
            f.write(str(len(p))+'\n')
            for m,c in p:f.write(f'{m[0]} {c}\n')
        def split(s):
            for k in ['full','primitive','unit_polynomial']:poly(s[k])
            f.write(s['QQ_unit']+'\n'+' '.join(map(str,s['five_unit_multiplicities']))+'\n')
        for a,c in [(1,0),(1,20),(1,60),(203,2700),(1,12)]:poly([[[1],str(a)]]+([[[0],str(c)]] if c else []))
        poly(old['determinant_nonzero'])
        for M in [old['selected_square'],old['final_adjugate']]:
            for rr in M:
                for p in rr:poly(p)
        for ex in powers:f.write(' '.join(map(str,ex))+'\n')
        for rr in rows:
            for i in columns:poly(rr['values98'][i]['numerator']);poly(rr['values98'][i]['denominator'])
        split(native['old_minor_strip']);f.write('15\n')
        for cand in native['candidates']:
            f.write(f"{cand['K2_light_position']} {cand['replaced_old_selected_row_position']}\n");split(cand['split'])
        poly(native['final_gcd']);f.write(' '.join(map(str,native['common_actual_minor_unit_powers5']))+'\n');poly(native['common_actual_minor_unit_polynomial'])
        for field in ['literal_primitive_Bezout_contributors','actual_minor_polynomial_Bezout_contributors']:
            f.write('5\n')
            for rec in native[field]:f.write(str(rec['basis_index'])+'\n');poly(rec['coefficient'])
    assert inputfile.stat().st_size<32*1024*1024
    exe=work/'verify-fixed-cover'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    with (out/'compile.log').open('w') as log:
        subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-started)))
    with inputfile.open() as f:
        proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,timeout=max(1,55-(time.monotonic()-started)))
    (out/'native.log').write_text(proc.stdout);(out/'native.stderr').write_text(proc.stderr)
    if proc.returncode:
        (out/'failure.json').write_text(json.dumps(dict(all_passed=False,reason=proc.stderr,returncode=proc.returncode,no_automatic_larger_retry=True),indent=2)+'\n');raise SystemExit(proc.returncode)
    result=json.loads(proc.stdout);assert result['all_passed'] and result['peak_with_reserve']<=150000
    proof=dict(all_passed=True,accepted_actual_G18_q5_source8_and_corrected_K2_source7_reused=True,
        old_square_and_adjugate_same_accepted_inputs=True,old_source_target_order_preserved=True,
        all56_new_RF_knownunit_LCD_clearings=True,all128_old_adj_products=True,
        exactly15_fixed_replacement8dots=True,all16_full_knownunit_primitive_factorizations=True,
        primitive_literal_Bezout_equals1=True,actual_minor_literal_Bezout_equals_E=True,
        E='z^21*(z+20)^93*(z+60)^6*(z+12)',
        all_polynomial_scalar_and_row_LCD_factors_are_units_on_entire_P100_chart=True,
        proof='The accepted actual pure source rows give the same fixed8 targets. Old and new LCD factors are whole-chart units. At every point the literal polynomial minor identity sums five fixed cleared8-minors to the unit E, so at least one is nonzero. Hence eight actual pure directions survive; together with the independently accepted141 block this proves original pointwise quartic rank at least149 on the entire permitted P100 chart.',
        whole_P100_original_quartic_pointwise_rank_lower149=True,
        no_fixed_original149minor_column_completion_full239_global_FULLU_or_exhaustion=True,
        no_gcd_rank_search_unknown_polynomial_inverse_sourcequery_point_prime_or_GB=True,checker=result)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report=dict(id=ID,all_passed=True,whole_P100_quartic_pointwise_lower149=True,
        no_full239_fixed149minor_or_global_FULLU_claim=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},ignored_native_input_sha256=sha(inputfile),
        ignored_executable_sha256=sha(exe),seconds=time.monotonic()-started,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
