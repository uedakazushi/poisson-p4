#!/usr/bin/env python3
"""P100: exact named-unit strip/factor and row-column replacement-solve input ONLY."""
import hashlib, importlib.util, json, math, resource, signal, subprocess, time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_THETA_P100_REPLACEMENT_PREFLIGHT_20261104_v63'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text())

def main(out):
    start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(exist_ok=False,parents=True)
    hp=ROOT/'code/s61_nonsquare_Theta_P80_replacement_preflight_20261104_v63.py'
    spec=importlib.util.spec_from_file_location('P80_polynomial_helpers_only',hp)
    ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
    poly,ser,mul,div,ppow=ar.poly,ar.ser,ar.multiply,ar.divide,ar.powpoly
    base=ROOT/'reports/S61_NONSQUARE_THETA_P100_CUBIC_GENERIC_20261104_v63'
    report=load(base/'report.json');ip=base/'substituted17x24-before-rank.json';np=base/'native-rank.json'
    assert sha(ip)==report['evidence_sha256'][ip.name] and sha(np)==report['evidence_sha256'][np.name]
    inp=load(ip);raw=[[poly(p)for p in row]for row in inp['integer_matrix17x24']]
    assert len(raw)==17 and all(len(row)==24 for row in raw)
    nr=load(np);rows=nr['selected_local_rows'];cols=nr['selected_local_columns'];extras=[j for j in range(24)if j not in cols]
    assert cols==list(range(17)) and extras==list(range(17,24)) and set(rows)==set(range(17))
    olddet=poly(nr['nonzero_selected_minor']);assert max(olddet)==447
    cb=ROOT/'reports/S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63';cp=cb/'certificate.json'
    assert sha(cp)==load(cb/'report.json')['certificate_sha256']
    curve=next(x for x in load(cp)['curve_profiles']if x['m']==100)
    units=[poly([[int(m[0]),v]for m,v in p])for p in curve['exact_univariate_open_factors']]
    assert [ser(p)for p in units]==[[[1,'1']],[[0,'20'],[1,'1']],[[0,'60'],[1,'1']],[[0,'2700'],[1,'203']],[[0,'12'],[1,'1']]]
    residual=dict(olddet);strip=[];product={0:Q(1)}
    for L in units:
        power=ar.linval(residual,L);factor=ppow(L,power);residual,rem=div(residual,factor);assert not rem
        strip.append({'unit':ser(L),'power':power});product=mul(product,factor)
    den=math.lcm(*(c.denominator for c in residual.values()))
    content=math.gcd(*(int(c*den)for c in residual.values()))
    qcontent=Q(content,den);primitive={m:c/qcontent for m,c in residual.items()}
    assert mul(product,{m:qcontent*c for m,c in primitive.items()})==olddet
    rbits=max(abs(c.numerator).bit_length()for c in primitive.values())
    assert max(primitive)<=100 and rbits<=2000
    stripcert={'source_selected_minor':ser(olddet),'source_minor_degree':447,'all5_named_unit_power_divisions':strip,
      'QQ_content':str(qcontent),'remaining_primitive_polynomial':ser(primitive),'remaining_degree':max(primitive),
      'remaining_integer_coefficient_bits':rbits,'whole_source_minor_reconstruction_checked':True}
    (out/'exact-known-unit-strip.json').write_text(json.dumps(stripcert,indent=2)+'\n')
    # Only this residual factorization, with literal expansion checked; no rank or solve.
    fi=out/'factor-input.txt'
    fi.write_text(str(len(primitive))+''.join(' '+str(c)+' '+str(m)+' 0 0'for m,c in sorted(primitive.items(),reverse=True))+'\n')
    cpp=ROOT/'code/s61_nonsquare_raw_H_B_common_curve_factor_20261104_v63.cpp';exe=ROOT/'work/P100-preflight-factor'
    with (out/'compile.log').open('w')as log:
        subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    exesha=sha(exe)
    with(out/'factor.log').open('w')as log:
        subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','P100-preflight',str(exe),str(fi),str(out/'factor-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    factors=load(out/'factor-certificate.json');expanded={0:Q(factors['factor_constant'])};inventory=[]
    for rec in factors['factors']:
        assert all(m[1:]==[0,0]for m,c in rec['factor'])
        p=poly([[m[0],c]for m,c in rec['factor']])
        for _ in range(rec['exponent']):expanded=mul(expanded,p)
        inventory.append({'degree':max(p),'exponent':rec['exponent'],'factor':ser(p)})
    assert expanded==primitive
    M=[[dict(p)for p in row]for row in raw];rf=[{0:Q(1)}for _ in range(17)];cf=[{0:Q(1)}for _ in range(24)]
    rexp=[[0]*5 for _ in range(17)];cexp=[[0]*5 for _ in range(24)]
    for i in range(17):
        for k,L in enumerate(units):
            v=min(ar.linval(p,L)for p in M[i]);assert v<10**9
            factor=ppow(L,v);rf[i]=mul(rf[i],factor);rexp[i][k]=v
            for j in range(24):M[i][j],rem=div(M[i][j],factor);assert not rem
    for j in range(24):
        for k,L in enumerate(units):
            v=min(ar.linval(M[i][j],L)for i in range(17));assert v<10**9
            factor=ppow(L,v);cf[j]=mul(cf[j],factor);cexp[j][k]=v
            for i in range(17):M[i][j],rem=div(M[i][j],factor);assert not rem
    rowqq=[]
    for i in range(17):
        den=math.lcm(*(v.denominator for p in M[i]for v in p.values()))
        content=math.gcd(*(int(v*den)for p in M[i]for v in p.values()));assert content>0
        scale=Q(content,den);rowqq.append(scale);rf[i]={m:v*scale for m,v in rf[i].items()}
        M[i]=[{m:v/scale for m,v in p.items()}for p in M[i]]
    for i in range(17):
        for j in range(24):assert mul(mul(rf[i],cf[j]),M[i][j])==raw[i][j]
    unitdet={0:Q(1)}
    for i in rows:unitdet=mul(unitdet,rf[i])
    for j in cols:unitdet=mul(unitdet,cf[j])
    expected,rem=div(olddet,unitdet);assert not rem and all(v.denominator==1 for v in expected.values())
    assert mul(unitdet,expected)==olddet
    removed=[sum(rexp[i][k]for i in rows)+sum(cexp[j][k]for j in cols)for k in range(5)]
    oldpowers=[x['power']for x in strip];remaining=[a-b for a,b in zip(oldpowers,removed)];assert min(remaining)>=0
    rcost=ar.costs(raw,rows,cols,extras);pcost=ar.costs(M,rows,cols,extras)
    assert pcost['matrix_terms']<150000 and pcost['all119_replacement_coeff_slot_upper']<150000
    cert={'variable':'zeta','whole_chart':'P100 intersect FULLU only; no compactification points',
      'source_matrix17x24':str(ip.relative_to(ROOT)),'source_row_order_for_A_B':rows,'A_source_columns':cols,'B_source_columns':extras,
      'five_known_chart_linear_units':[ser(p)for p in units],'row_knownunit_exponents':rexp,'column_knownunit_exponents':cexp,
      'row_QQ_units':list(map(str,rowqq)),'full_expanded_row_units':[ser(p)for p in rf],'full_expanded_column_units':[ser(p)for p in cf],
      'all408_source_entry_reconstruction_checked':True,'prepared_integer_matrix17x24':[[ser(p)for p in row]for row in M],
      'prepared_A17x17':[[ser(M[i][j])for j in cols]for i in rows],'prepared_B17x7':[[ser(M[i][j])for j in extras]for i in rows],
      'expected_prepared_det_from_saved_source_NO_newdet':ser(expected),'whole_det_unit_multiplier':ser(unitdet),
      'source_det_equals_unitdet_times_expected_checked':True,'old_known_unit_det_powers':oldpowers,
      'removed_row_column_det_unit_powers':removed,'remaining_known_unit_det_powers':remaining,
      'residual_irreducible_factor_inventory':inventory,'whole_residual_factor_product_independently_checked':True,
      'raw_cost':rcost,'prepared_cost':pcost,
      'future_gate':'A X=delta B with delta/expected a nonzero QQ unit; replacement cofactors need coprime cover of EVERY residual factor, plus explicit source-vector and Bezout checks. No fullcurve promotion from this preflight.',
      'no_solve_rank_rightinverse_GB_point_prime_jobs':True}
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    exe.unlink(missing_ok=True)
    paths=[Path(__file__),hp,cpp,ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip,np,
      ROOT/'reports/S61_NONSQUARE_THETA_P100_CUBIC_GENERIC_20261104_v63_REVIEW/report.json',cb/'report.json',cp]
    r={'id':ID,'all_passed':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
      'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_factor_executable_SHA_not_committed':exesha,
      'selected_minor_knownunit_powers':oldpowers,'residual_degree':max(primitive),'residual_coefficient_bits':rbits,
      'residual_factor_degrees_exponents':[[p['degree'],p['exponent']]for p in inventory],
      'prepared_expected_det_degree':max(expected),'raw_cost':rcost,'prepared_cost':pcost,
      'new_solve_rank_GB_point_prime_jobs':0,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    signal.alarm(0);(out/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
