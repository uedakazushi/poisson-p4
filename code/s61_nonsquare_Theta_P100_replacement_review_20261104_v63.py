#!/usr/bin/env python3
"""Independent dense-integer polynomial audit, with no solve or new rank test."""
import argparse,hashlib,json,math,resource,signal,sys,time
from fractions import Fraction as F
from pathlib import Path
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1];STEM='S61_NONSQUARE_THETA_P100_REPLACEMENT_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def trim(p):
    while p and p[-1]==0:p.pop()
    return p
def decode(p):
    q=[0]*(1+max([int(i) for i,c in p]+[-1]));seen=set()
    for i,c in p:
        i=int(i);assert i not in seen;seen.add(i);v=F(c);assert v.denominator==1;q[i]=v.numerator
    return trim(q)
def mul(a,b):
    if not a or not b:return []
    p=[0]*(len(a)+len(b)-1);aa=[(i,c) for i,c in enumerate(a) if c];bb=[(j,d) for j,d in enumerate(b) if d]
    for i,c in aa:
        for j,d in bb:p[i+j]+=c*d
    return trim(p)
def power(p,n):
    q=[1]
    while n:
        if n&1:q=mul(q,p)
        n//=2
        if n:p=mul(p,p)
    return q
def scale(a,c):return trim([x*c for x in a])
def sum_products(pairs):
    total=[]
    for a,b in pairs:
        p=mul(a,b)
        if len(total)<len(p):total +=[0]*(len(p)-len(total))
        for i,c in enumerate(p):total[i]+=c
    return trim(total)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    pb=ROOT/'reports/S61_NONSQUARE_THETA_P100_REPLACEMENT_PREFLIGHT_20261104_v63';sb=ROOT/'reports/S61_NONSQUARE_THETA_P100_REPLACEMENT_SOLVE_20261104_v63';gb=ROOT/'reports/S61_NONSQUARE_THETA_P100_CUBIC_GENERIC_20261104_v63';ub=pb
    inputs=[]
    for folder in [pb,sb,gb,ub]:
        r=read(folder/'report.json');assert r['frozen']
        for name,h in r['input_sha256'].items():assert sha(ROOT/name)==h,name
        for name,h in r['evidence_sha256'].items():assert sha(folder/name)==h,name
        inputs.append(folder/'report.json')
    prior=gb.with_name(gb.name+'_REVIEW')/'report.json';assert read(prior)['all_passed'];inputs.append(prior)
    for p,h in read(prior)['input_sha256'].items():assert sha(ROOT/p)==h,p
    c=read(pb/'certificate.json');n=read(sb/'native-certificate.json');u=read(ub/'exact-known-unit-strip.json');old=read(gb/'native-rank.json');source=read(gb/'substituted17x24-before-rank.json')
    original=[[decode(p) for p in row] for row in source['integer_matrix17x24']];M=[[decode(p) for p in row] for row in c['prepared_integer_matrix17x24']]
    assert len(M)==17 and all(len(row)==24 for row in M)
    linears=[[0,1],[20,1],[60,1],[2700,203],[12,1]];assert [decode(p) for p in c['five_known_chart_linear_units']]==linears
    rf=[decode(p) for p in c['full_expanded_row_units']];cf=[decode(p) for p in c['full_expanded_column_units']]
    for i in range(17):
        expected=[1]
        for p,k in zip(linears,c['row_knownunit_exponents'][i]):assert isinstance(k,int) and k>=0;expected=mul(expected,power(p,k))
        scalar=F(c['row_QQ_units'][i]);assert scalar and scalar.denominator==1;assert scale(expected,scalar.numerator)==rf[i]
    for j in range(24):
        expected=[1]
        for p,k in zip(linears,c['column_knownunit_exponents'][j]):assert isinstance(k,int) and k>=0;expected=mul(expected,power(p,k))
        assert expected==cf[j]
    for i in range(17):
        for j in range(24):assert mul(mul(rf[i],cf[j]),M[i][j])==original[i][j],(i,j)
    rows=c['source_row_order_for_A_B'];cols=c['A_source_columns'];extras=c['B_source_columns']
    assert rows==old['selected_local_rows'] and set(rows)==set(range(17));assert cols==old['selected_local_columns']==list(range(17)) and extras==list(range(17,24))
    A=[[M[i][j] for j in cols] for i in rows];B=[[M[i][j] for j in extras] for i in rows]
    assert A==[[decode(p) for p in row] for row in c['prepared_A17x17']] and B==[[decode(p) for p in row] for row in c['prepared_B17x7']]
    unitdet=[1]
    for i in rows:unitdet=mul(unitdet,rf[i])
    for j in cols:unitdet=mul(unitdet,cf[j])
    expected=decode(c['expected_prepared_det_from_saved_source_NO_newdet']);originaldet=decode(old['nonzero_selected_minor'])
    assert unitdet==decode(c['whole_det_unit_multiplier']) and mul(unitdet,expected)==originaldet and len(expected)-1==134
    removed=[sum(c['row_knownunit_exponents'][i][k] for i in rows)+sum(c['column_knownunit_exponents'][j][k] for j in cols) for k in range(5)]
    assert removed==c['removed_row_column_det_unit_powers']==[242,37,34,0,0]
    current=decode(u['remaining_primitive_polynomial']);assert len(current)-1==46
    assert decode(u['source_selected_minor'])==originaldet
    allpowers=[];whole=scale(current,int(u['QQ_content']))
    for p,step in zip(linears,u['all5_named_unit_power_divisions']):
        assert p==decode(step['unit']);k=step['power'];whole=mul(whole,power(p,k));allpowers.append(k)
    assert whole==originaldet
    rest=[old-new for old,new in zip(allpowers,removed)];assert rest==c['remaining_known_unit_det_powers']==[54,18,7,9,0]
    reconstructed=scale(current,int(u['QQ_content']))
    for p,k in zip(linears,rest):reconstructed=mul(reconstructed,power(p,k))
    assert reconstructed==expected
    content=math.gcd(*current);R=[x//content for x in current];assert scale(R,content)==current
    ratio=F(n['solver_denominator_over_expected_QQunit']);assert ratio and ratio.denominator==1;delta=decode(n['solver_denominator']);assert delta==scale(expected,ratio.numerator)
    selected=n['selected_Bcolumn'];which=n['selected_row'];assert selected==which==0 and n['Bezout_saved']
    X=[decode(p) for p in n['selected17_vector']];assert len(X)==17
    eqstats=[]
    for i in range(17):
        lhs=sum_products([(A[i][j],X[j]) for j in range(17)]);rhs=mul(delta,B[i][selected]);assert lhs==rhs,i
        eqstats.append({'row':i,'exact_residual_zero':True,'product_degree':len(lhs)-1})
    S=decode(n['Bezout_S_integer']);T=decode(n['Bezout_T_integer']);constant=int(n['Bezout_nonzero_integer_constant']);assert constant
    assert sum_products([(S,X[which]),(T,R)])==[constant]
    evidence={'independent_engine':'dense univariate integer coefficient arrays; no source-solver imports',
      'accepted_original_P100_source_and_minor_receipt':str(prior.relative_to(ROOT)),
      'all408_entry_unit_transports_checked':True,'all_row_and_column_unit_power_formulas_checked':True,
      'old447minor_all5_known_power_full_reconstruction_checked':True,'removed_known_powers':removed,'remaining_known_powers':rest,
      'prepared134det_is_actual_old_minor_divided_only_known_units':True,'solver_delta_to_actual_det_QQ_ratio':str(ratio),
      'selected_vector17_linear_identities':eqstats,'selected_replacement_row':which,'selected_extra_column':selected,
      'integer_Bezout_all_coefficients_zero_except_constant':True,'nonzero_Bezout_constant_bits':abs(constant).bit_length(),
      'R46_irreducibility_not_needed_for_coverage_proof':True,
      'two_actual_source_minors_cover_every_geometric_P100_point':True,'ALL_P100_full_cubic57_accepted':True,
      'all119_vector_dump_or_AX_claim_not_independently_replayed':True,
      'no_new_solve_ranksearch_point_prime_GB_P100_or_quartic239_conclusion':True}
    ep=out/'independent-unit-transport-AX-and-Bezout.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    inputs +=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pb/'certificate.json',sb/'native-certificate.json',gb/'substituted17x24-before-rank.json',gb/'native-rank.json',ub/'exact-known-unit-strip.json']
    report={'id':STEM,'all_passed':True,'frozen':True,'independent_actual408_transports_144det_17AX_and_Bezout':True,'ALL_P100_full_cubic57_proved':True,'no_author_solve_rerun_or_other_rank_query':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
