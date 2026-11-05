#!/usr/bin/env python3
"""Independent fixed nine cache projections, after SOURCE acceptance only."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_TAUA_CACHED9_VALUES_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_TAUA_CACHED9_VALUES_20261105_v63'
SOURCE='AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63'
PRIMARY='0bc5fcb9d6e4343a4436192a298d79cd2843e82281564c7da3637ab0acfb65a3'
SRCPEER='b7cd179090ff8a22e20de1fa007ffbe3d54f25bd7168bd374dd6d14f8bbf594b'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        pp=ROOT/'reports'/PRI/'report.json';assert sha(pp)==PRIMARY;pr=read(pp);assert pr['all_passed']
        srp=ROOT/'reports'/(SOURCE+'_REVIEW')/'report.json';assert sha(srp)==SRCPEER;sr=read(srp);assert sr['all_passed']
        pins={str(pp.relative_to(ROOT)):PRIMARY,str(srp.relative_to(ROOT)):SRCPEER}
        for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        for p,h in pr['evidence_sha256'].items():
            ep=ROOT/'reports'/PRI/p;assert sha(ep)==h;pins[str(ep.relative_to(ROOT))]=h
        assert pr['code_sha256']=='5891d81746e4eb9ad9df0e5d93879c0d8246cadfc496e03a17e79ccf42249888'
        ev=read(ROOT/'reports'/PRI/'all9-cached-cubic-dots-and-termwise-source-products.json')
        D=read(ROOT/'reports'/SOURCE/'nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json')
        c2=read(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
        K=[[Q(c) for c in r] for r in D['complete_constant_kernel12_by9']];free=D['free_columns9']
        # The SOURCE proves this literal simplification; do not pair any unselected pivot leaf.
        assert all(K[i][j]==int(i==free[j]) for i in range(12) for j in range(9))
        ratio=Q(D['gA_reduction_ratio']);assert ratio==Q(8,25)
        actual=[];height=0
        for j,i in enumerate(free):
            row=c2['physical12_coefficients'][i];v=dict(row['delta']);p=Q(v.get(29,'0'));q=Q(v.get(30,'0'))
            for x in [p,q,ratio]:assert max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=18
            pb=max(abs(ratio.numerator).bit_length()+abs(p.numerator).bit_length(),ratio.denominator.bit_length()+p.denominator.bit_length());assert pb<=2048
            # Reverse the prepared adjoint term order and use literal source indices29/30.
            prod=-ratio*p
            ab=max(max(abs(prod.numerator).bit_length()+q.denominator.bit_length(),abs(q.numerator).bit_length()+prod.denominator.bit_length())+1,prod.denominator.bit_length()+q.denominator.bit_length());assert ab<=2048
            value=prod+q;height=max(height,pb,ab)
            old=ev['all9_actual_cached_eA_cubic_dots'][j]
            assert old['source_index']==j and old['free_physical_row']==row['physical_row']
            expected=[[i,30,'1',str(q),str(q)],[i,29,'-8/25',str(p),str(prod)]]
            assert old['all_term_products']==expected and Q(old['value'])==value==0
            actual.append({'source_index':j,'physical_row':row['physical_row'],'delta_r3':str(p),'delta_r2ell':str(q),
              'reversed_weighted_r3_term':str(prod),'independent_value':str(value)})
        assert ev['nonzero_source_indices']==[] and not pr['actual_C3_length_upper19_proved']
        result={'all9_own_source_values':actual,'all18_saved_term_products_verified':True,
          'complete_kernel_literal_free_embedding_checked_before_simplification':True,
          'all9_exact_zero':True,'actual_C3_length_upper19_NOT_inferred':True,
          'scope':'Only the nine cached tauA eA-r2ell projections are zero; other cubic coordinates and quartic corrections remain open.',
          'maximum_operational_height_bound_used':height,'hard_bits':2048,
          'no_author_producer_imports_or_new_source_gauge_jet_solve':True,
          'static_numeric_upper':18000,'hard_numeric_cap':30000,'hard_RAW_bytes':400000}
        save(out/'independent-nine-cache-values-and-term-checks.json',result)
        raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in out.iterdir())
        assert raw<400000 and time.monotonic()-start<60
        report={'id':ID,'all_passed':True,'frozen':True,'primary_report_sha256':PRIMARY,'source_peer_sha256':SRCPEER,
          'input_sha256':pins,'code_sha256':sha(__file__),
          'evidence_sha256':{'independent-nine-cache-values-and-term-checks.json':sha(out/'independent-nine-cache-values-and-term-checks.json')},
          'all9_exact_zero_verified':True,'partial_negative_only':True,
          'no_C3_length_or_exponent_improvement_claim':True,'RAW_input_plus_output_bytes':raw,'elapsed_seconds':time.monotonic()-start}
        save(out/'report.json',report);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
    except Exception as e:
        save(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_positive_claim':True});raise
if __name__=='__main__':main()
