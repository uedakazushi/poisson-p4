#!/usr/bin/env python3
"""Independent Fraction proof check of exact strip and two Theta curves."""
import argparse,collections,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
from s61_nonsquare_p_branch_chain_review_20261104_v63 import ROOT,sha,read,dec,add,mul,sc,pw,shift,substitute,monic_reduce
STEM='S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW'
def enc(p):return [[list(reversed(m)),str(c)] for m,c in sorted(p.items())]
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    names=['THETA_D14_UNIT_STRIP','THETA_RANK_DROP_CURVES','THETA_DELTA_LAMBDA_PROFILES','THETA_REMAINING_DEFECTS','THETA_RATIONAL_COORDINATES','RAW_H_THETA_UNIT_PROFILE','RAW_H_GENERIC_RANK']
    inputs=[];cs=[]
    for n in names:
        directory=ROOT/'reports'/f'S61_NONSQUARE_{n}_20261104_v63';r=read(directory/'report.json');cs.append(read(directory/'certificate.json'))
        assert sha(directory/'certificate.json')==r['certificate_sha256']
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
        inputs.extend([directory/'report.json',directory/'certificate.json'])
    strip,curve,d,remaining,chart,t,g=cs
    finaldir=ROOT/'reports/S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63';operatorpath=finaldir/'operator-recipe.json';operator=read(operatorpath)
    assert sha(operatorpath)==read(finaldir/'report.json')['operator_recipe_sha256'];inputs.append(operatorpath)
    for n in ['RAW_H_P_BRANCH_CHAIN','THETA_FIXED_DEFECTS']:
        p=ROOT/'reports'/f'S61_NONSQUARE_{n}_20261104_v63_REVIEW/report.json';r=read(p);assert r['all_passed']
        for q,h in r['input_sha256'].items():assert sha(ROOT/q)==h,q
        inputs.append(p)
    O={(0,0):F(1)};D={(1,0):F(1)};z={(0,1):F(1)};zp=add(z,sc(O,20));L=add(pw(z,2,2),sc(D,50))
    localizers={'zeta':z,'zeta_plus20':zp,'L':L}
    for key in ['K','Fnum','Jnum']:localizers[key]=dec(chart['target_coordinates'][key])
    assert {k:dec(v) for k,v in strip['only_certified_unit_factors']}==localizers
    actualDelta=next(dec(x['polynomial']) for x in d['actual_images'] if x['name']=='Delta')
    actualD14=next(dec(x['defect']) for x in remaining['two_fixed_defects'] if x['raw_row']==14)
    expected_sources={'Delta':actualDelta,'D14':actualD14};results={}
    for result in strip['native_exact_strip']['results']:
        name=result['name'];source=dec(result['source']);assert source==expected_sources[name]
        content=F(result['QQ_content']);off=result['d_unit_exponent'];assert content and off==min(m[0] for m in source)
        current=shift(source,(-off,0),1/content);counts=collections.Counter()
        for step in result['steps']:
            nameunit=step['localizer'];factor=dec(step['gcd']);q=dec(step['source_quotient']);u=dec(step['localizer_quotient'])
            assert factor and mul(factor,q)==current and mul(factor,u)==localizers[nameunit]
            # The original localizer is a unit, so each factor is a unit.
            assert factor==localizers[nameunit];current=q;counts[nameunit]+=1
        rest=dec(result['remaining_unknown']);assert current==rest
        results[name]=(content,off,counts,rest)
    assert results['Delta'][1:3]==(-7,collections.Counter({'zeta':8,'zeta_plus20':7,'L':6}))
    assert results['Delta'][3]==sc(O,-1)
    assert results['D14'][1:3]==(-16,collections.Counter({'zeta':15,'zeta_plus20':14,'L':13}))
    P={m:add(sc(mul(z,zp),3),sc(D,m)) for m in [80,100]}
    assert results['D14'][3]==sc(mul(P[80],P[100]),-1)
    assert add(P[100],P[80],F(-1))==sc(D,20)
    assert dec(curve['P80'])==P[80] and dec(curve['P100'])==P[100]
    reconstructed={}
    for name,(content,off,counts,rest) in results.items():
        q=rest
        for key,n in counts.items():q=mul(q,pw(localizers[key],n,2))
        q=shift(q,(off,0),content);assert q==expected_sources[name];reconstructed[name]=True
    # Preserve the actual determinant ordering, without recomputing huge minors.
    assert dec(strip['literal_original_Delta'])==dec(g['nonzero_selected_minor_Delta'])
    bm=strip['bordered_minor_binding'];assert bm['row_order_raw']==[7,13,10,14] and bm['column_order']==[2,4,3,'B4_cleared_lastSchur']
    assert bm['original_raw_heavy5_row_order']==[8,7,13,10,14] and bm['original_raw_heavy5_column_order']==[0,2,4,3,5]
    rawdir=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';rawpeer=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW/report.json'
    rr=read(rawpeer);assert rr['all_passed'] and rr['independent_frame_Cdelta_eta_gamma_and_full_source_clearings_passed'] and rr['all90_final_checkpoint_coefficients_and_fixed_row_operations_match']
    pre=read(rawdir/'preflight.json');checkpoint=read(rawdir/'checkpoint.json');rawreport=read(rawdir/'report.json')
    for name in ['preflight.json','checkpoint.json']:assert sha(rawdir/name)==rawreport['evidence_sha256'][name]
    H=[[F(x) for x in row] for row in pre['constant_coefficient_column_frame7x7']];Hi=[[F(x) for x in row] for row in pre['constant_inverse7x7']]
    for a,b in [(H,Hi),(Hi,H)]:assert all(sum(a[i][k]*b[k][j] for k in range(7))==int(i==j) for i in range(7) for j in range(7))
    step=checkpoint['steps'][0];assert len(checkpoint['steps'])==1 and (step['pivot_row'],step['pivot_column'])==(1,1)
    specs=checkpoint['initial_row_divisors']+[step['pivot_unit'],step['row_scale_before']]+[q['row_divisor_after'] for q in step['eliminations']]
    assert specs==operator['first_two_raw_H_units']['all_unit_specifications']
    for spec in specs:assert F(spec['constant']) and len(spec['a_d_F_B_C_J_powers'])==6 and all(isinstance(x,int) for x in spec['a_d_F_B_C_J_powers'])
    assert step['pivot_unit']=={'constant':'-1/14062500','a_d_F_B_C_J_powers':[0,0,0,1,1,0]}
    assert operator['B_unit3']['selected_raw_rows']==[7,13,10] and operator['B_unit3']['selected_middle_column_order']==[0,2,1]
    assert operator['last_scalar']['zero_raw_rows']==[9,11,12] and operator['last_scalar']['nonzero_raw_row']==14
    inputs +=[rawpeer,rawdir/'preflight.json',rawdir/'checkpoint.json',rawdir/'report.json']
    # Universal bordered determinant: det[P,s;b,t]=det(P)*t-b*adj(P)*s.
    # Previously checked det(P)=Delta, Delta*b=lambda*P and the integral
    # source ring give b*adj(P)=lambda. This is polynomial equality before
    # localization. The actual raw8 unit pivot contributes f8/B^4.
    profiles=[]
    def lift_unary(p):return {(0,m[0]):F(c) for m,c in p if F(c)}
    source_profiles={key:dec(chart['target_coordinates'][key]) for key in ['a','e','B','C']}
    source_profiles.update(d=D,F=shift(localizers['Fnum'],(-2,0),F(1,40)),J=shift(localizers['Jnum'],(-1,0),F(1,100)))
    for m in [80,100]:
        dimage=sc(mul(z,zp),F(-3,m))
        assert not substitute(P[m],[dimage,z],2)
        mapped={key:substitute(value,[dimage,z],2) for key,value in localizers.items()}
        mapped['d']=dimage
        assert all(value and all(e[0]==0 for e in value) for value in mapped.values())
        product=O
        for value in mapped.values():product=mul(product,value)
        assert product
        saved=next(p for p in curve['curve_profiles'] if p['m']==m)
        assert lift_unary(saved['d_graph'])==dimage
        assert all(lift_unary(p)==mapped[key] for key,p in saved['L_K_Fnum_Jnum_exact_profiles'].items())
        for key,p in source_profiles.items():
            n=max([-ex[0] for ex in p]+[0]);poly=shift(p,(n,0));num=substitute(poly,[dimage,z],2);den=pw(dimage,n,2)
            sp=saved['original_a_e_d_B_C_F_J_rationals'][key];sn,sd=lift_unary(sp['numerator']),lift_unary(sp['denominator'])
            assert num and den and sn and sd and mul(num,sd)==mul(sn,den)
        open_factors=[lift_unary(p) for p in saved['exact_univariate_open_factors']];op=O;seen=[0]*len(open_factors);factorization=[]
        for factor in open_factors:assert factor and max(ex[1] for ex in factor)==1;op=mul(op,factor)
        assert op==lift_unary(saved['nonzero_open_product']) and op
        for key,value in mapped.items():
            rest=value;counts=[]
            for j,factor in enumerate(open_factors):
                lead=factor[(0,1)];monic=sc(factor,1/lead);count=0
                while rest and max(ex[1] for ex in rest)>0:
                    q,r=monic_reduce(rest,monic,1,1)
                    if r:break
                    rest=q;count+=1;seen[j]+=1
                counts.append(count)
            assert rest and set(rest)=={(0,0)};factorization.append({'localizer':key,'linear_multiplicities':counts,'nonzero_QQ_remainder':str(rest[(0,0)])})
        assert all(seen)
        # Every allowed localizer remains a nonzero polynomial in z, so
        # quotient is a nonempty localization of QQ[z], hence a domain.
        profiles.append({'m':m,'d_graph':enc(dimage),'localizer_images':{k:enc(v) for k,v in mapped.items()},'localizer_product':enc(product),'exact_open_factorizations':factorization,'all7_original_source_coordinate_rationals_checked':True})
    # Final primary-specific source profiles and operator recipe are bound
    # below after the independent fixed-form calculation.
    ep=out/'independent-strip-and-two-curves.json'
    evidence={'all63_saved_step_source_and_localizer_quotients_exact':True,
      'global_Delta_unit_known_factors_only':True,'D14_unknown_remainder_exact_negative_P80_P100':True,
      'all_source_reconstructions':reconstructed,'two_curve_profiles':profiles,
      'P100_minus_P80_is_20d_unit':True,'bordered_minor_orientation_and_raw8_scale_accepted':True,
      'actual_initial_two_units_frame49_both_inverses_and_all_specifications_rebound':True,
      'Theta_transpose_heavy_and_original_rawH_cokernel':'R/(P80*P100)',
      'rank_original_raw_H_on_two_curves':6,'rank_original_raw_H_off_two_curves':7,
      'nonempty_disjoint_reduced_rational_curves':True,
      'uniform_rawH_rank7_false_but_full_cubic57_or_quartic239_not_decided':True}
    ep.write_text(json.dumps(evidence,indent=2)+'\n')
    inputs +=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']
    report={'id':STEM,'all_passed':True,'independent_Fraction_all_saved_strip_divisions_and_source_reconstruction':True,
      'exact_two_curve_localizations_and_actual_cokernel_proof':True,'no_global_unknown_inverse_point_prime_ranksearch_GB_or_full_cubic_negative':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
