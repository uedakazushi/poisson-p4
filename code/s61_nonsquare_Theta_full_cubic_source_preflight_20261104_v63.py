#!/usr/bin/env python3
"""Read-only exact source/bit/support preflight. No coefficient substitution."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_p_branch_chain_review_20261104_v63 import ROOT,sha,read
STEM='S61_NONSQUARE_THETA_FULL_CUBIC_SOURCE_PREFLIGHT_20261104_v63'
def empty():return {'entries':0,'nonzero_entries':0,'terms':0,'max_entry_terms':0,'max_core_total_degree':0,'max_a_e_d_degrees':[0,0,0],'max_QQ_numerator_bits':0,'max_QQ_denominator_bits':0,'minimum_possible_z_exponent':0,'maximum_possible_z_exponent':0,'max_individual_z_cleared_degree':0,'curve_integer_coefficient_bit_bounds':{'80':0,'100':0}}
def record(stats,poly):
    stats['entries']+=1;stats['terms']+=len(poly);stats['max_entry_terms']=max(stats['max_entry_terms'],len(poly))
    if not poly:return
    stats['nonzero_entries']+=1;amax=[max(m[j] for m,c in poly) for j in range(3)];common=1
    low=min(d-2*a-e for (a,e,d),c in poly);high=max(2*d-a for (a,e,d),c in poly);clear=max(0,-low)
    stats['minimum_possible_z_exponent']=min(stats['minimum_possible_z_exponent'],low);stats['maximum_possible_z_exponent']=max(stats['maximum_possible_z_exponent'],high)
    stats['max_individual_z_cleared_degree']=max(stats['max_individual_z_cleared_degree'],high+clear)
    for j in range(3):stats['max_a_e_d_degrees'][j]=max(stats['max_a_e_d_degrees'][j],amax[j])
    for m,c in poly:
        stats['max_core_total_degree']=max(stats['max_core_total_degree'],sum(m));stats['max_QQ_numerator_bits']=max(stats['max_QQ_numerator_bits'],abs(c.numerator).bit_length());stats['max_QQ_denominator_bits']=max(stats['max_QQ_denominator_bits'],c.denominator.bit_length());common=math.lcm(common,c.denominator)
    for curve,den,lnorm in [(80,(45,9,80),(696,82,63)),(100,(9,9,100),(210,101,63))]:
        factors=[[den[j]**(amax[j]-k)*lnorm[j]**k for k in range(amax[j]+1)] for j in range(3)]
        bound=sum(abs(c.numerator)*(common//c.denominator)*factors[0][m[0]]*factors[1][m[1]]*factors[2][m[2]] for m,c in poly)
        stats['curve_integer_coefficient_bit_bounds'][str(curve)]=max(stats['curve_integer_coefficient_bit_bounds'][str(curve)],bound.bit_length())
def polyjson(p):return [(tuple(m),F(c)) for m,c in p if F(c)]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    raw=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt'
    sourcep=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';bp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';basis=read(bp)
    assert sha(raw)==basis['remaining_matrix_sha256'];source=read(sourcep);mons=source['normal_monomials'];heavy=[j for j,m in enumerate(mons) if m[5]];assert len(heavy)==7 and len(mons)==91
    profiles={n:empty() for n in ['full57x91','product41x84','genuine16x91','genuine16_pure84','genuine16_heavy7']}
    with raw.open() as f:
        assert f.readline().split()==['57','91','57']
        for i in range(57):
            for j in range(91):
                line=f.readline().split();count=int(line[0]);assert len(line)==1+4*count
                p=[(tuple(map(int,line[2+4*k:5+4*k])),F(line[1+4*k])) for k in range(count)];assert all(c and min(m)>=0 for m,c in p)
                record(profiles['full57x91'],p)
                if i<41:
                    if j in heavy:assert not p
                    else:record(profiles['product41x84'],p)
                else:
                    record(profiles['genuine16x91'],p);record(profiles['genuine16_heavy7' if j in heavy else 'genuine16_pure84'],p)
        assert not f.read().strip()
    augdir=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63';aug=read(augdir/'input.json');ar=read(augdir/'report.json');assert sha(augdir/'input.json')==ar['input_json_sha256']
    profiles['projected_augmented25x32']=empty()
    for row in aug['matrix25x32']:
        for p in row:record(profiles['projected_augmented25x32'],polyjson(p))
    cpdir=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63';cp=read(cpdir/'checkpoint.json');cr=read(cpdir/'report.json');assert sha(cpdir/'checkpoint.json')==cr['evidence_sha256']['checkpoint.json'] and cp['completed_pivots']==8
    rr,cc=cp['remaining_rows'],cp['remaining_columns'];assert(len(rr),len(cc))==(17,24)
    for n in ['residual17x24_numerators','residual17x24_all_denominators','residual17x24_nonzero_entry_denominators']:profiles[n]=empty()
    for i in rr:
        for j in cc:
            p,q=cp['matrix25x32'][i][j];assert q
            record(profiles['residual17x24_numerators'],polyjson(p));record(profiles['residual17x24_all_denominators'],polyjson(q))
            if p:record(profiles['residual17x24_nonzero_entry_denominators'],polyjson(q))
    assert profiles['full57x91']['terms']==164333 and profiles['projected_augmented25x32']['terms']==60380 and profiles['residual17x24_numerators']['terms']==84253
    gates=[ROOT/'reports'/n/'report.json' for n in ['S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW','S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63','S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW','S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW']]
    assert all(read(p)['all_passed'] for p in gates)
    criterion={
      'source':str(raw.relative_to(ROOT)),'source_rows':57,'target_columns':91,'pure_target_columns':84,'heavy_target_columns':7,
      'curve_generic_field':'QQ(zeta), separately P80 and P100; no numerical parameter evaluation',
      'whole_curve_known_ranks':{'product_image':41,'genuine_heavy_image':6,'combined':47},
      'exact_full_residual_shape':[10,43],
      'construction':'Select a nonzero pure41 product minor over QQ(zeta). Use six genuine rows with independent actual heavy6 coefficients. The resulting47 pivot rows/columns are independent. Every remaining row has zero heavy class since heavy image rank is exactly6; the10 residual rows are pure modulo product41, in rank43. Each entry can instead be specified by its actual48x48 bordered minor divided by the47pivot determinant.',
      'identity':'rank(full57x91)=47+rank(residual10x43)',
      'full57_iff':'residual row rank10',
      'fixed_old_DeltaP_or_Delta34_not_assumed_nonzero_on_either_curve':True,
      'cheaper_sufficient_test_shape':[17,24],
      'cheaper_identity':'rank(projected augmented25x32)=8+rank(this17x24)',
      'cheaper_full57_gate':'rank17 => full57; any smaller rank is inconclusive for full source because discarded pure directions can rescue',
      'smallest_formal_projected_residual_after_known_heavy4_and_quad3':[10,16],
      'smallest_formal_projected_residual_warning':'Not materialized without the old degree207 completion; not an exact full57 criterion.',
      'no_curve_coefficient_substitution_Gaussian_minor_or_POWER_coefficient_evaluation_performed':True}
    ep=out/'source-costs-and-exact-curve-criterion.json';ep.write_text(json.dumps({'profiles':profiles,'criterion':criterion,'curve_support_bound':'For a^i e^j d^k the Laurent zeta support lies between k-2i-j and2k-i. All curve coefficient bit values are conservative cleared-entry bounds, not substituted coefficients.','integer_bound_clear':'Per entry use source coefficient denominator LCM times45^imax*9^jmax*80^kmax (P80), or9^imax*9^jmax*100^kmax (P100), plus enough zeta power. The coefficient l1 bound is computed without polynomial expansion.'},indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',raw,sourcep,bp,augdir/'report.json',augdir/'input.json',cpdir/'report.json',cpdir/'checkpoint.json',*gates]
    report={'id':STEM,'all_passed':True,'actual_input_terms_degrees_and_rational_bits_counted':True,'exact10x43_full_curve_criterion_and17x24_sufficient_gate':True,'no_rank_substitution_new_point_prime_GB_or_power207_evaluation':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'all_passed':True,'seconds':report['seconds'],'max_rss_kib':report['max_rss_kib'],'profiles':profiles},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
