#!/usr/bin/env python3
"""ONE authorized P80 symbolic substitution+generic17x24 sufficient test."""
import argparse,json,os,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_p_branch_chain_review_20261104_v63 import ROOT,read,sha
STEM='S61_NONSQUARE_THETA_P80_CUBIC_GENERIC_20261104_v63'
def main(out):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    def timeout(*_):raise TimeoutError('one combined60s P80 source+rank cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    cpdir=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63';cp=read(cpdir/'checkpoint.json');cr=read(cpdir/'report.json');assert sha(cpdir/'checkpoint.json')==cr['evidence_sha256']['checkpoint.json']
    gates=[ROOT/'reports'/n/'report.json' for n in ['S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW','S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW','S61_NONSQUARE_THETA_FULL_CUBIC_SOURCE_PREFLIGHT_20261104_v63']];assert all(read(p)['all_passed'] for p in gates)
    rows,cols=cp['remaining_rows'],cp['remaining_columns'];assert len(rows)==17 and len(cols)==24 and cp['completed_pivots']==8
    B={(1,0,1):F(50),(0,2,1):F(18),(0,1,1):F(-9),(0,1,0):F(60),(0,0,1):F(-2),(0,0,0):F(-40)}
    classes=[];lines=['17 24'];inputterms=0
    for i in rows:
        for j in cols:
            p,q=cp['matrix25x32'][i][j];qd={tuple(m):F(c) for m,c in q if F(c)}
            if len(qd)==1 and (0,0,0) in qd:kind=0;scale=qd[(0,0,0)]
            elif len(qd)==1 and (1,0,0) in qd:kind=1;scale=qd[(1,0,0)]
            else:
                kind=2;assert set(qd)==set(B);scale=qd[(1,0,1)]/50;assert qd=={m:scale*c for m,c in B.items()}
            assert scale and scale.denominator==1
            classes.append({'source_row25':i,'target_column32':j,'known_unit_denominator_kind':['1','a','B'][kind],'nonzero_QQ_scale':str(scale)})
            lines.append(f'{kind} {scale.numerator} {len(p)}');inputterms+=len(p)
            for m,c in p:
                coefficient=F(c);assert coefficient.denominator==1;lines.append(str(coefficient.numerator)+' '+' '.join(map(str,m)))
    assert inputterms==84253
    # Before evaluation, source support guarantees degree<=34, hence at most
    # 408*35=14280 ordinary coefficient slots, below the150k field-term guard.
    assert 408*35<150000
    ip=out/'source17x24-input.txt';ip.write_text('\n'.join(lines)+'\n');bp=out/'source-denominator-and-operator-binding.json';bp.write_text(json.dumps({'curve':'P80','a':'-8*(7*zeta+80)/(45*zeta^2)','e':'2*(zeta+40)/(9*zeta)','d':'-3*zeta*(zeta+20)/80','row_ids25':rows,'target_ids32':cols,'known_unit_common_scalar':'zeta^24*a*B','all408_denominator_classifications':classes,'expanded_field_term_upper':14280,'raw_source_numerator_terms':inputterms,'eight_known_unit_checkpoint_equivalence_persists_on_P80':True,'projected_rank17_suffices_full57_even_when_full_heavy_rank6':True},separators=(',',':'))+'\n')
    cpp=Path(__file__).with_suffix('.cpp');exe=ROOT/'work/Theta-P80-cubic-generic-native';prefix=Path('/workspace/cloud-setup/prefix/usr');lib=prefix/'lib/x86_64-linux-gnu';env=os.environ.copy();env['LD_LIBRARY_PATH']=str(lib)+':'+env.get('LD_LIBRARY_PATH','')
    cmd=['g++','-O2','-std=c++17',str(cpp),'-I'+str(prefix/'include'),'-L'+str(lib),'-lflint','-lmpfr','-lgmp','-o',str(exe)]
    subprocess.run(cmd,check=True,timeout=min(20,59-time.monotonic()+started),env=env)
    checkpoint=out/'substituted17x24-before-rank.json';native=out/'native-rank.json';log=out/'native.log'
    failure=None
    try:
        with log.open('w') as f:subprocess.run([str(exe),str(ip),str(checkpoint),str(native)],check=True,timeout=59-time.monotonic()+started,stdout=f,stderr=subprocess.STDOUT,env=env)
    except (subprocess.CalledProcessError,subprocess.TimeoutExpired,TimeoutError) as e:failure=str(e)
    total=sum(p.stat().st_size for p in out.iterdir() if p.is_file());assert total<32*1024*1024
    inputs=[Path(__file__),cpp,ROOT/'research/notes'/f'{STEM}.md',cpdir/'report.json',cpdir/'checkpoint.json',*gates]
    report={'id':STEM,'all_passed':failure is None,'frozen':True,'independent_review_pending':True,'only_P80_generic_projected_test':True,'source_numerator_terms':inputterms,'input_field_coefficient_slot_bound':14280,'no_POWER207_P100_point_prime_or_full57_expansion':True,'failure_or_cap':failure,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'evidence_sha256':{p.name:sha(p) for p in [ip,bp,checkpoint,native,log] if p.exists()},'output_bytes':total,'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    if failure is None:
        result=read(native);data=read(checkpoint);assert data['ordinary_polynomial_terms']<=150000 and data['max_degree']<=34
        r=result['exact_QQ_zeta_rank'];assert r in range(18);report.update(exact_projected_remaining_rank=r,exact_augmented25x32_rank=8+r,full_cubic57_generic_P80_proved=r==17,all_P80_points_full57_proved=False,smaller_projected_rank_does_not_decide_full_cubic=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
