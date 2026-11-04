#!/usr/bin/env python3
"""Independent fixed source substitution and selected polynomial determinant."""
import argparse,hashlib,json,resource,signal,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_THETA_P100_CUBIC_GENERIC_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);work=ROOT/'work/P100-cubic-independent-review';work.mkdir(parents=True,exist_ok=True);hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def gate(p):
        r=read(p);assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        return r
    reports=ROOT/'reports';base=reports/'S61_NONSQUARE_THETA_P100_CUBIC_GENERIC_20261104_v63';primary=gate(base/'report.json')
    for stem in ['S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW','S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW','S61_NONSQUARE_THETA_FULL_CUBIC_SOURCE_PREFLIGHT_20261104_v63']:
        gate(reports/stem/'report.json')
    cp=read(reports/'S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json')
    binding=read(base/'source-denominator-and-operator-binding.json');matrix=read(base/'substituted17x24-before-rank.json');rank=read(base/'native-rank.json')
    rows=cp['remaining_rows'];cols=cp['remaining_columns'];assert len(rows)==17 and len(cols)==24 and cp['completed_pivots']==8
    assert rows==binding['row_ids25'] and cols==binding['target_ids32']
    B={(1,0,1):F(50),(0,2,1):F(18),(0,1,1):F(-9),(0,1,0):F(60),(0,0,1):F(-2),(0,0,0):F(-40)}
    # Independently bind the three P100 rational functions to the accepted Theta chart.
    import sympy as sy
    zz=sy.symbols('zeta')
    aa=-10*(zz+20)/(9*zz**2);ee=(zz+100)/(9*zz);dd=-3*zz*(zz+20)/100
    assert sy.cancel(ee-sy.Rational(2,3)-(zz**2-400)/(60*dd))==0
    assert sy.cancel(500*dd**2*aa-(zz+20)*((zz+20)*(zz-10)+50*dd))==0
    assert sy.cancel(100*dd+3*zz*(zz+20))==0
    entries=[];classes=[];rawterms=0;exponents=set()
    for i in rows:
        for j in cols:
            p,q=cp['matrix25x32'][i][j];pd={tuple(m):F(v)for m,v in p};qd={tuple(m):F(v)for m,v in q};assert all(v for v in pd.values())and all(v for v in qd.values())
            if set(qd)=={(0,0,0)}:kind=0;c=qd[(0,0,0)]
            elif set(qd)=={(1,0,0)}:kind=1;c=qd[(1,0,0)]
            else:
                kind=2;c=qd[(1,0,1)]/50;assert qd=={m:c*v for m,v in B.items()}
            assert c and c.denominator==1 and all(v.denominator==1 for v in pd.values())
            entries.append((kind,c,pd));classes.append(dict(source_row25=i,target_column32=j,known_unit_denominator_kind=['1','a','B'][kind],nonzero_QQ_scale=str(c)))
            rawterms+=len(pd);exponents.update(pd)
    assert classes==binding['all408_denominator_classifications'] and rawterms==84253 and len(exponents)==1113
    assert primary['all_P100_points_full57_proved'] is False
    assert sorted(rank['selected_local_rows'])==list(range(17)) and len(set(rank['selected_local_columns']))==17 and rank['exact_QQ_zeta_rank']==17
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_Theta_P80_cubic_generic_review_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']
    for p in paths:hashes[str(p.relative_to(ROOT))]=sha(p)
    inp=work/'fixed-source-and-minor.txt'
    with inp.open('w')as f:
        def uwrite(p):
            p={int(e):F(v)for e,v in p}if isinstance(p,list)else p
            f.write(str(len(p))+''.join(f' {e} {v}'for e,v in sorted(p.items()))+'\n')
        def swrite(p):f.write(str(len(p))+''.join(f' {a} {e} {d} {v}'for (a,e,d),v in sorted(p.items()))+'\n')
        # Own direct QQ Laurent numerators: a=an/zeta²,e=en/zeta,d=zeta*dn.
        uwrite({0:F(-200,9),1:F(-10,9)});uwrite({0:F(100,9),1:F(1,9)});uwrite({0:F(-3,5),1:F(-3,100)})
        swrite(B);uwrite({1:F(5),2:F(1,12)})
        for kind,c,p in entries:f.write(f'{kind} {c}\n');swrite(p)
        for i in range(17):
            f.write(f"{matrix['row_LCMs'][i]} {matrix['row_integer_contents'][i]}\n")
            for p in matrix['integer_matrix17x24'][i]:uwrite(p)
        f.write(' '.join(map(str,rank['selected_local_rows']))+'\n');f.write(' '.join(map(str,rank['selected_local_columns']))+'\n')
        uwrite(rank['nonzero_selected_minor']);uwrite(rank['fflu_denominator'])
    assert inp.stat().st_size<16*1024*1024
    exe=work/'verify-native'
    with(out/'compile.log').open('w')as log:
        subprocess.run(['c++','-O2','-std=c++17',str(paths[1]),'-o',str(exe),*subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
    with inp.open()as f:proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,check=True,timeout=max(1,55-(time.monotonic()-start)))
    (out/'native.log').write_text(proc.stdout);assert not proc.stderr;result=json.loads(proc.stdout);assert result['all_passed']
    proof=dict(accepted_actual_checkpoint_unit8_source_gate=True,all408_literal_denominator_classifications=True,
        actual_P100_parametrization_bound_to_accepted_curve_review=True,common_scalar_zeta24_a_B_unit_on_permitted_open=True,
        all_source_entries_and_row_constant_clearings_independently_checked=True,
        fixed_selected17_minor_recomputed_with_every_division_remultiplied=True,
        exact_projected_remaining_rank=17,augmented_rank=25,full_cubic_generic_rank=57,
        full_cubic_possible_defect_is_finite_on_permitted_P100=True,
        whole_P100_uniform_cubic57_not_proved=True,no_other_curve_quartic239_or_exhaustion_claim=True,
        no_POWER_cofactor_expansion_newpoint_prime_ranksearch_author_solver_or_GB=True,checker=result)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},
        ignored_verifier_input_sha256=sha(inp),ignored_executable_sha256=sha(exe),
        full_cubic57_generic_P100_independently_proved=True,uniform_whole_curve_cubic57_proved=False,
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
