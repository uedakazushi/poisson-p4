#!/usr/bin/env python3
"""Five-unknown scalar-offset test, not a dense1332-column kernel solve."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
from ambient_B18_comparison_T_next3_source_preflight_20261104_v63 import receipt_combos
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_QUADRATIC_L_SYZYGY_R4_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'rewrite':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'pure':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63/actual-linear-syzygies-and-pure-r3-source-values.json',
      'HF3peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
      'known40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_known_arrays_before_scalar_offset_test','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};NF,L,L1,H,HP,K40=[D[k] for k in paths]
    assert HP['all_passed'] and sha(ROOT/paths['HF3peer'])=='808637882f3d8d28f393a7b0f4abfeec59e076bd1269034bd67d8ac6dd8c706a'
    linear=list(map(rec,NF['linear_basis37']));assert len(linear)==37
    normalize,zeros,labels=receipt_combos(L['complete_linear_rank_certificate'],657)
    hsource={r['minimal_cycle_index657']:F(r['actual_e5_r2_coefficient']) for r in H['actual37_e5_r2_source_coefficients']}
    assert len(hsource)==37 and all(set(c)<=set(hsource) for c in normalize)
    hnorm=[sum((a*hsource[j] for j,a in c.items()),F(0)) for c in normalize]
    mons2=list(itertools.combinations_with_replacement(range(8),2));m2={m:i for i,m in enumerate(mons2)}
    mons3=[tuple(m) for m in NF['monomials3_120']];m3={m:i for i,m in enumerate(mons3)}
    steps={i:(step,rec(v)) for i,step,v in NF['all600_rewrite_steps']};standard=NF['standard_global600_positions'];assert len(standard)==5
    dep={}
    for i in reversed(range(600)):
        s,nf=steps[i]
        if s.get('standard'):dep[i]=set()
        else:
            t,m=divmod(i,120);rest=list(mons3[m]);rest.remove(s['leading_variable']);p=s['linear_relation_basis_index']
            dep[i]={36*p+m2[tuple(rest)]}
            for child,c in s['children']:
                assert child>i;dep[i].update(dep[child])
    bound=sum(map(len,dep.values()))+5*1332+1332*7+600*5+5000;assert bound<150000
    dump(out,'preflight.json',{'source_map_dimensions':[600,1332],'standard_unknowns5':standard,
      'structural_source_lift_support_upper':sum(map(len,dep.values())),'TOTAL_QQ_bound':bound,'cap':150000,
      'no_full_kernel_matrix_or_source_jet':True})
    lifts={};offset={}
    for i in reversed(range(600)):
        s,nf=steps[i]
        if s.get('standard'):c={}
        else:
            t,m=divmod(i,120);rest=list(mons3[m]);rest.remove(s['leading_variable']);p=s['linear_relation_basis_index'];c={36*p+m2[tuple(rest)]:F(1)}
            for child,a in s['children']:addscaled(c,lifts[child],F(a))
        assert set(c)<=dep[i];lifts[i]=c;offset[i]=sum((a*hnorm[j//36] for j,a in c.items() if mons2[j%36]==(7,7)),F(0))
    def source_map(c):
        v={}
        for index,a in c.items():
            p,m=divmod(index,36)
            for row,b in linear[p].items():
                t,x=divmod(row,8);j=120*t+m3[tuple(sorted(mons2[m]+(x,)))];z=v.get(j,F(0))+a*b
                if z:v[j]=z
                else:v.pop(j,None)
        return v
    for i,c in lifts.items():
        expected={i:F(1)}
        for j,a in steps[i][1].items():
            k=standard[j];v=expected.get(k,F(0))-a
            if v:expected[k]=v
            else:expected.pop(k,None)
        assert source_map(c)==expected
    basis=[];atrow={};sourcebasis=[];fbasis=[];rows=[];witness=None
    for p in range(37):
        for m,mon in enumerate(mons2):
            index=36*p+m;raw=source_map({index:F(1)});nf={};sub={};k=F(0)
            for i,a in raw.items():addscaled(nf,steps[i][1],a);addscaled(sub,lifts[i],a);k+=a*offset[i]
            lifted={index:F(1)};addscaled(lifted,sub,F(-1));fv=(hnorm[p] if mon==(7,7) else F(0))-k
            assert source_map(lifted)=={standard[i]:a for i,a in nf.items()}
            v=dict(nf);actual=dict(lifted);value=fv;trace=[]
            while v and min(v) in atrow:
                j=atrow[min(v)];c=v[min(v)];addscaled(v,basis[j],-c);addscaled(actual,sourcebasis[j],-c);value-=c*fbasis[j];trace.append([j,str(c)])
            row={'source1332_column':index,'normal_form5':encode(nf),'scalar_offset_functional':str(fv),'subtract':trace}
            if v:
                i=min(v);c=v[i];atrow[i]=len(basis);basis.append({j:a/c for j,a in v.items()});sourcebasis.append({j:a/c for j,a in actual.items()});fbasis.append(value/c)
                row.update(pivot=i,divide=str(c),normalized=encode(basis[-1]),normalized_scalar=str(value/c))
            else:
                row.update(zero=True,residual_scalar=str(value))
                assert not source_map(actual)
                if value:
                    witness={'source1332_inconsistent_constraint':index,'actual_normalized37_quadratic_coefficients':encode(actual),
                      'actual_pure_e5_r4_relation_value':str(value),'quadratic_coefficient_linear_syzygy_exact':True}
            rows.append(row)
            if witness:break
        if witness:break
    if witness:
        c=rec(witness['actual_normalized37_quadratic_coefficients']);value=sum((a*hnorm[j//36] for j,a in c.items() if mons2[j%36]==(7,7)),F(0));assert value==F(witness['actual_pure_e5_r4_relation_value'])
        translated={}
        for index,a in c.items():
            p,m=divmod(index,36)
            for source,b in normalize[p].items():
                j=36*source+m;v=translated.get(j,F(0))+a*b
                if v:translated[j]=v
                else:translated.pop(j,None)
        original_linear=list(map(rec,L1['actual40_by657_linear_relation_columns']));check={}
        for index,a in translated.items():
            source,m=divmod(index,36)
            for row,b in original_linear[source].items():
                t,v=divmod(row,8);j=120*t+m3[tuple(sorted(mons2[m]+(v,)))];z=check.get(j,F(0))+a*b
                if z:check[j]=z
                else:check.pop(j,None)
        assert not check;witness['actual_original657_quadratic_coefficients']=encode(translated)
    ev={'normalized_actual37_linear_relations':list(map(encode,linear)),
      'normalized_actual37_source_combinations_in657':list(map(encode,normalize)),
      'normalized_actual37_e5_r2_values':list(map(str,hnorm)),
      'all600_chosen_source_lifts_and_scalar_offsets':[[i,encode(lifts[i]),str(offset[i])] for i in range(600)],
      'all_checked_source_constraints':rows,'small_standard_rank':len(basis),
      'normalized_standard_basis_columns':list(map(encode,basis)),
      'normalized_actual1332_source_basis':list(map(encode,sourcebasis)),
      'normalized_scalar_functional_values':list(map(str,fbasis)),
      'nonzero_quartic_relation_witness':witness,
      'this_family_negative_if_no_witness':witness is None,
      'actual_HF4_computed':bool(witness),'actual_HF4_if_proved':0 if witness else None,
      'actual_T_total_length_if_proved':11 if witness else None,
      'no_source_X3_X4_sixth_jet_or_new_point':True}
    dump(out,'actual-quadratic-L-syzygy-scalar-offset-certificate.json',ev)
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL','checked_constraints':len(rows),'new_quartic_relation':bool(witness)})
    report={'id':STEM,'all_passed':True,'small_standard_unknowns':5,'checked_actual_source_constraints':len(rows),
      'actual_nonzero_pure_quartic_relation_found':bool(witness),'actual_T_total_length':11 if witness else None,
      'no_actual_HF4_claim_when_consistent':witness is None,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','checked_actual_source_constraints','actual_nonzero_pure_quartic_relation_found','actual_T_total_length','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
