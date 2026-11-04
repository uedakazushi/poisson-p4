#!/usr/bin/env python3
"""Seven zi*qB columns in the fixed5-unknown quartic scalar-offset system."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_ZI_QB_QUARTIC_EXTENSION_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'offset':'reports/AMBIENT_B18_T_QUADRATIC_L_SYZYGY_R4_20261104_v63/actual-quadratic-L-syzygy-scalar-offset-certificate.json',
      'rewrite':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
      'known40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
      'quad_peer':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json',
      'HF3peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_offsets_and_trueqB_before7columns','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};O,NF,K,QP,HP=[D[k] for k in paths];assert QP['all_passed'] and HP['all_passed']
    assert not O['nonzero_quartic_relation_witness'] and len(O['all_checked_source_constraints'])==1332
    linear=list(map(rec,O['normalized_actual37_linear_relations']));hnorm=list(map(F,O['normalized_actual37_e5_r2_values']))
    qb={(t,tuple(m)):F(c) for t,m,c in K['homogeneous_module_generators40'][38]}
    def ex(v):return tuple(v.count(i) for i in range(8))
    assert qb=={(3,ex([7,7])):F(1),(4,ex([5,7])):F(11767,3652)}
    mons2=list(itertools.combinations_with_replacement(range(8),2));mons3=[tuple(m) for m in NF['monomials3_120']];m3={m:i for i,m in enumerate(mons3)}
    standard=NF['standard_global600_positions'];nforms={i:rec(p) for i,s,p in NF['all600_rewrite_steps']}
    lifts={i:rec(p) for i,p,k in O['all600_chosen_source_lifts_and_scalar_offsets']};offset={i:F(k) for i,p,k in O['all600_chosen_source_lifts_and_scalar_offsets']}
    basis=list(map(rec,O['normalized_standard_basis_columns']));srcbasis=list(map(rec,O['normalized_actual1332_source_basis']));fbasis=list(map(F,O['normalized_scalar_functional_values']));at={min(v):j for j,v in enumerate(basis)}
    assert len(basis)==2
    qbcols=[]
    for k in range(7):
        v={}
        for (t,m),a in qb.items():
            inds=[i for i,n in enumerate(m) for q in range(n)]+[k];v[120*t+m3[tuple(sorted(inds))]]=a
        qbcols.append(v)
    def source_map(c):
        v={}
        for index,a in c.items():
            if index>=1332:addscaled(v,qbcols[index-1332],a);continue
            p,m=divmod(index,36)
            for row,b in linear[p].items():
                t,x=divmod(row,8);j=120*t+m3[tuple(sorted(mons2[m]+(x,)))];z=v.get(j,F(0))+a*b
                if z:v[j]=z
                else:v.pop(j,None)
        return v
    rows=[];witness=None
    for k,raw in enumerate(qbcols):
        nf={};lift={1332+k:F(1)};fv=F(0)
        for i,a in raw.items():addscaled(nf,nforms[i],a);addscaled(lift,lifts[i],-a);fv-=a*offset[i]
        assert source_map(lift)=={standard[j]:a for j,a in nf.items()}
        v=dict(nf);actual=dict(lift);value=fv;trace=[]
        while v and min(v) in at:
            j=at[min(v)];a=v[min(v)];addscaled(v,basis[j],-a);addscaled(actual,srcbasis[j],-a);value-=a*fbasis[j];trace.append([j,str(a)])
        item={'variable_z_index':k,'normal_form5':encode(nf),'source_lift1339':encode(lift),'scalar_offset':str(fv),'subtract':trace}
        if v:
            j=min(v);a=v[j];at[j]=len(basis);basis.append({i:c/a for i,c in v.items()});srcbasis.append({i:c/a for i,c in actual.items()});fbasis.append(value/a)
            item.update(pivot=j,divide=str(a))
        else:
            assert not source_map(actual);item.update(zero=True,residual_scalar=str(value))
            if value:
                check=sum((a*hnorm[i//36] for i,a in actual.items() if i<1332 and mons2[i%36]==(7,7)),F(0));assert check==value
                witness={'source1339_combination':encode(actual),'actual_pure_e5_r4_relation_value':str(value),
                  'cubic_initial_coefficients_exactly_zero':True,'source_new7_values_identically_zero_because_zi_non_r':True}
        rows.append(item)
        if witness:break
    proof={'actual_qB_initial_form':[[t,list(m),str(c)] for (t,m),c in qb.items()],
      'actual_Q_lift_exists_via_accepted_rank1_quadratic_source':True,
      'new_columns_source_index1332_plus_i_means_zi_times_Q':True,
      'Q_third_or_higher_values_not_used':'pure e5*r4 of zi*Q is zero fori0..6',
      'all_checked_new7constraints':rows,'nonzero_quartic_witness':witness,
      'limited_family_consistent_if_no_witness':not witness,
      'actual_HF4_if_proved':0 if witness else None,'actual_T_total_length_if_proved':11 if witness else None,
      'no_new_sourcejet_X3_X4_Q3_or_sixth':True}
    dump(out,'actual-seven-zi-qB-scalar-offset-certificate.json',proof)
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL','checked7':len(rows),'new_actual_quartic_relation':bool(witness)})
    report={'id':STEM,'all_passed':True,'new_zi_qB_columns_checked':len(rows),'nonzero_actual_quartic_relation_found':bool(witness),
      'actual_T_total_length':11 if witness else None,'no_actual_HF4_when_family_consistent':not witness,
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','new_zi_qB_columns_checked','nonzero_actual_quartic_relation_found','actual_T_total_length','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
