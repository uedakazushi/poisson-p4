#!/usr/bin/env python3
"""Source-only native/old first separation for118 cubic-kernel candidates."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
from ambient_B18_comparison_T_next3_source_preflight_20261104_v63 import receipt_combos
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_QUARTIC118_NATIVE_FIRST_GATE_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def put(v,i,a):
    c=v.get(i,F(0))+a
    if c:v[i]=c
    else:v.pop(i,None)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'cubicrank':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63/actual-linear-syzygy-cubic-source-values.json',
      'c1':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase2-linear-syzygy-source-root-contractions.json',
      'a2':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
      'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
      'native128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
      'Xp':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
      'N0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_first_arrays_before118_native_gate','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};R,C,A,F0,N,X,N0=[D[k] for k in paths]
    normalized,kernel,labels=receipt_combos(R['complete_fixed_rank_certificate'],119);assert len(kernel)==118
    c119=[[rec(v) for v in rows] for rows in C['all119_linear_coefficient_original768_combinations']]
    a119=list(map(rec,A['all119_actual_a2_sparse27_times36']))
    x2119=list(map(rec,A['all119_known_x2=sum_linear_l_times_Xj1_sparse56_times36']))
    raw=F0['source_columns640_rows56_sparse_dual9']+N['source_columns128_rows48_sparse_dual9'];assert len(raw)==768
    cols=[{i:rec(v) for i,v in c} for c in raw]
    xp=[rec(v) for v in X['normalized_Xp27_original768_combinations']]
    xp_native=[];xp_old=[]
    for c in xp:
        n={};o={}
        for j,a in c.items():
            for row,p in cols[j].items():
                b=p.get(0,F(0))
                if b:put(n if j>=640 else o,row,a*b)
        xp_native.append(n);xp_old.append(o)
    mons=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(mons)}
    cases=[];structural=0
    for number,c in enumerate(kernel):
        c1=[{} for k in range(8)];a2={};knownfull={}
        for j,a in c.items():
            for k in range(8):addscaled(c1[k],c119[j][k],a)
            addscaled(a2,a119[j],a);addscaled(knownfull,x2119[j],a)
        # Distinct original-source pieces, not the mixed normalizedXp label.
        ns=set();os=set()
        for k,cs in enumerate(c1):
            origin={}
            for j,a in cs.items():
                for row,p in cols[j].items():
                    if j>=640 and p.get(0):put(origin,row,a*p[0])
                    for d,b in p.items():
                        if d: (ns if j>=640 else os).add(36*row+mi[tuple(sorted((k,d-1)))])
            assert not origin,'native origin must vanish separately for each coefficient variable'
        for index,a in a2.items():
            p,m=divmod(index,36);ns.update(36*row+m for row in xp_native[p]);os.update(36*row+m for row in xp_old[p])
        structural+=len(ns)+len(os)
        cases.append((c1,a2,knownfull,ns,os))
    bound=structural+sum(len(v) for c1,a2,x,ns,os in cases for v in c1)+sum(len(a2) for c1,a2,x,ns,os in cases)+10000
    assert bound<150000,('SPARSE_SOURCE_BOUND',bound)
    dump(out,'preflight.json',{'source_kernel_count':118,'native_and_old_residual_structural_slots':structural,
      'TOTAL_sparse_QQ_bound':bound,'cap':150000,'dense118_by48_by36_forbidden':True,'no_new_source_coefficients':True})
    matrix=list(map(lambda row:[list(map(F,p)) for p in row],N0['N']));inverse=list(map(lambda row:[list(map(F,p)) for p in row],N0['N_inverse']))
    assert len(matrix)==len(inverse)==12
    def C4apply(v,mat):
        z={}
        for index,a in v.items():
            row,m=divmod(index,36);j,b=divmod(row,4)
            for i in range(12):
                for k,c in enumerate(mat[j][i]):
                    if c and not(b&k):put(z,36*(4*i+(b|k))+m,a*c)
        return z
    output=[];nonzero=0;support=0;r2only=[]
    for number,(c1,a2,knownfull,ns,os) in enumerate(cases):
        n={};o={};total={}
        for k,cs in enumerate(c1):
            for j,a in cs.items():
                for row,p in cols[j].items():
                    for d,b in p.items():
                        if d:
                            q=36*row+mi[tuple(sorted((k,d-1)))];put(n if j>=640 else o,q,a*b);put(total,q,a*b)
        assert total==knownfull
        for index,a in a2.items():
            p,m=divmod(index,36)
            for row,b in xp_native[p].items():put(n,36*row+m,-a*b)
            for row,b in xp_old[p].items():put(o,36*row+m,-a*b)
        full=dict(n);addscaled(full,o,F(1));assert not full
        assert set(n)<=ns and set(o)<=os and all(index//36<48 for index in n)
        proper=C4apply(n,inverse);assert C4apply(proper,matrix)==n
        if n:nonzero+=1
        support+=sum(map(len,c1))+len(a2)+len(n)+len(proper);assert support<150000
        rr={i:c for i,c in proper.items() if i%36==35}
        if rr:r2only.append(number)
        output.append({'kernel118_index':number,'fixed119_zero_receipt':labels[number],
          'combination_in119':encode(kernel[number]),'linear_original768_combinations_by8variables':list(map(encode,c1)),
          'known_normalized_Xp27_quadratic_a2':encode(a2),
          'native_only_X1_minus_native_Xp0_a2_residual48_times36':encode(n),
          'proper_S1_c1_minus_S0_XpNative_a2':encode(proper),
          'old_only_residual_is_exact_negative':True,'full56_residual_exact_zero':True})
    dump(out,'actual118-old-native-first-separation.json',{'monomials36':mons,'all118':output,
      'all_native_origin_coefficient_vectors_separately_zero':True,
      'nonzero_native_separation_count':nonzero,'nonzero_proper_pure_r2_candidate_indices':r2only,
      'N0_inverse_applied_in_true_C4_tauA2_tauB2_zero_ring':True,
      'unadjusted_L2_sensitivity_gate_only':True,
      'lower_cubic_corrections_and_final_adjoint_may_still_cancel_residual':True,
      'N2_needed_or_native_L2_cancels_not_claimed':True,'no_new_source_N2_X3_coefficient':True,'expanded_saved_QQ_support':support})
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL','nonzero_native_residuals':nonzero,'pure_r2_residuals':len(r2only)})
    report={'id':STEM,'all_passed':True,'actual118_unadjusted_native_separation_nonzero_count':nonzero,
      'actual118_unadjusted_native_pure_r2_nonzero_count':len(r2only),
      'new_N2_or_X3_values_computed':False,'final_N2_dependence_or_HF4_not_decided':True,
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual118_unadjusted_native_separation_nonzero_count','actual118_unadjusted_native_pure_r2_nonzero_count','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
