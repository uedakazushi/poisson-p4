#!/usr/bin/env python3
"""Whole619 quartic source-order adjoint; no new scalar value evaluation.

The high-order linear part is folded exactly. Lower X2 products are retained
as typed, source-bound operators, not replaced by incomplete B11 values.
"""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
from ambient_B18_comparison_T_next3_source_preflight_20261104_v63 import receipt_combos
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63'
def dump(out,name,x):
    (out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
def put(v,i,a):
    c=v.get(i,F(0))+a
    if c:v[i]=c
    else:v.pop(i,None)
def enc(v):return [[list(i) if isinstance(i,tuple) else i,str(a)] for i,a in sorted(v.items()) if a]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'source620':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
      'rank1':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
      'gamma3':'reports/AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63/phase1-Gamma-and-b3-source.json',
      'b3':'reports/AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63/actual-b3-and-all118-pure-quartic-values.json',
      'gammaC':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase1-cubic-Pi2-and-source-adjoint.json',
      'c4':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_VALUES_20261104_v63/actual-c4-and168-corrected-d3-quartic-values.json',
      'c4_peer':'reports/AMBIENT_B18_T_CORRECTED168_QUARTIC_VALUES_20261104_v63_REVIEW/report.json',
      'gates':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'selected37':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'all640first':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
      'all128first':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
      'native2':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'known619':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
      'known619_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
      'minimal':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
      'source_bridge':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'source_manifests_persisted_before619_source_adjoint','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};S,R,GA,B3,GC,C4,CP,G,B,H,A,L,X1,N1,N2,K,KP,MIN=[D[k] for k in list(paths)[:18]]
    assert D['source_bridge']['all_passed']
    assert sha(ROOT/paths['source_bridge'])=='6b608c0249efc99d52514285bd388112dc0a4e640075d91c9af389c7e12bf80a'
    assert CP['all_passed'] and KP['all_passed']
    assert sha(ROOT/paths['c4_peer'])=='88e314495185407ed08ea52cd493d557528bb3f74275f4e2134a2eb5897d6a27'
    assert all(not F(a) for a in K['all619_actual_pure_cubic_relation_values'])
    mons2=list(itertools.combinations_with_replacement(range(8),2));mons3=list(itertools.combinations_with_replacement(range(8),3))
    b3=F(B3['canonical_qB_pure_cubic_b3']);c4=F(C4['canonical_qC_pure_quartic_c4'])
    gamma=rec(GA['Gamma0_56_times120']);addscaled(gamma,rec(GA['Gammab_56_times120']),b3);addscaled(gamma,rec(GC['two_Xsource3_adjoint_rows56_times120'][0]),c4)
    ops=B['operator_scalar_nodes'];bc=list(map(dict,B['B48_by640_sparse_columns']));nodes=H['DAG_scalar_nodes']
    def xkey(j,src,m,order):
        if j>=640:return (2,j-640,src,order,m) if src<48 else None
        if src<48:
            root=bc[j].get(src)
            if root is None or ops[root][0]=='c':return None
            n=ops[root];return (1,n[1],order,m) if n[0]=='actual_C0_scalar_function' else (0,root,order,m)
        root=H['honest_D0_relation_columns'][j][src-48+1]
        return None if nodes[root][0]=='c' else (1,root,order,m)
    def rawP(c):
        v={}
        for j,a in c.items():
            if j<640:
                root=H['honest_D0_relation_columns'][j][8]
                if nodes[root][0]!='c':put(v,root,a)
        return v
    norm,kernel,labels=receipt_combos(R['complete_fixed_rank_certificate'],620);assert len(kernel)==619 and len(norm)==1
    c620=[rec(r['combination_in_actual768']) for r in S['all620']];a620=[rec(r['linear_Xp27_correction_row8j_plus_variable']) for r in S['all620']]
    xp=list(map(rec,S['normalized_Xp27_original768_combinations']))
    high=[];rootlabels=set();rawroots=set();xproots=set();support=len(gamma);a1oldr=0
    for i,combo in enumerate(kernel):
        c={};a1={}
        for j,a in combo.items():addscaled(c,c620[j],a);addscaled(a1,a620[j],a)
        p4=rawP(c);xp3={};xf3={};xf2={}
        for p in range(27):
            r=a1.get(8*p+7,F(0))
            if r:
                for root,a in G['normalized_Xp27_source55_rr'][p]['rawP8_root_coefficients']:put(xp3,root,-r*F(a))
        if xp3:a1oldr+=1
        # High-order contribution: -Gamma4(Xc3-Xp2*a1).
        # These are ACTUAL native coefficients; no effective-S2 alias is made.
        for index,g in gamma.items():
            src,m=divmod(index,120);mon=mons3[m]
            for j,a in c.items():
                label=xkey(j,src,m,3)
                if label is not None:put(xf3,label,-g*a)
            for v in set(mon):
                pair=list(mon);pair.remove(v);m2=mons2.index(tuple(pair))
                for p in range(27):
                    a=a1.get(8*p+v,F(0))
                    if a:
                        for j,b in xp[p].items():
                            label=xkey(j,src,m2,2)
                            if label is not None:put(xf2,label,g*a*b)
        rootlabels.update(xf3);rootlabels.update(xf2);rawroots.update(p4);xproots.update(xp3)
        support+=len(combo)+len(c)+len(a1)+len(p4)+len(xp3)+len(xf3)+len(xf2)
        assert support<150000,('HIGH_ADJOINT_SUPPORT_GUARD',support)
        high.append({'kernel619_index':i,'fixed620_zero_receipt':labels[i],'combination_in620':encode(combo),
          'actual_original768_constant_combination':encode(c),'actual_Xp27_linear_a1':encode(a1),
          'rawP55_pure_r4_roots':encode(p4),'Xp55_pure_r3_times_a1r_roots':encode(xp3),
          'actual_X3_high_order_source_functional':enc(xf3),'actual_Xp2_times_a1_source_functional':enc(xf2),
          'lower_source_operators_not_dropped':['b1','l1','a2','linear_oldP3_corrections','Gamma4_Xstar2_b1','Gamma4_XL2_l1','Gamma4_Xp1_a2','Xp55_2_a2rr']})
    orders={o:sorted({tuple((mons3 if o==3 else mons2)[label[-1]]) for label in rootlabels if label[-2]==o}) for o in [2,3]}
    native3={label[1] for label in rootlabels if label[0]==2 and label[-2]==3};native2={label[1] for label in rootlabels if label[0]==2 and label[-2]==2}
    dump(out,'all619-high-order-root-adjoints.json',{'Gamma4_56_times120':encode(gamma),'monomials2_36':mons2,'monomials3_120':mons3,
      'root_kind0_typedB_kind1_C0_kind2_ACTUAL_native':True,'all619':high,'all_values_still_symbolic':True})
    typed={
      'input_source':'c0 from each actual619 constant620-kernel combination',
      'a1':'frozen actual X1 cancellation Xc1=Xp0*a1',
      'C2':'Xc2-Xp1*a1 (full actual second coefficients)',
      'b1':'U0^-1*kappa_degree2(m0*C2)',
      'q2':'pi*C2-Ystar1*b1 in known177 span (from actual619 kernel)',
      'l1':'fixed first296 elimination receipt lift of q2 through selectedL37*variables',
      'a2':'fixed27 pivot lift of C2-Xstar1*b1-XL1*l1 (full56 equality required)',
      'source_z':'Xc-Xp*a1-Xstar*b1-XL*l1-Xp*a2, starts degree3',
      'C3':'Xc3-Xp2*a1-Xstar2*b1-XL2*l1-Xp1*a2',
      'beta2':'U0^-1*kappa_degree3(m0*C3)',
      'Qproxy3':'pi*C3-Ystar1*beta2',
      'Theta4':'Theta3+c4*PiC, b3=-42215/3456, c4=4565/68832',
      'corrected_pure_r4':'rawP4+linear_oldP3-Xp55_2*a2rr-Ystar55_2*beta2rr-Theta4(Qproxy3)',
      'expanded_linear_oldP3':'-Xp55_3*a1r-Xstar55_3*b1r-XL55_3*l1r',
      'm1_scope':'z starts3, so(mz)3=m0*C3. No claim that individual lower coefficient functions are cached.',
      'Euler_boundary':'difference d3*h with h in m has puredegree4 zero, using accepted corrected168 family',
      'native_source3':'L0*S3+L1*S2+L2*S1+L3*S0',
      'native_actual_source2':'L0*S2+L1*S1+L2*S0',
      'native_L3_scope':'S0*c0=0 separately; L3*S0*c0 cancels exactly',
      'native_L2_remaining':'L2*(S1*c0+S0*c1), c1=-Xp*a1-Xstar*b1-XL*l1; not assumed zero',
      'effective_S2_warning':'DO NOT replace native Xp2, Xstar2 or XL2 separately by effectiveS2; their L2*S0 terms participate in the remaining sensitivity',
      'lower_operator_computations_executed':False,'new_scalar_function_values_or_source_jets_computed':False}
    dump(out,'whole619-typed-source-and-order-profile.json',typed)
    profile={'actual619_sources':619,'high_order_folded_QQ_support':support,'high_order_scalar_labels':len(rootlabels),
      'high_order_X3_monomials':orders[3],'lower_Xp2_a1_monomials':orders[2],
      'rawP4_root_count':len(rawroots),'Xp3_root_count':len(xproots),'sources_with_explicit_oldg3_factor':a1oldr,
      'native_X3_source_columns':len(native3),'native_ACTUAL_X2_source_columns':len(native2),
      'new_highest_source_gates':['rawP4 common-source/pullback fourth sensitivity','oldXp3 honest quotient numerator6 if sensitivity survives',
        'X3 zi*r² source quotient derivatives if sensitivity survives','actual native L2 sensitivity after full c1 corrections',
        'lower b1/l1/a2 second-pair sensitivity before B11 substitution'],
      'source6_is_sufficient_or_necessary_proved':False,'no_actual619_quartic_values_or_HF4_computed':True,
      'no_new_source6_N2_common4_value_solve':True,'cap':150000,'high_order_only_total_bound':support+len(rootlabels)+10000}
    assert profile['high_order_only_total_bound']<150000
    dump(out,'preflight.json',profile);dump(out,'phase.json',{'phase':'FINAL_SOURCE_ONLY','high_support':support})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'source_only_not_actual_values':True,'actual619_sources':619,
      'high_order_source_support':support,'whole_lower_coefficient_expansion_NOT_performed':True,
      'new_source_coefficient_jet_or_solve':False,'actual_HF4_or_total_length':None,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump(out,'report.json',report);print(json.dumps(profile))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
