#!/usr/bin/env python3
"""Source-only full-c1 nativeL2 sensitivity of619 quartic channels."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_L2_ADJOINT_PREFLIGHT_20261104_v63'
def dump(o,n,v):
    (o/n).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in o.iterdir())<32*1024*1024
def put(v,k,a):
    z=v.get(k,F(0))+a
    if z:v[k]=z
    else:v.pop(k,None)
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def main(out):
    st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
      'high_report':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/report.json',
      'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
      'X1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
      'N0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
      'xp':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
      'star':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json',
      'L37':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
      'bridge':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/report.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_native_sources_before_full_c1_fold','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};H,HR,X,N,N0,XP,SS,LL,BR=[D[k] for k in paths]
    assert HR['all_passed'] and BR['all_passed']
    origin=[[F(x) for x in r] for r in X['full_X0_56_by768']]
    first=[{row:{k-1:F(c) for k,c in p if k and F(c)} for row,p in col} for col in N['source_columns128_rows48_sparse_dual9']]
    xp=list(map(rec,XP['normalized_Xp27_original768_combinations']))
    star=list(map(rec,SS['actual_unitstar84_original768_combinations']));L=list(map(rec,LL['actual_selected37_original768_constant_combinations']))
    nullchecks=[]
    for kind,cols in [('unitstar84',star),('selectedL37',L)]:
        for j,c in enumerate(cols):
            old={};native={}
            for row in range(56):
                a=sum((v*origin[row][i] for i,v in c.items() if i<640),F(0));b=sum((v*origin[row][i] for i,v in c.items() if i>=640),F(0))
                if a:old[row]=a
                if b:native[row]=b
            assert not old and not native,(kind,j)
            nullchecks.append([kind,j,True])
    mat=[[[F(a) for a in p] for p in row] for row in N0['N']];inv=[[[F(a) for a in p] for p in row] for row in N0['N_inverse']]
    def C4apply(v,M):
        z={}
        for index,a in v.items():
            row,var=divmod(index,8);i,b=divmod(row,4)
            for j in range(12):
                for k,c in enumerate(M[i][j]):
                    if c and not b&k:put(z,8*(4*j+(b|k))+var,a*c)
        return z
    xp0=[]
    for c in xp:
        n={}
        for i,a in c.items():
            if i>=640:
                for row in range(48):
                    if origin[row][i]:put(n,row,a*origin[row][i])
        xp0.append(n)
    mons2=list(itertools.combinations_with_replacement(range(8),2));mi2={m:i for i,m in enumerate(mons2)}
    mons3=list(itertools.combinations_with_replacement(range(8),3));gamma=rec(H['Gamma4_56_times120'])
    rows=[];support=0;nonzero=0;n2nonzero=0;requested=set();nativepairs=set()
    for case in H['all619']:
        c=rec(case['actual_original768_constant_combination']);a1=rec(case['actual_Xp27_linear_a1']);v={}
        for i,a in c.items():
            if i>=640:
                for row,p in first[i-640].items():
                    for var,b in p.items():put(v,8*row+var,a*b)
        for index,a in a1.items():
            p,var=divmod(index,8)
            for row,b in xp0[p].items():put(v,8*row+var,-a*b)
        proper=C4apply(v,inv);assert C4apply(proper,mat)==v
        # General ACTUAL L2 adjoint, including coefficient-basis effects.
        adj={}
        for index,a in gamma.items():
            outrow,m=divmod(index,120)
            if outrow>=48:continue
            mon=mons3[m]
            for var in set(mon):
                pair=list(mon);pair.remove(var);m2=mi2[tuple(pair)]
                for rv,b in proper.items():
                    inrow,w=divmod(rv,8)
                    if w==var:put(adj,(outrow,inrow,m2),-a*b)
        # Sub-adjoint for a free proper N2 direction, with basis multiplication
        # frozen at origin. Other Res4/basis L2 terms are not discarded.
        n2={}
        for (o,i,m),a in adj.items():
            j,cbit=divmod(o,4);h,bbit=divmod(i,4)
            if bbit&cbit==bbit:put(n2,(h,j,cbit^bbit,m),a)
        if adj:nonzero+=1
        if n2:n2nonzero+=1
        requested.update(mons2[k[2]] for k in adj);nativepairs.update(mons2[k[3]] for k in n2)
        support+=len(c)+len(a1)+len(v)+len(proper)+len(adj)+len(n2);assert support<150000,('NATIVE_ADJOINT_SUPPORT',support)
        rows.append({'kernel619_index':case['kernel619_index'],'actual_native_first_residual48_times8':encode(v),
          'proper_S1c0_minus_S0Xp_a1_48_times8':encode(proper),'actual_L2_source_adjoint':enc(adj),
          'free_proper_N2_direction_sub_adjoint':enc(n2)})
    dump(out,'all619-full-c1-native-L2-adjoints.json',{'all619':rows,'monomials2_36':mons2,
      'Gamma4_56_times120':encode(gamma),'all121_separate_origin_checks':nullchecks,
      'full_c1':'-Xp*a1-Xstar*b1-XL*l1','S0_Xstar_and_XL_zero':True,
      'a2_does_not_enter_degree1_source_multiplier':True,
      'actual_L2_remaining':'L2*(S1*c0-S0*Xp*a1)',
      'proper_residual_computed_as':'(N0^t)^-1*(nativeXc1-nativeXp0*a1)',
      'N2_sub_adjoint_scope':'Only free proper N2 directions under frozen-origin C4 multiplication; other basis L2 terms and oldB dependence remain separate',
      'no_actual_L2_N2_values_assigned':True})
    pre={'actual619_source_count':619,'native_L2_nonzero_adjoint_count':nonzero,'free_N2_nonzero_sub_adjoint_count':n2nonzero,
      'actual_L2_pairs':sorted(requested),'free_N2_pairs':sorted(nativepairs),'whole_saved_QQ_support':support,
      'all84_unitstar_and37_selectedL_separate_origins_zero':True,
      'lower_b1_l1_cannot_change_this_L2_sensitivity':True,'no_need_to_expand_b1_l1_a2_for_this_native_sensitivity':True,
      'new_N2_necessity_NOT_proved':'OldB source dependence may cancel coherently; actual L2 coefficients not evaluated',
      'no_new_source_function_jet_coefficient_solve':True,'cap':150000}
    dump(out,'preflight.json',pre);dump(out,'phase.json',{'phase':'FINAL_SOURCE_ONLY','support':support})
    note=ROOT/'research/notes'/f'{STEM}.md'
    r={'id':STEM,'all_passed':True,**pre,'actual619_quartic_values_or_HF4_computed':False,'input_sha256':hs,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump(out,'report.json',r);print(json.dumps(pre))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
