#!/usr/bin/env python3
"""NEW nine source-fourth functionals from cached≤3 and fixed105 adjoints."""
import argparse,collections,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63'
Z=(F(0),)*5
def mul(a,b):return tuple(sum((a[i]*b[k-i] for i in range(k+1)),F(0)) for k in range(5))
def inv(a):
    assert a[0];v=[1/a[0]]
    for k in range(1,5):v.append(-sum((a[i]*v[k-i] for i in range(1,k+1)),F(0))/a[0])
    assert mul(a,tuple(v))==(F(1),)+(F(0),)*4;return tuple(v)
def main(out):
    st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    def dump(n,x):
        (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths={
      'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
      'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
      'cache3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'cache2':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
      'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'higher_targets':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/phase4-all619-high-old-source-sensitivities.json',
      'copy_chain':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/phase1-source-closure-and-bindings.json',
      'ownCA6':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63/CA-source6-and-regular-quotient-third.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};dump('phase.json',{'phase':'NEW_NINE_FOURTH_FUNCTIONAL_INPUTS_BOUND','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};J,JP,C3,C2,S,H,COPY,CA=[D[k] for k in paths]
    assert JP['all_passed'] and not COPY['missing_source_leaves']
    unknown=J['all105_unknown_original_node_order'];eq=J['all105_equation_original_node_order'];index={i:k for k,i in enumerate(unknown)}
    targets=sorted({leaf for row in H['all619'] for (kind,leaf),a in row['old_high_source4_sensitivity'] if F(a) and kind=='physical'})
    assert len(targets)==9 and set(targets)<=index.keys()
    assert all(kind=='physical' for row in H['all619'] for (kind,leaf),a in row['old_high_source4_sensitivity'])
    bind={i:tuple(map(F,p))+(F(0),) for i,p in C3['actual_common105_rline_through3']};assert set(bind)==set(unknown)
    mon20=[tuple(m) for m in C2['base20_monomials_original8']];im={m:i for i,m in enumerate(mon20)}
    for i,p in C2['all105_actual_source2_values']:
        q={j:F(a) for j,a in p};assert bind[i][:3]==(q.get(0,F(0)),q.get(im[(7,)],F(0)),q.get(im[(7,7)],F(0)))
    nodes=S['DAG_scalar_nodes'];needed=set();todo=eq[:]
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i)
        if i in bind:continue
        n=nodes[i]
        if n[0] in ['a','m','n','i']:todo.extend(n[1:])
        else:assert n[0]=='c',('UNBOUND_SOURCE_VARIABLE',i,n)
    uses=collections.Counter()
    for i in needed:
        if i not in bind and nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
    uses.update(eq);live=set();peak=0
    for i in sorted(needed):
        live.add(i);peak=max(peak,len(live))
        if i not in bind and nodes[i][0] in ['a','m','n','i']:
            for j in nodes[i][1:]:
                uses[j]-=1
                if not uses[j]:live.discard(j)
    forecast=5*(peak+105)+2*105*105+9*105+10000;assert forecast<150000
    dump('preflight.json',{'kind':'NEW HIGHER SOURCE FUNCTIONAL RESEARCH','targets9':targets,'source_DAG_closure':len(needed),'peak_live_series':peak,
      'whole_QQ_bound':forecast,'cap':150000,'actual_all105_fourth_coefficients_solved':False,
      'only_cached_source0_through3_and_fixed105_inverse':True,'no_mixedthird_or_gauge3':True})
    jac=[list(map(F,p)) for p in J['fixed105_Jacobian']];ji=[list(map(F,p)) for p in J['fixed105_Jacobian_inverse']]
    rho=[ji[index[t]] for t in targets]
    for row,t in zip(rho,targets):
        assert all(sum((a*jac[i][j] for i,a in enumerate(row)),F(0))==F(j==index[t]) for j in range(105))
    uses=collections.Counter()
    for i in needed:
        if i not in bind and nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
    uses.update(eq);v={};maxbits=0
    for i in sorted(needed):
        n=nodes[i];o=n[0]
        if i in bind:p=bind[i]
        elif o=='c':p=(F(n[1],n[2]),)+(F(0),)*4
        elif o=='a':p=tuple(a+b for a,b in zip(v[n[1]],v[n[2]]))
        elif o=='n':p=tuple(-a for a in v[n[1]])
        elif o=='m':p=mul(v[n[1]],v[n[2]])
        elif o=='i':p=inv(v[n[1]])
        else:raise AssertionError((i,n))
        v[i]=p;maxbits=max(maxbits,max(max(a.numerator.bit_length(),a.denominator.bit_length()) for a in p));assert maxbits<100000
        if i not in bind:
            for j in n[1:] if o in ['a','m','n','i'] else []:
                uses[j]-=1
                if not uses[j]:del v[j]
    targets_eq=[Z for i in eq];targets_eq[55]=targets_eq[56]=(F(0),F(1),F(0),F(0),F(0))
    assert all(v[i][:4]==q[:4] for i,q in zip(eq,targets_eq))
    rhs=[v[i][4] for i in eq];values=[-sum((a*b for a,b in zip(row,rhs)),F(0)) for row in rho]
    own={i:list(map(F,p)) for i,p in CA['all57_source_jets_through6']};cross=[]
    for t,a in zip(targets,values):
        if t in own:assert own[t][4]==a;cross.append(t)
    assert cross,'at least the CA target has an independent saved order6 readback'
    dump('actual-nine-fourth-source-functionals.json',{'source_unknown_order105':unknown,'source_equation_order105':eq,
      'target_nodes9':targets,'rho9_by105':list(map(lambda p:list(map(str,p)),rho)),'all945_rhoJ_equal_target_entries_checked':True,
      'all105_source_equation_coefficients0_through3_checked':True,'known_F4_residual_with_unknown_fourths_zero':list(map(str,rhs)),
      'actual_nine_fourth_values':[[t,str(a)] for t,a in zip(targets,values)],
      'saved_ownCA6_readback_crosschecks':cross,'full105_fourth_vector_not_constructed':True,
      'proof':'Coefficient4 equation J0*u4+F4(cached_u<=3)=0; rhoJ0=e_target, so requestedu4=-rhoF4',
      'coherent_copied57_targets':'allcopy4 dependence already folded to these physical targets in pinnedhighJac; no independentcopycoeff is requested'})
    dump('phase.json',{'phase':'FINAL_NEW_NINE_ADJOINT_FUNCTIONALS','targets':targets})
    note=ROOT/'research/notes'/f'{STEM}.md'
    r={'id':STEM,'all_passed':True,'kind':'NEW RESTRICTED HIGHER SOURCE FUNCTIONAL RESEARCH','actual_fourth_functionals':9,
      'full105_fourth_solve_or_vector_constructed':False,'all_rhoJ_and_lower_source_equations_verified':True,
      'source_mixedthird_gauge3_or_source6_computed':False,'actual619_values_or_HF4_computed':False,
      'whole_QQ_bound':forecast,'max_QQ_bits':maxbits,'input_sha256':hs,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-st,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'independent_review':'PENDING'}
    dump('report.json',r);print(json.dumps({k:r[k] for k in ['id','all_passed','actual_fourth_functionals','whole_QQ_bound','max_QQ_bits','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
