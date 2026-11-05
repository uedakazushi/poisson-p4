#!/usr/bin/env python3
"""NEW restricted own-family order6:114 coefficients, no common/mixed solve."""
import argparse,collections,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63';Z=(F(0),)*7
def mul(a,b):return tuple(sum((a[i]*b[k-i] for i in range(k+1)),F(0)) for k in range(7))
def inv(a):
    assert a[0];v=[1/a[0]]
    for k in range(1,7):v.append(-sum((a[i]*v[k-i] for i in range(1,k+1)),F(0))/a[0])
    assert mul(a,tuple(v))==(F(1),)+(F(0),)*6;return tuple(v)
def bits(a):return max(a.numerator.bit_length(),a.denominator.bit_length())
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    def dump(n,x):
        (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths={
      'input5':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63/exact-two57-univariate5-input.json',
      'input5_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63_REVIEW/report.json',
      'CA5':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/CA-source-and-regular-quotient-second.json',
      'CB5':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/CB-source-and-regular-quotient-second.json',
      'actual5_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
      'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'source3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'own4':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
      'highJac':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/report.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};dump('phase.json',{'phase':'NEW_RESTRICTED_SOURCE6_INPUTS_PERSISTED','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};I,IP,CA,CB,VP,S,C3,O4,HJ=[D[k] for k in paths]
    assert IP['all_passed'] and VP['all_passed'] and HJ['all_passed']
    assert sha(ROOT/paths['actual5_peer'])=='808637882f3d8d28f393a7b0f4abfeec59e076bd1269034bd67d8ac6dd8c706a'
    assert I['physical_common_rline_binding']['actual_common_YB7_8_quadratic_r2_zero']
    assert I['physical_common_rline_binding']['quotient_origin0_from_source3']
    nodes=S['DAG_scalar_nodes'];planned=0;profiles=[]
    for f in I['families']:
        fam=f['family'];old=CA if fam=='CA' else CB
        assert [i for i,p in old['all57_source_jets_through5']]==f['original_unknown57_order']
        assert len(f['fixed_J57_inverse_sparse_rows'])==57
        order=f['source_needed_original_DAG_nodes'];roots=f['original_equations57']+[f['actual_g_original_node']]+f['actual_bbar4_original_nodes']
        uses=collections.Counter()
        for i in order:
            n=nodes[i]
            if n[0] in ['a','m','n','i']:uses.update(n[1:])
        uses.update(roots);live=set();peak=0
        for i in order:
            live.add(i);peak=max(peak,len(live));n=nodes[i]
            if n[0] in ['a','m','n','i']:
                for j in n[1:]:
                    uses[j]-=1
                    if not uses[j]:live.discard(j)
        bound=7*(peak+57+62)+20000;planned=max(planned,bound);assert planned<150000
        profiles.append({'family':fam,'source_nodes':len(order),'peak_live_nodes':peak,'whole_QQ_bound':bound,
          'saved_J57_determinant':f['fixed_J57_determinant'],'new_unknowns':57,'only_degree6_is_new':True,
          'last2_target_slopes':f['last2_target_slopes']})
    dump('preflight.json',{'kind':'NEW HIGHER PHYSICAL RESEARCH CALCULATION','families':profiles,'new_QQ_source_coefficients':114,
      'whole_QQ_bound':planned,'cap':150000,'max_bits_guard':100000,'no_common105_sixth_or_mixedthird_or_gauge3':True,
      'quotient_formula':'h3=(b6-g5*h1)/g3, sinceh0=g4=0','CB_actual_pullback':'add linear_h1 dot savedYB3, sinceYB2=0',
      'h0_global_source_and_J57_ancestry_pinned':True})
    maxbits=0;peakslots=0;outputs={};h3common={}
    def evaluate(order,roots,bind):
        nonlocal maxbits,peakslots
        uses=collections.Counter()
        for i in order:
            n=nodes[i]
            if n[0] in ['a','m','n','i']:uses.update(n[1:])
        uses.update(roots);v={}
        for i in order:
            n=nodes[i];o=n[0]
            if o=='v':p=bind[i]
            elif o=='c':p=(F(n[1],n[2]),)+(F(0),)*6
            elif o=='a':p=tuple(a+b for a,b in zip(v[n[1]],v[n[2]]))
            elif o=='n':p=tuple(-a for a in v[n[1]])
            elif o=='m':p=mul(v[n[1]],v[n[2]])
            elif o=='i':p=inv(v[n[1]])
            else:raise AssertionError((i,n))
            v[i]=p;peakslots=max(peakslots,7*len(v));assert peakslots+20000<150000
            maxbits=max(maxbits,max(bits(a) for a in p));assert maxbits<100000
            if o in ['a','m','n','i']:
                for j in n[1:]:
                    uses[j]-=1
                    if not uses[j]:del v[j]
        return {i:v[i] for i in roots}
    yb=[list(map(F,p)) for i,p in C3['actual_YB_last2_through3']];assert all(not p[2] for p in yb)
    cbfirst=[list(map(F,p)) for p in O4['families']['CB']['quotient_firsts4_by2']]
    for f in I['families']:
        fam=f['family'];old=CA if fam=='CA' else CB
        unknown=f['original_unknown57_order'];eq=f['original_equations57'];g=f['actual_g_original_node'];nums=f['actual_bbar4_original_nodes'];roots=eq+[g]+nums;order=f['source_needed_original_DAG_nodes']
        bind={i:tuple(map(F,p))+(F(0),) for i,p in old['all57_source_jets_through5']};before=evaluate(order,roots,bind)
        slopes=list(map(F,f['last2_target_slopes']));targets=[Z for i in eq];targets[55]=(F(0),slopes[0])+(F(0),)*5;targets[56]=(F(0),slopes[1])+(F(0),)*5
        assert all(before[i][:6]==q[:6] for i,q in zip(eq,targets))
        rhs=[before[i][6] for i in eq];u6=[-sum((F(c)*rhs[j] for j,c in row),F(0)) for row in f['fixed_J57_inverse_sparse_rows']]
        afterbind={i:bind[i][:6]+(c,) for i,c in zip(unknown,u6)}
        dump(f'{fam}-new-sixth-coefficient-checkpoint.json',{'family':fam,'source_unknown_order57':unknown,
          'degree6_residual57_before':list(map(str,rhs)),'new_degree6_corrections57':list(map(str,u6)),
          'saved_J57_determinant':f['fixed_J57_determinant'],'cached0_through5_source_path':paths['CA5' if fam=='CA' else 'CB5']})
        after=evaluate(order,roots,afterbind);assert all(after[i]==q for i,q in zip(eq,targets))
        gp=after[g];assert gp[:3]==(F(0),)*3 and gp[3] and not gp[4]
        own=[];actual=[];quotients=[]
        for j,node in enumerate(nums):
            h1=F(old['regular_quotient_origin0_and_first'][j][1]);h2=F(old['new_four_regular_quotient_second_line'][j]);bp=after[node]
            assert bp[:4]==(F(0),)*4 and bp[4]==gp[3]*h1 and bp[5]==gp[3]*h2
            h3=(bp[6]-gp[5]*h1)/gp[3];h=(F(0),h1,h2,h3,F(0),F(0),F(0));assert mul(gp,h)==bp
            bend=sum((cbfirst[j][k]*yb[k][3] for k in range(2)),F(0)) if fam=='CB' else F(0)
            own.append(h3);actual.append(h3+bend);quotients.append({'bbar_node':node,'bbar_through6':list(map(str,bp)),
              'known_h1':str(h1),'known_h2':str(h2),'new_own_h3':str(h3),'known_common_bending':str(bend),'actual_common_h3':str(h3+bend),
              'all7_coefficients_of_bbar_equals_g_times_h_checked':True})
        ids=[566400,566401,566402,566403] if fam=='CA' else [330164,330289,330524,330594]
        h3common.update(dict(zip(ids,actual)))
        outputs[fam]={'all57_source_equations_through6':[[i,list(map(str,after[i]))] for i in eq],
          'all57_source_jets_through6':[[i,list(map(str,afterbind[i]))] for i in unknown],
          'physical_x9_readbacks_through6':[[i,list(map(str,afterbind[i]))] for i in range(2,11)],
          'g_source_node':g,'g_through6':list(map(str,gp)),'four_quotient_records':quotients,
          'actual_four_common_third_values':list(map(str,actual)),
          'actual_common_pullback_uses_YB2_zero_and_saved_YB3':True}
        dump(f'{fam}-source6-and-regular-quotient-third.json',outputs[fam])
    dump('actual-eight-g-third-values.json',{'actual_eight_common_regular_g_third_values':[[i,str(a)] for i,a in sorted(h3common.items())],
      'actual_four_oldCB_minus_newCA_third_differences':[str(h3common[a]-h3common[b]) for a,b in zip([330164,330289,330524,330594],[566400,566401,566402,566403])],
      'new_restricted_own_source_coefficients':114,'no_common105_sixth_or_general_mixedthird_or_gauge3':True,
      'actual619_quartic_values_or_HF4_computed':False})
    dump('phase.json',{'phase':'FINAL_NEW_RESTRICTED_SOURCE6','new_coefficients':114})
    note=ROOT/'research/notes'/f'{STEM}.md'
    r={'id':STEM,'all_passed':True,'kind':'NEW RESTRICTED OWN-FAMILY HIGHER PHYSICAL RESEARCH','new_own_univariate_order6_coefficients':114,
      'all114_source_equations_checked_through6':True,'all_eight_g_third_values_fixed':True,
      'new_common105_sixth_mixedthird_or_gauge3':False,'actual619_values_HF4_or_length_computed':False,
      'max_live_series_QQ':peakslots,'max_QQ_bitheight':maxbits,'input_sha256':hs,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'independent_source_review':'PENDING'}
    dump('report.json',r);print(json.dumps({k:r[k] for k in ['id','all_passed','new_own_univariate_order6_coefficients','all_eight_g_third_values_fixed','max_live_series_QQ','max_QQ_bitheight','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
