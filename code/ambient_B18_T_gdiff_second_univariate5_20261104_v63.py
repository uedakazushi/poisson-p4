#!/usr/bin/env python3
"""ONE fifth-order univariate correction per own-family57, after root GO.

No full105 fifth solve, mixed jet, inverse search, or quotient zero-unit.
"""
import argparse,collections,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63'
Z=(F(0),)*6
def dump(out,name,x):
    p=out/name;p.write_text(json.dumps(x,separators=(',',':'))+'\n')
    assert sum(q.stat().st_size for q in out.iterdir())<32*1024*1024
def mul(a,b):return tuple(sum((a[i]*b[k-i] for i in range(k+1)),F(0)) for k in range(6))
def inv(a):
    assert a[0];q=[1/a[0]]
    for k in range(1,6):q.append(-sum((a[i]*q[k-i] for i in range(1,k+1)),F(0))/a[0])
    return tuple(q)
def bits(q):return max(q.numerator.bit_length(),q.denominator.bit_length())
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'input':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63/exact-two57-univariate5-input.json',
      'input_report':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63/report.json',
      'input_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63_REVIEW/report.json',
      'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'known':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase2-known-parts-and-missing-independent-projection.json',
      'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
      'completeness':'reports/AMBIENT_B18_T_PURE_CUBIC619_COMPLETENESS_PREFLIGHT_20261104_v63/report.json',
      'lower':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json',
      'upper':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'all_sources_bound_before_new114_coefficients','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};I,IR,IP,S,K,M,C,L,U=[D[k] for k in paths]
    assert all(d['all_passed'] for d in [IR,IP,C,L,U]);assert IR['new_source5_unknown_coefficients_needed']==114
    assert sha(ROOT/paths['input_peer'])=='917f2492740394b7ebbdeca99886cfa4cbbe18018c2171eaa5788af0c463fdd5'
    assert I['physical_common_rline_binding']['actual_common_YB7_8_quadratic_r2_zero']
    assert I['physical_common_rline_binding']['quotient_origin0_from_source3']
    nodes=S['DAG_scalar_nodes'];families={};g2={};maxbits=0;peakslots=0
    def evaluate(order,roots,bind):
        nonlocal maxbits,peakslots
        uses=collections.Counter()
        for i in order:
            n=nodes[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
        uses.update(roots);v={}
        for i in order:
            n=nodes[i];op=n[0]
            if op=='v':q=bind[i]
            elif op=='c':q=(F(n[1],n[2]),F(0),F(0),F(0),F(0),F(0))
            elif op=='a':q=tuple(a+b for a,b in zip(v[n[1]],v[n[2]]))
            elif op=='n':q=tuple(-a for a in v[n[1]])
            elif op=='m':q=mul(v[n[1]],v[n[2]])
            elif op=='i':q=inv(v[n[1]])
            else:raise AssertionError((i,n))
            v[i]=q;peakslots=max(peakslots,6*len(v));assert peakslots+2000<150000
            maxbits=max(maxbits,max(bits(a) for a in q));assert maxbits<100000
            for j in n[1:3] if op in ['a','m'] else n[1:] if op in ['n','i'] else []:
                uses[j]-=1
                if not uses[j]:v.pop(j,None)
        return {i:v[i] for i in roots}
    for f in I['families']:
        family=f['family'];unknown=f['original_unknown57_order'];eq=f['original_equations57'];g=f['actual_g_original_node'];nums=f['actual_bbar4_original_nodes']
        roots=eq+[g]+nums;order=f['source_needed_original_DAG_nodes'];slopes=list(map(F,f['last2_target_slopes']))
        bind={i:tuple(map(F,p))+(F(0),) for i,p in f['cached_univariate_source_through4']};assert list(bind)==unknown
        before=evaluate(order,roots,bind)
        targets=[Z for i in eq];targets[55]=(F(0),slopes[0],F(0),F(0),F(0),F(0));targets[56]=(F(0),slopes[1],F(0),F(0),F(0),F(0))
        assert all(before[i][:5]==t[:5] for i,t in zip(eq,targets))
        residual=[before[i][5] for i in eq];delta=[]
        for row in f['fixed_J57_inverse_sparse_rows']:delta.append(-sum((F(c)*residual[j] for j,c in row),F(0)))
        afterbind={i:bind[i][:5]+(c,) for i,c in zip(unknown,delta)}
        # Save the new114coefficients and residuals BEFORE any assertion or
        # quotient extraction; this is the true new math checkpoint.
        dump(out,f'{family}-fifth-coefficient-checkpoint.json',{'family':family,'actual_unknown_order57':unknown,
          'cached0_through4':f['cached_univariate_source_through4'],'degree5_residual57_before':list(map(str,residual)),
          'new_degree5_correction57':list(map(str,delta)),'unit_Jacobian_determinant':f['fixed_J57_determinant']})
        after=evaluate(order,roots,afterbind);assert all(after[i]==t for i,t in zip(eq,targets))
        gp=after[g];assert gp[:3]==(F(0),)*3 and gp[3] and gp[4]==0
        h1=list(map(F,f['cached_four_quotient_first_line']));values=[]
        for node,c1 in zip(nums,h1):
            p=after[node];assert p[:4]==(F(0),)*4 and p[4]==gp[3]*c1
            c2=(p[5]-gp[4]*c1)/gp[3]
            h=(F(0),c1,c2,F(0),F(0),F(0));assert mul(gp,h)==p
            values.append(c2)
        families[family]={'all57_source_equations_through5':[[i,list(map(str,after[i]))] for i in eq],
          'all57_source_jets_through5':[[i,list(map(str,afterbind[i]))] for i in unknown],
          'actual_physical_x9_readbacks_through5':[[i,list(map(str,afterbind[i]))] for i in range(2,11)],
          'g_source_node':g,'g_through5':list(map(str,gp)),
          'four_bbar_source_nodes':nums,'four_bbar_through5':[list(map(str,after[i])) for i in nums],
          'regular_quotient_origin0_and_first':[[str(F(0)),str(c)] for c in h1],
          'new_four_regular_quotient_second_line':list(map(str,values)),
          'division_is_r3_then_nonzero_g3_unit':True,'all6_coefficients_of_bbar_equal_g_times_h_verified':True}
        for i,c in zip(([566400,566401,566402,566403] if family=='CA' else [330164,330289,330524,330594]),values):g2[i]=c
        dump(out,f'{family}-source-and-regular-quotient-second.json',families[family])
        dump(out,'phase.json',{'phase':f'{family}_source5_and_g2_complete','new_coefficients_so_far':57*len(families)})
    assert len(g2)==8
    actual=[]
    for j,(known,row) in enumerate(zip(K['all619_known_part_missing_coefficients_split'],M['all619'])):
        assert row['kernel619_index']==j and not row['opaque_third_r3_sensitivity']
        value=F(known)+sum((F(c)*g2[i] for i,c in row['opaque_second_r2_sensitivity']),F(0));actual.append(value)
    nonzero=[j for j,v in enumerate(actual) if v];zero=not nonzero
    final={'actual_eight_regular_g_second_values':[[i,str(c)] for i,c in sorted(g2.items())],
      'actual_four_oldCB_minus_newCA_second_differences':[str(g2[i]-g2[j]) for i,j in zip([330164,330289,330524,330594],[566400,566401,566402,566403])],
      'all619_actual_pure_cubic_relation_values':list(map(str,actual)),
      'nonzero_pure_cubic_witness_indices':nonzero,'all_source_and_coordinate_bending_conditions_verified':True,
      'actual_m3_T_zero':bool(nonzero),'actual_T_length':10 if nonzero else None,
      'actual_HF_degree0_through3':[5,3,2,0 if nonzero else 1],
      'actual_T_length_interval':[10,10] if nonzero else [11,17],
      'full_nil_cube_claimed':False,'source_review_status':'PENDING'}
    dump(out,'actual-gdiff-second-and-complete619-pure-cubic-values.json',final)
    dump(out,'phase.json',{'phase':'FINAL','all619_pure_values_zero':zero,'T_length':10 if nonzero else None,'HF':final['actual_HF_degree0_through3']})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'new_source_coefficients_computed':114,'source_parameter_directions':1,
      'all114_fifth_coefficients_checked_by114_source_equations':True,'actual_nonzero_pure619_cubic_values':len(nonzero),
      'actual_HF_through_degree3':final['actual_HF_degree0_through3'],'actual_T_length':final['actual_T_length'],
      'actual_T_length_interval':final['actual_T_length_interval'],'full_nil_cube_not_inferred':True,
      'new_common105_order5_or_fullmixedjet_run':False,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'max_live_series_QQ_slots':peakslots,
      'max_QQ_bitheight':maxbits,'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      'independent_source_and_completeness_review_status':'PENDING'}
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','new_source_coefficients_computed','actual_nonzero_pure619_cubic_values','actual_HF_through_degree3','actual_T_length','actual_T_length_interval','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
