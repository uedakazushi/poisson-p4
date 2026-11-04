#!/usr/bin/env python3
"""SOURCE ONLY: force -t0 first in the ORIGINAL49039-term matrix.

Support/height/serialization forecast only. No coefficient update is built.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_REORDERED_START_PREFLIGHT_20261104_v63'
BASE='S61_SPARSE_GENERIC5_CONSTRAINED83_RAW_BLOCK_20261104_v63'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT_PREFIX_20261104_v63'
LOW='S61_SPARSE_GENERIC5_CONSTRAINED83_LOW_GROWTH_CONTINUATION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def size(x):return len(json.dumps(x,separators=(',',':')))
def support(p):return [] if p is None else [tuple(t[v+1]for v in range(5))for t in p['terms']]
def stat(p):
    if p is None:return 0,Q(0),0
    M=max(abs(int(t[0]))for t in p['terms']);return len(p['terms']),Q(p['content']),(1<<M.bit_length())-1
def lp_bytes(p,content=None,shift=None):
    return size({'content':p['content']if content is None else str(content),'shift':[0]*5 if shift is None else shift,'terms':p['terms']})
def bound_bytes(S,base,M,D):
    if not S or not M:return 0
    return size({'content':'','shift':base,'terms':[]})+len(S)*(len(str(M))+11)+len(str(D))+16+sum(len(str(e[v]-base[v]))for e in S for v in range(5))
def run():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports'/BASE;rp=base/'report.json';r=read(rp);assert r['all_passed'] and sha(rp)=='3b9d32cf070644c7651adf382155f2379702faa967d4e12047a0d730befdc303'
    mp=base/'actual-236x83-matrix.json';assert sha(mp)==r['evidence_sha256'][mp.name]
    matrix=read(mp);A={(i,j):p for i,j,p in matrix['sparse_entries']};assert matrix['total_matrix_terms']==49039
    oldrp=ROOT/'reports'/OLD/'report.json';lowrp=ROOT/'reports'/LOW/'report.json';old=read(oldrp);low=read(lowrp);assert old['all_passed'] and low['all_passed']
    oldfirst=read(oldrp.parent/'pivot-00.json');assert(oldfirst['pivot_original_row'],oldfirst['pivot_original_column'])==(188,37)
    row,col=192,22;pr=A[row,col];assert pr['content']=='-1' and pr['terms']==[['1',0,0,0,0,1]]
    pe=(0,0,0,0,1);unit=Q(-1);live=49039;rr=range(236);cc=range(83)
    pcols=[j for j in cc if j!=col and(row,j)in A];targets=[i for i in rr if i!=row and(i,col)in A]
    R=sum(len(p['terms'])for(i,j),p in A.items()if i==row);C=sum(len(A[i,col]['terms'])for i in targets)
    # Source nonzero coupling proves the first old normalized row changes.
    assert A[188,22]['content']=='-30' and A[188,22]['terms']==[['1',0,1,0,0,0]] and(192,37)not in A
    witness_col=next(j for j in pcols if j!=37)
    distinct={'old_first_pivot':[188,37],'new_first_pivot':[192,22],
              'old_first_row_at_new_pivot_column':A[188,22],
              'new_pivot_row_at_old_pivot_column':'exact_zero',
              'new_pivot_row_other_nonzero_column':witness_col,
              'new_pivot_row_other_source_value':A[192,witness_col],
              'proof':'The new first subtraction changes old row188 by -A188,22*(row192/(-t0)); its witness component is nonzero in the parameter domain. Thus this is not the original literal first instruction. It may represent the same eventual row space, as every correct elimination does; no distinct final algebra is asserted.',
              'existing_24_or29_steps_replayed':False,'future_min_forecast_path_not_computed':True}
    peak=live;logical=live;receipt=4096;residual=8192+32*235*82;updates=[]
    for(i,j),p in A.items():
        if i!=row and j!=col:residual+=lp_bytes(p)
    for j in pcols:
        p=A[row,j];receipt+=lp_bytes(p,Q(p['content'])/unit,[0,0,0,0,-1])+128
    for i in targets:
        factor=A[i,col];nf,cf,mf=stat(factor);fs=support(factor);receipt+=lp_bytes(factor)+512
        for j in pcols:
            oldp=A.get((i,j));a,ca,ma=stat(oldp);np,cp,mpm=stat(A[row,j]);cp/=unit
            ps=[tuple(e[v]-pe[v]for v in range(5))for e in support(A[row,j])]
            prod={tuple(x[v]+y[v]for v in range(5))for x in fs for y in ps};union=set(support(oldp))|prod;k=len(prod);b=len(union)
            support_scratch=len(fs)+len(ps)+a+k+b
            bound=max(logical+2*k+3*b+64,live+support_scratch+64);peak=max(peak,bound)
            cprod=cf*cp;D=ca.denominator*cprod.denominator//gcd(ca.denominator,cprod.denominator)
            M=abs(ca.numerator)*(D//ca.denominator)*ma+abs(cprod.numerator)*(D//cprod.denominator)*min(nf,np)*mf*mpm
            base_shift=[0,0,0,0,-1]
            residual+=bound_bytes(union,base_shift,M,D)-(lp_bytes(oldp)if oldp is not None else 0)
            meta=[i,j,a,k,b,logical,bound,M.bit_length(),D.bit_length()];updates.append(meta);receipt+=size(meta)+1
            logical+=b-a
        logical-=nf
    logical-=R
    future_frontier_reserve=128*235*82
    total_raw=receipt+residual+4096+1024*1024+future_frontier_reserve
    assert peak<=150000 and total_raw<=32*1024*1024
    certificate={'first_original_pivot':[row,col],'pivot':'-t0','parameter_exponent':list(pe),
                 'initial_terms':49039,'pivot_row_terms':R,'factor_column_terms_excluding_pivot':C,
                 'target_rows':len(targets),'other_pivot_columns':len(pcols),'coefficient_updates_forecast_only':len(updates),
                 'whole_pivot_coefficient_and_scratch_bound':peak,'final_support_union_terms_upper':logical,
                 'tight_RAW_bound':{'normalized_row_and_factors_and_update_metadata':receipt,'exact_residual_upper':residual,'selection_metadata':4096,'ordinary_metadata_reserve':1024*1024,'worst_next_selection_frontier_reserve':future_frontier_reserve,'total':total_raw},
                 'all_per_update_source_forecasts':updates,'no_actual_new_coefficients_or_row_operations':True}
    (out/'first-minus-t0-source-forecast.json').write_text(json.dumps(certificate,indent=2)+'\n')
    (out/'different-order-source-witness.json').write_text(json.dumps(distinct,indent=2)+'\n')
    plan={'source':'original236x83 only; no parsed old124840/143215 matrix or source C/N/T/P arrays',
          'policy':'Force original(192,22)=-t0 first. Subsequently minimum completed whole coefficient/scratch forecast among ALL eligible QQ constants or parameter Laurent monomials; source row/column ties; monotone domination allowed. No pivot executes after incomplete selection.',
          'exact_values':'content*x^shift*integerPolynomial with all five parameters; t0 inverse is a parameter-field unit. Signed exponent and representation guards before arithmetic.',
          'scratch':'Whole pivot support preforecast before normalization; L+2k+3b+64 includes original matrix, product, aligned add/sub copies, new/back-add and scalar LCM bounds; support-set scratch separately counted.',
          'RAW':'Proven numeric denominator LCM, primitive-content digit inequality and actual exponent digits. Reserve additionally128*(remainingRows)*(remainingCols) bytes for an uncompleted future selection frontier, beyond1MiB metadata.',
          'time':'47s soft stop before next pivot or unfinished selection;60s hard including final exact residual output,2GB. No guaranteed number of future pivots.',
          'inductive_whole_worker_bound':'Start49039; before EACH proposed whole pivot require all current+scratch and complete RAW bounds. Release retired row and scratch, then apply the same guard to the new actual state. Stop with exact residual, without increasing limits or changing to nonmonomial pivots.',
          'not_authorized_or_executed_here':['arithmetic branch','gcd','series','nonmonomial pivot','T/P expansion','Gram','dense inverse','new parameter point or prime'],
          'first_step_safe_only_not_a_promise_for_all83':True}
    (out/'whole-worker-branch-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    paths=[rp,mp,oldrp,lowrp,oldrp.parent/'pivot-00.json',ROOT/f'code/{Path(__file__).name}',ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_constrained83_low_growth_continuation_20261104_v63.cpp']
    result={'id':ID,'all_passed':True,'source_only_preflight':True,'branch_is_not_a_literal_existing_path_replay':True,
            'first_pivot':'original(192,22)=-t0','first_step_live_and_scratch_bound':peak,'first_step_final_terms_upper':logical,
            'first_step_complete_RAW_bound':total_raw,'first_step_safe':True,
            'all_five_parameters_constraints_and_old_results_preserved':True,'arithmetic_worker_or_series_not_run':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':run()
