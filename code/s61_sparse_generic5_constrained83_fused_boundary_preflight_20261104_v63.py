#!/usr/bin/env python3
"""Read-only logical-final vs scratch diagnosis; structural content screen.

No polynomial coefficient multiplication/update, pivot, gcd, or series.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FUSED_BOUNDARY_PREFLIGHT_20261104_v63'
PARENT='S61_SPARSE_GENERIC5_CONSTRAINED83_LOW_GROWTH_CONTINUATION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def support(p):
    return [] if p is None else [tuple(p['shift'][v]+t[v+1]for v in range(5))for t in p['terms']]
def run():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PARENT;rp=folder/'report.json';r=read(rp)
    assert r['all_passed'] and sha(rp)=='0b4bd1b5530d13e58ea0674753c7c60777ce3278c5022c269a5af0038dcd86c7'
    mp=folder/'residual-matrix.json';assert sha(mp)==r['evidence_sha256'][mp.name]
    m=read(mp);A={(i,j):p for i,j,p in m['sparse_entries']};rr=m['row_original_ids'];cc=m['column_original_ids'];live=m['total_terms']
    assert m['shape']==[207,54] and live==143215
    front=read(folder/'selection-frontier.json');assert len(front['records'])==front['eligible_candidates']==467 and all(x[2]=='terms'for x in front['records'])
    rowterms={i:sum(len(A[i,j]['terms'])for j in cc if(i,j)in A)for i in rr}
    colterms={j:sum(len(A[i,j]['terms'])for i in rr if(i,j)in A)for j in cc}
    monomials={(i,j)for(i,j),p in A.items()if len(p['terms'])==1}
    assert len(monomials)==467
    mono_rows={i for i,j in monomials};mono_cols={j for i,j in monomials}
    screen={
        'UFD_reason':'After clearing any Laurent coordinate powers, a slice containing a nonzero parameter monomial has no common nonmonomial factor. Monomial rescaling changes no coefficient term count.',
        'rows_with_monomial':sorted(mono_rows),'columns_with_monomial':sorted(mono_cols),
        'nonzero_rows_without_monomial':[{'row':i,'terms':rowterms[i],'entries':sum((i,j)in A for j in cc)}for i in rr if rowterms[i]and i not in mono_rows],
        'nonzero_columns_without_monomial':[{'column':j,'terms':colterms[j],'entries':sum((i,j)in A for i in rr)}for j in cc if colterms[j]and j not in mono_cols],
        'zero_rows':[i for i in rr if not rowterms[i]],
        'all467_eligible_pivot_rows_and_columns_exclude_nonmonomial_content':all(i in mono_rows and j in mono_cols for i,j in monomials),
        'no_gcd_or_content_division_computed':True,
        'unscreened_slices_are_only_opportunities_not_nontrivial_factor_claims':True}
    (out/'structural-content-screen.json').write_text(json.dumps(screen,indent=2)+'\n')
    def forecast(row,col):
        pivot=A[row,col];assert len(pivot['terms'])==1;pe=support(pivot)[0]
        pcols=[j for j in cc if j!=col and(row,j)in A]
        targets=[i for i in rr if i!=row and(i,col)in A]
        R=rowterms[row];C=colterms[col]-1
        changes=[];support_peak=live;old_peak=live;fused_peak=live
        logical=live;max_logical=live;perrow=[]
        for i in targets:
            fs=support(A[i,col]);updates=[]
            for j in pcols:
                if time.monotonic()-started>45:return {'row':row,'column':col,'complete':False,'reason':'read_only_time_guard'}
                ps=[tuple(e[v]-pe[v]for v in range(5))for e in support(A[row,j])]
                old=support(A.get((i,j)));a=len(old)
                prod=set();united=set(old)
                for f in fs:
                    for p in ps:
                        e=tuple(f[v]+p[v]for v in range(5));prod.add(e);united.add(e)
                        scratch=len(fs)+len(ps)+len(old)+len(prod)+len(united)
                        support_peak=max(support_peak,live+scratch+96)
                        if support_peak>150000:return {'row':row,'column':col,'complete':False,'reason':'read_only_support_scratch_guard','support_forecast_peak':support_peak}
                k=len(prod);b=len(united)
                old_peak=max(old_peak,logical+2*k+3*b+64)
                fused_peak=max(fused_peak,logical+b+96)
                updates.append((i,j,a,k,b,len(ps)));changes.append((i,j,a,k,b,len(ps)))
                logical+=b-a;max_logical=max(max_logical,logical)
            perrow.append((i,len(fs),updates));logical-=len(fs)
        final=logical-R
        # External normalized pivot row, read one entry at a time. This is a
        # different representation, assessed only; no value is serialized here.
        external_rowmajor=live-R;ext_r_peak=external_rowmajor
        for i,nf,updates in perrow:
            for ii,j,a,k,b,np in updates:
                ext_r_peak=max(ext_r_peak,external_rowmajor+np+b+96)
                external_rowmajor+=b-a
            external_rowmajor-=nf
        assert external_rowmajor==final
        bycolumn={j:[]for j in pcols}
        for change in changes:bycolumn[change[1]].append(change)
        external_colmajor=live-R;ext_c_peak=external_colmajor
        for j in pcols:
            np=len(A[row,j]['terms'])
            for i,jj,a,k,b,_ in bycolumn[j]:
                ext_c_peak=max(ext_c_peak,external_colmajor+np+b+96)
                external_colmajor+=b-a
        external_colmajor-=C
        assert external_colmajor==final
        value=Q(pivot['content'])*int(pivot['terms'][0][0])
        return {'row':row,'column':col,'complete':True,'pivot_QQ':str(value),'pivot_exponent':list(pe),
                'pivot_row_terms_R':R,'target_factor_column_terms_C':C,'updates':len(changes),
                'initial_logical_terms':live,'final_logical_support_union_upper':final,
                'max_processing_logical_terms_with_pivot_row_resident':max_logical,
                'old_full_polynomial_copy_scratch_peak':old_peak,
                'fused_rowmajor_pivot_resident_peak':fused_peak,
                'fused_external_pivot_row_rowmajor_peak':ext_r_peak,
                'fused_external_pivot_row_columnmajor_peak':ext_c_peak,
                'read_only_support_scratch_peak':support_peak,
                'minimum_proposed_whole_forecast':max(support_peak,min(fused_peak,ext_r_peak,ext_c_peak)),
                'final_union_bound_fits150k':final<=150000,
                'new_actual_coefficients_computed':False}
    # All old blocked candidates are considered in their deterministic cheap
    # growth order until a read-only resource stop. No pivot is selected/applied.
    results=[]
    for row,col,status,peak,b in front['records']:
        if time.monotonic()-started>45:break
        z=forecast(row,col);results.append(z)
        if not z['complete'] and z['reason']=='read_only_time_guard':break
    complete=[z for z in results if z['complete']]
    completed_all=len(results)==467 and all(z['complete']for z in results)
    best=sorted(complete,key=lambda z:(z['minimum_proposed_whole_forecast'],z['row'],z['column']))
    (out/'logical-final-and-fused-forecasts.json').write_text(json.dumps({'all467_completed':completed_all,'completed_count':len(complete),'examined_count':len(results),'best_forecast_among_completed':best[:12],'all_results':results,'no_actual_row_update_or_pivot':True},indent=2)+'\n')
    report={'id':ID,'all_passed':True,'read_only_forecast_and_structural_screen_only':True,
            'source143215_term_residual_unchanged':True,'all467_full_support_forecasts_completed':completed_all,
            'completed_candidate_forecasts':len(complete),
            'best_completed_candidate':None if not best else {k:v for k,v in best[0].items()if k not in('new_actual_coefficients_computed',)},
            'all_candidate_pivot_rows_columns_have_no_nonmonomial_content':screen['all467_eligible_pivot_rows_and_columns_exclude_nonmonomial_content'],
            'no_new_pivot_gcd_or_series_worker':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [rp,mp,folder/'selection-frontier.json',ROOT/f'code/{Path(__file__).name}',ROOT/f'research/notes/{ID}.md']},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':run()
