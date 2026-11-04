#!/usr/bin/env python3
"""READ ONLY: inspect alternative proven forecasts, never update a matrix cell.

Only support Minkowski sums, existing scalar contents and integer height
bounds are computed. No polynomial coefficient product/subtraction, pivot,
series evaluation or source recomputation occurs.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_TIGHT_FORECAST_PREFLIGHT_20261104_v63'
PARENT='S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT_PREFIX_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def encoded_bytes(x):return len(json.dumps(x,separators=(',',':')))
def digits(x):return len(str(abs(x)))
def stat(p):
    if p is None:return 0,Q(0),0,0,0
    c=Q(p['content']);m=max(abs(int(t[0]))for t in p['terms'])
    return len(p['terms']),c,m,m.bit_length()+c.numerator.bit_length(),c.denominator.bit_length()
def supp(p):
    if p is None:return []
    return [tuple(p['shift'][v]+t[v+1]for v in range(5))for t in p['terms']]
def old_cell_bytes(n,N,D):
    def d(bits):return (bits*30103+99999)//100000+2
    return 160+n*(d(N)+d(D)+5*22+48) if n else 0
def common_denominator_primitive_bytes(support,base,M,D):
    """A rigorous byte bound for content*gcd-normalized integer serialization.

    For unreduced integer numerators |H_i|<=M and common D, write H_i=g Z_i.
    digits(g)+sum digits(Z_i)<=n*digits(M)+1; a sign per term and content,
    denominator once, exact support exponents, and punctuation are included.
    """
    if not support or not M:return 0
    n=len(support)
    header=encoded_bytes({'content':'','shift':base,'terms':[]})
    expdigits=sum(digits(e-base[v])for e0 in support for v,e in enumerate(e0))
    assert all(e-base[v]>=0 for e0 in support for v,e in enumerate(e0))
    return header+n*(digits(M)+11)+digits(D)+16+expdigits
def run():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PARENT;reportp=folder/'report.json';r=read(reportp)
    assert sha(reportp)=='31ba1bf17d5f78015125e070147b2250353c3054504ef16f1e34e12f12c14e92' and r['all_passed']
    matrixp=folder/'residual-matrix.json';assert sha(matrixp)==r['evidence_sha256'][matrixp.name]
    matrix=read(matrixp);A={(i,j):p for i,j,p in matrix['sparse_entries']};rr=matrix['row_original_ids'];cc=matrix['column_original_ids']
    live=matrix['total_terms'];assert live==124840
    actual_size={key:encoded_bytes(p) for key,p in A.items()}
    heights={key:(s[0],s[3],s[4]) for key,p in A.items() for s in [stat(p)]}
    rowterms={i:sum(len(A[i,j]['terms']) for j in cc if(i,j)in A)for i in rr}
    colterms={j:sum(len(A[i,j]['terms']) for i in rr if(i,j)in A)for j in cc}
    rowcols={i:{j for j in cc if(i,j)in A}for i in rr}
    base_files=[folder/'preflight.json',folder/'native.log',*[folder/f'pivot-{i:02d}.json'for i in range(24)]]
    old_base_bytes=sum(p.stat().st_size for p in base_files)
    coefficients={
        'max_content_numerator_bits':max(Q(p['content']).numerator.bit_length()for p in A.values()),
        'max_content_denominator_bits':max(Q(p['content']).denominator.bit_length()for p in A.values()),
        'max_primitive_integer_bits':max(abs(int(t[0])).bit_length()for p in A.values()for t in p['terms']),
        'maximum_actual_exponents':[max(e[v] for p in A.values()for e in supp(p))for v in range(5)]}
    def forecast(row,col):
        pivot=A[row,col];assert len(pivot['terms'])==1
        pc=Q(pivot['content'])*int(pivot['terms'][0][0]);pe=supp(pivot)[0]
        uN=pc.numerator.bit_length();uD=pc.denominator.bit_length()
        pcols=[j for j in cc if j!=col and (row,j)in A]
        targets=[i for i in rr if i!=row and (i,col)in A]
        old_res=8192+32*(len(rr)-1)*(len(cc)-1)
        new_res=old_res
        for (i,j),p in A.items():
            if i!=row and j!=col:
                n,N,D=heights[i,j];old_res+=old_cell_bytes(n,N,D);new_res+=actual_size[i,j]
        old_receipt=4096;new_receipt=4096;reference_receipt=4096+64*(len(pcols)+len(targets))
        norm_stats={}
        for j in pcols:
            p=A[row,j];n,c,m,N,D=stat(p);cs=c/pc
            old_receipt+=old_cell_bytes(n,N+uD,D+uN)+128
            sh=[p['shift'][v]-pe[v]for v in range(5)]
            new_receipt+=actual_size[row,j]-len(p['content'])+len(str(cs))+encoded_bytes(sh)-encoded_bytes(p['shift'])+128
            norm_stats[j]=(n,cs,m,N+uD,D+uN,[tuple(e[v]-pe[v]for v in range(5))for e in supp(p)],sh)
        upper=live;peak=live;max_scratch=0;updates=0;detailed=[]
        maxima={'old_denominator_bits':0,'new_common_denominator_bits':0,'new_integer_numerator_bound_bits':0}
        for i in targets:
            factor=A[i,col];nf,cf,mf,Nf,Df=stat(factor);sf=supp(factor)
            old_receipt+=old_cell_bytes(nf,Nf,Df)+512
            new_receipt+=actual_size[i,col]+512
            perrow=[]
            for j in pcols:
                old=A.get((i,j));a,ca,ma,Na,Da=stat(old)
                np,cp,mp,Np,Dp,sp,pshift=norm_stats[j]
                product_support={tuple(x[v]+y[v]for v in range(5))for x in sf for y in sp}
                union=set(supp(old))|product_support;k=len(product_support);b=len(union)
                scratch=len(sf)+len(sp)+a+k+b
                max_scratch=max(max_scratch,scratch)
                bound=max(upper+2*k+3*b+16,live+scratch+16);peak=max(peak,bound)
                lg=(min(nf,np)-1).bit_length();Nprod=Nf+Np+lg+1;Dprod=Df+Dp
                Nold=max(Na+Dprod,Nprod+Da)+1 if a else Nprod
                Dold=Da+Dprod if a else Dprod
                old_res+=old_cell_bytes(b,Nold,Dold)-old_cell_bytes(a,Na,Da)
                # This is only a scalar common-content/height calculation.
                # No coefficient of the polynomial product is computed.
                cprod=cf*cp
                D=ca.denominator*cprod.denominator//gcd(ca.denominator,cprod.denominator)
                M=abs(ca.numerator)*(D//ca.denominator)*ma+abs(cprod.numerator)*(D//cprod.denominator)*min(nf,np)*mf*mp
                old_shift=old['shift'] if old is not None else [0]*5
                prodshift=[factor['shift'][v]+pshift[v]for v in range(5)]
                base=[min(old_shift[v],prodshift[v])for v in range(5)]
                size=common_denominator_primitive_bytes(union,base,M,D)
                new_res+=size-(actual_size.get((i,j),0))
                meta=[j,a,k,b,upper,bound,M.bit_length(),D.bit_length()]
                new_receipt+=encoded_bytes(meta)+1
                old_receipt+=320
                reference_receipt+=encoded_bytes([j,k,b,bound])+1
                maxima['old_denominator_bits']=max(maxima['old_denominator_bits'],Dold)
                maxima['new_common_denominator_bits']=max(maxima['new_common_denominator_bits'],D.bit_length())
                maxima['new_integer_numerator_bound_bits']=max(maxima['new_integer_numerator_bound_bits'],M.bit_length())
                perrow.append(meta);upper+=b-a;updates+=1
            detailed.append({'row':i,'entry_forecast_metadata':perrow})
            upper-=nf
        upper-=rowterms[row]
        reserve=1024*1024
        return {'row':row,'column':col,'pivot_QQ':str(pc),'pivot_exponent':list(pe),
                'target_rows':len(targets),'pivot_other_columns':len(pcols),'updates':updates,
                'whole_pivot_live_bound':peak,'max_support_scratch_slots':max_scratch,
                'final_residual_terms_upper':upper,
                'old_bound_breakdown':{'existing_pre_stop_files':old_base_bytes,'next_receipt':old_receipt,'residual':old_res,'reserve':reserve},
                'old_complete_output_bytes_upper':old_base_bytes+old_receipt+old_res+reserve,
                'tight_full_coefficient_receipt_breakdown':{'existing_pre_stop_files':old_base_bytes,'next_receipt':new_receipt,'residual':new_res,'reserve':reserve},
                'tight_complete_output_bytes_upper':old_base_bytes+new_receipt+new_res+reserve,
                'tight_state_reference_receipt_bytes_upper':reference_receipt,
                'tight_complete_with_reference_bytes_upper':old_base_bytes+reference_receipt+new_res+reserve,
                'new_continuation_folder_output_upper':new_receipt+new_res+reserve,
                'height_bounds':maxima,'new_candidate_not_executed':True,
                'all_per_update_support_only_forecasts':detailed}
    stopped=forecast(23,65)
    assert stopped['old_complete_output_bytes_upper']==37830400,stopped['old_complete_output_bytes_upper']
    assert stopped['whole_pivot_live_bound']==146462
    (out/'stopped-candidate-tight-forecast.json').write_text(json.dumps(stopped,indent=2)+'\n')
    candidates=[]
    for (i,j),p in A.items():
        if len(p['terms'])!=1:continue
        e=supp(p)[0];phase=int(any(e))
        fill=sum(len(rowcols[i]-rowcols[k])for k in rr if k!=i and(k,j)in A)
        pair=(rowterms[i]-1)*(colterms[j]-1)
        candidates.append({'row':i,'column':j,'phase':'constant'if not phase else'monomial','term_pair_majorant':pair,'structural_fill':fill})
    candidates.sort(key=lambda c:(c['term_pair_majorant'],c['structural_fill'],c['row'],c['column']))
    selected=[]
    for phase in ('constant','monomial'):
        selected.extend([c for c in candidates if c['phase']==phase][:3])
    comparisons=[]
    for c in selected:
        f=stopped if(c['row'],c['column'])==(23,65)else forecast(c['row'],c['column'])
        comparisons.append({**c,**{k:f[k]for k in ['whole_pivot_live_bound','final_residual_terms_upper','tight_complete_output_bytes_upper','tight_complete_with_reference_bytes_upper','updates']},'no_pivot_executed':True})
    (out/'coefficient-growth-candidate-comparison.json').write_text(json.dumps({'candidate_scores_no_elimination':candidates,'selected_readonly_forecasts':comparisons,'selection_changed_or_executed':False},indent=2)+'\n')
    report={'id':ID,'all_passed':True,'read_only_source_forecast_only':True,
            'old_stopped_32MiB_forecast_reproduced_exactly':True,
            'current_residual_coefficient_terms_unchanged':124840,
            'same_candidate_tight_live_bound':stopped['whole_pivot_live_bound'],
            'same_candidate_tight_output_bytes_upper':stopped['tight_complete_output_bytes_upper'],
            'same_candidate_state_reference_output_bytes_upper':stopped['tight_complete_with_reference_bytes_upper'],
            'same_candidate_fits_tight_forecast':stopped['whole_pivot_live_bound']<=150000 and stopped['tight_complete_output_bytes_upper']<=32*1024*1024,
            'no_row_update_Gaussian_field_pivot_series_or_new_source_computed':True,
            'coefficient_source_statistics':coefficients,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [reportp,matrixp,folder/'stopped-candidate.json',ROOT/f'code/{Path(__file__).name}',ROOT/f'research/notes/{ID}.md']},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256','coefficient_source_statistics')}))
if __name__=='__main__':run()
