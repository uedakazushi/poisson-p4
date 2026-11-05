#!/usr/bin/env python3
"""Next ONE source-only unit forecast on the exact compact generation39."""
from pathlib import Path
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from fractions import Fraction as Q
from math import gcd
import hashlib,json,resource,signal,time
from s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63 import stat,supp,encoded_bytes,common_denominator_primitive_bytes
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT43_UNIT_GATE_20261105_v64'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM43_MIXED_BINARY_COMPACTION_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);sd=ROOT/'reports'/SOURCE
    sr=sd/'report.json';ip=sd/'literal-file-index.json';matrix=sd/'current-matrix.bin';receipt=sd/'all-cell-row-roundtrip.json'
    assert sha(sr)=='64fd53c47d02b160aef71fe31e899a0d1872d7f0fb6e58ea6692cbd8d28dd046'
    r=read(sr);assert r['all_passed']and r['exact_current_generation_compacted']==43
    for p in(ip,matrix,receipt):assert sha(p)==r['evidence_sha256'][p.name]
    rr=read(receipt);assert rr['all_cells_equal']and len(rr['all_rows'])==193
    index=read(ip);assert index['logical_shape']==[193,40]and index['logical_coefficient_terms']==400727
    cells=index['all_literal_cells'];rows=index['all_original_row_ids'];cols=index['all_original_column_ids'];rc=index['row_cell_indices'];cc=index['column_cell_indices']
    rt={i:sum(cells[k]['terms']for k in rc.get(str(i),[]))for i in rows};rb={i:sum(cells[k]['payload_bytes']for k in rc.get(str(i),[]))for i in rows}
    monos=[k for k,c in enumerate(cells)if'monomial_QQ'in c]
    def score(k):
        c=cells[k];return(rt[c['row']],sum(rb[cells[z]['row']]for z in cc[str(c['column'])]if cells[z]['row']!=c['row']),c['row'],c['column'])
    shortlist=sorted(monos,key=score)[:16];dump(out/'metadata-shortlist.json',[{'cell_index':k,**cells[k],'score':score(k)}for k in shortlist])
    io=0
    def loadrow(i):
        nonlocal io
        z={}
        with matrix.open('rb')as f:
            for k in rc.get(str(i),[]):
                c=cells[k];f.seek(c['payload_offset']);b=f.read(c['payload_bytes']);io+=len(b)
                assert hashlib.sha256(b).hexdigest()==c['payload_sha256'];z[c['column']]=decode(b)
        return z
    helper=ROOT/'code/s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63.py'
    inputs=[sr,ip,matrix,receipt,helper,Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_primitive_binary_native_20261105_v64.hpp']
    active=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in inputs]
    base=sum(p['bytes']for p in active);completed=[];chosen=None
    for ordinal,k in enumerate(shortlist):
        if time.monotonic()-start>43:break
        c=cells[k];ri,pc=c['row'],c['column'];prow=loadrow(ri);R=rt[ri];unit=Q(c['monomial_QQ']);pe=c['monomial_exponent']
        pcols=[j for j in cols if j!=pc and j in prow];targets=[cells[z]['row']for z in cc[str(pc)]if cells[z]['row']!=ri]
        sm=len(monos);owned=sm+R+4*max(v['terms']for v in cells)+256;supportpeak=0;logical=400727;newbytes=0;plans=[];normbytes=512;pairs=0
        for p in prow.values():
            sh=[p['shift'][v]-pe[v]for v in range(5)]
            normbytes+=encoded_bytes(p)-len(p['content'])+len(str(Q(p['content'])/unit))+encoded_bytes(sh)-encoded_bytes(p['shift'])+32
        for i in targets:
            if time.monotonic()-start>45:raise RuntimeError('source45s stop; no arithmetic')
            cur=loadrow(i);factor=cur[pc];nf,cf,mf,_,_=stat(factor);mf=(1<<mf.bit_length())-1;fs=supp(factor)
            resident=sm+R+rt[i];owned=max(owned,resident+4*max(len(p['terms'])for p in cur.values())+256)
            rowbound=512+sum(encoded_bytes(p)+32 for j,p in cur.items()if j!=pc);updates=[]
            for j in pcols:
                old=cur.get(j);a,ca,ma,_,_=stat(old);ma=(1<<ma.bit_length())-1
                np,cp,mp,_,_=stat(prow[j]);mp=(1<<mp.bit_length())-1;cp/=unit
                ps=[tuple(e[v]-pe[v]for v in range(5))for e in supp(prow[j])]
                prod={tuple(f[v]+p[v]for v in range(5))for f in fs for p in ps};union=set(supp(old))|prod;kk,b=len(prod),len(union)
                peak=resident+2*kk+3*b+96;owned=max(owned,peak);supportpeak=max(supportpeak,resident+len(fs)+len(ps)+a+kk+b+96)
                cfcp=cf*cp;D=ca.denominator*cfcp.denominator//gcd(ca.denominator,cfcp.denominator)
                M=abs(ca.numerator)*(D//ca.denominator)*ma+abs(cfcp.numerator)*(D//cfcp.denominator)*min(nf,np)*mf*mp
                os=old['shift']if old is not None else[0]*5
                bs=[min(os[v],factor['shift'][v]+prow[j]['shift'][v]-pe[v])for v in range(5)]
                rowbound+=common_denominator_primitive_bytes(union,bs,M,D)+b+1-(encoded_bytes(old)if old is not None else 0)
                updates.append([j,a,kk,b,resident,peak,M.bit_length(),D.bit_length()]);resident+=b-a;logical+=b-a;pairs+=nf*np
            resident-=nf;logical-=nf;newbytes+=rowbound
            plans.append({'row':i,'source_terms':rt[i],'new_row_terms_upper':resident-sm-R,'new_immutable_row_file_bytes_upper':rowbound,
                          'factor_reference_cell_index':next(k for k in rc[str(i)]if cells[k]['column']==pc),'updates':updates})
            del cur
        logical-=R;owned=max(owned,supportpeak);nu=sum(len(p['updates'])for p in plans)
        metacount=5*len(cells)+len(cells)+len(rc)+len(rows)+len(cols)+5*len(plans)+8*nu+3*len(plans)+512
        reserve=3*1024*1024;raw=base+normbytes+newbytes+reserve
        gate={'ordinal':ordinal,'source_indexed_pivot':c,'pivot_row_terms':R,'target_rows':len(plans),'charged_candidate_scalar_metadata':sm,
              'owned_algebraic_and_scratch_bound':owned,'index_plan_receipt_integer_slots_upper':metacount,
              'whole_owned_scalar_and_scratch_upper':metacount+owned,'logical_after_support_union_upper':logical,
              'normalized_row_RAW_bound':normbytes,'all_new_rows_RAW_bound':newbytes,
              'actual_READ_inputs':active,'actual_READ_input_RAW_bytes':base,'RAW_future_all_extra_gate_worker_code_control_executable_receipts_manifest_reserve':reserve,
              'RAW_TOTAL_bound':raw,'integer_coefficient_pair_operations_upper':pairs,'all_row_forecasts':plans,
              'FULL_safe':metacount+owned<=150000 and raw<=32*1024*1024,'source_only':True,
              'retired_history_retained_not_read_or_deleted':True}
        gp=out/f'candidate-{ordinal:02d}-full-gate.json';dump(gp,gate);completed.append({j:v for j,v in gate.items()if j not in('all_row_forecasts','actual_READ_inputs')})
        if gate['FULL_safe']:chosen=(ordinal,gate);break
    dump(out/'completed-forecasts-summary.json',completed)
    report={'id':ID,'all_passed':True,'source_only':True,'generation43_not_replayed':True,'shape':[193,40],'logical_input_terms':400727,
            'FULL_safe_found':chosen is not None,'selected_original_row_column':[chosen[1]['source_indexed_pivot']['row'],chosen[1]['source_indexed_pivot']['column']]if chosen else None,
            'selected_gate_filename':f'candidate-{chosen[0]:02d}-full-gate.json'if chosen else None,'completed_forecasts':len(completed),'shortlist_size':16,
            'not_global_minimum_or_all_units_unsafe_claim':True,'no_pivot_field_GCD_point_or_prime':True,
            'input_sha256':{f['path']:f['sha256']for f in active},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
            'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
