#!/usr/bin/env python3
"""Source-only next-ONE row-stream gate from exact generation38."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json,resource,signal,time
from s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63 import stat,supp,encoded_bytes,common_denominator_primitive_bytes
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM38_SOURCE_GATE_20261104_v63'
PREV='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PIVOT_20261104_v63'
OLDINDEX='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PREFLIGHT_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def index_row(path,row,fileid):
    cells=[]
    with path.open('rb')as f:
        first=f.read(65536);at=first.index(b'"entries":[')+len(b'"entries":[');f.seek(at)
        buf=bytearray();depth=0;quoted=False;escaped=False;absolute=at;done=False
        while not done:
            block=f.read(65536);assert block
            for x in block:
                if not depth:
                    if x in b' \n\r\t,':absolute+=1;continue
                    if x==93:done=True;break
                    assert x==91;start=absolute;buf=bytearray()
                buf.append(x)
                if quoted:
                    if escaped:escaped=False
                    elif x==92:escaped=True
                    elif x==34:quoted=False
                elif x==34:quoted=True
                elif x==91:depth+=1
                elif x==93:depth-=1
                absolute+=1
                if depth==0:
                    col,p=json.loads(buf);local=buf.index(b'{');payload=bytes(buf[local:-1]);n=len(p['terms'])
                    assert n and col!=24 and len({tuple(z[1:])for z in p['terms']})==n
                    c={'file_index':fileid,'row':row,'column':col,'payload_offset':start+local,
                       'payload_bytes':len(payload),'payload_sha256':hashlib.sha256(payload).hexdigest(),'terms':n}
                    if n==1:
                        c['monomial_QQ']=str(Q(p['content'])*int(p['terms'][0][0]))
                        c['monomial_exponent']=[p['shift'][v]+p['terms'][0][v+1]for v in range(5)]
                    cells.append(c);del p,payload;buf=bytearray()
    return cells
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    pd=ROOT/'reports'/PREV;pr=pd/'report.json';mp=pd/'generation-1.json';ledger=pd/'final-RAW-ledger.json'
    assert sha(pr)=='72764ac409a435dd69689eaa2668c4c5ab2e1a3531d0a8019d1199eba09175a6'
    assert sha(mp)=='e2a4db6346a784e150d77c50fa7f801baba476a7d35422c4112960f874d2716a'
    r=read(pr);assert r['all_passed']and r['completed_total_on_reordered_path']==38
    m=read(mp);assert m['shape']==[198,45]and m['logical_total_terms']==174569
    oldip=ROOT/'reports'/OLDINDEX/'literal-file-index.json';old=read(oldip)
    assert sha(oldip)==m['index_sha256']
    oldfile=ROOT/old['source_file'];assert sha(oldfile)==old['source_file_sha256']
    files=[{'path':old['source_file'],'sha256':old['source_file_sha256'],'bytes':oldfile.stat().st_size}]
    cells=[];rowcells={};colcells={};maxcell=0
    for row in m['rows']:
        i=row['original_row'];assert i!=222
        if'literal_row_file'in row:
            p=pd/row['literal_row_file'];assert sha(p)==row['sha256'];fid=len(files)
            files.append({'path':str(p.relative_to(ROOT)),'sha256':row['sha256'],'bytes':p.stat().st_size})
            new=index_row(p,i,fid);assert sum(c['terms']for c in new)==row['terms']
        else:
            new=[{'file_index':0,**old['all_literal_cells'][k]}for k in row['indexed_source_cell_indices']]
            assert all(c['row']==i and c['column']!=24 for c in new)
        for c in new:
            assert c['column']in m['column_original_ids'];k=len(cells);cells.append(c)
            rowcells.setdefault(i,[]).append(k);colcells.setdefault(c['column'],[]).append(k);maxcell=max(maxcell,c['terms'])
    assert sum(c['terms']for c in cells)==174569
    assert len({(c['row'],c['column'])for c in cells})==len(cells)
    index={'generation':38,'source_manifest':str(mp.relative_to(ROOT)),'source_manifest_sha256':sha(mp),'files':files,
           'all_literal_cells':cells,'row_cell_indices':rowcells,'column_cell_indices':colcells,
           'all_original_row_ids':[v['original_row']for v in m['rows']],'all_original_column_ids':m['column_original_ids'],
           'logical_shape':[198,45],'logical_coefficient_terms':174569,'max_cell_terms':maxcell}
    ip=out/'mixed-literal-file-index.json';dump(ip,index)
    rowterms={i:sum(cells[k]['terms']for k in rowcells.get(i,[]))for i in index['all_original_row_ids']}
    rowbytes={i:sum(cells[k]['payload_bytes']for k in rowcells.get(i,[]))for i in rowterms}
    monos=[k for k,c in enumerate(cells)if'monomial_QQ'in c]
    def score(k):
        c=cells[k];targets=[cells[z]['row']for z in colcells[c['column']]if cells[z]['row']!=c['row']]
        return(rowterms[c['row']],sum(rowbytes[i]for i in targets),c['row'],c['column'])
    shortlist=sorted(monos,key=score)[:16];dump(out/'metadata-shortlist.json',[{'cell_index':k,**cells[k],'score':score(k)}for k in shortlist])
    del old,m,r
    io=0
    def loadrow(i):
        nonlocal io
        a={}
        # Open only the file needed by the current cell, never all polynomials.
        for k in rowcells.get(i,[]):
            c=cells[k];p=ROOT/files[c['file_index']]['path']
            with p.open('rb')as f:f.seek(c['payload_offset']);b=f.read(c['payload_bytes'])
            io+=len(b);assert hashlib.sha256(b).hexdigest()==c['payload_sha256'];a[c['column']]=json.loads(b)
        return a
    base=read(ledger)['TOTAL_RAW_bytes'];assert base==12681245
    source_bytes=Path(__file__).stat().st_size+(ROOT/f'research/notes/{ID}.md').stat().st_size
    results=[];chosen=None
    for ordinal,ci in enumerate(shortlist):
        if time.monotonic()-start>43:break
        c=cells[ci];ri,pc=c['row'],c['column'];prow=loadrow(ri);R=rowterms[ri]
        unit=Q(c['monomial_QQ']);pe=c['monomial_exponent'];pcols=[j for j in index['all_original_column_ids']if j!=pc and j in prow]
        targets=[cells[k]['row']for k in colcells[pc]if cells[k]['row']!=ri]
        scalar_meta=len(monos);rowforecasts=[];owned=0;supportpeak=0;parsepeak=scalar_meta+R+maxcell*2+96
        logical=174569;newrows=0;pairs=0;normalized=512
        for p in prow.values():
            sh=[p['shift'][v]-pe[v]for v in range(5)]
            normalized+=encoded_bytes(p)-len(p['content'])+len(str(Q(p['content'])/unit))+encoded_bytes(sh)-encoded_bytes(p['shift'])+32
        for i in targets:
            if time.monotonic()-start>45:raise RuntimeError('source45s stop; no arithmetic')
            cur=loadrow(i);factor=cur[pc];nf,cf,mf,_,_=stat(factor);mf=(1<<mf.bit_length())-1;fs=supp(factor)
            resident=scalar_meta+R+rowterms[i];parsepeak=max(parsepeak,resident+2*max((len(p['terms'])for p in cur.values()),default=0)+96)
            rowbound=512+sum(encoded_bytes(p)+32 for j,p in cur.items()if j!=pc);updates=[]
            for j in pcols:
                oldp=cur.get(j);a,ca,ma,_,_=stat(oldp);ma=(1<<ma.bit_length())-1
                np,cp,mpm,_,_=stat(prow[j]);mpm=(1<<mpm.bit_length())-1;cp/=unit
                ps=[tuple(e[v]-pe[v]for v in range(5))for e in supp(prow[j])]
                prod={tuple(f[v]+p[v]for v in range(5))for f in fs for p in ps};union=set(supp(oldp))|prod
                k,b=len(prod),len(union);peak=resident+2*k+3*b+96;owned=max(owned,peak)
                supportpeak=max(supportpeak,resident+len(fs)+len(ps)+a+k+b+96)
                cc=cf*cp;D=ca.denominator*cc.denominator//gcd(ca.denominator,cc.denominator)
                M=abs(ca.numerator)*(D//ca.denominator)*ma+abs(cc.numerator)*(D//cc.denominator)*min(nf,np)*mf*mpm
                osh=oldp['shift']if oldp is not None else[0]*5
                bsh=[min(osh[v],factor['shift'][v]+prow[j]['shift'][v]-pe[v])for v in range(5)]
                nb=common_denominator_primitive_bytes(union,bsh,M,D)+b+1
                rowbound+=nb-(encoded_bytes(oldp)if oldp is not None else 0)
                updates.append([j,a,k,b,resident,peak,M.bit_length(),D.bit_length()]);resident+=b-a;logical+=b-a;pairs+=nf*np
            resident-=nf;logical-=nf;newrows+=rowbound
            rowforecasts.append({'row':i,'source_terms':rowterms[i],'new_row_terms_upper':resident-scalar_meta-R,
                                 'new_immutable_row_file_bytes_upper':rowbound,'updates':updates,
                                 'factor_reference_cell_index':next(k for k in rowcells[i]if cells[k]['column']==pc)})
            del cur
        logical-=R;owned=max(owned,parsepeak,supportpeak)
        nu=sum(len(f['updates'])for f in rowforecasts)
        idxslots=6*len(cells)+len(cells)+len(rowcells)+sum(len(index[k])for k in('all_original_row_ids','all_original_column_ids'))+5*len(targets)+8*nu+3*len(targets)+4*len(files)+512
        # The old complete history stays charged. New immutable gate/index,
        # code and ONE future generation are added; no previous file removed.
        currentgate=sum(p.stat().st_size for p in out.iterdir()if p.is_file())
        reserve=3*1024*1024
        raw=base+currentgate+source_bytes+normalized+newrows+reserve
        record={'ordinal':ordinal,'source_indexed_pivot':c,'pivot_row_terms':R,'target_rows':len(targets),
                'charged_candidate_scalar_metadata':scalar_meta,'owned_algebraic_and_scratch_bound':owned,
                'index_plan_receipt_integer_slots_upper':idxslots,'whole_owned_scalar_and_scratch_upper':idxslots+owned,
                'parse_phase_bound':parsepeak,'support_phase_bound':supportpeak,'logical_after_support_union_upper':logical,
                'normalized_row_RAW_bound':normalized,'all_new_rows_RAW_bound':newrows,
                'RAW_history_base':base,'RAW_current_gate_index_source_bytes':currentgate+source_bytes,
                'RAW_future_worker_control_executable_all_receipts_manifest_and_extra_gate_reserve':reserve,
                'RAW_TOTAL_bound':raw,'integer_coefficient_pair_operations_upper':pairs,
                'all_row_forecasts':rowforecasts,'FULL_safe':idxslots+owned<=150000 and raw<=32*1024*1024,
                'source_only_no_coefficient_update':True}
        rp2=out/f'candidate-{ordinal:02d}-full-gate.json';dump(rp2,record);results.append({k:v for k,v in record.items()if k!='all_row_forecasts'})
        if record['FULL_safe']:chosen=(ordinal,record);break
        del prow,rowforecasts
    dump(out/'completed-forecasts-summary.json',results)
    report={'id':ID,'all_passed':True,'source_only':True,'source_generation38_not_replayed':True,
            'all_current_literal_cells_hash_bound':True,'logical_input_terms':174569,'source_shape':[198,45],
            'completed_candidate_forecasts':len(results),'shortlist_size':len(shortlist),'not_all_unit_candidates_claimed_examined':True,
            'FULL_safe_found':chosen is not None,'selected_original_row_column':[chosen[1]['source_indexed_pivot']['row'],chosen[1]['source_indexed_pivot']['column']]if chosen else None,
            'selected_gate_filename':f'candidate-{chosen[0]:02d}-full-gate.json'if chosen else None,
            'no_arithmetic_pivot_gcd_point_prime_or_nonmonomial_inverse':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[pr,mp,ledger,oldip,Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63.py']},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
            'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump(out/'report.json',report);signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
