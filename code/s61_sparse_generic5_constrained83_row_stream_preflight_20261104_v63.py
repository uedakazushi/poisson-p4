#!/usr/bin/env python3
"""SOURCE ONLY indexed row-stream representation; no coefficient update."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json,resource,signal,time
from s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63 import stat,supp,encoded_bytes,common_denominator_primitive_bytes
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PREFLIGHT_20261104_v63'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_FRONTIER_SAFE_CONTINUATION_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def index_cells(path):
    """Read raw bytes in chunks and parse only ONE literal coefficient cell."""
    cells=[];rowcells={};colcells={};max_cell_terms=0;max_buffer=0
    with path.open('rb')as f:
        prefix=f.read(65536);token=b'"sparse_entries":[';where=prefix.index(token)
        before=prefix[:where];assert before.endswith(b',')
        meta=json.loads(before[:-1]+b'}');start=where+len(token);f.seek(start)
        depth=0;quoted=False;escaped=False;buf=bytearray();entry_start=None;absolute=start;array_end=None
        while array_end is None:
            block=f.read(65536);assert block
            for x in block:
                if not depth:
                    if x in b' \r\n\t,':absolute+=1;continue
                    if x==93:array_end=absolute;break
                    assert x==91;entry_start=absolute;buf=bytearray();depth=0;quoted=False;escaped=False
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
                    assert buf.endswith(b'}]')
                    row,col,p=json.loads(buf);assert row in meta['row_original_ids']and col in meta['column_original_ids']
                    n=len(p['terms']);assert n and Q(p['content']) and len(p['shift'])==5
                    assert all(len(t)==6 and int(t[0]) and min(t[1:])>=0 for t in p['terms'])
                    assert len({tuple(t[1:])for t in p['terms']})==n
                    payload_local=buf.index(b'{');payload=bytes(buf[payload_local:-1])
                    info={'row':row,'column':col,'payload_offset':entry_start+payload_local,
                          'payload_bytes':len(payload),'payload_sha256':hashlib.sha256(payload).hexdigest(),'terms':n}
                    if n==1:
                        info['monomial_QQ']=str(Q(p['content'])*int(p['terms'][0][0]))
                        info['monomial_exponent']=[p['shift'][v]+p['terms'][0][v+1]for v in range(5)]
                    pos=len(cells);cells.append(info);rowcells.setdefault(row,[]).append(pos);colcells.setdefault(col,[]).append(pos)
                    max_cell_terms=max(max_cell_terms,n);max_buffer=max(max_buffer,len(buf))
                    del p,payload;buf=bytearray()
        f.seek(array_end+1);tail=f.read();assert tail.startswith(b',');summary=json.loads(b'{'+tail[1:])
    assert sum(c['terms']for c in cells)==summary['total_terms']==147370
    assert len(cells)==summary['nonzero_entries']==5358
    assert len({(c['row'],c['column'])for c in cells})==len(cells)
    return meta,cells,rowcells,colcells,max_cell_terms,max_buffer
def run():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/SOURCE;rp=folder/'report.json';r=read(rp)
    assert r['all_passed'] and sha(rp)=='1bc5d3f4867426f286ebe8c2d91f79a910dc3a97d4a693a765612d8291a61ca4'
    mp=folder/'residual-matrix.json';sourcehash=sha(mp);assert sourcehash==r['evidence_sha256'][mp.name]
    meta,cells,rowcells,colcells,maxcell,maxbuf=index_cells(mp)
    assert meta['shape']==[199,46]
    rowterms={i:sum(cells[k]['terms']for k in rowcells.get(i,[]))for i in meta['row_original_ids']}
    rowbytes={i:sum(cells[k]['payload_bytes']for k in rowcells.get(i,[]))for i in meta['row_original_ids']}
    monos=[k for k,c in enumerate(cells)if'monomial_QQ'in c]
    # Metadata-only deterministic choice: small resident pivot row, then
    # target-read volume, source row/column. No matrix arithmetic/rank search.
    def score(k):
        c=cells[k];i,j=c['row'],c['column']
        target=[cells[z]['row']for z in colcells[j]if cells[z]['row']!=i]
        return(rowterms[i],sum(rowbytes[x]for x in target),i,j)
    chosen=min(monos,key=score);pivot=cells[chosen];row,col=pivot['row'],pivot['column']
    index={'source_file':str(mp.relative_to(ROOT)),'source_file_sha256':sourcehash,'logical_shape':[199,46],
           'logical_coefficient_terms':147370,'all_literal_cells':cells,
           'row_cell_indices':rowcells,'column_cell_indices':colcells,
           'all_original_row_ids':meta['row_original_ids'],'all_original_column_ids':meta['column_original_ids'],
           'max_single_cell_terms_parsed_during_index':maxcell,'max_entry_raw_buffer_bytes':maxbuf,
           'no_full_coefficient_array_materialized':True}
    ip=out/'literal-file-index.json';ip.write_text(json.dumps(index,indent=2)+'\n')
    io_read=0;owner_peak=maxcell+len(monos)+96
    def load_row(i):
        nonlocal io_read,owner_peak
        a={}
        with mp.open('rb')as f:
            for k in rowcells.get(i,[]):
                c=cells[k];f.seek(c['payload_offset']);b=f.read(c['payload_bytes']);io_read+=len(b)
                assert hashlib.sha256(b).hexdigest()==c['payload_sha256']
                p=json.loads(b);assert len(p['terms'])==c['terms'];a[c['column']]=p
        return a
    prow=load_row(row);R=rowterms[row];unit=Q(pivot['monomial_QQ']);pe=tuple(pivot['monomial_exponent'])
    pcols=[j for j in meta['column_original_ids']if j!=col and j in prow]
    targets=[cells[k]['row']for k in colcells[col]if cells[k]['row']!=row]
    scalar_meta=len(monos);logical=147370;resident_peak=0;resident_support_peak=0
    normalized_bytes=512
    for j,p in prow.items():
        sh=[p['shift'][v]-pe[v]for v in range(5)]
        normalized_bytes+=encoded_bytes(p)-len(p['content'])+len(str(Q(p['content'])/unit))+encoded_bytes(sh)-encoded_bytes(p['shift'])+32
    newrows_bytes=0;receipt_bytes=1024;source_row_read_bytes=0;operation_pair_bound=0;rowforecasts=[]
    known_global_old_peak=147370
    for i in targets:
        if time.monotonic()-started>45:raise RuntimeError('source_preflight45s guard; no arithmetic')
        current=load_row(i);source_row_read_bytes+=rowbytes[i];factor=current[col];nf,cf,mf,_,_=stat(factor);mf=(1<<mf.bit_length())-1;fs=supp(factor)
        resident=scalar_meta+R+rowterms[i];owner_peak=max(owner_peak,resident+96)
        rowbound=512+sum(encoded_bytes(p)+32 for j,p in current.items()if j!=col);row_updates=[]
        for j in pcols:
            old=current.get(j);a,ca,ma,_,_=stat(old);ma=(1<<ma.bit_length())-1
            np,cp,mpm,_,_=stat(prow[j]);mpm=(1<<mpm.bit_length())-1;cp/=unit
            ps=[tuple(e[v]-pe[v]for v in range(5))for e in supp(prow[j])]
            product={tuple(f[v]+p[v]for v in range(5))for f in fs for p in ps};union=set(supp(old))|product
            k=len(product);b=len(union);scratch=len(fs)+len(ps)+a+k+b
            peak=resident+2*k+3*b+96
            resident_peak=max(resident_peak,peak)
            resident_support_peak=max(resident_support_peak,scalar_meta+R+rowterms[i]+scratch+96)
            assert peak<=150000 and resident_support_peak<=150000
            known_global_old_peak=max(known_global_old_peak,logical+2*k+3*b+64)
            cprod=cf*cp;D=ca.denominator*cprod.denominator//gcd(ca.denominator,cprod.denominator)
            M=abs(ca.numerator)*(D//ca.denominator)*ma+abs(cprod.numerator)*(D//cprod.denominator)*min(nf,np)*mf*mpm
            oldshift=old['shift']if old is not None else[0]*5
            base=[min(oldshift[v],factor['shift'][v]+prow[j]['shift'][v]-pe[v])for v in range(5)]
            newbytes=common_denominator_primitive_bytes(union,base,M,D)+b+1
            rowbound+=newbytes-(encoded_bytes(old)if old is not None else 0)
            row_updates.append([j,a,k,b,resident,peak,M.bit_length(),D.bit_length()])
            resident+=b-a;logical+=b-a;operation_pair_bound+=nf*np
        resident-=nf;logical-=nf
        newrows_bytes+=rowbound
        receipt_bytes+=1024+256*len(row_updates)
        rowforecasts.append({'row':i,'source_terms':rowterms[i],'new_row_terms_upper':resident-scalar_meta-R,
                             'new_immutable_row_file_bytes_upper':rowbound,'updates':row_updates,
                             'factor_reference_cell_index':next(k for k in rowcells[i]if cells[k]['column']==col)})
        del current
    logical-=R
    input_bytes=mp.stat().st_size;index_bytes=ip.stat().st_size
    rawtotal=input_bytes+index_bytes+normalized_bytes+newrows_bytes+receipt_bytes+2*1024*1024
    owned=max(owner_peak,resident_peak,resident_support_peak)
    assert owned<=150000 and rawtotal<=32*1024*1024
    full={'source_indexed_pivot':pivot,'pivot_selection':'minimum pivot-row terms, then target read bytes,source row,column',
          'logical_input_terms':147370,'logical_after_support_union_upper':logical,
          'all_resident_copy_forecast_peak_for_same_path':known_global_old_peak,
          'new_row_resident_peak':resident_peak,'source_support_forecast_peak':resident_support_peak,
          'actual_source_parsed_owner_peak':owner_peak,'charged_candidate_scalar_metadata':scalar_meta,
          'owned_live_with_all_scratch_bound':owned,'pivot_row_terms':R,'target_rows':len(targets),
          'max_input_row_terms':max(rowterms.values()),'max_input_cell_terms':maxcell,
          'RAW_ledger':{'frozen_input':input_bytes,'new_literal_index':index_bytes,'new_normalized_pivot_row':normalized_bytes,
                        'all_new_immutable_updated_rows':newrows_bytes,'all_receipts_and_forecasts':receipt_bytes,
                        'manifest_cursor_other_metadata_reserve':2*1024*1024,'TOTAL':rawtotal},
          'I_O':{'source_input_bytes':input_bytes,'source_row_payload_bytes_actually_read':io_read,
                 'first_arithmetic_read_payload_upper':rowbytes[row]+source_row_read_bytes,
                 'new_generation_write_bytes_upper':normalized_bytes+newrows_bytes+receipt_bytes,
                 'source_hash_plus_index_scan_read_upper':2*input_bytes+65536},
          'integer_coefficient_pair_operations_upper_for_one_product_pass':operation_pair_bound,
          'all_row_forecasts':rowforecasts,'no_actual_coefficient_update':True,
          'gate_passed_owned_and_RAW':True}
    ep=out/'first-indexed-row-stream-gate.json';ep.write_text(json.dumps(full,indent=2)+'\n')
    helper=ROOT/'code/s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63.py'
    policy={'resident_definition':'Only the scalar candidate index, one full current row, one fixed pivot row and one entry update scratch. All other literal coefficients remain in immutable source files and are represented by verified offset/length/hash references.',
            'not_a_cap_increase':'Logical full-matrix term count is separately reported and may exceed150000. OWNED resident coefficients plus every arithmetic/support scratch still must be<=150000 before every row update.',
            'exact_arithmetic_kernel':'Reuse the already-verified Laurent product/sub/back-add kernel. Fused new map is unnecessary for this first row-bound gate; no unreviewed fused arithmetic is assumed.',
            'logical_preservation':'Write each changed full row immutably and verify every back-add. Unchanged rows reuse exact indexed input references. Publish the next complete manifest only after ALL row files/receipts are verified; remove only the proved unit row/column.',
            'storage':'Charge frozen input once plus EVERY new row version, normalization receipt, index, manifest and operation record generated by the job. Keep all generated versions; no deletion or relabeling as external to evade32MiB.',
            'failure_atomicity':'Until the whole pivot commits, old manifest remains current. Preserve completed pending row files/identities on timeout/failure without claiming a completed new pivot.',
            'temporary_release':'Release parsed rows/cells from RAM after immutable write. Disk versions may be removed only with a separately verified exact reconstructing program; no deletion is needed or assumed here.',
            'time':'Source gate hard60s/2GB. A future arithmetic job needs its own deadline/soft guard and atomic-generation save. The byte-volume/term-pair bounds are not a promise of operating-system throughput.',
            'no_arithmetic_GO_here':True,'no_nonmonomial_point_prime_gcd_or_series':True}
    (out/'representation-and-atomic-generation-policy.json').write_text(json.dumps(policy,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'source_only':True,'all5358_literal_cells_indexed_and_hash_bound':True,
            'all147370_logical_coefficients_preserved_externally':True,
            'first_candidate':[row,col],'owned_live_and_scratch_bound':owned,
            'logical_new_support_terms_upper':logical,'RAW_input_plus_all_new_data_upper':rawtotal,
            'no_actual_row_update_pivot_gcd_series_or_nonmonomial_inverse':True,'arithmetic_requires_parent_review':True,
            'input_sha256':{str(p.relative_to(ROOT)):(sourcehash if p==mp else sha(p))for p in[rp,mp,Path(__file__),helper,ROOT/f'research/notes/{ID}.md']},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert input_bytes+sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
if __name__=='__main__':run()
