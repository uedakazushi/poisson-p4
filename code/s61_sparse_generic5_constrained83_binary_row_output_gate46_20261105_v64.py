#!/usr/bin/env python3
"""One changed binary-row representation forecast, no coefficient update."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,uint,encode,decode
from s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63 import stat,supp
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_BINARY_ROW_OUTPUT_GATE46_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def zbytes(n):
    k=(int(n).bit_length()+7)//8;return len(uint(2*k+1))+k
def binbound(union,base,M,D):
    if not union or not M:return 0
    return zbytes(M)+zbytes(D)+sum(len(uint(2*x if x>=0 else-2*x-1))for x in base)+len(uint(len(union)))+len(union)*zbytes(M)+sum(len(uint(e[v]-base[v]))for e in union for v in range(5))
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    stop=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT46_SOURCE_STOP_20261105_v64/report.json'
    frontier=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT46_UNIT_FRONTIER7_8_20261105_v64/report.json'
    assert sha(stop)=='ace48d24b242e6c5e309511aceff5de52543c227185a4a230cd23cfdad026764'
    assert sha(frontier)=='4ca50803c1b10151b7516e13d0a4d66ce69334938a0688fd312605ec3fc09811'
    f=read(frontier);assert f['all_passed']and not f['FULL_safe_found']
    gp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT46_UNIT_GATE_20261105_v64/candidate-06-full-gate.json';g=read(gp)
    assert g['ordinal']==6 and g['source_indexed_pivot']['row']==76 and g['source_indexed_pivot']['column']==66
    assert sha(gp)==read(stop)['input_sha256'][str(gp.relative_to(ROOT))]
    ip=next(ROOT/x['path']for x in g['actual_READ_inputs']if x['path'].endswith('/literal-file-index.json'));index=read(ip);matrix=ROOT/index['source_file']
    assert index['encoding']=='primitive-int-exponent-v1'and sha(matrix)==index['source_file_sha256']
    cells=index['all_literal_cells'];rc=index['row_cell_indices'];cols=index['all_original_column_ids'];io=0
    def loadrow(r):
        nonlocal io
        z={}
        with matrix.open('rb')as f:
            for k in rc.get(str(r),[]):
                c=cells[k];f.seek(c['payload_offset']);p=f.read(c['payload_bytes']);io+=len(p);assert hashlib.sha256(p).hexdigest()==c['payload_sha256'];z[c['column']]=decode(p)
        return z
    pivot=g['source_indexed_pivot'];pr,pc=pivot['row'],pivot['column'];q=Q(pivot['monomial_QQ']);pe=pivot['monomial_exponent'];prow=loadrow(pr);R=g['pivot_row_terms']
    assert sum(len(p['terms'])for p in prow.values())==R
    plans=[];totalbin=0;totalidx=0;maxwriter=0;nbounds=0;maxcell=0
    for old in g['all_row_forecasts']:
        assert time.monotonic()-start<48
        r=old['row'];cur=loadrow(r);factor=cur[pc];nf,cf,mf,_,_=stat(factor);mf=(1<<mf.bit_length())-1;fs=supp(factor);bounds={}
        for j,p in cur.items():
            if j==pc:continue
            n=len(p['terms']);ex=supp(p);mono=str(Q(p['content'])*int(p['terms'][0][0]))if n==1 else''
            bounds[j]={'column':j,'binary_upper':len(encode(p)),'terms_upper':n,'QQ_chars_upper':len(mono),'exponent_chars_upper':sum(max(len(str(e[v]))for e in ex)for v in range(5))}
        for u in old['updates']:
            j,a,k,b,resident,peak,Mbits,Dbits=u;p=cur.get(j);aa,ca,ma,_,_=stat(p);ma=(1<<ma.bit_length())-1
            np,cp,mp,_,_=stat(prow[j]);mp=(1<<mp.bit_length())-1;cp/=q
            ps=[tuple(e[v]-pe[v]for v in range(5))for e in supp(prow[j])]
            prod={tuple(x[v]+y[v]for v in range(5))for x in fs for y in ps};union=set(supp(p))|prod
            assert(aa,len(prod),len(union))==(a,k,b)
            cprod=cf*cp;D=ca.denominator*cprod.denominator//gcd(ca.denominator,cprod.denominator)
            M=abs(ca.numerator)*(D//ca.denominator)*ma+abs(cprod.numerator)*(D//cprod.denominator)*min(nf,np)*mf*mp
            assert(M.bit_length(),D.bit_length())==(Mbits,Dbits)
            os=p['shift']if p is not None else[0]*5;base=[min(os[v],factor['shift'][v]+prow[j]['shift'][v]-pe[v])for v in range(5)]
            bounds[j]={'column':j,'binary_upper':binbound(union,base,M,D),'terms_upper':b,'QQ_chars_upper':len(str(M))+len(str(D))+4,
                       'exponent_chars_upper':sum(max(len(str(e[v]))for e in union)for v in range(5))}
        ordered=[bounds[j]for j in cols if j in bounds];nupper=sum(x['terms_upper']for x in ordered);assert nupper==old['new_row_terms_upper']
        binary=len(MAGIC)+sum(x['binary_upper']for x in ordered);idx=512
        for x in ordered:
            x['index_entry_upper']=128+len(str(binary))+len(str(x['binary_upper']))+len(str(x['terms_upper']))+x['QQ_chars_upper']+x['exponent_chars_upper']
            idx+=x['index_entry_upper']+1
        cmax=max((x['terms_upper']for x in ordered),default=0);writer=g['charged_candidate_scalar_metadata']+R+nupper+4*cmax+256
        maxwriter=max(maxwriter,writer);maxcell=max(maxcell,cmax);totalbin+=binary;totalidx+=idx;nbounds+=len(ordered)
        plans.append({'row':r,'binary_file_upper':binary,'paired_index_upper':idx,'writer_owned_upper':writer,
                      'cells':[[x['column'],x['binary_upper'],x['terms_upper'],x['index_entry_upper']]for x in ordered]})
        del cur
    active={x['path']:x for x in g['actual_READ_inputs']}
    for p in(gp,stop,frontier,Path(__file__),ROOT/f'research/notes/{ID}.md'):
        active[str(p.relative_to(ROOT))]={'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}
    bp=out/'binary-row-writer-bounds.json';bp.write_bytes(enc({'source_gate':str(gp.relative_to(ROOT)),'source_gate_sha256':sha(gp),'rows':plans}))
    active[str(bp.relative_to(ROOT))]={'path':str(bp.relative_to(ROOT)),'sha256':sha(bp),'bytes':bp.stat().st_size}
    for x in active.values():assert sha(ROOT/x['path'])==x['sha256']
    meta=g['index_plan_receipt_integer_slots_upper']+4*nbounds+4*len(plans)+512
    alg=max(g['owned_algebraic_and_scratch_bound'],maxwriter);owned=meta+alg;inputbytes=sum(x['bytes']for x in active.values());aux=3145728
    total=inputbytes+totalbin+totalidx+g['normalized_row_RAW_bound']+aux
    proof={'all_passed':True,'source_only':True,'changed_representation':'New literal target rows use exact uncompressed primitive-int/exponent binary plus complete paired indices.',
           'source_generation':46,'candidate':[pr,pc],'unit_QQ':str(q),'unit_exponent':pe,'actual_READ_inputs':list(active.values()),'actual_READ_input_RAW_bytes':inputbytes,
           'source_text_gate':str(gp.relative_to(ROOT)),'source_text_gate_sha256':sha(gp),'binary_bounds_file':str(bp.relative_to(ROOT)),'binary_bounds_sha256':sha(bp),
           'new_binary_row_files_RAW_upper':totalbin,'new_paired_indices_RAW_upper':totalidx,'normalized_pivot_RAW_upper':g['normalized_row_RAW_bound'],
           'other_control_code_executable_events_receipts_manifests_upper':aux,'fresh_RAW_upper':total,
           'algebraic_and_writer_scratch_upper':alg,'integer_index_plan_writer_metadata_slots':meta,'whole_owned_upper':owned,
           'FULL_safe':owned<=150000 and total<=32*1024*1024,'no_coefficient_product_expansion_GCD_point_prime_or_pivot':True,
           'M_D_primitive_integer_bound':'common integer numerator magnitudes<=M and denominator D imply primitive content numerator<=M, denominator<=D and every primitive integer coefficient<=M. No content factor is assumed.',
           'writer_roundtrip':'Each produced cell is independently binary-decoded to a FLINT polynomial and compared to the computed QQ cell before atomic publication.'}
    (out/'whole-budget.json').write_bytes(enc(proof))
    report={'id':ID,'all_passed':True,'source_only':True,'FULL_safe_binary_row_pivot':proof['FULL_safe'],'candidate':[pr,pc],'fresh_RAW_upper':total,'owned_upper':owned,
            'binary_row_upper':totalbin,'paired_index_upper':totalidx,'rows':len(plans),'writer_scratch_upper':alg,'input_sha256':{k:x['sha256']for k,x in active.items()if k!=str(bp.relative_to(ROOT))},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start}
    (out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
