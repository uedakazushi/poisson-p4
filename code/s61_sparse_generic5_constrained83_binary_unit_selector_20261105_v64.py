#!/usr/bin/env python3
"""Direct binary-output unit source gate; exact known-content byte bounds."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,uint,encode,decode
from s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63 import stat,supp
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(numeric(x)for x in v)
    if isinstance(v,dict):return sum(numeric(x)for x in v.values())
    return 0
def v2(n):
    n=abs(n);return(n&-n).bit_length()-1 if n else 10**9
def zbytes(n):
    k=(int(n).bit_length()+7)//8;return len(uint(2*k+1))+k
def binbound(support,base,integerM,contentN,D):
    if not support or not integerM:return 0
    return zbytes(contentN)+zbytes(D)+sum(len(uint(2*x if x>=0 else-2*x-1))for x in base)+len(uint(len(support)))+len(support)*zbytes(integerM)+sum(len(uint(e[v]-base[v]))for e in support for v in range(5))
def main(a):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    id=f'S61_SPARSE_GENERIC5_CONSTRAINED83_BINARY_UNIT_SELECTOR{a.generation}_20261105_v64';out=ROOT/'reports'/id;out.mkdir(parents=True,exist_ok=False)
    sd=ROOT/'reports'/a.source;sr=sd/'report.json';assert sha(sr)==a.report_sha;r=read(sr);assert r['all_passed']and r['exact_current_generation_compacted']==a.generation
    ip=sd/r['current_index'];matrix=sd/r['current_snapshot'];receipt=sd/'all-cell-row-roundtrip.json'
    for p in(ip,matrix,receipt):assert sha(p)==r['evidence_sha256'][p.name]
    index=read(ip);assert index['encoding']=='primitive-int-exponent-v1';cells=index['all_literal_cells'];rc=index['row_cell_indices'];cols=index['all_original_column_ids'];rows=index['all_original_row_ids'];logical=index['logical_coefficient_terms']
    mono=[k for k,c in enumerate(cells)if'monomial_QQ'in c];rt={i:sum(cells[k]['terms']for k in rc.get(str(i),[]))for i in rows};rb={i:sum(cells[k]['payload_bytes']for k in rc.get(str(i),[]))for i in rows}
    def score(k):
        c=cells[k];return(rt[c['row']],sum(rb[cells[z]['row']]for z in index['column_cell_indices'][str(c['column'])]if cells[z]['row']!=c['row']),c['row'],c['column'])
    shortlist=[];seen=set()
    for k in sorted(mono,key=score):
        if cells[k]['column']in seen:continue
        seen.add(cells[k]['column']);shortlist.append(k)
        if len(shortlist)==16:break
    (out/'metadata-shortlist.json').write_bytes(enc([{'cell_index':k,**cells[k],'score':score(k)}for k in shortlist]))
    paths=[sr,ip,matrix,receipt,Path(__file__),ROOT/f'research/notes/{id}.md',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_primitive_binary_native_20261105_v64.hpp',ROOT/'code/s61_sparse_generic5_constrained83_tight_forecast_preflight_20261104_v63.py']
    def loadrow(i):
        z={}
        with matrix.open('rb')as f:
            for k in rc.get(str(i),[]):
                c=cells[k];f.seek(c['payload_offset']);p=f.read(c['payload_bytes']);assert hashlib.sha256(p).hexdigest()==c['payload_sha256'];z[c['column']]=decode(p)
        return z
    summaries=[];selected=None;sourcepeak=0
    for ordinal in range(a.start,min(len(shortlist),a.start+a.count)):
        if time.monotonic()-start>42:break
        c=cells[shortlist[ordinal]];pr,pc=c['row'],c['column'];q=Q(c['monomial_QQ']);pe=c['monomial_exponent'];prow=loadrow(pr);R=rt[pr];pcols=[j for j in cols if j!=pc and j in prow];sm=len(mono)
        normupper=512+sum(len(json.dumps(p,separators=(',',':')))-len(p['content'])+len(str(Q(p['content'])/q))+256 for p in prow.values())
        # 256 per cell covers all signed-shift integer lengths (source exponents
        # are machine integers); actual normalization is checked by the worker.
        targets=[cells[k]['row']for k in index['column_cell_indices'][str(pc)]if cells[k]['row']!=pr]
        supportplans=[];writerplans=[];proofrows=[];alg=sm+R+4*max(x['terms']for x in cells)+256;newlogical=logical;totalbin=0;totalidx=0;pairs=0
        for i in targets:
            if time.monotonic()-start>45:raise RuntimeError('source45s stop; complete prior candidates retained; no arithmetic')
            cur=loadrow(i);factor=cur[pc];nf,cf,mf,_,_=stat(factor);mf=(1<<mf.bit_length())-1;fs=supp(factor);resident=sm+R+rt[i]
            alg=max(alg,resident+4*max(len(p['terms'])for p in cur.values())+256);updates=[];bounds={};twoproofs=[]
            for j,p in cur.items():
                if j==pc:continue
                ex=supp(p);monoqq=str(Q(p['content'])*int(p['terms'][0][0]))if len(p['terms'])==1 else''
                bounds[j]=[len(encode(p)),len(p['terms']),len(monoqq),sum(max(len(str(e[v]))for e in ex)for v in range(5))]
            for j in pcols:
                old=cur.get(j);nold,ca,ma,_,_=stat(old);ma=(1<<ma.bit_length())-1;np,cp,mp,_,_=stat(prow[j]);mp=(1<<mp.bit_length())-1;cp/=q
                ps=[tuple(e[v]-pe[v]for v in range(5))for e in supp(prow[j])];prod={tuple(x[v]+y[v]for v in range(5))for x in fs for y in ps};union=set(supp(old))|prod;k,b=len(prod),len(union)
                peak=resident+2*k+3*b+96;alg=max(alg,peak,resident+len(fs)+len(ps)+nold+k+b+96)
                cc=cf*cp;D=ca.denominator*cc.denominator//gcd(ca.denominator,cc.denominator);A=ca.numerator*(D//ca.denominator);B=cc.numerator*(D//cc.denominator)
                M=abs(A)*ma+abs(B)*min(nf,np)*mf*mp
                # Guaranteed common power of two of the TWO KNOWN integer
                # amplitudes, read directly by bit operations; no polynomial
                # factor/GCD search and no characteristic specialization.
                h=min(v2(A),v2(B));assert h<10**9;dh=min(h,v2(D));I=M//(1<<h);N=M//(1<<dh);DD=D//(1<<dh)
                os=old['shift']if old is not None else[0]*5;base=[min(os[v],factor['shift'][v]+prow[j]['shift'][v]-pe[v])for v in range(5)]
                conservative=binbound(union,base,M,M,D);refined=binbound(union,base,I,N,DD)
                assert refined<=conservative
                bounds[j]=[refined,b,len(str(N))+len(str(DD))+4,sum(max(len(str(e[v]))for e in union)for v in range(5))]
                twoproofs.append([j,h,M.bit_length(),I.bit_length(),N.bit_length(),DD.bit_length(),conservative,refined])
                updates.append([j,nold,k,b,resident,peak,M.bit_length(),D.bit_length()]);resident+=b-nold;newlogical+=b-nold;pairs+=nf*np
                src=numeric(index)+len(mono)+numeric(supportplans)+numeric(writerplans)+numeric(proofrows)+numeric(updates)+numeric(twoproofs)+rt[i]+R+len(fs)+len(ps)+nold+k+b+4096
                sourcepeak=max(sourcepeak,src);assert src<=150000
            resident-=nf;newlogical-=nf;ordered=[(j,bounds[j])for j in cols if j in bounds];nupper=sum(x[1][1]for x in ordered);assert nupper==resident-sm-R
            binupper=len(MAGIC)+sum(x[1][0]for x in ordered);idxupper=512;bc=[]
            for j,(size,n,qlen,elen)in ordered:
                entry=128+len(str(binupper))+len(str(size))+len(str(n))+qlen+elen;idxupper+=entry+1;bc.append([j,size,n,entry])
            writerpeak=sm+R+nupper+4*max((x[2]for x in bc),default=0)+256;alg=max(alg,writerpeak)
            factorindex=next(k for k in rc[str(i)]if cells[k]['column']==pc)
            supportplans.append({'row':i,'source_terms':rt[i],'new_row_terms_upper':nupper,'new_immutable_row_file_bytes_upper':binupper+idxupper,'factor_reference_cell_index':factorindex,'updates':updates})
            writerplans.append({'row':i,'binary_file_upper':binupper,'paired_index_upper':idxupper,'writer_owned_upper':writerpeak,'cells':bc});proofrows.append({'row':i,'known_two_power_content_bounds':twoproofs});totalbin+=binupper;totalidx+=idxupper
            del cur
        newlogical-=R;nu=sum(len(p['updates'])for p in supportplans);nbc=sum(len(p['cells'])for p in writerplans)
        meta=5*len(cells)+len(cells)+len(rc)+len(rows)+len(cols)+5*len(targets)+8*nu+3*len(targets)+512+4*nbc+4*len(targets)+512;wholeowned=meta+alg
        gp=out/f'candidate-{ordinal:02d}-full-gate.json';wp=out/f'candidate-{ordinal:02d}-binary-writer-bounds.json'
        g={'ordinal':ordinal,'source_indexed_pivot':c,'pivot_row_terms':R,'target_rows':len(targets),'charged_candidate_scalar_metadata':sm,'all_row_forecasts':supportplans,
           'logical_after_support_union_upper':newlogical,'owned_algebraic_and_scratch_bound':alg,'index_plan_receipt_integer_slots_upper':meta,
           'whole_owned_scalar_and_scratch_upper':wholeowned,'normalized_row_RAW_bound':normupper,'integer_coefficient_pair_operations_upper':pairs,
           'source_only':True,'row_storage':'exact binary paired indices; known two-power canonical content bound'}
        gp.write_bytes(enc(g));wp.write_bytes(enc({'source_gate':str(gp.relative_to(ROOT)),'source_gate_sha256':sha(gp),'rows':writerplans,'coefficient_content_proofs':proofrows}))
        active={str(p.relative_to(ROOT)):{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in paths+[gp,wp]};base=sum(x['bytes']for x in active.values());aux=3145728;total=base+normupper+totalbin+totalidx+aux
        safe=wholeowned<=150000 and total<=32*1024*1024
        summary={'ordinal':ordinal,'candidate':[pr,pc],'binary_upper':totalbin,'index_upper':totalidx,'fresh_RAW_upper':total,'owned_upper':wholeowned,'source_owned_upper':sourcepeak,'FULL_safe':safe};summaries.append(summary)
        if safe:
            proof={'all_passed':True,'source_only':True,'source_generation':a.generation,'candidate':[pr,pc],'unit_QQ':str(q),'unit_exponent':pe,
                   'source_text_gate':str(gp.relative_to(ROOT)),'source_text_gate_sha256':sha(gp),'binary_bounds_file':str(wp.relative_to(ROOT)),'binary_bounds_sha256':sha(wp),
                   'actual_READ_inputs':list(active.values()),'actual_READ_input_RAW_bytes':base,'new_binary_row_files_RAW_upper':totalbin,'new_paired_indices_RAW_upper':totalidx,
                   'normalized_pivot_RAW_upper':normupper,'other_control_code_executable_events_receipts_manifests_upper':aux,'fresh_RAW_upper':total,
                   'algebraic_and_writer_scratch_upper':alg,'integer_index_plan_writer_metadata_slots':meta,'whole_owned_upper':wholeowned,'FULL_safe':True,
                   'known_content_proof':'Each integer numerator is A*U+B*V, so divisible by2^min(v2(A),v2(B)). Canonical primitive integer coefficients<=floor(M/2^h); content numerator<=floor(M/2^min(h,v2(D))), denominator<=D/2^min(h,v2(D)). Fixed bit operations only, no polynomial GCD or source prime.',
                   'no_coefficient_product_expansion_or_pivot':True}
            (out/'whole-budget.json').write_bytes(enc(proof));selected=(ordinal,gp,wp);break
        del prow,supportplans,writerplans,proofrows
    (out/'completed-forecasts-summary.json').write_bytes(enc(summaries))
    report={'id':id,'all_passed':True,'source_only':True,'source_generation':a.generation,'FULL_safe_binary_row_pivot':selected is not None,
            'selected_gate_file':selected[1].name if selected else None,'selected_bounds_file':selected[2].name if selected else None,
            'completed_candidates':len(summaries),'start_ordinal':a.start,'shortlist_size':len(shortlist),'source_owned_peak':sourcepeak,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,
            'no_general_field_elimination_polynomial_GCD_point_or_prime':True,'not_all_units_unsafe_claim':True}
    (out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--generation',type=int,required=True);p.add_argument('--source',required=True);p.add_argument('--report-sha',required=True);p.add_argument('--start',type=int,default=0);p.add_argument('--count',type=int,default=16);main(p.parse_args())
