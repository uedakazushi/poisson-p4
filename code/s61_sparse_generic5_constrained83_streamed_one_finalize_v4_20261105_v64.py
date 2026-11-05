#!/usr/bin/env python3
"""Metadata-only finalizer for one source-gated streamed pivot, no native run."""
from pathlib import Path
from fractions import Fraction
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,decode,encode
from fractions import Fraction
import argparse,hashlib,json,resource
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,indent=2)+'\n').encode()
def tree(p):return sum(q.stat().st_size for q in p.rglob('*')if q.is_file())
def main(id):
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/id;work=ROOT/'work'/id;rp=out/'report.json';assert not rp.exists()
    pre=read(out/'preflight.json');ex=read(out/'execution-status.json');n=read(out/'native-summary.json')
    assert pre['id']==id and ex['one_worker_finished_verified']and ex['total_seconds']<60 and n['all_passed']and n['completed_pivot']
    for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
    assert sha(work/'control.txt')==pre['control_sha256']
    assert n['owned_including_scratch_peak']<=pre['owned_algebraic_and_scratch_upper']
    assert n['index_plan_receipt_integer_slots_upper']==pre['index_plan_receipt_integer_slots_upper']
    assert n['whole_owned_scalar_and_scratch_upper']<=pre['whole_owned_upper']<150000
    assert n['no_full_coefficient_matrix']and n['no_second_pivot']and n[f"old{pre['prior_completed_steps']}_steps_not_replayed"]
    gp=next(ROOT/p for p in pre['input_sha256']if '/candidate-'in p and p.endswith('-full-gate.json'));g=read(gp)
    if pre.get('new_binary_row_output'):
        bp=ROOT/pre['binary_output_proof_path'];bb=read(bp);assert bb['all_passed']and bb['FULL_safe']
        assert bb['source_text_gate_sha256']==sha(gp)and bb['fresh_RAW_upper']==pre['RAW_TOTAL_bound']
        assert bb['new_binary_row_files_RAW_upper']+bb['new_paired_indices_RAW_upper']==pre['RAW_rows_cap']
        assert bb['algebraic_and_writer_scratch_upper']==pre['owned_algebraic_and_scratch_upper']
        assert bb['integer_index_plan_writer_metadata_slots']==pre['index_plan_receipt_integer_slots_upper']
        wb=read(ROOT/pre['binary_output_bounds_path']);binaryplans={p['row']:p for p in wb['rows']}
    else:
        assert g['FULL_safe']
        if g['RAW_TOTAL_bound']!=pre['RAW_TOTAL_bound']:
            bp=next(ROOT/p for p in pre['input_sha256']if p.endswith('/whole-budget.json'));b=read(bp)
            assert b['all_passed']and b['source_only']and b['whole_RAW_bound']==pre['RAW_TOTAL_bound']
            assert b['row_RAW_cap']==pre['RAW_rows_cap']and b['normalization_RAW_cap']==pre['RAW_normalization_cap']
            assert b['new_other_artifact_cap']==pre['RAW_other_cap']and b['whole_owned_upper']==pre['whole_owned_upper']
    m=read(out/'generation-1.json');m0=read(out/'generation-0.json')
    assert m['complete_pivot']and m['shape']==n['shape']and m['logical_total_terms']==n['logical_after_terms']
    assert not m0['complete_pivot']and m0['shape']==[n['shape'][0]+1,n['shape'][1]+1]and m0['logical_total_terms']==n['logical_input_terms']
    assert n['updated_rows']==g['target_rows']and n['back_add_identities_checked']==pre['expected_back_add_checks']
    assert m['normalized_pivot_sha256']==sha(out/'normalized-pivot-row.json')
    plans={r['row']:r for r in g['all_row_forecasts']};seen=set();writer_checks=0
    for r in m['rows']:
        if'literal_row_file'in r:
            p=out/r['literal_row_file'];q=read(out/f"receipts/row-{r['original_row']}.json");f=plans[r['original_row']]
            assert sha(p)==r['sha256']==q['literal_row_sha256']and q['all_cell_back_add_identities_checked']and q['pivot_column_cleared']
            assert q['terms']==r['terms']<=f['new_row_terms_upper']and q['RAW_bytes']==p.stat().st_size<=f['new_immutable_row_file_bytes_upper']
            seen.add(r['original_row'])
        elif 'literal_binary_row_file'in r:
            p=out/r['literal_binary_row_file'];ip=out/r['paired_index_file'];q=read(out/f"receipts/row-{r['original_row']}.json");f=plans[r['original_row']];bw=binaryplans[r['original_row']]
            assert sha(p)==r['sha256']==q['literal_row_sha256']and sha(ip)==r['paired_index_sha256']==q['paired_index_sha256']
            assert q['all_cell_back_add_identities_checked']and q['pivot_column_cleared']and q['binary_cell_value_roundtrips']
            ix=read(ip);assert ix['encoding']=='primitive-int-exponent-v1'and ix['source_binary_sha256']==r['sha256']and ix['original_row']==r['original_row']
            assert ix['source_binary_file']==p.name and ix['pivot_original_column']==n['pivot'][1]
            cb={x[0]:x for x in bw['cells']};totalterms=0;end=len(MAGIC);seen_cols=[]
            with p.open('rb')as stream:
                assert stream.read(len(MAGIC))==MAGIC
                for e in ix['entries']:
                    assert len(e)==7;col,offset,length,h,count,mono,exps=e
                    assert col in m['column_original_ids']and col in cb and offset==end and 0<count<=cb[col][2]and length<=cb[col][1]
                    stream.seek(offset);payload=stream.read(length);assert len(payload)==length and hashlib.sha256(payload).hexdigest()==h
                    value=decode(payload);assert encode(value)==payload and len(value['terms'])==count
                    if count==1:
                        assert str(Fraction(value['content'])*int(value['terms'][0][0]))==mono
                        assert [value['shift'][v]+value['terms'][0][v+1]for v in range(5)]==exps
                    else:assert mono is None and exps is None
                    end+=length;totalterms+=count;seen_cols.append(col);writer_checks+=1
                    del payload,value
            assert len(set(seen_cols))==len(seen_cols)and seen_cols==[j for j in m['column_original_ids']if j in seen_cols]
            assert end==p.stat().st_size<=bw['binary_file_upper']and ip.stat().st_size<=bw['paired_index_upper']
            assert totalterms==ix['total_terms']==q['terms']==r['terms']<=f['new_row_terms_upper']
            assert q['RAW_bytes']==p.stat().st_size+ip.stat().st_size
            seen.add(r['original_row'])
    assert seen==set(plans)
    if pre.get('new_binary_row_output'):assert writer_checks==n['writer_value_roundtrips']
    events=[json.loads(z)for z in(out/'arithmetic-events.jsonl').read_text().splitlines()]
    expected={(i,u[0])for i,p in plans.items()for u in p['updates']}
    if pre.get('compact_event_schema'):
        ctx=read(out/'event-context.json');assert ctx['pivot']==n['pivot']and ctx['normalized_pivot_sha256']==m['normalized_pivot_sha256']
        assert ctx['source_gate_sha256']==sha(gp)and dict(ctx['factor_source_cell_indices'])=={i:p['factor_reference_cell_index']for i,p in plans.items()}
        assert len(events)==len(expected)and {(e[0],e[1])for e in events}==expected and all(len(e)==4 and e[3]is True for e in events)
        bounds={(i,u[0]):u[3]for i,p in plans.items()for u in p['updates']};assert all(0<=e[2]<=bounds[e[0],e[1]]for e in events)
    else:
        assert len(events)==len(expected)and {(e['row'],e['column'])for e in events}==expected
        assert all(e['back_add_checked']for e in events)
    assert not list(out.rglob('*.partial'))
    rb=tree(out/'rows');nb=(out/'normalized-pivot-row.json').stat().st_size
    assert rb<=pre['RAW_rows_cap']and nb<=pre['RAW_normalization_cap']
    idx=next(p for p in pre['input_sha256']if p.endswith('/literal-file-index.json'))
    report={'id':id,'all_passed':True,'completed_exactly_one_pivot':True,'new_pivot':n['pivot'],
            'completed_total_on_reordered_path':pre['prior_completed_steps']+1,'shape':n['shape'],
            'logical_total_coefficients':n['logical_after_terms'],'owned_algebraic_and_scratch_peak':n['owned_including_scratch_peak'],
            'index_plan_receipt_integer_scalar_upper':n['index_plan_receipt_integer_slots_upper'],
            'whole_owned_including_index_scalars_upper':n['whole_owned_scalar_and_scratch_upper'],
            'all_back_adds_verified':n['back_add_identities_checked'],'all_changed_rows_hash_bound':len(seen),
            'all_five_parameters_and_constraints_preserved':True,'new_binary_row_output':pre.get('new_binary_row_output',False),'independent_stored_QQ_cell_roundtrips':writer_checks,'literal_input_index_path':idx,'output_generation_manifest':'generation-1.json',
            'no_full_matrix_in_RAM':True,'all_prior_generations_preserved':True,'no_second_pivot_nonmonomial_point_prime_GCD_or_exhaustion':True,
            'input_sha256':{**pre['input_sha256'],str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},
            'work_control_sha256':pre['control_sha256'],'ignored_executable_sha256':sha(work/'native'),
            'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},
            'total_worker_seconds':ex['total_seconds'],'native_seconds':n['seconds'],'hard_address_space_limit_bytes':2000000000,
            'empirical_worker_RSS_not_recorded':True,'scope':'One exact parameter-unit deletion; known rank is not a new result; no complete inverse or fiber exhaustion.'}
    lp=out/'final-RAW-ledger.json';assert not lp.exists();report['evidence_sha256'][lp.name]='0'*64
    ledger={'actual_READ_input_RAW_bytes':pre['RAW_actual_read_input_bytes'],'all_new_row_versions_bytes':rb,
            'normalized_pivot_bytes':nb,'other_artifacts_work_external_metadata_bytes':0,'TOTAL_RAW_bytes':0,
            'upper_bytes':pre['RAW_TOTAL_bound'],'all_generated_versions_retained':True,'no_gzip_or_artificial_deletion':True}
    finalizer_extra=0 if pre.get('finalizer_already_charged_in_external')else Path(__file__).stat().st_size
    existing=tree(out)+tree(work)+pre['RAW_other_external_charge']+finalizer_extra
    for _ in range(10):
        total=existing+pre['RAW_actual_read_input_bytes']+len(enc(ledger))+len(enc(report))
        if total==ledger['TOTAL_RAW_bytes']:break
        ledger['TOTAL_RAW_bytes']=total;ledger['other_artifacts_work_external_metadata_bytes']=total-pre['RAW_actual_read_input_bytes']-rb-nb
    else:raise AssertionError('RAW size fixed point')
    assert ledger['other_artifacts_work_external_metadata_bytes']<=pre['RAW_other_cap']
    assert total<=pre['RAW_TOTAL_bound']<=32*1024*1024
    lp.write_bytes(enc(ledger));report['evidence_sha256'][lp.name]=sha(lp);rp.write_bytes(enc(report))
    assert tree(out)+tree(work)+pre['RAW_other_external_charge']+finalizer_extra+pre['RAW_actual_read_input_bytes']==total
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
    print('RAW',total,'report_sha256',sha(rp),'manifest_sha256',sha(out/'generation-1.json'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--id',required=True);main(ap.parse_args().id)
