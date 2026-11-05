#!/usr/bin/env python3
"""Metadata-only finalizer for one source-gated streamed pivot, no native run."""
from pathlib import Path
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
    gp=next(ROOT/p for p in pre['input_sha256']if p.endswith('/candidate-00-full-gate.json'));g=read(gp)
    assert g['FULL_safe']and g['RAW_TOTAL_bound']==pre['RAW_TOTAL_bound']
    m=read(out/'generation-1.json');m0=read(out/'generation-0.json')
    assert m['complete_pivot']and m['shape']==n['shape']and m['logical_total_terms']==n['logical_after_terms']
    assert not m0['complete_pivot']and m0['shape']==[n['shape'][0]+1,n['shape'][1]+1]and m0['logical_total_terms']==n['logical_input_terms']
    assert n['updated_rows']==g['target_rows']and n['back_add_identities_checked']==pre['expected_back_add_checks']
    assert m['normalized_pivot_sha256']==sha(out/'normalized-pivot-row.json')
    plans={r['row']:r for r in g['all_row_forecasts']};seen=set()
    for r in m['rows']:
        if'literal_row_file'in r:
            p=out/r['literal_row_file'];q=read(out/f"receipts/row-{r['original_row']}.json");f=plans[r['original_row']]
            assert sha(p)==r['sha256']==q['literal_row_sha256']and q['all_cell_back_add_identities_checked']and q['pivot_column_cleared']
            assert q['terms']==r['terms']<=f['new_row_terms_upper']and q['RAW_bytes']==p.stat().st_size<=f['new_immutable_row_file_bytes_upper']
            seen.add(r['original_row'])
    assert seen==set(plans)
    events=[json.loads(z)for z in(out/'arithmetic-events.jsonl').read_text().splitlines()]
    expected={(i,u[0])for i,p in plans.items()for u in p['updates']}
    assert len(events)==len(expected)and {(e['row'],e['column'])for e in events}==expected
    assert all(e['back_add_checked']for e in events)and not list(out.rglob('*.partial'))
    rb=tree(out/'rows');nb=(out/'normalized-pivot-row.json').stat().st_size
    assert rb<=pre['RAW_rows_cap']and nb<=pre['RAW_normalization_cap']
    idx=next(p for p in pre['input_sha256']if p.endswith('/literal-file-index.json'))
    report={'id':id,'all_passed':True,'completed_exactly_one_pivot':True,'new_pivot':n['pivot'],
            'completed_total_on_reordered_path':pre['prior_completed_steps']+1,'shape':n['shape'],
            'logical_total_coefficients':n['logical_after_terms'],'owned_algebraic_and_scratch_peak':n['owned_including_scratch_peak'],
            'index_plan_receipt_integer_scalar_upper':n['index_plan_receipt_integer_slots_upper'],
            'whole_owned_including_index_scalars_upper':n['whole_owned_scalar_and_scratch_upper'],
            'all_back_adds_verified':n['back_add_identities_checked'],'all_changed_rows_hash_bound':len(seen),
            'all_five_parameters_and_constraints_preserved':True,'literal_input_index_path':idx,'output_generation_manifest':'generation-1.json',
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
    existing=tree(out)+tree(work)+pre['RAW_other_external_charge']+Path(__file__).stat().st_size
    for _ in range(10):
        total=existing+pre['RAW_actual_read_input_bytes']+len(enc(ledger))+len(enc(report))
        if total==ledger['TOTAL_RAW_bytes']:break
        ledger['TOTAL_RAW_bytes']=total;ledger['other_artifacts_work_external_metadata_bytes']=total-pre['RAW_actual_read_input_bytes']-rb-nb
    else:raise AssertionError('RAW size fixed point')
    assert ledger['other_artifacts_work_external_metadata_bytes']<=pre['RAW_other_cap']
    assert total<=pre['RAW_TOTAL_bound']<=32*1024*1024
    lp.write_bytes(enc(ledger));report['evidence_sha256'][lp.name]=sha(lp);rp.write_bytes(enc(report))
    assert tree(out)+tree(work)+pre['RAW_other_external_charge']+Path(__file__).stat().st_size+pre['RAW_actual_read_input_bytes']==total
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}))
    print('RAW',total,'report_sha256',sha(rp),'manifest_sha256',sha(out/'generation-1.json'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--id',required=True);main(ap.parse_args().id)
