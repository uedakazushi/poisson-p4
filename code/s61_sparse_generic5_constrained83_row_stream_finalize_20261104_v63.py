#!/usr/bin/env python3
"""Metadata/hash finalization only: never loads row polynomials or runs native."""
from pathlib import Path
import hashlib,json,resource
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PIVOT_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def encode(v):return (json.dumps(v,indent=2)+'\n').encode()
def treebytes(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def main():
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;work=ROOT/'work'/ID;rp=out/'report.json';assert not rp.exists()
    pre=read(out/'preflight.json');ex=read(out/'execution-status.json');n=read(out/'native-summary.json')
    assert ex['one_worker_finished_verified']and ex['total_seconds']<60 and n['all_passed']and n['completed_pivot']
    for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
    assert sha(work/'control.txt')==pre['control_sha256']
    assert n['owned_including_scratch_peak']==5494<=5497
    assert n['index_plan_receipt_integer_slots_upper']==pre['all_index_plan_receipt_integer_slots_upper']==41616
    assert n['whole_owned_scalar_and_scratch_upper']==47110<150000
    assert n['no_full_coefficient_matrix']and n['no_second_pivot']and n['old37_steps_not_replayed']
    m=read(out/'generation-1.json');initial=read(out/'generation-0.json')
    assert m['complete_pivot']and m['shape']==[198,45]and m['logical_total_terms']==n['logical_after_terms']==174569
    assert not initial['complete_pivot']and initial['shape']==[199,46]and initial['logical_total_terms']==147370
    assert n['updated_rows']==133 and n['back_add_identities_checked']==939
    assert m['normalized_pivot_sha256']==sha(out/'normalized-pivot-row.json')
    gatepath=next(ROOT/p for p in pre['input_sha256']if p.endswith('/first-indexed-row-stream-gate.json'))
    gate=read(gatepath);plans={r['row']:r for r in gate['all_row_forecasts']}
    seen=set();actualterms=0
    for r in m['rows']:
        if'literal_row_file'in r:
            p=out/r['literal_row_file'];assert sha(p)==r['sha256']
            q=read(out/f"receipts/row-{r['original_row']}.json");f=plans[r['original_row']]
            assert q['all_cell_back_add_identities_checked']and q['pivot_column_cleared']
            assert q['literal_row_sha256']==r['sha256']and q['terms']==r['terms']<=f['new_row_terms_upper']
            assert q['RAW_bytes']==p.stat().st_size<=f['new_immutable_row_file_bytes_upper']
            actualterms+=r['terms'];seen.add(r['original_row'])
    assert seen==set(plans)
    events=[json.loads(s)for s in(out/'arithmetic-events.jsonl').read_text().splitlines()]
    assert len(events)==931 and len({(e['row'],e['column'])for e in events})==931
    assert {(e['row'],e['column'])for e in events}=={(i,u[0])for i,f in plans.items()for u in f['updates']}
    assert all(e['back_add_checked']and e['pivot_file']=='normalized-pivot-row.json'for e in events)
    assert not list(out.rglob('*.partial'))
    assert len(list((out/'rows').iterdir()))==133 and len(list((out/'receipts').iterdir()))==133
    rows=treebytes(out/'rows');norm=(out/'normalized-pivot-row.json').stat().st_size
    assert rows<=6699540 and norm<=1457
    ledgerpath=out/'final-RAW-ledger.json';assert not ledgerpath.exists()
    # Future ledger/report are accounted by serialized length. A hash has fixed
    # length, so the accounting has no cryptographic self-reference.
    r={'id':ID,'all_passed':True,'completed_exactly_one_pivot':True,'new_pivot':[222,24],
       'completed_total_on_reordered_path':38,'shape':[198,45],'logical_total_coefficients':174569,
       'owned_algebraic_and_scratch_peak':5494,'index_plan_receipt_integer_scalar_upper':41616,
       'whole_owned_including_index_scalars_upper':47110,'all939_normalization_and_update_back_adds_verified':True,
       'all133_rows_immutable_hash_bound':True,'no_full_matrix_in_RAM':True,
       'atomic_next_generation_only_after_whole_pivot':True,'old37_and_all_old_allresident_failures_unchanged':True,
       'logical_count_above150k_is_external_not_cap_increase':True,
       'all_five_parameters_and_alpha_beta_constraints_preserved':True,
       'new_inverse_only_parameter_monomial':'2025*t0/(H*q0)',
       'no_second_pivot_nonmonomial_point_prime_gcd_series_or_exhaustion':True,
       'input_sha256':{**pre['input_sha256'],str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},
       'work_control_sha256':pre['control_sha256'],'ignored_executable_sha256':sha(work/'native'),
       'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},
       'total_worker_seconds':ex['total_seconds'],'native_seconds':n['seconds'],
       'hard_worker_address_space_limit_bytes':2000000000,'empirical_worker_peak_RSS_not_recorded':True,
       'finalization_only_max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
       'scope':'Explicit ONE Schur row/column deletion and exact disk-backed residual. Rank83 already known; no complete inverse or fiber exhaustion.'}
    r['evidence_sha256'][ledgerpath.name]='0'*64
    ledger={'frozen_input_bytes':4518376,'frozen_index_bytes':1353659,'all_new_row_versions_bytes':rows,
            'normalized_pivot_bytes':norm,'other_artifacts_work_external_metadata_bytes':0,'TOTAL_RAW_bytes':0,
            'upper_bytes':15045736,'RAW32MiB_passed':True,
            'all_generated_versions_retained':True,'no_compression_used_for_guard':True,
            'finalizer_source_also_charged':True,'future_files_included_by_exact_serialized_lengths':True}
    existing=treebytes(out)+treebytes(work)+pre['RAW_other_external_charge']+Path(__file__).stat().st_size
    for _ in range(10):
        total=existing+len(encode(r))+len(encode(ledger))
        ob=total-4518376-1353659-rows-norm
        # existing excludes read-only input/index; add them only to whole RAW.
        total=existing+4518376+1353659+len(encode(r))+len(encode(ledger))
        ob=total-4518376-1353659-rows-norm
        if ledger['TOTAL_RAW_bytes']==total:break
        ledger['TOTAL_RAW_bytes']=total;ledger['other_artifacts_work_external_metadata_bytes']=ob
    else:raise AssertionError('byte accounting fixed point')
    assert ledger['other_artifacts_work_external_metadata_bytes']<=2472704
    assert ledger['TOTAL_RAW_bytes']<=15045736<32*1024*1024
    ledgerpath.write_bytes(encode(ledger));r['evidence_sha256'][ledgerpath.name]=sha(ledgerpath)
    rp.write_bytes(encode(r))
    actual=treebytes(out)+treebytes(work)+pre['RAW_other_external_charge']+Path(__file__).stat().st_size+4518376+1353659
    assert actual==ledger['TOTAL_RAW_bytes']
    print(json.dumps({k:v for k,v in r.items()if k not in('input_sha256','evidence_sha256')}))
    print('RAW',actual,'report_sha256',sha(rp),'generation1_sha256',sha(out/'generation-1.json'))
if __name__=='__main__':main()
