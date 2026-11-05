#!/usr/bin/env python3
"""Read-only source/checkpoint summary; no polynomial arithmetic."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT_SEQUENCE46_49_SUMMARY_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for z in iter(lambda:f.read(65536),b''):h.update(z)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(x):return(json.dumps(x,indent=2)+'\n').encode()
def main():
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);paths={Path(__file__),ROOT/f'research/notes/{ID}.md'}
    items=[(46,'COMPACT45_PIVOT','03d0a6449c9a7349554e84f47f8d4c19f28f6c8fa6a7737fb183128a5b0d9d50'),
           (47,'BINARY_ROW_PIVOT46','55735c6d3f59e8af76dceeed488c0ef068312b8e963b7fbed64aac2bbac7abbb'),
           (48,'BINARY_ROW_UNIT48','36b28357a603103998820a89e0fcebcf77c61d6553f3302eaf22398865ee45c3'),
           (49,'BINARY_ROW_UNIT49','9786d0303b545d05c48519b041ecc5c1a78d1d0569b18ac3c2913c0820f089e5')]
    records=[];logical=467009
    for k,stem,h in items:
        d=ROOT/f'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_{stem}_20261105_v64';rp=d/'report.json';assert sha(rp)==h
        r=read(rp);assert r['all_passed']and r['completed_total_on_reordered_path']==k and r['completed_exactly_one_pivot']
        np=d/'native-summary.json';lp=d/'final-RAW-ledger.json';mp=d/'generation-1.json';pp=d/'preflight.json'
        for p in(np,lp,mp,pp):assert sha(p)==r['evidence_sha256'][p.name]
        n,l,m,p=map(read,[np,lp,mp,pp]);assert p['prior_completed_steps']==k-1
        assert n['logical_input_terms']==logical and n['logical_after_terms']==m['logical_total_terms']==r['logical_total_coefficients']
        assert r['shape']==m['shape']==[236-k,83-k]and m['complete_pivot']
        assert r['whole_owned_including_index_scalars_upper']<=150000 and l['TOTAL_RAW_bytes']<=32*1024*1024 and r['total_worker_seconds']<60
        assert n['all_passed']and n['completed_pivot']and n[f'old{k-1}_steps_not_replayed']
        assert n['back_add_identities_checked']==r['all_back_adds_verified']
        logical=r['logical_total_coefficients'];paths.update([rp,np,lp,mp,pp])
        records.append({'generation':k,'pivot':r['new_pivot'],'shape':r['shape'],'logical_terms':logical,
                        'owned_with_metadata':r['whole_owned_including_index_scalars_upper'],'RAW':l['TOTAL_RAW_bytes'],
                        'seconds':r['total_worker_seconds'],'back_adds':r['all_back_adds_verified'],'changed_rows':r['all_changed_rows_hash_bound'],
                        'binary_output':r.get('new_binary_row_output',False),'report':str(rp.relative_to(ROOT)),'report_sha256':h,'manifest_sha256':sha(mp)})
    gates={
      'EXACT49_CONSTRAINT_PROFILE':'1c97a068769399f85198e5d47943c8f925bfaf7757accee5fb874259ab0c9a62',
      'EXACT49_FIXED_CONSTRAINT_MINOR_PREFLIGHT':'15ec0fa2870a5985257a5a9e62500141ade5619b2ba1afe7b50a9bddcafea06e'}
    for stem,h in gates.items():
        p=ROOT/f'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_{stem}_20261105_v64/report.json';assert sha(p)==h and read(p)['all_passed'];paths.add(p)
    snap=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM49_ACTIVE_FRAGMENT_MIXED_BINARY_COMPACTION_20261105_v64/report.json'
    assert sha(snap)=='62a95473e52e2a85054586bd63e28b79ae58850bebf49d6e5a0a4daeb69e7970';s=read(snap);assert s['all_passed']and s['exact_current_generation_compacted']==49 and s['logical_terms']==logical;paths.add(snap)
    evidence={'all_passed':True,'records':records,'current_exact_snapshot_report':str(snap.relative_to(ROOT)),'current_snapshot_report_sha256':sha(snap),
              'all_old_generations_unsafe_forecasts_and_back_adds_preserved':True,'no50_arithmetic_executed':True,
              'frontier':'Five complete selected unit-column forecasts RAW-unsafe, among27 remaining monomial entries. This is not an all-entry barrier.',
              'next_source_method':'Two actual38x7 constraint blocks,487+663 terms. Fixed nonzero5x5 minors at oldL rows6,14,16,41,42 follow from accepted row-basis/kernel algebra. Coefficient/adjugate worker NOT authorized or executed in this summary.'}
    ep=out/'four-complete-pivots-and-next-source-gate.json';ep.write_bytes(enc(evidence))
    report={'id':ID,'all_passed':True,'start_generation':45,'completed_generation':49,'new_completed_pivots':4,
            'final_shape':[187,34],'final_logical_terms':logical,'changed_rows_checked':sum(v['changed_rows']for v in records),
            'normalization_and_back_add_identities':sum(v['back_adds']for v in records),
            'max_worker_owned_including_metadata':max(v['owned_with_metadata']for v in records),'max_worker_RAW':max(v['RAW']for v in records),'max_worker_seconds':max(v['seconds']for v in records),
            'all_five_parameters_and_alpha_beta_constraints_retained':True,'logical_vs_owned_counts_explicit':True,
            'no_nonmonomial_arithmetic_GCD_point_prime_Gram_or_full_inverse':True,'no50_arithmetic':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{ep.name:sha(ep)}}
    rp=out/'report.json';rp.write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
if __name__=='__main__':main()
