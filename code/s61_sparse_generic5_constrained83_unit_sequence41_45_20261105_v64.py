#!/usr/bin/env python3
"""Read-only five-completed-step cursor/hash receipt; no coefficient work."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT_SEQUENCE41_45_20261105_v64'
SPECS=[
 (41,'COMPACT40','9f1b6b8f3f9a33c195beda99f71f9d511ccb3c5206760a9863394eb153f5efe2','49bf160f9a3d869f1506cb05b7538ad4653c2616d81709c2983b478fcf167dff'),
 (42,'COMPACT41','597450c47a58f78e916aaf94a628d4b6599d9446a7ab6114f44294d21af31686','6badb17ab11dab9718bf1370b0369555de6a02428f00a07b3f37f54ad256cd35'),
 (43,'COMPACT42','76b37cd1c5d08fe6ccb39b6ffb0e55502e1436382f5dcf856cbeada94eaaefb7','394721829a4a38c61eae760ec77465593cdf106458f049356067499f65ddeadb'),
 (44,'COMPACT43','0395c01d8e50d09b09dfb4b52a6116cec0e75428c0173061e896a3e90e924698','7ddb65d968d62a574544b39229cf6671490eee0d2857f6d7a87c38df1b57efab'),
 (45,'COMPACT44','4359317000435b6f50771e95154a5dca5898602177c3b520a336ef52598cbe1a','59fdd8a5cfd733a611c6fa58debc752db5f76d3bc48888a6f4f7f9b8cdf87ef2')]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,indent=2)+'\n').encode()
def main():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);items=[];inputs={};back=0;rows=0
    for generation,stem,rh,mh in SPECS:
        folder=ROOT/'reports'/f'S61_SPARSE_GENERIC5_CONSTRAINED83_{stem}_PIVOT_20261105_v64';rp=folder/'report.json';mp=folder/'generation-1.json'
        assert sha(rp)==rh and sha(mp)==mh;r=read(rp);m=read(mp)
        assert r['all_passed']and r['completed_total_on_reordered_path']==generation and m['complete_pivot']
        assert m['shape']==[236-generation,83-generation]==r['shape']
        n=read(folder/'native-summary.json');ledger=read(folder/'final-RAW-ledger.json');pre=read(folder/'preflight.json')
        assert n['all_passed']and n['completed_pivot']and n['old'+str(generation-1)+'_steps_not_replayed']
        assert n['whole_owned_scalar_and_scratch_upper']<150000 and ledger['TOTAL_RAW_bytes']<=32*1024*1024 and r['total_worker_seconds']<60
        assert n['back_add_identities_checked']==r['all_back_adds_verified']and n['updated_rows']==r['all_changed_rows_hash_bound']
        assert m['logical_total_terms']==r['logical_total_coefficients']
        for f in('native-summary.json','final-RAW-ledger.json','preflight.json','normalized-pivot-row.json','arithmetic-events.jsonl'):
            p=folder/f;assert sha(p)==r['evidence_sha256'][f];inputs[str(p.relative_to(ROOT))]=sha(p)
        inputs[str(rp.relative_to(ROOT))]=rh;inputs[str(mp.relative_to(ROOT))]=mh
        context=folder/'event-context.json'
        if context.exists():assert sha(context)==r['evidence_sha256'][context.name];inputs[str(context.relative_to(ROOT))]=sha(context)
        gp=next(ROOT/p for p in pre['input_sha256']if'/candidate-'in p and p.endswith('-full-gate.json'));g=read(gp);assert sha(gp)==pre['input_sha256'][str(gp.relative_to(ROOT))]
        back+=n['back_add_identities_checked'];rows+=n['updated_rows']
        items.append({'generation':generation,'pivot':n['pivot'],'QQ_unit_coefficient':g['source_indexed_pivot']['monomial_QQ'],
                      'unit_parameter_exponent_H_L_s_q_t':g['source_indexed_pivot']['monomial_exponent'],
                      'shape':m['shape'],'logical_terms':m['logical_total_terms'],'owned_including_metadata':n['whole_owned_scalar_and_scratch_upper'],
                      'changed_rows':n['updated_rows'],'exact_back_adds':n['back_add_identities_checked'],'fresh_RAW_bytes':ledger['TOTAL_RAW_bytes'],
                      'compile_and_worker_seconds':r['total_worker_seconds'],'native_binary_value_roundtrips':n.get('native_binary_value_roundtrips',0),
                      'report':str(rp.relative_to(ROOT)),'report_sha256':rh,'complete_manifest':str(mp.relative_to(ROOT)),'manifest_sha256':mh,
                      'source_gate':str(gp.relative_to(ROOT)),'source_gate_sha256':sha(gp),'normalized_row':str((folder/'normalized-pivot-row.json').relative_to(ROOT)),
                      'source_index':r['literal_input_index_path'],'event_schema':'compact tuples + shared source/factor context'if context.exists()else'full per-update source references'})
    assert back==14967 and rows==487
    cols=m['column_original_ids'];assert len(cols)==38
    cursor={'completed_generation':45,'complete_manifest':items[-1]['complete_manifest'],'manifest_sha256':items[-1]['manifest_sha256'],
            'shape':[191,38],'logical_terms':467009,'original_columns':cols,
            'remaining_column_labels':{'scalar_c':[j for j in cols if j==0],'alpha20':[j for j in cols if 1<=j<=20],
                                       'beta20':[j for j in cols if 21<=j<=40],'normal42':[j for j in cols if 41<=j<=82]},
            'original_row_ids':m['rows']and[x['original_row']for x in m['rows']],
            'next_reconstruction':'Use exact generation45 row references: new JSON row files plus the frozen44 binary snapshot/index; resolve by encoding before any new whole-worker gate.',
            'no_generation46_launched':True,'stop_reason':'all five additionally authorized pivots completed',
            'next_method':'A fresh mixed-source exact binary compaction gate, then a separately bounded monomial/constant candidate gate. If active-input RAW grows too large, source-gate exact binary row outputs with paired indices; no automatic cap increase or nonmonomial pivot.',
            'not_a_claim_that_remaining_units_are_unsafe':True}
    (out/'five-step-program-and-cursor.json').write_bytes(enc({'steps':items,'cursor':cursor}))
    inputs[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/f'research/notes/{ID}.md';inputs[str(note.relative_to(ROOT))]=sha(note)
    for p in(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_primitive_binary_native_20261105_v64.hpp'):
        inputs[str(p.relative_to(ROOT))]=sha(p)
    result={'id':ID,'all_passed':True,'start_generation':40,'completed_generation':45,'new_completed_pivots':5,
            'changed_rows_checked':rows,'normalization_and_back_add_identities':back,'final_shape':[191,38],'final_logical_terms':467009,
            'max_coefficient_worker_owned_including_metadata':max(x['owned_including_metadata']for x in items),
            'max_coefficient_worker_fresh_RAW':max(x['fresh_RAW_bytes']for x in items),'max_coefficient_worker_seconds':max(x['compile_and_worker_seconds']for x in items),
            'all_five_parameters_alpha_beta_constraints_retained':True,'all_old_generations_and_unsafe_forecasts_preserved':True,
            'binary_is_exact_uncompressed_integer_exponent_format':True,'no_field_gaussian_Gram_point_prime_GCD_or_nonmonomial_pivot':True,
            'no_complete_inverse_new_rank_or_fiber_exhaustion_claim':True,'generation46_not_started':True,
            'input_sha256':inputs,'evidence_sha256':{'five-step-program-and-cursor.json':sha(out/'five-step-program-and-cursor.json')},
            'seconds':time.monotonic()-started}
    (out/'report.json').write_bytes(enc(result));signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
