#!/usr/bin/env python3
"""Read-only exact49 constraint blocks and selected-unit frontier receipt."""
from pathlib import Path
import hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode,encode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_CONSTRAINT_PROFILE_20261105_v64'
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
    if isinstance(v,list):return sum(map(numeric,v))
    if isinstance(v,dict):return sum(map(numeric,v.values()))
    return 0
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    sd=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM49_ACTIVE_FRAGMENT_MIXED_BINARY_COMPACTION_20261105_v64'
    sr=sd/'report.json';assert sha(sr)=='62a95473e52e2a85054586bd63e28b79ae58850bebf49d6e5a0a4daeb69e7970'
    r=read(sr);assert r['all_passed']and r['exact_current_generation_compacted']==49
    ip=sd/r['current_index'];mp=sd/r['current_snapshot'];assert sha(ip)==r['current_index_sha256']and sha(mp)==r['current_snapshot_sha256']
    oldrank=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW/report.json'
    ok=read(oldrank);assert ok['all_passed']and ok['rank_L18_plane3_normal42']
    raw=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_RAW_BLOCK_20261104_v63/report.json';assert read(raw)['all_passed']
    ix=read(ip);cs=ix['all_literal_cells'];rows=ix['all_original_row_ids'];cols=ix['all_original_column_ids']
    paths={sr,ip,mp,oldrank,raw,Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py'}
    assert ix['logical_shape']==[187,34]and ix['logical_coefficient_terms']==809647
    bdir=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_REORDERED_BRANCH_20261104_v63'
    brep=read(bdir/'report.json');assert brep['all_passed']and brep['unit_pivots']==30;paths.add(bdir/'report.json')
    sdir=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SHORTLIST_CONTINUATION_20261104_v63'
    srep=read(sdir/'report.json');assert srep['all_passed']and srep['total_on_reordered_path']==36;paths.add(sdir/'report.json')
    pairs=[];max_receipt_slots=0
    for pd,rr in [(bdir,brep),(sdir,srep)]:
        for p in sorted(pd.glob('pivot-*.json')):
            assert sha(p)==rr['evidence_sha256'][p.name];q=read(p);paths.add(p);max_receipt_slots=max(max_receipt_slots,numeric(q))
            if q['pivot_original_row']>=134:
                row,col=q['pivot_original_row'],q['pivot_original_column'];lo,hi=(1,21)if row<185 else(21,41)
                assert lo<=col<hi and q['completed']and q['all_row_update_back_add_identities_checked']
                assert all(lo<=j<hi for j,p in q['normalized_pivot_row'])
                pairs.append([row,col,str(p.relative_to(ROOT)),sha(p)])
            del q
    pd=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PIVOT_20261104_v63'
    pr=pd/'report.json';q=read(pr);assert q['all_passed']and q['completed_total_on_reordered_path']==38 and q['new_pivot']==[222,24]
    p=pd/'normalized-pivot-row.json';assert sha(p)==q['evidence_sha256'][p.name];z=read(p);assert all(21<=j<41 for j,p in z['entries']);paths.update([pr,p]);pairs.append([222,24,str(p.relative_to(ROOT)),sha(p)])
    assert len(pairs)==26 and len({q[0]for q in pairs})==26 and len({q[1]for q in pairs})==26
    assert {q[0]for q in pairs}==set(range(134,236))-set(rows)
    assert {q[1]for q in pairs}==set(range(1,41))-set(cols)
    groups=[];saved_slots=0;maxcell=0
    with mp.open('rb')as f:
        for name,rl,rh,cl,ch in [('alpha',134,185,1,21),('beta',185,236,21,41)]:
            rr=[x for x in rows if rl<=x<rh];cc=[x for x in cols if cl<=x<ch];entries=[]
            for k,c in enumerate(cs):
                if not rl<=c['row']<rh:continue
                assert cl<=c['column']<ch and 'monomial_QQ'not in c
                f.seek(c['payload_offset']);b=f.read(c['payload_bytes']);assert hashlib.sha256(b).hexdigest()==c['payload_sha256']
                p=decode(b);assert encode(p)==b and len(p['terms'])==c['terms'];entries.append([c['row'],c['column'],k,c['payload_sha256'],p]);maxcell=max(maxcell,c['terms'])
            assert len(rr)==38 and len(cc)==7 and sum(rl<=q[0]<rh for q in pairs)==13
            groups.append({'name':name,'rows':rr,'columns':cc,'shape':[38,7],'entries':entries,'terms':sum(len(e[4]['terms'])for e in entries),'rank5_from_accepted_L18_and13_pure_block_unit_deletions':True,'monomial_entries':0})
            saved_slots+=numeric(groups[-1])
    assert [g['terms']for g in groups]==[487,663]
    own=numeric(ix)+saved_slots+max_receipt_slots+24*maxcell+8192;assert own<=150000
    frontier=[]
    for k in range(5):
        fd=ROOT/f'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_BINARY_UNIT_SELECTOR49_FRONTIER{k}_20261105_v64';fr=fd/'report.json';fp=fd/'completed-forecasts-summary.json';q=read(fr)
        assert q['all_passed']and q['source_only']and q['source_generation']==49 and not q['FULL_safe_binary_row_pivot']
        assert sha(fp)==q['evidence_sha256'][fp.name];a=read(fp);assert len(a)==1 and a[0]['ordinal']==k and not a[0]['FULL_safe']
        assert a[0]['owned_upper']<=150000 and a[0]['fresh_RAW_upper']>32*1024*1024
        frontier.append({**a[0],'report':str(fr.relative_to(ROOT)),'report_sha256':sha(fr)});paths.update([fr,fp])
    mono=[{'row':c['row'],'column':c['column'],'QQ':c['monomial_QQ'],'exponent':c['monomial_exponent']}for c in cs if'monomial_QQ'in c];assert len(mono)==27
    evidence={'all_passed':True,'source_only':True,'source_generation':49,'actual_constraint_blocks':groups,
              'all26_pure_constraint_unit_deletion_sources':pairs,'all27_current_unit_entries':mono,'five_complete_selected_forecasts':frontier,
              'no_all27_unit_unsafe_claim':True,'remaining_unknown_partition':{'c':[0],'alpha':groups[0]['columns'],'beta':groups[1]['columns'],'normal':[c for c in cols if c>=41]},
              'next_source_method':'The two literal38x7 blocks have only1150 total terms and rank5 each. A future separately authorized source-bound nonmonomial5x5 instruction or reuse of the accepted L-field operator could remove ten constraint variables without expanding the two-plane kernel. No such elimination, determinant, polynomial GCD or coefficient product is executed here.',
              'rank83_already_known_not_recomputed':True,'full_inverse_and_exhaustion_uncomputed':True}
    ep=out/'actual-constraint-blocks-and-frontier.json';ep.write_bytes(enc(evidence))
    base=sum(p.stat().st_size for p in paths);upper=base+ep.stat().st_size+1024*1024;assert upper<32*1024*1024
    proof={'all_passed':True,'source_only':True,'owned_upper':own,'RAW_actual_READ_plus_output_reserve_upper':upper,'max_single_receipt_numeric_slots':max_receipt_slots,'no_arithmetic_worker':True}
    bp=out/'whole-budget.json';bp.write_bytes(enc(proof))
    report={'id':ID,'all_passed':True,'source_only':True,'source_generation':49,'shape':[187,34],'logical_terms':809647,
            'constraint_shapes':[[38,7],[38,7]],'constraint_terms':[487,663],'constraint_ranks_from_old_L18':[5,5],
            'selected_unit_frontier_count':5,'all_selected_frontier_RAW_unsafe':True,'not_all27_entries_claimed_unsafe':True,
            'owned_upper':own,'RAW_upper':upper,'no50_pivot_or_nonmonomial_arithmetic':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{ep.name:sha(ep),bp.name:sha(bp)},'seconds':time.monotonic()-start}
    rp=out/'report.json';rp.write_bytes(enc(report));assert base+sum(p.stat().st_size for p in out.iterdir())<=upper
    signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
if __name__=='__main__':main()
