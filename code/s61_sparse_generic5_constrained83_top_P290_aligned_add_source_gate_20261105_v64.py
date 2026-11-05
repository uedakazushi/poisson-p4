#!/usr/bin/env python3
"""SOURCE only: refine saved support bounds for one fixed-base align/add."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_TOP_P290_BATCH4_ALIGNED_SOURCE_GATE_20261105_v64'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_TOP_P290_BATCH4_STREAMCACHE_SOURCE_GATE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for z in iter(lambda:f.read(65536),b''):h.update(z)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def main():
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);od=ROOT/'reports'/OLD;orp=od/'report.json'
    assert sha(orp)=='46a501f5e4e81461c9033b92d46509aed913948ceb0a63ad27b623473770d16f';r=read(orp);assert r['all_passed']and r['source_only']and not r['FULL_safe_coefficient_job']
    op=od/'sixteen-source-support-height-cell-plans.json';obp=od/'whole-budget.json';assert sha(op)==r['evidence_sha256'][op.name]and sha(obp)==r['evidence_sha256'][obp.name]
    frame=read(op);old=read(obp);peak=old['future_preparation_owned_upper'];derivations=[]
    for p in frame['all16_cell_plans']:
        a=0;row=p['source_row_input_terms']
        for s in p['steps']:
            k=s['product_support_upper'];b=s['union_support_upper'];assert a<=b and k<=b
            cached=0 if s['kind']=='P_times_oldfree'else s['n_b'];base=290+row+cached+25000
            align=base+a+2*k+256;forward=base+a+k+b+256;back=forward
            # Includes ONE additional b coefficient-copy margin beyond all
            # explicit P objects. No in-place coefficient alias is assumed.
            aligned=base+a+k+2*b+256;parser=s['cached_factor_parser_upper'];bound=max(aligned,parser)
            derivations.append([p['row'],p['column'],s['step'],a,k,b,s['whole_owned_scratch_upper'],align,forward,back,aligned,parser,bound])
            s['previous_union_support_upper']=a;s['old_copy_adder_bound']=s['whole_owned_scratch_upper'];s['whole_owned_scratch_upper']=bound;s['fixed_common_base_alignment_bound']=aligned
            peak=max(peak,bound);a=b
        peak=max(peak,p['binary_roundtrip_phase_upper']);assert a==p['output_support_upper']
    assert peak<=56473 and peak==52712
    frame['addition_implementation']='Align product once to actual_accumulator_shift, verify actual QQ value hash, release unaligned product; direct same-base FLINT add/sub without two operand clones. Back-add target is the saved old accumulator, support<=a.'
    frame['old_SOURCE_guard_stop_unchanged']=str(orp.relative_to(ROOT));ep=out/op.name;ep.write_bytes(enc(frame))
    dp=out/'all96-fixed-base-owner-derivations.json';dp.write_bytes(enc({'all_passed':True,'source_only':True,'fields':['row','col','step','a','k','b','oldbound','align','forward','back','with_extra_b_margin','parser','newbound'],'derivations':derivations,'all_values_from_saved_source_supports':True,'support_or_coefficient_products_recomputed':False}))
    # Explicit coefficient READ ledger includes the original preserved support
    # files, but no duplicate support generation or old cell-plan copies.
    paths={ROOT/p for p in r['input_sha256']};paths.update([orp,op,obp,Path(__file__),ROOT/f'research/notes/{ID}.md'])
    for p in frame['all16_cell_plans']:
        sp=ROOT/p['support_file'];assert sha(sp)==p['support_sha256'];paths.add(sp)
    active=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)]
    raw=sum(p['bytes']for p in active)+ep.stat().st_size+dp.stat().st_size+old['all16_output_RAW_upper']+3*1024*1024
    # Charge future source budget/report bytes inside an explicit64KiB buffer.
    raw+=65536;assert raw<=32*1024*1024
    b={**old,'FULL_safe_sixteen_cell_coefficient_job':True,'actual_READ_inputs':active,'actual_READ_input_RAW_bytes':sum(p['bytes']for p in active),
       'all_source_plan_and_support_files_RAW_bytes':ep.stat().st_size+dp.stat().st_size+65536,'fresh_RAW_upper':raw,'whole_owned_with_fixed_metadata_upper':peak,
       'no_changed_parameter_point_prime_GCD_or_expression':True,'fixed_base_aligned_addition':True,'old_RAW_owner_guard_report_preserved':str(orp.relative_to(ROOT))}
    bp=out/'whole-budget.json';bp.write_bytes(enc(b))
    report={'id':ID,'all_passed':True,'source_only':True,'selected_rows':frame['selected_rows'],'cells':16,'FULL_safe_coefficient_job':True,
       'future_owned_upper':peak,'future_RAW_upper':raw,'coefficient_pair_upper':old['coefficient_term_pair_operations_upper'],'source_supports_recomputed':False,
       'input_sha256':{p['path']:p['sha256']for p in active},'evidence_sha256':{ep.name:sha(ep),dp.name:sha(dp),bp.name:sha(bp)},'seconds':time.monotonic()-begin}
    rp=out/'report.json';rp.write_bytes(enc(report));assert bp.stat().st_size+rp.stat().st_size<=65536
    signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
if __name__=='__main__':main()
