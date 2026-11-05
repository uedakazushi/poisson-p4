#!/usr/bin/env python3
"""Literal byte-custody repair only; no source arithmetic is rerun."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FOURROW_SOURCE_CUSTODY_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,indent=2)+'\n').encode()
def main():
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    od=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_TOP_P290_FOURROW_SOURCE_GATE_20261105_v64';orp=od/'report.json';oh=sha(orp)
    assert oh=='05817f4704c3b38930625249c8dea26205201f5f17bb83c8e106ab9b80b3dafc';old=read(orp);assert not old['all_passed']and old['honest_SOURCE_guard_stop']
    td=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_TOP_P290_FOURROW_THIN_SOURCE_GATE_20261105_v64';trp=td/'report.json';assert sha(trp)=='3b64458d6a4a6362f1234660e1aff73d9edfb176aa30dcfe31d589125789c83c';thin=read(trp);assert thin['all_passed']
    aliases=[];inputs={orp,trp,Path(__file__),ROOT/f'research/notes/{ID}.md'}
    for name,h in old['evidence_sha256'].items():
        advertised=od/name
        if advertised.exists():assert sha(advertised)==h;inputs.add(advertised);continue
        physical=od/'supports'/name;assert physical.is_file()and sha(physical)==h
        payload=physical.read_bytes();assert hashlib.sha256(payload).hexdigest()==h
        # New exact literal aliases, never modify any prior byte or report.
        with advertised.open('xb')as f:f.write(payload)
        assert advertised.read_bytes()==payload and sha(advertised)==h;inputs.add(physical)
        aliases.append({'declared_path':str(advertised.relative_to(ROOT)),'preserved_original_path':str(physical.relative_to(ROOT)),'sha256':h,'bytes':len(payload),'representation':'byte-identical literal alias','no_recomputation':True})
    assert len(aliases)==12 and sha(orp)==oh
    for name,h in old['evidence_sha256'].items():assert (od/name).is_file()and sha(od/name)==h
    all16=[]
    for row in[52,11,44,51]:
        for col in[7,13,27,33]:
            name=f'row-{row}-col-{col}.json';p=td/'supports'/name;key='supports/'+name;assert sha(p)==thin['evidence_sha256'][key];inputs.add(p)
            predecessor=od/'supports'/name
            if predecessor.exists():assert p.read_bytes()==predecessor.read_bytes();inputs.add(predecessor)
            all16.append({'row':row,'column':col,'thin_source_path':str(p.relative_to(ROOT)),'sha256':sha(p),
                          'old_failed_attempt_source_path':str(predecessor.relative_to(ROOT))if predecessor.exists()else None,
                          'scope':'first completed in successful thin SOURCE'if row==51 else'byte-identical to preserved old SOURCE support'})
    receipt={'all_passed':True,'no_arithmetic_rerun':True,'old_report_sha256_unchanged':oh,'old_failed_guard_status_unchanged':True,
        'all12_declared_paths_now_literal_and_hash_valid':aliases,'strict_original_all_evidence_path_checks_pass':True,'all16_successful_thin_source_supports':all16,
        'clarification':'The failed source had12 completed support arrays and stopped before row51. Its basename evidence keys omitted supports/. The successful thin source has all16. No original coefficient or SOURCE success claim is changed.'}
    ep=out/'literal-alias-and-all16-custody-receipt.json';ep.write_bytes(enc(receipt))
    report={'id':ID,'all_passed':True,'literal_aliases_created':12,'successful_thin_supports_verified':16,'old_report_unchanged':True,
        'no_coefficient_or_support_recomputation':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-begin}
    rp=out/'report.json';rp.write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
if __name__=='__main__':main()
