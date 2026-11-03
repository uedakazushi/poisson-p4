#!/usr/bin/env python3
"""One60s/2GB exact existing-field rank probe; no automatic QQ promotion."""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import time
from s61_nonsquare_structural_reductions_review_20261103_v63 import read,sha

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_FINITE_SIGMA_QUARTIC_PROBE_20261103_v63'


def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports/S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63_REVIEW';review=read(parent/'report.json');assert review['all_passed']
    for rel,h in review['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in review['evidence_sha256'].items():assert sha(parent/rel)==h,rel
    ops,raw=parent/'own-operators.txt',parent/'own-raw-modp.txt';cpp=ROOT/'code/s61_nonsquare_finite_sigma_quartic_probe_20261103_v63.cpp';exe=out/'native-bin'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    subprocess.run(['c++','-O3','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,timeout=20,capture_output=True)
    block=out/'actual-block.u32le';native=out/'native.json';completed=False
    with(out/'native.log').open('w')as log:
        try:
            subprocess.run([str(exe),str(ops),str(raw),str(block),str(native)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,59-(time.monotonic()-start)))
            completed=True
        except subprocess.TimeoutExpired:
            pass
    block_hash=sha(block)if block.exists()else None;block_bytes=block.stat().st_size if block.exists()else None
    result=dict(id=STEM,native_completed=completed,one_rank_call_only=True,QQ_promotion_not_performed=True,second_rank_inverse_or_additional_field_not_launched=True,
        compact_inputs_regenerate_full_actual_block=True,actual_block_u32le_sha256=block_hash,actual_block_bytes=block_bytes,
        resource_cap_seconds=60,resource_cap_bytes=2000000000)
    if completed:
        n=read(native);result.update(n);result['exact_GF_diagnostic_only']=True
        result['rank_matches_expected']=n['rank']==2184
        result['negative_QQ_inference_from_low_GF_rank']=False
    paths=[Path(__file__),cpp,ROOT/f'research/notes/{STEM}.md',parent/'report.json',ops,raw]
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths}
    # The33MB block is deterministic generated data, retained in ignored work;
    # publish only its SHA plus the compact verified inputs and native receipt.
    result['seconds']=time.monotonic()-start;(out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    exe.unlink(missing_ok=True)
    print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
