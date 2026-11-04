#!/usr/bin/env python3
"""One fixed native QQ factorization of accepted R434, and factor-list p tests only."""
import argparse,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_R434_FACTOR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert out.is_dir()and not(out/'report.json').exists();phase=read(out/'phase.json');pins={}
    def load(st,name):
        p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    sr='S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63';R=load(sr,'native.json')['remaining_polynomial'];assert len(R)==434 and max(sum(m)for m,c in R)==20
    assert load('S61_NONSQUARE_PREPOWER_H22_FIXED_MINORS_20261104_v63_REVIEW','report.json')['all_passed']
    p=load('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63','certificate.json')['p'];assert len(p)==21
    cs=[Q(c)for m,c in R];den=math.lcm(*(c.denominator for c in cs));zs=[c.numerator*(den//c.denominator)for c in cs];co=math.gcd(*zs)
    lc=zs[max(range(len(R)),key=lambda i:tuple(R[i][0]))];co=co if lc>0 else -co;scalar=Q(co,den);primitive=[[m,str(z//co)]for (m,c),z in zip(R,zs)]
    forecast=dict(input_R_terms=434,p_terms=21,total_degree=20,primitive_integer_coefficient_bits=max(abs(z//co).bit_length()for z in zs),signed_QQ_content=str(scalar),
      factor_output_support_upper=35420,accessible_support_upper=434+434+21+2*35420+3*1771+15000,
      source_factor_products_preserve_degree_at_most20=True,factor_internal_support_NOT_bounded=True,
      external_60s_2GB_guard=True,fullfactor_calls=1,no_Saux_factor=True,raw_artifact_forecast_bytes=30000000)
    assert forecast['accessible_support_upper']<150000;dump(out/'preflight.json',forecast)
    with(out/'input.txt').open('w')as f:
        for poly in [R,primitive,p]:
            f.write(str(len(poly))+'\n')
            for m,c in poly:f.write(' '.join([c]+list(map(str,m)))+'\n')
        f.write(str(scalar)+'\n')
    phase.update(phase='source and primitive QQ content pinned; compile before one call');dump(out/'phase.json',phase)
    cpp=ROOT/'code/s61_nonsquare_prePOWER_R434_factor_20261104_v63.cpp';exe=out/'native-worker'
    cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
    compile_run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,59-(time.monotonic()-start)));(out/'compile.log').write_text(compile_run.stdout);assert compile_run.returncode==0
    phase.update(phase='one native factor call about to run');dump(out/'phase.json',phase)
    run=subprocess.run([str(exe),str(out/'input.txt'),str(out/'native.json')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,59-(time.monotonic()-start)));(out/'native.log').write_text(run.stdout);assert run.returncode==0
    d=read(out/'native.json');assert d['all_passed']and d['single_factor_call']and d['all434_original_coefficients_reproduct_verified']
    assert d['max_accessible_terms_including_source_factors_scratch_reserve']<150000
    note=ROOT/f'research/notes/{ID}.md'
    for pth in [Path(__file__),cpp,note]:pins[str(pth.relative_to(ROOT))]=sha(pth)
    size=sum(p.stat().st_size for p in out.iterdir()if p!=exe);assert size<32*1024*1024
    phase.update(phase='one factor job complete; frozen result; no retry/next factor');dump(out/'phase.json',phase)
    result=dict(id=ID,all_passed=True,frozen=True,independent_review_pending=True,factor_count=len(d['factors']),
      factor_degrees=[x['total_degree']for x in d['factors']],factor_terms=[x['terms']for x in d['factors']],factor_multiplicities=[x['multiplicity']for x in d['factors']],
      gcd_with_named_p_is_constant=len(d['factor_list_gcd_with_p'])==1 and d['factor_list_gcd_with_p'][0][0]==[0,0,0],
      all_original_R_coefficients_reproduct_verified=True,no_full57_or_complement_removal_claim=True,
      max_accessible_terms=d['max_accessible_terms_including_source_factors_scratch_reserve'],factor_internal_support_uninstrumented=True,
      raw_artifact_bytes_before_report_excluding_native_binary=size,seconds=time.monotonic()-start,
      input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe})
    dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
