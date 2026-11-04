#!/usr/bin/env python3
"""Single bounded independent Brown-PRS identity audit, preserving every completed prefix."""
import argparse,gzip,hashlib,json,math,os,resource,shlex,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_SUBRESULTANT_PRS_20261104_v63_REVIEW'
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def main(out,work):
    started=time.monotonic();deadline=started+120
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RESULTANT_20261104_v63';ar=read(author/'report.json');assert ar['all_passed']
    for rel,h in ar['input_sha256'].items():assert sha(ROOT/rel)==h,rel
    for rel,h in ar['evidence_sha256'].items():assert sha(author/rel)==h,rel
    oldgate=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63_REVIEW/report.json';assert read(oldgate)['all_passed']
    source=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';sd=read(source)
    pre=read(author/'content-preflight.json');assert all(x['deleted_product']=='1'for x in pre['inputs'])
    binary_input=work/'integer-input.txt'
    with gzip.open(author/'integer-input.txt.gz','rb')as f:binary_input.write_bytes(f.read())
    words=iter(binary_input.read_text().split());assert next(words)=='6'and next(words)=='0'
    for index,name in enumerate(['Z1','Z2']):
        count=int(next(words));got={}
        for _ in range(count):
            c=int(next(words));e=tuple(int(next(words))for _ in range(6));assert e not in got;got[e]=Q(c)
        c=sd['localized_constraints'][name];expected={};scalar=Q(pre['inputs'][index]['clearing_rational'])
        for j,key in enumerate(['coefficient1','coefficient2','denominator']):
            for a,b,r,s,R in c[key]['coefficients']:
                e=[s,r,R,0,0,0]
                if j==1:e[3]=1
                if j==2:e[4+index]=1
                ee=tuple(e);expected[ee]=expected.get(ee,Q(0))+scalar*Q(int(a),int(b))*(-1 if j==2 else 1)
        expected={e:c for e,c in expected.items()if c};assert got==expected
    assert next(words,None)is None
    (out/'input-binding.json').write_text(json.dumps({'all_passed':True,'all284_and1076_input_coefficients_rebuilt_from_accepted_constraints':True,'integer_scales':[200,400],'deleted_parameter_contents':[1,1],'variable_order':['s','r','R','Z0','Z1','Z2']},indent=2)+'\n')
    native=ROOT/'code/q_new_family_two_shape_prs_review_20261104_v63.cpp';binary=work/'independent-prs'
    error=None;ret=None
    try:
        flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint','openssl'],text=True))
        subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(.1,deadline-time.monotonic()-2))
        remaining=deadline-time.monotonic()-2;assert remaining>0
        def limits():
            resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));n=max(1,math.ceil(remaining));resource.setrlimit(resource.RLIMIT_CPU,(n,n))
        with(out/'native.log').open('w')as log:
            ret=subprocess.run([str(binary),str(binary_input),str(author/'resultant.txt'),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=remaining,preexec_fn=limits).returncode
        if ret:raise RuntimeError('independent PRS native exit '+str(ret))
    except Exception as e:error=repr(e)
    status=read(out/'native-report.json')if(out/'native-report.json').exists()else None
    frontier=read(out/'native-frontier.json')if(out/'native-frontier.json').exists()else None
    # Compress completed coefficient arrays. A stopped worker still keeps
    # immutable completed stages and the exact last frontier.
    for p in list(out.glob('*.txt')):
        with p.open('rb')as f,gzip.open(p.with_suffix('.txt.gz'),'wb')as g:
            for b in iter(lambda:f.read(1<<20),b''):g.write(b)
        p.unlink()
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    paths=[Path(__file__),native,ROOT/'research/notes'/(STEM+'.md'),author/'report.json',author/'integer-input.txt.gz',author/'resultant.txt',author/'content-preflight.json',author/'determinant-recipe.json',source,oldgate]
    report={'id':STEM,'all_passed':error is None and bool(status and status['all_passed']),'error':error,'native_exit_code':ret,
            'status':status,'last_frontier':frontier,'one_total_wall_seconds_limit':120,'address_space_bytes':2_000_000_000,'artifact_bytes_limit':32*1024*1024,
            'independent_Brown_PRS_not_a_second_resultant_call':True,'no_retry_or_alternate_algorithm':True,
            'raw_r_projection_only_not_NEW_equation_or_actual_target_index_one':True,
            's_linear_recovery_only_if_degree1_subresultant_leading_coefficient_is_unit':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-started,'frozen':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','error','last_frontier','seconds']},indent=2))
    if error:raise SystemExit(2)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    a=p.parse_args();main(a.output.resolve(),a.work.resolve())
