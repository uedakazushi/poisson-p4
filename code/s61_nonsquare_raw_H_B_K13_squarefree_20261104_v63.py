#!/usr/bin/env python3
"""ONE fixed derivative-unit certificate; no factorization or field assumption."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from pathlib import Path
from s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63 import pa,pm
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_SQUAREFREE_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(i):int(v)for i,v in p}
def ser(p):return [[i,str(v)]for i,v in sorted(p.items(),reverse=True)]
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    db=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63';ab=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63';ap=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63_REVIEW'
    assert sha(ab/'report.json')=='1c8c9d9f4640d18f51f85e5e5accf6d31f4dd9c28a0f86292268c5c102970c66'
    assert read(ap/'report.json')['all_passed']
    g=ip(read(db/'checkpoint.json')['g13']);assert max(g)==13
    derivative={i-1:i*v for i,v in g.items()if i};content=math.gcd(*(abs(v)for v in derivative.values()));p={i:v//content for i,v in derivative.items()}
    bits=max(abs(v).bit_length()for v in list(g.values())+list(p.values()));hadamard=26*(bits+5)+math.factorial(26).bit_length()+64
    bound=27*(200+(hadamard*30103+99999)//100000)
    (out/'preflight.json').write_text(json.dumps(dict(degrees=[13,12],input_coefficient_bits=bits,primitive_derivative_content=str(content),Bezout_RAW_bytes_upper=bound,
        hard_seconds=60,hard_AS_bytes=2000000000,hard_artifact_bytes=32*1024*1024,no_factorization_or_quotient_field_assumption=True),indent=2)+'\n')
    inp=out/'native-input.txt'
    with inp.open('w')as f:
        for z in[g,p]:
            f.write(str(len(z))+'\n')
            for i,v in sorted(z.items(),reverse=True):f.write(str(i)+' '+str(v)+'\n')
    cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_named_unit_20261104_v63.cpp';exe=work/'squarefree-unit'
    with(out/'compile.log').open('w')as log:
        subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    answer=out/'native-certificate.json'
    subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','derivative-unit',str(exe),str(inp),str(answer)],check=True,timeout=max(1,57-(time.monotonic()-start)))
    n=read(answer);S=ip(n.get('S',[]));T=ip(n.get('T',[]));squarefree=bool(n['unit'])
    cert=dict(g13=ser(g),derivative=ser(derivative),primitive_derivative=ser(p),derivative_content=str(content),native_gcd=n['gcd'],squarefree=squarefree)
    if squarefree:
        R=int(n['R']);assert pa(pm(S,p),pm(T,g))=={0:R}
        Tg={i:content*v for i,v in T.items()};Rfull=content*R;assert pa(pm(S,derivative),pm(Tg,g))=={0:Rfull}
        cert.update(derivative_cofactor=ser(S),g13_cofactor=ser(Tg),nonzero_QQ_integer_constant=str(Rfull),
            whole_integer_identity='S*g13_derivative+T*g13=R_nonzero',all_coefficients_checked=True,
            actual_Bbad_geometrically_reduced_13points_given_accepted_scheme_iso=True,
            no_QQ_irreducibility_or_rational_point_claim=True,no_Poisson_fullrank_or_nil_claim=True)
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths=[Path(__file__),cpp,ROOT/'code/s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63.py',ROOT/f'research/notes/{ID}.md',db/'report.json',db/'checkpoint.json',ab/'report.json',ap/'report.json']
    report=dict(id=ID,all_passed=True,gcd_g13_derivative_one=squarefree,actual_Bbad_geometrically_reduced_13points_given_accepted_scheme_iso=squarefree,
        no_factorization_newpoint_prime_rank_GB_or_Poisson_defect_inference=True,independent_review_pending=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
