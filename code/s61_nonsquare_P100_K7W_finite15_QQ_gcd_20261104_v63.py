#!/usr/bin/env python3
"""ONE exact QQ gcd/xgcd of saved H15/R14; no source/minor reevaluation."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction
from pathlib import Path
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_K7W_FINITE15_QQ_GCD_20261104_v63'
def read(p):
    with p.open() as f:return json.load(f)
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out,run,resume=False):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57)
    if resume:
        assert out.exists()and (out/'compile-frontier.json').exists()and not(out/'native.json').exists()and not(out/'report.json').exists()
    else:out.mkdir(parents=True,exist_ok=False)
    b=ROOT/'reports/S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63'
    r=read(b/'report.json');assert sha(b/'report.json')=='6c6bf67338f63d991114572847886be2118003ce247739bc6beb0e89943ab7ae'
    assert sha(b/'native.json')==r['evidence_sha256']['native.json'];n=read(b/'native.json')
    assert r['complete']and n['all_passed']and n['literal_N_equals_HQ_plus_R_verified']and n['entire_special_fiber_quotient_dimension']==2
    keys=['monic_H15','finite_remainder_R'];polys=[n[k]for k in keys];del n
    profiles=[]
    for k,p in zip(keys,polys):
        cs=[Fraction(c)for m,c in p];den=math.lcm(*(c.denominator for c in cs));ns=[c.numerator*(den//c.denominator)for c in cs];cont=math.gcd(*ns)
        profiles.append({'key':k,'degree':max(m[0]for m,c in p),'terms':len(p),'max_raw_numerator_bits':max(abs(c.numerator).bit_length()for c in cs),'max_raw_denominator_bits':max(c.denominator.bit_length()for c in cs),'primitive_ZZ_height_bits':max(abs(x//cont).bit_length()for x in ns),'LCD_bits':den.bit_length()})
    assert [p['degree']for p in profiles]==[15,14]
    bitbound=14*(profiles[0]['primitive_ZZ_height_bits']+2)+15*(profiles[1]['primitive_ZZ_height_bits']+2)
    assert bitbound<=120754
    pre={'id':ID,'source_only_before_math':True,'source_keys':keys,'source_report_sha256':sha(b/'report.json'),'source_native_sha256':sha(b/'native.json'),'profiles':profiles,'primitive_Sylvester_size':29,'determinant_Hadamard_bit_bound':bitbound,'conservative_witness_serialized_bytes_bound':10*1024*1024,'conservative_live_coefficient_slots_bound':4000,'runtime_guards':{'seconds':60,'memory_bytes':2000000000,'total_output_bytes':32*1024*1024,'live_terms':150000,'coefficient_bits':300000},'one_QQ_xgcd_only':True,'no_factor_newprime_source_minor_or_NF':True}
    dump(out/'preflight.json',pre)
    ip=out/'native-input.txt';lines=[]
    for p in polys:lines.append(str(len(p)));lines.extend(f'{m[0]} {c}'for m,c in p)
    ip.write_text('\n'.join(lines)+'\n')
    phase=ROOT/'work/marked-square-v63/K7W-finite15-QQgcd-phase.json';phase.parent.mkdir(parents=True,exist_ok=True)
    dump(phase,{'status':'source_and_preflight_frozen_before_math','output':str(out.relative_to(ROOT)),'input_sha256':sha(ip),'preflight_sha256':sha(out/'preflight.json')})
    if not run:signal.alarm(0);print(json.dumps(pre,indent=2));return
    cpp=Path(__file__).with_suffix('.cpp');flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'QQ-gcd-native'
    subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,51-(time.monotonic()-start)))
    with(out/'native.log').open('w')as log:ret=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
    d=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':ret.returncode}
    passed=ret.returncode==0 and d['all_passed'];g=d.get('exact_QQ_gcd_degree')
    cert={'exact_chosen_fixed_minor_cover_generator_degree':g,'whole_P100_original_quartic_pointwise_lower150':passed and g==0,'proof':['Previous five actual9-minor normalized ideal=(H), and new actualK7W9-minor N=HQ+R.','Both saved contents times primitiveZZpolynomials equal H,R.','Exact literal g=S*primitiveH+T*primitiveR and H=g*QH,R=g*QR certify the QQ gcd, including multiplicity.','If g=1 the six fixed actual9-minors have no common complex zero on the whole permitted P100 curve, giving lower141+9=150.','Nonconstant g is only a fixed-minor frontier, not actual quartic rank failure.'],'no_new_source_minor_or_NF':True,'no_full239_or_global9_right_inverse_claim':True}
    dump(out/'certificate.json',cert)
    scalar={k:v for k,v in d.items()if not isinstance(v,(list,dict))}
    paths=[Path(__file__),cpp,b/'report.json',b/'native.json']
    report={'id':ID,'complete':passed,'exact_QQ_gcd_degree':g,'whole_P100_original_quartic_pointwise_lower150':passed and g==0,'native_audit':scalar,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);dump(out/'report.json',report)
    dump(phase,{'status':'one_QQ_gcd_complete','complete':passed,'output':str(out.relative_to(ROOT)),'report_sha256':sha(out/'report.json'),'exact_gcd_degree':g,'no_automatic_next_job':True})
    print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--run',action='store_true');p.add_argument('--resume-compile-frontier',action='store_true');a=p.parse_args();main(a.output.resolve(),a.run,a.resume_compile_frontier)
