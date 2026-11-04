#!/usr/bin/env python3
"""One fixed symbolic F32003 divrem, no factor/rank/point search."""
import argparse,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_R434_SAUX_DIVISION_20261104_v63';P=32003
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
    def load(st,name):
        p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
    dump(out/'phase.json',dict(phase='persisted before one modular division',single_prime=P,single_division_authorized=True,no_point_or_factor_search=True))
    R=load('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json')['remaining_polynomial']
    S=load('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63','native.json')['exact_Saux_Laurent_polynomial'];assert len(R)==434 and len(S)==7022
    assert load('S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW','report.json')['all_passed']
    denominators=set()
    def reduce(poly):
        out=[]
        for m,c in poly:
            x=Q(c);assert x.denominator%P;denominators.add(x.denominator);y=x.numerator%P*pow(x.denominator,-1,P)%P
            if y:out.append([m,y])
        return out
    rp,sp=reduce(R),reduce(S);assert max(sum(m)for m,c in rp)==20 and max(sum(m)for m,c in sp)==61
    den_gate=dict(all_QQ_denominators_p_units=True,distinct_denominators=[str(x)for x in sorted(denominators)],
      denominator_inverse_modp=[[str(x),pow(x,-1,P)]for x in sorted(denominators)],
      R_has_coefficient_p_unit=True,R_degree_preserved=20,Saux_degree_preserved=61,
      full_reduced_R_terms=len(rp),full_reduced_Saux_terms=len(sp))
    dump(out/'coefficient-reduction-gates.json',den_gate)
    forecast=dict(input_modslots=len(rp)+len(sp),graded_quotient_degree_at_most=41,graded_quotient_terms_at_most=math.comb(44,3),
      remainder_standard_monomial_bound=math.comb(64,3)-math.comb(44,3),difference_terms_at_most=35442,product_terms_at_most=41664,
      live_upper_including_10000_reserve=149470,raw_artifact_bytes_upper=6000000,
      Python_parent_polynomial_arrays_released_before_native=True,internal_division_support_NOT_instrumented=True,
      prime=P,symbolic_polynomial_not_parameter_point=True,single_divrem_call=True)
    assert forecast['live_upper_including_10000_reserve']<150000;dump(out/'preflight.json',forecast)
    with(out/'input.txt').open('w')as f:
        for poly in [sp,rp]:
            f.write(str(len(poly))+'\n')
            for m,c in poly:f.write(' '.join(map(str,[c]+m))+'\n')
    del R,S,rp,sp
    cpp=ROOT/'code/s61_nonsquare_prePOWER_R434_Saux_division_20261104_v63.cpp';exe=out/'native-worker'
    compile_cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
    run=subprocess.run(['bash','-lc',compile_cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,59-(time.monotonic()-start)));(out/'compile.log').write_text(run.stdout);assert run.returncode==0
    dump(out/'phase.json',dict(phase='one divrem call about to execute; all coefficient gates passed',single_prime=P))
    cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(out/'input.txt')+' '+str(out/'native.json')
    run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,59-(time.monotonic()-start)));(out/'native.log').write_text(run.stdout);assert run.returncode==0
    n=read(out/'native.json');assert n['all_passed']and n['single_divrem_call']and n['whole_fixed_polynomial_identity_verified']
    assert n['max_accessible_coefficient_slots_with_scratch_reserve']<150000
    note=ROOT/f'research/notes/{ID}.md'
    for p in [Path(__file__),cpp,note]:pins[str(p.relative_to(ROOT))]=sha(p)
    size=sum(p.stat().st_size for p in out.iterdir()if p!=exe);assert size<32*1024*1024
    dump(out/'phase.json',dict(phase='fixed polynomial division complete; FINAL; no next scalar job',single_prime=P))
    result=dict(id=ID,all_passed=True,frozen=True,independent_review_pending=True,
      QQ_R_does_not_divide_Saux_proved=n['Rbar_does_not_divide_Sauxbar'],QQ_gcd1_ONLY_if_independent_R_irreducibility_gate=True,
      quotient_terms=n['quotient_terms'],remainder_terms=n['remainder_terms'],max_accessible_slots=n['max_accessible_coefficient_slots_with_scratch_reserve'],
      no_unit_ideal_or_empty_intersection_or_full57_claim=True,no_new_point_prime_or_factor=True,
      raw_artifact_bytes_before_report_excluding_binary=size,seconds=time.monotonic()-start,
      input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe})
    dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
