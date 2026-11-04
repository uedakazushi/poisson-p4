#!/usr/bin/env python3
"""ONE exact h22^2 division and existing named-unit stripping; no factor search."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);assert not(out/'report.json').exists();phase=read(out/'phase.json')
 base=ROOT/'reports/S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63';assert sha(base/'report.json')==phase['source_fixed3_report_SHA256']and sha(base/'native.json')==phase['fixed3_native_SHA256'];n=read(base/'native.json');assert n['all_passed']and not n['determinant_is_zero']
 hp=ROOT/'reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json';h=read(hp);h22=h['source_h9'][h['original_target_labels'].index(22)]
 factors=[[[[1,0,0],'1']],[[[0,0,1],'1']],[[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']],[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']],[[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']],[[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']]]
 polys=[n['exact_fixed3_determinant'],h22]+factors;reserve=sum(map(len,polys));forecast=dict(input_literal_terms=reserve,source_terms=1895,degree=40,max_source_QQ_bits=166,temporary_polynomial_support_upper=12000,conservative_total_live_terms_upper=130000,hard_live_terms=150000,hard_seconds=60,hard_output_bytes=33554432)
 (out/'preflight.json').write_text(json.dumps(forecast,indent=2)+'\n');text=[str(reserve)]
 for p in polys:text.append(str(len(p)));text.extend(' '.join(map(str,m))+' '+c for m,c in p)
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');cpp=Path(__file__).with_suffix('.cpp');exe=out/'named-strip-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-start)))
 d=read(out/'native.json')if(out/'native.json').exists()else dict(all_passed=False,returncode=r.returncode);passed=r.returncode==0 and d['all_passed']
 if passed:assert d['peak_live_terms_with_input_reserve']<=150000 and d['full_QQ_unit_product_reconstruction_verified']
 proof=dict(structural_identity='det(hP_i[j]-P_i[q]Q[j])=h^2 det([Q;P1;P2;P3] columns[q,j1,j2,j3])',
  proof='Over universal integer polynomial ring localize at h=Q[q]. SubtractP_i[q]/h timesfirstrow, givingdiagonalblock(h,M/h). Determinant=h*det(M)/h^3=det(M)/h^2. Bothsidespolynomial; universalidentitydescendsbeforelocalization.',
  h22_only_chartunit_not_globalFULLUunit=True,source_rowLCD_units_not_changed=True,
  if_constant='Only ifremainingisconstant,pure4unscaledminorisFULLUunit; thenno h22 restrictionisneededforthis4rowpureblock. Thisneverprovesfullkernel9.',
  if_nonconstant='Unknownremainingpolynomialretained; doNOTfactororremoveitszero-locus.',
  no_new_minor_general_factorsearch_point_prime=True)
 (out/'symbolic-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md',base/'report.json',base/'native.json',hp]
 report=dict(id=ID,all_passed=passed,frozen=True,independent_review_pending=True,exact_h22_squared_division=passed,
  named_unit_powers=d.get('named_unit_powers_a_d_F_B_C_J'),pure4_minor_terms=d.get('pure4_minor_terms'),remaining_terms=d.get('remaining_terms'),remaining_is_nonzero_QQ_constant=d.get('remaining_is_nonzero_QQ_constant'),
  peak_live_terms_with_input_reserve=d.get('peak_live_terms_with_input_reserve'),no_full_kernel9_or_globalFULLU57_claim=True,no_general_factorsearch_or_new_minor=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
