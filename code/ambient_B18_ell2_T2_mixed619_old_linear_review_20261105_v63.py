#!/usr/bin/env python3
"""Read-only fixed702-factor formal-result check, no function values."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 bp='reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63';p=read(bp+'/report.json');assert p['all_passed'] and hs[bp+'/report.json'].startswith('83f5ab73')
 ep=bp+'/actual-formal619-oldsecond-weights-and-symbol-support.json';e=read(ep);assert hs[ep]==p['evidence_sha256'][Path(ep).name]
 gp='reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_GATE_20261105_v63';g=read(gp+'/report.json');d=read(gp+'/factored619-old-honest-second-linear-input-and-positive-bound.json');pr=read(gp+'_REVIEW/report.json');assert pr['all_passed'] and hs[gp+'_REVIEW/report.json']=='7d2918be6513ad05e830a449af7419301025af2132c51e1cf073cd96835629e0'
 for path,h in p['input_sha256'].items():assert sha(ROOT/path)==h
 # Group first, then sum each independent factor family. No callback value
 # participates in these rational coefficients.
 families={}
 for r in d['all619_four_channel_old_second_factor_records']:
  k=(r['case619'],r['channel4'],r['ordinary_pair8'],r['original_C0_symbol']);families.setdefault(k,[]).append(tuple(F(x) for x in r['rational_factors']))
 own={};maximum=0
 for k,terms in families.items():
  vals=[math.prod(v) for v in terms];q=sum(vals,F(0));maximum=max(maximum,*[max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in vals+[q]]);assert maximum<=46
  if q:own[k]=q
 saved={tuple(k):F(v) for k,v in e['exact_nonzero_old_second_symbol_coefficients']};assert own==saved and len(own)==576
 for key,ds,us in d['all_positive_symbol_coefficient_envelopes']:
  x=own.get(tuple(key),F(0));D,U=int(ds),int(us);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
 cases=sorted({k[0] for k in own});channels=sorted({(k[0],k[1]) for k in own});labels=sorted({k[3] for k in own});assert len(cases)==42 and len(channels)==110 and labels==e['nonzero_original_C0_symbols'] and cases==e['nonzero_case619_indices'] and [list(v) for v in channels]==e['nonzero_case_channel_pairs']
 assert all(0<=c<619 and 0<=ch<4 and 0<=pair<8 for c,ch,pair,label in own)
 # Absent factor families are zero FORMAL WEIGHTS, not relation values.
 for c in range(619):
  for ch in range(4):assert ((c,ch) in channels)==any(k[0]==c and k[1]==ch for k in own)
 ev=out/'independent702-factor-and-all-formal-slot-identity-checks.json';ev.write_text(json.dumps({'all702_factor_terms_regrouped_independently':True,'all576_coefficients_and_positive_envelopes_match':True,'nonzero_cases42':cases,'nonzero_case_channels110':[list(v) for v in channels],'all13_original_C0_symbols':labels,'all619_four_channel_coverage':True,'absent_formal_weight_not_zero_relation_value':True,'no_honestsecond_native_gauge_firstcross_or_carry_response_evaluated':True},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'formal_nonzero_slots':576,'formal_nonzero_cases':42,'formal_nonzero_case_channels':110,'no_new_function_or_source_values':True,'source_prior_bits':46,'actual_maximum_bits':maximum,'input_sha256':hs,'evidence_sha256':{ev.name:sha(ev)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','formal_nonzero_slots','formal_nonzero_cases','formal_nonzero_case_channels','actual_maximum_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
