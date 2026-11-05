#!/usr/bin/env python3
"""Publish T0/C4 actual result and SOURCE bounds for a distinct fixed integer long-division peer."""
import argparse,hashlib,json,time
from pathlib import Path
from math import comb
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T0_C4_LONGDIVISION_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 t=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 actual=load('S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63')
 assert actual['all_passed']and actual['result']['quotient_terms']==4500 and actual['result']['full_C4_multiplyback_zero']
 custody=load('S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63');assert custody['all_passed']
 peer=load('S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63');assert peer['all_passed']
 source=load('S61_NONSQUARE_PREPOWER_T0_C4_NEXT_SOURCE_20261105_v63','fixed-C4-integer-quotient-SOURCE.json')
 inp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63/native-input.txt')
 value=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63/quotient-integer-common-D.txt')
 assert sha(inp)==custody['evidence_sha256'][inp.name]and sha(value)==actual['evidence_sha256'][value.name]
 # Metadata-only bounds: no quotient, product polynomial, or long division is evaluated here.
 S=sum(comb(h+3,3)*29**h*100**(16-h)for h in range(17))
 qbound=(((1<<326)-1)*S+100**20-1)//100**20
 rembound=((1<<326)-1)+qbound*129**4
 assert qbound.bit_length()<=303 and rembound.bit_length()<=333
 owned=18000+2*(14592+9690+25)+2*9690
 assert owned==85994 and owned<100000<150000
 proof=dict(SOURCE_ONLY_NO_LONGDIVISION_ARITHMETIC=True,
  independent_literal_input='T0C4Q0017208/amax20/shift0/qmax16; immutable D197 and body; NO producer import or formula engine',
  fixed_divisor='expand four copies of[(a,100),(e2,18),(e,-9),(1,-2)] by independent256 ordered multinomial enumeration; combine to25 terms; leadinga4 coefficient1004',
  algorithm='R=literalI0. While leadingkey has adegree>=4, require leadingcoef divisible1004, q=coef/1004, subtract q*fixedC4; no finite-series/binomial evaluation.',
  order='packedkey=a*1048576+e*1024+d is injective over the stated bounds and induces lex leading monomial',
  integrality='accepted Cvaluation4 and primitiveGauss ensure actualJ0 is integer; leadingterm elimination is unique polynomial division, so encountered q terms are exactly J0, without reading expected coefficients',
  quotient_support_bound9690=9690,remainder_support_bound14592=14592,
  why_partial_remainder_bounds='R=I0-C4*(already removed terms of trueJ0); both summands satisfy a<=20/e<=40/d<=18,total<=57',
  coefficient_height_SOURCE=dict(original_l1_bits326=326,quotient_l1_upper_bits=qbound.bit_length(),all_partial_remainder_l1_bits=rembound.bit_length(),stored_hard518=518),
  height_formula_bound_only='ceil((2^326-1)*sum_h binom(h+3,3)29^h100^(16-h)/10020); using bound not engine; ||R||1<=||I||1+1294||J||1',
  temporary_source_bound='q*fixedC4 coefficient <=qbound*1004; any before-cancel sum <=rembound; conservative operandbit+operandbit precheck1038 and stored518 suffice',
  whole_owned85994=owned,proposed_guard100000=100000,
  owned_equation='18000metadata+2*(R14592+Q9690+fixed25)+2*Q9690 sorting/serializationalias; input streamed directly intoR, no secondI/expectedQ map',
  no_duplicate_expected_polynomial=True,
  product_upper270000=270000,source_main_product_bound25_times9690=25*9690,
  scalar_exactdivision_reversechecks9690=9690,
  completeness='must reach R=0 and exactly4500 Q terms; stream ALL4500 expected originalD197 coefficients and compare; denominator identical; no expectedlookup affects division',
  result_required='ALL4500 exact coefficients plus fullzero remainder; fixedC4 integer longdivision, not finite-series engine',
  RAW_proposed8000000=8000000,soft45_hard60_RAM2GB=True,
  no_originala4_doubleoffset=True,no_GCD_otherfactor_tail_or_badlocus=True)
 dump('fixed-integer-longdivision-SOURCE.json',proof)
 binding=dict(primary_report_sha256=sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63/report.json'),
  actual_terms4500=4500,original_D197=True,full_primary_C4multiplyback_zero=True,
  independent_longdivision_NOT_run_by_this_program=True,
  accepted_Horner_Cvaluation4_sha256=sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63/report.json'),
  mathematical_scope='Q0=T0/(a4 C4); namedunits only. FULLU sufficient cover unchanged; no GCD/radicalcoprimality/strictcover/complement result.')
 dump('actual-result-publication-binding.json',binding)
 note=pin(ROOT/'research/notes/S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_RESULT_20261105_v63.md');pin(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,actual_result_publication_ONLY=True,
  independent_longdivision_SOURCE_only=True,independent_quotient_value_peer_NOT_run=True,
  actual_terms4500=4500,proposed_owned100000=100000,hard_stored518=518,proposed_RAW8000000=8000000,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-t)
 dump('report.json',r);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
