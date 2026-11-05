#!/usr/bin/env python3
"""Known-input denominator/path bound for alternate native C72 arithmetic."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_ONE554_COMMON_DEN_HEIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def vals(v):
 if isinstance(v,(tuple,list)):
  for x in v:yield from vals(x)
 else:yield F(v)
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 ip='reports/AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63/ONE554-positive-four-native-and-nonzero644-outer-input.json';b=json.loads((ROOT/ip).read_text());p=b['full_primitive_input_J18_N_G_nu_GammaD']
 A=[q for row in p['Gamma_two56'] for i,q in row]+[q for ch in p['D_two_by10_by56'] for row in ch for q in row]
 classes={'known_partial_N9':p['N_known12_by12_by9_by4'],'canonical_G8':p['G8_C72'],'monics':p['nu_two_J18'],
  'Gamma_and_D_scalar_coefficients':A,'native_c0_constants':[v['constant_source_coefficient'] for v in b['parameterized_four_native_source_records']]}
 info={}
 for k,v in classes.items():
  xs=list(vals(v));d=1
  for x in xs:d=math.lcm(d,x.denominator);assert d.bit_length()<=384
  u=max([d]+[abs(x.numerator)*(d//x.denominator) for x in xs]);assert u.bit_length()<=384
  info[k]={'D':str(d),'U_including_scalar_one':str(u),'D_bits':d.bit_length(),'U_bits':u.bit_length(),'operands':len(xs)}
 # For ONE coefficient: four source terms*two Koszul rows*12 physical output
 # rows*16 N/G masks*9 N base cells*18 G base cells*19 optionalνbase choices.
 # At most11 Γ/D coefficient selectors are added. Inputmask is already fixed
 # by each source term. Alldoubleνpaths are structurally unsupported inJ18.
 path=4*2*12*16*9*18*19*11;assert path<2**28
 den=math.prod(int(v['D']) for v in info.values());u=math.prod(int(v['U_including_scalar_one']) for v in info.values())*(2**28)
 prior=max(den.bit_length(),u.bit_length());temp=3*prior+1;assert prior<=384 and temp<=1153
 ep='native-C72-known-input-common-denominator-and-direct-path-height-proof.json'
 e={'scalar_operand_classes':info,'common_D':str(den),'common_absolute_scaled_U_with_path2pow28':str(u),
  'explicit_native_direct_paths_upper':path,'rounded_paths2pow28':2**28,'committed_prior_bits':prior,'temporary_prior_bits':temp,
  'path_inputs_at_most_one_each':'constant cN, native N cell, G8 cell, oneνcell, oneΓ/Dcell',
  'double_nu_paths_structurally_absent':'J18 hasallνiνj supportsoutsidebasis; base exponents tested before coefficient products. Twoνhistories cannot contribute to a nonzero C72 table route.',
  'proper_intermediate_subsets_covered':'EachU includesD (scalar1), soSc/partialNSc/intermediatecoefficients withmissingclasses are dominated bythe sameD/U.',
  'requirements_for_value_algorithm':['Useonlytheknown9 partialNbasecells, keepingunknowns labelledformal; actualprojectedindependence followsfromweightedS0proof.','UsefixedC72 productroutes whichdiscardunsupportedbase exponents BEFOREcoefficientarithmetic.','Do notassertthepartialNSc matrixisactualnativeX inunrequestedslots; onlytheprojectedΓcubic+Dsecondisequal.','Verifyfinefinal152/111derivedcomponentenvelopesbeforewholeOLD/FIRST/outer/Qassembly.'],
  'joint_OLD_FIRST_outer_Q_sum_height_and_wholeowner_STILL_SEPARATE':True,'no_signed_native_or_outer_pairing_evaluated':True}
 (out/ep).write_text(json.dumps(e,separators=(',',':'))+'\n');rep={'id':STEM,'all_passed':True,'known_input_integer_bound_only':True,'native_alternate_C72_committed_prior':prior,'temporary_prior':temp,
  'input_sha256':{ip:sha(ROOT/ip)},'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'prior':prior,'temporary':temp,'paths':path}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
