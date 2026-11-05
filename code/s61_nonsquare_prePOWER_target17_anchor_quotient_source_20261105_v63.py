#!/usr/bin/env python3
"""SOURCE ONLY exact anchor ratio/bare12 theorem and finite polynomial-ratio gate."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_TARGET17_ANCHOR_QUOTIENT_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 en='S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63';r=load(en);assert r['all_passed']and r['Cramer_is_zero']==[True,True,False]
 yy=load(en,'Cramer2.json');DY=int(yy['denominator']);Y={tuple(m):int(c)for m,c in yy['integer_coefficients']};assert len(Y)==1768
 kr='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63';kp=load(kr,'root237.json')['polynomial'];DK=math.lcm(*(Q(c).denominator for m,c in kp));K={tuple(m):int(Q(c)*DK)for m,c in kp}
 GK=math.gcd(*(abs(z)for z in K.values()));GY=math.gcd(*(abs(z)for z in Y.values()));Kp={m:z//GK for m,z in K.items()};Hp={m:z//GY for m,z in Y.items()}
 assert len(K)==1174 and all(Y.values()) and GK>0 and GY>0
 assert max(sum(m)for m in K)==30 and max(sum(m)for m in Y)==34
 for j in range(3):assert min(m[j]for m in Y)-min(m[j]for m in K)>=0
 bounds=[max(m[j]for m in Y)-max(m[j]for m in K)for j in range(3)];assert bounds==[2,2,1]
 mon=sorted([m for m in itertools.product(range(3),range(3),range(2))if sum(m)<=4],key=lambda m:(sum(m),m),reverse=True);assert len(mon)==17
 lm=max(Kp,key=lambda m:(sum(m),m));lc=Kp[lm];KL=sum(map(abs,Kp.values()));HL=sum(map(abs,Hp.values()));Dcap=abs(lc)**17;Ncap=HL*(abs(lc)+KL+1)**17
 verifyL=Dcap*HL+17*KL*Ncap;terms=math.comb(37,3);owned=2*(len(Kp)+len(Hp))+2*terms+2*17+15000+256
 proof=load('S61_NONSQUARE_PREPOWER_K92_PIVOT14_NONZERO_GATE_20261105_v63');assert proof['all_passed'];p14=load('S61_NONSQUARE_PREPOWER_K92_PIVOT14_SCHUR_20261105_v63','operator.json')
 k11=load('S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63','operator.json');assert k11['identities'][0]=='kappa3*detA8=cprod*Saux*T'
 src=load('S61_NONSQUARE_PREPOWER_K92_TROPICAL_SOURCE_PREFLIGHT_V2_20261105_v63','source-binding.json');assert src['physical_targets']==[22,11,13,25,15,19,14,18,28,30,29]
 raw=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];assert raw[0][0] and all(Q(c)for m,c in raw[0][0])
 # The primitive/content normalizations are only fixed QQ scalings of already
 # stored coefficients. The prospective17 ratio coefficients are NOT solved.
 scalar=-Q(GY*DK,DY*GK)
 raw_upper=2*1024*1024+17*((Ncap.bit_length()+Dcap.bit_length())*30103//100000+100)+terms*((verifyL.bit_length()*30103)//100000+80)
 gate=dict(known_lambda='-Y2/kappa3',primitive_ratio='lambda=QQscalar*(Hprimitive/Kprimitive)',QQscalar=str(scalar),K_common_D=str(DK),K_integer_content=str(GK),Y_common_D=str(DY),Y_integer_content=str(GY),
  candidate_monomials=mon,leading_K_monomial=lm,leading_K_integer=str(lc),monomial_order='graded lex (totaldegree,a,e,d)',
  exact_scalar_rule='rho_m=(H_(m+LMK)-sum(previous_n rho_n*K_(m+LMK-n)))/LC_K; future monomials cannot contribute by admissible order',
  numerator_upper=str(Ncap),common_denominator_upper=str(Dcap),integer_verification_L1_upper=str(verifyL),
  temporary_integer_bits_upper=2*max(Ncap.bit_length(),Dcap.bit_length(),verifyL.bit_length())+32,
  coefficient_products_upper=len(Kp)*17+len(Hp),all_support_positions_upper=terms,all_keys_coefficients_metadata_owned_upper=owned,
  proposed_RAW_upper=raw_upper,slot_gate_passes=owned<=150000,RAW_gate_passes=raw_upper<32*1024*1024,
  no_polynomial_ratio_coefficients_computed=True,no_polynomial_division_factor_GCD_point_prime_or_newtarget=True,
  success_requires_ALL_coefficients_identity=True,failure_keeps_original_rational_lambda=True)
 dump(out/'finite-lambda-ratio-gate.json',gate)
 source=dict(exact_lambda=dict(numerator=f'reports/{en}/Cramer2.json',denominator=f'reports/{kr}/root237.json',sign=-1),
  pair_source_current15=[2,4,5,8],true_raw_proportionality_target17_primary=True,independent_threeCramer_peer_pending=True,
  corrected_target17='oldcolumn17-lambda*oldcolumn22; allfourpairentriesZERO',
  corrected_lastrow='(v1-lambda*c_anchor)-corrected_row10*S10inverse*(v10-lambda*C_anchor)',
  bare12=dict(rows=['Q9','G18/a','G39/a','G68','G20','G115/c10','(G120-G99)/c11','G92/c8','G85/c6','aG129/c12','G99/c9','G134 current15row14'],
    columns_original32=src['physical_targets']+[17],column_normalization='last3 H3 columns28,30,29 dividedby L=aFC²; column17 and all9 pure columns scale1',
    base_M11_reference='reports/S61_NONSQUARE_PREPOWER_K92_PIVOT14_SCHUR_20261105_v63/operator.json#/exports/P_K92',
    full_source_input_reference='reports/S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63/operator.json#/exports/input15x22',
    actual_same_entry='r17=detM12/P_K92',no_scalar_determinant_or_lastrow_value_computed=True),
  descent='temporary anchor-component localization only in proof; injective localization of domain; finaltwoexpressionsregularonoldOmega+P; NOnewanchorchart',
  literal_anchor_component_NOT_claimed_unit=True,all_H3_and_pairtail2_and_ten_corrections_retained=True,no_full57_or_globalclaim=True)
 dump(out/'exact-anchor-and-bare12-source.json',source)
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 rep=dict(id=ID,all_passed=True,frozen=True,source_only=True,exact_lambda_source_formula=True,full_same_target17_normalized_Schur_quotient=True,
  no_new_anchor_unit=True,no_new_scalar_ratio_or_lastrow_value=True,independent_Cramer_peer_pending=True,
  candidate_polynomial_ratio_dimension=17,prospective_owned=owned,prospective_RAW=raw_upper,prospective_temp_bits=gate['temporary_integer_bits_upper'],
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
