#!/usr/bin/env python3
"""Exact coefficient/Nakayama consequence of the frozen joint105 jet; no jet job."""
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    t=time.monotonic();folder=ROOT/'reports'/STEM;folder.mkdir(exist_ok=False)
    jp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    jr=jp.parent/'report.json'
    symmetric=ROOT/'research/notes/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63.md'
    d=json.loads(jp.read_text());assert json.loads(jr.read_text())['all_passed']
    mon=d['monomials'];qcols=[mon.index([2,0]),mon.index([1,1]),mon.index([0,2])]
    jets={x['CB_physical_zero_row']:list(map(F,x['all_degree_le3_coefficients'])) for x in d['all8_CB_zero_restrictions']}
    assert all(v==0 for jet in jets.values() for a,v in zip(mon,jet) if sum(a)<2)
    M=[[jets[row][i] for i in qcols] for row in [18,60,67]]
    a,b,c=M[0];e,f,g=M[1];h,i,j=M[2]
    det=a*(f*j-g*i)-b*(e*j-g*h)+c*(e*i-f*h);assert det==F(625,768)
    gj=list(map(F,d['g_CA_all10_coefficients']));assert all(v==0 for m,v in zip(mon,gj) if sum(m)<3)
    ev={'three_actual_CB_restriction_rows':[18,60,67],'monomial_order':['r^2','r*s','s^2'],
        'exact_quadratic_coefficient_matrix':[[str(v) for v in r] for r in M],
        'determinant':str(det),'all8_restrictions_have_zero_constant_and_linear_terms':True,
        'actual_ideal_F_subset_m_squared':True,'F_plus_m_cubed_equals_m_squared':True,
        'Nakayama_module':'m^2/F, finite overQQ[[r,s]], equalsm*(m^2/F), hencezero',
        'exact_actual_CB_restriction_ideal':'(r,s)^2',
        'g_CA_order3_hence_redundant_mod_actual_CB_restrictions':True,
        'exact_same_physical_base_quotient':'R0/(I0_CA+I0_CB)=QQ[[r,s]]/(r,s)^2',
        'quotient_length':3,'quotient_Hilbert_function':[1,2],
        'T3_annihilated_by_z7_and_r_s_squared':True,'T3_is_finite_length':True,
        'T3_length_upper_bound':36,'T3_nonzero_or_length_or_socle_not_computed':True}
    (folder/'certificate.json').write_text(json.dumps(ev,indent=2)+'\n')
    rep={'id':STEM,'all_passed':True,'status':'exact_length3_double_annihilator_not_T3_nonzero',
         'actual_I0_CA_plus_I0_CB_quotient_length':3,'T3_finite_length_proved':True,
         'T3_length_upper_bound':36,'maximal_ideal_squared_annihilates_T3':True,
         'T3_nonzero_depth19_colon_or_actual_length_claimed':False,
         'no_new_jets_GB_points_source_regeneration_or_formal_coordinate_alignment_assumption':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),jp,jr,symmetric]},
         'evidence_sha256':{'certificate.json':sha(folder/'certificate.json')},'seconds':time.monotonic()-t}
    (folder/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':main()
