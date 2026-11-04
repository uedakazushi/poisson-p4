#!/usr/bin/env python3
"""Fixed two-curve identity, complete unit profiles and typed raw-H cokernel."""
import argparse,ast,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as Q
from s61_nonsquare_Theta_rational_coordinates_20261104_v63 import sha,dec,enc,add,sm,mul,power,shift
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63'
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';dirs=[base/n for n in['S61_NONSQUARE_THETA_RATIONAL_COORDINATES_20261104_v63','S61_NONSQUARE_THETA_D14_UNIT_STRIP_20261104_v63','S61_NONSQUARE_THETA_REMAINING_DEFECTS_20261104_v63','S61_NONSQUARE_THETA_RAW12_DEFECT_20261104_v63','S61_NONSQUARE_RAW_H_THETA_ONE_DEFECT_20261104_v63','S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63']]
 rawdir=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';peer=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW/report.json'
 paths=[Path(__file__),Path(__file__).with_name('s61_nonsquare_Theta_rational_coordinates_20261104_v63.py')]+[p/f for p in dirs for f in['report.json','certificate.json']]+[rawdir/f for f in['report.json','preflight.json','checkpoint.json']]+[peer]
 rd=lambda p:json.loads(p.read_text());data=[]
 for p in dirs:
  r,c=rd(p/'report.json'),rd(p/'certificate.json');assert sha(p/'certificate.json')==r['certificate_sha256'];data.append(c)
 c,strip,remaining,raw12,old,theta=data;rawreport,pre,check=[rd(rawdir/f)for f in['report.json','preflight.json','checkpoint.json']];review=rd(peer)
 assert review['all_passed']and review['independent_frame_Cdelta_eta_gamma_and_full_source_clearings_passed']and review['all90_final_checkpoint_coefficients_and_fixed_row_operations_match']
 assert all(sha(rawdir/f)==rawreport['evidence_sha256'][f]for f in['preflight.json','checkpoint.json'])
 # Two original raw16x7 unit pivots: constant frame and literal unit functionals.
 H,Hi=pre['constant_coefficient_column_frame7x7'],pre['constant_inverse7x7']
 assert all(sum(Q(H[i][k])*Q(Hi[k][j])for k in range(7))==int(i==j)for i in range(7)for j in range(7))
 assert all(sum(Q(Hi[i][k])*Q(H[k][j])for k in range(7))==int(i==j)for i in range(7)for j in range(7))
 a={(1,0,0):Q(1)};e={(0,1,0):Q(1)};d={(0,0,1):Q(1)};one={(0,0,0):Q(1)}
 C=add(add(add(sm(a,100),sm(power(e,2,3),18)),sm(e,-9)),sm(one,-2));F=add(add(add(sm(a,500),sm(power(e,2,3),90)),sm(e,3)),sm(one,-24))
 B=add(mul(add(add(add(sm(a,50),sm(power(e,2,3),18)),sm(e,-9)),sm(one,-2)),d),add(sm(e,60),sm(one,-40)))
 delta=sm(power(C,2,3),Q(-64000));gamma=sm(mul(mul(F,B),power(C,2,3)),Q(-4,625));eta=sm(mul(mul(mul(a,F),power(B,2,3)),power(C,3,3)),Q(192,125))
 assert dec(pre['source_delta_change']['delta'])==delta
 pv=pre['pair_values'];assert dec(pv['gamma'])==gamma and dec(pv['eta'])==eta
 assert add(mul(dec(pv['cZ']),dec(pv['delta_z'])),mul(dec(pv['cW']),dec(pv['delta_w'])))==eta
 assert ast.literal_eval(pre['unit_Cdelta_row'])==(0,0,0,0,0,1,0)
 assert pre['remaining15_genuine_indices']==[i for i in range(16)if i!=6]
 assert pre['remaining6_coefficient_columns']==[0,1,2,3,4,6]
 assert check['completed_pivots']==1 and len(check['steps'])==1
 st=check['steps'][0];assert(st['pivot_row'],st['pivot_column'])==(1,1)
 assert st['pivot_unit']=={'constant':'-1/14062500','a_d_F_B_C_J_powers':[0,0,0,1,1,0]}
 specs=check['initial_row_divisors']+[st['pivot_unit'],st['row_scale_before']]+[q['row_divisor_after']for q in st['eliminations']]
 for q in specs:assert Q(q['constant'])and len(q['a_d_F_B_C_J_powers'])==6 and all(isinstance(v,int)for v in q['a_d_F_B_C_J_powers'])
 assert check['remaining_columns']==[0,2,3,4,5]and check['remaining_rows']==[0,2,3,4,5,6,7,8,9,10,11,12,13,14]
 # Fixed two-factor identity in the actual rational Theta ring.
 z={(1,0):Q(1)};D={(0,1):Q(1)};O={(0,0):Q(1)};zp20=add(z,sm(O,20));u=mul(z,zp20)
 P80=add(sm(u,3),sm(D,80));P100=add(sm(u,3),sm(D,100));nr=strip['native_exact_strip']['results'];Dstrip=next(x for x in nr if x['name']=='D14');Deltastrip=next(x for x in nr if x['name']=='Delta')
 assert dec(Dstrip['remaining_unknown'])==sm(mul(P80,P100),Q(-1));assert add(P100,P80,Q(-1))==sm(D,20)
 assert Deltastrip['remaining_is_nonzero_QQ_constant']and strip['Delta_global_unit_proved']
 assert raw12['defect_zero']and old['defect_identically_zero']and remaining['old_V_nonzero']
 assert next(q for q in remaining['two_fixed_defects']if q['raw_row']==9)['defect_zero']
 assert not next(q for q in remaining['two_fixed_defects']if q['raw_row']==14)['defect_zero']
 # Univariate graph profiles, all coefficients kept over QQ; no test point.
 Z={(1,):Q(1)};UO={(0,):Q(1)};Zp20=add(Z,sm(UO,20));u1=mul(Z,Zp20)
 def polynomial_on_graph(p,dm):
  ans={}
  for(i,j),v in p.items():
   assert j>=0;term=shift(power(dm,j,1),(i,),v);ans=add(ans,term)
  return ans
 def rational_on_graph(p,dm):
  clear=max([max(0,-m[1])for m in p]+[0]);num={}
  for(i,j),v in p.items():num=add(num,shift(power(dm,j+clear,1),(i,),v))
  return num,power(dm,clear,1)
 def lin(b,c0=0):return add(sm(Z,Q(b)),sm(UO,Q(c0)))
 profiles=[];target={k:dec(v)for k,v in c['target_coordinates'].items()}
 for m in[80,100]:
  dm=sm(u1,Q(-3,m));actual={name:polynomial_on_graph(target[name],dm)for name in['L','K','Fnum','Jnum']}
  if m==80:
   expected={'L':sm(mul(Z,lin(7,300)),Q(-1,8)),'K':sm(mul(Zp20,lin(7,80)),Q(-1,8)),
    'Fnum':sm(mul(mul(Z,power(Zp20,2,1)),lin(17,216)),Q(-1,16)),
    'Jnum':sm(mul(mul(Z,Zp20),lin(7,60)),Q(-1,8))}
   candidates={'a':(sm(lin(7,80),Q(-8)),sm(power(Z,2,1),Q(45))),
    'e':(sm(lin(1,40),Q(2)),sm(Z,Q(9))),'d':(dm,UO),'B':(mul(Z,lin(7,300)),sm(UO,Q(60))),
    'C':(sm(lin(7,300),Q(-4)),sm(Z,Q(9))),'F':(sm(lin(17,216),Q(-10)),sm(Z,Q(9))),
    'J':(lin(7,60),sm(UO,Q(30)))}
   open_factors=[Z,Zp20,lin(7,300),lin(7,80),lin(17,216),lin(7,60)]
  else:
   expected={'L':sm(mul(Z,lin(1,60)),Q(-1,2)),'K':sm(power(Zp20,2,1),Q(-1,2)),
    'Fnum':sm(mul(mul(Z,power(Zp20,2,1)),lin(203,2700)),Q(-1,250)),
    'Jnum':sm(mul(mul(Z,Zp20),lin(1,12)),Q(-1,2))}
   candidates={'a':(sm(Zp20,Q(-10)),sm(power(Z,2,1),Q(9))),
    'e':(lin(1,100),sm(Z,Q(9))),'d':(dm,UO),'B':(mul(Z,lin(1,60)),sm(UO,Q(12))),
    'C':(sm(lin(1,60),Q(-25)),sm(Z,Q(9))),'F':(sm(lin(203,2700),Q(-1)),sm(Z,Q(9))),
    'J':(lin(1,12),sm(UO,Q(6)))}
   open_factors=[Z,Zp20,lin(1,60),lin(203,2700),lin(1,12)]
  assert actual==expected and all(actual.values())
  # Original source coordinates, with rational full-localizer checks.
  source_rationals={name:rational_on_graph(target[name],dm)for name in['a','e','B','C']};source_rationals['d']=(dm,UO)
  source_rationals['F']=(actual['Fnum'],sm(power(dm,2,1),Q(40)));source_rationals['J']=(actual['Jnum'],sm(dm,Q(100)))
  for name,(n,q)in source_rationals.items():
   en,eq=candidates[name];assert n and q and mul(n,eq)==mul(en,q)
  op=UO
  for f in open_factors:assert f;op=mul(op,f)
  assert op
  profiles.append({'m':m,'d_graph':enc(dm),'L_K_Fnum_Jnum_exact_profiles':{name:enc(p)for name,p in actual.items()},
   'original_a_e_d_B_C_F_J_rationals':{name:{'numerator':enc(n),'denominator':enc(q)}for name,(n,q)in candidates.items()},
   'exact_univariate_open_factors':[enc(f)for f in open_factors],'nonzero_open_product':enc(op),
   'all_original_FULLU_localizers_nonzero_polynomials':True,'nonempty_rational_open_curve':True})
 operator={'orientation':'Original raw H coefficient row matrix presents cokernel R^7/(its16rows). All following are invertible row/target operations over R_Theta.',
  'first_two_raw_H_units':{'constant_GL7_and_inverse_checked':True,'source_C6_replacement':'delta*(C6-C1), delta=-64000C^2',
   'first_GL2_coefficient_column':'[[cZ/eta,cW/eta],[-delta_w/gamma,delta_z/gamma]], determinant1/gamma',
   'eta':'192*a*F*B^2*C^3/125','gamma':'-4*F*B*C^2/625','first_unit_row':'e5 from Cdelta',
   'known_unit_source_clearings':'Ci-(numerator_i/eta)*Cdelta',
   'second_fixed_unit':'checkpoint row1,col1, -B*C/14062500 after known-unit row clearings',
   'all_unit_specifications':specs,'accepted_full_coefficient_gate':str(peer.relative_to(ROOT)),
   'result':'coker(raw16x7) is coker(actualcheckpoint14x5) over FULLU and afterTheta basechange'},
  'Theta_unit_row8':{'known_inverse':'37500*C^2*g1/B^4','middle3zero':True,'last_target_scale':'B^4','zero_remaining_rows':[0,2,3,4,5,6],
   'result':'coker14x5 is coker7x4 [B|S]'},
  'B_unit3':{'selected_raw_rows':[7,13,10],'selected_middle_column_order':[0,2,1],'Delta_inverse':'inverse of explicit QQ*d^-7 product of certified unit gcd factors',
   'operator':'Use adj(Bpivot)/Delta, subtract the three pivot source rows and clear their S entries by a target shear; remaining scalar entries are the four defects/Delta.'},
  'last_scalar':{'zero_raw_rows':[9,11,12],'nonzero_raw_row':14,'scalar':'D14/Delta=unit*P80*P100',
   'unit_inverse_recipe':'reciprocal QQ content/d monomial; each gcd g has inverse (localizer/g)/localizer; Delta similarly. All localizers certifiedunits.'},
  'actual_cokernel':'coker(raw16x7) = coker(heavy14x5) = R_Theta/(P80*P100)',
  'CRT':'P100-P80=20d unit, so quotient is the product of its two reduced rational curve rings',
  'rank':'originalrawH rank6 exactly on either curve, rank7 off; heavy14x5 rank4 oncurves, rank5 off',
  'no_full_cubic57_or_quartic239_inference':True}
 (out/'operator-recipe.json').write_text(json.dumps(operator,separators=(',',':'))+'\n')
 cert={'P80':enc(P80),'P100':enc(P100),'literal_Q6_equals_minus_product':True,'P100_minus_P80_equals20d':True,
       'comaximal_certificate':'1=(P100-P80)/(20d)','curve_profiles':profiles,
       'two_initial_raw_H_unit_binding_verified':True,'raw_H_accepted_review_pinned':str(peer.relative_to(ROOT)),
       'exact_cokernel':'R_Theta/(P80*P100)','sharp_rank_on_and_off_curves':True,
       'uniform_raw_H7_is_false_on_two_nonempty_rational_curves':True,
       'no_point_GB_solver_new_minor_or_full_cubic_quartic_negative':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'fixed_P80_P100_identity_verified':True,'exact_comaximality_verified':True,
         'all8_translated_localizer_curve_profiles_verified':True,'all14_original_source_rational_profiles_verified':True,
         'two_initial_raw_H_unit_cokernel_binding_checked':True,'exact_Theta_heavy_and_rawH_cokernel':'R_Theta/(P80*P100)',
         'Theta_heavy_rank4_on_curves_5_off':True,'Theta_rawH_rank6_on_curves_7_off':True,
         'two_nonempty_rational_curves_no_newpoint':True,'no_full_cubic57_quartic239_negative':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'operator_recipe_sha256':sha(out/'operator-recipe.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
