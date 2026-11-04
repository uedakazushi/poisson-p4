#!/usr/bin/env python3
"""Only raw9/raw14 defects, plus the necessary old V nonzero composition."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from fractions import Fraction as Q
from s61_nonsquare_Theta_raw12_defect_20261104_v63 import sha,dec,enc,add,mul,shift
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_REMAINING_DEFECTS_20261104_v63'
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';gdir=base/'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';cdir=base/'S61_NONSQUARE_THETA_RATIONAL_COORDINATES_20261104_v63';tdir=base/'S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63';ddir=base/'S61_NONSQUARE_THETA_DELTA_LAMBDA_PROFILES_20261104_v63';rdir=base/'S61_NONSQUARE_THETA_RAW12_DEFECT_20261104_v63';odir=base/'S61_NONSQUARE_RAW_H_THETA_ONE_DEFECT_20261104_v63'
 dirs=[gdir,cdir,tdir,ddir,rdir,odir];cpp=Path(__file__).with_suffix('.cpp');paths=[Path(__file__),cpp,Path(__file__).with_name('s61_nonsquare_Theta_raw12_defect_20261104_v63.py')]+[p/f for p in dirs for f in['report.json','certificate.json']]
 rd=lambda p:json.loads(p.read_text());data=[]
 for p in dirs:
  report,cert=rd(p/'report.json'),rd(p/'certificate.json');assert sha(p/'certificate.json')==report['certificate_sha256'];data.append(cert)
 g,c,t,d,r,old=data;B=[[dec(p)for p in row]for row in g['actual_block7x3']];Delta=dec(g['nonzero_selected_minor_Delta']);piv=g['selected_rows7'];assert piv==[0,5,2]
 for i in[1,6]:
  for j in range(3):
   rhs={}
   for k in range(3):rhs=add(rhs,mul(dec(g['cleared_upper_row_coordinates7xr'][i][k]),B[piv[k]][j]))
   assert mul(Delta,B[i][j])==rhs
 A=dec(c['target_coordinates']['a']);E=dec(c['target_coordinates']['e']);An=shift(A,(0,2),Q(500));En=shift(E,(0,1),Q(60));text=[]
 for q in[An,En]:text.append(str(len(q)));text.extend(f'{v} {m[0]} {m[1]}'for m,v in sorted(q.items()))
 items=[('S_raw9',t['fixed_7x4_matrix'][1][3]),('S_raw14',t['fixed_7x4_matrix'][6][3]),('S_raw13',t['fixed_7x4_matrix'][5][3]),('old_V',g['actual_block7x3'][0][0])]
 text.append('4')
 for name,p in items:
  text.append(f'{name} {len(p)}');text.extend(f'{v} {m[0]} {m[1]} {m[2]}'for m,v in p)
 nativein=out/'native-input.txt';nativein.write_text('\n'.join(text)+'\n');flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'substitute-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20);subprocess.run([str(exe),str(nativein),str(out/'native.json')],check=True,timeout=40)
 nv=rd(out/'native.json');new={x['name']:dec(x['polynomial'])for x in nv['images']};assert len(new)==4 and new['old_V']
 S={x['name']:dec(x['polynomial'])for x in r['actual_transformed_S']};S.update({name:p for name,p in new.items()if name.startswith('S_')});images={x['name']:dec(x['polynomial'])for x in d['actual_images']};delta=images['Delta'];assert delta
 checks=[];stored=sum(len(p)for p in S.values())+len(new['old_V'])
 for i,raw in[(1,9),(6,14)]:
  products=[mul(delta,S[f'S_raw{raw}'])]+[mul(images[f'lambda_{i}_{j}'],S[f'S_raw{t["full_source_rows"][piv[j]]}'])for j in range(3)]
  defect=products[0]
  for p in products[1:]:defect=add(defect,p,Q(-1))
  stored+=sum(len(p)for p in products)+len(defect);assert stored<=150000
  checks.append({'raw_row':raw,'all3_original_B_coordinate_identities_verified':True,'transformed_product_polynomials':[enc(p)for p in products],'defect':enc(defect),'defect_zero':not defect})
 allzero=all(x['defect_zero']for x in checks)
 assert r['defect_zero']and old['defect_identically_zero']
 # The old V-cleared zero implies its Delta-coordinate version by cancellation
 # in the ACTUAL rational Theta domain. V is nonzero, never declared a unit.
 cert={'source_rows_checked_before_chart':[9,14],'all6_full_original_B_coordinate_identities_verified':True,
       'new_actual_source_images':nv['images'],'reused_S_images_from_raw12_certificate':[x for x in r['actual_transformed_S']if x['name']in['S_raw7','S_raw10']],
       'two_fixed_defects':checks,'old_V_nonzero':True,'old_V_actual_image':enc(new['old_V']),
       'old_raw11_domain_cancellation':'V*(Delta*S11-lambda*S7)=Delta*(V*S11-U*S7)+(Delta*U-lambda*V)*S7=0. V!=0 and domain imply Delta-coordinate defect0. No V inverse globally.',
       'all_four_nonpivot_defects_zero':allzero,
       'sharp_consequence_if_allzero':{'Theta_ring':'actual localization of QQ[zeta,d], integral, nonempty',
        'Delta_nonzero':'pinned exact56term polynomial; no inversion as globalunit',
        'generic_augmented_7x4_rank':3,'all_point_augmented_7x4_rank_upper':3,
        'generic_heavy14x5_rank':4,'all_point_heavy14x5_rank_upper':4,
        'generic_original_raw16x7_rank':6,'all_point_original_raw16x7_rank_upper':6,
        'proof':'Over domainfractionfield nonzeroDelta gives B3 and all4 rowcoordinates extend to S. Thus augmentedrank3; every4minor is0 in domain, so rank<=3 everyresiduefield. Acceptedunitrow8 adds1; acceptedtwoFULLUrawHunitpivots add2. Nonzero translatedFULLUlocalizer product gives nonemptyTheta. Hence uniformrawH7 false.'}if allzero else None,
       'no_full_cubic57_or_quartic239_negative_inference':True,'no_other_minor_solver_GB_point_prime_or_defect_evaluated':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'all6_original_source_B_coordinate_identities_verified':True,
         'new_source_image_terms':{name:len(p)for name,p in new.items()},'defect_termcounts':{str(x['raw_row']):len(x['defect'])for x in checks},
         'both_remaining_defects_zero':allzero,'old_V_nonzero_checked':True,'all_four_fixed_defects_zero':allzero,
         'sharp_all_Theta_rawH_rank_upper6_if_allzero':allzero,'generic_Theta_rawH_rank6_if_allzero':allzero,
         'no_cubic57_quartic239_negative_claim':True,'stored_polynomial_terms':stored,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'native_input_sha256':sha(nativein),'native_output_sha256':sha(out/'native.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
