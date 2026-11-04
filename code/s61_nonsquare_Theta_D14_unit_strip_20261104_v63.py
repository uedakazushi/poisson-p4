#!/usr/bin/env python3
"""Strip D14 and optional Delta by ONLY translated FULLU unit localizers."""
import argparse,itertools,json,resource,signal,subprocess,time
from pathlib import Path
from fractions import Fraction as Q
from s61_nonsquare_Theta_raw12_defect_20261104_v63 import sha,dec,enc,add,mul
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_D14_UNIT_STRIP_20261104_v63'
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';dirs=[base/n for n in['S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63','S61_NONSQUARE_THETA_RATIONAL_COORDINATES_20261104_v63','S61_NONSQUARE_THETA_DELTA_LAMBDA_PROFILES_20261104_v63','S61_NONSQUARE_THETA_REMAINING_DEFECTS_20261104_v63','S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63']]
 cpp=Path(__file__).with_suffix('.cpp');paths=[Path(__file__),cpp,Path(__file__).with_name('s61_nonsquare_Theta_raw12_defect_20261104_v63.py')]+[p/f for p in dirs for f in['report.json','certificate.json']]
 rd=lambda p:json.loads(p.read_text());data=[]
 for p in dirs:
  r,c=rd(p/'report.json'),rd(p/'certificate.json');assert sha(p/'certificate.json')==r['certificate_sha256'];data.append(c)
 g,c,d,r,t=data;B=[[dec(p)for p in row]for row in g['actual_block7x3']];rows=g['selected_rows7'];cols=g['selected_cols3'];assert rows==[0,5,2]and cols==[0,2,1]
 determinant={}
 for perm in itertools.permutations(range(3)):
  p={(0,0,0):Q(1)};inv=sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3))
  for i in range(3):p=mul(p,B[rows[i]][cols[perm[i]]])
  determinant=add(determinant,p,Q(-1)if inv%2 else Q(1))
 assert determinant==dec(g['nonzero_selected_minor_Delta'])
 D14=next(x for x in r['two_fixed_defects']if x['raw_row']==14);assert not D14['defect_zero'];delta=next(x['polynomial']for x in d['actual_images']if x['name']=='Delta')
 factors=[('zeta',[[[1,0],'1']]),('zeta_plus20',[[[1,0],'1'],[[0,0],'20']])]+[(name,c['target_coordinates'][name])for name in['L','K','Fnum','Jnum']]
 inputs=[('D14',D14['defect']),('Delta',delta)];text=['2']
 for name,p in inputs:text.append(name);text.append(str(len(p)));text.extend(f'{v} {m[0]} {m[1]}'for m,v in p)
 text.append(str(len(factors)))
 for name,p in factors:text.append(name);text.append(str(len(p)));text.extend(f'{v} {m[0]} {m[1]}'for m,v in p)
 nativein=out/'native-input.txt';nativein.write_text('\n'.join(text)+'\n');flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'strip-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20);subprocess.run([str(exe),str(nativein),str(out/'native.json')],check=True,timeout=40)
 nv=rd(out/'native.json');resultmap={x['name']:x for x in nv['results']}
 cert={'frozen_Delta_source_3x3_determinant_replayed':True,'literal_original_Delta':enc(determinant),
       'bordered_minor_binding':{'row_order_raw':[7,13,10,14],'column_order':[2,4,3,'B4_cleared_lastSchur'],
        'identity':'D14=det([B_pivots,S_pivots;B14,S14]) in this order. Follows from detB=Delta and literal Delta*B14=lambda*Bpivots by adjugate/bordered determinant identity; no extra minor selected.',
        'original_raw_heavy5_row_order':[8,7,13,10,14],'original_raw_heavy5_column_order':[0,2,4,3,5],
        'known_unit_scale':'det(rawheavy5)=(original_row8_first/B^4)*D14 moduloTheta; the row8 firstentry has the accepted known-unit inverse.'},
       'only_certified_unit_factors':factors,'native_exact_strip':nv,
       'D14_global_unit_proved':resultmap['D14']['remaining_is_nonzero_QQ_constant'],
       'Delta_global_unit_proved':resultmap['Delta']['remaining_is_nonzero_QQ_constant'],
       'no_unknown_factor_inverted_no_GB_point_prime_or_other_minor_search':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 report={'id':ID,'source_Delta3_literal_replay_verified':True,'D14_bordered_minor_identity_bound':True,
         'known_localizer_division_steps':{k:len(v['steps'])for k,v in resultmap.items()},
         'remaining_unknown_terms':{k:len(v['remaining_unknown'])for k,v in resultmap.items()},
         'D14_global_unit_proved':cert['D14_global_unit_proved'],'Delta_global_unit_proved':cert['Delta_global_unit_proved'],
         'no_unknown_factor_inversion':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
         'native_input_sha256':sha(nativein),'native_output_sha256':sha(out/'native.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
