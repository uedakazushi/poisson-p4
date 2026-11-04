#!/usr/bin/env python3
"""One whole P80 dual/gate and, only for a certified unit, explicit H1 right inverse."""
import argparse,gzip,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import sha,read,native_matrix
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P80_H1_RESCUE_20261104_v63'
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';raw=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';gdir=base/'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';compact=base/'S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';joint=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';hdir=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';curves=base/'S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63';cp=base/'S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz'
 kernels=[Path(__file__).with_name('s61_nonsquare_P80_H1_dual_kernel_20261104_v63.cpp'),Path(__file__).with_name('s61_nonsquare_P80_H1_right_inverse_kernel_20261104_v63.cpp')]
 paths=[Path(__file__),*kernels,Path(__file__).with_name('s61_nonsquare_P100_H1_right_inverse_20261104_v63.py'),raw/'report.json',raw/'preflight.json',raw/'checkpoint.json',gdir/'report.json',gdir/'certificate.json',compact/'report.json',compact/'six-original-quadrics.txt',compact/'source-metadata.txt',joint/'report.json',joint/'circuit.json',joint/'coordinate-forward7x7.txt',joint/'six-prime-quadrics-numerators.txt',joint/'six-prime-quadrics-denominators.txt',joint/'monomial-input.txt',hdir/'report.json',hdir/'certificate.json',curves/'report.json',curves/'certificate.json',cp]
 rr,pre,ch=map(read,paths[4:7]);gr,g=map(read,paths[7:9]);cm=read(paths[9]);jr=read(paths[12]);h=read(paths[19]);cr,cc=map(read,paths[20:22])
 assert sha(paths[5])==rr['evidence_sha256']['preflight.json']and sha(paths[6])==rr['evidence_sha256']['checkpoint.json']and sha(paths[8])==gr['certificate_sha256']
 assert sha(paths[21])==cr['certificate_sha256']and any(q['m']==80 and q['nonempty_rational_open_curve']for q in cc['curve_profiles'])
 for p in paths[10:12]:assert sha(p)==cm['artifact_sha256'][p.name]
 for p in paths[14:18]:assert sha(p)==jr['evidence_sha256'][p.name]
 _,_,q=native_matrix(paths[10]);_,_,forward=native_matrix(paths[14]);n,m,prime=native_matrix(paths[15]);assert(n,m)==(6,29);n,m,den=native_matrix(paths[16]);assert(n,m)==(6,1)
 oldmons=[list(map(int,s.split()))for s in paths[11].read_text().splitlines()[1:30]];newmons=[list(map(int,s.split()))for s in paths[17].read_text().splitlines()[30:59]]
 hcols=h['target122_Hdegree_classes']['1'];inds=h['parent_target122_monomial_indices'];hmons=[h['parent_joint_target239_monomials'][inds[j]]for j in hcols];assert len(hmons)==23
 defs=h['all112_cubic_source319_definitions'];physical=[0,1,2,3,4,6,7]
 assert all([defs[7*i+j][3]for j in range(7)]==physical for i in range(16))
 with gzip.open(cp,'rt')as f:pc=json.load(f)['raw_row_preconditioner']
 assert pc['replace_rows']==list(range(43,50))and pc['formula']=='row(43+j)=delta*(oldrow(43+j)-oldrow(8+j))'
 assert {tuple(m):v for m,v in pc['delta'][1]}=={(0,0,0):'1'}and{tuple(m):v for m,v in pc['delta'][0]}=={tuple(m):v for m,v in pre['source_delta_change']['delta']}
 ids=[1,7,8,10,13];flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();nativefiles=[]
 def run(kernel,lines,name):
  p=out/(name+'-input.txt');p.write_text('\n'.join(lines)+'\n');assert p.stat().st_size<8*1024*1024;exe=out/(name+'-native');subprocess.run(['g++','-O2','-std=c++17',str(kernel),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-started)));target=out/(name+'.json');subprocess.run([str(exe),str(p),str(target)],check=True,timeout=max(1,55-(time.monotonic()-started)));nativefiles.extend([p,target]);return read(target)
 def putpoly(lines,p):lines.append(str(len(p)));lines.extend(f'{v} {m[0]} {m[1]} {m[2]}'for m,v in p)
 lines=['16 7']
 for row in pre['raw16x7']:
  for p in row:putpoly(lines,p)
 lines.append('1 7 8 10 13')
 for i in ids:
  for p in ch['matrix15x6'][i]:putpoly(lines,p)
 putpoly(lines,g['nonzero_selected_minor_Delta'])
 for row in pre['constant_coefficient_column_frame7x7']:lines.append(' '.join(map(str,row)))
 for m in oldmons:lines.append(' '.join(map(str,m)))
 for p in q:putpoly(lines,p)
 (out/'preflight.json').write_text(json.dumps({'curve':'P80','raw_H_core_terms':sum(len(p)for row in pre['raw16x7']for p in row),'quad6_core_terms':sum(map(len,q)),'conditional_rightinverse_shape':[112,23],'actual_source_cache_rows':[1,112],'no_q5H_or_quartic_cachequery':True},indent=2)+'\n')
 dv=run(kernels[0],lines,'dual');assert dv['v3_is_one']and dv['all16_actual_rawH_annihilations_exact']
 units=dv['quadratic_unit_flags'];candidates=[1,0,2,3,4];chosen=next((i for i in candidates if units[i]),None);nv=None
 if chosen is not None:
  lines=[]
  for row in pre['raw16x7']:
   for p in row:putpoly(lines,p)
  for i in ids:
   for p in ch['matrix15x6'][i]:putpoly(lines,p)
  putpoly(lines,g['nonzero_selected_minor_Delta'])
  for H in[pre['constant_coefficient_column_frame7x7'],pre['constant_inverse7x7']]:
   for row in H:lines.append(' '.join(map(str,row)))
  for p in[pre['pair_values'][k]for k in['cZ','cW','eta']]+[pre['source_delta_change']['delta']]:putpoly(lines,p)
  def putr(r):
   for k in['numerator','denominator']:lines.append(str(len(r[k])));lines.extend(f'{m[0]} {v}'for m,v in r[k])
  for r in dv['eta_physical7']:putr(r)
  def unit(u):lines.append(u['constant']+' '+' '.join(map(str,u['a_d_F_B_C_J_powers'])))
  for u in ch['initial_row_divisors']:unit(u)
  el={v['row']:v for v in ch['steps'][0]['eliminations']}
  for i in ids:
   lines.append('1'if i in el else'0')
   if i in el:putpoly(lines,el[i]['coefficient_before']);unit(el[i]['row_divisor_after'])
  for m in oldmons:lines.append(' '.join(map(str,m)))
  for p in q:putpoly(lines,p)
  for p in forward:putpoly(lines,p)
  for m in newmons:lines.append(' '.join(map(str,m)))
  for p in prime:putpoly(lines,p)
  for p in den:putpoly(lines,p)
  for m in hmons:lines.append(' '.join(map(str,m)))
  lines.append(str(chosen));putr(dv['quadratic_evaluations'][chosen]);nv=run(kernels[1],lines,'right-inverse')
  assert nv['all7_fixed_source_linear_lifts_verified']and nv['all5_original_qH_NF_relations_verified']and nv['all23_source_identities_both_source_frames_verified']and nv['all_cofactor_denominators_are_whole_curve_units']and nv['physical_quadratic_coefficient_checks']==1288
 cert={'curve':'wholeP80open','dual_results_file':'dual.json','raw16_annihilations_verified':True,'six_q_evaluations_sourcebound':True,
       'chosen_eliminated_unit_quadric':chosen,'rightinverse_file':'right-inverse.json'if nv else None,'actual_rightinverse_rows112':'cache rows1..112, q5Hrow0 notused',
       'source_preconditioner':pc,'source_definitions112':defs,
       'cofactor_identity':'Q_j=sum112c*(rawHlinear_s)*v_k+rho*q_chosen; all23*28 checked beforeNF, both original anddelta7frames. Fullcorrected fiveqH NFrelations checked. H1projectionONLY.',
       'whole_curve_H1_surjective':nv is not None,'whole_curve_fullquartic_rank_lower141':nv is not None,
       'pure98_and_all21quartics_retained_lazy':True,'no_full134Gaussian_POWER_newpoint_prime_GB_or_full57_full239':True,
       'if_no_unit_scope':'No whole-curve inverse claimed; exact q values remain and finite nonunit zero sets require further work.'}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'all16_raw_H_annihilations_verified':True,'all6_original_q_evaluations_verified':True,
         'quadratic_zero_flags':dv['quadratic_zero_flags'],'quadratic_unit_flags':units,'chosen_eliminated_unit_quadric':chosen,
         'wholeP80_H1_rightinverse112x23_verified':nv is not None,'wholeP80_fullquartic_rank_lower141':nv is not None,
         'native_stored_terms':dv['stored_output_terms']+(nv['stored_terms']if nv else 0),'native_peak_poly_length':max(dv['peak_polynomial_length'],nv['peak_polynomial_length']if nv else 0),
         'no_cubic57_full239_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'native_evidence_sha256':{p.name:sha(p)for p in nativefiles},'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert result['native_stored_terms']<=150000 and sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
