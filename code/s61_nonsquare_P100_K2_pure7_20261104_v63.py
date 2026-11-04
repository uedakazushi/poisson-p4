#!/usr/bin/env python3
"""ONE cached K2(actual nonpivot23) times7 light block; no rank/series."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_K2_PURE7_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 pf=base/'S61_NONSQUARE_P100_PURE70_SOURCE_PREFLIGHT_20261104_v63';pr=read(pf/'report.json');pc=read(pf/'certificate.json');assert sha(pf/'certificate.json')==pr['certificate_sha256']
 forecast=next(r for r in pc['blocks10']if r['kernel_nonpivot_genuine']==2);assert forecast['all7_retained_output_RF_upper']==45384 and forecast['max_live_upper']==53925
 kb=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63';kp=kb/'pure-cubic84-2.json';kr=read(kb/'report.json');k=read(kp);binding=read(kb/'source-binding.json')
 assert kr['complete']and sha(kp)==kr['evidence_sha256'][kp.name]and k['both_full_source_frame_identities_verified']and k['all7_Hlinear_coefficients_zero']
 peer=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63_REVIEW/report.json';assert sha(peer)=='aa5dd03b5317acbd2ba35f2b671db55611301d5dcf3c44c713f7b6cb3d06e7c6'
 gp=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/joint-cubic-frame.json';B=read(gp)['physical_light_inverse7x7']
 fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';mp=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json'
 fd,h,m=map(read,[fdpath,hp,mp]);cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp'
 paths=[Path(__file__),cpp,utility,pf/'report.json',pf/'certificate.json',kb/'report.json',kb/'source-binding.json',kp,peer,gp,fdpath,hp,mp]
 for p in[gp,fdpath,hp]:assert sha(p)==pr['input_sha256'][str(p.relative_to(ROOT))]
 vals=k['pure84_values'];assert len(vals)==84 and k['nonpivot_genuine_index']==2;sa=fd['scalar_circuit']['atoms_num_den'];sn=fd['scalar_circuit']['nodes'];ln=fd['linear_nodes'];nf=fd['full_target_projection239_to122']
 reserve=sum(len(v['numerator'])+len(v['denominator'])for v in vals)+sum(len(v['numerator'])+len(v['denominator'])for row in B for v in row)+sum(len(sa[i][0])+len(sa[i][1])for i in m['NF_needed_atom_ids'])
 text=[str(reserve)]
 def rf(v):
  for part in['numerator','denominator']:text.append(str(len(v[part])));text.extend(f'{mm[0]} {c}'for mm,c in v[part])
 def poly(p):text.append(str(len(p)));text.extend(f'{c} {a} {e} {d}'for(a,e,d),c in p)
 def nodes(items):
  text.append(str(len(items)))
  for i,n in items:text.append(f'{i} l {len(n[1])} '+' '.join(f'{c} {j}'for c,j in n[1])if n[0]=='l'else f'{i} '+' '.join(map(str,n)))
 for v in vals:rf(v)
 for mon in binding['pure_monomials84']:text.append(' '.join(str(mon[j])for j in[0,1,2,3,4,6,7]))
 for row in B:
  for v in row:rf(v)
 text.append(str(len(m['NF_needed_atom_ids'])))
 for i in m['NF_needed_atom_ids']:text.append(str(i));poly(sa[i][0]);poly(sa[i][1])
 nodes([(i,sn[i])for i in m['NF_needed_scalar_ids']]);nodes([(i,ln[i])for i in m['NF_needed_linear_ids']])
 for ids in[m['NF_pure98_linear_roots'],m['pure98_target122_columns'],nf['input_variables']]:text.append(' '.join(map(str,ids)))
 for mon in h['parent_joint_target239_monomials']:text.append(' '.join(map(str,mon)))
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<2*1024*1024
 sb={'kernel_label':2,'actual_nonpivot135':23,'actual_full_source_combination':'K2=C23-sum_p alpha2,p*C_p in saved normalized Tn frame; it is not the literal uncorrected C23.',
     'old16_source_cofactors':k['old16_source_cofactors'],'actual_delta7_source_cofactors':k['actual_delta7_source_cofactors'],
     'physical_lights':[0,1,2,3,4,6,7],'original319_source_product_recipe':'For lightpositionj use old16 cofactors at source319 indices174+7*s+j. For preconditionedH1cache only use actual_delta7 cofactors at cacheindices1+7*s+j.',
     'pure_monomials84':binding['pure_monomials84'],'parent_whole_curve_kernel_peer':str(peer.relative_to(ROOT)),'all_Hlinear_H2_zero_and_normalization_rowunit_DIVISION_already_verified':True,
     'inverse_light_frame_source':str(gp.relative_to(ROOT)),'NF117_Hdegree_closure_source':str(hp.relative_to(ROOT)),
     'preflight_7outputs_upper':45384,'preflight_peak_live_upper':53925,'actual_literal_input_reserve':reserve,'stream_each_row_and_release':True,'no_other_kernel_quartic_series_rank_or_oracle':True}
 (out/'source-binding.json').write_text(json.dumps(sb,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'K2-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 a=read(out/'audit.json')if(out/'audit.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':r.returncode};complete=r.returncode==0 and a['all_passed']
 if complete:
  assert a['all1470_joint_combination_identities_checked']and a['all686_NF_outputs_evaluated']and a['peak_live_RF_terms_with_inputreserve']<=53925
  assert read(out/'physical-cubic84.json')['values84']==vals
  rowhashes={f'pure98-row{i}.json':sha(out/f'pure98-row{i}.json')for i in range(7)}
 else:rowhashes=None
 report={'id':ID,'complete':complete,'only_cached_K2_times7':True,'parent_kernel_source_and_84_values_bound':True,'native_audit':a,'seven_row_sha256':rowhashes,'no_rank_series_or_automatic_nextjob':True,'whole142_generic149_bounds_unchanged_by_this_coefficient_receipt':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({kk:v for kk,v in report.items()if kk not in['input_sha256','evidence_sha256','seven_row_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
