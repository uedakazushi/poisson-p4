#!/usr/bin/env python3
"""ONE literalcachedK7(actual68)*physicalW row; no other coefficient query."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_K7W_PURE98_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 pb=base/'S61_NONSQUARE_P100_PURE70_SOURCE_PREFLIGHT_20261104_v63';pr=read(pb/'report.json');pc=read(pb/'certificate.json');assert sha(pb/'certificate.json')==pr['certificate_sha256'];forecast=next(r for r in pc['blocks10']if r['kernel_nonpivot_genuine']==7);wprofile=next(r for r in forecast['seven_output_rows']if r['physical_light']==3);assert wprofile['pure98_output_RF_upper']==1730 and wprofile['total_live_upper_including_output_copy_and10000temps']==31782
 kb=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63';kp=kb/'pure-cubic84-7.json';kr=read(kb/'report.json');k=read(kp);binding=read(kb/'source-binding.json');assert kr['complete']and sha(kp)==kr['evidence_sha256'][kp.name]
 one={'numerator':[[[0],'1']],'denominator':[[[0],'1']]};zero={'numerator':[],'denominator':[[[0],'1']]};assert k['nonpivot_genuine_index']==7 and k['old16_source_cofactors']==[one if i==7 else zero for i in range(16)]and k['actual_delta7_source_cofactors']==k['old16_source_cofactors']
 peer=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63_REVIEW/report.json';assert sha(peer)=='aa5dd03b5317acbd2ba35f2b671db55611301d5dcf3c44c713f7b6cb3d06e7c6'
 gp=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/joint-cubic-frame.json';B=read(gp)['physical_light_inverse7x7'];assert B[3]==[one if j==6 else zero for j in range(7)]
 fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';mp=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json';fd,h,m=map(read,[fdpath,hp,mp]);assert h['all112_cubic_source319_definitions'][7*7+3]==['cubic_normalized',50,68,3]
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp';paths=[Path(__file__),cpp,utility,pb/'report.json',pb/'certificate.json',kb/'report.json',kb/'source-binding.json',kp,peer,gp,fdpath,hp,mp]
 vals=k['pure84_values'];assert len(vals)==84;sa=fd['scalar_circuit']['atoms_num_den'];sn=fd['scalar_circuit']['nodes'];ln=fd['linear_nodes'];nf=fd['full_target_projection239_to122']
 reserve=sum(len(v['numerator'])+len(v['denominator'])for v in vals)+sum(len(v['numerator'])+len(v['denominator'])for row in B for v in row)+sum(len(sa[i][0])+len(sa[i][1])for i in m['NF_needed_atom_ids']);text=[str(reserve)]
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
 source={'literal_source':'K7isexactlynormalizedC68,notacorrectednonliteralcombination; original16andactualdelta7cofactorsareunitvector7.',
         'original319_source_index':226,'source_label':['cubic_normalized',50,68,3],'physical_v3_equals_joint_W_exactly_verified':True,'joint_W_position':6,
         'old16_source_cofactors':k['old16_source_cofactors'],'actual_delta7_source_cofactors':k['actual_delta7_source_cofactors'],'pure_monomials84':binding['pure_monomials84'],
         'normalization_ROWUNIT_DIVISION_already_applied':True,'parent_kernel_peer':str(peer.relative_to(ROOT)),'coefficient_input657_RFterms':True,'output_upper1730':True,'live_upper31782':True,
         'no_H1_H2_cleanup_or_truequartic_Qhat_rescale':True,'no_other_row_source_or_Hfactor_inventory':True}
 (out/'source-binding.json').write_text(json.dumps(source,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'K7W-native';subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 a=read(out/'audit.json')if(out/'audit.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':r.returncode};complete=r.returncode==0 and a['all_passed']
 if complete:assert a['all210_joint_combination_identities_checked']and a['all98_NF_outputs_evaluated']and a['peak_live_RF_terms_with_inputreserve']<=31782 and read(out/'physical-cubic84.json')['values84']==vals
 report={'id':ID,'complete':complete,'one_literal_K7W_row_only':True,'physical_v3_to_joint_W_exact':True,'actual_source319_index226':True,'native_audit':a,'no_rank_otherNF_row_query_or_Hfactor_inventory':True,'existing149_bound_unchanged_by_coefficient_receipt':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
