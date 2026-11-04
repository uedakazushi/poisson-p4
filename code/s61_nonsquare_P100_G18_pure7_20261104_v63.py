#!/usr/bin/env python3
"""ONE seven-row G18*light pure98 block on P100; no rank or Gaussian."""
import argparse,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha,native_matrix
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_G18_PURE7_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';pf=base/'S61_NONSQUARE_P100_PURE21_SOURCE_PREFLIGHT_20261104_v63';pr=read(pf/'report.json');src=read(pf/'source-input.json');pc=read(pf/'certificate.json')
 assert sha(pf/'source-input.json')==pr['source_input_sha256']and sha(pf/'certificate.json')==pr['certificate_sha256']
 assert pc['full_source_and_all_Hzeros_verified']and pc['NF117_preserves_Hdegree_by_literal_linear_dependency']
 one=src['three_full_sources'][0];assert one['genuine_index']==0 and one['original135_equation']==18 and one['seven_full_Hlinear_coefficients_literal_zero']
 fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json'
 inv=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/coordinate-inverse-numerators7x7.txt';mp=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json'
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp'
 paths=[Path(__file__),cpp,utility,pf/'report.json',pf/'source-input.json',pf/'certificate.json',fdpath,hp,inv,mp]
 for p in[fdpath,hp,mp]:assert sha(p)==pr['input_sha256'][str(p.relative_to(ROOT))]
 fd,h,m=map(read,[fdpath,hp,mp]);_,_,n49=native_matrix(inv)
 sa=fd['scalar_circuit']['atoms_num_den'];sn=fd['scalar_circuit']['nodes'];ln=fd['linear_nodes'];nf=fd['full_target_projection239_to122']
 ps=one['full_pure84_coefficients'];assert len(ps)==84
 reserve=sum(map(len,ps))+sum(map(len,n49))+sum(len(sa[i][0])+len(sa[i][1])for i in m['NF_needed_atom_ids'])
 # Conservative live envelope: all seven outputs are NEVER retained together.
 boundparts={'literal_input_reserve':8746,'original_cubic84':2100,'joint_cubic84':7392,'inverse49':1127,'NFscalar':5488,'current_joint210':22890,'NFnode_peak_old_overestimate':26157,'current_pure98':14972,'temporaries':10000}
 assert reserve<=boundparts['literal_input_reserve']and sum(boundparts.values())==98872<150000
 serialized_bound=7*14972+2100+7392+1127;assert serialized_bound<150000
 text=[str(reserve)]
 def poly(p):text.append(str(len(p)));text.extend(f'{v} {a} {e} {d}'for(a,e,d),v in p)
 def nodes(items):
  text.append(str(len(items)))
  for i,n in items:
   text.append(f'{i} l {len(n[1])} '+' '.join(f'{c} {j}'for c,j in n[1])if n[0]=='l'else f'{i} '+' '.join(map(str,n)))
 for p in ps:poly(p)
 physical=src['physical_light_order']
 for mon in src['physical_pure_cubic_monomials84']:text.append(' '.join(str(mon[i])for i in physical))
 for p in n49:poly(p)
 text.append(str(len(m['NF_needed_atom_ids'])))
 for i in m['NF_needed_atom_ids']:text.append(str(i));poly(sa[i][0]);poly(sa[i][1])
 nodes([(i,sn[i])for i in m['NF_needed_scalar_ids']]);nodes([(i,ln[i])for i in m['NF_needed_linear_ids']])
 for ids in[m['NF_pure98_linear_roots'],m['pure98_target122_columns'],nf['input_variables']]:text.append(' '.join(map(str,ids)))
 for mon in h['parent_joint_target239_monomials']:text.append(' '.join(map(str,mon)))
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<2*1024*1024
 binding={'actual_full_source':one,'source319_rows':src['twentyone_actual_source319_products'][:7],
          'normalization':'Literal saved normalized fullTn cubics, already divided by row_unit; no extra scaling in this job.',
          'no_H1_H2_cleanup_required_for_these_literal_sources':True,'pure_input_remains_pure_by_accepted_NF117':True,
          'live_envelope_parts':boundparts,'live_envelope_sum':sum(boundparts.values()),'actual_literal_input_term_reserve':reserve,
          'serialized_terms_upper':serialized_bound,'stream_release':'Write/hash each pure98 row separately; release its210 input and98 output before next row. Only one joint cubic84 and its inverse49 are retained.',
          'no_other_two_pure_cubics_or_rankjob':True}
 (out/'source-binding.json').write_text(json.dumps(binding,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'seven-row-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 audit=read(out/'audit.json')if(out/'audit.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':r.returncode}
 complete=r.returncode==0 and audit['all_passed']
 if complete:
  assert audit['all1470_joint_combination_identities_checked']and audit['all686_NF_outputs_evaluated']and audit['peak_live_RF_terms_with_inputreserve']<=150000
  rowhashes={f'pure98-row{i}.json':sha(out/f'pure98-row{i}.json')for i in range(7)}
 else:rowhashes=None
 result={'id':ID,'complete':complete,'actual_saved_G18_times_all7_lights':True,'no_new_oracle_coefficient_query':True,'native_audit':audit,
         'seven_row_sha256':rowhashes,'no_rank_independence_full239_Gaussian_or_automatic_nextjob':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},
         'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
 signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256','seven_row_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
