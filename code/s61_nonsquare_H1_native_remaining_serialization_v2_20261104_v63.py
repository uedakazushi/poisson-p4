#!/usr/bin/env python3
"""Serialization repair ONLY: retain nine complete rows; query missing11 with actual byte accounting."""
import argparse,gzip,json,resource,signal,subprocess,time
from pathlib import Path
from s61_nonsquare_H1_native_remaining_quartics_20261104_v63 import read,sha,parse_factored,load_completed_raw_cache
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H1_NATIVE_REMAINING_SERIALIZATION_V2_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(49);out.mkdir(parents=True,exist_ok=False)
 first=ROOT/'reports/S61_NONSQUARE_H1_NATIVE_REMAINING_QUARTICS_20261104_v63';fr=read(first/'report.json');fc=read(first/'certificate.json');assert sha(first/'certificate.json')==fr['evidence_sha256']['certificate.json'];assert fr['new_completed_quartic_rows']==list(range(114,123))and fr['native_audit']['stop_reason']=='30MB_output_artifact_guard'
 for p,h in fr['evidence_sha256'].items():assert sha(first/p)==h
 ub=ROOT/'reports/S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63';up=ub/'certificate.json';assert sha(up)==read(ub/'report.json')['evidence_sha256'][up.name];units=read(up)
 ip=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';assert sha(ip)==units['source_input_sha256'];x=read(ip);nodes=x['scalar_circuit']['nodes'];atoms=x['scalar_circuit']['atoms_num_den'];block=x['complete_H1_image_block_scalar_roots'];missing=list(range(123,134));roots=[r for i in missing for r in block[i]]
 cp=ROOT/fc['completed_cache_descriptor']['old113_cache_path'];assert sha(cp)==fc['completed_cache_descriptor']['old113_cache_sha256'];cache=json.loads(gzip.decompress(cp.read_bytes()))
 needed=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if nodes[i][0]in['a','s','m']:todo.extend(nodes[i][1:])
 completed9roots={r for i in range(114,123)for r in block[i]};seedable=needed&completed9roots;assert not seedable
 completedroots={r for i in range(123)if i not in range(43,50)for r in block[i]};assert needed&completedroots<={0}
 aids=sorted({nodes[i][1]for i in needed if nodes[i][0]=='c'});cpp=ROOT/'code/s61_nonsquare_H1_native_remaining_serialization_v2_20261104_v63.cpp';oldcpp=ROOT/'code/s61_nonsquare_H1_native_remaining_quartics_20261104_v63.cpp';exe=ROOT/'work/H1-native-remaining-serialization-v2'
 kernel=lambda t:t.split('#include',1)[1].split('long byteswritten',1)[0]
 assert kernel(cpp.read_text())==kernel(oldcpp.read_text())
 def polyline(p):return str(len(p))+''.join(' '+c+' '+' '.join(map(str,m))for m,c in p)+'\n'
 lines=[]
 for p in units['known_factor_polynomials']:lines.append(polyline(p))
 for ai in aids:
  z=units['literal_constant_and_powers_by_denominator'][units['atomic_denominator_id_by_atom'][ai]];lines.append(str(ai)+' '+z['constant']+' '+' '.join(map(str,z['a_d_F_B_C_J_powers']))+'\n');lines.append(polyline(atoms[ai][0]))
 for i in sorted(needed):lines.append(str(i)+' '+' '.join(map(str,nodes[i]))+'\n')
 lines.append(' '.join(map(str,roots))+'\n')
 for st in cache['unit_pivots']:
  lines.append(str(st['pivot_column'])+'\n')
  for n,d in st['normalized23']:lines.extend([polyline(n),polyline(d)])
 body=''.join(lines);native=out/'native-input.txt';limit=32*1024*1024-len(body.encode())-1024*1024-256;assert limit>0
 native.write_text(f'6 {len(aids)} {len(needed)} 23 7 {len(nodes)} 11 123 {limit}\n'+body)
 command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,44-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as f:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','H1-serialization-v2',str(exe),str(native),str(out)],stdout=f,stderr=subprocess.STDOUT,timeout=max(1,45-(time.monotonic()-start)))
 audit=read(out/'native-audit.json')if(out/'native-audit.json').exists()else{'all_passed':False,'stop_reason':'native_abnormal_exit_without_audit','native_returncode':run.returncode};done=read(out/'completed-rows.json')if(out/'completed-rows.json').exists()else[];assert done==missing[:len(done)];complete=run.returncode==0 and audit['all_passed']
 combined=[]
 for item in fc['completed_cache_descriptor']['new_completed_rows']:
  r=dict(item)
  for k in['raw23_file','schur23_file','schur7coeff_file']:r[k]='../'+first.name+'/'+r[k]
  combined.append(r)
 for i in done:
  item={'source_row':i,'actual_quartic_label':x['all_remaining21_quartic_H2free_labels'][i-113]}
  for key,suffix,count in[('raw23','raw23.txt',23),('schur23','schur23.txt',23),('schur7coeff','schur7coeff.txt',7)]:
   p=out/f'row{i}-{suffix}';v=parse_factored(p);assert len(v)==count;item[key+'_file']=p.name;item[key+'_sha256']=sha(p);item[key+'_profile']={'nonzero':sum(bool(z['numerator'])for z in v),'numerator_terms':sum(len(z['numerator'])for z in v),'max_numerator_terms':max(len(z['numerator'])for z in v)}
   if key=='schur23':
    for st in cache['unit_pivots']:assert not v[st['pivot_column']]['numerator']
  combined.append(item)
 if complete:assert done==missing and audit['final_schur_identities_checked']==253 and audit['evaluated_DAG_nodes']==len(needed)
 descriptor=dict(fc['completed_cache_descriptor']);descriptor['new_completed_rows']=combined
 cert={'all134_H1_exact_raw_cache_complete':complete,'completed_exact_row_count':123+len(done),'completed_cache_descriptor':descriptor,'first_partial_packet_retained_immutable':str(first.relative_to(ROOT)),
  'first_partial_report_sha256':sha(first/'report.json'),'old9_completed_rows_recomputed':0,'available_old9_outputroot_ids_checked_for_seeding':sorted(completed9roots),'seedable_old9_outputroot_ids_in_missing11_dependencies':sorted(seedable),
  'seed_note':'All nine completed raw-output scalar IDs were checked. None occur in the missing11 dependency closure, so no nonzero seed is applicable. Shared internal values were not serialized and are explicitly recomputed for missing11.',
  'missing11_input_rows':missing,'needed_dependency_nodes':len(needed),'needed_atomic_inputs':len(aids),'all_new_dependency_nodes_were_previously_evaluated_but_not_serialized':True,
  'native_rational_arithmetic_kernel_identical_before_serialization':True,'serialization_only_change':'Count actual stream positions and next token bytes; reserve metadata/input bytes inside hard32MiB total. No 100-byte-per-term estimate.',
  'native_audit':audit,'new11_final_schur_identity_checks':23*len(done),'prior9_final_schur_identity_checks':207,'oldpsi23_final_schur_identity_checks':23,
  'full_source_schur_recipe':fc['full_source_schur_recipe'],'source134_convention':fc['source134_convention'],'pure98_queries':0,'Gaussian_jobs':0,'existing7_H1_units_and125_lower_unchanged':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');loaded=load_completed_raw_cache(out/'certificate.json');assert len(loaded)==123+len(done)
 if complete:assert set(loaded)==set(range(134))
 exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,oldcpp,ROOT/'code/s61_nonsquare_H1_native_remaining_quartics_20261104_v63.py',ROOT/'research/notes'/(ID+'.md'),first/'report.json',first/'certificate.json',up,ip,cp]
 result={'id':ID,'serialization_only_repair_all_missing11_complete':complete,'all134_H1_rows_exact_cache_complete':complete,'new_completed_rows_only':done,'old9_and_oldpsi_rows_recomputed':0,'shared_uncached_internal_dependencies_honestly_recomputed':len(needed),'native_audit':audit,
  'actual_byte_receipt':read(out/'serialization-byte-receipt.json')if(out/'serialization-byte-receipt.json').exists()else None,'Gaussian_jobs':0,'pure98_queries':0,'global7_units125_lower_and_known38_unchanged':True,'native_executable_sha256_not_committed':exesha,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'artifact_bytes_before_report':sum(p.stat().st_size for p in out.iterdir()if p.is_file()),'seconds':time.monotonic()-start,'first_plus_repair_seconds':fr['seconds']+time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert result['artifact_bytes_before_report']<32*1024*1024-30000;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
