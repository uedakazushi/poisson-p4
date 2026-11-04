#!/usr/bin/env python3
"""ONE shared query of remaining20 quartics; source-root cache only for completed rows."""
import argparse,gzip,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H1_NATIVE_REMAINING_QUARTICS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse_factored(p):
 ls=p.read_text().splitlines();ans=[]
 for line in ls[1:]:
  t=line.split();v=list(map(int,t[:6]));n=int(t[6]);assert len(t)==7+4*n
  ans.append({'signed_a_d_F_B_C_J_denominator_powers':v,'numerator':[[list(map(int,t[k+1:k+4])),t[k]]for k in range(7,len(t),4)]})
 assert len(ans)==int(ls[0]);return ans
def load_completed_raw_cache(path):
 """Read saved exact134 rows without evaluating any scalar/source DAG."""
 c=read(path)['completed_cache_descriptor'];p=ROOT/c['old113_cache_path'];assert sha(p)==c['old113_cache_sha256'];x=json.loads(gzip.decompress(p.read_bytes()))
 rows={i:[{'representation':'QQ_polynomial_fraction','numerator':n,'denominator':d}for n,d in row]for i,row in x['materialized_raw_rows']}
 p=ROOT/c['psi113_raw23_path'];assert sha(p)==c['psi113_raw23_sha256'];rows[113]=[dict(v,representation='named_unit_fraction')for v in parse_factored(p)]
 for r in c['new_completed_rows']:
  p=path.parent/r['raw23_file'];assert sha(p)==r['raw23_sha256'];rows[r['source_row']]=[dict(v,representation='named_unit_fraction')for v in parse_factored(p)]
 return rows
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 ub=ROOT/'reports/S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63';up=ub/'certificate.json';assert sha(up)==read(ub/'report.json')['evidence_sha256'][up.name];units=read(up)
 base=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';ip=base/'certificate.json';assert sha(ip)==units['source_input_sha256'];x=read(ip);nodes=x['scalar_circuit']['nodes'];atoms=x['scalar_circuit']['atoms_num_den'];block=x['complete_H1_image_block_scalar_roots'];newrows=list(range(114,134));roots=[r for i in newrows for r in block[i]]
 cp=ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz';assert sha(cp)==read(cp.parent/'report.json')['evidence_sha256'][cp.name];cache=json.loads(gzip.decompress(cp.read_bytes()));assert len(cache['materialized_raw_rows'])==113 and cache['completed_pivots']==7
 rv=ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63_REVIEW/report.json';assert read(rv)['all_passed']
 psi=ROOT/'reports/S61_NONSQUARE_H1_NATIVE_PSI_QUERY_20261104_v63';pr=read(psi/'report.json');assert pr['completed_psi113_query_only'];psip=psi/'raw23-factored.txt';assert sha(psip)==pr['evidence_sha256'][psip.name]
 def closure(rs):
  needed=set();todo=list(rs)
  while todo:
   i=todo.pop()
   if i in needed:continue
   needed.add(i)
   if nodes[i][0]in['a','s','m']:todo.extend(nodes[i][1:])
  return needed
 needed=closure(roots);aids=sorted({nodes[i][1]for i in needed if nodes[i][0]=='c'});assert len(needed)<=119695 and len(aids)<=1198
 oldroots={r for i in range(113)if i not in range(43,50)for r in block[i]}|set(block[113]);boundary=needed&oldroots;assert boundary<={0}
 native=out/'native-input.txt';cpp=ROOT/'code/s61_nonsquare_H1_native_remaining_quartics_20261104_v63.cpp';oldcpp=ROOT/'code/s61_nonsquare_H1_native_psi_query_20261104_v63.cpp';exe=ROOT/'work/H1-native-remaining-quartics'
 kernel=lambda t:t.split('#include',1)[1].split('int main(',1)[0]
 assert kernel(cpp.read_text())==kernel(oldcpp.read_text())
 def polyline(p):return str(len(p))+''.join(' '+c+' '+' '.join(map(str,m))for m,c in p)+'\n'
 with native.open('w')as f:
  f.write(f'6 {len(aids)} {len(needed)} 23 7 {len(nodes)} 20\n')
  for p in units['known_factor_polynomials']:f.write(polyline(p))
  for ai in aids:
   z=units['literal_constant_and_powers_by_denominator'][units['atomic_denominator_id_by_atom'][ai]];f.write(str(ai)+' '+z['constant']+' '+' '.join(map(str,z['a_d_F_B_C_J_powers']))+'\n');f.write(polyline(atoms[ai][0]))
  for i in sorted(needed):f.write(str(i)+' '+' '.join(map(str,nodes[i]))+'\n')
  f.write(' '.join(map(str,roots))+'\n')
  for st in cache['unit_pivots']:
   f.write(str(st['pivot_column'])+'\n')
   for n,d in st['normalized23']:f.write(polyline(n));f.write(polyline(d))
 command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as f:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','H1-native',str(exe),str(native),str(out)],stdout=f,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-start)))
 audit=read(out/'native-audit.json')if(out/'native-audit.json').exists()else{'all_passed':False,'stop_reason':'native_abnormal_exit_without_complete_audit','native_returncode':run.returncode}
 completed=read(out/'completed-rows.json')if(out/'completed-rows.json').exists()else[];assert completed==newrows[:len(completed)];complete=run.returncode==0 and audit['all_passed'];files=[]
 for i in completed:
  item={'source_row':i,'actual_quartic_label':x['all_remaining21_quartic_H2free_labels'][i-113]}
  for key,ending,count in[('raw23','raw23.txt',23),('schur23','schur23.txt',23),('schur7coeff','schur7coeff.txt',7)]:
   p=out/f'row{i}-{ending}';v=parse_factored(p);assert len(v)==count;item[key+'_file']=p.name;item[key+'_sha256']=sha(p);item[key+'_profile']={'nonzero':sum(bool(z['numerator'])for z in v),'numerator_terms':sum(len(z['numerator'])for z in v),'max_numerator_terms':max(len(z['numerator'])for z in v)}
   if key=='schur23':
    for st in cache['unit_pivots']:assert not v[st['pivot_column']]['numerator']
  files.append(item)
 if complete:assert completed==newrows and audit['evaluated_DAG_nodes']==len(needed)and audit['old7_input_fraction_crossproducts']==161 and audit['final_schur_identities_checked']==460
 descriptor={'old113_cache_path':str(cp.relative_to(ROOT)),'old113_cache_sha256':sha(cp),'psi113_raw23_path':str(psip.relative_to(ROOT)),'psi113_raw23_sha256':sha(psip),'new_completed_rows':files}
 cert={'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','all134_H1_exact_raw_cache_complete':complete,'completed_exact_row_count':114+len(completed),'completed_cache_descriptor':descriptor,
  'new_native_input_rows':newrows,'new_input_scalar_root_rows20x23':[block[i]for i in newrows],'needed_DAG_nodes':len(needed),'needed_atoms':len(aids),'completed_output_root_intersection_with_new_dependencies':sorted(boundary),'shared_internal_ancestors_with_prior_psi_graph':len(needed&closure(block[113])),
  'shared_ancestor_note':'No completed nonzero source output root is evaluated. Shared internal ancestors without serialized values are explicit dependencies of the new20 expressions, not claimed unavailable cache.',
  'native_arithmetic_kernel_identical_to_frozen_psi':True,'native_audit':audit,
  'source134_convention':'Same saved delta7 source change on cubic rows43..49; quartic113..133 unaffected. Cache0..112 is saved113 matrix,113 is frozennativepsi,114..133 are new values.',
  'full_source_schur_recipe':'Each reducedfullsource = originalfullsource - sum_i c_i*oldnormalizedfullpivot_i. Saved7 cofs and all23 restored identities bind the operation; pure98 remains the same lazy full-source arithmetic.',
  'old113_source_coefficient_root_queries':0,'psi113_coefficient_root_queries':0,'Gaussian_jobs':0,'pure98_queries':0,'all21_quartic_sources_retained':True,'known38_missing42_and7_H1_units_unchanged':True,'no_full121_full239_or_new_unit_rank_claim':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');loaded=load_completed_raw_cache(out/'certificate.json');assert len(loaded)==114+len(completed)
 if complete:assert set(loaded)==set(range(134))
 exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,oldcpp,ROOT/'research/notes'/(ID+'.md'),ub/'report.json',up,base/'report.json',ip,cp,rv,psi/'report.json',psip]
 result={'id':ID,'all_remaining20_quartics_complete':complete,'all134_exact_H1_raw_rows_complete':complete,'new_completed_quartic_rows':completed,'native_audit':audit,'no_old113_or_psi_outputroot_requery':True,'shared_internal_ancestors_explicitly_recorded':True,'Gaussian_jobs':0,'pure98_queries':0,'global7_H1_units_and125_lower_unchanged':True,'complete134_cache_reader_api_verified':complete,
  'native_executable_sha256_not_committed':exesha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
