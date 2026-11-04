#!/usr/bin/env python3
"""ONE full pure98 row, actual q5H cache0 minus its cubic H1 lift."""
import argparse,gzip,json,resource,signal,subprocess,time
from pathlib import Path
from fractions import Fraction as Q
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import sha,read,native_matrix
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63'
def mul3(p,q):
 r={}
 for m,a in p.items():
  for n,b in q.items():
   k=tuple(x+y for x,y in zip(m,n));r[k]=r.get(k,Q())+a*b
 return{k:v for k,v in r.items()if v}
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';pf=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63';pfr=read(pf/'report.json');manifest=read(pf/'manifest.json');ci=read(pf/'cubic-coefficient-input.json');assert sha(pf/'manifest.json')==pfr['manifest_sha256']and sha(pf/'cubic-coefficient-input.json')==pfr['cubic_input_sha256']
 fdpath=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';rip=base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json';cp=base/'S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz';inv=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/coordinate-inverse-numerators7x7.txt';qpair=base/'S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63/checkpoint.json';adapter=base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';cpp=Path(__file__).with_suffix('.cpp')
 paths=[Path(__file__),cpp,ROOT/'code/s61_nonsquare_P100_H1_dual_20261104_v63.cpp',pf/'report.json',pf/'manifest.json',pf/'cubic-coefficient-input.json',fdpath,hp,rip,cp,inv,qpair,adapter]
 for p in paths[6:11]:assert sha(p)==pfr['input_sha256'][str(p.relative_to(ROOT))]
 fd=read(fdpath);h=read(hp);ri=read(rip);qs=read(qpair)
 assert h['complete_H1_image_block_scalar_roots'][0]==h['remaining_q5H_H1_scalar_roots']
 with gzip.open(cp,'rt')as zf:old=json.load(zf)
 q5=next(row for i,row in old['materialized_raw_rows']if i==0);del old
 dec=lambda p:{tuple(m):Q(c)for m,c in p if Q(c)}
 classes=h['target122_Hdegree_classes'];remaining=qs['remaining_columns']
 assert all(not qs['matrix58x136'][57][remaining[j]][0]for j in classes['0']+classes['2'])
 for j,col in enumerate(classes['1']):
  n,d=qs['matrix58x136'][57][remaining[col]];nn,dd=q5[j];assert mul3(dec(n),dec(dd))==mul3(dec(nn),dec(d))
 assert read(base/'S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63/remaining-pair-block.json')['source_rows'][57]==[5,[0,0,0,0,0,0,0,1]]
 active=[i for i,row in enumerate(ri['right_inverse112x23'])if any(v['numerator']for v in row)];assert len(active)==42 and sorted({i//7 for i in active})==[1,6,8,9,11,14]
 _,_,n49=native_matrix(inv);ar=ci['needed_atoms'];on=ci['needed_nodes'];Satoms=fd['scalar_circuit']['atoms_num_den'];Snodes=fd['scalar_circuit']['nodes'];Lnodes=fd['linear_nodes'];NF=fd['full_target_projection239_to122'];text=[]
 inputterms=sum(len(n)+len(d)for _,(n,d)in ar)+sum(len(Satoms[i][0])+len(Satoms[i][1])for i in manifest['NF_needed_atom_ids'])+sum(len(n)+len(d)for n,d in q5)+sum(map(len,n49));text.append(str(inputterms))
 def poly(p):text.append(str(len(p)));text.extend(f'{c} {m[0]} {m[1]} {m[2]}'for m,c in p)
 def nodes(items):
  text.append(str(len(items)))
  for i,n in items:
   if n[0]=='l':text.append(f'{i} l {len(n[1])} '+ ' '.join(f'{a} {j}'for a,j in n[1]))
   else:text.append(f'{i} '+' '.join(map(str,n)))
 def rf(r):
  for k in['numerator','denominator']:text.append(str(len(r[k])));text.extend(f'{m[0]} {c}'for m,c in r[k])
 text.append(str(len(ar)))
 for ai,(n,d)in ar:text.append(str(ai));poly(n);poly(d)
 nodes(on)
 for r in ci['sources6']:text.append(' '.join(map(str,r['roots84'])))
 physical=[0,1,2,3,4,6,7]
 for m in ci['old_physical_monomials84']:text.append(' '.join(str(m[i])for i in physical))
 for n,d in q5:poly(n);poly(d)
 text.append(str(len(active)))
 for i in active:
  text.append(str(i))
  for r in ri['right_inverse112x23'][i]:rf(r)
 for n in n49:poly(n)
 text.append(str(len(manifest['NF_needed_atom_ids'])))
 for i in manifest['NF_needed_atom_ids']:text.append(str(i));poly(Satoms[i][0]);poly(Satoms[i][1])
 nodes([(i,Snodes[i])for i in manifest['NF_needed_scalar_ids']]);nodes([(i,Lnodes[i])for i in manifest['NF_needed_linear_ids']])
 text.append(' '.join(map(str,manifest['NF_pure98_linear_roots'])));text.append(' '.join(map(str,manifest['pure98_target122_columns'])));text.append(' '.join(map(str,NF['input_variables'])))
 for m in h['parent_joint_target239_monomials']:text.append(' '.join(map(str,m)))
 nativein=out/'native-input.txt';nativein.write_text('\n'.join(text)+'\n');assert nativein.stat().st_size<8*1024*1024
 binding={'actual_q5H_cache_row':0,'actual_qpair_source_row':57,'literal_q5H23_cache_matches_producer':True,'producer_q5H_pure98_and_H2_coefficients_identically_zero':True,
  'H2_deletion_binding':'q5H and every cubic*light source have Hdegree<=1, so separatechiH2source-deletion leaves them unchanged. Pure98 projection uses accepted NF117.',
  'sources6':ci['sources6'],'FrameAdapter_normalization':'F_actual(D*v)/basis57_source_unit (DIVISION)','root_count':504,'active_rightinverse_source_rows42':active,
  'core_input_term_reserve':inputterms,'preflight_full98_output_upper':manifest['bounds']['all98_expanded_output_terms_upper'],'no_all33quartic_or134_query':True}
 (out/'source-binding.json').write_text(json.dumps(binding,separators=(',',':'))+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'one-row-native';subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:run=subprocess.run([str(exe),str(nativein),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 audit=read(out/'audit.json')if(out/'audit.json').exists()else(read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'native_returncode':run.returncode})
 complete=run.returncode==0 and audit['all_passed'];values=read(out/'pure98.json')if complete else None
 if complete:assert audit['all210_combination_identities_checked']and audit['all98_NF117_outputs_evaluated']and audit['peak_live_RF_terms_with_inputreserve']<=150000
 result={'id':ID,'complete':complete,'actual_producer_q5H_ROW0_and_H2free_binding_verified':True,'normalized_source_cubic_coefficients504':True,'native_audit':audit,
  'nonzero_pure_positions':values['nonzero_positions']if complete else None,'known_unit_pure_positions':values['known_unit_positions']if complete else None,
  'one_source_row_only_no_rank98_full239_or_automatic_nextjob':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
  'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
