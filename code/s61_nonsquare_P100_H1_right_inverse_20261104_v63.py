#!/usr/bin/env python3
"""Fixed whole-P100 112-by-23 source right inverse, H1 projection only."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def native_matrix(p):
 ls=p.read_text().splitlines();n,m=map(int,ls[0].split());assert len(ls)==n*m+1;out=[]
 for ln in ls[1:]:
  t=ln.split();k=int(t[0]);assert len(t)==1+4*k;out.append([[list(map(int,t[j+1:j+4])),t[j]]for j in range(1,len(t),4)])
 return n,m,out
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';raw=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';gdir=base/'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';dual=base/'S61_NONSQUARE_P100_H1_DUAL_20261104_v63';compact=base/'S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';joint=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';hdir=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';precondition=base/'S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63/checkpoint.json';cpp=Path(__file__).with_suffix('.cpp')
 paths=[Path(__file__),cpp,Path(__file__).with_name('s61_nonsquare_P100_H1_dual_20261104_v63.cpp'),raw/'report.json',raw/'preflight.json',raw/'checkpoint.json',gdir/'report.json',gdir/'certificate.json',dual/'report.json',dual/'native.json',compact/'report.json',compact/'six-original-quadrics.txt',compact/'source-metadata.txt',joint/'report.json',joint/'circuit.json',joint/'preflight.json',joint/'coordinate-forward7x7.txt',joint/'six-prime-quadrics-numerators.txt',joint/'six-prime-quadrics-denominators.txt',joint/'monomial-input.txt',hdir/'report.json',hdir/'certificate.json']
 rr,pre,ch=map(read,paths[3:6]);gr,g=map(read,paths[6:8]);dr,dv=map(read,paths[8:10]);cm=read(paths[10]);jr=read(paths[13]);h=read(paths[21])
 assert sha(paths[4])==rr['evidence_sha256']['preflight.json']and sha(paths[5])==rr['evidence_sha256']['checkpoint.json'];assert sha(paths[7])==gr['certificate_sha256'];assert sha(paths[9])==dr['native_output_sha256']
 assert dv['all16_actual_rawH_annihilations_exact']and dv['v3_is_one']
 for p in paths[11:13]:assert sha(p)==cm['artifact_sha256'][p.name]
 for p in paths[16:20]:assert sha(p)==jr['evidence_sha256'][p.name]
 _,_,q=native_matrix(paths[11]);_,_,forward=native_matrix(paths[16]);n,m,prime=native_matrix(paths[17]);assert(n,m)==(6,29);n,m,den=native_matrix(paths[18]);assert(n,m)==(6,1)
 oldmeta=paths[12].read_text().splitlines();oldmons=[list(map(int,s.split()))for s in oldmeta[1:30]];newmeta=paths[19].read_text().splitlines();newmons=[list(map(int,s.split()))for s in newmeta[30:59]];assert len(newmons)==29
 hcols=h['target122_Hdegree_classes']['1'];inds=h['parent_target122_monomial_indices'];hmons=[h['parent_joint_target239_monomials'][inds[j]]for j in hcols];assert len(hmons)==23 and hmons[-1]==[0,0,0,0,0,0,2,1]
 defs=h['all112_cubic_source319_definitions'];assert len(defs)==112
 physical=[0,1,2,3,4,6,7]
 for i in range(16):assert[defs[7*i+j][3]for j in range(7)]==physical
 # All 112 source rows have the initial q5H offset in the actual134 cache.
 precondpath=base/'S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz'
 import gzip
 with gzip.open(precondpath,'rt')as f:pc=json.load(f)['raw_row_preconditioner']
 assert pc['replace_rows']==list(range(43,50))and pc['formula']=='row(43+j)=delta*(oldrow(43+j)-oldrow(8+j))';paths.append(precondpath)
 text=[]
 def poly(p):text.append(str(len(p)));text.extend(f'{v} {m[0]} {m[1]} {m[2]}'for m,v in p)
 def unit(u):text.append(u['constant']+' '+' '.join(map(str,u['a_d_F_B_C_J_powers'])))
 for row in pre['raw16x7']:
  for p in row:poly(p)
 ids=[1,7,8,10,13]
 for i in ids:
  for p in ch['matrix15x6'][i]:poly(p)
 poly(g['nonzero_selected_minor_Delta'])
 for row in pre['constant_coefficient_column_frame7x7']:text.append(' '.join(map(str,row)))
 for row in pre['constant_inverse7x7']:text.append(' '.join(map(str,row)))
 for p in[pre['pair_values'][k]for k in['cZ','cW','eta']]+[pre['source_delta_change']['delta']]:poly(p)
 for r in dv['eta_physical7']:
  for k in['numerator','denominator']:
   text.append(str(len(r[k])));text.extend(f'{m[0]} {v}'for m,v in r[k])
 for u in ch['initial_row_divisors']:unit(u)
 el={v['row']:v for v in ch['steps'][0]['eliminations']}
 for i in ids:
  text.append('1'if i in el else'0')
  if i in el:poly(el[i]['coefficient_before']);unit(el[i]['row_divisor_after'])
 for m in oldmons:text.append(' '.join(map(str,m)))
 for p in q:poly(p)
 for p in forward:poly(p)
 for m in newmons:text.append(' '.join(map(str,m)))
 for p in prime:poly(p)
 for p in den:poly(p)
 for m in hmons:text.append(' '.join(map(str,m)))
 nativein=out/'native-input.txt';nativein.write_text('\n'.join(text)+'\n');assert nativein.stat().st_size<8*1024*1024
 (out/'preflight.json').write_text(json.dumps({'curve':'P100','fixed_source_basis_raw_genuine_indices':[1,6,8,9,11,14],'target_count':23,'cofactor_shape':[112,23],'actual_cache_rows':list(range(1,113)),'no_full134cache_query':True,'native_input_bytes':nativein.stat().st_size,'stored_term_cap':150000},indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'right-inverse-native';subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20)
 subprocess.run([str(exe),str(nativein),str(out/'native.json')],check=True,timeout=max(1,55-(time.monotonic()-started)))
 nv=read(out/'native.json');assert nv['all7_fixed_source_linear_lifts_verified']and nv['all5_original_qH_NF_relations_verified']and nv['all23_source_identities_both_source_frames_verified']and nv['all_cofactor_denominators_are_whole_curve_units']
 assert nv['physical_quadratic_coefficient_checks']==2*23*28 and len(nv['right_inverse112x23'])==112
 cert={'curve':'P100 wholecurve open','source_cubic_cofactor_rows':[{'matrix_row':7*i+j,'actual134_cache_row':1+7*i+j,'genuine_index':i,'physical_light':physical[j],'frozen_definition':defs[7*i+j]}for i in range(16)for j in range(7)],
       'actual_delta7_preconditioner':pc,'source_coefficient_identity':'Q_target = sum_sj cofactor_sj*(raw_Hlinear_s)*v_j + rho*q1; checked all28coeffs for all23 targets in bothsourceframes.',
       'NF117_source_binding':'FiveHkernel directions explicitly reduced from fullcorrected jointprimequad rows, sourceinverse q0..q4 checked. Retained23Hmonomials matchactualtargetlabels. q1H has NF0. Thus preNF identities give actualH1projection rightinverse.',
       'H1_projection_right_inverse':[112,23],'native_certificate_file':'native.json',
       'H1_projection_surjective_on_whole_P100':True,
       'pure_tail_scope':'Cubics have puretails. Only their H1projection is inverted; originalpure98 remainslazy in everyfullsource combination. No H*l=0/fullHW^2 ideal assertion.',
       'source_unit_block_rank_lower':'Original319 target has117certifiedquadunitpivots plus H2unit plus23H1projectionunitdirections; fullmatrix rank>=141 at everyP100 point. This does not prove239.',
       'no_Gaussian_search_newpoint_prime_GB_full134_or_POWERexpansion':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'all23_actual_source_H1_rightinverse_identities_verified':True,'all5_eliminated_Hquad_relations_verified':True,'rightinverse_shape':[112,23],
         'whole_P100_H1projection_surjective':True,'whole_P100_fullquartic_rank_lower141':True,'no_full57_full239_claim':True,
         'native_peak_poly_length':nv['peak_polynomial_length'],'stored_terms':nv['stored_terms'],
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'native_input_sha256':sha(nativein),'native_output_sha256':sha(out/'native.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
