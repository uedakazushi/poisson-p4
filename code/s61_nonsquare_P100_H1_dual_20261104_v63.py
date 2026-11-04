#!/usr/bin/env python3
"""One P100 raw-H dual and six q evaluations; at most one conditional H1 row."""
import argparse,gzip,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_H1_DUAL_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def native_matrix(p):
 ls=p.read_text().splitlines();n,m=map(int,ls[0].split());assert len(ls)==n*m+1;out=[]
 for ln in ls[1:]:
  v=ln.split();k=int(v[0]);assert len(v)==1+4*k;out.append([[list(map(int,v[j+1:j+4])),v[j]]for j in range(1,len(v),4)])
 return n,m,out
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';rawdir=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';gdir=base/'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';curves=base/'S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63';compact=base/'S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';joint=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';hdir=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';cachepath=base/'S61_NONSQUARE_H1_NATIVE_REMAINING_SERIALIZATION_V2_20261104_v63/certificate.json'
 cpp=Path(__file__).with_suffix('.cpp');paths=[Path(__file__),cpp,rawdir/'report.json',rawdir/'preflight.json',rawdir/'checkpoint.json',gdir/'report.json',gdir/'certificate.json',curves/'report.json',curves/'certificate.json',curves/'operator-recipe.json',compact/'report.json',compact/'six-original-quadrics.txt',compact/'source-metadata.txt',joint/'report.json',joint/'circuit.json',joint/'coordinate-forward7x7.txt',hdir/'report.json',hdir/'certificate.json',cachepath]
 rr,pre,check=map(read,paths[2:5]);gr,g=map(read,paths[5:7]);cr,curve,op=map(read,paths[7:10]);cm=read(paths[10]);jr=read(paths[13]);h=read(paths[17]);cache=read(cachepath)
 assert sha(paths[3])==rr['evidence_sha256']['preflight.json']and sha(paths[4])==rr['evidence_sha256']['checkpoint.json']
 assert sha(paths[6])==gr['certificate_sha256']and sha(paths[8])==cr['certificate_sha256']and sha(paths[9])==cr['operator_recipe_sha256']
 for p in paths[11:13]:assert sha(p)==cm['artifact_sha256'][p.name]
 assert sha(paths[15])==jr['evidence_sha256'][paths[15].name]
 n,m,q=native_matrix(paths[11]);assert(n,m)==(6,29)
 n,m,forward=native_matrix(paths[15]);assert(n,m)==(7,7)
 meta=paths[12].read_text().splitlines();assert meta[0]=='29 91 239';mons=[list(map(int,s.split()))for s in meta[1:30]]
 inds=h['parent_target122_monomial_indices'];hm=[h['parent_joint_target239_monomials'][inds[j]]for j in h['target122_Hdegree_classes']['1']];assert len(hm)==23 and hm[22]==[0,0,0,0,0,0,2,1]
 descriptor=cache['completed_cache_descriptor'];oldcache=ROOT/descriptor['old113_cache_path'];assert sha(oldcache)==descriptor['old113_cache_sha256'];paths.append(oldcache)
 with gzip.open(oldcache,'rt')as f:old=json.load(f)
 q5=next(row for i,row in old['materialized_raw_rows']if i==112);assert len(q5)==23;del old
 entry=next(r for r in descriptor['new_completed_rows']if r['source_row']==118);assert entry['actual_quartic_label']==55
 quarticpath=(cachepath.parent/entry['raw23_file']).resolve();assert sha(quarticpath)==entry['raw23_sha256'];paths.append(quarticpath)
 lines=quarticpath.read_text().splitlines();assert int(lines[0])==23;quartic=[]
 for line in lines[1:]:
  t=line.split();powers=list(map(int,t[:6]));n=int(t[6]);assert len(t)==7+4*n
  p=[[list(map(int,t[i+1:i+4])),t[i]]for i in range(7,len(t),4)];quartic.append((powers,p))
 assert len(quartic)==23 and sum(len(p)for _,p in quartic)==425
 counts={'raw16_H_coeff_terms':sum(len(p)for row in pre['raw16x7']for p in row),'quad6_core_terms':sum(map(len,q)),
         'selected_checkpoint5_rows_terms':sum(len(p)for i in[1,7,8,10,13]for p in check['matrix15x6'][i]),
         'q5H_row112_numden_terms':sum(len(n)+len(d)for n,d in q5),'quartic118_actual55_terms':425,
         'one_curve_only':'P100','optional_rows_are_conditional_only':'Use cachedq5H112 if directq5 evaluation nonzero; otherwise onlyH2freequartic118. IfNFdescent fails neitherquery is evaluated.'}
 assert sum(v for v in counts.values()if isinstance(v,int))<25000
 (out/'preflight.json').write_text(json.dumps(counts,indent=2)+'\n')
 text=[]
 def poly(p):
  text.append(str(len(p)));text.extend(f'{v} {m[0]} {m[1]} {m[2]}'for m,v in p)
 text.append('16 7')
 for row in pre['raw16x7']:
  for p in row:poly(p)
 text.append('1 7 8 10 13')
 for i in[1,7,8,10,13]:
  for p in check['matrix15x6'][i]:poly(p)
 poly(g['nonzero_selected_minor_Delta'])
 for row in pre['constant_coefficient_column_frame7x7']:text.append(' '.join(map(str,row)))
 for v in mons:text.append(' '.join(map(str,v)))
 for p in q:poly(p)
 for p in forward:poly(p)
 for v in hm:text.append(' '.join(map(str,v)))
 for n,d in q5:poly(n);poly(d)
 for powers,p in quartic:text.append(' '.join(map(str,powers)));poly(p)
 nativein=out/'native-input.txt';nativein.write_text('\n'.join(text)+'\n');flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'dual-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20)
 subprocess.run([str(exe),str(nativein),str(out/'native.json')],check=True,timeout=max(1,55-(time.monotonic()-started)))
 nv=read(out/'native.json');assert nv['all16_actual_rawH_annihilations_exact']and nv['v3_is_one']
 cert={'source_binding':'All112 rawH entries from accepted actualTn preflight; originalq6 is the exact compact/source319 matrix also used by jointLell; zeroHquadratic terms checkednative.',
       'one_curve':'P100','dual_recipe':'Theta row8/B3knownunit inverse, checkpointrow1 lift, CdeltaGL2 lift and constantGL7; normalizephysicalv3unit.',
       'five_NF117_Hquad_directions':'Their span is originalq0..q4: fourprime spanq1..q4 via accepted jointsourcechange; additionalH pivot isq0H. Allfive are evaluated explicitly before23dualisused.',
       'actual_native_results':nv,'all_source_counts':counts,
       'source_H1_cache_row118':'Actual135genuinequartic55 afteracceptedH2unit sourcechange andfullNF117, no pure98query; all21quarticsretained inparentcache.',
       'no_full134_Gaussian_POWERexpansion_point_prime_GB_or_rank_probe':True,
       'no_full57_or_full239_claim':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'one_curve':'P100','actual_rawH_dual_all16_verified':True,
         'all6_original_quadric_evaluations_verified':True,'NF117_H5_descent':nv['all5_eliminated_Hquad_directions_annihilated'],
         'quadratic_zero_flags':nv['quadratic_zero_flags'],'optional_branch':nv['optional_branch'],
         'optional_scalar_nonzero':nv['optional_scalar_nonzero'],'optional_scalar_unit_entire_curve':nv['optional_scalar_unit_on_entire_curve'],
         'native_peak_poly_length':nv['peak_polynomial_length'],'stored_output_terms':nv['stored_output_terms'],
         'no_full57_full239_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
         'native_input_sha256':sha(nativein),'native_output_sha256':sha(out/'native.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
