#!/usr/bin/env python3
"""Stream10 actual pure cubic kernels using saved504 pivots; no quartic job."""
import argparse,gzip,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 pf=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_PREFLIGHT_20261104_v63';pr=read(pf/'report.json');pc=read(pf/'certificate.json');assert sha(pf/'certificate.json')==pr['certificate_sha256']
 matrix=base/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/remaining.txt.gz';bp=matrix.with_name('report.json');sp=base/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';adp=base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json'
 rawp=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json';sigp=base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json';etap=base/'S61_NONSQUARE_P100_H1_DUAL_20261104_v63/native.json'
 cip=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/cubic-coefficient-input.json';pivp=base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/cubic504.json';litp=base/'S61_NONSQUARE_P100_PURE21_SOURCE_PREFLIGHT_20261104_v63/source-input.json';g18p=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/physical-cubic84.json'
 cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp'
 paths=[Path(__file__),cpp,utility,pf/'report.json',pf/'certificate.json',matrix,bp,sp,adp,rawp,sigp,etap,cip,pivp,litp,g18p]
 for p in paths[5:]:assert sha(p)==pr['input_sha256'][str(p.relative_to(ROOT))]
 d=read(rawp);raw=d['raw16x7'];del d;d=read(sigp);sig=d['linear_lift_coefficients7x16'];assert d['all7_fixed_source_linear_lifts_verified'];del d
 d=read(etap);eta=d['eta_physical7'];del d;d=read(cip);oldmons=d['old_physical_monomials84'];pivotmap=d['sources6'];del d
 piv=read(pivp)['values6x84'];g18=read(g18p)['values84'];lit=read(litp);g18=[g18[lit['physical_pure_cubic_monomials84'].index(m)]for m in oldmons];del lit
 b,s,ad=map(read,[bp,sp,adp]);mons=s['normal_monomials'];monindex={tuple(m):j for j,m in enumerate(mons)};purecols=[monindex[tuple(m)]for m in oldmons];hcols=[]
 for v in[0,1,2,3,4,6,7]:
  m=[0]*8;m[5]=1;m[v]+=1;hcols.append(monindex[tuple(m)])
 P=[1,6,8,9,11,14];I=[i for i in range(16)if i not in P];assert [r['genuine_index']for r in pivotmap]==P
 for p,r in zip(P,pivotmap):assert ad['basis57_source_definitions'][41+p]['original_equation']==r['actual135_row']
 deg=lambda p:max((m[0]for m,c in p),default=-1)
 def rfwidth(v):return max(0,deg(v['numerator'])+1)+deg(v['denominator'])+1
 def cb(p):
  if not p:return(-1,0)
  clear=max(2*m[0]+m[1]for m,c in p);return clear+max(2*m[2]-m[0]for m,c in p),clear
 def addbs(v):
  v=[b for b in v if b[0]>=0]
  if not v:return(-1,0)
  dd=sum(b[1]for b in v);return max(b[0]+dd-b[1]for b in v),dd
 def ww(v):return 1 if v[0]<0 else v[0]+v[1]+2
 def rfb(v):return(-1,0)if not v['numerator']else(deg(v['numerator']),deg(v['denominator']))
 def alpha_bounds(i):return[addbs([(cb(raw[i][j])[0]+rfb(sig[j][p])[0],cb(raw[i][j])[1]+rfb(sig[j][p])[1])for j in range(7)if raw[i][j]and sig[j][p]['numerator']])for p in P]
 coflive=0
 for i in I:
  ab=alpha_bounds(i);old=sum(ww(v)for v in ab)+10;now=old+4+ww(addbs([ab[0],ab[1]]));coflive=max(coflive,sum(ww(v)for v in ab)+old+now)
 body=out/'native-body.tmp';rawterms=sum(len(p)for row in raw for p in row);maxcore=0;current84bound=0;total_new_core=0;source_records=[]
 literalrf=sum(len(v['numerator'])+len(v['denominator'])for row in piv for v in row)+sum(len(v['numerator'])+len(v['denominator'])for row in sig for v in row)+sum(len(v['numerator'])+len(v['denominator'])for v in eta+g18)
 canonical=lambda p:sorted((tuple(m),str(Q(c)))for m,c in p if Q(c))
 with body.open('w')as f:
  def poly(p):f.write(str(len(p))+'\n');f.writelines(f'{c} {m[0]} {m[1]} {m[2]}\n'for m,c in p)
  def rf(v):
   for k in['numerator','denominator']:f.write(str(len(v[k]))+'\n');f.writelines(f'{m[0]} {c}\n'for m,c in v[k])
  for row in raw:
   for p in row:poly(p)
  for row in sig:
   for v in row:rf(v)
  for v in eta:rf(v)
  for row in piv:
   for v in row:rf(v)
  f.write('10\n');hh=hashlib.sha256()
  with gzip.open(matrix,'rb')as zf:
   header=zf.readline();hh.update(header);assert header.strip()==b'57 91 57'
   for pos in range(57):
    i=pos-41;store=i in I and i!=0;rr={}
    for col in range(91):
     line=zf.readline();assert line;hh.update(line)
     if 0<=i<16 and(col in hcols or(store and col in purecols)):
      t=line.decode().split();n=int(t[0]);assert len(t)==1+4*n;p=[[list(map(int,t[2+4*k:5+4*k])),str(Q(t[1+4*k]))]for k in range(n)]
      if col in hcols:assert canonical(p)==canonical(raw[i][hcols.index(col)])
      if store and col in purecols:rr[col]=p;maxcore=max(maxcore,n)
    if i in I:
     f.write(f'{i} {int(i==0)}\n')
     if i==0:
      for v in g18:rf(v)
      current84bound=max(current84bound,sum(rfwidth(v)for v in g18))
     else:
      assert len(rr)==84
      for c in purecols:poly(rr[c])
      total_new_core+=sum(len(rr[c])for c in purecols);current84bound=max(current84bound,sum(ww(cb(rr[c]))for c in purecols))
     source_records.append({'genuine_index':i,'basis57_position':pos,'actual135':ad['basis57_source_definitions'][pos]['original_equation'],'original63':b['remaining_original_row_indices'][pos],'mode':'cached_G18'if i==0 else'streamed_complete_normalized_core_polynomials','row_unit_division_already_in_saved_source':ad['basis57_encoded_source_units'][pos]})
    rr.clear()
   assert not zf.read()
  assert hh.hexdigest()==b['remaining_matrix_sha256']
 # Retain no large row/core dictionary or full-source JSON while native runs.
 del rr
 if't'in locals():del t,p
 reserve=literalrf+rawterms+maxcore+256
 shared=sum(rfwidth(v)for row in piv for v in row)+sum(rfwidth(v)for row in sig for v in row)+sum(rfwidth(v)for v in eta)+sum(ww(cb(p))for row in raw for p in row)
 kernelmax=pc['bounds']['max_single_kernel84_output_terms_upper']+84
 live=reserve+shared+current84bound+kernelmax+coflive+2000+15000
 serialized_upper=pc['bounds']['all10_retained_output_terms_upper']+840+10*coflive
 assert live<150000 and serialized_upper<150000
 pre={'stream_literal_reserve':reserve,'cached_shared_dense_RF_upper':shared,'max_current_nonpivot84_RF_upper':current84bound,'max_kernel84_RF_upper':kernelmax,'current_source_cofactor_RF_upper':coflive,'power_cache_upper':2000,'temporary_upper':15000,'total_live_upper':live,'serialized_upper':serialized_upper,'max_streamed_core_coefficient_terms':maxcore,'newly_specialized_nine_nonpivot_core_terms':total_new_core,'whole_139746_source_JSON_never_loaded':True,'six_pivot504_values_reused_no_revaluation':True,'G18_84_values_reused_no_revaluation':True}
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');print(json.dumps(pre,indent=2),flush=True)
 ip=out/'native-input.txt'
 with ip.open('wb')as f:
  f.write((str(reserve)+'\n').encode())
  with body.open('rb')as bf:
   for chunk in iter(lambda:bf.read(1024*1024),b''):f.write(chunk)
 body.unlink();assert ip.stat().st_size<8*1024*1024
 bindings={'nonpivot_sources10':source_records,'six_pivot_bindings':pivotmap,'six_pivot_source_proof':'Accepted FrameAdapter full57 identity and independently accepted504 curve cache; positions/actual135/ROWUNIT DIVISION/mon84 checked. No repeated504 evaluation.','pure_monomials84':oldmons,'full112_Hcoefficients_match_complete_normalized57_matrix':True,'original_sigma_frame':True,'actual_delta7_cofactor_transport':['new6=old6/delta','new1=old1+old6'],'source_basis_determinant_original':1,'source_basis_determinant_actual_delta7':'1/delta, known whole-curve unit','no70products_NF_rank_or_new_oracle':True}
 (out/'source-binding.json').write_text(json.dumps(bindings,indent=2)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'kernel-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:run=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 audit=read(out/'audit.json')if(out/'audit.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':run.returncode}
 complete=run.returncode==0 and audit['all_passed']
 if complete:assert audit['pure_coefficient_identity_checks']==1680 and audit['peak_live_RF_terms_with_stream_input']<=150000
 report={'id':ID,'complete':complete,'ten_whole_curve_actual_pure_cubic_source_kernels':complete,'formal_kernel_basis_size10_NOT_image_rank':True,'native_audit':audit,'no70quartic_NF_rank_GB_newpoint_or_new_oracle':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
