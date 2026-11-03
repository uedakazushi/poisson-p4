#!/usr/bin/env python3
"""Independent compact audit of transported scalar witness and dual-column change."""
import argparse,gzip,json,resource,shlex,subprocess,time
from pathlib import Path
import sympy as sp
from s61_nonsquare_scalar_product_rows_review_20261103_v63 import ROOT,R,a,e,P,sha,decode,parse,dot
NAME='S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62'
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();folder=ROOT/'reports'/NAME
 with gzip.open(folder/'certificate.json.gz','rt')as f:d=json.load(f)
 first=ROOT/'reports/S61_BINARY_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v62/certificate.json.gz'
 with gzip.open(first,'rt')as f:old=json.load(f)
 parent=json.loads((folder/'report.json').read_text());assert parent['all_passed']
 source=ROOT/d['source'];assert sha(source)==d['source_sha256']==old['native_outputs']['row7_replacement']['sha256'];assert d['target_row']==6
 row=list(map(decode,d['row_sparse']));cofs=list(map(decode,d['cofactor_sparse']));assert len(row)==len(cofs)==59
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[25,59,1];assert f.readline().split()==['43'];ids=list(map(int,f.readline().split()));assert ids[6]==0
  for i in range(25):
   for j in range(59):
    s=f.readline();assert s
    if i==6:assert parse(s)==row[j]
  assert not f.read().strip()
 A0=list(map(decode,old['rows_sparse']['0']));A7=list(map(decode,old['rows_sparse']['7']));c7=list(map(decode,old['cofactors_sparse']['7']));T=[[decode(p)for p in line]for line in old['two_row_right_inverse_sparse']]
 assert c7[43].is_ground and c7[43];ratio=T[43][0].exquo(c7[43]);transport=[T[j][0]-ratio*c7[j]for j in range(60)if j!=43];assert transport==cofs and sum(map(len,cofs))==138
 cross=dot(A0,c7);raw=[a*A0[j]-cross*A7[j]for j in range(60)if j!=43];scale=R(sp.QQ(parent['row_cleaning_unit_monomial']['coefficient']));assert all(p==scale*a*q for p,q in zip(raw,row))
 target=a.exquo(scale);assert dot(row,cofs)==target
 alpha={int(j):sp.QQ(v)for j,v in d['constant_combo'].items()};assert set(alpha)=={50,53,54,55,57}and alpha[50]
 assert sum((cofs[j]*c for j,c in alpha.items()),R.zero)==R.one
 # Full literal target-row check for the dual QQ column/cofactor change.
 transformed=[]
 for j in range(59):
  if j==50:transformed.append(row[j].exquo(R(alpha[50])))
  else:transformed.append(row[j]-row[50]*R(alpha.get(j,sp.QQ.zero)/alpha[50]))
 newcofs=cofs[:];newcofs[50]=R.one;assert dot(transformed,newcofs)==target
 # Fresh DIFFERENT nonsquare2x3 fixture for the newly extended native API.
 # q=(e,e+1,0), alpha=(-1,1,0), source rows(a,-a,0),(1,e,P).
 # Its completed residual must be (e+1,P), keeping the unknown e+1 factor.
 fixture=out/'fixture';fixture.mkdir();(fixture/'input.txt').write_text('2 3 0\n1 1 1 0 0\n1 -1 1 0 0\n0\n1 1 0 0 0\n1 1 0 1 0\n1 1 0 0 1\n');(fixture/'witness.txt').write_text('3\n1 1 0 1 0\n2 1 0 1 0 1 0 0 0\n0\n');(fixture/'alpha.txt').write_text('3\n-1\n1\n0\n');(fixture/'units.txt').write_text('a\n')
 cpp=ROOT/'code/polynomial_cofactor_constant_span_pivot_20261103_v62.cpp';flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout);exe=out/'fixture-helper'
 subprocess.run(['c++','-std=c++17','-O2',str(cpp),*flags,'-lgmp','-o',str(exe)],check=True)
 def limits():resource.setrlimit(resource.RLIMIT_AS,(1000000000,1000000000));resource.setrlimit(resource.RLIMIT_CPU,(10,11))
 with(fixture/'log.txt').open('w')as f:p=subprocess.run([str(exe),str(fixture/'input.txt'),str(fixture/'output.txt'),str(fixture/'units.txt'),'0',str(fixture/'witness.txt'),'0',str(fixture/'alpha.txt')],stdout=f,stderr=subprocess.STDOUT,timeout=12,preexec_fn=limits)
 assert p.returncode==0;log=(fixture/'log.txt').read_text();assert'EXACT_COFACTOR_CONSTANT_COMBINATION 2 EXPOSED_COLUMN 0'in log and'EXACT_WITNESS_REPLACE 0 TARGET 0'in log and'PASS 1'in log
 with(fixture/'output.txt').open()as f:
  assert list(map(int,f.readline().split()))==[1,2,1];assert f.readline().split()==['0'];assert f.readline().split()==['1'];assert parse(f.readline())==e+1 and parse(f.readline())==P
 exe.unlink()
 output=ROOT/d['output']['path'];assert sha(output)==d['output']['sha256']==parent['generated_output']['sha256']and output.stat().st_size==d['output']['bytes'];assert d['output']['shape']==[24,58]
 native=(folder/'native.log').read_text();assert'EXACT_COFACTOR_CONSTANT_COMBINATION 5 EXPOSED_COLUMN 50'in native and'EXACT_WITNESS_REPLACE 50 TARGET 6'in native and'TERMS 2755382'in native and'PASS 1'in native
 for path,wanted in parent['input_sha256'].items():assert sha(ROOT/path)==wanted
 inputs=[Path(__file__),ROOT/'code/s61_nonsquare_scalar_product_rows_review_20261103_v63.py',folder/'report.json',folder/'certificate.json.gz',first,cpp,folder/'native.log',ROOT/'research/notes'/f'{NAME}.md']
 report={'id':'S61_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_V63_REVIEW','all_passed':True,'independent_algebra':'SymPy QQ[a,e,P], not producer Fraction code','all59_actual_source_row_entries_compared':True,'metadata_target_row6_is_original_row0':True,'138term_transport_independently_rebuilt':True,'actual_row_cleaning_unit_independently_checked_all59':'(422053134351/5039050240)*a','constant_combo5_components_equals1':True,'literal_dual_column_change_and_target_identity_checked':True,'new_native_API_fresh_2x3_fixture_pass':True,'fixture_unknown_eplus1_factor_retained':True,'full_growing_native_output_not_recomputed':True,'generated_output_sha256_checked':d['output']['sha256'],'accepted_residual_shape':[24,58],'accepted_constructive_product_unit_pivots':33,'allpoint_product41_fullcubic57_or_ambient_claim':False,'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','138term_transport_independently_rebuilt','new_native_API_fresh_2x3_fixture_pass','accepted_residual_shape','seconds']},indent=2))
if __name__=='__main__':main()
