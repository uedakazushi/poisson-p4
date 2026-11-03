#!/usr/bin/env python3
"""Independent QQ identities for the two modular-guided Dixon J0 units.

The prime used by the search is never used as mathematical evidence here.
"""
import argparse,gzip,json,shutil,time
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
from s61_nonsquare_J0_degree1_identity_20261103_v63 import sha,poly,encoded,decoded,identity
ROOT=Path(__file__).resolve().parents[1]
R,a,e=ring('a,e',sp.QQ)
FACTORS=[a,500*a+90*e**2+3*e-24,300*a*e+54*e**3+9*e**2-24*e-4,100*a+18*e**2-9*e-2,3*e+2]
NAMES=['a','F','Bs','C','3e+2']
def validate(p):
 rows=list(map(decoded,p['row_polynomials']));cofs=list(map(decoded,p['cofactor_polynomials']));target=identity(rows,cofs)
 assert target==decoded(p['target_polynomial']) and target
 q=R.from_dict({ex[:2]:sp.QQ(c.numerator,c.denominator)for ex,c in target.items()});exponents=[]
 for factor in FACTORS:
  n=0
  while True:
   quotient,remainder=divmod(q,factor)
   if remainder:break
   q=quotient;n+=1
  exponents.append(n)
 assert q and q.is_ground
 assert any(len(q)==1 and(0,0,0)in q for q in cofs)
 result={'active_cofactors':sum(bool(q)for q in cofs),'cofactor_terms':sum(map(len,cofs)),'constant_nonzero_cofactor_exists':True,'unit_factor_exponents':dict(zip(NAMES,exponents)),'nonzero_QQ_scalar':str(q)}
 assert result==p['independent_identity_result'];return result
def payload(row,witness,shape,targetrow,source):
 with row.open()as f:
  assert list(map(int,f.readline().split()))==[1,shape[1],0];ss=[f.readline()for _ in range(shape[1])];assert not f.read().strip()
 # Verify that the saved search row is the actual complete pinned matrix row.
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[*shape,1];f.readline();f.readline();actual=[]
  for i in range(shape[0]):
   for j in range(shape[1]):
    s=f.readline()
    if i==targetrow:actual.append(s)
 assert ss==actual;rows=list(map(poly,ss))
 with witness.open()as f:
  assert int(f.readline())==shape[1];cofs=[poly(f.readline())for _ in range(shape[1])];assert not f.read().strip()
 target=identity(rows,cofs);p={'source_shape':shape,'source_sha256':sha(source),'target_row':targetrow,'row_polynomials':list(map(encoded,rows)),'cofactor_polynomials':list(map(encoded,cofs)),'target_polynomial':encoded(target)}
 # Produce the result by the same simple independent division, then replay it.
 q=R.from_dict({ex[:2]:sp.QQ(c.numerator,c.denominator)for ex,c in target.items()});powers=[]
 for factor in FACTORS:
  n=0
  while True:
   quotient,remainder=divmod(q,factor)
   if remainder:break
   q=quotient;n+=1
  powers.append(n)
 assert q and q.is_ground
 p['independent_identity_result']={'active_cofactors':sum(bool(q)for q in cofs),'cofactor_terms':sum(map(len,cofs)),'constant_nonzero_cofactor_exists':True,'unit_factor_exponents':dict(zip(NAMES,powers)),'nonzero_QQ_scalar':str(q)}
 validate(p);return p
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path);ap.add_argument('--replay',type=Path);args=ap.parse_args()
 if args.replay:
  with gzip.open(args.replay,'rt')as f:ps=json.load(f)
  print(json.dumps({'all_passed':True,'identities':[validate(p)for p in ps]},indent=2));return
 assert args.output;out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();base=ROOT/'work/nonsquare-J0-pid-v63'
 previous=ROOT/'reports/S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_v63/report.json';s0=base/'degree1-row4-v1/scalar-primitive/residual.txt';w0=base/'degree1-dixon-v1';s1=w0/'residual.txt';w1=base/'dixon-chain-v1/stage-1'
 assert sha(s0)==json.loads(previous.read_text())['normalized_residual_sha256'];assert sha(s1)=='aea93a891e98c3cd22639ad47f04d59d8b8bdcbbb6135355ea719d6be7f194d6'
 ps=[payload(base/'degree1-chain-v1/stage-1/row1.txt',w0/'row1.witness',[20,54],1,s0),payload(w1/'row0.txt',w1/'row0.witness',[19,53],0,s1)]
 with gzip.open(out/'identities.json.gz','wt',compresslevel=9)as f:json.dump(ps,f,separators=(',',':'))
 # Verify the ONLY successful new QQ normalization independently.  The later
 # 18x52 normalization timed out and its zero-byte output is never consumed.
 with(w0/'QQ-factors.txt').open()as f:
  assert list(map(int,f.readline().split()))==[19,53];rf=[];cf=[]
  for i in range(19):s=f.readline().split();assert s[:2]==['ROW',str(i)];rf.append(F(s[2]))
  for j in range(53):s=f.readline().split();assert s[:2]==['COLUMN',str(j)];cf.append(F(s[2]))
 assert all(rf+cf);terms=0;t=time.monotonic()
 with(w0/'raw-residual.txt').open()as f,s1.open()as g:
  for _ in range(3):assert f.readline()==g.readline()
  for i in range(19):
   for j in range(53):
    p=poly(f.readline());q=poly(g.readline());scale=rf[i]*cf[j]
    assert q=={ex:c*scale for ex,c in p.items()};assert all(c.denominator==1 for c in q.values());terms+=len(p)
  assert not f.read().strip()and not g.read().strip()
 assert'PASS 1'in(w0/'pivot.log').read_text()and'PASS 1'in(w1/'pivot.log').read_text()
 assert'EXACT_WITNESS_REPLACE'in(w0/'pivot.log').read_text()and'EXACT_WITNESS_REPLACE'in(w1/'pivot.log').read_text()
 raw=w1/'raw-residual.txt';assert sha(raw)=='2fb8b2e8636ecb2f68ebe7c58acad75c5cac6a55a63198f936b42d92ab4d1486'
 for folder,prefix,names in [(w0,'unit38',['row1.witness','row1.log','pivot.log','normalize.log','QQ-factors.txt','receipt.json','pivot-receipt.json','normalize-receipt.json']),(w1,'unit39',['row0.witness','row0.log','pivot.log','normalize.log','receipt.json'])]:
  for n in names:shutil.copyfile(folder/n,out/f'{prefix}-{n}')
 shutil.copyfile(base/'units.txt',out/'units.txt');shutil.copyfile(base/'dixon-chain-v1/report.json',out/'chain-receipt.json')
 inputs=[Path(__file__),ROOT/'code/s61_nonsquare_J0_degree1_identity_20261103_v63.py',ROOT/'code/s61_nonsquare_J0_degree1_dixon_20261103_v63.cpp',ROOT/'code/s61_nonsquare_J0_dixon_chain_20261103_v63.py',ROOT/'code/polynomial_one_witness_pivot_20261103_v62.cpp',previous,s0,s1,w0/'row1.witness',w1/'row0.witness']
 report={'id':'S61_NONSQUARE_J0_DIXON_UNITS_20261103_V63','all_passed':True,'base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]','J_not_inverted':True,'unknown_core_factor_inverted':False,'modular_data_used_only_to_select_QQ_solve':True,'independent_full_QQ_polynomial_identities':[p['independent_identity_result']for p in ps],'unit38_search_seconds':0.2158732540001438,'unit39_controller_seconds_to_witness':5.340495431999443,'total_safe_pivots_from57basis':39,'new_residual_shape':[18,52],'new_residual_sha256':sha(raw),'new_residual_bytes':raw.stat().st_size,'last_confirmed_state':str(raw.relative_to(ROOT)),'unit38_all_QQ_diagonal_terms_independently_verified':terms,'QQ_diagonal_check_seconds':time.monotonic()-t,'unit39_optional_normalization_failed_and_excluded':True,'chain_stop':'raw matrix exceeds1.2GB guard; optional QQ normalization timed out at62s; confirmed raw state retained','no_full_coker_unit_claim':True,'no_finite_J0_defect_claim':True,'large_outputs_omitted_as_reproducible':True,'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','total_safe_pivots_from57basis','new_residual_shape','independent_full_QQ_polynomial_identities','seconds']},indent=2))
if __name__=='__main__':main()
