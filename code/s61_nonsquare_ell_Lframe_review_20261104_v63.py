#!/usr/bin/env python3
"""Independent inverse-polynomial audit of frozen L-frame inputs; no rank job."""
import argparse,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def mons(n):return[tuple(c.count(i)for i in range(5))for c in itertools.combinations_with_replacement(range(5),n)]
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63';r=read(folder/'report.json');c=read(folder/'certificate.json');n=read(folder/'native.json')
 assert r['all_source_bound_coordinate_and_polynomial_input_checks_passed']and r['frozen']and r['no_new_rank_or_determinant_computation']
 for p,h in r['source_and_recipe_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in r['artifact_sha256'].items():assert sha(folder/p)==h,p
 circ=read(ROOT/'reports/S61_NONSQUARE_ELL_COLON_CUBIC11_20261104_v63/circuit.json')
 old=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63';raw=old/'frozen57-source.txt';assert sha(raw)==circ['source_uncompressed_sha256']
 source=read(ROOT/circ['normal_monomials_source']);mon8=source['normal_monomials'];assert len(mon8)==91
 assert c['new_five_lights']==['L','x0','y1','y2','y3'] and c['new_y1_y2_y3']==['old_x2','old_x3','old_x4']
 assert c['actual_genuine_rows']==list(range(41,57))and c['nine_monomial_pivot_columns']==list(range(9))
 keep=list(range(15,35))+list(range(9,15));assert c['remaining26_full35_column_indices']==keep
 assert c['full35_cubic_monomials_new_light_order']==[list(x)for x in mons(3)]
 assert c['pair_column_reversal']==list(reversed(range(6)))
 # Independently derive each original monomial restriction by repeated
 # integer convolution in (old_x0,...,old_x4,H), not selection.txt.
 z=(0,)*6;unit=[tuple(int(i==j)for i in range(6))for j in range(6)]
 images=[{unit[0]:1},{unit[1]:1},{unit[2]:1},{},{unit[3]:1},{unit[5]:1},{unit[1]:1,unit[0]:-3},{unit[4]:1}]
 cache={}
 def restrict(m):
  m=tuple(m)
  if m not in cache:
   p={z:1}
   for power,image in zip(m,images):
    for _ in range(power):
     q={}
     for e,a in p.items():
      for f,b in image.items():
       v=tuple(x+y for x,y in zip(e,f));q[v]=q.get(v,0)+a*b
     p={e:a for e,a in q.items()if a}
   cache[m]=p
  return cache[m]
 i3={m:i for i,m in enumerate(mons(3))};restriction=[]
 for m in mon8:
  row={}
  for v,a in restrict(m).items():
   if not v[5]:j=i3[v[:5]]
   else:assert v[5]==1 and sum(v[:5])==1;j=35+v[:5].index(1)
   row[j]=row.get(j,0)+a
  restriction.append(sorted(row.items()))
 with(out/'independent-meta.txt').open('w')as f:
  f.write('91\n')
  for row in restriction:f.write(str(len(row))+' '+ ' '.join(str(v)for pair in row for v in pair)+'\n')
  f.write(str(n['q5_C_factors_cancelled'])+' '+ ' '.join(map(str,n['genuine_C_factors_cancelled']))+'\n')
 qpath=ROOT/circ['quad_source'];qd=read(qpath);i2={m:i for i,m in enumerate(mons(2))};Q4=[[{}for _ in range(15)]for _ in range(4)]
 for row,qi in enumerate([0,1,3,5]):
  for m,terms in qd['quadratic_rows_sparse'][qi]:
   for v,b in restrict(m).items():
    assert not v[5];j=i2[v[:5]]
    for e,a in terms:
     e=tuple(e);Q4[row][j][e]=Q4[row][j].get(e,Q())+Q(a)*b
 with(out/'original-qfour4x15.txt').open('w')as f:
  f.write('4 15\n')
  for row in Q4:
   for p in row:
    terms=sorted((e,a)for e,a in p.items()if a);f.write(str(len(terms))+''.join(' '+str(a)+' '+' '.join(map(str,e))for e,a in terms)+'\n')
 cpp=ROOT/'code/s61_nonsquare_ell_Lframe_review_20261104_v63.cpp';exe=out/'review-native'
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 def remaining():return max(.1,60-(time.monotonic()-started))
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=remaining())
 pair=ROOT/c['pair_source'];cmd=[str(exe),str(out),str(raw),str(old/'quad20x35.txt'),str(pair),str(folder)]
 completed=subprocess.run(cmd,check=True,capture_output=True,text=True,timeout=remaining());exe.unlink()
 (out/'native.log').write_text(completed.stdout+completed.stderr);own=json.loads(completed.stdout);assert own['all_passed']
 assert own['quad10_terms']==n['new_quad10x26']['terms']==317 and own['genuine16x31_terms']==n['new_genuine16x31']['terms']==50085
 deps=[ROOT/'reports'/s/'report.json'for s in['S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW','S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in deps)
 note=ROOT/'research/notes'/(ID+'.md');inputs=[Path(__file__),cpp,note,folder/'report.json',folder/'certificate.json',folder/'native.json',ROOT/circ['normal_monomials_source'],qpath,raw,pair,*deps]
 result=dict(id=ID,all_passed=True,actual_qfour_and_genuine16_source_frame_bound=True,
  inverse_coordinate_replay_all60_quad_and640_genuine_coefficients=True,original_quad20x35_all700_entries_checked=True,
  actual_q1_q3_normalizations_and_L_derivative_checked=True,only_known_a_C_unit_normalizations=True,
  exactly9_monomial_unit_pivots_and_all26_plus5_heavy_columns_retained=True,
  actual_quad10x26_block_and_pair_all12_coefficients_checked=True,
  input_only_no_rank_border_identity_generic_upper10_or_full239_claim=True,
  no_new_point_prime_GB_field_solver_or_parent_rerun=True,native=own,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
