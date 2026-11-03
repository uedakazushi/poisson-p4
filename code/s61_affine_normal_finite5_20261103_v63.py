#!/usr/bin/env python3
"""Full S61 algebra finite on five generators over its42normal coordinates.

Native exact inverse produces the finite quadratic cross-reduction recipe.
The derived107MB expansion stays in ignoredwork and is identified by its hash.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_affine_fiber_macaulay import load_fiber
ROOT=Path(__file__).resolve().parents[1];PRIME=32003

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(poly,mon,c):
 v=poly.get(mon,Fraction(0))+c
 if v:poly[mon]=v
 else:poly.pop(mon,None)
def mod(c):return int(s.numer(c)%PRIME)*pow(int(s.denom(c)%PRIME),-1,PRIME)%PRIME
def rows_mod(M):
 A=[[mod(c)for c in row]for row in M.tolist()];ids=list(range(M.rows));rank=0
 for j in range(M.cols):
  row=next((i for i in range(rank,M.rows)if A[i][j]),None)
  if row is None:continue
  A[rank],A[row]=A[row],A[rank];ids[rank],ids[row]=ids[row],ids[rank];inv=pow(A[rank][j],-1,PRIME);A[rank]=[c*inv%PRIME for c in A[rank]]
  for i in range(rank+1,M.rows):
   c=A[i][j]
   if c:A[i]=[(a-c*b)%PRIME for a,b in zip(A[i],A[rank])]
  rank+=1
 return ids[:rank]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();out=args.output.resolve();work=args.work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 dep=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/report.json';old=json.loads(dep.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(dep.parent/rel)==want,rel
 source=dep.parent/'linear129x42-QQ.json.gz';data=json.loads(gzip.decompress(source.read_bytes()));A=[[{tuple(m):Fraction(c)for m,c in poly}for poly in row]for row in data['all129_linear_factor_rows']]
 normal_monomials=[(j,)for j in range(42)]+list(itertools.combinations_with_replacement(range(42),2));assert len(normal_monomials)==945
 cross=[[A[i][j].get((k,),Fraction(0))for k in range(3)for j in range(42)]for i in range(129)];rest=[]
 for row in A:
  poly={}
  for j,column in enumerate(row):
   for mon,c in column.items():
    if mon and mon[0]<3:continue
    add(poly,tuple(sorted((j,)+tuple(k-3 for k in mon))),c)
  rest.append([poly.get(m,Fraction(0))for m in normal_monomials])
 assert all(c.denominator%PRIME for row in cross+rest for c in row)
 raw=('129 126 945\n'+'\n'.join(' '.join(map(str,row))for row in cross+rest)+'\n').encode();inputpath=work/'cross-input.txt';inputpath.write_bytes(raw);(out/'cross-input.txt.gz').write_bytes(gzip.compress(raw,mtime=0))
 cpp=ROOT/'code/s61_affine_normal_cross_inverse_20261103_v63.cpp';native=work/'native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();include=subprocess.check_output(['pkg-config','--variable=includedir','flint'],text=True).strip()
 subprocess.run(['g++','-O2','-std=c++17','-I'+include,str(cpp),'-o',str(native),*flags],check=True)
 derived=work/'cross-output.txt';native_started=time.monotonic()
 result=subprocess.run([str(native),str(inputpath),str(derived)],capture_output=True,text=True,timeout=60,check=True,preexec_fn=lambda:resource.setrlimit(resource.RLIMIT_AS,(2200*1024*1024,2200*1024*1024)))
 assert 'all126x945 mixed-coefficient identities PASS'in result.stdout
 (out/'native.log').write_text(result.stdout+result.stderr)
 with derived.open()as handle:
  assert handle.readline().split()==['126','945'];selected=list(map(int,handle.readline().split()))
 assert len(selected)==len(set(selected))==126
 # Independently reconstruct the actual five plane restrictions and their
 # homogeneous quadratic/cubic spanning blocks. Mixed reduction removes onlyn*z.
 F,qids,center,original,pointpath=load_fiber();embedding=[[Fraction(c)for c in row]for row in data['embedding']];plane=[]
 for polynomial in F[:5]:
  outpoly={}
  for mon,c in polynomial.items():
   if not mon:add(outpoly,(),c)
   elif len(mon)==1:
    for k in range(3):add(outpoly,(k,),c*embedding[mon[0]][k])
   else:
    for k in range(3):
     for l in range(3):add(outpoly,tuple(sorted((k,l))),c*embedding[mon[0]][k]*embedding[mon[1]][l])
  plane.append(outpoly)
 quads=list(itertools.combinations_with_replacement(range(3),2));cubics=list(itertools.combinations_with_replacement(range(3),3));Q=s.Matrix([[s.Rational(poly.get(m,0))for m in quads]for poly in plane]);pivots=rows_mod(Q.T);assert len(pivots)==5
 qi=DomainMatrix.from_Matrix(Q[:,pivots]).inv().to_Matrix();assert Q[:,pivots]*qi==s.eye(5);assert all(s.denom(c)%PRIME for c in qi)
 # Multiplication byz_k simply removes one occurrence; no multiplicityfactor.
 cubic=s.Matrix([[s.Rational(poly.get(tuple(list(mon[:mon.index(k)])+list(mon[mon.index(k)+1:])),0))if k in mon else 0 for mon in cubics]for k in range(3)for poly in plane])
 cubicrows=rows_mod(cubic);assert len(cubicrows)==10;ci=DomainMatrix.from_Matrix(cubic[cubicrows,:]).inv().to_Matrix();assert cubic[cubicrows,:]*ci==s.eye(10);assert all(s.denom(c)%PRIME for c in ci)
 basis=[[],[0],[1],[2],list(next(m for j,m in enumerate(quads)if j not in pivots))]
 cert=dict(normal_base_variables=42,plane_variables=3,cross_matrix_shape=[129,126],exact_cross_rank=126,selected_cross_rows=selected,normal_monomials=normal_monomials,
  all126_QQ_cross_identities='z_k*n_j=P_kj(n), degreeP<=2;P=-(C_selected)^-1*D_selected',
  inverse_and_all126x945_coefficients_verified_QQ=True,all_selected_inverse_denominators_prime_units=True,prime=PRIME,
  polynomial_cofactor_recipe='P coefficientmatrix is nativeoutput after126inverse rows; missing3G are D_other+C_other*P, pure-normalquadrics.',
  plane_quadratic_monomials=quads,plane_quadratic_matrix=[[str(c)for c in row]for row in Q.tolist()],plane_quadratic_pivots=pivots,plane_quadratic_inverse=[[str(c)for c in row]for row in qi.tolist()],
  plane_cubic_matrix=[[str(c)for c in row]for row in cubic.tolist()],plane_cubic_rows=cubicrows,plane_cubic_inverse=[[str(c)for c in row]for row in ci.tolist()],
  fixed_module_generators_z_monomials=basis,full_actual_algebra_finite_over_normal_base=True,module_generators_atmost5=True,module_freeness_not_claimed=True,
  residual_normal_projection_closed_immersion=True,
  exact_scalar_normal_ideal_recipe=dict(pure_normal_quadrics=3,compatibility_cubics=2583,cleared_plane_quartics=210,total_generators=2796,formula='J=(H3,n_i*P_kj-n_j*P_ki, n_i^2*fbar_l(P_0i/n_i,P_1i/n_i,P_2i/n_i,n)); I_res_normal=J:(n_0,...,n_41)^infinity'),
  actual_full_beta_finiteness_or_residual_emptiness_not_claimed=True,
  native_output_sha256=sha(derived),native_output_bytes=derived.stat().st_size,native_input_sha256=hashlib.sha256(raw).hexdigest(),native_seconds=time.monotonic()-native_started,
  large_expansion_policy='Reproducible exactinverse/126x945 derivedcoefficients stayignoredwork; source,input,outputhashandexactverificationlogareversioned.')
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),cpp,source,dep,original,pointpath,ROOT/'research/notes/S61_AFFINE_NORMAL_FINITE5_20261103_v63.md']
 report=dict(id='S61_AFFINE_NORMAL_FINITE5_20261103_v63',all_passed=True,new_actual_finite_module_over_42normal_base_generators5=True,exact_mixed_cross_rank126=True,p_integral_cross_and_degree23_reductions_at_frozen_prime=True,residual_closed_embedding_into_normal_affine42=True,low_degree_normal_saturation_criterion_retains_full_residual=True,finite_over_beta_or_unit_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('ACTUAL FULL S61 FINITE5 OVER NORMAL42; RESIDUAL CLOSED NORMAL EMBEDDING PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
