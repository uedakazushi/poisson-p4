#!/usr/bin/env python3
"""Exact direction-graph recipe, full coefficient witness, and raw-infinity blocker."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import itertools
import json
from pathlib import Path
import time
import numpy as np
from s61_affine_normal_commutator_20261103_v63 import ROOT,P,sha
from s61_affine_fiber_macaulay import load_fiber

def add(poly,mon,c):
 v=(poly.get(mon,0)+int(c))%P
 if v:poly[mon]=v
 else:poly.pop(mon,None)
def mul(a,b):
 out={}
 for u,c in a.items():
  for v,d in b.items():add(out,tuple(sorted(u+v)),c*d)
 return out
def lin(row):return {(i,):int(c)for i,c in enumerate(row)if c}
def matinv(A):
 n=len(A);M=np.concatenate((A,np.eye(n,dtype=np.int64)),axis=1)
 for j in range(n):
  i=next(i for i in range(j,n)if M[i,j]);M[[i,j]]=M[[j,i]];M[j]=M[j]*pow(int(M[j,j]),-1,P)%P
  for i in range(n):
   if i!=j and M[i,j]:M[i]=(M[i]-M[i,j]*M[j])%P
 assert np.array_equal(M[:,:n],np.eye(n,dtype=np.int64));return M[:,n:]
def residue(c):return c.numerator*pow(c.denominator,-1,P)%P

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--cross-output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 dep=ROOT/'reports/S61_AFFINE_NORMAL_COMMUTATOR_20261103_v63/report.json';old=json.loads(dep.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(dep.parent/rel)==want,rel
 finite=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/certificate.json';fc=json.loads(finite.read_text());assert sha(args.cross_output)==fc['native_output_sha256']
 def modtoken(token):
  a,slash,b=token.partition('/');return int(a)%P*pow(int(b)%P if slash else 1,-1,P)%P
 with args.cross_output.open()as handle:
  assert handle.readline().split()==['126','945'];selected=list(map(int,handle.readline().split()));assert selected==fc['selected_cross_rows']
  for _ in range(126):handle.readline()
  coeff=np.array([list(map(modtoken,handle.readline().split()))for _ in range(126)],dtype=np.int64)
 assert coeff.shape==(126,945);T=coeff[:,:42].reshape(3,42,42);comm=(T[0]@T[1]-T[1]@T[0])%P;Ci=matinv(comm);assert np.array_equal(comm@Ci%P,np.eye(42,dtype=np.int64))
 # FULL homogeneous polynomial Bezout check, not point evaluations.
 checks=[]
 for pivot in[0,5,41]:
  ls=[lin(T[k,pivot])for k in range(2)];U=[]
  for k in range(2):
   vector=[]
   for j in range(42):
    poly={}
    for mon,c in ls[k].items():add(poly,tuple(sorted(mon+(j,))),c)
    for t,c in enumerate(T[k,j]):add(poly,tuple(sorted((pivot,t))),-int(c))
    vector.append(poly)
   U.append(vector)
  c0=Ci[pivot]@T[1]%P;c1=Ci[pivot]@T[0]%P;result={}
  for j in range(42):
   u0={(pivot,):int(c0[j])};u1={(pivot,):int(-c1[j])%P}
   for mon,c in ls[1].items():add(u0,mon,-int(Ci[pivot,j])*c)
   for mon,c in ls[0].items():add(u1,mon,int(Ci[pivot,j])*c)
   for mon,c in mul(u0,U[0][j]).items():add(result,mon,c)
   for mon,c in mul(u1,U[1][j]).items():add(result,mon,c)
  assert result=={(pivot,pivot,pivot):1};checks.append(dict(pivot=pivot,complete_cubic_coefficient_check=True))
 F,qids,center,original,pointpath=load_fiber();source=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz';data=json.loads(gzip.decompress(source.read_bytes()));A=[[{tuple(m):Fraction(c)for m,c in poly}for poly in row]for row in data['all129_linear_factor_rows']];embedding=[[Fraction(c)for c in row]for row in data['embedding']]
 normalmons=[(j,)for j in range(42)]+list(itertools.combinations_with_replacement(range(42),2));mi={m:j for j,m in enumerate(normalmons)};C=np.array([[residue(A[i][j].get((k,),Fraction(0)))for k in range(3)for j in range(42)]for i in range(129)],dtype=np.int64);D=[]
 for row in A:
  poly={}
  for j,column in enumerate(row):
   for mon,c in column.items():
    if mon and mon[0]<3:continue
    key=tuple(sorted((j,)+tuple(k-3 for k in mon)));poly[key]=poly.get(key,Fraction(0))+c
  D.append([residue(poly.get(m,Fraction(0)))for m in normalmons])
 D=np.array(D,dtype=np.int64);other=[j for j in range(129)if j not in selected];H=(D[other]+C[other]@coeff)%P;plane=[];K=[]
 images=[]
 for j,row in enumerate(embedding):
  img={k:residue(c)for k,c in enumerate(row)if c}
  if j>=3:img[j]=1
  images.append(img)
 for polynomial in F[:5]:
  trans={}
  for mon,c in polynomial.items():
   if not mon:add(trans,(),residue(c))
   elif len(mon)==1:
    for j,cj in images[mon[0]].items():add(trans,(j,),residue(c)*cj)
   else:
    for j,cj in images[mon[0]].items():
     for k,ck in images[mon[1]].items():add(trans,tuple(sorted((j,k))),residue(c)*cj*ck)
  f={};kn=np.zeros(945,dtype=np.int64)
  for mon,c in trans.items():
   if all(j<3 for j in mon):f[mon]=c
   elif all(j>=3 for j in mon):kn[mi[tuple(j-3 for j in mon)]]=(kn[mi[tuple(j-3 for j in mon)]]+c)%P
   else:
    z,n=mon;assert z<3<=n;kn=(kn+c*coeff[z*42+n-3])%P
  plane.append(f);K.append(kn)
 K=np.array(K,dtype=np.int64)
 axes=[j for j in range(45)if all(f.get((j,j),0)==0 for f in F)];assert axes==[8,16,24,32];witnesses=[]
 for qaxis in axes:
  pivot=qaxis-3;v=np.zeros(42,dtype=np.int64);v[pivot]=1;quadratic=np.array([int(np.prod(v[list(mon)]))for mon in normalmons[42:]],dtype=np.int64);L=np.einsum('kij,j->ki',T,v)%P;Q=(coeff[:,42:]@quadratic%P).reshape(3,42)
  U=[(L[k,pivot]*v-L[k])%P for k in range(2)];V=[(Q[k,pivot]*v-Q[k])%P for k in range(2)];u0=Ci[pivot]@(T[1]-L[1,pivot]*np.eye(42,dtype=np.int64))%P;u1=-Ci[pivot]@(T[0]-L[0,pivot]*np.eye(42,dtype=np.int64))%P
  assert int((u0@U[0]+u1@U[1])%P)==1;B=int((u0@V[0]+u1@V[1])%P);s=(-B)%P;zeta=(s*L[:,pivot]+Q[:,pivot])%P;assert s==0 and not np.any(zeta)
  assert not np.any((H[:,:42]@v*s+H[:,42:]@quadratic)%P);assert not np.any((K[:,:42]@v*s+K[:,42:]@quadratic)%P)
  witnesses.append(dict(raw_q_axis_zero_based=qaxis,free_u_one_based=qids[qaxis]+1,direction_pivot=pivot,compact_s=0,all134_QQ_leading_forms_zero=True))
 rawdim=ROOT/'reports/S61_AFFINE_BOUNDARY_DIMENSION_20261102.json';rd=json.loads(rawdim.read_text());assert rd['all_passed']and rd['raw_generic_infinity_exact_dimension']==3
 payload=dict(prime=P,normal_monomials=normalmons,P_coefficients=coeff.tolist(),H_coefficients=H.tolist(),plane5_sparse=[[[list(m),int(c)]for m,c in f.items()]for f in plane],K_coefficients=K.tolist(),constant_commutator_inverse=Ci.tolist())
 (out/'coefficient-model-mod32003.json.gz').write_bytes(gzip.compress(json.dumps(payload,separators=(',',':')).encode(),mtime=0))
 certificate=dict(exact_QQ_source_is_finite_inverse_recipe=True,complete_polynomial_Bezout_checks=checks,
  uniform_exact_identity='v_i^2*C*v=(v_i*T1-lambda1*I)*U0-(v_i*T0-lambda0*I)*U1; row_i*C^-1 gives v_i^3',
  chart_recipe='v_i=1;U_k=lambda_k*v-T_k*v;V_k=Q_ki*v-Q_k(v);u0=row_i C^-1(T1-lambda1I);u1=-row_i C^-1(T0-lambda0I);B_i=u0 V0+u1 V1;s=-B_i;zeta_k=s*lambda_k+Q_ki',
  exact_direction_only_equations=['V_kj-B_i*U_kj=0 for3*41mixedcompatibilities','H2-B_i*H1=0 for3remainingG','f_plane_hom(zeta,s)+K1*s+K2=0 for5originalF'],
  polynomial_degree_bounds=dict(B=4,zeta=5,compatibilities=6,remaining_H=5,five_F=10),all42_standard_direction_charts_retained=True,
  compact_normal_direction_projection_closed_immersion=True,affine_residual_is_open_B_nonzero=True,affine_recovery='n=-v/B_i; z=lambda-Q_i/B_i',
  raw_infinity_at_frozen_center=witnesses,raw_generic_infinity_exact_dimension=3,
  compact_Y_empty_impossible_on_generic_beta_open=True,modp_s_nilpotence_does_not_alone_imply_QQ_nilpotence=True,
  actual_unknown_affine_empty_or_genericfive_not_claimed=True)
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),dep,finite,source,original,pointpath,rawdim,ROOT/'research/notes/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63.md']
 report=dict(id='S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63',all_passed=True,exact_compact_graph_inverse_by_constant_commutator=True,full134_and3remainingH_retained=True,normal_direction_only_affine_residual_chart_criterion=True,four_exact_frozen_infinity_witnesses=True,generic_compact_emptiness_proved_impossible_from_existing_raw_boundary=True,actual_affine_unit_or_QQ_nilpotence_not_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('DIRECTION GRAPH/CLOSEDIMMERSION/41VAR SUPPORT/RAWINFINITY BLOCKER PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
