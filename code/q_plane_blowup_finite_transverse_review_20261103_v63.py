#!/usr/bin/env python3
"""Independent exact review of finite transverse actual-blowup certificates.

Direct sparse differentiation/polarization rebuilds all350 linear and all350
corrected quadratic rows without importing the producer or its unpack helper.
"""
import argparse,gzip,hashlib,json,subprocess,sys,time
from fractions import Fraction
from pathlib import Path
from math import prod
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads(gzip.decompress(path.read_bytes()))
def cv(v):
 v=Fraction(v);return s.QQ(v.numerator,v.denominator)
def dm(rows):return DomainMatrix([[cv(v) for v in row] for row in rows],(len(rows),len(rows[0])),s.QQ)
def saved(rows):return dm(rows)
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 producer=ROOT/'code/q_plane_blowup_finite_transverse_20261103_v63.py';reportpath=ROOT/'reports/Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_v63/report.json';certpath=reportpath.parent/'certificate.json.gz'
 report=json.loads(reportpath.read_text());assert report['all_passed'];assert hashlib.sha256(certpath.read_bytes()).hexdigest()==report['certificate_sha256']
 for path,h in report['input_sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
 fresh=subprocess.run([sys.executable,str(producer),'--output',str(out/'producer-replay')],cwd=ROOT,check=True,capture_output=True,text=True);replay=json.loads(fresh.stdout);assert replay['all_passed'] and replay['branches']==report['branches']
 qpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz';cpath=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/cubic-initials.json.gz';ppath=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz'
 q=read(qpath);cub=read(cpath);planes=read(ppath);cert=read(certpath);assert q['variables']==cub['variables']==cert['variables'];assert len(q['variables'])==35
 qs=[[(i,j,Fraction(a)) for i,j,a in row] for row in q['all350_quadratic_rows']];cs=[[(tuple(mon),Fraction(a)) for mon,a in row] for row in cub['raw350_cubic_terms']];assert len(qs)==len(cs)==350
 evidence=[]
 for plane,stored in zip(planes['branches'],cert['branches']):
  assert plane['branch']==stored['branch'];free=[q['variables'].index(v) for v in plane['free_variables']];normal=[i for i in range(35) if i not in free+[23]];assert len(free)==14 and len(normal)==20
  assert free==stored['parameter_indices'] and normal==stored['normal_indices'];p={q['variables'].index(k):Fraction(v) for k,v in stored['parameter_point'].items()};p[23]=Fraction(1)
  v=[sum((Fraction(a)*prod(p[i] for i in mon) for mon,a in row),Fraction(0)) for row in plane['all35_graph']];assert list(map(str,v))==stored['direction_v35'];assert v[23]==1 and all(v[i]==p[i] for i in free)
  assert all(sum(a*v[i]*v[j] for i,j,a in row)==0 for row in qs)
  Lrows=[]
  for qr,cr in zip(qs,cs):
   deriv=[Fraction(0)]*35
   for i,j,a in qr:deriv[i]+=a*v[j];deriv[j]+=a*v[i]
   radial=sum((a*prod(v[i] for i in mon) for mon,a in cr),Fraction(0));Lrows.append([deriv[i] for i in normal]+[radial])
  L=dm(Lrows);assert (L-saved(stored['linear_matrix_350x21'])).is_zero_matrix
  pc=stored['linear_pivot_columns'];pr=stored['linear_pivot_rows'];assert len(pc)==len(pr)==19 and 20 in pc;A=L.extract(pr,pc);det=A.det();assert str(s.QQ.to_sympy(det))==stored['linear_minor'] and det!=0
  fc=[i for i in range(21) if i not in pc];assert fc==stored['free_normal_columns'];pv=-A.inv()*L.extract(pr,fc);Krows=[[s.QQ.zero,s.QQ.zero] for _ in range(21)]
  for i,row in zip(pc,pv.to_list()):Krows[i]=row
  for j,i in enumerate(fc):Krows[i][j]=s.QQ.one
  K=DomainMatrix(Krows,(21,2),s.QQ);assert (L*K).is_zero_matrix and (K-saved(stored['normalized_kernel_21x2'])).is_zero_matrix;assert all(a==0 for a in Krows[20])
  source=[[Fraction(0)]*35 for _ in range(2)]
  for j in range(2):
   for k,i in enumerate(normal):source[j][i]=Fraction(s.QQ.to_sympy(Krows[k][j]))
  B=[]
  for row in qs:
   left=sum((a*source[0][i]*source[0][j] for i,j,a in row),Fraction(0));cross=sum((a*(source[0][i]*source[1][j]+source[1][i]*source[0][j]) for i,j,a in row),Fraction(0));right=sum((a*source[1][i]*source[1][j] for i,j,a in row),Fraction(0));B.append([left,cross,right])
  BM=dm(B);assert (BM-saved(stored['raw_quadratic_350x3'])).is_zero_matrix;P=BM-L.extract(range(350),pc)*A.inv()*BM.extract(pr,range(3));assert (P-saved(stored['implicit_projected_quadratic_350x3'])).is_zero_matrix
  assert P.extract(pr,range(3)).is_zero_matrix;qr=stored['quadratic_minor_rows'];qdet=P.extract(qr,range(3)).det();assert qdet!=0 and str(s.QQ.to_sympy(qdet))==stored['quadratic_minor']
  evidence.append(dict(branch=plane['branch'],all350_linear_rows_rebuilt=True,entire_normalized_kernel_rebuilt=True,linear_rank_upper_by_two_independent_kernel_vectors=19,linear_minor=str(s.QQ.to_sympy(det)),radial_kernel_zero=True,all350_raw_quadratic_rows_rebuilt=True,all350_implicit_corrected_rows_rebuilt=True,quadratic_rank=3,quadratic_minor=str(s.QQ.to_sympy(qdet)),closed_fiber_length=3))
 inputs=[Path(__file__),producer,reportpath,certpath,qpath,cpath,ppath,ROOT/'research/notes/Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_v63.md']
 result=dict(id='Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_v63_REVIEW',all_passed=True,field='QQ',producer_fresh_replay=True,branches=evidence,formal_proof_accepted=True,complete_finiteness='closed fiber m² plus equivalent maximal/base-adic topologies gives finite <=3-generator R-module',saturation_completion='Noetherian principal colons stabilize; flat completion commutes with each multiplication kernel',dimension_argument='r-regular saturated finite R-module => dimT<=14 and dimT/r<=13',accepted_result='Each specified15plane meets actual tangent cone support in proper closed cone of dimension<=14',scope='Does not remove every plane point or the certificate point, and gives no global full35 dimension bound; high-degree>=5 interpretation additionally uses separate complete-quartic result.',input_sha256={str(a.relative_to(ROOT)):hashlib.sha256(a.read_bytes()).hexdigest() for a in inputs},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('BLOWUP FINITE TRANSVERSE INDEPENDENT REVIEW PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
