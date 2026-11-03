#!/usr/bin/env python3
"""Exact finite transverse certificates in actual blowup charts over15-planes.

For G_i=F_i(r*(v(p)+n))/r^2, the19-rank linear matrix is[dQ_v,C(v)].
Its2-dimensional kernel hasr=0.  Thus restricted quadratics are exactly
Q_i(kernel), independently of D_i and higher terms.  Full rank3 after
implicit projection proves transverse length3 and finiteness over14 base
parameters.  Saturation/dimension consequences are proved in the note.
"""
import argparse,gzip,hashlib,json,random,time
from fractions import Fraction as F
from math import prod
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from q_true_initial_eta_open_20261103_v63 import unpack
ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 qpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz';cpath=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/cubic-initials.json.gz'
 ppath=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz'
 q=json.loads(gzip.decompress(qpath.read_bytes()));c=json.loads(gzip.decompress(cpath.read_bytes()));planes=json.loads(gzip.decompress(ppath.read_bytes()))
 Q=sp.QQ;qs=[{(i,j):F(a)for i,j,a in row}for row in q['all350_quadratic_rows']];cs=[unpack(r)for r in c['raw350_cubic_terms']]
 ev=lambda poly,v:sum((a*prod(v[i]for i in m)for m,a in poly.items()),F(0));cv=lambda a:Q(a.numerator,a.denominator)
 stringify=lambda M:[[str(Q.to_sympy(a))for a in row]for row in M.to_list()]
 evidence=[];summaries=[]
 for branch in planes['branches']:
  graph=[{(23,)if not m else tuple(m):F(v)for m,v in d}for d in branch['all35_graph']]
  free=[q['variables'].index(n)for n in branch['free_variables']];normal=[i for i in range(35)if i not in free+[23]]
  rng=random.Random(20261103);point={i:F(rng.randrange(1,8))for i in free+[23]};point[23]=F(1)
  v=[ev(p,point)for p in graph];rows=[]
  for poly,cb in zip(qs,cs):
   r=[F(0)]*35
   for(i,j),a in poly.items():r[i]+=a*v[j];r[j]+=a*v[i]
   rows.append([cv(r[i])for i in normal]+[cv(ev(cb,v))])
  L=DomainMatrix(rows,(350,21),Q);_,cp=L.rref();assert len(cp)==19
  _,rp=L.extract(range(350),cp).transpose().rref();minor=L.extract(rp,cp);det=minor.det();assert det
  null=L.nullspace().transpose();assert null.shape==(21,2)and all(not a for a in null.to_list()[-1])
  nonpiv=[i for i in range(21)if i not in cp];assert len(nonpiv)==2 and 20 not in nonpiv
  null=null*null.extract(nonpiv,[0,1]).inv();assert (L*null).is_zero_matrix
  assert (null.extract(nonpiv,[0,1])-DomainMatrix.eye(2,Q)).is_zero_matrix
  vectors=[]
  for col in range(2):
   vec=[F(0)]*35
   for j,i in enumerate(normal):vec[i]=F(Q.to_sympy(null.to_list()[j][col]))
   vectors.append(vec)
  rawcols=[]
  for direction in [vectors[0],[a+b for a,b in zip(*vectors)],vectors[1]]:
   rawcols.append(DomainMatrix([[cv(ev(poly,direction))]for poly in qs],(350,1),Q))
  rawcols[1]=rawcols[1]-rawcols[0]-rawcols[2];quadratic=rawcols[0].hstack(rawcols[1],rawcols[2])
  # Actual implicit correction kills the19 selected equations.  This is the
  # full QQ coefficient matrix of the remaining three quadratic monomials.
  projected=quadratic-L.extract(range(350),cp)*minor.inv()*quadratic.extract(rp,range(3))
  _,piv=projected.rref();assert len(piv)==3
  _,qr=projected.transpose().rref();qdet=projected.extract(qr,range(3)).det();assert qdet
  evidence.append({'branch':branch['branch'],'parameter_indices':free,'normal_indices':normal,
      'parameter_point':{q['variables'][i]:str(point[i])for i in free},'direction_v35':list(map(str,v)),
      'linear_matrix_columns':[q['variables'][i]for i in normal]+['radial_r'],
      'linear_matrix_350x21':stringify(L),'linear_pivot_columns':list(cp),'linear_pivot_rows':list(rp),'linear_minor':str(Q.to_sympy(det)),
      'normalized_kernel_21x2':stringify(null),'free_normal_columns':nonpiv,
      'raw_quadratic_350x3':stringify(quadratic),'implicit_projected_quadratic_350x3':stringify(projected),
      'quadratic_minor_rows':list(qr),'quadratic_minor':str(Q.to_sympy(qdet))})
  summaries.append({'branch':branch['branch'],'normal_variables':20,'base_parameters':14,'radial_variables':1,
      'linear_rank':19,'kernel_dimension':2,'kernel_radial_coordinate_zero':True,'projected_quadratic_rank':3,
      'linear_minor':str(Q.to_sympy(det)),'quadratic_minor':str(Q.to_sympy(qdet)),
      'actual_weak_blowup_transverse_fiber_length':3,'weak_blowup_completion_finite_over14_base':True,
      'saturated_exceptional_completed_local_dimension_upper_bound':13,
      'entire15plane_not_contained_actual_tangent_cone':True})
 cpfile=args.output/'certificate.json.gz'
 with gzip.open(cpfile,'wt')as out:json.dump({'variables':q['variables'],'branches':evidence},out,separators=(',',':'))
 inputs=[Path(__file__),qpath,cpath,ppath,ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/report.json',ROOT/'reports/Q_ORDINARY_QUARTIC_INITIAL_20261103_v63/report.json']
 report={'id':'Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_V63','all_passed':True,'coefficient_field':'QQ','branches':summaries,
      'method':'Actual blowup weak transform F(r*(v(p)+n))/r^2; implicit19; transverse quadratic initial ideal m2; finite-over14base then r-saturation and dimension.',
      'Puiseux_directions_excluded_by_ordinary_arc_failure_alone':False,
      'full_plane_generic_QJac_rank18_claimed':False,'all15plane_points_excluded':False,
      'local_full_R35_dimension_upper12_claimed':False,'actual_new_branch_dimension_claimed':False,
      'certificate_sha256':hashlib.sha256(cpfile.read_bytes()).hexdigest(),
      'input_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},'seconds':time.monotonic()-started}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
