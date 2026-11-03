#!/usr/bin/env python3
"""Independent whole-family geometry and shear-normalization audit."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63';r=json.loads((folder/'report.json').read_text());assert r['all_passed']
 for name,want in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want
 source=ROOT/'reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/certificate.json.gz';old=json.loads(gzip.decompress(source.read_bytes()))
 n0,n1,n2,h,x,y,z,w,t,rr=s.symbols('n0 n1 n2 h x y z w t rr');V=s.Matrix([h,x,y,z,w]);coords=[x,y,z,w];pairs=list(itertools.combinations(range(4),2))
 P=s.Matrix(old['positive_P']);graph=list(map(s.sympify,old['exact99_polynomial_graph']));mons=list(map(s.sympify,old['coefficient_monomials']))
 for k,(i,j)in enumerate(pairs):
  value=sum(c*m for c,m in zip(graph[15*k:15*k+15],mons));P[i,j]+=value;P[j,i]-=value
 P=P.applyfunc(s.expand);assert P==s.Matrix(r['affine_P'])
 for i,j,k in itertools.combinations(range(4),3):assert s.expand(sum(P[i,a]*s.diff(P[j,k],coords[a])+P[j,a]*s.diff(P[k,i],coords[a])+P[k,a]*s.diff(P[i,j],coords[a])for a in range(4)))==0
 pf=s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]);assert s.Poly(pf,*coords).total_degree()==2
 Q=s.expand(h*h*pf.subs(dict(zip(coords,[v/h for v in coords])),simultaneous=True));A=s.hessian(Q,V)/2
 assert s.factor(A.det())==s.Rational(1024,531441)*n0**3*n1*(32*n0+n1*n2)
 assert A.subs(h,0)==A
 restricted=A.extract([1,2,3,4],[1,2,3,4]);assert all(restricted[2,i]==0 for i in range(4))
 assert s.factor(restricted.extract([0,1,3],[0,1,3]).det())==-s.Rational(4,729)*n0*n1*(32*n0+n1*n2)
 assert [s.expand(s.diff(Q,v).subs({h:0,x:0,y:0,z:1,w:0}))for v in V]==[-s.Rational(32,27)*n0,0,0,0,0]
 raw=s.zeros(5)
 for i,j in pairs:
  raw[i+1,j+1]=s.expand(h*h*P[i,j].subs(dict(zip(coords,[v/h for v in coords])),simultaneous=True));raw[j+1,i+1]=-raw[i+1,j+1]
 div=s.Matrix([s.expand(sum(s.diff(P[i,j],coords[i])for i in range(4)))for j in range(4)]);N=s.zeros(5);N[1:,0]=div.subs(dict.fromkeys(coords,0));N[1:,1:]=div.jacobian(coords)
 assert N==s.Matrix(r['homogeneous_modular_N'])and N**3==s.zeros(5)
 assert (N*N)[3,0]==4*n0*n0*n2/3
 assert s.factor(N.extract([2,3],[0,1]).det())==4*n0*n0/3
 G=s.eye(5)+t*N+t*t*N*N/2;Gi=s.eye(5)-t*N+t*t*N*N/2
 assert (G*raw.subs(dict(zip(V,Gi*V)),simultaneous=True)*G.T-raw).applyfunc(s.expand)==s.zeros(5)
 assert (N.T*A+A*N).applyfunc(s.expand)==s.zeros(5)
 # The z-h shear is ESSENTIAL. A bare diagonal fails the constant bracket.
 for eps in [1,-1]:
  D=s.diag(1,rr,rr*rr,rr**3,eps*rr);D[3,0]=(rr**3-eps)/(6*n0)
  moved=eps/rr**2*D*raw.subs(dict(zip(V,D.inv()*V)),simultaneous=True)*D.T
  expected=raw.subs({n0:eps*n0,n1:eps*rr*n1,n2:n2/rr},simultaneous=True)
  assert (moved-expected).applyfunc(s.cancel)==s.zeros(5)
 for key in ['orbit_sample','orbit_plus_scalar_sample','family_source_sample']:
  cert=r[key];M=s.Matrix(cert['matrix']);C=s.Matrix(cert['all_column_coefficients']);rows=cert['rows'];columns=cert['columns']
  assert M[:,columns]*C==M and M[rows,columns].det()==s.sympify(cert['determinant'])!=0
 assert [r[k]['rank']for k in ['orbit_sample','orbit_plus_scalar_sample','family_source_sample']]==[23,24,26]
 paths=[Path(__file__),folder/'report.json',source,ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',ROOT/'research/notes/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63.md',ROOT/'research/notes/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63_REVIEW.md']
 result=dict(id='R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63_REVIEW',all_passed=True,
  full_three_parameter_source_reconstructed_without_producer_imports=True,
  generic_divisor='3H+Q smooth rank5, H tangent at [0:0:0:1:0]',
  generic_modular_Jordan=[3,1,1],exact_finite_Ga_flow=True,
  both_signs_all150_homogeneous_shear_parameter_identities=True,
  z_h_shear_is_essential_and_no_bare_diagonal_claim=True,
  two_parameter_dominant_normal_form_cover_on_n0_n1_unit=True,
  generic_projective_orbit_dimension=23,connected_stabilizer='Ga',
  complete_orbit_equivalence_classification_or_generic_scheme_thickness_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT NEW25 GEOMETRY/FULL SHEAR REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
