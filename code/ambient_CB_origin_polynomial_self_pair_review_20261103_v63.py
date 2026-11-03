#!/usr/bin/env python3
"""Independent physical-source and full58-tangent audit of the new CB pair."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1];PARENT='AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63';STEM=PARENT+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports'/PARENT;report=load(folder/'report.json');assert report['all_passed']
 for n,h in report['input_sha256'].items():assert sha(ROOT/n)==h,n
 for n,h in report['evidence_sha256'].items():assert sha(folder/n)==h,n
 data=load(folder/'certificate.json.gz');h,x,y,z,w,b=s.symbols('h x y z w b');X=s.Matrix([h,x,y,z,w]);loc={str(v):v for v in [*X,b]}
 parse=lambda q:s.sympify(q,locals=loc)
 mat=lambda q:s.Matrix([[parse(v) for v in row] for row in q])
 col=lambda q:s.Matrix([parse(v) for v in q])
 K=s.QQ.frac_field(b);dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(K)
 same=lambda A,B: A.shape==B.shape and (A-B).applyfunc(s.cancel)==s.zeros(*A.shape)
 pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
 def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[c]) for i,j in pairs for a,c in mons])
 def tensor(v):
  P=s.zeros(5)
  for k,(i,j) in enumerate(pairs):P[i,j]=sum(v[15*k+m]*X[a]*X[c] for m,(a,c) in enumerate(mons));P[j,i]=-P[i,j]
  return P
 def canonical(P):
  d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
  return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
 def rawCB(f,g,H):
  P=s.zeros(5);P[1,2]=f;P[3,4]=g;P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z)
  return P-P.T
 def rawCA(f=0,g=0,H=0):
  P=s.zeros(5);P[0,2]=f;P[3,4]=g*x*x;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z);return P-P.T
 def lie(i,j,P):
  E=s.zeros(5);E[i,j]=1
  return (E*P+P*E.T-X[j]*P.diff(X[i])).applyfunc(s.expand)
 def push(G,P):
  inverse=G.inv();assert G*inverse==inverse*G==s.eye(5)
  return (G*P.subs(dict(zip(X,inverse*X)),simultaneous=True)*G.T).applyfunc(s.cancel)
 # Reconstruct actual P0 and fixed CA20; Q_b's independent polynomial formula.
 cubicCA=[s.prod((x,h,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
 ca=s.Matrix.hstack(*[vec(canonical(rawCA(f=f))) for f in [x*x,h*x,h*h]],vec(canonical(rawCA(g=1))),*[vec(canonical(rawCA(H=H))) for H in cubicCA])
 capath=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';old=load(capath)
 assert same(ca,s.Matrix(old['fixed_fiber_linear_basis150_by20']))
 p0=ca*s.Matrix(old['center_fixed20']);q=p0-b*(2*ca[:,1]+ca[:,9]);assert same(q,col(data['full_polynomial_Q150']))
 V=s.Matrix([x*x+s.Rational(8,5)*b*h*x,-s.Rational(2,5)*b*x*x,0,h*x+s.Rational(3,5)*b*x*z,h*h+z*z-s.Rational(7,5)*b*x*w]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
 assert same(vec(P),q) and canonical(P)==P
 jacobi=[]
 for i,j,k in itertools.combinations(range(5),3):
  J=s.expand(sum(P[i,l]*s.diff(P[j,k],X[l])+P[j,l]*s.diff(P[k,i],X[l])+P[k,l]*s.diff(P[i,j],X[l]) for l in range(5)))
  jacobi.extend(s.Poly(J,*X).coeff_monomial(s.prod(X[v] for v in mon)) for mon in itertools.combinations_with_replacement(range(5),3))
 assert len(jacobi)==350 and not any(jacobi)
 assert s.Poly(h*P[1,2]+x*P[2,0]+y*P[0,1],*X).coeff_monomial(x**3)==1
 G1=s.Matrix([[0,1,0,0,0],[b,1,0,-b*b,0],[0,0,1,0,0],[1,0,0,0,0],[0,0,0,0,1]])
 G2=s.eye(5);G2[0,0]=2*b;G2[0,1]=1
 assert G1.det()==-b*b and G2.det()==2*b
 sources=data['actual_two_CB_sources'];assert len(sources)==2
 Hs=[-3*b*h*w*z/2-h*h*w+z**3/3+z**3/(3*b*b)+(h-x)*z*z/b**3+(h-x)**2*z/b**4,
     -b*w*x*z+z**3/3-h*w*x/(2*b)+w*x*x/(2*b)+(h-x)**2*z/(4*b*b)]
 fs=[-b*h*x,2*b*x*x]
 cbglobalpath=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json';globaldata=load(cbglobalpath)
 para=[tuple(v) for v in globaldata['parabolic_GL_positions']];flags=[pos for pos in itertools.product(range(5),repeat=2) if pos not in para];assert len(flags)==9
 cubic=[s.prod((h,x,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
 binary=[h*h,h*x,x*x]
 F=s.Matrix.hstack(*[vec(canonical(rawCB(f,0,0))) for f in binary]);GG=s.Matrix.hstack(*[vec(canonical(rawCB(0,f,0))) for f in binary]);HH=s.Matrix.hstack(*[vec(canonical(rawCB(0,0,H))) for H in cubic])
 fixed=F.row_join(HH);rr=list(dm(fixed.T).rref()[1]);assert len(rr)==19;left=fixed[rr,:].inv()
 project=lambda M:(M-fixed*left*M[rr,:]).applyfunc(s.cancel)
 # Construct the FULL mu*f map first; its mixed partials are retained.
 u=s.symbols('u0:19');mu=s.Symbol('mu');universal=F*s.Matrix(u[:3])+mu*GG*s.Matrix(u[:3])+HH*s.Matrix(u[3:]);J20=universal.jacobian([*u,mu])
 for i in range(3):assert same(universal.diff(mu).diff(u[i]),GG[:,i])
 actualsources=[];fibers=[];full=[];minors=[]
 for index,(G,f,H) in enumerate(zip([G1,G2],fs,Hs)):
  witness=sources[index];assert same(G,mat(witness['G'])) and same(G.inv(),mat(witness['G_inverse']))
  assert s.cancel(f-parse(witness['f']))==s.cancel(H-parse(witness['H']))==0
  S=push(G,P);expected=canonical(rawCB(f,0,H));assert same(S,expected)
  assert same(vec(S),col(witness['full150_source'])) and same(vec(S),col(witness['full150_pushed_Q']))
  values=[*[s.Poly(f,h,x).coeff_monomial(mon) for mon in binary],*[s.Poly(H,h,x,z,w).coeff_monomial(mon) for mon in cubic],0]
  record=data['source_rank_certificates'][index];assert all(s.cancel(a-parse(c))==0 for a,c in zip(values,record['source20_values']))
  sub=dict(zip([*u,mu],values));fiber=J20.subs(sub).applyfunc(s.cancel);assert same(universal.subs(sub),vec(S))
  assert same(fiber[:,19],GG*s.Matrix(values[:3]))
  whole=fiber.row_join(s.Matrix.hstack(*[vec(lie(i,j,S)) for i,j in flags]));assert same(whole,mat(record['full150_by29']))
  quot=project(whole[:,19:]);assert same(quot,mat(record['quotient150_by10']))
  wit=record['quotient_witness'];assert wit['rank']==10
  sq=dm(quot.extract(wit['rows'],wit['columns']));det=s.factor(K.to_sympy(sq.det()));assert det==parse(wit['minor'])
  assert s.cancel(det-[432*b**5/15625,786432*b**10/15625][index])==0
  assert dm(whole).rank()==29
  actualsources.append(S);fibers.append(fiber);full.append(whole);minors.append(str(det))
 # Distinct ORIGINAL marked lines, not just labels of two coordinate charts.
 assert G1[0,:]==s.Matrix([[0,1,0,0,0]]) and G2[0,:]==s.Matrix([[2*b,1,0,0,0]])
 orbit=s.Matrix.hstack(*[vec(lie(i,j,P)) for i,j in itertools.product(range(5),repeat=2)])
 assert same(orbit,mat(data['full_GL150_by25']))
 witnesses=[('GL_rank_witness',orbit,19),('GL_scalar_rank_witness',orbit.row_join(q),19),('GL_scalar_parameter_rank_witness',orbit.row_join(q).row_join(q.diff(b)),20)]
 for label,M,rank in witnesses:
  w0=data[label];assert w0['rank']==rank and dm(M).rank()==rank
  assert s.cancel(K.to_sympy(dm(M.extract(w0['rows'],w0['columns'])).det())-parse(w0['minor']))==0
 # Stronger independent check: pull BOTH FULL29 tangent matrices to originalQ.
 back=[]
 for G,J in zip([G1,G2],full):back.append(s.Matrix.hstack(*[vec(push(G.inv(),tensor(J[:,j]))) for j in range(29)]))
 assert dm(back[0].row_join(back[1])).rank()==37
 relative=G1*G2.inv();second20=s.Matrix.hstack(*[vec(push(relative,tensor(fibers[1][:,j]))) for j in range(20)])
 cross=project(fibers[0][:,19].row_join(s.Matrix.hstack(*[vec(lie(i,j,actualsources[0])) for i,j in flags])).row_join(second20))
 assert same(cross,mat(data['same_point_cross_quotient150_by30']))
 cw=data['cross_quotient_witness'];assert cw['rank']==18 and dm(cross).rank()==18
 cdet=s.factor(K.to_sympy(dm(cross.extract(cw['rows'],cw['columns'])).det()));assert s.cancel(cdet-82944*b**25/15625)==0
 assert data['branch_tangent_union_rank']==37 and data['branch_tangent_intersection_dimension']==21
 proper=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';fiberproof=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json'
 assert load(proper)['all_passed'] and load(fiberproof)['all_passed']
 files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',folder/'report.json',folder/'certificate.json.gz',capath,cbglobalpath,proper,fiberproof]
 files += [ROOT/n for n in report['input_sha256']]
 result=dict(id=STEM,all_passed=True,field='QQ(b)',new_family_distinct_from_Pa=True,
  own_physical_Qb_and_full_CA20_binding=True,all350_strict_Jacobi_and25_divergence_zero=True,
  both_finite5x5_actual_CB_source_identities_all150=True,actual_flags_distinct_for_every_b_nonzero=True,
  universal_mu_times_f_differentiated_before_mu0=True,all_three_mixed_partials150_retained=True,
  both_actual_source150_by29_reconstructed=True,both_source_ranks=[29,29],both_quotient10_minors=minors,
  own_GL25_rank=19,own_GL_scalar_rank=19,own_GL_scalar_parameter_rank=20,projective_saturation_dimension=19,
  both_full29_tangents_pushed_to_original_Q=True,all58_union_rank=37,tangent_intersection_dimension=21,cross_minor=str(cdet),
  generic_whole_fiber_exactly_two_reduced_points_by_proper_open=True,
  generic_reduced_CB_completion='T/(I1 intersection I2), T/Ii = C[[29]]',generic_cone_depth_at_most=22,generic_projective_depth_at_most=21,
  meeting_scheme_edim=21,meeting_scheme_dimension_bounds=[20,21],
  all_nonzero_whole_flag_fiber_or_origin_conductor_multiplicity_claimed=False,
  full_CB_origin_odd_ideal_or_actual_CA_CB_intersection_D_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(files)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',STEM,result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
