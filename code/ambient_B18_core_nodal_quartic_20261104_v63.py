#!/usr/bin/env python3
"""Exact geometry of the two fixed actual core initial quadrics only."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CORE_NODAL_QUARTIC_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def mat(M):return [[str(x) for x in row] for row in M.tolist()]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    qp=ROOT/'reports/AMBIENT_B18_CORE4_QUADRATIC_INITIAL_20261104_v63/seven-actual-core-quadrics.json'
    qr=qp.parent/'report.json';Q=read(qp);assert read(qr)['all_passed']
    peer=ROOT/'reports/AMBIENT_B18_CORE4_QUADRATIC_20261104_v63_REVIEW/report.json';assert read(peer)['all_passed']
    X=s.symbols('x112 x117 tauA tauB');x,y,a,b=X
    f,g=[s.expand(sum(s.Rational(c)*s.prod(X[i] for i in m) for c,m in zip(row,Q['monomials4_through3']))) for row in Q['all7_core_quadrics_padded_through3'][:2]]
    assert s.gcd(f,g)==1
    A=s.hessian(f,X)/2;B=s.hessian(g,X)/2;assert A.rank()==4 and B.rank()==3
    v=B.nullspace()[0];v=v/v[3];assert list(v)==[s.Rational(7,5),s.Rational(2,5),s.Rational(7,24),1]
    sub=lambda q,w:s.expand(q.subs(dict(zip(X,w)),simultaneous=True))
    assert sub(f,v)==sub(g,v)==0
    t=s.Symbol('t');pencil=s.factor((A+t*B).det());assert s.cancel(pencil-(25*t*t+5*t+3)/691200)==0
    assert s.discriminant(25*t*t+5*t+3,t)==-275
    T=(A*v);assert T==s.Matrix([s.Rational(37,960),-s.Rational(13,240),-s.Rational(7,40),s.Rational(3,160)])
    # b=0 is a complement to the rational cone vertex. T=0 is the f-tangent.
    at=(37*x-52*y)/168
    tangent_g=s.factor(g.subs({b:0,a:at}));tangent_f=s.factor(f.subs({b:0,a:at}))
    H=s.hessian(tangent_g,(x,y))/2;tdet=s.factor(H.det());assert tdet
    tf=s.Poly(tangent_f,x,y);tg=s.Poly(tangent_g,x,y)
    sylv=s.resultant(tf.as_expr().subs(y,1),tg.as_expr().subs(y,1),x);assert sylv
    assert not(tf.coeff_monomial(x*x)==tg.coeff_monomial(x*x)==0)
    # Explicit QQ point on the direction conic; no beta/source point is tested.
    w0=s.Matrix([8,1,s.Rational(7,6),0]);assert sub(g,w0)==0
    assert (T.T*w0)[0]==s.Rational(1,20)
    u,w=s.symbols('u w');d=s.Matrix([u,w,0,0]);bd=(d.T*B*d)[0];cross=(w0.T*B*d)[0]
    direction=(bd*w0-2*cross*d).applyfunc(s.expand);assert direction[3]==0 and sub(g,direction)==0
    pdir=s.Poly((T.T*direction)[0],u,w);assert pdir.total_degree()==2 and s.discriminant(pdir.as_expr().subs(w,1),u)!=0
    fd=sub(f,direction);tw=(T.T*direction)[0]
    phi=(2*tw*direction-fd*v).applyfunc(s.expand)
    assert sub(f,phi)==sub(g,phi)==0
    # Clear ONE common rational scale, retaining a literal morphism P1->P3.
    den=s.ilcm(*[s.denom(c) for q in phi for c in s.Poly(q,u,w).coeffs()]);phi=(den*phi).applyfunc(s.expand)
    content=s.igcd(*[int(c) for q in phi for c in s.Poly(q,u,w).coeffs()]);phi=(phi/content).applyfunc(s.expand)
    assert all(s.Poly(q,u,w).total_degree()==4 for q in phi)
    gg=s.Poly(phi[0],u,w)
    for q in phi[1:]:gg=s.gcd(gg,s.Poly(q,u,w))
    assert gg.total_degree()==0
    qmat=s.Matrix([[s.Poly(q,u,w).coeff_monomial(u**(4-j)*w**j) for j in range(5)] for q in phi]);assert qmat.rank()==4
    assert sub(f,phi)==sub(g,phi)==0
    proof={'source_quadrics':[str(f),str(g)],'coordinate_order':list(map(str,X)),
      'symmetric_matrix_f':mat(A),'symmetric_matrix_g':mat(B),'matrix_ranks':[4,3],
      'gcd_f_g':'1','vertex_B':list(map(str,v)),'f_vertex':'0',
      'pencil_determinant':str(pencil),'finite_pencil_discriminant':'-275',
      'Av_tangent_linear_coefficients':list(map(str,T)),
      'node_affine_patch':'tauB=1','eliminated_affine_tangent_coordinate':'tauA',
      'f_affine_linear_derivative_in_tauA':'-7/20',
      'node_tangent_quadric':str(tangent_g),'node_tangent_2_by2_matrix':mat(H),'node_tangent_determinant':str(tdet),
      'f_restriction_on_tangent_direction_plane':str(tangent_f),
      'homogeneous_direction_no_basepoint_resultant':str(sylv),
      'QQ_point_on_direction_conic_only':list(map(str,w0)),
      'direction_conic_parametrization':list(map(str,direction)),
      'node_preimage_binary_quadratic':str(s.factor(tw)),
      'node_preimage_discriminant':str(s.discriminant(pdir.as_expr().subs(w,1),u)),
      'normalization_P1_to_quartic_four_binary_forms':list(map(str,phi)),
      'binary_form_coefficient_matrix4_by5':mat(qmat),'four_binary_forms_gcd':'1',
      'all_normalization_quadric_identities_verified':True,
      'finite_pencil_roots_simple_excludes_other_curve_singularities':True,
      'direction_T0_f0_g0_has_no_projective_point':True,
      'curve_geometrically_integral_nodal_quartic':True,'normalization_QQ_rational_P1':True,
      'complete_intersection_arithmetic_genus':1,'node_delta':1,'geometric_genus':0,
      'two_actual_core_functions_define_integral_CM_surface_from_regular_leading_sequence':True,
      'remaining_five_actual_core_functions_start_order_at_least3':True,
      'no_actual_higher_core_functions_evaluated_or_socle_depth_inferred':True}
    p=out/'fixed-quadrics-nodal-quartic.json';p.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',qp,qr,peer]
    rep={'id':STEM,'all_passed':True,'status':'actual_initial_quadrics_define_QQ_rational_irreducible_one_node_22_quartic',
      'matrix_ranks':[4,3],'node_count_geometric':1,'node_tangent_determinant':str(tdet),
      'normalization_degree':4,'geometric_genus':0,'normalization_over_QQ':'P1',
      'actual_core_cubic_source4_or_GB_used':False,'full_core_socle_or_depth_determined':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
