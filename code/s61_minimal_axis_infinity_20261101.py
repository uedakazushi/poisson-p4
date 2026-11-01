#!/usr/bin/env python3
"""A dominant sixteen-dimensional component of isolated M0 infinity points."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J5
from s61_nil22_infinity_20261101 import tensor,sm
from s61_known_fiber_audit import equation_jacobian
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def run():
    m=S61();weights=m.reduced['parameter_weights'];one=[i for i,t in enumerate(weights) if t==1];two=[i for i,t in enumerate(weights) if t==2];allu=s.symbols('u1:62');u=[allu[i] for i in one];x,y,z,w=X[1:];g=s.symbols('g0:9');k,L,U,V=s.symbols('k L U V');D=(3*k+2)*L-(6*k+8)*U
    G=sum(a*b for a,b in zip(g,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3],strict=True));Q=tensor([x*x,0,0,-s.diff(G,w),s.diff(G,z),0]);B=X[0]**2*J5;B[1:,1:]+=Q;normalized=m.normalize(B)['parameters'];values=[allu[i] if i in one else normalized[i] for i in range(61)]
    beta_sub={allu[3]:0,allu[11]:0,allu[4]:0,allu[12]:0,allu[13]:k,allu[48]:2*k+1,allu[25]:U,allu[34]:V,allu[49]:-(9*k+7)*V/(3*k+1),allu[50]:5*U-L}
    gv={g[0]:-2*V*(D-8*U)/((3*k+1)*D),g[1]:2*U*(L-2*U)/D,g[2]:-4*V*(3*k+2)/((3*k+1)*D),g[3]:-8*U/D,g[4]:0,g[5]:0,g[6]:2*(3*k+2)/D,g[7]:0,g[8]:0}
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';extra=json.loads(extra_path.read_text())['universal_new_equation'];top=[[term for term in f if sum(i in two for i in term[0])==2] for f in m.reduced['equations']+[extra]]
    restrictions=[s.expand(sum(s.Rational(n,d)*s.prod(values[i] for i in mon) for mon,n,d in f)) for f in top];assert restrictions[:70]==[0]*70
    assert all(s.cancel(f.subs(beta_sub).subs(gv,simultaneous=True))==0 for f in restrictions)
    vertex={allu[3]:0,allu[11]:0,allu[4]:0,allu[12]:0,allu[13]:k,allu[48]:2*k+1};assert s.expand(restrictions[-1].subs(vertex)-9*(2*k+1)*g[6]*g[7])==0
    sample={k:1,U:2,V:3,L:5,allu[24]:2,allu[26]:3,allu[33]:5,allu[35]:7,allu[36]:11,allu[47]:13};beta_map=s.Matrix([s.sympify(beta_sub.get(ui,ui)) for ui in u]);beta_values=beta_map.subs(sample);sample_g={a:s.sympify(b).subs(sample) for a,b in gv.items()};point=[s.cancel(s.sympify(v).subs(beta_sub).subs(gv,simultaneous=True).subs(sample)) for v in values];assert all(m.evaluate(f,point)==0 for f in top)
    zz,uu,beta=primitive_beta(m);assert uu==u;vector=s.Matrix([s.diff(beta[i,j],zz[h]) for (i,j),h in COORDS]);coef=vector.jacobian(u);rows=list(coef.T.rref()[1]);left=coef[rows,:].inv();primitive=coef*beta_values
    variables=[k,U,V,L]+[ui for ui in u if ui not in beta_sub];assert len(variables)==10
    param=beta_map.jacobian(variables).subs(sample);action=s.Matrix.hstack(*[left*infinitesimal_action(primitive,A)[rows,:] for A in sp_basis()]);differential=param.row_join(action);dm=minor_witness(differential);assert dm['rank']==16
    v=s.Matrix([0,1,0,0]);line=s.Matrix.hstack(*[(-A*v)[[0,2,3],:] for A in sp_basis()]);joint=differential.col_join(s.zeros(3,10).row_join(line));jm=minor_witness(joint);assert jm['rank']==16
    leading=equation_jacobian(top,point)[:,two];w134=minor_witness(leading[:134,:]);w135=minor_witness(leading);assert w134['rank']==w135['rank']==44
    fixed=s.Matrix(restrictions).jacobian(g).subs(dict(zip(u,beta_values,strict=True))).subs(sample_g);fw=minor_witness(fixed);assert fw['rank']==9
    c4,c6,delta,j=curve_invariants_from_beta(primitive);assert delta==4267428033666687552
    local_path=ROOT/'reports/S61_MINIMAL_AXIS_LOCAL_POINT_20261101.json';local=json.loads(local_path.read_text());assert point==list(map(s.Rational,local['point61']))
    fixed_q=next(j for j,i in enumerate(two) if point[i]!=0);chart=[j for j in range(45) if j!=fixed_q]
    tcol=s.Matrix([m.evaluate([term for term in f if sum(i in two for i in term[0])==1],point) for f in m.reduced['equations']+[extra]])
    full=leading[:,chart].row_join(tcol);rr=local['minor_rows_zero_based'];minor_dm=full[rr,:].to_DM();det=minor_dm.domain.to_sympy(minor_dm.det());assert det==s.Rational(local['minor_determinant']) and det!=0
    inputs=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_known_fiber_audit.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/elliptic_curve_invariants_20261101.py',ROOT/'code/weight1_structure.py',m.source/'symplectic_slice_reduced_equations.json',extra_path,local_path]
    return {'id':'S61_MINIMAL_AXIS_INFINITY_20261101','all_passed':True,'field':'QQ','G':str(G),'Q':sm(Q),'g_formula':{str(a):str(b) for a,b in gv.items()},'beta_formula':{str(a):str(b) for a,b in beta_sub.items()},'denominator_D':str(D),'source_open':'(3*k+1)*(3*k+2)*D != 0','all_135_identities':True,'G6_on_full_nine_G_vertex_chart':'9*(2*k+1)*g6*g7','source_variables':list(map(str,variables)),'sample_beta_values':list(map(str,beta_values)),'sample_g':{str(a):str(b) for a,b in sample_g.items()},'sample_S61_coordinates':list(map(str,point)),'weight_two_u_indices_zero_based':two,'dominance_matrix':sm(differential),'dominance_minor':dm,'joint_beta_vertex_minor':jm,'leading134_minor':w134,'leading135_minor':w135,'fixed_vertex_minor':fw,'c4':str(c4),'c6':str(c6),'Delta':str(delta),'independent_full_homogenized_minor':{'rows':rr,'determinant':str(det),'rank':45},'image_component_dimension':16,'fixed_beta_projective_tangent_dimension_at_sample':0,'scope':'A dominant generically finite component of the UNSATURATED leading incidence. The sample is a reduced isolated leading point. The full-homogenized rank45 minor is independently verified, so this component is removed by saturation over sufficiently general beta; no affine-fiber upper bound or minimal-stratum exhaustion asserted here.','proof_note':'research/notes/S61_MINIMAL_AXIS_INFINITY_20261101.md','input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'beta_rank':16,'leading_rank':44,'fixed_vertex_rank':9,'Delta_nonzero':True,'output':str(args.output)}))
if __name__=='__main__':main()
