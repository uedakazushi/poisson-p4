#!/usr/bin/env python3
"""The full nonzero Qa line stays in the common Ac8/Ac7 pinched strata.

Exact polynomial all-a projective flag ideals and existing whole-image odd
module certificates distinguish this from the ambient resonance at a=2/9.
"""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,pushforward,coefficients
ROOT=Path(__file__).resolve().parents[1]
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();a=s.Symbol('a');sourcepath=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json';body=json.loads(sourcepath.read_text());base=list(map(s.Symbol,body['base_parameters']));fiber=list(map(s.Symbol,body['fiber_parameters']));point=dict(zip(base,[0,0,s.Rational(4,3),s.Rational(3,2),0,0]))|dict.fromkeys(fiber,0);point[fiber[1]]=9*a/2
 assert str(fiber[1])=='k0';uCA=s.Matrix(list(map(s.sympify,body['families']['CA']['normalized_61_coordinates'])));uQ=uCA.subs(point);assert [(i+1,v)for i,v in enumerate(uQ)if v]!=[];assert [(i+1,v)for i,v in enumerate(uQ)if v]==[(25,-s.Rational(4,3)),(36,s.Integer(3)),(38,9*a/2)]
 model=S61();P=model.reconstruct(uQ);verify=model.verify_point(dict(parameters=list(uQ),coefficients=coefficients(P)),check_equations=True)
 # Membership in both actual closed families, including N, is exact before
 # transporting the line flag.
 for name,family in body['families'].items():
  ub=s.Matrix(list(map(s.sympify,family['normalized_61_coordinates']))).subs(point);assert ub==uQ,name
 r,z,w=s.symbols('r z w');g=s.Matrix([[1,r,-w,z],[0,1,0,0],[0,z,1,0],[0,w,0,1]]);full=s.eye(5);full[1:5,1:5]=g.inv();normal=model.normalize(pushforward(P,full));u=s.Matrix(normal['parameters']);bi=body['beta_indices_zero_based'];qi=body['q_indices_zero_based'];records={};constraints={}
 for name,family in body['families'].items():
  params=s.Matrix(family['beta_left_inverse'])*u[bi,0];shift=s.Matrix(list(map(s.sympify,family['q_shift']))).subs(dict(zip(base,params)),simultaneous=True);bc=s.Matrix(family['beta_closed_constraints'])*u[bi,0];qc=s.Matrix(family['fiber_closed_constraints'])*(u[qi,0]-shift);rows=[s.expand(v)for v in [*bc,*qc]if s.expand(v)!=0];constraints[name]=rows
 constraints['N']=constraints['CA']+constraints['CB']
 for name,rows in constraints.items():
  gb=s.groebner(rows,r,z,w,a,order='lex',domain=s.QQ);assert list(gb)==[r,z,w*w]
  records[name]=dict(all_a_constraints=list(map(str,rows)),entire_QQ_a_flag_ideal=list(map(str,gb)),whole_geometric_fiber_points=1,whole_scheme_fiber_length=2)
 # No other projective base flag: beta is constant, Kbeta has kernel span(e_y,e_w),
 # and the sole off-chart point e_w fails an a-independent base equation.
 off=s.eye(5);off[1:5,1:5]=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]);uoff=s.Matrix(model.normalize(pushforward(P,off))['parameters']);witness=s.Matrix(body['families']['CA']['beta_closed_constraints'])*uoff[bi,0];assert any(v!=0 for v in witness) and not any(v.has(a)for v in witness)
 oldpath=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json';old=json.loads(oldpath.read_text());odd={}
 for name,fa in old['families'].items():
  pars=list(map(s.Symbol,fa['source_fiber_parameters']));subs=dict.fromkeys(pars,0);subs[s.Symbol('k0')]=a;values=[s.expand(s.sympify(v).subs(subs))for v in fa['first_odd_coefficient_at_fixed_target_K']];assert all(v==0 for v in values);odd[name]=dict(free_varying_coordinate='k0',all45_first_odd_coefficients_zero=list(map(str,values)),whole_image_missing_ideal=fa['missing_module_ideal_in_target_residue_coordinates'],pinching_codimension=fa['pinching_codimension'])
 paths=[Path(__file__),sourcepath,oldpath,ROOT/'reports/Q_CA_CB_N_SOURCE_SMITH_20261103_v63/report.json',ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/certificate.json',ROOT/'research/notes/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63.md',ROOT/'code/s61_model.py',ROOT/'research/notes/AMBIENT_N_GLOBAL_R1_20261103_v63.md',ROOT/'research/notes/AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63.md']
 result=dict(id='Q_CA_CB_N_ALL_NONZERO_COMPLETED_MODELS_20261103_v63',all_passed=True,field='QQ[a]',S61_line=list(map(str,uQ)),exact134_150_Jacobi_verification=verify,physical_Qa_identification_scope='a!=0; S61 extension at a0 is not the physical rank2 Qa0',base_fixed_point={str(k):str(v)for k,v in point.items()if k in base},varying_fiber_coordinate='k0=9a/2',all_a_flag_fibers=records,no_off_chart_base_flags=True,off_chart_a_independent_witness=list(map(str,witness)),whole_image_odd_ideals_on_line=odd,all_nonzero_physical_Qa_models=dict(CA='A8 times smooth9',CB='A8 times smooth9',N='A7 times smooth9'),projective_completed_invariants=dict(CA=dict(dimension=28,embedding_dimension=37,depth=20,multiplicity=2),CB=dict(dimension=28,embedding_dimension=37,depth=20,multiplicity=2),N=dict(dimension=27,embedding_dimension=35,depth=20,multiplicity=2)),only_a1_pinch_claimed=False,physical_a0_completed_type_claimed=False,physical_a0_nonnormal_nonCM_by_closedness_and_global_R1=True,scope='Reduced CA/CB/N components; not whole Poisson ring, whose resonance at2/9 remains independent.',input_sha256={str(v.relative_to(ROOT)):hashlib.sha256(v.read_bytes()).hexdigest()for v in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ALL NONZERO Qa SAME FULL CA/CB/N PINCHED MODELS PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
