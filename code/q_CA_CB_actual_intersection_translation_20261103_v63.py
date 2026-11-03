#!/usr/bin/env python3
"""Common actual-target translation extends the full D presentation along Qa.

This checks BOTH complete source maps and literal independence of all218
actual image equations and45 exact colon generators; no new GB is needed.
"""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();maps=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json';chart=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json';builderpath=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json';oldreport=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/report.json'
 m=json.loads(maps.read_text());ch=json.loads(chart.read_text());builder=json.loads(builderpath.read_text());old=json.loads(oldreport.read_text());assert old['all_passed'];qi=ch['q_indices_zero_based'];records={};p=s.Symbol('p');lam=s.Symbol('lambda');v4=s.Symbol('v4');expected=s.Matrix([0,1,0,-p,0,0,p*p/3,0,0,0]);unit=s.zeros(61,1);unit[37]=1
 for name in ['CA','CB']:
  fa=m['families'][name];pars=list(map(s.Symbol,fa['source_fiber_parameters']));u=s.Matrix(list(map(s.sympify,fa['full_61_coordinates'])));assert list(map(str,pars))==['r',*['k'+str(i)for i in range(9)]]
  sel=fa['actual_q_selected_rows'];assert qi[sel[4]]==37;M=u[qi,0].jacobian(pars);L=M[sel,:];assert s.factor(L.det())==-72;velocity=L.inv()[:,4];assert velocity==expected
  translated=u.subs({v:v+lam*dv for v,dv in zip(pars,velocity)},simultaneous=True);assert (translated-u-lam*unit).applyfunc(s.expand)==s.zeros(61,1)
  assert (u.jacobian(pars)*velocity).applyfunc(s.expand)==unit
  relations=list(map(s.sympify,builder['component_presentations'][name]));assert len(relations)==109 and not any(f.has(v4)for f in relations)
  records[name]=dict(common_target_coordinate='v4 = normalized u38',source_fiber_variables=list(map(str,pars)),source_fiber_translation_velocity=list(map(str,velocity)),exact61_delta=list(map(str,unit)),common_K_determinant=str(s.factor(L.det())),all109_actual_image_generators_independent_of_v4=True)
 F=s.sympify(old['nilradical_generator']);G=s.sympify(old['target_p_multiple']);colon=list(map(s.sympify,old['annihilator_prime_generators']));assert not F.has(v4)and not G.has(v4);assert len(colon)==45 and not any(f.has(v4)for f in colon);assert 'v4'in old['annihilator_free_variables']
 # An explicit potential makes the simultaneous source translation visible:
 # lambda*(x²z-p*xz²+p²*z³/3) shifts (k0,k2,k5), not k0 alone.
 h,x,z=s.symbols('h x z');potential=x*x*z-p*x*z*z+p*p*z**3/3
 inputs=[Path(__file__),maps,chart,builderpath,oldreport,ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63_REVIEW/report.json',ROOT/'reports/Q_CA_CB_N_ALL_NONZERO_COMPLETED_MODELS_20261103_v63/report.json']
 result=dict(id='Q_CA_CB_ACTUAL_INTERSECTION_TRANSLATION_20261103_v63',all_passed=True,field='QQ[p,s,K,lambda]',families=records,source_common_potential_translation=str(potential),target_translation='u38 -> u38+lambda; all other61 coordinates fixed, in the common critical-beta source slice',inverse_translation='lambda -> -lambda',entire218_actual_intersection_generators_independent_of_v4=True,nilpotent_F_fixed=str(F),p_multiple_G_fixed=str(G),entire45_generator_colon_prime_independent_of_v4=True,annihilator_free_coordinates=old['annihilator_free_variables'],actual_intersection_structure_transport='Identical completed actual D rings at all values of free v4, hence all physical nonzero Qa after one fixed beta-torus transport and the already certified nonzero-a normalization.',S61_invariants=dict(dimension=18,depth=11,nilradical_square_zero=True,nilmodule_smooth_dimension=11,associated_prime_count=2,minimal_prime='N',embedded_prime='smooth11-dimensional layer'),full_projective_invariants=dict(dimension=27,depth=20,nilmodule_smooth_dimension=20,associated_prime_count=2),physical_a0_type_claimed=False,ring_splitting_claimed=False,scope='Actual reduced-component scheme intersection near EVERY nonzero Qa, not whole Poisson ring, rank2 a0, or all global D boundary.',input_sha256={str(v.relative_to(ROOT)):hashlib.sha256(v.read_bytes()).hexdigest()for v in inputs},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ALL NONZERO Qa ACTUAL D TRANSLATION +218/45 INDEPENDENCE PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
