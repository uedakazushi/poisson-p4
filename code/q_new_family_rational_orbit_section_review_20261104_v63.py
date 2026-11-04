#!/usr/bin/env python3
"""Independent finite tensor torus weights and cross-multiplied quadratic folding."""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];STEM='Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
 d=read(p);assert d['all_passed']
 for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
 for n,h in d.get('evidence_sha256',{}).items():
  q=ROOT/n if (ROOT/n).exists() else p.parent/n
  assert sha(q)==h,n
 return d
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 assert not out.exists();out.mkdir(parents=True)
 base=ROOT/'reports/Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63';ap=base/'report.json';accepted(ap);cp=base/'certificate.json';cert=read(cp)
 src=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json';source=read(src)
 repath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json';rep=read(repath)
 refs=[ROOT/f'reports/{n}/report.json' for n in ['Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW',
  'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW','Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW','Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63_REVIEW']]
 for p in refs:accepted(p)
 H=[1,-1,0,0];weights=[w[0] for w in rep['Cartan_weights']]
 basis=[[Q(v) for v in row] for row in rep['primitive24_by16_basis']]
 coords=[(i,j,k) for i,j in itertools.combinations(range(4),2) for k in range(4)]
 for r,(i,j,k) in enumerate(coords):
  for c in range(16):assert (H[i]+H[j]+H[k])*basis[r][c]==weights[c]*basis[r][c]
 rho=[[Q(v) for v in row] for row in rep['all_ten_representation_matrices'][1]]
 assert rho==[[weights[i] if i==j else 0 for j in range(16)] for i in range(16)]
 J=[[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]
 assert all(not J[i][j] or H[i]+H[j]==0 for i in range(4) for j in range(4))
 xi,r,t,q,B,u,v,c,X,W=s.symbols('xi r s q B u v c X W');h=t*t-25*r*q*q
 R,*_=ring('xi,r,s,q,B,u',s.QQ);F=R.to_field()
 substitutions={v:-32*u*t/h,c:xi*B*u/40,X:t/(48*u),W:(q+1)/3072}
 def fold(p):
  parts=[{},{}]
  for ex,cof in p.items():
   k=ex[5];e=list(ex[:5]);e[1]+=k//2;e=tuple(e)
   parts[k%2][e]=parts[k%2].get(e,Q(0))+Q(str(cof))*1600**(k//2)
  return [{e:c for e,c in a.items() if c} for a in parts]
 def zero_in_extension(expression):
  z=F.from_expr(expression);assert not any(fold(z.numer));assert any(fold(z.denom));return z
 results=[];nonzero=0
 for i,expr in enumerate(source['primitive_beta16']):
  bar=s.sympify(cert['beta_bar16'][i]);assert bar.free_symbols <= {xi,r,t,q,B}
  difference=s.sympify(expr).subs(substitutions,simultaneous=True)*u**weights[i]-bar
  z=zero_in_extension(difference);value=F.from_expr(bar);nonzero+=bool(value)
  results.append(dict(index=i,Cartan_weight=weights[i],zero_even_and_odd_numerator_remainders=True,
   denominator_nonzero_in_extension=True,beta_bar_numerator_terms=len(value.numer),beta_bar_denominator_terms=len(value.denom)))
 assert nonzero==10
 oldv=-98304*X/(147456*X**2-(3072*W-1)**2)
 zero_in_extension(oldv.subs({X:t/(48*u),W:(q+1)/3072},simultaneous=True)-substitutions[v])
 # Explicit transversality of the ACTUAL t=2 source section to the acceptedflow.
 flowp=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json';flow=read(flowp)
 rf,sf=s.symbols('r s');rr=s.sympify(flow['normalized_coordinates']['r']);ss=s.sympify(flow['normalized_coordinates']['s'])
 dt=s.sympify(flow['delta_t']).subs({rf:rr,sf:ss},simultaneous=True)
 tt=s.sympify(flow['fibre_parameter_t']).subs({rf:rr,sf:ss},simultaneous=True)
 sec=s.sympify(flow['section_v']);K=s.QQ.frac_field(X,W)
 dtsec=K.from_sympy(s.cancel(dt.subs(v,sec)));assert dtsec
 assert s.cancel(tt.subs(v,sec)-2)==0 and c not in sec.free_symbols
 assert not cert['full_hypersurface_rationality_or_generic_stabilizer_trivial_proved']
 evidence=dict(all16_inverse_source_fold_checks=results,nonzero_beta_coordinates=10,
  finite_actual24_cocharacter_all384_basis_entries=True,rho_generator1_exact_diagonal=True,
  source_t2_v_equals_minus32us_over_h=True,actual_delta_t_on_t2_section=str(dtsec),
  source_section_transverse_and_same_dense_Sp4_saturation_proof_in_note=True,
  index_one_generic_stabilizer_or_full_hypersurface_rationality_proved=False)
 ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ap,cp,src,repath,*refs,flowp,
  ROOT/'code/q_new_family_rational_orbit_section_20261104_v63.py',ROOT/'research/notes/Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63.md']
 report=dict(id=STEM,all_passed=True,independent_cross_multiplied_even_odd_monomial_fold_not_author_norm_pair=True,
  all16_rational_orbit_representative_identities_and_nonzero_denominators=True,
  actual_finiteH1_cocharacter384_basis_and_original_rho256_coefficients=True,
  source_c_independent_t2_section_and_delta_transversality_checked=True,
  same_dense_Sp4_saturation_on_generic_domain_proved=True,
  rational_five_parameter_family_not_one_orbit_per_parameter_or_stabilizer_claim=True,
  h14_direct_xi1_evaluation_recipe_only_no_expansion=True,
  actual_target_index_one_or_full_B15_rationality_or_source26_claim=False,
  no_new_point_invariant_expansion_GB_or_author_rerun=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k] for k in ['all_passed','all16_rational_orbit_representative_identities_and_nonzero_denominators','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
