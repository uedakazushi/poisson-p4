#!/usr/bin/env python3
"""Assemble the height-one transverse-bound certificates for the normalized binary H8."""
import argparse,hashlib,json,time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def certificate():
 st=time.monotonic();hs={};raw={};data={};names=[]
 def load(short,stem):
  name=f'reports/{stem}.json';p=ROOT/name;j=json.loads(p.read_text());assert j.get('all_passed',True)
  dep=j.get('input_sha256',j.get('source_sha256',{}))
  for f,h in dep.items():
   assert not Path(f).is_absolute() and digest(ROOT/f)==h,f
   assert f not in raw or raw[f]==h,f
   raw[f]=h
  hs[name]=digest(p);names.append(name);data[short]=j;return j
 for short,prefix in [('horizontal','BFREE_HORIZONTAL'),('eo','E_ONLY_SUPPORT'),('vert','VERTICAL_SLICE_CERTIFICATE'),('original','ORIGINAL_ROOT_QUARTICS'),('alt','C_TRANSVERSE_LENGTH'),('axes','AXES_TRANSVERSE_LENGTH'),('quad','QUADRATIC_CODIMENSION'),('cubic','CUBIC_CODIMENSION'),('heavy','HEAVY_SQUARE_PIVOT'),('frame','BFREE_FRAME'),('generic','GENERIC_RANKS'),('old','WEIGHTED_CHART')]:
  load(short,'S61_BINARY_HYPERPLANE_'+prefix+'_20261102')
 load('action','S61_BINARY_NUISANCE_CONTRACTION_20261102');load('ordinary','S61_BINARY_FORMAL_RING_20261101')
 h=data['horizontal'];assert h['horizontal_row_module_surjective'] and h['weighted_algebra_free_rank_over_Kd']==65 and h['weighted_Hilbert_function']==[1,7,23,34]
 ee=data['eo'];assert (ee['original_degree'],ee['leading_coefficient_mod'],ee['modular_gcd_degree'])==(942,2514,0)
 ve=data['vert'];assert ve['selected_minor_exact_a_degree_QQ_e']==575 and ve['direct_modular_minor']['leading_coefficient']==16008
 mo=ve['full_modular_row_modules']['3'];assert mo['rank']==238 and mo['degree']==287 and mo['known_orders']=={'a':26,'F':64,'D0':66,'C':131}
 assert data['quad']['bad_locus_dimension_upper_on_three_core_open']==1
 assert data['cubic']['cubic_rank_drop_dimension_upper']==1 and data['cubic']['generic_cubic_rank']==57
 org=data['original'];assert org['cuts']==[5,19,23]and org['original_equations_reconstructed']==135
 p=32003;factors={'a':(1,0),'F':(500,342),'D0':(50,36),'C':(100,52)}
 roots={k:(-c*pow(b,-1,p))%p for k,(b,c)in factors.items()};assert len(set(roots.values()))==4
 certs={}
 for k,v in [('a',924),('F',10766),('D0',28884)]:
  w=org['points'][k];assert w['a_modp']==roots[k]and w['minor_det_modp']==v and w['rank_134']==2814;certs[k]=v
 ct=data['alt'];assert ct['cuts']==[20,19,23]and Fraction(ct['new_cut_determinant_QQ'])==Fraction(79,25)and ct['weighted_Hilbert_function_QQ']==[1,7,23,34,0]and ct['length_QQ']==65
 q=Fraction(-13,25);assert q.numerator*pow(q.denominator,-1,p)%p==roots['C'];assert ct['all_weights_at_least_four_zero']
 # The two axis charts specialize the raw equations before changing the number of normals.
 axis=data['axes']
 for k in ('a','F'):
  a=axis['points'][k];assert a['length']==65 and a['weighted_Hilbert_function_QQ']==a['weighted_Hilbert_function_modp']==[1,7,23,34,0] and a['all_weights_at_least_four_zero']
  assert a['cuts']==org['cuts'] and a['transverse_determinant_modp']%p and a['weight4_minor_modp']==certs[k]
 assert all(x>0 for x in data['action']['nuisance_weights']) and data['action']['core_weights']==[0,0,0]
 assert data['ordinary']['uniform_normal_length_lower_bound']==65 and data['ordinary']['homogenizer_square_zero']
 extras=['research/notes/S61_BINARY_DEEPER_FULL_CORE_REMOVAL_20261102.md','research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md','research/notes/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.md']
 for x in extras:hs[x]=digest(ROOT/x)
 src=str(Path(__file__).relative_to(ROOT));note=f'research/notes/{STEM}.md';hs[src]=digest(ROOT/src);hs[note]=digest(ROOT/note)
 return dict(id=STEM,all_passed=True,scope='For the normalized H8 presentation there is a rationally defined open in all three core coordinates with complementary dimension at most one, with an actual transverse length bound 65 using one of finitely many cuts. All five nuisance parameters there have relative t squared zero; the corresponding ambient conclusion holds over sufficiently general beta. The remaining points over general beta are at most finite, not proved absent.',core_dimension=3,exceptional_core_dimension_upper=1,explicit_exceptional_curve_equations=False,checked_cases=['positive d degree','e only','positive a degree independent of d','four known a-positive divisors'],special_prime=p,monic_a_degree=575,monic_a_leading_residue=16008,special_a_roots=roots,old_cut_quartic_residues=certs,alternate_C_cut=ct['cuts'],alternate_C_cut_determinant=ct['new_cut_determinant_QQ'],coordinate_free_existence_of_transverse=True,old_cut_full_rank_everywhere_claimed=False,e_only_degree=942,e_only_LC_mod=2514,e_only_gcd_mod=0,relative_completion_base_dimension=11,relative_module_rank_on_zero_nuisance_section=65,all_nuisance_t_square_zero_on_open=True,residual_source_dimension_bound=1+5+10,residual_general_fiber_dimension_upper=0,residual_general_fiber_empty_claimed=False,full_S61_exhaustion_claimed=False,primary_inputs=names,primary_dependency_hash_count=len(raw),input_sha256=hs,primary_dependency_sha256=raw,seconds=time.monotonic()-st)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert not a.output.exists()
 r=certificate()
 with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
 print(f"PASS {len(r['primary_inputs'])} primary reports, {r['primary_dependency_hash_count']} dependencies; {r['seconds']:.3f} s")
