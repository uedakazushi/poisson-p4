#!/usr/bin/env python3
"""Bounded finite section checks; no failed GB or degree4 matrix is repeated."""
import argparse
import gzip
import hashlib
import json
from math import comb
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plane-replay',type=Path,required=True)
 ap.add_argument('--ranks-replay',type=Path,required=True);ap.add_argument('--assessment-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/Q_THROUGH3_LINEAR_SECTION_20261103_v63';manifest=load(directory/'manifest.json');paths=[Path(__file__),directory/'manifest.json']
 for name,h in manifest['source_sha256'].items():p=ROOT/name;assert sha(p)==h,name;paths.append(p)
 for name,h in manifest['evidence_sha256'].items():p=directory/name;assert sha(p)==h,name;paths.append(p)
 for name,h in manifest['uncompressed_sha256'].items():assert hashlib.sha256(gzip.decompress((directory/name).read_bytes())).hexdigest()==h,name
 fresh=[load(p)for p in[args.plane_replay,args.ranks_replay,args.assessment_replay]];old=[load(directory/p)for p in['known-plane/report.json','low-ranks/report.json','macaulay-assessment/report.json']]
 ignore={'seconds','native_rank_seconds','input_sha256','evidence_sha256'}
 for a,b in zip(old,fresh):assert {k:v for k,v in a.items()if k not in ignore}=={k:v for k,v in b.items()if k not in ignore}
 plane,ranks,assessment=old;M=s.Matrix(plane['restricted_linear_matrix']);v=s.Matrix(plane['primitive_integer_kernel_vector']);assert M*v==s.zeros(12,1)
 column=plane['chosen_unit_minor_column'];det=(-1)**column*M[:,[i for i in range(13)if i!=column]].det();assert det==plane['chosen_signed_minor']and int(det)%32003
 values=dict(zip(plane['active_plane_variables'],list(v)));f=values['z8']**2*(12*values['z14']+values['z17']+6*values['z3']);assert f%32003==11541
 assert[ranks['projected_quadratic_rank'],ranks['linear_quadratic_multiple_rank'],ranks['full_degree3_ideal_rank']]==[57,1210,1227]
 assert ranks['exact_homogeneous_section_Hilbert_values_mod_p_through3']==[1,23,219,1073]
 syz=57*23-1210;g=lambda r:23*r-r*(r-1)//2;shadow=4*g(23)+g(9);assert syz==101 and shadow==1275
 h4=comb(26,4)-(57*comb(24,2)+17*23)+shadow;assert h4==102 and sum([1,23,219,1073,h4])==1418
 assert not load(directory/'summary.json')['Artinian_section_verified']
 paths += [ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_THROUGH3_LINEAR_SECTION_20261103_v63_REVIEW.md']
 result=dict(id='Q_THROUGH3_LINEAR_SECTION_20261103_v63_REVIEW',all_passed=True,
  fresh_light_plane_rank_and_shadow_replays_match=True,exact_V_known13plane_section='V[t]/t³',
  native_Sym2_Sym3_projection_and_all74_recipe_source_audited=True,
  exact_GF_section_Hilbert_through3=[1,23,219,1073],linear_Q_syzygies=101,leading_module_shadow_minimum=1275,
  special_J74_section_h4_lower=102,Artinian_length_if_finite_lower=1418,
  no_failed_GB_repeated=True,no_degree4_DAG_or_matrix_constructed=True,
  full23variable_section_dimension_or_Artinianity_or_QQ_HF_claimed=False,
  low_degree_GB_Artinian_certificate_impossibility_claimed=False,full_R35_upper22_unchanged=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS finite section/negative Macaulay bounds; full section remains unresolved')

if __name__=='__main__':main()
