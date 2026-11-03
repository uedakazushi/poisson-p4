#!/usr/bin/env python3
"""Finite15-chart scalar criterion retaining the full generic affine residual."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
 fit=ROOT/'reports/S61_AFFINE_NORMAL_SCALAR_FITTING_20261103_v63/report.json';fr=json.loads(fit.read_text());assert fr['all_passed']
 for rel,want in fr['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in fr['evidence_sha256'].items():assert sha(fit.parent/rel)==want,rel
 finite=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/report.json';assert json.loads(finite.read_text())['all_passed']
 projection=ROOT/'reports/S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102.json';pr=json.loads(projection.read_text());assert pr['all_passed']and pr['generic_affine_projection_is_finite']and pr['dimension_upper']==2
 rowsets=[(1,2,3,4)];columnsets=list(itertools.combinations(range(6),4));cover=[dict(rows=list(r),columns=list(c))for r in rowsets for c in columnsets];assert len(cover)==15
 recipe=dict(base='K=C(beta16),intersectestablishednonemptybetaopens;Rprime=K[n42,S1,S2,H]',faithful_constant_extension='L=K(lambda_i_j:i0..669,j0..5),4020algebraicallyindependentconstants',matrix='Phi5x670fromglobalquintic;Psi=Phi*Lambda5x6',rank_on_unknown_B=4,dimension_B_atmost=2,
  bad_incidence=dict(parameter_dimension=4020,fiber_bad_codimension=(4-3)*(6-3),base_dimension=2,total_bad_dimension_atmost=4019),
  cover=cover,cover_basis='All15fourbyfourminorsdeltaofPsi_fixedbottom4rows',
  scalar_equations='g_j=det[Psi_(all5rows,Cfourcolumns),Phi_(all5rows,j)] forj0..669',
  actual_chart_identity='B_L localizeddelta = Spec Rprime_L[delta^-1]/(g_0,...,g_669)',scalar_generators_atmost=670,
  parameter_degree_bounds=dict(delta=8,g_j=10),known_origin_excluded_by_all_deltas=True,
  exact_exhaustion_criterion='B=empty iff all15scalarlocalizedidealsareunitsoverL',
  no_unknown_det_zero_complement_discarded=True,all_unknown_nilpotents_retained=True,
  not_claimed=['specificQQconstantLambdaverified','anychartunitcomputed','beta-onlyfiniteness','actualgenericfiveexhaustion'])
 (out/'cover-recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
 paths=[Path(__file__),fit,finite,projection,ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz',ROOT/'research/notes/S61_AFFINE_NORMAL_SCHUR_COVER_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_SCHUR_COVER_20261103_v63',all_passed=True,proof_only_finite_cover_receipt=True,charts=15,bad_incidence_codimension3_exceeds_base_dimension2=True,full_generic_unknown_residual_scalarized=True,actual_units_not_computed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'cover-recipe.json':sha(out/'cover-recipe.json')})
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('EXACT15SCHURCHART COVER/ALLUNKNOWN SUPPORT RETAINED;UNITSUNCOMPUTED')
if __name__=='__main__':main()
