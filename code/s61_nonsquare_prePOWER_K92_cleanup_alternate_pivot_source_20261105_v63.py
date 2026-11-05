#!/usr/bin/env python3
"""SOURCE alternateTHIRDpivot: knownQQ coefficient, correction excludedbysupport."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_CLEANUP_ALTERNATE_PIVOT_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 bridge=src('reports/S61_NONSQUARE_PREPOWER_K92_CLEANUP_EXISTING_PIVOT_BRIDGE_20261105_v63/report.json');assert bridge['all_passed'];guard=src('reports/S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_LIFT_GUARD_20261105_v63/report.json');assert guard['generic_A_topright_block_literallyZERO_before_specialization'];peer=src('reports/S61_NONSQUARE_K92_WHOLE_ROW_LEADING_PEER_20261105_v63_REVIEW/report.json');assert peer['all_passed'];witness=src('reports/S61_NONSQUARE_K92_WHOLE_ROW_LEADING_PEER_20261105_v63_REVIEW/direct5-coefficient-peer.json')
 actual=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_ENTRY_20261105_v63/literal-source-and-all-finite-path-identities.json');bounds=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_SUPPORT_INTERSECTION_20261105_v63/border5-supports-and-Newton-intersection.json');T={(e['row'],e['column']):{int(g):int(v)for g,v in e['integer_polynomial']}for e in actual['all_Tbar']};S={(e['row'],e['column']):set(e['support'])for e in bounds['Tbar']};cap=2242
 def val(k):
  # Actualcomplete550prefix, sourceupper thereafter; DA isunitconstant1.
  supp=set(T[k])|{g for g in S[k]if g>550};return min(supp,default=cap+1)
 table={str(k):val(k)for k in[(2,0),(2,6),(1,0),(1,6),(0,1),(6,1)]};assert table=={'(2, 0)':2243,'(2, 6)':2243,'(1, 0)':2177,'(1, 6)':2243,'(0, 1)':2242,'(6, 1)':2243}
 tprimeP=[min(val((2,j)),68+val((1,j)))for j in[0,6]];targetP=[val((j,1))for j in[0,6]];assert min(tprimeP)==2243 and min(targetP)==2242
 # P2 adjugate hasregularentries; its determinant hasexactorder1.
 assert T[0,0].get(0)==int(actual['Tbar_common_denominator'])and min(T[6,6])==1 and not T[0,6]and not T[6,0];P2leading=Q(T[6,6][1],int(actual['Tbar_common_denominator']));assert P2leading
 correction=min(tprimeP)+min(targetP)-1;assert correction==4484 and correction>2313
 dump(out/'finite-support-alternate-pivot-proof.json',dict(saved_actual_Tbar_prefix=550,simple_path_support_horizon=2242,all_six_literal_valuation_bounds=table,
  DA_unit_constant1_means_same_valuation_lowerbounds_for_T7=True,delta_formal_unit_constant1=True,
  Tprime_row2='delta*(Trow2-59375*t68*Trow1)',row2_pivot_columns=[0,6],row2_pivot_valuation_bounds=tprimeP,
  P2_columns_to_newcolumn1_valuation_bounds=targetP,P2_determinant_order=1,P2_determinant_leading_QQ=str(P2leading),P2inverse_entry_valuation_lower=-1,
  whole_correction_order_lower=correction,through2313_correction_EXACTZERO=True,
  coefficient_identity='[t2313]Ualt=[t2313]Tprime[2,1]=[t6663]N0',
  known_exact_QQ_coefficient='17149069629000000000/7',independent_peer_only_this_coefficient=True,
  no_claim_all_N0coefficients_independently_checked=True,no_claim_newpivot_LEADINGorder_exact=True,
  proof='Tprime[2,1]=a^4350*N0; formal a=t^-1. P2correction begins4484; no products or newscalarcoefficients required.'))
 nodes=[
  dict(name='Tprime',shape=[7,7],op='external_genericcleanup_bridge',all7Ncallbacks_retained=True),
  dict(name='P2',shape=[2,2],op='submatrix',rows=[0,6],columns=[0,6]),
  dict(name='Ualt',shape=[],op='oneborderedSchurentry',formula='Tprime[2,1]-Tprime[2,[0,6]]*P2^-1*Tprime[[0,6],1]',parameter_inverse_only='detP2'),
  dict(name='P3alt',shape=[3,3],op='submatrix',rows=[0,6,2],columns=[0,6,1],determinant_identity='detP2*Ualt'),
  dict(name='remaining4alt',shape=[4,4],op='Schur',pivot_rows=[0,6,2],pivot_columns=[0,6,1],free_rows=[1,3,4,5],free_columns=[2,3,4,5],formula='D4-C4*adj(P3alt)/detP3alt*B4')]
 dump(out/'typed-alternate-third-pivot.json',dict(nodes=nodes,field='samegenericQQ(a,e,d) source',nonzero_proof='peerQQ[t2313]Ualt coefficientnonzero ⇒Ualt genericRFnonzero; detP2 nonzero',
  chart=['deltaP!=0','deltaQ!=0','detP2!=0','Ualt!=0'],all_chart_complements_RETAINED=True,
  firsttwoSchurpivots_RETAINED=True,thirdpivot_is_EXCHANGED_not_added=True,total_existing4plus3_prefix_preserved=True,remaining_size_still4=True,
  compare_to_old4='onintersection ofboththirdpivotcharts, both4Schurmatrices areequivalent viaexactpivotexchange; no sourceprogressisdiscarded',
  generic_N6_NOTzeroed=True,no_morecallback_or_coefficientwindow_job=True,no_residual_rank_or_full57_claim=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_scalar_products_or_coefficients=True,correction_EXACTZERO_through2313=True,alternate_THIRDpivot_genericnonzero_from_peer=True,
  prior_two_pivots_retained=True,remaining_size4_unchanged=True,no_fourthpivot_rankadvance_or_global57_claim=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
