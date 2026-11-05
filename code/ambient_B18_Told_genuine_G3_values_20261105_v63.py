#!/usr/bin/env python3
"""One five-term cached scalar; no natural-map derivative is recomputed."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOLD_GENUINE_G3_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_TOLD_GENUINE_G3_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def numeric(v):
 if isinstance(v,(int,Q)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 if isinstance(v,(list,tuple)):return sum(map(numeric,v))
 return 0
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-peer-report',required=True);a=ap.parse_args()
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 sp=SOURCE+'/report.json';ep=SOURCE+'/actual-old-to-genuine-edge-one-cache-witness-and-diagram.json';rep=json.loads((ROOT/sp).read_text());src=json.loads((ROOT/ep).read_text());peer=json.loads((ROOT/a.source_peer_report).read_text())
 assert rep['all_passed']and peer['all_passed']and peer['input_sha256'][sp]==sha(ROOT/sp)
 w=src['one_cached_witness'];jp=w['natural_head0_r_first_callback_path'];old=json.loads((ROOT/jp).read_text());assert sha(ROOT/jp)==rep['input_sha256'][jp]
 assert src['all_five_head0_origin_entries_zero']and w['natural_origin_column5'][0]=='0'
 products=[]
 for j in range(4):
  b=Q(w['genuine_edge_b1_b2_b3_b4'][j]);m=Q(w['natural_origin_column5'][j+1]);products.append(b*m)
 xi=Q(w['genuine_edge_rh0_coefficient_xi']);Jr=Q(old[w['natural_head0_r_first_callback_key']]);products.append(xi*Jr)
 value=sum(products,Q(0));positive=bool(value)
 assert all(max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=512 for v in products+[value])
 data={'old_generator_index':4,'old_actual_free56_label':55,'all_four_origin_products_bj_Mj4':[str(v)for v in products[:4]],'cached_Jr':str(Jr),'known_xi':str(xi),'xi_times_Jr':str(products[4]),'actual_old_head4_genuine_G3_image':['0',str(value)],'canonical_head0_beta_not_used':True,'actual_edge_image_rank':1 if positive else None,'witness_zero_is_partial_negative_only':not positive,
  'actual_intersection_Told_Nbar_length':10 if positive else None,'actual_normalized_nil_image_Nbar_dimension':4 if positive else None,'actual_Ibar_to_NEbar_rank':4 if positive else None,'actual_L_DVR_rank':2 if positive else None,'actual_L_DVR_torsion_min_generator_count':3 if positive else None,
  'actual_Q_vs_Lbar':'Q=Tnil/Told is a length5 submodule ofNEbar; Lbar is a length5 quotient ofNEbar when positive. They are different modules/maps and no fullS-splitting or allz-action assertion is made.',
  'diagram_proof_if_positive':'Actualsource triangle hasker(theta:Tnil→NEbar)=Told. Compatibility with genuineNbar saysker(theta|Nbar)=Nbar∩Told. epsilonoldrank1 andToldlength11 giveintersection10; Nbar14givesimage4. Nbar→Ibar surjective, I Dfree7/NE Dfree9, soTor1^D(L,k)=ker(Ibar→NEbar)hasdimension3 andL hasDgenericrank2.',
  'actual_L_torsion_length_exponents_full_S_action_original_H0_and_nilcube_unknown':True,'new_source_gauge_jet_or_first_coefficient_solve':False}
 owned=numeric([rep,src,peer,old,data]);assert owned<15000,owned
 pins={p:sha(ROOT/p)for p in [sp,ep,jp,a.source_peer_report]};ev=out/'actual-five-cached-products-and-normalizer-image-rank-four.json';blob=json.dumps(data,indent=2)+'\n';assert sum((ROOT/p).stat().st_size for p in pins)+len(blob.encode())<2<<20;ev.write_text(blob)
 r={'id':STEM,'all_passed':True,'status':'PROVED_actual_Told_genuine_edge_rank1_and_normalizer_fiber_image4'if positive else'PARTIAL_NEGATIVE_single_cached_witness_zero','actual_edge_rank':1 if positive else None,'actual_L_generic_DVR_rank':2 if positive else None,'actual_L_torsion_mu':3 if positive else None,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'owned_numeric_slots':owned,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));signal.alarm(0)
if __name__=='__main__':main()
