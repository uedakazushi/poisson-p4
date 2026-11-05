#!/usr/bin/env python3
"""Prior fixed-membership product margin; no coefficient construction."""
import argparse,json,math,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_VERIFY_HEIGHT_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 p='reports/AMBIENT_B18_ELL2_T2_SOURCE8_RESOLUTION_PRIOR_HEIGHT_20261105_v63';prior=read(p+'/report.json');d=read(p+'/prior-input-receipt-product-and-partialsum-height-envelopes.json');t=read(H.PATHS['transport']);classes={}
 for kind,key in [('G2','active_F1_fixed67_columns'),('G3','active_F2_fixed182_columns')]:
  values=[F(x) for col in t[key] for i,x in col];D=math.lcm(*(x.denominator for x in values));U=max(abs(x.numerator)*(D//x.denominator) for x in values);classes[kind]=[D,U]
 membership=[];maximum=prior['factor_construction_prior_height_bound_bits']
 for r in d['all56_span_coordinate_positive_L1_bounds']:
  kind='G2' if r['kind'].startswith('G2') else 'G3';D,U=map(int,r['D_U']);Dm,Um=classes[kind];b=max((D*Dm).bit_length(),(U*Um).bit_length());membership.append(b);maximum=max(maximum,b)
 for kind,k in [('G2','repr67_L1_height_bound_bits'),('G3','repr182_L1_height_bound_bits')]:
  # Each representation has a single commonD/L1U of at most2^b.
  maximum=max(maximum,d[k]+max(x.bit_length() for x in classes[kind]))
 assert maximum<=384
 report={'id':STEM,'all_passed':True,'prior_factor_height_bits':prior['factor_construction_prior_height_bound_bits'],
  'fixed_membership_partialsum_height_bits':max(membership),'all_factor_and_membership_commit_height_bits':maximum,
  'uncommitted_scalar_product_sum_height_bits':3*maximum+1,'prepared384guard_covered':True,
  'fixed_matrix_common_D_U':{k:list(map(str,v)) for k,v in classes.items()},
  'proof':'ForvectorcoefficientswithcommonD andL1clearedboundU, eachpartialmatrixproductsum hasdenominator dividingD*Dm andabsnumerator<=(U*Um). No termcount neededbecauseU boundsSUMabsolutevectorcoefficients. Representation-membership usesfrozenreceiptL1bitboundsandthesameclass. No signednewfactororresponsecomputed.',
  'new_factor_response_or_physical_values':False,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','all_factor_and_membership_commit_height_bits','uncommitted_scalar_product_sum_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
