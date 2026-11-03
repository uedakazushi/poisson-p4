#!/usr/bin/env python3
"""Independent full-QQ proof of the new four-constant J0 unit identity."""
import argparse,gzip,hashlib,json,shutil,time
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);t=time.monotonic()
 source=ROOT/'work/nonsquare-J0-pid-v63/J0-residual.txt';work=ROOT/'work/nonsquare-J0-pid-v63/primitive-span-v1'
 first=ROOT/'reports/S61_NONSQUARE_J0_PID_20261103_v63/report.json';frozen=json.loads(first.read_text());assert sha(source)==frozen['omitted_easy_to_regenerate_outputs'][str(source.relative_to(ROOT))]['sha256']
 search=json.loads((work/'report.json').read_text());assert search['new_original_row_unit_witness_verified']and search['new_witness_row']==7
 assert search['pivot_replay']['status']=='completed'and search['pivot_replay']['residual_shape']==[22,56]and search['pivot_replay']['additional_known_unit_pivots']==1
 R,a,e=ring('a,e',sp.QQ);selected={};raw={}
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[23,57,13];f.readline();f.readline()
  for i in range(23):
   for j in range(57):
    z=f.readline().split()
    if i==7 and j in [4,15,25,30]:
     n=int(z[0]);assert len(z)==1+4*n and all(int(z[4+4*k])==0 for k in range(n))
     selected[j]=R.from_dict({tuple(map(int,z[2+4*k:4+4*k])):sp.QQ(z[1+4*k])for k in range(n)})
     raw[str(j)]=[[list(ex),str(c)]for ex,c in sorted(selected[j].items())]
 delta=selected[4]-selected[15]-selected[25]+selected[30]
 C=100*a+18*e**2-9*e-2;unit=a*C*(3*e+2)**3
 c=delta.exquo(unit);assert c and c.is_ground and delta==c*unit
 coefficients=[]
 with (work/'row7.witness').open()as f:
  assert int(f.readline())==57
  for j in range(57):
   z=f.readline().split();n=int(z[0])
   if n:assert n==1 and z[2:]==['0','0','0'];coefficients.append((j,sp.QQ(z[1])))
 assert [j for j,v in coefficients]==[4,15,25,30]
 assert [v for j,v in coefficients]==[1/c,-1/c,-1/c,1/c]
 result=work/'residual.txt';assert sha(result)==search['pivot_replay']['residual_sha256']
 payload={'source_shape':[23,57],'target_row':7,'columns':[4,15,25,30],'constant_coefficients':[1,-1,-1,1],
  'source_four_polynomials':raw,'nonzero_QQ_scalar':str(c),'target_factorization':'a*C*(3e+2)^3','C':'100a+18e^2-9e-2',
  'full_delta_terms':[[list(ex),str(v)]for ex,v in sorted(delta.items())]}
 cp=args.output/'identity.json.gz'
 with gzip.open(cp,'wt')as f:json.dump(payload,f,separators=(',',':'))
 for name in ['primitive-search.cpp','row7.witness','row7.log','pivot.log']:shutil.copyfile(work/name,args.output/name)
 shutil.copyfile(ROOT/'work/nonsquare-J0-pid-v63/units.txt',args.output/'units.txt')
 shutil.copyfile(work/'report.json',args.output/'search-and-pivot.json')
 inputs=[Path(__file__),ROOT/'code/s61_nonsquare_J0_primitive_span_20261103_v63.py',source,first,work/'report.json',work/'primitive-search.cpp',ROOT/'code/polynomial_unit_span_pivots_20261103_v62.cpp']
 report={'id':'S61_NONSQUARE_J0_PRIMITIVE_UNIT_20261103_V63','all_passed':True,'coefficient_field':'QQ',
  'base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]','J_not_inverted':True,'unknown_core_factor_inverted':False,
  'independent_full_polynomial_identity_verified':True,'target_row':7,'columns':[4,15,25,30],'coefficients':[1,-1,-1,1],
  'nonzero_QQ_scalar':str(c),'unit_factorization':'a*C*(3e+2)^3','unit_factor_count':5,
  'source_shape':[23,57],'new_residual_shape':[22,56],'additional_known_unit_pivots':1,'total_safe_pivots_from57basis':35,
  'constant_column_replacement_unimodular':True,'original_and_new_cokernels_isomorphic_on_stated_locus':True,
  'full_coker_unit_or_rank57_allpoints_proved':False,'J0_defect_support_finite_proved':False,
  'search_and_exact_pivot_seconds':search['seconds'],'independent_identity_seconds':time.monotonic()-t,
  'new_residual_sha256':sha(result),'new_residual_bytes':result.stat().st_size,
  'new_residual_omitted_as_reproducible':'Regenerate from archived J0 input and saved witness; approximately5seconds for this pivot after the prior10second13-pivot replay.',
  'identity_sha256':sha(cp),'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs}}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','unit_factorization','new_residual_shape','independent_identity_seconds']},indent=2))
if __name__=='__main__':main()
