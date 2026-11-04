#!/usr/bin/env python3
"""One seven-coefficient symbolic test of the projected15 cubic criterion."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 sp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';src=read(sp);mons=src['normal_monomials'];light=[0,1,2,3,4,6,7];hcols=[]
 for j in light:
  m=[0]*8;m[5]=1;m[j]=1;hcols.append(mons.index(m))
 raw=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt';bp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';basis=read(bp);assert sha(raw)==basis['remaining_matrix_sha256']
 assert basis['remaining_original_row_indices'][42]==43 and basis['remaining_original_row_indices'][47]==49 and src['labels'][43]==['cubic',20]and src['labels'][49]==['cubic',61]
 R,a,e,d=ring('a,e,d',s.QQ);qq=lambda terms:R.from_dict({tuple(m):s.QQ(v.numerator,v.denominator)for m,v in terms.items()})
 saved={};
 with raw.open()as f:
  assert f.readline().split()==['57','91','57']
  for i in range(57):
   for j in range(91):
    line=f.readline();assert line
    if i not in[42,47]or j not in hcols:continue
    t=line.split();n=int(t[0]);assert len(t)==1+4*n;poly={}
    for k in range(n):
     co=Q(t[1+4*k]);ex=tuple(map(int,t[2+4*k:5+4*k]));poly[ex]=poly.get(ex,Q())+co
    saved[i,j]=qq(poly)
  assert not f.read().strip()
 diff={j:saved[47,col]-saved[42,col]for j,col in zip(light,hcols)}
 oldfive=[diff[0]-3*diff[6],diff[1]+diff[6],diff[2],diff[4],diff[7]];assert all(not v for v in oldfive)
 qp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';qd=read(qp)
 assert all(not m[5]or not any(Q(v)for _,v in terms)for row in qd['quadratic_rows_sparse']for m,terms in row)
 C=100*a+18*e**2-9*e-2;F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;J=5*a*d-6*e-4;delta=-64000*C**2;known=[a,d,F,B,C,J]
 def factorunit(p):
  if not p:return None
  powers=[]
  for f in known:
   n=0
   while not p.rem(f):p=p.exquo(f);n+=1
   powers.append(n)
  if not p.is_ground:return None
  return{'constant':str(p),'a_d_F_B_C_J_powers':powers}
 enc=lambda p:[[list(m),str(v)]for m,v in p.items()]
 values={'Z':diff[6],'W':diff[3]};unitproof={j:factorunit(delta*v)for j,v in values.items()};allzero=all(not v for v in values.values())
 deps=[ROOT/'reports'/p/'report.json'for p in['S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW','S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW','S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW']];assert all(read(p)['all_passed']for p in deps)
 proof={'actual_source_rows':{'remaining57':[42,47],'normalized63':[43,49],'actual135':[20,61]},'physical_H_times_light_columns':dict(zip(map(str,light),hcols)),
 'unscaled_G61_minus_G20_heavy_coefficients':{str(j):enc(v)for j,v in diff.items()},
 'oldfive_after_ZW_substitution_all_zero':True,'remaining_heavy_form':'kZ*Z+kW*W',
 'known_unit_delta':enc(delta),'delta_scaled_HZ_HW_coefficients':{j:enc(delta*v)for j,v in values.items()},'known_unit_proofs':unitproof,
 'all_H_coefficients_zero':allzero,'atleastone_known_FULLU_unit_coefficient':any(v is not None for v in unitproof.values()),
 'projected15_criterion':{'domain':'FULLU D(Delta34), whereDelta34 is the accepted saved nonzero34minor of restricted36x40. NoU78/CI assumption needed for this degree3 quotient.',
 'exact_sequence':'0->T29->B3_weighted50->V21->0, from actual product41 and actualZW product19',
 'source_kernel':'G61-G20 primitive rank1 kernelofprojectedgenuine16 onDelta34; compactrelationdelta*diff=8100q1E+ZA+WB provesrestrictionzero',
 'full57_equivalence':'OnDelta34, actualfullcubicrank57 iffclassdelta*(G61-G20) is nonzeroineachB3fiber; a unit HZ/HW coefficient would give uniformpositivegate sincequadproductshaveNOH'},
 'no_new_point_prime_rank_GB_field_solver_quartic_expansion_or_source26':True,'no_allFULLU57_or_full239_claim':True}
 (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 result={'id':ID,'all_source_and_polynomial_checks_passed':True,'all7_heavy_coefficient_differences_extracted_exactly':True,'heavy_coefficient_unit_gate_succeeds':proof['atleastone_known_FULLU_unit_coefficient'],'all_heavy_coefficients_zero':allzero,
 'no_new_point_or_rank_job':True,'no_generic_oldpoint_repetition_or_full57_allFULLU_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),sp,bp,raw,qp,*deps]},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-began,'frozen':False}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
