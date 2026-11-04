#!/usr/bin/env python3
"""Source-only cost for56 fixed one-row replacement minors. No dot/minor/gcd."""
import argparse,json,math,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE15_COVER_PREFLIGHT_20261104_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports';rb=base/'S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63';rank=read(rb/'report.json');a=read(rb/'native.json');assert rank['generic_rank_saved_pure8x98']==8
 kb=base/'S61_NONSQUARE_P100_K2_PURE7_20261104_v63';kr=read(kb/'report.json');assert kr['complete'];rp=[kb/f'pure98-row{i}.json'for i in range(7)];data=[read(p)for p in rp];rows=[d['values98']for d in data]
 paths=[Path(__file__),rb/'report.json',rb/'native.json',rb/'source-labels.json',kb/'report.json',kb/'source-binding.json']+rp
 for p in rp:assert sha(p)==kr['seven_row_sha256'][p.name]
 cs=a['selected_columns'];H=a['final_adjugate'];D=a['determinant_nonzero'];factors=[(1,0),(1,20),(1,60),(203,2700),(1,12)];cache={}
 def degree(p):return max((m[0]for m,c in p),default=-1)
 def denpowers(poly):
  key=tuple((m[0],c)for m,c in poly)
  if key in cache:return cache[key]
  r={m[0]:Q(c)for m,c in poly};assert r[max(r)]==1;es=[]
  for aa,bb in factors:
   n=0
   while r and max(r)>=1:
    v=r.copy();q={}
    while v and max(v)>=1:
     j=max(v);c=v[j]/aa;q[j-1]=q.get(j-1,Q())+c
     for k,c0 in[(j,c*aa),(j-1,c*bb)]:
      v[k]=v.get(k,Q())-c0
      if not v[k]:del v[k]
    if v:break
    r=q;n+=1
   es.append(n)
  assert r and max(r)==0;cache[key]=es;return es
 def primitive_height(p):
  co=[Q(c)for m,c in p];den=1
  for c in co:den=math.lcm(den,c.denominator)
  return den,max([abs(c.numerator*(den//c.denominator)).bit_length()for c in co]+[0])
 hp=[[primitive_height(p)for p in row]for row in H];profiles=[];allpairs=allterms=bytes_upper=0
 for light,row in enumerate(rows):
  ex=[denpowers(v['denominator'])for v in row];mx=[max(v[j]for v in ex)for j in range(5)];dd=sum(mx);vv=[]
  for col in cs:
   v=row[col];n=v['numerator'];dn=degree(v['denominator']);nd=degree(n)+dd-dn if n else-1;den,nb=primitive_height(n)
   # LCD/monicDi is integer with this L1-height bound; no polynomial product.
   qb=sum((mx[j]-ex[col][j])*max(1,(abs(factors[j][0])+abs(factors[j][1])).bit_length())+ex[col][j]*max(0,factors[j][0].bit_length())for j in range(5))
   vv.append({'degree':nd,'terms_upper':max(0,nd+1),'common_QQ_denominator':den,'integer_numerator_height_upper':nb+qb,'source_RF_terms':len(n)+len(v['denominator'])})
  candidates=[]
  for replace in range(8):
   active=[k for k in range(8)if vv[k]['degree']>=0 and H[k][replace]]
   deg=max((vv[k]['degree']+degree(H[k][replace])for k in active),default=-1);terms=max(0,deg+1)
   pairs=sum(vv[k]['terms_upper']*len(H[k][replace])for k in active);common=1
   for k in active:common=math.lcm(common,vv[k]['common_QQ_denominator']*hp[k][replace][0])
   nbits=max([vv[k]['integer_numerator_height_upper']+hp[k][replace][1]+max(1,min(vv[k]['degree']+1,degree(H[k][replace])+1)).bit_length()+(common//(vv[k]['common_QQ_denominator']*hp[k][replace][0])).bit_length()for k in active]+[0])+4
   byte=terms*(32+math.ceil((nbits+common.bit_length())*math.log10(2)))
   allpairs+=pairs;allterms+=terms;bytes_upper+=byte;candidates.append({'replace_selected_row_position':replace,'replacement_poly_degree_upper':deg,'polynomial_terms_upper':terms,'naive_multiplication_pair_upper':pairs,'coefficient_num_bit_upper':nbits,'coefficient_common_den_bit':common.bit_length(),'raw_JSON_bytes_upper':byte})
  profiles.append({'K2_light':data[light]['physical_light'],'new_known_unit_LCD_powers5':mx,'selected8_row_degree_upper':max(v['degree']for v in vv),'fixed8_replacements':candidates})
 cert={'source_shape15x98':True,'compact_original_row_sources':['q5cleared','G18 times7','K2 times7'],'fixed_target_local98_columns8':cs,'old_selected_row_order':a['selected_rows'],
       'exact_future_formula':'For each newK2 rowv clear its known-unitLCD, then D_replace(j,i)=sum_k v_clear[k]*adj(S)[k,i]. This is a literal one-row replacement minor with oldrow/targetorder; no unknown denominator.',
       'positive_uniform_gate':'If gcd(oldD,all56replacements) becomes a known five-factor curve unit, selected8targetcolumns have rank8 at EVERY P100point, so originalquartic>=149 there. Literal QQgcd/Bezout or equivalent exact certificate would still be needed.',
       'negative_scope':'A remaining unknown gcd only marks the unresolved finite cover. One-row replacements all vanish when oldS has corank>=2 even if newrows can rescue; no actual failure follows.',
       'profiles7':profiles,'total_naive_pair_upper':allpairs,'total56_polynomial_terms_upper':allterms,'raw56_JSON_bytes_upper':bytes_upper,'oldD_degree':degree(D),'oldAdj_max_degree':max(degree(p)for row in H for p in row),
       'next_safeguard':'Before any dot/gcd authorization, known-unit stripping of oldD may reduce the univariate degree and witness size. Do not automatically evaluate56dots, gcd, higher replacements or a15x98rank.',
       'no_polynomial_dot_minor_factor_strip_gcd_or_rank_performed':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 report={'id':ID,'source_only':True,'fixed_replacement_candidates56':True,'total56_polynomial_terms_upper':allterms,'total_naive_pair_upper':allpairs,'raw_JSON_bytes_upper':bytes_upper,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True,'no_dot_gcd_or_rank_job':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
