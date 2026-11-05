#!/usr/bin/env python3
"""SOURCE only: saved T0/root239, fixed a4 removal and next known-C Hasse budget."""
import argparse,hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T0_C_HASSE_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  return json.loads(pin(ROOT/'reports'/st/name).read_text())
 old=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63');assert old['all_passed']and old['output_terms']==7208
 peer=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW');assert peer['all_passed']
 pre=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_PREFLIGHT_20261104_v63','operator-and-cost.json');D=int(pre['nodes'][pre['output']]['QQ_denominator']);bound=int(pre['nodes'][pre['output']]['integer_l1_upper']);assert D.bit_length()==197
 p=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root239.poly');assert sha(p)==old['evidence_sha256'][p.name]
 jp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root239.json');assert sha(jp)==old['evidence_sha256'][jp.name]
 n=L=word=numword=denword=0;box=[[999,-999]for _ in range(3)];deg=[999,-999];seen=set()
 with p.open()as f:
  assert int(f.readline())==4;assert int(f.readline())==7208
  for line in f:
   aa,ee,dd,c=line.split();a,e,d=map(int,[aa,ee,dd]);cc=Q(c);assert D%cc.denominator==0
   scale=D//cc.denominator;assert scale*cc.denominator==D
   v=cc.numerator*scale;assert v;n+=1;L+=abs(v);word=max(word,abs(v).bit_length())
   numword=max(numword,abs(cc.numerator).bit_length());denword=max(denword,cc.denominator.bit_length())
   k=a*1048576+e*1024+d;assert k not in seen;seen.add(k)
   for i,w in enumerate([a,e,d]):box[i]=[min(box[i][0],w),max(box[i][1],w)]
   deg=[min(deg[0],a+e+d),max(deg[1],a+e+d)]
 assert n==7208 and L<=bound and L.bit_length()==326 and word==320
 assert box==[[0,20],[0,40],[0,18]]and deg==[0,57]
 owned=2*(7208+441+1539)+16000;assert owned==34376
 inner=7208*441+2*7208*41
 allprod=inner+2*7208*23+1200+20+1+2*7208+100
 header=list(map(int,[4,7208]));assert header==[4,7208]
 d=dict(SOURCE_ONLY=True,T0_definition='the existing T=B*kappa*tail2, KAPPA_EXTRA2/root239, not a new source point',
  source_form='root239.poly firstline4 is the aOFFSET; its7208 bodyexponents already describe U0=T0/a^4',
  input_a_shift_for_future_worker='0 if using .poly body; subtract4 ONLY if instead decoding .json originalexponents; NEVER both',
  all7208_distinct_monomials=True,common_D197=D,normalized_integer_maxbits320=word,normalized_integer_l1bits326=L.bit_length(),
  original_QQ_numerator_bits235=numword,original_QQ_denominator_bits145=denword,
  constant_conversion='d divides D;scale=D/d, verify scale*d=D;I=n*scale. No I*d versus n*D cross-product.',
  integer_normalization_file_NOT_created=True,only_scalar_conversion_for_bounds=True,
  U0_bounds=dict(a=[0,20],e=[0,40],d=[0,18],total_degree=[0,57]),
  C='100*a+18*e^2-9*e-2',a0='(-18*e^2+9*e+2)/100',root_integer_low_to_high=[2,9,-18],
  Hasse_formula='D*100^20*[z^k]U0(a0+z,e,d)=sum I_ijl*choose(i,k)*100^(20-i+k)*(2+9e-18e^2)^(i-k)*e^j*d^l',
  orders='k0..20 firstNONZERO only; exact v_C(U0)=v_C(T0) because C does not divide a',
  source_input_l1_hard326=326,HH_integer_partial_bits481=326+20+(100**20).bit_length()+2,
  proposed_unchanged_hard_integer518=518,QQ_denominator_clear_bits_upper331=331,
  cache441=441,one_accumulator1539=1539,support='e<=40+40=80,d<=18=>81*19=1539',
  whole_owned34376=owned,ownership='two cells per inputcoeff/key7208,cache441,acc1539 plus ALLmetadata/scratch16000',
  prepare_streaming='plain .poly one line at a time; keyset7208 only if needed for duplicate check; never decode full .json together with allintegerpayload',
  no_opaque_poly_products_or_unknown_inverses=True,
  native_inner_products_upper=inner,ALL_conversion_scaled_root_products_upper=allprod,
  fixed_scalar_conversion_divisibility_checks7208=7208,
  RAW_upper6581960=(9188+16000)*170+100000+2000000+200000,
  proposed_RAW_hard8000000=8000000,soft45_hard60_RAM2GB=True,
  no_C_power_Hasse_or_valuation_evaluated=True,
  fresh_source_code_independent_audit_and_parent_GO_required=True,
  no_GCD_other_factor_tail_or_bad_locus_inference=True,
  purpose='compact the OLD T0 alongside normalizedT1 for later complement analysis; same FULLU sufficient-only D12 scope')
 assert d['HH_integer_partial_bits481']==481 and d['RAW_upper6581960']==6581960
 (out/'T0-literal-frame-and-preflight.json').write_text(json.dumps(d,indent=2)+'\n')
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,old_T0_7208_input_bound=True,
  aoffset4_literal_scope_pinned=True,C_valuation_NOT_computed=True,worker_NOT_built_or_run=True,
  normalized_scalar_file_NOT_created=True,whole_owned34376=34376,HH_bits481=481,hard518=518,
  ALL_products_upper=allprod,RAW_upper6581960=6581960,input_sha256=pins,
  evidence_sha256={'T0-literal-frame-and-preflight.json':sha(out/'T0-literal-frame-and-preflight.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
