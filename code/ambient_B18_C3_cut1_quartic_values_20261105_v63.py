#!/usr/bin/env python3
"""ONE cached CUT1 residual, with the other coordinate kept formal."""
import argparse,hashlib,json,re,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C3_CUT1_QUARTIC_SOURCE_20261105_v63/report.json'
SOURCE_SHA='ed1678417ff8e459de67d38ced0c8c835aab2a7ea8a1eeb315867a7570f9b633'
EVIDENCE='reports/AMBIENT_B18_C3_CUT1_QUARTIC_SOURCE_20261105_v63/literal-cut-quartic-source-and-three-unpaired-cache-operators.json'
C2='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
MAX_PRE=MAX_VALUE=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def guard(n):
 global MAX_PRE
 MAX_PRE=max(MAX_PRE,n);assert n<=2048
def note(q):
 global MAX_VALUE
 MAX_VALUE=max(MAX_VALUE,abs(q.numerator).bit_length(),q.denominator.bit_length());assert MAX_VALUE<=2048
 return q
def parse(s):
 guard(max(abs(int(v)).bit_length()for v in str(s).split('/')));return note(Q(s))
def add(a,b):
 if not a:return b
 if not b:return a
 guard(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def mul(a,b):
 if not a or not b:return Q(0)
 guard(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def neg(a):guard(max(abs(a.numerator).bit_length(),a.denominator.bit_length()));return note(-a)
def div(a,b):
 assert b;guard(max(abs(b.numerator).bit_length(),b.denominator.bit_length()))
 return mul(a,note(Q(b.denominator,b.numerator)))
def padd(a,b):return[add(x,y)for x,y in zip(a,b)]
def pmul(a,b):
 c=[Q(0)]*3
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<3:c[i+j]=add(c[i+j],mul(x,y))
 return c
def mm(a,b):
 out=[]
 for row in a:
  rr=[]
  for col in zip(*b):
   p=[Q(0)]*3
   for x,y in zip(row,col):p=padd(p,pmul(x,y))
   rr.append(p)
  out.append(rr)
 return out
def encode(a):return[[list(map(str,p))for p in row]for row in a]


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-peer-report',required=True);args=parser.parse_args()
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648));started=time.monotonic()
 assert sha(ROOT/SOURCE)==SOURCE_SHA
 source=json.loads((ROOT/SOURCE).read_text());peer_path=Path(args.source_peer_report)
 if not peer_path.is_absolute():peer_path=ROOT/peer_path
 peer=json.loads(peer_path.read_text());assert peer['all_passed']and SOURCE_SHA in set(peer.get('input_sha256',{}).values())
 assert sha(ROOT/EVIDENCE)==source['evidence_sha256'][Path(EVIDENCE).name]
 for path,want in source['input_sha256'].items():assert sha(ROOT/path)==want
 ev=json.loads((ROOT/EVIDENCE).read_text());c2=json.loads((ROOT/C2).read_text())
 owned=count(source)+count(peer)+count(ev)+count(c2)+6000;assert owned<=18000<30000
 raw_paths=set(source['input_sha256'])|{SOURCE,EVIDENCE,str(peer_path.relative_to(ROOT))}
 raw=sum((ROOT/path).stat().st_size for path in raw_paths);assert raw+100000<=400000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 (out/'pre-arithmetic-source-and-budget.json').write_text(json.dumps({'source_sha256':SOURCE_SHA,'worker_sha256':sha(Path(__file__)),'peer_sha256':sha(peer_path),'owned_whole_bound':owned,'input_RAW_bytes':raw,'cached_dot_not_yet_evaluated':True},indent=2)+'\n')
 vals={r['physical_row']:dict(r['delta'])for r in c2['physical12_coefficients']};d=Q(0);products=[]
 for key,w in ev['residual_B_minus_121A_over48_operator']:
  tag,physical,index=key;assert tag=='cache';weight=parse(w);leaf=parse(vals[physical].get(index,'0'));p=mul(weight,leaf);d=add(d,p)
  products.append([key,str(weight),str(leaf),str(p)])
 a=parse(ev['known_quartic_vector'][0]);b=parse(ev['known_quartic_vector'][1]);slope=parse(ev['residual_dual_slope'])
 assert a==Q(-3,5)and add(b,neg(mul(slope,a)))==0
 determinant=mul(a,d)
 if d:
  # In targetcoordinates(A,L), actualCUT=(t,d), t unpaired/formal.
  zero=[Q(0)]*3;one=[Q(1),Q(0),Q(0)]
  M=[[[a,Q(0),Q(0)],[Q(0),Q(1),Q(0)]],[zero,[d,Q(0),Q(0)]]]
  inv=[[[div(Q(1),a),Q(0),Q(0)],[Q(0),neg(div(Q(1),determinant)),Q(0)]],
       [zero,[div(Q(1),d),Q(0),Q(0)]]]
  I=[[one,zero],[zero,one]];left,right=mm(M,inv),mm(inv,M);assert left==I==right
  proof={'actual_target_coordinate_change':'(A4,B4)→(A4,L=B4-121A4/48)',
         'unpaired_actual_CUT_A4':'formal t; not assignedzero or evaluated',
         'universal_matrix_over_QQ_t':encode(M),'inverse_over_QQ_t':encode(inv),
         'both_universal_products':{'M_inverse':encode(left),'inverse_M':encode(right)},
         'nonzero_actual_minor_independent_of_t':str(determinant)}
 else:
  proof={'actual_CUT_column_in_known_HB_line':True,'column_factor':'formal t/(-3/5), where t is unevaluated actual A4; no point value assigned',
         'full_quartic_rank_or_C3_length20_not_inferred':True}
 assert MAX_PRE<=ev['whole_fold_minor_inverse_rank_operation_bound']<=2048
 data={'all24_termwise_cached_products':products,'actual_residual':str(d),'actual_known_HB_CUT_minor':str(determinant),
       'CUT_A4_B4_operators_remain_unpaired':True,'rank_source_certificate':proof,
       'actual_C3_length19_proved':bool(d),'actual_C3_HF_if_positive':[3,6,6,4,0]if d else None,
       'cut3_torsion_pattern_if_positive':[1,1,2,'e>=3']if d else None,
       'longer_exponent_genuine_C_ell_original_H0_nil_not_determined':True,
       'owned_numeric_bound_including_metadata_and_output':owned,'input_RAW_bytes':raw,
       'maximum_preoperation_bits':MAX_PRE,'maximum_value_bits':MAX_VALUE,'new_source_gauge_physical_IFT_or_jet':False}
 ep=out/'one-CUT1-residual-and-formal-coordinate-rank-certificate.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'POSITIVE_CUT1_EXACT_C3_LENGTH19'if d else'ZERO_CUT1_SUBSET_ONLY',
    'actual_C3_length19_proved':bool(d),'longer_exponent_or_original_H0_nil_determined':False,
    'input_sha256':{SOURCE:SOURCE_SHA,EVIDENCE:sha(ROOT/EVIDENCE),C2:sha(ROOT/C2),str(peer_path.relative_to(ROOT)):sha(peer_path)},
    'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
    'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert raw+sum(p.stat().st_size for p in out.iterdir())<=400000
 print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()
