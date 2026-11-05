#!/usr/bin/env python3
"""Four cachedold-C2 contractions for the actualDquotient7 genuineell-torsiontarget."""
import argparse,hashlib,json,re,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_GENUINE_T7_CACHED_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63/report.json'
SOURCE_SHA='a55309b43bd58e43f9a8516fca219213e38881b6142e66f16a1b09c8505c4983'
EVIDENCE='reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63/actual-filtered23-Dfree16-quotient7-and-genuine4-source.json'
GATE_REPORT='reports/AMBIENT_B18_C_GENUINE_T7_VALUE_GATE_20261105_v63/report.json'
GATE='reports/AMBIENT_B18_C_GENUINE_T7_VALUE_GATE_20261105_v63/four-cached-values-code-and-whole-positive-custody.json'
C2='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
MAX_PRECHECK=MAX_VALUE=0

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def check(n):
 global MAX_PRECHECK
 MAX_PRECHECK=max(MAX_PRECHECK,n);assert n<=2048
def note(q):
 global MAX_VALUE
 MAX_VALUE=max(MAX_VALUE,abs(q.numerator).bit_length(),q.denominator.bit_length());assert MAX_VALUE<=2048
 return q
def parse(x):
 p=str(x).split('/');assert len(p)in[1,2]
 check(max(abs(int(v)).bit_length()for v in p));return note(Q(str(x)))
def add(a,b):
 if not a:return b
 if not b:return a
 check(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def times(a,b):
 if not a or not b:return Q(0)
 check(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def neg(a):return note(-a)
def conv(a,b):
 v=[Q(0)]*4
 for n in range(4):
  for k in range(n+1):v[n]=add(v[n],times(a[k],b[n-k]))
 return v
def inverse(x):
 assert x[0]==1
 v=[Q(1),Q(0),Q(0),Q(0)]
 for n in range(1,4):
  for k in range(1,n+1):v[n]=add(v[n],neg(times(x[k],v[n-k])))
 return v

def mm(A,B):
 out=[]
 for row in A:
  outrow=[]
  for j in range(len(B[0])):
   v=Q(0)
   for k,a in enumerate(row):v=add(v,times(a,B[k][j]))
   outrow.append(v)
  out.append(outrow)
 return out

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-peer-report',required=True);args=parser.parse_args()
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648));started=time.monotonic()
 assert sha(ROOT/SOURCE)==SOURCE_SHA
 src=json.loads((ROOT/SOURCE).read_text());d=json.loads((ROOT/EVIDENCE).read_text());gr=json.loads((ROOT/GATE_REPORT).read_text());gate=json.loads((ROOT/GATE).read_text())
 pp=Path(args.source_peer_report)
 if not pp.is_absolute():pp=ROOT/pp
 peer=json.loads(pp.read_text());assert peer['all_passed']
 for h in[SOURCE_SHA,sha(ROOT/GATE_REPORT),sha(Path(__file__))]:assert h in peer['input_sha256'].values()
 assert sha(ROOT/EVIDENCE)==src['evidence_sha256'][Path(EVIDENCE).name]
 assert sha(ROOT/GATE)==gr['evidence_sha256'][Path(GATE).name]and sha(ROOT/C2)==gr['input_sha256'][C2]
 doc=json.loads((ROOT/C2).read_text())
 owned=count(src)+count(d)+count(gr)+count(gate)+count(peer)+count(doc)+12000
 assert owned<=25000<30000<150000
 raw=sum(p.stat().st_size for p in[ROOT/SOURCE,ROOT/EVIDENCE,ROOT/GATE_REPORT,ROOT/GATE,ROOT/C2,pp]);assert raw+100000<=600000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 (out/'pre-arithmetic-source-and-custody.json').write_text(json.dumps({'source_sha256':SOURCE_SHA,'peer_sha256':sha(pp),'gate_sha256':sha(ROOT/GATE_REPORT),
  'owned_numeric_upper':owned,'RAW_input_bytes':raw,'all_four_response_values_unpaired':True},indent=2)+'\n')
 dots={};products=[]
 try:
  byrow={r['physical_row']:dict(r['delta'])for r in doc['physical12_coefficients']}
  old={tuple(k):v for k,v in d['known_old27_leaf_values_copied']};keys={tuple(k)for op in d['four_ell_coefficients_UNPAIRED']for k,w in op['operator_old_C2_only']}
  assert len(keys)==24 and all(k[0]=='cache'for k in keys)
  leaves={}
  for k in sorted(keys):
   q=parse(byrow[k[1]].get(k[2],'0'));assert q==parse(old[k]);leaves[k]=q
  for op in d['four_ell_coefficients_UNPAIRED']:
   v=Q(0)
   for key,w in op['operator_old_C2_only']:
    q=leaves[tuple(key)];p=times(parse(w),q);v=add(v,p)
    products.append([op['label'],key,w,str(q),str(p)])
   dots[op['label']]=v
  a,b,c,e=[dots[k]for k in['alpha_A','alpha_B','beta_A','beta_B']]
  E=[[Q(0)]*7 for _ in range(7)];E[2][0]=1;E[5][1]=a;E[6][1]=b;E[5][2]=c;E[6][2]=e;E[5][3]=1;E[6][4]=1
  K=[[Q(0)]*4 for _ in range(7)];K[1][0]=1;K[3][0]=neg(a);K[4][0]=neg(b);K[2][1]=1;K[3][1]=neg(c);K[4][1]=neg(e);K[5][2]=1;K[6][3]=1
  assert mm(E,K)==[[Q(0)]*4 for _ in range(7)]
  assert[[E[i][j]for j in[0,3,4]]for i in[2,5,6]]==[[int(i==j)for j in range(3)]for i in range(3)]
  assert[[K[i][j]for j in range(4)]for i in[1,2,5,6]]==[[int(i==j)for j in range(4)]for i in range(4)]
  E2=mm(E,E);assert c or e
  assert E2==[[c if i==5 and j==0 else e if i==6 and j==0 else Q(0)for j in range(7)]for i in range(7)]
  assert mm(E2,E)==[[Q(0)]*7 for _ in range(7)]
  assert MAX_PRECHECK<=gate['whole_before_operation_bound']
 except BaseException as err:
  (out/'ROOT_ONLY_STOP_CHECKPOINT.json').write_text(json.dumps({'status':'STOP_no_cap_raise','error':repr(err),'completed_response_values':{k:str(v)for k,v in dots.items()},
   'max_precheck_bits':MAX_PRECHECK,'max_value_bits':MAX_VALUE},indent=2)+'\n');raise
 data={'all4_actual_cached_coefficients':{k:str(v)for k,v in dots.items()},'all60_termwise_products':products,
  'all24_original_C2_leaf_bindings':[[list(k),str(v)]for k,v in sorted(leaves.items())],
  'actual_E7':[[str(v)for v in row]for row in E],'actual_kernel7_by4':[[str(v)for v in row]for row in K],
  'all28_EK_products_zero':True,'fixed_rank3_I3_rows':[2,5,6],'fixed_rank3_I3_columns':[0,3,4],
  'kernel_freeI4_rows':[1,2,5,6],'actual_E2':[[str(v)for v in row]for row in E2],'beta_vector_nonzero_verified':True,'actual_E3_zero':True,
  'actual_E7_ranks':[7,3,1,0],'actual_genuine_C_ell_dimension':4,
  'basis_scope':'TheseareactualgenuineC[ell]coordinatesviatheSOURCEisomorphismkerE23/imE³→kerE7, relative toactualprimitiveD-freeF. NoexpandedC4torsionrepresentativesarecomputed.',
  'R0_transport_still_separate':'FNOTR0stable; doNOTreadr/zactionsfromthenaiquotient. No connectingRegDivell ororiginalH0values.',
  'no_new_source_physical_gauge_IFT_jet':True,'whole_max_precheck_bits':MAX_PRECHECK,'whole_max_value_bits':MAX_VALUE,'owned_numeric_upper':owned,'RAW_input_bytes':raw}
 ep=out/'four-cached-coefficients-actual-E7-and-genuine-four-target-coordinates.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'ACTUAL_D_TARGET_E7_RANKS7_3_1_0_GENUINE4_COORDINATES',
  'original_H0_nil_or_connecting_values_computed':False,'input_sha256':{SOURCE:SOURCE_SHA,EVIDENCE:sha(ROOT/EVIDENCE),GATE_REPORT:sha(ROOT/GATE_REPORT),GATE:sha(ROOT/GATE),C2:sha(ROOT/C2),str(pp.relative_to(ROOT)):sha(pp)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert raw+sum(p.stat().st_size for p in out.iterdir())<=600000
 print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
