#!/usr/bin/env python3
"""One gated954-node cached full-tauAB scalar circuit, no root/source recomputation."""
import argparse,hashlib,json,re,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63/report.json'
SOURCE_SHA='2debebb789b978069fc3dd5863cb5ff5f834366add91716c34065209f3dd852c'
CIRCUIT='reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63/full-tauAB12-quartic-scalar-circuit-and-source-proof.json'
BOUNDS='reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63/all-scalar-majorants-no-value-evaluation.json'
GATE_REPORT='reports/AMBIENT_B18_C_FULL_TAUAB12_WHOLE_GATE_20261105_v63/report.json'
GATE='reports/AMBIENT_B18_C_FULL_TAUAB12_WHOLE_GATE_20261105_v63/whole-source-values-rank-and-C4-positive-profile.json'
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

def divide(a,b):
 assert b
 check(max(abs(b.numerator).bit_length(),b.denominator.bit_length()))
 return times(a,note(Q(b.denominator,b.numerator)))

def determinant2(a):return add(times(a[0][0],a[1][1]),neg(times(a[0][1],a[1][0])))

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-peer-report',required=True);args=parser.parse_args()
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648));started=time.monotonic()
 import ambient_B18_C_full_tauAB12_input_interface_20261105_v63 as iface
 assert sha(ROOT/SOURCE)==SOURCE_SHA
 source=json.loads((ROOT/SOURCE).read_text());c=json.loads((ROOT/CIRCUIT).read_text());bounds=json.loads((ROOT/BOUNDS).read_text())
 gr=json.loads((ROOT/GATE_REPORT).read_text());gate=json.loads((ROOT/GATE).read_text())
 pp=Path(args.source_peer_report)
 if not pp.is_absolute():pp=ROOT/pp
 peer=json.loads(pp.read_text());assert peer['all_passed']
 for h in[SOURCE_SHA,sha(ROOT/GATE_REPORT),sha(Path(__file__))]:assert h in peer['input_sha256'].values()
 for p in[CIRCUIT,BOUNDS]:assert sha(ROOT/p)==source['evidence_sha256'][Path(p).name]
 assert sha(ROOT/GATE)==gr['evidence_sha256'][Path(GATE).name]
 ifacepath=ROOT/'code/ambient_B18_C_full_tauAB12_input_interface_20261105_v63.py'
 assert sha(ifacepath)==gr['input_sha256'][str(ifacepath.relative_to(ROOT))]
 assert len(c['nodes'])==954 and len(c['provided_inputs'])==39
 owned=count(source)+count(c)+count(bounds)+count(gr)+count(gate)+count(peer)+954+39+gate['largest_single_parsed_source_numeric_count']+14000
 assert owned<=45000<60000<150000
 raw=sum(p.stat().st_size for p in[ROOT/SOURCE,ROOT/CIRCUIT,ROOT/BOUNDS,ROOT/GATE_REPORT,ROOT/GATE,pp])+gate['literal_source_RAW_bytes']
 assert raw+180000<=1000000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 (out/'pre-arithmetic-source-and-whole-custody.json').write_text(json.dumps({'source_sha256':SOURCE_SHA,'source_peer_sha256':sha(pp),'gate_sha256':sha(ROOT/GATE_REPORT),
  'code_sha256':sha(Path(__file__)),'owned_numeric_upper':owned,'RAW_input_bytes':raw,'all_tauAB_responses_not_yet_evaluated':True},indent=2)+'\n')
 values=[];i=-1
 try:
  inputs,checks,metrics=iface.read39(source,c,parse)
  assert metrics['largest_single_parsed_source_numeric_count']==gate['largest_single_parsed_source_numeric_count']
  for i,node in enumerate(c['nodes']):
   op=node[0]
   if op=='c':v=parse(str(node[1])+'/'+str(node[2]))
   elif op=='input':v=inputs[i]
   elif op=='a':v=add(values[node[1]],values[node[2]])
   elif op=='m':v=times(values[node[1]],values[node[2]])
   else:raise AssertionError(op)
   U,D,h=bounds['bounds_U_D_height'][i];U,D=int(U),int(D)
   check(U.bit_length()+v.denominator.bit_length())
   assert abs(v.numerator)<=U*v.denominator and D%v.denominator==0
   values.append(v)
  circuit_pre,circuit_height=MAX_PRECHECK,MAX_VALUE
  assert circuit_pre<=gate['circuit_before_operation_bits']
  norm=[ [values[n]for n in r['NF4_full6_nodes']]for r in c['outputs3']]
  residual=[ [values[n]for n in r['residual5_nodes']]for r in c['outputs3']]
  v=list(map(parse,c['known_full_HB_vector']));slopes=list(map(parse,c['five_residual_dual_slopes_vs_coordinate0']))
  assert v[0]==parse('-3/5')
  for col,res in zip(norm,residual):
   assert res==[add(col[k],neg(times(slopes[k-1],col[0])))for k in range(1,6)]
  positive=[(j,k)for j in range(3)for k in range(1,6)if residual[j][k-1]]
  cut_positive=[(j,k)for j in range(3)for k in[1,3,4]if residual[j][k-1]]
  if positive:
   j,k=positive[0];W=norm[j];block=[[v[0],W[0]],[v[k],W[k]]];det=determinant2(block)
   assert det==times(v[0],residual[j][k-1])and det
   inv=[[divide(block[1][1],det),divide(neg(block[0][1]),det)],
        [divide(neg(block[1][0]),det),divide(block[0][0],det)]]
   for A,B in[(block,inv),(inv,block)]:
    assert[[add(times(A[a][0],B[0][b]),times(A[a][1],B[1][b]))for b in range(2)]for a in range(2)]==[[1,0],[0,1]]
   coefficients=[]
   for X in[v]+norm:
    # Cramer ratios give the proved8H+2 stored bound; no unboundedsolve.
    n0=determinant2([[X[0],W[0]],[X[k],W[k]]]);n1=determinant2([[v[0],X[0]],[v[k],X[k]]])
    ab=[divide(n0,det),divide(n1,det)]
    assert all(add(times(v[a],ab[0]),times(W[a],ab[1]))==X[a]for a in range(6))
    coefficients.append(list(map(str,ab)))
   rank=2;cert={'source_pivot_index':j,'NF_rows':[0,k],'minor':str(det),'inverse':[[str(x)for x in row]for row in inv],
    'both_inverse_products':True,'all4_actual_source_columns_factor':coefficients}
  else:
   coefficients=[]
   for X in norm:
    a=divide(X[0],v[0]);assert all(times(a,x)==y for x,y in zip(v,X));coefficients.append(str(a))
   rank=1;cert={'nonzero_HB_anchor':'-3/5','all3_tauAB_columns_proportional_to_HB':coefficients}
  c4cert=None
  if cut_positive:
   j,k=cut_positive[0];W=norm[j];minor=times(v[0],residual[j][k-1]);assert minor
   c4cert={'pivot_source_index':j,'rows_ordered':[0,k,2,5],
    'columns':['HB0','tauAB_pivot','literalell4A','literalell4B'],'universal_block_minor':str(minor),
    'matrix':[[str(v[a]),str(W[a]),str(int(a==2)),str(int(a==5))]for a in[0,k,2,5]],
    'known_C4_length_upper23_equals_lower23':True}
  assert MAX_PRECHECK<=gate['rank_and_factorization_before_operation_bits']
 except BaseException as err:
  (out/'ROOT_ONLY_STOP_CHECKPOINT.json').write_text(json.dumps({'status':'STOP_no_cap_raise','error':repr(err),'completed_scalar_node':i,
   'available_output_nodes':[[n,str(values[n])]for r in c['outputs3']for n in r['NF4_full6_nodes']+r['residual5_nodes']if n<len(values)],
   'max_precheck_bits':MAX_PRECHECK,'max_value_bits':MAX_VALUE},indent=2)+'\n');raise
 data={'literal39_cached_input_checks':checks,'full6_by3_tauAB_columns_as3columns':[[str(x)for x in col]for col in norm],
  'all15_residual_values_as3_by5':[[str(x)for x in col]for col in residual],
  'actual_known_HB_plus_tauAB_rank':rank,'fixed_rank_and_all4column_certificate':cert,'C4_rank4_certificate_if_positive':c4cert,
  'actual_full_C_HF4_if_rank2':4 if rank==2 else'OPEN_4_OR5','exact_comparison_C_e3_proved':bool(cut_positive),
  'long_C_exponent_scope':'3'if cut_positive else('3_OR4'if rank==2 else'e>=3_OPEN'),
  'rank1_scope':'OnlyHB+alltauAB12submoduleknownrank1; nofullquarticrank/eclaim.',
  'no_new_physical_gauge_IFT_jet_delta4':True,'original_H0_nilcube_or_exact_nilindex_not_determined':True,
  'circuit_max_precheck_bits':circuit_pre,'circuit_max_value_bits':circuit_height,'whole_max_precheck_bits':MAX_PRECHECK,'whole_max_value_bits':MAX_VALUE,
  'owned_numeric_upper':owned,'RAW_input_bytes':raw}
 ep=out/'all-full-tauAB12-values-and-fixed-known-rank-C4-certificates.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'KNOWN_TAUAB_POSITIVE_C4_LENGTH23_E3'if cut_positive else('KNOWN_TAUAB_POSITIVE_FULL_HF4_4_E3_OR4'if rank==2 else'ALL_TAUAB12_RESIDUALS_ZERO_KNOWN_SUBMODULE_RANK1_ONLY'),
  'actual_known_HB_plus_tauAB_rank':rank,'exact_comparison_C_e3_proved':bool(cut_positive),
  'input_sha256':{SOURCE:SOURCE_SHA,CIRCUIT:sha(ROOT/CIRCUIT),BOUNDS:sha(ROOT/BOUNDS),GATE_REPORT:sha(ROOT/GATE_REPORT),GATE:sha(ROOT/GATE),str(pp.relative_to(ROOT)):sha(pp)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert raw+sum(p.stat().st_size for p in out.iterdir())<=1000000
 print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
