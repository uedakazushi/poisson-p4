#!/usr/bin/env python3
"""Four unpaired oldC2 contractions: exact code/input/whole positive custody."""
import ast,hashlib,importlib.util,json,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_GENUINE_T7_VALUE_GATE_20261105_v63'
SOURCE='reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63/report.json'
EVIDENCE='reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63/actual-filtered23-Dfree16-quotient7-and-genuine4-source.json'
WORKER='code/ambient_B18_C_genuine_T7_cached_values_20261105_v63.py'
C2='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
ORIGINAL='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/report.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def parse(s):
 a=str(s).split('/');assert max(abs(int(x)).bit_length()for x in a)<=2048
 return Q(str(s))
def degree(x):return max((sum(mon)for mon,q in x),default=-1)
def pathstats(A,B):
 maximum=total=0
 for i in range(len(A)):
  for j in range(len(B[0])):
   terms=0
   for k in range(len(B)):
    if A[i][k]and B[k][j]:assert degree(A[i][k])+degree(B[k][j])<=1;terms+=1
   maximum=max(maximum,terms);total+=terms
 return maximum,total
def main():
 started=time.monotonic();signal.alarm(60);src=json.loads((ROOT/SOURCE).read_text());d=json.loads((ROOT/EVIDENCE).read_text())
 doc=json.loads((ROOT/C2).read_text());original=json.loads((ROOT/ORIGINAL).read_text())
 assert src['all_passed']and sha(ROOT/C2)==original['input_sha256'][C2]
 byrow={r['physical_row']:dict(r['delta'])for r in doc['physical12_coefficients']};old={tuple(k):s for k,s in d['known_old27_leaf_values_copied']}
 keys=sorted({tuple(k)for op in d['four_ell_coefficients_UNPAIRED']for k,w in op['operator_old_C2_only']});assert len(keys)==24
 checks=[]
 for k in keys:
  assert k[0]=='cache';a=parse(byrow[k[1]].get(k[2],'0'));assert a==parse(old[k]);checks.append([list(k),str(a)])
 assert sum(len(x['operator_old_C2_only'])for x in d['four_ell_coefficients_UNPAIRED'])==60
 E=d['symbolic_E7_rows7_columns7'];K=d['symbolic_kernel7_by4']
 EK=pathstats(E,K);EE=pathstats(E,E);assert EK[0]<=2 and EE[0]<=1
 H=d['future_coefficient_height_bound'];op=2*H+1;assert H==43 and op==87<2048
 owned=count(src)+count(d)+count(doc)+5000+11000+1500;assert owned<25000<30000<150000
 raw=sum((ROOT/p).stat().st_size for p in[SOURCE,EVIDENCE,C2]);assert raw+80000+30000+100000<600000
 worker=ROOT/WORKER;tree=ast.parse(worker.read_text());spec=importlib.util.spec_from_file_location('unrun_T7_fourcache',worker);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 assert m.SOURCE_SHA==sha(ROOT/SOURCE)and not(ROOT/'reports'/m.STEM).exists()
 x=list(map(m.parse,['1','2/3','-3/7','4/11']));iv=m.inverse(x);assert m.conv(x,iv)==[1,0,0,0]==m.conv(iv,x)
 try:m.times(m.Q(1<<2040),m.Q(1<<20))
 except AssertionError:pass
 else:raise AssertionError('before-operation guard failed')
 result=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,text=True,timeout=10)
 assert result.returncode==0 and'--source-peer-report'in result.stdout
 data={
  'scope':'Fouralreadycachedold-C2scalarcontractionsandDtargetE7/kernelcoordinatesonly. NotR0actiontransportororiginalconnectingvalues.',
  'all24_original_C2_leaf_copies_verified':checks,'four_response_values_not_evaluated':True,'weighted_terms':60,
  'source_response_height_bound':H,'whole_before_operation_bound':op,'hard_bits':2048,
  'all_matrix_nonzero_path_products_degree_atmost1':True,'EK_max_nonzero_terms_and_total':EK,'E2_max_nonzero_terms_and_total':EE,
  'operational_proof':'Fourdotpartials boundedbySOURCEH43; allmatrixnonzeroproductshaveatmostoneunknownresponsefactor andonefixedunit, withatmost2nonzerotermsperEKentry. Thusstored<=H43andsumpre<=2H+1=87, includingE²/E³checks. Noalpha*beta productorunknowninverse.',
  'owned_numeric_upper_with5000_peer_and_output_reserves':owned,'static_numeric':25000,'hard_numeric':30000,'global_numeric':150000,
  'existing_RAW_bytes':raw,'peer_RAW_reserve':80000,'gate_RAW_reserve':30000,'output_RAW_reserve':100000,'RAW_cap':600000,'global_RAW':33554432,'seconds_cap':60,'memory_cap':2147483648,
  'worker_sha256':sha(worker),'import_only_help_passed':True,'unrelated_inverse_and_overflow_guard_fixtures_passed':True,
  'AST_functions':[n.name for n in tree.body if isinstance(n,ast.FunctionDef)],
  'future_mandatory':['all24literalC2leaves/SHA','all60products/fourdotvalues','all28EKzero/fixedI3/freeI4','betaAorBnonzero/E²rank1/E³zero','scopeFonlyDstable/noR0quotientactions'],
  'future_command':'.venv/bin/python '+WORKER+' --source-peer-report reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63_REVIEW/report.json',
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'four-cached-values-code-and-whole-positive-custody.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'UNRUN_DTARGET_E7_FOUR_OLD_CACHE_VALUES_CODE_WHOLE_GATE',
  'input_sha256':{p:sha(ROOT/p)for p in[SOURCE,EVIDENCE,WORKER,C2,ORIGINAL]},'code_sha256':sha(Path(__file__)),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
