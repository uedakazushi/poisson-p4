#!/usr/bin/env python3
"""Freeze unrun954-node values, literal39 inputs and positive whole budgets."""
import ast,hashlib,importlib.util,json,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_full_tauAB12_input_interface_20261105_v63 as iface

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_TAUAB12_WHOLE_GATE_20261105_v63'
SOURCE='reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63/report.json'
CIRCUIT='reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63/full-tauAB12-quartic-scalar-circuit-and-source-proof.json'
BOUNDS='reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63/all-scalar-majorants-no-value-evaluation.json'
WORKER='code/ambient_B18_C_full_tauAB12_values_20261105_v63.py'
IFACE='code/ambient_B18_C_full_tauAB12_input_interface_20261105_v63.py'
C4='reports/AMBIENT_B18_C4_HB_RESIDUAL_SOURCE_20261105_v63/report.json'
C4PEER='reports/AMBIENT_B18_C4_HB_RESIDUAL_SOURCE_20261105_v63_REVIEW/report.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(s):
 a=str(s).split('/');assert max(abs(int(x)).bit_length()for x in a)<=2048
 return Q(str(s))
def main():
 started=time.monotonic();signal.alarm(60)
 src=json.loads((ROOT/SOURCE).read_text());c=json.loads((ROOT/CIRCUIT).read_text());b=json.loads((ROOT/BOUNDS).read_text())
 c4=json.loads((ROOT/C4).read_text());p4=json.loads((ROOT/C4PEER).read_text());assert c4['all_passed']and p4['all_passed']
 assert src['all_passed']and len(c['nodes'])==954 and len(c['provided_inputs'])==39
 vals,checks,metrics=iface.read39(src,c,parse)
 # The known inputs may be copied and checked, but NO SLP output is paired.
 assert len(vals)==39
 H=b['max_height'];circuit_op=max(b['max_precheck'],H);rank_op=18*H+5
 assert H==86 and circuit_op==171 and rank_op==1553<2048
 owned=(iface.count(src)+iface.count(c)+iface.count(b)+5000+954+39+
        metrics['largest_single_parsed_source_numeric_count']+14000+1000)
 assert owned<=45000<60000<150000
 raw=sum((ROOT/p).stat().st_size for p in[SOURCE,CIRCUIT,BOUNDS])+metrics['literal_source_RAW_bytes']
 assert raw+80000+50000+180000<1000000
 worker=ROOT/WORKER;tree=ast.parse(worker.read_text())
 spec=importlib.util.spec_from_file_location('unrun_full_tauAB12_values',worker);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 assert m.SOURCE_SHA==sha(ROOT/SOURCE)and not(ROOT/'reports'/m.STEM).exists()
 x=list(map(m.parse,['1','2/3','-3/7','4/11']));inv=m.inverse(x);assert m.conv(x,inv)==[1,0,0,0]==m.conv(inv,x)
 assert m.divide(m.parse('2/7'),m.parse('5/11'))==Q(22,35)
 try:m.times(m.Q(1<<2040),m.Q(1<<20))
 except AssertionError:pass
 else:raise AssertionError('before-operation guard failed')
 help_result=subprocess.run([str(ROOT/'.venv/bin/python'),str(worker),'--help'],capture_output=True,text=True,timeout=10)
 assert help_result.returncode==0 and'--source-peer-report'in help_result.stdout
 data={
  'scope':'SAME ONE954-node cachedscalar evaluation andfixedrankcertificates; no roots/source/gauge/IFT/jet replay.',
  **metrics,'all39_literal_inputs_match_pinned_source_locators':checks,
  'circuit_scalar_nodes':954,'circuit_max_stored_height_bound':H,'circuit_before_operation_bits':circuit_op,
  'rank_and_factorization_before_operation_bits':rank_op,'hard_bits':2048,
  'positive_postrank_proof':'AlltargetentriesandknownHBconstants haveheight<=H. A2x2det hasheight<=4H+1; eachinverseentry<=5H+1 andbothproductsoperationally<=12H+3. Cramernumerators<=4H+1, ratios<=8H+2; multiplyingbackwithbasisentriesandaddingtwoproductsisboundedbeforeoperationby18H+5. No lowerboundondeterminantsisassumed; nonzeroactualdet swapsnumerator/denominator withoutheightgrowth.',
  'generalized_C4_source_proof':'ForANYnewactualquarticW, residualLk=Wk-(vk/v0)W0. Onorderedrows[0,k,2,5], columns[HB0,W,ell4A,ell4B] haveblockmatrix[[v0,W0,0,0],[vk,Wk,0,0],[p,q,1,0],[s,t,0,1]], determinantv0*Wk-vk*W0=v0*Lk independentofp,q,s,t. Fork1,3or4 theliteralcutentriesattoprowsare0. CombinedwithreviewedC4m5zero/prefix21/lower23, Lk!=0⇒C4length23/e3.',
  'C4_negative_scope':'All3 residualcutcoordinates0 iscriterionnegativeonly, not e>=4 orC4length24.',
  'owned_numeric_upper_with5000_peer_and_output_reserves':owned,
  'value_static_numeric':45000,'value_hard_numeric':60000,'global_hard_numeric':150000,
  'existing_RAW_bytes_including_all5_literal_input_files':raw,
  'future_peer_RAW_reserve':80000,'gate_custody_RAW_reserve':50000,'future_output_RAW_reserve':180000,
  'whole_RAW_cap':1000000,'global_RAW_cap':33554432,'seconds_cap':60,'memory_cap':2147483648,
  'future_worker_sha256':sha(worker),'import_only_help_passed':True,'unrelated_inverse_division_and_guard_fixtures_passed':True,
  'AST_functions':[n.name for n in tree.body if isinstance(n,ast.FunctionDef)],
  'all_full_tauAB_responses_uncomputed':True,
  'positive_scope':'Nonzeroanyresidual⇒knownrank2/fullHF4=4,b2,e3or4. Ifcutresidualk1,3or4nonzero⇒C4rank4/L23/e3. No originalH0/nilcube claim.',
  'allzero_scope':'HB+all12tauABknownsubmodulerank1only. Completequarticspace/b/eopen.',
  'future_command':'.venv/bin/python '+WORKER+' --source-peer-report reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63_REVIEW/report.json',
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'whole-source-values-rank-and-C4-positive-profile.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'UNRUN_FULL_TAUAB12_VALUES_CODE_WHOLE_POSITIVE_GATE_NO_RESPONSES',
  'input_sha256':{p:sha(ROOT/p)for p in[SOURCE,CIRCUIT,BOUNDS,WORKER,IFACE,C4,C4PEER]},
  'code_sha256':sha(Path(__file__)),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
