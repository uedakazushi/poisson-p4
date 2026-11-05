#!/usr/bin/env python3
"""Integer-only prior bound for independent scalar-symbol source expansion."""
import hashlib,json,math,pathlib
from collections import Counter
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
 out=ROOT/'reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63_REVIEW'
 sourcepath=ROOT/'reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63/full-tauAB12-quartic-scalar-circuit-and-source-proof.json'
 data=json.loads(sourcepath.read_text());nodes=data['nodes'];roots=[]
 for row in data['outputs3']:roots+=row['NF4_full6_nodes']+row['residual5_nodes']
 for row in data['quartic_polynomial_nodes']:
  roots+=[node for key in ['WA4_polynomial','WB4_polynomial'] for mon,node in row[key]]
 uses=Counter()
 for node in nodes:
  if node[0]in ['a','m']:uses.update(node[1:])
 uses.update(roots);bounds=[];degrees=[];live={};peak=0
 for i,node in enumerate(nodes):
  if node[0]=='c':value=int(node[1]!=0);degree=0
  elif node[0]=='input':value=1;degree=1
  else:
   left,right=node[1:]
   degree=max(degrees[left],degrees[right]) if node[0]=='a' else degrees[left]+degrees[right]
   value=min(bounds[left]+bounds[right] if node[0]=='a' else bounds[left]*bounds[right],math.comb(39+degree,degree))
  bounds.append(value);degrees.append(degree);live[i]=value;peak=max(peak,sum(live.values()))
  if node[0]in ['a','m']:
   for child in node[1:]:
    uses[child]-=1
    if not uses[child]:del live[child]
  if not uses[i]:del live[i]
 base=data['source_bound_metadata']['source_owned_upper'];whole=base+8*peak+24000
 assert whole<150000 and max(degrees)==2
 proof={'status':'PRE_EXPANSION_SOURCE_BOUND_PASS_NO_RESPONSE_VALUES',
  'input_sha256':hashlib.sha256(sourcepath.read_bytes()).hexdigest(),'nodes':len(nodes),'symbols':39,
  'max_symbol_degree':2,'all_node_support_upper':bounds,'all_node_degrees':degrees,'retained_root_ids':roots,
  'lastuse_peak_support_upper':peak,'input_original_source_whole_upper':base,
  'whole_symbolic_proof_upper':whole,'hard_numeric':150000,
  'reason':'8 slots per term cover circuit plus independent coefficient-polynomial copies;24000 additional scratch/controller/publication reserve',
  'no_signed_response_or_mu_pairing':True,'code_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
 dest=out/'pre-expansion-support-and-whole-bound.json';assert not dest.exists();dest.write_text(json.dumps(proof,indent=2)+'\n')
 print(json.dumps({k:v for k,v in proof.items() if k not in ['all_node_support_upper','all_node_degrees','retained_root_ids']}))
if __name__=='__main__':main()
