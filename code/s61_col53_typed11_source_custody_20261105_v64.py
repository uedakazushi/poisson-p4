"""Structural SOURCE/custody audit of immutable typed11; no scalar coefficient evaluation."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_TYPED11_SOURCE_CUSTODY_20261105_v64'
PRIMARY='reports/S61_COL53_TYPED11_CALLBACK_SOURCE_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
if __name__=='__main__':
 p=ROOT/PRIMARY/'report.json';assert sha(p)=='cecda91623ef6010af0dceba045038339ec910806bddeb794bcd547d7248f583'
 r=json.loads(p.read_text());assert r['all_passed']and r['SOURCE_typed_callbacks_only_no_matrix_RHS_VALUE']
 pins={str(p.relative_to(ROOT)):sha(p)}
 for f,h in r['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
 for f,h in r['evidence_sha256'].items():assert sha(p.parent/f)==h;pins[str((p.parent/f).relative_to(ROOT))]=h
 O=json.loads((p.parent/'actual-full166x13-and-RHS-typed11-operator.json').read_text())
 ns=O['scalar_and_RHS_nodes'];types=[]
 for k,n in enumerate(ns):
  assert n['id']==k and all(0<=j<k for j in n['args'])
  op=n['op'];args=n['args']
  if op in['QQ_constant','known_parameter_polynomial','external_known_polynomial','external_full_SOURCE_DAG','external_current_cleared_numerator']:
   assert args==[];kind='parameter_scalar'
  elif op=='external_current_RHS':assert args==[];kind='RHS'
  elif op in['parameter_multiply','parameter_subtract']:
   assert len(args)==2 and all(types[j]=='parameter_scalar'for j in args);kind='parameter_scalar'
  elif op=='certified_parameter_unit_inverse':
   assert len(args)==1 and types[args[0]]=='parameter_scalar'and n['certificate']=='CHART_PEER';kind='parameter_scalar'
  elif op=='parameter_times_RHS':
   assert len(args)==2 and [types[j]for j in args]==['parameter_scalar','RHS'];kind='RHS'
  else:assert op=='RHS_subtract'and len(args)==2 and all(types[j]=='RHS'for j in args);kind='RHS'
  types.append(kind)
 inv=[n for n in ns if n['op']=='certified_parameter_unit_inverse'];assert len(inv)==1
 inv=inv[0];np=inv['args'][0];assert ns[np]['op']=='parameter_multiply'
 pref,f11=ns[np]['args'];assert ns[f11]['op']=='external_full_SOURCE_DAG'and ns[f11]['root_index']==1202
 assert O['owners']['FULL_F11']['all_nodes']==1203 and O['owners']['FULL_F11']['not_the_Hminus3_face_polynomial']
 assert ns[pref]['op']=='parameter_multiply'
 p7,f9=ns[pref]['args'];assert ns[p7]['formula']=='Q*R'and ns[f9]['owner']=='F9'
 A=O['remaining_shapes_and_ids'];rows=A['full_row_ids'];cols=A['column_ids'];zeros=set(A['zero66_consistency_RHS'])
 roots={(i,j):k for i,j,k in O['full166x13_cleared_numerator_roots']};rhs=dict(O['full166_RHS_roots']);lams=dict(O['all166_lambda_roots']);pivotrow=dict(O['pivot_row13_numerator_nodes'])
 assert len(roots)==166*13 and len(rhs)==len(lams)==166 and len(pivotrow)==13
 assert set(roots)=={(i,j)for i in rows for j in cols}
 def ext(k,i,j):
  n=ns[k];assert n['op']=='external_current_cleared_numerator'and n['owner']=='CURRENT_M10'and(n['row'],n['column'])==(i,j)
 def hr(k,i):
  n=ns[k];assert n['op']=='external_current_RHS'and n['owner']=='CURRENT_RHS10'and n['row']==i
 for j,k in pivotrow.items():ext(k,6,j)
 hp=O['pivot_RHS_node'];hr(hp,6)
 checks=0
 for i in rows:
  if i in zeros:
   assert ns[lams[i]]['op']=='QQ_constant'and ns[lams[i]]['value']=='0';hr(rhs[i],i)
  else:
   ln=ns[lams[i]];assert ln['op']=='parameter_multiply'and ln['args'][1]==inv['id'];ext(ln['args'][0],i,53)
   hn=ns[rhs[i]];assert hn['op']=='RHS_subtract';h,t=hn['args'];hr(h,i)
   assert ns[t]['op']=='parameter_times_RHS'and ns[t]['args']==[lams[i],hp]
  for j in cols:
   k=roots[i,j];assert types[k]=='parameter_scalar'
   if i in zeros:assert ns[k]['op']=='QQ_constant'and ns[k]['value']=='0'
   else:
    n=ns[k];assert n['op']=='parameter_subtract';a,b=n['args'];assert ns[a]['op']==ns[b]['op']=='parameter_multiply'
    assert ns[a]['args'][0]==np;ext(ns[a]['args'][1],i,j)
    assert ns[b]['args'][0]==ns[lams[i]]['args'][0]and ns[b]['args'][1]==pivotrow[j]
   checks+=1
 assert all(types[k]=='RHS'for k in rhs.values())
 assert len(A['extra87_RHS'])==87 and len(zeros)==66 and len(A['top_row_ids'])==100 and len(A['selected_row_ids'])==13
 pins['code/s61_col53_typed11_source_custody_20261105_v64.py']=sha(__file__)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 d={'id':ID,'all_passed':True,'structural_SOURCE_and_byte_custody_only':True,'not_independent_agent_or_scalar_VALUE_replay':True,'nodes_typed_acyclic':len(ns),'all_current_numerator_and_pivot_row_refs_checked':checks,'RHS_and_lambda_checked':len(rhs),'exact_full_F11_not_face':True,'ONE_parameter_unit_inverse_CERTIFIED':True,'source_INPUT_RAW_bytes':sum((ROOT/f).stat().st_size for f in pins),'all_history_preserved':True,'no_matrix_VALUE':True,'input_sha256':pins}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'));print('nodes/cells/RHS',len(ns),checks,len(rhs))
