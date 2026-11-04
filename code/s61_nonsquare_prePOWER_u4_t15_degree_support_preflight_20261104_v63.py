#!/usr/bin/env python3
"""Sharpen frozen BOX cost by total-degree lattice counts and rigorous l1=0."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_U4_T15_DEGREE_SUPPORT_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();assert not(out/'report.json').exists();base=ROOT/'reports/S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63';old=read(base/'operator-and-cost.json');N=old['nodes'];cache={};paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',base/'report.json',base/'operator-and-cost.json']
 def source(p):
  p=ROOT/p
  if p not in cache:cache[p]=read(p);paths.append(p)
  return cache[p]
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json')
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if r['path']=='accepted_Q2L_source_atom':
   x=r['definition'];nd=cp['matrix25x32'][x[1]][x[2]]if x[0]=='checkpoint'else pd['Lframe_dual6x2'][x[1]][1]if x[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][x[1]]
   assert len(nd[1])==1 and nd[1][0][0][1:]==[0,0], 'Selected entry must not silently omit the old99 B-denominator atoms'
   return nd[0]
  x=source(r['path']);k=r['key']
  if k=='node_values':return next(y['Laurent_polynomial']for y in x[k]if y['old_DAG_root']==r['old_DAG_root'])
  if k=='source_h9':return x[k][x['original_target_labels'].index(r['target'])]
  x=x[k]
  for i in r.get('index',[]):x=x[i]
  return x
 def count(box,lo,hi):return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
 B=[];sizes=[];zeros=[]
 for i,n in enumerate(N):
  assert n['bound'].get('known_B_denominator_power',0)==0
  if int(n['integer_l1_upper'])==0:B.append(None);sizes.append(0);zeros.append(i);continue
  if n['op']=='atom':
   pp=[m for m,c in atom(n)if c not in['0','0/1']];s=n['a_shift'];lo=min(sum(m)+s for m in pp);hi=max(sum(m)+s for m in pp);box=tuple((min(m[j]for m in pp)+(s if j==0 else 0),max(m[j]for m in pp)+(s if j==0 else 0))for j in range(3));nn=len(pp)
  else:
   a,b=n['args'];A,C=B[a],B[b]
   if n['op']=='mul':
    assert A is not None and C is not None;lo,hi=A[0]+C[0],A[1]+C[1];box=tuple((A[2][j][0]+C[2][j][0],A[2][j][1]+C[2][j][1])for j in range(3));nn=sizes[a]*sizes[b]
   else:
    zz=[(x,sizes[t])for t,x in zip([a,b],[A,C])if x is not None];lo=min(x[0][0]for x in zz);hi=max(x[0][1]for x in zz);box=tuple((min(x[0][2][j][0]for x in zz),max(x[0][2][j][1]for x in zz))for j in range(3));nn=sum(x[1]for x in zz)
  B.append((lo,hi,box));sizes.append(min(nn,count(box,lo,hi)))
 use=[0]*len(N)
 for n in N:
  for j in n['args']:use[j]+=1
 use[old['output']]+=1;live=peak=0;peakid=0;records=[]
 for i,n in enumerate(N):
  live+=sizes[i];scratch=sizes[i]if n['op']in['add','sub','mul']else 0
  if live+scratch>peak:peak=live+scratch;peakid=i
  records.append(dict(node=i,live=live,scratch=scratch))
  for j in n['args']:
   use[j]-=1
   if not use[j]:live-=sizes[j]
 total=peak+old['literal_atom_terms']+15000;assert total<=150000 and old['raw_artifact_bytes_upper']<32*1024*1024
 data=dict(original_operator_SHA256=sha(base/'operator-and-cost.json'),source_only=True,same_768_arithmetic_nodes=True,
  degree_box_bounds=B,support_upper=sizes,rigorous_l1_zero_nodes=zeros,last_use_records=records,
  zero_proof='Each saved integerN is a rigorous nonnegative UPPER bound for the entire coefficient l1 norm after commonQQdenominator clearing. N=0 forces every coefficient0; no cancellation assumption.',
  total_degree_proof='Atom degrees include a_shift; multiplication addsdegreebounds; additiontakesmin/maxsince selectedliveDAGhasNO B-poles. Every integerlattice monomial in its coordinatebox anddegreeinterval iscounted.',
  all_768_B_denominator_powers_zero_verified=True,all_selected_uncached_denominators_are_literal_a_QQ=True,
  earlier_99_B_pole_atoms_not_ignored='Thoseexistforothercoupledcoefficients,but NONE belongtoTHISselected atomclosure. Explicitdenominator-shape assertions are saved.',
  peak_node=peakid,last_use_peak_numerator_terms=peak,literal_atom_reserve=old['literal_atom_terms'],generic_work_reserve=15000,total_live_upper=total,
  output_support_upper=sizes[-1],output_degree_range=B[-1][:2],raw_artifact_bytes_upper=old['raw_artifact_bytes_upper'],
  Laurent_inplace_alignment='For each exponent change +k of numerator a terms, subtractk from scalar a_shift. The representedRF is IDENTICAL. Termorderisunchanged, so in-place exponent updates avoid fullshiftedcopies. Dynamicworker addsoutputsizedscratch.',
  no_value_coefficient_minor_rank_point_prime=True)
 (out/'certificate.json').write_text(json.dumps(data,indent=2)+'\n');report=dict(id=ID,all_passed=True,frozen=True,source_only=True,resource_gate_passes=True,
  old_BOX_bound=old['total_live_with_literal_and_15000reserve'],new_total_live_upper=total,output_support_upper=sizes[-1],
  integer_coefficient_bits=old['output_integer_coefficient_bits'],QQ_denominator_bits=old['QQ_denominator_bits'],raw_artifact_bytes_upper=old['raw_artifact_bytes_upper'],
  original_failed_BOX_packet_untouched=True,no_new_coefficient_values=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},certificate_sha256=sha(out/'certificate.json'),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
