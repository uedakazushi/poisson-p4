#!/usr/bin/env python3
"""Gated six-root cache composition and coherent RIGHT2 quartic rank."""
import argparse
import hashlib
import json
import re
import resource
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/report.json'
SOURCE_SHA='ae2b9c8162f53a396dcfb52288edbb53035569742b384de9b34286cf8cdc0d53'
COMPILED='reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/compiled-six-even-roots-and-provided-joint105-plane3-cache.json'
PROFILE='reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/lambda3-cache-and-right2-whole-source-profile.json'
Q20='reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json'
OPERATOR='reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json'
MATH_PEER='reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json'
MATH_PEER_SHA='05d8c4eb681961627e8d3d6210856052b52a6c6945efebfc055e2621b0b975af'
GATE='reports/AMBIENT_B18_C3_RIGHT2_WHOLE_VALUE_GATE_20261105_v63/report.json'
GATE_SHA='9fc7c8fffd965ebd5fc0cbaac7399d203e960b671fb69707c802849295f1a0bb'
C2='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
MON=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)]
MI={m:i for i,m in enumerate(MON)}
PAIRS=[[(a,b)for a,u in enumerate(MON)for b,v in enumerate(MON)
        if(u[0]+v[0],u[1]+v[1])==m]for m in MON]
MAX_PRECHECK=0
MAX_VALUE=0


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def field(text,key):
 token=json.dumps(key);assert text.count(token)==1;i=text.index(':',text.index(token))+1
 while text[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(text,i)[0]
def checked_height(n):
 global MAX_PRECHECK
 MAX_PRECHECK=max(MAX_PRECHECK,n);assert n<=2048
def note(q):
 global MAX_VALUE
 MAX_VALUE=max(MAX_VALUE,abs(q.numerator).bit_length(),q.denominator.bit_length())
 assert MAX_VALUE<=2048
 return q
def parse(s):
 p=str(s).split('/');checked_height(max(abs(int(v)).bit_length()for v in p))
 return note(Q(s))
def add(a,b):
 if not a:return b
 if not b:return a
 checked_height(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),
                         abs(b.numerator).bit_length()+a.denominator.bit_length())+1,
                    a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def times(a,b):
 if not a or not b:return Q(0)
 checked_height(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),
                    a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def neg(a):return note(-a)
def conv(a,b):
 v=[Q(0)]*6
 for m,pairs in enumerate(PAIRS):
  for i,j in pairs:v[m]=add(v[m],times(a[i],b[j]))
 return v
def inverse(x):
 assert x[0]==1
 v=[Q(1)]+[Q(0)]*5
 for m in range(1,6):
  for a,b in PAIRS[m]:
   if a:v[m]=add(v[m],neg(times(x[a],v[b])))
 return v
def divide(a,b):
 assert b
 checked_height(max(abs(b.numerator).bit_length(),b.denominator.bit_length()))
 return times(a,note(Q(b.denominator,b.numerator)))
def mm(a,b):
 result=[]
 for row in a:
  out=[]
  for col in zip(*b):
   q=Q(0)
   for x,y in zip(row,col):q=add(q,times(x,y))
   out.append(q)
  result.append(out)
 return result
def eval_form(terms,lam,cached):
 q=Q(0);records=[]
 for key,w in terms:
  tag,physical,index=key
  if tag=='lambda':v=lam[physical]
  elif tag=='cache':v=parse(cached[physical].get(index,'0'))
  elif tag=='constant':v=Q(1)
  else:raise AssertionError(tag)
  weight=parse(w);p=times(weight,v);q=add(q,p)
  records.append([key,str(weight),str(v),str(p)])
 return q,records
def quartic_result(operator,gate,lam,cached):
 matrix=[];terms=[]
 for row in operator['matrix2_rows_A4_B4_columns_HB0_HB1_linear_leaf_operators']:
  rr=[];tt=[]
  for form in row:
   value,products=eval_form(form,lam,cached);rr.append(value);tt.append(products)
  matrix.append(rr);terms.append(tt)
 a0,a1=matrix[0];b0,b1=matrix[1]
 c0=add(b0,times(parse('1/4'),a0));c1=add(b1,times(parse('1/4'),a1))
 f=add(times(parse('10'),a0),a1)
 factors=[eval_form(form,lam,cached)for form in gate['two_affine_factor_forms']]
 assert c0==factors[0][0]and f==factors[1][0]
 assert add(c1,times(parse('10'),c0))==0
 determinant=add(times(a0,b1),neg(times(a1,b0)))
 assert determinant==neg(times(c0,f))
 if determinant:
  rank=2;iv=[[divide(b1,determinant),neg(divide(a1,determinant))],
             [neg(divide(b0,determinant)),divide(a0,determinant)]]
  eye=[[Q(1),Q(0)],[Q(0),Q(1)]]
  left,right=mm(matrix,iv),mm(iv,matrix);assert left==eye==right
  proof={'fixed_full_minor':str(determinant),'inverse':[[str(v)for v in row]for row in iv],
         'matrix_times_inverse':[[str(v)for v in row]for row in left],
         'inverse_times_matrix_and_all_two_source_factors':[[str(v)for v in row]for row in right]}
 else:
  pivot=next(((i,j)for i in range(2)for j in range(2)if matrix[i][j]),None)
  if pivot is None:rank=0;proof={'all_four_entries_zero':True}
  else:
   rank=1;i,j=pivot;column=[matrix[t][j]for t in range(2)]
   coefficients=[divide(matrix[i][k],column[i])for k in range(2)]
   assert all(times(column[t],coefficients[k])==matrix[t][k]for t in range(2)for k in range(2))
   proof={'pivot_row_column':[i,j],'pivot_column':list(map(str,column)),
          'all_two_source_column_factors':list(map(str,coefficients)),'all4_factorizations_verified':True}
 return {'actual_quartic_matrix_rows_A4_B4_columns_HB0_HB1':[[str(v)for v in row]for row in matrix],
     'all_termwise_affine_products':terms,'two_affine_factors':[str(c0),str(f)],
     'two_factor_termwise_products':[p for v,p in factors],
     'C1_plus10C0_zero':True,'determinant':str(determinant),'det_equals_minus_factor_product':True,
     'actual_RIGHT2_rank':rank,'rank_certificate':proof,
     'actual_C3_length19_proved':rank==2,'actual_C3_HF_if_rank2':[3,6,6,4,0]if rank==2 else None,
     'C_DVR_cut3_torsion_pattern_if_rank2':[1,1,2,'e>=3']if rank==2 else None,
     'rank_below2_only_source_subset_inconclusive':rank<2,
     'longer_exponent_full_torsion_genuine_C_ell_H0_or_nil_not_determined':True}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--source-peer-report',required=True)
 args=parser.parse_args();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))
 started=time.monotonic();assert sha(ROOT/SOURCE)==SOURCE_SHA
 source=json.loads((ROOT/SOURCE).read_text())
 peer_path=Path(args.source_peer_report)
 if not peer_path.is_absolute():peer_path=ROOT/peer_path
 peer=json.loads(peer_path.read_text());assert peer['all_passed']
 assert SOURCE_SHA in set(peer.get('input_sha256',{}).values())
 for p in[COMPILED,PROFILE]:assert sha(ROOT/p)==source['evidence_sha256'][Path(p).name]
 assert sha(ROOT/Q20)==source['input_sha256'][Q20]
 c=json.loads((ROOT/COMPILED).read_text());profile=json.loads((ROOT/PROFILE).read_text())
 assert sha(ROOT/GATE)==GATE_SHA and sha(ROOT/MATH_PEER)==MATH_PEER_SHA
 gate=json.loads((ROOT/GATE).read_text());math_peer=json.loads((ROOT/MATH_PEER).read_text())
 assert gate['all_passed']and math_peer['all_passed']
 assert sha(ROOT/OPERATOR)==source['input_sha256'][OPERATOR]==gate['right2_operator_sha256']
 assert sha(ROOT/C2)==source['input_sha256'][C2]==gate['C2_input_sha256']
 operator=json.loads((ROOT/OPERATOR).read_text());c2full=json.loads((ROOT/C2).read_text())
 text=(ROOT/Q20).read_text();qm=field(text,'base27_monomials');qvals=dict(field(text,'all105_cached_composite_values'));del text
 qi=[qm.index(m)for m in[[],[7],[8],[7,7],[7,8],[8,8]]]
 records=c['operator_records_originalID_op_args'];assert len(records)==7095
 input_count=count(source)+count(peer)+count(c)+count(profile)+count(gate)+count(math_peer)+count(operator)+count(c2full)
 # Includes the largest temporary selected Q20 field (2940), converted
 # inputs, dependency counts, every last-use buffer, output copies and
 # serialization/controller metadata. Dense None controls are not scalars.
 owned_bound=(input_count+2940+735+7095+984+22000)
 assert owned_bound<=75000<100000
 raw_inputs=sum(p.stat().st_size for p in[ROOT/SOURCE,ROOT/COMPILED,ROOT/PROFILE,ROOT/Q20,peer_path,ROOT/GATE,ROOT/MATH_PEER,ROOT/OPERATOR,ROOT/C2])
 assert raw_inputs+150000<=1000000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 (out/'pre-arithmetic-source-and-budget-custody.json').write_text(json.dumps({
     'source_report_sha256':SOURCE_SHA,'worker_sha256':sha(Path(__file__)),
     'source_peer_sha256':sha(peer_path),'owned_numeric_bound':owned_bound,
     'raw_input_bytes':raw_inputs,'hard_bits':2048,'hard_seconds':60,
     'root_values_not_yet_evaluated':True},indent=2)+'\n')
 phi={};overlap_count=0
 for i,v in c['provided105_originalID_r_s_coefficients10']:
  a=list(map(parse,v));lower=[a[0],add(a[1],a[2]),times(a[2],parse('1/5')),
      add(add(a[3],a[4]),a[5]),times(add(a[4],times(parse('2'),a[5])),parse('1/5')),
      times(a[5],parse('1/25'))]
  assert lower==[parse(qvals[i][k])for k in qi];overlap_count+=6
  # Ordinary Taylor coefficients under s=r+ell/5; B6 discards r^2.
  phi[i]=[lower[0],lower[1],lower[2],lower[4],lower[5],
          times(add(a[8],times(parse('3'),a[9])),parse('1/25'))]
 assert overlap_count==630 and set(phi)==set(qvals)
 del qvals,qm
 max_id=max(r[0]for r in records)
 uses=[None]*(max_id+1);cache=[None]*(max_id+1)
 for r in records:uses[r[0]]=0
 roots=[r[n]for r in c['root_pairs_physical3']for n in['aA_node','aB_node']]
 assert [r['physical_row']for r in c['root_pairs_physical3']]==[18,60,67]
 assert roots==[20114,40892,20138,40916,20153,40931]
 for r in records:
  if r[1]not in['input','c']:
   for j in r[2:]:assert uses[j]is not None;uses[j]+=1
 for i in roots:uses[i]+=1
 live=peak=0;unit_checks=[]
 try:
  for r in records:
   i,op=r[:2]
   if op=='input':v=phi[i]
   elif op=='c':v=[parse(str(r[2])+'/'+str(r[3]))]+[Q(0)]*5
   elif op=='a':v=[add(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
   elif op=='n':v=[neg(a)for a in cache[r[2]]]
   elif op=='m':v=conv(cache[r[2]],cache[r[3]])
   elif op=='i':
    v=inverse(cache[r[2]]);one=[Q(1)]+[Q(0)]*5
    assert conv(v,cache[r[2]])==one==conv(cache[r[2]],v)
    unit_checks.append([i,r[2],True])
   else:raise AssertionError(op)
   cache[i]=v;live+=1;peak=max(peak,live)
   if op not in['input','c']:
    for j in r[2:]:
     uses[j]-=1
     if not uses[j]:cache[j]=None;live-=1
   if not uses[i]:cache[i]=None;live-=1
 except BaseException as err:
  checkpoint={'status':'STOP_no_cap_raise','error':repr(err),'completed_original_node':i,
      'known_root_values':[[n,list(map(str,cache[n]))]for n in roots if cache[n]is not None],
      'max_precheck_bits':MAX_PRECHECK,'max_value_bits':MAX_VALUE}
  (out/'ROOT_ONLY_STOP_CHECKPOINT.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
  raise
 assert peak==162 and unit_checks==[[a,b,True]for a,b,z in c['both_original_inverse_centers']]
 expected=dict(c['expected6_lower_B6_coefficients_order_1_r_ell_rell_ell2'])
 for i in roots:assert cache[i][:5]==list(map(parse,expected[i]))
 actual=[]
 for r in c['root_pairs_physical3']:
  aa,ab=cache[r['aA_node']],cache[r['aB_node']]
  delta=[add(x,neg(y))for x,y in zip(aa,ab)]
  assert delta[:5]==[Q(0)]*5
  actual.append({'physical_row':r['physical_row'],'aA_node':r['aA_node'],'aB_node':r['aB_node'],
      'aA_B6':list(map(str,aa)),'aB_B6':list(map(str,ab)),
      'delta_B6':list(map(str,delta)),'delta_r_ell2':str(delta[5])})
 assert MAX_PRECHECK<=79
 root_phase={'all6_even_roots_B6':[[i,list(map(str,cache[i]))]for i in roots],
    'all3_lambda_readbacks':actual,'all630_input_overlaps_verified':True,
    'all30_prior_root_lower_overlaps_verified':True,'both_unit_inverse_products':unit_checks,
    'root_phase_max_precheck_bits':MAX_PRECHECK,'root_phase_max_value_bits':MAX_VALUE,
    'no_new_physical_source_gauge_IFT_or_jet':True}
 phase1=out/'phase1-six-root-existing-cache-readbacks.json';phase1.write_text(json.dumps(root_phase,indent=2)+'\n')
 lam={r['physical_row']:parse(r['delta_r_ell2'])for r in actual}
 cached={r['physical_row']:dict(r['delta'])for r in c2full['physical12_coefficients']}
 try:
  post=quartic_result(operator,gate,lam,cached)
  assert MAX_PRECHECK<=gate['whole_root_and_post_operation_bound']<=2048
 except BaseException as err:
  checkpoint={'status':'POSTROOT_STOP_no_readback_rerun','error':repr(err),'phase1_sha256':sha(phase1),
      'lambda_values':{str(k):str(v)for k,v in lam.items()},'max_precheck_bits':MAX_PRECHECK}
  (out/'ROOT_ONLY_POST_STOP_CHECKPOINT.json').write_text(json.dumps(checkpoint,indent=2)+'\n');raise
 post.update({'owned_numeric_bound_all_metadata_copies':owned_bound,'input_RAW_bytes':raw_inputs,
     'whole_max_precheck_bits':MAX_PRECHECK,'whole_max_value_bits':MAX_VALUE,
     'all_source_and_post_operations_prechecked_at2048':True,'new_physical_gauge_source_or_jet':False})
 ep=out/'phase2-actual-RIGHT2-quartics-two-factors-and-rank.json';ep.write_text(json.dumps(post,indent=2)+'\n')
 rank=post['actual_RIGHT2_rank']
 report={'id':STEM,'all_passed':True,'status':'RIGHT2_RANK2_EXACT_C3_LENGTH19'if rank==2 else'RIGHT2_RANK_BELOW2_SUBSET_INCONCLUSIVE',
   'actual_RIGHT2_rank':rank,'actual_C3_length19_proved':rank==2,
   'longer_C_torsion_exponent_or_original_H0_nil_not_determined':True,
   'input_sha256':{SOURCE:SOURCE_SHA,COMPILED:sha(ROOT/COMPILED),PROFILE:sha(ROOT/PROFILE),
       Q20:sha(ROOT/Q20),GATE:GATE_SHA,MATH_PEER:MATH_PEER_SHA,OPERATOR:sha(ROOT/OPERATOR),C2:sha(ROOT/C2),str(peer_path.relative_to(ROOT)):sha(peer_path)},
   'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
   'evidence_sha256':{p.name:sha(p)for p in[phase1,ep]},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 assert raw_inputs+sum(p.stat().st_size for p in out.iterdir())<=1000000
 print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()
