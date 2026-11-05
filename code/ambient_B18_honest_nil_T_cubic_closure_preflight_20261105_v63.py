#!/usr/bin/env python3
"""Known cubic-input construction only, no cubic rank/GB/source computation.
Reduce the600x1284 known closure to an exactly equivalent24x1308 input.
"""
import collections,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63'
PATHS={
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'first_peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW/report.json',
 'quadbasis':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'actual_quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63/phase3-all661-actual-six-quadratic-values-and-rank.json',
 'actual_quad_report':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_VALUES_20261105_v63/report.json',
 'length11':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'Tnil_length':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def put(d,k,c):
 if c:
  d[k]=d.get(k,Q(0))+c
  if not d[k]:del d[k]
def sp(d):return [[i,c]for i,c in sorted(d.items())if c]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists()
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert D['first_peer']['all_passed']and D['actual_quad_report']['all_passed']and D['length11']['all_passed']and D['Tnil_length']['all_passed']
 A,F,B=D['actual_quad'],D['first'],D['quadbasis'];assert A['actual_operator_rank']==3 and F['actual_first_relation_rank']==35 and B['known_quadratic_rank']==174
 m2=list(itertools.combinations_with_replacement(range(8),2));m3=list(itertools.combinations_with_replacement(range(8),3));i2={m:i for i,m in enumerate(m2)};i3={m:i for i,m in enumerate(m3)};assert len(m2)==36 and len(m3)==120
 pi6=[]
 for row in B['exact_quotient_duals_rows_sparse']:pi6.append(dict((i,Q(c))for i,c in row))
 nf6=[36*target+i2[tuple(mon)]for target,mon in B['quotient_coordinate_labels']]
 assert [[pi6[i].get(j,Q(0))for j in nf6]for i in range(6)]==[[Q(i==j)for j in range(6)]for i in range(6)]
 cp=A['fixed_rank_columns'];pr=A['fixed_rank_rows'];Jinv=[[Q(c)for c in row]for row in A['fixed_minor_inverse']];Qcols=[[Q(c)for c in A['all661_actual_quadratic_six_values'][j]]for j in cp];free=[i for i in range(6)if i not in pr];assert free==[2,4,5]
 p36=[]
 for f in free:
  row=[Q(i==f)for i in range(6)]
  for k,p in enumerate(pr):row[p]-=sum((Qcols[j][f]*Jinv[j][k]for j in range(3)),Q(0))
  p36.append(row)
 assert all(sum((p36[i][j]*Qcols[k][j]for j in range(6)),Q(0))==0 for i in range(3)for k in range(3))
 pi3=[]
 for row in p36:
  v={}
  for j,c in enumerate(row):
   if c:
    for k,a in pi6[j].items():put(v,k,c*a)
  pi3.append(v)
 nf3=[nf6[i]for i in free]
 assert [[pi3[i].get(j,Q(0))for j in nf3]for i in range(3)]==[[Q(i==j)for j in range(3)]for i in range(3)]
 lin=[dict((i,Q(c))for i,c in F['all849_actual_minimal_first_columns40_sparse'][j])for j in F['rank_minor_columns']]
 quadrics=[{nf6[i]:c for i,c in enumerate(col)if c}for col in Qcols]
 assert sum(map(len,lin))==197 and sum(map(len,quadrics))==9
 raw=[];labels=[]
 for j,col in enumerate(lin):
  for mon in m2:
   v={}
   for row,c in col.items():target,var=divmod(row,8);put(v,120*target+i3[tuple(sorted(mon+(var,)))],c)
   raw.append(v);labels.append(['first35',F['rank_minor_columns'][j],list(mon)])
 for j,col in enumerate(quadrics):
  for var in range(8):
   v={}
   for row,c in col.items():target,mi=divmod(row,36);put(v,120*target+i3[tuple(sorted(m2[mi]+(var,)))],c)
   raw.append(v);labels.append(['quadratic3',cp[j],var])
 assert len(raw)==1284 and sum(map(len,raw))==7164
 # q:F3->V3 tensor S1. Canonical divisor is the first two sorted variables.
 qcols=[];divisor=[]
 for row in range(600):
  target,mi=divmod(row,120);mon=m3[mi];mon2=mon[:2];var=mon[2];k=36*target+i2[mon2]
  qcols.append({8*f+var:pi3[f][k]for f in range(3)if k in pi3[f]});divisor.append([target,list(mon2),var])
 def project(col):
  v={}
  for row,c in col.items():
   for j,a in qcols[row].items():put(v,j,c*a)
  return v
 reduced=[project(c)for c in raw];E=[];cons=[]
 for f,quadrow in enumerate(nf3):
  target,mi=divmod(quadrow,36)
  for var in range(8):
   row=120*target+i3[tuple(sorted(m2[mi]+(var,)))];E.append(row);v={8*f+var:Q(1)}
   for j,c in qcols[row].items():put(v,j,-c)
   cons.append(v)
 reduced+=cons;assert len(reduced)==1308
 den=[1]*24
 for col in reduced:
  for row,c in col.items():den[row]=math.lcm(den[row],c.denominator)
 integer=[{row:int(c*den[row])for row,c in col.items()}for col in reduced]
 assert all(Q(v[row],den[row])==c for col,v in zip(reduced,integer)for row,c in col.items())
 bits=[max([1]+[abs(col.get(row,0)).bit_length()for col in integer])for row in range(24)]
 detbits=sum(bits)+math.factorial(24).bit_length()+1
 # Future rank job stores minor/inverse only; all-column factors stream to a hash.
 inverse_raw_upper=576*(2*((detbits*30103+99999)//100000)+32)
 terms={'source600x1284':sum(map(len,raw)),'reduced24x1308':sum(map(len,reduced)),'projection600':sum(map(len,qcols)),'linear_and_quadratic_input':206,'future_basis_and_inverse':1152,'future_one_column_and_factor':48,'misc_receipts':3000}
 aggregate=sum(terms.values());assert aggregate<150000
 data={'same_actual_S_variable_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],'source_actual_HF_through2':[5,5,3],
  'all35_actual_linear_relations40_sparse':[sp(v)for v in lin],'all3_actual_quadratic_NF_relations180_sparse':[sp(v)for v in quadrics],
  'actual_quadratic_source_columns661':cp,'quadratic_NF_basis180_indices':nf3,'quadratic_projection3_by180_sparse':[sp(v)for v in pi3],
  'all1284_known_cubic_columns600_sparse':[sp(v)for v in raw],'all1284_source_labels':labels,
  'canonical_q24_by600_columns_sparse':[sp(v)for v in qcols],'all600_divisor_witnesses_for_Id_minus_Eq_in_known_cubic':divisor,
  'E600_by24_unit_column_rows':E,'all24_consistency_columns_Id_minus_qE_sparse':[sp(v)for v in cons],
  'equivalent_known_cubic24_by1308_sparse_columns':[sp(v)for v in reduced],'row_denominator_units24':den,'equivalent_integer24_by1308_sparse_columns':[sp(v)for v in integer],
  'all_row_scaling_source_identities_verified':True,'complete_known_degree3_source_count':1284,
  'equivalence_proof':'pi3 U3=I and ker(pi3)=known degree2 span177. For each cubic monomial e choose divisor e2*x: e-Eq(e)=x*(e2-U3pi3(e2)) belongs known degree3 K. Therefore F600/K is isomorphic F24/(qK+(I-qE)F24), retaining all1284 source paths and24 consistency relations.',
  'degree3_completeness_scope':'all multiples of the35 actual linear and3 actual quadratic NF relations; NEW actual cubic relations may remain',
  'known_cubic_rank_or_actual_HF3_not_computed':True,'pending_cubic_rank_formula':'known upper for actual HF3 is24-rank(integer24x1308)',
  'budget_terms':terms,'aggregate_QQ_slot_upper':aggregate,'integer_row_entry_bitbounds':bits,'any24minor_bitbound':detbits,'future_inverse_raw_byte_upper':inverse_raw_upper,
  'future_factorizations_must_stream_hash_not_save_dense24x1308':True,
  'conditional_total_length16_tail_options':[[5,5,3,3],[5,5,3,2,1],[5,5,3,1,1,1]],
  'tail_options_proof':'length16 and first553 leave3. HF3 cannot0 by Nakayama; if HF3=1 every later value<=1 by the unique-standard-monomial divisor lemma; HF3=2 leaves one further class; HF3=3 exhausts length.'}
 p=out/'literal-known-cubic-input-and-exact24-row-reduction.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())+inverse_raw_upper<32<<20;p.write_text(blob)
 r={'id':STEM,'all_passed':True,'status':'SOURCE_INPUT_AND_BUDGET_ONLY_no_cubic_rank_run','original_known_cubic_shape':[600,1284],'equivalent_shape':[24,1308],'all1284_paths_and24_consistency_relations_retained':True,'aggregate_QQ_upper':aggregate,'literal_input_raw_bytes':len(blob.encode()),'future_inverse_raw_byte_upper':inverse_raw_upper,'cubic_rank_or_HF3_computed':False,'new_source_gauge_jet':False,'quadratic_independent_peer_required_before_rank_job':True,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
