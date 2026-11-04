#!/usr/bin/env python3
"""SOURCE ONLY: one correct raw10 candidate retaining all four pair source rows."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K120_T14_SOURCE_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
BP=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
@lru_cache(None)
def lattice(box,lo,hi):return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
def atom(p):
 if not p:return None
 D=math.lcm(*(Q(c).denominator for m,c in p));L=sum(abs(Q(c).numerator)*(D//Q(c).denominator)for m,c in p)
 return dict(support=len(p),box=[[min(m[k]for m,c in p),max(m[k]for m,c in p)]for k in range(3)],degree=[min(sum(m)for m,c in p),max(sum(m)for m,c in p)],D=D,L=L,Bpower=0)
ONE=atom([[[0,0,0],'1']]);B=atom(BP)
def size(a):return a['support']if a else 0
def shift(a,k):
 if a is None:return None
 return dict(a,box=[[x+k for x in a['box'][0]]]+a['box'][1:],degree=[x+k for x in a['degree']])
def mul(a,b):
 if a is None or b is None:return None
 box=[[a['box'][k][0]+b['box'][k][0],a['box'][k][1]+b['box'][k][1]]for k in range(3)];lo=a['degree'][0]+b['degree'][0];hi=a['degree'][1]+b['degree'][1];L=a['L']*b['L'];return dict(support=min(a['support']*b['support'],lattice(tuple(map(tuple,box)),lo,hi),L),box=box,degree=[lo,hi],D=a['D']*b['D'],L=L,Bpower=a['Bpower']+b['Bpower'])
def timesB(a,k):
 if a is None:return None
 bp=a['Bpower'];a=dict(a,Bpower=0)
 for _ in range(k):a=mul(a,B)
 return dict(a,Bpower=bp)
def add(a,b):
 if a is None:return b
 if b is None:return a
 bp=max(a['Bpower'],b['Bpower']);a=timesB(a,bp-a['Bpower']);b=timesB(b,bp-b['Bpower']);box=[[min(a['box'][k][0],b['box'][k][0]),max(a['box'][k][1],b['box'][k][1])]for k in range(3)];lo=min(a['degree'][0],b['degree'][0]);hi=max(a['degree'][1],b['degree'][1]);D=math.lcm(a['D'],b['D']);L=a['L']*(D//a['D'])+b['L']*(D//b['D']);return dict(support=min(a['support']+b['support'],lattice(tuple(map(tuple,box)),lo,hi),L),box=box,degree=[lo,hi],D=D,L=L,Bpower=bp)
def scale(a,c):return mul(a,atom([[[0,0,0],str(Q(c))]]))
def enc(a):
 if a is None:return dict(zero=True,support=0)
 return dict(support=a['support'],box=a['box'],degree=a['degree'],QQ_denominator_bits=a['D'].bit_length(),integer_l1_bits=a['L'].bit_length(),known_B_denominator_power=a['Bpower'])
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');h9=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');fx=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');fixed={r['old_DAG_root']:r['Laurent_polynomial']for r in fx['node_values']};bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];ke=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63','source-binding.json');cols=q['new22_target_original32_columns'];cache={};inventory=[]
 def imp(i):
  if i in cache:return cache[i]
  a=q['scalar_DAG']['nodes'][i];op=a[0]
  if i in fixed:v=atom(fixed[i])
  elif op=='z':v=None
  elif op=='one':v=ONE
  elif op=='m':v=mul(imp(a[1]),imp(a[2]))
  elif op in['a','s']:v=None if op=='s'and a[1]==a[2]else add(imp(a[1]),imp(a[2]))
  else:
   num,den=cp['matrix25x32'][a[1]][a[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   if len(den)==1:assert den[0][0][1:]==[0,0];unit=Q(den[0][1]);v=scale(shift(atom(num),-den[0][0][0]),1/unit);bp=0
   else:
    dd={tuple(m):Q(c)for m,c in den};bb={tuple(m):Q(c)for m,c in BP};assert dd.keys()==bb.keys();unit=dd[next(iter(bb))]/bb[next(iter(bb))];assert unit and all(dd[m]==unit*bb[m]for m in bb);v=scale(atom(num),1/unit);v=dict(v,Bpower=1)if v else None;bp=1
   inventory.append(dict(old_root=i,denominator_known_B_power=bp))
  cache[i]=v;return v
 def row(r,labels):return[imp(q['remaining15x22_scalar_roots'][r][cols.index(j)])for j in labels]
 pure=[22,11,13,25,15,2,19,14];labels=pure+[27,31];six=[22,11,13,25,15,19];matrix=[]
 for r in[0,1,3,6]:
  rr=row(r,pure)
  if r==0:rr=[atom(h9['source_h9'][h9['original_target_labels'].index(j)])if j in h9['original_target_labels']else x for j,x in zip(pure,rr)]
  matrix.append(rr+[None,None])
 for k,r in enumerate([2,4,5,8]):
  rr=row(r,pure)
  for j in six:rr[pure.index(j)]=atom(bot[k][six.index(j)])
  # Saved extra2 column was B-cleared. Retain its precise B denominator.
  cle=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63',f"root{ke['cleared_column_roots'][k]}.json")['polynomial'];rr[pure.index(2)]=dict(atom(cle),Bpower=1)
  comp=[1,3,4,7][k];c=Q(hc['known_row_units_rho14'][comp]['constant']);matrix.append(rr+[scale(atom(hc['raw_compact14x5'][comp][j]),c)for j in[0,4]])
 rho=hc['known_row_units_rho14'];c6=Q(rho[6]['constant']);c10=Q(rho[10]['constant']);c11=Q(rho[11]['constant']);assert rho[11]==rho[9]
 g=atom(hc['raw_compact14x5'][6][1]);h=atom(hc['raw_compact14x5'][10][1]);dif={tuple(m):Q(c)for m,c in hc['raw_compact14x5'][11][1]}
 for m,c in hc['raw_compact14x5'][9][1]:dif[tuple(m)]=dif.get(tuple(m),Q(0))-Q(c)
 hp=[[list(m),str(c)]for m,c in sorted(dif.items())if c];assert len(hp)==35;hh=atom(hp)
 r85=row(7,pure);r115=row(11,pure);r120=row(12,pure);r99=row(10,pure)
 K115=[add(scale(mul(g,x),1/c10),scale(mul(h,y),-1/c6))for x,y in zip(r115,r85)];K120=[add(scale(mul(g,add(x,scale(y,-1))),1/c11),scale(mul(hh,z),-1/c6))for x,y,z in zip(r120,r99,r85)]
 H=hc['raw_compact14x5'];K115H=[add(mul(g,atom(H[10][j])),scale(mul(h,atom(H[6][j])),-1))for j in[0,4]];K120H=[add(mul(g,add(atom(H[11][j]),scale(atom(H[9][j]),-1))),scale(mul(hh,atom(H[6][j])),-1))for j in[0,4]];matrix+=[K115+K115H,K120+K120H]
 # Only known FULLU column clears; no p/R/h22 or unknown coefficient inverse.
 bp=[max((r[j]['Bpower']if r[j]else 0)for r in matrix)for j in range(10)];matrix=[[dict(timesB(x,bp[j]-x['Bpower']),Bpower=0)if x else None for j,x in enumerate(r)]for r in matrix];rowa=[min(x['box'][0][0]for x in r if x)for r in matrix];matrix=[[shift(x,-a)for x in r]for r,a in zip(matrix,rowa)]
 # Determinant support propagation only. No coefficient multiplication occurs.
 dp={0:ONE};maxstates=0
 for i,r in enumerate(matrix):
  nex={}
  for mask,x in dp.items():
   for j,y in enumerate(r):
    if not(mask>>j&1)and y:nex[mask|1<<j]=add(nex.get(mask|1<<j),mul(x,y))
  dp=nex;maxstates=max(maxstates,len(dp))
  dump(out/'phase.json',dict(phase='source-only determinant support propagation',completed_rows=i+1,states=len(dp)))
 result=dp[1023];dense=lattice(tuple(map(tuple,result['box'])),*result['degree']);inputs=sum(size(x)for r in matrix for x in r);minimum_resident=inputs+dense+15000
 cert=dict(matrix10x10_bounds=[[enc(x)for x in r]for r in matrix],output_bound=enc(result),matrix_input_support_sum=inputs,output_lattice_slots=dense,minimum_input_output_work_reserve=minimum_resident,
  lastuse_height_scratch_fullworker_NOT_evaluated=True,no_worker_gate_pass_claim=True,output_and_inputs_lowerbound_already_exceeds150k=minimum_resident>150000,
  source_raw_heavy_columns_divided_only_FULLU_L_and_Lgamma=True,known_B_column_clearings=bp,known_a_row_factors=rowa,
  K120_h_literal_existing_coefficient_difference=hp,all_other_parameter_coefficients_unevaluated=True,source_atom_denominator_inventory=inventory)
 dump(out/'cost.json',cert)
 recipe=dict(current15_rows='P4[0,1,3,6],pair4[2,4,5,8],K115small,K120small',original32_columns=labels,
  K115small='g*current11/c10-h*current7/c6',K120small='g*(current12-current10)/c11-(N11,1-N9,1)*current7/c6',
  raw_H3zero='All selected rows kill middleheavy3, after B3 row subtraction; pair4 andP4 have literalzero there. K120 follows actual Cramer[lambda442,0,Delta] and equalrho11=rho9.',
  append_B3_pivots_current15=[7,13,10],append_B3_targets=[28,30,29],selected_rowtransform_factor='g^2/(c10*c11), g alreadyunit onD(N6)',
  actual_selected13minor_recipe='det(raw10 afterliteral source/knownunitcolumnnormalization) timesactual B3det; no claimitsvalue/nonzero',
  why_bare7_is_NOT_used='The six purepair columnsrank1 isproved only on[22,11,13,25,15,19]. Target14 isunknown andextra2isnonzeroκ direction; retainALL4pairrows andALL8tail effects.',
  scope='ONE remainingK120 row/target14 source-only. No scalar/minor/rank/point/prime/GB. No incorrect rank1 extrapolation.')
 dump(out/'recipe.json',recipe);dump(out/'phase.json',dict(phase='FINAL source-only; no coefficient worker launched'));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));sizeout=sum(p.stat().st_size for p in out.iterdir());assert sizeout<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,one_remaining_source_current15=12,one_target_original32=14,source_DAG_nodes_profiled=len(cache),matrix_input_support_sum=inputs,determinant_output_lattice_slots=dense,output_totaldegree=result['degree'],minimum_input_output_reserve=minimum_resident,wholeworker_gate_NOT_passed=True,no_scalar_job=True,all_four_pair_rows_and_tail2_retained=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
