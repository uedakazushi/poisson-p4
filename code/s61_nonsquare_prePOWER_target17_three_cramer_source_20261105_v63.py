#!/usr/bin/env python3
"""SOURCE ONLY three fixed Cramer numerators, with full pair/source normalization."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def lat(b,d):return sum(max(0,min(b[2][1],d[1]-a-e)-max(b[2][0],d[0]-a-e)+1)for a in range(b[0][0],b[0][1]+1)for e in range(b[1][0],b[1][1]+1))
def shape(poly,scale=Q(1)):
 v=[(m,Q(c)*scale)for m,c in poly if Q(c)];D=math.lcm(*(c.denominator for m,c in v));L=sum(abs(c.numerator)*(D//c.denominator)for m,c in v)
 return dict(support=len(v),box=[[min(m[k]for m,c in v),max(m[k]for m,c in v)]for k in range(3)],degree=[min(sum(m)for m,c in v),max(sum(m)for m,c in v)],D=D,L1=L)
def comb(a,b,mul=False):
 if mul:
  z=[[x[0]+y[0],x[1]+y[1]]for x,y in zip(a['box'],b['box'])];d=[a['degree'][0]+b['degree'][0],a['degree'][1]+b['degree'][1]];D=a['D']*b['D'];L=a['L1']*b['L1'];s=a['support']*b['support']
 else:
  z=[[min(x[0],y[0]),max(x[1],y[1])]for x,y in zip(a['box'],b['box'])];d=[min(a['degree'][0],b['degree'][0]),max(a['degree'][1],b['degree'][1])];D=math.lcm(a['D'],b['D']);L=a['L1']*(D//a['D'])+b['L1']*(D//b['D']);s=a['support']+b['support']
 return dict(support=min(s,lat(z,d)),box=z,degree=d,D=D,L1=L)
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,fn='report.json'):
  p=ROOT/'reports'/st/fn;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 en='S61_NONSQUARE_PREPOWER_KAPPA17_ENTRY_20261105_v63';r=load(en);assert r['all_passed']and r['kappa17_is_zero']
 hst='S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63';h=load(hst,'certificate.json')
 botst='S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63';bot=load(botst,'native.json')['bottom4x6']
 k11=load('S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63','operator.json');assert k11['identities'][0]=='kappa3*detA8=cprod*Saux*T'
 assert load('S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63_REVIEW')['all_passed']
 staged=load('S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63','operator.json');assert staged['nodes'][22]['last_row_equals']=='minus kappa'
 for st in ['S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW','S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63_REVIEW','S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63_REVIEW']:assert load(st)['all_passed']
 A=[];V=[];sources=[]
 for i,compact in enumerate([1,3,4]):
  rho=h['known_row_units_rho14'][compact];assert rho['a_d_F_B_C_J_powers']==[1,0,1,0,2,0];c=Q(rho['constant'])
  row=[]
  for col in [0,4]:
   row.append(shape(h['raw_compact14x5'][compact][col],c));sources.append(dict(kind='A3heavy',row=i,col=0 if col==0 else 1,source=f'reports/{hst}/certificate.json#/raw_compact14x5/{compact}/{col}',QQ_scale=str(c)))
  row.append(shape(bot[i][0]));sources.append(dict(kind='A3anchor',row=i,col=2,source=f'reports/{botst}/native.json#/bottom4x6/{i}/0',QQ_scale='1'));A.append(row)
  v=load(en,f'source-root{i}.json');V.append(shape([[m,str(Q(int(z),int(v['denominator'])))]for m,z in v['integer_coefficients']]));sources.append(dict(kind='v17',row=i,source=f'reports/{en}/source-root{i}.json',QQ_scale='1'))
 assert [x['support']for x in V]==[314,314,300]
 original_input_terms=sum(x['support']for row in A for x in row)+sum(x['support']for x in V)
 steps=[];outputs=[];held=0;peak=0;maxbits=0
 # Each minor/cofactor is produced once for one Cramer term, consumed and
 # released; no nine-cofactor cache or matrix inverse is formed.
 for j in range(3):
  acc=None
  for i in range(3):
   rr=[k for k in range(3)if k!=i];cc=[k for k in range(3)if k!=j]
   p=comb(A[rr[0]][cc[0]],A[rr[1]][cc[1]],True);q=comb(A[rr[0]][cc[1]],A[rr[1]][cc[0]],True);co=comb(p,q);term=comb(co,V[i],True);na=term if acc is None else comb(acc,term)
   phase=2*original_input_terms+held+(0 if acc is None else acc['support'])+p['support']+q['support']+co['support']+term['support']+na['support']
   # Both sparse keys and coefficient cells are counted. This deliberately
   # permits a conservative two-slot sparse representation in a future job.
   peak=max(peak,2*phase+15000);maxbits=max(maxbits,*[x['L1'].bit_length()+x['D'].bit_length()for x in [p,q,co,term,na]])
   steps.append(dict(output=j,term=i,cofactor_sign=(-1)**(i+j),minor_rows=rr,minor_columns=cc,first_product=p,second_product=q,cofactor=co,term_product=term,accumulator=na,sparse_all_live_owned_upper=2*phase+15000));acc=na
  outputs.append(acc);held+=acc['support']
 raw=2*1024*1024+sum(x['support']*((maxbits*30103+99999)//100000+100)for x in outputs)+original_input_terms*200
 external=sum((ROOT/p).stat().st_size for p in pins);forecast=dict(original_input_terms=original_input_terms,output_bounds=outputs,all_sparse_coefficients_AND_keys_owned_upper=peak,temporary_integer_bits_upper=2*maxbits+16,TOTAL_RAW_upper=raw+external,
  slot_gate_passes=peak<=150000,RAW_gate_passes=raw+external<32*1024*1024,actual_polynomial_job_NOT_run=True)
 proof=dict(detA3='-kappa3',inverse_kappa3='detA8/(cprod*Saux*T)',cprod='-5120000000000000/45927',no_new_chart=True,
  cramer_formula='Y_j=sum_i Cofactor(A3)_(i,j)*v17_i; x_first3=-Y/kappa3, x4=0',fourth_row_recovery='kappa dot v17=0; kappa annihilates first3 A4 columns; kappa3 unit implies fourth-row residual zero',
  ten_columns_with_x4_zero=[11,13,25,15,28,30,29,19,14,18],reason='seven fixed rawpair8 proportional columns plus three literal H3-zero columns; nonzero anchor polynomial cancels only in a polynomial identity, never inverted as a new chart',
  full_source_normalization='current rows[2,4,5] map compact[1,3,4],rho=c_i*aFC²; normalized H27/H31 remove L/Lgamma; anchor22 scale1',
  independent_kappa17_peer_pending=True,old4plus10_and_all_tail_effects_retained=True,no_actual_Cramer_or_lastrow_value=True)
 dump(out/'source-and-unit-proof.json',proof);dump(out/'three-cramer-schedule.json',dict(A3_bounds=A,v17_bounds=V,literal_sources=sources,steps=steps,forecast=forecast))
 code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 rep=dict(id=ID,all_passed=True,frozen=True,source_only=True,all_three_pair_corrections_retained=True,no_new_unit_chart=True,conditional_on_kappa17_independent_zero_peer=True,
  no_new_Cramer_or_lastrow_scalar_value=True,prospective_slot_gate_passes=forecast['slot_gate_passes'],prospective_RAW_gate_passes=forecast['RAW_gate_passes'],forecast=forecast,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
