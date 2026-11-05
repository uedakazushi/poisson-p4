#!/usr/bin/env python3
"""Independent Laurent source parser + all100 ZZ inverse identities.
No determinant/cofactor/rank/GCD/point algorithm is called by this peer.
"""
import hashlib,json,math,os,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_FIXED_CONSTRAINT_ADJUGATES_20261105_v64_REVIEW'
OUT=ROOT/'reports'/ID;WORK=ROOT/'work'/ID
PRIMARY='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED_CONSTRAINT_ADJUGATES_20261105_v64'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_FIXED_CONSTRAINT_MINOR_PREFLIGHT_20261105_v64'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def poly(p):
 if p is None:return{}
 c=F(p['content']);sh=p.get('shift',p.get('shift5'));assert len(sh)==5;out={}
 for t in p['terms']:
  e=t[1]if len(t)==2 and isinstance(t[1],list)else t[1:];assert len(e)==5
  key=tuple(int(x)+int(y)for x,y in zip(sh,e));assert key not in out
  a=c*int(t[0]);assert a;out[key]=a
 return out
def scaled(p,rho,mon):return{tuple(a+b for a,b in zip(e,mon)):c*rho for e,c in p.items()}
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);WORK.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2)+'\n')
 rep=load(PRIMARY);assert pins[f'reports/{PRIMARY}/report.json']=='5c45544b64ce273de96ff56fbe26ae5e14ecc56c363b4d2e0b41b443f7d09b2a'and rep['all_passed']
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h,p
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/PRIMARY/n)==h,n
 source=load(SOURCE);assert pins[f'reports/{SOURCE}/report.json']=='15ec0fa2870a5985257a5a9e62500141ade5619b2ba1afe7b50a9bddcafea06e'
 for p,h in source['input_sha256'].items():assert sha(ROOT/p)==h,p
 for n,h in source['evidence_sha256'].items():assert sha(ROOT/'reports'/SOURCE/n)==h,n
 inputs=load(SOURCE,'two-actual-fixed-constraint-minors.json');assert inputs['variable_order']==['H','L','s0','q0','t0'] and inputs['original_free_coordinates']==[6,12]
 profile=load('S61_SPARSE_GENERIC5_CONSTRAINED83_EXACT49_CONSTRAINT_PROFILE_20261105_v64','actual-constraint-blocks-and-frontier.json');groups={g['name']:g for g in profile['actual_constraint_blocks']}
 old=load('S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63');assert old['all_passed']and old['kernel_complete']
 assert load('S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW')['all_passed']
 # Hash the large old kernel, but do NOT parse159249 numeric metadata.
 kp=ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63/generic-field-kernel.json'
 with kp.open('rb')as f:header=''.join(f.read(240).decode().split())
 assert '"full_L_rank":18'in header and '"free_original_columns":[6,12]'in header
 save('phase.json',dict(phase='all immutable actual source/profile and output hashes pinned; no rank or determinant replay'))
 blocks=[];checks=[];qq_terms=0
 for block in inputs['blocks']:
  name=block['name'];group=groups[name];off=134 if name=='alpha'else 185;cells={(r,c):v for r,c,k,h,v in group['entries']};rows=block['rows'];cols=block['pivot_columns'];free=block['free_columns']
  assert rows==[off+i for i in[6,14,16,41,42]]and len(cols)==5 and len(free)==2
  H=[]
  for i,row in enumerate(rows):
   literal=block['literal_H5x5'][i]+block['literal_B5x2'][i];assert literal==[cells.get((row,c))for c in cols+free]
   pr=[poly(v)for v in literal];den=1
   for p in pr:
    for c in p.values():den=math.lcm(den,c.denominator)
   mins=[min(e[j]for p in pr for e in p)for j in range(5)];scale=block['row_clearings'][i]
   assert scale['row']==row and den==int(scale['positive_integer_denominator'])and [-v for v in mins]==scale['parameter_monomial_exponent']
   for j,p in enumerate(pr):
    q=scaled(p,den,scale['parameter_monomial_exponent']);assert all(min(e)>=0 and c.denominator==1 for e,c in q.items())
    expected=poly(load(PRIMARY,f'{name}-cleared-{i}-{j}.json'));assert expected==q
    if j<5:H.append(q)
    checks.append([name,row,cols[j]if j<5 else free[j-5],len(q)])
  A=[poly(load(PRIMARY,f'{name}-adj-{i}-{j}.json'))for i in range(5)for j in range(5)];D=poly(load(PRIMARY,f'{name}-det.json'));assert D
  assert all(min(e)>=0 and c.denominator==1 for p in H+A+[D]for e,c in p.items());qq_terms+=sum(map(len,H+A+[D]));blocks.append((name,H,A,D))
 assert len(checks)==70 and qq_terms<150000
 control=WORK/'independent-five-variable-ZZ-input.txt'
 with control.open('w')as f:
  f.write('2\n')
  for name,H,A,D in blocks:
   f.write(name+'\n')
   for p in H+A+[D]:
    f.write(str(len(p))+'\n')
    for e,c in sorted(p.items()):f.write(str(c.numerator)+' '+' '.join(map(str,e))+'\n')
 cpp=ROOT/'code/s61_sparse_generic5_fixed_constraint_adjugates_review_20261105_v64.cpp';exe=WORK/'checker';native=OUT/'independent-ZZ-left-right-identities.json'
 comp=subprocess.run(['bash','-lc',f'g++ -O2 -std=c++17 {cpp} -o {exe} $(pkg-config --cflags --libs flint)'],capture_output=True,text=True,timeout=15);(OUT/'compile.log').write_text(comp.stdout+comp.stderr);assert comp.returncode==0
 run=subprocess.run([str(exe),str(control),str(native)],capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-st)));(OUT/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0
 N=json.loads(native.read_text());assert N['all_passed']and N['all_left_right_identities']==100 and N['both_nonzero_denominators']==2 and N['no_determinant_cofactor_rank_or_GCD_algorithm']
 # A source/code audit of the primary fixed-mask recurrence: the k-th
 # row expansion sign is (-1)^(k-1+rank(column in mask)); the saved adj
 # labels transpose the deleted row/column and apply (-1)^(i+j).
 # The new arithmetic certificate independently proves the required
 # inverse formulas; it does not recompute those52 fixed determinants.
 proof=dict(variable_order=inputs['variable_order'],all70_source_clearings=checks,all50_actual_matrix_entries_bound_to_profile=True,actual_nonzero_integer_inverse_denominators=True,
  inverse_formula='x_p=Adj(Htilde)*(rho*rhs-Btilde*x_f)/D; BOTH Htilde*Adj=Adj*Htilde=D*I5 proved in ZZ[H,L,s0,q0,t0]',
  signed_cofactor_source_audit='primary fixed subset DP is the last-row Laplace recurrence, with sign k-1+column-rank; adj-i-j deletes source row j and column i',
  independent_arithmetic_scope='All100 literal inverse identities and both nonzero supplied denominator polynomials; no second determinant/cofactor/rank computation.',
  five_rows_are_actual_constraint_basis='Accepted old18 row basis modulo13 retired pure unit rows has five nonzero survivors6,14,16,41,42. Free original6,12 have identity kernel projection. Therefore the fixed five equations control all38 constraint rows overK5.',
  localization='K5=QQ(H,L,s0,q0,t0). For specialization retain old monomial denominators and both supplied D nonzero; no complementary divisor coverage/global unit/full83 rank claim.',
  no_normal_unknown_parameters=True,no_new_source_point_prime_or_GCD=True,new_polynomial_products=N['fixed_polynomial_products'],owned_native_slots=N['owned_slots_including25000_metadata_upper'])
 save('independent-literal-source-and-two-inverse-operators.json',proof)
 for p in[Path(__file__),cpp,control,exe,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 raw=sum(p.stat().st_size for p in OUT.rglob('*')if p.is_file())+sum(p.stat().st_size for p in WORK.rglob('*')if p.is_file());assert raw<32<<20 and time.monotonic()-st<60
 result=dict(id=ID,all_passed=True,frozen=True,all70_literal_source_clearings=True,all100_left_right_inverse_identities=True,both_nonzero_base_denominators=True,only_K5_and_stated_specialization_chart=True,no_determinant_rank_GCD_or_new_point=True,actual_fresh_RAW=raw,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 save('report.json',result);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
