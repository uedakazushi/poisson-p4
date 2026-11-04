#!/usr/bin/env python3
"""SOURCE degree/support majorants only; no parameter products or minors evaluated."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_H22_KERNEL9_SUPPORT_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
Z=None
def vol(b):return 0 if b is None else __import__('math').prod(hi-lo+1 for lo,hi in b[1])
def atom(p):
 if not p:return Z
 return (len(p),tuple((min(m[j]for m,c in p),max(m[j]for m,c in p))for j in range(3)),0)
BP=([[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']])
def aa(a,b):
 if a is None:return b
 if b is None:return a
 power=max(a[2],b[2]);a=numB(a,power-a[2]);b=numB(b,power-b[2]);box=tuple((min(a[1][j][0],b[1][j][0]),max(a[1][j][1],b[1][j][1]))for j in range(3));r=(a[0]+b[0],box);return(min(r[0],vol(r)),box,power)
def mm(a,b):
 if a is None or b is None:return Z
 box=tuple((a[1][j][0]+b[1][j][0],a[1][j][1]+b[1][j][1])for j in range(3));r=(a[0]*b[0],box);return(min(r[0],vol(r)),box,a[2]+b[2])
def numB(a,k):
 if a is None:return Z
 power=a[2];b=(a[0],a[1],0)
 for _ in range(k):b=mm(b,atom(BP))
 return(b[0],b[1],power)
ONE=atom([[[0,0,0],'1']])
def shift(a,k):
 if a is None:return Z
 return(a[0],((a[1][0][0]+k,a[1][0][1]+k),a[1][1],a[1][2]),a[2])
def enc(a):return dict(zero=True,support_upper=0)if a is None else dict(zero=False,support_upper=a[0],Laurent_exponent_box=a[1],clear_a_power=max(0,-a[1][0][0]),known_B_denominator_power=a[2])
def size(a):return 0 if a is None else a[0]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists();phase=read(out/'phase.json');paths=[]
 def load(stem,file='certificate.json'):
  p=ROOT/'reports'/stem/file;paths.append(p);return read(p)
 cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63')
 hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');rank=load('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63')
 tr=load('S61_NONSQUARE_PREPOWER_HEAVY_Q_CHART_20261104_v63','report.json');assert sha(paths[-1])=='45022e18b75a98f93f4e8076a3432944a5f1ddc13ef5d0d7fbd9522ef5e36ca7'
 nodes=q['scalar_DAG']['nodes'];v=[];raw_atom_terms=0;literal_denominators=[]
 def rf(nd,ref):
  nonlocal raw_atom_terms
  n,d=nd;raw_atom_terms+=len(n)+len(d)
  if len(d)==1:
   assert d[0][0][1:]==[0,0]and Q(d[0][1]);literal_denominators.append(dict(source=ref,a_power=d[0][0][0],B_power=0,QQ_unit=d[0][1]));return shift(atom(n),-d[0][0][0])
  dd={tuple(m):Q(c)for m,c in d};bb={tuple(m):Q(c)for m,c in BP};assert dd.keys()==bb.keys();ratio=dd[next(iter(bb))]/bb[next(iter(bb))];assert ratio and all(dd[m]==ratio*bb[m]for m in bb)
  literal_denominators.append(dict(source=ref,a_power=0,B_power=1,QQ_unit=str(ratio)));b=atom(n);return Z if b is None else(b[0],b[1],1)
 for i,n in enumerate(nodes):
  if n[0]=='z':b=Z
  elif n[0]=='one':b=ONE
  elif n[0]=='checkpoint':b=rf(cp['matrix25x32'][n[1]][n[2]],['checkpoint',n[1],n[2]])
  elif n[0]=='accepted_pair_second':b=rf(pd['Lframe_dual6x2'][n[1]][1],['accepted_pair_second',n[1]])
  elif n[0]=='new_unit_atom':b=rf(q['scalar_DAG']['new_unit_atoms_num_den'][n[1]],['new_unit_atom',n[1]])
  elif n[0]=='m':b=mm(v[n[1]],v[n[2]])
  elif n[0]in['a','s']:b=Z if n[0]=='s'and n[1]==n[2]else aa(v[n[1]],v[n[2]])
  else:raise AssertionError(n)
  v.append(b)
 cols=q['new22_target_original32_columns'][:17];matrix=[[v[i]for i in row[:17]]for row in q['remaining15x22_scalar_roots']]
 for k,label in enumerate(h['original_target_labels']):matrix[0][cols.index(label)]=atom(h['source_h9'][k])
 pivot=cols.index(22);h22=matrix[0][pivot];assert len(h['source_h9'][7])==14
 pure3ids=[1,3,6];pure3=[]
 for i in pure3ids:pure3.append([aa(mm(h22,matrix[i][j]),mm(matrix[i][pivot],matrix[0][j]))for j in range(17)if j!=pivot])
 targets=[j for j in cols if j!=22]
 candidates=[]
 for C in itertools.combinations(range(16),3):
  b=Z;terms=0
  for perm in itertools.permutations(range(3)):
   t=ONE
   for i in range(3):t=mm(t,pure3[i][C[perm[i]]])
   b=aa(b,t)
   if t is not None:terms+=1
  if terms:
   candidates.append(dict(targets=[targets[j]for j in C],support_upper=size(b),term_count_nonzero_support_matchings=terms,
    maximum_input_entry_support=max(size(pure3[i][j])for i in range(3)for j in C),sum_input_entry_support=sum(size(pure3[i][j])for i in range(3)for j in C),bound=enc(b)))
 candidates.sort(key=lambda x:(x['support_upper'],x['sum_input_entry_support'],x['targets']))
 # Block-triangular2+3 decomposition, using ONLY already materialized cofactors.
 N=[[atom(p)for p in row]for row in hc['raw_compact14x5']];pair=[1,3];paircols=[0,4];brows=[6,8,9,10,11,12,13];Bsel=[6,12,9];Brest=[8,10,11,13]
 A2=[[N[i][j]for j in paircols]for i in pair];adj2=[[A2[1][1],A2[0][1]],[A2[1][0],A2[0][0]]]
 det2=atom(rank['literal_pair_determinant']);det3=atom(rank['nonzero_selected_minor_Delta']);lc=rank['cleared_upper_row_coordinates7xr']
 # rho common-denominator factors: exact known-unit monomial LCM, no gcd.
 named=[atom([[[1,0,0],'1']]),atom([[[0,0,1],'1']]),atom([[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']]),atom([[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]),atom([[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']]),atom([[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']])]
 def unitbound(powers):
  b=ONE
  for x,n in zip(named,powers):
   for _ in range(n):b=mm(b,x)
  return b
 coupled=[]
 for u in [4,7]+Brest:
  source={}
  if u in[4,7]:source[u]=ONE;c=[N[u][j]for j in paircols];pureheavy_scalar=det2
  else:
   lam=[atom(p)for p in lc[brows.index(u)]];source[u]=det3
   for j,i in enumerate(Bsel):source[i]=lam[j]
   c=[aa(mm(det3,N[u][j]),aa(mm(lam[0],N[Bsel[0]][j]),aa(mm(lam[1],N[Bsel[1]][j]),mm(lam[2],N[Bsel[2]][j]))))for j in paircols];pureheavy_scalar=det2
  source={i:mm(pureheavy_scalar,b)for i,b in source.items()}
  for j,i in enumerate(pair):source[i]=aa(mm(c[0],adj2[0][j]),mm(c[1],adj2[1][j]))
  supportrows=[i for i,b in source.items()if b is not None];powers=[[x for x in hc['known_row_units_rho14'][i]['a_d_F_B_C_J_powers']]for i in supportrows];lcm=[max(row[j]for row in powers)for j in range(6)]
  src={i:mm(b,unitbound([lcm[j]-hc['known_row_units_rho14'][i]['a_d_F_B_C_J_powers'][j]for j in range(6)]))for i,b in source.items()}
  kp=[]
  for j in range(17):
   b=Z
   for i,a in src.items():b=aa(b,mm(a,matrix[i+1][j]))
   kp.append(b)
  image=[aa(mm(h22,kp[j]),mm(kp[pivot],matrix[0][j]))for j in range(17)if j!=pivot]
  coupled.append(dict(compact_source_row=u,kind='pair_only'if u in[4,7]else'B3_kernel_then_pair',
   source_cofactor_support_bounds={str(i):enc(b)for i,b in src.items()},row_known_rho_LCM_powers=lcm,
   pure16_support_bounds=[enc(b)for b in image],row_total_support_upper=sum(map(size,image)),max_entry_support_upper=max(map(size,image))))
 proof=dict(source_only=True,polynomial_coefficients_or_parameter_products_NOT_evaluated=True,
  support_model='Laurent a exponent boxes with min(convolution-support bound,box volume); additive cancellation NOT assumed except literalzero/same-root subtraction and accepted nine actualh aliases.',
  all_source_DAG_atoms_have_only_a_B_and_QQ_denominators=True,
  first_a_only_denominator_assumption_failed_preserved='first-attempt-frontier.json; 99needed checkpoint atoms have literalQQ*Bdenominators, nowtracked explicitly withoutgcd.',
  actual_h22_terms14=True,all8_Q9_tails_retained=True,
  pure3_source_labels=[10,13,16],pure3_actual135=[18,39,68],pure3_Hadjugate_NOT_needed=True,
  staged_Hchart=dict(pair_rows=pair,pair_columns=paircols,det2_formula='-(6/5)*B*p^2',B3_rows=Bsel,B3_columns=[1,3,2],
   det3_terms=len(rank['nonzero_selected_minor_Delta']),det5_terms=len(rank['nonzero_minor_of14x5']),saved_B7_Cramer_cofactor_terms=sum(len(p)for row in lc for p in row),
   exact_source_formula='For four B-kernel rows, first det3*G_u-sumL_ui G_Bi, then det2*that - (remainingHpair)*adj(A2)*G_pair. The two pair-only rows use the second operation directly. Three Hzero rows use neither. This is the same D(det2*det3) chart; allrho factors are knownFULLU units.'),
  selected_chart='FULLU intersect D(existingDelta5 * actualh22); extra complement retained, includingDelta0,p0 andotherh charts.',
  no_new_3minor_or_rankclaim=True,cheapest3minor_is_only_source_support_candidate=True,
  whole_FULLU57_or_quartic239_NOT_proved=True)
 data=dict(proof=proof,pure17_target_labels=cols,pure16_target_labels=targets,Q9_h22_bound=enc(h22),
  pure3x16_bounds=[[enc(b)for b in row]for row in pure3],pure3_row_total_support_upper=[sum(map(size,row))for row in pure3],
  cheapest12_pure3_fixed_minor_candidates=candidates[:12],six_staged_coupled_row_bounds=coupled,
  all9x16_support_upper=sum(size(b)for row in pure3 for b in row)+sum(r['row_total_support_upper']for r in coupled),
  retained_source_atom_terms_counting_references=raw_atom_terms,
  denominator_inventory=literal_denominators)
 (out/'certificate.json').write_text(json.dumps(data,indent=2)+'\n');phase.update(phase='source_support_preflight_complete_no_value_job',all9_support_upper=data['all9x16_support_upper']);(out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
 paths+=[Path(__file__),ROOT/f'research/notes/{ID}.md'];report=dict(id=ID,all_passed=True,source_only=True,frozen=True,
  pure3_row_support_upper=data['pure3_row_total_support_upper'],six_coupled_row_support_upper=[r['row_total_support_upper']for r in coupled],
  all9x16_support_upper=data['all9x16_support_upper'],cheapest3minor_support_upper=candidates[0]['support_upper'],cheapest3minor_targets=candidates[0]['targets'],
  no_value_minor_rank_point_prime_GB=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},certificate_sha256=sha(out/'certificate.json'),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
