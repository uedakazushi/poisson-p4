#!/usr/bin/env python3
"""Independent positive cached right85 initial-source and rank proof."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def plus(v,w,c=F(1)):
 z=dict(v)
 for i,a in w.items():
  z[i]=z.get(i,F(0))+c*a
  if not z[i]:del z[i]
 return z
def prod(A,B):
 cols=[{i:c for i,c in enumerate(col)if c}for col in zip(*B)];rows=[{i:c for i,c in enumerate(row)if c}for row in A]
 return [[sum((c*col.get(i,F(0))for i,c in row.items()),F(0))for col in cols]for row in rows]
def identity(n):return [[F(i==j)for j in range(n)]for i in range(n)]
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63';r=load(primary,h='eb28ca1caa331bdd0a82d86decd157dd7dce73e53acb2fd749817093ed55e183');assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 e=load(primary,'all85-actual-cached-cubic-values-and-source-membership.json',h='2ab3b00768dbf01cc640718cbbadad511584b0db769ce70eb3e4aabd02223f74')
 assert load('AMBIENT_B18_CPRIME_RIGHT85_SOURCE_20261105_v63_REVIEW',h='8caeb044c5e71ba6a88e633ca2d80f4120fc9f977b304761c5ad37022ca6f9ba')['complete_right85_kernel_source_and_adjoint_gate_passes']
 source=load('AMBIENT_B18_CPRIME_RIGHT85_SOURCE_PREFLIGHT_20261105_v63','complete-right85-and-cached-plane-source-adjoints.json');qv=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63','phase4-all-twenty-actual-quadratic-values.json');assert load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW',h='5faf10001ed74c2b9e107d7fa1f4fe873600423affa75305a75a8f1d3fad5bb3')['actual_Cprime_HF_through2']==[3,3,3]
 assert load('AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW',h='481862051b9e604ac92163b452f88e0cb46682d043393668b6648a89610e3559')['known_cubic_rank']==26
 old=load('AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63','phase2-complete-streamed-cubic-upper-certificate.json');dualdata=load('AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW','independent-known-cubic-fixed-certificate-and-open-dual.json');dual={int(i):F(c)for i,c in dualdata['known_quotient_dual495'].items()};q27={int(i):F(c)for i,c in dualdata['known_quotient_dual27'].items()}
 cp=source['actual_V1_fixed24_source_columns'];raw=[list(map(F,col))for col in qv['all44_actual_V2_pureplane_nine_values']];K=[pairs(k)for k in source['all85_syzygy_columns216_sparse']];adj=[pairs(a)for a in source['all85_cached_plane_Q2_adjoints396_sparse']];m3=list(itertools.combinations_with_replacement(range(9),3));i3={m:i for i,m in enumerate(m3)};polys=[];values=[]
 for k,a in zip(K,adj):
  p={}
  for j,c in k.items():
   sourcej,x=divmod(j,9)
   if x<7:continue
   for target in range(3):
    for pair,mon in enumerate([(7,7),(7,8),(8,8)]):p=plus(p,{165*target+i3[tuple(sorted(mon+(x,)))]:c*raw[cp[sourcej]][3*target+pair]})
  value=sum((c*p.get(i,F(0))for i,c in dual.items()),F(0));folded=sum((c*raw[j//9][j%9]for j,c in a.items()),F(0));assert folded==value;polys.append(p);values.append(value)
 assert len(values)==85 and values==list(map(F,e['all85_actual_residual_cubic_values']))and polys==[pairs(p)for p in e['all85_actual_plane_cubic_initial_representatives_sparse']]
 nz=[i for i,v in enumerate(values)if v];assert nz==e['nonzero_indices']==[49,57,59,77,84]and values[49]==F(99,64)
 # Explicit actual source witness, with the degree2 zero checked in ALL135
 # coordinates. No higher source basis changes or M3 inputs occur.
 first=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase1-actual-fullnine-first-V-columns.json');L=[pairs(v)for v in first['all44_actual_linear_columns27_sparse']];m2=list(itertools.combinations_with_replacement(range(9),2));i2={m:i for i,m in enumerate(m2)};zero={};lift=[]
 for j,c in K[49].items():
  sj,x=divmod(j,9);original=cp[sj];lift.append([original,x,c])
  for targeti,a in L[original].items():t,y=divmod(targeti,9);zero=plus(zero,{45*t+i2[tuple(sorted((x,y)))]:c*a})
 assert not zero
 proj=[pairs(v)for v in old['all495_source_projection_columns27_sparse']];new={}
 for i,c in polys[49].items():new=plus(new,proj[i],c)
 assert sum((c*new.get(i,F(0))for i,c in q27.items()),F(0))==F(99,64)
 rr=old['fixed_rows'];basis=[pairs(v)for v in old['selected_original_compressed_columns']];A=[[p.get(i,F(0))for p in basis]for i in rr];Ai=[[F(c)for c in row]for row in old['fixed_inverse']];assert prod(A,Ai)==prod(Ai,A)==identity(26)
 top=[new.get(i,F(0))for i in rr];bottom=[p.get(25,F(0))for p in basis];d=new.get(25,F(0));u=[sum((c*v for c,v in zip(row,top)),F(0))for row in Ai];v=[sum((bottom[i]*Ai[i][j]for i in range(26)),F(0))for j in range(26)];s=d-sum((c*x for c,x in zip(bottom,u)),F(0));assert s==F(99,64)
 inverse=[[Ai[i][j]+u[i]*v[j]/s for j in range(26)]+[-u[i]/s]for i in range(26)]+[[-c/s for c in v]+[1/s]];matrix=[A[i]+[top[i]]for i in range(26)]+[bottom+[d]];assert prod(matrix,inverse)==prod(inverse,matrix)==identity(27)
 det=F(old['fixed_minor_determinant'])*s;assert det==F(-671,128000)
 assert e['source_kernel_sha256']==source['kernel_sha256']and r['actual_right85_rank']==1 and r['positive_Cprime_ell3_length9_and_long_exponent2']
 rank=load('AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63');assert rank['actual_generic_DVR_rank_Cprime']==2 and rank['actual_DVR_torsion_minimal_generators']==2
 assert load('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW',h='92d5c943b12da089eb504862043a6e811ad17bef9d07e2b693e4eaf07a539b30')['actual_Cprime2_length']==7
 save('independent-right85-positive-source-and-full27-certificate.json',dict(all85_actual_residual_values=values,all85_pureplane_projections_of_actual_initial_forms=polys,all85_direct_and_folded_pairings_match=True,nonzero_indices=nz,positive_witness_index=49,positive_original_V_source_linear_coefficients=lift,positive_source_V1_all135_coefficients_zero=True,positive_actual_cubic_projection27=new,positive_residual_scalar=s,new27_rows=rr+[25],new27_source_columns=old['fixed_columns']+['actualV_u1_witness49'],new27_matrix=matrix,new27_inverse=inverse,both27_inverse_products_verified=True,new27_determinant=det,full_actual_cubic_rank27_by_known26_plus_genuine_initial=True,positive_completion='Actualinitialrelations spanall degree3 in Cprime/ell3. Nakayama givesm3=0; actualHF333 giveslength9. DVRrank2/torsmu2/Cprime2length7 imply{1,e>=2}; lengthell3=7+min(e,3)=9 forces e2.',actual_Cprime_ell3_length=9,actual_Cprime_DVR_torsion_exponents=[1,2],no_original_nilideal_index_or_generators_claim=True,no_genuine_connecting_target_basis_or_actions_built=True,no_new_source_gauge_jet_or_rank_search=True))
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 total=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir()if p.is_file());assert total<2_000_000 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,actual_right85_rank=1,all85_source_cached_values_verified=True,explicit_actual27_initial_minor_and_both_inverse_products=True,actual_Cprime_ell3_length=9,actual_Cprime_DVR_torsion_exponents=[1,2],positive_long_exponent2_proved=True,no_new_source_gauge_jet_or_rank_search=True,no_original_nilindex_or_connecting_basis_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=total,seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
 print(json.dumps(dict(all_passed=True,first_nonzero='99/64',actual_Cprime_ell3_length=9,torsion_exponents=[1,2],seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
