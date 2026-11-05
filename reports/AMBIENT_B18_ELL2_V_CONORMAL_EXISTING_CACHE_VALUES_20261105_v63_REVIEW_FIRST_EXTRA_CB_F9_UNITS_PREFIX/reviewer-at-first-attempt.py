#!/usr/bin/env python3
"""Independent cache/source substitution and actual V conormal postcheck.

No implicit coefficient is solved and no producer module is imported.
The second source transport is reconstructed from the raw physical maps.
"""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,itertools,json,resource,signal,time

ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63'
STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json',
 'values':f'reports/{A}/phase1-actual105-B56-cache-composition.json',
 'sourcevalues':f'reports/{A}/phase2-all105-equations-and-primitive56-source-columns.json',
 'matrix':f'reports/{A}/actual-V-conormal26-by56-unit-minor-and-source-factorization.json',
 'profile':'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63/actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json',
 'source_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'physical':'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'catalog':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):
 (o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in o.iterdir() if p.is_file())<32*1024*1024
def put(p,m,c):
 v=p.get(m,Q(0))+c
 if v:p[m]=v
 else:p.pop(m,None)
def add(*ps):
 p={}
 for a in ps:
  for m,c in a.items():put(p,m,c)
 return p
def sc(p,c):return {m:c*v for m,v in p.items() if c*v}
def sp(p):return {int(i):Q(c) for i,c in p if Q(c)}
def enc(p):return [[i,str(c)] for i,c in sorted(p.items()) if c]
def mv(cols,v):return add(*(sc(cols[j],c) for j,c in v.items()))
class SparseJets:
 def __init__(self):
  self.mons=[()]+[(i,) for i in range(9)]
  self.mons+=list(m for m in itertools.combinations_with_replacement(range(9),2) if m!=(8,8))
  self.mons += [(7,7,7),(7,7,8)]
  self.allowed=set(self.mons);assert len(self.allowed)==56
  self.one={():Q(1)}
 def mul(self,a,b,order=3):
  p={}
  for m,c in a.items():
   for n,d in b.items():
    k=tuple(sorted(m+n))
    if len(k)<=order and k in self.allowed:put(p,k,c*d)
  return p
 def inv(self,a):
  c=a.get((),Q(0));assert c
  h=sc(add(self.one,sc(a,-1/c)),Q(1))
  v=dict(self.one);term=dict(self.one)
  for _ in range(3):term=self.mul(term,h);v=add(v,term)
  assert not self.mul(term,h)
  v=sc(v,1/c);assert self.mul(a,v)==self.one;return v
 def comp(self,coefficients,mons,args,order=2):
  out={}
  for mon,c in zip(mons,coefficients):
   if not Q(c):continue
   p={():Q(c)}
   for i in mon:p=self.mul(p,args[i],order)
   out=add(out,p)
  return out

def main():
 t=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()}
 dump(o,'phase.json',{'phase':'INDEPENDENT_CACHE_TRANSPORT_AND_ORIGINAL_SOURCE','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert D['primary']['all_passed'] and D['source_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 R=SparseJets();mons=D['FJ']['source2_monomials'];assert mons==D['physical']['source2_monomials']
 variables=[{(i,):Q(1)} for i in range(9)]
 YA=variables[:8]+[add(variables[7],sc(variables[8],Q(1,5)))]
 pa,pb=[D['physical']['families'][f] for f in ['CA','CB']]
 XA=[R.comp(p,mons,YA) for p in pa['physical_R0_readback_jets_in_own_F9']]
 Binv=[[Q(c) for c in row] for row in pb['fixed_inverse9']]
 Blin=[[Q(c) for c in row] for row in pb['R0_by_F9_unit_Jacobian']]
 assert all(sum(Binv[i][k]*Blin[k][j] for k in range(9))==Q(i==j) for i in range(9) for j in range(9))
 lin=lambda p:{m:c for m,c in p.items() if len(m)==1}
 Y1=[add(*(sc(lin(p),c) for p,c in zip(XA,row))) for row in Binv]
 XB2=[R.comp([c if len(m)==2 else 0 for m,c in zip(mons,p)],mons,Y1) for p in pb['physical_R0_readback_jets_in_own_F9']]
 YB=[add(y,*(sc(add(p,sc(lin(p),-1),sc(q,-1)),c) for p,q,c in zip(XA,XB2,row))) for y,row in zip(Y1,Binv)]
 assert [R.comp(p,mons,YB) for p in pb['physical_R0_readback_jets_in_own_F9']]==XA
 assert [R.comp(p,mons,XA) for p in pa['inverse_F9_in_common_R0_through2']]==YA
 own=D['FJ']['both57_source2']
 ca={i:R.comp(p,mons,YA) for i,p in own['CA']['all57_source2jets']}
 cb={i:R.comp(p,mons,YB) for i,p in own['CB']['all57_source2jets']}
 shared=set(ca)&set(cb);assert shared==set(range(2,11)) and all(ca[i]==cb[i] for i in shared)
 phi=ca|cb;assert len(phi)==105
 J=D['joint'];plane={}
 for i,coeff in J['all105_unknown_jets']:
  p={}
  for (a,b),c in zip(J['monomials'],coeff):
   z={():Q(c)} if Q(c) else {}
   for _ in range(a):z=R.mul(z,variables[7])
   for _ in range(b):z=R.mul(z,YA[8])
   p=add(p,z)
  plane[i]=p
  assert {m:c for m,c in phi[i].items() if all(k>=7 for k in m)}=={m:c for m,c in p.items() if len(m)<=2},('COMMON_PLANE_OVERLAP',i)
  phi[i]=add(phi[i],{m:c for m,c in p.items() if len(m)==3})
 basis=[tuple(i for i,n in enumerate(m) for _ in range(n)) for m in D['values']['basis56']]
 assert set(basis)==R.allowed;idx={m:i for i,m in enumerate(basis)}
 encode=lambda p:enc({idx[m]:c for m,c in p.items()})
 ids=J['all105_unknown_original_node_order'];assert ids==D['values']['unknown105_order']
 assert [encode(phi[i]) for i in ids]==D['values']['all105_cached_values_B56']
 # All original equations are evaluated directly at the provided cached
 # functions, including cubic coefficients from the common105 chart.
 nodes=D['source']['DAG_scalar_nodes'];roots=D['catalog']['requested_root_nodes'];roots=list(dict.fromkeys(roots+own['CA']['F9_nodes']+own['CB']['F9_nodes']))
 needed=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if i in phi:continue
  n=nodes[i];assert n[0] in ['a','m','n','i','c'],('UNBOUND',i,n)
  if n[0]!='c':todo.extend(n[1:])
 uses=collections.Counter(roots)
 for i in needed:
  if i not in phi and nodes[i][0]!='c':uses.update(nodes[i][1:])
 v={};inverses=[];peak=0
 for i in sorted(needed):
  n=nodes[i];op=n[0]
  if i in phi:p=phi[i]
  elif op=='c':p={():Q(n[1],n[2])} if n[1] else {}
  elif op=='a':p=add(v[n[1]],v[n[2]])
  elif op=='n':p=sc(v[n[1]],-1)
  elif op=='m':p=R.mul(v[n[1]],v[n[2]])
  elif op=='i':p=R.inv(v[n[1]]);inverses.append([i,n[1],encode(v[n[1]]),encode(p)])
  else:raise AssertionError(n)
  v[i]=p;peak=max(peak,sum(map(len,v.values())));assert peak<150000
  if i not in phi:
   for j in n[1:] if op!='c' else []:
    uses[j]-=1
    if not uses[j]:del v[j]
 eq=J['all105_equation_original_node_order']
 assert all(v[i]==(YA[j-48] if 48<=j<57 else {}) for j,i in enumerate(eq))
 assert all({m:c for m,c in v[i].items() if len(m)<=2}==p for i,p in zip(own['CB']['F9_nodes'],YB))
 SV=D['sourcevalues'];assert len(inverses)==6 and inverses==SV['six_unit_inverse_products']
 assert [[i,encode(v[i])] for i in D['catalog']['requested_root_nodes']]==SV['all164_original_source_roots']
 dump(o,'phase1-independent-common105-B56-source.json',{'all105_cache_values_reconstructed':True,'nonlinear_CB_transport_rebuilt_from_raw_physical_maps':True,'all105_original_equations_in_B56_zero':True,'all164_roots_and6unit_products_verified':True,'source_live_QQ':peak})
 nuA,nuB=[v[i] for i in D['catalog']['nu_roots_CA_CB']]
 V=[];fdata=[]
 for q in D['catalog']['physical12_coefficient_roots']:
  aA,aB,bA,bB=[v[q[k]] for k in ['aA_node','aB_node','bA_node','bB_node']];delta=add(aA,sc(aB,-1))
  V.extend([[bA,sc(bB,-1),{}],[delta,{},sc(bB,-1)],[{},delta,bA],[sc(R.mul(nuB,bB),-1),R.mul(nuA,bA),delta]])
  fdata.append({'row':q['physical_row'],'roots':q,'delta':encode(delta),'bA':encode(bA),'bB':encode(bB)})
 V.extend([[sc(variables[i],-1),{},{}] for i in range(7)])
 gA={}
 for (a,b),c in zip(J['monomials'],J['g_CA_all10_coefficients']):
  p={():Q(c)} if Q(c) else {}
  for _ in range(a):p=R.mul(p,variables[7])
  for _ in range(b):p=R.mul(p,YA[8])
  gA=add(gA,p)
 V.append([sc(gA,-1),{},{}]);assert len(V)==56
 assert fdata==SV['primitive_F12'] and encode(gA)==SV['actual_new_gA_plane_through3']
 assert [[encode(p) for p in col] for col in V]==SV['all56_primitive_V_columns']
 library=[]
 for d in D['profile']['scalar_monomial_adjoint_library']:
  mon=tuple(i for i,n in enumerate(d['monomial9']) for _ in range(n));library.append((d['physical_target'],mon,sp(d['conormal26_value'])))
 cols=[add(*(sc(p,col[target].get(mon,Q(0))) for target,mon,p in library)) for col in V]
 MD=D['matrix'];assert list(map(enc,cols))==MD['actual26_by56_columns']
 selected=MD['fixed_selected26_columns'];assert len(selected)==26
 fixed=[cols[j] for j in selected]
 mat=[[p.get(i,Q(0)) for p in fixed] for i in range(26)]
 assert [[str(c) for c in row] for row in mat]==MD['fixed26_minor']
 # Independent exact determinant, using only the already selected minor.
 det=Q(1)
 for j in range(26):
  p=next(i for i in range(j,26) if mat[i][j])
  if p!=j:mat[j],mat[p]=mat[p],mat[j];det=-det
  d=mat[j][j];det*=d
  for i in range(j+1,26):
   c=mat[i][j]/d
   if c:
    for k in range(j+1,26):mat[i][k]-=c*mat[j][k]
   mat[i][j]=Q(0)
 assert det and str(det)==MD['fixed26_minor_determinant']
 assert all(mv(fixed,sp(p))==col for p,col in zip(MD['all56_columns_in_selected26_basis'],cols))
 ep=o/'independent-common105-primitiveV-and26-minor-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all105_original_cache_source_functions_reconstructed':True,'all105_source_equations_and6unit_products_verified':True,'all164_roots_checked':True,'primitive56_V_columns_and_new_gA_plane_exact':True,'all26_by56_entries_recomputed':True,'selected26_columns':selected,'independent_selected26_determinant':str(det),'all56_column_factorizations_verified':True,'source_unit_section_mn_I26_accepted':True,'no_implicit_source_solve_or_new_physical_coefficients':True,'actual_T2_or_H0_length_not_computed':True})
 report={'id':STEM,'all_passed':True,'actual_V_conormal_rank26_accepted':True,'formal_unit_section_n_and_h0_source_gate_accepted':True,'actual_T2_or_H0_length_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep),'phase1-independent-common105-B56-source.json':sha(o/'phase1-independent-common105-B56-source.json')},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
