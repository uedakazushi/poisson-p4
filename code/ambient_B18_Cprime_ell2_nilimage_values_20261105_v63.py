#!/usr/bin/env python3
"""Finite33 canonical unit gauge and right-exact Cprime/ell² only.
Reuses actual105 physical functions; does not solve new physical equations.
"""
import collections,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_Cprime_ell2_nilimage_source_preflight_20261105_v63 import PATHS
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63'
GATE='AMBIENT_B18_CPRIME_ELL2_NILIMAGE_SOURCE_PREFLIGHT_20261105_v63'
EXTRA={
 'sourcepeer':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_SOURCE_20261105_v63_REVIEW/report.json',
 'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
 'B11gauge':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase2-gauge30-B11-correction.json',
 'B11proper':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'oldCbar':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-physical12-actions-and-C-special.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def sp(v):return [[i,c]for i,c in sorted(v.items())if c]
def put(v,i,c):
 if c:
  v[i]=v.get(i,Q(0))+c
  if not v[i]:del v[i]
def add(v,w,c=Q(1)):
 z=dict(v)
 for i,a in w.items():put(z,i,c*a)
 return z
def mv(A,v):
 z={}
 for j,c in v.items():z=add(z,A[j],c)
 return z
def bits(v):
 assert all(max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=2048 for c in v.values()),'DERIVED_2048_BIT_GUARD'
 return v
class Base:
 def __init__(self,d):self.n=33;self.table=[[rec(v)for v in r]for r in d['table_columns_for_each_left_basis']];self.deg=d['basis_degrees']
 def mul(self,v,w):
  z={}
  for i,c in v.items():
   for j,a in w.items():
    if self.deg[i]+self.deg[j]<5:z=add(z,self.table[i][j],c*a)
  return bits(z)
 def inv(self,v):
  c=v.get(0,Q(0));assert c
  h={i:-a/c for i,a in v.items()if i};z={0:Q(1)};p=z
  for _ in range(4):p=self.mul(p,h);z=add(z,p)
  z={i:a/c for i,a in z.items()};assert self.mul(v,z)=={0:Q(1)};return bits(z)
 def order(self,v):return min((self.deg[i]for i,c in v.items()if c),default=5)
class Trace:
 def __init__(self,B,m=0,nus=()):
  self.B=B;self.m=m;self.b=1<<m;self.n=33*self.b;self.nu=[rec(v)if isinstance(v,list)and v and isinstance(v[0],list)else {i:Q(c)for i,c in enumerate(v)if Q(c)}for v in nus];self.Z=(Q(0),)*self.n;self.ONE=self.c(1)
 def c(self,c):return (Q(c),)+(Q(0),)*(self.n-1)
 def add(self,*ps):
  z=tuple(sum((p[i]for p in ps),Q(0))for i in range(self.n));bits({i:c for i,c in enumerate(z)if c});return z
 def neg(self,p):return tuple(-c for c in p)
 def scale(self,p,c):
  z=tuple(c*a for a in p);bits({i:a for i,a in enumerate(z)if a});return z
 def mul(self,*ps):
  v=self.ONE
  for p in ps:v=self.mul2(v,p)
  return v
 def mul2(self,p,q):
  if p==self.Z or q==self.Z:return self.Z
  if p==self.ONE:return q
  if q==self.ONE:return p
  pp=[{i:c for i,c in enumerate(p[k::self.b])if c}for k in range(self.b)];qq=[{i:c for i,c in enumerate(q[k::self.b])if c}for k in range(self.b)];z=[{}for _ in range(self.b)]
  for i,x in enumerate(pp):
   for j,y in enumerate(qq):
    if not x or not y:continue
    v=self.B.mul(x,y);carry=i&j
    for k in range(self.m):
     if carry&(1<<k):v=self.B.mul(v,self.nu[k])
    z[i^j]=add(z[i^j],v)
  for v in z:bits(v)
  return tuple(z[k].get(i,Q(0))for i in range(33)for k in range(self.b))
 def inv(self,p):
  assert p[0];h=self.add(self.ONE,self.scale(p,-1/p[0]));v=self.ONE;z=self.ONE
  for _ in range(10):v=self.mul(v,h);z=self.add(z,v)
  z=self.scale(z,1/p[0]);assert self.mul(p,z)==self.ONE;return z
 def power(self,p,n):
  if n<0:return self.power(self.inv(p),-n)
  z=self.ONE
  while n:
   if n&1:z=self.mul(z,p)
   n//=2
   if n:p=self.mul(p,p)
  return z
 def dual(self,p):return tuple(Q(p[i])if k==0 else Q(0)for i in range(33)for k in range(self.b))
 def tau(self,k):
  v=list(self.Z);v[1<<k]=Q(1);return tuple(v)
 def install(self):
  for k in ['c','add','neg','scale','mul','inv','power']:setattr(H,k,getattr(self,k))
  H.Z=self.Z;H.ONE=self.ONE;H.TA=self.tau(0)if self.m else self.Z;H.TB=self.tau(1)if self.m>1 else self.Z
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 def save(name,d):
  blob=json.dumps(enc(d),separators=(',',':'))+'\n';assert len(blob.encode())+sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32<<20;(out/name).write_text(blob);(out/'phase.json').write_text(json.dumps({'phase':name,'seconds':time.monotonic()-t})+'\n')
 gp=ROOT/'reports'/GATE/'report.json';gr=json.loads(gp.read_text());assert gr['all_passed']and gr['whole_saved_live_QQ_upper']==118256
 for p,h in gr['input_sha256'].items():assert sha(ROOT/p)==h,p
 paths={**PATHS,**EXTRA};pins={p:sha(ROOT/p)for p in paths.values()};pins[str(gp.relative_to(ROOT))]=sha(gp)
 save('preflight.json',{'input_sha256':pins,'hard_saved_live_QQ':118256,'hard_derived_bits':2048,'new_physical_source_values':0,'mandatory_extra_overlap_inputs':EXTRA,'no_ell_division':True})
 D={k:json.loads((ROOT/p).read_text())for k,p in paths.items()};assert D['sourcepeer']['all_passed'];B=Base(D['table33']);b=Trace(B);S=D['source'];T=D['source33'];first={i:tuple(rec(v).get(k,Q(0))for k in range(33))for i,v in zip(T['joint105_unknown_order'],T['all105_unknown_final_values'])};nodes=S['DAG_scalar_nodes']
 roots=[i for f in S['families'].values()for p in f['complete_normalized_physical_w_coefficients'][0]for i in p]+[f['actual_tracefree_nu']for f in S['families'].values()];seen=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in first:continue
  n=nodes[i]
  if n[0]in ['a','m','n','i']:todo.extend(n[1:])
  else:assert n[0]=='c'
 uses=collections.Counter(roots)
 for i in seen:
  if i not in first and nodes[i][0]in ['a','m','n','i']:uses.update(nodes[i][1:])
 vals={};live=peak=0
 for i in sorted(seen):
  n=nodes[i];o=n[0]
  if i in first:v=first[i]
  elif o=='c':v=b.c(Q(n[1],n[2]))
  elif o=='a':v=b.add(vals[n[1]],vals[n[2]])
  elif o=='m':v=b.mul(vals[n[1]],vals[n[2]])
  elif o=='n':v=b.neg(vals[n[1]])
  else:assert o=='i';v=b.inv(vals[n[1]])
  vals[i]=v;live+=33;peak=max(peak,live);assert peak<=13299
  if i not in first and o in ['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j];live-=33
 overlap=0
 for i,v in T['all_requested164_root_values']:
  if i in vals:assert vals[i]==tuple(rec(v).get(k,Q(0))for k in range(33));overlap+=1
 pairs={f:[[vals[x],vals[y]]for x,y in data['complete_normalized_physical_w_coefficients'][0]]for f,data in S['families'].items()};nu=[vals[f['actual_tracefree_nu']]for f in S['families'].values()];assert all(not v[0]for v in nu)
 save('phase1-source150-readback-fingerprints.json',{'all600_physical_root_readbacks_sha256':{f:hashlib.sha256(json.dumps(enc(v),separators=(',',':')).encode()).hexdigest()for f,v in pairs.items()},'source164_exact_overlap_count':overlap,'source_stream_peak_QQ':peak,'nu_A_B33':nu,'original105_values_not_recomputed':True});del vals
 I,P,GL=D['interface'],D['proper'],D['firstpreflight']['GL18_source_parameter_positions'];pos=I['actual_group15_positions'];rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices'];centers=list(map(Q,P['actual_CA_source19_center']))
 def frame_family(a,family,mask):
  a.install();f=S['families'][family];p=H.add(a.tau(mask),H.scale(a.dual(first[f['p_quadratic_relation']['s']]),Q(1,2)));g=H.ident(5)
  for name,(x,y)in f['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);g[i][j]=H.add(g[i][j],a.dual(first[x]),H.mul(a.dual(first[y]),p))
  core=H.ident(5)
  if family=='CA':core[0][3]=H.scale(p,-2);core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
  else:core[0][3]=H.scale(p,Q(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),Q(1,4))
  return H.mm(g,core)
 def physical(a,family,mask):
  a.install();return [H.add(a.dual(x),H.mul(a.dual(y),a.tau(mask)))for x,y in pairs[family]]
 ev=[tuple([Q(v[0])]+[Q(0)]*32)for v in D['gauge1']['even15_first9']];od=[tuple([Q(v[0])]+[Q(0)]*32)for v in D['gauge1']['odd15_first9']]
 def gauge(a):
  a.install();g=H.ident(5)
  for(i,j),x,y in zip(pos,ev,od):g[i][j]=H.add(g[i][j],a.dual(x),H.mul(a.dual(y),a.tau(0)))
  return g
 aa=Trace(B,1,[nu[0]]);HA=frame_family(aa,'CA',0);PA=physical(aa,'CA',0)
 def constraints(g):
  aa.install();fa=H.frame(H.mm(g,HA));tr=H.mm(fa,g);moved=H.physical_push(PA,tr,H.mi(tr));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in ar]);own=[v[0]for v in own]
  return [fa[1][0],fa[1][3],fa[1][4],fa[0][3],fa[0][4]]+[H.add(own[j],H.c(-centers[j]))for j in fixed]
 Ji=[[(j,Q(c))for j,c in enumerate(row)if Q(c)]for row in D['J30']['J30_inverse']]
 for stage in range(1,5):
  rr=constraints(gauge(aa));res=[{k:c for k,c in enumerate(v[parity::2])if c}for parity in range(2)for v in rr];assert min(B.order(v)for v in res)>=stage
  corr=[]
  for row in Ji:
   v={}
   for j,c in row:v=add(v,res[j],-c)
   corr.append(bits(v))
  for i in range(15):ev[i]=tuple(Q(ev[i][k])+corr[i].get(k,Q(0))for k in range(33));od[i]=tuple(Q(od[i][k])+corr[15+i].get(k,Q(0))for k in range(33))
  save(f'stage{stage}-canonical30-residual-correction.json',{'minimum_residual_order':min(B.order(v)for v in res),'residual30':list(map(sp,res)),'correction30':list(map(sp,corr))})
 assert all(v==aa.Z for v in constraints(gauge(aa)))
 # Exact source-frame overlap to the independently accepted B11 gauge.
 mons=[tuple(m)for m in D['B11']['B11_basis_monomials8']];mi={m:i for i,m in enumerate(mons)};base=D['table33']['basis33'];maps=[]
 for z,r,e in base:
  mon=[0]*8;mon[7]=r
  if z:mon[z-1]+=1
  maps.append({mi[tuple(mon)]:Q(1)}if not e and tuple(mon)in mi else{})
 def proj11(v):return mv(maps,{i:c for i,c in enumerate(v)if c})
 assert all(proj11(v)=={i:Q(c)for i,c in enumerate(w)if Q(c)}for v,w in zip(ev,D['B11gauge']['even15_B11']))
 assert all(proj11(v)=={i:Q(c)for i,c in enumerate(w)if Q(c)}for v,w in zip(od,D['B11gauge']['odd15_B11']))
 save('phase2-canonical30-exact-values.json',{'actual_even15_B33':ev,'actual_odd15_B33':od,'all30_constraints_exact':True,'all30_gauge_B11_projections_exact':True,'fixed_J30_determinant':'4/9'});del HA,PA
 a=Trace(B,2,nu);gg=gauge(a);HB=frame_family(a,'CB',1);FB=H.frame(H.mm(gg,HB));del HB
 HA=frame_family(a,'CA',0);assert H.frame(H.mm(gg,HA))==H.ident(5);del HA
 PA=physical(a,'CA',0);moved=H.physical_push(PA,gg,H.mi(gg));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in ar]);own=[v[0]for v in own];assert all(own[j]==H.c(centers[j])for j in fixed);del PA,moved
 actual=[H.add(own[j],H.c(-centers[j]))for j in free]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j]for j in I['retained_source9_indices']]
 for v,w in zip(kappa,D['B11proper']['kappa9_C_B44']):assert [proj11(v[mask::4]).get(i,Q(0))for i in range(11)for mask in range(4)]==list(map(Q,w))
 core=[kappa[i]for i in[0,1,2,6]];av,d,x,tt=core;hh=H.add(av,H.neg(d),x,H.scale(tt,-2));cc=H.scale(H.add(d,x),Q(1,2));ll=H.scale(H.add(av,d,H.neg(x)),Q(1,3));nil=[H.mul(d,hh),H.mul(d,cc,ll)];nil8=[H.mul(v,a.tau(0)if mask==1 else a.tau(1)if mask==2 else H.mul(a.tau(0),a.tau(1))if mask==3 else a.ONE)for v in nil for mask in range(4)]
 save('phase3-actual-kappa9-and-normalized-nilframe.json',{'actual14_C132':actual,'kappa9_C132':kappa,'actual_nil2_C132':nil,'actual_nil8_C132':nil8,'all_kappa_B11_projections_exact':True,'actual_same_CA_frame_and_fiber_constraints_exact':True})
 C2=D['C2'];raw=D['raw33'];q6=[[rec(v)for v in b['quotient_map6_by33_columns']]for b in raw['canonical_branch_quotients6']];proj=[rec(v)for v in C2['all45_projection_columns']]
 def raw45(v):
  a0=mv(q6[0],{i:c for i,c in enumerate(v[1::4])if c});b0=mv(q6[1],{i:c for i,c in enumerate(v[2::4])if c});h0={i:c for i,c in enumerate(v[3::4])if c};return dict(list(a0.items())+[(i+6,c)for i,c in b0.items()]+[(i+12,c)for i,c in h0.items()])
 images=[mv(proj,raw45(v))for v in nil8];acts=[[rec(v)for v in A]for A in C2['C2_coordinate_actions9_sparse_columns']]
 def basact(v,mon):
  z,r,e=mon
  if z:v=mv(acts[z-1],v)
  for _ in range(r):v=mv(acts[7],v)
  for _ in range(e):v=mv(acts[8],v)
  return v
 closure=[basact(v,m)for v in images for m in base];assert len(closure)==264
 piv={};sel=[]
 def reduce(v):
  for p,row in sorted(piv.items()):v=add(v,row,-v.get(p,Q(0)))
  return v
 for j,v in enumerate(closure):
  z=reduce(v)
  if z:
   p=min(z);c=z[p];piv[p]={i:a/c for i,a in z.items()};sel.append(j)
 assert all(not reduce(v)for v in closure)
 assert all(not reduce(mv(A,v))for A in acts for v in piv.values())
 free=[i for i in range(14)if i not in piv];q=len(free);assert 6<=q<=8
 def qp(v):return {i:reduce(v).get(k,Q(0))for i,k in enumerate(free)if reduce(v).get(k,Q(0))}
 qacts=[[qp(mv(A,{j:Q(1)}))for j in free]for A in acts];assert all(mv(A,mv(B,{j:Q(1)}))==mv(B,mv(A,{j:Q(1)}))for A in qacts for B in qacts for j in range(q))
 # Mandatory exact ell0 nil-frame and quotient/action readback.
 old26=[rec(v)for v in D['oldCbar']['all26_projection_columns']];to26=[rec(v)for v in C2['raw45_to_raw26_ell0_columns']];oldnil=D['Cprime']['all8_actual_nil_image_vectors_in_Cbar8']
 for v,w in zip(nil8,oldnil):assert mv(old26,mv(to26,raw45(v)))=={i:Q(c)for i,c in enumerate(w)if Q(c)}
 cbpiv={p:{i:Q(c)for i,c in enumerate(v)if Q(c)}for p,v in D['Cprime']['normalized_triangular_image_basis']};cbfree=D['Cprime']['quotient_original_basis_indices']
 def cbred(v):
  for p,w in sorted(cbpiv.items()):v=add(v,w,-v.get(p,Q(0)))
  return {i:v[k]for i,k in enumerate(cbfree)if v.get(k,Q(0))}
 C2rawbasis=C2['C2_basis_raw45_indices'];barcols=[cbred(mv(old26,to26[j]))for j in C2rawbasis];qbar=[barcols[j]for j in free];cbacts=[[{i:Q(c)for i,c in enumerate(v)if Q(c)}for v in A]for A in D['Cprime']['all8_quotient_coordinate_actions_columns']]+[[{}for _ in range(4)]]
 assert all(mv(qbar,mv(A,{j:Q(1)}))==mv(B,mv(qbar,{j:Q(1)}))for A,B in zip(qacts,cbacts)for j in range(q))
 def rank(cols):
  ps={}
  for v in cols:
   z=dict(v)
   for p,w in sorted(ps.items()):z=add(z,w,-z.get(p,Q(0)))
   if z:
    p=min(z);c=z[p];ps[p]={i:a/c for i,a in z.items()}
  return len(ps)
 assert rank(qbar)==4 and rank(qacts[8])==q-4 and all(not mv(qacts[8],v)for v in qacts[8])
 data={'all8_nil_images_C2':list(map(sp,images)),'all264_R0_coefficient_closure_columns':list(map(sp,closure)),'fixed_image_pivot_columns':sel,'fixed_image_normalized_basis':[[p,sp(v)]for p,v in sorted(piv.items())],'image_rank':len(piv),'actual_Cprime2_length':q,'quotient_original_C2_basis_indices':free,'actual_Cprime2_coordinate_actions9':[list(map(sp,A))for A in qacts],'actual_Cprime2_to_Cprimebar4_columns':list(map(sp,qbar)),'all9_action_intertwiners_and_nilframe8_readbacks_exact':True,'actual_ell_action_rank':q-4,'both_actual_DVR_torsion_exponents_equal1':q==6,'genuine_ell_torsion_target_identification_only_if6':q==6,'no_ell_division_or_H1_basechange':True,'actual_connecting_values_or_Nbar14_cycles_computed':False}
 save('phase4-actual-Cprime2-quotient-and-actions.json',data)
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_right_exact_Cprime_ell2_nilimage_quotient','actual_Cprime2_length':q,'actual_image_rank_in_C2':len(piv),'actual_Cprime2_ell_rank':q-4,'both_torsion_exponents_equal1':q==6,'actual_connecting_kernel_or_nil_generators_computed':False,'new_physical_source_or_jet_values':0,'new_derived_gauge_coefficient_upper':990,'whole_QQ_bound':118256,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name not in ['report.json','phase.json']},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']}))
if __name__=='__main__':main()
