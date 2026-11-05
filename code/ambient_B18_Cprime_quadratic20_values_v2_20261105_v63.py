#!/usr/bin/env python3
"""Prepared cache-only B27/30-gauge/complete4x20 job; no physical solve."""
import argparse,collections,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_honest_nil_T_gauge120_cached_20261105_v63 as G
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_Cprime_quadratic20_source_preflight_20261105_v63 import PATHS
from ambient_B18_Cprime_linear_initial_values_20261105_v63 import fixed,sp,rec
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,polynomial
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63';GATE='AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE_PREFLIGHT_20261105_v63';SUB='AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE2_SUBSET_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def qq(v):
 if isinstance(v,Q):return 1
 if isinstance(v,dict):return sum(map(qq,v.values()))
 if isinstance(v,(list,tuple)):return sum(map(qq,v))
 return 0
def guard(v):
 if isinstance(v,Q):assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=2048
 elif isinstance(v,dict):
  for x in v.values():guard(x)
 elif isinstance(v,(tuple,list)):
  for x in v:guard(x)
 return v
def algebra(m=0,nus=()):
 a=G.Finite(m,nus)
 for name in ['add','neg','scale','mul2','dual']:
  original=getattr(a,name)
  def checked(*args,_f=original):return guard(_f(*args))
  setattr(a,name,checked)
 return a
def main(peerpath):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={};custody={};saved=0;peak=0
 def record(path):
  p=ROOT/path if not Path(path).is_absolute()else Path(path);pins[str(p.relative_to(ROOT))]=sha(p);return p
 def read(path,keys=None):
  p=record(path)
  if keys is None:return json.loads(p.read_text())
  text=p.read_text();decoder=json.JSONDecoder();d={};offsets={}
  for k in keys:
   token=json.dumps(k);assert text.count(token)==1,(p,k);i=text.index(':',text.index(token))+1
   while text[i].isspace():i+=1
   d[k],j=decoder.raw_decode(text,i);offsets[k]=[i,j]
  custody[str(p.relative_to(ROOT))]={'parsed_fields_only':offsets,'whole_JSON_not_decoded':True};return d
 def save(name,d,live=(),envelope=0):
  nonlocal saved,peak
  guard(d);n=qq(d);peak=max(peak,saved+n+qq(live)+envelope);assert peak<=118472,('WHOLE_OWNED_QQ_GUARD',name,peak)
  b=json.dumps(enc(d),separators=(',',':'))+'\n';raw=sum((ROOT/p).stat().st_size for p in pins);assert raw+len(b.encode())+sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32<<20
  (out/name).write_text(b);saved+=n;(out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ':saved,'owned_peak_with_named_live_and_envelope':peak,'seconds':time.monotonic()-t})+'\n')
 report=read(f'reports/{GATE}/report.json');assert report['all_passed']
 gate=read(f'reports/{GATE}/complete-four-by20-and-ell-square-cache-source-gate.json');assert gate['whole_saved_live_QQ_upper']==118472
 peer=read(peerpath);assert peer['all_passed'];assert peer['input_sha256'][f'reports/{GATE}/report.json']==sha(ROOT/f'reports/{GATE}/report.json')
 for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h,p;record(p)
 subrep=read(f'reports/{SUB}/report.json');assert subrep['all_passed']and subrep['whole_FULL9_JSON_never_decoded'];subpath=f'reports/{SUB}/minimal-two57-source2-and-labels.json';assert sha(ROOT/subpath)==subrep['evidence_sha256']['minimal-two57-source2-and-labels.json'];J=read(subpath)
 # Never decode the310866-leaf FULL9 object; only its immutable subset isused.
 keys={
  'source':['DAG_scalar_nodes','families'],
  'interface':['actual_group15_positions','retained_source9_indices'],
  'proper':['zero_coordinate_rows36','CA_constant_left_inverse_rows19','CA_fiber10_fixed_coefficient_indices','CA_free_fiber9_coefficient_indices','actual_CA_source19_center','CA_left_inverse19'],
  'second5':['source9_variable_order','selected5_pivot_order','five_actual_rational_roots'],
  'J30':['J30_inverse','det_J30'],
  'firstpreflight':['GL18_source_parameter_positions'],
  'source33':['joint105_unknown_order','all105_unknown_final_values','all_requested164_root_values'],
  'Vfirst':['all44_actual_linear_columns27_sparse'],
  'gauge33':['actual_even15_B33','actual_odd15_B33'],
  'proper33':['actual14_C132','kappa9_C132'],
  'adapter2':['YB_in_YA_linear','YB_in_YA_quadratic'],
  'radical':['R9_generators_polynomials9','d_polynomial']}
 D={k:read(PATHS[k],fields)for k,fields in keys.items()}
 save('preflight.json',{'input_sha256':pins,'selective_field_custody':custody,'FULL9_source_custody_report_sha256':sha(ROOT/f'reports/{SUB}/report.json'),'owned_QQ_cap':118472,'raw_total_cap':32<<20,'new_physical_source_solve':False})
 MON=[tuple(m)for m in gate['base27_monomials']];MI={m:i for i,m in enumerate(MON)};G.MON=MON;G.IDX=MI;b=algebra();S=D['source'];nodes=S['DAG_scalar_nodes'];y=[tuple(Q(i==MI[(v,)])for i in range(27))for v in range(9)];YA=y[:8]+[b.add(y[7],b.scale(y[8],Q(1,5)))];YB=[]
 def compose(values,args):
  z=b.Z
  for mon,c in zip(J['source2_monomials'],values):
   if Q(c):z=b.add(z,b.scale(b.mul(*(args[i]for i in mon)),Q(c)))
  return z
 for lin,quad in zip(D['adapter2']['YB_in_YA_linear'],D['adapter2']['YB_in_YA_quadratic']):
  z=b.Z
  for i,c in enumerate(lin):
   if Q(c):z=b.add(z,b.scale(YA[i],Q(c)))
  for exp,c in quad:
   if Q(c):z=b.add(z,b.scale(b.mul(*(YA[i]for i,e in enumerate(exp)for _ in range(e))),Q(c)))
  YB.append(z)
 own={f:{i:compose(v,YA if f=='CA'else YB)for i,v in J['both57_source2'][f]['all57_source2jets']}for f in['CA','CB']};shared=set(own['CA'])&set(own['CB']);assert len(shared)==9 and all(own['CA'][i]==own['CB'][i]for i in shared);first=own['CA']|own['CB'];assert len(first)==105
 maps=[rec(v)for v in gate['actual33_to26_columns']]
 def old26(v):
  z=[Q(0)]*27
  for i,c in v.items():
   for j,d in maps[i].items():z[j]+=c*d
  return tuple(z)
 old105={i:rec(v)for i,v in zip(D['source33']['joint105_unknown_order'],D['source33']['all105_unknown_final_values'])};assert set(first)==set(old105)and all(first[i][:26]==old26(old105[i])[:26]for i in first)
 roots=[i for f in S['families'].values()for p in f['complete_normalized_physical_w_coefficients'][0]for i in p]+[f['actual_tracefree_nu']for f in S['families'].values()]+[i for f in['CA','CB']for i in J['both57_source2'][f]['equations57']]+[7790,28539];needed=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if i in first:continue
  n=nodes[i]
  if n[0]in ['a','m','n','i']:todo.extend(n[1:])
  else:assert n[0]=='c'
 uses=collections.Counter(roots)
 for i in needed:
  if i not in first and nodes[i][0]in ['a','m','n','i']:uses.update(nodes[i][1:])
 vals={};live=streampeak=0
 for i in sorted(needed):
  n=nodes[i];op=n[0]
  if i in first:v=first[i]
  elif op=='c':v=b.c(Q(n[1],n[2]))
  elif op=='a':v=b.add(vals[n[1]],vals[n[2]])
  elif op=='m':v=b.mul(vals[n[1]],vals[n[2]])
  elif op=='n':v=b.neg(vals[n[1]])
  else:assert op=='i';v=b.inv(vals[n[1]])
  guard(v);vals[i]=v;live+=27;streampeak=max(streampeak,live);assert streampeak<=gate['budget']['stream_actual_physical_DAG_peak']
  if i not in first and op in ['a','m','n','i']:
   for j in n[1:]:
    uses[j]-=1
    if not uses[j]:del vals[j];live-=27
 residual114=[]
 for f,args in [('CA',YA),('CB',YB)]:
  eq=[vals[i]for i in J['both57_source2'][f]['equations57']];res=list(eq[:48])+[b.add(v,b.neg(w))for v,w in zip(eq[48:],args)];assert res==[b.Z]*57;residual114.extend([[i,[]]for i in J['both57_source2'][f]['equations57']])
 overlap=0
 for i,v in D['source33']['all_requested164_root_values']:
  if i in vals:assert vals[i][:26]==old26(rec(v))[:26];overlap+=1
 assert overlap==164
 pairs={f:[[vals[x],vals[y]]for x,y in data['complete_normalized_physical_w_coefficients'][0]]for f,data in S['families'].items()};nu=[vals[f['actual_tracefree_nu']]for f in S['families'].values()];gA=vals[7790]
 save('phase1-common105-and-114-source-equations.json',{'base27_monomials':MON,'all105_cached_composite_values':[[i,v]for i,v in sorted(first.items())],'all114_original_residuals':residual114,'all9_shared_point_coefficients_match':True,'all105_B26_source_cache_coefficients_match':True,'all_old164_root_overlap_count':overlap,'nu_A_B27':nu,'gA_source7790_B27':gA,'all600_physical_root_fingerprints':{f:hashlib.sha256(json.dumps(enc(v),separators=(',',':')).encode()).hexdigest()for f,v in pairs.items()},'source_stream_peak_QQ':streampeak},live=[first,pairs],envelope=8000);del vals,own,YB,J,old105;del D['source33'],D['adapter2']
 ev=[old26({i:Q(c)for i,c in enumerate(v)if Q(c)})for v in D['gauge33']['actual_even15_B33']];od=[old26({i:Q(c)for i,c in enumerate(v)if Q(c)})for v in D['gauge33']['actual_odd15_B33']];I,P=D['interface'],D['proper'];pos=I['actual_group15_positions'];GL=D['firstpreflight']['GL18_source_parameter_positions'];rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];fixedfiber=P['CA_fiber10_fixed_coefficient_indices'];freefiber=P['CA_free_fiber9_coefficient_indices'];centers=list(map(Q,P['actual_CA_source19_center']))
 def family_frame(a,f,mask):
  a.install();s=S['families'][f];p=H.add(a.tau(mask),H.scale(a.dual(first[s['p_quadratic_relation']['s']]),Q(1,2)));gg=H.ident(5)
  for name,(x,z)in s['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':i,j=divmod(GL[int(name[5:])-1],5);gg[i][j]=H.add(gg[i][j],a.dual(first[x]),H.mul(a.dual(first[z]),p))
  core=H.ident(5)
  if f=='CA':core[0][3]=H.scale(p,-2);core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
  else:core[0][3]=H.scale(p,Q(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),Q(1,4))
  return H.mm(gg,core)
 def physical(a,f,mask):a.install();return [H.add(a.dual(x),H.mul(a.dual(z),a.tau(mask)))for x,z in pairs[f]]
 def gauge(a):
  a.install();gg=H.ident(5)
  for(i,j),x,z in zip(pos,ev,od):gg[i][j]=H.add(gg[i][j],a.dual(x),H.mul(a.dual(z),a.tau(0)))
  return gg
 del D['gauge33']
 a54=algebra(1,[nu[0]]);HA=family_frame(a54,'CA',0);PA=physical(a54,'CA',0)
 def constraints():
  a54.install();gg=gauge(a54);fa=H.frame(H.mm(gg,HA));tr=H.mm(fa,gg);moved=H.physical_push(PA,tr,H.mi(tr));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in ar]);return [fa[1][0],fa[1][3],fa[1][4],fa[0][3],fa[0][4]]+[H.add(own[j][0],H.c(-centers[j]))for j in fixedfiber]
 residual=constraints();assert all(not v[2*i+parity]for v in residual for i in range(26)for parity in range(2));rhs=[v[2*26+parity]for parity in range(2)for v in residual];Ji=[[Q(c)for c in row]for row in D['J30']['J30_inverse']];correction=[-sum((c*v for c,v in zip(row,rhs)),Q(0))for row in Ji]
 for i in range(15):v=list(ev[i]);v[26]+=correction[i];ev[i]=tuple(v);v=list(od[i]);v[26]+=correction[15+i];od[i]=tuple(v)
 assert all(v==a54.Z for v in constraints());save('phase2-thirty-ell-square-gauge-values.json',{'RHS30':rhs,'correction30':correction,'actual_even15_B27':ev,'actual_odd15_B27':od,'all30_constraints_exact':True,'all780_old_B26_gauge_coefficients_preserved':True,'only_new_gauge_pair':[8,8]},live=[pairs,first,ev,od],envelope=43200);del HA,PA,residual
 a=algebra(2,nu);gg=gauge(a);HA=family_frame(a,'CA',0);assert H.frame(H.mm(gg,HA))==H.ident(5);HB=family_frame(a,'CB',1);FB=H.frame(H.mm(gg,HB));PA=physical(a,'CA',0);moved=H.physical_push(PA,gg,H.mi(gg));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in ar]);own=[v[0]for v in own];assert all(own[j]==H.c(centers[j])for j in fixedfiber)
 actual=[H.add(own[j],H.c(-centers[j]))for j in freefiber]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j]for j in I['retained_source9_indices']]
 for key,vs in [('actual14_C132',actual),('kappa9_C132',kappa)]:
  for old,new in zip(D['proper33'][key],vs):
   for mask in range(4):assert tuple(new[4*i+mask]for i in range(26))==old26({i:Q(c)for i,c in enumerate(old[mask::4])if Q(c)})[:26]
 save('phase3-actual-proper-B27-values.json',{'actual14_C108':actual,'kappa9_C108':kappa,'all_old_B26_proper_coefficients_exact':True,'all_CA_frame_and_fiber_constraints_exact':True},live=[pairs,first,ev,od,kappa],envelope=43200);del gg,HA,HB,FB,PA,moved,own,pairs;del D['proper33']
 # Literal NFsupport removes every nonplane coefficient before this source
 # projection; no discarded coefficient is filled withzero in the full ring.
 plane=[tuple(m)for m in gate['required_plane6_monomials']];ix=[MI[m]for m in plane];G.MON=plane;G.IDX={m:i for i,m in enumerate(plane)};nu6=[[v[i]for i in ix]for v in nu];a6=algebra(2,nu6);kap=[tuple(v[4*i+mask]for i in ix for mask in range(4))for v in kappa];core=[kap[i]for i in[0,1,2,6]];av,d,x,tt=core;K=D['second5'];du=a6.add(a6.ONE,a6.neg(a6.mul(d,d)));pu=a6.add(a6.c(5),a6.mul(d,d));u5=[]
 for name in K['selected5_pivot_order']:
  r=K['five_actual_rational_roots'][name];terms=[[m,Q(n,z)]for m,n,z in r['numerator']];v=a6.mul(polynomial(a6,terms,core),a6.power(du,-r['D_power']),a6.power(pu,-r['P_power']));j=K['source9_variable_order'].index(name);u5.append(a6.add(kap[j],a6.neg(v)))
 hh=a6.add(av,a6.neg(d),x,a6.scale(tt,-2));cc=a6.scale(a6.add(d,x),Q(1,2));rr=a6.scale(a6.add(a6.scale(av,-2),d,a6.neg(x)),Q(1,6));ll=a6.scale(a6.add(av,d,a6.neg(x)),Q(1,3));vs=u5+[hh,cc,rr,ll];R=D['radical'];dd=polynomial(a6,R['d_polynomial'],vs);assert dd==d
 def unit(mask):return tuple(Q(i==mask)for i in range(24))
 V=[]
 for pol in R['R9_generators_polynomials9']:
  v=a6.mul(dd,polynomial(a6,pol,vs))
  for mask in range(4):V.append(a6.mul(v,unit(mask)))
 for i in range(8):
  v=[Q(0)]*24
  if i==7:
   for j,k in enumerate(ix):v[4*j+1]=-gA[k]
  V.append(tuple(v))
 first=[rec(v)for v in D['Vfirst']['all44_actual_linear_columns27_sparse']];raw=[];mons9=list(__import__('itertools').combinations_with_replacement(range(9),2));p4=[rec(v)for v in gate['actual_NF4_rows_by135_sparse']];nf=[]
 for col,v in enumerate(V):
  gam=[[v[4*j+1],v[4*j+3],v[4*j+2]-Q(11,12)*v[4*j+1]]for j in range(6)];assert gam[0]==[Q(0)]*3
  assert all(gam[j][tgt]==first[col].get(9*tgt+var,Q(0))for j,var in[(1,7),(2,8)]for tgt in range(3));raw.append([gam[j][tgt]for tgt in range(3)for j in[3,4,5]])
  poly={45*tgt+mons9.index(plane[j]):gam[j][tgt]for tgt in range(3)for j in[3,4,5]if gam[j][tgt]};nf.append([sum((c*poly.get(i,Q(0))for i,c in row.items()),Q(0))for row in p4])
 kernel=[rec(v)for v in gate['all20_constant_kernel_columns44_sparse']];values=[[sum((c*nf[j][f]for j,c in col.items()),Q(0))for f in range(4)]for col in kernel];cvals=[{i:c for i,c in enumerate(v)if c}for v in values];cp,rr,iv,det,digest=fixed(cvals,4);rank=len(cp);hf2=4-rank
 save('phase4-all-twenty-actual-quadratic-values.json',{'all44_actual_V2_pureplane_nine_values':raw,'all44_projected_NF4_values':nf,'all20_actual_quadratic_values4':values,'source_constant_kernel20':gate['all20_constant_kernel_columns44_sparse'],'fixed_rank_columns':cp,'fixed_rank_rows':rr,'fixed_minor_inverse':iv,'fixed_minor_determinant':det,'all20_factorization_sha256':digest,'actual_complete_rank':rank,'actual_Cprime_HF2':hf2,'actual_HF_through2':[3,3,hf2],'no_long_exponent_from_this_rank_alone':True},live=[V,nf,values],envelope=8000)
 result={'id':STEM,'all_passed':True,'status':'PROVED_complete_actual4x20_cached_quadratic_values','actual_complete_rank':rank,'actual_Cprime_HF_through2':[3,3,hf2],'new_derived_gauge_coefficients':30,'new_physical_source_or_source3_values':0,'all114_source_and30_gauge_constraints_exact':True,'all_old_B26_source_and_proper_coefficients_match':True,'actual_long_exponent_or_nilcube_decided':False,'saved_QQ':saved,'owned_QQ_peak_bound':peak,'input_sha256':pins,'selective_field_custody':custody,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name not in['report.json','phase.json']},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256','selective_field_custody']}))
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--source-peer-report',required=True);main(parser.parse_args().source_peer_report)
