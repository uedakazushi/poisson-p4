#!/usr/bin/env python3
"""Cached actual full-nine Vfirst and known-initial positive upper test."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_honest_nil_T_gauge120_cached_20261105_v63 as G
from ambient_B18_Cprime_linear_initial_source_preflight_20261105_v63 import PATHS
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,polynomial,invq,mv
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63';GATE='AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_PREFLIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(v,i,c):
 if c:
  v[i]=v.get(i,Q(0))+c
  if not v[i]:del v[i]
def add(v,w,c=Q(1)):
 v=dict(v)
 for i,a in w.items():put(v,i,c*a)
 return v
def sp(v):return [[i,c]for i,c in sorted(v.items())if c]
def pivot(cols):
 p={};chosen=[];rows=[]
 for j,v in enumerate(cols):
  z=dict(v)
  for i,w in sorted(p.items()):z=add(z,w,-z.get(i,Q(0)))
  if z:
   i=min(z);c=z[i];p[i]={k:a/c for k,a in z.items()};chosen.append(j);rows.append(i)
 return chosen,rows
def fixed(cols,n):
 cp,rr=pivot(cols);rank=len(cp);A=[[cols[j].get(i,Q(0))for j in cp]for i in rr];iv,det=invq(A)if rank else([],Q(1));ivsp=[[(i,c)for i,c in enumerate(r)if c]for r in iv];digest=hashlib.sha256()
 for v in cols:
  f=[sum((c*v.get(rr[i],Q(0))for i,c in row),Q(0))for row in ivsp];z={}
  for k,c in enumerate(f):
   if c:z=add(z,cols[cp[k]],c)
  assert z==v;digest.update(json.dumps(enc(f),separators=(',',':')).encode()+b'\n')
 return cp,rr,iv,det,digest.hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=True);assert not(out/'report.json').exists()
 def save(name,d):
  blob=json.dumps(enc(d),separators=(',',':'))+'\n';assert len(blob.encode())+sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32<<20;(out/name).write_text(blob);(out/'phase.json').write_text(json.dumps({'phase':name,'seconds':time.monotonic()-t})+'\n')
 gp=ROOT/'reports'/GATE/'report.json';gate=json.loads(gp.read_text());assert gate['all_passed']and gate['whole_saved_live_QQ_upper']==85555
 for p,h in gate['input_sha256'].items():assert sha(ROOT/p)==h,p
 peer=ROOT/'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_20261105_v63_REVIEW/report.json';assert json.loads(peer.read_text())['all_passed']
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};pins={**gate['input_sha256'],str(gp.relative_to(ROOT)):sha(gp),str(peer.relative_to(ROOT)):sha(peer)};save('preflight.json',{'input_sha256':pins,'whole_QQ_bound':85555,'no_new_source_gauge_or_higherjet':True})
 B=D['table'];base=B['basis33'];MON=[()]+[(i,)for i in range(9)];MI={m:i for i,m in enumerate(MON)};mapping={0:0}
 for i,(z,r,e)in enumerate(base):
  if B['basis_degrees'][i]!=1:continue
  var=z-1 if z else 7 if r else 8;mapping[i]=var+1
 assert set(mapping.values())==set(range(10));G.MON=MON;G.IDX=MI
 def pr33(v):
  z=[Q(0)]*10
  for i,c in enumerate(v):
   if i in mapping:z[mapping[i]]+=Q(c)
  return z
 rawroots={i:{k:Q(c)for k,c in v}for i,v in D['source33']['all_requested164_root_values']};nus=[]
 for f in D['source']['families'].values():nus.append(pr33([rawroots[f['actual_tracefree_nu']].get(i,Q(0))for i in range(33)]))
 a=G.Finite(2,nus)
 kappa=[]
 for v in D['kappa']['kappa9_C132']:
  parts=[pr33(v[mask::4])for mask in range(4)];kappa.append(tuple(parts[mask][j]for j in range(10)for mask in range(4)))
 core=[kappa[i]for i in[0,1,2,6]];av,d,x,tt=core;K=D['core5'];Dunit=a.add(a.ONE,a.neg(a.mul(d,d)));Punit=a.add(a.c(5),a.mul(d,d));units=[]
 for name in K['selected5_pivot_order']:
  root=K['five_actual_rational_roots'][name];root_terms=[[m,Q(num,den)]for m,num,den in root['numerator']];v=a.mul(polynomial(a,root_terms,core),a.power(Dunit,-root['D_power']),a.power(Punit,-root['P_power']));idx=K['source9_variable_order'].index(name);units.append(a.add(kappa[idx],a.neg(v)))
 hh=a.add(av,a.neg(d),x,a.scale(tt,-2));cc=a.scale(a.add(d,x),Q(1,2));rr=a.scale(a.add(a.scale(av,-2),d,a.neg(x)),Q(1,6));ll=a.scale(a.add(av,d,a.neg(x)),Q(1,3));variables=units+[hh,cc,rr,ll];R=D['radical'];dd=polynomial(a,R['d_polynomial'],variables);assert dd==d
 def traceunit(mask):return tuple(Q(i==mask)for i in range(40))
 vals=[]
 for pol in R['R9_generators_polynomials9']:
  v=a.mul(dd,polynomial(a,pol,variables))
  for mask in range(4):vals.append(a.mul(v,traceunit(mask)))
 g=pr33([rawroots[7790].get(i,Q(0))for i in range(33)])
 for i in range(8):
  p=[Q(0)]*40
  if i<7:p[4*(i+1)+1]=Q(-1)
  else:
   for j,c in enumerate(g):p[4*j+1]=-c
  vals.append(tuple(p))
 columns=[]
 for v in vals:
  gam=[[v[4*j+1],v[4*j+3],v[4*j+2]-Q(11,12)*v[4*j+1]]for j in range(10)];assert gam[0]==[Q(0)]*3
  columns.append({9*t+k:gam[k+1][t]for t in range(3)for k in range(9)if gam[k+1][t]})
 assert len(columns)==44;cp1,rr1,iv1,det1,hash1=fixed(columns,27);h1=27-len(cp1)
 save('phase1-actual-fullnine-first-V-columns.json',{'coordinate_order9':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],'all44_actual_linear_columns27_sparse':list(map(sp,columns)),'all44_V_origins_zero':True,'fixed_pivot_columns':cp1,'fixed_pivot_rows':rr1,'fixed_minor_inverse':iv1,'fixed_minor_determinant':det1,'all44_factorization_sha256':hash1,'actual_HF1_Cprime':h1,'same_trace_units_D_P_and_d_identity_exact':True})
 m2=list(itertools.combinations_with_replacement(range(9),2));im2={m:i for i,m in enumerate(m2)};qcols=[]
 for v in columns:
  for var in range(9):
   z={}
   for k,c in v.items():target,i=divmod(k,9);put(z,45*target+im2[tuple(sorted((i,var)))],c)
   qcols.append(z)
 cp,rr,iv,det,hsh=fixed(qcols,135);free=[i for i in range(135)if i not in rr];h2=len(free)
 save('phase2-known-quadratic-closure.json',{'all396_known_quadratic_columns135_sparse':list(map(sp,qcols)),'fixed_pivot_columns':cp,'fixed_pivot_rows':rr,'fixed_minor_inverse':iv,'fixed_minor_determinant':det,'all396_factorization_sha256':hsh,'known_quadratic_quotient_dimension':h2,'not_actual_HF2_claimed':True})
 positive=False;cubic=None
 if h2<=3:
  q=[]
  for f in free:
   v={f:Q(1)}
   for k,r in enumerate(rr):put(v,r,-sum((qcols[cp[j]].get(f,Q(0))*iv[j][k]for j in range(len(cp))),Q(0)))
   q.append(v)
  m3=list(itertools.combinations_with_replacement(range(9),3));im3={m:i for i,m in enumerate(m3)};q3=[]
  for row in range(495):
   target,mi=divmod(row,165);mon=m3[mi];divisor=45*target+im2[mon[:2]];var=mon[2];q3.append({9*f+var:v[divisor]for f,v in enumerate(q)if divisor in v})
  def project(v):
   z={}
   for row,c in v.items():z=add(z,q3[row],c)
   return z
  # Release the135-row rank scratch before the genuinely streaming phase.
  # At most27 independent compressed input columns andtheir fixedminor
  # survive; all1980+3+consistency columns are regenerated forthe factor pass.
  del qcols,iv
  def cubic_stream():
   for j,col in enumerate(columns):
    for mon in m2:
     raw={}
     for k,c in col.items():target,i=divmod(k,9);put(raw,165*target+im3[tuple(sorted(mon+(i,)))],c)
     yield ['first44_times_quadratic',j,list(mon)],raw,project(raw)
   for target in range(3):
    raw={165*target+im3[(8,8,8)]:Q(1)};yield ['ell3',target],raw,project(raw)
   for f,k in enumerate(free):
    target,mi=divmod(k,45)
    for var in range(9):yield ['consistency',f,var],{},add({9*f+var:Q(1)},q3[165*target+im3[tuple(sorted(m2[mi]+(var,)))]],-1)
  sourcehash=hashlib.sha256();piv3={};p3=[];r3=[];basis3=[];count3=0
  for index,(label,raw,col)in enumerate(cubic_stream()):
   sourcehash.update(json.dumps(enc([label,sp(raw),sp(col)]),separators=(',',':')).encode()+b'\n');count3+=1;v=dict(col)
   for i,w in sorted(piv3.items()):v=add(v,w,-v.get(i,Q(0)))
   if v:
    i=min(v);c=v[i];piv3[i]={k:a/c for k,a in v.items()};p3.append(index);r3.append(i);basis3.append(dict(col))
   if count3%256==0:save('phase3-stream-progress.json',{'columns_completed':count3,'pivot_columns':p3,'pivot_rows':r3,'selected_original_compressed_columns':list(map(sp,basis3)),'rank':len(p3)})
  assert count3==1983+9*h2 and len(basis3)<=27
  del piv3
  # Small inverse in-place, without simultaneously retaining a separate
  # minor andnormalized pivot matrix. Live cubic scratch stays<5000 QQ.
  rank3=len(p3);aug=[[v.get(i,Q(0))for v in basis3]+[Q(irow==j)for j in range(rank3)]for irow,i in enumerate(r3)];det3=Q(1)
  for j in range(rank3):
   k=next(k for k in range(j,rank3)if aug[k][j])
   if k!=j:aug[k],aug[j]=aug[j],aug[k];det3=-det3
   c=aug[j][j];det3*=c;aug[j]=[v/c for v in aug[j]]
   for i in range(rank3):
    if i!=j and aug[i][j]:
     c=aug[i][j];aug[i]=[v-c*w for v,w in zip(aug[i],aug[j])]
  assert all(aug[i][j]==Q(i==j)for i in range(rank3)for j in range(rank3))
  iv3=[row[rank3:]for row in aug];del aug;factorhash=hashlib.sha256()
  for label,raw,v in cubic_stream():
   factors=mv(iv3,[v.get(i,Q(0))for i in r3]);back={}
   for col,c in zip(basis3,factors):back=add(back,col,c)
   assert back==v;factorhash.update(json.dumps(enc(factors),separators=(',',':')).encode()+b'\n')
  zero3=len(p3)==9*h2;upper=3+h1+h2
  if zero3:assert upper>=9,'contradiction to independent genuine length>=9'
  positive=zero3 and upper==9
  cubic={'quotient_projection_rows3_by135_sparse':list(map(sp,q)),'all495_canonical_projection_columns_sparse':list(map(sp,q3)),'all_cubic_source_and_reduction_stream_sha256':sourcehash.hexdigest(),'source_count1980_plus3ell3_plus_consistency':count3,'compressed_target_dimension':9*h2,'fixed_pivot_columns':p3,'fixed_pivot_rows':r3,'selected_original_compressed_columns':list(map(sp,basis3)),'fixed_minor_inverse':iv3,'fixed_minor_determinant':det3,'all_compressed_column_factorization_sha256':factorhash.hexdigest(),'known_third_quotient_zero':zero3,'cumulative_known_upper_through2':upper,'positive_length9_bound_proved':positive,'all_compressed_columns_never_materialized':True}
  save('phase3-known-cubic-ell3-upper-test.json',cubic)
 report={'id':STEM,'all_passed':True,'status':'PROVED_cached_Vfirst_and_known_initial_bounds','actual_Cprime_HF1':h1,'known_quadratic_quotient_dimension':h2,'positive_Cprime_ell3_length9_proved':positive,'actual_long_torsion_exponent_equals2_if_positive':positive,'no_negative_exponent_or_actual_HF2_from_failed_upper':True,'new_source_gauge_or_higherjet':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name not in['report.json','phase.json']},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']}))
if __name__=='__main__':main()
