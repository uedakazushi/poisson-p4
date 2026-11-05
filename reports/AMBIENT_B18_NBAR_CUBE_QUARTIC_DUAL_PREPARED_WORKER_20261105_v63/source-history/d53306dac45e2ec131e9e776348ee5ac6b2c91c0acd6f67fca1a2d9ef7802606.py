#!/usr/bin/env python3
"""Preparedtiny3x8cachedquarticrelationrank; no cubeproducts."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];GATE='AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_PREFLIGHT_20261105_v63';AUG='AMBIENT_B18_NBAR_CUBE_QUARTIC_ZERO3_SOURCE_20261105_v63';STEM='AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def count(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def main(peer_path):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={};custody={}
 def read(path,keys=None):
  p=ROOT/path;pins[path]=sha(p)
  if keys is None:return json.loads(p.read_text())
  text=p.read_text();dec=json.JSONDecoder();result={}
  for key in keys:
   token=json.dumps(key);assert text.count(token)==1;i=text.index(':',text.index(token))+1
   while text[i].isspace():i+=1
   result[key],j=dec.raw_decode(text,i);custody[key]=[i,j]
  return result
 g=f'reports/{GATE}/report.json';h=f'reports/{AUG}/report.json';r=read(g);s=read(h);peer=read(peer_path);assert all(d['all_passed']for d in[r,s,peer]);assert peer['input_sha256'][g]==sha(ROOT/g)and peer['input_sha256'][h]==sha(ROOT/h)
 for rep in[r,s]:
  for p,v in rep['input_sha256'].items():assert sha(ROOT/p)==v;pins[p]=v
 ep=f'reports/{GATE}/tiny-five-relations-and-constructive-quartic-dual-source.json';D=read(ep,['all1650_constructive_quartic_reduction_columns3_sparse','source_q3_restricted_to_five_pure_r3_heads','source_actual_natural_origin_columns5','known_pure_r4_basis_rows1650']);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 ap=f'reports/{AUG}/three-actual-zero-heads-and-corrected-first-source-relations.json';Z=read(ap);assert sha(ROOT/ap)==s['evidence_sha256'][Path(ap).name]
 A=[[Q(c)for c in row]for row in D['source_q3_restricted_to_five_pure_r3_heads']];M=[[Q(c)for c in col]for col in D['source_actual_natural_origin_columns5']];p4=[{i:Q(c)for i,c in col}for col in D['all1650_constructive_quartic_reduction_columns3_sparse']]
 values=[[sum((A[i][k]*col[k]for k in range(5)),Q(0))for i in range(3)]for col in M];m4=list(itertools.combinations_with_replacement(range(8),4));idx={m:i for i,m in enumerate(m4)}
 for col in Z['three_corrected_actual_M1_columns40']:
  v=[Q(0)]*3
  for key,c in col.items():target,var=divmod(int(key),8);source=330*target+idx[tuple(sorted((var,7,7,7)))];v=[x+Q(c)*p4[source].get(i,Q(0))for i,x in enumerate(v)]
  values.append(v)
 assert len(values)==8;piv={};chosen=[]
 for j,col in enumerate(values):
  v=list(col)
  for row,w in sorted(piv.items()):c=v[row];v=[a-c*b for a,b in zip(v,w)]
  if any(v):row=next(i for i,c in enumerate(v)if c);c=v[row];piv[row]=[a/c for a in v];chosen.append(j)
 rank=len(chosen);assert rank<=2,('CONTRADICTS_PROVED_ACTUAL_HF4_1',values)
 dual=None;phi=None;det=None;factor=[];rows=[]
 if rank==2:
  u,v=[values[j]for j in chosen];cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];f=next(i for i,c in enumerate(cross)if c);dual=[c/cross[f]for c in cross];assert all(sum((a*b for a,b in zip(dual,col)),Q(0))==0 for col in values)
  rows=list(next((i,j)for i,j in itertools.combinations(range(3),2)if u[i]*v[j]-u[j]*v[i]));i,j=rows;det=u[i]*v[j]-u[j]*v[i];inverse=[[v[j]/det,-v[i]/det],[-u[j]/det,u[i]/det]];minor=[[u[i],v[i]],[u[j],v[j]]]
  assert all(sum((minor[i][k]*inverse[k][j]for k in range(2)),Q(0))==Q(i==j)for i in range(2)for j in range(2))and all(sum((inverse[i][k]*minor[k][j]for k in range(2)),Q(0))==Q(i==j)for i in range(2)for j in range(2))
  for col in values:
   co=[sum((c*col[rows[k]]for k,c in enumerate(row)),Q(0))for row in inverse];assert [co[0]*u[k]+co[1]*v[k]for k in range(3)]==col;factor.append(co)
  phi=[[j,sum((dual[i]*col.get(i,Q(0))for i in range(3)),Q(0))]for j,col in enumerate(p4)];phi=[p for p in phi if p[1]];assert [dict(phi).get(j,Q(0))for j in D['known_pure_r4_basis_rows1650']]==dual
 data={'all_eight_actual_quartic_relation_vectors3':values,'first5_from_M0_r4_last3_from_corrected_ZERO3_M1_r3':True,'actual_source_relation_rank':rank,'selected_independent_source_columns':chosen,'fixed2minor_rows':rows,'fixed2minor_determinant':det,'all8_fixed_factorizations':factor,
  'actual_complete_quartic_dual_on_three_r4_representatives':dual,'actual_complete_quartic_dual1650_sparse':phi,'actual_quartic_dual_proved':rank==2,
  'proof':'Constructivesourcep4givesactualsurjectionfromthreepure-r4classes. Rank2actualrelationslieinkernel;knownactualHF4=1forcesEXACTkernel2. Thusdual*p4isactualcompletequarterdual. No assumptionp4iswell-definedbeforethisdimensionargument.',
  'rank_less2_is_inconclusive':True,'no_cube_product_values_computed':True,'selective_field_custody':custody}
 owned=count([r,s,peer,D,Z,A,M,p4,values,piv,dual,phi,data]);assert owned<=24000,owned
 def guard(x):
  if isinstance(x,Q):assert max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=2048
  elif isinstance(x,(list,tuple)):
   for y in x:guard(y)
  elif isinstance(x,dict):
   for y in x.values():guard(y)
 guard(data);blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert sum((ROOT/p).stat().st_size for p in pins)+len(blob.encode())<32<<20;ev=out/'actual-eight-relations-and-complete-quartic-dual.json';ev.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_cached_tiny_quartic_relation_test','actual_relation_rank':rank,'actual_complete_quartic_dual_proved':rank==2,'no_cube_values_or_new_source_gauge_jet':True,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'owned_numeric_slots':owned,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source-peer-report',required=True);main(p.parse_args().source_peer_report)
