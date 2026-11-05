#!/usr/bin/env python3
"""Independent target17 Cramer coefficients, source and unit/sign review."""
import argparse,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_prePOWER_kappa17_packed_review_20261105_v63 import sha,shape,combine,pack
ROOT=Path(__file__).resolve().parents[1]
A='S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63';GATE=A+'_PACKED_REVIEW_SOURCE_GATE';REVIEW=A+'_REVIEW';CAP=73322
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,separators=(',',':'))+'\n')
def cv(b):return {('denominator' if k=='D' else 'integer_l1' if k=='L1' else k):v for k,v in b.items()}
def skip(s,i):
 while s[i].isspace():i+=1
 if s[i] not in '[{':
  if s[i]=='"':
   j=i+1
   while True:
    if s[j]=='\\':j+=2
    elif s[j]=='"':return j+1
    else:j+=1
  j=i
  while j<len(s) and s[j] not in ',]}':j+=1
  return j
 stack=[s[i]];j=i+1
 while stack:
  if s[j]=='"':j=skip(s,j);continue
  if s[j] in '[{':stack.append(s[j])
  elif s[j] in ']}':stack.pop()
  j+=1
 return j
def leaf(s,key,indices):
 tag=json.dumps(key)+':';assert s.count(tag)==1;i=s.index(tag)+len(tag)
 for wanted in indices:
  while s[i].isspace():i+=1
  assert s[i]=='[';i+=1
  for k in range(wanted):
   i=skip(s,i)
   while s[i].isspace():i+=1
   assert s[i]==',';i+=1
 return json.JSONDecoder().raw_decode(s,i)[0]
def asdict(records):
 out={tuple(m):Q(c) for m,c in records};assert len(out)==len(records) and all(out.values());return out
def universal():
 zero=(0,)*9
 def add(a,b,c=1):
  d=dict(a)
  for m,v in b.items():
   q=d.get(m,0)+c*v
   if q:d[m]=q
   else:d.pop(m,None)
  return d
 def mul(a,b):
  d={}
  for m,v in a.items():
   for n,w in b.items():d=add(d,{tuple(x+y for x,y in zip(m,n)):v*w})
  return d
 M=[]
 for i in range(3):
  row=[]
  for j in range(3):m=list(zero);m[3*i+j]=1;row.append({tuple(m):1})
  M.append(row)
 cof=[]
 for i in range(3):
  row=[]
  for j in range(3):
   rr=[k for k in range(3) if k!=i];cc=[k for k in range(3) if k!=j];v=add(mul(M[rr[0]][cc[0]],M[rr[1]][cc[1]]),mul(M[rr[0]][cc[1]],M[rr[1]][cc[0]]),-1);row.append({m:(-1)**(i+j)*c for m,c in v.items()})
  cof.append(row)
 determinant={}
 for i in range(3):determinant=add(determinant,mul(M[i][0],cof[i][0]))
 for i in range(3):
  for j in range(3):
   v={}
   for k in range(3):v=add(v,mul(M[i][k],cof[j][k]))
   assert v==(determinant if i==j else {})
 return True
def custody():
 pp=ROOT/'reports'/A/'report.json';P=read(pp);assert P['all_passed'];pins={str(pp.relative_to(ROOT)):sha(pp)}
 for p,h in P['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 for p,h in P['evidence_sha256'].items():
  if p.startswith('Cramer') or p=='full-source-and-adjugate-proof.json':path=pp.parent/p;assert sha(path)==h;pins[str(path.relative_to(ROOT))]=h
 sp='reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_SOURCE_20261105_v63/three-cramer-schedule.json';cp='reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_WORKER_GATE_20261105_v63/literal-source-custody.json';schedule=read(ROOT/sp);C=read(ROOT/cp);pins[sp]=sha(ROOT/sp);pins[cp]=sha(ROOT/cp)
 assert C['physical_row_order']==[2,4,5] and C['physical_columns']==[27,31,22] and C['target']==17 and C['root_order']==[25,39,53]
 A3=C['input_A3_and_v17']['A3'];V=C['input_A3_and_v17']['v17'];assert len(A3)==len(V)==3
 # Decode only six heavy polynomials and three row units, not whole H5 arrays.
 hp='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json';s=(ROOT/hp).read_text();scales=[Q(40000,9),Q(160000,63),Q(-800000,81)]
 for i,row in enumerate([1,3,4]):
  rho=leaf(s,'known_row_units_rho14',[row]);assert Q(rho['constant'])==scales[i] and rho['a_d_F_B_C_J_powers']==[1,0,1,0,2,0]
  assert C['source_normalizations'][i]['current_row']==[2,4,5][i] and C['source_normalizations'][i]['compact_row']==row
  for j,col in enumerate([0,4]):assert asdict(A3[i][j])=={m:q*scales[i] for m,q in asdict(leaf(s,'raw_compact14x5',[row,col])).items()}
 s=None
 bp='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json';s=(ROOT/bp).read_text()
 for i in range(3):assert asdict(A3[i][2])==asdict(leaf(s,'bottom4x6',[i,0]))
 s=None
 for i in range(3):
  p=f'reports/S61_NONSQUARE_PREPOWER_KAPPA17_ENTRY_20261105_v63/source-root{i}.json';r=read(ROOT/p);assert asdict(V[i])=={tuple(m):Q(int(c),int(r['denominator'])) for m,c in r['integer_coefficients']}
 # Freeze existing inherited unit and the sign already carried by kappa.
 for stem in ['S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63_REVIEW','S61_NONSQUARE_PREPOWER_KAPPA17_ENTRY_20261105_v63_REVIEW']:
  p=f'reports/{stem}/report.json';r=read(ROOT/p);assert r['all_passed'];pins[p]=sha(ROOT/p)
  if stem.startswith('S61_NONSQUARE_PREPOWER_KAPPA11'):assert r['no_extra_kappa3_p_R_h22_Delta5_inverse']
 op='reports/S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63/operator.json';k11=read(ROOT/op);assert k11['identities'][0]=='kappa3*detA8=cprod*Saux*T';pins[op]=sha(ROOT/op)
 staged='reports/S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63/operator.json';st=read(ROOT/staged);assert st['nodes'][22]['last_row_equals']=='minus kappa';pins[staged]=sha(ROOT/staged)
 assert scales[0]*scales[1]*scales[2]==Q(-5120000000000000,45927)
 u='reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_UNIVERSAL_SOURCE_PROOF_20261105_v63/report.json';ur=read(ROOT/u);assert ur['all_passed'] and ur['free_nine_variable_polynomial_identity_Aadj_equals_detI'];pins[u]=sha(ROOT/u);assert universal()
 return P,schedule,C,A3,V,pins
def source():
 t=time.monotonic();P,S,C,A3,V,pins=custody();ab=[[shape(p) for p in row] for row in A3];vb=list(map(shape,V));assert [[cv(b) for b in row] for row in S['A3_bounds']]==ab and list(map(cv,S['v17_bounds']))==vb
 allbounds=[b for row in ab for b in row]+vb;ops=0;outputs=[];maxco=0
 for j in range(3):
  acc=None
  for i in range(3):
   step=S['steps'][3*j+i];rr=[k for k in range(3) if k!=i];cc=[k for k in range(3) if k!=j];assert step['minor_rows']==rr and step['minor_columns']==cc and step['cofactor_sign']==(-1)**(i+j)
   p=combine(ab[rr[0]][cc[0]],ab[rr[1]][cc[1]],True);q=combine(ab[rr[0]][cc[1]],ab[rr[1]][cc[0]],True);co=combine(p,q);term=combine(co,vb[i],True);acc=combine(acc,term)
   for k,b in [('first_product',p),('second_product',q),('cofactor',co),('term_product',term),('accumulator',acc)]:assert b==cv(step[k]);allbounds.append(b)
   maxco=max(maxco,co['support']);ops+=ab[rr[0]][cc[0]]['support']*ab[rr[1]][cc[1]]['support']+ab[rr[0]][cc[1]]['support']*ab[rr[1]][cc[0]]['support']+co['support']*vb[i]['support']
  outputs.append(acc)
 assert outputs==list(map(cv,S['forecast']['output_bounds']))
 assert all(b['box'][0][1]-b['box'][0][0]<64 and b['box'][1][1]<64 and min(b['box'][1][0],b['box'][2][0])>=0 for b in allbounds)
 bits=max(max(b['denominator'].bit_length(),b['integer_l1'].bit_length()) for b in allbounds);assert bits<434
 rawinput=4*sum(len(p) for row in A3 for p in row)+4*sum(map(len,V));parsed=2*S['forecast']['original_input_terms']+24;meta=15000;largest=max(b['support'] for b in outputs)
 # Packed sparse cofactor and direct output accumulator; one comparison JSON.
 peaks={'custody_literal_plus_parsed':rawinput+parsed+meta+4000,'source_heavy_oneleaf':rawinput+4*157+meta+4000,'direct_cofactor_accumulator':parsed+2*maxco+2+2*largest+2+meta,'readback_one_output':parsed+2*largest+2+4*1768+2+meta}
 whole=max(peaks.values());assert whole<=CAP<150000 and ops<4_000_000
 external=sum((ROOT/p).stat().st_size for p in pins);raw=external+300000;assert raw<5820000
 out=ROOT/'reports'/GATE;out.mkdir(exist_ok=False);ep=out/'independent-three-Cramer-custody-and-packed-source-bound.json';dump(ep,{'all12_literal_sources_and_rowunits_verified':True,'54_typed_formula_structure':'12atoms;9*(2products+cofactor difference+termproduct)+6accumulator adds, with roots25,39,53',
 'universal9_variable_Aadj_identity_independently_verified':True,'A3_bounds':ab,'v17_bounds':vb,'output_bounds':outputs,'integer_bits_bound':bits,'integer_products_bound':ops,'ownership_phases':peaks,'whole_numeric_bound':whole,'old_cap_not_raised':CAP,'TOTAL_RAW_upper':raw,
 'sign_unit':'kappa=-lastadjrow(A4), hence detA3=-kappa3. Accepted kappa3*detA8=cprod*Saux*T with FULLU/Delta407/Saux/T units, no new inverse. Cramers adjA3*v have sign(-1)^(i+j) only; x=-Y/kappa3.',
 'no_old_Q2_or_checkpoint_read':True,'signed_Cramer_values_not_computed':True})
 report={'id':GATE,'all_passed':True,'source_only':True,'whole_numeric_bound':whole,'integer_bits_bound':bits,'TOTAL_RAW_upper':raw,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{GATE}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(out/'report.json',report);print(json.dumps({'source_gate_passed':True,'owned':whole,'bits':bits,'products_bound':ops,'report_sha256':sha(out/'report.json')}))
def verify():
 t=time.monotonic();gp=ROOT/'reports'/GATE/'report.json';g=read(gp);assert g['all_passed'] and g['code_sha256']==sha(Path(__file__))
 for p,h in g['input_sha256'].items():assert sha(ROOT/p)==h
 D=read(gp.parent/'independent-three-Cramer-custody-and-packed-source-bound.json');P,S,C,A3,V,pins=custody();polys=[]
 def parse(records,b):
  den=b['denominator'];v={}
  for m,c in records:q=Q(c);assert den%q.denominator==0;k=pack(m);assert k not in v;v[k]=q.numerator*(den//q.denominator)
  return v,den
 mat=[[parse(A3[i][j],D['A3_bounds'][i][j]) for j in range(3)] for i in range(3)];vec=[parse(V[i],D['v17_bounds'][i]) for i in range(3)];products=0;peak=0;maxbits=0;held=sum(2*len(p)+2 for row in mat for p,d in row)+sum(2*len(p)+2 for p,d in vec)
 def put(v,k,c):
  nonlocal maxbits
  if c:maxbits=max(maxbits,abs(c).bit_length());assert abs(c).bit_length()<=D['integer_bits_bound'];v[k]=c
  else:v.pop(k,None)
 results=[]
 for j in range(3):
  acc={};denout=D['output_bounds'][j]['denominator']
  for i in range(3):
   step=S['steps'][3*j+i];rr=step['minor_rows'];cc=step['minor_columns'];den=step['cofactor']['D'];co={}
   for (r0,c0,r1,c1,sign) in [(rr[0],cc[0],rr[1],cc[1],step['cofactor_sign']),(rr[0],cc[1],rr[1],cc[0],-step['cofactor_sign'])]:
    left,dl=mat[r0][c0];right,dr=mat[r1][c1];factor=den//(dl*dr);assert factor*dl*dr==den
    for ka,a in left.items():
     for kb,b in right.items():put(co,ka+kb,co.get(ka+kb,0)+sign*a*b*factor);products+=1
   v,dv=vec[i];factor=denout//(den*dv);assert factor*den*dv==denout
   for ka,a in co.items():
    for kb,b in v.items():put(acc,ka+kb,acc.get(ka+kb,0)+a*b*factor);products+=1
   peak=max(peak,held+2*len(co)+2+2*len(acc)+2+15000);assert peak<=g['whole_numeric_bound']<=CAP
  out=read(ROOT/'reports'/A/f'Cramer{j}.json');provided=out['integer_coefficients'];d2=int(out['denominator']);assert len(provided)==len(acc)
  peak=max(peak,held+2*len(acc)+2+4*len(provided)+2+15000);assert peak<=g['whole_numeric_bound']
  for m,c in provided:assert acc.get(pack(m),0)*d2==int(c)*denout
  assert len(acc)==[0,0,1768][j];results.append({'output':j,'all_coefficients_exact':len(acc),'denominator':str(denout),'provided_sha256':sha(ROOT/'reports'/A/f'Cramer{j}.json')})
 assert products<=D['integer_products_bound'] and maxbits<434
 out=ROOT/'reports'/REVIEW;out.mkdir(exist_ok=False);ep=out/'independent-three-Cramer-fullQQ-and-unit-transport-proof.json';dump(ep,{'all12_source_polynomials_and_rowunit_scales_verified':True,'all_three_output_fullQQ_matches':results,'Cramer_terms':[0,0,1768],
 'universal_adjugate_and_sign_independently_verified':True,'unit_and_fourth_row':'detA3=-kappa3, inheritedunit inverse detA8/(cprod*Saux*T); xfirst3=-Y/kappa3,x4=0. Independentkappa17zero andkappa-annihilator first3columns force fourth residual0.',
 'actual_pair_consequence':'Onlytarget17: first2andfourthpaircoordinateszero; coefficientofanchor22 is -Y2/kappa3, nonzero as a rational function. Thus pairv17 is proportional to actualanchor22 on the existing chart; no anchor inversion is introduced.',
 'all_ten_block_and_tail_corrections_retained':True,'no_whole_lastrow_global_rank_or57_239_claim':True,'producer_native_or_old_checkpoint_used':False,'actual_integer_products':products,'actual_integer_bits':maxbits,'peak_numeric_with_metadata':peak})
 report={'id':REVIEW,'all_passed':True,'Cramer_fullQQ_terms_verified':[0,0,1768],'same_target17_pair_proportional_anchor22_accepted':True,'no_whole_lastrow_or_global_rank_claim':True,
 'input_sha256':{**pins,str(gp.relative_to(ROOT)):sha(gp)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{REVIEW}.md'),'evidence_sha256':{ep.name:sha(ep)},'actual_integer_products':products,'actual_integer_bits':maxbits,'peak_numeric_with_metadata':peak,'seconds':time.monotonic()-t};dump(out/'report.json',report);print(json.dumps({'all_passed':True,'Cramer_terms':[0,0,1768],'peak':peak,'products':products,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['source','verify'],required=True);a=p.parse_args();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));source() if a.mode=='source' else verify()
