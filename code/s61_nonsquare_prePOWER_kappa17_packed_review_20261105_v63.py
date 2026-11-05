#!/usr/bin/env python3
"""Independent packed-exponent integer Laurent check; no producer imports."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_PREPOWER_KAPPA17_ENTRY_20261105_v63'
GATE='S61_NONSQUARE_PREPOWER_KAPPA17_PACKED_REVIEW_SOURCE_GATE_20261105_v63'
REVIEW=BASE+'_REVIEW'
PLAN='reports/S61_NONSQUARE_PREPOWER_KAPPA17_SOURCE_GATE_20261105_v63/literal-source-plan.json'
LITERAL='reports/S61_NONSQUARE_PREPOWER_KAPPA17_LITERAL_INPUT_GATE_V2_20261105_v63/literal-atoms20.json'
CAP=75410;RADIX=64
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,d):p.write_text(json.dumps(d,separators=(',',':'))+'\n')
def pack(m):return int(m[0])+64*int(m[1])+4096*int(m[2])
def lattice(box,degree):
 return sum(max(0,min(box[2][1],degree[1]-a-e)-max(box[2][0],degree[0]-a-e)+1) for a in range(box[0][0],box[0][1]+1) for e in range(box[1][0],box[1][1]+1))
def shape(records,den=Q(1),ashift=0):
 vals=[((int(m[0])-ashift,int(m[1]),int(m[2])),Q(c)/den) for m,c in records if Q(c)]
 assert len({m for m,c in vals})==len(vals)
 if not vals:return None
 D=math.lcm(*(q.denominator for m,q in vals));U=sum(abs(q.numerator)*(D//q.denominator) for m,q in vals)
 return {'support':len(vals),'box':[[min(m[i] for m,q in vals),max(m[i] for m,q in vals)] for i in range(3)],'degree':[min(sum(m) for m,q in vals),max(sum(m) for m,q in vals)],'denominator':D,'integer_l1':U}
def combine(a,b,multiply=False):
 if multiply and (a is None or b is None):return None
 if a is None:return b
 if b is None:return a
 if multiply:
  box=[[u[0]+v[0],u[1]+v[1]] for u,v in zip(a['box'],b['box'])];degree=[a['degree'][0]+b['degree'][0],a['degree'][1]+b['degree'][1]];D=a['denominator']*b['denominator'];U=a['integer_l1']*b['integer_l1'];support=a['support']*b['support']
 else:
  box=[[min(u[0],v[0]),max(u[1],v[1])] for u,v in zip(a['box'],b['box'])];degree=[min(a['degree'][0],b['degree'][0]),max(a['degree'][1],b['degree'][1])];D=math.lcm(a['denominator'],b['denominator']);U=a['integer_l1']*(D//a['denominator'])+b['integer_l1']*(D//b['denominator']);support=a['support']+b['support']
 return {'support':min(support,lattice(box,degree)),'box':box,'degree':degree,'denominator':D,'integer_l1':U}
def inputs():
 paths=[f'reports/{BASE}/report.json',PLAN,LITERAL,'reports/S61_NONSQUARE_PREPOWER_KAPPA17_SOURCE_GATE_20261105_v63/report.json','reports/S61_NONSQUARE_PREPOWER_KAPPA17_LITERAL_INPUT_GATE_V2_20261105_v63/report.json','reports/S61_NONSQUARE_PREPOWER_KAPPA17_DENSE_OWNERSHIP_GATE_20261105_v63/report.json']
 P=read(ROOT/paths[0]);assert P['all_passed'] and P['kappa17_is_zero']
 # Hash only the explicit small frozen custody/worker inputs. In particular,
 # never open the old25x32 checkpoint mentioned inside source references.
 for p,h in P['input_sha256'].items():assert sha(ROOT/p)==h
 plan=read(ROOT/PLAN);literal=read(ROOT/LITERAL)
 assert plan['variable_order']==['a','e','d'] and plan['actual_current15_pair_rows']==[2,4,5,8] and plan['actual_original32_target']==17
 assert literal['source49_definitions']==[{'node':r['source_node'],'definition':r['definition']} for r in plan['source_nodes']]
 assert literal['source_order']==plan['source_root_order']==[2611,2821,2926,3241]
 assert literal['kappa_sources']==plan['kappa_coefficients_in_exact_order']
 refs={r['node']:r for r in plan['source_atom_references']}
 assert len(refs)==len(literal['atom_payloads'])==20
 for r in literal['atom_payloads']:
  assert r['source_reference']==refs[r['node']]['source'] and r['source_sha256']==refs[r['node']]['source_reference_sha256']
 paths.extend(r['path'] for r in plan['kappa_coefficients_in_exact_order'])
 paths.extend(f'reports/{BASE}/source-root{i}.json' for i in range(4))
 for p,h in P['evidence_sha256'].items():
  if p.startswith('source-root'):assert sha(ROOT/'reports'/BASE/p)==h
 return P,plan,literal,{p:sha(ROOT/p) for p in dict.fromkeys(paths)}
def source_gate():
 t=time.monotonic();P,plan,literal,pins=inputs();atoms={}
 for rec in literal['atom_payloads']:
  numerator,denominator=rec['original_num_den'];assert len(denominator)==1
  m,c=denominator[0];assert m[1:]==[0,0] and Q(c)
  atoms[rec['node']]=shape(numerator,Q(c),m[0])
 bounds={};ops=0
 for r in plan['source_nodes']:
  i=r['source_node'];node=r['definition'];op=node[0]
  if op=='z':b=None
  elif op in ('a','s','m'):
   b=None if op=='s' and node[1]==node[2] else combine(bounds[node[1]],bounds[node[2]],op=='m')
   if op=='m' and bounds[node[1]] and bounds[node[2]]:ops+=bounds[node[1]]['support']*bounds[node[2]]['support']
  else:b=atoms[i]
  assert b==r['bound'];bounds[i]=b
 kappas=[]
 for ref,n in zip(plan['kappa_coefficients_in_exact_order'],[1257,1257,1232,1174],strict=True):
  p=ROOT/ref['path'];assert sha(p)==ref['sha256'];records=read(p)['polynomial'];b=shape(records);assert b['support']==n;kappas.append(b)
 accumulator=None;Dstar=1;allbounds=[b for b in bounds.values() if b]+kappas;products=[]
 for j,root in enumerate(plan['source_root_order']):
  product=combine(kappas[j],bounds[root],True);accumulator=combine(accumulator,product);assert product==plan['dot_product_schedule'][j]['product'] and accumulator==plan['dot_product_schedule'][j]['updated_accumulator'];products.append(product);Dstar=math.lcm(Dstar,product['denominator']);ops+=kappas[j]['support']*bounds[root]['support'];allbounds.extend([product,accumulator])
 amin=min(b['box'][0][0] for b in allbounds);amax=max(b['box'][0][1] for b in allbounds);emax=max(b['box'][1][1] for b in allbounds);dmax=max(b['box'][2][1] for b in allbounds)
 assert amax-amin<64 and emax<64 and all(b['box'][1][0]>=0 and b['box'][2][0]>=0 for b in allbounds)
 bits=max(max(b['integer_l1'].bit_length(),b['denominator'].bit_length()) for b in allbounds);assert bits==201<548
 # Own representation charges one packed key and one integer coefficient;
 # direct multiplication adds into the old accumulator, never a product copy.
 payload=4*sum(len(r['original_num_den'][0])+len(r['original_num_den'][1]) for r in literal['atom_payloads'])
 atomcharge=sum(2*b['support']+2 for b in atoms.values());metadata=15000
 users={r['source_node']:0 for r in plan['source_nodes']}
 for r in plan['source_nodes']:
  if r['definition'][0] in ('a','s','m'):
   for c in r['definition'][1:]:users[c]+=1
 for root in plan['source_root_order']:users[root]+=1
 cache=dict(atoms);sourcepeak=payload+atomcharge+metadata
 def charge(c):return sum(2*b['support']+2 for b in c.values() if b)
 for r in plan['source_nodes']:
  i=r['source_node'];op=r['definition'][0];b=bounds[i]
  if op in ('a','s','m','z'):
   sourcepeak=max(sourcepeak,charge(cache)+(0 if b is None else 2*b['support']+2)+metadata);cache[i]=b
  if op in ('a','s','m'):
   for c in r['definition'][1:]:
    users[c]-=1
    if not users[c]:cache.pop(c,None)
 rootcharge=sum(2*bounds[i]['support']+2 for i in plan['source_root_order']);onejson=max(4*b['support']+2 for b in kappas);onesparse=max(2*b['support']+2 for b in kappas);acccharge=2*accumulator['support']+2
 dotpeak=rootcharge+onejson+onesparse+acccharge+metadata
 rootcomparepeak=rootcharge+4*max(bounds[i]['support'] for i in plan['source_root_order'])+2+metadata
 gatepeak=payload+onejson+metadata+2000;whole=max(sourcepeak,dotpeak,rootcomparepeak,gatepeak)
 assert whole<=CAP<150000 and ops<4_000_000
 external=sum((ROOT/p).stat().st_size for p in pins);raw=external+500000;assert raw<32*1024*1024
 evidence={'representation':'one integer numerator dictionary keyed by a+64e+4096d, with fixed common denominator per node; direct four-term accumulator, no separate product dictionary',
  'all49_source_bounds_and20_literal_custody_verified':True,'old_checkpoint_opened':False,'global_exponent_bounds':[amin,amax,emax,dmax],'radix_injective_proof':'difference in a has absolute value<64; then e difference<64, so packed equality forces all three exponents equal. Packed addition represents multiplication even for negative a.',
  'source49_integer_l1_bounds':{str(i):None if b is None else b for i,b in bounds.items()},'kappa4_bounds':kappas,
  'common_accumulator_denominator':str(Dstar),'integer_partialsum_bound_bits':201,'original_integer_cap_548_not_raised':True,
  'numeric_ownership_bounds':{'raw_literal_payload':payload,'all_parsed_atoms':atomcharge,'source_or_parse_phase':sourcepeak,'direct_dot_phase':dotpeak,'one_root_comparison':rootcomparepeak,'source_gate_phase':gatepeak,'all_metadata_and_temporary_reserve':metadata,'whole':whole},
  'integer_product_bound':ops,'external_RAW':external,'TOTAL_RAW_upper':raw,'signed_source_or_product_values_evaluated':False}
 out=ROOT/'reports'/GATE;out.mkdir(exist_ok=False);ep=out/'independent-packed-Laurent-source-and-resource-gate.json';dump(ep,evidence)
 report={'id':GATE,'all_passed':True,'source_only':True,'whole_numeric_ownership_bound':whole,'authorized_old_cap':CAP,'integer_bits_bound':201,'TOTAL_RAW_upper':raw,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{GATE}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(out/'report.json',report);print(json.dumps({'source_gate_passed':True,'owned_bound':whole,'bits':201,'products_bound':ops,'report_sha256':sha(out/'report.json')}))
def verify():
 t=time.monotonic();gp=ROOT/'reports'/GATE/'report.json';g=read(gp);assert g['all_passed'] and g['code_sha256']==sha(Path(__file__))
 for p,h in g['input_sha256'].items():assert sha(ROOT/p)==h
 D=read(gp.parent/'independent-packed-Laurent-source-and-resource-gate.json');P,plan,literal,pins=inputs();bounds={int(i):b for i,b in D['source49_integer_l1_bounds'].items()};cache={};peak=0;products=0;maxbits=0
 def live(extra=0):
  nonlocal peak
  q=sum(2*len(v)+2 for v,d in cache.values())+extra+15000;peak=max(peak,q);assert q<=g['whole_numeric_ownership_bound']<=CAP
 def put(d,k,c):
  nonlocal maxbits
  if c:
   maxbits=max(maxbits,abs(c).bit_length());assert abs(c).bit_length()<=201
   d[k]=c
  else:d.pop(k,None)
 for rec in literal['atom_payloads']:
  i=rec['node'];num,den=rec['original_num_den'];m,c=den[0];denom=Q(c);bound=bounds[i];targetD=bound['denominator'];poly={}
  for exp,cs in num:
   value=Q(cs)/denom;assert targetD%value.denominator==0;key=pack(exp)-m[0];assert key not in poly;put(poly,key,value.numerator*(targetD//value.denominator))
  cache[i]=(poly,targetD);live(D['numeric_ownership_bounds']['raw_literal_payload'])
 literal.clear();rec=None;num=None;den=None
 nodes=plan['source_nodes'];roots=plan['source_root_order'];users={r['source_node']:0 for r in nodes}
 for r in nodes:
  if r['definition'][0] in ('a','s','m'):
   for child in r['definition'][1:]:users[child]+=1
 for root in roots:users[root]+=1
 observed=[]
 for rec in nodes:
  i=rec['source_node'];op=rec['definition'][0];bound=bounds[i]
  if op=='z':cache[i]=({},1)
  elif op in ('a','s','m'):
   left,right=[cache[j] for j in rec['definition'][1:]];out={};targetD=bound['denominator'] if bound else 1
   if op=='m':
    scaleD=targetD//(left[1]*right[1]);assert targetD==scaleD*left[1]*right[1]
    for ka,a in left[0].items():
     for kb,b in right[0].items():put(out,ka+kb,out.get(ka+kb,0)+a*b*scaleD);products+=1
   else:
    for poly,den,sign in [(left[0],left[1],1),(right[0],right[1],-1 if op=='s' else 1)]:
     assert targetD%den==0
     for k,a in poly.items():put(out,k,out.get(k,0)+sign*a*(targetD//den))
   live(2*len(out)+2);cache[i]=(out,targetD)
  poly,den=cache[i]
  assert bound is None and not poly or bound is not None and len(poly)<=bound['support'] and sum(map(abs,poly.values()))<=bound['integer_l1'] and den==bound['denominator']
  observed.append([i,len(poly),den])
  if op in ('a','s','m'):
   for child in rec['definition'][1:]:
    users[child]-=1
    if not users[child]:del cache[child]
 assert set(cache)==set(roots)
 rootmatches=[]
 for n,root in enumerate(roots):
  records=read(ROOT/'reports'/BASE/f'source-root{n}.json');poly,den=cache[root];den2=int(records['denominator']);terms=records['integer_coefficients'];assert len(terms)==len(poly);live(4*len(terms)+2)
  for exp,value in terms:assert poly.get(pack(exp),0)*den2==int(value)*den
  rootmatches.append({'source_root':root,'all_terms_exact':len(poly),'provided_hash':sha(ROOT/'reports'/BASE/f'source-root{n}.json')});records=None;terms=None
 total={};Dstar=int(D['common_accumulator_denominator'])
 for n,ref in enumerate(plan['kappa_coefficients_in_exact_order']):
  records=read(ROOT/ref['path'])['polynomial'];b=D['kappa4_bounds'][n];den=b['denominator'];kappa={}
  for exp,value in records:
   q=Q(value);assert den%q.denominator==0;k=pack(exp);assert k not in kappa;put(kappa,k,q.numerator*(den//q.denominator))
  root,rd=cache[roots[n]];factor=Dstar//(rd*den);assert factor*rd*den==Dstar
  live(4*len(records)+2+2*len(kappa)+2+2*D['source49_integer_l1_bounds'][str(roots[n])]['support']+2)
  for ka,a in root.items():
   for kb,bv in kappa.items():put(total,ka+kb,total.get(ka+kb,0)+a*bv*factor);products+=1
  live(4*len(records)+2+2*len(kappa)+2+2*len(total)+2)
  assert products<=D['integer_product_bound'];del cache[roots[n]];records=None;kappa=None
 assert not total and maxbits<=201 and peak<=CAP
 out=ROOT/'reports'/REVIEW;out.mkdir(exist_ok=False);ep=out/'independent49-source-four-root-and-kappa17-zero.json'
 dump(ep,{'all49_nodes_reconstructed_from_literal20':True,'source_node_observed_termcounts_and_denominators':observed,'all4_source_root_QQ_equalities':rootmatches,
          'all4_kappa_fixed_order_no_extra_sign':True,'kappa17_exact_zero':True,'result_scope':'Only x4=kappa17/T=0 for original32 target17. No wholelastrow/globalrank/full57/239 inference; other3pair and10block corrections retained.',
          'old_checkpoint_read':False,'producer_import_or_native_reuse':False,'actual_integer_products':products,'actual_integer_bits':maxbits,'peak_owned_numeric_with_metadata':peak})
 report={'id':REVIEW,'all_passed':True,'kappa17_exact_zero_independently_verified':True,'whole_lastrow_or_global_rank_claimed':False,
  'input_sha256':{**pins,str(gp.relative_to(ROOT)):sha(gp)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{REVIEW}.md'),'evidence_sha256':{ep.name:sha(ep)},
  'actual_integer_products':products,'actual_integer_bits':maxbits,'peak_owned_numeric_with_metadata':peak,'seconds':time.monotonic()-t};dump(out/'report.json',report);print(json.dumps({'all_passed':True,'kappa17_zero':True,'products':products,'peak':peak,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['source','verify'],required=True);a=p.parse_args();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));source_gate() if a.mode=='source' else verify()
