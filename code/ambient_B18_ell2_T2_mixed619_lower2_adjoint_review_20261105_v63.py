#!/usr/bin/env python3
"""Independent fixed-certificate replay; no producer imports or rank search."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as CACHE
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(v):return {int(i):Q(x) for i,x in v if Q(x)}
def record(v):return [[i,str(x)] for i,x in sorted(v.items())]
def count(v):
 if isinstance(v,(int,Q)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(tuple,list)):return sum(map(count,v))
 return 0
def main(out):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={};peakbits=0
 def load(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=sha(ROOT/p);return CACHE.selected(p,k)
 def guard(x):
  nonlocal peakbits
  peakbits=max(peakbits,abs(x.numerator).bit_length(),x.denominator.bit_length());assert peakbits<=384;return x
 def plus(v,w,c=Q(1)):
  q=dict(v)
  for i,x in w.items():
   y=guard(q.get(i,Q(0))+guard(c*x))
   if y:q[i]=y
   else:q.pop(i,None)
  return q
 def image(b,v):
  q={}
  for i,x in v.items():q=plus(q,b[i],x)
  return q
 def dot(a,b):
  q=Q(0)
  for i,x in a.items():
   if i in b:q=guard(q+guard(x*b[i]))
  return q
 base='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63'
 primary=load(base+'/report.json');assert hs[base+'/report.json']=='0fbbef09a102e9bfad24eea83665057b1a462c7db68bed1d2d000e6ac1e451ff' and primary['all_passed']
 evpath=base+'/actual-two-mixed-output-adjoint-certificates-or-obstruction.json';cert=load(evpath);assert hs[evpath]==primary['evidence_sha256'][Path(evpath).name]
 custody=load('reports/AMBIENT_B18_ELL2_T2_MIXED619_ADJOINT_PARENT_CUSTODY_20261105_v63/report.json');assert custody['all_passed'] and custody['same_exact_bytes']
 for p,h in custody['input_sha256'].items():assert sha(ROOT/p)==h
 for p,r in primary['input_sha256'].items():assert sha(ROOT/r['path'])==r['sha256']
 inp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_MIXED_ADJOINT_INPUT_20261105_v63/actual136-mixed-output210-slot-adjoint-source-input.json';a=load(inp)
 covpath='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_COVECTORS_20261105_v63/actual-four-full918-source-covectors-and52-plane-adjoint-input.json';savedcov=load(covpath)
 pre='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';G3=list(map(sp,field(pre,'actual_full_G3linear102_by418_columns')));PY=list(map(sp,field(pre,'PiY3_from102_columns')))
 rank25=load('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json')
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';B=list(map(sp,field(mu,'actual_selected_Y0_columns102')));bn={i:sp(v) for i,v in field(mu,'normalized_basis_rows')};br={i:sp(v) for i,v in field(mu,'normalized_basis_original_selected_source_representations')}
 hsfirst=load('reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63/exact171-decomposition-162-zero-9-lower-first-source-input.json');Hr=list(map(sp,hsfirst['Hstar1_r_selected99_columns']));He=list(map(sp,hsfirst['Hstar1_ell_selected99_columns']))
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';G2=list(map(sp,field(mug,'full_G2_linear234_by102_columns')))
 input_owned=count([primary,cert,custody,a,savedcov,G3,PY,rank25,B,bn,br,hsfirst,G2]);assert input_owned+30000<150000
 Hphase={'phase':'independent_fixed_receipts_and_original_source_only','owned_input_plus_scratch':input_owned+30000,'no_rank_search_or_unknown_values':True};(out/'phase.json').write_text(json.dumps(Hphase)+'\n')
 # Rebuild all25 selected full/source projections, then replay ONLY the saved
 # pivot sequence. No basis or pivot is searched.
 chosen=rank25['selected_columns'];J=[G3[j] for j in chosen];JP=[]
 for col in J:
  p={}
  for k,x in col.items():
   row,var=divmod(k,9)
   for s,y in PY[row].items():p=plus(p,{9*s+var:y},x)
  JP.append(p)
 norms={};reps={}
 for n,r in enumerate(rank25['fixed_rank_receipts']):
  assert r['column']==chosen[n] and r['basis_column']==n;v=dict(JP[n]);w={n:Q(1)}
  for p,cs in r['subtract']:
   c=Q(cs);assert v.get(p,Q(0))==c;v=plus(v,norms[p],-c);w=plus(w,reps[p],-c)
  p=r['pivot'];d=Q(r['divide_by']);assert d and v.get(p)==d and min(v)==p;norms[p]={i:guard(x/d) for i,x in v.items()};reps[p]={i:guard(x/d) for i,x in w.items()};assert image(JP,reps[p])==norms[p]
 assert norms=={i:sp(v) for i,v in rank25['normalized_basis']}
 def extract(q,norm,rep):
  v=dict(q);w={}
  for p in sorted(norm):
   c=v.get(p,Q(0))
   if c:v=plus(v,norm[p],-c);w=plus(w,rep[p],c)
  return w,v
 F25=[extract({i:Q(1)},norms,reps)[0] for i in range(27)]
 # Recheck the old99 point-unit rows, then pull the two literal Hstar duals
 # through their SAME fixed triangular extraction. This reconstructs alpha
 # from original point/first arrays, rather than trusting supplied targets.
 for p in bn:assert min(bn[p])==p and image(B,br[p])==bn[p]
 duals=[[Q(0),Q(1),Q(4)],[Q(11),Q(0),Q(270)]];sr=[];se=[]
 for ell in duals:
  hr=[sum((ell[k]*v.get(k,Q(0)) for k in range(3)),Q(0)) for v in Hr];he=[sum((ell[k]*v.get(k,Q(0)) for k in range(3)),Q(0)) for v in He];rrow={};erow={}
  for i in range(102):
   c,_=extract({i:Q(1)},bn,br);r=dot({j:x for j,x in enumerate(hr) if x},c);e=dot({j:x for j,x in enumerate(he) if x},c)
   if r:rrow[i]=r
   if e:erow[i]=e
  sr.append(rrow);se.append(erow)
 alpha=[]
 for r,e in zip(sr,se):
  raw=[{9*i+7:x for i,x in r.items()},plus({9*i+7:x for i,x in e.items()},{9*i+8:x for i,x in r.items()})]
  for L in raw:
   LJ={j:dot(L,col) for j,col in enumerate(J)};v={}
   for i in range(102):
    for var in range(9):
     py={9*j+var:x for j,x in PY[i].items()};f=image(F25,py);q=guard(L.get(9*i+var,Q(0))-dot(LJ,f))
     if q:v[9*i+var]=q
   assert all(dot(v,col)==0 for col in J);alpha.append(v)
 assert alpha==list(map(sp,savedcov['full918_covectors4']));assert not alpha[0] and not alpha[1]
 # Rebuild all136 output covectors directly from the ORIGINAL G2linear,
 # using ordinary monomial coefficients: a square contributes once.
 pairs=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)];assert a['monomial_output_order8']==[list(v) for v in pairs] and a['active_target_order17']==list(range(17));outputs=[{} for _ in range(136)]
 for j,col in enumerate(G2):
  for k,c in col.items():
   target,var=divmod(k,9)
   if target>=17:continue
   for m,(u,v) in enumerate(pairs):
    if var==u:outputs[8*target+m]=plus(outputs[8*target+m],{9*j+v:c})
    if u!=v and var==v:outputs[8*target+m]=plus(outputs[8*target+m],{9*j+u:c})
 nonzero=[i for i,v in enumerate(outputs) if v];assert len(nonzero)==91 and nonzero==a['nonzero_output_indices']
 # The declared210-slot rectangle includes zero padding. It need not equal
 # the union of nonzero supports; every one of the other708 slots must vanish.
 domain=[9*j+v for j in range(42) for v in (0,4,5,7,8)];assert len(domain)==210 and domain==a['full918_domain_slots210'];assert all(set(v)<=set(domain) for v in outputs+alpha);pos={i:j for j,i in enumerate(domain)}
 assert [{pos[k]:x for k,x in outputs[i].items()} for i in nonzero]==list(map(sp,a['matrix_columns_on210_domain']))
 assert [{pos[k]:x for k,x in row.items()} for row in alpha[2:]]==list(map(sp,a['target_covectors2_on210_domain']))
 R91=list(map(sp,cert['solutions_on91_nonzero_output_columns']));R136=list(map(sp,cert['solutions_on136_literal_outputs']));res=[]
 for row,(r,rr,target) in enumerate(zip(R91,R136,alpha[2:])):
  assert rr=={nonzero[j]:x for j,x in r.items()} and len(r)==[21,42][row];assert all(i%8>=2 for i in rr)
  value=image(outputs,rr);assert value==target;assert all(value.get(i,Q(0))-target.get(i,Q(0))==0 for i in range(918));res.append({i:value.get(i,Q(0))-target.get(i,Q(0)) for i in range(918) if value.get(i,Q(0))!=target.get(i,Q(0))})
 owned=input_owned+count([J,JP,norms,reps,F25,sr,se,alpha,outputs,R91,R136])+10000;assert owned<150000
 ep=out/'independent-original-G2-alpha-and-all918-adjoint-checks.json';ep.write_text(json.dumps({'original_output_frame17_by8':[list(v) for v in pairs],'all136_output_rows_rebuilt':True,'nonzero_output_indices91':nonzero,'literal_domain_slots210':domain,'all_other708_slots_zero':True,'own_original_alpha4':[record(v) for v in alpha],'fixed_rho136':[record(v) for v in R136],'all2_times918_residuals':[record(v) for v in res],
  'fixed25_source_receipts_replayed_no_rank_search':True,'point99_and_Hstar_rebuilt_targets':True,'chiA_zero_scope':'Only its two cofactor alpha rows arezero; fullchiA lower relation retainsdirectchiA*Pi3 Z2 andisNOTevaluated.',
  'all_four_full_relation_source_formula':'Gamma=chi*Pi3_plane-rho*M0 onZ2=Xgamma2-Xp1*a1; chiA rho=0 stillhasdirectPi3 term.',
  'no_physical_gauge_X2_or619_relation_values':True},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'fixed_adjoint_nonzeros':[21,42],'all_literal_domain_equations_verified':1836,'original136_and_pruned91_and_domain210_verified':True,'chiA_direct_term_not_zero_filled':True,'no_rank_search_or_new_source_values':True,'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'owned_numeric_bound':owned,'actual_maximum_normalized_bits':peakbits,'preserved_first_reader_prefix':'reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER2_ADJOINT_20261105_v63_REVIEW_DOMAIN_READER_PREFIX','first_reader_issue':'overstrictly required fixed210 rectangle to equal nonzero support; zero padding is allowed and all708 outside slots are checked','seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','fixed_adjoint_nonzeros','all_literal_domain_equations_verified','owned_numeric_bound','actual_maximum_normalized_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
