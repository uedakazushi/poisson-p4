#!/usr/bin/env python3
"""Independent fullC40 source and PROVIDED first/quad certificates only."""
import gc,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
from ambient_B18_honest_nil_T_gauge120_review_20261105_v63 import Base,Cover
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID;WORK=ROOT/'work'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def construct(D):
 B=Base([()]+[(i,)for i in range(9)]);mapping={0:0}
 for j,(z,r,e)in enumerate(D['table']['basis33']):
  if D['table']['basis_degrees'][j]==1:mapping[j]=1+(z-1 if z else 7 if r else 8)
 assert set(mapping.values())==set(range(10))
 roots={i:Q.pairs(v)for i,v in D['T']['all_requested164_root_values']}
 def project(v):return {mapping[i]:F(c)for i,c in enumerate(v)if i in mapping and F(c)}
 nu=[project([roots[f['actual_tracefree_nu']].get(i,F(0))for i in range(33)])for f in D['source']['families'].values()];ring=Cover(B,nu)
 kappa=[{(mapping[i],mask):F(c)for i in range(33)for mask in range(4)if i in mapping and(c:=F(v[4*i+mask]))}for v in D['kappa']['kappa9_C132']]
 core=[kappa[i]for i in[0,1,2,6]];a,d,x,t=core
 def poly(p,vs):
  out={}
  for term in p:
   e=term[0];c=F(term[1],term[2])if len(term)==3 else F(term[1]);v=ring.c(c)
   for z,n in zip(vs,e):
    if n:v=ring.mul(v,ring.power(z,n))
   out=ring.add(out,v)
  return out
 Dunit=ring.add(ring.c(1),ring.neg(ring.mul(d,d)));Punit=ring.add(ring.c(5),ring.mul(d,d));assert Dunit.get((0,0),F(0))==1 and Punit.get((0,0),F(0))==5
 units=[];C5=D['core5']
 for name in C5['selected5_pivot_order']:
  p=C5['five_actual_rational_roots'][name];root=ring.tm(poly(p['numerator'],core),ring.power(Dunit,-p['D_power']),ring.power(Punit,-p['P_power']));units.append(ring.add(kappa[C5['source9_variable_order'].index(name)],ring.neg(root)))
 h=ring.add(a,ring.neg(d),x,ring.scale(t,-2));c=ring.scale(ring.add(d,x),F(1,2));r=ring.scale(ring.add(ring.scale(a,-2),d,ring.neg(x)),F(1,6));l=ring.scale(ring.add(a,d,ring.neg(x)),F(1,3));vs=units+[h,c,r,l];R=D['rad'];dd=poly(R['d_polynomial'],vs);assert dd==d
 columns=[];originchecks=0
 def append(p):
  nonlocal originchecks
  zero=[p.get((0,1),F(0)),p.get((0,3),F(0)),p.get((0,2),F(0))-F(11,12)*p.get((0,1),F(0))];assert not any(zero);originchecks+=1;v={}
  for var in range(9):
   vals=[p.get((var+1,1),F(0)),p.get((var+1,3),F(0)),p.get((var+1,2),F(0))-F(11,12)*p.get((var+1,1),F(0))]
   for target,c in enumerate(vals):
    if c:v[9*target+var]=c
  columns.append(v)
 for p in R['R9_generators_polynomials9']:
  v=ring.mul(dd,poly(p,vs))
  for mask in range(4):append(ring.mul(v,{(0,mask):F(1)}))
 g=project([roots[7790].get(i,F(0))for i in range(33)]);assert not g
 for i in range(8):append({(i+1,1):F(-1)}if i<7 else{(j,1):-c for j,c in g.items()})
 assert len(columns)==44 and originchecks==44
 m2=list(itertools.combinations_with_replacement(range(9),2));index={m:i for i,m in enumerate(m2)};qcols=[]
 for v in columns:
  for var in range(9):
   z={}
   for i,c in v.items():t,k=divmod(i,9);Q.put(z,45*t+index[tuple(sorted((k,var)))],c)
   qcols.append(z)
 return columns,qcols,dict(all44_fullC40_origins_zero=True,all_nine_variables_and_monic_carries=True,root_rational_numerators_parsed_exactly_as_num_den_triples=True,proper_polynomials_not_linearized_beforetrace=True,actual_d_polynomial_equals_kappa1=True,actual_gA_first_zero=True,unit_inverses_checked=ring.inverse_count)
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);WORK.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63';rep=load(primary);assert pins[f'reports/{primary}/report.json']=='e2bbf49070e0f4d25b637a793d52bc438df6144533175617e84f0701a4c2e7d1'
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 history=load(primary,'source-history-lossless-parser-receipt.json');assert pins[f'reports/{primary}/source-history-lossless-parser-receipt.json']=='1ff90e249aa7743e40a595a81308ada5fe60d0e6c094ca5a0c30b73496726bde'
 sourcepeer=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_20261105_v63_REVIEW');oldhash=history['prepared_source_review_sha256'];snapshot=ROOT/history['exact_preserved_prepared_source'];assert sha(snapshot)==oldhash=='1be2cda944e1dfc0a608e0ac03035168438f4240abc279ad34739146c474b028'
 worker=ROOT/'code/ambient_B18_Cprime_linear_initial_values_20261105_v63.py';assert sha(worker)==history['executed_final_source_sha256']==rep['code_sha256']
 for p,h in sourcepeer['input_sha256'].items():
  if p==str(worker.relative_to(ROOT)):assert h==oldhash and sha(snapshot)==h
  else:assert sha(ROOT/p)==h
 p1=load(primary,'phase1-actual-fullnine-first-V-columns.json');p2=load(primary,'phase2-known-quadratic-closure.json')
 names={'table':('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63','exact-existing-raw33-multiplication-table.json'),'T':('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-source105-ell2-raw33-values.json'),'source':('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json'),'kappa':('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63','phase3-actual-kappa9-and-normalized-nilframe.json'),'core5':('AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63','actual4variable7core-polynomials.json'),'rad':('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json')};D={k:load(*p)for k,p in names.items()}
 assert load('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW')['all_passed'];save('phase.json',dict(phase='immutable actualC40source and exact prepared-code history pinned'))
 L,Q2,proof=construct(D);assert L==[Q.pairs(v)for v in p1['all44_actual_linear_columns27_sparse']]and Q2==[Q.pairs(v)for v in p2['all396_known_quadratic_columns135_sparse']]
 save('independent-fullC40-and-all-known-quadratic-source-columns.json',dict(**proof,all44_actual_linear_columns27=L,all396_known_quadratic_columns135=Q2,no_actual_quadraticV2_or_newsource=True))
 control=WORK/'two-fixed-rational-certificate-input.txt';expected=[]
 with control.open('w')as f:
  f.write('2\n')
  def col(v):
   f.write(str(len(v))+'\n')
   for i,c in sorted(v.items()):f.write(str(i)+' '+str(c)+'\n')
  for n,cols,cert,rank,hashkey in[(27,L,p1,24,'all44_factorization_sha256'),(135,Q2,p2,131,'all396_factorization_sha256')]:
   cp,rr=cert['fixed_pivot_columns'],cert['fixed_pivot_rows'];assert len(cp)==len(rr)==rank and len(set(cp))==len(set(rr))==rank;iv=cert['fixed_minor_inverse'];assert len(iv)==rank and all(len(row)==rank for row in iv)
   f.write(f'{n} {len(cols)} {rank}\n'+' '.join(map(str,rr))+'\n'+' '.join(map(str,cp))+'\n')
   for j in cp:col(cols[j])
   for row in iv:f.write(' '.join(str(F(c))for c in row)+'\n')
   for v in cols:col(v)
   expected.append(cert[hashkey])
 # Release source/certificate coefficient arrays before the native fixed
 # matrix products; peak native slots incl12000metadata is<=85555.
 del D,L,Q2,p1,p2,sourcepeer;gc.collect()
 cpp=ROOT/'code/ambient_B18_Cprime_linear_initial_values_review_20261105_v63.cpp';exe=WORK/'checker';native=OUT/'native-fixed-products-and-streamed-factorizations.json'
 comp=subprocess.run(['bash','-lc',f'g++ -O2 -std=c++17 {cpp} -o {exe} $(pkg-config --cflags --libs flint openssl)'],capture_output=True,text=True,timeout=15);(OUT/'compile.log').write_text(comp.stdout+comp.stderr);assert comp.returncode==0
 run=subprocess.run([str(exe),str(control),str(native)],capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-st)));(OUT/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0
 N=json.loads(native.read_text());assert N['all_passed']and N['all_streamed_columns']==440 and[N['blocks'][i]['factorization_sha256']for i in range(2)]==expected
 assert rep['actual_Cprime_HF1']==3 and rep['known_quadratic_quotient_dimension']==4 and not rep['positive_Cprime_ell3_length9_proved']
 assert all('phase3-known-cubic'not in n for n in rep['evidence_sha256'])
 save('independent-ranks-source-history-and-positive-only-scope.json',dict(actual_fullnine_V1_rank=24,actual_Cprime_HF1=3,known_quadratic_closure_rank=131,known_quadratic_quotient_upper4=True,not_actual_HF2=True,conditional_cubic_not_run=True,no_long_exponent_claim=True,provided_fixed_inverses_both_products=True,all440_column_span_checks_and_factor_hashes=True,source_history=history,source_history_scope='Only exacthistorical prepared-code pin is resolved throughits preserved1be2 snapshot; all other hashes are literal.',no_new_source_gauge_or_value_solve=True))
 for p in[Path(__file__),Path(Q.__file__),cpp,control,exe,snapshot,worker,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 raw=sum(p.stat().st_size for d in[OUT,WORK]for p in d.rglob('*')if p.is_file());assert raw<32<<20 and time.monotonic()-st<60
 final=dict(id=ID,all_passed=True,frozen=True,actual_fullnine_V1_rank=24,actual_Cprime_HF1=3,known_quadratic_closure_rank=131,known_quadratic_quotient_dimension4=True,cubic_and_long_torsion_exponent_not_decided=True,no_new_source_gauge_or_solve=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=raw,seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 save('report.json',final);print(json.dumps({k:v for k,v in final.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
