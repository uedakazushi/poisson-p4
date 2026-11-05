#!/usr/bin/env python3
"""Independent scalar-symbol polynomial proof of FULL tauAB quartic sources."""
import hashlib,json,math,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
SOURCE='2debebb789b978069fc3dd5863cb5ff5f834366add91716c34065209f3dd852c'
MAXPRE=MAXCOF=0
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def guard(n):
 global MAXPRE
 MAXPRE=max(MAXPRE,n);assert n<=2048
def note(q):
 global MAXCOF
 MAXCOF=max(MAXCOF,abs(q.numerator).bit_length(),q.denominator.bit_length());assert MAXCOF<=2048
 return q
def qa(a,b):
 if not a:return b
 if not b:return a
 guard(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def qm(a,b):
 if not a or not b:return Q(0)
 guard(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def const(q):q=note(Q(q));return {():q} if q else {}
def plus(a,b):
 d=dict(a)
 for k,v in b.items():
  d[k]=qa(d.get(k,Q(0)),v)
  if not d[k]:del d[k]
 return d
def scale(a,c):return {k:qm(v,c) for k,v in a.items() if v and c}
def times(a,b):
 d={}
 for x,u in a.items():
  for y,v in b.items():
   k=tuple(sorted(x+y));assert len(k)<=2
   d[k]=qa(d.get(k,Q(0)),qm(u,v))
 return {k:v for k,v in d.items() if v}
def rpplus(a,b):
 d={m:dict(v) for m,v in a.items()}
 for m,v in b.items():d[m]=plus(d.get(m,{}),v)
 return {m:v for m,v in d.items() if v}
def rpscale(a,c):return {m:scale(v,c) for m,v in a.items() if v and c}
def rptimes(a,b):
 d={}
 for x,u in a.items():
  for y,v in b.items():
   m=tuple(i+j for i,j in zip(x,y));assert sum(m)<=4
   d[m]=plus(d.get(m,{}),times(u,v))
 return {m:v for m,v in d.items() if v}
def rpsum(seq):
 d={}
 for p in seq:d=rpplus(d,p)
 return d
def encoded(p):return [[list(k),str(v)] for k,v in sorted(p.items())]
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;assert out.exists() and not (out/'independent-full-tauAB-scalar-symbol-polynomial-proof.json').exists();pins={}
 def pin(p,h=None):
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  pins[str(p.relative_to(ROOT))]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pre=pin(out/'pre-expansion-support-and-whole-bound.json');assert pre['whole_symbolic_proof_upper']==113222<150000
 pr=pin(ROOT/'reports'/PRI/'report.json',SOURCE);docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 data=pin(ROOT/'reports'/PRI/'full-tauAB12-quartic-scalar-circuit-and-source-proof.json',
  'bd4ae81f549e2463fc6a972bca845b653719a69cff3c119ab5dd57b2279dfb1c')
 assert data['source_bound_metadata']['source_owned_upper']==36870
 nodes=data['nodes'];provided=data['provided_inputs'];assert len(nodes)==954 and len(provided)==39
 symbols={tuple(rec['key']):j for j,rec in enumerate(provided)}
 var=lambda key:{(symbols[tuple(key)],):Q(1)}
 frame=docs['actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json']
 c2=docs['actual-physical12-and-C2-right-exact-quotient.json'];raw=docs['actual-ell2-product33-and-raw45-actions.json'];basis=raw['actual33_basis']
 C={r['physical_row']:r for r in c2['physical12_coefficients']};rows=list(C);piv=[18,60,67];active=[0,4,5]
 first=[basis.index([i+1,0,0]) for i in range(7)];ir2,irl=[basis.index([0,*m]) for m in [(2,0),(1,1)]]
 ir3,ir2l=[basis.index([0,*m]) for m in [(3,0),(2,1)]]
 izr=[basis.index([i+1,1,0]) for i in range(7)];izl=[basis.index([i+1,0,1]) for i in range(7)]
 oldnu={i:dict(v) for i,v in docs['actual-source105-ell2-raw33-values.json']['all_requested164_root_values']}
 ell=docs['all164-cached-pure-ell-readbacks-and-physical12-second-values.json']
 lam={r['physical_row']:Q(r['delta_r_ell2']) for r in docs['phase1-six-root-existing-cache-readbacks.json']['all3_lambda_readbacks']}
 mu={r['physical_row']:Q(r['delta_ell3']) for r in docs['phase1-six-cached-pivot-pure-ell4-values.json']['all3_pivot_mu03_readbacks']}
 for rec in provided:
  k=rec['key'];kind=k[0]
  if kind=='delta':q=Q(dict(C[k[1]]['delta']).get(k[2],'0'))
  elif kind=='nu_sparse33':q=Q(oldnu[k[1]].get(k[2],'0'))
  elif kind=='nu_ell2':q=Q(ell['actual_nuA_B_pure_ell3'][k[1]][2])
  elif kind=='lambda':q=lam[k[1]]
  elif kind=='mu03':q=mu[k[1]]
  else:raise AssertionError(kind)
  assert q==Q(rec['value']) and tuple(nodes[rec['node']][1])==tuple(k)
 # Derive beta2 and quadratic normalizations from the already accepted
 # physical data; every requested response below remains a polynomial in39
 # INDEPENDENT scalar symbols, never their actual values.
 betaQ=[[Q(q) for q in r] for r in docs['four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json']['beta2_plane_r2_rell_ell2']]
 mon2=[(2,0),(1,1),(0,2)];mon3=[(3,0),(2,1),(1,2),(0,3)];mon4=[(4,0),(3,1),(2,2),(1,3),(0,4)]
 beta=[{m:const(q) for m,q in zip(mon2,row) if q} for row in betaQ]
 assert [i for i,b in enumerate(beta) if b]==active
 fm={r['physical_row']:r for r in frame['all12_verified_source_records']}
 r={(1,0):const(1)};l={(0,1):const(1)}
 def dd(p,ix):return var(('delta',p,ix))
 nu=[]
 for side,rootid in enumerate([20086,40865]):
  assert Q(oldnu[rootid].get(0,'0'))==Q(oldnu[rootid].get(basis.index([0,1,0]),'0'))==Q(ell['actual_nuA_B_pure_ell3'][side][1])==0
  nu.append({mon2[0]:var(('nu_sparse33',rootid,ir2)),mon2[1]:var(('nu_sparse33',rootid,irl)),mon2[2]:var(('nu_ell2',side))})
 nuAB=rpplus(nu[0],rpscale(rpsum([{m:times(v,var(('nu_sparse33',20086,first[i]))) for m,v in beta[i].items()} for i in active]),-Q(1)))
 delta=[];E=[]
 for p in piv:
  dp={mon3[0]:dd(p,ir3),mon3[1]:dd(p,ir2l),mon3[2]:var(('lambda',p)),mon3[3]:var(('mu03',p))};delta.append(dp)
  correction=rpsum([rptimes({(1,0):dd(p,izr[i]),(0,1):dd(p,izl[i])},beta[i]) for i in active])
  E.append(rpplus(dp,rpscale(correction,-Q(1))))
 oldHB=docs['two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json']
 T=[[Q(q) for q in row] for row in oldHB['constant_source_normalization_T']]
 H=[[Q(q) for q in row] for row in oldHB['source_first7_h_rows']]
 PA=[];PB=[]
 for k in range(3):
  PA.append(rpsum([rpscale(delta[j],T[j][k]) for j in range(3)]))
  PB.append(rpsum([rpscale(E[q],sum(T[j][k]*H[j][i] for j in range(3))) for q,i in enumerate(active)]))
 nf=[[Q(q) for q in row] for row in docs['full-HB6-by2-affine-source-operator-and-old27-leaves.json']['full_NF4']]
 knownv=[Q(row[0]) for row in docs['phase2-full-six-coordinate-HB-matrix-and-fixed-rank-certificate.json']['actual_full_HB6_by2_matrix']]
 assert knownv[0]==-Q(3,5)
 slopes=[q/knownv[0] for q in knownv[1:]]
 assert list(map(str,slopes))==data['five_residual_dual_slopes_vs_coordinate0']
 expected_polys=[];expected_nf=[];expected_res=[]
 for p in piv:
  pos=rows.index(p);bb={m:const(q) for m,q in zip(mon2,map(Q,fm[p]['bB_second_plane'])) if q}
  u=rpsum([rpscale(beta[i],Q(fm[p]['bA_first_z7'][i])) for i in range(7)])
  A=rpscale(rptimes(nu[1],bb),-Q(1));B=rpscale(rptimes(nuAB,u),-Q(1))
  B=rpplus(B,rpscale(rpsum([rptimes({(1,0):dd(p,izr[i]),(0,1):dd(p,izl[i])},E[q]) for q,i in enumerate(active)]),-Q(1)))
  clear=[{(1,0):dd(p,ir3),(0,1):dd(p,ir2l)},{(0,1):var(('lambda',p))},{(0,1):var(('mu03',p))}]
  A=rpplus(A,rpscale(rpsum([rptimes(clear[k],PA[k]) for k in range(3)]),-Q(1)))
  B=rpplus(B,rpscale(rpsum([rptimes(clear[k],PB[k]) for k in range(3)]),-Q(1)))
  assert all(sum(m)==4 for m in list(A)+list(B))
  full=[]
  for side,P in enumerate([A,B]):
   for row in nf[3*side:3*side+3]:
    v={}
    for m,c in zip(mon4,row[5*side:5*side+5]):v=plus(v,scale(P.get(m,{}),c))
    full.append(v)
  res=[plus(full[i],scale(full[0],-slopes[i-1])) for i in range(1,6)]
  expected_polys.append((A,B));expected_nf.append(full);expected_res.append(res)
 # Independent source polynomial interpretation with last-use release;
 # expand variables, NOT their cached numerical values.
 roots=pre['retained_root_ids'];uses=Counter()
 for node in nodes:
  if node[0]in['a','m']:uses.update(node[1:])
 uses.update(roots);cache={};peak=0
 for i,node in enumerate(nodes):
  if node[0]=='c':p=const(Q(node[1],node[2]))
  elif node[0]=='input':p=var(node[1])
  elif node[0]=='a':p=plus(cache[node[1]],cache[node[2]])
  elif node[0]=='m':p=times(cache[node[1]],cache[node[2]])
  else:raise AssertionError(node[0])
  assert len(p)<=pre['all_node_support_upper'][i];cache[i]=p;peak=max(peak,sum(map(len,cache.values())))
  if node[0]in['a','m']:
   for child in node[1:]:
    uses[child]-=1
    if not uses[child]:del cache[child]
  if not uses[i]:del cache[i]
 assert peak<=6544
 for j,(outrow,polyrow) in enumerate(zip(data['outputs3'],data['quartic_polynomial_nodes'])):
  assert outrow['physical_row']==polyrow['physical_row']==piv[j]
  assert [cache[n] for n in outrow['NF4_full6_nodes']]==expected_nf[j]
  assert [cache[n] for n in outrow['residual5_nodes']]==expected_res[j]
  for side,key in enumerate(['WA4_polynomial','WB4_polynomial']):
   assert {tuple(m):cache[n] for m,n in polyrow[key]}==expected_polys[j][side]
 # Literal source zeros for the other9; all inactive mixed coefficients
 # remain unspecified and multiply E_i3=0, never guessed zero.
 zero_rows=[]
 for p in rows:
  if p in piv:continue
  pos=rows.index(p);ds=dict(C[p]['delta'])
  assert not any(map(Q,fm[p]['bB_second_plane']))
  assert not any(betaQ[i][k]*Q(fm[p]['bA_first_z7'][i]) for i in range(7) for k in range(3))
  assert all(not Q(ds.get(ix,'0')) for ix in [ir3,ir2l]+[ix for i in active for ix in [izr[i],izl[i]]])
  zero_rows.append(p)
 assert zero_rows==data['all9_nonpivot_tauAB_quartic_classes_zero_source_rows']
 proof={'input_sha256':pins,'all39_actual_cached_leaf_locators_verified':True,
  'formal_symbols_independent_no_numerical_pairing':True,
  'independent_primitive_quartic_formulas':['A=-nuB2*bB2','B=-(nuA2-nuA_z*beta2)*u2-delta2_mixed*E3'],
  'full_H3_correction_retained':['(delta30*r+delta21*ell)*P0,3','lambda12*ell*P1,3','mu03*ell*P2,3'],
  'all_mu_squared_and_beta_squared_carry_retained':True,
  'all3_WA_WB_polynomials_18fullNF_15residual_polynomials_match':True,
  'all9_nonpivot_source_zero_rows':zero_rows,
  'order_proof':[
   'nuA/B firstplane and origins vanish; firstframe bA2plane=delta2plane=0. Hence graph-evaluated primitive quartics need only nu/b secondjets plus b firstz. nuA_B2=nuA2-nuA_z beta2 retains beta² carry.',
   'H2mixed clears by plane-linear*RH_i, producing -delta_j2mixed E_i3 in B4. H2zz clears have at least2 z factors in their cubic B corrections; graph substitution raises them to>=5.',
   'H3mixed clears by quadratic*RH_i. Its B2plane=0 and graph-substituted B2mixed starts3, so all A/B effects start>=5. H4, including unknowndelta4, clears by quadratic*P whose A/B starts3, again>=5.',
   'Every planeH3 term is canceled using actual P0/P1/P2 with ALLfour coefficients. In particular mu03ellP2 is retained, not a quotient cut.',
   'After coherent second-order z normalization, pureplane rewrite matrices are diag(0,-beta_i2,0) and commute. Terms containing z in these second-order matrices can contribute only>=5 to their quartic ambiguity. Other H3 decompositions differ by actual HB sources, killed by the five residual duals at accepted full HB values.'
  ],
  'expanded_symbolic_support_peak':peak,'prior_whole_bound':pre['whole_symbolic_proof_upper'],
  'max_formal_QQ_coefficient_bits':MAXCOF,'max_before_formal_operation_bits':MAXPRE,
  'all3_fullNF_symbol_polynomials':[[encoded(p) for p in row] for row in expected_nf],
  'all15_residual_symbol_polynomials':[[encoded(p) for p in row] for row in expected_res],
  'code_sha256':sha(__file__),'response_values_not_evaluated':True,'elapsed_seconds':time.monotonic()-start}
 ep=out/'independent-full-tauAB-scalar-symbol-polynomial-proof.json';save(ep,proof)
 save(out/'formal-phase.json',{'status':'FORMAL_SOURCE_PASS_WAIT_WHOLE_WORKER_GATE','proof_sha256':sha(ep),'response_values_not_evaluated':True})
 signal.alarm(0);print(json.dumps({'formal_source_passed':True,'proof_sha256':sha(ep),'seconds':proof['elapsed_seconds'],'peak':peak,'maxpre':MAXPRE}))
if __name__=='__main__':main()
