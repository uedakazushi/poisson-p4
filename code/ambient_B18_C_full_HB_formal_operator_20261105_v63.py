#!/usr/bin/env python3
"""Literal full-HB affine coefficients; three pivot ell3 values remain formal."""
import hashlib,json,math,re
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATHS={
 'HB':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json',
 'HB_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json',
 'old_value':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase1-six-root-existing-cache-readbacks.json',
 'old_value_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'beta':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'reduction':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json',
 'eta':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63/report.json',
 'eta_peer':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'graph':'reports/AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63/report.json',
 'graph_peer':'reports/AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63_REVIEW/report.json',
}
PIV=[18,60,67]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def add(a,b):
 v=a.copy()
 for k,c in b.items():
  v[k]=v.get(k,Q(0))+c
  if not v[k]:del v[k]
 return v
def scale(a,c):return{k:v*c for k,v in a.items()if v*c}
def leaf(kind,row,index):return{(kind,row,index):Q(1)}
def encode(a):return[[list(k),str(v)]for k,v in sorted(a.items())]

def build():
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};pins={p:sha(ROOT/p)for p in PATHS.values()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['eta']['nonzero_eta_count']==0 and d['graph']['all_five_graph03_coefficients_zero_source_proved']
 hb=d['HB'];T=[[Q(x)for x in row]for row in hb['constant_source_normalization_T']]
 H=[[Q(x)for x in row]for row in hb['source_first7_h_rows']]
 raw=d['raw33']['actual33_basis'];deltas={r['physical_row']:dict(r['delta'])for r in d['C2']['physical12_coefficients']}
 beta=[]
 for row in d['beta']['beta2_before_ell']:
  coeff=[Q(0)]*3
  for mon,s in row:
   assert not any(mon[:7])and sum(mon)==2
   nr,ns=mon[7:]
   for k in range(ns+1):coeff[k]+=Q(s)*math.comb(ns,k)/5**k
  beta.append(coeff)
 assert[i for i,row in enumerate(beta)if any(row)]==[0,4,5]
 D=[];E=[];known={}
 for row in PIV:
  ds=[leaf('cache',row,29),leaf('cache',row,30),leaf('lambda12',row,-1),leaf('mu03',row,-1)];D.append(ds)
  for ix in[29,30]:known[('cache',row,ix)]=Q(deltas[row].get(ix,'0'))
  gamma=[{}for _ in range(4)]
  for z in[0,4,5]:
   ri,ei=raw.index([z+1,1,0]),raw.index([z+1,0,1]);kr,ke=('cache',row,ri),('cache',row,ei)
   known[kr]=Q(deltas[row].get(ri,'0'));known[ke]=Q(deltas[row].get(ei,'0'))
   for m in range(3):
    gamma[m]=add(gamma[m],scale(leaf(*kr),beta[z][m]))
    gamma[m+1]=add(gamma[m+1],scale(leaf(*ke),beta[z][m]))
  E.append([add(a,scale(b,Q(-1)))for a,b in zip(ds,gamma)])
 for row in d['old_value']['all3_lambda_readbacks']:known[('lambda12',row['physical_row'],-1)]=Q(row['delta_r_ell2'])
 # All other first-seven E_i3 vanish: inactive beta2=0; free delta3=0;
 # active mixed leaves of free rows are literal0, even on pureell3.
 for row in[19,28,29,72]:
  for z in[0,4,5]:
   for mon in[[z+1,1,0],[z+1,0,1]]:assert Q(deltas[row].get(raw.index(mon),'0'))==0
 PA=[[{}for _ in range(4)]for _ in range(3)];PB=[[{}for _ in range(4)]for _ in range(3)]
 for k in range(3):
  for p in range(3):
   for m in range(4):PA[k][m]=add(PA[k][m],scale(D[p][m],T[p][k]))
  for q,z in enumerate([0,4,5]):
   coeff=sum((T[p][k]*H[p][z]for p in range(3)),Q(0))
   for m in range(4):PB[k][m]=add(PB[k][m],scale(E[q][m],coeff))
 f3=[[Q(x)for x in row]for row in d['reduction']['full_no_cut_cubic_NF6_rows12_columns']]
 f4=[[Q(0)]*15 for _ in range(6)]
 for side in range(2):
  a1,a2,a3=[-f3[3*side+k][4*side]for k in range(3)]
  block=[[a1*a1-a2,-a1,1,0,0],[a1*a2-a3,-a2,0,1,0],[a1*a3,-a3,0,0,1]]
  for k in range(3):f4[3*side+k][5*side:5*side+5]=list(map(Q,block[k]))
 matrix=[[{}for _ in range(2)]for _ in range(6)]
 for j in range(2):
  for side,P in enumerate([PA,PB]):
   poly=[{}for _ in range(5)]
   for m in range(4):poly[m+1]=add(poly[m+1],P[j][m]);poly[m]=add(poly[m],scale(P[j+1][m],Q(-1)))
   for k in range(3):
    for m in range(5):matrix[3*side+k][j]=add(matrix[3*side+k][j],scale(poly[m],f4[3*side+k][5*side+m]))
 # Full source projection is literally the older complete fixed affine
 # operator: compare symbolic coefficients, not just cached numbers.
 for side,k in enumerate([0,3]):
  for j in range(2):
   old={('lambda12'if key[0]=='lambda'else key[0],key[1],key[2]):Q(v)for key,v in hb['matrix2_rows_A4_B4_columns_HB0_HB1_linear_leaf_operators'][side][j]}
   assert matrix[k][j]==old
 mu={(kind,row,index)for cols in matrix for a in cols for kind,row,index in a if kind=='mu03'}
 assert mu=={('mu03',r,-1)for r in PIV}
 # Only supplied OLD leaf values are copied. The three mu03 remain formal.
 assert len(known)==27 and all(k[0]!='mu03'for k in known)
 data={
  'scope':'TWO genuine full-C known quartic HB sources, all SIX NF coordinates retained. Not complete quartic source space.',
  'rows':['A*r²ell²','A*r ell³','A*ell4','B*r²ell²','B*r ell³','B*ell4'],
  'columns':['ell*P0-r*P1','ell*P1-r*P2'],
  'full_NF4':[[str(x)for x in row]for row in f4],
  'P3_A':[[encode(a)for a in row]for row in PA],'P3_B':[[encode(a)for a in row]for row in PB],
  'full_HB_6_by2_affine_operators':[[encode(a)for a in cols]for cols in matrix],
  'known27_old_cached_leaf_values':[[list(k),str(v)]for k,v in sorted(known.items())],
  'new_three_mu03_callbacks':[[list(k),'UNREAD_FORMAL']for k in sorted(mu)],
  'old_projection_symbolic_operator_identity_verified':True,
  'full_source_carry_proof':'ActualPk secondterms are exactly r²H,r ell H,ell²H, so both HB have FULL order3 zero beforeanycut. Theirorder4 onlyuseslinear*Pk3. Allz-containingthirdpartsretainz; zframecarrybegin5. Hquarticsclearbydegree2*actualP2 withA/Bcarry>=5. Negative RB beta2mixed carry is inE_i3; no beta3/unknownmixedsource3/δ4/newIFT.',
  'no_literal_ell3H_cut_used':True,'all_r_ell3_and_ell4_coordinates_retained':True,
  'source_reduction_to_three_pivots':'Allfree/graphdelta plane3 nowzero INCLUDINGpureell3. Inactivebeta2zero; all24freeactiveδ2mixed0. Henceonlyi0,4,5 (physical18/60/67) remain inHᵗT*E3.',
  'rank_rule':'Row0 is(-3/5,6). Computew=col1+10col0. Ifanyw!=0, rank2andminor=(-3/5)*w_j; ifall0, theseTWOcolumnsrank1 only. No wholequarticrank claimifrank1.',
  'rank2_consequence':'ActualfullHF4=4 byknownupper6-2 andD-ranklower4; Tor(m²C)=(2,2), so C long exponent e<=4 and>=3. NOTautomatic e3.',
  'rank1_scope':'OnlyHBknownsubmodule rank1; actualfullHF4 remains4or5. No exactb/e/totalCtorsion result.',
  'new_values_computed':False,
 }
 return data,pins,count(d)
